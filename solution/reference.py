#!/usr/bin/env python3
"""Recover missing watershed-sample zones and publish a design-based report."""

from __future__ import annotations

import hashlib
import itertools
import json
import math
import os
import re
import sys
import tempfile
from pathlib import Path


TOP_KEYS = {
    "watershed_id",
    "campaign_id",
    "samples",
    "zones",
    "calibration",
    "design_confidence",
}
PANEL_KEYS = {"sample_ref", "aliquot_count", "dry_mass_g", "concentrations_ug_l"}
STRATUM_KEYS = {
    "zone_id",
    "catchment_area_ha",
    "sample_count",
    "dry_mass_total",
    "analyte_total",
    "analyte_sq_total",
}
CALIBRATION_KEYS = {"z_score", "tie_break"}
LEDGER_KEYS = {"watershed_id", "input_sha256", "result_sha256", "basin_mean_concentration_ug_l"}
PANEL_PATTERN = re.compile(r"^[A-Z][A-Z0-9_-]{2,15}$")
STRATUM_PATTERN = re.compile(r"^[A-Z][A-Z0-9_-]{0,23}$")
HEX_DIGITS = set("0123456789abcdef")


class Refusal(Exception):
    """A user-facing refusal with a stable business code."""

    def __init__(self, kind: str) -> None:
        super().__init__(kind)
        self.kind = kind


def fail(kind: str) -> None:
    sys.stderr.write(f"运行拒收: {kind}\n")
    raise SystemExit(1)


def finite(value: object) -> bool:
    return (
        isinstance(value, (int, float))
        and not isinstance(value, bool)
        and math.isfinite(float(value))
    )


def integer(value: object, low: int, high: int) -> bool:
    return (
        isinstance(value, int)
        and not isinstance(value, bool)
        and low <= value <= high
    )


def strict_pairs(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise Refusal("BAD_INPUT")
        result[key] = value
    return result


def ledger_pairs(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise Refusal("BAD_LEDGER")
        result[key] = value
    return result


def visible_label(value: object, low: int, high: int) -> bool:
    return (
        isinstance(value, str)
        and low <= len(value) <= high
        and value == value.strip()
        and all(ord(character) >= 32 and ord(character) != 127 for character in value)
    )


def validate(payload):
    if not isinstance(payload, dict) or set(payload) != TOP_KEYS:
        raise Refusal("BAD_INPUT")
    watershed_id = payload["watershed_id"]
    campaign_id = payload["campaign_id"]
    if not visible_label(watershed_id, 6, 48) or not watershed_id.startswith("RIVER-"):
        raise Refusal("BAD_INPUT")
    if not visible_label(campaign_id, 2, 48):
        raise Refusal("BAD_INPUT")
    samples = payload["samples"]
    if not isinstance(samples, list) or not 12 <= len(samples) <= 24:
        raise Refusal("BAD_INPUT")
    seen_samples = set()
    clean_samples = []
    for row in samples:
        if not isinstance(row, dict) or set(row) != PANEL_KEYS:
            raise Refusal("BAD_INPUT")
        reference = row["sample_ref"]
        if (
            not isinstance(reference, str)
            or PANEL_PATTERN.fullmatch(reference) is None
            or reference in seen_samples
        ):
            raise Refusal("BAD_INPUT")
        seen_samples.add(reference)
        count = row["aliquot_count"]
        if not integer(count, 10, 5000):
            raise Refusal("BAD_INPUT")
        weight = row["dry_mass_g"]
        if not finite(weight) or not 0 < float(weight) <= 100:
            raise Refusal("BAD_INPUT")
        concentrations_ug_l = row["concentrations_ug_l"]
        if not isinstance(concentrations_ug_l, list) or len(concentrations_ug_l) != count:
            raise Refusal("BAD_INPUT")
        if any(not integer(concentration, 0, 1000) for concentration in concentrations_ug_l):
            raise Refusal("BAD_INPUT")
        numeric_concentrations_ug_l = [int(concentration) for concentration in concentrations_ug_l]
        clean_samples.append(
            {
                "sample_ref": reference,
                "aliquot_count": count,
                "dry_mass_g": float(weight),
                "concentrations_ug_l": numeric_concentrations_ug_l,
            }
        )
    zones = payload["zones"]
    if not isinstance(zones, list) or not 3 <= len(zones) <= 5:
        raise Refusal("BAD_INPUT")
    seen_zones = set()
    clean_zones = []
    for row in zones:
        if not isinstance(row, dict) or set(row) != STRATUM_KEYS:
            raise Refusal("BAD_INPUT")
        identifier = row["zone_id"]
        if (
            not isinstance(identifier, str)
            or STRATUM_PATTERN.fullmatch(identifier) is None
            or identifier in seen_zones
        ):
            raise Refusal("BAD_INPUT")
        seen_zones.add(identifier)
        catchment_area_ha = row["catchment_area_ha"]
        if not integer(catchment_area_ha, 2, 1_000_000_000):
            raise Refusal("BAD_INPUT")
        sample_count = row["sample_count"]
        if not integer(sample_count, 1, len(clean_samples)):
            raise Refusal("BAD_INPUT")
        numeric = {}
        for key in (
            "dry_mass_total",
            "analyte_total",
            "analyte_sq_total",
        ):
            value = row[key]
            if not finite(value) or float(value) < 0:
                raise Refusal("BAD_INPUT")
            numeric[key] = float(value)
        if numeric["dry_mass_total"] <= 0 or numeric["analyte_total"] <= 0:
            raise Refusal("BAD_INPUT")
        clean_zones.append(
            {
                "zone_id": identifier,
                "catchment_area_ha": catchment_area_ha,
                "sample_count": sample_count,
                **numeric,
            }
        )
    if sum(row["sample_count"] for row in clean_zones) != len(clean_samples):
        raise Refusal("BAD_INPUT")
    calibration = payload["calibration"]
    if not isinstance(calibration, dict) or set(calibration) != CALIBRATION_KEYS:
        raise Refusal("BAD_INPUT")
    if not finite(calibration["z_score"]) or not 0 < float(calibration["z_score"]) <= 10:
        raise Refusal("BAD_INPUT")
    if calibration["tie_break"] != "lexicographic":
        raise Refusal("BAD_INPUT")
    confidence = payload["design_confidence"]
    if not finite(confidence) or not 0.8 <= float(confidence) <= 0.99:
        raise Refusal("BAD_INPUT")
    return (
        watershed_id,
        campaign_id,
        sorted(clean_samples, key=lambda item: item["sample_ref"]),
        sorted(clean_zones, key=lambda item: item["zone_id"]),
        float(calibration["z_score"]),
        float(confidence),
    )


def sample_signature(sample):
    count = sample["aliquot_count"]
    weight = sample["dry_mass_g"]
    concentration_total = sum(sample["concentrations_ug_l"])
    square_total = sum(value * value for value in sample["concentrations_ug_l"])
    return (
        1.0,
        weight * count,
        weight * concentration_total,
        weight * square_total,
    )


def zone_signature(zone):
    return (
        float(zone["sample_count"]),
        zone["dry_mass_total"],
        zone["analyte_total"],
        zone["analyte_sq_total"],
    )


def same_signature(left, right, tolerance=1e-7):
    return all(
        abs(left[index] - right[index]) <= tolerance
        for index in range(len(left))
    )


def enumerate_candidates(samples, zone):
    target = zone_signature(zone)
    count = zone["sample_count"]
    candidates = []
    for indexes in itertools.combinations(range(len(samples)), count):
        total = [0.0, 0.0, 0.0, 0.0]
        for index in indexes:
            signature = sample_signature(samples[index])
            for offset in range(4):
                total[offset] += signature[offset]
        if same_signature(total, target):
            candidates.append(indexes)
    candidates.sort(
        key=lambda row: tuple(samples[index]["sample_ref"] for index in row)
    )
    return candidates


def recover_partition(samples, zones):
    by_zone = {
        row["zone_id"]: enumerate_candidates(samples, row)
        for row in zones
    }
    if any(not choices for choices in by_zone.values()):
        raise Refusal("BAD_MODEL")

    def search(position, used, assignment):
        if position == len(zones):
            return dict(assignment)
        zone = zones[position]
        for candidate in by_zone[zone["zone_id"]]:
            if any(index in used for index in candidate):
                continue
            updated = set(used)
            updated.update(candidate)
            next_assignment = dict(assignment)
            next_assignment[zone["zone_id"]] = candidate
            found = search(position + 1, updated, next_assignment)
            if found is not None:
                return found
        return None

    recovered = search(0, set(), {})
    if recovered is None:
        raise Refusal("BAD_MODEL")
    return recovered


def calculate_report(watershed_id, samples, zones, z_score):
    assignment = recover_partition(samples, zones)
    total_catchment_area_ha = sum(row["catchment_area_ha"] for row in zones)
    means = []
    variances = []
    references = []
    for zone in zones:
        chosen = [samples[index] for index in assignment[zone["zone_id"]]]
        chosen.sort(key=lambda item: item["sample_ref"])
        mean = zone["analyte_total"] / zone["dry_mass_total"]
        dry_mass_total = zone["dry_mass_total"]
        weight_square_sum = sum(
            sample["aliquot_count"] * sample["dry_mass_g"] ** 2
            for sample in chosen
        )
        effective_n = dry_mass_total * dry_mass_total / weight_square_sum
        centered = sum(
            sample["dry_mass_g"] * (concentration - mean) ** 2
            for sample in chosen
            for concentration in sample["concentrations_ug_l"]
        )
        denominator = dry_mass_total - weight_square_sum / dry_mass_total
        sample_variance = centered / denominator
        share = zone["catchment_area_ha"] / total_catchment_area_ha
        fpc = max(0.0, 1.0 - effective_n / zone["catchment_area_ha"])
        variance = share * share * fpc * sample_variance / effective_n
        means.append(mean)
        variances.append(variance)
        references.append([sample["sample_ref"] for sample in chosen])
    combined = sum(
        zone["catchment_area_ha"] / total_catchment_area_ha * mean
        for zone, mean in zip(zones, means)
    )
    standard_error_ug_l = math.sqrt(sum(variances))
    return {
        "watershed_id": watershed_id,
        "status": "complete",
        "zone_order": [row["zone_id"] for row in zones],
        "zone_samples": references,
        "zone_mean_concentration_ug_l": means,
        "basin_mean_concentration_ug_l": combined,
        "standard_error_ug_l": standard_error_ug_l,
        "ci_low_ug_l": combined - z_score * standard_error_ug_l,
        "ci_high_ug_l": combined + z_score * standard_error_ug_l,
    }


def valid_sha(value):
    return (
        isinstance(value, str)
        and len(value) == 64
        and all(character in HEX_DIGITS for character in value)
    )


def read_ledger(path):
    if not path.exists():
        return []
    try:
        raw = path.read_bytes()
        text = raw.decode("utf-8")
    except (OSError, UnicodeError):
        raise Refusal("BAD_LEDGER")
    if not text:
        return []
    if not text.endswith("\n"):
        raise Refusal("BAD_LEDGER")
    records = []
    for raw_line in text.splitlines():
        try:
            record = json.loads(raw_line, object_pairs_hook=ledger_pairs)
        except (json.JSONDecodeError, UnicodeError, ValueError):
            raise Refusal("BAD_LEDGER")
        if not isinstance(record, dict) or list(record) != [
            "watershed_id",
            "input_sha256",
            "result_sha256",
            "basin_mean_concentration_ug_l",
        ]:
            raise Refusal("BAD_LEDGER")
        canonical = json.dumps(record, ensure_ascii=False, separators=(",", ":"))
        if raw_line != canonical:
            raise Refusal("BAD_LEDGER")
        if (
            not visible_label(record["watershed_id"], 6, 48)
            or not valid_sha(record["input_sha256"])
            or not valid_sha(record["result_sha256"])
            or not finite(record["basin_mean_concentration_ug_l"])
        ):
            raise Refusal("BAD_LEDGER")
        records.append(record)
    return records


def atomic_write(path, payload):
    descriptor = None
    stage = None
    try:
        descriptor, name = tempfile.mkstemp(prefix=".watershed-", dir=path.parent)
        stage = Path(name)
        with os.fdopen(descriptor, "wb") as handle:
            descriptor = None
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(stage, path)
        stage = None
    finally:
        if descriptor is not None:
            os.close(descriptor)
        if stage is not None:
            stage.unlink(missing_ok=True)


def publish_result(target, payload):
    if (
        not target.parent.is_dir()
        or target.is_dir()
        or not os.access(target.parent, os.W_OK)
    ):
        raise Refusal("IO_FAILURE")
    try:
        atomic_write(target, payload)
    except OSError:
        raise Refusal("IO_FAILURE")


def publish_pair(source, target, ledger, input_bytes, output_bytes, combined):
    if (
        not target.parent.is_dir()
        or not ledger.parent.is_dir()
        or target.is_dir()
        or ledger.is_dir()
        or not os.access(target.parent, os.W_OK)
        or not os.access(ledger.parent, os.W_OK)
    ):
        raise Refusal("IO_FAILURE")
    if ledger in (source, target):
        raise Refusal("BAD_LEDGER")
    records = read_ledger(ledger)
    input_sha = hashlib.sha256(input_bytes).hexdigest()
    result_sha = hashlib.sha256(output_bytes).hexdigest()
    for record in records:
        if record["input_sha256"] == input_sha:
            if (
                record["result_sha256"] == result_sha
                and record["watershed_id"] == json.loads(output_bytes)["watershed_id"]
                and record["basin_mean_concentration_ug_l"] == combined
            ):
                publish_result(target, output_bytes)
                return
            raise Refusal("BAD_LEDGER")
        if (
            record["result_sha256"] == result_sha
            and record["input_sha256"] != input_sha
        ):
            raise Refusal("BAD_LEDGER")
    new_record = {
        "watershed_id": json.loads(output_bytes)["watershed_id"],
        "input_sha256": input_sha,
        "result_sha256": result_sha,
        "basin_mean_concentration_ug_l": combined,
    }
    new_line = (
        json.dumps(new_record, ensure_ascii=False, separators=(",", ":")) + "\n"
    ).encode("utf-8")
    old_target = target.read_bytes() if target.exists() else None
    old_ledger = ledger.read_bytes() if ledger.exists() else b""
    ledger_bytes = old_ledger + new_line
    try:
        atomic_write(ledger, ledger_bytes)
        try:
            atomic_write(target, output_bytes)
        except OSError:
            if ledger.exists():
                if old_ledger:
                    atomic_write(ledger, old_ledger)
                else:
                    ledger.unlink(missing_ok=True)
            raise
    except OSError:
        if old_target is None:
            target.unlink(missing_ok=True)
        else:
            try:
                atomic_write(target, old_target)
            except OSError:
                pass
        raise Refusal("IO_FAILURE")


def main(argv):
    if len(argv) not in (3, 4):
        fail("BAD_INPUT")
    source = Path(argv[1])
    target = Path(argv[2])
    try:
        input_bytes = source.read_bytes()
        payload = json.loads(
            input_bytes.decode("utf-8"),
            object_pairs_hook=strict_pairs,
        )
        watershed_id, _, samples, zones, z_score, _ = validate(payload)
        report = calculate_report(watershed_id, samples, zones, z_score)
        encoded = (
            json.dumps(report, ensure_ascii=False, separators=(",", ":")) + "\n"
        ).encode("utf-8")
        if len(argv) == 4:
            publish_pair(
                source,
                target,
                Path(argv[3]),
                input_bytes,
                encoded,
                report["basin_mean_concentration_ug_l"],
            )
        else:
            publish_result(target, encoded)
        return 0
    except Refusal as exc:
        fail(exc.kind)
    except (OSError, UnicodeError, json.JSONDecodeError, ValueError):
        fail("BAD_INPUT")


if __name__ == "__main__":
    sys.exit(main(sys.argv))
