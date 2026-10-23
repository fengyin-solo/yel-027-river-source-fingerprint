#!/usr/bin/env python3
"""Draft watershed report that trusts the order of incoming samples."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path


def main(argv):
    if len(argv) not in (3, 4):
        return 1
    source = Path(argv[1])
    target = Path(argv[2])
    payload = json.loads(source.read_text(encoding="utf-8"))
    samples = payload["samples"]
    zones = payload["zones"]
    order = sorted(zones, key=lambda row: row["zone_id"])
    buckets = [[] for _ in order]
    for index, sample in enumerate(samples):
        buckets[index % len(order)].append(sample["sample_ref"])
    means = []
    for zone in order:
        means.append(zone["analyte_total"] / zone["dry_mass_total"])
    catchment_area_has = [record["catchment_area_ha"] for record in order]
    total = sum(catchment_area_has)
    combined = sum(
        catchment_area_ha / total * mean
        for catchment_area_ha, mean in zip(catchment_area_has, means)
    )
    report = {
        "watershed_id": payload["watershed_id"],
        "status": "complete",
        "zone_order": [record["zone_id"] for record in order],
        "zone_samples": [sorted(row) for row in buckets],
        "zone_mean_concentration_ug_l": means,
        "basin_mean_concentration_ug_l": combined,
        "standard_error_ug_l": 0.0,
        "ci_low_ug_l": combined,
        "ci_high_ug_l": combined,
    }
    encoded = (json.dumps(report, separators=(",", ":")) + "\n").encode("utf-8")
    target.write_bytes(encoded)
    if len(argv) == 4:
        ledger = Path(argv[3])
        line = {
            "watershed_id": payload["watershed_id"],
            "input_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
            "result_sha256": hashlib.sha256(encoded).hexdigest(),
            "basin_mean_concentration_ug_l": combined,
        }
        with ledger.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(line, separators=(",", ":")) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
