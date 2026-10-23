"""Acceptance tests for missing-zone watershed reconstruction."""

from __future__ import annotations

import json
import math
import os
import subprocess
import sys
from pathlib import Path

import pytest


ENTRY = os.environ.get("RIVER_SOURCE_ENTRY", "/app/sediment_source_fingerprint.py")


def sample(reference, count, weight, start):
    concentrations_ug_l = [((start + index * 3) % 91) + 5 for index in range(count)]
    return {
        "sample_ref": reference,
        "aliquot_count": count,
        "dry_mass_g": weight,
        "concentrations_ug_l": concentrations_ug_l,
    }


def make_document(sample_total=12, zone_total=3, shuffle=False):
    if (sample_total, zone_total) == (12, 3):
        groups = [
            [
                sample("SED-A01", 10, 1.25, 2),
                sample("SED-A02", 14, 2.5, 13),
                sample("SED-A03", 11, 3.25, 31),
                sample("SED-A04", 17, 1.75, 47),
            ],
            [
                sample("SED-B01", 13, 1.5, 7),
                sample("SED-B02", 16, 2.25, 23),
                sample("SED-B03", 12, 4.0, 39),
                sample("SED-B04", 19, 0.75, 55),
            ],
            [
                sample("SED-C01", 15, 1.1, 11),
                sample("SED-C02", 18, 2.1, 29),
                sample("SED-C03", 10, 5.0, 43),
                sample("SED-C04", 20, 1.4, 61),
            ],
        ]
        catchment_area_has = [1200, 2600, 4100]
        identifiers = ["A", "B", "C"]
    else:
        groups = []
        per_group = sample_total // zone_total
        catchment_area_has = []
        identifiers = []
        ordinal = 0
        for group_index in range(zone_total):
            rows = []
            for member_index in range(per_group):
                ordinal += 1
                rows.append(
                    sample(
                        f"P-{chr(65 + group_index)}{member_index:02d}",
                        10 + (ordinal % 9),
                        0.5 + ((ordinal * 7) % 19) / 4,
                        3 + ordinal * 11 + group_index * 5,
                    )
                )
            groups.append(rows)
            catchment_area_has.append(900 + group_index * 1100)
            identifiers.append(chr(65 + group_index))
    truth = {
        identifier: [row["sample_ref"] for row in rows]
        for identifier, rows in zip(identifiers, groups)
    }
    samples = [row for group in groups for row in group]
    if shuffle:
        samples = samples[6:9] + samples[:6] + samples[9:]
    zones = []
    for identifier, rows, catchment_area_ha in zip(identifiers, groups, catchment_area_has):
        dry_mass_total = sum(row["aliquot_count"] * row["dry_mass_g"] for row in rows)
        concentration_total = sum(
            row["dry_mass_g"] * sum(row["concentrations_ug_l"]) for row in rows
        )
        square_total = sum(
            row["dry_mass_g"] * sum(value * value for value in row["concentrations_ug_l"])
            for row in rows
        )
        zones.append(
            {
                "zone_id": identifier,
                "catchment_area_ha": catchment_area_ha,
                "sample_count": len(rows),
                "dry_mass_total": dry_mass_total,
                "analyte_total": concentration_total,
                "analyte_sq_total": square_total,
            }
        )
    document = {
        "watershed_id": "RIVER-BASIN-04",
        "campaign_id": "CAMPAIGN-AUTUMN",
        "samples": samples,
        "zones": zones,
        "calibration": {"z_score": 1.959963984540054, "tie_break": "lexicographic"},
        "design_confidence": 0.95,
    }
    return document, truth


def recalibrate(document):
    groups = {}
    for sample_row in document["samples"]:
        prefix = sample_row["sample_ref"].split("-")[1][0]
        groups.setdefault(prefix, []).append(sample_row)
    for zone in document["zones"]:
        rows = groups[zone["zone_id"]]
        zone["sample_count"] = len(rows)
        zone["dry_mass_total"] = sum(
            row["aliquot_count"] * row["dry_mass_g"] for row in rows
        )
        zone["analyte_total"] = sum(
            row["dry_mass_g"] * sum(row["concentrations_ug_l"]) for row in rows
        )
        zone["analyte_sq_total"] = sum(
            row["dry_mass_g"] * sum(value * value for value in row["concentrations_ug_l"])
            for row in rows
        )
    return document


def write_document(path, document, raw=None):
    if raw is None:
        raw = json.dumps(document, ensure_ascii=False, separators=(",", ":")) + "\n"
    path.write_bytes(raw.encode("utf-8"))
    return path


def invoke(tmp_path, document, raw=None, with_ledger=False):
    source = write_document(tmp_path / "watershed.json", document, raw)
    target = tmp_path / "report.json"
    ledger = tmp_path / "delivery.jsonl"
    argv = [sys.executable, ENTRY, str(source), str(target)]
    if with_ledger:
        argv.append(str(ledger))
    result = subprocess.run(argv, capture_output=True, text=True, timeout=30)
    return result, source, target, ledger


def expected_values(document, truth):
    by_reference = {row["sample_ref"]: row for row in document["samples"]}
    zones = sorted(document["zones"], key=lambda row: row["zone_id"])
    total_catchment_area_ha = sum(row["catchment_area_ha"] for row in zones)
    means = []
    variances = []
    for zone in zones:
        rows = [by_reference[value] for value in truth[zone["zone_id"]]]
        mean = zone["analyte_total"] / zone["dry_mass_total"]
        weight_square_sum = sum(
            row["aliquot_count"] * row["dry_mass_g"] ** 2 for row in rows
        )
        effective_n = zone["dry_mass_total"] ** 2 / weight_square_sum
        centered = sum(
            row["dry_mass_g"] * (concentration - mean) ** 2
            for row in rows
            for concentration in row["concentrations_ug_l"]
        )
        sample_variance = centered / (
            zone["dry_mass_total"] - weight_square_sum / zone["dry_mass_total"]
        )
        share = zone["catchment_area_ha"] / total_catchment_area_ha
        fpc = max(0.0, 1.0 - effective_n / zone["catchment_area_ha"])
        variances.append(share * share * fpc * sample_variance / effective_n)
        means.append(mean)
    combined = sum(
        row["catchment_area_ha"] / total_catchment_area_ha * mean
        for row, mean in zip(zones, means)
    )
    standard_error_ug_l = math.sqrt(sum(variances))
    z_score = document["calibration"]["z_score"]
    return means, combined, standard_error_ug_l, combined - z_score * standard_error_ug_l, combined + z_score * standard_error_ug_l


def read_report(path):
    return json.loads(path.read_text(encoding="utf-8"))


def test_success_is_quiet_and_creates_report(tmp_path):
    """题面：「正常返回码为 0 时，标准输出与标准错误均为空；」；"""
    document, _ = make_document()
    result, _, target, _ = invoke(tmp_path, document)
    assert result.returncode == 0 and result.stdout == "" and result.stderr == ""
    assert target.is_file()


def test_existing_compact_json_suffix_is_preserved(tmp_path):
    """题面：「结果文件只允许一个 JSON 对象和结尾换行」。"""
    document, _ = make_document()
    result, _, target, _ = invoke(tmp_path, document)
    raw = target.read_bytes()
    assert result.returncode == 0 and raw.endswith(b"\n") and raw.count(b"\n") == 1


def test_report_has_exact_ordered_keys(tmp_path):
    """题面：「结果文件只允许一个 JSON 对象和结尾换行」。"""
    document, _ = make_document()
    result, _, target, _ = invoke(tmp_path, document)
    assert result.returncode == 0
    assert list(read_report(target)) == [
        "watershed_id",
        "status",
        "zone_order",
        "zone_samples",
        "zone_mean_concentration_ug_l",
        "basin_mean_concentration_ug_l",
        "standard_error_ug_l",
        "ci_low_ug_l",
        "ci_high_ug_l",
    ]


def test_hidden_partition_is_recovered(tmp_path):
    """题面：「所有样本必须完整划分到各来源区，任何一份都不得重复或遗漏；」；"""
    document, truth = make_document()
    result, _, target, _ = invoke(tmp_path, document)
    references = read_report(target)["zone_samples"]
    assert result.returncode == 0
    assert references == [sorted(truth[key]) for key in ("A", "B", "C")]


def test_input_order_does_not_change_business_bytes(tmp_path):
    """题面：「归属并列的规范化顺序由区编号序列和区内 `sample_ref` 序列共同组成，取字典序最小者；」"""
    first, _ = make_document()
    shuffled, _ = make_document(shuffle=True)
    result_one, _, target_one, _ = invoke(tmp_path, first)
    write_document(tmp_path / "watershed.json", shuffled)
    result_two = subprocess.run(
        [sys.executable, ENTRY, str(tmp_path / "watershed.json"), str(tmp_path / "second.json")],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result_one.returncode == 0 and result_two.returncode == 0
    assert target_one.read_bytes() == (tmp_path / "second.json").read_bytes()


def test_zone_mean_concentration_ug_l_use_published_dry_mass_totals(tmp_path):
    """题面：「zone_mean_h = analyte_total_h / dry_mass_total_h」；"""
    document, truth = make_document()
    result, _, target, _ = invoke(tmp_path, document)
    expected = expected_values(document, truth)[0]
    assert result.returncode == 0
    assert read_report(target)["zone_mean_concentration_ug_l"] == pytest.approx(expected)


def test_basin_mean_concentration_ug_l_uses_catchment_area_ha_shares(tmp_path):
    """题面：「basin_mean_concentration_ug_l = sum(share_h * zone_mean_h)」；"""
    document, truth = make_document()
    result, _, target, _ = invoke(tmp_path, document)
    expected = expected_values(document, truth)[1]
    assert result.returncode == 0 and read_report(target)["basin_mean_concentration_ug_l"] == pytest.approx(expected)


def test_standard_error_ug_l_uses_effective_sample_size(tmp_path):
    """题面：「standard_error_ug_l = sqrt(sum(variance_term_h))」；"""
    document, truth = make_document()
    result, _, target, _ = invoke(tmp_path, document)
    expected = expected_values(document, truth)[2]
    assert result.returncode == 0 and read_report(target)["standard_error_ug_l"] == pytest.approx(expected)


def test_confidence_interval_uses_declared_z_score(tmp_path):
    """题面：「置信下界与上界分别为均值减、加 `z_score * standard_error_ug_l`；」；"""
    document, truth = make_document()
    result, _, target, _ = invoke(tmp_path, document)
    expected = expected_values(document, truth)
    report = read_report(target)
    assert result.returncode == 0 and report["ci_low_ug_l"] == pytest.approx(expected[3])
    assert report["ci_high_ug_l"] == pytest.approx(expected[4])


def test_design_confidence_is_validated_but_not_used_as_z(tmp_path):
    """题面：「`design_confidence` 参与字段校验，但不替代已锁定的 `z_score`；」；"""
    document, _ = make_document()
    changed = json.loads(json.dumps(document))
    changed["design_confidence"] = 0.9
    original, _, original_target, _ = invoke(tmp_path, document)
    altered, _, altered_target, _ = invoke(tmp_path, changed)
    assert original.returncode == 0 and altered.returncode == 0
    assert original_target.read_bytes() == altered_target.read_bytes()


def test_minimum_sample_boundary_is_accepted(tmp_path):
    """题面：「`samples` 有 12 到 24 份样本；」；"""
    document, _ = make_document(12, 3)
    result, _, target, _ = invoke(tmp_path, document)
    assert result.returncode == 0 and target.is_file()


def test_maximum_sample_boundary_is_accepted(tmp_path):
    """题面：「`samples` 有 12 到 24 份样本；」；"""
    document, _ = make_document(24, 4)
    result, _, target, _ = invoke(tmp_path, document)
    assert result.returncode == 0 and len(read_report(target)["zone_order"]) == 4


def test_eleven_samples_are_refused(tmp_path):
    """题面：「`samples` 有 12 到 24 份样本；」；"""
    document, _ = make_document()
    document["samples"] = document["samples"][:-1]
    result, _, target, _ = invoke(tmp_path, document)
    assert result.returncode == 1 and result.stderr == "运行拒收: BAD_INPUT\n" and not target.exists()


def test_extra_top_level_key_is_refused(tmp_path):
    """题面：「六个顶层名称组成唯一合法根形状；」"""
    document, _ = make_document()
    document["operator"] = "night"
    result, _, target, _ = invoke(tmp_path, document)
    assert result.returncode == 1 and result.stderr == "运行拒收: BAD_INPUT\n" and not target.exists()


def test_watershed_identifier_prefix_is_refused(tmp_path):
    """题面：「前者以 `RIVER-` 起头」"""
    document, _ = make_document()
    document["watershed_id"] = "HHS-2026-04"
    result, _, target, _ = invoke(tmp_path, document)
    assert result.returncode == 1 and result.stderr == "运行拒收: BAD_INPUT\n" and not target.exists()


def test_invalid_zone_identifier_is_refused(tmp_path):
    """题面：「`zone_id` 匹配 `[A-Z][A-Z0-9_-]{0,23}` 且不得重复；」；"""
    document, _ = make_document()
    document["zones"][0]["zone_id"] = "a"
    result, _, target, _ = invoke(tmp_path, document)
    assert result.returncode == 1 and result.stderr == "运行拒收: BAD_INPUT\n" and not target.exists()


def test_campaign_identifier_with_surrounding_space_is_refused(tmp_path):
    """题面：「`watershed_id` 和 `campaign_id` 的字符数应落在 2–48」"""
    document, _ = make_document()
    document["campaign_id"] = " campaign "
    result, _, target, _ = invoke(tmp_path, document)
    assert result.returncode == 1 and result.stderr == "运行拒收: BAD_INPUT\n" and not target.exists()
    minimum, _ = make_document()
    minimum["campaign_id"] = "AB"
    minimum_result, _, minimum_target, _ = invoke(tmp_path, minimum)
    assert minimum_result.returncode == 0 and json.loads(minimum_target.read_text(encoding="utf-8"))["watershed_id"] == minimum["watershed_id"]
    spaced, _ = make_document()
    spaced["campaign_id"] = "A B"
    spaced_result, _, spaced_target, _ = invoke(tmp_path, spaced)
    assert spaced_result.returncode == 0 and json.loads(spaced_target.read_text(encoding="utf-8"))["watershed_id"] == spaced["watershed_id"]
    spaced_bytes = spaced_target.read_bytes()
    control, _ = make_document()
    control["campaign_id"] = "A\tB"
    control_result, _, control_target, _ = invoke(tmp_path, control)
    assert control_result.returncode == 1 and control_result.stderr == "运行拒收: BAD_INPUT\n" and control_target.read_bytes() == spaced_bytes


def test_duplicate_sample_reference_is_refused(tmp_path):
    """题面：「`sample_ref` 的格式为“大写字母 + 2 至 15 个字母数字/下划线/连字符”，同批不可复现；」"""
    document, _ = make_document()
    document["samples"][1]["sample_ref"] = document["samples"][0]["sample_ref"]
    result, _, target, _ = invoke(tmp_path, document)
    assert result.returncode == 1 and result.stderr == "运行拒收: BAD_INPUT\n" and not target.exists()


def test_calibration_extra_key_is_refused(tmp_path):
    """题面：「`calibration` 仅有两个键：」；"""
    document, _ = make_document()
    document["calibration"]["method"] = "manual"
    result, _, target, _ = invoke(tmp_path, document)
    assert result.returncode == 1 and result.stderr == "运行拒收: BAD_INPUT\n" and not target.exists()


def test_missing_sample_key_is_refused(tmp_path):
    """题面：「每份样本只含四个键：」；"""
    document, _ = make_document()
    document["samples"][0].pop("dry_mass_g")
    result, _, target, _ = invoke(tmp_path, document)
    assert result.returncode == 1 and result.stderr == "运行拒收: BAD_INPUT\n" and not target.exists()


def test_duplicate_json_key_is_refused(tmp_path):
    """题面：「UTF-8 解码过程中若出现重复名称即拒收」。"""
    document, _ = make_document()
    raw = json.dumps(document, ensure_ascii=False, separators=(",", ":"))[:-1] + ',"watershed_id":"RIVER-X"}\n'
    result, _, target, _ = invoke(tmp_path, document, raw=raw)
    assert result.returncode == 1 and result.stderr == "运行拒收: BAD_INPUT\n" and not target.exists()


def test_boolean_aliquot_count_is_refused(tmp_path):
    """题面：「真/假值不属于这里可接受的整数或小数；」；"""
    document, _ = make_document()
    document["samples"][0]["aliquot_count"] = True
    result, _, target, _ = invoke(tmp_path, document)
    assert result.returncode == 1 and result.stderr == "运行拒收: BAD_INPUT\n" and not target.exists()


def test_concentration_length_mismatch_is_refused(tmp_path):
    """题面：「`concentrations_ug_l` 的条目数应与 `aliquot_count` 完全相等」"""
    document, _ = make_document()
    document["samples"][0]["concentrations_ug_l"].pop()
    result, _, target, _ = invoke(tmp_path, document)
    assert result.returncode == 1 and result.stderr == "运行拒收: BAD_INPUT\n" and not target.exists()


def test_boolean_concentration_is_refused(tmp_path):
    """题面：「真/假值不属于这里可接受的整数或小数」。"""
    document, _ = make_document()
    document["samples"][0]["concentrations_ug_l"][0] = False
    result, _, target, _ = invoke(tmp_path, document)
    assert result.returncode == 1 and result.stderr == "运行拒收: BAD_INPUT\n" and not target.exists()


def test_zone_sample_count_must_sum_to_sample_total(tmp_path):
    """题面：「每区的 `sample_count` 是实际样本数」；"""
    document, _ = make_document()
    document["zones"][0]["sample_count"] += 1
    result, _, target, _ = invoke(tmp_path, document)
    assert result.returncode == 1 and result.stderr == "运行拒收: BAD_INPUT\n" and not target.exists()


def test_inconsistent_aggregate_reports_bad_model(tmp_path):
    """题面：「闭合失败时不生成结果，退出 1，stderr 固定 `运行拒收: BAD_MODEL`；」；"""
    document, _ = make_document()
    document["zones"][0]["analyte_total"] += 9.0
    result, _, target, _ = invoke(tmp_path, document)
    assert result.returncode == 1 and result.stderr == "运行拒收: BAD_MODEL\n" and not target.exists()


def test_wrong_tie_break_is_refused(tmp_path):
    """题面：「`tie_break` 固定为字符串 `lexicographic`；」；"""
    document, _ = make_document()
    document["calibration"]["tie_break"] = "arrival"
    result, _, target, _ = invoke(tmp_path, document)
    assert result.returncode == 1 and result.stderr == "运行拒收: BAD_INPUT\n" and not target.exists()


def test_ledger_has_exact_four_key_line(tmp_path):
    """题面：「每行四键，紧凑 JSON，行尾单个 LF：」；"""
    document, _ = make_document()
    result, _, target, ledger = invoke(tmp_path, document, with_ledger=True)
    line = ledger.read_text(encoding="utf-8")
    record = json.loads(line)
    assert result.returncode == 0 and line.endswith("\n") and list(record) == [
        "watershed_id",
        "input_sha256",
        "result_sha256",
        "basin_mean_concentration_ug_l",
    ]
    assert record["result_sha256"] == __import__("hashlib").sha256(target.read_bytes()).hexdigest()


def test_repeated_success_keeps_ledger_bytes(tmp_path):
    """题面：「相同输入成功重跑时，台账一个字节都不变；」；"""
    document, _ = make_document()
    first, _, _, ledger = invoke(tmp_path, document, with_ledger=True)
    before = ledger.read_bytes()
    second = subprocess.run(
        [sys.executable, ENTRY, str(tmp_path / "watershed.json"), str(tmp_path / "report.json"), str(ledger)],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert first.returncode == 0 and second.returncode == 0 and ledger.read_bytes() == before


def test_distinct_input_appends_in_success_order(tmp_path):
    """题面：「复核单摘要变化时，流水在其末尾登记新行；」；"""
    document, _ = make_document()
    first, _, _, ledger = invoke(tmp_path, document, with_ledger=True)
    changed = json.loads(json.dumps(document))
    changed["watershed_id"] = "RIVER-BASIN-05"
    changed["campaign_id"] = "CAMPAIGN-FOLLOWUP"
    changed["samples"][0]["concentrations_ug_l"][0] += 1
    recalibrate(changed)
    second, _, _, _ = invoke(tmp_path, changed, with_ledger=True)
    assert first.returncode == 0 and second.returncode == 0
    assert [json.loads(line)["watershed_id"] for line in ledger.read_text(encoding="utf-8").splitlines()] == [
        document["watershed_id"],
        changed["watershed_id"],
    ]


def test_bad_ledger_is_refused_without_replacing_report(tmp_path):
    """题面：「流水重放或格式异常统一以 `BAD_LEDGER` 拒绝；」；"""
    document, _ = make_document()
    result, source, target, ledger = invoke(tmp_path, document)
    target.write_bytes(b"old bytes\n")
    ledger.write_bytes(b"{}\n")
    attempted = subprocess.run(
        [sys.executable, ENTRY, str(source), str(target), str(ledger)],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0 and attempted.returncode == 1
    assert attempted.stderr == "运行拒收: BAD_LEDGER\n" and target.read_bytes() == b"old bytes\n"


def test_ledger_parent_failure_preserves_old_result(tmp_path):
    """题面：「复核单和台账应同时发布或同时保持旧字节；」；"""
    document, _ = make_document()
    result, source, target, _ = invoke(tmp_path, document)
    target.write_bytes(b"previous report\n")
    absent_ledger = tmp_path / "nested" / "delivery.jsonl"
    attempted = subprocess.run(
        [sys.executable, ENTRY, str(source), str(target), str(absent_ledger)],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0 and attempted.returncode == 1
    assert attempted.stderr == "运行拒收: IO_FAILURE\n" and target.read_bytes() == b"previous report\n"


def test_input_file_is_not_modified(tmp_path):
    """题面：「输入材料本身不可改写；」；"""
    document, _ = make_document()
    result, source, _, _ = invoke(tmp_path, document, with_ledger=True)
    assert result.returncode == 0 and source.read_bytes() == (
        json.dumps(document, ensure_ascii=False, separators=(",", ":")) + "\n"
    ).encode("utf-8")


def test_success_leaves_no_staging_files(tmp_path):
    """题面：「成功发布后不得残留临时文件」；"""
    document, _ = make_document()
    result, _, _, _ = invoke(tmp_path, document, with_ledger=True)
    names = sorted(path.name for path in tmp_path.iterdir())
    assert result.returncode == 0 and names == ["delivery.jsonl", "report.json", "watershed.json"]


def test_repeat_without_ledger_is_byte_stable(tmp_path):
    """题面：「结果文件只允许一个 JSON 对象和结尾换行」。"""
    document, _ = make_document()
    first, _, target, _ = invoke(tmp_path, document)
    before = target.read_bytes()
    second = subprocess.run(
        [sys.executable, ENTRY, str(tmp_path / "watershed.json"), str(target)],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert first.returncode == 0 and second.returncode == 0 and target.read_bytes() == before


def test_missing_input_is_bad_input(tmp_path):
    """题面：「读取或校验材料失败，退出码 1，标准错误仅含 `运行拒收: BAD_INPUT`；」；"""
    target = tmp_path / "report.json"
    result = subprocess.run(
        [sys.executable, ENTRY, str(tmp_path / "absent.json"), str(target)],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 1 and result.stderr == "运行拒收: BAD_INPUT\n" and not target.exists()


def test_missing_result_parent_reports_io_failure(tmp_path):
    """题面：「读写位置异常时终端仅显示 `运行拒收: IO_FAILURE`；」"""
    document, _ = make_document()
    source = write_document(tmp_path / "batch.json", document)
    target = tmp_path / "absent" / "report.json"
    result = subprocess.run(
        [sys.executable, ENTRY, str(source), str(target)],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 1 and result.stderr == "运行拒收: IO_FAILURE\n" and not target.exists()


def test_failure_does_not_emit_traceback(tmp_path):
    """题面：「异常路径禁止向终端展示 Python `Traceback`」；"""
    document, _ = make_document()
    document["calibration"]["z_score"] = -1
    result, _, target, _ = invoke(tmp_path, document)
    assert result.returncode == 1 and "Traceback" not in result.stderr and not target.exists()


def test_maximum_group_shape_finishes_within_limit(tmp_path):
    """题面：「所有样本必须完整划分到各来源区」；"""
    document, _ = make_document(24, 4)
    result, _, target, _ = invoke(tmp_path, document)
    assert result.returncode == 0 and target.is_file()


def test_three_source_materials_produce_report(tmp_path):
    """题面：「实验室合并批次时丢失了每份泥沙样本所属的集水区编号，只留下样本干重、各等分样浓度和三个来源区的汇总值；」"""
    document, _ = make_document()
    result, _, target, _ = invoke(tmp_path, document)
    report = read_report(target)
    assert result.returncode == 0 and report["watershed_id"] == document["watershed_id"]
    assert report["status"] == "complete" and report["zone_order"] == ["A", "B", "C"]


def test_five_zones_are_accepted(tmp_path):
    """题面：「`zones` 有 3 到 5 个来源区，每个来源区的成员固定为：」；"""
    document, _ = make_document(20, 5)
    result, _, target, _ = invoke(tmp_path, document)
    assert result.returncode == 0 and read_report(target)["zone_order"] == ["A", "B", "C", "D", "E"]


def test_zone_order_follows_identifier_order(tmp_path):
    """题面：「`zone_order` 按 `zone_id` 升序；」；"""
    document, _ = make_document()
    result, _, target, _ = invoke(tmp_path, document)
    report = read_report(target)
    assert result.returncode == 0 and report["zone_order"] == sorted(report["zone_order"])


def test_lexicographic_normalization_is_stable(tmp_path):
    """题面：「取字典序最小者」。"""
    document, _ = make_document()
    shuffled, _ = make_document(shuffle=True)
    first, _, first_target, _ = invoke(tmp_path, document)
    write_document(tmp_path / "watershed.json", shuffled)
    second = subprocess.run(
        [sys.executable, ENTRY, str(tmp_path / "watershed.json"), str(tmp_path / "shuffled.json")],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert first.returncode == 0 and second.returncode == 0
    assert first_target.read_bytes() == (tmp_path / "shuffled.json").read_bytes()


def test_zone_arrays_follow_identifier_order(tmp_path):
    """题面：「`zone_samples` 和 `zone_mean_concentration_ug_l` 的元素位置与 `zone_order` 同步」"""
    document, _ = make_document()
    result, _, target, _ = invoke(tmp_path, document)
    report = read_report(target)
    count = len(report["zone_order"])
    assert result.returncode == 0 and len(report["zone_samples"]) == count
    assert len(report["zone_mean_concentration_ug_l"]) == count


def test_member_references_are_sorted_within_zone(tmp_path):
    """题面：「区内 `sample_ref` 仍按升序展开」。"""
    document, _ = make_document()
    result, _, target, _ = invoke(tmp_path, document)
    report = read_report(target)
    assert result.returncode == 0 and all(row == sorted(row) for row in report["zone_samples"])


def test_ledger_is_append_only_jsonl(tmp_path):
    """题面：「台账采用只追加的 UTF-8 JSONL；」；"""
    document, _ = make_document()
    first, _, _, ledger = invoke(tmp_path, document, with_ledger=True)
    before = ledger.read_bytes()
    second, _, _, _ = invoke(tmp_path, document, with_ledger=True)
    assert first.returncode == 0 and second.returncode == 0 and ledger.read_bytes() == before
    assert before.decode("utf-8").endswith("\n")
