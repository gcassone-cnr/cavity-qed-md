#!/usr/bin/env python3
"""Summarize eT force-call timing and SCF iterations from driver logs."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import re
import statistics


RECORD = re.compile(
    r"force evaluation\s+(?P<step>\d+):.*?"
    r"SCF iterations\s*=\s*(?P<iterations>\d+|None),\s*"
    r"restart\s*=\s*(?P<restart>True|False),\s*"
    r"wall\s*=\s*(?P<wall>[0-9.eE+-]+)\s*s"
)


def percentile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    position = fraction * (len(ordered) - 1)
    low = math.floor(position)
    high = math.ceil(position)
    weight = position - low
    return ordered[low] * (1.0 - weight) + ordered[high] * weight


def parse_log(path: Path) -> list[dict[str, float | int | bool | None]]:
    records: list[dict[str, float | int | bool | None]] = []
    for match in RECORD.finditer(path.read_text(encoding="utf-8", errors="replace")):
        iteration_text = match.group("iterations")
        records.append(
            {
                "step": int(match.group("step")),
                "scf_iterations": (
                    None if iteration_text == "None" else int(iteration_text)
                ),
                "restart_loaded": match.group("restart") == "True",
                "wall_seconds": float(match.group("wall")),
            }
        )
    return records


def summarize(
    path: Path, records: list[dict[str, float | int | bool | None]], discard: int,
    timestep_fs: float,
) -> dict[str, object]:
    retained = records[discard:]
    if not retained:
        raise ValueError(
            f"{path}: no force records remain after discarding {discard}"
        )
    walls = [float(record["wall_seconds"]) for record in retained]
    iterations = [
        int(record["scf_iterations"])
        for record in retained
        if record["scf_iterations"] is not None
    ]
    median_wall = statistics.median(walls)
    mean_wall = statistics.fmean(walls)
    return {
        "path": str(path.resolve()),
        "records_total": len(records),
        "records_used": len(retained),
        "discarded": discard,
        "wall_seconds_mean": mean_wall,
        "wall_seconds_median": median_wall,
        "wall_seconds_p10": percentile(walls, 0.10),
        "wall_seconds_p90": percentile(walls, 0.90),
        "scf_iterations_mean": statistics.fmean(iterations) if iterations else None,
        "scf_iterations_median": statistics.median(iterations) if iterations else None,
        "restart_loaded_fraction": statistics.fmean(
            float(bool(record["restart_loaded"])) for record in retained
        ),
        "estimated_fs_per_day_from_mean": timestep_fs * 86400.0 / mean_wall,
        "estimated_fs_per_day_from_median": timestep_fs * 86400.0 / median_wall,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("logs", nargs="+", help="et_driver_*.log files")
    parser.add_argument(
        "--timestep-fs", type=float, required=True, help="i-PI MD timestep in fs"
    )
    parser.add_argument(
        "--discard", type=int, default=1,
        help="discard this many initial force calls from each log (default: 1)",
    )
    parser.add_argument("--json", dest="json_path", help="optional JSON output")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.timestep_fs <= 0.0:
        raise SystemExit("--timestep-fs must be positive")
    if args.discard < 0:
        raise SystemExit("--discard must be non-negative")
    summaries = []
    for name in args.logs:
        path = Path(name)
        records = parse_log(path)
        if not records:
            raise SystemExit(f"{path}: no driver force-evaluation records found")
        summaries.append(summarize(path, records, args.discard, args.timestep_fs))

    for item in summaries:
        print(item["path"])
        print(
            f"  calls used: {item['records_used']}/{item['records_total']}  "
            f"wall median: {item['wall_seconds_median']:.2f} s  "
            f"p10-p90: {item['wall_seconds_p10']:.2f}-"
            f"{item['wall_seconds_p90']:.2f} s"
        )
        print(
            f"  SCF iterations median: {item['scf_iterations_median']}  "
            f"restart loaded: {100.0 * item['restart_loaded_fraction']:.1f}%  "
            f"estimated throughput: "
            f"{item['estimated_fs_per_day_from_median']:.2f} fs/day"
        )

    if args.json_path:
        output = Path(args.json_path)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(
            json.dumps(
                {"timestep_fs": args.timestep_fs, "logs": summaries},
                indent=2,
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )
        print(f"JSON: {output.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
