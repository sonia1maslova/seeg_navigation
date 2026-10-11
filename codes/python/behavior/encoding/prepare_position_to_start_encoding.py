"""Create isolated current-position-to-start bearing tables; never edit sources."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path

TAU = 2 * math.pi


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def exclusive_or_identical(path: Path, payload: bytes) -> None:
    if path.exists():
        if path.read_bytes() != payload:
            raise FileExistsError(f"Refusing to overwrite different existing file: {path}")
        return
    with path.open("xb") as handle:
        handle.write(payload)


def build_table(source: Path, target: Path, keys: list[str], angles: list[str]) -> dict:
    source_bytes = source.read_bytes()
    reader = csv.DictReader(io.StringIO(source_bytes.decode("utf-8-sig")))
    original_columns = reader.fieldnames
    if not original_columns or not set(keys + angles).issubset(original_columns):
        raise ValueError(f"Unexpected schema in {source}")
    rows = list(reader)
    old_columns = [name + "_start_to_position" for name in angles]
    if set(old_columns).intersection(original_columns):
        raise ValueError("Input already appears to have been recoded")
    seen = set()
    counts = {name: 0 for name in angles}
    for row in rows:
        key = tuple(row[name] for name in keys)
        if key in seen:
            raise ValueError(f"Duplicate key {key} in {source}")
        seen.add(key)
        for name in angles:
            original = row[name]
            row[name + "_start_to_position"] = original
            if original.strip().lower() in {"", "na", "nan"}:
                continue
            theta = float(original)
            if not math.isfinite(theta) or not 0 <= theta < TAU:
                raise ValueError(f"Invalid angle {name}={original}, key={key}")
            reversed_theta = (theta + math.pi) % TAU
            error = abs(math.atan2(math.sin(reversed_theta - theta - math.pi),
                                   math.cos(reversed_theta - theta - math.pi)))
            if error > 2e-12:
                raise AssertionError("Bearing inversion failed")
            row[name] = repr(reversed_theta)
            counts[name] += 1
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=original_columns + old_columns,
                            lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    payload = output.getvalue().encode("utf-8")
    reread = list(csv.DictReader(io.StringIO(payload.decode())))
    original_rows = list(csv.DictReader(io.StringIO(source_bytes.decode("utf-8-sig"))))
    if len(original_rows) != len(reread):
        raise AssertionError("Row count changed")
    for before, after in zip(original_rows, reread):
        for col in original_columns:
            if col not in angles and before[col] != after[col]:
                raise AssertionError(f"Unrequested change: {col}")
        for col in angles:
            if before[col] != after[col + "_start_to_position"]:
                raise AssertionError("Original angle not preserved")
    exclusive_or_identical(target, payload)
    if source.read_bytes() != source_bytes:
        raise AssertionError("Original encoding file changed during generation")
    return {"source": str(source), "target": str(target), "rows": len(rows),
            "unique_key": keys, "recoded_angles": counts,
            "source_sha256": hashlib.sha256(source_bytes).hexdigest(),
            "target_sha256": sha256(target), "other_columns_unchanged": True}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    args = parser.parse_args()
    root = args.project_root.resolve()
    encoding = root / "data" / "encoding"
    if not encoding.is_dir():
        raise FileNotFoundError(encoding)
    records = []
    for stem, keys, angles in [
        ("corridor_subsegi_encoding_table", ["Trial", "Corridor", "Half"],
         ["Angle2Xaxis", "Angle2orientation"]),
        ("turn_encoding_table", ["Trial", "Turn"], ["Angle2Xaxis"]),
    ]:
        records.append(build_table(encoding / f"{stem}.csv",
                                  encoding / f"{stem}_position_to_start.csv", keys, angles))
    manifest = {"branch": "bearing_position_to_start", "direction": "current_position_to_start",
                "transform": "(original_angle + pi) modulo (2*pi)",
                "unit": "radian", "range": "[0, 2*pi)",
                "source_direction": "start_to_current_position", "tables": records}
    path = encoding / "bearing_position_to_start_encoding_manifest.json"
    exclusive_or_identical(path, (json.dumps(manifest, ensure_ascii=False, indent=2) + "\n").encode())
    print(json.dumps(manifest, ensure_ascii=False))


if __name__ == "__main__":
    main()
