"""Run all shared plotting scripts in batch mode.

Usage from the project root:

    python Codes/run_all_figures.py
    python Codes/run_all_figures.py --only main
    python Codes/run_all_figures.py --only sup
    python Codes/run_all_figures.py --stop-on-error
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from datetime import datetime
from pathlib import Path

from project_paths import CODES_DIR, PROJECT_ROOT, LOG_DIR


def script_sort_key(path: Path):
    name = path.stem
    is_sup = name.lower().startswith("sup_")
    return (1 if is_sup else 0, name.lower())


def collect_scripts(only: str):
    scripts = []
    if only in {"all", "main"}:
        scripts.extend(CODES_DIR.glob("Figure*.py"))
    if only in {"all", "sup"}:
        scripts.extend(CODES_DIR.glob("Sup_fig*.py"))
    return sorted(scripts, key=script_sort_key)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run all figure plotting scripts.")
    parser.add_argument("--only", choices=["all", "main", "sup"], default="all")
    parser.add_argument("--stop-on-error", action="store_true")
    args = parser.parse_args()

    LOG_DIR.mkdir(parents=True, exist_ok=True)
    scripts = collect_scripts(args.only)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = LOG_DIR / f"run_all_figures_{timestamp}.log"

    print(f"Project root: {PROJECT_ROOT}")
    print(f"Scripts     : {len(scripts)}")
    print(f"Log file    : {log_file}")
    print("=" * 80)

    env = os.environ.copy()
    env["MPLBACKEND"] = "Agg"
    env["PYTHONPATH"] = str(CODES_DIR) + os.pathsep + env.get("PYTHONPATH", "")

    success = []
    failed = []
    with log_file.open("w", encoding="utf-8") as log:
        log.write(f"Project root: {PROJECT_ROOT}\n")
        log.write(f"Python      : {sys.executable}\n")
        log.write(f"Scripts     : {len(scripts)}\n")
        log.write("=" * 80 + "\n")

        for script in scripts:
            rel = script.relative_to(PROJECT_ROOT)
            print(f"Running: {rel}", flush=True)
            log.write(f"\n>>> Running: {rel}\n")
            log.flush()

            proc = subprocess.run(
                [sys.executable, str(script)],
                cwd=str(PROJECT_ROOT),
                env=env,
                text=True,
                capture_output=True,
            )

            if proc.stdout:
                log.write("--- stdout ---\n")
                log.write(proc.stdout)
                if not proc.stdout.endswith("\n"):
                    log.write("\n")
            if proc.stderr:
                log.write("--- stderr ---\n")
                log.write(proc.stderr)
                if not proc.stderr.endswith("\n"):
                    log.write("\n")
            log.write(f"--- return code: {proc.returncode} ---\n")
            log.flush()

            if proc.returncode == 0:
                success.append(rel)
                print(f"  OK: {rel}")
            else:
                failed.append(rel)
                print(f"  FAILED: {rel}")
                if args.stop_on_error:
                    break

    print("=" * 80)
    print("SUMMARY")
    print(f"  Success: {len(success)}")
    print(f"  Failed : {len(failed)}")
    if failed:
        print("Failed scripts:")
        for item in failed:
            print(f"  - {item}")
    print(f"Log file: {log_file}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
