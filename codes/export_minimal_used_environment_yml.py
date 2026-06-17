#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
export_minimal_used_environment_yml.py
=====================================

Purpose
-------
Create a MINIMAL Conda environment YAML file that contains only the third-party
packages actually imported by the current plotting code, instead of exporting
the whole Conda environment.

Recommended location:
    D:/下载/project_shared_plotting_updated/project_shared_plotting/codes/export_minimal_used_environment_yml.py

Run in your current Conda environment:
    conda activate DiversityMaping
    cd /d "D:/下载/project_shared_plotting_updated/project_shared_plotting"
    python ".\\codes\\export_minimal_used_environment_yml.py"

Main output:
    environment_minimal_used.yml

Other outputs:
    environment_minimal_used_requirements.txt
    environment_minimal_used_packages.csv
    environment_minimal_imports_detected.csv
    environment_minimal_missing_imports.txt

Notes
-----
1. This script scans Python import statements in the code folder.
2. It excludes standard-library modules and local project modules.
3. It checks package versions in your CURRENT active environment.
4. It writes only directly used third-party packages. Conda will resolve their
   dependencies when recreating the environment.
"""

from __future__ import annotations

import argparse
import ast
import csv
import json
import os
import platform
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

try:
    from importlib import metadata as importlib_metadata
except Exception:  # pragma: no cover
    import importlib_metadata  # type: ignore


# Import module name -> likely distribution/package names.
# This helps when import name and install name differ.
MODULE_TO_PACKAGE_HINTS: Dict[str, List[str]] = {
    "PIL": ["pillow", "Pillow"],
    "cv2": ["opencv-python", "opencv"],
    "sklearn": ["scikit-learn"],
    "mpl_toolkits": ["matplotlib"],
    "yaml": ["pyyaml", "PyYAML"],
    "dateutil": ["python-dateutil"],
    "statsmodels": ["statsmodels"],
    "seaborn": ["seaborn"],
    "matplotlib": ["matplotlib"],
    "numpy": ["numpy"],
    "pandas": ["pandas"],
    "scipy": ["scipy"],
    "nilearn": ["nilearn"],
    "networkx": ["networkx"],
    "tqdm": ["tqdm"],
    "h5py": ["h5py"],
    "pymer4": ["pymer4"],
    "IPython": ["ipython"],
}


# Local modules in this shared plotting project.
LOCAL_MODULE_NAMES = {
    "project_paths",
    "functions",
    "figure_helpers",
    "my_plot_setup",
}


@dataclass
class ImportRecord:
    file: str
    line: int
    module: str
    statement: str


def normalize_name(name: str) -> str:
    """Normalize package names for comparison."""
    return re.sub(r"[-_.]+", "-", name).lower()


def infer_project_root() -> Path:
    """
    Infer project root if this script is placed under project_root/codes or project_root/Codes.
    """
    this_file = Path(__file__).resolve()
    if this_file.parent.name.lower() in {"codes", "code"}:
        return this_file.parent.parent
    return Path.cwd().resolve()


def default_code_dir(project_root: Path) -> Path:
    for name in ("codes", "Codes", "code", "Code"):
        candidate = project_root / name
        if candidate.exists():
            return candidate.resolve()
    return (project_root / "codes").resolve()


def get_stdlib_modules() -> Set[str]:
    stdlib = set(sys.builtin_module_names)
    if hasattr(sys, "stdlib_module_names"):
        stdlib.update(sys.stdlib_module_names)
    else:
        stdlib.update(
            {
                "abc", "argparse", "ast", "collections", "copy", "csv", "datetime",
                "functools", "glob", "gzip", "importlib", "inspect", "itertools",
                "json", "logging", "math", "os", "pathlib", "pickle", "platform",
                "re", "shutil", "site", "statistics", "subprocess", "sys",
                "tempfile", "time", "typing", "warnings", "zipfile",
            }
        )
    return stdlib


def top_level_name(module_name: str) -> str:
    return module_name.split(".")[0]


def iter_python_files(code_dir: Path) -> Iterable[Path]:
    skip_dirs = {"__pycache__", ".git", ".ipynb_checkpoints", "logs", "environment"}
    for path in sorted(code_dir.rglob("*.py")):
        if any(part in skip_dirs for part in path.parts):
            continue
        # Avoid scanning this exporter itself and older environment exporters if desired.
        if path.name in {
            "export_conda_environment_yml.py",
            "export_minimal_used_environment_yml.py",
            "get_used_package_versions.py",
        }:
            continue
        yield path


def parse_imports(py_file: Path, project_root: Path) -> List[ImportRecord]:
    try:
        source = py_file.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        source = py_file.read_text(encoding="utf-8-sig")

    rel_file = str(py_file.relative_to(project_root))
    records: List[ImportRecord] = []

    try:
        tree = ast.parse(source, filename=str(py_file))
    except SyntaxError as exc:
        records.append(
            ImportRecord(
                file=rel_file,
                line=getattr(exc, "lineno", 0) or 0,
                module="<SYNTAX_ERROR>",
                statement=f"SyntaxError: {exc}",
            )
        )
        return records

    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                records.append(
                    ImportRecord(
                        file=rel_file,
                        line=node.lineno,
                        module=alias.name,
                        statement=f"import {alias.name}",
                    )
                )

        elif isinstance(node, ast.ImportFrom):
            # Relative imports, e.g. from .ExtractFeature import seg_name, are local.
            if node.level and node.level > 0:
                continue
            if node.module:
                records.append(
                    ImportRecord(
                        file=rel_file,
                        line=node.lineno,
                        module=node.module,
                        statement=f"from {node.module} import ...",
                    )
                )

    return records


def get_local_top_level_modules(code_dir: Path) -> Set[str]:
    local = set(LOCAL_MODULE_NAMES)

    for path in code_dir.glob("*.py"):
        if path.name.startswith("_"):
            continue
        local.add(path.stem)

    for path in code_dir.iterdir():
        if path.is_dir():
            local.add(path.name)
            if (path / "__init__.py").exists():
                local.add(path.name)

    functions_dir = code_dir / "functions"
    if functions_dir.exists():
        for path in functions_dir.glob("*.py"):
            if path.name != "__init__.py":
                local.add(path.stem)

    return local


def is_third_party_module(module: str, code_dir: Path) -> bool:
    top = top_level_name(module)
    if module == "<SYNTAX_ERROR>":
        return False
    if top in get_stdlib_modules():
        return False
    if top in get_local_top_level_modules(code_dir):
        return False
    if top.startswith("_"):
        return False
    return True


def build_module_to_distribution_map() -> Dict[str, List[str]]:
    mapping: Dict[str, List[str]] = {}

    try:
        pkg_map = importlib_metadata.packages_distributions()
        for module_name, dist_names in pkg_map.items():
            mapping[module_name] = sorted(set(dist_names), key=str.lower)
    except Exception:
        pass

    for module_name, dist_names in MODULE_TO_PACKAGE_HINTS.items():
        existing = mapping.get(module_name, [])
        mapping[module_name] = sorted(set(existing + dist_names), key=str.lower)

    return mapping


def get_installed_distribution_version(dist_name: str) -> Optional[str]:
    try:
        return importlib_metadata.version(dist_name)
    except Exception:
        return None


def run_command(cmd: List[str]) -> Tuple[int, str, str]:
    result = subprocess.run(
        cmd,
        text=True,
        capture_output=True,
        check=False,
        encoding="utf-8",
        errors="replace",
    )
    return result.returncode, result.stdout, result.stderr


def load_conda_list() -> Dict[str, Dict[str, str]]:
    """
    Return conda list package info keyed by normalized package name.
    """
    code, stdout, stderr = run_command(["conda", "list", "--json"])
    if code != 0:
        print("WARNING: failed to run `conda list --json`.")
        print(stderr.strip())
        return {}

    try:
        records = json.loads(stdout)
    except json.JSONDecodeError:
        print("WARNING: failed to parse `conda list --json` output.")
        return {}

    out: Dict[str, Dict[str, str]] = {}
    for rec in records:
        name = rec.get("name", "")
        if not name:
            continue
        out[normalize_name(name)] = {
            "name": name,
            "version": rec.get("version", ""),
            "channel": rec.get("channel", ""),
            "build_string": rec.get("build_string", ""),
        }
    return out


def resolve_package_for_module(
    module_name: str,
    module_to_dist: Dict[str, List[str]],
    conda_pkgs: Dict[str, Dict[str, str]],
) -> Dict[str, str]:
    """
    Resolve import module to package name/version/source.

    source:
        conda  -> available in conda list, channel not pypi
        pip    -> available in conda list with channel pypi, or dist metadata only
        missing
    """
    top = top_level_name(module_name)

    candidates: List[str] = []
    candidates.extend(module_to_dist.get(top, []))
    candidates.extend(MODULE_TO_PACKAGE_HINTS.get(top, []))
    candidates.append(top)

    # Deduplicate while preserving order.
    seen = set()
    candidates_unique = []
    for c in candidates:
        key = normalize_name(c)
        if key not in seen:
            candidates_unique.append(c)
            seen.add(key)

    for candidate in candidates_unique:
        key = normalize_name(candidate)
        if key in conda_pkgs:
            rec = conda_pkgs[key]
            channel = rec.get("channel", "")
            source = "pip" if channel == "pypi" else "conda"
            return {
                "import_module": top,
                "package_name": rec["name"],
                "version": rec["version"],
                "source": source,
                "channel": channel,
                "status": "installed",
            }

    # Fall back to importlib.metadata. This usually means pip metadata exists
    # but conda list matching failed.
    for candidate in candidates_unique:
        version = get_installed_distribution_version(candidate)
        if version:
            return {
                "import_module": top,
                "package_name": candidate,
                "version": version,
                "source": "pip",
                "channel": "",
                "status": "installed_metadata_only",
            }

    return {
        "import_module": top,
        "package_name": top,
        "version": "",
        "source": "",
        "channel": "",
        "status": "missing",
    }


def write_csv(path: Path, rows: List[Dict[str, str]], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def yaml_quote(value: str) -> str:
    # Simple safe YAML quoting for environment name.
    if re.match(r"^[A-Za-z0-9_.-]+$", value):
        return value
    return '"' + value.replace('"', '\\"') + '"'


def make_environment_yml(
    env_name: str,
    python_version: str,
    conda_deps: List[Tuple[str, str]],
    pip_deps: List[Tuple[str, str]],
    channels: List[str],
) -> str:
    lines: List[str] = []
    lines.append(f"name: {yaml_quote(env_name)}")
    lines.append("channels:")
    for ch in channels:
        lines.append(f"  - {ch}")
    lines.append("dependencies:")
    lines.append(f"  - python={python_version}")

    for name, version in conda_deps:
        if version:
            lines.append(f"  - {name}={version}")
        else:
            lines.append(f"  - {name}")

    if pip_deps:
        # Ensure pip exists before the pip subsection.
        if not any(normalize_name(name) == "pip" for name, _ in conda_deps):
            lines.append("  - pip")
        lines.append("  - pip:")
        for name, version in pip_deps:
            if version:
                lines.append(f"      - {name}=={version}")
            else:
                lines.append(f"      - {name}")

    return "\n".join(lines) + "\n"


def main() -> int:
    project_root_default = infer_project_root()

    parser = argparse.ArgumentParser(
        description="Export a minimal environment_minimal_used.yml based on imports used by project code."
    )
    parser.add_argument(
        "--project-root",
        type=Path,
        default=project_root_default,
        help="Project root. Default: inferred from this script location.",
    )
    parser.add_argument(
        "--code-dir",
        type=Path,
        default=None,
        help="Code directory. Default: project_root/codes or project_root/Codes.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Output directory. Default: project root.",
    )
    parser.add_argument(
        "--env-name",
        type=str,
        default=None,
        help="Name written into environment_minimal_used.yml. Default: current_env_minimal.",
    )
    parser.add_argument(
        "--channels",
        type=str,
        default="conda-forge,defaults",
        help="Comma-separated channels for the YAML. Default: conda-forge,defaults.",
    )
    parser.add_argument(
        "--no-version-pin",
        action="store_true",
        help="Do not pin package versions. Only write package names.",
    )
    args = parser.parse_args()

    project_root = args.project_root.resolve()
    code_dir = args.code_dir.resolve() if args.code_dir else default_code_dir(project_root)
    out_dir = args.out_dir.resolve() if args.out_dir else project_root
    out_dir.mkdir(parents=True, exist_ok=True)

    if not code_dir.exists():
        print(f"ERROR: code directory not found: {code_dir}")
        return 1

    current_env = os.environ.get("CONDA_DEFAULT_ENV", "plotting_env")
    env_name = args.env_name or f"{current_env}_minimal"
    py_version = f"{sys.version_info.major}.{sys.version_info.minor}"
    channels = [x.strip() for x in args.channels.split(",") if x.strip()]

    print(f"Project root       : {project_root}")
    print(f"Code folder        : {code_dir}")
    print(f"Output folder      : {out_dir}")
    print(f"Current conda env  : {current_env}")
    print(f"Python executable  : {sys.executable}")
    print("-" * 80)

    all_imports: List[ImportRecord] = []
    for py_file in iter_python_files(code_dir):
        all_imports.extend(parse_imports(py_file, project_root))

    import_rows: List[Dict[str, str]] = []
    third_party_modules: Set[str] = set()

    for rec in all_imports:
        top = top_level_name(rec.module)
        is_third = is_third_party_module(rec.module, code_dir)
        if is_third:
            third_party_modules.add(top)
        import_rows.append(
            {
                "file": rec.file,
                "line": str(rec.line),
                "module": rec.module,
                "top_level_module": top,
                "third_party": str(is_third),
                "statement": rec.statement,
            }
        )

    conda_pkgs = load_conda_list()
    module_to_dist = build_module_to_distribution_map()

    resolved: List[Dict[str, str]] = []
    for module in sorted(third_party_modules, key=str.lower):
        resolved.append(resolve_package_for_module(module, module_to_dist, conda_pkgs))

    # Deduplicate packages, prefer conda records over pip metadata if duplicate.
    package_map: Dict[str, Dict[str, str]] = {}
    for rec in resolved:
        key = normalize_name(rec["package_name"])
        old = package_map.get(key)
        if old is None:
            package_map[key] = rec
        else:
            if old.get("source") == "pip" and rec.get("source") == "conda":
                package_map[key] = rec

    package_rows = sorted(package_map.values(), key=lambda x: normalize_name(x["package_name"]))

    conda_deps: List[Tuple[str, str]] = []
    pip_deps: List[Tuple[str, str]] = []
    missing: List[str] = []

    for rec in package_rows:
        name = rec["package_name"]
        version = "" if args.no_version_pin else rec["version"]
        if rec["status"].startswith("installed"):
            if rec["source"] == "conda":
                conda_deps.append((name, version))
            elif rec["source"] == "pip":
                pip_deps.append((name, version))
            else:
                pip_deps.append((name, version))
        else:
            missing.append(f"{rec['import_module']}\t{name}\t{rec['status']}")

    # Sort and remove accidental python/pip duplication in conda deps.
    conda_deps = sorted(
        [(n, v) for n, v in conda_deps if normalize_name(n) not in {"python"}],
        key=lambda x: normalize_name(x[0]),
    )
    pip_deps = sorted(pip_deps, key=lambda x: normalize_name(x[0]))

    yml_text = make_environment_yml(
        env_name=env_name,
        python_version=py_version,
        conda_deps=conda_deps,
        pip_deps=pip_deps,
        channels=channels,
    )

    env_yml_path = out_dir / "environment_minimal_used.yml"
    env_yml_path.write_text(yml_text, encoding="utf-8")

    req_lines = []
    for name, version in conda_deps + pip_deps:
        req_lines.append(f"{name}=={version}" if version else name)
    (out_dir / "environment_minimal_used_requirements.txt").write_text(
        "\n".join(req_lines) + ("\n" if req_lines else ""),
        encoding="utf-8",
    )

    write_csv(
        out_dir / "environment_minimal_used_packages.csv",
        package_rows,
        ["import_module", "package_name", "version", "source", "channel", "status"],
    )
    write_csv(
        out_dir / "environment_minimal_imports_detected.csv",
        import_rows,
        ["file", "line", "module", "top_level_module", "third_party", "statement"],
    )
    (out_dir / "environment_minimal_missing_imports.txt").write_text(
        "\n".join(missing) + ("\n" if missing else ""),
        encoding="utf-8",
    )

    print("Minimal packages detected:")
    for rec in package_rows:
        print(
            f"  {rec['package_name']:<22} "
            f"{rec['version']:<12} "
            f"{rec['source']:<6} "
            f"from import `{rec['import_module']}`"
        )

    print("-" * 80)
    print("Saved:")
    print(f"  {env_yml_path}")
    print(f"  {out_dir / 'environment_minimal_used_requirements.txt'}")
    print(f"  {out_dir / 'environment_minimal_used_packages.csv'}")
    print(f"  {out_dir / 'environment_minimal_imports_detected.csv'}")
    print(f"  {out_dir / 'environment_minimal_missing_imports.txt'}")
    print("-" * 80)
    print("Recreate on another machine:")
    print("  conda env create -f environment_minimal_used.yml")
    print()
    print("This file contains only packages directly imported by the plotting code.")
    print("Conda will install required dependencies automatically.")

    if missing:
        print()
        print("WARNING: Some imports could not be resolved. Check:")
        print(f"  {out_dir / 'environment_minimal_missing_imports.txt'}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
