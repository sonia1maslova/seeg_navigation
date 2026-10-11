"""Plot saved 18-bin simulation inputs and method A/B comparison results.

Run after test_circular_vector_averaging_equivalence.py:
    python codes/python/pac/test/plot_circular_vector_averaging_equivalence.py
Requires numpy and matplotlib. This script does not simulate new data or
recalculate complex vectors. Use --overwrite explicitly to replace a PNG.
"""

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


def load_saved_data(result_dir):
    """读取已保存的归一化分布和复向量；不导入或调用计算脚本。"""
    cases = {}
    with (result_dir / "simulated_distributions.csv").open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            samples = cases.setdefault(row["case"], {})
            bins = samples.setdefault(int(row["sample"]), {})
            bin_number = int(row["phase_bin"])
            if bin_number in bins:
                raise ValueError("Duplicate sample/bin record")
            bins[bin_number] = (float(row["phase_center_deg"]), float(row["probability"]))
    if len(cases) != 6:
        raise ValueError("The demo figure requires the six saved simulation cases")
    profiles = {}
    for name, samples in cases.items():
        probability, centers = [], None
        for sample in sorted(samples):
            bins = samples[sample]
            if set(bins) != set(range(1, 19)):
                raise ValueError(f"Incomplete 18-bin profile: {name}, sample {sample}")
            theta = np.array([bins[k][0] for k in range(1, 19)])
            values = np.array([bins[k][1] for k in range(1, 19)])
            if not np.isfinite(theta).all() or not np.isfinite(values).all() or np.any(values < 0):
                raise ValueError("Nonfinite or negative saved profile values")
            if not np.isclose(values.sum(), 1, rtol=0, atol=1e-12):
                raise ValueError("Saved probabilities must sum to one")
            if centers is not None and not np.array_equal(theta, centers):
                raise ValueError("Samples have inconsistent phase-bin centers")
            centers = theta
            probability.append(values)
        profiles[name] = dict(theta_deg=centers, probability=np.array(probability))
    results = {}
    with (result_dir / "equivalence_results.csv").open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            name = row["case"]
            if name not in profiles:
                continue  # 随机压力测试只保存在结果表中，不逐组绘图。
            if name in results:
                raise ValueError(f"Duplicate comparison result: {name}")
            if int(row["n_samples"]) != len(profiles[name]["probability"]):
                raise ValueError(f"Sample count mismatch: {name}")
            values = {key: float(row[key]) for key in
                      ("a_real", "a_imag", "a_length", "b_real", "b_imag", "b_length")}
            if not all(np.isfinite(value) for value in values.values()):
                raise ValueError(f"Nonfinite saved vector: {name}")
            if row["phase_defined"] not in ("True", "False"):
                raise ValueError(f"Invalid phase-defined flag: {name}")
            values["phase_defined"] = row["phase_defined"] == "True"
            results[name] = values
    if set(results) != set(profiles):
        raise ValueError("Missing comparison results for a demo case")
    return profiles, results


def plot_demo(profiles, results, path):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "axes.labelweight": "bold", "axes.titleweight": "bold"})
    fig, axes = plt.subplots(3, 4, figsize=(11, 7.2))
    for i, (name, profile_data) in enumerate(profiles.items()):
        profile_ax, vector_ax = axes[i // 2, (i % 2) * 2:(i % 2) * 2 + 2]
        theta = profile_data["theta_deg"]
        probability = profile_data["probability"]
        result = results[name]
        for profile in probability:
            profile_ax.plot(theta, profile, color="#aaaaaa", lw=0.7, alpha=0.65)
        profile_ax.plot(theta, probability.mean(axis=0), color="black", lw=1.5)
        profile_ax.set(title=name, xlim=(-180, 180), ylim=(0, None),
                       xticks=[-180, 0, 180], xlabel="Phase-bin center (deg)", ylabel="Probability")
        for key, color, style, width in [("a", "#49a5cf", "-", 3),
                                          ("b", "black", "--", 1.2)]:
            x, y = result[key + "_real"], result[key + "_imag"]
            vector_ax.plot([0, x], [0, y], color=color, ls=style, lw=width)
            vector_ax.plot(x, y, "o", color=color, ms=3)
        bound = max(0.1, 1.2 * result["a_length"])
        vector_ax.axhline(0, color="#dddddd", lw=0.6)
        vector_ax.axvline(0, color="#dddddd", lw=0.6)
        vector_ax.set(xlim=(-bound, bound), ylim=(-bound, bound),
                      xlabel="Real part", ylabel="Imaginary part", aspect="equal")
        vector_ax.set_title(f"R = {result['a_length']:.5f}", fontsize=10)
        if not result["phase_defined"]:
            vector_ax.text(0.5, 0.12, "Phase undefined (R ~ 0)", ha="center",
                           transform=vector_ax.transAxes, fontsize=7)
        for ax in (profile_ax, vector_ax):
            ax.spines[["top", "right"]].set_visible(False)
            ax.tick_params(length=2)
    handles = [Line2D([0], [0], color="#aaaaaa", lw=1, label="Individual distributions"),
               Line2D([0], [0], color="black", lw=1.5, label="Mean distribution"),
               Line2D([0], [0], color="#49a5cf", lw=3, label="A: distribution first"),
               Line2D([0], [0], color="black", ls="--", label="B: vector first")]
    fig.legend(handles=handles, loc="upper center", ncol=4, frameon=False, fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.94), w_pad=1.5, h_pad=1.8)
    fig.savefig(path, dpi=200, facecolor="white")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[4])
    parser.add_argument("--results-dir", type=Path, help="Read an alternative saved results directory")
    parser.add_argument("--output", type=Path, help="Save the PNG at an alternative path")
    parser.add_argument("--overwrite", action="store_true", help="Allow replacing the output PNG")
    args = parser.parse_args()
    relative = Path("pac/test/circular_vector_averaging_equivalence")
    result_dir = args.results_dir if args.results_dir is not None else args.project_root / "results" / relative
    output = args.output if args.output is not None else args.project_root / "figures" / relative / "equivalence_demo.png"
    if output.suffix.lower() != ".png":
        parser.error("--output must be a PNG path")
    if output.exists() and not args.overwrite:
        raise FileExistsError(f"PNG already exists; use --overwrite explicitly: {output}")
    profiles, results = load_saved_data(result_dir)
    output.parent.mkdir(parents=True, exist_ok=True)
    plot_demo(profiles, results, output)
    print(f"Saved: {output}")


if __name__ == "__main__":
    main()
