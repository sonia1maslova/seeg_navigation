"""Compare averaging normalized 18-bin profiles with averaging full vectors.

Only simulated data are used. Run from the project with:
    python codes/python/pac/test/test_circular_vector_averaging_equivalence.py
Requires numpy only. Saves simulation inputs, comparison results and a JSON
summary. Existing output files are never overwritten. Plot the saved results
separately with plot_circular_vector_averaging_equivalence.py.
"""

import argparse
import csv
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform

import numpy as np


N_BINS = 18
TOLERANCE = 1e-14
ZERO_THRESHOLD = 1e-12
EDGES = np.linspace(-np.pi, np.pi, N_BINS + 1)
CENTERS = (EDGES[:-1] + EDGES[1:]) / 2
BASIS = np.exp(1j * CENTERS)
SOURCE = ("https://github.com/sonia1maslova/ai-assisted-learning-notes/blob/main/"
          "circular-statistics/circular_distribution_vector_averaging_equivalence_zh.md")


def normalize_rows(raw):
    """每个样本先独立归一化，不能先合并原始分布再归一化。"""
    raw = np.asarray(raw, dtype=np.float64)
    if raw.ndim != 2 or raw.shape[0] == 0 or raw.shape[1] != N_BINS:
        raise ValueError("Expected a nonempty samples x 18 array")
    if not np.isfinite(raw).all() or np.any(raw < 0):
        raise ValueError("All bin values must be finite and nonnegative")
    totals = raw.sum(axis=1, keepdims=True)
    if not np.isfinite(totals).all() or np.any(totals <= 0):
        raise ValueError("Every sample must have a finite positive total")
    return raw / totals


def method_a(raw):
    """方法 A：平均各样本归一化的分布，再计算合成复向量。"""
    probability = normalize_rows(raw)
    mean_probability = probability.mean(axis=0)
    return np.sum(mean_probability * BASIS)


def method_b(raw):
    """方法 B：各样本先计算复向量，再等权平均；保留各自长度。"""
    probability = normalize_rows(raw)
    sample_vectors = np.sum(probability * BASIS, axis=1)
    return sample_vectors.mean()


def exact_coefficients_match(raw):
    """整数输入用 Fraction 核验两表达式的 bin 系数严格相等。

    这不是用有理数精确计算 exp(i*theta)，而是核验它的线性系数。
    """
    n = len(raw)
    totals = [sum(int(value) for value in row) for row in raw]
    for k in range(N_BINS):
        a = sum((Fraction(int(row[k]), total)
                 for row, total in zip(raw, totals)), Fraction(0)) / n
        b = sum((Fraction(int(row[k]), total) / n
                 for row, total in zip(raw, totals)), Fraction(0))
        if a != b:
            return False
    return True


def demo_cases(rng):
    def peaked(means, kappas):
        shape = np.exp(np.asarray(kappas)[:, None] *
                       (np.cos(CENTERS - np.asarray(means)[:, None]) - 1))
        scales = np.array([60, 150, 400, 900, 2000, 5000])[:, None]
        return np.maximum(1, np.rint(scales * shape)).astype(np.int64)

    opposite = np.zeros((6, N_BINS), dtype=np.int64)
    opposite[:3, 3] = [10, 100, 1000]
    opposite[3:, 12] = [20, 200, 2000]
    return {
        "Random profiles": rng.integers(1, 1001, (6, N_BINS)),
        "Aligned peaks": peaked(np.deg2rad([20, 30, 40, 25, 35, 30]), [3] * 6),
        "Different concentrations": peaked(np.deg2rad([-60, -20, 30, 70, 110, 150]),
                                            [0.1, 0.4, 1, 2, 4, 8]),
        "Across +/-180 degrees": peaked(np.deg2rad([170, -170, 175, -175, 165, -165]),
                                        [4] * 6),
        "Uniform profiles": np.arange(1, 7)[:, None] * np.ones((6, N_BINS), dtype=int),
        "Opposite directions": opposite,
    }


def compare_case(name, raw, rng, exact=False):
    a, b = method_a(raw), method_b(raw)
    difference = float(abs(a - b))
    defined = min(abs(a), abs(b)) > ZERO_THRESHOLD
    # 用圆周角差，避免 +180/-180 的分界造成虚假的大差异。
    phase_difference = float(abs(np.rad2deg(np.angle(a * np.conj(b))))) if defined else None
    scales = 10.0 ** rng.uniform(-8, 8, (len(raw), 1))
    scale_error = float(max(abs(method_a(raw * scales) - a),
                            abs(method_b(raw * scales) - b)))
    exact_match = exact_coefficients_match(raw) if exact else None
    passed = difference <= TOLERANCE and scale_error <= TOLERANCE
    passed = passed and (exact_match is not False)
    return dict(case=name, n_samples=len(raw),
                a_real=float(a.real), a_imag=float(a.imag), a_length=float(abs(a)),
                a_phase_deg=float(np.rad2deg(np.angle(a))) if defined else None,
                b_real=float(b.real), b_imag=float(b.imag), b_length=float(abs(b)),
                b_phase_deg=float(np.rad2deg(np.angle(b))) if defined else None,
                complex_abs_difference=difference,
                length_abs_difference=float(abs(abs(a) - abs(b))),
                circular_phase_difference_deg=phase_difference,
                phase_defined=bool(defined), bitwise_equal=bool(a == b),
                sample_scaling_max_abs_error=scale_error,
                exact_bin_coefficients_equal=exact_match, passed=bool(passed))


def check_invalid_inputs():
    good = np.ones((2, N_BINS))
    invalid = [np.ones((2, 17)), np.ones(18), np.empty((0, 18)),
               np.vstack([good[0], np.zeros(18)]), -good,
               np.full((2, 18), np.nan), np.full((2, 18), np.inf)]
    rejected = 0
    for raw in invalid:
        for method in (method_a, method_b):
            try:
                method(raw)
            except ValueError:
                rejected += 1
            else:
                raise AssertionError("An invalid input was accepted")
    return rejected


def write_csv(path, rows):
    with path.open("x", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[4])
    parser.add_argument("--seed", type=int, default=20261008)
    parser.add_argument("--random-tests", type=int, default=1000)
    args = parser.parse_args()
    if args.random_tests < 1:
        parser.error("--random-tests must be positive")
    relative = Path("pac/test/circular_vector_averaging_equivalence")
    result_dir = args.project_root / "results" / relative
    if result_dir.exists() and any(result_dir.iterdir()):
        raise FileExistsError(f"Output directory is not empty: {result_dir}")
    rng = np.random.default_rng(args.seed)
    cases = demo_cases(rng)
    rows = [compare_case(name, raw, rng, exact=True) for name, raw in cases.items()]
    for i in range(args.random_tests):
        raw = rng.integers(0, 100001, (int(rng.integers(2, 65)), N_BINS))
        rows.append(compare_case(f"random_{i + 1:04d}", raw, rng))
    invalid_rejected = check_invalid_inputs()
    if not all(row["passed"] for row in rows):
        raise AssertionError("At least one equivalence or scaling test failed")
    inputs = []
    for name, raw in cases.items():
        for i, (values, probability) in enumerate(zip(raw, normalize_rows(raw)), 1):
            for k in range(N_BINS):
                inputs.append(dict(case=name, sample=i, phase_bin=k + 1,
                                   phase_center_deg=float(np.rad2deg(CENTERS[k])),
                                   raw_bin_value=int(values[k]), probability=float(probability[k])))
    summary = dict(seed=args.seed, n_bins=N_BINS, demo_cases=len(cases),
                   random_tests=args.random_tests, all_passed=True,
                   absolute_tolerance=TOLERANCE, zero_phase_threshold=ZERO_THRESHOLD,
                   max_complex_abs_difference=max(r["complex_abs_difference"] for r in rows),
                   max_length_abs_difference=max(r["length_abs_difference"] for r in rows),
                   max_sample_scaling_error=max(r["sample_scaling_max_abs_error"] for r in rows),
                   bitwise_equal_count=sum(r["bitwise_equal"] for r in rows),
                   exact_coefficient_pass_count=sum(r["exact_bin_coefficients_equal"] is True for r in rows),
                   invalid_method_calls_rejected=invalid_rejected, source=SOURCE,
                   python=platform.python_version(), numpy=np.__version__,
                   script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                   interpretation="Mathematically identical; floating-point order can change last bits. "
                                  "Exact Fraction checks compare bin coefficients, not complex exponentials. "
                                  "Phase is not interpreted when length is near zero.")
    result_dir.mkdir(parents=True, exist_ok=True)
    write_csv(result_dir / "simulated_distributions.csv", inputs)
    write_csv(result_dir / "equivalence_results.csv", rows)
    with (result_dir / "summary.json").open("x", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, ensure_ascii=False)
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
