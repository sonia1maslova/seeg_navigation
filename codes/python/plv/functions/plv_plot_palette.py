"""Shared categorical palette for PLV figures.

The values freeze the first five colors of Matplotlib's default ``tab10``
cycle, matching the historical/current Figure-6 PLV significance color order.
"""

PLV_BAND_COLORS = {
    (0.5, 4): "#1f77b4",
    (4, 12): "#ff7f0e",
    (12, 30): "#2ca02c",
    (30, 70): "#d62728",
    (70, 150): "#9467bd",
}

PLV_BAND_LABEL_COLORS = {
    "0.5-4 Hz": PLV_BAND_COLORS[(0.5, 4)],
    "4-12 Hz": PLV_BAND_COLORS[(4, 12)],
    "12-30 Hz": PLV_BAND_COLORS[(12, 30)],
    "30-70 Hz": PLV_BAND_COLORS[(30, 70)],
    "70-150 Hz": PLV_BAND_COLORS[(70, 150)],
}
