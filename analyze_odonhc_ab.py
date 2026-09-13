"""Compare experiment A (GNSS/INS) and B (GNSS/INS + ODO/NHC).

The script aligns both navigation solutions with dataset3/truth.nav, computes
NED position, velocity, and attitude errors, and writes figures plus a CSV
summary to dataset3/results/analysis_ab.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
RESULTS_DIR = ROOT / "dataset3" / "results"
OUTPUT_DIR = RESULTS_DIR / "analysis_ab"

A_NAV_PATH = RESULTS_DIR / "NavResult_A_GNSS_INS.nav"
B_NAV_PATH = RESULTS_DIR / "NavResult_B_GNSS_INS_ODONHC.nav"
TRUTH_PATH = ROOT / "dataset3" / "truth.nav"
A_STD_PATH = RESULTS_DIR / "NavSTD_A_GNSS_INS.txt"
B_STD_PATH = RESULTS_DIR / "NavSTD_B_GNSS_INS_ODONHC.txt"

WGS84_A = 6378137.0
WGS84_F = 1.0 / 298.257223563

COLORS = {"A": "#d95f02", "B": "#1b9e77", "truth": "#303030"}


def load_navigation(path: Path) -> np.ndarray:
    """Return [time, lat, lon, h, vn, ve, vd, roll, pitch, yaw]."""
    if not path.exists():
        raise FileNotFoundError(f"Missing navigation file: {path}")
    data = np.loadtxt(path, usecols=range(1, 11))
    if data.ndim != 2 or data.shape[1] != 10:
        raise ValueError(f"Unexpected navigation format: {path}")
    if not np.all(np.diff(data[:, 0]) > 0):
        raise ValueError(f"Time must increase strictly: {path}")
    return data


def interpolate_navigation(data: np.ndarray, target_time: np.ndarray) -> np.ndarray:
    result = np.empty((target_time.size, data.shape[1]), dtype=float)
    result[:, 0] = target_time
    for column in range(1, 9):
        result[:, column] = np.interp(target_time, data[:, 0], data[:, column])

    # Heading crosses 0/360 degrees, so unwrap it before interpolation.
    yaw_unwrapped = np.unwrap(np.deg2rad(data[:, 9]))
    result[:, 9] = np.rad2deg(np.interp(target_time, data[:, 0], yaw_unwrapped))
    return result


def wrap_degrees(angle: np.ndarray) -> np.ndarray:
    return (angle + 180.0) % 360.0 - 180.0


def radii_of_curvature(latitude_rad: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    eccentricity_sq = WGS84_F * (2.0 - WGS84_F)
    denominator = np.sqrt(1.0 - eccentricity_sq * np.sin(latitude_rad) ** 2)
    rn = WGS84_A / denominator
    rm = WGS84_A * (1.0 - eccentricity_sq) / denominator**3
    return rm, rn


def position_error_ned(solution: np.ndarray, truth: np.ndarray) -> np.ndarray:
    truth_lat = np.deg2rad(truth[:, 1])
    delta_lat = np.deg2rad(solution[:, 1] - truth[:, 1])
    delta_lon = np.deg2rad(solution[:, 2] - truth[:, 2])
    rm, rn = radii_of_curvature(truth_lat)
    north = (rm + truth[:, 3]) * delta_lat
    east = (rn + truth[:, 3]) * np.cos(truth_lat) * delta_lon
    down = -(solution[:, 3] - truth[:, 3])
    return np.column_stack((north, east, down))


def navigation_errors(solution: np.ndarray, truth: np.ndarray) -> dict[str, np.ndarray]:
    attitude_error = solution[:, 7:10] - truth[:, 7:10]
    attitude_error = wrap_degrees(attitude_error)
    return {
        "position": position_error_ned(solution, truth),
        "velocity": solution[:, 4:7] - truth[:, 4:7],
        "attitude": attitude_error,
    }


def axis_rmse(error: np.ndarray) -> np.ndarray:
    return np.sqrt(np.mean(error**2, axis=0))


def vector_rmse(error: np.ndarray, columns: slice | tuple[int, ...]) -> float:
    selected = error[:, columns]
    return float(np.sqrt(np.mean(np.sum(selected**2, axis=1))))


def metric_rows(errors_a: dict[str, np.ndarray], errors_b: dict[str, np.ndarray]):
    rows: list[tuple[str, str, float, float]] = []
    axis_names = {
        "position": ("Position North", "Position East", "Position Down"),
        "velocity": ("Velocity North", "Velocity East", "Velocity Down"),
        "attitude": ("Roll", "Pitch", "Yaw"),
    }
    units = {"position": "m", "velocity": "m/s", "attitude": "deg"}

    for group in ("position", "velocity", "attitude"):
        values_a = axis_rmse(errors_a[group])
        values_b = axis_rmse(errors_b[group])
        for label, value_a, value_b in zip(axis_names[group], values_a, values_b):
            rows.append((label, units[group], float(value_a), float(value_b)))

    rows.extend(
        [
            (
                "Position horizontal",
                "m",
                vector_rmse(errors_a["position"], (0, 1)),
                vector_rmse(errors_b["position"], (0, 1)),
            ),
            (
                "Position 3D",
                "m",
                vector_rmse(errors_a["position"], slice(None)),
                vector_rmse(errors_b["position"], slice(None)),
            ),
            (
                "Velocity 3D",
                "m/s",
                vector_rmse(errors_a["velocity"], slice(None)),
                vector_rmse(errors_b["velocity"], slice(None)),
            ),
            (
                "Attitude 3D",
                "deg",
                vector_rmse(errors_a["attitude"], slice(None)),
                vector_rmse(errors_b["attitude"], slice(None)),
            ),
        ]
    )
    return rows


def save_metrics(rows: list[tuple[str, str, float, float]]) -> None:
    csv_path = OUTPUT_DIR / "rmse_metrics.csv"
    with csv_path.open("w", newline="", encoding="utf-8-sig") as file:
        writer = csv.writer(file)
        writer.writerow(["Metric", "Unit", "A_RMSE", "B_RMSE", "Improvement_percent"])
        for name, unit, value_a, value_b in rows:
            improvement = 100.0 * (value_a - value_b) / value_a if value_a else np.nan
            writer.writerow([name, unit, f"{value_a:.9g}", f"{value_b:.9g}", f"{improvement:.4f}"])

    summary_path = OUTPUT_DIR / "summary.txt"
    with summary_path.open("w", encoding="utf-8") as file:
        file.write("A: GNSS/INS\nB: GNSS/INS + ODO/NHC\n")
        file.write("Positive improvement means B has lower RMSE than A.\n\n")
        file.write(f"{'Metric':<24} {'A RMSE':>12} {'B RMSE':>12} {'Improve':>10}\n")
        file.write("-" * 64 + "\n")
        for name, unit, value_a, value_b in rows:
            improvement = 100.0 * (value_a - value_b) / value_a if value_a else np.nan
            file.write(
                f"{name:<24} {value_a:>10.5f} {unit:<5} "
                f"{value_b:>10.5f} {unit:<5} {improvement:>8.2f}%\n"
            )


def plot_error_series(
    elapsed: np.ndarray,
    error_a: np.ndarray,
    error_b: np.ndarray,
    labels: tuple[str, str, str],
    unit: str,
    title: str,
    filename: str,
) -> None:
    fig, axes = plt.subplots(3, 1, figsize=(12, 8), sharex=True, constrained_layout=True)
    for index, axis in enumerate(axes):
        axis.plot(elapsed, error_a[:, index], color=COLORS["A"], linewidth=0.8, label="A: GNSS/INS")
        axis.plot(elapsed, error_b[:, index], color=COLORS["B"], linewidth=0.8, label="B: + ODO/NHC")
        axis.axhline(0.0, color="0.45", linewidth=0.6)
        axis.set_ylabel(f"{labels[index]} ({unit})")
        axis.grid(True, alpha=0.25)
    axes[0].legend(ncol=2)
    axes[-1].set_xlabel("Elapsed time (s)")
    fig.suptitle(title)
    fig.savefig(OUTPUT_DIR / filename, dpi=220)
    plt.close(fig)


def relative_ne(position: np.ndarray, origin: np.ndarray) -> np.ndarray:
    origin_lat = np.deg2rad(origin[0])
    rm, rn = radii_of_curvature(np.array([origin_lat]))
    north = (rm[0] + origin[2]) * np.deg2rad(position[:, 0] - origin[0])
    east = (rn[0] + origin[2]) * np.cos(origin_lat) * np.deg2rad(position[:, 1] - origin[1])
    return np.column_stack((north, east))


def plot_trajectory(solution_a: np.ndarray, solution_b: np.ndarray, truth: np.ndarray) -> None:
    origin = truth[0, 1:4]
    ne_a = relative_ne(solution_a[:, 1:4], origin)
    ne_b = relative_ne(solution_b[:, 1:4], origin)
    ne_truth = relative_ne(truth[:, 1:4], origin)

    step = max(1, len(truth) // 5000)
    fig, axis = plt.subplots(figsize=(9, 8), constrained_layout=True)
    axis.plot(ne_truth[::step, 1], ne_truth[::step, 0], color=COLORS["truth"], linewidth=2.0, label="Truth")
    axis.plot(ne_a[::step, 1], ne_a[::step, 0], color=COLORS["A"], linewidth=1.0, label="A: GNSS/INS")
    axis.plot(ne_b[::step, 1], ne_b[::step, 0], color=COLORS["B"], linewidth=1.0, label="B: + ODO/NHC")
    axis.set_xlabel("East (m)")
    axis.set_ylabel("North (m)")
    axis.set_title("Horizontal trajectory comparison")
    axis.axis("equal")
    axis.grid(True, alpha=0.25)
    axis.legend()
    fig.savefig(OUTPUT_DIR / "trajectory_comparison.png", dpi=220)
    plt.close(fig)


def plot_rmse(rows: list[tuple[str, str, float, float]]) -> None:
    row_map = {row[0]: row for row in rows}
    groups = [
        ("Position RMSE", ("Position North", "Position East", "Position Down", "Position horizontal"), "m"),
        ("Velocity RMSE", ("Velocity North", "Velocity East", "Velocity Down", "Velocity 3D"), "m/s"),
        ("Attitude RMSE", ("Roll", "Pitch", "Yaw", "Attitude 3D"), "deg"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True)
    for axis, (title, names, unit) in zip(axes, groups):
        x = np.arange(len(names))
        values_a = [row_map[name][2] for name in names]
        values_b = [row_map[name][3] for name in names]
        width = 0.38
        axis.bar(x - width / 2, values_a, width, color=COLORS["A"], label="A")
        axis.bar(x + width / 2, values_b, width, color=COLORS["B"], label="B")
        axis.set_xticks(x, [name.replace("Position ", "").replace("Velocity ", "") for name in names], rotation=25)
        axis.set_ylabel(unit)
        axis.set_title(title)
        axis.grid(axis="y", alpha=0.25)
    axes[0].legend(title="Experiment")
    fig.suptitle("RMSE comparison (lower is better)")
    fig.savefig(OUTPUT_DIR / "rmse_comparison.png", dpi=220)
    plt.close(fig)


def load_std(path: Path) -> np.ndarray:
    if not path.exists():
        raise FileNotFoundError(f"Missing STD file: {path}")
    return np.loadtxt(path)


def plot_std() -> None:
    std_a = load_std(A_STD_PATH)
    std_b = load_std(B_STD_PATH)
    start = max(std_a[0, 0], std_b[0, 0])
    end = min(std_a[-1, 0], std_b[-1, 0])
    time = std_a[(std_a[:, 0] >= start) & (std_a[:, 0] <= end), 0]
    elapsed = time - time[0]

    b_velocity = np.column_stack([np.interp(time, std_b[:, 0], std_b[:, i]) for i in range(4, 7)])
    b_attitude = np.column_stack([np.interp(time, std_b[:, 0], std_b[:, i]) for i in range(7, 10)])
    a_mask = (std_a[:, 0] >= start) & (std_a[:, 0] <= end)
    a_velocity = std_a[a_mask, 4:7]
    a_attitude = std_a[a_mask, 7:10]

    fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True, constrained_layout=True)
    axes[0].plot(elapsed, np.linalg.norm(a_velocity, axis=1), color=COLORS["A"], linewidth=0.9, label="A")
    axes[0].plot(elapsed, np.linalg.norm(b_velocity, axis=1), color=COLORS["B"], linewidth=0.9, label="B")
    axes[0].set_ylabel("Velocity STD norm (m/s)")
    axes[1].plot(elapsed, np.linalg.norm(a_attitude, axis=1), color=COLORS["A"], linewidth=0.9, label="A")
    axes[1].plot(elapsed, np.linalg.norm(b_attitude, axis=1), color=COLORS["B"], linewidth=0.9, label="B")
    axes[1].set_ylabel("Attitude STD norm (deg)")
    axes[1].set_xlabel("Elapsed time (s)")
    for axis in axes:
        axis.grid(True, alpha=0.25)
        axis.legend()
    fig.suptitle("Filter uncertainty comparison")
    fig.savefig(OUTPUT_DIR / "std_comparison.png", dpi=220)
    plt.close(fig)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    solution_a_raw = load_navigation(A_NAV_PATH)
    solution_b_raw = load_navigation(B_NAV_PATH)
    truth_raw = load_navigation(TRUTH_PATH)

    common_start = max(solution_a_raw[0, 0], solution_b_raw[0, 0], truth_raw[0, 0])
    common_end = min(solution_a_raw[-1, 0], solution_b_raw[-1, 0], truth_raw[-1, 0])
    mask_a = (solution_a_raw[:, 0] >= common_start) & (solution_a_raw[:, 0] <= common_end)
    target_time = solution_a_raw[mask_a, 0]

    solution_a = solution_a_raw[mask_a]
    solution_b = interpolate_navigation(solution_b_raw, target_time)
    truth = interpolate_navigation(truth_raw, target_time)

    errors_a = navigation_errors(solution_a, truth)
    errors_b = navigation_errors(solution_b, truth)
    rows = metric_rows(errors_a, errors_b)
    save_metrics(rows)

    elapsed = target_time - target_time[0]
    plot_trajectory(solution_a, solution_b, truth)
    plot_error_series(elapsed, errors_a["position"], errors_b["position"], ("North", "East", "Down"), "m", "Position errors", "position_errors.png")
    plot_error_series(elapsed, errors_a["velocity"], errors_b["velocity"], ("North", "East", "Down"), "m/s", "Velocity errors", "velocity_errors.png")
    plot_error_series(elapsed, errors_a["attitude"], errors_b["attitude"], ("Roll", "Pitch", "Yaw"), "deg", "Attitude errors", "attitude_errors.png")
    plot_rmse(rows)
    plot_std()

    print(f"Analysis complete: {OUTPUT_DIR}")
    for name, unit, value_a, value_b in rows:
        improvement = 100.0 * (value_a - value_b) / value_a if value_a else np.nan
        print(f"{name:24s} A={value_a:10.5f} {unit:5s} B={value_b:10.5f} {unit:5s} improvement={improvement:8.2f}%")


if __name__ == "__main__":
    main()
