#!/usr/bin/env python3
"""Plot C++ StochasticHeatEquation ASCII diagnostics.
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path

import numpy as np


KB = 1.380649e-23

REQUIRED_INPUTS = {
    "npts",
    "xlen",
    "rho",
    "cv",
    "crossA",
    "icor",
    "iper",
    "uinit",
    "jmidl",
    "jmidr",
    "midfact",
    "iresl",
    "iresr",
    "uleft",
    "uright",
}


def parse_value(text: str):
    token = text.strip().split()[0]
    if token.lower() in {"true", "false"}:
        return int(token.lower() == "true")
    try:
        if any(c in token for c in ".eE"):
            return float(token)
        return int(token)
    except ValueError:
        return token.strip("\"'")


def read_inputs(path: Path) -> dict:
    values = {}
    for line in path.read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if not line or "=" not in line:
            continue
        key, raw = line.split("=", 1)
        values[key.strip()] = parse_value(raw)
    return values


def require_inputs(params: dict, path: Path, keys: set[str]) -> None:
    missing = sorted(key for key in keys if key not in params)
    if missing:
        names = ", ".join(missing)
        raise ValueError(f"{path}: missing required input parameter(s): {names}")


def normalize_tag(tag: str) -> str:
    return tag.strip("_")


def label_sort_key(label: str) -> tuple[int, int]:
    ens, step = label.split("_", 1)
    return int(ens), int(step)


def find_labels(data_dir: Path, tag: str) -> list[str]:
    pattern = re.compile(rf"^mean_{re.escape(tag)}_(\d{{6}}_\d{{9}})\.dat$")
    labels = []
    for path in data_dir.glob(f"mean_{tag}_*.dat"):
        match = pattern.match(path.name)
        if match:
            labels.append(match.group(1))
    return sorted(labels, key=label_sort_key)


def read_xy_blocks(path: Path) -> list[np.ndarray]:
    blocks: list[list[tuple[float, float]]] = [[]]
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line:
            if blocks[-1]:
                blocks.append([])
            continue
        if line.startswith("#"):
            continue
        cols = line.split()
        if len(cols) < 2:
            raise ValueError(f"{path}: expected at least two columns in line {raw!r}")
        blocks[-1].append((float(cols[0]), float(cols[1])))

    return [np.asarray(block, dtype=float) for block in blocks if block]


def read_xy(path: Path) -> np.ndarray:
    blocks = read_xy_blocks(path)
    if not blocks:
        raise ValueError(f"{path}: no data found")
    return blocks[0]


def read_field(data_dir: Path, name: str, tag: str, label: str) -> np.ndarray:
    path = data_dir / f"{name}_{tag}_{label}.dat"
    if not path.exists():
        raise FileNotFoundError(path)
    return read_xy(path)


def reference_profile(
    npts: int, uinit: float, midfact: float, jmidl: int, jmidr: int
) -> np.ndarray:
    profile = np.full(npts, uinit, dtype=float)
    for i in range(npts):
        fortran_cell = i + 1
        if jmidl <= fortran_cell <= jmidr:
            profile[i] *= midfact
    return profile


def deterministic_profile(x: np.ndarray, params: dict) -> np.ndarray:
    if (
        int(params["iper"]) == 0
        and int(params["iresl"]) == 1
        and int(params["iresr"]) == 1
    ):
        return float(params["uleft"]) + (
            float(params["uright"]) - float(params["uleft"])
        ) * x / float(params["xlen"])

    return reference_profile(
        int(params["npts"]),
        float(params["uinit"]),
        float(params["midfact"]),
        int(params["jmidl"]),
        int(params["jmidr"]),
    )


def is_nonequilibrium(params: dict) -> bool:
    return (
        int(params["iper"]) == 0
        and int(params["iresl"]) == 1
        and int(params["iresr"]) == 1
        and not math.isclose(float(params["uleft"]), float(params["uright"]))
    )


def theory_variance(x: np.ndarray, steady: np.ndarray, params: dict, kb: float) -> np.ndarray:
    npts = int(params["npts"])
    dx = float(params["xlen"]) / npts
    dvol = float(params["crossA"]) * dx

    if is_nonequilibrium(params):
        profile = steady
    else:
        profile = np.full(npts, float(params["uinit"]))

    var = kb * profile**2 / (float(params["rho"]) * float(params["cv"]) * dvol)
    if int(params["iper"]) == 1:
        var *= 1.0 - 1.0 / npts
    return var


def theory_correlation(x: np.ndarray, steady: np.ndarray, params: dict, kb: float) -> np.ndarray:
    npts = int(params["npts"])
    icor = int(params["icor"]) - 1
    dx = float(params["xlen"]) / npts
    dvol = float(params["crossA"]) * dx

    if is_nonequilibrium(params):
        tref = steady[icor]
    else:
        tref = float(params["uinit"])

    var_ref = kb * tref**2 / (float(params["rho"]) * float(params["cv"]) * dvol)
    if int(params["iper"]) == 1:
        corr = -var_ref / npts * np.ones(npts)
        corr[icor] = var_ref * (1.0 - 1.0 / npts)
    else:
        corr = np.zeros(npts)
        corr[icor] = var_ref

    if is_nonequilibrium(params):
        length = float(params["xlen"])
        prefactor = (
            kb
            * (float(params["uright"]) - float(params["uleft"])) ** 2
            / (float(params["rho"]) * float(params["cv"]) * float(params["crossA"]) * length**3)
        )
        xstar = x[icor]
        for i, xi in enumerate(x):
            if i < icor:
                corr[i] += prefactor * xi * (length - xstar)
            else:
                corr[i] += prefactor * xstar * (length - xi)

    return corr


STATE_HEADER = re.compile(
    r"step,time\s*=\s*(\d+)\s+[-+0-9.eE]+.*ensemble\s*=\s*(\d+)"
)


def parse_state_from_log(path: Path, ensemble: int, step: int, npts: int) -> np.ndarray | None:
    lines = path.read_text(errors="replace").splitlines()
    for iline, line in enumerate(lines):
        match = STATE_HEADER.search(line)
        if not match:
            continue
        if int(match.group(1)) != step or int(match.group(2)) != ensemble:
            continue

        rows: list[tuple[float, float]] = []
        for raw in lines[iline + 1 :]:
            cols = raw.split()
            if len(cols) >= 2:
                try:
                    rows.append((float(cols[0]), float(cols[1])))
                except ValueError:
                    if rows:
                        break
                    continue
            elif rows:
                break
            if len(rows) == npts:
                return np.asarray(rows, dtype=float)
    return None


def find_state_log(data_dir: Path, requested: Path | None, ensemble: int, step: int, npts: int):
    candidates: list[Path]
    if requested is not None:
        candidates = [requested]
    else:
        candidates = sorted(data_dir.glob("*.out"), key=lambda p: p.stat().st_mtime, reverse=True)

    for path in candidates:
        state = parse_state_from_log(path, ensemble, step, npts)
        if state is not None:
            return state, path
    return None, None


def setup_matplotlib(show: bool):
    if not show:
        import matplotlib

        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def save_figure(fig, outdir: Path, stem: str, fmt: str) -> Path:
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / f"{stem}.{fmt}"
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    return path


def plot_initial(plt, x: np.ndarray, initial: np.ndarray, steady: np.ndarray):
    fig, ax = plt.subplots()
    ax.plot(x, initial, "g*", label="Initial state")
    ax.plot(x, steady, "b--", label="Deterministic")
    ax.set_xlabel("Position (m)")
    ax.set_ylabel("Temperature (K)")
    ax.set_title("Initial state (stars), steady state (dashed)")
    ax.legend()
    return fig


def plot_temperature(plt, x: np.ndarray, mean: np.ndarray, steady: np.ndarray, final_state):
    fig, ax = plt.subplots()
    if final_state is not None:
        ax.plot(final_state[:, 0], final_state[:, 1], "g*", label="Final state")
    ax.plot(x, mean, "ro", label="Average")
    ax.plot(x, steady, "b--", label="Deterministic")
    ax.set_xlabel("Position (m)")
    ax.set_ylabel("Temperature (K)")
    if final_state is not None:
        ax.set_title("Final state (stars), Average (circles), Deterministic (line)")
    else:
        ax.set_title("Average (circles), Deterministic (line)")
    ax.legend()
    return fig


def plot_variance(plt, x: np.ndarray, var: np.ndarray, var_th: np.ndarray):
    fig, ax = plt.subplots()
    ax.plot(x, var, "ro", label="Average")
    ax.plot(x, var_th, "b--", label="Equilibrium theory")
    ax.set_xlabel("Position (m)")
    ax.set_ylabel("Variance < dT**2 >")
    ax.set_title("Average (circles), Equilibrium theory (dashed)")
    ax.legend()
    return fig


def plot_correlation(plt, x: np.ndarray, corr: np.ndarray, corr_th: np.ndarray):
    fig, ax = plt.subplots()
    ax.plot(x, corr, "ro", label="Average")
    ax.plot(x, corr_th, "b--", label="Theory")
    ax.set_xlabel("Position (m)")
    ax.set_ylabel("Correlation < dT(x) dT(x*) >")
    ax.set_title("Average (circles), Theory (dashed)")
    ax.legend()
    return fig


def plot_correlation_without_delta(
    plt, x: np.ndarray, corr: np.ndarray, corr_th: np.ndarray, icor: int
):
    corr_masked = corr.copy()
    theory_masked = corr_th.copy()
    corr_masked[icor] = np.nan
    theory_masked[icor] = np.nan

    fig, ax = plt.subplots()
    ax.plot(x, corr_masked, "ro", label="Average")
    ax.plot(x, theory_masked, "b--", label="Theory")
    ax.set_xlabel("Position (m)")
    ax.set_ylabel("< dT(x) dT(x*) > [delta-function removed]")
    ax.set_title("Average (circles); Theory (dashed)")
    ax.legend()
    return fig


def structure_factor_from_correlation(corr: np.ndarray, icor: int) -> np.ndarray:
    corr_by_displacement = np.roll(corr, -icor)
    # The exact periodic equilibrium covariance is even in displacement. A
    # finite-sample row is not exactly even, so symmetrize before the DFT.
    corr_even = corr_by_displacement.copy()
    npts = len(corr_even)
    for r in range(1, (npts + 1) // 2):
        partner = npts - r
        avg = 0.5 * (corr_even[r] + corr_even[partner])
        corr_even[r] = avg
        corr_even[partner] = avg
    return len(corr_even) * np.real(np.fft.fft(corr_even))


def plot_structure_factor(plt, corr: np.ndarray, params: dict, kb: float, icor: int):
    npts = int(params["npts"])
    dx = float(params["xlen"]) / npts
    dvol = float(params["crossA"]) * dx
    var_ref = kb * float(params["uinit"]) ** 2 / (
        float(params["rho"]) * float(params["cv"]) * dvol
    )

    sk = structure_factor_from_correlation(corr, icor)
    nyq = npts // 2
    ki = np.arange(1, nyq)
    sk_eq = var_ref * npts * np.ones_like(ki, dtype=float)

    fig, ax = plt.subplots()
    ax.plot(ki, sk[1:nyq], "ro", label="Average")
    ax.plot(ki, sk_eq, "b--", label="Theory")
    ax.set_xlabel("Wavenumber index")
    ax.set_ylabel("Structure factor S(k)")
    ax.set_title("Average (circles); Theory (dashed)")
    ax.legend()
    return fig


def parse_args(argv: list[str]) -> argparse.Namespace:
    default_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(
        description="Generate notebook-style plots from StochasticHeatEquation ASCII output."
    )
    parser.add_argument("--data-dir", type=Path, default=default_dir)
    parser.add_argument(
        "--inputs",
        type=Path,
        default=None,
        help="AMReX inputs file. Defaults to <data-dir>/inputs.",
    )
    parser.add_argument(
        "--tag", default=None, help="ASCII tag, e.g. fv or ga. Defaults from inputs."
    )
    parser.add_argument("--label", default=None, help="Run label, e.g. 000002_022000000.")
    parser.add_argument(
        "--stdout-log",
        type=Path,
        default=None,
        help="Optional captured stdout file. If omitted, matching *.out files are searched.",
    )
    parser.add_argument("--outdir", type=Path, default=None)
    parser.add_argument("--format", default="png", help="Output image format.")
    parser.add_argument(
        "--kb", type=float, default=KB, help="Boltzmann constant for theory curves."
    )
    parser.add_argument("--show", action="store_true", help="Show figures interactively.")
    parser.add_argument(
        "--skip-structure",
        action="store_true",
        help="Do not derive the equilibrium structure-factor plot from cor11.",
    )
    return parser.parse_args(argv)


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    data_dir = args.data_dir.resolve()
    inputs_path = args.inputs.resolve() if args.inputs else data_dir / "inputs"
    params = read_inputs(inputs_path)
    required_inputs = set(REQUIRED_INPUTS)
    require_inputs(params, inputs_path, required_inputs)

    tag = normalize_tag(args.tag if args.tag is not None else "ga")
    labels = find_labels(data_dir, tag)
    if not labels:
        print(f"No mean_{tag}_*.dat files found in {data_dir}", file=sys.stderr)
        return 2

    label = args.label or labels[-1]
    if label not in labels:
        print(f"Label {label!r} was not found for tag {tag!r}. Available labels:", file=sys.stderr)
        for item in labels:
            print(f"  {item}", file=sys.stderr)
        return 2

    ensemble, step = label_sort_key(label)
    outdir = args.outdir.resolve() if args.outdir else data_dir / f"plots_{tag}_{label}"

    mean_xy = read_field(data_dir, "mean", tag, label)
    var_xy = read_field(data_dir, "var", tag, label)
    corr_xy = read_field(data_dir, "cor11", tag, label)

    x = mean_xy[:, 0]
    mean = mean_xy[:, 1]
    var = var_xy[:, 1]
    corr = corr_xy[:, 1]

    npts = int(params["npts"])
    icor = int(params["icor"]) - 1
    if len(x) != npts:
        print(f"Warning: inputs npts={npts}, but mean file has {len(x)} rows", file=sys.stderr)
        npts = len(x)
    if len(var) != npts or len(corr) != npts:
        raise ValueError("mean, var and first cor11 block must have matching row counts")

    initial = reference_profile(
        npts,
        float(params["uinit"]),
        float(params["midfact"]),
        int(params["jmidl"]),
        int(params["jmidr"]),
    )
    steady = deterministic_profile(x, params)
    var_th = theory_variance(x, steady, params, args.kb)
    corr_th = theory_correlation(x, steady, params, args.kb)

    final_state, log_path = find_state_log(data_dir, args.stdout_log, ensemble, step, npts)
    if args.stdout_log and final_state is None:
        print(
            f"Warning: {args.stdout_log} has no state for ensemble={ensemble}, step={step}",
            file=sys.stderr,
        )

    plt = setup_matplotlib(args.show)
    written: list[Path] = []
    written.append(
        save_figure(plot_initial(plt, x, initial, steady), outdir, "initial_state", args.format)
    )
    written.append(
        save_figure(
            plot_temperature(plt, x, mean, steady, final_state),
            outdir,
            "temperature_average",
            args.format,
        )
    )
    written.append(save_figure(plot_variance(plt, x, var, var_th), outdir, "variance", args.format))
    written.append(
        save_figure(plot_correlation(plt, x, corr, corr_th), outdir, "correlation", args.format)
    )
    written.append(
        save_figure(
            plot_correlation_without_delta(plt, x, corr, corr_th, icor),
            outdir,
            "correlation_delta_removed",
            args.format,
        )
    )

    if not args.skip_structure:
        if int(params["iper"]) == 1 and not is_nonequilibrium(params):
            fig = plot_structure_factor(plt, corr, params, args.kb, icor)
            written.append(save_figure(fig, outdir, "structure_factor", args.format))
        else:
            print("Structure factor not computed for non-periodic or non-equilibrium inputs")

    if args.show:
        plt.show()
    else:
        plt.close("all")

    print(f"Using inputs: {inputs_path}")
    print(f"Using diagnostics: tag={tag}, label={label}")
    if log_path is not None:
        print(f"Using final state from: {log_path}")
    else:
        print("No matching stdout state found; final-state markers were omitted")
    print(f"Correlation for x* = {x[icor]:.12e} (m)")
    for path in written:
        print(f"Wrote {path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
