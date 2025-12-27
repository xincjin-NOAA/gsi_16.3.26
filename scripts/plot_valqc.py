#!/usr/bin/env python3

import argparse
from dataclasses import dataclass
from typing import Iterable, Tuple

import numpy as np


@dataclass(frozen=True)
class ValqcParams:
    ratio_errors: float = 1.0
    cvar_pg: float = 0.0
    cvar_b: float = 8.0
    cg_term: float = 1.0
    tiny: float = 1.0e-12


def compute_valqc_from_val(val: np.ndarray, p: ValqcParams) -> np.ndarray:
    """Compute valqc using the same formulas as setupgnssrspd.f90.

    Here val == error * ddiff (i.e., normalized departure used in the Fortran).

    Fortran excerpt:
      exp_arg  = -half*val2
      rat_err2 = ratio_errors**2
      if (cvar_pg(ikx) > tiny_r_kind .and. error > tiny_r_kind) then
         arg  = exp(exp_arg)
         wnotgross= one-cvar_pg(ikx)
         cg_gnssrspd=cvar_b(ikx)
         wgross = cg_term*cvar_pg(ikx)/(cg_gnssrspd*wnotgross)
         term = log((arg+wgross)/(one+wgross))
      else
         term = exp_arg
      endif
      valqc = -two*rat_err2*term

    Notes:
    - This function uses p.tiny in place of tiny_r_kind.
    - We treat "error > tiny" as always true because val is already formed;
      if you want to mimic the branch exactly, pass cvar_pg=0 to disable varqc.
    """

    val = np.asarray(val, dtype=float)

    exp_arg = -0.5 * (val * val)
    rat_err2 = p.ratio_errors * p.ratio_errors

    if p.cvar_pg > p.tiny:
        # guard against pg->1
        wnotgross = 1.0 - p.cvar_pg
        if wnotgross <= p.tiny:
            raise ValueError(f"cvar_pg too close to 1 (cvar_pg={p.cvar_pg}); 1-cvar_pg must be > tiny")
        if p.cvar_b <= p.tiny:
            raise ValueError(f"cvar_b must be > tiny (cvar_b={p.cvar_b})")

        arg = np.exp(exp_arg)
        wgross = p.cg_term * p.cvar_pg / (p.cvar_b * wnotgross)
        term = np.log((arg + wgross) / (1.0 + wgross))
    else:
        term = exp_arg

    valqc = -2.0 * rat_err2 * term
    return valqc


def _try_import_matplotlib():
    try:
        import matplotlib.pyplot as plt  # noqa: WPS433

        return plt
    except Exception as e:  # pragma: no cover
        raise RuntimeError(
            "matplotlib is required for plotting. Install with: python -m pip install matplotlib numpy"
        ) from e


def _parse_csv_floats(s: str) -> Tuple[float, ...]:
    parts = [p.strip() for p in s.split(",") if p.strip()]
    if not parts:
        raise ValueError("Expected a comma-separated list of floats")
    return tuple(float(x) for x in parts)


def plot_valqc(
    val: np.ndarray,
    sweep: Iterable[Tuple[str, Iterable[ValqcParams]]],
    out: str,
    title: str,
    ylim: Tuple[float, float] | None,
):
    plt = _try_import_matplotlib()

    sweep = list(sweep)
    n = len(sweep)
    fig, axes = plt.subplots(1, n, figsize=(6.5 * n, 5.0), sharey=True)
    if n == 1:
        axes = [axes]

    for ax, (panel_title, param_list) in zip(axes, sweep, strict=True):
        for p in param_list:
            y = compute_valqc_from_val(val, p)
            label = f"pg={p.cvar_pg:g}, b={p.cvar_b:g}, ratio={p.ratio_errors:g}, cg_term={p.cg_term:g}"
            ax.plot(val, y, linewidth=2.0, label=label)

        ax.set_title(panel_title)
        ax.set_xlabel("val = error * ddiff")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)

    axes[0].set_ylabel("valqc")
    if ylim is not None:
        axes[0].set_ylim(*ylim)

    fig.suptitle(title)
    fig.tight_layout()

    fig.savefig(out, dpi=150)
    print(f"Wrote: {out}")


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Plot valqc vs normalized departure val (= error * ddiff) using the same formulas as GSI setupgnssrspd.f90.\n\n"
            "Dependencies: numpy, matplotlib"
        )
    )
    ap.add_argument("--out", default="valqc_sweeps.png", help="Output image file (png/pdf/svg).")
    ap.add_argument("--xmax", type=float, default=8.0, help="Plot range for val: [-xmax, xmax].")
    ap.add_argument("--npoints", type=int, default=1001, help="Number of points along x.")

    ap.add_argument("--base_ratio", type=float, default=1.0, help="Baseline ratio_errors.")
    ap.add_argument("--base_pg", type=float, default=0.0, help="Baseline cvar_pg.")
    ap.add_argument("--base_b", type=float, default=8.0, help="Baseline cvar_b.")
    ap.add_argument("--cg_term", type=float, default=1.0, help="Gross-error constant cg_term.")

    ap.add_argument(
        "--pg_list",
        default="0,0.01,0.05,0.1",
        help="Comma-separated sweep values for cvar_pg (panel 1).",
    )
    ap.add_argument(
        "--b_list",
        default="2,4,8,16",
        help="Comma-separated sweep values for cvar_b (panel 2).",
    )
    ap.add_argument(
        "--ratio_list",
        default="0.5,1,1.5,2",
        help="Comma-separated sweep values for ratio_errors (panel 3).",
    )

    ap.add_argument(
        "--ylim",
        default="",
        help="Optional y-limits as 'ymin,ymax' (e.g. 0,30). Leave empty for autoscale.",
    )

    args = ap.parse_args()

    val = np.linspace(-abs(args.xmax), abs(args.xmax), int(args.npoints))

    base = ValqcParams(
        ratio_errors=float(args.base_ratio),
        cvar_pg=float(args.base_pg),
        cvar_b=float(args.base_b),
        cg_term=float(args.cg_term),
    )

    pg_list = _parse_csv_floats(args.pg_list)
    b_list = _parse_csv_floats(args.b_list)
    ratio_list = _parse_csv_floats(args.ratio_list)

    sweeps = [
        (
            "Sweep cvar_pg (prob. gross)",
            [ValqcParams(base.ratio_errors, pg, base.cvar_b, base.cg_term) for pg in pg_list],
        ),
        (
            "Sweep cvar_b (gross scale)",
            [ValqcParams(base.ratio_errors, base.cvar_pg, b, base.cg_term) for b in b_list],
        ),
        (
            "Sweep ratio_errors",
            [ValqcParams(r, base.cvar_pg, base.cvar_b, base.cg_term) for r in ratio_list],
        ),
    ]

    ylim = None
    if args.ylim.strip():
        yparts = _parse_csv_floats(args.ylim)
        if len(yparts) != 2:
            raise ValueError("--ylim must be 'ymin,ymax'")
        ylim = (float(yparts[0]), float(yparts[1]))

    title = (
        "valqc vs val (GSI varQC formula)\n"
        f"baseline: pg={base.cvar_pg:g}, b={base.cvar_b:g}, ratio={base.ratio_errors:g}, cg_term={base.cg_term:g}"
    )

    plot_valqc(val=val, sweep=sweeps, out=args.out, title=title, ylim=ylim)


if __name__ == "__main__":
    main()
