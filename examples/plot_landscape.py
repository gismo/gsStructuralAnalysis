#!/usr/bin/env python3
"""Plot the equilibrium diagram of a ``gsALMLandscape`` CSV export.

The ``example_BratuExploration`` and ``example_ModifiedBratuExploration``
drivers of ``gsStructuralAnalysis`` write ``<ResultsDir>/landscape.csv`` via
``gsALMLandscape::writeCsv``.  This script renders the classical
:math:`\\lambda`--:math:`\\|u\\|` equilibrium diagram from such a file:

* one poly-line per ``curve``, ordered by ``point``;
* solid where the branch is stable (``stability == +1``), dashed where it is
  unstable (``stability == -1``), dotted where stability was never set
  (``stability == 0``) -- a single curve is split at every sign change;
* ``negatives`` (tangent inertia / negative-pivot count) is parsed and
  reported per curve; ``-1`` means "not recorded";
* bifurcation points (``isBifurcation == 1``) are marked;
* the :math:`\\lambda`-maximum (fold / limit point) of every curve is marked
  and annotated;
* child curves are joined to the parent point they branched from,
  ``(parentCurve, parentPointIdx)``.

Notes
-----
``gsALMLandscape::writeCsv`` sets no stream precision, so the floating point
columns of the input carry only about six significant figures.  Nothing
printed or drawn by this script is more accurate than that, and the summary
deliberately reports six significant digits only.

Only the standard library, NumPy and Matplotlib are used.  The ``Agg``
backend is selected so the script runs headless.

Examples
--------
Render the Bratu landscape and compare its fold against the analytic value::

    python3 plot_landscape.py BratuExplorationResults/landscape.csv \\
        --reference 6.808124423 -o bratu.png

Run the built-in self test::

    python3 plot_landscape.py --selftest

Compare the merged loci against the thesis figure (needs ``h5py`` and the
``landscape.h5`` the driver writes next to the CSV)::

    python3 plot_landscape.py ModifiedBratuResults/landscape.csv \\
        --merge-loci --window 0.01 --thesis-metrics --norm l2

Locus-level thesis-comparison mode
----------------------------------
``--merge-loci`` reclassifies every point by its solution *profile*, read
from the landscape HDF5 export (same basename as the CSV, ``.h5``; datasets
``curve_%04d_U`` of shape ``(N_dof, n_points)``), and merges all traced
segments of one solution *locus* into a single polyline.  ``--window``
clips to the thesis's continuation window first, solver-jump links are
removed from a measured step-length criterion, and ``--thesis-metrics``
prints the per-locus deviation from the analytic constant branch and from
polylines digitized out of the thesis figure render (``--thesis-fig``).
This mode needs ``h5py``; without it (or without the ``.h5`` file) the
script refuses ``--merge-loci`` with a clean message, because a norm-only
classification cannot separate the loci.

One source of the segmentation this flag was written to repair is gone:
the explorer now sweeps both arc-length directions of a seed into **one**
curve, so a locus is no longer split merely because it was traced
forwards and backwards.  Merging is still needed for the rest -- see
:func:`merge_loci`.

State norm
----------
The comparison quantity is a norm of the solution *function* psi, and
``--norm`` selects its convention.  The default ``l2`` is the continuous
norm ``sqrt(int_0^1 psi^2)``, evaluated exactly from the coefficients and
the mass matrix of the (reconstructed) spline basis.  ``discrete``
reproduces the project's historical ``||U||_2/sqrt(N_dof)``, which weights
the Greville abscissae uniformly although the two end spacings are half
the interior spacing: it is exact on constant states and biased on every
non-constant one, by up to 3 % on the traced mode-2 states.  ``h1`` adds
the derivative term, ``thesis`` mimics the thesis's own discrete
convention (see :func:`state_norm`).  The basis is reconstructed from
``N_dof`` because the HDF5 export stores no basis description; override
the assumed degree with ``--basis-degree``, which is checked against the
knot structure the export can have come from
(:func:`check_basis_degree`).
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import io
import os
import sys
import tempfile
from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402  (must follow matplotlib.use)
import numpy as np  # noqa: E402

try:  # optional: only the --merge-loci mode needs it
    import h5py  # noqa: E402

    _HAS_H5PY = True
except ImportError:  # pragma: no cover - depends on the environment
    _HAS_H5PY = False

__all__ = [
    "LandscapeError",
    "Curve",
    "read_landscape",
    "summarize",
    "plot_landscape",
    "read_solution_vectors",
    "open_uniform_knots",
    "bspline_basis",
    "spline_gram_matrices",
    "check_basis_degree",
    "state_norm",
    "norm_label",
    "count_sign_changes",
    "classify_locus",
    "derive_jump_cut",
    "merge_loci",
    "summarize_loci",
    "plot_loci",
    "locus_a_deviation",
    "digitize_thesis_fig",
    "thesis_metrics",
]

#: Columns that must be present in the CSV header.
REQUIRED_COLUMNS: Tuple[str, ...] = ("curve", "point", "L", "normU")

#: Columns that may be absent; the value used when they are.
OPTIONAL_COLUMNS: Dict[str, int] = {
    "stability": 0,
    "negatives": -1,
    "isBifurcation": 0,
    "parentCurve": -1,
    "parentPointIdx": -1,
    "equilibrium": 1,
}

#: Line style per stability flag.
_STABILITY_STYLE: Dict[int, str] = {1: "-", -1: "--", 0: ":"}

#: Human readable name per stability flag.
_STABILITY_NAME: Dict[int, str] = {1: "stable", -1: "unstable", 0: "unknown"}

#: Number of significant digits the input actually carries.
_SIGFIG = 6


class LandscapeError(Exception):
    """Raised when the input file cannot be interpreted as a landscape CSV.

    The message is intended to be shown to the user verbatim, without a
    traceback.
    """


class Curve:
    """A single continuation branch of the landscape.

    Parameters
    ----------
    cid : int
        Curve identifier, the value of the ``curve`` column.
    lam : numpy.ndarray
        Load / continuation parameter :math:`\\lambda` of each point, ordered
        by the ``point`` column.  Shape ``(n,)``.
    norm_u : numpy.ndarray
        Solution norm :math:`\\|u\\|` of each point.  Shape ``(n,)``.
    stability : numpy.ndarray
        Stability flag of each point, in ``{-1, 0, +1}``.  Shape ``(n,)``.
    bifurcation : numpy.ndarray
        Boolean mask of the bifurcation points.  Shape ``(n,)``.
    parent_curve : int
        Identifier of the curve this branch emanates from, or ``-1`` for a
        root curve.
    parent_point : int
        Index (``point`` value) on ``parent_curve`` this branch emanates
        from, or ``-1``.
    negatives : numpy.ndarray, optional
        Tangent inertia (negative-eigenvalue/pivot count) of each point.
        Shape ``(n,)``.  ``-1`` means "not recorded"; defaults to an
        all-``-1`` array when omitted.

    Attributes
    ----------
    size : int
        Number of points on the curve.
    """

    def __init__(
        self,
        cid: int,
        lam: np.ndarray,
        norm_u: np.ndarray,
        stability: np.ndarray,
        bifurcation: np.ndarray,
        parent_curve: int,
        parent_point: int,
        negatives: Optional[np.ndarray] = None,
    ) -> None:
        self.cid = cid
        self.lam = lam
        self.norm_u = norm_u
        self.stability = stability
        self.bifurcation = bifurcation
        self.parent_curve = parent_curve
        self.parent_point = parent_point
        self.negatives = (
            negatives if negatives is not None else np.full(lam.size, -1, dtype=int)
        )
        # Values of the ``point`` column; overwritten by read_landscape with
        # the values actually found in the file.
        self._point_ids = np.arange(lam.size, dtype=int)

    @property
    def size(self) -> int:
        """int: Number of points on this curve."""
        return int(self.lam.size)

    def lambda_max_index(self) -> int:
        """Index of the point with the largest :math:`\\lambda`.

        Returns
        -------
        int
            Position of the :math:`\\lambda`-maximum within the curve arrays.
            The first occurrence is returned when the maximum is attained
            several times.
        """
        return int(np.argmax(self.lam))

    def fold_is_interior(self) -> bool:
        """Whether the :math:`\\lambda`-maximum is a genuine interior fold.

        Returns
        -------
        bool
            ``True`` when the :math:`\\lambda`-maximum lies strictly inside
            the traced arc, i.e. the branch actually turns.  ``False`` when it
            sits on an end point, in which case the maximum is only where the
            continuation happened to stop.
        """
        idx = self.lambda_max_index()
        return 0 < idx < self.size - 1

    def stability_transitions(self) -> List[Tuple[int, int, int]]:
        """Locate the changes of the stability flag along the curve.

        Returns
        -------
        list of tuple of int
            One ``(index, before, after)`` triple per transition, where
            ``index`` is the position of the first point carrying the new
            flag.
        """
        out: List[Tuple[int, int, int]] = []
        for i in range(1, self.size):
            if self.stability[i] != self.stability[i - 1]:
                out.append((i, int(self.stability[i - 1]), int(self.stability[i])))
        return out

    def point_at(self, point_index: int) -> Optional[Tuple[float, float]]:
        """Coordinates of the point whose ``point`` column equals a value.

        Parameters
        ----------
        point_index : int
            Value of the ``point`` column being looked up.

        Returns
        -------
        tuple of float or None
            ``(lambda, normU)`` of that point, or ``None`` when the curve has
            no such point.
        """
        where = np.flatnonzero(self._point_ids == point_index)
        if where.size == 0:
            return None
        j = int(where[0])
        return float(self.lam[j]), float(self.norm_u[j])


def _to_float(text: str, column: str, line: int) -> float:
    """Convert a CSV field to ``float`` or raise :class:`LandscapeError`.

    Parameters
    ----------
    text : str
        Raw field content.
    column : str
        Column name, used in the error message.
    line : int
        1-based line number in the file, used in the error message.

    Returns
    -------
    float
        The parsed value.

    Raises
    ------
    LandscapeError
        If the field is empty or not a floating point literal.
    """
    try:
        return float(text)
    except (TypeError, ValueError):
        raise LandscapeError(
            "line {}: column '{}' is not a number (got {!r})".format(
                line, column, text
            )
        )


def _to_int(text: str, column: str, line: int) -> int:
    """Convert a CSV field to ``int`` or raise :class:`LandscapeError`.

    Values written as floats (``"3"`` or ``"3.0"``) are both accepted.

    Parameters
    ----------
    text : str
        Raw field content.
    column : str
        Column name, used in the error message.
    line : int
        1-based line number in the file, used in the error message.

    Returns
    -------
    int
        The parsed value.

    Raises
    ------
    LandscapeError
        If the field cannot be read as an integer.
    """
    value = _to_float(text, column, line)
    if value != int(value):
        raise LandscapeError(
            "line {}: column '{}' must be an integer (got {!r})".format(
                line, column, text
            )
        )
    return int(value)


def read_landscape(csv_path: str) -> List[Curve]:
    """Read a ``gsALMLandscape`` CSV export.

    Parsing is keyed on the header names, never on column position, so an
    extra column added to ``writeCsv`` in the future does not shift the
    interpretation of the existing ones.

    Parameters
    ----------
    csv_path : str
        Path to ``landscape.csv``.

    Returns
    -------
    list of Curve
        One entry per distinct value of the ``curve`` column, sorted by that
        value.  Points within a curve are sorted by the ``point`` column.

    Raises
    ------
    LandscapeError
        If the file is missing, empty, has no header, lacks one of
        ``curve``, ``point``, ``L``, ``normU``, or contains a row that cannot
        be parsed.
    """
    try:
        with open(csv_path, "r", newline="") as handle:
            reader = csv.DictReader(handle)
            fieldnames = reader.fieldnames
            if not fieldnames:
                raise LandscapeError(
                    "'{}' is empty or has no header row".format(csv_path)
                )
            header = [name.strip() for name in fieldnames]
            missing = [c for c in REQUIRED_COLUMNS if c not in header]
            if missing:
                raise LandscapeError(
                    "'{}' is missing required column(s): {}. Found: {}".format(
                        csv_path, ", ".join(missing), ", ".join(header)
                    )
                )
            reader.fieldnames = header

            records: List[Dict[str, float]] = []
            for line, row in enumerate(reader, start=2):
                if row.get(None) is not None:
                    raise LandscapeError(
                        "line {}: more fields than the header declares".format(line)
                    )
                # Skip fully blank trailing lines.
                if all((v is None or str(v).strip() == "") for v in row.values()):
                    continue
                rec: Dict[str, float] = {
                    "curve": _to_int(row["curve"], "curve", line),
                    "point": _to_int(row["point"], "point", line),
                    "L": _to_float(row["L"], "L", line),
                    "normU": _to_float(row["normU"], "normU", line),
                }
                for name, default in OPTIONAL_COLUMNS.items():
                    raw = row.get(name)
                    if raw is None or str(raw).strip() == "":
                        rec[name] = default
                    else:
                        rec[name] = _to_int(raw, name, line)
                records.append(rec)
    except LandscapeError:
        raise
    except FileNotFoundError:
        raise LandscapeError("no such file: '{}'".format(csv_path))
    except IsADirectoryError:
        raise LandscapeError("'{}' is a directory, not a CSV file".format(csv_path))
    except OSError as exc:
        raise LandscapeError("cannot read '{}': {}".format(csv_path, exc))
    except UnicodeDecodeError:
        raise LandscapeError(
            "'{}' is not a text file (binary content)".format(csv_path)
        )

    if not records:
        raise LandscapeError("'{}' contains a header but no data rows".format(csv_path))

    by_id: Dict[int, List[Dict[str, float]]] = {}
    for rec in records:
        by_id.setdefault(int(rec["curve"]), []).append(rec)

    curves: List[Curve] = []
    for cid in sorted(by_id):
        rows = sorted(by_id[cid], key=lambda r: r["point"])
        curve = Curve(
            cid=cid,
            lam=np.array([r["L"] for r in rows], dtype=float),
            norm_u=np.array([r["normU"] for r in rows], dtype=float),
            stability=np.array([int(r["stability"]) for r in rows], dtype=int),
            bifurcation=np.array(
                [int(r["isBifurcation"]) == 1 for r in rows], dtype=bool
            ),
            parent_curve=int(rows[0]["parentCurve"]),
            parent_point=int(rows[0]["parentPointIdx"]),
            negatives=np.array([int(r["negatives"]) for r in rows], dtype=int),
        )
        curve._point_ids = np.array([int(r["point"]) for r in rows], dtype=int)
        curves.append(curve)
    return curves


def _stability_runs(stability: np.ndarray) -> List[Tuple[int, int, int]]:
    """Split an index range into maximal runs of constant stability.

    Parameters
    ----------
    stability : numpy.ndarray
        Per-point stability flags.

    Returns
    -------
    list of tuple of int
        ``(start, stop, flag)`` triples with ``stop`` exclusive, covering the
        whole array.
    """
    n = int(stability.size)
    runs: List[Tuple[int, int, int]] = []
    if n == 0:
        return runs
    start = 0
    for i in range(1, n):
        if stability[i] != stability[start]:
            runs.append((start, i, int(stability[start])))
            start = i
    runs.append((start, n, int(stability[start])))
    return runs


def summarize(
    curves: Sequence[Curve], reference: Optional[float] = None
) -> str:
    """Build the textual report of a landscape.

    Parameters
    ----------
    curves : sequence of Curve
        Curves as returned by :func:`read_landscape`.
    reference : float, optional
        Reference value of :math:`\\lambda` (e.g. the analytic Bratu fold).
        When given, the signed deviation of the largest traced
        :math:`\\lambda` from it is reported.

    Returns
    -------
    str
        Multi-line human readable summary.
    """
    lines: List[str] = []
    lines.append("Landscape summary ({} curve(s))".format(len(curves)))
    lines.append(
        "  note: writeCsv writes ~{} significant figures; "
        "values below are not more precise.".format(_SIGFIG)
    )
    lines.append("")

    total_bif = 0
    global_max = -np.inf
    global_max_curve = -1
    for curve in curves:
        idx = curve.lambda_max_index()
        lam_max = float(curve.lam[idx])
        if lam_max > global_max:
            global_max = lam_max
            global_max_curve = curve.cid
        n_bif = int(np.count_nonzero(curve.bifurcation))
        total_bif += n_bif

        parent = (
            "root"
            if curve.parent_curve < 0
            else "child of curve {} at point {}".format(
                curve.parent_curve, curve.parent_point
            )
        )
        lines.append("curve {} ({})".format(curve.cid, parent))
        lines.append("  points            : {}".format(curve.size))
        lines.append(
            "  lambda range      : [{:.6g}, {:.6g}]".format(
                float(np.min(curve.lam)), float(np.max(curve.lam))
            )
        )
        lines.append(
            "  lambda_max        : {:.6g} at point {} ({})".format(
                lam_max,
                int(curve._point_ids[idx]),
                "interior fold" if curve.fold_is_interior() else "end of branch",
            )
        )
        lines.append("  |u| at lambda_max : {:.6g}".format(float(curve.norm_u[idx])))

        flags = sorted(set(int(s) for s in curve.stability))
        lines.append(
            "  stability         : {}".format(
                ", ".join(
                    "{} x {}".format(
                        int(np.count_nonzero(curve.stability == f)),
                        _STABILITY_NAME.get(f, "flag {}".format(f)),
                    )
                    for f in flags
                )
            )
        )
        neg_flags = sorted(set(int(v) for v in curve.negatives))
        lines.append(
            "  negatives         : {}".format(
                ", ".join(
                    "{} x {}".format(
                        int(np.count_nonzero(curve.negatives == f)),
                        "not recorded" if f == -1 else str(f),
                    )
                    for f in neg_flags
                )
            )
        )
        transitions = curve.stability_transitions()
        if transitions:
            lines.append(
                "  transitions       : {} ({})".format(
                    len(transitions),
                    "; ".join(
                        "point {}: {} -> {}".format(
                            int(curve._point_ids[i]),
                            _STABILITY_NAME.get(a, a),
                            _STABILITY_NAME.get(b, b),
                        )
                        for i, a, b in transitions
                    ),
                )
            )
        else:
            lines.append("  transitions       : 0")
        lines.append("  bifurcations      : {}".format(n_bif))
        if n_bif:
            lines.append(
                "    at lambda = {}".format(
                    ", ".join(
                        "{:.6g}".format(v) for v in curve.lam[curve.bifurcation]
                    )
                )
            )
        lines.append("")

    lines.append(
        "total: {} point(s), {} bifurcation(s)".format(
            sum(c.size for c in curves), total_bif
        )
    )
    lines.append(
        "global lambda_max: {:.6g} (curve {})".format(global_max, global_max_curve)
    )
    if reference is not None:
        deviation = global_max - reference
        rel = deviation / reference if reference != 0.0 else float("nan")
        lines.append("reference lambda : {:.10g}".format(reference))
        lines.append(
            "deviation        : {:+.6g} (relative {:+.3g})".format(deviation, rel)
        )
    return "\n".join(lines)


def plot_landscape(
    curves: Sequence[Curve],
    output: str,
    title: Optional[str] = None,
    reference: Optional[float] = None,
    dpi: int = 150,
) -> str:
    """Render the :math:`\\lambda`--:math:`\\|u\\|` equilibrium diagram.

    Parameters
    ----------
    curves : sequence of Curve
        Curves as returned by :func:`read_landscape`.
    output : str
        Path of the image file to write.  The extension selects the format.
    title : str, optional
        Figure title.
    reference : float, optional
        If given, a horizontal reference line is drawn at this
        :math:`\\lambda`.
    dpi : int, optional
        Output resolution, by default 150.

    Returns
    -------
    str
        The ``output`` path, for convenience.

    Raises
    ------
    LandscapeError
        If the figure cannot be written.
    """
    by_id: Dict[int, Curve] = {c.cid: c for c in curves}
    cmap = plt.get_cmap("tab10")

    fig, ax = plt.subplots(figsize=(8.0, 6.0))

    # Branch connectors first, so they sit underneath the curves.
    connector_drawn = False
    for curve in curves:
        if curve.parent_curve < 0 or curve.parent_curve not in by_id:
            continue
        anchor = by_id[curve.parent_curve].point_at(curve.parent_point)
        if anchor is None or curve.size == 0:
            continue
        ax.plot(
            [anchor[0], float(curve.lam[0])],
            [anchor[1], float(curve.norm_u[0])],
            color="0.55",
            linewidth=0.9,
            linestyle=(0, (2, 2)),
            zorder=1,
            label="branch connector" if not connector_drawn else None,
        )
        connector_drawn = True

    # The curve carrying the overall largest lambda; its maximum is always
    # annotated even when it only marks where the continuation stopped.
    global_max_cid = max(
        (c for c in curves if c.size), key=lambda c: float(np.max(c.lam)), default=None
    )
    global_max_cid = global_max_cid.cid if global_max_cid is not None else None

    for k, curve in enumerate(curves):
        color = cmap(k % 10)
        if curve.size == 1:
            ax.plot(
                curve.lam,
                curve.norm_u,
                marker="o",
                markersize=4,
                linestyle="none",
                color=color,
                zorder=3,
                label="curve {}".format(curve.cid),
            )
        else:
            labelled = False
            for start, stop, flag in _stability_runs(curve.stability):
                # Extend one point to the left so consecutive runs touch and
                # the poly-line shows no gap at a stability change.
                lo = max(start - 1, 0)
                ax.plot(
                    curve.lam[lo:stop],
                    curve.norm_u[lo:stop],
                    color=color,
                    linestyle=_STABILITY_STYLE.get(flag, ":"),
                    linewidth=1.6,
                    zorder=3,
                    label=None if labelled else "curve {}".format(curve.cid),
                )
                labelled = True

        # lambda-maximum (fold / limit point).
        idx = curve.lambda_max_index()
        interior = curve.fold_is_interior()
        ax.plot(
            [curve.lam[idx]],
            [curve.norm_u[idx]],
            marker="v" if interior else "|",
            markersize=9 if interior else 10,
            markerfacecolor="none",
            markeredgecolor=color,
            markeredgewidth=1.6,
            linestyle="none",
            zorder=5,
        )
        # Annotating every curve maximum clutters the figure; a maximum that
        # only sits at the end of the traced arc is not a fold, so it is
        # labelled only when it is also the overall maximum.
        if interior or curve.cid == global_max_cid:
            # Flip the label to the left near the right edge so it stays
            # inside the figure.
            lam_lo, lam_hi = ax.get_xlim()
            span = lam_hi - lam_lo
            on_right = span > 0 and (curve.lam[idx] - lam_lo) / span > 0.7
            ax.annotate(
                "{}$\\lambda$={:.6g}".format(
                    "fold " if interior else "max ", curve.lam[idx]
                ),
                xy=(float(curve.lam[idx]), float(curve.norm_u[idx])),
                xytext=(-8 if on_right else 8, 10 if k % 2 == 0 else -14),
                textcoords="offset points",
                horizontalalignment="right" if on_right else "left",
                fontsize=8,
                color=color,
                bbox=dict(
                    boxstyle="round,pad=0.2",
                    facecolor="white",
                    edgecolor="none",
                    alpha=0.7,
                ),
                arrowprops=dict(arrowstyle="-", color=color, linewidth=0.6),
            )

    # Bifurcation points, drawn on top and with a single shared legend entry.
    bif_lam = np.concatenate(
        [c.lam[c.bifurcation] for c in curves] + [np.empty(0)]
    )
    bif_u = np.concatenate(
        [c.norm_u[c.bifurcation] for c in curves] + [np.empty(0)]
    )
    if bif_lam.size:
        ax.plot(
            bif_lam,
            bif_u,
            marker="*",
            markersize=13,
            linestyle="none",
            color="k",
            zorder=6,
            label="bifurcation",
        )

    if reference is not None:
        ax.axvline(
            reference,
            color="crimson",
            linewidth=1.0,
            linestyle="-.",
            zorder=2,
            label="reference $\\lambda$={:.6g}".format(reference),
        )

    # Style legend entries (proxy artists, no data).
    for flag in (1, -1, 0):
        if any(np.any(c.stability == flag) for c in curves):
            ax.plot(
                [],
                [],
                color="0.25",
                linestyle=_STABILITY_STYLE[flag],
                linewidth=1.6,
                label=_STABILITY_NAME[flag],
            )

    ax.set_xlabel(r"load parameter $\lambda$")
    ax.set_ylabel(r"$\|u\|$")
    ax.set_title(title if title else "Equilibrium landscape")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=8, framealpha=0.9)
    fig.tight_layout()

    try:
        fig.savefig(output, dpi=dpi)
    except (OSError, ValueError) as exc:
        plt.close(fig)
        raise LandscapeError("cannot write '{}': {}".format(output, exc))
    plt.close(fig)
    return output


def _default_output(csv_path: str) -> str:
    """Derive the default image path from the CSV path.

    Parameters
    ----------
    csv_path : str
        Path to the input CSV.

    Returns
    -------
    str
        ``<csv without extension>.png``, beside the CSV.
    """
    base, _ = os.path.splitext(csv_path)
    return base + ".png"


# ---------------------------------------------------------------------------
# Locus-level thesis-comparison mode (--merge-loci / --thesis-metrics)
# ---------------------------------------------------------------------------

#: Colour per solution locus, following the thesis figure (A navy,
#: B yellow-green, C pink/violet); unclassified points are grey.
_LOCUS_COLOR: Dict[str, str] = {
    "A": "navy",
    "B": "yellowgreen",
    "C": "orchid",
    "other": "0.5",
}

#: Deterministic locus iteration order.
_LOCUS_ORDER: Tuple[str, ...] = ("A", "B", "C", "other")

#: Relative spread (``(max(U)-min(U)) / max|U|`` over the landscape) below
#: which a state counts as spatially constant, i.e. locus A.  Derived from
#: the measured ModifiedBratu landscapes: constant-branch spreads are at
#: most 1.4e-10 of the state scale while the smallest non-constant
#: (emerging mode-1) spread is 6.6e-2 of it; 1e-5 sits more than three
#: decades away from both measured bounds (log-midpoint 3e-6).
_SPREAD_REL_TOL: float = 1e-5

#: Relative magnitude below which a coefficient of ``U - mean(U)`` is
#: treated as numerical noise when counting sign changes.  On constant
#: states the residual is pure round-off (measured <= 1.4e-10 of the state
#: scale); on the smallest traced mode-1 state the residual amplitude is
#: 6.6e-2 of it.  1e-6 sits >= 3.9 decades from either measured bound and
#: is the value the task-59/60/61 review used.
_SIGN_REL_TOL: float = 1e-6

#: Minimum ratio across a gap in the sorted step-length distribution for
#: that gap to count as separating solver jumps from genuine continuation
#: steps.  Measured over four landscapes: the largest consecutive ratio
#: between genuine (upper-tail) steps is 1.81, the smallest jump gap ratio
#: is 4.47; 3.0 separates the two with margin on both sides.
_JUMP_GAP_RATIO: float = 3.0

#: State-space distance below which two merged points count as duplicates.
#: The smallest measured distance between *distinct* landscape points is
#: 1.94e-4 (the stalled-fold points of the default run's curve 1), so 1e-5
#: can only ever drop genuinely repeated points.
_DEDUP_TOL: float = 1e-5

#: Tick-label values of Wouters (2019) fig. 9.3, read from the rendered
#: figure itself (400 dpi render of physical page 284 of the thesis PDF).
_FIG93_XTICK_VALUES: Tuple[float, ...] = (
    0.01, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40,
)
#: The axis title (a mu glyph) sits between the '0.15' and '0.2' labels and
#: shows up as one extra cluster at this position in the label row.
_FIG93_XTICK_MU_TEXT_INDEX: int = 4
_FIG93_YTICK_VALUES: Tuple[float, ...] = (6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0)

#: Spline degree of the driver's discretization.  The landscape HDF5 stores
#: only the coefficient matrices -- it carries NO basis description -- so the
#: basis has to be reconstructed.  ``example_ModifiedBratuExploration``
#: builds a degree-1 unit-interval geometry and elevates it once
#: (``numElevate = 1``), i.e. degree 2, then refines it (``numHref = 5``) to
#: 32 uniform elements of an open knot vector on [0, 1], giving
#: ``N_dof = n_elem + degree = 34`` -- which is the value measured in every
#: landscape export, and which the driver itself prints
#: (``Patches: 1, degree: 2, elements: 32`` / ``Number of free DoFs: 34``;
#: the Neumann problem eliminates no DoF, so the exported vector *is* the
#: coefficient vector).  ``--basis-degree`` overrides this assumption; the
#: element count is always derived as ``N_dof - degree``.
_BASIS_DEGREE: int = 2

#: State-norm conventions selectable with ``--norm``.
_STATE_NORMS: Tuple[str, ...] = ("l2", "h1", "discrete", "thesis")

#: Number of finite-difference grid points of the thesis's own
#: discretization (Wouters 2019, printed p. 258), used by the ``thesis``
#: norm convention: ``<.,.>_h = (1/(n-1)) <.,.>_2`` (eq. 2.8, printed p. 17).
_THESIS_FD_POINTS: int = 100


def open_uniform_knots(degree: int, n_elem: int) -> np.ndarray:
    """Open uniform knot vector of ``n_elem`` elements on [0, 1].

    Parameters
    ----------
    degree : int
        Spline degree.
    n_elem : int
        Number of (uniform) elements.

    Returns
    -------
    numpy.ndarray
        Knot vector of length ``n_elem + 2*degree + 1``; the first and last
        knot are repeated ``degree + 1`` times.
    """
    inner = np.linspace(0.0, 1.0, int(n_elem) + 1)
    return np.concatenate(
        [np.zeros(int(degree)), inner, np.ones(int(degree))]
    )


def bspline_basis(
    knots: np.ndarray, degree: int, x: np.ndarray, deriv: int = 0
) -> np.ndarray:
    """Evaluate all B-spline basis functions (or their first derivative).

    Plain Cox--de Boor recursion; no SciPy, keeping the module's stated
    dependency set (standard library, NumPy, Matplotlib).  The result was
    checked against ``scipy.interpolate.BSpline`` on this space and agrees
    to 0.0 (bit-identical) for both ``deriv = 0`` and ``deriv = 1``.

    Parameters
    ----------
    knots : numpy.ndarray
        Knot vector, non-decreasing.
    degree : int
        Spline degree ``p``.
    x : numpy.ndarray
        Evaluation points, inside ``[knots[0], knots[-1]]``.
    deriv : int, optional
        ``0`` for the basis functions, ``1`` for their first derivative.

    Returns
    -------
    numpy.ndarray
        Shape ``(x.size, N_dof)`` with ``N_dof = len(knots) - degree - 1``.
    """
    x = np.atleast_1d(np.asarray(x, dtype=float))
    n_dof = int(knots.size) - int(degree) - 1
    # Degree 0: the indicator of every non-empty knot span (the last span
    # is closed on the right so x = 1 is covered).
    basis = np.zeros((x.size, int(knots.size) - 1))
    for i in range(int(knots.size) - 1):
        if knots[i + 1] > knots[i]:
            inside = (x >= knots[i]) & (x < knots[i + 1])
            if knots[i + 1] == knots[-1]:
                inside = inside | (x == knots[-1])
            basis[:, i] = inside
    lower = basis
    for p in range(1, int(degree) + 1):
        raised = np.zeros((x.size, int(knots.size) - 1 - p))
        for i in range(raised.shape[1]):
            d1 = knots[i + p] - knots[i]
            d2 = knots[i + p + 1] - knots[i + 1]
            left = ((x - knots[i]) / d1) * basis[:, i] if d1 > 0.0 else 0.0
            right = (
                ((knots[i + p + 1] - x) / d2) * basis[:, i + 1]
                if d2 > 0.0
                else 0.0
            )
            raised[:, i] = left + right
        lower = basis
        basis = raised
    if deriv == 0:
        return basis[:, :n_dof]
    # B'_{i,p} = p * ( B_{i,p-1}/(t_{i+p} - t_i)
    #                  - B_{i+1,p-1}/(t_{i+p+1} - t_{i+1}) )
    out = np.zeros((x.size, n_dof))
    p = int(degree)
    for i in range(n_dof):
        d1 = knots[i + p] - knots[i]
        d2 = knots[i + p + 1] - knots[i + 1]
        a = (lower[:, i] / d1) if d1 > 0.0 else 0.0
        b = (lower[:, i + 1] / d2) if d2 > 0.0 else 0.0
        out[:, i] = p * (a - b)
    return out


def spline_gram_matrices(
    n_dof: int, degree: int = _BASIS_DEGREE
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Mass and stiffness Gram matrices of the driver's spline space.

    ``M[i, j] = \\int_0^1 B_i B_j`` and ``K[i, j] = \\int_0^1 B_i' B_j'``
    are assembled element by element with Gauss--Legendre quadrature of
    ``degree + 1`` points, which is exact for polynomials up to degree
    ``2*degree + 1``.  The integrand ``B_i B_j`` is a polynomial of degree
    ``2*degree`` on every element (the basis is piecewise polynomial with
    breakpoints exactly at the element boundaries), so **the rule is exact
    and both matrices are exact to round-off** -- no quadrature error
    enters the norm.  ``B_i' B_j'`` has degree ``2*degree - 2`` and is
    covered a fortiori.

    Parameters
    ----------
    n_dof : int
        Number of coefficients per solution vector, read from the HDF5.
        The element count is ``n_dof - degree`` (open uniform knots).
    degree : int, optional
        Spline degree; see :data:`_BASIS_DEGREE` for the assumption.

    Returns
    -------
    tuple of numpy.ndarray
        ``(M, K, knots)``.

    Raises
    ------
    LandscapeError
        If the requested degree cannot describe ``n_dof`` coefficients.
    """
    degree = int(degree)
    n_dof = int(n_dof)
    if degree < 1:
        raise LandscapeError(
            "basis degree must be >= 1, got {}".format(degree)
        )
    n_elem = n_dof - degree
    if n_elem < 1:
        raise LandscapeError(
            "N_dof = {} is too small for a degree-{} basis (needs at least "
            "{} coefficients) -- wrong --basis-degree?".format(
                n_dof, degree, degree + 1
            )
        )
    knots = open_uniform_knots(degree, n_elem)
    nodes, weights = np.polynomial.legendre.leggauss(degree + 1)
    edges = np.linspace(0.0, 1.0, n_elem + 1)
    mass = np.zeros((n_dof, n_dof))
    stiff = np.zeros((n_dof, n_dof))
    for e in range(n_elem):
        a, b = float(edges[e]), float(edges[e + 1])
        xq = 0.5 * (b - a) * nodes + 0.5 * (a + b)
        wq = 0.5 * (b - a) * weights
        bq = bspline_basis(knots, degree, xq, 0)
        dq = bspline_basis(knots, degree, xq, 1)
        mass += bq.T @ (wq[:, None] * bq)
        stiff += dq.T @ (wq[:, None] * dq)
    return mass, stiff, knots


def _admissible_degrees(n_dof: int) -> List[int]:
    """Spline degrees whose element count is a power of two.

    See :func:`check_basis_degree` for the derivation.

    Parameters
    ----------
    n_dof : int
        Number of coefficients per solution vector.

    Returns
    -------
    list of int
        Ascending; empty only for ``n_dof < 2``.
    """
    out: List[int] = []
    for p in range(1, int(n_dof)):
        n_elem = int(n_dof) - p
        if n_elem >= 1 and (n_elem & (n_elem - 1)) == 0:
            out.append(p)
    return out


def check_basis_degree(n_dof: int, degree: int) -> List[str]:
    """Validate an assumed spline degree against the knot structure.

    The degree is the one *un-measured* input of the ``l2``/``h1``/
    ``thesis`` norms -- the landscape HDF5 stores no basis description --
    and it is not harmless: reading a real locus-C state with degree 3
    instead of degree 2 moves its L2 norm from 6.390882 to 6.311828,
    i.e. by 1.24 % (measured).  The data does constrain it, though.
    ``example_ModifiedBratuExploration`` adds a single-element
    ``BSplineUnitInterval``, raises it to degree ``1 + numElevate`` and
    then calls ``dbasis.uniformRefine()`` ``numHref`` times, and every
    uniform refinement *halves* each element, so the export always
    satisfies

    ``N_dof = 2**numHref + degree``

    -- the element count ``N_dof - degree`` must be a **power of two**.
    That bound is derived from the driver's own construction, not tuned:
    it rejects exactly the degrees the export cannot have come from.  It
    catches both hazards seen in practice: degree 3 on a 34-dof export
    (31 elements, the 1.24 % case above) and the default degree 2 on a
    35-dof export produced with ``-e 2`` (33 elements), which is
    otherwise perfectly legal arithmetic.

    Several degrees can remain admissible (``N_dof = 34`` admits 2, 18,
    26, 30, 32 and 33, i.e. ``numHref = 5, 4, 3, 2, 1, 0``).  The
    remaining ambiguity is *not* decidable from the export, so a
    non-smallest choice is reported as a note rather than refused: a
    landscape genuinely run with few refinements and a high degree is
    legitimate.

    The rule also refuses the *two-dimensional* landscape of
    ``example_BratuExploration`` (unit square, Dirichlet-eliminated
    16 x 16 = 256 free dofs -- measured on its own export), which has no
    one-dimensional reconstruction to begin with; there not one
    admissible degree leaves more elements than its own order, which is
    what the message then says.  ``--norm discrete`` needs no basis and
    stays available on such an export.

    Parameters
    ----------
    n_dof : int
        Number of coefficients per solution vector, read from the HDF5.
    degree : int
        Assumed spline degree (``--basis-degree``).

    Returns
    -------
    list of str
        Diagnostic note lines (empty when the degree is the smallest
        admissible one, which is the driver's default configuration).

    Raises
    ------
    LandscapeError
        If the degree is below 1, leaves fewer than one element, or
        implies an element count that is not a power of two.
    """
    degree = int(degree)
    n_dof = int(n_dof)
    if degree < 1:
        raise LandscapeError(
            "basis degree must be >= 1, got {}".format(degree)
        )
    n_elem = n_dof - degree
    if n_elem < 1:
        raise LandscapeError(
            "N_dof = {} is too small for a degree-{} basis (needs at least "
            "{} coefficients) -- wrong --basis-degree?".format(
                n_dof, degree, degree + 1
            )
        )
    admissible = _admissible_degrees(n_dof)
    if n_elem & (n_elem - 1):
        # When not one admissible degree leaves more elements than its own
        # order, no 1D reconstruction of this export is credible at all.
        hint = ""
        if not [p for p in admissible if (n_dof - p) > p]:
            hint = (
                "; no admissible degree leaves more elements than its own "
                "order, so this is most likely not a 1D export (the 2D "
                "example_BratuExploration landscape has 16x16 = 256 free "
                "dofs) -- use --norm discrete, which needs no basis"
            )
        listed = ", ".join(str(p) for p in admissible)
        raise LandscapeError(
            "basis degree {} implies {} elements for N_dof = {}, but the "
            "driver refines a single element uniformly (numHref doublings), "
            "so N_dof - degree must be a power of two; this export admits "
            "degree(s) {} -- and the degree is not cosmetic: a wrong one "
            "shifted a measured locus-C state by 1.24 % (6.390882 vs "
            "6.311828){}".format(
                degree, n_elem, n_dof, listed if listed else "none", hint
            )
        )
    notes: List[str] = []
    if admissible and degree != admissible[0]:
        notes.append(
            "NOTE: degree {} means {} element(s) (numHref = {}); the "
            "smallest degree N_dof = {} admits is {} ({} elements), which "
            "is the driver's default configuration -- the export cannot "
            "distinguish them, so check it against the run log".format(
                degree,
                n_elem,
                int(n_elem).bit_length() - 1,
                n_dof,
                admissible[0],
                n_dof - admissible[0],
            )
        )
    return notes


def state_norm(
    u: np.ndarray, kind: str = "l2", degree: int = _BASIS_DEGREE
) -> np.ndarray:
    """State-space norm of one or more solution coefficient vectors.

    The comparison quantity of the thesis figures is a norm of the
    *solution function* ``psi``, not of its coefficient vector.  Four
    conventions are available:

    ``l2``
        The continuous norm ``sqrt(int_0^1 psi^2) = sqrt(U^T M U)`` with
        the exact mass matrix of :func:`spline_gram_matrices`.  This is
        the default of the thesis-comparison mode.
    ``h1``
        ``sqrt(U^T (M + K) U) = sqrt(int psi^2 + int psi'^2)``, the full
        H1 norm (not the semi-norm).
    ``discrete``
        The project's historical quantity ``||U||_2 / sqrt(N_dof)``.  It
        is a quadrature rule with *uniform* weights on the Greville
        abscissae, whose two end spacings are half the interior spacing,
        so the two endpoints are over-weighted by about 2x.  It is exact
        on constant states and biased on every non-constant one; kept for
        back-comparison with earlier reports.
    ``thesis``
        The thesis's own discrete convention ``sqrt((1/(n-1)) sum psi_i^2)``
        on ``n =`` :data:`_THESIS_FD_POINTS` uniform points (eq. 2.8), but
        evaluated on *our* spline field.  It quantifies the residual
        convention mismatch against the thesis's y-axis; it is **not** the
        thesis's own solution (a different discretization), so it bounds
        the convention gap only.

    On a spatially constant state ``psi = c`` the first three agree
    exactly with ``c``; ``thesis`` gives ``c * sqrt(n/(n-1))``.

    Parameters
    ----------
    u : numpy.ndarray
        Coefficient vector of shape ``(N_dof,)`` or a matrix of shape
        ``(N_dof, n_points)`` (the HDF5 layout).
    kind : str, optional
        One of :data:`_STATE_NORMS`.
    degree : int, optional
        Spline degree; see :data:`_BASIS_DEGREE`.

    Returns
    -------
    numpy.ndarray
        One norm value per column (shape ``(n_points,)``; shape ``(1,)``
        for a single vector).

    Raises
    ------
    LandscapeError
        On an unknown ``kind`` or an unusable ``(N_dof, degree)`` pair.
    """
    if kind not in _STATE_NORMS:
        raise LandscapeError(
            "unknown state norm '{}' -- expected one of {}".format(
                kind, ", ".join(_STATE_NORMS)
            )
        )
    u = np.asarray(u, dtype=float)
    if u.ndim == 1:
        u = u[:, None]
    if u.ndim != 2:
        raise LandscapeError(
            "state_norm expects a vector or an (N_dof, n) matrix"
        )
    n_dof = int(u.shape[0])
    if kind == "discrete":
        return np.linalg.norm(u, axis=0) / np.sqrt(float(n_dof))
    if kind == "thesis":
        _mass, _stiff, knots = spline_gram_matrices(n_dof, degree)
        xs = np.linspace(0.0, 1.0, _THESIS_FD_POINTS)
        field = bspline_basis(knots, degree, xs, 0) @ u
        return np.sqrt(
            np.sum(field * field, axis=0) / float(_THESIS_FD_POINTS - 1)
        )
    mass, stiff, _knots = spline_gram_matrices(n_dof, degree)
    gram = mass if kind == "l2" else mass + stiff
    quad = np.einsum("ij,ik,kj->j", u, gram, u)
    return np.sqrt(np.maximum(quad, 0.0))


def norm_label(kind: str) -> str:
    """One-line description of a state-norm convention, for axes/reports.

    Parameters
    ----------
    kind : str
        One of :data:`_STATE_NORMS`.

    Returns
    -------
    str
        Human readable label.

    Raises
    ------
    LandscapeError
        On an unknown ``kind``.
    """
    labels = {
        "l2": "||psi||_L2 = sqrt(int_0^1 psi^2)",
        "h1": "||psi||_H1 = sqrt(int psi^2 + int psi'^2)",
        "discrete": "||U||_2 / sqrt(N_dof)  (biased off constant states)",
        "thesis": (
            "sqrt((1/(n-1)) sum psi(x_i)^2), n = {} uniform "
            "points".format(_THESIS_FD_POINTS)
        ),
    }
    if kind not in labels:
        raise LandscapeError(
            "unknown state norm '{}' -- expected one of {}".format(
                kind, ", ".join(_STATE_NORMS)
            )
        )
    return labels[kind]


def _norm_axis_label(kind: str) -> str:
    """Matplotlib y-axis label (TeX) for a state-norm convention.

    Parameters
    ----------
    kind : str
        One of :data:`_STATE_NORMS`.

    Returns
    -------
    str
        TeX label string.
    """
    return {
        "l2": r"$\|\psi\|_{L^2} = \sqrt{\int_0^1 \psi^2}$",
        "h1": r"$\|\psi\|_{H^1}$",
        "discrete": r"$\|\psi\| = \|U\|_2/\sqrt{N_\mathrm{dof}}$",
        "thesis": r"$\|\psi\|_h$ (thesis convention, $n=100$)",
    }.get(kind, r"$\|\psi\|$")


def _const_exactness_note(kind: str) -> str:
    """Whether ``kind`` reproduces a constant state exactly, as a phrase.

    ``l2``, ``h1`` and ``discrete`` all return ``c`` on ``psi = c``,
    which is why the locus-A comparison is norm-independent; ``thesis``
    does **not** -- it returns ``c*sqrt(n/(n-1))``.  Stating this
    unconditionally would contradict the warning printed above the same
    table.

    Parameters
    ----------
    kind : str
        One of :data:`_STATE_NORMS`.

    Returns
    -------
    str
        Parenthetical phrase, without surrounding brackets.
    """
    if kind == "thesis":
        return (
            "this norm is NOT exact on constant states: it reads "
            "c*sqrt(n/(n-1))"
        )
    return "the '{}' norm is exact on constant states".format(kind)


def _ref_precision_note(kind: str) -> str:
    """Trailing qualifier for a deviation against the digitized figure.

    The ``~0.1`` read-off floor of the figure is a statement about an
    L2-type y-axis.  Under ``h1`` the rows carry the extra
    ``int psi'^2`` term and are not a figure comparison at all, so the
    floor is meaningless there and must not be printed.

    Parameters
    ----------
    kind : str
        One of :data:`_STATE_NORMS`.

    Returns
    -------
    str
        Text appended to the deviation line (may be empty).
    """
    if kind == "h1":
        return (
            "(no read-off floor applies: H1 is not the figure's axis, "
            "see the WARNING above)"
        )
    return "(reference precision ~0.1)"


def read_solution_vectors(
    h5_path: str, curves: Sequence[Curve]
) -> Dict[int, np.ndarray]:
    """Read the per-curve solution vectors from a landscape HDF5 export.

    The export stores one dataset ``curve_%04d_U`` of shape
    ``(N_dof, n_points)`` per curve (plus per-point scalars and a
    ``structure`` table, which are not needed here -- the CSV carries them).

    Parameters
    ----------
    h5_path : str
        Path to ``landscape.h5`` (same basename as the CSV).
    curves : sequence of Curve
        Curves as returned by :func:`read_landscape`; used to check that
        every curve has a matching dataset of the right length.

    Returns
    -------
    dict of int to numpy.ndarray
        ``curve id -> U`` with ``U`` of shape ``(N_dof, n_points)``.

    Raises
    ------
    LandscapeError
        If ``h5py`` is unavailable, the file cannot be read, a dataset is
        missing, or a dataset's point count disagrees with the CSV.
    """
    if not _HAS_H5PY:
        raise LandscapeError(
            "h5py is not installed; --merge-loci needs the HDF5 solution "
            "vectors (a norm-only classification cannot separate the loci)"
        )
    try:
        handle = h5py.File(h5_path, "r")
    except OSError as exc:
        raise LandscapeError("cannot read '{}': {}".format(h5_path, exc))
    vectors: Dict[int, np.ndarray] = {}
    with handle:
        for curve in curves:
            name = "curve_{:04d}_U".format(curve.cid)
            if name not in handle:
                raise LandscapeError(
                    "'{}' has no dataset '{}' -- not a landscape HDF5 "
                    "export matching the CSV".format(h5_path, name)
                )
            u = np.asarray(handle[name][...], dtype=float)
            if u.ndim != 2 or u.shape[1] != curve.size:
                raise LandscapeError(
                    "dataset '{}' has shape {} but the CSV curve has {} "
                    "point(s) -- .csv/.h5 pair out of sync?".format(
                        name, u.shape, curve.size
                    )
                )
            vectors[curve.cid] = u
    return vectors


def count_sign_changes(u: np.ndarray, rel_tol: float = _SIGN_REL_TOL) -> int:
    """Count the sign changes of ``u - mean(u)``, ignoring round-off noise.

    Coefficients with ``|u - mean(u)| <= rel_tol * max|u|`` are skipped
    before counting, so a numerically constant vector reports zero sign
    changes instead of the one-to-three a naive counter returns on pure
    round-off.  Note the count runs over B-spline *coefficients*, which by
    the variation-diminishing property only bound the function's sign
    changes from above; at the traced mode amplitudes this is decisive.

    Parameters
    ----------
    u : numpy.ndarray
        Solution coefficient vector.  Shape ``(N_dof,)``.
    rel_tol : float, optional
        Noise threshold relative to ``max|u|``; see :data:`_SIGN_REL_TOL`
        for the measured derivation of the default.

    Returns
    -------
    int
        Number of sign changes of the thresholded residual.
    """
    u = np.asarray(u, dtype=float)
    if u.size == 0:
        return 0
    resid = u - float(u.mean())
    scaleu = float(np.abs(u).max())
    if scaleu == 0.0:
        return 0
    signs = np.sign(np.where(np.abs(resid) > rel_tol * scaleu, resid, 0.0))
    signs = signs[signs != 0.0]
    if signs.size < 2:
        return 0
    return int(np.count_nonzero(signs[1:] != signs[:-1]))


def classify_locus(
    u: np.ndarray,
    scale: float,
    spread_rel_tol: float = _SPREAD_REL_TOL,
    sign_rel_tol: float = _SIGN_REL_TOL,
) -> str:
    """Classify one solution vector into a locus by its spatial profile.

    The test is the one measured decisive for the ModifiedBratu landscape
    (norm-only tests were shown to fail there):

    * ``spread = max(u) - min(u)`` at round-off level -> constant state,
      locus ``"A"``;
    * one sign change of ``u - mean(u)`` -> mode 1, locus ``"B"``;
    * two sign changes -> mode 2, locus ``"C"``;
    * anything else -> ``"other"``.

    Parameters
    ----------
    u : numpy.ndarray
        Solution coefficient vector.  Shape ``(N_dof,)``.
    scale : float
        ``max|U|`` over the whole landscape; referencing the *global* scale
        keeps near-zero states on the trivial branch classified as
        constant.  Guarded below by 1 so an all-zero landscape does not
        produce a zero threshold.
    spread_rel_tol, sign_rel_tol : float, optional
        See :data:`_SPREAD_REL_TOL` and :data:`_SIGN_REL_TOL`.

    Returns
    -------
    str
        One of ``"A"``, ``"B"``, ``"C"``, ``"other"``.
    """
    u = np.asarray(u, dtype=float)
    spread = float(u.max() - u.min()) if u.size else 0.0
    if spread <= spread_rel_tol * max(scale, 1.0):
        return "A"
    n = count_sign_changes(u, rel_tol=sign_rel_tol)
    if n == 1:
        return "B"
    if n == 2:
        return "C"
    return "other"


def derive_jump_cut(
    steps: np.ndarray,
) -> Tuple[Optional[float], str]:
    """Derive the solver-jump step-length cut from the data.

    A solver jump is a single continuation step whose state-space length is
    an outlier of the step distribution.  The cut is placed at the
    geometric mean across the largest consecutive gap (in ratio) of the
    sorted upper tail (steps at or above the median), and only accepted
    when that ratio is at least :data:`_JUMP_GAP_RATIO` -- below it the
    distribution shows no jump and nothing is cut.

    Parameters
    ----------
    steps : numpy.ndarray
        Pooled state-space step lengths of all curves.

    Returns
    -------
    tuple of (float or None, str)
        The cut (``None`` when the distribution shows no jump gap) and a
        one-line description of how it was derived.
    """
    d = np.sort(np.asarray(steps, dtype=float))
    d = d[d > 0.0]
    if d.size < 4:
        return None, "none (fewer than 4 positive steps)"
    upper = d[d >= float(np.median(d))]
    if upper.size < 2:
        return None, "none (degenerate step distribution)"
    ratios = upper[1:] / upper[:-1]
    k = int(np.argmax(ratios))
    if float(ratios[k]) < _JUMP_GAP_RATIO:
        return None, (
            "none (largest upper-tail gap ratio {:.3g} < {:g})".format(
                float(ratios[k]), _JUMP_GAP_RATIO
            )
        )
    cut = float(np.sqrt(upper[k] * upper[k + 1]))
    return cut, (
        "{:.6g} (geometric mean across the largest upper-tail step gap, "
        "{:.6g} -> {:.6g}, ratio {:.3g})".format(
            cut, float(upper[k]), float(upper[k + 1]), float(ratios[k])
        )
    )


def _chain_order(mu: np.ndarray, state: np.ndarray) -> np.ndarray:
    """Order pooled points along their branch by nearest-neighbour chaining.

    Ordering by ``mu`` alone is not enough -- locus A folds -- so the chain
    is grown greedily in the ``(mu, state)`` plane, starting from the point
    farthest from the centroid (a branch extremity).

    Parameters
    ----------
    mu, state : numpy.ndarray
        Point coordinates.  Shape ``(n,)`` each.

    Returns
    -------
    numpy.ndarray
        Permutation of ``arange(n)`` giving the chained order.
    """
    n = int(mu.size)
    if n <= 2:
        return np.arange(n, dtype=int)
    cx, cy = float(mu.mean()), float(state.mean())
    start = int(np.argmax(np.hypot(mu - cx, state - cy)))
    order = [start]
    remaining = np.ones(n, dtype=bool)
    remaining[start] = False
    cur = start
    for _ in range(n - 1):
        rem = np.flatnonzero(remaining)
        d = np.hypot(mu[rem] - mu[cur], state[rem] - state[cur])
        cur = int(rem[int(np.argmin(d))])
        order.append(cur)
        remaining[cur] = False
    return np.asarray(order, dtype=int)


def merge_loci(
    curves: Sequence[Curve],
    vectors: Dict[int, np.ndarray],
    window: Optional[float] = None,
    norm: str = "l2",
    degree: int = _BASIS_DEGREE,
) -> Tuple[Dict[str, Dict[str, object]], List[str]]:
    """Merge all traced segments of each solution locus into one polyline.

    Pipeline: (1) clip to ``mu >= window``; (2) split every curve at
    solver-jump links (cut derived by :func:`derive_jump_cut` from the
    pooled step distribution); (3) classify every point via
    :func:`classify_locus`; (4) pool per locus, drop duplicates within
    :data:`_DEDUP_TOL`, chain via :func:`_chain_order`, and re-split the
    chained polyline at links longer than the jump cut (or, when no jump
    was detected, twice the largest surviving step -- adjacent samples of
    one branch cannot be farther apart than that).

    What is left to merge
    ---------------------
    Part of the segmentation this function was written for is now done
    upstream: the explorer sweeps **both arc-length directions of a seed
    into one curve**, so a locus traced forwards and backwards from the
    same seed no longer arrives as two curves (a pristine rest state is
    still swept in one direction only).  What still arrives split, and is
    what merging is now for:

    * **branch legs of the same locus from different branch points** --
      one locus is reached by several branch curves (one job per mode
      sign, and several branch points along the parent), each exported as
      its own curve;
    * **solver jumps** -- a single traced curve can step across a gap in
      state space, which is cut here and never drawn as a link;
    * the **overlaps** those legs leave behind, which are dropped as
      duplicates within :data:`_DEDUP_TOL`.

    Parameters
    ----------
    curves : sequence of Curve
        Curves as returned by :func:`read_landscape`.
    vectors : dict of int to numpy.ndarray
        Solution vectors per curve id, from :func:`read_solution_vectors`.
    window : float, optional
        Keep only points with ``lambda >= window`` (the thesis's stopping
        condition, eq. 9.1); ``None`` keeps everything.
    norm : str, optional
        State-norm convention, see :func:`state_norm`; ``"l2"`` (the
        continuous norm of the spline solution) by default.
    degree : int, optional
        Spline degree of the reconstructed basis; see
        :data:`_BASIS_DEGREE`.

    Returns
    -------
    tuple of (dict, list of str)
        ``loci`` maps a locus label to a dict with the chained arrays
        ``mu``, ``state`` (in the requested norm), ``stability``,
        ``bif``, the polyline ``parts`` (index arrays), and ``n_dupes``;
        the second element is a list of diagnostic note lines.

    Raises
    ------
    LandscapeError
        On missing/mismatched vectors, or when nothing survives the clip.
    """
    notes: List[str] = []
    if not curves:
        raise LandscapeError("no curves to merge")
    n_dof: Optional[int] = None
    for curve in curves:
        u = vectors.get(curve.cid)
        if u is None:
            raise LandscapeError(
                "no solution vectors for curve {}".format(curve.cid)
            )
        if n_dof is None:
            n_dof = int(u.shape[0])
        elif int(u.shape[0]) != n_dof:
            raise LandscapeError(
                "inconsistent N_dof across curves ({} vs {})".format(
                    int(u.shape[0]), n_dof
                )
            )
    assert n_dof is not None
    scale = max(
        float(np.abs(vectors[c.cid]).max()) if vectors[c.cid].size else 0.0
        for c in curves
    )
    notes.append(
        "N_dof = {} (from the HDF5 solution vectors), state scale "
        "max|U| = {:.6g}".format(n_dof, scale)
    )
    notes.append("state norm: {} [{}]".format(norm, norm_label(norm)))
    if norm != "discrete":
        # The HDF5 carries no basis description; state the reconstruction
        # -- and refuse a degree the export cannot have come from.
        degree_notes = check_basis_degree(n_dof, degree)
        notes.append(
            "basis (assumed, reconstructed from N_dof -- the HDF5 stores no "
            "basis): degree {}, {} uniform elements, open knots on "
            "[0,1]".format(int(degree), n_dof - int(degree))
        )
        notes.extend(degree_notes)
    elif int(degree) != _BASIS_DEGREE:
        notes.append(
            "note: --basis-degree {} is unused by the 'discrete' norm, "
            "which needs no basis".format(int(degree))
        )

    # Per-point state coordinate and locus label; window clip into
    # contiguous runs so no step is taken across a removed point.
    runs: List[Tuple[np.ndarray, ...]] = []
    n_clipped = 0
    for curve in curves:
        u_all = vectors[curve.cid]
        state = state_norm(u_all, kind=norm, degree=degree)
        # The CSV's normU is ||U||_2 at ~6 significant figures; disagreement
        # beyond that means the .csv/.h5 pair is out of sync.  This check is
        # deliberately independent of the plotted norm.
        denom = np.maximum(np.abs(curve.norm_u), 1.0)
        mismatch = np.abs(np.linalg.norm(u_all, axis=0) - curve.norm_u) / denom
        if float(mismatch.max()) > 1e-4:
            raise LandscapeError(
                "curve {}: CSV normU and HDF5 ||U||_2 disagree by {:.3g} "
                "(relative) -- .csv/.h5 pair out of sync?".format(
                    curve.cid, float(mismatch.max())
                )
            )
        labels = [
            classify_locus(u_all[:, j], scale) for j in range(curve.size)
        ]
        keep = (
            np.ones(curve.size, dtype=bool)
            if window is None
            else curve.lam >= window
        )
        n_clipped += int(np.count_nonzero(~keep))
        j = 0
        while j < curve.size:
            if not keep[j]:
                j += 1
                continue
            k = j
            while k < curve.size and keep[k]:
                k += 1
            runs.append(
                (
                    curve.lam[j:k],
                    state[j:k],
                    curve.stability[j:k],
                    curve.bifurcation[j:k],
                    labels[j:k],
                )
            )
            j = k
    if window is not None:
        notes.append(
            "window: clipped {} point(s) below mu = {:g}".format(
                n_clipped, window
            )
        )
    if not runs:
        raise LandscapeError("no points remain after the window clip")

    steps = np.concatenate(
        [np.hypot(np.diff(r[0]), np.diff(r[1])) for r in runs]
        + [np.empty(0)]
    )
    cut, cut_note = derive_jump_cut(steps)
    notes.append("jump cut: " + cut_note)

    subruns: List[Tuple[np.ndarray, ...]] = []
    n_jump_links = 0
    for mu, st, stab, bif, labels in runs:
        if cut is None or mu.size < 2:
            subruns.append((mu, st, stab, bif, labels))
            continue
        d = np.hypot(np.diff(mu), np.diff(st))
        breaks = np.flatnonzero(d > cut)
        n_jump_links += int(breaks.size)
        j = 0
        for b in breaks:
            subruns.append(
                (
                    mu[j : b + 1],
                    st[j : b + 1],
                    stab[j : b + 1],
                    bif[j : b + 1],
                    labels[j : b + 1],
                )
            )
            j = int(b) + 1
        subruns.append((mu[j:], st[j:], stab[j:], bif[j:], labels[j:]))
    if cut is not None:
        notes.append(
            "dropped {} solver-jump link(s) longer than the cut".format(
                n_jump_links
            )
        )

    pools: Dict[str, List[Tuple[float, float, int, bool]]] = {}
    for mu, st, stab, bif, labels in subruns:
        for j, lab in enumerate(labels):
            pools.setdefault(lab, []).append(
                (float(mu[j]), float(st[j]), int(stab[j]), bool(bif[j]))
            )

    kept_steps = steps[steps <= cut] if cut is not None else steps
    max_step = float(kept_steps.max()) if kept_steps.size else 0.0
    link_cut = cut if cut is not None else (
        2.0 * max_step if max_step > 0.0 else None
    )
    if cut is None and link_cut is not None:
        notes.append(
            "polyline link cut: {:.6g} (no jump detected; twice the "
            "largest measured step)".format(link_cut)
        )

    loci: Dict[str, Dict[str, object]] = {}
    for lab in _LOCUS_ORDER:
        if lab not in pools:
            continue
        pts = pools[lab]
        kept: List[Tuple[float, float, int, bool]] = []
        for p in pts:
            dup = False
            for q in kept:
                if float(np.hypot(p[0] - q[0], p[1] - q[1])) <= _DEDUP_TOL:
                    dup = True
                    break
            if not dup:
                kept.append(p)
        n_dupes = len(pts) - len(kept)
        mu = np.array([p[0] for p in kept], dtype=float)
        st = np.array([p[1] for p in kept], dtype=float)
        stab = np.array([p[2] for p in kept], dtype=int)
        bif = np.array([p[3] for p in kept], dtype=bool)
        order = _chain_order(mu, st)
        mu, st, stab, bif = mu[order], st[order], stab[order], bif[order]
        if mu.size <= 1 or link_cut is None:
            parts = [np.arange(mu.size, dtype=int)]
        else:
            d = np.hypot(np.diff(mu), np.diff(st))
            breaks = np.flatnonzero(d > link_cut)
            parts = []
            j = 0
            for b in breaks:
                parts.append(np.arange(j, int(b) + 1, dtype=int))
                j = int(b) + 1
            parts.append(np.arange(j, mu.size, dtype=int))
        loci[lab] = {
            "mu": mu,
            "state": st,
            "stability": stab,
            "bif": bif,
            "parts": parts,
            "n_dupes": n_dupes,
        }
    return loci, notes


def summarize_loci(
    loci: Dict[str, Dict[str, object]],
    notes: Sequence[str],
    norm: str = "l2",
) -> str:
    """Build the textual report of a locus-merged landscape.

    Parameters
    ----------
    loci : dict
        Output of :func:`merge_loci`.
    notes : sequence of str
        Diagnostic notes from :func:`merge_loci`.
    norm : str, optional
        State-norm convention the states were computed in; reported with
        the state ranges so no number is printed without its convention.

    Returns
    -------
    str
        Multi-line human readable summary.
    """
    lines: List[str] = []
    lines.append("Locus-merged landscape ({} locus/loci)".format(len(loci)))
    for note in notes:
        lines.append("  " + note)
    lines.append("")
    for lab in _LOCUS_ORDER:
        if lab not in loci:
            continue
        entry = loci[lab]
        mu = entry["mu"]
        st = entry["state"]
        stab = entry["stability"]
        lines.append(
            "locus {}: {} point(s), {} polyline part(s), "
            "{} duplicate(s) dropped".format(
                lab, int(mu.size), len(entry["parts"]), entry["n_dupes"]
            )
        )
        if mu.size:
            lines.append(
                "  mu range    : [{:.6g}, {:.6g}]".format(
                    float(mu.min()), float(mu.max())
                )
            )
            lines.append(
                "  state range : [{:.6g}, {:.6g}]   (state = {})".format(
                    float(st.min()), float(st.max()), norm_label(norm)
                )
            )
            flags = sorted(set(int(s) for s in stab))
            lines.append(
                "  stability   : {}".format(
                    ", ".join(
                        "{} x {}".format(
                            int(np.count_nonzero(stab == f)),
                            _STABILITY_NAME.get(f, "flag {}".format(f)),
                        )
                        for f in flags
                    )
                )
            )
            lines.append(
                "  singular pts: {}".format(
                    int(np.count_nonzero(entry["bif"]))
                )
            )
    return "\n".join(lines)


def plot_loci(
    loci: Dict[str, Dict[str, object]],
    output: str,
    title: Optional[str] = None,
    dpi: int = 150,
    norm: str = "l2",
) -> str:
    """Render the locus-merged landscape in the thesis's line style.

    Solid where stable, dashed where unstable (dotted where unknown), one
    colour per locus, cyan dots on the marked singular points; the
    unclassified remainder is drawn dashed grey and labelled.

    Parameters
    ----------
    loci : dict
        Output of :func:`merge_loci`.
    output : str
        Path of the image file to write.
    title : str, optional
        Figure title.
    dpi : int, optional
        Output resolution, by default 150.
    norm : str, optional
        State-norm convention of ``loci``; sets the y-axis label so the
        axis always names the convention it is drawn in.

    Returns
    -------
    str
        The ``output`` path, for convenience.

    Raises
    ------
    LandscapeError
        If the figure cannot be written.
    """
    fig, ax = plt.subplots(figsize=(8.0, 6.0))
    for lab in _LOCUS_ORDER:
        if lab not in loci:
            continue
        entry = loci[lab]
        mu = entry["mu"]
        st = entry["state"]
        stab = entry["stability"]
        color = _LOCUS_COLOR[lab]
        labelled = False
        legend_name = (
            "other (unclassified)" if lab == "other" else "locus {}".format(lab)
        )
        for part in entry["parts"]:
            if part.size == 0:
                continue
            if part.size == 1:
                ax.plot(
                    mu[part],
                    st[part],
                    marker=".",
                    markersize=5,
                    linestyle="none",
                    color=color,
                    zorder=3,
                    label=None if labelled else legend_name,
                )
                labelled = True
                continue
            if lab == "other":
                ax.plot(
                    mu[part],
                    st[part],
                    color=color,
                    linestyle="--",
                    linewidth=1.4,
                    zorder=3,
                    label=None if labelled else legend_name,
                )
                labelled = True
                continue
            for start, stop, flag in _stability_runs(stab[part]):
                lo = max(start - 1, 0)
                idx = part[lo:stop]
                ax.plot(
                    mu[idx],
                    st[idx],
                    color=color,
                    linestyle=_STABILITY_STYLE.get(flag, ":"),
                    linewidth=1.8,
                    zorder=3,
                    label=None if labelled else legend_name,
                )
                labelled = True
    bif_mu = np.concatenate(
        [loci[lab]["mu"][loci[lab]["bif"]] for lab in loci] + [np.empty(0)]
    )
    bif_st = np.concatenate(
        [loci[lab]["state"][loci[lab]["bif"]] for lab in loci] + [np.empty(0)]
    )
    if bif_mu.size:
        ax.plot(
            bif_mu,
            bif_st,
            marker="o",
            markersize=6,
            linestyle="none",
            color="cyan",
            markeredgecolor="teal",
            zorder=6,
            label="singular point",
        )
    for flag in (1, -1):
        if any(
            np.any(loci[lab]["stability"] == flag)
            for lab in loci
            if lab != "other"
        ):
            ax.plot(
                [],
                [],
                color="0.25",
                linestyle=_STABILITY_STYLE[flag],
                linewidth=1.6,
                label=_STABILITY_NAME[flag],
            )
    ax.set_xlabel(r"$\mu$")
    ax.set_ylabel(_norm_axis_label(norm))
    ax.set_title(title if title else "Locus-merged landscape")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=8, framealpha=0.9)
    fig.tight_layout()
    try:
        fig.savefig(output, dpi=dpi)
    except (OSError, ValueError) as exc:
        plt.close(fig)
        raise LandscapeError("cannot write '{}': {}".format(output, exc))
    plt.close(fig)
    return output


def locus_a_deviation(
    mu: np.ndarray, state: np.ndarray
) -> Dict[str, float]:
    """Deviation of locus-A points from the analytic constant branch.

    On a spatially constant state ``psi = c`` the ModifiedBratu equation
    reduces to ``mu = c * exp(-c)`` exactly.  The state coordinate equals
    ``c`` in the ``l2``, ``h1`` and ``discrete`` conventions alike (a
    constant is the one state on which they agree), so this comparison is
    unaffected by the norm choice; in the ``thesis`` convention it carries
    the factor ``sqrt(n/(n-1))`` and the caller is warned.  Both the
    vertical residual ``|mu - c e^{-c}|`` and the orthogonal distance to
    the curve ``t -> (t e^{-t}, t)`` are reported; the latter is computed
    on a dense parameter grid with a parabolic refinement of the minimum.

    Parameters
    ----------
    mu, state : numpy.ndarray
        Locus-A point coordinates.  Shape ``(n,)`` each.

    Returns
    -------
    dict of str to float
        ``n``, ``res_max``, ``res_mean``, ``orth_max``, ``orth_mean``.
    """
    mu = np.asarray(mu, dtype=float)
    c = np.asarray(state, dtype=float)
    res = np.abs(mu - c * np.exp(-c))
    lo = min(float(c.min()), 0.0) - 1.0
    hi = float(c.max()) + 1.0
    tgrid = np.linspace(lo, hi, 20001)
    gx = tgrid * np.exp(-tgrid)
    orth = np.empty(c.size, dtype=float)
    for i in range(c.size):
        d2 = (gx - mu[i]) ** 2 + (tgrid - c[i]) ** 2
        j = int(np.argmin(d2))
        if 0 < j < tgrid.size - 1:
            # Parabolic vertex of d^2 over the three bracketing samples.
            y0, y1, y2 = d2[j - 1], d2[j], d2[j + 1]
            denom = y0 - 2.0 * y1 + y2
            shift = 0.5 * (y0 - y2) / denom if denom > 0.0 else 0.0
            t = tgrid[j] + shift * (tgrid[1] - tgrid[0])
            orth[i] = float(np.hypot(t * np.exp(-t) - mu[i], t - c[i]))
        else:
            orth[i] = float(np.sqrt(d2[j]))
    return {
        "n": int(c.size),
        "res_max": float(res.max()) if res.size else 0.0,
        "res_mean": float(res.mean()) if res.size else 0.0,
        "orth_max": float(orth.max()) if orth.size else 0.0,
        "orth_mean": float(orth.mean()) if orth.size else 0.0,
    }


def _cluster_1d(
    indices: np.ndarray, weights: np.ndarray, gap: int
) -> List[float]:
    """Group sorted pixel indices into clusters and return their centroids.

    Parameters
    ----------
    indices : numpy.ndarray
        Sorted pixel indices with nonzero weight.
    weights : numpy.ndarray
        Weight per pixel index (full-axis array, indexed by ``indices``).
    gap : int
        Two indices farther apart than this start a new cluster.

    Returns
    -------
    list of float
        Weighted centroid per cluster, in ascending order.
    """
    groups: List[List[int]] = []
    for idx in indices:
        if groups and int(idx) - groups[-1][-1] <= gap:
            groups[-1].append(int(idx))
        else:
            groups.append([int(idx)])
    return [
        float(np.average(gr, weights=weights[gr])) for gr in groups
    ]


def digitize_thesis_fig(
    png_path: str,
) -> Tuple[Dict[str, np.ndarray], List[str]]:
    """Digitize curves B and C of Wouters (2019) fig. 9.3 from a render.

    Expects the 400 dpi ``pdftoppm`` render of physical page 284 of the
    thesis PDF (any render containing fig. 9.3 as the topmost axes frame
    works).  The plot frame is found as the first pair of horizontal rules
    joined by two vertical rules; the pixel-to-data calibration is fitted
    to the tick-label centroids (values :data:`_FIG93_XTICK_VALUES` /
    :data:`_FIG93_YTICK_VALUES`, read from the figure itself); curve
    pixels are selected by colour (B yellow-green, C pink/violet) and
    reduced to one point per pixel column.  Locus A is not digitized: it
    is double-valued in ``mu`` (it folds) and is compared analytically by
    :func:`locus_a_deviation` instead.

    The two polylines are written as ``<png basename>_locusB.csv`` and
    ``..._locusC.csv`` next to the PNG so the metric is reproducible.

    Parameters
    ----------
    png_path : str
        Path to the page render.

    Returns
    -------
    tuple of (dict, list of str)
        ``label -> (n, 2)`` array of ``(mu, psi)`` sorted by ``mu``, and
        diagnostic note lines (calibration figures, output files).

    Raises
    ------
    LandscapeError
        If the image cannot be read, the frame or the labels cannot be
        located, or the calibration residuals are out of tolerance.
    """
    try:
        img = plt.imread(png_path)
    except (OSError, ValueError) as exc:
        raise LandscapeError("cannot read '{}': {}".format(png_path, exc))
    if img.ndim != 3 or img.shape[2] < 3:
        raise LandscapeError(
            "'{}' is not an RGB image (shape {})".format(png_path, img.shape)
        )
    rgb = np.asarray(img[:, :, :3], dtype=float)
    if rgb.max() > 1.5:
        rgb = rgb / 255.0
    gray = rgb.sum(axis=2)
    dark = gray < 0.6
    height, width = gray.shape

    rowcount = dark.sum(axis=1)
    hlines = _cluster_1d(
        np.flatnonzero(rowcount > 0.35 * width), rowcount, gap=3
    )
    frame = None
    for i in range(len(hlines)):
        for j in range(i + 1, len(hlines)):
            top, bot = int(round(hlines[i])), int(round(hlines[j]))
            if bot - top < 100:
                continue
            colcount = dark[top:bot, :].sum(axis=0)
            vcols = _cluster_1d(
                np.flatnonzero(colcount > 0.95 * (bot - top)), colcount, gap=3
            )
            if len(vcols) >= 2:
                frame = (top, bot, int(round(vcols[0])), int(round(vcols[-1])))
                break
        if frame is not None:
            break
    if frame is None:
        raise LandscapeError(
            "cannot locate the fig. 9.3 plot frame in '{}'".format(png_path)
        )
    top, bot, left, right = frame
    notes: List[str] = [
        "thesis frame: rows {}..{}, cols {}..{}".format(top, bot, left, right)
    ]

    # x calibration from the tick-label centroids below the frame.  The
    # two edge labels can be clipped/shifted, so the fit uses the interior
    # labels only and the edges are merely detected.
    band = (gray[bot + 12 : bot + 80, :] < 1.2).sum(axis=0)
    xcents = _cluster_1d(np.flatnonzero(band > 0), band, gap=15)
    if len(xcents) != len(_FIG93_XTICK_VALUES) + 1:
        raise LandscapeError(
            "expected {} x-label clusters below the frame, found {}".format(
                len(_FIG93_XTICK_VALUES) + 1, len(xcents)
            )
        )
    del xcents[_FIG93_XTICK_MU_TEXT_INDEX]
    xs = np.asarray(xcents[1:-1], dtype=float)
    xvals = np.asarray(_FIG93_XTICK_VALUES[1:-1], dtype=float)
    xcoef = np.polyfit(xs, xvals, 1)
    xresid = xvals - np.polyval(xcoef, xs)
    if float(np.abs(xresid).max()) > 0.005:
        raise LandscapeError(
            "x calibration residual {:.3g} exceeds 0.005 in mu".format(
                float(np.abs(xresid).max())
            )
        )
    notes.append(
        "x calibration: {:.4g} px per 0.05 in mu, max label residual "
        "{:.2g} mu".format(
            0.05 / xcoef[0], float(np.abs(xresid).max())
        )
    )

    # y calibration: label digits left of the frame; the axis title sits
    # farther left and can bleed clusters into the band, so the two
    # extreme clusters anchor a provisional map and only clusters within
    # 12 px of a predicted tick row are kept.
    r0 = max(top - 30, 0)
    band = (gray[r0 : bot + 30, max(left - 56, 0) : left - 6] < 1.2).sum(axis=1)
    ycents = [
        c + r0 for c in _cluster_1d(np.flatnonzero(band > 0), band, gap=15)
    ]
    yvals_all = np.asarray(_FIG93_YTICK_VALUES, dtype=float)
    if len(ycents) < len(yvals_all):
        raise LandscapeError(
            "expected at least {} y-label clusters, found {}".format(
                len(yvals_all), len(ycents)
            )
        )
    row_first, row_last = ycents[0], ycents[-1]
    pred = row_first + (yvals_all - yvals_all[0]) * (
        (row_last - row_first) / (yvals_all[-1] - yvals_all[0])
    )
    inl_rows: List[float] = []
    inl_vals: List[float] = []
    cents_arr = np.asarray(ycents, dtype=float)
    for pr, v in zip(pred, yvals_all):
        k = int(np.argmin(np.abs(cents_arr - pr)))
        if abs(cents_arr[k] - pr) <= 12.0:
            inl_rows.append(float(cents_arr[k]))
            inl_vals.append(float(v))
    if len(inl_rows) < 5:
        raise LandscapeError(
            "only {} y-label clusters match the tick grid".format(
                len(inl_rows)
            )
        )
    ycoef = np.polyfit(inl_rows, inl_vals, 1)
    yresid = np.asarray(inl_vals) - np.polyval(ycoef, inl_rows)
    if float(np.abs(yresid).max()) > 0.06:
        raise LandscapeError(
            "y calibration residual {:.3g} exceeds 0.06 in ||psi||".format(
                float(np.abs(yresid).max())
            )
        )
    notes.append(
        "y calibration: {:.4g} px per 1.0 in ||psi|| ({} label(s), max "
        "residual {:.2g})".format(
            1.0 / abs(ycoef[0]), len(inl_rows), float(np.abs(yresid).max())
        )
    )

    sub = rgb[top + 3 : bot - 2, left + 3 : right - 2]
    red, grn, blu = sub[:, :, 0], sub[:, :, 1], sub[:, :, 2]
    masks = {
        "B": (grn > 0.7) & (blu < 0.65) & (grn > blu + 0.2) & (red > 0.4),
        "C": (red > 0.6) & (blu > 0.6) & (grn < np.minimum(red, blu) - 0.08),
    }
    refs: Dict[str, np.ndarray] = {}
    base, _ = os.path.splitext(png_path)
    rows_idx = np.arange(sub.shape[0], dtype=float)
    for lab in ("B", "C"):
        mask = masks[lab]
        cols = np.flatnonzero(mask.any(axis=0))
        if cols.size == 0:
            raise LandscapeError(
                "no pixels matched the colour of thesis locus {} in "
                "'{}'".format(lab, png_path)
            )
        pts = []
        for col in cols:
            wcol = mask[:, col].astype(float)
            rmean = float(np.average(rows_idx, weights=wcol))
            pts.append(
                (
                    float(np.polyval(xcoef, col + left + 3)),
                    float(np.polyval(ycoef, rmean + top + 3)),
                )
            )
        arr = np.asarray(pts, dtype=float)
        arr = arr[np.argsort(arr[:, 0])]
        refs[lab] = arr
        out_csv = "{}_locus{}.csv".format(base, lab)
        try:
            with open(out_csv, "w") as handle:
                handle.write("mu,psi\n")
                for p in arr:
                    handle.write("{:.6g},{:.6g}\n".format(p[0], p[1]))
        except OSError as exc:
            raise LandscapeError(
                "cannot write '{}': {}".format(out_csv, exc)
            )
        notes.append(
            "thesis locus {}: {} column sample(s), mu in "
            "[{:.4g}, {:.4g}], wrote {}".format(
                lab, arr.shape[0], float(arr[0, 0]), float(arr[-1, 0]),
                out_csv,
            )
        )
    return refs, notes


def thesis_metrics(
    loci: Dict[str, Dict[str, object]],
    thesis_fig: Optional[str] = None,
    norm: str = "l2",
) -> str:
    """Quantitative per-locus deviation table (``--thesis-metrics``).

    Parameters
    ----------
    loci : dict
        Output of :func:`merge_loci`.
    thesis_fig : str, optional
        Path to the thesis figure render; enables the digitized B/C
        comparison via :func:`digitize_thesis_fig`.
    norm : str, optional
        State-norm convention of ``loci``; reported, and used to phrase
        the residual convention mismatch against the thesis's own y-axis.

    Returns
    -------
    str
        Multi-line metrics report.

    Raises
    ------
    LandscapeError
        Propagated from :func:`digitize_thesis_fig`.
    """
    lines: List[str] = []
    lines.append("Thesis-comparison metrics")
    lines.append("  state norm: {} [{}]".format(norm, norm_label(norm)))
    lines.append(
        "  NOTE: the B/C reference is a read-off of the thesis figure "
        "render; its precision is about +/-0.1 in ||psi||.  Deviations at "
        "or below that scale bottom out at the thesis's own resolution."
    )
    lines.append(
        "  NOTE: the thesis's own y-axis is ITS discrete norm "
        "(1/(n-1)) sum psi_i^2 on n = {} FD points, not a continuous "
        "L2 norm; run --norm thesis to measure that convention gap on "
        "these states (it does NOT vanish).".format(_THESIS_FD_POINTS)
    )
    if "A" in loci:
        dev = locus_a_deviation(loci["A"]["mu"], loci["A"]["state"])
        if norm == "thesis":
            lines.append(
                "  WARNING: in the 'thesis' convention the constant state c "
                "reads c*sqrt(n/(n-1)) = 1.00504*c, so the locus-A residual "
                "below is dominated by that convention factor, not by the "
                "solution."
            )
        lines.append(
            "locus A vs analytic mu = c*exp(-c), c = the state coordinate "
            "({} point(s); {}):".format(
                dev["n"], _const_exactness_note(norm)
            )
        )
        lines.append(
            "  |mu - c e^-c|       : max {:.3g}, mean {:.3g}".format(
                dev["res_max"], dev["res_mean"]
            )
        )
        lines.append(
            "  orthogonal distance : max {:.3g}, mean {:.3g}".format(
                dev["orth_max"], dev["orth_mean"]
            )
        )
    else:
        lines.append("locus A: absent from this landscape")
    if thesis_fig is None:
        lines.append(
            "locus B/C: no --thesis-fig given; digitized comparison skipped"
        )
        return "\n".join(lines)
    if norm == "h1":
        lines.append(
            "  WARNING: the digitized reference is an L2-type norm of the "
            "thesis's solution, so the B/C rows below are NOT a thesis "
            "comparison in the H1 convention -- they carry the extra "
            "int psi'^2 term (measured: 1.4 to 5.6 on these landscapes).  "
            "Use --norm l2 (or --norm thesis) to compare against the figure."
        )
    refs, dnotes = digitize_thesis_fig(thesis_fig)
    for note in dnotes:
        lines.append("  " + note)
    for lab in ("B", "C"):
        if lab not in loci:
            lines.append(
                "locus {}: absent from this landscape (nothing to "
                "compare)".format(lab)
            )
            continue
        ref = refs[lab]
        mu = loci[lab]["mu"]
        st = loci[lab]["state"]
        inr = (mu >= ref[0, 0]) & (mu <= ref[-1, 0])
        if not np.any(inr):
            lines.append(
                "locus {}: no points inside the digitized mu range "
                "[{:.4g}, {:.4g}]".format(
                    lab, float(ref[0, 0]), float(ref[-1, 0])
                )
            )
            continue
        interp = np.interp(mu[inr], ref[:, 0], ref[:, 1])
        devs = np.abs(st[inr] - interp)
        lines.append(
            "locus {} vs digitized thesis polyline ({} of {} point(s) "
            "inside its mu range):".format(
                lab, int(np.count_nonzero(inr)), int(mu.size)
            )
        )
        lines.append(
            "  |state - psi_fig(mu)| : max {:.3g}, mean {:.3g}   "
            "{}".format(
                float(devs.max()),
                float(devs.mean()),
                _ref_precision_note(norm),
            )
        )
    return "\n".join(lines)


_SELFTEST_CSV = """curve,point,L,normU,stability,negatives,isBifurcation,parentCurve,parentPointIdx,equilibrium
0,0,0.0,0.0,1,0,0,-1,-1,1
0,1,0.5,0.4,1,0,0,-1,-1,1
0,2,0.9,0.9,1,-1,1,-1,-1,1
0,3,1.1,1.6,1,0,0,-1,-1,1
0,4,1.0,2.3,-1,1,0,-1,-1,1
0,5,0.7,2.9,-1,1,0,-1,-1,1
1,0,0.88,1.0,-1,1,0,0,2,1
1,1,0.80,1.3,-1,1,0,0,2,1
1,2,0.70,1.7,-1,1,0,0,2,1
2,0,0.30,3.5,0,-1,0,-1,-1,0
"""


def _make_constant_curve(
    cid: int, cs: np.ndarray, n_dof: int = 34
) -> Tuple[Curve, np.ndarray]:
    """Build a synthetic constant-state curve on the analytic branch.

    Parameters
    ----------
    cid : int
        Curve identifier.
    cs : numpy.ndarray
        Constant state values ``c`` per point; ``mu = c * exp(-c)``.
    n_dof : int, optional
        Number of coefficients per solution vector.

    Returns
    -------
    tuple of (Curve, numpy.ndarray)
        The curve and its ``(n_dof, n)`` solution-vector matrix.
    """
    cs = np.asarray(cs, dtype=float)
    n = int(cs.size)
    curve = Curve(
        cid=cid,
        lam=cs * np.exp(-cs),
        norm_u=cs * np.sqrt(float(n_dof)),
        stability=np.ones(n, dtype=int),
        bifurcation=np.zeros(n, dtype=bool),
        parent_curve=-1,
        parent_point=-1,
    )
    return curve, np.tile(cs, (n_dof, 1))


def _selftest_norms() -> None:
    """Self-check the spline basis, the Gram matrices and the state norms.

    The Gram matrix is checked against *identities*, not against an
    eyeballed value: every row of the mass matrix must integrate its own
    basis function, ``sum_j M[i, j] = (t_{i+p+1} - t_i) / (p+1)``, the
    total mass must be exactly 1 (partition of unity on a unit domain),
    ``M`` must be symmetric positive definite and banded with half
    bandwidth ``p``.  Aggregate checks alone would pass with a mis-scaled
    row.

    Raises
    ------
    AssertionError
        If any check fails.
    """
    p, n_dof = _BASIS_DEGREE, 34
    mass, stiff, knots = spline_gram_matrices(n_dof, p)
    n_elem = n_dof - p
    assert knots.size == n_dof + p + 1
    # Per-row integral identity -- catches a mis-scaled or shifted row.
    rows = mass.sum(axis=1)
    exact = np.array(
        [(knots[i + p + 1] - knots[i]) / (p + 1.0) for i in range(n_dof)]
    )
    assert float(np.abs(rows - exact).max()) < 1e-14, (
        "mass-matrix row sums must equal (t_{i+p+1} - t_i)/(p+1)"
    )
    assert abs(float(mass.sum()) - 1.0) < 1e-14, (
        "the total mass must be exactly 1 on the unit domain"
    )
    assert float(np.abs(mass - mass.T).max()) < 1e-15, "M must be symmetric"
    assert float(np.linalg.eigvalsh(mass).min()) > 0.0, "M must be SPD"
    band = np.abs(np.subtract.outer(np.arange(n_dof), np.arange(n_dof)))
    assert float(np.abs(mass[band > p]).max()) == 0.0, (
        "M must be banded with half bandwidth p"
    )
    assert float(np.abs(stiff.sum(axis=1)).max()) < 1e-12, (
        "the stiffness matrix must annihilate constants"
    )

    # Constant state: l2 == h1 == discrete == c exactly; the thesis
    # convention carries the documented sqrt(n/(n-1)) factor.
    c = 3.7
    const = np.full(n_dof, c)
    for kind in ("l2", "h1", "discrete"):
        assert abs(float(state_norm(const, kind)[0]) - c) < 1e-12, (
            "norm '{}' must be exact on a constant state".format(kind)
        )
    factor = np.sqrt(
        _THESIS_FD_POINTS / float(_THESIS_FD_POINTS - 1)
    )
    assert abs(float(state_norm(const, "thesis")[0]) - c * factor) < 1e-12, (
        "the thesis convention must give c*sqrt(n/(n-1)) on a constant"
    )

    # A hand-computable non-constant state: psi(x) = x lies in the space
    # (degree >= 1), its coefficients are the Greville abscissae, and
    # int_0^1 x^2 = 1/3, int_0^1 (x^2 + 1) = 4/3.
    greville = np.array(
        [float(knots[i + 1 : i + p + 1].mean()) for i in range(n_dof)]
    )
    ramp_l2 = float(state_norm(greville, "l2")[0])
    ramp_h1 = float(state_norm(greville, "h1")[0])
    ramp_disc = float(state_norm(greville, "discrete")[0])
    assert abs(ramp_l2 - 1.0 / np.sqrt(3.0)) < 1e-12, (
        "the L2 norm of psi(x) = x must be exactly sqrt(1/3)"
    )
    assert abs(ramp_h1 - 2.0 / np.sqrt(3.0)) < 1e-12, (
        "the H1 norm of psi(x) = x must be exactly sqrt(4/3)"
    )
    # ... and the historical norm is measurably wrong on it.
    assert abs(ramp_disc - 1.0 / np.sqrt(3.0)) > 1e-3, (
        "the discrete norm is expected to be biased off constant states"
    )

    # The pitchfork pair: psi = cbar +/- a cos(2 pi x) are two DIFFERENT
    # functions with the SAME continuous L2 norm (the cross term integrates
    # to zero).  The discrete norm splits them, because the 32 interior
    # Greville cosines cancel while the two endpoints contribute 1 each,
    # leaving sum cos = 2 and a squared gap of 4*cbar*a*2/N_dof.
    cbar, amp = 6.0, 1.556
    wave = np.cos(2.0 * np.pi * greville)
    well, bump = cbar + amp * wave, cbar - amp * wave
    l2w = float(state_norm(well, "l2")[0])
    l2b = float(state_norm(bump, "l2")[0])
    assert abs(l2w - l2b) / l2w < 1e-4, (
        "the two pitchfork sides must agree in the continuous L2 norm"
    )
    dw = float(state_norm(well, "discrete")[0])
    db = float(state_norm(bump, "discrete")[0])
    assert dw - db > 0.17, "the discrete norm must show the known artefact"
    predicted = 4.0 * cbar * amp * 2.0 / n_dof
    assert abs((dw * dw - db * db) - predicted) < 1e-9, (
        "the squared discrete gap must equal 4*cbar*a*sum(cos)/N_dof"
    )
    assert abs(float(wave.sum()) - 2.0) < 1e-12, (
        "the Greville cosines must sum to 2 (only the endpoints survive)"
    )

    # Matrix input: one norm per column.
    both = np.stack([well, bump], axis=1)
    vals = state_norm(both, "l2")
    assert vals.shape == (2,)
    assert abs(vals[0] - l2w) < 1e-14 and abs(vals[1] - l2b) < 1e-14

    # Clean refusals.
    try:
        state_norm(const, "nope")
    except LandscapeError as exc:
        assert "unknown state norm" in str(exc)
    else:  # pragma: no cover - defensive
        raise AssertionError("an unknown norm was not rejected")
    try:
        spline_gram_matrices(2, 5)
    except LandscapeError as exc:
        assert "too small" in str(exc)
    else:  # pragma: no cover - defensive
        raise AssertionError("an impossible (N_dof, degree) was not rejected")
    assert n_elem == 32, "the reconstructed element count must be N_dof - p"

    # --- the assumed degree is checked against the knot structure -----
    # N_dof = 2**numHref + degree, so N_dof - degree must be a power of
    # two.  On the shipped 34-dof export only 2, 18, 26, 30, 32, 33 are
    # reachable; the default 2 is the smallest and needs no note.
    assert check_basis_degree(34, 2) == [], (
        "the driver's own (N_dof, degree) must pass without a note"
    )
    assert _admissible_degrees(34) == [2, 18, 26, 30, 32, 33], (
        "admissible degrees must be exactly those with 2^k elements"
    )
    for bad in (3, 4, 5, 17):  # 31 / 30 / 29 / 17 elements
        try:
            check_basis_degree(34, bad)
        except LandscapeError as exc:
            assert "power of two" in str(exc), (
                "a non-refinable element count must say so"
            )
            assert "1.24" in str(exc), (
                "the refusal must name the measured consequence"
            )
        else:  # pragma: no cover - defensive
            raise AssertionError(
                "degree {} was accepted on 34 dofs".format(bad)
            )
    try:
        check_basis_degree(34, 3)
    except LandscapeError as exc:
        assert "not a 1D export" not in str(exc), (
            "a 34-dof export must not be blamed on dimensionality"
        )
    # The 2D landscape (256 free dofs) has no credible 1D reconstruction:
    # no admissible degree leaves more elements than its own order.
    try:
        check_basis_degree(256, _BASIS_DEGREE)
    except LandscapeError as exc:
        assert (
            "not a 1D export" in str(exc) and "--norm discrete" in str(exc)
        ), "a 256-dof export must be diagnosed as not 1D, with the way out"
    else:  # pragma: no cover - defensive
        raise AssertionError("a 2D export was read as a 1D spline")
    # A landscape run with -e 2 exports 35 dofs: the DEFAULT degree is
    # then the wrong one, and must be refused with the right suggestion.
    try:
        check_basis_degree(35, _BASIS_DEGREE)
    except LandscapeError as exc:
        assert "degree(s) 3," in str(exc), (
            "the refusal must name the degrees the export can have"
        )
    else:  # pragma: no cover - defensive
        raise AssertionError("degree 2 was accepted on a 35-dof export")
    # Admissible but not the driver's configuration: a note, not a refusal
    # (the export cannot decide between them).
    for alt in (18, 33):
        alt_notes = check_basis_degree(34, alt)
        assert len(alt_notes) == 1 and alt_notes[0].startswith("NOTE:"), (
            "a non-default admissible degree must be reported"
        )
    assert "numHref = 0" in check_basis_degree(34, 33)[0], (
        "the note must state the implied refinement count"
    )
    for degree, n_dofs, expect in ((0, 34, ">= 1"), (5, 2, "too small")):
        try:
            check_basis_degree(n_dofs, degree)
        except LandscapeError as exc:
            assert expect in str(exc)
        else:  # pragma: no cover - defensive
            raise AssertionError("an impossible degree was accepted")


def _selftest_loci() -> None:
    """Self-check the locus classifier, merge, window clip and jump cut.

    Raises
    ------
    AssertionError
        If any check fails.
    """
    # --- classifier on synthetic constant / mode-1 / mode-2 vectors ---
    x = np.linspace(0.0, 1.0, 34)
    const = np.full(34, 2.16)
    mode1 = 2.16 + 1.3 * np.cos(np.pi * x)
    mode2 = 2.16 + 1.5 * np.cos(2.0 * np.pi * x)
    mode3 = 2.16 + 1.5 * np.cos(3.0 * np.pi * x)
    scale = 4.0
    assert classify_locus(const, scale) == "A"
    # Round-off-level noise must not break the constant classification.
    assert classify_locus(const + 1e-9 * np.sin(7.0 * x), scale) == "A"
    assert classify_locus(mode1, scale) == "B"
    assert classify_locus(mode2, scale) == "C"
    assert classify_locus(mode3, scale) == "other"
    assert count_sign_changes(mode1 + 1e-8 * np.sin(13.0 * x)) == 1
    assert count_sign_changes(mode2 + 1e-8 * np.sin(13.0 * x)) == 2
    assert count_sign_changes(const) == 0, "constant state must count 0"

    # --- merge: two overlapping A-segments -> ONE ordered polyline ---
    c1, u1 = _make_constant_curve(0, np.linspace(0.10, 1.00, 10))
    c2, u2 = _make_constant_curve(1, np.linspace(0.95, 2.00, 12))
    loci, _ = merge_loci([c1, c2], {0: u1, 1: u2})
    assert set(loci) == {"A"}, "two constant segments must merge into A only"
    entry = loci["A"]
    assert len(entry["parts"]) == 1, "overlapping segments must give 1 part"
    assert entry["mu"].size == 22
    dstate = np.diff(entry["state"])
    assert np.all(dstate > 0.0) or np.all(dstate < 0.0), (
        "chained polyline must be ordered along the branch"
    )

    # --- duplicate dropping within the state-space tolerance ---
    c3, u3 = _make_constant_curve(2, np.array([1.0, 1.5]))
    loci_dup, _ = merge_loci([c1, c2, c3], {0: u1, 1: u2, 2: u3})
    assert loci_dup["A"]["n_dupes"] == 1, "the repeated c=1 point must drop"
    assert loci_dup["A"]["mu"].size == 23

    # --- window clip removes everything below MU_MIN before merging ---
    loci_win, _ = merge_loci([c1, c2], {0: u1, 1: u2}, window=0.2)
    assert float(loci_win["A"]["mu"].min()) >= 0.2
    n_expect = int(
        np.count_nonzero(np.concatenate([c1.lam, c2.lam]) >= 0.2)
    )
    assert loci_win["A"]["mu"].size == n_expect

    # --- solver-jump removal splits the polyline at the outlier link ---
    c4, u4 = _make_constant_curve(
        0, np.concatenate([np.linspace(0.10, 0.40, 16), [5.0]])
    )
    loci_jump, _ = merge_loci([c4], {0: u4})
    assert len(loci_jump["A"]["parts"]) == 2, (
        "a solver jump must split the merged polyline"
    )
    c5, u5 = _make_constant_curve(0, np.linspace(0.10, 0.40, 16))
    loci_nojump, _ = merge_loci([c5], {0: u5})
    assert len(loci_nojump["A"]["parts"]) == 1, (
        "without a jump the polyline must stay in one part"
    )
    cut, _note = derive_jump_cut(
        np.hypot(np.diff(c5.lam), np.diff(np.linspace(0.10, 0.40, 16)))
    )
    assert cut is None, "a smooth step distribution must yield no cut"

    # --- merge_loci honours --norm ------------------------------------
    # On the constant fixtures every convention but 'thesis' agrees, and
    # 'thesis' differs by exactly the documented sqrt(n/(n-1)) factor.
    loci_l2, _ = merge_loci([c1, c2], {0: u1, 1: u2}, norm="l2")
    loci_disc, _ = merge_loci([c1, c2], {0: u1, 1: u2}, norm="discrete")
    assert float(
        np.abs(loci_l2["A"]["state"] - loci_disc["A"]["state"]).max()
    ) < 1e-12, "on constant states l2 and discrete must coincide"
    loci_th, _ = merge_loci([c1, c2], {0: u1, 1: u2}, norm="thesis")
    ratio = loci_th["A"]["state"] / loci_l2["A"]["state"]
    assert float(np.abs(ratio - np.sqrt(100.0 / 99.0)).max()) < 1e-9, (
        "the thesis convention must scale a constant state by sqrt(n/(n-1))"
    )
    try:
        merge_loci([c1, c2], {0: u1, 1: u2}, norm="nope")
    except LandscapeError:
        pass
    else:  # pragma: no cover - defensive
        raise AssertionError("merge_loci accepted an unknown norm")

    # --- the degree is validated wherever the basis is actually used ---
    try:
        merge_loci([c1, c2], {0: u1, 1: u2}, degree=3)
    except LandscapeError as exc:
        assert "power of two" in str(exc)
    else:  # pragma: no cover - defensive
        raise AssertionError("merge_loci accepted an unreachable degree")
    _loci_d3, notes_d3 = merge_loci(
        [c1, c2], {0: u1, 1: u2}, norm="discrete", degree=3
    )
    assert any("unused by the 'discrete' norm" in n for n in notes_d3), (
        "the 'discrete' norm must say that --basis-degree did nothing"
    )
    _loci_d2, notes_d2 = merge_loci([c1, c2], {0: u1, 1: u2})
    assert not any(n.startswith("NOTE:") for n in notes_d2), (
        "the driver's own degree must not raise a degree note"
    )

    with tempfile.TemporaryDirectory(prefix="plot_landscape_loci_") as tmp:
        # The locus plot must render.
        out = os.path.join(tmp, "loci.png")
        plot_loci(loci, out, title="selftest loci", dpi=80)
        assert os.path.getsize(out) > 0, "locus plot was not written"

        # The HDF5 read path (round trip) -- only when h5py is present.
        if _HAS_H5PY:
            h5_path = os.path.join(tmp, "landscape.h5")
            with h5py.File(h5_path, "w") as handle:
                handle.create_dataset("curve_0000_U", data=u1)
            vectors = read_solution_vectors(h5_path, [c1])
            assert vectors[0].shape == u1.shape
            assert float(np.abs(vectors[0] - u1).max()) == 0.0
            # A dataset/CSV length mismatch must be rejected cleanly.
            c_short, _u_short = _make_constant_curve(
                0, np.linspace(0.95, 2.00, 12)
            )
            try:
                read_solution_vectors(h5_path, [c_short])
            except LandscapeError as exc:
                assert "out of sync" in str(exc)
            else:  # pragma: no cover - defensive
                raise AssertionError("h5/CSV mismatch was not detected")

        # A missing .h5 must be refused with a clean message.
        try:
            read_solution_vectors(os.path.join(tmp, "nope.h5"), [c1])
        except LandscapeError:
            pass
        else:  # pragma: no cover - defensive
            raise AssertionError("missing .h5 was not refused")


def _selftest_reports() -> None:
    """Self-check the norm-dependent claims printed by ``--thesis-metrics``.

    Both strings checked here used to be printed unconditionally: the
    "exact on constant states" parenthetical (false, and self-contradicted
    by the WARNING above it, under ``--norm thesis``) and the "~0.1"
    read-off floor (meaningless on H1 rows).  They are exercised through
    :func:`thesis_metrics` itself -- the helpers alone would not notice a
    hard-coded string creeping back into the report.

    Raises
    ------
    AssertionError
        If any check fails.
    """
    # The phrase helpers, per convention.
    assert "NOT exact" in _const_exactness_note("thesis")
    for kind in ("l2", "h1", "discrete"):
        note = _const_exactness_note(kind)
        assert "NOT exact" not in note and "exact on constant" in note
    assert "~0.1" not in _ref_precision_note("h1")
    for kind in ("l2", "discrete", "thesis"):
        assert _ref_precision_note(kind) == "(reference precision ~0.1)"

    # A constant-branch locus A, plus a B locus inside a stub reference.
    cs = np.linspace(0.30, 1.20, 12)
    loci: Dict[str, Dict[str, object]] = {
        "A": {"mu": cs * np.exp(-cs), "state": cs},
        "B": {
            "mu": np.linspace(0.06, 0.18, 5),
            "state": np.full(5, 5.5),
        },
    }
    a_thesis = thesis_metrics(loci, norm="thesis")
    assert "WARNING:" in a_thesis, "the thesis warning must still fire"
    assert "every norm is exact" not in a_thesis, (
        "the locus-A header must not claim exactness for every norm"
    )
    assert "NOT exact on constant states" in a_thesis, (
        "under --norm thesis the header must say the norm is not exact"
    )
    a_l2 = thesis_metrics(loci, norm="l2")
    assert "the 'l2' norm is exact on constant states" in a_l2
    assert "NOT exact" not in a_l2

    # The digitized B/C rows, through a stub digitizer so no PNG (and no
    # figure render) is needed.
    ref = np.array([[0.05, 5.0], [0.20, 6.0]])

    def _stub_digitizer(path: str):
        """Stand-in for :func:`digitize_thesis_fig` (selftest only)."""
        return {"B": ref, "C": ref}, ["stub reference"]

    real = globals()["digitize_thesis_fig"]
    globals()["digitize_thesis_fig"] = _stub_digitizer
    try:
        h1_text = thesis_metrics(loci, thesis_fig="stub.png", norm="h1")
        l2_text = thesis_metrics(loci, thesis_fig="stub.png", norm="l2")
    finally:
        globals()["digitize_thesis_fig"] = real
    assert "locus B vs digitized" in h1_text, (
        "the stubbed B row must actually be produced"
    )
    assert "(reference precision ~0.1)" not in h1_text, (
        "the read-off floor must not be printed on H1 rows"
    )
    assert "no read-off floor applies" in h1_text
    assert "locus B vs digitized" in l2_text
    assert "(reference precision ~0.1)" in l2_text, (
        "the read-off floor must survive on L2 rows"
    )


def _selftest_cli() -> None:
    """Self-check the command line guards that no value test can express.

    ``--norm l2`` and ``--basis-degree 2`` must be refused without
    ``--merge-loci`` exactly like every other value of those options, and
    an unreachable ``--basis-degree`` must come back as one clean error
    line, not a traceback.

    Raises
    ------
    AssertionError
        If any check fails.
    """

    def _run(argv: Sequence[str]) -> Tuple[int, str, str]:
        """Run :func:`main`, capturing the exit code and both streams."""
        out, err = io.StringIO(), io.StringIO()
        try:
            with contextlib.redirect_stdout(out):
                with contextlib.redirect_stderr(err):
                    code = main(list(argv))
        except SystemExit as exc:  # argparse's parser.error()
            code = int(exc.code) if exc.code is not None else 0
        return code, out.getvalue(), err.getvalue()

    for opt in (["--norm", "l2"], ["--norm", "discrete"]):
        code, _out, err = _run(["nowhere.csv"] + opt)
        assert code == 2 and "--norm requires --merge-loci" in err, (
            "'{}' must be refused without --merge-loci".format(" ".join(opt))
        )
    for val in ("2", "3"):
        code, _out, err = _run(["nowhere.csv", "--basis-degree", val])
        assert code == 2 and "--basis-degree requires --merge-loci" in err, (
            "--basis-degree {} must be refused without --merge-loci".format(
                val
            )
        )

    if not _HAS_H5PY:  # pragma: no cover - depends on the environment
        return
    with tempfile.TemporaryDirectory(prefix="plot_landscape_cli_") as tmp:
        csv_path = os.path.join(tmp, "landscape.csv")
        curve, u = _make_constant_curve(0, np.linspace(0.30, 1.20, 12))
        with open(csv_path, "w") as handle:
            handle.write(
                "curve,point,L,normU,stability,isBifurcation,parentCurve,"
                "parentPointIdx,equilibrium\n"
            )
            for j in range(curve.size):
                handle.write(
                    "0,{},{!r},{!r},1,0,-1,-1,1\n".format(
                        j, float(curve.lam[j]), float(curve.norm_u[j])
                    )
                )
        with h5py.File(os.path.join(tmp, "landscape.h5"), "w") as handle:
            handle.create_dataset("curve_0000_U", data=u)
        png = os.path.join(tmp, "loci.png")

        code, _out, err = _run(
            [csv_path, "-o", png, "--merge-loci", "--basis-degree", "3"]
        )
        assert code == 2, "an unreachable --basis-degree must fail"
        assert err.startswith("error: basis degree 3"), (
            "the refusal must be one clean 'error:' line, got {!r}".format(
                err
            )
        )
        assert len(err.strip().splitlines()) == 1, "no traceback, one line"

        code, out, _err = _run(
            [csv_path, "-o", png, "--merge-loci", "--basis-degree", "33"]
        )
        assert code == 0, "an admissible --basis-degree must be accepted"
        assert "NOTE: degree 33" in out, (
            "a non-default admissible degree must be noted, not silent"
        )
        code, out, _err = _run([csv_path, "-o", png, "--merge-loci"])
        assert code == 0 and "NOTE: degree" not in out, (
            "the default degree must stay quiet"
        )


def _selftest() -> int:
    """Build a synthetic landscape, render it and check the output.

    The synthetic file has two traced curves plus a one-point curve, a
    branch, a stability flip and an interior fold.  The locus-mode
    machinery is checked by :func:`_selftest_loci`.

    Returns
    -------
    int
        ``0`` on success.

    Raises
    ------
    AssertionError
        If any check on the parsed data or the rendered file fails.
    """
    with tempfile.TemporaryDirectory(prefix="plot_landscape_selftest_") as tmp:
        csv_path = os.path.join(tmp, "landscape.csv")
        with open(csv_path, "w") as handle:
            handle.write(_SELFTEST_CSV)

        curves = read_landscape(csv_path)
        assert len(curves) == 3, "expected 3 curves, got {}".format(len(curves))
        assert curves[0].size == 6
        assert curves[1].parent_curve == 0 and curves[1].parent_point == 2
        assert curves[2].size == 1, "single-point curve must survive parsing"
        assert curves[0].fold_is_interior(), "curve 0 must have an interior fold"
        assert len(curves[0].stability_transitions()) == 1
        assert int(np.count_nonzero(curves[0].bifurcation)) == 1
        assert int(curves[2].stability[0]) == 0

        # negatives column: parsed from the header-keyed 10-column CSV, not a
        # constant (a constant-column assertion would pass even if the parser
        # mapped the wrong field).
        assert list(int(v) for v in curves[0].negatives) == [0, 0, -1, 0, 1, 1]
        assert any(v == -1 for v in curves[0].negatives), "expected a sentinel value"
        assert any(v >= 0 for v in curves[0].negatives), "expected a recorded value"

        text = summarize(curves, reference=1.0)
        assert "curve 1 (child of curve 0 at point 2)" in text
        assert "deviation" in text
        assert "negatives" in text
        assert "1 x not recorded, 3 x 0, 2 x 1" in text, (
            "curve 0's negatives summary line did not match the expected counts"
        )

        # Legacy-format regression: a 9-column CSV using the OLD header must
        # still parse, with every negatives entry defaulting to -1 (proves the
        # OPTIONAL_COLUMNS default path, not just that parsing survives an
        # inserted column).
        legacy_csv = os.path.join(tmp, "legacy.csv")
        with open(legacy_csv, "w") as handle:
            handle.write(
                "curve,point,L,normU,stability,isBifurcation,parentCurve,"
                "parentPointIdx,equilibrium\n"
                "0,0,0.0,0.0,1,0,-1,-1,1\n"
                "0,1,0.5,0.4,1,0,-1,-1,1\n"
                "0,2,0.9,0.9,1,1,-1,-1,1\n"
            )
        legacy_curves = read_landscape(legacy_csv)
        assert len(legacy_curves) == 1
        assert all(int(v) == -1 for v in legacy_curves[0].negatives), (
            "legacy 9-column CSV must default every negatives entry to -1"
        )

        out = os.path.join(tmp, "landscape.png")
        plot_landscape(curves, out, title="selftest", reference=1.0, dpi=80)
        assert os.path.exists(out), "output file was not created"
        size = os.path.getsize(out)
        assert size > 0, "output file is empty"

        # A file missing a required column must fail cleanly.
        bad = os.path.join(tmp, "bad.csv")
        with open(bad, "w") as handle:
            handle.write("curve,point,L\n0,0,1.0\n")
        try:
            read_landscape(bad)
        except LandscapeError as exc:
            assert "normU" in str(exc)
        else:  # pragma: no cover - defensive
            raise AssertionError("missing column was not detected")

        # A non-numeric field must fail cleanly too.
        bad2 = os.path.join(tmp, "bad2.csv")
        with open(bad2, "w") as handle:
            handle.write("curve,point,L,normU\n0,0,nonsense,1.0\n")
        try:
            read_landscape(bad2)
        except LandscapeError as exc:
            assert "line 2" in str(exc)
        else:  # pragma: no cover - defensive
            raise AssertionError("malformed value was not detected")

        print(text)
        print("selftest: wrote {} bytes to {}".format(size, out))

    _selftest_norms()
    print("selftest: spline basis / Gram matrices / state norms OK")
    _selftest_loci()
    print("selftest: locus classifier / merge / window / jump cut OK")
    _selftest_reports()
    print("selftest: norm-dependent metric claims OK")
    _selftest_cli()
    print("selftest: command line guards / basis-degree refusal OK")
    print("selftest: OK")
    return 0


class _RecordExplicit(argparse.Action):
    """Store a value and record that the option was actually given.

    ``--norm l2`` and ``--basis-degree 2`` are indistinguishable from
    their defaults *by value*, so a "requires --merge-loci" guard written
    as ``args.norm != 'l2'`` lets exactly those two through while
    refusing every other value -- an inconsistent guard.  This action
    records the option in ``namespace._explicit`` instead, so the guard
    can test whether the user typed it.
    """

    def __call__(self, parser, namespace, values, option_string=None):
        """Set ``dest`` and add it to the explicit-options set."""
        explicit = getattr(namespace, "_explicit", None)
        if explicit is None:
            explicit = set()
            setattr(namespace, "_explicit", explicit)
        explicit.add(self.dest)
        setattr(namespace, self.dest, values)


def _build_parser() -> argparse.ArgumentParser:
    """Construct the command line parser.

    Returns
    -------
    argparse.ArgumentParser
        The configured parser.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Plot the lambda-|u| equilibrium diagram of a gsALMLandscape "
            "landscape.csv export."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "csv_path",
        nargs="?",
        help="path to landscape.csv (omit only with --selftest)",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="output image path (default: the CSV path with a .png suffix)",
    )
    parser.add_argument(
        "--reference",
        type=float,
        default=None,
        help=(
            "reference lambda to draw and compare against, e.g. 6.808124423 "
            "for the Bratu fold"
        ),
    )
    parser.add_argument("--title", default=None, help="figure title")
    parser.add_argument("--dpi", type=int, default=150, help="output resolution")
    parser.add_argument(
        "--selftest",
        action="store_true",
        help="render a small synthetic landscape in a temp dir and self-check",
    )
    parser.add_argument(
        "--merge-loci",
        action="store_true",
        help=(
            "classify every point by its solution profile (from the "
            "landscape .h5 next to the CSV; needs h5py) and merge all "
            "segments of one solution locus into a single polyline, "
            "thesis line-style"
        ),
    )
    parser.add_argument(
        "--window",
        type=float,
        default=None,
        metavar="MU_MIN",
        help=(
            "with --merge-loci: clip to mu >= MU_MIN before merging (the "
            "thesis's stopping condition, eq. 9.1; e.g. 0.01)"
        ),
    )
    parser.add_argument(
        "--norm",
        choices=list(_STATE_NORMS),
        default="l2",
        action=_RecordExplicit,
        help=(
            "with --merge-loci: state-norm convention of the y-axis and of "
            "every printed deviation.  'l2' is the continuous norm of the "
            "spline solution, sqrt(int psi^2) (default); 'h1' adds the "
            "derivative term, sqrt(int psi^2 + int psi'^2) -- a "
            "mode-content diagnostic, NOT a comparison against the "
            "thesis figure; 'discrete' is the historical "
            "||U||_2/sqrt(N_dof), which over-weights the two end "
            "coefficients and is exact on constant states only; 'thesis' is "
            "the thesis's own discrete convention on 100 uniform points"
        ),
    )
    parser.add_argument(
        "--basis-degree",
        type=int,
        default=_BASIS_DEGREE,
        metavar="P",
        action=_RecordExplicit,
        help=(
            "with --merge-loci: spline degree assumed when reconstructing "
            "the basis for the L2/H1/thesis norms (the landscape HDF5 "
            "stores no basis); the element count is N_dof - P and must be "
            "a power of two, or the value is refused"
        ),
    )
    parser.add_argument(
        "--thesis-metrics",
        action="store_true",
        help=(
            "with --merge-loci: print per-locus deviations (A vs the "
            "analytic constant branch; B/C vs the digitized thesis figure "
            "when --thesis-fig is given)"
        ),
    )
    parser.add_argument(
        "--thesis-fig",
        default=None,
        metavar="PNG",
        help=(
            "with --thesis-metrics: page render containing Wouters (2019) "
            "fig. 9.3; curves B/C are digitized from it and written as "
            "CSV next to it"
        ),
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Command line entry point.

    Parameters
    ----------
    argv : sequence of str, optional
        Argument vector; ``sys.argv[1:]`` when omitted.

    Returns
    -------
    int
        ``0`` on success, ``2`` on a user-facing error.
    """
    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.selftest:
        try:
            return _selftest()
        except AssertionError as exc:
            print("selftest FAILED: {}".format(exc), file=sys.stderr)
            return 2
        except LandscapeError as exc:
            print("selftest FAILED: {}".format(exc), file=sys.stderr)
            return 2

    if not args.csv_path:
        parser.error("csv_path is required unless --selftest is given")

    if args.dpi <= 0:
        parser.error("--dpi must be positive")

    if args.window is not None and not args.merge_loci:
        parser.error("--window requires --merge-loci")
    if args.thesis_metrics and not args.merge_loci:
        parser.error("--thesis-metrics requires --merge-loci")
    if args.thesis_fig is not None and not args.thesis_metrics:
        parser.error("--thesis-fig requires --thesis-metrics")
    # Tested by "was it given", not by "does it differ from the default":
    # --norm l2 and --basis-degree 2 would slip through a value test.
    explicit = getattr(args, "_explicit", set())
    if "norm" in explicit and not args.merge_loci:
        parser.error("--norm requires --merge-loci")
    if "basis_degree" in explicit and not args.merge_loci:
        parser.error("--basis-degree requires --merge-loci")
    if args.basis_degree < 1:
        parser.error("--basis-degree must be at least 1")

    try:
        curves = read_landscape(args.csv_path)
        output = args.output if args.output else _default_output(args.csv_path)
        title = args.title if args.title else os.path.basename(args.csv_path)
        if args.merge_loci:
            base, _ = os.path.splitext(args.csv_path)
            h5_path = base + ".h5"
            if not os.path.exists(h5_path):
                raise LandscapeError(
                    "'{}' not found; --merge-loci needs the HDF5 solution "
                    "vectors next to the CSV (a norm-only classification "
                    "cannot separate the loci)".format(h5_path)
                )
            vectors = read_solution_vectors(h5_path, curves)
            loci, notes = merge_loci(
                curves,
                vectors,
                window=args.window,
                norm=args.norm,
                degree=args.basis_degree,
            )
            plot_loci(loci, output, title=title, dpi=args.dpi, norm=args.norm)
            print(summarize_loci(loci, notes, norm=args.norm))
            if args.thesis_metrics:
                print("")
                print(
                    thesis_metrics(
                        loci, thesis_fig=args.thesis_fig, norm=args.norm
                    )
                )
        else:
            plot_landscape(
                curves,
                output,
                title=title,
                reference=args.reference,
                dpi=args.dpi,
            )
            print(summarize(curves, reference=args.reference))
    except LandscapeError as exc:
        print("error: {}".format(exc), file=sys.stderr)
        return 2

    print("")
    print("wrote {}".format(output))
    return 0


if __name__ == "__main__":
    sys.exit(main())
