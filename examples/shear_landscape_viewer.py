#!/usr/bin/env python3
"""Interactive viewer for a ``gsALMLandscape`` shear-wrinkling export.

``example_ShearLandscape`` writes, into ``ShearLandscapeResults/``:

``landscape.csv``
    one row per equilibrium point: ``curve, point, L, normU, stability,
    negatives, isBifurcation, parentCurve, parentPointIdx, ...``
``landscape_amplitude.csv``
    the same rows plus the wrinkling scalars ``maxAbsZ`` and ``halfWaves``
``landscape.h5``
    per point, the deformed shell as a ``gsMultiPatch`` -- a bicubic
    ``TensorBSpline2`` control net (``coefs``) plus both knot vectors.

This script draws the equilibrium diagram of every traced branch and, beside
it, the deformed surface at the selected point, with sliders over branches and
over points along the selected branch.

The surface is evaluated directly from the control net: for a tensor-product
B-spline of degree (p0, p1),

    S(u,v) = sum_i sum_j N_{i,p0}(u) N_{j,p1}(v) C[j,i]

with the coefficients stored direction-0-fastest, i.e. flat index
``k = i + n0 * j``.  Evaluation is O(n_samples^2 * (p0+1) * (p1+1)) per point
via ``scipy.interpolate.BSpline``; the control nets here are 7x7, so a redraw
is dominated by matplotlib, not by the spline.

Two different ``max|z|`` values are shown, and they are *supposed* to differ:
the CSV column ``maxAbsZ`` is measured over the **control net**, while the
panel title reports the maximum of the **evaluated surface**.  A B-spline is
contained in the convex hull of its control polygon, so the surface value is
the smaller of the two (for the r=2 landscape, 0.0838 against 0.1022).  Neither
is wrong; they are different measurements.

Only NumPy, SciPy, Matplotlib and h5py are used.

Examples
--------
Interactive::

    python3 shear_landscape_viewer.py ShearLandscapeResults

Headless contact sheet (one PNG per branch, no GUI needed)::

    python3 shear_landscape_viewer.py ShearLandscapeResults --save out/
"""

from __future__ import annotations

import argparse
import csv
import pathlib
import sys

import numpy as np

# ----------------------------------------------------------------------------
# data loading
# ----------------------------------------------------------------------------


def read_landscape_csv(path):
    """Return ``{column: np.ndarray}`` from a gsALMLandscape CSV export.

    ``writeCsv`` sets no stream precision, so the float columns carry only
    about six significant figures.  Nothing derived from them here is more
    accurate than that.
    """
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        raise SystemExit(f"{path}: no data rows")
    cols = {}
    for key in rows[0]:
        raw = [r[key] for r in rows]
        try:
            cols[key] = np.array([float(v) for v in raw])
        except ValueError:
            cols[key] = np.array(raw, dtype=object)
    return cols


def load_points(results_dir):
    """Prefer landscape_amplitude.csv; fall back to landscape.csv.

    The amplitude file is a superset (it adds ``maxAbsZ`` / ``halfWaves``), but
    only the shear driver writes it -- the Bratu drivers write the plain file.
    """
    results_dir = pathlib.Path(results_dir)
    amp = results_dir / "landscape_amplitude.csv"
    base = results_dir / "landscape.csv"
    if amp.is_file():
        return read_landscape_csv(amp), True
    if base.is_file():
        return read_landscape_csv(base), False
    raise SystemExit(
        f"{results_dir}: neither landscape_amplitude.csv nor landscape.csv found.\n"
        "Run example_ShearLandscape first (see --help)."
    )


class Geometry:
    """Lazy accessor for the per-point deformed patches in landscape.h5."""

    def __init__(self, h5path):
        self.path = pathlib.Path(h5path)
        self._f = None
        self.available = self.path.is_file()
        if self.available:
            try:
                import h5py  # noqa: F401
            except ImportError:
                self.available = False
                self.reason = "h5py is not installed"
                return
        else:
            self.reason = f"{self.path.name} not found"

    def _file(self):
        if self._f is None:
            import h5py

            self._f = h5py.File(self.path, "r")
        return self._f

    def patch(self, curve, point):
        """Control net + knots for one point, or None when not stored.

        ``example_ShearLandscape`` does not attach geometry to every point; the
        landscape records that per point in ``curve_XXXX_hasgeom``.
        """
        if not self.available:
            return None
        f = self._file()
        key = f"curve_{int(curve):04d}_point_{int(point):04d}_mp/patches/patch_0"
        if key not in f:
            return None
        g = f[key]
        return (
            np.asarray(g["coefs"]),
            np.asarray(g["knots_0"]),
            np.asarray(g["knots_1"]),
            int(g.attrs["degree_0"]),
            int(g.attrs["degree_1"]),
        )


# ----------------------------------------------------------------------------
# tensor B-spline evaluation
# ----------------------------------------------------------------------------


def eval_surface(coefs, kv0, kv1, p0, p1, n=41):
    """Evaluate a tensor-product B-spline patch on an n x n parametric grid.

    Returns ``(X, Y, Z)``, each ``(n, n)``.  ``coefs`` is ``(n0*n1, 3)`` stored
    direction-0-fastest (flat index ``k = i + n0 * j``), which is G+Smo's
    ``gsTensorBSpline`` convention.
    """
    from scipy.interpolate import BSpline

    n0 = len(kv0) - p0 - 1
    n1 = len(kv1) - p1 - 1
    if n0 * n1 != coefs.shape[0]:
        raise ValueError(
            f"control net {coefs.shape[0]} != n0*n1 = {n0}*{n1}; "
            "knot vectors and coefficient count disagree"
        )
    # C[j, i, :] -- direction 0 is the fast index, so the (n1, n0) C-order
    # reshape puts j (direction 1) outermost, which is what the flat index
    # k = i + n0*j means.
    C = coefs.reshape(n1, n0, 3)

    u = np.linspace(kv0[p0], kv0[-p0 - 1], n)
    v = np.linspace(kv1[p1], kv1[-p1 - 1], n)
    # Collocation matrices: B0 is (n, n0), B1 is (n, n1).
    B0 = np.empty((n, n0))
    for i in range(n0):
        e = np.zeros(n0)
        e[i] = 1.0
        B0[:, i] = BSpline(kv0, e, p0, extrapolate=False)(u)
    B1 = np.empty((n, n1))
    for j in range(n1):
        e = np.zeros(n1)
        e[j] = 1.0
        B1[:, j] = BSpline(kv1, e, p1, extrapolate=False)(v)
    B0 = np.nan_to_num(B0)
    B1 = np.nan_to_num(B1)

    # S[a, b, :] = sum_j B1[b, j] * sum_i B0[a, i] * C[j, i, :]
    tmp = np.einsum("ai,jic->ajc", B0, C)   # (n, n1, 3)
    S = np.einsum("bj,ajc->abc", B1, tmp)   # (n, n,  3)
    return S[..., 0], S[..., 1], S[..., 2]


# ----------------------------------------------------------------------------
# figure
# ----------------------------------------------------------------------------


def curve_slices(cols):
    """Return ``[(curve_id, row_indices_ordered_by_point), ...]``."""
    cid = cols["curve"].astype(int)
    pid = cols["point"].astype(int)
    out = []
    for c in sorted(set(cid.tolist())):
        idx = np.where(cid == c)[0]
        out.append((c, idx[np.argsort(pid[idx])]))
    return out


def build_figure(cols, geom, has_amp, results_dir):
    import matplotlib.pyplot as plt
    from matplotlib.widgets import Slider

    slices = curve_slices(cols)
    fig = plt.figure(figsize=(13, 6.2))
    fig.suptitle(f"shear-wrinkling landscape — {results_dir}", fontsize=11)
    ax_d = fig.add_subplot(2, 2, 1)
    ax_s = fig.add_subplot(1, 2, 2, projection="3d")
    ax_z = fig.add_subplot(2, 2, 3)
    fig.subplots_adjust(bottom=0.22, top=0.90, wspace=0.15, hspace=0.45)

    # static backdrop: every branch, faint
    for c, idx in slices:
        ax_d.plot(cols["normU"][idx], cols["L"][idx], "-", color="0.82",
                  lw=1.0, zorder=1)
    ax_d.set_xlabel(r"$\|u\|$", labelpad=1)
    ax_d.set_ylabel(r"$\lambda$  (engineering shear strain)")
    ax_d.grid(alpha=0.3)

    hi, = ax_d.plot([], [], "-", color="C0", lw=2.0, zorder=3)
    bif = ax_d.plot([], [], "o", color="C3", ms=6, mfc="none", zorder=4)[0]
    here, = ax_d.plot([], [], "o", color="C1", ms=9, zorder=5)
    readout = ax_d.text(0.02, 0.98, "", transform=ax_d.transAxes, va="top",
                        family="monospace", fontsize=8.5)

    ax_c = fig.add_axes([0.13, 0.09, 0.33, 0.03])
    ax_p = fig.add_axes([0.13, 0.03, 0.33, 0.03])
    s_c = Slider(ax_c, "branch", 0, len(slices) - 1, valinit=0, valstep=1)
    s_p = Slider(ax_p, "point", 0, max(len(slices[0][1]) - 1, 1), valinit=0,
                 valstep=1)

    state = {"warned": False}

    def draw(_=None):
        ci = int(s_c.val)
        c, idx = slices[ci]
        n = len(idx)
        if s_p.valmax != n - 1:
            s_p.valmax = max(n - 1, 1)
            s_p.ax.set_xlim(s_p.valmin, s_p.valmax)
            if s_p.val > n - 1:
                s_p.set_val(n - 1)
                return
        pi = min(int(s_p.val), n - 1)
        row = idx[pi]

        hi.set_data(cols["normU"][idx], cols["L"][idx])
        b = idx[cols["isBifurcation"][idx] > 0.5]
        bif.set_data(cols["normU"][b], cols["L"][b])
        here.set_data([cols["normU"][row]], [cols["L"][row]])

        txt = [f"branch {c}   point {pi+1}/{n}",
               f"lambda  {cols['L'][row]:.6g}",
               f"|u|     {cols['normU'][row]:.6g}",
               f"stab    {int(cols['stability'][row]):+d}"]
        if has_amp:
            txt.append(f"max|z|  {cols['maxAbsZ'][row]:.6g}  (net)")
            txt.append(f"half-w  {int(cols['halfWaves'][row])}")
        pc = int(cols["parentCurve"][row])
        txt.append("seed curve" if pc < 0
                   else f"from curve {pc} pt {int(cols['parentPointIdx'][row])}")
        readout.set_text("\n".join(txt))
        ax_d.relim()
        ax_d.autoscale_view()

        ax_s.clear()
        ax_z.clear()
        ax_z.set_xlabel("x  (mid-width station)", fontsize=8)
        ax_z.set_ylabel("z", fontsize=8)
        ax_z.tick_params(labelsize=7)
        ax_z.grid(alpha=0.3)
        ax_z.axhline(0.0, color="0.6", lw=0.8)
        patch = geom.patch(c, cols["point"][row])
        if patch is None:
            ax_s.set_axis_off()
            msg = ("no geometry stored for this point"
                   if geom.available else f"no surface: {geom.reason}")
            ax_s.text2D(0.5, 0.5, msg, ha="center", transform=ax_s.transAxes,
                        fontsize=9, color="0.4")
        else:
            X, Y, Z = eval_surface(*patch, n=81)
            zmax = float(np.abs(Z).max())
            # Mid-width profile: the driver counts half-waves as (interior sign
            # changes + 1) on exactly this kind of station, ignoring samples too
            # small to carry signal relative to the profile's own amplitude.
            prof = Z[:, Z.shape[1] // 2]
            ax_z.plot(X[:, X.shape[1] // 2], prof, "-", color="C0", lw=1.4)
            sig = prof[np.abs(prof) > 0.02 * max(np.abs(prof).max(), 1e-30)]
            nsc = int((np.diff(np.sign(sig)) != 0).sum()) if sig.size else 0
            ax_z.set_title(f"mid-width profile — {nsc} interior sign change(s)"
                           f"  =>  {nsc + 1} half-wave(s)", fontsize=8)
            lim = zmax if zmax > 1e-12 else 1.0
            ax_s.plot_surface(X, Y, Z, cmap="coolwarm", vmin=-lim, vmax=lim,
                              linewidth=0, antialiased=True, rcount=41,
                              ccount=41)
            ax_s.set_zlim(-max(lim, 1e-3) * 1.2, max(lim, 1e-3) * 1.2)
            ax_s.set_title(f"deformed shell   max|z| = {zmax:.4g}  (surface)",
                           fontsize=9)
            ax_s.set_xlabel("x"); ax_s.set_ylabel("y")
        fig.canvas.draw_idle()

    s_c.on_changed(draw)
    s_p.on_changed(draw)
    draw()
    return fig, s_c, s_p, draw, slices


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__.split("\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Run example_ShearLandscape first; it writes ShearLandscapeResults/.")
    ap.add_argument("results_dir", nargs="?", default="ShearLandscapeResults",
                    help="directory written by example_ShearLandscape "
                         "(default: %(default)s)")
    ap.add_argument("--save", metavar="DIR",
                    help="headless: write one PNG per branch to DIR and exit")
    ap.add_argument("--point", type=int, default=None,
                    help="with --save, the point index to render per branch "
                         "(default: the point of largest max|z|)")
    args = ap.parse_args(argv)

    import matplotlib
    if args.save:
        matplotlib.use("Agg")

    cols, has_amp = load_points(args.results_dir)
    geom = Geometry(pathlib.Path(args.results_dir) / "landscape.h5")
    if not geom.available:
        print(f"note: surfaces unavailable ({geom.reason}); "
              "the diagram still works", file=sys.stderr)

    fig, s_c, s_p, draw, slices = build_figure(cols, geom, has_amp,
                                               args.results_dir)

    if args.save:
        out = pathlib.Path(args.save)
        out.mkdir(parents=True, exist_ok=True)
        key = "maxAbsZ" if has_amp else "normU"
        for ci, (c, idx) in enumerate(slices):
            s_c.set_val(ci)
            pick = (args.point if args.point is not None
                    else int(np.argmax(cols[key][idx])))
            s_p.set_val(min(pick, len(idx) - 1))
            draw()
            p = out / f"branch_{c:02d}.png"
            fig.savefig(p, dpi=130)
            print(f"wrote {p}")
        return 0

    import matplotlib.pyplot as plt
    plt.show()
    return 0


if __name__ == "__main__":
    sys.exit(main())
