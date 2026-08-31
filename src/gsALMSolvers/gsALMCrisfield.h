 /** @file gsALMCrisfield.h

    @brief Performs the arc length method to solve a nonlinear equation system.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#pragma once

#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMBase.h>

namespace gismo
{

/**
    @brief Performs the Crisfield arc length method to solve a nonlinear equation system.

    == The BORDERED-SOLVE CORRECTOR (option \c BorderedMode, default "Off") ==

    \c BorderedMode selects, at RUNTIME, which corrector chart \a step() uses and in which
    order: \c "Off" (default) never touches the bordered chart; \c "Fallback" attempts the
    elimination first and retries with the bordered chart on failure; \c "Primary" attempts
    the bordered chart first (the references' own practice -- both AUTO and pde2path use the
    bordered (n+1) solve as THE corrector, not a fallback) and retries with the elimination.
    The old switch \c BorderedFallback is a DEPRECATED ALIAS for \c BorderedMode="Fallback",
    kept working and warning once per solver object; see \a getOptions().

    (a) WHY. Near a fold the corrector's failure is a CANCELLATION failure of the
    \f$\delta u_t/\bar{u}\f$ ELIMINATION, and NOT a singular-matrix failure. Measured on the
    2-DOF fold fixture F at \c Scaling=0 and \f$ds=0.05\f$, where the predictor of step 19
    lands on \f$u_1 = 1+2.2\cdot10^{-16}\f$, i.e. \f$K = \mathrm{diag}(-4.4\cdot10^{-16},1)\f$:
    \c SimplicialLDLT factorizes that tangent WITHOUT complaint (the reported status is
    \c gsStatus::NotConverged = 1, never \c SolverError = 3), and
    \f$|\delta u_t| = 2.25\cdot10^{15}\f$, \f$|\bar u| = 5.63\cdot10^{12}\f$ are both
    individually "correct". What dies is their combination: with
    \f$a_0 = |\delta u_t|^2\f$, \f$b_1 = 2\bar u\cdot\delta u_t\f$, \f$c_2 = |\bar u|^2\f$
    the leading part of the constraint quadratic is the PERFECT SQUARE
    \f$(\alpha\,\delta\Lambda+\beta)^2/\varepsilon^2\f$, so its discriminant
    \f$\alpha_2^2-4\alpha_1\alpha_3\f$ is a difference of two \f$O(\varepsilon^{-4})\f$
    numbers with true value \f$O(\varepsilon^{-2})\f$ -- a relative cancellation of
    \f$\varepsilon^2 \approx 10^{-32}\f$, far below double precision. MEASURED: the
    discriminant comes out EXACTLY 0, both roots collapse onto the degenerate double root,
    the resulting point violates \f$|\Delta U| = ds\f$ by three orders of magnitude, and the
    corrector settles into a bounded PERIOD-2 limit cycle (\f$u_1\f$ alternating
    \f$1.0 \leftrightarrow 0.9502\f$, \f$|\delta U|\f$ and the residual constant to all
    printed digits) until \c MaxIter is exhausted. By contrast, the bordered matrix \a B
    below is NONSINGULAR at exactly that point -- with \f$K=\mathrm{diag}(0,1)\f$,
    \f$f\notin\mathrm{range}(K)\f$ and border row \f$w = 2\Delta U\f$, \a B has full rank.

    (b) WHAT IT IS. An EXACT CHANGE OF BASIS of the same affine solution set of the
    underdetermined \f$K\,\delta U - \delta\Lambda f = -r\f$, not an inexact-Newton
    perturbation and not a regularisation. With
    \f[ B = \begin{bmatrix} K & -f \\ w^T & \gamma\end{bmatrix},\quad
        w = 2\Delta U,\ \gamma = 2A_0\Delta\Lambda \f]
    (i.e. \a B is the exact Jacobian of \f$\{R=0,\ c=0\}\f$, the very system the corrector
    solves, so "\a B nonsingular" IS "this Crisfield step is well posed") we solve
    \f$B[p_t;q_t]=[0;1]\f$ and \f$B[p;q]=[-r;0]\f$ and parameterise
    \f$\delta U = \eta p + t\,p_t\f$, \f$\delta\Lambda = \eta q + t\,q_t\f$. Whenever \a K is
    invertible, \f$q_t = 1/(w^T\delta u_t+\gamma)\f$, \f$p_t = q_t\,\delta u_t\f$,
    \f$q = -q_t (w^T\bar u)\f$, \f$p = \bar u + q\,\delta u_t\f$, so the two charts describe
    the same set and produce the same accepted point in exact arithmetic. The difference is
    only that \f$(p_t,q_t)\f$ is \f$(\delta u_t,1)\f$ NORMALISED by the \f$O(1/\varepsilon)\f$
    factor \f$w^T\delta u_t+\gamma\f$, so both bordered solutions stay \f$O(1)\f$ and nothing
    cancels.

    (c) THE COEFFICIENT MAP is term-for-term (see \a _computeLambdasBordered), and
    \a computeLambdasModified() and \a computeLambdasEta() are INVARIANT and deliberately
    left unchanged: the \a t-quadratic is the \f$\delta\Lambda\f$-quadratic composed with the
    affine map \f$\delta\Lambda=\eta q+t q_t\f$, hence \f$D_t = q_t^2 D_\Lambda\f$ as
    polynomials in \f$\eta\f$ -- the sign of the discriminant, the \f$\eta\f$ roots and
    \f$-\alpha_2/\alpha_1\f$ are the SAME decisions.

    (d) ⚠ TRAP: the INDEX ORDER FLIPS when \f$q_t<0\f$. Both quadratics have a positive
    leading coefficient, so root 0 is the larger in each parameter, but
    \f$\delta\Lambda=\eta q+t q_t\f$ reverses the order for negative \f$q_t\f$. Selection is
    by comparison and never by index, so this is harmless HERE -- it is fatal for anyone who
    assumes \c m_deltaLs[0] came from \c t[0].

    (e) THIS PATH NEVER EMITS \c SolverError. The step-fail retry in
    \a gsALMExploration::traceCurve() retries only on \c NotConverged || \c AssemblyError, so
    a \c SolverError falls through to its "record the point" branch. Every failure here -- LU failure, singular \a B, non-finite
    \f$q_t\f$, a complex root -- warns and \c throw 1, and \a step() additionally clamps the
    reported status so attempt 2 can never be worse than attempt 1.

    (f) THE OPTION IS OPT-IN: \c BorderedMode DEFAULTS \c "Off", and \a defaultOptions() is
    the authority on that literal -- read it there, not this sentence. It was briefly
    defaulted to a bordered-first mode on 2026-07-30 as the package partner
    of \c Scaling=0; BOTH of those default changes were reverted the same day, so the shipped
    pair is \c Scaling=-1 with \c BorderedMode="Off". Turning it ON (\c "Fallback" or
    \c "Primary") remains a good idea for a caller that hits fold-step failures.
    \c "Fallback" attempts the elimination first and retries with the bordered chart only
    after the elimination fails; \c "Primary" attempts the bordered chart FIRST and retries
    with the elimination -- the references' own practice, see the class-doc heading. Under
    \c "Primary" the retry (the elimination) CAN emit \c SolverError even though the bordered
    chart itself never does (e); \a step() clamps the reported status to attempt 1's own
    status either way, so this asymmetry is invisible to the caller.
    The deprecated switch \c BorderedFallback remains registered and, when \c BorderedMode is
    left at \c "Off", is mapped onto \c "Fallback" with a once-per-object warning; an
    explicitly non-\c "Off" \c BorderedMode always wins over the alias (see \a getOptions()).
    The trigger for the ONE retry is a step-level failure -- \c NotConverged or
    \c SolverError, see \a step() -- so no new code runs on any path that succeeds under
    \c "Off". \c BorderedMode != "Off" is NOT bit-neutral on runs that DO fail steps, and
    that is deliberate: an \c example_BratuExploration Crisfield run halves the arc length
    ~28 times at the \f$\lambda\approx5.91\f$ barrier and every one of those halvings is a step that
    returned \c NotConverged, so the retry fires dozens of times per run. MEASURED WITH THE
    SWITCH ON, on the three Crisfield landscape CSVs: the \c -m 1 landscape loses exactly one
    point -- a successful retry replacing a halving -- with \f$\lambda_{\max}=5.90762\f$
    BIT-IDENTICAL, the \a gsALMConsistentCrisfield landscape is bit-identical (containment by
    class), and no acceptance check changes verdict. Why the two go together: with
    \c Scaling=0 and this switch OFF, the 2-DOF fold fixture F dies at step 19 with
    \c NotConverged and \f$\max u_1 = 0.95\f$ (measured, at all three load scales); with both
    on it rounds the fold to \f$\max u_1 = 1.8\f$ at \f$ds = 0.05\f$ unchanged.

    (g) NOT COVERED. Branch points (\f$K\f$ singular AND \f$f\in\mathrm{range}(K)\f$): there
    the null space is 2-dimensional and \a B is singular for EVERY \a w; that case belongs to
    \a gsALMBase::_extendedSystemSolve, and here it correctly fails as \c NotConverged. The
    three sibling solvers \a gsALMRiks, \a gsALMConsistentCrisfield and \a gsALMLoadControl
    carry the identical elimination and the identical exposure and are
    deliberately untouched. The complex-root fallback \a computeLambdasComplex() is not
    re-derived in the bordered chart: it is guarded by a warn+throw, on the measured evidence
    that it does not fire anywhere on the step-19 trace at any of the three load scales.

    (h) ★ THE TRANSVERSALITY GUARD of \a _borderedSolve. The identity in (b),
    \f$q_t = 1/(w^T\delta u_t+\gamma)\f$, says that the bordered chart itself degenerates
    when \f$b_0 := w^T\delta u_t+\gamma\f$ vanishes: \f$b_0\f$ is the SCHUR COMPLEMENT of
    \a K in \a B, hence \f$\det B = \det(K)\,b_0\f$, and with \a K nonsingular \a B is
    singular \b iff \f$b_0=0\f$. This is exactly the case (2.20a)/(2.20b) split of
    Riks, E. (1984), CMAME 47, p. 228, whose eq. (2.31) reads
    \f$J = (D/\lambda')(n^t\cdot x')\f$ with \f$D=\det K\f$: our \f$b_0\f$ IS his
    \f$(n^t\cdot x')\f$, our \a w his \a n, our \f$\det B\f$ his \a J. MEASURED:
    \f$K=I_2\f$, \f$f=(1,0)\f$ and a corrector state \f$\Delta U=(1.4\cdot10^{-9},0.13)\f$
    reached by this class's own predictor + one bordered iteration gives
    \f$q_t=3.65\cdot10^8\f$, which \c SparseLU returns without any complaint and both
    existing guards accept -- and the point the corrector then lands on misses
    \f$|\Delta U| = ds\f$ by a factor TWO.
    - DETECT on the \b scale-free cosine of the bordering direction \f$d=2(\Delta U,
      \Delta\Lambda)\f$ against the augmented tangent \f$x'=(\delta u_t,1)\f$, in the
      constraint metric \f$\langle(a,\alpha),(b,\beta)\rangle_A = a\cdot b + A_0\alpha\beta\f$
      in which the corrector's own quadratic is written. Since the last row of \a B gives
      \f$\langle d,(p_t,q_t)\rangle_A = 1\f$ identically, that cosine is
      \f$1/(\|d\|_A\,\|(p_t,q_t)\|_A)\f$ -- see \a _chartTransversality: no cancellation,
      no division by the vanishing quantity, and NO second copy of \f$b_0\f$.
    - DISTINGUISH Riks's (2.20a) (\a G singular, \a G* full rank -- a CHART artifact,
      \a removable) from (2.20b) (\a G and \a G* both rank deficient -- a BIFURCATION,
      which a "re-border and continue" would silently walk past). His summary (2.29) is
      the discriminator: \f$b_0\neq0\f$ at a limit point, so \f$b_0\approx0\f$ TOGETHER
      with a rank-deficient \a K is the bifurcation signature. THREE tests, in this order:
      -# RE-BORDER, which is Riks's own remedy -- "it can always be avoided if we take
         another \f$f^*_{N+1}\f$ so that \f$n^*\neq n\f$" -- and therefore also his own
         TEST: a GENERIC bordering is transversal with probability one (Govaerts, W. J. F.
         (2000), \a Numerical \a Methods \a for \a Bifurcations \a of \a Dynamical
         \a Equilibria, SIAM, Prop. 3.2.1 p. 51: "a generic choice of B, C, D will do ...
         could be generated by a random number generator"), so a degeneracy that SURVIVES
         it is not a chart artifact. Deterministic: the fixed-seed LCG of
         \a gsALMBase::_computeCriticalMode -- see \a _genericBorderDirection.
      -# Is \a K numerically RANK DEFICIENT (Riks's \f$D=0\f$)? \a _chartSingularityProbe,
         run THROUGH the re-bordered chart because it is the only one accurate enough.
      -# If it is, the second sourced criterion decides what to say: \f$\varphi^T P\f$ of
         Wriggers & Simo (1990) eqs. (9a)/(9b). AGREEMENT with (2.29) is reported as a
         candidate bifurcation; DISAGREEMENT is reported as unclassifiable. Neither
         re-borders.
    - RESPOND: only a confirmed (2.20a) -- degeneracy removed by the generic bordering AND
      a well-conditioned tangent -- continues, on the re-bordered chart, which describes
      the SAME affine solution set. Every other outcome warns and \c throw 1, i.e. reports
      the step rather than handing back a chart it cannot compute in.
    - WHY MONITOR AT ALL is Govaerts, W. (1991), SIAM J. Matrix Anal. Appl. 12(3),
      469-483, Prop. 3.1 / Cor. 3.2: the bordered solve's accuracy is governed by
      \f$\kappa(M)\f$ of the BORDERED matrix. ⚠ Note the orientation: that paper's setting
      is \a K nearly singular with \a M well conditioned -- the ordinary fold, our (b).
      The case here is the COMPLEMENTARY one, \a K perfectly conditioned and \a M singular
      through \f$b_0\to0\f$, which it does not study; what it supplies is the criterion.
    - ⚠ THE TWO NUMERICAL THRESHOLDS ARE OURS AND UNSOURCED. Riks's condition is exact
      (\f$=0\f$) and Govaerts's criterion is \f$\kappa(M)\f$; no surveyed source gives a
      cutoff. See \a m_chartCosTol / \a m_chartRankTol.

    \tparam T coefficient type

    \ingroup gsALMSolvers
*/
template <class T>
class gsALMCrisfield : public gsALMBase<T>
{

    typedef gsALMBase<T> Base;

    typedef typename Base::ALResidual_t  ALResidual_t;
    typedef typename Base::Jacobian_t    Jacobian_t;
    typedef typename Base::dJacobian_t   dJacobian_t;

public:

    using Base::setLength;
    using Base::computeStability;

    /**
     * @brief      One arc-length step, with the RUNTIME-selectable BORDERED-SOLVE corrector.
     *
     * \c BorderedMode (default \c "Off", see \a defaultOptions()) selects which chart attempt
     * 1 uses and whether there is an attempt 2 at all:
     * - \c "Off" -- attempt 1 is \a gsALMBase::step() verbatim, using the elimination chart,
     *   and \a step() returns immediately: no snapshot, no retry, no new code on any path
     *   (bit-neutral by construction).
     * - \c "Fallback" -- attempt 1 is the elimination; if it fails (see the gate below),
     *   attempt 2 retries once with the bordered chart documented on this class.
     * - \c "Primary" -- attempt 1 is the BORDERED chart (the references' own practice, see
     *   the class-doc heading); if it fails, attempt 2 retries once with the elimination.
     *
     * The retry gate is \c (st == gsStatus::NotConverged || st == gsStatus::SolverError):
     * widened past \c NotConverged alone because a failed factorization at a near-singular
     * corrector-iterate tangent is exactly the case the bordered chart exists for (audit
     * B.6), and today's elimination-only path never reaches a retry for it. \a predictor()
     * factorizes \a K on BOTH charts regardless of \c BorderedMode (unchanged, out of scope),
     * so a tangent that is already singular at the step's SEED state still cannot be rescued;
     * the case this gate DOES rescue is a factorization failure inside a corrector iteration
     * (\a quasiNewtonIteration()), where the bordered chart never factorizes \a K at all.
     *
     * A failed retry leaves no trace: \a m_DeltaUold, \a m_DeltaLold, \a m_arcLength,
     * \a m_arcLength_prev, \a m_stepTaken, \a m_stabilityPrev, \a m_phi, and the four
     * stability-bookkeeping members \a m_negatives, \a m_indicator, \a m_stabilityVec and
     * \a m_stability are snapshotted beforehand and restored (see C6 / the fourth member is
     * required because \a gsALMBase::_computeStability() derives \a m_stability from
     * \a m_indicator in the very same call that a failed attempt 2 also runs, so restoring
     * the other three without it would leave the pair mutually inconsistent). The reported
     * status is clamped to attempt 1's OWN status, so attempt 2 can never be reported as
     * WORSE than attempt 1 (see the \c !! CORRECTED comment in the
     * implementation for why this still matters under \c "Primary", where attempt 2 -- the
     * elimination -- CAN emit \c SolverError even though the bordered chart never does, (e)).
     *
     * @note With \c AdaptiveLength on, a SUCCESSFUL retry makes \a computeLength() consume
     *       the RETRY's iteration count, so every subsequent step gets a different arc
     *       length. That is inherent to retrying and is not guarded here.
     */
    gsStatus step();

protected:

    using Base::computeJacobian;
    using Base::computeResidual;
    using Base::computeResidualNorms;
    using Base::computeUt;
    using Base::computeUbar;
    using Base::computeLength;

public:

    /// Constructor
    gsALMCrisfield( const Jacobian_t  &Jacobian,
                    const ALResidual_t&ALResidual,
                    const gsVector<T> &Force )
    : Base(Jacobian,ALResidual,Force)
    {
        defaultOptions();
        getOptions();

        initMethods();
    }

    /// Constructor using the jacobian that takes the solution and the solution step
    gsALMCrisfield( const dJacobian_t &dJacobian,
                    const ALResidual_t&ALResidual,
                    const gsVector<T> &Force )
    : Base(dJacobian,ALResidual,Force)
    {
        defaultOptions();
        getOptions();

        initMethods();
    }

public:
    /// Distance in the (U,L) plane, measured in the constraint metric of the CURRENT
    /// step; the \f$\|f\|^2\f$ scaling is the frozen step forcing (see
    /// gsALMBase::stepForcing), so that this is exactly the metric in which the
    /// corrector's quadratic constraint is written.
    T distance(const gsVector<T>& DeltaU, const T DeltaL) const
    {
        T A0 = math::pow(m_phi,2)*this->stepForcing().dot(this->stepForcing());
        return math::pow(DeltaU.dot(DeltaU) + A0*math::pow(DeltaL,2.0),0.5);
    }

protected:

    /// See gsALMBase
    void initMethods();
    /// See gsALMBase
    void initiateStep();
    /// See gsALMBase
    void iterationFinish();

    /// See gsALMBase
    void quasiNewtonPredictor();
    /// See gsALMBase
    void quasiNewtonIteration();

    /// See gsALMBase
    void predictor();
    void predictorGuess();
    /// See gsALMBase
    void iteration();

    /// See gsALMBase
    void initOutput();
    /// See gsALMBase
    void stepOutput();

    /// See gsALMBase
    void defaultOptions();
    /// See gsALMBase
    void getOptions();

    /// Compute the load factors
    void computeLambdas();
    /// Compute the load factors
    void computeLambdasSimple();
    /// Compute the load factors
    void computeLambdasModified();
    /// Compute the load factors
    void computeLambdasComplex();
    /// Compute the load factors
    void computeLambdasEta();
    /// Compute the load factors
    void computeLambdaDET();
    /// Compute the load factors
    void computeLambdaDOT();
    /// Compute the load factors
    void computeLambdaMU();

    // --- bordered-solve fallback (see the class documentation) -----------------------

    /// Assembles the tangent at the current iterate WITHOUT factorizing it. The ordinary
    /// path cannot be reused: \a gsALMBase::_computeJacobian(), which every
    /// \a computeJacobian() overload delegates to, calls \a factorizeMatrix() internally and
    /// would \c throw 3 on an exactly singular tangent -- the very tangent the bordered chart
    /// exists to handle. Same shape as the \c factorizeShifted() ladder in
    /// \a gsALMBase::_extendedSystemIteration().
    gsSparseMatrix<T> _assembleJacobianUnfactorized();

    /// Builds \f$B=[[K,-f],[w^T,\gamma]]\f$ and performs the TWO bordered solves, filling
    /// \a m_deltaUt = \f$p_t\f$, \a m_deltaLt = \f$q_t\f$, \a m_deltaUbar = \f$p\f$,
    /// \a m_deltaLbar = \f$q\f$. Warns and \c throw 1 (never 3, see (e)) on any failure.
    /// Carries the transversality guard of (h).
    void _borderedSolve();

    // --- the transversality guard, see (h) --------------------------------------------

    /// Assembles \f$B=[[K,-f],[w^T,\gamma]]\f$ from the CURRENT \a m_jacMat and performs the
    /// two solves \f$B[p_t;q_t]=[0;1]\f$, \f$B[p;q]=[-r;0]\f$, writing \a m_deltaUt,
    /// \a m_deltaLt, \a m_deltaUbar, \a m_deltaLbar. Returns \c false -- and never throws --
    /// when the factorization fails or the chart comes back non-finite, so that the CALLER
    /// decides between the historic warn+\c throw 1 (the constraint bordering) and the
    /// (2.20b) verdict (a re-bordering that also fails).
    bool _borderedChartSolve(const gsVector<T> & w, const T gamma, const gsVector<T> & f);

    /// The transversality cosine of (h):
    /// \f$|\langle d,x'\rangle_A|/(\|d\|_A\|x'\|_A) = 1/(\|d\|_A\|(p_t,q_t)\|_A)\f$, using
    /// \f$\langle d,(p_t,q_t)\rangle_A = 1\f$ (the last row of \a B) and
    /// \f$(p_t,q_t) = q_t\,x'\f$. \a du / \a dl are the bordering DIRECTION, whose metric
    /// dual \f$(d_u, A_0 d_\lambda)\f$ is the border row actually used --
    /// \f$d = 2(\Delta U,\Delta\Lambda)\f$ for Crisfield's own constraint. Reads the chart
    /// members, so it is valid only right after \a _borderedChartSolve returned \c true.
    T _chartTransversality(const gsVector<T> & du, const T dl, const T A0) const;

    /// The TWO signals of the (2.20a)/(2.20b) split, both read off the ALREADY factorized
    /// \a B with one extra solve and no factorization of \a K (which the bordered path
    /// deliberately never performs). It solves \f$B[z;\zeta]=[v;0]\f$ for a fixed-seed
    /// generic unit \a v and removes the bordering's own contribution exactly,
    /// \f$z^* = z - (\zeta/q_t)p_t\f$, which satisfies \f$K z^* = v\f$ -- ONE step of inverse
    /// iteration on \a K, obtained through \a B. ⚠ Its accuracy is the chart's: the
    /// subtraction cancels a term \f$1/\cos\f$ larger than the answer, so it must be run on a
    /// TRANSVERSAL chart (measured in (h): through the degenerate chart it reports
    /// \c rank \c = \c 1 on a tangent with \f$\det K = 10^{-12}\f$).
    /// \a rank  = \f$(1/\|z^*\|)/\|K\|_F\f$, an UPPER bound on
    ///            \f$\sigma_{\min}(K)/\|K\|_F\f$: small means \a K is numerically rank
    ///            deficient (Riks's \f$D=0\f$), \b in \b any direction, including one the
    ///            load never excites.
    /// \a modeCos = \f$|\varphi\cdot f|/\|f\|\f$ with \f$\varphi = z^*/\|z^*\|\f$, i.e. the
    ///            limit-vs-branch criterion \f$\varphi^T P\f$ of Wriggers & Simo (1990),
    ///            IJNME 30, p. 159, eqs. (9a)/(9b) (= Wriggers, Wagner & Miehe 1988, eq. (25);
    ///            Spence & Jepson 1984) -- the same statistic
    ///            \a gsALMBase::_testSingularPoint thresholds, recomputed here from a mode
    ///            this class can obtain WITHOUT overwriting \a m_V or refactorizing
    ///            \a m_jacMat. Two honest limitations: \f$\varphi\f$ is ONE step of inverse
    ///            iteration (against up to \c SingularPointTestIt sweeps there), and the
    ///            equivalence uses the LEFT null vector, so it needs the symmetry the default
    ///            \c SimplicialLDLT path already assumes.
    void _chartSingularityProbe(const gsVector<T> & f, T & rank, T & modeCos);

    /// A GENERIC bordering direction (Govaerts 2000, Prop. 3.2.1: "could be generated by a
    /// random number generator"), scaled to \f$\|d\|_A = \a scaleA\f$ so that the re-bordered
    /// \a B keeps the row scale of the one it replaces. DETERMINISM: a fixed-seed 32-bit LCG
    /// with the Numerical Recipes constants of \a gsALMBase::_computeCriticalMode, re-seeded
    /// on every call -- no global RNG state, no dependence on the path taken to this state,
    /// so a restart re-borders identically.
    void _genericBorderDirection(gsVector<T> & du, T & dl, const T scaleA, const T A0) const;

    /// The bordered counterpart of \a computeLambdas(): the SAME six constraint
    /// coefficients under the substitution \f$(\delta u_t,1)\to(p_t,q_t)\f$,
    /// \f$(\bar u,0)\to(p,q)\f$, the same branch structure, and the same root DECISION
    /// (computed on bounded quantities). Stores \f$\delta\Lambda\f$ -- not \a t -- in
    /// \a m_deltaLs, which is a base member with existing readers.
    void _computeLambdasBordered();

protected:

    // Number of degrees of freedom
    using Base::m_numDof;

    using Base::m_jacobian;
    using Base::m_djacobian;
    using Base::m_residualFun;
    using Base::m_forcing;

    /// Solver options
    using Base::m_options;

    /// Number of Arc Length iterations performed
    using Base::m_numIterations;

    /// Maximum number of Arc Length iterations allowed
    using Base::m_maxIterations;

    /// Length of the step in the u,f plane
    using Base::m_arcLength;

    /// Output verbosity
    using Base::m_verbose;

    /// Note
    using Base::m_note;

    /// Convergence result
    using Base::m_converged;

    /// Force residuum
    using Base::m_residueF;

    /// Displacement residuum
    using Base::m_residueU;

    /// Load residuum
    using Base::m_residueL;

    /// Indicator for bifurcation
    using Base::m_indicator;
    using Base::m_negatives;

    /// Option \c SingularPointTestTol, the tolerance \a gsALMBase::_testSingularPoint
    /// thresholds \f$|V\cdot f|/\|f\|\f$ with. The guard of (h) reuses it -- the USER's own
    /// limit-vs-branch angle tolerance -- for the \f$\varphi^T P\f$ signal of
    /// \a _chartSingularityProbe, rather than inventing a second one for the same criterion.
    using Base::m_SPTestTol;

    /// Relaxation factor
    using Base::m_relax;

    // Previous update
    using Base::m_DeltaUold;
    using Base::m_DeltaLold;
    /// Displacement vector (present, at previously converged point)
    using Base::m_U;
    using Base::m_Uprev;
    using Base::m_Uguess;
    /// Update of displacement vector
    using Base::m_DeltaU;
    /// u_bar
    using Base::m_deltaUbar;
    /// u_t
    using Base::m_deltaUt;
    /// Update of update of displacement vector
    using Base::m_deltaU;

    /// Lambda (present, at previously converged point)
    using Base::m_L;
    using Base::m_Lprev;
    using Base::m_Lguess;
    /// Update of lambdaGeneralizedSelfAdjointEigenSolver
    using Base::m_DeltaL;
    /// Update of update of lambda
    using Base::m_deltaL;
    /// Vector with lambda updates
    using Base::m_deltaLs;

    /// Jacobian matrix
    using Base::m_jacMat;
    using Base::m_detKT;

    /// Value of the residual function (right-hand side of the second bordered solve)
    using Base::m_resVec;

    // Angle determination method: 0: determine based on previous load step. 1: determine based on previous iteration
    index_t m_angleDetermine;

    /// Scaling parameter
    T m_phi;
    bool m_phi_user;

    // MODIFIED ARC LENGTH METHOD
    /// factor (modified arc length method)
    T m_eta;

    // discriminant
    T m_discriminant;

    T m_alpha1;
    T m_alpha2;
    T m_alpha3;

    T m_a0;
    T m_b0,m_b1;
    T m_c0,m_c1,m_c2;

    // --- bordered-solve corrector state (see the class documentation) -----------------

    /// Option \c BorderedMode, resolved to one of \a borderedmode's values (same idiom as
    /// \a m_angleDetermine / \a angmethod: stored as \c index_t, compared against
    /// \a borderedmode::Off / \a Fallback / \a Primary). Default \c Off -- see (f).
    /// \a defaultOptions() holds the authoritative literal ("Off").
    index_t m_borderedMode;
    /// Whether the deprecated \c BorderedFallback alias has already warned for THIS solver
    /// object. Seeded \c false in \a defaultOptions() (NOT \a initMethods() and NOT
    /// \a getOptions()) so it survives the repeated \a getOptions() calls every
    /// \a gsALMBase::applyOptions() makes -- see the class \c C1 rationale.
    bool m_borderedAliasWarned;
    /// True ONLY inside whichever attempt of \a step() is running the bordered chart. The
    /// single switch every bordered branch tests; when false, not one line of bordered code
    /// runs. \a step()'s \c Off branch never assigns it (re-asserted \c false by
    /// \a getOptions() instead), which is what lets the \c forceBordered white-box test idiom
    /// keep working.
    bool m_useBordered;
    /// Option \c BorderedSolver (default "LU" = Eigen \c SparseLU, which pivots). Kept
    /// SEPARATE from \a gsALMBase::m_solver, which must keep \a K's factors: the
    /// Determinant branch of \a gsALMBase::_computeStability() recovers them by \c dynamic_cast
    /// and \a computeUbar() reuses them. Besides, \a B is nonsymmetric whenever
    /// \f$w\neq-f\f$, and at \c Scaling=0 its corner \f$\gamma=2A_0\Delta\Lambda\f$ is a
    /// STRUCTURAL ZERO on the diagonal -- no non-pivoting factorization can be used.
    std::string m_borderedSolverName;
    typename gsSparseSolver<T>::uPtr m_borderedSolver;

    // --- transversality-guard state and thresholds, see (h) ---------------------------

    /// The transversality cosine of the chart the last \a _borderedSolve settled on (after a
    /// re-bordering, the RE-BORDERED one). Diagnostic only: nothing in the corrector reads
    /// it. BORDERED-ONLY, poisoned with NaN in \a initMethods() for the same reason
    /// \a m_deltaLt is.
    T m_chartCosine;
    /// How many times the guard has re-bordered since construction. Cumulative on purpose:
    /// a per-solve flag cannot answer "did this ever fire anywhere in the analysis", which is
    /// the question both the unit tests and a user chasing a moved trajectory actually ask.
    index_t m_chartReborderings;

    /// ⚠ OURS AND UNSOURCED (see (h)); NOT options, because no source offers a value to
    /// expose. \a m_chartCosTol: below this cosine the chart is treated as degenerate. The
    /// FORM is the point and follows the C1 precedent -- there an absolute \c 1e-6 on
    /// \f$|V\cdot f|\f$ was replaced by the dimensionless \f$|V\cdot f|/\|f\|\f$, because a
    /// threshold on an unnormalised quantity silently means a relative tolerance that depends
    /// on the model's units; the same applies here, and a EUCLIDEAN cosine would NOT do: it
    /// mixes displacement with load units, so on a stiff model (\f$\|\delta u_t\|\ll1\f$) it
    /// would fire on healthy charts. The VALUE is calibrated, not derived: MEASURED,
    /// the cosine is \b 1 at every step of the fixture-F fold traces this fallback exists for
    /// (all three load scales) and \f$\geq 0.966\f$ at every corrector iteration of a
    /// two-dimensional trace on fixture P's bifurcated branch at a deliberately coarse
    /// \f$ds = 0.5\f$ (both Scalings), while the counterexample of (h) sits at
    /// \f$1.05\cdot10^{-8}\f$ -- 1e-6 is placed in that gap, nearly six orders below the
    /// healthy minimum and two above the defect. Govaerts's \f$\kappa(M)\f$ criterion says
    /// only that the bordered solve is outside every stated guarantee once
    /// \f$1/\cos\gtrsim1/\epsilon_{\mathrm{mach}}\f$; it gives no cutoff.
    T m_chartCosTol;
    /// ⚠ OURS AND UNSOURCED. \a m_chartRankTol: below this value of
    /// \f$\sigma_{\min}(K)/\|K\|_F\f$ the tangent counts as numerically rank deficient, i.e.
    /// Riks's \f$D=0\f$. Deliberately SEVERE (1e-10, ~6 orders below the counterexample's
    /// 0.71): the asymmetry is that a false "rank deficient" only costs a reported step and an
    /// arc-length halving, whereas a false "chart artifact" re-borders past a bifurcation --
    /// the exact failure this guard exists to prevent. A genuinely ill-conditioned tangent on
    /// a (2.20a) chart is therefore REPORTED rather than rescued, on purpose.
    T m_chartRankTol;

    /// \f$q_t\f$ and \f$q\f$, the load-factor components of the two bordered solutions.
    /// ⚠ BORDERED-ONLY: written and read exclusively inside the bordered branch. Their
    /// "natural" non-bordered values would be 1 and 0 (the chart the elimination implies)
    /// and they are deliberately NOT seeded with those: a member initialised as if it were
    /// always meaningful but maintained on only one path is precisely the shape of the
    /// m_DeltaV defect of the ALM round-5 work. They are poisoned with NaN in
    /// \a initMethods() instead, so that any future non-bordered reader fails loudly rather
    /// than silently consuming a plausible-looking 1/0.
    T m_deltaLt, m_deltaLbar;

protected:
    /// Angle determination method option
    struct angmethod
    {
        enum type
        {
            Step = 0,
            Iteration  = 1,
            Predictor  = 2,
        };
    };

    /// \c BorderedMode option: which corrector chart(s) \a step() uses, and in which order --
    /// see the class documentation and \a step().
    struct borderedmode
    {
        enum type
        {
            Off      = 0,
            Fallback = 1,
            Primary  = 2,
        };
    };

};


} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsALMCrisfield.hpp)
#endif
