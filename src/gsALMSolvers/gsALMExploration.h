 /** @file gsALMExploration.h

    @brief Landscape-explorer orchestrator (thesis Alg. 8.1 ConstructLandscape,
    multiplicity-1 only) built on top of a gsALMBase arc-length solver.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#pragma once

#include <algorithm>
#include <functional>
#include <queue>
#include <vector>
#include <utility>

#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMBase.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMLandscape.h>

namespace gismo
{

/**
    @brief Explores the solution landscape of a nonlinear equilibrium problem.

    Starting from one converged seed solution, gsALMExploration traces
    equilibrium curves with a caller-supplied gsALMBase solver, detects singular
    points, classifies turning (limit) versus branch points, queues branch-switch
    jobs (one per \f$\pm\f$ critical mode), de-duplicates retraced curves, and
    assembles everything into a gsALMLandscape.

    ONE JOB PRODUCES ONE CURVE, TRACED IN BOTH ARC-LENGTH DIRECTIONS. On the
    TANGENT path (the common case: \a gsALMBase::computeBranchTangent
    resolved the pde2path swibra ABE emanating tangent \f$\tau_1\f$) each job's
    two directions are literally \f$\pm\texttt{sign}\cdot\tau_1\f$ ALONG the
    emanating branch, for a FIXED sign carried by that job -- the side is no
    longer ambiguous, only the two symmetric halves of a (generically pitchfork)
    bifurcation remain. A SINGLE
    job whose two sweeps are \f$+\tau_1\f$/\f$-\tau_1\f$ (arc-length direction
    tied to the tangent sign) is only correct when the corrector's predictor
    actually reads the installed secant beyond the first step (Riks,
    (Consistent)Crisfield); \c gsALMLoadControl's predictor does NOT -- it
    marches \f$\Lambda\f$ by the raw sign of the arc length for every subsequent
    step -- so at a (near-)symmetric bifurcation (\f$\dot\alpha_1\approx0\f$,
    BOTH mode signs require \f$\Lambda\f$ to increase) tying the sign to the arc-
    length direction makes only ONE sign reachable, MEASURED on the toy pitchfork
    fixture. The traceCurve/traceSweep caller therefore spawns one job PER SIGN
    (\c BranchPoints of them, mirroring the FALLBACK-NUDGE path's own loop and
    its documented driver default of 1 = "+ only"): each job's FORWARD sweep
    (arc length increasing) carries that job's fixed \f$\pm\tau_1\f$ as the
    predictor offset and is the "correct-side" sweep for its sign; each job's
    BACKWARD sweep carries the opposite offset with \f$\Lambda\f$ decreasing and
    is the "wrong-side" sweep, expected to fail or retrace -- see traceSweep and
    the branch-job creation block in traceCurve's caller. On the
    FALLBACK-NUDGE path (the ABE could not be resolved) the same reasoning
    still applies: a subcritical branch lives on the other side of
    \f$\Lambda^*\f$ and the side cannot be read off the critical mode alone (see
    traceCurve), so both directions must be swept. Either way EACH JOB's two
    directions are swept into a SINGLE curve object -- thesis Alg. 8.1 lines
    17--23: the second sweep starts from the same state with the opposite arc
    length (or opposite tangent sign), and the first sweep's points are reverted
    so the assembled point list runs from one far end through the start state to
    the other. Consequently the two halves of a seed, and the two halves of ONE
    JOB's branch, are no longer two curves that could be duplicates of each
    other -- but at \c BranchPoints \c >= \c 2 on the TANGENT path, the two SIGNS
    are now, once again, two separate jobs (hence two curves), because they are
    the only way to reach both physical halves when the corrector cannot.

    STORAGE ORDER OF A BOTH-DIRECTIONS CURVE, for consumers that do index
    arithmetic on gsALMLandscape<T>::Curve::points. The list reads FAR END ->
    START STATE -> FAR END, i.e. the FIRST sweep's block is stored REVERSED (it
    was traced outward from the start state, and traceCurve reverts it when the
    second sweep contributes anything). Storage order is therefore NOT trace
    order, the start state sits in the MIDDLE and not at index 0, and
    gsALMLandscape<T>::Curve::parentPointIdx of every child spawned inside that
    block is RENUMBERED by the reversal (i -> nFirst-1-i, applied before any
    child observes it). There is no public marker for the junction between the
    two sweeps: a consumer that needs one leg on its own must select it by a
    problem quantity (a monotone load, say), never by "the first half". A curve
    swept in ONE direction only -- a pristine rest-state seed, or a job whose
    second sweep contributed nothing -- keeps plain trace order.

    Points carry their own equilibrium PROVENANCE; the landscape is not uniformly
    certified. gsALMLandscape<T>::Point::equilibrium is true iff the row satisfies
    the residual test that was used to accept it: curve points from a converged
    corrector always have it true, and so does a refined singular point whose
    extended-system solve converged (\c gsStatus::Success) -- REGARDLESS of the
    solver's \c SingularPointComposite switch, which instead controls whether that
    termination test is the STRICTER composite one (\f$\|K_T V\|\f$ AND the
    equilibrium residuals TolF/TolU) or \f$\|K_T V\|\f$ alone; see the constructor
    note. A singular-point solve that did NOT converge stores its traced
    post-crossing point with \c equilibrium \c == \c false AND
    \c gsALMLandscape<T>::Point::unresolved \c == \c true, so it is never
    indistinguishable from a certified singular point in memory, in the CSV or in
    the HDF5 checkpoint. A consumer that needs the stronger composite certificate
    cannot read it off a per-row flag: that is a run-level property of the caller's
    own \c SingularPointComposite setting.

    On the FALLBACK-NUDGE path, of the two arc-length directions swept at a
    singular point, the one pointing away from the emanating branch falls back
    onto the parent branch and retraces it. Because both directions are swept
    into ONE curve, that fall-back is a SWEEP and not a curve: the C_start test
    (isRetrace, at \c RetraceTol) rewinds the offending sweep and leaves the
    genuine half that the other sweep found. Rewinding the offending sweep
    therefore does not cost the curve. The test runs at EVERY accepted step from
    \c StartSteps onward (not once), because a fall-back sweep can read above
    threshold at exactly \c StartSteps and only collapse onto the parent's own
    sampling a few steps later -- see retraceThreshold()'s implementation doc for
    the measured trajectory that motivated this. On the TANGENT path each spawned job's TWO
    swept directions are \f$+\f$ and \f$-\f$ its OWN fixed sign of \f$\tau_1\f$
    (see the ONE-JOB-PRODUCES-ONE-CURVE paragraph above for why the sign is now
    fixed per job rather than tied to the arc-length direction) -- for a
    secant-reading corrector (Riks, (Consistent)Crisfield) this fall-back IS
    expected to fire ~never, because the job's own FORWARD sweep already reaches
    the physical branch and its BACKWARD sweep is a short-lived excursion the
    corrector itself declines to continue past the singular point in the wrong
    load direction; for \c gsALMLoadControl it fires ROUTINELY instead, because
    that corrector's BACKWARD sweep marches
    \f$\Lambda\f$ away from the branch by construction (see the ONE-JOB
    paragraph) and is exactly the "wrong-side" sweep this net exists to remove,
    just as on the fallback path. The retrace safety net (\c m_retraceTests /
    \c m_retraceFired, reported by \a solve()) runs unconditionally on both
    paths and both correctors, as a net rather than as something either path or
    corrector is assumed to make unreachable.

    \c RetraceTol is a comparison scope of its own, deliberately looser than
    \c DedupTol and unrelated to \c MoveTol -- three tolerances, three jobs, and
    three defaults that are meant to be read together: \c DedupTol \c = \c 1e-4
    (branch-JOB de-duplication, isDuplicateJob), \c MoveTol \c = \c 1e-6 (the
    no-progress guard on an accepted step) and \c RetraceTol \c = \c 8e-2 (this
    test). \c MoveTol was split out of \c DedupTol, which used to serve both, and
    its default reproduces the old derived value \c DedupTol*1e-2 EXACTLY at the
    default \c DedupTol -- so introducing it changed no behaviour, while raising
    \c DedupTol no longer tightens the progress guard (which was MEASURED to
    delete points of an unrelated curve). isRetrace compares a SAMPLE of the
    sweep being traced against the stored SAMPLES of other curves, and two
    different samplings of one locus never coincide -- when SwitchLength and
    Length put parent and child on different grids the offset was MEASURED at
    2.99e-2 relative, i.e. far above any tolerance a job- or progress-comparison
    could use. It is NOT compared directly: the threshold is calibrated at the
    arc length \c StartSteps*SwitchLength the sweep has travelled from the branch
    point after \c StartSteps accepted steps -- retraceThreshold() scales
    \c RetraceTol by
    \f$\sqrt{\texttt{StartSteps}\cdot\texttt{SwitchLength}/\texttt{Length}}\f$ --
    but the TEST ITSELF now runs at every accepted step from
    \c StartSteps onward against that SAME fixed threshold, not once: NEAR the
    branch point, a genuine emanating branch departs from its parent monotonically
    (like the SQUARE ROOT of the arc length travelled, hence typically only
    growing past this calibration point OVER THE CONFIGURATIONS MEASURED -- see
    the caveat immediately below for when that stops holding),
    while a sweep that has fallen back onto the parent can clear the threshold at
    \c StartSteps and still collapse under it a few steps later, once it reaches
    the parent's own dense sampling.
    A single below-threshold reading no longer discards the sweep by
    itself past \c StartSteps: \c RetraceHits CONSECUTIVE below-threshold
    evaluations are required there, so that this widened window catches a
    persistent fall-back like the one just measured without also catching a
    genuine branch that merely crosses another stored curve transversally (a
    single, isolated below-threshold reading); a hit AT exactly \c StartSteps
    still discards on its own, unconditionally, as it always has. A problem whose genuinely
    distinct branches approach each other more closely than that needs
    \c RetraceTol lowered; a sweep that retraces on a coarser parent grid than
    this one needs it raised. Neither knob can recover a branch the tracer never
    reached: a nudge (the solver's \c Perturbation) that lets the corrector fall
    back onto the parent produces a sweep that IS a retrace, and it is removed as
    one -- a \a gsWarn says so whenever a branch curve disappears with a rewind
    among its causes: one message when EVERY sweep was rewound, a distinct one
    for the mixed case in which a sweep was rewound while the surviving sweep was
    kept but stored no point at all, and a third when a discarded sweep
    dropped more than \c StartSteps accepted points regardless of whether its
    sibling survives. A curve whose sweeps merely all failed is
    silent: no dedup decision was taken there.

    \c MaxPointsPerCurve is a per-SWEEP budget: a curve traced in both arc-length
    directions runs it once per direction and may hold up to twice that many
    accepted points -- plus any singular points stored on it, which the budget
    never counted. That keeps each leg's reach exactly what it was when the two
    directions were two curves.

    This is library code: it carries no shell / KLShell dependency. Everything
    problem-specific enters through the solver and an optional solution-constructor
    callback that maps a free-DOF vector to a deformed gsMultiPatch.

    The per-curve loop mirrors the proven de-risk driver
    (examples/example_ShearExploration.cpp): step + halve-arc-length-on-failure
    retry, perturbation-free stability detection, and extended-system singular
    point computation seeded from the LOCALIZED pre-crossing-side endpoint of the
    bisected crossing interval (see \c BisecMax / \c BisecLengthFloor and
    \a _localizeCrossing below) -- not, as before, from the pre-crossing ACCEPTED
    point itself.

    ARC-LENGTH CONTROL is halve-on-failure + snap-back-to-base-length by default.
    \c StepGrowth is an OPT-IN alternative (off by default, so the default control
    above is unchanged): an accepted step whose corrector converged within
    \c StepGrowthIter iterations doubles the arc length instead of snapping back,
    clamped in magnitude at the sweep's base length (pde2path \c sscontrol.m
    shape); see \c defaultOptions() and \c traceSweep().

    LOCALIZING A DETECTED CROSSING (\c BisecMax, \c BisecLengthFloor). A detected
    inertia flip is first BRACKETED between the two consecutive accepted points
    (pde2path \c bifdetec.m shape) and then bisected on the tangent inertia -- not
    on the residual or the indicator -- so the classification below runs on a
    tangent assembled close to the crossing rather than on the accepted arc-length
    grid, which may be much coarser. \c BisecMax caps the number of bisection
    probes; \c BisecLengthFloor is a tolerance RELATIVE to the sweep's base arc
    length (\c Length on a seed curve, \c SwitchLength on a branch curve) below
    which the bracket is considered resolved. Whichever bound is hit first stops
    the bisection; at the library defaults (\c Length \c = \c 1e-2, \c BisecMax \c
    = \c 10) that is \c BisecMax. Exhausting either bound is a NORMAL termination;
    only a bracket that no single probe ever refined is reported as a failed
    localization.

    RETRYING A FAILED LOCALIZATION (\c LocalizeRetries). A failed localization
    (no single bisection probe converged) is permanent by default: the crossing is
    marked unresolved on the landscape and abandoned. \a _localizeCrossing restarts
    its bracket from scratch on every call and is otherwise deterministic, so a
    repeated attempt at UNCHANGED knobs could never reach a different outcome.
    \c LocalizeRetries (default 1) instead re-runs \a _localizeCrossing with
    \c BisecMax and \c BisecLengthFloor escalated, COMPOUNDING across attempts --
    attempt N runs at \c BisecMax doubled N times and \c BisecLengthFloor halved N
    times relative to the caller-set values, saturating rather than
    overflowing/underflowing once the ladder runs far enough. Both knobs are
    restored to their caller-set values once, after the whole retry sequence ends
    (including when an exception propagates out of it) -- no caller ever observes
    an escalated value. Relaxing both is required: the two failure exits of
    \a _localizeCrossing are gated by different knobs (exhausting \c BisecMax, and
    retreating below \c BisecLengthFloor), so loosening only one leaves the other
    exit reachable. \c LocalizeRetries \c = \c 0 disables retrying, reproducing the
    pre-retry behaviour bit-for-bit. If every retry also fails, the honest
    unresolved marking stands unchanged; the failed bracket (lambda bounds and
    total probe count) is recorded on the landscape \c Point as a hook for a later
    resume, not consumed here.

    \tparam T coefficient type

    \ingroup gsALMSolvers
*/
template <class T>
class gsALMExploration
{
public:

    /// Constructs a free-DOF vector into a deformed geometry (optional; enables
    /// per-point gsMultiPatch storage in the landscape). Returns true on success.
    typedef std::function<bool(const gsVector<T> &, gsMultiPatch<T> &)> SolutionConstructor_t;

    /// @brief Why one arc-length sweep (traceSweep) stopped.
    ///
    /// The sweep loop has exactly four exits and this enumerates all four. The two
    /// Underflow* values are the same 1e-6-of-base-length arc-length floor reached
    /// from two different causes: a corrector that never converged, and a corrector
    /// that converged onto the point it started from.
    enum struct SweepTermination
    {
        PointBudget,          ///< the per-sweep MaxPointsPerCurve budget was exhausted
        UnderflowStepFailed,  ///< a non-Success step halved the arc length below the floor
        UnderflowNoProgress,  ///< a converged but non-moving step halved it below the floor
        RetraceRewound        ///< the C_start dedup rewound the sweep (traceSweep returned false)
    };

public:

    /**
     * @brief Constructor.
     *
     * @param solver caller-configured ALM method (non-owning pointer). The caller
     *               is responsible for the singular-point regime: bisection first
     *               stage on (SingularPointComputeTolB != 0), a problem-scaled
     *               Perturbation tau -- GOVERNING THE FALLBACK-NUDGE PATH ONLY
     *               (the common path now seeds from the ABE emanating
     *               tangent \f$\tau_1\f$ via \a gsALMBase::computeBranchTangent
     *               instead). When that resolution fails, \c switchBranch()'s
     *               nudge \f$U^*+V/\tau\f$ is used as before, so \f$1/\tau\f$ must
     *               still be comparable to the amplitude the emanating branch has
     *               one SwitchLength away from the branch point -- too small a
     *               nudge simply falls back onto the parent branch) and a
     *               TolF/MaxIter regime in which the extended system can reach the
     *               equilibrium manifold. \c BranchTangentPerturbation and
     *               \c BranchTangentTol (on the solver) govern the tangent path.
     *
     * @note The explorer modifies NO solver option. \c SingularPointComposite in
     *       particular is only READ, never set: it is the caller's switch and the
     *       explorer merely records its value on the singular point it publishes.
     *       Forcing it on was measured to drive the extended Newton -- an inexact
     *       Newton on a singular tangent -- past the point where it stops improving
     *       and into divergence, so the decision is deliberately left with the caller.
     *
     *       CONSEQUENCE, and it is actionable. A converged refined singular point is
     *       ALWAYS stored with gsALMLandscape<T>::Point::equilibrium \c == \c true --
     *       it met its own termination test either way. If you want that test to be
     *       the STRICTER one (the equilibrium residuals TolF/TolU, not just
     *       \f$\|K_T V\|\to 0\f$), you must BOTH
     *       \c options().setSwitch("SingularPointComposite",true) AND
     *       \c applyOptions() before handing the solver over -- the solve branches on
     *       the applied value, while the explorer reads the option list, and there is
     *       no public accessor for the applied state. Otherwise the point is still
     *       stored (equilibrium = true), but a \a gsWarn says its termination test did
     *       NOT include the equilibrium residuals, so a downstream consumer that
     *       specifically needs THAT stronger certificate cannot get it from the stored
     *       point and must track its own \c SingularPointComposite setting. \c TolF
     *       stays the caller's either way, so how strong the certificate is remains
     *       configurable.
     */
    explicit gsALMExploration(gsALMBase<T> * solver);

    /// Access the options (see defaultOptions()).
    gsOptionList & options() { return m_options; }

    /// Sets the optional solution-constructor callback.
    void setSolutionConstructor(const SolutionConstructor_t & fun) { m_solutionConstructor = fun; }

    /**
     * @brief Explores the landscape (Alg. 8.1).
     *
     * Seeds the job queue with (U0,L0) -- which must be a converged equilibrium;
     * use U0=0, L0=0 for an undeformed rest state -- and explores until the queue
     * is empty or MaxCurves is reached.
     *
     * @return Success unless the FIRST curve cannot take a single step (empty
     *         landscape). Individual curve failures are logged and tolerated.
     */
    gsStatus solve(const gsVector<T> & U0, T L0);

    /**
     * @brief Explores the landscape from several seed states (Alg. 8.1).
     *
     * Grows ONE landscape -- with shared curve/job de-duplication and a single
     * branch-job queue -- from multiple seeds. Each seed must be a converged
     * equilibrium; use U0=0, L0=0 for an undeformed rest state. Seed jobs are
     * queued in the given order, forward before backward per seed (the backward
     * job is skipped at a pristine rest state). Explores until the queue is empty
     * or MaxCurves is reached.
     *
     * @param seeds converged (U0,L0) equilibria to seed the exploration from.
     *
     * @return Success unless the FIRST curve cannot take a single step (empty
     *         landscape). Individual curve failures are logged and tolerated.
     */
    gsStatus solve(const std::vector<std::pair<gsVector<T>,T> > & seeds);

    /// Const access to the assembled landscape.
    const gsALMLandscape<T> & landscape() const { return m_landscape; }

    /// @brief Termination reason of every arc-length sweep executed by the last
    /// solve() call, paired with the curve index the sweep was tracing, in
    /// execution order (a bothDirections curve contributes TWO entries).
    ///
    /// Cleared at the start of every solve(). NOTE the curve index is the index at
    /// the time of the sweep: a curve every one of whose sweeps produced no point
    /// is removed by traceCurve (gsALMLandscape<T>::removeLastCurve), so an entry
    /// may name a curve that is no longer in the landscape.
    const std::vector<std::pair<index_t,SweepTermination> > & sweepTerminations() const
    { return m_sweepTermination; }

protected:

    /// Sets the default options.
    void defaultOptions();

    /// @brief One unit of work: trace one equilibrium curve from a seed state.
    ///
    /// A job produces exactly ONE curve. When \a bothDirections is true the curve
    /// is traced in BOTH arc-length directions from the same start state and the
    /// two sweeps are assembled into a single ordered point list (thesis Alg. 8.1
    /// lines 17--23); see traceCurve.
    struct Job
    {
        gsVector<T> U;              ///< converged parent state (displacement)
        T           L;              ///< converged parent state (load)
        gsVector<T> Unudged;        ///< nudged start (FALLBACK branch jobs); empty for seed
                                     ///< jobs and for TANGENT-path branch jobs
        gsVector<T> tangentU;       ///< emanating-branch tangent (pde2path swibra
                                     ///< ABE), displacement part; EMPTY on the fallback-nudge
                                     ///< path and on seed jobs
        T           tangentL;       ///< emanating-branch tangent, load part
        bool        backward;       ///< arc-length direction of the FIRST sweep
        bool        bothDirections; ///< also sweep the opposite direction into the SAME curve
        index_t     parentCurve;    ///< parent curve index, -1 for seed jobs
        index_t     parentPointIdx; ///< branch-point index in the parent curve, -1 for seed jobs

        /// Precedence rule: \c tangentU non-empty => TANGENT path (Euler
        /// predictor step along +-tau1, see traceSweep); else \c Unudged non-empty =>
        /// FALLBACK-NUDGE path (\f$U^*+V/\tau\f$); else a seed job.
    };

    /// @brief Traces a single curve for \a job, appending points to the landscape
    /// and pushing new branch jobs onto \a queue. Returns true if at least one
    /// point was accepted (a non-empty curve was produced).
    bool traceCurve(const Job & job, std::queue<Job> & queue);

    /// @brief ONE arc-length sweep of curve \a cid, from \a job's start state in
    /// direction \a backward, appending accepted points to that curve.
    ///
    /// Branch jobs discovered by the sweep are NOT pushed onto the queue directly:
    /// they are collected in \a pending together with the index (in \a cid) of the
    /// branch point they emanate from, in \a pendingPt. traceCurve may revert the
    /// first sweep's point block when assembling the two sweeps, which renumbers
    /// those indices; deferring the push lets it remap them before any child can
    /// observe a stale parent index. \c m_jobKeys is still updated immediately, so
    /// isDuplicateJob sees exactly the jobs it saw before.
    ///
    /// \returns false if the sweep was rewound as a retrace (its own points and
    /// pending jobs removed); true otherwise. Points contributed by an EARLIER
    /// sweep of the same curve are never touched.
    bool traceSweep(const Job & job, index_t cid, bool backward,
                    std::vector<Job> & pending, std::vector<index_t> & pendingPt);

    /// @brief Localizes a detected inertia crossing by bisecting the arc-length
    /// interval between two consecutive ACCEPTED points on the tangent INERTIA
    /// (pde2path bifdetec shape), and returns the pre-crossing-side endpoint of the
    /// final bracket in (\a Uloc,\a Lloc).
    ///
    /// \param Uold,Lold,negOld  previous accepted point and its inertia (bracket low end)
    /// \param Ucur,Lcur,negCur  just-accepted point and its inertia (bracket high end);
    ///                          negCur != negOld is the caller's detection condition
    /// \param dLb   SIGNED arc length of the accepted step that produced (Ucur,Lcur)
    /// \param dLb0  the sweep's base (unbisected) signed arc length; scales BisecLengthFloor
    /// \param cid   curve index, for Verbose reporting only
    /// \param[out] Uloc,Lloc  the localized point (pre-crossing side of the bracket)
    /// \param[out] bracket    the final bracket WIDTH in arc length (always set)
    /// \param[out] probesOut  if non-null, receives the number of bisection probes
    ///                        spent by THIS call, set on every return path; \c nullptr
    ///                        (default) is accepted and simply skips the write.
    /// \returns false ONLY when not a single bisection probe converged, i.e. the interval
    ///          was never refined; exhausting BisecMax or reaching BisecLengthFloor are
    ///          NORMAL terminations and return true.
    bool _localizeCrossing(const gsVector<T> & Uold, T Lold, index_t negOld,
                           const gsVector<T> & Ucur, T Lcur, index_t negCur,
                           T dLb, T dLb0, index_t cid,
                           gsVector<T> & Uloc, T & Lloc, T & bracket,
                           index_t * probesOut = nullptr);

    /// @brief Branch-job dedup: true if a job with load \a Lstar and signed
    /// normalized nudge direction \a dir already exists (queued or executed).
    ///
    /// The arc-length direction is NOT part of the key: since both directions are
    /// traced into one curve by a single job (see traceCurve), one (Lstar,dir) pair
    /// corresponds to exactly one job -- there is no same-key/opposite-direction pair
    /// left to distinguish.
    bool isDuplicateJob(T Lstar, const gsVector<T> & dir) const;

    /// @brief C_start dedup: true if (U,L) coincides with a stored point of a
    /// curve OTHER than \a selfCurve (simplified thesis eq. 8.1), at the
    /// scale-aware tolerance retraceThreshold(). Delegates to retraceDistance():
    /// \c isRetrace(U,L,selfCurve) \c == \c retraceDistance(...) \c <
    /// \c retraceThreshold(), bit-identical to the predicate this replaces --
    /// \c max(ratioU,ratioL) \c < \c tol is the same test as the old
    /// \c sameU \c && \c sameL, and a MINIMUM ratio below tol exists iff some
    /// point's ratio is below tol.
    bool isRetrace(const gsVector<T> & U, T L, index_t selfCurve) const;

    /// @brief Diagnostic: the isRetrace predicate's own distance, exposed
    /// so callers (and traceSweep's Verbose print) can see the number instead of
    /// only the boolean. Returns the MINIMUM over all OTHER curves' stored points
    /// \a pt of \f$\max(\|U-U_p\|/\max(1,\|U_p\|),\,|L-L_p|/\max(1,|L_p|))\f$, and
    /// reports which curve/point attained it in \a matchCurve / \a matchPoint
    /// (-1/-1 when no comparable point exists -- e.g. the very first branch curve
    /// -- in which case a large sentinel is returned so the caller need not
    /// special-case it). \a selfCurve is excluded exactly as in isRetrace().
    ///
    /// Complexity: O(sum of stored points over all curves != selfCurve), same as
    /// isRetrace() before it (this does not add a pass: isRetrace() now calls this
    /// once instead of running its own loop).
    T retraceDistance(const gsVector<T> & U, T L, index_t selfCurve,
                      index_t & matchCurve, index_t & matchPoint) const;

    /// @brief The distance below which isRetrace calls a point a retrace:
    /// \c RetraceTol scaled by \f$\sqrt{\texttt{StartSteps}\cdot
    /// \texttt{SwitchLength}/\texttt{Length}}\f$.
    ///
    /// This VALUE is calibrated at the arc length \c StartSteps*SwitchLength
    /// travelled from the branch point, and a GENUINE emanating branch is by then
    /// only \f$O(\sqrt{\cdot})\f$ away from its parent. The caller (the
    /// C_start block in traceSweep) now tests every accepted step from
    /// \c StartSteps onward against this SAME fixed value, not a single step --
    /// NEAR this calibration point, over the configurations measured, a genuine
    /// branch's distance typically only grows, so widening the test window costs
    /// it little (see the caveat in the class doxygen's \c RetraceTol
    /// comparison-scope paragraph for the shape -- distinct branches that run
    /// closer together than this -- where that stops holding),
    /// while a sweep that falls back
    /// onto the parent can (MEASURED) clear the threshold once and collapse under
    /// it a few steps later; \c RetraceHits CONSECUTIVE below-threshold
    /// evaluations past \c StartSteps are what the caller now requires before
    /// discarding on that basis, precisely so that "costs it little" does not have
    /// to mean "costs it nothing": a single below-threshold reading, on its own,
    /// no longer discards a sweep past \c StartSteps. A FIXED tolerance is therefore a different criterion at every
    /// \c SwitchLength -- MEASURED to delete a genuine branch at
    /// \c SwitchLength <= \c Length/40. In the scaled variable the retrace and
    /// the genuine populations are 2.3x away from the default on both sides
    /// over the MEASURED range \c SwitchLength in [5e-4, 5e-3] at
    /// \c Length \c = \c 2e-2; the calibration table is in the implementation.
    ///
    /// WORKED VALUES, so that "scaled" is a number and not an adjective:
    /// * calibration configuration (\c StartSteps \c 3, \c SwitchLength \c 5e-3,
    ///   \c Length \c 2e-2): scale \f$\sqrt{3\cdot 5\!\cdot\!10^{-3}/2\!\cdot\!10^{-2}}
    ///   = 0.866\f$, threshold \f$8\cdot 10^{-2}\cdot 0.866 = 6.93\cdot 10^{-2}\f$
    ///   -- within 1% of the fixed 7e-2 the criterion replaced, which is why that
    ///   landscape did not move;
    /// * THIS CLASS'S OWN DEFAULTS (\c Length \c 1e-2 AND \c SwitchLength \c 1e-2,
    ///   \c StartSteps \c 3): the ratio is 3, so the threshold is
    ///   \f$8\cdot 10^{-2}\sqrt{3} = 0.1386\f$ -- 1.73x the bare \c RetraceTol and
    ///   1.98x the fixed 7e-2 it replaced, i.e. at library defaults this test is
    ///   markedly LOOSER, not tighter. No benchmark of this library exercises the
    ///   default configuration (every driver sets both lengths), so that regime is
    ///   covered only indirectly, by a probe at the same ratio 3
    ///   (\c SwitchLength \c 0.02 against \c Length \c 2e-2), where the genuine
    ///   population still cleared the threshold by 2.2x.
    T retraceThreshold() const;

protected:

    /// Non-owning ALM solver.
    gsALMBase<T> * m_solver;

    /// Exploration options.
    gsOptionList m_options;

    /// Optional deformed-geometry constructor.
    SolutionConstructor_t m_solutionConstructor;

    /// The assembled landscape.
    gsALMLandscape<T> m_landscape;

    /// @brief Dedup key of a created branch job.
    struct JobKey
    {
        T           L;          ///< load at the singular point
        gsVector<T> dir;        ///< signed normalized nudge direction
    };

    /// Dedup keys of created branch jobs.
    std::vector<JobKey> m_jobKeys;

    /// isRetrace safety-net telemetry: evaluations attempted / times the
    /// predicate FIRED (i.e. resulted in \c retrace \c == \c true, the two-clause
    /// decision -- see the C_start block in traceSweep). Reset per solve() call; printed
    /// unconditionally under Verbose at the end of solve(). The DENOMINATOR is
    /// no longer one evaluation per sweep -- the C_start block now evaluates at every
    /// accepted step from StartSteps onward, so a branch sweep contributes up to
    /// MaxPointsPerCurve-StartSteps+1 evaluations instead of exactly one; the printed
    /// "fired N of M evaluations" rate therefore moves at every configuration that traces
    /// a branch curve (reported, not pinned by any test -- see the implementation's
    /// grounding note). The NUMERATOR counts a predicate
    /// FIRING at most ONCE per sweep -- a sweep the spawnedJobs guard keeps (see the
    /// C_start block) would otherwise re-fire, hence re-count, at every remaining step
    /// once \c retrace latches true, inflating the printed rate; the guard against that
    /// is the same \c spawnedJobsNoticed latch that silences the repeated Verbose notice.
    /// Expected to fire ~never on the TANGENT path for a
    /// secant-reading corrector (Riks, (Consistent)Crisfield); fires ROUTINELY there
    /// too for gsALMLoadControl (round-1 respec: each spawned job's own BACKWARD
    /// sweep is by construction the "wrong-side" one for that corrector -- see the
    /// class doxygen). Meaningful on the FALLBACK-NUDGE path regardless of corrector.
    index_t m_retraceTests, m_retraceFired;

    /// Per-sweep termination reasons of the current solve(), see sweepTerminations().
    std::vector<std::pair<index_t,SweepTermination> > m_sweepTermination;

}; // class gsALMExploration


} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsALMExploration.hpp)
#endif
