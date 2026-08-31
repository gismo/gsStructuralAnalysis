 /** @file gsALMExploration.hpp

    @brief Implementation of gsALMExploration.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#pragma once

#include <limits>

#include <gsCore/gsMemory.h>

namespace gismo
{

template <class T>
gsALMExploration<T>::gsALMExploration(gsALMBase<T> * solver)
:
m_solver(solver), m_retraceTests(0), m_retraceFired(0)
{
    GISMO_ENSURE(m_solver != nullptr, "gsALMExploration: solver pointer is null.");
    defaultOptions();
}

template <class T>
void gsALMExploration<T>::defaultOptions()
{
    m_options.addInt   ("MaxCurves","Hard cap on landscape curves",8);
    // PER-SWEEP budget, not per-curve: a curve traced in both arc-length
    // directions runs the budget once per direction and may hold up to
    // 2*MaxPointsPerCurve ACCEPTED points, plus the singular points stored on it,
    // which the budget never counted (MEASURED on the calibration benchmark: a
    // merged curve of 31 points = 10 + 20 accepted + 1 stored singular point
    // against a cap of 20, and 11 = 5 + 5 + 1 against a cap of 5). The budget per
    // DIRECTION -- i.e. how far each leg of a locus reaches -- is what the option
    // controlled before the two directions shared one curve object, and it is
    // unchanged; a shared budget would have halved each leg's reach.
    m_options.addInt   ("MaxPointsPerCurve","Maximum accepted points per arc-length SWEEP; a curve traced in both directions may hold up to twice this",50);
    m_options.addReal  ("Length","Per-curve initial arc length",1e-2);
    m_options.addReal  ("SwitchLength","Arc length for the first steps of a branch curve",1e-2);
    // The comparison against existing curves now runs at EVERY accepted step
    // from StartSteps onward (was: exactly once, at StartSteps), so this is a FLOOR on
    // when scrutiny starts, not a single evaluation point. Consequence for the two
    // degenerate settings this inverts nobody's existing usage of (no driver or
    // unittest in this module sets StartSteps outside {3}), but is worth stating: with
    // the old `==` predicate, StartSteps <= 0 disabled the test entirely (stepsTaken is
    // >= 1 by the time it is checked, so `stepsTaken == 0` never matched) -- a usable
    // "turn the dedup off" idiom. With `>=` that is no longer true: StartSteps <= 0 now
    // makes the test fire from the FIRST accepted step, i.e. maximally aggressive, the
    // opposite of disabled. A caller that wants the dedup off must set RetraceTol to an
    // unreachable value (e.g. 0) instead.
    // With the two-clause predicate below (clause A at exactly StartSteps,
    // clause B = RetraceHits consecutive hits later), StartSteps <= 0 has a second
    // consequence: clause A (`stepsTaken == startSteps`) becomes UNREACHABLE, because
    // this block only ever runs once stepsTaken >= 1, so `stepsTaken == 0` (or any
    // negative StartSteps) never matches. The predicate then degenerates to clause B
    // alone -- pure persistence, from the first accepted step -- which is a different,
    // not merely earlier, criterion than the StartSteps >= 1 configurations this module
    // is calibrated on.
    m_options.addInt   ("StartSteps","C_start dedup: FIRST accepted step (of a branch sweep) at which comparison against existing curves begins; every accepted step from here onward is tested, not a single evaluation point. NOTE: <= 0 now means \"test from step 1\", the opposite of \"disabled\" under the predicate this replaced -- see the comment above.",3);
    m_options.addReal  ("DedupTol","Duplicate tolerance for branch-JOB de-duplication (isDuplicateJob)",1e-4);
    // Part B: the no-progress guard used to derive its tolerance as
    // DedupTol*1e-2, i.e. one number served two unrelated semantic jobs -- a
    // comparison scope and a step-length/progress control. Raising the dedup
    // tolerance therefore also tightened the progress guard, which was MEASURED
    // (63-report §3) to delete the two fold-bracketing rows of an unrelated curve.
    // The default below is EXACTLY the old effective value (DedupTol 1e-4 * 1e-2),
    // so introducing it changes nothing; the two are now independent.
    m_options.addReal  ("MoveTol","No-progress guard: relative state change below which an accepted step is rejected as a stall",1e-6);
    // Third scope, third tolerance (PyNCT carries nine, and derives its step
    // control from none of them). isRetrace compares a SAMPLE of the sweep being
    // traced against the stored SAMPLES of other curves, so it must tolerate the
    // offset between two different samplings of the same locus -- it can never be
    // as tight as a job- or progress-comparison.
    //
    // IT IS NOT COMPARED DIRECTLY: retraceThreshold() scales it by
    // sqrt(StartSteps*SwitchLength/Length), and the derivation lives there. In
    // short: the distance the test sees for a GENUINE emanating branch scales as
    // the SQUARE ROOT of the arc length the sweep has travelled AT StartSteps
    // accepted steps -- that is the arc length the calibration below is anchored
    // to, and it stays fixed even though the predicate that compares
    // against it now runs at every accepted step from StartSteps onward, not once
    // -- so a FIXED tolerance is a different criterion at every SwitchLength
    // -- it was MEASURED to delete a genuine branch at SwitchLength <= Length/40.
    // In the scaled variable the two populations are 2.3x apart on both sides at
    // EVERY SwitchLength, and the default below is the geometric centre of the
    // MEASURED window (3.450e-2, 1.851e-1); see retraceThreshold().
    // At the calibration configuration (StartSteps 3, SwitchLength 5e-3,
    // Length 2e-2) the effective threshold is 8e-2*0.866 = 6.93e-2, i.e. within
    // 1% of the fixed 7e-2 this option shipped with, which is why the calibration
    // landscape is unchanged. It is calibrated on ONE benchmark: a problem whose
    // distinct branches approach each other more closely than this needs it
    // lowered; a sweep that retraces on a coarser parent grid needs it raised.
    m_options.addReal  ("RetraceTol","C_start dedup: relative distance -- scaled by sqrt(StartSteps*SwitchLength/Length), see retraceThreshold() -- below which a traced point counts as lying on an already-stored curve",8e-2);
    // RetraceTol's calibration is (by design, see the comment above and
    // retraceThreshold()) only 1.2% apart between the genuine and the phantom sweep AT
    // StartSteps -- MEASURED on the modified-Bratu A-C crossing (factors 2.318 vs
    // 2.346). Retuning the threshold to separate them is arithmetically excluded: any
    // threshold tight enough to catch the phantom's step-3 reading also deletes the
    // genuine curve's own 2.318. The signal that DOES separate the two populations is
    // PERSISTENCE: the phantom's kept sweep stays below threshold at every step from
    // its first hit (MEASURED: step 8) through the end of the sweep, while a genuine
    // sweep that merely crosses another stored curve transversally produces a SINGLE
    // below-threshold reading before resuming its climb (the discriminator the earlier
    // measurement identified but did not act on). RetraceHits is the number of
    // CONSECUTIVE below-threshold evaluations of ACCEPTED steps -- see the C_start
    // block's own comment for exactly what "consecutive" means -- required to discard a
    // sweep once PAST StartSteps; a hit AT exactly StartSteps keeps its own,
    // unconditional single-evaluation semantics (clause A) regardless of this option.
    // Three regimes from one knob:
    //   <= 0        : clause B disabled -> the predicate degenerates to clause A alone,
    //                 i.e. a single evaluation at StartSteps.
    //   1           : EXACTLY the original one-hit predicate (one hit anywhere
    //                 >= StartSteps discards the sweep).
    //   2 (default) : the shipped predicate -- a hit at StartSteps discards on its own
    //                 (clause A), or RetraceHits consecutive hits anywhere later in the
    //                 sweep discard it (clause B).
    // This option alone cannot disable clause A; a caller wanting C_start off entirely
    // must still set RetraceTol to an unreachable value (e.g. 0), same idiom as before.
    // Default 2, not 3: a transversal crossing (the case this knob exists to survive)
    // is, by construction, a SINGLE below-threshold hit, which RetraceHits = 2 already
    // rejects; no oracle in this tree can measure a benefit of 3 over 2 (none has a
    // transversal-crossing fixture -- that fixture remains to be built), while the
    // COST of a larger RetraceHits is measurable and one-sided: every extra step before
    // the discard is one more step in which singular-point detection (which runs AFTER
    // this block) can queue a branch job and thereby SUPPRESS the discard via the
    // spawnedJobs guard below. Do not ship RetraceHits > 3.
    m_options.addInt   ("RetraceHits","C_start dedup: number of CONSECUTIVE below-threshold accepted-step evaluations, past StartSteps, required to discard a sweep (clause B); a hit AT exactly StartSteps still discards unconditionally (clause A). <= 0 disables clause B (pre-01 predicate: clause A only); 1 reproduces the original one-hit predicate (one hit anywhere discards); 2 (default) is the shipped persistence-aware predicate. Never > 3.",2);
    // Number of critical-mode/tangent SIGNS explored at a singular point
    // (1 = +mode/+tau1 only, 2 = +/- both). Each sign spawns exactly ONE job, and
    // that job sweeps BOTH arc-length directions into a single curve (see
    // traceCurve for why a single direction is not enough), so a singular point
    // creates at most BranchPoints jobs -- not 2*BranchPoints, as it did while the
    // two directions were queued as two independent jobs.
    //
    // Round 0 governed the FALLBACK-NUDGE path only, spawning exactly
    // ONE job on the TANGENT path regardless of this option. That was MEASURED wrong for gsALMLoadControl --
    // its predictor does not read the installed secant beyond the first step, so
    // tying +tau1/-tau1 to the sweep's arc-length direction (as a single job's two
    // sweeps do) makes only ONE sign of a (near-)symmetric bifurcation reachable
    // (see the class doxygen's ONE-JOB-PRODUCES-ONE-CURVE paragraph). BranchPoints
    // now governs BOTH paths identically: on the TANGENT path it selects how many
    // SIGNED tau1 jobs are spawned (mirroring the fallback loop exactly), and at
    // the driver default of 1 this reproduces the round-0 single-job behaviour
    // bit-for-bit (MEASURED on example_ModifiedBratuExploration).
    m_options.addInt   ("BranchPoints","Critical-mode/tangent signs explored per singular point (+/- mode or +/- tau1 = 2); each spawns ONE job, swept in both arc-length directions. Governs BOTH the FALLBACK-NUDGE and the TANGENT path (round 1) identically.",2);
    // tau1 is a DIRECTION -- the ABE's genuine
    // contribution is WHICH WAY the emanating branch leaves, not HOW FAR to step. The
    // first attempt at this seeding scaled the Euler predictor offset by SwitchLength
    // (the arc-length base), which is wrong by a POWER: at a pitchfork the emanating
    // branch leaves the singular point with lambdadot=0 and grows like C*sqrt(s), so the
    // offset needed to clear the trivial branch's Newton basin scales like O(sqrt(ds)),
    // not O(ds). MEASURED consequence: at SwitchLength=5e-3 (~200x smaller than the
    // tuned fallback offset 1/Perturbation=1.0 that example_ModifiedBratuExploration's
    // own comment states is required), BOTH +-tau1 sweeps retraced and curve B was
    // dropped entirely -- strictly worse than the fallback nudge. Sweeping SwitchLength
    // up does NOT recover it (it collapses the corrector into arc-length halving
    // instead): this is a scaling error, not a tuning problem.
    //
    // BranchEscape decouples MAGNITUDE from DIRECTION. -1 (default) resolves to the
    // UNFLOORED fallback-nudge scale 1/Perturbation, read from the solver's own options
    // -- the scale that is already tuned and proven on these drivers, see
    // switchBranch()/the Unudged construction; any positive value is used directly as
    // the offset magnitude. tau1 is re-normalised to unit Euclidean length before being
    // scaled by this magnitude (see traceSweep), so the escape distance no longer
    // depends on SwitchLength or on whatever metric computeBranchTangent used to
    // normalise tau1U/tau1L internally.
    //
    // The FALLBACK-NUDGE path's OWN magnitude may now be floored above
    // 1/Perturbation by RetraceBallFloor (see its doc immediately below) -- the two
    // paths can therefore DISAGREE in magnitude even though both read the same -1
    // sentinel semantics. BranchEscape = -1 always reuses the RAW, unfloored
    // 1/Perturbation; it does NOT pick up the fallback path's floor. A caller wanting
    // the TANGENT path to match a floored fallback magnitude must set BranchEscape
    // explicitly.
    m_options.addReal  ("BranchEscape","Euler-predictor escape magnitude for the TANGENT path: -1 (default) reuses the UNFLOORED fallback-nudge scale 1/Perturbation (see RetraceBallFloor, which floors the fallback path's OWN magnitude but is NOT read here); any positive value is used directly. Inert on the FALLBACK-NUDGE path.",-1);
    // Floors the FALLBACK-NUDGE branch-job's absolute magnitude 1/Perturbation
    // at margin*retraceThreshold()*max(1,||Ustar||) -- the SAME normalization
    // isRetrace's own ball uses (see retraceDistance()) -- so a caller-tuned
    // 1/Perturbation that clears the retrace ball at a LOW-state-norm locus (where it
    // was calibrated) does not silently fall short at a HIGHER-state-norm locus, where
    // the ball radius has grown but a fixed absolute nudge has not. <= 0 (default)
    // DISABLES the floor entirely: the fallback path is then bit-identical to
    // the unfloored nudge, and a caller-set Perturbation is never overridden. When
    // active (> 0), this OVERRIDES a caller-set Perturbation on the fallback path at
    // whichever locus the floor binds -- it is a floor, not a replacement: any locus
    // whose existing 1/Perturbation nudge already clears margin*ball is left
    // bit-identical (max() picks 1/Perturbation). There is no library-derived default:
    // the margin is a driver-specific calibration (an equality at ONE locus of
    // example_ModifiedBratuExploration, see that driver's own setReal call), so this
    // ships disabled and each driver that needs it derives and sets its own value.
    m_options.addReal  ("RetraceBallFloor","Floor on the FALLBACK-NUDGE branch-job's absolute magnitude 1/Perturbation, as a multiple (margin) of isRetrace's own ball radius retraceThreshold()*max(1,||Ustar||); <= 0 (default) disables the floor and leaves the fallback path, and any caller-set Perturbation, untouched. No library-derived default exists -- each driver derives and sets its own margin.",0);
    m_options.addInt ("BisecMax","Singular-point localization: maximum bisection probes on the "
                                 "detection interval (pde2path bisecmax)",10);
    m_options.addReal("BisecLengthFloor","Singular-point localization: bisection stops when the "
                                 "bracket width drops below this FRACTION of the sweep's base arc "
                                 "length (Length on seed curves, SwitchLength on branch curves)",1e-6);
    m_options.addInt ("LocalizeRetries","Singular-point localization: additional "
                                 "_localizeCrossing attempts after the first one fails to refine "
                                 "the interval; each retry doubles BisecMax and halves "
                                 "BisecLengthFloor for that call only. 0 disables retrying.",1);
    m_options.addSwitch("Verbose","Verbose output",false);
    m_options.addString("OutputPrefix","When non-empty, writeCsv(prefix+\".csv\") after every completed curve","");
    // Iteration-count step GROWTH (pde2path sscontrol.m:37-43, "very good step ⇒
    // increase ds by dsincfac"). OFF by default: the shipped control is
    // halve-on-failure + snap back to the base length, which keeps the sampling
    // density of a run reproducible and every CSV oracle untouched. With the switch
    // ON the snap-back is replaced by a x2 ladder, clamped in MAGNITUDE at the base
    // length (Length for a seed curve, SwitchLength for a branch curve), so the arc
    // length can never EXCEED what the default control uses -- it only approaches it
    // more gradually after a halving.
    m_options.addSwitch("StepGrowth","Opt-in iteration-count arc-length growth after an accepted step (pde2path sscontrol); off = halve-and-reset, the shipped control",false);
    // pde2path: dsinciter = imax/2 (stanparam.m:50), i.e. HALF the corrector's own
    // iteration cap. -1 reproduces that literally against the solver's MaxIter
    // (gsALMBase.hpp:26, default 100); a positive value overrides it, which matters
    // because MaxIter is a generous CAP here, not a target (contrast the solver-side
    // AdaptiveIterations, default 10).
    m_options.addInt   ("StepGrowthIter","Corrector-iteration threshold for StepGrowth: a step converged in <= this many iterations doubles the arc length. -1 = MaxIter/2 of the solver (pde2path dsinciter)",-1);
}

template <class T>
bool gsALMExploration<T>::isDuplicateJob(T Lstar, const gsVector<T> & dir) const
{
    const T tol = m_options.getReal("DedupTol");
    for (size_t i = 0; i != m_jobKeys.size(); ++i)
    {
        const JobKey & k = m_jobKeys[i];
        const bool sameLoad = math::abs(Lstar - k.L) < tol * math::max( (T)1, math::abs(Lstar) );
        const bool sameDir  = (dir.size() == k.dir.size()) && ( (dir - k.dir).norm() < tol );
        if (sameLoad && sameDir)
            return true;
    }
    return false;
}

template <class T>
T gsALMExploration<T>::retraceThreshold() const
{
    // WHY THE THRESHOLD IS NOT A FIXED NUMBER.
    //
    // isRetrace is now evaluated at EVERY accepted step from
    // StartSteps onward (traceSweep's C_start block, `stepsTaken >= startSteps`),
    // not once. This function's own CALIBRATION, below, is still anchored at
    // exactly StartSteps accepted steps -- the threshold magnitude does not grow
    // with stepsTaken, only the CALLER's evaluation cadence changed. That anchoring
    // is deliberately kept: the sweep has then travelled  s = StartSteps *
    // SwitchLength  in arc length from the branch point. NEAR the branch point, at a
    // bifurcation, the emanating branch leaves the parent like the SQUARE ROOT of that
    // distance -- so, OVER THE CONFIGURATIONS MEASURED (this class's calibration
    // benchmark), a GENUINE sweep's distance from its parent typically GROWS past this
    // point (see the caveat in gsALMExploration.h's class doxygen, RetraceTol
    // comparison-scope paragraph: a problem whose genuinely distinct branches run
    // closer together than that needs RetraceTol lowered, and that is exactly the
    // shape this is a LOCAL asymptotic about, not a general property). A sweep that has fallen back onto the parent can (and, MEASURED on the
    // modified-Bratu benchmark's A-C crossing, does) drop back under the very same
    // threshold a few steps after clearing it once -- see the C_start block's own
    // comment for the measured trajectory and for RetraceHits, the k-consecutive-hits
    // predicate that distinguishes that persistent fall-back from a genuine
    // branch's single, transversal below-threshold reading. Re-deriving the threshold
    // at every
    // stepsTaken (scaling by stepsTaken instead of StartSteps) is NOT done here:
    // it would loosen the test as the sweep gets longer, defeating the very
    // mechanism that catches the late fall-back. The toy fixture of
    // gsALMExploration_test has it in closed form (u2 = +-sqrt(2(lambda-1))), and
    // example_ModifiedBratuExploration MEASURED the same law over three orders of
    // magnitude in SwitchLength (spread ~ C*sqrt(nPts*switchLen), C = 7.2..8.7).
    // The distance the predicate sees for a GENUINE sweep therefore scales with
    // sqrt(s), while the distance it sees for a sweep that FELL BACK onto the
    // parent is set by the parent's own sampling. A fixed tolerance is thus a
    // different criterion at every SwitchLength, and it was MEASURED to delete a
    // genuine branch once SwitchLength dropped to Length/40.
    //
    // Dividing by sqrt(Length) -- the arc length of the parent curves the sweep is
    // compared against -- makes the ratio dimensionless. MEASURED on the
    // modified-Bratu benchmark (OMP_NUM_THREADS=1, --sptestit 7, the predicate
    // value at the test step over a factor 10 in SwitchLength, each sweep
    // identified by its PROFILE, spread ~ 1e-13 = fell back onto the constant
    // branch, spread = O(1) = genuine branch B):
    //
    //   SwitchLength   genuine   /scale     fell back   /scale
    //   5e-3           1.606e-1  1.8545e-1  2.988e-2    3.450e-2
    //   2e-3           1.015e-1  1.8532e-1  1.400e-2    2.556e-2
    //   1e-3           7.168e-2  1.8507e-1  2.735e-3    7.062e-3
    //   5e-4           5.071e-2  1.8514e-1  2.893e-3    1.056e-2
    //
    // The genuine column collapses to 1.851e-1 within 0.2%: the sqrt law is the
    // right normalizer, and it is measured, not fitted. The window is therefore
    // (3.450e-2, 1.851e-1) at EVERY SwitchLength, and RetraceTol's default is its
    // geometric centre, sqrt(3.450e-2*1.851e-1) = 7.99e-2 -> 8e-2, i.e. 2.32x
    // above the retrace population and 2.32x below the genuine one -- the same
    // geometric-centre derivation the fixed default used, now applied in the
    // variable that does not move with the option values.
    //
    // SwitchLength, not the sweep's baseLen: isRetrace is only ever reached under
    // isBranchCurve (see the C_start block in traceSweep), and for a branch sweep
    // traceSweep's baseLen IS SwitchLength. Anyone who makes the test reachable
    // from a seed sweep must pass that sweep's own arc length in here instead --
    // otherwise the threshold would be scaled by a length the sweep never used.
    //
    // Complexity: O(1) per call. Called at every accepted step from
    // StartSteps onward (was: once per sweep) -- up to MaxPointsPerCurve-StartSteps+1
    // times per sweep (steps StartSteps..MaxPointsPerCurve inclusive).
    const T tol       = m_options.getReal("RetraceTol");
    const T switchLen = m_options.getReal("SwitchLength");
    const T length    = m_options.getReal("Length");
    const T steps     = static_cast<T>(m_options.getInt("StartSteps"));
    // Degenerate configurations must not turn the threshold into 0 or infinity:
    // either would classify EVERY sweep (or none) as a retrace, which is the
    // failure class this scaling exists to close. Fall back to the unscaled
    // tolerance instead.
    const T ratio = steps * math::abs(switchLen) / math::abs(length);
    if (!(ratio > (T)0) || !(ratio < (T)1e30))
        return tol;
    return tol * math::sqrt(ratio);
}

template <class T>
T gsALMExploration<T>::retraceDistance(const gsVector<T> & U, T L, index_t selfCurve,
                                       index_t & matchCurve, index_t & matchPoint) const
{
    // Diagnostic: exposes the isRetrace predicate's own distance instead
    // of only its boolean. min over all OTHER curves' stored points pt of
    // max( ||U-U_p||/max(1,||U_p||), |L-L_p|/max(1,|L_p|) ) -- EXACTLY the
    // quantity the old isRetrace body compared per-point via sameU && sameL
    // (max(ratioU,ratioL) < tol is the same test as (ratioU < tol && ratioL < tol),
    // and a minimum below tol exists iff some point's ratio is below tol -- see
    // isRetrace()'s doc for why this makes the refactor behaviour-preserving).
    //
    // Complexity: O(sum of stored points over all curves != selfCurve), same
    // single pass isRetrace() ran before it.
    matchCurve = -1;
    matchPoint = -1;
    T best = std::numeric_limits<T>::max(); // sentinel: no comparable point found
    for (index_t c = 0; c != static_cast<index_t>(m_landscape.nCurves()); ++c)
    {
        if (c == selfCurve)
            continue;
        const typename gsALMLandscape<T>::Curve & curve = m_landscape.curve(c);
        for (size_t p = 0; p != curve.points.size(); ++p)
        {
            const typename gsALMLandscape<T>::Point & pt = curve.points[p];
            if (pt.U.size() != U.size())
                continue;
            const T ratioU = (U - pt.U).norm() / math::max( (T)1, pt.U.norm() );
            const T ratioL = math::abs(L - pt.L) / math::max( (T)1, math::abs(pt.L) );
            const T ratio  = math::max(ratioU, ratioL);
            if (ratio < best)
            {
                best       = ratio;
                matchCurve = c;
                matchPoint = static_cast<index_t>(p);
            }
        }
    }
    return best;
}

template <class T>
bool gsALMExploration<T>::isRetrace(const gsVector<T> & U, T L, index_t selfCurve) const
{
    index_t matchCurve, matchPoint; // unused here; the caller that wants them uses retraceDistance directly
    const T dist = retraceDistance(U, L, selfCurve, matchCurve, matchPoint);
    return dist < retraceThreshold();
}

template <class T>
bool gsALMExploration<T>::traceCurve(const Job & job, std::queue<Job> & queue)
{
    const bool verbose = m_options.getSwitch("Verbose");

    const index_t cid = m_landscape.addCurve(job.parentCurve, job.parentPointIdx);

    // Branch jobs discovered while sweeping, with the index (within cid) of the
    // branch point each emanates from. Deferred to the end of this function
    // because assembling the two sweeps can renumber those indices; see traceSweep.
    std::vector<Job>     pending;
    std::vector<index_t> pendingPt;

    // --- Sweep 1: the job's own arc-length direction ---------------------------
    // false means the sweep was rewound as a retrace, which is a normal outcome
    // (the wrong-side half of a branch) as long as the OTHER sweep survives. What
    // the curve is worth is decided by the points it ends up holding, tested at
    // the end of this function -- but a curve that ends up EMPTY after a sweep was
    // rewound disappears, and a caller that queued a branch job at a correctly
    // classified branch point would then get neither a curve nor a diagnostic.
    // BOTH tallies are kept for exactly those two warnings: anyKept separates the
    // all-rewound case from the mixed one (a sweep rewound while a kept sweep
    // stored nothing), anyRewound separates BOTH of them from a pure solver
    // failure, in which the dedup made no decision and there is nothing to report.
    const bool kept1   = traceSweep(job, cid, job.backward, pending, pendingPt);
    bool       anyKept = kept1;
    bool    anyRewound = !kept1;
    const size_t nFirst = m_landscape.curve(cid).points.size();

    // --- Sweep 2: the OPPOSITE direction, into the SAME curve ------------------
    // Thesis Alg. 8.1 lines 17-23 / PyNCT SolutionDiagram: one (branch point,
    // tangent) yields ONE curve, swept from the same start state with the negated
    // initial tangent, the first sweep's points reverted and the second appended.
    // The historic code queued the two directions as two independent JOBS, so one
    // seed became two curve objects that were halves of a single locus, and one
    // branch mode became two curves of which the wrong-side one could only be
    // recognised after the fact by a point-coincidence test.
    if (job.bothDirections)
    {
        const bool kept2 = traceSweep(job, cid, !job.backward, pending, pendingPt);
        anyKept    = kept2 || anyKept;
        anyRewound = !kept2 || anyRewound;
    }

    // --- Assemble: revert the first sweep's block ------------------------------
    // Both sweeps run OUTWARD from the same start state, so appending them leaves
    // the point list running start->far, start->far (a jump at the junction).
    // Reverting the first block makes the assembled curve read far -> start ->
    // far, i.e. one ordered polyline, which is what every consumer (CSV row order,
    // plots, connectivity) assumes. Only done when the second sweep actually
    // contributed: otherwise the curve keeps exactly its historic point order.
    if (m_landscape.curve(cid).points.size() > nFirst && nFirst > 1)
    {
        typename gsALMLandscape<T>::Curve & c = m_landscape.curve(cid);
        std::reverse(c.points.begin(), c.points.begin() + nFirst);
        // Renumber the branch points of the reverted block: i -> nFirst-1-i.
        for (size_t i = 0; i != pendingPt.size(); ++i)
            if (static_cast<size_t>(pendingPt[i]) < nFirst)
                pendingPt[i] = static_cast<index_t>(nFirst) - 1 - pendingPt[i];
    }

    // --- Flush the deferred branch jobs, with their final parent indices -------
    for (size_t i = 0; i != pending.size(); ++i)
    {
        Job nj = pending[i];
        nj.parentPointIdx = pendingPt[i];
        queue.push(give(nj));
    }

    if (m_landscape.curve(cid).points.empty())
    {
        // A curve that vanishes because EVERY sweep was rewound as a retrace is
        // not a solver failure and not a plain empty curve: it is the dedup test
        // itself deciding that a branch the caller asked for is already in the
        // landscape. That decision must be visible OUTSIDE Verbose -- it is how a
        // too-large RetraceTol (or a nudge that never left the parent branch)
        // silently costs a consumer a branch it was entitled to. Not gsInfo: the
        // sweep-level notices are Verbose-gated on purpose, this one is not.
        //
        // THE MIXED CASE gets its own diagnostic (second branch below). If one
        // sweep is rewound as a retrace while the other is KEPT but stores zero
        // points -- every one of its steps failed until the arc length underflowed
        // -- then anyKept is true, yet the curve disappears just the same and the
        // dedup DID make a decision that contributed to the loss. Reporting it as
        // "every sweep retraced" would be false, and reporting nothing (the
        // behaviour before this branch existed) hides a dedup decision behind a
        // solver failure; the two causes are named separately instead, because
        // they have different remedies. A curve whose sweeps ALL merely failed
        // (nothing rewound) stays silent: no dedup decision was taken there, it is
        // an ordinary failed branch job, and saying otherwise would put a dedup
        // warning on every such run.
        if (!anyKept)
            gsWarn << "gsALMExploration: the branch curve emanating from point "
                   << job.parentPointIdx << " of curve " << job.parentCurve
                   << " (lambda = " << job.L << ") was DROPPED: every arc-length sweep "
                      "retraced an already-stored curve at RetraceTol = "
                   << m_options.getReal("RetraceTol") << " (effective threshold "
                   << retraceThreshold() << " at SwitchLength = "
                   << m_options.getReal("SwitchLength") << "). If a distinct branch is "
                      "expected here, lower RetraceTol, or check that the branch-switch "
                      "nudge (the solver's \"Perturbation\") actually leaves the parent "
                      "branch -- a sweep that falls back onto its parent is correctly "
                      "removed here, and no dedup setting can recover the branch it "
                      "never reached.\n";
        else if (anyRewound)
            gsWarn << "gsALMExploration: the branch curve emanating from point "
                   << job.parentPointIdx << " of curve " << job.parentCurve
                   << " (lambda = " << job.L << ") was DROPPED for TWO reasons: one "
                      "arc-length sweep retraced an already-stored curve at RetraceTol = "
                   << m_options.getReal("RetraceTol") << " (effective threshold "
                   << retraceThreshold() << " at SwitchLength = "
                   << m_options.getReal("SwitchLength") << ") and was rewound, while "
                      "the opposite sweep was kept but stored NO point (every step "
                      "failed until the arc length underflowed). Neither cause alone "
                      "would have emptied the curve. If a distinct branch is expected "
                      "here, check the solver side first (the kept sweep never took a "
                      "step), then RetraceTol and the branch-switch nudge (the solver's "
                      "\"Perturbation\").\n";
        m_landscape.removeLastCurve();
        return false;
    }

    if (verbose)
        gsInfo << "  [curve " << cid << "] " << m_landscape.curve(cid).points.size()
               << " points (" << (job.bothDirections ? 2 : 1) << " sweep(s)).\n";

    return true;
}

template <class T>
bool gsALMExploration<T>::traceSweep(const Job & job, index_t cid, bool backward,
                                     std::vector<Job> & pending,
                                     std::vector<index_t> & pendingPt)
{
    const bool    verbose          = m_options.getSwitch("Verbose");
    const bool    isBranchCurve    = (job.parentCurve != -1);
    const index_t maxPoints        = m_options.getInt("MaxPointsPerCurve");
    const index_t startSteps       = m_options.getInt("StartSteps");
    // Clause B of the C_start predicate below -- see RetraceHits' own doc in
    // defaultOptions() for the three regimes (<= 0 / 1 / 2) this single knob spans.
    const index_t retraceHits      = m_options.getInt("RetraceHits");
    // RetraceHits' own doc says "Never > 3" (the cost argument: every
    // extra step before the discard is one more step in which singular-point detection
    // can queue a branch job and suppress the discard via the spawnedJobs guard), but
    // that was prose only -- options() is public, so nothing stopped a caller from
    // setting a larger value and getting silent nonsense. GISMO_ENSURE, not
    // GISMO_ASSERT: this build configuration is Release (-O3 -DNDEBUG), where
    // GISMO_ASSERT is decorative.
    GISMO_ENSURE(retraceHits <= 3, "gsALMExploration: RetraceHits = " << retraceHits
                 << " is out of range; the C_start dedup predicate (clause B) is only "
                    "calibrated up to 3 consecutive below-threshold hits -- see the "
                    "RetraceHits option doc in defaultOptions().");
    const index_t branchPoints     = m_options.getInt("BranchPoints");
    // StepGrowth (pde2path sscontrol.m:37-43): opt-in x2 arc-length growth after an
    // accepted step, clamped at the base length. Reading the solver's MaxIter here is a
    // READ, not a write -- the explorer's "modifies NO solver option" contract
    // (gsALMExploration.h:157) is about writes.
    const bool    stepGrowth       = m_options.getSwitch("StepGrowth");
    const index_t growIterOpt      = m_options.getInt("StepGrowthIter");
    const index_t growIter         = (growIterOpt >= 0) ? growIterOpt
                                                         : m_solver->options().getInt("MaxIter")/2;

    // Base arc length: seed curves use Length, branch curves use SwitchLength.
    const T baseLen = isBranchCurve ? m_options.getReal("SwitchLength")
                                    : m_options.getReal("Length");
    const T len     = backward ? -baseLen : baseLen;

    // Diagnostic (alpha): the live solver state AS INHERITED FROM THE PREVIOUS
    // SWEEP, i.e. before ANY of the four seeding calls below runs. Verbose-gated, read-only
    // (gsALMBase.h's stabilityPrev()/solutionUPrev()/solutionLPrev()/getLengthPrev()/
    // stepTaken() accessors), no effect on any oracle's output. Compare against the beta
    // dump after the seeding block: the difference is exactly what the four setters reset.
    if (verbose)
        gsInfo << "  [curve " << cid << "] alpha (pre-seed): |Uprev|=" << m_solver->solutionUPrev().norm()
               << ", Lprev=" << m_solver->solutionLPrev() << ", lengthPrev=" << m_solver->getLengthPrev()
               << ", stepTaken=" << m_solver->stepTaken() << ", indicator=" << m_solver->indicator()
               << ", negatives=" << m_solver->negatives() << ", stabilityPrev=" << m_solver->stabilityPrev()
               << ".\n";

    // Seed the solver from the start state and reset the stability memory (as the
    // proven driver does before its loop). EVERY sweep re-seeds from the start
    // state: the second sweep of a curve must not inherit the first sweep's
    // far-end state, its inertia memory or its bisected arc length, or it would
    // trace on from the wrong end and compare its first point's inertia against a
    // point at the other end of the curve (a spurious singular-point detection at
    // the junction).
    //
    // TWO start-state kinds now exist, precedence job.tangentU (non-empty)
    // > job.Unudged (non-empty) > seed job (neither) -- see the Job doxygen.
    gsVector<T> Ustart;
    T           Lstart;
    // Diagnostic: the escape/nudge MAGNITUDE in force for this sweep,
    // reported by the C_start Verbose print below. Declared here (top of the
    // function) because the TANGENT branch's own `escape` local is scoped inside
    // that branch and would not otherwise survive to the C_start block further
    // down. For a FALLBACK-NUDGE job (job.Unudged non-empty), read the ACTUAL
    // magnitude back out of job.Unudged (== ||job.Unudged-job.U||, since the
    // branch-job creation block's `dir` is unit-normalised) rather than
    // recomputing 1/Perturbation here -- the two can now differ (the
    // ballFloorMargin floor at the branch-job creation site), so recomputing
    // would silently print a stale number whenever the floor is active. Seed
    // jobs (neither Unudged nor tangentU set) never reach the isBranchCurve-
    // gated print below, so their fallback value here is never displayed.
    T escapeUsed = (job.Unudged.size() != 0)
                 ? (job.Unudged - job.U).norm()
                 : (T)1 / m_solver->options().getReal("Perturbation");
    if (job.tangentU.size() != 0)
    {
        // TANGENT path: an EULER PREDICTOR STEP of one ESCAPE MAGNITUDE along s*tauHat
        // (s = +1 for the FORWARD sweep, -1 for backward), with the outward SECANT
        // installed so both families of predictor pick the outward direction.
        //
        // The magnitude is DECOUPLED from
        // baseLen/SwitchLength -- see the BranchEscape option doc in defaultOptions().
        // Tying the offset to the arc length (baseLen) is wrong by a power: the emanating branch leaves a singular
        // point like C*sqrt(s), not O(s)). tau1 is therefore only a DIRECTION here: it
        // is re-normalised to unit Euclidean length (job.tangentU/tangentL as stored may
        // already be close to unit length -- computeBranchTangent normalises by
        // distance(), a solver-specific metric that need not be Euclidean -- so this
        // re-normalisation makes the direction used here independent of that choice).
        const T tauNorm = math::sqrt(job.tangentU.squaredNorm()
                                     + job.tangentL*job.tangentL);
        gsVector<T> tauHatU = job.tangentU;
        T           tauHatL = job.tangentL;
        if (tauNorm > (T)0 && math::isfinite(tauNorm))
        {
            tauHatU /= tauNorm;
            tauHatL /= tauNorm;
        }
        // BranchEscape: -1 (default) reuses the UNFLOORED 1/Perturbation magnitude
        // (see switchBranch() / the Unudged construction in the branch-job creation
        // block); any positive value is used directly. Reading it here rather than
        // caching it keeps this in sync with a caller that changes the solver's
        // Perturbation option between solve() calls.
        //
        // This is the RAW 1/Perturbation, never RetraceBallFloor's floored
        // magnitude -- the fallback path's OWN nudge (see the branch-job creation
        // block) may be floored above 1/Perturbation while this TANGENT-path escape is
        // not, so the two paths can disagree in magnitude even at BranchEscape's -1
        // default. A caller wanting them to match must set BranchEscape explicitly.
        const T branchEscapeOpt = m_options.getReal("BranchEscape");
        const T escape = (branchEscapeOpt > (T)0)
                        ? branchEscapeOpt
                        : (T)1 / m_solver->options().getReal("Perturbation");
        escapeUsed = escape; // diagnostic: overwrite the fallback default above
        //  - setLength() BEFORE setPrevious() is the ORDER RULE documented at
        //    gsALMBase.h (setPrevious sets m_arcLength_prev = m_arcLength), and
        //    only that order makes the Riks / ConsistentCrisfield secant
        //    predictor's divisor equal the step, i.e.
        //    delta U = (U-Uprev)/Delta s_prev * Delta s = s*escape*tauHatU
        //    exactly (gsALMRiks.hpp's predictor). len KEEPS ITS OWN FORMULA (baseLen
        //    above): the arc length for the SUBSEQUENT sweep is unrelated to the
        //    predictor's escape offset.
        //  - gsALMCrisfield does not extrapolate the secant; it uses
        //    (m_DeltaUold,m_DeltaLold) -- which setPrevious() sets to the same
        //    outward secant -- to select the root of its quadratic, and
        //    initiateStep() does NOT clear them, so the same seeding serves
        //    both families.
        //  - With a nonzero secant, gsALMCrisfield::predictor() takes its else
        //    branch (not the fresh-start branch) -- intended.
        const T s = backward ? (T)(-1) : (T)(+1);
        Ustart = job.U + s*escape*tauHatU;
        Lstart = job.L + s*escape*tauHatL;
        m_solver->setSolution(Ustart, Lstart);
        m_solver->setLength(len);                 // BEFORE setPrevious -- see above
        m_solver->setPrevious(job.U, job.L);       // secant = s*escape*tauHat, m_arcLength_prev = len
        m_solver->setIndicator(0);
    }
    else
    {
        // FALLBACK-NUDGE path: EXACTLY today's four calls in today's order (bit-
        // identical seeding to the pre-tangent-path behaviour). Unudged non-empty for branch
        // jobs, job.U for seed jobs.
        Ustart = (job.Unudged.size() != 0) ? job.Unudged : job.U;
        Lstart = job.L;
        m_solver->setSolution(Ustart, Lstart);
        // Seed the previous state to the start state so every curve begins in
        // the fresh-start predictor branch (prev==current => zero secant). Without
        // this a non-rest seed (backward/branch job) inherits a stale m_Uprev and the
        // secant predictor degenerates, freezing the trace at the seed.
        m_solver->setPrevious(Ustart, Lstart);
        m_solver->setIndicator(0);
        m_solver->setLength(len);
    }

    // Diagnostic (beta): the live solver state immediately AFTER the seeding
    // block. The alpha/beta difference is exactly what the four setters reset.
    if (verbose)
        gsInfo << "  [curve " << cid << "] beta (post-seed): |Uprev|=" << m_solver->solutionUPrev().norm()
               << ", Lprev=" << m_solver->solutionLPrev() << ", lengthPrev=" << m_solver->getLengthPrev()
               << ", stepTaken=" << m_solver->stepTaken() << ", indicator=" << m_solver->indicator()
               << ", negatives=" << m_solver->negatives() << ", stabilityPrev=" << m_solver->stabilityPrev()
               << ".\n";

    // Diagnostic: the STEP-0 (predictor-state) half of the
    // discriminating measurement -- ||Ustart-job.U|| and the retrace distance of the
    // predictor state itself, BEFORE the first corrector step, for a branch sweep.
    // Verbose-gated (the explorer's own local, already ON in the ModifiedBratu
    // drivers); O(sum of stored points over all curves != selfCurve), same cost class
    // as the per-step C_start probe below, paid once per sweep.
    if (verbose && isBranchCurve)
    {
        index_t step0MatchCurve = -1, step0MatchPoint = -1;
        const T step0Dist = retraceDistance(Ustart, Lstart, cid, step0MatchCurve, step0MatchPoint);
        gsInfo << "  [curve " << cid << "] step-0 predictor state: ||Ustart-job.U||="
               << (Ustart - job.U).norm() << ", retrace distance=" << step0Dist
               << " vs threshold " << retraceThreshold() << ", matched curve "
               << step0MatchCurve << " point " << step0MatchPoint << ".\n";
    }

    // Where THIS sweep's points and jobs start: a retrace rewind returns here and
    // never touches what an earlier sweep of the same curve contributed.
    const size_t sweepStartPt  = m_landscape.curve(cid).points.size();
    const size_t sweepStartJob = pending.size();

    // Pre-crossing point memory: seeds the extended singular-point solve.
    // job.L must NOT survive as the start load: on the TANGENT path a tangent
    // moves Lambda too (Lstart != job.L in general), and leaving Lold = job.L
    // would seed the next singular-point solve with an inconsistent pre-crossing
    // pair.
    gsVector<T> Uold = Ustart;
    T           Lold = Lstart;

    T       dLb         = len;
    const T dLb0        = len;
    bool    bisected    = false;
    bool    spawnedJobs = false; // true once this sweep has queued a branch job
    // Latches the "keeping sweep to preserve parent references" Verbose
    // notice to ONE print per sweep. With the C_start test now running at every
    // step from StartSteps onward (not once), a persistently-retracing sweep with
    // spawnedJobs already true would otherwise reprint the same notice every
    // remaining step.
    bool    spawnedJobsNoticed = false;
    // Length of the current CONSECUTIVE run of below-threshold C_start
    // evaluations, i.e. clause B's counter -- see the C_start block's own comment for
    // the exact semantics (reset on any at-or-above-threshold evaluation, never reset
    // by the spawnedJobs guard, never carried across sweeps: a fresh traceSweep call is
    // a fresh sweep and a fresh counter).
    index_t belowRun = 0;
    // Counts ACCEPTED points of THIS sweep. MaxPointsPerCurve is therefore a
    // per-sweep budget: a two-sweep curve may hold up to 2*MaxPointsPerCurve
    // points. Deliberate -- a shared budget would halve each direction's reach and
    // truncate rows that the one-curve-per-direction code traced in full.
    index_t stepsTaken  = 0;

    // Falling out of the while condition below is the budget exit; the two
    // arc-length-underflow breaks overwrite this before breaking.
    SweepTermination reason = SweepTermination::PointBudget;

    // STORED landscape stability of consecutive ACCEPTED points (not
    // solver-side corrector-iterate transients). 0 = no prior accepted point yet.
    // This value drives the stored `stab` semantics only, NOT the detection trigger.
    index_t stabPrevAccepted = 0;

    // Inertia-based singular-point detection: the trigger fires on a
    // change of the NEGATIVE-EIGENVALUE COUNT (tangent inertia) between consecutive
    // ACCEPTED points. This strictly generalizes the old min-eigenvalue sign flip (a
    // 0<->1 count change), and additionally catches SECONDARY bifurcations on
    // already-unstable branches where the 2nd/3rd/... eigenvalue crosses zero while
    // the minimum stays negative. -1 = no prior accepted point yet.
    index_t negPrevAccepted = -1;

    while (stepsTaken < maxPoints)
    {
        gsStatus status;
        try { status = m_solver->step(); }
        catch (...) { status = gsStatus::AssemblyError; } // NaN metric on a bad predictor

        // --- Step-fail retry: halve the arc length and re-seed from (Uold,Lold) ---
        //
        // The predicate is a CATCH-ALL, not an enumeration of failure statuses
        // Three things motivate that:
        //
        //  * WHAT A FAILED STEP WOULD OTHERWISE RECORD. A step that does not
        //    return Success commits nothing: m_U += m_DeltaU happens in
        //    iterationFinish() only (gsALMLoadControl<T>::iterationFinish and its
        //    siblings), which gsALMBase<T>::_step() reaches only after the
        //    convergence test. So solutionU()/solutionL() below would return the
        //    PREVIOUS state bit-exactly, and the "converged step" branch would
        //    store the curve's own seed as if it were a traced point -- carrying a
        //    `stab` read off an m_jacMat belonging to some other state.
        //  * WHY THE PROGRESS GUARD BELOW IS NOT A SUBSTITUTE. It rejects exactly
        //    that zero-length move, but it is skipped while stepsTaken == 0, so it
        //    never sees a curve's FIRST step. Measured on
        //    example_ModifiedBratuExploration --sptestit 7 -m 2, where both legs of
        //    the c = 1.5 seed returned SolverError on their first step and each
        //    stored one phantom row equal to the seed.
        //  * WHY NOT `|| status == gsStatus::SolverError`. Over gsStatus
        //    (gsStructuralAnalysisTypes.h) the catch-all additionally admits only
        //    OtherError -- an unknown exception escaping step(), which must
        //    certainly not be recorded -- and NotStarted, unreachable after step().
        //    It is therefore behaviourally identical to the enumerated form today,
        //    while being closed against the defect's own class: an enumerated predicate
        //    silently drifts out of sync whenever a status is added to the enum without
        //    updating it. It also agrees with the catch(...) above, which already
        //    funnels ANY escaping exception into this same retry.
        //
        // A curve whose every step fails now stores nothing and is dropped by the
        // points.empty() branch at the end of this function, so no empty curve
        // reaches the landscape, the CSV or the HDF5 file.
        if (status != gsStatus::Success)
        {
            dLb = dLb / (T)2;
            // NEGATED form of `< 1e-6` -- same ratio test, same constant, identical for
            // every finite value, but it also fires on NaN. dLb0 = 0 (a curve traced at
            // zero length) makes dLb/dLb0 = 0/0 = NaN, and `NaN < 1e-6` is FALSE, which
            // would spin this retry loop forever. Matches the form used at the
            // two corresponding sites in gsAPALM<T>::_initiation / _correction, which were
            // copied FROM here; this is the original catching up with its own copies.
            if (!(math::abs(dLb / dLb0) >= (T)1e-6))
            {
                // Graceful termination: keep the partial curve.
                if (verbose)
                    gsInfo << "  [curve " << cid << "] arc length underflow; keeping partial curve ("
                           << stepsTaken << " points).\n";
                reason = SweepTermination::UnderflowStepFailed;
                break;
            }
            m_solver->setLength(dLb);
            m_solver->setSolution(Uold, Lold);
            bisected = true;
            if (verbose)
                // Carry the status: "did not converge" misattributes a SolverError
                // (a failed factorization/solve, i.e. a step never attempted to
                // completion) to the Newton corrector. See gsStatus in
                // gsStructuralAnalysisTypes.h for the value.
                gsInfo << "  [curve " << cid << "] step failed (status " << (int)status
                       << "); halving arc length to " << dLb << ".\n";
            continue; // retry without consuming the point budget
        }

        // --- Converged step: record the point ---------------------------------
        // State fidelity: jacobian=TRUE. With jacobian=false this
        // call reused the corrector's stale m_jacMat -- a tangent belonging to the LAST
        // CORRECTOR ITERATE, not to the committed point (m_U,m_L) whose stability is
        // about to be stored. Every stored +-1 `stab` and every `negatives()` count read
        // below was therefore classified PRE-COMMIT. Re-assembling here costs one extra
        // Jacobian per ACCEPTED point (not per corrector iteration).
        // MEASURED NEUTRAL: applying exactly this edit, all six CSV
        // oracles came back byte-identical -- so those runs are ALSO useless as a
        // regression gate for it. Note gsALMBase<T>::_step()'s per-iteration
        // computeStability(false) calls are deliberately left alone; they are a separate,
        // UNMEASURED question.
        m_solver->computeStability(true);

        const gsVector<T> Ucur = m_solver->solutionU();
        const T           Lcur = m_solver->solutionL();

        // --- Progress guard: reject a converged step that did not move ---------
        // A step that returns a point bit-identical (within moveTol) to the
        // previous accepted state (Uold,Lold) is treated EXACTLY like a
        // non-converged step: halve the arc length and re-seed. This prevents a
        // stalled solver from re-accepting the same point up to MaxPointsPerCurve;
        // the arc-length underflow check then terminates the curve gracefully.
        // Only applies once the curve has >=1 accepted point: a branch-job's first
        // corrected point can legitimately sit close to its nudged seed.
        if (stepsTaken > 0)
        {
            // Part B: own option, no longer DedupTol*1e-2. Default 1e-6
            // reproduces the old effective value at the default DedupTol = 1e-4.
            const T moveTol = m_options.getReal("MoveTol");
            const bool noProgress =
                ((Ucur - Uold).norm() <= moveTol * math::max((T)1, Uold.norm())) &&
                (math::abs(Lcur - Lold) <= moveTol * math::max((T)1, math::abs(Lold)));
            if (noProgress)
            {
                dLb = dLb / (T)2;
                // Negated form, see the identical break in the step-fail retry above: it
                // fires on NaN too (dLb0 = 0 gives dLb/dLb0 = NaN).
                if (!(math::abs(dLb / dLb0) >= (T)1e-6))
                {
                    // Graceful termination: keep the partial curve.
                    if (verbose)
                        gsInfo << "  [curve " << cid << "] arc length underflow; keeping partial curve ("
                               << stepsTaken << " points).\n";
                    reason = SweepTermination::UnderflowNoProgress;
                    break;
                }
                m_solver->setLength(dLb);
                m_solver->setSolution(Uold, Lold);
                bisected = true;
                if (verbose)
                    gsInfo << "  [curve " << cid << "] no progress; halving arc length to " << dLb << ".\n";
                continue; // retry without consuming the point budget
            }
        }

        const T           curIndicator = m_solver->indicator();
        const index_t     stab = (curIndicator > 0) ? +1 : -1;
        // Tangent inertia (negative-eigenvalue/pivot count) from computeStability
        // above; drives inertia-based singular-point detection below.
        const index_t     negCur = m_solver->negatives();
        // Captured HERE, not in the arc-length tail below: computeSingularPoint's
        // bisection runs its own step()s and overwrites m_numIterations.
        const index_t     itersCur = m_solver->numIterations();

        // AUTO-07p's fold test function at the accepted point (gsALMBase::foldTestFunction).
        // INSTRUMENT ONLY: it is printed, never tested -- detection stays on the inertia
        // count above. NaN means "not available at this state" (see the accessor's doc).
        if (verbose && !math::isnan(m_solver->foldTestFunction()))
            gsInfo << "  [curve " << cid << "] fold test function = "
                   << m_solver->foldTestFunction() << ".\n";

        gsMultiPatch<T> deformed;
        bool haveGeometry = false;
        if (m_solutionConstructor)
            haveGeometry = m_solutionConstructor(Ucur, deformed);

        // /*equilibrium=*/true is a spelling change, not a semantic one: it
        // reproduces the declaration's existing default, which this call was
        // already getting. Appending a trailing negatives argument simply
        // forces the intervening arguments to be spelled out; the VALUE here
        // is untouched.
        m_landscape.addPoint(cid, Ucur, Lcur, stab, haveGeometry ? &deformed : nullptr,
                             /*isBifurcation=*/false, /*equilibrium=*/true,
                             /*negatives=*/negCur);
        ++stepsTaken;

        // --- C_start dedup (simplified eq. 8.1): rewind a retracing SWEEP ------
        // Under one-curve-per-direction this removed the whole curve, because the
        // curve WAS the sweep. Now it rewinds only this sweep: the wrong-side half
        // of a branch is dropped while the genuine half traced by the other sweep
        // of the same curve is kept. That is the mechanism which stops the
        // fall-back half from surviving as an appendage of the child curve.
        // This splits the isRetrace evaluation from the predicate's result so the
        // safety-net telemetry (m_retraceTests/m_retraceFired) counts every EVALUATION,
        // not just every FIRING -- counting inside the old short-circuited `&&` condition
        // would give a 100%-or-undefined rate (isRetrace is only ever reached, hence only
        // ever counted, when it fires).
        //
        // An earlier measurement identified the discriminator; the two-clause predicate
        // below implements it. A SINGLE evaluation at stepsTaken == startSteps missed a phantom sweep that
        // reads ABOVE threshold at StartSteps but falls back onto the parent branch a
        // few steps later (MEASURED on the modified-Bratu A-C crossing at bare default
        // flags, --sptestit 7: curve 3's surviving sweep reads factor 2.35 at step 3,
        // rises to 2.44 (step 4) then 3.52 (step 7) mimicking a genuine branch, and only
        // at step 8 collapses to factor 0.45 -- below threshold -- as it lands back on
        // curve 1's own dense sampling; every step from 8 to 20 stays below 1). The
        // GENUINE sweep (curve 2) was confirmed, over its own full sweep, to be safe
        // against this widened window: its ratio/threshold factor is 2.32 at step 3 and
        // rises, NEAR the branch point and OVER THE CONFIGURATIONS MEASURED, to 5.60 at
        // step 20 -- see the caveat in gsALMExploration.h's class doxygen, RetraceTol
        // comparison-scope paragraph, for the shape (closely approaching distinct
        // branches) that would need RetraceTol revisited instead.
        //
        // What was originally shipped from that measurement was WEAKER than the signal it
        // found: "one hit anywhere >= StartSteps discards the sweep" -- any below-
        // threshold reading, however isolated, deletes the whole sweep. The measurement
        // itself identified PERSISTENCE as the discriminator: the phantom's kept sweep
        // stays below threshold at EVERY step from its first hit (8) through the end of
        // the sweep (20), whereas a genuine sweep that merely crosses another stored
        // curve transversally is expected to produce a SINGLE below-threshold reading
        // before resuming its climb (untested in this tree -- no oracle has that shape).
        // A one-hit predicate cannot tell those two apart; a run-length
        // predicate can.
        //
        // The shipped predicate is therefore TWO clauses (RetraceHits, see its own doc
        // in defaultOptions() for the <=0 / 1 / 2 regimes this one knob spans):
        //   clause A: stepsTaken == startSteps, below-threshold discards unconditionally
        //             -- exactly the original decision at that one step. This
        //             is NOT a hedge, it is a bit-identical-to-the-original
        //             GUARANTEE: two sweeps in oracle b fire here (both children's
        //             wrong-side sweep), and nothing in this tree has ever measured step
        //             startSteps+1 of those sweeps -- the original predicate always truncated
        //             the trace the moment step StartSteps fired. Dropping clause A would
        //             make that decision depend on an unmeasured next step.
        //   clause B: RetraceHits CONSECUTIVE below-threshold evaluations past
        //             StartSteps (belowRun >= RetraceHits, RetraceHits > 0).
        //             "Consecutive" means consecutive EVALUATIONS of this very block,
        //             which -- because the block runs immediately after
        //             addPoint()/++stepsTaken above -- is the same thing as consecutive
        //             ACCEPTED steps: a failed/retried corrector step or a
        //             MoveTol-rejected step never reaches here, so it cannot break a run.
        // belowRun (declared with the sweep locals above) is reset to 0 by ANY evaluation
        // that reads dist >= thr, and is never carried across sweeps (a new traceSweep
        // call is a new sweep, hence a fresh counter starting at 0).
        //
        // Interaction with the spawnedJobs guard below, stated so nobody reads it as an
        // oversight: a sweep the guard keeps (spawnedJobs true) neither resets belowRun
        // nor returns, so belowRun keeps accumulating and `retrace` stays true at every
        // subsequent step for as long as the sweep keeps reading below threshold. That is
        // intended -- the guard decides whether to ACT on `retrace`, the counter does not
        // decide whether `retrace` fires -- and the telemetry latch (spawnedJobsNoticed)
        // keeps m_retraceFired at one per sweep regardless. Do NOT add a belowRun reset in
        // the guard branch: it would change the step at which a later discard fires.
        //
        // `stepsTaken >= startSteps` (the evaluation guard, unchanged) is
        // still the floor on when scrutiny starts; only the DECISION inside it changed,
        // from the original one-hit rule to the two-clause rule above.
        // Complexity: up to MaxPointsPerCurve-StartSteps+1 evaluations per sweep (steps
        // StartSteps..MaxPointsPerCurve inclusive; was exactly one originally), each
        // O(sum of stored points over all curves != selfCurve) -- see retraceDistance().
        bool retrace = false;
        if (isBranchCurve && stepsTaken >= startSteps)
        {
            ++m_retraceTests;
            // Diagnostic: call retraceDistance() directly (instead of the
            // isRetrace() boolean wrapper) so the min ratio, the matched
            // curve/point and the threshold are all available to print below --
            // isRetrace() itself is unchanged and still used everywhere else.
            index_t matchCurve = -1, matchPoint = -1;
            const T    dist  = retraceDistance(Ucur, Lcur, cid, matchCurve, matchPoint);
            const T    thr   = retraceThreshold();
            const bool below = dist < thr;
            belowRun = below ? belowRun + 1 : 0;
            retrace = below && (  stepsTaken == startSteps                        // clause A
                                || (retraceHits > 0 && belowRun >= retraceHits) ); // clause B
            // Count a PREDICATE FIRING at most once per sweep, guarded by
            // !spawnedJobsNoticed: without that guard ++m_retraceFired would run at EVERY
            // remaining step of a sweep the spawnedJobs guard below keeps, so one
            // persistently-retracing sweep could dominate the printed "fired N of M" rate.
            // m_retraceTests still counts every EVALUATION, unchanged.
            if (retrace && !spawnedJobsNoticed) ++m_retraceFired;
            if (verbose)
                gsInfo << "  [curve " << cid << "] retrace probe: step " << stepsTaken
                       << ", ratio " << dist << " vs threshold " << thr << " (factor "
                       << (thr > (T)0 ? dist / thr : dist) << "), matched curve "
                       << matchCurve << " point " << matchPoint << ", escape = "
                       << escapeUsed << ", SwitchLength = "
                       << m_options.getReal("SwitchLength")
                       << ", below-threshold run = " << belowRun << " of " << retraceHits << ".\n";
        }
        if (retrace)
        {
            // If a branch job was already queued from this sweep (possible once the
            // test above runs at every step past StartSteps, not only once), rewinding
            // now would orphan it. Keep the sweep in that case (harmless duplicate) to
            // preserve parents; the notice is latched to fire ONCE per sweep (not once
            // per remaining step) via spawnedJobsNoticed below.
            if (spawnedJobs)
            {
                if (verbose && !spawnedJobsNoticed)
                    gsInfo << "  [curve " << cid << "] retrace detected but branch jobs already "
                           << "queued from it; keeping sweep to preserve parent references.\n";
                spawnedJobsNoticed = true;
            }
            else
            {
                if (verbose)
                    gsInfo << "  [curve " << cid << "] duplicate retrace detected after "
                           << stepsTaken << " steps; discarding this sweep ("
                           << stepsTaken << " points).\n";
                // Make the data loss visible: a discard can now drop up to
                // MaxPointsPerCurve accepted points (a PER-SWEEP budget, see its own doc in
                // defaultOptions()) rather than just <= StartSteps; in the MIXED case
                // -- this sweep's sibling survives WITH points of its own -- nothing else
                // warns (traceCurve's two gsWarn's are gated on !anyKept and on "surviving
                // sweep stored no point at all" respectively, neither of which is this
                // case). A singular point MARKED during this window --
                // singular-point detection runs AFTER this block, below -- sets
                // points.back().isBifurcation = true (a confirmed limit/fold point via
                // markBifurcation(), OR a classification-FAILED point via
                // markUnresolvedSingular(), see the singular-point-detection block below and
                // gsALMLandscape.hpp) and stores NO new point, so it is erased here together
                // with the rest of the discarded sweep, leaving a Verbose/gsWarn line that
                // named a point that no longer exists. A BRANCH point cannot do this: it sets
                // spawnedJobs, which routes detection here through the spawnedJobs-guard
                // branch above instead, never this one. droppedSingular counts how many such
                // marked points this discard erased. MEASURED NONZERO on oracle
                // b's own shipping gate (bare default flags): curve 3's kept sweep detects an
                // inertia flip at step 8 (its first below-threshold hit), classification
                // FAILS ("critical mode unresolved"), markUnresolvedSingular() marks that same
                // step's point, and the clause-B discard at step 9 erases it --
                // droppedSingular == 1, dropped == 9. This path cannot fire under a predicate
                // that discards at the FIRST below-threshold step (here, step 8): such a
                // predicate would return before this detection block ever ran; RetraceHits'
                // extra step is what exposes it. A
                // zero value on a given run means no marked point fell inside that run's
                // discard window -- not a guarantee for every configuration; either way it
                // must be reported, not left implicit.
                typename gsALMLandscape<T>::Curve & c = m_landscape.curve(cid);
                const index_t dropped = static_cast<index_t>(c.points.size() - sweepStartPt);
                index_t droppedSingular = 0;
                for (size_t i = sweepStartPt; i != c.points.size(); ++i)
                    if (c.points[i].isBifurcation) ++droppedSingular;
                if (dropped > startSteps)
                    gsWarn << "gsALMExploration: curve " << cid << "'s sweep discarded at step "
                           << stepsTaken << " by the C_start dedup dropped " << dropped
                           << " accepted point(s) (" << droppedSingular << " of them carrying a "
                           << "singular-point mark, isBifurcation) at RetraceTol = "
                           << m_options.getReal("RetraceTol")
                           << " (effective threshold " << retraceThreshold() << "); this "
                           << "curve's sibling sweep may still survive, with or without "
                           << "points of its own.\n";
                c.points.erase(c.points.begin() + sweepStartPt, c.points.end());
                pending.erase  (pending.begin()   + sweepStartJob, pending.end());
                pendingPt.erase(pendingPt.begin() + sweepStartJob, pendingPt.end());
                m_sweepTermination.push_back(std::make_pair(cid, SweepTermination::RetraceRewound));
                return false;
            }
        }

        // --- Singular point detection -----------------------------------------
        // Inertia-based trigger: fire whenever the negative-eigenvalue count of the
        // tangent changes between the previous and current ACCEPTED point. This
        // subsumes the old min-eigenvalue sign flip (a 0<->1 count change) and also
        // catches secondary bifurcations on unstable branches (e.g. 1->2 negatives).
        // Skip detection for the first two steps of a branch curve: the switched
        // state starts near-singular and would re-detect the parent bifurcation.
        const bool flipped = (negPrevAccepted != -1) && (negCur != negPrevAccepted);
        const bool detect  = flipped && (!isBranchCurve || stepsTaken > 2);

        // A flip that the branch-curve guard suppresses must be DEFERRED, not absorbed:
        // while suppression is active the reference inertia is HELD at its pre-switch
        // value (see the memory update at the end of the loop), so an inertia change
        // that PERSISTS is detected at the first unsuppressed step, while one that has
        // already reverted is correctly not detected.
        const bool suppressed = flipped && !detect;
        if (suppressed && verbose)
            gsInfo << "  [curve " << cid << "] inertia flip (" << negPrevAccepted
                   << " -> " << negCur << ") suppressed by the branch-curve guard at step "
                   << stepsTaken << "; holding the reference inertia (deferred).\n";

        if (detect)
        {
            // A change of MORE than one in the accepted-point inertia is a multi-step
            // crossing (e.g. 1->3): _localizeCrossing narrows its bracket to keep
            // negLo != negHi and returns exactly ONE crossing, so at most one branch
            // is followed and the reduction is otherwise silent. Explicit integer
            // comparison, not math::abs -- every other use of math::abs in this file
            // is on a value of type T, not on an index_t inertia difference.
            const index_t dneg = negCur - negPrevAccepted;
            const bool multReduced = (dneg > 1) || (dneg < -1);
            if (multReduced)
                gsWarn << "gsALMExploration: multi-step inertia crossing on curve " << cid
                       << " over lambda in [" << Lold << "," << Lcur << "] (inertia "
                       << negPrevAccepted << " -> " << negCur << ") is reduced to ONE "
                          "simple crossing: at most one branch is followed. "
                          "gsALMBase::computeBranchTangent resolves the emanating "
                          "tangent only at a simple (multiplicity-1) branch point. The "
                          "stored singular point is flagged multiplicityReduced = true.\n";

            // --- Localize: bisect the accepted-point interval on the tangent -------
            // inertia (pde2path bifdetec.m:42-99 shape) BEFORE classifying, instead of
            // classifying at the (possibly coarse) accepted-point grid on a tangent
            // assembled at the wrong end of it (see the class doxygen).
            gsVector<T> Uloc; T Lloc; T bracket;
            index_t probes = 0, probesTotal = 0;
            bool localized = _localizeCrossing(Uold, Lold, negPrevAccepted,
                                               Ucur, Lcur, negCur, dLb, dLb0, cid,
                                               Uloc, Lloc, bracket, &probes);
            probesTotal = probes;

            // A failed localization loses the crossing permanently unless retried. Retry with
            // both bisection bounds relaxed: BisecMax doubled reopens the "BisecMax"
            // exit, BisecLengthFloor halved reopens the retreat-below-floor exit, and
            // the routine is deterministic, so relaxing BOTH is what makes a second
            // attempt able to reach a different outcome. AUTO discards-and-re-probes;
            // pde2path saves-or-discards on |Re mu|.
            const index_t nRetries = m_options.getInt("LocalizeRetries");
            index_t attemptsUsed = 0;
            for (index_t attempt = 1; !localized && attempt <= nRetries; ++attempt)
            {
                attemptsUsed = attempt;
                const index_t bisecMax0 = m_options.getInt ("BisecMax");
                const T       floor0    = m_options.getReal("BisecLengthFloor");
                m_options.setInt ("BisecMax",         bisecMax0 * 2);
                m_options.setReal("BisecLengthFloor", floor0 / (T)2);
                try
                {
                    localized = _localizeCrossing(Uold, Lold, negPrevAccepted,
                                                  Ucur, Lcur, negCur, dLb, dLb0, cid,
                                                  Uloc, Lloc, bracket, &probes);
                }
                catch (...)
                {
                    m_options.setInt ("BisecMax",         bisecMax0);
                    m_options.setReal("BisecLengthFloor", floor0);
                    throw;
                }
                m_options.setInt ("BisecMax",         bisecMax0);
                m_options.setReal("BisecLengthFloor", floor0);
                probesTotal += probes;
                if (verbose)
                    gsInfo << "  [curve " << cid << "] localization retry " << attempt
                           << "/" << nRetries << " (BisecMax " << bisecMax0*2
                           << ", BisecLengthFloor " << floor0/(T)2 << "): "
                           << (localized ? "localized." : "still not localized.") << "\n";
            }

            if (!localized)
            {
                // LOCALIZATION-FAILURE. Not a single bisection probe
                // converged, so no point near the crossing is available to classify
                // or store; the honest "unresolved singular point" marking attaches
                // here. Fall back to the coarse post-crossing point, as today -- it is
                // flagged isBifurcation=true, equilibrium=false, unresolved=true, so it
                // is never indistinguishable from a certified bifurcation. The failed
                // bracket (Lold,Lcur) and the total probe count over the first attempt
                // and every retry are recorded on the landscape point as a hook for a
                // later resume; nothing consumes them yet.
                m_landscape.markUnresolvedSingular(cid, Lold, Lcur, probesTotal);
                gsWarn << "gsALMExploration: possible singular point near lambda = " << Lcur
                       << " on curve " << cid << " -- NOT refined (bisection localization "
                          "failed after " << (1 + attemptsUsed) << " attempt(s), "
                       << probesTotal << " probe(s) total: no probe converged in the "
                          "interval lambda in [" << Lold
                       << "," << Lcur << "]); no branch followed. The traced point is "
                          "stored with unresolved = true and equilibrium = false; it is "
                          "NOT a certified bifurcation.\n";
                if (verbose)
                    gsInfo << "  [curve " << cid << "] singular point LOCALIZATION FAILED "
                           << "after " << (1 + attemptsUsed) << " attempt(s), " << probesTotal
                           << " probe(s) total (no bisection probe converged, interval lambda in ["
                           << Lold << "," << Lcur << "]); marking post-crossing point as "
                           << "unresolved; no branch jobs.\n";
            }
            else
            {
                // Classify at the LOCALIZED point on a tangent freshly assembled
                // there (jacobian=true). This is also the single classification call
                // left in the file, and it seeds m_V for the extended solve's initial
                // null-vector guess below -- nothing may touch the solver between
                // this call and computeSingularPoint on the branch-point path.
                //
                // Fix 4: restore the SWEEP's
                // arc length here, BEFORE seeding the singular-point solve.
                // _localizeCrossing's bisection probes leave m_solver's arc length at
                // whatever the LAST probe used -- a probe-sized length (bracket width
                // ~4.9e-5 against a base dLb0 of 0.05, ~1000x smaller, MEASURED both on
                // the pitchfork fixture and on example_BratuExploration) -- and
                // computeSingularPoint's own bisection stage (SingularPointComputeTolB)
                // steps at, and its internal step()s therefore move by, whatever arc
                // length is installed when it starts. Leaving a probe-sized length in
                // place made that stage's own steps tiny and, MEASURED on the pitchfork
                // unit-test fixtures (8 branch-point solves, MaxIter=50), it consumed 22
                // iterations before this restore and 11 after (identical across all 8;
                // the unrelated riks_fold_extended_solve fixture stayed at 49 either
                // way, see the task report) -- a real, reproducible ~2x reduction, NOT
                // full-budget exhaustion on this toy problem. No correctness impact was
                // observed -- the worst case leaves (U,L) at the already-localized point
                // -- this is a wasted-work / robustness fix only.
                m_solver->setLength(dLb);
                m_solver->setSolution(Uloc, Lloc);
                bool isBranch = false;
                bool classificationThrew = false;
                try { isBranch = m_solver->isBifurcation(/*jacobian=*/true); }
                catch (...) { classificationThrew = true; } // -> honest-failure handling below

                if (!isBranch)
                {
                    // The "mode unresolved" contract: isBifurcation(true) already
                    // returns false (indistinguishable from a LIMIT verdict) when the
                    // critical mode failed to resolve; singularPointVerdict() is the
                    // machine-readable discriminator, queried right after the call.
                    const bool unresolved = !classificationThrew &&
                        (m_solver->singularPointVerdict() == gsALMBase<T>::SPverdict::Unresolved);

                    if (classificationThrew || unresolved)
                    {
                        // The SAME honest-failure marking as the localization
                        // failure above, for the two other ways the classification at
                        // the localized point can fail to produce a verdict. This also
                        // routes the "mode unresolved" classification outcome
                        // (SPverdict::Unresolved) to the same honest path: it must NOT
                        // fall through to the limit-point branch below, which would
                        // silently assert a fold verdict the classifier explicitly
                        // refused to give.
                        m_landscape.markUnresolvedSingular(cid);
                        gsWarn << "gsALMExploration: possible singular point near lambda = " << Lcur
                               << " on curve " << cid << " -- NOT refined (classification at "
                                  "the localized point lambda=" << Lloc << " "
                               << (classificationThrew ? "threw" : "left the critical mode "
                                  "unresolved")
                               << "); no branch followed. The traced point is stored with "
                                  "unresolved = true and equilibrium = false; it is NOT a "
                                  "certified bifurcation.\n";
                        if (verbose)
                            gsInfo << "  [curve " << cid << "] singular point classification "
                                   << "FAILED at the localized point (lambda=" << Lloc << "; "
                                   << (classificationThrew ? "threw" : "critical mode unresolved")
                                   << "); marking post-crossing point as unresolved; no "
                                      "branch jobs.\n";
                    }
                    else
                    {
                        // Limit / turning point, classified at the LOCALIZED point with
                        // a fresh tangent.
                        if (verbose)
                            gsInfo << "  [curve " << cid << "] limit point at lambda=" << Lloc << ".\n";
                        m_landscape.markBifurcation(cid);
                    }
                }
                else
                {
                    // Branch point: extended solve seeded from the LOCALIZED point.
                    //
                    // EQUILIBRIUM PROVENANCE of the point stored below. The extended system's
                    // BASIC termination test is ||K_T.V|| < SingularPointComputeTolE, which
                    // certifies SINGULARITY only: the set on which it holds is in general a
                    // whole CURVE through the branch point -- on the modified-Bratu benchmark
                    // the tangent is singular along mu*exp(c) = c_1 for EVERY c -- so the test
                    // bounds nothing at all about ||R(U*,L*)|| unless SingularPointComposite is
                    // on. Measured there (SingularPointComposite off), the returned point sits
                    // at ||R||/(1e-6 max(1,|L| |f|)) = 168, while every point produced by the
                    // ordinary corrector sits at ~1e-7 of the same bound.
                    //
                    // Semantics (locked): Point::equilibrium is true iff the row
                    // satisfies the residual test that was USED TO ACCEPT it. A converged
                    // extended solve (spStatus == Success) met exactly that -- its own
                    // termination test, whichever one is active -- so it is stored
                    // equilibrium = true REGARDLESS of spCertified below; a solve that did NOT
                    // converge is stored equilibrium = false (see the failure branch, which
                    // also sets unresolved = true). spCertified therefore no longer decides the
                    // STORED flag; it keeps controlling whether the termination test just met
                    // is the STRICTER composite one (||K_T.V|| AND the TolF/TolU equilibrium
                    // residuals) or ||K_T.V|| alone -- a RUN-LEVEL property of the solver
                    // option, not a per-row one. A consumer that needs the stronger certificate
                    // must read the caller's own SingularPointComposite setting; the
                    // replacement filter for "give me a certified singular point" is
                    // isBifurcation == true && stability == 0 (this point), and "give me only
                    // the ones a solve actually resolved" is !unresolved -- see
                    // gsALMLandscape<T>::Point::equilibrium / Point::unresolved.
                    //
                    // Only gsALMBase's "SingularPointComposite" adds the equilibrium residuals
                    // to the termination test. It is READ, not overridden: the singular-point
                    // regime belongs to the caller (see the constructor doxygen), and forcing
                    // it on has been MEASURED to push the extended Newton -- an inexact Newton
                    // on a singular tangent -- past the point where it stops improving and into
                    // divergence (L -> 1e208 on that same crossing).
                    //
                    // @note This reads the OPTION LIST, while _extendedSystemSolve branches on
                    //       the solver's cached m_SPComposite, refreshed only by applyOptions()
                    //       / setOptions(). A caller that sets the switch and never applies it
                    //       therefore gets the WEAKER (||K_T.V||-only) termination test even
                    //       though it asked for the composite one -- silently, because there is
                    //       no public accessor for the applied state -- while the stored point
                    //       is still (correctly, under the new semantics) equilibrium = true:
                    //       it met the test that actually ran. The contract is unchanged: apply
                    //       your options before handing the solver over if you want the
                    //       stronger test to be the one that ran.
                    const bool spCertified =
                        m_solver->options().getSwitch("SingularPointComposite");

                    // testPoint=false: the classification just ran, above, on this very
                    // point and this very tangent (isBifurcation(true)), and it left m_V
                    // as the extended solve's initial null vector. jacobian is INERT when
                    // testPoint==false (gsALMBase.hpp:524: the tangent is only
                    // re-assembled inside the testPoint branch), so passing false here
                    // reuses that same tangent rather than re-assembling it for nothing.
                    gsStatus spStatus = gsStatus::NotConverged;
                    try {
                        spStatus = m_solver->computeSingularPoint(Uloc, Lloc,
                                                                  /*switchBranch=*/false,
                                                                  /*jacobian=*/false,
                                                                  /*testPoint=*/false);
                    } catch (...) { spStatus = gsStatus::OtherError; }

                    if (spStatus != gsStatus::Success)
                    {
                        // Honest-failure path (distinct from the localization/classification
                        // failures above): singular point solve did not converge. The
                        // extended (U*,L*) is not available, so fall back to marking the
                        // post-crossing traced point as the approximate bifurcation -- flagged
                        // unresolved, never indistinguishable from a certified one.
                        m_landscape.markUnresolvedSingular(cid);
                        gsWarn << "gsALMExploration: possible singular point near lambda = " << Lcur
                               << " on curve " << cid << " -- NOT refined (singular-point solve "
                                  "returned status " << (int)spStatus << "); no branch followed. "
                                  "The traced point is stored with unresolved = true and "
                                  "equilibrium = false; it is NOT a certified bifurcation.\n";
                        if (verbose)
                            gsInfo << "  [curve " << cid << "] singular point solve did NOT converge; "
                                   << "marking post-crossing point as unresolved; no branch jobs.\n";
                    }
                    else
                    {
                        const gsVector<T> Ustar = m_solver->solutionU();
                        const T           Lstar = m_solver->solutionL();
                        const gsVector<T> V     = m_solver->solutionV();
                        const T           tau   = m_solver->options().getReal("Perturbation");
                        const T           vnorm = V.norm();

                        // The converged extended-solve singular
                        // point IS the branch point. Store it as a flagged landscape
                        // point (stability neutral: it is the crossing itself) and use
                        // its index as the parent of the branch jobs. Do NOT mark the
                        // post-crossing traced point in this case.
                        gsMultiPatch<T> spDeformed;
                        bool spHaveGeometry = false;
                        if (m_solutionConstructor)
                            spHaveGeometry = m_solutionConstructor(Ustar, spDeformed);
                        // The extended-system solve produces no certified inertia count at
                        // (Ustar,Lstar) -- the tangent there is singular by construction --
                        // so the point records the "not recorded" sentinel rather than the
                        // pre- or post-crossing count.
                        //
                        // equilibrium = true unconditionally, NOT
                        // spCertified. spStatus == Success means the extended solve met ITS
                        // OWN termination test (see the provenance note above); that is the
                        // new definition of equilibrium. SingularPointComposite keeps
                        // controlling whether that test was the stricter composite one, but
                        // no longer decides the stored flag -- it stays in use below only to
                        // pick the right gsWarn.
                        m_landscape.addPoint(cid, Ustar, Lstar, /*stability=*/0,
                                             spHaveGeometry ? &spDeformed : nullptr,
                                             /*isBifurcation=*/true,
                                             /*equilibrium=*/true, /*negatives=*/-1);
                        if (!spCertified)
                            gsWarn << "gsALMExploration: the singular point at lambda = " << Lstar
                                   << " was accepted on ||K_T.V|| alone; its termination test did "
                                      "NOT include the equilibrium residuals. Enable the solver "
                                      "option \"SingularPointComposite\" if downstream use "
                                      "requires that stronger certificate.\n";
                        const index_t spPointIdx =
                            static_cast<index_t>(m_landscape.curve(cid).points.size()) - 1;

                        // --- pde2path swibra ABE emanating tangent ------------
                        // tau0 is the ACCEPTED-point secant (Ucur-Uold, Lcur-Lold): Uold/Lold
                        // is traceSweep's pre-crossing memory, a plain loop local untouched by
                        // the computeSingularPoint() call above (that call only moves the
                        // SOLVER's internal state -- m_U/m_L/m_V/... -- not these locals), so
                        // it remains the genuine incoming branch direction regardless of when
                        // this call happens relative to computeSingularPoint. Do NOT use
                        // m_solver->solutionDU()/solutionDL() here: by this point they hold the
                        // EXTENDED SOLVE's own increment, not the pre-crossing secant.
                        // negPrevAccepted != -1 (required for `detect`, hence for reaching this
                        // block at all) guarantees Uold/Lold really is a prior accepted point --
                        // on a sweep's very first crossing that is the sweep's own start state,
                        // still a valid incoming tangent.
                        gsVector<T> tau1U; T tau1L = 0;
                        const typename gsALMBase<T>::branchTangent::type bt =
                            m_solver->computeBranchTangent(Ustar, Lstar, V,
                                                           Ucur - Uold, Lcur - Lold,
                                                           tau1U, tau1L);

                        // Diagnostic: the explorer-side half of the
                        // discriminating measurement -- everything computeBranchTangent's
                        // out-parameters and return code make available at THIS call site,
                        // gated on the explorer's own `verbose` local (already ON in the
                        // ModifiedBratu drivers), independent of the solver's own m_verbose
                        // (which the drivers never set). See gsALMBase.hpp's
                        // computeBranchTangent for the internal al1/a1/b1/al1b/lam0 line,
                        // which needs the solver's Verbose to fire.
                        if (verbose)
                            gsInfo << "  [curve " << cid << "] computeBranchTangent at lambda="
                                   << Lstar << ": bt=" << (int)bt << ", tau1L=" << tau1L
                                   << ", ||tau1U||=" << tau1U.norm() << ".\n";

                        if (bt == gsALMBase<T>::branchTangent::Success ||
                            bt == gsALMBase<T>::branchTangent::TrivialBranch)
                        {
                            // TANGENT path (swibra ABE resolved). A design of "ONE job regardless of
                            // BranchPoints, +tau1/-tau1 are the two sweeps of that one job"
                            // rests on a premise that is FALSE for gsALMLoadControl and
                            // MEASURED to be so: traceSweep ties the predictor OFFSET sign to
                            // the sweep's ARC-LENGTH direction (backward), but LoadControl's own
                            // predictor() never reads the installed secant beyond that first
                            // step -- it marches lambda by the raw SIGN of `len` (== the
                            // arc-length direction) for every subsequent step. At a
                            // (near-)symmetric bifurcation (alpha1 ~ 0) BOTH mode signs need
                            // lambda to INCREASE to reach the real branch (fixture P's closed
                            // form: u2^2 = 2(lambda-1), valid at lambda>1 for EITHER sign of
                            // u2), so tying the sign to the arc-length direction makes exactly
                            // ONE of the two signs reachable; the other's sweep decreases lambda
                            // and collapses back onto the trivial solution before it ever
                            // reaches the branch. MEASURED: unittest fixture
                            // symmetric_branches_or_dedup (gsALMLoadControl, an internal
                            // scratch check) -- one child curve, u2 range clamped to [0,+X],
                            // haveNeg never true, independent of BranchEscape.
                            //
                            // Fix: mirror the FALLBACK path's structure -- one job PER SIGN
                            // (BranchPoints signs, matching the fallback path's own loop and its
                            // documented driver default of 1 = "+ only"), each job's tangent
                            // FIXED to sign*tau1. traceSweep is untouched: within one job, its
                            // FORWARD sweep (backward=false, s=+1) carries +sign*tau1 as the
                            // predictor offset with lambda INCREASING -- the "correct-side"
                            // sweep for that sign -- and its BACKWARD sweep carries -sign*tau1
                            // with lambda decreasing -- the "wrong-side" sweep, expected to
                            // fail/retrace, exactly as the fallback path's own wrong-side sweep
                            // is cheap by comparison (see its comment below). At the library
                            // default BranchPoints=2 this reproduces the OLD two-job dedup
                            // behaviour for a (near-)symmetric bifurcation; at BranchPoints=1
                            // (both example drivers' own default) this spawns exactly the SAME
                            // single job (sign=+1) as before -- MEASURED bit-identical on
                            // example_ModifiedBratuExploration (see the task report's A/B).
                            //
                            // Dedup key: normalised sign*tau1U, or sign*Vhat if tau1U happens to
                            // have zero norm (unreachable in practice -- computeBranchTangent
                            // normalises tau1U/tau1L by distance() before returning Success, and
                            // TrivialBranch sets tau1U := phi1, which step 1 there already
                            // guarded to have positive norm -- kept defensively).
                            const T tau1norm = tau1U.norm();
                            for (index_t s = 0; s < branchPoints; ++s)
                            {
                                const T sign = (s == 0) ? (T)(+1) : (T)(-1);
                                const gsVector<T> signedTauU = sign * tau1U;
                                const T           signedTauL = sign * tau1L;

                                gsVector<T> dir;
                                if (tau1norm > 0)
                                    dir = signedTauU / tau1norm;
                                else if (vnorm > 0)
                                    dir = sign * V / vnorm;
                                else
                                    dir = signedTauU;

                                if (isDuplicateJob(Lstar, dir))
                                {
                                    if (verbose)
                                        gsInfo << "  [curve " << cid << "] duplicate branch job "
                                                  "(tangent path) discarded.\n";
                                    continue;
                                }

                                Job nj;
                                nj.U              = Ustar;
                                nj.L              = Lstar;
                                // nj.Unudged left EMPTY: the tangent path seeds from
                                // nj.tangentU/tangentL (see traceSweep), not a nudged state.
                                nj.tangentU       = signedTauU;
                                nj.tangentL       = signedTauL;
                                nj.backward       = false;  // sweep 1 = +sign*tau1, sweep 2 = -sign*tau1
                                nj.bothDirections = true;
                                nj.parentCurve    = cid;
                                nj.parentPointIdx = spPointIdx; // provisional; remapped on flush
                                pending.push_back( give(nj) );
                                pendingPt.push_back(spPointIdx);

                                JobKey key;
                                key.L = Lstar; key.dir = dir;
                                m_jobKeys.push_back( give(key) );
                                spawnedJobs = true;

                                if (verbose)
                                    gsInfo << "  [curve " << cid << "] queued branch job at lambda=" << Lstar
                                           << " (tangent path, sign " << sign << ", outcome "
                                           << (bt == gsALMBase<T>::branchTangent::Success
                                               ? "Success" : "TrivialBranch")
                                           << ", both directions along +-sign*tau1).\n";
                            }
                        }
                        else
                        {
                            // FALLBACK-NUDGE path: the ABE emanating tangent could not be
                            // resolved (ModeUnresolved / Degenerate / AssemblyFailed). Fall back
                            // to switchBranch()'s nudge predictor U*+-V/Perturbation, unchanged.
                            const char * btName =
                                (bt == gsALMBase<T>::branchTangent::ModeUnresolved) ? "ModeUnresolved" :
                                (bt == gsALMBase<T>::branchTangent::Degenerate)     ? "Degenerate"     :
                                (bt == gsALMBase<T>::branchTangent::AssemblyFailed) ? "AssemblyFailed" :
                                                                                       "Unknown";
                            gsWarn << "gsALMExploration: computeBranchTangent could not resolve "
                                      "the emanating branch tangent at lambda=" << Lstar
                                   << " (outcome " << btName << "); falling back to the nudge "
                                      "predictor U* +- V/Perturbation.\n";

                            if (vnorm > 0)
                            {
                                const gsVector<T> Vhat = V / vnorm;
                                // ONE job per mode sign; that job sweeps BOTH arc-length
                                // directions into one curve (see traceCurve).
                                //
                                // The mode sign selects WHICH of the two symmetric halves of the
                                // emanating branch is followed; it says nothing about the side of
                                // Lstar on which that branch LIVES. Tracing forward only
                                // would therefore trace away from the branch the job
                                // was created for whenever the bifurcation is subcritical -- e.g.
                                // the modified-Bratu A-B pitchfork, whose curve B exists only at
                                // mu < mu*, so every branch job was doomed by construction.
                                //
                                // Both directions are swept rather than one being CHOSEN from the
                                // emanating tangent: on THIS path the ABE above could not resolve
                                // that tangent in the first place (that is exactly why we are
                                // here), so there is no first-order sign to read either. The
                                // wrong-side SWEEP is cheap by comparison: it either fails to
                                // converge (arc-length underflow ends it) or falls back onto the
                                // parent branch, where the C_start test rewinds it -- and because
                                // it is a sweep and not a curve, rewinding it cannot cost the
                                // genuine half that the other sweep found.
                                //
                                // BACKWARD IS SWEPT FIRST: it is the direction a forward-only sweep
                                // could never take.
                                //
                                // OPTIONALLY floor the absolute nudge 1/Perturbation
                                // at a magnitude derived from the SAME normalization isRetrace's
                                // own ball uses (retraceThreshold() * max(1,||U_p||), see
                                // retraceDistance()) evaluated at THIS branch point's own state
                                // Ustar. 1/Perturbation is tuned once, driver-wide, at whichever
                                // locus the caller calibrated it on -- example_ModifiedBratu-
                                // Exploration's own header comment documents that calibration at
                                // the (low-state-norm) A-B point, nudge 1.0 clearing the retrace
                                // ball there with a 2.3x margin. The retrace ball's radius SCALES
                                // with ||U_p|| while a fixed absolute nudge does not, so the same
                                // margin does not survive to a HIGHER-state-norm locus: MEASURED
                                // on this driver's own A-C point (mu=0.03507, ||Ustar||=28.8504),
                                // the unfloored nudge 1.0 gives a retrace-predicate ratio of
                                // 0.0626 against threshold 0.0693 (factor 0.90 -- the genuine
                                // branch is one arc-length sweep away from being wrongly captured
                                // as a retrace and IS captured on this problem, dropping the whole
                                // curve), while 4.0 gives 0.1625 (factor 2.35, clear).
                                //
                                // A hardcoded margin here
                                // was calibrated as an EQUALITY at exactly ONE driver's ONE locus
                                // (example_ModifiedBratuExploration's A-B point, margin 1.2458 =
                                // 1.0/(0.069282*11.5858)) and applied UNCONDITIONALLY to every
                                // caller. At the LIBRARY default Perturbation=1e3 (gsALMBase.hpp,
                                // nudgeBase=1e-3) that margin floors the nudge to
                                // 1.2458*retraceThreshold()*max(1,||Ustar||), which for an ordinary
                                // PDE state norm is ~1e2-1e4x the documented 1/Perturbation scale
                                // -- i.e. the floor would be ACTIVE, and Perturbation effectively
                                // INERT, for essentially every caller that does not override it,
                                // with no way to opt out. RetraceBallFloor (defaultOptions(),
                                // <= 0 = disabled) fixes this: the floor is now driver-selected,
                                // not library-imposed, and no library-derived margin is asserted --
                                // see the option's own doc for why. margin <= 0 makes this block
                                // bit-identical to the unfloored nudge (nudgeMag == nudgeBase
                                // always, no gsInfo line, Perturbation never overridden).
                                //
                                // max(...) is a FLOOR, never a ceiling: any locus whose EXISTING
                                // 1/Perturbation nudge already clears the ball at the caller's
                                // margin is left BIT-IDENTICAL (nudgeBase wins the max); only a
                                // locus where it does not is affected. Perturbation therefore keeps
                                // its documented inverse-scale role wherever it already suffices,
                                // and this is a floor on top of it, not a replacement for it.
                                const T ballFloorMargin = m_options.getReal("RetraceBallFloor");
                                const T nudgeBase  = (T)1 / tau;
                                const T ballRadius = retraceThreshold() * math::max((T)1, Ustar.norm());
                                const T nudgeMag   = (ballFloorMargin > (T)0)
                                                    ? math::max(nudgeBase, ballFloorMargin * ballRadius)
                                                    : nudgeBase;
                                if (verbose && nudgeMag > nudgeBase)
                                    gsInfo << "  [curve " << cid << "] fallback-nudge magnitude "
                                              "floored from " << nudgeBase << " (1/Perturbation) to "
                                           << nudgeMag << " (" << ballFloorMargin
                                           << " * retraceThreshold() * max(1,||Ustar||) = "
                                           << ballFloorMargin << " * " << retraceThreshold() << " * "
                                           << math::max((T)1, Ustar.norm()) << ") at lambda=" << Lstar
                                           << ".\n";
                                for (index_t s = 0; s < branchPoints; ++s)
                                {
                                    const T sign = (s == 0) ? (T)(+1) : (T)(-1);
                                    const gsVector<T> dir     = sign * Vhat;       // signed nudge direction
                                    const gsVector<T> Unudged = Ustar + nudgeMag * dir;

                                    if (isDuplicateJob(Lstar, dir))
                                    {
                                        if (verbose)
                                            gsInfo << "  [curve " << cid << "] duplicate branch job discarded.\n";
                                        continue;
                                    }

                                    Job nj;
                                    nj.U              = Ustar;
                                    nj.L              = Lstar;
                                    nj.Unudged        = Unudged;
                                    nj.backward       = true;   // backward swept first
                                    nj.bothDirections = true;
                                    nj.parentCurve    = cid;
                                    nj.parentPointIdx = spPointIdx; // provisional; remapped on flush
                                    pending.push_back( give(nj) );
                                    pendingPt.push_back(spPointIdx);

                                    JobKey key;
                                    key.L = Lstar; key.dir = dir;
                                    m_jobKeys.push_back( give(key) );
                                    spawnedJobs = true;

                                    if (verbose)
                                        gsInfo << "  [curve " << cid << "] queued branch job at lambda=" << Lstar
                                               << " (sign " << sign << ", both directions).\n";
                                }
                            }
                        }
                    }
                }
            }

            // Landscape bookkeeping for the multi-step-crossing reduction warned about
            // above, placed here (after all five storage paths, before the solver-state
            // restore below) because every one of them writes or relabels exactly the
            // LAST point of curve cid.
            if (multReduced)
            {
                typename gsALMLandscape<T>::Curve & c = m_landscape.curve(cid);
                GISMO_ASSERT(!c.points.empty(),
                             "gsALMExploration: curve "<<cid<<" has no points to flag as "
                             "multiplicity-reduced.");
                c.points.back().multiplicityReduced = true;
            }

            // --- ONE restore, reached by EVERY exit path out of this block --------
            // (success, limit, singular-solve failure, classification failure and
            // localization failure): this restore is load-bearing on every path, not
            // just the success path, because the bisection probes in
            // _localizeCrossing and the classification/solve calls above all move the
            // solver (m_U, m_L, m_Uprev, the arc length, the tangent and the
            // stability data), not just the success path.
            //
            // ORDER matters: setPrevious derives m_DeltaUold/m_DeltaLold from the
            // CURRENT m_U/m_L and copies m_arcLength into m_arcLength_prev, so it
            // must run AFTER setSolution/setLength restore the current point.
            m_solver->setSolution(Ucur, Lcur);
            m_solver->setLength(dLb);
            // setPrevious(Uold,Lold) restores EXACTLY the value the accepted step
            // itself left (every solver's iterationFinish() sets m_Uprev=m_U,
            // m_Lprev=m_L BEFORE m_U+=m_DeltaU) -- the truthful predictor history,
            // which the localization probes and the base's own bisection stage would
            // otherwise leave pointing at a probe state (see the seeding block
            // above: a stale m_Uprev degenerates the secant predictor).
            m_solver->setPrevious(Uold, Lold);
            // Re-sync the stability memory to the restored post-crossing side:
            // computeSingularPoint / the probes leave m_stability synced to whatever
            // they last touched, so without this the next step's stabilityChange()
            // would spuriously re-detect the SAME crossing (inserting a duplicate
            // singular point). Use the RELIABLE post-crossing indicator captured
            // before any of this ran -- NOT a recompute at a near-singular point,
            // whose sign is numerically unreliable. Trivially revertible.
            m_solver->setIndicator(curIndicator);
            // m_deltaU is NOT restored here: it is polluted by the probes exactly as
            // it is already polluted by the corrector, and gsALMBase.hpp's own note
            // in _computeSingularPoint records that staleness as pre-existing and out
            // of scope for this fix -- unchanged in kind by this restore.
        }

        // Advance the pre-crossing memory and restore the base arc length.
        Uold = Ucur;
        Lold = Lcur;
        stabPrevAccepted = stab;    // track last accepted stored stability sign
        // Track last accepted tangent inertia -- but NOT while a flip is being
        // suppressed by the branch-curve guard: holding the
        // reference here is what lets a PERSISTING flip be caught at the first
        // unsuppressed step, instead of being silently absorbed as the new baseline.
        if (!suppressed)
            negPrevAccepted = negCur;

        if (stepGrowth)
        {
            // pde2path sscontrol.m:39, p.sol.ds = ds*p.nc.dsincfac with dsincfac = 2,
            // fired on a "very good step" (iter < dsinciter = imax/2, stanparam.m:50).
            // The clamp is on |dLb|: dLb < 0 on a backward sweep and dLb0 is the
            // SIGNED base length, so a magnitude comparison is the only correct one.
            const T dLbPrev = dLb;
            if (itersCur <= growIter)
                dLb *= (T)2;
            if (math::abs(dLb) > math::abs(dLb0))
                dLb = dLb0;                      // ceiling = Length / SwitchLength
            m_solver->setLength(dLb);
            if (verbose && dLb != dLbPrev)       // only when dLb actually changed
                gsInfo << "  [curve " << cid << "] step converged in " << itersCur
                       << " iterations; growing arc length to " << dLb << ".\n";
        }
        else if (!bisected)
        {
            dLb = dLb0;
            m_solver->setLength(dLb);
        }
        bisected = false;
    }

    m_sweepTermination.push_back(std::make_pair(cid, reason));
    return true;
}

template <class T>
bool gsALMExploration<T>::_localizeCrossing(const gsVector<T> & Uold, T Lold, index_t negOld,
                                            const gsVector<T> & Ucur, T Lcur, index_t negCur,
                                            T dLb, T dLb0, index_t cid,
                                            gsVector<T> & Uloc, T & Lloc, T & bracket,
                                            index_t * probesOut)
{
    const bool    verbose  = m_options.getSwitch("Verbose");
    const index_t bisecMax = m_options.getInt("BisecMax");
    // Relative floor: an ABSOLUTE ds would be a moving target under a
    // state-dependent arc-length metric (Scaling=-1), see the option doxygen.
    const T       floor    = m_options.getReal("BisecLengthFloor") * math::abs(dLb0);

    // Bracket on the arc-length parameter s in [0,dLb], measured from (Uold,Lold):
    // (slo,negLo) is the pre-crossing side, (shi,negHi) the post-crossing side.
    // Invariant: negLo != negHi. The probe is the position currently being
    // stepped to; it is NOT one of the two bracket endpoints.
    T           slo   = (T)0, shi = dLb;
    index_t     negLo = negOld, negHi = negCur;
    gsVector<T> Ulo   = Uold;
    T           Llo   = Lold;

    index_t      probes  = 0;
    bool         refined = false; // true once a single probe has certified a half-interval
    const char * reason  = "no probe converged";
    T            sprobe  = (slo + shi) / (T)2;

    while (true)
    {
        // Termination, checked BEFORE spending a probe: BisecMax exhausted, or the
        // bracket already at/below the floor. Both are NORMAL terminations.
        if (probes >= bisecMax)            { reason = "BisecMax";         break; }
        if (math::abs(shi - slo) <= floor) { reason = "BisecLengthFloor"; break; }

        // Each probe: re-seed from the PRE-crossing accepted point and step to the
        // current midpoint, exactly like the sweep's own step-fail retry.
        m_solver->setSolution(Uold, Lold);
        m_solver->setLength(sprobe);
        gsStatus status;
        try { status = m_solver->step(); }
        catch (...) { status = gsStatus::AssemblyError; } // NaN metric on a bad predictor
        ++probes;

        if (status == gsStatus::Success)
        {
            // Converged probe: fresh tangent and inertia at the probe, exactly the
            // call the sweep makes on an accepted point.
            m_solver->computeStability(true);
            const index_t     negProbe = m_solver->negatives();
            const gsVector<T> Uprobe   = m_solver->solutionU();
            const T           Lprobe   = m_solver->solutionL();

            if (negProbe == negLo)
            {
                slo = sprobe; Ulo = Uprobe; Llo = Lprobe; negLo = negProbe;
            }
            else
            {
                // Keeps the invariant negLo != negHi when a multi-step crossing
                // (e.g. 1->3) is probed at an intermediate inertia (e.g. 2).
                shi = sprobe; negHi = negProbe;
            }
            refined = true;
            sprobe = (slo + shi) / (T)2;
        }
        else
        {
            // Non-converged probe: a probe with no inertia cannot certify a
            // half-interval, so neither bracket endpoint moves. Retreat the probe
            // toward the low end instead (the same halving the sweep itself uses
            // on a failed step) and count the attempt against BisecMax.
            sprobe = (slo + sprobe) / (T)2;
            if (math::abs(sprobe - slo) < floor)
            {
                reason = refined ? "BisecLengthFloor" : "no probe converged";
                break;
            }
        }
    }
    GISMO_UNUSED(negHi); // tracked only to keep the documented bracket invariant

    bracket = math::abs(shi - slo);
    Uloc    = Ulo;
    Lloc    = Llo;
    if (probesOut != nullptr) *probesOut = probes;

    if (verbose)
        gsInfo << "  [curve " << cid << "] localized crossing: interval lambda in ["
               << Lold << "," << Lcur << "], " << probes << " probe(s), stopped on "
               << reason << ", bracket width = " << bracket
               << ", localized lambda = " << Lloc << ".\n";

    return refined;
}

template <class T>
gsStatus gsALMExploration<T>::solve(const gsVector<T> & U0, T L0)
{
    // Delegate to the multi-seed overload with a single seed (bit-identical to
    // the previous single-seed implementation).
    std::vector<std::pair<gsVector<T>,T> > seeds(1, std::make_pair(U0, L0));
    return solve(seeds);
}

template <class T>
gsStatus gsALMExploration<T>::solve(const std::vector<std::pair<gsVector<T>,T> > & seeds)
{
    m_landscape = gsALMLandscape<T>();
    m_jobKeys.clear();
    m_sweepTermination.clear();
    m_retraceTests = 0;
    m_retraceFired = 0;

    const index_t maxCurves = m_options.getInt("MaxCurves");
    const bool    verbose   = m_options.getSwitch("Verbose");
    const std::string prefix = m_options.getString("OutputPrefix");

    std::queue<Job> queue;

    // ONE seed job per seed (U0,L0), parent -1, tracing FORWARD first and then
    // BACKWARD into the same curve -- the two halves of a seed state are two ends
    // of ONE locus, not two curves (thesis Alg. 8.1 lines 17-23; see traceCurve).
    // At a pristine rest state (U0=0, L0=0) only the forward direction is swept:
    // there is nothing behind the undeformed configuration, so a backward sweep
    // only re-explores the forward branch mirrored. That is the same condition
    // that suppresses the backward JOB here, so a rest-state seed yields
    // exactly the points a forward-only sweep would produce.
    for (size_t i = 0; i != seeds.size(); ++i)
    {
        const gsVector<T> & U0 = seeds[i].first;
        const T             L0 = seeds[i].second;

        const bool restState = (math::abs(L0) < 1e-14) && (U0.size() == 0 || U0.norm() == 0);

        Job sj;
        sj.U = U0; sj.L = L0; sj.backward = false;
        sj.bothDirections = !restState;
        sj.parentCurve = -1; sj.parentPointIdx = -1;
        sj.tangentL = 0; // tangentU stays default-empty: a seed job has no tangent
        queue.push(give(sj));
    }

    while (!queue.empty() &&
           static_cast<index_t>(m_landscape.nCurves()) < maxCurves)
    {
        Job job = queue.front();
        queue.pop();

        traceCurve(job, queue);

        // Crash checkpoint: cheap re-write after every completed curve.
        if (!prefix.empty())
        {
            m_landscape.writeCsv(prefix + ".csv");
#ifdef gsHDF5_ENABLED
            // Whole-file HDF5 checkpoint. A storage failure must not
            // kill a long exploration run, so swallow it with a warning here;
            // the direct saveHDF5 API stays throwing.
            try { m_landscape.saveHDF5(prefix + ".h5"); }
            catch (...)
            {
                gsWarn << "gsALMExploration: HDF5 checkpoint to " << prefix
                       << ".h5 failed; continuing.\n";
            }
#endif
        }
    }

    if (verbose)
        gsInfo << "gsALMExploration: " << m_landscape.nCurves() << " curves, "
               << m_landscape.nPoints() << " points.\n";

    // Unconditional (under Verbose) report of the isRetrace safety-net rate,
    // printed EXACTLY ONCE per solve() call, including the M=0 case -- a 0% rate is an
    // expected, acceptable result (the tangent path is expected to fire it ~never; see
    // the class doxygen).
    if (verbose)
        gsInfo << "gsALMExploration: retrace safety net (isRetrace) fired "
               << m_retraceFired << " of " << m_retraceTests << " evaluations ("
               << (m_retraceTests > 0 ? (100.0*m_retraceFired)/m_retraceTests : 0.0)
               << "%).\n";

    // The only hard failure is an empty landscape: the first curve could not take
    // a single step. All other curve-level failures are logged and tolerated.
    if (m_landscape.nPoints() == 0)
        return gsStatus::NotConverged;

    return gsStatus::Success;
}


} // namespace gismo
