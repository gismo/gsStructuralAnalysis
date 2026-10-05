/** @file gsALMSolver_test.cpp

    @brief Regression recordings of the four arc-length presets (load control,
    Riks, Crisfield, consistent Crisfield) on the two-dof von Mises truss of
    gsStage0Truss.h.

    * Reference gates: for each preset and each combination of AdaptiveLength
      and Quasi, the per-step (lambda, |U|, status, iterations) trajectory is
      compared with a recording embedded below. Steps whose iteration count is
      not robust against a change of the convergence tolerance are marked as
      boundary steps; their iteration count is not gated.
    * Negative controls: load control cannot pass the limit point, Riks and
      Crisfield do.

    Setting the environment variable GSSA_STAGE0_PRINT prints all reference
    arrays as C++ initialisers (through gsInfo) and skips the comparison with
    the embedded arrays.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
 **/

#include "gismo_unittest.h"

#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <string>
#include <vector>

#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMLoadControl.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMRiks.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMCrisfield.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMConsistentCrisfield.h>
#include "gsStage0Truss.h"

SUITE(gsALMSolver_test)
{
    /// One recorded step: load factor, norm of the solution, status, number of
    /// iterations and whether the iteration count is tolerance-sensitive
    struct T1Step
    {
        real_t     L;
        real_t     normU;
        gsStatus   status;
        index_t    iterations;
        bool       boundary;
    };

    /// A recorded step together with the apex deflection w = U[0] and the
    /// load factor the step aimed at (previous lambda plus the arc length;
    /// for load control the target of the load increment)
    struct T1Rec
    {
        T1Step step;
        real_t w;
        real_t target;
    };

    /// Base convergence tolerance (TolF = TolU). In a non-adaptive trajectory a
    /// step whose iteration count flips moves lambda and |U| by about one last
    /// Newton correction, which stays below the value gate of 1e-8. With
    /// AdaptiveLength the iteration count sets the next arc length, so a flipped
    /// count changes all later steps by far more than 1e-8.
    static const real_t T1_tol = 1e-10;

    /// Arc-length parameters and number of steps of one trajectory
    struct T1Params
    {
        real_t  ds;
        index_t N;
    };

    enum T1Preset { LoadControl = 0, Riks = 1, Crisfield = 2, ConsistentCrisfield = 3 };

    /// Step size and step count of trajectory (preset, adaptive, quasi).
    /// Load control uses ds as load increment; the other presets as arc length.
    static T1Params t1Params(const T1Preset preset, const bool adaptive, const bool quasi)
    {
        GISMO_UNUSED(adaptive); GISMO_UNUSED(quasi);
        T1Params p;
        if (preset == LoadControl) { p.ds = 0.07; p.N = 20; }
        else                       { p.ds = 0.05; p.N = 40; }
        return p;
    }

    /// Truss, its operators and a solver of type S, configured with
    /// the options of one trajectory
    template<class S>
    struct T1Setup
    {
        typedef typename gsStructuralAnalysisOps<real_t>::Jacobian_t   Jacobian_t;
        typedef typename gsStructuralAnalysisOps<real_t>::ALResidual_t ALResidual_t;

        gsStage0Truss<real_t> truss;
        Jacobian_t            J;
        ALResidual_t          R;
        gsVector<real_t>      F;
        S                     solver;

        T1Setup(real_t ds, bool adaptive, bool quasi, real_t tolScale)
        : J(truss.jacobian()), R(truss.alResidual()), F(truss.force()), solver(J,R,F)
        {
            // setOptions() re-reads "Length", so all options are set at once
            gsOptionList opts = solver.options();
            opts.setReal("Length", ds);
            opts.setSwitch("AdaptiveLength", adaptive);
            opts.setInt("AdaptiveIterations", 4);
            opts.setSwitch("Quasi", quasi);
            opts.setReal("TolF", T1_tol*tolScale);
            opts.setReal("TolU", T1_tol*tolScale);
            solver.setOptions(opts);
            solver.initialize();
        }
    };

    /// Runs N steps (stopping after the first failed step) of a fresh solver
    template<class S>
    std::vector<T1Rec> t1Run(real_t ds, index_t N, bool adaptive, bool quasi, real_t tolScale = 1)
    {
        T1Setup<S> setup(ds, adaptive, quasi, tolScale);
        std::vector<T1Rec> rec;
        for (index_t k = 0; k < N; ++k)
        {
            const real_t target = setup.solver.solutionL() + setup.solver.getLength();
            const gsStatus st = setup.solver.step();
            T1Rec r;
            r.target          = target;
            r.step.L          = setup.solver.solutionL();
            r.step.normU      = setup.solver.solutionU().norm();
            r.step.status     = st;
            r.step.iterations = setup.solver.numIterations();
            r.step.boundary   = false;
            r.w               = setup.solver.solutionU()[0];
            rec.push_back(r);
            if (st != gsStatus::Success)
                break;
        }
        return rec;
    }

    /// Per-trajectory gating data derived from the tolerance-perturbed runs
    struct T1Meta
    {
        index_t gated;          ///< number of leading steps whose values and status are gated
        bool    fragile;        ///< the perturbed runs disagree on status or length from step \a gated on
        bool    gateIterations; ///< iteration counts of non-boundary steps are gated
    };

    /// Result of the three tolerance-perturbed runs of one trajectory
    struct T1Marked
    {
        std::vector<T1Rec> rec;   ///< the run at scale 1, boundary flags set
        T1Meta             meta;
        bool               targetPrefixEnd;   ///< the prefix was cut at the first load-control target near lambda_max
        bool               adaptivePrefixEnd; ///< the prefix was cut at the first boundary step (AdaptiveLength)
    };

    /// Runs a trajectory at the tolerance scales 0.7, 1 and 1.4. A step is a
    /// boundary step when its iteration count differs between the runs or the
    /// step is missing in one of them. The gated prefix ends at the first
    /// index where the runs disagree on status or on the existence of a step
    /// with AdaptiveLength at the first boundary step, and for load control at
    /// the first step whose target lies within 1e-2 of lambda_max. Iteration counts are gated only if at most 20 % of the gated steps are
    /// boundary steps.
    template<class S>
    T1Marked t1RunMarked(real_t ds, index_t N, bool adaptive, bool quasi, bool loadControl)
    {
        const std::vector<T1Rec> lo = t1Run<S>(ds, N, adaptive, quasi, 0.7);
        T1Marked m;
        m.adaptivePrefixEnd = false;
        m.targetPrefixEnd   = false;
        m.rec = t1Run<S>(ds, N, adaptive, quasi, 1.0);
        const std::vector<T1Rec> hi = t1Run<S>(ds, N, adaptive, quasi, 1.4);
        std::vector<T1Rec> & mid = m.rec;

        const size_t maxLen = std::max(lo.size(), std::max(mid.size(), hi.size()));
        size_t K = maxLen;
        for (size_t i = 0; i < maxLen; ++i)
        {
            const bool inLo = i < lo.size(), inMid = i < mid.size(), inHi = i < hi.size();
            if (!(inLo && inMid && inHi) || lo[i].step.status != mid[i].step.status
                || hi[i].step.status != mid[i].step.status)
            {
                K = i;
                break;
            }
        }
        for (size_t i = 0; i < mid.size(); ++i)
        {
            const bool inLo = i < lo.size(), inHi = i < hi.size();
            mid[i].step.boundary = !inLo || !inHi
                || lo[i].step.iterations != mid[i].step.iterations
                || hi[i].step.iterations != mid[i].step.iterations;
        }
        // With AdaptiveLength the iteration count sets the next arc length, so a
        // boundary step changes every later step: the gated prefix ends there.
        if (adaptive)
            for (size_t i = 0; i < K && i < mid.size(); ++i)
                if (mid[i].step.boundary)
                {
                    K = i;
                    m.adaptivePrefixEnd = true;
                    break;
                }
        // Whether a load-control step with a target close to lambda_max
        // converges, fails, or jumps to the far branch is decided by round-off.
        if (loadControl)
            for (size_t i = 0; i < K && i < mid.size(); ++i)
                if (math::abs(mid[i].target - gsStage0Truss<real_t>().lambdaMax()) <= 1e-2)
                {
                    K = i;
                    m.targetPrefixEnd = true;
                    break;
                }
        index_t b = 0;
        for (size_t i = 0; i < K && i < mid.size(); ++i)
            b += mid[i].step.boundary ? 1 : 0;
        m.meta.gated          = static_cast<index_t>(std::min(K, mid.size()));
        m.meta.fragile        = K < maxLen;
        m.meta.gateIterations = 5*b <= m.meta.gated;
        return m;
    }

    static bool t1PrintMode() { return std::getenv("GSSA_STAGE0_PRINT") != nullptr; }

    static const char * t1StatusName(const gsStatus s)
    {
        switch (s)
        {
        case gsStatus::Success:       return "gsStatus::Success";
        case gsStatus::NotConverged:  return "gsStatus::NotConverged";
        case gsStatus::AssemblyError: return "gsStatus::AssemblyError";
        case gsStatus::SolverError:   return "gsStatus::SolverError";
        case gsStatus::NotStarted:    return "gsStatus::NotStarted";
        default:                      return "gsStatus::OtherError";
        }
    }

    static index_t t1CountBoundary(const std::vector<T1Rec> & rec, size_t n)
    {
        index_t b = 0;
        for (size_t i = 0; i < rec.size() && i < n; ++i)
            b += rec[i].step.boundary ? 1 : 0;
        return b;
    }

    /// True when a successful step beyond the limit point has a lower load
    /// factor than its predecessor
    static bool t1PassedLimit(const std::vector<T1Rec> & rec, const real_t wLimit)
    {
        for (size_t i = 1; i < rec.size(); ++i)
            if (rec[i].step.status == gsStatus::Success && rec[i].w > wLimit
                && rec[i].step.L < rec[i-1].step.L)
                return true;
        return false;
    }

    /// Prints the recording and its gating data as ready-to-paste initialisers
    static void t1Print(const std::string & name, const T1Params & p, const bool adaptive,
                        const bool quasi, const T1Marked & m)
    {
        const std::vector<T1Rec> & rec = m.rec;
        char buf[512];
        const index_t b = t1CountBoundary(rec, static_cast<size_t>(m.meta.gated));
        std::snprintf(buf, sizeof(buf), "// %s: ds=%g N=%d adaptive=%d quasi=%d boundary %d / %d steps, limit point passed %d, final status %s at lambda %.6g\n",
                      name.c_str(), double(p.ds), int(p.N), int(adaptive), int(quasi),
                      int(b), int(m.meta.gated), int(t1PassedLimit(rec, gsStage0Truss<real_t>().wLimit())),
                      rec.empty() ? "none" : t1StatusName(rec.back().step.status),
                      rec.empty() ? 0.0 : double(rec.back().step.L));
        gsInfo << buf;
        if (m.meta.fragile)
            gsInfo << "// " << name << ": fragile from step " << m.meta.gated << "\n";
        if (m.targetPrefixEnd)
            gsInfo << "// " << name << ": load-control prefix ends at step " << m.meta.gated << "\n";
        if (m.adaptivePrefixEnd)
            gsInfo << "// " << name << ": adaptive prefix ends at step " << m.meta.gated << "\n";
        if (!m.meta.gateIterations)
        {
            std::snprintf(buf, sizeof(buf), "// %s: iterations not gated: %d / %d boundary\n",
                          name.c_str(), int(b), int(m.meta.gated));
            gsInfo << buf;
        }
        gsInfo << "static const std::vector<T1Step> " << name << " = {\n";
        for (size_t i = 0; i < rec.size(); ++i)
        {
            const T1Step & s = rec[i].step;
            std::snprintf(buf, sizeof(buf), "    { %.17g, %.17g, %s, %d, %s },\n",
                          double(s.L), double(s.normU), t1StatusName(s.status),
                          int(s.iterations), s.boundary ? "true" : "false");
            gsInfo << buf;
        }
        gsInfo << "};\n";
        gsInfo << "static const T1Meta " << name << "_meta = { " << m.meta.gated << ", "
               << (m.meta.fragile ? "true" : "false") << ", "
               << (m.meta.gateIterations ? "true" : "false") << " };\n";
    }

    /// Checks the recording of one trajectory against the reference \a ref
    /// and its gating data \a meta. The problem is nondimensional
    /// (lambda_max = 1, h = 1), so the value tolerance 1e-8 is relative to
    /// max(|ref|, 1). Only the first meta.gated steps are compared; the length
    /// of the recording is compared unless the trajectory is fragile. In print
    /// mode the recording is printed instead of compared; the checks on the
    /// fresh recording (every gated load-control
    /// step target stays clear of lambda_max, the
    /// non-vacuity of iteration gating of the plain cells) always run.
    template<class S>
    void t1Gate(const std::string & name, const std::vector<T1Step> & ref, const T1Meta & meta,
                const T1Preset preset, const bool adaptive, const bool quasi)
    {
        const T1Params p = t1Params(preset, adaptive, quasi);
        const gsStage0Truss<real_t> truss;
        const bool print = t1PrintMode();
        const bool plain = !adaptive && !quasi;

        std::vector<T1Rec> rec;
        T1Meta fresh = { 0, false, false };
        if (print)
        {
            const T1Marked m = t1RunMarked<S>(p.ds, p.N, adaptive, quasi, preset == LoadControl);
            t1Print(name, p, adaptive, quasi, m);
            rec = m.rec;
            fresh = m.meta;
        }
        else
            rec = t1Run<S>(p.ds, p.N, adaptive, quasi);

        if (preset == LoadControl)
        {
            const size_t n = std::min(rec.size(), static_cast<size_t>(print ? fresh.gated : meta.gated));
            for (size_t i = 0; i < n; ++i)
                CHECK(math::abs(rec[i].target - truss.lambdaMax()) > 1e-2);
        }

        if (print)
        {
            if (plain)
                CHECK(fresh.gateIterations);
            return;
        }

        CHECK(!ref.empty());
        if (plain)
            CHECK(meta.gateIterations);
        {
            index_t b = 0;
            for (size_t i = 0; i < static_cast<size_t>(meta.gated) && i < ref.size(); ++i)
                b += ref[i].boundary ? 1 : 0;
            CHECK_EQUAL(meta.gateIterations, 5*b <= meta.gated);
            if (adaptive)
                CHECK_EQUAL(0, b);
        }
        CHECK(meta.gated <= static_cast<index_t>(ref.size()));
        if (meta.fragile)
            CHECK(rec.size() >= static_cast<size_t>(meta.gated));
        else
            CHECK_EQUAL(ref.size(), rec.size());

        for (size_t i = 0; i < static_cast<size_t>(meta.gated) && i < ref.size() && i < rec.size(); ++i)
        {
            CHECK(rec[i].step.status == ref[i].status);
            CHECK_CLOSE(ref[i].L,     rec[i].step.L,     1e-8*math::max(math::abs(ref[i].L), real_t(1)));
            CHECK_CLOSE(ref[i].normU, rec[i].step.normU, 1e-8*math::max(math::abs(ref[i].normU), real_t(1)));
            if (meta.gateIterations && !ref[i].boundary)
                CHECK_EQUAL(ref[i].iterations, rec[i].step.iterations);
        }
    }

/* Reference recordings (16 trajectories). Provenance:
   generator command, run from the worktree root:
     bash /home/hverhelst/.claude/plugins/cache/gsagents/gismo/0.8.0/skills/evidence-log/scripts/run_logged.sh .claude/plans/gssa-redesign/tasks 01-r1d-generate env GSSA_STAGE0_PRINT=1 OMP_NUM_THREADS=1 bash /home/hverhelst/.claude/plugins/cache/gsagents/gismo/0.8.0/skills/run-tests/scripts/run_unittests.sh gsALMSolver_test
   generator log: .claude/plans/gssa-redesign/tasks/logs/063-01-r1d-generate.log, md5 67cbd4e834f56e02ad016a59b5551998
   gsStructuralAnalysis commit b5ec853; md5 of its uncommitted diff 8e48e0dff51a2811698d1c2a08f6406d
   date 2026-10-03
   truss: c = 3 sqrt(3)/2, h = 1, k = 5; TolF = TolU = 1e-10; AdaptiveIterations = 4;
   all other options at their defaults.
   Step size ds and number of steps N per trajectory: load control ds = 0.07 (load
   increment), N = 20; Riks, Crisfield, consistent Crisfield ds = 0.05, N = 40;
   a trajectory stops after its first failed step.
   Each array is followed by its T1Meta {gated steps, fragile, iterations gated}.
   With AdaptiveLength the gated prefix ends at the first boundary step; for load
   control it ends at the first step whose target is within 1e-2 of lambda_max.
   A boundary step has an iteration count that differs between the runs at tolerance
   scales 0.7, 1 and 1.4. */
// T1_LoadControl_A0_Q0: ds=0.07 N=20 adaptive=0 quasi=0 boundary 0 / 20 steps, limit point passed 0, final status gsStatus::Success at lambda 1.4
static const std::vector<T1Step> T1_LoadControl_A0_Q0 = {
    { 0.070000000000000007, 0.030975053280708726, gsStatus::Success, 3, false },
    { 0.14000000000000001, 0.062767990466369133, gsStatus::Success, 3, false },
    { 0.21000000000000002, 0.095485647222139508, gsStatus::Success, 3, false },
    { 0.28000000000000003, 0.12925955248539236, gsStatus::Success, 3, false },
    { 0.35000000000000003, 0.16425465591286578, gsStatus::Success, 3, false },
    { 0.42000000000000004, 0.20068245427560297, gsStatus::Success, 3, false },
    { 0.49000000000000005, 0.23882156737663307, gsStatus::Success, 3, false },
    { 0.56000000000000005, 0.27905167803791892, gsStatus::Success, 4, false },
    { 0.63000000000000012, 0.32191324971498469, gsStatus::Success, 4, false },
    { 0.70000000000000018, 0.36822188968365632, gsStatus::Success, 4, false },
    { 0.77000000000000024, 0.41931450588763941, gsStatus::Success, 4, false },
    { 0.8400000000000003, 0.47767975029250997, gsStatus::Success, 4, false },
    { 0.91000000000000036, 0.54913022011864188, gsStatus::Success, 4, false },
    { 0.98000000000000043, 0.65851078197277269, gsStatus::Success, 5, false },
    { 1.0500000000000005, 3.2081433983711385, gsStatus::Success, 17, false },
    { 1.1200000000000006, 3.2308812754676088, gsStatus::Success, 3, false },
    { 1.1900000000000006, 3.253402737551506, gsStatus::Success, 3, false },
    { 1.2600000000000007, 3.2757192492134988, gsStatus::Success, 3, false },
    { 1.3300000000000007, 3.2978413624934477, gsStatus::Success, 3, false },
    { 1.4000000000000008, 3.3197788136095094, gsStatus::Success, 3, false },
};
static const T1Meta T1_LoadControl_A0_Q0_meta = { 20, false, true };
// T1_LoadControl_A0_Q1: ds=0.07 N=20 adaptive=0 quasi=1 boundary 3 / 15 steps, limit point passed 0, final status gsStatus::NotConverged at lambda 0.98
static const std::vector<T1Step> T1_LoadControl_A0_Q1 = {
    { 0.070000000000000007, 0.030975053280708618, gsStatus::Success, 4, false },
    { 0.14000000000000001, 0.062767990466368925, gsStatus::Success, 4, false },
    { 0.21000000000000002, 0.09548564722213905, gsStatus::Success, 4, false },
    { 0.28000000000000003, 0.12925955248539134, gsStatus::Success, 4, false },
    { 0.35000000000000003, 0.1642546559128632, gsStatus::Success, 4, false },
    { 0.42000000000000004, 0.20068245427559597, gsStatus::Success, 4, true },
    { 0.49000000000000005, 0.23882156737663299, gsStatus::Success, 5, false },
    { 0.56000000000000005, 0.27905167803791864, gsStatus::Success, 5, false },
    { 0.63000000000000012, 0.32191324971498331, gsStatus::Success, 5, false },
    { 0.70000000000000018, 0.36822188968364389, gsStatus::Success, 5, false },
    { 0.77000000000000024, 0.41931450588763775, gsStatus::Success, 6, false },
    { 0.8400000000000003, 0.47767975029241211, gsStatus::Success, 6, true },
    { 0.91000000000000036, 0.54913022011857215, gsStatus::Success, 8, false },
    { 0.98000000000000043, 0.65851078197091495, gsStatus::Success, 16, true },
    { 0.98000000000000043, 0.65851078197091495, gsStatus::NotConverged, 99, false },
};
static const T1Meta T1_LoadControl_A0_Q1_meta = { 15, false, true };
// T1_LoadControl_A1_Q0: ds=0.07 N=20 adaptive=1 quasi=0 boundary 0 / 8 steps, limit point passed 0, final status gsStatus::NotConverged at lambda 0.91
// T1_LoadControl_A1_Q0: fragile from step 8
// T1_LoadControl_A1_Q0: load-control prefix ends at step 8
static const std::vector<T1Step> T1_LoadControl_A1_Q0 = {
    { 0.070000000000000007, 0.030975053280708726, gsStatus::Success, 3, false },
    { 0.16333333333333333, 0.073565211902333671, gsStatus::Success, 3, false },
    { 0.2877777777777778, 0.13308418379077394, gsStatus::Success, 4, false },
    { 0.41222222222222227, 0.19655622486626986, gsStatus::Success, 4, false },
    { 0.53666666666666674, 0.26538017521821894, gsStatus::Success, 4, false },
    { 0.6611111111111112, 0.34200611859560642, gsStatus::Success, 4, false },
    { 0.78555555555555567, 0.43153037869404243, gsStatus::Success, 4, false },
    { 0.91000000000000014, 0.54913022011864188, gsStatus::Success, 5, false },
    { 0.91000000000000014, 0.54913022011864188, gsStatus::NotConverged, 99, false },
};
static const T1Meta T1_LoadControl_A1_Q0_meta = { 8, true, true };
// T1_LoadControl_A1_Q1: ds=0.07 N=20 adaptive=1 quasi=1 boundary 0 / 5 steps, limit point passed 0, final status gsStatus::Success at lambda 0.882294
// T1_LoadControl_A1_Q1: fragile from step 5
// T1_LoadControl_A1_Q1: adaptive prefix ends at step 5
static const std::vector<T1Step> T1_LoadControl_A1_Q1 = {
    { 0.070000000000000007, 0.030975053280708618, gsStatus::Success, 4, false },
    { 0.14000000000000001, 0.062767990466368925, gsStatus::Success, 4, false },
    { 0.21000000000000002, 0.09548564722213905, gsStatus::Success, 4, false },
    { 0.28000000000000003, 0.12925955248539134, gsStatus::Success, 4, false },
    { 0.35000000000000003, 0.1642546559128632, gsStatus::Success, 4, false },
    { 0.42000000000000004, 0.20068245427559597, gsStatus::Success, 4, true },
    { 0.49000000000000005, 0.23882156737663299, gsStatus::Success, 5, true },
    { 0.54600000000000004, 0.27081543603609459, gsStatus::Success, 5, true },
    { 0.59079999999999999, 0.29754653147492099, gsStatus::Success, 4, true },
    { 0.63559999999999994, 0.32547763068878532, gsStatus::Success, 5, true },
    { 0.67143999999999993, 0.34884223242158524, gsStatus::Success, 4, true },
    { 0.70727999999999991, 0.37328266668093696, gsStatus::Success, 5, false },
    { 0.73595199999999994, 0.39374810032443647, gsStatus::Success, 4, false },
    { 0.76462399999999997, 0.41517560576307438, gsStatus::Success, 4, true },
    { 0.793296, 0.43775090051494531, gsStatus::Success, 5, true },
    { 0.8162336, 0.45680464211514132, gsStatus::Success, 5, true },
    { 0.83458368000000005, 0.47280853840321269, gsStatus::Success, 4, true },
    { 0.8529337600000001, 0.48961616201650782, gsStatus::Success, 5, false },
    { 0.86761382400000009, 0.50374546117734098, gsStatus::Success, 4, false },
    { 0.88229388800000008, 0.51859609255617123, gsStatus::Success, 5, true },
};
static const T1Meta T1_LoadControl_A1_Q1_meta = { 5, true, true };
// T1_Riks_A0_Q0: ds=0.05 N=40 adaptive=0 quasi=0 boundary 1 / 40 steps, limit point passed 1, final status gsStatus::Success at lambda -0.243381
static const std::vector<T1Step> T1_Riks_A0_Q0 = {
    { 0.064672981509771882, 0.028590303647158165, gsStatus::Success, 6, false },
    { 0.12908899508374924, 0.057754503561862354, gsStatus::Success, 4, false },
    { 0.19321512499062035, 0.087549640124625863, gsStatus::Success, 4, false },
    { 0.25701168331909807, 0.11804223466296489, gsStatus::Success, 4, false },
    { 0.3204302368885546, 0.14931055822798403, gsStatus::Success, 4, false },
    { 0.38341087083102743, 0.18144762188101965, gsStatus::Success, 4, false },
    { 0.44587831169035352, 0.21456517748638257, gsStatus::Success, 4, false },
    { 0.50773630103733924, 0.24879915718764339, gsStatus::Success, 4, false },
    { 0.56885920485783104, 0.2843171949472566, gsStatus::Success, 4, false },
    { 0.62907910786752486, 0.32132920631229156, gsStatus::Success, 4, false },
    { 0.68816525327252709, 0.36010250586085019, gsStatus::Success, 4, false },
    { 0.74578995903994483, 0.40098364181154189, gsStatus::Success, 4, false },
    { 0.80146955842746737, 0.44442982495881855, gsStatus::Success, 4, false },
    { 0.85445714904830727, 0.491052215940533, gsStatus::Success, 4, false },
    { 0.90353937796696715, 0.54166521993954331, gsStatus::Success, 4, false },
    { 0.94664589344980266, 0.59729390406887906, gsStatus::Success, 5, false },
    { 0.98017091996545092, 0.65892003688561296, gsStatus::Success, 5, false },
    { 0.99840668124471099, 0.72630142733759673, gsStatus::Success, 5, false },
    { 0.99560949347742145, 0.7956561336616812, gsStatus::Success, 5, false },
    { 0.97189679062091283, 0.86067635631308514, gsStatus::Success, 5, false },
    { 0.93325122792738979, 0.91812945186232831, gsStatus::Success, 5, false },
    { 0.88551061017800814, 0.96844742851747057, gsStatus::Success, 5, false },
    { 0.8322826503422629, 1.0131130888055768, gsStatus::Success, 5, false },
    { 0.7755886966621619, 1.0534693002277824, gsStatus::Success, 4, false },
    { 0.71658690816318471, 1.0905304699900669, gsStatus::Success, 4, false },
    { 0.65597671695986526, 1.1250396893416463, gsStatus::Success, 4, false },
    { 0.59420428180324247, 1.1575469992969767, gsStatus::Success, 4, false },
    { 0.53156872531403743, 1.1884689245712463, gsStatus::Success, 4, false },
    { 0.46827963609979173, 1.218129235136701, gsStatus::Success, 4, false },
    { 0.40448979062263712, 1.2467863941676258, gsStatus::Success, 4, false },
    { 0.34031468781681429, 1.2746522518931631, gsStatus::Success, 4, false },
    { 0.27584474495788341, 1.301905065469829, gsStatus::Success, 4, false },
    { 0.21115324482323289, 1.3286988260299373, gsStatus::Success, 4, false },
    { 0.14630174376054733, 1.3551701647470402, gsStatus::Success, 4, false },
    { 0.081343931259869315, 1.3814436697688626, gsStatus::Success, 3, false },
    { 0.0163285445278325, 1.407636176423523, gsStatus::Success, 3, false },
    { -0.048698270585885522, 1.4338604304184781, gsStatus::Success, 3, false },
    { -0.11369087508343914, 1.4602284298287127, gsStatus::Success, 3, false },
    { -0.17860210528568929, 1.4868547052728232, gsStatus::Success, 3, false },
    { -0.24338108086515753, 1.5138597885966416, gsStatus::Success, 3, true },
};
static const T1Meta T1_Riks_A0_Q0_meta = { 40, false, true };
// T1_Riks_A0_Q1: ds=0.05 N=40 adaptive=0 quasi=1 boundary 5 / 40 steps, limit point passed 1, final status gsStatus::Success at lambda -0.243381
static const std::vector<T1Step> T1_Riks_A0_Q1 = {
    { 0.064672981509758309, 0.028590303647188866, gsStatus::Success, 8, false },
    { 0.12908899508374916, 0.057754503561863353, gsStatus::Success, 4, false },
    { 0.19321512499061996, 0.087549640124627584, gsStatus::Success, 4, false },
    { 0.25701168331909696, 0.11804223466296802, gsStatus::Success, 4, false },
    { 0.32043023688855488, 0.14931055822798425, gsStatus::Success, 5, true },
    { 0.38341087083102776, 0.18144762188101982, gsStatus::Success, 5, false },
    { 0.44587831169035375, 0.21456517748638285, gsStatus::Success, 5, false },
    { 0.50773630103733947, 0.2487991571876437, gsStatus::Success, 5, false },
    { 0.56885920485783104, 0.28431719494725738, gsStatus::Success, 5, false },
    { 0.62907910786752408, 0.32132920631229356, gsStatus::Success, 5, false },
    { 0.68816525327252309, 0.3601025058608569, gsStatus::Success, 5, false },
    { 0.74578995903994538, 0.4009836418115425, gsStatus::Success, 6, true },
    { 0.80146955842746714, 0.44442982495882005, gsStatus::Success, 6, false },
    { 0.85445714904830128, 0.49105221594054077, gsStatus::Success, 6, false },
    { 0.9035393779669042, 0.54166521993960415, gsStatus::Success, 6, true },
    { 0.94664589344979744, 0.5972939040688976, gsStatus::Success, 7, false },
    { 0.98017091996535333, 0.65892003688567868, gsStatus::Success, 7, false },
    { 0.99840668124459409, 0.72630142733766756, gsStatus::Success, 7, false },
    { 0.99560949347735372, 0.79565613366175569, gsStatus::Success, 6, false },
    { 0.97189679062086909, 0.86067635631316919, gsStatus::Success, 6, false },
    { 0.93325122792733128, 0.91812945186240269, gsStatus::Success, 7, false },
    { 0.88551061017793953, 0.96844742851753551, gsStatus::Success, 7, false },
    { 0.83228265034219606, 1.0131130888056441, gsStatus::Success, 6, false },
    { 0.77558869666208463, 1.0534693002278348, gsStatus::Success, 6, false },
    { 0.71658690816310444, 1.0905304699901146, gsStatus::Success, 6, true },
    { 0.65597671695978432, 1.1250396893416934, gsStatus::Success, 5, false },
    { 0.59420428180315932, 1.1575469992970195, gsStatus::Success, 5, false },
    { 0.53156872531395305, 1.1884689245712867, gsStatus::Success, 5, false },
    { 0.46827963609970658, 1.2181292351367399, gsStatus::Success, 5, false },
    { 0.40448979062255197, 1.2467863941676651, gsStatus::Success, 4, true },
    { 0.34031468781672825, 1.2746522518932002, gsStatus::Success, 4, false },
    { 0.27584474495779693, 1.3019050654698652, gsStatus::Success, 4, false },
    { 0.2111532448231461, 1.3286988260299728, gsStatus::Success, 4, false },
    { 0.14630174376046037, 1.3551701647470755, gsStatus::Success, 4, false },
    { 0.081343931259782246, 1.3814436697688974, gsStatus::Success, 4, false },
    { 0.016328544527745417, 1.4076361764235579, gsStatus::Success, 3, false },
    { -0.048698270585972592, 1.433860430418513, gsStatus::Success, 3, false },
    { -0.1136908750835261, 1.4602284298287482, gsStatus::Success, 3, false },
    { -0.17860210528577614, 1.4868547052728591, gsStatus::Success, 4, false },
    { -0.24338108086524426, 1.5138597885966782, gsStatus::Success, 4, false },
};
static const T1Meta T1_Riks_A0_Q1_meta = { 40, false, true };
// T1_Riks_A1_Q0: ds=0.05 N=40 adaptive=1 quasi=0 boundary 0 / 40 steps, limit point passed 1, final status gsStatus::Success at lambda 0.907123
static const std::vector<T1Step> T1_Riks_A1_Q0 = {
    { 0.064672981509771882, 0.028590303647158165, gsStatus::Success, 6, false },
    { 0.10764735613317858, 0.047965949697939575, gsStatus::Success, 4, false },
    { 0.15049818327825329, 0.06761294755946623, gsStatus::Success, 4, false },
    { 0.19321480086207604, 0.087549487446565363, gsStatus::Success, 4, false },
    { 0.23578510263380814, 0.10779575772555394, gsStatus::Success, 4, false },
    { 0.27819526782660103, 0.12837425416806983, gsStatus::Success, 4, false },
    { 0.32042942533669799, 0.14931015141951207, gsStatus::Success, 4, false },
    { 0.36246923273788328, 0.17063175211343481, gsStatus::Success, 4, false },
    { 0.40429334328776217, 0.19237103358833577, gsStatus::Success, 4, false },
    { 0.44587672385931632, 0.21456431819155208, gsStatus::Success, 4, false },
    { 0.48718977190320001, 0.23725310121802137, gsStatus::Success, 4, false },
    { 0.52819715770289732, 0.26048508133349596, gsStatus::Success, 4, false },
    { 0.56885628547363554, 0.28431545276487802, gsStatus::Success, 4, false },
    { 0.60911521701413551, 0.30880853763007499, gsStatus::Success, 4, false },
    { 0.64890982429601807, 0.33403986141420577, gsStatus::Success, 4, false },
    { 0.68815981519282077, 0.3600988046717396, gsStatus::Success, 4, false },
    { 0.72676307998463108, 0.38709199604093902, gsStatus::Success, 4, false },
    { 0.76458748489439887, 0.41514763226436024, gsStatus::Success, 4, false },
    { 0.80145870740915326, 0.44442088003609043, gsStatus::Success, 4, false },
    { 0.83714182850868502, 0.47510032105267441, gsStatus::Success, 4, false },
    { 0.87131297603152502, 0.50741474857440816, gsStatus::Success, 4, false },
    { 0.90351521541020863, 0.54163773777655211, gsStatus::Success, 4, false },
    { 0.93309067396767598, 0.57808237658356409, gsStatus::Success, 4, false },
    { 0.95908249669459567, 0.61706621949011375, gsStatus::Success, 4, false },
    { 0.98012156301772624, 0.65880168110643356, gsStatus::Success, 5, false },
    { 0.99217506171754932, 0.69407951275447666, gsStatus::Success, 5, false },
    { 0.9980154585662705, 0.72324173946309833, gsStatus::Success, 5, false },
    { 0.99992587717613368, 0.74690370377248616, gsStatus::Success, 5, false },
    { 0.99958487839124877, 0.76586001132266046, gsStatus::Success, 4, false },
    { 0.99756713917365203, 0.78468763492622495, gsStatus::Success, 4, false },
    { 0.99391378920888318, 0.80324837827695106, gsStatus::Success, 4, false },
    { 0.98871205422846775, 0.82142288626876159, gsStatus::Success, 4, false },
    { 0.98208274816907948, 0.83911932798855127, gsStatus::Success, 4, false },
    { 0.97416541245230071, 0.85627623720334267, gsStatus::Success, 4, false },
    { 0.96510482512104079, 0.8728604066356731, gsStatus::Success, 4, false },
    { 0.95504097840509683, 0.88886186759244568, gsStatus::Success, 4, false },
    { 0.94410292619656588, 0.90428794215554764, gsStatus::Success, 4, false },
    { 0.93240586022915128, 0.91915769706911188, gsStatus::Success, 4, false },
    { 0.92005042898952849, 0.9334974198029895, gsStatus::Success, 4, false },
    { 0.90712339410306686, 0.94733724337495584, gsStatus::Success, 4, false },
};
static const T1Meta T1_Riks_A1_Q0_meta = { 40, false, true };
// T1_Riks_A1_Q1: ds=0.05 N=40 adaptive=1 quasi=1 boundary 0 / 18 steps, limit point passed 0, final status gsStatus::Success at lambda 0.933768
// T1_Riks_A1_Q1: fragile from step 18
// T1_Riks_A1_Q1: adaptive prefix ends at step 18
static const std::vector<T1Step> T1_Riks_A1_Q1 = {
    { 0.064672981509758309, 0.028590303647188866, gsStatus::Success, 8, false },
    { 0.09691493249251773, 0.043097258400547414, gsStatus::Success, 4, false },
    { 0.12908879710331092, 0.057754412792756255, gsStatus::Success, 4, false },
    { 0.16119028523643369, 0.072569148094405525, gsStatus::Success, 4, false },
    { 0.19321468731540267, 0.087549433961352494, gsStatus::Success, 4, false },
    { 0.2251568182762167, 0.10270389385506414, gsStatus::Success, 4, false },
    { 0.25701095199163654, 0.11804187982772105, gsStatus::Success, 4, false },
    { 0.2887707441364159, 0.13357355831122469, gsStatus::Success, 4, false },
    { 0.32042914099464009, 0.14931000888675491, gsStatus::Success, 4, false },
    { 0.35197827106641244, 0.16526333842790186, gsStatus::Success, 4, false },
    { 0.38340931549357904, 0.18144681352578937, gsStatus::Success, 4, false },
    { 0.41471235222617181, 0.19787501474435076, gsStatus::Success, 4, false },
    { 0.44587616739989105, 0.21456401704971595, gsStatus::Success, 4, false },
    { 0.47688802545989595, 0.23153160174872625, gsStatus::Success, 4, false },
    { 0.50773338696286019, 0.24879750650537999, gsStatus::Success, 4, false },
    { 0.53839555945379924, 0.26638372153593198, gsStatus::Success, 4, false },
    { 0.56885526196572678, 0.28431484197366858, gsStatus::Success, 4, false },
    { 0.59909007697384475, 0.30261848869903685, gsStatus::Success, 4, false },
    { 0.62907375423598422, 0.32132581267956384, gsStatus::Success, 4, true },
    { 0.65877531765854691, 0.34047210101182507, gsStatus::Success, 5, true },
    { 0.68230880127004645, 0.35613190949844858, gsStatus::Success, 4, true },
    { 0.70561606846084879, 0.37212138380057458, gsStatus::Success, 4, false },
    { 0.72867126244268732, 0.38846636220173009, gsStatus::Success, 5, true },
    { 0.74691363638171404, 0.40181755417042925, gsStatus::Success, 4, true },
    { 0.76495738018767567, 0.41543108744968338, gsStatus::Success, 4, false },
    { 0.78278184424234654, 0.42932495311121233, gsStatus::Success, 4, true },
    { 0.80036321343088601, 0.44351890955056711, gsStatus::Success, 5, true },
    { 0.81423479232208251, 0.4551045858281918, gsStatus::Success, 4, true },
    { 0.82791700088156039, 0.46690823986078361, gsStatus::Success, 4, false },
    { 0.84139146749336247, 0.47894287579033101, gsStatus::Success, 4, false },
    { 0.854637284734036, 0.49122246554958648, gsStatus::Success, 5, true },
    { 0.86505330930187552, 0.50123245937137495, gsStatus::Success, 4, true },
    { 0.87529397914200757, 0.51141689330822315, gsStatus::Success, 4, false },
    { 0.88534387886823163, 0.52178428170498048, gsStatus::Success, 4, false },
    { 0.89518575757183783, 0.53234350280191767, gsStatus::Success, 4, false },
    { 0.90480028063790718, 0.54310374504142922, gsStatus::Success, 5, false },
    { 0.91231364725971575, 0.55186291976205215, gsStatus::Success, 4, false },
    { 0.91965534567209628, 0.56076162951493225, gsStatus::Success, 4, false },
    { 0.92681185477814076, 0.5698047252359143, gsStatus::Success, 4, false },
    { 0.93376834234243389, 0.57899696201150797, gsStatus::Success, 4, false },
};
static const T1Meta T1_Riks_A1_Q1_meta = { 18, true, true };
// T1_Crisfield_A0_Q0: ds=0.05 N=40 adaptive=0 quasi=0 boundary 0 / 40 steps, limit point passed 1, final status gsStatus::Success at lambda 0.441009
static const std::vector<T1Step> T1_Crisfield_A0_Q0 = {
    { 0.080311815895432251, 0.035604821679769909, gsStatus::Success, 3, false },
    { 0.15885630407663703, 0.071485363644873276, gsStatus::Success, 3, false },
    { 0.23552638360449152, 0.1076715075278195, gsStatus::Success, 3, false },
    { 0.31020179299310968, 0.14419743057189255, gsStatus::Success, 3, false },
    { 0.38274649344076228, 0.18110239471714906, gsStatus::Success, 3, false },
    { 0.45300536501358024, 0.21843170193109937, gsStatus::Success, 3, false },
    { 0.52079995302724558, 0.2562378465175259, gsStatus::Success, 3, false },
    { 0.58592292401472756, 0.29458189136671442, gsStatus::Success, 3, false },
    { 0.64813075047431712, 0.33353507718062608, gsStatus::Success, 3, false },
    { 0.70713394546694941, 0.37318062298750454, gsStatus::Success, 3, false },
    { 0.7625838988292557, 0.41361555192819421, gsStatus::Success, 3, false },
    { 0.81405503616968122, 0.45495208908424095, gsStatus::Success, 3, false },
    { 0.86102073082100039, 0.49731753417352198, gsStatus::Success, 3, false },
    { 0.90282154352245414, 0.54085010046100912, gsStatus::Success, 3, false },
    { 0.938626208297897, 0.58568522141409829, gsStatus::Success, 3, false },
    { 0.96739288872502482, 0.63192096535846864, gsStatus::Success, 3, false },
    { 0.98785909400086613, 0.67954187824440349, gsStatus::Success, 3, false },
    { 0.9986364359039096, 0.72827490715326015, gsStatus::Success, 3, false },
    { 0.99855028458108663, 0.77738397178726915, gsStatus::Success, 3, false },
    { 0.98730630089402227, 0.82554538970087832, gsStatus::Success, 3, false },
    { 0.96613320862725316, 0.871098925166087, gsStatus::Success, 3, false },
    { 0.93762625338037831, 0.91269846180486125, gsStatus::Success, 3, false },
    { 0.90472978192481968, 0.94978992459062894, gsStatus::Success, 3, false },
    { 0.86982820347362355, 0.98252580306267945, gsStatus::Success, 3, false },
    { 0.83449815893791701, 1.0114093857415567, gsStatus::Success, 3, false },
    { 0.79966625622938181, 1.0370202560095503, gsStatus::Success, 3, false },
    { 0.76583666023938624, 1.0598882169992045, gsStatus::Success, 3, false },
    { 0.73326257201798906, 1.0804574432881104, gsStatus::Success, 3, false },
    { 0.70205343502380746, 1.0990882062431608, gsStatus::Success, 3, false },
    { 0.67223763915717505, 1.1160699527634612, gsStatus::Success, 3, false },
    { 0.64379848004981233, 1.1316355680581822, gsStatus::Success, 3, false },
    { 0.61669476938072576, 1.1459736664209832, gsStatus::Success, 3, false },
    { 0.59087271717699741, 1.1592383617164304, gsStatus::Success, 3, false },
    { 0.56627282041970495, 1.1715567951070776, gsStatus::Success, 3, false },
    { 0.5428338600274718, 1.1830348748656072, gsStatus::Success, 3, false },
    { 0.52049520123110327, 1.1937616511288913, gsStatus::Success, 3, false },
    { 0.49919808670388566, 1.2038126656439416, gsStatus::Success, 3, false },
    { 0.47888632648774776, 1.2132525350382861, gsStatus::Success, 3, false },
    { 0.45950662517980934, 1.222136959671152, gsStatus::Success, 3, false },
    { 0.44100869146463567, 1.2305142995989196, gsStatus::Success, 3, false },
};
static const T1Meta T1_Crisfield_A0_Q0_meta = { 40, false, true };
// T1_Crisfield_A0_Q1: ds=0.05 N=40 adaptive=0 quasi=1 boundary 6 / 40 steps, limit point passed 1, final status gsStatus::Success at lambda 0.441009
static const std::vector<T1Step> T1_Crisfield_A0_Q1 = {
    { 0.080311815895432279, 0.035604821679769902, gsStatus::Success, 4, false },
    { 0.15885630407663712, 0.071485363644873262, gsStatus::Success, 4, false },
    { 0.23552638360449163, 0.1076715075278195, gsStatus::Success, 4, false },
    { 0.31020179299310985, 0.14419743057189258, gsStatus::Success, 4, false },
    { 0.3827464934407625, 0.18110239471714906, gsStatus::Success, 4, false },
    { 0.45300536501358069, 0.21843170193109934, gsStatus::Success, 4, false },
    { 0.52079995302724624, 0.25623784651752585, gsStatus::Success, 4, false },
    { 0.58592292401472856, 0.2945818913667142, gsStatus::Success, 4, false },
    { 0.6481307504743189, 0.33353507718062564, gsStatus::Success, 4, false },
    { 0.70713394546695263, 0.37318062298750365, gsStatus::Success, 4, false },
    { 0.76258389882926159, 0.41361555192819249, gsStatus::Success, 4, false },
    { 0.81405503616969166, 0.45495208908423784, gsStatus::Success, 4, false },
    { 0.8610207308210196, 0.49731753417351632, gsStatus::Success, 4, true },
    { 0.90282154352245458, 0.54085010046100934, gsStatus::Success, 5, true },
    { 0.93862620829789767, 0.58568522141409851, gsStatus::Success, 5, false },
    { 0.96739288872502538, 0.63192096535846898, gsStatus::Success, 5, false },
    { 0.98785909400086647, 0.67954187824440393, gsStatus::Success, 5, true },
    { 0.99863643590392603, 0.72827490715325949, gsStatus::Success, 4, false },
    { 0.99855028458108719, 0.77738397178726781, gsStatus::Success, 4, false },
    { 0.98730630089402283, 0.82554538970087699, gsStatus::Success, 3, false },
    { 0.96613320862725394, 0.87109892516608578, gsStatus::Success, 4, false },
    { 0.93762625338038597, 0.9126984618048638, gsStatus::Success, 4, false },
    { 0.90472978192482179, 0.94978992459062705, gsStatus::Success, 5, false },
    { 0.86982820347362555, 0.98252580306267767, gsStatus::Success, 5, false },
    { 0.83449815893791912, 1.0114093857415551, gsStatus::Success, 5, true },
    { 0.79966625622938459, 1.0370202560095505, gsStatus::Success, 4, false },
    { 0.76583666023938846, 1.0598882169992037, gsStatus::Success, 4, false },
    { 0.73326257201799105, 1.0804574432881093, gsStatus::Success, 4, false },
    { 0.70205343502380935, 1.0990882062431595, gsStatus::Success, 4, false },
    { 0.67223763915717694, 1.1160699527634603, gsStatus::Success, 4, false },
    { 0.6437984800498141, 1.1316355680581813, gsStatus::Success, 4, false },
    { 0.61669476938072743, 1.1459736664209821, gsStatus::Success, 4, false },
    { 0.59087271717699896, 1.1592383617164297, gsStatus::Success, 4, false },
    { 0.56627282041970639, 1.1715567951070769, gsStatus::Success, 4, true },
    { 0.54283386002747325, 1.1830348748656068, gsStatus::Success, 3, true },
    { 0.52049520123110471, 1.1937616511288909, gsStatus::Success, 3, false },
    { 0.49919808670388699, 1.203812665643941, gsStatus::Success, 3, false },
    { 0.47888632648774904, 1.2132525350382857, gsStatus::Success, 3, false },
    { 0.45950662517981056, 1.2221369596711515, gsStatus::Success, 3, false },
    { 0.44100869146463684, 1.2305142995989191, gsStatus::Success, 3, false },
};
static const T1Meta T1_Crisfield_A0_Q1_meta = { 40, false, true };
// T1_Crisfield_A1_Q0: ds=0.05 N=40 adaptive=1 quasi=0 boundary 0 / 31 steps, limit point passed 1, final status gsStatus::Success at lambda -264.334
// T1_Crisfield_A1_Q0: fragile from step 31
// T1_Crisfield_A1_Q0: adaptive prefix ends at step 31
static const std::vector<T1Step> T1_Crisfield_A1_Q0 = {
    { 0.080311815895432251, 0.035604821679769909, gsStatus::Success, 3, false },
    { 0.18474910087809604, 0.083569205308693828, gsStatus::Success, 3, false },
    { 0.31886791692414351, 0.14852774166111879, gsStatus::Success, 3, false },
    { 0.48730912976424207, 0.23731965041645334, gsStatus::Success, 4, false },
    { 0.64116943795319059, 0.32904469880755671, gsStatus::Success, 4, false },
    { 0.77691910185353497, 0.4247028662384113, gsStatus::Success, 4, false },
    { 0.88903586342123198, 0.52569531312308448, gsStatus::Success, 4, false },
    { 0.96832024773775971, 0.63368544544676975, gsStatus::Success, 4, false },
    { 0.99996724134144577, 0.74879613481271567, gsStatus::Success, 4, false },
    { 0.97080383637901646, 0.86273135458159023, gsStatus::Success, 3, false },
    { 0.86578021655124815, 0.98601911934727271, gsStatus::Success, 4, false },
    { 0.75035143411542149, 1.0698244945573019, gsStatus::Success, 4, false },
    { 0.6485466886337734, 1.1290757983946473, gsStatus::Success, 4, false },
    { 0.56193384508794808, 1.1736998582526545, gsStatus::Success, 3, false },
    { 0.46386479447540652, 1.2201489125713758, gsStatus::Success, 3, false },
    { 0.35874119537171362, 1.2667357752937056, gsStatus::Success, 3, false },
    { 0.25351309053421817, 1.3112080021751515, gsStatus::Success, 3, false },
    { 0.15729365776666926, 1.3507042766894095, gsStatus::Success, 3, false },
    { 0.079843053315429457, 1.3820489864745233, gsStatus::Success, 3, false },
    { 0.028564175004507864, 1.402708715853797, gsStatus::Success, 3, false },
    { 0.0044576909371206556, 1.4124176490098268, gsStatus::Success, 2, false },
    { -0.0030149226815992473, 1.4154283703884956, gsStatus::Success, 2, false },
    { -0.013101508815731818, 1.4194936779270844, gsStatus::Success, 2, false },
    { -0.10051340834977256, 1.4548630935963454, gsStatus::Success, 3, false },
    { -0.97214117495173946, 1.9840422730114133, gsStatus::Success, 6, false },
    { -4.748138201398648, 1.5228844668599848, gsStatus::Success, 10, false },
    { -15.067940673095523, 4.1140117388354458, gsStatus::Success, 5, false },
    { -24.455645618743105, 6.2895081798528478, gsStatus::Success, 4, false },
    { -34.306839484454443, 8.5000537215280776, gsStatus::Success, 4, false },
    { -44.454163506371138, 10.733384054961437, gsStatus::Success, 4, false },
    { -54.809441658377779, 12.982959074496211, gsStatus::Success, 4, false },
    { -65.320066690705858, 15.244823606189104, gsStatus::Success, 4, true },
    { -75.951901278470018, 17.516378548204415, gsStatus::Success, 3, true },
    { -90.263384908904158, 20.554350153320947, gsStatus::Success, 3, true },
    { -109.53707727792114, 24.618786295325467, gsStatus::Success, 4, true },
    { -128.99997045486987, 28.699749149568131, gsStatus::Success, 3, true },
    { -155.15508328413296, 34.156508003842355, gsStatus::Success, 4, true },
    { -181.51325382443437, 39.631562310191633, gsStatus::Success, 3, true },
    { -216.87782454660666, 46.94891286046137, gsStatus::Success, 3, true },
    { -264.33442827909852, 56.730076291075655, gsStatus::Success, 4, true },
};
static const T1Meta T1_Crisfield_A1_Q0_meta = { 31, true, true };
// T1_Crisfield_A1_Q1: ds=0.05 N=40 adaptive=1 quasi=1 boundary 0 / 12 steps, limit point passed 1, final status gsStatus::Success at lambda 0.485453
// T1_Crisfield_A1_Q1: fragile from step 12
// T1_Crisfield_A1_Q1: adaptive prefix ends at step 12
static const std::vector<T1Step> T1_Crisfield_A1_Q1 = {
    { 0.080311815895432279, 0.035604821679769902, gsStatus::Success, 4, false },
    { 0.15885630407663712, 0.071485363644873262, gsStatus::Success, 4, false },
    { 0.23552638360449163, 0.1076715075278195, gsStatus::Success, 4, false },
    { 0.31020179299310985, 0.14419743057189258, gsStatus::Success, 4, false },
    { 0.3827464934407625, 0.18110239471714906, gsStatus::Success, 4, false },
    { 0.45300536501358069, 0.21843170193109934, gsStatus::Success, 4, false },
    { 0.52079995302724624, 0.25623784651752585, gsStatus::Success, 4, false },
    { 0.58592292401472856, 0.2945818913667142, gsStatus::Success, 4, false },
    { 0.6481307504743189, 0.33353507718062564, gsStatus::Success, 4, false },
    { 0.70713394546695263, 0.37318062298750365, gsStatus::Success, 4, false },
    { 0.76258389882926159, 0.41361555192819249, gsStatus::Success, 4, false },
    { 0.81405503616969166, 0.45495208908423784, gsStatus::Success, 4, false },
    { 0.8610207308210196, 0.49731753417351632, gsStatus::Success, 4, true },
    { 0.90282154352245458, 0.54085010046100934, gsStatus::Success, 5, true },
    { 0.93197057359325197, 0.57658007526026023, gsStatus::Success, 4, true },
    { 0.95677397542338127, 0.61320119226448022, gsStatus::Success, 4, false },
    { 0.97661738916098262, 0.65073277284088926, gsStatus::Success, 4, false },
    { 0.99080578368158312, 0.68911781742778977, gsStatus::Success, 4, false },
    { 0.99862373998210352, 0.72816173752969249, gsStatus::Success, 4, false },
    { 0.99947894929957692, 0.76745871831646262, gsStatus::Success, 4, false },
    { 0.99314070462455195, 0.80634891951283916, gsStatus::Success, 3, true },
    { 0.97421310264262928, 0.85618166360692982, gsStatus::Success, 4, true },
    { 0.9458276405132453, 0.90196531278252035, gsStatus::Success, 4, false },
    { 0.91155296979715739, 0.94271216629025356, gsStatus::Success, 5, false },
    { 0.88215161922476903, 0.97154083116343726, gsStatus::Success, 4, true },
    { 0.85212579549562517, 0.99742943883004398, gsStatus::Success, 4, false },
    { 0.82218900577285237, 1.0207379899268521, gsStatus::Success, 4, false },
    { 0.792777249326022, 1.0418199206970513, gsStatus::Success, 4, false },
    { 0.76414386411236235, 1.0609893301474869, gsStatus::Success, 4, false },
    { 0.73642691824446449, 1.0785132138145348, gsStatus::Success, 4, false },
    { 0.70969309065333841, 1.094614276725616, gsStatus::Success, 4, false },
    { 0.68396525855535584, 1.1094771718525385, gsStatus::Success, 4, false },
    { 0.65923965105126114, 1.1232551181152854, gsStatus::Success, 4, false },
    { 0.63549651557967879, 1.1360757968076156, gsStatus::Success, 4, false },
    { 0.61270678029026171, 1.1480462536589184, gsStatus::Success, 4, false },
    { 0.5908362368288621, 1.1592568503922844, gsStatus::Success, 4, true },
    { 0.56984817513474983, 1.1697844082947868, gsStatus::Success, 3, true },
    { 0.54298599784887114, 1.1829611325843283, gsStatus::Success, 4, true },
    { 0.517559884593771, 1.1951567787346289, gsStatus::Success, 3, true },
    { 0.48545251456508859, 1.2102157665444804, gsStatus::Success, 4, true },
};
static const T1Meta T1_Crisfield_A1_Q1_meta = { 12, true, true };
// T1_ConsistentCrisfield_A0_Q0: ds=0.05 N=40 adaptive=0 quasi=0 boundary 2 / 40 steps, limit point passed 1, final status gsStatus::Success at lambda 0.0698236
static const std::vector<T1Step> T1_ConsistentCrisfield_A0_Q0 = {
    { 0.080311815895432251, 0.035604821679769909, gsStatus::Success, 5, false },
    { 0.158302485065131, 0.071228353776569928, gsStatus::Success, 4, false },
    { 0.23386076037132991, 0.10687195182666948, gsStatus::Success, 4, false },
    { 0.30686854986347417, 0.1425371997069137, gsStatus::Success, 4, false },
    { 0.37720051614459493, 0.17822597878665805, gsStatus::Success, 4, false },
    { 0.44472367270261798, 0.21394056352045102, gsStatus::Success, 4, false },
    { 0.50929698271388335, 0.24968375532455817, gsStatus::Success, 4, false },
    { 0.57077096643069014, 0.28545907287004385, gsStatus::Success, 4, false },
    { 0.62898732305057947, 0.3212710272520079, gsStatus::Success, 4, false },
    { 0.68377857096505534, 0.35712552802135317, gsStatus::Success, 4, false },
    { 0.73496770461818739, 0.39303049698238252, gsStatus::Success, 4, false },
    { 0.78236785311897294, 0.42899682372351194, gsStatus::Success, 4, false },
    { 0.82578189744075781, 0.46503990798331102, gsStatus::Success, 4, false },
    { 0.86500194226911697, 0.50118226482669559, gsStatus::Success, 4, false },
    { 0.8998084040746045, 0.537458187672056, gsStatus::Success, 4, false },
    { 0.92996815909180608, 0.57392275299833051, gsStatus::Success, 4, true },
    { 0.95523035808614642, 0.61067109122444729, gsStatus::Success, 5, false },
    { 0.97531591438531251, 0.64788611413728503, gsStatus::Success, 5, false },
    { 0.98988614377532624, 0.68598644070408066, gsStatus::Success, 5, false },
    { 0.99840717110596766, 0.72630547902605846, gsStatus::Success, 5, false },
    { 0.99878145081398373, 0.77532468254569242, gsStatus::Success, 10, false },
    { 0.99433266784917906, 0.80148465992039497, gsStatus::Success, 8, false },
    { 0.98526631574995194, 0.8311342999091903, gsStatus::Success, 6, false },
    { 0.97105102178056768, 0.86227007797653632, gsStatus::Success, 5, false },
    { 0.95142925768646214, 0.89414870442323524, gsStatus::Success, 5, false },
    { 0.92626518778100431, 0.92643782806549779, gsStatus::Success, 5, false },
    { 0.89550215989137116, 0.95896454576061474, gsStatus::Success, 5, false },
    { 0.85914934666065113, 0.99162924923438811, gsStatus::Success, 5, false },
    { 0.81727755979900407, 1.0243710274892386, gsStatus::Success, 4, true },
    { 0.77001804635297288, 1.0571519496634896, gsStatus::Success, 4, false },
    { 0.71756189388186986, 1.0899492205007713, gsStatus::Success, 4, false },
    { 0.66015899944649259, 1.1227509469853088, gsStatus::Success, 4, false },
    { 0.59811613709077627, 1.1555536681723872, gsStatus::Success, 4, false },
    { 0.53179397124916805, 1.1883607703121519, gsStatus::Success, 4, false },
    { 0.46160305704772653, 1.2211813404424063, gsStatus::Success, 4, false },
    { 0.38799898941429162, 1.2540292213816757, gsStatus::Success, 4, false },
    { 0.31147692648937042, 1.2869221419540753, gsStatus::Success, 4, false },
    { 0.23256572685884358, 1.3198808592478855, gsStatus::Success, 4, false },
    { 0.15182191348448906, 1.3529282873979274, gsStatus::Success, 4, false },
    { 0.069823621784810691, 1.3860886106454593, gsStatus::Success, 4, false },
};
static const T1Meta T1_ConsistentCrisfield_A0_Q0_meta = { 40, false, true };
// T1_ConsistentCrisfield_A0_Q1: ds=0.05 N=40 adaptive=0 quasi=1 boundary 2 / 21 steps, limit point passed 0, final status gsStatus::NotConverged at lambda 0.998549
static const std::vector<T1Step> T1_ConsistentCrisfield_A0_Q1 = {
    { 0.080311815895446226, 0.035604821679763893, gsStatus::Success, 7, false },
    { 0.15785626348948922, 0.071021319925430093, gsStatus::Success, 4, false },
    { 0.23344565986427313, 0.10667278697148895, gsStatus::Success, 5, true },
    { 0.30650037102178018, 0.14235399619216516, gsStatus::Success, 5, false },
    { 0.37687356582325193, 0.17805670644251345, gsStatus::Success, 5, false },
    { 0.4444389974302414, 0.21378664171152642, gsStatus::Success, 5, false },
    { 0.50905605847611191, 0.24954712745636631, gsStatus::Success, 5, false },
    { 0.5705750390936577, 0.28534195446508798, gsStatus::Success, 5, false },
    { 0.62883744912540052, 0.32117604015438245, gsStatus::Success, 5, false },
    { 0.6836755560719131, 0.35705582402761193, gsStatus::Success, 5, false },
    { 0.73491200213050112, 0.39298992056348508, gsStatus::Success, 5, false },
    { 0.78235943279673081, 0.42899015256326378, gsStatus::Success, 5, false },
    { 0.82582007094543064, 0.46507322337389778, gsStatus::Success, 5, false },
    { 0.86508512482514699, 0.50126355296850778, gsStatus::Success, 5, false },
    { 0.89993377473387925, 0.53759838302206175, gsStatus::Success, 5, true },
    { 0.93013113951465998, 0.57413770724945856, gsStatus::Success, 6, false },
    { 0.95542370554900613, 0.61098571774119448, gsStatus::Success, 6, false },
    { 0.97552777398397084, 0.64834449942556194, gsStatus::Success, 6, false },
    { 0.99009418864696974, 0.68668223676913998, gsStatus::Success, 7, false },
    { 0.99854856820694649, 0.72750202097344852, gsStatus::Success, 9, false },
    { 0.99854856820694649, 0.72750202097344852, gsStatus::NotConverged, 99, false },
};
static const T1Meta T1_ConsistentCrisfield_A0_Q1_meta = { 21, false, true };
// T1_ConsistentCrisfield_A1_Q0: ds=0.05 N=40 adaptive=1 quasi=0 boundary 0 / 20 steps, limit point passed 0, final status gsStatus::Success at lambda 0.999881
// T1_ConsistentCrisfield_A1_Q0: fragile from step 20
// T1_ConsistentCrisfield_A1_Q0: adaptive prefix ends at step 20
static const std::vector<T1Step> T1_ConsistentCrisfield_A1_Q0 = {
    { 0.080311815895432251, 0.035604821679769909, gsStatus::Success, 5, false },
    { 0.14280345254770851, 0.064059748432613015, gsStatus::Success, 4, false },
    { 0.20375206427395956, 0.09252460771559555, gsStatus::Success, 4, false },
    { 0.26309846408409482, 0.12099990221337967, gsStatus::Success, 4, false },
    { 0.32078060778039458, 0.14948620618082417, gsStatus::Success, 4, false },
    { 0.3767334701668511, 0.17798418474276881, gsStatus::Success, 4, false },
    { 0.43088892263523654, 0.20649461932590213, gsStatus::Success, 4, false },
    { 0.48317561372104939, 0.23501844140174408, gsStatus::Success, 4, false },
    { 0.53351885443011904, 0.26355677761148905, gsStatus::Success, 4, false },
    { 0.58184051027979233, 0.2921110106661538, gsStatus::Success, 4, false },
    { 0.62805890200472969, 0.32068286241771388, gsStatus::Success, 4, false },
    { 0.67208871662262282, 0.3492745085981912, gsStatus::Success, 4, false },
    { 0.71384092982734715, 0.37788873965372938, gsStatus::Success, 4, false },
    { 0.7532227391062295, 0.4065291901732821, gsStatus::Success, 4, false },
    { 0.79013750390321036, 0.43520067308459981, gsStatus::Success, 4, false },
    { 0.82448468334769731, 0.46390967887018542, gsStatus::Success, 4, false },
    { 0.85615975117402365, 0.49266514443293563, gsStatus::Success, 4, false },
    { 0.88505404663066534, 0.52147968251714083, gsStatus::Success, 4, false },
    { 0.91105447882749035, 0.55037164141731443, gsStatus::Success, 4, false },
    { 0.93404291484690005, 0.57936876553081074, gsStatus::Success, 4, false },
    { 0.95389488213103402, 0.60851521833115174, gsStatus::Success, 5, true },
    { 0.9673195424098221, 0.63178243813560286, gsStatus::Success, 5, false },
    { 0.97644495071171766, 0.65035130230760374, gsStatus::Success, 5, false },
    { 0.98269180724257132, 0.66516714513181252, gsStatus::Success, 5, false },
    { 0.98700539536318221, 0.67698625075327246, gsStatus::Success, 5, false },
    { 0.99001404956620143, 0.68641337820414239, gsStatus::Success, 4, true },
    { 0.99264482350520156, 0.69587855592286441, gsStatus::Success, 4, true },
    { 0.99489224420757816, 0.70539711249815051, gsStatus::Success, 5, false },
    { 0.99639712784705103, 0.71299634245201127, gsStatus::Success, 5, false },
    { 0.99741245476338058, 0.71906207929825861, gsStatus::Success, 5, false },
    { 0.99810358610292182, 0.72390281818087865, gsStatus::Success, 5, false },
    { 0.99857879772863889, 0.72776523274439031, gsStatus::Success, 5, false },
    { 0.99890919041565507, 0.7308465597124808, gsStatus::Success, 4, true },
    { 0.9991973566875314, 0.73394130438806937, gsStatus::Success, 4, true },
    { 0.99944295434765318, 0.73705449831505598, gsStatus::Success, 5, true },
    { 0.99960708287094568, 0.73954043984799223, gsStatus::Success, 5, false },
    { 0.99971767472025652, 0.74152506255428041, gsStatus::Success, 5, false },
    { 0.99979289730209886, 0.74310910523234119, gsStatus::Success, 5, true },
    { 0.99984460212547677, 0.74437315567723517, gsStatus::Success, 5, true },
    { 0.99988055025779166, 0.74538167698984459, gsStatus::Success, 4, true },
};
static const T1Meta T1_ConsistentCrisfield_A1_Q0_meta = { 20, true, true };
// T1_ConsistentCrisfield_A1_Q1: ds=0.05 N=40 adaptive=1 quasi=1 boundary 0 / 17 steps, limit point passed 0, final status gsStatus::Success at lambda 0.973598
// T1_ConsistentCrisfield_A1_Q1: fragile from step 17
// T1_ConsistentCrisfield_A1_Q1: adaptive prefix ends at step 17
static const std::vector<T1Step> T1_ConsistentCrisfield_A1_Q1 = {
    { 0.080311815895446226, 0.035604821679763893, gsStatus::Success, 7, false },
    { 0.12476993369373189, 0.055776038894288164, gsStatus::Success, 4, false },
    { 0.16871424372058061, 0.076070183725724327, gsStatus::Success, 4, false },
    { 0.21186232255635915, 0.096369863461395297, gsStatus::Success, 4, false },
    { 0.25418604650763715, 0.11667219673920028, gsStatus::Success, 4, false },
    { 0.29566512406171191, 0.13697830635279812, gsStatus::Success, 4, false },
    { 0.33627637095379259, 0.15728833361084299, gsStatus::Success, 4, false },
    { 0.37599572395750003, 0.17760238480860169, gsStatus::Success, 4, false },
    { 0.41479834051645598, 0.19792059261060288, gsStatus::Success, 4, false },
    { 0.45265856100658114, 0.21824311165639462, gsStatus::Success, 4, false },
    { 0.48954988852609088, 0.2385701231623901, gsStatus::Success, 4, false },
    { 0.52544497148832947, 0.25890184166587082, gsStatus::Success, 4, false },
    { 0.56031558754062638, 0.27923852340760946, gsStatus::Success, 4, false },
    { 0.59413262921771104, 0.29958047690199452, gsStatus::Success, 4, false },
    { 0.62686609163267704, 0.31992807634734188, gsStatus::Success, 4, false },
    { 0.65848506247840821, 0.34028177873643828, gsStatus::Success, 4, false },
    { 0.68895771458451782, 0.3606421458357465, gsStatus::Success, 4, false },
    { 0.71825130121427228, 0.38100987264177266, gsStatus::Success, 4, true },
    { 0.74633215417610399, 0.40138582456252248, gsStatus::Success, 5, true },
    { 0.76785222467692094, 0.4176560581268402, gsStatus::Success, 4, true },
    { 0.78855403095184706, 0.43392838642093085, gsStatus::Success, 4, false },
    { 0.80842232816896775, 0.45020566388324007, gsStatus::Success, 4, false },
    { 0.82743828605794856, 0.46648845740552775, gsStatus::Success, 4, true },
    { 0.84558286229444435, 0.48277759888257876, gsStatus::Success, 5, true },
    { 0.85942394782829901, 0.49578138464454963, gsStatus::Success, 4, true },
    { 0.87268411544054258, 0.50878680840239809, gsStatus::Success, 4, false },
    { 0.88535549219693777, 0.52179649398081929, gsStatus::Success, 4, false },
    { 0.89742767783848532, 0.53481094748851876, gsStatus::Success, 4, true },
    { 0.90889036947974466, 0.54783104870162946, gsStatus::Success, 5, true },
    { 0.91759107531545059, 0.5582230501495411, gsStatus::Success, 4, false },
    { 0.92588942769271709, 0.56861685628842129, gsStatus::Success, 4, false },
    { 0.93378161756702371, 0.57901492170862723, gsStatus::Success, 4, false },
    { 0.9412620368489758, 0.58941773881128767, gsStatus::Success, 5, true },
    { 0.94693037652222278, 0.59772021944820475, gsStatus::Success, 4, true },
    { 0.9523288819653466, 0.60602441249037098, gsStatus::Success, 4, false },
    { 0.95745579279066806, 0.6143323291437035, gsStatus::Success, 4, false },
    { 0.96230814579524748, 0.62264434764454102, gsStatus::Success, 4, false },
    { 0.96688314389255126, 0.63096127212972153, gsStatus::Success, 5, false },
    { 0.97033098445122046, 0.63759811139879441, gsStatus::Success, 4, false },
    { 0.97359815113255177, 0.64423713592380105, gsStatus::Success, 4, false },
};
static const T1Meta T1_ConsistentCrisfield_A1_Q1_meta = { 17, true, true };

    TEST(T1_truss_selfcheck)
    {
        const gsStage0Truss<real_t> truss;
        const real_t ws[4] = { 0.1, truss.wLimit(), 0.8, 1.5 };
        for (index_t i = 0; i < 4; ++i)
        {
            gsVector<real_t> R;
            CHECK(truss.residual(truss.pathPoint(ws[i]), truss.lambdaOf(ws[i]), R));
            CHECK(R.norm() <= 1e-12);
        }

        // Tangent against a central finite difference of the residual
        const real_t eps = 1e-6;
        for (index_t i = 0; i < 4; ++i)
        {
            const gsVector<real_t> U = truss.pathPoint(ws[i]);
            gsSparseMatrix<real_t> K;
            CHECK(truss.jacobianMatrix(U, K));
            for (index_t j = 0; j < 2; ++j)
            {
                gsVector<real_t> Up = U, Um = U, Rp, Rm;
                Up[j] += eps;
                Um[j] -= eps;
                truss.residual(Up, 0.3, Rp);
                truss.residual(Um, 0.3, Rm);
                const gsVector<real_t> fd = (Rp - Rm)/(2*eps);
                for (index_t r = 0; r < 2; ++r)
                    CHECK_CLOSE(K.coeff(r,j), fd[r], 1e-6);
            }
        }

        CHECK_CLOSE(truss.lambdaMax(), truss.lambdaOf(truss.wLimit()), 1e-14);
        CHECK_CLOSE(0.0, truss.dfint(truss.wLimit()), 1e-12);
    }

    TEST(T1_LoadControl_reference)
    {
        UNITTEST_TIME_CONSTRAINT(1000);
        t1Gate<gsALMLoadControl<real_t> >("T1_LoadControl_A0_Q0", T1_LoadControl_A0_Q0, T1_LoadControl_A0_Q0_meta, LoadControl, false, false);
        t1Gate<gsALMLoadControl<real_t> >("T1_LoadControl_A0_Q1", T1_LoadControl_A0_Q1, T1_LoadControl_A0_Q1_meta, LoadControl, false, true);
        t1Gate<gsALMLoadControl<real_t> >("T1_LoadControl_A1_Q0", T1_LoadControl_A1_Q0, T1_LoadControl_A1_Q0_meta, LoadControl, true,  false);
        t1Gate<gsALMLoadControl<real_t> >("T1_LoadControl_A1_Q1", T1_LoadControl_A1_Q1, T1_LoadControl_A1_Q1_meta, LoadControl, true,  true);
    }

    TEST(T1_Riks_reference)
    {
        UNITTEST_TIME_CONSTRAINT(1000);
        t1Gate<gsALMRiks<real_t> >("T1_Riks_A0_Q0", T1_Riks_A0_Q0, T1_Riks_A0_Q0_meta, Riks, false, false);
        t1Gate<gsALMRiks<real_t> >("T1_Riks_A0_Q1", T1_Riks_A0_Q1, T1_Riks_A0_Q1_meta, Riks, false, true);
        t1Gate<gsALMRiks<real_t> >("T1_Riks_A1_Q0", T1_Riks_A1_Q0, T1_Riks_A1_Q0_meta, Riks, true,  false);
        t1Gate<gsALMRiks<real_t> >("T1_Riks_A1_Q1", T1_Riks_A1_Q1, T1_Riks_A1_Q1_meta, Riks, true,  true);
    }

    TEST(T1_Crisfield_reference)
    {
        UNITTEST_TIME_CONSTRAINT(1000);
        t1Gate<gsALMCrisfield<real_t> >("T1_Crisfield_A0_Q0", T1_Crisfield_A0_Q0, T1_Crisfield_A0_Q0_meta, Crisfield, false, false);
        t1Gate<gsALMCrisfield<real_t> >("T1_Crisfield_A0_Q1", T1_Crisfield_A0_Q1, T1_Crisfield_A0_Q1_meta, Crisfield, false, true);
        t1Gate<gsALMCrisfield<real_t> >("T1_Crisfield_A1_Q0", T1_Crisfield_A1_Q0, T1_Crisfield_A1_Q0_meta, Crisfield, true,  false);
        t1Gate<gsALMCrisfield<real_t> >("T1_Crisfield_A1_Q1", T1_Crisfield_A1_Q1, T1_Crisfield_A1_Q1_meta, Crisfield, true,  true);
    }

    TEST(T1_ConsistentCrisfield_reference)
    {
        UNITTEST_TIME_CONSTRAINT(1000);
        t1Gate<gsALMConsistentCrisfield<real_t> >("T1_ConsistentCrisfield_A0_Q0", T1_ConsistentCrisfield_A0_Q0, T1_ConsistentCrisfield_A0_Q0_meta, ConsistentCrisfield, false, false);
        t1Gate<gsALMConsistentCrisfield<real_t> >("T1_ConsistentCrisfield_A0_Q1", T1_ConsistentCrisfield_A0_Q1, T1_ConsistentCrisfield_A0_Q1_meta, ConsistentCrisfield, false, true);
        t1Gate<gsALMConsistentCrisfield<real_t> >("T1_ConsistentCrisfield_A1_Q0", T1_ConsistentCrisfield_A1_Q0, T1_ConsistentCrisfield_A1_Q0_meta, ConsistentCrisfield, true,  false);
        t1Gate<gsALMConsistentCrisfield<real_t> >("T1_ConsistentCrisfield_A1_Q1", T1_ConsistentCrisfield_A1_Q1, T1_ConsistentCrisfield_A1_Q1_meta, ConsistentCrisfield, true,  true);
    }

    /// Steps a fresh non-adaptive, non-quasi solver of type S through a
    /// trajectory and checks that every step is Success and lies on the
    /// analytic equilibrium path
    template<class S>
    std::vector<T1Rec> t1RunOnPath(const T1Preset preset, const bool requireAllSuccess)
    {
        const T1Params p = t1Params(preset, false, false);
        const gsStage0Truss<real_t> truss;
        std::vector<T1Rec> rec = t1Run<S>(p.ds, p.N, false, false);
        for (size_t i = 0; i < rec.size(); ++i)
        {
            if (requireAllSuccess)
                CHECK(rec[i].step.status == gsStatus::Success);
            if (rec[i].step.status == gsStatus::Success)
                CHECK(math::abs(rec[i].step.L - truss.lambdaOf(rec[i].w))
                      <= 3*T1_tol*math::max(real_t(1), math::abs(rec[i].step.L)));
        }
        return rec;
    }

    /// Load control: steps below lambda_max are on the path; the first step
    /// whose target exceeds lambda_max either fails and leaves the state
    /// unchanged, or jumps to the far stable branch
    template<class S>
    void t1LoadControlLimit()
    {
        const T1Params p = t1Params(LoadControl, false, false);
        const gsStage0Truss<real_t> truss;
        T1Setup<S> setup(p.ds, false, false, 1);

        real_t Lprev = 0, wprev = 0;
        bool crossed = false;
        for (index_t k = 0; k < p.N && !crossed; ++k)
        {
            const real_t target = (k + 1)*p.ds;
            const gsVector<real_t> Ubefore = setup.solver.solutionU();
            const real_t Lbefore = setup.solver.solutionL();
            const gsStatus st = setup.solver.step();
            if (target < truss.lambdaMax())
            {
                CHECK(st == gsStatus::Success);
                CHECK(math::abs(setup.solver.solutionL() - truss.lambdaOf(setup.solver.solutionU()[0]))
                      <= 3*T1_tol*math::max(real_t(1), math::abs(setup.solver.solutionL())));
                CHECK(setup.solver.solutionL() >= Lprev);
            }
            else
            {
                crossed = true;
                const real_t jump = setup.solver.solutionU()[0] - wprev;
                const bool failedUnchanged = st != gsStatus::Success
                    && setup.solver.solutionL() == Lbefore
                    && (setup.solver.solutionU() - Ubefore).norm() == 0;
                const bool jumped = st == gsStatus::Success && jump > 0.5;
                gsInfo << "T1 load control crossing step " << k << ": status "
                       << t1StatusName(st) << (failedUnchanged ? " (state unchanged)" : "")
                       << ", w jump " << jump << "\n";
                CHECK(failedUnchanged || jumped);
            }
            Lprev = setup.solver.solutionL();
            wprev = setup.solver.solutionU()[0];
        }
        CHECK(crossed);
    }

    TEST(T1_LoadControl_fails_at_limit_point)
    {
        t1LoadControlLimit<gsALMLoadControl<real_t> >();
    }

    /// Riks/Crisfield-type solvers pass the limit point: all steps succeed, lie
    /// on the path, reach lambda_max and later have a decreasing lambda
    template<class S>
    void t1PassesLimit(const T1Preset preset)
    {
        const gsStage0Truss<real_t> truss;
        const std::vector<T1Rec> rec = t1RunOnPath<S>(preset, true);
        real_t Lmax = 0;
        for (size_t i = 0; i < rec.size(); ++i)
            Lmax = math::max(Lmax, rec[i].step.L);
        CHECK(Lmax >= 0.98*truss.lambdaMax());
        CHECK(t1PassedLimit(rec, truss.wLimit()));
    }

    TEST(T1_Riks_passes_limit_point)
    {
        t1PassesLimit<gsALMRiks<real_t> >(Riks);
    }

    TEST(T1_Crisfield_passes_limit_point)
    {
        t1PassesLimit<gsALMCrisfield<real_t> >(Crisfield);
    }

    /// One converged step of a stability run: apex deflection, analytic
    /// stability sign from f'(w), and what the solver reports
    struct StabRec
    {
        real_t  w;          ///< apex deflection U[0]
        real_t  dfint;      ///< f'(w); K is positive definite iff f'(w) > 0
        index_t expected;   ///< +1 stable (f' > 0), -1 unstable
        index_t stability;  ///< solver.stability() after computeStability(true)
        bool    change;     ///< solver.stabilityChange()
        real_t  indicator;  ///< solver.indicator() (smallest eigenvalue or LDL^T pivot)
        gsStatus stepStatus, stabStatus;
    };

    /// Riks (full Newton, fixed arc length 0.05) on the T1 truss with the given
    /// BifurcationMethod; computeStability(true) at every converged point.
    /// Stops after the first failed step or the first step with w > 1, which lies
    /// between the limit point w* = 1 - 1/sqrt(3) and the second critical point
    /// w2 = 1 + 1/sqrt(3), so the window contains exactly one stability change.
    /// \a initStability receives solver.stability() after initialize() (U = 0).
    static std::vector<StabRec> stabRun(index_t bifMethod, index_t & initStability)
    {
        typedef gsStructuralAnalysisOps<real_t>::Jacobian_t   Jacobian_t;
        typedef gsStructuralAnalysisOps<real_t>::ALResidual_t ALResidual_t;

        gsStage0Truss<real_t> truss;
        Jacobian_t       J = truss.jacobian();
        ALResidual_t     R = truss.alResidual();
        gsVector<real_t> F = truss.force();
        gsALMRiks<real_t> riks(J, R, F);
        // gsALMRiks re-declares computeStability() as protected; the public
        // interface is that of the base class
        gsALMBase<real_t> & solver = riks;

        // setOptions() re-reads "Length", so all options are set at once
        gsOptionList opts = solver.options();
        opts.setReal("Length", 0.05);
        opts.setSwitch("AdaptiveLength", false);
        opts.setSwitch("Quasi", false);
        opts.setReal("TolF", 1e-10);
        opts.setReal("TolU", 1e-10);
        opts.setInt("BifurcationMethod", bifMethod);
        solver.setOptions(opts);
        solver.initialize();
        initStability = solver.stability();

        std::vector<StabRec> rec;
        for (index_t k = 0; k < 40; ++k)
        {
            StabRec r;
            r.stepStatus = solver.step();
            r.stabStatus = solver.computeStability(true);
            r.w          = solver.solutionU()[0];
            r.dfint      = truss.dfint(r.w);
            r.expected   = (r.dfint > 0) ? 1 : -1;
            r.stability  = solver.stability();
            r.change     = solver.stabilityChange();
            r.indicator  = solver.indicator();
            rec.push_back(r);
            if (r.stepStatus != gsStatus::Success || r.w > 1)
                break;
        }
        return rec;
    }

    static void stabSignCheck(index_t bifMethod, const char * name)
    {
        index_t init = 0;
        const std::vector<StabRec> rec = stabRun(bifMethod, init);
        CHECK_EQUAL(1, init);
        for (size_t i = 0; i < rec.size(); ++i)
            CHECK_EQUAL(rec[i].expected, rec[i].stability);
        gsInfo << name << ": sign checks " << rec.size() + 1 << "\n";
    }

    static void stabChangeCheck(index_t bifMethod)
    {
        index_t init = 0;
        const std::vector<StabRec> rec = stabRun(bifMethod, init);
        CHECK(!rec.empty());

        index_t nPlus = 0, nMinus = 0, nFired = 0, jFired = -1, jExpected = -1;
        index_t prev = 1; // U = 0 is stable: f'(0) = 2c > 0
        for (size_t i = 0; i < rec.size(); ++i)
        {
            CHECK(rec[i].stepStatus == gsStatus::Success);
            CHECK(rec[i].stabStatus == gsStatus::Success);
            CHECK(math::abs(rec[i].dfint) > 1e-6);
            CHECK((rec[i].indicator > 0) == (rec[i].expected == 1));
            (rec[i].expected == 1 ? nPlus : nMinus)++;
            if (rec[i].expected != prev)
                jExpected = static_cast<index_t>(i);
            prev = rec[i].expected;
            if (rec[i].change)
            {
                ++nFired;
                jFired = static_cast<index_t>(i);
            }
        }
        CHECK(nPlus >= 3);
        CHECK(nMinus >= 3);
        if (!rec.empty())
            CHECK(rec.back().w > 1);
        CHECK_EQUAL(1, nFired);
        CHECK_EQUAL(jExpected, jFired);
    }

    TEST(StabilitySign_Riks_Eigenvalue)
    {
        stabSignCheck(gsALMBase<real_t>::bifmethod::Eigenvalue, "StabilitySign_Riks_Eigenvalue");
    }

    TEST(StabilitySign_Riks_Determinant)
    {
        stabSignCheck(gsALMBase<real_t>::bifmethod::Determinant, "StabilitySign_Riks_Determinant");
    }

    TEST(StabilitySign_ChangeOnce_Riks_Eigenvalue)
    {
        stabChangeCheck(gsALMBase<real_t>::bifmethod::Eigenvalue);
    }

    TEST(StabilitySign_ChangeOnce_Riks_Determinant)
    {
        stabChangeCheck(gsALMBase<real_t>::bifmethod::Determinant);
    }
}
