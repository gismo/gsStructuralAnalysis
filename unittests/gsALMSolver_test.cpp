/** @file gsALMSolver_test.cpp

    @brief Regression recordings of the four arc-length presets (load control,
    Riks, Crisfield, consistent Crisfield) on the two-dof von Mises truss of
    gsStage0Truss.h.

    * Reference gates: for each preset and each combination of AdaptiveLength
      and Quasi, the per-step (lambda, |U|, status, iterations) trajectory is
      compared with a recording stored in gsStage0TrussReferences.h. Steps
      whose iteration count is not robust against a change of the convergence
      tolerance are marked as boundary steps; their iteration count is not
      gated.
    * Negative controls: load control cannot pass the limit point; Riks,
      Crisfield and consistent Crisfield do.

    Setting the environment variable GSSA_STAGE0_PRINT prints all reference
    arrays as C++ initialisers for gsStage0TrussReferences.h (through gsInfo)
    and skips the comparison with the stored arrays.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
 **/

#include "gismo_unittest.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <algorithm>
#include <string>
#include <vector>

#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMLoadControl.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMRiks.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMCrisfield.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMConsistentCrisfield.h>
#include "gsStage0Truss.h"
#include "gsStage0TrussReferences.h"

SUITE(gsALMSolver_test)
{
    using namespace gismo::stage0;

    /// A recorded step together with the apex deflection w = U[0] and the
    /// load factor the step aimed at (previous lambda plus the arc length;
    /// for load control the target of the load increment)
    struct T1Rec
    {
        T1Step step;
        real_t w;
        real_t target;
    };

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

    /// Consistent Crisfield solver exposing its scaling factor phi
    class T1CCProbe : public gsALMConsistentCrisfield<real_t>
    {
    public:
        using gsALMConsistentCrisfield<real_t>::gsALMConsistentCrisfield;
        real_t phi() const { return this->m_cc->scaling(); }
    };

    /// Steps a fresh non-adaptive solver of type S with arc length ds for N
    /// steps; checks that every step is Success (when requested) and, for a
    /// non-quasi run, that every successful step lies on the analytic path
    template<class S>
    std::vector<T1Rec> t1RunOnPathAt(real_t ds, index_t N, bool quasi, bool requireAllSuccess)
    {
        const gsStage0Truss<real_t> truss;
        std::vector<T1Rec> rec = t1Run<S>(ds, N, false, quasi);
        for (size_t i = 0; i < rec.size(); ++i)
        {
            if (requireAllSuccess)
                CHECK(rec[i].step.status == gsStatus::Success);
            if (!quasi && rec[i].step.status == gsStatus::Success)
                CHECK(math::abs(rec[i].step.L - truss.lambdaOf(rec[i].w))
                      <= 3*T1_tol*math::max(real_t(1), math::abs(rec[i].step.L)));
        }
        return rec;
    }

    /// All N steps succeed, lambda_max is reached and the limit point is passed
    template<class S>
    void t1PassesLimitAt(real_t ds, index_t N, bool quasi)
    {
        const gsStage0Truss<real_t> truss;
        const std::vector<T1Rec> rec = t1RunOnPathAt<S>(ds, N, quasi, true);
        CHECK_EQUAL(size_t(N), rec.size());
        real_t Lmax = 0;
        for (size_t i = 0; i < rec.size(); ++i)
            Lmax = math::max(Lmax, rec[i].step.L);
        CHECK(Lmax >= 0.98*truss.lambdaMax());
        CHECK(t1PassedLimit(rec, truss.wLimit()));
    }

    TEST(T1_ConsistentCrisfield_passes_limit_point)
    {
        UNITTEST_TIME_CONSTRAINT(1000);
        t1PassesLimitAt<gsALMConsistentCrisfield<real_t> >(0.1,  40, false);
        t1PassesLimitAt<gsALMConsistentCrisfield<real_t> >(0.05, 40, true);
    }

    /// Tangent ratio |u_t|/|F| with u_t = K(U)^{-1} F, from a dense solve
    static real_t t1TangentPhi(const gsStage0Truss<real_t> & truss, const gsVector<real_t> & F, const gsVector<real_t> & U)
    {
        gsSparseMatrix<real_t> K;
        truss.jacobian()(U, K);
        const gsMatrix<real_t> Kd = K.toDense();
        const gsVector<real_t> ut = Kd.partialPivLu().solve(F);
        return std::sqrt(ut.squaredNorm() / F.squaredNorm());
    }

    /// With automatic scaling, phi is one constant per run: the tangent ratio
    /// |u_t|/|F| at the start point, whatever the step does afterwards
    TEST(T1_ConsistentCrisfield_scaling_is_frozen)
    {
        UNITTEST_TIME_CONSTRAINT(1000);
        struct Cell { real_t ds; bool adaptive; bool quasi; };
        const Cell cells[] = { {0.05, false, false}, {0.05, false, true},
                               {0.05, true,  false}, {0.05, true,  true},
                               {0.1,  false, false} };
        index_t checked = 0;
        for (size_t c = 0; c < sizeof(cells)/sizeof(cells[0]); ++c)
        {
            T1Setup<T1CCProbe> setup(cells[c].ds, cells[c].adaptive, cells[c].quasi, 1);
            real_t phi0 = 0;
            for (index_t k = 0; k < 40; ++k)
            {
                const gsStatus st = setup.solver.step();
                if (st != gsStatus::Success)
                    break;
                if (k == 0)
                {
                    phi0 = setup.solver.phi();
                    const real_t tangent = t1TangentPhi(setup.truss, setup.F, gsVector<real_t>::Zero(2));
                    CHECK_CLOSE(tangent, phi0, 1e-12*phi0);
                }
                else
                {
                    CHECK(setup.solver.phi() == phi0);
                    ++checked;
                }
            }
        }
        CHECK(checked > 40);
    }

    /// The path crosses lambda = 0 forward at w > 1: every step is a Success
    /// with w increasing (no reversal), also for the five steps after the
    /// crossing
    TEST(T1_ConsistentCrisfield_crosses_zero_load_forward)
    {
        UNITTEST_TIME_CONSTRAINT(1000);
        const real_t steps[] = { 0.2, 0.4 };
        for (size_t c = 0; c < 2; ++c)
            for (index_t q = 0; q < 2; ++q)
            {
                const bool quasi = (q == 1);
                T1Setup<T1CCProbe> setup(steps[c], false, quasi, 1);
                const gsStage0Truss<real_t> & truss = setup.truss;
                real_t wprev = 0;
                index_t crossing = -1, after = 0;
                for (index_t k = 0; k < 30; ++k)
                {
                    const gsStatus st = setup.solver.step();
                    CHECK(st == gsStatus::Success);
                    if (st != gsStatus::Success)
                        break;
                    const real_t L = setup.solver.solutionL();
                    const real_t w = setup.solver.solutionU()[0];
                    CHECK(w > wprev);
                    if (!quasi)
                        CHECK(math::abs(L - truss.lambdaOf(w)) <= 3*T1_tol*math::max(real_t(1), math::abs(L)));
                    wprev = w;
                    if (crossing < 0)
                    {
                        if (L < 0 && w > 1)
                            crossing = k;
                    }
                    else if (++after == 5)
                        break;
                }
                gsInfo << "T1 zero-load crossing: ds " << steps[c] << ", quasi " << quasi
                       << ", crossing step " << crossing << ", steps after " << after << "\n";
                CHECK(crossing >= 0);
                CHECK_EQUAL(5, after);
            }
    }

    /// A mid-run applyOptions() resets the Scaling option, not the scaling of
    /// the run: the trajectory is unchanged
    TEST(T1_ConsistentCrisfield_scaling_survives_applyOptions)
    {
        UNITTEST_TIME_CONSTRAINT(1000);
        const index_t N = 40;
        T1Setup<T1CCProbe> R(0.05, false, false, 1);
        std::vector<T1Step> ref;
        for (index_t k = 0; k < N; ++k)
        {
            T1Step s;
            s.status     = R.solver.step();
            s.L          = R.solver.solutionL();
            s.normU      = R.solver.solutionU().norm();
            s.iterations = R.solver.numIterations();
            s.boundary   = false;
            ref.push_back(s);
        }

        T1Setup<T1CCProbe> X(0.05, false, false, 1);
        real_t phi0 = 0;
        index_t nSame = 0;
        for (index_t k = 0; k < N; ++k)
        {
            if (k == 10)
                X.solver.applyOptions();
            const gsStatus st = X.solver.step();
            CHECK(st == ref[k].status);
            CHECK_CLOSE(ref[k].L,     X.solver.solutionL(),         1e-14*math::abs(ref[k].L));
            CHECK_CLOSE(ref[k].normU, X.solver.solutionU().norm(),  1e-14*math::abs(ref[k].normU));
            if (st == gsStatus::Success)
            {
                if (k == 0)
                    phi0 = X.solver.phi();
                CHECK(X.solver.phi() == phi0);
                ++nSame;
            }
        }
        CHECK(phi0 > 0);
        CHECK(nSame > 20);
    }

    /// A run that does not start at the origin takes its scaling from the
    /// first step of that run: after setSolution on a fresh solver, after
    /// initialize() on a used solver, and through an initial guess
    TEST(T1_ConsistentCrisfield_scaling_on_restart)
    {
        UNITTEST_TIME_CONSTRAINT(1000);
        T1Setup<T1CCProbe> A(0.05, false, false, 1);
        gsVector<real_t> U9, U10;
        real_t L9 = 0, L10 = 0, phi0A = 0;
        for (index_t k = 0; k < 10; ++k)
        {
            A.solver.step();
            if (k == 0)
                phi0A = A.solver.phi();
            if (k == 8)
            {
                U9 = A.solver.solutionU();
                L9 = A.solver.solutionL();
            }
        }
        U10 = A.solver.solutionU();
        L10 = A.solver.solutionL();
        const real_t phiRestart = t1TangentPhi(A.truss, A.F, U10);

        // (b1) fresh solver restarted at step 10
        T1Setup<T1CCProbe> Y(0.05, false, false, 1);
        Y.solver.setSolution(U10, L10);
        Y.solver.setPrevious(U9, L9);
        const gsStatus stY = Y.solver.step();
        CHECK_CLOSE(phiRestart, Y.solver.phi(), 1e-10*phiRestart);
        CHECK(Y.solver.phi() != phi0A);

        // (b2) the used solver, re-initialised and restarted at step 10
        A.solver.initialize();
        A.solver.setSolution(U10, L10);
        A.solver.setPrevious(U9, L9);
        const gsStatus stA = A.solver.step();
        CHECK(A.solver.phi() == Y.solver.phi());
        CHECK_CLOSE(Y.solver.solutionL(),        A.solver.solutionL(),        1e-14*math::abs(Y.solver.solutionL()));
        CHECK_CLOSE(Y.solver.solutionU().norm(), A.solver.solutionU().norm(), 1e-14*Y.solver.solutionU().norm());

        // (b3) fresh solver restarted through an initial guess
        T1Setup<T1CCProbe> Z(0.05, false, false, 1);
        Z.solver.setSolution(U10, L10);
        Z.solver.setPrevious(U9, L9);
        Z.solver.setInitialGuess(2*U10 - U9, 2*L10 - L9);
        const gsStatus stZ = Z.solver.step();
        CHECK_CLOSE(phiRestart, Z.solver.phi(), 1e-10*phiRestart);

        gsInfo << "T1 restart: phi0 of the run " << phi0A << ", tangent phi at restart " << phiRestart
               << "; statuses b1 " << t1StatusName(stY) << ", b2 " << t1StatusName(stA)
               << ", b3 " << t1StatusName(stZ) << "\n";
    }

    /// Putting the solver back at the start of its last converged step, as the
    /// drivers do for a redo and for computeSingularPoint(), keeps the scaling
    /// of the run
    TEST(T1_ConsistentCrisfield_scaling_survives_redo)
    {
        UNITTEST_TIME_CONSTRAINT(1000);
        {
            T1Setup<T1CCProbe> setup(0.05, false, false, 1);
            gsVector<real_t> Uold;
            real_t Lold = 0, phi0 = 0;
            for (index_t k = 0; k < 10; ++k)
            {
                Uold = setup.solver.solutionU();
                Lold = setup.solver.solutionL();
                setup.solver.step();
                if (k == 0)
                    phi0 = setup.solver.phi();
            }
            CHECK(phi0 > 0);
            setup.solver.setSolution(Uold, Lold);
            for (index_t j = 0; j < 2; ++j)
            {
                const gsStatus st = setup.solver.step();
                CHECK(st == gsStatus::Success);
                CHECK(setup.solver.phi() == phi0);
            }
        }
        {
            T1Setup<T1CCProbe> setup(0.05, false, false, 1);
            gsALMBase<real_t> & base = setup.solver;
            gsVector<real_t> Uold;
            real_t Lold = 0, phi0 = 0;
            bool found = false;
            for (index_t k = 0; k < 40 && !found; ++k)
            {
                Uold = setup.solver.solutionU();
                Lold = setup.solver.solutionL();
                const gsStatus st = setup.solver.step();
                if (k == 0)
                    phi0 = setup.solver.phi();
                if (st != gsStatus::Success)
                    break;
                base.computeStability(true);
                found = base.stabilityChange();
            }
            gsInfo << "T1 redo: stability change found " << found << "\n";
            CHECK(found);
            if (found)
            {
                base.computeSingularPoint(Uold, Lold, false);
                for (index_t j = 0; j < 3; ++j)
                {
                    setup.solver.step();
                    CHECK(setup.solver.phi() == phi0);
                }
            }
        }
    }

    /// The only exercise of predictorGuess(): an initial guess on the straight
    /// continuation of the last two converged points reproduces the plain run
    TEST(T1_ConsistentCrisfield_predictorGuess)
    {
        UNITTEST_TIME_CONSTRAINT(1000);
        typedef gsALMConsistentCrisfield<real_t> S;
        const real_t ds = 0.05;
        const index_t N = 40;
        const std::vector<T1Rec> ref = t1Run<S>(ds, N, false, true);

        T1Setup<S> setup(ds, false, true, 1);
        std::vector<T1Rec> rec;
        gsVector<real_t> Uprev = setup.solver.solutionU();
        real_t Lprev = setup.solver.solutionL();
        for (index_t k = 0; k < N; ++k)
        {
            const gsVector<real_t> U = setup.solver.solutionU();
            const real_t           L = setup.solver.solutionL();
            if (k > 0)
                setup.solver.setInitialGuess(2*U - Uprev, 2*L - Lprev);
            const gsStatus st = setup.solver.step();
            T1Rec r;
            r.step.L          = setup.solver.solutionL();
            r.step.normU      = setup.solver.solutionU().norm();
            r.step.status     = st;
            r.step.iterations = setup.solver.numIterations();
            r.step.boundary   = false;
            r.w               = setup.solver.solutionU()[0];
            r.target          = 0;
            rec.push_back(r);
            Uprev = U;
            Lprev = L;
            if (st != gsStatus::Success)
                break;
        }
        gsInfo << "T1 predictorGuess: plain run " << ref.size() << " steps, guess run "
               << rec.size() << " steps, last guess-run status "
               << t1StatusName(rec.back().step.status) << "\n";
        CHECK_EQUAL(size_t(N), ref.size());
        CHECK_EQUAL(size_t(N), rec.size());
        for (size_t i = 0; i < std::min(ref.size(), rec.size()); ++i)
        {
            CHECK(ref[i].step.status == gsStatus::Success);
            CHECK(rec[i].step.status == gsStatus::Success);
            CHECK_CLOSE(ref[i].step.L,     rec[i].step.L,     1e-8*math::max(math::abs(ref[i].step.L), real_t(1)));
            CHECK_CLOSE(ref[i].step.normU, rec[i].step.normU, 1e-8*math::max(math::abs(ref[i].step.normU), real_t(1)));
        }
        CHECK(t1PassedLimit(rec, gsStage0Truss<real_t>().wLimit()));
    }

    /// A user-set "Scaling" is never overwritten by the predictor
    TEST(T1_ConsistentCrisfield_user_scaling)
    {
        UNITTEST_TIME_CONSTRAINT(1000);
        typedef gsStructuralAnalysisOps<real_t>::Jacobian_t   Jacobian_t;
        typedef gsStructuralAnalysisOps<real_t>::ALResidual_t ALResidual_t;
        gsStage0Truss<real_t> truss;
        Jacobian_t       J = truss.jacobian();
        ALResidual_t     R = truss.alResidual();
        gsVector<real_t> F = truss.force();
        T1CCProbe solver(J, R, F);
        gsOptionList opts = solver.options();
        opts.setReal("Length", 0.05);
        opts.setSwitch("AdaptiveLength", false);
        opts.setInt("AdaptiveIterations", 4);
        opts.setSwitch("Quasi", false);
        opts.setReal("TolF", T1_tol);
        opts.setReal("TolU", T1_tol);
        opts.setReal("Scaling", 0.5);
        solver.setOptions(opts);
        solver.initialize();
        for (index_t k = 0; k < 40; ++k)
        {
            const gsStatus st = solver.step();
            CHECK_EQUAL(real_t(0.5), solver.phi());
            if (st != gsStatus::Success)
                break;
        }
    }

    /// Truss, its operators and a solver of type S for the acceptance-gate
    /// tests. A \a scaling > 0 sets "Scaling" (only Crisfield and consistent
    /// Crisfield have this option); a \a constraintTol >= 0 sets "ConstraintTol".
    /// AdaptiveLength and Quasi are off.
    template<class S>
    struct GateSetup
    {
        typedef typename gsStructuralAnalysisOps<real_t>::Jacobian_t   Jacobian_t;
        typedef typename gsStructuralAnalysisOps<real_t>::ALResidual_t ALResidual_t;

        gsStage0Truss<real_t> truss;
        Jacobian_t            J;
        ALResidual_t          R;
        gsVector<real_t>      F;
        S                     solver;

        GateSetup(real_t ds, real_t tolFU, index_t maxIter, real_t scaling = -1, real_t constraintTol = -1)
        : J(truss.jacobian()), R(truss.alResidual()), F(truss.force()), solver(J,R,F)
        {
            // setOptions() re-reads "Length", so all options are set at once
            gsOptionList opts = solver.options();
            opts.setReal("Length", ds);
            opts.setSwitch("AdaptiveLength", false);
            opts.setSwitch("Quasi", false);
            opts.setInt("MaxIter", maxIter);
            opts.setReal("TolF", tolFU);
            opts.setReal("TolU", tolFU);
            if (scaling > 0)
                opts.setReal("Scaling", scaling);
            if (constraintTol >= 0)
                opts.setReal("ConstraintTol", constraintTol);
            solver.setOptions(opts);
            solver.initialize();
        }
    };

    /// |g|/s^2 of the path constraint of \a preset, evaluated from the accepted
    /// increment (DU, DL) and the arc length s. Scaling is 1 for Crisfield and
    /// consistent Crisfield; Riks uses phi = 1/numDof = 1/2 on the 2-dof truss.
    static real_t gateResidual(const T1Preset preset, const real_t s, const gsVector<real_t> & DU, const real_t DL)
    {
        switch (preset)
        {
        case Riks:
            return math::abs(0.5*DU.dot(DU) + 0.5*DL*DL - s*s) / (s*s);
        case LoadControl:
            return math::abs(DL/s - 1);
        default:
            return math::abs(DU.dot(DU) + DL*DL - s*s) / (s*s);
        }
    }

    /// Three steps with acceptance tolerances so loose that every finite
    /// iterate passes TolF and TolU; each accepted step must still satisfy the
    /// path constraint to 1e-6.
    template<class S>
    void gateHolds(const T1Preset preset, const char * name, const real_t ds, const real_t scaling)
    {
        GateSetup<S> setup(ds, 10, 100, scaling);
        for (index_t k = 0; k < 3; ++k)
        {
            const gsStatus st = setup.solver.step();
            const real_t g2 = gateResidual(preset, setup.solver.getLength(),
                                           setup.solver.solutionDU(), setup.solver.solutionDL());
            char line[256];
            std::snprintf(line, sizeof(line), "%s step %d: status %d, iterations %d, |g|/s^2 = %.3e\n",
                          name, static_cast<int>(k), static_cast<int>(st),
                          static_cast<int>(setup.solver.numIterations()), static_cast<double>(g2));
            gsInfo << line;
            CHECK(st == gsStatus::Success);
            CHECK(g2 < 1e-6);
        }
    }

    TEST(ALMGate_LoadControl)
    {
        gateHolds<gsALMLoadControl<real_t> >(LoadControl, "ALMGate_LoadControl", 0.07, -1);
    }

    TEST(ALMGate_Riks)
    {
        gateHolds<gsALMRiks<real_t> >(Riks, "ALMGate_Riks", 0.05, -1);
    }

    TEST(ALMGate_Crisfield)
    {
        gateHolds<gsALMCrisfield<real_t> >(Crisfield, "ALMGate_Crisfield", 0.05, 1);
    }

    TEST(ALMGate_ConsistentCrisfield)
    {
        gateHolds<gsALMConsistentCrisfield<real_t> >(ConsistentCrisfield, "ALMGate_ConsistentCrisfield", 0.05, 1);
    }

    /// With TolF = TolU = 0 the strict residue test never holds, so a step
    /// must run exactly MaxIter corrector iterations and report NotConverged.
    template<class S>
    void maxIterIsIterationCount(const real_t ds)
    {
        const index_t budgets[] = {1, 5};
        for (index_t i = 0; i < 2; ++i)
        {
            GateSetup<S> setup(ds, 0, budgets[i]);
            CHECK(setup.solver.step() == gsStatus::NotConverged);
            CHECK_EQUAL(budgets[i], setup.solver.numIterations());
        }
    }

    TEST(ALMMaxIter_LoadControl)
    {
        maxIterIsIterationCount<gsALMLoadControl<real_t> >(0.07);
    }

    TEST(ALMMaxIter_Riks)
    {
        maxIterIsIterationCount<gsALMRiks<real_t> >(0.05);
    }

    TEST(ALMMaxIter_Crisfield)
    {
        maxIterIsIterationCount<gsALMCrisfield<real_t> >(0.05);
    }

    TEST(ALMMaxIter_ConsistentCrisfield)
    {
        maxIterIsIterationCount<gsALMConsistentCrisfield<real_t> >(0.05);
    }

    template<class S>
    void constraintTolDefault(const real_t ds)
    {
        GateSetup<S> setup(ds, 1e-3, 100);
        CHECK_EQUAL(1e-6, setup.solver.options().getReal("ConstraintTol"));
    }

    TEST(ALMGate_ConstraintTol_default)
    {
        constraintTolDefault<gsALMLoadControl<real_t> >(0.07);
        constraintTolDefault<gsALMRiks<real_t> >(0.05);
        constraintTolDefault<gsALMCrisfield<real_t> >(0.05);
        constraintTolDefault<gsALMConsistentCrisfield<real_t> >(0.05);
    }

    /// A huge ConstraintTol leaves only the residue tests in the gate: Riks
    /// then accepts after the first iteration although the predictor is off the
    /// constraint.
    TEST(ALMGate_ConstraintTol_honoured)
    {
        GateSetup<gsALMRiks<real_t> > setup(0.05, 10, 100, -1, 1e300);
        CHECK(setup.solver.step() == gsStatus::Success);
        CHECK_EQUAL(1, setup.solver.numIterations());
    }

    /// Pins the iteration count of the first Riks step on a problem with
    /// |u_t| = |K^{-1}F| << 1 (stiff truss, |u_t| = 4.4e-7): the first-step
    /// predictor must lie on the Riks constraint, so the corrector converges in a
    /// few iterations, and the converged step satisfies
    /// phi |U|^2 + (1-phi) lambda^2 = s^2 with phi = 1/2.
    TEST(T1_Riks_first_step_on_constraint)
    {
        typedef gsStructuralAnalysisOps<real_t>::Jacobian_t   Jacobian_t;
        typedef gsStructuralAnalysisOps<real_t>::ALResidual_t ALResidual_t;

        gsStage0Truss<real_t> truss(real_t(1.5) * math::sqrt(real_t(3)) * 1e6, 1, 5e6);
        Jacobian_t       J = truss.jacobian();
        ALResidual_t     R = truss.alResidual();
        gsVector<real_t> F = truss.force();
        gsALMRiks<real_t> solver(J, R, F);

        const real_t ds = 0.05;
        gsOptionList opts = solver.options();
        opts.setReal("Length", ds);
        opts.setSwitch("AdaptiveLength", false);
        opts.setSwitch("Quasi", false);
        opts.setReal("TolF", 1e-10);
        opts.setReal("TolU", 1e-10);
        solver.setOptions(opts);
        solver.initialize();

        const gsStatus st = solver.step();
        CHECK(st == gsStatus::Success);
        gsInfo << "T1_Riks_first_step_on_constraint: iterations "
               << solver.numIterations() << "\n";
        CHECK(solver.numIterations() <= 5);

        const gsVector<real_t> U = solver.solutionU();
        const real_t L = solver.solutionL();
        const real_t g = 0.5 * U.dot(U) + 0.5 * L * L;
        CHECK_CLOSE(ds * ds, g, 1e-8 * ds * ds);
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

    /// Runs \a n steps through the base interface; (lambda, |U|) after every
    /// successful step is appended to \a rec
    template<class S>
    void frozenScalingRun(T1Setup<S> & setup, index_t n, std::vector<std::pair<real_t,real_t> > & rec)
    {
        gsALMBase<real_t> & b = setup.solver;
        for (index_t k = 0; k < n; ++k)
        {
            if (b.step() != gsStatus::Success)
                break;
            rec.push_back(std::make_pair(b.solutionL(), b.solutionU().norm()));
        }
    }

    /// frozenScaling, hasFrozenScaling and setFrozenScaling are virtual on
    /// gsALMBase: consistent Crisfield overrides them, the other presets ignore them
    TEST(ALMFrozenScaling_base_hook)
    {
        UNITTEST_TIME_CONSTRAINT(1000);

        // Automatic scaling: frozen by the first step
        real_t p = 0;
        {
            T1Setup<gsALMConsistentCrisfield<real_t> > A(0.05, false, false, 1);
            gsALMBase<real_t> & b = A.solver;
            CHECK(!b.hasFrozenScaling());
            CHECK_EQUAL(real_t(0), b.frozenScaling());
            CHECK(b.step() == gsStatus::Success);
            p = b.frozenScaling();
            CHECK(b.hasFrozenScaling());
            CHECK(p > 0);
            CHECK_EQUAL(p, A.solver.frozenScaling());
        }

        // Imposed through the base interface: kept by the first step
        {
            T1Setup<gsALMConsistentCrisfield<real_t> > B(0.05, false, false, 1);
            gsALMBase<real_t> & b = B.solver;
            b.setFrozenScaling(real_t(1.5)*p);
            CHECK(B.solver.hasFrozenScaling());
            CHECK_EQUAL(real_t(1.5)*p, B.solver.frozenScaling());
            CHECK(b.step() == gsStatus::Success);
            CHECK_EQUAL(real_t(1.5)*p, b.frozenScaling());
        }

        // User scaling: never frozen, so nothing is exchanged
        {
            typedef gsStructuralAnalysisOps<real_t>::Jacobian_t   Jacobian_t;
            typedef gsStructuralAnalysisOps<real_t>::ALResidual_t ALResidual_t;
            gsStage0Truss<real_t> truss;
            Jacobian_t       J = truss.jacobian();
            ALResidual_t     R = truss.alResidual();
            gsVector<real_t> F = truss.force();
            gsALMConsistentCrisfield<real_t> solver(J, R, F);
            gsOptionList opts = solver.options();
            opts.setReal("Length", 0.05);
            opts.setSwitch("AdaptiveLength", false);
            opts.setInt("AdaptiveIterations", 4);
            opts.setSwitch("Quasi", false);
            opts.setReal("TolF", T1_tol);
            opts.setReal("TolU", T1_tol);
            opts.setReal("Scaling", 0.5);
            solver.setOptions(opts);
            solver.initialize();
            gsALMBase<real_t> & b = solver;
            for (index_t k = 0; k < 3; ++k)
                CHECK(b.step() == gsStatus::Success);
            CHECK(!b.hasFrozenScaling());
            CHECK_EQUAL(real_t(0), b.frozenScaling());
        }

        // Presets without a frozen scaling ignore the setter
        {
            T1Setup<gsALMRiks<real_t> >      Rx(0.05, false, false, 1), Ry(0.05, false, false, 1);
            T1Setup<gsALMCrisfield<real_t> > Cx(0.05, false, false, 1), Cy(0.05, false, false, 1);
            std::vector<std::pair<real_t,real_t> > rx, ry, cx, cy;
            gsALMBase<real_t> & bRx = Rx.solver;
            gsALMBase<real_t> & bCx = Cx.solver;
            CHECK(!bRx.hasFrozenScaling());
            CHECK_EQUAL(real_t(0), bRx.frozenScaling());
            CHECK(!bCx.hasFrozenScaling());
            CHECK_EQUAL(real_t(0), bCx.frozenScaling());
            bRx.setFrozenScaling(1.0);
            bCx.setFrozenScaling(1.0);
            frozenScalingRun(Rx, 5, rx);
            frozenScalingRun(Ry, 5, ry);
            frozenScalingRun(Cx, 5, cx);
            frozenScalingRun(Cy, 5, cy);
            CHECK(!rx.empty());
            CHECK(!cx.empty());
            CHECK_EQUAL(ry.size(), rx.size());
            CHECK_EQUAL(cy.size(), cx.size());
            for (size_t k = 0; k < rx.size() && k < ry.size(); ++k)
            {
                CHECK_EQUAL(ry[k].first,  rx[k].first);
                CHECK_EQUAL(ry[k].second, rx[k].second);
            }
            for (size_t k = 0; k < cx.size() && k < cy.size(); ++k)
            {
                CHECK_EQUAL(cy[k].first,  cx[k].first);
                CHECK_EQUAL(cy[k].second, cx[k].second);
            }
            CHECK(!bRx.hasFrozenScaling());
            CHECK(!bCx.hasFrozenScaling());
        }
    }

    // ------------------------------------------------------------------
    // Gates of the presets on the continuation driver.
    // ------------------------------------------------------------------

    /// The truss in the driver convention (R = lambda F - F_int, K = dF_int/dU)
    /// on a caller's problem and state space. Every tangent assembly records
    /// its load factor and the working state it sees.
    struct FacadeTrussRig
    {
        gsStage0Truss<real_t> truss;
        gsStructuralAnalysisOps<real_t>::Jacobian_t   jac;
        gsStructuralAnalysisOps<real_t>::ALResidual_t alres;
        gsVectorStateSpace<real_t> space;
        gsContinuationProblem<real_t> P;
        gsVector<real_t> Ralm;
        std::vector<real_t> jacLambda;
        std::vector<gsVector<real_t> > jacState;

        FacadeTrussRig() : truss(), jac(truss.jacobian()), alres(truss.alResidual()), space(2)
        {
            P.F = truss.force();
            P.residual = [this](real_t l, gsVector<real_t> & R, gsEvalMode) -> bool
            {
                if (!alres(space.working(), l, Ralm))
                    return false;
                R = -Ralm;
                return true;
            };
            P.jacobian = [this](real_t l, gsSparseMatrix<real_t> & K) -> bool
            {
                jacLambda.push_back(l);
                jacState.push_back(space.working());
                return jac(space.working(), K);
            };
        }
    private:
        FacadeTrussRig(const FacadeTrussRig &);
        FacadeTrussRig & operator=(const FacadeTrussRig &);
    };

    /// Truss, a Jacobian operator that counts its calls, and a non-adaptive
    /// solver of type S (arc length \a ds, TolF = TolU = \a tol, the given
    /// BifurcationMethod). \a tweak may change further options before they are applied.
    template<class S>
    struct FacadeSetup
    {
        typedef typename gsStructuralAnalysisOps<real_t>::Jacobian_t   Jacobian_t;
        typedef typename gsStructuralAnalysisOps<real_t>::ALResidual_t ALResidual_t;
        typedef std::function<void(gsOptionList &)> Tweak;

        gsStage0Truss<real_t> truss;
        index_t               count;   ///< Jacobian assemblies so far
        Jacobian_t            J;
        ALResidual_t          R;
        gsVector<real_t>      F;
        S                     solver;

        FacadeSetup(real_t ds, bool quasi, index_t bifMethod, real_t tol = T1_tol, const Tweak & tweak = Tweak())
        : truss(), count(0),
          J([this](gsVector<real_t> const & U, gsSparseMatrix<real_t> & K) -> bool
            { ++count; return truss.jacobianMatrix(U, K); }),
          R(truss.alResidual()), F(truss.force()), solver(J, R, F)
        {
            gsOptionList opts = solver.options();
            opts.setReal("Length", ds);
            opts.setSwitch("AdaptiveLength", false);
            opts.setInt("AdaptiveIterations", 4);
            opts.setSwitch("Quasi", quasi);
            opts.setReal("TolF", tol);
            opts.setReal("TolU", tol);
            opts.setInt("BifurcationMethod", bifMethod);
            if (tweak)
                tweak(opts);
            solver.setOptions(opts);
            solver.initialize();
        }
    private:
        FacadeSetup(const FacadeSetup &);
        FacadeSetup & operator=(const FacadeSetup &);
    };

    /// The row gate of the T1 recordings: status exact, lambda and |U| to
    /// 1e-8 (relative to max(|ref|,1)), iterations exact unless the row is a boundary step
    static void facadeGateRow(const T1Step & ref, const T1Meta & meta,
                              const gsStatus st, const real_t L, const real_t normU, const index_t it)
    {
        CHECK(st == ref.status);
        CHECK_CLOSE(ref.L,     L,     1e-8*math::max(math::abs(ref.L), real_t(1)));
        CHECK_CLOSE(ref.normU, normU, 1e-8*math::max(math::abs(ref.normU), real_t(1)));
        if (meta.gateIterations && !ref.boundary)
            CHECK_EQUAL(ref.iterations, it);
    }

    // ---- distance -------------------------------------------------------

    /// After three accepted steps the distance of the accepted increment is
    /// its arc length: |d^2 - s^2| <= ConstraintTol s^2 by the acceptance gate
    /// (load control: |DL| = s). gsAPALM relies on this to compare the pieces of one path.
    template<class S>
    void facadeDistanceAfterSteps(T1Setup<S> & s, const bool loadControl, real_t & worst)
    {
        gsALMBase<real_t> & b = s.solver;
        for (index_t k = 0; k < 3; ++k)
        {
            CHECK(b.step() == gsStatus::Success);
            const real_t len = b.getLength();
            const real_t d = b.distance(b.solutionDU(), b.solutionDL());
            CHECK((d - d) == 0);
            CHECK(b.solutionDL() != 0);
            CHECK(b.solutionDU().norm() > 0);
            const real_t err = loadControl ? math::abs(d - len)/len : math::abs(d*d - len*len)/(len*len);
            CHECK(err <= 1e-6);
            worst = math::max(worst, err);
        }
    }

    /// distance() is defined before the first step and has the closed form of the
    /// constraint: |DL| (load control), sqrt(phi |U|^2 + (1-phi) DL^2) with phi = 1/n
    /// (Riks), sqrt(|U|^2 + phi^2 |F|^2 DL^2) with the sentinel phi = -1 of the
    /// automatic scaling (Crisfield, consistent Crisfield)
    template<class S>
    void facadeDistance(const T1Preset preset)
    {
        const T1Params p = t1Params(preset, false, false);
        T1Setup<S> s(p.ds, false, false, 1);
        gsALMBase<real_t> & b = s.solver;
        const gsVector<real_t> v = s.F/s.F.norm()*0.1;
        const real_t DL = 0.2, FF = s.F.dot(s.F);
        real_t ref = 0;
        if (preset == LoadControl)
            ref = DL;
        else if (preset == Riks)
            ref = math::sqrt(0.5*v.dot(v) + 0.5*DL*DL);
        else
            ref = math::sqrt(v.dot(v) + FF*DL*DL);
        const real_t d = b.distance(v, DL);
        CHECK((d - d) == 0);
        CHECK(d > 0);
        CHECK_CLOSE(ref, d, 1e-14*ref);
        real_t worst = 0;
        facadeDistanceAfterSteps(s, preset == LoadControl, worst);
        gsInfo << "Facade distance preset " << int(preset) << ": before the first step " << d << " (closed form "
               << ref << "), max |d^2 - s^2|/s^2 over 3 accepted steps " << worst << "\n";
    }

    /// A frozen scaling set before the first step enters the distance
    template<class S>
    void facadeDistanceFrozenScaling()
    {
        T1Setup<S> s(0.05, false, false, 1);
        gsALMBase<real_t> & b = s.solver;
        const gsVector<real_t> v = s.F/s.F.norm()*0.1;
        const real_t DL = 0.2, p = 0.37;
        b.setFrozenScaling(p);
        const real_t ref = math::sqrt(v.dot(v) + p*p*s.F.dot(s.F)*DL*DL);
        const real_t d = b.distance(v, DL);
        CHECK((d - d) == 0);
        CHECK_CLOSE(ref, d, 1e-14*ref);
    }

    /// A user "Scaling" enters the distance before the first step
    template<class S>
    void facadeDistanceUserScaling()
    {
        GateSetup<S> s(0.05, T1_tol, 100, 0.5);
        const gsVector<real_t> v = s.F/s.F.norm()*0.1;
        const real_t DL = 0.2;
        const real_t ref = math::sqrt(v.dot(v) + 0.25*s.F.dot(s.F)*DL*DL);
        const real_t d = s.solver.distance(v, DL);
        CHECK((d - d) == 0);
        CHECK_CLOSE(ref, d, 1e-14*ref);
    }

    TEST(Facade_distance_all_presets)
    {
        UNITTEST_TIME_CONSTRAINT(1000);
        facadeDistance<gsALMLoadControl<real_t> >(LoadControl);
        facadeDistance<gsALMRiks<real_t> >(Riks);
        facadeDistance<gsALMCrisfield<real_t> >(Crisfield);
        facadeDistance<gsALMConsistentCrisfield<real_t> >(ConsistentCrisfield);
        facadeDistanceFrozenScaling<gsALMConsistentCrisfield<real_t> >();
        facadeDistanceUserScaling<gsALMCrisfield<real_t> >();
    }

    // ---- restarting a path piece (APALM) --------------------------------

    /// A fresh solver started with setLength, setSolution, setPrevious (and the
    /// frozen scaling of the run that produced the point) continues a path
    /// exactly like the run that produced it: rows k+1 and k+2 of the recording
    /// \a ref. With \a pastLimit the window lies behind the limit point, where
    /// a solver that lost its history takes the first-step predictor, which
    /// points toward increasing lambda and lands on another point.
    template<class S>
    void facadeApalm(const T1Preset preset, const std::vector<T1Step> & ref, const T1Meta & meta,
                     const index_t k, const bool pastLimit)
    {
        const T1Params p = t1Params(preset, false, false);
        CHECK(ref.size() > static_cast<size_t>(k + 2));
        if (pastLimit)
        {
            CHECK(ref[k].L < ref[k-1].L);
            CHECK(ref[k+1].L < ref[k].L);
            CHECK(ref[k+2].L < ref[k+1].L);
        }

        gsVector<real_t> Uk, Ukm1;
        real_t Lk = 0, Lkm1 = 0, phi = 0;
        bool frozen = false;
        {
            T1Setup<S> A(p.ds, false, false, 1);
            for (index_t i = 0; i <= k; ++i)
            {
                CHECK(A.solver.step() == gsStatus::Success);
                if (i == k - 1)
                {
                    Ukm1 = A.solver.solutionU();
                    Lkm1 = A.solver.solutionL();
                }
            }
            Uk = A.solver.solutionU();
            Lk = A.solver.solutionL();
            CHECK_CLOSE(ref[k].L, Lk, 1e-8*math::max(math::abs(ref[k].L), real_t(1)));
            const gsALMBase<real_t> & a = A.solver;
            frozen = a.hasFrozenScaling();
            phi = a.frozenScaling();
            if (frozen)
                CHECK(phi > 0);
        }
        CHECK_EQUAL(2, Uk.size());
        CHECK_EQUAL(2, Ukm1.size());

        T1Setup<S> B(p.ds, false, false, 1);
        gsALMBase<real_t> & b = B.solver;
        b.setLength(p.ds);
        b.setSolution(Uk, Lk);
        b.setPrevious(Ukm1, Lkm1);
        if (frozen)
            b.setFrozenScaling(phi);
        real_t dL = 0;
        for (index_t j = 1; j <= 2; ++j)
        {
            const gsStatus st = b.step();
            facadeGateRow(ref[k + j], meta, st, b.solutionL(), b.solutionU().norm(), b.numIterations());
            dL = math::max(dL, math::abs(b.solutionL() - ref[k + j].L));
        }
        if (frozen)
            CHECK_EQUAL(phi, b.frozenScaling());

        // Without the frozen scaling the restarted run freezes its own, taken at the
        // restart point, which differs from that of the run it continues
        real_t dLfree = 0;
        if (frozen)
        {
            T1Setup<S> C(p.ds, false, false, 1);
            gsALMBase<real_t> & c = C.solver;
            c.setLength(p.ds);
            c.setSolution(Uk, Lk);
            c.setPrevious(Ukm1, Lkm1);
            CHECK(c.step() == gsStatus::Success);
            CHECK(c.hasFrozenScaling());
            const real_t own = t1TangentPhi(C.truss, C.F, Uk);
            CHECK_CLOSE(own, c.frozenScaling(), 1e-10*own);
            CHECK(c.frozenScaling() != phi);
            dLfree = math::abs(c.solutionL() - ref[k + 1].L);
        }
        gsInfo << "Facade apalm preset " << int(preset) << " window " << k << ": max |lambda - reference| over 2 steps " << dL;
        if (frozen)
            gsInfo << ", first step without the frozen scaling " << dLfree;
        gsInfo << "\n";
    }

    TEST(Facade_apalm_sequence)
    {
        UNITTEST_TIME_CONSTRAINT(1000);
        facadeApalm<gsALMLoadControl<real_t> >(LoadControl, T1_LoadControl_A0_Q0, T1_LoadControl_A0_Q0_meta, 4, false);
        facadeApalm<gsALMRiks<real_t> >(Riks, T1_Riks_A0_Q0, T1_Riks_A0_Q0_meta, 22, true);
        facadeApalm<gsALMCrisfield<real_t> >(Crisfield, T1_Crisfield_A0_Q0, T1_Crisfield_A0_Q0_meta, 22, true);
        facadeApalm<gsALMConsistentCrisfield<real_t> >(ConsistentCrisfield, T1_ConsistentCrisfield_A0_Q0,
                                                       T1_ConsistentCrisfield_A0_Q0_meta, 22, true);
    }

    // ---- one factorization per Jacobian ---------------------------------

    /// A preset that exposes its driver
    template<class P>
    class FacadeProbe : public P
    {
    public:
        template<class A, class B, class V>
        FacadeProbe(A & J, B & R, V & F) : P(J, R, F) { }
        using P::driver;
    };

    /// On the preset path every Jacobian assembly is followed by exactly one
    /// factorization of the driver's tangent solver, also when the stability
    /// is evaluated, in full Newton and in Quasi mode. Only this solver is
    /// counted: a refactorization by another solver (the one behind
    /// computeStability(false) when no committed tangent is held) is not visible here.
    template<class S>
    void facadeFactorizations(const T1Preset preset, const index_t bif, const bool quasi)
    {
        const T1Params p = t1Params(preset, false, quasi);
        FacadeSetup<S> s(p.ds, quasi, bif);
        const index_t jac0 = s.count;
        const index_t fact0 = s.solver.driver().linearSolver().numFactorizations();
        index_t ok = 0;
        for (index_t k = 0; k < p.N; ++k)
        {
            if (s.solver.step() != gsStatus::Success)
                break;
            ++ok;
        }
        const index_t jac  = s.count - jac0;
        const index_t fact = s.solver.driver().linearSolver().numFactorizations() - fact0;
        CHECK(jac > 0);
        CHECK(ok >= 3);
        CHECK_EQUAL(jac, fact);
        gsInfo << "Facade factorizations preset " << int(preset) << (bif == 0 ? " Determinant" : " Eigenvalue")
               << (quasi ? " Quasi" : "") << ": " << ok << " steps, " << jac << " Jacobians, " << fact << " factorizations\n";
    }

    template<template<class> class W>
    void facadeFactorizationsAll()
    {
        for (index_t bif = 0; bif < 2; ++bif)
            for (index_t q = 0; q < 2; ++q)
            {
                facadeFactorizations<W<gsALMLoadControl<real_t> > >(LoadControl, bif, q == 1);
                facadeFactorizations<W<gsALMRiks<real_t> > >(Riks, bif, q == 1);
                facadeFactorizations<W<gsALMCrisfield<real_t> > >(Crisfield, bif, q == 1);
                facadeFactorizations<W<gsALMConsistentCrisfield<real_t> > >(ConsistentCrisfield, bif, q == 1);
            }
    }

    TEST(Facade_one_factorization_per_jacobian)
    {
        UNITTEST_TIME_CONSTRAINT(1000);
        facadeFactorizationsAll<FacadeProbe>();
    }

    // ---- stability at the committed point -------------------------------

    /// After every accepted step the stored indicator and sign equal those of a
    /// fresh evaluation at the committed point, exactly. The stability of the
    /// last Newton iterate's tangent differs from it by ~1e-11 in full Newton
    /// and ~1e-2 in Quasi mode, so a tolerance would miss the first.
    template<class S>
    void facadeStoredFresh(const index_t bif, const bool quasi)
    {
        FacadeSetup<S> s(0.05, quasi, bif);
        gsALMBase<real_t> & base = s.solver;
        index_t ok = 0, nStable = 0, nUnstable = 0;
        real_t worst = 0;
        for (index_t k = 0; k < 40; ++k)
        {
            const gsStatus st = base.step();
            CHECK(st == gsStatus::Success);
            if (st != gsStatus::Success)
                break;
            ++ok;
            const real_t ind = base.indicator();
            const bool stable = base.isStable();
            CHECK(base.computeStability(true) == gsStatus::Success);
            CHECK_EQUAL(ind, base.indicator());
            CHECK_EQUAL(stable, base.isStable());
            worst = math::max(worst, math::abs(ind - base.indicator()));
            (stable ? nStable : nUnstable)++;
        }
        CHECK(ok >= 30);
        CHECK(nStable >= 3);
        CHECK(nUnstable >= 3);
        gsInfo << "Facade stored vs fresh" << (bif == 0 ? " Determinant" : " Eigenvalue") << (quasi ? " Quasi" : "")
               << ": " << ok << " steps, " << nStable << " stable, " << nUnstable << " unstable, max |stored - fresh| "
               << worst << "\n";
    }

    /// The tangent of an explicit stability or bifurcation evaluation is assembled
    /// at the committed point, load factor solutionL() and state solutionU(), not
    /// at the load factor of a trial increment. The legacy operators ignore the
    /// load factor, so this runs on a caller's problem that records it.
    template<class S>
    void facadeTangentLoadFactor()
    {
        FacadeTrussRig rig;
        S solver(rig.P, rig.space);
        gsOptionList opts = solver.options();
        opts.setReal("Length", 0.05);
        opts.setSwitch("AdaptiveLength", false);
        opts.setReal("TolF", T1_tol);
        opts.setReal("TolU", T1_tol);
        solver.setOptions(opts);
        solver.initialize();
        gsALMBase<real_t> & base = solver;

        real_t Lmax = 0;
        index_t evaluations = 0;
        for (index_t k = 0; k < 25; ++k)
        {
            CHECK(base.step() == gsStatus::Success);
            CHECK(base.solutionDL() != 0);
            Lmax = math::max(Lmax, base.solutionL());
            for (index_t pass = 0; pass < 2; ++pass)
            {
                rig.jacLambda.clear();
                rig.jacState.clear();
                if (pass == 0)
                    base.computeStability(true);
                else
                    base.isBifurcation(true);
                CHECK(!rig.jacLambda.empty());
                CHECK_EQUAL(rig.jacLambda.size(), rig.jacState.size());
                for (size_t i = 0; i < rig.jacLambda.size() && i < rig.jacState.size(); ++i)
                {
                    CHECK(rig.jacLambda[i] == base.solutionL());
                    CHECK_EQUAL(base.solutionU().size(), rig.jacState[i].size());
                    CHECK((rig.jacState[i] - base.solutionU()).norm() == 0);
                    ++evaluations;
                }
            }
        }
        CHECK(base.solutionL() < Lmax);
        gsInfo << "Facade tangent load factor: " << evaluations << " tangents at the committed point, lambda "
               << base.solutionL() << " behind the limit load " << Lmax << "\n";
    }

    TEST(Facade_stability_at_committed_point)
    {
        UNITTEST_TIME_CONSTRAINT(1000);
        for (index_t bif = 0; bif < 2; ++bif)
            for (index_t q = 0; q < 2; ++q)
            {
                facadeStoredFresh<gsALMRiks<real_t> >(bif, q == 1);
                facadeStoredFresh<gsALMCrisfield<real_t> >(bif, q == 1);
            }
        facadeTangentLoadFactor<gsALMRiks<real_t> >();
    }

    // ---- stability after switchBranch and setSolution --------------------

    /// After a branch switch the stored stability is unknown: no change is
    /// reported right after the switch, nor by the first step after it, whose
    /// accepted value only seeds the stored one. The premises that make both
    /// checks able to fail are asserted: the step that finds the change leaves a
    /// stored pair (unstable, stable); the point of the switch is stable, so a
    /// previous stability taken from it would report a change; the first step
    /// after the switch ends unstable.
    template<class S>
    void facadeSwitchBranch(const index_t bif)
    {
        FacadeSetup<S> s(0.05, false, bif);
        gsALMBase<real_t> & base = s.solver;

        gsVector<real_t> Uold;
        real_t Lold = 0;
        bool found = false;
        for (index_t k = 0; k < 40 && !found; ++k)
        {
            Uold = base.solutionU();
            Lold = base.solutionL();
            const gsStatus st = base.step();
            CHECK(st == gsStatus::Success);
            if (st != gsStatus::Success)
                break;
            found = base.stabilityChange();
        }
        CHECK(found);
        // the step that found the change ended unstable, after a stable step
        CHECK(base.stabilityChange());
        CHECK(!base.isStable());

        // the branch switch needs the eigenvector estimate that isBifurcation(true) stores;
        // storing it leaves the reported change in place
        base.isBifurcation(true);
        CHECK(base.stabilityChange());

        base.computeSingularPoint(Uold, Lold, true, false, false);
        CHECK(!base.stabilityChange());
        CHECK(!base.isStable());

        // the point of the switch itself is stable: a previous stability taken from it
        // would report a change in the first step after the switch
        {
            FacadeSetup<S> f(0.05, false, bif);
            f.solver.setSolution(base.solutionU(), base.solutionL());
            CHECK(f.solver.isStable());
        }

        CHECK(base.step() == gsStatus::Success);
        // the first step after the switch ends unstable, so a change would be reported
        // if the previous stability were known
        CHECK(!base.isStable());
        CHECK(!base.stabilityChange());
        gsInfo << "Facade switchBranch" << (bif == 0 ? " Determinant" : " Eigenvalue")
               << ": premises hold, no change reported after the switch or in the first step after it\n";
    }

    /// setSolution on an initialized solver evaluates the stability once at
    /// the set point (one Jacobian assembly), stored as the previous stability of
    /// the next step: a job set to just before a stability change reports the
    /// change in its first step, also when the solver was beyond it before.
    template<class S>
    void facadeSetSolution(const index_t bif)
    {
        FacadeSetup<S> A(0.05, false, bif);
        gsALMBase<real_t> & a = A.solver;
        std::vector<gsVector<real_t> > Us;
        std::vector<real_t> Ls;
        Us.push_back(a.solutionU());
        Ls.push_back(a.solutionL());
        index_t j = -1;
        for (index_t k = 0; k < 40 && j < 0; ++k)
        {
            CHECK(a.step() == gsStatus::Success);
            Us.push_back(a.solutionU());
            Ls.push_back(a.solutionL());
            if (a.stabilityChange())
                j = k + 1;
        }
        CHECK(j >= 2);
        if (j < 2)
            return;
        for (index_t k = 0; k < 3; ++k)
        {
            CHECK(a.step() == gsStatus::Success);
            Us.push_back(a.solutionU());
            Ls.push_back(a.solutionL());
        }

        // (a) stored == fresh, with a single assembly
        {
            FacadeSetup<S> B(0.05, false, bif);
            gsALMBase<real_t> & b = B.solver;
            CHECK(b.isStable());
            const index_t c0 = B.count;
            b.setSolution(Us[j + 3], Ls[j + 3]);
            CHECK_EQUAL(1, B.count - c0);
            CHECK(!b.isStable());
            CHECK(!b.stabilityChange());
            const real_t ind = b.indicator();
            b.computeStability(true);
            CHECK_EQUAL(ind, b.indicator());
        }

        // (b) a solver beyond the change, set back to just before it
        CHECK(!a.isStable());
        a.setSolution(Us[j - 1], Ls[j - 1]);
        a.setPrevious(Us[j - 2], Ls[j - 2]);
        CHECK(a.isStable());
        CHECK(a.step() == gsStatus::Success);
        CHECK(a.stabilityChange());
        CHECK(!a.isStable());
        gsInfo << "Facade setSolution" << (bif == 0 ? " Determinant" : " Eigenvalue")
               << ": change found at step " << j << ", reported again by the first step of the restarted job\n";
    }

    TEST(Facade_no_change_after_switchBranch)
    {
        UNITTEST_TIME_CONSTRAINT(1000);
        facadeSwitchBranch<gsALMRiks<real_t> >(0);
        facadeSwitchBranch<gsALMRiks<real_t> >(1);
    }

    TEST(Facade_setSolution_evaluates_stability)
    {
        UNITTEST_TIME_CONSTRAINT(1000);
        facadeSetSolution<gsALMRiks<real_t> >(0);
        facadeSetSolution<gsALMRiks<real_t> >(1);
    }

    // ---- frozen scaling of consistent Crisfield --------------------------

    /// initialize() starts a new run and forgets a frozen scaling, also one set
    /// before it, so the run freezes its own at its first step; the setter
    /// rejects values that are not positive and finite
    template<class S>
    void facadeFrozenScalingReset()
    {
        T1Setup<S> A(0.05, false, false, 1);
        gsALMBase<real_t> & b = A.solver;
        const real_t phi0 = t1TangentPhi(A.truss, A.F, gsVector<real_t>::Zero(2));
        const real_t q = real_t(1.5)*phi0;
        b.setFrozenScaling(q);
        b.applyOptions();
        CHECK(b.hasFrozenScaling());
        CHECK_EQUAL(q, b.frozenScaling());
        b.initialize();
        CHECK(!b.hasFrozenScaling());
        CHECK_EQUAL(real_t(0), b.frozenScaling());
        CHECK(b.step() == gsStatus::Success);
        CHECK_CLOSE(phi0, b.frozenScaling(), 1e-10*phi0);

        T1Setup<S> B(0.05, false, false, 1);
        gsALMBase<real_t> & c = B.solver;
        c.setFrozenScaling(q);
        c.initialize();
        CHECK(!c.hasFrozenScaling());

        CHECK_THROW(c.setFrozenScaling(0), std::runtime_error);
        CHECK_THROW(c.setFrozenScaling(-1), std::runtime_error);
        CHECK_THROW(c.setFrozenScaling(std::numeric_limits<real_t>::quiet_NaN()), std::runtime_error);
        CHECK_THROW(c.setFrozenScaling(std::numeric_limits<real_t>::infinity()), std::runtime_error);
        CHECK(!c.hasFrozenScaling());
    }

    TEST(Facade_CC_frozen_scaling_reset_and_guards)
    {
        UNITTEST_TIME_CONSTRAINT(1000);
        facadeFrozenScalingReset<gsALMConsistentCrisfield<real_t> >();
    }

    // ---- fallback roots of the Crisfield constraint ----------------------

    /// Crisfield preset that exposes the root counters of its constraint
    class FacadeCris : public gsALMCrisfield<real_t>
    {
    public:
        template<class A, class B, class V>
        FacadeCris(A & J, B & R, V & F) : gsALMCrisfield<real_t>(J, R, F) { }
        index_t numModified()   const { return this->m_cris->numModifiedRoots(); }
        index_t numSimplified() const { return this->m_cris->numSimplifiedRoots(); }
        gsCrisfieldConstraint<real_t>::Root lastRoot() const { return this->m_cris->lastRoot(); }
    };

    /// The Crisfield constraint on the driver directly, with the acceptance gate
    /// and the settings the preset maps its options to
    struct FacadeDirectCris
    {
        FacadeTrussRig rig;
        memory::shared_ptr<gsCrisfieldConstraint<real_t> > cris;
        gsContinuationDriver<real_t> drv;

        FacadeDirectCris(real_t ds, real_t scaling, real_t tol)
        : rig(), cris(new gsCrisfieldConstraint<real_t>()), drv(rig.P, rig.space)
        {
            cris->setScaling(scaling);
            cris->setAngleMethod(0);
            cris->setRelaxation(1);
            drv.setConstraint(cris);
            drv.setTest(gsConvergenceTest<real_t>::makeArcLengthGate(tol, tol, 1e-6, 100));
            drv.controller().setLength(ds);
            drv.controller().setAdaptive(false, 4);
            drv.corrector().setQuasi(false, -1);
            drv.corrector().setRelaxation(1);
            drv.linearSolver().setSolver("SimplicialLDLT");
            drv.setup();
        }
    };

    /// |g|/s^2 of the accepted increment of the last step
    static real_t facadeConstraintResidual(gsALMBase<real_t> & b)
    {
        const real_t s = b.getLength();
        const real_t d = b.distance(b.solutionDU(), b.solutionDL());
        return math::abs(d*d - s*s)/(s*s);
    }

    /// Steps the preset and the constraint run directly on the driver in
    /// parallel, at most \a n steps, until the preset has taken the fallback root
    /// \a modified (modified eta) or simplified for the first time; their trajectories,
    /// iteration counts and root counters agree bitwise. Returns the index of that
    /// step (-1 if none), with the preset left on it.
    template<class S>
    index_t facadeCrisfieldRoots(S & solver, FacadeDirectCris & direct, const index_t n, const bool modified)
    {
        index_t first = -1;
        for (index_t k = 0; k < n && first < 0; ++k)
        {
            const index_t before = modified ? solver.numModified() : solver.numSimplified();
            const gsStatus st = solver.step();
            const gsStatus sd = direct.drv.step();
            CHECK(st == gsStatus::Success);
            CHECK(sd == st);
            CHECK_EQUAL(direct.drv.numIterations(), solver.numIterations());
            CHECK_EQUAL(direct.drv.lambda(), solver.solutionL());
            CHECK_EQUAL(solver.solutionU().size(), direct.rig.space.committed().size());
            CHECK((solver.solutionU() - direct.rig.space.committed()).norm() == 0);
            CHECK_EQUAL(direct.cris->numModifiedRoots(), solver.numModified());
            CHECK_EQUAL(direct.cris->numSimplifiedRoots(), solver.numSimplified());
            if ((modified ? solver.numModified() : solver.numSimplified()) > before)
                first = k;
        }
        return first;
    }

    /// Scaling 2 with a long first step makes the discriminant of the corrector
    /// quadratic negative, and the modified-eta double root (Lam and Morley) is
    /// taken; the step still converges onto the arc.
    template<class S>
    void facadeCrisfieldModifiedRoot()
    {
        const real_t ds = 1.2;
        FacadeSetup<S> s(ds, false, 1, T1_tol, [](gsOptionList & o) { o.setReal("Scaling", 2); });
        FacadeDirectCris direct(ds, 2, T1_tol);
        const index_t step = facadeCrisfieldRoots(s.solver, direct, 1, true);
        CHECK_EQUAL(0, step);
        CHECK(s.solver.numModified() > 0);
        CHECK_EQUAL(0, s.solver.numSimplified());
        CHECK(s.solver.numIterations() > 1);
        const real_t g = facadeConstraintResidual(s.solver);
        CHECK(g <= 1e-6);
        gsInfo << "Facade Crisfield modified root: " << s.solver.numModified() << " modified root(s) in "
               << s.solver.numIterations() << " iterations, |g|/s^2 of the end point " << g << "\n";
    }

    /// Automatic scaling, a step of 1.6 and acceptance tolerances of 10: the
    /// step is accepted after the first corrector update, whose negative
    /// discriminant (no modified eta exists) makes the simplified rescaling of
    /// the Newton point onto the arc the last correction, so the end point lies
    /// on the arc to round-off.
    template<class S>
    void facadeCrisfieldSimplifiedRoot()
    {
        const real_t ds = 1.6;
        FacadeSetup<S> s(ds, false, 1, 10);
        FacadeDirectCris direct(ds, -1, 10);
        const index_t step = facadeCrisfieldRoots(s.solver, direct, 8, false);
        CHECK(step >= 0);
        CHECK_EQUAL(1, s.solver.numSimplified());
        CHECK_EQUAL(0, s.solver.numModified());
        CHECK_EQUAL(1, s.solver.numIterations());
        CHECK(s.solver.lastRoot() == gsCrisfieldConstraint<real_t>::Root::Simplified);
        const real_t g = facadeConstraintResidual(s.solver);
        CHECK(g <= 1e-12);
        gsInfo << "Facade Crisfield simplified root: step " << step << ", " << s.solver.numSimplified()
               << " simplified root(s), |g|/s^2 of the end point " << g << "\n";
    }

    TEST(Facade_Crisfield_modified_root)
    {
        UNITTEST_TIME_CONSTRAINT(1000);
        facadeCrisfieldModifiedRoot<FacadeCris>();
    }

    TEST(Facade_Crisfield_simplified_root)
    {
        UNITTEST_TIME_CONSTRAINT(1000);
        facadeCrisfieldSimplifiedRoot<FacadeCris>();
    }
}
