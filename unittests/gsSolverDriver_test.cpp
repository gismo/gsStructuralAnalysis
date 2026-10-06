/** @file gsSolverDriver_test.cpp

    @brief Tests of gsContinuationDriver and the parts it composes (state space,
    path constraints, step controller, Newton corrector, tangent solver), through
    their public interfaces, on the two-dof von Mises truss of gsStage0Truss.h.

    * State-space contract: a counting state space records every state
      operation and every driver event in one string, and a checker validates
      that string against the contract (no copies on the stateful path,
      rollback first in every attempt and in every rejected attempt, commit
      directly before every acceptance). Every test that steps a driver on the
      counting space checks its log with it (the pure-step variant after each
      pure step); the setup and copyability tests assert the absence of state
      operations by character counts instead.
    * Property gates: failure reports and rollback, MaxIter = N gives N
      corrector updates, the path constraint holds at acceptance (recomputed
      from the step history, not read from the gate), the Newton order is
      quadratic and the quasi-Newton order is not, the step-length bound and
      the termination of a never-converging step, one factorization per
      Jacobian, the secant retry policy.
    * Reference trajectories: the 16 trajectories of gsStage0TrussReferences.h
      (load control, Riks, Crisfield, consistent Crisfield, each with and
      without AdaptiveLength and Quasi) are reproduced through the driver.
    * Displacement control passes the load limit point; consistent Crisfield
      crosses lambda = 0 forward.
    * In-suite negative controls show that the central checks can fail:
      checker self-tests, the wrong residual sign, a non-binding constraint
      tolerance, a repeated factorization version, a quasi-Newton run and the
      uninjected recovery step.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
 **/

#include "gismo_unittest.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include <gsStructuralAnalysis/src/gsSolverDriver/gsStateSpace.h>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsSolverDriverTypes.h>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsContinuationProblem.h>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsConvergenceTest.h>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsTangentSolver.h>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsPathConstraint.h>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsPathConstraintCrisfield.h>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsStepController.h>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsNewtonCorrector.h>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsContinuationDriver.h>
#include "gsStage0Truss.h"
#include "gsStage0TrussReferences.h"

SUITE(gsSolverDriver_test)
{
    using namespace gismo::stage0;

    typedef gsPathConstraint<real_t>::Ptr  CPtr;
    typedef gsConvergenceTest<real_t>::Ptr TPtr;
    typedef gsContinuationDriver<real_t>   Driver;

    // ------------------------------------------------------------------
    // Counting state space and contract checker
    // ------------------------------------------------------------------

    /// Forwards to a gsVectorStateSpace and appends one character per call to a
    /// shared event log: t tangentSize, y copyable, m magnitude, a apply,
    /// c commit, r rollback, e exportWorking, i importCommitted,
    /// v inverseApplyExact, x chartJacobian. The event observers of the rigs
    /// below append S/I/A/J for StepStarted/IterationDone/StepAccepted/
    /// StepRejected, so one string records the interleaving of state operations
    /// and driver events. Not copyable. The test side reads the state through
    /// inner(), which is not logged.
    class CountingSpace : public gsStateSpace<real_t>
    {
    public:
        CountingSpace(index_t n, std::string & log)
        : m_in(n), m_log(&log), m_copyable(true), m_failApplyAt(-1), m_failCommit(false), m_applies(0)
        {}

        index_t tangentSize() const { put('t'); return m_in.tangentSize(); }
        bool apply(const gsVector<real_t> & dU, real_t s)
        {
            put('a');
            if (m_applies++ == m_failApplyAt) return false;
            return m_in.apply(dU, s);
        }
        bool commit()   { put('c'); return m_failCommit ? false : m_in.commit(); }
        bool rollback() { put('r'); return m_in.rollback(); }
        bool copyable() const { put('y'); return m_copyable; }
        bool exportWorking(gsVector<real_t> & U) const { put('e'); return m_copyable && m_in.exportWorking(U); }
        bool importCommitted(const gsVector<real_t> & U) { put('i'); return m_copyable && m_in.importCommitted(U); }
        bool magnitude(gsVector<real_t> & x) const { put('m'); return m_in.magnitude(x); }
        bool inverseApplyExact() const { put('v'); return false; }
        bool chartJacobian(const gsVector<real_t> & d, gsSparseMatrix<real_t> & K) const
        { put('x'); return m_in.chartJacobian(d, K); }

        gsVectorStateSpace<real_t> & inner() { return m_in; }
        const gsVectorStateSpace<real_t> & inner() const { return m_in; }

        /// Makes copyable() false and exportWorking()/importCommitted() fail
        void setCopyable(bool c) { m_copyable = c; }
        /// The k-th apply() call from now on (0-based) fails; -1 disarms
        void failApplyAt(index_t k) { m_failApplyAt = k; m_applies = 0; }
        /// Makes commit() fail
        void failCommit(bool f) { m_failCommit = f; }

    private:
        CountingSpace(const CountingSpace &);
        CountingSpace & operator=(const CountingSpace &);
        void put(char ch) const { m_log->push_back(ch); }

        gsVectorStateSpace<real_t> m_in;
        std::string * m_log;
        bool m_copyable;
        index_t m_failApplyAt;
        bool m_failCommit;
        index_t m_applies;
    };

    /// The log without the query characters t, y and m (they do not change the
    /// state, so they do not count as "immediately preceding")
    static std::string withoutQueries(const std::string & log)
    {
        std::string q;
        for (size_t k = 0; k < log.size(); ++k)
            if (log[k] != 't' && log[k] != 'y' && log[k] != 'm')
                q.push_back(log[k]);
        return q;
    }

    static std::string violation(int rule, size_t pos, const std::string & q, const char * what)
    {
        std::ostringstream os;
        os << "clause " << rule << " violated at position " << pos << " of '" << q << "': " << what;
        return os.str();
    }

    /// "" when @a log obeys the state-space contract, else a description of the
    /// first violation. Rules, on the log without queries:
    /// 1. stateful path: no e, i, v, x (the driver never copies a state); pure
    ///    step: exactly one i, before the first S, at most one e, only right
    ///    after the A, no v, x;
    /// 2. every S is immediately preceded by r (every attempt begins with a
    ///    rollback);
    /// 3. every J is immediately preceded by r (every rejected attempt rolls
    ///    back);
    /// 4. every A is immediately preceded by c, and every c is immediately
    ///    followed by A (accepted) or by r J (a failed commit, rolled back);
    /// 5. every a and c lies inside an open attempt;
    /// 6. every S is closed by exactly one A or J before the next S or the end.
    static std::string contractViolation(const std::string & log, bool pureStep)
    {
        const std::string q = withoutQueries(log);
        size_t nImport = 0, nExport = 0, firstS = std::string::npos;
        bool open = false;
        for (size_t k = 0; k < q.size(); ++k)
        {
            const char ch = q[k];
            if (ch == 'v' || ch == 'x')
                return violation(1, k, q, "inverseApplyExact/chartJacobian is never called by the driver");
            if (ch == 'e' || ch == 'i')
            {
                if (!pureStep)
                    return violation(1, k, q, "export/import on the stateful path");
                if (ch == 'i')
                {
                    if (++nImport > 1)
                        return violation(1, k, q, "second import in one pure step");
                    if (firstS != std::string::npos)
                        return violation(1, k, q, "import after the first S");
                }
                else
                {
                    if (++nExport > 1)
                        return violation(1, k, q, "second export in one pure step");
                    if (k == 0 || q[k-1] != 'A')
                        return violation(1, k, q, "export not right after A");
                }
            }
            if (ch == 'S')
            {
                if (firstS == std::string::npos) firstS = k;
                if (k == 0 || q[k-1] != 'r')
                    return violation(2, k, q, "S not preceded by r");
                if (open)
                    return violation(6, k, q, "S while the previous attempt is not closed");
                open = true;
            }
            else if (ch == 'A' || ch == 'J')
            {
                if (!open)
                    return violation(6, k, q, "A/J without an open attempt");
                open = false;
                if (ch == 'A' && (k == 0 || q[k-1] != 'c'))
                    return violation(4, k, q, "A not preceded by c");
                if (ch == 'J' && (k == 0 || q[k-1] != 'r'))
                    return violation(3, k, q, "J not preceded by r");
            }
            else if (ch == 'a')
            {
                if (!open)
                    return violation(5, k, q, "apply outside an attempt");
            }
            else if (ch == 'c')
            {
                if (!open)
                    return violation(5, k, q, "commit outside an attempt");
                const bool accepted = (k + 1 < q.size() && q[k+1] == 'A');
                const bool failed   = (k + 2 < q.size() && q[k+1] == 'r' && q[k+2] == 'J');
                if (!accepted && !failed)
                    return violation(4, k, q, "commit followed by neither A nor r J");
            }
        }
        if (open)
            return violation(6, q.size(), q, "attempt not closed at the end of the log");
        if (pureStep && nImport != 1)
            return violation(1, q.size(), q, "a pure step needs exactly one import");
        return std::string();
    }

    static index_t countChar(const std::string & s, char ch)
    {
        return static_cast<index_t>(std::count(s.begin(), s.end(), ch));
    }

    static char eventChar(gsStepEvent<real_t>::Kind k)
    {
        switch (k)
        {
        case gsStepEvent<real_t>::StepStarted:   return 'S';
        case gsStepEvent<real_t>::IterationDone: return 'I';
        case gsStepEvent<real_t>::StepAccepted:  return 'A';
        case gsStepEvent<real_t>::StepRejected:  return 'J';
        }
        return '?';
    }

    // ------------------------------------------------------------------
    // The truss behind a gsContinuationProblem, with fault injection
    // ------------------------------------------------------------------

    /// The truss of gsStage0Truss.h as a gsContinuationProblem on the working
    /// vector @a work of a state space (set with setWork()).
    ///
    /// The truss residual is R_ALM = F_int - lambda F; the driver residual is
    /// R = lambda F - F_int, hence R = -R_ALM (signR = -1). Call counters,
    /// per-attempt Jacobian counts and fault injection are driven by the step
    /// events (feed every event to onEvent()): the residual closure fails, or
    /// returns a vector of the wrong size, when inject/badSize(step, attempt,
    /// call) holds, where attempt counts the attempts of a step from 1 and
    /// call counts the residual calls of the attempt from 0 (call 0 is at the
    /// predicted point, call k after update k).
    class TrussProblem
    {
    public:
        typedef std::function<bool(index_t, index_t, index_t)> Hook_t;

        TrussProblem()
        : truss(), alres(truss.alResidual()), jac(truss.jacobian()),
          work(NULL), lin(NULL), signR(-1), failJacobianAt(-1), recordOrder(false),
          resCalls(0), jacCalls(0), callInAttempt(0), attemptInStep(0), curStep(-1),
          nStarted(0), nAccepted(0), nRejected(0), factStart(0)
        {
            P.F = truss.force();
            P.residual = [this](real_t l, gsVector<real_t> & R, gsEvalMode) -> bool
            {
                const index_t call = callInAttempt++;
                ++resCalls;
                if (inject && inject(curStep, attemptInStep, call))
                    return false;
                if (badSize && badSize(curStep, attemptInStep, call))
                {
                    R = gsVector<real_t>::Zero(3);
                    return true;
                }
                if (!alres(*work, l, Ralm))
                    return false;
                R = signR*Ralm;
                if (recordOrder)
                {
                    // round-off floor eps * || |K(U)| |U| ||, the RoundoffFloor
                    // formula with c = 1
                    gsSparseMatrix<real_t> K;
                    truss.jacobianMatrix(*work, K);
                    const gsMatrix<real_t> Kd = K.toDense();
                    const gsVector<real_t> kabs = Kd.cwiseAbs() * work->cwiseAbs();
                    curR.push_back(R.norm());
                    curF.push_back(std::numeric_limits<real_t>::epsilon() * kabs.norm());
                }
                return true;
            };
            P.jacobian = [this](real_t, gsSparseMatrix<real_t> & K) -> bool
            {
                const index_t j = jacCalls++;
                if (!attemptJac.empty())
                    ++attemptJac.back();
                if (j == failJacobianAt)
                    return false;
                return jac(*work, K);
            };
        }

        void setWork(const gsVector<real_t> & w) { work = &w; }

        /// Bookkeeping of one step event
        void onEvent(const gsStepEvent<real_t> & ev)
        {
            switch (ev.kind)
            {
            case gsStepEvent<real_t>::StepStarted:
                ++nStarted;
                if (ev.ctx.step != curStep) { curStep = ev.ctx.step; attemptInStep = 0; }
                ++attemptInStep;
                callInAttempt = 0;
                attemptJac.push_back(0);
                startedSizes.push_back(ev.ctx.stepSize);
                factStart = lin ? lin->numFactorizations() : 0;
                curR.clear();
                curF.clear();
                break;
            case gsStepEvent<real_t>::IterationDone:
                break;
            case gsStepEvent<real_t>::StepAccepted:
                ++nAccepted;
                attemptFact.push_back(lin ? lin->numFactorizations() - factStart : 0);
                seqR.push_back(curR);
                seqF.push_back(curF);
                break;
            case gsStepEvent<real_t>::StepRejected:
                ++nRejected;
                attemptFact.push_back(lin ? lin->numFactorizations() - factStart : 0);
                break;
            }
        }

        gsStage0Truss<real_t> truss;
        gsStructuralAnalysisOps<real_t>::ALResidual_t alres;
        gsStructuralAnalysisOps<real_t>::Jacobian_t   jac;
        gsContinuationProblem<real_t> P;
        const gsVector<real_t> * work;
        const gsTangentSolver<real_t> * lin;
        gsVector<real_t> Ralm;
        real_t signR;
        Hook_t inject, badSize;
        index_t failJacobianAt;
        bool recordOrder;
        index_t resCalls, jacCalls, callInAttempt, attemptInStep, curStep;
        index_t nStarted, nAccepted, nRejected;
        index_t factStart;
        std::vector<index_t> attemptJac, attemptFact;   ///< Jacobian assemblies and factorizations per attempt
        std::vector<real_t> startedSizes;               ///< step size of every attempt
        std::vector<real_t> curR, curF;
        std::vector<std::vector<real_t> > seqR, seqF;   ///< |R| and floor per residual call, per accepted attempt

    private:
        TrussProblem(const TrussProblem &);
        TrussProblem & operator=(const TrussProblem &);
    };

    /// One record per step of a run: load factor, |U|, apex deflection w,
    /// status, corrector updates
    struct DRec
    {
        real_t L;
        real_t normU;
        real_t w;
        gsStatus status;
        index_t iterations;
    };

    /// The truss problem, a CountingSpace and a gsContinuationDriver on both,
    /// with an observer that logs S/I/A/J and checks that observers see the
    /// space the driver was given, not a copy.
    class TrussRig : public TrussProblem
    {
    public:
        TrussRig()
        : TrussProblem(), log(), space(2, log), d(P, space)
        {
            setWork(space.inner().working());
            lin = &d.linearSolver();
            d.addObserver([this](const gsStepEvent<real_t> & ev)
            {
                log.push_back(eventChar(ev.kind));
                CHECK(&ev.space == &space);
                onEvent(ev);
            });
        }

        void configure(const CPtr & c, const TPtr & t, real_t ds, bool adaptive = false,
                       bool quasi = false, bool doSetup = true)
        {
            d.setConstraint(c);
            d.setTest(t);
            d.controller().setLength(ds);
            d.controller().setAdaptive(adaptive, 4);
            d.corrector().setQuasi(quasi, -1);
            d.corrector().setRelaxation(1);
            if (doSetup)
                d.setup();
        }

        /// The configuration of the reference recordings: convergence tolerance
        /// T1_tol on residual and update, ConstraintTol 1e-6, MaxIter 100,
        /// AdaptiveIterations 4
        void configureT1(const CPtr & c, real_t ds, bool adaptive, bool quasi, bool doSetup = true)
        {
            configure(c, gsConvergenceTest<real_t>::makeArcLengthGate(T1_tol, T1_tol, 1e-6, 100),
                      ds, adaptive, quasi, doSetup);
        }

        /// Steps until N steps are done or a step fails (that step is recorded)
        std::vector<DRec> run(index_t N)
        {
            std::vector<DRec> out;
            for (index_t k = 0; k < N; ++k)
            {
                const gsStatus st = d.step();
                const gsVector<real_t> & U = space.inner().committed();
                const DRec r = { d.lambda(), U.norm(), U[0], st, d.numIterations() };
                out.push_back(r);
                if (st != gsStatus::Success)
                    break;
            }
            return out;
        }

        std::string log;
        CountingSpace space;
        Driver d;
    };

    static TPtr arcLengthGate(real_t tolF, real_t tolU, real_t tolC, index_t maxIter)
    {
        return gsConvergenceTest<real_t>::makeArcLengthGate(tolF, tolU, tolC, maxIter);
    }

    /// Report check by fields (the note is checked by the caller where it matters)
    static void checkReport(const gsSolverReport & rep, gsStatus status, gsDriverOp op,
                            index_t step, index_t iteration)
    {
        CHECK(rep.status == status);
        CHECK(rep.failedOp == op);
        CHECK_EQUAL(step, rep.step);
        CHECK_EQUAL(iteration, rep.iteration);
    }

    // ------------------------------------------------------------------
    // Predicates and the gate, each self-tested in section 1
    // ------------------------------------------------------------------

    /// True when a successful step beyond the limit point has a lower load
    /// factor than its predecessor
    static bool passedLimit(const std::vector<DRec> & rec, real_t wLimit)
    {
        for (size_t i = 1; i < rec.size(); ++i)
            if (rec[i].status == gsStatus::Success && rec[i].w > wLimit && rec[i].L < rec[i-1].L)
                return true;
        return false;
    }

    /// True when every step is a Success with w strictly increasing from 0 (no
    /// reversal), some step has L < 0 at w > 1, and at least 5 steps follow it
    static bool crossesForward(const std::vector<DRec> & rec)
    {
        real_t wprev = 0;
        index_t crossing = -1;
        for (size_t i = 0; i < rec.size(); ++i)
        {
            if (rec[i].status != gsStatus::Success || !(rec[i].w > wprev))
                return false;
            wprev = rec[i].w;
            if (crossing < 0 && rec[i].L < 0 && rec[i].w > 1)
                crossing = static_cast<index_t>(i);
        }
        return crossing >= 0 && static_cast<index_t>(rec.size()) - 1 - crossing >= 5;
    }

    /// Smallest order estimate p = log(r_{j+1}/r_j) / log(r_j/r_{j-1}) over all
    /// triples of the residual sequences @a r (one sequence per step) whose last
    /// value exceeds 10 times its round-off floor in @a f. A triple that is not
    /// strictly contracting counts as p = 0. For Newton's method r_{k+1} =
    /// C r_k^2 gives p = 2 whatever C; for a chord iteration (frozen K)
    /// r_{k+1} = C e_0 r_k gives p = 1. NaN if no triple is assessed.
    static real_t minNewtonOrder(const std::vector<std::vector<real_t> > & r,
                                 const std::vector<std::vector<real_t> > & f,
                                 index_t * count, std::vector<real_t> * all = NULL)
    {
        real_t pmin = std::numeric_limits<real_t>::quiet_NaN();
        index_t n = 0;
        for (size_t s = 0; s < r.size(); ++s)
            for (size_t j = 1; j + 1 < r[s].size(); ++j)
            {
                const real_t a = r[s][j-1], b = r[s][j], c = r[s][j+1];
                if (!(c > 10*f[s][j+1]))
                    continue;
                const real_t p = (c < b && b < a) ? std::log(c/b)/std::log(b/a) : real_t(0);
                pmin = (n == 0) ? p : std::min(pmin, p);
                ++n;
                if (all) all->push_back(p);
            }
        if (count) *count = n;
        return pmin;
    }

    /// The gate of one recorded trajectory against its reference, without
    /// CHECKs; the first differing row is described in @a firstDiff.
    static bool gateWouldPass(const std::vector<T1Step> & ref, const T1Meta & meta,
                              const std::vector<DRec> & rec, std::string * firstDiff = NULL)
    {
        bool pass = !ref.empty();
        if (meta.fragile) pass = pass && rec.size() >= static_cast<size_t>(meta.gated);
        else              pass = pass && ref.size() == rec.size();
        index_t compared = 0;
        for (size_t i = 0; i < static_cast<size_t>(meta.gated) && i < ref.size() && i < rec.size(); ++i)
        {
            ++compared;
            bool ok = rec[i].status == ref[i].status
                && math::abs(ref[i].L - rec[i].L) <= 1e-8*math::max(math::abs(ref[i].L), real_t(1))
                && math::abs(ref[i].normU - rec[i].normU) <= 1e-8*math::max(math::abs(ref[i].normU), real_t(1));
            if (meta.gateIterations && !ref[i].boundary)
                ok = ok && ref[i].iterations == rec[i].iterations;
            if (!ok && firstDiff && firstDiff->empty())
            {
                std::ostringstream os;
                os << "row " << i << ": ref (L " << ref[i].L << ", normU " << ref[i].normU << ", it "
                   << ref[i].iterations << ") vs rec (L " << rec[i].L << ", normU " << rec[i].normU
                   << ", it " << rec[i].iterations << ")";
                *firstDiff = os.str();
            }
            pass = pass && ok;
        }
        return pass && compared == meta.gated;
    }

    /// Checks a driver recording against the reference of one trajectory:
    /// status and the values L and |U| of the first meta.gated rows to
    /// 1e-8*max(|ref|,1) (the problem is nondimensional, lambda_max = 1), the
    /// iteration count of every non-boundary row when meta.gateIterations, the
    /// recording length unless the trajectory is fragile. Prints one summary line.
    static void driverGate(const std::string & name, const std::vector<T1Step> & ref, const T1Meta & meta,
                           const std::vector<DRec> & rec, index_t jacCalls, index_t factorizations,
                           index_t exports, index_t imports)
    {
        CHECK(!ref.empty());
        if (meta.fragile)
            CHECK(rec.size() >= static_cast<size_t>(meta.gated));
        else
            CHECK_EQUAL(ref.size(), rec.size());

        index_t compared = 0, mmIt = 0, mmItB = 0;
        real_t dL = 0, dU = 0;
        for (size_t i = 0; i < static_cast<size_t>(meta.gated) && i < ref.size() && i < rec.size(); ++i)
        {
            ++compared;
            CHECK(rec[i].status == ref[i].status);
            CHECK_CLOSE(ref[i].L,     rec[i].L,     1e-8*math::max(math::abs(ref[i].L), real_t(1)));
            CHECK_CLOSE(ref[i].normU, rec[i].normU, 1e-8*math::max(math::abs(ref[i].normU), real_t(1)));
            dL = math::max(dL, math::abs(ref[i].L - rec[i].L)/math::max(math::abs(ref[i].L), real_t(1)));
            dU = math::max(dU, math::abs(ref[i].normU - rec[i].normU)/math::max(math::abs(ref[i].normU), real_t(1)));
            if (ref[i].iterations != rec[i].iterations)
                (ref[i].boundary ? mmItB : mmIt) += 1;
            if (meta.gateIterations && !ref[i].boundary)
                CHECK_EQUAL(ref[i].iterations, rec[i].iterations);
        }
        CHECK_EQUAL(meta.gated, compared);
        gsInfo << "Driver_T1 " << name << ": rows compared " << compared << " / gated " << meta.gated
               << " (recorded " << rec.size() << ", ref " << ref.size() << "), max rel diff L " << dL
               << ", normU " << dU << ", iteration mismatches non-boundary " << mmIt
               << ", boundary " << mmItB << "; jacobians " << jacCalls << ", factorizations "
               << factorizations << ", export/import " << exports << "/" << imports << "\n";
    }

    /// One reference cell through the driver on a counting space: the gate
    /// plus the contract checks of the whole run.
    static void runT1Cell(const std::string & name, const CPtr & c,
                          const std::vector<T1Step> & ref, const T1Meta & meta,
                          const T1Preset preset, const bool adaptive, const bool quasi)
    {
        const T1Params p = t1Params(preset, adaptive, quasi);
        TrussRig r;
        r.configureT1(c, p.ds, adaptive, quasi);
        const std::vector<DRec> rec = r.run(p.N);

        driverGate(name, ref, meta, rec, r.jacCalls, r.d.linearSolver().numFactorizations(),
                   countChar(r.log, 'e'), countChar(r.log, 'i'));

        index_t successes = 0;
        for (size_t i = 0; i < rec.size(); ++i)
            successes += (rec[i].status == gsStatus::Success) ? 1 : 0;
        CHECK_EQUAL(std::string(), contractViolation(r.log, false));
        CHECK_EQUAL(0, countChar(r.log, 'e'));
        CHECK_EQUAL(0, countChar(r.log, 'i'));
        CHECK_EQUAL(successes, countChar(r.log, 'A'));
        CHECK_EQUAL(r.jacCalls, r.d.linearSolver().numFactorizations());
    }

    // ------------------------------------------------------------------
    // 1. Self-tests of the checkers (negative controls of the suite)
    // ------------------------------------------------------------------

    TEST(Checker_contract_rejects_bad_logs)
    {
        CHECK_EQUAL(std::string(), contractViolation("trSaaacA", false));
        CHECK_EQUAL(std::string(), contractViolation("rSaarJrSacA", false));
        CHECK_EQUAL(std::string(), contractViolation("rSaacArSaarJ", false));
        CHECK_EQUAL(std::string(), contractViolation("rSacrJ", false));

        const char * bad[] = { "SaacA",     // no rollback before S
                               "rSaaJ",     // reject without rollback
                               "rSaaA",     // accept without commit
                               "rSacaA",    // commit followed by neither A nor r J
                               "rSaacAe",   // export on the stateful path
                               "rSaacAi",   // import on the stateful path
                               "arSacA",    // apply outside an attempt
                               "rSaarSacA", // an attempt that is never closed
                               "rSacAx" };  // chartJacobian
        for (size_t k = 0; k < sizeof(bad)/sizeof(bad[0]); ++k)
            CHECK(!contractViolation(bad[k], false).empty());

        CHECK_EQUAL(std::string(), contractViolation("irSacAe", true));
        CHECK(!contractViolation("rSacAe", true).empty());     // no import
        CHECK(!contractViolation("irSacAee", true).empty());   // two exports
    }

    TEST(Checker_order_estimator)
    {
        const real_t tol = 1e-12;
        index_t n = 0;
        std::vector<std::vector<real_t> > r(1), f(1);
        const real_t quad[] = { 1e-1, 1e-2, 1e-4, 1e-8 };
        r[0].assign(quad, quad + 4);
        f[0].assign(4, real_t(0));
        CHECK_CLOSE(2, minNewtonOrder(r, f, &n), tol);
        CHECK_EQUAL(2, n);

        const real_t lin[] = { 1e-1, 1e-2, 1e-3, 1e-4 };
        r[0].assign(lin, lin + 4);
        CHECK_CLOSE(1, minNewtonOrder(r, f, &n), tol);
        CHECK_EQUAL(2, n);

        // every value below 10 times its floor: nothing is assessed
        f[0].assign(4, real_t(1));
        minNewtonOrder(r, f, &n);
        CHECK_EQUAL(0, n);
    }

    TEST(Checker_passed_limit_predicate)
    {
        const real_t wLim = 0.5;
        const DRec up[]    = { {0.2, 0, 0.1, gsStatus::Success, 1}, {0.6, 0, 0.3, gsStatus::Success, 1},
                               {0.9, 0, 0.6, gsStatus::Success, 1}, {0.8, 0, 0.9, gsStatus::Success, 1} };
        std::vector<DRec> rec(up, up + 4);
        CHECK(passedLimit(rec, wLim));

        const DRec mono[]  = { {0.2, 0, 0.1, gsStatus::Success, 1}, {0.6, 0, 0.3, gsStatus::Success, 1},
                               {0.9, 0, 0.6, gsStatus::Success, 1}, {0.95, 0, 0.9, gsStatus::Success, 1} };
        rec.assign(mono, mono + 4);
        CHECK(!passedLimit(rec, wLim));

        // the load drops, but only at a deflection below the limit point
        const DRec early[] = { {0.2, 0, 0.1, gsStatus::Success, 1}, {0.1, 0, 0.2, gsStatus::Success, 1},
                               {0.6, 0, 0.6, gsStatus::Success, 1}, {0.9, 0, 0.9, gsStatus::Success, 1} };
        rec.assign(early, early + 4);
        CHECK(!passedLimit(rec, wLim));
    }

    TEST(Checker_crossing_predicate)
    {
        std::vector<DRec> rec;
        for (index_t k = 0; k < 9; ++k)
        {
            // w = 0.5, 0.75, 1, 1.25, ...; L < 0 first at w = 1.25 (row 3), 5 rows follow
            const DRec x = { (k == 3 ? real_t(-0.1) : real_t(0.1)), 0, real_t(0.5) + real_t(0.25)*k,
                             gsStatus::Success, 1 };
            rec.push_back(x);
        }
        CHECK(crossesForward(rec));

        std::vector<DRec> rev = rec;            // a reversal at one row
        rev[5].w = rev[4].w - 0.01;
        CHECK(!crossesForward(rev));

        std::vector<DRec> none = rec;           // no row with L < 0 at w > 1
        none[3].L = 0.1;
        CHECK(!crossesForward(none));

        std::vector<DRec> shortTail = rec;      // only 4 rows after the crossing
        shortTail.pop_back();
        CHECK(!crossesForward(shortTail));

        std::vector<DRec> failed = rec;         // a failed step
        failed[6].status = gsStatus::NotConverged;
        CHECK(!crossesForward(failed));
    }

    // ------------------------------------------------------------------
    // 2. State-space contract
    // ------------------------------------------------------------------

    /// Each setup-time requirement, violated on its own, makes setup() throw;
    /// the unmodified configuration does not. A driver without setup() refuses
    /// to step and touches the state space not at all.
    TEST(Contract_setup_ensures)
    {
        const CPtr riks(new gsRiksConstraint<real_t>);
        const TPtr gate = arcLengthGate(T1_tol, T1_tol, 1e-6, 100);
        {
            TrussRig r;
            r.configure(riks, gate, 0.05, false, false, false);
            CHECK(r.d.setup());
        }
        {   // (a) F of the wrong size
            TrussRig r;
            r.configure(riks, gate, 0.05, false, false, false);
            r.P.F = gsVector<real_t>::Zero(3);
            CHECK_THROW(r.d.setup(), std::runtime_error);
        }
        {   // (b) no constraint
            TrussRig r;
            r.d.setTest(gate);
            CHECK_THROW(r.d.setup(), std::runtime_error);
        }
        {   // (c) no test
            TrussRig r;
            r.d.setConstraint(riks);
            CHECK_THROW(r.d.setup(), std::runtime_error);
        }
        {   // (d) a test tree without MaxIter: a leaf and a composite
            TrussRig r;
            r.d.setConstraint(riks);
            r.d.setTest(memory::make_shared(new gsAbsResidualTest<real_t>(1e-10)));
            CHECK_THROW(r.d.setup(), std::runtime_error);

            gsAllTest<real_t> * all = new gsAllTest<real_t>;
            all->add(TPtr(new gsRelResidualTest<real_t>(1e-10, gsResidualReference::LoadNorm)));
            r.d.setTest(TPtr(all));
            CHECK_THROW(r.d.setup(), std::runtime_error);
        }
        {   // (e) the secant predictor needs a load-increment constraint
            TrussRig r;
            r.d.setConstraint(riks);
            r.d.setTest(gate);
            r.d.setSecantPredictor(true);
            CHECK_THROW(r.d.setup(), std::runtime_error);
        }
        {   // (f) no residual
            TrussRig r;
            r.d.setConstraint(riks);
            r.d.setTest(gate);
            r.P.residual = gsContinuationProblem<real_t>::Residual_t();
            CHECK_THROW(r.d.setup(), std::runtime_error);
        }
        {   // no setup(): OtherError, no state operation
            TrussRig r;
            r.configure(riks, gate, 0.05, false, false, false);
            CHECK(r.d.step() == gsStatus::OtherError);
            CHECK(r.d.report().failedOp == gsDriverOp::None);
            CHECK(r.log.find_first_of("acrei") == std::string::npos);
        }
    }

    /// A pure step needs a copyable space: without export/import it fails with
    /// {OtherError, Import} before any state operation.
    TEST(Contract_pure_step_needs_copyable)
    {
        const gsVector<real_t> zero = gsVector<real_t>::Zero(2);
        gsPathPoint<real_t> from, to;
        from.U = zero; from.lambda = 0;
        gsStepHistory<real_t> hist, histOut;
        gsSolverReport rep;
        {
            gsVectorStateSpace<real_t> inner(2);
            index_t applies = 0;
            gsFunctionalStateSpace<real_t> sp = gsStateSpace<real_t>::fromOps(2,
                [&inner, &applies](const gsVector<real_t> & dU, real_t s) { ++applies; return inner.apply(dU, s); },
                [&inner]() { return inner.commit(); },
                [&inner]() { return inner.rollback(); });
            CHECK(!sp.copyable());
            gsVector<real_t> U, m;
            CHECK(!sp.exportWorking(U));
            CHECK(sp.magnitude(m));
            CHECK_EQUAL(0, m.size());

            TrussProblem prob;
            prob.setWork(inner.working());
            Driver d(prob.P, sp);
            d.setConstraint(CPtr(new gsRiksConstraint<real_t>));
            d.setTest(arcLengthGate(T1_tol, T1_tol, 1e-6, 100));
            d.controller().setLength(0.05);
            CHECK(d.setup());

            CHECK(d.step(from, hist, 0.05, to, histOut, &rep) == gsStatus::OtherError);
            CHECK(rep.status == gsStatus::OtherError);
            CHECK(rep.failedOp == gsDriverOp::Import);
            CHECK_EQUAL(0, applies);
            CHECK(inner.committed() == zero);
            CHECK(inner.working() == zero);
        }
        {
            TrussRig r;
            r.configureT1(CPtr(new gsRiksConstraint<real_t>), 0.05, false, false);
            r.space.setCopyable(false);
            r.log.clear();
            CHECK(r.d.step(from, hist, 0.05, to, histOut, &rep) == gsStatus::OtherError);
            CHECK(rep.status == gsStatus::OtherError);
            CHECK(rep.failedOp == gsDriverOp::Import);
            CHECK_EQUAL(0, countChar(r.log, 'a'));
            CHECK(r.space.inner().committed() == zero);
        }
    }

    /// One successful pure step from {0, 0} with no history: one import before
    /// the attempt, one export right after it, nothing else copies the state.
    TEST(Contract_pure_step_log)
    {
        TrussRig r;
        r.configureT1(CPtr(new gsRiksConstraint<real_t>), 0.05, false, false);
        gsPathPoint<real_t> from, to;
        from.U = gsVector<real_t>::Zero(2); from.lambda = 0;
        gsStepHistory<real_t> hist, histOut;
        r.log.clear();
        CHECK(r.d.step(from, hist, 0.05, to, histOut) == gsStatus::Success);
        CHECK_EQUAL(std::string(), contractViolation(r.log, true));
        const std::string q = withoutQueries(r.log);
        CHECK(!q.empty() && q[0] == 'i');
        CHECK(q.size() >= 3 && q.substr(q.size() - 3) == "cAe");
    }

    // ------------------------------------------------------------------
    // 3. Status model
    // ------------------------------------------------------------------

    /// The residual fails at call 3 of the attempt of step 2 (after update 3):
    /// the report names the operation, the step and the update, and the full
    /// rollback leaves no trace: committed state, working state, exported state
    /// and lambda are bitwise those after step 1, and the next step reproduces
    /// the uninjected run.
    TEST(Status_assembly_error_step2_iter3)
    {
        UNITTEST_TIME_CONSTRAINT(1000);
        CHECK(T1_Riks_A0_Q0[2].iterations >= 4);   // the attempt of step 2 reaches update 3

        real_t lambda2 = 0;
        {
            TrussRig ref;
            ref.configureT1(CPtr(new gsRiksConstraint<real_t>), 0.05, false, false);
            for (index_t k = 0; k < 3; ++k)
                CHECK(ref.d.step() == gsStatus::Success);
            lambda2 = ref.d.lambda();
        }

        TrussRig r;
        r.configureT1(CPtr(new gsRiksConstraint<real_t>), 0.05, false, false);
        r.inject = [](index_t step, index_t, index_t call) { return step == 2 && call == 3; };
        CHECK(r.d.step() == gsStatus::Success);
        CHECK(r.d.step() == gsStatus::Success);
        const gsVector<real_t> U1 = r.space.inner().committed();
        const real_t L1 = r.d.lambda();
        gsVector<real_t> Uexp1;
        r.space.inner().exportWorking(Uexp1);

        CHECK(r.d.step() == gsStatus::AssemblyError);
        checkReport(r.d.report(), gsStatus::AssemblyError, gsDriverOp::Residual, 2, 3);
        gsVector<real_t> Uexp;
        r.space.inner().exportWorking(Uexp);
        CHECK(Uexp == Uexp1);
        CHECK(r.space.inner().committed() == U1);
        CHECK(r.space.inner().working() == r.space.inner().committed());
        CHECK(r.d.lambda() == L1);
        CHECK_EQUAL(2, r.d.context().step);
        const std::string q = withoutQueries(r.log);
        CHECK(q.size() >= 2 && q.substr(q.size() - 2) == "rJ");

        // control: without the injection the same step succeeds and equals the uninjected run
        r.inject = TrussProblem::Hook_t();
        CHECK(r.d.step() == gsStatus::Success);
        CHECK(r.d.lambda() == lambda2);
        CHECK_EQUAL(std::string(), contractViolation(r.log, false));
    }

    /// One row of the failure-to-report mapping of the driver per sub-case, on
    /// a fresh Riks driver at step 0; the committed state and lambda are
    /// bitwise unchanged after every failed step.
    TEST(Status_mapping_table)
    {
        UNITTEST_TIME_CONSTRAINT(1000);
        const gsVector<real_t> zero = gsVector<real_t>::Zero(2);
        {   // jacobian false at the predictor
            TrussRig r;
            r.configureT1(CPtr(new gsRiksConstraint<real_t>), 0.05, false, false);
            r.failJacobianAt = 0;
            CHECK(r.d.step() == gsStatus::AssemblyError);
            checkReport(r.d.report(), gsStatus::AssemblyError, gsDriverOp::Jacobian, 0, 0);
            CHECK(r.space.inner().committed() == zero);
            CHECK(r.d.lambda() == 0);
            CHECK_EQUAL(std::string(), contractViolation(r.log, false));
        }
        {   // apply false at the predictor: not retried although the controller would allow it
            TrussRig r;
            r.configureT1(CPtr(new gsRiksConstraint<real_t>), 0.05, false, false);
            r.d.controller().setMaxBisect(1);
            r.space.failApplyAt(0);
            CHECK(r.d.step() == gsStatus::OtherError);
            checkReport(r.d.report(), gsStatus::OtherError, gsDriverOp::Apply, 0, 0);
            CHECK_EQUAL(1, countChar(r.log, 'S'));
            CHECK_EQUAL(0, r.d.controller().depth());
            CHECK(r.space.inner().committed() == zero);
            CHECK(r.d.lambda() == 0);
            CHECK_EQUAL(std::string(), contractViolation(r.log, false));
        }
        {   // commit false: the step converged, so k = the recorded iteration count
            TrussRig r;
            r.configureT1(CPtr(new gsRiksConstraint<real_t>), 0.05, false, false);
            r.space.failCommit(true);
            CHECK(r.d.step() == gsStatus::OtherError);
            CHECK_EQUAL(T1_Riks_A0_Q0[0].iterations, r.d.numIterations());
            checkReport(r.d.report(), gsStatus::OtherError, gsDriverOp::Commit, 0, T1_Riks_A0_Q0[0].iterations);
            CHECK(r.space.inner().committed() == zero);
            CHECK(r.d.lambda() == 0);
            CHECK_EQUAL(std::string(), contractViolation(r.log, false));
        }
        {   // post-iteration hook aborts after update 1
            TrussRig r;
            r.configureT1(CPtr(new gsRiksConstraint<real_t>), 0.05, false, false);
            r.d.setPostIteration([](const gsStepContext<real_t> & ctx)
            { return ctx.iteration == 1 ? gsHookDecision::Abort : gsHookDecision::Continue; });
            CHECK(r.d.step() == gsStatus::OtherError);
            checkReport(r.d.report(), gsStatus::OtherError, gsDriverOp::Hook, 0, 1);
            CHECK(r.space.inner().committed() == zero);
            CHECK(r.d.lambda() == 0);
            CHECK_EQUAL(std::string(), contractViolation(r.log, false));
        }
        {   // pre-step hook rejects; the default controller grants no retry
            TrussRig r;
            r.configureT1(CPtr(new gsRiksConstraint<real_t>), 0.05, false, false);
            r.d.setPreStep([](real_t, real_t &) { return gsHookDecision::Reject; });
            CHECK(r.d.step() == gsStatus::NotConverged);
            checkReport(r.d.report(), gsStatus::NotConverged, gsDriverOp::Hook, 0, 0);
            CHECK(r.space.inner().committed() == zero);
            CHECK(r.d.lambda() == 0);
            CHECK_EQUAL(std::string(), contractViolation(r.log, false));
        }
        {   // a residual of the wrong size at the predicted point
            TrussRig r;
            r.configureT1(CPtr(new gsRiksConstraint<real_t>), 0.05, false, false);
            r.badSize = [](index_t, index_t, index_t call) { return call == 0; };
            CHECK(r.d.step() == gsStatus::AssemblyError);
            checkReport(r.d.report(), gsStatus::AssemblyError, gsDriverOp::Residual, 0, 0);
            CHECK(r.space.inner().committed() == zero);
            CHECK(r.d.lambda() == 0);
            CHECK_EQUAL(std::string(), contractViolation(r.log, false));
        }
    }

    // ------------------------------------------------------------------
    // 4. Property gates
    // ------------------------------------------------------------------

    /// With TolF = TolU = 0 the strict residue test never holds, so a step runs
    /// exactly MaxIter = N corrector updates and ends NotConverged with the
    /// iteration-limit report. The predictor (all four constraints need the
    /// tangent at the first step) and each update assemble one K and evaluate
    /// one residual, so both counts are N + 1.
    TEST(MaxIter_is_iteration_count)
    {
        UNITTEST_TIME_CONSTRAINT(1000);
        const index_t budgets[] = { 1, 5, 0 };
        for (int which = 0; which < 4; ++which)
            for (size_t b = 0; b < 3; ++b)
            {
                const index_t N = budgets[b];
                CPtr c;
                real_t ds = 0.05;
                if (which == 0)      { c = CPtr(new gsLoadConstraint<real_t>); ds = 0.07; }
                else if (which == 1) c = CPtr(new gsRiksConstraint<real_t>);
                else if (which == 2) c = CPtr(new gsCrisfieldConstraint<real_t>);
                else                 c = CPtr(new gsConsistentCrisfieldConstraint<real_t>);
                TrussRig r;
                r.configure(c, arcLengthGate(0, 0, 1e-6, N), ds);
                CHECK(r.d.step() == gsStatus::NotConverged);
                CHECK_EQUAL(N, r.d.numIterations());
                checkReport(r.d.report(), gsStatus::NotConverged, gsDriverOp::None, 0, N);
                CHECK(r.d.report().note.find("iteration limit") != std::string::npos);
                CHECK_EQUAL(N + 1, r.resCalls);
                if (N > 0)
                    CHECK_EQUAL(N + 1, r.jacCalls);
                CHECK_EQUAL(std::string(), contractViolation(r.log, false));
            }
    }

    /// The path constraint g holds at every accepted step, recomputed from the
    /// step history with closed forms (never read from the gate, whose
    /// PathConstraint test compares the same value). With TolF = TolU = 10
    /// every finite iterate passes the residue tests, and MinIter(1) still
    /// forces one update, so the constraint tolerance is the binding criterion.
    TEST(Constraint_holds_at_acceptance)
    {
        UNITTEST_TIME_CONSTRAINT(1000);
        enum Kind { KRiks, KCrisfield, KConsistent, KLoad, KDisplacement };
        const Kind kinds[] = { KRiks, KCrisfield, KConsistent, KLoad, KDisplacement };
        const char * names[] = { "Riks", "Crisfield", "ConsistentCrisfield", "Load", "Displacement" };
        for (size_t kk = 0; kk < 5; ++kk)
        {
            CPtr c;
            real_t ds = 0.05;
            if (kinds[kk] == KRiks)             c = CPtr(new gsRiksConstraint<real_t>);
            else if (kinds[kk] == KCrisfield)   { gsCrisfieldConstraint<real_t> * p = new gsCrisfieldConstraint<real_t>; p->setScaling(1); c = CPtr(p); }
            else if (kinds[kk] == KConsistent)  { gsConsistentCrisfieldConstraint<real_t> * p = new gsConsistentCrisfieldConstraint<real_t>; p->setScaling(1); c = CPtr(p); }
            else if (kinds[kk] == KLoad)        { c = CPtr(new gsLoadConstraint<real_t>); ds = 0.07; }
            else                                c = CPtr(new gsDisplacementConstraint<real_t>(1));
            TrussRig r;
            r.configure(c, arcLengthGate(10, 10, 1e-6, 100), ds);
            for (index_t k = 0; k < 3; ++k)
            {
                CHECK(r.d.step() == gsStatus::Success);
                const gsStepHistory<real_t> & h = r.d.history();
                const real_t s = h.stepSize;
                real_t g2 = 0;   // |g|/s^2
                switch (kinds[kk])
                {
                case KRiks:        // phi = 1/n = 1/2 on the truss
                    g2 = math::abs(0.5*h.DeltaU.dot(h.DeltaU) + 0.5*h.DeltaL*h.DeltaL - s*s)/(s*s);
                    break;
                case KCrisfield:
                case KConsistent:  // phi = 1, |F| = 1
                    g2 = math::abs(h.DeltaU.dot(h.DeltaU) + h.DeltaL*h.DeltaL - s*s)/(s*s);
                    break;
                case KLoad:
                    CHECK(h.DeltaL == h.stepSize);
                    CHECK(h.stepSize == real_t(0.07));
                    g2 = math::abs(h.DeltaL - h.stepSize)/(s*s);
                    break;
                case KDisplacement:
                    g2 = math::abs(h.DeltaU[1] - s)/(s*s);
                    break;
                }
                gsInfo << "Constraint_holds " << names[kk] << " step " << k << ": iterations "
                       << r.d.numIterations() << ", |g|/s^2 = " << g2 << "\n";
                CHECK(g2 < 1e-6);
            }
            CHECK_EQUAL(std::string(), contractViolation(r.log, false));
        }

        // control: with ConstraintTol 1e300 the gate does not bind. MinIter(1)
        // forces exactly one update, which leaves the quadratic remainder
        // phi |dU|^2 + (1-phi) dl^2 of a single linearised Riks update in g.
        TrussRig r;
        r.configure(CPtr(new gsRiksConstraint<real_t>), arcLengthGate(10, 10, 1e300, 100), 0.4);
        CHECK(r.d.step() == gsStatus::Success);
        CHECK_EQUAL(1, r.d.numIterations());
        const gsStepHistory<real_t> & h = r.d.history();
        const real_t g2 = math::abs(0.5*h.DeltaU.dot(h.DeltaU) + 0.5*h.DeltaL*h.DeltaL - h.stepSize*h.stepSize)
                          /(h.stepSize*h.stepSize);
        gsInfo << "Constraint_holds control (ConstraintTol 1e300, ds 0.4): |g|/s^2 = " << g2 << "\n";
        CHECK(g2 > 1e-6);
        CHECK_EQUAL(std::string(), contractViolation(r.log, false));
    }

    /// Runs @a steps steps of the Load or Riks driver with the order recording
    /// on and returns the smallest order estimate over the accepted attempts.
    static real_t newtonOrder(bool riks, bool quasi, index_t steps, index_t * count)
    {
        TrussRig r;
        r.recordOrder = true;
        r.configureT1(riks ? CPtr(new gsRiksConstraint<real_t>) : CPtr(new gsLoadConstraint<real_t>),
                      riks ? 0.05 : 0.07, false, quasi);
        for (index_t k = 0; k < steps; ++k)
            CHECK(r.d.step() == gsStatus::Success);
        CHECK_EQUAL(std::string(), contractViolation(r.log, false));
        std::vector<real_t> all;
        const real_t p = minNewtonOrder(r.seqR, r.seqF, count, &all);
        gsInfo << "Newton order " << (riks ? "Riks" : "Load") << (quasi ? " Q1" : " Q0") << ": min p " << p
               << " over " << *count << " triples; p =";
        for (size_t i = 0; i < all.size(); ++i)
            gsInfo << " " << all[i];
        gsInfo << "\n";
        return p;
    }

    /// Newton's method converges quadratically: the order estimate over the
    /// residual sequences of the accepted attempts is at least 1.8 (below the
    /// limit point, steps 0-9 of load control and of Riks).
    TEST(Newton_order_full)
    {
        UNITTEST_TIME_CONSTRAINT(1000);
        index_t n = 0;
        const real_t pLoad = newtonOrder(false, false, 10, &n);
        CHECK(n >= 8);
        CHECK(pLoad >= 1.8);
        const real_t pRiks = newtonOrder(true, false, 10, &n);
        CHECK(n >= 10);
        CHECK(pRiks >= 1.8);
    }

    /// Control: with K refactorized only at update 1 (quasi-Newton, a chord
    /// iteration) the order estimate is about 1, which the estimator detects.
    TEST(Newton_order_quasi_is_linear)
    {
        UNITTEST_TIME_CONSTRAINT(1000);
        index_t n = 0;
        const real_t p = newtonOrder(false, true, 10, &n);
        CHECK(n >= 5);
        CHECK(p < 1.8);
    }

    /// A problem that never converges: R = (1,1) whatever U and lambda, K = I.
    /// The relative residual is at least sqrt(2)/|lambda + Delta lambda|, never
    /// below the tolerance, so every attempt is exhausted after 2 updates.
    static void makeNeverConverging(TrussRig & r)
    {
        r.P.residual = [](real_t, gsVector<real_t> & R, gsEvalMode) -> bool
        { R.resize(2); R << 1, 1; return true; };
        r.P.jacobian = [](real_t, gsSparseMatrix<real_t> & K) -> bool
        {
            K.resize(2, 2);
            K.setIdentity();
            return true;
        };
    }

    /// Arc length: with MaxBisect = 3 a step makes exactly 4 attempts, at s0,
    /// s0/2, s0/4 and s0/8 (halving is exact in binary), so the shortest
    /// attempted length is s0/2^MaxBisect. The next step makes one attempt, at
    /// the length the controller was left with.
    TEST(Controller_bisection_bound)
    {
        UNITTEST_TIME_CONSTRAINT(1000);
        const real_t s0 = 0.05;
        TrussRig r;
        makeNeverConverging(r);
        r.configure(CPtr(new gsRiksConstraint<real_t>), arcLengthGate(T1_tol, T1_tol, 1e-6, 2), s0);
        r.d.controller().setMaxBisect(3);

        CHECK(r.d.step() == gsStatus::NotConverged);
        CHECK_EQUAL(size_t(4), r.startedSizes.size());
        CHECK(r.startedSizes[0] == s0);
        CHECK(r.startedSizes[1] == s0/2);
        CHECK(r.startedSizes[2] == s0/4);
        CHECK(r.startedSizes[3] == s0/8);
        CHECK(*std::min_element(r.startedSizes.begin(), r.startedSizes.end()) == s0/8);
        CHECK(r.d.report().status == gsStatus::NotConverged);
        CHECK(r.d.report().failedOp == gsDriverOp::None);
        CHECK(r.d.report().note.find("; no retry (controller)") != std::string::npos);
        CHECK_EQUAL(3, r.d.controller().depth());
        CHECK(r.d.controller().length() == s0/8);

        CHECK(r.d.step() == gsStatus::NotConverged);
        CHECK_EQUAL(size_t(5), r.startedSizes.size());
        CHECK(r.startedSizes[4] == s0/8);
        CHECK_EQUAL(5, countChar(r.log, 'J'));
        CHECK_EQUAL(0, countChar(r.log, 'A'));
        CHECK_EQUAL(std::string(), contractViolation(r.log, false));
    }

    /// Load control at lambda = 1: the attempts halve s = 0.07 = 1.12 * 2^-4
    /// (exact) until 1 + s == 1 in floating point, which holds iff s <= 2^-digits.
    /// The attempt at depth k has s = 1.12 * 2^-(4+k), so depth digits-4 still
    /// moves lambda and depth digits-3 does not: exactly digits-3 attempts
    /// (50 for double), and then the load-stall guard returns without a
    /// further attempt.
    TEST(Controller_terminates_on_load_stall)
    {
        UNITTEST_TIME_CONSTRAINT(1000);
        const index_t attempts = std::numeric_limits<real_t>::digits - 3;
        TrussRig r;
        makeNeverConverging(r);
        r.configure(CPtr(new gsLoadConstraint<real_t>), arcLengthGate(T1_tol, T1_tol, 1e-6, 2), 0.07);
        r.d.setLambda(1);
        r.d.controller().setMaxBisect(1000);

        CHECK(r.d.step() == gsStatus::NotConverged);
        CHECK_EQUAL(size_t(attempts), r.startedSizes.size());
        CHECK_EQUAL(attempts, r.d.controller().depth());
        checkReport(r.d.report(), gsStatus::NotConverged, gsDriverOp::Controller, 0, -1);
        CHECK_EQUAL(std::string("load step does not change lambda"), r.d.report().note);
        CHECK_EQUAL(std::string(), contractViolation(r.log, false));
    }

    // ------------------------------------------------------------------
    // 5. Reference trajectories through the driver
    // ------------------------------------------------------------------

    TEST(Driver_T1_LoadControl)
    {
        UNITTEST_TIME_CONSTRAINT(1000);
        runT1Cell("T1_LoadControl_A0_Q0", CPtr(new gsLoadConstraint<real_t>), T1_LoadControl_A0_Q0, T1_LoadControl_A0_Q0_meta, LoadControl, false, false);
        runT1Cell("T1_LoadControl_A0_Q1", CPtr(new gsLoadConstraint<real_t>), T1_LoadControl_A0_Q1, T1_LoadControl_A0_Q1_meta, LoadControl, false, true);
        runT1Cell("T1_LoadControl_A1_Q0", CPtr(new gsLoadConstraint<real_t>), T1_LoadControl_A1_Q0, T1_LoadControl_A1_Q0_meta, LoadControl, true,  false);
        runT1Cell("T1_LoadControl_A1_Q1", CPtr(new gsLoadConstraint<real_t>), T1_LoadControl_A1_Q1, T1_LoadControl_A1_Q1_meta, LoadControl, true,  true);
    }

    TEST(Driver_T1_Riks)
    {
        UNITTEST_TIME_CONSTRAINT(1000);
        runT1Cell("T1_Riks_A0_Q0", CPtr(new gsRiksConstraint<real_t>), T1_Riks_A0_Q0, T1_Riks_A0_Q0_meta, Riks, false, false);
        runT1Cell("T1_Riks_A0_Q1", CPtr(new gsRiksConstraint<real_t>), T1_Riks_A0_Q1, T1_Riks_A0_Q1_meta, Riks, false, true);
        runT1Cell("T1_Riks_A1_Q0", CPtr(new gsRiksConstraint<real_t>), T1_Riks_A1_Q0, T1_Riks_A1_Q0_meta, Riks, true,  false);
        runT1Cell("T1_Riks_A1_Q1", CPtr(new gsRiksConstraint<real_t>), T1_Riks_A1_Q1, T1_Riks_A1_Q1_meta, Riks, true,  true);
    }

    /// The Crisfield cells also show that the simplified (complex-root)
    /// fallback is never taken on this path.
    TEST(Driver_T1_Crisfield)
    {
        UNITTEST_TIME_CONSTRAINT(1000);
        const bool adaptive[] = { false, false, true, true };
        const bool quasi[]    = { false, true,  false, true };
        const std::vector<T1Step> * refs[] = { &T1_Crisfield_A0_Q0, &T1_Crisfield_A0_Q1, &T1_Crisfield_A1_Q0, &T1_Crisfield_A1_Q1 };
        const T1Meta * metas[] = { &T1_Crisfield_A0_Q0_meta, &T1_Crisfield_A0_Q1_meta, &T1_Crisfield_A1_Q0_meta, &T1_Crisfield_A1_Q1_meta };
        const char * names[] = { "T1_Crisfield_A0_Q0", "T1_Crisfield_A0_Q1", "T1_Crisfield_A1_Q0", "T1_Crisfield_A1_Q1" };
        for (size_t k = 0; k < 4; ++k)
        {
            memory::shared_ptr<gsCrisfieldConstraint<real_t> > c(new gsCrisfieldConstraint<real_t>);
            runT1Cell(names[k], c, *refs[k], *metas[k], Crisfield, adaptive[k], quasi[k]);
            CHECK_EQUAL(0, c->numSimplifiedRoots());
            gsInfo << "Driver_T1 " << names[k] << ": modified roots " << c->numModifiedRoots()
                   << ", relaxed " << c->numRelaxed() << "\n";
        }
    }

    TEST(Driver_T1_ConsistentCrisfield)
    {
        UNITTEST_TIME_CONSTRAINT(1000);
        runT1Cell("T1_ConsistentCrisfield_A0_Q0", CPtr(new gsConsistentCrisfieldConstraint<real_t>), T1_ConsistentCrisfield_A0_Q0, T1_ConsistentCrisfield_A0_Q0_meta, ConsistentCrisfield, false, false);
        runT1Cell("T1_ConsistentCrisfield_A0_Q1", CPtr(new gsConsistentCrisfieldConstraint<real_t>), T1_ConsistentCrisfield_A0_Q1, T1_ConsistentCrisfield_A0_Q1_meta, ConsistentCrisfield, false, true);
        runT1Cell("T1_ConsistentCrisfield_A1_Q0", CPtr(new gsConsistentCrisfieldConstraint<real_t>), T1_ConsistentCrisfield_A1_Q0, T1_ConsistentCrisfield_A1_Q0_meta, ConsistentCrisfield, true,  false);
        runT1Cell("T1_ConsistentCrisfield_A1_Q1", CPtr(new gsConsistentCrisfieldConstraint<real_t>), T1_ConsistentCrisfield_A1_Q1, T1_ConsistentCrisfield_A1_Q1_meta, ConsistentCrisfield, true,  true);
    }

    /// Control of the residual sign: the driver solves K dU = R with
    /// R = lambda F - F_int. Handing it the ALM residual F_int - lambda F makes
    /// every Newton update point away from the solution, and the gate rejects it.
    TEST(Driver_T1_sign_control)
    {
        UNITTEST_TIME_CONSTRAINT(1000);
        const T1Params p = t1Params(Riks, false, false);
        TrussRig r;
        r.signR = 1;
        r.configureT1(CPtr(new gsRiksConstraint<real_t>), p.ds, false, false);
        const std::vector<DRec> rec = r.run(p.N);
        std::string firstDiff;
        const bool wouldPass = gateWouldPass(T1_Riks_A0_Q0, T1_Riks_A0_Q0_meta, rec, &firstDiff);
        CHECK(!wouldPass);
        CHECK_EQUAL(std::string(), contractViolation(r.log, false));
        gsInfo << "Driver_T1 sign control: ALM sign would pass the gate: " << (wouldPass ? "YES" : "no")
               << "; recorded " << rec.size() << " rows; first differing " << firstDiff << "\n";
    }

    // ------------------------------------------------------------------
    // 6. Pure step, factorization versions
    // ------------------------------------------------------------------

    /// The pure step chained on its own space equals the stateful driver bitwise
    TEST(PureStep_equals_stateful)
    {
        UNITTEST_TIME_CONSTRAINT(1000);
        TrussRig r1, r2;
        r1.configureT1(CPtr(new gsRiksConstraint<real_t>), 0.05, false, false);
        r2.configureT1(CPtr(new gsRiksConstraint<real_t>), 0.05, false, false);

        gsPathPoint<real_t> from, to;
        from.U = gsVector<real_t>::Zero(2); from.lambda = 0;
        gsStepHistory<real_t> hist, histOut;
        index_t steps = 0, bitwise = 0;
        for (index_t k = 0; k < 40; ++k)
        {
            const gsStatus s1 = r1.d.step();
            gsSolverReport rep;
            r2.log.clear();
            const gsStatus s2 = r2.d.step(from, hist, 0.05, to, histOut, &rep);
            CHECK(s1 == s2);
            CHECK(rep.status == r1.d.report().status);
            CHECK_EQUAL(std::string(), contractViolation(r2.log, true));
            if (s1 != gsStatus::Success || s2 != gsStatus::Success)
                break;
            ++steps;
            const bool eqL = (to.lambda == r1.d.lambda());
            const bool eqU = ((to.U - r1.space.inner().committed()).norm() == 0);
            CHECK(eqL);
            CHECK(eqU);
            CHECK_EQUAL(r1.d.history().iterations, histOut.iterations);
            if (eqL && eqU) ++bitwise;
            from = to;
            hist = histOut;
        }
        CHECK_EQUAL(40, steps);
        CHECK_EQUAL(40, bitwise);
        CHECK_EQUAL(std::string(), contractViolation(r1.log, false));
        // the pure step leaves the driver's own state untouched
        CHECK(r2.d.lambda() == 0);
        CHECK(!r2.d.history().valid);
        CHECK_EQUAL(0, r2.d.context().step);
        gsInfo << "PureStep_equals_stateful: " << steps << " steps, bitwise-equal steps " << bitwise << "\n";
    }

    /// Quasi-Newton Riks with a valid history needs no predictor tangent and
    /// factorizes only before update 1, so every step after the first makes
    /// exactly one Jacobian assembly and one factorization. Were the version
    /// counter of a pure-step context to restart at 0, two consecutive pure
    /// steps would both factorize under version 1; the second factorize()
    /// would do nothing and the second pure step would run its corrector on the
    /// first one's K. The deltas of both counters per step, and the bitwise
    /// agreement with an independent stateful run, exclude that.
    TEST(PureStep_after_stateful_refactorizes)
    {
        UNITTEST_TIME_CONSTRAINT(1000);
        std::vector<real_t> refL;
        std::vector<gsVector<real_t> > refU;
        {
            TrussRig ref;
            ref.configureT1(CPtr(new gsRiksConstraint<real_t>), 0.05, false, true);
            for (index_t k = 0; k < 7; ++k)
            {
                CHECK(ref.d.step() == gsStatus::Success);
                refL.push_back(ref.d.lambda());
                refU.push_back(ref.space.inner().committed());
            }
        }

        TrussRig r;
        r.configureT1(CPtr(new gsRiksConstraint<real_t>), 0.05, false, true);
        index_t jac0 = 0, fact0 = 0;
        for (index_t k = 0; k < 3; ++k)
        {
            r.log.clear();
            CHECK(r.d.step() == gsStatus::Success);
            CHECK_EQUAL(std::string(), contractViolation(r.log, false));
            const index_t dj = r.jacCalls - jac0, df = r.d.linearSolver().numFactorizations() - fact0;
            CHECK_EQUAL(dj, df);
            CHECK_EQUAL(k == 0 ? 2 : 1, dj);
            jac0 = r.jacCalls;
            fact0 = r.d.linearSolver().numFactorizations();
            CHECK(r.d.lambda() == refL[k]);
        }

        gsPathPoint<real_t> from, to;
        from.U = r.space.inner().committed();
        from.lambda = r.d.lambda();
        gsStepHistory<real_t> hist = r.d.history(), histOut;
        for (index_t k = 3; k < 6; ++k)
        {
            r.log.clear();
            CHECK(r.d.step(from, hist, 0.05, to, histOut) == gsStatus::Success);
            CHECK_EQUAL(std::string(), contractViolation(r.log, true));
            const index_t dj = r.jacCalls - jac0, df = r.d.linearSolver().numFactorizations() - fact0;
            CHECK_EQUAL(1, dj);
            CHECK_EQUAL(1, df);
            jac0 = r.jacCalls;
            fact0 = r.d.linearSolver().numFactorizations();
            CHECK(to.lambda == refL[k]);
            CHECK((to.U - refU[k]).norm() == 0);
            from = to;
            hist = histOut;
        }

        r.d.setLambda(from.lambda);
        r.d.setHistory(hist);
        r.log.clear();
        CHECK(r.d.step() == gsStatus::Success);
        CHECK_EQUAL(std::string(), contractViolation(r.log, false));
        CHECK_EQUAL(1, r.jacCalls - jac0);
        CHECK_EQUAL(1, r.d.linearSolver().numFactorizations() - fact0);
        CHECK(r.d.lambda() == refL[6]);
        CHECK((r.space.inner().committed() - refU[6]).norm() == 0);
    }

    /// The tangent solver does no work for a repeated version: the solve then
    /// uses the stale factor. A new version recomputes. (Control of the trap
    /// described above.)
    TEST(TangentSolver_repeated_version_is_stale)
    {
        gsSparseMatrix<real_t> K1(2, 2), K2(2, 2);
        K1.insert(0, 0) = 2; K1.insert(1, 1) = 2; K1.makeCompressed();
        K2.insert(0, 0) = 4; K2.insert(1, 1) = 4; K2.makeCompressed();
        gsVector<real_t> b(2), x;
        b << 1, 1;

        gsTangentSolver<real_t> s;
        CHECK(s.factorize(K1, 7));
        CHECK(s.factorize(K2, 7));
        CHECK_EQUAL(1, s.numFactorizations());
        CHECK(s.solve(b, x));
        CHECK_CLOSE(0.5, x[0], 1e-14);
        CHECK_CLOSE(0.5, x[1], 1e-14);

        CHECK(s.factorize(K2, 8));
        CHECK_EQUAL(2, s.numFactorizations());
        CHECK(s.solve(b, x));
        CHECK_CLOSE(0.25, x[0], 1e-14);
        CHECK_CLOSE(0.25, x[1], 1e-14);
    }

    // ------------------------------------------------------------------
    // 7. Secant retry and drop, step history, solve
    // ------------------------------------------------------------------

    /// Load control with the secant predictor and MaxBisect 2; the residual
    /// fails at call 0 of selected attempts. Jacobian assemblies and
    /// factorizations are counted per attempt (a secant attempt needs no
    /// predictor tangent, a tangent attempt one), rejections by StepRejected.
    ///  (a) step 0 uses the tangent, steps 1-2 the secant;
    ///  (b) step 3: the secant attempt fails, the tangent retry at the same
    ///      length succeeds (depth stays 0);
    ///  (c) step 4: secant and tangent retry fail, the halved attempt succeeds
    ///      (depth 1, length and history step size 0.07/2); the failed tangent
    ///      retry assembled one K before its residual failed;
    ///  (d) step 5 uses the secant again;
    ///  (e) step 6: every attempt fails (secant, tangent retry, halved), the
    ///      step ends AssemblyError after 3 rejections; the next step, without
    ///      injection, uses the tangent: the drop of the secant held.
    TEST(Secant_retry_and_drop)
    {
        UNITTEST_TIME_CONSTRAINT(1000);
        TrussRig r;
        r.configure(CPtr(new gsLoadConstraint<real_t>), arcLengthGate(T1_tol, T1_tol, 1e-6, 100), 0.07);
        r.d.setSecantPredictor(true);
        r.d.controller().setMaxBisect(2);
        CHECK(r.d.setup());
        r.inject = [](index_t step, index_t attempt, index_t call)
        {
            return call == 0 && ((step == 3 && attempt == 1) || (step == 4 && attempt <= 2) || step == 6);
        };

        std::vector<std::vector<index_t> > jacs(8), facts(8);
        std::vector<index_t> its(8), rej(8);
        std::vector<gsStatus> sts(8);
        for (index_t k = 0; k < 8; ++k)
        {
            if (k == 7) r.inject = TrussProblem::Hook_t();   // step 6 keeps its index after the failure
            const size_t a0 = r.attemptJac.size();
            const index_t r0 = r.nRejected;
            sts[k] = r.d.step();
            its[k] = r.d.numIterations();
            rej[k] = r.nRejected - r0;
            for (size_t a = a0; a < r.attemptJac.size(); ++a)
            {
                jacs[k].push_back(r.attemptJac[a]);
                facts[k].push_back(r.attemptFact[a]);
            }
            if (k == 3) { CHECK_EQUAL(0, r.d.controller().depth()); }
            if (k == 4)
            {
                CHECK_EQUAL(1, r.d.controller().depth());
                CHECK(r.d.controller().length() == real_t(0.07)/2);
                CHECK(r.d.history().stepSize == real_t(0.07)/2);
            }
            if (k == 6)
            {
                CHECK(r.d.report().failedOp == gsDriverOp::Residual);
                CHECK_EQUAL(2, r.d.controller().depth());
                CHECK(r.d.controller().length() == real_t(0.07)/4);
            }
        }
        for (index_t k = 0; k < 8; ++k)
            CHECK(jacs[k] == facts[k]);   // one factorization per Jacobian, per attempt

        // (a)
        CHECK(sts[0] == gsStatus::Success); CHECK_EQUAL(0, rej[0]);
        CHECK_EQUAL(size_t(1), jacs[0].size()); CHECK_EQUAL(its[0] + 1, jacs[0][0]);
        for (index_t k = 1; k < 3; ++k)
        {
            CHECK(sts[k] == gsStatus::Success); CHECK_EQUAL(0, rej[k]);
            CHECK_EQUAL(size_t(1), jacs[k].size()); CHECK_EQUAL(its[k], jacs[k][0]);
        }
        // (b)
        CHECK(sts[3] == gsStatus::Success); CHECK_EQUAL(1, rej[3]);
        CHECK_EQUAL(size_t(2), jacs[3].size()); CHECK_EQUAL(0, jacs[3][0]); CHECK_EQUAL(its[3] + 1, jacs[3][1]);
        // (c)
        CHECK(sts[4] == gsStatus::Success); CHECK_EQUAL(2, rej[4]);
        CHECK_EQUAL(size_t(3), jacs[4].size());
        CHECK_EQUAL(0, jacs[4][0]); CHECK_EQUAL(1, jacs[4][1]); CHECK_EQUAL(its[4] + 1, jacs[4][2]);
        // (d)
        CHECK(sts[5] == gsStatus::Success); CHECK_EQUAL(0, rej[5]);
        CHECK_EQUAL(size_t(1), jacs[5].size()); CHECK_EQUAL(its[5], jacs[5][0]);
        // (e)
        CHECK(sts[6] == gsStatus::AssemblyError); CHECK_EQUAL(3, rej[6]);
        CHECK_EQUAL(size_t(3), jacs[6].size());
        CHECK_EQUAL(0, jacs[6][0]); CHECK_EQUAL(1, jacs[6][1]); CHECK_EQUAL(1, jacs[6][2]);
        CHECK(sts[7] == gsStatus::Success); CHECK_EQUAL(0, rej[7]);
        CHECK_EQUAL(size_t(1), jacs[7].size()); CHECK_EQUAL(its[7] + 1, jacs[7][0]);

        gsInfo << "Secant_retry_and_drop: jacobians per attempt";
        for (index_t k = 0; k < 8; ++k)
        {
            gsInfo << " | step " << k << " (it " << its[k] << ", rej " << rej[k] << "):";
            for (size_t a = 0; a < jacs[k].size(); ++a) gsInfo << " " << jacs[k][a];
        }
        gsInfo << "\n";
        CHECK_EQUAL(std::string(), contractViolation(r.log, false));
    }

    /// The step size recorded in the history is the size the attempt used,
    /// after a pre-step resize; solve(lambdaEnd) steps up to lambdaEnd, and
    /// with a step length that points away from lambdaEnd it returns Success
    /// without a step.
    TEST(Driver_pins_resize_and_solve)
    {
        UNITTEST_TIME_CONSTRAINT(1000);
        {
            TrussRig r;
            r.configureT1(CPtr(new gsLoadConstraint<real_t>), 0.07, false, false);
            r.d.setPreStep([](real_t, real_t & s) { s = 0.01; return gsHookDecision::Continue; });
            CHECK(r.d.step() == gsStatus::Success);
            CHECK(r.d.history().stepSize == real_t(0.01));
            CHECK(r.d.history().DeltaL == real_t(0.01));
            CHECK(r.d.controller().length() == real_t(0.07));
            CHECK_EQUAL(std::string(), contractViolation(r.log, false));
        }
        {
            TrussRig r;
            r.configureT1(CPtr(new gsLoadConstraint<real_t>), 0.07, false, false);
            CHECK(r.d.solve(0.14) == gsStatus::Success);
            CHECK_CLOSE(0.14, r.d.lambda(), 1e-12);
            CHECK_EQUAL(2, countChar(r.log, 'A'));

            const real_t L = r.d.lambda();
            const index_t started = countChar(r.log, 'S');
            const gsVector<real_t> U = r.space.inner().committed();
            CHECK(r.d.solve(-1) == gsStatus::Success);       // the step length is +0.07
            CHECK_EQUAL(-1, r.d.report().step);
            CHECK_EQUAL(started, countChar(r.log, 'S'));
            CHECK(r.d.lambda() == L);
            CHECK(r.space.inner().committed() == U);
            CHECK_EQUAL(std::string(), contractViolation(r.log, false));
        }
    }

    // ------------------------------------------------------------------
    // 8. Displacement control and the zero-load crossing
    // ------------------------------------------------------------------

    /// Prescribing the loaded dof v: v(w) = w + f(w)/k is strictly increasing,
    /// so displacement control follows the path through the limit point
    /// v* = w* + 1/5 = 0.6226 up to v = 0.8, where lambda = 0.73 and the
    /// LoadNorm basis of the residual test stays O(1).
    TEST(Displacement_passes_limit_point)
    {
        UNITTEST_TIME_CONSTRAINT(1000);
        const index_t N = 16;
        TrussRig r;
        r.configureT1(CPtr(new gsDisplacementConstraint<real_t>(1)), 0.05, false, false);
        std::vector<DRec> rec;
        real_t Lmax = 0;
        for (index_t k = 0; k < N; ++k)
        {
            const gsStatus st = r.d.step();
            CHECK(st == gsStatus::Success);
            const gsVector<real_t> & U = r.space.inner().committed();
            const DRec x = { r.d.lambda(), U.norm(), U[0], st, r.d.numIterations() };
            rec.push_back(x);
            // the prescribed dof is exact after every step
            CHECK(math::abs(U[1] - (k + 1)*real_t(0.05)) <= 1e-12);
            // and the state lies on the analytic path
            CHECK(math::abs(x.L - r.truss.lambdaOf(x.w)) <= 3*T1_tol*math::max(real_t(1), math::abs(x.L)));
            Lmax = math::max(Lmax, x.L);
            if (st != gsStatus::Success)
                break;
        }
        CHECK_EQUAL(size_t(N), rec.size());
        CHECK(Lmax >= 0.98*r.truss.lambdaMax());
        CHECK(passedLimit(rec, r.truss.wLimit()));
        CHECK_EQUAL(std::string(), contractViolation(r.log, false));
    }

    /// The path crosses lambda = 0 forward at w > 1: every step is a Success
    /// with w increasing (no reversal), also for the five steps after the
    /// crossing. The scaling of consistent Crisfield is frozen at its first
    /// step to |u_t|/|F| at U = 0 (0.4370973689763234, from the first-step
    /// tangent K(0)^{-1}F of the truss, logged in
    /// .claude/plans/gssa-redesign/tasks/logs/792-B00b2a-cells.log:350) and
    /// stays bitwise constant.
    TEST(CC_crosses_zero_load_forward)
    {
        UNITTEST_TIME_CONSTRAINT(1000);
        const real_t steps[] = { 0.2, 0.4 };
        const real_t phi0 = 0.4370973689763234;
        for (size_t c = 0; c < 2; ++c)
            for (index_t q = 0; q < 2; ++q)
            {
                const bool quasi = (q == 1);
                memory::shared_ptr<gsConsistentCrisfieldConstraint<real_t> > cc(new gsConsistentCrisfieldConstraint<real_t>);
                TrussRig r;
                r.configureT1(cc, steps[c], false, quasi);
                real_t wprev = 0;
                index_t crossing = -1, after = 0;
                std::vector<DRec> rec;
                for (index_t k = 0; k < 30; ++k)
                {
                    const gsStatus st = r.d.step();
                    CHECK(st == gsStatus::Success);
                    if (st != gsStatus::Success)
                        break;
                    const real_t L = r.d.lambda();
                    const real_t w = r.space.inner().committed()[0];
                    const DRec x = { L, r.space.inner().committed().norm(), w, st, r.d.numIterations() };
                    rec.push_back(x);
                    CHECK(w > wprev);
                    if (!quasi)
                        CHECK(math::abs(L - r.truss.lambdaOf(w)) <= 3*T1_tol*math::max(real_t(1), math::abs(L)));
                    wprev = w;
                    if (k == 0)
                    {
                        CHECK(cc->hasFrozenScaling());
                        CHECK_CLOSE(phi0, cc->scaling(), 1e-12*phi0);
                    }
                    CHECK(cc->scaling() == cc->frozenScaling());
                    if (crossing < 0)
                    {
                        if (L < 0 && w > 1)
                            crossing = k;
                    }
                    else if (++after == 5)
                        break;
                }
                gsInfo << "CC zero-load crossing: ds " << steps[c] << ", quasi " << quasi
                       << ", crossing step " << crossing << ", steps after " << after << "\n";
                CHECK(crossing >= 0);
                CHECK_EQUAL(5, after);
                CHECK(crossesForward(rec));
                CHECK_EQUAL(std::string(), contractViolation(r.log, false));
            }

        // information only: Crisfield's automatic scaling follows the secant ratio |U|/(|lambda||F|),
        // which is not defined at lambda = 0; no claim is made about this run
        {
            TrussRig r;
            r.configureT1(CPtr(new gsCrisfieldConstraint<real_t>), 0.2, false, false);
            const std::vector<DRec> rec = r.run(30);
            gsInfo << "Crisfield (automatic scaling) ds 0.2 Q0: " << rec.size() << " steps recorded, crossesForward "
                   << (crossesForward(rec) ? "yes" : "no") << "\n";
        }
    }

    // ------------------------------------------------------------------
    // 9. Constraint and corrector units (no driver)
    // ------------------------------------------------------------------

    /// Crisfield with complex roots reads the residual R. n = 3, phi = 0 (so
    /// psi = 0), s = 1, F = (0,0,1), lambda = 0.5; the iterate is Delta U =
    /// (2,0,0), Delta lambda = 0.1; dUt = (0,0,1), dUbar = (0,1,0). The
    /// quadratic for the load update has the discriminant -16 (no real roots)
    /// and the modified one -192, so the simplified rescaling is
    /// taken: the Newton point (Delta U + dUbar, Delta lambda_cr) is scaled onto
    /// the sphere |Delta U| = s. With F_int = (lambda + Delta lambda) F - R and
    /// lambda_cr = F_int.F/F.F:
    ///  (i) R = (0,0,0.3): F_int = (0,0,0.3), lambda_cr = 0.3, Delta lambda_cr
    ///      = -0.2, the Newton point is ((2,1,0), -0.2) with |(2,1,0)| = sqrt 5,
    ///      so mu = 1/sqrt 5, the new increment (2,1,0)/sqrt 5, -0.2/sqrt 5 and
    ///      deltaU = (2/sqrt 5 - 2, 1/sqrt 5, 0), deltaL = -0.2/sqrt 5 - 0.1;
    ///  (ii) the same start with R = 0: F_int = 0.6 F, lambda_cr = 0.6,
    ///      Delta lambda_cr = 0.1, deltaL = 0.1/sqrt 5 - 0.1.
    /// The two deltaL differ, so the fallback depends on R.
    TEST(Crisfield_complex_fallback_reads_residual)
    {
        const real_t tol = 1e-14;
        const real_t r5 = std::sqrt(real_t(5));
        gsCrisfieldConstraint<real_t> c;
        c.setScaling(0);
        gsVector<real_t> dUt(3), dUbar(3);
        dUt << 0, 0, 1;
        dUbar << 0, 1, 0;

        for (int which = 0; which < 2; ++which)
        {
            gsStepContext<real_t> ctx;
            ctx.F.resize(3);      ctx.F << 0, 0, 1;
            ctx.DeltaU.resize(3); ctx.DeltaU << 2, 0, 0;
            ctx.deltaU = gsVector<real_t>::Zero(3);
            ctx.R = gsVector<real_t>::Zero(3);
            if (which == 0) ctx.R << 0, 0, 0.3;
            ctx.lambda = 0.5; ctx.DeltaL = 0.1; ctx.deltaL = 0.1;
            ctx.stepSize = 1;
            ctx.prev.valid = false;
            ctx.relaxation = 1;

            CHECK(c.correct(ctx, dUbar, dUt));
            CHECK(c.lastRoot() == gsCrisfieldConstraint<real_t>::Root::Simplified);
            if (which == 0)
            {
                CHECK_CLOSE(-0.2/r5 - 0.1, ctx.deltaL, tol);
                CHECK_CLOSE(2/r5 - 2, ctx.deltaU[0], tol);
                CHECK_CLOSE(1/r5,     ctx.deltaU[1], tol);
                CHECK_CLOSE(0,        ctx.deltaU[2], tol);
            }
            else
                CHECK_CLOSE(0.1/r5 - 0.1, ctx.deltaL, tol);
            ctx.DeltaU += ctx.deltaU;
            ctx.DeltaL += ctx.deltaL;
            CHECK(math::abs(c.residual(ctx)) < 1e-14);   // s = 1
        }
        CHECK_EQUAL(2, c.numSimplifiedRoots());
    }

    /// Conditional relaxation: n = 2, phi = 0, s = 1, F = (0,1), Delta U =
    /// (0.6, 0), Delta lambda = 0, R = 0, previous increment (0,-1), dUt =
    /// (0,1), dUbar = 0, previous load update lamold = 1. The real roots of the
    /// quadratic are delta lambda = +-0.8 (|(0.6, delta lambda)| = 1); the
    /// angle criterion takes the root whose update has the larger inner
    /// product with the previous increment (0,-1): -0.8 gives +0.8, +0.8 gives
    /// -0.8, so -0.8 is taken. Its sign differs from lamold without growth in
    /// magnitude, so an oscillation is detected: with relaxation omega = 0.5
    /// the constraint overwrites ctx.relaxation by 0.5, with omega = 1 it
    /// leaves it alone.
    TEST(Crisfield_conditional_relaxation)
    {
        const real_t tol = 1e-14;
        for (int which = 0; which < 2; ++which)
        {
            gsCrisfieldConstraint<real_t> c;
            c.setScaling(0);
            c.setRelaxation(which == 0 ? real_t(0.5) : real_t(1));
            gsStepContext<real_t> ctx;
            ctx.F.resize(2);      ctx.F << 0, 1;
            ctx.DeltaU.resize(2); ctx.DeltaU << 0.6, 0;
            ctx.deltaU = gsVector<real_t>::Zero(2);
            ctx.R = gsVector<real_t>::Zero(2);
            ctx.DeltaL = 0; ctx.lambda = 0; ctx.stepSize = 1;
            ctx.prev.DeltaU.resize(2); ctx.prev.DeltaU << 0, -1;
            ctx.prev.DeltaL = 0; ctx.prev.stepSize = 1; ctx.prev.valid = true;
            ctx.deltaL = 1;
            ctx.relaxation = 0.7;
            gsVector<real_t> dUt(2), dUbar(2);
            dUt << 0, 1;
            dUbar << 0, 0;

            CHECK(c.correct(ctx, dUbar, dUt));
            CHECK(c.lastRoot() == gsCrisfieldConstraint<real_t>::Root::Real);
            CHECK_CLOSE(which == 0 ? 0.5 : 0.7, ctx.relaxation, tol);
            CHECK_CLOSE(-0.8, ctx.deltaL, tol);
            CHECK_CLOSE(0,    ctx.deltaU[0], tol);
            CHECK_CLOSE(-0.8, ctx.deltaU[1], tol);
            CHECK_EQUAL(which == 0 ? 1 : 0, c.numRelaxed());
            if (which == 1)
            {
                ctx.DeltaU += ctx.deltaU;
                ctx.DeltaL += ctx.deltaL;
                CHECK(math::abs(c.residual(ctx)) < 1e-14);   // |(0.6, -0.8)| = 1
            }
        }
    }

    /// K is assembled before update k iff !quasi || k <= 1 || (interval > 0 &&
    /// k % interval == 0); always for k <= 0
    TEST(Corrector_refactorize_truth_table)
    {
        struct Row { bool quasi; index_t interval; const char * expect; };   // k = 1..5
        const Row rows[] = { { false, -1, "TTTTT" }, { true, -1, "TFFFF" }, { true, 0, "TFFFF" },
                             { true,   1, "TTTTT" }, { true,  2, "TTFTF" }, { true, 3, "TFTFF" } };
        for (size_t i = 0; i < sizeof(rows)/sizeof(rows[0]); ++i)
            for (index_t k = 1; k <= 5; ++k)
                CHECK_EQUAL(rows[i].expect[k-1] == 'T',
                            gsNewtonCorrector<real_t>::refactorizeAt(k, rows[i].quasi, rows[i].interval));
        CHECK(gsNewtonCorrector<real_t>::refactorizeAt(0, true, -1));
    }
}
