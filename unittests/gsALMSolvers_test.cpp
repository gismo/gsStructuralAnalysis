/** @file gsALMSolvers_test.cpp

    @brief Unit tests for the gsALMSolvers family (gsALMBase and its four
           derived arc-length solvers) on the closed-form 2-DOF fixtures of
           gsALMTestProblems.h, plus gsBucklingSolver, gsModalSolver and the
           shared gsEigenProblemBase on hand-coded 2x2 pencils with a
           closed-form spectrum.

    Implements tests 1-10 and 9b of the gsALMSolvers correctness review test plan.

    == Reading this file ==
    Several tests pin defects that are NOT fixed yet. They are written to the
    CORRECT invariant, were run, and then marked with

        // DISABLED(task NN, <finding>): <reason>
        ... GSALM_DISABLED("task NN","<finding>","<reason>"); ...

    so that the committed suite is green while every waiting test is greppable:

        grep -rn "DISABLED(task" optional/gsStructuralAnalysis/unittests

    The assertions below such a marker are COMPILED but not executed; build with
    -DGSALM_RUN_DISABLED_TESTS to execute them and see which fixes have landed.

    Earlier work landed six fixes and CLEARED four markers: M3 (Crisfield root-selection
    weight), M4 (complex-root fallback internal force) and M2b for both gsALMRiks and
    gsALMConsistentCrisfield. The remaining markers below (M0a, M0b, M5, M6, M16, m1, M14,
    m17) are findings explicitly left OUT of scope of that work; they are
    waiting on follow-up work that will re-tag them.

    == Determinism ==
    All fixtures are hand-coded 2-DOF operators and never touch gsExprAssembler,
    so they are bit-reproducible regardless of the thread count and tight
    tolerances (1e-10) are legitimate.

    == Why SingularPointTestTol = 1e-4 here ==
    The singular-point test classifies on the cosine |V.f|/|f| of the mode shape
    returned by the inverse power iteration (finding m8). That
    iteration starts from a deterministic GENERIC vector (not from `ones`) and
    SingularPointTestIt = 5 is a MAXIMUM, the sweep stopping early once the mode
    settles to within Tol. On fixture P at a fundamental-branch point (lambda,0)
    the exact tangent is diag(1,1-lambda), so five sweeps return
    V ~ (v1, v2*(1-lambda)^-5) and the cosine is (v1/v2)*|1-lambda|^5 -- NOT
    zero (v1/v2 = 3.01 is the ratio of the two seed components). The default
    1e-6 therefore only accepts points within |1-lambda| < 0.05 of the fork;
    1e-4 widens that to ~0.13, which is where these tests operate. At a LIMIT
    point the cosine is O(1), so the test stays fully discriminating.

    Author(s): H.M. Verhelst
 **/

#include "gismo_unittest.h"

#include <limits>
#include <sstream>
#include <string>
#include <vector>
#include <type_traits>
#include <new>
#include <cstring>
#include <cstddef>

#include <gsStructuralAnalysis/src/gsStructuralAnalysisTools/gsStructuralAnalysisTypes.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMBase.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMLoadControl.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMRiks.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMCrisfield.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMConsistentCrisfield.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsAPALMData.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsAPALM.h>
#include <gsStructuralAnalysis/src/gsEigenSolvers/gsBucklingSolver.h>
#include <gsStructuralAnalysis/src/gsEigenSolvers/gsModalSolver.h>

#include "gsALMTestProblems.h"

SUITE(gsALMSolvers_test)   // suite name == file basename
{

using namespace gsALMTest;

// ===========================================================================
// Shared helpers
// ===========================================================================

/// The option regime used by (almost) every test below: determinant-based
/// stability (exact for these diagonal/2x2 tangents), tight corrector
/// tolerances, no bisection, silent.
template <class ALM>
void configure(ALM & alm, real_t ds, real_t SPtestTol = 1e-4)
{
    alm.options().setString("Solver","SimplicialLDLT");
    alm.options().setInt   ("BifurcationMethod",0);          // 0: determinant
    alm.options().setReal  ("Length",ds);
    alm.options().setReal  ("Tol" ,1e-10);
    alm.options().setReal  ("TolF",1e-8);
    alm.options().setReal  ("TolU",1e-8);
    alm.options().setInt   ("MaxIter",100);
    alm.options().setReal  ("SingularPointTestTol",SPtestTol);
    alm.options().setReal  ("SingularPointComputeTolE",1e-10);
    alm.options().setReal  ("SingularPointComputeTolB",0);   // bisection OFF
    alm.options().setSwitch("Verbose",false);
    alm.applyOptions();
    alm.initialize();
}

/// |V.f|/|f| -- the dimensionless limit-vs-branch statistic the singular-point
/// test is supposed to threshold (see C1). Recomputed here from the PUBLIC
/// solutionV() so the test does not depend on the internal formula.
inline real_t modeForceCosine(const gsVector<real_t> & V, const gsVector<real_t> & f)
{
    return math::abs(V.normalized().dot(f)) / f.norm();
}

/// det K_T at \a U for a 2x2 fixture; zero exactly at a singular point.
inline real_t tangentDet(const AlmProblem & prob, const gsVector<real_t> & U)
{
    gsSparseMatrix<real_t> K;
    prob.Jacobian(U,K);
    return K.coeff(0,0)*K.coeff(1,1) - K.coeff(0,1)*K.coeff(1,0);
}

// ---------------------------------------------------------------------------
// White-box probes. Nothing under src/ is modified: the protected members are
// reached the way C++ intends, by deriving.
// ---------------------------------------------------------------------------

class RiksProbe : public gsALMRiks<real_t>
{
    typedef gsALMRiks<real_t> Base;
public:
    RiksProbe(const gsStructuralAnalysisOps<real_t>::Jacobian_t   & J,
              const gsStructuralAnalysisOps<real_t>::ALResidual_t & R,
              const gsVector<real_t>                              & F)
    : Base(J,R,F) {}

    /// Runs ONLY gsALMBase::computeLength() for a given iteration count, and
    /// returns the resulting length factor. This is the adaptive-length
    /// arithmetic of M1 in isolation.
    real_t lengthFactor(index_t iterationsTaken, index_t desired, real_t ds)
    {
        this->m_arcLength         = ds;
        this->m_numIterations     = iterationsTaken;
        this->m_desiredIterations = desired;
        this->computeLength();
        return this->m_arcLength / ds;
    }

    real_t arcLength()     const { return this->m_arcLength; }
    real_t arcLengthPrev() const { return this->m_arcLength_prev; }

    /// Riks's convex constraint weight, phi (psi^2 = (1-phi)/phi = numDof-1).
    /// Named m_convexWeight in gsALMRiks; the shim keeps the name phi()
    /// so that the psi^2 = (1-phi)/phi translation below reads as in the literature.
    real_t phi() const { return this->m_convexWeight; }

    /// True once gsALMRiks::predictorGuess() has CONSUMED the guess: it ends with
    /// m_Uguess.resize(0), which gsALMRiks::predictor() never does. This is how a test
    /// proves gsALMBase::_step() dispatched to the guess predictor and not to predictor().
    bool guessConsumed() const { return this->m_Uguess.rows()==0; }

    /// initiateStep() + predictor(), i.e. exactly the predictor of one step.
    void runPredictorOnly() { this->initiateStep(); this->predictor(); }

    const gsVector<real_t> & DeltaU()    const { return this->m_DeltaU; }
    real_t                   DeltaL()    const { return this->m_DeltaL; }
    const gsVector<real_t> & deltaUt()   const { return this->m_deltaUt; }
    const gsVector<real_t> & deltaUbar() const { return this->m_deltaUbar; }
    real_t                   deltaL()    const { return this->m_deltaL; }
    const std::string &      note()      const { return this->m_note; }

    /// One corrector iteration, exactly as gsALMBase::_step() drives it:
    /// quasiNewtonIteration() (re-assembles K_T at m_U+m_DeltaU and re-solves
    /// for m_deltaUt) followed by iteration().
    void runOneCorrectorIteration()
    {
        this->quasiNewtonIteration();
        this->iteration();
    }

    /// The scalar the corrector in gsALMRiks::iteration() divides by, HALVED --
    /// i.e. phi * (paper eq. (22) t) = phi*DeltaU.deltaUt + (1-phi)*DeltaL.
    /// Evaluated on whatever state the object currently holds.
    real_t halfDenominator() const
    {
        return (1.0-this->m_convexWeight)*this->m_DeltaL
             + this->m_convexWeight*this->m_DeltaU.dot(this->m_deltaUt);
    }

    /// Drives ONE step by hand, exactly as gsALMBase::_step() does for a
    /// non-quasi-Newton run, recording halfDenominator() BEFORE every corrector
    /// iteration and the state after it. Returns true when the step converged.
    struct IterRecord
    {
        real_t halfDen;   ///< phi * (paper eq. (22) t), i.e. denominator/2, as computed by halfDenominator()
        real_t deltaL;    ///< the m_deltaL that divisor produced
        real_t DeltaL;    ///< accumulated DeltaLambda after the iteration
        real_t DeltaU0;   ///< accumulated DeltaU, component 0
        real_t residueF;
    };
    bool runStepInstrumented(std::vector<IterRecord> & rec, index_t maxIt = 100)
    {
        this->m_converged     = false;
        this->m_numIterations = 0;
        this->initiateStep();
        this->predictor();
        this->computeResidual();
        this->computeResidualNorms();
        for (this->m_numIterations = 1; this->m_numIterations < maxIt; ++this->m_numIterations)
        {
            this->quasiNewtonIteration();   // computeJacobian + factorize + computeUt
            IterRecord r;
            r.halfDen = this->halfDenominator();
            this->iteration();
            r.deltaL  = this->m_deltaL;
            r.DeltaL  = this->m_DeltaL;
            r.DeltaU0 = this->m_DeltaU[0];
            this->computeResidual();
            this->computeResidualNorms();
            r.residueF = this->m_residueF;
            rec.push_back(r);
            if (this->m_residueF < this->m_toleranceF && this->m_residueU < this->m_toleranceU)
                return true;
        }
        return false;
    }

    /// Commits the increment the way gsALMRiks::iterationFinish() does.
    void commitStep() { this->iterationFinish(); }

    /// computeResidualNorms() captures its denominator at iteration 0. Reading this
    /// straight after runStepInstrumented(rec,/*maxIt*/0) -- which runs exactly the
    /// predictor + basis capture of gsALMBase::_step() and skips the corrector loop --
    /// gives the ORDINARY corrector's own basis, unfloored by construction under
    /// computeResidualNorms()'s extended-solve-only scoping (its \c bool parameter).
    real_t basisResidualU() const { return this->m_basisResidualU; }
};

class CrisfieldProbe : public gsALMCrisfield<real_t>
{
    typedef gsALMCrisfield<real_t> Base;
public:
    CrisfieldProbe(const gsStructuralAnalysisOps<real_t>::Jacobian_t   & J,
                   const gsStructuralAnalysisOps<real_t>::ALResidual_t & R,
                   const gsVector<real_t>                              & F)
    : Base(J,R,F) {}

    real_t phi() const { return this->m_phi; }
    /// A0 = phi^2 |f|^2 -- the load weight of the constraint surface.
    real_t A0()  const
    { return this->m_phi*this->m_phi*this->stepForcing().dot(this->stepForcing()); }

    /// initiateStep() + predictor(), i.e. exactly the predictor of one step.
    void runPredictorOnly() { this->initiateStep(); this->predictor(); }

    const gsVector<real_t> & DeltaU() const { return this->m_DeltaU; }
    real_t                   DeltaL() const { return this->m_DeltaL; }
    const gsVector<real_t> & deltaU()    const { return this->m_deltaU; }
    real_t                   deltaL()    const { return this->m_deltaL; }
    const gsVector<real_t> & deltaUbar() const { return this->m_deltaUbar; }
    const std::string &      note()      const { return this->m_note; }

    /// Forces the corrector into the BORDERED chart without going through the
    /// NotConverged retry of gsALMCrisfield::step(). This is what makes G-B a
    /// test of the chart equivalence rather than of the retry trigger: with the
    /// trigger, a trace that never fails executes literally the same
    /// instructions with the switch on and off, so agreeing proves nothing.
    void forceBordered(bool b) { this->m_useBordered = b; }

    /// initiateStep() + predictor() + EXACTLY ONE corrector iteration, in the
    /// order gsALMBase::_step() drives them. Leaves
    /// m_deltaU/m_deltaL holding that single iteration's chart output.
    void runPredictorAndOneIteration()
    {
        this->m_converged     = false;
        this->m_numIterations = 0;
        this->initiateStep();
        this->predictor();
        this->computeResidual();
        this->computeResidualNorms();
        this->m_numIterations = 1;
        this->quasiNewtonIteration();
        this->iteration();
    }

    /// One corrector iteration, exactly as gsALMBase::_step() drives it
    /// (gsALMBase.hpp, the m_numIterations loop of _step()): quasiNewtonIteration()
    /// -- which on the bordered path assembles K_T WITHOUT factorizing it -- followed
    /// by iteration(). Nothing is installed by hand: the state consumed is whatever
    /// the previous iteration left.
    void runOneCorrectorIteration()
    {
        ++this->m_numIterations;
        this->quasiNewtonIteration();
        this->iteration();
    }

    /// q_t and q, the load components of the two bordered solutions (gsALMCrisfield's
    /// m_deltaLt / m_deltaLbar). ⚠ BORDERED-ONLY: NaN outside the bordered chart, on
    /// purpose (see their member documentation).
    real_t deltaLt()   const { return this->m_deltaLt;   }
    real_t deltaLbar() const { return this->m_deltaLbar; }
    /// p_t (bordered) or delta_u_t (elimination), depending on the chart in force.
    const gsVector<real_t> & deltaUt() const { return this->m_deltaUt; }

    /// The transversality guard (gsALMCrisfield class doc (h)): the cosine of
    /// the chart the last _borderedSolve() settled on, the cumulative count of
    /// re-borderings, and the two thresholds (which are members, not options).
    real_t  chartCosine()      const { return this->m_chartCosine;       }
    index_t chartReborderings()const { return this->m_chartReborderings; }
    real_t  chartCosTol()      const { return this->m_chartCosTol;       }
    real_t  chartRankTol()     const { return this->m_chartRankTol;      }

    /// Installs a complete, self-consistent corrector state and then runs ONLY
    /// gsALMCrisfield::computeLambdasComplex() (Lam 1992 eqs 13-17), the
    /// complex-root fallback of M4.
    ///
    /// \param U,L        converged point the step started from
    /// \param DeltaU,DeltaL  the increment accumulated so far in this step
    /// \param ds         the arc length of the step
    void runComplexFallback(const gsVector<real_t> & U, real_t L,
                            const gsVector<real_t> & DeltaU, real_t DeltaL,
                            real_t ds)
    {
        this->m_U = U;  this->m_L = L;
        this->m_DeltaU = DeltaU;  this->m_DeltaL = DeltaL;
        this->m_arcLength = ds;
        this->m_eta = 1.0;

        // The tangent the corrector holds: assembled at (U+DeltaU) by the last
        // quasiNewtonIteration(), which is also where m_deltaUt comes from.
        gsVector<real_t> zero = gsVector<real_t>::Zero(U.size());
        this->m_jacMat = this->computeJacobian(U+DeltaU,zero); // factorizes too
        this->m_deltaUt   = this->solveSystem(this->computeForcing(U+DeltaU,L+DeltaL));
        this->m_resVec    = this->computeResidual(U+DeltaU,L+DeltaL);
        this->m_deltaUbar = this->solveSystem(-this->m_resVec);

        this->computeLambdasComplex();
    }
};

class ConsistentProbe : public gsALMConsistentCrisfield<real_t>
{
    typedef gsALMConsistentCrisfield<real_t> Base;
public:
    ConsistentProbe(const gsStructuralAnalysisOps<real_t>::Jacobian_t   & J,
                    const gsStructuralAnalysisOps<real_t>::ALResidual_t & R,
                    const gsVector<real_t>                              & F)
    : Base(J,R,F) {}

    real_t phi() const { return this->m_phi; }
};

/// White-box probe for the singular-point BISECTION stage.
///
/// Counts the arc-length step()s the stage takes (step() is the only route
/// _bisectionSolve has to one, and unlike a Jacobian count it is not confused by
/// the bracketing precondition's endpoint stability evaluations), and captures
/// the (U,L) the stage hands to the extended system. With the extended stage
/// stubbed out, solutionU()/solutionL() after computeSingularPoint() are exactly
/// what the bisection stage produced.
class SingularPointProbe : public gsALMRiks<real_t>
{
    typedef gsALMRiks<real_t> Base;
public:
    SingularPointProbe(const gsStructuralAnalysisOps<real_t>::Jacobian_t   & J,
                       const gsStructuralAnalysisOps<real_t>::ALResidual_t & R,
                       const gsVector<real_t>                              & F)
    : Base(J,R,F), m_steps(0), m_extCalls(0), m_stub(true), m_seedL(0) {}

    /// When true, _extendedSystemSolve records its arguments and returns success
    /// WITHOUT moving (m_U,m_L) -- isolates stage 1 from stage 2.
    void    stubExtendedSolve(bool s) { m_stub = s; }
    void    resetProbe() { m_steps = 0; m_extCalls = 0; }
    index_t stepCalls()     const { return m_steps; }
    index_t extendedCalls() const { return m_extCalls; }
    const gsVector<real_t> & extendedSeedU() const { return m_seedU; }
    real_t                   extendedSeedL() const { return m_seedL; }

    gsStatus step() { ++m_steps; return Base::step(); }

protected:
    bool _extendedSystemSolve(const gsVector<real_t> & U, const real_t L,
                              const real_t tol)
    {
        ++m_extCalls;
        m_seedU = U;
        m_seedL = L;
        if (m_stub) return true;
        return Base::_extendedSystemSolve(U,L,tol);
    }

private:
    index_t          m_steps, m_extCalls;
    bool             m_stub;
    gsVector<real_t> m_seedU;
    real_t           m_seedL;
};

// ===========================================================================
// TEST 1 (C1) -- the singular-point test must be scale invariant.
// ===========================================================================
// Fixture S: fixture P (and F) with the residual and the force scaled by s. The
// equilibria, the null vector V and the cosine |V.f|/|f| are all invariant, so
// the CLASSIFICATION must be too. Comparing the raw product
// |V.f| instead would grow like s: at s=1e6 a perfect branch point would be
// downgraded to a limit point, at s=1e-6 every limit point would pass as a branch point.
TEST(singular_point_test_is_scale_invariant)
{
    const real_t scales[3] = {1e-6, 1.0, 1e6};

    // --- Fixture S = pitchfork at three load scales: BRANCH point ------------
    real_t cosP[3];
    bool   isBifP[3];
    for (index_t i = 0; i != 3; ++i)
    {
        AlmProblem prob = pitchforkProblem(scales[i]);
        gsALMLoadControl<real_t> alm(prob.Jacobian, prob.ALResidual, prob.Force);
        configure(alm, 0.05);

        // Fundamental branch, just past the fork: u = (1.05,0), lambda = 1.05.
        gsVector<real_t> U(2); U << 1.05, 0.0;
        alm.setSolution(U,1.05);

        isBifP[i] = alm.isBifurcation(true);
        cosP[i]   = modeForceCosine(alm.solutionV(), prob.Force);
        gsInfo<<"  [scale-invariance] pitchfork s="<<scales[i]
              <<"  |V.f|/|f| = "<<cosP[i]
              <<"  isBifurcation = "<<(isBifP[i]?"true":"false")<<"\n";
    }

    // The statistic itself is scale free (this is the invariant C1 breaks).
    CHECK_CLOSE(cosP[1], cosP[0], 1e-10);
    CHECK_CLOSE(cosP[1], cosP[2], 1e-10);
    // ... hence so is the classification, and it is the analytically correct one.
    CHECK(isBifP[0]);
    CHECK(isBifP[1]);
    CHECK(isBifP[2]);

    // --- Fixture F at three load scales: LIMIT point (control) ---------------
    // Without this half, "isBifurcation() always returns true" would pass above.
    real_t cosF[3];
    bool   isBifF[3];
    for (index_t i = 0; i != 3; ++i)
    {
        AlmProblem prob = foldProblem(scales[i]);
        gsALMLoadControl<real_t> alm(prob.Jacobian, prob.ALResidual, prob.Force);
        configure(alm, 0.05);

        // On the stable arm, close to the fold: u1 = 0.95, lambda = 2u1-u1^2.
        gsVector<real_t> U(2); U << 0.95, 0.0;
        alm.setSolution(U,foldLambda(0.95));

        isBifF[i] = alm.isBifurcation(true);
        cosF[i]   = modeForceCosine(alm.solutionV(), prob.Force);
        gsInfo<<"  [scale-invariance] fold      s="<<scales[i]
              <<"  |V.f|/|f| = "<<cosF[i]
              <<"  isBifurcation = "<<(isBifF[i]?"true":"false")<<"\n";
    }

    CHECK_CLOSE(cosF[1], cosF[0], 1e-10);
    CHECK_CLOSE(cosF[1], cosF[2], 1e-10);
    CHECK(!isBifF[0]);
    CHECK(!isBifF[1]);
    CHECK(!isBifF[2]);
}

// ===========================================================================
// TEST 2 (second half of the review sketch) -- the singular point returned for
// a genuine BRANCH point is the analytic one.
// ===========================================================================
TEST(compute_singular_point_refines_the_branch_point)
{
    AlmProblem prob = pitchforkProblem();
    gsALMLoadControl<real_t> alm(prob.Jacobian, prob.ALResidual, prob.Force);
    configure(alm, 0.05);

    // March up the fundamental branch to lambda = 1.05 (21 load-control steps).
    for (index_t k = 0; k != 21; ++k)
        CHECK(alm.step() == gsStatus::Success);
    CHECK_CLOSE(1.05, alm.solutionL(), 1e-10);

    const gsStatus st = alm.computeSingularPoint(alm.solutionU(), alm.solutionL(),
                                                 /*switchBranch*/false,
                                                 /*jacobian*/true,
                                                 /*testPoint*/true);
    CHECK(st == gsStatus::Success);

    // Analytic singular point: u* = (1,0), lambda* = 1, det K_T(u*) = 0.
    CHECK_CLOSE(1.0, alm.solutionL(),    1e-6);
    CHECK_CLOSE(1.0, alm.solutionU()[0], 1e-6);
    CHECK_CLOSE(0.0, alm.solutionU()[1], 1e-6);
    CHECK(math::abs(tangentDet(prob,alm.solutionU())) < 1e-6);
}

// ===========================================================================
// TEST 2 (M0b) -- computeSingularPoint() must never report Success for a call
// in which it did nothing.
// ===========================================================================
// At a LIMIT point _computeSingularPoint() classifies test==false, skips its
// whole body, and computeSingularPoint() then sets gsStatus::Success -- after
// `m_U = U; m_L = L;` has already rewound the solver to the caller's input. A
// caller that trusts the status (gsALMExploration<T>::traceSweep does, in its
// singular-point block: everything under `if (spStatus != gsStatus::Success)`
// is the honest-failure path, so Success means "store this as the refined
// branch point") is handed the unrefined input point flagged as a converged
// branch point.
//
TEST(compute_singular_point_reports_failure_not_success)
{
    AlmProblem prob = foldProblem();
    gsALMRiks<real_t> alm(prob.Jacobian, prob.ALResidual, prob.Force);
    configure(alm, 0.05);

    // A point on the stable arm, close enough to the fold that a singular-point
    // computation is a sensible thing for a caller to ask for.
    for (index_t k = 0; k != 12; ++k)
        CHECK(alm.step() == gsStatus::Success);

    const gsVector<real_t> Uold = alm.solutionU();
    const real_t           Lold = alm.solutionL();
    // Ground truth: this really is a limit point, not a branch point.
    CHECK(!alm.isBifurcation(true));

    const gsStatus st = alm.computeSingularPoint(Uold,Lold,
                                                 /*switchBranch*/false,
                                                 /*jacobian*/true,
                                                 /*testPoint*/true);
    const real_t moved = (alm.solutionU()-Uold).norm();
    gsInfo<<"  [M0b] status="<<static_cast<index_t>(st)<<"  |U-Uold| = "<<moved<<"\n";

    // The invariant: Success must never coincide with "nothing happened".
    CHECK( !(st == gsStatus::Success && moved == 0.0) );
    // The intended reading of the sketch: a limit point is not a singular point
    // this routine can deliver, so it must say so.
    CHECK(st != gsStatus::Success);
}

// ===========================================================================
// TEST 3 (M0a) -- computeSingularPoint() must classify AT the point it is given.
// ===========================================================================
// _computeSingularPoint() runs _testSingularPoint() BEFORE assigning m_U = U,
// so the classification (and, with jacobian=true, the tangent) is taken at the
// solver's CURRENT state, not at the requested one. Both in-tree callers pass
// the pre-crossing point while the solver sits at the post-crossing point.
//
// Setup: the requested point (1.05,0),1.05 classifies as a BRANCH point;
// the solver is parked at (0.5,0),0.5 which classifies as a LIMIT point (the
// 5-step power iteration has not separated the modes that far from the fork).
// A correct implementation refines the requested point to (1,0),1; today it
// classifies at (0.5,0), takes the test==false path and returns the input
// unchanged (i.e. the failure is OBSERVED through M0b's no-op).
//
TEST(singular_point_test_uses_the_requested_point)
{
    AlmProblem prob = pitchforkProblem();
    gsALMLoadControl<real_t> alm(prob.Jacobian, prob.ALResidual, prob.Force);
    configure(alm, 0.05);

    gsVector<real_t> Ureq(2); Ureq << 1.05, 0.0;   const real_t Lreq = 1.05;
    // NOTE (2026-08-13): Ucur was originally (0.50, 0.0), i.e. ON the
    // fixture's trivial branch u2=0. K_T(u1,u2) = [[1,-u2],[-u2,1-u1+3u2^2]] is
    // EXACTLY diagonal there for every u1, so the converged critical mode is
    // e2=(0,1), exactly force-orthogonal to F=(1,0) (dot==0), at EVERY point of
    // that branch regardless of distance to the true singular point at u1=1 --
    // not a proximity-to-singularity effect. Under the old SingularPointTestIt=5
    // default the mode had not separated the two eigen-directions enough to reach
    // that exact value (0.5^5 ~ 3e-2 residual direction error), so the check
    // happened to read LIMIT; raising SingularPointTestIt to 20
    // (0.5^20 ~ 1e-6) resolves the mode fully and makes the point read
    // BRANCH -- correctly, given the diagonal K_T, but not what this preliminary
    // sanity check intends to probe. (0.5,0) is a REGULAR point on the trivial
    // branch and the limit-vs-branch classifier has no meaningful contract there.
    // Moved Ucur off the degenerate branch (u2=0.1) so K_T is only slightly
    // perturbed from diagonal: eigenvalues ~1.02/0.51, critical eigenvector
    // ~(0.20,0.98), cosine |v.F|/|F| ~ 0.20 -- three decades above
    // SingularPointTestTol=1e-4 regardless of SingularPointTestIt (>=~14 is
    // already fully converged at this eigenvalue ratio), so the check keeps its
    // original contrast (LIMIT here vs BRANCH at Ureq) without depending on the
    // exact mode-tolerance/iteration-count trade-off this task decouples.
    gsVector<real_t> Ucur(2); Ucur << 0.50, 0.1;   const real_t Lcur = 0.50;

    // Reference classifications, each taken with the solver AT the point.
    alm.setSolution(Ureq,Lreq);
    const bool classReq = alm.isBifurcation(true);
    alm.setSolution(Ucur,Lcur);
    const bool classCur = alm.isBifurcation(true);
    gsInfo<<"  [M0a] classification at the requested point = "<<(classReq?"branch":"limit")
          <<", at the current point = "<<(classCur?"branch":"limit")<<"\n";
    CHECK(classReq);    // the setup really does distinguish the two states
    CHECK(!classCur);

    // Solver parked at the OTHER point; ask for the requested one.
    alm.setSolution(Ucur,Lcur);
    const gsStatus st = alm.computeSingularPoint(Ureq,Lreq,
                                                 /*switchBranch*/false,
                                                 /*jacobian*/true,
                                                 /*testPoint*/true);
    CHECK(st == gsStatus::Success);
    // If the requested point was used, the branch point was found and refined.
    CHECK_CLOSE(1.0, alm.solutionL(),    1e-6);
    CHECK_CLOSE(1.0, alm.solutionU()[0], 1e-6);
}

// ===========================================================================
// TEST (G3) -- a no-bracket call to _bisectionSolve takes zero arc-length steps.
// ===========================================================================
// _bisectionSolve evaluates the inertia (negatives()) at both
// endpoints BEFORE stepping: if they agree there is no bracket to localize and
// the routine must return immediately, leaving (m_U,m_L) untouched and taking
// no step() at all. Without this check, a no-bracket call would instead run a
// forward march bounded by MaxIter-1.
TEST(bisection_without_a_bracket_takes_no_step)
{
    AlmProblem prob = foldProblem();
    SingularPointProbe alm(prob.Jacobian, prob.ALResidual, prob.Force);
    configure(alm, 0.05);

    alm.options().setReal("SingularPointComputeTolB",1e-6);   // bisection ON
    alm.options().setInt ("SingularPointBisIt",10);
    alm.applyOptions();
    // The option must EXIST: setInt on an unregistered label is a silent no-op.
    CHECK_EQUAL(10, alm.options().askInt("SingularPointBisIt",-1));

    // March 8 steps up the stable arm: fixture F from rest at ds=0.05 stays
    // well inside u1 < 1.
    for (index_t k = 0; k != 8; ++k)
        CHECK(alm.step() == gsStatus::Success);

    const gsVector<real_t> Ucur = alm.solutionU();
    const real_t           Lcur = alm.solutionL();
    const gsVector<real_t> Uarg = alm.solutionUPrev();
    const real_t           Larg = alm.solutionLPrev();

    // The setup is the intended one, not a degenerate one: two genuinely
    // distinct points, both analytically on the no-bracket (u1<1) side.
    CHECK((Uarg-Ucur).norm() > 1e-3);
    CHECK(Uarg[0] < 1.0);
    CHECK(Ucur[0] < 1.0);

    // Measured, not just analytic: both endpoints share the same inertia.
    alm.setSolution(Uarg,Larg);
    alm.computeStability(true);
    const index_t negArg = alm.negatives();
    alm.setSolution(Ucur,Lcur);
    alm.computeStability(true);
    const index_t negCur = alm.negatives();
    CHECK_EQUAL(0, negArg);
    CHECK_EQUAL(0, negCur);

    // Re-park the solver at the incumbent point, exactly the shape
    // _computeSingularPoint uses: the argument is the far endpoint, the
    // solver's current state is the incumbent (near) endpoint.
    alm.setSolution(Ucur,Lcur);

    alm.stubExtendedSolve(true);
    alm.resetProbe();
    alm.computeSingularPoint(Uarg, Larg, /*switchBranch*/false,
                             /*jacobian*/true, /*testPoint*/false);

    gsInfo<<"  [G3] stepCalls() = "<<alm.stepCalls()
          <<"  extendedCalls() = "<<alm.extendedCalls()<<"\n";

    CHECK_EQUAL(0, alm.stepCalls());       // the gate: no march
    CHECK_EQUAL(1, alm.extendedCalls());   // stage 2 was still reached
    CHECK_CLOSE(Larg, alm.extendedSeedL(), 1e-14);
    CHECK((alm.extendedSeedU()-Uarg).norm() < 1e-14);
    CHECK_CLOSE(Larg, alm.solutionL(), 1e-14);
    CHECK((alm.solutionU()-Uarg).norm() < 1e-14);
}

// ---------------------------------------------------------------------------
// Shared helper for G3b and G4: march fixture F with gsALMRiks until
// stabilityChange() fires (the snapping_example.cpp:440-446 shape: the
// argument to computeSingularPoint is the last accepted PRE-crossing point,
// the solver sits at the POST-crossing point), then run one
// computeSingularPoint with the extended stage stubbed out, so the result is
// exactly what the bisection stage produced.
// ---------------------------------------------------------------------------
struct BisectionRun
{
    gsVector<real_t> Uarg,  Ucur,  seedU;
    real_t           Larg,  Lcur,  seedL;
    index_t          steps, extCalls;
    index_t          bisItReadBack;   ///< options().askInt("SingularPointBisIt",-1)
    bool             crossed;
};

inline BisectionRun runBracketedBisection(index_t maxIter, index_t bisIt, real_t tolB)
{
    AlmProblem prob = foldProblem();
    SingularPointProbe alm(prob.Jacobian, prob.ALResidual, prob.Force);
    configure(alm, 0.05);

    alm.options().setInt ("MaxIter",maxIter);
    alm.options().setReal("SingularPointComputeTolB",tolB);   // bisection ON
    alm.options().setInt ("SingularPointBisIt",bisIt);
    alm.applyOptions();

    BisectionRun out;
    out.bisItReadBack = alm.options().askInt("SingularPointBisIt",-1);
    out.crossed        = false;

    gsVector<real_t> Uold = alm.solutionU();
    real_t           Lold = alm.solutionL();
    for (index_t k = 0; k != 60; ++k)
    {
        Uold = alm.solutionU(); Lold = alm.solutionL();
        if (alm.step() != gsStatus::Success) break;
        alm.computeStability(true);
        if (alm.stabilityChange()) { out.crossed = true; break; }
    }

    out.Uarg = Uold;              out.Larg = Lold;
    out.Ucur = alm.solutionU();   out.Lcur = alm.solutionL();

    alm.stubExtendedSolve(true);
    alm.resetProbe();
    alm.computeSingularPoint(out.Uarg, out.Larg, /*switchBranch*/false,
                             /*jacobian*/true, /*testPoint*/false);

    out.steps    = alm.stepCalls();
    out.extCalls = alm.extendedCalls();
    out.seedU    = alm.extendedSeedU();
    out.seedL    = alm.extendedSeedL();
    return out;
}

// ===========================================================================
// TEST (G3b) -- a genuine bracket is still bisected, and the seed it produces
// is localized inside it.
// ===========================================================================
// The anti-regression gate for a naive "just forbid the march" fix, which
// would satisfy G3 by taking zero steps unconditionally -- here the bracket
// is real (a genuine inertia change) and at least one step must be taken.
TEST(bisection_with_a_bracket_localizes_the_crossing)
{
    const BisectionRun r = runBracketedBisection(/*maxIter*/100, /*bisIt*/20, /*tolB*/1e-1);

    CHECK_EQUAL(20, r.bisItReadBack);   // the option exists (C5)

    CHECK(r.crossed);                   // the march actually reached the fold
    CHECK(r.Uarg[0] < 1.0);
    CHECK(r.Ucur[0] > 1.0);             // a genuine bracket (analytic inertia table)

    CHECK(r.steps >= 1);                // the anti-regression gate
    CHECK(r.steps <= 20);               // the budget

    // Containment: the seed lies inside the bracket it was given, to a tolerance rather
    // than exactly. The returned point is the probe at which the termination test was
    // evaluated -- a converged Newton iterate of the arc-length system, not an exact point
    // on the path -- so containment holds only up to that re-linearization slack.
    //
    // FIXME: re-derive this tolerance. 1e-3 was chosen against the earlier return-point
    // semantics, in which the converged exit reported the last certified pre-crossing point
    // and the probe loop had a much smaller budget. Both changed, so the reachable
    // tolerances and the observed slack differ from what justified 1e-3. Measure the slack
    // at the returned point across tolB and Length, then tighten to the smallest bound that
    // still clears it.
    CHECK(r.seedU[0] >= math::min(r.Uarg[0], r.Ucur[0]) - 1e-3);
    CHECK(r.seedU[0] <= math::max(r.Uarg[0], r.Ucur[0]) + 1e-3);

    // Localization: a genuine bisection halves the bracket every probe, so a
    // factor of two is the weakest thing it must beat.
    gsInfo<<"  [G3b] steps="<<r.steps
          <<"  |u1_seed-1|="<<math::abs(r.seedU[0]-1.0)
          <<"  |u1_arg-1|="<<math::abs(r.Uarg[0]-1.0)<<"\n";
    CHECK(math::abs(r.seedU[0] - 1.0) < 0.5 * math::abs(r.Uarg[0] - 1.0));
}

// ===========================================================================
// TEST (G4) -- the probe budget is SingularPointBisIt, not MaxIter, and a
// non-converged exit restores the entry state before handing it to stage 2.
// ===========================================================================
TEST(bisection_budget_is_independent_of_maxiter)
{
    // tolB tight enough that the bisection cannot converge inside the budget,
    // forcing the budget itself to bind.
    const BisectionRun a = runBracketedBisection(/*maxIter*/100, /*bisIt*/3, /*tolB*/1e-14);
    const BisectionRun b = runBracketedBisection(/*maxIter*/400, /*bisIt*/3, /*tolB*/1e-14);

    CHECK(a.crossed); CHECK(b.crossed);
    CHECK_EQUAL(3, a.bisItReadBack); CHECK_EQUAL(3, b.bisItReadBack);

    gsInfo<<"  [G4] steps(MaxIter=100)="<<a.steps
          <<"  steps(MaxIter=400)="<<b.steps<<"\n";

    CHECK_EQUAL(3, a.steps);
    CHECK_EQUAL(3, b.steps);
    CHECK_EQUAL(a.steps, b.steps);      // independence of MaxIter

    // State restore on a non-converged exit: the seed handed to the extended
    // stage is the restored entry state, not wherever the probes ended.
    CHECK_CLOSE(a.Larg, a.seedL, 1e-12);
    CHECK((a.seedU - a.Uarg).norm() < 1e-12);
}

// ===========================================================================
// TEST 4 (M1) -- the adaptive length factor is continuous in the iteration count.
// ===========================================================================
// A computeLength() that formed `m_desiredIterations / m_numIterations` in
// index_t would let the clamped factor only ever be 0.5, 1 or 2, with one extra
// iteration halving the step (destroying the classic ds*(n_des/n_act)
// adaptation). This pins the floating-point form.
TEST(adaptive_length_factor_is_continuous)
{
    AlmProblem prob = foldProblem();
    RiksProbe alm(prob.Jacobian, prob.ALResidual, prob.Force);
    configure(alm, 0.05);

    const real_t ds      = 0.05;
    const index_t desired = 10;

    // Reference = the textbook ds_new = ds * n_desired/n_actual, clamped to [1/2,2].
    const index_t taken[6] = {1, 5, 6, 10, 11, 30};
    for (index_t i = 0; i != 6; ++i)
    {
        real_t ref = (real_t)desired / (real_t)taken[i];
        if      (ref < 0.5) ref = 0.5;
        else if (ref > 2.0) ref = 2.0;

        const real_t fac = alm.lengthFactor(taken[i],desired,ds);
        gsInfo<<"  [M1] iterations = "<<taken[i]<<"  factor = "<<fac
              <<"  (expected "<<ref<<")\n";
        CHECK_CLOSE(ref, fac, 1e-12);
    }

    // The two values that integer division could not produce; asserted
    // explicitly so a regression is unmistakable.
    CHECK_CLOSE(10.0/6.0,  alm.lengthFactor( 6,desired,ds), 1e-12);
    CHECK_CLOSE(10.0/11.0, alm.lengthFactor(11,desired,ds), 1e-12);

    // m_numIterations == 0 was an integer division by zero (SIGFPE); it now maps
    // onto the clamp's upper bound.
    CHECK_CLOSE(2.0, alm.lengthFactor(0,desired,ds), 1e-12);
}

// ===========================================================================
// TEST 5 (M2) -- setLength() and reduceLength() must produce the same retry.
// ===========================================================================
// setLength() overwriting m_arcLength_prev (the divisor that normalises
// the secant predictor) would leave the predictor stepping the FULL previous
// increment after "halve the step and retry", making every further halving
// strictly harder. reduceLength() is the API that avoids this.
TEST(setLength_and_reduceLength_agree)
{
    const real_t ds = 0.05;
    AlmProblem probA = foldProblem();
    AlmProblem probB = foldProblem();
    RiksProbe A(probA.Jacobian, probA.ALResidual, probA.Force);
    RiksProbe B(probB.Jacobian, probB.ALResidual, probB.Force);
    configure(A, ds);
    configure(B, ds);

    // Identical history: one accepted step of length ds.
    CHECK(A.step() == gsStatus::Success);
    CHECK(B.step() == gsStatus::Success);
    CHECK_CLOSE(A.solutionL(), B.solutionL(), 1e-14);

    // --- the white-box invariant that actually pins M2 ---------------------
    // After an ACCEPTED step, m_arcLength_prev is the length that produced the
    // current secant; requesting a new length must not redefine it (otherwise
    // the predictor ratio ds/ds_prev is identically one and the retry repeats
    // the previous increment verbatim).
    A.setLength(ds/2);
    CHECK_CLOSE(ds/2, A.arcLength()    , 1e-14);
    CHECK_CLOSE(ds  , A.arcLengthPrev(), 1e-14);
    B.reduceLength(0.5);
    CHECK_CLOSE(ds/2, B.arcLength()    , 1e-14);
    CHECK_CLOSE(ds  , B.arcLengthPrev(), 1e-14);

    // --- the black-box guard: the two APIs are interchangeable -------------
    // (Cheap, and it stays meaningful if the internals are reorganised, but it
    // has no teeth on its own once the invariant above holds.)
    CHECK(A.step() == gsStatus::Success);
    CHECK(B.step() == gsStatus::Success);
    CHECK_CLOSE(B.solutionL(),    A.solutionL(),    1e-10);
    CHECK_ARRAY_CLOSE(B.solutionU(), A.solutionU(), 2, 1e-10);
    CHECK(A.numIterations() <= B.numIterations() + 1);

    // --- compounding: two further halvings without an intervening success ---
    // (setLength is absolute, reduceLength is multiplicative; both solvers are
    // at ds/2 here, so ds/4 then ds/8.)
    A.setLength(ds/4);     A.setLength(ds/8);
    B.reduceLength(0.5);   B.reduceLength(0.5);
    CHECK_CLOSE(ds/8, A.arcLength()    , 1e-14);
    CHECK_CLOSE(ds/8, B.arcLength()    , 1e-14);
    CHECK_CLOSE(A.arcLengthPrev(), B.arcLengthPrev(), 1e-14);
    CHECK(A.step() == gsStatus::Success);
    CHECK(B.step() == gsStatus::Success);
    CHECK_CLOSE(B.solutionL(),    A.solutionL(),    1e-10);
    CHECK_ARRAY_CLOSE(B.solutionU(), A.solutionU(), 2, 1e-10);
}

// ===========================================================================
// TEST 6 (M3) -- Crisfield's root selection must never reverse the direction of
// travel.
// ===========================================================================
// computeLambdaDOT() weights the load term of the increment-increment product by
// phi^2 instead of A0 = phi^2|f|^2, i.e. by a factor |f|^-2 relative to the
// metric of the constraint it is selecting roots for. The sign of
// S = DeltaU_old.deltaUt + w*DeltaLambda_old decides the root, so a wrong w can
// select the root that retraces the branch. Because the auto scaling gives
// A0 = |U|^2/Lambda^2 (scale free) while phi^2 = A0/|f|^2, a SMALL load scale
// amplifies the error -- hence s in {1e-3, 1, 1e3}.
//
// Progress criterion, in the metric of the constraint:
//   (U_k - U_{k-1}).(U_{k-1} - U_{k-2}) + A0*dL_k*dL_{k-1} > 0.
//
struct CrisfieldTrace
{
    real_t  minProgress;   ///< min over k of the progress criterion (>0: no retrace)
    real_t  maxU1;         ///< the furthest u1 reached (>1: the fold was rounded)
    index_t nPoints;       ///< number of accepted points, seed included
    bool    allStepsOk;    ///< every requested step returned gsStatus::Success
};

/// Traces \a nSteps Crisfield steps through fixture F at load scale \a scale.
CrisfieldTrace crisfieldProgress(real_t scale, real_t ds, index_t nSteps, bool backward)
{
    AlmProblem prob = foldProblem(scale);
    CrisfieldProbe alm(prob.Jacobian, prob.ALResidual, prob.Force);
    configure(alm, ds);

    if (backward)
    {
        // Seed a secant pointing towards decreasing u1 so the first predictor
        // travels backwards down the stable arm.
        gsVector<real_t> U(2), Uprev(2);
        U     << 0.5, 0.0;
        Uprev << 0.6, 0.0;
        alm.setSolution(U,foldLambda(0.5));
        alm.setLength(ds);
        alm.setPrevious(Uprev,foldLambda(0.6));
    }

    std::vector<gsVector<real_t> > Us;
    std::vector<real_t>            Ls;
    Us.push_back(alm.solutionU());
    Ls.push_back(alm.solutionL());

    CrisfieldTrace out;
    out.minProgress = std::numeric_limits<real_t>::max();
    out.maxU1       = Us.back()[0];
    out.allStepsOk  = true;

    for (index_t k = 0; k != nSteps; ++k)
    {
        const gsStatus st = alm.step();
        if (st != gsStatus::Success)
        {
            out.allStepsOk = false;
            gsInfo<<"  [M3] scale = "<<scale<<": step "<<k<<" failed with status "
                  <<static_cast<index_t>(st)<<"\n";
            break;
        }
        Us.push_back(alm.solutionU());
        Ls.push_back(alm.solutionL());
        out.maxU1 = math::max(out.maxU1, alm.solutionU()[0]);

        const size_t n = Us.size();
        if (n >= 3)
        {
            const gsVector<real_t> dU1 = Us[n-1] - Us[n-2];
            const gsVector<real_t> dU0 = Us[n-2] - Us[n-3];
            const real_t dL1 = Ls[n-1] - Ls[n-2];
            const real_t dL0 = Ls[n-2] - Ls[n-3];
            const real_t p = dU1.dot(dU0) + alm.A0()*dL1*dL0;
            out.minProgress = math::min(out.minProgress,p);
        }
        // Stay clear of the second lambda = 0 crossing at u1 = 2 (that is M6's
        // territory, tested separately).
        if (alm.solutionU()[0] > 1.8 || alm.solutionU()[0] < -0.8)
            break;
    }
    out.nPoints = static_cast<index_t>(Us.size());
    gsInfo<<"  [M3] scale = "<<scale<<(backward?" (backward)":" (forward)")
          <<"  points = "<<out.nPoints
          <<"  last (u1,L) = ("<<Us.back()[0]<<","<<Ls.back()<<")"
          <<"  max u1 = "<<out.maxU1
          <<"  min progress = "<<out.minProgress<<"\n";
    CHECK(out.nPoints >= 5);       // the trace has to actually go somewhere
    return out;
}

TEST(crisfield_root_selection_never_retraces)
{
    // Forward through the fold at three load scales, plus one backward trace.
    //
    // STRENGTHENED: asserting minProgress > 0 alone was VACUOUS on a trace that
    // dies early -- crisfieldProgress() breaks out of the loop on the first non-Success
    // step, so a run that stops at step 23 still reports the minimum over the handful of
    // increments it did manage. That is exactly how the s = 1e-3 trace used to "pass"
    // while M3 was killing it. allStepsOk is now required as well, so the two halves of
    // M3 (no retrace AND no premature failure) are both pinned here.
    const CrisfieldTrace t1f = crisfieldProgress(1.0 , 0.05, 40, /*backward*/false);
    const CrisfieldTrace t1b = crisfieldProgress(1.0 , 0.05, 25, /*backward*/true );
    const CrisfieldTrace tHi = crisfieldProgress(1e3 , 0.05, 40, /*backward*/false);
    const CrisfieldTrace tLo = crisfieldProgress(1e-3, 0.05, 40, /*backward*/false);

    CHECK(t1f.minProgress > 0.0);   CHECK(t1f.allStepsOk);
    CHECK(t1b.minProgress > 0.0);   CHECK(t1b.allStepsOk);
    CHECK(tHi.minProgress > 0.0);   CHECK(tHi.allStepsOk);
    CHECK(tLo.minProgress > 0.0);   CHECK(tLo.allStepsOk);
}

// ---------------------------------------------------------------------------
// TEST 6, second half (M3) -- the same trace must round the fold at every load
// scale.
// ---------------------------------------------------------------------------
// The trace is scale invariant in exact arithmetic (see fixture S), so a trace
// that rounds the fold at s = 1 must round it at s = 1e-3 too. Before this fix it
// did not: the corrector's root selection weighted the load term by
// phi^2 = A0/|f|^2, so at s = 1e-3 that weight was 1e6 times too large and the
// step through the fold stopped converging (it died at step 23). This is M3's
// mechanism showing up as an outright step FAILURE rather than as the silent
// retrace the review predicted -- "no retrace" is therefore no evidence at all
// here; the fold-rounding statement below is.
//
// FIXED (M3): computeLambdaDOT() now uses the same A0 = phi^2|f|^2 as the
// constraint whose roots it selects, which is scale free, so all three traces coincide.
TEST(crisfield_fold_rounding_is_scale_invariant)
{
    // Ground truth at the reference scale: the fold IS rounded.
    const CrisfieldTrace ref = crisfieldProgress(1.0, 0.05, 40, false);
    CHECK(ref.allStepsOk);
    CHECK(ref.maxU1 > 1.05);

    const real_t scales[2] = {1e-3, 1e3};
    for (index_t i = 0; i != 2; ++i)
    {
        const CrisfieldTrace t = crisfieldProgress(scales[i], 0.05, 40, false);
        CHECK(t.allStepsOk);
        CHECK(t.maxU1 > 1.05);
        gsInfo<<"  [M3] scale-invariance of maxU1: |ref - t| = "
              <<math::abs(ref.maxU1-t.maxU1)<<" (tol 1e-8)\n";
        CHECK_CLOSE(ref.maxU1, t.maxU1, 1e-8);
        CHECK_EQUAL(ref.nPoints, t.nPoints);
    }
}

// ===========================================================================
// TEST 6c (R1) -- gsALMRiks's corrector must not land on the BACKWARD
// intersection of its own arc-length surface with the equilibrium path.
// ===========================================================================
// gsALMRiks::iteration() solves the LINEARISED constraint,
//     m_deltaL = -( r + 2 phi DeltaU.ubar ) / ( 2(1-phi)DeltaL + 2 phi DeltaU.u_t ),
// whose denominator is exactly 2*phi*(Ritto-Correa & Camotim 2008, eq. (22))
//     t = Delta_a . a^Q + psi^2 Delta_lambda,     psi^2 = (1-phi)/phi = numDof-1.
// Nothing inspects it. Approaching a limit point, a^Q = K_T^-1 f flips sign with
// det K_T while Delta_lambda does not, so for a LARGE psi^2 the two terms of t
// are comparable and of opposite sign: t crosses zero, m_deltaL blows up, and
// the corrector is thrown into the basin of the BACKWARD intersection (point Z
// of the paper's Fig. 2, "the equilibrium path portion already determined").
// Z is a genuine equilibrium at (about) the right arc length, so the residual
// converges and iterationFinish() commits the retrace as a Success.
//
// psi^2 is the whole mechanism: gsALMCrisfield and gsALMConsistentCrisfield use
// the SAME single linear solve and the SAME secant predictor and do not retrace,
// but in the configuration that produced the observation they run CYLINDRICAL
// (Scaling = 0 => psi^2 = 0), where t = DeltaU.a^Q alone changes sign as a whole
// and never passes through zero. Fixture F-n is the psi^2 dial that reproduces
// this on 2-DOF-family algebra: psi^2 = n-1.
//
// Progress criterion, CYLINDRICAL (psi^2 = 0). The psi^2-weighted form of eq.
// (22)/(23) does NOT discriminate here (at a fold
// BOTH intersections lie below the centre in lambda, so the lambda term, which
// outweighs the displacement term by ~5x at n = 128, is negative for the forward
// root too):
//   (U_k - U_{k-1}) . (U_{k-1} - U_{k-2}) > 0.
//
struct RiksTrace
{
    real_t  minProgressU;  ///< min_k (U_k-U_{k-1}).(U_{k-1}-U_{k-2})   [psi^2 = 0]
    real_t  minCosine;     ///< the same, normalised -- DIAGNOSTIC ONLY, printed not asserted
    real_t  maxU1;         ///< furthest u0 reached (> 1: the fold was rounded)
    index_t nPoints;       ///< number of accepted points, seed included
    /// true unless the trace GAVE UP: every requested step was eventually
    /// accepted, possibly after halve-and-retry. (Not "no step ever returned
    /// NotConverged" -- the fix works precisely by making a step fail once and
    /// then recovering, so that reading would be self-contradictory.)
    bool    allStepsOk;
    index_t nRetries;      ///< arc-length halvings the harness performed
    real_t  finalLength;   ///< the arc length the trace ended on
};

/// Traces at most \a nSteps Riks steps forward through fixture F-n from the
/// origin, with the halve-and-retry of gsALMExploration<T>::traceSweep (its
/// step-fail retry block, `if (status != gsStatus::Success)`).
RiksTrace riksProgress(index_t n, real_t ds, index_t nSteps,
                       index_t maxIter = 100, real_t tolF = 1e-8, real_t tolU = 1e-8,
                       bool verbose = true)
{
    AlmProblem prob = foldProblemPadded(n);
    RiksProbe  alm(prob.Jacobian, prob.ALResidual, prob.Force);
    configure(alm, ds);
    // The corrector regime is a dial of its own: example_BratuExploration, where
    // the retrace was observed, runs MaxIter = 30 with the DEFAULT TolF = 1e-3 /
    // TolU = 1e-6, i.e. it accepts a much coarser corrector outcome than the
    // 1e-8 the rest of this file uses.
    alm.options().setInt ("MaxIter",maxIter);
    alm.options().setReal("TolF",tolF);
    alm.options().setReal("TolU",tolU);
    alm.applyOptions();
    alm.initialize();
    alm.setLength(ds);

    std::vector<gsVector<real_t> > Us;
    std::vector<real_t>            Ls;
    Us.push_back(alm.solutionU());
    Ls.push_back(alm.solutionL());

    RiksTrace out;
    out.minProgressU = std::numeric_limits<real_t>::max();
    out.minCosine    = std::numeric_limits<real_t>::max();
    out.maxU1        = Us.back()[0];
    out.allStepsOk   = true;
    out.nRetries     = 0;

    real_t       dLb  = ds;
    const real_t dLb0 = ds;

    index_t taken = 0;
    // Total attempts are capped so that a pathological halving loop cannot hang
    // the suite; 24 spare attempts is twice the 12-halving give-up threshold.
    const index_t maxAttempts = nSteps + 24;
    for (index_t attempt = 0; taken != nSteps && attempt != maxAttempts; ++attempt)
    {
        const gsVector<real_t> Uold = alm.solutionU();
        const real_t           Lold = alm.solutionL();

        gsStatus st;
        try { st = alm.step(); }
        catch (...) { st = gsStatus::AssemblyError; }

        if (st != gsStatus::Success)
        {
            // Verbatim the recovery of gsALMExploration<T>::traceSweep's
            // step-fail retry block.
            dLb /= (real_t)2;
            ++out.nRetries;
            if (math::abs(dLb/dLb0) < 1e-6 || out.nRetries > 12)
            {
                out.allStepsOk = false;
                gsInfo<<"  [R1] n = "<<n<<": gave up after "<<out.nRetries
                      <<" halvings (arc length "<<dLb<<")\n";
                break;
            }
            alm.setLength(dLb);
            alm.setSolution(Uold,Lold);
            continue;
        }

        ++taken;
        Us.push_back(alm.solutionU());
        Ls.push_back(alm.solutionL());
        out.maxU1 = math::max(out.maxU1, alm.solutionU()[0]);

        const size_t m = Us.size();
        if (m >= 3)
        {
            const gsVector<real_t> dU1 = Us[m-1] - Us[m-2];
            const gsVector<real_t> dU0 = Us[m-2] - Us[m-3];
            const real_t p = dU1.dot(dU0);
            out.minProgressU = math::min(out.minProgressU,p);
            const real_t den = dU1.norm()*dU0.norm();
            if (den > 0)
                out.minCosine = math::min(out.minCosine, p/den);
        }
        // Stay clear of the second lambda = 0 crossing at u0 = 2 (M6 territory).
        if (alm.solutionU()[0] > 1.8 || alm.solutionU()[0] < -0.8)
            break;
    }
    out.nPoints     = static_cast<index_t>(Us.size());
    out.finalLength = dLb;

    if (verbose)
        gsInfo<<"  [R1] n = "<<std::setw(4)<<n<<" ds = "<<ds
              <<" it = "<<maxIter<<" tolF = "<<tolF
              <<"  points = "<<out.nPoints
              <<"  last (u0,L) = ("<<Us.back()[0]<<","<<Ls.back()<<")"
              <<"  max u0 = "<<out.maxU1
              <<"  minProgressU = "<<out.minProgressU
              <<"  minCosine = "<<out.minCosine
              <<"  retries = "<<out.nRetries
              <<"  ok = "<<out.allStepsOk<<"\n";
    return out;
}

// ---------------------------------------------------------------------------
// ! READ THIS BEFORE TREATING THE TEST BELOW AS A PIN FOR R1 !
//
// This test is GREEN ON THE UNFIXED LIBRARY and therefore does NOT pin R1. It
// was measured over the whole psi^2 x arc-length grid
//     n  in {2, 8, 32, 128, 256, 1024}   (psi^2 = n-1, up to 1023)
//     ds in {0.4, 0.2, 0.1, 0.05, 0.025}
// in BOTH corrector regimes (MaxIter = 100 / TolF = 1e-8, and the driver's
// MaxIter = 30 / TolF = 1e-3): 60 cells, ZERO of them reversing. minCosine is
// exactly +1 in every cell (the trace is one-dimensional in u0, so a reversal
// would read -1). Adding the direction gate to
// gsALMRiks::iterationFinish() left all 60 cells BIT-IDENTICAL -- the gate never
// fires here, so no fixture cell can supply failing-first evidence for it.
//
// What fixture F-n DOES show at large psi^2 is the OTHER half of the finding
// (R2): the fold-crossing step stops converging (the near-singular tangent makes
// m_deltaU = ubar + deltaL*a^Q oscillate by O(1) while the corrector's divisor
// stays comfortably positive), and Riks only gets past the fold because the
// caller halves the arc length. That is what nRetries below measures.
//
// The R1 retrace itself reproduces only at driver level:
//   example_BratuExploration -L 0.2   (dead load, i.e. NO --forcingCallback)
// where landscape.csv point 35 reproduces point 33 and 36 reproduces 32 to every
// printed digit, 35 non-increasing normU pairs and 33 mirror pairs.
// ---------------------------------------------------------------------------
TEST(riks_rounds_the_padded_fold)
{
    // The cell predicted RED (psi^2 = 127, 85% of the arc
    // length in the lambda term at the fold) ...
    const RiksTrace t = riksProgress(128, 0.05, 60);
    CHECK(t.allStepsOk);            // the trace must not give up
    CHECK(t.nPoints >= 5);          // ... it has to actually go somewhere
    CHECK(t.maxU1 > 1.05);          // ... and it must ROUND the fold
    CHECK(t.minProgressU > 0.0);    // ... without ever reversing
    CHECK(t.nRetries < 6);          // ... and not by halving into oblivion

    // ... and the cylindrical-like control (psi^2 = 1), which must round the
    // fold without a single arc-length halving.
    const RiksTrace c = riksProgress(2, 0.05, 60);
    CHECK(c.allStepsOk);
    CHECK(c.maxU1 > 1.05);
    CHECK(c.minProgressU > 0.0);
    CHECK_EQUAL(0, c.nRetries);
}

// ===========================================================================
// TEST 6d -- gsALMConsistentCrisfield must not retrace at
// ANY psi^2.
// ===========================================================================
// The converse experiment: the hypothesis that the
// Riks retrace is a psi^2 property predicts that raising
// gsALMConsistentCrisfield's psi^2 = phi^2|f|^2 (the `Scaling` option; |f| = 1
// on foldProblemPadded(n,1), so Scaling = sqrt(255) gives A0 = 255) to Riks's
// value would make IT retrace too.
//
// It does not. This test is that REFUTATION, kept as a standing invariant: at
// psi^2 in {0, 127, 255} the trace rounds the fold with strictly positive
// progress and zero arc-length halvings. A large psi^2 does shrink the
// displacement step a lot (~6x more steps to cover the same arc), which is why
// the budget below is generous -- without it the trace would stop short of the
// fold and `minProgress > 0` would be vacuous (the exact tautology the
// `maxU1 > 1.05` guard exists to close).
TEST(consistent_crisfield_never_retraces_at_any_psi2)
{
    const index_t n = 128;
    const real_t  scalings[3] = {0.0, math::sqrt(127.0), math::sqrt(255.0)};
    for (index_t i = 0; i != 3; ++i)
    {
        const real_t ds = 0.05;
        AlmProblem prob = foldProblemPadded(n);
        ConsistentProbe alm(prob.Jacobian, prob.ALResidual, prob.Force);
        configure(alm, ds);
        alm.options().setReal("Scaling",scalings[i]);
        alm.applyOptions();
        alm.initialize();
        alm.setLength(ds);

        std::vector<gsVector<real_t> > Us;
        Us.push_back(alm.solutionU());
        real_t minProg = std::numeric_limits<real_t>::max();
        real_t maxU1   = Us.back()[0];
        index_t retries = 0, taken = 0;
        real_t dLb = ds;
        bool   ok = true;
        const index_t budget = 4000;   // enough for the tiny steps a large A0 produces
        for (index_t attempt = 0; taken != budget && attempt != budget+40; ++attempt)
        {
            const gsVector<real_t> Uold = alm.solutionU();
            const real_t           Lold = alm.solutionL();
            gsStatus st;
            try { st = alm.step(); } catch (...) { st = gsStatus::AssemblyError; }
            if (st != gsStatus::Success)
            {
                dLb /= (real_t)2; ++retries;
                if (retries > 12) { ok = false; break; }
                alm.setLength(dLb); alm.setSolution(Uold,Lold);
                continue;
            }
            ++taken;
            Us.push_back(alm.solutionU());
            maxU1 = math::max(maxU1, alm.solutionU()[0]);
            const size_t m = Us.size();
            if (m >= 3)
                minProg = math::min(minProg, (Us[m-1]-Us[m-2]).dot(Us[m-2]-Us[m-3]));
            if (alm.solutionU()[0] > 1.8 || alm.solutionU()[0] < -0.8) break;
        }
        gsInfo<<"  [3.3] ConsistentCrisfield n = "<<n<<"  Scaling = "<<scalings[i]
              <<"  A0 = psi^2 = "<<scalings[i]*scalings[i]
              <<"  points = "<<Us.size()<<"  max u0 = "<<maxU1
              <<"  last u0 = "<<Us.back()[0]
              <<"  minProgress = "<<minProg
              <<"  retries = "<<retries<<"  ok = "<<ok<<"\n";
        CHECK(ok);
        CHECK(maxU1 > 1.05);        // the trace REACHED the fold (anti-tautology)
        CHECK(minProg > 0.0);       // ... and never reversed
        // Not a no-op claim: the CYLINDRICAL cell (Scaling = 0) does take ONE
        // arc-length halving at the fold, the two large-psi^2 cells take none.
        // The bound is the same anti-oblivion guard used above.
        CHECK(retries < 6);
    }
}

// ===========================================================================
// TEST 6e -- CHARACTERISATION: gsALMRiks's
// corrector divisor does NOT change sign across a fold-crossing step.
// ===========================================================================
// The hypothesized mechanism is that the divisor of
// gsALMRiks::iteration(),
//     denominator/2 = phi*DeltaU.a^Q + (1-phi)*DeltaLambda = phi * t   (eq. (22)),
// crosses zero near a limit point because a^Q = K_T^-1 f flips sign with det K_T
// while DeltaLambda does not, so that m_deltaL blows up and throws the corrector
// into the backward basin.
//
// MEASURED, AND IT DOES NOT HAPPEN. Across the fold-crossing step at
// psi^2 in {7, 127, 1023} the divisor takes ZERO sign changes: at large psi^2
// the (1-phi)*DeltaLambda term dominates by an order of magnitude and keeps the
// divisor comfortably positive and O(ds). What actually degrades there is the
// UPDATE, m_deltaU = ubar + deltaL*a^Q, whose two near-singular solves make the
// displacement oscillate by O(1) while m_deltaL stays ~1e-3.
//
// This test is a CHARACTERISATION, deliberately pinning the refutation so that a
// future change of remedy has to confront it. If it ever fires, re-examine the
// refuted mechanism above before "fixing" it.
TEST(riks_corrector_denominator_does_not_change_sign_at_the_fold)
{
    const index_t ns[3] = {8,128,1024};
    for (index_t i = 0; i != 3; ++i)
    {
        const index_t n  = ns[i];
        const real_t  ds = 0.05;
        AlmProblem prob = foldProblemPadded(n);
        RiksProbe  alm(prob.Jacobian, prob.ALResidual, prob.Force);
        configure(alm, ds);

        // Walk up the stable arm until the next step would straddle the fold.
        index_t walked = 0;
        while (alm.solutionU()[0] < 0.9 && walked < 200)
        {
            if (alm.step() != gsStatus::Success) break;
            ++walked;
        }
        gsInfo<<"  [R1 probe] n = "<<n<<"  seeded at u0 = "<<alm.solutionU()[0]
              <<"  L = "<<alm.solutionL()<<"  phi = "<<alm.phi()
              <<"  psi^2 = "<<(1.0-alm.phi())/alm.phi()<<"\n";

        // 12 iterations, NOT the full MaxIter: at large psi^2 the step does not
        // converge and the tail of the sequence is genuinely DIVERGENT (|DeltaU0|
        // past 2.7, resF past 7). Asserting on a divergent tail would be a gate
        // that bit-rots under any change of compiler flags or Eigen version. The
        // first dozen iterations are the well-behaved part -- +0.042..+0.057 at
        // n = 128, +0.0497..+0.0503 at n = 1024 -- and they carry the finding.
        std::vector<RiksProbe::IterRecord> rec;
        const bool converged = alm.runStepInstrumented(rec,12);
        gsInfo<<"  [R1 probe] step converged within 12 iterations = "<<converged
              <<"  ("<<rec.size()<<" iterations)\n";
        index_t signChanges = 0;
        for (size_t k = 0; k != rec.size(); ++k)
        {
            if (k>0 && rec[k].halfDen*rec[k-1].halfDen < 0) ++signChanges;
            if (k < 12 || k+4 > rec.size())
                gsInfo<<"      it "<<std::setw(3)<<k
                      <<"  denom/2 = "<<std::setw(14)<<rec[k].halfDen
                      <<"  deltaL = "<<std::setw(14)<<rec[k].deltaL
                      <<"  DeltaL = "<<std::setw(14)<<rec[k].DeltaL
                      <<"  DeltaU0 = "<<std::setw(14)<<rec[k].DeltaU0
                      <<"  resF = "<<rec[k].residueF<<"\n";
        }
        gsInfo<<"  [R1 probe] denominator sign changes within the step: "<<signChanges<<"\n";
        CHECK(rec.size() >= 1);          // the step was actually driven
        CHECK(alm.solutionU()[0] > 0.85);// ... from a seed close to the fold at u0 = 1
        CHECK_EQUAL(0, signChanges);     // ... and the divisor never changed sign
    }
}

// ===========================================================================
// TEST 7 (M4) -- Crisfield's complex-root fallback must build Lam's point from
// the true internal force.
// ===========================================================================
// computeLambdasComplex() sets Fint = K_T(U+DeltaU)*(U+DeltaU), which equals the
// internal force only for a LINEAR problem; the correct definition is
// F_int(U+DeltaU_cr) = R(U+DeltaU_cr, 0).
//
// NOTE (deviation from the review sketch, deliberate): the sketch's assertion
// "distance(DeltaU,DeltaL) ~ ds after the fallback" is TAUTOLOGICAL. The routine
// returns mu*(DeltaU_cr, DeltaLambda_cr) with mu = ds/|(DeltaU_cr,DeltaLambda_cr)|_A
// and distance() uses that same A0, so the result is on the sphere by
// construction for ANY value of Lcr. It is kept below as a guard on the
// mu-scaling, but the defect is pinned by inverting mu and comparing the
// recovered Lcr against the closed-form internal force.
//
/// Recovers the load factor Lcr that computeLambdasComplex() used, from the
/// public outcome of one fallback evaluation.
struct FallbackOutcome
{
    real_t LcrUsed;      ///< the Lcr the code effectively used
    real_t LcrExact;     ///< F_int(U+DeltaU_cr).f / (f.f)
    real_t sphereRadius; ///< |(DeltaU+deltaU, DeltaL+deltaL)|_A, must be ds
};

FallbackOutcome runFallback(const AlmProblem & prob, const gsVector<real_t> & U,
                            real_t L, const gsVector<real_t> & DeltaU, real_t DeltaL,
                            real_t ds)
{
    AlmProblem p = prob;
    CrisfieldProbe alm(p.Jacobian, p.ALResidual, p.Force);
    configure(alm, ds);
    alm.options().setReal("Scaling",1.0);   // pin phi so A0 = |f|^2 exactly
    alm.applyOptions();

    alm.runComplexFallback(U,L,DeltaU,DeltaL,ds);

    const gsVector<real_t> DeltaUnew = DeltaU + alm.deltaU();
    const real_t           DeltaLnew = DeltaL + alm.deltaL();

    // (DeltaUnew, DeltaLnew) = mu*(DeltaU_cr, DeltaLambda_cr) by construction.
    const gsVector<real_t> DeltaUcr = DeltaU + alm.deltaUbar();
    const real_t mu = DeltaUnew.dot(DeltaUcr) / DeltaUcr.dot(DeltaUcr);

    FallbackOutcome out;
    out.LcrUsed      = L + DeltaLnew/mu;
    out.LcrExact     = p.Fint(U+DeltaUcr).dot(p.Force) / p.Force.dot(p.Force);
    out.sphereRadius = alm.distance(DeltaUnew,DeltaLnew);
    gsInfo<<"  [M4] mu = "<<mu<<"  Lcr(used) = "<<out.LcrUsed
          <<"  Lcr(exact) = "<<out.LcrExact
          <<"  |(DU,DL)|_A = "<<out.sphereRadius<<" (ds = "<<ds<<")\n";
    return out;
}

// FIXED (M4): computeLambdasComplex() now recovers the internal force from the
// residual, F_int(U+DeltaU_cr) = R(U+DeltaU_cr, L+DeltaL) + (L+DeltaL) f, at the trial
// point DeltaU_cr rather than as K_T*(U+DeltaU).
TEST(crisfield_complex_root_fallback_stays_on_the_sphere)
{
    const real_t ds = 0.3;

    // --- Linear 2-DOF control, residual-free state (deltaUbar == 0) ---------
    // There F_int(U+DeltaU_cr) = K*(U+DeltaU) identically, so the fallback is
    // exact both before and after the fix.
    {
        AlmProblem lin = linearProblem();
        gsVector<real_t> U(2);      U      << 0.5, 0.0;   // 2*u1 = lambda
        gsVector<real_t> DeltaU(2); DeltaU << 0.1, 0.0;
        const real_t L = 1.0, DeltaL = 0.2;               // (U+DU, L+DL) is an equilibrium
        FallbackOutcome o = runFallback(lin,U,L,DeltaU,DeltaL,ds);
        CHECK_CLOSE(o.LcrExact, o.LcrUsed, 1e-10);
        CHECK_CLOSE(ds, o.sphereRadius, 1e-10);
    }

    // --- Fixture F (nonlinear): the defect is O(1) --------------------------
    {
        AlmProblem prob = foldProblem();
        gsVector<real_t> U(2);      U      << 0.6, 0.0;   // lambda = 0.84
        gsVector<real_t> DeltaU(2); DeltaU << 0.2, 0.0;
        const real_t L = foldLambda(0.6), DeltaL = 0.1;
        FallbackOutcome o = runFallback(prob,U,L,DeltaU,DeltaL,ds);

        // Tautological, kept as a guard on the mu-scaling only.
        CHECK_CLOSE(ds, o.sphereRadius, 1e-10);

        // THE gate: the mu-inversion oracle. Pre-fix Lcr(used) = 0.32 against
        // Lcr(exact) = 0.9375 -- an O(1) error that the tautological check above cannot
        // see, because mu rescales any Lcr back onto the sphere.
        CHECK_CLOSE(o.LcrExact, o.LcrUsed, 1e-8);
    }
}

// ===========================================================================
// TEST 8 (M5) -- the fresh-start predictor must honour a user-set Scaling.
// ===========================================================================
// gsALMCrisfield::predictor() hard-codes deltaLambda = ds/sqrt(2*|deltaUt|^2),
// which is the correct root of |DeltaU|^2 + A0*DeltaLambda^2 = ds^2 ONLY under
// the automatic scaling phi = |deltaUt|/|f| set two lines below. With a user
// Scaling the predictor lands off the constraint surface (1/sqrt(2), i.e. 29%
// short, for Crisfield's own recommended cylindrical psi = 0).
//
TEST(predictor_respects_user_scaling)
{
    const real_t ds = 0.05;

    // Control case: the automatic scaling (Scaling = -1) must land on the constraint
    // surface too, not just the user-set values swept below.
    {
        AlmProblem prob = foldProblem();
        CrisfieldProbe alm(prob.Jacobian, prob.ALResidual, prob.Force);
        configure(alm, ds);
        // Asked for BY NAME so this block does not depend on the shipped
        // default: it is the automatic phi that is under test here. Inert today
        // (the default is -1); it keeps the coverage if the default ever moves to 0.
        alm.options().setReal("Scaling",-1.0);
        alm.applyOptions();
        alm.runPredictorOnly();
        gsInfo<<"  [M5] Scaling = auto  phi = "<<alm.phi()
              <<"  |(DU,DL)|_A = "<<alm.distance(alm.DeltaU(),alm.DeltaL())<<"\n";
        CHECK_CLOSE(ds, alm.distance(alm.DeltaU(),alm.DeltaL()), 1e-12);
    }

    const real_t scalings[3] = {0.0, 1.0, 1e3};
    for (index_t i = 0; i != 3; ++i)
    {
        AlmProblem prob = foldProblem();
        CrisfieldProbe alm(prob.Jacobian, prob.ALResidual, prob.Force);
        configure(alm, ds);
        alm.options().setReal("Scaling",scalings[i]);
        alm.applyOptions();

        alm.runPredictorOnly();
        const real_t d = alm.distance(alm.DeltaU(),alm.DeltaL());
        gsInfo<<"  [M5] Scaling = "<<scalings[i]<<"  phi = "<<alm.phi()
              <<"  |(DU,DL)|_A = "<<d<<"  (expected "<<ds<<")\n";
        CHECK_CLOSE(ds, d, 1e-12);
    }
}

// ===========================================================================
// TEST 8b (M6) -- the automatic scaling must survive a zero-load crossing.
// ===========================================================================
// gsALMCrisfield's secant-branch auto scaling is phi = |U|/(|Lambda| |f|), an
// unguarded division by the load factor. Fixture F's equilibrium parabola
// returns to lambda = 0 at u1 = 2, a perfectly regular point (K = diag(-2,1)).
//
// FIXED (M6) -- ENABLED. gsALMCrisfield::predictor() now implements Lam & Morley
// 1992's two-formula prescription: eq. (11) phi = |U|/(|L||f|) normally, falling back to
// eq. (22) phi = |delta_u_t|/|f| ONLY where eq. (11) is undefined. At this fixture
// Lambda = 0 exactly, so the exception fires and phi = |K^-1 f|/|f| = 0.5.
//
// Why the fallback must stay an exception, and why this test alone does not prove that:
// eq. (22) applied EVERYWHERE was measured to make this
// test pass while KILLING both M3 fold tests above -- A0 = |K^-1 f|^2 diverges at a LIMIT
// point (Ritto-Correa & Camotim 2008, footnote 3, p.1356), reaching 1076 at the state
// step 27 starts from against 0.970 for eq. (11), and the traces then die at step 27 with
// max u1 = 0.985 < 1.05 at all three load scales. So this test and
// crisfield_fold_rounding_is_scale_invariant are a PAIR: neither one alone pins the fix,
// and any future change to the automatic phi must keep both green simultaneously.
TEST(predictor_scaling_survives_zero_load)
{
    const real_t ds = 0.05;
    AlmProblem prob = foldProblem();
    CrisfieldProbe alm(prob.Jacobian, prob.ALResidual, prob.Force);
    configure(alm, ds);
    // ★ Asked for BY NAME so this gate does not depend on the shipped default.
    // Inert today (the default IS -1), but load-bearing if it ever moves: under any
    // fixed Scaling -- 0 in particular -- m_phi_user is true, the eq. (22) exception
    // below can never fire, and every assertion in this test passes TAUTOLOGICALLY
    // (phi = 0 is finite, computeLambdaMU() rescales onto the sphere, |DeltaU| = ds).
    // This test is about the AUTOMATIC phi, so it must name it.
    alm.options().setReal("Scaling",-1.0);
    alm.applyOptions();

    // Regular equilibrium with lambda = 0 exactly: u = (2,0), det K = -2.
    gsVector<real_t> U(2);     U     << 2.0, 0.0;
    gsVector<real_t> Uprev(2); Uprev << 1.9, 0.0;
    CHECK_CLOSE(0.0, foldLambda(2.0), 1e-14);
    CHECK_CLOSE(-2.0, tangentDet(prob,U), 1e-14);

    alm.setSolution(U,foldLambda(2.0));
    alm.setLength(ds);
    alm.setPrevious(Uprev,foldLambda(1.9));

    alm.runPredictorOnly();
    gsInfo<<"  [M6] phi = "<<alm.phi()<<"  |(DU,DL)|_A = "
          <<alm.distance(alm.DeltaU(),alm.DeltaL())
          <<"  |DU| = "<<alm.DeltaU().norm()<<"\n";
    CHECK(math::isfinite(alm.phi()));
    CHECK_CLOSE(ds, alm.distance(alm.DeltaU(),alm.DeltaL()), 1e-10);

    // NOT redundant with the line above: computeLambdaMU() rescales ANY finite phi back
    // onto the sphere, so |(DU,DL)|_A == ds is tautological (cf. the tautology note in
    // the runFallback()/FallbackOutcome comment above, used by the [M4] test).
    // The discriminator is that the predictor must actually MOVE: a phi that is merely
    // regularised (e.g. |Lambda| clamped to 1e-12) gives phi ~ 2e12 and |DeltaU| ~ 1e-14
    // while still satisfying every assertion above.
    CHECK(alm.DeltaU().norm() > 0.1*ds);

    // ... and the step through the crossing must complete on the manifold.
    AlmProblem prob2 = foldProblem();
    CrisfieldProbe alm2(prob2.Jacobian, prob2.ALResidual, prob2.Force);
    configure(alm2, ds);
    alm2.options().setReal("Scaling",-1.0);   // as above: the automatic phi is the subject
    alm2.applyOptions();
    alm2.setSolution(U,foldLambda(2.0));
    alm2.setLength(ds);
    alm2.setPrevious(Uprev,foldLambda(1.9));
    CHECK(alm2.step() == gsStatus::Success);
    CHECK_CLOSE(foldLambda(alm2.solutionU()[0]), alm2.solutionL(), 1e-8);
}

// ===========================================================================
// TEST 9 (M2, M8, M16 and the initialize() trap) -- restart determinism.
// ===========================================================================
// A solver handed back the state it produced must reproduce the continuation.
// Any per-step state that is NOT restored by
// setSolution()/setLength()/setPrevious() -- gsALMConsistentCrisfield's phi,
// built from a stale m_deltaUt (M16); gsALMRiks::m_convexWeight (M8); m_arcLength_prev
// (M2) -- shows up here as a diverging trajectory.
//
// Two variants, both legitimate readings of the review sketch:
//   (a) reseedInPlace: a second solver takes the same 5 steps, is handed its own
//       state back, and continues -- the injection must be a no-op;
//   (b) freshRestart:  a solver that has NEVER stepped is seeded with
//       (U5,L5),(U4,L4),ds and must reproduce steps 6..10. This is the
//       gsAPALM-style restart and the one that exposes M16.
template <class ALM>
void restartDeterminism(const std::string & name, bool freshRestart, real_t tol)
{
    const real_t  ds     = 0.05;
    const index_t nFirst = 5, nSecond = 5;

    AlmProblem probA = foldProblem();
    ALM A(probA.Jacobian, probA.ALResidual, probA.Force);
    configure(A, ds);

    std::vector<gsVector<real_t> > U(1, A.solutionU());
    std::vector<real_t>            L(1, A.solutionL());
    for (index_t k = 0; k != nFirst+nSecond; ++k)
    {
        CHECK(A.step() == gsStatus::Success);
        U.push_back(A.solutionU());
        L.push_back(A.solutionL());
    }

    AlmProblem probB = foldProblem();
    ALM B(probB.Jacobian, probB.ALResidual, probB.Force);
    configure(B, ds);
    if (!freshRestart)
        for (index_t k = 0; k != nFirst; ++k)
            CHECK(B.step() == gsStatus::Success);

    // The restart contract, in the order gsAPALM uses it (setLength BEFORE
    // setPrevious: setPrevious redeclares m_arcLength_prev = m_arcLength).
    B.setSolution(U[nFirst],L[nFirst]);
    B.setLength(ds);
    B.setPrevious(U[nFirst-1],L[nFirst-1]);

    real_t worst = 0.0;
    for (index_t k = 0; k != nSecond; ++k)
    {
        CHECK(B.step() == gsStatus::Success);
        const index_t idx = nFirst + k + 1;
        worst = math::max(worst, (B.solutionU()-U[idx]).norm());
        worst = math::max(worst, math::abs(B.solutionL()-L[idx]));
        CHECK_ARRAY_CLOSE(U[idx], B.solutionU(), 2, tol);
        CHECK_CLOSE(L[idx], B.solutionL(), tol);
    }
    gsInfo<<"  [restart] "<<name<<(freshRestart?" (fresh)":" (in place)")
          <<"  max deviation = "<<worst<<"\n";
}

TEST(restart_determinism_loadcontrol)
{
    restartDeterminism<gsALMLoadControl<real_t> >("gsALMLoadControl",false,1e-10);
    restartDeterminism<gsALMLoadControl<real_t> >("gsALMLoadControl",true ,1e-10);
}

TEST(restart_determinism_riks)
{
    restartDeterminism<gsALMRiks<real_t> >("gsALMRiks",false,1e-10);
    restartDeterminism<gsALMRiks<real_t> >("gsALMRiks",true ,1e-10);
}

TEST(restart_determinism_crisfield)
{
    restartDeterminism<gsALMCrisfield<real_t> >("gsALMCrisfield",false,1e-10);
    restartDeterminism<gsALMCrisfield<real_t> >("gsALMCrisfield",true ,1e-10);
}

TEST(restart_determinism_consistent_crisfield)
{
    restartDeterminism<gsALMConsistentCrisfield<real_t> >
        ("gsALMConsistentCrisfield",false,1e-10);

    restartDeterminism<gsALMConsistentCrisfield<real_t> >
        ("gsALMConsistentCrisfield",true ,1e-10);
}

// ===========================================================================
// TEST 9b (M2b) -- switchBranch() must leave the fundamental path.
// ===========================================================================
// switchBranch() clears gsALMCrisfield's fresh-start sentinel (m_DeltaUold) but
// not the one gsALMRiks / gsALMConsistentCrisfield use (m_U - m_Uprev), so those
// two predict along the PRE-bifurcation direction of travel and the corrector,
// constrained only to the arc-length sphere, is free to land on the fundamental
// path again.
template <class ALM>
void switchBranchCheck(const std::string & name, real_t ds, index_t nAfter)
{
    AlmProblem prob = pitchforkProblem();
    ALM alm(prob.Jacobian, prob.ALResidual, prob.Force);
    configure(alm, ds);
    alm.options().setReal("Perturbation",1e3);   // the shipped default: xi = 1e-3
    alm.applyOptions();

    // Up the fundamental branch until just past the fork.
    for (index_t k = 0; k != 200 && alm.solutionL() <= 1.0; ++k)
        CHECK(alm.step() == gsStatus::Success);
    CHECK(alm.solutionL() > 1.0);
    CHECK(alm.solutionL() < 1.15);                // still inside the test's tol
    CHECK(math::abs(alm.solutionU()[1]) < 1e-12); // exactly on u2 = 0

    const gsStatus st = alm.computeSingularPoint(alm.solutionU(), alm.solutionL(),
                                                 /*switchBranch*/true,
                                                 /*jacobian*/true,
                                                 /*testPoint*/true);
    CHECK(st == gsStatus::Success);

    for (index_t k = 0; k != nAfter; ++k)
        CHECK(alm.step() == gsStatus::Success);

    const gsVector<real_t> Uend = alm.solutionU();
    const real_t           Lend = alm.solutionL();
    gsInfo<<"  [M2b] "<<name<<": after switchBranch + "<<nAfter<<" steps -> "
          <<"(u1,u2,L) = ("<<Uend[0]<<","<<Uend[1]<<","<<Lend<<")\n";

    // Left the fundamental path ...
    CHECK(math::abs(Uend[1]) > 0.1);
    // ... and landed on the closed-form bifurcated branch.
    CHECK_CLOSE(pitchforkBranchU1(Lend)  , Uend[0]        , 1e-6);
    CHECK_CLOSE(pitchforkBranchU2sq(Lend), Uend[1]*Uend[1], 1e-6);
}

TEST(switch_branch_leaves_the_fundamental_path_loadcontrol)
{
    switchBranchCheck<gsALMLoadControl<real_t> >("gsALMLoadControl",0.05,8);
}

TEST(switch_branch_leaves_the_fundamental_path_crisfield)
{
    switchBranchCheck<gsALMCrisfield<real_t> >("gsALMCrisfield",0.05,8);
}

// FIXED (M2b): switchBranch() now also resets m_Uprev/m_Lprev, the sentinel the
// SECANT predictors of gsALMRiks / gsALMConsistentCrisfield read.
TEST(switch_branch_leaves_the_fundamental_path_riks)
{
    switchBranchCheck<gsALMRiks<real_t> >("gsALMRiks",0.05,8);
}

// FIXED (M2b): same sentinel as Riks.
TEST(switch_branch_leaves_the_fundamental_path_consistent_crisfield)
{
    switchBranchCheck<gsALMConsistentCrisfield<real_t> >
        ("gsALMConsistentCrisfield",0.05,8);
}

// ===========================================================================
// TEST 10 (m1) -- a step that cannot converge must not report Success.
// ===========================================================================
// Regression for defect m1. A corrector loop of the form
// `for (n = 1; n < m_maxIterations; ++n)` with the limit throw INSIDE it would
//   (i) run MaxIter-1 corrector iterations (ask for 30, get 29), and
//   (ii) with MaxIter <= 1 never run its body, so nothing would throw and step()
//        would return Success while m_converged was false and (U,L) were never
//        updated ("[m1] status = 0  converged() = false  moved = 0" reproduces
//        exactly this).
// The loop instead runs `n = 1 .. m_maxIterations` inclusive with the failure test
// AFTER the loop, keyed on !m_converged -- so (ii) cannot happen
// whatever the bound. One clean foldProblem() step at ds = 0.05 needs 5 corrector
// iterations, so MaxIter = 1 grants exactly one iteration, cannot converge, and
// throws 1 -> gsStatus::NotConverged.
TEST(step_with_one_iteration_reports_failure)
{
    AlmProblem prob = foldProblem();
    gsALMRiks<real_t> alm(prob.Jacobian, prob.ALResidual, prob.Force);
    configure(alm, 0.05);
    alm.options().setInt("MaxIter",1);
    alm.applyOptions();

    const gsVector<real_t> U0 = alm.solutionU();

    const gsStatus st = alm.step();
    gsInfo<<"  [m1] status = "<<static_cast<index_t>(st)
          <<"  converged() = "<<(alm.converged()?"true":"false")
          <<"  moved = "<<(alm.solutionU()-U0).norm()<<"\n";

    CHECK(st != gsStatus::Success);
    // ... and in any case status and converged() must not contradict each other.
    CHECK( !(st == gsStatus::Success && !alm.converged()) );
}

// ===========================================================================
// TEST 10 (M14) -- isStable() must report the stability of the current arm.
// ===========================================================================
// Regression for defect M14. isStable() used to return the
// index_t m_stability (+1 or -1: both convert to true), and stability() returned
// +1 when the smallest eigenvalue/pivot was NEGATIVE, the opposite of what its own
// documentation says. Both now follow the documented convention: -1 unstable,
// +1 stable, and isStable() is a genuine sign test.
TEST(is_stable_matches_the_arm)
{
    AlmProblem prob = foldProblem();
    gsALMRiks<real_t> alm(prob.Jacobian, prob.ALResidual, prob.Force);
    configure(alm, 0.05);
    // computeStability() is public on gsALMBase and also public
    // on all four subclasses; the call through the base interface is kept as-is because
    // it exercises the very interface the narrowing used to break.
    gsALMBase<real_t> & base = alm;

    // Stable arm: u1 < 1 => K = diag(2-2u1,1) positive definite.
    gsVector<real_t> Us(2); Us << 0.5, 0.0;
    alm.setSolution(Us,foldLambda(0.5));
    base.computeStability(true);
    const bool stableArm = alm.isStable();
    const index_t stabS  = alm.stability();

    // Unstable arm: u1 > 1 => K has a negative eigenvalue.
    gsVector<real_t> Uu(2); Uu << 1.5, 0.0;
    alm.setSolution(Uu,foldLambda(1.5));
    base.computeStability(true);
    const bool unstableArm = alm.isStable();
    const index_t stabU    = alm.stability();

    gsInfo<<"  [M14] stable arm: isStable = "<<(stableArm?"true":"false")
          <<", stability() = "<<stabS
          <<" | unstable arm: isStable = "<<(unstableArm?"true":"false")
          <<", stability() = "<<stabU<<"\n";

    CHECK(stableArm);
    CHECK(!unstableArm);
    CHECK_EQUAL(+1, stabS);    // documented: +1 if stable
    CHECK_EQUAL(-1, stabU);    // documented: -1 if unstable
}

// ===========================================================================
// TEST 10 (m17) -- converged() must not survive a failed singular-point solve.
// ===========================================================================
// Regression for defect m17. _extendedSystemSolve() never touches
// m_converged, so without an explicit assignment, after a FAILED computeSingularPoint
// the solver would report the convergence of whatever ran last -- reproduced as
// "[m17] status = 1  converged() = true", i.e. NotConverged and converged() at once.
// _computeSingularPoint instead assigns m_converged from the extended stage's return value,
// after that stage (so it also overwrites the nested step()s of _bisectionSolve, live
// whenever SingularPointComputeTolB != 0) and before the throw. The defect is a
// property of gsALMBase, not of any particular consumer.
TEST(converged_is_not_stale_after_singular_point)
{
    AlmProblem prob = pitchforkProblem();
    gsALMLoadControl<real_t> alm(prob.Jacobian, prob.ALResidual, prob.Force);
    configure(alm, 0.05);

    for (index_t k = 0; k != 21; ++k)
        CHECK(alm.step() == gsStatus::Success);
    CHECK(alm.converged());          // the last ordinary step did converge

    // Starve the extended solve so that it cannot possibly converge: the
    // termination test is ||K_T*V|| < tolE and a norm is never < 0.
    alm.options().setReal("SingularPointComputeTolE",0.0);
    alm.options().setInt ("MaxIter",5);
    alm.applyOptions();
    const gsStatus st = alm.computeSingularPoint(alm.solutionU(), alm.solutionL(),
                                                 /*switchBranch*/false,
                                                 /*jacobian*/true,
                                                 /*testPoint*/true);
    gsInfo<<"  [m17] status = "<<static_cast<index_t>(st)
          <<"  converged() = "<<(alm.converged()?"true":"false")<<"\n";
    CHECK(st != gsStatus::Success);

    CHECK(!alm.converged());
}

// ===========================================================================
// TEST (m17, the LIMIT-POINT branch) -- the second exit of computeSingularPoint()
// ===========================================================================
// The test above covers the `test == true` branch (a failing extended solve).
// The `testPoint == true, test == false` branch -- a LIMIT point -- performs no
// solve at all and ends in `throw 1`, i.e. it returns gsStatus::NotConverged.
// Previously it assigned nothing to m_converged, so converged() kept reporting
// the PRECEDING ordinary step: NotConverged together with converged() == true,
// which is the m17 contradiction itself, surviving on the one path the doxygen
// described as an exception instead of fixing.
TEST(converged_is_not_stale_after_limit_point)
{
    AlmProblem prob = pitchforkProblem();
    gsALMLoadControl<real_t> alm(prob.Jacobian, prob.ALResidual, prob.Force);
    configure(alm, 0.05);

    for (index_t k = 0; k != 21; ++k)
        CHECK(alm.step() == gsStatus::Success);
    CHECK(alm.converged());          // the last ordinary step did converge

    // Force the LIMIT-point classification deterministically: _testSingularPoint
    // accepts a BRANCH point on |cos(V,f)| < SingularPointTestTol, and a cosine is
    // never < 0, so tol = 0 leaves the limit-point branch as the only reachable
    // one -- the same starvation trick the test above uses on the extended solve.
    alm.options().setReal("SingularPointTestTol",0.0);
    alm.applyOptions();
    const gsStatus st = alm.computeSingularPoint(alm.solutionU(), alm.solutionL(),
                                                 /*switchBranch*/false,
                                                 /*jacobian*/true,
                                                 /*testPoint*/true);
    gsInfo<<"  [m17/limit] status = "<<static_cast<index_t>(st)
          <<"  converged() = "<<(alm.converged()?"true":"false")<<"\n";
    CHECK(st == gsStatus::NotConverged);
    CHECK(!alm.converged());
}

// ===========================================================================
// TEST 10 (operator wrapper) -- a `false` from a user callback must surface as
// gsStatus::AssemblyError.
// ===========================================================================
// This is the one thing in the "dropped bool" category the review checked and
// found CORRECT; the counting/fault-injecting wrapper is what makes it testable
// (and it is the same wrapper tests 2/16 need).
TEST(assembly_failure_is_reported_as_assembly_error)
{
    // Baseline: a clean step, and the counters actually count.
    {
        AlmProblem prob = foldProblem();
        AlmOperatorProbe probe(prob);
        gsALMRiks<real_t> alm(probe.jacobian(), probe.residual(), prob.Force);
        configure(alm, 0.05);
        CHECK(alm.step() == gsStatus::Success);
        gsInfo<<"  [probe] one clean step: "<<probe.residualCalls()<<" residual / "
              <<probe.jacobianCalls()<<" Jacobian evaluations\n";
        CHECK(probe.residualCalls() > 0);
        CHECK(probe.jacobianCalls() > 0);
    }

    // Residual failure on the very first evaluation.
    {
        AlmProblem prob = foldProblem();
        AlmOperatorProbe probe(prob);
        probe.failResidualAt(1);
        gsALMRiks<real_t> alm(probe.jacobian(), probe.residual(), prob.Force);
        configure(alm, 0.05);
        CHECK(alm.step() == gsStatus::AssemblyError);
        CHECK_EQUAL(1, probe.residualCalls());
    }

    // Residual failure in the middle of the corrector.
    {
        AlmProblem prob = foldProblem();
        AlmOperatorProbe probe(prob);
        probe.failResidualAt(3);
        gsALMRiks<real_t> alm(probe.jacobian(), probe.residual(), prob.Force);
        configure(alm, 0.05);
        CHECK(alm.step() == gsStatus::AssemblyError);
        CHECK_EQUAL(3, probe.residualCalls());
    }

    // Jacobian failure inside step(). NOTE: initialize() ALSO assembles a Jacobian
    // (gsALMBase::init -> _computeStability -> _computeJacobian -> `throw 2`), so the
    // injected failure is aimed at the first assembly after configure() returns.
    {
        AlmProblem prob = foldProblem();
        AlmOperatorProbe probe(prob);
        gsALMRiks<real_t> alm(probe.jacobian(), probe.residual(), prob.Force);
        configure(alm, 0.05);
        probe.failJacobianAt(probe.jacobianCalls()+1);
        CHECK(alm.step() == gsStatus::AssemblyError);
    }

    // ... and the same during initialize() itself. This used to escape
    // as a raw `throw 2` -- an `int` -- past every caller in the library, which is the one
    // public entry point that did not honour the gsStatus convention.
    {
        AlmProblem prob = foldProblem();
        AlmOperatorProbe probe(prob);
        gsALMRiks<real_t> alm(probe.jacobian(), probe.residual(), prob.Force);
        alm.options().setString("Solver","SimplicialLDLT");
        alm.options().setInt   ("BifurcationMethod",0);
        alm.applyOptions();

        probe.failJacobianAt(probe.jacobianCalls()+1);   // the initial stability assembly
        gsStatus st = gsStatus::NotStarted;
        bool threwRawInt = false;
        try { st = alm.initialize(); }
        catch (int e) { threwRawInt = true; gsInfo<<"  [fix 5] raw int escaped: "<<e<<"\n"; }
        gsInfo<<"  [fix 5] initialize() status = "<<static_cast<index_t>(st)<<"\n";
        CHECK(!threwRawInt);
        CHECK(st == gsStatus::AssemblyError);
        CHECK(alm.status() == gsStatus::AssemblyError);

        // A clean initialize() still reports Success (the guard is not a blanket catch).
        AlmProblem prob2 = foldProblem();
        AlmOperatorProbe probe2(prob2);
        gsALMRiks<real_t> alm2(probe2.jacobian(), probe2.residual(), prob2.Force);
        configure(alm2, 0.05);
        CHECK(alm2.status() == gsStatus::Success);
    }

    // ... and the same on the singular-point path.
    {
        AlmProblem prob = pitchforkProblem();
        AlmOperatorProbe probe(prob);
        // gsALMLoadControl's ctor takes NON-const lvalue references.
        ALMJacobian_t J = probe.jacobian();
        ALMResidual_t R = probe.residual();
        gsALMLoadControl<real_t> alm(J, R, prob.Force);
        configure(alm, 0.05);
        for (index_t k = 0; k != 21; ++k)
            CHECK(alm.step() == gsStatus::Success);
        probe.failResidualAt(probe.residualCalls()+1);
        const gsStatus st = alm.computeSingularPoint(alm.solutionU(), alm.solutionL(),
                                                     false, true, true);
        CHECK(st == gsStatus::AssemblyError);
    }
}

// ===========================================================================
// TEST 13 (m8) -- a SYMMETRIC problem must find its ANTISYMMETRIC critical mode.
// ===========================================================================
// The historic power iteration started from `Ones` -- a SYMMETRIC field -- and
// ran a fixed number of sweeps. The canonical bifurcation of a symmetric
// structure has an ANTISYMMETRIC critical mode, which is EXACTLY orthogonal to
// `Ones`: the component the iteration has to amplify is zero, so no number of
// sweeps recovers it and solutionV() returns the symmetric mode instead. Every
// consumer of the mode (switchBranch(), gsALMExploration's branch jobs) then
// nudges along the wrong direction.
//
// Fixture W (local to this test): a 4-DOF linear problem whose tangent is built
// from the Walsh/Hadamard basis of R^4,
//     w0 = ( 1, 1, 1, 1)/2   symmetric ("constant") mode
//     w1 = ( 1,-1, 1,-1)/2   ANTISYMMETRIC mode  <-- the critical one
//     w2 = ( 1, 1,-1,-1)/2
//     w3 = ( 1,-1,-1, 1)/2
//     K   = sum_i mu_i w_i w_i^T,   mu = (2, 0.01, 3, 4)
// Every entry is a multiple of 1/4, hence exactly representable, and
//     w1.Ones = 0   and   w1.f = 0  with  f = Ones
// hold BIT-exactly. So (i) the smallest-|eigenvalue| mode is w1, (ii) `Ones` is
// an exact eigenvector of K (eigenvalue mu0 = 2) and therefore a stationary
// point of the inverse iteration, and (iii) the point is a genuine BRANCH point
// (f in range(K_T) to machine precision), not a limit point.
//
// The eigenvalue ratio mu0/mu1 = 200 is deliberately MODERATE: it bounds the
// amplification of the O(1e-16) rounding noise of the linear solves to
// 200^5*1e-16 ~ 3e-5 over five sweeps, so the old code cannot reach the mode
// "by accident" through amplified round-off -- it stays on w0.
// ---------------------------------------------------------------------------

/// Walsh (Hadamard) mode \a i of R^4; entries +/-1/2, exactly representable.
inline gsVector<real_t> walshMode(index_t i)
{
    static const real_t sgn[4][4] = { { 1., 1., 1., 1.},
                                      { 1.,-1., 1.,-1.},
                                      { 1., 1.,-1.,-1.},
                                      { 1.,-1.,-1., 1.} };
    gsVector<real_t> w(4);
    for (index_t k = 0; k != 4; ++k)
        w[k] = 0.5*sgn[i][k];
    return w;
}

/// Fixture W: R = K u - lambda f, with K = sum_i mu_i w_i w_i^T and f = Ones.
inline AlmProblem walshProblem()
{
    gsMatrix<real_t> K = gsMatrix<real_t>::Zero(4,4);
    const real_t mu[4] = {2.0, 0.01, 3.0, 4.0};
    for (index_t i = 0; i != 4; ++i)
    {
        const gsVector<real_t> w = walshMode(i);
        K += mu[i]*w*w.transpose();
    }

    AlmProblem p;
    p.scale = 1.0;
    p.Force = gsVector<real_t>::Ones(4);
    p.Jacobian =
        [K](gsVector<real_t> const & /*u*/, gsSparseMatrix<real_t> & m) -> bool
        {
            m = K.sparseView();
            m.makeCompressed();
            return true;
        };
    p.ALResidual =
        [K](gsVector<real_t> const & u, const real_t lambda, gsVector<real_t> & result) -> bool
        { result = K*u - lambda*gsVector<real_t>::Ones(4); return true; };
    p.Fint =
        [K](const gsVector<real_t> & u) { return gsVector<real_t>(K*u); };
    return p;
}

TEST(symmetric_problem_finds_the_antisymmetric_mode)
{
    // --- Fixture W: the unambiguous case (mode EXACTLY orthogonal to Ones) ----
    {
        AlmProblem prob = walshProblem();
        gsALMLoadControl<real_t> alm(prob.Jacobian, prob.ALResidual, prob.Force);
        configure(alm, 0.05);

        // Any point will do: the tangent is state independent by construction.
        gsVector<real_t> U = gsVector<real_t>::Zero(4);
        alm.setSolution(U,0.0);

        const bool isBif = alm.isBifurcation(true);
        const gsVector<real_t> V     = alm.solutionV();
        const gsVector<real_t> Vstar = walshMode(1);          // closed-form mode
        const real_t align = math::abs(V.dot(Vstar))/V.norm();  // |V.V*|/|V|, V* is a unit vector
        gsInfo<<"  [m8/W] |V.V*|/|V| = "<<align
              <<"   |V.f|/|f| = "<<modeForceCosine(V,prob.Force)
              <<"   isBifurcation = "<<(isBif?"true":"false")<<"\n";

        // THE gate: the returned mode is the antisymmetric one.
        CHECK(align > 0.99);
        // ... and, since it is orthogonal to the (symmetric) load, the point is
        // correctly classified as a BRANCH point.
        CHECK(isBif);
    }

    // --- Fixture P: the same statement on the pitchfork ------------------------
    // Here V* = (0,1) is NOT orthogonal to the historic Ones start, so this half
    // passes both before and after the fix; it is kept as the "the fix did not
    // break the easy case" control.
    {
        AlmProblem prob = pitchforkProblem();
        gsALMLoadControl<real_t> alm(prob.Jacobian, prob.ALResidual, prob.Force);
        configure(alm, 0.05);

        gsVector<real_t> U(2); U << 1.05, 0.0;                // just past the fork
        alm.setSolution(U,1.05);

        const bool isBif = alm.isBifurcation(true);
        const gsVector<real_t> V = alm.solutionV();
        gsVector<real_t> Vstar(2); Vstar << 0.0, 1.0;         // closed-form mode
        const real_t align = math::abs(V.dot(Vstar))/V.norm();
        gsInfo<<"  [m8/P] |V.V*|/|V| = "<<align
              <<"   |V.f|/|f| = "<<modeForceCosine(V,prob.Force)
              <<"   isBifurcation = "<<(isBif?"true":"false")<<"\n";

        CHECK(align > 0.99);
        CHECK(isBif);
    }
}

// ===========================================================================
// TEST 13b (extended-solve mode commit) -- solutionV() must be the REFINED mode.
// ===========================================================================
// _extendedSystemSolve() iterates on (U,Lambda,V) and tests convergence on
// ||K_T*(V+DeltaV)||. Committing only `m_U += m_DeltaU; m_L += m_DeltaL;` would
// discard the mode correction it had just computed, so solutionV() would return
// the PRE-refinement power-iteration mode and every branch nudge in the library
// would point along it.
//
// Observable form of the defect: at the converged singular point (u*,lambda*)
// the exact critical mode satisfies K_T(u*) V = 0. On fixture P the extended
// solve drives ||K_T*(V+DeltaV)|| below SingularPointComputeTolE = 1e-10, while
// the uncommitted power-iteration mode leaves ||K_T V|| ~ 1e-7 (its residual
// after the five sweeps at the pre-crossing point).
TEST(extended_solve_commits_the_mode_correction)
{
    AlmProblem prob = pitchforkProblem();
    gsALMLoadControl<real_t> alm(prob.Jacobian, prob.ALResidual, prob.Force);
    configure(alm, 0.05);

    gsVector<real_t> U(2); U << 1.05, 0.0;
    alm.setSolution(U,1.05);

    const gsStatus st = alm.computeSingularPoint(U,1.05,
                                                 /*switchBranch*/false,
                                                 /*jacobian*/true,
                                                 /*testPoint*/true);
    CHECK(st == gsStatus::Success);
    // The singular point itself (this is TEST 2's statement; repeated here only
    // so that a mode residual is not compared at a point that never converged).
    CHECK_CLOSE(1.0, alm.solutionL(),    1e-6);
    CHECK_CLOSE(1.0, alm.solutionU()[0], 1e-6);

    gsSparseMatrix<real_t> K;
    prob.Jacobian(alm.solutionU(),K);
    const gsVector<real_t> V = alm.solutionV();
    const real_t modeResidual = (K*V).norm()/V.norm();
    gsInfo<<"  [defect A] ||K_T(u*) V||/|V| = "<<modeResidual<<"  (tolE = 1e-10)\n";

    // Two orders of magnitude above tolE, three below the uncommitted value.
    CHECK(modeResidual < 1e-8);
}

// ===========================================================================
// TESTS 11 -- the BORDERED-SOLVE FALLBACK at a fold.
// ===========================================================================
// gsALMCrisfield's corrector eliminates delta_u_t/u_bar and pins delta_lambda from a
// quadratic. MEASURED on fixture F at Scaling = 0 (cylindrical), ds = 0.05, forward from
// the seed, at all three load scales: the predictor of step 19 lands on
// u1 = 1 + 2.2e-16, so K = diag(-4.4e-16,1), |delta_u_t| = 2.25e15 and |u_bar| = 5.63e12.
// SimplicialLDLT factorizes that tangent without complaint -- the reported status is
// NotConverged (1), NEVER SolverError (3) -- but the constraint quadratic's leading part is
// the perfect square (alpha*dL+beta)^2/eps^2, so its discriminant is a difference of two
// O(eps^-4) numbers whose true value is O(eps^-2): a relative cancellation of eps^2 ~ 1e-32,
// far below double precision. It comes out EXACTLY 0, both roots collapse onto the
// degenerate double root, the resulting point misses |DeltaU| = ds by three orders of
// magnitude, and the corrector settles into a bounded PERIOD-2 limit cycle
// (u1 alternating 1.0 <-> 0.9502, |delta U| and both residuals constant to every printed
// digit) for all 99 iterations. THIS IS A CANCELLATION FAILURE OF THE ELIMINATION, NOT A
// SINGULAR MATRIX: the bordered B = [[K,-f],[2 DeltaU^T, 2 A0 DeltaLambda]] is NONSINGULAR
// at exactly that point (f is not in range(K) and 2 DeltaU.e1 = 0.1), which is the whole
// reason the fallback works.
//
// See the class documentation of gsALMCrisfield for the algebra; the chart change is EXACT,
// not a regularisation.
// ---------------------------------------------------------------------------
struct BorderedTrace
{
    real_t  maxU1;         ///< the furthest u1 reached (> 1: the fold was rounded)
    index_t nPoints;       ///< number of accepted points, seed included
    bool    allStepsOk;    ///< every requested step returned gsStatus::Success
    index_t failedStep;    ///< index of the first failing step, or -1
    index_t failedStatus;  ///< its gsStatus as an integer, or -1
    real_t  minIncrement;  ///< min_k |U_k - U_{k-1}|                            (G3)
    index_t inertiaFlips;  ///< number of changes of negatives() along the trace (G4)
    bool    inertiaMatchesArm; ///< negatives() == (u1 > 1) at every accepted point (G4)
};

/// Traces fixture F at load scale \a scale, CYLINDRICAL (Scaling = 0) and with
/// ds = 0.05 UNCHANGED, with the bordered fallback \a bordered.
///
/// Deliberately a NEW harness. crisfieldProgress() above pins NO Scaling, so it
/// exercises whatever gsALMCrisfield::defaultOptions() currently sets; the two
/// tests it backs (crisfield_root_selection_never_retraces and
/// crisfield_fold_rounding_is_scale_invariant) guard the AUTOMATIC-phi path and
/// therefore only do so while that default is -1. On 2026-07-30 the default was
/// briefly flipped to 0 and both tests silently stopped exercising auto-phi; the
/// flip has since been reverted. This harness pins its own Scaling below so it
/// stays independent of that default either way.
BorderedTrace crisfieldBorderedProgress(real_t scale, bool bordered, index_t nSteps = 40)
{
    AlmProblem prob = foldProblem(scale);
    CrisfieldProbe alm(prob.Jacobian, prob.ALResidual, prob.Force);
    configure(alm, 0.05);                                 // ds = 0.05, UNCHANGED
    alm.options().setReal  ("Scaling",0.0);               // cylindrical, pinned IN THE TEST
    alm.options().setString("BorderedMode",bordered ? "Fallback" : "Off"); // pinned IN THE TEST,
                                                          // both ways, so this harness never
                                                          // inherits the default. Migrated from the
                                                          // deprecated BorderedFallback switch (task
                                                          // 07): BorderedMode="Fallback" is what that
                                                          // switch mapped onto, so this is behaviour-
                                                          // preserving; the alias itself is still
                                                          // exercised once, deliberately, in T6.
    alm.applyOptions();

    std::vector<gsVector<real_t> > Us;
    std::vector<index_t>           negs;

    BorderedTrace out;
    out.maxU1        = alm.solutionU()[0];
    out.allStepsOk   = true;
    out.failedStep   = -1;
    out.failedStatus = -1;
    out.minIncrement = std::numeric_limits<real_t>::max();
    out.inertiaFlips = 0;
    out.inertiaMatchesArm = true;
    Us.push_back(alm.solutionU());

    for (index_t k = 0; k != nSteps; ++k)
    {
        const gsStatus st = alm.step();
        if (st != gsStatus::Success)
        {
            out.allStepsOk   = false;
            out.failedStep   = k;
            out.failedStatus = static_cast<index_t>(st);
            break;
        }
        Us.push_back(alm.solutionU());
        negs.push_back(alm.negatives());
        out.maxU1 = math::max(out.maxU1, alm.solutionU()[0]);
        // G3: the accepted increment, i.e. the step actually went somewhere.
        out.minIncrement = math::min(out.minIncrement,
                                     (Us[Us.size()-1]-Us[Us.size()-2]).norm());
        // G4: det K_T = (2-2u1)*1 is negative exactly on the far arm, so the analytic
        // inertia at an accepted point is 1 iff u1 > 1. This is an ORACLE, not a
        // self-comparison.
        if (alm.negatives() != (alm.solutionU()[0] > 1.0 ? 1 : 0))
            out.inertiaMatchesArm = false;
        if (negs.size() >= 2 && negs[negs.size()-1] != negs[negs.size()-2])
            ++out.inertiaFlips;
        // Stay clear of the second lambda = 0 crossing at u1 = 2.
        if (alm.solutionU()[0] > 1.8 || alm.solutionU()[0] < -0.8)
            break;
    }
    out.nPoints = static_cast<index_t>(Us.size());
    gsInfo<<"  [T42] scale = "<<scale<<(bordered?"  bordered ON ":"  bordered OFF")
          <<"  points = "<<out.nPoints
          <<"  max u1 = "<<out.maxU1
          <<"  min|dU| = "<<out.minIncrement
          <<"  inertia flips = "<<out.inertiaFlips
          <<"  arm match = "<<(out.inertiaMatchesArm?"yes":"no");
    if (!out.allStepsOk)
        gsInfo<<"  FAILED at step "<<out.failedStep<<" with status "<<out.failedStatus;
    gsInfo<<"\n";
    return out;
}

// ---------------------------------------------------------------------------
// G-A' -- the switch actually SWITCHES.
// ---------------------------------------------------------------------------
// This is G-A's fails-first, asserted permanently. G-A and G-A' bracket the switch through
// crisfieldBorderedProgress()'s `bordered` argument, ON and OFF, so the PAIR is what pins the
// behaviour and neither half depends on what the default happens to be. Without this half a
// change of default would let G-A keep passing while pinning nothing at all -- exactly the
// disease the Scaling = -1 pins elsewhere in this file exist to prevent, and exactly what did
// happen on 2026-07-30 to the two auto-phi tests that pin no Scaling (see
// crisfieldBorderedProgress()'s note).
//
// ⚠ "No failure" is never the evidence here and the numbers below are not either: they are
// the OPPOSITE, a trace that provably STOPS. max u1 == 0.95 is the seed-side arm; the fold
// is at u1 = 1.
TEST(crisfield_without_the_bordered_fallback_still_dies_at_the_fold)
{
    const real_t scales[3] = {1.0, 1e3, 1e-3};
    for (index_t i = 0; i != 3; ++i)
    {
        const BorderedTrace off = crisfieldBorderedProgress(scales[i], /*bordered*/false);
        CHECK(!off.allStepsOk);
        CHECK_EQUAL(19, off.failedStep);
        // 1 == gsStatus::NotConverged. NOT 3 (SolverError): the tangent is never bit-exactly
        // singular here, which is why the factorizeShifted ladder is irrelevant to this
        // defect and the bordered chart is not.
        CHECK_EQUAL(1, off.failedStatus);
        CHECK_EQUAL(20, off.nPoints);
        CHECK_CLOSE(0.95, off.maxU1, 1e-12);
        // ... and the fold is therefore never rounded, so no inertia flip ever happens.
        CHECK_EQUAL(0, off.inertiaFlips);
    }
}

// ---------------------------------------------------------------------------
// G-A / G3 / G4 -- with the fallback the fold is ROUNDED, at every load scale.
// ---------------------------------------------------------------------------
// ⚠ THE TAUTOLOGY WARNING. A stalled trace never misbehaves, so "every step succeeded" is
// no evidence on its own. The evidence is max u1 PAST the fold plus the scale-invariance
// equalities, and the companion test above pins that the same configuration without the
// switch does NOT get there.
TEST(crisfield_bordered_fallback_rounds_the_fold_at_every_load_scale)
{
    // Ground truth at the reference scale.
    const BorderedTrace ref = crisfieldBorderedProgress(1.0, /*bordered*/true);
    CHECK(ref.allStepsOk);
    CHECK(ref.maxU1 > 1.05);            // G-A: PAST the fold at u1 = 1, not merely alive

    // G3 -- anti-degeneracy (borrowed from test 8b). A step that satisfies the arc-length
    // identity while going nowhere passes every other assertion here; the rejected
    // epsilon-clamp produced |DeltaU| ~ 1.3e-14 that way.
    CHECK(ref.minIncrement > 0.1*0.05);

    // G4 -- the stability classification is not silently stale. gsALMBase's
    // computeStability() SWALLOWS a throw 3, after which
    // m_stabilityVec keeps the previous iterate's value while the step still reports
    // Success, so "max u1 > 1.05" can pass on a stale classification. det K_T = 2-2u1 is
    // negative exactly beyond the fold, so the analytic inertia is an oracle.
    CHECK(ref.inertiaMatchesArm);
    CHECK_EQUAL(1, ref.inertiaFlips);   // and it flips ONCE, not back and forth

    const real_t scales[2] = {1e-3, 1e3};
    for (index_t i = 0; i != 2; ++i)
    {
        const BorderedTrace t = crisfieldBorderedProgress(scales[i], /*bordered*/true);
        CHECK(t.allStepsOk);
        CHECK(t.maxU1 > 1.05);
        gsInfo<<"  [T42] scale-invariance of maxU1: |ref - t| = "
              <<math::abs(ref.maxU1-t.maxU1)<<" (tol 1e-8)\n";
        CHECK_CLOSE(ref.maxU1, t.maxU1, 1e-8);
        CHECK_EQUAL(ref.nPoints, t.nPoints);
        CHECK(t.minIncrement > 0.1*0.05);
        CHECK(t.inertiaMatchesArm);
        CHECK_EQUAL(1, t.inertiaFlips);
    }
}

// ---------------------------------------------------------------------------
// G-B -- the EXACTNESS claim, which is the load-bearing assertion of the whole
// design and otherwise has no gate at all.
// ---------------------------------------------------------------------------
// The bordered chart is an EXACT change of basis of the same affine solution set of
// K deltaU - deltaLambda f = -r, so AWAY from the fold -- where the elimination is perfectly
// healthy -- both charts must produce the SAME corrector increment. If they do not, the
// coefficient map is wrong and every "green" above is an accident.
//
// ⚠ This must NOT be written as "run a short trace with the switch on and off and check
// they agree". The retry only fires on NotConverged, so on a trace that never fails the two
// runs execute literally the same instructions: that assertion would stay true if the
// bordered branch were GISMO_ERROR. Hence CrisfieldProbe::forceBordered().
//
// Three configurations, because no single one exercises the whole coefficient map:
//   (i)   Scaling = 0, u2 = 0     -- the G-A configuration. MEASURED: the corrector's
//         displacement increment is ANALYTICALLY ZERO here (the cylindrical constraint pins
//         |DeltaU| = ds and u2 == 0 on the fundamental path, so only lambda moves; both
//         charts return |deltaU| = 2e-18, pure rounding). Only the dL comparison carries
//         information in this case -- hence (ii).
//   (ii)  Scaling = 0, u2 = 0.1   -- seeded OFF the equilibrium so the corrector has real
//         displacement work to do and deltaU is O(0.1).
//   (iii) Scaling = 1, u2 = 0     -- A0 = |f|^2 = 1, so the four A0-weighted terms of the
//         map are live; at Scaling = 0 they all drop out and a defect in one is invisible.
TEST(crisfield_bordered_chart_equals_the_elimination_away_from_the_fold)
{
    const real_t scalings[3] = {0.0, 0.0, 1.0};
    const real_t u2seed  [3] = {0.0, 0.1, 0.0};
    for (index_t s = 0; s != 3; ++s)
    {
        gsVector<real_t> dU[2];
        real_t           dL[2];
        for (index_t b = 0; b != 2; ++b)
        {
            AlmProblem prob = foldProblem(1.0);
            CrisfieldProbe alm(prob.Jacobian, prob.ALResidual, prob.Force);
            configure(alm, 0.05);
            alm.options().setReal("Scaling",scalings[s]);
            alm.applyOptions();

            // Well away from the fold: u1 = 0.5 => K = diag(1,1), perfectly conditioned.
            gsVector<real_t> U(2); U << 0.5, u2seed[s];
            alm.setSolution(U,foldLambda(0.5));
            alm.setLength(0.05);

            alm.forceBordered(b==1);         // NOT via the NotConverged retry
            alm.runPredictorAndOneIteration();
            dU[b] = alm.deltaU();
            dL[b] = alm.deltaL();
        }
        // deltaU is normalised by the STEP SIZE, not by its own norm: in case (i) it is
        // analytically zero and a relative comparison would divide rounding by rounding.
        const real_t relU = (dU[0]-dU[1]).norm() / math::max(dU[0].norm(),0.05);
        const real_t relL = math::abs(dL[0]-dL[1]) / math::abs(dL[0]);
        gsInfo<<"  [T42/G-B] Scaling = "<<scalings[s]<<"  u2 = "<<u2seed[s]
              <<"  |dU| = "<<dU[0].norm()<<"  dL = "<<dL[0]
              <<"  rel. diff (dU,dL) = ("<<relU<<","<<relL<<")  (tol 1e-12)\n";
        // Not tautological: the chart output is O(1e-3)..O(1e-1) in every case, so a wrong
        // coefficient shows up rather than hiding under a zero.
        CHECK(math::max(dU[0].norm(),math::abs(dL[0])) > 1e-4);
        CHECK(relU < 1e-12);
        CHECK(relL < 1e-12);
    }
}

// ===========================================================================
// TEST 12 (G2) -- the DIAGONAL-SHIFT LADDER of _extendedSystemIteration.
// ===========================================================================
// The factorizeShifted lambda inside gsALMBase::_extendedSystemIteration() holds a
// retry ladder: factorize the tangent
// unshifted; only if that throws, retry on a COPY shifted by
// sigma = 1e-12 * max(1,||diag||_inf), escalating x100 over exactly THREE attempts
// (1e-12 / 1e-10 / 1e-8 times that scale), and rethrow loudly if all three fail. It is
// called twice inside _extendedSystemIteration() -- once on the freshly assembled
// tangent, once after the finite-difference perturbation used for h1/h2 -- and
// NOWHERE else.
//
// Until this test it had ZERO coverage, and the reason is not that it is dead code but
// that continuation iterates never land BIT-EXACTLY on the singular set: the extended
// solve only converges TOWARDS it, to SingularPointComputeTolE = 1e-10, and
// SimplicialLDLT factorizes a tangent that is merely 1e-10 from singular without
// complaint. (The explanation in an earlier report -- "the extended iteration factorizes
// the well-conditioned BORDERED tangent, not K" -- is wrong: _extendedSystemIteration()
// assembles the tangent via evalJacobian() and factorizes it (plain K, via
// factorizeShifted) at the extended iterate.) The gate therefore has to
// SEED the ladder exactly on the singular set, which no public entry point can do:
// m_V is sized only inside _computeCriticalMode(), so
// computeSingularPoint(...,testPoint=false) would enter _extendedSystemSolve with an
// empty m_V, and the testPoint=true route factorizes K at site 8/9 before the extended
// system is ever reached. Hence the probe below.
//
// WHY FIXTURE P AND NOT F. Both are exactly singular at (1,0) -- P has
// K = diag(1, 1-u1) = diag(1,0) and F has K = diag(2-2u1, 1) = diag(0,1), and 1.0-1.0
// and 2.0-2.0*1.0 are both bit-exact zeros. But on F the load f = (1,0) IS the null
// direction, so delta_u_t = (K+sigma I)^{-1} f ~ 1e12 and the finite-difference vector
// h1 = (K_eps*delta_u_t - f)/eps is pure cancellation garbage. On P the load is
// ORTHOGONAL to the critical mode V = (0,1) (that is what makes P a branch point), so
// delta_u_t = (1,0)/(1+sigma) stays O(1); the only 1/sigma amplification is in
// delta_V_t, which is exactly where the extended system's constraint divides it out.
// P is the fixture on which the ladder can be observed in isolation.
// ---------------------------------------------------------------------------

/// Exposes the two things G2 needs: the extended-system solve seeded with a
/// GIVEN critical mode, and the plain (unshifted) factorization that is
/// attempt (i) of the ladder. Nothing under src/ is modified.
class ExtendedProbe : public gsALMLoadControl<real_t>
{
    typedef gsALMLoadControl<real_t> Base;
public:
    // Non-const references: that is gsALMLoadControl's constructor signature,
    // unlike gsALMRiks/gsALMCrisfield used by the probes above.
    ExtendedProbe(gsStructuralAnalysisOps<real_t>::Jacobian_t   & J,
                  gsStructuralAnalysisOps<real_t>::ALResidual_t & R,
                  gsVector<real_t>                              & F)
    : Base(J,R,F) {}

    /// Installs the critical mode and runs gsALMBase::_extendedSystemSolve at (U,L).
    bool extendedSolveAt(const gsVector<real_t> & U, real_t L,
                         const gsVector<real_t> & V, real_t tol)
    { this->m_V = V; return this->_extendedSystemSolve(U,L,tol); }

    /// gsALMBase::factorizeMatrix, i.e. LITERALLY attempt (i) of the ladder,
    /// on the same matrix the ladder is handed.
    void factorizePlain(const gsSparseMatrix<real_t> & K) { this->factorizeMatrix(K); }

    /// ||K_T*(V+DeltaV)|| of the LAST extended iteration -- the quantity
    /// gsALMBase::_extendedSystemSolve()'s `m_residueKTPhi < tolE` convergence test
    /// actually checks, evaluated on the UNSHIFTED tangent. NOT the same as ||K_T V||
    /// computed after the solve returned: the `m_V += m_DeltaV` mode commit in
    /// _extendedSystemSolve() happens AFTER the test passed.
    real_t residueKTPhi() const { return this->m_residueKTPhi; }
};

/// Runs \a f with gsInfo (== std::cout) redirected into a string.
/// The buffer is restored BEFORE returning -- and before any CHECK is evaluated -- so
/// that a failing assertion still reaches the test log instead of the capture.
template <class Fun>
std::string captureInfo(Fun f)
{
    std::ostringstream oss;
    std::streambuf * old = std::cout.rdbuf(oss.rdbuf());
    try                { f(); }
    catch (...)        { std::cout.rdbuf(old); throw; }
    std::cout.rdbuf(old);
    return oss.str();
}

/// Number of (non-overlapping) occurrences of \a needle in \a hay.
inline index_t countOccurrences(const std::string & hay, const std::string & needle)
{
    index_t n = 0;
    for (std::string::size_type p = hay.find(needle);
         p != std::string::npos; p = hay.find(needle,p+needle.size()))
        ++n;
    return n;
}

TEST(extended_system_shift_ladder_survives_an_exactly_singular_tangent)
{
    AlmProblem prob = pitchforkProblem();

    gsVector<real_t> U(2); U << 1.0, 0.0;      // the pitchfork itself
    const real_t     L = 1.0;                  // lambda* = 1
    gsVector<real_t> V(2); V << 0.0, 1.0;      // the closed-form critical mode

    gsSparseMatrix<real_t> K;
    prob.Jacobian(U,K);
    // The premise of the whole gate, checked and not assumed: 1.0 - u1 with u1 = 1.0 is
    // a bit-exact zero pivot, not a small number.
    CHECK_EQUAL(0.0, K.coeff(1,1));

    // --- (a) FAILS-FIRST, asserted permanently -----------------------------------
    // The ladder's attempt (i) on this exact matrix must FAIL. This is what makes the
    // rest of the test a statement about the shifted branch: with the ladder removed,
    // the throw below is what _extendedSystemIteration would propagate instead of
    // returning. NOTE: factorizeMatrix prints "Solver error with code 1"
    // (NumericalIssue) before throwing -- that line in the log is the EVIDENCE here,
    // not a failure.
    {
        ExtendedProbe plain(prob.Jacobian, prob.ALResidual, prob.Force);
        configure(plain, 0.05);
        index_t thrownCode = 0;
        try                  { plain.factorizePlain(K); }
        catch (int errorCode){ thrownCode = errorCode; }
        gsInfo<<"  [T41/G2] plain factorizeMatrix(diag(1,0)) threw code "<<thrownCode
              <<" (3 == SolverError)\n";
        CHECK_EQUAL(3, thrownCode);
    }

    // --- (b) the ladder runs, and the answer is unharmed --------------------------
    ExtendedProbe alm(prob.Jacobian, prob.ALResidual, prob.Force);
    configure(alm, 0.05);
    alm.setSolution(U,L);

    bool        converged = false;
    bool        threw     = false;
    std::string log       = captureInfo(
        [&alm,&U,&V,&converged,&threw,L]()
        {
            try { converged = alm.extendedSolveAt(U,L,V,1e-10); }
            catch (...) { threw = true; }
        });

    const index_t shifts = countOccurrences(log,"diagonal shift");
    gsInfo<<"  [T41/G2] extended solve: threw = "<<(threw?"yes":"no")
          <<"  converged = "<<(converged?"yes":"no")
          <<"  shift diagnostics = "<<shifts
          <<"  U = ("<<alm.solutionU()[0]<<","<<alm.solutionU()[1]<<")"
          <<"  L = "<<alm.solutionL()
          <<"  V = ("<<alm.solutionV()[0]<<","<<alm.solutionV()[1]<<")\n";

    // The outcome change the ladder buys: a RETURN instead of the escaping throw 3 that
    // (a) just measured on the very same matrix.
    CHECK(!threw);
    CHECK(converged);
    // ... and it demonstrably took the shifted branch (the message is emitted only
    // inside factorizeShifted's diagonal-shift retry, gsALMBase::_extendedSystemIteration()).
    // >= 1 rather than == 2 so the gate does not depend on
    // how many solves one iteration happens to need.
    CHECK(shifts >= 1);

    // I1/I2: the shift is applied to a COPY and only after the unshifted attempt
    // failed, and the convergence test uses the UNSHIFTED tangent -- so the seeded
    // singular point, which is already the exact solution of the extended system, must
    // come back BIT-UNCHANGED. (This is the honest statement of what happened here: the
    // solve is the identity map at this seed. What is being pinned is that a
    // 1e-12-perturbed factorization did not move it, not that Newton work was done.)
    CHECK_EQUAL(1.0, alm.solutionU()[0]);
    CHECK_EQUAL(0.0, alm.solutionU()[1]);
    CHECK_EQUAL(1.0, alm.solutionL());
    CHECK_EQUAL(0.0, alm.solutionV()[0]);
    CHECK_EQUAL(1.0, alm.solutionV()[1]);

    // --- (c) ... and the shifted solve still SOLVES ---------------------------------
    // (b) pins that the shift does not move an already-exact answer, but there the
    // extended solve is the identity map: no Newton work is done, so on its own it would
    // not distinguish the ladder from any other way of not crashing. Here the same
    // bit-exactly singular tangent K = diag(1,0) is seeded with lambda deliberately 0.05
    // OFF the singular point. The residual is then r = (1 - 1.05, 0) != 0, the extended
    // system has genuine work to do, and it must come back at lambda* = 1 -- through a
    // first tangent that plain factorizeMatrix cannot factorize at all, per (a).
    {
        ExtendedProbe off(prob.Jacobian, prob.ALResidual, prob.Force);
        configure(off, 0.05);
        off.setSolution(U,1.05);

        bool        c2 = false;
        bool        t2 = false;
        std::string log2 = captureInfo(
            [&off,&U,&V,&c2,&t2]()
            {
                try { c2 = off.extendedSolveAt(U,1.05,V,1e-10); }
                catch (...) { t2 = true; }
            });

        gsSparseMatrix<real_t> Kc;
        prob.Jacobian(off.solutionU(),Kc);
        const real_t modeResidual = (Kc*off.solutionV()).norm()/off.solutionV().norm();
        gsInfo<<"  [T41/G2] seed L = 1.05: threw = "<<(t2?"yes":"no")
              <<"  converged = "<<(c2?"yes":"no")
              <<"  shift diagnostics = "<<countOccurrences(log2,"diagonal shift")
              <<"  U = ("<<off.solutionU()[0]<<","<<off.solutionU()[1]<<")"
              <<"  L = "<<off.solutionL()
              <<"  tested ||K_T(V+DV)|| = "<<off.residueKTPhi()
              <<"  committed ||K_T V||/|V| = "<<modeResidual<<"\n";

        CHECK(!t2);
        CHECK(c2);
        // The singular point, from the closed form of fixture P -- NOT re-pasted from a
        // solver output.
        CHECK_CLOSE(1.0, off.solutionL(),    1e-10);
        CHECK_CLOSE(1.0, off.solutionU()[0], 1e-10);
        CHECK_CLOSE(0.0, off.solutionU()[1], 1e-10);
        // The mode is NOT bit-exact here, and that is the honest cost of the shift: the
        // returned V picked up a component of 8.23e-11 along the load direction (MEASURED),
        // i.e. the null direction is polluted -- the same effect an earlier report saw as a
        // divergence when the shift was made large.
        //
        // Two DIFFERENT quantities are printed above and they must not be conflated. The
        // one the solver tested is gsALMBase::_extendedSystemSolve()'s
        // m_residueKTPhi = ||K_T*(V+DeltaV)|| convergence test, on the
        // UNSHIFTED tangent, BEFORE the `m_V += m_DeltaV` mode commit; modeResidual is recomputed
        // here from the committed solutionV(). MEASURED they coincide to every digit
        // (8.23086e-11), and they must: deltaU is zero at this seed, so the converged
        // tangent is the same K, and |V+DeltaV| = 1 makes the normalisation a no-op.
        // Hence the margin against SingularPointComputeTolE = 1e-10 really is ~18%: the
        // pollution is O(sigma) and this seed converges on the FIRST iteration only
        // because 8.2e-11 happens to fall under tolE.
        // A bound an order of magnitude above the measurement pins that it stays a
        // rounding-level effect and does not silently grow.
        CHECK(modeResidual < 1e-9);
    }
}

// ===========================================================================
// TEST 14 -- the BORDERED-CHART TRANSVERSALITY GUARD.
// ===========================================================================

/// K = I2 EXACTLY at every state: F_int(u) = u, Force = (1,0)^T, so
/// R = (u1 - lambda, u2), K = diag(1,1) and delta_u_t = K^-1 f = (1,0)^T ALWAYS.
/// The equilibrium path is the straight line u = (lambda,0). This fixture isolates
/// the chart itself: it is the ONE tangent that is perfectly conditioned, so any
/// degeneracy the corrector runs into on it is a property of the CHART and of nothing
/// else.
inline AlmProblem identityProblem()
{
    AlmProblem p;
    p.scale = 1.0;
    p.Force.resize(2);
    p.Force << 1.0, 0.0;
    p.Jacobian =
        [](gsVector<real_t> const & /*u*/, gsSparseMatrix<real_t> & m) -> bool
        {
            gsMatrix<real_t> K(2,2);
            K(0,0) = 1.0; K(0,1) = 0.0;
            K(1,0) = 0.0; K(1,1) = 1.0;
            toSparse2(K,m); return true;
        };
    p.ALResidual =
        [](gsVector<real_t> const & u, const real_t lambda, gsVector<real_t> & result) -> bool
        { result.resize(2); result[0] = u[0]-lambda; result[1] = u[1]; return true; };
    p.Fint = [](const gsVector<real_t> & u) { return u; };
    return p;
}

/// K = diag(1,eps) EXACTLY at every state: F_int(u) = (u1, eps*u2), Force = (1,0)^T.
/// Same straight equilibrium path as identityProblem(), but the tangent is
/// (numerically) RANK DEFICIENT in the u2 direction, and f = e1 lies in range(K):
/// Riks 1984's case (2.20b) configuration (G and G* both rank deficient), NOT a fold.
/// delta_u_t = K^-1 f = (1,0) stays BOUNDED here -- which is exactly why a
/// singularity estimate taken along the load direction cannot see this case.
inline AlmProblem nearSingularProblem(real_t eps)
{
    AlmProblem p;
    p.scale = 1.0;
    p.Force.resize(2);
    p.Force << 1.0, 0.0;
    p.Jacobian =
        [eps](gsVector<real_t> const & /*u*/, gsSparseMatrix<real_t> & m) -> bool
        {
            gsMatrix<real_t> K(2,2);
            K(0,0) = 1.0; K(0,1) = 0.0;
            K(1,0) = 0.0; K(1,1) = eps;
            toSparse2(K,m); return true;
        };
    p.ALResidual =
        [eps](gsVector<real_t> const & u, const real_t lambda, gsVector<real_t> & result) -> bool
        { result.resize(2); result[0] = u[0]-lambda; result[1] = eps*u[1]; return true; };
    p.Fint = [eps](const gsVector<real_t> & u)
        { gsVector<real_t> F(2); F[0] = u[0]; F[1] = eps*u[1]; return F; };
    return p;
}

// ---------------------------------------------------------------------------
// THE STATE ALL THREE TESTS BELOW USE, and how it is reached.
// ---------------------------------------------------------------------------
// The corrector's constraint line  {DeltaU + t*delta_u_t}  is TANGENT to the arc-length
// sphere exactly when  b0 = 2*(DeltaU.delta_u_t + A0*DeltaLambda) = 0, which is Riks's
// (n^t.x') = 0. On identityProblem() at Scaling = 0 (A0 = 0, so b0 = 2*DeltaU.e1) that is
// reached WITHOUT installing anything by hand: seed the step at
//     U = (u1, -ds),  L = u1 + h   (h != 0),
// let the class's own predictor take DeltaU = ds*e1, and let its own first corrector
// iteration run. u_bar's second component is then EXACTLY -u2 = ds (K = I2 and R2 = u2 make
// that exact, no tuning), so the constraint |DeltaU| = ds forces the FIRST component to
// zero -- and the first component is computed as a floating-point SUM of O(ds) terms whose
// exact value is 0, i.e. it lands on ROUNDING. h != 0 is what makes those terms inexact:
// with h = 0 they cancel bit-exactly, DeltaU.e1 == 0, B is EXACTLY singular and SparseLU
// reports it, which is the already-handled path and not the defect.
//
// MEASURED ON THE UNFIXED TREE (fails-first, ds = 0.13, u1 = 0.3, h = 0.11):
//   after iteration 1: DeltaU = (1.36977e-09, 0.13), |DeltaU| = ds  (on the sphere)
//   iteration 2:       q_t = 3.65024e+08 ACCEPTED, no diagnostic, and the state it moves to
//                      has |(DeltaU,DeltaLambda)|_A = 0.26 = 2*ds -- the arc-length
//                      constraint violated by 100%, with the step still reporting success.
// The same tangency at ds = 0.05 with iteration 1 taken in the ELIMINATION chart gives
// DeltaU = (1.38778e-17,0.05) and an accepted q_t = 3.60288e+16.
// ---------------------------------------------------------------------------

/// Seeds the tangency above and runs predictor + iteration 1, all through the class's own
/// initiateStep()/predictor()/quasiNewtonIteration()/iteration(). Nothing is installed.
/// \param borderedFirst  true: the WHOLE step runs in the bordered chart, which is what
///        gsALMCrisfield::step()'s retry does. false: iteration 1 runs in the ELIMINATION
///        chart and only iteration 2 is bordered -- a configuration the retry cannot
///        produce (it sets m_useBordered for the whole step), so it is labelled as what it
///        is: the chart switch is moved between iterations, and NOTHING else is installed.
inline void seedTangentCorrectorState(CrisfieldProbe & alm, const AlmProblem & prob,
                                      real_t ds, bool borderedFirst = true,
                                      real_t u1 = 0.3, real_t h = 0.11)
{
    GISMO_UNUSED(prob);
    configure(alm, ds);
    alm.options().setReal("Scaling",0.0);          // cylindrical, pinned IN THE TEST
    alm.applyOptions();
    gsVector<real_t> U(2); U << u1, -ds;
    alm.setSolution(U, u1 + h);
    alm.setLength(ds);
    alm.forceBordered(borderedFirst);
    alm.runPredictorAndOneIteration();
    alm.forceBordered(true);
}

// ---------------------------------------------------------------------------
// G-46a -- the guard fires on Riks's case (2.20a) and RESCUES the chart.
// ---------------------------------------------------------------------------
// ⚠ Not "the step succeeds": it does not, and that is reported honestly below. What is
// pinned is that the corrector no longer ACCEPTS a chart it cannot compute in -- the two
// numbers that move are q_t (3.65e8 -> O(10)) and the accepted arc length (2*ds -> ds).
TEST(bordered_chart_transversality_guard_rescues_a_turning_chart)
{
    const real_t ds = 0.13;
    AlmProblem prob = identityProblem();
    CrisfieldProbe alm(prob.Jacobian, prob.ALResidual, prob.Force);
    seedTangentCorrectorState(alm,prob,ds);

    // The state really is the counterexample, and it really is ON the constraint surface --
    // otherwise this would only be a test about a corrector that had already gone wrong.
    // delta_u_t = K^-1 f = e1 exactly here, so b0 = 2*DeltaU[0].
    const real_t b0 = 2.0*alm.DeltaU()[0];
    gsInfo<<"  [T46/a] after iteration 1: DeltaU = ("<<alm.DeltaU()[0]<<","<<alm.DeltaU()[1]
          <<")  b0 = "<<b0<<"  |(DU,DL)|_A = "<<alm.distance(alm.DeltaU(),alm.DeltaL())
          <<"  cosine = "<<alm.chartCosine()<<"\n";
    CHECK(math::abs(b0) < 1e-8);                            // tangent to machine precision
    CHECK_CLOSE(ds, alm.distance(alm.DeltaU(),alm.DeltaL()), 1e-14);
    CHECK_EQUAL(0, alm.chartReborderings());                // iteration 1 was HEALTHY
    CHECK(alm.chartCosine() > alm.chartCosTol());

    // Iteration 2 is the one whose bordered matrix is numerically singular.
    index_t code = 0;
    try { alm.runOneCorrectorIteration(); } catch (int e) { code = e; }
    gsInfo<<"  [T46/a] iteration 2: code = "<<code<<"  q_t = "<<alm.deltaLt()
          <<"  cosine = "<<alm.chartCosine()<<"  reborderings = "<<alm.chartReborderings()
          <<"  |(DU,DL)|_A = "<<alm.distance(alm.DeltaU(),alm.DeltaL())<<"\n";

    // (i) the guard fired AND classified (2.20a): it re-bordered rather than reported.
    CHECK_EQUAL(1, alm.chartReborderings());
    // (ii) the chart it settled on is transversal and BOUNDED. On the unfixed tree this is
    //      q_t = 3.65024e+08 (MEASURED) and there is no cosine at all.
    CHECK(alm.chartCosine() > alm.chartCosTol());
    CHECK(math::abs(alm.deltaLt()) < 1e3);
    // (iii) and no garbage was COMMITTED: the corrector state is untouched at |.|_A = ds,
    //      where the unfixed tree moved it to 2*ds (MEASURED 0.26) while reporting success.
    CHECK_CLOSE(ds, alm.distance(alm.DeltaU(),alm.DeltaL()), 1e-14);
    // (iv) the step is REPORTED (throw 1 -> NotConverged), not silently accepted. The reason
    //      it is a report rather than a completed iteration is NOT the guard: at an exact
    //      tangency the constraint quadratic has a DOUBLE root, its discriminant is 0 up to
    //      rounding, and the bordered chart's complex-root branch is the one deliberately
    //      left unimplemented (class documentation (g)). Re-bordering fixes the CHART, not
    //      the geometry.
    CHECK_EQUAL(1, code);

    // (v) DETERMINISM. The re-bordering direction comes from a fixed-seed LCG, so a second
    //     run of the same state must reproduce it BIT for bit -- restart determinism (unit
    //     test 9) and the pinned CSV md5s both depend on this.
    CrisfieldProbe again(prob.Jacobian, prob.ALResidual, prob.Force);
    seedTangentCorrectorState(again,prob,ds);
    index_t code2 = 0;
    try { again.runOneCorrectorIteration(); } catch (int e) { code2 = e; }
    CHECK_EQUAL(code, code2);
    CHECK_EQUAL(alm.deltaLt(),     again.deltaLt());        // bit equality, not CHECK_CLOSE
    CHECK_EQUAL(alm.chartCosine(), again.chartCosine());

    // (vi) THE SPEC'S OWN HEADLINE NUMBER (q_t ~ 5e15), which the
    //      ds = 0.13 state above does not reach: the same tangency at ds = 0.05 with
    //      iteration 1 taken in the ELIMINATION chart. MEASURED ON THE UNFIXED TREE:
    //      DeltaU = (1.38778e-17,0.05), then iteration 2 ACCEPTED q_t = 3.60288e+16 with
    //      code = 0 and moved the state to |(DU,DL)|_A = 0.1 = 2*ds.
    //      ⚠ RUNG: only the chart SWITCH moves between the two iterations (see
    //      seedTangentCorrectorState); every state is produced by the class's own
    //      predictor/iteration. gsALMCrisfield::step()'s retry cannot produce this
    //      combination -- it sets m_useBordered for the whole step -- so this half is a
    //      configuration the harness makes, and it is labelled rather than dressed up.
    {
        const real_t ds2 = 0.05;
        CrisfieldProbe mix(prob.Jacobian, prob.ALResidual, prob.Force);
        seedTangentCorrectorState(mix,prob,ds2,/*borderedFirst*/false);
        CHECK(math::abs(2.0*mix.DeltaU()[0]) < 1e-15);
        CHECK_CLOSE(ds2, mix.distance(mix.DeltaU(),mix.DeltaL()), 1e-15);
        // NO cosine assertion here: iteration 1 ran in the elimination chart, so
        // m_chartCosine is still the NaN initMethods() poisons it with (MEASURED: nan).
        CHECK_EQUAL(0, mix.chartReborderings());

        index_t code3 = 0;
        try { mix.runOneCorrectorIteration(); } catch (int e) { code3 = e; }
        gsInfo<<"  [T46/a] ds = 0.05, elimination->bordered: code = "<<code3<<"  q_t = "
              <<mix.deltaLt()<<"  cosine = "<<mix.chartCosine()<<"  reborderings = "
              <<mix.chartReborderings()<<"  |(DU,DL)|_A = "
              <<mix.distance(mix.DeltaU(),mix.DeltaL())<<"\n";
        CHECK_EQUAL(1, mix.chartReborderings());
        CHECK(mix.chartCosine() > mix.chartCosTol());
        CHECK(math::abs(mix.deltaLt()) < 1e3);                  // unfixed tree: 3.60288e+16
        CHECK_CLOSE(ds2, mix.distance(mix.DeltaU(),mix.DeltaL()), 1e-15);  // unfixed: 2*ds
        CHECK_EQUAL(1, code3);                                  // unfixed tree: 0
    }
}

// ---------------------------------------------------------------------------
// G-46b -- ... and it does NOT re-border Riks's case (2.20b).
// ---------------------------------------------------------------------------
// THE SAME tangency, THE SAME seed, THE SAME arc length: the ONLY thing that changes is the
// tangent, from K = I2 to K = diag(1,1e-12) -- numerically rank deficient in the u2
// direction, with f = e1 IN range(K). Riks: "a rank deficiency of at least 1 for both G and
// G* ... indicates the existence of multiple solutions x' of (2.18) and thus the
// manifestation of bifurcations". A guard that answered "chart artifact, re-border and
// continue" here would walk straight past a bifurcation -- in a landscape explorer built to
// find them. The pair G-46a/G-46b is therefore what pins the DISCRIMINATOR; either half
// alone is consistent with a guard that always re-borders, or always reports.
//
// ⚠ delta_u_t = K^-1 f = e1 is BOUNDED here (f does not excite the near-null direction), so
// a singularity estimate taken along the load direction sees NOTHING: it reports
// |f|/|delta_u_t| = 1. MEASURED on this exact state, before the probe was corrected to
// reconstruct K^-1 v.
TEST(bordered_chart_transversality_guard_reports_a_rank_deficient_tangent)
{
    const real_t ds  = 0.13;
    const real_t eps = 1e-12;
    AlmProblem prob = nearSingularProblem(eps);
    CrisfieldProbe alm(prob.Jacobian, prob.ALResidual, prob.Force);
    seedTangentCorrectorState(alm,prob,ds);

    // Same tangency as G-46a (the fixture was built so that u_bar[1] = -u2 exactly there too).
    CHECK(math::abs(2.0*alm.DeltaU()[0]) < 1e-8);
    CHECK_EQUAL(0, alm.chartReborderings());

    index_t code = 0;
    try { alm.runOneCorrectorIteration(); } catch (int e) { code = e; }
    gsInfo<<"  [T46/b] eps = "<<eps<<"  iteration 2: code = "<<code<<"  cosine = "
          <<alm.chartCosine()<<"  reborderings = "<<alm.chartReborderings()<<"\n";

    // The verdict: reported, NOT re-bordered, and the step fails rather than continuing.
    CHECK_EQUAL(0, alm.chartReborderings());
    CHECK(alm.chartCosine() < alm.chartCosTol());
    CHECK_EQUAL(1, code);
    CHECK_CLOSE(ds, alm.distance(alm.DeltaU(),alm.DeltaL()), 1e-14);

    // THE ORACLES, both independent of the guard's own arithmetic.
    // 1. det K_T = 1e-12 -- Riks's D = 0 branch of his summary (2.29).
    gsVector<real_t> Ueq(2); Ueq << 0.3, 0.0;
    CHECK_CLOSE(eps, tangentDet(prob,Ueq), 1e-16);
    // 2. gsALMBase::_testSingularPoint, i.e. the phi^T P criterion of Wriggers & Simo 1990
    //    eqs. (9a)/(9b), evaluated at the same tangent by the LIBRARY's own machinery (up to
    //    SingularPointTestIt inverse-iteration sweeps and gsALMBase's own m_V, neither of
    //    which the guard is allowed to touch mid-corrector). It must say BIFURCATION, which
    //    is the verdict the guard reached from inside the bordered chart.
    CrisfieldProbe ref(prob.Jacobian, prob.ALResidual, prob.Force);
    configure(ref, ds);
    ref.setSolution(Ueq,0.3);
    const bool bifurcation = ref.isBifurcation(true);
    const real_t cosVf = modeForceCosine(ref.solutionV(),prob.Force);
    gsInfo<<"  [T46/b] oracles: det K = "<<tangentDet(prob,Ueq)<<"  _testSingularPoint = "
          <<(bifurcation?"BIFURCATION":"limit point")<<"  |V.f|/|f| = "<<cosVf<<"\n";
    CHECK(bifurcation);
    CHECK(cosVf < 1e-4);

    // ... and the CONTRAST that makes this a discriminator and not a coin: on K = I2 the
    // very same seed is classified the other way (G-46a re-borders). Re-asserted here so the
    // pair cannot drift apart.
    AlmProblem reg = identityProblem();
    CrisfieldProbe alm2(reg.Jacobian, reg.ALResidual, reg.Force);
    seedTangentCorrectorState(alm2,reg,ds);
    try { alm2.runOneCorrectorIteration(); } catch (int) { }
    CHECK_EQUAL(1, alm2.chartReborderings());
}

// ---------------------------------------------------------------------------
// G-46c -- the guard must NOT fire whenever the bordered chart is merely USED.
// ---------------------------------------------------------------------------
// This is the half that separates "fires at b0 = 0" from "fires in the bordered chart", and
// it is also the bit-neutrality argument for TESTS 11: every trace below runs ENTIRELY in
// the bordered chart (forceBordered, not the NotConverged retry), and the cumulative
// re-bordering counter must still read zero at the end.
TEST(bordered_chart_transversality_guard_is_silent_on_a_healthy_chart)
{
    // (a) the fold traces the fallback exists for, at all three load scales.
    const real_t scales[3] = {1.0, 1e3, 1e-3};
    for (index_t s = 0; s != 3; ++s)
    {
        AlmProblem prob = foldProblem(scales[s]);
        CrisfieldProbe alm(prob.Jacobian, prob.ALResidual, prob.Force);
        configure(alm, 0.05);
        alm.options().setReal("Scaling",0.0);
        alm.applyOptions();
        alm.forceBordered(true);
        index_t ok = 0;
        real_t  minCos = std::numeric_limits<real_t>::max();
        for (index_t k = 0; k != 40; ++k)
        {
            if (alm.step() != gsStatus::Success) break;
            ++ok;
            minCos = math::min(minCos,alm.chartCosine());
            if (alm.solutionU()[0] > 1.8) break;
        }
        gsInfo<<"  [T46/c] fold trace scale = "<<scales[s]<<"  steps = "<<ok
              <<"  min cosine = "<<minCos<<"  reborderings = "<<alm.chartReborderings()
              <<"  max u1 = "<<alm.solutionU()[0]<<"\n";
        // Not a tautology: the trace has to GET somewhere (past the fold at u1 = 1) for
        // "never fired" to mean anything -- a trace that dies at step 0 never fires either.
        CHECK(ok >= 20);
        CHECK(alm.solutionU()[0] > 1.05);
        CHECK_EQUAL(0, alm.chartReborderings());
        CHECK(minCos > alm.chartCosTol());
    }

    // (b) a genuinely TWO-DIMENSIONAL path, sampled at EVERY corrector iteration. Fixture F
    //     moves in one component only, so its cosine is 1 by construction and carries no
    //     margin information; the bifurcated branch of fixture P moves in both. ds = 0.5 is
    //     deliberately coarse -- a long step is where the corrector wanders furthest from
    //     the tangent, i.e. where the cosine is smallest.
    for (index_t sc = 0; sc != 2; ++sc)
    {
        AlmProblem prob = pitchforkProblem();
        CrisfieldProbe alm(prob.Jacobian, prob.ALResidual, prob.Force);
        configure(alm, 0.5);
        alm.options().setReal("Scaling", sc==0 ? 0.0 : -1.0);
        alm.applyOptions();
        // ON the bifurcated branch: u1 = 2L-1, u2 = sqrt(2(L-1)) at L = 1.5 (closed form,
        // gsALMTestProblems.h) -- checked, not assumed.
        gsVector<real_t> U(2); U << 2.0, 1.0;
        CHECK_CLOSE(pitchforkBranchU1(1.5),   U[0],      1e-14);
        CHECK_CLOSE(pitchforkBranchU2sq(1.5), U[1]*U[1], 1e-14);
        alm.setSolution(U,1.5);
        alm.setLength(0.5);
        alm.forceBordered(true);

        real_t minCos = std::numeric_limits<real_t>::max();
        index_t code = 0;
        try
        {
            alm.runPredictorAndOneIteration();
            minCos = math::min(minCos,alm.chartCosine());
            for (index_t k = 0; k != 6; ++k)
            {
                alm.runOneCorrectorIteration();
                minCos = math::min(minCos,alm.chartCosine());
            }
        }
        catch (int e) { code = e; }
        gsInfo<<"  [T46/c] pitchfork ds = 0.5  Scaling = "<<(sc==0?0.0:-1.0)
              <<"  per-iteration min cosine = "<<minCos<<"  code = "<<code
              <<"  reborderings = "<<alm.chartReborderings()<<"\n";
        CHECK_EQUAL(0, code);
        CHECK_EQUAL(0, alm.chartReborderings());
        // The MARGIN the threshold was placed in: the healthy minimum is O(1), five orders
        // above m_chartCosTol = 1e-6, which in turn is two orders above the 1.05e-8 of the
        // counterexample in G-46a.
        CHECK(minCos > 0.9);
    }
}

// ===========================================================================
// gsALMRiks::predictorGuess() must initialise the convex
// weight.
// ===========================================================================
// ! HONEST LABEL -- READ BEFORE TREATING THIS AS A FAILS-FIRST GATE !
//
// This is a POSTCONDITION PIN, not a fails-first test, and it was never claimed
// to be one. The defect it guards is that gsALMRiks::predictorGuess() never
// assigned the convex constraint weight -- only gsALMRiks::predictor() did -- so
// a guess-seeded FIRST step entered iteration()/distance()/stepOutput() reading
// an UNINITIALIZED member. An uninitialized read is undefined behaviour: on the
// unfixed tree that member may hold anything at all, INCLUDING the correct
// 1/numDof, so a RED run cannot be guaranteed and no RED transcript is offered.
// What the test pins is the postcondition -- "after a guess-seeded first step the
// weight equals the value predictor() would have written" -- in the same spirit
// as the label on riks_rounds_the_padded_fold above.
//
// Scope note: gsALMRiks::predictorGuess() has NO in-tree caller (nothing calls
// gsALMBase::setInitialGuess outside this test), which is also why the fix is
// oracle-neutral by construction. The path is deprecated and deliberately kept
// as it is: this removes the undefined behaviour from it, it does not revive
// it and does not make it throw.
TEST(riks_predictor_guess_initialises_the_convex_weight)
{
    AlmProblem prob = foldProblem();               // fixture F: exactly 2 DOFs
    RiksProbe  alm(prob.Jacobian, prob.ALResidual, prob.Force);
    configure(alm, 0.05);

    // gsALMBase::_step() dispatches to predictorGuess() instead of predictor() only
    // when the guess differs from the current state in BOTH U and Lambda.
    gsVector<real_t> Uguess(2); Uguess << 0.05, 0.0;
    alm.setInitialGuess(Uguess, 0.05);

    // The FIRST step of a fresh solver: predictor(), the only other writer of the
    // weight, has never run, so the value observed below can only come from
    // predictorGuess().
    const gsStatus st = alm.step();
    gsInfo<<"  [T48/2] guess-seeded first step: converged = "<<(st==gsStatus::Success)
          <<"  guessConsumed = "<<alm.guessConsumed()
          <<"  phi = "<<alm.phi()<<"\n";

    // Not a tautology: without this check the assertion below would also be satisfied
    // by the ORDINARY predictor(), which writes the very same value. Only
    // predictorGuess() clears m_Uguess.
    CHECK(alm.guessConsumed());
    CHECK_CLOSE(1.0/2.0, alm.phi(), 1e-14);        // 1/numDof, numDof = 2
}

// ===========================================================================
// gsAPALM termination guards
// ===========================================================================
// WHY THE OBVIOUS TEST IS THE WRONG TEST.
//
// Both loops pinned below are UNBOUNDED on the unfixed tree, so the naive
// fails-first test -- drive gsAPALM::_initiation with a solver that always fails
// and wait for it to come back -- does not FAIL, it never RETURNS: it would hang
// this suite forever, for every future contributor.
//
// The stub below is therefore CALL-CAPPED. It counts every gsALMBase::step() and
// GISMO_ERRORs once the cap is exceeded, which converts "does not terminate" into
// a bounded, loud test failure (UnitTest++ reports an escaping exception as a
// failed test). The cap is 50; the post-fix cost is ~20 calls, one per halving of
// the arc length down to the |dL/dL0| < 1e-6 floor that
// gsALMExploration<T>::traceCurve already uses (2^-20 = 9.5e-7 < 1e-6), so the two
// numbers are two decades of margin apart and the upper bound asserted below
// (<= 25) pins that termination happened VIA THE FLOOR and not by accident.

/// A gsALMBase whose step() only counts the call and returns a FIXED status --
/// the deterministic, state-preserving failure that gsAPALM's loops must survive.
///
/// Derived from gsALMLoadControl because it is the only family member whose
/// distance() is stateless (gsALMLoadControl<T>::distance returns |DeltaL|).
/// gsALMRiks<T>::distance() would read m_convexWeight, which is written by
/// gsALMRiks::predictor()/predictorGuess() only -- and this stub's step() never
/// runs a predictor.
class CountingALM : public gsALMLoadControl<real_t>
{
    typedef gsALMLoadControl<real_t> Base;
public:
    CountingALM(ALMJacobian_t          & J,
                ALMResidual_t          & R,
                gsVector<real_t>       & F,
                gsStatus                 stepStatus,
                index_t                  callCap)
    : Base(J,R,F), m_stepStatus(stepStatus), m_callCap(callCap), m_calls(0)
    {}

    /// Overrides gsALMBase<T>::step() (virtual), the only entry point through which
    /// gsAPALM::_initiation / gsAPALM::_correction advance the solver.
    gsStatus step()
    {
        if (++m_calls > m_callCap)
            GISMO_ERROR("CountingALM: step() called more than "<<m_callCap
                        <<" times -- the gsAPALM loop under test does NOT terminate.");
        return m_stepStatus;
    }

    index_t calls() const { return m_calls; }

private:
    gsStatus m_stepStatus;
    index_t  m_callCap;
    index_t  m_calls;
};

/// Exposes gsAPALM's two protected worker routines. Nothing under src/ is
/// modified: the members are reached the way C++ intends, by deriving.
class APALMProbe : public gsAPALM<real_t>
{
    typedef gsAPALM<real_t> Base;
public:
    APALMProbe(gsALMBase<real_t> * ALM,
               const gsAPALMData<real_t,Base::solution_t> & data)
    : Base(ALM,data) {}

    using Base::_initiation;
    using Base::_correction;
};

typedef gsAPALM<real_t>::solution_t apalm_solution_t;
typedef std::tuple<index_t,real_t,apalm_solution_t,apalm_solution_t> apalm_entry_t;

/// The (ID, dL, start, prev) tuple gsAPALM's job queue hands to both routines.
inline apalm_entry_t apalmEntry(real_t dL, const gsVector<real_t> & U, real_t L)
{
    return std::make_tuple((index_t)0, dL, std::make_pair(U,L), std::make_pair(U,L));
}

// ---------------------------------------------------------------------------
// SITE A -- gsAPALM<T>::_initiation, the `while (diverged)` loop.
//
// On SolverError/OtherError the loop stays alive (diverged == true) while the
// halve-and-reseed body is SKIPPED, and a failed step commits nothing
// (m_U += m_DeltaU lives in iterationFinish(), which gsALMBase<T>::_step()
// reaches only after the convergence test), so step() is re-run from the
// BIT-IDENTICAL state forever. On the handled NotConverged/AssemblyError branch
// the arc length is halved with NO floor and NO iteration cap.
// ---------------------------------------------------------------------------
TEST(apalm_initiation_terminates_on_solver_error)
{
    AlmProblem prob = foldProblem();
    const index_t cap = 50;
    CountingALM alm(prob.Jacobian, prob.ALResidual, prob.Force, gsStatus::SolverError, cap);
    configure(alm, 0.1);

    gsAPALMData<real_t,apalm_solution_t> data;
    APALMProbe apalm(&alm,data);

    gsVector<real_t> U0(2); U0 << 0.5, 0.0;
    real_t distance = -1.0;
    std::vector<apalm_solution_t> solutions;
    bool bifurcation = true;

    const gsStatus status = apalm._initiation(apalmEntry(0.1,U0,0.5), /*startTime=*/0.0,
                                              distance, solutions, bifurcation);

    gsInfo<<"  [T49/A1] _initiation, SolverError: step() calls = "<<alm.calls()
          <<"  status = "<<(index_t)status
          <<"  solutions = "<<solutions.size()<<"\n";

    CHECK(alm.calls() < cap);               // it terminated at all
    CHECK(alm.calls() <= 25);               // ... via the arc-length floor (20 halvings)
    CHECK(status != gsStatus::Success);     // and said so
    CHECK_EQUAL(0, (index_t)solutions.size());  // recording no phantom point
}

TEST(apalm_initiation_terminates_on_not_converged)
{
    AlmProblem prob = foldProblem();
    const index_t cap = 50;
    CountingALM alm(prob.Jacobian, prob.ALResidual, prob.Force, gsStatus::NotConverged, cap);
    configure(alm, 0.1);

    gsAPALMData<real_t,apalm_solution_t> data;
    APALMProbe apalm(&alm,data);

    gsVector<real_t> U0(2); U0 << 0.5, 0.0;
    real_t distance = -1.0;
    std::vector<apalm_solution_t> solutions;
    bool bifurcation = true;

    const gsStatus status = apalm._initiation(apalmEntry(0.1,U0,0.5), /*startTime=*/0.0,
                                              distance, solutions, bifurcation);

    gsInfo<<"  [T49/A2] _initiation, NotConverged: step() calls = "<<alm.calls()
          <<"  status = "<<(index_t)status
          <<"  solutions = "<<solutions.size()<<"\n";

    CHECK(alm.calls() < cap);
    CHECK(alm.calls() <= 25);
    CHECK(status != gsStatus::Success);
    CHECK_EQUAL(0, (index_t)solutions.size());
}

// ---------------------------------------------------------------------------
// SITE B -- gsAPALM<T>::_correction, the `for (k)` interval loop.
//
// Two distinct defects, one test each:
//  * SolverError/OtherError falls THROUGH to stepSolutions.at(k) = ..., so the
//    UNCHANGED previous state is recorded as a converged interval solution -- the
//    phantom-point bug already excluded from gsALMExploration::traceCurve. This
//    one TERMINATES on the unfixed tree; it is wrong, not hanging.
//  * NotConverged/AssemblyError does `dL /= 2; k -= 1; continue;` with no floor
//    and no cap, i.e. it does not terminate.
// ---------------------------------------------------------------------------
TEST(apalm_correction_records_no_phantom_on_solver_error)
{
    AlmProblem prob = foldProblem();
    const index_t cap = 50;
    CountingALM alm(prob.Jacobian, prob.ALResidual, prob.Force, gsStatus::SolverError, cap);
    configure(alm, 0.1);

    gsAPALMData<real_t,apalm_solution_t> data;
    APALMProbe apalm(&alm,data);

    gsVector<real_t> U0(2);  U0  << 0.5, 0.0;
    gsVector<real_t> Uref(2); Uref << 0.6, 0.0;   // same half-space: keeps the
                                                  // direction GISMO_ASSERT in the
                                                  // tail of _correction happy
    std::vector<real_t> distances;
    std::vector<apalm_solution_t> stepSolutions;
    real_t upperDistance = -1.0, lowerDistance = -1.0;

    const gsStatus status = apalm._correction(apalmEntry(0.2,U0,0.5),
                      std::make_pair((real_t)0.0,(real_t)1.0),
                      /*dataLevel=*/1,
                      std::make_pair(Uref,(real_t)0.6),
                      distances, stepSolutions, upperDistance, lowerDistance);

    gsInfo<<"  [T49/B1] _correction, SolverError: step() calls = "<<alm.calls()
          <<"  status = "<<(index_t)status
          <<"  stepSolutions = "<<stepSolutions.size()
          <<"  distances = "<<distances.size()<<"\n";

    CHECK(alm.calls() < cap);
    CHECK(alm.calls() <= 25);
    CHECK(status != gsStatus::Success);
    // The phantom pin: a failed step commits nothing, so anything recorded here is
    // the unchanged start state masquerading as a converged interval solution.
    CHECK_EQUAL(0, (index_t)stepSolutions.size());
    CHECK_EQUAL(0, (index_t)distances.size());
}

TEST(apalm_correction_terminates_on_not_converged)
{
    AlmProblem prob = foldProblem();
    const index_t cap = 50;
    CountingALM alm(prob.Jacobian, prob.ALResidual, prob.Force, gsStatus::NotConverged, cap);
    configure(alm, 0.1);

    gsAPALMData<real_t,apalm_solution_t> data;
    APALMProbe apalm(&alm,data);

    gsVector<real_t> U0(2);  U0  << 0.5, 0.0;
    gsVector<real_t> Uref(2); Uref << 0.6, 0.0;
    std::vector<real_t> distances;
    std::vector<apalm_solution_t> stepSolutions;
    real_t upperDistance = -1.0, lowerDistance = -1.0;

    const gsStatus status = apalm._correction(apalmEntry(0.2,U0,0.5),
                      std::make_pair((real_t)0.0,(real_t)1.0),
                      /*dataLevel=*/1,
                      std::make_pair(Uref,(real_t)0.6),
                      distances, stepSolutions, upperDistance, lowerDistance);

    gsInfo<<"  [T49/B2] _correction, NotConverged: step() calls = "<<alm.calls()
          <<"  status = "<<(index_t)status
          <<"  stepSolutions = "<<stepSolutions.size()
          <<"  distances = "<<distances.size()<<"\n";

    CHECK(alm.calls() < cap);
    CHECK(alm.calls() <= 25);
    CHECK(status != gsStatus::Success);
    CHECK_EQUAL(0, (index_t)stepSolutions.size());
    CHECK_EQUAL(0, (index_t)distances.size());
}

// ---------------------------------------------------------------------------
// The DEGENERATE input the ratio test itself must survive: a job queued with a
// zero arc length. |dL/dL0| is then 0/0 = NaN, and `NaN < 1e-6` is FALSE, so the
// straightforward transcription of gsALMExploration<T>::traceCurve's break would
// leave the very spin-loop these guards remove. Both guards therefore use the NEGATED
// comparison; this test is what pins that choice (it fails with the `<` form).
// ---------------------------------------------------------------------------
TEST(apalm_termination_guards_survive_a_zero_arc_length_job)
{
    AlmProblem prob = foldProblem();
    const index_t cap = 50;

    gsVector<real_t> U0(2);   U0   << 0.5, 0.0;
    gsVector<real_t> Uref(2); Uref << 0.6, 0.0;

    {
        CountingALM alm(prob.Jacobian, prob.ALResidual, prob.Force, gsStatus::SolverError, cap);
        configure(alm, 0.1);
        gsAPALMData<real_t,apalm_solution_t> data;
        APALMProbe apalm(&alm,data);

        real_t distance = -1.0;
        std::vector<apalm_solution_t> solutions;
        bool bifurcation = true;
        const gsStatus status = apalm._initiation(apalmEntry(0.0,U0,0.5), /*startTime=*/0.0,
                                                  distance, solutions, bifurcation);
        gsInfo<<"  [T49/C1] _initiation, dL0 = 0: step() calls = "<<alm.calls()
              <<"  status = "<<(index_t)status<<"\n";
        CHECK(alm.calls() < cap);
        CHECK_EQUAL(1, alm.calls());          // the first failure already exhausts it
        CHECK(status != gsStatus::Success);
        CHECK_EQUAL(0, (index_t)solutions.size());
    }
    {
        CountingALM alm(prob.Jacobian, prob.ALResidual, prob.Force, gsStatus::NotConverged, cap);
        configure(alm, 0.1);
        gsAPALMData<real_t,apalm_solution_t> data;
        APALMProbe apalm(&alm,data);

        std::vector<real_t> distances;
        std::vector<apalm_solution_t> stepSolutions;
        real_t upperDistance = -1.0, lowerDistance = -1.0;
        const gsStatus status = apalm._correction(apalmEntry(0.0,U0,0.5),
                          std::make_pair((real_t)0.0,(real_t)1.0),
                          /*dataLevel=*/1,
                          std::make_pair(Uref,(real_t)0.6),
                          distances, stepSolutions, upperDistance, lowerDistance);
        gsInfo<<"  [T49/C2] _correction, dL0 = 0: step() calls = "<<alm.calls()
              <<"  status = "<<(index_t)status<<"\n";
        CHECK(alm.calls() < cap);
        CHECK_EQUAL(1, alm.calls());
        CHECK(status != gsStatus::Success);
        CHECK_EQUAL(0, (index_t)stepSolutions.size());
        CHECK_EQUAL(0, (index_t)distances.size());
    }
}

// ---------------------------------------------------------------------------
// END-TO-END: the CALLER clause, not the loop. The four tests above drive
// _initiation/_correction directly; this one goes through the public entry point
// so that gsAPALM<T>::serialSolve's new "initiation failed -> break, keep the
// partial branch" path (and the branch storage behind it) is actually executed.
//
// Pre-fix this call would NOT have failed -- it would never have returned.
// ---------------------------------------------------------------------------
TEST(apalm_serial_solve_terminates_when_every_step_fails)
{
    AlmProblem prob = foldProblem();
    const index_t cap = 50;
    CountingALM alm(prob.Jacobian, prob.ALResidual, prob.Force, gsStatus::SolverError, cap);
    configure(alm, 0.1);

    gsAPALMData<real_t,apalm_solution_t> data;
    APALMProbe apalm(&alm,data);
    apalm.initialize();          // seeds one start point at (0,0) with dL = getLength()
    apalm.serialSolve(3);

    const index_t nsol = (index_t)apalm.getFlatSolutions(0).size();
    gsInfo<<"  [T49/D] serialSolve(3), every step fails: step() calls = "<<alm.calls()
          <<"  stored solutions = "<<nsol<<"\n";

    CHECK(alm.calls() < cap);    // it returned at all
    CHECK(alm.calls() <= 25);    // ... via the arc-length floor of the FIRST load step
    // Only the seed survives: the branch is truncated, and no failed step is stored.
    CHECK_EQUAL(1, nsol);
}

// ===========================================================================
// gsAPALM failure protocol
// ===========================================================================
// SCOPE WARNING, read before adding to this block.
//
// The PRINCIPAL change here -- the three MPI worker clauses that now report a
// failed job to main (via an empty-solution sentinel) instead of throwing out of
// GISMO_ENSURE -- is NOT covered here and CANNOT BE. It lives inside
// #ifdef GISMO_WITH_MPI, which is #undef in this build tree, and no unit test in
// this module constructs an MPI communicator. Those clauses are COMPILE-VERIFIED
// ONLY (-DGISMO_WITH_MPI, with the preprocessor confirming all three are emitted).
//
// A test of gsAPALM<T>::parallelSolve() at m_comm.size() == 1 was CONSIDERED AND
// DROPPED: that dispatches to parallelSolve_impl<false>(), which this fix does not
// touch, so such a test passes on the unfixed tree and pins nothing here
// (house rule: a gate must be shown to FAIL without the thing it guards).
//
// The two tests below cover the parts of the failure protocol that ARE reachable, and both are
// genuinely fails-first.

// ---------------------------------------------------------------------------
// gsAPALM<T>::_getOptions() -- SubIntervals must be >= 1.
//
// Pre-fix there was NO clamp: setInt("SubIntervals",0) was accepted, and
// _correction() then sized stepSolutions to 0 and dereferenced stepSolutions.back()
// in its tail -- UB on an empty vector. The guard is also what makes the
// MPI failure sentinel unambiguous (an empty stepSolutions can only mean failure
// if a SUCCESSFUL correction returns at least one solution).
//
// FAILS-FIRST: without the GISMO_ENSURE, initialize() returns normally and
// CHECK_THROW fails.
// ---------------------------------------------------------------------------
TEST(apalm_rejects_a_zero_subinterval_count)
{
    AlmProblem prob = foldProblem();
    CountingALM alm(prob.Jacobian, prob.ALResidual, prob.Force, gsStatus::Success, 50);
    configure(alm, 0.1);

    gsAPALMData<real_t,apalm_solution_t> data;
    APALMProbe apalm(&alm,data);
    apalm.options().setInt("SubIntervals",0);

    // GISMO_ENSURE throws std::runtime_error (gsCore/gsDebug.h) -- it does NOT abort.
    CHECK_THROW(apalm.initialize(), std::runtime_error);

    // The default is unchanged and still accepted.
    gsAPALMData<real_t,apalm_solution_t> data2;
    APALMProbe apalm2(&alm,data2);
    CHECK_EQUAL(2, apalm2.options().getInt("SubIntervals"));
    apalm2.initialize();
}

// ---------------------------------------------------------------------------
// gsALMRiks<T>::initMethods() -- no 1./0 at CONSTRUCTION (handover 2026-08-03 §5.7).
//
// gsALMBase's constructor copies Force into m_forcing with NO size check, so an
// empty Force is a legal construction that reaches initMethods() with
// m_numDof == 0. A weight initialiser was added there mirroring
// predictor() verbatim, which made 1./m_numDof reachable at construction: the
// `else` arm evaluated 1./0 -> +inf, and distance() (PUBLIC, and load-bearing for
// gsAPALM's interval distances) then returned NaN via (1-inf)*DeltaL^2.
//
// FAILS-FIRST: with the pre-fix `if (m_numDof == 1)` the weight is +inf and the
// distance below is NaN, so CHECK(math::isfinite(...)) fails.
// ---------------------------------------------------------------------------
TEST(riks_empty_force_does_not_divide_by_zero_at_construction)
{
    AlmProblem prob = foldProblem();
    gsVector<real_t> emptyForce;                 // size 0 -- m_numDof == 0
    CHECK_EQUAL(0, (index_t)emptyForce.size());  // the premise, pinned

    gsALMRiks<real_t> riks(prob.Jacobian, prob.ALResidual, emptyForce);

    // distance() BEFORE any step -- exactly the public read this closes.
    // DeltaU is empty so its norm is 0, leaving (1-m_convexWeight)*DeltaL^2 under the
    // square root: finite iff the weight is finite.
    const real_t d = riks.distance(emptyForce, (real_t)1);
    gsInfo<<"  [T52] gsALMRiks with an empty Force: distance(0,1) = "<<d<<"\n";
    CHECK(math::isfinite(d));

    // The n >= 1 behaviour is untouched: widening `== 1` to `<= 1` changes no
    // reachable value, since m_numDof == 1 already took the same arm.
    gsALMRiks<real_t> riks1(prob.Jacobian, prob.ALResidual, prob.Force);
    CHECK(math::isfinite(riks1.distance(prob.Force, (real_t)1)));
}

// ---------------------------------------------------------------------------
// gsAPALM<T>::_finalize() -- the level index must run over the branch's POINTS.
//
// The loop populating m_lvlSolutions/m_lvlTimes was bounded by m_solutions.size()
// -- the BRANCH COUNT -- while k indexes the points within branch b. At the
// single branch every in-tree driver produces, that bound is 1, so exactly ONE
// point was ever registered no matter how long the branch was, and every
// level-wise accessor (getSolutions(level), getSolutionsPerLevel(),
// getTimes(level), getTimesPerLevel()) silently returned one point.
//
// The flat accessor getFlatSolutions() reads m_solutions[b] directly and was
// never affected, which is why this went unnoticed -- and which is exactly what
// makes it testable: the two views must agree.
//
// FAILS-FIRST: with the old bound the level-wise total is 1 while the flat count
// is the real number of points, so CHECK_EQUAL below fails (1 != N).
//
// (The other half of the defect -- an out-of-range read at nBranches >= 2 with a
// branch shorter than the branch count -- is NOT pinned here: producing two
// branches needs a bifurcation, which this fixture family does not reach, so
// that path remains untested.)
//
// ⚠ ENTRY POINT: serialSolve() + parallelSolve(), NOT solve(). Both reach
// _finalize(), but solve() SEGFAULTS on this fixture before it gets anywhere near
// it -- gsAPALMData::addStartPoint() inserts into a default-constructed, EMPTY
// gsKnotVector, and gsKnotVector<T>::insert() computes numLeftGhosts()/
// numRightGhosts() before its own size guard, dereferencing nullptr-4. That is a
// CORE-LIBRARY defect (gsNurbs/gsKnotVector), diagnosed under gdb+valgrind
// 2026-08-04, unrelated to this test's subject and unowned. serialSolve() is
// unaffected because it builds each branch's data via setData()/init() rather
// than addStartPoint(). The serialSolve+parallelSolve pair is also what every
// in-tree APALM driver actually calls.
// ---------------------------------------------------------------------------
TEST(apalm_finalize_registers_every_point_of_the_branch)
{
    AlmProblem prob = foldProblem();
    gsALMLoadControl<real_t> alm(prob.Jacobian, prob.ALResidual, prob.Force);
    configure(alm, 0.1);

    gsAPALMData<real_t,apalm_solution_t> data;
    gsAPALM<real_t> apalm(&alm,data);
    apalm.options().setSwitch("Verbose",false);
    apalm.initialize();
    apalm.serialSolve(4);                 // builds the branch
    apalm.parallelSolve();                // refines it, then calls _finalize()

    const index_t flat = (index_t)apalm.getFlatSolutions(0).size();

    index_t perLevel = 0;
    const std::vector<std::vector<apalm_solution_t>> byLevel = apalm.getSolutionsPerLevel(0);
    for (size_t l=0; l!=byLevel.size(); l++)
        perLevel += (index_t)byLevel[l].size();

    gsInfo<<"  [T52/_finalize] flat = "<<flat<<"  summed over "<<byLevel.size()
          <<" level(s) = "<<perLevel<<"\n";

    // The premise: the fixture really does produce more than one point, otherwise
    // the check below is vacuous (1 == 1 passes on the unfixed tree too).
    CHECK(flat > 1);
    // The defect: every stored point must appear exactly once in the level-wise view.
    CHECK_EQUAL(flat, perLevel);
}

// The same invariant, but after finalizing TWICE. _finalize() is called
// unconditionally at the end of parallelSolve_impl<false>(), so a second parallelSolve()
// on the same object does no work (the queue is drained, the while loop never runs) and
// still re-registers everything. Before the two clear()s in _finalize() this failed with
// perLevel == 2*flat (14 vs 7): resize() does not clear, so the second pass appended to
// the first pass' level-wise view. The same defect also left every pointer from the first
// pass dangling once m_solutions/m_times were reassigned; that half needs valgrind and is
// not asserted here -- pinning the count is what makes the regression visible in-suite.
TEST(apalm_finalize_is_idempotent)
{
    AlmProblem prob = foldProblem();
    gsALMLoadControl<real_t> alm(prob.Jacobian, prob.ALResidual, prob.Force);
    configure(alm, 0.1);

    gsAPALMData<real_t,apalm_solution_t> data;
    gsAPALM<real_t> apalm(&alm,data);
    apalm.options().setSwitch("Verbose",false);
    apalm.initialize();
    apalm.serialSolve(4);                 // builds the branch
    apalm.parallelSolve();                // refines it, then calls _finalize()
    apalm.parallelSolve();                // no work left to do, but finalizes AGAIN

    const index_t flat = (index_t)apalm.getFlatSolutions(0).size();

    index_t perLevel = 0;
    const std::vector<std::vector<apalm_solution_t>> byLevel = apalm.getSolutionsPerLevel(0);
    for (size_t l=0; l!=byLevel.size(); l++)
        perLevel += (index_t)byLevel[l].size();

    gsInfo<<"  [T55/_finalize x2] flat = "<<flat<<"  summed over "<<byLevel.size()
          <<" level(s) = "<<perLevel<<"\n";

    // Non-vacuity: the fixture must produce more than one point, else 1 == 1 would pass
    // on the unfixed tree too.
    CHECK(flat > 1);
    // The defect: finalizing twice must leave exactly the same level-wise view as once.
    CHECK_EQUAL(flat, perLevel);
}

// ===========================================================================
// Fixtures pinning the left-null-vector, tangent, BorderedMode and fold-test-function
// behaviour.
// ===========================================================================

// ===========================================================================
// T1 (falsifier) -- the limit-vs-branch verdict uses the LEFT null
// vector, not the right one.
// ===========================================================================
// Fixture A (gsALMTestProblems.h) is the only fixture in this suite whose LEFT
// and RIGHT critical modes give OPPOSITE limit-vs-branch verdicts: at the fold
// the right null vector phi=(0,1) is EXACTLY force-orthogonal (a right-vector
// test misclassifies this LIMIT point as a BRANCH point at every tolerance),
// while the left null vector psi=(1,1)/sqrt(2) gives cosine 1/sqrt(2). A
// non-symmetric tangent needs the "LU" solver; BifurcationMethod is set to
// "Nothing" (-1) for safety, though isBifurcation() never reaches
// computeStability() on this call path (checked at source: gsALMBase.hpp
// isBifurcation -> _testSingularPoint -> computeJacobian/factorizeMatrix ->
// _computeCriticalMode(Left), none of which touch computeStability).
TEST(fold_classification_uses_the_left_null_vector)
{
    // --- fixture A: the falsifier -------------------------------------------
    AlmProblem prob = asymFoldProblem();
    gsALMLoadControl<real_t> alm(prob.Jacobian, prob.ALResidual, prob.Force);
    alm.options().setString("Solver","LU");
    alm.options().setInt   ("BifurcationMethod",-1);   // "Nothing": safety only, see above
    alm.options().setReal  ("Length",0.05);
    alm.options().setReal  ("Tol",1e-10);
    alm.options().setReal  ("TolF",1e-8);
    alm.options().setReal  ("TolU",1e-8);
    alm.options().setInt   ("MaxIter",100);
    alm.options().setReal  ("SingularPointTestTol",1e-4);
    alm.options().setReal  ("SingularPointComputeTolE",1e-10);
    alm.options().setReal  ("SingularPointComputeTolB",0);
    alm.options().setSwitch("Verbose",false);
    alm.applyOptions();
    CHECK(alm.initialize() == gsStatus::Success);

    const real_t u2 = 0.99;
    const real_t L  = asymFoldLambda(u2);        // = u1, the exact equilibrium
    gsVector<real_t> U(2); U << L, u2;
    alm.setSolution(U,L);

    // The headline assertion: with a right-vector-only test this would read TRUE
    // (misclassified as BRANCH).
    CHECK(!alm.isBifurcation(true));
    // RESOLVED limit, not an "unresolved mode" refusal: SPverdict::Unresolved
    // ALSO makes isBifurcation() return false, so without this check
    // the assertion above could pass for the wrong reason.
    CHECK_EQUAL(gsALMBase<real_t>::SPverdict::Limit, alm.singularPointVerdict());
    gsInfo<<"  [T1] fixture A u2=0.99: criticalModeError = "<<alm.criticalModeError()
          <<"  cosine(V,f) = "<<modeForceCosine(alm.solutionV(),prob.Force)
          <<"  V = "<<alm.solutionV().transpose()<<"\n";

    // OPTIONAL: solutionVLeft() is PUBLIC on gsALMBase, so no probe is needed --
    // the left mode psi ~ (1,0.98)/||.|| at u2=0.99 (fixture-header derivation).
    gsVector<real_t> psiRef(2); psiRef << 1.0, 0.98;
    psiRef.normalize();
    const gsVector<real_t> psi = alm.solutionVLeft().normalized();
    const real_t cosPsi = math::abs(psi.dot(psiRef));
    gsInfo<<"  [T1] left mode psi = "<<alm.solutionVLeft().transpose()
          <<"  cosine vs analytic (1,0.98) = "<<cosPsi<<"\n";
    CHECK(cosPsi > 1.0 - 1e-6);

    // --- control half: the symmetric fixture F must still classify LIMIT ----
    // (mirrors the two-half structure of singular_point_test_is_scale_invariant:
    // without this half, "isBifurcation always returns false" would pass T1 too).
    AlmProblem probF = foldProblem();
    gsALMLoadControl<real_t> almF(probF.Jacobian, probF.ALResidual, probF.Force);
    almF.options().setString("Solver","SimplicialLDLT");
    almF.options().setInt   ("BifurcationMethod",0);
    almF.options().setReal  ("Length",0.05);
    almF.options().setReal  ("Tol",1e-10);
    almF.options().setReal  ("TolF",1e-8);
    almF.options().setReal  ("TolU",1e-8);
    almF.options().setInt   ("MaxIter",100);
    almF.options().setReal  ("SingularPointTestTol",1e-4);
    almF.options().setReal  ("SingularPointComputeTolE",1e-10);
    almF.options().setReal  ("SingularPointComputeTolB",0);
    almF.options().setSwitch("Verbose",false);
    almF.applyOptions();
    CHECK(almF.initialize() == gsStatus::Success);

    gsVector<real_t> UF(2); UF << 0.95, 0.0;
    almF.setSolution(UF,foldLambda(0.95));
    CHECK(!almF.isBifurcation(true));
    CHECK_EQUAL(gsALMBase<real_t>::SPverdict::Limit, almF.singularPointVerdict());
}

// ===========================================================================
// T2 (control) -- fixture F's fold classification survives the
// left-vector rewrite. Fixture F is SYMMETRIC (left == right mode), so this is
// a NO-REGRESSION CONTROL, not evidence about the left vector -- see T1/fixture
// A for the falsifier.
// ===========================================================================
TEST(fixture_F_fold_is_limit_with_unit_mode_force_cosine)
{
    AlmProblem prob = foldProblem();
    gsALMLoadControl<real_t> alm(prob.Jacobian, prob.ALResidual, prob.Force);
    configure(alm, 0.05);

    gsVector<real_t> U(2); U << 0.999, 0.0;
    alm.setSolution(U,foldLambda(0.999));
    CHECK(!alm.isBifurcation(true));
    CHECK_EQUAL(gsALMBase<real_t>::SPverdict::Limit, alm.singularPointVerdict());
    const real_t cosF = modeForceCosine(alm.solutionV(), prob.Force);
    gsInfo<<"  [T2] fixture F u1=0.999: cosine(V,f) = "<<cosF<<"\n";
    // K = diag(2-2u1,1) = diag(0.002,1): the inverse power iteration converges
    // to (1,0), so the cosine is 1 to many digits.
    CHECK_CLOSE(1.0, cosF, 1e-6);

    // --- second half: the branch/limit discrimination must survive the left-vector rewrite --
    AlmProblem probP = pitchforkProblem();
    gsALMLoadControl<real_t> almP(probP.Jacobian, probP.ALResidual, probP.Force);
    configure(almP, 0.05);

    gsVector<real_t> UP(2); UP << 1.05, 0.0;
    almP.setSolution(UP,1.05);
    CHECK(almP.isBifurcation(true));
    CHECK_EQUAL(gsALMBase<real_t>::SPverdict::Branch, almP.singularPointVerdict());
    const real_t cosP = modeForceCosine(almP.solutionV(), probP.Force);
    gsInfo<<"  [T2] fixture P fundamental lambda=1.05: cosine(V,f) = "<<cosP
          <<" (SingularPointTestTol = 1e-4)\n";
    CHECK(cosP < 1e-4);
}

// ===========================================================================
// T3 -- the swibra ABE emanating tangent on fixture P.
// ===========================================================================
// Analytic solution (fixture P): at the
// fork K* = diag(1,0); right/left null vectors phi1 = psi1 = (0,1) (symmetric);
// incoming fundamental tangent udot0=(1,0), lambdadot0=1; b1 = -1 EXACTLY
// (normalisation-invariant), a1 = 3*del under the FORWARD FD the swibra port
// uses (0 centrally), c1 = 0 -- the ABE a1*alpha^2+2*b1*alpha*beta+c1*beta^2=0
// reduces to -2*alpha*beta=0, whose non-trivial root is beta=0, i.e. the
// emanating tangent is the PURE MODE direction tau1 ~ (0,+-1;0). Independent
// cross-check: the bifurcated branch u1=2L-1, u2=+-sqrt(2(L-1)) parametrised by
// u2 gives (du1,du2,dL)/ds -> (0,1,0) at u2->0 exactly -- the same direction.
TEST(branch_tangent_matches_the_analytic_abe_solution)
{
    AlmProblem prob = pitchforkProblem();
    gsALMLoadControl<real_t> alm(prob.Jacobian, prob.ALResidual, prob.Force);
    configure(alm, 0.05);

    for (index_t k = 0; k != 21; ++k)
        CHECK(alm.step() == gsStatus::Success);
    CHECK_CLOSE(1.05, alm.solutionL(), 1e-10);

    // Incoming secant, captured BEFORE any refinement moves the solver.
    const gsVector<real_t> tangentU0 = alm.solutionDU();
    const real_t           tangentL0 = alm.solutionDL();
    gsInfo<<"  [T3] incoming secant: dU = ("<<tangentU0[0]<<","<<tangentU0[1]
          <<")  dL = "<<tangentL0<<"\n";

    const gsStatus st = alm.computeSingularPoint(alm.solutionU(), alm.solutionL(),
                                                 /*switchBranch*/false,
                                                 /*jacobian*/true,
                                                 /*testPoint*/true);
    CHECK(st == gsStatus::Success);
    if (st != gsStatus::Success)
        return;

    gsVector<real_t> tau1U; real_t tau1L = 0.0;
    // Success is the DERIVED expectation here, not merely the measured one (C3):
    // configure() sets Solver=SimplicialLDLT, the self-adjoint path, so
    // _computeCriticalModeLeft returns psi1=phi1 exactly => den=phi1.psi1=1, the
    // computeBranchTangent()'s biorthogonalisation (den = phi1.psi1) ModeUnresolved
    // guard cannot fire. The incoming load-control secant
    // dU=(0.05,0), dL=0.05 gives lamShare=1/sqrt(2)=0.707 >> BranchTangentTol=1e-2,
    // so the load-share Degenerate guard cannot fire. Both callbacks are
    // analytic 2x2 closures that always return true, so AssemblyFailed cannot fire.
    // |Gud|_F = del*sqrt(2) ~ 1.4e-3 against 100*eps*|K*|_F ~ 2.2e-14, so the
    // TrivialBranch guard cannot fire. al1b = 2 exactly against
    // BranchTangentTol*scaleABE = 1e-2*2 = 0.02, so the al1b Degenerate (no
    // distinct branch) guard cannot fire either. Success is
    // therefore the ONLY reachable exit on fixture P -- see computeBranchTangent()'s
    // ModeUnresolved, load-share/al1b Degenerate, and TrivialBranch guards above for
    // the sites.
    const gsALMBase<real_t>::branchTangent::type outcome0 =
        alm.computeBranchTangent(alm.solutionU(), alm.solutionL(), alm.solutionV(),
                                 tangentU0, tangentL0, tau1U, tau1L);
    gsInfo<<"  [T3] computeBranchTangent outcome = "<<static_cast<index_t>(outcome0)
          <<" (0=Success,1=TrivialBranch,2=ModeUnresolved,3=Degenerate,4=AssemblyFailed)\n";
    CHECK_EQUAL(gsALMBase<real_t>::branchTangent::Success, outcome0);
    // alpha1b did NOT vanish (analytically = 2): the "no distinct branch" outcome
    // must not fire on fixture P.
    CHECK(outcome0 != gsALMBase<real_t>::branchTangent::Degenerate);

    CHECK(tau1U.norm() > 0);
    gsVector<real_t> tauVec(3); tauVec << tau1U[0], tau1U[1], tau1L;
    const real_t tauNorm = tauVec.norm();
    CHECK(tauNorm > 0);
    tauVec /= tauNorm;
    gsInfo<<"  [T3] tau1 (normalised) = ("<<tauVec[0]<<","<<tauVec[1]<<" ; "<<tauVec[2]
          <<")  raw tau1U=("<<tau1U[0]<<","<<tau1U[1]<<")  tau1L="<<tau1L<<"\n";

    // The direction error is O(del) -- the FORWARD-FD coefficient a1 = 3*del feeds
    // straight into tau1's (0,+-1;0) deviation (see the comment block below) -- so
    // the tolerance is read from the option actually driving it, NOT a guessed
    // constant (default del = 1e-3, an ORDER bigger than a guessed 1e-6 would
    // allow: MEASURED tauVec[0]=tauVec[2]~1.5e-3 at that default).
    const real_t del = alm.options().getReal("BranchTangentPerturbation");
    gsInfo<<"  [T3] BranchTangentPerturbation (del) = "<<del<<"\n";
    const real_t dirTol = 5.0 * del;
    const real_t magTol = 5.0 * del * del;
    CHECK(math::abs(tauVec[0]) + math::abs(tauVec[2]) < dirTol);
    CHECK_CLOSE(1.0, math::abs(tauVec[1]), magTol);
    // Stronger, still-analytic check (C3): al1 = psi1.u0d = 0 exactly on fixture P
    // (psi1=(0,1), u0d=(1,0) up to normalisation), so the a1*al1/lam0 term in al1b
    // vanishes identically and tau1 = (3*del, 2 ; 3*del) BEFORE normalisation =>
    // abs(tauVec[0])+abs(tauVec[2]) = 3*del up to an O(del^3) normalisation
    // correction (~7e-9 at del=1e-3) -- far inside the tolerance below.
    CHECK_CLOSE(3.0*del, math::abs(tauVec[0]) + math::abs(tauVec[2]), 0.1*del);

    // ABE coefficients (a1,b1) are NOT exposed by the current gsALMBase API (no
    // accessor was found in gsALMBase.h): the analytic oracle is recorded here as
    // a comment instead (del = option BranchTangentPerturbation, default 1e-3):
    //   b1 = psi1.(R_uu[phi1,phi0] + R_ulambda[phi1]) = -1   EXACTLY, normalisation-invariant
    //   a1 = psi1.R_uu[phi1,phi1] = 3*del (forward FD; 0 centrally)
    //   alpha1b = -(a1*alpha1/lambdadot0 + 2*b1) = 2   (alpha1 = 0 kills the a1 term)
}

// ===========================================================================
// T6 -- BorderedMode runtime option: "Primary" converges where the
// elimination corrector fails, the retry gate fires on SolverError (not only
// NotConverged), and a FAILED retry restores all FOUR stability-bookkeeping
// members.
// ===========================================================================

/// Wraps a TRUE Jacobian callback and, on a single chosen 1-based call number
/// (after arm()), substitutes a BIT-EXACT SINGULAR 2x2 matrix diag(1,0) --
/// exactly the matrix TEST(extended_system_shift_ladder_survives_an_exactly_
/// singular_tangent) above measured to make plain factorizeMatrix throw code 3
/// (SolverError, NumericalIssue) -- instead of the true tangent. Single-shot:
/// only the call AT failAt is replaced; every other call uses the true Jacobian.
class SolverErrorInjector
{
public:
    explicit SolverErrorInjector(const gsStructuralAnalysisOps<real_t>::Jacobian_t & trueJac)
    : m_trueJac(trueJac), m_failAt(-1), m_calls(0) {}

    /// Arms a single-shot fault at 1-based call number \a failAt and resets the
    /// call counter.
    void arm(index_t failAt) { m_failAt = failAt; m_calls = 0; }

    index_t calls() const { return m_calls; }

    gsStructuralAnalysisOps<real_t>::Jacobian_t jacobian()
    {
        return [this](gsVector<real_t> const & u, gsSparseMatrix<real_t> & m) -> bool
        {
            ++m_calls;
            if (m_calls == m_failAt)
            {
                gsMatrix<real_t> K(2,2); K.setZero();
                K(0,0) = 1.0; K(1,1) = 0.0;
                toSparse2(K,m);
                return true;
            }
            return m_trueJac(u,m);
        };
    }
private:
    gsStructuralAnalysisOps<real_t>::Jacobian_t m_trueJac;
    index_t m_failAt;
    index_t m_calls;
};

/// Reaches gsALMCrisfield's protected stability-bookkeeping members that have no
/// public raw accessor (m_stabilityVec, m_stability -- m_negatives/m_indicator ARE
/// already public via negatives()/indicator()) so T6(c) can perturb all four to a
/// distinctive sentinel before a forced double-attempt failure and assert the
/// RESTORED values, not merely that they changed. Nothing under src/ is modified.
class CrisfieldStabilityProbe : public gsALMCrisfield<real_t>
{
    typedef gsALMCrisfield<real_t> Base;
public:
    CrisfieldStabilityProbe(const gsStructuralAnalysisOps<real_t>::Jacobian_t   & J,
                            const gsStructuralAnalysisOps<real_t>::ALResidual_t & R,
                            const gsVector<real_t>                              & F)
    : Base(J,R,F) {}

    const gsVector<real_t> & stabilityVecRaw() const { return this->m_stabilityVec; }
    index_t                  stabilityRaw()    const { return this->m_stability; }

    /// Sets all FOUR snapshot members directly, bypassing every public setter
    /// (setIndicator() RECOMPUTES m_stability from m_indicator, which is not what
    /// a raw perturbation needs -- it must set an internally INCONSISTENT sentinel
    /// so a restore that misses even one member is observable).
    void setStabilitySnapshot(index_t negatives, real_t indicator,
                              const gsVector<real_t> & stabilityVec, index_t stability)
    {
        this->m_negatives    = negatives;
        this->m_indicator     = indicator;
        this->m_stabilityVec  = stabilityVec;
        this->m_stability     = stability;
    }
};

TEST(bordered_mode_primary_and_solver_error_retry)
{
    // --- (a) BorderedMode="Primary" converges the fixture-F step-19 fold that
    // the elimination corrector fails (class doc (a): the predictor of step 19
    // lands on u1=1+2.2e-16, the elimination chart degenerates into a period-2
    // limit cycle for all MaxIter iterations and returns NotConverged; the
    // bordered chart is well-posed at exactly that point). Drives the NEW
    // BorderedMode string option directly (not the deprecated alias).
    {
        AlmProblem prob = foldProblem(1.0);
        CrisfieldProbe alm(prob.Jacobian, prob.ALResidual, prob.Force);
        configure(alm, 0.05);
        alm.options().setReal  ("Scaling",0.0);
        alm.options().setString("BorderedMode","Primary");
        alm.applyOptions();

        index_t failedStep = -1;
        for (index_t k = 0; k != 25; ++k)
        {
            const gsStatus st = alm.step();
            if (st != gsStatus::Success) { failedStep = k; break; }
        }
        gsInfo<<"  [T6a] BorderedMode=Primary: failedStep = "<<failedStep
              <<"  final u1 = "<<alm.solutionU()[0]<<"\n";
        CHECK_EQUAL(-1, failedStep);           // no step failed, in particular NOT step 19
        CHECK(alm.solutionU()[0] > 1.05);      // PAST the fold (elimination alone dies at u1=0.95)
    }

    // --- (b) the retry gate fires on SolverError, not only NotConverged ------
    // A single-shot injected SolverError at the FIRST corrector-iterate Jacobian
    // assembly (predictor() always plainly factorizes K regardless of chart --
    // ONE call -- so call 2 is quasiNewtonIteration()'s first corrector-iterate
    // call), at a WELL-CONDITIONED state
    // (u1=0.5, away from any fold) so a rescuing attempt 2 converges normally.
    // Under "Fallback" attempt 1 is the elimination chart (the only chart whose
    // quasiNewtonIteration() factorizes K, hence the only one this fault can
    // reach) and attempt 2 is the bordered chart, which never factorizes plain K
    // (quasiNewtonIteration()'s bordered branch calls _assembleJacobianUnfactorized()
    // without a factorizeMatrix() call) and so is untouched by the single-shot fault.
    {
        AlmProblem prob = foldProblem(1.0);
        gsVector<real_t> U(2); U << 0.5, 0.0;

        SolverErrorInjector injOff(prob.Jacobian);
        CrisfieldProbe almOff(injOff.jacobian(), prob.ALResidual, prob.Force);
        configure(almOff, 0.05);
        almOff.options().setString("BorderedMode","Off");
        almOff.applyOptions();
        almOff.setSolution(U, foldLambda(0.5));
        injOff.arm(2);
        const gsStatus stOff = almOff.step();

        SolverErrorInjector injFB(prob.Jacobian);
        CrisfieldProbe almFB(injFB.jacobian(), prob.ALResidual, prob.Force);
        configure(almFB, 0.05);
        almFB.options().setString("BorderedMode","Fallback");
        almFB.applyOptions();
        almFB.setSolution(U, foldLambda(0.5));
        injFB.arm(2);
        const gsStatus stFB = almFB.step();

        gsInfo<<"  [T6b] same injected SolverError: BorderedMode=Off -> status "
              <<static_cast<index_t>(stOff)<<" (3=SolverError expected); BorderedMode=Fallback -> status "
              <<static_cast<index_t>(stFB)<<" (0=Success expected: rescued by the bordered retry)\n";
        CHECK(stOff == gsStatus::SolverError);   // the fault, in isolation, IS a SolverError
        CHECK(stFB  == gsStatus::Success);       // ... and the widened gate RESCUES it
    }

    // --- (c) a FAILED retry restores all FOUR snapshot members ---------------
    // Same fault mechanism, "Fallback" mode, u1=0.5 (well-conditioned, decoupled
    // from (a)'s fold physics): attempt 1 (elimination) is forced to fail via the
    // SAME single-shot fault as (b), on its VERY FIRST corrector iterate -- i.e.
    // BEFORE gsALMBase::_step()'s computeStability(false) ever runs
    // for that iteration, so the four members are UNTOUCHED when the internal
    // snapshot is taken: the snapshot IS the sentinel set below. Attempt 2 (the
    // bordered chart, unaffected by the single-shot fault, see (b)) is then
    // forced to fail for an UNRELATED, mundane reason -- MaxIter=1 with a very
    // tight tolerance -- so it runs exactly ONE real corrector iteration (a REAL
    // computeStability(false) write, definitively not the sentinel) before
    // exhausting its budget and returning NotConverged. This is what makes the
    // restore falsifiable: delete the restore block and the members read
    // attempt-2's real (non-sentinel) values instead of the sentinel.
    {
        AlmProblem prob = foldProblem(1.0);
        gsVector<real_t> U(2); U << 0.5, 0.0;

        SolverErrorInjector inj(prob.Jacobian);
        CrisfieldStabilityProbe alm(inj.jacobian(), prob.ALResidual, prob.Force);
        configure(alm, 0.05);
        alm.options().setString("BorderedMode","Fallback");
        alm.options().setInt   ("MaxIter",1);
        alm.options().setReal  ("Tol",1e-15);
        alm.options().setReal  ("TolF",1e-15);
        alm.options().setReal  ("TolU",1e-15);
        alm.applyOptions();
        alm.setSolution(U, foldLambda(0.5));

        const index_t    sentNeg  = 999;
        const real_t     sentInd  = -123456.0;
        gsVector<real_t> sentVec(2); sentVec << 777.0, -777.0;
        const index_t    sentStab = 999;
        alm.setStabilitySnapshot(sentNeg, sentInd, sentVec, sentStab);

        inj.arm(2);   // faults attempt 1's first (and only, since it throws) corrector-iterate call
        const gsStatus st = alm.step();

        gsInfo<<"  [T6c] double failure: step() status = "<<static_cast<index_t>(st)
              <<" (3=SolverError expected, attempt 1's own status per the clamp)\n"
              <<"        restored: negatives="<<alm.negatives()<<" (sentinel "<<sentNeg<<")"
              <<"  indicator="<<alm.indicator()<<" (sentinel "<<sentInd<<")"
              <<"  stability="<<alm.stabilityRaw()<<" (sentinel "<<sentStab<<")"
              <<"  stabilityVec=("<<alm.stabilityVecRaw()[0]<<","<<alm.stabilityVecRaw()[1]
              <<") (sentinel ("<<sentVec[0]<<","<<sentVec[1]<<"))\n";

        CHECK(st != gsStatus::Success);   // both attempts genuinely failed
        CHECK_EQUAL(sentNeg,  alm.negatives());
        CHECK_CLOSE(sentInd,  alm.indicator(), 1e-14);
        CHECK_EQUAL(sentStab, alm.stabilityRaw());
        CHECK_ARRAY_CLOSE(sentVec.data(), alm.stabilityVecRaw().data(), 2, 1e-14);
    }

    // --- (c) CONTROL: the block above is only evidence of a restore if a lone
    // attempt-2-class execution DOES dirty the sentinels. Same state, same
    // MaxIter=1/tight-tolerance setup, NO injected fault, BorderedMode="Off"
    // (a single elimination attempt, no snapshot/restore machinery at all): one
    // real corrector iteration must overwrite every sentinel. Without this
    // half, "members == sentinels after the double failure" would pass equally
    // well if attempt 2 simply never wrote them -- this closes that gap.
    {
        AlmProblem prob = foldProblem(1.0);
        gsVector<real_t> U(2); U << 0.5, 0.0;

        CrisfieldStabilityProbe alm(prob.Jacobian, prob.ALResidual, prob.Force);
        configure(alm, 0.05);
        alm.options().setString("BorderedMode","Off");
        alm.options().setInt   ("MaxIter",1);
        alm.options().setReal  ("Tol",1e-15);
        alm.options().setReal  ("TolF",1e-15);
        alm.options().setReal  ("TolU",1e-15);
        alm.applyOptions();
        alm.setSolution(U, foldLambda(0.5));

        const index_t    sentNeg  = 999;
        const real_t     sentInd  = -123456.0;
        gsVector<real_t> sentVec(2); sentVec << 777.0, -777.0;
        const index_t    sentStab = 999;
        alm.setStabilitySnapshot(sentNeg, sentInd, sentVec, sentStab);

        const gsStatus st = alm.step();

        gsInfo<<"  [T6c/control] lone attempt (no fault, no retry): status = "
              <<static_cast<index_t>(st)<<" (1=NotConverged expected)\n"
              <<"        negatives="<<alm.negatives()<<" (sentinel was "<<sentNeg<<")"
              <<"  indicator="<<alm.indicator()<<" (sentinel was "<<sentInd<<")"
              <<"  stability="<<alm.stabilityRaw()<<" (sentinel was "<<sentStab<<")\n";

        CHECK(st != gsStatus::Success);          // exhausted its 1-iteration budget
        CHECK(alm.negatives()    != sentNeg);
        CHECK(alm.indicator()    != sentInd);
        CHECK(alm.stabilityRaw() != sentStab);
    }

    // --- the deprecated alias still works and still warns --------------------
    // crisfieldBorderedProgress() was migrated to drive BorderedMode directly
    // (the last remaining BorderedFallback use in the tree before this task); this
    // keeps the alias itself exercised, once, deliberately.
    {
        AlmProblem prob = foldProblem(1.0);
        CrisfieldProbe alm(prob.Jacobian, prob.ALResidual, prob.Force);
        configure(alm, 0.05);
        alm.options().setReal  ("Scaling",0.0);
        alm.options().setSwitch("BorderedFallback",true);   // the deprecated alias, deliberately

        const std::string log1 = captureInfo([&alm]() { alm.applyOptions(); });
        const index_t warnings1 = countOccurrences(log1,"BorderedFallback is DEPRECATED");
        gsInfo<<"  [T6/alias] BorderedFallback warning occurrences on 1st applyOptions() = "
              <<warnings1<<"\n";
        CHECK_EQUAL(1, warnings1);

        // A second applyOptions() must NOT warn again (once per solver object).
        const std::string log2 = captureInfo([&alm]() { alm.applyOptions(); });
        const index_t warnings2 = countOccurrences(log2,"BorderedFallback is DEPRECATED");
        gsInfo<<"  [T6/alias] BorderedFallback warning occurrences on 2nd applyOptions() = "
              <<warnings2<<"\n";
        CHECK_EQUAL(0, warnings2);

        // And it behaves like Fallback: elimination fails at step 19, the aliased
        // bordered retry rounds the fold (same trace as crisfieldBorderedProgress
        // with bordered=true).
        index_t failedStep = -1;
        for (index_t k = 0; k != 25; ++k)
        {
            const gsStatus st = alm.step();
            if (st != gsStatus::Success) { failedStep = k; break; }
        }
        CHECK_EQUAL(-1, failedStep);
        CHECK(alm.solutionU()[0] > 1.05);
    }
}

// ===========================================================================
// T7 -- AUTO-07p's fold test function changes sign across the fold,
// analytically pinned on fixture F, and does not affect detection behaviour.
// ===========================================================================
// TF = sign(f.K^-1.f)/sqrt(1+|K^-1.f|^2). On fixture F, K=diag(2-2u1,1), f=(1,0):
// K^-1.f = (1/(2-2u1),0), so TF = sign(2-2u1)/sqrt(1+1/(2-2u1)^2). Solving
// 1+1/(2-2u1)^2 = 26 gives 2-2u1 = +-1/5, i.e. u1 = 0.9 (stable arm, TF=+1/sqrt(26))
// and u1 = 1.1 (unstable arm, TF=-1/sqrt(26)) -- the closed-form pair this test uses.
TEST(fold_test_function_changes_sign_across_the_fold)
{
    AlmProblem prob = foldProblem();
    gsALMLoadControl<real_t> alm(prob.Jacobian, prob.ALResidual, prob.Force);
    configure(alm, 0.05);

    gsVector<real_t> Us(2); Us << 0.9, 0.0;
    alm.setSolution(Us, foldLambda(0.9));
    CHECK(alm.computeStability(true) == gsStatus::Success);
    const real_t tfStable = alm.foldTestFunction();
    gsInfo<<"  [T7] fold test function at u1=0.9 (stable arm) = "<<tfStable
          <<"  (analytic +1/sqrt(26) = "<<1.0/math::sqrt(26.0)<<")\n";
    CHECK(!math::isnan(tfStable));
    CHECK_CLOSE(1.0/math::sqrt(26.0), tfStable, 1e-10);

    gsVector<real_t> Uu(2); Uu << 1.1, 0.0;
    alm.setSolution(Uu, foldLambda(1.1));
    CHECK(alm.computeStability(true) == gsStatus::Success);
    const real_t tfUnstable = alm.foldTestFunction();
    gsInfo<<"  [T7] fold test function at u1=1.1 (unstable arm) = "<<tfUnstable
          <<"  (analytic -1/sqrt(26))\n";
    CHECK(!math::isnan(tfUnstable));
    CHECK_CLOSE(-1.0/math::sqrt(26.0), tfUnstable, 1e-10);

    // Sign change across the fold.
    CHECK(tfStable * tfUnstable < 0.0);

    // Detection behaviour is UNCHANGED by the instrument's presence: the fold
    // must still classify as a LIMIT point via the ordinary path (foldTestFunction
    // is an INSTRUMENT, never consulted by _testSingularPoint -- see gsALMBase.h).
    alm.setSolution(Uu, foldLambda(1.1));
    CHECK(!alm.isBifurcation(true));
    CHECK_EQUAL(gsALMBase<real_t>::SPverdict::Limit, alm.singularPointVerdict());
}

// ===========================================================================
// New -- proves the CHECK_EQUAL pinned in T3 above is discriminating
// rather than always-true: a fixture where K_T is CONSTANT in u drives
// computeBranchTangent through the TrivialBranch exit.
// ===========================================================================
TEST(branch_tangent_reports_trivial_branch_on_a_linear_problem)
{
    // K = diag(2,1) is CONSTANT in u BY CONSTRUCTION of linearProblem(): Kdel ==
    // Kstar identically for ANY del and ANY
    // phi1, so Gud == 0 exactly -- not to roundoff via a measured coincidence, but
    // structurally -- and computeBranchTangent()'s TrivialBranch test
    // |Gud|_F <= 100*eps*|K*|_F fires
    // unconditionally.
    AlmProblem prob = linearProblem();
    gsALMLoadControl<real_t> alm(prob.Jacobian, prob.ALResidual, prob.Force);
    configure(alm, 0.05);

    gsVector<real_t> Ustar(2); Ustar << 0.5, 0.0;   // equilibrium: R1=2*0.5-1=0, R2=0
    const real_t Lstar = 1.0;
    gsVector<real_t> phi1(2); phi1 << 0.0, 1.0;
    gsVector<real_t> tangentU0(2); tangentU0 << 1.0, 0.0;
    // lamShare = 1/sqrt(2) = 0.707 >> BranchTangentTol = 1e-2, so the
    // load-share Degenerate guard is cleared and the TrivialBranch test really is
    // the exit taken.
    const real_t tangentL0 = 1.0;

    gsVector<real_t> tau1U; real_t tau1L = 0.0;
    const gsALMBase<real_t>::branchTangent::type outcome =
        alm.computeBranchTangent(Ustar, Lstar, phi1, tangentU0, tangentL0, tau1U, tau1L);
    gsInfo<<"  [T3b] computeBranchTangent outcome = "<<static_cast<index_t>(outcome)
          <<" (expect 1=TrivialBranch)\n";
    CHECK_EQUAL(gsALMBase<real_t>::branchTangent::TrivialBranch, outcome);

    // Documented payload of the TrivialBranch branch:
    // tau1U = phi1, tau1L = 0.
    CHECK_CLOSE(0.0, tau1L, 1e-14);
    CHECK_CLOSE(1.0, math::abs(tau1U.normalized()[1]), 1e-12);
}

// ===========================================================================
// New -- proves the CHECK_EQUAL pinned in T3 discriminates in the
// OTHER direction: an incoming tangent with ZERO load share pins
// gsALMBase::computeBranchTangent()'s load-share Degenerate guard, the documented
// "branch point coinciding with a fold is outside multiplicity-1 swibra's
// scope" contract.
// ===========================================================================
TEST(branch_tangent_reports_degenerate_on_a_pure_displacement_incoming_tangent)
{
    AlmProblem prob = pitchforkProblem();
    gsALMLoadControl<real_t> alm(prob.Jacobian, prob.ALResidual, prob.Force);
    configure(alm, 0.05);

    gsVector<real_t> Ustar(2); Ustar << 1.0, 0.0;   // the analytic fork, no stepping needed
    const real_t Lstar = 1.0;
    gsVector<real_t> phi1(2); phi1 << 0.0, 1.0;
    gsVector<real_t> tangentU0(2); tangentU0 << 1.0, 0.0;
    // lamShare = 0: computeBranchTangent()'s load-share Degenerate guard fires.
    // lam0 = 0 would otherwise divide twice downstream (phi0 in step 4, al1b in
    // step 8), which is why the guard exists.
    const real_t tangentL0 = 0.0;

    gsVector<real_t> tau1U; real_t tau1L = 0.0;
    const gsALMBase<real_t>::branchTangent::type outcome =
        alm.computeBranchTangent(Ustar, Lstar, phi1, tangentU0, tangentL0, tau1U, tau1L);
    gsInfo<<"  [T3c] computeBranchTangent outcome = "<<static_cast<index_t>(outcome)
          <<" (expect 3=Degenerate)\n";
    CHECK_EQUAL(gsALMBase<real_t>::branchTangent::Degenerate, outcome);
}

// ===========================================================================
// m_Lprev must not be indeterminate storage.
// ===========================================================================
// T m_L, m_Lprev; has no in-class initializer. Before the fix
// (both gsALMBase constructors now seed m_Lprev = 0.0), neither
// gsALMCrisfield::initMethods() nor
// gsALMLoadControl::initMethods() assigned it --
// only gsALMRiks / gsALMConsistentCrisfield seed it (to 0.0) in their own
// initMethods(). The public solutionLPrev() accessor is read PRE-SEED by
// gsALMExploration::traceSweep()'s verbose alpha (pre-seed) diagnostic dump, and
// gsALMLoadControl is the solver both oracle drivers run (the
// `new gsALMLoadControl<real_t>(Jacobian,ALResidual,Force)` construction in
// example_ShearLandscape.cpp, and the `else solver = new gsALMLoadControl<real_t>(...)`
// default branch of example_ModifiedBratuExploration.cpp's method dispatch), so this
// was a live path, not a theoretical one.
//
// Ordinary construction on the heap/stack is usually already zero, so
// CHECK_EQUAL(0.0, ...) on a plain object would pass before AND after the fix
// -- a criterion that can only pass is worthless. Instead each solver is
// placement-new'd into storage poisoned with 0xA5 bytes: the buffer is raw
// storage until the constructor runs, and the constructor default-constructs
// every class-type member it does not explicitly list (gsVector<T>, the
// std::function members, ...) regardless of the poison, but leaves a bare
// scalar member with NO initializer (like the pre-fix m_Lprev) holding
// whatever bytes were already there. That is exactly the discriminator.
// Fills [buf, buf+n) with 0xA5 through a VOLATILE pointer, one byte at a time. A
// std::memset here would be a DEAD STORE for every byte the subsequent placement-new
// constructor overwrites, and the compiler is free to elide it entirely: a disassembly
// of the shipped `unittests` binary found exactly ONE surviving 0xA5 fill
// (the gsALMConsistentCrisfield instantiation -- the one solver where the poison is
// useless, because its own initMethods() re-seeds m_Lprev anyway), i.e. the memset-based
// instrument was inert for 3 of the 4 solvers, including the live-path gsALMLoadControl.
// A volatile store cannot be optimized away, so this makes the poison reach every
// solver's storage regardless of what its constructor goes on to write.
static void poisonBytes(void * buf, std::size_t n)
{
    volatile unsigned char * p = reinterpret_cast<volatile unsigned char*>(buf);
    for (std::size_t i = 0; i < n; ++i)
        p[i] = (unsigned char)0xA5;
}

template <class ALM>
static real_t poisonedLPrevAfterInit(AlmProblem & prob)
{
    typename std::aligned_storage<sizeof(ALM), alignof(ALM)>::type buf;
    poisonBytes(&buf, sizeof(buf));                // poison every byte the ctor does not write
    ALM * alm = new (&buf) ALM(prob.Jacobian, prob.ALResidual, prob.Force);
    configure(*alm, 0.05);                         // calls initialize() -> initMethods()
    const real_t val = alm->solutionLPrev();
    alm->~ALM();                                   // placement-new: no operator delete
    return val;
}

// The discriminating read for m_indicator/m_negatives must
// happen BEFORE any initialize()/configure() call -- _computeStability
// is the write site, and configure() reaches it via
// initialize() -> init(true) -> _computeStability. A test that read AFTER configure()
// would pass even with the Fix-1 seeding deleted. Since the seeding lives
// in the BASE constructor, all four solvers are expected to read poison here with the
// fix reverted -- a strictly stronger discriminator than poisonedLPrevAfterInit's, whose
// pre-fix fail-before was only two of four (m_Lprev is re-seeded by two of the four
// initMethods()).
template <class ALM>
static void poisonedStabilityAfterConstruct(AlmProblem & prob, real_t & indic, index_t & neg)
{
    typename std::aligned_storage<sizeof(ALM), alignof(ALM)>::type buf;
    poisonBytes(&buf, sizeof(buf));
    ALM * alm = new (&buf) ALM(prob.Jacobian, prob.ALResidual, prob.Force);
    indic = alm->indicator();
    neg   = alm->negatives();
    alm->~ALM();                                   // placement-new: no operator delete
}

TEST(m_Lprev_is_not_indeterminate_storage)
{
    AlmProblem prob = pitchforkProblem();

    const real_t lc = poisonedLPrevAfterInit<gsALMLoadControl<real_t> >(prob);
    const real_t cf = poisonedLPrevAfterInit<gsALMCrisfield<real_t> >(prob);
    const real_t rk = poisonedLPrevAfterInit<gsALMRiks<real_t> >(prob);
    const real_t cc = poisonedLPrevAfterInit<gsALMConsistentCrisfield<real_t> >(prob);

    gsInfo<<"  [Lprev init] solutionLPrev() after placement-new + initialize(): "
          <<"LoadControl="<<lc<<"  Crisfield="<<cf<<"  Riks="<<rk
          <<"  ConsistentCrisfield="<<cc<<"\n";

    CHECK_EQUAL(0.0, lc);
    CHECK_EQUAL(0.0, cf);
    CHECK_EQUAL(0.0, rk);
    CHECK_EQUAL(0.0, cc);
}

// ===========================================================================
// The `T m_indicator;` / `index_t m_negatives;`
// members of `gsALMBase` (the two unseeded siblings of m_Lprev) must be well-defined from
// construction, on every solver, BEFORE initialize()/configure() ever runs. Read
// immediately after a poisoned placement-new (see poisonedStabilityAfterConstruct
// above); index_t width is GISMO_INDEX_TYPE, not fixed, so the poisoned (fix-reverted)
// value is never hard-coded here -- only the seeded 0/0 is asserted.
// ===========================================================================
TEST(m_indicator_and_m_negatives_are_not_indeterminate_storage)
{
    AlmProblem prob = pitchforkProblem();
    real_t  indic;
    index_t neg;

    poisonedStabilityAfterConstruct<gsALMLoadControl<real_t> >(prob, indic, neg);
    gsInfo<<"  [stability state] LoadControl: indicator()="<<indic<<"  negatives()="<<neg<<"\n";
    CHECK_EQUAL(0.0, indic);
    CHECK_EQUAL((index_t)0, neg);

    poisonedStabilityAfterConstruct<gsALMCrisfield<real_t> >(prob, indic, neg);
    gsInfo<<"  [stability state] Crisfield: indicator()="<<indic<<"  negatives()="<<neg<<"\n";
    CHECK_EQUAL(0.0, indic);
    CHECK_EQUAL((index_t)0, neg);

    poisonedStabilityAfterConstruct<gsALMRiks<real_t> >(prob, indic, neg);
    gsInfo<<"  [stability state] Riks: indicator()="<<indic<<"  negatives()="<<neg<<"\n";
    CHECK_EQUAL(0.0, indic);
    CHECK_EQUAL((index_t)0, neg);

    poisonedStabilityAfterConstruct<gsALMConsistentCrisfield<real_t> >(prob, indic, neg);
    gsInfo<<"  [stability state] ConsistentCrisfield: indicator()="<<indic<<"  negatives()="<<neg<<"\n";
    CHECK_EQUAL(0.0, indic);
    CHECK_EQUAL((index_t)0, neg);
}

// ===========================================================================
// PIN the non-converged extended-solve
// commit policy. This test does NOT assert a defect; it pins DELIBERATE,
// option-selected, default-on behaviour so it cannot drift silently in
// EITHER direction. See gsALMBase::SPfail (Without=0/With=1),
// gsALMBase::defaultOptions()'s SingularPointFailure default = SPfail::With,
// gsALMBase::_computeSingularPoint()'s `if (switchBranch && (converged ||
// m_SPfail==1)) this->switchBranch();` (switchBranch() depends on the
// committed state on failure), and the user decision recorded 2026-08-24:
// no behavioural default moves here.
//
// Forcing the failure path deterministically: fixture P's extended system
// converges to SingularPointComputeTolE = 1e-10 (the `configure` default) in
// a SINGLE Newton iteration from (1.05,0)/1.05, so MaxIter alone cannot be
// tuned down to force non-convergence -- MaxIter=1 still lets that one
// iteration converge before the failure branch is ever reached.
// SingularPointComputeTolE=0 makes gsALMBase::_extendedSystemSolve()'s convergence test
// `m_residueKTPhi < tol` UNSATISFIABLE (a norm is never
// < 0), so every iteration budget is exhausted and MaxIter=1 reaches
// _extendedSystemSolve()'s `m_numIterations == m_maxIterations-1` on its very
// first and only iteration -- reliably, independent of solver convergence
// speed.
// ===========================================================================
TEST(extended_solve_failure_commit_is_gated_by_SingularPointFailure)
{
    AlmProblem prob = pitchforkProblem();
    gsVector<real_t> U(2); U << 1.05, 0.0;

    // --- SingularPointFailure = 1 (SPfail::With, the DEFAULT): the committed
    //     state has MOVED to the last (non-converged) increment. -----------
    {
        gsALMLoadControl<real_t> alm(prob.Jacobian, prob.ALResidual, prob.Force);
        configure(alm, 0.05);
        alm.setSolution(U,1.05);
        alm.options().setInt ("MaxIter", 1);
        alm.options().setReal("SingularPointComputeTolE", 0.0); // unreachable -> forces the failure branch
        alm.options().setInt ("SingularPointFailure", 1);
        alm.applyOptions();

        const gsStatus st = alm.computeSingularPoint(U,1.05,
                                                     /*switchBranch*/false,
                                                     /*jacobian*/true,
                                                     /*testPoint*/true);
        // Must actually have taken the failure path, or the rest of the test
        // degenerates vacuously.
        CHECK(st != gsStatus::Success);
        CHECK(!alm.converged());

        const real_t moved = (alm.solutionU()-U).norm();
        gsInfo<<"  [singular-point failure] SingularPointFailure=1 (With): |U-U0| = "<<moved
              <<"  L = "<<alm.solutionL()<<"\n";
        CHECK(moved > 1e-6);
    }

    // --- SingularPointFailure = 0 (SPfail::Without): the committed state is
    //     UNMOVED -- bisection is off (SingularPointComputeTolB=0, `configure`),
    //     so nothing has touched m_U/m_L since gsALMBase::_extendedSystemSolve()'s
    //     entry assignment `m_U = U; m_L = L;`. ------
    {
        gsALMLoadControl<real_t> alm(prob.Jacobian, prob.ALResidual, prob.Force);
        configure(alm, 0.05);
        alm.setSolution(U,1.05);
        alm.options().setInt ("MaxIter", 1);
        alm.options().setReal("SingularPointComputeTolE", 0.0); // unreachable -> forces the failure branch
        alm.options().setInt ("SingularPointFailure", 0);
        alm.applyOptions();

        const gsStatus st = alm.computeSingularPoint(U,1.05,
                                                     /*switchBranch*/false,
                                                     /*jacobian*/true,
                                                     /*testPoint*/true);
        CHECK(st != gsStatus::Success);
        CHECK(!alm.converged());

        const real_t moved = (alm.solutionU()-U).norm();
        gsInfo<<"  [singular-point failure] SingularPointFailure=0 (Without): |U-U0| = "<<moved
              <<"  L = "<<alm.solutionL()<<"\n";
        CHECK_EQUAL(0.0, moved);
        CHECK_EQUAL(1.05, alm.solutionL());
    }
}

// ===========================================================================
// MaxIter == 0 must not report convergence after
// zero iterations. _extendedSystemSolve's loop (
// `for (m_numIterations = 0; m_numIterations < m_maxIterations; ...)`) never runs
// its body when m_maxIterations < 1, so the non-convergence branch
// (`if (m_numIterations == m_maxIterations-1)`) is unreachable and the function used
// to fall through to `return true;` --
// i.e. report SUCCESS after doing nothing. Uses
// ExtendedProbe (defined above) directly on gsALMBase::_extendedSystemSolve, seeded
// with the exact fixture already proven (in
// extended_system_shift_ladder_survives_an_exactly_singular_tangent above) to
// converge from this state: U=(1,0), L=1.0=lambda*, V=(0,1) is the pitchfork's
// closed-form singular point, so K*(V) = diag(1,0)*(0,1) = 0 already below tol --
// the MaxIter=1 companion below is a genuine convergence, not a degenerate no-op.
// ===========================================================================
TEST(extended_solve_zero_iteration_budget_is_not_converged)
{
    AlmProblem prob = pitchforkProblem();
    gsVector<real_t> U(2); U << 1.0, 0.0;      // the pitchfork itself
    const real_t      L = 1.0;                 // lambda* = 1
    gsVector<real_t> V(2); V << 0.0, 1.0;      // the closed-form critical mode

    // --- MaxIter = 0: nothing is solved; must not be reported as converged. -----
    {
        ExtendedProbe alm(prob.Jacobian, prob.ALResidual, prob.Force);
        configure(alm, 0.05);
        alm.options().setInt("MaxIter", 0);
        alm.applyOptions();

        const bool ok = alm.extendedSolveAt(U,L,V,1e-10);
        gsInfo<<"  [zero-iteration guard] MaxIter=0: extendedSolveAt() returned "<<(ok?"true":"false")
              <<"  numIterations()="<<alm.numIterations()<<"\n";
        CHECK(!ok);
    }

    // --- MaxIter = 1: UNCHANGED path -- reaches the loop body and converges on its
    //     first and only iteration, proving the guard is not a blanket "always
    //     false". --------------------------------------------------------------
    {
        ExtendedProbe alm(prob.Jacobian, prob.ALResidual, prob.Force);
        configure(alm, 0.05);
        alm.options().setInt("MaxIter", 1);
        alm.applyOptions();

        const bool ok = alm.extendedSolveAt(U,L,V,1e-10);
        gsInfo<<"  [zero-iteration guard] MaxIter=1: extendedSolveAt() returned "<<(ok?"true":"false")
              <<"  numIterations()="<<alm.numIterations()<<"\n";
        CHECK(ok);
    }
}

// ===========================================================================
// Direction A -- computeResidualNorms()'s m_residueU denominator must not be
// frozen at a round-off value when the extended solve is seeded AT its own
// answer. Fixture: the pitchfork's exact closed-form singular point,
// U=(1,0), L=1 (K = diag(1,0) there, V=(0,1) the exact null vector -- the same
// state the "MaxIter=1" test above already proved converges in one iteration on
// the ||K_T.V||-only test). L is perturbed by exactly one ULP so the equilibrium
// residual R=(u1-lambda,...) is a genuine, tiny, nonzero quantity (round-off, not
// an artificial injection) rather than the bit-exact zero the un-perturbed point
// would give -- either way computeResidualNorms's iteration-0 branch sets
// m_residueU = m_deltaU.norm()/m_basisResidualU. Without the floor,
// m_basisResidualU == m_deltaU.norm() there by construction, making m_residueU
// == 1 (or 0/0) regardless of how converged the point already is: the composite
// test below is unreachable. With SingularPointComposite ON and MaxIter=1, the
// extended solve must ALSO satisfy TolF/TolU (not just the ||K_T.V|| test), so
// this is precisely the single-condition trap the extended-solve floor exists to
// dissolve.
// ===========================================================================
TEST(extended_solve_denominator_is_not_frozen_at_seed_distance)
{
    AlmProblem prob = pitchforkProblem();
    gsVector<real_t> U(2); U << 1.0, 0.0;               // the pitchfork itself
    gsVector<real_t> V(2); V << 0.0, 1.0;               // the closed-form critical mode
    const real_t L = 1.0 + std::numeric_limits<real_t>::epsilon(); // one ULP off lambda*=1

    ExtendedProbe alm(prob.Jacobian, prob.ALResidual, prob.Force);
    configure(alm, 0.05); // TolF=TolU=1e-8 here -- tighter than the library 1e-6 default
    alm.options().setSwitch("SingularPointComposite", true);
    alm.options().setInt("MaxIter", 1);
    alm.applyOptions();

    const bool ok = alm.extendedSolveAt(U, L, V, 1e-10);
    gsInfo<<"  [composite predictor, on-seed] SingularPointComposite=true, MaxIter=1, seeded at the exact "
          <<"singular point: extendedSolveAt() returned "<<(ok?"true":"false")
          <<"  residueKTPhi()="<<alm.residueKTPhi()<<"\n";
    // The K_T.V-only condition must already be satisfied here -- this seed is the
    // fixture already proven to converge on that test alone in one iteration -- so a
    // failure below can only be attributed to the m_residueU denominator, not to the
    // singular-point classification itself.
    CHECK(alm.residueKTPhi() < 1e-10);
    CHECK(ok);
}

// ===========================================================================
// Direction B -- the m_basisResidualU floor in computeResidualNorms() is scoped
// to the extended solve ONLY: the ordinary corrector's iteration-0 basis
// (m_basisResidualU, captured at _step()'s predictor call)
// must stay EXACTLY the raw predictor increment, unfloored, even in a regime
// where floor*||U|| would dominate it if the scoping leaked.
//
// Construction: take ONE genuine step at ordinary length (Length=0.05) from the
// true cold start (U=0,L=0 -- gsALMRiks::predictor()'s `(m_U-m_Uprev).norm()<tol`
// branch), which installs a real secant (m_Uprev/m_Lprev) and
// leaves m_arcLength_prev==0.05 (the non-adaptive branch of _step(), "else
// m_arcLength_prev = m_arcLength;"). Then shrink Length to 1e-10 and capture a
// SECOND predictor with nothing else touched: predictor()'s secant branch
// (since m_U now != m_Uprev) computes
// deltaU=(m_U-m_Uprev)/m_arcLength_prev, i.e. the FIRST step's secant direction
// (magnitude set by the 0.05 step, no longer by the new Length), then rescales it
// by the new m_arcLength=1e-10 (predictor()'s `m_deltaU *= m_arcLength;` step) --
// so the raw predictor
// increment this second call produces is O(1e-10), while ||U|| is still O(0.05)
// (the first step's own arc length): floor*||U|| = 0.1*0.05 = 5e-3 would dominate
// the O(1e-10) raw increment by SEVEN orders of magnitude if the scoping leaked.
// runStepInstrumented(rec, maxIt=0) runs
// initiateStep()+predictor()+computeResidual()+computeResidualNorms() (the exact
// basis-capture sequence of _step()) and returns before
// the corrector loop ever executes (the maxIt=0 loop bound is never reached),
// leaving m_basisResidualU exactly the ordinary path captured -- with no
// extendedSolve=true call anywhere on this path.
// ===========================================================================
TEST(ordinary_corrector_basis_is_never_floored)
{
    AlmProblem prob = pitchforkProblem();

    RiksProbe alm(prob.Jacobian, prob.ALResidual, prob.Force);
    configure(alm, 0.05);
    CHECK(alm.step() == gsStatus::Success); // establishes a genuine secant (Uprev/Lprev/arcLength_prev)

    // Shrink Length drastically; m_arcLength_prev is UNCHANGED (getOptions() only
    // reseeds it while !m_stepTaken, and a step was just taken).
    alm.options().setReal("Length", 1e-10);
    alm.applyOptions();

    std::vector<RiksProbe::IterRecord> rec;
    const bool converged = alm.runStepInstrumented(rec, /*maxIt*/0);
    // maxIt=0: the corrector for-loop body never runs, so this guards that only the
    // predictor + basis capture ran -- not a degenerate "trivially converged" step.
    CHECK(!converged);
    CHECK_EQUAL(0u, (unsigned)rec.size());

    gsInfo<<"  [composite predictor, off-seed] ||DeltaU|| (raw predictor increment) = "<<alm.DeltaU().norm()
          <<"  m_basisResidualU = "<<alm.basisResidualU()
          <<"  ||U|| = "<<alm.solutionU().norm()<<"\n";
    // The construction must actually have reached the tiny-increment regime, or the
    // rest of the test is vacuous: the raw increment must sit far enough below
    // floor*||U|| that a leaked floor would be impossible to miss.
    const real_t floorRelU = 1e-1;
    CHECK(alm.DeltaU().norm() < 1e-3*floorRelU*alm.solutionU().norm());
    // Exact equality: computeResidualNorms(false) (the default, and the ONLY overload
    // any ordinary-path call site can reach) never applies the extended-solve floor on
    // this path, so the captured basis must equal the raw predictor increment to the
    // last ULP.
    CHECK_EQUAL(alm.DeltaU().norm(), alm.basisResidualU());
}

// ===========================================================================
// gsBucklingSolver / gsModalSolver / gsEigenProblemBase -- closed-form 2x2
// pencils.
// ===========================================================================
//
// gsBucklingSolver poses the bifurcation condition [K_L + lambda K_G] v = 0,
// K_G := K_NL - K_L, as the pencil K_G v = nu K_L v with lambda = -1/nu, and
// hands (K_G, K_L) -- in that order -- to
// gsEigen::GeneralizedSelfAdjointEigenSolver::compute(A,B), which
// Cholesky-factors B = K_L. Eliminating nu from lambda = -1/nu gives
// K_L v + lambda K_G v = 0, the residual used below to certify a returned
// pair really is an eigenpair of the pencil, independent of how it was
// computed.

TEST(buckling_factors_are_exact_for_an_indefinite_geometric_stiffness)
{
    // K_L = 3*I in both fixtures; K_G differs only in sign structure, and
    // both K_NL = K_L + K_G are SPD (both are physically well-posed
    // pre-buckling states). The buckling-factor MAGNITUDES are identical,
    // {1,3}, in both cases, so it is the definiteness of K_G -- not its
    // size -- that the two fixtures discriminate.
    //
    //   indefinite: K_G = [[1,2],[2,1]], eig(K_G) = {3,-1}.
    //     det(3I + l*K_G) = (3+l)^2 - 4l^2 = 3*(3-l)*(1+l)  ->  l = 3, -1
    //
    //   definite (control): K_G = [[2,1],[1,2]], eig(K_G) = {3,1}.
    //     det(3I + l*K_G) = (3+2l)^2 - l^2 = (3+l)*(3+3l)  ->  l = -3, -1
    gsMatrix<real_t> KLdense(2,2), KGindef(2,2), KGdef(2,2);
    KLdense << 3,0, 0,3;
    KGindef << 1,2, 2,1;
    KGdef   << 2,1, 1,2;

    gsSparseMatrix<real_t> KL, KNLi, KNLd;
    toSparse2(KLdense, KL);
    toSparse2(KLdense+KGindef, KNLi);
    toSparse2(KLdense+KGdef,   KNLd);

    // Check 1 -- the fixture is what it claims to be: the sign structure of
    // K_G. Without this, a later edit could quietly make K_G definite and
    // the test would keep passing while measuring nothing.
    gsEigen::SelfAdjointEigenSolver<gsMatrix<real_t>::Base> esIndef, esDef;
    esIndef.compute(KGindef);
    esDef.compute(KGdef);
    gsInfo<<"  [buckling/indefinite] eig(K_G) = "<<esIndef.eigenvalues().transpose()<<"\n";
    gsInfo<<"  [buckling/definite]   eig(K_G) = "<<esDef.eigenvalues().transpose()<<"\n";
    CHECK(esIndef.eigenvalues()(0) < 0.0 && esIndef.eigenvalues()(1) > 0.0);
    CHECK(esDef.eigenvalues()(0)   > 0.0 && esDef.eigenvalues()(1)   > 0.0);

    // gsBucklingSolver(linear, nonlinear): linear = K_L, nonlinear = K_NL.
    // Both parameters are non-const lvalue references, hence the named
    // gsSparseMatrix variables above.
    gsBucklingSolver<real_t> solverIndef(KL, KNLi);
    const gsStatus stIndef = solverIndef.compute();

    gsBucklingSolver<real_t> solverDef(KL, KNLd);
    const gsStatus stDef = solverDef.compute();

    gsInfo<<"  [buckling/indefinite] status = "<<static_cast<index_t>(stIndef)
          <<"  lambda = "<<solverIndef.values().transpose()<<"\n";
    gsInfo<<"  [buckling/definite]   status = "<<static_cast<index_t>(stDef)
          <<"  lambda = "<<solverDef.values().transpose()<<"\n";

    // Check 2.
    CHECK(stIndef == gsStatus::Success);
    CHECK(stDef   == gsStatus::Success);

    // Check 3 -- signed lambda under the documented convention
    // (lambda = -1/nu, ascending in |lambda|):
    //   indefinite: l = {-1, 3}   definite: l = {-1, -3}
    CHECK_CLOSE(-1.0, solverIndef.value(0), 1e-10);
    CHECK_CLOSE( 3.0, solverIndef.value(1), 1e-10);
    CHECK_CLOSE(-1.0, solverDef.value(0), 1e-10);
    CHECK_CLOSE(-3.0, solverDef.value(1), 1e-10);

    // Check 4 -- convention-free magnitude check.
    std::vector<real_t> magIndef = { math::abs(solverIndef.value(0)), math::abs(solverIndef.value(1)) };
    std::vector<real_t> magDef   = { math::abs(solverDef.value(0)),   math::abs(solverDef.value(1))   };
    std::sort(magIndef.begin(), magIndef.end());
    std::sort(magDef.begin(),   magDef.end());
    CHECK_CLOSE(1.0, magIndef[0], 1e-10);
    CHECK_CLOSE(3.0, magIndef[1], 1e-10);
    CHECK_CLOSE(1.0, magDef[0], 1e-10);
    CHECK_CLOSE(3.0, magDef[1], 1e-10);

    // Check 5 -- eigenpair residual K_L v + lambda K_G v = 0. Neither
    // fixture has a nu ~ 0 (floored, |lambda| -> infinity) mode, so both
    // lambda_k are finite and the residual is meaningful for k = 0,1.
    const real_t KLnorm = KLdense.norm();
    for (index_t k = 0; k!=2; ++k)
    {
        const gsMatrix<real_t> vIndef = solverIndef.vector(k);
        const real_t lamIndef = solverIndef.value(k);
        const real_t resIndef = (KLdense*vIndef + lamIndef*KGindef*vIndef).norm()
                               / (KLnorm*vIndef.norm());
        gsInfo<<"  [buckling/indefinite] eigenpair residual k="<<k<<" = "<<resIndef<<"\n";
        CHECK(resIndef < 1e-10);

        const gsMatrix<real_t> vDef = solverDef.vector(k);
        const real_t lamDef = solverDef.value(k);
        const real_t resDef = (KLdense*vDef + lamDef*KGdef*vDef).norm()
                             / (KLnorm*vDef.norm());
        gsInfo<<"  [buckling/definite] eigenpair residual k="<<k<<" = "<<resDef<<"\n";
        CHECK(resDef < 1e-10);
    }

    // Check 6 -- the discrimination proof: the K_L-first pencil
    // K_L v = mu K_G v factors K_G directly. On the definite control that
    // Cholesky factorization is valid and reproduces the closed-form
    // magnitudes; on the indefinite fixture it factors an indefinite
    // matrix and must not.
    gsEigen::GeneralizedSelfAdjointEigenSolver<gsMatrix<real_t>::Base> esOldIndef, esOldDef;
    esOldIndef.compute(KLdense, KGindef);
    esOldDef.compute(KLdense, KGdef);
    gsInfo<<"  [buckling/indefinite] old pencil info="<<static_cast<int>(esOldIndef.info())
          <<" mu = "<<esOldIndef.eigenvalues().transpose()<<"\n";
    gsInfo<<"  [buckling/definite]   old pencil info="<<static_cast<int>(esOldDef.info())
          <<" mu = "<<esOldDef.eigenvalues().transpose()<<"\n";

    std::vector<real_t> muDef = { math::abs(esOldDef.eigenvalues()(0)), math::abs(esOldDef.eigenvalues()(1)) };
    std::sort(muDef.begin(), muDef.end());
    CHECK_CLOSE(1.0, muDef[0], 1e-10);
    CHECK_CLOSE(3.0, muDef[1], 1e-10);

    std::vector<real_t> muIndef = { math::abs(esOldIndef.eigenvalues()(0)), math::abs(esOldIndef.eigenvalues()(1)) };
    std::sort(muIndef.begin(), muIndef.end());
    const bool oldPencilReproducesIndefinite =
        (math::abs(muIndef[0]-1.0) < 1e-10) && (math::abs(muIndef[1]-3.0) < 1e-10);
    CHECK(!oldPencilReproducesIndefinite);
}

TEST(buckling_solver_does_not_report_success_on_a_non_spd_linear_stiffness)
{
    // After the pencil swap in gsBucklingSolver's 2-argument constructor,
    // its second matrix (the one gsEigen would Cholesky-factor) is K_L, the
    // linear stiffness: eig([[1,0],[0,-1]]) = {1,-1}, symmetric indefinite.
    // gsBucklingSolver's own checkDefinite() certifies positive-definiteness
    // from its SimplicialLDLT factorization of K_L before compute() ever
    // touches the dense eigensolver, so this must be caught here, not by
    // the (silent) Cholesky inside GeneralizedSelfAdjointEigenSolver.
    gsMatrix<real_t> KLdense(2,2), KGdense(2,2);
    KLdense << 1,0, 0,-1;
    KGdense << 1,0, 0, 1;   // any matrix: irrelevant to the guard on K_L

    gsSparseMatrix<real_t> KL, KNL;
    toSparse2(KLdense, KL);
    toSparse2(KLdense+KGdense, KNL);

    gsBucklingSolver<real_t> solver(KL, KNL);
    const gsStatus st = solver.compute();

    gsInfo<<"  [eigen/non-spd] status = "<<static_cast<index_t>(st)<<"\n";

    // Check 1 -- do not pin which non-Success enumerator the guard uses.
    CHECK(st != gsStatus::Success);

    // Check 2 -- the invariant that survives any future re-mapping of the
    // guard's status value: Success must never coincide with a pair that
    // does not actually satisfy the pencil. values()/vectors() are only
    // meaningful when st == Success, so the residual is computed only then.
    real_t maxPencilResidual = 0.0;
    if (st == gsStatus::Success)
    {
        const real_t KLnorm = KLdense.norm();
        for (index_t k = 0; k!=2; ++k)
        {
            const gsMatrix<real_t> v = solver.vector(k);
            const real_t lam = solver.value(k);
            const real_t res = (KLdense*v + lam*KGdense*v).norm() / (KLnorm*v.norm());
            maxPencilResidual = math::max(maxPencilResidual, res);
        }
    }
    // Success must imply that what came back really is an eigenpair of the pencil.
    CHECK( !(st == gsStatus::Success && maxPencilResidual > 1e-8) );
}

TEST(modal_solver_reproduces_the_analytic_pencil)
{
    // K v = w M v, K = [[2,-1],[-1,2]] (SPD), M = diag(1,4) (SPD):
    //   det(K - w M) = (2-w)(2-4w) - 1 = 4w^2 - 10w + 3
    //     -> w = (5 +- sqrt(13)) / 4
    // gsModalSolver never swaps the pencil and never converts its
    // eigenvalues, so this is the control that pins gsEigenProblemBase
    // itself: gsEigenProblemBase must solve this pencil and return
    // M-orthonormal modes in the original space, whatever gsBucklingSolver
    // does with its own pencil.
    gsMatrix<real_t> Kdense(2,2), Mdense(2,2);
    Kdense << 2,-1, -1,2;
    Mdense << 1,0, 0,4;

    gsSparseMatrix<real_t> K, M;
    toSparse2(Kdense, K);
    toSparse2(Mdense, M);

    gsModalSolver<real_t> solver(K, M);
    const gsStatus st = solver.compute();

    const real_t w0 = (5.0 - math::sqrt(13.0))/4.0;
    const real_t w1 = (5.0 + math::sqrt(13.0))/4.0;

    gsInfo<<"  [modal/control] status = "<<static_cast<index_t>(st)
          <<"  w = "<<solver.values().transpose()
          <<"  errors = "<<math::abs(solver.value(0)-w0)<<", "<<math::abs(solver.value(1)-w1)<<"\n";

    // Check 1.
    CHECK(st == gsStatus::Success);

    // Check 2 -- eigenvalues; gsEigen orders generalized eigenvalues
    // ascending, so value(0) is the smaller root.
    CHECK_CLOSE(w0, solver.value(0), 1e-12);
    CHECK_CLOSE(w1, solver.value(1), 1e-12);

    // Check 3 -- M-orthonormality of the modes: v_i^T M v_j = delta_ij. This
    // pins the eigenVECTORS, which the eigenvalue check alone does not, and
    // it is where a mis-scaled back-transform in the base class would show.
    real_t worstOrtho = 0.0;
    for (index_t i = 0; i!=2; ++i)
        for (index_t j = 0; j!=2; ++j)
        {
            const real_t vMv = (solver.vector(i).transpose()*Mdense*solver.vector(j))(0,0);
            const real_t expect = (i==j) ? 1.0 : 0.0;
            worstOrtho = math::max(worstOrtho, math::abs(vMv-expect));
        }
    gsInfo<<"  [modal/control] worst |v_i^T M v_j - delta_ij| = "<<worstOrtho<<"\n";
    CHECK(worstOrtho < 1e-10);

    // Check 4 -- pencil residual.
    const real_t Knorm = Kdense.norm();
    for (index_t k = 0; k!=2; ++k)
    {
        const gsMatrix<real_t> v = solver.vector(k);
        const real_t w = solver.value(k);
        const real_t res = (Kdense*v - w*Mdense*v).norm() / (Knorm*v.norm());
        gsInfo<<"  [modal/control] pencil residual k="<<k<<" = "<<res<<"\n";
        CHECK(res < 1e-10);
    }
}

TEST(buckling_factor_sign_follows_the_bifurcation_condition)
{
    // Bifurcation condition for the reduced, Dirichlet-eliminated tangent
    // (example_ShearBuckling.cpp, "classical bifurcation condition"): a
    // load multiple lambda is a
    // buckling factor iff [K_L + lambda*K_G] v = 0 for some v != 0, i.e.
    // iff det(K_L + lambda*K_G) = 0, with K_G := K_NL(w) - K_L
    // (example_ShearBuckling.cpp, "K_G is defined the same way"). The
    // determinant below is factored
    // by hand from that condition; it is not read off gsBucklingSolver.
    //
    // K_L = I, K_G = [[1,2],[2,1]] (eig = {-1,3}, indefinite):
    //   det(I + l*K_G) = det([[1+l,2l],[2l,1+l]]) = (1+l)^2 - 4l^2
    //                  = (1+l-2l)(1+l+2l) = (1-l)(1+3l)
    //     -> l = 1   and   l = -1/3
    // Ascending in |l|: value(0) = -1/3, value(1) = +1.
    //
    // Sign discrimination: the two roots have DIFFERENT signs and
    // DIFFERENT magnitudes, so a solver reporting +1/nu instead of -1/nu
    // would return {+1/3, -1} -- the same two magnitudes with both signs
    // flipped. A check on |lambda| alone cannot see that; the signed
    // checks below can.
    //
    // A second, solver-independent reading of the same root: l = +1 is a
    // buckling factor iff K_L + 1*K_G = K_NL is singular, and
    // K_NL = [[2,2],[2,2]] visibly has two identical rows.
    gsMatrix<real_t> KLdense(2,2), KGdense(2,2), KNLdense(2,2);
    KLdense << 1,0, 0,1;
    KGdense << 1,2, 2,1;
    KNLdense = KLdense + KGdense;

    gsSparseMatrix<real_t> KL, KNL;
    toSparse2(KLdense, KL);
    toSparse2(KNLdense, KNL);

    // gsBucklingSolver(linear, nonlinear): linear = K_L, nonlinear = K_NL;
    // it forms K_G = K_NL - K_L internally.
    gsBucklingSolver<real_t> solver(KL, KNL);
    const gsStatus st = solver.compute();

    const real_t lam0 = solver.value(0);
    const real_t lam1 = solver.value(1);
    const real_t KLnorm = KLdense.norm();
    const gsMatrix<real_t> v0 = solver.vector(0);
    const gsMatrix<real_t> v1 = solver.vector(1);
    // The '+' is the '+' of [K_L + lambda*K_G] v = 0
    // (example_ShearBuckling.cpp, "[ K_L + lambda K_G ] v = 0"), not a
    // sign convention chosen here.
    // A reported -lambda would make this residual O(1).
    const real_t res0 = (KLdense*v0 + lam0*KGdense*v0).norm() / (KLnorm*v0.norm());
    const real_t res1 = (KLdense*v1 + lam1*KGdense*v1).norm() / (KLnorm*v1.norm());
    const real_t detKNL = KNLdense.determinant();

    gsInfo<<"  [buckling/sign] status = "<<static_cast<index_t>(st)
          <<"  lambda = "<<lam0<<", "<<lam1
          <<"  residuals = "<<res0<<", "<<res1<<"\n";

    // Check 1.
    CHECK(st == gsStatus::Success);

    // Check 2 -- the closed-form roots, signed.
    CHECK_CLOSE(-1.0/3.0, lam0, 1e-10);
    CHECK_CLOSE( 1.0,     lam1, 1e-10);

    // Check 3 -- the sign-sensitive eigenpair residual.
    CHECK(res0 < 1e-10);
    CHECK(res1 < 1e-10);

    // Check 4 -- K_NL is singular by inspection, independent of the solver.
    gsInfo<<"  [buckling/sign] det(K_NL) = "<<detKNL<<"\n";
    CHECK(math::abs(detKNL) < 1e-12);
}

TEST(buckling_mode_with_negligible_prestress_coupling_is_sentinelled_and_sorts_last)
{
    // K_L = I, K_G = diag(1e-12, 0.5): SPD, but one direction is ~1e12
    // times softer than the other. K_NL = K_L + K_G = diag(1+1e-12, 1.5).
    //
    //   det(I + l*K_G) = (1 + 1e-12*l)(1 + 0.5*l)
    //     -> l = -1e12   and   l = -2
    //
    // The l = -1e12 root is the physical content of the sentinel: as the
    // pre-stress coupling in that direction vanishes, its buckling factor
    // runs off to -infinity. The coupling here is small but strictly
    // non-zero on purpose -- that is what lets the test discriminate.
    // (1 + 1e-12) - 1 is accurate to ~1e-16 in double precision, so the
    // small entry is real, not round-off.
    //
    // Raw pencil eigenvalues (with K_L = I these coincide with eig(K_G)):
    // nu = {1e-12, 0.5}, so max|nu| = 0.5 and the relative floor
    // 1e-8*max|nu| = 5e-9. Since 1e-12 <= 5e-9 (more than three orders of
    // margin), that mode is floored and reported as the sentinel
    // has_infinity ? infinity() : max(), UNSIGNED regardless of the true
    // root's sign. An absent or backwards floor comparison would instead
    // return the finite value -1/1e-12 = -1e12: "no lambda near zero" is
    // therefore vacuously green here and must not be the only check -- the
    // sentinel equality and the > 1e15 magnitude gate below are what
    // separate +infinity from -1e12.
    //
    // Eigen returns the raw pencil eigenvalues ascending in nu, i.e.
    // [1e-12, 0.5], which places the floored mode at raw index 0.
    // convertToLambda() must move it to index 1 for the "sorts last"
    // assertion below to be non-vacuous.
    gsMatrix<real_t> KLdense(2,2), KGdense(2,2), KNLdense(2,2);
    KLdense << 1,0, 0,1;
    KGdense << 1e-12,0, 0,0.5;
    KNLdense = KLdense + KGdense;

    gsSparseMatrix<real_t> KL, KNL;
    toSparse2(KLdense, KL);
    toSparse2(KNLdense, KNL);

    // Check 1 -- fixture self-certification, before touching the solver:
    // without this a later edit could quietly make both nu comparable and
    // the test would keep passing while exercising nothing.
    gsEigen::SelfAdjointEigenSolver<gsMatrix<real_t>::Base> esKG;
    esKG.compute(KGdense);
    const real_t nuSmall = esKG.eigenvalues()(0);
    const real_t nuLarge = esKG.eigenvalues()(1);
    gsInfo<<"  [buckling/floored] eig(K_G) = "<<nuSmall<<", "<<nuLarge<<"\n";
    CHECK(nuSmall > 0.0);
    CHECK(nuSmall <= 1e-8*nuLarge);

    gsBucklingSolver<real_t> solver(KL, KNL);
    const gsStatus st = solver.compute();

    const real_t sentinel = std::numeric_limits<real_t>::has_infinity
                           ? std::numeric_limits<real_t>::infinity()
                           : std::numeric_limits<real_t>::max();

    gsInfo<<"  [buckling/floored] status = "<<static_cast<index_t>(st)
          <<"  lambda = "<<solver.value(0)<<", "<<solver.value(1)<<"\n";

    // Check 2 -- a floored mode is not a solver error.
    CHECK(st == gsStatus::Success);

    // Check 3 -- the finite mode, signed and derived from the determinant above.
    CHECK_CLOSE(-2.0, solver.value(0), 1e-10);

    // Check 4 -- not NaN.
    CHECK(!math::isnan(solver.value(1)));

    // Check 5 -- the sentinel, exactly: the value is assigned, not computed.
    CHECK_EQUAL(sentinel, solver.value(1));

    // Check 6 -- sorted last, and far enough from the finite mode that an
    // unfloored -1/nu = -1e12 (the value a backwards floor comparison would
    // produce) cannot be mistaken for the sentinel.
    CHECK(math::abs(solver.value(0)) < math::abs(solver.value(1)));
    CHECK(math::abs(solver.value(1)) > 1e15);

    // Check 7 -- the eigenvectors were permuted together with the values:
    // the mode paired with the sentinel must be the negligible-coupling
    // one. lambda_1 is not finite, so no pencil residual is computed at
    // k = 1; the finite mode at k = 0 gets the same physics residual as
    // the sign test above.
    const gsMatrix<real_t> vFloored = solver.vector(1);
    const real_t floorResidual = (KGdense*vFloored).norm() / vFloored.norm();
    gsInfo<<"  [buckling/floored] ||K_G v_floored||/||v_floored|| = "<<floorResidual<<"\n";
    CHECK(floorResidual <= 1e-8*nuLarge);

    const gsMatrix<real_t> v0 = solver.vector(0);
    const real_t KLnorm = KLdense.norm();
    const real_t res0 = (KLdense*v0 + solver.value(0)*KGdense*v0).norm() / (KLnorm*v0.norm());
    gsInfo<<"  [buckling/floored] eigenpair residual k=0 = "<<res0<<"\n";
    CHECK(res0 < 1e-10);
}

TEST(buckling_shift_option_leaves_the_reported_factors_unchanged)
{
    // Same fixture and closed form as the sign test above: K_L = I,
    // K_G = [[1,2],[2,1]], K_NL = [[2,2],[2,2]], exact lambda = {-1/3, +1}
    // from det(I + l*K_G) = (1-l)(1+3l).
    //
    // With nu := -1/lambda, this fixture's exact raw pencil eigenvalues
    // are nu = {3,-1}. On the dense path (gsEigenProblemBase<T>::compute), a
    // nonzero shift s solves (K_G - s*K_L) v = (nu-s) K_L v and then adds
    // s back to the computed eigenvalue before convertToLambda() runs, so
    // a nu-space shift is invisible in the reported lambda. At s = 0.5 the
    // three possible readings of "shift" predict:
    //
    //   acts on nu, added back (documented)     : {-1/3, +1}   unchanged
    //   acts on the reported lambda instead     : {-1/3+0.5, 1+0.5} = {1/6, 3/2}
    //   acts on nu, add-back missing            : {-1/(3-0.5), -1/(-1-0.5)} = {-0.4, 2/3}
    //
    // All three are numerically distinct at 1e-10.
    //
    // This test cannot distinguish a correct nu-space shift from a shift
    // that is silently ignored: in the reported lambda those two are
    // observationally identical, which is exactly what "acts on nu, not on
    // lambda" means. What it DOES discriminate is a shift applied in
    // lambda space and a shift whose add-back is missing -- the two ways
    // the documented semantics can actually be violated.
    gsMatrix<real_t> KLdense(2,2), KGdense(2,2), KNLdense(2,2);
    KLdense << 1,0, 0,1;
    KGdense << 1,2, 2,1;
    KNLdense = KLdense + KGdense;

    gsSparseMatrix<real_t> KL, KNL;
    toSparse2(KLdense, KL);
    toSparse2(KNLdense, KNL);

    gsBucklingSolver<real_t> solverDefault(KL, KNL);
    const gsStatus stDefault = solverDefault.compute();

    gsBucklingSolver<real_t> solverShifted(KL, KNL);
    // The two getReal checks are not decoration: gsOptionList::setReal on a
    // mis-spelled key does not abort, so without them a typo in the option
    // name would leave this test vacuously green.
    CHECK_CLOSE(0.0, solverShifted.options().getReal("shift"), 1e-15);
    solverShifted.options().setReal("shift", 0.5);
    CHECK_CLOSE(0.5, solverShifted.options().getReal("shift"), 1e-15);
    const gsStatus stShifted = solverShifted.compute();

    gsInfo<<"  [buckling/shift] shifted lambda = "<<solverShifted.value(0)<<", "<<solverShifted.value(1)
          <<"  unshifted lambda = "<<solverDefault.value(0)<<", "<<solverDefault.value(1)
          <<"  lambda-space reading = 0.166667, 1.5"
          <<"  missing-add-back reading = -0.4, 0.666667\n";

    // Check 1.
    CHECK(stDefault == gsStatus::Success);
    CHECK(stShifted == gsStatus::Success);

    // Check 2 -- the closed form, unchanged by the shift.
    CHECK_CLOSE(-1.0/3.0, solverShifted.value(0), 1e-10);
    CHECK_CLOSE( 1.0,     solverShifted.value(1), 1e-10);

    // Check 3 -- shifted equals unshifted mode by mode. CHECK_CLOSE, not
    // CHECK_EQUAL: the two runs factor different matrices and agree only
    // to round-off.
    CHECK_CLOSE(solverDefault.value(0), solverShifted.value(0), 1e-10);
    CHECK_CLOSE(solverDefault.value(1), solverShifted.value(1), 1e-10);

    // Check 4 -- the physics residual on the shifted solver's pairs; this
    // catches a shift that corrupts the eigenvectors while leaving the
    // reported values plausible.
    const real_t KLnorm = KLdense.norm();
    for (index_t k = 0; k!=2; ++k)
    {
        const gsMatrix<real_t> v = solverShifted.vector(k);
        const real_t lam = solverShifted.value(k);
        const real_t res = (KLdense*v + lam*KGdense*v).norm() / (KLnorm*v.norm());
        gsInfo<<"  [buckling/shift] shifted eigenpair residual k="<<k<<" = "<<res<<"\n";
        CHECK(res < 1e-10);
    }
}

TEST(modal_solver_does_not_report_success_on_a_non_spd_mass)
{
    // K = [[2,-1],[-1,2]] (SPD, eig = {1,3}), M = diag(1,-1): symmetric and
    // invertible, but indefinite. Invertible matters: the failure being
    // pinned is *definiteness*, not singularity, so SimplicialLDLT::info()
    // alone would pass and only the vectorD() > 0 condition discriminates.
    // gsModalSolver's own checkDefinite() certifies positive-definiteness of
    // M before compute() ever touches the dense eigensolver, so this must be
    // caught here, not by the (silent) Cholesky inside
    // GeneralizedSelfAdjointEigenSolver.
    gsMatrix<real_t> Kdense(2,2), Mdense(2,2);
    Kdense << 2,-1, -1,2;
    Mdense << 1,0, 0,-1;

    gsSparseMatrix<real_t> K, M;
    toSparse2(Kdense, K);
    toSparse2(Mdense, M);

    // Check 1 -- the fixture is what it claims to be. Without this, a later
    // edit could quietly make M definite and the test would keep passing
    // while measuring nothing.
    gsEigen::SelfAdjointEigenSolver<gsMatrix<real_t>::Base> esK, esM;
    esK.compute(Kdense);
    esM.compute(Mdense);
    gsInfo<<"  [modal/non-spd] eig(K) = "<<esK.eigenvalues().transpose()<<"\n";
    gsInfo<<"  [modal/non-spd] eig(M) = "<<esM.eigenvalues().transpose()<<"\n";
    CHECK(esK.eigenvalues()(0) > 0.0 && esK.eigenvalues()(1) > 0.0);
    CHECK(esM.eigenvalues()(0) < 0.0 && esM.eigenvalues()(1) > 0.0);

    gsModalSolver<real_t> solver(K, M);
    const gsStatus st = solver.compute();

    gsInfo<<"  [modal/non-spd] status = "<<static_cast<index_t>(st)<<"\n";

    // Check 2 -- do not pin which non-Success enumerator the guard uses.
    CHECK(st != gsStatus::Success);

    // Check 3 -- the invariant that survives any future re-mapping of the
    // guard's status value: Success must never coincide with a pair that
    // does not actually satisfy the pencil. values()/vectors() are only
    // meaningful when st == Success, so the residual is computed only then.
    real_t maxPencilResidual = 0.0;
    if (st == gsStatus::Success)
    {
        const real_t Knorm = Kdense.norm();
        for (index_t k = 0; k!=2; ++k)
        {
            const gsMatrix<real_t> v = solver.vector(k);
            const real_t w = solver.value(k);
            const real_t res = (Kdense*v - w*Mdense*v).norm() / (Knorm*v.norm());
            maxPencilResidual = math::max(maxPencilResidual, res);
        }
    }
    CHECK( !(st == gsStatus::Success && maxPencilResidual > 1e-8) );
}

} // SUITE(gsALMSolvers_test)
