/** @file example_ModifiedBratuExploration.cpp

    @brief MT benchmark driver: 1D MODIFIED-BRATU landscape traced with
    gsALMRiks + gsALMExploration, validated against EXACT discrete oracles.

    Reproduces the Wouters-thesis benchmark landscape (Sec. 9.4.1, fig. 9.3) --
    the modified Bratu equation from pde2path --

        -psi'' + lambda * (psi - mu * exp(psi)) = 0   on  Omega = [0,1],
        NEUMANN (homogeneous natural) BCs,  lambda = 10 FIXED,

    with the continuation parameter mu playing the role of gsALMBase's load "L".
    The thesis domain is [-0.5,0.5]; the problem is translation-invariant, so
    [0,1] is spectrally identical (same Neumann pencil eigenvalues) -- only the
    coordinate origin shifts, nothing measured here depends on it.

    == Closed-form structure (derive-checked; drives every acceptance oracle) ==
    Weak residual (test/trial space = the B-spline space, no elimination):
        R(psi,mu) = (Klap + lambda*M) * psi  -  mu*lambda * b(psi),
        b(psi)_i  = integral( exp(psi) * phi_i ),
        Me(psi)_ij= integral( exp(psi) * phi_i * phi_j ),   Jac = KplusLM - mu*lambda*Me.

    Curve A -- spatially CONSTANT states psi == c:
      * K*1 = 0 (partition of unity => the constant is in the stiffness null
        space) and exp(c) integrates EXACTLY on the spline space at Gauss p+1, so
        the constant c is an EXACT DISCRETE equilibrium with
              mu(c) = c * exp(-c).
      * Along A the reaction mass collapses (mu*exp(c) = c), so the tangent is
              J = Klap + lambda*(1-c)*M,
        whose eigenvalues in the pencil (K phi = theta M phi) are theta_k+lambda(1-c):
          - constant mode theta_0 = 0  -> zero at c = 1     => FOLD, mu* = 1/e,
          - mode k>=1  theta_k        -> zero at c = 1+theta_k/lambda => BRANCH pt.

    Fold of A at c = 1:   mu* = e^{-1} = 0.36787944117144233 (EXACT, discrete too).
      The singular mode is the CONSTANT => classified LIMIT point (folds are only
      approached by Riks, never traversed).

    A-B branch point (k=1):  c1 = 1 + theta_1/lambda ~ 1.987 (theta_1 ~ pi^2),
      mu_BP = c1*exp(-c1) ~ 0.2725  (thesis: "~0.27").  The mode is cos-like with
      ZERO MEAN => the reduced force projection is small => classified BRANCH point
      (gsALMBase's classification maps exactly onto this problem).

    Curve B -- nonconstant one-sided states bifurcating at mu_BP (a mirror pair;
      the explorer's +/- mode jobs give BOTH mirrors, which coincide in ||psi||).
      Thesis: B continues onward to curve C at mu ~ 0.03 (qualitative at our scope).

    == Round 5: the STATE-DEPENDENT forcing callback (opt-in, --forcingCallback) ==
    gsALMBase assumes a DEAD load: f = -dR/dLambda is the CONSTANT vector handed to the
    constructor. This problem violates that assumption exactly:
        R(psi,mu) = KplusLM*psi - mu*lambda*b(psi)   =>   -dR/dmu = lambda*b(psi),
    which DEPENDS on psi. The constructor argument Force = lambda*b(0) is therefore only a
    frozen representative; on a constant state psi == c the true load derivative is
    kappa*Force with kappa = exp(c). gsALMBase::setForcingFunction accepts the
    consistent f(psi,mu) = lambda*b(psi); it is byte-identical when unset and --forcingCallback
    is its first functional exercise anywhere in this plan.

    EXACT oracle for the callback itself (check 11, runs unconditionally, hard gate):
    on a constant state psi == c the partition of unity gives b(c*ones) = exp(c)*b(0)
    EXACTLY, hence  f(c*ones,mu) == exp(c)*Force  for EVERY mu. This is an identity, not an
    FD check. It exists because gsStructuralAnalysisOps<T>::ALForce_t and ALResidual_t are
    the SAME C++ type, so setForcingFunction(ALResidual) compiles silently and would
    produce garbage that nothing downstream could distinguish from physics.

    Gate (a) (check 12, only under --forcingCallback) asks the falsifiable question the
    callback exists for: does a plain Riks corrector ROUND the fold at mu* = 1/e with the
    consistent f, where the frozen-Force corrector cannot? The documented dead-load
    divergence criterion (see setForcingFunction's theory block)
        phi*u*kappa*(kappa-2) > 1-phi,   phi = 1/numDof,  u = ||K_T^{-1} Force||^2,
    is evaluated in CLOSED FORM along curve A first, as the quantitative prediction of
    where a frozen-Force corrector must break down; K_T = Klap + lambda*(1-c)*M is singular
    at c=1, so u blows up and the criterion must trip strictly BELOW mu*.
    Gate (b) is the full check list re-run with the callback installed: movement there is
    PREDICTED (the extended system's FD block is missing d(K_T V)/dLambda once the load
    stiffness makes K_T depend on Lambda), and is characterised, not chased.
    Callback OFF -- the default -- is byte-identical to round 4.

    Skeleton: examples/example_BratuExploration.cpp (the 2D Bratu-Gelfand fold
    driver); this is its 1D / Neumann / fixed-lambda modified-Bratu sibling.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gismo.h>

#include <gsStructuralAnalysis/src/gsStructuralAnalysisTools/gsStructuralAnalysisTypes.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMBase.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMLoadControl.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMRiks.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMCrisfield.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMExploration.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMLandscape.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    // ------------------------------------------------------------------ CLI
    index_t numHref    = 5;      // 32 elements on the default degree-2 basis
    index_t numElevate = 1;      // degree 1 + 1 = degree 2
    real_t  lambdaFixed= 10.0;   // FIXED reaction stiffness (thesis lambda)
    real_t  dMu        = 0.02;   // arc length in (psi,mu) space (mu-scale ~ 0.4,
                                 // so much finer than the Gelfand driver's 0.1)
    index_t maxPoints  = 20;     // MaxPointsPerCurve, a PER-SWEEP budget in the explorer:
                                 // the c=1.5 seed is swept in both arc-length directions
                                 // into ONE curve, so -N 20 allows up to 40 accepted
                                 // points on it (plus the stored singular point; MEASURED
                                 // at the default: 20 / 31 / 20 over the three curves).
                                 // 20 captures both fold sides (mu -> 1/e) and the A-B
                                 // branch traversal in ~3 s; the mu -> 0 tail of that
                                 // curve's backward LEG (the B-C region) is a slow
                                 // exp-blowup zone, opt-in via a larger -N.
    index_t maxCurves  = 8;      // MaxCurves
    // Branch-switch nudge scale: the explorer starts a branch job at U* + V/tau with
    // |V| = 1. It must be comparable to the amplitude curve B actually has one
    // SwitchLength away from the branch point, otherwise the corrector simply falls back
    // onto the parent (constant) branch A and the child curve is a retrace of A.
    // MEASURED with the branch-direction probe below (fixed-mu Newton from U*+eps*V at
    // mu = mu* - 0.005, printing the converged spread for a whole eps ladder):
    //   eps <= 0.15 -> spread ~ 1e-13 (back onto A);  eps = 0.4 -> Newton diverges;
    //   eps >= 0.8  -> spread 0.505   (curve B).
    // The pitchfork amplitude grows like sqrt(mu*-mu), so the requirement is a nudge of
    // the order of B's own amplitude. tau = 1 (nudge 1.0) sits inside that basin with
    // margin; the historic tau = 10 (nudge 0.1) is an order of magnitude too small -- it
    // was never exposed before because no correctly-directed branch job ever ran.
    real_t  tau        = 1.0;    // branch-switch nudge scale (Perturbation)
    // Root cause (see header note): the unstable branch's tangent
    // J = K - 5M is banded INDEFINITE, and SimplicialLDLT is
    // UNPIVOTED LDLT -> it silently mis-factorizes it, so the corrector diverges
    // from the EXACT c=1.5 seed. The pivoted "LU" (gsEigen::SparseLU) is the
    // default; determinant-stability HARD-REQUIRES SimplicialLDLT, so LU is
    // paired with the Eigenvalue stability method (gsSpectra), bifMethod=1.
    // gsALMLoadControl is the DEFAULT corrector. LoadControl fixes mu per
    // step and does a PLAIN Newton on U (exactly the check-2b corrector that
    // converges quadratically on the indefinite tangent J = K-5M). The multi-seed
    // bypass only needs to APPROACH the fold (LoadControl dies AT folds, which the
    // explorer handles by graceful truncation) -- it never rounds one.
    index_t method     = 0;        // 0=LoadControl (default), 1=Riks, 2=Crisfield
    std::string linSolver = "LU";  // "LU" (pivoted, default) or "SimplicialLDLT"
    index_t bifMethod  = 1;        // 0=Determinant (needs SimplicialLDLT), 1=Eigenvalue
    // SingularPointTestIt: inverse-power sweeps used to isolate the critical mode for the
    // limit-vs-branch test. The A-B crossing mode is the ZERO-MEAN cos mode; the classifier
    // needs enough sweeps to isolate it, otherwise |V.f|/|f| stays above the tolerance and
    // the crossing is misread as a LIMIT point (no branch job, check 5 FAILs).
    // MEASURED on the current pipeline: the A-B crossing classifies BRANCH at the library
    // default (20 sweeps), and even at --sptestit 5/7 the crossing reads BRANCH here (see
    // the table) -- the sweep count is not the lever that determines the classification.
    // These locals used to hard-code a LITERAL mirror of the
    // library defaults (5 and 1e-6), which silently went stale the moment the library default
    // changed and made every "default flags" run of this driver blind to a
    // library-default change. They now start at the sentinel -1 (index_t) / -1.0 (real_t) --
    // values outside each option's valid range ([1,.] resp. [0,1]) -- and `solver->options()`
    // below is set ONLY when `cmd.getValues` actually overwrote the local away from -1, i.e.
    // only when --sptestit / --sptesttol was passed explicitly. With no flag, the
    // freshly-constructed `solver` already carries the LIBRARY default from its own
    // `defaultOptions()` call, so that default governs unmodified. The `>= 0` sentinel test
    // (not `> 0`) matters because `SingularPointTestTol = 0` is itself a meaningful, valid
    // value elsewhere in this codebase (see gsALMBase.hpp's `_testSingularPoint` unresolved-
    // band guard) and must not be swallowed if a caller ever passes `--sptesttol 0`.
    index_t spTestIt  = -1;        // sentinel: unset => library SingularPointTestIt governs
    real_t  spTestTol = -1.0;      // sentinel: unset => library SingularPointTestTol governs
    real_t  switchLen = 0.005;     // arc length for branch (child) curves
    // Critical-mode SIGNS explored per singular point (explorer option "BranchPoints";
    // the LIBRARY default is 2, this driver hard-codes 1). See the long note
    // at the expl.options().setInt("BranchPoints",...) call for what each sign costs and
    // buys, and why 1 is NOT information-neutral at the A-C (even, mode-2) branch point.
    // The default landscape is unchanged at 1;
    // 2 is opt-in and fits inside the default --maxCurves 8 (MEASURED; the job
    // queue is FIFO and gated by MaxCurves, and the budget ladder is in the same note).
    index_t branchPoints = 1;      // BranchPoints (1 = + mode only; 2 = both signs)
    bool    plot       = false;
    // Verbose exploration output is ON by default (unchanged shipped behaviour; ctest runs
    // these drivers with NO arguments). gsCmdLine::addSwitch XORs the bound variable
    // (gsCmdLine.cpp:296-298), so a switch variable must start false -- binding the
    // already-true `verbose` to --verbose is what used to make the flag SILENCE the driver.
    bool    verbose    = true;   // resolved from the two switches below, after getValues
    bool    verboseFlag = false; // --verbose : explicit ON (wins over --quiet)
    bool    quietFlag   = false; // --quiet   : turn the default verbose output OFF
    // State-dependent forcing callback. OPT-IN; the default path must stay
    // bit-identical to the dead-load path (its CSV md5 is this driver's only bit-level
    // regression net).
    bool    useForcingCallback = false; // install f(psi,mu) = -dR/dmu = lambda*b(psi)
    index_t gateSteps  = 200;    // max ATTEMPTED steps in the gate-(a) fold-rounding loop
    // Diagnostic additions for measuring SingularPointComposite's never-converged
    // behaviour. All three
    // default to the shipped behaviour: solverVerbose/spComposite start false (gsCmdLine::
    // addSwitch XORs the bound bool, gsCmdLine.cpp:296-298), and spTolB's default is the
    // SAME 1e-4 literal it replaces below, so a bare run is unchanged.
    bool    solverVerbose = false; // print gsALMBase's per-iteration extended-solve residual table (Verbose)
    bool    spComposite   = false; // SingularPointComposite ON for the EXPLORATION-path extended solve
    real_t  spTolB     = 1e-4;     // SingularPointComputeTolB (bisection tol, 0=off)

    gsCmdLine cmd("MT benchmark: 1D modified-Bratu landscape with gsALMRiks + gsALMExploration.");
    cmd.addInt ("r","hRefine",        "Number of uniform h-refinement steps", numHref);
    cmd.addInt ("e","degreeElevation","Number of degree elevation steps", numElevate);
    cmd.addReal("l","lambda",         "Fixed reaction stiffness lambda", lambdaFixed);
    cmd.addReal("L","dMu",            "Arc length (mu-scale)", dMu);
    cmd.addInt ("N","maxPoints",      "Maximum number of accepted points per arc-length SWEEP "
                                      "(the explorer's MaxPointsPerCurve); a curve swept in both "
                                      "directions -- the c=1.5 seed here -- may hold twice this, "
                                      "plus its stored singular points", maxPoints);
    cmd.addInt ("","maxCurves",       "Maximum number of landscape curves", maxCurves);
    cmd.addReal("","tau",             "Branch-switch nudge scale (Perturbation)", tau);
    cmd.addInt ("m","method",         "Corrector: 0=LoadControl (default), 1=Riks, 2=Crisfield", method);
    cmd.addString("","solver",        "Linear solver: LU (pivoted, default) or SimplicialLDLT", linSolver);
    cmd.addInt ("","bifmethod",       "Stability method: 0=Determinant (needs SimplicialLDLT), 1=Eigenvalue", bifMethod);
    cmd.addInt ("","sptestit",        "SingularPointTestIt: inverse-power iterations for limit-vs-branch classification (unset/-1: library SingularPointTestIt default governs)", spTestIt);
    cmd.addReal("","sptesttol",       "SingularPointTestTol: |V.Force| threshold below which a crossing is a BRANCH (unset/negative: library SingularPointTestTol default governs)", spTestTol);
    cmd.addReal("","switchLen",       "Arc length for branch (child) curves", switchLen);
    cmd.addInt ("","branchPoints",    "BranchPoints: critical-mode signs explored per singular point (1 = + only, the driver default; 2 = both signs, i.e. 2 jobs per point, each swept in both arc-length directions -- needs --maxCurves >= 6 on this problem, which the default 8 already gives)", branchPoints);
    cmd.addSwitch("plot",   "Plot the branch-point solution field in ParaView format", plot);
    cmd.addSwitch("verbose","Verbose exploration output. It is ON by default; this flag "
                            "requests it explicitly and overrides --quiet", verboseFlag);
    cmd.addSwitch("quiet",  "Turn OFF the verbose exploration output that is on by default "
                            "(--verbose wins if both are given)", quietFlag);
    cmd.addSwitch("forcingCallback",
                  "Install the CONSISTENT state-dependent load derivative f(psi,mu) = lambda*b(psi) "
                  "via gsALMBase::setForcingFunction (default OFF = the dead-load path, byte-identical "
                  "to round 4). Enables gate (a), check 12.", useForcingCallback);
    cmd.addInt ("","gateSteps",  "Maximum attempted steps in the gate-(a) fold-rounding Riks loop", gateSteps);
    cmd.addSwitch("","solverVerbose","Print gsALMBase's per-iteration extended-solve residual table "
                                     "(sets the Verbose option; diagnostic switch)", solverVerbose);
    cmd.addSwitch("","spComposite","SingularPointComposite ON for the EXPLORATION-path extended solve "
                                    "(diagnostic switch; interacts with check (7)'s spCompositeAtExploration "
                                    "gate, see the comment above that capture)", spComposite);
    cmd.addReal("","spTolB","SingularPointComputeTolB (bisection tol, 0=off); default identical to the "
                            "literal it replaces (diagnostic default)", spTolB);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // --verbose wins over --quiet; neither flag => the shipped default (verbose ON).
    verbose = verboseFlag || !quietFlag;

    // Guard the library constraint: Determinant stability only works with
    // SimplicialLDLT (gsALMBase._computeStability dynamic_casts to it). If the
    // caller pairs LU with Determinant, promote to the Eigenvalue method.
    if (linSolver != "SimplicialLDLT" && bifMethod == 0)
    {
        gsWarn << "solver=" << linSolver << " is incompatible with bifmethod=0 "
                  "(Determinant needs SimplicialLDLT); switching to bifmethod=1 (Eigenvalue).\n";
        bifMethod = 1;
    }

    const real_t muFold = std::exp(-1.0); // 0.36787944117144233 -- EXACT discrete fold

    // Overall exit status: any hard [FAIL] flips this.
    bool allOk = true;
    auto report = [&allOk](bool ok, const std::string & msg)
    {
        gsInfo << (ok ? "[ OK ] " : "[FAIL] ") << msg << "\n";
        if (!ok) allOk = false;
    };
    // Soft check: reported as a FINDING but does NOT flip the process exit code.
    // Granted (by the task) ONLY to the extended on-PDE fold solve (check 4).
    auto softReport = [](bool ok, const std::string & msg)
    {
        gsInfo << (ok ? "[ OK ] " : "[FIND] ") << msg
               << (ok ? "" : "  (soft: documented finding, does not fail the run)") << "\n";
    };

    // ---------------------------------------------------- Discretisation setup
    gsStopwatch clock; clock.restart();

    // 1D unit interval geometry (degree 1), then elevate the basis to degree 1+e.
    gsMultiPatch<> mp;
    mp.addPatch( gsNurbsCreator<>::BSplineUnitInterval(1) );
    gsMultiBasis<> dbasis(mp, true);
    dbasis.setDegree( dbasis.maxCwiseDegree() + numElevate );
    for (index_t r = 0; r < numHref; ++r)
        dbasis.uniformRefine();

    gsInfo << "Patches: " << mp.nPatches()
           << ", degree: " << dbasis.minCwiseDegree()
           << ", elements: " << dbasis.totalElements() << "\n";

    // NEUMANN: homogeneous natural BCs = an EMPTY gsBoundaryConditions. No DoF is
    // eliminated, so the free-DoF vector IS the coefficient vector.
    gsBoundaryConditions<> bc;
    bc.setGeoMap(mp);

    gsExprAssembler<> A(1,1);
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    A.setIntegrationDomain(dbasis.domain());
    geometryMap G = A.getMap(mp);
    space u = A.getSpace(dbasis);
    u.setup(bc, dirichlet::l2Projection, 0); // no elimination (empty bc)

    gsMatrix<> solVector;
    solution u_sol = A.getSolution(u, solVector);

    A.initSystem();
    const index_t numDof = A.numDofs();
    gsInfo << "Number of free DoFs: " << numDof << "\n";

    // Laplacian stiffness Klap and mass Mmat, assembled ONCE (geometry is fixed).
    A.assemble( igrad(u,G) * igrad(u,G).tr() * meas(G) );
    gsSparseMatrix<> Klap = A.matrix();

    A.clearMatrix();
    A.assemble( u * u.tr() * meas(G) );
    gsSparseMatrix<> Mmat = A.matrix();

    gsSparseMatrix<> KplusLM = Klap + lambdaFixed * Mmat;
    KplusLM.makeCompressed();

    // Representative constant Force = lambda * b(0) = lambda * integral(phi_i).
    solVector.setZero(numDof);
    A.clearRhs();
    A.assemble( u * u_sol.val().exp() * meas(G) ); // exp(0) = 1
    gsVector<real_t> Force = A.rhs();
    Force *= lambdaFixed;                            // NON-CONST lvalue (ctor takes gsVector<T>&)
    const real_t ForceNorm = Force.norm();
    gsInfo << "ForceNorm = ||lambda*integral(phi_i)|| = " << ForceNorm << "\n";
    // Partition-of-unity sanity check: Force is ALREADY multiplied by lambdaFixed
    // above, so the basis-sum identity is Force.sum()/lambdaFixed == 1.0 (the unit-interval
    // length), NOT Force.sum() == 1.0. A mismatch would mean the basis count or the assembly is
    // not what the ForceNorm derivation below assumes. Measurement only -- this is intentionally
    // NOT a report()/softReport() clause (this check must not add an acceptance criterion).
    const real_t partitionOfUnity = Force.sum() / lambdaFixed;
    gsInfo << "Partition-of-unity check: Force.sum()/lambdaFixed = " << partitionOfUnity
           << " (deviation from 1.0 = " << (partitionOfUnity - (real_t)1) << ")\n";

    gsInfo << "Setup / one-off assembly: " << clock.stop() << " s\n";

    // -------------------------------------------------- theta_1 eigenvalue oracle
    // Dense generalized self-adjoint eigensolve of the Neumann pencil
    // (Klap, Mmat). theta_0 = 0 (constant mode); theta_1 = smallest eigenvalue
    // above 1e-8. Gives the semi-analytic branch-point prediction.
    real_t theta1 = 0.0, c1 = 0.0, muBP_pred = 0.0;
    gsVector<real_t> cosMode;   // discrete cos mode = theta_1 generalized eigenvector
    {
        gsStopwatch eigClock; eigClock.restart();
        // Use the plain Eigen matrix type: GeneralizedSelfAdjointEigenSolver's
        // internal evaluators do not specialize for gsMatrix.
        typedef gsEigen::Matrix<real_t,gsEigen::Dynamic,gsEigen::Dynamic> EigMat;
        EigMat Kd = Klap.toDense();
        EigMat Md = Mmat.toDense();
        gsEigen::GeneralizedSelfAdjointEigenSolver<EigMat> ges(Kd, Md);
        const gsVector<real_t> ev = ges.eigenvalues();
        for (index_t i = 0; i < ev.size(); ++i)
            if (ev[i] > 1e-8) { theta1 = ev[i]; cosMode = ges.eigenvectors().col(i); break; }
        c1        = 1.0 + theta1 / lambdaFixed;
        muBP_pred = c1 * std::exp(-c1);
        gsInfo << "theta_1 = " << std::setprecision(10) << theta1
               << " (pi^2 = " << EIGEN_PI*EIGEN_PI << ")"
               << "  ->  c1 = " << c1 << ",  mu_BP_pred = " << muBP_pred
               << "  (" << eigClock.stop() << " s)\n";
    }

    // ------------------------------------------------------------- ALM operators
    // mu side channel: the residual callback writes mu_current before the Jacobian
    // reads it (the ALM always evaluates the residual before the Jacobian here).
    real_t mu_current = 0.0;
    real_t asmTime = 0.0;
    gsStopwatch asmClock;

    // R(psi,mu) = KplusLM*psi - mu*lambda*b(psi),  b(psi)_i = integral(exp(psi) phi_i).
    gsStructuralAnalysisOps<real_t>::ALResidual_t ALResidual =
        [&](gsVector<real_t> const & x, real_t mu, gsVector<real_t> & result) -> bool
    {
        asmClock.restart();
        mu_current = mu;
        // Clamp psi for exp() only: legitimate states have psi <= ~6; a doomed
        // forward branch job (curve B lies at mu < mu*) runs psi off to infinity and
        // would overflow exp() to inf and segfault the sparse solve. Clamping keeps it
        // finite so the step fails gracefully. Untouched for all valid states.
        solVector = x.cwiseMin((real_t)30);
        A.clearRhs();
        A.assemble( u * u_sol.val().exp() * meas(G) ); // b(psi)
        result = KplusLM * x - mu * lambdaFixed * A.rhs();
        asmTime += asmClock.stop();
        return true;
    };

    // Jac(psi) = KplusLM - mu_current*lambda*Me(psi),  Me_ij = integral(exp(psi) phi_i phi_j).
    gsStructuralAnalysisOps<real_t>::Jacobian_t Jacobian =
        [&](gsVector<real_t> const & x, gsSparseMatrix<real_t> & m) -> bool
    {
        asmClock.restart();
        solVector = x.cwiseMin((real_t)30);   // exp overflow guard (see ALResidual)
        A.clearMatrix();
        A.assemble( u_sol.val().exp() * u * u.tr() * meas(G) ); // Me(psi)
        m = KplusLM - mu_current * lambdaFixed * A.matrix();
        m.makeCompressed();
        asmTime += asmClock.stop();
        return true;
    };

    // Round-5 CONSISTENT load derivative f(psi,mu) = -dR/dmu = lambda*b(psi).
    // b(psi) is the very vector ALResidual already assembles, so this reuses that single
    // expression (and its cwiseMin(30) exp-overflow guard) rather than opening a second
    // code path for b. It is deliberately mu-INDEPENDENT: R is affine in mu.
    // NOTE it must NOT touch the mu_current side channel -- that channel belongs to the
    // residual/Jacobian pair, and the solver calls the forcing at points in the sequence
    // where a stray write would silently mis-linearize the tangent (asserted in check 11).
    gsStructuralAnalysisOps<real_t>::ALForce_t ALForcing =
        [&](gsVector<real_t> const & x, real_t /*mu*/, gsVector<real_t> & result) -> bool
    {
        asmClock.restart();
        solVector = x.cwiseMin((real_t)30);   // exp overflow guard (see ALResidual)
        A.clearRhs();
        A.assemble( u * u_sol.val().exp() * meas(G) ); // b(psi)
        result = lambdaFixed * A.rhs();
        asmTime += asmClock.stop();
        return true;
    };

    // --------------------------------- (1) Mandatory FD self-consistency gate ----
    // Central-FD the residual callback column-wise against the Jacobian callback
    // at a nonzero state (x = 0.1*ones, mu = 0.5) -- guards a silent sign / mu
    // side-channel mistake.
    {
        gsStopwatch fdClock; fdClock.restart();
        const real_t eps = 1e-6;
        gsVector<real_t> x0 = gsVector<real_t>::Constant(numDof, 0.1);
        gsVector<real_t> R0(numDof);
        ALResidual(x0, 0.5, R0);          // sets mu_current = 0.5
        gsSparseMatrix<real_t> K0;
        Jacobian(x0, K0);                 // uses mu_current = 0.5
        gsMatrix<real_t> Kd = K0.toDense();

        real_t maxRelErr = 0.0;
        for (index_t j = 0; j < numDof; ++j)
        {
            gsVector<real_t> xp = x0, xm = x0;
            xp[j] += eps; xm[j] -= eps;
            gsVector<real_t> Rp(numDof), Rm(numDof);
            ALResidual(xp, 0.5, Rp);
            ALResidual(xm, 0.5, Rm);
            gsVector<real_t> col = (Rp - Rm) / (2.0*eps);
            const real_t denom = math::max( (real_t)1, Kd.col(j).norm() );
            maxRelErr = math::max(maxRelErr, (col - Kd.col(j)).norm() / denom);
        }
        gsInfo << "FD Jacobian max relative column error = " << maxRelErr
               << " (" << fdClock.stop() << " s)\n";
        report(maxRelErr <= 1e-4, "FD Jacobian consistency");
        if (!allOk) { gsInfo << "Aborting: FD gate failed.\n"; return EXIT_FAILURE; }
    }

    // ------------- (11) MANDATORY forcing-callback identity gate (round 5) -------
    // ANTI-FOOTGUN. gsStructuralAnalysisOps<T>::ALForce_t and ALResidual_t are the SAME
    // C++ type, so setForcingFunction(ALResidual) compiles silently and produces garbage.
    // This gate is the only thing between that mistake and every number gate (a)/(b)
    // reports, so it runs UNCONDITIONALLY (it only CALLS the lambda, it never installs it,
    // so the dead-load default path is untouched) and aborts exactly as gate (1) does.
    // EXACT oracle -- an identity, not an FD check: on a constant state psi == c the
    // partition of unity gives b(c*ones) = exp(c)*b(0), hence f(c*ones,mu) = exp(c)*Force
    // for every mu. Asserted at four c and two mu, plus mu-independence, plus the absence
    // of a mu_current side-channel write (the copy-paste failure mode that the returned
    // vector alone cannot catch).
    {
        gsStopwatch fcClock; fcClock.restart();
        real_t worstRel = 0.0, worstMuDep = 0.0;
        bool   bitwiseMuIndep = true, muChannelClean = true;
        gsVector<real_t> f1(numDof), f2(numDof);
        gsInfo << std::setprecision(4);
        for (real_t c : {(real_t)0.0,(real_t)0.5,(real_t)1.0,(real_t)1.5})
        {
            const gsVector<real_t> Uc = gsVector<real_t>::Constant(numDof, c);
            const real_t muBefore = mu_current;
            ALForcing(Uc, (real_t)0.3, f1);
            ALForcing(Uc, (real_t)1.7, f2);
            if (mu_current != muBefore) muChannelClean = false;
            const real_t rel = (f1 - std::exp(c)*Force).norm() / ForceNorm;
            const real_t dep = (f1 - f2).norm() / ForceNorm;
            if ((f1 - f2).cwiseAbs().maxCoeff() != (real_t)0) bitwiseMuIndep = false;
            worstRel   = math::max(worstRel,   rel);
            worstMuDep = math::max(worstMuDep, dep);
            gsInfo << "  forcing identity: c = " << c
                   << "  ||f(c*ones,mu) - exp(c)*Force||/||Force|| = " << rel
                   << ",  mu-dependence = " << dep << "\n";
        }
        gsInfo << "  worst relative error = " << worstRel
               << ", worst mu-dependence = " << worstMuDep
               << " (bitwise mu-independent: " << (bitwiseMuIndep ? "yes" : "no")
               << "; mu_current side channel untouched: " << (muChannelClean ? "yes" : "no")
               << "; " << fcClock.stop() << " s)\n";
        report(worstRel <= 1e-12 && worstMuDep <= 1e-14 && muChannelClean,
               "forcing callback matches exp(c)*Force on constant states");
        if (!allOk)
        {
            gsInfo << "Aborting: forcing-callback identity gate failed -- the callback is "
                      "mis-wired and nothing downstream would mean anything.\n";
            return EXIT_FAILURE;
        }
    }

    // ------------------------------- Exact constant-state seeds (curve A) --------
    // A constant psi == c has coefficient vector c*ones (partition of unity); with
    // no elimination this equals the free-DoF vector. mu(c) = c*exp(-c).
    const real_t c15   = 1.5;
    const real_t mu15  = c15 * std::exp(-c15);
    gsVector<real_t> U15 = gsVector<real_t>::Constant(numDof, c15);

    // (2) Seed sanity: the constant c=1.5 must be an EXACT discrete equilibrium.
    {
        gsVector<real_t> R(numDof);
        ALResidual(U15, mu15, R);
        gsInfo << "||R(U15, mu15)|| = " << std::setprecision(4) << R.norm()
               << " (mu15 = " << std::setprecision(10) << mu15 << ")\n";
        report(R.norm() <= 1e-10, "constant c=1.5 is an exact discrete equilibrium (||R|| <= 1e-10)");
    }

    // (2b) Corrector-isolation probe (round-2 root cause): a PLAIN fixed-mu Newton
    // solve from a PERTURBED c=1.5 state, using a direct pivoted LU factorization
    // of the INDEFINITE tangent J = K - 5M. If plain Newton converges quadratically
    // back to the constant state, the linear solver AND the Jacobian are sound on
    // the unstable branch -- so the arc-length TRACING divergence (curves 1-2 below)
    // is isolated to gsALMBase's arc-length CONSTRAINT update, not the algebra.
    {
        gsVector<real_t> Un = U15 + gsVector<real_t>::Constant(numDof, 1e-3); // perturb c=1.5
        gsVector<real_t> R(numDof);
        gsInfo << "[probe] plain fixed-mu Newton from c=1.5+1e-3 (direct LU on J=K-5M):\n";
        bool newtonOk = false;
        for (index_t it = 0; it < 20; ++it)
        {
            ALResidual(Un, mu15, R);               // sets mu_current = mu15
            const real_t rn = R.norm();
            if (it < 6 || rn < 1e-12)
                gsInfo << "         it " << it << "  ||R|| = " << std::setprecision(4) << rn << "\n";
            if (rn < 1e-12) { newtonOk = true; break; }
            gsSparseMatrix<real_t> J;
            Jacobian(Un, J);                        // J = KplusLM - mu_current*lambda*Me
            typename gsSparseSolver<real_t>::LU lu; // pivoted (COLAMD) -- exact on indefinite
            lu.compute(J);
            Un -= lu.solve(R);
        }
        const real_t backErr = (Un - U15).norm();
        gsInfo << "         converged=" << (newtonOk?"yes":"no")
               << ", ||U_newton - 1.5*ones|| = " << std::setprecision(4) << backErr << "\n";
        report(newtonOk && backErr < 1e-8,
               "plain fixed-mu Newton converges on the INDEFINITE unstable tangent "
               "(linear solver + Jacobian sound; arc-length divergence is a corrector issue)");
    }

    // ===================== (12) GATE (a): does Riks + the CONSISTENT forcing =====
    // ===================== ROUND the fold, where frozen Force cannot?        =====
    // The CORRECTOR path only -- this touches no extended solve and no explorer
    // heuristic: a direct gsALMRiks stepping loop from rest (psi = 0, mu = 0), run TWICE
    // with everything identical except setForcingFunction. Rounding (not merely
    // approaching) mu* = 1/e is the claim; the callback-off run is the causal control,
    // without which a "pass" would only say that this driver configuration rounds folds.
    if (useForcingCallback)
    {
        gsStopwatch gateClock; gateClock.restart();
        // gsALMRiks's constraint scaling for numDof > 1 -- the convex-weight branch of
        // gsALMRiks<T>::predictor() (member gsALMRiks::m_convexWeight).
        const real_t phi = 1.0 / (real_t)numDof;

        // -- PREDICTION FIRST (closed form, no solve of the ALM involved). On a constant
        // state c the true load derivative is kappa*Force, kappa = exp(c), while a frozen
        // -Force corrector uses Force; gsALMBase::setForcingFunction documents that Riks
        // then degenerates to a linear fixed-point iteration DIVERGING when
        // phi*u*kappa*(kappa-2) > 1-phi with u = ||K_T^{-1}Force||^2. K_T = Klap +
        // lambda*(1-c)*M is SPD for c<1 and SINGULAR at c=1, so u blows up as the fold is
        // approached and the criterion MUST trip strictly below mu* = 1/e. c_crit below is
        // where the frozen-Force run is predicted to break down.
        gsInfo << "\n[gate a] dead-load divergence criterion along curve A (phi = 1/numDof = "
               << std::setprecision(6) << phi << "):\n";
        gsInfo << "      c       mu(c)      kappa=e^c   u=||Kt^-1 f||^2   phi*u*k*(k-2)   1-phi    diverges?\n";
        real_t cCrit = 0.0; bool haveCrit = false;
        {
            real_t prevC = 0.0, prevGap = 0.0; bool first = true;
            for (real_t c : {(real_t)0.2,(real_t)0.4,(real_t)0.5,(real_t)0.6,(real_t)0.7,
                             (real_t)0.8,(real_t)0.9,(real_t)0.95,(real_t)0.99})
            {
                const gsVector<real_t> Uc = gsVector<real_t>::Constant(numDof, c);
                const real_t muc = c * std::exp(-c);
                gsVector<real_t> dummy(numDof);
                ALResidual(Uc, muc, dummy);                    // sets mu_current
                gsSparseMatrix<real_t> Kt; Jacobian(Uc, Kt);   // K_T = Klap + lambda(1-c)M
                typename gsSparseSolver<real_t>::LU lu; lu.compute(Kt);
                const gsVector<real_t> ut = lu.solve(Force);
                const real_t uu = ut.dot(ut), kap = std::exp(c);
                const real_t lhs = phi*uu*kap*(kap-2.0), rhs = 1.0 - phi;
                const real_t gap = lhs - rhs;
                gsInfo << "   " << std::setw(6) << c << "  " << std::setw(10) << muc
                       << "  " << std::setw(10) << kap << "  " << std::setw(14) << uu
                       << "  " << std::setw(14) << lhs << "  " << std::setw(8) << rhs
                       << "   " << (gap > 0 ? "YES" : "no") << "\n";
                if (!first && !haveCrit && prevGap <= 0 && gap > 0)
                { cCrit = prevC + (c-prevC)*(-prevGap)/(gap-prevGap); haveCrit = true; }
                prevC = c; prevGap = gap; first = false;
            }
        }
        if (haveCrit)
            gsInfo << "   => predicted frozen-Force breakdown at c_crit ~ " << cCrit
                   << ", i.e. mu ~ " << cCrit*std::exp(-cCrit)
                   << " (1/e = " << muFold << ")\n";
        else
            gsInfo << "   => criterion never trips on the sampled ladder\n";

        struct GateRun
        {
            index_t nAccepted, idxMax, nPostFold;
            real_t  muMax, cAtMax, minLen, worstRes;
            bool    underflow, stalled;
            std::vector<real_t>  mu, cEst, res;
            std::vector<index_t> its;
        };

        // One direct Riks trace from rest. Fresh solver object per run: no option, no
        // factorization and no secant is shared between the two arms of the contrast.
        auto runFromRest = [&](bool useCb, const char* tag) -> GateRun
        {
            GateRun R;
            R.nAccepted = 0; R.idxMax = -1; R.nPostFold = 0;
            R.muMax = -1e30; R.cAtMax = 0.0; R.minLen = dMu; R.worstRes = 0.0;
            R.underflow = false; R.stalled = false;

            mu_current = 0.0;                        // reset the residual/Jacobian side channel
            gsALMRiks<real_t> riks(Jacobian, ALResidual, Force);
            riks.options().setString("Solver",linSolver);
            riks.options().setInt   ("BifurcationMethod",1);   // LU pairs with Eigenvalue
            riks.options().setReal  ("Length",dMu);
            riks.options().setReal  ("Tol",1e-8);
            riks.options().setReal  ("TolF",1e-7);
            riks.options().setInt   ("MaxIter",30);
            riks.options().setSwitch("Verbose",false);
            riks.applyOptions();
            if (useCb) riks.setForcingFunction(ALForcing);
            riks.initialize();
            riks.setLength(dMu);                                  // BEFORE setPrevious
            riks.setSolution(gsVector<real_t>::Zero(numDof), (real_t)0);
            riks.setPrevious(gsVector<real_t>::Zero(numDof), (real_t)0);

            gsVector<real_t> Uold = gsVector<real_t>::Zero(numDof);
            real_t Lold = 0.0, dL = dMu;
            bool firstFailShown = false;
            const index_t postFoldCap = 15;  // gate (a) is a FOLD test: stop well before the
                                             // A-B branch point at c ~ 1.987
            for (index_t k = 0; k < gateSteps; ++k)
            {
                gsStatus st;
                try { st = riks.step(); } catch (...) { st = gsStatus::AssemblyError; }

                if (st != gsStatus::Success)
                {
                    if (!firstFailShown)
                    {
                        // A failed step leaves (m_U,m_L,m_Uprev,m_Lprev,arc length) untouched
                        // (iterationFinish() is only reached on convergence), so replaying it
                        // with Verbose ON reproduces it EXACTLY and prints the per-iteration
                        // Riks table -- the evidence the honest-BLOCKED clause asks for.
                        firstFailShown = true;
                        gsInfo << "  [" << tag << "] FIRST non-converged step from mu = "
                               << std::setprecision(10) << Lold << ", c_est = "
                               << std::setprecision(6) << (Uold.norm()/math::sqrt((real_t)numDof))
                               << ", arc length = " << dL
                               << " -- exact replay with Verbose ON:\n";
                        riks.options().setSwitch("Verbose",true);
                        riks.applyOptions(); riks.setLength(dL);
                        try { riks.step(); } catch (...) {}
                        riks.options().setSwitch("Verbose",false);
                        riks.applyOptions(); riks.setLength(dL);
                        riks.setSolution(Uold, Lold);
                    }
                    dL *= 0.5;
                    R.minLen = math::min(R.minLen, dL);
                    if (dL < 1e-8) { R.underflow = true; break; }
                    riks.setLength(dL);
                    riks.setSolution(Uold, Lold);
                    continue;                      // retry, no point consumed
                }

                const gsVector<real_t> Ucur = riks.solutionU();
                const real_t           Lcur = riks.solutionL();
                if (R.nAccepted > 0 &&
                    (Ucur-Uold).norm() <= 1e-12*math::max((real_t)1, Uold.norm()) &&
                    math::abs(Lcur-Lold) <= 1e-12*math::max((real_t)1, math::abs(Lold)))
                { R.stalled = true; break; }       // accepted point repeats its predecessor

                gsVector<real_t> Rv(numDof); ALResidual(Ucur, Lcur, Rv);
                R.mu  .push_back(Lcur);
                R.cEst.push_back(Ucur.norm()/math::sqrt((real_t)numDof));
                R.res .push_back(Rv.norm());
                R.its .push_back(riks.numIterations());
                R.worstRes = math::max(R.worstRes, Rv.norm());
                if (Lcur > R.muMax) { R.muMax = Lcur; R.idxMax = R.nAccepted; R.cAtMax = R.cEst.back(); }
                ++R.nAccepted;
                Uold = Ucur; Lold = Lcur;
                if (R.idxMax >= 0 && R.nAccepted-1-R.idxMax >= postFoldCap) break;
            }
            // Post-fold run: CONSECUTIVE accepted points past the mu-peak with mu DECREASING
            // and ||psi|| INCREASING -- i.e. the c > 1 upper branch, which is what "rounded"
            // means (as opposed to "approached and turned back").
            for (index_t i = R.idxMax+1; i >= 1 && i < R.nAccepted; ++i)
            {
                if (R.mu[i] < R.mu[i-1] && R.cEst[i] > R.cEst[i-1]) ++R.nPostFold;
                else break;
            }
            return R;
        };

        gsInfo << "\n[gate a] direct gsALMRiks trace of curve A from rest (Length = "
               << dMu << ", <= " << gateSteps << " steps, halve-and-retry on failure):\n";
        const GateRun on  = runFromRest(true , "callback ON ");
        const GateRun off = runFromRest(false, "callback OFF");

        auto printRun = [&](const char* tag, const GateRun & R)
        {
            gsInfo << "  [" << tag << "] accepted = " << R.nAccepted
                   << ", mu_max = " << std::setprecision(10) << R.muMax
                   << " (|mu_max-1/e| = " << math::abs(R.muMax-muFold) << ")"
                   << ", c_est at peak = " << std::setprecision(6) << R.cAtMax
                   << ", post-fold pts = " << R.nPostFold
                   << ", min arc length = " << R.minLen
                   << ", worst ||R|| = " << R.worstRes
                   << (R.underflow ? ", ARC-LENGTH UNDERFLOW" : "")
                   << (R.stalled   ? ", STALLED (repeated point)" : "") << "\n";
            gsInfo << "        pt          mu        c_est   it        ||R||\n";
            const index_t lo = math::max((index_t)0, R.idxMax-3);
            const index_t hi = math::min(R.nAccepted, R.idxMax+16);
            for (index_t i = lo; i < hi; ++i)
                gsInfo << "      " << std::setw(4) << i
                       << "  " << std::setprecision(10) << std::setw(14) << R.mu[i]
                       << "  " << std::setprecision(6)  << std::setw(11) << R.cEst[i]
                       << "  " << std::setw(3) << R.its[i]
                       << "  " << std::setw(12) << R.res[i]
                       << (i == R.idxMax ? "   <== mu peak" : "") << "\n";
        };
        printRun("callback ON ", on);
        printRun("callback OFF", off);

        // The clauses are reported SEPARATELY: (i) is a sampling-quality statement (the
        // trace must place a point within dc ~ 0.023 of c=1 to see 1/e to 1e-4), (ii)
        // "rounding, not approaching" is the actual claim, (iii) is the no-stall/no-
        // arc-length-underflow requirement -- a trace that only gets there by collapsing
        // the step is not a rounded fold. All three enter the pass boolean below.
        // muTol (clause (i)): the apex UNDERSHOOT of an arc-length trace that steps ACROSS
        // the fold, DERIVED rather than chosen. mu(c) = c*exp(-c) has mu'(1) = 0 and
        // mu''(1) = -exp(-1) = -muFold, so mu ~ muFold - (muFold/2)*(c-1)^2. gsALMRiks's
        // metric is phi*||dpsi||^2 + (1-phi)*dmu^2 with phi = 1/numDof; on the near-constant
        // states psi ~ c*ones we have ||dpsi|| = sqrt(numDof)*dc, hence phi*||dpsi||^2 = dc^2,
        // and at the apex dmu -> 0, so the achieved step IN c equals the arc length dMu. The
        // nearest accepted point therefore sits at most dMu/2 from the apex and
        //     |mu_max - 1/e|  <=  (muFold/2)*(dMu/2)^2  =  muFold*dMu^2/8.
        // At the default dMu = 0.02 that is 1.84e-5 -- consistent with the "dc ~ 0.023 to see
        // 1/e to 1e-4" statement above, which is the same relation read the other way round.
        //
        // THE BOUND SCALES AS dMu^2; THE CONSTANT 1e-4 DID NOT -- it froze the default step
        // size into the gate and trips at dMu = sqrt(8e*1e-4) = 0.0466. MEASURED on this
        // driver (OMP_NUM_THREADS=1, --forcingCallback --sptestit 7):
        //     -L 0.08 : ON dev 1.205e-4, post-fold pts 15, min arc length 0.08 (never
        //               halved), no stall -- the callback arm demonstrably ROUNDS the fold,
        //               the dead-load control does not (dev 0.02643, 0 post-fold points),
        //               and yet the gate printed "FAILED (the consistent forcing does NOT
        //               round the fold)". The bound at that step is 2.943e-4, so the
        //               measurement respects the derivation and only the frozen constant
        //               was wrong.
        //     -L 0.06 : ON dev 3.77e-5 against a bound of 1.655e-4 -- passes either way.
        // The 1e-4 floor is retained so the default path and every finer step are unchanged.
        // WIDENING IS SAFE: a false [ OK ] needs roundedOn AND !roundedOff and clause (i)
        // enters BOTH arms, so loosening it can only ADD to roundedOff -- turning a PASS
        // into a softReport'ed INCONCLUSIVE, never manufacturing attribution. Clause (ii)
        // keeps discriminating (the dead-load arm has 0 post-fold points at every dMu
        // measured, and its dev is two orders above the bound at -L 0.06 and 0.08).
        const real_t muTol = math::max( (real_t)1e-4, muFold*dMu*dMu/8 );
        const bool onClause1  = math::abs(on .muMax - muFold) <= muTol;
        const bool onClause2  = on .nPostFold >= 5;
        const bool onClause3  = !on.underflow && !on.stalled;
        const bool offClause1 = math::abs(off.muMax - muFold) <= muTol;
        const bool offClause2 = off.nPostFold >= 5;
        const bool roundedOn  = onClause1  && onClause2 && onClause3;
        const bool roundedOff = offClause1 && offClause2;
        gsInfo << "  clause (i)   |mu_max - 1/e| <= " << muTol << "  : ON = "
               << std::setprecision(4) << math::abs(on.muMax-muFold) << " (" << (onClause1?"yes":"no")
               << "), OFF = " << math::abs(off.muMax-muFold) << " (" << (offClause1?"yes":"no") << ")\n";
        gsInfo << "  clause (ii)  >= 5 post-fold points     : ON = " << on.nPostFold
               << " (" << (onClause2?"yes":"no") << "), OFF = " << off.nPostFold
               << " (" << (offClause2?"yes":"no") << ")\n";
        gsInfo << "  clause (iii) no stall, no arc-length underflow (ON arm only): "
               << (onClause3?"yes":"no") << " (min arc length = " << on.minLen
               << ", initial = " << dMu << ")\n";
        // THREE-WAY outcome, encoded in the boolean handed to the reporter -- NOT just in
        // the message. PASS requires the callback arm to round the fold AND the frozen
        // -Force control NOT to: if both round, the outcome is not attributable to the
        // callback and is INCONCLUSIVE, which must never be printed as [ OK ]. Dormant as
        // long as the control stalls, but a future session that changes a default must not
        // be told "OK" by a check whose premise has silently dissolved.
        const bool gatePass = roundedOn && !roundedOff;
        const bool gateInconclusive = roundedOn && roundedOff;
        const std::string verdict =
            gatePass          ? std::string("PASSED (fold rounding is attributable to the callback)")
          : (gateInconclusive ? std::string("INCONCLUSIVE (the callback-OFF control ALSO rounds the "
                                            "fold -- the outcome is NOT attributable to the callback)")
                              : std::string("FAILED (the consistent forcing does NOT round the fold)"));
        gsInfo << "  GATE (a) VERDICT: " << verdict << "\n";
        gsInfo << "  (gate a: " << gateClock.stop() << " s)\n";
        const std::string gateMsg = "gate (a): Riks + consistent state-dependent forcing "
                                    "ROUNDS the fold at 1/e -- VERDICT: " + verdict;
        // PASSED -> hard [ OK ]; FAILED -> hard [FAIL] (exit 1); INCONCLUSIVE -> [FIND], a
        // documented non-pass that does not fail the run (the task grants all three as
        // acceptable OUTCOMES, but only the first is a pass).
        if (gateInconclusive) softReport(false, gateMsg);
        else                  report(gatePass, gateMsg);
    }

    // ----------------------------------------------------------- Solver + explorer
    const char* methodName[3] = {"LoadControl","Riks","Crisfield"};
    GISMO_ENSURE(method>=0 && method<=2, "method must be 0 (LoadControl), 1 (Riks) or 2 (Crisfield)");
    gsInfo << "Corrector: " << methodName[method] << ", linear solver: " << linSolver
           << ", stability method: " << (bifMethod==0 ? "Determinant" : "Eigenvalue") << "\n";
    gsALMBase<real_t>* solver;
    if      (method==1) solver = new gsALMRiks<real_t>     (Jacobian, ALResidual, Force);
    else if (method==2) solver = new gsALMCrisfield<real_t>(Jacobian, ALResidual, Force);
    else                solver = new gsALMLoadControl<real_t>(Jacobian, ALResidual, Force);
    solver->options().setString("Solver",linSolver);
    solver->options().setInt   ("BifurcationMethod",bifMethod);   // 0=det, 1=eigenvalue
    solver->options().setReal  ("Length",dMu);
    solver->options().setReal  ("Tol",1e-8);
    // TolF is the corrector's RELATIVE force tolerance, ||R||/(|mu| ||f||). Check (7) below
    // accepts ||R|| <= 1e-6*max(1,|mu| ||f||) at every stored point, and here
    // |mu| ||f|| ~ 0.48 < 1, so that acceptance bound is ||R|| <= 1e-6, i.e. a relative
    // 2.1e-6 -- three orders TIGHTER than the library default TolF = 1e-3 the driver used
    // to configure. Asking a solver for 1e-3 and then asserting 1e-6 is an internally
    // inconsistent configuration, so TolF is set here to 1e-7 (an order of margin on
    // 2.1e-6) and check (7) now rests on a configured tolerance rather than on how far
    // past its tolerance Newton happens to overshoot.
    // CONTROL RUN: this is hygiene, not a result. With TolF left at 1e-3 the landscape,
    // the marked branch mu and every residual ratio below are IDENTICAL -- RE-MEASURED
    // against the shipped configuration at --sptestit 7: 3 curves / 71 points,
    // mu = 0.2724240767 (|dev| = 1.113e-05), certified worst ratio 1.460508833e-07, and
    // the landscape CSV comes out BYTE-IDENTICAL to the shipped run's, not merely equal
    // in these summary numbers. Plain Newton on this problem lands at ~1e-13 relative,
    // far past either threshold. The CSV md5 moving is therefore attributable to tau and
    // the branch-job direction alone.
    solver->options().setReal  ("TolF",1e-7);
    solver->options().setInt   ("MaxIter",30);
    solver->options().setReal  ("SingularPointComputeTolE",1e-8);
    solver->options().setReal  ("SingularPointComputeTolB",spTolB); // bisection ON
    // (C.1b unpinning): only override when the flag was actually given (local left
    // at its sentinel), so a bare invocation is governed by the library default set inside
    // solver's own constructor -- see the sentinel note at the spTestIt/spTestTol declaration.
    if (spTestTol >= (real_t)0) solver->options().setReal("SingularPointTestTol",spTestTol);
    if (spTestIt  >= 0)         solver->options().setInt ("SingularPointTestIt",spTestIt);  // round-4: isolate the cos mode
    solver->options().setReal  ("Perturbation",tau);
    if (method==2) { solver->options().setInt("AngleMethod",0); solver->options().setReal("Scaling",0.0); }
    // Diagnostic (SingularPointComposite never-converged measurement): Verbose defaults
    // false, reproducing the library default on a bare run (getOptions() reads Verbose at
    // gsALMBase.hpp:154, so a later set would be a no-op -- this is set BEFORE applyOptions()).
    // spComposite defaults false (guarded `if`, mirroring the sentinel pattern above), so a
    // bare run writes nothing here and the library default (off) governs. Check (4) below
    // re-applies options at its own call site but does not reset Verbose, so the switch stays
    // live through it; spComposite here governs only the EXPLORATION-path extended solve
    // (check (4)'s own probe hard-codes SingularPointComposite=true independently, see :1103).
    solver->options().setSwitch("Verbose",solverVerbose);
    if (spComposite) solver->options().setSwitch("SingularPointComposite",true);
    solver->applyOptions();
    // Check (7) below: capture the SingularPointComposite setting that governs the
    // exploration run started a few lines down (expl.solve(seeds)) BEFORE the later check-(4)
    // probe (~:1041) flips this switch to true on the SAME `solver` object for its own,
    // unrelated extended-solve demonstration. Reading options().getSwitch(...) from inside
    // check (7) itself would silently pick up the probe's leftover value instead of the value
    // that actually governed the stored landscape points -- see gismo:advisor's finding on
    // this task.
    const bool spCompositeAtExploration = solver->options().getSwitch("SingularPointComposite");
    // GATE (b): install the CONSISTENT load derivative on the explorer's solver. Movement
    // in the singular-point refinement below is PREDICTED, not a defect: with a callback
    // set, K_T depends on Lambda through the load stiffness, so _extendedSystemIteration's
    // FD block is missing d(K_T V)/dLambda (documented @warning on setForcingFunction).
    // The limit-vs-branch CLASSIFICATION is scale-invariant (it compares the cosine
    // |V.f|/||f||) and f = exp(c)*Force on curve A, so it is predicted NOT to move.
    if (useForcingCallback)
    {
        solver->setForcingFunction(ALForcing);
        gsInfo << "State-dependent forcing callback INSTALLED on the explorer's solver "
                  "(f = lambda*b(psi)); the extended system is INEXACT Newton in this mode.\n";
    }
    solver->initialize();

    std::string dirname = "ModifiedBratuResults";
    gsFileManager::mkdir(dirname);

    gsALMExploration<real_t> expl(solver);
    expl.options().setInt   ("MaxCurves",maxCurves);
    expl.options().setInt   ("MaxPointsPerCurve",maxPoints);
    expl.options().setReal  ("Length",dMu);
    expl.options().setReal  ("SwitchLength",switchLen);
    // C.7: RetraceBallFloor floors the FALLBACK-NUDGE branch-job's absolute
    // magnitude 1/Perturbation at margin*retraceThreshold()*max(1,||Ustar||) -- the
    // SAME normalization isRetrace's own ball uses -- so a nudge tuned at one (low
    // state-norm) branch point still clears the retrace ball at a HIGHER state-norm
    // one (this driver's A-C point). The library ships this option DISABLED (<= 0):
    // there is no problem-independent default, because the margin below is an
    // EQUALITY calibration at THIS driver's own A-B point with THIS driver's own
    // tau=1.0 (Perturbation), 1.0/(0.069282*11.5858) = 1.2458 -- i.e. the largest
    // margin that provably leaves the A-B nudge unmoved (nudgeBase already wins the
    // max there). A library default would either (a) be this same driver-tuned number
    // silently applied to every OTHER caller regardless of their Perturbation/state
    // scale (measured: at the library default Perturbation=1e3 this
    // floors the nudge ~1e2-1e4x above 1/Perturbation for an ordinary PDE state norm,
    // making Perturbation effectively inert everywhere with no opt-out), or (b) some
    // other number tuned on a different benchmark that has no more claim to
    // universality. Setting it HERE, next to this driver's own SwitchLength/Length/tau,
    // keeps the derivation and its justification together at the one locus it is
    // actually valid for; the same unpinning approach is the precedent for a driver
    // supplying its own calibrated value rather than leaning on an unjustified library
    // default.
    expl.options().setReal  ("RetraceBallFloor",1.2458);
    // BranchPoints counts critical-mode SIGNS, and each sign spawns exactly ONE job, which
    // gsALMExploration<T>::traceCurve sweeps in BOTH arc-length directions into a single
    // curve. So 1 here means ONE job at the A-B point, with two LEGS: mu down (where the
    // branch-direction probe below finds curve B) and mu up (which lands back on the
    // constant branch A).
    //   ** The mu-up leg is REWOUND, not stored. ** It is a SWEEP of the child curve now,
    // not a curve of its own, so the C_start test can discard it without costing the genuine
    // leg the other sweep found -- which is exactly why it can act at all. MEASURED on the
    // default run (--sptestit 7, OMP_NUM_THREADS=1, verbose), the measured baseline
    // (post the `stepsTaken >= startSteps` fix): the landscape is 3
    // curves / 72 points (0(20)/1(32)/2(20)) where the one-curve-per-direction code stored 5
    // curves / 92 -- the 20 removed rows ARE that leg.
    //   WHY it now fires, from the library's own contract. gsALMExploration<T>::isRetrace is
    // a point-to-STORED-POINT test against a DISCRETE sample of the other curves, not a
    // distance-to-locus test. The same fix changed the CALLER's evaluation cadence:
    // it is no longer evaluated once, at StartSteps, but at EVERY accepted step from StartSteps
    // onward against the same StartSteps-anchored threshold (the C_start block just above
    // does not stop testing after the first pass) -- see the RetraceTol paragraph of
    // gsALMExploration's class documentation and gsALMExploration<T>::retraceThreshold(),
    // which are the current statements of this contract. The knob it reads is RetraceTol (via
    // retraceThreshold()), NOT DedupTol: DedupTol is the branch-JOB dedup scope only. The
    // threshold is not the bare tolerance either -- retraceThreshold() scales it by
    // sqrt(StartSteps*SwitchLength/Length), which here (3, 0.005, 0.02) is 0.866, so the
    // EFFECTIVE threshold is 8e-2*0.866 = 6.93e-02.
    //   MEASURED (from landscape.h5, back when this leg was kept for its full budget): for
    // the mu-up leg, min over the OTHER curves' stored points of
    // max( ||U-U_p||/max(1,||U_p||), |L-L_p|/max(1,|L_p|) ) -- the isRetrace predicate itself
    // -- is 2.99e-02 at the StartSteps point, i.e. a factor 2.3 BELOW the 6.93e-02 threshold,
    // while the genuine curve-B leg of the same job sits at 1.61e-01, a factor 2.3 ABOVE it.
    // The two populations are separated with margin on both sides, and the default sits at
    // the geometric centre of that window by construction. (The child's points sit BETWEEN
    // the parent's, roughly half a parent grid step away, because SwitchLength = Length/4 --
    // which is why the tolerance has to be this loose and cannot be a job-comparison scope.)
    //   THE COUPLING THAT WOULD HAVE BLOCKED THIS IS GONE. Raising the tolerance from the
    // driver alone would be impossible if the explorer derived its no-progress guard as
    // moveTol = DedupTol*1e-2: at DedupTol = 3.5e-02 that guard would also reject the
    // near-stalled steps with which the upper-fold seed curve brackets mu = 1/e, and that
    // curve would lose its last TWO points (the rows mu = 0.367879 at ||U|| = 5.83463 and
    // 5.83350; at 7e-02 it would lose four). MoveTol is instead
    // its own option (default 1e-6 = the same derived value), the retrace test has its
    // own RetraceTol, and both fold-bracketing rows are present in the shipped landscape.
    //   WHY THE DEFAULT IS 1, AND EXACTLY HOW FAR THAT IS JUSTIFIED. At the A-B point the
    // critical mode is mode 1, which is ODD under the domain reflection, so curve B is a
    // MIRROR pair: the -mode leg is the reflection of the +mode leg and coincides with it
    // in ||psi||. There, one sign really is information-neutral and the second only costs
    // two more curves against MaxCurves plus an extra exp-blowup job.
    //   ** That argument does NOT extend to the A-C point. **
    // The A-C critical mode is mode 2, which is
    // EVEN under the reflection (64-report: parity ||v - rev(v)|| = 1.9e-14, MEASURED), so
    // the two sides of that pitchfork are NOT reflections of one another: one is a WELL
    // (centre 4.46, ends 7.57) and the other a BUMP (centre 7.57, ends 4.46), i.e. genuinely
    // different solution FIELDS. ** They do, however, agree in ||psi||: ** they are half-period
    // translates, psi_bump(x) = psi_well(x+1/2 mod 1) to 1.3e-04 on a field range of 3.11,
    // so their L2 norms coincide to 2e-07 relative (MEASURED). Only the PROFILE
    // separates them -- no norm can. With BranchPoints=1
    // only the + sign is ever nudged (gsALMExploration<T>::traceCurve loops
    // `for (s = 0; s < branchPoints; ++s)` with `sign = (s == 0) ? +1 : -1`, symbol anchor),
    // so exactly one side of the A-C pitchfork is reachable at the default. --branchPoints 2
    // exposes the other.
    //   ** BUDGET NOTE for --branchPoints 2. ** gsALMExploration<T>::solve consumes a FIFO
    // job queue under `while (!queue.empty() && nCurves() < maxCurves)` (symbol anchor), and
    // each sign queues ONE job -- swept in both arc-length directions -- i.e. 2 jobs per
    // singular point, in the order (+) (-). The merged seed curve crosses BOTH the A-B and
    // the A-C point in one trace, so the configuration that exposes A-C
    // (--sptesttol 1e-2 --tau 0.25 --branchPoints 2) queues 4 branch jobs in total and its
    // complete landscape is 6 curves / 132 points.
    //   MEASURED --maxCurves ladder for that configuration: 4 -> 4 curves / 92 points
    // (starved), 5 -> 5 / 112 (starved), 6 -> 6 / 132, and 7/8/10/12/20 all reproduce
    // 6 / 132 byte-for-byte. So --maxCurves >= 6 is the requirement and the DEFAULT 8
    // already meets it; a requirement of >= 16 would belong to the era of
    // 4 jobs per singular point and would only waste queue headroom now. (At the plain
    // default tolerances --branchPoints 2 finds only the A-B point and saturates at
    // 4 curves / 91 points, unchanged from --maxCurves 8 up to 20.)
    expl.options().setInt   ("BranchPoints",branchPoints);
    expl.options().setSwitch("Verbose",verbose);
    expl.options().setString("OutputPrefix",dirname + "/landscape");

    // Multi-seed: rest state traces A from below; the exact upper
    // constant c=1.5 traces A from above BOTH ways -- forward (c up: crosses the
    // A-B branch point) and backward (c down: descends toward the fold from above).
    std::vector<std::pair<gsVector<real_t>,real_t> > seeds;
    seeds.push_back( std::make_pair(gsVector<real_t>::Zero(numDof), 0.0) );
    seeds.push_back( std::make_pair(U15, mu15) );

    gsStopwatch exploreClock; exploreClock.restart();
    gsStatus estatus = expl.solve(seeds);
    gsInfo << "Exploration: " << exploreClock.stop() << " s\n";
    report(estatus == gsStatus::Success, "exploration returned Success");

    const gsALMLandscape<real_t> & ls = expl.landscape();
    gsInfo << "Landscape: " << ls.nCurves() << " curves, " << ls.nPoints() << " points.\n";
    if (ls.nCurves() == 0) { gsInfo << "No curve traced; aborting.\n"; delete solver; return EXIT_FAILURE; }

    // ---------------------------------------------------------- Curve identification
    // Identify curves by CONTENT (not index): the solver's forward direction sign
    // is not controlled, so the c=1.5 forward/backward legs can map to either
    // curve. mean(U) ~ c along a constant state.
    auto meanU = [](const gsVector<real_t> & U) -> real_t { return U.mean(); };
    auto curveMaxMu = [&](index_t c) -> real_t
    {
        const auto & pts = ls.curve(c).points;
        real_t m = pts.front().L;
        for (const auto & p : pts) m = math::max(m, p.L);
        return m;
    };
    auto curveMinMu = [&](index_t c) -> real_t
    {
        const auto & pts = ls.curve(c).points;
        real_t m = pts.front().L;
        for (const auto & p : pts) m = math::min(m, p.L);
        return m;
    };

    // Under LoadControl (mu-parametrised), the two c=1.5 seed legs are:
    //   * FORWARD  (mu up):   A upper toward the fold from above (c down -> 1),
    //     max mu -> 1/e  == the upper side of the two-sided fold bracket;
    //   * BACKWARD (mu down): A upper away from the fold (c up), THROUGH the A-B
    //     branch point (c~1.987, mu~0.2724) and on toward mu -> 0.
    // Discriminate by max mu: only the forward/upper-fold leg approaches 1/e.
    index_t idxRest = -1, idxUpperFold = -1, idxBackward = -1;
    {
        std::vector<index_t> seedCurves;
        for (index_t c = 0; c < (index_t)ls.nCurves(); ++c)
        {
            const gsALMLandscape<real_t>::Curve & cv = ls.curve(c);
            if (cv.parentCurve != -1) continue;          // seed curves only
            if (math::abs(meanU(cv.points.front().U)) < 0.25) idxRest = c; // starts at rest
            else seedCurves.push_back(c);
        }
        // Among the non-rest seed curves, the upper-fold leg has the larger max mu.
        for (index_t c : seedCurves)
        {
            if (idxUpperFold < 0 || curveMaxMu(c) > curveMaxMu(idxUpperFold))
            { idxBackward = idxUpperFold; idxUpperFold = c; }
            else idxBackward = c;
        }
        // ONE-CURVE-PER-SEED: the two arc-length directions are no longer two curve
        // objects. gsALMExploration<T>::traceCurve owns the curve and calls traceSweep
        // once per direction (symbol anchor), reverting the first sweep's block, so the
        // c=1.5 seed produces a SINGLE curve reading far(mu up) -> seed -> far(mu down).
        // The loop above then finds only ONE non-rest seed curve and leaves idxBackward
        // at -1, which silently skipped the inertia scan below and blanked two table rows.
        // The backward leg is now a point RANGE inside the merged curve, not a curve --
        // but BOTH consumers of idxBackward are whole-curve reductions that select that
        // range for free (see the justifications at each of them), so pointing idxBackward
        // at the merged curve is EXACT, not an approximation, and no per-half index
        // arithmetic is needed. This fires only when the loop above found no second seed
        // curve: at size 0 it assigns -1 to -1 (no-op, the -m 2 landscape), at size >= 2
        // idxBackward is already >= 0, so a genuine second seed curve is never overwritten.
        if (idxBackward < 0) idxBackward = idxUpperFold;
    }
    gsInfo << "Seed curves: rest=" << idxRest << ", upper-fold(mu up)=" << idxUpperFold
           << ", backward(mu down, through A-B)=" << idxBackward << "\n";

    // (3) Fold, two-sided: both the from-rest (lower) and the c-down (upper)
    // constant-state curves approach mu* = 1/e from BELOW (mu(c) < 1/e on both
    // sides of c=1). max of the two approached mu must sit just under 1/e and be
    // within 5e-3 of it.
    {
        real_t muLower = (idxRest >= 0)      ? curveMaxMu(idxRest)      : 0.0;
        real_t muUpper = (idxUpperFold >= 0) ? curveMaxMu(idxUpperFold) : 0.0;
        real_t muApproach = math::max(muLower, muUpper);
        gsInfo << "Fold approach: mu_max(lower/from-rest) = " << std::setprecision(10) << muLower
               << ", mu_max(upper/from-above) = " << muUpper << ", 1/e = " << muFold << "\n";

        // Parabolic-fit mu* from each side (report only), vertex of a 3-point
        // parabola bracketing the peak in cumulative-chord abscissa.
        auto parabolaVertex = [&](index_t c) -> real_t
        {
            const auto & pts = ls.curve(c).points;
            const index_t n = (index_t)pts.size();
            if (n < 3) return (n ? curveMaxMu(c) : 0.0);
            index_t kmax = 0; for (index_t p = 1; p < n; ++p) if (pts[p].L > pts[kmax].L) kmax = p;
            index_t i1 = kmax; if (i1 < 1) i1 = 1; if (i1 > n-2) i1 = n-2;
            const index_t i0 = i1-1, i2 = i1+1;
            const real_t invSqrtN = 1.0 / math::sqrt((real_t)numDof);
            auto chord = [&](index_t a, index_t b) -> real_t
            {
                const real_t du = ((pts[b].U - pts[a].U) * invSqrtN).norm();
                const real_t dl = pts[b].L - pts[a].L;
                return math::sqrt(du*du + dl*dl);
            };
            const real_t s0 = 0.0, s1 = chord(i0,i1), s2 = s1 + chord(i1,i2);
            gsMatrix<real_t> M(3,3); gsVector<real_t> rhs(3);
            M << s0*s0,s0,1.0, s1*s1,s1,1.0, s2*s2,s2,1.0;
            rhs << pts[i0].L, pts[i1].L, pts[i2].L;
            gsVector<real_t> abc = M.partialPivLu().solve(rhs);
            const real_t a = abc[0], b = abc[1], cc = abc[2];
            const real_t vertex = cc - b*b/(4.0*a);
            if (math::abs(a) > 1e-30 && vertex==vertex && math::abs(vertex-curveMaxMu(c)) < 1.0)
                return vertex;
            return curveMaxMu(c);
        };
        if (idxRest >= 0)
            gsInfo << "  parabolic mu*(lower) = " << parabolaVertex(idxRest)
                   << " (|dev| = " << math::abs(parabolaVertex(idxRest)-muFold) << ")\n";
        if (idxUpperFold >= 0)
            gsInfo << "  parabolic mu*(upper) = " << parabolaVertex(idxUpperFold)
                   << " (|dev| = " << math::abs(parabolaVertex(idxUpperFold)-muFold) << ")\n";

        report(idxRest >= 0 && idxUpperFold >= 0 &&
               muApproach <= muFold + 1e-8 && muApproach >= muFold - 5e-3,
               "two-sided fold approach: max(mu_lower,mu_upper) in [1/e-5e-3, 1/e+1e-8]");
    }

    // (4) Extended on-PDE fold solve, seeded EXACTLY at c=0.95 (a near-fold exact
    // equilibrium). bisection OFF, composite ON, m_V seeded via computeStability
    // +isBifurcation, then computeSingularPoint(...,false,false,false).
    {
        gsStopwatch extClock; extClock.restart();
        const real_t c095  = 0.95;
        const real_t mu095 = c095 * std::exp(-c095);
        gsVector<real_t> U095 = gsVector<real_t>::Constant(numDof, c095);

        solver->options().setReal  ("SingularPointComputeTolB",0);     // no bisection
        solver->options().setSwitch("SingularPointComposite",true);
        solver->options().setInt   ("MaxIter",100);
        solver->applyOptions();

        mu_current = mu095;                    // side channel for computeStability
        solver->setLength(dMu);
        solver->setSolution(U095, mu095);
        solver->setPrevious(U095, mu095);
        solver->computeStability(true);        // factorize J at (U095, mu095)
        solver->isBifurcation(false);          // power iteration seeds m_V

        gsStatus sp = gsStatus::NotConverged;
        try {
            sp = solver->computeSingularPoint(U095, mu095, /*switchBranch=*/false,
                                              /*jacobian=*/false, /*testPoint=*/false);
        } catch (...) { sp = gsStatus::OtherError; }

        bool extOk = false; real_t muExt = 0.0, uInfErr = 0.0;
        if (sp == gsStatus::Success)
        {
            muExt = solver->solutionL();
            gsVector<real_t> Uext = solver->solutionU();
            uInfErr = (Uext - gsVector<real_t>::Ones(numDof)).cwiseAbs().maxCoeff();
            extOk = true;
        }
        gsInfo << "extended fold solve status = "
               << (sp == gsStatus::Success ? "Success" : "NOT converged")
               << " -> mu = " << std::setprecision(10) << muExt
               << ", ||U-1||_inf = " << uInfErr << " (" << extClock.stop() << " s)\n";
        if (extOk)
            report(math::abs(muExt - muFold) <= 1e-6 && uInfErr <= 1e-6,
                   "extended on-PDE fold: |mu*-1/e|<=1e-6 and ||U*-ones||_inf<=1e-6");
        else
            softReport(false, "extended on-PDE fold solve did not converge (recipe followed verbatim)");
    }

    // (5-probe) Branch-direction probe: at the A-B branch point (U* = c1*ones,
    // mu* = mu_BP), fixed-mu Newton (LoadControl's corrector) seeded from U* + eps*V
    // (V = discrete cos mode) at mu around mu* -- report the mu-direction in which the
    // converged state is NON-CONSTANT (= curve B). This tells us whether the explorer's
    // forward (mu-up) branch job can reach curve B under LoadControl.
    {
        gsVector<real_t> V = cosMode; V.normalize();
        const gsVector<real_t> Ustar = gsVector<real_t>::Constant(numDof, c1);
        auto newtonSpread = [&](real_t muT, real_t eps) -> real_t
        {
            gsVector<real_t> Un = Ustar + eps * V;
            for (index_t it = 0; it < 40; ++it)
            {
                if (!Un.allFinite() || Un.maxCoeff() > 30) return 0.0; // exp overflow guard
                gsVector<real_t> R(numDof); ALResidual(Un, muT, R);
                if (!R.allFinite()) return 0.0;
                if (R.norm() < 1e-11) break;
                gsSparseMatrix<real_t> J; Jacobian(Un, J);
                typename gsSparseSolver<real_t>::LU lu; lu.compute(J);
                if (lu.info() != gsEigen::Success) return 0.0;
                Un -= lu.solve(R);
            }
            if (!Un.allFinite()) return 0.0;
            return Un.maxCoeff() - Un.minCoeff();   // 0 for a constant state
        };
        gsInfo << "Branch-direction probe (U*+eps*cosMode, fixed-mu Newton): max spread found per mu\n";
        for (real_t dmu : {(real_t)-0.08,(real_t)-0.05,(real_t)-0.02,(real_t)-0.005,
                           (real_t)0.005,(real_t)0.02,(real_t)0.05,(real_t)0.08})
        {
            real_t best = 0.0, bestEps = 0.0;
            for (real_t eps : {(real_t)-1.5,(real_t)-0.8,(real_t)-0.4,(real_t)-0.15,(real_t)-0.05,
                               (real_t)0.05,(real_t)0.15,(real_t)0.4,(real_t)0.8,(real_t)1.5})
            {
                const real_t s = newtonSpread(muBP_pred+dmu, eps);
                if (s > best) { best = s; bestEps = eps; }
            }
            gsInfo << "   mu = mu* + " << std::setprecision(3) << dmu
                   << "  max spread = " << std::setprecision(4) << best
                   << (best > 0.05 ? "  <== NON-CONSTANT (curve B), eps=" + std::to_string(bestEps) : "") << "\n";
        }
    }

    // (5-pre) Inertia scan of the backward leg: locate the A-B branch point by the
    // NEGATIVE-eigenvalue COUNT of the tangent going 1 -> 2 (the cos mode crossing at
    // c~1.987). This documents WHETHER the branch point was TRACED THROUGH, using an
    // INDEPENDENT computation -- the dense eigensolve below is this driver's own, not the
    // explorer's.
    //   IT IS A CORROBORATION, NOT A COMPLAINT. Round 3 recorded that the explorer's
    // detection fired on a min-eigenvalue SIGN flip, which cannot see a bifurcation on an
    // already-unstable branch (min eigenvalue = the constant mode, negative throughout;
    // only the COUNT changes). That is history: gsALMExploration<T>::traceSweep detects on
    // a change of the negative-eigenvalue count -- the same inertia signal this scan uses
    // -- and on the default run it DOES mark the A-B crossing (check (5) below reports the
    // marked mu). So the scan and the explorer agree here, and the printed lines must not
    // claim otherwise.
    //
    // Consumer 1 of idxBackward. Under one-curve-per-seed idxBackward is the MERGED
    // curve, and scanning it WHOLE is exact rather than a compromise: mu decreases
    // monotonically along the merged storage order (MEASURED on the default run, curve 1:
    // 0.367879 -> 0.354695 across the reverted mu-up block, then 0.314695 -> -0.0653048
    // across the mu-down block), so the merged curve IS the A-upper branch swept once from
    // the fold down to low mu. The tangent's negative count therefore goes 1 -> 2 exactly
    // ONCE along it, at the cos-mode crossing c ~ 1.987 -- which lies in the backward half,
    // so the loop breaks there and never enters the rest of the curve. That is why no
    // per-half point-range arithmetic is needed here, and it is checked by the value, not
    // assumed: a spurious crossing inside the mu-up block would report mu ~ 0.36 and the
    // softReport below would stay [FIND]. (The one point that is OUT of mu order is the
    // stored singular point, index 14 at mu 0.272424, inserted between 13 and 15; the break
    // fires at 13, before it.)
    real_t muInertiaChange = 0.0; bool tracedThroughBP = false;
    if (idxBackward >= 0)
    {
        // Return (#negatives, 2nd-smallest eigenvalue = the cos mode that crosses).
        auto probe = [&](const gsVector<real_t> & U, real_t mu, real_t & eig2nd) -> index_t
        {
            gsVector<real_t> dummy(numDof);
            ALResidual(U, mu, dummy);              // sets mu_current
            gsSparseMatrix<real_t> J; Jacobian(U, J);
            gsEigen::SelfAdjointEigenSolver<gsMatrix<real_t>> es(J.toDense(), gsEigen::EigenvaluesOnly);
            const gsVector<real_t> ev = es.eigenvalues(); // ascending
            eig2nd = ev[1];
            return (ev.array() < -1e-9).count();
        };
        const auto & pts = ls.curve(idxBackward).points;
        index_t prevNeg = -1; real_t muLo = 0, e2Lo = 0;
        for (const auto & p : pts)
        {
            real_t e2; const index_t neg = probe(p.U, p.L, e2);
            if (prevNeg == 1 && neg == 2)
            {
                // Linear interpolation of the crossing eigenvalue e2 (>0 at muLo, <0 here)
                // pins the branch mu (where det J = 0) between the two stored points.
                muInertiaChange = muLo + (p.L - muLo) * e2Lo / (e2Lo - e2);
                tracedThroughBP = true; break;
            }
            prevNeg = neg; muLo = p.L; e2Lo = e2;
        }
        // "backward LEG", not "backward curve": under one-curve-per-seed the id below is
        // the merged both-directions curve, and calling it the backward curve would assert
        // something false about it.
        //   BOTH lines are GUARDED on tracedThroughBP. Printing "tangent negatives 1->2"
        // and "TRACED THROUGH" unconditionally asserted the scan's conclusion whenever the
        // classifier had merely found a curve to scan: MEASURED at --sptestit 7 -N 3, where
        // the pre-guard code printed both of them next to "interpolated A-B branch mu = 0"
        // while the softReport below correctly read [FIND] -- three lines, two of which
        // contradicted the third.
        if (tracedThroughBP)
        {
            gsInfo << "Inertia scan (backward leg, curve " << idxBackward << "): tangent negatives 1->2; "
                      "interpolated A-B branch mu = " << std::setprecision(10) << muInertiaChange
                   << " (mu_BP_pred = " << muBP_pred << ", |dev| = "
                   << math::abs(muInertiaChange - muBP_pred) << ")\n";
            gsInfo << "  => the A-B branch point IS TRACED THROUGH on this curve. Detection in "
                      "gsALMExploration is inertia-based too (a change of the negative-eigenvalue "
                      "count), so this independent scan CORROBORATES it rather than substituting "
                      "for it; check (5) below reports what the explorer marked.\n";
        }
        else
            gsInfo << "Inertia scan (backward leg, curve " << idxBackward << "): NO 1->2 "
                      "tangent-inertia change along the stored points -- the A-B branch point "
                      "was not traced through on this curve (nothing is claimed about whether "
                      "it exists; the trace may simply not have reached it).\n";
    }
    // Soft record: the branch point was traversed at the predicted mu (physics correct),
    // even though the explorer could not mark it (a detection-criterion library gap).
    softReport(tracedThroughBP && math::abs(muInertiaChange - muBP_pred) <= 5e-3,
               "A-B branch point traced through at mu within 5e-3 of mu_BP_pred (inertia 1->2)");

    // (5) Branch point: a bifurcation-marked point near mu_BP_pred exists AND the
    // explorer classified it as a BRANCH point (child curves parented at it).
    {
        const std::vector<std::pair<index_t,index_t> > conn = ls.connectivity();
        bool markNearBP = false; real_t bestMu = 0.0, bestDev = 1e30;
        for (index_t c = 0; c < (index_t)ls.nCurves(); ++c)
        {
            const std::vector<index_t> bidx = ls.bifurcationIndices(c);
            for (index_t bi : bidx)
            {
                const real_t muB = ls.curve(c).points[bi].L;
                const real_t dev = math::abs(muB - muBP_pred);
                if (dev < bestDev) { bestDev = dev; bestMu = muB; }
                if (dev <= 5e-3) markNearBP = true;
            }
        }
        gsInfo << "Branch point: nearest marked mu = " << std::setprecision(10) << bestMu
               << " (|dev vs mu_BP_pred " << muBP_pred << "| = " << bestDev
               << "), child curves = " << conn.size() << "\n";
        report(markNearBP && !conn.empty(),
               "A-B branch point marked within 5e-3 of mu_BP_pred AND classified BRANCH (children exist)");
    }

    // (6) Curve B exists: a child curve (parent != -1) with >=5 points whose last
    // point is NON-constant (max_i U - min_i U >= 0.1).
    //
    // BOTH numbers were frozen constants compared against quantities that -N and
    // --switchLen set, and BOTH were reachable defects:
    //   * ">= 5 points": MaxPointsPerCurve caps children too, so at -N < 5 no branch job,
    //     however correct, could satisfy it. Fixed by capping the requirement at -N.
    //   * "spread >= 0.1": the last point's amplitude is the PITCHFORK amplitude, which
    //     grows as sqrt(mu* - mu), and the branch job travels nPts*switchLen in arc length
    //     from the branch point, so
    //           spread ~ C * sqrt(nPts * switchLen).
    //     MEASURED on this driver (OMP_NUM_THREADS=1, --sptestit 7), C = 7.2 .. 8.7 across
    //     3 orders in switchLen and a factor 5 in -N:
    //        switchLen 5e-3 (default) nPts 20 -> 2.5156      -N 12 -> 1.8523
    //        switchLen 1e-3           nPts 20 -> 1.0249      -N  8 -> 1.4792
    //        switchLen 2e-4           nPts 20 -> 0.4519      -N  6 -> 1.2677
    //        switchLen 5e-5           nPts 20 -> 0.2262      -N  4 -> 1.0249
    //        switchLen 1e-5           nPts 20 -> 0.1033
    //        switchLen 5e-6           nPts 20 -> 0.0749  <== FAILED the frozen 0.1
    //        switchLen 1e-5, -N 6             -> 0.0571  <== FAILED the frozen 0.1
    //        switchLen 1e-6           nPts 20 -> 0.0390  <== FAILED the frozen 0.1
    //     In each of those three the child curve IS curve B -- its spread is 6 to 12 orders
    //     above the OTHER population, the child that fell back onto the constant branch A
    //     (measured 3.6e-15 .. 9.6e-10 in the same runs) -- and yet check (6) hard-FAILed
    //     and the process exited 1.
    // The threshold therefore scales with the reach the branch job was GIVEN, per curve and
    // from the points actually traced (not from -N, which is only an upper bound):
    //       spreadTol = 0.1 * sqrt( nPts*switchLen / 0.1 ),  capped at the historical 0.1,
    // where the reference reach 0.1 = 20*0.005 is exactly the DEFAULT -N * --switchLen at
    // which the 0.1 was calibrated -- so the default is unchanged, by construction. The
    // discrimination survives: the ratio of a genuine B spread to spreadTol is
    // C*sqrt(0.1)/0.1 ~ 23 at EVERY reach (both sides scale as sqrt), while the fell-back
    // population stays 6+ orders below it. And it does not paper over real failures:
    // --tau 10 (the historic, too-small nudge) still yields two constant children
    // (spread 6.7e-15, 8.0e-15) and check (6) still hard-FAILs.
    //   FLOOR 1e-6, so the check can never become UNFAILABLE: at --switchLen 0 the scaled
    // tolerance would be 0 and "spread >= 0" would accept the constant-branch child this
    // check exists to reject. 1e-6 is DERIVED from the two measured populations, not chosen:
    // it sits ~3 orders above the LARGEST fell-back spread in the table above (9.64e-10,
    // which is corrector round-off on an O(1) state) and ~4.6 orders below the SMALLEST
    // genuine curve-B spread observed (0.0390). It can never mis-reject a real B leg:
    // C*sqrt(nPts*switchLen) drops below 1e-6 only for nPts*switchLen < 2e-14.
    // The two curve-B criteria are declared HERE, at function scope, because the ParaView
    // export at the end of this file must select on EXACTLY the same criteria that check (6)
    // accepts on. They were two independent copies of the same two literals, so the
    // scaling above left the export behind: at --switchLen <= 1e-5 --plot check
    // (6) passed while the literal 0.1 in the export block silently skipped every curve and
    // nothing was written (MEASURED before this change: --switchLen 5e-6 --plot -> check (6)
    // [ OK ] with spread 0.0749 >= 0.00316, 0 files exported). One definition, two uses.
    //
    // MaxPointsPerCurve (-N) caps EVERY curve, children included, so a fixed ">= 5 points"
    // is structurally impossible at -N < 5: no branch job, however correct, could ever
    // satisfy it. Cap the requirement by what -N permits.
    const index_t minChildPts = math::min((index_t)5, maxPoints);
    // Per-child non-constancy threshold; the derivation, the measured C = 7.15..8.72 law and
    // the two floors are documented in the block above.
    //   reachRef = 20*0.005 = the DEFAULT -N * --switchLen, the reach at which the historical
    //   0.1 was calibrated -- so the default is unchanged by construction.
    auto curveBSpreadTol = [&switchLen](size_t nPts) -> real_t
    {
        const real_t reachRef = 0.1;
        const real_t reachB   = (real_t)nPts*switchLen;
        return math::max( (real_t)1e-6,
               math::min( (real_t)0.1, (real_t)0.1*math::sqrt(reachB/reachRef) ) );
    };
    index_t nBranchCurves = 0, nBranchB = 0, nFellBack = 0; real_t muBmin = 1e30;
    std::string fellBackIds;
    {
        for (index_t c = 0; c < (index_t)ls.nCurves(); ++c)
        {
            const gsALMLandscape<real_t>::Curve & cv = ls.curve(c);
            if (cv.parentCurve == -1) continue;
            ++nBranchCurves;
            // points.back() below is unconditional, so an EMPTY child (a branch job whose
            // very first step failed) must be dropped here: it is counted as a branch curve
            // and rejected as curve B, which is what a 0-point child is.
            if (cv.points.empty())
            { gsInfo << "  child curve " << c << ": 0 points (empty branch job)\n"; continue; }
            const gsVector<real_t> & Ulast = cv.points.back().U;
            const real_t spread = Ulast.maxCoeff() - Ulast.minCoeff();
            const real_t spreadTol = curveBSpreadTol(cv.points.size());
            // Printed for every child, accepted or not: the accept/reject decision below
            // is a two-population discrimination and the populations must stay visible.
            gsInfo << "  child curve " << c << ": " << cv.points.size()
                   << " points, last-point spread = " << spread
                   << " (non-constant if >= " << spreadTol << ")\n";
            // (6a) A child whose last point is CONSTANT to within the SAME two-population
            // threshold that check (6) uses to accept curve B is, by that same
            // discrimination, a leg that fell back onto the constant branch A -- i.e. a
            // RETRACE of its parent that the C_start dedup did not remove. No new constant
            // is introduced: this is the SPREAD HALF of check (6)'s accept predicate,
            // inverted, read off the identical spread/spreadTol pair printed on the line
            // above. It is deliberately NOT gated on minChildPts as check (6) is -- a SHORT
            // constant child is still not a distinct solution locus -- and EMPTY children
            // (a branch job whose first step failed) never reach here at all, having been
            // skipped by the points.empty() branch above; they are neither a retrace nor a
            // curve-B leg.
            if (spread < spreadTol)
            {
                ++nFellBack;
                fellBackIds += (fellBackIds.empty() ? "" : ",") + std::to_string(c);
            }
            if ((index_t)cv.points.size() < minChildPts) continue;
            if (spread >= spreadTol)
            {
                ++nBranchB;
                for (const auto & p : cv.points) muBmin = math::min(muBmin, p.L);
            }
        }
        gsInfo << "Branch curves = " << nBranchCurves << ", non-constant curve-B legs = "
               << nBranchB << (nBranchB ? " (min mu on B = " + std::to_string(muBmin) + ")" : "") << "\n";
        report(nBranchB >= 1, "curve B exists: >=1 non-constant child curve with >="
               + std::to_string(minChildPts) + " points");

        // (6a) Surviving retraces -- PROMOTED to a HARD check (was softReport). The
        // DedupTol/moveTol coupling that made a fell-back child unremovable from this driver is
        // GONE (see the BranchPoints note above -- MoveTol and RetraceTol are separate options
        // now, and the C_start `>=` fix rewinds every fell-back leg at the default). A
        // survivor here now means gsALMExploration<T>::retraceThreshold() did not fire on that
        // leg -- either it genuinely left the parent branch far enough to clear the threshold,
        // or the run's Perturbation/RetraceTol put it outside the calibrated window -- and
        // either way it is a distinct-solution-locus claim the dedup exists to retire; with
        // that fix in place a survivor is a driver defect, not a configuration finding.
        //
        // One predicate, not two: "no surviving retrace" and "no child below its own spread
        // threshold" are the SAME predicate. nFellBack counts exactly the children with
        // spread < spreadTol -- the identical spread/spreadTol pair printed on the per-child
        // line above and used by check (6)'s own curve-B accept test just below. There is no
        // second, independently introduced threshold here.
        //
        // Gate on the NUMERATOR (nFellBack) only, never on the denominator (nBranchCurves).
        // "0 of 0" (a configuration that grows no child curve at all, e.g. -m 2, oracle f) must
        // print [ OK ] without asserting anything about a dedup that never ran there -- the
        // vacuity clause below exists for exactly that reading, and is why the gate does NOT
        // require nBranchCurves > 0. This is not a hole: the OTHER way "0 of 0" could arise --
        // the dedup destroying the genuine curve B along with the phantom -- is caught
        // independently by the existing HARD check just above
        // (report(nBranchB >= 1, "curve B exists: ...")), which fails on its own if that
        // happens. So nFellBack == 0 needs no denominator gate; adding one would put a
        // spurious extra [FAIL] into a landscape that legitimately has zero child curves.
        gsInfo << "Retrace children surviving the C_start dedup = " << nFellBack
               << " of " << nBranchCurves << " child curve(s)"
               << (nFellBack ? " (curve ids " + fellBackIds + "; each is a retrace of the "
                               "constant branch A, NOT a distinct solution locus)" : "") << "\n";
        std::string dedupMsg;
        if (nFellBack == 0)
        {
            // Positive polarity (requirement (b)): must not claim more than the clause tests
            // (requirement (d)) -- no claim that curve B survived, no claim the dedup ran
            // unless nBranchCurves > 0.
            dedupMsg = "no child curve survived as a retrace of the constant branch A (0 of "
                     + std::to_string(nBranchCurves) + " child curve(s)"
                     + (nBranchCurves == 0
                        ? ")"
                        : ", last-point spread >= its per-curve threshold on every child)");
            if (nBranchCurves == 0)
                // Requirement (c): kept verbatim -- this exact substring is what keeps oracle
                // f's "0 of 0" from reading as an assertion about a dedup that never ran.
                dedupMsg += " (vacuous: this configuration grew no child curve)";
        }
        else
        {
            dedupMsg = std::to_string(nFellBack) + " of " + std::to_string(nBranchCurves)
                     + " child curve(s) fell back onto the constant branch A and were NOT removed "
                       "by the C_start dedup (curve ids " + fellBackIds + "; a survivor means"
                       " gsALMExploration<T>::retraceThreshold() did not fire on that leg --"
                       " check the solver's Perturbation and the explorer's RetraceTol)";
        }
        report(nFellBack == 0, dedupMsg);

        // An earlier diagnosis found oracle c (--forcingCallback) BLOCKED --
        // basin selection at the branch-switch predictor state -- is a KNOWN, WRITTEN-UP hold, not
        // an unexplained failure. This marker is documentation only: it does NOT flip allOk and
        // is NOT routed through softReport (whose "does not fail the run" phrasing this retires).
        // The promoted report() call above is NOT conditionalized on
        // useForcingCallback -- per that earlier diagnosis, c's only child WAS the retrace
        // survivor, so once fixed it is removed, nFellBack == 0, and the promoted clause PASSES on
        // c already; what keeps c red is the pre-existing hard "curve B exists" check above,
        // a separate concern from this gate's.
        if (useForcingCallback)
        {
            gsInfo << "[HOLD] --forcingCallback (oracle c): basin selection at the branch-switch"
                      " predictor state (Euler predictor Ustart = job.U + s*escape*tauHatU,"
                      " gsALMExploration.hpp:621-622, escape = 1/Perturbation,"
                      " gsALMExploration.hpp:600-603 -- pinned to this driver's own"
                      " RetraceBallFloor calibration radius, 1.2458*retraceThreshold() = 0.08631)"
                      " is decided by corrector-level solver state that diagnosis did not identify.\n"
                      "       c's [FAIL] curve B exists is a KNOWN, WRITTEN-UP hold pending that"
                      " root cause, not an unexplained failure. The dedup-survivor check above is"
                      " unaffected, but note it passes VACUOUSLY on c (0 of 0 child curves: c"
                      " produces no children at all, so that check has nothing to discriminate"
                      " here and its pass is not evidence about c).\n";
        }
    }

    // (7) True residual within tolerance at every EQUILIBRIUM-CERTIFIED stored point.
    // ★ (orchestrator addition 2026-08-13): INVERTED --
    // corrected -- the meaning of Point::equilibrium: it is now true unconditionally for
    // every converged extended-system singular-point solve, whether or not that solve's
    // termination test included the equilibrium residuals (SingularPointComposite); it is
    // false only for a markUnresolvedSingular fallback row (see gsALMLandscape.h:64-95).
    // A per-row `!p.equilibrium` filter can therefore no longer tell "this row's ||R|| was
    // never checked" from "this row is a genuine converged corrector point that merely isn't
    // a certified bifurcation" -- exactly backwards for what this check needs. Per
    // gsALMLandscape.h's documented replacement filter (~:81-87) and 03b-report.md:305-325,
    // whether a refined singular-point row's ||R|| is actually bounded is a RUN-LEVEL
    // property of the SingularPointComposite option, not a per-row flag; a stored row IS a
    // refined singular-point row iff isBifurcation && stability==0 (the ONLY addPoint call
    // site that passes stability=0 is gsALMExploration.hpp's refined-singular-point branch;
    // every ordinary corrector/limit-point row gets stability = +-1, gsALMExploration.hpp's
    // indicator sign, so stability==0 cannot false-positive on those). `solver` never turns
    // SingularPointComposite on before `expl.solve(seeds)` runs (that option is opt-in
    // per-run: this driver enables it only for the UNRELATED check-(4) fold probe, further
    // below) -- `spCompositeAtExploration`, captured right after `solver->applyOptions()`
    // and before that probe mutates the switch, is exactly the value that governed the
    // stored landscape's singular-point solves.
    {
        gsStopwatch resClock; resClock.restart();
        real_t worst = 0.0, worstExcluded = 0.0; bool resOk = true;
        index_t nExcluded = 0;
        for (index_t c = 0; c < (index_t)ls.nCurves(); ++c)
        {
            index_t pi = 0;
            for (const auto & p : ls.curve(c).points)
            {
                gsVector<real_t> R(numDof);
                ALResidual(p.U, p.L, R);
                const real_t bound = 1e-6 * math::max( (real_t)1, math::abs(p.L) * ForceNorm );
                const real_t ratio = R.norm() / bound;
                const bool notCertified = p.isBifurcation && p.stability==0 && !spCompositeAtExploration;
                if (notCertified)
                {
                    ++nExcluded;
                    worstExcluded = math::max(worstExcluded, ratio);
                    // Diagnostic: |L|*ForceNorm and which regime max(1,|L|*ForceNorm) lands
                    // in, settling whether normalization (the clamp-to-1 regime) or a genuine
                    // difference in ||R|| explains the ratio gap between excluded points. This
                    // print is a measurement only -- it feeds no report()/softReport() clause.
                    const real_t lForceNorm = math::abs(p.L) * ForceNorm;
                    const real_t regime = math::max((real_t)1, lForceNorm);
                    gsInfo << "  [excluded, not equilibrium-certified] curve " << c
                           << " point " << pi << "  mu = " << std::setprecision(10) << p.L
                           << "  ||R||/bound = " << ratio
                           << "  |L|*ForceNorm = " << lForceNorm
                           << "  max(1,|L|*ForceNorm) = " << regime
                           << (lForceNorm < (real_t)1 ? "  (the `1` clamp wins)"
                                                       : "  (the |L|*ForceNorm term wins)")
                           << "  bound = " << bound << "\n";
                }
                else
                {
                    worst = math::max(worst, ratio);
                    if (R.norm() > bound) resOk = false;
                }
                ++pi;
            }
        }
        gsInfo << "true-residual worst ratio (||R||/bound) over certified points = " << worst
               << " (" << nExcluded << " point(s) excluded, worst excluded ratio "
               << worstExcluded << "; " << resClock.stop() << " s)\n";
        report(resOk, "true residual within tolerance at every equilibrium-certified stored point");
    }

    // (8) Thesis-comparison table (report / findings).
    gsInfo << "\n=================== Thesis (fig. 9.3) comparison ===================\n";
    gsInfo << std::setprecision(6);
    gsInfo << "  feature | computed        | semi-analytic / exact | thesis\n";
    gsInfo << "  --------+-----------------+-----------------------+-------\n";
    const real_t muUpperFold = (idxUpperFold>=0)?curveMaxMu(idxUpperFold):0.0;
    gsInfo << "  fold    | " << math::max((idxRest>=0)?curveMaxMu(idxRest):0.0, muUpperFold)
           << " (2-sided) | 1/e = " << muFold << "        | ~0.37\n";
    gsInfo << "  A-B     | " << (tracedThroughBP ? muInertiaChange : 0.0)
           << " (traced) | mu_BP_pred = " << muBP_pred << " | ~0.27\n";
    // Consumer 2 of idxBackward. A MIN over the merged curve is automatically the backward
    // leg's min: mu decreases monotonically along the merged storage order (see the inertia
    // scan), so the mu-up half contributes nothing below the seed's mu. No point-range
    // arithmetic needed; the value is the pre-merge one (MEASURED, -0.0653048).
    gsInfo << "  B-C     | min mu on backward = " << (idxBackward>=0 ? curveMinMu(idxBackward) : 0.0)
           << " | (qualitative)          | ~0.03\n";
    gsInfo << "====================================================================\n\n";

    // ------------------------------------------------------------------- Outputs
    ls.writeCsv(dirname + "/landscape.csv");
    gsInfo << "Landscape CSV written to " << dirname << "/landscape.csv\n";

    if (plot)
    {
        // Export the last non-constant curve-B solution field, if any.
        // SELECTION ONLY -- this loop feeds no report()/softReport() and must not start to:
        // it decides what is written, never whether the run passes. It selects on EXACTLY
        // the criteria check (6) accepts on (minChildPts and curveBSpreadTol, both declared
        // at check (6)), so the two can no longer disagree. With the literals they used to
        // duplicate, --switchLen <= 1e-5 --plot passed check (6) and exported NOTHING.
        for (index_t c = 0; c < (index_t)ls.nCurves(); ++c)
        {
            const gsALMLandscape<real_t>::Curve & cv = ls.curve(c);
            if (cv.parentCurve == -1 || cv.points.empty()) continue;
            if ((index_t)cv.points.size() < minChildPts) continue;
            const gsVector<real_t> & Ulast = cv.points.back().U;
            if (Ulast.maxCoeff() - Ulast.minCoeff() < curveBSpreadTol(cv.points.size())) continue;
            solVector = Ulast;
            gsExprEvaluator<> ev(A);
            gsParaviewCollection collection(dirname + "/branchBSolution", &ev);
            collection.newTimeStep(&mp);
            collection.addField(u_sol, "psi");
            collection.saveTimeStep();
            collection.save();
            gsInfo << "Curve-B solution written to " << dirname << "/branchBSolution.pvd\n";
            break;
        }
    }

    gsInfo << "\nTotal assembly time in callbacks: " << asmTime << " s\n";
    gsInfo << (allOk ? "\nAll acceptance checks passed.\n"
                     : "\nOne or more acceptance checks FAILED.\n");

    delete solver;
    return allOk ? EXIT_SUCCESS : EXIT_FAILURE;
}
