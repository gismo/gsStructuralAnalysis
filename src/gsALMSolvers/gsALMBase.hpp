/** @file gsALMBase.hpp

    @brief Performs the arc length method to solve a nonlinear equation system.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#pragma once

#include <typeinfo>
#include <limits>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMHelper.h>

namespace gismo
{

template <class T>
void gsALMBase<T>::defaultOptions()
{
    m_options.addInt ("MaxIter","Maximum iterations",100);
    m_options.addReal("Tol","Tolerance",1e-6);
    m_options.addReal("TolF","Tolerance",1e-3);
    m_options.addReal("TolU","Tolerance",1e-6);
    m_options.addReal("Perturbation","Set Perturbation factor Tau",1e3);

    m_options.addReal("Length","Arclength",1e-2);

    m_options.addSwitch("AdaptiveLength","Adaptive length",false);
    m_options.addInt ("AdaptiveIterations","Desired iterations for adaptive length",10);

    m_options.addSwitch("Quasi","Use Quasi Newton method",false);
    m_options.addInt ("QuasiIterations","Number of iterations for quasi newton method",-1);

    m_options.addInt ("BifurcationMethod","Bifurcation Identification based on: -1: nothing, 0: Determinant;  1: Eigenvalue",bifmethod::Eigenvalue);

    m_options.addInt ("SingularPointFailure","What to do when a singular point determination fails?: 0 = Proceed without singular point (keep the original solution); 1 = Apply the last (non-converged) solution anyway",SPfail::With);
    m_options.addReal("SingularPointTestTol", "Detection tolerance for singular points. The cosine |V.f|/|f| of the (normalised) mode shape and the forcing should be below this tolerance for a branch point (dimensionless, in [0,1])", 1e-4);
    m_options.addReal("SingularPointModeTol","Convergence tolerance for the critical-mode "
                      "inverse power iteration (direction change per sweep, measured up to "
                      "sign). <= 0 (default): derive it as SingularPointTestTol*1e-2.",-1);
    m_options.addInt ("SingularPointTestIt" , "Maximum number of inverse-power iterations for the critical-mode (singular point) test; the iteration also stops early once the mode changes by less than SingularPointModeTol",20);

    m_options.addReal("SingularPointComputeTolE", "Tolerance for the extended iterations to compute a bifurcation point", 1e-10);
    m_options.addReal("SingularPointComputeTolB", "Tolerance for the bisection iterations to compute a bifurcation point. If tol = 0, no bi-section method is used.", 0);
    m_options.addInt ("SingularPointBisIt", "Bisection stage's OWN probe budget (number of "
      "step() calls _bisectionSolve may spend localizing the crossing), independent of "
      "MaxIter (pde2path bifdetec's bisecmax). Bisection halves the bracket per probe, so "
      "the budget needed scales as ~log2(dLb/desired resolution); raise it for a larger "
      "dLb or a tighter SingularPointComputeTolB.", 40);
    m_options.addSwitch("SingularPointComposite","Extended-system solve also requires the equilibrium residuals (TolF/TolU) at termination",false);

    m_options.addReal("BranchTangentPerturbation",
      "Finite-difference perturbation 'del' of the swibra ABE emanating-tangent computation "
      "(computeBranchTangent). pde2path's own value and warning: 'don't choose del too small' "
      "-- the two K_T assemblies are differenced, so del below the assembler's noise floor "
      "returns noise.", 1e-3);
    m_options.addReal("BranchTangentTol",
      "RELATIVE floor below which an ABE quantity counts as zero in computeBranchTangent "
      "(alpha1bar against the magnitude of its own two terms; the load share of the incoming "
      "tangent; the biorthogonalisation cosine |phi1.psi1|/(|phi1||psi1|)). Default 1e-2 = 10x "
      "the O(del) relative resolution of the finite-differenced coefficients, so the verdict "
      "sits ABOVE the resolution of the quantity it judges (cf. SingularPointTestTol vs the "
      "critical-mode tolerance).", 1e-2);

    m_options.addString("Solver","Sparse linear solver", "SimplicialLDLT");

    m_options.addSwitch ("Verbose","Verbose output",false);

    m_options.addReal("Relaxation","Set Relaxation factor alpha",1.0);

}

template <class T>
void gsALMBase<T>::getOptions()
{
    m_maxIterations       = m_options.getInt ("MaxIter");
    m_tolerance           = m_options.getReal("Tol");
    m_toleranceF          = m_options.getReal("TolF");
    m_toleranceU          = m_options.getReal("TolU");

    m_tau                 = m_options.getReal("Perturbation");

    m_adaptiveLength      = m_options.getSwitch("AdaptiveLength");
    m_desiredIterations   = m_options.getInt ("AdaptiveIterations");

    m_quasiNewton         = m_options.getSwitch("Quasi");
    m_quasiNewtonInterval = m_options.getInt ("QuasiIterations");

    m_bifurcationMethod   = m_options.getInt ("BifurcationMethod");
    m_solver = gsSparseSolver<T>::get( m_options.getString("Solver") );
    if  (!dynamic_cast<typename gsSparseSolver<T>::SimplicialLDLT*>(m_solver.get()) && m_bifurcationMethod==bifmethod::Determinant)
    {
        gsWarn<<"Determinant method cannot be used with solvers other than LDLT. Bifurcation method will be set to 'Eigenvalue'.\n";
        m_bifurcationMethod = bifmethod::Eigenvalue;
    }

    m_verbose             = m_options.getSwitch ("Verbose");
    m_relax               = m_options.getReal("Relaxation");

    // Same rule as setLength(): m_arcLength_prev is the arc length that PRODUCED the
    // current secant, so it must survive a (re)read of the options once a step has been
    // accepted - otherwise an applyOptions() inside a stepping or failure-recovery loop
    // silently re-declares the secant normaliser and neuters the reduction. Before the
    // first accepted step there is no secant yet, so it is still seeded here (this is the
    // path every constructor takes, hence bit-identical initial behaviour).
    m_arcLength = m_arcLength_ori = m_options.getReal("Length");
    if (!m_stepTaken)
      m_arcLength_prev = m_arcLength;

    m_SPfail = m_options.getInt ("SingularPointFailure");
    m_SPTestTol = m_options.getReal("SingularPointTestTol");
    m_SPTestIt  = m_options.getInt ("SingularPointTestIt");

    // SingularPointModeTol: sentinel default (see defaultOptions doc). Resolved here,
    // immediately after m_SPTestTol/m_SPTestIt, so every classification downstream sees a
    // consistent pair. Do NOT clamp an explicit user value -- warn and obey.
    m_SPModeTol = m_options.getReal("SingularPointModeTol");
    if (m_SPModeTol <= (T)0)
      m_SPModeTol = (m_SPTestTol > (T)0) ? m_SPTestTol*(T)1e-2 : m_tolerance;
    else if (m_SPModeTol >= m_SPTestTol && m_SPTestTol > (T)0)
      gsWarn<<"SingularPointModeTol ("<<m_SPModeTol<<") is not below SingularPointTestTol ("
            <<m_SPTestTol<<"): the limit-vs-branch verdict is then evaluated at or below the "
              "resolution of the mode it thresholds.\n";

    m_SPCompTolE = m_options.getReal("SingularPointComputeTolE");
    m_SPCompTolB = m_options.getReal("SingularPointComputeTolB");
    m_SPBisIt    = m_options.getInt ("SingularPointBisIt");
    m_SPComposite = m_options.getSwitch("SingularPointComposite");

    m_branchTangentDel = m_options.getReal("BranchTangentPerturbation");
    m_branchTangentTol = m_options.getReal("BranchTangentTol");

}

template <class T>
void gsALMBase<T>::init(bool stability)
{
  if (stability)
  {
    this->_computeStability(m_U);
    m_stability = this->stability(); // requires Jabobian!!
  }
}

template <class T>
void gsALMBase<T>::computeLength()
{
  // m_desiredIterations := optimal number of iterations
  // Both counters are index_t: the ratio MUST be formed in T, otherwise integer division
  // truncates and the clamp below leaves only the values {0.5,1,2}, i.e. no proportional
  // response at all.
  // m_numIterations==0 would be an integer division by zero (UB). _step() itself never
  // reaches computeLength() with a zero count (its corrector loop starts at 1), but
  // resetLength() also calls it and may run before any step has been taken. A step that
  // "converged in zero iterations" is the strongest possible argument for a larger arc
  // length, so it maps onto the clamp's upper bound 2.
  T fac = (0 == m_numIterations) ? (T)2
                                 : (T)m_desiredIterations / (T)m_numIterations;
  if (fac < 0.5)
    fac = 0.5;
  else if (fac > 2.0)
    fac = 2.0;

  m_arcLength_prev = m_arcLength;
  m_arcLength = m_arcLength*fac;
}

template <class T>
T gsALMBase<T>::reduceLength(T fac)
{
  m_arcLength *= fac;
  return m_arcLength;
}

template <class T>
T gsALMBase<T>::resetLength()
{
  // Same rule as setLength()/getOptions(): restoring the ORIGINAL arc length says nothing
  // about the secant that is currently installed, so m_arcLength_prev is only (re)seeded
  // while no step has been accepted yet.
  // @note In the adaptive branch the guard is inert: computeLength() re-establishes
  //       m_arcLength_prev = m_arcLength itself, which is the accepted-step semantics
  //       introduced by the setLength() fix and is deliberately left alone here.
  m_arcLength = m_arcLength_ori;
  if (!m_stepTaken)
    m_arcLength_prev = m_arcLength;
  if (m_adaptiveLength)
    computeLength();
  return m_arcLength;
}

template <class T>
gsVector<T> gsALMBase<T>::computeResidual(const gsVector<T> & U, const T & L)
{
  gsVector<T> resVec;
  if (!m_residualFun(U, L, resVec))
    throw 2;
  return resVec;
}

template <class T>
void gsALMBase<T>::computeResidual()
{
  m_resVec = this->computeResidual(m_U + m_DeltaU, m_L + m_DeltaL);
}

template <class T>
void gsALMBase<T>::computeResidualNorms(bool extendedSolve)
{
  if (m_numIterations ==0 ) // then define the residual
  {
    // Residual scaling ||(L+DL)*f||: evaluated with the forcing at the current iterate
    // (m_U+m_DeltaU, m_L+m_DeltaL), i.e. the state the residual m_resVec itself was
    // evaluated at. For a dead load this is the stored constant, unchanged.
    const gsVector<T> & F = this->computeForcing();
    const T scaleF = ((m_L+m_DeltaL) * F).norm();

    // m_basisResidualF = m_resVec.norm();
    m_basisResidualF = scaleF;
    m_basisResidualU = m_DeltaU.norm();

    // The extended-system solve (_extendedSystemSolve) is the one caller whose
    // iteration-0 increment can already be at the answer: it is seeded from an
    // accurate localizeCrossing output, so ||DeltaU|| there measures round-off, not
    // a genuine step. A denominator frozen AT that round-off value makes the ratio
    // below O(1) forever, regardless of how converged the iterate already is --
    // TolU becomes unreachable rather than merely tight. Floor the basis at
    // floorRelU*||U+DeltaU|| (the state m_resVec/scaleF above were evaluated at, for
    // the same reason m_basisResidualF is evaluated there) so a genuinely-small
    // first increment cannot make the relative measure meaningless.
    //
    // floorRelU = 1e-1 is derived from the FULL per-iteration trajectory of every
    // measured extended-system crossing, NOT from the row-0 ratio alone, which only
    // says whether the floor changes m_basisResidualU, not whether the resulting bound
    // is ever small enough to reach. Because the frozen denominator applies to EVERY
    // iteration (not just iteration 0), what determines reachability is the SMALLEST
    // ||deltaU|| the trajectory visits at any row where the other two sub-tests (TolF,
    // the K_T.V test) already hold: floorRelU must exceed that row's
    // ||deltaU||/(TolU*||U||) for the solve to be able to stop there. The binding
    // requirement across the measured traps is 0.0417 (the shear crossing whose seed is
    // farthest from round-off); the ceiling is 0.4246 (the one case,
    // B6/riks_fold_extended_solve, whose basis already reflects a genuine,
    // non-degenerate first correction and must stay bit-identical -- measured against
    // the SAME ||U+DeltaU|| the code floors against here, 0.424606/1.0000). 1e-1 sits
    // in that (0.0417, 0.4246) window, 0.38 decades above the binding requirement and
    // 0.63 decades below the ceiling.
    // One consequence worth stating plainly: a row-0-only "binds/inert" criterion --
    // which only asks whether floorRelU exceeds the row-0 ratio -- systematically
    // UNDERSTATES the constant this needs, because later rows can visit a smaller
    // ||deltaU|| than row 0 even while the trajectory is not monotonically converging
    // (measured on the ModifiedBratu composite crossings: row 0 is not where the
    // minimum sits). This floor does not, and cannot, fix every single-condition-trap
    // case: where the extended corrector's raw increment never dips below the bound at
    // ANY iteration before the corrector destabilises (a genuine ill-conditioning of the
    // extended system near a critical tangent, not a normalisation artefact), no floor
    // consistent with the ceiling above will make it converge.
    //
    // Scoped to the extended solve only (extendedSolve==true, set by no caller but
    // _extendedSystemSolve): the ordinary corrector's own iteration-0 basis is the
    // PREDICTOR increment, a genuine step by construction, and must not be floored.
    if (extendedSolve)
    {
      const T floorRelU = static_cast<T>(1e-1);
      m_basisResidualU = math::max(m_basisResidualU, floorRelU*(m_U+m_DeltaU).norm());
    }

    // m_residueF = m_resVec.norm() / (((m_L+m_DeltaL) * m_forcing).norm());
    m_residueF = m_resVec.norm() / scaleF;
    // Uses m_basisResidualU rather than m_DeltaU.norm() directly, so that iteration 0
    // is measured against the SAME denominator as every later iteration (see the else
    // branch below). Whenever the floor is inert (extendedSolve==false, or
    // floorRelU*||U+DeltaU|| < ||DeltaU||) m_basisResidualU == m_DeltaU.norm() and the
    // ratio is 1 by construction, as it is for every ordinary-corrector step.
    m_residueU = m_deltaU.norm()/m_basisResidualU;
    m_residueL = m_DeltaL;
  }
  else
  {
    // m_residueF = m_resVec.norm() / (((m_L+m_DeltaL) * m_forcing).norm());
    m_residueF = m_resVec.norm() / m_basisResidualF;
    m_residueU = m_deltaU.norm() / m_basisResidualU;
    m_residueL = m_deltaL / m_DeltaL;
  }
  // gsInfo<<"m_DeltaU = "<<m_DeltaU.norm()
  //       <<"\t m_deltaU = "<<m_deltaU.norm()
  //       <<"\t m_deltaUbar = "<<m_deltaUbar.norm()
  //       <<"\t m_deltaUt = "<<m_deltaUt.norm()<<"\n";
  // gsInfo<<"m_resVec.norm() = "<<m_resVec.norm()
  //       <<"\t ((m_L+m_DeltaL) * m_forcing).norm() = "<<((m_L+m_DeltaL) * m_forcing).norm()
  //       <<"\t m_deltaU.norm() = "<<m_deltaUbar.norm()
  //       <<"\t m_basisResidualU = "<<m_basisResidualU
  //       <<"\t m_deltaL = "<<m_deltaL
  //       <<"\t m_DeltaL = "<<m_DeltaL<<"\n";

}

template <class T>
void gsALMBase<T>::factorizeMatrix(const gsSparseMatrix<T> & M)
{
  m_solver->compute(M);
  if (m_solver->info()!=gsEigen::ComputationInfo::Success)
  {
    gsInfo<<"Solver error with code "<<m_solver->info()<<". See Eigen documentation on ComputationInfo \n"
                                                                  <<gsEigen::ComputationInfo::Success<<": Success"<<"\n"
                                                                  <<gsEigen::ComputationInfo::NumericalIssue<<": NumericalIssue"<<"\n"
                                                                  <<gsEigen::ComputationInfo::NoConvergence<<": NoConvergence"<<"\n"
                                                                  <<gsEigen::ComputationInfo::InvalidInput<<": InvalidInput"<<"\n";
    throw 3;
  }
}

template <class T>
gsVector<T> gsALMBase<T>::solveSystem(const gsVector<T> & F)
{
  try
  {
    return m_solver->solve(F);
  }
  catch (...)
  {
    throw 3;
  }
}

template <class T>
gsSparseMatrix<T> gsALMBase<T>::_computeJacobian(const gsVector<T> & U, const gsVector<T> & deltaU)
{
  // Compute Jacobian
  gsSparseMatrix<T> m;
  m_note += "J";
  if (!m_djacobian(U,deltaU,m))
    throw 2;
  this->factorizeMatrix(m);
  return m;
}

template <class T>
gsSparseMatrix<T> gsALMBase<T>::computeJacobian(const gsVector<T> & U, const gsVector<T> & deltaU)
{
  return this->_computeJacobian(U,deltaU);
}

template <class T>
gsSparseMatrix<T> gsALMBase<T>::computeJacobian(const gsVector<T> & U)
{
  gsVector<T> DeltaU(U.rows());
  DeltaU.setZero();
  return this->computeJacobian(U,DeltaU);
}

template <class T>
gsSparseMatrix<T> gsALMBase<T>::computeJacobian()
{
  // Compute Jacobian
  if (m_deltaU.rows() == 0)
    m_deltaU = gsVector<T>::Zero(m_DeltaU.rows());
  return this->computeJacobian(m_U + m_DeltaU, m_deltaU);
}

template <class T>
void gsALMBase<T>::computeUbar()
{
  m_deltaUbar = this->solveSystem(-m_resVec);
}

template <class T>
void gsALMBase<T>::computeUt()
{
  // delta_u_t = K_T^{-1} f, with f = -dR/dL. The forcing is taken at the very state at
  // which the tangent currently held by m_solver was assembled (m_U+m_DeltaU,
  // m_L+m_DeltaL), so that delta_u_t is the exact bordered-Newton tangent also for
  // state-dependent loads. Without a callback f is the stored dead load (unchanged).
  m_deltaUt = this->solveSystem(this->computeForcing());
}

// template <class T>
// void gsALMBase<T>::initiateStep()
// {
//   m_converged = false;
//   m_numIterations = 0;
//   if (m_method == method::Riks)
//     initiateStepRiks();
//   else if (m_method == method::ConsistentCrisfield)
//     initiateStepConsistentCrisfield();
//   else if (m_method == method::ExplicitIterations)
//     initiateStepExplicitIterations();
//   else if (m_method == method::Crisfield)
//     initiateStepCrisfield();
//   else if (m_method == method::LoadControl)
//     initiateStepLC();
//   else
//   {
//     gsInfo<<"Error: Method unknown...\n Terminating process...\n";
//     std::terminate();
//   }
// }

// template <class T>
// void gsALMBase<T>::predictor()
// {
//   if (m_method == method::Riks)
//     predictorRiks();
//   else if (m_method == method::ConsistentCrisfield)
//     predictorConsistentCrisfield();
//   else if (m_method == method::ExplicitIterations)
//     predictorExplicitIterations();
//   else if (m_method == method::Crisfield)
//     predictorCrisfield();
//   else if (m_method == method::LoadControl)
//     predictorLC();
//   else
//   {
//     gsInfo<<"Error: Method unknown...\n Terminating process...\n";
//     std::terminate();
//   }
// }

// template <class T>
// void gsALMBase<T>::iterationFinish()
// {
//   if (m_method == method::Riks)
//     iterationFinishRiks();
//   else if (m_method == method::ConsistentCrisfield)
//     iterationFinishConsistentCrisfield();
//   else if (m_method == method::ExplicitIterations)
//     iterationFinishExplicitIterations();
//   else if (m_method == method::Crisfield)
//     iterationFinishCrisfield();
//   else if (m_method == method::LoadControl)
//     iterationFinishLC();
//   else
//   {
//     gsInfo<<"Error: Method unknown...\n Terminating process...\n";
//     std::terminate();
//   }
// }


// HV:
// to do: make a stand-alone (static? const?) Unew,Lnew = step(Uold,Lold) function, which can be used in the bisection functions without problems
// Some ideas: just make a (static) step function inside the classes and get rid of virtual sub-functions
//
template <class T>
gsStatus gsALMBase<T>::step()
{
  try
  {
    _step();
    m_status = gsStatus::Success;
  }
  catch (int errorCode)
  {
    m_status = _statusFromCode(errorCode);
  }
  catch (...)
  {
    m_status = gsStatus::OtherError;
  }
  return m_status;
}

template <class T>
void gsALMBase<T>::_step()
{
  GISMO_ASSERT(m_initialized,"Arc-Length Method is not initialized! Call initialize()");

  if (m_verbose)
    initOutput();

  m_converged = false;
  m_numIterations = 0;
  initiateStep();

  if (m_Uguess.rows()!=0 && m_Uguess.cols()!=0 && (m_Uguess-m_U).norm()!=0 && (m_Lguess-m_L)!=0)
    predictorGuess();
  else
    predictor();

  computeResidual();
  computeResidualNorms();

  if (m_verbose)
     stepOutput();

  if (m_quasiNewton)
  {
    quasiNewtonPredictor();
  }

  m_stabilityPrev = m_stability;
  // The counter is deliberately 1-BASED and the bound is INCLUSIVE, so the corrector runs
  // n = 1 .. m_maxIterations, i.e. exactly the MaxIter iterations the user asked for. The
  // bound `n < m_maxIterations` would run MaxIter-1 of them (ask for 30, get 29).
  // Re-basing the counter to 0 is NOT an option: computeLength() reads m_numIterations
  // as "the number of corrector iterations the accepted step took" and forms
  // m_desiredIterations/m_numIterations, so a re-base would shift the adaptive arc length of
  // every step that already converges. With `<=` the value of m_numIterations on the success
  // path is unchanged for every step that converges strictly before the limit.
  for (m_numIterations = 1; m_numIterations <= m_maxIterations; ++m_numIterations)
  {
    if ( (!m_quasiNewton) || ( ( m_quasiNewtonInterval>0 ) && ( m_numIterations % m_quasiNewtonInterval) < 1e-10 ) )
    {
      quasiNewtonIteration();
    }

    iteration();
    computeStability(false);

    computeResidual();
    computeResidualNorms();
    if (m_verbose)
       stepOutput();

    if ( m_residueF < m_toleranceF && m_residueU < m_toleranceU )
    {
      iterationFinish();
      // Change arc length. Both branches (re)establish the invariant that
      // m_arcLength_prev is the arc length that produced the secant (m_U-m_Uprev)
      // just written by iterationFinish(); m_stepTaken records that the invariant
      // now holds, so setLength() stops seeding m_arcLength_prev (see setLength()).
      if (m_adaptiveLength)
        computeLength();
      else
        m_arcLength_prev = m_arcLength;
      m_stepTaken = true;

      break;
    }
  }

  // The failure test lives OUTSIDE the loop and keys on m_converged, not on the counter.
  // An in-loop test would never fire when the body does not run, so at MaxIter <= 1 it
  // would fall out of _step() without throwing, and step() would then report Success with
  // m_converged == false and (U,L) uncorrected. Testing !m_converged instead keeps the
  // reported status and converged() in lockstep by construction, whatever the loop bound.
  if (!m_converged)
  {
    // Hygiene, not behaviour: the loop leaves the counter at m_maxIterations+1 on the
    // failure path, so numIterations() would report more iterations than were allowed.
    // Nothing in-tree reads it on a failed step (computeLength() runs on the success
    // path only), but the public accessor should not lie.
    m_numIterations = m_maxIterations;
    gsInfo<<"maximum iterations reached. Solution did not converge\n";
    throw 1;
  }
}

// ------------------------------------------------------------------------------------------------------------
// ---------------------------------------Singular point methods-----------------------------------------------
// ------------------------------------------------------------------------------------------------------------
template <class T>
void gsALMBase<T>::_computeSingularPoint(const gsVector<T> & U, const T & L, bool switchBranch, bool jacobian, bool testPoint)
{
  // THE REQUESTED STATE IS INSTALLED FIRST, BEFORE ANYTHING READS IT.
  // _testSingularPoint below evaluates the tangent (when jacobian==true, via
  // computeJacobian(m_U,m_deltaU)), the critical mode and the forcing
  // computeForcing(m_U,m_L) at the solver's CURRENT state; the same is true of the two
  // solve stages, which are handed m_U/m_L. Assigning after the test would silently
  // classify the point the solver happened to sit at rather than the point (U,L) the
  // caller asked about. Both non-trivial in-tree callers
  // (gsALMExploration, example_ShearExploration) pass the PRE-crossing point while the
  // solver sits at the POST-crossing one, so the two states genuinely differ.
  //
  // Consequence, deliberate: if _testSingularPoint throws (computeJacobian or a forcing
  // callback -> throw 2), the solver is left at (U,L) and not at its prior state. The
  // arguments are the requested state, so leaving the solver there is the consistent
  // contract; see the doxygen on computeSingularPoint.
  //
  // NOT changed here: m_deltaU is left alone, so computeJacobian(m_U,m_deltaU) still uses
  // whatever increment the last step left. That staleness is pre-existing (it was equally
  // stale before, only at a different point) and is out of scope for this fix.
  //
  // Captured BEFORE the assignments below overwrite it: the solver's own incumbent state is
  // the far bracket endpoint _bisectionSolve needs (see the comment block above this one --
  // the two states genuinely differ for the non-trivial callers).
  const gsVector<T> Uincumbent = m_U;
  const T           Lincumbent = m_L;

  m_U = U;
  m_L = L;

  // Localization stage: bisection, run before classification so that
  // _testSingularPoint below acts on the localized point and on the tangent this
  // stage last factorized, rather than on the caller's raw (U,L) and a possibly stale
  // tangent. Its return bool is not a reliable convergence signal (the trailing
  // eigenvector power-iteration loop can set it true independently of whether the
  // bisection loop converged), so ignore it.
  if (m_SPCompTolB != 0)
    this->_bisectionSolve(m_U,m_L,m_SPCompTolB,Uincumbent,Lincumbent);

  // Determines if the point is a bifurcation point. If not, it assumes it is
  bool test;
  if (testPoint)
  {
    test = this->_testSingularPoint(jacobian);
  }
  else
  {
    // No classification is performed on this path, so no verdict was formed either --
    // same staleness rule as _testSingularPoint's own entry (do not let a PREVIOUS
    // classification's verdict survive a call that skipped classification).
    m_SPverdict = SPverdict::NotTested;
    test = true;
  }

  if (test)
  {
    bool converged = false;
    // Refinement stage: extended system, run after classification on the point
    // localization left. Returns false when the iterations do not converge (it does
    // not throw and does not set m_status).
    converged = this->_extendedSystemSolve(m_U, m_L, m_SPCompTolE);
    // Publish the outcome of THIS solve. converged() must describe the most
    // recent solve the object performed, never a previous one. Placement is load-bearing:
    //  (1) AFTER the extended stage, so it overwrites whatever the nested step() calls
    //      inside _bisectionSolve left in m_converged -- those are live whenever
    //      SingularPointComputeTolB != 0 (default 0, but gsALMExploration_test's
    //      RiksFoldFixture sets 1e-6), and they would otherwise leave a stale `true`;
    //  (2) BEFORE the `if (!converged) throw 1;` below, so the failing path -- the one
    //      that matters -- reports correctly too;
    //  (3) switchBranch() below does not touch m_converged (verified), so the assignment
    //      survives the branch switch either way.
    m_converged = converged;

    if (switchBranch && (converged || m_SPfail==1))
      this->switchBranch();

    // here we assume that the stability of the singular point is
    // equal to the one of the previous point...
    // to avoid the algorithm to find a singular point again
    m_stability = m_stabilityPrev;

    // Propagate non-convergence of the singular-point solve to the caller via the
    // existing int-error-code convention (mapped to gsStatus::NotConverged in the
    // computeSingularPoint wrapper). Thrown *after* the switchBranch decision so the
    // "switch anyway on failure" behavior (SingularPointFailure==With) is unchanged;
    // only the reported status differs.
    if (!converged)
      throw 1;
  }
  else
  {
    // LIMIT POINT or UNRESOLVED MODE: the routine has NOTHING to deliver. The body above
    // is the entire singular-point computation, so falling through here without reporting
    // anything would return gsStatus::Success while only (U,L) had been copied into the
    // solver - i.e. "Success" and "nothing happened" would be indistinguishable, and a
    // caller that trusts the status (gsALMExploration) would store the UNREFINED input
    // point flagged as a converged bifurcation. Report it through the existing int-code
    // convention instead.
    //
    // Thrown AFTER m_U = U; m_L = L; on purpose: the ~13 in-tree callers that discard the
    // status keep exactly the state contract they have today.
    //
    // The code 1 (-> gsStatus::NotConverged) is shared with the non-convergence throw
    // above AND with the (new) Unresolved outcome. Since gsStatus is shared with
    // gsDynamicSolvers and gsAPALM, no enumerator is added for either cause; the MESSAGE
    // (unconditional, not gated on m_verbose) and singularPointVerdict() are the
    // discriminators -- see the caller contract on _testSingularPoint / singularPointVerdict.
    if (m_SPverdict == SPverdict::Unresolved)
      gsInfo<<"Singular point computation: possible singular point -- the critical mode "
              "could not be resolved (|dV| = "<<m_SPModeError<<" >= SingularPointTestTol = "
            <<m_SPTestTol<<"), so the point was NOT classified, NOT refined and no branch "
              "was followed. Raise SingularPointTestIt or loosen SingularPointModeTol.\n";
    else
      gsInfo<<"Singular point computation: the point classifies as a LIMIT point, which "
              "this routine cannot deliver; no extended solve was attempted and the "
              "solution was not refined by the extended system; it is the requested "
              "(U,L), localized first when the bisection stage is enabled "
              "(SingularPointComputeTolB != 0).\n";
    // PUBLISH the outcome of THIS call. Without this
    // assignment the branch returned gsStatus::NotConverged (throw 1, see below) while
    // converged() still read the flag of the PRECEDING ordinary step -- typically true --
    // i.e. exactly the status/flag contradiction this rule exists to remove, surviving on the one
    // path that was documented as an exception instead of being fixed. Nothing was solved
    // here, so `false` is the accurate report.
    // Placement is load-bearing and mirrors the `test` branch above: it sits AFTER
    // m_U = U; m_L = L; and BEFORE the throw, so the ~13 in-tree callers that discard the
    // status keep exactly the state contract they have today.
    m_converged = false;
    throw 1;
  }
}

// tolB and switchBranch will be defaulted
template <class T>
gsStatus gsALMBase<T>::computeSingularPoint(const gsVector<T> & U, const T & L, bool switchBranch, bool jacobian, bool testPoint)
{
  try
  {
    this->_computeSingularPoint(U,L,switchBranch,jacobian, testPoint);
    // Reached only when _computeSingularPoint did not throw, i.e. the singular-point
    // solve converged. Set Success explicitly here (rather than relying on a stale
    // m_status, which internal step() calls from the bisection stage may have left
    // non-Success) so a clean run reports Success and any thrown code overrides it.
    m_status = gsStatus::Success;
  }
  catch (int errorCode)
  {
    m_status = _statusFromCode(errorCode);
  }
  catch (...)
  {
    m_status = gsStatus::OtherError;
  }
  return m_status;
}

template <class T>
bool gsALMBase<T>::_inversePowerSweeps(gsVector<T> & V, index_t maxIt, T tol, T & dist)
{
  // Shared loop body of the inverse power iteration V <- K_T^{-1}V / ||K_T^{-1}V|| (or its
  // transpose analogue), on the ALREADY factorized operator held by m_solver. V must
  // already hold the (deterministic) start vector on entry.
  //
  // TERMINATION. maxIt is a MAXIMUM; the sweep stops as soon as the mode direction
  // settles. The distance is measured UP TO SIGN because the critical eigenvalue is
  // NEGATIVE on the far side of every crossing, where the iterate flips sign at every
  // sweep and ||V_k - V_{k-1}|| would never tend to zero.
  bool converged = false;
  dist = std::numeric_limits<T>::max();
  gsVector<T> Vprev;
  for (index_t k = 0; k != maxIt; k++)
  {
    Vprev = V;
    V = this->solveSystem(V);
    const T nrm = V.norm();
    if (!(nrm > (T)0) || !math::isfinite(nrm))
    {
      // The solve annihilated the iterate or blew up: keep the last usable direction
      // rather than handing a zero/NaN vector to switchBranch(), which normalises it.
      // dist is set to the WORST value (not left at the previous, optimistic sweep's
      // value): reaching this path means the iterate is unusable, and an unusable mode
      // must never be reported as merely imprecise.
      gsWarn<<"_inversePowerSweeps: the inverse iteration returned a vector of norm "<<nrm
            <<"; keeping the previous iterate.\n";
      V = Vprev;
      dist = std::numeric_limits<T>::max();
      break;
    }
    V /= nrm;
    dist = (V.dot(Vprev) < (T)0) ? (V+Vprev).norm() : (V-Vprev).norm();
    if (dist < tol)
    {
      converged = true;
      break;
    }
  }
  return converged;
}

template <class T>
bool gsALMBase<T>::_solverIsSelfAdjoint() const
{
  // SimplicialLDLT and SimplicialLLT are the SAME type (gsSparseSolver.h ~78-79), so one
  // dynamic_cast covers both. Self-adjoint families beyond that: CGDiagonal, CGIdentity
  // (Conjugate Gradient assumes SPD/self-adjoint by construction) and, where enabled,
  // PardisoLDLT/PardisoLLT. An unrecognised solver returns false: the transpose path in
  // _computeCriticalModeLeft is always correct, the free path only under symmetry.
  if (dynamic_cast<typename gsSparseSolver<T>::SimplicialLDLT*>(m_solver.get()))
    return true;
  if (dynamic_cast<typename gsSparseSolver<T>::CGDiagonal*>(m_solver.get()))
    return true;
  if (dynamic_cast<typename gsSparseSolver<T>::CGIdentity*>(m_solver.get()))
    return true;
#ifdef GISMO_WITH_PARDISO
  if (dynamic_cast<typename gsSparseSolver<T>::PardisoLDLT*>(m_solver.get()))
    return true;
  if (dynamic_cast<typename gsSparseSolver<T>::PardisoLLT*>(m_solver.get()))
    return true;
#endif
  // MINRES is self-adjoint by construction too, but its Eigen adaptor is not compiled in
  // this tree (gsSparseSolver.h ~212: the GISMO_EIGEN_SPARSE_SOLVER(gsEigenMINRES,...)
  // instantiation is commented out, so gsEigenMINRES<T> stays an INCOMPLETE type and a
  // dynamic_cast against it does not compile). Omitted here; the cost of the omission is
  // never wrong -- an unrecognised self-adjoint solver simply pays for a transpose.
  return false;
}

template <class T>
bool gsALMBase<T>::_computeCriticalMode(index_t maxIt, T tol)
{
  // Inverse power iteration V <- K_T^{-1}V / ||K_T^{-1}V||, which converges to the
  // eigenvector of the SMALLEST-MAGNITUDE eigenvalue of K_T, i.e. to the critical mode at
  // a (near-)singular point. Complexity: one factorization plus at most maxIt
  // back-substitutions, so O(maxIt) triangular solves — the assembly is NOT repeated.
  //
  // START VECTOR. A start vector of Ones -- a SYMMETRIC field -- is the obvious choice and
  // the wrong one: the canonical bifurcation of a symmetric structure has an ANTISYMMETRIC critical mode,
  // EXACTLY orthogonal to Ones (the zero-mean cos mode of a 1D reaction-diffusion
  // problem, the alternating Walsh mode of the unit test, an out-of-phase wrinkle): the
  // component the iteration must amplify is then zero — or, worse, pure rounding noise —
  // and NO number of sweeps recovers it. A GENERIC start has an O(1) component along
  // every mode with probability one, which is the textbook remedy.
  //
  // DETERMINISM. The start is drawn from a FIXED-SEED 32-bit LCG (Numerical Recipes
  // constants). It is pseudo-random but bit-reproducible: unsigned arithmetic is exactly
  // modular and the mapping to [-1/2,1/2) uses only powers of two, so the same start
  // vector is produced run to run, build to build and machine to machine. No global RNG
  // state is read or written, and nothing depends on the solver history — the returned
  // mode is a pure function of (K_T, maxIt, tol). This matters: the alternative
  // "warm-start from the previously converged mode" would make a CLASSIFICATION depend on
  // the path taken to reach the point.
  //
  // The sweep LOOP itself is shared with the left-mode iteration via _inversePowerSweeps
  // (see that function for the termination argument); this function only seeds m_V and
  // reports the outcome. The helper is called from two places whose tangent has a
  // different provenance (_testSingularPoint assembles it, _bisectionSolve inherits
  // whatever the last stability evaluation left in m_jacMat), so state the coupling it
  // relies on.
  GISMO_ASSERT(m_numDof == m_jacMat.cols(),
               "The tangent has "<<m_jacMat.cols()<<" columns but the solver has "
               <<m_numDof<<" degrees of freedom.");
  m_V.resize(m_numDof);
  unsigned seed = 2463534242u;                     // arbitrary but FIXED
  for (index_t i = 0; i != m_numDof; ++i)
  {
    seed = 1664525u*seed + 1013904223u;            // exact modulo 2^32
    // Top 24 bits -> a multiple of 2^-24 in [-1/2,1/2): exactly representable.
    m_V[i] = (T)( (index_t)(seed >> 8) ) / (T)16777216 - (T)0.5;
  }
  m_V.normalize();

  // The caller may have left another matrix in the solver (jacobian==false paths).
  this->factorizeMatrix(m_jacMat);

  T dist;
  bool converged = this->_inversePowerSweeps(m_V, maxIt, tol, dist);
  // Property of THIS mode iteration; both call sites (_testSingularPoint, the seeding call
  // inside _bisectionSolve) may update it, only _testSingularPoint turns it into a verdict.
  m_SPModeError = dist;

  // Warn only when the remaining direction error could actually MATTER: the mode enters
  // the classification through the cosine |psi.f|/|f|, so a direction error well below
  // m_SPTestTol can change neither the limit-vs-branch verdict nor, at that size, the
  // branch-switch direction. Warning there would make this message per-step noise instead
  // of a signal. Written !(dist < .) so that a non-finite dist also warns.
  if (!converged && !(dist < m_SPTestTol))
    gsWarn<<"_computeCriticalMode: the critical mode did not converge in "<<maxIt
          <<" iterations (SingularPointTestIt); |dV| = "<<dist<<" (target SingularPointModeTol = "<<tol
          <<") still exceeds SingularPointTestTol = "<<m_SPTestTol
          <<", so the limit-vs-branch classification and the branch-switch direction may be "
          <<"inaccurate. Increase SingularPointTestIt.\n";
  return converged;
}

template <class T>
bool gsALMBase<T>::_computeCriticalModeLeft(index_t maxIt, T tol)
{
  // Precondition: _computeCriticalMode has already run for THIS tangent (m_jacMat), i.e.
  // m_solver currently holds m_jacMat's factors and m_SPModeError holds the right mode's
  // direction error. The only caller, _testSingularPoint, calls the two back to back.
  if (this->_solverIsSelfAdjoint())
  {
    // K_T^T = K_T under the configured solver, so the left null vector IS the right one:
    // no second factorization, no extra triangular solve. m_SPModeError is left UNTOUCHED
    // -- the right-mode value already describes both vectors. The returned flag
    // reconstructs exactly what the preceding _computeCriticalMode call returned: both
    // calls are made with the SAME tol (see _testSingularPoint), so
    // m_SPModeError < tol <=> that call's "dist < tol".
    GISMO_UNUSED(maxIt);
    m_Vleft = m_V;
    return (m_SPModeError < tol);
  }

  // Non-symmetric (or unrecognised) solver: the transpose path is always correct. K_T is
  // an EXPRESSION TEMPLATE's .transpose(), not storage, so materialize it into a real
  // matrix before factorizing.
  gsSparseMatrix<T> Kt = m_jacMat.transpose();
  this->factorizeMatrix(Kt);

  // Same fixed-seed LCG start and determinism argument as the right mode (see
  // _computeCriticalMode); duplicated rather than shared because it seeds a DIFFERENT
  // vector (m_Vleft) sized from the same m_numDof.
  m_Vleft.resize(m_numDof);
  unsigned seed = 2463534242u;
  for (index_t i = 0; i != m_numDof; ++i)
  {
    seed = 1664525u*seed + 1013904223u;
    m_Vleft[i] = (T)( (index_t)(seed >> 8) ) / (T)16777216 - (T)0.5;
  }
  m_Vleft.normalize();

  T dist;
  bool converged = this->_inversePowerSweeps(m_Vleft, maxIt, tol, dist);
  m_SPModeError = math::max(m_SPModeError, dist);

  // RESTORE: callers of _testSingularPoint (_bisectionSolve, _extendedSystemSolve, the
  // correctors) assume m_solver holds m_jacMat's factors, not K_T^T's.
  this->factorizeMatrix(m_jacMat);

  if (!converged && !(dist < m_SPTestTol))
    gsWarn<<"_computeCriticalModeLeft: the left critical mode did not converge in "<<maxIt
          <<" iterations (SingularPointTestIt); |dV| = "<<dist<<" (target SingularPointModeTol = "<<tol
          <<") still exceeds SingularPointTestTol = "<<m_SPTestTol
          <<", so the limit-vs-branch classification may be inaccurate. Increase "
          <<"SingularPointTestIt.\n";
  return converged;
}

template <class T>
bool gsALMBase<T>::_testSingularPoint(bool jacobian)
{
  // Staleness rule (same precedent as m_converged): no verdict has
  // been formed yet on THIS call. Every return path below -- including the throw paths
  // reached from computeJacobian()/computeForcing() (throw 2) -- must either assign the
  // final verdict or leave this at NotTested; a throw leaves it here, at NotTested, on
  // purpose (no classification happened).
  m_SPverdict = SPverdict::NotTested;

  if (jacobian)
  {
    m_jacMat = this->computeJacobian(m_U,m_deltaU);
    this->factorizeMatrix(m_jacMat);
  }

  // Approximate the RIGHT critical mode phi of the tangent (deterministic generic start,
  // at most m_SPTestIt sweeps, resolved to m_SPModeTol -- decoupled from the verdict
  // threshold m_SPTestTol; see the SingularPointModeTol option), then the LEFT critical
  // mode psi (free on the default self-adjoint SimplicialLDLT path: see
  // _computeCriticalModeLeft). Both use the SAME tol so m_SPModeError below is consistent
  // across the two.
  this->_computeCriticalMode(m_SPTestIt, m_SPModeTol);
  this->_computeCriticalModeLeft(m_SPTestIt, m_SPModeTol);

  // C3(d) bands, taken from the pre-existing warning gate (historically
  // !converged && !(dist < m_SPTestTol)) and its reasoning: a direction error well below
  // the verdict threshold can change neither the limit-vs-branch verdict nor the
  // branch-switch direction, so only a mode that has NOT reached that resolution is a
  // problem.
  //   resolved   (m_SPModeError <  m_SPModeTol):                     classify normally
  //   imprecise  (m_SPModeTol <= m_SPModeError < m_SPTestTol):       classify normally,
  //                                                                  SILENTLY (see below)
  //   unresolved (m_SPTestTol > 0 && !(m_SPModeError < m_SPTestTol)): REFUSE
  // The m_SPTestTol > 0 guard is load-bearing, not cosmetic:
  // unittests/gsALMSolvers_test.cpp:1620 sets SingularPointTestTol = 0.0 deliberately, to
  // starve the branch verdict BY CONSTRUCTION (a cosine is never < 0); without the guard
  // !(d < 0) is always true and that configuration would refuse every point regardless of
  // mode quality. Triggering the refusal on !converged alone would be WRONG for the same
  // reason: with the tightened SingularPointModeTol a mode that lands at, say, 1e-7 would
  // be refused although its error cannot move a 1e-6 verdict.
  if (m_SPTestTol > (T)0 && !(m_SPModeError < m_SPTestTol))
  {
    m_SPverdict = SPverdict::Unresolved;
    if (m_verbose)
      gsInfo<<"\t Unresolved singular point (|dV| = "<<m_SPModeError
            <<" >= SingularPointTestTol = "<<m_SPTestTol<<")\n";
    return false;
  }
  // Band 2 (imprecise) is deliberately NOT warned about here, unconditionally: the
  // reasoning above exists precisely to suppress a message at this size, or it would be
  // per-step noise instead of a signal.

  // Classification of the singular point: the (approximate) LEFT null vector psi of K_T is
  // orthogonal to the load direction iff f is in the range of K_T, i.e. iff the point is
  // a branch point rather than a limit point. By the fundamental theorem of linear
  // algebra, range(K_T) is the orthogonal complement of ker(K_T^T) -- NOT of ker(K_T) --
  // so the vector that must be tested against f is the LEFT null vector psi, not the right
  // one phi (pde2path bifdetec.m:120-131 uses the same left-vector criterion, correct
  // under follower / state-dependent loads where K_T need not be symmetric). The two
  // coincide only when K_T is symmetric, which is exactly the concession
  // gsALMCrisfield.h's _chartSingularityProbe doxygen (~350-359) already makes for the
  // default SimplicialLDLT path; see _computeCriticalModeLeft / _solverIsSelfAdjoint for
  // when that free identity applies. The forcing is taken at the TESTED state (m_U,m_L) —
  // the same state the tangent above was evaluated at — which makes the test exact for
  // state-dependent loads as well.
  //
  // The comparison is made on the COSINE |psi.f|/||f|| (m_Vleft is unit-norm on both paths:
  // each sweep normalises, and the symmetric path copies the already-normalised m_V), NOT
  // on the raw product |psi.f|. The raw product carries the physical magnitude of the load
  // vector, so a fixed threshold on it really means the RELATIVE tolerance
  // m_SPTestTol/||f||: unreachable on any model with a physical load (||f|| >> 1), where
  // every branch point is then silently downgraded to a limit point, and far too permissive
  // on a normalised model (||f|| << 1), where even genuine limit points pass. The cosine is
  // dimensionless, so m_SPTestTol is a genuine angle tolerance: |cos| = O(1) at a limit
  // point (f not in range(K_T)) and |cos| -> 0 at a branch point (f in range(K_T)).
  //
  // Degenerate ||f|| = 0: the DIRECTION of f is undefined, so no value is "the limit" of
  // the cosine (which is scale-invariant and does NOT tend to 0 as f -> 0 along a fixed
  // direction). We adopt |psi.f|/||f|| := 0, i.e. BRANCH point, by convention: f = 0 lies in
  // range(K_T) trivially, so the branch-point criterion holds exactly. The choice is inert
  // in practice - with f = 0 the ALM predictor deltaUt = K_T^-1 f vanishes and the
  // continuation is degenerate at that point whichever way it is classified - hence the
  // warning rather than a hard error. No absolute floor is placed on ||f||: for a small
  // BUT NONZERO load (e.g. a Dirichlet-driven gsALMLoadControl problem) the cosine is
  // already scale-free, and a floor would be the only thing reintroducing dimensionality.
  //
  // Bind the forcing once: with a forcing callback (setForcingFunction) every call
  // overwrites m_forcingEval, which the returned reference aliases.
  const gsVector<T> & F = this->computeForcing(m_U,m_L);
  const T Fnorm = F.norm();
  // The test on the WARNING is !(Fnorm>0), not (Fnorm==0), so that a non-finite ||f|| -
  // reachable through a diverging corrector or a forcing callback evaluated at a
  // diverged state - is never silent either. Its CLASSIFICATION is deliberately left
  // alone: only ||f|| == 0 takes the convention below, so a NaN keeps propagating into
  // dot and fails the comparison, i.e. it still yields "limit point" exactly as it did
  // before this criterion was normalised.
  if (!(Fnorm > (T)0))
    gsWarn<<"_testSingularPoint: the forcing at the tested state has |f| = "<<Fnorm
          <<" (vanishing or not finite); the limit-vs-branch test is degenerate.\n";
  const T dot = (Fnorm == (T)0) ? (T)0 : math::abs(m_Vleft.dot(F)) / Fnorm;
  if ( (dot / m_SPTestTol > 1e-1) && (dot / m_SPTestTol < 10) )
  {
    gsInfo<<"Warning: the singular point test is close to its tolerance. |V.f|/(|f| tol) = "<<dot/m_SPTestTol<<" |V.f|/|f| = "<<dot<<"\t |f| = "<<Fnorm<<"\t tolerance = "<<m_SPTestTol<<"\n";
  }
  if (dot < m_SPTestTol)    // Bifurcation point
  {
    m_SPverdict = SPverdict::Branch;
    if (m_verbose) {gsInfo<<"\t Bifurcation point\n";}
    return true;
  }
  else        // Limit point
  {
    m_SPverdict = SPverdict::Limit;
    if (m_verbose) {gsInfo<<"\t Limit point\n";}
    return false;
  }
}

template <class T>
bool gsALMBase<T>::isBifurcation(bool jacobian)
{
  // Controls the singular point test with the first two arguments
  return this->_testSingularPoint(jacobian);
}

template <class T>
gsStatus gsALMBase<T>::computeStability(bool jacobian, T shift)
{
  try
  {
    _computeStability(m_U,jacobian,shift);
    m_status = gsStatus::Success;
  }
  catch (int errorCode)
  {
    m_status = _statusFromCode(errorCode);
  }
  catch (...)
  {
    m_status = gsStatus::OtherError;
  }
  return m_status;
}

template <class T>
void gsALMBase<T>::_computeStability(const gsVector<T> & x, bool jacobian, T shift)
{
  if (jacobian)
  {
    gsVector<T> dx = gsVector<T>::Zero(x.size());
    m_jacMat = this->computeJacobian(x,dx);
    this->factorizeMatrix(m_jacMat);
  } // otherwise the jacobian is already computed (on m_U+m_DeltaU)

  // gsInfo<<"x = \n"<<x.transpose()<<"\n";
  if (m_bifurcationMethod == bifmethod::Determinant)
  {
    if ( auto * s = dynamic_cast<typename gsSparseSolver<T>::SimplicialLDLT*>(m_solver.get()) )
    {
      factorizeMatrix(m_jacMat);
      m_stabilityVec = s->vectorD();
    }
    else
    {
      gsWarn<<"Determinant stability method only works with SimplicialLDLT solver, current solver is "<<m_options.getString("Solver")<<"\n";
      throw 3;
    }
  }
  else if (m_bifurcationMethod == bifmethod::Eigenvalue)
  {
    #ifdef gsSpectra_ENABLED
    index_t number = std::min(static_cast<index_t>(std::floor(m_jacMat.cols()/5.)),10);
    /*
    // Without shift!
    // This one can sometimes not converge, because spectra is better at finding large values.
      gsSpectraSymSolver<gsSparseMatrix<T>> es(m_jacMat,number,5*number);
      es.init();
      es.compute(Spectra::SortRule::SmallestAlge,1000,1e-6,Spectra::SortRule::SmallestAlge);
      GISMO_ASSERT(es.info()==Spectra::CompInfo::Successful,"Spectra did not converge!"); // Reason for not converging can be due to the value of ncv (last input in the class member), which is too low.
    */

    // With shift!
    // This one converges easier. However, a shift must be provided!
    gsSpectraSymShiftSolver<gsSparseMatrix<T>> es(m_jacMat,number,5*number,shift);
    es.init();
    es.compute(Spectra::SortRule::LargestAlge,1000,1e-6,Spectra::SortRule::SmallestAlge);
    if (es.info()!=Spectra::CompInfo::Successful)
    {
      gsWarn<<"Spectra did not converge!\n"; // Reason for not converging can be due to the value of ncv (last input in the class member), which is too low.
      throw 3;
    }

    // if (es.info()==Spectra::CompInfo::NotComputed)
    // if (es.info()==Spectra::CompInfo::NotConverging)
    // if (es.info()==Spectra::CompInfo::NumericalIssue)
    // gsEigen::SelfAdjointEigenSolver< gsMatrix<T> > es(m_jacMat);
    m_stabilityVec = es.eigenvalues();
    #else
    GISMO_UNUSED(shift);
    gsEigen::SelfAdjointEigenSolver<gsMatrix<T>> es2(m_jacMat);
    m_stabilityVec = es2.eigenvalues();
    #endif
  }
  else if (m_bifurcationMethod == bifmethod::Nothing)
  {
    m_stabilityVec = gsVector<T>::Zero(x.size());
  }
  else
    gsInfo<<"bifurcation method unknown!";

  m_negatives = countNegatives(m_stabilityVec);
  m_indicator = m_stabilityVec.colwise().minCoeff()[0]; // This is required since D does not necessarily have one column.
  // stability() is the SINGLE source of the +1/-1 convention: it reads m_indicator,
  // which was just assigned above. Duplicating the ternary here is what let this
  // write and stability() drift into opposite polarities.
  m_stability = this->stability();

  // ---- AUTO-07p fold test function (toolboxae.f90:703-747); see foldTestFunction().
  // Only on the jacobian=true path: there m_solver holds the factorization of the
  // tangent AT x, which is the only state for which this number means anything. On the
  // jacobian=false path (the per-iteration calls from _step()) it is set to NaN rather
  // than left stale, so a reader can never mistake an old point's value for this one's.
  m_foldTF = std::numeric_limits<T>::quiet_NaN();
  if (jacobian)
  {
    // The FROZEN step forcing, exactly as gsALMCrisfield's constraint metric uses it
    // (gsALMCrisfield.hpp:765). NOT computeForcing(): re-entering the forcing callback
    // here would overwrite m_forcingEval, whose reference semantics are documented at
    // gsALMBase.h (computeForcing: "overwritten by the next call"), and would add a
    // throw-2 path to an INSTRUMENT. Without a callback this IS m_forcing, bit-exactly.
    const gsVector<T> & f = this->stepForcing();
    try
    {
      const gsVector<T> zt  = this->solveSystem(f);   // one back-substitution
      const T           num = f.dot(zt);
      const T           den = math::sqrt((T)1 + zt.dot(zt));
      if (math::isfinite(num) && math::isfinite(den) && den > (T)0)
        m_foldTF = ((num < (T)0) ? (T)(-1) : (T)1) / den;
    }
    catch (...) { /* instrument: a failed back-substitution leaves NaN, nothing else */ }
  }
}

template <class T>
index_t gsALMBase<T>::stability() const
{
  // Convention (documented in gsALMBase.h): -1 unstable, +1 stable.
  //
  // m_indicator is in ALL cases the MINIMUM entry of m_stabilityVec (_computeStability
  // takes minCoeff), never a product -- the name "Determinant" below is the enum's, not
  // the quantity's:
  //   BifurcationMethod = Eigenvalue   m_stabilityVec = eigenvalues of K_T
  //                                    => m_indicator = smallest eigenvalue
  //   BifurcationMethod = Determinant  m_stabilityVec = the LDLT pivot vector D
  //                                    => m_indicator = smallest pivot
  //   BifurcationMethod = Nothing      m_stabilityVec = 0 => m_indicator = 0 => stable
  // In the Determinant case only the SIGN of the minimum pivot is meaningful (its
  // magnitude is not an eigenvalue), and that sign is enough: by Sylvester's law of
  // inertia the LDLT pivots carry the same inertia as the eigenvalues, so a negative
  // minimum pivot holds iff K_T has a negative eigenvalue.
  //
  // Hence under BOTH methods m_indicator < 0 iff the tangent is indefinite, i.e.
  // UNSTABLE. Note this relies on the quantity being a MINIMUM: an actual determinant
  // (product of pivots) is POSITIVE for two negative eigenvalues and would misreport
  // such a point as stable.
  return (m_indicator < 0) ? -1 : 1;
}

template <class T>
bool gsALMBase<T>::stabilityChange() const
{
  return (m_stability*m_stabilityPrev < 0) ? true : false;
}


template <class T>
bool gsALMBase<T>::_extendedSystemSolve(const gsVector<T> & U, const T L, const T tol)
{
  m_U = U;
  m_L = L;
  gsInfo<<"Extended iterations --- Starting with U.norm = "<<m_U.norm()<<" and L = "<<m_L<<"\n";

  // Evaluate the tangent WITHOUT factorizing it: at the converged singular
  // point the tangent is (bit-exactly) singular, so factorizing it would abort
  // (SimplicialLDLT NumericalIssue) even though we only need the matrix here for
  // the residual norm (m_jacMat*V). computeJacobian() factorizes internally, so
  // it cannot be used at a singular iterate; call the assembly functor directly.
  auto evalJacobian = [this](const gsVector<T> & Uarg, const gsVector<T> & dUarg)
  {
    gsSparseMatrix<T> J;
    m_note += "J";
    if (!m_djacobian(Uarg,dUarg,J))
      throw 2;
    return J;
  };

  m_jacMat = evalJacobian(m_U,m_deltaU); // Jacobian evaluated on m_U (unfactorized)
  m_basisResidualKTPhi = (m_jacMat*m_V).norm();

  m_DeltaV = gsVector<T>::Zero(m_numDof);
  m_DeltaU.setZero();
  m_DeltaL = 0.0;

  m_deltaV = gsVector<T>::Zero(m_numDof);
  m_deltaU.setZero();
  m_deltaL = 0.0;
  // With m_maxIterations < 1 the loop below never
  // runs its body, and the non-convergence branch at m_numIterations==m_maxIterations-1 is
  // unreachable, so without this guard the function would fall through to `return true;`
  // -- reporting CONVERGENCE after doing nothing. A zero (or negative) iteration budget is never a
  // converged solve; every increment above is already zero, so nothing else here needs
  // to change. This return flows to _computeSingularPoint's `converged` (gsALMBase.hpp
  // ~609): `m_converged=false` and `throw 1` -> gsStatus::NotConverged, which is the
  // correct "nothing was solved" report; the `switchBranch()` gate there
  // (`converged || m_SPfail==1`) still fires whenever SingularPointFailure==With (the
  // default), unchanged for every caller that leaves the option at its default.
  if (m_maxIterations < 1)
    return false;
  if (m_verbose)
      _initOutputExtended();
  for (m_numIterations = 0; m_numIterations < m_maxIterations; ++m_numIterations)
  {
    _extendedSystemIteration();
    m_DeltaU += m_deltaU;
    m_DeltaL += m_deltaL;
    m_DeltaV += m_deltaV;
    // m_V.normalize();

    // m_resVec = m_residualFun(m_U, m_L, m_forcing);
    // m_residue = m_resVec.norm() / ( m_L * m_forcing.norm() );
    // m_residue = (m_jacobian(m_U).toDense()*m_V).norm() / refError;
    m_jacMat = evalJacobian(m_U+m_DeltaU,m_deltaU); // unfactorized: norm-only use
    m_residueKTPhi = (m_jacMat*(m_V+m_DeltaV)).norm(); // /m_basisResidualKTPhi;
    m_resVec = this->computeResidual(m_U+m_DeltaU,m_L+m_DeltaL);
    computeResidualNorms(true); // extended solve: see the floor derivation in computeResidualNorms
    if (m_verbose)
      _stepOutputExtended();

    // termination criteria
    // Fix C (opt-in via "SingularPointComposite"): when enabled, the converged
    // singular point must also lie on the equilibrium manifold (TolF/TolU), not
    // merely where ||K.V||->0. That set is generally a whole CURVE through the singular
    // point -- on the modified-Bratu benchmark the tangent is singular along
    // mu*exp(c) = c_1 for EVERY c -- so the ||K.V||-only test constrains the singularity
    // but places NO bound at all on ||R(U,L)||; a caller that PUBLISHES the returned
    // point (gsALMExploration does) has to know which of the two it got. Default OFF
    // => the test reduces to the plain ||K.V||-only criterion, unchanged for every caller
    // that leaves SingularPointComposite at its default.
    if ( m_residueKTPhi < tol &&
         (!m_SPComposite || (m_residueF < m_toleranceF && m_residueU < m_toleranceU)) )
    {
      m_U += m_DeltaU;
      m_L += m_DeltaL;
      // Commit the MODE correction as well. The extended system iterates on
      // (U,Lambda,V) and the convergence test two lines up is ||K_T*(V+DeltaV)||, i.e. it
      // is the CORRECTED mode that was certified; committing only (U,Lambda) would make
      // solutionV() return the pre-refinement power-iteration mode, so switchBranch() and
      // gsALMExploration's branch jobs would nudge along an unrefined direction. |V+DeltaV| = 1
      // to the same tolerance (the extended constraint enforces it), so the consumers'
      // normalisation is unaffected.
      m_V += m_DeltaV;
      gsInfo<<"Iterations finished. U.norm() = "<<m_U.norm()<<"\t L = "<<m_L<<"\n";
      break;
    }
    if (m_numIterations == m_maxIterations-1)
    {
      gsInfo<<"Warning: Extended iterations did not converge! \n";
      // The two branches below are DELIBERATE,
      // option-selected, default-ON policy, not an accidental leak of a failed solve's
      // state. "SingularPointFailure" == SPfail::With (== 1, the default, gsALMBase.hpp
      // ~line 42) commits m_U/m_L/m_V from the last (non-converged) increment;
      // == SPfail::Without (== 0) discards the increment and keeps the entry state.
      // _computeSingularPoint (gsALMBase.hpp ~line 622) relies on the committed state:
      // it calls switchBranch() whenever (converged || m_SPfail==1), i.e. the "switch
      // anyway on failure" behaviour is intentional and documented there. The caller
      // already receives the failure through m_converged=false and the thrown status
      // (gsStatus::NotConverged); gsALMExploration additionally restores its own U/L on
      // that status but does NOT undo m_V. Changing this block (guard or save/restore)
      // is therefore a behavioural change reaching every default-configuration ALM
      // consumer, and must not be made without a deliberate decision to do so --
      // see TEST(...) in gsALMSolvers_test.cpp for the pinning test that keeps this
      // legible and non-drifting. Do not "fix" this without a fresh decision.
      if (m_SPfail==1)
      {
        m_U += m_DeltaU;
        m_L += m_DeltaL;
        // Same mode commit as the converged branch above. GUARDED here because
        // this branch is reached on a DIVERGED extended solve, where nothing constrains
        // |V+DeltaV| to 1: it may be ~0 or non-finite, and both consumers of solutionV()
        // (switchBranch, gsALMExploration) divide by that norm. When the correction is
        // unusable the pre-refinement mode - the best direction available - is kept.
        const T Vnorm = (m_V+m_DeltaV).norm();
        if (Vnorm > (T)0 && math::isfinite(Vnorm))
          m_V += m_DeltaV;
        else
          gsWarn<<"Extended system: the mode correction gives |V+DV| = "<<Vnorm
                <<"; keeping the unrefined mode.\n";
        gsInfo<<"Iterations finished; continuing with last solution U.norm() = "<<m_U.norm()<<"\t L = "<<m_L<<"\n";
      }
      else
      {
        m_DeltaV = gsVector<T>::Zero(m_numDof);
        m_DeltaU.setZero();
        m_DeltaL = 0.0;

        m_deltaV = gsVector<T>::Zero(m_numDof);
        m_deltaU.setZero();
        m_deltaL = 0.0;
        gsInfo<<"Iterations finished. continuing with original solution U.norm() = "<<m_U.norm()<<"\t L = "<<m_L<<"\n";
      }
      return false;
    }
  }
  return true;
}

// TODO: Optimize memory
template <class T>
void gsALMBase<T>::_extendedSystemIteration()
{
  // Evaluate the tangent WITHOUT factorizing (computeJacobian() factorizes
  // internally and aborts on a singular tangent). Matrices are only assembled
  // here; factorization is done separately below via factorizeShifted().
  auto evalJacobian = [this](const gsVector<T> & Uarg, const gsVector<T> & dUarg)
  {
    gsSparseMatrix<T> J;
    m_note += "J";
    if (!m_djacobian(Uarg,dUarg,J))
      throw 2;
    return J;
  };

  // Load-bearing factorization of the (possibly exactly-singular) extended-system
  // tangent. Near the singular point the *bordered* system is well conditioned;
  // a tiny diagonal shift on the inner solves acts as an inexact-Newton
  // perturbation and keeps the converged singular point within tolerance (the
  // convergence test in _extendedSystemSolve uses the UNshifted tangent). We
  // factorize a shifted COPY so m_jacMat itself stays unshifted. On any tangent
  // that factorizes cleanly (the common/shell path) no shift is applied and the
  // result is bit-identical to before.
  auto factorizeShifted = [this](const gsSparseMatrix<T> & M)
  {
    try { this->factorizeMatrix(M); return; }
    catch (...) { /* singular: fall through to shifted retry */ }

    T dmax = 1;
    for (index_t k = 0; k < M.rows(); ++k)
      dmax = math::max(dmax, math::abs(M.coeff(k,k)));
    T sigma = static_cast<T>(1e-12) * dmax;

    gsSparseMatrix<T> Id(M.rows(), M.cols());
    Id.setIdentity();
    for (index_t attempt = 0; attempt < 3; ++attempt, sigma *= static_cast<T>(100))
    {
      try
      {
        gsSparseMatrix<T> Ms = M + sigma*Id;
        this->factorizeMatrix(Ms);
        gsInfo<<"Extended system: (near-)singular tangent; factorized with diagonal shift sigma = "<<sigma<<"\n";
        return;
      }
      catch (...) { /* escalate shift */ }
    }
    // Exhausted the shift attempts: rethrow the original failure loudly.
    this->factorizeMatrix(M);
  };

  m_resVec = this->computeResidual(m_U+m_DeltaU, m_L+m_DeltaL);
  m_jacMat = evalJacobian(m_U+m_DeltaU,m_deltaU); // unfactorized tangent
  factorizeShifted(m_jacMat);

  // m_jacMat = m_jacobian(m_U);

  // Forcing at the state the extended-system tangent was assembled at
  // (m_U+m_DeltaU, m_L+m_DeltaL). NOTE: the SAME f must be reused in h1 below, and no
  // other forcing evaluation may intervene (computeForcing() returns a reference into a
  // single scratch vector).
  const gsVector<T> & F = this->computeForcing();
  m_deltaUt = this->solveSystem(F); // DeltaV1
  m_deltaUbar = this->solveSystem(-m_resVec); // DeltaV2

  real_t eps = 1e-8;
  // jacMatEps is only used in matvecs (h1,h2); evaluate it without factorizing
  // (this also avoids clobbering m_solver, which still holds m_jacMat's factors).
  gsVector<T> zeroDelta = gsVector<T>::Zero((m_U+m_DeltaU).size());
  gsSparseMatrix<T> jacMatEps = evalJacobian((m_U+m_DeltaU) + eps*(m_V+m_DeltaV), zeroDelta);
  // h1 approximates d(K_T*delta_u_t)/dU . V by finite differences. The subtracted term
  // stands in for K_T(U)*delta_u_t, which equals f EXACTLY by the solve above — hence f
  // here must be the forcing of that solve, NOT f(U+eps*V); using the perturbed forcing
  // would corrupt the difference quotient.
  gsVector<T> h1 = 1/eps * ( jacMatEps * m_deltaUt ) - 1/eps * F;
  gsVector<T> h2 = m_jacMat * (m_V+m_DeltaV) + 1/eps * ( jacMatEps * m_deltaUbar + m_resVec );

  factorizeShifted(m_jacMat);

  m_deltaVt = this->solveSystem(-h1); // DeltaV1
  m_deltaVbar = this->solveSystem(-h2); // DeltaV2

  m_deltaL = -( ( (m_V+m_DeltaV)/(m_V+m_DeltaV).norm() ).dot(m_deltaVbar)  + (m_V+m_DeltaV).norm() - 1) / ( (m_V+m_DeltaV)/(m_V+m_DeltaV).norm() ).dot( m_deltaVt );

  m_deltaU = m_deltaL * m_deltaUt + m_deltaUbar;
  // gsInfo<<"m_DeltaU = \n"<<m_DeltaU<<"\n";

  m_deltaV = m_deltaL * m_deltaVt + m_deltaVbar;
  // gsInfo<<"m_DeltaV = \n"<<m_DeltaV<<"\n";
}

template <class T>
index_t gsALMBase<T>::_bisectionObjectiveFunction(const gsVector<T> & x, bool jacobian)
{
  this->_computeStability(x,jacobian);
  m_note += "(" + std::to_string(m_negatives) + ")";
  return m_negatives;
}

template <class T>
T gsALMBase<T>::_bisectionTerminationFunction(const gsVector<T> & x, bool jacobian)
{
  this->_computeStability(x,jacobian);
  // return m_stabilityVec.colwise().minCoeff()[0]; // This is required since D does not necessarily have one column.
  return m_indicator;
}

template <class T>
bool gsALMBase<T>::_bisectionSolve(const gsVector<T> & U, const T L, const T tol,
                                    const gsVector<T> & Ufar, const T Lfar)
{
  m_U = U;
  m_L = L;

  // Entry state, restored on every non-converged exit (see the doxygen restore guarantee).
  const gsVector<T> U_entry     = U;
  const T           L_entry     = L;
  const gsVector<T> Uprev_entry = m_Uprev;
  const T           Lprev_entry = m_Lprev;
  // Store original arc length
  const real_t dL = m_arcLength; // NOTE (pre-existing): real_t, not T - a narrowing for non-double T
  // The step()s below overwrite m_arcLength_prev with a bisection-INTERNAL length. Before
  // m_arcLength_prev became persistent (see setLength()), the consumer's next setLength()
  // happened to repair that; now nothing does, so save and restore it here as well.
  const T dLprev = m_arcLength_prev;
  const bool adaptiveBool = m_adaptiveLength;
  m_adaptiveLength = false;
  bool converged = false;
  m_SPBisProbes = 0;

  // Guard for a missing incumbent (e.g. the solver's state was never set: m_U starts as a
  // size-0 vector) -- there is nothing to evaluate Ufar's stability against, so this counts
  // as no bracket without touching Ufar at all.
  const bool haveIncumbent = (Ufar.size() == U.size());

  // Endpoint evaluation order is load-bearing (see the m_V postcondition in the doxygen):
  // _computeCriticalMode below reuses the tangent/factorization left resident by the LAST
  // _computeStability call, which must describe (U,L), not (Ufar,Lfar). Evaluate the far
  // endpoint first (when there is one), the near endpoint (U,L) last.
  index_t negHi = 0;
  if (haveIncumbent)
    negHi = _bisectionObjectiveFunction(Ufar, true);

  const T referenceError = _bisectionTerminationFunction(U, true);
  if (m_verbose)
    gsInfo<<"Bisection iterations --- Starting with U.norm = "<<m_U.norm()<<" and L = "<<m_L<<"; Reference error = "<<referenceError<<"\n";
  const index_t negLo = _bisectionObjectiveFunction(U, false); // jacobian already computed above

  const bool haveBracket = haveIncumbent && (negLo != negHi);

  if (haveBracket)
  {
    // Bracket on the arc-length parameter s in [0,dLb], measured from (U,L) (the
    // pre-crossing side) to (Ufar,Lfar) (the post-crossing side, s = dLb = m_arcLength at
    // entry). Invariant: negLo != negHi. In the shape of
    // gsALMExploration<T>::_localizeCrossing.
    const T dLb = dL;
    T slo = (T)0, shi = dLb;
    T sprobe = (slo + shi) / (T)2;
    index_t negLoLoop = negLo;

    for (index_t k = 1; k <= m_SPBisIt; ++k)
    {
      if (m_verbose) gsInfo<<"\t bisection iteration "<<k<<"\t arc length = "<<sprobe<<"; ";

      // Reset start point: re-seed from the FIXED near endpoint and step to the current
      // midpoint, exactly like _localizeCrossing's own step-fail retry
      // (m_solver->setSolution(Uold,Lold) on every probe there). Re-seed from U_entry/
      // L_entry, NOT from the U/L parameters: the one in-tree call site
      // (_computeSingularPoint) passes m_U/m_L themselves, so U/L ALIAS m_U/m_L and no
      // longer denote the entry point once the first step() below has moved m_U.
      // U_entry/L_entry are the value copies taken before that first mutation, so sprobe
      // is thereby a genuine coordinate in [0,dLb] measured from the entry point; slo/shi
      // are bracket COORDINATES only, used to pick the next probe -- the point actually
      // returned on convergence is the final probe itself, (m_U,m_L), not a bracket
      // endpoint reconstructed from slo.
      m_U = U_entry;
      m_L = L_entry;
      m_arcLength = sprobe;

      if (m_verbose) gsInfo<<"From U.norm = "<<m_U.norm()<<" and L = "<<m_L<<"\n";

      gsStatus status = step();
      ++m_SPBisProbes;

      // Fresh tangent and inertia at the probe (jacobian=true), regardless of status: the
      // termination check below reuses this SAME factorization via jacobian=false.
      const index_t negProbe = _bisectionObjectiveFunction(m_U, true);

      if (status == gsStatus::Success)
      {
        if (negProbe == negLoLoop)
        {
          slo = sprobe; negLoLoop = negProbe;
        }
        else
        {
          // Keeps the invariant negLo != negHi when a multi-step crossing (e.g. 1->3) is
          // probed at an intermediate inertia (e.g. 2).
          shi = sprobe;
        }
        sprobe = (slo + shi) / (T)2;
      }
      else
      {
        // Non-converged probe: a probe with no inertia cannot certify a half-interval, so
        // neither bracket endpoint moves. Retreat the probe toward the low end instead.
        sprobe = (slo + sprobe) / (T)2;
      }

      // termination criteria: relative error
      T term = _bisectionTerminationFunction(m_U, false); // jacobian on the new point already computed

      if (m_verbose) gsInfo<<"\t Finished. Relative error = "<< abs(term/referenceError)<<"\t"<<" obj.value = "<<term<<"\n";

      // A non-converged probe leaves m_U at an unconverged Newton iterate; its indicator
      // can be small by accident, so it must not be allowed to certify convergence (that
      // would skip the restore below and hand _extendedSystemSolve a spurious point).
      if ( status == gsStatus::Success && abs(term/referenceError) < tol )
      {
          converged = true;
          break;
      }
    }
  }

  // Reset arc length (and the secant divisor the internal step()s overwrote), for every
  // path -- no-bracket, budget-exhausted or converged.
  m_arcLength = dL;
  m_arcLength_prev = dLprev;
  m_adaptiveLength = adaptiveBool;

  if (!converged)
  {
    m_U = U_entry;
    m_L = L_entry;
    m_Uprev = Uprev_entry;
    m_Lprev = Lprev_entry;
    // On the no-bracket path (haveBracket == false) nothing moved and the tangent left
    // resident by the endpoint evaluations above already describes (U,L), by construction
    // of the evaluation order. On a budget-exhausted bracket attempt the tangent instead
    // describes the last probed point, not the just-restored (U,L) -- refresh it so the
    // m_V postcondition below holds here too.
    if (haveBracket)
      this->_computeStability(m_U, true);
  }

  // Critical mode at the point left in (m_U,m_L). On the testPoint=false path (no
  // classification) this is the mode that seeds the extended system (and hence every
  // branch nudge); on the testPoint=true path, classification runs after this call and
  // _testSingularPoint recomputes and OVERWRITES it instead, so it MUST use the same
  // iteration: this call routes through _computeCriticalMode rather than re-implementing
  // it, which avoids duplicating the deterministic generic start (see START VECTOR
  // above), hardcoding a separate sweep count, using a sign-blind convergence test, and
  // refactorizing the same matrix redundantly inside the loop.
  // Tolerance is m_SPModeTol (the mode's OWN resolution), not m_tolerance, matching the
  // _testSingularPoint call sites -- see the SingularPointModeTol option. This call SEEDS
  // the extended system only; it does not classify, so it must NOT write m_SPverdict
  // (only _testSingularPoint does). It may still update m_SPModeError, which is a property
  // of the last mode iteration rather than of a classification.
  // NOTE: `converged` keeps the value the BISECTION loop above set. The eigenvector
  // iteration's own convergence flag is discarded above (its return value is not
  // captured): that flag can settle to true independently of whether the bisection
  // converged, so folding it into `converged` here would let a converged mode mask a
  // failed bisection. This is why the only caller (_computeSingularPoint) documents that
  // it ignores this function's return value.
  this->_computeCriticalMode(m_SPTestIt, m_SPModeTol);
  return converged;
}

template <class T>
void gsALMBase<T>::switchBranch()
{
  // FALLBACK branch-switch predictor: used only when computeBranchTangent()
  // does not return branchTangent::Success/TrivialBranch, i.e. when the true swibra ABE
  // emanating tangent cannot be resolved. Nudge the singular point along the UNIT
  // critical mode, U <- U* + xi*V, with xi = 1/m_tau. Note: an expression of the form
  //     m_V.normalize();  real_t lenPhi = m_V.norm();  real_t xi = lenPhi/m_tau;
  // would read the norm of a just-normalised vector, i.e. IDENTICALLY 1 -- so the apparent
  // "scale the nudge with the mode length" would do nothing and xi would always be
  // 1/m_tau regardless. That is the source of the project's tau-sensitivity folklore: the
  // nudge magnitude is set by Perturbation (m_tau) ALONE. Writing xi directly as 1/m_tau,
  // as below, carries no loss of information: every producer of m_V returns it normalised
  // (the power iteration normalises each sweep, the extended system constrains |V| = 1),
  // and it matches gsALMExploration, which builds its branch jobs as U* + (1/tau)*V/|V|.
  m_V.normalize();
  const T xi = (T)1/m_tau;
  m_DeltaU = xi*m_V;
  m_U = m_U + m_DeltaU;

  // Both staleness sentinels must be cleared, not just one. gsALMCrisfield selects its
  // fresh-start predictor on (m_DeltaUold,m_DeltaLold); gsALMRiks and
  // gsALMConsistentCrisfield select theirs on the SECANT (m_U-m_Uprev, m_L-m_Lprev).
  // Clearing only the first left those two extrapolating along the PRE-bifurcation
  // direction of travel after a branch switch, so the corrector - constrained only to the
  // arc-length sphere - was free to fall back onto the fundamental path. Assigning the
  // (already nudged) current point to the previous one selects their fresh-start branch,
  // which is exactly what zeroing m_DeltaUold does for gsALMCrisfield.
  // @note A caller that issues setPrevious() AFTER switchBranch() re-installs a secant and
  //       defeats this reset; the in-tree callers only follow with setLength(), which does
  //       not touch (m_Uprev,m_Lprev) and (once a step has been taken) not m_arcLength_prev
  //       either, so the fix composes with the setLength()/setPrevious() rules.
  m_DeltaLold = 0.0;
  m_DeltaUold.setZero();
  m_Uprev = m_U;
  m_Lprev = m_L;
  // m_DeltaLold = m_DeltaL;
  // m_DeltaUold = m_DeltaU;
}

// TRUE emanating branch tangent at a simple (multiplicity-1) branch point, by the
// algebraic branching equation (ABE) of pde2path's swibra.m (lines 31,41-72). Ported
// verbatim, including the two things the summary that motivated this task omitted:
//  (1) the "if any(Gud)" TRIVIAL-branch case (swibra.m:61,68-70): K_T not varying along
//      phi1 (a LINEAR problem) is not a failure, tau1 = [phi1;0];
//  (2) del = 1e-3 with pde2path's own warning "don't choose del too small!".
//
// Sign convention: pde2path's residual G is our R; Gu is our K_T; Glam = dR/dLambda = -f,
// where f = computeForcing() (see its doxygen / setForcingFunction). Hence
// Glamd = -( f(U*+del*phi1,L*) - f(U*,L*) ), and it is IDENTICALLY ZERO on the dead-load path.
template <class T>
typename gsALMBase<T>::branchTangent::type
gsALMBase<T>::computeBranchTangent(const gsVector<T> & Ustar, const T & Lstar,
                                    const gsVector<T> & phi1in,
                                    const gsVector<T> & tangentU0, const T & tangentL0,
                                    gsVector<T> & tau1U, T & tau1L)
{
  // Euclidean fallback for the solver-specific distance() metric, used throughout this
  // function wherever a normalisation is needed: distance() is GISMO_NO_IMPLEMENTATION on
  // the base class (throws std::runtime_error) and is degenerate (identically zero) for
  // gsALMLoadControl on a pure-displacement vector, so every normalisation here must
  // tolerate both a throw and a non-positive/non-finite result.
  //
  // DEVIATION from the spec's literal guard ("fall back ... when it throws or returns a
  // value that is not > 0 / not finite"), MEASURED necessary: gsALMLoadControl's distance()
  // is |DeltaL| ALONE (gsALMLoadControl.h:84-91, its OWN corrector-constraint metric, by
  // design), which is POSITIVE and FINITE yet can be orders of magnitude smaller than the
  // vector's Euclidean length whenever the load component is small relative to the
  // displacement one -- exactly the pitchfork oracle's tau1L = 3e-3 against a raw vector of
  // Euclidean length ~2 (ratio 1.5e-3). Passing the literal guard there still normalised a
  // UNIT STEP DIRECTION (consumed by gsALMExploration::traceSweep as
  // job.U + s*baseLen*job.tangentU) by 3e-3, inflating the displacement component to ~667 --
  // MEASURED -- turning a small Euler predictor step into an O(1) state jump instead. The
  // spec calls distance() "the analogue of pde2path's xinorm", a full-state metric; a metric
  // that measures one near-zero component of an otherwise O(1) vector is not commensurate
  // with it. ratioFloor = 1e-2 stays clear of legitimate WEIGHTED metrics (Crisfield's
  // sqrt(dU.dot(dU)+A0*dL^2), Riks' sqrt(w*|dU|^2+(1-w)*dL^2)), which track the Euclidean
  // norm within an order of magnitude for any non-degenerate A0/w, while catching a metric
  // that measures a near-zero component of an O(1) vector.
  auto safeDistance = [this](const gsVector<T> & dU, const T dL) -> T
  {
    const T euclid = math::sqrt(dU.squaredNorm() + dL*dL);
    T d;
    try { d = this->distance(dU, dL); }
    catch (...) { d = std::numeric_limits<T>::quiet_NaN(); }
    const T ratioFloor = (T)1e-2;
    if (!(d > (T)0) || !math::isfinite(d) || !(d >= ratioFloor*euclid))
      d = euclid;
    return d;
  };

  // ---- 1. phi1: LOCAL, normalised copy (never read from m_V: the caller already holds
  // it, e.g. solutionV(), so this contract carries no hidden ordering requirement). ------
  gsVector<T> phi1 = phi1in;
  {
    const T nrm = phi1.norm();
    if (!(nrm > (T)0) || !math::isfinite(nrm))
    {
      if (m_verbose)
        gsInfo<<"computeBranchTangent: ModeUnresolved (|phi1| = "<<nrm<<", non-positive or non-finite).\n";
      return branchTangent::ModeUnresolved;
    }
    phi1 /= nrm;
  }

  // Assemble K_T(Ustar) ONCE, through m_djacobian DIRECTLY -- NOT computeJacobian(), which
  // factorizes and throw-3's on the tangent that is singular by construction at U*. Reused
  // below both for the left-mode computation (step 2) and for Gud (step 5).
  gsSparseMatrix<T> Kstar;
  {
    const gsVector<T> zeroDelta = gsVector<T>::Zero(Ustar.size());
    m_note += "J";
    // computeBranchTangent is PUBLIC and must never let a raw int (the in-tree
    // assembly-failure convention, e.g. throw 2) escape -- see round-1 review fix 3.
    bool ok = false;
    try { ok = m_djacobian(Ustar, zeroDelta, Kstar); }
    catch (...) { ok = false; }
    if (!ok)
    {
      if (m_verbose) gsInfo<<"computeBranchTangent: AssemblyFailed (K_T(U*)).\n";
      return branchTangent::AssemblyFailed;
    }
  }
  // ||K_T(U*)||, computed once and reused by the TrivialBranch test's epsTriv comparison
  // (step 5).
  const T KstarNorm = Kstar.norm();

  // ---- 2. psi1: LEFT critical mode of K_T(Ustar). ---------------------------
  //
  // FRESHNESS: the critical-mode
  // pair the caller may otherwise be tempted to reuse is NOT classification-fresh --
  // _bisectionSolve tail-calls _computeCriticalMode, which overwrites m_V/m_SPModeError
  // for a DIFFERENT point than the one this function is asked about. phi1 sidesteps that
  // entirely (it is the caller's own, explicitly passed, externally certified mode -- the
  // converged extended-system solve's ||K_T.phi1|| < SingularPointComputeTolE). For psi1 we
  // seed m_V with THIS phi1 before calling _computeCriticalModeLeft: on the default
  // SELF-ADJOINT (SimplicialLDLT) path that routine performs NO factorization at all and
  // returns m_Vleft = m_V exactly (psi1=phi1 is then a mathematical IDENTITY under
  // symmetry, taken openly and not silently -- we are
  // calling the library's own routine, which decides that for itself via
  // _solverIsSelfAdjoint()), so its own convergence check is the ONLY gate left on that
  // path -- and that check reads m_SPModeError. See the VACUOUS-GUARD note below for why,
  // after two attempts at a phi1-quality-derived seed both measured to misfire, this
  // function deliberately does NOT try to make that check a function of phi1's quality.
  //
  // VACUOUS-GUARD. A left-mode freshness guard on the self-adjoint path needs a
  // residual to gate on, and neither direct candidate works:
  //  - An ABSOLUTE residual ||Kstar*phi1|| is reachable but not SCALE-INVARIANT
  //    against m_SPModeTol's documented DIMENSIONLESS meaning (gsALMBase.hpp:69-73,
  //    gsALMBase.h:675-690) -- MEASURED to flip to ModeUnresolved on fixture P
  //    rescaled by s>=10 alone, phi1's direction error fixed.
  //  - A RELATIVE residual ||Kstar*phi1||/||Kstar|| IS scale-invariant (re-verified:
  //    fixture P s-sweep unchanged) but is compared against m_SPModeTol -- an option
  //    UNRELATED to whatever tolerance actually certified phi1 upstream
  //    (SingularPointComputeTolE/B, an extended-system or bisection-level ABSOLUTE
  //    test). MEASURED on example_ModifiedBratuExploration --solver SimplicialLDLT:
  //    the accepted phi1 has ||Kstar||=225.2, relative residual 1.312e-07,
  //    comfortably ABOVE m_SPModeTol=1e-08 -- but the point was accepted upstream on
  //    a LOOSER test ("accepted on ||K_T.V|| alone", the driver's own log), not on
  //    SingularPointComputeTolE=1e-8: back-computing, the ABSOLUTE residual actually
  //    certified is ~3e-5, i.e. m_SPModeTol demands ~3 ORDERS OF MAGNITUDE better than
  //    anything upstream ever promised for this phi1. The relative form is therefore
  //    STILL a knife-edge against an unrelated certificate -- just a differently-shaped
  //    one -- and turns the tangent path off on the project's own self-adjoint-solver
  //    driver just as an absolute-residual gate would. Deriving a floor that is
  //    genuinely commensurate with SingularPointComputeTolE/B (which are extended-
  //    system/bisection options this function does not otherwise read, and are
  //    themselves ABSOLUTE tests against different quantities -- ||K_T.V|| vs a
  //    bracket width) is a larger redesign than this function's current scope; see the
  //    FOLLOW-UP below.
  //
  //  DECISION: state precisely why the check cannot fail and remove the dead check
  //  rather than leave a test that only pretends to run. m_SPModeError is seeded to
  //  (T)0: on the SELF-ADJOINT fast path this makes
  //  ModeUnresolved-from-the-left-mode STRUCTURALLY UNREACHABLE, by design, not by
  //  accident -- _computeCriticalModeLeft's self-adjoint branch performs NO
  //  factorization and returns m_Vleft = m_V exactly (psi1 = phi1, an IDENTITY under
  //  symmetry, not a shortcut this function takes -- the library's own routine
  //  decides that via _solverIsSelfAdjoint()), so there is no independent residual
  //  this function can responsibly gate on without either (a) inventing a threshold
  //  with no principled tie to what actually certified phi1 (measured above to be
  //  actively harmful), or (b) reading SingularPointComputeTolE/B directly, out of
  //  scope here. A phi1 of poor ACTUAL quality on the self-adjoint path is therefore
  //  not caught by this guard; it is caught downstream, if at all, by the
  //  biorthogonalisation cosine (den = phi1.psi1, degenerate to exactly 1 when
  //  psi1=phi1, so also no help here) and by the ABE's OWN degeneracy tests (steps 5,
  //  8), which operate on independently assembled K_T(U*+del*phi1) and are NOT
  //  subject to this same knife-edge.
  //
  //  On the NON-self-adjoint path this value is only a PRIOR: _computeCriticalModeLeft
  //  there recomputes its own dist from a fresh transpose power iteration and folds it
  //  in via math::max(m_SPModeError,dist), so seeding 0 can only ever be overridden
  //  upward by that fresh measurement -- never wrong, never vacuous there.
  //
  //  FOLLOW-UP: a principled self-adjoint-path freshness guard needs a floor read from
  //  SingularPointComputeTolE/B directly (whichever certified this specific phi1), not
  //  from SingularPointModeTol, which governs a DIFFERENT quantity (the LEFT mode's
  //  OWN inverse-power convergence, never exercised on the self-adjoint fast path in
  //  the first place).
  //
  // The non-self-adjoint path recomputes psi1 from scratch on K_T(Ustar), independent of
  // any staleness concern -- but K_T(Ustar) is singular by construction, so its OWN
  // restore (factorizeMatrix(m_jacMat) at the end of the non-self-adjoint branch) is
  // expected to throw under plain SimplicialLDLT/LU there too; this function then
  // correctly reports ModeUnresolved (via the outer try/catch below) rather than
  // propagating that throw. A non-self-adjoint solver therefore falls back to
  // switchBranch()'s nudge here by default; making it resolvable would need
  // _computeCriticalModeLeft's non-self-adjoint path to tolerate a singular tangent,
  // which it does not today.
  //
  // State/purity contract: save what we are about to touch and restore it exactly,
  // whatever happens (including a throw): m_jacMat, m_V, m_Vleft, m_SPModeError and the
  // solver's factorization. The non-self-adjoint path of _computeCriticalModeLeft is
  // documented as NOT exception-safe on a throw mid-path (its own mandatory
  // factorizeMatrix(m_jacMat) restore can be skipped there); redo the restore here
  // defensively rather than trust it, since our own FD probes below perturb K_T by
  // construction and must not inherit a factorization pointing at the wrong matrix.
  const gsSparseMatrix<T> savedJacMat  = m_jacMat;
  const gsVector<T>       savedV       = m_V;
  const gsVector<T>       savedVleft   = m_Vleft;
  const T                 savedModeErr = m_SPModeError;

  gsVector<T> psi1;
  bool haveLeftMode = false;
  // Captured for diagnostics only (the ModeUnresolved message below): the value
  // _computeCriticalModeLeft actually tested m_SPModeError against, read back BEFORE the
  // restore two lines down overwrites it. NaN if the try block threw before reaching it.
  // On the self-adjoint fast path this is always the seeded 0 (see the VACUOUS-GUARD note
  // above -- deliberately unreachable there); on the non-self-adjoint path it is the REAL
  // dist from a fresh transpose power iteration, genuine diagnostic information.
  T reportedModeError = std::numeric_limits<T>::quiet_NaN();
  try
  {
    m_jacMat      = Kstar;
    m_V           = phi1;
    // Seeded 0: see the VACUOUS-GUARD note above for why this is the
    // deliberate, documented choice rather than an oversight.
    m_SPModeError = (T)0;
    haveLeftMode  = this->_computeCriticalModeLeft(m_SPTestIt, m_SPModeTol);
    psi1          = m_Vleft;
    reportedModeError = m_SPModeError;
  }
  catch (...)
  {
    haveLeftMode = false;
  }

  m_jacMat      = savedJacMat;
  m_V           = savedV;
  m_Vleft       = savedVleft;
  m_SPModeError = savedModeErr;
  // The SELF-ADJOINT fast path of _computeCriticalModeLeft (the default
  // SimplicialLDLT configuration) touches NEITHER m_jacMat NOR the solver's
  // factorization -- it only copies m_V into m_Vleft -- so there is nothing to
  // restore there, and re-factorizing savedJacMat unconditionally would fire a
  // false alarm on EVERY branch point: K_T(U*) is singular by construction there,
  // and _extendedSystemSolve's own m_jacMat is deliberately left UNFACTORIZED for
  // exactly that reason (see its comment), so savedJacMat generally cannot be
  // factorized at all, regardless of whether this function touched anything. Only
  // the NON-self-adjoint path factorizes (Kstar's transpose, then Kstar itself to
  // restore its own m_solver state); only there is a real restore owed.
  if (!this->_solverIsSelfAdjoint())
  {
    try { this->factorizeMatrix(m_jacMat); }
    catch (...)
    {
      gsWarn<<"computeBranchTangent: could not re-factorize the solver's entry tangent "
              "after the non-self-adjoint left-mode computation; the solver's cached "
              "factorization may be stale until the next assembly.\n";
    }
  }

  if (!haveLeftMode)
  {
    // Unreachable-by-design on the self-adjoint fast path (m_SPModeError seeded 0, see the
    // VACUOUS-GUARD note above), so a ModeUnresolved here is always the NON-self-adjoint
    // path: reportedModeError is the real dist from its fresh transpose power iteration.
    if (m_verbose)
      gsInfo<<"computeBranchTangent: ModeUnresolved (the left critical mode of K_T(U*) did not "
              "resolve; |dV| = "<<reportedModeError<<" vs SingularPointModeTol = "<<m_SPModeTol
            <<").\n";
    return branchTangent::ModeUnresolved;
  }

  // Biorthogonalise: den = phi1.psi1. swibra.m:57 leaves this unguarded; at a near-double
  // mode (audit B.7, the shell-wrinkling case) it explodes silently -- guard the COSINE
  // against BranchTangentTol instead.
  T den = phi1.dot(psi1);
  {
    const T scale = phi1.norm() * psi1.norm();
    if (!(math::abs(den) > m_branchTangentTol * scale))
    {
      if (m_verbose)
        gsInfo<<"computeBranchTangent: ModeUnresolved (|phi1.psi1| = "<<math::abs(den)
              <<" not above BranchTangentTol*|phi1||psi1| = "<<m_branchTangentTol*scale<<").\n";
      return branchTangent::ModeUnresolved;
    }
  }
  psi1 /= den;

  // ---- 3. Incoming tangent tau0 = (u0d,lam0), normalised to unit Euclidean length. The
  // result tau1 is invariant under a joint rescaling of tau0 (see the task's self-check),
  // so this normalisation only makes the load-share guard below a stateable number. --------
  gsVector<T> u0d  = tangentU0;
  T           lam0 = tangentL0;
  {
    const T nrm = math::sqrt(u0d.squaredNorm() + lam0*lam0);
    if (!(nrm > (T)0) || !math::isfinite(nrm))
    {
      if (m_verbose)
        gsInfo<<"computeBranchTangent: Degenerate (incoming tangent (tangentU0,tangentL0) has "
                "zero or non-finite norm).\n";
      return branchTangent::Degenerate;
    }
    u0d  /= nrm;
    lam0 /= nrm;
  }

  // Load-share guard: lambdadot0 (lam0) divides twice below (phi0 here, al1b in step 8), so
  // a branch point whose incoming tangent carries essentially no load component is a branch
  // point coinciding with a fold -- outside multiplicity-1 swibra's scope.
  const gsVector<T> zeroU = gsVector<T>::Zero(u0d.size());
  const T lamShare = safeDistance(zeroU, lam0) / safeDistance(u0d, lam0);
  if (!(lamShare > m_branchTangentTol))
  {
    if (m_verbose)
      gsInfo<<"computeBranchTangent: Degenerate (incoming tangent's load share "<<lamShare
            <<" not above BranchTangentTol = "<<m_branchTangentTol<<").\n";
    return branchTangent::Degenerate;
  }

  // ---- 4. al1, phi0. ---------------------------------------------------------------------
  const T al1 = psi1.dot(u0d);
  const gsVector<T> phi0 = (u0d - al1*phi1) / lam0;

  // ---- 5. Kdel = K_T(U*+del*phi1), Gud = Kdel-Kstar; TRIVIAL-branch test (swibra.m:61). --
  gsSparseMatrix<T> Kdel;
  {
    const gsVector<T> zeroDelta = gsVector<T>::Zero(Ustar.size());
    const gsVector<T> Udel      = Ustar + m_branchTangentDel*phi1;
    m_note += "J";
    // Same escape-hatch as the Kstar assembly above: never let a raw int escape this
    // PUBLIC method (round-1 review fix 3).
    bool ok = false;
    try { ok = m_djacobian(Udel, zeroDelta, Kdel); }
    catch (...) { ok = false; }
    if (!ok)
    {
      if (m_verbose) gsInfo<<"computeBranchTangent: AssemblyFailed (K_T(U*+del*phi1)).\n";
      return branchTangent::AssemblyFailed;
    }
  }
  const gsSparseMatrix<T> Gud = Kdel - Kstar;

  // Roundoff-scale test ("the two assemblies returned the same matrix"), NOT a tunable:
  // KNOWN LIMITATION -- both norms are Frobenius over the WHOLE matrix, so on a large
  // sparse problem where only a few entries vary along phi1 the ratio can approach
  // 100*eps for a genuinely nonlinear problem, silently returning TrivialBranch with the
  // plausible-looking [phi1;0]. Exposure is low (TrivialBranch and Success are routed
  // identically downstream), which is exactly why the outcome name is always logged below.
  const T epsTriv = (T)100 * std::numeric_limits<T>::epsilon();
  if (!(Gud.norm() > epsTriv * KstarNorm))
  {
    tau1U = phi1;
    tau1L = (T)0;
    const T nrm = safeDistance(tau1U, tau1L);
    if (nrm > (T)0 && math::isfinite(nrm)) { tau1U /= nrm; tau1L /= nrm; }
    if (m_verbose)
      gsInfo<<"computeBranchTangent: al1="<<al1<<" [TrivialBranch: |Gud|="<<Gud.norm()
            <<" <= "<<epsTriv<<"*|K_T(U*)|="<<epsTriv*KstarNorm
            <<" -- K_T does not vary along phi1 to roundoff; tau1 = [phi1;0]].\n";
    return branchTangent::TrivialBranch;
  }

  // ---- 6. Glamd = -( f(U*+del*phi1,L*) - f(U*,L*) ); zero on the dead-load path. ---------
  gsVector<T> Glamd;
  if (!m_forcingFun)
  {
    Glamd = gsVector<T>::Zero(Ustar.size());
  }
  else
  {
    gsVector<T> fStar, fDel;
    try { fStar = this->computeForcing(Ustar, Lstar); } // COPY: computeForcing() aliases a
                                                          // single scratch vector
    catch (...) // widened from catch(int): never let a raw int (or anything else) escape
                // this PUBLIC method (round-1 review fix 3).
    {
      if (m_verbose) gsInfo<<"computeBranchTangent: AssemblyFailed (forcing at U*).\n";
      return branchTangent::AssemblyFailed;
    }
    try { fDel = this->computeForcing(Ustar + m_branchTangentDel*phi1, Lstar); }
    catch (...)
    {
      if (m_verbose) gsInfo<<"computeBranchTangent: AssemblyFailed (forcing at U*+del*phi1).\n";
      return branchTangent::AssemblyFailed;
    }
    Glamd = -(fDel - fStar);
  }

  // ---- 7. a1, b1, al1b. -------------------------------------------------------------------
  const T a1   = psi1.dot(Gud*phi1) / m_branchTangentDel;
  const T b1   = psi1.dot(Gud*phi0 + Glamd) / m_branchTangentDel;
  const T al1b = -(a1*al1/lam0 + (T)2*b1);

  // ---- 8. Degeneracy of al1b ("== 0" of swibra.m:66, made numerical). A RELATIVE
  // threshold, not exact equality: a1 and b1 carry an O(del) relative FD error, so a
  // cancellation below that is not resolvable, and the default BranchTangentTol = 1e-2 =
  // 10*del keeps the verdict above its own resolution.
  const T scaleABE = math::abs(a1*al1/lam0) + (T)2*math::abs(b1);
  if (!(math::abs(al1b) > m_branchTangentTol * scaleABE))
  {
    if (m_verbose)
      gsInfo<<"computeBranchTangent: al1="<<al1<<", a1="<<a1<<", b1="<<b1<<", al1b="<<al1b
            <<" [Degenerate: |al1b| not above BranchTangentTol*scaleABE = "
            <<m_branchTangentTol*scaleABE<<" -- no distinct branch to switch to].\n";
    return branchTangent::Degenerate;
  }

  // ---- 9. tau1 = al1b*phi1 + a1*phi0 (displacement), a1 (load); normalise. --------------
  tau1U = al1b*phi1 + a1*phi0;
  tau1L = a1;
  const T tau1nrm = safeDistance(tau1U, tau1L);
  if (tau1nrm > (T)0 && math::isfinite(tau1nrm)) { tau1U /= tau1nrm; tau1L /= tau1nrm; }

  // ---- 10. Verbose ABE line (mirrors swibra.m:65), mandatory observable of this method. -
  if (m_verbose)
    gsInfo<<"computeBranchTangent: al1="<<al1<<", a1="<<a1<<", b1="<<b1<<", al1b="<<al1b
          <<" [Success]\n";

  return branchTangent::Success;
}

// ------------------------------------------------------------------------------------------------------------
// ---------------------------------------Output functions-----------------------------------------------------
// ------------------------------------------------------------------------------------------------------------

template <class T>
void gsALMBase<T>::_initOutputExtended()
{
  gsInfo<<"\t";
  gsInfo<<std::setw(4)<<std::left<<"It.";
  gsInfo<<std::setw(17)<<std::left<<"Res. F";
  gsInfo<<std::setw(17)<<std::left<<"|dU|/|Du|";
  gsInfo<<std::setw(17)<<std::left<<"dL/DL";
  gsInfo<<std::setw(17)<<std::left<<"K_T * φ";
  gsInfo<<std::setw(17)<<std::left<<"|U|";
  gsInfo<<std::setw(17)<<std::left<<"|φ|";
  gsInfo<<std::setw(17)<<std::left<<"L";
  gsInfo<<std::setw(17)<<std::left<<"|DU|";
  gsInfo<<std::setw(17)<<std::left<<"|Dφ|";
  gsInfo<<std::setw(17)<<std::left<<"DL";
  gsInfo<<std::setw(17)<<std::left<<"|dU|";
  gsInfo<<std::setw(17)<<std::left<<"|dφ|";
  gsInfo<<std::setw(17)<<std::left<<"dL";
  gsInfo<<std::setw(17)<<std::left<<"Dmin";
  gsInfo<<std::setw(17)<<std::left<<"note";
  gsInfo<<"\n";

  m_note = "";
}

template <class T>
void gsALMBase<T>::_stepOutputExtended()
{
  gsInfo<<"\t";
  gsInfo<<std::setw(4)<<std::left<<m_numIterations;
  gsInfo<<std::setw(17)<<std::left<<m_residueF;
  gsInfo<<std::setw(17)<<std::left<<m_residueU;
  gsInfo<<std::setw(17)<<std::left<<m_residueL;
  gsInfo<<std::setw(17)<<std::left<<m_residueKTPhi;
  gsInfo<<std::setw(17)<<std::left<<(m_U).norm();
  gsInfo<<std::setw(17)<<std::left<<(m_V).norm();
  gsInfo<<std::setw(17)<<std::left<<(m_L);
  gsInfo<<std::setw(17)<<std::left<<m_DeltaU.norm();
  gsInfo<<std::setw(17)<<std::left<<m_DeltaV.norm();
  gsInfo<<std::setw(17)<<std::left<<m_DeltaL;
  gsInfo<<std::setw(17)<<std::left<<m_deltaU.norm();
  gsInfo<<std::setw(17)<<std::left<<m_deltaV.norm();
  gsInfo<<std::setw(17)<<std::left<<m_deltaL;
  gsInfo<<std::setw(17)<<std::left<<_bisectionTerminationFunction(m_U,false);
  gsInfo<<std::setw(17)<<std::left<<m_note;
  gsInfo<<"\n";

  m_note = "";
}

} // namespace gismo
