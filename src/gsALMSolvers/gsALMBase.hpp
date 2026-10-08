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
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMHelper.h>

namespace gismo
{

template <class T>
void gsALMBase<T>::defaultOptions()
{
    m_options.addInt ("MaxIter","Maximum number of corrector iterations per step",100);
    m_options.addReal("Tol","Tolerance",1e-6);
    m_options.addReal("TolF","Tolerance",1e-3);
    m_options.addReal("TolU","Tolerance",1e-6);
    m_options.addReal("ConstraintTol","Tolerance on the normalised path-constraint residual at acceptance (|g|/s^2; load control: |DL - DL_target|/|DL_target|)",1e-6);
    m_options.addReal("Perturbation","Set Perturbation factor Tau",1e3);

    m_options.addReal("Length","Arclength",1e-2);

    m_options.addSwitch("AdaptiveLength","Adaptive length",false);
    m_options.addInt ("AdaptiveIterations","Desired iterations for adaptive length",10);

    m_options.addSwitch("Quasi","Use Quasi Newton method",false);
    m_options.addInt ("QuasiIterations","Number of iterations for quasi newton method",-1);

    m_options.addInt ("BifurcationMethod","Bifurcation Identification based on: -1: nothing, 0: Determinant;  1: Eigenvalue",bifmethod::Eigenvalue);

    m_options.addInt ("SingularPointFailure","What to do when a singular point determination fails?: 0 = Apply solution anyways; 1 = Proceed without singular point",SPfail::With);
    m_options.addReal("SingularPointTestTol", "Detection tolerance for singular points. Product of the mode shape and the forcing should be below this tolerance", 1e-6);
    m_options.addInt ("SingularPointTestIt" , "Number of iterations of the power method for the singular point test",5);

    m_options.addReal("SingularPointComputeTolE", "Tolerance for the extended iterations to compute a bifurcation point", 1e-10);
    m_options.addReal("SingularPointComputeTolB", "Tolerance for the bisection iterations to compute a bifurcation point. If tol = 0, no bi-section method is used.", 0);

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
    m_toleranceC          = m_options.getReal("ConstraintTol");

    m_tau                 = m_options.getReal("Perturbation");

    m_adaptiveLength      = m_options.getSwitch("AdaptiveLength");
    m_desiredIterations   = m_options.getInt ("AdaptiveIterations");

    m_quasiNewton         = m_options.getSwitch("Quasi");
    m_quasiNewtonInterval = m_options.getInt ("QuasiIterations");

    m_bifurcationMethod   = m_options.getInt ("BifurcationMethod");
    if (m_bifurcationMethod!=bifmethod::Nothing && m_bifurcationMethod!=bifmethod::Determinant && m_bifurcationMethod!=bifmethod::Eigenvalue)
        GISMO_ERROR("Unknown BifurcationMethod "<<m_bifurcationMethod<<"; use -1 (nothing), 0 (determinant) or 1 (eigenvalue).");
    m_solver = gsSparseSolver<T>::get( m_options.getString("Solver") );
    if  (!dynamic_cast<typename gsSparseSolver<T>::SimplicialLDLT*>(m_solver.get()) && m_bifurcationMethod==bifmethod::Determinant)
    {
        gsWarn<<"Determinant method cannot be used with solvers other than LDLT. Bifurcation method will be set to 'Eigenvalue'.\n";
        m_bifurcationMethod = bifmethod::Eigenvalue;
    }

    m_verbose             = m_options.getSwitch ("Verbose");
    m_relax               = m_options.getReal("Relaxation");

    m_arcLength = m_arcLength_prev = m_arcLength_ori = m_options.getReal("Length");

    m_SPfail = m_options.getInt ("SingularPointFailure");
    m_SPTestTol = m_options.getReal("SingularPointTestTol");
    m_SPTestIt  = m_options.getInt ("SingularPointTestIt");

    m_SPCompTolE = m_options.getReal("SingularPointComputeTolE");
    m_SPCompTolB = m_options.getReal("SingularPointComputeTolB");

    this->_configureDriver();
}

template <class T>
void gsALMBase<T>::init(bool stability)
{
  if (stability && !m_constraint)
  {
    gsWarn<<"gsALMBase::initialize: no path constraint is installed; the stability is not computed.\n";
    m_stability = m_stabilityPrev = 0;
  }
  else if (stability)
  {
    this->_computeStability(m_U);
    m_stability = this->stability(); // requires Jabobian!!
    m_stabilityPrev = m_stability;
  }
  else
    m_stability = m_stabilityPrev = 0;
}

template <class T>
void gsALMBase<T>::computeLength()
{
  // Adaptive arc length: dl_{n+1} = dl_n * clamp(N_desired / N_n, 1/2, 2), N_n the
  // corrector iterations of the last converged step. N_n == 0 (no step yet, e.g.
  // resetLength() on a fresh solver) carries no information: keep the length.
  T fac = (m_numIterations > 0)
        ? static_cast<T>(m_desiredIterations) / static_cast<T>(m_numIterations)
        : T(1);
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
  m_driver->controller().reduceLength(fac);
  if (m_constraint)
    m_arcLength = m_driver->controller().length();
  return m_arcLength;
}

template <class T>
T gsALMBase<T>::resetLength()
{
  m_arcLength = m_arcLength_prev = m_arcLength_ori;
  if (m_adaptiveLength)
    computeLength();
  m_driver->controller().resetLength(m_numIterations);
  if (m_constraint)
  {
    m_arcLength      = m_driver->controller().length();
    m_arcLength_prev = m_driver->controller().previousLength();
    this->_syncHistoryStep();
  }
  return m_arcLength;
}

template <class T>
T gsALMBase<T>::getLength()
{
  return m_constraint ? m_driver->controller().length() : m_arcLength;
}

template <class T>
void gsALMBase<T>::setLength(T length)
{
  m_options.setReal("Length",length);
  m_arcLength = m_arcLength_prev = m_arcLength_ori = m_options.getReal("Length");
  m_driver->controller().setLength(m_arcLength);
  this->_syncHistoryStep();
}

template <class T>
void gsALMBase<T>::setLength(T length, bool adaptive)
{
  this->setLength(length);
  m_adaptiveLength = adaptive;
  m_desiredIterations = 10; // number of desired iterations defaults to 10
  m_driver->controller().setAdaptive(m_adaptiveLength, m_desiredIterations);
}

template <class T>
void gsALMBase<T>::setLength(T length, bool adaptive, index_t iterations)
{
  this->setLength(length);
  m_adaptiveLength = adaptive;
  m_desiredIterations = iterations;
  m_driver->controller().setAdaptive(m_adaptiveLength, m_desiredIterations);
}

template <class T>
void gsALMBase<T>::setLength(T length, index_t iterations)
{
  this->setLength(length);
  m_adaptiveLength = true;
  m_desiredIterations = iterations;
  m_driver->controller().setAdaptive(m_adaptiveLength, m_desiredIterations);
}

template <class T>
gsVector<T> gsALMBase<T>::computeResidual(const gsVector<T> & U, const T & L)
{
  gsVector<T> resVec;
  if (m_residualFun)
  {
    if (!m_residualFun(U, L, resVec))
      throw 2;
  }
  else if (!this->_evaluateAt(U, L, &resVec, NULL))
    throw 2;
  return resVec;
}

template <class T>
void gsALMBase<T>::computeResidual()
{
  m_resVec = this->computeResidual(m_U + m_DeltaU, m_L + m_DeltaL);
}

template <class T>
void gsALMBase<T>::computeResidualNorms()
{
  if (m_numIterations ==0 ) // then define the residual
  {
    // m_basisResidualF = m_resVec.norm();
    m_basisResidualF = ((m_L+m_DeltaL) * m_forcing).norm();
    m_basisResidualU = m_DeltaU.norm();

    // m_residueF = m_resVec.norm() / (((m_L+m_DeltaL) * m_forcing).norm());
    m_residueF = m_resVec.norm() / (((m_L+m_DeltaL) * m_forcing).norm());
    m_residueU = m_deltaU.norm()/m_DeltaU.norm();
    // m_residueU = m_basisResidualU;
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
  return this->_computeJacobianAt(U, deltaU, m_L + m_DeltaL);
}

template <class T>
gsSparseMatrix<T> gsALMBase<T>::_computeJacobianAt(const gsVector<T> & U, const gsVector<T> & deltaU, T L)
{
  // Compute Jacobian
  gsSparseMatrix<T> m;
  m_note += "J";
  const bool ok = m_djacobian ? m_djacobian(U,deltaU,m)
                              : this->_evaluateAt(U, L, NULL, &m);
  if (!ok)
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
  m_deltaUt = this->solveSystem(m_forcing);
}

template <class T>
gsStatus gsALMBase<T>::step()
{
  return this->_stepDriver();
}

// ------------------------------------------------------------------------------------------------------------
// ---------------------------------------Driver mode----------------------------------------------------------
// ------------------------------------------------------------------------------------------------------------

template <class T>
void gsALMBase<T>::_constructOwned(const gsVector<T> & Force)
{
  m_ownSpace.reset(new gsVectorStateSpace<T>(Force.size()));
  m_ownProblem.F = Force;

  // The legacy ALResidual_t is F_int - lambda*F; the driver needs
  // R = lambda*F - F_int = -ALResidual. The negation is exact.
  m_ownProblem.residual = [this](T l, gsVector<T> & R, gsEvalMode /*mode*/) -> bool
  {
    if (!m_residualFun(m_ownSpace->working(), l, R))
      return false;
    R = -R;
    return true;
  };

  // The tangent is passed through. The step dx handed to the legacy
  // dJacobian_t is the latest applied corrector update (zero right after the
  // predictor).
  m_ownProblem.jacobian = [this](T /*lambda*/, gsSparseMatrix<T> & K) -> bool
  {
    return m_djacobian(m_ownSpace->working(), m_driver->context().deltaU, K);
  };

  m_problem = &m_ownProblem;
  m_space   = m_ownSpace.get();
}

template <class T>
void gsALMBase<T>::_initMembers(bool readOptions)
{
  m_driver.reset(new gsContinuationDriver<T>(*m_problem, *m_space));
  m_driver->addObserver([this](const gsStepEvent<T> & ev) { this->_observe(ev); });

  m_constraint.reset();
  m_driverDirty = true;
  m_inBisection = false;
  m_prevKnown   = false;
  m_numAccepted = 0;
  m_initialized = false;
  m_verbose     = false;
  m_numDof      = m_space->tangentSize();

  // initialize variables
  m_numIterations = 0;
  this->defaultOptions();
  if (readOptions)
    gsALMBase<T>::getOptions();
  this->setLength(1e-2);
  m_converged = false;
  m_indicator = 0;
  m_negatives = 0;
  m_stability = 0;
  m_stabilityPrev = 0;

  // initialize errors
  m_basisResidualF = 0.0;
  m_basisResidualU = 0.0;

  m_status = gsStatus::NotStarted;

  // zero state, so that every accessor is defined before initialize()
  this->_refreshMirror();
  m_L = m_Lprev = m_DeltaL = m_deltaL = 0;
  m_Uprev = m_U;
  m_DeltaU = m_deltaU = gsVector<T>::Zero(m_numDof);
  m_DeltaUold = gsVector<T>::Zero(m_numDof);
  m_DeltaLold = 0;
}

template <class T>
void gsALMBase<T>::_configureDriver()
{
  m_driver->setTest(gsConvergenceTest<T>::makeArcLengthGate(m_toleranceF, m_toleranceU, m_toleranceC, m_maxIterations));
  m_driver->controller().setLength(m_arcLength);
  m_driver->controller().setAdaptive(m_adaptiveLength, m_desiredIterations);
  m_driver->corrector().setQuasi(m_quasiNewton, m_quasiNewtonInterval);
  m_driver->corrector().setRelaxation(T(1));
  m_driver->linearSolver().setSolver(m_options.getString("Solver"));
  m_driverDirty = true;
  this->_syncHistoryStep();
}

template <class T>
void gsALMBase<T>::setConstraint(const typename gsPathConstraint<T>::Ptr & c)
{
  GISMO_ENSURE(c, "gsALMBase::setConstraint: null constraint");
  m_constraint = c;
  m_driver->setConstraint(c);
  m_driverDirty = true;
}

template <class T>
void gsALMBase<T>::initialize(bool stability)
{
  m_initialized = true;
  this->initMethods();
  m_numAccepted = 0;
  if (m_constraint)
    this->_initializeDriver();
  this->init(stability);
}

template <class T>
void gsALMBase<T>::_initializeDriver()
{
  m_numDof = m_space->tangentSize();
  if (m_ownSpace)
    m_ownSpace->importCommitted(gsVector<T>::Zero(m_numDof));
  m_driver->setLambda(T(0));
  m_driver->resetHistory();
  m_prevKnown = false;

  this->_refreshMirror();
  m_L     = 0;
  m_Uprev = m_U;
  m_Lprev = 0;
  m_DeltaU = gsVector<T>::Zero(m_numDof);
  m_deltaU = gsVector<T>::Zero(m_numDof);
  m_DeltaL = 0;
}

// Cost: one O(n) export per accepted step (a copy for a vector space).
template <class T>
void gsALMBase<T>::_refreshMirror()
{
  if (!(m_space->copyable() && m_space->exportWorking(m_U)))
    m_U.resize(0);
}

template <class T>
bool gsALMBase<T>::_vectorLike() const
{
  if (!m_space->copyable())
    return false;
  gsVector<T> u;
  return m_space->exportWorking(u) && u.size() == m_space->tangentSize();
}

template <class T>
bool gsALMBase<T>::_importMirror()
{
  if (!m_space->importCommitted(m_U))
    return false;
  m_driver->setLambda(m_L);
  this->_setHistoryFromMirror();
  return true;
}

template <class T>
void gsALMBase<T>::_setHistory(const gsVector<T> & DU, T DL)
{
  gsStepHistory<T> h;
  h.DeltaU     = DU;
  h.DeltaL     = DL;
  h.stepSize   = m_driver->controller().previousLength();
  h.iterations = 0;
  h.valid      = true;
  m_driver->setHistory(h);
}

template <class T>
void gsALMBase<T>::_setHistoryFromMirror()
{
  if (m_prevKnown && m_U.size() != 0 && m_U.size() == m_Uprev.size())
  {
    const gsVector<T> DU = m_U - m_Uprev;
    this->_setHistory(DU, m_L - m_Lprev);
  }
  else
    m_driver->resetHistory();
}

template <class T>
void gsALMBase<T>::_syncHistoryStep()
{
  if (m_driver->history().valid)
  {
    gsStepHistory<T> h = m_driver->history();
    h.stepSize = m_driver->controller().previousLength();
    m_driver->setHistory(h);
  }
}

template <class T>
void gsALMBase<T>::setSolution(const gsVector<T> & U, const T & L)
{
  if (!m_constraint)
  {
    m_L = L; m_U = U;
    return;
  }

  const bool vectorLike = this->_vectorLike() && U.size() == m_space->tangentSize();
  GISMO_ENSURE(vectorLike, "gsALMBase::setSolution needs a copyable state space and a vector of the tangent size.");
  const bool imported = m_space->importCommitted(U);
  GISMO_ENSURE(imported, "gsALMBase::setSolution: the state space rejected the solution.");
  m_driver->setLambda(L);
  m_U = U;
  m_L = L;
  this->_setHistoryFromMirror();

  if (m_initialized && !m_inBisection)
    this->_stabilityAtSolution();
}

/*
  One stability evaluation at the committed point (tangent at m_L), made when
  the solution is set from outside. The stored stability of a step is that of
  the point it starts from, so a bifurcation in the first step after the
  external set is detected. If the evaluation fails the stability is unknown (0).
  The tangent is the driver's, assembled and factorized once by
  _stabilityAtCommitted(); the predictor of the next step reuses it when its
  constraint needs the load tangent. Cost: that of _stabilityAtCommitted().
*/
template <class T>
void gsALMBase<T>::_stabilityAtSolution()
{
  const gsStatus status = m_status;
  bool ok = true;
  try
  {
    if (m_driverDirty)
    {
      m_driver->setup();
      m_driverDirty = false;
    }
  }
  catch (...)
  {
    ok = false;
  }
  if (ok)
    ok = this->_stabilityAtCommitted();

  if (ok)
    m_stabilityPrev = m_stability;
  else
    m_stability = m_stabilityPrev = 0;
  m_status = status;
}

template <class T>
void gsALMBase<T>::setPrevious(const gsVector<T> & Uprev, const T & Lprev)
{
  m_Uprev = Uprev;
  m_Lprev = Lprev;
  if (!m_constraint || m_U.size() == m_Uprev.size())
  {
    m_DeltaUold = m_U - m_Uprev;
    m_DeltaLold = m_L - m_Lprev;
  }
  if (m_constraint)
  {
    m_prevKnown = true;
    this->_setHistoryFromMirror();
  }
}

template <class T>
void gsALMBase<T>::setSolutionStep(const gsVector<T> & DU, const T & DL)
{
  m_DeltaUold = DU;
  m_DeltaLold = DL;
  if (m_constraint)
    this->_setHistory(DU, DL);
}

template <class T>
void gsALMBase<T>::resetStep()
{
  m_DeltaUold.setZero();
  m_DeltaLold = 0;
  if (m_constraint)
    m_driver->resetHistory();
}

template <class T>
void gsALMBase<T>::applyInitialGuess(const gsVector<T> & Uguess, T Lguess)
{
  if (Uguess.size() != m_U.size())
    return;
  const gsVector<T> DU = Uguess - m_U;
  this->_setHistory(DU, Lguess - m_L);
}

template <class T>
T gsALMBase<T>::distance(const gsVector<T> & DeltaU, const T DeltaL) const
{
  if (m_constraint)
    return m_constraint->distance(DeltaU, DeltaL);
  GISMO_NO_IMPLEMENTATION;
}

/*
  Evaluates the residual (with the sign of the legacy ALResidual_t, i.e. the
  negative of the driver residual) or the tangent at the state U through the
  state space: the working state is saved, replaced by U, restored afterwards.
  Cost: three O(n) vector copies plus the evaluation. Used off the step path
  only (stability, singular points) and only if the solver was built on a
  caller's problem; the legacy operators are evaluated directly otherwise.
  A state space that is not copyable cannot be set to U: the tangent is then
  evaluated at the committed state and the residual is unavailable.
*/
template <class T>
bool gsALMBase<T>::_evaluateAt(const gsVector<T> & U, T L, gsVector<T> * R, gsSparseMatrix<T> * K)
{
  if (!m_space->rollback())
    return false;
  const bool copyable = m_space->copyable();
  if (R && !copyable)
    return false;

  gsVector<T> saved;
  bool ok = true;
  if (copyable)
    ok = m_space->exportWorking(saved) && m_space->importCommitted(U);

  if (ok)
  {
    if (R)
    {
      ok = m_problem->residual(L, *R, gsEvalMode::Probe) && R->size() == m_space->tangentSize();
      if (ok)
        *R = -*R;
    }
    else
      ok = m_problem->jacobian(L, *K);
  }

  if (copyable && saved.size() != 0)
    m_space->importCommitted(saved);
  m_driver->invalidateCommittedTangent();
  return ok;
}

template <class T>
gsStatus gsALMBase<T>::_stepDriver()
{
  if (!m_initialized)
  {
    m_note = "Arc-Length Method is not initialized! Call initialize()";
    m_status = gsStatus::OtherError;
    return m_status;
  }

  m_stabilityPrev = m_stability;
  if (!m_constraint)
  {
    m_note = "No path constraint is installed";
    m_status = gsStatus::OtherError;
    return m_status;
  }
  m_note.clear();

  gsStatus st = gsStatus::OtherError;
  bool failedByException = false;
  try
  {
    if (m_driverDirty)
    {
      m_driver->setup();
      m_driverDirty = false;
    }

    if (m_Uguess.rows()!=0 && m_Uguess.cols()!=0 && m_Uguess.size()==m_U.size()
        && (m_Uguess-m_U).norm()!=0 && (m_Lguess-m_L)!=0)
      this->applyInitialGuess(m_Uguess, m_Lguess);

    st = m_driver->step();
  }
  catch (...)
  {
    st = gsStatus::OtherError;
    failedByException = true;
    m_note = "exception in the problem operators or in the driver setup";
  }

  const gsStepContext<T> & ctx = m_driver->context();
  m_numIterations = m_driver->numIterations();
  m_DeltaU = ctx.DeltaU;
  m_DeltaL = ctx.DeltaL;
  m_deltaU = ctx.deltaU;
  m_converged = (st == gsStatus::Success);
  // The tangent of the last attempt, which is the failed-step tangent when the
  // step failed; after an accepted step the stability evaluation replaces it by
  // the tangent at the committed point (with the method Nothing it stays the
  // last assembled tangent). Cost: one O(nnz) copy per step.
  if (ctx.K.rows() == m_numDof)
    m_jacMat = ctx.K;
  if (!failedByException && st != gsStatus::Success)
    m_note = m_driver->report().note;

  if (st == gsStatus::Success)
  {
    m_Uprev = m_U;
    m_Lprev = m_L;
    this->_refreshMirror();
    m_L = m_driver->lambda();

    m_arcLength      = m_driver->controller().length();
    m_arcLength_prev = m_driver->controller().previousLength();
    m_prevKnown = true;

    // The stability is evaluated once per accepted step, on the tangent at the
    // committed point; the cost is one Jacobian assembly and factorization plus
    // one stability analysis (O(n^3) for a dense eigenvalue solve).
    this->_stabilityAfterStep();

    if (!m_inBisection)
    {
      if (m_onAccepted)
        m_onAccepted(m_numAccepted, m_L);
      ++m_numAccepted;
    }
  }

  m_status = st;
  return m_status;
}

/*
  The stability of an accepted point is that of the tangent there, as in
  pde2path (cont.m:77) and BifurcationKit: the driver assembles and factorizes K
  at the committed point (gsContinuationDriver::assembleCommittedTangent), the
  tangent is copied to m_jacMat for the singular-point tests, and the Determinant
  method classifies the LDL^T pivots of that very factorization, so it needs no
  second one. The Eigenvalue method (and a Determinant method without pivots)
  classifies m_jacMat through computeStability(false). Cost: one Jacobian
  assembly, one factorization, plus the analysis; the next predictor reuses the
  factorization when its constraint needs the load tangent.
*/
template <class T>
bool gsALMBase<T>::_stabilityAtCommitted()
{
  bool ok;
  try
  {
    ok = m_driver->assembleCommittedTangent();
  }
  catch (...)
  {
    ok = false;
  }
  if (!ok)
    return false;

  m_jacMat = m_driver->context().K;

  if (m_bifurcationMethod == bifmethod::Determinant && m_driver->hasCommittedTangent()
      && m_driver->linearSolver().pivots(m_stabilityVec))
  {
    m_stabilityAnalysis.setMethod(m_bifurcationMethod);
    if (m_stabilityAnalysis.classify(m_stabilityVec))
    {
      m_negatives = m_stabilityAnalysis.negatives();
      m_indicator = m_stabilityAnalysis.indicator();
      m_stability = m_stabilityAnalysis.stability();
      return true;
    }
  }

  return this->computeStability(false) == gsStatus::Success;
}

template <class T>
void gsALMBase<T>::_stabilityAfterStep()
{
  const bool ok = (m_bifurcationMethod == bifmethod::Nothing)
                ? this->computeStability(false) == gsStatus::Success
                : this->_stabilityAtCommitted();
  if (!ok)
    m_note = "stability evaluation failed; stability() keeps its previous value";
}

template <class T>
void gsALMBase<T>::_observe(const gsStepEvent<T> & ev) const
{
  if (!m_verbose)
    return;

  const gsStepContext<T> & c = ev.ctx;
  switch (ev.kind)
  {
    case gsStepEvent<T>::StepStarted:
      gsInfo<<"\t";
      gsInfo<<std::setw(4)<<std::left<<"It.";
      gsInfo<<std::setw(17)<<std::left<<"|R|";
      gsInfo<<std::setw(17)<<std::left<<"|dU|";
      gsInfo<<std::setw(17)<<std::left<<"|DU|";
      gsInfo<<std::setw(17)<<std::left<<"dL";
      gsInfo<<std::setw(17)<<std::left<<"DL";
      gsInfo<<std::setw(17)<<std::left<<"L+DL";
      gsInfo<<std::setw(17)<<std::left<<"g/s^2";
      gsInfo<<"\n";
      break;
    case gsStepEvent<T>::IterationDone:
      gsInfo<<"\t";
      gsInfo<<std::setw(4)<<std::left<<c.iteration;
      gsInfo<<std::setw(17)<<std::left<<c.R.norm();
      gsInfo<<std::setw(17)<<std::left<<c.deltaU.norm();
      gsInfo<<std::setw(17)<<std::left<<c.DeltaU.norm();
      gsInfo<<std::setw(17)<<std::left<<c.deltaL;
      gsInfo<<std::setw(17)<<std::left<<c.DeltaL;
      gsInfo<<std::setw(17)<<std::left<<(c.lambda + c.DeltaL);
      gsInfo<<std::setw(17)<<std::left<<(c.constraintResidual / (c.stepSize*c.stepSize));
      gsInfo<<"\n";
      break;
    case gsStepEvent<T>::StepAccepted:
      gsInfo<<"\t step accepted after "<<c.iteration<<" iterations; lambda = "<<c.lambda<<"\n";
      break;
    case gsStepEvent<T>::StepRejected:
      gsInfo<<"\t step rejected after "<<c.iteration<<" iterations\n";
      break;
  }
}

// ------------------------------------------------------------------------------------------------------------
// ---------------------------------------Singular point methods-----------------------------------------------
// ------------------------------------------------------------------------------------------------------------
template <class T>
void gsALMBase<T>::_computeSingularPoint(const gsVector<T> & U, const T & L, bool switchBranch, bool jacobian, bool testPoint)
{
  // Determines if the point is a bifurcation point. If not, it assumes it is
  bool test = (testPoint) ? this->_testSingularPoint(jacobian) : true;

  m_U = U;
  m_L = L;
  if (test)
  {
    bool converged = false;
    // First stage: bisection method
    if (m_SPCompTolB != 0)
      this->_bisectionSolve(m_U,m_L,m_SPCompTolB);

    converged = this->_extendedSystemSolve(m_U, m_L, m_SPCompTolE);

    if (switchBranch && (converged || m_SPfail==1))
    {
      this->switchBranch(); // leaves the stability unknown
      return;
    }

    // The stability of the singular point is taken equal to the one of the
    // previous point, so that the singular point is not detected again.
    m_stability = m_stabilityPrev;
  }
}

// tolB and switchBranch will be defaulted
template <class T>
gsStatus gsALMBase<T>::computeSingularPoint(const gsVector<T> & U, const T & L, bool switchBranch, bool jacobian, bool testPoint)
{
  if (m_constraint && !this->_vectorLike())
  {
    gsWarn<<"computeSingularPoint needs a copyable state space that exports a vector of the tangent size.\n";
    m_status = gsStatus::OtherError;
    return m_status;
  }

  try
  {
    this->_computeSingularPoint(U,L,switchBranch,jacobian, testPoint);
  }
  catch (int errorCode)
  {
    if      (errorCode==1)
      m_status = gsStatus::NotConverged;
    else if (errorCode==2)
      m_status = gsStatus::AssemblyError;
    else if (errorCode==3)
      m_status = gsStatus::SolverError;
    else
      m_status = gsStatus::OtherError;
  }
  catch (...)
  {
    m_status = gsStatus::OtherError;
  }

  if (m_constraint && !this->_importMirror())
    m_status = gsStatus::OtherError;
  return m_status;
}

template <class T>
bool gsALMBase<T>::_testSingularPoint(bool jacobian)
{
  // First, approximate the eigenvector of the Jacobian by a few arc length iterations
  // Initiate m_V and m_DeltaVDET

  if (jacobian)
  {
    m_jacMat = this->_computeJacobianAt(m_U, m_deltaU, m_L);
    this->factorizeMatrix(m_jacMat);
  }

  m_V = gsVector<T>::Ones(m_numDof);
  m_V.normalize();
  factorizeMatrix(m_jacMat);

  for (index_t k = 0; k<m_SPTestIt; k++)
  {
    m_V = this->solveSystem(m_V);
    m_V.normalize();
  }

  T dot = (abs(m_V.dot(m_forcing)));
  if ( (abs(dot) / m_SPTestTol > 1e-1) && (abs(dot) / m_SPTestTol < 10) )
  {
    gsInfo<<"Warning: the singular point test is close to its tolerance. dot/tol = "<<abs(dot)/m_SPTestTol<<" dot = "<<dot<<"\t tolerance = "<<m_SPTestTol<<"\n";
  }
  if (dot < m_SPTestTol)    // Bifurcation point
  {
    if (m_verbose) {gsInfo<<"\t Bifurcation point\n";}
    return true;
  }
  else        // Limit point
  {
    if (m_verbose) {gsInfo<<"\t Limit point\n";}
    return false;
  }
}

template <class T>
bool gsALMBase<T>::isBifurcation(bool jacobian)
{
  if (m_constraint && !this->_vectorLike())
  {
    gsWarn<<"isBifurcation needs a copyable state space that exports a vector of the tangent size.\n";
    return false;
  }
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
    if      (errorCode==1)
      m_status = gsStatus::NotConverged;
    else if (errorCode==2)
      m_status = gsStatus::AssemblyError;
    else if (errorCode==3)
      m_status = gsStatus::SolverError;
    else
      m_status = gsStatus::OtherError;
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
    m_jacMat = this->_computeJacobianAt(x, dx, m_L);
    this->factorizeMatrix(m_jacMat);
  } // otherwise m_jacMat holds the stored tangent (see step())

  m_stabilityAnalysis.setMethod(m_bifurcationMethod);

  bool ok = false;
  if (m_bifurcationMethod == bifmethod::Determinant)
  {
    if ( auto * s = dynamic_cast<typename gsSparseSolver<T>::SimplicialLDLT*>(m_solver.get()) )
    {
      factorizeMatrix(m_jacMat);
      m_stabilityVec = s->vectorD();
      ok = m_stabilityAnalysis.classify(m_stabilityVec);
    }
    else
    {
      gsWarn<<"Determinant stability method only works with SimplicialLDLT solver, current solver is "<<m_options.getString("Solver")<<"\n";
      throw 3;
    }
  }
  else if (m_bifurcationMethod == bifmethod::Eigenvalue)
    ok = m_stabilityAnalysis.compute(m_jacMat, shift, m_stabilityVec);
  else if (m_bifurcationMethod == bifmethod::Nothing)
  {
    m_stabilityVec.setZero(m_numDof);
    ok = m_stabilityAnalysis.classify(m_stabilityVec);
  }

  if (!ok)
    throw 3;

  m_negatives = m_stabilityAnalysis.negatives();
  m_indicator = m_stabilityAnalysis.indicator();
  m_stability = m_stabilityAnalysis.stability();
}

template <class T>
index_t gsALMBase<T>::stability() const
{
  return gsStabilityAnalysis<T>::sign(m_indicator);
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

  m_jacMat = this->_computeJacobianAt(m_U, m_deltaU, m_L); // Jacobian evaluated on (m_U, m_L)
  this->factorizeMatrix(m_jacMat);
  m_basisResidualKTPhi = (m_jacMat*m_V).norm();

  m_DeltaV = gsVector<T>::Zero(m_numDof);
  m_DeltaU.setZero();
  m_DeltaL = 0.0;

  m_deltaV = gsVector<T>::Zero(m_numDof);
  m_deltaU.setZero();
  m_deltaL = 0.0;
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
    m_jacMat = this->computeJacobian(m_U+m_DeltaU,m_deltaU);
    this->factorizeMatrix(m_jacMat);
    m_residueKTPhi = (m_jacMat*(m_V+m_DeltaV)).norm(); // /m_basisResidualKTPhi;
    m_resVec = this->computeResidual(m_U+m_DeltaU,m_L+m_DeltaL);
    computeResidualNorms();
    if (m_verbose)
      _stepOutputExtended();

    // termination criteria
    // if ( m_residueF < m_toleranceF && m_residueU < m_toleranceU &&  m_residueKTPhi< tol )
    if ( m_residueKTPhi< tol )
    {
      m_U += m_DeltaU;
      m_L += m_DeltaL;
      gsInfo<<"Iterations finished. U.norm() = "<<m_U.norm()<<"\t L = "<<m_L<<"\n";
      break;
    }
    if (m_numIterations == m_maxIterations-1)
    {
      gsInfo<<"Warning: Extended iterations did not converge! \n";
      if (m_SPfail==1)
      {
        m_U += m_DeltaU;
        m_L += m_DeltaL;
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
  m_resVec = this->computeResidual(m_U+m_DeltaU, m_L+m_DeltaL);
  m_jacMat = this->computeJacobian(m_U+m_DeltaU,m_deltaU);
  this->factorizeMatrix(m_jacMat);

  // m_jacMat = m_jacobian(m_U);

  m_deltaUt = this->solveSystem(m_forcing); // DeltaV1
  m_deltaUbar = this->solveSystem(-m_resVec); // DeltaV2

  real_t eps = 1e-8;
  gsSparseMatrix<T> jacMatEps = this->computeJacobian((m_U+m_DeltaU) + eps*(m_V+m_DeltaV));
  m_note += "J"; // mark jacobian computation
  gsVector<T> h1 = 1/eps * ( jacMatEps * m_deltaUt ) - 1/eps * m_forcing;
  gsVector<T> h2 = m_jacMat * (m_V+m_DeltaV) + 1/eps * ( jacMatEps * m_deltaUbar + m_resVec );

  this->factorizeMatrix(m_jacMat);

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
bool gsALMBase<T>::_bisectionSolve(const gsVector<T> & U, const T L, const T tol)
{
  m_U = U;
  m_L = L;

  gsVector<T> U_old  = m_U;
  T L_old = L;

  // Stores the arc length and the controller; the trial steps neither adapt the
  // length nor call the accepted-step callback. The state is restored on every
  // exit, also when a helper throws (computeSingularPoint catches it).
  struct Restore
  {
    gsALMBase<T> &            s;
    const T                   arcLength;
    const bool                adaptive;
    const gsStepController<T> controller;
    bool                      done;

    explicit Restore(gsALMBase<T> & a)
    : s(a), arcLength(a.m_arcLength), adaptive(a.m_adaptiveLength),
      controller(a.m_driver->controller()), done(false)
    {
      s.m_inBisection = true;
      s.m_adaptiveLength = false;
      if (s.m_constraint)
        s.m_driver->controller().setAdaptive(false, s.m_desiredIterations);
    }

    void release()
    {
      if (done)
        return;
      done = true;
      s.m_inBisection = false;
      s.m_arcLength = arcLength;
      s.m_adaptiveLength = adaptive;
      if (s.m_constraint)
      {
        s.m_driver->controller() = controller;
        s.m_arcLength_prev = s.m_driver->controller().previousLength();
      }
    }

    ~Restore() { release(); }
  } restore(*this);
  bool converged = false;

  T referenceError = _bisectionTerminationFunction(m_U, true);

  gsInfo<<"Bisection iterations --- Starting with U.norm = "<<m_U.norm()<<" and L = "<<m_L<<"; Reference error = "<<referenceError<<"\n";

  index_t fa, fb;
  fa = _bisectionObjectiveFunction(m_U, false); // jacobian is already computed in the termination function

  // m_arcLength = m_arcLength/2.0;

  // Store termination function value of initial solution U
  // T termU =  _bisectionTerminationFunction(m_U);
  for (index_t k = 1; k < m_maxIterations; ++k)
  {
    if (m_verbose) gsInfo<<"\t bisection iteration "<<k<<"\t arc length = "<<m_arcLength<<"; ";

    // Reset start point
    this->setSolution(U_old, L_old);

    gsInfo<<"From U.norm = "<<m_U.norm()<<" and L = "<<m_L<<"\n";

    // Make an arc length step; m_U and U_old and m_L and L_old are different
    gsStatus status = step();

    // Objective function on new point
    fb = _bisectionObjectiveFunction(m_U, true); // m_u is the new solution

    if (fa==fb && status==gsStatus::Success)
    {
      U_old = m_U;
      L_old = m_L;
      // Objective function on previous point
      fa = fb;
    }
    else
      this->reduceLength(0.5); // ONLY IF fa==fb??

    // termination criteria
    // relative error
    T term = _bisectionTerminationFunction(m_U, false); // jacobian on the new point already computed

    if (m_verbose) gsInfo<<"\t Finished. Relative error = "<< abs(term/referenceError)<<"\t"<<" obj.value = "<<term<<"\n";

    // if (  m_arcLength < tol && abs(term) < tol )
    if (  abs(term/referenceError) < tol )
    {
        converged = true;
        break;
    }

  }
  restore.release();

  // Compute eigenvector
  m_V = gsVector<T>::Ones(m_numDof);
  m_V.normalize();
  gsVector<T> Vold = gsVector<T>::Zero(m_numDof);
  for (index_t k = 0; k<5; k++)
  {
    this->factorizeMatrix(m_jacMat);

    m_V = this->solveSystem(m_V);
    m_V.normalize();
    if (  (m_V-Vold).norm() < m_tolerance )
    {
        converged = true;
        break;
    }
    Vold = m_V;
  }
  return converged;
}

template <class T>
void gsALMBase<T>::switchBranch()
{
  if (m_constraint && !this->_vectorLike())
  {
    gsWarn<<"switchBranch needs a copyable state space that exports a vector of the tangent size.\n";
    return;
  }

  m_V.normalize();
  real_t lenPhi = m_V.norm();
  real_t xi = lenPhi/m_tau;
  // gsInfo<<xi<<"\n";
  m_DeltaU = xi*m_V;
  m_U = m_U + m_DeltaU;

  m_DeltaLold = 0.0;
  m_DeltaUold.setZero();
  // m_DeltaLold = m_DeltaL;
  // m_DeltaUold = m_DeltaU;

  if (m_constraint)
  {
    m_prevKnown = false;
    if (!this->_importMirror())
      gsWarn<<"switchBranch: the new point could not be imported into the state space.\n";
    m_driver->resetHistory();
  }
  m_stability = m_stabilityPrev = 0;
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
