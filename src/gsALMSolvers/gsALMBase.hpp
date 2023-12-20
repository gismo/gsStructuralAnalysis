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

    m_arcLength = m_arcLength_prev = m_arcLength_ori = m_options.getReal("Length");

    m_SPfail = m_options.getInt ("SingularPointFailure");
    m_SPTestTol = m_options.getReal("SingularPointTestTol");
    m_SPTestIt  = m_options.getInt ("SingularPointTestIt");

    m_SPCompTolE = m_options.getReal("SingularPointComputeTolE");
    m_SPCompTolB = m_options.getReal("SingularPointComputeTolB");


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
  T fac = m_desiredIterations / m_numIterations;
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
  m_arcLength = m_arcLength_prev = m_arcLength_ori;
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
  m_deltaUt = this->solveSystem(m_forcing);
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
  for (m_numIterations = 1; m_numIterations < m_maxIterations; ++m_numIterations)
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
      // Change arc length
      if (m_adaptiveLength)
        computeLength();
      else
        m_arcLength_prev = m_arcLength;

      break;
    }
    else if (m_numIterations == m_maxIterations-1)
    {
      gsInfo<<"maximum iterations reached. Solution did not converge\n";
      throw 1;
    }

      // GISMO_ERROR("maximum iterations reached. Solution did not converge");
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
      this->switchBranch();

    // here we assume that the stability of the singular point is
    // equal to the one of the previous point...
    // to avoid the algorithm to find a singular point again
    m_stability = m_stabilityPrev;
  }
}

// tolB and switchBranch will be defaulted
template <class T>
gsStatus gsALMBase<T>::computeSingularPoint(const gsVector<T> & U, const T & L, bool switchBranch, bool jacobian, bool testPoint)
{
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
  return m_status;
}

template <class T>
bool gsALMBase<T>::_testSingularPoint(bool jacobian)
{
  // First, approximate the eigenvector of the Jacobian by a few arc length iterations
  // Initiate m_V and m_DeltaVDET

  if (jacobian)
  {
    m_jacMat = this->computeJacobian(m_U,m_deltaU);
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
  m_stability = (m_indicator < 0) ? 1 : -1;
}

template <class T>
index_t gsALMBase<T>::stability() const
{
  return (m_indicator < 0) ? 1 : -1;
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

  m_jacMat = this->computeJacobian(m_U,m_deltaU); // Jacobian evaluated on m_U
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
  // Store original arc length
  real_t dL = m_arcLength;
  bool adaptiveBool = m_adaptiveLength;
  m_adaptiveLength = false;
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
    m_U = U_old;
    m_L = L_old;

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
      m_arcLength = m_arcLength/2.0; // ONLY IF fa==fb??

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
  // Reset arc length
  m_arcLength = dL;
  m_adaptiveLength = adaptiveBool;

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
