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
#include <gsStructuralAnalysis/gsALMHelper.h>

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

    m_options.addInt ("BifurcationMethod","Bifurcation Identification based on: 0: Determinant;  1: Eigenvalue",bifmethod::Eigenvalue);

    m_options.addInt ("SingularPointFailure","What to do wne a singular point determination fails?: 0 = Apply solution anyways; 1 = Proceed without singular point",SPfail::With);

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

    m_arcLength = m_arcLength_prev = m_options.getReal("Length");

    m_SPfail = m_options.getInt ("SingularPointFailure");
}

template <class T>
void gsALMBase<T>::init(bool stability)
{
  if (stability)
  {
    this->computeStability(m_U);
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
void gsALMBase<T>::computeResidual()
{
  m_resVec = m_residualFun(m_U + m_DeltaU, m_L + m_DeltaL, m_forcing);
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
  GISMO_ASSERT(m_solver->info()==Eigen::ComputationInfo::Success,"Solver error with code "<<m_solver->info()<<". See Eigen documentation on ComputationInfo \n"
                                                              <<Eigen::ComputationInfo::Success<<": Success"<<"\n"
                                                              <<Eigen::ComputationInfo::NumericalIssue<<": NumericalIssue"<<"\n"
                                                              <<Eigen::ComputationInfo::NoConvergence<<": NoConvergence"<<"\n"
                                                              <<Eigen::ComputationInfo::InvalidInput<<": InvalidInput"<<"\n");
}

template <class T>
gsVector<T> gsALMBase<T>::solveSystem(const gsVector<T> & F)
{
  return m_solver->solve(F);
}

template <class T>
gsSparseMatrix<T> gsALMBase<T>::_computeJacobian(const gsVector<T> U, const gsVector<T> deltaU)
{
  // Compute Jacobian
  gsSparseMatrix<T> m;
  m = m_djacobian(U,deltaU);
  this->factorizeMatrix(m);
  m_note += "J";
  return m;
}

template <class T>
void gsALMBase<T>::computeJacobian(const gsVector<T> U, const gsVector<T> deltaU)
{
  m_jacMat = _computeJacobian(U,deltaU);
}

template <class T>
void gsALMBase<T>::computeJacobian()
{
  // Compute Jacobian
  if (m_deltaU.rows() == 0)
    m_deltaU = gsVector<T>::Zero(m_DeltaU.rows());
  this->computeJacobian(m_U + m_DeltaU, m_deltaU);
}

template <class T> void gsALMBase<T>::computeUbar() { m_deltaUbar = this->solveSystem(-m_resVec); }
template <class T> void gsALMBase<T>::computeUt()   { m_deltaUt = this->solveSystem(m_forcing); }

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

template <class T>
void gsALMBase<T>::step()
{
  GISMO_ASSERT(m_initialized,"Arc-Length Method is not initialized! Call initialize()");

  if (m_verbose)
    initOutput();

  m_converged = false;
  m_numIterations = 0;
  initiateStep();
  predictor();

  computeResidual();
  computeResidualNorms();

  if (m_verbose)
     stepOutput();

  if (m_quasiNewton)
  {
    quasiNewtonPredictor();
  }

  for (m_numIterations = 1; m_numIterations < m_maxIterations; ++m_numIterations)
  {
    if ( (!m_quasiNewton) || ( ( m_quasiNewtonInterval>0 ) && ( m_numIterations % m_quasiNewtonInterval) < 1e-10 ) )
    {
      quasiNewtonIteration();
    }

    iteration();

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
      gsInfo<<"maximum iterations reached. Solution did not converge\n";

      // GISMO_ERROR("maximum iterations reached. Solution did not converge");
  }
}

// // ------------------------------------------------------------------------------------------------------------
// // ---------------------------------------Load Control method--------------------------------------------------------
// // ------------------------------------------------------------------------------------------------------------

// template <class T>
// void gsALMBase<T>::iterationLC()
// {
//   m_deltaU = m_deltaUbar;
//   m_DeltaU += m_deltaU;
// }

// template <class T>
// void gsALMBase<T>::initiateStepLC()
// {
//   // (m_U,m_L) is the present solution (iteratively updated)
//   // (m_Uprev,m_Lprev) is the previously converged solution before (m_Lprev,m_Uprev)

//   // Reset step
//   m_DeltaU = m_deltaU = gsVector<T>::Zero(m_numDof);
//   m_DeltaL = m_deltaL = 0.0;

//   // Initiate verbose
//   if (m_verbose)
//     initOutputLC();
// }

// template <class T>
// void gsALMBase<T>::predictorLC()
// {
//   m_DeltaL = m_deltaL = m_arcLength;
//   m_deltaUt = this->solveSystem(m_forcing);
//   m_DeltaU = m_deltaL*m_deltaUt;

//   // gsDebugVar(m_DeltaU.norm());
//   m_note += "predictor\t";
// }

// template <class T>
// void gsALMBase<T>::iterationFinishLC()
// {
//   m_converged = true;
//   m_Uprev = m_U;
//   m_Lprev = m_L;
//   m_U += m_DeltaU;
//   m_L += m_DeltaL;
// }

// // ------------------------------------------------------------------------------------------------------------
// // ---------------------------------------Riks's method--------------------------------------------------------
// // ------------------------------------------------------------------------------------------------------------

// // template <class T>
// // void gsALMBase<T>::iterationRiks()
// // {
// //   m_deltaUbar = this->solveSystem(-m_resVec);
// //   m_deltaUt = this->solveSystem(m_forcing); // note the minus!!

// //   m_deltaL = - ( (m_DeltaU).dot(m_deltaUbar) ) / ( (m_DeltaL)*m_forcing.dot(m_forcing) + (m_DeltaU).dot(m_deltaUt) );
// //   m_deltaU = m_deltaUbar + m_deltaL*m_deltaUt;

// //   m_DeltaU += m_deltaU;
// //   m_DeltaL += m_deltaL;
// // }

// template <class T>
// void gsALMBase<T>::iterationRiks()
// {
//   // Residual function
//   // r = m_phi * m_DeltaU.dot(m_DeltaU)  + (1.0-m_phi) * m_DeltaL*m_DeltaL - m_arcLength*m_arcLength;
//   T r = m_phi*(m_U + m_DeltaU - m_U).dot(m_U + m_DeltaU - m_U) + (1-m_phi)*math::pow(m_L + m_DeltaL - m_L,2.0) - m_arcLength*m_arcLength;

//   m_deltaL = - ( r + 2*m_phi*(m_DeltaU).dot(m_deltaUbar) ) / ( 2*(1-m_phi)*(m_DeltaL) + 2*m_phi*(m_DeltaU).dot(m_deltaUt) );
//   m_deltaU = m_deltaUbar + m_deltaL*m_deltaUt;

//   m_DeltaU += m_deltaU;
//   m_DeltaL += m_deltaL;
// }

// template <class T>
// void gsALMBase<T>::initiateStepRiks()
// {
//   // (m_U,m_L) is the present solution (iteratively updated)
//   // (m_Uprev,m_Lprev) is the previously converged solution before (m_Lprev,m_Uprev)

//   // Reset step
//   m_DeltaU = m_deltaU =  gsVector<T>::Zero(m_numDof);
//   m_DeltaL = m_deltaL = 0.0;

//   // Initiate verbose
//   if (m_verbose)
//     initOutputRiks();
// }

// template <class T>
// void gsALMBase<T>::predictorRiks()
// {
//   // Define scaling
//   if (m_numDof ==1)
//     m_phi = 0.99999999999;
//   else
//     m_phi = 1./m_numDof;

//   // Check if the solution on start and prev are similar.
//   // Then compute predictor of the method
//   T tol = 1e-10;

//   if ( ((m_U-m_Uprev).norm() < tol) && ((m_L - m_Lprev) * (m_L - m_Lprev) < tol ) )
//   {
//     m_note += "predictor\t";
//     T DL = 1.;
//     m_deltaUt = this->solveSystem(m_forcing);
//     m_deltaU = m_deltaUt / math::sqrt( m_deltaUt.dot(m_deltaUt) + m_DeltaL*DL );
//     m_deltaL = DL / math::sqrt( m_deltaUt.dot(m_deltaUt) + m_DeltaL*DL );
//   }
//   else
//   {
//     m_deltaL = 1./m_arcLength_prev*(m_L - m_Lprev);
//     m_deltaU = 1./m_arcLength_prev*(m_U - m_Uprev);
//   }

//   // Update iterative step
//   m_deltaL *= m_arcLength;
//   m_deltaU *= m_arcLength;

//   // Update load step
//   m_DeltaU += m_deltaU;
//   m_DeltaL += m_deltaL;
// }

// // template <class T>
// // void gsALMBase<T>::predictorRiks()
// // {
// //   // Check if the solution on start and prev are similar.
// //   // Then compute predictor of the method
// //   T tol = 1e-10;

// //   if ( ((m_U-m_Uprev).norm() < tol) && ((m_L - m_Lprev) * (m_L - m_Lprev) < tol ) )
// //   {
// //     m_note += "predictor\t";
// //     T DL = 1.;
// //     computeJacobian();
// //     m_deltaUt = this->solveSystem(m_forcing);
// //     m_deltaU = m_deltaUt / math::sqrt( m_deltaUt.dot(m_deltaUt) + m_DeltaL*DL );
// //     m_deltaL = DL / math::sqrt( m_deltaUt.dot(m_deltaUt) + m_DeltaL*DL );
// //   }
// //   else
// //   {
// //     m_deltaL = 1./m_arcLength*(m_L - m_Lprev);
// //     m_deltaU = 1./m_arcLength*(m_U - m_Uprev);
// //   }

// //   // Update iterative step
// //   m_deltaL *= m_arcLength;
// //   m_deltaU *= m_arcLength;

// //   // Update load step
// //   m_DeltaU += m_deltaU;
// //   m_DeltaL += m_deltaL;
// // }

// template <class T>
// void gsALMBase<T>::iterationFinishRiks()
// {
//   m_converged = true;
//   m_Uprev = m_U;
//   m_Lprev = m_L;
//   m_U += m_DeltaU;
//   m_L += m_DeltaL;
// }

// // ------------------------------------------------------------------------------------------------------------
// // ---------------------------------------Consistent Crisfield's method----------------------------------------
// // ------------------------------------------------------------------------------------------------------------

// template <class T>
// void gsALMBase<T>::iterationConsistentCrisfield() // see Carrera1994 eqs 20 and 21
// {
//   T cres = m_DeltaU.dot(m_DeltaU) + m_phi*m_phi*m_DeltaL*m_DeltaL * m_forcing.dot(m_forcing) - m_arcLength*m_arcLength;
//   T num = cres + (2*m_DeltaU).dot(m_deltaUbar);
//   T denum = (2*m_DeltaU).dot(m_deltaUt) + m_phi*m_phi*2*m_DeltaL*m_forcing.dot(m_forcing);
//   m_deltaL = - num / denum;
//   m_deltaU = m_deltaL * m_deltaUt + m_deltaUbar;

//   m_DeltaU += m_deltaU;
//   m_DeltaL += m_deltaL;
// }

// template <class T>
// void gsALMBase<T>::initiateStepConsistentCrisfield()
// {
//   // (m_U,m_L) is the present solution (iteratively updated)
//   // (m_Uprev,m_Lprev) is the previously converged solution before (m_Lprev,m_Uprev)

//   // Reset step
//   m_DeltaU = m_deltaU =  gsVector<T>::Zero(m_numDof);
//   m_DeltaL = m_deltaL = 0.0;

//   // Initiate verbose
//   if (m_verbose)
//     initOutputRiks();
// }

// template <class T>
// void gsALMBase<T>::predictorConsistentCrisfield()
// {
//   // Check if the solution on start and prev are similar.
//   // Then compute predictor of the method
//   T tol = 1e-10;

//   if ( ((m_U-m_Uprev).norm() < tol) && ((m_L - m_Lprev) * (m_L - m_Lprev) < tol ) )
//   {
//     m_note += "predictor\t";
//     T DL = 1.;
//     m_deltaUt = this->solveSystem(m_forcing);
//     m_deltaU = m_deltaUt / math::sqrt( m_deltaUt.dot(m_deltaUt) + m_DeltaL*DL );
//     m_deltaL = DL / math::sqrt( m_deltaUt.dot(m_deltaUt) + m_DeltaL*DL );

//     if (!m_phi_user)
//       m_phi = math::pow( m_deltaUt.dot(m_deltaUt) / m_forcing.dot(m_forcing),0.5);
//   }
//   else
//   {
//     m_deltaL = 1./m_arcLength_prev*(m_L - m_Lprev);
//     m_deltaU = 1./m_arcLength_prev*(m_U - m_Uprev);

//     if (!m_phi_user)
//       m_phi = math::pow( m_deltaUt.dot(m_deltaUt) / m_forcing.dot(m_forcing),0.5);
//   }

//   // Update iterative step
//   m_deltaL *= m_arcLength;
//   m_deltaU *= m_arcLength;

//   // Update load step
//   m_DeltaU += m_deltaU;
//   m_DeltaL += m_deltaL;
// }

// template <class T>
// void gsALMBase<T>::iterationFinishConsistentCrisfield()
// {
//   m_converged = true;
//   m_Uprev = m_U;
//   m_Lprev = m_L;
//   m_U += m_DeltaU;
//   m_L += m_DeltaL;
//   m_DeltaUold = m_DeltaU;
//   m_DeltaLold = m_DeltaL;
// }

// // ------------------------------------------------------------------------------------------------------------
// // ---------------------------------------ExplicitIterations Method--------------------------------------------
// // ------------------------------------------------------------------------------------------------------------

// template <class T>
// void gsALMBase<T>::iterationExplicitIterations() // see Forde 1987 page 627
// {
//   // first do the Riks-Wempner-Tangent Method (Carrera1994 section 4.1)
//   T num = (m_DeltaU).dot(m_deltaUbar);
//   T denum = (m_DeltaU).dot(m_deltaUt) + m_phi*m_phi*m_DeltaL*m_forcing.dot(m_forcing);
//   m_deltaL = - num / denum;
//   m_deltaU = m_deltaL * m_deltaUt + m_deltaUbar;

//   // find the length of the tangent in the potential configuration
//   m_tangentLength = math::sqrt( m_tangentLength * m_tangentLength + m_deltaU.dot(m_deltaU) + m_phi*m_phi*m_deltaL*m_deltaL );

//   // residual of tangents
//   T res = ( m_arcLength*m_arcLength / ( m_tangentLength*m_tangentLength ) ) * ( m_tangentLength - m_arcLength );

//   m_deltaL = - ( res + m_DeltaU.dot(m_deltaUbar) ) / (m_DeltaU).dot(m_deltaUt) + m_phi*m_phi*m_DeltaL*m_forcing.dot(m_forcing);
//   m_deltaU = m_deltaL * m_deltaUt + m_deltaUbar;

//   m_DeltaU += m_deltaU;
//   m_DeltaL += m_deltaL;
// }

// template <class T>
// void gsALMBase<T>::initiateStepExplicitIterations()
// {
//   // (m_U,m_L) is the present solution (iteratively updated)
//   // (m_Uprev,m_Lprev) is the previously converged solution before (m_Lprev,m_Uprev)

//   // Reset step
//   m_DeltaU = m_deltaU =  gsVector<T>::Zero(m_numDof);
//   m_DeltaL = m_deltaL = 0.0;

//   // Initiate verbose
//   if (m_verbose)
//     initOutputRiks();
// }

// template <class T>
// void gsALMBase<T>::predictorExplicitIterations()
// {
//   // Check if the solution on start and prev are similar.
//   // Then compute predictor of the method
//   T tol = 1e-10;

//   if ( ((m_U-m_Uprev).norm() < tol) && ((m_L - m_Lprev) * (m_L - m_Lprev) < tol ) )
//   {
//     m_note += "predictor\t";
//     T DL = 1.;
//     m_deltaUt = this->solveSystem(m_forcing);
//     m_deltaU = m_deltaUt / math::sqrt( m_deltaUt.dot(m_deltaUt) + m_DeltaL*DL );
//     m_deltaL = DL / math::sqrt( m_deltaUt.dot(m_deltaUt) + m_DeltaL*DL );

//     if (!m_phi_user)
//       m_phi = math::pow( m_deltaUt.dot(m_deltaUt) / m_forcing.dot(m_forcing),0.5);
//   }
//   else
//   {
//     m_deltaL = 1./m_arcLength_prev*(m_L - m_Lprev);
//     m_deltaU = 1./m_arcLength_prev*(m_U - m_Uprev);

//     if (!m_phi_user)
//       m_phi = math::pow( m_deltaUt.dot(m_deltaUt) / m_forcing.dot(m_forcing),0.5);
//   }

//   // Update iterative step
//   m_deltaL *= m_arcLength;
//   m_deltaU *= m_arcLength;

//   // Update load step
//   m_DeltaU += m_deltaU;
//   m_DeltaL += m_deltaL;

//   m_tangentLength = math::sqrt( m_DeltaU.dot(m_DeltaU) + m_phi*m_phi*m_DeltaL*m_DeltaL );
// }

// template <class T>
// void gsALMBase<T>::iterationFinishExplicitIterations()
// {
//   m_converged = true;
//   m_Uprev = m_U;
//   m_Lprev = m_L;
//   m_U += m_DeltaU;
//   m_L += m_DeltaL;
//   m_DeltaUold = m_DeltaU;
//   m_DeltaLold = m_DeltaL;
// }

// // ------------------------------------------------------------------------------------------------------------
// // ---------------------------------------Crisfield's method---------------------------------------------------
// // ------------------------------------------------------------------------------------------------------------
// template <class T>
// void gsALMBase<T>::iterationCrisfield()
// {
//   m_eta = 1.0;

//   T lamold = m_deltaL;
//   computeLambdas();

//   // Relaxation against oscillating load factor
//   if (( (lamold*m_deltaL < 0) && (abs(m_deltaL) <= abs(lamold) ) ) && m_relax != 1.0 )
//   {
//      m_note += "\t relaxated solution!";
//      m_deltaU = m_relax * (m_deltaL*m_deltaUt + m_eta*m_deltaUbar);
//      m_deltaL = m_relax * m_deltaL;
//   }

//   m_DeltaU += m_deltaU;
//   m_DeltaL += m_deltaL;

//   if (m_angleDetermine == angmethod::Iteration)
//   {
//     m_DeltaUold = m_DeltaU;
//     m_DeltaLold = m_DeltaL;
//   }
// }

// template <class T>
// void gsALMBase<T>::initiateStepCrisfield()
// {
//   if (m_verbose)
//     initOutputCrisfield();

//   m_DeltaU = m_deltaUbar = m_deltaUt = gsVector<T>::Zero(m_numDof);
//   m_DeltaL = m_deltaL = 0.0;
//   m_eta = 1.0;
// }

// template <class T>
// void gsALMBase<T>::predictorCrisfield()
// {
//   m_deltaUt = this->solveSystem(m_forcing);

//   // Choose Solution
//   if (m_DeltaUold.dot(m_DeltaUold) == 0 && m_DeltaLold*m_DeltaLold == 0) // no information about previous step.
//   // {
//   //  m_deltaL = m_arcLength / math::pow(2*( m_deltaUt.dot(m_deltaUt) ) , 0.5);
//   //  m_deltaU = m_deltaUbar + m_deltaL*m_deltaUt;
//   //  m_phi = math::pow( m_deltaUt.dot(m_deltaUt) / m_forcing.dot(m_forcing),0.5);

//   // }
//   {
//     m_note += "predictor\t";
//     T DL = 1.;
//     // m_deltaL = m_arcLength * DL / math::sqrt( 2*( m_deltaUt.dot(m_deltaUt) + m_DeltaL*DL ) );
//     m_deltaL = m_arcLength * DL / math::sqrt( ( m_deltaUt.dot(m_deltaUt) + m_DeltaL*DL ) );

//     // m_deltaU = m_arcLength * m_deltaUt / math::sqrt( m_deltaUt.dot(m_deltaUt) + m_DeltaL*DL );
//     m_deltaU = m_deltaL*m_deltaUt;

//     if (!m_phi_user)
//       m_phi = math::pow( m_deltaUt.dot(m_deltaUt) / m_forcing.dot(m_forcing),0.5);

//     m_DeltaUold = m_deltaU;
//     m_DeltaLold = m_deltaL;
//   }
//   else // previous point is not in the origin
//   {
//     if (!m_phi_user)
//       m_phi = math::pow(m_U.dot(m_U)/( math::pow(m_L,2) * m_forcing.dot(m_forcing) ),0.5);
//     m_note += " phi=" + std::to_string(m_phi);
//     computeLambdaMU();
//   }

//   // Compute Temporary updates of DeltaL and DeltaU
//   m_DeltaU += m_deltaU;
//   m_DeltaL += m_deltaL;

//   if (m_angleDetermine == angmethod::Iteration || m_angleDetermine == angmethod::Predictor)
//   {
//    m_DeltaUold = m_DeltaU;
//    m_DeltaLold = m_DeltaL;
//   }
// }

// template <class T>
// void gsALMBase<T>::iterationFinishCrisfield()
// {
//   m_converged = true;
//   m_Uprev = m_U;
//   m_Lprev = m_L;
//   m_U += m_DeltaU;
//   m_L += m_DeltaL;
//   if (m_angleDetermine == angmethod::Step)
//   {
//     m_DeltaUold = m_DeltaU;
//     m_DeltaLold = m_DeltaL;
//   }
// }

// template <class T>
// void gsALMBase<T>::computeLambdasSimple() //Ritto-Corrêa et al. 2008
// {
//   T A0 = math::pow(m_phi,2)* m_forcing.dot(m_forcing); // see Lam et al. 1991,

//   m_a0 = m_deltaUt.dot(m_deltaUt) + A0;
//   m_b0 = 2*( m_deltaUt.dot(m_DeltaU) + m_DeltaL * A0 );
//   m_b1 = 2*( m_deltaUbar.dot(m_deltaUt) );
//   m_c0 = m_DeltaU.dot(m_DeltaU) + m_DeltaL*m_DeltaL * A0 - math::pow(m_arcLength,2);
//   m_c1 = 2*( m_DeltaU.dot(m_deltaUbar) );
//   m_c2 = m_deltaUbar.dot(m_deltaUbar);

//   /// Calculate the coefficients of the polynomial
//   m_alpha1 = m_a0;
//   m_alpha2 = m_b0 + m_eta*m_b1;
//   m_alpha3 = m_c0 + m_eta*m_c1 + m_eta*m_eta*m_c2;

//   m_discriminant = math::pow(m_alpha2 ,2) - 4 * m_alpha1 * m_alpha3;
//   // m_note += "\t D = " + std::to_string(m_discriminant);
// }

// template <class T>
// void gsALMBase<T>::computeLambdasEta()
// {
//   m_alpha1 = m_a0;
//   m_alpha2 = m_b0 + m_eta*m_b1;
//   m_alpha3 = m_c0 + m_eta*m_c1 + m_eta*m_eta*m_c2;

//   // Lam 1992
//   m_discriminant = math::pow(m_alpha2 ,2) - 4 * m_alpha1 * m_alpha3;
//   m_deltaLs[0] = (-m_alpha2 + math::sqrt(m_discriminant))/(2*m_alpha1);
//   m_deltaLs[1] = (-m_alpha2 - math::sqrt(m_discriminant))/(2*m_alpha1);

//   // Zhou 1995
//   // m_deltaLs[0] = (-m_alpha2 )/(2*m_alpha1);
//   // m_deltaLs[1] = (-m_alpha2 )/(2*m_alpha1);
// }

// template <class T>
// void gsALMBase<T>::computeLambdasModified()
// {
//   m_alpha1 = m_b1*m_b1 - 4.0*m_a0*m_c2;
//   m_alpha2 = 2.0*m_b0*m_b1 - 4.0*m_a0*m_c1;
//   m_alpha3 = m_b0*m_b0 - 4.0*m_a0*m_c0;

//   m_discriminant = math::pow(m_alpha2 ,2.0) - 4.0 * m_alpha1 * m_alpha3;

//   gsVector<T> etas(2);
//   etas.setZero();
//   if (m_discriminant >= 0)
//     {
//       etas[0] = (-m_alpha2 + math::pow(m_discriminant,0.5))/(2.0*m_alpha1);
//       etas[1] = (-m_alpha2 - math::pow(m_discriminant,0.5))/(2.0*m_alpha1);

//       T eta1 = std::min(etas[0],etas[1]);
//       T eta2 = std::max(etas[0],etas[1]);
//       gsInfo<<"eta 1 = "<<eta1<<"\t eta2 = "<<eta2<<"\n";

//       // Approach of Zhou 1995
//       // m_eta = std::min(1.0,eta2);
//       // if (m_eta <= 0)
//       //   gsInfo<<"Warning: both etas are non-positive!\n";
//       // if (m_eta <=0.5)
//       // {
//       //   gsInfo<<"Warning: eta is small; bisecting step length\n";
//       //   m_arcLength=m_arcLength/2.;
//       // }

//       // Approach of Lam 1992
//       T xi = 0.05*abs(eta2-eta1);
//       if (eta2<1.0)
//         m_eta = eta2-xi;
//       else if ( (eta2 > 1.0) && (-m_alpha2/m_alpha1 < 1.0) )
//         m_eta = eta2+xi;
//       else if ( (eta1 < 1.0) && (-m_alpha2/m_alpha1 > 1.0) )
//         m_eta = eta1-xi;
//       else if ( eta1 > 1.0 )
//         m_eta = eta1 + xi;

//       if (eta2<1.0)
//         m_note = note + " option 1";
//       else if ( (eta2 > 1.0) && (-m_alpha2/m_alpha1 < 1.0) )
//         m_note = note + " option 2";
//       else if ( (eta1 < 1.0) && (-m_alpha2/m_alpha1 > 1.0) )
//         m_note = note + " option 3";
//       else if ( eta1 > 1.0 )
//         m_note = note + " option 4";

//       // m_eta = eta2;
//     }
//   else
//   {
//     gsInfo<<"Discriminant was negative in modified method\n";
//   }
// }
// template <class T>
// void gsALMBase<T>::computeLambdasComplex()
// {
//   T A0 = math::pow(m_phi,2)* m_forcing.dot(m_forcing); // see Lam et al. 1991
//   gsVector<T> DeltaUcr = m_DeltaU + m_deltaUbar;

//   // Compute internal loads from residual and
//   // gsVector<T> R = m_residualFun(m_U + DeltaUcr, m_L + m_DeltaL , m_forcing);
//   // gsVector<T> Fint = R + (m_L + m_DeltaL) * m_forcing;
//   gsVector<T> Fint = m_jacMat*(m_U+m_DeltaU);
//   T Lcr = Fint.dot(m_forcing)/m_forcing.dot(m_forcing);
//   T DeltaLcr = Lcr - m_L;

//   T arcLength_cr = math::pow( DeltaUcr.dot(DeltaUcr) + A0 * math::pow(DeltaLcr, 2.0) ,0.5);
//   T mu = m_arcLength/arcLength_cr;

//   m_deltaL = mu*DeltaLcr - m_DeltaL;
//   m_deltaU = mu*DeltaUcr - m_DeltaU;
// }

// template <class T>
// void gsALMBase<T>::computeLambdas()
// {
//   m_deltaLs.setZero(2);
//   computeLambdasSimple();
//   if (m_discriminant >= 0)
//   {
//     m_eta = 1.0;
//     m_deltaLs[0] = (-m_alpha2 + math::pow(m_discriminant,0.5))/(2*m_alpha1);
//     m_deltaLs[1] = (-m_alpha2 - math::pow(m_discriminant,0.5))/(2*m_alpha1);
//     computeLambdaDOT();
//   }
//   else
//   {
//     m_note += "\tC";
//     // Compute eta
//     computeLambdasModified();
//     if ((m_discriminant >= 0) && m_eta > 0.5)
//     {
//       // recompute lambdas with new eta
//       computeLambdasEta();
//       // gsInfo<<"2: dL1 = "<<m_deltaLs[0]<<"\tdL2 = "<<m_deltaLs[1]<<"\t eta = "<<m_eta<<"\n";
//       computeLambdaDOT();
//       // gsInfo<<"2: dL1 = "<<m_deltaL<<"\t m_deltaU.norm = "<<m_deltaU.norm()<<"\t eta = "<<m_eta<<"\n";
//       gsInfo<<"Modified Complex Root Solve\n";
//     }
//     else
//     {
//       // if the roots of the modified method are still complex, we use the following function (see Lam 1992, eq 13-17)
//       m_eta = 1.0;
//       computeLambdasComplex();
//       // gsInfo<<"3: dL1 = "<<m_deltaL<<"\t m_deltaU.norm = "<<m_deltaU.norm()<<"\t eta = "<<m_eta<<"\n";
//       gsInfo<<"Simplified Complex Root Solve\n";
//       // Note: no selection of roots is needed
//     }
//   }
// }

// template <class T>
// void gsALMBase<T>::computeLambdaDET()
// {
//     computeLambdas();

//     if (sign(m_DeltaL + m_deltaLs[0]) == sign(m_detKT))
//       m_deltaL = m_deltaLs[0];
//     else
//       m_deltaL = m_deltaLs[1];

//     // Compute update of U (NOTE: m_eta=1.0)
//     m_deltaU = m_deltaUbar + m_deltaL*m_deltaUt;

//     // gsInfo<<"\t\t Choice based on DETERMINANT. Options:\n";
//     // gsInfo<<"\t\t DeltaL = "<<m_DeltaL+m_deltaLs[0]<<" DeltaU.norm = "<<(m_DeltaU + m_deltaUbar + m_deltaLs[0]*m_deltaUt).norm()<<"\n";
//     // gsInfo<<"\t\t DeltaL = "<<m_DeltaL+m_deltaLs[1]<<" DeltaU.norm = "<<(m_DeltaU + m_deltaUbar + m_deltaLs[1]*m_deltaUt).norm()<<"\n";
// }

// template <class T>
// void gsALMBase<T>::computeLambdaMU()
// {
//     T A0 = math::pow(m_phi,2)* m_forcing.dot(m_forcing); // see Lam et al. 1991
//     index_t dir = sign(m_DeltaUold.dot(m_deltaUt) + A0*m_DeltaLold); // Feng et al. 1995 with H = \Psi^2
//     T denum = ( math::pow( m_deltaUt.dot(m_deltaUt) + A0 ,0.5) ); // Feng et al. 1995 with H = \Psi^2

//     T mu;
//     if (denum==0)
//       mu = m_arcLength;
//     else
//       mu = m_arcLength / denum;

//     m_deltaL = dir*mu;
//     m_deltaU = m_deltaL*m_deltaUt;
// }

// template <class T>
// void gsALMBase<T>::computeLambdaDOT()
// {
//     gsVector<T> deltaU1, deltaU2;
//     deltaU1 = m_eta*m_deltaUbar + m_deltaUt*m_deltaLs[0];
//     deltaU2 = m_eta*m_deltaUbar + m_deltaUt*m_deltaLs[1];

//     // ---------------------------------------------------------------------------------
//     // Method by Ritto-Corea et al. 2008
//     T DOT1,DOT2;
//     DOT1 = m_deltaLs[0]*(m_DeltaUold.dot(m_deltaUt) + math::pow(m_phi,2)*m_DeltaLold);
//     DOT2 = m_deltaLs[1]*(m_DeltaUold.dot(m_deltaUt) + math::pow(m_phi,2)*m_DeltaLold);

//     if (DOT1 > DOT2)
//     {
//       m_deltaL = m_deltaLs[0];
//       m_deltaU = deltaU1;
//     }
//     else if (DOT1 < DOT2)
//     {
//       m_deltaL = m_deltaLs[1];
//       m_deltaU = deltaU2;
//     }
//     else
//     {
//       m_deltaL = m_deltaLs[0];
//       m_deltaU = deltaU1;
//     }

//     // ---------------------------------------------------------------------------------
//     // // Method by Crisfield 1981
//     // T DOT1,DOT2;
//     // DOT1 = (m_DeltaUold+deltaU1).dot(m_DeltaUold);
//     // DOT2 = (m_DeltaUold+deltaU2).dot(m_DeltaUold);

//     // DOT1 = (m_DeltaU+deltaU1).dot(m_DeltaUold);
//     // DOT2 = (m_DeltaU+deltaU2).dot(m_DeltaUold);

//     // if ((DOT1 > DOT2) && DOT2 <= 0)
//     // {
//     //   m_deltaL = m_deltaLs[0];
//     //   m_deltaU = deltaU1;
//     // }
//     // else if ((DOT1 < DOT2) && DOT1 <= 0)
//     // {
//     //   m_deltaL = m_deltaLs[1];
//     //   m_deltaU = deltaU2;
//     // }
//     // else if ((DOT1 >=0) && (DOT2 >=0))
//     // {
//     //   T linsol = -m_alpha3/m_alpha2;
//     //   T diff1 = abs(m_deltaLs[0]-linsol);
//     //   T diff2 = abs(m_deltaLs[1]-linsol);
//     //   m_note += "\t linear solution!";
//     //   // m_note += "Linsol\t" + std::to_string(diff1) + "\t" + std::to_string(diff2) + "\t" + std::to_string(linsol) + "\t" + std::to_string(m_deltaLs[0]) + "\t" + std::to_string(m_deltaLs[1]) + "\n";
//     //   if (diff1 > diff2)
//     //   {
//     //     m_deltaL = m_deltaLs[1];
//     //     m_deltaU = deltaU2;
//     //   }
//     //   else
//     //   {
//     //     m_deltaL = m_deltaLs[0];
//     //     m_deltaU = deltaU1;
//     //   }
//     // }
//     // else
//     // {
//     //   m_deltaL = m_deltaLs[0];
//     //   m_deltaU = deltaU1;
//     // }

// }

// ------------------------------------------------------------------------------------------------------------
// ---------------------------------------Singular point methods-----------------------------------------------
// ------------------------------------------------------------------------------------------------------------
template <class T>
void gsALMBase<T>::computeSingularPoint(T singTol, index_t kmax, gsVector<T> U, T L, T tolE, T tolB, bool switchBranch, bool jacobian)
{
  // Controls the singular point test with the first two arguments
  bool test = this->testSingularPoint(singTol, kmax,jacobian);
  // if (m_verbose)
  // {
  //  T value = this->bisectionTerminationFunction(m_U,false);
  //  gsInfo<<"\t Singular point details:";
  //  gsInfo<<"\tvalue: "<<value<<"\n";
  // }
  if (test)
  {
    if (m_verbose) {gsInfo<<"\t Bifurcation point\n";}

    // First stage: bisection method
    if (tolB != 0)
    {
      this->bisectionSolve(U,L,tolB);
      this->extendedSystemSolve(m_U, m_L, tolE);
    }
    else
      this->extendedSystemSolve(U, L, tolE);

    if (switchBranch && (m_converged || m_SPfail==1))
      this->switchBranch();

    // here we assume that the stability of the singular point is
    // equal to the one of the previous point...
    // to avoid the algorithm to find a singular point again
    m_stability = m_stabilityPrev;
  }
  else
    gsInfo<<"\t Limit point\n";
}

template <class T>
gsMatrix<T> gsALMBase<T>::computeModes(gsVector<T> x, bool jacobian, T shift)
{
  if (jacobian)
  {
    gsVector<T> dx = gsVector<T>::Zero(x.size());
    this->computeJacobian(x,dx);
  } // otherwise the jacobian is already computed (on m_U+m_DeltaU)

  gsMatrix<T> result;
#ifdef GISMO_WITH_SPECTRA
  index_t number = std::min(static_cast<index_t>(std::floor(m_jacMat.cols()/3.)),10);
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
  GISMO_ENSURE(es.info()==Spectra::CompInfo::Successful,"Spectra did not converge!"); // Reason for not converging can be due to the value of ncv (last input in the class member), which is too low.

  // if (es.info()==Spectra::CompInfo::NotComputed)
  // if (es.info()==Spectra::CompInfo::NotConverging)
  // if (es.info()==Spectra::CompInfo::NumericalIssue)
  // Eigen::SelfAdjointEigenSolver< gsMatrix<T> > es(m_jacMat);
  m_stabilityVec = es.eigenvalues();

  result.resize(m_numDof,countNegatives(m_stabilityVec));
  if (result.cols() == number)
    gsWarn<<"Number of computed eigenvectors is too small!";

  for (index_t k=0; k!=result.cols(); ++k)
  {
    if (m_stabilityVec.at(k) < 0)
      result.col(k) = es.eigenvectors().col(k).normalized();
    else
      break;
  }

#else
  Eigen::SelfAdjointEigenSolver<gsMatrix<T>> es2(m_jacMat);
  m_stabilityVec = es2.eigenvalues();

  result.resize(m_numDof,countNegatives(m_stabilityVec));
  for (index_t k=0; k!=result.cols(); ++k)
  {
    if (m_stabilityVec.at(k) < 0)
      result.col(k) = es2.eigenvectors().col(k).normalized();
    else
      break;
  }

#endif

  return result;

}

// tolB and switchBranch will be defaulted
template <class T>
void gsALMBase<T>::computeSingularPoint(gsVector<T> U, T L, T tolE, T tolB, bool switchBranch, bool jacobian)
{ this->computeSingularPoint(1e-6, 5, U, L, tolE, tolB, switchBranch,jacobian); }

// template <class T>
// void gsALMBase<T>::computeSingularPoint(gsVector<T> U, T L, T tolE, bool switchBranch)
// { this->computeSingularPoint(1e-6, 5, U, L, tolE, 0, switchBranch); }

// template <class T>
// void gsALMBase<T>::computeSingularPoint(T singTol, index_t kmax, gsVector<T> U, T L, T tolE, bool switchBranch)
// { this->computeSingularPoint(singTol, kmax, U, L, tolE, 0, switchBranch); }

template <class T>
void gsALMBase<T>::computeSingularPoint(T singTol, index_t kmax, T tolE, T tolB, bool switchBranch, bool jacobian)
{ this->computeSingularPoint(singTol, kmax, m_U, m_L, tolE, tolB, switchBranch,jacobian); }

// template <class T>
// void gsALMBase<T>::computeSingularPoint(T singTol, index_t kmax, T tolE, bool switchBranch)
// { this->computeSingularPoint(singTol, kmax, m_U, m_L, tolE, 0, switchBranch); }

// template <class T>
// void gsALMBase<T>::computeSingularPoint(T tolE, bool switchBranch)
// { this->computeSingularPoint(1e-6, 5, m_U, m_L, tolE, 0, switchBranch); }

// template <class T>
// void gsALMBase<T>::computeSingularPoint(T tolE, T tolB, bool switchBranch)
// { this->computeSingularPoint(1e-6, 5, m_U, m_L, tolE, tolB, switchBranch); }

// to do: ADD POINT FROM WHICH TO TEST
template <class T>
bool gsALMBase<T>::testSingularPoint(gsVector<T> U, T L, T tol, index_t kmax, bool jacobian)
{
  // First, approximate the eigenvector of the Jacobian by a few arc length iterations
  // Initiate m_V and m_DeltaVDET

  if (jacobian)
    this->computeJacobian(m_U,m_deltaU);

  m_V = gsVector<T>::Ones(m_numDof);
  m_V.normalize();
  factorizeMatrix(m_jacMat);

  for (index_t k = 0; k<kmax; k++)
  {
    m_V = this->solveSystem(m_V);
    m_V.normalize();
  }

  T dot = (abs(m_V.dot(m_forcing)));
  if ( (abs(dot) / tol > 1e-1) && (abs(dot) / tol < 10) )
  {
    gsInfo<<"Warning: the singular point test is close to its tolerance. dot/tol = "<<abs(dot)/tol<<" dot = "<<dot<<"\t tolerance = "<<tol<<"\n";
  }
  if (dot < tol)    // Bifurcation point
    return true;
  else        // Limit point
    return false;
}

template <class T>
bool gsALMBase<T>::testSingularPoint(T tol, index_t kmax, bool jacobian)
{
  return testSingularPoint(m_U,m_L,tol,kmax,jacobian);
}


template <class T>
void gsALMBase<T>::computeStability(gsVector<T> x, bool jacobian, T shift)
{
  if (jacobian)
  {
    gsVector<T> dx = gsVector<T>::Zero(x.size());
    this->computeJacobian(x,dx);
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
      GISMO_ERROR("Determinant stability method only works with SimplicialLDLT solver, current solver is "<<m_options.getString("Solver"));
  }
  else if (m_bifurcationMethod == bifmethod::Eigenvalue)
  {
    #ifdef GISMO_WITH_SPECTRA
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
    GISMO_ENSURE(es.info()==Spectra::CompInfo::Successful,"Spectra did not converge!"); // Reason for not converging can be due to the value of ncv (last input in the class member), which is too low.

    // if (es.info()==Spectra::CompInfo::NotComputed)
    // if (es.info()==Spectra::CompInfo::NotConverging)
    // if (es.info()==Spectra::CompInfo::NumericalIssue)
    // Eigen::SelfAdjointEigenSolver< gsMatrix<T> > es(m_jacMat);
    m_stabilityVec = es.eigenvalues();
    #else
    Eigen::SelfAdjointEigenSolver<gsMatrix<T>> es2(m_jacMat);
    m_stabilityVec = es2.eigenvalues();
    #endif
  }
  else
    gsInfo<<"bifurcation method unknown!";

  m_negatives = countNegatives(m_stabilityVec);
  m_indicator = m_stabilityVec.colwise().minCoeff()[0]; // This is required since D does not necessarily have one column.
}

template <class T>
index_t gsALMBase<T>::stability()
{
  index_t tmp = 1;
  if (m_indicator < 0)
    tmp =  -1;
  return tmp;
}

template <class T>
index_t gsALMBase<T>::stability(gsVector<T> x, bool jacobian)
{
  this->computeStability(x, jacobian);
  index_t tmp = this->stability();
  return tmp;
}

template <class T>
bool gsALMBase<T>::stabilityChange()
{
  m_stabilityPrev = m_stability;
  m_stability = this->stability();
  if (m_stability*m_stabilityPrev < 0) // then singular point passed
    return true;
  else
    return false;
}


template <class T>
void gsALMBase<T>::extendedSystemSolve(gsVector<T> U, T L, T tol)
{
  m_U = U;
  m_L = L;
  gsInfo<<"Extended iterations --- Starting with U.norm = "<<m_U.norm()<<" and L = "<<m_L<<"\n";

  this->computeJacobian(m_U,m_deltaU); // Jacobian evaluated on m_U
  m_basisResidualKTPhi = (m_jacMat*m_V).norm();

  m_DeltaV = gsVector<T>::Zero(m_numDof);
  m_DeltaU.setZero();
  m_DeltaL = 0.0;

  m_deltaV = gsVector<T>::Zero(m_numDof);
  m_deltaU.setZero();
  m_deltaL = 0.0;
  m_converged = false;
  if (m_verbose)
      initOutputExtended();
  for (m_numIterations = 0; m_numIterations < m_maxIterations; ++m_numIterations)
  {
    extendedSystemIteration();
    m_DeltaU += m_deltaU;
    m_DeltaL += m_deltaL;
    m_DeltaV += m_deltaV;
    // m_V.normalize();

    // m_resVec = m_residualFun(m_U, m_L, m_forcing);
    // m_residue = m_resVec.norm() / ( m_L * m_forcing.norm() );
    // m_residue = (m_jacobian(m_U).toDense()*m_V).norm() / refError;
    this->computeJacobian(m_U+m_DeltaU,m_deltaU);
    m_residueKTPhi = (m_jacMat*(m_V+m_DeltaV)).norm(); // /m_basisResidualKTPhi;
    m_resVec = m_residualFun(m_U+m_DeltaU,m_L+m_DeltaL,m_forcing);
    computeResidualNorms();
    if (m_verbose)
      stepOutputExtended();

    // termination criteria
    // if ( m_residueF < m_toleranceF && m_residueU < m_toleranceU &&  m_residueKTPhi< tol )
    if ( m_residueKTPhi< tol )
    {
      m_converged = true;
      m_U += m_DeltaU;
      m_L += m_DeltaL;
      gsInfo<<"Iterations finished. U.norm() = "<<m_U.norm()<<"\t L = "<<m_L<<"\n";

      break;
    }
    if ( (m_numIterations == m_maxIterations-1) && (!m_converged) )
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
    }

  }
}

// TODO: Optimize memory
template <class T>
void gsALMBase<T>::extendedSystemIteration()
{
  m_resVec = m_residualFun(m_U+m_DeltaU, m_L+m_DeltaL, m_forcing);
  this->computeJacobian(m_U+m_DeltaU,m_deltaU);
  // m_jacMat = m_jacobian(m_U);

  m_deltaUt = this->solveSystem(m_forcing); // DeltaV1
  m_deltaUbar = this->solveSystem(-m_resVec); // DeltaV2

  real_t eps = 1e-8;
  gsSparseMatrix<T> jacMatEps = m_jacobian((m_U+m_DeltaU) + eps*(m_V+m_DeltaV));
  m_note += "J"; // mark jacobian computation
  gsVector<T> h1 = 1/eps * ( jacMatEps * m_deltaUt ) - 1/eps * m_forcing;
  gsVector<T> h2 = m_jacMat * (m_V+m_DeltaV) + 1/eps * ( jacMatEps * m_deltaUbar + m_resVec );

  factorizeMatrix(m_jacMat);

  m_deltaVt = this->solveSystem(-h1); // DeltaV1
  m_deltaVbar = this->solveSystem(-h2); // DeltaV2

  m_deltaL = -( ( (m_V+m_DeltaV)/(m_V+m_DeltaV).norm() ).dot(m_deltaVbar)  + (m_V+m_DeltaV).norm() - 1) / ( (m_V+m_DeltaV)/(m_V+m_DeltaV).norm() ).dot( m_deltaVt );

  m_deltaU = m_deltaL * m_deltaUt + m_deltaUbar;
  // gsInfo<<"m_DeltaU = \n"<<m_DeltaU<<"\n";

  m_deltaV = m_deltaL * m_deltaVt + m_deltaVbar;
  // gsInfo<<"m_DeltaV = \n"<<m_DeltaV<<"\n";
}

template <class T>
index_t gsALMBase<T>::bisectionObjectiveFunction(const gsVector<T> & x, bool jacobian)
{
  this->computeStability(x,jacobian);
  index_t negs = countNegatives(m_stabilityVec);
  m_note += "(" + std::to_string(negs) + ")";
  return negs;
}

template <class T>
T gsALMBase<T>::bisectionTerminationFunction(const gsVector<T> & x, bool jacobian)
{
  this->computeStability(x,jacobian);
  // return m_stabilityVec.colwise().minCoeff()[0]; // This is required since D does not necessarily have one column.
  return m_indicator;
}

template <class T>
void gsALMBase<T>::bisectionSolve(gsVector<T> U, T L, T tol)
{
  m_U = U;
  m_L = L;

  gsVector<T> U_old  = U;
  T L_old = L;
  // Store original arc length
  real_t dL = m_arcLength;
  bool adaptiveBool = m_adaptiveLength;
  m_adaptiveLength = false;
  m_converged = false;

  T referenceError = bisectionTerminationFunction(U, true);

  gsInfo<<"Bisection iterations --- Starting with U.norm = "<<m_U.norm()<<" and L = "<<m_L<<"; Reference error = "<<referenceError<<"\n";


  index_t fa, fb;
  fa = bisectionObjectiveFunction(m_U, false); // jacobian is already computed in the termination function

  // m_arcLength = m_arcLength/2.0;

  // Store termination function value of initial solution U
  // T termU =  bisectionTerminationFunction(m_U);
  for (index_t k = 1; k < m_maxIterations; ++k)
  {
    if (m_verbose) gsInfo<<"\t bisection iteration "<<k<<"\t arc length = "<<m_arcLength<<"; ";

    // Reset start point
    m_U = U_old;
    m_L = L_old;

    gsInfo<<"From U.norm = "<<m_U.norm()<<" and L = "<<m_L<<"\n";

    // Make an arc length step; m_U and U_old and m_L and L_old are different
    step();

    // Objective function on new point
    fb = bisectionObjectiveFunction(m_U, true); // m_u is the new solution

    if (fa==fb && m_converged)
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
    T term = bisectionTerminationFunction(m_U, false); // jacobian on the new point already computed

    if (m_verbose) gsInfo<<"\t Finished. Relative error = "<< abs(term/referenceError)<<"\t"<<" obj.value = "<<term<<"\n";

    // if (  m_arcLength < tol && abs(term) < tol )
    if (  abs(term/referenceError) < tol )
    {
        m_converged = true;
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
    factorizeMatrix(m_jacMat);

    m_V = this->solveSystem(m_V);
    m_V.normalize();
    if (  (m_V-Vold).norm() < m_tolerance )
    {
        m_converged = true;
        break;
    }
    Vold = m_V;
  }
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
void gsALMBase<T>::initOutputExtended()
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
void gsALMBase<T>::stepOutputExtended()
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
  gsInfo<<std::setw(17)<<std::left<<bisectionTerminationFunction(m_U,false);
  gsInfo<<std::setw(17)<<std::left<<m_note;
  gsInfo<<"\n";

  m_note = "";
}

} // namespace gismo