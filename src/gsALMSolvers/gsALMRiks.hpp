/** @file gsALMRiks.hpp

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
void gsALMRiks<T>::initMethods()
{
  m_numDof = m_forcing.size();
  m_DeltaU = m_U = gsVector<T>::Zero(m_numDof);
  m_DeltaL = m_L = 0.0;

  m_Uprev = gsVector<T>::Zero(m_numDof);
  m_Lprev = 0.0;
}

// ------------------------------------------------------------------------------------------------------------
// ---------------------------------------Riks's method--------------------------------------------------------
// ------------------------------------------------------------------------------------------------------------

template <class T>
void gsALMRiks<T>::quasiNewtonPredictor()
{
  m_jacMat = computeJacobian();
  this->factorizeMatrix(m_jacMat);
  computeUt(); // rhs does not depend on solution
  computeUbar(); // rhs contains residual and should be computed every time

}

template <class T>
void gsALMRiks<T>::quasiNewtonIteration()
{
  m_jacMat = computeJacobian();
  this->factorizeMatrix(m_jacMat);
  computeUt(); // rhs does not depend on solution
}

// template <class T>
// void gsALMRiks<T>::iteration()
// {
//   m_deltaUbar = this->solveSystem(-m_resVec);
//   m_deltaUt = this->solveSystem(m_forcing); // note the minus!!

//   m_deltaL = - ( (m_DeltaU).dot(m_deltaUbar) ) / ( (m_DeltaL)*m_forcing.dot(m_forcing) + (m_DeltaU).dot(m_deltaUt) );
//   m_deltaU = m_deltaUbar + m_deltaL*m_deltaUt;

//   m_DeltaU += m_deltaU;
//   m_DeltaL += m_deltaL;
// }

template <class T>
void gsALMRiks<T>::iteration()
{
  computeUbar(); // rhs contains residual and should be computed every time

  // Compute next solution
  //
  // Residual function
  // r = m_phi * m_DeltaU.dot(m_DeltaU)  + (1.0-m_phi) * m_DeltaL*m_DeltaL - m_arcLength*m_arcLength;
  T r = m_phi*(m_U + m_DeltaU - m_U).dot(m_U + m_DeltaU - m_U) + (1-m_phi)*math::pow(m_L + m_DeltaL - m_L,2.0) - m_arcLength*m_arcLength;

  m_deltaL = - ( r + 2*m_phi*(m_DeltaU).dot(m_deltaUbar) ) / ( 2*(1-m_phi)*(m_DeltaL) + 2*m_phi*(m_DeltaU).dot(m_deltaUt) );
  m_deltaU = m_deltaUbar + m_deltaL*m_deltaUt;

  m_DeltaU += m_deltaU;
  m_DeltaL += m_deltaL;
}

template <class T>
void gsALMRiks<T>::initiateStep()
{
  // (m_U,m_L) is the present solution (iteratively updated)
  // (m_Uprev,m_Lprev) is the previously converged solution before (m_Lprev,m_Uprev)

  // Reset step
  m_DeltaU = m_deltaU =  gsVector<T>::Zero(m_numDof);
  m_DeltaL = m_deltaL = 0.0;
}

template <class T>
void gsALMRiks<T>::predictor()
{
  m_jacMat = computeJacobian();
  this->factorizeMatrix(m_jacMat);

  // Define scaling
  if (m_numDof ==1)
    m_phi = 0.99999999999;
  else
    m_phi = 1./m_numDof;

  // Check if the solution on start and prev are similar.
  // Then compute predictor of the method
  T tol = 1e-10;

  if ( ((m_U-m_Uprev).norm() < tol) && ((m_L - m_Lprev) * (m_L - m_Lprev) < tol ) )
  {
    m_note+= "predictor\t";
    T DL = 1.;
    m_deltaUt = this->solveSystem(m_forcing);
    m_deltaU = m_deltaUt / math::sqrt( m_deltaUt.dot(m_deltaUt) + m_DeltaL*DL );
    m_deltaL = DL / math::sqrt( m_deltaUt.dot(m_deltaUt) + m_DeltaL*DL );
  }
  else
  {
    m_deltaL = 1./m_arcLength_prev*(m_L - m_Lprev);
    m_deltaU = 1./m_arcLength_prev*(m_U - m_Uprev);
  }

  // Update iterative step
  m_deltaL *= m_arcLength;
  m_deltaU *= m_arcLength;

  // Update load step
  m_DeltaU += m_deltaU;
  m_DeltaL += m_deltaL;
}

template <class T>
void gsALMRiks<T>::predictorGuess()
{
  GISMO_ASSERT(m_Uguess.rows()!=0 && m_Uguess.cols()!=0,"Guess is empty");

  m_jacMat = computeJacobian();
  this->factorizeMatrix(m_jacMat);

  m_deltaUt = this->solveSystem(m_forcing);

  T tol = 1e-10;

  if ( ((m_Uguess-m_U).norm() < tol) && ((m_Lguess - m_L) * (m_Lguess - m_L) < tol ) )
  {
    m_note+= "predictor\t";
    T DL = 1.;
    m_deltaUt = this->solveSystem(m_forcing);
    m_deltaU = m_deltaUt / math::sqrt( m_deltaUt.dot(m_deltaUt) + m_DeltaL*DL );
    m_deltaL = DL / math::sqrt( m_deltaUt.dot(m_deltaUt) + m_DeltaL*DL );
  }
  else
  {
    m_deltaL = 1./m_arcLength_prev*(m_Lguess - m_L);
    m_deltaU = 1./m_arcLength_prev*(m_Uguess - m_U);
  }

  // Update iterative step
  m_deltaL *= m_arcLength;
  m_deltaU *= m_arcLength;

  m_DeltaU = m_deltaU;
  m_DeltaL = m_deltaL;

  m_Uguess.resize(0);
}

// template <class T>
// void gsALMRiks<T>::predictor()
// {
//   // Check if the solution on start and prev are similar.
//   // Then compute predictor of the method
//   T tol = 1e-10;

//   if ( ((m_U-m_Uprev).norm() < tol) && ((m_L - m_Lprev) * (m_L - m_Lprev) < tol ) )
//   {
//     m_note+= "predictor\t";
//     T DL = 1.;
//     m_jacMat = computeJacobian();
//     this->factorizeMatrix(m_jacMat);
//     m_deltaUt = this->solveSystem(m_forcing);
//     m_deltaU = m_deltaUt / math::sqrt( m_deltaUt.dot(m_deltaUt) + m_DeltaL*DL );
//     m_deltaL = DL / math::sqrt( m_deltaUt.dot(m_deltaUt) + m_DeltaL*DL );
//   }
//   else
//   {
//     m_deltaL = 1./m_arcLength*(m_L - m_Lprev);
//     m_deltaU = 1./m_arcLength*(m_U - m_Uprev);
//   }

//   // Update iterative step
//   m_deltaL *= m_arcLength;
//   m_deltaU *= m_arcLength;

//   // Update load step
//   m_DeltaU += m_deltaU;
//   m_DeltaL += m_deltaL;
// }

template <class T>
void gsALMRiks<T>::iterationFinish()
{
  m_converged = true;
  m_Uprev = m_U;
  m_Lprev = m_L;
  m_U += m_DeltaU;
  m_L += m_DeltaL;
}

// ------------------------------------------------------------------------------------------------------------
// ---------------------------------------Output functions-----------------------------------------------------
// ------------------------------------------------------------------------------------------------------------

template <class T>
void gsALMRiks<T>::initOutput()
{
  gsInfo<<"\t";
  gsInfo<<std::setw(4)<<std::left<<"It.";
  gsInfo<<std::setw(17)<<std::left<<"Res. F";
  gsInfo<<std::setw(17)<<std::left<<"|dU|/|Du|";
  gsInfo<<std::setw(17)<<std::left<<"dL/DL";
  gsInfo<<std::setw(17)<<std::left<<"|U|";
  gsInfo<<std::setw(17)<<std::left<<"L";
  gsInfo<<std::setw(17)<<std::left<<"|DU|";
  gsInfo<<std::setw(17)<<std::left<<"DL";
  gsInfo<<std::setw(17)<<std::left<<"|dU|";
  gsInfo<<std::setw(17)<<std::left<<"dL";
  gsInfo<<std::setw(17)<<std::left<<"ds²";
  gsInfo<<std::setw(17)<<std::left<<"|dU|²";
  gsInfo<<std::setw(17)<<std::left<<"dL²";
  gsInfo<<std::setw(17)<<std::left<<"Dmin";
  gsInfo<<std::setw(17)<<std::left<<"m_note";
  gsInfo<<"\n";

  m_note = "";
}

template <class T>
void gsALMRiks<T>::stepOutput()
{
  computeStability(false);

  gsInfo<<"\t";
  gsInfo<<std::setw(4)<<std::left<<m_numIterations;
  gsInfo<<std::setw(17)<<std::left<<m_residueF;
  gsInfo<<std::setw(17)<<std::left<<m_residueU;
  gsInfo<<std::setw(17)<<std::left<<m_residueL;
  gsInfo<<std::setw(17)<<std::left<<(m_U+m_DeltaU).norm();
  gsInfo<<std::setw(17)<<std::left<<(m_L + m_DeltaL);
  gsInfo<<std::setw(17)<<std::left<<m_DeltaU.norm();
  gsInfo<<std::setw(17)<<std::left<<m_DeltaL;
  gsInfo<<std::setw(17)<<std::left<<m_deltaU.norm();
  gsInfo<<std::setw(17)<<std::left<<m_deltaL;
  gsInfo<<std::setw(17)<<std::left<<m_phi * math::pow(m_DeltaU.norm(),2.0) + (1.0-m_phi) * math::pow(m_DeltaL,2.0);
  gsInfo<<std::setw(17)<<std::left<<m_phi * math::pow(m_DeltaU.norm(),2.0);
  gsInfo<<std::setw(17)<<std::left<<(1-m_phi) * math::pow(m_DeltaL,2.0);
  gsInfo<<std::setw(17)<<std::left<<m_indicator;
  gsInfo<<std::setw(17)<<std::left<<m_note;
  gsInfo<<"\n";

  m_note = "";
}

} // namespace gismo