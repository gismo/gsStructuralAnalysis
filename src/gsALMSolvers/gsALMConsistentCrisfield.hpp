/** @file gsALMConsistentCrisfield.hpp

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
void gsALMConsistentCrisfield<T>::defaultOptions()
{
    Base::defaultOptions();
    m_options.addReal("Scaling","Set Scaling factor Phi",-1);
}

template <class T>
void gsALMConsistentCrisfield<T>::getOptions()
{
    Base::getOptions();
    m_phi                 = m_options.getReal("Scaling");
    m_phi_user = m_phi == -1 ? false : true;
}

/// What to initialize for this method
template <class T>
void gsALMConsistentCrisfield<T>::initMethods()
{
  m_numDof = m_forcing.size();
  m_DeltaU = m_U = gsVector<T>::Zero(m_numDof);
  m_DeltaL = m_L = 0.0;

  m_Uprev = gsVector<T>::Zero(m_numDof);
  m_Lprev = 0.0;
}

// ------------------------------------------------------------------------------------------------------------
// ---------------------------------------Consistent Crisfield's method----------------------------------------
// ------------------------------------------------------------------------------------------------------------

/// Define what is computed for quasi newton iterations before iterations start
template <class T>
void gsALMConsistentCrisfield<T>::quasiNewtonPredictor()
{
  m_jacMat = computeJacobian();
  this->factorizeMatrix(m_jacMat);
  computeUt(); // rhs does not depend on solution
  computeUbar(); // rhs contains residual and should be computed every time

}

/// Define what is computed it no quasi newton iterations during iterations
template <class T>
void gsALMConsistentCrisfield<T>::quasiNewtonIteration()
{
  m_jacMat = computeJacobian();
  this->factorizeMatrix(m_jacMat);
  computeUt(); // rhs does not depend on solution
}


template <class T>
void gsALMConsistentCrisfield<T>::iteration() // see Carrera1994 eqs 20 and 21
{
  computeUbar(); // rhs contains residual and should be computed every time

  T cres = m_DeltaU.dot(m_DeltaU) + m_phi*m_phi*m_DeltaL*m_DeltaL * m_forcing.dot(m_forcing) - m_arcLength*m_arcLength;
  T num = cres + (2*m_DeltaU).dot(m_deltaUbar);
  T denum = (2*m_DeltaU).dot(m_deltaUt) + m_phi*m_phi*2*m_DeltaL*m_forcing.dot(m_forcing);
  m_deltaL = - num / denum;
  m_deltaU = m_deltaL * m_deltaUt + m_deltaUbar;

  m_DeltaU += m_deltaU;
  m_DeltaL += m_deltaL;
}

template <class T>
void gsALMConsistentCrisfield<T>::initiateStep()
{
  // (m_U,m_L) is the present solution (iteratively updated)
  // (m_Uprev,m_Lprev) is the previously converged solution before (m_Lprev,m_Uprev)

  // Reset step
  m_DeltaU = m_deltaU =  gsVector<T>::Zero(m_numDof);
  m_DeltaL = m_deltaL = 0.0;

}

template <class T>
void gsALMConsistentCrisfield<T>::predictor()
{
  m_jacMat = computeJacobian();
  this->factorizeMatrix(m_jacMat);

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

    if (!m_phi_user)
      m_phi = math::pow( m_deltaUt.dot(m_deltaUt) / m_forcing.dot(m_forcing),0.5);
  }
  else
  {
    m_deltaL = 1./m_arcLength_prev*(m_L - m_Lprev);
    m_deltaU = 1./m_arcLength_prev*(m_U - m_Uprev);

    if (!m_phi_user)
      m_phi = math::pow( m_deltaUt.dot(m_deltaUt) / m_forcing.dot(m_forcing),0.5);
  }

  // Update iterative step
  m_deltaL *= m_arcLength;
  m_deltaU *= m_arcLength;

  // Update load step
  m_DeltaU += m_deltaU;
  m_DeltaL += m_deltaL;
}

template <class T>
void gsALMConsistentCrisfield<T>::predictorGuess()
{
  m_jacMat = computeJacobian();
  this->factorizeMatrix(m_jacMat);

  // Check if the solution on start and prev are similar.
  // Then compute predictor of the method
  T tol = 1e-10;

  if ( ((m_Uguess-m_U).norm() < tol) && ((m_Lguess - m_L) * (m_Lguess - m_L) < tol ) )
  {
    m_note+= "predictor\t";
    T DL = 1.;
    m_deltaUt = this->solveSystem(m_forcing);
    m_deltaU = m_deltaUt / math::sqrt( m_deltaUt.dot(m_deltaUt) + m_DeltaL*DL );
    m_deltaL = DL / math::sqrt( m_deltaUt.dot(m_deltaUt) + m_DeltaL*DL );

    if (!m_phi_user)
      m_phi = math::pow( m_deltaUt.dot(m_deltaUt) / m_forcing.dot(m_forcing),0.5);
  }

  else
  {
    m_deltaL = 1./m_arcLength_prev*(m_Lguess - m_L);
    m_deltaU = 1./m_arcLength_prev*(m_Uguess - m_U);

    if (!m_phi_user)
      m_phi = math::pow( m_deltaUt.dot(m_deltaUt) / m_forcing.dot(m_forcing),0.5);
  }

  // Update iterative step
  m_deltaL *= m_arcLength;
  m_deltaU *= m_arcLength;

  // Update load step
  m_DeltaU += m_deltaU;
  m_DeltaL += m_deltaL;
}

template <class T>
void gsALMConsistentCrisfield<T>::iterationFinish()
{
  m_converged = true;
  m_Uprev = m_U;
  m_Lprev = m_L;
  m_U += m_DeltaU;
  m_L += m_DeltaL;
  m_DeltaUold = m_DeltaU;
  m_DeltaLold = m_DeltaL;
}


// ------------------------------------------------------------------------------------------------------------
// ---------------------------------------Output functions-----------------------------------------------------
// ------------------------------------------------------------------------------------------------------------

template <class T>
void gsALMConsistentCrisfield<T>::initOutput()
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
void gsALMConsistentCrisfield<T>::stepOutput()
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