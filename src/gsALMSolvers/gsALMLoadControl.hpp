/** @file gsALMLoadControl.hpp

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
void gsALMLoadControl<T>::quasiNewtonPredictor()
{
  m_jacMat = computeJacobian();
  this->factorizeMatrix(m_jacMat);
  computeUbar(); // rhs contains residual and should be computed every time

}

template <class T>
void gsALMLoadControl<T>::quasiNewtonIteration()
{
  m_jacMat = computeJacobian();
  this->factorizeMatrix(m_jacMat);
}

template <class T>
void gsALMLoadControl<T>::initMethods()
{
  m_numDof = m_forcing.size();
  m_DeltaU = m_U = gsVector<T>::Zero(m_numDof);
  m_DeltaL = m_L = 0.0;
}

// ------------------------------------------------------------------------------------------------------------
// ---------------------------------------Load Control method--------------------------------------------------------
// ------------------------------------------------------------------------------------------------------------

template <class T>
void gsALMLoadControl<T>::iteration()
{
  computeUbar(); // rhs contains residual and should be computed every time

  // Compute next solution
  m_deltaU = m_deltaUbar;
  m_DeltaU += m_deltaU;
}

template <class T>
void gsALMLoadControl<T>::initiateStep()
{
  // (m_U,m_L) is the present solution (iteratively updated)
  // (m_Uprev,m_Lprev) is the previously converged solution before (m_Lprev,m_Uprev)

  // Reset step
  m_DeltaU = m_deltaU = gsVector<T>::Zero(m_numDof);
  m_DeltaL = m_deltaL = 0.0;
}

template <class T>
void gsALMLoadControl<T>::predictor()
{
  m_jacMat = computeJacobian();
  this->factorizeMatrix(m_jacMat);

  m_DeltaL = m_deltaL = m_arcLength;
  m_deltaUt = this->solveSystem(m_forcing);
  m_DeltaU = m_deltaL*m_deltaUt;

  // gsDebugVar(m_DeltaU.norm());
  m_note+= "predictor\t";
}

template <class T>
void gsALMLoadControl<T>::predictorGuess()
{
  m_jacMat = computeJacobian();
  this->factorizeMatrix(m_jacMat);

  m_DeltaL = m_deltaL = m_Lguess - m_L;
  m_deltaUt = this->solveSystem(m_forcing);
  m_DeltaU = m_deltaL*m_deltaUt;

  // gsDebugVar(m_DeltaU.norm());
  m_note+= "predictor\t";
}

template <class T>
void gsALMLoadControl<T>::iterationFinish()
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
void gsALMLoadControl<T>::initOutput()
{
  gsInfo<<"\t";
  gsInfo<<std::setw(4)<<std::left<<"It.";
  gsInfo<<std::setw(17)<<std::left<<"Res. F";
  gsInfo<<std::setw(17)<<std::left<<"|dU|/|Du|";
  gsInfo<<std::setw(17)<<std::left<<"|U|";
  gsInfo<<std::setw(17)<<std::left<<"L";
  gsInfo<<std::setw(17)<<std::left<<"|DU|";
  gsInfo<<std::setw(17)<<std::left<<"DL";
  gsInfo<<std::setw(17)<<std::left<<"|dU|";
  gsInfo<<std::setw(17)<<std::left<<"dL";
  gsInfo<<std::setw(17)<<std::left<<"Dmin";
  gsInfo<<std::setw(17)<<std::left<<"m_note";
  gsInfo<<"\n";

  m_note = "";
}

template <class T>
void gsALMLoadControl<T>::stepOutput()
{
  computeStability(false);

  gsInfo<<"\t";
  gsInfo<<std::setw(4)<<std::left<<m_numIterations;
  gsInfo<<std::setw(17)<<std::left<<m_residueF;
  gsInfo<<std::setw(17)<<std::left<<m_residueU;
  gsInfo<<std::setw(17)<<std::left<<(m_U+m_DeltaU).norm();
  gsInfo<<std::setw(17)<<std::left<<(m_L + m_DeltaL);
  gsInfo<<std::setw(17)<<std::left<<m_DeltaU.norm();
  gsInfo<<std::setw(17)<<std::left<<m_DeltaL;
  gsInfo<<std::setw(17)<<std::left<<m_deltaU.norm();
  gsInfo<<std::setw(17)<<std::left<<m_deltaL;
  gsInfo<<std::setw(17)<<std::left<<m_indicator;
  gsInfo<<std::setw(17)<<std::left<<m_note;
  gsInfo<<"\n";

  m_note = "";
}

} // namespace gismo