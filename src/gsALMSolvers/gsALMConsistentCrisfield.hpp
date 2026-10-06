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

  m_phi0    = 0;
  m_phi0Set = false;
}

// Constraint |DeltaU|^2 + phi^2 DeltaL^2 |F|^2 = s^2 (see iteration()).
//
// With automatic scaling, phi is one constant per run, phi0 = |u_t,0| / |F|
// with u_t,0 = K^{-1} F at the start point of the run. The load term
// phi |F| DeltaL = |u_t,0| DeltaL is the displacement of the initial tangent
// for the load increment DeltaL, so both terms of the constraint have the units
// of a displacement, and the weight psi = phi^2 |F|^2 of the load term is the
// same at every step.
//
// A scaling that follows the path changes the metric of the constraint from
// step to step. The tangent ratio |u_t|/|F| is unbounded at a limit point,
// where K is singular. The secant ratio |U| / (|lambda| |F|) is unbounded as
// lambda -> 0 at a regular point, which bounds |DeltaL| <= s |lambda| / |U|, so
// that lambda = 0 is approached only geometrically. A constant phi has neither
// singularity.
//
// A run starts at initialize(). Its phi0 is taken from u_t at the first step of
// the run, wherever that step starts (the origin, or a restart through
// setSolution), and is not changed by later steps, option calls or by putting
// the solver back at an earlier converged point.
template <class T>
void gsALMConsistentCrisfield<T>::setAutoScaling(const bool utSolved)
{
  if (m_phi_user)
    return;

  if (!m_phi0Set)
  {
    if (!utSolved)
      m_deltaUt = this->solveSystem(m_forcing);
    m_phi0    = math::pow( m_deltaUt.dot(m_deltaUt) / m_forcing.dot(m_forcing),0.5);
    m_phi0Set = true;
  }
  m_phi = m_phi0;
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
T gsALMConsistentCrisfield<T>::constraintResidual() const
{
  // Spherical constraint |DeltaU|^2 + phi^2*DeltaL^2*|F|^2 = s^2, relative to s^2
  const T s2 = m_arcLength*m_arcLength;
  return math::abs(m_DeltaU.dot(m_DeltaU) + m_phi*m_phi*m_DeltaL*m_DeltaL*m_forcing.dot(m_forcing) - s2) / s2;
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

    // Automatic scaling phi: constant per run, see setAutoScaling()
    this->setAutoScaling(true);
  }
  else
  {
    m_deltaL = 1./m_arcLength_prev*(m_L - m_Lprev);
    m_deltaU = 1./m_arcLength_prev*(m_U - m_Uprev);

    // Automatic scaling phi: constant per run, see setAutoScaling()
    this->setAutoScaling(false);
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

    // Automatic scaling phi: constant per run, see setAutoScaling()
    this->setAutoScaling(true);
  }

  else
  {
    m_deltaL = 1./m_arcLength_prev*(m_Lguess - m_L);
    m_deltaU = 1./m_arcLength_prev*(m_Uguess - m_U);

    // Automatic scaling phi: constant per run, see setAutoScaling()
    this->setAutoScaling(false);
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