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
  computeUt(); // rhs is state-independent for a dead load; with a forcing callback set it IS state-dependent (see gsALMBase::setForcingFunction)
  computeUbar(); // rhs contains residual and should be computed every time

}

/// Define what is computed it no quasi newton iterations during iterations
template <class T>
void gsALMConsistentCrisfield<T>::quasiNewtonIteration()
{
  m_jacMat = computeJacobian();
  this->factorizeMatrix(m_jacMat);
  computeUt(); // rhs is state-independent for a dead load; with a forcing callback set it IS state-dependent (see gsALMBase::setForcingFunction)
}


template <class T>
void gsALMConsistentCrisfield<T>::iteration() // see Carrera1994 eqs 20 and 21
{
  computeUbar(); // rhs contains residual and should be computed every time

  // Constraint metric: the ||f||^2 scaling is FROZEN per step (see
  // gsALMBase::stepForcing), so cres and its derivative denum refer to one and the same
  // fixed quadratic constraint throughout this step's corrector. The consistent
  // (state-dependent) forcing enters instead through m_deltaUt, see computeUt().
  const gsVector<T> & Fref = this->stepForcing();
  T cres = m_DeltaU.dot(m_DeltaU) + m_phi*m_phi*m_DeltaL*m_DeltaL * Fref.dot(Fref) - m_arcLength*m_arcLength;
  T num = cres + (2*m_DeltaU).dot(m_deltaUbar);
  T denum = (2*m_DeltaU).dot(m_deltaUt) + m_phi*m_phi*2*m_DeltaL*Fref.dot(Fref);
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

  // Freeze the constraint-metric forcing at the predictor state (m_U,m_L): the ||f||^2
  // scaling of the constraint (and phi) must stay fixed during the corrector. A no-op
  // without a forcing callback.
  this->_freezeStepForcing();
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
    // Predictor tangent at the Jacobian's linearization state (m_U,m_L).
    m_deltaUt = this->solveSystem(this->computeForcing(m_U,m_L));
    m_deltaU = m_deltaUt / math::sqrt( m_deltaUt.dot(m_deltaUt) + m_DeltaL*DL );
    m_deltaL = DL / math::sqrt( m_deltaUt.dot(m_deltaUt) + m_DeltaL*DL );

    // Constraint-metric scaling: frozen step forcing (see gsALMBase::stepForcing)
    if (!m_phi_user)
      m_phi = math::pow( m_deltaUt.dot(m_deltaUt) / this->stepForcing().dot(this->stepForcing()),0.5);
  }
  else
  {
    m_deltaL = 1./m_arcLength_prev*(m_L - m_Lprev);
    m_deltaU = 1./m_arcLength_prev*(m_U - m_Uprev);

    // Constraint-metric scaling: frozen step forcing (see gsALMBase::stepForcing).
    // m_deltaUt is NOT computed on this branch, so the historic phi read whatever the
    // PREVIOUS step's corrector left behind -- and on a solver restarted through
    // setSolution()/setLength()/setPrevious() it has never been sized at all (size 0 =>
    // phi = 0). Recompute the predictor tangent at (m_U,m_L), the same state the
    // fresh-start branch linearises at, so that phi is a function of the injected state
    // alone and the restart contract holds. The Jacobian is already factorized above, so
    // this costs one back-substitution.
    // Scoped under !m_phi_user: with a user Scaling nothing here reads m_deltaUt, and
    // keeping the branch untouched in that case makes the change bit-neutral for every
    // in-tree driver (all of which pin Scaling = 0).
    if (!m_phi_user)
    {
      m_deltaUt = this->solveSystem(this->computeForcing(m_U,m_L));
      m_phi = math::pow( m_deltaUt.dot(m_deltaUt) / this->stepForcing().dot(this->stepForcing()),0.5);
    }
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
    // Predictor tangent at the Jacobian's linearization state (m_U,m_L) - NOT at the
    // guess, which only supplies a secant direction.
    m_deltaUt = this->solveSystem(this->computeForcing(m_U,m_L));
    m_deltaU = m_deltaUt / math::sqrt( m_deltaUt.dot(m_deltaUt) + m_DeltaL*DL );
    m_deltaL = DL / math::sqrt( m_deltaUt.dot(m_deltaUt) + m_DeltaL*DL );

    // Constraint-metric scaling: frozen step forcing (see gsALMBase::stepForcing)
    if (!m_phi_user)
      m_phi = math::pow( m_deltaUt.dot(m_deltaUt) / this->stepForcing().dot(this->stepForcing()),0.5);
  }

  else
  {
    m_deltaL = 1./m_arcLength_prev*(m_Lguess - m_L);
    m_deltaU = 1./m_arcLength_prev*(m_Uguess - m_U);

    // Constraint-metric scaling: frozen step forcing (see gsALMBase::stepForcing).
    // Same defect and same fix as the secant branch of predictor() above: m_deltaUt is
    // not computed on this branch, so recompute the predictor tangent at (m_U,m_L) --
    // NOT at the guess, which only supplies a secant direction -- to make phi a function
    // of the injected state alone. Scoped under !m_phi_user, so it is bit-neutral for a
    // user Scaling. This branch has no unit-test coverage.
    if (!m_phi_user)
    {
      m_deltaUt = this->solveSystem(this->computeForcing(m_U,m_L));
      m_phi = math::pow( m_deltaUt.dot(m_deltaUt) / this->stepForcing().dot(this->stepForcing()),0.5);
    }
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
  // "ds", not "ds²": the column below prints distance(DeltaU,DeltaL), which is a SQUARE
  // ROOT. The two columns after it remain squares, so the constraint identity reads
  // ds² = |dU|² + dL², not ds = |dU|² + dL². Header corrected rather than the value, to
  // keep this class's three columns identical to gsALMCrisfield::stepOutput()'s, as the
  // comment in stepOutput() asserts.
  gsInfo<<std::setw(17)<<std::left<<"ds";
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

  // Output only; uses the same (frozen) constraint metric as distance() -- i.e. exactly the
  // cres/denum pair of iteration(), |DeltaU|^2 + A0*DeltaLambda^2 = ds^2 with
  // A0 = phi^2*|f|^2. Same three columns, same functional form, as
  // gsALMCrisfield::stepOutput().
  //
  // ⚠ these three columns
  // used to print gsALMRiks's CONVEX form phi*|DeltaU|^2 + (1-phi)*DeltaLambda^2. That form
  // belongs to a different constraint and to a different meaning of phi -- here phi is the
  // metric scaling psi of the Scaling option, not a convex weight in (0,1) -- so at any
  // psi > 1 the printed "dL^2" column went NEGATIVE and the "ds^2" column bore no relation to
  // the constraint this class actually enforces.
  T A0 = math::pow(m_phi,2)*this->stepForcing().dot(this->stepForcing());

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
  gsInfo<<std::setw(17)<<std::left<<this->distance(m_DeltaU,m_DeltaL);
  gsInfo<<std::setw(17)<<std::left<<math::pow(m_DeltaU.norm(),2.0);
  gsInfo<<std::setw(17)<<std::left<<A0*math::pow(m_DeltaL,2.0);
  gsInfo<<std::setw(17)<<std::left<<m_indicator;
  gsInfo<<std::setw(17)<<std::left<<m_note;
  gsInfo<<"\n";

  m_note = "";
}

} // namespace gismo