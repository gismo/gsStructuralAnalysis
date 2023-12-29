/** @file gsDynamicExplicitEuler.hpp

    @brief Performs the arc length method to solve a nonlinear equation system.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#pragma once

namespace gismo
{

template <class T>
gsStatus gsDynamicExplicitEuler<T>::_step(const T dt)
{
  gsVector<T> R;
  gsSparseMatrix<T> M, C;
  _computeResidual(m_U,m_time,R);
  _computeDamping(m_U,m_time,C);
  _computeMass(m_time,M);

  m_Uold = m_U;
  m_Vold = m_V;
  m_Aold = m_A;

  m_U = m_Uold + m_dt * m_Vold;
  // if (m_lumped)
  // {
  //     gsVector<T> M = _computeMassLumped(m_t);
  //     m_DeltaV = m_dt * M.cwiseProduct(R - C * m_Vold);// element-wise;
  // }
  // else
  // {
      // gsSparseMatrix<T> M = _computeMass(m_t);
      m_solver->compute(M);
      m_V = m_Vold + m_dt * m_solver->solve( R - C * m_Vold );
  // }

  m_time += m_dt;

  return gsStatus::Success;
}

template <class T>
void gsDynamicExplicitEuler<T>::_initOutput()
{
  gsInfo<<"\t";
  gsInfo<<std::setw(4)<<std::left<<"It.";
  gsInfo<<std::setw(17)<<std::left<<"|R|/|R0|";
  gsInfo<<std::setw(17)<<std::left<<"|dU|/|U0|";
}

template <class T>
void gsDynamicExplicitEuler<T>::_stepOutput(const index_t it, const index_t resnorm, const index_t updatenorm)
{
  gsInfo<<"\t";
  gsInfo<<std::setw(4)<<std::left<<it;
  gsInfo<<std::setw(17)<<std::left<<resnorm;
  gsInfo<<std::setw(17)<<std::left<<updatenorm;
}

} // namespace gismo