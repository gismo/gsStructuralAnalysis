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

// template <class T, bool _NL>
// gsStatus gsDynamicExplicitEuler<T,_NL>::_step(const T dt)
// {
//   gsVector<T> R;
//   gsSparseMatrix<T> M, C;
//   _computeResidual(m_U,m_time,R);
//   _computeDamping(m_U,m_time,C);
//   _computeMass(m_time,M);

//   m_Uold = m_U;
//   m_Vold = m_V;
//   m_Aold = m_A;

//   m_U = m_Uold + m_dt * m_Vold;
//   // if (m_lumped)
//   // {
//   //     gsVector<T> M = _computeMassLumped(m_t);
//   //     m_DeltaV = m_dt * M.cwiseProduct(R - C * m_Vold);// element-wise;
//   // }
//   // else
//   // {
//       // gsSparseMatrix<T> M = _computeMass(m_t);
//       m_solver->compute(M);
//       m_V = m_Vold + m_dt * m_solver->solve( R - C * m_Vold );
//   // }

//   m_time += m_dt;

//   return gsStatus::Success;
// }

template <class T, bool _NL>
template <bool _nonlinear>
typename std::enable_if<(_nonlinear==false), gsStatus>::type
gsDynamicExplicitEuler<T,_NL>::_step_impl(const T t, const T dt, gsVector<T> & U, gsVector<T> & V, gsVector<T> & A) const
{
  gsVector<T> Uold = U;
  gsVector<T> Vold = V;
  gsVector<T> Aold = A;

  index_t N = U.rows();
  gsVector<T> sol(2*N);
  sol.topRows(N) = Uold;
  sol.bottomRows(N) = Vold;

  gsVector<T> F;
  gsSparseMatrix<T> M, Minv, C, K;

  // Computed at t=t0
  this->_computeMass(t,M);
  this->_computeMassInverse(M,Minv);
  this->_computeForce(t,F);
  this->_computeDamping(U,t,C);
  this->_computeJacobian(U,t,K);

  this->_initOutput();
  sol.topRows(N) += dt * Vold;
  sol.bottomRows(N) += dt * Minv * (F - K * Uold - C * Vold);
  this->_stepOutput(0,sol.norm(),0.);
  gsDebugVar(sol.transpose());

  U = sol.topRows(N);
  V = sol.bottomRows(N);

  if (math::isinf(sol.norm()) || math::isnan(sol.norm()))
    return gsStatus::NotConverged;
  else
    return gsStatus::Success;
}


template <class T, bool _NL>
template <bool _nonlinear>
typename std::enable_if<(_nonlinear==true), gsStatus>::type
gsDynamicExplicitEuler<T,_NL>::_step_impl(const T t, const T dt, gsVector<T> & U, gsVector<T> & V, gsVector<T> & A) const
{
  gsVector<T> Uold = U;
  gsVector<T> Vold = V;
  gsVector<T> Aold = A;

  index_t N = U.rows();
  gsVector<T> sol(2*N);
  sol.topRows(N) = Uold;
  sol.bottomRows(N) = Vold;

  gsVector<T> R;
  gsSparseMatrix<T> M, Minv, C, K;

  // Computed at t=t0
  this->_computeMass(t,M);
  this->_computeMassInverse(M,Minv);
  this->_computeDamping(Uold,t,C);
  this->_computeResidual(Uold,t,R);

  this->_initOutput();
  sol.topRows(N) += dt * Vold;
  sol.bottomRows(N) += dt * Minv * ( - R - C * Vold);
  this->_stepOutput(0,sol.norm(),0.);

  U = sol.topRows(N);
  V = sol.bottomRows(N);

  if (math::isinf(sol.norm()) || math::isnan(sol.norm()))
    return gsStatus::NotConverged;
  else
    return gsStatus::Success;
}

template <class T, bool _NL>
gsStatus gsDynamicExplicitEuler<T,_NL>::_step(const T t, const T dt,
                                        gsVector<T> & U, gsVector<T> & V,
                                        gsVector<T> & A) const
{
    gsStatus status = gsStatus::NotStarted;
    status = _step_impl<_NL>(t,dt,U,V,A);
    return status;
}

template <class T, bool _NL>
void gsDynamicExplicitEuler<T,_NL>::_initOutput() const
{
  if (m_options.getSwitch("Verbose"))
  {
    gsInfo<<"\t";
    gsInfo<<std::setw(4)<<std::left<<"It.";
    gsInfo<<std::setw(17)<<std::left<<"|R|";
    gsInfo<<"\n";    
  }
}

template <class T, bool _NL>
void gsDynamicExplicitEuler<T,_NL>::_stepOutput(const index_t it, const T resnorm, const T updatenorm) const
{
  if (m_options.getSwitch("Verbose"))
  {
    gsInfo<<"\t";
    gsInfo<<std::setw(4)<<std::left<<it;
    gsInfo<<std::setw(17)<<std::left<<resnorm;
    gsInfo<<"\n";    
  }
}

} // namespace gismo