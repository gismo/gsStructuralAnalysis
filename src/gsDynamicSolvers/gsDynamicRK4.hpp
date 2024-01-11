/** @file gsDynamicRK4.hpp

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
// gsStatus gsDynamicRK4<T,_NL>::_step(const T dt)
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
gsDynamicRK4<T,_NL>::_step_impl(const T t, const T dt, gsVector<T> & U, gsVector<T> & V, gsVector<T> & A) const
{
  /*
  */
  gsVector<T> Uold = U;
  gsVector<T> Vold = V;
  gsVector<T> Aold = A;

  index_t N = U.rows();
  gsVector<T> sol(2*N);
  sol.topRows(N) = Uold;
  sol.bottomRows(N) = Vold;

  gsVector<T> F, R; //R is the residual vector, and the new _computeForce give the inline factor
  gsSparseMatrix<T> M, Minv, C, K;

  // Computed at t=t0
  this->_computeMass(t,M);
  this->_computeMassInverse(M,Minv);
  // this->_computeForce(t,F);
  this->_computeDamping(U,t,C); //C is damping
  this->_computeJacobian(U,t,K);

  // this->_initOutput();
  // Initialize parameters for RK4
  gsVector<T> k1(2*N), k2(2*N), k3(2*N), k4(2*N);
  gsVector<T> Utmp, Vtmp;

  //Step1 (calculate k1)
  _computeForce(t, F);
  R = F - K * Uold;
  k1.topRows(N) = Vold;
  k1.bottomRows(N) = Minv * (R - C * Vold);

  //Step2 (calculate k2)
  Utmp = Uold + dt/2. * k1.topRows(N);
  Vtmp = Vold + dt/2. * k1.bottomRows(N);
  _computeForce(t + dt/2.,F);
  R = F - K * Utmp;
  k2.topRows(N) = Vtmp;
  k2.bottomRows(N) = Minv * ( R - C * Vtmp);

  //Step3 (calculate k3)
  Utmp = Uold + m_dt/2. * k2.topRows(N);
  Vtmp = Vold + m_dt/2. * k2.bottomRows(N);
  R = _computeForce(t + dt/2.) - K * Utmp;
  k3.topRows(N) = Vtmp;
  k3.bottomRows(N) = Minv * ( R - C * Vtmp);

  //Step4 (calculate k4)
  Utmp = Uold + m_dt/2. * k3.topRows(N);
  Vtmp = Vold + m_dt/2. * k3.bottomRows(N);
  R = _computeForce(t + dt/2.) - K * Utmp;
  k4.topRows(N) = Vtmp;
  k4.bottomRows(N) = Minv * ( R - C * Vtmp);

  sol += 1./6 * dt * (k1 + 2.*k2 + 2.*k3 + k4);

  U = std::move(sol.topRows(N)); //Use std move to update vectors to avoid unnecessary copying
  V = std::move(sol.bottomRows(N));

  if (math::isinf(sol.norm()) || math::isnan(sol.norm()))
    return gsStatus::NotConverged;
  else
    return gsStatus::Success;
}


template <class T, bool _NL>
template <bool _nonlinear>
typename std::enable_if<(_nonlinear==true), gsStatus>::type
gsDynamicRK4<T,_NL>::_step_impl(const T t, const T dt, gsVector<T> & U, gsVector<T> & V, gsVector<T> & A) const
{
  /*
  */
  gsVector<T> Uold = U;
  gsVector<T> Vold = V;
  gsVector<T> Aold = A;

  index_t N = U.rows();
  gsVector<T> sol(2*N);
  sol.topRows(N) = Uold;
  sol.bottomRows(N) = Vold;

  gsVector<T> F, R; //R is the residual vector, and the new _computeForce give the inline factor
  gsSparseMatrix<T> M, Minv, C, K;

  // Computed at t=t0
  this->_computeMass(t,M);
  this->_computeMassInverse(M,Minv);
  // this->_computeForce(t,F);
  this->_computeDamping(U,t,C); //C is damping
  this->_computeJacobian(U,t,K);

  // this->_initOutput();
  // Initialize parameters for RK4
  gsVector<T> k1(2*N), k2(2*N), k3(2*N), k4(2*N);
  gsVector<T> Utmp, Vtmp;

  //Step1 (calculate k1)
  R = _computeResidual(Uold, t);
  k1.topRows(N) = Vold;
  k1.bottomRows(N) = Minv * (R - C * Vold);

  //Step2 (calculate k2)
  Utmp = Uold + dt/2. * k1.topRows(N);
  Vtmp = Vold + dt/2. * k1.bottomRows(N);
  R =_computeResidual(Utmp,t + dt/2.);
  k2.topRows(N) = Vtmp;
  k2.bottomRows(N) = Minv * ( R - C * Vtmp);

  //Step3 (calculate k3)
  Utmp = Uold + m_dt/2. * k2.topRows(N);
  Vtmp = Vold + m_dt/2. * k2.bottomRows(N);
  R =_computeResidual(Utmp,t + dt/2.);
  k3.topRows(N) = Vtmp;
  k3.bottomRows(N) = Minv * ( R - C * Vtmp);

  //Step4 (calculate k4)
  Utmp = Uold + m_dt/2. * k3.topRows(N);
  Vtmp = Vold + m_dt/2. * k3.bottomRows(N);
  R =_computeResidual(Utmp,t + dt/2.);
  k4.topRows(N) = Vtmp;
  k4.bottomRows(N) = Minv * ( R - C * Vtmp);

  sol += 1./6 * dt * (k1 + 2.*k2 + 2.*k3 + k4);

  U = std::move(sol.topRows(N)); //Use std move to update vectors to avoid unnecessary copying
  V = std::move(sol.bottomRows(N));

  if (math::isinf(sol.norm()) || math::isnan(sol.norm()))
    return gsStatus::NotConverged;
  else
    return gsStatus::Success;
}
}

template <class T, bool _NL>
void gsDynamicRK4<T,_NL>::_initOutput() const
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
void gsDynamicRK4<T,_NL>::_stepOutput(const index_t it, const T resnorm, const T updatenorm) const
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