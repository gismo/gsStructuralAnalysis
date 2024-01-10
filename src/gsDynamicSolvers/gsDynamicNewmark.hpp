/** @file gsDynamicNewmark.hpp

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

template <class T, bool _NL>
void gsDynamicNewmark<T,_NL>::defaultOptions()
{
  Base::defaultOptions();
  m_options.addReal("alpha","Beta parameter for Newmark's method, such that 0 =< 2 beta =< 1",0.25);
  m_options.addReal("delta","Beta parameter for Newmark's method, such that 0 =< gamma =< 1",0.50);
}


template <class T, bool _NL>
template <bool _nonlinear>
typename std::enable_if<(_nonlinear==false), gsStatus>::type
gsDynamicNewmark<T,_NL>::_step_impl(const T t, const T dt, gsVector<T> & U, gsVector<T> & V, gsVector<T> & A) const
{
  gsVector<T> Uold = U;
  gsVector<T> Vold = V;

  index_t N = U.rows();
  gsVector<T> &sol(2*N);
  sol.topRows(N) = Uold;
  sol.bottomRows(N) = Vold;

  gsVector<T> F, R;
  gsSparseMatrix<T> M, Minv, C, K;

  // Computed at t=t0
  this->_computeMass(t,M);
  this->_computeMassInverse(M,Minv);
  this->_computeForce(t,F);
  this->_computeDamping(U,t,C);
  this->_computeJacobian(U,t,K);

  gsVector<T> k1(2*N), k2(2*N), k3(2*N), k4(2*N);
  gsVector<T> Utmp, Vtmp;

  // Stage 1
  R = _computeForce(t) - K * Uold;
  k1.topRows(N) = Vold;
  k1.bottomRows(N) = Minv * ( R - C * Vold);

  // Stage 2
  // NOTE: Compute K, C, M on t+dt/2
  Utmp = Uold + dt/2. * k1.topRows(N);
  Vtmp = Vold + dt/2. * k1.bottomRows(N);
  R = _computeForce(t + dt/2.) - K * Utmp;
  k2.topRows(N) = Vtmp;
  k2.bottomRows(N) = Minv * ( R - C * Vtmp);

  // Stage 3
  Utmp = Uold + dt/2. * k2.topRows(N);
  Vtmp = Vold + dt/2. * k2.bottomRows(N);
  R = _computeForce(t + dt/2.) - K * Utmp;
  k3.topRows(N) = Vtmp;
  k3.bottomRows(N) = Minv * ( R - C * Vtmp);

  // Stage 4
  Utmp = Uold + dt * k3.topRows(N);
  Vtmp = Vold + dt * k3.bottomRows(N);
  R = _computeForce(t + dt) - K * Utmp;
  k4.topRows(N) = Vtmp;
  k4.bottomRows(N) = Minv * ( R - C * Vtmp);

  sol += 1./6. * dt *  ( k1 + 2.*k2 + 2.*k3 + k4 );
  U = sol.topRows(N);
  V = sol.bottomRows(N);

  if (math::isinf(U.norm()) || math::isnan(U.norm()))
    return gsStatus::NotConverged;
  else
    return gsStatus::Success;
}


template <class T, bool _NL>
template <bool _nonlinear>
typename std::enable_if<(_nonlinear==true), gsStatus>::type
gsDynamicNewmark<T,_NL>::_step_impl(const T t, const T dt, gsVector<T> & U, gsVector<T> & V, gsVector<T> & A) const
{
  //https://ethz.ch/content/dam/ethz/special-interest/baug/ibk/structural-mechanics-dam/education/femII/presentation_05_dynamics_v3.pdf
  gsVector<T> Uold = U;
  gsVector<T> Vold = V;
  gsVector<T> Aold = A;

  gsVector<T> R;
  gsSparseMatrix<T> M, C, K;

  T alpha = m_options.getReal("alpha");
  T delta = m_options.getReal("delta");

  gsSparseMatrix<T> lhs;
  A.setZero();
  U = Uold + dt*Vold + Aold*(0.5 - alpha)*dt*dt + alpha*dt*dt*A;
  V = Vold + Aold*(1 - delta)*dt + delta*dt*A;
  // Computed at t=t0+dt
  this->_computeMass(t+dt,M);
  this->_computeDamping(U,t+dt,C);
  this->_computeResidual(U,t+dt,R);

  gsMatrix<T> rhs = R - C * V - M*(A);

  T tolU = m_options.getReal("TolU");
  T tolF = m_options.getReal("TolF");
  T updateNorm   = 10.0*tolU;
  T residualNorm  = rhs.norm();
  T residualNorm0 = residualNorm;
  gsVector<T> dA;
  this->_initOutput();
  typename gsSparseSolver<T>::uPtr solver;  
  for (index_t numIterations = 0; numIterations < m_options.getInt("MaxIter"); ++numIterations)
  {
    if ((!m_options.getSwitch("Quasi")) || ((numIterations==0) || (numIterations % m_options.getInt("QuasiIterations") == 0)))
    {
      // Computed at t=t0+dt
      this->_computeDamping(U,t+dt,C);
      this->_computeJacobian(U,t+dt,K);
    }

    lhs = M + delta*dt*C + dt*dt*alpha*K;

    solver = gsSparseSolver<T>::get( m_options.getString("Solver") ); 
    solver->compute(lhs);

    dA = solver->solve(rhs);
    A += dA;
    V += dA*delta*dt;
    U += dA*alpha*dt*dt;

    this->_computeResidual(U,t+dt,R);
    rhs = R - C*V - M*A;
    residualNorm = rhs.norm() / residualNorm0;
    updateNorm = dA.norm() / A.norm();

    this->_stepOutput(numIterations,residualNorm,updateNorm);

    if ( (updateNorm<tolU && residualNorm<tolF) )
    {
        return gsStatus::Success;
    }
  }

  gsInfo<<"maximum iterations reached. Solution did not converge\n";
  return gsStatus::NotConverged;
}

template <class T, bool _NL>
void gsDynamicNewmark<T,_NL>::_initOutput() const
{
  if (m_options.getSwitch("Verbose"))
  {
    gsInfo<<"\t";
    gsInfo<<std::setw(4)<<std::left<<"It.";
    gsInfo<<std::setw(17)<<std::left<<"|R|/|R0|";
    gsInfo<<std::setw(17)<<std::left<<"|dU|/|U0|"<<"\n";
  }
}

template <class T, bool _NL>
void gsDynamicNewmark<T,_NL>::_stepOutput(const index_t it, const T resnorm, const T updatenorm) const
{
  if (m_options.getSwitch("Verbose"))
  {
    gsInfo<<"\t";
    gsInfo<<std::setw(4)<<std::left<<it;
    gsInfo<<std::setw(17)<<std::left<<resnorm;
    gsInfo<<std::setw(17)<<std::left<<updatenorm<<"\n";
  }
}

} // namespace gismo