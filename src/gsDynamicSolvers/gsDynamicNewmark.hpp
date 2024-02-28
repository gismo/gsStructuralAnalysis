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
  gsVector<T> Aold = A;

  gsVector<T> F;
  gsSparseMatrix<T> M, C, K;

  T alpha = m_options.getReal("alpha");
  T delta = m_options.getReal("delta");

  // Computed at t=t0+dt
  this->_computeMass(t+dt,M);
  this->_computeForce(t+dt,F);

  gsSparseMatrix<T> lhs;
  gsMatrix<T> rhs;
  // predictors
  U = Uold + dt*Vold + Aold*(0.5 - alpha)*dt*dt;
  V = Vold + Aold*(1 - delta)*dt;

  this->_computeDamping(U,t+dt,C);
  this->_computeJacobian(U,t+dt,K);

  // lhs and rhs
  lhs = M + delta*dt*C + dt*dt*alpha*K;
  rhs = F - K*U - C * V;

  this->_initOutput();
  typename gsSparseSolver<T>::uPtr solver;
  solver = gsSparseSolver<T>::get( m_options.getString("Solver") ); 
  solver->compute(lhs);

  A = solver->solve(rhs);
  V += A*delta*dt;
  U += A*alpha*dt*dt;
  
  this->_stepOutput(0,(F - K*U - M * A - C * V).norm(),U.norm());
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
  T residualNorm0 = (residualNorm!=0) ? residualNorm : 1;
  gsVector<T> dA;
  T Anorm, dAnorm;
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

    Anorm = A.norm();
    dAnorm = dA.norm();
    updateNorm = (Anorm!=0) ? dAnorm/Anorm : dAnorm;

    this->_computeResidual(U,t+dt,R);
    rhs = R - C*V - M*A;
    residualNorm = rhs.norm() / residualNorm0;

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