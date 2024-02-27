/** @file gsDynamicBathe.hpp

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
void gsDynamicBathe<T,_NL>::defaultOptions()
{
  Base::defaultOptions();
  m_options.addReal("alpha","Beta parameter for Bathe's method, such that 0 =< 2 beta =< 1",0.25);
  m_options.addReal("delta","Beta parameter for Bathe's method, such that 0 =< delta =< 1",0.50);
  m_options.addReal("gamma","Beta parameter for Bathe's method, such that 0 =< gamma =< 1",0.50);
}


template <class T, bool _NL>
template <bool _nonlinear>
typename std::enable_if<(_nonlinear==false), gsStatus>::type
gsDynamicBathe<T,_NL>::_step_impl(const T t, const T dt, gsVector<T> & U, gsVector<T> & V, gsVector<T> & A) const
{
  T gamma = m_options.getReal("gamma");

  gsVector<T> Uold = U;
  gsVector<T> Vold = V;
  gsVector<T> Aold = A;

  gsVector<T> Ustep = U;
  gsVector<T> Vstep = V;
  gsVector<T> Astep = A;

  /// Stage 1: A Newmark step with DT=gamma*dt
  this->_stageOutput(1);
  Base::_step(t,gamma*dt,Ustep,Vstep,Astep);

  /// Stage 2: A Newmark step with DT=gamma*dt
  gsVector<T> F;
  gsSparseMatrix<T> M, C, K;

  // Computed at t=t0+dt
  this->_computeMass(t+dt,M);
  this->_computeForce(t+dt,F);
  this->_computeDamping(U,t+dt,C);
  this->_computeJacobian(U,t+dt,K);

  T c1 = (1-gamma)/(gamma*dt);
  T c2 = (-1)/((1-gamma)*gamma*dt);
  T c3 = (2-gamma)/((1-gamma)*dt);

  gsSparseMatrix<T> lhs = K + c3*c3*M + c3*C;
  gsMatrix<T> rhs = F - M*(c1*Vold+c2*Vstep+c1*c3*Uold+c3*c2*Ustep) - C*(c2*Ustep+c1*Uold);

  this->_stageOutput(2);
  this->_initOutput();
  typename gsSparseSolver<T>::uPtr solver;
  solver = gsSparseSolver<T>::get( m_options.getString("Solver") ); 
  solver->compute(lhs);
  U = solver->solve(rhs);
  V = c1*Uold + c2*Ustep + c3*U;
  A = c1*Vold + c2*Vstep + c3*V;
  
  this->_stepOutput(0,(F - K*U - M * A - C * V).norm(),U.norm());
  if (math::isinf(U.norm()) || math::isnan(U.norm()))
    return gsStatus::NotConverged;
  else
    return gsStatus::Success;
}


template <class T, bool _NL>
template <bool _nonlinear>
typename std::enable_if<(_nonlinear==true), gsStatus>::type
gsDynamicBathe<T,_NL>::_step_impl(const T t, const T dt, gsVector<T> & U, gsVector<T> & V, gsVector<T> & A) const
{
  T gamma = m_options.getReal("gamma");

  T c1 = (1-gamma)/(gamma*dt);
  T c2 = (-1)/((1-gamma)*gamma*dt);
  T c3 = (2-gamma)/((1-gamma)*dt);


  gsVector<T> Uold = U;
  gsVector<T> Vold = V;
  gsVector<T> Aold = A;

  gsVector<T> Ustep = U;
  gsVector<T> Vstep = V;
  gsVector<T> Astep = A;

  /// Stage 1: A Newmark step with DT=gamma*dt
  this->_stageOutput(1);
  Base::_step(t,gamma*dt,Ustep,Vstep,Astep);

  /// Stage 2: A Newmark step with DT=gamma*dt
  gsVector<T> R;
  gsSparseMatrix<T> M, C, K;
  gsSparseMatrix<T> lhs;
  gsMatrix<T> rhs;
  gsMatrix<T> dU;

  // Computed at t=t0+dt
  this->_computeMass(t+dt,M);
  this->_computeDamping(U,t+dt,C);
  this->_computeResidual(U,t+dt,R);

  U = Ustep;
  V = Vstep;
  A = Astep;

  rhs = R - M*(c1*Vold+c2*Vstep+c1*c3*Uold+c3*c2*Ustep+c3*c3*U) - C*(c1*Uold+c2*Ustep+c3*U);

  T tolU = m_options.getReal("TolU");
  T tolF = m_options.getReal("TolF");
  T updateNorm   = 10.0*tolU;
  T residualNorm  = rhs.norm();
  T residualNorm0 = (residualNorm!=0) ? residualNorm : residualNorm;
  this->_initOutput();
  T Unorm, dUnorm;
  for (index_t numIterations = 0; numIterations < m_options.getInt("MaxIter"); ++numIterations)
  {
    if ((!m_options.getSwitch("Quasi")) || ((numIterations==0) || (numIterations % m_options.getInt("QuasiIterations") == 0)))
    {
      // Computed at t=t0+dt      
      this->_computeDamping(U,t+dt,C);
      this->_computeJacobian(U,t+dt,K);
    }

    lhs = K + c3*c3*M + c3*C;

    typename gsSparseSolver<T>::uPtr solver;
    solver = gsSparseSolver<T>::get( m_options.getString("Solver") ); 
    solver->compute(lhs);
    // U = solver->solve(rhs);
    // V = delta/( dt*alpha ) * (U-Uold) - Vold;
    // A = 1/( dt*dt*alpha ) * (U-Uold-Vold*dt) - Aold;

    dU = solver->solve(rhs);
    U += dU;

    Unorm = U.norm();
    dUnorm = dU.norm();
    updateNorm = (Unorm!=0) ? dUnorm/Unorm : dUnorm;

    this->_computeResidual(U,t+dt,R);
    rhs = R - M*(c1*Vold+c2*Vstep+c1*c3*Uold+c3*c2*Ustep+c3*c3*U) - C*(c1*Uold+c2*Ustep+c3*U);
    residualNorm = rhs.norm() / residualNorm0;
    updateNorm = dU.norm() / U.norm();

    this->_stepOutput(numIterations,residualNorm,updateNorm);

    if ( (updateNorm<tolU && residualNorm<tolF) )
    {
        V = c1*Uold + c2*Ustep + c3*U;
        A = c1*Vold + c2*Vstep + c3*V;
        return gsStatus::Success;
    }
  }

  V = c1*Uold + c2*Ustep + c3*U;
  A = c1*Vold + c2*Vstep + c3*V;
  gsInfo<<"maximum iterations reached. Solution did not converge\n";
  return gsStatus::NotConverged;
}

template <class T, bool _NL>
void gsDynamicBathe<T,_NL>::_stageOutput(index_t stage) const
{
  if (m_options.getSwitch("Verbose"))
  {
    gsInfo<<"Stage "<<stage<<":";
    gsInfo<<"\n";
  }
}

template <class T, bool _NL>
void gsDynamicBathe<T,_NL>::_initOutput() const
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
void gsDynamicBathe<T,_NL>::_stepOutput(const index_t it, const T resnorm, const T updatenorm) const
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