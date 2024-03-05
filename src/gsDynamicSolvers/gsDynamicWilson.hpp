/** @file gsDynamicWilson.hpp

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
void gsDynamicWilson<T,_NL>::defaultOptions()
{
  Base::defaultOptions();
  m_options.addReal("alpha","Beta parameter for Wilson's method, such that 0 =< 2 beta =< 1",0.25);
  m_options.addReal("delta","Beta parameter for Wilson's method, such that 0 =< delta =< 1",0.50);
  m_options.addReal("gamma","Beta parameter for Wilson's method, such that 0 =< gamma =< 1",1.4);
}


template <class T, bool _NL>
template <bool _nonlinear>
typename std::enable_if<(_nonlinear==false), gsStatus>::type
gsDynamicWilson<T,_NL>::_step_impl(const T t, const T dt, gsVector<T> & U, gsVector<T> & V, gsVector<T> & A) const
{
  T alpha = m_options.getReal("alpha");
  T delta = m_options.getReal("delta");
  T gamma = m_options.getReal("gamma");

  gsVector<T> Uold = U;
  gsVector<T> Vold = V;
  gsVector<T> Aold = A;

  gsVector<T> Utmp = U;
  gsVector<T> Vtmp = V;
  gsVector<T> Atmp = A;

  // Perform a Newmark step with DT = gamma*dt
  Base::_step(t,gamma*dt,Utmp,Vtmp,Atmp);

  // Compute solutions
  A = (1-1/gamma)*Aold + 1/gamma*Atmp;
  V += A*delta*dt;
  U += A*alpha*dt*dt;
  
  return gsStatus::Success;
}


template <class T, bool _NL>
template <bool _nonlinear>
typename std::enable_if<(_nonlinear==true), gsStatus>::type
gsDynamicWilson<T,_NL>::_step_impl(const T t, const T dt, gsVector<T> & U, gsVector<T> & V, gsVector<T> & A) const
{
  T alpha = m_options.getReal("alpha");
  T delta = m_options.getReal("delta");
  T gamma = m_options.getReal("gamma");

  gsVector<T> Uold = U;
  gsVector<T> Vold = V;
  gsVector<T> Aold = A;

  gsVector<T> Utmp = U;
  gsVector<T> Vtmp = V;
  gsVector<T> Atmp = A;

  /// Perform a Newmark step with DT = gamma*dt
  Base::_step(t,gamma*dt,Utmp,Vtmp,Atmp);

  // Compute solutions (see Wilson et al. (1973) NONLINEAR DYNAMIC ANALYSIS OF COMPLEX STRUCTURES)
  A = (1-1/gamma)*Aold + 1/gamma*Atmp;
  U = Uold + dt*Vold + Aold*(0.5 - alpha)*dt*dt + alpha*dt*dt*A;
  V = Vold + Aold*(1 - delta)*dt + delta*dt*A;

  return gsStatus::Success;
}

template <class T, bool _NL>
void gsDynamicWilson<T,_NL>::_stageOutput(index_t stage) const
{
  if (m_options.getSwitch("Verbose"))
  {
    gsInfo<<"Stage "<<stage<<":";
    gsInfo<<"\n";
  }
}

template <class T, bool _NL>
void gsDynamicWilson<T,_NL>::_initOutput() const
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
void gsDynamicWilson<T,_NL>::_stepOutput(const index_t it, const T resnorm, const T updatenorm) const
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