/** @file gsALMCrisfield.hpp

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
void gsALMCrisfield<T>::_install()
{
    m_cris.reset(new gsCrisfieldConstraint<T>());
    this->setConstraint(m_cris);
    m_cris->setForce(this->m_forcing);
    defaultOptions();
    getOptions();
}

template <class T>
void gsALMCrisfield<T>::defaultOptions()
{
    Base::defaultOptions();
    this->m_options.addReal("Scaling","Set Scaling factor Phi",-1);
    this->m_options.addInt ("AngleMethod","Angle determination method: 0 = Previous step; 1 = Previous iteration",angmethod::Step);
}

template <class T>
void gsALMCrisfield<T>::getOptions()
{
    Base::getOptions();
    m_cris->setScaling(this->m_options.getReal("Scaling"));
    m_cris->setAngleMethod(this->m_options.getInt("AngleMethod"));
    m_cris->setRelaxation(this->m_relax);
}

template <class T>
void gsALMCrisfield<T>::initMethods()
{
    m_cris->resetCounters();
}

// The previous increment of the predictor is -(Uguess - U, Lguess - L), the
// increment from the guess back to the present point: the predictor takes the
// sign of its inner product with (u_t, 1), so a guess ahead on the path turns
// the step around.
template <class T>
void gsALMCrisfield<T>::applyInitialGuess(const gsVector<T> & Uguess, T Lguess)
{
    if (Uguess.size() != this->m_U.size())
        return;

    gsStepHistory<T> h;
    h.DeltaU     = -(Uguess - this->m_U);
    h.DeltaL     = -(Lguess - this->m_L);
    h.stepSize   = this->driver().controller().previousLength();
    h.iterations = 0;
    h.valid      = true;
    this->driver().setHistory(h);
    m_cris->requestTangentScaling();
    this->m_Uguess.resize(0);
}

} // namespace gismo
