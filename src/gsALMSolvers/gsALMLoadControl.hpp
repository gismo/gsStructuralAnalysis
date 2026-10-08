/** @file gsALMLoadControl.hpp

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
void gsALMLoadControl<T>::_install()
{
    m_load.reset(new gsLoadConstraint<T>());
    this->setConstraint(m_load);
    defaultOptions();
    getOptions();
}

// The load constraint takes the load increment from the step length of the
// controller, so the guess step replaces the length for one step.
template <class T>
void gsALMLoadControl<T>::applyInitialGuess(const gsVector<T> & /*Uguess*/, T Lguess)
{
    if (!m_guessStep)
    {
        m_savedController = this->driver().controller();
        m_guessStep = true;
    }
    this->driver().controller().setLength(Lguess - this->m_L);
}

// The controller adapts its length to the iterations of the step: the guess
// step is therefore accounted to the length it replaced.
template <class T>
gsStatus gsALMLoadControl<T>::step()
{
    const gsStatus st = Base::step();
    if (m_guessStep)
    {
        m_guessStep = false;
        this->driver().controller() = m_savedController;
        if (st == gsStatus::Success)
        {
            this->driver().controller().accepted(this->m_numIterations);
            this->m_Uguess.resize(0);
        }
        this->m_arcLength      = this->driver().controller().length();
        this->m_arcLength_prev = this->driver().controller().previousLength();
    }
    return st;
}

} // namespace gismo
