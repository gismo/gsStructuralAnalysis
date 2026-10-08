/** @file gsALMConsistentCrisfield.hpp

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
void gsALMConsistentCrisfield<T>::_install()
{
    m_cc.reset(new gsConsistentCrisfieldConstraint<T>());
    this->setConstraint(m_cc);
    m_cc->setForce(this->m_forcing);
    defaultOptions();
    getOptions();
}

template <class T>
void gsALMConsistentCrisfield<T>::defaultOptions()
{
    Base::defaultOptions();
    this->m_options.addReal("Scaling","Set Scaling factor Phi",-1);
}

template <class T>
void gsALMConsistentCrisfield<T>::getOptions()
{
    Base::getOptions();
    m_cc->setScaling(this->m_options.getReal("Scaling"));
}

template <class T>
void gsALMConsistentCrisfield<T>::initMethods()
{
    m_cc->resetFrozenScaling();
}

} // namespace gismo
