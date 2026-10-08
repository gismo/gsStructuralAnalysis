/** @file gsALMRiks.hpp

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
void gsALMRiks<T>::_install()
{
    m_riks.reset(new gsRiksConstraint<T>());
    this->setConstraint(m_riks);
    defaultOptions();
    getOptions();
}

template <class T>
void gsALMRiks<T>::applyInitialGuess(const gsVector<T> & Uguess, T Lguess)
{
    Base::applyInitialGuess(Uguess, Lguess);
    this->m_Uguess.resize(0);
}

} // namespace gismo
