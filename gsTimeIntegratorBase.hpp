/** @file gsTimeIntegrator.hpp

    @brief Provides temporal solvers for structural analysis problems

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
    void gsTimeIntegratorBase<T>::defaultOptions()
    {
        m_options.addInt ("MaxIter","Maximum iterations",100);
        m_options.addReal("Tol","Tolerance",1e-6);
        m_options.addSwitch ("Verbose","Verbose output",false);
        m_options.addSwitch("Quasi","Use Quasi Newton method",false);
    }


    template <class T>
    void gsTimeIntegratorBase<T>::getOptions()
    {
        m_maxIterations       = m_options.askInt("MaxIter",100);
        m_tolerance           = m_options.askReal("Tol",1e-6);
        m_verbose             = m_options.askSwitch ("Verbose",false);
        m_quasiNewton         = m_options.askSwitch("Quasi",false);
    }

} // namespace gismo