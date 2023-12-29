/** @file gsDynamicBase.hpp

    @brief Base class to perform time integration of second-order structural dynamics systems

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#pragma once

#include <typeinfo>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMHelper.h>

namespace gismo
{

template <class T>
void gsDynamicBase<T>::defaultOptions()
{
    m_options.addInt ("MaxIter","Maximum iterations",100);
    m_options.addReal("Tol","Tolerance",1e-6);
    m_options.addReal("TolF","Tolerance",1e-3);
    m_options.addReal("TolU","Tolerance",1e-6);

    m_options.addReal("DT","Time step",1e-2);

    m_options.addSwitch("Quasi","Use Quasi Newton method",false);
    m_options.addInt ("QuasiIterations","Number of iterations for quasi newton method",-1);

    m_options.addString("Solver","Sparse linear solver", "SimplicialLDLT");

    m_options.addSwitch ("Verbose","Verbose output",false);
}

template <class T>
void gsDynamicBase<T>::getOptions()
{
    m_maxIterations       = m_options.getInt ("MaxIter");
    m_tolerance           = m_options.getReal("Tol");
    m_toleranceF          = m_options.getReal("TolF");
    m_toleranceU          = m_options.getReal("TolU");

    m_dt = m_options.getReal("DT");

    m_quasiNewton         = m_options.getSwitch("Quasi");
    m_quasiNewtonInterval = m_options.getInt ("QuasiIterations");

    m_solver = gsSparseSolver<T>::get( m_options.getString("Solver") );

    m_verbose             = m_options.getSwitch ("Verbose");
}

} // namespace gismo
