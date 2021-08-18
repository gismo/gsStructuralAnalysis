/** @file gsALMBase.hpp

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
void gsALMExplorer<T>::defaultOptions()
{
    m_options.addInt ("maxBranches","Maximum number of branches to follow",100);
    m_options.addInt ("maxSteps","Maximum number of steps per branch",100);

};

template <class T>
void gsALMExplorer<T>::getOptions() const
{
    // m_maxIterations = m_options.getInt("MaxIterations");
    // m_toleranceF = m_options.getReal("ToleranceF")!=-1 ? m_options.getReal("ToleranceF") : m_options.getReal("Tolerance");
    // m_toleranceU = m_options.getReal("ToleranceU")!=-1 ? m_options.getReal("ToleranceU") : m_options.getReal("Tolerance");
    // m_verbose = m_options.getInt("Verbose");
    // m_relax = m_options.getReal("Relaxation");

    // m_bifurcationMethod   = m_options.getInt ("BifurcationMethod");
    // m_solverType          = m_options.getInt ("Solver");
    // if (m_solverType!=solver::LDLT && m_bifurcationMethod==bifmethod::Determinant)
    // {
    //   gsWarn<<"Determinant method cannot be used with solvers other than LDLT. Bifurcation method will be set to 'Eigenvalue'.\n";
    //   m_bifurcationMethod = bifmethod::Eigenvalue;
    // }
};

} // namespace gismo