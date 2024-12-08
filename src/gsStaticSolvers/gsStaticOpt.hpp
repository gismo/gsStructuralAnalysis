 /** @file gsStaticNewton.hpp

    @brief Static solver using optimization methods

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#include <gsOptimizer/gsGradientDescent.h>

#ifdef gsHLBFGS_ENABLED
#include <gsHLBFGS/gsHLBFGS.h>
#endif

#ifdef gsOptim_ENABLED
#include <gsOptim/gsOptim.h>
#endif

#ifdef gsIpOpt_ENABLED
#include <gsIpOpt/gsIpOpt.h>
#endif

#pragma once

namespace gismo
{

template <class T, class Optimizer>
void gsStaticOpt<T,Optimizer>::defaultOptions()
{
    Base::defaultOptions();
}

template <class T, class Optimizer>
void gsStaticOpt<T,Optimizer>::_init()
{
    m_U.setZero(m_dofs);
    m_DeltaU.setZero(m_dofs);

    defaultOptions();

    m_status = gsStatus::NotStarted;
    m_numIterations = 0;

    // m_optProblem = gsOptProblemStatic<T>(m_residualFun, m_dofs);
    // m_optimizer = Optimizer(&m_optProblem);
}

template <class T, class Optimizer>
void gsStaticOpt<T,Optimizer>::getOptions()
{
    Base::getOptions();

    m_optimizer.options().setInt("MaxIterations",m_maxIterations);
    m_optimizer.options().setInt("Verbose",m_verbose);
    if (dynamic_cast<gsGradientDescent<T>*>(&m_optimizer))
    {
        m_optimizer.options().setReal("MinGradientLength",m_tolF);
        m_optimizer.options().setReal("MinStepLength",m_tolU);
    }
#ifdef gsHLBFGS_ENABLED
    else if (dynamic_cast<gsHLBFGS<T>*>(&m_optimizer))
    {
        m_optimizer.options().setReal("MinGradientLength",m_tolF);
        m_optimizer.options().setReal("MinStepLength",m_tolU);
    }
#endif
#ifdef gsOptim_ENABLED
    else if (dynamic_cast<gsOptim<T>*>(&m_optimizer))
    {
        m_optimizer.options().setReal("RelObjFnChangeTol",m_tolF);
        m_optimizer.options().setReal("RelSolChangeTol",m_tolU);
    }
#endif
#ifdef gsIpOpt_ENABLED
    else if (dynamic_cast<gsIpOpt<T>*>(&m_optimizer))
    {
        m_optimizer.options().setReal("MinGradientLength",m_tolF);
        m_optimizer.options().setReal("MinStepLength",m_tolU);
    }
#endif
    else
        GISMO_ERROR("Optimizer not recognized\n");
}

template <class T, class Optimizer>
gsStatus gsStaticOpt<T,Optimizer>::solve()
{
    this->getOptions();
    m_optimizer.solve(m_U+m_DeltaU);

    m_DeltaU = m_optimizer.currentDesign() - m_U;
    m_U += m_DeltaU;

    m_numIterations = m_optimizer.iterations();
    return gsStatus::Success;
}
} // namespace gismo
