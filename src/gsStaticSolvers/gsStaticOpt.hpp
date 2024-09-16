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



//! [OptProblem]

// template <class T, class Optimizer>
// void gsStaticOpt<T,Optimizer>::getOptions()
// {
//     Base::getOptions();
// }

// template <class T, class Optimizer>
// gsStatus gsStaticOpt<T,Optimizer>::solve()
// {
//     try
//     {
//         _solve();
//         m_status = gsStatus::Success;
//     }
//     catch (int errorCode)
//     {
//         if      (errorCode==1)
//             m_status = gsStatus::NotConverged;
//         else if (errorCode==2)
//             m_status = gsStatus::AssemblyError;
//         else if (errorCode==3)
//             m_status = gsStatus::SolverError;
//         else
//             m_status = gsStatus::OtherError;
//     }
//     catch (...)
//     {
//         m_status = gsStatus::OtherError;
//     }
//     return m_status;
// }

// template <class T, class Optimizer>
// void gsStaticOpt<T,Optimizer>::_solve()
// {
//     m_Eks.clear();
//     m_Eks.reserve(m_maxIterations);

//     if (m_verbose) initOutput();

//     _start();

//     m_Ek0 = m_Ek;
//     m_Eks.push_back(m_Ek);

//     if (m_verbose != 0) stepOutput(0);
//     index_t resetIt = 0;
//     for (m_numIterations=1; m_numIterations!=m_maxIterations; m_numIterations++, resetIt++)
//     {
//         _iteration();
//         if ((m_c==0 && m_Ek_prev > m_Ek) || resetIt==m_resetIterations)// || (m_Ek/m_Ek_prev > 1/m_tolE && m_Ek_prev!=0))
//         {
//             resetIt = 0;
//             _peak();
//         }

//         if (m_verbose!=0)
//             if (m_numIterations % m_verbose == 0 || m_verbose==-1 ) stepOutput(m_numIterations);

//         m_Eks.push_back(m_Ek);

//         m_residualOld = m_residual;

//         if (m_residual/m_residualIni < m_tolF && m_Ek/m_Ek0 < m_tolE && m_deltaU.norm()/m_DeltaU.norm() < m_tolU)
//         {
//             m_U += m_DeltaU;
//             gsDebug <<"Converged: \n";
//             gsDebug <<"\t |R|/|R0| = "<<m_residual/m_residualIni<<" < tolF = "<<m_tolF<<"\n";
//             gsDebug <<"\t |E|/|E0| = "<<m_Ek/m_Ek0              <<" < tolE = "<<m_tolE<<"\n";
//             gsDebug <<"\t |U|/|U0| = "<<m_deltaU.norm()/m_DeltaU.norm()<<" < tolF = "<<m_tolU<<"\n";
//             break;
//         }
//         if (m_numIterations==m_maxIterations-1)
//         {
//             m_U += m_DeltaU;
//             gsInfo<<"Maximum iterations reached. Solution did not converge\n";
//             throw 1;
//         }
//     }
//     gsInfo<<"\n";
// };

// template <class T, class Optimizer>
// void gsStaticOpt<T,Optimizer>::initialize()
// {
//     this->reset();
//     getOptions();
// }

// template <class T, class Optimizer>
// void gsStaticOpt<T,Optimizer>::reset()
// {
//     Base::reset();
//     m_status = gsStatus::NotStarted;
// }

// template <class T, class Optimizer>
// void gsStaticOpt<T,Optimizer>::_init()
// {
//     this->reset();
//     defaultOptions();
// }

// template <class T, class Optimizer>
// gsVector<T> gsStaticOpt<T,Optimizer>::_computeResidual(const gsVector<T> & U)
// {
//   gsVector<T> resVec;
//   if (!m_residualFun(U, resVec))
//     throw 2;
//   return resVec;
// }

// template <class T, class Optimizer>
// gsOptProblemStatic<T>::gsOptProblemStatic()
// {

//     // Number of design variables
//     m_numDesignVars  = 2;
//     // Number of constraints
//     m_numConstraints = 0;
//     // Number of non-zeros in the Jacobian of the constraints
//     m_numConJacNonZero = 0;

//     m_desLowerBounds.resize(m_numDesignVars);
//     m_desUpperBounds.resize(m_numDesignVars);

//     m_desLowerBounds.setConstant(-1e19);
//     m_desUpperBounds.setConstant( 1e19);

//     // we initialize x in bounds, in the upper right quadrant
//     m_curDesign.resize(m_numDesignVars,1);
//     m_curDesign.setZero();
// }

} // namespace gismo
