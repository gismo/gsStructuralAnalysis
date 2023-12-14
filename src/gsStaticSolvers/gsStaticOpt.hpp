 /** @file gsStaticNewton.hpp

    @brief Static solver using optimization methods

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#include <gsOptimizer/gsOptProblem.h>
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

template <class T>
void gsStaticOpt<T>::defaultOptions()
{
    Base::defaultOptions();
}

template <class T>
void gsStaticOpt<T>::getOptions()
{
    Base::getOptions();
}

template <class T>
void gsStaticOpt<T>::initOutput()
{
    gsInfo<<"\t";
    gsInfo<<std::setw(4)<<std::left<<"It.";
    gsInfo<<std::setw(17)<<std::left<<"|R|";
    gsInfo<<std::setw(17)<<std::left<<"|R|/|R0|";
    gsInfo<<std::setw(17)<<std::left<<" Ek";
    gsInfo<<std::setw(17)<<std::left<<" Ek/Ek0";
    gsInfo<<std::setw(17)<<std::left<<"|dU|";
    gsInfo<<std::setw(17)<<std::left<<"|dU|/|DU|";
    gsInfo<<std::setw(17)<<std::left<<"|dU|/|U+DU|";
    gsInfo<<std::setw(17)<<std::left<<"|dV|";
    gsInfo<<"\n";
}

template <class T>
void gsStaticOpt<T>::stepOutput(index_t k)
{
    gsInfo<<"\t";
    gsInfo<<std::setw(4)<<std::left<<k;
    gsInfo<<std::setw(17)<<std::left<<m_residual;
    gsInfo<<std::setw(17)<<std::left<<m_residual/m_residualIni;
    gsInfo<<std::setw(17)<<std::left<<m_Ek;
    gsInfo<<std::setw(17)<<std::left<<m_Ek/m_Ek0;
    gsInfo<<std::setw(17)<<std::left<<m_deltaU.norm();
    gsInfo<<std::setw(17)<<std::left<<m_deltaU.norm()/m_DeltaU.norm();
    gsInfo<<std::setw(17)<<std::left<<m_deltaU.norm()/(m_U+m_DeltaU).norm();
    gsInfo<<std::setw(17)<<std::left<<m_v.norm();
    gsInfo<<"\n";
}

template <class T>
gsStatus gsStaticOpt<T>::solve()
{
    try
    {
        _solve();
        m_status = gsStatus::Success;
    }
    catch (int errorCode)
    {
        if      (errorCode==1)
            m_status = gsStatus::NotConverged;
        else if (errorCode==2)
            m_status = gsStatus::AssemblyError;
        else if (errorCode==3)
            m_status = gsStatus::SolverError;
        else
            m_status = gsStatus::OtherError;
    }
    catch (...)
    {
        m_status = gsStatus::OtherError;
    }
    return m_status;
}

template <class T>
void gsStaticOpt<T>::_solve()
{
    m_Eks.clear();
    m_Eks.reserve(m_maxIterations);

    if (m_verbose) initOutput();

    _start();

    m_Ek0 = m_Ek;
    m_Eks.push_back(m_Ek);

    if (m_verbose != 0) stepOutput(0);
    index_t resetIt = 0;
    for (m_numIterations=1; m_numIterations!=m_maxIterations; m_numIterations++, resetIt++)
    {
        _iteration();
        if ((m_c==0 && m_Ek_prev > m_Ek) || resetIt==m_resetIterations)// || (m_Ek/m_Ek_prev > 1/m_tolE && m_Ek_prev!=0))
        {
            resetIt = 0;
            _peak();
        }

        if (m_verbose!=0)
            if (m_numIterations % m_verbose == 0 || m_verbose==-1 ) stepOutput(m_numIterations);

        m_Eks.push_back(m_Ek);

        m_residualOld = m_residual;

        if (m_residual/m_residualIni < m_tolF && m_Ek/m_Ek0 < m_tolE && m_deltaU.norm()/m_DeltaU.norm() < m_tolU)
        {
            m_U += m_DeltaU;
            gsDebug <<"Converged: \n";
            gsDebug <<"\t |R|/|R0| = "<<m_residual/m_residualIni<<" < tolF = "<<m_tolF<<"\n";
            gsDebug <<"\t |E|/|E0| = "<<m_Ek/m_Ek0              <<" < tolE = "<<m_tolE<<"\n";
            gsDebug <<"\t |U|/|U0| = "<<m_deltaU.norm()/m_DeltaU.norm()<<" < tolF = "<<m_tolU<<"\n";
            break;
        }
        if (m_numIterations==m_maxIterations-1)
        {
            m_U += m_DeltaU;
            gsInfo<<"Maximum iterations reached. Solution did not converge\n";
            throw 1;
        }
    }
    gsInfo<<"\n";
};

template <class T>
void gsStaticOpt<T>::initialize()
{
    this->reset();
    getOptions();
}

template <class T>
void gsStaticOpt<T>::reset()
{
    Base::reset();
    m_status = gsStatus::NotStarted;
}

template <class T>
void gsStaticOpt<T>::_init()
{
    this->reset();
    defaultOptions();
}

template <class T>
gsVector<T> gsStaticOpt<T>::_computeResidual(const gsVector<T> & U)
{
  gsVector<T> resVec;
  if (!m_residualFun(U, resVec))
    throw 2;
  return resVec;
}

template <class T>
gsOptProblemStatic<T>::gsOptProblemStatic()
{

    // Number of design variables
    m_numDesignVars  = 2;
    // Number of constraints
    m_numConstraints = 0;
    // Number of non-zeros in the Jacobian of the constraints
    m_numConJacNonZero = 0;

    m_desLowerBounds.resize(m_numDesignVars);
    m_desUpperBounds.resize(m_numDesignVars);

    m_desLowerBounds.setConstant(-1e19);
    m_desUpperBounds.setConstant( 1e19);

    // we initialize x in bounds, in the upper right quadrant
    m_curDesign.resize(m_numDesignVars,1);
    m_curDesign.setZero();
}

} // namespace gismo
