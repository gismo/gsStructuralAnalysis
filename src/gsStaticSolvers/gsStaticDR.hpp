 /** @file gsStaticNewton.hpp

    @brief Performs linear modal analysis given a matrix or functions of a matrix

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#include <typeinfo>
#pragma once

namespace gismo
{

template <class T>
void gsStaticDR<T>::defaultOptions()
{
    Base::defaultOptions();
    m_options.addReal("damping","damping factor",1.0);
    m_options.addReal("alpha","mass coefficient",2.0);
    m_options.addReal("tolE","Kinetic energy tolerance",1e-6);
    m_options.addInt("ResetIt","Reset rate of velocities if damping is zero",-1);
}

template <class T>
void gsStaticDR<T>::getOptions()
{
    Base::getOptions();
    m_c = m_options.getReal("damping");
    m_alpha = m_options.getReal("alpha");
    m_tolE = m_options.getReal("tolE");
    m_resetIterations = m_options.getInt("ResetIt");
}

template <class T>
void gsStaticDR<T>::initOutput()
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
void gsStaticDR<T>::stepOutput(index_t k)
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
gsStatus gsStaticDR<T>::solve()
{
    // try
    // {
        _solve();
        m_status = gsStatus::Success;
    // }
    // catch (int errorCode)
    // {
    //     if      (errorCode==1)
    //         m_status = gsStatus::NotConverged;
    //     else if (errorCode==2)
    //         m_status = gsStatus::AssemblyError;
    //     else if (errorCode==3)
    //         m_status = gsStatus::SolverError;
    //     else
    //         m_status = gsStatus::OtherError;
    // }
    // catch (...)
    // {
    //     m_status = gsStatus::OtherError;
    // }
    return m_status;
}

template <class T>
void gsStaticDR<T>::_solve()
{
    m_Eks.clear();
    m_Eks.reserve(m_maxIterations);

    if (m_verbose) initOutput();

    _start(); // initializes m_Ek, etc

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
void gsStaticDR<T>::initialize()
{
    this->reset();
    getOptions();
    m_massInv *= 1./m_alpha;
    m_damp = m_c * m_mass;
}

template <class T>
void gsStaticDR<T>::reset()
{
    m_dofs = m_mass.rows();
    m_massInv = m_mass.array().inverse();
    m_Ek = m_Ek0 = 0.0;
    // resets m_U, m_DeltaU, m_deltaU, m_R, m_L, m_DeltaL, m_deltaL and m_headstart
    Base::reset();
    m_v.setZero(m_dofs);
    m_status = gsStatus::NotStarted;
}

template <class T>
void gsStaticDR<T>::_init()
{
    this->reset();
    m_dt = 1.0;
    defaultOptions();
}

template <class T>
gsVector<T> gsStaticDR<T>::_computeResidual(const gsVector<T> & U)
{
  gsVector<T> resVec;
  if (!m_residualFun(U, resVec))
    throw 2;
  return resVec;
}

template <class T>
void gsStaticDR<T>::_iteration()
{
    m_Ek_prev = m_Ek;
    m_R = _computeResidual(m_U+m_DeltaU) - m_damp.cwiseProduct(m_v);
    m_residual = m_R.norm();
//----------------------------------------------------------------------------------
    m_v += m_dt * m_massInv.cwiseProduct(m_R);                    // Velocities at t+dt/2
    m_deltaU = m_dt * m_v;               // Velocities at t+dt
//----------------------------------------------------------------------------------
    // m_deltaU = m_deltaU + m_dt*m_dt*m_massInv*m_R;
    // m_v = m_deltaU/m_dt;
//----------------------------------------------------------------------------------

    m_DeltaU += m_deltaU;               // Velocities at t+dt
    m_Ek = m_v.transpose() * m_mass.cwiseProduct(m_v);
}

template <class T>
void gsStaticDR<T>::_peak()
{
    m_R = _computeResidual(m_U+m_DeltaU)- m_damp.cwiseProduct(m_v);
    m_deltaU = - 1.5 * m_dt * m_v + m_dt*m_dt / 2. * m_massInv.cwiseProduct(m_R);
    m_DeltaU += m_deltaU;

    m_v.setZero();
    m_Ek = 0.0;
    // m_v = 0.5 * m_dt * m_massInv * m_R; // Velocities at dt/2
}

template <class T>
void gsStaticDR<T>::_start()
{
    m_v.setZero();

//----------------------------------------------------------------------------------
    // Define residual measures:
    // m_residual    = current residual
    // m_residualIni = residual on m_U
    // m_residualOld = residual in previous step
    if (!m_headstart) // no headstart
    {
        // We can reset the update to ensure we properly restart
        m_DeltaU.setZero(m_dofs);
        // Compute current residual and its norm
        m_R = _computeResidual(m_U);
        // m_residual = m_R.norm();
        m_residual = m_forcing.norm();
        // If the residual is 0 (e.g. with purely displacment loading), we set it to 1 to allow divisions
        if (m_residual==0) m_residual=1;
        // All residual norms are equal
        m_residualIni = m_residualOld = m_residual;
    }
    else
    {
        // Compute current residual and its norm
        m_R = _computeResidual(m_U + m_DeltaU);
        // m_residual = m_R.norm();
        m_residual = m_forcing.norm();
        // If the residual is 0 (e.g. with purely displacment loading), we set it to 1 to allow divisions
        if (m_residual==0) m_residual=1;
        // The previous step residual is the same as the residual
        m_residualOld = m_residual;
        // Residual0 is the residual without m_DeltaU
        m_residualIni = _computeResidual(m_U).norm();
        // If the residual is 0 (e.g. with purely displacment loading), we set it to 1 to allow divisions
        if (m_residualIni==0) m_residualIni=1;

        // Now we can reset the headstart
        m_headstart = false;
    }

//----------------------------------------------------------------------------------
    // m_R = _computeResidual(m_U+m_DeltaU);
    m_deltaU = 0.5*m_dt*m_dt*m_massInv.cwiseProduct(m_R);
    m_v = m_deltaU/m_dt;
//----------------------------------------------------------------------------------

    m_DeltaU += m_deltaU;

    m_Ek = m_v.transpose() * m_mass.cwiseProduct(m_v);
}

} // namespace gismo
