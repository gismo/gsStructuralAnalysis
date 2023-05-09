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
void gsStaticNewton<T>::defaultOptions()
{
    Base::defaultOptions();
    m_options.setString("Solver","CGDiagonal"); // The CG solver is robust for membrane models, where zero-blocks in the matrix might occur.
    m_options.addReal("Relaxation","Relaxation parameter",1);
};

template <class T>
void gsStaticNewton<T>::getOptions()
{
    Base::getOptions();

    // if (m_solverType!=solver::LDLT && m_stabilityMethod==stabmethod::Determinant)
    // {
    //   gsWarn<<"Determinant method cannot be used with solvers other than LDLT. Bifurcation method will be set to 'Eigenvalue'.\n";
    //   m_stabilityMethod = stabmethod::Eigenvalue;
    // }

    m_relax = m_options.getReal("Relaxation");
};

template <class T>
void gsStaticNewton<T>::initOutput()
{
    gsInfo<<"\t";
    gsInfo<<std::setw(4)<<std::left<<"It.";
    gsInfo<<std::setw(17)<<std::left<<"|R|";
    gsInfo<<std::setw(17)<<std::left<<"|R|/|R0|";
    gsInfo<<std::setw(17)<<std::left<<"|dU|";
    gsInfo<<std::setw(17)<<std::left<<"|dU|/|DU|";
    gsInfo<<std::setw(17)<<std::left<<"|dU|/|U+DU|";
    gsInfo<<std::setw(17)<<std::left<<"log(Ri/R0):";
    gsInfo<<std::setw(17)<<std::left<<"log(Ri+1/R0)";
    gsInfo<<"\n";
}

template <class T>
void gsStaticNewton<T>::stepOutput(index_t k)
{
    gsInfo<<"\t";
    gsInfo<<std::setw(4)<<std::left<<k;
    gsInfo<<std::setw(17)<<std::left<<m_residual;
    gsInfo<<std::setw(17)<<std::left<<m_residual/m_residualIni;
    gsInfo<<std::setw(17)<<std::left<<m_relax * m_deltaU.norm();
    gsInfo<<std::setw(17)<<std::left<<m_relax * m_deltaU.norm()/m_DeltaU.norm();
    gsInfo<<std::setw(17)<<std::left<<m_relax * m_deltaU.norm()/(m_U+m_DeltaU).norm();
    gsInfo<<std::setw(17)<<std::left<<math::log10(m_residualOld/m_residualIni);
    gsInfo<<std::setw(17)<<std::left<<math::log10(m_residual/m_residualIni);
    gsInfo<<"\n";

}

template <class T>
gsVector<T> gsStaticNewton<T>::solveLinear()
{
    this->getOptions();
    if (m_verbose==2)
    {
        gsInfo<<"Matrix: \n"<<m_linear.toDense()<<"\n";
        gsInfo<<"Vector: \n"<<m_force<<"\n";
    }
    _factorizeMatrix( m_linear );
    m_U = _solveSystem(m_force);
    m_converged = true;


    return m_U;
};

template <class T>
gsVector<T> gsStaticNewton<T>::solveNonlinear()
{
    this->getOptions();
    // m_start: true -> m_U given
    // m_headstart: true -> m_DeltaU given

    if (m_DeltaU.norm()==0 && m_DeltaU.rows()==0) ///
    {
        m_deltaU = m_DeltaU = this->solveLinear();
        m_U.setZero(); // Needed because linear solve modifies m_U.
        m_headstart = true; // due to this, the relative residual is based on the solution of the linear solve
    }
    _start();

    gsSparseMatrix<T> jacMat;

    if (m_verbose>0) { initOutput(); }

    for (m_numIterations = 0; m_numIterations != m_maxIterations; ++m_numIterations)
    {
        jacMat = m_dnonlinear(m_U+m_DeltaU,m_deltaU);
        if (m_verbose==2)
        {
            gsInfo<<"Matrix: \n"<<jacMat.toDense()<<"\n";
            gsInfo<<"Vector: \n"<<m_R<<"\n";
        }

        _factorizeMatrix(jacMat);
        m_deltaU = _solveSystem(m_R); // this is the UPDATE
        m_DeltaU += m_relax * m_deltaU;

        m_R = m_residualFun(m_U+m_DeltaU);
        m_residual = m_R.norm();

        if (m_verbose>0) { stepOutput(m_numIterations); }

        m_residualOld = m_residual;

        if (m_relax * m_deltaU.norm()/m_DeltaU.norm()  < m_tolU && m_residual/m_residualIni < m_tolF)
        {
            m_converged = true;
            m_U+=m_DeltaU;
            break;
        }
        else if (m_numIterations+1 == m_maxIterations)
        {
            m_converged = false;
            gsWarn<<"Maximum iterations reached!\n";
        }
    }

    return m_U;
};

template <class T>
void gsStaticNewton<T>::_factorizeMatrix(const gsSparseMatrix<T> & jacMat) const
{
    m_solver->compute(jacMat);
    // If 1: matrix is not SPD
    GISMO_ENSURE(m_solver->info()==Eigen::ComputationInfo::Success,"Solver error with code "<<m_solver->info()<<". See Eigen documentation on ComputationInfo \n"
                 <<Eigen::ComputationInfo::Success<<": Success"<<"\n"
                 <<Eigen::ComputationInfo::NumericalIssue<<": NumericalIssue"<<"\n"
                 <<Eigen::ComputationInfo::NoConvergence<<": NoConvergence"<<"\n"
                 <<Eigen::ComputationInfo::InvalidInput<<": InvalidInput"<<"\n");
}

template <class T>
gsVector<T> gsStaticNewton<T>::_solveSystem(const gsVector<T> & F)
{
    return m_solver->solve(F);
}

template <class T>
void gsStaticNewton<T>::_init()
{
    m_stabilityMethod = 0;
    m_converged = false;
    m_start = false;
    m_headstart = false;
    m_dofs = m_force.rows();

    if (m_dofs==0)
        gsWarn<<"The number of degrees of freedom is equal to zero. This can lead to bad initialization.\n";

    m_U.setZero(m_dofs);
    m_DeltaU.setZero(m_dofs);
    m_deltaU.setZero(m_dofs);
    m_R.setZero(m_dofs);

    m_residual = m_residualIni = m_residualOld = 0;

    defaultOptions();
}

template <class T>
void gsStaticNewton<T>::_start()
{
    m_converged = false;

    // Define residual measures:
    // residual    = current residual
    // residual0   = residual on m_U
    // residualOld = residual in previous step

    if (!m_headstart) // no headstart
    {
        // We can reset the update to ensure we properly restart
        if (m_dofs==0) gsWarn<<"The number of degrees of freedom is equal to zero. This can lead to bad initialization.\n";
        m_DeltaU.setZero(m_dofs);
        // Compute current residual and its norm
        m_R = m_residualFun(m_U);
        m_residual = m_R.norm();
        // If the residual is 0 (e.g. with purely displacment loading), we set it to 1 to allow divisions
        if (m_residual==0) m_residual=1;
        // All residual norms are equal
        m_residualIni = m_residualOld = m_residual;
    }
    else
    {
        // If we have a headstart, we need to compute Residual0 on the solution m_U
        // Compute current residual and its norm
        m_R = m_residualFun(m_U + m_DeltaU);
        m_residual = m_R.norm();
        // If the residual is 0 (e.g. with purely displacment loading), we set it to 1 to allow divisions
        if (m_residual==0) m_residual=1;
        // The previous step residual is the same as the residual
        m_residualOld = m_residual;
        // Residual0 is the residual without m_DeltaU
        m_residualIni = m_residualFun(m_U).norm();
        // If the residual is 0 (e.g. with purely displacment loading), we set it to 1 to allow divisions
        if (m_residualIni==0) m_residualIni=1;

        // Now we can reset the headstart
        m_headstart = false;
    }

    // Reset incremental update
    m_deltaU.setZero(m_dofs);

}

} // namespace gismo
