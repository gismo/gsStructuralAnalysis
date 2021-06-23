 /** @file gsStaticSolver.hpp

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
gsVector<T> gsStaticSolver<T>::solveLinear()
{
    this->getOptions();
    if (m_verbose==2)
    {
        gsInfo<<"Matrix: \n"<<m_linear.toDense()<<"\n";
        gsInfo<<"Vector: \n"<<m_force<<"\n";
    }
    m_solver.compute( m_linear );
    m_solVec = m_solver.solve(m_force);
    m_converged = true;
    return m_solVec;
};

template <class T>
gsVector<T> gsStaticSolver<T>::solveNonlinear()
{
    m_converged = false;
    this->getOptions();
    if (m_solVec.rows()==0)
        m_solVec = this->solveLinear();

    gsVector<T> resVec = m_residual(m_solVec);
    T residual = resVec.norm();
    if (residual==0) residual=1;
    T residual0 = residual;
    T residualOld = residual;
    m_DeltaU = gsVector<T>::Zero(m_solVec.rows());
    m_deltaU = gsVector<T>::Zero(m_solVec.rows());
    gsSparseMatrix<T> jacMat;


    if (m_verbose>0)
    {
        gsInfo<<"\t";
        gsInfo<<std::setw(4)<<std::left<<"It.";
        gsInfo<<std::setw(17)<<std::left<<"|R|";
        gsInfo<<std::setw(17)<<std::left<<"|R|/|R0|";
        gsInfo<<std::setw(17)<<std::left<<"|dU|";
        gsInfo<<std::setw(17)<<std::left<<"|dU|/|DU|";
        gsInfo<<std::setw(17)<<std::left<<"|dU|/|U|";
        gsInfo<<std::setw(17)<<std::left<<"log(Ri/R0):";
        gsInfo<<std::setw(17)<<std::left<<"log(Ri+1/R0)";
        gsInfo<<"\n";
    }

    for (m_iterations = 0; m_iterations != m_maxIterations; ++m_iterations)
    {
        jacMat = m_nonlinear(m_solVec+m_DeltaU);
        if (m_verbose==2)
        {
            gsInfo<<"Matrix: \n"<<jacMat.toDense()<<"\n";
            gsInfo<<"Vector: \n"<<resVec<<"\n";
        }
        m_solver.compute(jacMat);
        m_deltaU = m_solver.solve(resVec); // this is the UPDATE
        m_DeltaU += m_relax * m_deltaU;

        resVec = m_residual(m_solVec+m_DeltaU);
        residual = resVec.norm();

        if (m_verbose>0)
        {
            gsInfo<<"\t";
            gsInfo<<std::setw(4)<<std::left<<m_iterations;
            gsInfo<<std::setw(17)<<std::left<<residual;
            gsInfo<<std::setw(17)<<std::left<<residual/residual0;
            gsInfo<<std::setw(17)<<std::left<<m_relax * m_deltaU.norm();
            gsInfo<<std::setw(17)<<std::left<<m_relax * m_deltaU.norm()/m_DeltaU.norm();
            gsInfo<<std::setw(17)<<std::left<<m_relax * m_deltaU.norm()/m_solVec.norm();
            gsInfo<<std::setw(17)<<std::left<<math::log10(residualOld/residual0);
            gsInfo<<std::setw(17)<<std::left<<math::log10(residual/residual0);
            gsInfo<<"\n";
        }

        residualOld = residual;

        if (m_relax * m_deltaU.norm()/m_solVec.norm()  < m_toleranceU && residual/residual0 < m_toleranceF)
        {
            m_converged = true;
            m_solVec+=m_DeltaU;
            break;
        }
        else if (m_iterations+1 == m_maxIterations)
        {
            m_converged = false;
            gsWarn<<"Maximum iterations reached!\n";
        }


            // ADD DIRICHLET HOMOGENIZE
    }
    return m_solVec;
};

template <class T>
void gsStaticSolver<T>::_computeStability(const gsVector<T> x) const
{
    gsVector<T> stabilityVec;
    gsSparseMatrix<T> jacMat = m_nonlinear(x);

    #ifdef GISMO_WITH_SPECTRA
    index_t number = std::min(static_cast<index_t>(std::floor(jacMat.cols()/3.)),10);
    gsSpectraSymSolver<gsSparseMatrix<T>> es(jacMat,number,3*number);
    es.init();
    es.compute(Spectra::SortRule::LargestMagn,1000,1e-6);
    GISMO_ASSERT(es.info()==Spectra::CompInfo::Successful,"Spectra did not converge!"); // Reason for not converging can be due to the value of ncv (last input in the class member), which is too low.
    stabilityVec = es.eigenvalues();
    stabilityVec = stabilityVec.reverse();
    #else
    Eigen::SelfAdjointEigenSolver<gsMatrix<T>> es(jacMat);
    stabilityVec = es.eigenvalues();
    #endif

    m_indicator = stabilityVec.colwise().minCoeff()[0]; // This is required since D does not necessarily have one column.
}

template <class T>
void gsStaticSolver<T>::defaultOptions()
{
    m_options.addInt("Verbose","Verbose output",verbose::iterations);
    m_options.addInt("MaxIterations","Maximum number of iterations",25);
    m_options.addReal("Tolerance","Relative Tolerance Force",1e-6);
    m_options.addReal("ToleranceF","Relative Tolerance Force",-1);
    m_options.addReal("ToleranceU","Relative Tolerance Displacements",-1);
    m_options.addReal("Relaxation","Relaxation parameter",1);
};

template <class T>
void gsStaticSolver<T>::getOptions() const
{
    m_maxIterations = m_options.getInt("MaxIterations");
    m_toleranceF = m_options.getReal("ToleranceF")!=-1 ? m_options.getReal("ToleranceF") : m_options.getReal("Tolerance");
    m_toleranceU = m_options.getReal("ToleranceU")!=-1 ? m_options.getReal("ToleranceU") : m_options.getReal("Tolerance");
    m_verbose = m_options.getInt("Verbose");
    m_relax = m_options.getReal("Relaxation");
};

} // namespace gismo
