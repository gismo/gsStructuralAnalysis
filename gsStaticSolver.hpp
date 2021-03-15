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

    T residual = m_force.norm();
    T residual0 = residual;
    T residualOld = residual;
    gsVector<T> updateVector = m_solVec;
    gsVector<T> resVec = m_residual(m_solVec);
    gsSparseMatrix<T> jacMat;
    for (m_iterations = 0; m_iterations != m_maxIterations; ++m_iterations)
    {
        jacMat = m_nonlinear(m_solVec);
        if (m_verbose==2)
        {
            gsInfo<<"Matrix: \n"<<jacMat.toDense()<<"\n";
            gsInfo<<"Vector: \n"<<resVec<<"\n";
        }
        m_solver.compute(jacMat);
        updateVector = m_solver.solve(resVec); // this is the UPDATE
        m_solVec += updateVector;

        resVec = m_residual(m_solVec);
        residual = resVec.norm();

        if (m_verbose>0)
        {
            gsInfo<<"Iteration: "<< m_iterations
               <<", residue: "<< residual
               <<", update norm: "<<updateVector.norm()
               <<", log(Ri/R0): "<< math::log10(residualOld/residual0)
               <<", log(Ri+1/R0): "<< math::log10(residual/residual0)
               <<"\n";
        }

        residualOld = residual;

        if (updateVector.norm() < m_tolerance)
        {
            m_converged = true;
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

    index_t number = std::min(static_cast<index_t>(std::floor(jacMat.cols()/3.)),10);
    gsSpectraSymSolver<gsSparseMatrix<T>> es(jacMat,number,3*number);
    es.init();
    es.compute(Spectra::SortRule::LargestMagn,1000,1e-6);
    GISMO_ASSERT(es.info()==Spectra::CompInfo::Successful,"Spectra did not converge!"); // Reason for not converging can be due to the value of ncv (last input in the class member), which is too low.
    stabilityVec = es.eigenvalues();
    stabilityVec = stabilityVec.reverse();
    m_indicator = stabilityVec.colwise().minCoeff()[0]; // This is required since D does not necessarily have one column.
}

template <class T>
void gsStaticSolver<T>::defaultOptions()
{
    m_options.addInt("Verbose","Verbose output",verbose::iterations);
    m_options.addInt("MaxIterations","Maximum number of iterations",25);
    m_options.addReal("Tolerance","Relative Tolerance",1e-6);
};

template <class T>
void gsStaticSolver<T>::getOptions() const
{
    m_maxIterations = m_options.getInt("MaxIterations");
    m_tolerance = m_options.getReal("Tolerance");
    m_verbose = m_options.getInt("Verbose");
};

} // namespace gismo
