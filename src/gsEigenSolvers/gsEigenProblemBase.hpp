 /** @file gsEigenProblemBase.hpp

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
gsStatus gsEigenProblemBase<T>::compute()
{
    if (m_status==gsStatus::AssemblyError)
        return m_status;

    bool verbose = m_options.getSwitch("verbose");
    if (verbose) { gsInfo<<"Solving eigenvalue problem" ; }
    try
    {
        T shift = m_options.getReal("shift");
        if (shift!=0.0)
            m_eigSolver.compute(m_A-shift*m_B,m_B);
        else
            m_eigSolver.compute(m_A,m_B);

        if (verbose) { gsInfo<<"." ; }
        m_values  = m_eigSolver.eigenvalues();
        m_values.array() += shift;
        if (verbose) { gsInfo<<"." ; }
        m_vectors = m_eigSolver.eigenvectors();
        if (verbose) { gsInfo<<"." ; }
        if (verbose) { gsInfo<<"Finished\n" ; }
        m_status = gsStatus::Success;
    }
    catch (...)
    {
        m_status = gsStatus::SolverError;
    }
    return m_status;
};


template <class T>
gsStatus gsEigenProblemBase<T>::computeSparse(const index_t number)
{
    if (m_status==gsStatus::AssemblyError)
        return m_status;
    
    #ifdef gsSpectra_ENABLED
        if (m_options.getInt("solver")==0)
            return computeSparse_impl<Spectra::GEigsMode::Cholesky>(number);
        else if (m_options.getInt("solver")==1)
            return computeSparse_impl<Spectra::GEigsMode::RegularInverse>(number);
        else if (m_options.getInt("solver")==2)
            return computeSparse_impl<Spectra::GEigsMode::ShiftInvert>(number);
        else if (m_options.getInt("solver")==3)
            return computeSparse_impl<Spectra::GEigsMode::Buckling>(number);
        else if (m_options.getInt("solver")==4)
            return computeSparse_impl<Spectra::GEigsMode::Cayley>(number);
        else 
        {
            gsWarn<<"Solver index "<<m_options.getInt("solver")<<" unknown.\n";
            m_status=gsStatus::NotStarted;
            return m_status;
        }
    #else
        gsWarn<<"Sparse solver is not implemented without gsSpectra. Please compile gismo with Spectra.\n";
        m_status=gsStatus::NotStarted;
        return gsStatus::NotStarted;
    #endif
};

#ifdef gsSpectra_ENABLED
template <class T>
template< Spectra::GEigsMode _GEigsMode>
typename std::enable_if<_GEigsMode==Spectra::GEigsMode::Cholesky ||
                        _GEigsMode==Spectra::GEigsMode::RegularInverse
                        ,
                        gsStatus>::type
gsEigenProblemBase<T>::computeSparse_impl(index_t number)
{
    bool verbose = m_options.getSwitch("verbose");
    T shift = m_options.getReal("shift");

    Spectra::SortRule selectionRule = static_cast<Spectra::SortRule>(m_options.getInt("selectionRule"));
    Spectra::SortRule sortRule = static_cast<Spectra::SortRule>(m_options.getInt("sortRule"));

    index_t ncvFac = m_options.getInt("ncvFac");
    T tol = m_options.getReal("tolerance");
    if (verbose) { gsInfo<<"Solving eigenvalue problem" ; }

    gsSparseMatrix<T> Atmp;
    if (shift!=0.0)
        Atmp = m_A-shift*m_B;
    else
        Atmp = m_A;

    gsSpectraGenSymSolver<gsSparseMatrix<T>,_GEigsMode> solver(Atmp,m_B,number,ncvFac*number);

    if (verbose) { gsInfo<<"." ; }
    solver.init();
    if (verbose) { gsInfo<<"." ; }
    solver.compute(selectionRule,1000,tol,sortRule);

    if      (solver.info()==Spectra::CompInfo::Successful)
    {
        gsDebug<<"Spectra converged in "<<solver.num_iterations()<<" iterations and with "<<solver.num_operations()<<"operations. \n";
        if (verbose) { gsInfo<<"." ; }
        m_values  = solver.eigenvalues();
        m_values.array() += shift;
        if (verbose) { gsInfo<<"." ; }
        m_vectors = solver.eigenvectors();
        if (verbose) { gsInfo<<"Finished\n" ; }

        m_status = gsStatus::Success;
    }
    else if (solver.info()==Spectra::CompInfo::NotConverging) 
    {
        gsWarn<<"Spectra did not converge! Error code: NotConverging\n";
        m_status = gsStatus::NotConverged;
    }
    else if (solver.info()==Spectra::CompInfo::NumericalIssue)
    {
        gsWarn<<"Spectra did not converge! Error code: NumericalIssue\n";
        m_status = gsStatus::SolverError;
    }
    else if (solver.info()==Spectra::CompInfo::NotComputed)
    {
        gsWarn<<"Spectra did not converge! Error code: NotComputed\n";
        m_status = gsStatus::SolverError;
    }
    else
    {
        gsWarn<<"No error code known\n";
        m_status = gsStatus::OtherError;
    }
    return m_status;
}
#endif

#ifdef gsSpectra_ENABLED
template <class T>
template< Spectra::GEigsMode _GEigsMode>
typename std::enable_if<_GEigsMode==Spectra::GEigsMode::ShiftInvert ||
                        _GEigsMode==Spectra::GEigsMode::Buckling ||
                        _GEigsMode==Spectra::GEigsMode::Cayley
                        ,
                        gsStatus>::type
gsEigenProblemBase<T>::computeSparse_impl(index_t number)
{
    bool verbose = m_options.getSwitch("verbose");
    T shift = m_options.getReal("shift");

    Spectra::SortRule selectionRule = static_cast<Spectra::SortRule>(m_options.getInt("selectionRule"));
    Spectra::SortRule sortRule = static_cast<Spectra::SortRule>(m_options.getInt("sortRule"));

    if (selectionRule == Spectra::SortRule::SmallestMagn ) // =4
    gsWarn<<"Selection Rule 'SmallestMagn' is selected, but for ShiftInvert, Buckling and Cayley it is advised to use 'LargestMagn'!\n";
    if (selectionRule == Spectra::SortRule::SmallestAlge ) // =7
    gsWarn<<"Selection Rule 'SmallestAlge' is selected, but for ShiftInvert, Buckling and Cayley it is advised to use 'LargestMagn'!\n";

    index_t ncvFac = m_options.getInt("ncvFac");
    T tol = m_options.getReal("tolerance");
    if (verbose) { gsInfo<<"Solving eigenvalue problem" ; }
    gsSpectraGenSymShiftSolver<gsSparseMatrix<T>,_GEigsMode> solver(m_A,m_B,number,ncvFac*number,shift);
    if (verbose) { gsInfo<<"." ; }
    solver.init();
    if (verbose) { gsInfo<<"." ; }
    solver.compute(selectionRule,1000,tol,sortRule);

    if      (solver.info()==Spectra::CompInfo::Successful)
    {
        gsDebug<<"Spectra converged in "<<solver.num_iterations()<<" iterations and with "<<solver.num_operations()<<"operations. \n";
        if (verbose) { gsInfo<<"." ; }
        m_values  = solver.eigenvalues();
        // m_values.array() += shift;
        if (verbose) { gsInfo<<"." ; }
        m_vectors = solver.eigenvectors();
        if (verbose) { gsInfo<<"Finished\n" ; }

        m_status = gsStatus::Success;
    }
    else if (solver.info()==Spectra::CompInfo::NotConverging) 
    {
        gsWarn<<"Spectra did not converge! Error code: NotConverging\n";
        m_status = gsStatus::NotConverged;
    }
    else if (solver.info()==Spectra::CompInfo::NumericalIssue)
    {
        gsWarn<<"Spectra did not converge! Error code: NumericalIssue\n";
        m_status = gsStatus::SolverError;
    }
    else if (solver.info()==Spectra::CompInfo::NotComputed)
    {
        gsWarn<<"Spectra did not converge! Error code: NotComputed\n";
        m_status = gsStatus::SolverError;
    }
    else
    {
        gsWarn<<"No error code known\n";
        m_status = gsStatus::OtherError;
    }
    return m_status;
};
#endif

template <class T>
gsStatus gsEigenProblemBase<T>::computePower()
{
    if (m_status==gsStatus::AssemblyError)
        return m_status;

    bool verbose = m_options.getSwitch("verbose");
    if (verbose) { gsInfo<<"Solving eigenvalue problem" ; }
    gsMatrix<T> D = m_A.toDense().inverse() * (m_B);

    gsVector<T> v(D.cols());
    v.setOnes();
    gsVector<T> v_old(D.cols());
    v_old.setZero();

    index_t kmax = 100;
    real_t error,tol = m_options.getReal("tolerance");
    for (index_t k=0; k!=kmax; k++)
    {
        v = D*v;
        v.normalize();

        error = (v-v_old).norm();

        if ( error < tol )
        {
            m_status = gsStatus::Success;
            break;        
        }
        else if (k==kmax-1)
            m_status = gsStatus::NotConverged;

        v_old = v;
    }

    m_vectors = v;
    m_values =  (v.transpose() * v) / (v.transpose() * D * v);

    if (verbose) { gsInfo<<"Finished\n" ; }
    return m_status;
};

template <class T>
std::vector<std::pair<T,gsMatrix<T>> > gsEigenProblemBase<T>::makeMode(int k) const
{
    std::vector<std::pair<T,gsMatrix<T>> > mode;
    mode.push_back( std::make_pair( m_values.at(k), m_vectors.col(k) ) );
    return mode;
};

} // namespace gismo
