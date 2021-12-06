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

template <class T, Spectra::GEigsMode GEigsMode>
void gsEigenProblemBase<T,GEigsMode>::compute()
{
    if (m_verbose) { gsInfo<<"Solving eigenvalue problem" ; }
    m_eigSolver.compute(m_A,m_B);
    if (m_verbose) { gsInfo<<"." ; }
    m_values  = m_eigSolver.eigenvalues();
    if (m_verbose) { gsInfo<<"." ; }
    m_vectors = m_eigSolver.eigenvectors();
    if (m_verbose) { gsInfo<<"." ; }
    if (m_verbose) { gsInfo<<"Finished\n" ; }
};

template <class T, Spectra::GEigsMode GEigsMode>
void gsEigenProblemBase<T,GEigsMode>::compute(T shift)
{
    if (m_verbose) { gsInfo<<"Solving eigenvalue problem" ; }
    m_eigSolver.compute(m_A-shift*m_B,m_B);
    if (m_verbose) { gsInfo<<"." ; }
    m_values  = m_eigSolver.eigenvalues();
    m_values.array() += shift;
    if (m_verbose) { gsInfo<<"." ; }
    m_vectors = m_eigSolver.eigenvectors();
    if (m_verbose) { gsInfo<<"." ; }
    if (m_verbose) { gsInfo<<"Finished\n" ; }
};

template <class T, Spectra::GEigsMode GEigsMode>
template< Spectra::GEigsMode _GEigsMode>
typename std::enable_if<_GEigsMode==Spectra::GEigsMode::Cholesky ||
                        _GEigsMode==Spectra::GEigsMode::RegularInverse
                        ,
                        void>::type
gsEigenProblemBase<T,GEigsMode>::computeSparse_impl(T shift, index_t number, index_t ncvFac, Spectra::SortRule selectionRule, Spectra::SortRule sortRule)
{
#ifdef GISMO_WITH_SPECTRA
  if (m_verbose) { gsInfo<<"Solving eigenvalue problem" ; }
  gsSpectraGenSymSolver<gsSparseMatrix<T>,GEigsMode> solver(m_A-shift*m_B,m_B,number,ncvFac*number);
  if (m_verbose) { gsInfo<<"." ; }
  solver.init();
  if (m_verbose) { gsInfo<<"." ; }
  solver.compute(selectionRule,1000,1e-6,sortRule);

  if (solver.info()==Spectra::CompInfo::Successful)         { gsDebug<<"Spectra converged in "<<solver.num_iterations()<<" iterations and with "<<solver.num_operations()<<"operations. \n"; }
  else if (solver.info()==Spectra::CompInfo::NumericalIssue){ GISMO_ERROR("Spectra did not converge! Error code: NumericalIssue"); }
  else if (solver.info()==Spectra::CompInfo::NotConverging) { GISMO_ERROR("Spectra did not converge! Error code: NotConverging"); }
  else if (solver.info()==Spectra::CompInfo::NotComputed)   { GISMO_ERROR("Spectra did not converge! Error code: NotComputed");   }
  else                                                      { GISMO_ERROR("No error code known"); }

  if (m_verbose) { gsInfo<<"." ; }
  m_values  = solver.eigenvalues();
  m_values.array() += shift;
  if (m_verbose) { gsInfo<<"." ; }
  m_vectors = solver.eigenvectors();
  if (m_verbose) { gsInfo<<"Finished\n" ; }
#else
    GISMO_ERROR("Spectra not available. Tun CMake with -DGISMO_WITH_SPECTRA=ON");
#endif
}

template <class T, Spectra::GEigsMode GEigsMode>
template< Spectra::GEigsMode _GEigsMode>
typename std::enable_if<_GEigsMode==Spectra::GEigsMode::ShiftInvert ||
                        _GEigsMode==Spectra::GEigsMode::Buckling ||
                        _GEigsMode==Spectra::GEigsMode::Cayley
                        ,
                        void>::type
gsEigenProblemBase<T,GEigsMode>::computeSparse_impl(T shift, index_t number, index_t ncvFac, Spectra::SortRule selectionRule, Spectra::SortRule sortRule)
{
#ifdef GISMO_WITH_SPECTRA
  if (selectionRule == Spectra::SortRule::SmallestMagn )
    gsWarn<<"Selection Rule 'SmallestMagn' is selected, but for ShiftInvert, Buckling and Cayley it is advised to use 'LargestMagn'!\n";
  if (selectionRule == Spectra::SortRule::SmallestAlge )
    gsWarn<<"Selection Rule 'SmallestAlge' is selected, but for ShiftInvert, Buckling and Cayley it is advised to use 'LargestMagn'!\n";

  if (m_verbose) { gsInfo<<"Solving eigenvalue problem" ; }
  gsSpectraGenSymShiftSolver<gsSparseMatrix<T>,GEigsMode> solver(m_A,m_B,number,ncvFac*number,shift);
  if (m_verbose) { gsInfo<<"." ; }
  solver.init();
  if (m_verbose) { gsInfo<<"." ; }
  solver.compute(selectionRule,1000,1e-6,sortRule);

  if (solver.info()==Spectra::CompInfo::Successful)         { gsDebug<<"Spectra converged in "<<solver.num_iterations()<<" iterations and with "<<solver.num_operations()<<"operations. \n"; }
  else if (solver.info()==Spectra::CompInfo::NumericalIssue){ GISMO_ERROR("Spectra did not converge! Error code: NumericalIssue"); }
  else if (solver.info()==Spectra::CompInfo::NotConverging) { GISMO_ERROR("Spectra did not converge! Error code: NotConverging"); }
  else if (solver.info()==Spectra::CompInfo::NotComputed)   { GISMO_ERROR("Spectra did not converge! Error code: NotComputed");   }
  else                                                      { GISMO_ERROR("No error code known"); }

  if (m_verbose) { gsInfo<<"." ; }
  m_values  = solver.eigenvalues();
  if (m_verbose) { gsInfo<<"." ; }
  m_vectors = solver.eigenvectors();
  if (m_verbose) { gsInfo<<"Finished\n" ; }
#else
    GISMO_ERROR("Spectra not available. Tun CMake with -DGISMO_WITH_SPECTRA=ON");
#endif
};

template <class T, Spectra::GEigsMode GEigsMode>
void gsEigenProblemBase<T,GEigsMode>::computePower()
{
    if (m_verbose) { gsInfo<<"Solving eigenvalue problem" ; }
    gsMatrix<T> D = m_A.toDense().inverse() * (m_B);

    gsVector<T> v(D.cols());
    v.setOnes();
    gsVector<T> v_old(D.cols());
    v_old.setZero();

    index_t kmax = 100;
    real_t error,tol = 1e-5;
    for (index_t k=0; k!=kmax; k++)
    {
      v = D*v;
      v.normalize();

      error = (v-v_old).norm();

      if ( error < tol )
        break;

      v_old = v;
    }

    m_vectors = v;
    m_values =  (v.transpose() * v) / (v.transpose() * D * v);

    if (m_verbose) { gsInfo<<"Finished\n" ; }
};

template <class T, Spectra::GEigsMode GEigsMode>
std::vector<std::pair<T,gsMatrix<T>> > gsEigenProblemBase<T,GEigsMode>::makeMode(int k) const
{
    std::vector<std::pair<T,gsMatrix<T>> > mode;
    mode.push_back( std::make_pair( m_values.at(k), m_vectors.col(k) ) );
    return mode;
};

} // namespace gismo
