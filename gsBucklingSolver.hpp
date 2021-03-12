 /** @file gsBucklingSolver.hpp

    @brief Performs linear buckling analysis given a matrix or functions of a matrix

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#pragma once


//////////////////////////////////////////////////
//////////////////////////////////////////////////

namespace gismo
{

template <class T>
void gsBucklingSolver<T>::initializeMatrix()
{
  if (m_verbose) { gsInfo<<"Computing matrices" ; }
  m_solver.compute(m_linear);
  if (m_verbose) { gsInfo<<"." ; }
  m_solVec = m_solver.solve(m_scaling*m_rhs);
  if (m_verbose) { gsInfo<<"." ; }
  m_nonlinear = m_nonlinearFun(m_solVec);
  if (m_verbose) { gsInfo<<"." ; }
  if (m_verbose) { gsInfo<<"Finished\n" ; }
};

template <class T>
void gsBucklingSolver<T>::compute()
{
    if (m_verbose) { gsInfo<<"Solving eigenvalue problem" ; }
    m_eigSolver.compute(m_linear,m_nonlinear - m_linear);
    if (m_verbose) { gsInfo<<"." ; }
    m_values  = m_eigSolver.eigenvalues();
    if (m_verbose) { gsInfo<<"." ; }
    m_vectors = m_eigSolver.eigenvectors();
    if (m_verbose) { gsInfo<<"." ; }
    if (m_verbose) { gsInfo<<"Finished\n" ; }
};

template <class T>
void gsBucklingSolver<T>::computeSparse(T shift, index_t number)
{
    if (m_verbose) { gsInfo<<"Solving eigenvalue problem" ; }
    gsSparseMatrix<T> lhs = m_linear-shift*(m_nonlinear - m_linear);
    gsSpectraGenSymSolver<gsSparseMatrix<T>,Spectra::GEigsMode::Cholesky> solver(lhs,m_nonlinear - m_linear,number,2*number);
    if (m_verbose) { gsInfo<<"." ; }
    solver.init();
    if (m_verbose) { gsInfo<<"." ; }
    solver.compute(Spectra::SortRule::SmallestMagn,1000,1e-10,Spectra::SortRule::SmallestAlge);
    if (m_verbose) { gsInfo<<"." ; }
    m_values  = solver.eigenvalues();
    m_values.array() += shift;
    if (m_verbose) { gsInfo<<"." ; }
    m_vectors = solver.eigenvectors();
    if (m_verbose) { gsInfo<<"Finished\n" ; }
};

template <class T>
void gsBucklingSolver<T>::computePower()
{
    if (m_verbose) { gsInfo<<"Solving eigenvalue problem" ; }
    gsMatrix<T> D = m_linear.toDense().inverse() * (m_nonlinear.toDense() - m_linear.toDense());

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


template <class T>
std::vector<std::pair<T,gsMatrix<T>> > gsBucklingSolver<T>::makeMode(int k) const
{
    std::vector<std::pair<T,gsMatrix<T>> > mode;
    mode.push_back( std::make_pair( m_values.at(k), m_vectors.col(k) ) );
    return mode;
};

} // namespace gismo
