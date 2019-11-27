 /** @file gsBucklingSolver.h

    @brief Performs linear buckling analysis given a matrix or functions of a matrix

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <typeinfo>
#pragma once


namespace gismo
{

/**
    @brief Performs the arc length method to solve a nonlinear equation system.

    \tparam T coefficient type

    \ingroup ThinShell
*/
template <class T>
class gsBucklingSolver
{
protected:
    typedef std::vector<std::pair<T,gsMatrix<T>> > modes_t;

public:

  /// Constructor giving access to the gsShellAssembler object to create a linear system per iteration
  gsBucklingSolver(     gsSparseMatrix<T> &linear,
                        gsVector<T> &rhs,
                        std::function < gsSparseMatrix<T> ( gsVector<T> const & ) > &nonlinear,
                        T scaling = 1.0) :
    m_linear(linear),
    m_rhs(rhs),
    m_nonlinearFun(nonlinear),
    m_scaling(scaling)
  {
    this->initializeMatrix();
  }

  /// Constructor giving access to the gsShellAssembler object to create a linear system per iteration
  gsBucklingSolver(     gsSparseMatrix<T> &linear,
                        gsSparseMatrix<T> &nonlinear ) :
    m_linear(Jacobian),
    m_nonlinear(Residual)
  {

  }
public:

    void compute();

    gsMatrix<T> values() const { return m_values; };
    T value(int k) const { return m_values.at(k); };

    gsMatrix<T> vector() const { return m_vectors; };
    gsMatrix<T> vector(int k) const { return m_vectors.col(k); };

    modes_t mode(int k) const {makeMode(k); return m_mode; }

protected:

    const gsSparseMatrix<T> m_linear
    gsSparseMatrix<T> m_nonlinear;
    const gsVector<T> m_rhs;
    const std::function < gsSparseMatrix<T> ( gsVector<T> const & ) > m_nonlinearFun;
    T m_scaling;

    /// Linear solver employed
    gsSparseSolver<>::SimplicialLDLT  m_solver;
    Eigen::GeneralizedSelfAdjointEigenSolver< gsMatrix<real_t>::Base >  m_eigSolver;

    gsSparseMatrix<T> m_diff;
    gsVector<T> m_solVec;
    gsMatrix<T> m_values,m_vectors;

    modes_t m_mode;

protected:

    void initializeMatrix();
    void makeMode(int k);

};


} // namespace gismo

//////////////////////////////////////////////////
//////////////////////////////////////////////////

namespace gismo
{

template <class T>
void gsBucklingSolver<T>::initializeMatrix()
{
    solver.compute(m_linear);
    m_solVec = solver.solve(scaling*m_rhs);
    m_nonlinear = m_nonlinearFun(m_solVec);
};

template <class T>
void gsBucklingSolver<T>::compute()
{
    m_eigSolver.compute(m_linear,m_nonlinear - m_linear);

    m_values  = m_eigSolver.eigenvalues();
    m_vectors = m_eigSolver.eigenvectors();
};

template <class T>
void gsModalSolver<T>::makeMode(int k)
{
    m_mode.push_back( std::make_pair( m_values.at(k), m_vectors.col(k) ) );
};


} // namespace gismo
