 /** @file gsModalSolver.h

    @brief Performs linear modal analysis given a matrix or functions of a matrix

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
class gsModalSolver
{
protected:
    typedef std::vector<std::pair<T,gsMatrix<T>> > modes_t;
public:

  /// Constructor giving access to the gsShellAssembler object to create a linear system per iteration
  gsModalSolver(        gsSparseMatrix<T> &stiffness,
                        gsSparseMatrix<T> &mass     ) :
    m_stiffness(stiffness),
    m_mass(mass)
  {
    m_verbose = false;
  }

public:

    void compute();

    gsMatrix<T> values() const { return m_values; };
    T value(int k) const { return m_values.at(k); };

    gsMatrix<T> vectors() const { return m_vectors; };
    gsMatrix<T> vector(int k) const { return m_vectors.col(k); };

    modes_t mode(int k) const {makeMode(k); return m_mode; }

protected:

    const gsSparseMatrix<T> m_stiffness;
    const gsSparseMatrix<T> m_mass;

    Eigen::GeneralizedSelfAdjointEigenSolver< gsMatrix<real_t>::Base >  m_eigSolver;

    gsVector<T> m_solVec;
    gsMatrix<T> m_values,m_vectors;

    modes_t m_mode;

    bool m_verbose;

protected:

    void makeMode(int k);

};


} // namespace gismo

//////////////////////////////////////////////////
//////////////////////////////////////////////////

namespace gismo
{

template <class T>
void gsModalSolver<T>::compute()
{
    if (m_verbose) { gsInfo<<"Solving eigenvalue problem" ; }
    m_eigSolver.compute(m_stiffness,m_mass);

    if (m_verbose) { gsInfo<<"." ; }
    m_values  = m_eigSolver.eigenvalues();
    if (m_verbose) { gsInfo<<"." ; }
    m_vectors = m_eigSolver.eigenvectors();
    if (m_verbose) { gsInfo<<"." ; }
    if (m_verbose) { gsInfo<<"Finished\n" ; }
};

template <class T>
void gsModalSolver<T>::makeMode(int k)
{
    m_mode.push_back( std::make_pair( m_values.at(k), m_vectors.col(k) ) );
};

} // namespace gismo
