 /** @file gsEigenSolverBase.h

    @brief Performs linear buckling analysis given a matrix or functions of a matrix

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#include <typeinfo>
#include <gsSpectra/gsSpectra.h>
#pragma once


namespace gismo
{

/**
    @brief Performs the arc length method to solve a nonlinear equation system.

    \tparam T coefficient type

    \ingroup ThinShell
*/
template <class T>
class gsEigenSolverBase
{
protected:
    // typedef typename std::vector<std::pair<T,gsMatrix<T>> > modes_t;

public:

    virtual ~gsEigenSolverBase() {};

public:

    void verbose() {m_verbose=true; };

    void compute();
    void computeSparse(T shift = 0.0, index_t number = 10);
    void computePower();

    gsMatrix<T> values() const { return m_values; };
    T value(int k) const { return m_values.at(k); };

    gsMatrix<T> vectors() const { return m_vectors; };
    gsMatrix<T> vector(int k) const { return m_vectors.col(k); };

    std::vector<std::pair<T,gsMatrix<T>> > mode(int k) const {return makeMode(k); }

protected:

    const gsSparseMatrix<T> m_linear;
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

    bool m_verbose;

    index_t m_num;

protected:

    void initializeMatrix();
    std::vector<std::pair<T,gsMatrix<T>> > makeMode(int k) const;

};


} // namespace gismo


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsEigenSolverBase.hpp)
#endif