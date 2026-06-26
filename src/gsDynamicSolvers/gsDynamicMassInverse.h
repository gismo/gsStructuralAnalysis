/** @file gsDynamicMassInverse.h

    @brief Inverse mass function objects for dynamic solvers

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M.M. Hodzelmans
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>

#include <string>

namespace gismo
{

/// @brief Function object that computes and applies an inverse mass matrix.
template <class T>
class gsDynamicMassInverse
{
public:
    typedef memory::shared_ptr<gsDynamicMassInverse> Ptr;

    virtual ~gsDynamicMassInverse() {}

    /// Compute or update the inverse mass representation from \a mass.
    virtual void compute(const gsSparseMatrix<T> & mass) = 0;

    /// Apply the inverse mass representation, result = M^{-1} rhs.
    virtual void apply(const gsVector<T> & rhs, gsVector<T> & result) const = 0;
};

/// @brief Default inverse mass object that explicitly stores M^{-1}.
template <class T>
class gsDynamicExplicitMassInverse : public gsDynamicMassInverse<T>
{
public:
    void compute(const gsSparseMatrix<T> & mass)
    {
        gsSparseMatrix<T> eye(mass.rows(), mass.cols());
        eye.setIdentity();
        typename gsSparseSolver<T>::LU solver(mass);
        gsMatrix<T> inverse = solver.solve(eye);
        m_massInv = inverse.sparseView();
    }

    void apply(const gsVector<T> & rhs, gsVector<T> & result) const
    {
        result = m_massInv * rhs;
    }

private:
    gsSparseMatrix<T> m_massInv;
};

/// @brief Inverse mass object backed by a sparse factorization.
template <class T>
class gsDynamicSparseMassInverse : public gsDynamicMassInverse<T>
{
public:
    explicit gsDynamicSparseMassInverse(const std::string & solver = "SimplicialLDLT")
    :
    m_solver(gsSparseSolver<T>::get(solver))
    {}

    void compute(const gsSparseMatrix<T> & mass)
    {
        m_solver->compute(mass);
    }

    void apply(const gsVector<T> & rhs, gsVector<T> & result) const
    {
        result = m_solver->solve(rhs);
    }

private:
    typename gsSparseSolver<T>::uPtr m_solver;
};

} // namespace gismo
