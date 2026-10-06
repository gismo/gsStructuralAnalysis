 /** @file gsTangentSolver.hpp

    @brief Implementation of gsTangentSolver.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#pragma once
#include <gsStructuralAnalysis/src/gsSolverDriver/gsTangentSolver.h>

namespace gismo
{

namespace internal
{
/// True for the solvers of gsSparseSolver::get() that need a symmetric matrix
inline bool gsTangentSolverSymmetricOnly(const std::string & name)
{
    return name == "SimplicialLDLT" || name == "SimplicialLLT"
        || name == "CGDiagonal"     || name == "CGIdentity"
        || name == "PardisoLDLT"    || name == "PardisoLLT";
}
} // namespace internal

template<class T>
gsTangentSolver<T>::gsTangentSolver()
: m_n(0), m_version(0), m_factorized(false),
  m_equilibrate(false), m_symmetric(true), m_checkFinite(false),
  m_numFactorizations(0)
{
    setSolver("SimplicialLDLT", true);
}

template<class T>
gsTangentSolver<T>::gsTangentSolver(const gsOptionList & opt)
: m_n(0), m_version(0), m_factorized(false),
  m_equilibrate(false), m_symmetric(true), m_checkFinite(false),
  m_numFactorizations(0)
{
    setOptions(opt);
}

template<class T>
gsOptionList gsTangentSolver<T>::defaultOptions()
{
    gsOptionList opt;
    opt.addString("Solver", "Sparse linear solver for the tangent", "SimplicialLDLT");
    opt.addSwitch("Equilibrate", "Symmetric diagonal scaling of the tangent before factorization", false);
    opt.addSwitch("Symmetric", "The tangent is symmetric", true);
    opt.addSwitch("CheckFinite", "Report a non-finite solution as a failed solve", false);
    return opt;
}

template<class T>
void gsTangentSolver<T>::setOptions(const gsOptionList & opt)
{
    gsOptionList o = defaultOptions();
    o.update(opt, gsOptionList::ignoreIfUnknown);
    m_equilibrate = o.getSwitch("Equilibrate");
    m_checkFinite = o.getSwitch("CheckFinite");
    setSolver(o.getString("Solver"), o.getSwitch("Symmetric"));
}

template<class T>
void gsTangentSolver<T>::setSolver(const std::string & name, bool symmetric)
{
    GISMO_ENSURE(symmetric || !internal::gsTangentSolverSymmetricOnly(name),
                 "gsTangentSolver: solver '" << name << "' needs a symmetric tangent; "
                 "declare Symmetric = true or choose an unsymmetric solver such as LU");
    m_solver = gsSparseSolver<T>::get(name);
    m_name = name;
    m_symmetric = symmetric;
    invalidate();
}

template<class T>
void gsTangentSolver<T>::setSolver(const std::string & name)
{
    setSolver(name, m_symmetric);
}

template<class T>
void gsTangentSolver<T>::setEquilibrate(bool on)
{
    m_equilibrate = on;
    invalidate();
}

template<class T>
void gsTangentSolver<T>::setSymmetric(bool sym)
{
    GISMO_ENSURE(sym || !internal::gsTangentSolverSymmetricOnly(m_name),
                 "gsTangentSolver: solver '" << m_name << "' needs a symmetric tangent; "
                 "change the solver and the symmetry together with setSolver(name, symmetric)");
    m_symmetric = sym;
    invalidate();
}

template<class T>
void gsTangentSolver<T>::setCheckFinite(bool on)
{
    m_checkFinite = on;
}

template<class T>
bool gsTangentSolver<T>::factorize(const gsSparseMatrix<T> & K, index_t version)
{
    if (m_factorized && version == m_version)
        return true;

    m_factorized = false;
    m_n = K.rows();
    if (m_equilibrate)
    {
        m_d = K.diagonal().cwiseAbs();
        for (index_t i = 0; i != m_d.size(); ++i)
            m_d(i) = m_d(i) > 0 ? T(1) / math::sqrt(m_d(i)) : T(1);
        gsSparseMatrix<T> Ks = m_d.asDiagonal() * K * m_d.asDiagonal();
        m_solver->compute(Ks);
    }
    else
        m_solver->compute(K);
    ++m_numFactorizations;

    if (m_solver->info() != gsEigen::ComputationInfo::Success)
        return false;

    m_version = version;
    m_factorized = true;
    return true;
}

template<class T>
bool gsTangentSolver<T>::solve(const gsVector<T> & b, gsVector<T> & x) const
{
    if (!m_factorized || b.size() != m_n)
        return false;

    if (m_equilibrate)
    {
        gsMatrix<T> rhs = m_d.asDiagonal() * b;
        gsMatrix<T> y = m_solver->solve(rhs);
        x = m_d.asDiagonal() * y.col(0);
    }
    else
    {
        gsMatrix<T> y = m_solver->solve(b);
        x = y.col(0);
    }
    return m_checkFinite ? x.allFinite() : true;
}

template<class T>
typename gsSparseSolver<T>::SimplicialLDLT * gsTangentSolver<T>::ldlt() const
{
    if (!m_factorized || !m_symmetric)
        return NULL;
    return dynamic_cast<typename gsSparseSolver<T>::SimplicialLDLT*>(m_solver.get());
}

template<class T>
bool gsTangentSolver<T>::inertia(index_t & npos, index_t & nneg, index_t & nzero) const
{
    typename gsSparseSolver<T>::SimplicialLDLT * s = ldlt();
    if (!s)
        return false;
    const gsVector<T> D = s->vectorD();
    npos = nneg = nzero = 0;
    for (index_t i = 0; i != D.size(); ++i)
    {
        if      (D(i) > 0) ++npos;
        else if (D(i) < 0) ++nneg;
        else               ++nzero;
    }
    return true;
}

template<class T>
bool gsTangentSolver<T>::pivots(gsVector<T> & D) const
{
    typename gsSparseSolver<T>::SimplicialLDLT * s = ldlt();
    if (!s)
        return false;
    D = s->vectorD();
    return true;
}

} // namespace gismo
