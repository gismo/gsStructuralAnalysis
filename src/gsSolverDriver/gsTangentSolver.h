 /** @file gsTangentSolver.h

    @brief Tangent solver of the solver driver: the single owner of the
    factorization of the tangent stiffness.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#pragma once
#include <string>

#include <gsCore/gsLinearAlgebra.h>
#include <gsIO/gsOptionList.h>

namespace gismo
{

/**
   @brief Owns one gsSparseSolver (chosen by name) and the factorization of
   the tangent stiffness \f$K\f$ it holds.

   The factorization is computed once per version stamp: the driver passes
   gsStepContext::Kversion to factorize(), and a repeated call with the same
   stamp does no work. With equilibration the factorized matrix is
   \f$D K D\f$ with \f$D = \mathrm{diag}(1/\sqrt{|K_{ii}|})\f$ (1 where
   \f$K_{ii} = 0\f$), and solve() undoes the scaling, so callers always see
   the system \f$K x = b\f$. The solver does not know a sign convention.

   Options (see defaultOptions()): <tt>Solver</tt> (string, default
   "SimplicialLDLT"), <tt>Equilibrate</tt> (switch, false), <tt>Symmetric</tt>
   (switch, true), <tt>CheckFinite</tt> (switch, false).

   The class holds a unique pointer and is not copyable.

   @tparam T coefficient type

   \ingroup gsStructuralAnalysis
*/
template<class T>
class gsTangentSolver
{
public:
    /// Defaults: Solver "SimplicialLDLT", Equilibrate off, Symmetric on, CheckFinite off.
    gsTangentSolver();

    /// Defaults, then setOptions(@a opt)
    explicit gsTangentSolver(const gsOptionList & opt);

    /// The options with their defaults: <tt>Solver</tt>, <tt>Equilibrate</tt>,
    /// <tt>Symmetric</tt>, <tt>CheckFinite</tt>.
    static gsOptionList defaultOptions();

    /// Merge @a opt into defaultOptions() (unknown keys are ignored), read the
    /// four keys, and apply them. The (Solver, Symmetric) pair is checked once
    /// by setSolver(name, symmetric); throws on a symmetric-only solver with
    /// Symmetric false (GISMO_ENSURE) and on an unknown solver name
    /// (GISMO_ERROR in gsSparseSolver::get()). Invalidates the factorization.
    void setOptions(const gsOptionList & opt);

    /// Choose the solver by name and declare whether \f$K\f$ is symmetric.
    /// Throws if the name is unknown to gsSparseSolver::get(), or if
    /// @a symmetric is false and the solver is symmetric-only
    /// (SimplicialLDLT, SimplicialLLT, CGDiagonal, CGIdentity, PardisoLDLT,
    /// PardisoLLT). To change both the solver and the symmetry (for example
    /// to "LU", unsymmetric) use this overload or setOptions(), never
    /// setSymmetric(false) first. Invalidates the factorization.
    void setSolver(const std::string & name, bool symmetric);

    /// Choose the solver by name, keeping the current symmetry flag; the
    /// resulting pair is checked as in setSolver(name, symmetric).
    void setSolver(const std::string & name);

    /// Switch the symmetric diagonal equilibration on or off. Invalidates the factorization.
    void setEquilibrate(bool on);

    /// Declare whether \f$K\f$ is symmetric. Throws if @a sym is false and
    /// the current solver is symmetric-only; to leave a symmetric-only solver
    /// use setSolver(name, symmetric). inertia() and pivots() require true.
    /// Invalidates the factorization.
    void setSymmetric(bool sym);

    /// Check the finiteness of the solution in solve() (default off)
    void setCheckFinite(bool on);

    /// Factorize @a K unless a valid factorization of the same @a version is
    /// held. Without equilibration the solver factorizes @a K itself; with it,
    /// \f$D K D\f$ (extra O(nnz(K)) work and one copy of K). Every compute call
    /// is counted by numFactorizations(), failed ones included.
    /// @param[in] K tangent stiffness, n x n
    /// @param[in] version version stamp of @a K (any value, 0 included)
    /// @return false if the solver reports failure (the factorization is then invalid)
    /// Cost: zero for a repeated version, otherwise the cost of the chosen solver's factorization.
    bool factorize(const gsSparseMatrix<T> & K, index_t version);

    /// True iff a valid factorization is held
    bool factorized() const { return m_factorized; }

    /// Version stamp of the held factorization; meaningful iff factorized()
    index_t version() const { return m_version; }

    /// Drop the factorization: the next factorize() recomputes
    void invalidate() { m_factorized = false; }

    /// Number of factorizations computed since construction or resetCounter()
    index_t numFactorizations() const { return m_numFactorizations; }

    /// Set numFactorizations() to zero
    void resetCounter() { m_numFactorizations = 0; }

    /// Solve \f$K x = b\f$ with the held factorization.
    /// @param[in] b right-hand side, size n
    /// @param[out] x solution, size n
    /// @return false if no valid factorization is held, if b.size() differs
    /// from n, or if CheckFinite is on and @a x is not finite
    /// Cost: two triangular solves plus O(n).
    bool solve(const gsVector<T> & b, gsVector<T> & x) const;

    /// Inertia of \f$K\f$: the numbers of positive, negative and zero pivots of
    /// the LDL^T factorization. By Sylvester's law of inertia these are the
    /// counts of the eigenvalues of \f$K\f$, also under equilibration (\f$D K
    /// D\f$ with positive diagonal \f$D\f$ is a congruence) and under the
    /// fill-reducing permutation.
    /// @return false unless a factorization is held, Symmetric is on, and the
    /// solver is an LDLT solver (SimplicialLDLT or SimplicialLLT)
    /// Cost: O(n).
    bool inertia(index_t & npos, index_t & nneg, index_t & nzero) const;

    /// Pivots \f$D\f$ of the LDL^T factorization, size n. Without
    /// equilibration, D.minCoeff() is the determinant-stability indicator of
    /// gsALMBase. With equilibration the pivots belong to the scaled matrix:
    /// their signs are those of \f$K\f$, their values are not.
    /// @return false under the preconditions of inertia()
    /// Cost: O(n).
    bool pivots(gsVector<T> & D) const;

private:
    gsTangentSolver(const gsTangentSolver &);
    gsTangentSolver & operator=(const gsTangentSolver &);

    typename gsSparseSolver<T>::SimplicialLDLT * ldlt() const;

    typename gsSparseSolver<T>::uPtr m_solver;
    std::string m_name;
    gsVector<T> m_d;
    index_t m_n;
    index_t m_version;
    bool m_factorized;
    bool m_equilibrate;
    bool m_symmetric;
    bool m_checkFinite;
    index_t m_numFactorizations;
};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsTangentSolver.hpp)
#endif
