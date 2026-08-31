 /** @file gsModalSolver.h

    @brief Performs linear modal analysis given a matrix or functions of a matrix

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#include <typeinfo>
#include <gsStructuralAnalysis/src/gsEigenSolvers/gsEigenProblemBase.h>
#include <gsIO/gsOptionList.h>

#pragma once


namespace gismo
{

/**
    @brief Performs linear modal analysis given a matrix or functions of a matrix

    This class poses the generalized eigenvalue pencil K v = w M v (m_A == K,
    m_B == M, unswapped) and hands it to
    gsEigen::GeneralizedSelfAdjointEigenSolver::compute(A,B), which
    Cholesky-factors *B* == M. That factorization requires M symmetric
    positive definite, and Eigen's LLT does not report a failed factorization
    as such: it stops at the first non-positive pivot and returns, leaving
    the reduced problem silently wrong instead of caught, so the class
    certifies M's positive-definiteness itself (see checkDefinite()) before
    any solve and reports a non-Success gsStatus rather than returning
    plausible numbers for an unverified matrix.

    \tparam T           coefficient type

    \ingroup gsEigenSolvers
*/
template <class T>
class gsModalSolver : public gsEigenProblemBase<T>
{
protected:

    typedef gsEigenProblemBase<T> Base;

public:

  /**
   * @brief      Constructor
   *
   * @param      stiffness  The stiffness matrix
   * @param      mass       The mass matrix
   */
  gsModalSolver(    const gsSparseMatrix<T> &stiffness,
                    const gsSparseMatrix<T> &mass     )
  {
    m_A = stiffness;
    m_B = mass;
    m_massStatus = checkDefinite();
    if (m_massStatus != gsStatus::Success)
        m_status = m_massStatus;
  }

    gsStatus compute()
    {
        if (m_massStatus != gsStatus::Success)
        {
            m_status = m_massStatus;
            return m_status;
        }
        m_status = Base::compute();
        return m_status;
    }

    gsStatus computeSparse(const index_t number = 10)
    {
        if (m_massStatus != gsStatus::Success)
        {
            m_status = m_massStatus;
            return m_status;
        }
        m_status = Base::computeSparse(number);
        return m_status;
    }

    gsStatus computePower()
    {
        if (m_massStatus != gsStatus::Success)
        {
            m_status = m_massStatus;
            return m_status;
        }
        m_status = Base::computePower();
        return m_status;
    }

protected:

    /// Certifies that m_B (the mass matrix M) is symmetric positive
    /// definite via a sparse SimplicialLDLT factorization. LDLT accepts any
    /// invertible symmetric matrix and reports failure (a non-Success
    /// info()) only on an exact zero pivot, so info() alone does not
    /// certify positive-definiteness; the additional condition is that
    /// every entry of the factorization's diagonal (vectorD()) is strictly
    /// positive. Both reads are O(n) on data the factorization already
    /// produced, never a second, dense O(n^3) check.
    gsStatus checkDefinite()
    {
        typename gsSparseSolver<T>::SimplicialLDLT ldlt(m_B);

        const bool isSPD = (ldlt.info()==gsEigen::Success) &&
                           (ldlt.vectorD().array() > 0).all();

        if (!isSPD)
        {
            gsWarn<<"gsModalSolver: the mass matrix is not symmetric positive "
                    "definite; the modal eigenproblem cannot be solved.\n";
            return gsStatus::SolverError;
        }
        return gsStatus::Success;
    }

protected:

    using Base::m_A;
    using Base::m_B;
    using Base::m_options;
    using Base::m_status;

    /// Cached verdict of checkDefinite(), set once at construction. Kept
    /// separate from m_status (rather than latched into it the way
    /// gsBucklingSolver does) because the mass matrix cannot change after
    /// construction but the inherited "shift"/"tolerance" options can be
    /// changed between calls: a first NotConverged solve must remain
    /// retryable, which a shared latch on m_status would prevent.
    gsStatus m_massStatus;
};


} // namespace gismo
