 /** @file gsBucklingSolver.h

    @brief Performs linear buckling analysis given a matrix or functions of a matrix

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#include <gsStructuralAnalysis/src/gsEigenSolvers/gsEigenProblemBase.h>

#include <gsIO/gsFileData.h>
#include <gsIO/gsReadFile.h>

#include <algorithm>
#include <limits>
#include <vector>

#pragma once

namespace gismo
{

/**
    @brief Performs linear buckling analysis by posing and solving
    the *definite* pencil K_G v = nu K_L v.

    The bifurcation condition for a linearly pre-stressed structure is

        [ K_L + lambda K_G ] v = 0 ,   K_G := K_NL(w) - K_L ,

    with K_L the linear stiffness, w the pre-buckling state and K_G the
    geometric (initial-stress) stiffness evaluated at w. Rearranging as
    K_G v = nu K_L v with nu = -1/lambda gives the buckling factor
    lambda = -1/nu; nu < 0 iff lambda > 0.

    This class poses exactly that pencil -- K_G first, K_L second -- and
    hands it to gsEigen::GeneralizedSelfAdjointEigenSolver::compute(A,B),
    which Cholesky-factors *B*. Posing it the other way round
    (K_L v = mu K_G v, mu = 1/nu) would factor K_G, which is indefinite by
    construction under any pre-stress state that has both a tension and a
    compression principal direction (e.g. shear); that factorization is
    invalid, and Eigen's LLT does not report it as a failure -- it stops at
    the first non-positive pivot and returns, leaving the eigenproblem
    silently wrong instead of caught. Posing K_G v = nu K_L v instead means
    only the symmetric positive definite K_L is ever factored; this class
    additionally certifies that positive-definiteness before every solve
    (see checkDefinite()) and reports a non-Success gsStatus rather than
    solving against an unverified second matrix.

    The swapped pencil is also more accurate: a mode with negligible
    pre-stress coupling has nu ~ 0, i.e. |lambda| -> infinity, landing
    harmlessly at the far end of the buckling-factor spectrum. In the
    K_L-first formulation the same mode is a spurious near-zero eigenvalue --
    a fictitious low buckling load.

    values(), value(k), vectors(), vector(k) and mode(k) all report the
    buckling factor lambda = -1/nu, never the raw pencil eigenvalue nu,
    ascending in |lambda| (the most critical mode first); m_vectors is
    permuted identically, so mode(k) always pairs lambda_k with its own
    eigenvector. A mode whose |nu| is at or below a relative floor (1e-8
    times the largest |nu| of the solve) is the |lambda| -> infinity case
    above: it is assigned a sentinel value
    (std::numeric_limits<T>::infinity() where the scalar type has one, the
    finite std::numeric_limits<T>::max() otherwise, since infinity() is 0 for
    some gsMultiPrecision scalars) and always sorts last -- never as a
    spurious near-zero lambda.

    The `shift` option (inherited from gsEigenProblemBase) shifts the
    (K_G, K_L) pencil that is actually factored, i.e. it acts on the raw
    eigenvalue nu, not on the reported buckling factor lambda: shifting a
    pencil is not a shift of its reciprocal, so no lambda-space shift is
    offered.

    gsSpectra's sparse eigensolvers, when enabled, rank and select modes by
    nu, not by lambda. Because lambda = -1/nu, the buckling-critical modes
    (smallest |lambda|) are the *largest*-|nu| modes -- the opposite end of
    the spectrum from the inherited default selectionRule/sortRule = 4
    (SmallestMagn). A caller driving computeSparse() must therefore set
    "selectionRule" and "sortRule" to 0 (LargestMagn) to retrieve the
    physically relevant modes, and must use the Cholesky Spectra mode
    ("solver" option 0), which Cholesky-factors the *second* matrix of the
    pencil (here m_B == K_L, symmetric positive definite) -- exactly the
    matrix this class's own dense path also factors. The Buckling mode
    ("solver" option 3) must *not* be used with this pencil: it requires the
    *first* matrix to be positive semi-definite, and this class's first
    matrix is m_A == K_G, indefinite by construction. gsSpectra is not
    enabled in every build, so this path is documented here rather than
    exercised by this class's own tests.

    computePower() is refused rather than inherited: the base-class power
    method inverts the *first* matrix of the pencil, which for this class is
    K_G -- exactly the matrix that is singular for the negligible-pre-stress
    modes this pencil is designed to push to the far end of the spectrum.

    \tparam T           coefficient type

    \ingroup gsEigenSolvers
*/
template <class T>
class gsBucklingSolver : public gsEigenProblemBase<T>
{
protected:

    typedef gsEigenProblemBase<T> Base;

    typedef std::function < bool ( gsVector<T> const &,                         gsSparseMatrix<T> & ) > Jacobian_t;
    typedef std::function < bool ( gsVector<T> const &, gsVector<T> const &,    gsSparseMatrix<T> & ) > dJacobian_t;

public:

    /**
     * @brief      Constructor
     *
     * @param      linear     The linear stiffness matrix
     * @param      rhs        The external force vector for linearization
     * @param      nonlinear  The Jacobian
     * @param[in]  scaling    A scaling factor (optional)
     */
    gsBucklingSolver(   gsSparseMatrix<T> &linear,
                        gsVector<T> &rhs,
                        Jacobian_t  &nonlinear,
                        T scaling = 1.0) :
    m_rhs(rhs),
    m_nonlinear(nonlinear),
    m_scaling(scaling)
    {
        m_A = linear;
        m_dnonlinear = [this](gsVector<T> const & x, gsVector<T> const & /*dx*/, gsSparseMatrix<T> & m) -> bool
        {
            return m_nonlinear(x,m);
        };

        m_solver = gsSparseSolver<T>::get( "SimplicialLDLT" );

        m_status = this->initializeMatrix();
    }

  /**
   * @brief      Constructor
   *
   * @param      linear     The linear stiffness matrix
   * @param      rhs        The external force vector for linearization
   * @param      nonlinear  The Jacobian taking the solution and the update as argument
   * @param[in]  scaling    A scaling factor (optional)
   */
  gsBucklingSolver(     gsSparseMatrix<T> &linear,
                        gsVector<T> &rhs,
                        dJacobian_t &dnonlinear,
                        T scaling = 1.0) :
    m_rhs(rhs),
    m_dnonlinear(dnonlinear),
    m_scaling(scaling)
    {
        m_A = linear;

        m_solver = gsSparseSolver<T>::get( "SimplicialLDLT" );

        m_status = this->initializeMatrix();
    }


    /**
    * @brief      Constructor
    *
    * @param      linear     The linear stiffness matrix
    * @param      nonlinear  The Jacobian which has already been assembled
    */
    gsBucklingSolver(     gsSparseMatrix<T> &linear,
                        gsSparseMatrix<T> &nonlinear )
    {
        m_A = linear;
        m_B = nonlinear-m_A;               // m_B == K_G = K_NL - K_L
        m_A.swap(m_B);                     // m_A == K_G, m_B == K_L

        // This constructor performs no pre-buckling solve, so unlike
        // initializeMatrix() it has no factorization of K_L to reuse: it
        // must compute one here purely to certify positive-definiteness.
        m_solver = gsSparseSolver<T>::get( "SimplicialLDLT" );
        m_solver->compute(m_B);
        m_status = checkDefinite();
    }

  // todo: add solver option
  // gsOptionList defaultOptions()
  // {
  //   gsOptionList options;
  //   options = Base::defaultOptions()
  //   options.addString("Solver","Specify the sparse solver","SimplicialLDLT");
  //   return options;
  // }

    gsStatus compute()
    {
        if (m_status != gsStatus::Success)
            return m_status;

        m_status = Base::compute();
        if (m_status == gsStatus::Success)
            convertToLambda();
        return m_status;
    }

    gsStatus computeSparse(const index_t number = 10)
    {
        if (m_options.getInt("solver") != 0 )
            gsWarn<<"It is highly recommended to use the Cholesky solver (option 0): "
                    "it factors K_L, the second (positive definite) matrix of this "
                    "class's pencil. The Buckling solver (option 3) requires the "
                    "*first* matrix positive semi-definite, which K_G is not.\n";

        if (m_status != gsStatus::Success)
            return m_status;

        m_status = Base::computeSparse(number);
        if (m_status == gsStatus::Success)
            convertToLambda();
        return m_status;
    }

    /// Refused: the base-class power method inverts m_A, which for this
    /// class is K_G -- indefinite by construction and singular for exactly
    /// the negligible-pre-stress modes this pencil pushes to the far end of
    /// the spectrum (see the class documentation). Use compute() or
    /// computeSparse() instead.
    gsStatus computePower()
    {
        gsWarn<<"gsBucklingSolver::computePower() is not supported: it would "
                "invert K_G, the indefinite (and, for negligible-pre-stress "
                "modes, singular) matrix of the pencil. Use compute() or "
                "computeSparse() instead.\n";
        m_status = gsStatus::SolverError;
        return m_status;
    }


protected:

    gsStatus initializeMatrix()
    {
        bool verbose = m_options.getSwitch("verbose");
        if (verbose) { gsInfo<<"Computing matrices" ; }
        m_solver->compute(m_A);            // m_A == K_L here
        if (verbose) { gsInfo<<"." ; }
        m_solVec = m_solver->solve(m_scaling*m_rhs);
        if (verbose) { gsInfo<<"." ; }
        try
        {
            m_dnonlinear(m_solVec,gsVector<T>::Zero(m_solVec.rows()),m_B);
        }
        catch (...)
        {
            m_status = gsStatus::AssemblyError;
            return m_status;
        }
        m_B -= m_A;                       // m_B == K_G = K_NL - K_L
        if (verbose) { gsInfo<<"." ; }

        m_A.swap(m_B);                     // m_A == K_G, m_B == K_L

        // m_solver above already holds the LDLT factorization of the
        // pre-swap m_A (== K_L); checkDefinite() reuses it rather than
        // factoring K_L a second time.
        m_status = checkDefinite();

        if (verbose) { gsInfo<<"Finished\n" ; }

        return m_status;
    }

    /// Certifies that m_B (== K_L after the pencil swap) is symmetric
    /// positive definite, reusing the sparse LDLT factorization m_solver
    /// already holds. LDLT accepts any invertible symmetric matrix and
    /// reports failure (a non-Success info()) only on an exact zero pivot,
    /// so it does not by itself certify positive-definiteness; the
    /// additional condition is that every entry of the factorization's
    /// diagonal (vectorD()) is strictly positive. Both checks read data the
    /// factorization already produced -- O(n) beyond the sparse
    /// factorization itself, never a second, dense O(n^3) check.
    gsStatus checkDefinite()
    {
        typedef typename gsSparseSolver<T>::SimplicialLDLT LDLT;
        const LDLT * ldlt = static_cast<const LDLT*>(m_solver.get());

        const bool isSPD = (m_solver->info()==gsEigen::Success) &&
                            (ldlt->vectorD().array() > 0).all();

        if (!isSPD)
        {
            gsWarn<<"gsBucklingSolver: the linear stiffness matrix K_L (the "
                    "second matrix of the pencil) is not symmetric positive "
                    "definite; the buckling eigenproblem cannot be solved.\n";
            m_status = gsStatus::SolverError;
        }
        else
            m_status = gsStatus::Success;

        return m_status;
    }

    /// Converts the raw pencil eigenvalues nu held in m_values (from
    /// K_G v = nu K_L v) into buckling factors lambda = -1/nu, in place, and
    /// permutes m_values together with the columns of m_vectors into
    /// ascending |lambda| -- the order mode(k) relies on to pair lambda_k
    /// with its own eigenvector. Called once, at the end of a successful
    /// solve, never lazily inside values(): a lazy conversion would
    /// double-convert on a second call and would force values() to be
    /// non-const.
    ///
    /// A mode with |nu| at or below a relative floor has negligible
    /// pre-stress coupling: dividing it out would either overflow or, worse,
    /// round to a spurious near-zero lambda that looks like a critical mode.
    /// Such modes are assigned a sentinel instead and are kept at the end of
    /// the ordering by an explicit flag, never by a numeric comparison
    /// against the sentinel itself.
    void convertToLambda()
    {
        const index_t n = m_values.rows();
        if (n==0)
            return;

        const T nuFloor = static_cast<T>(1e-8) * m_values.cwiseAbs().maxCoeff();
        const T inf = std::numeric_limits<T>::has_infinity
                    ? std::numeric_limits<T>::infinity()
                    : std::numeric_limits<T>::max();

        std::vector<bool> floored(n);
        for (index_t i = 0; i!=n; ++i)
        {
            const T nu = m_values(i,0);
            floored[i] = ( math::abs(nu) <= nuFloor );
            m_values(i,0) = floored[i] ? inf : ( -static_cast<T>(1)/nu );
        }

        std::vector<index_t> perm(n);
        for (index_t i = 0; i!=n; ++i) perm[i] = i;
        std::sort(perm.begin(), perm.end(),
            [this,&floored](index_t a, index_t b) -> bool
            {
                if (floored[a] != floored[b])
                    return floored[b];             // floored modes sort last
                if (floored[a])
                    return a < b;                   // stable order among floored modes
                return math::abs(m_values(a,0)) < math::abs(m_values(b,0));
            });

        gsMatrix<T> values(n,1);
        gsMatrix<T> vectors(m_vectors.rows(),n);
        for (index_t k = 0; k!=n; ++k)
        {
            values(k,0)    = m_values(perm[k],0);
            vectors.col(k) = m_vectors.col(perm[k]);
        }
        m_values  = give(values);
        m_vectors = give(vectors);
    }

protected:

    using Base::m_A;
    const gsVector<T> m_rhs;
    const Jacobian_t m_nonlinear;
         dJacobian_t m_dnonlinear;
    T m_scaling;
    using Base::m_B;


    /// Linear solver employed
    mutable typename gsSparseSolver<T>::uPtr m_solver;
    gsVector<> m_solVec;

    using Base::m_options;

    using Base::m_status;

    using Base::m_values;
    using Base::m_vectors;
};

} // namespace gismo
