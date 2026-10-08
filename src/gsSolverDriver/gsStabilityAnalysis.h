 /** @file gsStabilityAnalysis.h

    @brief Stability analysis of a symmetric tangent stiffness: eigenvalues or
    LDL^T pivots, and their classification.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>
#ifdef gsSpectra_ENABLED
#include <gsSpectra/gsSpectra.h>
#endif

namespace gismo
{

/**
   @brief Computes and classifies the stability of a symmetric tangent
   stiffness matrix \f$K\f$.

   The stability of an equilibrium is read from the signs of the spectrum of
   \f$K\f$: it is stable if all eigenvalues are non-negative. Three methods are
   available (see method):
   - Eigenvalue: the eigenvalues of \f$K\f$. With gsSpectra_ENABLED, the
     min(floor(n/5),10) eigenvalues nearest above a given shift are computed
     by shift-invert Spectra iterations; they are the smallest ones only if
     the shift lies below the whole spectrum, and eigenvalues below the shift
     are not seen. Non-convergence reports false, but Spectra itself throws
     std::invalid_argument for n < 5, where no eigenvalue is requested.
     Without Spectra, the full spectrum is computed with a dense
     gsEigen::SelfAdjointEigenSolver, which costs O(n^3) time and O(n^2) memory
     for every call.
   - Determinant: the pivots \f$D\f$ of the factorization
     \f$K = P^T L D L^T P\f$. By Sylvester's law of inertia, \f$D\f$ has as many
     negative entries as \f$K\f$ has negative eigenvalues, so the minimum
     and the sign count of the pivots classify \f$K\f$, but the pivots are not
     the eigenvalues. The cost is one sparse LDL^T factorization per call.
   - Nothing: a zero vector of the size of \f$K\f$ (always "stable").

   The class holds no solver: the factorization of the Determinant method is
   a local variable, so no factor is retained between calls and the object is
   copyable. Failures are reported by the return value and the class throws
   nothing itself; with gsSpectra_ENABLED, the Eigenvalue method can propagate
   an exception thrown by Spectra (see above).

   @tparam T coefficient type

   \ingroup gsStructuralAnalysis
*/
template<class T>
class gsStabilityAnalysis
{
public:
    /// Methods; the values equal those of gsALMBase::bifmethod and gsStaticBase::stabmethod
    struct method
    {
        enum type { Nothing = -1, Determinant = 0, Eigenvalue = 1 };
    };

    /// Constructor with the \a method (a value of method::type)
    explicit gsStabilityAnalysis(index_t method = method::Eigenvalue);

    /// Stores the \a method (a value of method::type); its validity is checked by compute()
    void setMethod(index_t method);

    /// Returns the stored method
    index_t getMethod() const;

    /**
       @brief Fills \a values from \a K with the stored method and classifies them
       with classify().

       - Determinant: factorizes \a K with a local SimplicialLDLT (released on
         return) and copies the pivots; returns false if the factorization fails.
       - Eigenvalue: with gsSpectra_ENABLED, a shift-invert Spectra computation
         around \a shift, returning false (and warning) if Spectra does not
         converge; Spectra throws std::invalid_argument if K.rows() < 5.
         Otherwise a dense eigendecomposition, for which \a shift is unused.
       - Nothing: \a values is the zero vector of size K.rows().

       For an unknown method or an empty \a K, the function returns false and leaves
       \a values and the stored results untouched. If it returns false for another
       reason, \a values may be modified but the stored results are untouched.

       @param[in]  K       symmetric tangent stiffness
       @param[in]  shift   shift of the Spectra shift-invert iterations; it must lie
                           below the smallest eigenvalue of interest
       @param[out] values  eigenvalues or pivots of \a K
       @return true on success
    */
    bool compute(const gsSparseMatrix<T> & K, T shift, gsVector<T> & values);

    /**
       @brief Classifies a vector of values (eigenvalues, or LDL^T pivots computed by
       the caller): indicator() becomes the smallest value, negatives() the number
       of negative values and stability() its sign.

       Costs O(n). For an empty \a values, the function returns false and leaves
       the stored results untouched.
    */
    bool classify(const gsVector<T> & values);

    /// Smallest value of the last successful compute() or classify()
    T indicator() const;

    /// Number of negative values of the last successful compute() or classify().
    /// For the Determinant method (Sylvester's law of inertia) and the dense
    /// Eigenvalue method this is the number of negative eigenvalues of K; with
    /// gsSpectra_ENABLED it counts only the negatives among the computed
    /// eigenvalues (at most ten), a lower bound.
    index_t negatives() const;

    /// +1 if indicator() >= 0, -1 if indicator() < 0
    index_t stability() const;

    /// Returns -1 if \a indicator < 0 and +1 otherwise. Trap: NaN is not
    /// negative and therefore maps to +1.
    static index_t sign(T indicator);

private:
    index_t m_method;
    T       m_indicator;
    index_t m_negatives;
};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsStabilityAnalysis.hpp)
#endif
