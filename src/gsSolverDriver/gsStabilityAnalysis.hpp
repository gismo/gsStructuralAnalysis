 /** @file gsStabilityAnalysis.hpp

    @brief Implementation of gsStabilityAnalysis.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#pragma once
#include <gsStructuralAnalysis/src/gsSolverDriver/gsStabilityAnalysis.h>

#include <algorithm>
#include <cmath>

namespace gismo
{

template<class T>
gsStabilityAnalysis<T>::gsStabilityAnalysis(index_t method)
: m_method(method), m_indicator(0), m_negatives(0)
{ }

template<class T>
void gsStabilityAnalysis<T>::setMethod(index_t method)
{
    m_method = method;
}

template<class T>
index_t gsStabilityAnalysis<T>::getMethod() const
{
    return m_method;
}

template<class T>
bool gsStabilityAnalysis<T>::compute(const gsSparseMatrix<T> & K, T shift, gsVector<T> & values)
{
    if (K.rows()==0)
        return false;

    if (m_method == method::Determinant)
    {
        typename gsSparseSolver<T>::SimplicialLDLT ldlt;
        ldlt.compute(K);
        if (ldlt.info()!=gsEigen::ComputationInfo::Success)
            return false;
        values = ldlt.vectorD();
    }
    else if (m_method == method::Eigenvalue)
    {
#ifdef gsSpectra_ENABLED
        index_t number = std::min(static_cast<index_t>(std::floor(K.cols()/5.)),10);

        // With shift: shift-invert converges more easily than the unshifted
        // smallest-eigenvalue iteration, but requires a shift.
        gsSpectraSymShiftSolver<gsSparseMatrix<T>> es(K,number,5*number,shift);
        es.init();
        es.compute(Spectra::SortRule::LargestAlge,1000,1e-6,Spectra::SortRule::SmallestAlge);
        if (es.info()!=Spectra::CompInfo::Successful)
        {
            gsWarn<<"Spectra did not converge!\n"; // The number of Lanczos vectors (5*number) can be too low.
            return false;
        }
        values = es.eigenvalues();
#else
        GISMO_UNUSED(shift);
        gsEigen::SelfAdjointEigenSolver<gsMatrix<T>> es(K);
        values = es.eigenvalues();
#endif
    }
    else if (m_method == method::Nothing)
    {
        values = gsVector<T>::Zero(K.rows());
    }
    else
        return false;

    return classify(values);
}

template<class T>
bool gsStabilityAnalysis<T>::classify(const gsVector<T> & values)
{
    if (values.size()==0)
        return false;

    m_indicator = values.minCoeff();
    m_negatives = 0;
    for (index_t i = 0; i!=values.size(); ++i)
        if (values[i] < 0)
            ++m_negatives;
    return true;
}

template<class T>
T gsStabilityAnalysis<T>::indicator() const
{
    return m_indicator;
}

template<class T>
index_t gsStabilityAnalysis<T>::negatives() const
{
    return m_negatives;
}

template<class T>
index_t gsStabilityAnalysis<T>::stability() const
{
    return sign(m_indicator);
}

template<class T>
index_t gsStabilityAnalysis<T>::sign(T indicator)
{
    return (indicator < 0) ? -1 : 1;
}

} // namespace gismo
