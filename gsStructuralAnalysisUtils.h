/** @file gsStructuralAnalysisUtils.h

    @brief Provides utilities for the gsStructuralAnalysis class

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#pragma once

namespace gismo
{

//  _computeStability(const gsVector<T> x) const
// {
// gsVector<T> stabilityVec;
// gsSparseMatrix<T> jacMat = m_nonlinear(x);

// #ifdef GISMO_WITH_SPECTRA
// index_t number = std::min(static_cast<index_t>(std::floor(jacMat.cols() / 3.)), 10);
// gsSpectraSymSolver<gsSparseMatrix<T>> es(jacMat, number, 3 * number);
// es.init();
// es.compute(Spectra::SortRule::LargestMagn, 1000, 1e-6);
// GISMO_ASSERT(es.info() == Spectra::CompInfo::Successful, "Spectra did not converge!"); // Reason for not converging can be due to the value of ncv (last input in the class member), which is too low.
// stabilityVec = es.eigenvalues();
// stabilityVec = stabilityVec.reverse();
// #else
// Eigen::SelfAdjointEigenSolver<gsMatrix<T>> es(m_jacMat);
// m_stabilityVec = es.eigenvalues();
// #endif

// m_indicator = stabilityVec.colwise().minCoeff()[0]; // This is required since D does not necessarily have one column.

// }

} // namespace

