 /** @file gsModalSolver.hpp

    @brief Performs linear modal analysis given a matrix or functions of a matrix

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#include <typeinfo>
#pragma once

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
void gsModalSolver<T>::computeSparse(T shift, index_t number, Spectra::SortRule selectionRule, Spectra::SortRule sortRule)
{
    if (m_verbose) { gsInfo<<"Solving eigenvalue problem" ; }
    gsSparseMatrix<T> lhs = m_stiffness-shift*m_mass;
    gsSpectraGenSymSolver<gsSparseMatrix<T>,Spectra::GEigsMode::Cholesky> solver(lhs,m_mass,number,2*number);
    if (m_verbose) { gsInfo<<"." ; }
    solver.init();
    if (m_verbose) { gsInfo<<"." ; }
    /*
        Compute with Spectra
        first argument:     sorting rule for SELECTION of the eigenvalues
        second argument:    max iterations
        third argument:     tolerance
        fourth argument:    rule for SORTING the output
    */
    solver.compute(selectionRule,1000,1e-10,sortRule);
    if (m_verbose) { gsInfo<<"." ; }
    m_values  = solver.eigenvalues();
    m_values.array() += shift;
    if (m_verbose) { gsInfo<<"." ; }
    m_vectors = solver.eigenvectors();
    if (m_verbose) { gsInfo<<"Finished\n" ; }
};

template <class T>
std::vector<std::pair<T,gsMatrix<T>> > gsModalSolver<T>::makeMode(int k) const
{
    std::vector<std::pair<T,gsMatrix<T>> > mode;
    mode.push_back( std::make_pair( m_values.at(k), m_vectors.col(k) ) );
    return mode;
};

} // namespace gismo
