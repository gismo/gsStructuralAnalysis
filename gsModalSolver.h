 /** @file gsModalSolver.h

    @brief Performs linear modal analysis given a matrix or functions of a matrix

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#include <typeinfo>
#include <gsSpectra/gsSpectra.h>
#include <gsStructuralAnalysis/gsEigenProblemBase.h>
#pragma once

#include <gsSpectra/gsSpectra.h>
#include <gsIO/gsOptionList.h>

namespace gismo
{

/**
    @brief
    Performs linear modal analysis given a matrix or functions of a matrix

    Linear buckling analysis is performed by solving the following eigenvalue problem:

    \f{align*}{
        (K_L-\sigma K_{NL}(\mathbf{u}^L_h))\mathbf{v}_h = \lambda K_{NL}(\mathbf{u}_h) \mathbf{v}_h
    \f}

    Where \f$K_L\f$ is the linear stiffness matrix, \f$K_{NL}(\mathbf{u}_h)\f$ is the tangential
    stiffness matrix assembled around \f$\mathbf{u}^L_h\f$. The solution \f$\mathbf{u}^L_h\f$ is
    obtained by solving a linear problem \f$K_L\mathbf{u}^L_h = \mathbf{P}\f$. Furthermore,
    \f$\sigma\f$ is a shift and \f$(\lambda+\sigma\f)\mathbf{P}$ is the critical buckling load.
    The modeshape is represented by \f$\phi\f$.

    An example with the use of this class is in \ref gsThinShell_Buckling.cpp

    \tparam T           coefficient type

    \tparam GEigsMode   The mode for the Spectra solver

    \ingroup gsStructuralAnalysis
*/
template <class T, Spectra::GEigsMode GEigsMode = Spectra::GEigsMode::Cholesky>
class gsModalSolver : public gsEigenProblemBase<T,GEigsMode>
{
protected:

    typedef gsEigenProblemBase<T,GEigsMode> Base;

public:

  /**
   * @brief      Constructor
   *
   * @param      stiffness  The stiffness matrix
   * @param      mass       The mass matrix
   */
  gsModalSolver(        gsSparseMatrix<T> &stiffness,
                        gsSparseMatrix<T> &mass     )
  {
    m_A = stiffness;
    m_B = mass;
    m_verbose = false;
  }


protected:

    using Base::m_A;
    using Base::m_B;
    using Base::m_verbose;

};


} // namespace gismo
