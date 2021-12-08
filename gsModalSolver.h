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
    @brief Performs the arc length method to solve a nonlinear equation system.

    \tparam T coefficient type

    \ingroup ThinShell
*/
template <class T, Spectra::GEigsMode GEigsMode = Spectra::GEigsMode::Cholesky>
class gsModalSolver : public gsEigenProblemBase<T,GEigsMode>
{
protected:

    typedef gsEigenProblemBase<T,GEigsMode> Base;

public:

  /// Constructor giving access to the gsShellAssembler object to create a linear system per iteration
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
