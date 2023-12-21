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

    \tparam T           coefficient type

    \ingroup gsModalSolver
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
  }

protected:

    using Base::m_A;
    using Base::m_B;
    using Base::m_options;
};


} // namespace gismo
