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

#pragma once

namespace gismo
{

/**
    @brief Performs linear buckling analysis given a matrix or functions of a matrix

    \tparam T           coefficient type

    \ingroup gsBucklingSolver
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
        m_dnonlinear = [this](gsVector<T> const & x, gsVector<T> const & dx, gsSparseMatrix<T> & m) -> bool
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
        m_B = nonlinear-m_A;
    }

  // todo: add solver option
  // gsOptionList defaultOptions()
  // {
  //   gsOptionList options;
  //   options = Base::defaultOptions()
  //   options.addString("Solver","Specify the sparse solver","SimplicialLDLT");
  //   return options;
  // }

    gsStatus computeSparse(const index_t number = 10)
    {
        if (m_options.getInt("solver") != 3 )
            gsWarn<<"It is highly recommended to use the Buckling solver (option 4)\n";

        return Base::computeSparse(number);
    }


protected:

    gsStatus initializeMatrix()
    {
        bool verbose = m_options.getSwitch("verbose");
        if (verbose) { gsInfo<<"Computing matrices" ; }
        m_solver->compute(m_A);
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
        m_B -= m_A;
        if (verbose) { gsInfo<<"." ; }
        if (verbose) { gsInfo<<"Finished\n" ; }

        return m_status;
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
};

} // namespace gismo
