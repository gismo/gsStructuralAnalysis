 /** @file gsStaticSolver.h

    @brief Performs linear modal analysis given a matrix or functions of a matrix

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

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
template <class T>
class gsStaticSolver
{
protected:
    // typedef std::vector<std::pair<T,gsMatrix<T>> > modes_t;
public:

  /// Constructor giving access to the gsShellAssembler object to create a linear system per iteration
  gsStaticSolver(   const gsSparseMatrix<T> &linear,
                    const gsVector<T> &force
                    ) :
    m_linear(linear),
    m_force(force)
  {
   m_NL = false;
   defaultOptions();
  }

private:
  /**
   * @brief      helper constructor
   *
   * @param[in]  linear      The linear matrix
   * @param[in]  force       The force vector
   * @param[in]  residual    The residual function
   */
  gsStaticSolver(   const gsSparseMatrix<T> &linear,
                    const gsVector<T> &force,
                    const std::function < gsVector<T> ( gsVector<T> const & ) > &residual
                    ) :
    m_linear(linear),
    m_force(force),
    m_nonlinear(nullptr),
    m_dnonlinear(nullptr),
    m_residual(residual)
  {
    m_NL = true;
    m_converged = false;
    defaultOptions();
  }

public:
  /**
   * @brief      { function_description }
   *
   * @param[in]  linear      The linear matrix
   * @param[in]  force       The force vector
   * @param[in]  dnonlinear  The jacobian taking the solution x and the iterative update dx
   * @param[in]  residual    The residual function
   */
  gsStaticSolver(   const gsSparseMatrix<T> &linear,
                    const gsVector<T> &force,
                    const std::function < gsSparseMatrix<T> ( gsVector<T> const &, gsVector<T> const &) > &dnonlinear,
                    const std::function < gsVector<T> ( gsVector<T> const & ) > &residual
                    ) :
    gsStaticSolver(linear,force,residual)
  {
    m_dnonlinear = dnonlinear;
  }

  /**
   * @brief      { function_description }
   *
   * @param[in]  linear     XXXX
   * @param[in]  force      XXXX
   * @param[in]  nonlinear  XXXX
   * @param[in]  residual   XXXX
   */
  gsStaticSolver(   const gsSparseMatrix<T> &linear,
                    const gsVector<T> &force,
                    const std::function < gsSparseMatrix<T> ( gsVector<T> const & ) > &nonlinear,
                    const std::function < gsVector<T> ( gsVector<T> const & ) > &residual
                    ) :
    gsStaticSolver(linear,force,residual)
  {
    m_nonlinear  = nonlinear;
    m_dnonlinear = [this](gsVector<T> const & x, gsVector<T> const & dx)
    {
        return m_nonlinear(x);
    };
  }

public:

    gsVector<T> solveLinear();
    gsVector<T> solveNonlinear();

    gsOptionList & options() {return m_options;}
    void setOptions(gsOptionList opt) { m_options = opt; } // gsOptionList opt
    void defaultOptions();
    void getOptions() const;

    int iterations() {return m_iterations;}

    gsVector<T> solution() {return m_solVec ;}

    bool converged() { return m_converged; }

    T indicator() const
    {
        _computeStability(m_solVec);
        return m_indicator;
    }

protected:
    void _computeStability(const gsVector<T> x) const;

protected:

    const gsSparseMatrix<T> & m_linear;
    const gsVector<T> & m_force;
    std::function < gsSparseMatrix<T> ( gsVector<T> const & ) > m_nonlinear;
    std::function < gsSparseMatrix<T> ( gsVector<T> const &, gsVector<T> const & ) > m_dnonlinear;
    const std::function < gsVector<T> ( gsVector<T> const & ) > m_residual;

    gsVector<T> m_solVec;

    mutable index_t m_verbose;
    bool m_NL;

    mutable gsOptionList m_options;

    mutable T m_toleranceU, m_toleranceF;
    mutable T m_relax;
    mutable index_t m_maxIterations;
    index_t m_iterations;

    bool m_converged;

    /// Linear solver employed
    gsSparseSolver<>::CGDiagonal m_solver;

    /// Indicator for bifurcation
    mutable T m_indicator;


    /// @brief Specifies (in)compressibility
    struct verbose
    {
        enum type
        {
            off = 0,
            iterations = 1,
            full = 2
        };
    };

};


} // namespace gismo


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsStaticSolver.hpp)
#endif