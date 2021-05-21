 /** @file gsStaticSolver.h

    @brief Performs linear modal analysis given a matrix or functions of a matrix

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#include <typeinfo>
#include <gsSpectra/gsSpectra.h>
#pragma once


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
   m_bifurcationMethod = 0;
   m_solverType = 0;
   m_solVec.setZero(force.rows());
   defaultOptions();
  }

  gsStaticSolver(   const gsSparseMatrix<T> &linear,
                    const gsVector<T> &force,
                    const std::function < gsSparseMatrix<T> ( gsVector<T> const & ) > &nonlinear,
                    const std::function < gsVector<T> ( gsVector<T> const & ) > &residual
                    ) :
    m_linear(linear),
    m_force(force),
    m_nonlinear(nonlinear),
    m_residual(residual)
  {
    m_NL = true;
    m_bifurcationMethod = 0;
    m_solverType = 0;
    m_converged = false;
    m_headstart = false;
    // m_solVec.setZero(force.rows());
    defaultOptions();
  }

public:

    gsVector<T> solveLinear();
    gsVector<T> solveNonlinear();

    gsOptionList & options() {return m_options;}
    void setOptions(gsOptionList opt) { m_options.update(opt,gsOptionList::addIfUnknown); } // gsOptionList opt
    void defaultOptions();
    void getOptions() const;

    void factorizeMatrix(const gsSparseMatrix<T> & M) const;
    gsVector<T> solveSystem(const gsVector<T> & F);

    int iterations() {return m_iterations;}

    gsVector<T> solution()  { return m_solVec; }
    gsVector<T> increment() { return m_DeltaU; }

    bool converged() { return m_converged; }

    T indicator() const
    {
        _computeStability(m_solVec);
        return m_indicator;
    }

    void setSolution(const gsVector<T> & solution)
    {
        m_solVec = solution;
    }

    void setSolutionStep(const gsVector<T> & step)
    {
        m_DeltaU = step;
        m_headstart = true;
    }

    void setSolutionStep(const gsVector<T> & step, const gsVector<T> & solution)
    {
        m_DeltaU = step;
        m_solVec = solution;
        m_headstart = true;
    }


protected:
    void _computeStability(const gsVector<T> x) const;

protected:

    const gsSparseMatrix<T> & m_linear;
    const gsVector<T> & m_force;
    const std::function < gsSparseMatrix<T> ( gsVector<T> const & ) > m_nonlinear;
    const std::function < gsVector<T> ( gsVector<T> const & ) > m_residual;

    mutable gsVector<T> m_solVec;
    mutable gsVector<T> m_DeltaU, m_deltaU;

    mutable index_t m_verbose;
    bool m_NL;

    mutable gsOptionList m_options;

    mutable T m_toleranceU, m_toleranceF;
    mutable T m_relax;
    mutable index_t m_maxIterations;
    index_t m_iterations;

    bool m_converged;
    bool m_headstart;

    /// Linear solver employed
    mutable gsSparseSolver<>::SimplicialLDLT  m_LDLTsolver;   // Cholesky
    mutable gsSparseSolver<>::CGDiagonal      m_CGsolver;     // CG

    mutable index_t m_bifurcationMethod;
    mutable index_t m_solverType;

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

    struct solver
    {
        enum type
        {
            LDLT = 0,
            CG  = 1, // The CG solver is robust for membrane models, where zero-blocks in the matrix might occur.
        };
    };

    struct bifmethod
    {
        enum type
        {
            Determinant = 0,
            Eigenvalue  = 1,
        };
    };

};


} // namespace gismo


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsStaticSolver.hpp)
#endif