 /** @file gsStaticNewton.h

    @brief Performs linear modal analysis given a matrix or functions of a matrix

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#include <typeinfo>
#include <gsSpectra/gsSpectra.h>
#include <gsStructuralAnalysis/gsStaticBase.h>
#pragma once


namespace gismo
{

/**
    @brief Performs the arc length method to solve a nonlinear equation system.

    \tparam T coefficient type

    \ingroup ThinShell
*/
template <class T>
class gsStaticNewton : public gsStaticBase<T>
{
protected:

    typedef gsStaticBase<T> Base;

    typedef std::function<gsVector<T>(gsVector<T> const &)   >          Residual_t;
    typedef std::function<gsVector<T>(gsVector<T> const &, T)>          ALResidual_t;
    typedef std::function < gsSparseMatrix<T> ( gsVector<T> const & ) > Jacobian_t;

public:

  /// Constructor giving access to the gsShellAssembler object to create a linear system per iteration
  gsStaticNewton(   const gsSparseMatrix<T> &linear,
                    const gsVector<T> &force
                    ) :
    m_linear(linear),
    m_force(force),
    m_nonlinear(nullptr),
    m_residualFun(nullptr)
  {
    this->_init();
    m_NL = false;
    m_U.setZero(m_dofs);
  }

  gsStaticNewton(   const gsSparseMatrix<T> &linear,
                    const gsVector<T> &force,
                    const Jacobian_t &nonlinear,
                    const Residual_t &residual
                    ) :
    m_linear(linear),
    m_force(force),
    m_nonlinear(nonlinear),
    m_residualFun(residual),
    m_ALresidualFun(nullptr)
  {
    this->_init();
    m_NL = true;
  }

  gsStaticNewton(   const gsSparseMatrix<T> &linear,
                    const gsVector<T>  &force,
                    const Jacobian_t   &nonlinear,
                    const ALResidual_t &ALResidual
                    ) :
    m_linear(linear),
    m_force(force),
    m_nonlinear(nonlinear),
    m_ALresidualFun(ALResidual)
  {
    m_L = 1.0;
    m_residualFun = [this](gsVector<T> const & x)
    {
        return m_ALresidualFun(x,m_L);
    };
    this->_init();
    m_NL = true;
  }

public:

    void solve() { m_U = solveNonlinear(); }
    gsVector<T> solveLinear();
    gsVector<T> solveNonlinear();

    void initialize();

    void initOutput();
    void stepOutput(index_t k);

    void defaultOptions();
    void getOptions();


    using Base::indicator;
    using Base::stabilityVec;

    T indicator()
    {
        return indicator(m_nonlinear(m_U));
    }

    gsVector<T> stabilityVec()
    {
        return stabilityVec(m_nonlinear(m_U));
    }

protected:

    using Base::_computeStability;


    void _init();
    void _start();
    void _factorizeMatrix(const gsSparseMatrix<T> & jacMat) const;
    gsVector<T> _solveSystem(const gsVector<T> & F);

protected:

    const gsSparseMatrix<T> & m_linear;
    const gsVector<T> & m_force;
    const Jacobian_t m_nonlinear;
    Residual_t m_residualFun;
    const ALResidual_t m_ALresidualFun;

    using Base::m_R;

    using Base::m_U;
    using Base::m_DeltaU;
    using Base::m_deltaU;

    using Base::m_L;

    using Base::m_stabilityVec;

    using Base::m_verbose;

    bool m_NL;

    using Base::m_options;

    using Base::m_tolF;
    using Base::m_tolU;

    mutable T m_relax;

    using Base::m_maxIterations;
    using Base::m_numIterations;
    using Base::m_converged;

    // Residual norms
    using Base::m_residual;
    using Base::m_residualIni;
    using Base::m_residualOld;

    using Base::m_start;
    using Base::m_headstart;

    /// Linear solver employed
    using Base::m_LDLTsolver;   // Cholesky
    mutable gsSparseSolver<>::CGDiagonal      m_CGsolver;     // CG

    using Base::m_stabilityMethod;
    mutable index_t m_solverType;

    /// Indicator for bifurcation
    using Base::m_indicator;

    using Base::m_dofs;

    struct solver
    {
        enum type
        {
            LDLT = 0,
            CG  = 1, // The CG solver is robust for membrane models, where zero-blocks in the matrix might occur.
        };
    };



};


} // namespace gismo


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsStaticNewton.hpp)
#endif