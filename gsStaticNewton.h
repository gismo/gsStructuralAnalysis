 /** @file gsStaticNewton.h

    @brief Static solver using a newton method

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
    @brief Static solver using a newton method

    \tparam T coefficient type

    \ingroup gsStructuralAnalysis
*/
template <class T>
class gsStaticNewton : public gsStaticBase<T>
{
protected:

    typedef gsStaticBase<T> Base;

    typedef std::function<gsVector<T>(gsVector<T> const &)   >                                          Residual_t;
    typedef std::function<gsVector<T>(gsVector<T> const &, T)>                                          ALResidual_t;
    typedef std::function < gsSparseMatrix<T> ( gsVector<T> const & ) >                                 Jacobian_t;
    typedef std::function<gsSparseMatrix<real_t> (gsVector<real_t> const &,gsVector<real_t> const &)>   dJacobian_t;

public:

    /**
     * @brief      Constructor
     *
     * @param[in]  linear  The linear stiffness matrix
     * @param[in]  force   The external force
     */
    gsStaticNewton(     const gsSparseMatrix<T> &linear,
                        const gsVector<T> &force        )
    :
        m_linear(linear),
        m_force(force),
        m_nonlinear(nullptr),
        m_residualFun(nullptr)
    {
        this->_init();
        m_NL = false;
        m_U.setZero(m_dofs);
    }

    /**
     * @brief      Constructor
     *
     * @param[in]  linear     The linear stiffness matrix
     * @param[in]  force      The external force
     * @param[in]  nonlinear  The Jacobian
     * @param[in]  residual   The residual
     */
    gsStaticNewton( const gsSparseMatrix<T> &linear,
                    const gsVector<T> &force,
                    const Jacobian_t &nonlinear,
                    const Residual_t &residual      )
    :
        m_linear(linear),
        m_force(force),
        m_nonlinear(nonlinear),
        m_dnonlinear(nullptr),
        m_residualFun(residual),
        m_ALresidualFun(nullptr)
    {
        m_dnonlinear = [this](gsVector<T> const & x, gsVector<T> const & dx)
        {
            return m_nonlinear(x);
        };

        this->_init();
        m_NL = true;
    }

    /**
     * @brief      Constructor
     *
     * @param[in]  linear     The linear stiffness matrix
     * @param[in]  force      The external force
     * @param[in]  nonlinear  The Jacobian
     * @param[in]  residual   The residual as arc-length object
     */
    gsStaticNewton( const gsSparseMatrix<T> &linear,
                    const gsVector<T>  &force,
                    const Jacobian_t   &nonlinear,
                    const ALResidual_t &ALResidual  )
    :
        m_linear(linear),
        m_force(force),
        m_nonlinear(nonlinear),
        m_dnonlinear(nullptr),
        m_residualFun(nullptr),
        m_ALresidualFun(ALResidual)
    {
        m_L = 1.0;
        m_residualFun = [this](gsVector<T> const & x)
        {
            return m_ALresidualFun(x,m_L);
        };

        m_dnonlinear = [this](gsVector<T> const & x, gsVector<T> const & dx)
        {
            return m_nonlinear(x);
        };

        this->_init();
        m_NL = true;
    }

    /**
     * @brief      { function_description }
     *
     * @param[in]  linear      The linear matrix
     * @param[in]  force       The force vector
     * @param[in]  dnonlinear  The jacobian taking the solution x and the iterative update dx
     * @param[in]  residual    The residual function
     */
    gsStaticNewton(   const gsSparseMatrix<T> &linear,
                      const gsVector<T> &force,
                      const dJacobian_t &dnonlinear,
                      const Residual_t &residual
                      )
    :
        m_linear(linear),
        m_force(force),
        m_nonlinear(nullptr),
        m_dnonlinear(dnonlinear),
        m_residualFun(residual),
        m_ALresidualFun(nullptr)
    {
        this->_init();
        m_NL = true;
    }

public:

    /// See \ref gsStaticBase
    void solve() { m_U = solveNonlinear(); }
    /// Perform a linear solve
    gsVector<T> solveLinear();
    /// Perform a nonlinear solve
    gsVector<T> solveNonlinear();

    /// See \ref gsStaticBase
    void initialize() {_init(); };

    /// See \ref gsStaticBase
    void initOutput();
    /// See \ref gsStaticBase
    void stepOutput(index_t k);

    /// See \ref gsStaticBase
    void defaultOptions();
    /// See \ref gsStaticBase
    void getOptions();

    using Base::indicator;
    using Base::stabilityVec;

    /// See \ref gsStaticBase
    T indicator()
    {
        return indicator(m_nonlinear(m_U));
    }

    /// See \ref gsStaticBase
    gsVector<T> stabilityVec()
    {
        return stabilityVec(m_nonlinear(m_U));
    }

protected:

    using Base::_computeStability;

    /// Initializes the method
    void _init();
    /// Starts the method
    void _start();
    /// Factorizes the \a jacMat
    void _factorizeMatrix(const gsSparseMatrix<T> & jacMat) const;
    /// Solves the system with RHS \a F and LHS the Jacobian
    gsVector<T> _solveSystem(const gsVector<T> & F);

protected:

    const gsSparseMatrix<T> & m_linear;
    const gsVector<T> & m_force;
    const Jacobian_t m_nonlinear;
    dJacobian_t m_dnonlinear;
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