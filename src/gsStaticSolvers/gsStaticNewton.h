 /** @file gsStaticNewton.h

    @brief Static solver using a newton method

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#include <typeinfo>

#include <gsStructuralAnalysis/src/gsStaticSolvers/gsStaticBase.h>
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

    typedef typename Base::Residual_t    Residual_t;
    typedef typename Base::ALResidual_t  ALResidual_t;
    typedef typename Base::Jacobian_t    Jacobian_t;
    typedef typename Base::dJacobian_t   dJacobian_t;

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
        m_dnonlinear = [this](gsVector<T> const & x, gsVector<T> const & dx, gsSparseMatrix<T> & m) -> bool
        {
            return m_nonlinear(x,m);
        };

        this->_init();
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
        m_residualFun = [this](gsVector<T> const & x, gsVector<T> & result) -> bool
        {
            return m_ALresidualFun(x,m_L,result);
        };

        m_dnonlinear = [this](gsVector<T> const & x, gsVector<T> const & dx, gsSparseMatrix<T> & m) -> bool
        {
            return m_nonlinear(x,m);
        };

        this->_init();
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
    }

public:

    /// See \ref gsStaticBase
    gsStatus solve() override
    {
        if (m_NL)
            return solveNonlinear();
        else
            return solveLinear();
    }
    /// Perform a linear solve
    gsStatus solveLinear();
    gsStatus solveLinear(gsVector<T> & solution)
    {
        gsStatus status = this->solveLinear();
        solution = m_U;
        return status;
    }
    /// Perform a nonlinearg solve
    gsStatus solveNonlinear(gsVector<T> & solution)
    {
        gsStatus status = this->solveNonlinear();
        solution = m_U;
        return status;
    }
    gsStatus solveNonlinear();

    /// See \ref gsStaticBase
    void initOutput() override;

    /// See \ref gsStaticBase
    void stepOutput(index_t k) override;

    /// See \ref gsStaticBase
    void defaultOptions() override;

    /// See \ref gsStaticBase
    void reset() override;

    /// See \ref gsStaticBase
    void getOptions() override;

    using Base::indicator;
    using Base::stabilityVec;

    /// See \ref gsStaticBase
    T indicator()
    {
        gsSparseMatrix<T> m;
        GISMO_ENSURE(m_nonlinear(m_U, m),"Assembly failed");
        return indicator(m);
    }

    /// See \ref gsStaticBase
    gsVector<T> stabilityVec()
    {
        gsSparseMatrix<T> m;
        GISMO_ENSURE(m_nonlinear(m_U, m),"Assembly failed");
        return stabilityVec(m);
    }

protected:

    /// Perform a linear solve
    gsVector<T> _solveLinear();
    /// Perform a nonlinear solve
    gsVector<T> _solveNonlinear();

    using Base::_computeStability;

    /// Initializes the method
    void _init();
    /// Starts the method
    void _start();
    /// Factorizes the \a jacMat
    void _factorizeMatrix(const gsSparseMatrix<T> & jacMat) const;
    /// Solves the system with RHS \a F and LHS the Jacobian
    gsVector<T> _solveSystem(const gsVector<T> & F);

    gsVector<T> _computeResidual(const gsVector<T> & U);

    gsSparseMatrix<T> _computeJacobian(const gsVector<T> & U, const gsVector<T> & deltaU);

protected:

    const gsSparseMatrix<T> & m_linear;
    const gsVector<T> & m_force;
    const Jacobian_t      m_nonlinear;
          dJacobian_t     m_dnonlinear;
          Residual_t      m_residualFun;
    const ALResidual_t    m_ALresidualFun;

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

    // Residual norms
    using Base::m_residual;
    using Base::m_residualIni;
    using Base::m_residualOld;

    using Base::m_start;
    using Base::m_headstart;

    /// Linear solver employed
    using Base::m_solver;   // Cholesky by default
    
    using Base::m_stabilityMethod;

    /// Indicator for bifurcation
    using Base::m_indicator;

    using Base::m_dofs;

    // Solver status
    using Base::m_status;
};


} // namespace gismo


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsStaticNewton.hpp)
#endif
