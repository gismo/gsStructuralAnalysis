/** @file gsStaticOpt.cpp

    @brief Static solver using optimization methods

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#include <gsStructuralAnalysis/src/gsStaticSolvers/gsStaticBase.h>
#include <gsOptimizer/gsGradientDescent.h>

#pragma once
namespace gismo
{

/**
 * @file gsStaticOpt.h
 * @brief Header file for the gsStaticOpt and gsOptProblemStatic classes.
 *
 * This file contains the definitions for the gsStaticOpt and gsOptProblemStatic classes,
 * which are used for static optimization in structural analysis.
 */

/**
 * @class gsOptProblemStatic
 * @brief A class representing a static optimization problem.
 *
 * @tparam T The coefficient type.
 *
 * This class inherits from gsOptProblem and provides functionality for evaluating
 * the objective function and its gradient for static optimization problems.
 */
template <typename T>
class gsOptProblemStatic : public gsOptProblem<T>
{
protected:

    typedef gsOptProblem<T> Base;
    typedef typename gsStructuralAnalysisOps<T>::Residual_t   Residual_t;
    typedef typename gsStructuralAnalysisOps<T>::ALResidual_t ALResidual_t;

public:
    /**
     * @brief Constructor for gsOptProblemStatic.
     *
     * @param residualFun The residual function.
     * @param L A coefficient (unused but needs to be assigned).
     * @param numDesignVars The number of design variables.
     */
    gsOptProblemStatic(const typename gsStructuralAnalysisOps<T>::Residual_t & residualFun, const T & L, index_t numDesignVars)
    :
    m_residualFun(residualFun),
    m_L(L) // unused, in fact, but needs to be assigned
    {
        m_numDesignVars = numDesignVars;
        m_curDesign.resize(numDesignVars,1);
        m_curDesign.setZero();
    }

    /**
     * @brief Constructor for gsOptProblemStatic with arc-length residual function.
     *
     * @param ALresidualFun The arc-length residual function.
     * @param L A coefficient.
     * @param numDesignVars The number of design variables.
     */
    gsOptProblemStatic(const typename gsStructuralAnalysisOps<T>::ALResidual_t & ALresidualFun, const T & L, index_t numDesignVars)
    :
    m_ALresidualFun(ALresidualFun),
    m_L(L)
    {
        m_numDesignVars = numDesignVars;
        m_residualFun = [this](gsVector<T> const & x, gsVector<T> & result) -> bool
        {
            return m_ALresidualFun(x,m_L,result);
        };
        m_curDesign.resize(numDesignVars,1);
        m_curDesign.setZero();
    }

public:
    /**
     * @brief Evaluates the objective function.
     *
     * @param u The input vector.
     * @return The norm of the residual.
     */
    T evalObj( const gsAsConstVector<T> & u ) const
    {
        gsVector<T> result;
        m_residualFun(u,result);
        return result.norm();
    }

    /**
     * @brief Computes the gradient of the objective function.
     *
     * @param u The input vector.
     * @param result The output vector for the gradient.
     */
    void gradObj_into( const gsAsConstVector<T> & u, gsAsVector<T> & result) const
    {
        gsVector<T> tmp;
        result.resize(u.rows());
        m_residualFun(u,tmp);
        result = -tmp;
    }

private:

    Residual_t m_residualFun;
    ALResidual_t m_ALresidualFun;
    const T & m_L;
    // Lastly, we forward the memebers of the base clase gsOptProblem
    using Base::m_numDesignVars;
    using Base::m_numConstraints;
    using Base::m_numConJacNonZero;

    using Base::m_desLowerBounds;
    using Base::m_desUpperBounds;

    using Base::m_conLowerBounds;
    using Base::m_conUpperBounds;

    using Base::m_conJacRows;
    using Base::m_conJacCols;

    using Base::m_curDesign;
};

/**
 * @class gsStaticOpt
 * @brief Static solver using the Dynamic Relaxation method.
 *
 * @tparam T The coefficient type.
 * @tparam Optimizer The optimizer type (default is gsGradientDescent<T>).
 *
 * This class inherits from gsStaticBase and provides functionality for solving
 * static optimization problems using the Dynamic Relaxation method.

 * \ingroup gsStaticBase
 */
template <class T, class Optimizer = gsGradientDescent<T>>
class gsStaticOpt : public gsStaticBase<T>
{
protected:

    typedef gsStaticBase<T> Base;

    typedef typename Base::Residual_t    Residual_t;
    typedef typename Base::ALResidual_t  ALResidual_t;

public:

    /**
     * @brief Constructor for gsStaticOpt.
     *
     * @param Residual The residual function.
     * @param numDofs The number of degrees of freedom.
     */
    gsStaticOpt( const Residual_t &Residual, const index_t numDofs )
    :
    m_optimizer(&m_optProblem),
    m_optProblem(Residual,m_L,numDofs)
    {
        m_dofs = numDofs;
        this->_init();
    }

    /**
     * @brief Constructor for gsStaticOpt with arc-length residual function.
     *
     * @param ALResidual The arc-length residual function.
     * @param numDofs The number of degrees of freedom.
     */
    gsStaticOpt( const ALResidual_t  & ALResidual, const index_t numDofs )
    :
    m_optimizer(&m_optProblem),
    m_optProblem(ALResidual,m_L,numDofs)
    {
        m_L = 1.0;
        m_dofs = numDofs;
        this->_init();
    }

public:
    /**
     * @brief Returns the optimizer options.
     *
     * @return A reference to the optimizer options.
     */
    gsOptionList & optimizerOptions() { return m_optimizer.options(); }

/// gsStaticBase base functions
public:
    /// See \ref gsStaticBase
    gsStatus solve() override;

    // /// See \ref gsStaticBase
    // void initialize() override;

    // /// See \ref gsStaticBase
    // void initOutput() override;

    // /// See \ref gsStaticBase
    // void stepOutput(index_t k) override;

    /// See \ref gsStaticBase
    void defaultOptions() override;

    // /// See \ref gsStaticBase
    // void reset() override;

    /// See \ref gsStaticBase
    void getOptions() override;

    /// Returns the stability indicator
    T indicator(const gsSparseMatrix<T> &, T)
    {
        GISMO_NO_IMPLEMENTATION;
    }

    /// Returns the stability vector
    gsVector<T> stabilityVec(const gsSparseMatrix<T> &, T)
    {
        GISMO_NO_IMPLEMENTATION;
    }

protected:
    /// See \ref solve()
    void _solve();
    /// Initializes the method
    void _init();


protected:
    Optimizer m_optimizer;
    gsOptProblemStatic<T> m_optProblem;

    // // Solution
    using Base::m_U;
    using Base::m_DeltaU;
    // using Base::m_deltaU;

    using Base::m_L;
    // using Base::m_DeltaL;
    // using Base::m_deltaL;

    // Iterations
    using Base::m_numIterations;
    using Base::m_maxIterations;

    // // Residuals
    // using Base::m_R;

    // // Tolerances
    using Base::m_tolF;
    using Base::m_tolU;

    // Options
    using Base::m_options;
    using Base::m_verbose;

    // DoFs
    using Base::m_dofs;

    // Headstart
    using Base::m_headstart;

    // Solver status
    using Base::m_status;
};



} //namespace

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsStaticOpt.hpp)
#endif
