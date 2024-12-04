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

template <typename T> 
class gsOptProblemStatic : public gsOptProblem<T>
{
protected:

    typedef gsOptProblem<T> Base;
    typedef typename gsStructuralAnalysisOps<T>::Residual_t   Residual_t;
    typedef typename gsStructuralAnalysisOps<T>::ALResidual_t ALResidual_t;

public:

    gsOptProblemStatic(const typename gsStructuralAnalysisOps<T>::Residual_t & residualFun, const T & L, index_t numDesignVars)
    :
    m_residualFun(residualFun),
    m_L(L) // unused, in fact, but needs to be assigned
    {
        m_numDesignVars = numDesignVars;
        m_curDesign.resize(numDesignVars,1);
        m_curDesign.setZero();
    }

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

    T evalObj( const gsAsConstVector<T> & u ) const
    {
        gsVector<T> result;
        m_residualFun(u,result);
        return result.norm();
    }

    void gradObj_into( const gsAsConstVector<T> & u, gsAsVector<T> & result) const
    {
        gsVector<T> tmp;
        result.resize(u.rows());
        m_residualFun(u,tmp);
        result = -tmp;
    }

    // void evalCon_into( const gsAsConstVector<T> & u, gsAsVector<T> & result) const;

    // void jacobCon_into( const gsAsConstVector<T> & u, gsAsVector<T> & result) const;

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
 * @brief Static solver using the Dynamic Relaxation method
 *
 * @tparam     T     coefficient type
 *
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
     * @brief      Constructor
     *
     * @param[in]  Residual  The residual
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
     * @brief      Constructs a new instance.
     *
     * @param[in]  ALResidual  The residual as arc-length object
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

// Class-specific functions
// public:
//     void initialize();

protected:
    /// See \ref solve()
    void _solve();
    /// Initializes the method
    void _init();

//     gsVector<T> _computeResidual(const gsVector<T> & U);

// public:

//     /// Return the residual norm
//     T residualNorm() const { return m_R.norm(); }

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
