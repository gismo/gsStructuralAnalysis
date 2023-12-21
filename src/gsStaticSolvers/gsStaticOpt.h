/** @file gsStaticOpt.cpp

    @brief Static solver using optimization methods

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#include <gsStructuralAnalysis/src/gsStaticSolvers/gsStaticBase.h>

namespace gismo
{

/**
 * @brief Static solver using the Dynamic Relaxation method
 *
 * @tparam     T     coefficient type
 *
 * \ingroup gsStaticBase
 */
template <class T>
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
    gsStaticOpt( const Residual_t &Residual )
    :
    m_residualFun(Residual),
    m_ALresidualFun(nullptr),
    m_jacobian(nullptr)
    {
        this->_init();
    }

    /**
     * @brief      Constructor
     *
     * @param[in]  Residual  The residual
     * @param[in]  Jacobian  The jacobian
     */
    gsStaticOpt( const Residual_t &Residual,
                 const Jacobian_t &Jacobian )
    :
    m_residualFun(Residual),
    m_ALresidualFun(nullptr),
    m_jacobian(Jacobian)
    {
        this->_init();
    }


    /**
     * @brief      Constructs a new instance.
     *
     * @param[in]  ALResidual  The residual as arc-length object
     * @param[in]  Jacobian  The jacobian
     */
    gsStaticOpt( const ALResidual_t  & ALResidual )
    :
    m_ALresidualFun(ALResidual),
    m_jacobian(nullptr)
    {
        m_L = 1.0;
        m_residualFun = [this](gsVector<T> const & x, gsVector<T> & result) -> bool
        {
            return m_ALresidualFun(x,m_L,result);
        };
        this->_init();
    }

    /**
     * @brief      Constructs a new instance.
     *
     * @param[in]  ALResidual  The residual as arc-length object
     * @param[in]  Jacobian  The jacobian
     */
    gsStaticOpt( const ALResidual_t  & ALResidual,
                 const Jacobian_t &Jacobian )
    :
    m_ALresidualFun(ALResidual),
    m_jacobian(Jacobian)
    {
        m_L = 1.0;
        m_residualFun = [this](gsVector<T> const & x, gsVector<T> & result) -> bool
        {
            return m_ALresidualFun(x,m_L,result);
        };
        this->_init();
    }

/// gsStaticBase base functions
public:
    /// See \ref gsStaticBase
    gsStatus solve() override;

    /// See \ref gsStaticBase
    void initialize() override;

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

protected:
    /// See \ref solve()
    void _solve();
    /// Initializes the method
    void _init();

    gsVector<T> _computeResidual(const gsVector<T> & U);

public:

    /// Return the residual norm
    T residualNorm() const { return m_R.norm(); }

protected:
    Residual_t m_residualFun;
    const ALResidual_t m_ALresidualFun;
    const Jacobian_t m_jacobian;

    // Solution
    using Base::m_U;
    using Base::m_DeltaU;
    using Base::m_deltaU;

    using Base::m_L;
    using Base::m_DeltaL;
    using Base::m_deltaL;

    // Iterations
    using Base::m_numIterations;
    using Base::m_maxIterations;

    // Residuals
    using Base::m_R;

    // Tolerances
    using Base::m_tolF;
    using Base::m_tolU;

    // Residual norms
    using Base::m_residual;
    using Base::m_residualIni;
    using Base::m_residualOld;

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

//! [OptProblemExample Class]
template <typename T>
class gsOptProblemStatic : public gsOptProblem<T>
//! [OptProblemExample Class]
{
public:

    gsOptProblemStatic();

public:

    T evalObj( const gsAsConstVector<T> & u ) const;

    void gradObj_into( const gsAsConstVector<T> & u, gsAsVector<T> & result) const;

    void evalCon_into( const gsAsConstVector<T> & u, gsAsVector<T> & result) const;

    void jacobCon_into( const gsAsConstVector<T> & u, gsAsVector<T> & result) const;

private:

    // Lastly, we forward the memebers of the base clase gsOptProblem
    using gsOptProblem<T>::m_numDesignVars;
    using gsOptProblem<T>::m_numConstraints;
    using gsOptProblem<T>::m_numConJacNonZero;

    using gsOptProblem<T>::m_desLowerBounds;
    using gsOptProblem<T>::m_desUpperBounds;

    using gsOptProblem<T>::m_conLowerBounds;
    using gsOptProblem<T>::m_conUpperBounds;

    using gsOptProblem<T>::m_conJacRows;
    using gsOptProblem<T>::m_conJacCols;

    using gsOptProblem<T>::m_curDesign;
};
//! [OptProblem]

} //namespace


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsStaticOpt.hpp)
#endif
