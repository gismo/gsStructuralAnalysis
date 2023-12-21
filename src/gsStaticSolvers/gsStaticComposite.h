 /** @file gsStaticComposite.h

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
class gsStaticComposite : public gsStaticBase<T>
{
protected:

    typedef gsStaticBase<T> Base;

    typedef typename Base::Residual_t    Residual_t;
    typedef typename Base::ALResidual_t  ALResidual_t;
    typedef typename Base::Jacobian_t    Jacobian_t;
    typedef typename Base::dJacobian_t   dJacobian_t;

public:

    /**
     * @brief      Constructs a new instance.
     *
     * @param[in]  solvers  The solvers
     */
    gsStaticComposite(  std::vector<gsStaticBase<T> * > solvers )
    :
    m_solvers(solvers)
    {
        this->defaultOptions();
    }

public:

    /// See \ref gsStaticBase
    gsStatus solve() override;

    /// See \ref gsStaticBase
    void initialize() override;

    /// See \ref gsStaticBase
    void defaultOptions() override;

    /// See \ref gsStaticBase
    void getOptions() override;

    /// See \ref gsStaticBase
    void setOptions(gsOptionList & options) override;

    /// See \ref gsStaticBase
    void setDisplacement(const gsVector<T> & displacement) override;

    /// See \ref gsStaticBase
    void setLoad(const T L) override;

    /// See \ref gsStaticBase
    void setSolution(const gsVector<T> & displacement, const T L) override;

    /// See \ref gsStaticBase
    void setUpdate(const gsVector<T> & update) override;

    /// See \ref gsStaticBase
    index_t numDofs() override { return m_solvers.front()->numDofs(); }

    /// See \ref gsStaticBase
    void reset() override;

protected:

    std::vector<gsStaticBase<T> *> m_solvers;

    using Base::m_U;
    using Base::m_DeltaU;

    using Base::m_L;

    using Base::m_verbose;

    using Base::m_options;

    using Base::m_numIterations;

    // Solver status
    using Base::m_status;
};


} // namespace gismo


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsStaticComposite.hpp)
#endif
