 /** @file gsDynamicExplicitEuler.h

    @brief Class to perform time integration of second-order structural dynamics systems using the Explicit Euler method

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)

    TODO (June 2023):
    *    Change inputs to const references!
*/

#pragma once
#include <gsCore/gsLinearAlgebra.h>

#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicBase.h>
#include <gsIO/gsOptionList.h>

namespace gismo
{

/**
    @brief Performs the arc length method to solve a nonlinear system of equations.

    \tparam T coefficient type

    \ingroup gsStructuralAnalysis
*/
template <class T>
class gsDynamicExplicitEuler : public gsDynamicBase<T>
{
    typedef gsDynamicBase<T> Base;

protected:

    typedef typename gsStructuralAnalysisOps<T>::Force_t     Force_t;
    typedef typename gsStructuralAnalysisOps<T>::TForce_t    TForce_t;
    typedef typename gsStructuralAnalysisOps<T>::Residual_t  Residual_t;
    typedef typename gsStructuralAnalysisOps<T>::TResidual_t TResidual_t;
    typedef typename gsStructuralAnalysisOps<T>::Mass_t      Mass_t;
    typedef typename gsStructuralAnalysisOps<T>::TMass_t     TMass_t;
    typedef typename gsStructuralAnalysisOps<T>::Damping_t   Damping_t;
    typedef typename gsStructuralAnalysisOps<T>::TDamping_t  TDamping_t;
    typedef typename gsStructuralAnalysisOps<T>::Stiffness_t Stiffness_t;
    typedef typename gsStructuralAnalysisOps<T>::Jacobian_t  Jacobian_t;
    typedef typename gsStructuralAnalysisOps<T>::TJacobian_t TJacobian_t;

public:

    virtual ~gsDynamicExplicitEuler() {};

    /// Constructor
    gsDynamicExplicitEuler(
                    const Mass_t        & Mass,
                    const Damping_t     & Damping,
                    const Stiffness_t   & Stiffness,
                    const Force_t       & Force
                )
    : Base(Mass,Damping,Stiffness,Force)
    {}

    /// Constructor
    gsDynamicExplicitEuler(
                    const Mass_t        & Mass,
                    const Damping_t     & Damping,
                    const Stiffness_t   & Stiffness,
                    const TForce_t      & TForce
                )
    : Base(Mass,Damping,Stiffness,TForce)
    {}

    /// Constructor
    gsDynamicExplicitEuler(
                    const Mass_t        & Mass,
                    const Damping_t     & Damping,
                    const Jacobian_t    & Jacobian,
                    const Residual_t    & Residual
                )
    : Base(Mass,Damping,Jacobian,Residual)
    {}

    /// Constructor
    gsDynamicExplicitEuler(
                    const Mass_t        & Mass,
                    const Damping_t     & Damping,
                    const TJacobian_t   & TJacobian,
                    const TResidual_t   & TResidual
                )
    : Base(Mass,Damping,TJacobian,TResidual)
    {}

    /// Constructor
    gsDynamicExplicitEuler(
                    const TMass_t       & TMass,
                    const TDamping_t    & TDamping,
                    const TJacobian_t   & TJacobian,
                    const TResidual_t   & TResidual
                )
    : Base(TMass,TDamping,TJacobian,TResidual)
    {}

// General functions
public:

protected:

    gsStatus _step(const T dt);

    // static std::pair<gsvector,gsvector>

private:
    void _initOutput();
    void _stepOutput(const index_t it, const index_t resnorm, const index_t updatenorm);

    using Base::_computeResidual;
    using Base::_computeMass;
    using Base::_computeDamping;
    using Base::_computeJacobian;

protected:
    using Base::m_U;
    using Base::m_V;
    using Base::m_A;

    using Base::m_Uold;
    using Base::m_Vold;
    using Base::m_Aold;

    using Base::m_updateNorm;
    using Base::m_residualNorm;
    using Base::m_dt;
    using Base::m_time;

    using Base::m_tolerance;

    using Base::m_solver;

    using Base::m_maxIterations;
    using Base::m_converged;

};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsDynamicExplicitEuler.hpp)
#endif
