 /** @file gsALMConsistentCrisfield.h

    @brief Performs the arc length method to solve a nonlinear equation system.

    Schweizerhof, K. H., & Wriggers, P. (1986).
    Consistent linearization for path following methods in nonlinear FE analysis.
    Computer Methods in Applied Mechanics and Engineering, 59(3), 261–279.
    https://doi.org/10.1016/0045-7825(86)90001-0

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#pragma once

#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMBase.h>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsPathConstraintCrisfield.h>

namespace gismo
{

/**
    @brief Performs the Consistent Crisfield arc length method to solve a nonlinear equation system.

    The path constraint is gsConsistentCrisfieldConstraint: the spherical
    constraint \f$\Delta U\cdot\Delta U+\varphi^2\|F\|^2\Delta\lambda^2=s^2\f$
    with the consistently linearised corrector of Schweizerhof and Wriggers.

    <b>Scaling.</b> With automatic scaling, \f$\varphi\f$ is one constant per run,
    \f$\varphi_0=\|u_{t,0}\|/\|F\|\f$ with \f$u_{t,0}=K^{-1}F\f$ at the start
    point of the run. The load term \f$\varphi\|F\|\Delta\lambda=\|u_{t,0}\|
    \Delta\lambda\f$ is then the displacement of the initial tangent for the
    load increment \f$\Delta\lambda\f$, so both terms of the constraint have the
    units of a displacement and the weight of the load term is the same at every
    step. A scaling that follows the path changes the metric of the constraint
    from step to step: the tangent ratio \f$\|u_t\|/\|F\|\f$ is unbounded at a
    limit point, where \f$K\f$ is singular, and the secant ratio
    \f$\|U\|/(|\lambda|\|F\|)\f$ is unbounded as \f$\lambda\to0\f$ at a regular
    point, which bounds \f$|\Delta\lambda|\le s\,|\lambda|/\|U\|\f$, so that
    \f$\lambda=0\f$ is approached only geometrically. A constant \f$\varphi\f$
    has neither singularity.

    A run starts at initialize(). Its \f$\varphi_0\f$ is taken from \f$u_t\f$ at
    the first step of the run, wherever that step starts (the origin, or a
    restart through setSolution), and is not changed by later steps, by
    applyOptions() or by putting the solver back at an earlier converged
    point. A run restarted exactly at a limit point, where \f$K\f$ is singular,
    freezes a huge scaling.

    Options, in addition to those of gsALMBase: "Scaling" (-1), the scaling
    \f$\varphi\f$; -1 selects the automatic scaling.

    An initial guess (see setInitialGuess) is applied as the previous increment
    \f$(U_{guess}-U,\ \lambda_{guess}-\lambda)\f$ of the secant predictor and is
    not cleared: it is applied again at every step while it differs from the
    present point.

    \tparam T coefficient type

    \ingroup gsALMBase
*/
template <class T>
class gsALMConsistentCrisfield : public gsALMBase<T>
{

    typedef gsALMBase<T> Base;

    typedef typename Base::ALResidual_t  ALResidual_t;
    typedef typename Base::Jacobian_t    Jacobian_t;
    typedef typename Base::dJacobian_t   dJacobian_t;

public:

    using Base::setLength;

protected:

    using Base::computeJacobian;
    using Base::computeResidual;
    using Base::computeResidualNorms;
    using Base::computeUt;
    using Base::computeUbar;
    using Base::computeStability;
    using Base::computeLength;

    /// The consistent Crisfield constraint installed on the driver (same object as constraint())
    memory::shared_ptr<gsConsistentCrisfieldConstraint<T> > m_cc;

public:

    /// Constructor
    gsALMConsistentCrisfield(   const Jacobian_t  &Jacobian,
                                const ALResidual_t&ALResidual,
                                const gsVector<T> &Force )
    : Base(Jacobian,ALResidual,Force)
    {
        this->_install();
    }

    /// Constructor using the jacobian that takes the solution and the solution step
    gsALMConsistentCrisfield(   const dJacobian_t &dJacobian,
                                const ALResidual_t&ALResidual,
                                const gsVector<T> &Force )
    : Base(dJacobian,ALResidual,Force)
    {
        this->_install();
    }

    /// Constructor on a caller's problem (driver convention
    /// \f$R=\lambda F-F_{int}\f$, nothing is negated) and state space; both must outlive the solver
    gsALMConsistentCrisfield(   gsContinuationProblem<T> & problem,
                                gsStateSpace<T>          & space )
    : Base(problem,space)
    {
        this->_install();
    }

public:

    /// The frozen scaling \f$\varphi_0=\|u_{t,0}\|/\|F\|\f$ of the run; 0 while none is frozen.
    T frozenScaling() const override { return m_cc->frozenScaling(); }

    /// True iff \f$\varphi_0\f$ is frozen (by the first step of the run or by setFrozenScaling()).
    bool hasFrozenScaling() const override { return m_cc->hasFrozenScaling(); }

    /// Freezes \f$\varphi_0\f$ to @a phi0 (> 0, finite; GISMO_ENSURE): the run uses it instead of
    /// freezing its own at its first step, so secant steps skip the tangent solve that freezing needs.
    /// Call it after initialize(), which resets \f$\varphi_0\f$; applyOptions() keeps it. Ignored
    /// while a user "Scaling" (not -1) is set. Used to give several solvers one arc metric.
    void setFrozenScaling(T phi0) override { m_cc->setFrozenScaling(phi0); }

protected:

    /// Forgets the frozen scaling: initialize() starts a new run
    virtual void initMethods();

    /// See gsALMBase; adds "Scaling"
    virtual void defaultOptions();
    /// See gsALMBase; maps "Scaling" to the constraint without touching the frozen scaling
    virtual void getOptions();

private:

    /// Installs the constraint and reads the options
    void _install();
};


} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsALMConsistentCrisfield.hpp)
#endif
