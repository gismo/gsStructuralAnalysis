 /** @file gsALMCrisfield.h

    @brief Performs the arc length method to solve a nonlinear equation system.

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
    @brief Performs the Crisfield arc length method to solve a nonlinear equation system.

    The path constraint is gsCrisfieldConstraint: the spherical constraint
    \f$\Delta U\cdot\Delta U+\varphi^2\|F\|^2\Delta\lambda^2=s^2\f$ with the
    quadratic solved for the load increment (see there for the root choice and
    the fallbacks).

    Options, in addition to those of gsALMBase:
    <ul>
    <li>"Scaling" (-1): the scaling \f$\varphi\f$; -1 selects the automatic scaling.</li>
    <li>"AngleMethod" (0): root choice by the previous increment of 0 = the last
    accepted step, 1 = the previous iteration, 2 = the predictor. Any other value
    is rejected by a GISMO_ENSURE when the options are applied.</li>
    <li>"Relaxation" (1): damping of an oscillating load factor, see gsALMBase.</li>
    </ul>

    <b>Automatic scaling.</b> With \c Scaling = -1 the scaling is recomputed at
    every step from the secant ratio \f$\|U\|/(|\lambda|\|F\|)\f$; it is
    unbounded as \f$\lambda\to0\f$ and the method can stall near limit points.
    gsALMConsistentCrisfield freezes one \f$\varphi_0\f$ per run instead.

    An initial guess (see setInitialGuess) makes the previous increment of the
    predictor \f$-(U_{guess}-U,\ \lambda_{guess}-\lambda)\f$, with the tangent
    ratio \f$\|u_t\|/\|F\|\f$ as the automatic scaling of that step; it is
    cleared by the step that uses it.

    \tparam T coefficient type

    \ingroup gsALMBase
*/
template <class T>
class gsALMCrisfield : public gsALMBase<T>
{

    typedef gsALMBase<T> Base;

    typedef typename Base::ALResidual_t  ALResidual_t;
    typedef typename Base::Jacobian_t    Jacobian_t;
    typedef typename Base::dJacobian_t   dJacobian_t;

public:

    using Base::setLength;
    using Base::computeStability;

protected:

    using Base::computeJacobian;
    using Base::computeResidual;
    using Base::computeResidualNorms;
    using Base::computeUt;
    using Base::computeUbar;
    using Base::computeLength;

    /// The Crisfield constraint installed on the driver (same object as constraint())
    memory::shared_ptr<gsCrisfieldConstraint<T> > m_cris;

public:

    /// Constructor
    gsALMCrisfield( const Jacobian_t  &Jacobian,
                    const ALResidual_t&ALResidual,
                    const gsVector<T> &Force )
    : Base(Jacobian,ALResidual,Force)
    {
        this->_install();
    }

    /// Constructor using the jacobian that takes the solution and the solution step
    gsALMCrisfield( const dJacobian_t &dJacobian,
                    const ALResidual_t&ALResidual,
                    const gsVector<T> &Force )
    : Base(dJacobian,ALResidual,Force)
    {
        this->_install();
    }

    /// Constructor on a caller's problem (driver convention
    /// \f$R=\lambda F-F_{int}\f$, nothing is negated) and state space; both must outlive the solver
    gsALMCrisfield( gsContinuationProblem<T> & problem,
                    gsStateSpace<T>          & space )
    : Base(problem,space)
    {
        this->_install();
    }

protected:

    /// Resets the counters of the root choices of the constraint
    virtual void initMethods();

    /// See gsALMBase; adds "Scaling" and "AngleMethod"
    virtual void defaultOptions();
    /// See gsALMBase; maps "Scaling", "AngleMethod" and "Relaxation" to the constraint
    virtual void getOptions();

    /// See gsALMBase; the guess is cleared once it is used
    virtual void applyInitialGuess(const gsVector<T> & Uguess, T Lguess);

protected:
    /// Angle determination method option
    struct angmethod
    {
        enum type
        {
            Step = 0,
            Iteration  = 1,
            Predictor  = 2,
        };
    };

private:

    /// Installs the Crisfield constraint and reads the options
    void _install();
};


} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsALMCrisfield.hpp)
#endif
