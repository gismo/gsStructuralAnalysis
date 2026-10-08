 /** @file gsALMLoadControl.h

    @brief Performs the arc length method to solve a nonlinear equation system.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#pragma once

#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMBase.h>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsPathConstraint.h>

namespace gismo
{

/**
    @brief Performs the load-controlled arc length method to solve a nonlinear equation system.

    The path constraint is gsLoadConstraint: every step prescribes the load
    increment \f$\Delta\lambda=s\f$, the length \f$s\f$ of the step controller.
    The method cannot pass a limit point of the load.

    Options: those of gsALMBase; "ConstraintTol" is relative to the load
    increment.

    An initial guess (see setInitialGuess) is used for one load step to
    \a Lguess: the step is made with \f$\Delta\lambda=\lambda_{guess}-\lambda\f$
    and the guess is cleared once that step succeeds. The step length is that of
    the controller before the guess step and adapts to the iterations of the
    guess step as to those of a regular step.

    \tparam T coefficient type

    \ingroup gsALMBase
*/
template <class T>
class gsALMLoadControl : public gsALMBase<T>
{

    typedef gsALMBase<T> Base;

    typedef typename Base::ALResidual_t  ALResidual_t;
    typedef typename Base::Jacobian_t    Jacobian_t;
    typedef typename Base::dJacobian_t   dJacobian_t;

public:

    using Base::setLength;

protected:

    using Base::defaultOptions;
    using Base::getOptions;
    using Base::computeJacobian;
    using Base::computeResidual;
    using Base::computeResidualNorms;
    using Base::computeUt;
    using Base::computeUbar;
    using Base::computeStability;
    using Base::computeLength;

    /// The load constraint installed on the driver (same object as constraint())
    memory::shared_ptr<gsLoadConstraint<T> > m_load;

public:

    /// Constructor
    gsALMLoadControl(   Jacobian_t  &Jacobian,
                        ALResidual_t&ALResidual,
                        gsVector<T> &Force )
    : Base(Jacobian,ALResidual,Force), m_guessStep(false)
    {
        this->_install();
    }

    /// Constructor using the jacobian that takes the solution and the solution step
    gsALMLoadControl(   dJacobian_t &dJacobian,
                        ALResidual_t&ALResidual,
                        gsVector<T> &Force )
    : Base(dJacobian,ALResidual,Force), m_guessStep(false)
    {
        this->_install();
    }

    /// Constructor on a caller's problem (driver convention
    /// \f$R=\lambda F-F_{int}\f$, nothing is negated) and state space; both must outlive the solver
    gsALMLoadControl(   gsContinuationProblem<T> & problem,
                        gsStateSpace<T>          & space )
    : Base(problem,space), m_guessStep(false)
    {
        this->_install();
    }

    /// See gsALMBase. An initial guess makes one load step to \a Lguess; it is cleared once that step succeeds.
    virtual gsStatus step();

protected:

    /// Uses \a Lguess as the load increment of the next step: the step length
    /// of the controller is replaced by \f$\lambda_{guess}-\lambda\f$ until the step is over
    virtual void applyInitialGuess(const gsVector<T> & Uguess, T Lguess);

private:

    /// Installs the load constraint and reads the options
    void _install();

    /// True while the next step is made with the load increment of an initial guess
    bool m_guessStep;
    /// The controller before the guess step
    gsStepController<T> m_savedController;
};


} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsALMLoadControl.hpp)
#endif
