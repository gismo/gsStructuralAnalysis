 /** @file gsALMRiks.h

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
    @brief Performs the Riks arc length method to solve a nonlinear equation system.

    The path constraint is gsRiksConstraint: the quadratic constraint
    \f$\varphi\|\Delta U\|^2+(1-\varphi)\Delta\lambda^2=s^2\f$ with the
    constant weight \f$\varphi=1/n\f$ (\f$n\f$ the number of unknowns), solved
    by Newton's method.

    Options: those of gsALMBase.

    An initial guess (see setInitialGuess) replaces the previous increment of
    the secant predictor by \f$(U_{guess}-U,\ \lambda_{guess}-\lambda)\f$ and is
    cleared by the step that uses it.

    \tparam T coefficient type

    \ingroup gsALMBase
*/
template <class T>
class gsALMRiks : public gsALMBase<T>
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

    /// The Riks constraint installed on the driver (same object as constraint())
    memory::shared_ptr<gsRiksConstraint<T> > m_riks;

public:

    /// Constructor
    gsALMRiks(  const Jacobian_t  &Jacobian,
                const ALResidual_t&ALResidual,
                const gsVector<T> &Force )
    : Base(Jacobian,ALResidual,Force)
    {
        this->_install();
    }

    /// Constructor using the jacobian that takes the solution and the solution step
    gsALMRiks(  const dJacobian_t &dJacobian,
                const ALResidual_t&ALResidual,
                const gsVector<T> &Force )
    : Base(dJacobian,ALResidual,Force)
    {
        this->_install();
    }

    /// Constructor on a caller's problem (driver convention
    /// \f$R=\lambda F-F_{int}\f$, nothing is negated) and state space; both must outlive the solver
    gsALMRiks(  gsContinuationProblem<T> & problem,
                gsStateSpace<T>          & space )
    : Base(problem,space)
    {
        this->_install();
    }

protected:

    /// See gsALMBase; the guess is cleared once it is used
    virtual void applyInitialGuess(const gsVector<T> & Uguess, T Lguess);

private:

    /// Installs the Riks constraint and reads the options
    void _install();
};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsALMRiks.hpp)
#endif
