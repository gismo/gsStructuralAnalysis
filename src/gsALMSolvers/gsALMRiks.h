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

namespace gismo
{

/**
    @brief Performs the Riks arc length method to solve a nonlinear equation system.

    \tparam T coefficient type

    \ingroup gsALMSolvers
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
    // computeStability() is PUBLIC on gsALMBase; re-exporting it protected here would
    // narrow the base interface (callers holding a gsALMBase& could reach it, callers
    // holding the derived type could not). Kept public, as in gsALMCrisfield.
    using Base::computeStability;

protected:

    using Base::defaultOptions;
    using Base::getOptions;
    using Base::computeJacobian;
    using Base::computeResidual;
    using Base::computeResidualNorms;
    using Base::computeUt;
    using Base::computeUbar;
    using Base::computeLength;

public:

    /// Constructor
    gsALMRiks(  const Jacobian_t  &Jacobian,
                const ALResidual_t&ALResidual,
                const gsVector<T> &Force )
    : Base(Jacobian,ALResidual,Force)
    {
        defaultOptions();
        getOptions();

        initMethods();
    }

    /// Constructor using the jacobian that takes the solution and the solution step
    gsALMRiks(  const dJacobian_t &dJacobian,
                const ALResidual_t&ALResidual,
                const gsVector<T> &Force )
    : Base(dJacobian,ALResidual,Force)
    {
        defaultOptions();
        getOptions();

        initMethods();
    }

public:
    /// Distance in the (U,L) plane, measured in the CONVEX constraint metric this class
    /// enforces: \f$ \phi|\Delta U|^2 + (1-\phi)\Delta\Lambda^2 = \Delta s^2 \f$, with
    /// \f$\phi\f$ = \a m_convexWeight (see its declaration; it is NOT the metric scaling
    /// \f$\psi\f$ that gsALMCrisfield/gsALMConsistentCrisfield expose as their \c Scaling
    /// option).
    /// @note PUBLIC, and the value ESCAPES this class - it is load-bearing, not diagnostic.
    ///       Verified at source 2026-08-03, CORRECTED 2026-08-04: the live readers are
    ///       \a gsAPALM::_initiation and \a gsAPALM::_correction, which accumulate it into their
    ///       interval distances (including the \c upperDistance / \c lowerDistance of a bisected
    ///       interval). Outside those, the only in-tree readers are stepOutput() and the unit
    ///       tests.
    ///       !! The earlier version of this note also credited \c gsAPALMBase::parallelSolve with
    ///       two more consumer sites. That class was DEAD CODE and was DELETED 2026-08-04, so the
    ///       live consumer count is FIVE, not the seven recorded earlier. Anyone re-deriving
    ///       the read-site map should start from five.
    T distance(const gsVector<T>& DeltaU, const T DeltaL) const
    {
        return math::pow(m_convexWeight * math::pow(DeltaU.norm(),2.0) + (1.0-m_convexWeight) * math::pow(DeltaL,2.0),0.5);
    }

protected:

    /// See gsALMBase
    void initMethods();
    /// See gsALMBase
    void initiateStep();
    /// See gsALMBase
    void iterationFinish();

    /// See gsALMBase
    void quasiNewtonPredictor();
    /// See gsALMBase
    void quasiNewtonIteration();

    /// See gsALMBase
    void predictor();
    void predictorGuess();
    /// See gsALMBase
    void iteration();

    /// See gsALMBase
    void initOutput();
    /// See gsALMBase
    void stepOutput();

protected:

    // Number of degrees of freedom
    using Base::m_numDof;

    using Base::m_jacobian;
    using Base::m_djacobian;
    using Base::m_residualFun;
    using Base::m_forcing;

    /// Solver options
    using Base::m_options;

    /// Number of Arc Length iterations performed
    using Base::m_numIterations;

    /// Maximum number of Arc Length iterations allowed
    using Base::m_maxIterations;

    /// Length of the step in the u,f plane
    using Base::m_arcLength;
    using Base::m_arcLength_prev;

    /// Output verbosity
    using Base::m_verbose;

    /// Note
    using Base::m_note;

    /// Convergence result
    using Base::m_converged;

    /// Force residuum
    using Base::m_residueF;

    /// Displacement residuum
    using Base::m_residueU;

    /// Load residuum
    using Base::m_residueL;

    /// Indicator for bifurcation
    using Base::m_indicator;
    using Base::m_negatives;

    /// Relaxation factor
    using Base::m_relax;

    // Previous update
    using Base::m_DeltaUold;
    using Base::m_DeltaLold;
    /// Displacement vector (present, at previously converged point)
    using Base::m_U;
    using Base::m_Uprev;
    using Base::m_Uguess;
    /// Update of displacement vector
    using Base::m_DeltaU;
    /// u_bar
    using Base::m_deltaUbar;
    /// u_t
    using Base::m_deltaUt;
    /// Update of update of displacement vector
    using Base::m_deltaU;

    /// Lambda (present, at previously converged point)
    using Base::m_L;
    using Base::m_Lprev;
    using Base::m_Lguess;
    /// Update of lambdaGeneralizedSelfAdjointEigenSolver
    using Base::m_DeltaL;
    /// Update of update of lambda
    using Base::m_deltaL;
    /// Vector with lambda updates
    using Base::m_deltaLs;

    /// Jacobian matrix
    using Base::m_jacMat;
    using Base::m_detKT;

    /// Convex interpolation weight \f$\phi \in (0,1)\f$ of Riks's constraint
    /// \f$ \phi|\Delta U|^2 + (1-\phi)\Delta\Lambda^2 = \Delta s^2 \f$. Purely INTERNAL:
    /// this class registers no \c Scaling option and the weight is hard-set to
    /// \f$1/n_{dof}\f$ by \a predictor() / \a predictorGuess().
    ///
    /// @note Renamed out of the historical member name \c phi, which the siblings
    ///       gsALMCrisfield and gsALMConsistentCrisfield still carry for an INCOMPATIBLE
    ///       quantity - the metric scaling \f$\psi\f$ of \f$A_0 = \psi^2\|f\|^2\f$, exposed
    ///       as their \c Scaling option. The two are related by
    ///       \f$\psi^2 = (1-\phi)/\phi\f$; only the test harness translates between them.
    T m_convexWeight;
};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsALMRiks.hpp)
#endif
