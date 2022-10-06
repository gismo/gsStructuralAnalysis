 /** @file gsALMRiks.h

    @brief Performs the arc length method to solve a nonlinear equation system.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#pragma once

#ifdef GISMO_WITH_SPECTRA
#include <gsSpectra/gsSpectra.h>
#endif

#include <gsStructuralAnalysis/gsALMBase.h>

namespace gismo
{

/**
    @brief Performs the Riks arc length method to solve a nonlinear equation system.

    \tparam T coefficient type

    \ingroup gsALMBase
*/
template <class T>
class gsALMRiks : public gsALMBase<T>
{

    typedef gsALMBase<T> Base;

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

public:

    /// Constructor
    gsALMRiks(  std::function < gsSparseMatrix<T> ( gsVector<T> const & ) > &Jacobian,
                std::function < gsVector<T> ( gsVector<T> const &, T, gsVector<T> const & ) > &Residual,
                gsVector<T> &Force )
    : Base(Jacobian,Residual,Force)
    {
        defaultOptions();
        getOptions();

        initMethods();
    }

    /// Constructor using the jacobian that takes the solution and the solution step
    gsALMRiks(  std::function < gsSparseMatrix<T> ( gsVector<T> const &, gsVector<T> const & ) > &dJacobian,
                std::function < gsVector<T> ( gsVector<T> const &, T, gsVector<T> const & ) > &Residual,
                gsVector<T> &Force )
    : Base(dJacobian,Residual,Force)
    {
        defaultOptions();
        getOptions();

        initMethods();
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
    /// Update of lambdaGeneralizedSelfAdjointEigenSolver
    using Base::m_DeltaL;
    /// Update of update of lambda
    using Base::m_deltaL;
    /// Vector with lambda updates
    using Base::m_deltaLs;

    /// Jacobian matrix
    using Base::m_jacMat;
    using Base::m_detKT;

    /// Scaling parameter
    T m_phi;
};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsALMRiks.hpp)
#endif
