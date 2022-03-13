 /** @file gsALMCrisfield.h

    @brief Performs the arc length method to solve a nonlinear equation system.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#pragma once
#include <gsSpectra/gsSpectra.h>
#include <gsStructuralAnalysis/gsALMBase.h>

namespace gismo
{

/**
    @brief Performs the arc length method to solve a nonlinear equation system.

    \tparam T coefficient type

    \ingroup gsStructuralAnalysis
*/
template <class T>
class gsALMCrisfield : public gsALMBase<T>
{

    typedef gsALMBase<T> Base;

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

public:

    /// Constructor giving access to the gsShellAssembler object to create a linear system per iteration
    gsALMCrisfield( std::function < gsSparseMatrix<T> ( gsVector<T> const & ) > &Jacobian,
                    std::function < gsVector<T> ( gsVector<T> const &, T, gsVector<T> const & ) > &Residual,
                    gsVector<T> &Force )
    : Base(Jacobian,Residual,Force)
    {
        defaultOptions();
        getOptions();

        initMethods();
    }

    T distance(gsVector<T>& DeltaU, T DeltaL)
    {
        T A0 = math::pow(m_phi,2)*m_forcing.dot(m_forcing);
        return math::pow(DeltaU.dot(DeltaU) + A0*math::pow(DeltaL,2.0),0.5);
    }

protected:

// Implementations for virtual functions
    void initMethods();
    void initiateStep();
    void iterationFinish();

    void quasiNewtonPredictor();
    void quasiNewtonIteration();

    void predictor();
    void predictorGuess();
    void iteration();

    void initOutput();
    void stepOutput();

// Additional internal functions
    void defaultOptions();
    void getOptions();

    void computeLambdas();
    void computeLambdasSimple();
    void computeLambdasModified();
    void computeLambdasComplex();
    void computeLambdasEta();

    void computeLambdaDET();

    void computeLambdaDOT();

    void computeLambdaMU();

protected:

    // Number of degrees of freedom
    using Base::m_numDof;

    using Base::m_jacobian;
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

    // Angle determination method: 0: determine based on previous load step. 1: determine based on previous iteration
    index_t m_angleDetermine;

    /// Scaling parameter
    T m_phi;
    bool m_phi_user;

    // MODIFIED ARC LENGTH METHOD
    /// factor (modified arc length method)
    T m_eta;

    // discriminant
    T m_discriminant;

    T m_alpha1;
    T m_alpha2;
    T m_alpha3;

    T m_a0;
    T m_b0,m_b1;
    T m_c0,m_c1,m_c2;

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

};


} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsALMCrisfield.hpp)
#endif
