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


public:

    void initialize() {m_initialized = true; this -> initMethods(); this -> init();}

    void step();

protected:

    void init();
    void initMethods();

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

    void initOutput();

    void stepOutput();

    void iteration();

    void initiateStep();

    void predictor();

    void iterationFinish();

protected:

    // Number of degrees of freedom
    using Base::m_numDof;

    using Base::m_jacobian;
    using Base::m_residualFun;
    using Base::m_residualFunMod;
    using Base::m_forcing;

    /// Linear solver employed
    using Base::m_LDLTsolver;   // Cholesky
    using Base::m_CGsolver;     // CG

protected:

/// @brief Specifies the material law to use
    struct method
    {
        enum type
        {
            LoadControl  = 0,
            Riks  = 1,
            Crisfield = 2,
            ConsistentCrisfield = 3,
            ExplicitIterations = 4,
        };
    };

    struct bifmethod
    {
        enum type
        {
            Determinant = 0,
            Eigenvalue  = 1,
        };
    };

    struct angmethod
    {
        enum type
        {
            Step = 0,
            Iteration  = 1,
            Predictor  = 2,
        };
    };

    struct solver
    {
        enum type
        {
            LDLT = 0,
            CG  = 1, // The CG solver is robust for membrane models, where zero-blocks in the matrix might occur.
        };
    };

    struct SPfail
    {
        enum type
        {
            Without  = 0,
            With     = 1,
        };
    };


  protected:
    using Base::m_options;

    /// Number of Arc Length iterations performed
    using Base::m_numIterations;

    /// Maximum number of Arc Length iterations allowed
    using Base::m_maxIterations;

    /// Number of desired iterations
    using Base::m_desiredIterations;

    /// Length of the step in the u,f plane
    using Base::m_arcLength;
    using Base::m_arcLength_prev;
    using Base::m_adaptiveLength;

    /// Scaling parameter
    T m_phi;
    bool m_phi_user;

    /// Tolerance value to decide convergence
    using Base::m_tolerance;

    /// Tolerance value to decide convergence - Force criterion
    using Base::m_toleranceF;

    /// Tolerance value to decide convergence - Displacement criterion
    using Base::m_toleranceU;

    using Base::m_verbose;
    using Base::m_initialized;

    using Base::m_quasiNewton;
    using Base::m_quasiNewtonInterval;

    using Base::m_note;

    /// Convergence result
    using Base::m_converged;

    /// Force residuum
    using Base::m_residue;

    /// Force residuum
    using Base::m_residueF;
    using Base::m_basisResidualF;

    /// Displacement residuum
    using Base::m_residueU;
    using Base::m_basisResidualU;

    /// Load residuum
    using Base::m_residueL;

    /// Singular point
    using Base::m_residueKTPhi;
    using Base::m_basisResidualKTPhi;

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

    /// Value of residual function
    using Base::m_resVec;

    // eigenvector
    using Base::m_V;
    // step eigenvector
    using Base::m_deltaV;
    using Base::m_deltaVbar;
    using Base::m_deltaVt;
    using Base::m_DeltaV;

    // Stability indicator
    using Base::m_stabilityVec;

    // Integer if point is unstable (-1) or not (+1)
    using Base::m_stabilityPrev; //previous step
    using Base::m_stability; //current step
    // Method to check if a point is a bifurcation point
    using Base::m_bifurcationMethod;
    using Base::m_solverType;

    // Angle determination method: 0: determine based on previous load step. 1: determine based on previous iteration
    using Base::m_angleDetermine;

    // What to do after computeSingularPoint fails?
    using Base::m_SPfail;

    // Branch switch parameter
    using Base::m_tau;

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

    T m_tangentLength;
};


} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsALMCrisfield.hpp)
#endif
