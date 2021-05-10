 /** @file gsArcLengthIterator.h

    @brief Performs the arc length method to solve a nonlinear equation system.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#pragma once
#include <gsSpectra/gsSpectra.h>

namespace gismo
{

/**
    @brief Performs the arc length method to solve a nonlinear equation system.

    \tparam T coefficient type

    \ingroup gsStructuralAnalysis
*/
template <class T>
class gsArcLengthIterator
{
public:

    /// Constructor giving access to the gsShellAssembler object to create a linear system per iteration
    gsArcLengthIterator(  std::function < gsSparseMatrix<T> ( gsVector<T> const & ) > &Jacobian,
                          std::function < gsVector<T> ( gsVector<T> const &, T, gsVector<T> const & ) > &Residual,
                          gsVector<T> &Force )
    : m_jacobian(Jacobian),
      m_residualFun(Residual),
      m_forcing(Force)
    {
      this->defaultOptions();
      this->getOptions();

      // initialize variables
      m_numIterations = 0;
      m_arcLength = 1e-2;
      m_arcLength_prev = 1e-2;
      m_converged = false;

      // initialize errors
      m_basisResidualF = 0.0;
      m_basisResidualU = 0.0;

      m_initialized = false;

      // m_deltaLs = gsMatrix<T>::Zero(2,1);
    }


public:

    void initialize() {m_initialized = true; this -> initMethods(); this -> init();}

    void step();

    bool testSingularPoint(T tol, index_t kmax, bool jacobian=false);

    void computeSingularPoint(T singTol, index_t kmax, gsVector<T> U, T L, T tolE, T tolB=0, bool switchBranch=false);
    void computeSingularPoint(T singTol, index_t kmax, T tolE, T tolB=0, bool switchBranch=false);
    void computeSingularPoint(gsVector<T> U, T L, T tolE, T tolB=0, bool switchBranch=false);

    bool stabilityChange();

    void computeStability(gsVector<T> x, bool jacobian=true);
    index_t stability();
    index_t stability(gsVector<T> x, bool jacobian=true);

    void switchBranch();

public:

      void setLength(T length)
      {
        m_options.setReal("Length",length);
        m_arcLength = m_arcLength_prev = m_options.getReal("Length");
      }


    // Output
      /// Tells if the Arc Length method converged
      bool converged() const {return m_converged;}

      /// Returns the number of Newton iterations performed
      index_t numIterations() const { return m_numIterations;}

      /// Returns the tolerance value used
      T tolerance() const {return m_tolerance;}

      /// Returns the error after solving the nonlinear system
      T residue()   const {return m_residue;}

      /// Returns the value of the Determinant or other indicator
      T indicator() const {return m_indicator;}
      void setIndicator(T indicator) {m_indicator = indicator; m_stability = this->stability();}

      /// Return the solution vector and factor
      const gsVector<T> & solutionU() const {return  m_U;}
      const gsVector<T> & solutionDU() const {return  m_DeltaU;}
      T  solutionL() const {return  m_L;}
      T  solutionDL() const {return  m_DeltaL;}
      const gsVector<T> & solutionV() const {return  m_V;}

      // Returns if solution passed a bifurcation point
      bool isStable() const {return m_stability;}

      T determinant() const {return m_jacobian(m_U).toDense().determinant();}

    // Miscelaneous
      void resetStep() {m_DeltaUold.setZero();}

      // Set initial guess for solution
      void setInitialGuess(const gsVector<T> guess) {m_U = guess;}

      void setSolution(const gsVector<T> U, T L) {m_L = L; m_U = U; }// m_DeltaUold.setZero(); m_DeltaLold = 0;}

      gsOptionList & options() {return m_options;}

      void applyOptions() {this->getOptions(); }
      void resetOptions() {this->defaultOptions(); this->getOptions(); }

      void printSettings() { gsInfo<<m_options<<"\n"; }

public:
  // --------------------------------begin DEPRECIATED!!---------------------------------------------
  // Control parameters

      /// Set Arclength method
      void setMethod(std::string method)
      {
        if (method=="Riks")
          m_method = method::Riks;
        else if (method=="Crisfield")
          m_method = method::Crisfield;
        else if (method=="ConsistentCrisfield")
          m_method = method::ConsistentCrisfield;
        else if (method=="ExplicitIterations")
          m_method = method::ExplicitIterations;
        else if (method=="LoadControl")
          m_method = method::LoadControl;

        this -> initMethods();
      }


      /// Set the maximum number of Newton iterations allowed
      void setMaxIterations(index_t nIter) {m_maxIterations = nIter;}

      /// Set the tolerance for convergence
      void setTolerance(T tol) {m_tolerance = tol;}
      void setToleranceF(T tol) {m_toleranceF = tol;}
      void setToleranceU(T tol) {m_toleranceU = tol;}

      /// Set mode shape parameter tau
      void setTau(T tau) {m_tau = tau;}

      /// Set scaling factpr
      void setPhi(T phi) {m_phi = phi; m_phi_user = true;}

      /// Set arc length
      // void setLength(T length) {m_arcLength_prev = m_arcLength; m_arcLength = length;}
      void setLength(T length, bool adaptive)
      {
        m_arcLength_prev = m_arcLength;
        m_arcLength = length;
        m_adaptiveLength = adaptive;
        m_desiredIterations = 10;
      } // number of desired iterations defaults to 10
      void setLength(T length, bool adaptive, index_t iterations)
      {
        m_arcLength_prev = m_arcLength;
        m_arcLength = length;
        m_adaptiveLength = adaptive;
        m_desiredIterations = iterations;
      }
      void setLength(T length, index_t iterations)
      {
        m_arcLength_prev = m_arcLength;
        m_arcLength = length;
        m_adaptiveLength = true;
        m_desiredIterations = iterations;
      }

      // Use quasi-newton or not
      void quasiNewton() {m_quasiNewton = true; };
      void quasiNewton(index_t interval) {m_quasiNewton = true; m_quasiNewtonInterval = interval;}
    // Methods and procedures
      /// Set method for singular point detection
      void setBifurcationMethod(std::string method)
      {
        if (method=="Eigenvalue")
          m_bifurcationMethod = bifmethod::Eigenvalue;
        if (method=="Determinant")
          m_bifurcationMethod = bifmethod::Determinant;
      }

      /// Set AngleDetermination
      void setAngleDeterminationMethod(index_t method) {m_angleDetermine = method;}

      void setRelaxation(T relaxation) {m_relax = relaxation; }

    // Toggles
      // Verbose output
      void verbose() {m_verbose = true;}
  // --------------------------------end DEPRECIATED!!---------------------------------------------

protected:

    void init();
    void initMethods();

    void defaultOptions();
    void getOptions();

    void factorizeMatrix(const gsSparseMatrix<T> & M);
    gsVector<T> solveSystem(const gsVector<T> & F);

    void computeLambdas();
    void computeLambdasSimple();
    void computeLambdasModified();
    void computeLambdasComplex();
    void computeLambdasEta();

    void computeLambdaDET();

    void computeLambdaDOT();

    void computeLambdaMU();

    void extendedSystemIteration();

    index_t bisectionObjectiveFunction(const gsVector<T> & x, bool jacobian=true);

    T bisectionTerminationFunction(const gsVector<T> & x, bool jacobian=true);

    void computeResidual();
    void computeResidualNorms();
    void computeJacobian(gsVector<T> U);
    void computeJacobian();
    void computeLength();
    void computeUt();
    void computeUbar();

    // void extendedSystemSolve();
    // void extendedSystemSolve(T tol);
    // void extendedSystemSolve(const gsVector<T> U, T L);
    void extendedSystemSolve(const gsVector<T> U, T L, T tol);

    // void bisectionSolve();
    // void bisectionSolve(T tol);
    // void bisectionSolve(const gsVector<T> U, T L);
    void bisectionSolve(const gsVector<T> U, T L, T tol);


    void initOutput();
    void initOutputRiks();
    void initOutputCrisfield();
    void initOutputExtended();
    void initOutputLC();

    void stepOutput();
    void stepOutputRiks();
    void stepOutputCrisfield();
    void stepOutputExtended();
    void stepOutputLC();

    void iteration();
    void iterationRiks();
    void iterationConsistentCrisfield();
    void iterationExplicitIterations();
    void iterationCrisfield();
    void iterationLC();

    void initiateStep();
    void initiateStepRiks();
    void initiateStepConsistentCrisfield();
    void initiateStepExplicitIterations();
    void initiateStepCrisfield();
    void initiateStepLC();

    void predictor();
    void predictorRiks();
    void predictorConsistentCrisfield();
    void predictorExplicitIterations();
    void predictorCrisfield();
    void predictorLC();

    void iterationFinish();
    void iterationFinishRiks();
    void iterationFinishConsistentCrisfield();
    void iterationFinishExplicitIterations();
    void iterationFinishCrisfield();
    void iterationFinishLC();

protected:

    // Number of degrees of freedom
    index_t m_numDof;

    std::function < gsSparseMatrix<T> ( gsVector<T> const & ) > m_jacobian;
    std::function < gsMatrix<T> ( gsVector<T> const &, T, gsVector<T> const & ) > m_residualFun;
    std::function < gsMatrix<T> ( gsVector<T> const &, T, gsVector<T> const & ) > m_residualFunMod;
    gsVector<T> m_forcing;

    /// Linear solver employed
    gsSparseSolver<>::SimplicialLDLT  m_LDLTsolver;   // Cholesky
    gsSparseSolver<>::CGDiagonal      m_CGsolver;     // CG

protected:

    mutable gsOptionList m_options;

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


  protected:


    /// Number of Arc Length iterations performed
    index_t m_numIterations;

    /// Maximum number of Arc Length iterations allowed
    index_t m_maxIterations;

    /// Number of desired iterations
    index_t m_desiredIterations;

    /// Length of the step in the u,f plane
    T m_arcLength;
    T m_arcLength_prev;
    bool m_adaptiveLength;

    /// Scaling parameter
    T m_phi;
    bool m_phi_user;

    /// Tolerance value to decide convergence
    T m_tolerance;

    /// Tolerance value to decide convergence - Force criterion
    T m_toleranceF;

    /// Tolerance value to decide convergence - Displacement criterion
    T m_toleranceU;

    bool m_verbose;
    bool m_initialized;

    bool m_quasiNewton;
    index_t m_quasiNewtonInterval;

    std::string note;

protected:

    /// Convergence result
    bool m_converged;

    /// Force residuum
    T m_residue;

    /// Force residuum
    T m_residueF;
    T m_basisResidualF;

    /// Displacement residuum
    T m_residueU;
    T m_basisResidualU;

    /// Load residuum
    T m_residueL;

    /// Singular point
    T m_residueKTPhi;
    T m_basisResidualKTPhi;

    /// Indicator for bifurcation
    T m_indicator;
    index_t m_negatives;

    /// Relaxation factor
    T m_relax;

protected:

    // Arc length method (either "Crisfield" or "Riks")
    index_t m_method;

    // Previous update
    gsVector<T> m_DeltaUold;
    real_t m_DeltaLold;
    /// Displacement vector (present, at previously converged point)
    gsVector<T> m_U, m_Uprev;
    /// Update of displacement vector
    gsVector<T> m_DeltaU;
    /// u_bar
    gsVector<T> m_deltaUbar;
    /// u_t
    gsVector<T> m_deltaUt;
    /// Update of update of displacement vector
    gsVector<T> m_deltaU;

    /// Lambda (present, at previously converged point)
    T m_L, m_Lprev;
    /// Update of lambdaGeneralizedSelfAdjointEigenSolver
    T m_DeltaL;
    /// Update of update of lambda
    T m_deltaL;
    /// Vector with lambda updates
    gsVector<T> m_deltaLs;

    /// Jacobian matrix
    gsSparseMatrix<T> m_jacMat;
    T m_detKT;

    /// Value of residual function
    gsVector<T> m_resVec;

    // eigenvector
    gsVector<T> m_V;
    // step eigenvector
    gsVector<T> m_deltaV;
    gsVector<T> m_deltaVbar;
    gsVector<T> m_deltaVt;
    gsVector<T> m_DeltaV;

    // Stability indicator
    gsVector<T> m_stabilityVec;

    // Integer if point is unstable (-1) or not (+1)
    index_t m_stabilityPrev; //previous step
    index_t m_stability; //current step
    // Method to check if a point is a bifurcation point
    index_t m_bifurcationMethod;
    index_t m_solverType;

    // Angle determination method: 0: determine based on previous load step. 1: determine based on previous iteration
    index_t m_angleDetermine;

    // Branch switch parameter
    T m_tau;

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
#include GISMO_HPP_HEADER(gsArcLengthIterator.hpp)
#endif
