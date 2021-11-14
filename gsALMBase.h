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
class gsALMBase
{
public:

    virtual ~gsALMBase() {};

    /// Constructor giving access to the gsShellAssembler object to create a linear system per iteration
    gsALMBase(  std::function < gsSparseMatrix<T> ( gsVector<T> const & ) > &Jacobian,
                std::function < gsVector<T> ( gsVector<T> const &, T, gsVector<T> const & ) > &Residual,
                gsVector<T> &Force )
    : m_jacobian(Jacobian),
      m_residualFun(Residual),
      m_forcing(Force)
    {
        // initialize variables
        m_numIterations = 0;
        m_arcLength = m_arcLength_prev = 1e-2;
        m_converged = false;

        // initialize errors
        m_basisResidualF = 0.0;
        m_basisResidualU = 0.0;
    }

// General functions
public:
    virtual index_t numDofs() {return m_forcing.size();}
    virtual T getLength() {return m_arcLength; }

    virtual void step();

    virtual void initialize(bool stability = true)
    {
        m_initialized = true;
        this -> initMethods();
        this -> init(stability);
    }



    /// Set arc length
    virtual void setLength(T length)
    {
      m_options.setReal("Length",length);
      m_arcLength = m_arcLength_prev = m_options.getReal("Length");
    }

    virtual void setLength(T length, bool adaptive)
    {
      m_arcLength_prev = m_arcLength;
      m_arcLength = length;
      m_adaptiveLength = adaptive;
      m_desiredIterations = 10;
    } // number of desired iterations defaults to 10
    virtual void setLength(T length, bool adaptive, index_t iterations)
    {
      m_arcLength_prev = m_arcLength;
      m_arcLength = length;
      m_adaptiveLength = adaptive;
      m_desiredIterations = iterations;
    }
    virtual void setLength(T length, index_t iterations)
    {
      m_arcLength_prev = m_arcLength;
      m_arcLength = length;
      m_adaptiveLength = true;
      m_desiredIterations = iterations;
    }

    // Output
    /// Tells if the Arc Length method converged
    virtual bool converged() const {return m_converged;}

    /// Returns the number of Newton iterations performed
    virtual index_t numIterations() const { return m_numIterations;}

    /// Returns the tolerance value used
    virtual T tolerance() const {return m_tolerance;}

    /// Returns the error after solving the nonlinear system
    virtual T residue()   const {return m_residue;}

    /// Returns the value of the Determinant or other indicator
    virtual T indicator() const {return m_indicator;}
    virtual void setIndicator(T indicator) {m_indicator = indicator; m_stability = this->stability();}

    /// Return the solution vector and factor
    virtual const gsVector<T> & solutionU() const {return  m_U;}
    virtual const gsVector<T> & solutionDU() const {return  m_DeltaU;}
    virtual T  solutionL() const {return  m_L;}
    virtual T  solutionDL() const {return  m_DeltaL;}
    virtual const gsVector<T> & solutionV() const {return  m_V;}

    // Returns if solution passed a bifurcation point
    virtual bool isStable() const {return m_stability;}

    virtual T determinant() const {return m_jacobian(m_U).toDense().determinant();}

    // Miscelaneous
    virtual void resetStep() {m_DeltaUold.setZero(); m_DeltaLold = 0;}

    // Set initial guess for solution
    // virtual void setInitialGuess(const gsVector<T> guess) {m_U = guess;}
    virtual void setInitialGuess(const gsVector<T> Uguess, T Lguess) {m_Uguess = Uguess; m_Lguess = Lguess;}

    virtual void setSolution(const gsVector<T> U, T L) {m_L = L; m_U = U; }// m_DeltaUold.setZero(); m_DeltaLold = 0;}
    virtual void setSolutionStep(const gsVector<T> DU, T DL) {m_DeltaUold = DU; m_DeltaLold = DL;}// m_DeltaUold.setZero(); m_DeltaLold = 0;}

    virtual gsOptionList & options() {return m_options;};
    virtual void setOptions(gsOptionList options) { m_options = options; this->getOptions(); };
    virtual const void options_into(gsOptionList options) {options = m_options;};

    virtual void applyOptions() {this->getOptions(); }

    virtual T distance(gsVector<T>& DeltaU, T DeltaL)
    {
        GISMO_NO_IMPLEMENTATION;
    }

// ------------------------------------------------------------------------------------------------------------
// ---------------------------------------Singular point methods-----------------------------------------------
// ------------------------------------------------------------------------------------------------------------
public:
    virtual void computeSingularPoint(T singTol, index_t kmax, gsVector<T> U, T L, T tolE, T tolB=0, bool switchBranch=false, bool jacobian=false);
    virtual void computeSingularPoint(T singTol, index_t kmax, T tolE, T tolB=0, bool switchBranch=false, bool jacobian=false);
    virtual void computeSingularPoint(gsVector<T> U, T L, T tolE, T tolB=0, bool switchBranch=false, bool jacobian=false);

    virtual bool testSingularPoint(T tol, index_t kmax, bool jacobian=false);
    virtual bool testSingularPoint(gsVector<T> U, T L, T tol, index_t kmax, bool jacobian=false);

    virtual bool stabilityChange();

    // The shift is needed to ensure that a negative eigenvalue is found
    virtual void computeStability(gsVector<T> x, bool jacobian=true, T shift = -1e2);
    virtual gsMatrix<T> computeModes(gsVector<T> x, bool jacobian=true, T shift = -1e2);
    virtual index_t stability();
    virtual index_t stability(gsVector<T> x, bool jacobian=true);

    virtual void switchBranch();

protected:
    virtual void extendedSystemIteration();

    virtual index_t bisectionObjectiveFunction(const gsVector<T> & x, bool jacobian=true);
    virtual T bisectionTerminationFunction(const gsVector<T> & x, bool jacobian=true);

    // void extendedSystemSolve();
    // void extendedSystemSolve(T tol);
    // void extendedSystemSolve(const gsVector<T> U, T L);
    virtual void extendedSystemSolve(const gsVector<T> U, T L, T tol);

    // void bisectionSolve();
    // void bisectionSolve(T tol);
    // void bisectionSolve(const gsVector<T> U, T L);
    virtual void bisectionSolve(const gsVector<T> U, T L, T tol);

    virtual void initOutputExtended();
    virtual void stepOutputExtended();

// ------------------------------------------------------------------------------------------------------------
// ---------------------------------------Computations---------------------------------------------------------
// ------------------------------------------------------------------------------------------------------------
protected:

    virtual void defaultOptions();
    virtual void getOptions();

    virtual void init(bool stability);

    virtual void factorizeMatrix(const gsSparseMatrix<T> & M);
    virtual gsVector<T> solveSystem(const gsVector<T> & F);

    virtual void computeResidual();
    virtual void computeResidualNorms();
    virtual void computeJacobian(gsVector<T> U);
    virtual void computeJacobian();
    virtual void computeLength();
    virtual void computeUt();
    virtual void computeUbar();

// Purely virtual functions
protected:
    virtual void initMethods() = 0;
    virtual void initiateStep() = 0;
    virtual void iterationFinish() = 0;

    virtual void quasiNewtonPredictor() = 0;
    virtual void quasiNewtonIteration() = 0;

    virtual void predictor() = 0;
    virtual void predictorGuess() = 0;
    virtual void iteration() = 0;

    virtual void initOutput() = 0;
    virtual void stepOutput() = 0;

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

    struct bifmethod
    {
        enum type
        {
            Determinant = 0,
            Eigenvalue  = 1,
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

    std::string m_note;

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

    gsVector<T> m_Uguess;
    T m_Lguess;

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

    // What to do after computeSingularPoint fails?
    index_t m_SPfail;

    // Branch switch parameter
    T m_tau;
};


} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsALMBase.hpp)
#endif