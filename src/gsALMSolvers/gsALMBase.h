 /** @file gsALMBase.h

    @brief Base class to perform the arc length method to solve a nonlinear equation system.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)

    TODO (June 2023):
    *    Change inputs to const references!
*/

#pragma once
#include <gsCore/gsLinearAlgebra.h>

#ifdef gsSpectra_ENABLED
#include <gsSpectra/gsSpectra.h>
#endif
#include <gsIO/gsOptionList.h>
#include <gsStructuralAnalysis/src/gsStructuralAnalysisTools/gsStructuralAnalysisTypes.h>

namespace gismo
{

/**
    @brief Performs the arc length method to solve a nonlinear system of equations.

    \tparam T coefficient type

    \ingroup gsStructuralAnalysis
*/
template <class T>
class gsALMBase
{
protected:

    typedef typename gsStructuralAnalysisOps<T>::ALResidual_t    ALResidual_t;
    typedef typename gsStructuralAnalysisOps<T>::Jacobian_t      Jacobian_t;
    typedef typename gsStructuralAnalysisOps<T>::dJacobian_t     dJacobian_t;

public:

    virtual ~gsALMBase() {};

    /// Constructor
    gsALMBase(  const Jacobian_t   & Jacobian,
                const ALResidual_t & ALResidual,
                const gsVector<T>  & Force )
    : m_residualFun(ALResidual),
      m_forcing(Force)
    {
        m_jacobian  = Jacobian;
        m_djacobian = [this](gsVector<T> const & x, gsVector<T> const & dx, gsSparseMatrix<T> & m) -> bool
        {
            return m_jacobian(x,m);
        };

        // initialize variables
        m_numIterations = 0;
        this->defaultOptions();
        this->setLength(1e-2);
        m_converged = false;

        // initialize errors
        m_basisResidualF = 0.0;
        m_basisResidualU = 0.0;

        m_status = gsStatus::NotStarted;
    }

    /// Constructor using the jacobian that takes the solution and the solution step
    gsALMBase(  const dJacobian_t &  dJacobian,
                const ALResidual_t & Residual,
                const gsVector<T>  & Force )
    : m_residualFun(Residual),
      m_forcing(Force)
    {
        m_djacobian  = dJacobian;

        // initialize variables
        m_numIterations = 0;
        this->defaultOptions();
        this->setLength(1e-2);
        m_converged = false;

        // initialize errors
        m_basisResidualF = 0.0;
        m_basisResidualU = 0.0;
    }

// General functions
public:

    virtual gsStatus status() { return m_status; }

    // Returns the number of DoFs
    virtual index_t numDofs() {return m_forcing.size();}

    // Returns the current length
    virtual T getLength() {return m_arcLength; }

    /// Perform one arc-length step
    virtual gsStatus step();

    /// Initialize the arc-length method, computes the stability of the initial configuration if \a stability is true
    virtual void initialize(bool stability = true)
    {
        m_initialized = true;
        this -> initMethods();
        this -> init(stability);
    }

    /// Set arc length to \a length
    virtual void setLength(T length)
    {
      m_options.setReal("Length",length);
      m_arcLength = m_arcLength_prev = m_arcLength_ori = m_options.getReal("Length");
    }

    /// Set arc length to \a length, enables \a adaptive steps
    virtual void setLength(T length, bool adaptive)
    {
      this->setLength(length);
      m_adaptiveLength = adaptive;
      m_desiredIterations = 10; // number of desired iterations defaults to 10
    }

    /// Set arc length to \a length, enables \a adaptive steps aiming for \a iterations number of iterations per step
    virtual void setLength(T length, bool adaptive, index_t iterations)
    {
      this->setLength(length);
      m_adaptiveLength = adaptive;
      m_desiredIterations = iterations;
    }
    /// Set arc length to \a length, enables adaptive steps aiming for \a iterations number of iterations per step
    virtual void setLength(T length, index_t iterations)
    {
      this->setLength(length);
      m_adaptiveLength = true;
      m_desiredIterations = iterations;
    }

    // Output
    /// True if the Arc Length method converged
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

    /// Returns if solution passed a bifurcation point
    virtual bool isStable() const {return m_stability;}

    // /// Returns the value of the deterimant of the jacobian
    // virtual T determinant() const
    // {
    //     return m_jacobian(m_U).toDense().determinant();
    // }

    /// Resets the step
    virtual void resetStep() {m_DeltaUold.setZero(); m_DeltaLold = 0;}

    // Set initial guess for solution
    virtual void setInitialGuess(const gsVector<T> & Uguess, const T & Lguess) {m_Uguess = Uguess; m_Lguess = Lguess;}
    virtual void setPrevious(const gsVector<T> & Uprev, const T & Lprev)
    {
        m_Uprev = Uprev;
        m_Lprev = Lprev;
        m_DeltaUold = m_U - m_Uprev;
        m_DeltaLold = m_L - m_Lprev;
    }

    /// Sets the solution
    virtual void setSolution(const gsVector<T> & U, const T & L) {m_L = L; m_U = U; }// m_DeltaUold.setZero(); m_DeltaLold = 0;}

    /// Sets the solution step
    virtual void setSolutionStep(const gsVector<T> & DU, const T & DL) {m_DeltaUold = DU; m_DeltaLold = DL;}// m_DeltaUold.setZero(); m_DeltaLold = 0;}

    /// Access the options
    virtual gsOptionList & options() {return m_options;};

    /// Set the options to \a options
    virtual void setOptions(gsOptionList options) {m_options.update(options,gsOptionList::addIfUnknown); this->getOptions(); };

    /// Return the options into \a options
    virtual const void options_into(gsOptionList options) {options = m_options;};

    /// Apply the options
    virtual void applyOptions() {this->getOptions(); }

    virtual T distance(const gsVector<T>& DeltaU, const T DeltaL) const
    {
        GISMO_NO_IMPLEMENTATION;
    }

// ------------------------------------------------------------------------------------------------------------
// ---------------------------------------Singular point methods-----------------------------------------------
// ------------------------------------------------------------------------------------------------------------
public:


    /**
     * @brief      Calculates the singular point.
     *
     * @param[in]  singTol       The tolerance for convergence
     * @param[in]  kmax          The maximum number of iterations for the initial power method
     * @param[in]  U             The point to start (displacements)
     * @param[in]  L             The point to start (loads)
     * @param[in]  tolE          The tolerance for the extended iterations
     * @param[in]  tolB          The tolerance for the bisection method
     * @param[in]  switchBranch  Switches branch if true
     * @param[in]  jacobian      Evaluate the Jacobian?
     */
    virtual gsStatus computeSingularPoint(const gsVector<T> & U, const T & L, bool switchBranch=false, bool jacobian=false, bool testPoint=true);

    /// Returns true if the point is a bifurcation
    virtual bool isBifurcation(bool jacobian = false);

    /// Checks if the stability of the system changed since the previously known solution
    virtual bool stabilityChange() const;

    /**
     * @brief      Calculates the stability of the solution \a x
     *
     * note; The shift is needed to ensure that a negative eigenvalue is found
     *
     * @param[in]  x         Solution vector
     * @param[in]  jacobian  Compute the jacobian?
     * @param[in]  shift     The shift to apply
     */
    virtual gsStatus computeStability(bool jacobian=true, T shift = -1e2);

    /// Computes the stability: -1 if unstable, +1 if stable
    virtual index_t stability() const;

    /// Switches branches
    virtual void switchBranch();

    /// Reduce the length by multiplication with a factor \a fac
    virtual T reduceLength(T fac = 0.5);

    /// Reset the length
    virtual T resetLength();

protected:

    /// See \a computeStability
    virtual void _computeStability(const gsVector<T> & x, bool jacobian=true, T shift = -1e2);

    /// See \a computeSingularPoint
    virtual void _computeSingularPoint(const gsVector<T> & U, const T & L, bool switchBranch=false, bool jacobian=false, bool testPoint=true);

    /**
     * @brief      Tests if a point is a bifurcation point
     *
     * @param[in]  jacobian  Evaluate the Jacobian?
     *
     * @return     True if it is a bifurcation point
     */
    virtual bool _testSingularPoint(bool jacobian=false);

    /// Perform an extended system iteration
    virtual void _extendedSystemIteration();

    /// Returns the objective function for the bisection method given solution \a x
    virtual index_t _bisectionObjectiveFunction(const gsVector<T> & x, bool jacobian=true);
    /// Returns the termination function for the bisection method given solution \a x
    virtual T       _bisectionTerminationFunction(const gsVector<T> & x, bool jacobian=true);

    /// Perform an extended system solve to find a singular point
    virtual bool _extendedSystemSolve(const gsVector<T> & U, const T L, const T tol);

    /// Perform a bisection system solve to find a singular point
    virtual bool _bisectionSolve(const gsVector<T> & U, const T L, const T tol);

    /// Initialize the output for extended iterations
    virtual void _initOutputExtended();

    /// Step output for extended iterations
    virtual void _stepOutputExtended();

// ------------------------------------------------------------------------------------------------------------
// ---------------------------------------Computations---------------------------------------------------------
// ------------------------------------------------------------------------------------------------------------
protected:
    /// Implementation of step
    virtual void _step();

    /// Set default options
    virtual void defaultOptions();

    /// Apply options
    virtual void getOptions();

    /// Initialize the solver
    virtual void init(bool stability);

    /// Factorize the matrix \a M
    virtual void factorizeMatrix(const gsSparseMatrix<T> & M);

    /// Solve the system with right-hand side \a F
    virtual gsVector<T> solveSystem(const gsVector<T> & F);

    /// Compute the residual
    virtual gsVector<T> computeResidual(const gsVector<T> & U, const T & L);
    virtual void computeResidual();

    /// Compute the residual error norms
    virtual void computeResidualNorms();

    /// Compute the jacobian matrix
    virtual gsSparseMatrix<T> _computeJacobian(const gsVector<T> & U, const gsVector<T> & dU);
    virtual gsSparseMatrix<T> computeJacobian(const gsVector<T> & U, const gsVector<T> & dU);
    virtual gsSparseMatrix<T> computeJacobian(const gsVector<T> & U);
    virtual gsSparseMatrix<T> computeJacobian();

    /// Compute the adaptive arc-length
    virtual void computeLength();

    /// Compute \f$\u_t\f$
    virtual void computeUt();
    /// Compute \f$\bar_u\f$
    virtual void computeUbar();

// Purely virtual functions
protected:
    /// Initialize the ALM
    virtual void initMethods() = 0;
    /// Initiate the first iteration
    virtual void initiateStep() = 0;
    /// Finish the iterations
    virtual void iterationFinish() = 0;

    /// Provide a specialized predictor when using quasi newton methods
    virtual void quasiNewtonPredictor() = 0;
    /// Perform iteration using quasi-newton method
    virtual void quasiNewtonIteration() = 0;

    /// Step predictor
    virtual void predictor() = 0;
    virtual void predictorGuess() = 0;
    /// A single iteration
    virtual void iteration() = 0;

    /// Initialize the output
    virtual void initOutput() = 0;
    /// Provide step-wise output
    virtual void stepOutput() = 0;

protected:

    // Number of degrees of freedom
    index_t m_numDof;

    Jacobian_t      m_jacobian;
    dJacobian_t     m_djacobian;
    const ALResidual_t    m_residualFun;
    const gsVector<T>     m_forcing;

    mutable typename gsSparseSolver<T>::uPtr m_solver; // Cholesky by default

public:


    struct bifmethod
    {
        enum type
        {
            Nothing = -1,
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
    mutable gsOptionList m_options;


    /// Number of Arc Length iterations performed
    index_t m_numIterations;

    /// Maximum number of Arc Length iterations allowed
    index_t m_maxIterations;

    /// Number of desired iterations
    index_t m_desiredIterations;

    /// Length of the step in the u,f plane
    T m_arcLength;
    T m_arcLength_prev;
    T m_arcLength_ori;
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

    gsStatus m_status;

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

    // What to do after computeSingularPoint fails?
    index_t m_SPfail;

    // Number of iterations and the tolerance for the singular point test
    index_t m_SPTestIt;
    T m_SPTestTol;

    // Singular point computation tolerances
    T m_SPCompTolE; // extended iterations
    T m_SPCompTolB; // bisection method

    // Branch switch parameter
    T m_tau;
};


} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsALMBase.hpp)
#endif
