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
#include <functional>

#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsMemory.h>

#ifdef gsSpectra_ENABLED
#include <gsSpectra/gsSpectra.h>
#endif
#include <gsIO/gsOptionList.h>
#include <gsStructuralAnalysis/src/gsStructuralAnalysisTools/gsStructuralAnalysisTypes.h>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsStabilityAnalysis.h>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsContinuationDriver.h>

namespace gismo
{

/**
    @brief Performs the arc length method to solve a nonlinear system of equations.

    One step is one step of a gsContinuationDriver, which factorizes the
    tangent once per Jacobian and reports failure through a status instead of
    exceptions. The path constraint that defines the method is installed with
    \a setConstraint, normally in the constructor of a subclass (the presets
    gsALMLoadControl, gsALMRiks, gsALMCrisfield and gsALMConsistentCrisfield).
    Without a constraint, \a step returns gsStatus::OtherError and sets a note.

    The state lives in a gsStateSpace. The class either owns a
    gsVectorStateSpace and a gsContinuationProblem built from the Jacobian and
    residual operators of the first two constructors, or refers to a caller's
    problem and state space (non-owning constructor).

    <b>Sign convention.</b> The operators passed to the legacy constructors
    use the residual \f$F_{int}-\lambda F_{ext}\f$ and the Newton correction
    solves \f$K\,\delta U=-R\f$. The driver works with
    \f$R=\lambda F_{ext}-F_{int}\f$ and solves \f$K\,\delta U=R\f$, so the
    problem built from the legacy operators negates the residual (exactly) and
    passes the tangent through. The non-owning constructor expects a problem
    that is already in the driver convention.

    The class is not copyable: the operators of the owned problem refer to the
    object itself.

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

    /// Constructor. The class owns the state space and the problem that
    /// wrap \a Jacobian, \a ALResidual and \a Force.
    gsALMBase(  const Jacobian_t   & Jacobian,
                const ALResidual_t & ALResidual,
                const gsVector<T>  & Force )
    : m_residualFun(ALResidual),
      m_forcing(Force)
    {
        m_jacobian  = Jacobian;
        m_djacobian = [this](gsVector<T> const & x, gsVector<T> const & /*dx*/, gsSparseMatrix<T> & m) -> bool
        {
            return m_jacobian(x,m);
        };

        this->_constructOwned(Force);
        this->_initMembers(false);
    }

    /// Constructor using the jacobian that takes the solution and the solution step
    gsALMBase(  const dJacobian_t &  dJacobian,
                const ALResidual_t & Residual,
                const gsVector<T>  & Force )
    : m_residualFun(Residual),
      m_forcing(Force)
    {
        m_djacobian  = dJacobian;

        this->_constructOwned(Force);
        this->_initMembers(false);
    }

    /**
       @brief Constructor on a caller's problem and state space (non-owning).

       No Jacobian or residual operator of the legacy signature exists. A path
       constraint must be installed (see \a setConstraint) for \a step to
       run; without one, \a step returns gsStatus::OtherError.

       @param[in] problem the problem, already in the driver convention
                  \f$R=\lambda F-F_{int}\f$, \f$K=\partial F_{int}/\partial U\f$;
                  nothing is negated. Must outlive the solver.
       @param[in] space   the state space; must outlive the solver. The
                  solution vector \a solutionU() and the extended-system
                  functions (singular points, branch switching) need a
                  copyable space whose exported vector has the tangent size.
    */
    gsALMBase(  gsContinuationProblem<T> & problem,
                gsStateSpace<T>          & space )
    : m_forcing(problem.F)
    {
        m_problem = &problem;
        m_space   = &space;
        this->_initMembers(true);
    }

// General functions
public:

    virtual gsStatus status() { return m_status; }

    // Returns the number of DoFs
    virtual index_t numDofs() {return m_forcing.size();}

    /// Returns the current length (the length of the step controller)
    virtual T getLength();

    /**
       @brief Perform one arc-length step.

       Runs one gsContinuationDriver step: one factorization per assembled
       Jacobian and no exception on the step path. A failed step leaves the
       solution and the stability unchanged.

       With the \a BifurcationMethod Determinant or Eigenvalue the stability
       of the new point is computed once per accepted step, from one assembly
       and factorization of the tangent at the committed point (the converged
       point, not the last iterate). The next predictor reuses that
       factorization when its constraint needs the load tangent (load control,
       Crisfield), so from the second step on such a step costs no more than
       without the stability evaluation; the first step, and the first step
       after anything that discards the committed tangent (\a initialize,
       \a resetStep, \a switchBranch, \a computeSingularPoint, a change of
       options), pays one assembly and factorization more. The Riks and the
       consistent Crisfield predictors need no tangent, so they pay one
       assembly and factorization more per step. With
       the Determinant method the LDL^T pivots of that factorization are
       classified, with the Eigenvalue method the dense eigenvalue problem of
       that tangent is solved. A committed step cannot be undone, so a
       failing stability evaluation does not change the returned status
       (gsStatus::Success); \a note reports it and \a stability keeps its
       previous value. With the method Nothing no extra work is done.

       If \a setInitialGuess was called with a guess that differs from the
       current point, the guess is applied before the step by
       \a applyInitialGuess.

       If a callback was set with \a onStepAccepted, it is called after every
       accepted step, after the solution has been refreshed.

       @return the status, also available as \a status()
    */
    virtual gsStatus step();

    /// Initialize the arc-length method, computes the stability of the initial configuration if \a stability is true.
    /// The solution is reset to zero (an owned state space) or taken from the
    /// state space (a caller's state space), and the step history is cleared.
    /// Without an installed constraint the stability is not computed (a warning is printed).
    virtual void initialize(bool stability = true);

    /// Set arc length to \a length
    virtual void setLength(T length);

    /// Set arc length to \a length, enables \a adaptive steps
    virtual void setLength(T length, bool adaptive);

    /// Set arc length to \a length, enables \a adaptive steps aiming for \a iterations number of iterations per step
    virtual void setLength(T length, bool adaptive, index_t iterations);

    /// Set arc length to \a length, enables adaptive steps aiming for \a iterations number of iterations per step
    virtual void setLength(T length, index_t iterations);

    /// Sets the function that is called after every accepted step with the
    /// number of steps accepted before it (0, 1, ...; reset by \a initialize)
    /// and the load factor of the new point. An empty function clears it.
    /// It is not called for the trial steps of the bisection inside
    /// computeSingularPoint.
    virtual void onStepAccepted(const std::function<void(index_t, T)> & f) { m_onAccepted = f; }

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

    /// Return the solution vector and factor.
    /// \a solutionU is a member buffer, refreshed in place from
    /// gsStateSpace::exportWorking after every accepted step; it is empty for
    /// a state space that is not copyable.
    virtual const gsVector<T> & solutionU() const {return  m_U;}
    virtual const gsVector<T> & solutionDU() const {return  m_DeltaU;}
    virtual T  solutionL() const {return  m_L;}
    virtual T  solutionDL() const {return  m_DeltaL;}
    virtual const gsVector<T> & solutionV() const {return  m_V;}

    /// Returns true iff the last computed (or set, see \ref setIndicator) stability is +1.
    /// It is false if the solution is unstable, or if the stability has not been computed
    /// (after construction, or after \ref initialize with \a stability = false).
    virtual bool isStable() const {return m_stability > 0;}

    // /// Returns the value of the deterimant of the jacobian
    // virtual T determinant() const
    // {
    //     return m_jacobian(m_U).toDense().determinant();
    // }

    /// Resets the step. This clears the step history, so the next step uses
    /// the first-step predictor.
    virtual void resetStep();

    // Set initial guess for solution
    virtual void setInitialGuess(const gsVector<T> & Uguess, const T & Lguess) {m_Uguess = Uguess; m_Lguess = Lguess;}

    /// Sets the previous solution. Call it after \a setSolution: the previous step
    /// is the difference between the present and the previous solution.
    virtual void setPrevious(const gsVector<T> & Uprev, const T & Lprev);

    /**
       @brief Sets the solution.

       The solution is imported into the state space, which requires a
       copyable state space (GISMO_ENSURE) and a vector of the tangent size
       (without an installed constraint only the mirror is set). The step history becomes the difference between the new and the previous solution, with the previous
       step length, as long as a previous solution is known (after an
       accepted step or \a setPrevious); otherwise it is cleared.
       Note that after accepted steps, a call without a following \a setPrevious
       sets the history to \f$U-U_{prev}\f$ also for the Crisfield method.

       On an initialized solver the stability is evaluated once at the new
       point (once per call, not per step), by the same assembly and
       factorization of the tangent at the committed point as after an
       accepted step, with every \a BifurcationMethod. It is stored as the
       previous stability of the next step, so that a bifurcation in the
       first step after the call is detected, and the factorization is reused
       by the predictor of that step when its constraint needs the load
       tangent. Not evaluated for the trial
       steps of the bisection in \a computeSingularPoint.
    */
    virtual void setSolution(const gsVector<T> & U, const T & L);

    /// Sets the solution step
    virtual void setSolutionStep(const gsVector<T> & DU, const T & DL);

    /// Access the options
    virtual gsOptionList & options() {return m_options;};

    /// Set the options to \a options
    virtual void setOptions(gsOptionList options) {m_options.update(options,gsOptionList::addIfUnknown); this->getOptions(); };

    /// Return the options into \a options
    virtual void options_into(gsOptionList options) {options = m_options;};

    /// Apply the options
    virtual void applyOptions() {this->getOptions(); }

    /// Distance of the increment (\a DeltaU, \a DeltaL) in the metric of the
    /// installed path constraint; not implemented without one. It is defined
    /// before the first step.
    virtual T distance(const gsVector<T>& DeltaU, const T DeltaL) const;

    /// Frozen arc-length scaling \f$\varphi_0\f$ of the run, for methods that freeze one
    /// (gsALMConsistentCrisfield with automatic "Scaling"); 0 otherwise.
    /// Solvers that continue one path in pieces (e.g. the workers of gsAPALM) exchange it
    /// so that every piece uses the same arc metric.
    virtual T frozenScaling() const { return T(0); }

    /// True iff a frozen scaling is set (see \a frozenScaling). Default: false.
    virtual bool hasFrozenScaling() const { return false; }

    /// Imposes the frozen scaling \a phi0 (> 0, finite) obtained from another solver's
    /// \a frozenScaling. Methods without a frozen scaling ignore it (default).
    virtual void setFrozenScaling(T /*phi0*/) { }

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
     *
     * The state space must be copyable and export a vector of the
     * tangent size; otherwise a warning is printed and gsStatus::OtherError is returned.
     * The resulting point is imported into the state space and the driver
     * (load factor and step history, see \a setSolution). The trial steps of the
     * bisection do not call the \a onStepAccepted callback.
     */
    virtual gsStatus computeSingularPoint(const gsVector<T> & U, const T & L, bool switchBranch=false, bool jacobian=false, bool testPoint=true);

    /// Returns true if the point is a bifurcation. The state
    /// space must be copyable and export a vector of the tangent size; otherwise
    /// a warning is printed and false is returned. Without \a jacobian the stored
    /// tangent is used: after an accepted step with the \a BifurcationMethod
    /// Determinant or Eigenvalue (or after \a setSolution) the tangent at the
    /// committed point; with the method Nothing the last assembled tangent of
    /// that step; after a failed step the last tangent of the failed attempt.
    virtual bool isBifurcation(bool jacobian = false);

    /// Checks if the stability of the system changed since the previously known solution.
    /// It is false while either the current or the previous stability is unknown (0).
    ///
    /// The previous stability of a step is the one stored at the end of the previous
    /// accepted step (or evaluated by \a setSolution), at the committed point (from
    /// the tangent assembled there, never from the last Newton iterate); it is
    /// never recomputed at the start of a step and does not depend on the verbosity.
    /// After \a switchBranch it is unknown, so the first step after a branch switch
    /// reports no change and only seeds the stored value (as AUTO and pde2path do).
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

    /// Switches branches. The new point is imported into the
    /// state space and the step history is cleared, so the next step uses the
    /// first-step predictor of the constraint. The stored stability becomes unknown (0),
    /// so \a stabilityChange is false in the first step after the switch, whose
    /// accepted value seeds the stored stability. Needs a copyable state space.
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
    /// Set default options
    virtual void defaultOptions();

    /// Apply options; ends by mapping them to the driver (convergence test,
    /// step controller, corrector, tangent solver), so subclasses that override
    /// it call this version first
    virtual void getOptions();

    /// Initialize the solver
    virtual void init(bool stability);

    /// Factorize the matrix \a M
    virtual void factorizeMatrix(const gsSparseMatrix<T> & M);

    /// Solve the system with right-hand side \a F
    virtual gsVector<T> solveSystem(const gsVector<T> & F);

    /// Compute the residual (sign of ALResidual_t). Without legacy operators (solver
    /// built on a caller's problem) the residual is evaluated through the state
    /// space at \a U, which costs O(n) extra per call.
    virtual gsVector<T> computeResidual(const gsVector<T> & U, const T & L);
    virtual void computeResidual();

    /// Compute the residual error norms
    virtual void computeResidualNorms();

    /// Compute and factorize the jacobian matrix at the displacement \a U for
    /// the trial increment \a m_DeltaL, i.e. at the load factor \a m_L + \a m_DeltaL
    /// (use \a _computeJacobianAt for a given load factor). Without legacy
    /// operators (solver built on a caller's problem) the tangent is evaluated
    /// through the state space, which costs O(n) extra per call.
    virtual gsSparseMatrix<T> _computeJacobian(const gsVector<T> & U, const gsVector<T> & dU);
    /// Compute and factorize the jacobian matrix at the displacement \a U and
    /// the load factor \a L, with \a dU the latest corrector update. Operators
    /// built on the legacy signatures do not depend on \a L. Evaluations at a
    /// committed point (stability, singular points) pass the committed load
    /// factor \a m_L.
    gsSparseMatrix<T> _computeJacobianAt(const gsVector<T> & U, const gsVector<T> & dU, T L);
    virtual gsSparseMatrix<T> computeJacobian(const gsVector<T> & U, const gsVector<T> & dU);
    virtual gsSparseMatrix<T> computeJacobian(const gsVector<T> & U);
    virtual gsSparseMatrix<T> computeJacobian();

    /// Compute the adaptive arc-length
    virtual void computeLength();

    /// Compute \f$\u_t\f$
    virtual void computeUt();
    /// Compute \f$\bar_u\f$
    virtual void computeUbar();

// ------------------------------------------------------------------------------------------------------------
// ---------------------------------------Driver mode----------------------------------------------------------
// ------------------------------------------------------------------------------------------------------------
protected:

    /**
       @brief Installs the path constraint that defines the method.

       The constraint is handed to the driver, whose setup is marked dirty.
       Install it before \a initialize. Subclasses call it in their
       constructors and map their options to the constraint in \a getOptions
       (after calling the base \a getOptions, which configures the driver).

       @param[in] c the constraint; must not be null (GISMO_ENSURE)
    */
    void setConstraint(const typename gsPathConstraint<T>::Ptr & c);

    /// The installed path constraint; null until \a setConstraint is called
    const typename gsPathConstraint<T>::Ptr & constraint() const { return m_constraint; }

    /// The continuation driver; exists from construction on
    gsContinuationDriver<T> &       driver()       { return *m_driver; }
    /// The continuation driver; exists from construction on
    const gsContinuationDriver<T> & driver() const { return *m_driver; }

    /// The state space of the driver
    gsStateSpace<T> &       stateSpace()       { return *m_space; }
    /// The state space of the driver
    const gsStateSpace<T> & stateSpace() const { return *m_space; }

    /**
       @brief Applies the initial guess before a step.

       The default sets the step history to
       \f$(U_{guess}-U,\ \lambda_{guess}-\lambda)\f$ with the previous step
       length, so that a secant predictor scales it by \f$s/s_{prev}\f$ (the
       behaviour of the Riks and consistent Crisfield methods). The guess is
       not cleared by the default, so a guess that is not replaced is applied
       again at the next step while it differs from the current point.
       Subclasses refine it, e.g. load control uses the load difference as a
       one-off step.

       @param[in] Uguess guess of the solution
       @param[in] Lguess guess of the load factor
    */
    virtual void applyInitialGuess(const gsVector<T> & Uguess, T Lguess);

    /// Reset hook called first by \a initialize; subclasses reset their method state (e.g. the frozen scaling)
    virtual void initMethods() {}

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

    /// Tolerance on the normalised path-constraint residual (constraintResidual())
    /// that an iterate must satisfy, in addition to the force and displacement
    /// criteria, to be accepted
    T m_toleranceC;

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
    gsStabilityAnalysis<T> m_stabilityAnalysis;

    // Stability of the point: -1 unstable, +1 stable, 0 not computed (see isStable)
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

protected:

    /// Problem built from the legacy operators (legacy constructors)
    gsContinuationProblem<T>                  m_ownProblem;
    /// State space owned by the legacy constructors
    memory::unique_ptr<gsVectorStateSpace<T> > m_ownSpace;
    /// The problem the driver works on: \a m_ownProblem or the caller's
    gsContinuationProblem<T> *                m_problem;
    /// The state space the driver works on: \a m_ownSpace or the caller's
    gsStateSpace<T> *                         m_space;
    /// The driver, built on \a m_problem and \a m_space
    memory::unique_ptr<gsContinuationDriver<T> > m_driver;
    /// The installed path constraint; null until \a setConstraint is called
    typename gsPathConstraint<T>::Ptr         m_constraint;
    /// True if the driver configuration changed since its last setup()
    bool                                      m_driverDirty;
    /// True if \a m_Uprev, \a m_Lprev hold a previous solution
    bool                                      m_prevKnown;
    /// Number of accepted steps since \a initialize
    index_t                                   m_numAccepted;
    /// Callback of \a onStepAccepted
    std::function<void(index_t, T)>           m_onAccepted;

private:

    gsALMBase(const gsALMBase &);
    gsALMBase & operator=(const gsALMBase &);

    /// Builds the owned state space and problem from the legacy operators
    void _constructOwned(const gsVector<T> & Force);
    /// Initialises the members, builds the driver and its observer
    void _initMembers(bool readOptions);
    /// Maps the options to the driver, called at the end of \a getOptions
    void _configureDriver();
    /// One step
    gsStatus _stepDriver();
    /// Driver-mode part of \a initialize
    void _initializeDriver();
    /// Re-reads the solution mirror \a m_U from the state space
    void _refreshMirror();
    /// True if the state space is copyable and exports a vector of the tangent size
    bool _vectorLike() const;
    /// Imports \a m_U, \a m_L into the state space and the driver; false if the state space refuses \a m_U
    bool _importMirror();
    /// Writes the step history (\a DU, \a DL, previous length) into the driver
    void _setHistory(const gsVector<T> & DU, T DL);
    /// Sets the history from the mirrors: present minus previous solution, or cleared if the latter is unknown
    void _setHistoryFromMirror();
    /// Keeps the step size of a valid history equal to the previous controller length
    void _syncHistoryStep();
    /// Evaluates the residual (legacy sign) or Jacobian at \a U through the state space
    bool _evaluateAt(const gsVector<T> & U, T L, gsVector<T> * R, gsSparseMatrix<T> * K);
    /// Prints the verbose output of a driver event
    void _observe(const gsStepEvent<T> & ev) const;
    /// Stability of the point just accepted: \a _stabilityAtCommitted, or
    /// \a computeStability(false) on the last tangent with the method Nothing
    void _stabilityAfterStep();
    /// Stability at the committed point after an external \a setSolution:
    /// \a _stabilityAtCommitted, stored as current and previous stability
    void _stabilityAtSolution();
    /// Assembles and factorizes the tangent at the committed point through the
    /// driver (reused by the next predictor), stores it in \a m_jacMat and sets
    /// the stability from the LDL^T pivots (Determinant) or from \a m_jacMat
    /// (otherwise); false if the assembly or the analysis fails
    bool _stabilityAtCommitted();

    /// True while the bisection of \a computeSingularPoint makes its trial steps
    bool                                      m_inBisection;
};


} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsALMBase.hpp)
#endif
