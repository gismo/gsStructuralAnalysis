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
#include <limits>
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

    \ingroup gsALMSolvers

    \par Deferred: matrix-free / iterative continuation

    This class factorizes the tangent stiffness directly (\a factorizeMatrix /
    \a solveSystem, \c SimplicialLDLT by default, \c CG also selectable through the
    \c Solver option) rather than solving it matrix-free. A PyNCT-style architecture
    (secant tangents plus Newton-Krylov with deflated GMRES, where the tangent
    stiffness is never factorized) does not graft onto this class as an option
    toggle:

    -# **The delta.** PyNCT replaces the Newton corrector's linear solve with
       Newton-Krylov (deflated GMRES on secant-approximated tangents); no
       factorization of \f$K_T\f$ is ever formed.
    -# **Why it does not graft on.** This class's singular-point DETECTION signal
       IS the factorization's inertia (\a negatives(), filled from the LDLT
       \c vectorD for \c bifmethod::Determinant), which is free with a direct
       solve and expensive without. The iterative counterparts -- shift-invert
       eigensolves (which themselves want a factorization), or pde2path's
       \c bifcheck=2 with an iterative eigensolver -- mean a matrix-free backend
       does not merely swap the linear solver: it REMOVES the detection signal
       this class relies on and must replace it with something else.
    -# **Revisit triggers.** A target problem whose factorization dominates the
       explore loop (well beyond 1e5-1e6 DoFs, e.g. 3-D solids), or an M-series
       milestone that requires it.
    -# **What it would be.** A design study PLUS a solver backend, not an option
       toggle. Deferred by user decision (2026-08-13) after review of this
       trade-off; this note records that decision, it does not promise a design.

    \par Provisional: the SingularPointTestTol default

    The \c SingularPointTestTol default (\c 1e-4) is calibrated BY MEASUREMENT
    against the 2-DOF unit fixtures of \c unittests/gsALMSolvers_test.cpp
    (\c TEST(singular_point_test_is_scale_invariant)) and the 1-D Bratu
    continuation driver \c example_ModifiedBratuExploration ONLY, run with
    \c OMP_NUM_THREADS=1 and every option set explicitly (never via the
    library default). On the Modified-Bratu thesis benchmark, \c 1e-4 is
    the smallest of the four measured candidates
    \f$\{10^{-6},10^{-5},10^{-4},10^{-3}\}\f$ that classifies the curve-C branch
    point (measured cosine \c 5.293e-05, invariant across \c SingularPointTestIt
    and across the bracketing tolerances) as a BRANCH; \c 1e-3 gives curve C a
    numerically larger margin but is not chosen, because no PDE-scale LIMIT point
    in the measured evidence has an observed cosine against which to bound the
    false-branch risk. Its behaviour at PDE scale -- large 2-D/3-D shell
    or solid problems -- is UNEXAMINED: this is an open question, not a claim
    that the default is either known good or known bad there. A caller running
    such a problem should treat the default as provisional and cross-check any
    limit-vs-branch verdict against an independent signal, e.g. the tangent
    inertia count exposed by \a negatives() (the \c negatives column of the
    landscape CSV; see \ref gsALMSolvers_SingularPoints).

    \par Singular-point option semantics

    \c SingularPointTestTol thresholds a DIMENSIONLESS cosine
    \f$|V\cdot f|/\|f\|\in[0,1]\f$ between the (normalised) approximate null
    vector \f$V\f$ and the forcing \f$f\f$ (see \a _testSingularPoint): because
    the comparison is normalised by \f$\|f\|\f$, the effective threshold does not
    depend on the load magnitude or the unit system, so a caller who rescales the
    forcing does not have to retune the option. Compared with a hypothetical
    un-normalised test on the raw product \f$|V\cdot f|\f$ (whose effective
    tolerance would be \f$\texttt{tol}/\|f\|\f$), \f$\|f\|>1\f$ yields more
    BRANCH verdicts and \f$\|f\|<1\f$ more LIMIT verdicts; at \f$\|f\|=1\f$ the
    two coincide exactly (division by 1 is exact).

    \c SingularPointModeTol is the convergence tolerance for the critical-mode
    inverse power iteration ITSELF, for BOTH the right mode
    (\a _computeCriticalMode) and the left mode (\a _computeCriticalModeLeft);
    see \ref gsALMSolvers_SingularPoints for why it is decoupled from
    \c SingularPointTestTol. Two decades (\c 1e-2) is the smallest margin that
    turns the verdict into a statement about the operator rather than about the
    power iteration. At \c <= 0 (the default, \c -1) it is derived in
    \a getOptions as \c SingularPointTestTol*1e-2, or, when
    \c SingularPointTestTol is itself \c <= 0, as \c Tol -- a literal default
    here would silently stop tracking \c SingularPointTestTol should that
    default's own value ever change. An explicit value at or above
    \c SingularPointTestTol is honoured, not clamped, but warned about only
    when \c SingularPointTestTol itself is positive; at \c SingularPointTestTol
    \c <= 0 the guard is silent.

    \c SingularPointTestIt bounds the number of inverse-power sweeps employed
    to approximate the critical mode, both for the limit-vs-branch test
    (\a _testSingularPoint) and for the mode that seeds the extended system
    (\a _bisectionSolve); see \a _computeCriticalMode. It is an UPPER BOUND, not
    an exact count: the iteration stops as soon as the mode direction changes by
    less than \c SingularPointModeTol between two sweeps (measured up to sign).
    Each sweep is one back-substitution on the already factorized \f$K_T\f$ --
    \a _computeCriticalMode factorizes once, outside the loop -- so a bound of
    \c 20 costs little next to the assembly and factorization that
    precede it, and buys the two extra decades of resolution
    \c SingularPointModeTol demands of the mode. Reaching the bound without
    converging to \c SingularPointModeTol is not by itself fatal; see the three
    bands documented on \a _computeCriticalMode and \a singularPointVerdict.

    \c SingularPointComposite governs what an extended-system solve certifies:
    the set on which \f$\|K_T V\|\f$ vanishes is generally a curve through the
    singular point and therefore bounds nothing about the equilibrium residual
    \f$\|R\|\f$; see \a _extendedSystemSolve. With the option off, a converged
    solve certifies the point as singular but not as an equilibrium; see
    \ref gsALMSolvers_SingularPoints for how a consumer should read that
    distinction.
*/
template <class T>
class gsALMBase
{
protected:

    typedef typename gsStructuralAnalysisOps<T>::ALResidual_t    ALResidual_t;
    typedef typename gsStructuralAnalysisOps<T>::ALForce_t       ALForce_t;
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
        m_djacobian = [this](gsVector<T> const & x, gsVector<T> const & /*dx*/, gsSparseMatrix<T> & m) -> bool
        {
            return m_jacobian(x,m);
        };

        // initialize variables
        m_numIterations = 0;
        // m_initialized/m_stepTaken must be set BEFORE setLength(): the latter reads
        // m_stepTaken to decide whether m_arcLength_prev may still be seeded.
        m_initialized = false;
        m_stepTaken   = false;
        this->defaultOptions();
        this->setLength(1e-2);
        m_converged = false;
        // m_stability/m_stabilityPrev are otherwise written ONLY by _computeStability()
        // (through init(true)/computeStability()/setIndicator()) and by _step(). A caller
        // that uses initialize(false) reaches stabilityChange() -- which reads BOTH --
        // before either has been assigned. Seeding them EQUAL makes that read report "no
        // change" instead of comparing indeterminate values. +1 is the stability()
        // convention for STABLE, the natural default for the undeformed reference state
        // this solver is seeded from.
        m_stability = m_stabilityPrev = 1;

        // initialize errors
        m_basisResidualF = 0.0;
        m_basisResidualU = 0.0;

        // m_Lprev has no in-class initializer, and gsALMCrisfield::initMethods() /
        // gsALMLoadControl::initMethods() do not assign it (only gsALMRiks /
        // gsALMConsistentCrisfield seed it, to 0.0, in their own initMethods()). Left
        // unseeded it is indeterminate storage until the first accepted step or an
        // explicit setPrevious() call, and the public solutionLPrev() accessor would
        // report that indeterminate value. Seeded here, before any derived initMethods()
        // runs, for all four current solvers (and any future subclass).
        m_Lprev = 0.0;

        // m_indicator / m_negatives are the same defect class as m_Lprev above: no
        // in-class initializer, none of the four initMethods() assign them, and both have
        // public accessors (indicator(), negatives()) consumed by gsALMExploration's
        // landscape inertia column. Seeded 0/0 rather than left indeterminate: stability()
        // then agrees with the m_stability = m_stabilityPrev = 1 (STABLE) seed above, and 0
        // is exactly what _computeStability itself produces on the bifmethod::Nothing path
        // -- not an invented sentinel. m_negatives = 0 is the matching inertia ("no
        // negative eigenvalues/pivots yet").
        m_indicator = 0.0;
        m_negatives = 0;

        // m_SPBisProbes is the same defect class as m_indicator/m_negatives above: it has a
        // public accessor (bisectionProbes()) and its sole writer before any solve is
        // _bisectionSolve, so it must not read as indeterminate storage before the first
        // call.
        m_SPBisProbes = 0;

        m_status = gsStatus::NotStarted;
        m_foldTF = std::numeric_limits<T>::quiet_NaN();
        // Same precedent as m_converged above: an outcome flag must not survive a call
        // that did not produce it. Seed NotTested / "no mode measured yet" so the very
        // first classification cannot inherit a stale verdict or a stale (spuriously
        // small) direction error.
        m_SPverdict = SPverdict::NotTested;
        m_SPModeError = std::numeric_limits<T>::max();
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
        // See the Jacobian_t constructor: both flags precede setLength().
        m_initialized = false;
        m_stepTaken   = false;
        this->defaultOptions();
        this->setLength(1e-2);
        m_converged = false;
        // See the Jacobian_t constructor: both stability flags are seeded EQUAL so that a
        // stabilityChange() reached before the first _computeStability()/_step() reports
        // "no change" rather than reading indeterminate values.
        m_stability = m_stabilityPrev = 1;

        // initialize errors
        m_basisResidualF = 0.0;
        m_basisResidualU = 0.0;

        // See the Jacobian_t constructor: same m_Lprev-seeding precedent (uninitialized
        // otherwise for gsALMCrisfield / gsALMLoadControl until the first accepted step).
        m_Lprev = 0.0;

        // See the Jacobian_t constructor: same m_indicator/m_negatives-seeding precedent
        // -- uninitialized otherwise, publicly readable via indicator()/negatives(), and
        // consumed by gsALMExploration's landscape columns.
        m_indicator = 0.0;
        m_negatives = 0;

        // See the Jacobian_t constructor: same m_SPBisProbes-seeding precedent.
        m_SPBisProbes = 0;

        m_status = gsStatus::NotStarted;
        m_foldTF = std::numeric_limits<T>::quiet_NaN();
        // See the Jacobian_t constructor: same m_converged-precedent seeding.
        m_SPverdict = SPverdict::NotTested;
        m_SPModeError = std::numeric_limits<T>::max();
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

    /**
     * @brief      Initialize the arc-length method; computes the stability of the initial
     *             configuration if \a stability is true.
     *
     * With \a stability the initial tangent is ASSEMBLED here, so a `false` from the user's
     * Jacobian reaches this scope as the library-internal `throw 2`, mapped onto
     * \a gsStatus exactly as every other public entry point (\a step,
     * \a computeSingularPoint, \a computeStability) maps it, and reported through the
     * return value (and \a status()).
     *
     * @note \a m_initialized stays true even on a failure: \a initMethods() has already
     *       sized the state, only the (optional) stability assembly failed, and the caller
     *       is told about it through the returned status.
     */
    virtual gsStatus initialize(bool stability = true)
    {
        m_initialized = true;
        try
        {
            this -> initMethods();
            this -> init(stability);
            m_status = gsStatus::Success;
        }
        catch (int errorCode)
        {
            m_status = _statusFromCode(errorCode);
        }
        catch (...)
        {
            m_status = gsStatus::OtherError;
        }
        return m_status;
    }

    /**
     * @brief      Set arc length to \a length
     *
     * \a m_arcLength_prev is deliberately NOT overwritten once a step has been accepted:
     * it is the arc length that PRODUCED the current secant \f$(U-U_{prev},\ \Lambda-\Lambda_{prev})\f$,
     * and the secant predictors of gsALMRiks / gsALMConsistentCrisfield normalise with it
     * (\f$\delta U = (U-U_{prev})/\Delta s_{prev}\f$, then \f$\delta U \mathrel{*}= \Delta s\f$).
     * Assigning both to \a length would make that ratio identically one, so every
     * setLength-based failure-recovery loop ("halve the arc length and retry") would
     * reproduce the previous increment verbatim instead of a reduced one. The invariant is
     * maintained by _step() (and computeLength()) at the end of every accepted step.
     *
     * Before the first accepted step there is no secant yet, so \a m_arcLength_prev is
     * also seeded here, which is what makes setLength() usable directly to configure the
     * INITIAL step length.
     */
    virtual void setLength(T length)
    {
      m_options.setReal("Length",length);
      m_arcLength = m_arcLength_ori = m_options.getReal("Length");
      if (!m_stepTaken)
        m_arcLength_prev = m_arcLength;
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
    /// True if the MOST RECENT solve performed by this object converged.
    ///
    /// "Most recent solve" means the last call to either public solve entry point:
    ///  - \a step(), whose corrector sets the flag from \a iterationFinish() and clears
    ///    it at entry, so a step that hits the iteration limit reports \c false;
    ///  - \a computeSingularPoint(), which publishes the outcome of its extended-system
    ///    stage. On the `testPoint == true` branch that classifies the point as a LIMIT
    ///    point, no solve is attempted at all -- it throws code 1 -- and the flag is set to
    ///    \c false before the throw, so \a converged() and the returned
    ///    \a gsStatus::NotConverged never contradict each other on either path.
    ///
    /// One exception does remain, and it is a THROW path rather than a return path: when
    /// \a _testSingularPoint itself throws (a failing \a computeJacobian or forcing
    /// callback, code 2 -> \a gsStatus::AssemblyError), it does so BEFORE any assignment to
    /// the flag, so \a converged() still describes the previous solve. Nothing is solved on
    /// that path either; the status is the reliable signal there.
    virtual bool converged() const {return m_converged;}

    /// Returns the number of Newton iterations performed
    virtual index_t numIterations() const { return m_numIterations;}

    /// Returns the number of \a step() calls spent by the MOST RECENT \a _bisectionSolve.
    /// \c 0 is AMBIGUOUS: it is returned both when that call found no bracket to localize (or
    /// has not run yet) AND when a bracket WAS found but \c SingularPointBisIt was configured
    /// \c &lt;= 0, so the probe loop never executes a single iteration. Neither case can be
    /// told apart from this accessor alone; a caller that needs to distinguish "no bracket"
    /// from "budget configured away" must check \c SingularPointBisIt itself.
    virtual index_t bisectionProbes() const {return m_SPBisProbes;}

    /// Returns the tolerance value used
    virtual T tolerance() const {return m_tolerance;}

    /// Returns the error after solving the nonlinear system
    virtual T residue()   const {return m_residue;}

    /// Returns the value of the Determinant or other indicator. Explicitly seeded
    /// \c 0.0 in both \a gsALMBase constructors, so this returns a
    /// well-defined value from construction onwards on every subclass, not only
    /// after the first \a _computeStability() call has run.
    virtual T indicator() const {return m_indicator;}
    virtual void setIndicator(T indicator) {m_indicator = indicator; m_stability = this->stability();}

    /// Returns the number of negative entries in the stability vector, i.e. the
    /// inertia (count of negative eigenvalues or, for the Determinant method,
    /// negative pivots) of the tangent stiffness at the current point. Filled by
    /// _computeStability for BOTH the Determinant and Eigenvalue methods. A change
    /// of this count between accepted points signals a singular point, including
    /// SECONDARY bifurcations on already-unstable branches that a min-eigenvalue
    /// sign flip cannot see.
    /// \note For the Determinant method the count comes from the unpivoted LDLT
    /// vectorD of a possibly indefinite matrix and can be unreliable; the
    /// Eigenvalue method yields exact inertia counts.
    /// \sa indicator(). Explicitly seeded \c 0 in both \a gsALMBase constructors,
    /// so this returns a well-defined value from construction onwards
    /// on every subclass, not only after the first \a _computeStability() call.
    virtual index_t negatives() const {return m_negatives;}

    /**
     * @brief AUTO-07p's FOLD test function at the state whose tangent was last
     *        assembled AND factorized by \a _computeStability (i.e. by a
     *        \a computeStability(true) / \a init(true) call).
     *
     * Reference: auto-07p \c toolboxae.f90:703-747 (\c FNLPAE) -- the last component of
     * the NORMALISED solution of the bordered system
     * \f{align*}{
     *   \begin{bmatrix} K_T & -f \\ w^T & \gamma \end{bmatrix}
     *   \begin{bmatrix} z \\ \zeta \end{bmatrix} =
     *   \begin{bmatrix} 0 \\ 1 \end{bmatrix} .
     * \f}
     * Eliminating the top block gives \f$z = \zeta\,z_t\f$ with
     * \f$z_t = K_T^{-1} f\f$ and \f$\zeta = 1/(w\cdot z_t + \gamma)\f$, so after
     * normalisation the last component is
     * \f[ \mathrm{TF} = \frac{\mathrm{sign}(w\cdot z_t+\gamma)}{\sqrt{1+\|z_t\|^2}} . \f]
     * The MAGNITUDE does not depend on the bordering \f$(w,\gamma)\f$ at all; only the
     * sign does. We take \f$(w,\gamma)=(f,0)\f$, so
     * \f$\mathrm{TF}=\mathrm{sign}(f\cdot z_t)/\sqrt{1+\|z_t\|^2}\f$ with
     * \f$f\cdot z_t = f^T K_T^{-1} f\f$ the compliance in the load direction.
     * (AUTO borders with the arc-length constraint row; that row is method-specific
     * here and needs step state -- \a m_DeltaU is zero or stale when
     * \a computeStability(true) is called standalone -- so \f$(f,0)\f$ was chosen
     * instead: it needs no step state, and \f$\mathrm{sign}(f^TK_T^{-1}f)\f$ IS the
     * fold-vs-branch discriminator. The magnitude is unaffected by the choice.)
     *
     * At a simple LIMIT point the critical mode satisfies \f$\phi\cdot f\neq0\f$, so
     * \f$f^TK_T^{-1}f\sim(\phi\cdot f)^2/\mu\f$ changes sign with the critical
     * eigenvalue \f$\mu\f$ and TF crosses zero: **a sign change between two accepted
     * points brackets a fold.** At a BRANCH point \f$\phi\cdot f\approx0\f$ and it does
     * not change sign -- which is exactly what makes AUTO's LP function a FOLD test and
     * not a general singularity test.
     *
     * @note INSTRUMENT ONLY. Nothing in this library reads it: singular-point detection
     *       is the inertia count \a negatives(), and the limit-vs-branch verdict is
     *       \a _testSingularPoint. It costs ONE back-substitution on a factorization
     *       that has just been formed with a DIRECT solver (the default,
     *       \c SimplicialLDLT); with an ITERATIVE solver (\c CG, see the \c solver::type
     *       enum) it is a FULL solve per \c jacobian=true stability evaluation.
     * @note Returns NaN when no value is available: before the first stability
     *       evaluation, after a \a _computeStability(x,false) call (which cannot own the
     *       factorization it would need), and when the back-substitution itself fails.
     *       Check with \c math::isnan before using it.
     */
    virtual T foldTestFunction() const {return m_foldTF;}

    /// Return the solution vector and factor
    virtual const gsVector<T> & solutionU() const {return  m_U;}
    virtual const gsVector<T> & solutionDU() const {return  m_DeltaU;}
    virtual T  solutionL() const {return  m_L;}
    virtual T  solutionDL() const {return  m_DeltaL;}
    virtual const gsVector<T> & solutionV() const {return  m_V;}

    /// The verdict of the most recent limit-vs-branch classification (\a _testSingularPoint);
    /// see \a SPverdict.
    virtual index_t singularPointVerdict() const {return m_SPverdict;}
    /// Direction error |dV| achieved by the critical-mode iteration of that classification
    /// (max over the right and left modes). Compare against \c SingularPointModeTol /
    /// \c SingularPointTestTol. \c SingularPointTestTol's default is provisional at PDE
    /// scale -- see the class-level "Provisional: the SingularPointTestTol default" note.
    virtual T criticalModeError() const {return m_SPModeError;}
    /// The LEFT critical mode psi of the current tangent (equal to \a solutionV() whenever
    /// the configured solver is self-adjoint); see \a _computeCriticalModeLeft.
    /// @warning Valid only IMMEDIATELY AFTER a classification. \a _extendedSystemSolve
    ///          refines the RIGHT mode in place (m_V += m_DeltaV) and does NOT refine this
    ///          one, so after a converged computeSingularPoint() the two are inconsistent
    ///          even under symmetry. Re-syncing them is not implemented.
    virtual const gsVector<T> & solutionVLeft() const {return m_Vleft;}

    /// Returns true if the current solution point is STABLE, i.e. if stability()
    /// is +1. (Whether a bifurcation point was passed is stabilityChange(), not
    /// this function.) m_stability is an index_t holding +1 or -1, so this MUST be
    /// a sign test: returning it directly would make isStable() always true.
    virtual bool isStable() const {return m_stability > 0;}

    /// Diagnostic accessor (read-only): the STORED stability sign of the
    /// PREVIOUS accepted point, i.e. the second operand of \a stabilityChange(). Exposed
    /// so a caller can observe whether it survives a reseed (\a setSolution / \a
    /// setPrevious / \a setIndicator / \a setLength do not write it -- see \a
    /// gsALMExploration::traceSweep's :538-544 claim) without inferring it from source.
    virtual index_t stabilityPrev() const {return m_stabilityPrev;}

    /// Diagnostic accessor (read-only): the previously converged point injected
    /// by \a setPrevious (or the construction-time rest state). Exposed for the same
    /// reason as \a stabilityPrev(). \a m_Uprev is default-constructed (size 0, so
    /// \c .norm()==0) from construction until the first accepted step or \a setPrevious().
    virtual const gsVector<T> & solutionUPrev() const {return m_Uprev;}
    /// \sa solutionUPrev(). \a m_Lprev is explicitly seeded \c 0.0 in both \a gsALMBase
    /// constructors, so this returns a well-defined \c 0.0 from construction
    /// onwards on every subclass, not only the two (\a gsALMRiks, \a gsALMConsistentCrisfield)
    /// whose \a initMethods() re-seeds it.
    virtual T solutionLPrev() const {return m_Lprev;}

    /// Diagnostic accessor (read-only): the arc length that produced the CURRENT
    /// secant (see \a m_arcLength_prev's own doc for the invariant and the ORDER-
    /// DEPENDENT reseed rule in \a setPrevious()).
    virtual T getLengthPrev() const {return m_arcLength_prev;}

    /// Diagnostic accessor (read-only): true once at least one step has been
    /// ACCEPTED (see \a m_stepTaken's own doc -- it is deliberately never cleared by
    /// \a setLength once true, which gates whether \a m_arcLength_prev is reseeded).
    virtual bool stepTaken() const {return m_stepTaken;}

    // /// Returns the value of the deterimant of the jacobian
    // virtual T determinant() const
    // {
    //     return m_jacobian(m_U).toDense().determinant();
    // }

    /// Resets the step
    virtual void resetStep() {m_DeltaUold.setZero(); m_DeltaLold = 0;}

    // Set initial guess for solution
    virtual void setInitialGuess(const gsVector<T> & Uguess, const T & Lguess) {m_Uguess = Uguess; m_Lguess = Lguess;}
    /**
     * @brief      Injects the previously converged point, i.e. the secant
     *             \f$(U-U_{prev},\ \Lambda-\Lambda_{prev})\f$ the predictors extrapolate along.
     *
     * The injected secant was produced by an arc length this object never saw - \a setPrevious
     * carries none - so \a m_arcLength_prev is (re)declared to be the CURRENT arc length. That
     * is what every caller means, and it is what the pre-\a setLength-fix code effectively did
     * (it assigned m_arcLength_prev = m_arcLength on every setLength). Without this, a REUSED
     * solver object would normalise an externally injected secant by the step length of an
     * unrelated earlier interval - see gsAPALM, which drives every job through one m_ALM.
     *
     * @note The rule is ORDER-DEPENDENT: it yields the intended ratio one only when
     *       \a setLength() is called BEFORE \a setPrevious() - which is what gsAPALM does at
     *       both of its seeding sites (the setLength/setSolution/setPrevious triples at the
     *       head of \a gsAPALM::_initiation and of \a gsAPALM::_correction).
     *       Where the order is reversed - the \a setPrevious ... \a setLength
     *       pair in \a gsALMExploration::traceSweep(), which seeds the solver once
     *       per arc-length SWEEP - the divisor keeps whatever
     *       \a m_arcLength held at that moment; that is harmless there because the same call
     *       site sets
     *       \f$U_{prev} = U\f$, which selects the predictors' fresh-start branch (the divisor
     *       is not read), and \a _step() re-establishes it at the first accepted step.
     * @note \a m_stepTaken is deliberately NOT cleared: clearing it would let the NEXT
     *       \a setLength() re-seed \a m_arcLength_prev and so neuter the failure-recovery fix.
     */
    virtual void setPrevious(const gsVector<T> & Uprev, const T & Lprev)
    {
        m_Uprev = Uprev;
        m_Lprev = Lprev;
        m_DeltaUold = m_U - m_Uprev;
        m_DeltaLold = m_L - m_Lprev;
        m_arcLength_prev = m_arcLength;
    }

    /// Sets the solution
    virtual void setSolution(const gsVector<T> & U, const T & L) {m_L = L; m_U = U; }// m_DeltaUold.setZero(); m_DeltaLold = 0;}

    /// Sets the solution step
    virtual void setSolutionStep(const gsVector<T> & DU, const T & DL) {m_DeltaUold = DU; m_DeltaLold = DL;}// m_DeltaUold.setZero(); m_DeltaLold = 0;}

    /**
     * @brief      Sets an (OPT-IN) state-dependent forcing function
     *             \f$ f(U,\Lambda) = -\partial R/\partial \Lambda \f$.
     *
     * By default the arc-length methods assume a **dead load**: the residual is
     * \f$ R(U,\Lambda) = F_{int}(U) - \Lambda F_{ext}\f$ so that
     * \f$ -\partial R/\partial\Lambda = F_{ext}\f$ is the *constant* vector \a Force
     * handed to the constructor. Setting this callback tells the solver that the load
     * derivative depends on the state, and the correctors then evaluate it at the
     * current iterate instead of using the stored constant.
     *
     * ### Theory
     * An arc-length method is Newton's method on the bordered system
     * \f{align*}{
     *    R(U,\Lambda) &= 0, \\
     *    c(\Delta U,\Delta\Lambda) &= 0,
     * \f}
     * where \f$c\f$ is the arc-length constraint. Block elimination of the Newton step
     * gives, with \f$K_T = \partial R/\partial U\f$ and \f$R_\Lambda = \partial R/\partial\Lambda\f$,
     * \f{align*}{
     *    \delta U &= \delta\bar{u} + \delta\Lambda\,\delta u_t, \\
     *    \delta\bar{u} &= -K_T^{-1} R, \\
     *    \delta u_t &= -K_T^{-1} R_\Lambda \;=\; K_T^{-1} f(U,\Lambda), \qquad f \equiv -\partial R/\partial\Lambda .
     * \f}
     * The dead-load code path uses \f$f = \f$ \a m_forcing \f$ = \f$ const; with this
     * callback the *same* formulas are used but \f$f\f$ is evaluated at the CURRENT
     * iterate. The \f$\delta\Lambda\f$ constraint updates (Riks' hypersphere
     * linearization, Crisfield's quadratic root solve) are unchanged.
     *
     * ### Why it matters
     * If the true load derivative is \f$\kappa f\f$ while the solver uses \f$f\f$
     * (\f$\kappa\neq 1\f$), every corrector iteration injects a residual
     * \f$\approx\delta\Lambda(\kappa-1)f\f$. The locally *quadratic* bordered Newton
     * degenerates into a linear fixed-point iteration with multiplier \f$1-g\f$,
     * \f[
     *    g = \frac{(1-\varphi) + \varphi\,\kappa^2 u}{(1-\varphi) + \varphi\,\kappa\, u},
     *    \qquad u = \|K_T^{-1}f\|^2 ,
     * \f]
     * which **diverges** (sign-alternating, period-2) whenever
     * \f$\varphi\,u\,\kappa(\kappa-2) > 1-\varphi\f$. That criterion contains no arc
     * length, so step-size reduction can never repair it. For \f$\kappa=1\f$ (dead
     * load) \f$g=1\f$ and the multiplier is \f$0\f$ — which is why the dead-load path
     * is robust even on unstable branches and through limit points. Supplying the
     * consistent \f$f(U,\Lambda)\f$ restores quadratic convergence and exact fold
     * rounding.
     *
     * ### Contract
     * - \a fun must return \f$-\partial R/\partial\Lambda\f$ at the given state,
     *   consistently with the \a ALResidual passed to the constructor. It returns
     *   \c false on an assembly failure (same convention as the residual and Jacobian
     *   functors); the solver then reports gsStatus::AssemblyError.
     * - The Jacobian functor must supply the **consistent tangent**, i.e. it must
     *   already include the load-stiffness contribution \f$-\Lambda\,\partial b/\partial U\f$
     *   for a load \f$\Lambda b(U)\f$. This is the pre-existing Jacobian contract, but
     *   it only becomes load-bearing here.
     * - Load stiffness is generally **nonsymmetric** (e.g. follower pressure). In that
     *   case select a nonsymmetric solver and the eigenvalue-based bifurcation test:
     *   \c options().setString("Solver","LU") and
     *   \c options().setInt("BifurcationMethod",bifmethod::Eigenvalue), since the
     *   default SimplicialLDLT/Determinant path assumes symmetry.
     *
     * ### Where \f$f\f$ is (re)evaluated
     * - corrector tangent \f$\delta u_t\f$ (computeUt, and the extended-system solve):
     *   at the state the Jacobian of that solve was evaluated at,
     *   \f$(U+\Delta U,\ \Lambda+\Delta\Lambda)\f$;
     * - predictors: at the predictor linearization state \f$(U,\Lambda)\f$;
     * - residual-norm scaling (computeResidualNorms): at
     *   \f$(U+\Delta U,\ \Lambda+\Delta\Lambda)\f$;
     * - singular-point classification (\f$|V\cdot f|\f$): at the tested state
     *   \f$(U,\Lambda)\f$ — with the consistent \f$f\f$ this test becomes exact, since a
     *   branch point is characterised by \f$f\in\mathrm{range}(K_T)\f$;
     * - constraint metric scalings of Crisfield's quadratic (the \f$\|f\|^2\f$ terms and
     *   the scaling \f$\varphi\f$) are **frozen per step** at the predictor state, see
     *   \a stepForcing. This is a definitional choice of the constraint surface, not an
     *   approximation: it keeps \f$c\f$ a fixed quadratic during the corrector, which is
     *   what makes the \f$\delta\Lambda\f$ root solve well posed.
     *
     * @warning The **extended system** (\a _extendedSystemIteration) is exact Newton only
     *          for a dead load. Its finite-difference block differentiates \f$K_T\f$ in the
     *          \f$V\f$ direction alone, whereas a state-dependent load makes \f$K_T\f$
     *          depend on \f$\Lambda\f$ as well (through the load stiffness
     *          \f$-\Lambda\,\partial b/\partial U\f$), so the \f$\partial(K_TV)/\partial
     *          \Lambda\f$ block is missing. With a callback set, the extended solve is
     *          therefore an inexact Newton and may converge more slowly or not at all;
     *          the singular-point *classification* above is unaffected. Dead-load
     *          behaviour (callback unset) is exactly as before.
     *
     * @note Default (callback unset) behaviour is bit-identical to the dead-load code:
     *       all helpers then return a reference to the stored \a m_forcing itself.
     *       When the callback IS set, it costs one extra evaluation per iteration (an
     *       assembly of a load vector — the same cost class as the residual); the result
     *       is deliberately NOT cached across states.
     *
     * @param[in]  fun   The load-derivative functor \f$(U,\Lambda)\mapsto -\partial R/\partial\Lambda\f$
     */
    virtual void setForcingFunction(const ALForce_t & fun)
    {
        m_forcingFun = fun;
        // Keep the frozen step forcing sized (and equal to the dead load) until the
        // first initiateStep() refreshes it.
        m_forcingStep = m_forcing;
    }

    /// Access the options
    virtual gsOptionList & options() {return m_options;};

    /// Set the options to \a options
    virtual void setOptions(gsOptionList options) {m_options.update(options,gsOptionList::addIfUnknown); this->getOptions(); };

    /// Return the options into \a options
    virtual void options_into(gsOptionList options) {options = m_options;};

    /// Apply the options
    virtual void applyOptions() {this->getOptions(); }

    virtual T distance(const gsVector<T>& /*DeltaU*/, const T /*DeltaL*/) const
    {
        GISMO_NO_IMPLEMENTATION;
    }

// ------------------------------------------------------------------------------------------------------------
// ---------------------------------------Singular point methods-----------------------------------------------
// ------------------------------------------------------------------------------------------------------------
public:


    /**
     * @brief      Computes the singular point seeded from \a (U,L).
     *
     * The solver state is set to \a (U,L) FIRST, so nothing this routine evaluates is ever
     * evaluated at whatever state the solver happened to sit at on entry. From there the
     * pipeline is localize, then classify, then refine: when
     * \a SingularPointComputeTolB != 0, the bisection stage localizes the crossing
     * starting from the just-installed \a (U,L); the limit-vs-branch classification
     * (\a _testSingularPoint, only when \a testPoint is true), the tangent it uses (only
     * re-assembled when \a jacobian is true; otherwise the tangent the localization stage
     * left, or the caller's \a m_jacMat unchanged if that stage did not run), and
     * \a computeForcing all then act on that localized point; the extended system solve
     * refines it last. The arguments are the point the caller asks ABOUT, not merely a
     * starting guess appended after the fact.
     *
     * Return contract:
     *  - a converged branch-point refinement returns \a gsStatus::Success and leaves the
     *    solver at the refined singular point;
     *  - a non-converged extended solve returns \a gsStatus::NotConverged;
     *  - \a testPoint == true and the point classifying as a LIMIT point also returns
     *    \a gsStatus::NotConverged: the bisection stage (when
     *    \a SingularPointComputeTolB != 0) may already have run, so the solver is left at
     *    the localized point rather than necessarily at \a (U,L); no extended solve is
     *    attempted; and an unconditional message on \a gsInfo names the limit-point cause.
     *    This is the only way to tell that case apart from the non-converged one, since no
     *    dedicated \a gsStatus enumerator exists (the enum is shared with the dynamic and
     *    APALM solvers). \a Success therefore never coincides with "nothing happened".
     *  - an assembly / linear-solver failure returns \a gsStatus::AssemblyError /
     *    \a gsStatus::SolverError (see \a _statusFromCode).
     *
     * Throw semantics: if the classification throws (a failing \a computeJacobian or
     * forcing callback), the solver is left at the point the pipeline had reached when the
     * throw occurred -- the localized point when the bisection stage ran
     * (\a SingularPointComputeTolB != 0), \a (U,L) otherwise -- rather than at its state on
     * entry.
     *
     * @param[in]  U             The point to compute the singular point at (displacements)
     * @param[in]  L             The point to compute the singular point at (load factor)
     * @param[in]  switchBranch  Switch onto the bifurcated branch after a successful solve
     * @param[in]  jacobian      Re-assemble the tangent at the point being classified? If
     *                           false, the tangent left in the solver by the caller is
     *                           reused as-is only when the bisection stage did not run
     *                           (\a SingularPointComputeTolB == 0); when that stage did
     *                           run, the classification instead sees the tangent it last
     *                           factorized.
     * @param[in]  testPoint     Classify the point first? If false, the point is ASSUMED to
     *                           be a branch point and the solve stages run unconditionally
     *                           (so neither the classification nor the limit-point return
     *                           contract above applies).
     *
     * @return     The status of the computation, see above.
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
    ///
    /// FALLBACK branch-switch predictor: nudges the singular point along the
    /// unit critical mode, \f$U \leftarrow U^*+V/\tau\f$ (option \c Perturbation).
    /// Used when \a computeBranchTangent does not return
    /// \a branchTangent::Success or \a branchTangent::TrivialBranch -- i.e.
    /// whenever the true swibra ABE emanating tangent cannot be resolved.
    virtual void switchBranch();

    /// Outcome of \a computeBranchTangent.
    struct branchTangent
    {
        enum type
        {
            Success        = 0, ///< tau1 is the emanating tangent
            TrivialBranch  = 1, ///< K_T does not vary along phi1 (linear problem): tau1 = [phi1;0]
            ModeUnresolved = 2, ///< no usable left mode, or |phi1.psi1| below BranchTangentTol
            Degenerate     = 3, ///< alpha1bar ~ 0 ("no distinct branch"), or lambdadot0 ~ 0
            AssemblyFailed = 4  ///< the Jacobian or the forcing callback failed
        };
    };

    /**
     * @brief      Computes the TRUE emanating branch tangent \f$\tau_1\f$ at a
     *             simple (multiplicity-1) branch point, by the algebraic
     *             branching equation (ABE) of pde2path's \c swibra.
     *
     * Sign convention: pde2path's residual \c G is our \c R; \c Gu is our
     * \f$K_T\f$; \c Glam \f$=\partial R/\partial\Lambda = -f\f$, where \f$f\f$
     * is what \a computeForcing returns.
     *
     * @param[in]  Ustar      The branch point (displacements)
     * @param[in]  Lstar      The branch point (load factor)
     * @param[in]  phi1       The (right) critical mode at \a (Ustar,Lstar), e.g.
     *                        \a solutionV(). Passed explicitly rather than read
     *                        from \c m_V, so the contract carries no hidden
     *                        ordering requirement.
     * @param[in]  tangentU0  The incoming (pre-crossing) tangent/secant,
     *                        displacement part
     * @param[in]  tangentL0  The incoming (pre-crossing) tangent/secant, load part
     * @param[out] tau1U      The emanating branch tangent, displacement part
     * @param[out] tau1L      The emanating branch tangent, load part
     *
     * @return     The outcome, see \a branchTangent.
     *
     * @note State/purity contract: leaves \c m_U, \c m_L, \c m_V, \c m_jacMat
     *       and \c m_Vleft exactly as it found them (member values restored
     *       unconditionally). The solver's cached FACTORIZATION is untouched on
     *       the default self-adjoint path (the left mode there is a pure copy,
     *       no factorization is ever performed) and best-effort restored
     *       otherwise -- see the implementation note on why an unconditional
     *       restore is neither possible nor correct at a branch point, where
     *       \f$K_T(U^*)\f$ is singular by construction. Performs two Jacobian
     *       ASSEMBLIES (never through \a computeJacobian, which factorizes and
     *       would abort on that singular tangent) and, at most, whatever the
     *       LEFT critical-mode routine costs -- no solve of its own.
     */
    virtual typename branchTangent::type
    computeBranchTangent(const gsVector<T> & Ustar, const T & Lstar,
                         const gsVector<T> & phi1,
                         const gsVector<T> & tangentU0, const T & tangentL0,
                         gsVector<T> & tau1U, T & tau1L);

    /// Reduce the length by multiplication with a factor \a fac
    virtual T reduceLength(T fac = 0.5);

    /// Reset the length
    virtual T resetLength();

protected:

    /**
     * @brief      Maps the library-internal integer error codes onto \a gsStatus.
     *
     * The internal routines signal failure by `throw <int>` (1: not converged,
     * 2: assembly failure, 3: linear-solver failure). Every PUBLIC entry point that can
     * reach such a throw catches it and reports a \a gsStatus instead; this single
     * definition is what keeps those entry points from drifting apart.
     */
    static gsStatus _statusFromCode(int errorCode)
    {
        if      (errorCode==1) return gsStatus::NotConverged;
        else if (errorCode==2) return gsStatus::AssemblyError;
        else if (errorCode==3) return gsStatus::SolverError;
        else                   return gsStatus::OtherError;
    }

    /// See \a computeStability
    virtual void _computeStability(const gsVector<T> & x, bool jacobian=true, T shift = -1e2);

    /**
     * @brief      See \a computeSingularPoint, which wraps this and maps the int codes.
     *
     * Sets \a m_U = U and \a m_L = L BEFORE anything else, so that the localization stage,
     * the classification and the solve stages all act on the requested point (or on the
     * point the localization stage moved it to); may therefore leave the solver at the
     * localized point rather than at \a (U,L) when it throws. Signals a limit point
     * (nothing computed) and a non-converged extended solve alike with \a throw 1,
     * distinguished only by the message each one prints.
     */
    virtual void _computeSingularPoint(const gsVector<T> & U, const T & L, bool switchBranch=false, bool jacobian=false, bool testPoint=true);

    /**
     * @brief      Tests if a point is a bifurcation point
     *
     * @param[in]  jacobian  Evaluate the Jacobian?
     *
     * @return     True if it is a bifurcation point
     */
    virtual bool _testSingularPoint(bool jacobian=false);

    /**
     * @brief      Approximates the CRITICAL MODE of the current tangent \a m_jacMat into
     *             \a m_V by inverse power iteration.
     *
     * The critical mode is the eigenvector belonging to the smallest-magnitude eigenvalue
     * of \f$K_T\f$; it spans the null space at a singular point and is what the
     * limit-vs-branch test (\a _testSingularPoint), the extended system
     * (\a _extendedSystemSolve) and the branch switch (\a switchBranch) all consume.
     *
     * The iteration starts from a DETERMINISTIC pseudo-random vector (fixed-seed LCG) and
     * NOT from \f$\mathbf{1}\f$: a constant start is a symmetric field and is exactly
     * orthogonal to the antisymmetric critical mode of the canonical symmetric-structure
     * bifurcation, which no number of sweeps can then recover. Being a pure function of
     * \f$(K_T,\ \texttt{maxIt},\ \texttt{tol})\f$, it introduces no run-to-run or
     * history-dependent variability.
     *
     * \a maxIt is a MAXIMUM: the sweep stops as soon as the mode direction changes by less
     * than \a tol between two sweeps, measured UP TO SIGN (the critical eigenvalue is
     * negative past a crossing, where the iterate alternates sign).
     *
     * The achieved direction error is published through \a m_SPModeError (hence
     * \a criticalModeError()) for every call, whether or not it converged; \a
     * _testSingularPoint then reads it against three bands, from the tightest to the
     * loosest:
     *   -# **resolved**   (\c |dV| \< \c SingularPointModeTol): classify normally;
     *   -# **imprecise**  (\c SingularPointModeTol \f$\le\f$ \c |dV| \< \c
     *      SingularPointTestTol): classify normally, silently -- a direction error still
     *      below the verdict threshold cannot change the limit-vs-branch verdict;
     *   -# **unresolved** (\c |dV| \f$\ge\f$ \c SingularPointTestTol, when \c
     *      SingularPointTestTol \> 0): the classification REFUSES (\a SPverdict::Unresolved)
     *      instead of certifying a verdict evaluated at or below the resolution of the mode
     *      that feeds it -- see \a singularPointVerdict().
     *
     * @param[in]  maxIt  Maximum number of sweeps (option \c SingularPointTestIt)
     * @param[in]  tol    Convergence tolerance on the mode direction (option \c
     *                    SingularPointModeTol, NOT \c Tol -- decoupled from the verdict
     *                    threshold \c SingularPointTestTol so the classification is never
     *                    evaluated below the mode's own resolution)
     *
     * @return     True if the mode converged within \a maxIt sweeps to \a tol. A \c false
     *             is warned about here only when the achieved error also exceeds
     *             \c SingularPointTestTol (see the bands above); the caller decides whether
     *             that is fatal.
     */
    virtual bool _computeCriticalMode(index_t maxIt, T tol);

    /**
     * @brief      Approximates the LEFT critical mode \f$\psi\f$ of the current tangent
     *             \f$K_T\f$ (i.e. the right null vector of \f$K_T^{\mathsf T}\f$) into
     *             \a m_Vleft, by inverse power iteration on the transpose.
     *
     * Required precondition: \a _computeCriticalMode has already run for THIS tangent (the
     * only caller, \a _testSingularPoint, calls the two back to back).
     *
     * Self-adjoint solver (\a _solverIsSelfAdjoint): \f$K_T^{\mathsf T}=K_T\f$, so the left
     * mode IS the right mode. This branch is \c m_Vleft = m_V; and returns the flag the
     * right-mode iteration achieved, with NO call to \a factorizeMatrix and NO call to
     * \a solveSystem -- the identity and the zero extra cost are then true by construction.
     * It leaves \a m_SPModeError untouched: the right-mode value already describes both
     * vectors.
     *
     * Non-self-adjoint (or unrecognised) solver: the transpose path is always correct but
     * costs one extra factorization plus at most \a maxIt back-substitutions. \a m_jacMat is
     * materialized into a real transposed matrix (its \c .transpose() is an expression
     * template, not storage), factorized, swept with the same fixed-seed LCG start and
     * sign-blind stopping rule as the right mode, and \a m_SPModeError is raised to the max
     * of its previous value and this sweep's error. \a m_solver is then RESTORED to hold
     * \a m_jacMat's factors before returning -- callers of \a _testSingularPoint
     * (\a _bisectionSolve, \a _extendedSystemSolve, the correctors) assume that.
     *
     * @param[in]  maxIt  Maximum number of sweeps (option \c SingularPointTestIt)
     * @param[in]  tol    Convergence tolerance on the mode direction (option \c
     *                    SingularPointModeTol)
     *
     * @return     True if the left mode converged within \a maxIt sweeps to \a tol.
     */
    virtual bool _computeCriticalModeLeft(index_t maxIt, T tol);

    /// One inverse power SWEEP LOOP (up to \a maxIt sweeps) on the ALREADY factorized
    /// operator held by \a m_solver, applied to \a V in place -- shared by
    /// \a _computeCriticalMode and the non-self-adjoint path of
    /// \a _computeCriticalModeLeft so both use ONE copy of the sign-blind stopping rule.
    /// \a V must already hold the (deterministic fixed-seed LCG) start vector on entry.
    /// Writes the achieved direction error into \a dist (set to
    /// \c std::numeric_limits<T>::max() if the solve annihilates or blows up the iterate,
    /// so an unusable mode is never reported as merely imprecise) and returns true iff it
    /// fell below \a tol.
    virtual bool _inversePowerSweeps(gsVector<T> & V, index_t maxIt, T tol, T & dist);

    /// True when the CONFIGURED solver is self-adjoint, i.e. when the caller has asserted
    /// that K_T is symmetric (SimplicialLDLT/LLT are the same type; CG/Pardiso-LDLT/LLT are
    /// self-adjoint by construction). An unrecognised solver returns false: the transpose
    /// path in \a _computeCriticalModeLeft is always correct, the free path is only correct
    /// under symmetry.
    virtual bool _solverIsSelfAdjoint() const;

    /// Perform an extended system iteration
    virtual void _extendedSystemIteration();

    /// Returns the objective function for the bisection method given solution \a x
    virtual index_t _bisectionObjectiveFunction(const gsVector<T> & x, bool jacobian=true);
    /// Returns the termination function for the bisection method given solution \a x
    virtual T       _bisectionTerminationFunction(const gsVector<T> & x, bool jacobian=true);

    /// Perform an extended system solve to find a singular point
    virtual bool _extendedSystemSolve(const gsVector<T> & U, const T L, const T tol);

    /// Localizes a singular point by bisection on the arc-length parameter \c s between two
    /// bracket endpoints: \c (U,L) is the NEAR endpoint (\c s = 0) and \c (Ufar,Lfar) the FAR
    /// one (\c s = dLb), where \c dLb = \a m_arcLength at entry. The routine evaluates the
    /// inertia (\a _bisectionObjectiveFunction) at BOTH endpoints itself; when the two counts
    /// are equal (or \c Ufar carries no usable incumbent, i.e. \c Ufar.size() != U.size())
    /// there is no bracket to bisect and it returns \c false WITHOUT taking a single \a step()
    /// call, leaving \c (m_U,m_L) at the arguments. Otherwise it bisects, maintaining the
    /// invariant that the two endpoints' inertia counts differ, in the shape of
    /// \a gsALMExploration::_localizeCrossing.
    ///
    /// The probe budget is the option \c SingularPointBisIt (\a m_SPBisIt), independent of
    /// \c MaxIter; \a bisectionProbes() reports how many \a step() calls the most recent call
    /// spent (\c 0 on the no-bracket path).
    ///
    /// On every NON-CONVERGED exit \c (m_U,m_L) are restored to the arguments \c (U,L), and
    /// \c m_Uprev/m_Lprev/m_arcLength/m_arcLength_prev/m_adaptiveLength are restored to their
    /// values on entry. On a CONVERGED exit \c (m_U,m_L) is left at the final PROBE -- the
    /// point the indicator-ratio termination test was actually evaluated at and cleared
    /// \c tol on, NOT a reconstructed bracket endpoint -- so containment (\c s &lt;= dLb) holds
    /// up to corrector slack rather than exactly. This is the whole point of the routine for
    /// its \c snapping_* callers.
    ///
    /// Postcondition on EVERY return path, convergence or not: \a m_V holds the critical mode
    /// (\a _computeCriticalMode) at the point left in \c (m_U,m_L) -- this is what makes the
    /// no-bracket early return safe for a caller (\a _extendedSystemSolve) that seeds itself
    /// from \c m_V without recomputing it.
    ///
    /// \note \c dLb is inferred as \c m_arcLength at entry, not passed by the caller (unlike
    /// \a gsALMExploration::_localizeCrossing, which owns the sweep length). A \c dLb too SMALL
    /// silently breaks the contract, because \c s = dLb then does not reach the incumbent and
    /// the far endpoint's inertia is attributed to a position no probe can reach. A \c dLb too
    /// LARGE also breaks containment: the first probe sits at \c s = dLb/2, already past the
    /// incumbent, and when the indicator tolerance is loose enough to certify convergence off
    /// that single probe the loop never bisects back -- measured slack \c 0.0708 at an
    /// inflation factor of 4 and \c 0.455 at 16, on a bracket of width \c 0.0706 (a factor of 2
    /// was not itself observed to break containment). \a m_arcLength can be the too-large one
    /// under \c AdaptiveLength=true: \a computeLength() rescales it by a factor clamped to
    /// [0.5, 2.0] after every accepted step, so the entry \c dLb can be up to 2x the length
    /// that actually produced the incumbent. \a m_arcLength_prev is the length that produced
    /// the bracket -- i.e. the correct \c dLb -- but is not used here.
    virtual bool _bisectionSolve(const gsVector<T> & U, const T L, const T tol,
                                  const gsVector<T> & Ufar, const T Lfar);

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

    // On WHY these are direct factorizations and not a Krylov backend, see the class doc,
    // "Deferred: matrix-free / iterative continuation".
    /// Factorize the matrix \a M
    virtual void factorizeMatrix(const gsSparseMatrix<T> & M);

    /// Solve the system with right-hand side \a F
    virtual gsVector<T> solveSystem(const gsVector<T> & F);

    /// Compute the residual
    virtual gsVector<T> computeResidual(const gsVector<T> & U, const T & L);
    virtual void computeResidual();

    /**
     * @brief      Returns the forcing \f$f = -\partial R/\partial\Lambda\f$ at the state \f$(U,L)\f$.
     *
     * Without a forcing callback (see \a setForcingFunction) this returns a reference to
     * the stored dead load \a m_forcing ITSELF — not a copy — so every downstream
     * operation is bit-identical to the pre-callback code. With a callback it returns a
     * reference to the freshly evaluated \a m_forcingEval, which is overwritten by the
     * next call; hold the reference only for as long as no other evaluation intervenes.
     *
     * @param[in]  U     State (displacements)
     * @param[in]  L     State (load factor)
     */
    const gsVector<T> & computeForcing(const gsVector<T> & U, const T & L)
    {
      if (!m_forcingFun) return m_forcing;
      if (!m_forcingFun(U, L, m_forcingEval))
        throw 2;
      GISMO_ENSURE(m_forcingEval.size()==m_forcing.size(),
                   "The forcing function returned "<<m_forcingEval.size()<<" entries, expected "<<m_forcing.size()<<".");
      return m_forcingEval;
    }

    /// Returns the forcing at the current iterate \f$(U+\Delta U,\ \Lambda+\Delta\Lambda)\f$,
    /// i.e. at the state where computeJacobian() evaluates the tangent. Mirrors the
    /// argument-less computeResidual()/computeJacobian() overloads. Note that the state
    /// vector is only formed when a callback is set, so the dead-load path is free.
    const gsVector<T> & computeForcing()
    {
      if (!m_forcingFun) return m_forcing;
      return this->computeForcing(m_U + m_DeltaU, m_L + m_DeltaL);
    }

    /**
     * @brief      Returns the forcing used by the constraint metric of the current step.
     *
     * The \f$\|f\|^2\f$-type scalings entering Crisfield's quadratic constraint (and the
     * scaling parameter \f$\varphi\f$) are FROZEN per step, at the predictor state, by
     * \a _freezeStepForcing. This is a definitional choice: it keeps the constraint a
     * fixed quadratic in \f$(\Delta U,\Delta\Lambda)\f$ throughout the corrector, so the
     * \f$\delta\Lambda\f$ root solve stays consistent from iteration to iteration. It is
     * NOT the dead-load approximation, which concerns \f$\delta u_t\f$ (see
     * \a setForcingFunction).
     *
     * Without a callback this returns \a m_forcing itself, hence bit-identical results.
     */
    const gsVector<T> & stepForcing() const
    {
      if (!m_forcingFun) return m_forcing;
      GISMO_ASSERT(m_forcingStep.size()==m_forcing.size(),"The step forcing is not initialized; it is frozen in initiateStep().");
      return m_forcingStep;
    }

    /// Freezes the constraint-metric forcing for the current step at the predictor
    /// state \f$(U,\Lambda)\f$; a no-op (and free) when no callback is set.
    void _freezeStepForcing()
    {
      if (!m_forcingFun) return;
      // m_U is sized from construction on: every subclass constructor calls initMethods(),
      // and the only caller of initiateStep() is _step(), which asserts m_initialized.
      if (!m_forcingFun(m_U, m_L, m_forcingStep))
        throw 2;
      // Checked here, at the cause, rather than later at an arbitrary stepForcing() read.
      GISMO_ENSURE(m_forcingStep.size()==m_forcing.size(),
                   "The forcing function returned "<<m_forcingStep.size()<<" entries, expected "<<m_forcing.size()<<".");
    }

    /// Compute the residual error norms.
    ///
    /// \param extendedSolve Selects the denominator floor applied to \a m_basisResidualU
    ///        at \a m_numIterations==0 (see the .hpp for the derivation). \c true only at
    ///        the extended-system solve's own call site (\a _extendedSystemSolve): that is
    ///        the one caller whose iteration-0 increment can be seeded already at the
    ///        answer (an extended solve started from an accurate \c localizeCrossing
    ///        output), which otherwise freezes the denominator at a round-off value and
    ///        makes \a m_residueU unreachable regardless of how converged the iterate is.
    ///        The ordinary corrector (\a _step, both of its call sites) always takes the
    ///        default \c false: its first increment is a genuine predictor step and must
    ///        not be floored.
    /// \note This is a \c virtual with a defaulted argument: the default is resolved by the
    ///       STATIC type at the call site, not by the dynamic type. That is harmless here
    ///       because no class in this hierarchy overrides \a computeResidualNorms — every
    ///       derived solver only re-exports it (\c using \c Base::computeResidualNorms;),
    ///       so there is exactly one definition and the two call sites in \a gsALMBase.hpp
    ///       are the only place the default is ever resolved.
    virtual void computeResidualNorms(bool extendedSolve = false);

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
    /// Dead load: the constant \f$-\partial R/\partial\Lambda\f$. Also the fallback (and
    /// the reference size) whenever no forcing callback is set.
    const gsVector<T>     m_forcing;

    /// OPT-IN state-dependent load derivative, empty by default; see \a setForcingFunction
    ALForce_t             m_forcingFun;
    /// Scratch for the per-iteration evaluation of \a m_forcingFun (unused when empty)
    gsVector<T>           m_forcingEval;
    /// Constraint-metric forcing, frozen once per step by \a _freezeStepForcing
    gsVector<T>           m_forcingStep;

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

    /// Outcome of the most recent limit-vs-branch classification (\a _testSingularPoint).
    struct SPverdict
    {
        enum type
        {
            NotTested  = -2,  ///< no classification has been performed since the last reset
            Unresolved = -1,  ///< the critical mode could not be resolved: NO verdict was formed
            Limit      =  0,  ///< |psi.f|/|f| >= SingularPointTestTol: limit point
            Branch     =  1,  ///< |psi.f|/|f| <  SingularPointTestTol: branch point
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
    /// The arc length that produced the CURRENT secant (m_U-m_Uprev, m_L-m_Lprev); the
    /// divisor that normalises the secant predictors. Valid once \a m_stepTaken is true;
    /// before that it merely mirrors \a m_arcLength. See setLength().
    T m_arcLength_prev;
    T m_arcLength_ori;
    bool m_adaptiveLength;
    /// True once at least one step has been ACCEPTED, i.e. once m_arcLength_prev carries
    /// the secant-producing arc length rather than just a configured value.
    bool m_stepTaken;

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

    /// AUTO's fold test function at the last jacobian=true stability evaluation; NaN when
    /// unavailable. See foldTestFunction().
    T m_foldTF;

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
    /// LEFT critical mode (right null vector of K_T^T); see \a _computeCriticalModeLeft and
    /// \a solutionVLeft(). Equal to \a m_V whenever the configured solver is self-adjoint.
    gsVector<T> m_Vleft;
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
    /// Convergence tolerance for the critical-mode inverse power iteration itself
    /// (option \c SingularPointModeTol); decoupled from \a m_SPTestTol, the verdict
    /// threshold it feeds. See \a defaultOptions / \a getOptions.
    T m_SPModeTol;

    /// Finite-difference perturbation \c del of \a computeBranchTangent (option
    /// \c BranchTangentPerturbation). See \a defaultOptions.
    T m_branchTangentDel;
    /// Relative degeneracy floor of \a computeBranchTangent (option
    /// \c BranchTangentTol). See \a defaultOptions.
    T m_branchTangentTol;

    /// Outcome of the most recent classification; see \a SPverdict. Seeded to
    /// \c SPverdict::NotTested in both constructors and reset at every entry to
    /// \a _testSingularPoint (same staleness precedent as \a m_converged).
    index_t m_SPverdict;
    /// Direction error |dV| achieved by the critical-mode iteration of the most recent
    /// classification (max over the right and left modes); see \a criticalModeError().
    T m_SPModeError;

    // Singular point computation tolerances
    T m_SPCompTolE; // extended iterations
    T m_SPCompTolB; // bisection method

    /// Probe budget of \a _bisectionSolve (option \c SingularPointBisIt), independent of
    /// \a m_maxIterations. See \a defaultOptions.
    index_t m_SPBisIt;
    /// Number of \a step() calls spent by the most recent \a _bisectionSolve. \c 0 is
    /// ambiguous: it means no bracket was found, or the routine has not run, or
    /// \c SingularPointBisIt \c <= 0 so the probe loop never iterated. See
    /// \a bisectionProbes(), which documents how to tell these apart.
    index_t m_SPBisProbes;

    // Opt-in (option SingularPointComposite): also require the equilibrium residuals
    // (TolF/TolU) at extended-solve termination. Default false is the ||K.V||-only test.
    bool m_SPComposite;

    // Branch switch parameter
    T m_tau;
};


} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsALMBase.hpp)
#endif
