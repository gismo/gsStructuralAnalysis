 /** @file gsSolverDriverTypes.h

    @brief Enumerations and plain structs shared by the parts of the solver
    driver (step context, history, report, observer payload).

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#pragma once
#include <string>

#include <gsCore/gsLinearAlgebra.h>
#include <gsStructuralAnalysis/src/gsStructuralAnalysisTools/gsStructuralAnalysisTypes.h>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsStateSpace.h>

namespace gismo
{

/// Mode of a residual evaluation
enum class gsEvalMode
{
    Accept, ///< Normal, write-enabled evaluation; may directly precede gsStateSpace::commit()
    Probe   ///< Write-disabled trial evaluation (e.g. for a line search); never directly precedes commit()
};

/// The operation that failed, recorded in gsSolverReport::failedOp
enum class gsDriverOp
{
    None,          ///< No failure
    Residual,      ///< The problem's residual op returned false
    Jacobian,      ///< The problem's jacobian op returned false
    LoadDirection, ///< The problem's loadDirection op returned false
    Mass,          ///< The mass op returned false (time driver)
    Damping,       ///< The damping op returned false (time driver)
    Apply,         ///< gsStateSpace::apply returned false
    Commit,        ///< gsStateSpace::commit returned false
    Rollback,      ///< gsStateSpace::rollback returned false
    Import,        ///< gsStateSpace::importCommitted returned false
    Export,        ///< gsStateSpace::exportWorking returned false
    Factorize,     ///< Factorization of the tangent failed
    Solve,         ///< A linear solve failed
    Constraint,    ///< The constraint's predict or correct returned false
    Test,          ///< A convergence test reported Failed (non-finite or diverged)
    Controller,    ///< The step controller refused a retry (depth, budget or stall)
    Hook           ///< A hook returned gsHookDecision::Abort
};

/// Decision returned by a pre-step or post-iteration hook
enum class gsHookDecision
{
    Continue, ///< Proceed
    Reject,   ///< Roll back and treat the step as rejected, so the controller may retry
    Abort     ///< Roll back and return with failedOp = gsDriverOp::Hook
};

/// Outcome of a driver call
struct gsSolverReport
{
    gsStatus    status   = gsStatus::NotStarted; ///< Overall status
    gsDriverOp  failedOp = gsDriverOp::None;     ///< Operation that failed, None on success
    index_t     step      = -1;                  ///< Step of the failure (counted as gsStepContext::step), -1 if none
    index_t     iteration = -1;                  ///< Corrector update of the failure (counted as gsStepContext::iteration), -1 if none
    std::string note;                            ///< Free text, e.g. the reason given by a hook

    /// Restore the defaults
    void reset()
    {
        status    = gsStatus::NotStarted;
        failedOp  = gsDriverOp::None;
        step      = -1;
        iteration = -1;
        note.clear();
    }
};

/// The last accepted step, written by the driver at every commit
template<class T>
struct gsStepHistory
{
    gsVector<T> DeltaU;        ///< \f$\Delta U\f$ of the accepted step, tangent size
    T           DeltaL = 0;    ///< \f$\Delta\lambda\f$ of the accepted step
    T           stepSize = 0;  ///< Step size \f$s\f$ of the accepted step
    index_t     iterations = 0;///< Corrector updates the accepted step needed
    bool        valid = false; ///< False until the first step is accepted

    /// Restore the defaults (DeltaU empty)
    void reset()
    {
        DeltaU.resize(0);
        DeltaL     = 0;
        stepSize   = 0;
        iterations = 0;
        valid      = false;
    }
};

/// A point on the solution path in exported coordinates. Only meaningful for
/// spaces with gsStateSpace::copyable() true.
template<class T>
struct gsPathPoint
{
    gsVector<T> U;      ///< Exported state, of the <em>export</em> size of gsStateSpace::exportWorking (not necessarily the tangent size)
    T           lambda = 0; ///< Load factor
};

/**
   @brief Tangent-space scratch of one continuation step.

   Owned by the driver; the parts (constraint, corrector, controller) receive
   it by reference, tests and observers by const reference. It never holds
   the state, which lives in the gsStateSpace. All vectors have the tangent
   size \f$n\f$ unless stated otherwise.

   Semantics the parts rely on:
   - @b iteration counts corrector updates: 0 right after the predictor, k
     after the k-th application of an update. <tt>MaxIter = N</tt> means N
     corrector updates; the step fails with NotConverged only if the
     convergence test evaluated after the N-th update still reports
     Unconverged.
   - @b Iteration @b 0: after the predictor, DeltaU and DeltaL hold the
     predictor increment, deltaU = DeltaU, deltaL = DeltaL, and R and
     constraintResidual are evaluated at the predicted point.
   - @b constraintResidual is the value of the constraint's residual(ctx), so
     a test can form \f$|g|/s^2\f$ from constraintResidual and stepSize alone.
   - @b relaxation: the applied update is relaxation*(deltaU, deltaL). The
     corrector sets relaxation from its option before calling the
     constraint's correct(), and correct() may lower it (e.g. on an
     oscillating \f$\delta\lambda\f$) or reset it to 1. After the corrector
     applies the update, deltaU and deltaL hold the <em>applied (relaxed)</em>
     increment, which is what a relative-update test reads.

   @tparam T coefficient type

   \ingroup gsStructuralAnalysis
*/
template<class T>
struct gsStepContext
{
    T lambda = 0;              ///< \f$\lambda\f$ at the last commit
    T DeltaL = 0;              ///< \f$\Delta\lambda\f$ accumulated in this step (predictor and corrections)
    T deltaL = 0;              ///< \f$\delta\lambda\f$ of the latest update
    T stepSize = 0;            ///< Step size \f$s\f$: the \f$\Delta\lambda\f$ target (load/displacement control) or the arc length; the constraint's stepIsLoadIncrement() says which
    gsVector<T> DeltaU;        ///< \f$\Delta U\f$ accumulated in this step, tangent size
    gsVector<T> deltaU;        ///< \f$\delta U\f$ of the latest update, tangent size
    gsVector<T> R;             ///< Residual \f$\lambda_w F - F_{int}\f$ at the working state, \f$\lambda_w\f$ = lambda + DeltaL, tangent size
    gsVector<T> F;             ///< Reference external force (copy of the problem's F), tangent size
    gsVector<T> dRdl;          ///< \f$\partial R/\partial\lambda\f$ at the working state; equals F when the problem has no loadDirection
    gsSparseMatrix<T> K;       ///< Latest tangent stiffness \f$\partial F_{int}/\partial U\f$, n x n
    index_t Kversion = 0;      ///< Incremented each time K is reassembled; the tangent solver factorizes once per version
    index_t step = 0;          ///< Index of the current step (number of accepted steps before it)
    index_t iteration = 0;     ///< Corrector updates applied in this step: 0 right after the predictor, k after the k-th update
    T constraintResidual = 0;  ///< Constraint value \f$g(\Delta U,\Delta\lambda)\f$ at the working point, refreshed after the predictor and after every correction; 0 for fixed lambda / load control
    T relaxation = 1;          ///< Factor applied to (deltaU, deltaL) before the update is applied
    gsStepHistory<T> prev;     ///< The last accepted step
    gsSolverReport report;     ///< Status of the current or last step

    /// Prepare for a new step: zero DeltaL, deltaL, iteration and
    /// constraintResidual, set relaxation to 1 and zero DeltaU and deltaU
    /// (keeping their size). lambda, prev, F, K and Kversion are untouched.
    void resetStep()
    {
        DeltaL = 0;
        deltaL = 0;
        iteration = 0;
        constraintResidual = 0;
        relaxation = 1;
        DeltaU.setZero();
        deltaU.setZero();
    }
};

/**
   @brief Read-only payload handed to step observers.

   Observers return void and cannot fail a step. They see the state only
   through the const surface of gsStateSpace (exportWorking, magnitude, sizes).
   The event holds references and is only valid during the observer call.

   \ingroup gsStructuralAnalysis
*/
template<class T>
struct gsStepEvent
{
    /// What happened
    enum Kind
    {
        StepStarted,   ///< A step began (after the initial rollback)
        IterationDone, ///< A corrector update was applied and tested
        StepAccepted,  ///< The step converged and was committed
        StepRejected   ///< The step was rolled back
    };

    /// Constructor
    gsStepEvent(Kind k, const gsStepContext<T>& c, const gsStateSpace<T>& s)
    : kind(k), ctx(c), space(s)
    {}

    Kind kind;                     ///< Event kind
    const gsStepContext<T>& ctx;   ///< Step scratch
    const gsStateSpace<T>& space;  ///< State, const access only
};

} // namespace gismo
