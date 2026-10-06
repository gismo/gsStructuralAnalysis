 /** @file gsContinuationDriver.h

    @brief Predictor-corrector driver of the solver driver: one continuation
    step, a stateful step with retry policy, a load-interval solve and a pure
    step, over a gsStateSpace and a gsContinuationProblem.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#pragma once
#include <functional>
#include <string>
#include <vector>

#include <gsCore/gsLinearAlgebra.h>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsSolverDriverTypes.h>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsStateSpace.h>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsContinuationProblem.h>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsConvergenceTest.h>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsTangentSolver.h>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsPathConstraint.h>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsStepController.h>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsNewtonCorrector.h>

namespace gismo
{

/**
   @brief Continuation driver: advances the state of
   \f$R(U,\lambda) = \lambda F - F_{int}(U) = 0\f$ along a solution path, one
   constrained predictor-corrector step at a time.

   <b>Mathematics.</b> A step from the committed point \f$(U,\lambda)\f$
   seeks increments \f$(\Delta U,\Delta\lambda)\f$ with
   \f$R(U+\Delta U,\lambda+\Delta\lambda) = 0\f$ and one scalar constraint
   \f$g(\Delta U,\Delta\lambda;s) = 0\f$ (a gsPathConstraint). With the
   tangent \f$K = \partial F_{int}/\partial U = -\partial R/\partial U\f$ the
   Newton update solves the bordered system
   \f[
     \begin{pmatrix} K & -\partial R/\partial\lambda \\
                     \nabla_U g^T & \partial_\lambda g \end{pmatrix}
     \begin{pmatrix}\delta U\\ \delta\lambda\end{pmatrix}
     = \begin{pmatrix} R \\ -g \end{pmatrix},
   \f]
   because \f$R(U+\delta U) \approx R - K\,\delta U\f$: the correction solves
   \f$K\,\delta U = R\f$, never \f$K\,\delta U = -R\f$. Eliminating
   \f$\delta U = \delta\bar U + \delta\lambda\,u_t\f$ with
   \f$\delta\bar U = K^{-1}R\f$ and \f$u_t = K^{-1}\partial R/\partial\lambda\f$
   costs two back-substitutions with one factorization; the constraint owns
   the remaining scalar equation for \f$\delta\lambda\f$ (closed form or root
   choice).

   <b>Protocol of one attempt</b> (the constraint, test, controller and
   corrector are called exactly in this order):
   -# rollback() the state to the committed point, reset the step context;
   -# optional pre-step hook (may resize the step of this attempt);
   -# predictor: when the constraint needs the load tangent, assemble
      \f$K\f$ at the committed point, factorize, solve
      \f$K u_t = \partial R/\partial\lambda\f$; the constraint writes
      \f$(\Delta U,\Delta\lambda)\f$ and the driver applies \f$\Delta U\f$
      with scale 1;
   -# residual and constraint residual at the predicted point, then
      gsConvergenceTest::initStep() once and test() at iteration 0;
   -# while the test is Unconverged: assemble and factorize \f$K\f$ when the
      gsNewtonCorrector schedule asks, solve \f$K\delta\bar U = R\f$ and
      \f$K u_t = \partial R/\partial\lambda\f$, let the constraint choose
      \f$\delta\lambda\f$, apply the (relaxed) update, refresh the residuals,
      test, then call the optional post-iteration hook;
   -# Converged: commit() the state, add \f$\Delta\lambda\f$ to \f$\lambda\f$
      and write the step history; any failure: rollback().

   The committed state and \f$\lambda\f$ never change on failure. The state is
   reached only through gsStateSpace (accumulation is in place, no
   \f$U+\Delta U\f$ is ever formed); the problem only through
   gsContinuationProblem. Each Jacobian assembly gets a fresh version stamp
   from one driver-wide counter and is followed by exactly one
   gsTangentSolver::factorize(), so there is one factorization per Jacobian.

   <b>Failure policy.</b> No member throws. Every operation returns a status
   and fills a gsSolverReport; GISMO_ENSURE is used only in the setters and
   in setup(). The driver itself never inspects norms for
   non-finite values: a non-finite iterate runs on until a test ends it.
   Termination is guaranteed instead by setup() requiring a MaxIter guard
   (gsConvergenceTest::maxIterations() >= 0) in the test tree; non-finite
   iterates fail earlier only through an opt-in gsDivergedTest (status
   NotConverged, failedOp Test) or the tangent solver's CheckFinite (status
   SolverError, failedOp Solve).

   Report mapping (report.step is the number of accepted steps before the
   failing one; report.iteration the corrector updates applied at the failing
   operation, 0 for the predictor):

   | event | status | failedOp | retried by step() |
   |---|---|---|---|
   | residual false, or residual size differs from n | AssemblyError | Residual | yes |
   | jacobian false | AssemblyError | Jacobian | yes |
   | loadDirection false | AssemblyError | LoadDirection | yes |
   | factorize false | SolverError | Factorize | yes |
   | solve false (including CheckFinite) | SolverError | Solve | yes |
   | test tree Exhausted (iteration limit) | NotConverged | None, note "iteration limit" | yes |
   | test tree Failed (gsDivergedTest) | NotConverged | Test, note = test name | yes |
   | pre-step / post-iteration hook Reject | NotConverged | Hook | yes |
   | pre-step / post-iteration hook Abort | OtherError | Hook | no |
   | predict/correct false, or predicted DeltaU of size other than n | OtherError | Constraint | no |
   | apply false (predictor or corrector) | OtherError | Apply | no |
   | commit false | OtherError | Commit | no |
   | rollback false at the start of an attempt | OtherError | Rollback | no |
   | importCommitted false or space not copyable (pure step) | OtherError | Import | no |
   | exportWorking false (pure step, after commit) | OtherError | Export | no |
   | load stall or length stall (step()) | NotConverged | Controller | no |
   | setup() not called | OtherError | None | no |

   <b>Stateful step.</b> step() runs attempts until one is accepted or the
   policy gives up: a retryable failure is retried after the controller has
   halved the length (gsStepController::rejected()), and, when the Load secant
   predictor is on and was used, first once at the same length with the
   tangent predictor. With the default controller (no bisection) step() is
   exactly one attempt.

   <b>Load secant predictor.</b> With setSecantPredictor(true) a step whose
   constraint has stepIsLoadIncrement() extrapolates the previous increment,
   \f$\Delta U = s\,\Delta U_{prev}/s_{prev}\f$, \f$\Delta\lambda = s\f$
   (first-order in the load), and skips the Jacobian and factorization at the
   committed point. It is valid only where the step parameter is the load
   increment, so setup() requires stepIsLoadIncrement(). Arc-length
   constraints keep their own secant in gsStepContext::prev, which the driver
   never invalidates on a retry.

   <b>Pure step.</b> The overload step(from, hist, ds, to, histOut) imports
   @a from into the state space, runs one attempt without retries and exports
   the result. "Pure" means that the driver's own members (lambda, history,
   controller, secant flag, report) are unchanged; the committed state of the
   space becomes @a to on success and is @a from after a failed attempt (see
   the member's doc for the import/export failure branches), constraints may keep
   internal state, and the tests keep per-attempt bases. Callers mixing pure
   and stateful steps on one space resynchronise with setLambda() and
   setHistory().

   <b>Complexity.</b> Per corrector update: one residual assembly, at most one
   Jacobian assembly and one factorization (per the refactorization
   schedule), two triangular solves, O(n) vector work and one
   gsStateSpace::apply(). The predictor costs one Jacobian assembly,
   factorization and solve when the constraint needs the tangent, nothing in
   secant mode.

   <b>Thread safety.</b> The work vectors and the tangent solver are mutable
   members: a driver is used by one thread at a time, also through the pure
   step. The driver does not own the problem or the state space; both must
   outlive it. It is not copyable.

   @tparam T coefficient type

   \ingroup gsStructuralAnalysis
*/
template<class T>
class gsContinuationDriver
{
public:
    /// Observer of step events; cannot fail a step
    typedef std::function<void(const gsStepEvent<T> &)>              Observer_t;
    /// Pre-step hook: receives the committed lambda and the step of this
    /// attempt, which it may resize (for this attempt only)
    typedef std::function<gsHookDecision(T lambda, T & stepSize)>    PreStep_t;
    /// Post-iteration hook: called after every corrector update and its test
    typedef std::function<gsHookDecision(const gsStepContext<T> &)>  PostIteration_t;

    /// Non-owning: @a P and @a S must outlive the driver. Non-const references,
    /// so temporaries cannot bind. O(1).
    gsContinuationDriver(gsContinuationProblem<T> & P, gsStateSpace<T> & S);

    /// Sets the path constraint; invalidates setup(). Requires a non-null pointer.
    void setConstraint(const typename gsPathConstraint<T>::Ptr & c);
    /// The path constraint. O(1).
    const typename gsPathConstraint<T>::Ptr & constraint() const { return m_constraint; }

    /// Sets the convergence-test tree; invalidates setup(). Requires a non-null pointer.
    void setTest(const typename gsConvergenceTest<T>::Ptr & t);
    /// The convergence-test tree. O(1).
    const typename gsConvergenceTest<T>::Ptr & test() const { return m_test; }

    /// The step controller (length, adaptivity, retry budget). O(1).
    gsStepController<T> & controller() { return m_controller; }
    /// The step controller. O(1).
    const gsStepController<T> & controller() const { return m_controller; }

    /// The Newton corrector policy (refactorization schedule, relaxation). O(1).
    gsNewtonCorrector<T> & corrector() { return m_corrector; }
    /// The Newton corrector policy. O(1).
    const gsNewtonCorrector<T> & corrector() const { return m_corrector; }

    /// The tangent solver (solver choice, equilibration, factorization counter). O(1).
    gsTangentSolver<T> & linearSolver() { return m_lin; }
    /// The tangent solver. O(1).
    const gsTangentSolver<T> & linearSolver() const { return m_lin; }

    /// Switches the Load secant predictor; invalidates setup(). Default off. O(1).
    void setSecantPredictor(bool on);
    /// Whether the Load secant predictor is on. O(1).
    bool secantPredictor() const { return m_secant; }

    /// Appends an observer; observers are called in insertion order. O(1).
    void addObserver(const Observer_t & f);
    /// Removes all observers. O(number of observers).
    void clearObservers() { m_observers.clear(); }
    /// Sets the pre-step hook; an empty function means none. O(1).
    void setPreStep(const PreStep_t & f) { m_preStep = f; }
    /// Sets the post-iteration hook; an empty function means none. O(1).
    void setPostIteration(const PostIteration_t & f) { m_postIteration = f; }

    /**
       @brief Checks the configuration and sizes the step context.

       Requires (GISMO_ENSURE): the problem's residual and jacobian are set; a
       constraint and a test are set; <tt>P.F.size() == S.tangentSize()</tt>;
       the test tree contains a MaxIter guard (this is what guarantees
       termination of a step); the Load secant predictor is only on for a
       constraint with stepIsLoadIncrement(). Then copies P.F into the context.
       lambda, the step history and the step counter are not reset. Call
       again after changing the problem's F.

       Cost: O(n).
    */
    bool setup();

    /// One continuation step with the retry policy, see the class doc. Returns
    /// the status, also available in report(). Cost: see the class doc, per attempt.
    gsStatus step();

    /**
       @brief Steps until lambda reaches @a lambdaEnd (load-increment
       constraints only).

       The last step is clamped to land on @a lambdaEnd (the increment is
       @a lambdaEnd - lambda, exact whenever that subtraction is exact, see
       gsStepController::targetLoad()). Returns Success without a step when
       lambda is already at @a lambdaEnd, in the direction of the step length.
       A non-finite lambda, @a lambdaEnd or step length returns NotConverged
       with failedOp Controller. Stops at the first failing step.

       Cost: the sum of the steps taken, each as for step(); O(1) when no step
       is needed.
    */
    gsStatus solve(T lambdaEnd);

    /**
       @brief Pure step: one attempt from @a from with step size @a ds.

       Needs a copyable state space. The committed state of the space is set
       to @a from, the attempt is run once (no retry, no controller call), and
       on success @a to, @a histOut hold the new point and history. On a failed
       attempt @a to and @a histOut are untouched and the space's committed and
       working state equal @a from (the attempt rolled back). If exportWorking()
       fails after a converged attempt, the space has already committed the new
       point, @a histOut and @a to.lambda are untouched and @a to.U is
       unspecified (failedOp Export). If setup() was not called, the space is
       not copyable or importCommitted() fails, no attempt is run, @a to and
       @a histOut are untouched and the space is as the failed import (or the
       caller) left it. The report's step is 0. See the class doc for what "pure" covers.

       @param[in] from start point (exported coordinates, lambda)
       @param[in] hist step history of @a from
       @param[in] ds step size of this step
       @param[out] to new point
       @param[out] histOut new step history
       @param[out] report optional report

       Cost: one attempt (see the class doc) plus O(n) for the import and export.
    */
    gsStatus step(const gsPathPoint<T> & from, const gsStepHistory<T> & hist, T ds,
                  gsPathPoint<T> & to, gsStepHistory<T> & histOut,
                  gsSolverReport * report = NULL) const;

    /// Committed load factor. O(1).
    T lambda() const { return m_ctx.lambda; }
    /// Sets the committed load factor (the state is not touched). O(1).
    void setLambda(T lambda) { m_ctx.lambda = lambda; }

    /// History of the last accepted step. O(1).
    const gsStepHistory<T> & history() const { return m_ctx.prev; }
    /// Sets the step history; the Load secant predictor is armed iff it is usable. O(n).
    void setHistory(const gsStepHistory<T> & h);
    /// Clears the step history and disarms the Load secant predictor. O(1).
    void resetHistory();

    /// Corrector updates of the last attempt (N after an iteration-limit failure). O(1).
    index_t numIterations() const { return m_ctx.iteration; }
    /// The step context of the stateful path. O(1).
    const gsStepContext<T> & context() const { return m_ctx; }
    /// Report of the last stateful call. O(1).
    const gsSolverReport & report() const { return m_ctx.report; }

private:
    gsContinuationDriver(const gsContinuationDriver &);
    gsContinuationDriver & operator=(const gsContinuationDriver &);

    gsStatus _step(bool clampToEnd, T lambdaEnd);
    gsStatus _attempt(gsStepContext<T> & ctx, T s, bool useSecant) const;
    gsStatus _fail(gsStepContext<T> & ctx, gsStatus st, gsDriverOp op, const std::string & note) const;
    void _emit(typename gsStepEvent<T>::Kind kind, const gsStepContext<T> & ctx) const;
    bool _evaluate(gsStepContext<T> & ctx, index_t n) const;
    gsDriverOp _tangent(gsStepContext<T> & ctx, T lambdaK) const;
    bool _secantUsable(const gsStepHistory<T> & h) const;
    static gsStatus _statusOf(gsDriverOp op);
    static bool _notFinite(T x) { return !((x - x) == T(0)); }

    gsContinuationProblem<T> * m_P;
    gsStateSpace<T> * m_S;
    typename gsPathConstraint<T>::Ptr m_constraint;
    typename gsConvergenceTest<T>::Ptr m_test;
    gsStepController<T> m_controller;
    gsNewtonCorrector<T> m_corrector;
    mutable gsTangentSolver<T> m_lin;
    /// Driver-wide Jacobian version counter, shared by the stateful context and
    /// every context of a pure step: a context counting from zero would hand the
    /// tangent solver a version it already holds and reuse a stale factorization.
    mutable index_t m_Kstamp;
    mutable gsVector<T> m_dUt, m_dUbar;
    gsVector<T> m_none; ///< empty vector passed as u_t when the predictor needs none
    bool m_secant, m_secantOK, m_setup;
    std::vector<Observer_t> m_observers;
    PreStep_t m_preStep;
    PostIteration_t m_postIteration;
    gsStepContext<T> m_ctx;
};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsContinuationDriver.hpp)
#endif
