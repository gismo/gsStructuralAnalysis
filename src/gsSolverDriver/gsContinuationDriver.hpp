 /** @file gsContinuationDriver.hpp

    @brief Implementation of the continuation driver.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#pragma once
#include <cmath>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsContinuationDriver.h>

namespace gismo
{

template<class T>
gsContinuationDriver<T>::gsContinuationDriver(gsContinuationProblem<T> & P, gsStateSpace<T> & S)
: m_P(&P), m_S(&S), m_Kstamp(0), m_committedK(0), m_secant(false), m_secantOK(false), m_setup(false)
{ }

template<class T>
void gsContinuationDriver<T>::setConstraint(const typename gsPathConstraint<T>::Ptr & c)
{
    GISMO_ENSURE(c, "gsContinuationDriver::setConstraint: null constraint");
    m_constraint = c;
    m_setup = false;
}

template<class T>
void gsContinuationDriver<T>::setTest(const typename gsConvergenceTest<T>::Ptr & t)
{
    GISMO_ENSURE(t, "gsContinuationDriver::setTest: null convergence test");
    m_test = t;
    m_setup = false;
}

template<class T>
void gsContinuationDriver<T>::setSecantPredictor(bool on)
{
    m_secant = on;
    m_setup = false;
}

template<class T>
void gsContinuationDriver<T>::addObserver(const Observer_t & f)
{
    GISMO_ENSURE(bool(f), "gsContinuationDriver::addObserver: empty observer");
    m_observers.push_back(f);
}

template<class T>
bool gsContinuationDriver<T>::setup()
{
    GISMO_ENSURE(bool(m_P->residual), "gsContinuationDriver::setup: the problem has no residual");
    GISMO_ENSURE(bool(m_P->jacobian), "gsContinuationDriver::setup: the problem has no jacobian");
    GISMO_ENSURE(m_constraint, "gsContinuationDriver::setup: no path constraint set");
    GISMO_ENSURE(m_test, "gsContinuationDriver::setup: no convergence test set");
    const index_t n = m_S->tangentSize();
    GISMO_ENSURE(m_P->F.size() == n, "gsContinuationDriver::setup: F has " << m_P->F.size()
                 << " entries, the state space has tangent size " << n);
    GISMO_ENSURE(m_test->maxIterations() >= 0,
                 "gsContinuationDriver::setup: the convergence test has no MaxIter guard");
    GISMO_ENSURE(!m_secant || m_constraint->stepIsLoadIncrement(),
                 "gsContinuationDriver::setup: the secant predictor needs a load-increment constraint");

    m_ctx.F = m_P->F;
    m_ctx.DeltaU = gsVector<T>::Zero(n);
    m_ctx.deltaU = gsVector<T>::Zero(n);
    m_ctx.dRdl = m_ctx.F;
    m_committedK = 0;
    m_setup = true;
    return true;
}

template<class T>
bool gsContinuationDriver<T>::assembleCommittedTangent()
{
    m_committedK = 0;
    if (!m_setup || !m_S->rollback())
        return false;

    // The predictor of the next step sees deltaU = 0 (resetStep()); the
    // operators may read it, so it is zeroed for the evaluation and restored.
    m_deltaUsave = m_ctx.deltaU;
    m_ctx.deltaU.setZero(m_S->tangentSize());
    const gsDriverOp op = _tangent(m_ctx, m_ctx.lambda);
    m_ctx.deltaU = m_deltaUsave;
    if (op != gsDriverOp::None)
        return false;
    m_committedK = m_ctx.Kversion;
    return true;
}

template<class T>
void gsContinuationDriver<T>::setHistory(const gsStepHistory<T> & h)
{
    m_ctx.prev = h;
    m_secantOK = _secantUsable(h);
}

template<class T>
void gsContinuationDriver<T>::resetHistory()
{
    m_ctx.prev.reset();
    m_secantOK = false;
    m_committedK = 0;
}

template<class T>
bool gsContinuationDriver<T>::_secantUsable(const gsStepHistory<T> & h) const
{
    return h.valid && h.DeltaU.size() == m_S->tangentSize()
           && h.stepSize != T(0) && std::isfinite(h.stepSize);
}

template<class T>
gsStatus gsContinuationDriver<T>::_statusOf(gsDriverOp op)
{
    switch (op)
    {
    case gsDriverOp::Factorize:
    case gsDriverOp::Solve:
        return gsStatus::SolverError;
    default:
        return gsStatus::AssemblyError;
    }
}

template<class T>
void gsContinuationDriver<T>::_emit(typename gsStepEvent<T>::Kind kind, const gsStepContext<T> & ctx) const
{
    for (size_t i = 0; i < m_observers.size(); ++i)
    {
        const gsStepEvent<T> ev(kind, ctx, *m_S);
        m_observers[i](ev);
    }
}

// The first failure of an attempt is the one reported: every failure returns
// at once, after the rollback.
template<class T>
gsStatus gsContinuationDriver<T>::_fail(gsStepContext<T> & ctx, gsStatus st, gsDriverOp op,
                                        const std::string & note) const
{
    ctx.report.status    = st;
    ctx.report.failedOp  = op;
    ctx.report.step      = ctx.step;
    ctx.report.iteration = ctx.iteration;
    ctx.report.note      = note;
    if (!m_S->rollback())
        ctx.report.note += note.empty() ? "rollback failed" : "; rollback failed";
    _emit(gsStepEvent<T>::StepRejected, ctx);
    return st;
}

// Residual at the working state with lambda_w = lambda + DeltaL, then the
// constraint residual g(DeltaU, DeltaL; s).
template<class T>
bool gsContinuationDriver<T>::_evaluate(gsStepContext<T> & ctx, index_t n) const
{
    if (!m_P->residual(ctx.lambda + ctx.DeltaL, ctx.R, gsEvalMode::Accept))
        return false;
    if (ctx.R.size() != n)
        return false;
    ctx.constraintResidual = m_constraint->residual(ctx);
    return true;
}

// Assembles K (and dR/dlambda when the problem has a load direction) at the
// working state and factorizes it: one factorize() per assembly, with a
// version drawn from the driver-wide counter. Returns the failing operation,
// gsDriverOp::None on success.
template<class T>
gsDriverOp gsContinuationDriver<T>::_tangent(gsStepContext<T> & ctx, T lambdaK) const
{
    if (!m_P->jacobian(lambdaK, ctx.K))
        return gsDriverOp::Jacobian;
    ctx.Kversion = ++m_Kstamp;
    if (!m_lin.factorize(ctx.K, ctx.Kversion))
        return gsDriverOp::Factorize;
    if (m_P->loadDirection && !m_P->loadDirection(lambdaK, ctx.dRdl))
        return gsDriverOp::LoadDirection;
    return gsDriverOp::None;
}

template<class T>
gsStatus gsContinuationDriver<T>::_attempt(gsStepContext<T> & ctx, T s, bool useSecant) const
{
    const index_t n = m_S->tangentSize();

    // Every attempt starts from the committed state.
    ctx.resetStep();
    ctx.stepSize = s;
    if (!m_S->rollback())
    {
        ctx.report.status    = gsStatus::OtherError;
        ctx.report.failedOp  = gsDriverOp::Rollback;
        ctx.report.step      = ctx.step;
        ctx.report.iteration = ctx.iteration;
        ctx.report.note      = "rollback failed";
        _emit(gsStepEvent<T>::StepRejected, ctx);
        return gsStatus::OtherError;
    }
    _emit(gsStepEvent<T>::StepStarted, ctx);

    if (m_preStep)
    {
        const gsHookDecision d = m_preStep(ctx.lambda, s);
        if (d == gsHookDecision::Reject)
            return _fail(ctx, gsStatus::NotConverged, gsDriverOp::Hook, "pre-step hook rejected");
        if (d == gsHookDecision::Abort)
            return _fail(ctx, gsStatus::OtherError, gsDriverOp::Hook, "pre-step hook aborted");
    }
    ctx.stepSize = s;
    if (!m_P->loadDirection)
        ctx.dRdl = ctx.F;

    // Predictor. The Load secant direction is the previous increment per unit
    // load; otherwise u_t = K^{-1} dR/dlambda at the committed state when the
    // constraint needs it.
    const gsVector<T> * dUt = &m_none;
    if (useSecant)
    {
        m_dUt = ctx.prev.DeltaU / ctx.prev.stepSize;
        dUt = &m_dUt;
    }
    else if (m_constraint->predictorNeedsTangent(ctx))
    {
        // The factorization left by assembleCommittedTangent() is that of the
        // matrix assembled here: same state, lambda and deltaU = 0.
        const bool reuse = m_committedK != 0 && ctx.Kversion == m_committedK
                           && m_lin.factorized() && m_lin.version() == m_committedK;
        if (!reuse)
        {
            const gsDriverOp op = _tangent(ctx, ctx.lambda);
            if (op != gsDriverOp::None)
                return _fail(ctx, _statusOf(op), op, "");
        }
        if (!m_lin.solve(ctx.dRdl, m_dUt))
            return _fail(ctx, gsStatus::SolverError, gsDriverOp::Solve, "");
        dUt = &m_dUt;
    }
    if (!m_constraint->predict(ctx, *dUt, *m_S))
        return _fail(ctx, gsStatus::OtherError, gsDriverOp::Constraint, "predict failed");
    if (ctx.DeltaU.size() != n)
        return _fail(ctx, gsStatus::OtherError, gsDriverOp::Constraint, "predicted increment has the wrong size");
    if (!m_S->apply(ctx.DeltaU, T(1)))
        return _fail(ctx, gsStatus::OtherError, gsDriverOp::Apply, "predictor apply failed");

    if (!_evaluate(ctx, n))
        return _fail(ctx, gsStatus::AssemblyError, gsDriverOp::Residual, "");

    m_test->initStep(ctx, *m_S);
    gsConvergenceStatus st = m_test->test(ctx, *m_S);

    // Corrector: K dUbar = R (R = lambda F - F_int, K = dF_int/dU).
    while (st == gsConvergenceStatus::Unconverged)
    {
        if (m_corrector.refactorize(ctx))
        {
            const gsDriverOp op = _tangent(ctx, ctx.lambda + ctx.DeltaL);
            if (op != gsDriverOp::None)
                return _fail(ctx, _statusOf(op), op, "");
        }
        if (!m_lin.solve(ctx.R, m_dUbar))
            return _fail(ctx, gsStatus::SolverError, gsDriverOp::Solve, "");
        if (!m_lin.solve(ctx.dRdl, m_dUt))
            return _fail(ctx, gsStatus::SolverError, gsDriverOp::Solve, "");

        m_corrector.prepare(ctx);
        if (!m_constraint->correct(ctx, m_dUbar, m_dUt))
            return _fail(ctx, gsStatus::OtherError, gsDriverOp::Constraint, "correct failed");
        if (!m_corrector.apply(ctx, *m_S))
            return _fail(ctx, gsStatus::OtherError, gsDriverOp::Apply, "corrector apply failed");

        if (!_evaluate(ctx, n))
            return _fail(ctx, gsStatus::AssemblyError, gsDriverOp::Residual, "");
        st = m_test->test(ctx, *m_S);
        _emit(gsStepEvent<T>::IterationDone, ctx);

        if (m_postIteration)
        {
            const gsHookDecision d = m_postIteration(ctx);
            if (d == gsHookDecision::Reject)
                return _fail(ctx, gsStatus::NotConverged, gsDriverOp::Hook, "post-iteration hook rejected");
            if (d == gsHookDecision::Abort)
                return _fail(ctx, gsStatus::OtherError, gsDriverOp::Hook, "post-iteration hook aborted");
        }
    }

    if (st == gsConvergenceStatus::Exhausted)
        return _fail(ctx, gsStatus::NotConverged, gsDriverOp::None, "iteration limit");
    if (st == gsConvergenceStatus::Failed)
        return _fail(ctx, gsStatus::NotConverged, gsDriverOp::Test, m_test->name());

    // Converged: the last residual was evaluated in Accept mode at this point.
    if (!m_S->commit())
        return _fail(ctx, gsStatus::OtherError, gsDriverOp::Commit, "commit failed");
    // The committed point moves on, whatever the solver still holds.
    m_committedK = 0;
    ctx.lambda += ctx.DeltaL;
    ctx.prev.DeltaU     = ctx.DeltaU;
    ctx.prev.DeltaL     = ctx.DeltaL;
    ctx.prev.stepSize   = ctx.stepSize;
    ctx.prev.iterations = ctx.iteration;
    ctx.prev.valid      = true;
    ctx.report.status    = gsStatus::Success;
    ctx.report.failedOp  = gsDriverOp::None;
    ctx.report.step      = ctx.step;
    ctx.report.iteration = ctx.iteration;
    ctx.report.note.clear();
    _emit(gsStepEvent<T>::StepAccepted, ctx);
    ++ctx.step;
    return gsStatus::Success;
}

template<class T>
gsStatus gsContinuationDriver<T>::step()
{
    return _step(false, T(0));
}

// One step with the retry policy. A secant attempt that fails retryably is
// repeated once at the same length with the tangent predictor; after that the
// controller decides on halving. ctx.prev is only written at commit.
template<class T>
gsStatus gsContinuationDriver<T>::_step(bool clampToEnd, T lambdaEnd)
{
    if (!m_setup)
    {
        m_ctx.report.reset();
        m_ctx.report.status = gsStatus::OtherError;
        m_ctx.report.note   = "setup() not called";
        return gsStatus::OtherError;
    }
    m_ctx.report.reset();
    const bool load = m_constraint->stepIsLoadIncrement();
    bool trySecant = m_secant && m_secantOK && _secantUsable(m_ctx.prev);

    for (;;)
    {
        T s = m_controller.length();
        T target;
        if (clampToEnd)
        {
            target = gsStepController<T>::targetLoad(m_ctx.lambda, s, lambdaEnd);
            s = target - m_ctx.lambda;
        }
        else
            target = m_ctx.lambda + s;

        if (load && gsStepController<T>::loadStalled(m_ctx.lambda, target))
        {
            m_ctx.report.status    = gsStatus::NotConverged;
            m_ctx.report.failedOp  = gsDriverOp::Controller;
            m_ctx.report.step      = m_ctx.step;
            m_ctx.report.iteration = -1;
            m_ctx.report.note      = "load step does not change lambda";
            m_secantOK = false;
            return gsStatus::NotConverged;
        }

        const gsStatus st = _attempt(m_ctx, s, trySecant);
        if (st == gsStatus::Success)
        {
            m_controller.accepted(m_ctx.prev.iterations);
            m_secantOK = true;
            return st;
        }
        if (st == gsStatus::OtherError)
        {
            m_secantOK = false;
            return st;
        }
        if (trySecant)
        {
            trySecant = false;
            continue;
        }
        if (!m_controller.rejected())
        {
            m_ctx.report.note += m_ctx.report.note.empty() ? "no retry (controller)" : "; no retry (controller)";
            m_secantOK = false;
            return st;
        }
        m_secantOK = false;
        if (!load && gsStepController<T>::lengthStalled(m_controller.length()))
        {
            m_ctx.report.status    = gsStatus::NotConverged;
            m_ctx.report.failedOp  = gsDriverOp::Controller;
            m_ctx.report.step      = m_ctx.step;
            m_ctx.report.iteration = -1;
            m_ctx.report.note      = "step length stalled";
            return gsStatus::NotConverged;
        }
    }
}

template<class T>
gsStatus gsContinuationDriver<T>::solve(T lambdaEnd)
{
    m_ctx.report.reset();
    if (!m_setup)
    {
        m_ctx.report.status = gsStatus::OtherError;
        m_ctx.report.note   = "setup() not called";
        return gsStatus::OtherError;
    }
    if (!m_constraint->stepIsLoadIncrement())
    {
        m_ctx.report.status = gsStatus::OtherError;
        m_ctx.report.note   = "solve(lambdaEnd) needs a load-increment constraint";
        return gsStatus::OtherError;
    }
    // reachedEnd() reports a NaN load as reached, so the finiteness guard comes first
    if (_notFinite(m_ctx.lambda) || _notFinite(lambdaEnd) || _notFinite(m_controller.length()))
    {
        m_ctx.report.status   = gsStatus::NotConverged;
        m_ctx.report.failedOp = gsDriverOp::Controller;
        m_ctx.report.step     = m_ctx.step;
        m_ctx.report.note     = "non-finite lambda, lambdaEnd or step length";
        return gsStatus::NotConverged;
    }

    index_t taken = 0;
    while (!gsStepController<T>::reachedEnd(m_ctx.lambda, m_controller.length(), lambdaEnd))
    {
        const gsStatus st = _step(true, lambdaEnd);
        if (st != gsStatus::Success)
            return st;
        ++taken;
    }
    m_ctx.report.reset();
    m_ctx.report.status    = gsStatus::Success;
    m_ctx.report.step      = taken > 0 ? m_ctx.step - 1 : -1;
    m_ctx.report.iteration = taken > 0 ? m_ctx.prev.iterations : -1;
    return gsStatus::Success;
}

template<class T>
gsStatus gsContinuationDriver<T>::step(const gsPathPoint<T> & from, const gsStepHistory<T> & hist, T ds,
                                       gsPathPoint<T> & to, gsStepHistory<T> & histOut,
                                       gsSolverReport * report) const
{
    gsSolverReport rep;
    gsStatus st;
    m_committedK = 0; // the space and the tangent solver are about to move to @a from
    if (!m_setup)
    {
        rep.status = gsStatus::OtherError;
        rep.note   = "setup() not called";
        st = rep.status;
    }
    else if (!m_S->copyable())
    {
        rep.status   = gsStatus::OtherError;
        rep.failedOp = gsDriverOp::Import;
        rep.step     = 0;
        rep.note     = "pure step needs a copyable state space";
        st = rep.status;
    }
    else if (!m_S->importCommitted(from.U))
    {
        rep.status   = gsStatus::OtherError;
        rep.failedOp = gsDriverOp::Import;
        rep.step     = 0;
        rep.note     = "importCommitted failed";
        st = rep.status;
    }
    else
    {
        const index_t n = m_S->tangentSize();
        gsStepContext<T> ctx;
        ctx.F      = m_ctx.F;
        ctx.DeltaU = gsVector<T>::Zero(n);
        ctx.deltaU = gsVector<T>::Zero(n);
        ctx.dRdl   = ctx.F;
        ctx.lambda = from.lambda;
        ctx.prev   = hist;
        ctx.step   = 0;

        const bool useSecant = m_secant && _secantUsable(hist);
        st = _attempt(ctx, ds, useSecant);
        rep = ctx.report;
        if (st == gsStatus::Success)
        {
            if (m_S->exportWorking(to.U))
            {
                to.lambda = ctx.lambda;
                histOut   = ctx.prev;
            }
            else
            {
                rep.status   = gsStatus::OtherError;
                rep.failedOp = gsDriverOp::Export;
                rep.note     = "exportWorking failed";
                st = rep.status;
            }
        }
    }
    if (report)
        *report = rep;
    return st;
}

} // namespace gismo
