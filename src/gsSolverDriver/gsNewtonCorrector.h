 /** @file gsNewtonCorrector.h

    @brief Newton corrector policy of the continuation driver: refactorization
    schedule, relaxation and application of one update to the state space.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsSolverDriverTypes.h>

namespace gismo
{

/**
   @brief Policy of the (quasi-)Newton corrector: when to refactorize, the
   relaxation factor, and the application of one update.

   The class holds no matrix, solver, problem or constraint, and does not test
   convergence or bound the iteration count (that belongs to the convergence
   tests). The driver runs the following sequence for each corrector update
   \f$k = \f$ ctx.iteration + 1, with the residual
   \f$R = \lambda F - F_{int}\f$ and \f$K = \partial F_{int}/\partial U\f$:

   \code
   if corrector.refactorize(ctx):  problem.jacobian(lambda+DeltaL, ctx.K); ++ctx.Kversion; tangent solver factorizes ctx.K
   dUbar = K^{-1} R          // K delta = R, not -R
   dUt   = K^{-1} dRdl       // re-solved every update: one back-substitution, bitwise equal to a cached value
   corrector.prepare(ctx)                          // ctx.relaxation = omega
   constraint.correct(ctx, dUbar, dUt)             // writes ctx.deltaU, ctx.deltaL; may overwrite ctx.relaxation
   corrector.apply(ctx, space)                     // relaxed apply + accumulate + ++iteration
   evaluate R and constraintResidual at the new working point; run the convergence test
   \endcode

   <b>Refactorization schedule.</b> With \f$q\f$ = quasi and interval \f$m\f$,
   K is assembled and factorized before update \f$k\f$ iff
   \f$ \lnot q \;\lor\; k \le 1 \;\lor\; (m > 0 \land k \bmod m = 0)\f$.
   The first update always uses a factorization at the predicted point.

   <b>Relaxation.</b> \f$\omega\f$ is a generic, unconditional factor.
   prepare() writes ctx.relaxation = \f$\omega\f$ before every constraint
   correction. A path constraint may overwrite ctx.relaxation, but never
   multiplies it. apply() scales the update by whatever ctx.relaxation holds
   when it is called. A constraint that overwrites ctx.relaxation is
   responsible for any rebuild of ctx.deltaU that its rule needs.

   apply() is the only place where ctx.iteration is incremented during the
   corrector; the driver must not increment it. Refreshing ctx.R and
   ctx.constraintResidual after apply() is the driver's job.

   Every operation is O(1) except apply(), which is O(n) plus one
   space.apply().

   \tparam T coefficient type

   \ingroup gsStructuralAnalysis
*/
template<class T>
class gsNewtonCorrector
{
public:
    /// Default constructor: relaxation 1, quasi-Newton off, interval -1
    gsNewtonCorrector()
    : m_relaxation(T(1)), m_quasi(false), m_interval(-1)
    {}

    /// Sets the relaxation factor \f$\omega\f$ (default 1; any value accepted)
    void setRelaxation(T omega) { m_relaxation = omega; }

    /// Relaxation factor \f$\omega\f$
    T relaxation() const { return m_relaxation; }

    /// Sets quasi-Newton mode and the refactorization interval; a
    /// non-positive interval means: factorize only before the first update
    void setQuasi(bool quasi, index_t interval = -1)
    {
        m_quasi    = quasi;
        m_interval = interval;
    }

    /// Whether the factorization is reused between updates
    bool quasi() const { return m_quasi; }

    /// Refactorization interval of quasi-Newton mode
    index_t quasiInterval() const { return m_interval; }

    /**
       @brief The refactorization rule for the k-th update (1-based).

       True iff \f$ \lnot\mathrm{quasi} \lor k \le 1 \lor
       (\mathrm{interval} > 0 \land k \bmod \mathrm{interval} = 0)\f$; true
       for \f$k \le 0\f$. O(1).
    */
    static bool refactorizeAt(index_t k, bool quasi, index_t interval)
    {
        return !quasi || k <= 1 || (interval > 0 && k % interval == 0);
    }

    /// Whether K is to be reassembled and factorized before the next update,
    /// i.e. refactorizeAt(ctx.iteration + 1, quasi(), quasiInterval()). O(1).
    bool refactorize(const gsStepContext<T>& ctx) const
    {
        return refactorizeAt(ctx.iteration + 1, m_quasi, m_interval);
    }

    /// Writes ctx.relaxation = \f$\omega\f$ and nothing else. Call immediately
    /// before the path constraint's correction. O(1).
    void prepare(gsStepContext<T>& ctx) const
    {
        ctx.relaxation = m_relaxation;
    }

    /**
       @brief Applies one corrector update.

       Moves the working state by ctx.relaxation * ctx.deltaU, scales ctx.deltaU
       and ctx.deltaL by ctx.relaxation (only when it differs from 1), so that
       they hold the applied increment, adds them to ctx.DeltaU and ctx.DeltaL
       and increments ctx.iteration.

       Returns false with ctx untouched when ctx.deltaU or ctx.DeltaU does not
       have space.tangentSize() entries, or when space.apply() fails; the
       driver then rolls the step back.

       Complexity: O(n) plus the cost of one space.apply(), n = tangentSize().
    */
    bool apply(gsStepContext<T>& ctx, gsStateSpace<T>& space) const
    {
        const index_t n = space.tangentSize();
        if (ctx.deltaU.size() != n || ctx.DeltaU.size() != n)
            return false;
        if (!space.apply(ctx.deltaU, ctx.relaxation))
            return false;
        if (ctx.relaxation != T(1))
        {
            ctx.deltaU *= ctx.relaxation;
            ctx.deltaL *= ctx.relaxation;
        }
        ctx.DeltaU += ctx.deltaU;
        ctx.DeltaL += ctx.deltaL;
        ++ctx.iteration;
        return true;
    }

private:
    T m_relaxation;
    bool m_quasi;
    index_t m_interval;
};

} // namespace gismo
