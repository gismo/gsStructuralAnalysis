 /** @file gsStepController.h

    @brief Step-length policy of the continuation driver: adaptive growth,
    halving on rejection, regrow, failure budget, end clamp and stall guards.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#pragma once
#include <cmath>

#include <gsCore/gsLinearAlgebra.h>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsSolverDriverTypes.h>

namespace gismo
{

/**
   @brief Owns the step length \f$s\f$ of a continuation run.

   Three lengths are kept: the current length \f$s\f$, the length of the
   previous step \f$s_{prev}\f$ (the secant predictor scales its extrapolation
   with it) and the original length \f$s_{ori}\f$ to which resetLength()
   returns.

   <b>Adaptive length.</b> After an accepted step that needed \f$k\f$
   corrector updates, with desired count \f$N_d\f$,
   \f[
      s_{prev} \leftarrow s,\qquad
      s \leftarrow s\cdot\min\!\big(2,\max(\tfrac12, N_d/k)\big),
   \f]
   and \f$N_d/k := 1\f$ for \f$k = 0\f$. With adaptive length off,
   \f$s_{prev} \leftarrow s\f$ only.

   <b>Retry policy.</b> rejected() halves \f$s\f$ while the net bisection depth
   stays below MaxBisect and the number of failures below MaxFailures. With
   regrow on, an accepted step doubles \f$s\f$ again and lowers the depth.

   <b>Bound.</b> With adaptive length off and no reduceLength() call after the
   last setLength(), \f$s = s_{ori}\,2^{-\mathrm{depth}} \ge
   s_{ori}/2^{\mathrm{MaxBisect}}\f$; regrow keeps this identity because it
   doubles \f$s\f$ and decrements the depth together. With adaptive length on,
   the factor 1/2 of the growth clamp shrinks \f$s\f$ on accepted steps without
   any failure, and reduceLength() changes \f$s\f$ without touching the depth;
   in those cases the bound holds for the halvings only.

   <b>Default configuration.</b> MaxBisect 0, Regrow off and MaxFailures -1
   (unlimited) make rejected() always return false and the regrow rule
   unreachable: the length then changes only through setLength(),
   reduceLength(), resetLength() and the adaptive update.

   Defaults: \f$s = s_{prev} = s_{ori} = 10^{-2}\f$, adaptive off,
   \f$N_d = 10\f$. No value is validated: a non-positive \f$N_d\f$ yields a
   non-positive factor, which the growth clamp turns into 1/2, and the sign of
   the length is free (unloading uses a negative load increment).

   The controller assembles, solves and tests nothing. Every operation is
   O(1).

   \tparam T coefficient type

   \ingroup gsStructuralAnalysis
*/
template<class T>
class gsStepController
{
public:
    /// Default constructor: length 1e-2, adaptive off, 10 desired iterations,
    /// MaxBisect 0, Regrow off, MaxFailures -1
    gsStepController()
    : m_length(T(1e-2)), m_lengthPrev(T(1e-2)), m_lengthOri(T(1e-2)),
      m_adaptive(false), m_desired(10),
      m_maxBisect(0), m_regrow(false), m_maxFailures(-1),
      m_depth(0), m_failures(0)
    {}

    /// Constructor with initial length, see setLength()
    explicit gsStepController(T length)
    : m_adaptive(false), m_desired(10),
      m_maxBisect(0), m_regrow(false), m_maxFailures(-1),
      m_depth(0), m_failures(0)
    {
        setLength(length);
    }

    /// Current length \f$s\f$
    T length() const { return m_length; }

    /// Length \f$s_{prev}\f$ of the previous accepted step
    T previousLength() const { return m_lengthPrev; }

    /// Original length \f$s_{ori}\f$
    T originalLength() const { return m_lengthOri; }

    /// Sets \f$s = s_{prev} = s_{ori} = l\f$ and clears depth and failures
    void setLength(T l)
    {
        m_length = m_lengthPrev = m_lengthOri = l;
        resetCounters();
    }

    /// Sets the adaptive flag and the desired iteration count \f$N_d\f$
    void setAdaptive(bool adaptive, index_t desired)
    {
        m_adaptive = adaptive;
        m_desired  = desired;
    }

    /// Whether the length adapts to the iteration count of accepted steps
    bool adaptive() const { return m_adaptive; }

    /// Desired iteration count \f$N_d\f$
    index_t desiredIterations() const { return m_desired; }

    /// Maximal net bisection depth (default 0: no retry)
    void setMaxBisect(index_t n) { m_maxBisect = n; }
    /// Maximal net bisection depth
    index_t maxBisect() const { return m_maxBisect; }

    /// Whether an accepted step undoes one halving (default false)
    void setRegrow(bool r) { m_regrow = r; }
    /// Whether an accepted step undoes one halving
    bool regrow() const { return m_regrow; }

    /// Maximal number of rejections granted, counted over the run since the
    /// last setLength()/resetCounters(); negative means unlimited (default -1)
    void setMaxFailures(index_t n) { m_maxFailures = n; }
    /// Failure budget; negative means unlimited
    index_t maxFailures() const { return m_maxFailures; }

    /// Net bisection depth: halvings minus regrows
    index_t depth() const { return m_depth; }

    /// Number of rejections granted since the last reset
    index_t failures() const { return m_failures; }

    /// Sets depth and failures to zero
    void resetCounters() { m_depth = m_failures = 0; }

    /**
       @brief Growth factor of the adaptive update.

       \f$ N_d/k \f$ for \f$k > 0\f$ and 1 for \f$k = 0\f$, clamped to
       \f$[1/2, 2]\f$.

       @param[in] desired    desired iteration count \f$N_d\f$
       @param[in] iterations corrector updates \f$k\f$ of the accepted step
    */
    static T adaptiveFactor(index_t desired, index_t iterations)
    {
        T fac = (iterations > 0)
              ? static_cast<T>(desired) / static_cast<T>(iterations)
              : T(1);
        if (fac < 0.5)
            fac = 0.5;
        else if (fac > 2.0)
            fac = 2.0;
        return fac;
    }

    /// Adaptive update: \f$s_{prev} = s\f$, \f$s = s\cdot\f$adaptiveFactor().
    /// Independent of the adaptive flag.
    void computeLength(index_t iterations)
    {
        const T fac = adaptiveFactor(m_desired, iterations);
        m_lengthPrev = m_length;
        m_length = m_length*fac;
    }

    /// Multiplies \f$s\f$ by fac and returns it. \f$s_{prev}\f$, \f$s_{ori}\f$,
    /// depth and failures are unchanged.
    T reduceLength(T fac)
    {
        m_length *= fac;
        return m_length;
    }

    /**
       @brief Returns to the original length.

       Sets \f$s = s_{prev} = s_{ori}\f$ and, when adaptive, applies
       computeLength() with the iteration count passed by the caller. Depth
       and failures are unchanged.

       @param[in] iterations the caller's current iteration count (the count
       of the last step, or 0 when there is none)
    */
    T resetLength(index_t iterations)
    {
        m_length = m_lengthPrev = m_lengthOri;
        if (m_adaptive)
            computeLength(iterations);
        return m_length;
    }

    /**
       @brief Registers an accepted step.

       With adaptive length computeLength(), else \f$s_{prev} = s\f$. Then,
       with regrow on and a positive depth, doubles \f$s\f$ and decrements the
       depth.

       @param[in] iterations corrector updates of the accepted step
    */
    void accepted(index_t iterations)
    {
        if (m_adaptive)
            computeLength(iterations);
        else
            m_lengthPrev = m_length;

        if (m_regrow && m_depth > 0)
        {
            m_length = m_length*2;
            --m_depth;
        }
    }

    /**
       @brief Registers a rejected step and halves the length when a retry is
       granted.

       Returns false, leaving everything unchanged, when the failure budget is
       used up (MaxFailures non-negative and reached) or the depth has reached
       MaxBisect. Otherwise increments depth and failures, halves \f$s\f$ and
       returns true. \f$s_{prev}\f$ is not touched.
    */
    bool rejected()
    {
        if (m_maxFailures >= 0 && m_failures >= m_maxFailures)
            return false;
        if (m_depth >= m_maxBisect)
            return false;
        ++m_depth;
        ++m_failures;
        m_length = m_length/2;
        return true;
    }

    /**
       @brief Load value targeted by a step of increment ds, clamped at the end
       of the load interval.

       Returns \f$\lambda + ds\f$, replaced by lambdaEnd when it would pass
       lambdaEnd in the direction of ds. It is a target load, not a step. A
       load constraint applying \f$\Delta\lambda = \mathrm{target} - \lambda\f$
       reaches the target exactly whenever that subtraction is exact (Sterbenz:
       same sign and \f$\lambda \in [\mathrm{target}/2, 2\,\mathrm{target}]\f$,
       the usual case on the last clamped step); otherwise it may miss by one
       ulp.
    */
    static T targetLoad(T lambda, T ds, T lambdaEnd)
    {
        T lnew = lambda + ds;
        if ( (ds>0 && lnew>lambdaEnd) || (ds<0 && lnew<lambdaEnd) )
            lnew = lambdaEnd;
        return lnew;
    }

    /// True when lambda has reached lambdaEnd within the relative tolerance
    /// \f$10^{-12}|\lambda_{end}|\f$ in the direction of ds, or when ds = 0
    static bool reachedEnd(T lambda, T ds, T lambdaEnd)
    {
        return !( (ds>0 && lambda < lambdaEnd - T(1e-12)*std::abs(lambdaEnd)) ||
                  (ds<0 && lambda > lambdaEnd + T(1e-12)*std::abs(lambdaEnd)) );
    }

    /// True when the target load equals lambda (the step is below the
    /// resolution of lambda) or either value is non-finite
    static bool loadStalled(T lambda, T lambdaTarget)
    {
        return !std::isfinite(lambda) || !std::isfinite(lambdaTarget)
               || lambdaTarget == lambda;
    }

    /// True when the length is zero or non-finite: the stall guard when
    /// lambda is not controlled (arc length)
    static bool lengthStalled(T ds)
    {
        return ds == T(0) || !std::isfinite(ds);
    }

private:
    T m_length, m_lengthPrev, m_lengthOri;
    bool m_adaptive;
    index_t m_desired;
    index_t m_maxBisect;
    bool m_regrow;
    index_t m_maxFailures;
    index_t m_depth, m_failures;
};

} // namespace gismo
