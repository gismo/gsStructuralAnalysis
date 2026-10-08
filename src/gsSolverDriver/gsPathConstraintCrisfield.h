 /** @file gsPathConstraintCrisfield.h

    @brief Spherical path constraints of the solver driver: Crisfield's
    constraint with a quadratic root choice, and the consistently linearised
    variant with a run-constant scaling.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#pragma once
#include <string>
#include <limits>

#include <gsCore/gsLinearAlgebra.h>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsPathConstraint.h>

namespace gismo
{

/**
   @brief Crisfield's spherical arc-length constraint with the quadratic
   solved for the load increment (Crisfield 1981; Lam and Morley 1992;
   Ritto-Correa and Camotim 2008; Feng et al. 1995).

   <b>Constraint.</b> With the scaling \f$\varphi\f$ and
   \f$\psi=\varphi^2\|F\|^2\f$ the increment \f$(\Delta U,\Delta\lambda)\f$
   of a step lies on the sphere
   \f[ g=\Delta U\!\cdot\!\Delta U+\psi\,\Delta\lambda^2-s^2=0 .\f]
   Warning: the load term is weighted by \f$\psi=\varphi^2\|F\|^2\f$ in the
   constraint and the predictor, but by \f$\varphi^2\f$ alone in the DOT
   product of the root choice.

   <b>Scaling.</b> setScaling() with \f$\varphi\neq-1\f$ fixes \f$\varphi\f$.
   With \f$\varphi=-1\f$ (automatic), \f$\varphi\f$ is recomputed at every
   predict(): on the first step \f$\varphi=\|u_t\|/\|F\|\f$, which places the
   predictor on the sphere, and on every other step the secant ratio
   \f$\varphi=\|U\|/(|\lambda|\,\|F\|)\f$ of the committed state, taken
   without a guard. It is therefore infinite or NaN at \f$\lambda=0\f$ and
   unbounded as \f$\lambda\to0\f$; such a value is not an error here and
   reaches the convergence tests.

   <b>Predictor.</b> First step (no usable previous increment): with
   \f$u_t=K^{-1}F\f$, \f$\Delta\lambda=s/\sqrt{2\,u_t\!\cdot\!u_t}\f$ and
   \f$\Delta U=\Delta\lambda\,u_t\f$. Otherwise the MU predictor
   \f$\Delta\lambda=\pm s/\sqrt{u_t\!\cdot\!u_t+\psi}\f$,
   \f$\Delta U=\Delta\lambda\,u_t\f$, which lies on the sphere; the sign is
   that of \f$\Delta U_{old}\!\cdot\!u_t+\psi\,\Delta\lambda_{old}\f$, the
   inner product of the previous increment with \f$(u_t,1)\f$ in the
   \f$\psi\f$-weighted metric.

   <b>Corrector.</b> With \f$\delta U=\eta\,\delta\bar U+\delta\lambda\,u_t\f$
   (\f$\eta=1\f$ unless modified) the constraint at the updated increment is
   the quadratic \f$a_0\delta\lambda^2+(b_0+\eta b_1)\delta\lambda
   +(c_0+\eta c_1+\eta^2c_2)=0\f$ with
   \f$a_0=u_t\!\cdot\!u_t+\psi\f$,
   \f$b_0=2(u_t\!\cdot\!\Delta U+\psi\,\Delta\lambda)\f$,
   \f$b_1=2\,\delta\bar U\!\cdot\!u_t\f$,
   \f$c_0=g\f$,
   \f$c_1=2\,\Delta U\!\cdot\!\delta\bar U\f$,
   \f$c_2=\delta\bar U\!\cdot\!\delta\bar U\f$.
   <ul>
   <li><b>Real roots.</b> The root whose update has the larger inner product
   with the previous increment, \f$\delta\lambda\,(\Delta U_{old}\!\cdot\!u_t
   +\varphi^2\Delta\lambda_{old})\f$ (Ritto-Correa and Camotim 2008), is taken.</li>
   <li><b>Complex roots, modified \f$\eta\f$ (Lam and Morley 1992).</b>
   \f$\eta\f$ is chosen so that the discriminant vanishes: the roots in
   \f$\eta\f$ of the resulting quadratic are moved by 5% of their distance to
   the side of 1 they lie on, and \f$\delta\lambda=-\alpha_2/(2\alpha_1)\f$
   (a double root). Used when such \f$\eta\f$ exist and \f$\eta>0.05\f$.</li>
   <li><b>Complex roots, simplified rescaling.</b> Otherwise the Newton point
   \f$(\Delta U+\delta\bar U,\lambda_{cr}-\lambda)\f$ is scaled onto the
   sphere. \f$\lambda_{cr}\f$ makes the out-of-balance force orthogonal to
   \f$F\f$, with the internal force at the current iterate
   \f$F_{int}=(\lambda+\Delta\lambda)F-R\f$. The formula assumes a dead load:
   ctx.F is the load direction at every state.</li>
   </ul>

   <b>Angle methods.</b> Root choice in correct(): the previous increment is
   the last accepted step (Step; zero if there is none), the current
   increment of the running step (Iteration), or the predictor increment of
   the running step (Predictor). MU sign in predict(): the last accepted
   step (Step and Iteration), or the predictor increment of the previous
   predict() call (Predictor; the last accepted step if there is none).

   <b>Relaxation.</b> A damping factor \f$\omega\neq1\f$ (setRelaxation())
   is applied to an update whose load increment changes sign and does not
   grow in magnitude relative to the previous one (an oscillating load
   factor): correct() then rebuilds \f$\delta U=\delta\lambda\,u_t+\eta\,
   \delta\bar U\f$ and overwrites ctx.relaxation with \f$\omega\f$, so that
   the caller scales both updates by \f$\omega\f$.

   Vector sizes: ctx.F, ctx.R, ctx.DeltaU, ctx.deltaU, dUt, dUbar have the
   tangent size \f$n\f$.

   \ingroup gsStructuralAnalysis
*/
template<class T>
class gsCrisfieldConstraint : public gsPathConstraint<T>
{
public:
    /// Angle methods, see setAngleMethod()
    struct angle { enum type { Step = 0, Iteration = 1, Predictor = 2 }; };

    /// Branch taken by the last correct()
    enum class Root { None, Real, Modified, Simplified };

    /// Automatic scaling, angle method Step, relaxation 1. O(1).
    gsCrisfieldConstraint()
    : m_phi(-1), m_phiUser(false), m_tangentScalingNext(false), m_angle(angle::Step), m_relax(1),
      m_FF(std::numeric_limits<T>::quiet_NaN()), m_DeltaLpred(0), m_last(Root::None),
      m_numModified(0), m_numSimplified(0), m_numRelaxed(0)
    { }

    /// Fixed scaling @a phi; -1 selects the automatic scaling. O(1).
    void setScaling(T phi) { m_phi = phi; m_phiUser = (phi != -1); }

    /// The user value, or -1 before the first automatic predict(), then the
    /// last computed value. O(1).
    T scaling() const { return m_phi; }

    /// Caches \f$F\cdot F\f$ so that distance() is defined before the first
    /// predict(); predict() refreshes it from ctx.F. O(n).
    void setForce(const gsVector<T> & F) { m_FF = F.dot(F); }

    /// The next predict() with automatic scaling takes
    /// \f$\varphi=\|u_t\|/\|F\|\f$ instead of the secant ratio (an initial
    /// guess); consumed by that predict(). O(1).
    void requestTangentScaling() { m_tangentScalingNext = true; }

    /// 0 = Step, 1 = Iteration, 2 = Predictor (GISMO_ENSURE on any other value). O(1).
    void setAngleMethod(index_t m);

    /// The angle method. O(1).
    index_t angleMethod() const { return m_angle; }

    /// Damping factor of the oscillation relaxation; 1 switches it off. O(1).
    void setRelaxation(T relax) { m_relax = relax; }

    /// The relaxation factor. O(1).
    T relaxation() const { return m_relax; }

    /// Branch of the last correct() (None before the first). O(1).
    Root lastRoot() const { return m_last; }

    /// Number of correct() calls that used the modified-\f$\eta\f$ roots. O(1).
    index_t numModifiedRoots() const { return m_numModified; }

    /// Number of correct() calls that used the simplified rescaling. O(1).
    index_t numSimplifiedRoots() const { return m_numSimplified; }

    /// Number of correct() calls that set ctx.relaxation. O(1).
    index_t numRelaxed() const { return m_numRelaxed; }

    /// Zero the three counters and reset lastRoot(). O(1).
    void resetCounters();

    /// Always true: both predictor branches use \f$u_t\f$.
    bool predictorNeedsTangent(const gsStepContext<T> & ctx) const;

    /// First-step or MU predictor of length ctx.stepSize (no rescaling by
    /// \f$s\f$ afterwards); caches \f$\|F\|^2\f$ for distance(). False iff dUt
    /// or ctx.F differ in size from space.tangentSize(), or if the automatic
    /// scaling needs the state magnitude and @a space provides neither
    /// magnitude() nor exportWorking(), or the vector it provides is not of size n. O(n).
    bool predict(gsStepContext<T> & ctx, const gsVector<T> & dUt,
                 const gsStateSpace<T> & space);

    /// Root choice of the quadratic for \f$\delta\lambda\f$, with the complex
    /// fallbacks and the conditional relaxation. False iff dUbar, dUt, ctx.F
    /// or ctx.R differ in size from ctx.DeltaU. Writes ctx.deltaU,
    /// ctx.deltaL and, when the relaxation fires, ctx.relaxation. O(n).
    bool correct(gsStepContext<T> & ctx, const gsVector<T> & dUbar,
                 const gsVector<T> & dUt);

    /// Signed \f$\Delta U\!\cdot\!\Delta U+\varphi^2\|F\|^2\Delta\lambda^2-s^2\f$. O(n).
    T residual(const gsStepContext<T> & ctx) const;

    /// \f$\sqrt{\Delta U\!\cdot\!\Delta U+\varphi^2\|F\|^2\Delta\lambda^2}\f$
    /// with the \f$\|F\|^2\f$ of the last predict(), or of setForce() before
    /// the first one (NaN if neither happened). O(n).
    T distance(const gsVector<T> & DeltaU, T DeltaL) const;

    /// False: ctx.stepSize is an arc length.
    bool stepIsLoadIncrement() const;

    /// "Crisfield"
    std::string name() const;

private:
    /// True iff there is no usable previous increment in ctx.prev (invalid,
    /// of a size different from ctx.F, or exactly zero)
    bool firstStep(const gsStepContext<T> & ctx) const;

private:
    T m_phi;                    ///< scaling \f$\varphi\f$
    bool m_phiUser;             ///< true iff \f$\varphi\f$ was set by the user
    bool m_tangentScalingNext;  ///< the next automatic predict() takes the tangent ratio
    index_t m_angle;            ///< angle method
    T m_relax;                  ///< relaxation factor
    T m_FF;                     ///< \f$F\!\cdot\!F\f$ of setForce() or of the last predict()
    gsVector<T> m_DeltaUpred;   ///< displacement increment of the last predictor
    T m_DeltaLpred;             ///< load increment of the last predictor
    Root m_last;                ///< branch of the last correct()
    index_t m_numModified;      ///< count of modified-\f$\eta\f$ corrections
    index_t m_numSimplified;    ///< count of simplified rescalings
    index_t m_numRelaxed;       ///< count of relaxed corrections
};

/**
   @brief Consistently linearised spherical constraint with a run-constant
   scaling (Schweizerhof and Wriggers, CMAME 59(3) 261-279, 1986).

   <b>Constraint.</b> \f$g=\Delta U\!\cdot\!\Delta U+\varphi^2\|F\|^2
   \Delta\lambda^2-s^2=0\f$.

   <b>Scaling.</b> With automatic scaling \f$\varphi\f$ is one constant per
   run, \f$\varphi_0=\|u_{t,0}\|/\|F\|\f$ with \f$u_{t,0}=K^{-1}F\f$ at the
   first step. The load term \f$\varphi\|F\|\Delta\lambda=\|u_{t,0}\|
   \Delta\lambda\f$ then has the units of a displacement and the weight
   \f$\psi=\varphi^2\|F\|^2\f$ is the same at every step. A scaling that
   follows the path would change the metric from step to step: the tangent
   ratio \f$\|u_t\|/\|F\|\f$ is unbounded at a limit point, where \f$K\f$ is
   singular, and the secant ratio \f$\|U\|/(|\lambda|\|F\|)\f$ is unbounded as
   \f$\lambda\to0\f$ at a regular point. A run that does not start from a first
   step (a restart) freezes \f$\varphi_0\f$ at its own first step from a freshly
   solved \f$u_t\f$. \f$\varphi_0\f$ is assigned once: later first steps
   (redone steps), setScaling() and later steps leave it, and only
   resetFrozenScaling() starts a new run.

   <b>Predictor.</b> First step: \f$(\Delta U,\Delta\lambda)=(u_t,1)/\|u_t\|
   \cdot s\f$. Otherwise the secant: the previous increment rescaled by
   \f$s/s_{prev}\f$. The predicted increment does not depend on
   \f$\varphi\f$.

   <b>Corrector (consistent linearisation).</b> The linearisation
   \f$g+\partial g/\partial\Delta U\cdot\delta U+\partial g/\partial\Delta
   \lambda\,\delta\lambda=0\f$ with \f$\delta U=\delta\bar U+\delta\lambda\,
   u_t\f$ gives in closed form
   \f[ \delta\lambda=-\frac{g+2\,\Delta U\!\cdot\!\delta\bar U}
       {2\,\Delta U\!\cdot\!u_t+2\varphi^2\|F\|^2\Delta\lambda}. \f]
   A zero or non-finite denominator is not a failure; the NaN or infinite
   update reaches the convergence tests.

   Vector sizes: ctx.F, ctx.DeltaU, dUt, dUbar have the tangent size
   \f$n\f$.

   \ingroup gsStructuralAnalysis
*/
template<class T>
class gsConsistentCrisfieldConstraint : public gsPathConstraint<T>
{
public:
    /// Automatic scaling, no frozen \f$\varphi_0\f$. O(1).
    gsConsistentCrisfieldConstraint()
    : m_phi(-1), m_phiUser(false), m_phi0(0), m_phi0Set(false),
      m_FF(std::numeric_limits<T>::quiet_NaN())
    { }

    /// Fixed scaling @a phi; -1 selects the automatic scaling, which keeps a
    /// frozen \f$\varphi_0\f$ if there is one. Does not touch the frozen
    /// \f$\varphi_0\f$. O(1).
    void setScaling(T phi)
    {
        m_phiUser = (phi != -1);
        m_phi = (!m_phiUser && m_phi0Set) ? m_phi0 : phi;
    }

    /// The user value; with automatic scaling -1 while no \f$\varphi_0\f$ is
    /// frozen, then the frozen \f$\varphi_0\f$ (after resetFrozenScaling() the
    /// last value until the next predict()). O(1).
    T scaling() const { return m_phi; }

    /// Caches \f$F\cdot F\f$ so that distance() is defined before the first
    /// predict(); predict() refreshes it from ctx.F. O(n).
    void setForce(const gsVector<T> & F) { m_FF = F.dot(F); }

    /// Freezes \f$\varphi_0\f$ to @a phi0, which must be positive and finite
    /// (GISMO_ENSURE; NaN is rejected). The automatic scaling then uses it
    /// instead of freezing its own at the first predict(); a user scaling is
    /// not affected. O(1).
    void setFrozenScaling(T phi0)
    {
        GISMO_ENSURE(phi0 > 0 && phi0 < std::numeric_limits<T>::infinity(),
                     "gsConsistentCrisfieldConstraint::setFrozenScaling: the scaling must be positive and finite, got " << phi0);
        m_phi0 = phi0;
        m_phi0Set = true;
        if (!m_phiUser)
            m_phi = phi0;
    }

    /// Forget the frozen \f$\varphi_0\f$: the next automatic first step freezes a new one. O(1).
    void resetFrozenScaling() { m_phi0 = 0; m_phi0Set = false; }

    /// The frozen \f$\varphi_0\f$ (0 if none). O(1).
    T frozenScaling() const { return m_phi0; }

    /// True iff \f$\varphi_0\f$ has been frozen. O(1).
    bool hasFrozenScaling() const { return m_phi0Set; }

    /// True iff the first-step predictor is used, or the automatic scaling
    /// still has to freeze \f$\varphi_0\f$ from \f$u_t\f$. Reads ctx.prev and
    /// ctx.F.size(). O(n).
    bool predictorNeedsTangent(const gsStepContext<T> & ctx) const;

    /// First-step or secant predictor of length ctx.stepSize; freezes
    /// \f$\varphi_0\f$ when automatic and not yet frozen; caches
    /// \f$\|F\|^2\f$ for distance(). False iff ctx.F differs in size from
    /// space.tangentSize(), or dUt does while it is needed. O(n).
    bool predict(gsStepContext<T> & ctx, const gsVector<T> & dUt,
                 const gsStateSpace<T> & space);

    /// Closed-form consistently linearised update. False iff dUbar, dUt or
    /// ctx.F differ in size from ctx.DeltaU. Writes ctx.deltaU and ctx.deltaL
    /// only. O(n).
    bool correct(gsStepContext<T> & ctx, const gsVector<T> & dUbar,
                 const gsVector<T> & dUt);

    /// Signed \f$\Delta U\!\cdot\!\Delta U+\varphi^2\|F\|^2\Delta\lambda^2-s^2\f$. O(n).
    T residual(const gsStepContext<T> & ctx) const;

    /// \f$\sqrt{\Delta U\!\cdot\!\Delta U+\varphi^2\|F\|^2\Delta\lambda^2}\f$
    /// with the \f$\|F\|^2\f$ of the last predict(), or of setForce() before
    /// the first one (NaN if neither happened). O(n).
    T distance(const gsVector<T> & DeltaU, T DeltaL) const;

    /// False: ctx.stepSize is an arc length.
    bool stepIsLoadIncrement() const;

    /// "ConsistentCrisfield"
    std::string name() const;

private:
    /// True iff there is no usable secant in ctx.prev
    bool firstStep(const gsStepContext<T> & ctx) const;

private:
    T m_phi;        ///< scaling \f$\varphi\f$
    bool m_phiUser; ///< true iff \f$\varphi\f$ was set by the user
    T m_phi0;       ///< frozen \f$\varphi_0\f$
    bool m_phi0Set; ///< true iff \f$\varphi_0\f$ is frozen
    T m_FF;         ///< \f$F\!\cdot\!F\f$ of setForce() or of the last predict()
};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsPathConstraintCrisfield.hpp)
#endif
