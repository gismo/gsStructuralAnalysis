 /** @file gsPathConstraint.h

    @brief Path constraints of the solver driver: the predictor, the choice of
    the load increment from the two bordered solves, and the constraint
    residual of one continuation method.

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
#include <gsStructuralAnalysis/src/gsSolverDriver/gsStateSpace.h>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsSolverDriverTypes.h>

namespace gismo
{

/**
   @brief Interface of a path constraint.

   A continuation step solves \f$R(U,\lambda)=0\f$ together with a scalar
   constraint \f$g(\Delta U,\Delta\lambda;s)=0\f$ for the increments
   \f$(\Delta U,\Delta\lambda)\f$ from the last committed point. The residual
   is \f$R = \lambda F - F_{int}\f$ and the tangent \f$K=\partial F_{int}/\partial U\f$,
   so a Newton correction solves \f$K\,\delta = R\f$. Linearising the
   residual and the constraint about the current iterate gives the bordered
   system
   \f[
     \begin{pmatrix} K & -\partial R/\partial\lambda \\ \nabla_U g^T & \partial_\lambda g \end{pmatrix}
     \begin{pmatrix}\delta U\\ \delta\lambda\end{pmatrix}
     = \begin{pmatrix} R \\ -g \end{pmatrix}.
   \f]
   Eliminating \f$\delta U = \delta\bar U + \delta\lambda\, u_t\f$ with
   \f$\delta\bar U = K^{-1}R\f$ and \f$u_t = K^{-1}\partial R/\partial\lambda\f$
   leaves one scalar equation for \f$\delta\lambda\f$. A constraint owns that
   scalar equation (correct()), the predictor (predict()) and the value of
   \f$g\f$ (residual()); it does not factorize, assemble or touch the state.

   <b>Call protocol</b> (the driver follows it exactly). For each step
   attempt, with the state rolled back to the committed point:
   @code
   ctx.F set (tangent size n) before the first query; ctx.prev = the secant data
   ctx.resetStep(); ctx.stepSize = s;                // DeltaU/deltaU keep size n
   if constraint.predictorNeedsTangent(ctx):
       K at the committed state; ++ctx.Kversion; factorize
       ctx.dRdl at the committed state (F if the problem has no loadDirection)
       solve K dUt = ctx.dRdl                        // dUt of size n
   else: dUt is empty and is not read
   constraint.predict(ctx, dUt, space)   // writes DeltaU, DeltaL, deltaU = DeltaU, deltaL = DeltaL
   space.apply(ctx.DeltaU, 1)            // the driver applies the predictor, never the constraint
   R at lambda + DeltaL; ctx.constraintResidual = constraint.residual(ctx)
   convergence tests: initStep, test (iteration 0)
   per corrector update: [refactorize] -> dUbar = K^{-1} R, dUt = K^{-1} dRdl
       -> constraint.correct(ctx, dUbar, dUt) -> apply (relaxed) update
       -> R; ctx.constraintResidual = constraint.residual(ctx); test
   @endcode

   <b>Contract.</b>
   -# <em>Where the state is touched.</em> predict() reads the state only
      through the const @a space (tangentSize(), and for constraints that need
      it magnitude()/exportWorking()). The driver applies ctx.DeltaU with
      \f$s=1\f$ right after a successful predict(), so a failed apply is
      attributed to the apply operation and a false predict() to the
      constraint.
   -# <em>Same context for the query and the predict.</em>
      predictorNeedsTangent() and predict() see the same @a ctx (in
      particular the same ctx.prev), so the query decides the branch that
      predict() takes.
   -# <em>No stale factorization.</em> Skipping the factorization at the
      committed state when the predictor does not need it is safe, because
      the corrector refactorizes before its first update.
   -# <em>The step history and the secant.</em> ctx.prev is the secant data
      \f$\{\Delta U_{prev},\Delta\lambda_{prev},s_{prev},\text{valid}\}\f$. A
      secant predictor scales \f$(\Delta U_{prev},\Delta\lambda_{prev})\f$ by
      \f$s/s_{prev}\f$, so ctx.prev.stepSize is the length by which the
      previous increment is normalised; the driver stores the accepted step's
      \f$s\f$ there at commit. A caller that changes the step length between
      steps, or supplies its own previous point, sets it consistently.
      ctx.prev.valid false, or a ctx.prev.DeltaU whose size differs from
      \f$n\f$, means "no secant information": constraints with a secant fall
      back to their first-step predictor.
   -# <em>What correct() writes.</em> ctx.deltaU (size \f$n\f$) and
      ctx.deltaL only. It never touches DeltaU, DeltaL, iteration, lambda,
      prev, R, K or constraintResidual. The corrector accumulates. A
      constraint may overwrite ctx.relaxation but never multiplies it.
   -# <em>What predict() writes.</em> ctx.DeltaU (resized to \f$n\f$),
      ctx.DeltaL, ctx.deltaU = ctx.DeltaU and ctx.deltaL = ctx.DeltaL. It does
      not touch iteration, lambda, prev, R, K, Kversion, relaxation or
      constraintResidual.
   -# <em>Failure semantics.</em> predict() and correct() return false only on
      a structural error: a vector size mismatch, a missing dUt when
      predictorNeedsTangent() was true, or an invalid configuration. A NaN or
      infinite value computed from finite-size inputs is not a failure: it is
      written to @a ctx and left to the convergence tests, which must be able
      to iterate a non-finite iterate to their iteration limit. The only
      numerical exception is the exact zero pivot of the displacement
      constraint (gsDisplacementConstraint), which is reported as failure.
   -# <em>The residual and the convergence test.</em> residual() returns the
      signed, un-normalised \f$g\f$. The path-constraint test measures
      \f$|g|/s^2\f$ (0 when \f$g\f$ is exactly 0); each class states what
      that measure means for it.

   predict() and correct() are non-const because constraints may keep state.
   The driver holds a constraint through Ptr, so a const driver method can
   still call them.

   @tparam T coefficient type

   \ingroup gsStructuralAnalysis
*/
template<class T>
class gsPathConstraint
{
public:
    typedef memory::shared_ptr<gsPathConstraint> Ptr;

    virtual ~gsPathConstraint() {}

    /// True iff predict() needs the load tangent \f$u_t=K^{-1}\partial R/\partial\lambda\f$
    /// at the committed state. Same @a ctx as the following predict().
    virtual bool predictorNeedsTangent(const gsStepContext<T> & ctx) const = 0;

    /// Predictor: writes ctx.DeltaU, ctx.DeltaL, ctx.deltaU = ctx.DeltaU and
    /// ctx.deltaL = ctx.DeltaL. Does not touch the state.
    /// @param[in,out] ctx step context; ctx.stepSize is the step size
    /// @param[in] dUt \f$u_t\f$, size \f$n\f$; empty and unread if predictorNeedsTangent() was false
    /// @param[in] space state space, read-only
    /// @return false only on a structural error (see the class contract)
    virtual bool predict(gsStepContext<T> & ctx, const gsVector<T> & dUt,
                         const gsStateSpace<T> & space) = 0;

    /// Corrector: from \f$\delta\bar U=K^{-1}R\f$ and \f$u_t\f$ writes
    /// ctx.deltaU and ctx.deltaL.
    /// @param[in,out] ctx step context
    /// @param[in] dUbar \f$K^{-1}R\f$, size \f$n\f$
    /// @param[in] dUt \f$K^{-1}\partial R/\partial\lambda\f$, size \f$n\f$
    /// @return false only on a structural error (see the class contract)
    virtual bool correct(gsStepContext<T> & ctx, const gsVector<T> & dUbar,
                         const gsVector<T> & dUt) = 0;

    /// Signed constraint value \f$g(\Delta U,\Delta\lambda;s)\f$ at
    /// (ctx.DeltaU, ctx.DeltaL) with \f$s\f$ = ctx.stepSize; 0 on the constraint.
    virtual T residual(const gsStepContext<T> & ctx) const = 0;

    /// Path distance of an increment (@a DeltaU, @a DeltaL) in the
    /// constraint's own metric.
    virtual T distance(const gsVector<T> & DeltaU, T DeltaL) const = 0;

    /// True iff ctx.stepSize is a load-factor increment (false: arc length or displacement).
    virtual bool stepIsLoadIncrement() const = 0;

    /// Short identifier for verbose output, e.g. "Riks".
    virtual std::string name() const = 0;
};

/**
   @brief Fixed load: \f$\lambda\f$ is held at its committed value and the
   step is a plain Newton iteration at fixed load.

   The predictor is the committed state itself (\f$\Delta U=0\f$,
   \f$\Delta\lambda=0\f$, no extrapolation) and does not need the tangent.
   ctx.stepSize is ignored. A relative-update test whose basis is the
   predictor increment (\f$\|\Delta U\|=0\f$ here) cannot converge, so a
   fixed-load preset uses its own convergence tests.

   The bordered system degenerates: \f$\delta\lambda=0\f$ and
   \f$\delta U=K^{-1}R\f$. \f$g\equiv 0\f$.

   \ingroup gsStructuralAnalysis
*/
template<class T>
class gsFixedConstraint : public gsPathConstraint<T>
{
public:
    /// Always false.
    bool predictorNeedsTangent(const gsStepContext<T> & ctx) const;

    /// Zero predictor of size space.tangentSize(); @a dUt is not read. O(n).
    bool predict(gsStepContext<T> & ctx, const gsVector<T> & dUt,
                 const gsStateSpace<T> & space);

    /// deltaU = dUbar, deltaL = 0; false iff dUbar.size() differs from ctx.DeltaU.size(). O(n).
    bool correct(gsStepContext<T> & ctx, const gsVector<T> & dUbar,
                 const gsVector<T> & dUt);

    /// Always 0. O(1).
    T residual(const gsStepContext<T> & ctx) const;

    /// \f$\|\Delta U\|\f$: \f$\lambda\f$ does not move, so the displacement increment is the only path measure. O(n).
    T distance(const gsVector<T> & DeltaU, T DeltaL) const;

    /// False: ctx.stepSize is ignored.
    bool stepIsLoadIncrement() const;

    /// "Fixed"
    std::string name() const;
};

/**
   @brief Load control: the load increment \f$\Delta\lambda=s\f$ is prescribed.

   With \f$\Delta\lambda\f$ prescribed the bordered system degenerates:
   \f$\delta\lambda=0\f$ and \f$\delta U=K^{-1}R\f$. The predictor is the
   Euler tangent \f$\Delta U=s\,K^{-1}\partial R/\partial\lambda\f$.

   <b>The constraint is exactly satisfied.</b> correct() never changes
   \f$\Delta\lambda\f$, and predict() stores \f$\Delta\lambda\f$ as the
   target by assignment, so \f$g=\Delta\lambda-\Delta\lambda_{target}\f$
   (relative to the target if it is non-zero) is exactly 0 after predict() and
   after every update. The path-constraint test's \f$|g|/s^2\f$ therefore
   agrees with the relative measure
   \f$|\Delta\lambda-\Delta\lambda_{target}|/|\Delta\lambda_{target}|\f$
   (absolute for a zero target): the two scalings differ in general but
   coincide at 0, the only value \f$g\f$ takes. A zero target (\f$s=0\f$)
   gives \f$g=0\f$, never 0/0.

   \ingroup gsStructuralAnalysis
*/
template<class T>
class gsLoadConstraint : public gsPathConstraint<T>
{
public:
    /// Default constructor; the stored load-increment target is 0 until the first predict(). O(1).
    gsLoadConstraint() : m_target(0) {}

    /// Always true.
    bool predictorNeedsTangent(const gsStepContext<T> & ctx) const;

    /// \f$\Delta\lambda=\delta\lambda=s\f$, \f$\Delta U=\delta U=s\,u_t\f$,
    /// and the stored target is set to \f$s\f$. False iff dUt.size() differs
    /// from space.tangentSize(). O(n).
    bool predict(gsStepContext<T> & ctx, const gsVector<T> & dUt,
                 const gsStateSpace<T> & space);

    /// deltaU = dUbar, deltaL = 0; false iff dUbar.size() differs from ctx.DeltaU.size(). O(n).
    bool correct(gsStepContext<T> & ctx, const gsVector<T> & dUbar,
                 const gsVector<T> & dUt);

    /// Signed \f$\Delta\lambda-\Delta\lambda_{target}\f$, divided by
    /// \f$|\Delta\lambda_{target}|\f$ if the target is non-zero. O(1).
    T residual(const gsStepContext<T> & ctx) const;

    /// \f$|\Delta\lambda|\f$, the load-increment length. O(1).
    T distance(const gsVector<T> & DeltaU, T DeltaL) const;

    /// True.
    bool stepIsLoadIncrement() const;

    /// "Load"
    std::string name() const;

private:
    T m_target; ///< load increment prescribed by the last predict()
};

/**
   @brief Displacement control of one tangent dof \f$k\f$ (Batoz and Dhatt,
   Int. J. Numer. Meth. Eng. 14 (1979)): the increment \f$\Delta U_k=s\f$ is
   prescribed and \f$\lambda\f$ is the unknown.

   The constraint \f$g=\Delta U_k-s\f$ is linear, so the bordered system gives
   \f$\delta\lambda=-(g+\delta\bar U_k)/u_{t,k}\f$ and
   \f$\delta U=\delta\bar U+\delta\lambda\,u_t\f$, and after every update
   \f$\Delta U_k=s\f$ up to rounding. The predictor
   \f$\Delta\lambda=s/u_{t,k}\f$, \f$\Delta U=\Delta\lambda\,u_t\f$ satisfies
   \f$\Delta U_k=s\f$ up to rounding; the equilibrium is only approximated to first order.

   The path-constraint test's \f$|g|/s^2\f$ is not scale-invariant for a
   displacement \f$s\f$: it is \f$|\Delta U_k-s|/s^2\f$. A preset over this
   constraint chooses its constraint tolerance with that in mind.

   The exact zero pivot \f$u_{t,k}=0\f$ is a structural singularity of the
   bordered system (the dof is not controllable by \f$\lambda\f$ at this
   tangent) and is the one numerical case reported as failure.

   \ingroup gsStructuralAnalysis
*/
template<class T>
class gsDisplacementConstraint : public gsPathConstraint<T>
{
public:
    /// @param[in] dof index of the controlled tangent dof, non-negative (GISMO_ENSURE)
    explicit gsDisplacementConstraint(index_t dof);

    /// Set the controlled dof, non-negative (GISMO_ENSURE)
    void setDof(index_t dof);

    /// The controlled dof
    index_t dof() const { return m_dof; }

    /// Always true.
    bool predictorNeedsTangent(const gsStepContext<T> & ctx) const;

    /// \f$\delta\lambda=s/u_{t,k}\f$, \f$\delta U=\delta\lambda\,u_t\f$, copied
    /// to DeltaL and DeltaU. False iff dof() >= space.tangentSize(), dUt has
    /// the wrong size, or \f$u_{t,k}=0\f$. O(n).
    bool predict(gsStepContext<T> & ctx, const gsVector<T> & dUt,
                 const gsStateSpace<T> & space);

    /// \f$\delta\lambda=-(g+\delta\bar U_k)/u_{t,k}\f$, \f$\delta U=\delta\bar U+\delta\lambda\,u_t\f$.
    /// False on a size mismatch, dof() >= ctx.DeltaU.size(), or \f$u_{t,k}=0\f$. O(n).
    bool correct(gsStepContext<T> & ctx, const gsVector<T> & dUbar,
                 const gsVector<T> & dUt);

    /// Signed \f$\Delta U_k-s\f$; NaN if dof() is out of range. O(1).
    T residual(const gsStepContext<T> & ctx) const;

    /// \f$|\Delta U_k|\f$; NaN if dof() is out of range. O(1).
    T distance(const gsVector<T> & DeltaU, T DeltaL) const;

    /// False: @a s is a displacement increment.
    bool stepIsLoadIncrement() const;

    /// "Displacement"
    std::string name() const;

private:
    index_t m_dof;
};

/**
   @brief Riks' arc-length constraint with dof-averaged weighting.

   \f[ g = \varphi\,\Delta U^T\Delta U + (1-\varphi)\,\Delta\lambda^2 - s^2 = 0,
   \qquad \varphi = 1/n. \f]
   The weight \f$\varphi=1/n\f$ makes the displacement term the dof-average
   squared increment, so \f$\varphi\|\Delta U\|^2\f$ and
   \f$(1-\varphi)\Delta\lambda^2\f$ are commensurate for any mesh size
   (scaling()). For \f$n=1\f$ the value 0.99999999999 keeps the load term
   alive. The class is stateless: \f$\varphi\f$ is a function of the
   tangent size. The path-constraint test's \f$|g|/s^2\f$ is the relative
   violation of the constraint.

   <b>First-step predictor</b> (no secant). The tangent direction
   \f$(u_t,1)\f$, \f$u_t=K^{-1}\partial R/\partial\lambda\f$, is scaled by
   \f$\Delta\lambda\f$ so that
   \f$\varphi\Delta\lambda^2u_t\!\cdot\!u_t+(1-\varphi)\Delta\lambda^2=s^2\f$,
   i.e. \f$\Delta\lambda=s/\sqrt{\varphi\,u_t\!\cdot\!u_t+(1-\varphi)}\f$; the
   predicted point lies on the constraint.

   <b>Secant predictor.</b> The previous increment is rescaled by
   \f$s/s_{prev}\f$. It lies on the constraint iff the previous increment
   satisfied it with \f$s_{prev}\f$: ctx.prev.stepSize is the length the
   previous increment is normalised by (the driver stores the accepted step's
   \f$s\f$; a caller that changes the length between steps or supplies its own
   previous point sets it accordingly). Started off the constraint by a factor
   \f$\kappa\f$, Newton on this quadratic constraint roughly halves the excess
   per iteration, which costs about \f$\log_2\kappa\f$ iterations independent
   of \f$s\f$. The first-step predictor is used when ctx.prev is invalid, has
   a size different from ctx.F, or the previous increment is numerically zero.

   <b>Corrector.</b> Linearising \f$g\f$ about the current increment,
   \f$g+2\varphi\,\Delta U\cdot(\delta\bar U+\delta\lambda\,u_t)+2(1-\varphi)\Delta\lambda\,\delta\lambda=0\f$,
   gives
   \f[ \delta\lambda=-\frac{g+2\varphi\,\Delta U\cdot\delta\bar U}
       {2(1-\varphi)\Delta\lambda+2\varphi\,\Delta U\cdot u_t},\qquad
       \delta U=\delta\bar U+\delta\lambda\,u_t. \f]
   A zero or non-finite denominator is not a failure; the NaN or infinite
   update reaches the convergence tests.

   \ingroup gsStructuralAnalysis
*/
template<class T>
class gsRiksConstraint : public gsPathConstraint<T>
{
public:
    /// Weight \f$\varphi\f$ of the displacement term for tangent size @a n:
    /// 1/n for n >= 2, 0.99999999999 for n = 1, 0 for n <= 0. O(1).
    static T scaling(index_t n);

    /// True iff the first-step predictor is used, i.e. the tangent is needed.
    /// Reads ctx.prev and ctx.F.size(). O(n).
    bool predictorNeedsTangent(const gsStepContext<T> & ctx) const;

    /// First-step or secant predictor, scaled to length ctx.stepSize. False
    /// iff the first-step branch is taken and dUt.size() differs from
    /// space.tangentSize(). O(n).
    bool predict(gsStepContext<T> & ctx, const gsVector<T> & dUt,
                 const gsStateSpace<T> & space);

    /// Newton update of \f$\delta\lambda\f$ on the quadratic constraint. False
    /// iff dUbar, dUt and ctx.DeltaU differ in size. O(n).
    bool correct(gsStepContext<T> & ctx, const gsVector<T> & dUbar,
                 const gsVector<T> & dUt);

    /// Signed \f$\varphi\|\Delta U\|^2+(1-\varphi)\Delta\lambda^2-s^2\f$. O(n).
    T residual(const gsStepContext<T> & ctx) const;

    /// \f$\sqrt{\varphi\|\Delta U\|^2+(1-\varphi)\Delta\lambda^2}\f$ with
    /// \f$\varphi\f$ from DeltaU.size(). O(n).
    T distance(const gsVector<T> & DeltaU, T DeltaL) const;

    /// False: ctx.stepSize is an arc length.
    bool stepIsLoadIncrement() const;

    /// "Riks"
    std::string name() const;

private:
    /// True iff there is no usable secant in ctx.prev
    bool firstStep(const gsStepContext<T> & ctx) const;
};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsPathConstraint.hpp)
#endif
