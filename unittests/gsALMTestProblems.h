/** @file gsALMTestProblems.h

    @brief Shared closed-form 2-DOF test problems for the gsALMSolvers unit tests.

    Three fixtures, all pure algebra (no shell / KLShell / assembler), hence
    exactly reproducible run-to-run and safe with tight tolerances:

    == Fixture P (pitchfork) ==
    Potential  Pi(u,lambda) = 1/2 u1^2 + 1/2 u2^2 - 1/2 u1 u2^2 + 1/4 u2^4 - lambda u1.
    Residual   R = grad_u Pi = F_int(u) - lambda*Force, Force = (1,0)^T:
        R1 = u1 - 0.5*u2^2 - lambda
        R2 = u2 - u1*u2 + u2^3
    Symmetric tangent K = dR/du:
        K = [ 1      -u2              ]
            [ -u2    1 - u1 + 3*u2^2  ]
    Closed form:
      * fundamental branch  u2 = 0, u1 = lambda, K = diag(1, 1-lambda):
        stable for lambda<1, unstable for lambda>1;
      * pitchfork at lambda* = 1, u* = (1,0); critical mode V ~ (0,1)^T with
        V.Force = 0 EXACTLY  => BRANCH point;
      * bifurcated branch (supercritical, lambda>1):
        u1 = 2*lambda-1, u2 = +/- sqrt(2(lambda-1)), stable (det K = 2(lambda-1)>0).

    == Fixture F (fold) ==
    Same residual convention, a single limit point instead of a pitchfork:
        F_int(u) = ( 2*u1 - u1^2, u2 ),  Force = (1,0)^T
        R1 = 2*u1 - u1^2 - lambda
        R2 = u2
        K  = diag( 2 - 2*u1, 1 )
    Equilibria: u2 = 0, lambda = 2*u1 - u1^2 -- one parabola with a fold at
    u* = (1,0), lambda* = 1. Stable for u1<1, unstable for u1>1. The critical
    mode at the fold is V = (1,0)^T, so |V.Force|/|Force| = 1 >> tol => LIMIT
    point (contrast fixture P). The path crosses lambda = 0 a second time at
    u1 = 2, which is the M6 (phi = |U|/(|L| |f|)) lever.

    == Fixture F-n (padded fold) ==
    Fixture F with n-2 extra DECOUPLED unit DOFs:
        F_int(u) = ( 2*u0 - u0^2,  u1, u2, ..., u_{n-1} ),  Force = s*e0
        R0 = s*(2*u0 - u0^2 - lambda),   Ri = s*ui  (i >= 1)
        K  = s*diag( 2 - 2*u0, 1, 1, ..., 1 )
    The equilibrium manifold (the u0-parabola with ui = 0), the fold at
    (1,0,...,0) with lambda* = 1, the critical mode V = e0 and the tangent
    a^Q = K^-1 f = e0/(2-2*u0) are ALL IDENTICAL to fixture F for every n. The
    only thing n changes is m_numDof, hence gsALMRiks's phi = 1/n and its
    effective load scaling psi^2 = (1-phi)/phi = n-1 (the `m_convexWeight = 1./m_numDof;`
    assignment in gsALMRiks::predictor()).
    n is therefore a pure PSI^2 DIAL on an otherwise invariant problem -- the
    same kind of lever fixture S is for the load scale s. At n = 2 the fixture
    reduces to fixture F itself (psi^2 = 1), at n = 128 psi^2 = 127 and the
    constraint carries ~85% of its arc length in the lambda term even AT the
    fold, i.e. Riks is there very nearly load control.

    == Fixture S (scaled) ==
    NOT a separate problem: both fixtures above take a scale s. The residual is
    scaled to  R_s = s*R  and the force to  f_s = s*f, so that
        K_s = s*K,  equilibria (u,lambda) UNCHANGED,  null vector V UNCHANGED,
        delta_u_t = K_s^{-1} f_s = K^{-1} f  UNCHANGED,
        |V.f_s|/|f_s| = |V.f|/|f|            UNCHANGED.
    Every quantity a scale-invariant arc-length method may depend on is therefore
    invariant under s, which makes s the lever that exposes every scale-dependence
    bug (C1, M3, ...). s = 1 reproduces the historic fixtures BIT-EXACTLY
    (multiplication by 1.0 is exact in IEEE754 and the expression trees are
    unchanged).

    == Fixture A (asymmetric fold) ==
    The only fixture whose LEFT and RIGHT null vectors
    give OPPOSITE limit-vs-branch verdicts. A non-gradient system -- K is NOT
    symmetric, on purpose:
        R1 = s*( u1 - lambda )
        R2 = s*( -u1 + 2*u2 - u2*u2 )
        K  = s*[ [ 1,  0 ],
                 [ -1, 2 - 2*u2 ] ]                 (K12 = 0, K21 = -1: NON-symmetric)
        F_int(u) = R(u,0) = s*( u1, -u1 + 2*u2 - u2*u2 ),  Force = s*(1,0)^T.
    Equilibrium manifold: u1 = lambda, lambda = 2*u2 - u2^2 -- the same parabola as
    fixture F, now in u2. Fold at u2* = 1, u1* = 1, lambda* = 1; det K = s^2*(2-2*u2)
    changes sign there exactly as in fixture F. At the fold K* = [[1,0],[-1,0]]:
      * RIGHT null vector phi = (0,1)  =>  |phi.f|/|f| = 0 EXACTLY
        => the pre-task-01 (right-vector) test misclassifies this fold as a BRANCH
        point at every tolerance;
      * LEFT null vector psi = (1,1)/sqrt(2)  =>  |psi.f|/|f| = 1/sqrt(2) = 0.70710678...
        => the task-01 (left-vector) test classifies it correctly as a LIMIT point
        (Fredholm: the branch has dlambda/ds = 0 iff psi.Force = -psi.R_lambda != 0).
    A non-symmetric tangent needs a non-symmetric solver ("LU"); BifurcationMethod
    must be "Nothing" (-1) since neither Determinant (needs SimplicialLDLT) nor
    Eigenvalue (symmetric eigensolver) applies -- see gsALMTestProblems_test usage
    notes at the call site.

    == Fixture PF (subcritical pitchfork with a fold on the emanating branch) ==
    A gradient system (symmetric K, so the ordinary
    SimplicialLDLT/Determinant regime applies), obtained from fixture P by making
    the pitchfork SUBcritical and adding a u2^6 term that turns the emanating
    branch around:
        Pi = 1/2 u1^2 - lambda u1 + 1/2 u2^2 - 1/2 u1 u2^2 - 1/8 u2^4 + 1/6 u2^6
        R1 = s*( u1 - 0.5*u2^2 - lambda )
        R2 = s*( u2 - u1*u2 - 0.5*u2^3 + u2^5 )
        K  = s*[ [ 1,     -u2                        ],
                 [ -u2,   1 - u1 - 1.5*u2^2 + 5*u2^4  ] ],   Force = s*(1,0)^T.
    Fundamental branch u2=0, u1=lambda, K=diag(1,1-lambda): IDENTICAL to fixture P
    -- pitchfork at lambda*=1, u*=(1,0), critical mode (0,1), V.Force=0 => BRANCH.
    Emanating branch, with t = u2^2 > 0:
        lambda(t) = 1 - t + t^2,   u1(t) = lambda(t) + t/2 = 1 - t/2 + t^2,   u2 = +/- sqrt(t).
    SUBcritical (dlambda/dt = -1 < 0 at t=0). det K = 2t(2t-1) on that branch =>
    negatives = 1 for 0 < t < 1/2 (born unstable), negatives = 0 for t > 1/2 (the
    branch RE-STABILISES at a fold). Fold on the emanating branch at t = 1/2:
        u2* = +/- 1/sqrt(2) = +/-0.70710678118654752,  u1* = 1.0,  lambda* = 0.75.
    There K = [[1,-u2],[-u2,u2^2]], null vector V ~ (u2,1) = (0.7071,1),
    |V.f|/|f| = 0.57735... => classifies as a LIMIT point on the child curve, so
    the task-03b deferred-flip detection takes the markBifurcation path with no
    extended solve involved.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
**/

#pragma once

#include <functional>

#include <gsCore/gsLinearAlgebra.h>
#include <gsStructuralAnalysis/src/gsStructuralAnalysisTools/gsStructuralAnalysisTypes.h>

// ---------------------------------------------------------------------------
// Marker for a test that pins a defect which is KNOWN-UNFIXED today.
//
// The assertions below the marker are still COMPILED (so they cannot bit-rot)
// but are not EXECUTED, which keeps the committed suite green. Every use is
// accompanied by a
//      // DISABLED(task NN, <finding>): <one-line reason>
// comment above the TEST, so
//      grep -rn "DISABLED(task" optional/gsStructuralAnalysis/unittests
// lists everything that is waiting on a fix, with the task that will enable it.
//
// Build with -DGSALM_RUN_DISABLED_TESTS to run them all and see which fixes have
// landed.
// ---------------------------------------------------------------------------
#ifdef GSALM_RUN_DISABLED_TESTS
#  define GSALM_DISABLED(task, finding, reason) ((void)0)
#else
#  define GSALM_DISABLED(task, finding, reason)                                 \
    do {                                                                        \
        gsInfo << "  [DISABLED(" << task << ", " << finding << ")] "             \
               << reason                                                        \
               << " -- assertions below are compiled, not executed.\n";         \
        return;                                                                 \
    } while (false)
#endif

namespace gsALMTest
{

using namespace gismo;

typedef gsStructuralAnalysisOps<real_t>::Jacobian_t   ALMJacobian_t;
typedef gsStructuralAnalysisOps<real_t>::ALResidual_t ALMResidual_t;

// ---------------------------------------------------------------------------
// Dense 2x2 -> sparse, in the exact order the historic fixtures used.
// ---------------------------------------------------------------------------
inline void toSparse2(const gsMatrix<real_t> & K, gsSparseMatrix<real_t> & m)
{
    m.resize(2,2);
    m.setZero();
    m.coeffRef(0,0) = K(0,0);
    m.coeffRef(0,1) = K(0,1);
    m.coeffRef(1,0) = K(1,0);
    m.coeffRef(1,1) = K(1,1);
    m.makeCompressed();
}

// ---------------------------------------------------------------------------
// Fixture P: pitchfork. See the file header.
// ---------------------------------------------------------------------------
inline gsVector<real_t> pitchforkResidual(const gsVector<real_t> & u, real_t lambda,
                                          real_t s = 1.0)
{
    gsVector<real_t> R(2);
    const real_t u1 = u[0], u2 = u[1];
    R[0] = s*(u1 - 0.5*u2*u2 - lambda);
    R[1] = s*(u2 - u1*u2 + u2*u2*u2);
    return R;
}

inline gsMatrix<real_t> pitchforkTangent(const gsVector<real_t> & u, real_t s = 1.0)
{
    gsMatrix<real_t> K(2,2);
    const real_t u1 = u[0], u2 = u[1];
    K(0,0) = s*(1.0);          K(0,1) = s*(-u2);
    K(1,0) = s*(-u2);          K(1,1) = s*(1.0 - u1 + 3.0*u2*u2);
    return K;
}

// ---------------------------------------------------------------------------
// Fixture F: fold. See the file header.
// ---------------------------------------------------------------------------
inline gsVector<real_t> foldResidual(const gsVector<real_t> & u, real_t lambda,
                                     real_t s = 1.0)
{
    gsVector<real_t> R(2);
    const real_t u1 = u[0], u2 = u[1];
    R[0] = s*(2.0*u1 - u1*u1 - lambda);
    R[1] = s*(u2);
    return R;
}

inline gsMatrix<real_t> foldTangent(const gsVector<real_t> & u, real_t s = 1.0)
{
    gsMatrix<real_t> K(2,2);
    const real_t u1 = u[0];
    K(0,0) = s*(2.0 - 2.0*u1);   K(0,1) = s*(0.0);
    K(1,0) = s*(0.0);            K(1,1) = s*(1.0);
    return K;
}

// ---------------------------------------------------------------------------
// Closed-form oracles (independent of any solver output -- never re-paste a
// computed value as "truth").
// ---------------------------------------------------------------------------

/// Fixture P, bifurcated branch: u1(lambda)
inline real_t pitchforkBranchU1  (real_t lambda) { return 2.0*lambda - 1.0; }
/// Fixture P, bifurcated branch: u2(lambda)^2
inline real_t pitchforkBranchU2sq(real_t lambda) { return 2.0*(lambda - 1.0); }
/// Fixture F, equilibrium manifold: lambda(u1) (with u2 = 0)
inline real_t foldLambda(real_t u1) { return 2.0*u1 - u1*u1; }

// ---------------------------------------------------------------------------
// Bundled operators. gsALMBase copies the callbacks and the Force vector, so an
// AlmProblem may die before the solver it configured.
// ---------------------------------------------------------------------------
struct AlmProblem
{
    /// Load scale s (see the file header). All members below are consistent with it.
    real_t           scale;
    /// f = -dR/dLambda = s * (1,0)^T
    gsVector<real_t> Force;
    ALMJacobian_t    Jacobian;
    ALMResidual_t    ALResidual;
    /// Closed-form internal force F_int(u) = R(u,0); used as an oracle only.
    std::function<gsVector<real_t>(const gsVector<real_t> &)> Fint;
};

inline AlmProblem pitchforkProblem(real_t s = 1.0)
{
    AlmProblem p;
    p.scale = s;
    p.Force.resize(2);
    p.Force << s*(1.0), s*(0.0);
    p.Jacobian =
        [s](gsVector<real_t> const & u, gsSparseMatrix<real_t> & m) -> bool
        { toSparse2(pitchforkTangent(u,s), m); return true; };
    p.ALResidual =
        [s](gsVector<real_t> const & u, const real_t lambda, gsVector<real_t> & result) -> bool
        { result = pitchforkResidual(u,lambda,s); return true; };
    p.Fint =
        [s](const gsVector<real_t> & u) { return pitchforkResidual(u,0.0,s); };
    return p;
}

inline AlmProblem foldProblem(real_t s = 1.0)
{
    AlmProblem p;
    p.scale = s;
    p.Force.resize(2);
    p.Force << s*(1.0), s*(0.0);
    p.Jacobian =
        [s](gsVector<real_t> const & u, gsSparseMatrix<real_t> & m) -> bool
        { toSparse2(foldTangent(u,s), m); return true; };
    p.ALResidual =
        [s](gsVector<real_t> const & u, const real_t lambda, gsVector<real_t> & result) -> bool
        { result = foldResidual(u,lambda,s); return true; };
    p.Fint =
        [s](const gsVector<real_t> & u) { return foldResidual(u,0.0,s); };
    return p;
}

// ---------------------------------------------------------------------------
// Fixture F-n: the fold, padded with n-2 decoupled unit DOFs. See the file
// header -- n is a pure psi^2 dial for gsALMRiks (psi^2 = n-1).
// ---------------------------------------------------------------------------

/// Diagonal n x n -> sparse.
inline void toSparseDiag(const gsVector<real_t> & d, gsSparseMatrix<real_t> & m)
{
    const index_t n = d.size();
    m.resize(n,n);
    m.setZero();
    m.reserve(gsVector<index_t>::Constant(n,1));
    for (index_t i = 0; i != n; ++i)
        m.coeffRef(i,i) = d[i];
    m.makeCompressed();
}

/// Fixture F padded to \a n DOFs at load scale \a s. n = 2 reproduces
/// foldProblem(s) exactly (same expression trees, same coefficient order).
inline AlmProblem foldProblemPadded(index_t n, real_t s = 1.0)
{
    GISMO_ASSERT(n >= 2, "foldProblemPadded: need at least the 2 DOFs of fixture F");
    AlmProblem p;
    p.scale = s;
    p.Force = gsVector<real_t>::Zero(n);
    p.Force[0] = s*(1.0);
    p.Jacobian =
        [n,s](gsVector<real_t> const & u, gsSparseMatrix<real_t> & m) -> bool
        {
            gsVector<real_t> d = gsVector<real_t>::Constant(n, s*(1.0));
            d[0] = s*(2.0 - 2.0*u[0]);
            toSparseDiag(d,m);
            return true;
        };
    p.ALResidual =
        [s](gsVector<real_t> const & u, const real_t lambda, gsVector<real_t> & result) -> bool
        {
            result = s*u;
            result[0] = s*(2.0*u[0] - u[0]*u[0] - lambda);
            return true;
        };
    p.Fint =
        [s](const gsVector<real_t> & u)
        {
            gsVector<real_t> F = s*u;
            F[0] = s*(2.0*u[0] - u[0]*u[0]);
            return F;
        };
    return p;
}

/// A LINEAR 2-DOF control problem: F_int(u) = diag(2,1)*u, Force = s*(1,0)^T.
/// No fold, no bifurcation; every quantity that is exact only for a linear
/// internal force (e.g. gsALMCrisfield's complex-root fallback, see M4) must be
/// correct here both before and after a fix.
inline AlmProblem linearProblem(real_t s = 1.0)
{
    AlmProblem p;
    p.scale = s;
    p.Force.resize(2);
    p.Force << s*(1.0), s*(0.0);
    p.Jacobian =
        [s](gsVector<real_t> const & /*u*/, gsSparseMatrix<real_t> & m) -> bool
        {
            gsMatrix<real_t> K(2,2);
            K(0,0) = s*(2.0); K(0,1) = 0.0;
            K(1,0) = 0.0;     K(1,1) = s*(1.0);
            toSparse2(K,m);
            return true;
        };
    p.ALResidual =
        [s](gsVector<real_t> const & u, const real_t lambda, gsVector<real_t> & result) -> bool
        {
            result.resize(2);
            result[0] = s*(2.0*u[0] - lambda);
            result[1] = s*(u[1]);
            return true;
        };
    p.Fint =
        [s](const gsVector<real_t> & u)
        {
            gsVector<real_t> F(2);
            F[0] = s*(2.0*u[0]);
            F[1] = s*(u[1]);
            return F;
        };
    return p;
}

// ---------------------------------------------------------------------------
// Fixture A: asymmetric fold. See the file header. NOT a gradient system: K is
// NOT symmetric (K12=0, K21=-1) -- this is the whole point.
// ---------------------------------------------------------------------------
inline gsVector<real_t> asymFoldResidual(const gsVector<real_t> & u, real_t lambda,
                                         real_t s = 1.0)
{
    gsVector<real_t> R(2);
    const real_t u1 = u[0], u2 = u[1];
    R[0] = s*(u1 - lambda);
    R[1] = s*(-u1 + 2.0*u2 - u2*u2);
    return R;
}

inline gsMatrix<real_t> asymFoldTangent(const gsVector<real_t> & u, real_t s = 1.0)
{
    gsMatrix<real_t> K(2,2);
    const real_t u2 = u[1];
    K(0,0) = s*(1.0);   K(0,1) = s*(0.0);
    K(1,0) = s*(-1.0);  K(1,1) = s*(2.0 - 2.0*u2);
    return K;
}

/// Fixture A, equilibrium manifold: lambda(u2) (with u1 = lambda)
inline real_t asymFoldLambda(real_t u2) { return 2.0*u2 - u2*u2; }

inline AlmProblem asymFoldProblem(real_t s = 1.0)
{
    AlmProblem p;
    p.scale = s;
    p.Force.resize(2);
    p.Force << s*(1.0), s*(0.0);
    p.Jacobian =
        [s](gsVector<real_t> const & u, gsSparseMatrix<real_t> & m) -> bool
        { toSparse2(asymFoldTangent(u,s), m); return true; };
    p.ALResidual =
        [s](gsVector<real_t> const & u, const real_t lambda, gsVector<real_t> & result) -> bool
        { result = asymFoldResidual(u,lambda,s); return true; };
    p.Fint =
        [s](const gsVector<real_t> & u) { return asymFoldResidual(u,0.0,s); };
    return p;
}

// ---------------------------------------------------------------------------
// Fixture PF: subcritical pitchfork with a fold on the emanating branch. See the
// file header. A gradient system (symmetric K).
// ---------------------------------------------------------------------------
inline gsVector<real_t> subPitchforkResidual(const gsVector<real_t> & u, real_t lambda,
                                             real_t s = 1.0)
{
    gsVector<real_t> R(2);
    const real_t u1 = u[0], u2 = u[1];
    R[0] = s*(u1 - 0.5*u2*u2 - lambda);
    R[1] = s*(u2 - u1*u2 - 0.5*u2*u2*u2 + u2*u2*u2*u2*u2);
    return R;
}

inline gsMatrix<real_t> subPitchforkTangent(const gsVector<real_t> & u, real_t s = 1.0)
{
    gsMatrix<real_t> K(2,2);
    const real_t u1 = u[0], u2 = u[1];
    K(0,0) = s*(1.0);   K(0,1) = s*(-u2);
    K(1,0) = s*(-u2);   K(1,1) = s*(1.0 - u1 - 1.5*u2*u2 + 5.0*u2*u2*u2*u2);
    return K;
}

/// Fixture PF, emanating branch parametrised by t = u2^2 > 0: lambda(t)
inline real_t subPitchforkBranchLambda(real_t t) { return 1.0 - t + t*t; }
/// Fixture PF, emanating branch: u1(t)
inline real_t subPitchforkBranchU1(real_t t) { return subPitchforkBranchLambda(t) + 0.5*t; }
/// Fixture PF, fold on the emanating branch: t* = 1/2
static const real_t subPitchforkFoldT      = 0.5;
static const real_t subPitchforkFoldU2     = 0.70710678118654752;  // 1/sqrt(2)
static const real_t subPitchforkFoldU1     = 1.0;
static const real_t subPitchforkFoldLambda = 0.75;

inline AlmProblem subPitchforkProblem(real_t s = 1.0)
{
    AlmProblem p;
    p.scale = s;
    p.Force.resize(2);
    p.Force << s*(1.0), s*(0.0);
    p.Jacobian =
        [s](gsVector<real_t> const & u, gsSparseMatrix<real_t> & m) -> bool
        { toSparse2(subPitchforkTangent(u,s), m); return true; };
    p.ALResidual =
        [s](gsVector<real_t> const & u, const real_t lambda, gsVector<real_t> & result) -> bool
        { result = subPitchforkResidual(u,lambda,s); return true; };
    p.Fint =
        [s](const gsVector<real_t> & u) { return subPitchforkResidual(u,0.0,s); };
    return p;
}

// ---------------------------------------------------------------------------
// Call-counting / fault-injecting operator wrapper.
//
// Forwards to a wrapped AlmProblem while counting invocations, and can be told
// to return false on the n-th call. The `false` return is the documented
// assembly-failure convention of gsStructuralAnalysisOps; the solvers must map
// it to gsStatus::AssemblyError.
//
// The callbacks returned by residual()/jacobian() capture `this`, so the probe
// must OUTLIVE every solver it was handed to.
// ---------------------------------------------------------------------------
class AlmOperatorProbe
{
public:
    explicit AlmOperatorProbe(const AlmProblem & problem)
    :
    m_problem(problem),
    m_nResidual(0), m_nJacobian(0),
    m_failResidualAt(-1), m_failJacobianAt(-1)
    {}

    /// Force the \a n-th residual evaluation (1-based) to report failure; -1 disables.
    void failResidualAt(index_t n) { m_failResidualAt = n; }
    /// Force the \a n-th Jacobian evaluation (1-based) to report failure; -1 disables.
    void failJacobianAt(index_t n) { m_failJacobianAt = n; }

    index_t residualCalls() const { return m_nResidual; }
    index_t jacobianCalls() const { return m_nJacobian; }

    void resetCounters() { m_nResidual = 0; m_nJacobian = 0; }

    const gsVector<real_t> & force() const { return m_problem.Force; }

    ALMResidual_t residual()
    {
        return [this](gsVector<real_t> const & u, const real_t lambda,
                      gsVector<real_t> & result) -> bool
        {
            ++m_nResidual;
            if (m_nResidual == m_failResidualAt) return false;
            return m_problem.ALResidual(u,lambda,result);
        };
    }

    ALMJacobian_t jacobian()
    {
        return [this](gsVector<real_t> const & u, gsSparseMatrix<real_t> & m) -> bool
        {
            ++m_nJacobian;
            if (m_nJacobian == m_failJacobianAt) return false;
            return m_problem.Jacobian(u,m);
        };
    }

private:
    AlmProblem m_problem;
    index_t    m_nResidual;
    index_t    m_nJacobian;
    index_t    m_failResidualAt;
    index_t    m_failJacobianAt;
};

} // namespace gsALMTest
