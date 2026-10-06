 /** @file gsConvergenceTest.h

    @brief Convergence tests of the solver driver: leaf tests on the step
    context, All/Any composition, and the arc-length gate.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#pragma once
#include <limits>
#include <vector>
#include <string>
#include <cmath>

#include <gsCore/gsLinearAlgebra.h>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsSolverDriverTypes.h>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsStateSpace.h>

namespace gismo
{

/// Outcome of a convergence test at one iterate.
enum class gsConvergenceStatus
{
    Unconverged, ///< keep iterating
    Converged,   ///< accept the iterate
    Exhausted,   ///< iteration budget used up: the step ends NotConverged (not an error of the test)
    Failed       ///< non-finite or diverging iterate: the step ends NotConverged with failedOp = gsDriverOp::Test
};

/// Reference against which a relative residual test measures
enum class gsResidualReference
{
    R0,      ///< norm of the residual at the predicted point (iteration 0)
    LoadNorm ///< norm of the external load at the predicted load level, \f$\|(\lambda+\Delta\lambda)F\|\f$
};

/**
   @brief Base class of the convergence tests of the solver driver.

   A test reads only the step context and the const state space. Leaf tests
   are of two kinds: <em>tolerance tests</em> return only Converged or
   Unconverged; <em>guards</em> (gsMaxIterTest, gsDivergedTest) return only
   Exhausted/Failed or Unconverged and are never Converged. A tolerance test
   whose measure is NaN or infinite returns Unconverged. Only
   gsDivergedTest turns non-finite values into Failed.

   Tests are combined with gsAllTest and gsAnyTest, whose result follows a
   fixed precedence (see there).

   @b Call @b protocol. For each step attempt, retries included, the driver:
   1. predicts;
   2. evaluates the residual (and the constraint residual) at the predicted point;
   3. calls initStep() once;
   4. calls test() at iteration 0;
   5. after every corrector update, evaluates the residual and calls test().

   A test whose initStep() was not called in the current attempt holds stale
   reference values, so the driver must not skip step 3.

   initStep() and test() are non-const because the tests store their
   reference values and last measure. The driver holds the tests through
   Ptr, so a const driver method may still call them.

   @tparam T coefficient type

   \ingroup gsStructuralAnalysis
*/
template<class T>
class gsConvergenceTest
{
public:
    typedef memory::shared_ptr<gsConvergenceTest> Ptr;
    virtual ~gsConvergenceTest() {}

    /// Captures the per-step reference values (initial residual, load norm,
    /// predictor increment). Called once per step attempt at
    /// ctx.iteration == 0, before the first test(). Composites forward it to
    /// every child. Cost: O(n), \f$n\f$ = tangent size.
    virtual void initStep(const gsStepContext<T> & ctx, const gsStateSpace<T> & space)
    { GISMO_UNUSED(ctx); GISMO_UNUSED(space); }

    /// Evaluate the test at iteration 0 (after initStep()) and after every
    /// corrector update. Cost: see the leaf class.
    virtual gsConvergenceStatus test(const gsStepContext<T> & ctx, const gsStateSpace<T> & space) = 0;

    /// True for gsMaxIterTest and gsDivergedTest, and for a gsAnyTest
    /// containing a guard.
    virtual bool isGuard() const { return false; }

    /// Iteration limit \f$N\f$ of the tightest gsMaxIterTest in this tree;
    /// -1 if the tree has none.
    virtual index_t maxIterations() const { return -1; }

    /// Last measured quantity (ratio, norm, count); quiet NaN before the
    /// first test() and for composites.
    virtual T value() const { return std::numeric_limits<T>::quiet_NaN(); }

    /// Short identifier for verbose output, e.g. "RelResidual(LoadNorm)".
    virtual std::string name() const = 0;

    /**
       @brief The acceptance rule of the arc-length family (gsALMBase and
       its presets), as a test tree.

       Builds
       <pre>
       Any( All( MinIter(1),
                 RelResidual(tolF, LoadNorm),
                 RelUpdate(tolU),
                 PathConstraint(tolC) ),
            MaxIter(maxIter) )
       </pre>
       It reproduces the quantities of gsALMBase:
       - residueF = |R| / |(lambda + DeltaL) F| (basis fixed at iteration 0),
         tested strictly against tolF (option <tt>TolF</tt>);
       - residueU = |deltaU| / |DeltaU of the predictor|, tested strictly
         against tolU (option <tt>TolU</tt>);
       - the constraint measure |g| / s^2 against tolC (option
         <tt>ConstraintTol</tt>); g is 0 for load control, whose measure is
         then 0 for every step size;
       - no acceptance at iteration 0 (MinIter);
       - <tt>MaxIter = N</tt> gives exactly N corrector updates: an iterate
         that converges at update N is accepted, otherwise the result is
         Exhausted;
       - a NaN residual is never Converged and the loop runs on to N,
         because the tree contains no gsDivergedTest.

       Complexity of one test(): O(n).
       @param[in] tolF tolerance of the relative residual
       @param[in] tolU tolerance of the relative update
       @param[in] tolC tolerance of the path constraint
       @param[in] maxIter number of corrector updates, non-negative
    */
    static Ptr makeArcLengthGate(T tolF, T tolU, T tolC, index_t maxIter);
};

/**
   @brief Relative residual test: \f$\|R\|/\text{basis} < \text{tol}\f$.

   The basis is either the norm of the residual at the predicted point
   (gsResidualReference::R0) or the norm of the load at the predicted load
   level (gsResidualReference::LoadNorm). A zero basis makes the ratio NaN or
   infinite, so the test never converges; combine it with gsAbsResidualTest
   under gsAnyTest if that case matters.

   Cost of test(): O(n).
*/
template<class T>
class gsRelResidualTest : public gsConvergenceTest<T>
{
public:
    /// @param[in] tol tolerance
    /// @param[in] ref reference of the relative measure
    gsRelResidualTest(T tol, gsResidualReference ref)
    : m_tol(tol), m_ref(ref),
      m_basis(std::numeric_limits<T>::quiet_NaN()),
      m_value(std::numeric_limits<T>::quiet_NaN())
    {}

    void initStep(const gsStepContext<T> & ctx, const gsStateSpace<T> &)
    {
        if (m_ref == gsResidualReference::LoadNorm)
            m_basis = ((ctx.lambda + ctx.DeltaL) * ctx.F).norm();
        else
            m_basis = ctx.R.norm();
    }

    gsConvergenceStatus test(const gsStepContext<T> & ctx, const gsStateSpace<T> &)
    {
        m_value = ctx.R.norm() / m_basis;
        return (m_value < m_tol) ? gsConvergenceStatus::Converged : gsConvergenceStatus::Unconverged;
    }

    T value() const { return m_value; }

    std::string name() const
    { return m_ref == gsResidualReference::LoadNorm ? "RelResidual(LoadNorm)" : "RelResidual(R0)"; }

private:
    T m_tol;
    gsResidualReference m_ref;
    T m_basis, m_value;
};

/**
   @brief Absolute residual test: \f$\|R\| < \text{tol}\f$.

   Cost of test(): O(n).
*/
template<class T>
class gsAbsResidualTest : public gsConvergenceTest<T>
{
public:
    /// @param[in] tol tolerance
    explicit gsAbsResidualTest(T tol)
    : m_tol(tol), m_value(std::numeric_limits<T>::quiet_NaN())
    {}

    gsConvergenceStatus test(const gsStepContext<T> & ctx, const gsStateSpace<T> &)
    {
        m_value = ctx.R.norm();
        return (m_value < m_tol) ? gsConvergenceStatus::Converged : gsConvergenceStatus::Unconverged;
    }

    T value() const { return m_value; }
    std::string name() const { return "AbsResidual"; }

private:
    T m_tol, m_value;
};

/**
   @brief Relative update test: \f$\|\delta U\|/\|\Delta U_{pred}\| < \text{tol}\f$.

   The basis is the predictor increment \f$\Delta U\f$ at iteration 0;
   \f$\delta U\f$ is the applied (relaxed) update of the latest corrector
   iteration. A zero predictor increment makes the ratio NaN or infinite, so
   the test never converges.

   Cost of test(): O(n).
*/
template<class T>
class gsRelUpdateTest : public gsConvergenceTest<T>
{
public:
    /// @param[in] tol tolerance
    explicit gsRelUpdateTest(T tol)
    : m_tol(tol),
      m_basis(std::numeric_limits<T>::quiet_NaN()),
      m_value(std::numeric_limits<T>::quiet_NaN())
    {}

    void initStep(const gsStepContext<T> & ctx, const gsStateSpace<T> &)
    { m_basis = ctx.DeltaU.norm(); }

    gsConvergenceStatus test(const gsStepContext<T> & ctx, const gsStateSpace<T> &)
    {
        m_value = ctx.deltaU.norm() / m_basis;
        return (m_value < m_tol) ? gsConvergenceStatus::Converged : gsConvergenceStatus::Unconverged;
    }

    T value() const { return m_value; }
    std::string name() const { return "RelUpdate"; }

private:
    T m_tol, m_basis, m_value;
};

/**
   @brief Path-constraint test: \f$|g|/s^2 < \text{tol}\f$, with
   \f$g\f$ = gsStepContext::constraintResidual and \f$s\f$ =
   gsStepContext::stepSize.

   For \f$g = 0\f$ exactly (fixed lambda, load and displacement control) the
   measure is 0 for every \f$s\f$, including \f$s = 0\f$, where the quotient
   would be 0/0. For \f$g \ne 0\f$ and \f$s = 0\f$ the measure is infinite
   and the test does not converge.

   Cost of test(): O(1).
*/
template<class T>
class gsPathConstraintTest : public gsConvergenceTest<T>
{
public:
    /// @param[in] tol tolerance
    explicit gsPathConstraintTest(T tol)
    : m_tol(tol), m_value(std::numeric_limits<T>::quiet_NaN())
    {}

    gsConvergenceStatus test(const gsStepContext<T> & ctx, const gsStateSpace<T> &)
    {
        const T g = ctx.constraintResidual;
        const T s = ctx.stepSize;
        m_value = (g == 0) ? T(0) : math::abs(g) / (s*s);
        return (m_value < m_tol) ? gsConvergenceStatus::Converged : gsConvergenceStatus::Unconverged;
    }

    T value() const { return m_value; }
    std::string name() const { return "PathConstraint"; }

private:
    T m_tol, m_value;
};

/**
   @brief Minimum iteration test: Converged iff at least @a m corrector
   updates have been applied. Used inside gsAllTest to forbid acceptance
   at iteration 0.

   Cost of test(): O(1).
*/
template<class T>
class gsMinIterTest : public gsConvergenceTest<T>
{
public:
    /// @param[in] m minimum number of corrector updates
    explicit gsMinIterTest(index_t m)
    : m_min(m), m_value(std::numeric_limits<T>::quiet_NaN())
    {}

    gsConvergenceStatus test(const gsStepContext<T> & ctx, const gsStateSpace<T> &)
    {
        m_value = T(ctx.iteration);
        return (ctx.iteration >= m_min) ? gsConvergenceStatus::Converged : gsConvergenceStatus::Unconverged;
    }

    T value() const { return m_value; }
    std::string name() const { return "MinIter"; }

private:
    index_t m_min;
    T m_value;
};

/**
   @brief Guard on the number of corrector updates: Exhausted iff
   gsStepContext::iteration >= N, so that exactly N updates are applied
   before the step fails. Never Converged.

   Cost of test(): O(1).
*/
template<class T>
class gsMaxIterTest : public gsConvergenceTest<T>
{
public:
    /// @param[in] n number of corrector updates
    explicit gsMaxIterTest(index_t n)
    : m_max(n), m_value(std::numeric_limits<T>::quiet_NaN())
    {}

    gsConvergenceStatus test(const gsStepContext<T> & ctx, const gsStateSpace<T> &)
    {
        m_value = T(ctx.iteration);
        return (ctx.iteration >= m_max) ? gsConvergenceStatus::Exhausted : gsConvergenceStatus::Unconverged;
    }

    bool isGuard() const { return true; }
    index_t maxIterations() const { return m_max; }
    T value() const { return m_value; }
    std::string name() const { return "MaxIter"; }

private:
    index_t m_max;
    T m_value;
};

/**
   @brief Guard against a non-finite or diverging iterate: Failed iff the
   residual norm is not finite, or (after the first update) the update norm
   is not finite, or <tt>factor > 0</tt> and \f$\|R\| > \text{factor}\,\|R_0\|\f$
   with \f$R_0\f$ the residual at iteration 0. Never Converged.

   Cost of test(): O(n).
*/
template<class T>
class gsDivergedTest : public gsConvergenceTest<T>
{
public:
    /// @param[in] factor growth factor of the residual norm that counts as
    /// divergence; non-positive switches the growth check off
    explicit gsDivergedTest(T factor)
    : m_factor(factor), m_R0(std::numeric_limits<T>::quiet_NaN()),
      m_value(std::numeric_limits<T>::quiet_NaN())
    {}

    void initStep(const gsStepContext<T> & ctx, const gsStateSpace<T> &)
    { m_R0 = ctx.R.norm(); }

    gsConvergenceStatus test(const gsStepContext<T> & ctx, const gsStateSpace<T> &)
    {
        m_value = ctx.R.norm();
        const bool bad = !std::isfinite(m_value)
            || (ctx.iteration >= 1 && !std::isfinite(ctx.deltaU.norm()))
            || (m_factor > 0 && m_value > m_factor * m_R0);
        return bad ? gsConvergenceStatus::Failed : gsConvergenceStatus::Unconverged;
    }

    bool isGuard() const { return true; }
    T value() const { return m_value; }
    std::string name() const { return "Diverged"; }

private:
    T m_factor, m_R0, m_value;
};

/**
   @brief Round-off floor of the residual: Converged iff
   \f$\|R\| \le c\,\varepsilon\,\|\,|K|\,\hat x\,\|_2\f$.

   The right-hand side is the forward-error bound of the internal-force
   product, with \f$\hat x\f$ = gsStateSpace::magnitude() (|x| on additive
   dofs, 0 elsewhere) and \f$\varepsilon\f$ the machine epsilon. Residuals
   below it cannot be reduced further by Newton updates, so the test is
   meant to sit under a gsAnyTest next to a relative residual test.

   The floor is off (Unconverged, value NaN) at iteration 0, for
   <tt>c <= 0</tt>, when the space has no magnitude, and when the magnitude
   does not have the size of K. Before the first update only the relative
   test may accept, because the floor can exceed the initial residual.

   gsStepContext::K is the most recently @e assembled tangent. When the
   residual is evaluated and tested before the next Jacobian assembly, that
   is the tangent of the previous iterate, and older when the tangent is
   reused over several iterations; the floor is therefore an estimate. For
   a rod the floor was measured at about 1e-9 relative.

   Cost of test(): O(nnz(K)); the buffers for \f$\hat x\f$ and \f$|K|\hat x\f$ are reused.
*/
template<class T>
class gsRoundoffFloorTest : public gsConvergenceTest<T>
{
public:
    /// @param[in] c safety factor of the floor; non-positive switches the test off
    explicit gsRoundoffFloorTest(T c)
    : m_c(c), m_eps(std::numeric_limits<T>::epsilon()),
      m_value(std::numeric_limits<T>::quiet_NaN())
    {}

    gsConvergenceStatus test(const gsStepContext<T> & ctx, const gsStateSpace<T> & space)
    {
        m_value = std::numeric_limits<T>::quiet_NaN();
        if (m_c <= 0 || ctx.iteration < 1)
            return gsConvergenceStatus::Unconverged;
        if (!space.magnitude(m_xh) || m_xh.size() == 0 || m_xh.size() != ctx.K.cols())
            return gsConvergenceStatus::Unconverged;

        m_Kxh.noalias() = ctx.K.cwiseAbs() * m_xh;
        m_value = m_c * m_eps * m_Kxh.norm();
        if (!std::isfinite(m_value))
            return gsConvergenceStatus::Unconverged;
        return (ctx.R.norm() <= m_value) ? gsConvergenceStatus::Converged : gsConvergenceStatus::Unconverged;
    }

    T value() const { return m_value; }
    std::string name() const { return "RoundoffFloor"; }

private:
    T m_c, m_eps, m_value;
    gsVector<T> m_xh, m_Kxh;
};

/**
   @brief Conjunction of tests.

   Result: Converged iff every child is Converged; otherwise Failed if any
   child is Failed, otherwise Exhausted if any child is Exhausted, otherwise
   Unconverged. An empty conjunction is Unconverged. Every child is
   evaluated on every call (no short-circuit), so that each child's value()
   stays current. Guards are rejected as children, since a guard is never
   Converged and would block acceptance for good.

   Cost of test(): the sum of the children's costs.
*/
template<class T>
class gsAllTest : public gsConvergenceTest<T>
{
public:
    typedef typename gsConvergenceTest<T>::Ptr Ptr;

    /// Add a child; throws (GISMO_ENSURE) if @a child is null or a guard.
    /// @return *this
    gsAllTest& add(const Ptr & child)
    {
        GISMO_ENSURE(child, "gsAllTest::add: null child");
        GISMO_ENSURE(!child->isGuard(), "gsAllTest::add: guard '" << child->name() << "' can never be Converged");
        m_children.push_back(child);
        return *this;
    }

    void initStep(const gsStepContext<T> & ctx, const gsStateSpace<T> & space)
    {
        for (size_t i = 0; i != m_children.size(); ++i)
            m_children[i]->initStep(ctx, space);
    }

    gsConvergenceStatus test(const gsStepContext<T> & ctx, const gsStateSpace<T> & space)
    {
        bool allConverged = !m_children.empty();
        bool failed = false, exhausted = false;
        for (size_t i = 0; i != m_children.size(); ++i)
        {
            const gsConvergenceStatus s = m_children[i]->test(ctx, space);
            allConverged = allConverged && (s == gsConvergenceStatus::Converged);
            failed       = failed    || (s == gsConvergenceStatus::Failed);
            exhausted    = exhausted || (s == gsConvergenceStatus::Exhausted);
        }
        if (allConverged) return gsConvergenceStatus::Converged;
        if (failed)       return gsConvergenceStatus::Failed;
        if (exhausted)    return gsConvergenceStatus::Exhausted;
        return gsConvergenceStatus::Unconverged;
    }

    std::string name() const { return "All"; }

    /// The children, in insertion order
    const std::vector<Ptr> & children() const { return m_children; }

private:
    std::vector<Ptr> m_children;
};

/**
   @brief Disjunction of tests.

   Result: Converged if any child is Converged (so an iterate that converges
   at the iteration limit is accepted); otherwise Failed if any child is
   Failed, otherwise Exhausted if any child is Exhausted, otherwise
   Unconverged. An empty disjunction is Unconverged. Every child is
   evaluated on every call (no short-circuit).

   Cost of test(): the sum of the children's costs.
*/
template<class T>
class gsAnyTest : public gsConvergenceTest<T>
{
public:
    typedef typename gsConvergenceTest<T>::Ptr Ptr;

    /// Add a child; throws (GISMO_ENSURE) if @a child is null.
    /// @return *this
    gsAnyTest& add(const Ptr & child)
    {
        GISMO_ENSURE(child, "gsAnyTest::add: null child");
        m_children.push_back(child);
        return *this;
    }

    void initStep(const gsStepContext<T> & ctx, const gsStateSpace<T> & space)
    {
        for (size_t i = 0; i != m_children.size(); ++i)
            m_children[i]->initStep(ctx, space);
    }

    gsConvergenceStatus test(const gsStepContext<T> & ctx, const gsStateSpace<T> & space)
    {
        bool converged = false, failed = false, exhausted = false;
        for (size_t i = 0; i != m_children.size(); ++i)
        {
            const gsConvergenceStatus s = m_children[i]->test(ctx, space);
            converged = converged || (s == gsConvergenceStatus::Converged);
            failed    = failed    || (s == gsConvergenceStatus::Failed);
            exhausted = exhausted || (s == gsConvergenceStatus::Exhausted);
        }
        if (converged) return gsConvergenceStatus::Converged;
        if (failed)    return gsConvergenceStatus::Failed;
        if (exhausted) return gsConvergenceStatus::Exhausted;
        return gsConvergenceStatus::Unconverged;
    }

    /// True iff any child is a guard
    bool isGuard() const
    {
        for (size_t i = 0; i != m_children.size(); ++i)
            if (m_children[i]->isGuard()) return true;
        return false;
    }

    /// Smallest non-negative maxIterations() of the children, -1 if none
    index_t maxIterations() const
    {
        index_t n = -1;
        for (size_t i = 0; i != m_children.size(); ++i)
        {
            const index_t c = m_children[i]->maxIterations();
            if (c >= 0 && (n < 0 || c < n)) n = c;
        }
        return n;
    }

    std::string name() const { return "Any"; }

    /// The children, in insertion order
    const std::vector<Ptr> & children() const { return m_children; }

private:
    std::vector<Ptr> m_children;
};

template<class T>
typename gsConvergenceTest<T>::Ptr
gsConvergenceTest<T>::makeArcLengthGate(T tolF, T tolU, T tolC, index_t maxIter)
{
    GISMO_ENSURE(maxIter >= 0, "makeArcLengthGate: negative MaxIter " << maxIter);
    memory::shared_ptr< gsAllTest<T> > all(new gsAllTest<T>());
    all->add(Ptr(new gsMinIterTest<T>(1)))
        .add(Ptr(new gsRelResidualTest<T>(tolF, gsResidualReference::LoadNorm)))
        .add(Ptr(new gsRelUpdateTest<T>(tolU)))
        .add(Ptr(new gsPathConstraintTest<T>(tolC)));
    memory::shared_ptr< gsAnyTest<T> > any(new gsAnyTest<T>());
    any->add(Ptr(all)).add(Ptr(new gsMaxIterTest<T>(maxIter)));
    return Ptr(any);
}

} // namespace gismo
