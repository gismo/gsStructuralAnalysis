 /** @file gsStateSpace.h

    @brief Abstract state-space contract of the solver driver, with a plain
    vector space and a lambda-built implementation.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#pragma once
#include <functional>

#include <gsCore/gsLinearAlgebra.h>

namespace gismo
{

template<class T> class gsFunctionalStateSpace;

/**
   @brief Abstract state space: the state of a nonlinear problem that the
   solver driver advances, commits and rolls back, without ever seeing it.

   The state may live on a manifold (e.g. positions and quaternions of a rod).
   The driver works exclusively in the <em>tangent space</em>: increments
   \f$\delta U\f$ and \f$\Delta U\f$ have the <em>tangent size</em>
   \f$n\f$ = tangentSize(), the number of unknowns of the residual. The state
   itself is reached only through the members below.

   Contract:
   1. The driver calls nothing on the state except these members, and no
      \f$U\f$ argument exists in the driver.
   2. Every step and every rejected step begins with rollback(); the driver
      never copies a state.
   3. The state keeps a @e committed and a @e working copy. apply() acts on
      the working copy, commit() makes the working copy the committed one and
      rollback() restores the working copy from the committed one.
   4. tangentSize() == F.size() of the problem (and copyable() where states
      are stored) is checked once at setup.
   5. inverseApplyExact() means that apply(-d, s) directly after apply(d, s)
      restores the working state exactly. Finite differences use that when it
      is true, else export/import. On a non-vector space, differenced
      Jacobians are chart-biased unless pulled back with chartJacobian().
   6. A commit follows the last write-enabled residual evaluation at the
      working state. A gsEvalMode::Probe evaluation never directly precedes
      commit().

   Ops return false on failure instead of throwing.

   @tparam T coefficient type

   \ingroup gsStructuralAnalysis
*/
template<class T>
class gsStateSpace
{
public:
    /// Signature of an apply operation, see apply()
    typedef std::function<bool(const gsVector<T>&, T)> Apply_t;
    /// Signature of a commit operation, see commit()
    typedef std::function<bool()> Commit_t;
    /// Signature of a rollback operation, see rollback()
    typedef std::function<bool()> Rollback_t;
    /// Signature of an export operation, see exportWorking()
    typedef std::function<bool(gsVector<T>&)> Export_t;
    /// Signature of an import operation, see importCommitted()
    typedef std::function<bool(const gsVector<T>&)> Import_t;
    /// Signature of a magnitude operation, see magnitude()
    typedef std::function<bool(gsVector<T>&)> Magnitude_t;
    /// Signature of a chart-Jacobian operation, see chartJacobian()
    typedef std::function<bool(const gsVector<T>&, gsSparseMatrix<T>&)> ChartJacobian_t;

    virtual ~gsStateSpace() {}

    /// Number \f$n\f$ of tangent-space unknowns, the length of every
    /// increment passed to apply() and of the residual.
    virtual index_t tangentSize() const = 0;

    /// Update the working state in place, working <- working (+) s*dU.
    /// @param[in] dU increment of tangent size \f$n\f$
    /// @param[in] s scaling of the increment
    /// @return false on failure (e.g. a size mismatch)
    /// Cost: O(ambient size of the state).
    virtual bool apply(const gsVector<T>& dU, T s) = 0;

    /// Make the working state the committed state, including any
    /// material/history store the client keeps.
    virtual bool commit()   = 0;

    /// Restore the working state from the committed state.
    virtual bool rollback() = 0;

    /// True iff exportWorking() and importCommitted() are implemented.
    virtual bool copyable() const { return false; }

    /// Write the working state to @a U.
    /// @param[out] U resized to the <em>export size</em>, the length of the
    /// exported coordinates. It equals tangentSize() for vector spaces and
    /// may be larger for a manifold state (e.g. positions plus quaternions).
    /// @return false if the space is not copyable()
    virtual bool exportWorking(gsVector<T>& U) const { GISMO_UNUSED(U); return false; }

    /// Set the committed state from @a U and reset the working state to it,
    /// which makes the import idempotent with a following rollback().
    /// @param[in] U coordinates of the <em>export size</em>
    /// @return false if the space is not copyable() or @a U is invalid
    virtual bool importCommitted(const gsVector<T>& U) { GISMO_UNUSED(U); return false; }

    /// Magnitude of the working state, used for round-off estimates.
    /// @param[out] xhat either of length tangentSize() with non-negative
    /// entries (|x| on additive dofs, 0 on rotational ones), or of length 0
    /// meaning that no magnitude is available (round-off floor off).
    /// The default is length 0.
    virtual bool magnitude(gsVector<T>& xhat) const { GISMO_UNUSED(xhat); xhat.resize(0); return true; }

    /// True iff apply(-d, s) directly after apply(d, s) restores the working
    /// state exactly (bitwise or by construction).
    virtual bool inverseApplyExact() const { return false; }

    /// Pull a Jacobian back to the chart of the unperturbed working state.
    /// @param[in] d perturbation (tangent size) at which @a K was assembled
    /// @param[in,out] K Jacobian assembled at working (+) d; replaced by its
    /// pullback to the chart of the unperturbed working state
    /// @return false if the space cannot do it
    virtual bool chartJacobian(const gsVector<T>& d, gsSparseMatrix<T>& K) const
    { GISMO_UNUSED(d); GISMO_UNUSED(K); return false; }

    /// @brief Build a state space from client closures.
    ///
    /// Optional ops that are left empty take the conservative base defaults
    /// of this class: not copyable, no magnitude, no exact inverse, no chart
    /// pullback. They are not identity/vector-space defaults, because a
    /// lambda-built space is typically a manifold, and identity defaults
    /// would silently give chart-biased finite-difference Jacobians.
    /// Plain vector states use gsVectorStateSpace.
    ///
    /// Typical use: <tt>auto S = gsStateSpace<real_t>::fromOps(n, apply, commit, rollback);</tt>
    /// and pass @c S as a <tt>gsStateSpace<real_t>&</tt>.
    ///
    /// @param[in] n tangent size
    /// @param[in] apply, commit, rollback required ops
    /// @param[in] exportWorking, importCommitted either both or neither
    /// @param[in] magnitude optional
    /// @param[in] inverseApplyExact see inverseApplyExact()
    /// @param[in] chartJacobian optional
    static gsFunctionalStateSpace<T> fromOps(index_t n, Apply_t apply, Commit_t commit, Rollback_t rollback,
                                             Export_t exportWorking = Export_t(), Import_t importCommitted = Import_t(),
                                             Magnitude_t magnitude = Magnitude_t(), bool inverseApplyExact = false,
                                             ChartJacobian_t chartJacobian = ChartJacobian_t());
};

/**
   @brief The additive state space \f$\mathbb{R}^n\f$, where
   apply() is \f$U \leftarrow U + s\,dU\f$.

   Export size equals tangent size. Because floating-point
   \f$(U + s d) - s d\f$ is not bitwise \f$U\f$, inverseApplyExact() is false;
   the space is copyable, so finite differences take the exact export/import
   route. The chart of \f$\mathbb{R}^n\f$ is the identity, so chartJacobian()
   leaves the Jacobian unchanged.

   @tparam T coefficient type

   \ingroup gsStructuralAnalysis
*/
template<class T>
class gsVectorStateSpace : public gsStateSpace<T>
{
public:
    /// Zero state of size @a n (committed and working)
    explicit gsVectorStateSpace(index_t n)
    {
        GISMO_ENSURE(n >= 0, "gsVectorStateSpace: negative size " << n);
        m_work      = gsVector<T>::Zero(n);
        m_committed = gsVector<T>::Zero(n);
    }

    /// State initialised to @a U0 (committed and working)
    explicit gsVectorStateSpace(const gsVector<T>& U0)
    : m_work(U0), m_committed(U0)
    {}

    index_t tangentSize() const { return m_work.size(); }

    /// Returns false if dU.size() differs from tangentSize()
    bool apply(const gsVector<T>& dU, T s)
    {
        if (dU.size() != m_work.size())
            return false;
        m_work.noalias() += s * dU;
        return true;
    }

    bool commit()   { m_committed = m_work; return true; }
    bool rollback() { m_work = m_committed; return true; }

    bool copyable() const { return true; }

    bool exportWorking(gsVector<T>& U) const { U = m_work; return true; }

    /// Sets both committed and working state; returns false if the size
    /// differs from tangentSize()
    bool importCommitted(const gsVector<T>& U)
    {
        if (U.size() != m_work.size())
            return false;
        m_committed = U;
        m_work      = U;
        return true;
    }

    /// |working|, of length tangentSize()
    bool magnitude(gsVector<T>& xhat) const { xhat = m_work.cwiseAbs(); return true; }

    bool inverseApplyExact() const { return false; }

    /// The chart of \f$\mathbb{R}^n\f$ is the identity: @a K is unchanged
    bool chartJacobian(const gsVector<T>& d, gsSparseMatrix<T>& K) const
    { GISMO_UNUSED(d); GISMO_UNUSED(K); return true; }

    /// Working state
    const gsVector<T>& working() const { return m_work; }
    /// Committed state
    const gsVector<T>& committed() const { return m_committed; }

private:
    gsVector<T> m_work, m_committed;
};

/**
   @brief State space assembled from client closures, for manifold states
   driven by lambdas.

   Missing optional ops are bound once, at construction, to the base-class
   defaults of gsStateSpace (not copyable, no magnitude, no chart pullback);
   no member tests a closure for emptiness. See gsStateSpace::fromOps().

   @tparam T coefficient type

   \ingroup gsStructuralAnalysis
*/
template<class T>
class gsFunctionalStateSpace : public gsStateSpace<T>
{
public:
    typedef typename gsStateSpace<T>::Apply_t         Apply_t;
    typedef typename gsStateSpace<T>::Commit_t        Commit_t;
    typedef typename gsStateSpace<T>::Rollback_t      Rollback_t;
    typedef typename gsStateSpace<T>::Export_t        Export_t;
    typedef typename gsStateSpace<T>::Import_t        Import_t;
    typedef typename gsStateSpace<T>::Magnitude_t     Magnitude_t;
    typedef typename gsStateSpace<T>::ChartJacobian_t ChartJacobian_t;

    /// @brief Constructor; throws (GISMO_ENSURE) if @a n is negative, if
    /// @a apply, @a commit or @a rollback is empty, or if exactly one of
    /// @a exportWorking and @a importCommitted is given.
    /// @param[in] n tangent size
    /// @param[in] apply, commit, rollback required ops
    /// @param[in] exportWorking, importCommitted either both or neither; the
    /// space is copyable() iff both are given
    /// @param[in] magnitude optional; default is length 0 (no magnitude)
    /// @param[in] inverseApplyExact see gsStateSpace::inverseApplyExact()
    /// @param[in] chartJacobian optional; default is "not available"
    gsFunctionalStateSpace(index_t n, Apply_t apply, Commit_t commit, Rollback_t rollback,
                           Export_t exportWorking = Export_t(), Import_t importCommitted = Import_t(),
                           Magnitude_t magnitude = Magnitude_t(), bool inverseApplyExact = false,
                           ChartJacobian_t chartJacobian = ChartJacobian_t())
    : m_n(n),
      m_apply(apply), m_commit(commit), m_rollback(rollback),
      m_export(exportWorking), m_import(importCommitted),
      m_magnitude(magnitude), m_chart(chartJacobian),
      m_copyable(false), m_exactInverse(inverseApplyExact)
    {
        GISMO_ENSURE(n >= 0, "gsFunctionalStateSpace: negative size " << n);
        GISMO_ENSURE(bool(apply) && bool(commit) && bool(rollback),
                     "gsFunctionalStateSpace: apply, commit and rollback are required");
        GISMO_ENSURE(bool(exportWorking) == bool(importCommitted),
                     "gsFunctionalStateSpace: exportWorking and importCommitted must be given together");

        m_copyable = bool(exportWorking);
        if (!m_copyable)
        {
            m_export = [](gsVector<T>&) { return false; };
            m_import = [](const gsVector<T>&) { return false; };
        }
        if (!m_magnitude)
            m_magnitude = [](gsVector<T>& xhat) { xhat.resize(0); return true; };
        if (!m_chart)
            m_chart = [](const gsVector<T>&, gsSparseMatrix<T>&) { return false; };
    }

    index_t tangentSize() const { return m_n; }
    bool apply(const gsVector<T>& dU, T s) { return m_apply(dU, s); }
    bool commit()   { return m_commit(); }
    bool rollback() { return m_rollback(); }
    bool copyable() const { return m_copyable; }
    bool exportWorking(gsVector<T>& U) const { return m_export(U); }
    bool importCommitted(const gsVector<T>& U) { return m_import(U); }
    bool magnitude(gsVector<T>& xhat) const { return m_magnitude(xhat); }
    bool inverseApplyExact() const { return m_exactInverse; }
    bool chartJacobian(const gsVector<T>& d, gsSparseMatrix<T>& K) const { return m_chart(d, K); }

private:
    index_t m_n;
    Apply_t m_apply;
    Commit_t m_commit;
    Rollback_t m_rollback;
    Export_t m_export;
    Import_t m_import;
    Magnitude_t m_magnitude;
    ChartJacobian_t m_chart;
    bool m_copyable;
    bool m_exactInverse;
};

template<class T>
gsFunctionalStateSpace<T>
gsStateSpace<T>::fromOps(index_t n, Apply_t apply, Commit_t commit, Rollback_t rollback,
                         Export_t exportWorking, Import_t importCommitted,
                         Magnitude_t magnitude, bool inverseApplyExact,
                         ChartJacobian_t chartJacobian)
{
    return gsFunctionalStateSpace<T>(n, apply, commit, rollback, exportWorking, importCommitted,
                                     magnitude, inverseApplyExact, chartJacobian);
}

} // namespace gismo
