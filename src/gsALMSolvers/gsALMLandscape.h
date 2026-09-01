 /** @file gsALMLandscape.h

    @brief Container for the connected bifurcation diagram (solution landscape)
    produced by an arc-length exploration.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsMultiPatch.h>

namespace gismo
{

/**
    @brief Stores a connected bifurcation diagram as a set of curves.

    A landscape is a list of Curves. Each Curve is an ordered list of Points
    (a solution vector, its continuation parameter, stability and bifurcation
    flags, plus an optional deformed geometry). Curves carry parent-curve /
    branch-point connectivity so that the whole diagram forms a tree rooted at
    one (or more) seed curves.

    This is a pure in-memory container plus CSV/Paraview writers: it contains
    no solver logic.

    \tparam T coefficient type

    \ingroup gsALMSolvers
*/
template <class T>
class gsALMLandscape
{
public:

    /// @brief A single point on a curve of the landscape.
    struct Point
    {
        gsVector<T>     U;              ///< free-DOF solution vector
        T               L;              ///< load / continuation parameter
        index_t         stability;      ///< +1 stable, -1 unstable (caller's sign convention)
        /// @brief Tangent inertia at this point: the number of negative
        /// eigenvalues/pivots of \f$K_T\f$, as reported by the solver's
        /// `negatives()`. \c -1 means NOT RECORDED.
        ///
        /// The sentinel cannot be 0 the way it is for \c stability: 0 is a
        /// VALID inertia (a stable point has zero negative pivots), so the
        /// "unset" value has to sit outside the range of the quantity. \c -1
        /// is already the "no inertia recorded yet" sentinel used by
        /// gsALMExploration's `negPrevAccepted`.
        ///
        /// Storing the COUNT (not the +-1 sign in \c stability) is what lets a
        /// reader distinguish a 0->1 crossing from a 1->2 secondary one.
        index_t         negatives;
        bool            isBifurcation;  ///< marked singular point
        /// @brief True iff this ROW satisfies the residual test that was used to
        /// ACCEPT it. An ordinary corrector point, a limit point and a converged
        /// extended-system singular-point solve (status Success) all set this
        /// true -- each met its own termination test. A singular-point solve
        /// that did NOT converge is stored with this false (see \c unresolved,
        /// which then also reads true).
        ///
        /// This is NOT a per-row certificate of the STRENGTH of that test. The
        /// extended-system solve's basic termination test is
        /// \f$\|K_T V\|\to 0\f$, which certifies SINGULARITY only; the solver
        /// option \c SingularPointComposite adds the equilibrium residuals
        /// (TolF/TolU) to it. Whether a PARTICULAR solve used the stronger
        /// composite test is a RUN-LEVEL property of that option, not encoded
        /// per row here -- a consumer that needs the stronger certificate must
        /// read the caller's own \c SingularPointComposite setting, not this
        /// flag.
        ///
        /// The replacement filter for "give me only certified singular points":
        /// a stored extended-solve point is identified by
        /// \c isBifurcation \c == \c true \c && \c stability \c == \c 0;
        /// \c unresolved \c == \c true marks a point no singular-point solve
        /// resolved at all (its location is only the traced, post-crossing
        /// point). Consumers that only need the location of a crossing (plots,
        /// connectivity) may use every point regardless of this flag.
        bool            equilibrium;
        /// @brief True when this point was FLAGGED as a singular point that the
        /// singular-point solve could NOT resolve: its location is only the traced
        /// (post-crossing) point, no refined (U*,L*) exists, and no branch was followed.
        /// Such a point has isBifurcation == true and equilibrium == false. A CERTIFIED
        /// singular point has unresolved == false. Consumers must never treat the two
        /// as the same thing.
        bool            unresolved;
        /// @brief Lambda at the pre-crossing end of a FAILED localization's detection
        /// interval (\c Lold at the \c markUnresolvedSingular call site); meaningful
        /// only when \c bracketProbes \c >= \c 0, see there. \c 0 where the record is
        /// absent.
        T               bracketLo;
        /// @brief Lambda at the post-crossing end of a FAILED localization's detection
        /// interval (\c Lcur at the \c markUnresolvedSingular call site); meaningful
        /// only when \c bracketProbes \c >= \c 0, see there. \c 0 where the record is
        /// absent.
        T               bracketHi;
        /// @brief Total bisection probes spent, over the first localization attempt
        /// and every retry, at this crossing. \c -1 means NOT RECORDED -- either the
        /// point is not an unresolved localization failure, or it predates this field.
        /// \c bracketProbes \c >= \c 0 is the single discriminator that says
        /// \c bracketLo / \c bracketHi are meaningful.
        ///
        /// This is a HOOK for a later resume API to re-attempt the crossing instead of
        /// rediscovering it -- no such consumer exists yet.
        index_t         bracketProbes;
        /// @brief OBSERVATION, not inference: the tangent inertia reported by the
        /// solver (\c negatives()) changed by MORE than one between the two accepted
        /// points bracketing this singular point, and \c gsALMExploration's crossing
        /// detection reduced that multi-step change to a SINGLE simple crossing --
        /// one localized crossing, at most one branch followed (see
        /// \c _localizeCrossing's bracket-narrowing comment). The reported inertia
        /// count need not be exact: with the eigenvalue stability method it comes
        /// from an inverse-power iteration. The two usual causes are a genuine
        /// multiplicity >= 2 crossing, and two distinct simple crossings merged by a
        /// coarse arc-length step.
        ///
        /// Consequence for a consumer: the landscape is INCOMPLETE over the interval
        /// that produced this point -- one or more branches emanating somewhere in
        /// that interval (not necessarily at this point) are not present. Real
        /// multiplicity >= 2 branch switching is not implemented (
        /// \c gsALMBase::computeBranchTangent resolves the emanating tangent only at
        /// a simple, multiplicity-1, branch point).
        ///
        /// \c false means either a simple crossing or a point that was never a
        /// crossing -- it is NOT a certificate that the crossing was simple in the
        /// exact problem, only that the observed accepted-point inertia change was
        /// by one.
        ///
        /// Set in-class (default member initializer) rather than by \c addPoint,
        /// which this flag is not among the fields of: \c gsALMExploration writes it
        /// directly on the stored point once a multi-step crossing is detected.
        ///
        /// NOT written by \c writeCsv -- the CSV schema is frozen and carries no
        /// column for it, so a reader of the CSV cannot see this flag. It IS
        /// persisted by \c saveHDF5 / \c loadHDF5, following the same pattern as
        /// \c bracketLo / \c bracketHi / \c bracketProbes above.
        bool            multiplicityReduced = false;
        gsMultiPatch<T> deformed;       ///< OPTIONAL deformed geometry; empty (0 patches) if not provided
    };

    /// @brief An ordered list of points forming one branch of the landscape.
    struct Curve
    {
        std::vector<Point> points;      ///< ordered points along the curve
        index_t parentCurve;            ///< index of parent curve, -1 for the seed curve
        index_t parentPointIdx;         ///< index of the branch point in the parent curve, -1 for seed
    };

    /// @brief Default constructor: an empty landscape.
    gsALMLandscape() = default;

    /// @brief Starts a new curve; returns its id.
    index_t addCurve(index_t parentCurve = -1, index_t parentPointIdx = -1);

    /// @brief Appends a point to curve \a curveId. \a deformed may be nullptr (no geometry stored).
    /// \a equilibrium defaults to true: it means "this row satisfies the residual
    /// test that was used to accept it", true for the ordinary corrector path and
    /// for a converged extended-system singular-point solve alike; pass false only
    /// for a point no solve accepted at all (see Point::equilibrium). \a unresolved
    /// is always initialised false here -- use \c markUnresolvedSingular to
    /// relabel an already-added point as an unresolved singular point.
    /// \a negatives is the tangent inertia count; it defaults to -1 ("not
    /// recorded"). NOTE the parameter order deliberately does NOT match the CSV
    /// column order (where `negatives` sits right after `stability`): appending
    /// keeps every existing call site source-compatible. Do not "fix" this.
    void addPoint(index_t curveId, const gsVector<T> & U, T L, index_t stability,
                  const gsMultiPatch<T> * deformed = nullptr, bool isBifurcation = false,
                  bool equilibrium = true, index_t negatives = -1);

    /// @brief Marks the LAST point of \a curveId as a bifurcation point.
    void markBifurcation(index_t curveId);

    /// @brief Marks the LAST point of \a curveId as an UNRESOLVED singular point:
    /// isBifurcation = true, equilibrium = false, unresolved = true, atomically.
    /// Used when a singular-point solve failed, so the flagged point is the traced
    /// post-crossing point and no refined (U*,L*) exists. Its \c stability is left
    /// untouched (it is a genuine converged equilibrium whose inertia WAS measured).
    ///
    /// \a bracketLo, \a bracketHi, \a bracketProbes are OPTIONAL: they record a
    /// FAILED localization's detection interval and total probe count, as a hook for
    /// a later resume to re-attempt the crossing (see Point::bracketProbes). Left at
    /// their defaults, \a bracketProbes stays -1 ("not recorded") -- e.g. for the
    /// classification-failure and singular-solve-failure call sites, whose crossing
    /// DID localize and whose bracket is therefore not the failed one this hook is
    /// about.
    void markUnresolvedSingular(index_t curveId,
                                T bracketLo = (T)0, T bracketHi = (T)0,
                                index_t bracketProbes = -1);

    /// @brief Removes the last curve added (used by dedup abort).
    void removeLastCurve();

    /// @brief Number of curves in the landscape.
    size_t nCurves() const;

    /// @brief Total number of points over all curves.
    size_t nPoints() const;

    /// @brief Const access to curve \a i.
    const Curve & curve(index_t i) const;

    /// @brief Mutable access to curve \a i.
    Curve & curve(index_t i);

    /// @brief (childCurve, parentCurve) edges, skipping seed curves (parent==-1).
    std::vector<std::pair<index_t,index_t> > connectivity() const;

    /// @brief Indices of bifurcation points of curve \a i.
    /// @note The returned indices MIX certified singular points and unresolved
    /// ones (see Point::unresolved): a caller that needs only certified singular
    /// points must additionally filter on \c !points[idx].unresolved.
    std::vector<index_t> bifurcationIndices(index_t i) const;

    /// @brief Writes the landscape to a CSV file, one row per point:
    /// curve,point,L,normU,stability,negatives,isBifurcation,parentCurve,parentPointIdx,equilibrium,unresolved
    /// (header row first; parent columns repeated per row for greppability). The
    /// \c negatives column is the tangent inertia count (see Point::negatives); \c -1
    /// means not recorded. The \c equilibrium column is 0 for a row that did NOT
    /// satisfy the residual test used to accept it -- see Point::equilibrium for
    /// what that test was and was not. The trailing \c unresolved column is 1 only
    /// for a point flagged as a singular point that no singular-point solve could
    /// resolve (see Point::unresolved / markUnresolvedSingular); a reader that
    /// wants only certified singular points must filter on \c isBifurcation \c ==
    /// \c 1 \c && \c unresolved \c == \c 0.
    void writeCsv(const std::string & fname) const;

    /// @brief Writes one .vts per point that HAS a stored deformed geometry
    /// (basename_c<curve>_p<point>) plus a gsParaviewCollection <basename>.pvd over them.
    /// Points without geometry are skipped silently. No-op if nothing has geometry.
    void writeParaview(const std::string & basename, index_t npts = 1000) const;

#ifdef gsHDF5_ENABLED
    /// @brief Schema version of the HDF5 checkpoint layout written by saveHDF5().
    /// Stored as the integer attribute "gsALMLandscapeSchemaVersion" on the file's
    /// root group. loadHDF5() rejects any file whose marker differs or is absent.
    /// Bump this whenever a dataset in saveHDF5() is added, removed or renamed.
    static index_t hdf5SchemaVersion() { return 1; }

    /// @brief Saves the FULL landscape (all curves/points, including a stored
    /// deformed gsMultiPatch for every point that has one) to an HDF5 file.
    /// Whole-file rewrite (H5F_ACC_TRUNC); intended as a cheap crash checkpoint.
    /// A second, brief write phase reopens the file to stamp the root-group
    /// attribute "gsALMLandscapeSchemaVersion" (see hdf5SchemaVersion()) once the
    /// first phase's handle has closed. Only available in gsHDF5-enabled builds.
    /// Lets gsHDF5 errors propagate.
    void saveHDF5(const std::string & fname) const;

    /// @brief Replaces the current contents with a landscape previously written
    /// by saveHDF5. After load the container compares equal to the saved one
    /// (curves, points, U/L/stability/isBifurcation, parent links, geometries).
    /// Throws (GISMO_ENSURE) before any field is read if the file lacks the
    /// "gsALMLandscapeSchemaVersion" root-group attribute or carries a value
    /// other than hdf5SchemaVersion(). Only available in gsHDF5-enabled builds.
    void loadHDF5(const std::string & fname);
#endif

protected:

    /// The curves making up the landscape.
    std::vector<Curve> m_curves;

}; // class gsALMLandscape


} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsALMLandscape.hpp)
#endif
