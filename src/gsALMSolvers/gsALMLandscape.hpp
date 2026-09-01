 /** @file gsALMLandscape.hpp

    @brief Implementation of gsALMLandscape.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#pragma once

#include <fstream>

#include <gsCore/gsMemory.h>
#include <gsIO/gsFileManager.h>
#include <gsIO/gsWriteParaview.h>
#include <gsIO/gsParaviewCollection.h>

#ifdef gsHDF5_ENABLED
#include <cstdio>          // std::snprintf for 4-digit zero-padded dataset names
#include <cmath>           // std::lround
#include <gsHDF5/gsHDF5.h>
#endif

namespace gismo
{

template <class T>
index_t gsALMLandscape<T>::addCurve(index_t parentCurve, index_t parentPointIdx)
{
    Curve c;
    c.parentCurve    = parentCurve;
    c.parentPointIdx = parentPointIdx;
    m_curves.push_back(give(c));
    return static_cast<index_t>(m_curves.size()) - 1;
}

template <class T>
void gsALMLandscape<T>::addPoint(index_t curveId, const gsVector<T> & U, T L, index_t stability,
                                 const gsMultiPatch<T> * deformed, bool isBifurcation,
                                 bool equilibrium, index_t negatives)
{
    GISMO_ASSERT(curveId >= 0 && curveId < static_cast<index_t>(m_curves.size()),
                 "gsALMLandscape::addPoint: curveId "<<curveId<<" out of range [0,"<<m_curves.size()<<").");

    Point p;
    p.U             = U;
    p.L             = L;
    p.stability     = stability;
    p.negatives     = negatives;
    p.isBifurcation = isBifurcation;
    p.equilibrium   = equilibrium;
    p.unresolved    = false;   // set true only by markUnresolvedSingular
    p.bracketLo     = (T)0;    // set only by markUnresolvedSingular
    p.bracketHi     = (T)0;    // set only by markUnresolvedSingular
    p.bracketProbes = -1;      // NOT RECORDED; set only by markUnresolvedSingular
    if (deformed != nullptr)
        p.deformed = *deformed;   // copy the geometry; otherwise leave default (0 patches)

    m_curves[curveId].points.push_back(give(p));
}

template <class T>
void gsALMLandscape<T>::markBifurcation(index_t curveId)
{
    GISMO_ASSERT(curveId >= 0 && curveId < static_cast<index_t>(m_curves.size()),
                 "gsALMLandscape::markBifurcation: curveId "<<curveId<<" out of range [0,"<<m_curves.size()<<").");
    GISMO_ASSERT(!m_curves[curveId].points.empty(),
                 "gsALMLandscape::markBifurcation: curve "<<curveId<<" has no points.");
    m_curves[curveId].points.back().isBifurcation = true;
}

template <class T>
void gsALMLandscape<T>::markUnresolvedSingular(index_t curveId,
                                               T bracketLo, T bracketHi, index_t bracketProbes)
{
    GISMO_ASSERT(curveId >= 0 && curveId < static_cast<index_t>(m_curves.size()),
                 "gsALMLandscape::markUnresolvedSingular: curveId "<<curveId<<" out of range [0,"<<m_curves.size()<<").");
    GISMO_ASSERT(!m_curves[curveId].points.empty(),
                 "gsALMLandscape::markUnresolvedSingular: curve "<<curveId<<" has no points.");
    Point & p = m_curves[curveId].points.back();
    p.isBifurcation = true;
    p.equilibrium   = false;
    p.unresolved    = true;
    p.bracketLo     = bracketLo;
    p.bracketHi     = bracketHi;
    p.bracketProbes = bracketProbes;
}

template <class T>
void gsALMLandscape<T>::removeLastCurve()
{
    GISMO_ASSERT(!m_curves.empty(), "gsALMLandscape::removeLastCurve: no curves to remove.");
    m_curves.pop_back();
}

template <class T>
size_t gsALMLandscape<T>::nCurves() const
{
    return m_curves.size();
}

template <class T>
size_t gsALMLandscape<T>::nPoints() const
{
    size_t n = 0;
    for (const Curve & c : m_curves)
        n += c.points.size();
    return n;
}

template <class T>
const typename gsALMLandscape<T>::Curve & gsALMLandscape<T>::curve(index_t i) const
{
    GISMO_ASSERT(i >= 0 && i < static_cast<index_t>(m_curves.size()),
                 "gsALMLandscape::curve: index "<<i<<" out of range [0,"<<m_curves.size()<<").");
    return m_curves[i];
}

template <class T>
typename gsALMLandscape<T>::Curve & gsALMLandscape<T>::curve(index_t i)
{
    GISMO_ASSERT(i >= 0 && i < static_cast<index_t>(m_curves.size()),
                 "gsALMLandscape::curve: index "<<i<<" out of range [0,"<<m_curves.size()<<").");
    return m_curves[i];
}

template <class T>
std::vector<std::pair<index_t,index_t> > gsALMLandscape<T>::connectivity() const
{
    std::vector<std::pair<index_t,index_t> > edges;
    for (size_t i = 0; i != m_curves.size(); ++i)
        if (m_curves[i].parentCurve != -1)
            edges.push_back( std::make_pair(static_cast<index_t>(i), m_curves[i].parentCurve) );
    return edges;
}

template <class T>
std::vector<index_t> gsALMLandscape<T>::bifurcationIndices(index_t i) const
{
    GISMO_ASSERT(i >= 0 && i < static_cast<index_t>(m_curves.size()),
                 "gsALMLandscape::bifurcationIndices: index "<<i<<" out of range [0,"<<m_curves.size()<<").");
    std::vector<index_t> idx;
    const Curve & c = m_curves[i];
    for (size_t p = 0; p != c.points.size(); ++p)
        if (c.points[p].isBifurcation)
            idx.push_back(static_cast<index_t>(p));
    return idx;
}

template <class T>
void gsALMLandscape<T>::writeCsv(const std::string & fname) const
{
    std::ofstream file(fname.c_str());
    GISMO_ENSURE(file.is_open(), "gsALMLandscape::writeCsv: could not open "<<fname<<" for writing.");

    file << "curve,point,L,normU,stability,negatives,isBifurcation,parentCurve,parentPointIdx,equilibrium,unresolved\n";
    for (size_t i = 0; i != m_curves.size(); ++i)
    {
        const Curve & c = m_curves[i];
        for (size_t p = 0; p != c.points.size(); ++p)
        {
            const Point & pt = c.points[p];
            file << i << ","
                 << p << ","
                 << pt.L << ","
                 << pt.U.norm() << ","
                 << pt.stability << ","
                 << pt.negatives << ","
                 << (pt.isBifurcation ? 1 : 0) << ","
                 << c.parentCurve << ","
                 << c.parentPointIdx << ","
                 << (pt.equilibrium ? 1 : 0) << ","
                 << (pt.unresolved ? 1 : 0) << "\n";
        }
    }
    file.close();
}

template <class T>
void gsALMLandscape<T>::writeParaview(const std::string & basename, index_t npts) const
{
    // Bail out (no .pvd produced) if no point carries a stored geometry.
    bool anyGeometry = false;
    for (const Curve & c : m_curves)
        for (const Point & pt : c.points)
            if (pt.deformed.nPatches() != 0) { anyGeometry = true; break; }

    if (!anyGeometry)
        return;

    // Ensure the target directory exists (basename may contain a path).
    gsFileManager::mkdir(gsFileManager::getPath(basename));

    gsParaviewCollection collection(basename);

    for (size_t i = 0; i != m_curves.size(); ++i)
    {
        const Curve & c = m_curves[i];
        for (size_t p = 0; p != c.points.size(); ++p)
        {
            const Point & pt = c.points[p];
            if (pt.deformed.nPatches() == 0)
                continue;

            const std::string fnBase = basename + "_c" + util::to_string(i) + "_p" + util::to_string(p);
            const std::string fnBase_nopath = gsFileManager::getFilename(fnBase);

            // gsWriteParaview writes one file per patch as fnBase_<k>.<ext>.
            gsWriteParaview<T>(pt.deformed, fnBase, static_cast<unsigned>(npts));

            for (size_t k = 0; k != pt.deformed.nPatches(); ++k)
            {
                const std::string ext =
                    (pt.deformed.patch(k).domainDim() == 1) ? ".vtp" : ".vts";
                collection.addPart(fnBase_nopath + "_" + util::to_string(k) + ext);
            }
        }
    }

    collection.save();
}

#ifdef gsHDF5_ENABLED

namespace {
/// Formats a flat HDF5 dataset/group name with 4-digit zero-padded indices.
/// (kept file-local; header-inline, so the anonymous namespace avoids ODR clashes)
inline std::string almCurveName(index_t c, const char * suffix)
{
    char buf[64];
    std::snprintf(buf, sizeof(buf), "curve_%04d_%s", (int)c, suffix);
    return std::string(buf);
}
inline std::string almPointMpName(index_t c, index_t p)
{
    char buf[64];
    std::snprintf(buf, sizeof(buf), "curve_%04d_point_%04d_mp", (int)c, (int)p);
    return std::string(buf);
}
/// Name of the root-group attribute holding the checkpoint schema version.
const char * const almSchemaAttrName = "gsALMLandscapeSchemaVersion";
} // anonymous namespace

template <class T>
void gsALMLandscape<T>::saveHDF5(const std::string & fname) const
{
    {
        // Whole-file rewrite: gsHDF5Writer opens with H5F_ACC_TRUNC, there is no
        // append. All top-level names are created at the file root (flat scheme).
        gsHDF5Writer writer(fname);

        const index_t nc = static_cast<index_t>(m_curves.size());

        // structure: row c = [nPts_c, parentCurve_c, parentPointIdx_c]. Integer
        // values are stored as T; the reader rounds them back with std::lround.
        // An empty landscape is stored as a 0x3 matrix (verified accepted by the
        // writer/HDF5 as a zero-extent dataspace; loadHDF5 reads it back as 0 rows).
        gsMatrix<T> structure(nc, 3);
        for (index_t c = 0; c != nc; ++c)
        {
            structure(c,0) = static_cast<T>(m_curves[c].points.size());
            structure(c,1) = static_cast<T>(m_curves[c].parentCurve);
            structure(c,2) = static_cast<T>(m_curves[c].parentPointIdx);
        }
        writer.write("structure", structure);

        for (index_t c = 0; c != nc; ++c)
        {
            const Curve & cur = m_curves[c];
            const index_t nPts = static_cast<index_t>(cur.points.size());
            if (nPts == 0)
                continue; // structure records nPts==0; nothing else to store

            // U: numDof x nPts, column j = point j's free-DOF vector.
            const index_t numDof = cur.points[0].U.size();
            gsMatrix<T> U(numDof, nPts);
            gsMatrix<T> L(nPts, 1), stab(nPts, 1), neg(nPts, 1), bif(nPts, 1), hasgeom(nPts, 1), eq(nPts, 1), unres(nPts, 1);
            gsMatrix<T> bLo(nPts, 1), bHi(nPts, 1), bProbes(nPts, 1), multRed(nPts, 1);
            for (index_t p = 0; p != nPts; ++p)
            {
                const Point & pt = cur.points[p];
                if (pt.U.size() == numDof)
                    U.col(p) = pt.U;
                else
                    U.col(p).setZero(); // defensive: ragged curve (should not occur)
                L(p,0)       = pt.L;
                stab(p,0)    = static_cast<T>(pt.stability);
                neg(p,0)     = static_cast<T>(pt.negatives);
                bif(p,0)     = pt.isBifurcation ? static_cast<T>(1) : static_cast<T>(0);
                hasgeom(p,0) = (pt.deformed.nPatches() != 0) ? static_cast<T>(1)
                                                             : static_cast<T>(0);
                eq(p,0)      = pt.equilibrium ? static_cast<T>(1) : static_cast<T>(0);
                unres(p,0)   = pt.unresolved  ? static_cast<T>(1) : static_cast<T>(0);
                bLo(p,0)     = pt.bracketLo;
                bHi(p,0)     = pt.bracketHi;
                bProbes(p,0) = static_cast<T>(pt.bracketProbes);
                multRed(p,0) = pt.multiplicityReduced ? static_cast<T>(1) : static_cast<T>(0);
            }
            writer.write(almCurveName(c,"U"),         U);
            writer.write(almCurveName(c,"L"),         L);
            writer.write(almCurveName(c,"stability"), stab);
            writer.write(almCurveName(c,"negatives"), neg);
            writer.write(almCurveName(c,"bif"),       bif);
            writer.write(almCurveName(c,"hasgeom"),   hasgeom);
            writer.write(almCurveName(c,"equilibrium"), eq);
            writer.write(almCurveName(c,"unresolved"),  unres);
            writer.write(almCurveName(c,"bracketLo"),     bLo);
            writer.write(almCurveName(c,"bracketHi"),     bHi);
            writer.write(almCurveName(c,"bracketProbes"), bProbes);
            writer.write(almCurveName(c,"multiplicityReduced"), multRed);

            // One gsMultiPatch group per point that carries a stored geometry.
            for (index_t p = 0; p != nPts; ++p)
                if (cur.points[p].deformed.nPatches() != 0)
                    writer.write(almPointMpName(c,p), cur.points[p].deformed);
        }
    }   // writer's H5File closes here

    // The schema marker is stamped in a second pass: gsHDF5Writer keeps its
    // H5File private, so there is no way to attach a root attribute while it is
    // open, and HDF5 refuses a second conflicting open on a file still held for
    // writing. This is safe only because gsHDF5Writer's constructor opens with
    // H5F_ACC_TRUNC, recreating the file empty on every save, so the reopened
    // root group never already carries the attribute; if gsHDF5Writer ever
    // gains an append mode, this write needs an attrExists() + removeAttr()
    // guard first.
    H5::H5File file(fname, H5F_ACC_RDWR);
    H5::Group root = file.openGroup("/");
    internal::h5WriteAttr(root, almSchemaAttrName, (int)hdf5SchemaVersion());
}

template <class T>
void gsALMLandscape<T>::loadHDF5(const std::string & fname)
{
    {
        // Probe the schema marker before touching m_curves or opening a
        // gsHDF5Reader, so a rejected load leaves the landscape unchanged and
        // fails before any dataset read that could throw an opaque
        // H5::FileIException instead of an explicit message.
        H5::H5File file(fname, H5F_ACC_RDONLY);
        H5::Group root = file.openGroup("/");
        GISMO_ENSURE(root.attrExists(almSchemaAttrName),
            "gsALMLandscape::loadHDF5: \"" << fname << "\" has no "
            "\"" << almSchemaAttrName << "\" attribute, so it predates HDF5 "
            "checkpoint schema versioning and lacks the bracketLo, bracketHi, "
            "bracketProbes and multiplicityReduced datasets this build's "
            "loadHDF5 (schema version " << hdf5SchemaVersion() << ") requires. "
            "Regenerate the checkpoint with this build.");
        const int found = internal::h5ReadIntAttr(root, almSchemaAttrName);
        GISMO_ENSURE(found == (int)hdf5SchemaVersion(),
            "gsALMLandscape::loadHDF5: \"" << fname << "\" was written with "
            "schema version " << found << ", but this build's loadHDF5 requires "
            "version " << hdf5SchemaVersion() << ". Regenerate the checkpoint "
            "with this build.");
    }   // probe handle closes before the reader opens

    m_curves.clear();

    gsHDF5Reader reader(fname);

    gsMatrix<T> structure;
    reader.read("structure", structure);

    const index_t nc = structure.rows();
    m_curves.resize(nc);

    for (index_t c = 0; c != nc; ++c)
    {
        // Integer fields were stored as T; round back to index_t.
        const index_t nPts = static_cast<index_t>(std::lround(structure(c,0)));
        Curve & cur = m_curves[c];
        cur.parentCurve    = static_cast<index_t>(std::lround(structure(c,1)));
        cur.parentPointIdx = static_cast<index_t>(std::lround(structure(c,2)));

        if (nPts == 0)
            continue;

        gsMatrix<T> U, L, stab, neg, bif, hasgeom, eq, unres;
        gsMatrix<T> bLo, bHi, bProbes, multRed;
        reader.read(almCurveName(c,"U"),         U);
        reader.read(almCurveName(c,"L"),         L);
        reader.read(almCurveName(c,"stability"), stab);
        reader.read(almCurveName(c,"negatives"), neg);
        reader.read(almCurveName(c,"bif"),       bif);
        reader.read(almCurveName(c,"hasgeom"),   hasgeom);
        reader.read(almCurveName(c,"equilibrium"), eq);
        reader.read(almCurveName(c,"unresolved"),  unres);
        // Unconditional, like every field above: loadHDF5 validated the schema
        // version before reaching this point, so a file that lacks these datasets
        // has already been rejected with an explicit message. Adding, removing or
        // renaming any dataset here requires bumping hdf5SchemaVersion().
        reader.read(almCurveName(c,"bracketLo"),     bLo);
        reader.read(almCurveName(c,"bracketHi"),     bHi);
        reader.read(almCurveName(c,"bracketProbes"), bProbes);
        reader.read(almCurveName(c,"multiplicityReduced"), multRed);

        cur.points.resize(nPts);
        for (index_t p = 0; p != nPts; ++p)
        {
            Point & pt = cur.points[p];
            pt.U             = U.col(p);
            pt.L             = L(p,0);
            pt.stability     = static_cast<index_t>(std::lround(stab(p,0)));
            pt.negatives     = static_cast<index_t>(std::lround(neg(p,0)));
            pt.isBifurcation = (std::lround(bif(p,0)) != 0);
            pt.equilibrium   = (std::lround(eq(p,0)) != 0);
            pt.unresolved    = (std::lround(unres(p,0)) != 0);
            pt.bracketLo     = bLo(p,0);
            pt.bracketHi     = bHi(p,0);
            pt.bracketProbes = static_cast<index_t>(std::lround(bProbes(p,0)));
            pt.multiplicityReduced = (std::lround(multRed(p,0)) != 0);
            if (std::lround(hasgeom(p,0)) != 0)
                reader.read(almPointMpName(c,p), pt.deformed);
        }
    }
}

#endif // gsHDF5_ENABLED


} // namespace gismo
