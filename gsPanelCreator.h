/** @file gsPanelCreator.h

    @brief Provides declaration of the NurbsCreator struct.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/


// Note:
// will provide functions to make bsplines, nurbs etc
// like the identity=0, ruled, l spapes, rings, etc
// linear extrude .....
// sweep ( eg. half-cylinder by sweeping half-circle
// revolve operation

#pragma once

#include <gsCore/gsMultiPatch.h>

namespace gismo
{

/**
   @brief Class gsPanelCreator provides some simple examples of Nurbs Geometries

   \ingroup Nurbs
*/

template<class T>
struct gsPanelCreator
{
public:
    static gsMultiPatch<T> Plate(T const & Lp, T const & Wp, T const & x=0, T const & y=0, T const & z=0);
    static gsMultiPatch<T> Strip(T const &  Lb, T const &  Hw, T const & x=0, T const & y=0, T const & z=0);
    static gsMultiPatch<T> IBeam(T const & Lb, T const & Hw, T const & Wf, T const & x=0, T const & y=0, T const & z=0);
    static gsMultiPatch<T> TBeam(T const & Lb, T const & Hw, T const & Wf, T const & x=0, T const & y=0, T const & z=0);
    static gsMultiPatch<T> LBeam(T const & Lb, T const & Hw, T const & Wf, T const & x=0, T const & y=0, T const & z=0);
    static gsMultiPatch<T> PanelStrip(T const & Lp, T const & Wp, T const & Hw, T const & x=0, T const & y=0, T const & z=0);
    static gsMultiPatch<T> PanelT(T const & Lp, T const & Wp, T const & Hw, T const & Wf, T const & x=0, T const & y=0, T const & z=0);
    static gsMultiPatch<T> PanelL(T const & Lp, T const & Wp, T const & Hw, T const & Wf, T const & x=0, T const & y=0, T const & z=0);
    static gsMultiPatch<T> PlateGirderL(T const & Lp, T const & Wp, T const & Hwg, T const & Wfg, T const & Hws, T const & Wfs, T const & x=0, T const & y=0, T const & z=0);
}; // struct

#ifdef GISMO_WITH_PYBIND11

    /**
     * @brief Initializes the Python wrapper for the class: gsBoundaryConditions
     */
    void pybind11_init_gsPanelCreator(pybind11::module &m);

#endif

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsPanelCreator.hpp)
#endif
