/** @file gsPanelCreator.hpp

    @brief Provides implementation of the NurbsCreator struct.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsMultiPatch.h>
#include <gsNurbs/gsKnotVector.h>
#include <gsNurbs/gsTensorBSpline.h>

namespace gismo
{

/*
   @brief Class gsPanelCreator provides some simple examples of panels with B-splines
*/

template<class T>
gsMultiPatch<T> gsPanelCreator<T>::Plate(T const & Lp, T const &  Wp, T const & x, T const & y, T const & z)
{
   gsMultiPatch<T> result;

   gsKnotVector<T> KV (0,1,0,2) ;
   gsMatrix<T> C(4,3) ;

   C.row(0)<< x      ,y      ,z;
   C.row(1)<< Lp+x   ,y      ,z;
   C.row(2)<< x      ,Wp+y   ,z;
   C.row(3)<< Lp+x   ,Wp+y   ,z;

   result.addPatch(gsTensorBSpline<2,T>(KV,KV, give(C)));
   return result;
   // return TensorBSpline2Ptr(new gsTensorBSpline<2,T>(KV,KV, give(C)));
}

template<class T>
gsMultiPatch<T> gsPanelCreator<T>::Strip(T const & Lb, T const & Hw, T const & x, T const & y, T const & z)
{
   gsMultiPatch<T> result;

   gsKnotVector<T> KV (0,1,0,2) ;
   gsMatrix<T> C(4,3) ;

   C.row(0)<< x,   y   ,z;
   C.row(1)<< Lb+x,y   ,z;
   C.row(2)<< x   ,y   ,Hw+z;
   C.row(3)<< Lb+x,y   ,Hw+z;

   result.addPatch(gsTensorBSpline<2,T>(KV,KV, give(C)));
   return result;
   // return TensorBSpline2Ptr(new );
}

template <class T>
gsMultiPatch<T> gsPanelCreator<T>::IBeam(T const & Lb, T const & Hw, T const & Wf, T const & x, T const & y, T const & z)
{
   gsMultiPatch<T> result;
   gsMultiPatch<T> tmp;

   // Web
   tmp = Strip(Lb,Hw,0,0,-Hw/2.);
   result.addPatch(tmp.patch(0));

   tmp = Plate(Lb,Wf/2.,0,0,Hw/2.);
   result.addPatch(tmp.patch(0));

   tmp = Plate(Lb,Wf/2.,0,-Wf/2.,Hw/2.);
   result.addPatch(tmp.patch(0));

   tmp = Plate(Lb,Wf/2.,0,0,-Hw/2.);
   result.addPatch(tmp.patch(0));

   tmp = Plate(Lb,Wf/2.,0,-Wf/2.,-Hw/2.);
   result.addPatch(tmp.patch(0));

   for (size_t p = 0; p!=result.nPatches(); p++)
   {
      result.patch(p).coefs().col(0).array() += x;
      result.patch(p).coefs().col(1).array() += y;
      result.patch(p).coefs().col(2).array() += z;
   }

   result.computeTopology();
   result.addAutoBoundaries();
   return result;
}

template <class T>
gsMultiPatch<T> gsPanelCreator<T>::TBeam(T const & Lb, T const & Hw, T const & Wf, T const & x, T const & y, T const & z)
{
    gsMultiPatch<T> result;
    gsMultiPatch<T> tmp;

    // Web
    tmp = Strip(Lb,Hw,0,0,0);
    result.addPatch(tmp.patch(0));

    // Flange, left
    tmp = Plate(Lb,Wf/2,0,0,Hw);
    result.addPatch(tmp.patch(0));

    // Flange, right
    tmp = Plate(Lb,Wf/2,0,-Wf/2,Hw);
    result.addPatch(tmp.patch(0));

    for (size_t p = 0; p!=result.nPatches(); p++)
    {
        result.patch(p).coefs().col(0).array() += x;
        result.patch(p).coefs().col(1).array() += y;
        result.patch(p).coefs().col(2).array() += z;
    }

    result.computeTopology();

    // result.addInterface(&result.patch(0),4,&result.patch(1),1);
    // result.addInterface(&result.patch(0),4,&result.patch(2),2);
    // result.addInterface(&result.patch(1),1,&result.patch(2),2);

    result.addAutoBoundaries();
    return result;
}

template <class T>
gsMultiPatch<T> gsPanelCreator<T>::LBeam(T const & Lb, T const & Hw, T const & Wf, T const & x, T const & y, T const & z)
{
    gsMultiPatch<T> result;
    gsMultiPatch<T> tmp;

    // Web
    tmp = Strip(Lb,Hw,0,0,0);
    for (size_t p=0; p!=tmp.nPatches(); p++)
        result.addPatch(tmp.patch(p));

    // Flange, left
    tmp = Plate(Lb,Wf,0,0,Hw);
    for (size_t p=0; p!=tmp.nPatches(); p++)
        result.addPatch(tmp.patch(p));

    for (size_t p = 0; p!=result.nPatches(); p++)
    {
        result.patch(p).coefs().col(0).array() += x;
        result.patch(p).coefs().col(1).array() += y;
        result.patch(p).coefs().col(2).array() += z;
    }

    result.addInterface(&result.patch(0),4,&result.patch(1),1);

    result.addAutoBoundaries();

    return result;
}

template <class T>
gsMultiPatch<T> gsPanelCreator<T>::PanelT(T const & Lp, T const & Wp, T const & Hw, T const & Wf, T const & x, T const & y, T const & z)

{
    gsMultiPatch<T> result;
    gsMultiPatch<T> tmp;

    // Base plate, left
    // Flange, left
    tmp = Plate(Lp,Wp/2,0,0,0);
    for (size_t p=0; p!=tmp.nPatches(); p++)
        result.addPatch(tmp.patch(p));

    // Base plate, right
    tmp = Plate(Lp,Wp/2,0,-Wp/2,0);
    for (size_t p=0; p!=tmp.nPatches(); p++)
        result.addPatch(tmp.patch(p));

    // T-Beam
    gsMultiPatch<> beam = TBeam(Lp,Hw,Wf);

    for (size_t p=0; p!=beam.nPatches(); p++)
        result.addPatch(beam.patch(p));

    for (size_t p = 0; p!=result.nPatches(); p++)
    {
        result.patch(p).coefs().col(0).array() += x;
        result.patch(p).coefs().col(1).array() += y;
        result.patch(p).coefs().col(2).array() += z;
    }

    // result.addInterface(&result.patch(1),2,&result.patch(0),1);
    // result.addInterface(&result.patch(2),3,&result.patch(0),1);
    // result.addInterface(&result.patch(2),4,&result.patch(3),1);
    // result.addInterface(&result.patch(2),4,&result.patch(4),2);

    result.computeTopology();
    result.addAutoBoundaries();

    return result;
}

template <class T>
gsMultiPatch<T> gsPanelCreator<T>::PanelStrip(T const & Lp, T const & Wp, T const & Hw, T const & x, T const & y, T const & z)
{
    gsMultiPatch<T> result, tmp;
    std::vector<gsMultiPatch<T>> panels(3);
    panels.at(0) = Plate(Lp,Wp/2,0,0,0);
    panels.at(1) = Plate(Lp,Wp/2,0,-Wp/2.,0);
    panels.at(2) = Strip(Lp,Hw);

    for (typename std::vector<gsMultiPatch<T>>::iterator it = panels.begin(); it!=panels.end(); it++)
        for (size_t p = 0; p!=it->nPatches(); p++)
            result.addPatch(it->patch(p));

    for (size_t p = 0; p!=result.nPatches(); p++)
    {
        result.patch(p).coefs().col(0).array() += x;
        result.patch(p).coefs().col(1).array() += y;
        result.patch(p).coefs().col(2).array() += z;
    }

    result.computeTopology();
    result.addAutoBoundaries();
    return result;
}

template <class T>
gsMultiPatch<T> gsPanelCreator<T>::PanelL(T const & Lp, T const & Wp, T const & Hw, T const & Wf, T const & x, T const & y, T const & z)
{
    gsMultiPatch<T> result, tmp;
    std::vector<gsMultiPatch<T>> panels(4);
    panels.at(0) = Plate(Lp,Wp/2,0,-Wp/2.,0);
    panels.at(1) = Plate(Lp,Wf,0,0,0);
    panels.at(2) = Plate(Lp,Wp/2-Wf,0,Wf,0);

    // L-stiffener
    panels.at(3) = LBeam(Lp,Hw,Wf);
    for (typename std::vector<gsMultiPatch<T>>::iterator it = panels.begin(); it!=panels.end(); it++)
        for (size_t p = 0; p!=it->nPatches(); p++)
            result.addPatch(it->patch(p));

    for (size_t p = 0; p!=result.nPatches(); p++)
    {
        result.patch(p).coefs().col(0).array() += x;
        result.patch(p).coefs().col(1).array() += y;
        result.patch(p).coefs().col(2).array() += z;
    }

    result.computeTopology();
    result.addAutoBoundaries();
    return result;
}

template <class T>
gsMultiPatch<T> gsPanelCreator<T>::PlateGirderL(T const & Lp, T const & Wp, T const & Hwg, T const & Wfg, T const & Hws, T const & Wfs, T const & x, T const & y, T const & z)
{
    gsMultiPatch<T> result, tmp;

    // make sub panels
    std::vector<gsMultiPatch<T>> panels(14);
    panels.at(0) = Plate(Lp/2.,Wp/2.,                  0.,                 -Wp/2., 0.);
    panels.at(1) = Plate(Lp/2.,Wfs,                 0.,                 0.,             0.);
    panels.at(2) = Plate(Lp/2.,Wp/2.-Wfs,   0.,                 Wfs, 0.);
    panels.at(3) = Plate(Lp/2.,Wp/2.,                  -Lp/2.,    -Wp/2., 0.);
    panels.at(4) = Plate(Lp/2.,Wfs,                 -Lp/2.,    0.,             0.);
    panels.at(5) = Plate(Lp/2.,Wp/2.-Wfs,   -Lp/2.,    Wfs, 0.);

    panels.at(6) = LBeam(Lp/2.,Hws,Wfs,0);
    panels.at(7) = LBeam(Lp/2.,Hws,Wfs,-Lp/2.);

    panels.at(8) = TBeam(Wfs,              Hwg-Hws,Wfg,0,0,Hws);
    panels.at(9) = TBeam(Wp/2.-Wfs,Hwg-Hws,Wfg,0,0,Hws);
    panels.at(10) = TBeam(Wp/2.,              Hwg-Hws,Wfg,0,0,Hws);

    for (size_t p = 0; p!=panels[8].nPatches(); p++)
    {
        panels[8].patch(p).coefs().col(0).swap(panels[8].patch(p).coefs().col(1));
        // panels[8].patch(p).coefs().col(1).array() -= Wp / 2.;
    }
    for (size_t p = 0; p!=panels[9].nPatches(); p++)
    {
        panels[9].patch(p).coefs().col(0).swap(panels[9].patch(p).coefs().col(1));
        panels[9].patch(p).coefs().col(1).array() += Wfs;
    }
    for (size_t p = 0; p!=panels[10].nPatches(); p++)
    {
        panels[10].patch(p).coefs().col(0).swap(panels[10].patch(p).coefs().col(1));
        panels[10].patch(p).coefs().col(1).array() -= Wp / 2.;
    }


    panels.at(11) = Strip(Wfs,              Hws);
    panels.at(12) = Strip(Wp/2.-Wfs,Hws);
    panels.at(13) = Strip(Wp/2.,               Hws);
    for (size_t p = 0; p!=panels[11].nPatches(); p++)
    {
        panels[11].patch(p).coefs().col(0).swap(panels[11].patch(p).coefs().col(1));
        // panels[11].patch(p).coefs().col(1).array() -= Wp / 2.;
    }
    for (size_t p = 0; p!=panels[12].nPatches(); p++)
    {
        panels[12].patch(p).coefs().col(0).swap(panels[12].patch(p).coefs().col(1));
        panels[12].patch(p).coefs().col(1).array() += Wfs;
    }
    for (size_t p = 0; p!=panels[13].nPatches(); p++)
    {
        panels[13].patch(p).coefs().col(0).swap(panels[13].patch(p).coefs().col(1));
        panels[13].patch(p).coefs().col(1).array() -= Wp / 2.;
    }


    for (typename std::vector<gsMultiPatch<T>>::iterator it = panels.begin(); it!=panels.end(); it++)
        for (size_t p = 0; p!=it->nPatches(); p++)
            result.addPatch(it->patch(p));


    result.computeTopology();
    result.addAutoBoundaries();

    return result;
}

} // namespace gismo
