/** @file example_shell3D.cpp

    @brief Simple 3D examples for the shell class

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst (2019 - ..., TU Delft)
*/

#include <gismo.h>

#include <gsKLShell/src/gsThinShellAssembler.h>
#include <gsKLShell/src/gsMaterialMatrixContainer.h>
#include <gsKLShell/src/gsMaterialMatrixEval.h>
#include <gsKLShell/src/gsMaterialMatrixIntegrate.h>
#include <gsKLShell/src/getMaterialMatrix.h>
#include <gsCore/gsPiecewiseFunction.h>
#include <gsStructuralAnalysis/src/gsStructuralAnalysisTools/gsPanelCreator.h>

using namespace gismo;

template <class T>
gsMultiPatch<T> Strip(T Lp, T Wp, T x = 0, T y = 0, T z = 0);

template <class T>
gsMultiPatch<T> Plate(T Lp, T Wp, T alpha=0.5, T x = 0, T y = 0, T z = 0);

template <class T>
gsMultiPatch<T> Panel3(T Lp, T Wp, T Hw, T alpha=0.5, T x = 0, T y = 0, T z = 0);

template <class T>
gsMultiPatch<T> Panel4(T Lp, T Wp, T Hw, T alpha=0.5, T x = 0, T y = 0, T z = 0);


// Choose among various shell examples, default = Thin Plate
int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot  = false;
    bool stress= false;
    index_t numRefine  = 0;
    index_t numElevate = 1;
    index_t testCase = 1;
    bool nonlinear = false;

    bool membrane = false;
    bool composite = false;

    real_t E_modulus = 1.0;
    real_t PoissonRatio = 0.0;
    real_t Density = 1.0;
    real_t thickness = 1.0;

    real_t ifcDirichlet = 1.0;
    real_t ifcClamped = 1.0;

    real_t alpha = 0.5;
    real_t beta = 0.5;

    gsCmdLine cmd("2D shell example.");
    cmd.addReal( "D", "Dir", "Dirichlet penalty scalar",  ifcDirichlet );
    cmd.addReal( "C", "Cla", "Clamped penalty scalar",  ifcClamped );
    cmd.addReal( "a", "alpha", "Alpha parameter",  alpha );
    cmd.addReal( "b", "beta", "Beta parameter", beta );
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving",  numRefine );
    cmd.addInt( "t", "testCase", "Test case to run: 0 = square plate with pressure; 1 = Scordelis Lo Roof; 2 = quarter hemisphere; 3 = pinched cylinder",  testCase );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("stress", "Create a ParaView visualization file with the stresses", stress);
    cmd.addSwitch("membrane", "Use membrane model (no bending)", membrane);
    cmd.addSwitch("composite", "Composite material", composite);
    cmd.addSwitch( "nl", "Print information", nonlinear );

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Define material parameters and geometry per example]
    gsMultiPatch<> mp;
    gsMultiPatch<> mp_def;
    if (testCase == 0 )
    {
        thickness = 10.; // plate
        E_modulus = 206000;
        PoissonRatio = 0.3;
        mp = gsPanelCreator<real_t>::PanelT(3500.,3500.,200.,100.);
    }
    else if (testCase == 1 )
    {
        thickness = 10.; // plate
        E_modulus = 206000;
        PoissonRatio = 0.3;
        mp = gsPanelCreator<real_t>::PanelStrip(3500.,3500.,200.);
    }
    else if (testCase == 2 )
    {
        thickness = 10.; // plate
        E_modulus = 206000;
        PoissonRatio = 0.3;
        mp = gsPanelCreator<real_t>::PanelL(3500.,3500.,200.,38.);
    }
    else if (testCase == 3 )
    {
        thickness = 10.; // plate
        E_modulus = 206000;
        PoissonRatio = 0.3;
        mp = gsPanelCreator<real_t>::PlateGirderL(7020.,3500.,580.,250.,200.,38.);
    }
    else
        GISMO_ERROR("Testcase not found");
    //! [Define material parameters and geometry per example]

    //! [Refine and elevate]
    // p-refine
    if (numElevate!=0)
        mp.degreeElevate(numElevate);

    // h-refine
    for (int r =0; r < numRefine; ++r)
        mp.uniformRefine();

    mp_def = mp;
    if (plot) gsWriteParaview<>( mp_def    , "mp", 1000, true);
    //! [Refine and elevate]

    gsMultiBasis<> dbasis(mp);

    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
    gsInfo << dbasis.basis(0)<<"\n";


    //! [Set boundary conditions]
    gsBoundaryConditions<> bc;
    bc.setGeoMap(mp);

    gsPiecewiseFunction<> force(mp.nPatches());
    gsPiecewiseFunction<> t(mp.nPatches());
    // gsPiecewiseFunction<> nu(mp.nPatches());

    gsVector<> tmp(3), null(3);
    tmp << 0, 0, 1;
    null.setZero();

    gsMatrix<> refPoints;
    gsVector<index_t> refPatches;
    if (testCase == 0)
    {
        for (index_t d = 0; d!=3; d++)
        {
            bc.addCondition(0, boundary::east, condition_type::dirichlet, 0, 0, false, d);
            bc.addCondition(0, boundary::west, condition_type::dirichlet, 0, 0, false, d);
            bc.addCondition(0, boundary::north, condition_type::dirichlet, 0, 0, false, d);

            bc.addCondition(1, boundary::east, condition_type::dirichlet, 0, 0, false, d);
            bc.addCondition(1, boundary::west, condition_type::dirichlet, 0, 0, false, d);
            bc.addCondition(1, boundary::south, condition_type::dirichlet, 0, 0, false, d);

            bc.addCondition(2, boundary::east, condition_type::dirichlet, 0, 0, false, d);
            bc.addCondition(2, boundary::west, condition_type::dirichlet, 0, 0, false, d);

            bc.addCondition(3, boundary::east, condition_type::dirichlet, 0, 0, false, d);
            bc.addCondition(3, boundary::west, condition_type::dirichlet, 0, 0, false, d);

            bc.addCondition(4, boundary::east, condition_type::dirichlet, 0, 0, false, d);
            bc.addCondition(4, boundary::west, condition_type::dirichlet, 0, 0, false, d);

        }

        // Point loads
        tmp << 0,0,1e-2;
        gsConstantFunction<> piece0(tmp,3);
        force.addPiece(piece0);
        tmp << 0,0,1e-2;
        gsConstantFunction<> piece1(tmp,3);
        force.addPiece(piece1);
        tmp << 0,0,0;
        gsConstantFunction<> piece2(tmp,3);
        force.addPiece(piece2);
        tmp << 0,0,0;
        gsConstantFunction<> piece3(tmp,3);
        force.addPiece(piece3);
        tmp << 0,0,0;
        gsConstantFunction<> piece4(tmp,3);
        force.addPiece(piece4);

        // thickness
        gsFunctionExpr<> t0(std::to_string(thickness), 3);
        t.addPiece(t0);
        gsFunctionExpr<> t1(std::to_string(thickness), 3);
        t.addPiece(t1);
        gsFunctionExpr<> t2(std::to_string(thickness), 3);
        t.addPiece(t2);
        gsFunctionExpr<> t3(std::to_string(thickness), 3);
        t.addPiece(t3);
        gsFunctionExpr<> t4(std::to_string(thickness), 3);
        t.addPiece(t4);

        refPoints.resize(2,1);
        refPoints<<1.0,1.0;
        refPatches.resize(1);
        refPatches<<0;
    }
    else if (testCase == 1)
    {
        // Horizontal displacements
        bc.addCondition(2,boundary::east, condition_type::dirichlet, 0, 0, false, 0 );      // stiffener web
        bc.addCondition(2,boundary::west, condition_type::dirichlet, 0, 0, false, 0 );      // stiffener web
        // bc.addCondition(3,boundary::north, condition_type::dirichlet, 0, 0, false, 1 );     // stiffener flange
        // bc.addCondition(3,boundary::south, condition_type::dirichlet, 0, 0, false, 1 );     // stiffener flange

        // Vertical displacements
        bc.addCondition(0,boundary::east, condition_type::dirichlet, 0, 0, false, 2 );     // plate
        bc.addCondition(0,boundary::west, condition_type::dirichlet, 0, 0, false, 2 );     // plate
        bc.addCondition(0,boundary::north, condition_type::dirichlet, 0, 0, false, 2 );      // plate
        bc.addCondition(1,boundary::east, condition_type::dirichlet, 0, 0, false, 2 );     // plate
        bc.addCondition(1,boundary::west, condition_type::dirichlet, 0, 0, false, 2 );     // plate
        bc.addCondition(1,boundary::south, condition_type::dirichlet, 0, 0, false, 2 );      // plate
        bc.addCondition(2,boundary::east, condition_type::dirichlet, 0, 0, false, 2 );      // stiffener web
        bc.addCondition(2,boundary::west, condition_type::dirichlet, 0, 0, false, 2 );      // stiffener web
        // bc.addCondition(3,boundary::north, condition_type::dirichlet, 0, 0, false, 2 );     // stiffener flange
        // bc.addCondition(3,boundary::south, condition_type::dirichlet, 0, 0, false, 2 );     // stiffener flange


        // Clampings
        bc.addCondition(0,boundary::east, condition_type::clamped, 0, 0, false, 2 );
        bc.addCondition(0,boundary::west, condition_type::clamped, 0, 0, false, 2 );
        bc.addCondition(0,boundary::north, condition_type::clamped, 0, 0, false, 2 );
        bc.addCondition(1,boundary::east, condition_type::clamped, 0, 0, false, 2 );
        bc.addCondition(1,boundary::west, condition_type::clamped, 0, 0, false, 2 );
        bc.addCondition(1,boundary::south, condition_type::clamped, 0, 0, false, 2 );
        bc.addCondition(2,boundary::east, condition_type::clamped, 0, 0, false, 1 );
        bc.addCondition(2,boundary::west, condition_type::clamped, 0, 0, false, 1 );
        // bc.addCondition(3,boundary::north, condition_type::clamped, 0, 0, false, 2 );
        // bc.addCondition(3,boundary::south, condition_type::clamped, 0, 0, false, 2 );

        std::vector<index_t> plateIDs{0,1};
        std::vector<index_t> stiffenerWebIDs{2};
        // std::vector<index_t> stiffenerFlangeIDs{3};

        real_t pressure = 1e-3;
        // Surface loads
        gsVector<> pressvec(3);
        pressvec << 0,0,pressure;
        tmp << 0,0,0;
        gsConstantFunction<> forcep(pressvec,3);
        gsConstantFunction<> force0(tmp,3);

        typedef typename std::vector<gsFunction<real_t>*>       FunctionContainer;
        FunctionContainer forceContainer(mp.nPatches());
        // Plates
        for (typename std::vector<index_t>::iterator it = plateIDs.begin(); it!=plateIDs.end(); it++)
            forceContainer.at(*it) = &forcep;
        // Stiffener webs
        for (typename std::vector<index_t>::iterator it = stiffenerWebIDs.begin(); it!=stiffenerWebIDs.end(); it++)
            forceContainer.at(*it) = &force0;
        // Stiffener flanges
        // for (typename std::vector<index_t>::iterator it = stiffenerFlangeIDs.begin(); it!=stiffenerFlangeIDs.end(); it++)
        //     forceContainer.at(*it) = &force0;

        for (size_t p = 0; p != forceContainer.size(); p++)
            force.addPiece(*forceContainer.at(p));


        FunctionContainer tContainer(mp.nPatches());

        gsVector<> tmp2(1);
        gsFunctionExpr<> tp(std::to_string(1.), 3);
        gsFunctionExpr<> tsw(std::to_string(1.), 3);
        gsFunctionExpr<> tsf(std::to_string(1.), 3);

        // Plates
        for (typename std::vector<index_t>::iterator it = plateIDs.begin(); it!=plateIDs.end(); it++)
            tContainer.at(*it) = &tp;
        // Stiffener webs
        for (typename std::vector<index_t>::iterator it = stiffenerWebIDs.begin(); it!=stiffenerWebIDs.end(); it++)
            tContainer.at(*it) = &tsw;
        // Stiffener flanges
        // for (typename std::vector<index_t>::iterator it = stiffenerFlangeIDs.begin(); it!=stiffenerFlangeIDs.end(); it++)
        //     tContainer.at(*it) = &tsf;

        for (size_t p = 0; p != tContainer.size(); p++)
            t.addPiece(*tContainer.at(p));

        refPoints.resize(2,1);
        refPoints<<0,0.5;
        refPatches.resize(1);
        refPatches<<0;
    }
    else if (testCase == 2)
    {
        // Horizontal displacements
        bc.addCondition(3,boundary::east, condition_type::dirichlet, 0, 0, false, 1 );      // stiffener web
        bc.addCondition(3,boundary::west, condition_type::dirichlet, 0, 0, false, 1 );      // stiffener web
        bc.addCondition(4,boundary::east, condition_type::dirichlet, 0, 0, false, 1 );     // stiffener flange
        bc.addCondition(4,boundary::west, condition_type::dirichlet, 0, 0, false, 1 );     // stiffener flange

        // Vertical displacements
        bc.addCondition(0,boundary::east, condition_type::dirichlet, 0, 0, false, 2 );     // plate
        bc.addCondition(0,boundary::west, condition_type::dirichlet, 0, 0, false, 2 );     // plate
        bc.addCondition(0,boundary::south, condition_type::dirichlet, 0, 0, false, 2 );      // plate
        bc.addCondition(1,boundary::east, condition_type::dirichlet, 0, 0, false, 2 );     // plate
        bc.addCondition(1,boundary::west, condition_type::dirichlet, 0, 0, false, 2 );     // plate
        bc.addCondition(2,boundary::east, condition_type::dirichlet, 0, 0, false, 2 );     // plate
        bc.addCondition(2,boundary::west, condition_type::dirichlet, 0, 0, false, 2 );     // plate
        bc.addCondition(2,boundary::north, condition_type::dirichlet, 0, 0, false, 2 );      // plate
        bc.addCondition(3,boundary::east, condition_type::dirichlet, 0, 0, false, 2 );      // stiffener web
        bc.addCondition(3,boundary::west, condition_type::dirichlet, 0, 0, false, 2 );      // stiffener web
        bc.addCondition(4,boundary::east, condition_type::dirichlet, 0, 0, false, 2 );     // stiffener flange
        bc.addCondition(4,boundary::west, condition_type::dirichlet, 0, 0, false, 2 );     // stiffener flange

        // Clampings
        bc.addCondition(0,boundary::east, condition_type::clamped, 0, 0, false, 2 );
        bc.addCondition(0,boundary::west, condition_type::clamped, 0, 0, false, 2 );
        bc.addCondition(0,boundary::south, condition_type::clamped, 0, 0, false, 2 );
        bc.addCondition(1,boundary::east, condition_type::dirichlet, 0, 0, false, 2 );     // plate
        bc.addCondition(1,boundary::west, condition_type::dirichlet, 0, 0, false, 2 );     // plate
        bc.addCondition(2,boundary::east, condition_type::clamped, 0, 0, false, 2 );
        bc.addCondition(2,boundary::west, condition_type::clamped, 0, 0, false, 2 );
        bc.addCondition(2,boundary::north, condition_type::clamped, 0, 0, false, 2 );
        bc.addCondition(3,boundary::east, condition_type::clamped, 0, 0, false, 0 );
        bc.addCondition(3,boundary::west, condition_type::clamped, 0, 0, false, 0 );
        bc.addCondition(4,boundary::east, condition_type::clamped, 0, 0, false, 2 );
        bc.addCondition(4,boundary::west, condition_type::clamped, 0, 0, false, 2 );

        std::vector<index_t> plateIDs{0,1,2};
        std::vector<index_t> stiffenerWebIDs{3};
        std::vector<index_t> stiffenerFlangeIDs{4};

        real_t pressure = 1e-3;
        // Surface loads
        gsVector<> pressvec(3);
        pressvec << 0,0,pressure;
        tmp << 0,0,0;
        gsConstantFunction<> forcep(pressvec,3);
        gsConstantFunction<> force0(tmp,3);

        typedef typename std::vector<gsFunction<real_t>*>       FunctionContainer;
        FunctionContainer forceContainer(mp.nPatches());
        // Plates
        for (typename std::vector<index_t>::iterator it = plateIDs.begin(); it!=plateIDs.end(); it++)
            forceContainer.at(*it) = &forcep;
        // Stiffener webs
        for (typename std::vector<index_t>::iterator it = stiffenerWebIDs.begin(); it!=stiffenerWebIDs.end(); it++)
            forceContainer.at(*it) = &force0;
        // Stiffener flanges
        for (typename std::vector<index_t>::iterator it = stiffenerFlangeIDs.begin(); it!=stiffenerFlangeIDs.end(); it++)
            forceContainer.at(*it) = &force0;

        for (size_t p = 0; p != forceContainer.size(); p++)
            force.addPiece(*forceContainer.at(p));


        FunctionContainer tContainer(mp.nPatches());

        gsVector<> tmp2(1);
        gsFunctionExpr<> tp(std::to_string(1.), 3);
        gsFunctionExpr<> tsw(std::to_string(1.), 3);
        gsFunctionExpr<> tsf(std::to_string(1.), 3);

        gsDebugVar(tsf.targetDim());

        // Plates
        for (typename std::vector<index_t>::iterator it = plateIDs.begin(); it!=plateIDs.end(); it++)
            tContainer.at(*it) = &tp;
        // Stiffener webs
        for (typename std::vector<index_t>::iterator it = stiffenerWebIDs.begin(); it!=stiffenerWebIDs.end(); it++)
            tContainer.at(*it) = &tsw;
        // Stiffener flanges
        for (typename std::vector<index_t>::iterator it = stiffenerFlangeIDs.begin(); it!=stiffenerFlangeIDs.end(); it++)
            tContainer.at(*it) = &tsf;

        for (size_t p = 0; p != tContainer.size(); p++)
            t.addPiece(*tContainer.at(p));

        refPoints.resize(2,1);
        refPoints<<0.0,0.5;
        refPatches.resize(1);
        refPatches<<0;
    }
    else if (testCase == 3)
    {
        // Plate
        bc.addCondition(0,boundary::east, condition_type::dirichlet, 0, 0, false, 2 );
        bc.addCondition(0,boundary::south, condition_type::dirichlet, 0, 0, false, 2 );
        bc.addCondition(1,boundary::east, condition_type::dirichlet, 0, 0, false, 2 );
        bc.addCondition(2,boundary::east, condition_type::dirichlet, 0, 0, false, 2 );
        bc.addCondition(2,boundary::north, condition_type::dirichlet, 0, 0, false, 2 );
        bc.addCondition(3,boundary::west, condition_type::dirichlet, 0, 0, false, 2 );
        bc.addCondition(3,boundary::south, condition_type::dirichlet, 0, 0, false, 2 );
        bc.addCondition(4,boundary::west, condition_type::dirichlet, 0, 0, false, 2 );
        bc.addCondition(5,boundary::west, condition_type::dirichlet, 0, 0, false, 2 );
        bc.addCondition(5,boundary::north, condition_type::dirichlet, 0, 0, false, 2 );

        bc.addCondition(0,boundary::east, condition_type::clamped, 0, 0, false, 2 );
        bc.addCondition(0,boundary::south, condition_type::clamped, 0, 0, false, 2 );
        bc.addCondition(1,boundary::east, condition_type::clamped, 0, 0, false, 2 );
        bc.addCondition(2,boundary::east, condition_type::clamped, 0, 0, false, 2 );
        bc.addCondition(2,boundary::north, condition_type::clamped, 0, 0, false, 2 );
        bc.addCondition(3,boundary::west, condition_type::clamped, 0, 0, false, 2 );
        bc.addCondition(3,boundary::south, condition_type::clamped, 0, 0, false, 2 );
        bc.addCondition(4,boundary::west, condition_type::clamped, 0, 0, false, 2 );
        bc.addCondition(5,boundary::west, condition_type::clamped, 0, 0, false, 2 );
        bc.addCondition(5,boundary::north, condition_type::clamped, 0, 0, false, 2 );

        // Girder
        // Web
        bc.addCondition(13,boundary::east, condition_type::dirichlet, 0, 0, false, 2 );
        bc.addCondition(20,boundary::east, condition_type::dirichlet, 0, 0, false, 2 );
        bc.addCondition(16,boundary::west, condition_type::dirichlet, 0, 0, false, 2 );
        bc.addCondition(21,boundary::west, condition_type::dirichlet, 0, 0, false, 2 );

        bc.addCondition(13,boundary::east, condition_type::clamped, 0, 0, false, 1 );
        bc.addCondition(20,boundary::east, condition_type::clamped, 0, 0, false, 1 );
        bc.addCondition(21,boundary::west, condition_type::clamped, 0, 0, false, 1 );
        bc.addCondition(16,boundary::west, condition_type::clamped, 0, 0, false, 1 );

        //Flange
        bc.addCondition(14,boundary::north, condition_type::dirichlet, 0, 0, false, 1 );
        bc.addCondition(15,boundary::north, condition_type::dirichlet, 0, 0, false, 1 );
        bc.addCondition(17,boundary::south, condition_type::dirichlet, 0, 0, false, 1 );
        bc.addCondition(18,boundary::south, condition_type::dirichlet, 0, 0, false, 1 );

        bc.addCondition(14,boundary::north, condition_type::dirichlet, 0, 0, false, 2 );
        bc.addCondition(15,boundary::north, condition_type::dirichlet, 0, 0, false, 2 );
        bc.addCondition(17,boundary::south, condition_type::dirichlet, 0, 0, false, 2 );
        bc.addCondition(18,boundary::south, condition_type::dirichlet, 0, 0, false, 2 );

        bc.addCondition(14,boundary::north, condition_type::clamped, 0, 0, false, 2 );
        bc.addCondition(15,boundary::north, condition_type::clamped, 0, 0, false, 2 );
        bc.addCondition(17,boundary::south, condition_type::clamped, 0, 0, false, 2 );
        bc.addCondition(18,boundary::south, condition_type::clamped, 0, 0, false, 2 );

        // Stiffener
        // Web
        bc.addCondition(8,boundary::west, condition_type::dirichlet, 0, 0, false, 0 );
        bc.addCondition(6,boundary::east, condition_type::dirichlet, 0, 0, false, 0 );

        bc.addCondition(8,boundary::west, condition_type::dirichlet, 0, 0, false, 2 );
        bc.addCondition(6,boundary::east, condition_type::dirichlet, 0, 0, false, 2 );

        bc.addCondition(8,boundary::west, condition_type::clamped, 0, 0, false, 1 );
        bc.addCondition(6,boundary::east, condition_type::clamped, 0, 0, false, 1 );
        //Flange
        bc.addCondition(9,boundary::south, condition_type::dirichlet, 0, 0, false, 0 );
        bc.addCondition(7,boundary::north, condition_type::dirichlet, 0, 0, false, 0 );

        bc.addCondition(9,boundary::south, condition_type::dirichlet, 0, 0, false, 2 );
        bc.addCondition(7,boundary::north, condition_type::dirichlet, 0, 0, false, 2 );

        bc.addCondition(9,boundary::south, condition_type::clamped, 0, 0, false, 2 );
        bc.addCondition(7,boundary::north, condition_type::clamped, 0, 0, false, 2 );

        std::vector<index_t> plateIDs{0,1,2,3,4,5};
        std::vector<index_t> girderWebIDs{10,13,16,19,20,21};
        std::vector<index_t> girderFlangeIDs{11,12,14,15,17,18};
        std::vector<index_t> stiffenerWebIDs{6,8};
        std::vector<index_t> stiffenerFlangeIDs{7,9};

        real_t pressure = 1e-3;
        // Surface loads
        gsVector<> pressvec(3);
        pressvec << 0,0,pressure;
        tmp << 0,0,0;
        gsConstantFunction<> forcep(pressvec,3);
        gsConstantFunction<> force0(tmp,3);

        typedef typename std::vector<gsFunction<real_t>*>       FunctionContainer;
        FunctionContainer forceContainer(mp.nPatches());
        // Plates
        for (typename std::vector<index_t>::iterator it = plateIDs.begin(); it!=plateIDs.end(); it++)
            forceContainer.at(*it) = &forcep;
        // Girder webs
        for (typename std::vector<index_t>::iterator it = girderWebIDs.begin(); it!=girderWebIDs.end(); it++)
            forceContainer.at(*it) = &force0;
        // Girder flanges
        for (typename std::vector<index_t>::iterator it = girderFlangeIDs.begin(); it!=girderFlangeIDs.end(); it++)
            forceContainer.at(*it) = &force0;
        // Stiffener webs
        for (typename std::vector<index_t>::iterator it = stiffenerWebIDs.begin(); it!=stiffenerWebIDs.end(); it++)
            forceContainer.at(*it) = &force0;
        // Stiffener flanges
        for (typename std::vector<index_t>::iterator it = stiffenerFlangeIDs.begin(); it!=stiffenerFlangeIDs.end(); it++)
            forceContainer.at(*it) = &force0;

        for (size_t p = 0; p != forceContainer.size(); p++)
            force.addPiece(*forceContainer.at(p));

        // force = gsPiecewiseFunction<>(forceContainer);

        // Plate (patches 0, 1, 4, 5)
        // tp = 10
        // Girder web (patches 0, 1, 4, 5)
        // tgw = 15
        // Girder flange (patches 0, 1, 4, 5)
        // tgf = 30
        // Stiffener web (patches 0, 1, 4, 5)
        // tsw = 10
        // Stiffener flange (patches 0, 1, 4, 5)
        // tsf = 15

        FunctionContainer tContainer(mp.nPatches());

        gsVector<> tmp2(1);
        gsFunctionExpr<> tp(std::to_string(1.), 3);
        gsFunctionExpr<> tgw(std::to_string(1.), 3);
        gsFunctionExpr<> tgf(std::to_string(1.), 3);
        gsFunctionExpr<> tsw(std::to_string(1.), 3);
        gsFunctionExpr<> tsf(std::to_string(1.), 3);

        gsDebugVar(tsf.targetDim());

        // Plates
        for (typename std::vector<index_t>::iterator it = plateIDs.begin(); it!=plateIDs.end(); it++)
            tContainer.at(*it) = &tp;
        // Girder webs
        for (typename std::vector<index_t>::iterator it = girderWebIDs.begin(); it!=girderWebIDs.end(); it++)
            tContainer.at(*it) = &tgw;
        // Girder flanges
        for (typename std::vector<index_t>::iterator it = girderFlangeIDs.begin(); it!=girderFlangeIDs.end(); it++)
            tContainer.at(*it) = &tgf;
        // Stiffener webs
        for (typename std::vector<index_t>::iterator it = stiffenerWebIDs.begin(); it!=stiffenerWebIDs.end(); it++)
            tContainer.at(*it) = &tsw;
        // Stiffener flanges
        for (typename std::vector<index_t>::iterator it = stiffenerFlangeIDs.begin(); it!=stiffenerFlangeIDs.end(); it++)
            tContainer.at(*it) = &tsf;

        for (size_t p = 0; p != tContainer.size(); p++)
            t.addPiece(*tContainer.at(p));

        // t = gsPiecewiseFunction<>(tContainer);
    }
    else
        GISMO_ERROR("Test case not known");
    //! [Set boundary conditions]


    //! [Make material functions]
    // Linear isotropic material model
    gsFunctionExpr<> E(std::to_string(E_modulus),3);
    // gsConstantFunction<> t(thickness,3);
    gsFunctionExpr<> rho(std::to_string(Density),3);
    gsConstantFunction<> nu(PoissonRatio,3);
    // gsFunctionExpr<> force("0","0","0",3);

    // Linear anisotropic material model (only one layer for example purposes)
    index_t kmax = 1; // number of layers
    std::vector<gsFunctionSet<> * > Gs(kmax); // Material matrices
    std::vector<gsFunctionSet<> * > Ts(kmax); // Thickness per layer
    std::vector<gsFunctionSet<> * > Phis(kmax); // Fiber angle per layer

    // Make material matrix
    gsMatrix<> Gmat = gsCompositeMatrix(E_modulus,E_modulus,0.5 * E_modulus / (1+PoissonRatio),PoissonRatio,PoissonRatio);
    Gmat.resize(Gmat.rows()*Gmat.cols(),1);
    gsConstantFunction<> Gfun(Gmat,3);
    Gs[0] = &Gfun;

    // Define fiber angle
    gsConstantFunction<> phi;
    phi.setValue(0,3);
    Phis[0] = &phi;

    // Define thickness
    gsConstantFunction<> thicks(thickness/kmax,3);
    Ts[0] = &thicks;

    //! [Make assembler]
    std::vector<gsFunctionSet<>*> parameters;
    gsMaterialMatrixBase<real_t>* materialMatrix;
    gsOptionList options;
    // Make gsMaterialMatrix depending on the user-defined choices
    if (composite)
    {
        materialMatrix = new gsMaterialMatrixComposite<3,real_t>(mp,Ts,Gs,Phis);
    }
    else
    {
        parameters.resize(2);
        parameters[0] = &E;
        parameters[1] = &nu;
        options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
        options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
        materialMatrix = getMaterialMatrix<3,real_t>(mp,t,parameters,rho,options);
    }

    gsMaterialMatrixContainer<real_t> materialMats(mp.nPatches());
    for (size_t p = 0; p!=mp.nPatches(); p++)
        materialMats.add(materialMatrix);


    // gsMaterialMatrixContainer<real_t> materialMatsSingle(mp.nPatches());
    // gsMaterialMatrixBase<real_t> * mmtmp;
    // for (size_t p = 0; p!=mp.nPatches(); p++)
    // {
    //     parameters.resize(2);
    //     parameters[0] = const_cast<gsFunction<> *>(&(E.function(p)));
    //     parameters[1] = const_cast<gsFunction<> *>(&(nu.function(p)));
    //     options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
    //     options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
    //     mmtmp = getMaterialMatrix<3,real_t>(mp,t,parameters,rho,options);

    //     materialMatsSingle.add(mmtmp);
    // }



    // gsDebugVar(materialMats);
    // gsDebugVar(materialMatsSingle);

    // gsMatrix<> z(1,1);
    // z.setZero();
    // // gsMaterialMatrixEval2<real_t,MaterialOutput::VectorN> vectorN(materialMats,&mp,z);

    // gsMaterialMatrixEval<real_t,MaterialOutput::Thickness> Thickness(materialMatrix,&mp,z);
    // gsMaterialMatrixEval<real_t,MaterialOutput::Parameters> Parameters(materialMatrix,&mp,z);


    // gsMatrix<> result;
    // gsVector<> pt(2); pt.setConstant(0.25);

    // for (size_t p = 0; p!=mp.nPatches(); ++p)
    // {
    //     gsDebug<<"-----------Patch "<<p<<"\n";
    //     Thickness.piece(p).eval_into(pt,result);
    //     gsDebug<<"Thickness: "<<result.transpose()<<"\n";

    //     Parameters.piece(p).eval_into(pt,result);
    //     gsDebug<<"Parameters: "<<result.transpose()<<"\n";
    // }


    // gsMaterialMatrixEval<real_t,MaterialOutput::Thickness> Thickness2(materialMatsSingle,&mp,z);
    // gsMaterialMatrixEval<real_t,MaterialOutput::Parameters> Parameters2(materialMatsSingle,&mp,z);


    // for (size_t p = 0; p!=mp.nPatches(); ++p)
    // {
    //     gsDebug<<"-----------Patch "<<p<<"\n";
    //     Thickness2.piece(p).eval_into(pt,result);
    //     gsDebug<<"Thickness: "<<result.transpose()<<"\n";

    //     Parameters2.piece(p).eval_into(pt,result);
    //     gsDebug<<"Parameters: "<<result.transpose()<<"\n";
    // }

    // Construct the gsThinShellAssembler
    gsThinShellAssemblerBase<real_t>* assembler;
    if(membrane) // no bending term
        assembler = new gsThinShellAssembler<3, real_t, false>(mp,dbasis,bc,force,materialMatrix);
    else
        assembler = new gsThinShellAssembler<3, real_t, true >(mp,dbasis,bc,force,materialMatrix);

    // Set the penalty parameter for the interface C1 continuity
    assembler->options().setInt("Continuity",-1);
    assembler->options().setReal("IfcDirichlet",ifcDirichlet);
    assembler->options().setReal("IfcClamped",ifcClamped);
    assembler->addWeakC0(mp.topology().interfaces());
    assembler->addWeakC1(mp.topology().interfaces());
    assembler->initInterfaces();
    //! [Make assembler]

    // Set stopwatch
    gsStopwatch stopwatch,stopwatch2;
    real_t time = 0.0;
    real_t totaltime = 0.0;

    //! [Define jacobian and residual]
    // Function for the Jacobian
    typedef std::function<gsSparseMatrix<real_t> (gsVector<real_t> const &)>    Jacobian_t;
    typedef std::function<gsVector<real_t> (gsVector<real_t> const &) >         Residual_t;
    Jacobian_t Jacobian = [&time,&stopwatch,&assembler,&mp_def](gsVector<real_t> const &x)
    {
      stopwatch.restart();
      assembler->constructSolution(x,mp_def);
      assembler->assembleMatrix(mp_def);
      time += stopwatch.stop();
      gsSparseMatrix<real_t> m = assembler->matrix();
      return m;
    };
    // Function for the Residual
    Residual_t Residual = [&time,&stopwatch,&assembler,&mp_def](gsVector<real_t> const &x)
    {
      stopwatch.restart();
      assembler->constructSolution(x,mp_def);
      assembler->assembleVector(mp_def);
      time += stopwatch.stop();
      return assembler->rhs();
    };
    //! [Define jacobian and residual]

    stopwatch.restart();
    stopwatch2.restart();
    assembler->assemble();
    time += stopwatch.stop();

    //! [Assemble linear part]
    gsSparseMatrix<> matrix = assembler->matrix();
    gsVector<> vector = assembler->rhs();
    //! [Assemble linear part]

    //! [Solve linear problem]
    gsInfo<<"Solving system with "<<assembler->numDofs()<<" DoFs\n";
    gsVector<> solVector;
    gsSparseSolver<>::CGDiagonal solver;
    solver.compute( matrix );
    solVector = solver.solve(vector);
    //! [Solve linear problem]


    //! [Solve non-linear problem]
    if (nonlinear)
    {
        real_t residual = vector.norm();
        real_t residual0 = residual;
        real_t residualOld = residual;
        gsVector<real_t> updateVector = solVector;
        gsVector<real_t> resVec = Residual(solVector);
        gsSparseMatrix<real_t> jacMat;
        for (index_t it = 0; it != 100; ++it)
        {
            jacMat = Jacobian(solVector);
            solver.compute(jacMat);
            updateVector = solver.solve(resVec); // this is the UPDATE
            solVector += updateVector;

            resVec = Residual(solVector);
            residual = resVec.norm();

            gsInfo<<"Iteration: "<< it
               <<", residue: "<< residual
               <<", update norm: "<<updateVector.norm()
               <<", log(Ri/R0): "<< math::log10(residualOld/residual0)
               <<", log(Ri+1/R0): "<< math::log10(residual/residual0)
               <<"\n";

            residualOld = residual;

            if (updateVector.norm() < 1e-6)
                break;
            else if (it+1 == it)
                gsWarn<<"Maximum iterations reached!\n";
        }
    }
    //! [Solve non-linear problem]

    totaltime += stopwatch2.stop();

    //! [Construct and evaluate solution]
    mp_def = assembler->constructSolution(solVector);
    gsMultiPatch<> deformation = assembler->constructDisplacement(solVector);
    //! [Construct and evaluate solution]


    //! [Construct and evaluate solution]
    gsMatrix<> refVals(3,refPoints.cols());
    for (index_t k=0; k!=refPoints.cols(); k++)
        refVals.col(k) = deformation.patch(refPatches[k]).eval(refPoints.col(k));

    // gsInfo << "Displacement at reference point: "<<numVal<<"\n";
    gsInfo << "Displacement at reference point: "<<refVals<<"\n";
    //! [Construct and evaluate solution]

    // ! [Export visualization in ParaView]
    if (plot)
    {
        gsField<> solField(mp_def, deformation);
        // gsField<> solField(mp, deformation);
        gsInfo<<"Plotting in Paraview...\n";
        // gsWriteParaview<>( solField, "Deformation", 1000, true);
        gsWriteParaview<>( solField, "Deformation", 1000, false);
    }
    if (stress)
    {
        gsPiecewiseFunction<> membraneStresses;
        assembler->constructStress(mp_def,membraneStresses,stress_type::membrane);
        gsField<> membraneStress(mp_def,membraneStresses, true);

        gsPiecewiseFunction<> flexuralStresses;
        assembler->constructStress(mp_def,flexuralStresses,stress_type::flexural);
        gsField<> flexuralStress(mp_def,flexuralStresses, true);

        gsWriteParaview(membraneStress,"MembraneStress",1000);
        gsWriteParaview(flexuralStress,"FlexuralStress",1000);
    }
    // ! [Export visualization in ParaView]

    gsInfo<<"Total ellapsed assembly time: \t\t"<<time<<" s\n";
    gsInfo<<"Total ellapsed solution time (incl. assembly): \t"<<totaltime<<" s\n";

    delete assembler;
    delete materialMatrix;
    return EXIT_SUCCESS;

}// end main

template <class T>
gsMultiPatch<T> Plate(T Lp, T Wp, T alpha, T x, T y, T z)
{
    gsMultiPatch<T> result, tmp;

    // Base plate, left
    result.addPatch(gsNurbsCreator<>::BSplineSquare());
    result.patch(0).embed(3);
    result.patch(0).coefs().row(0)<< Wp*alpha,0,0;
    result.patch(0).coefs().row(1)<< Wp,0,0;
    result.patch(0).coefs().row(2)<< Wp*alpha,Lp,0;
    result.patch(0).coefs().row(3)<< Wp,Lp,0;

    // Base plate, right
    result.addPatch(gsNurbsCreator<>::BSplineSquare());
    result.patch(1).embed(3);
    result.patch(1).coefs().row(0)<< 0,0,0;
    result.patch(1).coefs().row(1)<< Wp*alpha,0,0;
    result.patch(1).coefs().row(2)<< 0,Lp,0;
    result.patch(1).coefs().row(3)<< Wp*alpha,Lp,0;

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
gsMultiPatch<T> Strip(T Lb, T Hw, T x, T y, T z)
{
    gsMultiPatch<T> result, tmp;

    // Web
    result.addPatch(gsNurbsCreator<>::BSplineSquare());
    result.patch(0).embed(3);
    result.patch(0).coefs().row(0)<< 0,0,0;
    result.patch(0).coefs().row(1)<< 0,Lb,0;
    result.patch(0).coefs().row(2)<< 0,0,Hw;
    result.patch(0).coefs().row(3)<< 0,Lb,Hw;

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
gsMultiPatch<T> Panel3(T Lp, T Wp, T Hw, T alpha, T x, T y, T z)

{
    gsMultiPatch<T> result, tmp;

    // Base plate, left
    result.addPatch(gsNurbsCreator<>::BSplineSquare());
    result.patch(0).embed(3);
    result.patch(0).coefs().row(0)<< Wp*alpha,0,0;
    result.patch(0).coefs().row(1)<< Wp,0,0;
    result.patch(0).coefs().row(2)<< Wp*alpha,Lp,0;
    result.patch(0).coefs().row(3)<< Wp,Lp,0;

    // Base plate, right
    result.addPatch(gsNurbsCreator<>::BSplineSquare());
    result.patch(1).embed(3);
    result.patch(1).coefs().row(0)<< 0,0,0;
    result.patch(1).coefs().row(1)<< Wp*alpha,0,0;
    result.patch(1).coefs().row(2)<< 0,Lp,0;
    result.patch(1).coefs().row(3)<< Wp*alpha,Lp,0;

    // T-Beam
    gsMultiPatch<> beam = Strip(Lp,Hw,Wp*alpha);

    for (size_t p=0; p!=beam.nPatches(); p++)
        result.addPatch(beam.patch(p));

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
gsMultiPatch<T> Panel4(T Lp, T Wp, T Hw, T alpha, T x, T y, T z)

{
    gsMultiPatch<T> result, tmp;

    // Base plate, left
    result.addPatch(gsNurbsCreator<>::BSplineSquare());
    result.patch(0).embed(3);
    result.patch(0).coefs().row(0)<< 0,0,0;
    result.patch(0).coefs().row(1)<< Wp/2,0,0;
    result.patch(0).coefs().row(2)<< 0,Lp,0;
    result.patch(0).coefs().row(3)<< Wp/2,Lp,0;

    // Base plate, middle
    result.addPatch(gsNurbsCreator<>::BSplineSquare());
    result.patch(1).embed(3);
    result.patch(1).coefs().row(0)<< Wp/2,0,0;
    result.patch(1).coefs().row(1)<< Wp*alpha,0,0;
    result.patch(1).coefs().row(2)<< Wp/2,Lp,0;
    result.patch(1).coefs().row(3)<< Wp*alpha,Lp,0;

    // Base plate, right
    result.addPatch(gsNurbsCreator<>::BSplineSquare());
    result.patch(2).embed(3);
    result.patch(2).coefs().row(0)<< Wp*alpha,0,0;
    result.patch(2).coefs().row(1)<< Wp,0,0;
    result.patch(2).coefs().row(2)<< Wp*alpha,Lp,0;
    result.patch(2).coefs().row(3)<< Wp,Lp,0;

    // T-Beam
    gsMultiPatch<> beam = Strip(Lp,Hw,Wp/2);

    for (size_t p=0; p!=beam.nPatches(); p++)
        result.addPatch(beam.patch(p));

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

