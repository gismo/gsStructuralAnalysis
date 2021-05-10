/** @file gsThinShell_Static.cpp

    @brief Static simulations of a shell

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#include <gismo.h>

#include <gsKLShell/gsThinShellAssembler.h>
#include <gsKLShell/getMaterialMatrix.h>

#include <gsElasticity/gsWriteParaviewMultiPhysics.h>

#include <gsStructuralAnalysis/gsStaticSolver.h>

//#include <gsThinShell/gsNewtonIterator.h>

using namespace gismo;

template <class T>
gsMultiPatch<T> RectangularDomain(int n, int m, int p, int q, T L, T B);
template <class T>
gsMultiPatch<T> RectangularDomain(int n, int p, T L, T B);

template <class T>
gsMultiPatch<T> Rectangle(T L, T B);

template <class T>
gsMultiPatch<T> RectangularDomainVert(int n, int m, int p, int q, T L, T B);
template <class T>
gsMultiPatch<T> RectangularDomainVert(int n, int p, T L, T B);

template <class T>
gsMultiPatch<T> RectangularDomain90(int n, int m, int p, int q, T L, T B);
template <class T>
gsMultiPatch<T> RectangularDomain90(int n, int p, T L, T B);

template <class T>
gsMultiPatch<T> FrustrumDomain(int n, int p, T R1, T R2, T h);

// Choose among various shell examples, default = Thin Plate
int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot  = false;
    bool stress= false;
    index_t numRefine  = 1;
    index_t numElevate = 1;
    index_t testCase = 1;

    index_t Compressibility = 0;
    index_t material = 0;
    bool composite = false;
    index_t impl = 1; // 1= analytical, 2= generalized, 3= spectral

    bool nonlinear = false;
    int verbose = 0;
    std::string fn;
    bool membrane = false;
    bool weak = false;

    std::string assemberOptionsFile("options/solver_options.xml");

    real_t E_modulus = 1.0;
    real_t PoissonRatio = 0.0;
    real_t Density = 1.0;
    real_t thickness = 1.0;
    real_t Ratio = 7.0;

    gsCmdLine cmd("Static analysis for thin shells.");
    cmd.addString( "f", "file", "Input XML file for assembler options", assemberOptionsFile );
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving",  numRefine );
    cmd.addReal( "R", "Ratio", "Mooney Rivlin Ratio",  Ratio );
    cmd.addInt( "t", "testCase", "Test case to run: 1 = unit square; 2 = Scordelis Lo Roof",  testCase );

    cmd.addInt( "m", "Material", "Material law",  material );
    cmd.addInt( "c", "Compressibility", "1: compressible, 0: incompressible",  Compressibility );
    cmd.addInt( "I", "Implementation", "Implementation: 1= analytical, 2= generalized, 3= spectral",  impl );
    cmd.addSwitch("composite", "Composite material", composite);

    cmd.addSwitch("nl", "Solve nonlinear problem", nonlinear);
    cmd.addInt("v","verbose", "0: no; 1: iteration output; 2: Full matrix and vector output", verbose);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("stress", "Create a ParaView visualization file with the stresses", stress);
    cmd.addSwitch("membrane", "Use membrane model (no bending)", membrane);
    cmd.addSwitch("weak", "Impose boundary conditions weakly", weak);


    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Read input file]
    gsFileData<> fd(assemberOptionsFile);
    gsOptionList opts;
    fd.getFirst<gsOptionList>(opts);

    gsMultiPatch<> mp;
    gsMultiPatch<> mp_def;
    if (testCase == 1 )
    {
        thickness = 0.25;
        E_modulus = 4.32E8;
        fn = "surface/scordelis_lo_roof.xml";
        gsReadFile<>(fn, mp);
        PoissonRatio = 0.0;
    }
    else if (testCase == 2)
    {
        thickness = 0.04;
        E_modulus = 6.825E7;
        PoissonRatio = 0.3;
        gsReadFile<>("surface/quarter_hemisphere.xml", mp);
    }
    else if (testCase == 3)
    {
        thickness = 3;
        E_modulus = 3E6;
        PoissonRatio = 0.3;
        gsReadFile<>("surface/pinched_cylinder.xml", mp);
    }
    else if (testCase == 4)
    {
        real_t L = 1.0;
        real_t B = 1.0;
        real_t mu = 1.5e6;
        thickness = 0.001;
        if (!Compressibility)
          PoissonRatio = 0.5;
        else
          PoissonRatio = 0.45;
        E_modulus = 2*mu*(1+PoissonRatio);
        // PoissonRatio = 0;
        mp = Rectangle(L, B);

        gsInfo<<"mu = "<<E_modulus / (2 * (1 + PoissonRatio))<<"\n";
    }
    else if (testCase == 5)
    {
        thickness = 1;
        E_modulus = 1;
        PoissonRatio = 0.5;
        gsReadFile<>("planar/unitcircle.xml", mp);
    }
    else if (testCase == 6)
    {
        thickness = 0.1;
        real_t mu = 4.225e5;
        PoissonRatio = 0.3;
        E_modulus = 2*mu*(1+PoissonRatio);
        gsReadFile<>("surface/quarter_sphere.xml", mp);
    }
    else if (testCase == 7)
    {
        thickness = 0.1;
        real_t mu = 4.225e5;
        PoissonRatio = 0.3;
        E_modulus = 2*mu*(1+PoissonRatio);
        gsReadFile<>("surface/quarter_frustrum.xml", mp);
    }
    // else if (testCase == 8)
    // {
    //     thickness = 1;
    //     PoissonRatio = 0.4998999999;
    //     E_modulus = 240.5653612;
    //     gsReadFile<>("surfaces/cooks_membrane.xml", mp);
    //     mp.embed(3);

    //     real_t mu = E_modulus/(2*(1+PoissonRatio));
    //     real_t K = E_modulus/(3-6*PoissonRatio);
    //     gsDebug<<"K = "<<K<<"\tmu = "<<mu<<"\n";
    // }
    else if (testCase == 8 || testCase == 9)
    {
        thickness = 0.1;
        real_t mu = 4.225;
        PoissonRatio = 0.5;
        E_modulus = 2*mu*(1+PoissonRatio);
        // gsReadFile<>("quarter_frustrum.xml", mp);

        // R1 is radius on bottom, R2 is radius on top
        mp = FrustrumDomain(numRefine,numElevate+2,2.0,1.0,1.0);
    }
    else if (testCase == 10)
    {
        E_modulus = 1;
        // thickness = 0.15;
        thickness = 1;
        if (!Compressibility)
          PoissonRatio = 0.499;
        else
          PoissonRatio = 0.45;

        E_modulus = 1;

        // real_t bDim = thickness / 1.9e-3;
        // real_t aDim = 2*bDim;

        real_t bDim = 1;
        real_t aDim = 1;


        // Ratio = 2.5442834138486314;
        // Ratio = 1e2;

        mp = Rectangle(aDim, bDim);

        // mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // degree
        // mp.addAutoBoundaries();
        // mp.embed(3);
        // E_modulus = 1e0;
        // thickness = 1e0;
        // // PoissonRatio = 0.0;
        // PoissonRatio = 0.4999;

    }
    else if (testCase == 17)
    {
        // Unit square
        mp = RectangularDomain(0,2,1.0,1.0);
        E_modulus = 4.497000000e6;
        thickness = 0.001;
        PoissonRatio = 0.4999;
    }
    else if (testCase == 18)
    {
        // Unit square
        // gsReadFile<>("planar/annulus_4p.xml", mp);
        gsReadFile<>("surface/frustrum.xml", mp);
        mp.computeTopology();

        E_modulus = 1;
        thickness = 1;
        PoissonRatio = 0;
    }
    else if (testCase == 20)
    {
        // Unit square
        // gsReadFile<>("planar/annulus_4p.xml", mp);
        gsMultiPatch<> mp_old;
        mp_old.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // degree
        mp = mp_old.uniformSplit();
        mp.computeTopology();
        mp.embed(3);
        E_modulus = 1;
        thickness = 1;
        PoissonRatio = 0;
    }
    else if (testCase==13)
    {
        // Unit square
        mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // degree
        mp.addAutoBoundaries();
        mp.embed(3);
        E_modulus = 1e6;
        thickness = 0.005;
        // PoissonRatio = 0.5;
        // PoissonRatio = 0.499;
        PoissonRatio = 0.3;
    }
    else
    {
        // Unit square
        mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // degree
        mp.addAutoBoundaries();
        mp.embed(3);
        E_modulus = 1e0;
        thickness = 1e0;
        // PoissonRatio = 0.5;
        // PoissonRatio = 0.499;
        PoissonRatio = 0.0;
        // PoissonRatio = 0.5;
    }
    //! [Read input file]
    // p-refine
    if (testCase != 8 && testCase != 9)
    {
        if (numElevate!=0)
            mp.degreeElevate(numElevate);

        // h-refine
        for (index_t r =0; r < numRefine; ++r)
            mp.uniformRefine();
    }
    mp_def = mp;
    gsWriteParaview<>( mp_def    , "mp", 1000, true);


    //! [Refinement]
    gsMultiBasis<> dbasis(mp);

    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
    gsInfo<<mp_def<<"\n";
    gsInfo << dbasis.basis(0)<<"\n";

    gsBoundaryConditions<> bc;
    bc.setGeoMap(mp);
    gsVector<> tmp(3);
    tmp << 0, 0, 0;

    gsVector<> neu(3);
    neu << 0, 0, 0;

    gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();

    gsConstantFunction<> displx(0.1,3);
    gsConstantFunction<> disply(0.25,3);

    gsFunctionExpr<> neuDataFun1;
    gsConstantFunction<> neuData(neu,3);
    gsVector<> weak_D(3);
    gsVector<> weak_C(3);
    weak_D.setZero();
    weak_C.setZero();
    gsConstantFunction<> weak_drch(weak_D,3);
    gsConstantFunction<> weak_clmp(weak_C,3);
    real_t pressure = 0.0;
    if (testCase == 0)
    {
        if (weak)
        {
            weak_drch.setValue(weak_D,3);
            bc.addCondition(boundary::north, condition_type::weak_dirichlet, &weak_drch, 0, false, -1 ); // unknown 2 - z
            bc.addCondition(boundary::east, condition_type::weak_dirichlet, &weak_drch, 0, false, -1 ); // unknown 2 - z
            bc.addCondition(boundary::south, condition_type::weak_dirichlet, &weak_drch, 0, false, -1); // unknown 2 - z
            bc.addCondition(boundary::west, condition_type::weak_dirichlet, &weak_drch, 0, false, -1 ); // unknown 2 - z
        }
        else
        {
            for (index_t i=0; i!=3; ++i)
            {
                bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, i ); // unknown 0 - x
                bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, i ); // unknown 1 - y
                bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, i ); // unknown 2 - z
                bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, i ); // unknown 2 - z
            }
        }

        tmp << 0,0,-1;
    }
    else if (testCase == 1)
    {
        if (weak)
        {
            GISMO_ERROR("Not implemented. Component-wise weak BCs need to be implemented first");

            // Diaphragm conditions
            bc.addCondition(boundary::west, condition_type::weak_dirichlet, 0, 0 ,false, 1 ); // unknown 1 - y
            bc.addCondition(boundary::west, condition_type::weak_dirichlet, 0, 0 ,false, 2 ); // unknown 2 - z

            bc.addCondition(boundary::east, condition_type::weak_dirichlet, 0, 0 ,false, 1 ); // unknown 1 - y
            bc.addCondition(boundary::east, condition_type::weak_dirichlet, 0, 0 ,false, 2 ); // unknown 2 - z
        }
        else
        {
            // Diaphragm conditions
            bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0 ,false, 1 ); // unknown 1 - y
            bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0 ,false, 2 ); // unknown 2 - z

            bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0 ,false, 1 ); // unknown 1 - y
            bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0 ,false, 2 ); // unknown 2 - z
        }

        bc.addCornerValue(boundary::southwest, 0.0, 0, 0); // (corner,value, patch, unknown)

        // Surface forces
        tmp << 0, 0, -90;
    }
    else if (testCase == 2)
    {
        if (weak)
        {
            GISMO_ERROR("Not implemented. Component-wise weak BCs need to be implemented first");

            bc.addCondition(boundary::north, condition_type::weak_dirichlet, 0, 0, false, 0 ); // unknown 0 - x
            bc.addCondition(boundary::north, condition_type::weak_dirichlet, 0, 0, false, 1 ); // unknown 1 - y
            bc.addCondition(boundary::north, condition_type::weak_dirichlet, 0, 0, false, 2 ); // unknown 2 - z

            // Symmetry in x-direction:
            bc.addCondition(boundary::east, condition_type::weak_dirichlet, 0, 0, false, 0 );
            bc.addCondition(boundary::east, condition_type::weak_clamped, 0, 0, false, 1 );
            bc.addCondition(boundary::east, condition_type::weak_clamped, 0, 0, false, 2 );

            // Symmetry in y-direction:
            bc.addCondition(boundary::west, condition_type::weak_clamped, 0, 0, false, 0 );
            bc.addCondition(boundary::west, condition_type::weak_dirichlet, 0, 0, false, 1 );
            bc.addCondition(boundary::west, condition_type::weak_clamped, 0, 0, false, 2 );

        }
        else
        {
            bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
            bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
            bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

            // Symmetry in x-direction:
            bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 );
            bc.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 1 );
            bc.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 2 );

            // Symmetry in y-direction:
            bc.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 0 );
            bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 );
            bc.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 2 );
        }

        // Surface forces
        tmp.setZero();

        // Point loads
        gsVector<> point(2);
        gsVector<> load (3);
        point<< 0.0, 0.0 ; load << 1.0, 0.0, 0.0 ;
        pLoads.addLoad(point, load, 0 );
        point<< 1.0, 0.0 ; load << 0.0, -1.0, 0.0 ;
        pLoads.addLoad(point, load, 0 );
    }
    else if (testCase == 3)
    {
        if (weak)
        {
            GISMO_ERROR("Not implemented. Component-wise weak BCs need to be implemented first");

            // Symmetry in y-direction for back side
            bc.addCondition(boundary::north, condition_type::weak_clamped, 0, 0, false, 0 );
            bc.addCondition(boundary::north, condition_type::weak_dirichlet, 0, 0, false, 1 );
            bc.addCondition(boundary::north, condition_type::weak_clamped, 0, 0, false, 2 );

            // Diaphragm conditions for left side
            bc.addCondition(boundary::west, condition_type::weak_dirichlet, 0, 0, false, 1 ); // unknown 1 - y
            bc.addCondition(boundary::west, condition_type::weak_dirichlet, 0, 0, false, 2 ); // unknown 2 - z

            // Symmetry in x-direction: for right side
            bc.addCondition(boundary::east, condition_type::weak_dirichlet, 0, 0, false, 0 );
            bc.addCondition(boundary::east, condition_type::weak_clamped, 0, 0, false, 1 );
            bc.addCondition(boundary::east, condition_type::weak_clamped, 0, 0, false, 2 );

            // Symmetry in z-direction:for the front side
            bc.addCondition(boundary::south, condition_type::weak_clamped, 0, 0, false, 0 );
            bc.addCondition(boundary::south, condition_type::weak_clamped, 0, 0, false, 1 );
            bc.addCondition(boundary::south, condition_type::weak_dirichlet, 0, 0, false, 2 );
        }
        else
        {
            // Symmetry in y-direction for back side
            bc.addCondition(boundary::north, condition_type::clamped, 0, 0, false, 0 );
            bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 );
            bc.addCondition(boundary::north, condition_type::clamped, 0, 0, false, 2 );

            // Diaphragm conditions for left side
            bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
            bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

            // Symmetry in x-direction: for right side
            bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 );
            bc.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 1 );
            bc.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 2 );

            // Symmetry in z-direction:for the front side
            bc.addCondition(boundary::south, condition_type::clamped, 0, 0, false, 0 );
            bc.addCondition(boundary::south, condition_type::clamped, 0, 0, false, 1 );
            bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 );
        }

        // Surface forces
        tmp.setZero();

        // Point loads
        gsVector<> point(2); point<< 1.0, 1.0 ;
        gsVector<> load (3); load << 0.0, 0.0, -0.25 ;
        pLoads.addLoad(point, load, 0 );
    }
    else if (testCase == 4) // Uniaxial tension; use with hyperelastic material model!
    {
        neu << 2625, 0, 0;
        neuData.setValue(neu,3);

        if (weak)
        {
            GISMO_ERROR("Not implemented. Component-wise weak BCs need to be implemented first");

            bc.addCondition(boundary::west, condition_type::weak_dirichlet, 0, 0, false, 0 ); // unknown 0 - x
            bc.addCondition(boundary::west, condition_type::weak_dirichlet, 0, 0, false, 2 ); // unknown 2 - z

            // bc.addCondition(boundary::east, condition_type::collapsed, 0, 0, false, 0 ); // unknown 1 - y
            bc.addCondition(boundary::east, condition_type::weak_dirichlet, 0, 0, false, 2 ); // unknown 2 - z.

            bc.addCondition(boundary::south, condition_type::weak_dirichlet, 0, 0, false, 1 ); // unknown 1 - y
            bc.addCondition(boundary::south, condition_type::weak_dirichlet, 0, 0, false, 2 ); // unknown 1 - y
            bc.addCondition(boundary::north, condition_type::weak_dirichlet, 0, 0, false, 2 ); // unknown 1 - y

        }
        else
        {
            bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
            bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

            // bc.addCondition(boundary::east, condition_type::collapsed, 0, 0, false, 0 ); // unknown 1 - y
            bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z.

            bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
            bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 1 - y
            bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 1 - y

        }

        bc.addCondition(boundary::east, condition_type::neumann, &neuData ); // unknown 1 - y

        gsVector<> point(2);
        gsVector<> load (3);
        point<< 1.0, 0.5 ;
        load << 1.0,0.0, 0.0;
        pLoads.addLoad(point, load, 0 );
    }
    else if (testCase == 5)
    {
        if (weak)
        {
            bc.addCondition(boundary::north, condition_type::weak_dirichlet, &weak_drch );
            bc.addCondition(boundary::east, condition_type::weak_dirichlet, &weak_drch );
            bc.addCondition(boundary::south, condition_type::weak_dirichlet, &weak_drch );
            bc.addCondition(boundary::west, condition_type::weak_dirichlet, &weak_drch );
        }
        else
        {
            for (index_t i = 0; i!=3; i++)
            {
                bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, i );
                bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, i );
                bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, i );
                bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, i );
            }
        }

        // Point loads
        gsVector<> point(2); point<< 0.5,0.5 ;
        gsVector<> load (3); load << 0.0, 0.0, -1 ;
        pLoads.addLoad(point, load, 0 );
    }
    else if (testCase == 6) // balloon
    {
        if (weak)
        {
            GISMO_ERROR("Not implemented. Component-wise weak BCs need to be implemented first");

            bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

            // Symmetry in x-direction:
            bc.addCondition(boundary::east, condition_type::weak_dirichlet, 0, 0, false, 0 );
            bc.addCondition(boundary::east, condition_type::weak_clamped, 0, 0, false, 1 );
            bc.addCondition(boundary::east, condition_type::weak_clamped, 0, 0, false, 2 );

            // Symmetry in y-direction:
            bc.addCondition(boundary::west, condition_type::weak_clamped, 0, 0, false, 0 );
            bc.addCondition(boundary::west, condition_type::weak_dirichlet, 0, 0, false, 1 );
            bc.addCondition(boundary::west, condition_type::weak_clamped, 0, 0, false, 2 );
        }
        else
        {
            bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

            // Symmetry in x-direction:
            bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 );
            bc.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 1 );
            bc.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 2 );

            // Symmetry in y-direction:
            bc.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 0 );
            bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 );
            bc.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 2 );
        }

        // Pressure
        pressure = 5e3;
    }
    else if (testCase == 7)
    {
        neu << 0, 0, -100.0;
        neuData.setValue(neu,3);

        bc.addCondition(boundary::north, condition_type::neumann, &neuData );

        if (weak)
        {
            GISMO_ERROR("Not implemented. Component-wise weak BCs need to be implemented first");

            bc.addCondition(boundary::south, condition_type::weak_dirichlet, 0, 0, false, 0 ); // unknown 0 - x
            bc.addCondition(boundary::south, condition_type::weak_dirichlet, 0, 0, false, 1 ); // unknown 1 - y
            bc.addCondition(boundary::south, condition_type::weak_dirichlet, 0, 0, false, 2 ); // unknown 2 - z

            // Symmetry in x-direction:
            bc.addCondition(boundary::east, condition_type::weak_dirichlet, 0, 0, false, 0 );
            bc.addCondition(boundary::east, condition_type::weak_clamped, 0, 0, false, 1 );
            bc.addCondition(boundary::east, condition_type::weak_clamped, 0, 0, false, 2 );

            // Symmetry in y-direction:
            bc.addCondition(boundary::west, condition_type::weak_clamped, 0, 0, false, 0 );
            bc.addCondition(boundary::west, condition_type::weak_dirichlet, 0, 0, false, 1 );
            bc.addCondition(boundary::west, condition_type::weak_clamped, 0, 0, false, 2 );
        }
        else
        {
            bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
            bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
            bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

            // Symmetry in x-direction:
            bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 );
            bc.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 1 );
            bc.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 2 );

            // Symmetry in y-direction:
            bc.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 0 );
            bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 );
            bc.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 2 );
        }

    }
    // else if (testCase == 8)
    // {
    //     bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
    //     bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y

    //     // Symmetry in x-direction:
    //     neu << 0, 6.25, 0;
    //     neuData.setValue(neu,3);
    //     bc.addCondition(boundary::east, condition_type::neumann, &neuData );

    //     bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 2 );
    //     bc.addCondition(boundary::east,  condition_type::dirichlet, 0, 0, false, 2 );
    //     bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 );
    //     bc.addCondition(boundary::west,  condition_type::dirichlet, 0, 0, false, 2 );
    // }
    else if (testCase == 8)
    {
    //     real_t Load = -0.0168288;
    //     neu << 0, 0, Load;
    //     neuData.setValue(neu,3);

    //     bc.addCondition(boundary::north, condition_type::neumann, &neuData );
    //     bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 2 - z
    //     bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 2 - z
    //     bc.addCondition(boundary::north, condition_type::collapsed, 0, 0, false, 2 ); // unknown 1 - y

        displx.setValue(-0.027815,3);
        if (weak)
        {
            GISMO_ERROR("Not implemented. Component-wise weak BCs need to be implemented first");

            bc.addCondition(boundary::north, condition_type::weak_dirichlet, 0, 0, false, 0 ); // unknown 2 - z
            bc.addCondition(boundary::north, condition_type::weak_dirichlet, 0, 0, false, 1 ); // unknown 2 - z
            bc.addCondition(boundary::north, condition_type::weak_dirichlet, &displx, 0, false, 2 ); // unknown 1 - y

            bc.addCondition(boundary::south, condition_type::weak_dirichlet, 0, 0, false, 0 ); // unknown 0 - x
            bc.addCondition(boundary::south, condition_type::weak_dirichlet, 0, 0, false, 1 ); // unknown 1 - y
            bc.addCondition(boundary::south, condition_type::weak_dirichlet, 0, 0, false, 2 ); // unknown 2 - z

            // Symmetry in x-direction:
            bc.addCondition(boundary::east, condition_type::weak_dirichlet, 0, 0, false, 0 );
            bc.addCondition(boundary::east, condition_type::weak_clamped, 0, 0, false, 1 );
            bc.addCondition(boundary::east, condition_type::weak_clamped, 0, 0, false, 2 );

            // Symmetry in y-direction:
            bc.addCondition(boundary::west, condition_type::weak_clamped, 0, 0, false, 0 );
            bc.addCondition(boundary::west, condition_type::weak_dirichlet, 0, 0, false, 1 );
            bc.addCondition(boundary::west, condition_type::weak_clamped, 0, 0, false, 2 );
        }
        else
        {
            bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 2 - z
            bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 2 - z
            bc.addCondition(boundary::north, condition_type::dirichlet, &displx, 0, false, 2 ); // unknown 1 - y

            bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
            bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
            bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

            // Symmetry in x-direction:
            bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 );
            bc.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 1 );
            bc.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 2 );

            // Symmetry in y-direction:
            bc.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 0 );
            bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 );
            bc.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 2 );
        }
    }
    else if (testCase == 9)
    {
        real_t  Load = -0.0168288;
        neu << 0, 0, Load;
        neuData.setValue(neu,3);

        bc.addCondition(boundary::north, condition_type::neumann, &neuData );
        bc.addCondition(boundary::north, condition_type::collapsed, 0, 0, false, 2 ); // unknown 1 - y

        if (weak)
        {
            GISMO_ERROR("Not implemented. Component-wise weak BCs need to be implemented first");

            bc.addCondition(boundary::south, condition_type::weak_dirichlet, 0, 0, false, 0 ); // unknown 0 - x
            bc.addCondition(boundary::south, condition_type::weak_dirichlet, 0, 0, false, 1 ); // unknown 1 - y
            bc.addCondition(boundary::south, condition_type::weak_dirichlet, 0, 0, false, 2 ); // unknown 2 - z

            // Symmetry in x-direction:
            bc.addCondition(boundary::east, condition_type::weak_dirichlet, 0, 0, false, 0 );
            bc.addCondition(boundary::east, condition_type::weak_clamped, 0, 0, false, 1 );
            bc.addCondition(boundary::east, condition_type::weak_clamped, 0, 0, false, 2 );

            // Symmetry in y-direction:
            bc.addCondition(boundary::west, condition_type::weak_clamped, 0, 0, false, 0 );
            bc.addCondition(boundary::west, condition_type::weak_dirichlet, 0, 0, false, 1 );
            bc.addCondition(boundary::west, condition_type::weak_clamped, 0, 0, false, 2 );
        }
        else
        {
            bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
            bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
            bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

            // Symmetry in x-direction:
            bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 );
            bc.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 1 );
            bc.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 2 );

            // Symmetry in y-direction:
            bc.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 0 );
            bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 );
            bc.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 2 );
        }
    }
    else if (testCase == 10)
    {
        if (weak)
        {
            GISMO_ERROR("Not implemented. Component-wise weak BCs need to be implemented first");

            for (index_t i=0; i!=3; ++i)
            {
                bc.addCondition(boundary::west, condition_type::weak_dirichlet, 0, 0, false, i ); // unknown 2 - z
            }
            bc.addCondition(boundary::north, condition_type::weak_dirichlet, 0, 0, false, 2 ); // unknown 2 - z
            bc.addCondition(boundary::south, condition_type::weak_dirichlet, 0, 0, false, 2 ); // unknown 2 - z

            bc.addCondition(boundary::east, condition_type::weak_dirichlet, 0, 0 ,false,1);
            bc.addCondition(boundary::east, condition_type::weak_dirichlet, 0, 0 ,false,2);
        }
        else
        {
            for (index_t i=0; i!=3; ++i)
            {
                bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, i ); // unknown 2 - z
            }
            bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
            bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

            bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0 ,false,1);
            bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0 ,false,2);
        }
        bc.addCondition(boundary::east, condition_type::collapsed, 0, 0 ,false,0);

        gsVector<> point(2); point<< 1.0, 0.5 ;
        gsVector<> load (3); load << 0.25, 0.0, 0.0 ;
        pLoads.addLoad(point, load, 0 );
    }
    else if (testCase == 11)
    {
        if (weak)
        {
            GISMO_ERROR("Not implemented. Component-wise weak BCs need to be implemented first");

            for (index_t i=0; i!=3; ++i)
            {
                bc.addCondition(boundary::east, condition_type::weak_dirichlet, 0, 0, false, i ); // unknown 1 - y
                bc.addCondition(boundary::west, condition_type::weak_dirichlet, 0, 0, false, i ); // unknown 2 - z
            }

            bc.addCondition(boundary::east, condition_type::weak_clamped, 0, 0 ,false,2);
            bc.addCondition(boundary::west, condition_type::weak_clamped, 0, 0 ,false,2);
        }
        else
        {
            for (index_t i=0; i!=3; ++i)
            {
                bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, i ); // unknown 1 - y
                bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, i ); // unknown 2 - z
            }

            bc.addCondition(boundary::east, condition_type::clamped, 0, 0 ,false,2);
            bc.addCondition(boundary::west, condition_type::clamped, 0, 0 ,false,2);
        }

        tmp<<0,0,-1;
    }
    else if (testCase == 12)
    {
        if (weak)
        {
            GISMO_ERROR("Not implemented. Component-wise weak BCs need to be implemented first");

            for (index_t i=0; i!=3; ++i)
            {
                bc.addCondition(boundary::west, condition_type::weak_dirichlet, 0, 0, false, i ); // unknown 2 - z
            }

            bc.addCondition(boundary::west, condition_type::weak_clamped, 0, 0 ,false,2);
        }
        else
        {
            for (index_t i=0; i!=3; ++i)
            {
                bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, i ); // unknown 2 - z
            }

            bc.addCondition(boundary::west, condition_type::clamped, 0, 0 ,false,2);
        }


        neu << 0, 0, -0.1;
        neuData.setValue(neu,3);

        bc.addCondition(boundary::east, condition_type::neumann, &neuData );
    }
    else if (testCase == 13)
    {
        if (weak)
        {
            weak_drch.setValue(weak_D,3);
            weak_clmp.setValue(weak_C,3);
            bc.addCondition(boundary::north, condition_type::weak_dirichlet, &weak_drch );
            bc.addCondition(boundary::east, condition_type::weak_dirichlet, &weak_drch );
            bc.addCondition(boundary::south, condition_type::weak_dirichlet, &weak_drch );
            bc.addCondition(boundary::west, condition_type::weak_dirichlet, &weak_drch );

            bc.addCondition(boundary::north, condition_type::weak_clamped, &weak_clmp);
            bc.addCondition(boundary::east, condition_type::weak_clamped, &weak_clmp);
            bc.addCondition(boundary::south, condition_type::weak_clamped, &weak_clmp);
            bc.addCondition(boundary::west, condition_type::weak_clamped, &weak_clmp);
        }
        else
        {
            for (index_t i=0; i!=3; ++i)
            {
                bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, i ); // unknown 0 - x
                bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, i ); // unknown 1 - y
                bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, i ); // unknown 2 - z
                bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, i ); // unknown 2 - z
            }

            bc.addCondition(boundary::north, condition_type::clamped, 0, 0 ,false,2);
            bc.addCondition(boundary::east, condition_type::clamped, 0, 0 ,false,2);
            bc.addCondition(boundary::south, condition_type::clamped, 0, 0 ,false,2);
            bc.addCondition(boundary::west, condition_type::clamped, 0, 0 ,false,2);
        }

        tmp << 0,0,-1;
    }
    else if (testCase == 14)
    {
        if (weak)
        {
            GISMO_ERROR("Not implemented. Component-wise weak BCs need to be implemented first");

            for (index_t i=0; i!=3; ++i)
            {
                bc.addCondition(boundary::west, condition_type::weak_dirichlet, 0, 0, false, i ); // unknown 2 - z
            }
            bc.addCondition(boundary::north, condition_type::weak_dirichlet, 0, 0, false, 2 ); // unknown 2 - z
            bc.addCondition(boundary::south, condition_type::weak_dirichlet, 0, 0, false, 2 ); // unknown 2 - z
        }
        else
        {
            for (index_t i=0; i!=3; ++i)
            {
                bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, i ); // unknown 2 - z
            }
            bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
            bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
        }
        tmp << 0,0,-1;

        bc.addCondition(boundary::east, condition_type::collapsed, 0, 0 ,false,0);
        bc.addCondition(boundary::east, condition_type::collapsed, 0, 0, false, 2 );
    }
    else if (testCase == 14)
    {
        if (weak)
        {
            GISMO_ERROR("Not implemented. Component-wise weak BCs need to be implemented first");

            for (index_t i=0; i!=3; ++i)
            {
                bc.addCondition(boundary::west, condition_type::weak_dirichlet, 0, 0, false, i ); // unknown 2 - z
            }
            bc.addCondition(boundary::north, condition_type::weak_dirichlet, 0, 0, false, 2 ); // unknown 2 - z
            bc.addCondition(boundary::south, condition_type::weak_dirichlet, 0, 0, false, 2 ); // unknown 2 - z
        }
        else
        {
            for (index_t i=0; i!=3; ++i)
            {
                bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, i ); // unknown 2 - z
            }
            bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
            bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
        }

        tmp << 0,0,-1;

        bc.addCondition(boundary::east, condition_type::collapsed, 0, 0 ,false,0);
        bc.addCondition(boundary::east, condition_type::collapsed, 0, 0, false, 2 );
    }
    else if (testCase == 15)
    {
        if (weak)
        {
            GISMO_ERROR("Not implemented. Component-wise weak BCs need to be implemented first");

            for (index_t i=0; i!=3; ++i)
            {
                bc.addCondition(boundary::west, condition_type::weak_dirichlet, 0, 0, false, i ); // unknown 2 - z
            }
            bc.addCondition(boundary::north, condition_type::weak_dirichlet, 0, 0, false, 2 ); // unknown 2 - z
            bc.addCondition(boundary::south, condition_type::weak_dirichlet, 0, 0, false, 2 ); // unknown 2 - z


            bc.addCondition(boundary::east, condition_type::weak_dirichlet, 0, 0 ,false,1);
            bc.addCondition(boundary::east, condition_type::weak_dirichlet, 0, 0 ,false,2);
        }
        else
        {
            for (index_t i=0; i!=3; ++i)
            {
                bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, i ); // unknown 2 - z
            }
            bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
            bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z


            bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0 ,false,1);
            bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0 ,false,2);
        }
        bc.addCondition(boundary::east, condition_type::collapsed, 0, 0 ,false,0);

        gsVector<> point(2); point<< 1.0, 0.5 ;
        gsVector<> load (3); load << 0.1, 0.0, 0.0 ;
        pLoads.addLoad(point, load, 0 );
    }
    else if (testCase == 16)
    {
        if (weak)
        {
            GISMO_ERROR("Not implemented. Component-wise weak BCs need to be implemented first");

            bc.addCondition(boundary::west, condition_type::weak_dirichlet, 0, 0, false, 0 ); // unknown 0 - x
            bc.addCondition(boundary::west, condition_type::weak_dirichlet, 0, 0, false, 1 ); // unknown 1 - y
            bc.addCondition(boundary::west, condition_type::weak_dirichlet, 0, 0, false, 2 ); // unknown 2 - z

            bc.addCondition(boundary::west, condition_type::weak_clamped, 0, 0, false, 2 );
        }
        else
        {
            bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
            bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
            bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

            bc.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 2 );
        }

        real_t Load = 1e-2;
        gsVector<> point(2);
        gsVector<> load (3);
        point<< 1.0, 0.5 ;
        load << 0.0, 0.0, Load ;
        pLoads.addLoad(point, load, 0 );
    }
    else if (testCase == 17)
    {
        if (weak)
        {
            GISMO_ERROR("Not implemented. Component-wise weak BCs need to be implemented first");

            bc.addCondition(boundary::west, condition_type::weak_dirichlet, 0, 0, false, 1 ); // unknown 1 - y
            bc.addCondition(boundary::west, condition_type::weak_dirichlet, 0, 0, false, 2 ); // unknown 2 - z

            bc.addCondition(boundary::east, condition_type::weak_dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        }
        else
        {
            bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
            bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

            bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        }
        real_t Load = 1e-1;
        neu << -Load, 0, 0;
        neuData.setValue(neu,3);

        bc.addCondition(boundary::west, condition_type::neumann, &neuData ); // unknown 0 - x
        bc.addCondition(boundary::west, condition_type::collapsed, 0, 0, false, 0 ); // unknown 0 - x
    }
    else if (testCase == 20)
    {
        real_t Load = 1e-1;
        neu << 0, 0, -Load;
        neuData.setValue(neu,3);

        real_t T1 = 0.5;
        std::string nx = std::to_string(-T1) + "*cos(atan2(y,x))";
        std::string ny = std::to_string(-T1) + "*sin(atan2(y,x))";
        std::string nz = "0";
        neuDataFun1 = gsFunctionExpr<>(nx,ny,nz,3);

        // for (index_t p=0; p!=mp.nPatches(); p++)
        // {
        //     // bc.addCondition(p,boundary::west, condition_type::neumann, &neuDataFun1 ); // unknown 0 - x
        //     bc.addCondition(p,boundary::west, condition_type::neumann, &neuData ); // unknown 0 - x
        //     bc.addCondition(p,boundary::east, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 1 - y
        //     bc.addCondition(p,boundary::east, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        //     bc.addCondition(p,boundary::east, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
        // }

        bc.addCondition(0,boundary::south, condition_type::neumann, &neuDataFun1 ); // unknown 0 - x
        // bc.addCondition(0,boundary::south, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        // bc.addCondition(0,boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 0 - x
        // bc.addCondition(0,boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 0 - x
        bc.addCondition(0,boundary::west, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        bc.addCondition(0,boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 0 - x
        bc.addCondition(0,boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 0 - x

        // bc.addCondition(1,boundary::west, condition_type::neumann, &neuDataFun1 ); // unknown 0 - x
        bc.addCondition(1,boundary::west, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        bc.addCondition(1,boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 0 - x
        bc.addCondition(1,boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 0 - x
        bc.addCondition(1,boundary::north, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        bc.addCondition(1,boundary::north, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 0 - x
        bc.addCondition(1,boundary::north, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 0 - x

        // bc.addCondition(2,boundary::south, condition_type::neumann, &neuDataFun1 ); // unknown 0 - x
        bc.addCondition(2,boundary::south, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        bc.addCondition(2,boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 0 - x
        bc.addCondition(2,boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 0 - x
        bc.addCondition(2,boundary::east, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        bc.addCondition(2,boundary::east, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 0 - x
        bc.addCondition(2,boundary::east, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 0 - x

        // bc.addCondition(3,boundary::east, condition_type::neumann, &neuDataFun1 ); // unknown 0 - x
        bc.addCondition(3,boundary::east, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        bc.addCondition(3,boundary::east, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 0 - x
        bc.addCondition(3,boundary::east, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 0 - x
        bc.addCondition(3,boundary::north, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        bc.addCondition(3,boundary::north, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 0 - x
        bc.addCondition(3,boundary::north, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 0 - x


        tmp << 0,0,-1e-6;
    }

    gsDebugVar(bc);
    //! [Refinement]

    // Linear isotropic material model
    gsConstantFunction<> force(tmp,3);
    gsFunctionExpr<> t(std::to_string(thickness), 3);
    gsFunctionExpr<> E(std::to_string(E_modulus),3);
    gsFunctionExpr<> nu(std::to_string(PoissonRatio),3);
    gsFunctionExpr<> rho(std::to_string(Density),3);
    gsConstantFunction<> ratio(Ratio,3);

    real_t mu = E_modulus / (2 * (1 + PoissonRatio));
    gsConstantFunction<> alpha1(1.3,3);
    gsConstantFunction<> mu1(6.3e5/4.225e5*mu,3);
    gsConstantFunction<> alpha2(5.0,3);
    gsConstantFunction<> mu2(0.012e5/4.225e5*mu,3);
    gsConstantFunction<> alpha3(-2.0,3);
    gsConstantFunction<> mu3(-0.1e5/4.225e5*mu,3);

    real_t pi = math::atan(1)*4;
    index_t kmax = 2;

    std::vector<gsFunctionSet<> * > Gs(kmax);
    std::vector<gsFunctionSet<> * > Ts(kmax);
    std::vector<gsFunctionSet<> * > Phis(kmax);

    gsMatrix<> Gmat = gsCompositeMatrix(E_modulus,E_modulus,0.5 * E_modulus / (1+PoissonRatio),PoissonRatio,PoissonRatio);
    Gmat.resize(Gmat.rows()*Gmat.cols(),1);
    gsConstantFunction<> Gfun(Gmat,3);
    Gs[0] = Gs[1] = &Gfun;

    gsConstantFunction<> phi1, phi2;
    phi1.setValue(0,3);
    phi2.setValue(pi/2.0,3);

    Phis[0] = &phi1;
    Phis[1] = &phi2;

    gsConstantFunction<> thicks(thickness/kmax,3);
    Ts[0] = Ts[1] = &thicks;

    std::vector<gsFunction<>*> parameters;
    if (material==0) // SvK & Composites
    {
      parameters.resize(2);
      parameters[0] = &E;
      parameters[1] = &nu;
    }
    else if (material==1 || material==2) // NH & NH_ext
    {
      parameters.resize(2);
      parameters[0] = &E;
      parameters[1] = &nu;
    }
    else if (material==3) // MR
    {
      parameters.resize(3);
      parameters[0] = &E;
      parameters[1] = &nu;
      parameters[2] = &ratio;
    }
    else if (material==4) // OG
    {
      parameters.resize(8);
      parameters[0] = &E;
      parameters[1] = &nu;
      parameters[2] = &mu1;
      parameters[3] = &alpha1;
      parameters[4] = &mu2;
      parameters[5] = &alpha2;
      parameters[6] = &mu3;
      parameters[7] = &alpha3;
    }

    gsMaterialMatrixBase<real_t>* materialMatrix;

    gsOptionList options;
    if      (material==0 && impl==1)
    {
        if (composite)
        {
            materialMatrix = new gsMaterialMatrixComposite<3,real_t>(mp,Ts,Gs,Phis);
        }
        else
        {
            parameters.resize(2);
            options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
            options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
            materialMatrix = getMaterialMatrix<3,real_t>(mp,t,parameters,rho,options);
        }
    }
    else
    {
        options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",material);
        options.addSwitch("Compressibility","Compressibility: (false): Imcompressible | (true): Compressible",Compressibility);
        options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",impl);
        materialMatrix = getMaterialMatrix<3,real_t>(mp,t,parameters,rho,options);
    }

    gsThinShellAssemblerBase<real_t>* assembler;
    if(membrane)
        assembler = new gsThinShellAssembler<3, real_t, false>(mp,dbasis,bc,force,materialMatrix);
    else
        assembler = new gsThinShellAssembler<3, real_t, true >(mp,dbasis,bc,force,materialMatrix);

    opts.addReal("WeakDirichlet","Penalty parameter weak dirichlet conditions",1e5);
    opts.addReal("WeakClamped","Penalty parameter weak clamped conditions",1e5);
    // Construct assembler object
    assembler->setOptions(opts);
    assembler->setPointLoads(pLoads);

    // gsVector<> found_vec(3);
    // found_vec<<0,0,2;
    // gsConstantFunction<> found(found_vec,3);
    // assembler->setFoundation(found);

    gsStopwatch stopwatch,stopwatch2;
    real_t time = 0.0;
    real_t totaltime = 0.0;

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

    // Define Matrices
    stopwatch.restart();
    stopwatch2.restart();
    assembler->assemble();
    time += stopwatch.stop();

    gsSparseMatrix<> matrix = assembler->matrix();
    gsVector<> vector = assembler->rhs();

    // Configure Structural Analsysis module
    gsStaticSolver<real_t> staticSolver(matrix,vector,Jacobian,Residual);
    gsOptionList solverOptions = staticSolver.options();
    solverOptions.setInt("Verbose",verbose);
    solverOptions.setInt("MaxIterations",10);
    solverOptions.setReal("Tolerance",1e-6);
    staticSolver.setOptions(solverOptions);

    // Solve linear problem
    gsVector<> solVector;
    solVector = staticSolver.solveLinear();
    if (nonlinear)
        solVector = staticSolver.solveNonlinear();

    totaltime += stopwatch2.stop();

    mp_def = assembler->constructSolution(solVector);

    gsMultiPatch<> deformation = mp_def;
    for (size_t k = 0; k != mp_def.nPatches(); ++k)
        deformation.patch(k).coefs() -= mp.patch(k).coefs();

    if (testCase==4)
    {
        gsVector<> pt(2);
        pt<<1,0;
      gsMatrix<> lambdas = assembler->computePrincipalStretches(pt,mp_def,0);
      real_t S = 2625 / 1e-3 / lambdas(0) / lambdas(2);
      real_t San = mu * (math::pow(lambdas(1),2)-1/lambdas(1));
      gsInfo<<"S = \t"<<S<<"\t San = \t"<<San<<"\t |S-San| = \t"<<abs(S-San)<<"\n";
      gsInfo<<"lambda = \t"<<lambdas(1)<<"\t 1/lambda = \t"<<1/lambdas(1)<<"\t lambda_0 = \t"<<lambdas(0)<<"\t lambda_2 = \t"<<lambdas(2)<<"\n";
    }


    // ! [Export visualization in ParaView]
    if (plot)
    {
        gsField<> solField(mp_def, deformation);
        gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview<>( solField, "solution", 1000, true);
        // ev.options().setSwitch("plot.elements", true);
        // ev.writeParaview( u_sol   , G, "solution");

        // gsFileManager::open("solution.pvd");

        gsInfo <<"Maximum deformation coef: "
               << deformation.patch(0).coefs().colwise().maxCoeff() <<".\n";
        gsInfo <<"Minimum deformation coef: "
               << deformation.patch(0).coefs().colwise().minCoeff() <<".\n";
    }
    if (stress)
    {

        gsPiecewiseFunction<> membraneStresses;
        assembler->constructStress(mp_def,membraneStresses,stress_type::membrane);
        gsField<> membraneStress(mp_def,membraneStresses, true);

        gsPiecewiseFunction<> flexuralStresses;
        assembler->constructStress(mp_def,flexuralStresses,stress_type::flexural);
        gsField<> flexuralStress(mp_def,flexuralStresses, true);

        gsPiecewiseFunction<> stretches;
        assembler->constructStress(mp_def,stretches,stress_type::principal_stretch);
        gsField<> Stretches(mp_def,stretches, true);

        // gsPiecewiseFunction<> membraneStresses_p;
        // assembler->constructStress(mp_def,membraneStresses_p,stress_type::principal_stress_membrane);
        // gsField<> membraneStress_p(mp_def,membraneStresses_p, true);

        // gsPiecewiseFunction<> flexuralStresses_p;
        // assembler->constructStress(mp_def,flexuralStresses_p,stress_type::principal_stress_flexural);
        // gsField<> flexuralStress_p(mp_def,flexuralStresses_p, true);

        gsPiecewiseFunction<> stretch1;
        assembler->constructStress(mp_def,stretch1,stress_type::principal_stretch_dir1);
        gsField<> stretchDir1(mp_def,stretch1, true);

        gsPiecewiseFunction<> stretch2;
        assembler->constructStress(mp_def,stretch2,stress_type::principal_stretch_dir2);
        gsField<> stretchDir2(mp_def,stretch2, true);

        gsPiecewiseFunction<> stretch3;
        assembler->constructStress(mp_def,stretch3,stress_type::principal_stretch_dir3);
        gsField<> stretchDir3(mp_def,stretch3, true);


        gsField<> solutionField(mp,deformation, true);


        // gsField<> stressField = assembler->constructStress(mp_def,stress_type::membrane_strain);

        #ifdef GISMO_ELASTICITY
        std::map<std::string,const gsField<> *> fields;
        fields["Deformation"] = &solutionField;
        fields["Membrane Stress"] = &membraneStress;
        fields["Flexural Stress"] = &flexuralStress;
        fields["Principal Stretch"] = &Stretches;
        // fields["Principal Membrane Stress"] = &membraneStress_p;
        // fields["Principal Flexural Stress"] = &flexuralStress_p;
        fields["Principal Direction 1"] = &stretchDir1;
        fields["Principal Direction 2"] = &stretchDir2;
        fields["Principal Direction 3"] = &stretchDir3;

        gsWriteParaviewMultiPhysics(fields,"stress",5000,true);
        #else
        gsWriteParaview(solutionField, "Deformation");
        gsWriteParaview(membraneStress, "MembraneStress");
        gsWriteParaview(flexuralStress, "FlexuralStress");
        gsWriteParaview(Stretches, "PrincipalStretch");
        gsWriteParaview(stretchDir1, "PrincipalDirection1");
        gsWriteParaview(stretchDir2, "PrincipalDirection2");
        gsWriteParaview(stretchDir3, "PrincipalDirection3");
        #endif

    }
    gsInfo<<"Total ellapsed assembly time: \t\t"<<time<<" s\n";
    gsInfo<<"Total ellapsed solution time (incl. assembly): \t"<<totaltime<<" s\n";

    return EXIT_SUCCESS;

}// end main

template <class T>
gsMultiPatch<T> RectangularDomain(int n, int p, T L, T B)
{
  gsMultiPatch<T> mp = RectangularDomain(n, n, p, p, L, B);
  return mp;
}

template <class T>
gsMultiPatch<T> RectangularDomain(int n, int m, int p, int q, T L, T B)
{
  // -------------------------------------------------------------------------
  // --------------------------Make beam geometry-----------------------------
  // -------------------------------------------------------------------------
  int dim = 3; //physical dimension
  gsKnotVector<> kv0;
  kv0.initUniform(0,1,0,p+1,1);
  gsKnotVector<> kv1;
  kv1.initUniform(0,1,0,q+1,1);

  for(index_t i = 0; i< n; ++i)
      kv0.uniformRefine();
  for(index_t i = 0; i< m; ++i)
      kv1.uniformRefine();

  // Make basis
  gsTensorBSplineBasis<2,T> basis(kv0,kv1);

  // Initiate coefficient matrix
  gsMatrix<> coefs(basis.size(),dim);
  // Number of control points needed per component
  size_t len0 = basis.component(0).size();
  size_t len1 = basis.component(1).size();
  gsVector<> coefvec0(len0);
  // Uniformly distribute control points per component
  coefvec0.setLinSpaced(len0,0.0,L);
  gsVector<> coefvec1(basis.component(1).size());
  coefvec1.setLinSpaced(len1,0.0,B);

  // Z coordinate is zero
  coefs.col(2).setZero();

  // Define a matrix with ones
  gsVector<> temp(len0);
  temp.setOnes();
  for (size_t k = 0; k < len1; k++)
  {
    // First column contains x-coordinates (length)
    coefs.col(0).segment(k*len0,len0) = coefvec0;
    // Second column contains y-coordinates (width)
    coefs.col(1).segment(k*len0,len0) = temp*coefvec1.at(k);
  }
  // Create gsGeometry-derived object for the patch
  gsTensorBSpline<2,real_t> shape(basis,coefs);

  gsMultiPatch<T> mp;
  mp.addPatch(shape);
  mp.addAutoBoundaries();

  return mp;
}

template <class T>
gsMultiPatch<T> Rectangle(T L, T B)
{
  // -------------------------------------------------------------------------
  // --------------------------Make beam geometry-----------------------------
  // -------------------------------------------------------------------------
  int dim = 3; //physical dimension
  gsKnotVector<> kv0;
  kv0.initUniform(0,1,0,2,1);
  gsKnotVector<> kv1;
  kv1.initUniform(0,1,0,2,1);

  // Make basis
  gsTensorBSplineBasis<2,T> basis(kv0,kv1);

  // Initiate coefficient matrix
  gsMatrix<> coefs(basis.size(),dim);
  // Number of control points needed per component
  size_t len0 = basis.component(0).size();
  size_t len1 = basis.component(1).size();
  gsVector<> coefvec0(len0);
  // Uniformly distribute control points per component
  coefvec0.setLinSpaced(len0,0.0,L);
  gsVector<> coefvec1(basis.component(1).size());
  coefvec1.setLinSpaced(len1,0.0,B);

  // Z coordinate is zero
  coefs.col(2).setZero();

  // Define a matrix with ones
  gsVector<> temp(len0);
  temp.setOnes();
  for (size_t k = 0; k < len1; k++)
  {
    // First column contains x-coordinates (length)
    coefs.col(0).segment(k*len0,len0) = coefvec0;
    // Second column contains y-coordinates (width)
    coefs.col(1).segment(k*len0,len0) = temp*coefvec1.at(k);
  }
  // Create gsGeometry-derived object for the patch
  gsTensorBSpline<2,real_t> shape(basis,coefs);

  gsMultiPatch<T> mp;
  mp.addPatch(shape);
  mp.addAutoBoundaries();

  return mp;
}


template <class T>
gsMultiPatch<T> RectangularDomainVert(int n, int p, T L, T B)
{
  gsMultiPatch<T> mp = RectangularDomainVert(n, n, p, p, L, B);
  return mp;
}

template <class T>
gsMultiPatch<T> RectangularDomainVert(int n, int m, int p, int q, T L, T B)
{
  // -------------------------------------------------------------------------
  // --------------------------Make beam geometry-----------------------------
  // -------------------------------------------------------------------------
  int dim = 3; //physical dimension
  gsKnotVector<> kv0;
  kv0.initUniform(0,1,0,p+1,1);
  gsKnotVector<> kv1;
  kv1.initUniform(0,1,0,q+1,1);

  for(index_t i = 0; i< n; ++i)
      kv0.uniformRefine();
  for(index_t i = 0; i< m; ++i)
      kv1.uniformRefine();

  // Make basis
  gsTensorBSplineBasis<2,T> basis(kv0,kv1);

  // Initiate coefficient matrix
  gsMatrix<> coefs(basis.size(),dim);
  // Number of control points needed per component
  size_t len0 = basis.component(0).size();
  size_t len1 = basis.component(1).size();
  gsVector<> coefvec0(len0);
  // Uniformly distribute control points per component
  coefvec0.setLinSpaced(len0,0.0,L);
  gsVector<> coefvec1(basis.component(1).size());
  coefvec1.setLinSpaced(len1,0.0,B);

  // Z coordinate is zero
  coefs.col(2).setZero();

  // Define a matrix with ones
  gsVector<> temp(len0);
  temp.setOnes();
  for (index_t k = 0; k < len1; k++)
  {
    // First column contains x-coordinates (length)
    coefs.col(0).segment(k*len0,len0) = coefvec0;
    // Second column contains z-coordinates (height)
    coefs.col(2).segment(k*len0,len0) = temp*coefvec1.at(k);
  }
  // Create gsGeometry-derived object for the patch
  gsTensorBSpline<2,real_t> shape(basis,coefs);

  gsMultiPatch<T> mp;
  mp.addPatch(shape);
  mp.addAutoBoundaries();

  return mp;
}

template <class T>
gsMultiPatch<T> RectangularDomain90(int n, int p, T L, T B)
{
  gsMultiPatch<T> mp = RectangularDomain90(n, n, p, p, L, B);
  return mp;
}

template <class T>
gsMultiPatch<T> RectangularDomain90(int n, int m, int p, int q, T L, T B)
{
  // -------------------------------------------------------------------------
  // --------------------------Make beam geometry-----------------------------
  // -------------------------------------------------------------------------
  int dim = 3; //physical dimension
  gsKnotVector<> kv0;
  kv0.initUniform(0,1,0,p+1,1);
  gsKnotVector<> kv1;
  kv1.initUniform(0,1,0,q+1,1);

  for(index_t i = 0; i< n; ++i)
      kv0.uniformRefine();
  for(index_t i = 0; i< m; ++i)
      kv1.uniformRefine();

  // Make basis
  gsTensorBSplineBasis<2,T> basis(kv0,kv1);

  // Initiate coefficient matrix
  gsMatrix<> coefs(basis.size(),dim);
  // Number of control points needed per component
  size_t len0 = basis.component(0).size();
  size_t len1 = basis.component(1).size();
  gsVector<> coefvec0(len0);
  // Uniformly distribute control points per component
  coefvec0.setLinSpaced(len0,0.0,L);
  gsVector<> coefvec1(basis.component(1).size());
  coefvec1.setLinSpaced(len1,0.0,B);

  // Z coordinate is zero
  coefs.col(2).setZero();

  // Define a matrix with ones
  gsVector<> temp(len0);
  temp.setOnes();
  for (index_t k = 0; k < len1; k++)
  {
    // First column contains x-coordinates (length)
    coefs.col(1).segment(k*len0,len0) = coefvec0;
    // Second column contains z-coordinates (height)
    coefs.col(0).segment(k*len0,len0) = temp*coefvec1.at(k);
  }
  // Create gsGeometry-derived object for the patch
  gsTensorBSpline<2,real_t> shape(basis,coefs);

  gsMultiPatch<T> mp;
  mp.addPatch(shape);
  mp.addAutoBoundaries();

  return mp;
}

template <class T>
gsMultiPatch<T> FrustrumDomain(int n, int p, T R1, T R2, T h)
{
  // -------------------------------------------------------------------------
  // --------------------------Make beam geometry-----------------------------
  // -------------------------------------------------------------------------
  // n = number of uniform refinements over the height; n = 0, only top and bottom part

  int dim = 3; //physical dimension
  gsKnotVector<> kv0;
  kv0.initUniform(0,1,0,3,1);
  gsKnotVector<> kv1;
  kv1.initUniform(0,1,0,3,1);

  // Refine n times
  for(index_t i = 0; i< n; ++i)
      kv1.uniformRefine();

  gsDebug<<kv1;

  // Make basis
  // gsTensorNurbsBasis<2,T> basis(kv0,kv1);

  // Initiate coefficient matrix
  index_t N = math::pow(2,n)+2;
  gsMatrix<> coefs(3*N,dim);
  gsMatrix<> tmp(3,3);
  T R,H;

  gsMatrix<> weights(3*N,1);
  for (index_t k=0; k!= N; k++)
  {
    R = k*(R2-R1)/(N-1) + R1;
    H = k*h/(N-1);
    tmp<< R,0,H,
          R,R,H,
          0,R,H;

    coefs.block(3*k,0,3,3) = tmp;

    weights.block(3*k,0,3,1) << 1,0.70711,1;
  }

  // Create gsGeometry-derived object for the patch
  gsTensorNurbs<2,real_t> shape(kv0,kv1,coefs,weights);

  gsMultiPatch<T> mp;
  mp.addPatch(shape);
  mp.addAutoBoundaries();

  // Elevate up to order p
  if (p>2)
  {
    for(index_t i = 2; i< p; ++i)
        mp.patch(0).degreeElevate();    // Elevate the degree
  }

  // // Refine n times
  // for(index_t i = 0; i< n; ++i)
  //     mp.patch(0).uniformRefine();

  return mp;
}


/*
    to do:
    =  make function for construction of the solution given the space and the mp
*/



/*
template<class T>
void gsShellAssembler<T>::applyLoads()
{
    gsMatrix<T>        bVals;
    gsMatrix<unsigned> acts,globalActs;

    for (size_t i = 0; i< m_pLoads.numLoads(); ++i )
    {
        if ( m_pLoads[i].parametric )
        {
            m_bases.front().basis(m_pLoads[i].patch).active_into( m_pLoads[i].point, acts );
            m_bases.front().basis(m_pLoads[i].patch).eval_into  ( m_pLoads[i].point, bVals);
        }
        else
        {
            gsMatrix<> forcePoint;
            m_patches.patch(m_pLoads[i].patch).invertPoints(m_pLoads[i].point,forcePoint);
            u.source().piece(m_pLoads[i].patch).active_into( forcePoint, acts );
            u.source().piece(m_pLoads[i].patch).active_into( forcePoint, bVals);
        }

        // translate patch-local indices to global dof indices
        for (size_t j = 0; j< 3; ++j)
        {
            if (m_pLoads[i].value[j] != 0.0)
            {
                u.dofMappers[j].localToGlobal(acts, m_pLoads[i].patch, globalActs);

                for (index_t k=0; k < globalActs.rows(); ++k)
                {
                    if (int(globalActs(k,0)) < m_dofs)
                        m_rhs(globalActs(k,0), 0) += bVals(k,0) * m_pLoads[i].value[j];
                }
            }
        }
    }
}
*/
