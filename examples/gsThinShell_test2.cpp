/** @file gsThinShell_test2.cpp

    @brief Example testing and debugging thin shell solver. Based on gsThinShell_test

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst
*/

#include <gismo.h>
#include <gismo_dev.h>

#include <gsThinShell2/gsThinShellAssembler.h>
#include <gsThinShell2/gsMaterialMatrix.h>

//#include <gsThinShell/gsNewtonIterator.h>

using namespace gismo;

template <class T>
gsMultiPatch<T> RectangularDomain(int n, int m, int p, int q, T L, T B);
template <class T>
gsMultiPatch<T> RectangularDomain(int n, int p, T L, T B);

template <class T>
gsMultiPatch<T> RectangularDomainVert(int n, int m, int p, int q, T L, T B);
template <class T>
gsMultiPatch<T> RectangularDomainVert(int n, int p, T L, T B);

template <class T>
gsMultiPatch<T> RectangularDomain90(int n, int m, int p, int q, T L, T B);
template <class T>
gsMultiPatch<T> RectangularDomain90(int n, int p, T L, T B);


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
    bool nonlinear = false;
    bool verbose = false;
    std::string fn;

    real_t E_modulus = 1.0;
    real_t PoissonRatio = 0.0;
    real_t Density = 1.0;
    real_t thickness = 1.0;

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving",  numRefine );
    cmd.addInt( "t", "testCase", "Test case to run: 1 = unit square; 2 = Scordelis Lo Roof",  testCase );
    cmd.addInt( "m", "Material", "Material law",  material );
    cmd.addInt( "c", "Compressibility", "1: compressible, 0: incompressible",  Compressibility );
    cmd.addString( "f", "file", "Input XML file", fn );
    cmd.addSwitch("nl", "Solve nonlinear problem", nonlinear);
    cmd.addSwitch("verbose", "Full matrix and vector output", verbose);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("stress", "Create a ParaView visualization file with the stresses", stress);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Read input file]
    gsMultiPatch<> mp;
    gsMultiPatch<> mp_def;
    if (testCase == 1 )
    {
        thickness = 0.25;
        E_modulus = 4.32E8;
        fn = "../extensions/unsupported/filedata/scordelis_lo_roof.xml";
        gsReadFile<>(fn, mp);
        PoissonRatio = 0.0;
    }
    else if (testCase == 2)
    {
        thickness = 0.04;
        E_modulus = 6.825E7;
        PoissonRatio = 0.3;
        gsReadFile<>("quarter_hemisphere.xml", mp);
    }
    else if (testCase == 3)
    {
        thickness = 3;
        E_modulus = 3E6;
        PoissonRatio = 0.3;
        gsReadFile<>("pinched_cylinder.xml", mp);
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
        E_modulus = (2+PoissonRatio)*mu;
        gsReadFile<>("quarter_sphere.xml", mp);
    }
    else if (testCase == 7)
    {
        thickness = 0.1;
        real_t mu = 4.225e5;
        PoissonRatio = 0.3;
        E_modulus = (2+PoissonRatio)*mu;
        gsReadFile<>("quarter_frustrum.xml", mp);
    }
    else if (testCase == 17)
    {
        // Unit square
        mp = RectangularDomain(0,2,1.0,1.0);
        E_modulus = 4.497000000e6;
        thickness = 0.001;
        PoissonRatio = 0.4999;
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
    }
    //! [Read input file]

    // p-refine
    if (numElevate!=0)
        mp.degreeElevate(numElevate);

    // h-refine
    for (int r =0; r < numRefine; ++r)
        mp.uniformRefine();

    mp_def = mp;

    //! [Refinement]
    gsMultiBasis<> dbasis(mp);

    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
    gsInfo << dbasis.basis(0)<<"\n";

    gsBoundaryConditions<> bc;
    bc.setGeoMap(mp);
    gsVector<> tmp(3);
    tmp << 0, 0, 0;

    gsVector<> neu(3);
    neu << 0, 0, 0;

    gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();

    gsFunctionExpr<> displx("1.0",3);
    gsFunctionExpr<> disply("0.5",3);

    gsConstantFunction<> neuData(neu,3);
    real_t pressure = 0.0;
    if (testCase == 0)
    {
        for (index_t i=0; i!=3; ++i)
        {
            bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0 ,false,i);
            bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0 ,false,i);
            bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0 ,false,i);
            bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0 ,false,i);
        }

        // bc.addCondition(boundary::north, condition_type::clamped, 0, 0 ,false,2);
        // bc.addCondition(boundary::east, condition_type::clamped, 0, 0 ,false,2);
        // bc.addCondition(boundary::south, condition_type::clamped, 0, 0 ,false,2);
        // bc.addCondition(boundary::west, condition_type::clamped, 0, 0 ,false,2);

        // tmp << 0,0,0;
        tmp << 0,0,-1;

        // Point loads
        // gsVector<> point(2);
        // gsVector<> load (3);
        // point<< 0.5, 0.5 ; load << 0.0, 1.0, 0.0 ;
        // pLoads.addLoad(point, load, 0 );
    }
    else if (testCase == 1)
    {
        // Diaphragm conditions
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0 ,false, 1 ); // unknown 1 - y
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0 ,false, 2 ); // unknown 2 - z

        bc.addCornerValue(boundary::southwest, 0.0, 0, 0); // (corner,value, patch, unknown)

        bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0 ,false, 1 ); // unknown 1 - y
        bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0 ,false, 2 ); // unknown 2 - z

        // Surface forces
        tmp << 0, 0, -90;
    }
    else if (testCase == 2)
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

        // Surface forces
        tmp.setZero();

        // Point loads
        gsVector<> point(2); point<< 1.0, 1.0 ;
        gsVector<> load (3); load << 0.0, 0.0, -0.25 ;
        pLoads.addLoad(point, load, 0 );
    }
    else if (testCase == 4)
    {
        for (index_t i = 0; i!=3; i++)
        {
            bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, i );
            bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, i );
            bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, i );
            bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, i );
        }
        pressure = -1.0;
    }
    else if (testCase == 5)
    {
        for (index_t i = 0; i!=3; i++)
        {
            bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, i );
            bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, i );
            bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, i );
            bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, i );
        }

        // Point loads
        gsVector<> point(2); point<< 0.5,0.5 ;
        gsVector<> load (3); load << 0.0, 0.0, -1 ;
        pLoads.addLoad(point, load, 0 );
    }
    else if (testCase == 6) // balloon
    {
        // bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        // bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        // bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

        // bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        // bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

        // Symmetry in x-direction:
        bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 );
        bc.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 1 );
        bc.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 2 );

        // Symmetry in y-direction:
        bc.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 0 );
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 );
        bc.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 2 );

        // Pressure
        pressure = 5e3;
    }
    else if (testCase == 7)
    {
        neu << 0, 0, -100.0;
        neuData.setValue(neu,3);

        bc.addCondition(boundary::north, condition_type::neumann, &neuData );
        // bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 2 - z
        // bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 2 - z

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
    else if (testCase == 10)
    {
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
        bc.addCondition(boundary::east, condition_type::dirichlet, &displx, 0, false,  0 ); // unknown 0 - x
    }
    else if (testCase == 11)
    {
        for (index_t i=0; i!=3; ++i)
        {
            bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, i ); // unknown 1 - y
            bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, i ); // unknown 2 - z
        }

        bc.addCondition(boundary::east, condition_type::clamped, 0, 0 ,false,2);
        bc.addCondition(boundary::west, condition_type::clamped, 0, 0 ,false,2);

        tmp<<0,0,-1;
    }
    else if (testCase == 12)
    {
        for (index_t i=0; i!=3; ++i)
        {
            bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, i ); // unknown 2 - z
        }

        bc.addCondition(boundary::west, condition_type::clamped, 0, 0 ,false,2);

        neu << 0, 0, -0.1;
        neuData.setValue(neu,3);

        bc.addCondition(boundary::east, condition_type::neumann, &neuData );
    }
    else if (testCase == 13)
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

        tmp << 0,0,-1;
    }
    else if (testCase == 14)
    {
        for (index_t i=0; i!=3; ++i)
        {
            bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, i ); // unknown 2 - z
        }
        bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
        bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

        tmp << 0,0,-1;

        bc.addCondition(boundary::east, condition_type::collapsed, 0, 0 ,false,0);
        bc.addCondition(boundary::east, condition_type::collapsed, 0, 0, false, 2 );
    }
    else if (testCase == 14)
    {
        for (index_t i=0; i!=3; ++i)
        {
            bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, i ); // unknown 2 - z
        }
        bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
        bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

        tmp << 0,0,-1;

        bc.addCondition(boundary::east, condition_type::collapsed, 0, 0 ,false,0);
        bc.addCondition(boundary::east, condition_type::collapsed, 0, 0, false, 2 );
    }
    else if (testCase == 15)
    {
        for (index_t i=0; i!=3; ++i)
        {
            bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, i ); // unknown 2 - z
        }
        bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
        bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

        bc.addCondition(boundary::east, condition_type::collapsed, 0, 0 ,false,0);

        gsVector<> point(2); point<< 1.0, 0.5 ;
        gsVector<> load (3); load << 0.1, 0.0, 0.0 ;
        pLoads.addLoad(point, load, 0 );
    }
    else if (testCase == 16)
    {
      bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
      bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
      bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

      bc.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 2 );

        real_t Load = 1e-2;
        gsVector<> point(2);
        gsVector<> load (3);
        point<< 1.0, 0.5 ;
        load << 0.0, 0.0, Load ;
        pLoads.addLoad(point, load, 0 );
    }
    else if (testCase == 17)
    {
        real_t Load = 1e-1;
        neu << -Load, 0, 0;
        neuData.setValue(neu,3);

        bc.addCondition(boundary::west, condition_type::neumann, &neuData ); // unknown 0 - x
        bc.addCondition(boundary::west, condition_type::collapsed, 0, 0, false, 0 ); // unknown 0 - x
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

        bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
    }
    //! [Refinement]

    // Linear isotropic material model
    gsConstantFunction<> force(tmp,3);
    gsConstantFunction<> pressFun(pressure,3);
    gsFunctionExpr<> t(std::to_string(thickness), 3);
    gsFunctionExpr<> E(std::to_string(E_modulus),3);
    gsFunctionExpr<> nu(std::to_string(PoissonRatio),3);
    gsFunctionExpr<> rho(std::to_string(Density),3);
    gsConstantFunction<> ratio(7.0,3);

    real_t mu = E_modulus / (2 * (1 + PoissonRatio));

    gsDebugVar(mu);
    gsConstantFunction<> alpha1(2.0,3);
    gsConstantFunction<> mu1(7.0*mu/8.0,3);
    gsConstantFunction<> alpha2(-2.0,3);
    gsConstantFunction<> mu2(-mu/8.0,3);
    // gsConstantFunction<> alpha3()
    // gsConstantFunction<> mu3()

    // gsMaterialMatrix materialMatrixNonlinear(mp,mp_def,t,E,nu,rho);
    std::vector<gsFunction<>*> parameters(3);
    parameters[0] = &E;
    parameters[1] = &nu;
    parameters[2] = &ratio;
    gsMaterialMatrix materialMatrixNonlinear(mp,mp_def,t,parameters,rho);

    std::vector<gsFunction<>*> parameters2(6);
    if (material==14)
    {
        parameters2[0] = &E;
        parameters2[1] = &nu;
        parameters2[2] = &mu1;
        parameters2[3] = &alpha1;

        parameters2[4] = &mu2;
        parameters2[5] = &alpha2;

        // parameters[6] = ;
        // parameters[7] = ;
        materialMatrixNonlinear.setParameters(parameters2);
    }


    materialMatrixNonlinear.options().setInt("MaterialLaw",material);
    materialMatrixNonlinear.options().setInt("Compressibility",Compressibility);


    // Linear anisotropic material model
    real_t pi = math::atan(1)*4;
    std::vector<std::pair<real_t,real_t>> Evec;
    std::vector<std::pair<real_t,real_t>> nuvec;
    std::vector<real_t> Gvec;
    std::vector<real_t> tvec;
    std::vector<real_t> phivec;

    index_t kmax = 2;
    for (index_t k=0; k != kmax; ++k)
    {
        Evec.push_back( std::make_pair(E_modulus,E_modulus) );
        nuvec.push_back( std::make_pair(PoissonRatio,PoissonRatio) );
        Gvec.push_back( 0.5 * E_modulus / (1+PoissonRatio) );
        tvec.push_back( thickness/kmax );
        phivec.push_back( k / kmax * pi/2.0);
    }
    gsMaterialMatrix materialMatrixComposite(mp,mp_def,tvec,Evec,Gvec,nuvec,phivec);

    gsThinShellAssembler assembler(mp,dbasis,bc,force,materialMatrixNonlinear);
    assembler.setPointLoads(pLoads);
    if (pressure!= 0.0)
        assembler.setPressure(pressFun);


    gsSparseSolver<>::CGDiagonal solver;

    // gsInfo<< A.numDofs() <<"\n"<<std::flush;

    // assembler.assembleMass();

    gsVector<> pt(2); pt.setConstant(0.25);
    // evaluateFunction(ev, u * ff * meas(G), pt); // evaluates an expression on a point

    // ! [Solve linear problem]

    // assemble system
    assembler.assemble();
    // solve system
    solver.compute( assembler.matrix() );

    if (verbose)
    {
        gsInfo<<"Matrix: \n"<<assembler.matrix().toDense()<<"\n";
        gsInfo<<"Vector: \n"<<assembler.rhs().transpose()<<"\n";
    }


    // // TEST MATRIX INTEGRATION
    gsMaterialMatrix materialMatrixNonlinear2(mp,mp_def,t,parameters,rho);
    if (material==14)
        materialMatrixNonlinear2.setParameters(parameters2);
    materialMatrixNonlinear2.options().setInt("MaterialLaw",material);
    materialMatrixNonlinear2.options().setInt("Compressibility",Compressibility);

    materialMatrixNonlinear2.info();

    gsMaterialMatrix materialMatrixTest(mp,mp_def,t,parameters,rho);
    materialMatrixTest.options().setInt("MaterialLaw",13);
    materialMatrixTest.options().setInt("Compressibility",Compressibility);
    gsVector<> testPt(2);
    testPt<<0.351135,0.85235;
    // testPt.setConstant(0.25);
    gsMatrix<> testResult1, testResult2;
    // materialMatrixTest.makeVector(0);
    materialMatrixTest.makeMatrix(0);
    materialMatrixTest.eval_into(testPt,testResult1);
    gsDebugVar(testResult1);

    // materialMatrixNonlinear2.makeVector(0);
    materialMatrixNonlinear2.makeMatrix(0);
    materialMatrixNonlinear2.eval_into(testPt,testResult2);
    gsDebugVar(testResult2);
    gsDebugVar(testResult1-testResult2);

    materialMatrixTest.makeVector(0);
    materialMatrixTest.eval_into(testPt,testResult1);
    gsDebugVar(testResult1);

    materialMatrixNonlinear2.makeVector(0);
    materialMatrixNonlinear2.eval_into(testPt,testResult2);
    gsDebugVar(testResult2);
    gsDebugVar(testResult1-testResult2);
    // // ! TEST MATRIX INTEGRATION

    gsVector<> solVector = solver.solve(assembler.rhs());

    mp_def = assembler.constructSolution(solVector);

    // mp_def = RectangularDomain(0,2,2.0,1.0);
    //     // p-refine
    // if (numElevate!=0)
    //     mp_def.degreeElevate(numElevate);

    // // h-refine
    // for (int r =0; r < numRefine; ++r)
    //     mp_def.uniformRefine();

    gsWriteParaview<>( mp    , "mp", 1000, true);
    gsWriteParaview<>( mp_def, "mp_def", 1000, true);


    if ((plot) && (!nonlinear))
    {

        gsMultiPatch<> deformation = mp_def;
        for (size_t k = 0; k != mp_def.nPatches(); ++k)
            deformation.patch(0).coefs() -= mp.patch(0).coefs();

        gsField<> solField(mp, deformation);
        gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview<>( solField, "solution", 1000, true);
    }
    if ((stress) && (!nonlinear))
    {
        gsPiecewiseFunction<> stresses;
        assembler.constructStress(mp_def,stresses,stress_type::principal_stretch);
        gsField<> stressField(mp,stresses, true);

        gsWriteParaview( stressField, "stress", 5000);

        gsMatrix<> stretch;
        stretch=assembler.computePrincipalStretches(pt,mp_def,0.0);
        gsDebugVar(stretch);
    }

    /*Something with Dirichlet homogenization*/

    // ! [Solve linear problem]

    // ! [Solve nonlinear problem]

    // assembler.assemble();
    // assembler.assemble(mp_def);

    real_t residual = assembler.rhs().norm();
    real_t residual0 = residual;
    real_t residualOld = residual;
    gsVector<> updateVector = solVector;
    if (nonlinear)
    {
        index_t itMax = 100;
        real_t tol = 1e-8;
        for (index_t it = 0; it != itMax; ++it)
        {
            assembler.assemble(mp_def);
            // solve system
            solver.compute( assembler.matrix() );

            if (verbose)
            {
                gsInfo<<"Matrix: \n"<<assembler.matrix().toDense()<<"\n";
                gsInfo<<"Vector: \n"<<assembler.rhs().transpose()<<"\n";
            }


            updateVector = solver.solve(assembler.rhs()); // this is the UPDATE
            residual = assembler.rhs().norm();

            solVector += updateVector;
            mp_def = assembler.constructSolution(solVector);


            gsInfo<<"Iteration: "<< it
                   <<", residue: "<< residual
                   <<", update norm: "<<updateVector.norm()
                   <<", log(Ri/R0): "<< math::log10(residualOld/residual0)
                   <<", log(Ri+1/R0): "<< math::log10(residual/residual0)
                   <<"\n";

            residualOld = residual;

            if (updateVector.norm() < tol)
                break;

            // ADD DIRICHLET HOMOGENIZE
        }
    }

    // TEST MATRIX INTEGRATION
    materialMatrixNonlinear2.options().setInt("MaterialLaw",material);
    materialMatrixNonlinear2.options().setInt("Compressibility",Compressibility);

    materialMatrixTest.options().setInt("MaterialLaw",13);
    materialMatrixTest.options().setInt("Compressibility",Compressibility);

    // materialMatrixTest.makeVector(0);
    materialMatrixTest.makeMatrix(0);
    materialMatrixTest.eval_into(testPt,testResult1);
    gsDebugVar(testResult1);

    // materialMatrixNonlinear2.makeVector(0);
    materialMatrixNonlinear2.makeMatrix(0);
    materialMatrixNonlinear2.eval_into(testPt,testResult2);
    gsDebugVar(testResult2);
    gsDebugVar(testResult1-testResult2);

    materialMatrixTest.makeVector(0);
    materialMatrixTest.eval_into(testPt,testResult1);
    gsDebugVar(testResult1);

    materialMatrixNonlinear2.makeVector(0);
    materialMatrixNonlinear2.eval_into(testPt,testResult2);
    gsDebugVar(testResult2);
    gsDebugVar(testResult1-testResult2);

    // // ! TEST MATRIX INTEGRATION

    // gsDebugVar(assembler.computePrincipalStretches(pts,mp_def));


    // ! [Solve nonlinear problem]

    // // For Neumann (same for Dirichlet/Nitche) conditions
    // variable g_N = assembler.getBdrFunction();
    // assembler.assembleRhsBc(u * g_N.val() * nv(G).norm(), bc.neumannSides() );

    // // Penalize the matrix? (we need values for the DoFs to be enforced..
    // // function/call:  penalize_matrix(DoF_indices, DoF_values)
    // // otherwise: should we tag the DoFs inside "u" ?

    // gsInfo<<"RHS rows = "<<assembler.rhs().rows()<<"\n";
    // gsInfo<<"RHS cols = "<<assembler.rhs().cols()<<"\n";
    // gsInfo<<"MAT rows = "<<assembler.matrix().rows()<<"\n";
    // gsInfo<<"MAT cols = "<<assembler.matrix().cols()<<"\n";

    // gsInfo<< assembler.rhs().transpose() <<"\n";
    // gsInfo<< assembler.matrix().toDense()<<"\n";

    // // ADD BOUNDARY CONDITIONS! (clamped will be tricky..............)

    mp_def = assembler.constructSolution(solVector);
    gsMultiPatch<> deformation = mp_def;
    for (size_t k = 0; k != mp_def.nPatches(); ++k)
        deformation.patch(0).coefs() -= mp.patch(0).coefs();

    // ! [Export visualization in ParaView]
    if ( (plot) && (nonlinear) )
    {

        gsField<> solField(mp, deformation);
        gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview<>( solField, "solution", 1000, true);
        // ev.options().setSwitch("plot.elements", true);
        // ev.writeParaview( u_sol   , G, "solution");

        // gsFileManager::open("solution.pvd");
    }
    if ((stress) && (nonlinear))
    {
        gsPiecewiseFunction<> stresses;
        assembler.constructStress(mp_def,stresses,stress_type::principal_stretch);
        gsField<> stressField(mp,stresses, true);

        gsWriteParaview( stressField, "stress", 5000);
    }

    gsInfo <<"Maximum deformation coef: "
           << deformation.patch(0).coefs().colwise().maxCoeff() <<".\n";
    gsInfo <<"Minimum deformation coef: "
           << deformation.patch(0).coefs().colwise().minCoeff() <<".\n";


    gsInfo <<"Area (undeformed) = "<<assembler.getArea(mp)<<"\tArea (deformed) = "<<assembler.getArea(mp_def)<<"\n";

    return EXIT_SUCCESS;

}// end main

template <class T>
void evaluateFunction(gsExprEvaluator<T> ev, auto expression, gsVector<T> pt)
{
    gsMatrix<T> evresult = ev.eval( expression,pt );
    gsInfo << "Eval on point ("<<pt.at(0)<<" , "<<pt.at(1)<<") :\n"<< evresult;
    gsInfo << "\nEnd ("<< evresult.rows()<< " x "<<evresult.cols()<<")\n";
};

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
