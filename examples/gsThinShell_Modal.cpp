/** @file gsThinShell_Modal.cpp

    @brief Example to compute eigenvalues and eigenmodes of a vibrating shell

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#include <gismo.h>

#ifdef gsKLShell_ENABLED
#include <gsKLShell/src/gsThinShellAssembler.h>
#include <gsKLShell/src/getMaterialMatrix.h>
#endif

#include <gsStructuralAnalysis/src/gsEigenSolvers/gsModalSolver.h>

using namespace gismo;

void writeToCSVfile(std::string name, gsMatrix<> matrix)
{
    std::ofstream file(name.c_str());
    for(int  i = 0; i < matrix.rows(); i++){
        for(int j = 0; j < matrix.cols(); j++){
           std::string str = std::to_string(matrix(i,j));
           if(j+1 == matrix.cols()){
               file<<std::setprecision(10)<<str;
           }else{
               file<<std::setprecision(10)<<str<<',';
           }
        }
        file<<'\n';
    }
  }

template <class T>
gsMultiPatch<T> Rectangle(T L, T B);

#ifdef gsKLShell_ENABLED
int main (int argc, char** argv)
{
    // Input options
    int numElevate  = 1;
    int numHref     = 1;
    int numKref     = 1;
    bool plot       = false;
    bool sparse     = false;
    bool nonlinear  = false;
    bool first  = false;
    int mode = 0;

    std::string assemberOptionsFile("options/solver_options.xml");

    real_t thickness     = 1;
    real_t width = 1; // Width of the strip is equal to 1.
    real_t length = 10; // Length of the strip is equal to 10.
    real_t Area = thickness*width;

    real_t E_modulus     = 1e0;
    real_t PoissonRatio = 0;
    real_t Density = 1e0;
    gsMultiPatch<> mp;

    real_t shift = 0.0;

    int testCase = 0;

    int result = 0;

    bool write = false;

    gsCmdLine cmd("Modal analysis for thin shells.");
    cmd.addString( "f", "file", "Input XML file for assembler options", assemberOptionsFile );
    cmd.addInt("r","hRefine",
               "Number of dyadic h-refinement (bisection) steps to perform before solving",
               numHref);
    cmd.addInt("t", "testcase",
                "Test case: 0: Beam - pinned-pinned, 1: Beam - fixed-fixed, 2: beam - fixed-free, 3: plate - fully pinned, 4: plate - fully fixed, 5: circle - fully pinned, 6: 5: circle - fully fixed",
               testCase);
    cmd.addInt("m", "nmode",
               "Mode shape number, starting from 0",
              mode);
    cmd.addInt("e","degreeElevation",
               "Number of degree elevation steps to perform on the Geometry's basis before solving",
               numElevate);
    cmd.addInt("k","continuityDecrease",
               "++",
               numKref);
    cmd.addReal("s","shift", "eigenvalue shift", shift);
    cmd.addSwitch("nl", "Nonlinear elasticity (otherwise linear)", nonlinear);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addSwitch("first", "Plot only first", first);
    cmd.addSwitch("write", "Write convergence data to file", write);
    cmd.addSwitch("sparse", "Use sparse solver", sparse);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }


    gsFileData<> fd(assemberOptionsFile);
    gsOptionList opts;
    fd.getFirst<gsOptionList>(opts);

    std::string fn;

    real_t EA(0),EI(0),r(0);
    if (testCase==0 || testCase==1 || testCase==2)
    {
        if (testCase==2)
        {
            width = 0.5;
            thickness = 0.5;
            E_modulus = 210e9;
            Density = 7800;
        }
        
        mp.addPatch(gsNurbsCreator<>::BSplineRectangle(0,0,length,width));
        mp.addAutoBoundaries();
        mp.embed(3);
        EA = E_modulus*Area;
        EI = 1.0/12.0*(width*math::pow(thickness,3))*E_modulus;
        r = math::sqrt(EI/EA);
        gsInfo<<"EI = "<<EI<<"; EA = "<<EA<<"; r = "<<r<<"\n";
    }
    else if (testCase==3)
    {
        mp.addPatch(gsNurbsCreator<>::BSplineRectangle(0,0,4.0,4.0));
        mp.addAutoBoundaries();
        mp.embed(3);
        E_modulus = 1e5;
        PoissonRatio = 0.3;
        Density = 1e0;
        thickness = 1e-3;
    }
    else if (testCase==4 || testCase==5)
    {
      mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // degree
      mp.addAutoBoundaries();
      mp.embed(3);
    }

    else if (testCase==6 || testCase==7)
    {
      fn = "planar/unitcircle.xml";
      gsReadFile<>(fn, mp);
    }
    else if (testCase==8)
    {
        mp.addPatch( gsNurbsCreator<>::BSplineTriangle(1,1) ); // degree
        mp.addAutoBoundaries();
        mp.embed(3);
    }
    else if (testCase==9)
    {
        std::string fn = "planar/weirdShape.xml";
        gsReadFile<>(fn, mp);
        mp.addAutoBoundaries();

        gsKnotVector<> kv1(0,1,12,3);
        gsKnotVector<> kv2(0,1,0,3);
        gsTensorBSplineBasis<2> basis(kv2,kv1);
        gsMatrix<> coefs;
        gsQuasiInterpolate<real_t>::localIntpl(basis, mp.patch(0), coefs);
        gsTensorBSpline<2> bspline(kv2,kv1,coefs);

        gsDebugVar(mp.patch(0).coefs());
        gsDebugVar(coefs);
        gsDebugVar(bspline);

        mp.patch(0) = bspline;
        mp.embed(3);
        mp.degreeElevate(1);

        thickness = 0.01;
        PoissonRatio = 0.3;
        E_modulus     = 1e0;
        Density = 1e0;
    }

    for(index_t i = 0; i< numElevate; ++i)
        mp.patch(0).degreeElevate();    // Elevate the degree

    // h-refine
    for(index_t i = 0; i< numHref; ++i)
        mp.patch(0).uniformRefine();

    for(index_t i = 0; i< numKref; ++i)
        mp.patch(0).degreeElevate();    // Elevate the degree

    gsMultiBasis<> dbasis(mp);


    gsInfo<<"Basis (patch 0): "<< mp.patch(0).basis() << "\n";

    gsWriteParaview<>(mp, "mp", 500,true,true);

    // Boundary conditions
    std::vector< std::pair<patchSide,int> > clamped;
    gsBoundaryConditions<> BCs;
    BCs.setGeoMap(mp);
    gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();
    gsPointLoads<real_t> pMass = gsPointLoads<real_t>();

    // Initiate Surface forces
    std::string tx("0");
    std::string ty("0");
    std::string tz("0");

    gsVector<> tmp(3);
    tmp << 0, 0, 0;

    std::vector<real_t> omegas;

    if (testCase == 0)
    {
        // Beam
        // Pinned-Pinned [no loads]
            // Left
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
            // Right
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

        // Restrain sides
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y

        for (index_t n=0; n!=10; n++)
          omegas.push_back(math::pow((n*3.1415926535/length),2)*math::pow(EI/(Density*Area),0.5));

    }
    else if (testCase == 1)
    {
        // Beam
        // Fixed-Fixed
            // Left
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
        BCs.addCondition(boundary::west, condition_type::clamped,0,0,false,2);
            // Right
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
        BCs.addCondition(boundary::east, condition_type::clamped,0,0,false,2);

        // Restrain sides
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y

        omegas.push_back(math::pow((4.73004/length),2)*math::pow(EI/(Density*Area),0.5));

    }
    else if (testCase == 2)
    {
        // Beam
        // Fixed-Free
            // Left
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
        BCs.addCondition(boundary::west, condition_type::clamped,0,0,false,2);

        // Restrain sides
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        // Restrain tip in horizontal and side direction
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y

        omegas.push_back(3.52*math::pow((1/length),2)*math::pow(EI/(Density*Area),0.5));

        gsVector<> point(2);
        gsVector<> mass (1);

        point<< 1.0, 0.5 ; mass << 100.0 ;
        pMass.addLoad(point, mass, 0 );
    }
    else if (testCase == 3)
    {
        // Plate
        // Free-Free-Free-Free
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
            // Right
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
            // Top
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
            // Bottom
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        omegas.push_back(0.0);
    }
    else if (testCase == 4)
    {
        thickness = 0.01;
        E_modulus = 1e0;
        Density = 1e0;
        PoissonRatio = 0.3;
        // Plate
        // Pinned-Pinned-Pinned-Pinned
            // Left
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
            // Right
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
            // Top
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
            // Bottom
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z


        real_t D = E_modulus*math::pow(thickness,3)/(12*(1-math::pow(PoissonRatio,2)));
        gsInfo<<"D = "<<D<<"\n";
        for (index_t m=1; m!=10; m++)
          for (index_t n=1; n!=10; n++)
            omegas.push_back((math::pow(m/1.0,2)+math::pow(n/1.0,2))*math::pow(3.1415926535,2)*math::sqrt(D / (Density * thickness)));

        std::sort(omegas.begin(),omegas.end());
    }
    else if (testCase == 5)
    {
        thickness = 0.01;
        PoissonRatio = 0.3;
        E_modulus = 1e0;
        Density = 1e0;
        // Plate
        // Clamped-Clamped-Clamped-Clamped
            // Left
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
        BCs.addCondition(boundary::west, condition_type::clamped,0,0,false,2);

            // Right
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
        BCs.addCondition(boundary::east, condition_type::clamped,0,0,false,2);

            // Top
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
        BCs.addCondition(boundary::north, condition_type::clamped,0,0,false,2);

            // Bottom
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
        BCs.addCondition(boundary::south, condition_type::clamped,0,0,false,2);


        // real_t D = E_modulus*math::pow(thickness,3)/(12*(1-math::pow(PoissonRatio,2)));
        // gsInfo<<"D = "<<D<<"\n";
        // for (index_t m=0; m!=10; m++)
        //   for (index_t n=0; n!=10; n++)
        //     omegas.push_back((math::pow(m/1.0,2)+math::pow(n/1.0,2))*math::pow(3.1415926535,2)*math::sqrt(D / (Density * thickness)));

        // std::sort(omegas.begin(),omegas.end());

        omegas.push_back(0);
    }
    else if (testCase == 6)
    {
        thickness = 0.01;
        PoissonRatio = 0.3;
        // Circle
        // Pinned-Pinned-Pinned-Pinned
            // Left
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

            // Right
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

            // Top
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

            // Bottom
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z


        omegas.push_back(0);
    }
    else if (testCase == 7)
    {
        thickness = 0.01;
        PoissonRatio = 0.3;
        // Circle
        // Pinned-Pinned-Pinned-Pinned
            // Left
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
        BCs.addCondition(boundary::west, condition_type::clamped,0,0,false,2);

            // Right
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
        BCs.addCondition(boundary::east, condition_type::clamped,0,0,false,2);

            // Top
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
        BCs.addCondition(boundary::north, condition_type::clamped,0,0,false,2);

            // Bottom
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
        BCs.addCondition(boundary::south, condition_type::clamped,0,0,false,2);

        real_t D = E_modulus*math::pow(thickness,3)/(12*(1-math::pow(PoissonRatio,2)));
        gsInfo<<"D = "<<D<<"\n";
        gsVector<> gammas(10);
        gammas<<3.1962206158252, 4.6108999, 5.9056782, 6.3064370, 7.1435310, 7.7992738, 8.3466059, 9.1968826, 9.4394991, 9.5257014;
        for (index_t n=0; n!=gammas.size(); n++)
          omegas.push_back(math::pow(math::pow(gammas[n],4)*D/(Density*thickness),0.5));
    }
    else if (testCase == 8)
    {
        thickness = 0.01;
        PoissonRatio = 0.3;
        // Circle
        // Pinned-Pinned-Pinned-Pinned
        for (index_t k = 0; k!=3; ++k)
        {
            BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, k ); // unknown 0 - x
            BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, k ); // unknown 0 - x
            BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, k ); // unknown 0 - x
            BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, k ); // unknown 0 - x
        }

        BCs.addCondition(boundary::west, condition_type::clamped,0,0,false,2);
        BCs.addCondition(boundary::east, condition_type::clamped,0,0,false,2);
        BCs.addCondition(boundary::north, condition_type::clamped,0,0,false,2);
        BCs.addCondition(boundary::south, condition_type::clamped,0,0,false,2);

        omegas.push_back(0);
    }
    else if (testCase==9)
    {
        for (index_t i = 0; i != 3; ++i)
        {
            BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, i);
        }
        BCs.addCondition(boundary::south, condition_type::clamped,0,0,false,2);

        omegas.push_back(0);
    }
    else
        GISMO_ERROR("TESTCASE UNKNOWN!");

    gsFunctionExpr<> surfForce(tx,ty,tz,3);
    // Initialise solution object
    gsMultiPatch<> mp_def = mp;

    // Linear isotropic material model
    gsConstantFunction<> force(tmp,3);
    gsFunctionExpr<> t(std::to_string(thickness), 3);
    gsFunctionExpr<> E(std::to_string(E_modulus),3);
    gsFunctionExpr<> nu(std::to_string(PoissonRatio),3);
    gsFunctionExpr<> rho(std::to_string(Density),3);

    std::vector<gsFunctionSet<>*> parameters(2);
    parameters[0] = &E;
    parameters[1] = &nu;

    gsMaterialMatrixBase<real_t>* materialMatrix;

    gsOptionList options;
    options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
    options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
    materialMatrix = getMaterialMatrix<3,real_t>(mp,t,parameters,rho,options);

    gsThinShellAssemblerBase<real_t>* assembler;
    assembler = new gsThinShellAssembler<3, real_t, true >(mp,dbasis,BCs,force,materialMatrix);

    assembler->setPointMass(pMass);
    assembler->setOptions(opts);
    // Initialise solution object
    gsMultiPatch<> solution = mp;

    assembler->assemble();
    gsSparseMatrix<> K = assembler->matrix();
    assembler->assembleMass();
    gsSparseMatrix<> M = assembler->massMatrix();

    gsModalSolver<real_t> modal(K,M);
    modal.options().setInt("solver",2);
    modal.options().setInt("selectionRule",0);
    modal.options().setInt("sortRule",4);
    modal.options().setSwitch("verbose",true);
    modal.options().setInt("ncvFac",2);
    modal.options().setReal("shift",shift);

    gsStatus status;
    if (!sparse)
        status = modal.compute();
    else
        status = modal.computeSparse(10);//,2,Spectra::SortRule::LargestMagn,Spectra::SortRule::SmallestMagn);

    GISMO_ENSURE(status == gsStatus::Success,"Buckling solver failed");

    gsMatrix<> values = modal.values();
    gsMatrix<> vectors = modal.vectors();

    gsInfo<<"First eigenfrequency: "<<"\t[Analytical]: "<< omegas[0]<<"\t[Numerical]: "<<math::sqrt(values.at(0))<<"\n";


    gsInfo<< "First 10 eigenvalues:\n";
    for (index_t k = 0; k<10; k++)
        gsInfo<<"\t"<<std::setprecision(20)<<values.at(k)<<"\n";
    gsInfo<<"\n";

    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";
        int systemRet = system("mkdir -p ModalResults");
        GISMO_ASSERT(systemRet!=-1,"Something went wrong with calling the system argument");

        gsMultiPatch<> deformation = solution;
        gsMatrix<> modeShape;
        gsParaviewCollection collection("ModalResults/modes");

        int N = 1;
        if (!first)
          N = vectors.cols();
        for (index_t m=0; m<N; m++)
        {

          // Compute solution based on eigenmode with number 'mode'
          modeShape = modal.vector(m);
          assembler->constructSolution(modeShape, solution);

          // compute the deformation spline
          deformation = solution;
          deformation.patch(0).coefs() -= mp.patch(0).coefs();// assuming 1 patch here

          // Normalize mode shape amplitude in z coordinate
          real_t maxAmpl = std::max(math::abs(deformation.patch(0).coefs().col(2).maxCoeff()),math::abs(deformation.patch(0).coefs().col(2).minCoeff()));
          if (maxAmpl!=0.0)
          {
            deformation.patch(0).coefs() = deformation.patch(0).coefs()/maxAmpl;
          }

          gsField<> solField(mp,deformation);
          std::string fileName = "ModalResults/modes" + util::to_string(m);
          gsWriteParaview<>(solField, fileName, 5000);
          fileName = "modes" + util::to_string(m) + "0";
          collection.addPart(fileName + ".vts",m);

        }

        gsFunctionExpr<> analytical("0","0","sin(3.1415926535*x)*sin(3.1415926535*y)",3);
        gsPiecewiseFunction<> func(analytical);
        gsField<> an(mp,func);
        gsWriteParaview(an,"analytical");

        collection.save();
    }

    if (write)
    {
        int systemRet = system("mkdir -p ModalResults");
        GISMO_ASSERT(systemRet!=-1,"Something went wrong with calling the system argument");

        std::string wnM = "ModalResults/eigenvalues.txt";
        writeToCSVfile(wnM,values);
    }

    delete materialMatrix;
    delete assembler;

    return result;
}
#else//gsKLShell_ENABLED
int main(int argc, char *argv[])
{
    gsWarn<<"G+Smo is not compiled with the gsKLShell module.";
    return EXIT_FAILURE;
}
#endif

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
