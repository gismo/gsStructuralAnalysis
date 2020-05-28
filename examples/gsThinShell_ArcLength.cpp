/** @file gsThinShell_BucklingArcLength.cpp

    @brief Code for the arc-length method of a shell based on loads

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gismo.h>
#include <gismo_dev.h>

#include <gsIO/gsMatrixToFile.h>
#include <gsThinShell2/gsThinShellAssembler.h>

#include <gsThinShell/gsArcLengthIterator.h>

using namespace gismo;

template <class T>
gsMultiPatch<T> RectangularDomain(int n, int m, int p, int q, T L, T B, bool clamped = false, T offset = 0.1);
template <class T>
gsMultiPatch<T> RectangularDomain(int n, int p, T L, T B, bool clamped = false, T offset = 0.1);
template <class T>
gsMultiPatch<T> AnnularDomain(int n, int p, T R1, T R2);
template <class T>
gsMultiPatch<T> FrustrumDomain(int n, int p, T R1, T R2, T h);

void writeToCSVfile(std::string name, gsMatrix<> matrix)
{
    std::ofstream file(name.c_str());
    for(int  i = 0; i < matrix.rows(); i++){
        for(int j = 0; j < matrix.cols(); j++){
           std::string str = std::to_string(matrix(i,j));
           if(j+1 == matrix.cols()){
               file<<str;
           }else{
               file<<str<<',';
           }
        }
        file<<'\n';
    }
  }

int main (int argc, char** argv)
{
    // Input options
    int numElevate  = 1;
    int numHref     = 1;
    bool plot       = false;
    bool stress       = false;
    bool nonlinear  = false;
    bool first  = false;
    bool SingularPoint = false;
    bool quasiNewton = false;
    int quasiNewtonInt = -1;
    bool adaptive = false;
    int step = 10;
    int method = 1; // (1: Crisfield's method; 2: Riks' method)

    real_t thickness = 1e-1;
    real_t width = 0.1; // Width of the strip is equal to 0.1.
    real_t length = 1; // Length of the strip is equal to 1.
    real_t Area = thickness*width;

    real_t E_modulus     = 1;
    real_t PoissonRatio = 0;
    real_t Density = 1e0;
    gsMultiPatch<> mp;
    real_t eta = 0;
    real_t tau = 1e4;

    index_t Compressibility = 0;
    index_t material = 0;
    real_t Ratio = 7.0;

    real_t relax = 1.0;

    int testCase = 0;

    int result = 0;

    bool write = false;

    // Arc length method options
    real_t dL = 0; // General arc length
    real_t dLb = 0.5; // Ard length to find bifurcation
    real_t tol = 1e-6;
    real_t tolU = 1e-6;
    real_t tolF = 1e-3;

    std::string wn("data.csv");

    gsCmdLine cmd("Thin shell plate example.");
    cmd.addInt("r","hRefine",
       "Number of dyadic h-refinement (bisection) steps to perform before solving",
       numHref);
    cmd.addInt("t", "testcase",
        "Test case: 0: clamped-clamped, 1: pinned-pinned, 2: clamped-free",
       testCase);
    cmd.addInt("N", "maxsteps",
       "Maximum number of steps",
      step);
    cmd.addInt("e","degreeElevation",
      "Number of degree elevation steps to perform on the Geometry's basis before solving",
      numElevate);
    cmd.addInt("m","Method",
      "Arc length method; 1: Crisfield's method; 2: RIks' method.",
      method);
    cmd.addInt( "M", "Material", "Material law",  material );
    cmd.addInt( "c", "Compressibility", "1: compressible, 0: incompressible",  Compressibility );
    cmd.addReal("T","thickness",
       "thickenss of the plate",
       thickness);
    cmd.addReal("L","dLb",
       "arc length",
       dLb);
    cmd.addReal("l","dL",
       "arc length after bifurcation",
       dL);
    cmd.addReal("R","relaxation",
    "Relaxation factor for arc length method",
       relax);
    cmd.addSwitch("B","BifurcationPath",
       "Compute singular points and bifurcation paths",
       SingularPoint);
    cmd.addSwitch("Q","QuasiNewton",
       "Use the Quasi Newton method",
       quasiNewton);
    cmd.addReal("f","factor",
       "factor for bifurcation perturbation",
       tau);
    cmd.addInt("q","QuasiNewtonInt","Use the Quasi Newton method every INT iterations",quasiNewtonInt);
    cmd.addSwitch("A","adaptive", "Adaptive length ", adaptive);
    cmd.addSwitch("nl", "Nonlinear elasticity (otherwise linear)", nonlinear);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addSwitch("stress", "Plot stress in ParaView format", stress);
    cmd.addSwitch("first", "Plot only first", first);
    cmd.addSwitch("write", "Write convergence data to file", write);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    real_t aDim = 1;
    real_t bDim = 1;

    if (dL==0)
    {
      dL = dLb;
    }


    if (testCase==0 || testCase==1)
    {
      E_modulus = 2.886751346e12;
      thickness = 0.3464101615e-3;
      PoissonRatio = 0.0;

      real_t L = 1.0;
      real_t B = 0.1;
      Area = B*thickness;
      mp = RectangularDomain(numHref, 0, numElevate+2, 2, L, B);
    }
    else if (testCase==2 || testCase==3)
    {
      E_modulus = 1;
      thickness = 0.01;
      PoissonRatio = 0.0;
      real_t L = 1.0;
      real_t B = 0.01;
      Area = B*thickness;
      mp = RectangularDomain(numHref, 0, numElevate+2, 2, L, B);
    }
    else if (testCase==4 || testCase==5)
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
      Area = B*thickness;
      mp = RectangularDomain(numHref, numElevate+2, L, B);
    }
    else if (testCase==6)
    {
      // real_t L = 10.0e-3;
      // real_t B = 10.0e-3;
      // real_t mu = 10e3;
      // thickness = 0.25e-3;
      // if ((!Compressibility) && (material!=0))
      //   PoissonRatio = 0.5;
      // else
      //   PoissonRatio = 0.45;

      // E_modulus = 2*mu*(1+PoissonRatio);
      // gsDebugVar(E_modulus);

      real_t L = 9e-3;
      real_t B = 3e-3;
      real_t mu = 10e3;
      thickness = 0.9e-3;
      if ((!Compressibility) && (material!=0))
        PoissonRatio = 0.5;
      else
        PoissonRatio = 0.45;

      E_modulus = 2*mu*(1+PoissonRatio);
      gsDebugVar(E_modulus);

      mp = RectangularDomain(numHref, numElevate+2, L, B);
    }
    else if (testCase == 7)
    {
        thickness = 0.1;
        real_t mu = 4.225e5;
        if (!Compressibility)
          PoissonRatio = 0.5;
        else
          PoissonRatio = 0.45;
        E_modulus = 2*mu*(1+PoissonRatio);
        gsReadFile<>("quarter_sphere.xml", mp);

        for(index_t i = 0; i< numElevate; ++i)
          mp.patch(0).degreeElevate();    // Elevate the degree

        // h-refine
        for(index_t i = 0; i< numHref; ++i)
          mp.patch(0).uniformRefine();
    }
    else if (testCase == 8 || testCase == 9)
    {
        thickness = 0.1;
        real_t mu = 4.225;
        PoissonRatio = 0.5;
        E_modulus = 2*mu*(1+PoissonRatio);
        // gsReadFile<>("quarter_frustrum.xml", mp);

        // R1 is radius on bottom, R2 is radius on top
        mp = FrustrumDomain(numHref,numElevate+2,2.0,1.0,1.0);
    }
    else if (testCase == 10)
    {
        thickness = 0.5;
        real_t mu = 200;
        PoissonRatio = 0.5;
        E_modulus = 2*mu*(1+PoissonRatio);
        gsReadFile<>("quarterCircle.xml", mp);
        Ratio = 4;

        for(index_t i = 0; i< numElevate; ++i)
          mp.patch(0).degreeElevate();    // Elevate the degree

        // h-refine
        for(index_t i = 0; i< numHref; ++i)
          mp.patch(0).uniformRefine();

    }
    else if (testCase == 11)
    {
        thickness = 2e-3;
        PoissonRatio = 0.4;
        E_modulus = 168e9; // GPa
        gsReadFile<>("half_cylinder.xml", mp);
        Ratio = 4;

        for(index_t i = 0; i< numElevate; ++i)
          mp.patch(0).degreeElevate();    // Elevate the degree

        // h-refine
        for(index_t i = 0; i< numHref; ++i)
          mp.patch(0).uniformRefine();
    }
    else if (testCase==12 || testCase==13 )
    {
      real_t L = 2.5;
      real_t B = 1.0;
      thickness = 1e-3;
        if ((!Compressibility) && (material!=0))
          PoissonRatio = 0.5;
        else
          PoissonRatio = 0.499;
      E_modulus = 1e6;

      // We model symmetry over the width axis
      mp = RectangularDomain(1,numHref, 4, numElevate+2, L, B/2., true, 0.05);
    }
    else if (testCase==14 || testCase==15 )
    {
      real_t L = 2.5;
      real_t B = 1.0;
      thickness = 1e-3;
        if ((!Compressibility) && (material!=0))
          PoissonRatio = 0.5;
        else
          PoissonRatio = 0.499;
      E_modulus = 1e6;

      // We model symmetry over the width axis
      mp = RectangularDomain(1,numHref, 4, numElevate+2, L/2., B/2., true, 0.05);
    }
    else if (testCase == 16)
    {
        thickness = 1e-3;
        real_t mu = 11538.46e6;
        PoissonRatio = 0.4999999;
        E_modulus = 2*mu*(1+PoissonRatio);
        gsReadFile<>("pinched_cylinder2.xml", mp);

        for(index_t i = 0; i< numElevate; ++i)
          mp.patch(0).degreeElevate();    // Elevate the degree

        // h-refine
        for(index_t i = 0; i< numHref; ++i)
          mp.patch(0).uniformRefine();
    }
    else if (testCase==21  )
    {
      std::string fn;
      // thickness = 0.5*2.286;
      // E_modulus = 3102.75e2;
      // PoissonRatio = 0.3;

      thickness = 6.35;
      // thickness = 0.5*12.7;
      // thickness = 0.5*16.75;
      E_modulus = 3102.75;
      PoissonRatio = 0.3;

      fn = "scordelis_lo_roof_shallow.xml";

      gsReadFile<>(fn, mp);

      for(index_t i = 0; i< numElevate; ++i)
        mp.patch(0).degreeElevate();    // Elevate the degree

      // h-refine
      for(index_t i = 0; i< numHref; ++i)
        mp.patch(0).uniformRefine();
    }

    gsMultiBasis<> dbasis(mp);


    real_t EA = E_modulus*Area;
    real_t EI = 1.0/12.0*(width*math::pow(thickness,3))*E_modulus;

    real_t r = math::sqrt(EI/EA);
    gsInfo<<"EI = "<<EI<<"; EA = "<<EA<<"; r = "<<r<<"\n";

    gsInfo<<"Basis (patch 0): "<< mp.patch(0).basis() << "\n";

    // Boundary conditions
    std::vector< std::pair<patchSide,int> > clamped;
    gsBoundaryConditions<> BCs;
    gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();

    // Initiate Surface forces
    std::string tx("0");
    std::string ty("0");
    std::string tz("0");

    real_t displ = 1;
    gsFunctionExpr<> displ1(std::to_string( displ),3);
    gsFunctionExpr<> displ2(std::to_string(-displ),3);

    gsVector<> tmp(3);
    gsVector<> neu(3);
    tmp << 0, 0, 0;
    neu << 0, 0, 0;
    gsConstantFunction<> neuData(neu,3);

    // Buckling coefficient
    real_t fac = 1;
    // Unscaled load
    real_t Load = 0;

    std::string output = "solution";
    std::string dirname = "ArcLengthResults";
    real_t pressure = 0.0;

    if (testCase == 0)
    {
        // Pinned-Pinned
        tmp << 1e-4, 0, 0;
        neuData.setValue(tmp,3);
        // // Clamped-Clamped
        BCs.addCondition(boundary::west, condition_type::neumann, &neuData ); // unknown 0 - x
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

        fac = 1;

        // dL =  1e-2;
        // dLb = 1e-4;
        // tol = 1e-3;

        output = "Case" + std::to_string(testCase) + "solution";
        wn = output + "data.txt";
        SingularPoint = true;
    }
    if (testCase == 1)
    {
        Load = EA/width*1e-1;
        tmp << Load, 0, 0;
        // tmp << 0, 0, Load;
        neuData.setValue(tmp,3);
        // // Clamped-Clamped
        BCs.addCondition(boundary::west, condition_type::neumann, &neuData ); // unknown 0 - x
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y

        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

        BCs.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 2 );
        BCs.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 2 );


        fac = 0.5;

        // dL =  1e-3;
        // dLb = 1e-3;

        output = "Case" + std::to_string(testCase) + "solution";
        wn = output + "data.txt";
        SingularPoint = true;
    }
    if (testCase == 2)
    {
      BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
      BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
      BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

      BCs.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 2 );

      Load = 1e-4;
      gsVector<> point(2);
      gsVector<> load (3);
      point<< 1.0, 0.5 ;
      load << 0.0, 0.0, Load ;
      pLoads.addLoad(point, load, 0 );

      // dL =  1e-3;
      // dLb = 2e0;

      output = "Case" + std::to_string(testCase) + "solution";
      wn = output + "data.txt";
      SingularPoint = false;
    }
    if (testCase == 3)
    {
      // Clamped-Clamped
      BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
      BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
      BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

      BCs.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 2 );

      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
      BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y

      Load = 1e-1;

      // dL =  3e-0;
      // dLb = 0.8e-4;

      Load = 1e-1;
      gsVector<> point(2);
      gsVector<> load (3);
      point<< 1.0, 0.5 ;
      load << -Load, 0.0, 0.0 ;
      pLoads.addLoad(point, load, 0 );

      output = "Case" + std::to_string(testCase) + "solution_r" + std::to_string(numHref) + "_e" + std::to_string(numElevate);
      wn = output + "data.txt";
      SingularPoint = true;
    }
    if (testCase == 4) // Uniaxial tension; use with hyperelastic material model!
    {
      BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
      BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

      BCs.addCondition(boundary::east, condition_type::collapsed, 0, 0, false, 0 ); // unknown 1 - y
      BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z.


      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 1 - y
      BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 1 - y

      // tmp << 1,0,0;
      // neuData.setValue(tmp,3);
      // // // Clamped-Clamped
      // BCs.addCondition(boundary::east, condition_type::neumann, &neuData ); // unknown 0 - x

      Load = 1e0;
      gsVector<> point(2);
      gsVector<> load (3);
      point<< 1.0, 0.5 ;
      load << Load,0.0, 0.0;
      pLoads.addLoad(point, load, 0 );

      output = "Case" + std::to_string(testCase) + "solution";
      wn = output + "data.txt";
      SingularPoint = true;
    }
    if (testCase == 5) // Bi-axial tension; use with hyperelastic material model!
    {
      BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
      BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

      BCs.addCondition(boundary::east, condition_type::collapsed, 0, 0, false, 0 ); // unknown 1 - y
      BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z.

      BCs.addCondition(boundary::north, condition_type::collapsed, 0, 0, false, 1 ); // unknown 1 - y
      BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z.

      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 1 - y

      real_t load_factor = 1;
      Load = 1e0;
      gsVector<> point(2);
      gsVector<> load (3);
      point<< 1.0, 0.5 ;
      load << Load,0.0, 0.0;
      pLoads.addLoad(point, load, 0 );

      point<< 0.5, 1.0 ;
      load << 0.0, Load/load_factor, 0.0;
      pLoads.addLoad(point, load, 0 );

      output = "Case" + std::to_string(testCase) + "solution";
      wn = output + "data.txt";
      SingularPoint = true;
    }
    // else if (testCase == 6) // Plate with pressure
    // {
    //   for (index_t k=0; k!=3; k++)
    //   {
    //     BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, k ); // unknown 1 - x
    //     BCs.addCondition(boundary::east,  condition_type::dirichlet, 0, 0, false, k ); // unknown 1 - x
    //     BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, k ); // unknown 1 - x
    //     BCs.addCondition(boundary::west,  condition_type::dirichlet, 0, 0, false, k ); // unknown 1 - x
    //   }

    //   BCs.addCondition(boundary::north, condition_type::clamped, 0, 0, false, 2 ); // unknown 2 - z.
    //   BCs.addCondition(boundary::east,  condition_type::clamped, 0, 0, false, 2 ); // unknown 2 - z.
    //   BCs.addCondition(boundary::south, condition_type::clamped, 0, 0, false, 2 ); // unknown 1 - y
    //   BCs.addCondition(boundary::west,  condition_type::clamped, 0, 0, false, 2 ); // unknown 1 - y

    //   pressure = 1.0;

    //   output = "Case" + std::to_string(testCase) + "solution";
    //   wn = output + "data.txt";
    //   SingularPoint = false;
    // }
    else if (testCase == 6) // Plate with pressure
    {
      BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 2 - z.
      BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 2 - z.

      BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 2 - z.
      BCs.addCondition(boundary::east, condition_type::collapsed, 0, 0, false, 0 ); // unknown 2 - z.


      BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z.
      BCs.addCondition(boundary::east,  condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z.
      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 1 - y
      BCs.addCondition(boundary::west,  condition_type::dirichlet, 0, 0, false, 2 ); // unknown 1 - y

      Load = 1e0;
      gsVector<> point(2);
      gsVector<> load (3);
      point<< 1.0, 0.5 ;
      load << Load,0.0, 0.0;
      pLoads.addLoad(point, load, 0 );

      output = "Case" + std::to_string(testCase) + "solution";
      wn = output + "data.txt";
      SingularPoint = false;


    }

    else if (testCase == 7)
    {
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

        // Symmetry in x-direction:
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 );
        BCs.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 1 );
        BCs.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 2 );

        // Symmetry in y-direction:
        BCs.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 0 );
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 );
        BCs.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 2 );

        // Pressure
        pressure = 1e3;

        output = "Balloon_solution";
        wn = dirname + "/" + output + "data.txt";

    }
    else if (testCase == 8)
    {
        Load = -1;
        neu << 0, 0, Load;
        neuData.setValue(neu,3);

        BCs.addCondition(boundary::north, condition_type::neumann, &neuData );
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 2 - z
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 2 - z
        BCs.addCondition(boundary::north, condition_type::collapsed, 0, 0, false, 2 ); // unknown 1 - y

        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

        // Symmetry in x-direction:
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 );
        BCs.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 1 );
        BCs.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 2 );

        // Symmetry in y-direction:
        BCs.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 0 );
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 );
        BCs.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 2 );

        dirname = dirname + "/" + "Frustrum_-r=" + std::to_string(numHref) + "-e" + std::to_string(numElevate) + "_solution";
        output = "solution";
        wn = dirname + "/" + output + "data.txt";
    }
    else if (testCase == 9)
    {
        Load = -1;
        neu << 0, 0, Load;
        neuData.setValue(neu,3);

        BCs.addCondition(boundary::north, condition_type::neumann, &neuData );
        // BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 2 - z
        // BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 2 - z
        BCs.addCondition(boundary::north, condition_type::collapsed, 0, 0, false, 2 ); // unknown 1 - y

        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

        // Symmetry in x-direction:
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 );
        BCs.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 1 );
        BCs.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 2 );

        // Symmetry in y-direction:
        BCs.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 0 );
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 );
        BCs.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 2 );

        dirname = dirname + "/" + "Frustrum2_-r=" + std::to_string(numHref) + "-e" + std::to_string(numElevate) + "_solution";
        output = "solution";
        wn = dirname + "/" + output + "data.txt";
    }
    else if (testCase == 10)
    {
        pressure = 1;

        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 2 - z
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 0 - x
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 1 - y

        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 );
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 );
        BCs.addCondition(boundary::west, condition_type::collapsed, 0, 0, false, 2 );
        BCs.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 2 );

        // Symmetry in x-direction:
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 0 );
        BCs.addCondition(boundary::north, condition_type::clamped, 0, 0, false, 1 );
        BCs.addCondition(boundary::north, condition_type::clamped, 0, 0, false, 2 );

        // Symmetry in y-direction:
        BCs.addCondition(boundary::south, condition_type::clamped, 0, 0, false, 0 );
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 );
        BCs.addCondition(boundary::south, condition_type::clamped, 0, 0, false, 2 );

        output = "Case" + std::to_string(testCase) + "solution";
        wn = output + "data.txt";
        SingularPoint = false;
    }
    else if (testCase == 11)
    {
      Load = -1;
      neu << 0, 0, Load;
      neuData.setValue(neu,3);

      BCs.addCondition(boundary::north, condition_type::neumann, &neuData );
      BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 );
      BCs.addCondition(boundary::north, condition_type::clamped, 0, 0, false, 2 );

      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 0 );
      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 );
      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 );
      BCs.addCondition(boundary::south, condition_type::clamped, 0, 0, false, 2 );

      output = "Case" + std::to_string(testCase) + "solution";
      wn = output + "data.txt";
      SingularPoint = false;
    }
    // Anti-symmetric
    else if (testCase == 12)
    {
      BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
      BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 2 - z.
      BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z.
      BCs.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 2 ); // unknown 2 - z.

      BCs.addCondition(boundary::east, condition_type::collapsed, 0, 0, false, 0 ); // unknown 1 - y
      BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 2 - z.
      BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z.
      BCs.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 2 ); // unknown 2 - z.

      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 2 - z.
      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z.

      Load = 1e0;
      gsVector<> point(2);
      gsVector<> load (3);
      point<< 1.0, 0.5 ;
      load << Load,0.0, 0.0;
      pLoads.addLoad(point, load, 0 );

      // dL =  750;
      // dLb = 750;

      dirname = dirname + "/" + "Sheet_Symm_Half_" + "-r" + std::to_string(numHref) + "-e" + std::to_string(numElevate) + "-M" + std::to_string(material) + "-c" + std::to_string(Compressibility);
      output = "solution";
      wn = dirname + "/" + output + "data.txt";
      SingularPoint = true;
    }
    // Symmetric
    else if (testCase == 13)
    {
      BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
      BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 2 - z.
      BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z.
      BCs.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 2 ); // unknown 2 - z.

      BCs.addCondition(boundary::east, condition_type::collapsed, 0, 0, false, 0 ); // unknown 1 - y
      BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 2 - z.
      BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z.
      BCs.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 2 ); // unknown 2 - z.

      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 2 - z.
      BCs.addCondition(boundary::south, condition_type::clamped, 0, 0, false, 2 ); // unknown 2 - z.

      Load = 1e0;
      gsVector<> point(2);
      gsVector<> load (3);
      point<< 1.0, 0.5 ;
      load << Load,0.0, 0.0;
      pLoads.addLoad(point, load, 0 );

      // dL =  750;
      // dLb = 750;

      dirname = dirname + "/" + "Sheet_Asymm_Half_" + "-r" + std::to_string(numHref) + "-e" + std::to_string(numElevate) + "-M" + std::to_string(material) + "-c" + std::to_string(Compressibility);
      output = "solution";
      wn = dirname + "/" + output + "data.txt";
      SingularPoint = true;
    }
    // Anti-symmetric
    else if (testCase == 14)
    {
      BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
      BCs.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 2 ); // unknown 2 - z.

      BCs.addCondition(boundary::east, condition_type::collapsed, 0, 0, false, 0 ); // unknown 1 - y
      BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 2 - z.
      BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z.
      BCs.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 2 ); // unknown 2 - z.

      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 2 - z.
      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z.

      Load = 1e0;
      gsVector<> point(2);
      gsVector<> load (3);
      point<< 1.0, 0.5 ;
      load << Load,0.0, 0.0;
      pLoads.addLoad(point, load, 0 );

      // dL =  750;
      // dLb = 750;

      dirname = dirname + "/" + "Sheet_Symm_Quarter_" + "-r" + std::to_string(numHref) + "-e" + std::to_string(numElevate) + "-M" + std::to_string(material) + "-c" + std::to_string(Compressibility);
      output = "solution";
      wn = dirname + "/" + output + "data.txt";
      SingularPoint = true;
    }
    // Symmetric
    else if (testCase == 15)
    {
      BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
      BCs.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 2 ); // unknown 2 - z.

      BCs.addCondition(boundary::east, condition_type::collapsed, 0, 0, false, 0 ); // unknown 1 - y
      BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 2 - z.
      BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z.
      BCs.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 2 ); // unknown 2 - z.

      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 2 - z.
      BCs.addCondition(boundary::south, condition_type::clamped, 0, 0, false, 2 ); // unknown 2 - z.

      Load = 1e0;
      gsVector<> point(2);
      gsVector<> load (3);
      point<< 1.0, 0.5 ;
      load << Load,0.0, 0.0;
      pLoads.addLoad(point, load, 0 );

      // dL =  750;
      // dLb = 750;

      dirname = dirname + "/" + "Sheet_Asymm_Quarter_" + "-r" + std::to_string(numHref) + "-e" + std::to_string(numElevate) + "-M" + std::to_string(material) + "-c" + std::to_string(Compressibility);
      output = "solution";
      wn = dirname + "/" + output + "data.txt";
      SingularPoint = true;
    }
    else if (testCase == 16)
    {
      // Symmetry in y-direction for back side
      BCs.addCondition(boundary::north, condition_type::clamped, 0, 0, false, 0 );
      BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 );
      BCs.addCondition(boundary::north, condition_type::clamped, 0, 0, false, 2 );

      // Diaphragm conditions for left side
      BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
      BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

      // Symmetry in x-direction: for right side
      BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 );
      BCs.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 1 );
      BCs.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 2 );

      // Symmetry in z-direction:for the front side
      BCs.addCondition(boundary::south, condition_type::clamped, 0, 0, false, 0 );
      BCs.addCondition(boundary::south, condition_type::clamped, 0, 0, false, 1 );
      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 );

      // Surface forces
      tmp.setZero();

      // Point loads
      gsVector<> point(2); point<< 1.0, 1.0 ;
      gsVector<> load (3); load << 0.0, 0.0, -0.25 ;
      pLoads.addLoad(point, load, 0 );

      output = "Case" + std::to_string(testCase) + "solution";
      wn = output + "data.txt";
      SingularPoint = true;
    }


    else if (testCase == 21)
    {
      // Diaphragm conditions
      BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
      BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
      BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
      // BCs.addCornerValue(boundary::southwest, 0.0, 0, 0); // (corner,value, patch, unknown)

      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

      Load = -1e1;
      // Point loads
      gsVector<> point(2);
      gsVector<> load (3);
      point<< 0.5, 0.5 ;
      load << 0.0, 0.0, Load ;
      pLoads.addLoad(point, load, 0 );

      dirname = dirname + "/" +  "Roof1_t="+ std::to_string(thickness) + "-r=" + std::to_string(numHref) + "-e" + std::to_string(numElevate) +"_solution";
      output = "solution";
      wn = dirname + "/" + output + "data.txt";
      SingularPoint = false;
    }

    std::string commands = "mkdir -p " + dirname;
    const char *command = commands.c_str();
    system(command);

    // plot geometry
    if (plot)
      gsWriteParaview(mp,"mp",1000,true);

    if (write)
    {
      std::ofstream file;
      file.open(wn,std::ofstream::out);
      file  << std::setprecision(20)
            << "Deformation norm" << ","
            << "Left end - x" << ","
            << "Left end - y" << ","
            << "Left end - z" << ","
            << "Mid point - x" << ","
            << "Mid point - y" << ","
            << "Mid point - z" << ","
            << "Right end - x" << ","
            << "Right end - y" << ","
            << "Right end - z" << ","
            << "Lambda" << ","
            << "Indicator"
            << "\n";
      file.close();
    }

    gsInfo<<"Results will be written in folder: "<<dirname<<"\n";

    gsFunctionExpr<> surfForce(tx,ty,tz,3);
    gsConstantFunction<> pressFun(pressure,3);
    // Initialise solution object
    gsMultiPatch<> mp_def = mp;
    gsSparseSolver<>::LU solver;

    // Linear isotropic material model
    gsConstantFunction<> force(tmp,3);
    gsFunctionExpr<> t(std::to_string(thickness), 3);
    gsFunctionExpr<> E(std::to_string(E_modulus),3);
    gsFunctionExpr<> nu(std::to_string(PoissonRatio),3);
    gsFunctionExpr<> rho(std::to_string(Density),3);
    gsConstantFunction<> ratio(7.0,3);

    real_t mu = E_modulus / (2 * (1 + PoissonRatio));
    gsConstantFunction<> alpha1(1.3,3);
    gsConstantFunction<> mu1(6.3e5/4.225e5*mu,3);
    gsConstantFunction<> alpha2(5.0,3);
    gsConstantFunction<> mu2(0.012e5/4.225e5*mu,3);
    gsConstantFunction<> alpha3(-2.0,3);
    gsConstantFunction<> mu3(-0.1e5/4.225e5*mu,3);

    // gsMaterialMatrix materialMatrixNonlinear(mp,mp_def,t,E,nu,rho);
    std::vector<gsFunction<>*> parameters(3);
    parameters[0] = &E;
    parameters[1] = &nu;
    parameters[2] = &ratio;
    gsMaterialMatrix materialMatrixNonlinear(mp,mp_def,t,parameters,rho);

    std::vector<gsFunction<>*> parameters2(8);
    if (material==14)
    {
        parameters2[0] = &E;
        parameters2[1] = &nu;
        parameters2[2] = &mu1;
        parameters2[3] = &alpha1;

        parameters2[4] = &mu2;
        parameters2[5] = &alpha2;

        parameters2[6] = &mu3;
        parameters2[7] = &alpha3;
        materialMatrixNonlinear.setParameters(parameters2);
    }


    materialMatrixNonlinear.options().setInt("MaterialLaw",material);
    materialMatrixNonlinear.options().setInt("Compressibility",Compressibility);

    // Construct assembler object
    gsThinShellAssembler assembler(mp,dbasis,BCs,surfForce,materialMatrixNonlinear);
    assembler.setPointLoads(pLoads);
    if (pressure!= 0.0)
        assembler.setPressure(pressFun);


    // Function for the Jacobian
    std::function<gsSparseMatrix<real_t> (gsVector<real_t> const &)> Jacobian;
    Jacobian = [&assembler,&mp_def](gsVector<real_t> const &x)
    {
      // assembler.constructSolution(x); // DOES NOT WORK!!
      assembler.constructSolution(x,mp_def);
      assembler.assembleMatrix(mp_def);
      gsSparseMatrix<real_t> m = assembler.matrix();
      // gsInfo<<"matrix = \n"<<m.toDense()<<"\n";
      return m;
    };
    // Function for the Residual
    std::function<gsVector<real_t> (gsVector<real_t> const &, real_t, gsVector<real_t> const &) > Residual;
    Residual = [&assembler,&mp_def](gsVector<real_t> const &x, real_t lam, gsVector<real_t> const &force)
    {
      // assembler.assembleVector(x); // DOES NOT WORK!!
      assembler.constructSolution(x,mp_def);
      assembler.assembleVector(mp_def);
      gsVector<real_t> Fint = -(assembler.rhs() - force);
      gsVector<real_t> result = Fint - lam * force;

      // gsDebugVar(Fint);
      // gsDebugVar(result);
      // The residual is now defined as the internal forces minus lam*force
      // gsVector<real_t> r =
      // gsVector<> result = Fint - lam * force;
      // gsDebugVar(lam * force);
      // gsDebugVar(assembler.rhs());
      // gsDebugVar(Fint);
      // gsDebugVar((lam*force).sum());
      return result; // - lam * force;
    };
    // Assemble linear system to obtain the force vector
    assembler.assemble();
    gsVector<> Force = assembler.rhs();

    // gsDebugVar(Force.sum());

    gsArcLengthIterator<real_t> arcLength(Jacobian, Residual, Force);
    arcLength.setLength(dLb); // dLb
    arcLength.setBifurcationMethod("Determinant");
    arcLength.setLength(dLb,adaptive,5); // dLb
    arcLength.setTau(tau);
    arcLength.setTolerance(tol); //tol
    arcLength.setToleranceU(tolU);
    arcLength.setToleranceF(tolF);
    arcLength.setMaxIterations(50);
    arcLength.verbose();
    arcLength.setAngleDeterminationMethod(0);
    arcLength.setPhi(0.0);
    if (testCase==4 || testCase==5 || testCase==8 || testCase==9 || testCase==21)
      arcLength.setPhi(0.0);
    // arcLength.setRelaxation(relax);

    if (method==1)
      arcLength.setMethod("Crisfield");
    else if (method==2)
      arcLength.setMethod("Riks");
    else if (method==3)
      arcLength.setMethod("ConsistentCrisfield");

    if (quasiNewton)
      arcLength.quasiNewton();

    if (quasiNewtonInt>0)
      arcLength.quasiNewton(quasiNewtonInt);

    gsParaviewCollection collection(dirname + "/" + output);
    gsMultiPatch<> deformation = mp;

    // Make objects for previous solutions
    real_t Lold = 0;
    gsMatrix<> Uold = Force;
    Uold.setZero();

    gsMatrix<> solVector;
    real_t indicator = 0.0;
    arcLength.setIndicator(indicator); // RESET INDICATOR
    for (index_t k=0; k<step; k++)
    {
      gsInfo<<"Load step "<< k<<"\n";
      // assembler.constructSolution(solVector,solution);
      arcLength.step();

      // gsInfo<<"m_U = "<<arcLength.solutionU()<<"\n";
      if (!(arcLength.converged()))
      {
        gsInfo<<"Error: Loop terminated, arc length method did not converge.\n";
        arcLength.setLength(dLb/2.);
        arcLength.setSolution(Uold,Lold);
        k -= 1;
        continue;
        // if (plot)
        // {
        //   solVector = arcLength.solutionU();
        //   Uold = solVector;
        //   Lold = arcLength.solutionL();
        //   assembler.constructSolution(solVector,mp_def);

        //   deformation = mp_def;
        //   deformation.patch(0).coefs() -= mp.patch(0).coefs();// assuming 1 patch here

        //   gsField<> solField(mp,deformation);
        //   std::string fileName = dirname + "/" + output + util::to_string(k);
        //   gsWriteParaview<>(solField, fileName, 5000);
        //   fileName = output + util::to_string(k) + "0";
        //   collection.addTimestep(fileName,k,".vts");
        // }
        // break;
      }

      if (SingularPoint)
      {
        arcLength.computeStability(arcLength.solutionU(),quasiNewton);
        if (arcLength.stabilityChange())
        {
          gsInfo<<"Bifurcation spotted!"<<"\n";
          arcLength.computeSingularPoint(1e-6, 5, Uold, Lold, 1e-10, 0, false);
          arcLength.switchBranch();
          arcLength.setLength(dL);
        }
      }
      indicator = arcLength.indicator();

      solVector = arcLength.solutionU();
      Uold = solVector;
      Lold = arcLength.solutionL();
      assembler.constructSolution(solVector,mp_def);

      gsMatrix<> pts(2,1);
      pts<<0.5,0.5;
      if (testCase==8 || testCase==9)
      {
        pts.resize(2,3);
        pts.col(0)<<0.0,1.0;
        pts.col(1)<<0.5,1.0;
        pts.col(2)<<1.0,1.0;
      }
      gsDebugVar(assembler.computePrincipalStretches(pts,mp_def,0));

      deformation = mp_def;
      deformation.patch(0).coefs() -= mp.patch(0).coefs();// assuming 1 patch here

      // gsDebugVar(mp_def.patch(0).coefs());

      if (stress)
      {
        gsPiecewiseFunction<> stresses;
        assembler.constructStress(mp_def,stresses,stress_type::principal_stretch);
        gsField<> stressField(mp,stresses, true);
        gsWriteParaview( stressField, "stress", 5000);
      }


      gsInfo<<"pressures:\n"<<pressure*arcLength.solutionL()<<"\n"
                            <<pressure*arcLength.solutionL() * assembler.getArea(mp) / assembler.getArea(mp_def)<<"\n";


      if (plot)
      {
        gsField<> solField(mp,deformation);
        std::string fileName = dirname + "/" + output + util::to_string(k);
        gsWriteParaview<>(solField, fileName, 5000);
        fileName = output + util::to_string(k) + "0";
        collection.addTimestep(fileName,k,".vts");
      }
      // gsDebugVar(mp_def.patch(0).coefs());


      if (write)
      {
        gsMatrix<> P(2,1);
        // Compute end point displacement
        if (testCase==8 || testCase==9)
          P<<0.0,1.0;
        else
          P<<0.0,0.5;

        gsMatrix<> left;
        deformation.patch(0).eval_into(P,left);

        if (testCase==8 || testCase==9)
          P<<0.5,1.0;
        else
          P<<0.5,0.5;

        gsMatrix<> mid;
        deformation.patch(0).eval_into(P,mid);

        if (testCase==8 || testCase==9)
          P<<1.0,1.0;
        else
          P<<1.0,0.5;

        gsMatrix<> right;
        deformation.patch(0).eval_into(P,right);

        std::ofstream file;
        file.open(wn,std::ofstream::out | std::ofstream::app);
        if (testCase!=12)
        {
          file  << std::setprecision(6)
                << arcLength.solutionU().norm() << ","
                << left.at(0) << ","
                << left.at(1) << ","
                << left.at(2) << ","
                << mid.at(0) << ","
                << mid.at(1) << ","
                << mid.at(2) << ","
                << right.at(0) << ","
                << right.at(1) << ","
                << right.at(2) << ","
                << -arcLength.solutionL() << ","
                << indicator << ","
                << "\n";
        }
        else
        {
          gsVector<> w(101);
          gsMatrix<> Q(2,1);
          for (int k = 0; k <= 100; k ++)
          {
            Q<<0.0,1.0*k/100;

            gsMatrix<> res;
            deformation.patch(0).eval_into(Q,res);
            w.at(k) = res.at(2);
            // gsInfo<<res.at(0)<<","<<res.at(1)<<","<<res.at(2)<<"\n";
          }
          file  << std::setprecision(6)
                << arcLength.solutionU().norm() << ","
                << left.at(0) << ","
                << left.at(1) << ","
                // << left.at(2) << ","
                << std::max(abs(w.maxCoeff()),abs(w.minCoeff())) << ","
                << mid.at(0) << ","
                << mid.at(1) << ","
                << mid.at(2) << ","
                << right.at(0) << ","
                << right.at(1) << ","
                << right.at(2) << ","
                << -arcLength.solutionL() << ","
                << indicator << ","
                << "\n";
        }
        file.close();

      }
    }
    if (plot)
      collection.save();

  return result;
}

template <class T>
gsMultiPatch<T> RectangularDomain(int n, int p, T L, T B, bool clamped, T clampoffset)
{
  gsMultiPatch<T> mp = RectangularDomain(n, n, p, p, L, B, clamped, clampoffset);
  return mp;
}

template <class T>
gsMultiPatch<T> RectangularDomain(int n, int m, int p, int q, T L, T B, bool clamped, T clampoffset)
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

  if (clamped)
  {
    T knotval;
    knotval = kv0.uValue(1);
    kv0.insert(std::min(clampoffset,knotval/2.));

    knotval = kv0.uValue(kv0.uSize()-2);
    kv0.insert(std::max(1-clampoffset,knotval/2.));
  }

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
gsMultiPatch<T> AnnularDomain(int n, int p, T R1, T R2)
{
  // -------------------------------------------------------------------------
  // --------------------------Make beam geometry-----------------------------
  // -------------------------------------------------------------------------
  int dim = 3; //physical dimension
  gsKnotVector<> kv0;
  kv0.initUniform(0,1,0,3,1);
  gsKnotVector<> kv1;
  kv1.initUniform(0,1,0,3,1);

  // Make basis
  // gsTensorNurbsBasis<2,T> basis(kv0,kv1);

  // Initiate coefficient matrix
  gsMatrix<> coefs(9,dim);

  coefs<<R1,0,0,
  (R1+R2)/2,0,0,
  R2,0,0,
  R1,R1,0,
  (R1+R2)/2,(R1+R2)/2,0,
  R2,R2,0,
  0,R1,0,
  0,(R1+R2)/2,0,
  0,R2,0;

  gsMatrix<> weights(9,1);
  weights<<1,1,1,
  0.707106781186548,0.707106781186548,0.707106781186548,
  1,1,1;

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

  // Refine n times
  for(index_t i = 0; i< n; ++i)
      mp.patch(0).uniformRefine();

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
