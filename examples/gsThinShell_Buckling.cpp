/** @file gsThinShell_Buckling.cpp

    @brief Example to compute eigenvalues and eigenmodes of a buckled shell

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#include <gismo.h>

#include <gsKLShell/gsThinShellAssembler.h>
#include <gsKLShell/getMaterialMatrix.h>
#include <gsStructuralAnalysis/gsBucklingSolver.h>

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
gsMultiPatch<T> RectangularDomain(int n, int m, int p, int q, T L, T B);
template <class T>
gsMultiPatch<T> RectangularDomain(int n, int p, T L, T B);
template <class T>
gsMultiPatch<T> Rectangle(T L, T B);
template <class T>
gsMultiPatch<T> AnnularDomain(int n, int p, T R1, T R2);

int main (int argc, char** argv)
{
    // Input options
    int numElevate  = 1;
    int numHref     = 1;
    int numElevateL = -1;
    int numHrefL    = -1;
    bool plot       = false;
    bool sparse     = false;
    bool nonlinear  = false;
    bool first  = false;
    int mode = 0;


    real_t E_modulus     = 1e8;
    real_t PoissonRatio = 0;
    real_t Density = 1e0;

    real_t thickness = 1e-3;
    real_t aDim = 1.0;
    real_t bDim = 1.0;

    index_t Compressibility = 0;
    index_t material = 0;
    real_t Ratio = 7.0;
    bool composite = false;
    index_t impl = 1; // 1= analytical, 2= generalized, 3= spectral

    real_t fac = 1;

    real_t shift = 0.0;

    int testCase = 0;

    int result = 0;

    bool write = false;

    std::string assemberOptionsFile("options/solver_options.xml");

    gsCmdLine cmd("Buckling analysis for thin shells.");
    cmd.addString( "f", "file", "Input XML file for assembler options", assemberOptionsFile );

    cmd.addInt("t", "testcase", "Test case: 0: clamped-clamped, 1: pinned-pinned, 2: clamped-free", testCase);

    cmd.addInt("r","hRefine", "Number of dyadic h-refinement (bisection) steps to perform before solving", numHref);
    cmd.addInt("e","degreeElevation", "Number of degree elevation steps to perform on the Geometry's basis before solving", numElevate);
    cmd.addInt("R","hRefine2", "Number of dyadic h-refinement (bisection) steps to perform before solving (secondary direction)", numHrefL);
    cmd.addInt("E","degreeElevation2", "Number of degree elevation steps to perform on the Geometry's basis before solving (secondary direction)", numElevateL);

    cmd.addInt( "M", "Material", "Material law",  material );
    cmd.addInt( "c", "Compressibility", "1: compressible, 0: incompressible",  Compressibility );
    cmd.addInt( "I", "Implementation", "Implementation: 1= analytical, 2= generalized, 3= spectral",  impl );
    cmd.addSwitch("composite", "Composite material", composite);

    cmd.addInt("m", "nmode",
               "Mode shape number, starting from 0",
              mode);

    cmd.addReal("T","hdim", "thickness of the plate", thickness);
    cmd.addReal("a","adim", "dimension a", aDim);
    cmd.addReal("b","bdim", "dimension b", bDim);

    cmd.addReal("F","fac", "factor linear problem", fac);

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

    gsMultiPatch<> mp;

    if (numHrefL==-1)
      numHrefL = numHref;
    if (numElevateL==-1)
      numElevateL = numElevate;

    real_t length,width;
    if (testCase==0 || testCase==1)
    {
      thickness = 2*0.005313292845913;
      length = 1.0;
      width = 0.1;
      E_modulus = 1e10;
      PoissonRatio = 0.0;

      mp = RectangularDomain(numHref,1,numElevate,2,length,width);
    }
    else if (testCase==2)
    {
      E_modulus = 1;
      thickness = 0.01;
      PoissonRatio = 0.0;
      length = 1.0;
      width = 0.01;

      mp = RectangularDomain(numHref,1,numElevate,1,length,width);
    }
    else if (testCase==3 || testCase==4 || testCase==5 || testCase==6)
    {
      length = aDim;
      width = bDim;
      thickness = 1e-2;
      E_modulus = 200e9;
      PoissonRatio = 0.3;

      // length = 1.0;
      // width = 1.0;
      // thickness = 0.5*0.01;
      // E_modulus = 200e9;
      // PoissonRatio = 0.3;

      mp = RectangularDomain(numHrefL,numHref,numElevateL,numElevate,length,width);
    }
    else if (testCase==8)
    {
      E_modulus = 1e6;
      PoissonRatio = 0.3;
      gsDebug<<"E = "<<E_modulus<<"; nu = "<<PoissonRatio<<"\n";

      aDim = 2;
      bDim = 1;
      thickness = 1e-3;

      mp = Rectangle(aDim,bDim);

      for(index_t i = 0; i< numElevate; ++i)
        mp.patch(0).degreeElevate();    // Elevate the degree

      // h-refine
      for(index_t i = 0; i< numHref; ++i)
        mp.patch(0).uniformRefine();

    }
    else if (testCase==10)
    {
      if ((!Compressibility) && (material!=0))
        PoissonRatio = 0.5;
      else
        PoissonRatio = 0.499;

      real_t mu;
      if (material==3||material==13||material==23)
        mu = 440200;
      else
        mu = 2*1.91e5;
      PoissonRatio = 0.5;
      E_modulus = 2*mu*(1+PoissonRatio);

      gsDebug<<"E = "<<E_modulus<<"; nu = "<<PoissonRatio<<"; mu = "<<mu<<"\n";

      aDim = 0.28;
      bDim = 0.14;
      thickness = 140e-6;

      Ratio = 2.5442834138486314;

      mp = RectangularDomain(numHrefL, numHref, numElevateL+2, numElevate + 2, aDim, bDim);
    }
    else if (testCase==16 )
    {
      if ((!Compressibility) && (material!=0))
        PoissonRatio = 0.5;
      else
        PoissonRatio = 0.499;

      real_t mu, C01,C10;
      if (material==3||material==13||material==23)
      {
        C10 = (0.5-1/22.)*1e6;      // c1/2
        C01 = (1/22.)*1e6;          // c2/2

        Ratio = C10/C01;
        mu = 2*(C01+C10);
      }
      else
      {
        C10 = (0.5)*1e6;
        mu = 2*C10;
      }
      E_modulus = 2*mu*(1+PoissonRatio);
      gsDebug<<"E = "<<E_modulus<<"; nu = "<<PoissonRatio<<"; mu = "<<mu<<"; ratio = "<<Ratio<<"\n";

      aDim = 0.28;
      bDim = 0.14;
      thickness = 0.14e-3;

      // We model symmetry over the width axis
      mp = Rectangle(aDim/2., bDim/2.);//, true, 0.001);

      for(index_t i = 0; i< numElevate; ++i)
        mp.patch(0).degreeElevate();    // Elevate the degree

      // h-refine
      for(index_t i = 0; i< numHref; ++i)
        mp.patch(0).uniformRefine();
    }
    gsMultiBasis<> dbasis(mp);


    gsInfo<<"Basis (patch 0): "<< mp.patch(0).basis() << "\n";

    // Boundary conditions
    gsBoundaryConditions<> BCs;
    gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();

    // Initiate Surface forces
    std::string tx("0");
    std::string ty("0");
    std::string tz("0");

    gsVector<> tmp(3);
    gsVector<> neu(3);
    tmp << 0, 0, 0;
    neu << 0, 0, 0;
    gsConstantFunction<> neuData(neu,3);
    gsConstantFunction<> neuData2(neu,3);

    // Buckling coefficient
    real_t Load = -1e4;
    real_t pressure = 0.0;

    if (testCase == 0)
    {
        neu << Load, 0, 0;
        neuData.setValue(neu,3);
        // // Clamped-Clamped
        BCs.addCondition(boundary::west, condition_type::neumann, &neuData ); // unknown 0 - x
        // BCs.addCondition(boundary::west, condition_type::dirichlet, &displ1, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0,false,1 ); // unknown 1 - y
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0,false,2 ); // unknown 2 - z

        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0,false,1 ); // unknown 1 - y
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0,false,1 ); // unknown 1 - y


        // BCs.addCondition(boundary::east, condition_type::dirichlet, &displ2, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0,false,0 ); // unknown 0 - x
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0,false,1 ); // unknown 1 - y
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0,false,2 ); // unknown 2 - z
    }
    else if (testCase == 1)
    {
        neu << Load, 0, 0;
        neuData.setValue(neu,3);
        // // Clamped-Clamped
        BCs.addCondition(boundary::west, condition_type::neumann, &neuData ); // unknown 0 - x
        // BCs.addCondition(boundary::west, condition_type::dirichlet, &displ1, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0,false,1 ); // unknown 1 - y
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0,false,2 ); // unknown 2 - z

        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0,false,1 ); // unknown 1 - y
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0,false,1 ); // unknown 1 - y


        // BCs.addCondition(boundary::east, condition_type::dirichlet, &displ2, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0,false,0 ); // unknown 0 - x
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0,false,1 ); // unknown 1 - y
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0,false,2 ); // unknown 2 - z

        BCs.addCondition(boundary::east, condition_type::clamped,0,0,false,2);
        BCs.addCondition(boundary::west, condition_type::clamped,0,0,false,2);
    }
    else if (testCase == 2)
    {
      // Clamped-Clamped
      BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0,false,0 ); // unknown 0 - x
      BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0,false,1 ); // unknown 1 - y
      BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0,false,2 ); // unknown 2 - z

      BCs.addCondition(boundary::west, condition_type::clamped,0,0,false,2);

      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0,false,1 ); // unknown 1 - y
      BCs.addCondition(boundary::south, condition_type::clamped, 0, 0,false,2 ); // unknown 1 - y
      BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0,false,1 ); // unknown 1 - y
      BCs.addCondition(boundary::north, condition_type::clamped, 0, 0,false,2 ); // unknown 1 - y

      BCs.addCondition(boundary::east, condition_type::collapsed, 0, 0,false,0 ); // unknown 1 - y
      BCs.addCondition(boundary::east, condition_type::collapsed, 0, 0,false,1 ); // unknown 1 - y
      BCs.addCondition(boundary::east, condition_type::collapsed, 0, 0,false,2 ); // unknown 1 - y


      Load = 1e-3;
      gsVector<> point(2);
      gsVector<> load (3);
      point<< 1.0, 0.5 ;
      load << -Load, 0.0, 0.0 ;
      pLoads.addLoad(point, load, 0 );
    }
    else if (testCase == 3)
    {
        // poisson_ratio = 0.0;
        Load = 1e8;
        neu << Load, 0, 0;
        neuData.setValue(neu,3);
        // // Clamped-Clamped
        BCs.addCondition(boundary::east, condition_type::neumann, &neuData ); // unknown 0 - x

        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 1); // unknown 1 - y
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 2); // unknown 2 - z

        // BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false,0 ); // unknown 0 - x
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1); // unknown 1 - y
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 2); // unknown 2 - z

        // BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false,0 ); // unknown 0 - x
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1); // unknown 1 - y
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2); // unknown 2 - z

        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0); // unknown 0 - x
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false,1 ); // unknown 1 - y
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false,2 ); // unknown 2 - z
    }
    else if (testCase == 4)
    {
        PoissonRatio = 0.0;
        Load = 1e5;
        neu << Load, 0, 0;
        neuData.setValue(neu,3);
        // // Clamped-Clamped
        BCs.addCondition(boundary::east, condition_type::neumann, &neuData ); // unknown 0 - x
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false,1 ); // unknown 1 - y
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false,2 ); // unknown 2 - z
        BCs.addCondition(boundary::east, condition_type::clamped, 0, 0, false,2 ); // unknown 2 - z

        // BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false,1 ); // unknown 1 - y
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false,2 ); // unknown 2 - z
        BCs.addCondition(boundary::north, condition_type::clamped, 0, 0, false,2 ); // unknown 2 - z

        // BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false,1 ); // unknown 1 - y
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false,2 ); // unknown 2 - z
        BCs.addCondition(boundary::south, condition_type::clamped, 0, 0, false,2 ); // unknown 2 - z


        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false,0 ); // unknown 0 - x
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false,1 ); // unknown 1 - y
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false,2 ); // unknown 2 - z
        BCs.addCondition(boundary::west, condition_type::clamped, 0, 0, false,2 ); // unknown 2 - z
    }
    else if (testCase == 5)
    {
        PoissonRatio = 0.3;
        Load = 1e-2;
        neu << Load, 0, 0;
        neuData.setValue(neu,3);
        // // Clamped-Clamped
        BCs.addCondition(boundary::east, condition_type::neumann, &neuData ); // unknown 0 - x
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false,1 ); // unknown 1 - y
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false,2 ); // unknown 2 - z

        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false,1 ); // unknown 1 - y
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false,2 ); // unknown 2 - z

        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false,1 ); // unknown 1 - y
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false,2 ); // unknown 2 - z


        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false,0 ); // unknown 0 - x
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false,1 ); // unknown 1 - y
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false,2 ); // unknown 2 - z

        BCs.addCondition(boundary::north, condition_type::clamped,0,0,false,2);
        BCs.addCondition(boundary::east, condition_type::clamped,0,0,false,2);
        BCs.addCondition(boundary::south, condition_type::clamped,0,0,false,2);
        BCs.addCondition(boundary::west, condition_type::clamped,0,0,false,2);

    }
    else if (testCase == 6)
    {
        PoissonRatio = 0.;
        Load = 1e-8;
        neu << Load, 0, 0;
        neuData.setValue(neu,3);
        neu << -Load, 0, 0;
        neuData2.setValue(neu,3);
        // // Clamped-Clamped
        BCs.addCondition(boundary::east, condition_type::neumann, &neuData ); // unknown 0 - x
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0,false,2 ); // unknown 2 - z

        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0,false,1 ); // unknown 1 - y
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0,false,2 ); // unknown 2 - z

        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0,false,2 ); // unknown 2 - z


        BCs.addCondition(boundary::west, condition_type::neumann, &neuData2 ); // unknown 0 - x
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0,false,2 ); // unknown 2 - z
    }
    else if (testCase == 8)
    {
      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0 ,false,0);
      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0 ,false,1);
      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0 ,false,2);

      BCs.addCondition(boundary::north, condition_type::collapsed, 0, 0 ,false,0);
      BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0 ,false,1);
      BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0 ,false,2);

      Load = 1e0;
      gsVector<> point(2); point<< 1.0, 1.0 ;
      gsVector<> load (3); load << Load,0.0, 0.0;
      pLoads.addLoad(point, load, 0 );
    }
    else if (testCase == 16)
    {
      BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
      BCs.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 2 ); // unknown 2 - z

      BCs.addCondition(boundary::east, condition_type::collapsed, 0, 0, false, 0 ); // unknown 1 - y
      BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 2 - z.
      BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z.

      BCs.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 0 ); // unknown 2 - z.
      // BCs.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 1 ); // unknown 2 - z
      BCs.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 2 ); // unknown 2 - z
      // BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z.

      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 2 - z.

      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z.

      Load = 1e0;
      gsVector<> point(2);
      gsVector<> load (3);
      point<< 1.0, 0.5 ;
      load << Load,0.0, 0.0;
      pLoads.addLoad(point, load, 0 );

    }

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
    gsConstantFunction<> ratio(Ratio,3);

    real_t mu = E_modulus / (2 * (1 + PoissonRatio));
    gsConstantFunction<> alpha1(1.3,3);
    gsConstantFunction<> mu1(6.3e5/4.225e5*mu,3);
    gsConstantFunction<> alpha2(5.0,3);
    gsConstantFunction<> mu2(0.012e5/4.225e5*mu,3);
    gsConstantFunction<> alpha3(-2.0,3);
    gsConstantFunction<> mu3(-0.1e5/4.225e5*mu,3);

    index_t kmax = 1;

    std::vector<gsFunctionSet<> * > Gs(kmax);
    std::vector<gsFunctionSet<> * > Ts(kmax);
    std::vector<gsFunctionSet<> * > Phis(kmax);

    gsMatrix<> Gmat = gsCompositeMatrix(E_modulus,E_modulus,0.5 * E_modulus / (1+PoissonRatio),PoissonRatio,PoissonRatio);
    Gmat.resize(Gmat.rows()*Gmat.cols(),1);
    gsConstantFunction<> Gfun(Gmat,3);
    Gs[0] = &Gfun;

    gsConstantFunction<> phi;
    phi.setValue(0,3);

    Phis[0] = &phi;

    gsConstantFunction<> thicks(thickness/kmax,3);
    Ts[0] = &thicks;

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
    assembler = new gsThinShellAssembler<3, real_t, true >(mp,dbasis,BCs,force,materialMatrix);


    // Construct assembler object
    assembler->setOptions(opts);
    assembler->setPointLoads(pLoads);

    // Initialise solution object
    gsMultiPatch<> solution = mp;

    assembler->assemble();
    gsSparseMatrix<> K_L =  assembler->matrix();
    gsVector<> rhs = assembler->rhs();

    typedef std::function<gsSparseMatrix<real_t> (gsVector<real_t> const &)>    Jacobian_t;
    Jacobian_t K_NL = [&assembler,&mp_def](gsVector<real_t> const &x)
    {
      assembler->constructSolution(x,mp_def);
      assembler->assemble(mp_def);
      gsSparseMatrix<real_t> m = assembler->matrix();
      return m;
    };

      gsBucklingSolver<real_t,Spectra::GEigsMode::ShiftInvert> buckling(K_L,rhs,K_NL);
      buckling.verbose();
      // buckling.computePower();

      if (!sparse)
        buckling.compute();
      else
        buckling.computeSparse(shift,10,2,Spectra::SortRule::LargestMagn,Spectra::SortRule::SmallestMagn);

      gsMatrix<> values = buckling.values();
      gsMatrix<> vectors = buckling.vectors();

      gsDebugVar(buckling.vectors().cols());

    gsInfo<< "First 10 eigenvalues:\n";
    for (index_t k = 0; k<10; k++)
        gsInfo<<"\t"<<std::setprecision(20)<<values.at(k)<<"\n";
    gsInfo<<"\n";

    for (index_t k = 0; k<10; k++)
    {
      gsInfo<<"\t"<<values.at(k)*Load<<"\n";
    }

    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";
        system("mkdir -p BucklingResults");
        gsMultiPatch<> deformation = solution;
        gsMatrix<> modeShape;
        gsParaviewCollection collection("BucklingResults/modes");

        int N = 1;
        if (!first)
          N = vectors.cols();
        for (index_t m=0; m<N; m++)
        {

          // Compute solution based on eigenmode with number 'mode'
          modeShape = vectors.col(m);//solver.solve( assembler->rhs() );
          assembler->constructSolution(modeShape, solution);

          // compute the deformation spline
          deformation = solution;
          deformation.patch(0).coefs() -= mp.patch(0).coefs();// assuming 1 patch here

          // gsField<> solField(mp,deformation);
          // gsField<> mpField(mp,mp);
          //
          // real_t norm = solField.distanceL2(mpField);

          // Normalize mode shape amplitude in z coordinate
          real_t maxAmpl = std::max(math::abs(deformation.patch(0).coefs().col(2).maxCoeff()),math::abs(deformation.patch(0).coefs().col(2).minCoeff()));
          if (maxAmpl!=0.0)
          {
            deformation.patch(0).coefs() = deformation.patch(0).coefs()/maxAmpl;
          }

          gsField<> solField(mp,deformation);
          std::string fileName = "BucklingResults/modes" + util::to_string(m);
          gsWriteParaview<>(solField, fileName, 5000);
          fileName = "modes" + util::to_string(m) + "0";
          collection.addTimestep(fileName,m,".vts");
        }
        collection.save();
    }

    if (write)
    {
        system("mkdir -p BucklingResults");
        std::string wnM = "BucklingResults/eigenvalues.txt";
        writeToCSVfile(wnM,values);
    }

    return result;
}

template <class T>
gsMultiPatch<T> RectangularDomain(int n, int p, T L, T B)
{
  int q = p;
  int m = n;
  gsMultiPatch<T> mp = RectangularDomain(n, m, p, q, L, B);
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

  // Refine n times
  for(index_t i = 0; i< n; ++i)
      mp.patch(0).uniformRefine();
  // Elevate up to order p
  if (p>2)
  {
    for(index_t i = 2; i< p; ++i)
        mp.patch(0).degreeElevate();    // Elevate the degree
  }

  return mp;
}
