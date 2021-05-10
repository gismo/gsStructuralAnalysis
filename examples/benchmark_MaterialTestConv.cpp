/** @file benchmark_MaterialTestConv.cpp

    @brief Stretches a planar material with both ends clamped; convergence test

    See fig 1 of Roohbakhshan and Sauer 2017
    Roohbakhshan, F., & Sauer, R. A. (2017). Efficient isogeometric thin shell formulations for soft biological materials. Biomechanics and Modeling in Mechanobiology, 16(5), 1569–1597. https://doi.org/10.1007/s10237-017-0906-6

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#include <gismo.h>

#include <gsKLShell/gsThinShellAssembler.h>
#include <gsKLShell/getMaterialMatrix.h>

// #include <gsThinShell/gsArcLengthIterator.h>
#include <gsStructuralAnalysis/gsStaticSolver.h>

using namespace gismo;

template <class T>
gsMultiPatch<T> Rectangle(T L, T B);

int main (int argc, char** argv)
{
    // Input options
    int numElevate  = -1;
    int numHref     = -1;
    bool plot       = false;
    bool last       = false;

    index_t Compressibility = 0;
    index_t material = 1;
    bool composite = false;
    index_t impl = 1; // 1= analytical, 2= generalized, 3= spectral

    int testCase = 0;

    int result = 0;

    index_t maxit = 50;

    real_t tol = 1e-6;

    std::string wn("data.csv");

    gsCmdLine cmd("Wrinkling analysis with thin shells.");

    cmd.addInt("t", "testcase", "Test case: 0: clamped-clamped, 1: pinned-pinned, 2: clamped-free", testCase);

    cmd.addInt("r","hRefine", "Number of dyadic h-refinement (bisection) steps to perform before solving", numHref);
    cmd.addInt("e","degreeElevation", "Number of degree elevation steps to perform on the Geometry's basis before solving", numElevate);

    cmd.addInt( "M", "Material", "Material law",  material );
    cmd.addInt( "c", "Compressibility", "1: compressible, 0: incompressible",  Compressibility );
    cmd.addInt( "I", "Implementation", "Implementation: 1= analytical, 2= generalized, 3= spectral",  impl );
    cmd.addSwitch("composite", "Composite material", composite);

    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addSwitch("last", "Only last refinement step", last);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    real_t aDim,bDim;
    real_t thickness;
    real_t E_modulus;
    real_t PoissonRatio;
    real_t Density = 1e0;
    real_t Ratio;
    real_t mu;

    real_t A1,A2,A3;
    real_t M1,M2,M3;

    real_t Load;

    /*
      Data from Roohbakhshan2017
    */
    aDim = 9e-3;
    bDim = 3e-3;
    thickness = 0.3e-3;

    mu = 10e3;
    if ((!Compressibility) && (material!=0))
      PoissonRatio = 0.5;
    else
      PoissonRatio = 0.45;

    if (material==1)
    {
      mu = 10e3;
      if (material==4)
      {
        A1 = 1.3;
        A2 = 5.0;
        A3 = -2.0;
        M1 = 6.3e5/4.225e5*mu;
        M2 = 0.012e5/4.225e5*mu;
        M3 = -0.1e5/4.225e5*mu;
      }
    }
    else if (material==3)
    {
      mu = 30e3;
    }

    E_modulus = 2*mu*(1+PoissonRatio);

    real_t EA = 2*mu*(1+PoissonRatio) * bDim*thickness;
    gsInfo<<"EA = "<<EA<<"\n";

    Ratio = 0.5;
    Load = 0.5*EA;

    gsMultiPatch<> mp,mp_def;

    mp = Rectangle(aDim   , bDim   );

    for(index_t i = 0; i< numElevate; ++i)
      mp.patch(0).degreeElevate();    // Elevate the degree

    if (last)
    {
      // h-refine
      for(index_t i = 0; i< numHref; ++i)
        mp.patch(0).uniformRefine();
      numHref = 0;
    }


    mp_def = mp;
    gsInfo<<"alpha = "<<aDim/bDim<<"; beta = "<<bDim/thickness<<"\n";


    gsMultiBasis<> dbasis(mp);
    gsInfo<<"Basis (patch 0): "<< mp.patch(0).basis() << "\n";

    // Boundary conditions
    gsBoundaryConditions<> BCs;
    gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();

    BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 2 - z
    BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 2 - z

    BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0 ,false,1);
    BCs.addCondition(boundary::east, condition_type::collapsed, 0, 0 ,false,0);

    gsVector<> point(2); point<< 1.0, 0.5 ;
    gsVector<> load (2); load << Load,0.0;
    pLoads.addLoad(point, load, 0 );

    // plot geometry
    if (plot)
      gsWriteParaview(mp,"mp",1000,true);

    // Linear isotropic material model
    gsFunctionExpr<> force("0","0",2);
    gsConstantFunction<> t(thickness,2);
    gsConstantFunction<> E(E_modulus,2);
    gsConstantFunction<> nu(PoissonRatio,2);
    gsConstantFunction<> rho(Density,2);
    gsConstantFunction<> ratio(Ratio,2);

    mu = E_modulus / (2 * (1 + PoissonRatio));
    gsConstantFunction<> alpha1(A1,2);
    gsConstantFunction<> mu1(M1,2);
    gsConstantFunction<> alpha2(A2,2);
    gsConstantFunction<> mu2(M2,2);
    gsConstantFunction<> alpha3(A3,2);
    gsConstantFunction<> mu3(M3,2);

    index_t kmax = 1;

    std::vector<gsFunctionSet<> * > Gs(kmax);
    std::vector<gsFunctionSet<> * > Ts(kmax);
    std::vector<gsFunctionSet<> * > Phis(kmax);

    gsMatrix<> Gmat = gsCompositeMatrix(E_modulus,E_modulus,0.5 * E_modulus / (1+PoissonRatio),PoissonRatio,PoissonRatio);
    Gmat.resize(Gmat.rows()*Gmat.cols(),1);
    gsConstantFunction<> Gfun(Gmat,2);
    Gs[0] = &Gfun;

    gsConstantFunction<> phi;
    phi.setValue(0,2);

    Phis[0] = &phi;

    gsConstantFunction<> thicks(thickness/kmax,2);
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
            materialMatrix = new gsMaterialMatrixComposite<2,real_t>(mp,Ts,Gs,Phis);
        }
        else
        {
            parameters.resize(2);
            options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
            options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
            materialMatrix = getMaterialMatrix<2,real_t>(mp,t,parameters,rho,options);
        }
    }
    else
    {
        options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",material);
        options.addSwitch("Compressibility","Compressibility: (false): Imcompressible | (true): Compressible",Compressibility);
        options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",impl);
        materialMatrix = getMaterialMatrix<2,real_t>(mp,t,parameters,rho,options);
    }

    gsThinShellAssemblerBase<real_t>* assembler;

    gsStopwatch stopwatch;
    real_t time = 0.0;

    typedef std::function<gsSparseMatrix<real_t> (gsVector<real_t> const &)>    Jacobian_t;
    typedef std::function<gsVector<real_t> (gsVector<real_t> const &) >         Residual_t;
    // Function for the Jacobian
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


    gsSparseMatrix<> matrix;
    gsVector<> vector, solVector;
    std::vector<real_t> strains(numHref+1);
    for (index_t k=0; k!=numHref+1; k++)
    {
      assembler = new gsThinShellAssembler<2,real_t,false>(mp,dbasis,BCs,force,materialMatrix);
      assembler->setPointLoads(pLoads);
      gsInfo<<"DoFs = "<<assembler->numDofs()<<"\n";

      assembler->assemble();
      matrix = assembler->matrix();
      vector = assembler->rhs();
      // Configure Structural Analsysis module
      gsStaticSolver<real_t> staticSolver(matrix,vector,Jacobian,Residual);
      gsOptionList solverOptions = staticSolver.options();
      solverOptions.setInt("Verbose",true);
      solverOptions.setInt("MaxIterations",maxit);
      solverOptions.setReal("Tolerance",tol);
      staticSolver.setOptions(solverOptions);
      solVector = staticSolver.solveNonlinear();

      mp_def = assembler->constructSolution(solVector);
      gsMultiPatch<> deformation = mp_def;
      for (size_t k = 0; k != mp_def.nPatches(); ++k)
          deformation.patch(k).coefs() -= mp.patch(k).coefs();

      gsVector<> u(2);
      u<<1,0.5;
      gsMatrix<> def = mp_def.patch(0).eval(u);
      strains.at(k) = (def.at(0)-aDim)/aDim;

      mp.uniformRefine();
      mp_def = mp;
      dbasis.uniformRefine();

    }

    for (size_t k=0; k!=strains.size(); k++)
      gsInfo<< std::setprecision(12) <<strains.at(k)<<",";
    gsInfo<<"\n";


  return result;
}

template <class T>
gsMultiPatch<T> Rectangle(T L, T B) //, int n, int m, std::vector<boxSide> sides, T offset)
{
  // -------------------------------------------------------------------------
  // --------------------------Make beam geometry-----------------------------
  // -------------------------------------------------------------------------
  int dim = 2; //physical dimension
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

