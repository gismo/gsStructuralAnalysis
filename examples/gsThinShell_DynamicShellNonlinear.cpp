/** @file gsThinShell_DynamicShellNonlinear.cpp

    @brief Example for nonlinear time integration of a nonlinear shell

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

/*
  To do:
  - Add initial displacement/velocity/acceleration by L2 projection
  - Add other test case
*/

#include <gismo.h>

#include <gsKLShell/gsThinShellAssembler.h>

#include <gsStructuralAnalysis/gsTimeIntegrator.h>

using namespace gismo;

// Choose among various shell examples, default = Thin Plate
int main (int argc, char** argv)
{
    // Input options
    int numElevate  = 1;
    int numHref     = 1;
    bool plot       = false;
    bool quasiNewton = false;

    real_t thickness = 0.01576;
    real_t width = 1; // Width of the strip is equal to 1.
    real_t Area = thickness*width;

    real_t E_modulus = 1e7;
    real_t PoissonRatio = 0.3;
    real_t Density = 0.000245;
    gsMultiPatch<> mp;

    real_t EA = E_modulus*Area;
    real_t EI = 1.0/12.0*(width*math::pow(thickness,3))*E_modulus;

    int testCase = 0;

    int method = 2; // 1: Explicit Euler, 2: Implicit Euler, 3: Newmark, 4: Bathe

    int result = 0;
    std::string fn("surface/sphericalCap.xml");

    bool write = false;
    std::string wn;

    int steps = 100;
    real_t tend = 1e-5;

    gsCmdLine cmd("Thin shell plate example.");
    cmd.addInt("r","hRefine",
               "Number of dyadic h-refinement (bisection) steps to perform before solving",
               numHref);
    cmd.addInt("t", "testcase",
                "Test case: 0: clamped-clamped, 1: pinned-pinned, 2: clamped-free",
               testCase);
    cmd.addInt("m", "method",
               "1: Explicit Euler, 2: Implicit Euler, 3: Newmark, 4: Bathe",
              method);
    cmd.addInt("s", "steps",
               "Number of time steps",
              steps);
    cmd.addReal("p", "endtime",
           "End time of simulation",
          tend);
//     cmd.addString("g","geometry","File containing Geometry (.xml, .axl, .txt)",fn);
    cmd.addInt("e","degreeElevation",
               "Number of degree elevation steps to perform on the Geometry's basis before solving",
               numElevate);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addSwitch("q","quasi", "Use quasi newton method", quasiNewton);
    cmd.addSwitch("write", "Write convergence data to file", write);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }


    gsInfo<<"Simulation time:\t"<<tend<<"\n";
    real_t dt = tend/((double)steps);
    gsInfo<<"EI:\t"<<EI<<"\n";
    gsInfo<<"EA:\t"<<EA<<"\n";
    gsInfo<<"Density:\t"<<Density<<"\n";

    std::string dirname;
    // dirname = "DynamicShellResults_NM_r" + std::to_string(numHref) + "e" + std::to_string(numElevate) + "_dt" + std::to_string(dt);
    dirname = "DynamicShellResults";
    system(("mkdir -p " + dirname).c_str());

    wn = dirname + "/output.csv";

    std::string methodName;
    if (method==1)
      methodName = "ExplEuler";
    else if (method==2)
      methodName = "ImplEuler";
    else if (method==3)
      methodName = "Newmark";
    else if (method==4)
      methodName = "Bathe";


    if (write)
    {
      std::ofstream file;
      file.open(wn,std::ofstream::out);
      file  << std::setprecision(20)
            << "time" << ","
            << "displacement"
            << "\n";
      file.close();
    }

//------------------------------------------------------------------------------
// Make domain and define basis
//------------------------------------------------------------------------------
    gsReadFile<>(fn, mp);

    for(index_t i = 0; i< numElevate; ++i)
        mp.patch(0).degreeElevate();    // Elevate the degree

    // h-refine
    for(index_t i = 0; i< numHref; ++i)
        mp.patch(0).uniformRefine();

    gsInfo<<"Basis (patch 0): "<< mp.patch(0).basis() << "\n";

//------------------------------------------------------------------------------
// Define Boundary conditions and Initial conditions
//------------------------------------------------------------------------------

    gsBoundaryConditions<> BCs;
    gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();

    BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
    BCs.addCondition(boundary::north, condition_type::clamped,0,0,false,2);

    BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
    BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
    BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
    BCs.addCondition(boundary::east, condition_type::clamped,0,0,false,2);

    BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
    BCs.addCondition(boundary::south, condition_type::clamped,0,0,false,2);

    BCs.addCondition(boundary::west, condition_type::clamped,0,0,false,2);

    gsVector<> point(2); point<< 0, 0 ; // Point where surface force is applied
    gsVector<> loadvec (3); loadvec << 0, 0, -25 ;
    pLoads.addLoad(point, loadvec, 0 );

    gsVector<> tmp(3);
    tmp<<0,0,0;

//------------------------------------------------------------------------------
// Define Beam Assembler and Assembler for L2-projection
//------------------------------------------------------------------------------

    gsMultiPatch<> mp_def = mp, solution = mp;
    gsMultiBasis<> dbasis(mp);
    // Linear isotropic material model
    gsConstantFunction<> force(tmp,3);
    gsFunctionExpr<> thick(std::to_string(thickness), 3);
    gsFunctionExpr<> E(std::to_string(E_modulus),3);
    gsFunctionExpr<> nu(std::to_string(PoissonRatio),3);
    gsFunctionExpr<> rho(std::to_string(Density),3);

    std::vector<gsFunction<>*> parameters(2);
    parameters[0] = &E;
    parameters[1] = &nu;
    gsMaterialMatrix materialMat(mp,mp_def,thick,parameters,rho);

    // Construct assembler object for dynamic computations
    gsThinShellAssembler assembler(mp,dbasis,BCs,force,materialMat);
    assembler.setPointLoads(pLoads);

    assembler.assemble();
    size_t N = assembler.numDofs();
    gsMatrix<> uOld(N,1);
    gsMatrix<> vOld(N,1);
    gsMatrix<> aOld(N,1);
    uOld.setZero();
    vOld.setZero();
    aOld.setZero();


//------------------------------------------------------------------------------
// Initiate mass and stiffness matrices and vectors for velocity, displacement and acceleration
//------------------------------------------------------------------------------

    gsMatrix<real_t> Minv;
    gsSparseMatrix<> M;
    gsSparseMatrix<> K;
    gsSparseMatrix<> K_T;

//------------------------------------------------------------------------------
// Nonlinear time integration
//------------------------------------------------------------------------------
gsParaviewCollection collection(dirname + "/solution");

typedef std::function<gsSparseMatrix<real_t> (gsVector<real_t> const &)>    Jacobian_t;
typedef std::function<gsVector<real_t> (gsVector<real_t> const &) >         Residual_t;
// Function for the Jacobian
Jacobian_t Jacobian = [&assembler,&solution](gsMatrix<real_t> const &x)
{
  // to do: add time dependency of forcing
  assembler.constructSolution(x,solution);
  assembler.assemble(solution);
  gsSparseMatrix<real_t> m = assembler.matrix();
  return m;
};

// Function for the Residual
Residual Residual = [&assembler,&solution](gsMatrix<real_t> const &x, real_t time)
{
  assembler.constructSolution(x,solution);
  assembler.assemble(solution);
  gsMatrix<real_t> r = assembler.rhs();
  return r;
};

// Compute mass matrix (since it is constant over time)
assembler.assembleMass();
M = assembler.matrix();
// pre-assemble system
assembler.assemble();

// // set damping Matrix (same dimensions as M)
// C.setZero(M.rows(),M.cols());

gsTimeIntegrator<real_t> timeIntegrator(M,Jacobian,Residual,dt);

timeIntegrator.verbose();
timeIntegrator.setTolerance(1e-6);
timeIntegrator.setMethod(methodName);
if (quasiNewton)
  timeIntegrator.quasiNewton();

//------------------------------------------------------------------------------
// Initial Conditions
//------------------------------------------------------------------------------
gsMatrix<> uNew,vNew,aNew;
uNew = uOld;
vNew = vOld;
aNew = aOld;
timeIntegrator.setDisplacement(uNew);
timeIntegrator.setVelocity(vNew);
timeIntegrator.setAcceleration(aNew);


for (index_t i=0; i<steps; i++)
{
  real_t time = timeIntegrator.currentTime();
  timeIntegrator.step();
  timeIntegrator.constructSolution();
  gsMatrix<> displacements = timeIntegrator.displacements();

  assembler.constructSolution(displacements,solution);

  solution.patch(0).coefs() -= mp.patch(0).coefs();// assuming 1 patch here
  gsField<> solField(mp,solution);
  std::string fileName = dirname + "/solution" + util::to_string(i);
  gsWriteParaview<>(solField, fileName, 500);
  fileName = "solution" + util::to_string(i) + "0";
  collection.addTimestep(fileName,i,".vts");

  if (write)
  {
    gsMatrix<> v(2,1);
    v<<  0.0,0.0;
    gsMatrix<> res2;
    solution.patch(0).eval_into(v,res2);

    std::ofstream file;
    file.open(wn,std::ofstream::out | std::ofstream::app);
    file  << std::setprecision(10)
          << time << ","
          << res2(2,0) <<"\n";
    file.close();
  }
  // Update solution with multipatch coefficients to generate geometry again
  solution.patch(0).coefs() += mp.patch(0).coefs();// assuming 1 patch here

  // gsInfo<<displacements.transpose()<<"\n";
}

collection.save();
return result;
}