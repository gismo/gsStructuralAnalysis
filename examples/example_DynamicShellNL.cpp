/** @file example_DynamicShellNL.cpp

    @brief Example for nonlinear time integration of a nonlinear shell

    Fig 12 of:
    Filho, L. A. D., & Awruch, A. M. (2004).
    Geometrically nonlinear static and dynamic analysis of shells and plates using the eight-node hexahedral element with one-point quadrature.
    Finite Elements in Analysis and Design. https://doi.org/10.1016/j.finel.2003.08.012

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

#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicExplicitEuler.h>
#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicImplicitEuler.h>
#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicNewmark.h>
#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicBathe.h>
#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicWilson.h>
using namespace gismo;

// Choose among various shell examples, default = Thin Plate
#ifdef gsKLShell_ENABLED
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
    real_t tend = 3e-4;

    std::string assemberOptionsFile("options/solver_options.xml");

    gsCmdLine cmd("Dynamics of a nonlinear spherical cap.");
    cmd.addString( "f", "file", "Input XML file for assembler options", assemberOptionsFile );
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

    gsFileData<> fd(assemberOptionsFile);
    gsOptionList opts;
    fd.getFirst<gsOptionList>(opts);

    gsInfo<<"Simulation time:\t"<<tend<<"\n";
    real_t dt = tend/((double)steps);
    gsInfo<<"EI:\t"<<EI<<"\n";
    gsInfo<<"EA:\t"<<EA<<"\n";
    gsInfo<<"Density:\t"<<Density<<"\n";

    std::string dirname;
    dirname = "DynamicShellResults";
    gsFileManager::mkdir(dirname);



    wn = dirname + "/output.csv";

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

    gsWriteParaview<>(mp, "mp", 500);

//------------------------------------------------------------------------------
// Define Boundary conditions and Initial conditions
//------------------------------------------------------------------------------

    gsBoundaryConditions<> BCs;
    BCs.setGeoMap(mp);
    gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();

    BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
    BCs.addCondition(boundary::north, condition_type::clamped,0,0,false,2);

    BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
    BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
    BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
    BCs.addCondition(boundary::east, condition_type::clamped,0,0,false,2);

    BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
    BCs.addCondition(boundary::south, condition_type::clamped,0,0,false,2);

    BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 1 - y
    BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
    BCs.addCondition(boundary::west, condition_type::clamped,0,0,false,2);
    BCs.addCondition(boundary::west, condition_type::collapsed,0,0,false,2);

    gsVector<> point(2); point<< 0, 0 ; // Point where surface force is applied
    gsVector<> loadvec (3); loadvec << 0, 0, -25 ;
    pLoads.addLoad(point, loadvec, 0 );

    gsVector<> tmp(3);
    tmp<<0,0,0;

//------------------------------------------------------------------------------
// Define Beam Assembler and Assembler for L2-projection
//------------------------------------------------------------------------------

    gsMultiPatch<> mp_def = mp;
    gsMultiBasis<> dbasis(mp);
    // Linear isotropic material model
    gsFunctionExpr<> force("0","0","0",3);
    gsFunctionExpr<> t(std::to_string(thickness), 3);
    gsFunctionExpr<> E(std::to_string(E_modulus),3);
    gsFunctionExpr<> nu(std::to_string(PoissonRatio),3);
    gsFunctionExpr<> rho(std::to_string(Density),3);

    gsMaterialMatrixLinear<3,real_t> materialMatrix(mp,t,E,nu,rho);
    gsThinShellAssemblerBase<real_t>* assembler;
    assembler = new gsThinShellAssembler<3, real_t, true >(mp,dbasis,BCs,force,&materialMatrix);

    assembler->setOptions(opts);
    assembler->setPointLoads(pLoads);

    assembler->assemble();
    size_t N = assembler->numDofs();
    gsMatrix<> uOld(N,1);
    gsMatrix<> vOld(N,1);
    gsMatrix<> aOld(N,1);
    uOld.setZero();
    vOld.setZero();
    aOld.setZero();


//------------------------------------------------------------------------------
// Initiate mass and stiffness matrices and vectors for velocity, displacement and acceleration
//------------------------------------------------------------------------------

    gsSparseMatrix<> M, C;
    gsSparseMatrix<> K;
    gsSparseMatrix<> K_T;

//------------------------------------------------------------------------------
// Nonlinear time integration
//------------------------------------------------------------------------------
gsParaviewCollection collection(dirname + "/solution");

// Function for the Jacobian
gsStructuralAnalysisOps<real_t>::Jacobian_t Jacobian = [&assembler,&mp_def](gsMatrix<real_t> const &x, gsSparseMatrix<real_t> & m)
{
    // to do: add time dependency of forcing
    ThinShellAssemblerStatus status;
    assembler->constructSolution(x,mp_def);
    status = assembler->assembleMatrix(mp_def);
    m = assembler->matrix();
    return status == ThinShellAssemblerStatus::Success;
};

// Function for the Residual
gsStructuralAnalysisOps<real_t>::TResidual_t Residual = [&assembler,&mp_def](gsMatrix<real_t> const &x, real_t /*time*/, gsVector<real_t> & result)
{
    ThinShellAssemblerStatus status;
    assembler->constructSolution(x,mp_def);
    status = assembler->assembleVector(mp_def);
    result = assembler->rhs();
    return status == ThinShellAssemblerStatus::Success;
  };

// // Function for the Residual (TEST FO CHANGE FUNCTION!)
// Residual_t Residual = [&force,&assembler,&solution](gsMatrix<real_t> const &x, real_t time)
// {
//   gsFunctionExpr<> force2("0","0",std::to_string(time),3);
//   force.swap(force2);
//   assembler->constructSolution(x,solution);
//   assembler->assembleVector(solution);
//   return assembler->rhs();
// };

// Compute mass matrix (since it is constant over time)
assembler->assembleMass();
M = assembler->matrix();
C = gsSparseMatrix<>(assembler->numDofs(),assembler->numDofs());
// pre-assemble system
assembler->assemble();

gsStructuralAnalysisOps<real_t>::Mass_t    Mass    = [&M](                          gsSparseMatrix<real_t> & m) { m = M; return true; };
gsStructuralAnalysisOps<real_t>::Damping_t Damping = [&C](const gsVector<real_t> &, gsSparseMatrix<real_t> & m) { m = C; return true; };

gsDynamicBase<real_t> * timeIntegrator;
if (method==1)
    timeIntegrator = new gsDynamicExplicitEuler<real_t,true>(Mass,Damping,Jacobian,Residual);
else if (method==2)
    timeIntegrator = new gsDynamicImplicitEuler<real_t,true>(Mass,Damping,Jacobian,Residual);
else if (method==3)
    timeIntegrator = new gsDynamicNewmark<real_t,true>(Mass,Damping,Jacobian,Residual);
else if (method==4)
    timeIntegrator = new gsDynamicBathe<real_t,true>(Mass,Damping,Jacobian,Residual);
else if (method==5)
{
    timeIntegrator = new gsDynamicWilson<real_t,true>(Mass,Damping,Jacobian,Residual);
    timeIntegrator->options().setReal("gamma",1.4);
}
else
    GISMO_ERROR("Method "<<method<<" not known");

timeIntegrator->options().setReal("DT",dt);
timeIntegrator->options().setReal("TolU",1e-3);
timeIntegrator->options().setSwitch("Quasi",quasiNewton);
timeIntegrator->options().setSwitch("Verbose",true);

//------------------------------------------------------------------------------
// Initial Conditions
//------------------------------------------------------------------------------
gsMatrix<> uNew,vNew,aNew;
uNew = uOld;
vNew = vOld;
aNew = aOld;
timeIntegrator->setU(uNew);
timeIntegrator->setV(vNew);
timeIntegrator->setA(aNew);

real_t time;
for (index_t i=0; i<steps; i++)
{
  gsStatus status = timeIntegrator->step();
  if (status!=gsStatus::Success)
    GISMO_ERROR("Time integrator did not succeed");

  gsMatrix<> displacements = timeIntegrator->solutionU();

  assembler->constructSolution(displacements,mp_def);

  mp_def.patch(0).coefs() -= mp.patch(0).coefs();// assuming 1 patch here
  gsField<> solField(mp,mp_def);
  std::string fileName = dirname + "/solution" + util::to_string(i);
  gsWriteParaview<>(solField, fileName, 500);
  fileName = "solution" + util::to_string(i) + "0";
  collection.addPart(fileName + ".vts",i);

  if (write)
  {
    gsMatrix<> v(2,1);
    v<<  0.0,0.0;
    gsMatrix<> res2;
    mp_def.patch(0).eval_into(v,res2);
    time = timeIntegrator->time();
    std::ofstream file;
    file.open(wn,std::ofstream::out | std::ofstream::app);
    file  << std::setprecision(10)
          << time << ","
          << res2(2,0) <<"\n";
    file.close();
  }
  // Update solution with multipatch coefficients to generate geometry again
  mp_def.patch(0).coefs() += mp.patch(0).coefs();// assuming 1 patch here

  // gsInfo<<displacements.transpose()<<"\n";
}

collection.save();

delete assembler;
delete timeIntegrator;


return result;
}
#else//gsKLShell_ENABLED
int main(int argc, char *argv[])
{
    gsWarn<<"G+Smo is not compiled with the gsKLShell module.";
    return EXIT_FAILURE;
}
#endif