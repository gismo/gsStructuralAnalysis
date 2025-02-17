/** @file gsThinShell_DynamicBeamNL.cpp

    @brief Computes nonlinear dynamic analysis of a beam

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
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
#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicRK4.h>

using namespace gismo;

template <class T>
gsMultiPatch<T> RectangularDomain(int n, int m, int p, int q, T L, T B);
template <class T>
gsMultiPatch<T> RectangularDomain(int n, int p, T L, T B);


// Choose among various shell examples, default = Thin Plate
#ifdef gsKLShell_ENABLED
int main (int argc, char** argv)
{
    // Input options
    int numElevate  = 1;
    int numHref     = 1;
    bool plot       = false;
    bool nonlinear  = false;

    real_t thickness     = 1.0;
    real_t width = 1; // Width of the strip is equal to 1.
    real_t length = 10; // Length of the strip is equal to 10.
    real_t Area = thickness*width;

    real_t E_modulus     = 1e0;
    real_t PoissonRatio = 0;
    real_t Density = 1e0;
    real_t load = 1e-3;
    gsMultiPatch<> mp;

    real_t EA = E_modulus*Area;
    real_t EI = 1.0/12.0*(width*math::pow(thickness*2,3))*E_modulus;

    int method = 2; // 1: Explicit Euler, 2: Implicit Euler, 3: Implicit Euler (gsBlockOp)

    int result = 0;
    std::string fn("planar/strip.xml");

    bool write = false;

    int steps = 100;
    real_t periods = 1;

    gsCmdLine cmd("Dynamics of a nonlinear beam.");
    cmd.addInt("r","hRefine",
               "Number of dyadic h-refinement (bisection) steps to perform before solving",
               numHref);
    cmd.addInt("m", "method",
               "1: Explicit Euler, 2: Implicit Euler, 3: Implicit Euler (gsBlockOp)",
              method);
    cmd.addInt("s", "steps",
               "Number of time steps",
              steps);
    cmd.addReal("p", "periods",
           "Number of periods",
          periods);
//     cmd.addString("g","geometry","File containing Geometry (.xml, .axl, .txt)",fn);
    cmd.addInt("e","degreeElevation",
               "Number of degree elevation steps to perform on the Geometry's basis before solving",
               numElevate);
    cmd.addSwitch("nl", "Nonlinear elasticity (otherwise linear)", nonlinear);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addSwitch("write", "Write convergence data to file", write);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    real_t omega = math::pow((3.1415926535/length),2)*math::pow(EI/(Density*Area),0.5); // vibration frequency (man. solution)

    omega = 1;
    gsInfo<<"eigenfrequency:\t"<<omega<<"\n";
    // real_t tend = 2*3.1415926535/omega*periods;
    real_t tend = 2*1/omega*periods;
    gsInfo<<"Simulation time:\t"<<tend<<"\n";
    real_t dt = tend/((double)steps);
    gsInfo<<"EI:\t"<<EI<<"\n";
    gsInfo<<"EA:\t"<<EA<<"\n";
    gsInfo<<"rho:\t"<<Density<<"\n";

//------------------------------------------------------------------------------
// Make domain and define basis
//------------------------------------------------------------------------------

    // -------------------------------------------------------------------------
    // --------------------------Make beam geometry-----------------------------
    // -------------------------------------------------------------------------
    mp = RectangularDomain(numHref, 0, numElevate, 1, length,width);

    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------

    gsInfo<<"Basis (patch 0): "<< mp.patch(0).basis() << "\n";

//------------------------------------------------------------------------------
// Define Boundary conditions and Initial conditions
//------------------------------------------------------------------------------

    // Boundary conditions
    gsBoundaryConditions<> BCs;
    BCs.setGeoMap(mp);

    // BCs.addCondition(boundary::east, condition_type::neumann,&neuData );
    // Left side is always restrained
    BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false,0 ); // unknown 0 - x
    BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
    BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

    // Initiate Surface forces
    std::string uinix("0");
    std::string uiniy("0");
    std::string uiniz("0");

    std::string ainix("0");
    std::string ainiy("0");
    std::string ainiz("0");

    // Pinned-Pinned
    BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
    BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
    BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

    real_t time = 0;
    char buffer_u_ini[200];
    sprintf(buffer_u_ini,"%e*x/(24*%e)*(%e^3-2*%e*x^2+x^3)*cos(%e*pi*%e)",load,EI,length,length,omega,time);
    uiniz = buffer_u_ini;

    // The initial velocity is zero
    // char buffer_ini[200];
    // sprintf(buffer_ini,"-1*x/(24*%e)*(%e^3-2*%e*x^2+x^3)*%e*pi*cos(%e*pi*%e)",EI,length,length,omega,omega,t)
    // viniz = buffer_ini;

    // Expression for initial acceleration
    char buffer_a_ini[200];
    sprintf(buffer_a_ini,"-%e*x/(24*%e)*(%e^3-2*%e*x^2+x^3)*(%e*pi)^2*cos(%e*pi*%e)",load,EI,length,length,omega,omega,time);
    ainiz = buffer_a_ini;

//------------------------------------------------------------------------------
// Define Beam Assembler and Assembler for L2-projection
//------------------------------------------------------------------------------

    gsMultiPatch<> mp_def = mp;
    gsMultiBasis<> dbasis(mp);
    // Linear isotropic material model
    gsFunctionExpr<> force("0","0","0",3);
    gsFunctionExpr<> initialShape(uinix,uiniy,uiniz,3);
    gsFunctionExpr<> initialAcc(ainix,ainiy,ainiz,3);

    gsFunctionExpr<> t(std::to_string(thickness), 3);
    gsFunctionExpr<> E(std::to_string(E_modulus),3);
    gsFunctionExpr<> nu(std::to_string(PoissonRatio),3);
    gsFunctionExpr<> rho(std::to_string(Density),3);

    std::vector<gsFunctionSet<>*> parameters(2);
    parameters[0] = &E;
    parameters[1] = &nu;

    gsMaterialMatrixBase<real_t>::uPtr materialMatrix;

    gsOptionList options;
    options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
    options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
    materialMatrix = getMaterialMatrix<3,real_t>(mp,t,parameters,rho,options);

    gsThinShellAssemblerBase<real_t>* assembler;
    assembler = new gsThinShellAssembler<3, real_t, true >(mp,dbasis,BCs,force,materialMatrix);

    // assembler->setOptions(opts);


//------------------------------------------------------------------------------
// Compute initial deformation and set initial Velocity to zero
//------------------------------------------------------------------------------

    size_t N = assembler->numDofs();
    gsVector<> U(N);
    gsVector<> V(N);
    gsVector<> A(N);
    // U.setZero();
    // V.setZero();
    // A.setZero();

    U = assembler->projectL2(initialShape);
    // Set initial velocity to zero
    V.setZero();
    A = assembler->projectL2(initialAcc);

    // Compute Error
    gsMultiPatch<> deformationIni = mp;
    assembler->constructSolution(U,deformationIni);
    deformationIni.patch(0).coefs() -= mp.patch(0).coefs();
    gsField<> inifield(mp,deformationIni);
    gsInfo<<"Initial error:\t"<<inifield.distanceL2(initialShape)<<"\n";

    gsWriteParaview<>( inifield, "Beam_Numerical", 1000);

    gsField<> u_anfield(mp,initialShape);
    gsWriteParaview( u_anfield, "Beam_Analytical", 1000);

//------------------------------------------------------------------------------
// Initiate mass and stiffness matrices
//------------------------------------------------------------------------------

    gsMatrix<real_t> Minv;
    gsSparseMatrix<> M;
    gsSparseMatrix<> K;
    gsSparseMatrix<> K_T;
    gsMatrix<> F;

//------------------------------------------------------------------------------
// Nonlinear time integration
//------------------------------------------------------------------------------
std::string dirname = "DynamicBeamResults";
gsFileManager::mkdir(dirname);

std::string wn = dirname + "/output.csv";

gsParaviewCollection collection("DynamicBeamResults/solution");
gsParaviewCollection collection_an("DynamicBeamResults/analytical");

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

// Function for the Residual (TEST TO CHANGE FUNCTION!)
gsStructuralAnalysisOps<real_t>::TResidual_t Residual = [&EA,&load,&omega,&length,&EI,&force,&Density,&Area,&assembler,&mp_def](gsMatrix<real_t> const &x, real_t time, gsVector<real_t> & result)
{
  /// Make time dependent forcing
  ThinShellAssemblerStatus status;
  char buffer_traction[200];
  sprintf(buffer_traction,"%e*%e^2*cos(%e*pi*%e)^2*(%e^3-6*%e*x^2+4*x^3)*x*(%e-x)/(48*%e^2)",EA,load,omega,time,length,length,length,EI);
  std::string traction = buffer_traction;

  char buffer_pressure[400];

  sprintf(buffer_pressure,"-%e*(-(1/32)*%e*%e^2*x*(%e-x)*(%e-2*x)^2*(%e^2+2*%e*x-2*x^2)^2*cos(%e*pi*%e)^2+%e^2*(%e*%e^3*pi^2*%e^2*x-2*%e*%e*pi^2*%e^2*x^3+%e*pi^2*%e^2*x^4-24*%e))*cos(%e*pi*%e)/(24*%e^3)"
  ,load,EA,load,length,length,length,length,omega,time,EI,Density*Area,length,omega,Density*Area,length,omega,Density*Area,omega,EI,omega,time,EI);
  std::string pressure = buffer_pressure;

  gsFunctionExpr<> surfForceTemp(traction,"0",pressure,3);
  force.swap(surfForceTemp);
  assembler->constructSolution(x,mp_def);
  status = assembler->assembleVector(mp_def);
  result = assembler->rhs();
  return status == ThinShellAssemblerStatus::Success;
};

// Compute mass matrix (since it is constant over time)
assembler->assembleMass();
M = assembler->matrix();
// pre-assemble system
assembler->assemble();

gsSparseMatrix<> C = gsSparseMatrix<>(assembler->numDofs(),assembler->numDofs());
gsStructuralAnalysisOps<real_t>::Damping_t Damping = [&C](const gsVector<real_t> &, gsSparseMatrix<real_t> & m) { m = C; return true; };
gsStructuralAnalysisOps<real_t>::Mass_t    Mass    = [&M](                          gsSparseMatrix<real_t> & m) { m = M; return true; };

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
else if (method==6)
    timeIntegrator = new gsDynamicRK4<real_t,true>(Mass,Damping,Jacobian,Residual);
else
    GISMO_ERROR("Method "<<method<<" not known");

timeIntegrator->options().setReal("DT",dt);
timeIntegrator->options().setReal("TolU",1e-3);
timeIntegrator->options().setSwitch("Verbose",true);

timeIntegrator->setU(U);
timeIntegrator->setV(V);
timeIntegrator->setA(A);

//------------------------------------------------------------------------------
// Initial Conditions
//------------------------------------------------------------------------------
for (index_t i=0; i<steps; i++)
{
  gsStatus status = timeIntegrator->step(time,dt,U,V,A);
  GISMO_ASSERT(status == gsStatus::Success,"Time integrator did not succeed");

  assembler->constructSolution(U,mp_def);
  time = timeIntegrator->time();

  // Define manufactured solution
  char buffer_w_an[200];
  sprintf(buffer_w_an,"%e*x/(24*%e)*(%e^3-2*%e*x^2+x^3)*cos(%e*pi*%e)",load,EI,length,length,omega,time);
  std::string w_an = buffer_w_an;
  gsFunctionExpr<> analytical("0","0",w_an,3);

  // Update solution and export
  assembler->constructSolution(U,mp_def);
  mp_def.patch(0).coefs() -= mp.patch(0).coefs();// assuming 1 patch here
  gsField<> solField(mp,mp_def);
  std::string fileName = dirname + "/solution" + util::to_string(i);
  gsWriteParaview<>(solField, fileName, 500);
  fileName = "solution" + util::to_string(i) + "0";
  collection.addPart(fileName + ".vts",i);

  gsField<> anField(mp,analytical);
  fileName = dirname + "/analytical" + util::to_string(i);
  gsWriteParaview<>(anField, fileName, 500);
  fileName = "analytical" + util::to_string(i) + "0";
  collection_an.addPart(fileName + ".vts",i);

  if (write)
  {
    gsMatrix<> v(2,1);
    v<<  0.0,0.0;
    gsMatrix<> res2;
    mp_def.patch(0).eval_into(v,res2);
    std::ofstream file;
    file.open(wn,std::ofstream::out | std::ofstream::app);
    file  << std::setprecision(10)
          << time << ","
          << res2(2,0) <<"\n";
    file.close();
  }

  if (i==steps-1)
  {
    // Compute Error
    gsField<> anField(mp,analytical);
    gsInfo<<"Error at t = "<<t<<":\t"<<solField.distanceL2(anField)<<"\n";
    gsMatrix<> u(3,1);
    u << 5,0.5,0;
    gsMatrix<> v(2,1);
    v<<  0.5,0.5;
    gsMatrix<> res1;
    gsMatrix<> res2;
    analytical.eval_into(u,res1);
    // gsInfo<<"Analytical: \n"<<res1<<"\n";
    mp_def.patch(0).eval_into(v,res2);
    // gsInfo<<"Numerical: "<<res2<<"\n";

    gsInfo<<"Mid-point error at t = "<<t<<":\t"<<math::abs(res1(2,0)-res2(2,0))<<"\n";

    gsInfo<<"L2-norm error at t = "<<t<<":\t"<<solField.distanceL2(anField)<<"\n";
    //gsInfo <<"Deformation norm       : "<< deformation.patch(0).coefs().norm() <<".\n";
    gsInfo <<"Maximum deformation coef: "
           << mp_def.patch(0).coefs().colwise().maxCoeff() <<".\n";
    gsInfo <<"Minimum deformation coef: "
           << mp_def.patch(0).coefs().colwise().minCoeff() <<".\n";
  }
  // Update solution with multipatch coefficients to generate geometry again
  mp_def.patch(0).coefs() += mp.patch(0).coefs();// assuming 1 patch here

  // gsInfo<<displacements.transpose()<<"\n";
}
collection.save();
collection_an.save();

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
