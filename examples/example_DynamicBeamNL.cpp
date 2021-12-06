/** @file gsThinShell_DynamicBeam.cpp

    @brief Example testing thin shell solver.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gismo.h>

#include <gsKLShell/gsThinShellAssembler.h>
#include <gsKLShell/getMaterialMatrix.h>
#include <gsStructuralAnalysis/gsTimeIntegrator.h>

using namespace gismo;

template <class T>
gsMultiPatch<T> RectangularDomain(int n, int m, int p, int q, T L, T B);
template <class T>
gsMultiPatch<T> RectangularDomain(int n, int p, T L, T B);


// Choose among various shell examples, default = Thin Plate
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

    gsCmdLine cmd("Thin shell plate example.");
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

    std::string methodName;
    if (method==1)
      methodName = "ExplEuler";
    else if (method==2)
      methodName = "ImplEuler";
    else if (method==3)
      methodName = "Newmark";
    else if (method==4)
      methodName = "Bathe";
    else if (method==5)
      methodName = "CentralDiff";
    else if (method==6)
      methodName = "RK4";


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

    std::vector<gsFunction<>*> parameters(2);
    parameters[0] = &E;
    parameters[1] = &nu;

    gsMaterialMatrixBase<real_t>* materialMatrix;

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
    gsMatrix<> uOld(N,1);
    gsMatrix<> vOld(N,1);
    gsMatrix<> aOld(N,1);

    uOld = assembler->projectL2(initialShape);
    // Set initial velocity to zero
    vOld.setZero();
    aOld = assembler->projectL2(initialAcc);

    // Compute Error
    gsMultiPatch<> deformationIni = mp;
    assembler->constructSolution(uOld,deformationIni);
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
system(("mkdir -p " + dirname).c_str());
std::string wn = dirname + "/output.csv";

gsParaviewCollection collection("DynamicBeamResults/solution");
gsParaviewCollection collection_an("DynamicBeamResults/analytical");

typedef std::function<gsSparseMatrix<real_t> (gsVector<real_t> const &)>        Jacobian_t;
typedef std::function<gsVector<real_t> (gsVector<real_t> const &, real_t time) >Residual_t;
// Function for the Jacobian
Jacobian_t Jacobian = [&assembler,&mp_def](gsMatrix<real_t> const &x)
{
  // to do: add time dependency of forcing
  assembler->constructSolution(x,mp_def);
  assembler->assembleMatrix(mp_def);
  gsSparseMatrix<real_t> m = assembler->matrix();
  return m;
};

// Function for the Residual (TEST TO CHANGE FUNCTION!)
Residual_t Residual = [&EA,&load,&omega,&length,&EI,&force,&Density,&Area,&assembler,&mp_def](gsMatrix<real_t> const &x, real_t time)
{
  /// Make time dependent forcing
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
  assembler->assembleVector(mp_def);
  return assembler->rhs();
};

// Compute mass matrix (since it is constant over time)
assembler->assembleMass();
M = assembler->matrix();
// pre-assemble system
assembler->assemble();

// // set damping Matrix (same dimensions as M)
// C.setZero(M.rows(),M.cols());

gsTimeIntegrator<real_t> timeIntegrator(M,Jacobian,Residual,dt);

timeIntegrator.verbose();
timeIntegrator.setTolerance(1e-6);
timeIntegrator.setMethod(methodName);
// if (quasiNewton)
//   timeIntegrator.quasiNewton();

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
  timeIntegrator.step();
  timeIntegrator.constructSolution();
  uNew = timeIntegrator.displacements();

  assembler->constructSolution(uNew,mp_def);

  time = timeIntegrator.currentTime();
  // Define manufactured solution
  char buffer_w_an[200];
  sprintf(buffer_w_an,"%e*x/(24*%e)*(%e^3-2*%e*x^2+x^3)*cos(%e*pi*%e)",load,EI,length,length,omega,time);
  std::string w_an = buffer_w_an;
  gsFunctionExpr<> analytical("0","0",w_an,3);

  // Update solution and export
  assembler->constructSolution(uNew,mp_def);
  mp_def.patch(0).coefs() -= mp.patch(0).coefs();// assuming 1 patch here
  gsField<> solField(mp,mp_def);
  std::string fileName = dirname + "/solution" + util::to_string(i);
  gsWriteParaview<>(solField, fileName, 500);
  fileName = "solution" + util::to_string(i) + "0";
  collection.addTimestep(fileName,i,".vts");

  gsField<> anField(mp,analytical);
  fileName = dirname + "/analytical" + util::to_string(i);
  gsWriteParaview<>(anField, fileName, 500);
  fileName = "analytical" + util::to_string(i) + "0";
  collection_an.addTimestep(fileName,i,".vts");

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
    return result;
}

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
