/** @file gsThinShell_ArcLength.cpp

    @brief Code for the arc-length method of a shell based on loads

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#include <gismo.h>

#include <gsKLShell/gsThinShellAssembler.h>
#include <gsKLShell/getMaterialMatrix.h>
#include <gsStructuralAnalysis/gsALMBase.h>
#include <gsStructuralAnalysis/gsALMCrisfield.h>
#include <gsStructuralAnalysis/gsSolutionFitter.h>
#include <gsStructuralAnalysis/gsSpaceTimeFitter.h>
#include <gsStructuralAnalysis/gsStructuralAnalysisUtils.h>

#include <gsHSplines/gsKdNode.h>


using namespace gismo;

int main (int argc, char** argv)
{
  // Input options
  int numElevate    = 1;
  int numHref       = 1;
  bool plot         = false;
  bool make_fit     = false;
  bool mesh         = false;
  bool stress       = false;
  bool membrane     = false;
  bool quasiNewton  = false;
  int quasiNewtonInt= -1;
  bool adaptive     = false;
  int step          = 10;
  int method        = 2; // (0: Load control; 1: Riks' method; 2: Crisfield's method; 3: consistent crisfield method; 4: extended iterations)
  bool deformed     = false;

  bool composite = false;

  real_t relax      = 1.0;

  int testCase      = 0;

  int result        = 0;

  bool write        = false;

  index_t maxit     = 20;
  // index_t iniLevels  = 2;
  // index_t maxLevels  = 4;
  index_t maxLevel  = 2;

  // Arc length method options
  real_t dL        = 0.5; // Ard length to find bifurcation
  real_t tol        = 1e-6;
  real_t tolU       = 1e-6;
  real_t tolF       = 1e-3;

  index_t deg_z = 2;

  std::string wn("data.csv");

  std::string assemberOptionsFile("options/solver_options.xml");

  gsCmdLine cmd("Arc-length analysis for thin shells.");
  cmd.addString( "f", "file", "Input XML file for assembler options", assemberOptionsFile );

  cmd.addInt("t", "testcase", "Test case: 0: clamped-clamped, 1: pinned-pinned, 2: clamped-free", testCase);

  cmd.addInt("r","hRefine", "Number of dyadic h-refinement (bisection) steps to perform before solving", numHref);
  cmd.addInt("e","degreeElevation", "Number of degree elevation steps to perform on the Geometry's basis before solving", numElevate);
  cmd.addSwitch("composite", "Composite material", composite);

  cmd.addInt("m","Method", "Arc length method; 1: Crisfield's method; 2: RIks' method.", method);
  cmd.addReal("L","dL", "arc length", dL);
  // cmd.addInt("I","inilvl", "Initial levels", iniLevels);
  // cmd.addInt("M","maxlvl", "Max levels", maxLevels);
  cmd.addInt("l","level", "Max level", maxLevel);
  cmd.addReal("A","relaxation", "Relaxation factor for arc length method", relax);

  cmd.addInt("q","QuasiNewtonInt","Use the Quasi Newton method every INT iterations",quasiNewtonInt);
  cmd.addInt("N", "maxsteps", "Maximum number of steps", step);
  cmd.addInt("z", "degz", "Degree of fitting splin", deg_z);

  cmd.addSwitch("adaptive", "Adaptive length ", adaptive);
  cmd.addSwitch("quasi", "Use the Quasi Newton method", quasiNewton);
  cmd.addSwitch("plot", "Plot result in ParaView format", plot);
  cmd.addSwitch("mesh", "Plot mesh?", mesh);
  cmd.addSwitch("stress", "Plot stress in ParaView format", stress);
  cmd.addSwitch("membrane", "Use membrane model (no bending)", membrane);
  cmd.addSwitch("deformed", "plot on deformed shape", deformed);
  cmd.addSwitch("write", "write to file", write);
  cmd.addSwitch("fit", "Make a fit and refer to that", make_fit);

  try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

  // GISMO_ASSERT(maxLevels>iniLevels,"Max levels must  be more than initial levels!");

  gsFileData<> fd(assemberOptionsFile);
  gsOptionList opts;
  fd.getFirst<gsOptionList>(opts);

  gsMultiPatch<> mp;

  real_t thickness;
  real_t Exx, Eyy, Gxy;
  real_t PoissonRatio = 0.3;
  real_t Density    = 1e0;

  if (composite)
  {
    Exx = 3300;
    Eyy = 1100;
    Gxy = 660;
    PoissonRatio = 0.25;
  }
  else
  {
    Exx  = 3102.75;
    PoissonRatio = 0.3;
  }


  if (testCase==1)
  {
    thickness = 6.35;
  }
  else if (testCase==2)
  {
    thickness = 12.7;
  }
  else if (testCase==3)
  {
    thickness = 16.75;
  }

  gsReadFile<>("surface/scordelis_lo_roof_shallow.xml", mp);

  for(index_t i = 0; i< numElevate; ++i)
    mp.patch(0).degreeElevate();    // Elevate the degree

  // h-refine
  for(index_t i = 0; i< numHref; ++i)
    mp.patch(0).uniformRefine();

  gsMultiBasis<> dbasis(mp);
  gsInfo<<"Basis (patch 0): "<< mp.patch(0).basis() << "\n";

  // Boundary conditions
  gsBoundaryConditions<> BCs;
  BCs.setGeoMap(mp);
  gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();

  // Initiate Surface forces
  std::string tx("0");
  std::string ty("0");
  std::string tz("0");

  gsVector<> tmp(mp.targetDim());
  gsVector<> neu(mp.targetDim());
  tmp.setZero();
  neu.setZero();
  gsConstantFunction<> neuData(neu,mp.targetDim());

  // Unscaled load
  real_t Load = 0;

  std::string output = "solution";
  std::string dirname = "ArcLengthResults";

  gsMatrix<> writePoints(2,3);
  writePoints.col(0)<< 0.0,0.5;
  writePoints.col(1)<< 0.5,0.5;
  writePoints.col(2)<< 1.0,0.5;

  GISMO_ASSERT(mp.targetDim()==3,"Geometry must be surface (targetDim=3)!");
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

  dirname = dirname + "/" +  "Roof_t="+ std::to_string(thickness) + "-r=" + std::to_string(numHref) + "-e" + std::to_string(numElevate) +"_solution";
  output =  "solution";
  wn = "data.txt";
  std::string line = "line.txt";

  std::string commands = "mkdir -p " + dirname;
  const char *command = commands.c_str();
  system(command);

  // plot geometry
  if (plot)
    gsWriteParaview(mp,dirname + "/" + "mp",1000,mesh);

  // Initialise solution object
  gsMultiPatch<> mp_def = mp;

  gsMaterialMatrixBase<real_t>* materialMatrix;

  real_t pi = math::atan(1)*4;
  index_t kmax = 3;

  std::vector<gsFunctionSet<> * > Gs(kmax);
  std::vector<gsFunctionSet<> * > Ts(kmax);
  std::vector<gsFunctionSet<> * > Phis(kmax);

  gsMatrix<> Gmat = gsCompositeMatrix(Exx,Eyy,Gxy,PoissonRatio,PoissonRatio*Eyy/Exx);
  Gmat.resize(Gmat.rows()*Gmat.cols(),1);
  gsConstantFunction<> Gfun(Gmat,3);
  Gs[0] = Gs[1] = Gs[2] = &Gfun;

  gsConstantFunction<> phi1,phi2,phi3;
  phi1.setValue(pi/2,3);
  phi2.setValue(0,3);
  phi3.setValue(pi/2,3);

  Phis[0] = &phi1;
  Phis[1] = &phi2;
  Phis[2] = &phi3;

  gsConstantFunction<> thicks(thickness/kmax,3);
  Ts[0] = Ts[1] = Ts[2] = &thicks;

  gsConstantFunction<> force(tmp,3);
  gsFunctionExpr<> t(std::to_string(thickness), 3);
  gsFunctionExpr<> E(std::to_string(Exx),3);
  gsFunctionExpr<> nu(std::to_string(PoissonRatio),3);
  gsFunctionExpr<> rho(std::to_string(Density),3);

  std::vector<gsFunction<>*> parameters;
  parameters.resize(2);
  parameters[0] = &E;
  parameters[1] = &nu;

  gsOptionList options;
  if (composite)
  {
      materialMatrix = new gsMaterialMatrixComposite<3,real_t>(mp,Ts,Gs,Phis);
  }
  else
  {
      options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
      options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
      materialMatrix = getMaterialMatrix<3,real_t>(mp,t,parameters,rho,options);
  }

  gsThinShellAssemblerBase<real_t>* assembler;
  assembler = new gsThinShellAssembler<3, real_t, true >(mp,dbasis,BCs,force,materialMatrix);

  // Construct assembler object
  assembler->setOptions(opts);
  assembler->setPointLoads(pLoads);

  typedef std::function<gsSparseMatrix<real_t> (gsVector<real_t> const &)>                                Jacobian_t;
  typedef std::function<gsVector<real_t> (gsVector<real_t> const &, real_t, gsVector<real_t> const &) >   ALResidual_t;
  // Function for the Jacobian
  Jacobian_t Jacobian = [&assembler](gsVector<real_t> const &x)
  {
    gsMultiPatch<> mp_def;
    assembler->constructSolution(x,mp_def);
    assembler->assembleMatrix(mp_def);
    gsSparseMatrix<real_t> m = assembler->matrix();
    return m;
  };
  // Function for the Residual
  ALResidual_t ALResidual = [&assembler](gsVector<real_t> const &x, real_t lam, gsVector<real_t> const &force)
  {
    gsMultiPatch<> mp_def;
    assembler->constructSolution(x,mp_def);
    assembler->assembleVector(mp_def);
    gsVector<real_t> Fint = -(assembler->rhs() - force);
    gsVector<real_t> result = Fint - lam * force;
    return result; // - lam * force;
  };
  // Assemble linear system to obtain the force vector
  assembler->assemble();
  gsVector<> Force = assembler->rhs();

  gsALMBase<real_t> * arcLength;
  arcLength = new gsALMCrisfield<real_t>(Jacobian, ALResidual, Force);

  arcLength->options().setString("Solver","SimplicialLDLT"); // LDLT solver
  arcLength->options().setInt("BifurcationMethod",0); // 0: determinant, 1: eigenvalue
  arcLength->options().setReal("Length",dL);
  // arcLength->options().setInt("AngleMethod",0); // 0: step, 1: iteration
  arcLength->options().setSwitch("AdaptiveLength",adaptive);
  arcLength->options().setInt("AdaptiveIterations",5);
  arcLength->options().setReal("Scaling",1.0);
  arcLength->options().setReal("Tol",tol);
  arcLength->options().setReal("TolU",tolU);
  arcLength->options().setReal("TolF",tolF);
  arcLength->options().setInt("MaxIter",maxit);
  arcLength->options().setSwitch("Verbose",true);
  arcLength->options().setReal("Relaxation",relax);
  if (quasiNewtonInt>0)
  {
    quasiNewton = true;
    arcLength->options().setInt("QuasiIterations",quasiNewtonInt);
  }
  arcLength->options().setSwitch("Quasi",quasiNewton);


  gsInfo<<arcLength->options();
  arcLength->applyOptions();
  arcLength->initialize();

  gsMultiPatch<> deformation = mp;

  // Make objects for previous solutions
  real_t Lguess,Lold, L0;
  gsMatrix<> Uguess,Uold, U0;
  Uold.setZero(Force.size(),1);
  U0.setZero(Force.size(),1);
  L0 = Lold = 0.0;

  gsMatrix<> solVector;
  real_t indicator = 0.0;
  arcLength->setIndicator(indicator); // RESET INDICATOR
  real_t dL0 = dL;


  typedef std::pair<gsVector<real_t>,real_t> solution_t;

  /*
    \a solutions is a container the for each level contains the solutions per point
    \a points contains all the points across levels in the format (level, U, lambda) -------------------------------> OVERKILL? WHY NEEDED?
    \a refPoints is a container that contains (level, U, lambda) of the points from which a refinement should START in level+1
    \a errors is a container that contains the error[l][i] e_i at the ith point of level l
  */

  std::vector<solution_t> solutions;
  std::vector<real_t> times;
  real_t s = 0;

  index_t level = 0;
  gsInfo<<"------------------------------------------------------------------------------------\n";
  gsInfo<<"\t\t\tLevel "<<level<<" (dL = "<<dL<<") -- Coarse grid \n";
  gsInfo<<"------------------------------------------------------------------------------------\n";

  bool bisected = false;

  index_t stepi = step; // number of steps for level i
  stepi = step * (math::pow(2,level));

  // Add the undeformed solution
  solutions.push_back({U0,L0});
  times.push_back(s);
  // Add other solutions
  arcLength->setSolution(U0,L0);

  std::vector<gsVector<real_t>> solutionContainer(solutions.size());
  index_t N = solutions[0].first.rows();
  for (size_t k=0; k!= solutions.size(); k++)
  {
    solutionContainer.at(k).resize(N+1);
    solutionContainer.at(k).head(N) = solutions[k].first;
    solutionContainer.at(k).at(N) = solutions[k].second;
  }

  gsSolutionFitter<real_t> cfitter(solutionContainer,
                                      gsAsVector<>(times),
                                      deg_z);

  for (index_t k=1; k<stepi+1; k++)
  {
    s+=dL;

    gsInfo<<"Load step "<< k<<"\t"<<"dL = "<<dL<<"; curve time = "<<s<<"\n";
    // arcLength->setLength(dL);
    // assembler->constructSolution(solVector,solution);
    arcLength->step();

    if (!(arcLength->converged()))
    {
      s -= dL;
      gsInfo<<"Error: Loop terminated, arc length method did not converge.\n";
      dL = dL / 2.;
      arcLength->setLength(dL);
      arcLength->setSolution(Uold,Lold);
      bisected = true;
      k -= 1;
      continue;
    }

    Uold = arcLength->solutionU();
    Lold = arcLength->solutionL();

    real_t lambda = arcLength->solutionL();
    solutions.push_back({arcLength->solutionU(),lambda});
    times.push_back(s);

    // add data point to fit
    index_t N = arcLength->solutionU().rows();
    gsVector<real_t> sol(N+1);
    sol.head(N) = arcLength->solutionU();
    sol.at(N) = lambda;
    cfitter.addDataPoint(sol,s);
    if (k>deg_z+1)
      cfitter.compute();

    // if (k>deg_z+1)
    // {
    //   cfitter.compute();

    //   gsVector<> xi;
    //   xi.setLinSpaced(100,times[0],times[times.size()-1]);

    //   if (write)
    //   {
    //     gsALMOutput<real_t> line(dirname + "/line.csv",refPoints);
    //     if (write)
    //         line.init(headers);

    //     gsField<> solField;
    //     gsMultiPatch<> mp_tmp;
    //     gsMatrix<> tmp;

    //     for (index_t k = 0; k!=xi.size(); k++)
    //     {
    //       tmp = cfitter.slice(xi.at(k));
    //       Lold = tmp(0,N);
    //       Uold = tmp.row(0).head(N);

    //       gsMatrix<> pointResults(mp.geoDim(),refPoints.cols());
    //       assembler->constructSolution(Uold,mp_tmp);
    //       deformation.patch(0) = mp_tmp.patch(0);
    //       deformation.patch(0).coefs() -= mp.patch(0).coefs();// assuming 1 patch here
    //       solField = gsField<>(mp,deformation);
    //       for (index_t p=0; p!=refPoints.cols(); p++)
    //           pointResults.col(p) = solField.value(refPoints.col(p),refPatches(0,p));
    //       line.add(pointResults,Lold);
    //     }

    //     gsALMOutput<real_t> coefs(dirname + "/coefs.csv",refPoints);
    //     if (write)
    //         coefs.init(headers);

    //     tmp = cfitter.fit().coefs();
    //     gsMatrix<> res(2,tmp.rows());
    //     for (index_t x = 0; x!=tmp.rows(); x++)
    //     {
    //       Lold = tmp(x,N);
    //       Uold = tmp.row(x).head(N);
    //       gsMatrix<> pointResults(mp.geoDim(),refPoints.cols());
    //       assembler->constructSolution(Uold,mp_tmp);
    //       deformation.patch(0) = mp_tmp.patch(0);
    //       deformation.patch(0).coefs() -= mp.patch(0).coefs();// assuming 1 patch here
    //       solField = gsField<>(mp,deformation);
    //       for (index_t p=0; p!=refPoints.cols(); p++)
    //           pointResults.col(p) = solField.value(refPoints.col(p),refPatches(0,p));
    //       coefs.add(pointResults,Lold);
    //     }
    //   }
    // }

    //////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////

    if (!bisected)
    {
      dL = dL0;
      arcLength->setLength(dL);
    }
    bisected = false;
  }

  if (plot || write)
  {
    gsVector<> refPoints(2);
    gsVector<index_t> refPatches(1);
    refPoints<<0.5,0.5;
    refPatches<<0;
    gsField<> solField;
    gsMultiPatch<> mp_tmp;
    std::vector<std::string> headers = {"u_x","u_y","u_z"};

    gsParaviewCollection dataCollection(dirname + "/" + "data");
    gsALMOutput<real_t> data(dirname + "/data.csv",refPoints);
    if (write)
        data.init(headers);

    for (size_t k=0; k!= solutions.size(); k++)
    {
      Lold = solutions[k].second;
      Uold = solutions[k].first;

      assembler->constructSolution(Uold,mp_tmp);
      deformation.patch(0) = mp_tmp.patch(0);
      deformation.patch(0).coefs() -= mp.patch(0).coefs();// assuming 1 patch here
      solField = gsField<>(mp,deformation);

      if (write)
        if (refPoints.cols()!=0)
        {
            gsMatrix<> pointResults(mp.geoDim(),refPoints.cols());
            for (index_t p=0; p!=refPoints.cols(); p++)
                pointResults.col(p) = solField.value(refPoints.col(p),refPatches(0,p));
            data.add(pointResults,Lold);
        }

      if (plot)
      {
        std::string fileName = dirname + "/" + "data" + util::to_string(k);
        gsWriteParaview<>(solField, fileName, 1000,mesh);
        fileName = "data" + util::to_string(k) + "0";
        dataCollection.addTimestep(fileName,times[k],".vts");
        if (mesh) dataCollection.addTimestep(fileName,times[k],"_mesh.vtp");
      }
    }

    gsParaviewCollection lineCollection(dirname + "/" + "line");
    gsALMOutput<real_t> line(dirname + "/line.csv",refPoints);
    if (write)
        line.init(headers);
    gsVector<> xi;
    xi.setLinSpaced(times.size()*10,times[0],times[times.size()-1]);
    gsVector<> sol;
    for (index_t k = 0; k!=xi.size(); k++)
    {
      sol = cfitter.slice(xi.at(k));
      Lold = sol.at(N);
      Uold = sol.head(N);

      assembler->constructSolution(Uold,mp_tmp);
      deformation.patch(0) = mp_tmp.patch(0);
      deformation.patch(0).coefs() -= mp.patch(0).coefs();// assuming 1 patch here
      solField = gsField<>(mp,deformation);

      if (write)
      {
        gsMatrix<> pointResults(mp.geoDim(),refPoints.cols());
        for (index_t p=0; p!=refPoints.cols(); p++)
            pointResults.col(p) = solField.value(refPoints.col(p),refPatches(0,p));
        line.add(pointResults,Lold);
      }

      if (plot)
      {
        std::string fileName = dirname + "/" + "line" + util::to_string(k);
        gsWriteParaview<>(solField, fileName, 1000,mesh);
        fileName = "line" + util::to_string(k) + "0";
        lineCollection.addTimestep(fileName,xi[k],".vts");
        if (mesh) lineCollection.addTimestep(fileName,xi[k],"_mesh.vtp");
      }
    }

    gsALMOutput<real_t> coefs(dirname + "/coefs.csv",refPoints);
    if (write)
        coefs.init(headers);

    gsMatrix<> solCoefs = cfitter.fit().coefs();
    gsMatrix<> res(2,sol.rows());
    for (index_t x = 0; x!=solCoefs.rows(); x++)
    {
      Lold = solCoefs(x,N);
      Uold = solCoefs.row(x).head(N);
      gsMatrix<> pointResults(mp.geoDim(),refPoints.cols());
      assembler->constructSolution(Uold,mp_tmp);
      deformation.patch(0) = mp_tmp.patch(0);
      deformation.patch(0).coefs() -= mp.patch(0).coefs();// assuming 1 patch here
      solField = gsField<>(mp,deformation);
      for (index_t p=0; p!=refPoints.cols(); p++)
        pointResults.col(p) = solField.value(refPoints.col(p),refPatches(0,p));
      coefs.add(pointResults,Lold);
    }

    if (plot)
    {
      lineCollection.save();
      dataCollection.save();
    }
  }

  // gsMatrix<> spacePoint(solutions[1].first.size()+1,1);
  // spacePoint.block(0,0,solutions[1].first.size(),1) = solutions[1].first;
  // spacePoint(solutions[1].first.size(),0) = solutions[1].second;
  // gsMatrix<> bla;

  // cfitter.nearestParam(spacePoint,bla);
  // gsDebugVar(bla);


  return 0;

  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////

  // Store solution coefficients in a matrix
  // index_t blocksize = mp.patch(0).coefs().rows();
  // gsMatrix<> solutionCoefs(solutions.size(),3*blocksize);
  std::vector<gsMatrix<real_t>> solutionCoefs(solutions.size());
  gsVector<> loads(solutions.size());

  for (size_t k=0; k!= solutions.size(); k++)
  {
    gsMultiPatch<> mp_tmp;
    assembler->constructSolution(solutions[k].first,mp_tmp);
    solutionCoefs.at(k) = mp_tmp.patch(0).coefs();

    loads.at(k) = solutions[k].second;
  }

  gsSpaceTimeFitter<2,real_t> fitter( solutionCoefs,
                                      loads,
                                      gsAsVector<>(times),
                                      dbasis,
                                      deg_z);

  fitter.compute();

  typename gsTensorBSpline<3,real_t>::BoundaryGeometryType target;



  delete arcLength;
  return result;
}