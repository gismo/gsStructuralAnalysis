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
#include <gsStructuralAnalysis/gsALMRiks.h>
#include <gsStructuralAnalysis/gsHierarchicalALM.h>
#include <gsStructuralAnalysis/gsSpaceTimeHierarchy.h>
#include <gsStructuralAnalysis/gsAdaptiveSpaceTime.h>
#include <gsStructuralAnalysis/gsSpaceTimeFitter.h>
#include <gsStructuralAnalysis/gsSolutionFitter.h>

#include <gsHSplines/gsKdNode.h>


using namespace gismo;

template <class T>
void initStepOutput( const std::string name, const gsMatrix<T> & points);

template <class T>
void writeStepOutput(const T lambda, const gsVector<T> & solution, const gsMultiPatch<T> & deformation, const std::string name, const gsMatrix<T> & points, const index_t extreme=-1, const index_t kmax=100);

template <class T>
void writeStepOutput2(const T lambda, const gsMultiPatch<T> & deformation, const std::string name, const gsMatrix<T> & points, const index_t extreme=-1, const index_t kmax=100);


template<index_t dim, class T>
gsTensorBSpline<dim,T> gsSpaceTimeFit(const std::vector<gsMatrix<T>> & solutionCoefs, const gsVector<T> & times, const gsVector<T> & ptimes, gsMultiBasis<T> & spatialBasis, index_t deg = 2);

template<class T>
gsBSpline<T> gsSolutionFit(const std::vector<gsVector<T>> & solutions, const gsVector<T> & times, const gsVector<T> & ptimes, index_t deg = 2);

template<class T>
gsMatrix<T> nearestPoint(const gsMatrix<T>& spacePoint, const gsBSpline<T> &bspline, index_t nTrialPoints = 20, index_t nIterations = 50);



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

  if (write)
  {
    initStepOutput(dirname + "/" + wn, writePoints);
    initStepOutput(dirname + "/" + line, writePoints);
  }

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
  // arcLength = new gsALMRiks<real_t>(Jacobian, ALResidual, Force);

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


  gsVector<> pt(2);
  pt<<0.5,0.5;
  gsMatrix<> pt_result;
  gsMultiPatch<> mp_tmp;
  /////////////////////////////////////////////////////////////////////////////////////////////
  std::ofstream file, fitfile;
  file.open("init.txt",std::ofstream::out);
  assembler->constructSolution(U0,mp_tmp);
  deformation.patch(0) = mp_tmp.patch(0);
  deformation.patch(0).coefs() -= mp.patch(0).coefs();// assuming 1 patch here
  deformation.patch(0).eval_into(pt,pt_result);
  file  << std::setprecision(20) << pt_result(2,0) << "," << L0<<"\n";
  //////////////////////////////////////////////////////////////////////////////////////////////


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

    //////////////////////////////////////////////////////////////////////////////////////////////
    assembler->constructSolution(Uold,mp_tmp);
    deformation.patch(0) = mp_tmp.patch(0);
    deformation.patch(0).coefs() -= mp.patch(0).coefs();// assuming 1 patch here
    deformation.patch(0).eval_into(pt,pt_result);
    file  << std::setprecision(20) << pt_result(2,0) << "," << Lold<<"\n";
    //////////////////////////////////////////////////////////////////////////////////////////////


    real_t lambda = arcLength->solutionL();
    solutions.push_back({arcLength->solutionU(),lambda});
    times.push_back(s);

    if (!bisected)
    {
      dL = dL0;
      arcLength->setLength(dL);
    }
    bisected = false;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////
  file.close();
  //////////////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////

  if (make_fit)
  {
    // Store solution coefficients in a matrix
    std::vector<gsMatrix<real_t>> solutionCoefs(solutions.size());
    std::vector<gsMatrix<real_t>> solutionContainer(solutions.size());
    gsVector<> loads(solutions.size());

    for (size_t k=0; k!= solutions.size(); k++)
    {
      assembler->constructSolution(solutions[k].first,mp_tmp);
      solutionCoefs.at(k) = mp_tmp.patch(0).coefs();
      solutionContainer.at(k) = solutions[k].first;

      loads.at(k) = solutions[k].second;
    }

    gsSpaceTimeFitter<2,real_t> fitter( solutionCoefs,
                                        loads,
                                        gsAsVector<>(times),
                                        dbasis,
                                        deg_z);
    fitter.compute();

    gsSolutionFitter<real_t> cfitter(solutionContainer,
                                        loads,
                                        gsAsVector<>(times),
                                        deg_z);
    cfitter.compute();

    typename gsTensorBSpline<3,real_t>::BoundaryGeometryType target;

    gsParaviewCollection collection(dirname + "/" + output);
    gsParaviewCollection datacollection(dirname + "/" + "data");

    gsField<> solField;

    gsVector<> xi;
    xi.setLinSpaced(100,times[0],times[times.size()-1]);
    if (plot || write)
    {
      if (write)
        initStepOutput(dirname + "/" + line, writePoints);

      gsMatrix<> Uold;
      real_t Lold;

      for (index_t k = 0; k!=xi.size(); k++)
      {
        std::tie(Lold,Uold) = cfitter.slice(xi.at(k));

        assembler->constructSolution(Uold,mp_tmp);
        deformation.patch(0) = mp_tmp.patch(0);
        deformation.patch(0).coefs() -= mp.patch(0).coefs();// assuming 1 patch here
        deformation.patch(0).eval_into(pt,pt_result);
        if (plot)
        {
          solField = gsField<>(mp,deformation);
          gsWriteParaview(solField,"slice");

          gsDebugVar(xi[k]);

          std::string fileName = dirname + "/" + output + util::to_string(k);
          gsWriteParaview<>(solField, fileName, 1000,mesh);
          fileName = output + util::to_string(k) + "0";
          collection.addTimestep(fileName,xi[k],".vts");
          if (mesh) collection.addTimestep(fileName,xi[k],"_mesh.vtp");
        }
        if (write)
        {
            writeStepOutput2(Lold,deformation, dirname + "/" + line, writePoints,1, 201);
        }
      }
    }

    for (size_t k=0; k!= solutions.size(); k++)
    {
      assembler->constructSolution(solutions[k].first,mp_tmp);

      real_t lambda = solutions[k].second;

      real_t Time = times[k];

      deformation.patch(0) = mp_tmp.patch(0);
      deformation.patch(0).coefs() -= mp.patch(0).coefs();// assuming 1 patch here

      if (plot)
      {
        solField = gsField<>(mp,deformation);
        std::string fileName = dirname + "/" + "data" + util::to_string(k);
        gsWriteParaview<>(solField, fileName, 1000,mesh);
        fileName = "data" + util::to_string(k) + "0";
        datacollection.addTimestep(fileName,Time,".vts");
        if (mesh) datacollection.addTimestep(fileName,Time,"_mesh.vtp");
      }
      if (write)
      {
          writeStepOutput(lambda,solutions[k].first,deformation, dirname + "/" + wn, writePoints,1, 201);
      }
    }

    if (plot)
    {
      collection.save();
      datacollection.save();
    }
  }
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////

  gsAdaptiveSpaceTime<real_t,solution_t> hierarchy(times,solutions);
  hierarchy.options().setInt("MaxLevel",2);
  hierarchy.options().setInt("Split",false);
  hierarchy.options().setReal("Tolerance",2e-1);
  hierarchy.init();
  hierarchy.printQueue();

  solution_t start, guess, reference;
  index_t ID;
  real_t tstart = 0;
  real_t tend = 0;
  real_t dt, dt0;
  index_t it = 0;
  index_t itmax = 100;
  real_t TOL = 1e-2;
  gsVector<> DeltaU;
  real_t DeltaL;
  index_t Nintervals;
  std::vector<solution_t> stepSolutions;
  std::vector<real_t> stepTimes;
  std::vector<real_t> distances;
  bisected = false;
  real_t dt_rem = 0;

  while (!hierarchy.empty() && it < itmax)
  {
    //////////////////////////////////////////////////////////////////////////////////////////////
    file.open("it_" + std::to_string(it) + ".txt",std::ofstream::out);
    //////////////////////////////////////////////////////////////////////////////////////////////
    gsDebugVar(it);

    Nintervals = 3;
    stepSolutions.resize(Nintervals);
    stepTimes.resize(Nintervals);
    distances.resize(Nintervals+1);

    std::tie(ID,tstart,tend,dt0,start,guess) = hierarchy.pop();
    std::tie(Uold,Lold) = start;
    std::tie(Uguess,Lguess) = guess;

    gsMatrix<> Uori = Uold;
    real_t Lori = Lold;

    gsDebugVar(dt0);

    dt0 = dt0 / Nintervals;
    dt = dt0;

    arcLength->setLength(dt);
    arcLength->setSolution(Uold,Lold);
    // arcLength->resetStep();
    arcLength->setPrevious(Uguess,Lguess);
    gsDebug<<"Start - ||u|| = "<<Uold.norm()<<", L = "<<Lold<<"\n";
    gsDebug<<"Guess - ||u|| = "<<Uguess.norm()<<", L = "<<Lguess<<"\n";

    gsVector<> tmpU = Uold-Uguess;
    real_t tmpL = Lold-Lguess;


    gsDebugVar(arcLength->distance(tmpU,tmpL));

    //////////////////////////////////////////////////////////////////////////////////////////////
    assembler->constructSolution(Uold,mp_tmp);
    deformation.patch(0) = mp_tmp.patch(0);
    deformation.patch(0).coefs() -= mp.patch(0).coefs();// assuming 1 patch here
    deformation.patch(0).eval_into(pt,pt_result);
    file  << std::setprecision(20) << pt_result(2,0) << "," << Lold<<"\n";
    //////////////////////////////////////////////////////////////////////////////////////////////

    real_t s = 0;

    gsInfo<<"Starting with ID "<<ID<<" from (|U|,L) = ("<<Uold.norm()<<","<<Lold<<"), curve time = "<<tstart<<"\n";
    for (index_t k = 0; k!=Nintervals; k++)
    {
      gsDebug<<"Interval "<<k+1<<" of "<<Nintervals<<"\n";
      gsDebug<<"Start - ||u|| = "<<Uold.norm()<<", L = "<<Lold<<"\n";
      arcLength->step();
      if (!(arcLength->converged()))
      {
        gsInfo<<"Error: Loop terminated, arc length method did not converge.\n";
        dt = dt / 2.;
        dt_rem += dt; // add the remainder of the interval to dt_rem
        arcLength->setLength(dt);
        arcLength->setSolution(Uold,Lold);
        bisected = true;
        k -= 1;
        continue;
      }
      GISMO_ENSURE(arcLength->converged(),"Loop terminated, arc length method did not converge.\n");

      stepSolutions.at(k) = std::make_pair(arcLength->solutionU(),arcLength->solutionL());
      stepTimes.at(k) = tstart + dt;
      DeltaU = arcLength->solutionU() - Uold;
      DeltaL = arcLength->solutionL() - Lold;

      distances.at(k) = arcLength->distance(DeltaU,DeltaL);
      gsDebugVar(arcLength->distance(DeltaU,DeltaL));

      s += distances.at(k);

      Uold = arcLength->solutionU();
      Lold = arcLength->solutionL();

      //////////////////////////////////////////////////////////////////////////////////////////////
      assembler->constructSolution(Uold,mp_tmp);
      deformation.patch(0) = mp_tmp.patch(0);
      deformation.patch(0).coefs() -= mp.patch(0).coefs();// assuming 1 patch here
      deformation.patch(0).eval_into(pt,pt_result);
      file  << std::setprecision(20) << pt_result(2,0) << "," << Lold<<"\n";
      //////////////////////////////////////////////////////////////////////////////////////////////

      if (!bisected) // if bisected = false
        dt = dt0;
      else
      {
        // if the current interval has finished, but was refined before.
        // The next interval should have the remaining length.
        // Also, Nintervals should increase
        //
        dt = dt_rem;
        Nintervals++;
        stepSolutions.resize(Nintervals);
        stepTimes.resize(Nintervals);
        distances.resize(Nintervals+1);
      }

      arcLength->setLength(dt);
      dt_rem = 0;
      bisected = false;
    }

    gsDebugVar(s);

    bool success = hierarchy.getReferenceByID(ID,reference);
    GISMO_ASSERT(success,"Reference not found");

    //////////////////////////////////////////////////////////////////////////////////////////////



    //////////////////////////////////////////////////////////////////////////////////////////////


    DeltaU = reference.first - arcLength->solutionU();
    DeltaL = reference.second - arcLength->solutionL();
    distances.back() = arcLength->distance(DeltaU,DeltaL);
    // gsDebugVar(arcLength->distance(DeltaU,DeltaL));

    // DeltaU = arcLength->solutionU() - Uori;
    // DeltaL = arcLength->solutionL() - Lori;
    // gsDebugVar(arcLength->distance(DeltaU,DeltaL));

    s += distances.back();
    // gsDebugVar(s);


    //////////////////////////////////////////////////////////////////////////////////////////////
    assembler->constructSolution(reference.first,mp_tmp);
    deformation.patch(0) = mp_tmp.patch(0);
    deformation.patch(0).coefs() -= mp.patch(0).coefs();// assuming 1 patch here
    deformation.patch(0).eval_into(pt,pt_result);
    file  << std::setprecision(20) << pt_result(2,0) << "," << reference.second<<"\n";
    //////////////////////////////////////////////////////////////////////////////////////////////
    file.close();



    fitfile.close();

    hierarchy.submit(ID,distances,stepTimes,stepSolutions);
    hierarchy.addJobs(ID);
    hierarchy.finishJob(ID);

    if (make_fit)
    {
      for (index_t k=0; k!=stepTimes.size(); k++)
        cfitter.addDataPoint(stepSolutions[k].first,stepSolutions[k].second,stepTimes[k],-1);

      cfitter.compute();

      if (plot || write)
      {
        if (write)
          initStepOutput(dirname + "/" + line + "_" + str::to_string(it), writePoints);

        gsMatrix<> Uold;
        real_t Lold;

        for (index_t k = 0; k!=xi.size(); k++)
        {
          std::tie(Lold,Uold) = cfitter.slice(xi.at(k));

          assembler->constructSolution(Uold,mp_tmp);
          deformation.patch(0) = mp_tmp.patch(0);
          deformation.patch(0).coefs() -= mp.patch(0).coefs();// assuming 1 patch here
          deformation.patch(0).eval_into(pt,pt_result);
          if (plot)
          {
            solField = gsField<>(mp,deformation);
            gsWriteParaview(solField,"slice");

            gsDebugVar(xi[k]);

            std::string fileName = dirname + "/" + output + util::to_string(k);
            gsWriteParaview<>(solField, fileName, 1000,mesh);
            fileName = output + util::to_string(k) + "0";
            collection.addTimestep(fileName,xi[k],".vts");
            if (mesh) collection.addTimestep(fileName,xi[k],"_mesh.vtp");
          }
          if (write)
          {
              writeStepOutput2(Lold,deformation, dirname + "/" + line, writePoints,1, 201);
          }
        }
      }
    }

    it++;
  }

  hierarchy.printQueue();
  std::tie(times,solutions) = hierarchy.getFlatSolution();

  hierarchy.printKnots();

  gsDebugVar(gsAsVector(times));


  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////

  gsParaviewCollection collection2(dirname + "/" + output + "_refit");
  gsParaviewCollection datacollection2(dirname + "/" + "data" + "_refit");


  if (plot || write)
  {
    for (size_t k=0; k!= solutions.size(); k++)
    {
      assembler->constructSolution(solutions[k].first,mp_tmp);

      real_t lambda = solutions[k].second;

      real_t Time = times[k];

      deformation.patch(0) = mp_tmp.patch(0);
      deformation.patch(0).coefs() -= mp.patch(0).coefs();// assuming 1 patch here

      gsField<> solField;
      if (plot)
      {
        solField = gsField<>(mp,deformation);
        std::string fileName = dirname + "/" + "data" + util::to_string(k);
        gsWriteParaview<>(solField, fileName, 1000,mesh);
        fileName = "data" + util::to_string(k) + "0";
        datacollection2.addTimestep(fileName,0,Time,".vts");
        if (mesh) datacollection2.addTimestep(fileName,0,Time,"_mesh.vtp");
      }
      if (write)
      {
          writeStepOutput(lambda,solutions[k].first,deformation, dirname + "/" + wn, writePoints,1, 201);
      }
    }
  }

  if (plot)
  {
    collection2.save();
    datacollection2.save();
  }


  delete arcLength;
  return result;
}

template <class T>
void initStepOutput(const std::string name, const gsMatrix<T> & points)
{
  std::ofstream file;
  file.open(name,std::ofstream::out);
  file  << std::setprecision(20)
        << "Deformation norm" << ",";
        for (index_t k=0; k!=points.cols(); k++)
        {
          file<< "point "<<k<<" - x" << ","
              << "point "<<k<<" - y" << ","
              << "point "<<k<<" - z" << ",";
        }

  file  << "Lambda" << ","
        << "Indicator"
        << "\n";
  file.close();

  gsInfo<<"Step results will be written in file: "<<name<<"\n";
}

template <class T>
void writeStepOutput(const T lambda, const gsVector<T> & solution, const gsMultiPatch<T> & deformation, const std::string name, const gsMatrix<T> & points, const index_t extreme, const index_t kmax) // extreme: the column of point indices to compute the extreme over (default -1)
{
  gsMatrix<T> P(2,1), Q(2,1);
  gsMatrix<T> out(3,points.cols());
  gsMatrix<T> tmp;

  for (index_t p=0; p!=points.cols(); p++)
  {
    P<<points.col(p);
    deformation.patch(0).eval_into(P,tmp);
    out.col(p) = tmp;
  }

  std::ofstream file;
  file.open(name,std::ofstream::out | std::ofstream::app);
  if (extreme==-1)
  {
    file  << std::setprecision(6)
          << solution.norm() << ",";
          for (index_t p=0; p!=points.cols(); p++)
          {
            file<< out(0,p) << ","
                << out(1,p) << ","
                << out(2,p) << ",";
          }

    file  << lambda << ","
          << "NA" << ","
          << "\n";
  }
  else if (extreme==0 || extreme==1)
  {
    gsMatrix<T> out2(kmax,points.cols()); // evaluation points in the rows, output (per coordinate) in columns
    for (int p = 0; p != points.cols(); p ++)
    {
      Q.at(1-extreme) = points(1-extreme,p);
      for (int k = 0; k != kmax; k ++)
      {
        Q.at(extreme) = 1.0*k/(kmax-1);
        deformation.patch(0).eval_into(Q,tmp);
        out2(k,p) = tmp.at(2); // z coordinate
      }
    }

    file  << std::setprecision(6)
          << solution.norm() << ",";
          for (index_t p=0; p!=points.cols(); p++)
          {
            file<< out(0,p) << ","
                << out(1,p) << ","
                << std::max(abs(out2.col(p).maxCoeff()),abs(out2.col(p).minCoeff())) << ",";
          }

    file  << lambda << ","
          << "NA" << ","
          << "\n";
  }
  else
    GISMO_ERROR("Extremes setting unknown");

  file.close();
}


template <class T>
void writeStepOutput2(const T lambda, const gsMultiPatch<T> & deformation, const std::string name, const gsMatrix<T> & points, const index_t extreme, const index_t kmax) // extreme: the column of point indices to compute the extreme over (default -1)
{
  gsMatrix<T> P(2,1), Q(2,1);
  gsMatrix<T> out(3,points.cols());
  gsMatrix<T> tmp;

  for (index_t p=0; p!=points.cols(); p++)
  {
    P<<points.col(p);
    deformation.patch(0).eval_into(P,tmp);
    out.col(p) = tmp;
  }

  std::ofstream file;
  file.open(name,std::ofstream::out | std::ofstream::app);
  if (extreme==-1)
  {
    file  << std::setprecision(6)
          << "NA" << ",";
          for (index_t p=0; p!=points.cols(); p++)
          {
            file<< out(0,p) << ","
                << out(1,p) << ","
                << out(2,p) << ",";
          }

    file  << lambda << ","
          << "NA" << ","
          << "\n";
  }
  else if (extreme==0 || extreme==1)
  {
    gsMatrix<T> out2(kmax,points.cols()); // evaluation points in the rows, output (per coordinate) in columns
    for (int p = 0; p != points.cols(); p ++)
    {
      Q.at(1-extreme) = points(1-extreme,p);
      for (int k = 0; k != kmax; k ++)
      {
        Q.at(extreme) = 1.0*k/(kmax-1);
        deformation.patch(0).eval_into(Q,tmp);
        out2(k,p) = tmp.at(2); // z coordinate
      }
    }

    file  << std::setprecision(6)
          << "NA" << ",";
          for (index_t p=0; p!=points.cols(); p++)
          {
            file<< out(0,p) << ","
                << out(1,p) << ","
                << std::max(abs(out2.col(p).maxCoeff()),abs(out2.col(p).minCoeff())) << ",";
          }

    file  << lambda << ","
          << "NA" << ","
          << "\n";
  }
  else
    GISMO_ERROR("Extremes setting unknown");

  file.close();
}

template<index_t dim, class T>
gsTensorBSpline<dim,T> gsSpaceTimeFit(const std::vector<gsMatrix<T>> & solutionCoefs, const gsVector<T> & times, const gsVector<T> & ptimes, gsMultiBasis<T> & spatialBasis, index_t deg)
{
  GISMO_ASSERT(solutionCoefs.size()==(size_t) times.rows(),"Number of time and solution steps should match! "<<solutionCoefs.size()<<"!="<<times.rows());
  GISMO_ASSERT(solutionCoefs.at(0).cols() == dim,"Is the dimension correct?"<<solutionCoefs.at(0).cols() <<"!="<<dim);
  index_t nsteps = times.rows();
  index_t bsize = solutionCoefs.at(0).rows();

  // Prepare fitting basis
  gsKnotVector<> kv(ptimes.minCoeff(),ptimes.maxCoeff(),nsteps-(deg+1),deg+1);
  gsBSplineBasis<T> lbasis(kv);


  //////// TO DO:
  //////// - Make multi-patch compatible
  //////// - Include dimension d

  // for (index_t p = 0; p!=dbasis.nBases(); ++p)
  // {
    gsTensorBSplineBasis<dim,T> tbasis(
                                            static_cast<gsBSplineBasis<T> *>(&spatialBasis.basis(0).component(0))->knots(),
                                            static_cast<gsBSplineBasis<T> *>(&spatialBasis.basis(0).component(1))->knots(),
                                            kv
                                            );
  // }

  gsMatrix<> rhs(times.size(),(dim+1)*bsize);
  gsVector<> ones; ones.setOnes(bsize);

  for (index_t lam = 0; lam!=nsteps; ++lam)
  {
    rhs.block(lam,0,1,dim * bsize) = solutionCoefs.at(lam).reshape(1,dim * bsize);
    rhs.block(lam,dim*bsize,1,bsize) = times.at(lam) * ones.transpose();
  }

  // get the Greville Abcissae (anchors)
  gsMatrix<> anchors = lbasis.anchors();

  // Get the collocation matrix at the anchors
  gsSparseMatrix<> C = lbasis.collocationMatrix(anchors);

  gsSparseSolver<>::LU solver;
  solver.compute(C);

  gsMatrix<> sol, coefs((nsteps)*bsize,dim+1);
  sol = solver.solve(rhs);

  for (index_t lam = 0; lam!=nsteps; ++lam)
  {
    gsMatrix<> tmp = sol.block(lam,0,1,dim * bsize);
    coefs.block(lam * bsize,0,bsize,dim) = tmp.reshape(bsize,dim);
    coefs.block(lam * bsize,dim,bsize,1) = sol.block(lam,dim*bsize,1,bsize).transpose();
  }

  // gsTensorBSpline<3,T> tspline = tbasis.makeGeometry(give(coefs)).release();
  gsTensorBSpline<dim,T> tspline(tbasis,give(coefs));
  return tspline;
}

template<index_t dim, class T>
class gsSpaceTimeFitter
{
public:
  gsSpaceTimeFitter  ( const std::vector<gsMatrix<T>> & solutionCoefs,
                      const gsVector<T> & times,
                      const gsVector<T> & ptimes,
                      const gsMultiBasis<T> & spatialBasis,
                      const index_t deg = 2)
  :
  m_data(solutionCoefs),
  m_times(times),
  m_ptimes(ptimes),
  m_bases(spatialBasis),
  m_deg(deg)
  {

  }

  gsTensorBSpline<dim,T> compute()
  {
    GISMO_ASSERT(m_data.size()==m_times.rows(),"Number of time and solution steps should match! "<<m_data.size()<<"!="<<m_times.rows());
    GISMO_ASSERT(m_data.at(0).cols()== dim,"Is the dimension correct?"<<m_data.at(0).cols()<<"!="<<dim);
    index_t nsteps = m_times.rows();
    index_t bsize = m_data.at(0).rows();

    // Prepare fitting basis
    gsKnotVector<> kv(m_ptimes.minCoeff(),m_ptimes.maxCoeff(),nsteps-(m_deg+1),m_deg+1);
    gsBSplineBasis<T> lbasis(kv);


    //////// TO DO:
    //////// - Make multi-patch compatible
    //////// - Include dimension d

    // for (index_t p = 0; p!=dbasis.nBases(); ++p)
    // {
      gsTensorBSplineBasis<dim,T> tbasis(
                                              static_cast<gsBSplineBasis<T> *>(&m_bases.basis(0).component(0))->knots(),
                                              static_cast<gsBSplineBasis<T> *>(&m_bases.basis(0).component(1))->knots(),
                                              kv
                                              );
    // }

    gsMatrix<> rhs(m_times.size(),(dim+1)*bsize);
    gsVector<> ones; ones.setOnes(bsize);

    for (index_t lam = 0; lam!=nsteps; ++lam)
    {
      rhs.block(lam,0,1,dim * bsize) = m_data.at(lam).reshape(dim * bsize);
      rhs.block(lam,dim*bsize,1,bsize) = m_times.at(lam) * ones.transpose();
    }

    // get the Greville Abcissae (anchors)
    gsMatrix<> anchors = lbasis.anchors();

    // Get the collocation matrix at the anchors
    gsSparseMatrix<> C;
    lbasis.collocationMatrix(anchors,C);

    gsSparseSolver<>::LU solver;
    solver.compute(C);

    gsMatrix<> sol, coefs((nsteps)*bsize,dim+1);
    sol = solver.solve(rhs);

    for (index_t lam = 0; lam!=nsteps; ++lam)
    {
      gsMatrix<> tmp = sol.block(lam,0,1,dim * bsize);
      coefs.block(lam * bsize,0,bsize,dim) = tmp.reshape(bsize,dim);
      coefs.block(lam * bsize,dim,bsize,1) = sol.block(lam,dim*bsize,1,bsize).transpose();
    }

    // gsTensorBSpline<3,T> tspline = tbasis.makeGeometry(give(coefs)).release();
    gsTensorBSpline<dim,T> tspline(tbasis,give(coefs));
    return tspline;
  }

protected:
  gsMatrix<T> m_data;
  gsVector<T> m_times;
  gsVector<T> m_ptimes;
  gsMultiBasis<T> m_bases;
  index_t m_deg;

};
