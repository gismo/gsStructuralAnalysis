/** @file gsThinShell_ArcLength.cpp

    @brief Code for the arc-length method of a shell based on loads

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

#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMBase.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMCrisfield.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMRiks.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsAPALMData.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsAPALM.h>
#include <gsStructuralAnalysis/src/gsStructuralAnalysisTools/gsStructuralAnalysisUtils.h>
using namespace gismo;

template <class T>
gsMultiPatch<T> RectangularDomain(int n, int m, int p, int q, T L, T B, bool clamped = false, T offset = 0.1);
template <class T>
gsMultiPatch<T> RectangularDomain(int n, int p, T L, T B, bool clamped = false, T offset = 0.1);

template <class T>
gsMultiPatch<T> Rectangle(T L, T B);

#ifdef gsKLShell_ENABLED
template<class T>
class gsAPALMBeam : public gsAPALM<T>
{

  using Base = gsAPALM<T>;
  typedef typename Base::solution_t solution_t;

public:
  gsAPALMBeam(const gsMpiComm & comm,
              gsALMBase<T> * ALM,
              const gsAPALMData<T,solution_t> & Data,
              const gsThinShellAssemblerBase<T> * assembler,
              std::string & dirname,
              const gsMatrix<T> & refPoints,
              const gsVector<index_t> & refPatches          )
  :
  Base(ALM,Data,comm),
  m_assembler(assembler),
  m_dirname(dirname),
  m_refPoints(refPoints),
  m_refPatches(refPatches)
  {

  }

  void parallelIntervalOutput(const std::vector<std::pair<gsVector<T>,T>> & stepSolutions, const std::vector<T> & stepTimes, index_t level, index_t ID)
  {
    gsStructuralAnalysisOutput<real_t> data(m_dirname + "/interval_"+std::to_string(ID)+".csv",m_refPoints);
    gsMultiPatch<T> deformation,mp_tmp, mp;
    deformation = mp = m_assembler->geometry();
    std::vector<std::string> pointheaders = {"u_x","u_y","u_z"};
    std::vector<std::string> otherheaders = {"U-norm","lambda","time","level"};

    data.init(pointheaders,otherheaders);

    for (size_t k=0; k!=stepSolutions.size(); k++)
    {
      m_assembler->constructSolution(stepSolutions[k].first,mp_tmp);
      deformation.patch(0) = mp_tmp.patch(0);
      deformation.patch(0).coefs() -= mp.patch(0).coefs();// assuming 1 patch here
      gsField<T> solField(mp,deformation);

      if (m_refPoints.cols()!=0)
      {
          gsVector<> otherData(4);
          otherData<<stepSolutions[k].first.norm(),stepSolutions[k].second,stepTimes[k],level;
          gsMatrix<> pointResults(mp.geoDim(),m_refPoints.cols());
          for (index_t p=0; p!=m_refPoints.cols(); p++)
              pointResults.col(p) = solField.value(m_refPoints.col(p),m_refPatches.at(p));
          data.add(pointResults,otherData);
      }
    }
  }

protected:
  const gsThinShellAssemblerBase<T> * m_assembler;

  const std::string m_dirname;

  const gsMatrix<T> m_refPoints;
  const gsVector<index_t> m_refPatches;
};

int main (int argc, char** argv)
{

  // Input options
  int numElevate    = 2;
  int numHref       = 5;
  bool plot         = false;
  bool mesh         = false;
  bool stress       = false;
  bool membrane     = false;
  bool quasiNewton  = false;
  int quasiNewtonInt= -1;
  bool adaptive     = false;
  int step          = 10;
  int SubIntervals  = 2;
  int method        = 2; // (0: Load control; 1: Riks' method; 2: Crisfield's method; 3: consistent crisfield method; 4: extended iterations)
  bool deformed     = false;
  bool sequential   = false;

  bool composite = false;

  real_t relax      = 1.0;

  int testCase      = 1;

  int result        = 0;

  bool write        = false;

  index_t maxit     = 20;
  // index_t iniLevels  = 2;
  // index_t maxLevels  = 4;
  index_t maxLevel  = 2;

  // Arc length method options
  real_t dL         = 0; // General arc length
  real_t dLb        = 0.5; // Arc length to find bifurcation
  real_t ALM_tol    = 1e-6;
  real_t ALM_tolU   = 1e-6;
  real_t ALM_tolF   = 1e-3;

  real_t APALM_tol  = 1e-3;

  index_t verbose = 0;

  std::string wn("data.csv");

  std::string assemberOptionsFile("options/solver_options.xml");

  gsCmdLine cmd("Arc-length analysis for thin shells.");
  cmd.addString( "f", "file", "Input XML file for assembler options", assemberOptionsFile );

  cmd.addInt("t", "testcase", "Test case: 0: clamped-free with vertical load, 1: clamped-free with horizontal load", testCase);

  cmd.addInt("r","hRefine", "Number of dyadic h-refinement (bisection) steps to perform before solving", numHref);
  cmd.addInt("e","degreeElevation", "Number of degree elevation steps to perform on the Geometry's basis before solving", numElevate);
  cmd.addSwitch("composite", "Composite material", composite);

  cmd.addInt("m","Method", "Arc length method; 1: Crisfield's method; 2: RIks' method.", method);
  cmd.addReal("L","dLb", "arc length", dLb);
  cmd.addReal("l","dL", "arc length after bifurcation", dL);
  cmd.addReal("T", "APALM_tol", "APALM Tolerance", APALM_tol);
  cmd.addInt("d","level", "Max level", maxLevel);
  cmd.addReal("A","relaxation", "Relaxation factor for arc length method", relax);

  cmd.addInt("q","QuasiNewtonInt","Use the Quasi Newton method every INT iterations",quasiNewtonInt);
  cmd.addInt("N", "maxsteps", "Maximum number of steps", step);
  cmd.addInt("n", "SubIntervals", "Number of steps in subintervals", SubIntervals);

  cmd.addInt("v", "verbose", "verbose", verbose);

  cmd.addSwitch("adaptive", "Adaptive length ", adaptive);
  cmd.addSwitch("quasi", "Use the Quasi Newton method", quasiNewton);
  cmd.addSwitch("plot", "Plot result in ParaView format", plot);
  cmd.addSwitch("mesh", "Plot mesh?", mesh);
  cmd.addSwitch("stress", "Plot stress in ParaView format", stress);
  cmd.addSwitch("membrane", "Use membrane model (no bending)", membrane);
  cmd.addSwitch("deformed", "plot on deformed shape", deformed);
  cmd.addSwitch("write", "write to file", write);
  cmd.addSwitch("sequential", "Solve sequential (serial -> parallel) instead of mixed", sequential);

  try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

  if (dL==0)
  {
    dL = dLb;
  }

  gsFileData<> fd(assemberOptionsFile);
  gsOptionList opts;
  fd.getFirst<gsOptionList>(opts);

  gsMultiPatch<> mp;
  real_t aDim = 0;
  real_t bDim = 0;

  real_t E_modulus, thickness, PoissonRatio;
  /*
    Case 2: Clamped beam (left) under vertical end load                   --- Validation settings: -L 5e-1 -M 0 -N 10 -r 2 -e 1 (--plot --write -q 5)
                                                                              Fig 3b from: Pagani, A., & Carrera, E. (2018). Unified formulation of geometrically nonlinear refined beam theories. Mechanics of Advanced Materials and Structures, 25(1), 15–31. https://doi.org/10.1080/15376494.2016.1232458
    Case 3: Clamped beam (left) under horizontal compressive end load     --- Validation settings: -L 5e-5 -l 1e-1 -M 0 -N 100 -r 4 -e 2
                                                                              Fig 5  from: Pagani, A., & Carrera, E. (2018). Unified formulation of geometrically nonlinear refined beam theories. Mechanics of Advanced Materials and Structures, 25(1), 15–31. https://doi.org/10.1080/15376494.2016.1232458
  */
  E_modulus = 75e6;
  thickness = 0.01;
  PoissonRatio = 0.0;
  aDim = 1.0;
  bDim = 0.01;
  mp = RectangularDomain(numHref, 0, numElevate+2, 2, aDim, bDim);

  gsMultiBasis<> dbasis(mp);
  gsInfo<<"Basis (patch 0): "<< mp.patch(0).basis() << "\n";

  gsDebugVar(mp.patch(0).coefs());

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
  std::string cores;
  real_t pressure = 0.0;

  gsMatrix<> writePoints(2,3);
  writePoints.col(0)<< 0.0,0.5;
  writePoints.col(1)<< 0.5,0.5;
  writePoints.col(2)<< 1.0,0.5;

  const gsMpi & mpi = gsMpi::init();
  gsMpiComm comm = mpi.worldComm();
  cores = "_ncores="+std::to_string(comm.size());


  if (testCase == 0)
  {
      GISMO_ASSERT(mp.targetDim()==3,"Geometry must be surface (targetDim=3)!");
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

    dirname = dirname + "/Beam_clamped-verticalLoad" + cores;
    output =  "solution";
    wn = output + "data.txt";
  }
  else if (testCase == 1)
  {
    GISMO_ASSERT(mp.targetDim()==3,"Geometry must be surface (targetDim=3)!");
    BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
    BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
    BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

    BCs.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 2 );

    BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
    BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y

    // dL =  3e-0;
    // dLb = 0.8e-4;

    Load = 1e-1;
    gsVector<> point(2);
    gsVector<> load (3);
    point<< 1.0, 0.5 ;
    load << -Load, 0.0, 0.0 ;
    pLoads.addLoad(point, load, 0 );

    dirname = dirname + "/Beam_clamped-horizontalLoad"  + cores;
    output =  "solution";
    wn = output + "data.txt";
  }

  // Prepare and create directory with dirname
  dirname = gsFileManager::getCurrentPath() + dirname;
  if (sequential)
    dirname = dirname + "_seq";

  GISMO_ENSURE(gsFileManager::mkdir(dirname),"Failed to create directory " + dirname);
  // Made directory

  gsConstantFunction<> pressFun(pressure,3);
  // Initialise solution object
  gsMultiPatch<> mp_def = mp;

  // Linear isotropic material model
  gsConstantFunction<> force(tmp,3);
  gsFunctionExpr<> t(std::to_string(thickness), 3);
  gsFunctionExpr<> E(std::to_string(E_modulus),3);
  gsFunctionExpr<> nu(std::to_string(PoissonRatio),3);

  std::vector<gsFunctionSet<>*> parameters(2);
  parameters[0] = &E;
  parameters[1] = &nu;

  gsMaterialMatrixBase<real_t>* materialMatrix;
  materialMatrix = new gsMaterialMatrixLinear<3,real_t>(mp,t,parameters);

  gsThinShellAssemblerBase<real_t>* assembler;
  assembler = new gsThinShellAssembler<3, real_t, true >(mp,dbasis,BCs,force,materialMatrix);

  // Construct assembler object
  assembler->setOptions(opts);
  assembler->setPointLoads(pLoads);

  // Assemble linear system to obtain the force vector
  assembler->assemble();
  gsVector<> Force = assembler->rhs();

  // Function for the Jacobian
  gsStructuralAnalysisOps<real_t>::Jacobian_t Jacobian = [&assembler,&mp_def](gsVector<real_t> const &x, gsSparseMatrix<real_t> & m)
  {
    ThinShellAssemblerStatus status;
    assembler->constructSolution(x,mp_def);
    status = assembler->assembleMatrix(mp_def);
    m = assembler->matrix();
    return status == ThinShellAssemblerStatus::Success;
  };
  // Function for the Residual
  gsStructuralAnalysisOps<real_t>::ALResidual_t ALResidual = [&assembler,&mp_def,&Force](gsVector<real_t> const &x, real_t lam, gsVector<real_t> & result)
  {
    ThinShellAssemblerStatus status;
    assembler->constructSolution(x,mp_def);
    status = assembler->assembleVector(mp_def);
    result = Force - lam * Force - assembler->rhs(); // assembler rhs - force = Finternal
    return status == ThinShellAssemblerStatus::Success;
  };

  gsALMBase<real_t> * arcLength;
  arcLength = new gsALMCrisfield<real_t>(Jacobian, ALResidual, Force);
  // arcLength = new gsALMRiks<real_t>(Jacobian, ALResidual, Force);

  arcLength->options().setString("Solver","SimplicialLDLT"); // LDLT solver
  arcLength->options().setInt("BifurcationMethod",0); // 0: determinant, 1: eigenvalue
  arcLength->options().setReal("Length",dLb);
  // arcLength->options().setInt("AngleMethod",0); // 0: step, 1: iteration
  arcLength->options().setSwitch("AdaptiveLength",adaptive);
  arcLength->options().setInt("AdaptiveIterations",5);
  arcLength->options().setReal("Scaling",0.0);
  arcLength->options().setReal("Tol",ALM_tol);
  arcLength->options().setReal("TolU",ALM_tolU);
  arcLength->options().setReal("TolF",ALM_tolF);
  arcLength->options().setInt("MaxIter",maxit);
  arcLength->options().setSwitch("Verbose",(verbose>0));
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
  real_t Lold;
  gsVector<> Uold;

  typedef std::pair<gsVector<real_t>,real_t> solution_t;

  /*
    \a solutions is a container the for each level contains the solutions per point
    \a points contains all the points across levels in the format (level, U, lambda) -------------------------------> OVERKILL? WHY NEEDED?
    \a refPoints is a container that contains (level, U, lambda) of the points from which a refinement should START in level+1
    \a errors is a container that contains the error[l][i] e_i at the ith point of level l
  */

  std::vector<solution_t> solutions;
  std::vector<real_t> times;
  std::vector<index_t> levels;

  index_t level = 0;
  gsInfo<<"------------------------------------------------------------------------------------\n";
  gsInfo<<"\t\t\tLevel "<<level<<" (dL = "<<dLb<<") -- Coarse grid \n";
  gsInfo<<"------------------------------------------------------------------------------------\n";

  gsVector<> refPoints(2);
  gsVector<index_t> refPatches(1);
  refPoints<<1.0,0.5;
  refPatches<<0;

  gsAPALMData<real_t,solution_t> apalmData;
  apalmData.options().setInt("MaxLevel",maxLevel);
  apalmData.options().setInt("Verbose",verbose);
  apalmData.options().setReal("Tolerance",APALM_tol);

  gsAPALMBeam<real_t> apalm(comm,arcLength,apalmData,assembler,dirname,refPoints,refPatches);
  apalm.options().setSwitch("Verbose",(verbose>0));
  apalm.options().setInt("SubIntervals",SubIntervals);
  apalm.options().setSwitch("SingularPoint",true);
  apalm.options().setReal("BranchLengthMultiplier",dL/dLb);
  apalm.options().setReal("BifLengthMultiplier",0.05);
  apalm.initialize();

  if (!sequential)
  {
    real_t time = mpi.wallTime();
    apalm.solve(step+1);
    time = mpi.wallTime() - time;
    if (apalm.isMain()) gsInfo<<"Time = "<<time<<"\n";

    if (apalm.isMain())
    {
      if (plot || write)
      {
        gsField<> solField;
        gsMultiPatch<> mp_tmp;
        std::vector<std::string> pointheaders = {"u_x","u_y","u_z"};
        std::vector<std::string> otherheaders = {"U-norm","lambda","time","level"};

        for (index_t b=0; b!=apalm.getHierarchy().nBranches(); b++)
        {
          solutions = apalm.getFlatSolutions(b);
          times     = apalm.getFlatTimes(b);
          levels    = apalm.getFlatLevels(b);

          gsParaviewCollection dataCollection(dirname + "/" + "data_branch"+std::to_string(b));
          gsStructuralAnalysisOutput<real_t> data(dirname + "/data_branch"+std::to_string(b)+".csv",refPoints);
          if (write)
              data.init(pointheaders,otherheaders);
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
                  gsVector<> otherData(4);
                  otherData<<Uold.norm(),Lold,times[k],levels[k];
                  gsMatrix<> pointResults(mp.geoDim(),refPoints.cols());
                  for (index_t p=0; p!=refPoints.cols(); p++)
                      pointResults.col(p) = solField.value(refPoints.col(p),refPatches(0,p));
                  data.add(pointResults,otherData);
              }

            if (plot)
            {
              std::string fileName = dirname + "/" + "data_branch"+std::to_string(b) + util::to_string(k);
              gsWriteParaview<>(solField, fileName, 1000,mesh);
              fileName = "data_branch"+std::to_string(b) + util::to_string(k) + "0";
              dataCollection.addPart(fileName + ".vts",times[k]);
              if (mesh) dataCollection.addPart(fileName + "_mesh.vtp",times[k]);
            }
          }
          if (plot)
          {
            dataCollection.save();
          }
        }
      }
    }
    if (apalm.isMain())
    {
      std::ofstream file;
      file.open(dirname + "/times.txt");
      file<<"solution time: "<<time<<" s\n";
      file.close();
    }
  }
  else
  {
    real_t serialTime = mpi.wallTime();
    apalm.serialSolve(step+1);
    serialTime = mpi.wallTime() - serialTime;
    if (apalm.isMain()) gsInfo<<"Serial time = "<<serialTime<<"\n";

    if (apalm.isMain())
    {
      if (plot || write)
      {
        gsField<> solField;
        gsMultiPatch<> mp_tmp;
        std::vector<std::string> pointheaders = {"u_x","u_y","u_z"};
        std::vector<std::string> otherheaders = {"U-norm","lambda","time","level"};

        for (index_t b=0; b!=apalm.getHierarchy().nBranches(); b++)
        {
          solutions = apalm.getFlatSolutions(b);
          times     = apalm.getFlatTimes(b);
          levels    = apalm.getFlatLevels(b);

          gsParaviewCollection dataCollection(dirname + "/" + "data_serial_branch"+std::to_string(b));
          gsStructuralAnalysisOutput<real_t> data(dirname + "/data_serial_branch"+std::to_string(b)+".csv",refPoints);
          if (write)
              data.init(pointheaders,otherheaders);
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
                  gsVector<> otherData(4);
                  otherData<<Uold.norm(),Lold,times[k],levels[k];
                  gsMatrix<> pointResults(mp.geoDim(),refPoints.cols());
                  for (index_t p=0; p!=refPoints.cols(); p++)
                      pointResults.col(p) = solField.value(refPoints.col(p),refPatches(0,p));
                  data.add(pointResults,otherData);
              }

            if (plot)
            {
              std::string fileName = dirname + "/" + "data_serial_branch"+std::to_string(b) + util::to_string(k);
              gsWriteParaview<>(solField, fileName, 1000,mesh);
              fileName = "data_serial_branch"+std::to_string(b) + util::to_string(k) + "0";
              dataCollection.addPart(fileName + ".vts",times[k]);
              if (mesh) dataCollection.addPart(fileName + "_mesh.vtp",times[k]);
            }
          }
          if (plot)
          {
            dataCollection.save();
          }
        }
      }
    }



    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    real_t parallelTime = mpi.wallTime();
    apalm.parallelSolve();
    parallelTime = mpi.wallTime() - parallelTime;
    if (apalm.isMain()) gsInfo<<"Parallel time = "<<parallelTime<<"\n";

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    if (apalm.isMain())
    {
      if (plot || write)
      {
        gsField<> solField;
        gsMultiPatch<> mp_tmp;
        std::vector<std::string> pointheaders = {"u_x","u_y","u_z"};
        std::vector<std::string> otherheaders = {"U-norm","lambda","time","level"};

        for (index_t b=0; b!=apalm.getHierarchy().nBranches(); b++)
        {
          solutions = apalm.getFlatSolutions(b);
          times     = apalm.getFlatTimes(b);
          levels    = apalm.getFlatLevels(b);

          gsParaviewCollection dataCollection(dirname + "/" + "data_parallel_branch"+std::to_string(b));
          gsStructuralAnalysisOutput<real_t> data(dirname + "/data_parallel_branch"+std::to_string(b)+".csv",refPoints);
          if (write)
              data.init(pointheaders,otherheaders);
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
                  gsVector<> otherData(4);
                  otherData<<Uold.norm(),Lold,times[k],levels[k];
                  gsMatrix<> pointResults(mp.geoDim(),refPoints.cols());
                  for (index_t p=0; p!=refPoints.cols(); p++)
                      pointResults.col(p) = solField.value(refPoints.col(p),refPatches(0,p));
                  data.add(pointResults,otherData);
              }

            if (plot)
            {
              std::string fileName = dirname + "/" + "data_parallel_branch"+std::to_string(b) + util::to_string(k);
              gsWriteParaview<>(solField, fileName, 1000,mesh);
              fileName = "data_parallel_branch"+std::to_string(b) + util::to_string(k) + "0";
              dataCollection.addPart(fileName + ".vts",times[k]);
              if (mesh) dataCollection.addPart(fileName + "_mesh.vtp",times[k]);
            }
          }
          if (plot)
          {
            dataCollection.save();
          }
        }
      }
    }
    if (apalm.isMain())
    {
      std::ofstream file;
      file.open(dirname + "/times.txt");
      file<<"serial   time: "<<serialTime<<" s\n";
      file<<"parallel time: "<<parallelTime<<" s\n";
      file.close();
    }
  }

  delete assembler;
  delete materialMatrix;
  delete arcLength;
  return EXIT_SUCCESS;
}
#else//gsKLShell_ENABLED
int main(int argc, char *argv[])
{
    gsWarn<<"G+Smo is not compiled with the gsKLShell module.";
    return EXIT_FAILURE;
}
#endif

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
void addClamping(gsMultiPatch<T>& mp, index_t patch, std::vector<boxSide> sides, T offset) //, std::vector<boxSide> sides, T offset)
{

    gsTensorBSpline<2,T> *geo = dynamic_cast< gsTensorBSpline<2,real_t> * > (&mp.patch(patch));

    T dknot0 = geo->basis().component(0).knots().minIntervalLength();
    T dknot1 = geo->basis().component(1).knots().minIntervalLength();

    gsInfo<<"sides.size() = "<<sides.size()<<"\n";

    index_t k =0;


    for (std::vector<boxSide>::iterator it = sides.begin(); it != sides.end(); it++)
    {
        gsInfo<<"side = "<<(*it)<<"\n";

      if (*it==boundary::west || *it==boundary::east) // west or east
      {
        if (*it==boundary::east) // east, val = 1
          geo->insertKnot(1 - std::min(offset, dknot0 / 2),0);
        else if (*it==boundary::west) // west
          geo->insertKnot(std::min(offset, dknot0 / 2),0);
      }
      else if (*it==boundary::south || *it==boundary::north) // west or east
      {
       if (*it==boundary::north) // north
         geo->insertKnot(1 - std::min(offset, dknot0 / 2),1);
       else if (*it==boundary::south) // south
         geo->insertKnot(std::min(offset, dknot0 / 2),1);
      }
      else if (*it==boundary::none)
        gsWarn<<*it<<"\n";
      else
        GISMO_ERROR("Side unknown, side = " <<*it);

        k++;
gsInfo<<"k = "<<k<<"\n";
    }
}

template <class T>
gsMultiPatch<T> Rectangle(T L, T B) //, int n, int m, std::vector<boxSide> sides, T offset)
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
