/** @file benchmark_Elasticity_Beam_APALM.cpp

    @brief Code for arc-length analysis of a solid beam using the APALM

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#include <gismo.h>

#ifdef gsElasticity_ENABLED
#include <gsElasticity/gsGeoUtils.h>
#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>
#endif

#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMBase.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMLoadControl.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMRiks.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMCrisfield.h>

#include <gsStructuralAnalysis/src/gsALMSolvers/gsAPALMData.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsAPALM.h>
#include <gsStructuralAnalysis/src/gsStructuralAnalysisTools/gsStructuralAnalysisUtils.h>

using namespace gismo;

template <class T>
gsMultiPatch<T> BrickDomain(int n, int m, int o, int p, int q ,int r, T L, T B, T H);
template <class T>
gsMultiPatch<T> BrickDomain(int n, int p, T L, T B, T H);

#ifdef gsElasticity_ENABLED
template<class T>
class gsAPALMBeam : public gsAPALM<T>
{

  using Base = gsAPALM<T>;
  typedef typename Base::solution_t solution_t;

public:
  gsAPALMBeam(
              const gsMpiComm & comm,
              gsALMBase<T> * ALM,
              const gsAPALMData<T,solution_t> & Data,
              const gsElasticityAssembler<T> & assembler,
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
    gsMultiPatch<T> deformation,mp_tmp, mp, mp_def;
    std::vector<std::string> pointheaders = {"u_x","u_y","u_z"};
    std::vector<std::string> otherheaders = {"U-norm","lambda","time","level"};

    data.init(pointheaders,otherheaders);

    mp = m_assembler.patches();
    std::vector<gsMatrix<> > fixedDofs = m_assembler.allFixedDofs();
    for (size_t k=0; k!=stepSolutions.size(); k++)
    {
      m_assembler.constructSolution(stepSolutions[k].first,fixedDofs,mp_def);
      gsField<T> solField(m_assembler.patches(),mp_def);

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
  const gsElasticityAssembler<T> m_assembler;

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
  bool make_fit     = false;
  bool mesh         = false;
  bool stress       = false;
  bool membrane     = false;
  bool quasiNewton  = false;
  int quasiNewtonInt= -1;
  bool adaptive     = false;
  int step          = 4;
  int SubIntervals  = 2;
  int method        = 2; // (0: Load control; 1: Riks' method; 2: Crisfield's method; 3: consistent crisfield method; 4: extended iterations)
  bool deformed     = false;

  bool composite = false;

  real_t relax      = 1.0;

  int testCase      = 3;

  int result        = 0;

  bool write        = false;

  index_t maxit     = 20;
  // index_t iniLevels  = 2;
  // index_t maxLevels  = 4;
  index_t maxLevel  = 6;

  // Arc length method options
  real_t dL         = 5e0 ; // General arc length
  real_t dLb        = 5e-1; // Arc length to find bifurcation
  real_t tol        = 1e-6;
  real_t tolU       = 1e-6;
  real_t tolF       = 1e-3;


  index_t verbose = 0;


  std::string wn("data.csv");

  std::string assemberOptionsFile("options/solver_options.xml");

  gsCmdLine cmd("Arc-length analysis for thin shells.");
  cmd.addString( "f", "file", "Input XML file for assembler options", assemberOptionsFile );

  cmd.addInt("t", "testcase", "Test case: 0: clamped-clamped, 1: pinned-pinned, 2: clamped-free", testCase);

  cmd.addInt("r","hRefine", "Number of dyadic h-refinement (bisection) steps to perform before solving", numHref);
  cmd.addInt("e","degreeElevation", "Number of degree elevation steps to perform on the Geometry's basis before solving", numElevate);
  cmd.addSwitch("composite", "Composite material", composite);

  cmd.addInt("m","Method", "Arc length method; 1: Crisfield's method; 2: RIks' method.", method);
  cmd.addReal("L","dLb", "arc length", dLb);
  cmd.addReal("l","dL", "arc length after bifurcation", dL);
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
  cmd.addSwitch("fit", "Make a fit and refer to that", make_fit);

  try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

  if (dL==0)
  {
    dL = dLb;
  }

  gsFileData<> fd(assemberOptionsFile);
  gsOptionList opts;
  fd.getFirst<gsOptionList>(opts);

  gsMultiPatch<> mp;
  real_t L,B,H;

  real_t E_modulus, PoissonRatio;
  /*
    Case 2: Clamped beam (left) under vertical end load                   --- Validation settings: -L 5e-1 -M 0 -N 10 -r 2 -e 1 (--plot --write -q 5)
                                                                              Fig 3b from: Pagani, A., & Carrera, E. (2018). Unified formulation of geometrically nonlinear refined beam theories. Mechanics of Advanced Materials and Structures, 25(1), 15–31. https://doi.org/10.1080/15376494.2016.1232458
    Case 3: Clamped beam (left) under horizontal compressive end load     --- Validation settings: -L 5e-5 -l 1e-1 -M 0 -N 100 -r 4 -e 2
                                                                              Fig 5  from: Pagani, A., & Carrera, E. (2018). Unified formulation of geometrically nonlinear refined beam theories. Mechanics of Advanced Materials and Structures, 25(1), 15–31. https://doi.org/10.1080/15376494.2016.1232458
  */
  E_modulus = 1.;
  H = 0.01;
  PoissonRatio = 0.0;
  L = 1.0;
  B = 0.01;
  mp = BrickDomain(numHref, 0, 0, numElevate, 1, 1, L, B,H);

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

  std::string dirname = "ArcLengthResultsElasticity";
  std::string output = "solution";

  gsMatrix<> writePoints(3,3);
  writePoints.col(0)<< 0.0,0.5,0.5;
  writePoints.col(1)<< 0.5,0.5,0.5;
  writePoints.col(2)<< 1.0,0.5,0.5;

  if (testCase == 2)
  {
    real_t Load = 1e1;
    tmp << 0, 0, Load/(B*H);
    neuData.setValue(tmp,3);

    BCs.addCondition(0,boundary::west,condition_type::dirichlet,nullptr,0); // last number is a component (coordinate) number
    BCs.addCondition(0,boundary::west,condition_type::dirichlet,nullptr,1);
    BCs.addCondition(0,boundary::west,condition_type::dirichlet,nullptr,2);


    // BCs.addCondition(0,boundary::east,condition_type::dirichlet,nullptr,0);
    // BCs.addCondition(0,boundary::east,condition_type::dirichlet,nullptr,1);
    // vertical load
    BCs.addCondition(0,boundary::east,condition_type::neumann,&neuData);

    // dL =  1e-3;
    // dLb = 2e0;

    dirname = dirname + "/Beam_clamped-verticalLoad";
    output =  "solution";
    wn = output + "data.txt";
  }
  else if (testCase == 3)
  {
    real_t Load = 1e1;
    tmp << -Load/(B*H), 0, 0;
    gsInfo<<"Applied "<<tmp.transpose()<<" as load vector\n";
    neuData.setValue(tmp,3);

    BCs.addCondition(0,boundary::west,condition_type::dirichlet,nullptr,0); // last number is a component (coordinate) number
    BCs.addCondition(0,boundary::west,condition_type::dirichlet,nullptr,1);
    BCs.addCondition(0,boundary::west,condition_type::dirichlet,nullptr,2);

    BCs.addCondition(0,boundary::north,condition_type::dirichlet,nullptr,1);
    BCs.addCondition(0,boundary::south,condition_type::dirichlet,nullptr,1);

    // BCs.addCondition(0,boundary::east,condition_type::dirichlet,nullptr,0);
    // BCs.addCondition(0,boundary::east,condition_type::dirichlet,nullptr,1);
    // vertical load
    BCs.addCondition(0,boundary::east,condition_type::neumann,&neuData);

    dirname = dirname + "/Beam_clamped-horizontalLoad";
    output =  "solution";
    wn = output + "data.txt";
  }

  std::string commands = "mkdir -p " + dirname;
  const char *command = commands.c_str();
  system(command);

  gsFunctionExpr<> g("0","0","0",3);

  gsElasticityAssembler<real_t> assembler(mp,dbasis,BCs,g);
  assembler.options().setReal("YoungsModulus",E_modulus);
  assembler.options().setReal("PoissonsRatio",PoissonRatio);
  assembler.options().setInt("MaterialLaw",material_law::hooke);
  gsInfo << "Initialized system with " << assembler.numDofs() << " dofs.\n";

  gsMultiPatch<> mp_def;

  gsStopwatch stopwatch;
  real_t time = 0.0;

  // Assemble linear system to obtain the force vector
  assembler.assemble();
  gsVector<> Force = assembler.rhs();

  std::vector<gsMatrix<> > fixedDofs = assembler.allFixedDofs();
  // Function for the Jacobian
  gsStructuralAnalysisOps<real_t>::Jacobian_t Jacobian = [&time,&stopwatch,&assembler,&fixedDofs](gsVector<real_t> const &x, gsSparseMatrix<real_t> &m)
  {
    stopwatch.restart();
    assembler.assemble(x,fixedDofs);
    time += stopwatch.stop();

    m = assembler.matrix();
    // gsInfo<<"matrix = \n"<<m.toDense()<<"\n";
    return true;
  };
  // Function for the Residual
  gsStructuralAnalysisOps<real_t>::ALResidual_t ALResidual = [&time,&stopwatch,&assembler,&fixedDofs,&Force](gsVector<real_t> const &x, real_t lam, gsVector<real_t> &result)
  {
    stopwatch.restart();
    assembler.assemble(x,fixedDofs);
    result = Force - lam * Force - assembler.rhs(); // assembler rhs - force = Finternal
    time += stopwatch.stop();
    return true;
  };

  // if (material==0)
    assembler.options().setInt("MaterialLaw",material_law::saint_venant_kirchhoff);
  // else if (material==2)
  //   assembler.options().setInt("MaterialLaw",material_law::neo_hooke_ln);
  // else if (material==22)
  //   assembler.options().setInt("MaterialLaw",material_law::neo_hooke_quad);

  gsALMBase<real_t> * arcLength;
  if (method==0)
    arcLength = new gsALMLoadControl<real_t>(Jacobian, ALResidual, Force);
  else if (method==1)
    arcLength = new gsALMRiks<real_t>(Jacobian, ALResidual, Force);
  else if (method==2)
    arcLength = new gsALMCrisfield<real_t>(Jacobian, ALResidual, Force);
  else
    GISMO_ERROR("Method "<<method<<" unknown");

  arcLength->options().setString("Solver","SimplicialLDLT"); // LDLT solver
  arcLength->options().setInt("BifurcationMethod",0); // 0: determinant, 1: eigenvalue
  arcLength->options().setReal("Length",dLb);
  // arcLength->options().setInt("AngleMethod",0); // 0: step, 1: iteration
  arcLength->options().setSwitch("AdaptiveLength",adaptive);
  arcLength->options().setInt("AdaptiveIterations",5);
  arcLength->options().setReal("Scaling",1.0);
  arcLength->options().setReal("Tol",tol);
  arcLength->options().setReal("TolU",tolU);
  arcLength->options().setReal("TolF",tolF);
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

  gsVector<> refPoints(3);
  gsVector<index_t> refPatches(1);
  refPoints<<1.0,0.5,0.5;
  refPatches<<0;

  gsAPALMData<real_t,solution_t> apalmData;
  apalmData.options().setInt("MaxLevel",maxLevel);
  apalmData.options().setInt("Verbose",verbose);
  apalmData.options().setReal("Tolerance",1e-3);

  const gsMpi & mpi = gsMpi::init(argc, argv);
  gsMpiComm comm = mpi.worldComm();

  gsAPALMBeam<real_t> apalm(comm,arcLength,apalmData,assembler,dirname,refPoints,refPatches);
  apalm.options().setSwitch("Verbose",(verbose>0));
  apalm.options().setInt("SubIntervals",SubIntervals);
  apalm.options().setSwitch("SingularPoint",true);
  apalm.options().setReal("BranchLengthMultiplier",dL/dLb);
  apalm.options().setReal("BifLengthMultiplier",0.05*dLb/dL);
  apalm.initialize();
  apalm.serialSolve(step+1);

  if (apalm.isMain())
  {
    if (plot || write)
    {
      gsField<> solField, stressField;
      gsMultiPatch<> mp_tmp;
      std::vector<std::string> pointheaders = {"u_x","u_y","u_z"};
      std::vector<std::string> otherheaders = {"U-norm","lambda","time","level"};
      gsPiecewiseFunction<> stresses;

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

          assembler.constructSolution(Uold,fixedDofs,mp_def);
          solField = gsField<>(assembler.patches(),mp_def);
          assembler.constructCauchyStresses(mp_def,stresses,stress_components::von_mises);
          stressField = gsField<>(assembler.patches(),stresses,true);

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
            std::map<std::string,const gsField<> *> fields;
            fields["Deformation"] = &solField;
            fields["Stress"] = &stressField;
            gsWriteParaviewMultiPhysics(fields,fileName,1000,true);
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

  apalm.parallelSolve();

  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////

  if (apalm.isMain())
  {
    if (plot || write)
    {
      gsField<> solField, stressField;
      gsMultiPatch<> mp_tmp;
      std::vector<std::string> pointheaders = {"u_x","u_y","u_z"};
      std::vector<std::string> otherheaders = {"U-norm","lambda","time","level"};
      gsPiecewiseFunction<> stresses;

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

          assembler.constructSolution(Uold,fixedDofs,mp_def);
          solField = gsField<>(assembler.patches(),mp_def);
          assembler.constructCauchyStresses(mp_def,stresses,stress_components::von_mises);
          stressField = gsField<>(assembler.patches(),stresses,true);

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
            std::map<std::string,const gsField<> *> fields;
            fields["Deformation"] = &solField;
            fields["Stress"] = &stressField;
            gsWriteParaviewMultiPhysics(fields,fileName,1000,true);
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

  delete arcLength;
  return result;
}
#else//gsElasticity_ENABLED
int main(int argc, char *argv[])
{
    gsWarn<<"G+Smo is not compiled with the gsElasticity module.";
    return EXIT_FAILURE;
}
#endif

template <class T>
gsMultiPatch<T> BrickDomain(int n, int p, T L, T B, T H)
{
  int q = p;
  int r = p;

  int m = n;
  int o = n;
  gsMultiPatch<T> mp = BrickDomain(n, m, o, p, q, r, L, B, H);
  return mp;
}

template <class T>
gsMultiPatch<T> BrickDomain(int n, int m, int o, int p, int q ,int r, T L, T B, T H)
{
  // -------------------------------------------------------------------------
  // --------------------------Make beam geometry-----------------------------
  // -------------------------------------------------------------------------
  int dim = 3; //physical dimension
  gsKnotVector<> kv0;
  kv0.initUniform(0,1,0,p+1,1);
  gsKnotVector<> kv1;
  kv1.initUniform(0,1,0,q+1,1);
  gsKnotVector<> kv2;
  kv2.initUniform(0,1,0,r+1,1);

  for(index_t i = 0; i< n; ++i)
      kv0.uniformRefine();
  for(index_t i = 0; i< m; ++i)
      kv1.uniformRefine();
  for(index_t i = 0; i< o; ++i)
      kv2.uniformRefine();

  // Make basis
  gsTensorBSplineBasis<3,T> basis(kv0,kv1,kv2);

  // Initiate coefficient matrix
  gsMatrix<> coefs(basis.size(),dim);
  // Number of control points needed per component
  size_t len0 = basis.component(0).size();
  size_t len1 = basis.component(1).size();
  size_t len2 = basis.component(2).size();
  // Uniformly distribute control points per component
  gsVector<> coefvec0(len0);
  coefvec0.setLinSpaced(len0,0.0,L);
  gsVector<> coefvec1(basis.component(1).size());
  coefvec1.setLinSpaced(len1,0.0,B);
  gsVector<> coefvec2(basis.component(2).size());
  coefvec2.setLinSpaced(len2,0.0,H);

  // Define a matrix with ones
  gsVector<> temp(len0);
  temp.setOnes();
  for (size_t l = 0; l < len2; l++)
    {
        for (size_t k = 0; k < len1; k++)
        {
            index_t offset = l*len0*len1;
            // First column contains x-coordinates (length)
            coefs.col(0).segment(k*len0+offset,len0) = coefvec0;
            // Second column contains y-coordinates (width)
            coefs.col(1).segment(k*len0+offset,len0) = temp*coefvec1.at(k);

            coefs.col(2).segment(k*len0+offset,len0) = temp*coefvec2.at(l);
        }
    }
  // gsInfo<<"\n"<<coefs<<"\n";
  // Create gsGeometry-derived object for the patch
  gsTensorBSpline<3,real_t> shape(basis,coefs);

  gsMultiPatch<T> mp;
  mp.addPatch(shape);
  mp.addAutoBoundaries();

  return mp;
}
