/** @file benchmark_Frustrum_APALM.cpp

    @brief Benchmark for the collapsing frustrum with the APALM

    Based on:
    Başar, Y., & Itskov, M. (1998). Finite element formulation of the Ogden material model with application to ruber-like shells.
    International Journal for Numerical Methods in Engineering. https://doi.org/10.1002/(SICI)1097-0207(19980815)42:7<1279::AID-NME437>3.0.CO;2-I

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
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMLoadControl.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsAPALMData.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsAPALM.h>
#include <gsStructuralAnalysis/src/gsStructuralAnalysisTools/gsStructuralAnalysisUtils.h>
using namespace gismo;

template <class T>
gsMultiPatch<T> FrustrumDomain(int n, int p, T R1, T R2, T h);

template <class T>
void initStepOutput( const std::string name, const gsMatrix<T> & points);

template <class T>
void writeStepOutput(const gsALMBase<T> * arcLength, const gsMultiPatch<T> & deformation, const std::string name, const gsMatrix<T> & points, const index_t extreme=-1, const index_t kmax=100);

#ifdef gsKLShell_ENABLED
template<class T>
class gsAPALMFrustrum : public gsAPALM<T>
{

  using Base = gsAPALM<T>;
  typedef typename Base::solution_t solution_t;

public:
  gsAPALMFrustrum(const gsMpiComm & comm,
                  gsALMBase<T> * ALM,
                  const gsAPALMData<T,solution_t> & Data,
                  const gsThinShellAssemblerBase<T> * assembler,
                  std::string dirname,
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
    std::vector<std::string> otherheaders = {"lambda","level"};

    data.init(pointheaders,otherheaders);

    for (size_t k=0; k!=stepSolutions.size(); k++)
    {
      m_assembler->constructSolution(stepSolutions[k].first,mp_tmp);
      deformation.patch(0) = mp_tmp.patch(0);
      deformation.patch(0).coefs() -= mp.patch(0).coefs();// assuming 1 patch here
      gsField<T> solField(mp,deformation);

      if (m_refPoints.cols()!=0)
      {
          gsVector<> otherData(2);
          otherData<<stepSolutions[k].second,level;
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
    int numElevate    = 1;
    int numHref       = 5;
    bool plot         = false;
    bool sequential   = false;
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

    bool composite = false;
    index_t material  = 3;
    index_t impl = 1; // 1= analytical, 2= generalized, 3= spectral

    real_t relax      = 1.0;

    int testCase      = 1;

    int result        = 0;

    bool write        = false;

    index_t maxit     = 20;
    // index_t iniLevels  = 2;
    // index_t maxLevels  = 4;
    index_t maxLevel  = 2;

    // Arc length method options
    real_t dL         = -1; // General arc length
    real_t ALM_tol    = 1e-6;
    real_t ALM_tolU   = 1e-6;
    real_t ALM_tolF   = 1e-3;

    real_t APALM_tol  = 1e-2;


    index_t verbose = 0;

    std::string wn("data.csv");

    std::string assemberOptionsFile("options/solver_options.xml");

    gsCmdLine cmd("APALM analysis of a collapsing frustrum.");
    cmd.addString( "f", "file", "Input XML file for assembler options", assemberOptionsFile );

    cmd.addInt("t", "testcase", "Test case: 0: Constant top radius, 1: variable top radius", testCase);

    cmd.addInt("r","hRefine", "Number of dyadic h-refinement (bisection) steps to perform before solving", numHref);
    cmd.addInt("e","degreeElevation", "Number of degree elevation steps to perform on the Geometry's basis before solving", numElevate);

    cmd.addInt( "M", "Material", "Material law",  material );
    cmd.addInt( "I", "Implementation", "Implementation: 1= analytical, 2= generalized, 3= spectral",  impl );
    cmd.addSwitch("composite", "Composite material", composite);

    cmd.addInt("m","Method", "Arc length method; 1: Crisfield's method; 2: Riks' method.", method);
    cmd.addReal("L","dL", "arc length", dL);
    cmd.addReal("T", "APALM_tol", "APALM Tolerance", APALM_tol);
    // cmd.addInt("I","inilvl", "Initial levels", iniLevels);
    // cmd.addInt("M","maxlvl", "Max levels", maxLevels);
    cmd.addInt("l","level", "Max level", maxLevel);
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
    cmd.addSwitch("sequential", "Solve sequentially, i.e. serial --> parallel", sequential);


    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    real_t mu = 4.225;
    real_t thickness = 0.1;
    real_t PoissonRatio;
    if (material==0)
      PoissonRatio = 0.499;
    else
      PoissonRatio = 0.5;
    real_t E_modulus = 2*mu*(1+PoissonRatio);
    real_t Density    = 1e0;
    real_t Ratio      = 7.0;

    gsMultiPatch<> mp,mp_def;
    mp = FrustrumDomain(numHref,numElevate+2,2.0,1.0,1.0);
    // Initialise solution object
    mp_def = mp;

    gsMultiBasis<> dbasis(mp,true);
    gsInfo<<"Basis (patch 0): "<< mp.patch(0).basis() << "\n";

    // Boundary conditions
    gsBoundaryConditions<> BCs;
    BCs.setGeoMap(mp);

    real_t Load = -1;
    gsVector<> neu(3);
    neu << 0, 0, Load;
    gsConstantFunction<> neuData(neu,3);

    std::string output = "solution";
    std::string dirname = "ArcLengthResults";
    std::string cores;

    gsMatrix<> writePoints(2,3);
    writePoints.col(0)<<0.0,1.0;
    writePoints.col(1)<<0.5,1.0;
    writePoints.col(2)<<1.0,1.0;

    /*
      Case 0: Frustrum with constrained top boundary        (Complete cycle with -N 1200)
      Case 1: Frustrum with unconstrained top boundary      (Complete cycle with -N 850)
    */
    const gsMpi & mpi = gsMpi::init();
    gsMpiComm comm = mpi.worldComm();
    cores = "_ncores="+std::to_string(comm.size());

    if (testCase == 0)
    {
      BCs.addCondition(boundary::north, condition_type::neumann, &neuData );
      BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 2 - z
      BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 2 - z
      BCs.addCondition(boundary::north, condition_type::collapsed, 0, 0, false, 2 ); // unknown 1 - y

      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

      // Symmetry in x-direction:
      BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 );
      BCs.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 1 );
      BCs.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 2 );

      // Symmetry in y-direction:
      BCs.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 0 );
      BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 );
      BCs.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 2 );

      dirname = dirname + "/" + "Frustrum_-r=" + std::to_string(numHref) + "-e" + std::to_string(numElevate) + "-M" + std::to_string(material) + "_solution" + cores;
      if (dL == -1 && material!=3) { dL = 5e-2; }
      if (dL == -1 && material==3) { dL = 4e-2; }
    }
    else if (testCase == 1)
    {
      BCs.addCondition(boundary::north, condition_type::neumann, &neuData );
      BCs.addCondition(boundary::north, condition_type::collapsed, 0, 0, false, 2 ); // unknown 1 - y

      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

      // Symmetry in x-direction:
      BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 );
      BCs.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 1 );
      BCs.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 2 );

      // Symmetry in y-direction:
      BCs.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 0 );
      BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 );
      BCs.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 2 );

      dirname = dirname + "/" + "Frustrum2_-r=" + std::to_string(numHref) + "-e" + std::to_string(numElevate) + "-M" + std::to_string(material) + "_solution" + cores;
      if (dL == -1) { dL = 4e-2; }
    }
    else
      GISMO_ERROR("Test case" << testCase << " does not exist!");

    wn = output + "data.txt";

    // Prepare and create directory with dirname
    dirname = gsFileManager::getCurrentPath() + dirname;
    if (sequential)
      dirname = dirname + "_seq";

    GISMO_ENSURE(gsFileManager::mkdir(dirname),"Failed to create directory " + dirname);
    // Made directory

    // Linear isotropic material model
    gsFunctionExpr<> force("0","0","0",3);
    gsConstantFunction<> t(thickness,3);
    gsConstantFunction<> E(E_modulus,3);
    gsConstantFunction<> nu(PoissonRatio,3);
    gsConstantFunction<> rho(Density,3);
    gsConstantFunction<> ratio(Ratio,3);

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

    std::vector<gsFunctionSet<>*> parameters;
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
            options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
            options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
            materialMatrix = getMaterialMatrix<3,real_t>(mp,t,parameters,rho,options);
        }
    }
    else
    {
        options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",material);
        options.addSwitch("Compressibility","Compressibility: (false): Imcompressible | (true): Compressible",false);
        options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",impl);
        materialMatrix = getMaterialMatrix<3,real_t>(mp,t,parameters,rho,options);
    }

    gsThinShellAssemblerBase<real_t>* assembler;
    assembler = new gsThinShellAssembler<3, real_t, true >(mp,dbasis,BCs,force,materialMatrix);


    gsStopwatch stopwatch;
    real_t time = 0.0;

    // Assemble linear system to obtain the force vector
    assembler->assemble();
    gsVector<> Force = assembler->rhs();
    
    // Function for the Jacobian
    gsStructuralAnalysisOps<real_t>::Jacobian_t Jacobian = [&time,&stopwatch,&assembler,&mp_def](gsVector<real_t> const &x, gsSparseMatrix<real_t> & m)
    {
      ThinShellAssemblerStatus status;
      stopwatch.restart();
      assembler->constructSolution(x,mp_def);
      status = assembler->assembleMatrix(mp_def);
      m = assembler->matrix();
      time += stopwatch.stop();
      return status == ThinShellAssemblerStatus::Success;
    };
    // Function for the Residual
    gsStructuralAnalysisOps<real_t>::ALResidual_t ALResidual = [&time,&stopwatch,&assembler,&mp_def,&Force](gsVector<real_t> const &x, real_t lam, gsVector<real_t> & result)
    {
        ThinShellAssemblerStatus status;
        stopwatch.restart();
        assembler->constructSolution(x,mp_def);
        status = assembler->assembleVector(mp_def);
        result = Force - lam * Force - assembler->rhs(); // assembler rhs - force = Finternal
        time += stopwatch.stop();
        return status == ThinShellAssemblerStatus::Success;
    };

    gsALMBase<real_t> * arcLength;
    if (method==0)
      arcLength = new gsALMLoadControl<real_t>(Jacobian, ALResidual, Force);
    else if (method==1)
      arcLength = new gsALMRiks<real_t>(Jacobian, ALResidual, Force);
    else if (method==2)
      arcLength = new gsALMCrisfield<real_t>(Jacobian, ALResidual, Force);
    else
      GISMO_ERROR("Method "<<method<<" unknown");

    arcLength->options().setString("Solver","CGDiagonal"); // LDLT solver
    arcLength->options().setInt("BifurcationMethod",1); // 0: determinant, 1: eigenvalue
    arcLength->options().setReal("Length",dL);
    arcLength->options().setInt("AngleMethod",0); // 0: step, 1: iteration
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


    gsDebug<<arcLength->options();
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

    gsVector<> refPoints(2);
    gsVector<index_t> refPatches(1);
    refPoints<<1.0,1.0;
    refPatches<<0;

    gsAPALMData<real_t,solution_t> apalmData;
    apalmData.options().setInt("MaxLevel",maxLevel);
    apalmData.options().setInt("Verbose",verbose);
    apalmData.options().setReal("Tolerance",APALM_tol);

    gsAPALMFrustrum<real_t> apalm(comm,arcLength,apalmData,assembler,dirname,refPoints,refPatches);

#   ifdef _OPENMP
    const int nt  = omp_get_num_threads();
    gsInfo << "Gismo was compiled with MPI support with "<<apalm.size()<<" ranks and with OpenMP support with "<<nt<<" threads per rank.\n";
#   else
    gsInfo << "Gismo was compiled with MPI support with "<<apalm.size()<<" cores.\n";
#   endif

    apalm.options().setSwitch("Verbose",(verbose>0));
    apalm.options().setInt("SubIntervals",SubIntervals);
    apalm.options().setInt("MaxIt",1000);
    apalm.initialize();


    if(!sequential)
    {
      real_t time = mpi.wallTime();
      apalm.solve(step+1);
      time = mpi.wallTime() - time;
      if (apalm.isMain()) gsInfo<<"Time = "<<time<<"\n";

      // plot geometry
      if (apalm.isMain())
        if (plot)
          gsWriteParaview(mp,dirname + "/" + "mp",1000,true);

      if (apalm.isMain())
        if (write)
          initStepOutput(dirname + "/" + wn, writePoints);

      if (apalm.isMain())
      {
        solutions = apalm.getFlatSolutions();
        times     = apalm.getFlatTimes();
        levels    = apalm.getFlatLevels();

        if (plot || write)
        {
          gsField<> solField;
          gsMultiPatch<> mp_tmp;
          std::vector<std::string> pointheaders = {"u_x","u_y","u_z"};
          std::vector<std::string> otherheaders = {"lambda","level"};

          gsParaviewCollection dataCollection(dirname + "/" + "data");
          gsStructuralAnalysisOutput<real_t> data(dirname + "/data.csv",refPoints);
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
                  gsVector<> otherData(2);
                  otherData<<Lold,levels[k];
                  gsMatrix<> pointResults(mp.geoDim(),refPoints.cols());
                  for (index_t p=0; p!=refPoints.cols(); p++)
                      pointResults.col(p) = solField.value(refPoints.col(p),refPatches(0,p));
                  data.add(pointResults,otherData);
              }

            if (plot)
            {
              std::string fileName = dirname + "/" + "data" + util::to_string(k);
              gsWriteParaview<>(solField, fileName, 1000,mesh);
              fileName = "data" + util::to_string(k) + "0";
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
      if (apalm.isMain())
      {
        gsInfo<<"------------------------------------------------------------------------------------\n";
        gsInfo<<"\t\t\tLevel "<<level<<" (dL = "<<dL<<") -- Coarse grid \n";
        gsInfo<<"------------------------------------------------------------------------------------\n";
      }
      apalm.serialSolve(step+1);
      serialTime = mpi.wallTime() - serialTime;
      if (apalm.isMain()) gsInfo<<"Serial time = "<<serialTime<<"\n";

      // plot geometry
      if (apalm.isMain())
        if (plot)
          gsWriteParaview(mp,dirname + "/" + "mp",1000,true);

      if (apalm.isMain())
        if (write)
          initStepOutput(dirname + "/" + wn, writePoints);

      if (apalm.isMain())
      {
        solutions = apalm.getFlatSolutions();
        times     = apalm.getFlatTimes();
        levels    = apalm.getFlatLevels();

        if (plot || write)
        {
          gsField<> solField;
          gsMultiPatch<> mp_tmp;
          std::vector<std::string> pointheaders = {"u_x","u_y","u_z"};
          std::vector<std::string> otherheaders = {"lambda","level"};

          gsParaviewCollection dataCollection(dirname + "/" + "data_serial");
          gsStructuralAnalysisOutput<real_t> data(dirname + "/data_serial.csv",refPoints);
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
                  gsVector<> otherData(2);
                  otherData<<Lold,levels[k];
                  gsMatrix<> pointResults(mp.geoDim(),refPoints.cols());
                  for (index_t p=0; p!=refPoints.cols(); p++)
                      pointResults.col(p) = solField.value(refPoints.col(p),refPatches(0,p));
                  data.add(pointResults,otherData);
              }

            if (plot)
            {
              std::string fileName = dirname + "/" + "data_serial" + util::to_string(k);
              gsWriteParaview<>(solField, fileName, 1000,mesh);
              fileName = "data_serial" + util::to_string(k) + "0";
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
        solutions = apalm.getFlatSolutions();
        times     = apalm.getFlatTimes();
        levels    = apalm.getFlatLevels();

        if (plot || write)
        {
          gsField<> solField;
          gsMultiPatch<> mp_tmp;
          std::vector<std::string> pointheaders = {"u_x","u_y","u_z"};
          std::vector<std::string> otherheaders = {"lambda","level"};

          gsParaviewCollection dataCollection2(dirname + "/" + "data_parallel");
          gsStructuralAnalysisOutput<real_t> data2(dirname + "/data_parallel.csv",refPoints);
          if (write)
              data2.init(pointheaders,otherheaders);

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
                gsVector<> otherData(2);
                otherData<<Lold,levels[k];
                gsMatrix<> pointResults(mp.geoDim(),refPoints.cols());
                for (index_t p=0; p!=refPoints.cols(); p++)
                    pointResults.col(p) = solField.value(refPoints.col(p),refPatches(0,p));
                data2.add(pointResults,otherData);
              }

            if (plot)
            {
              std::string fileName = dirname + "/" + "data_parallel" + util::to_string(k);
              gsWriteParaview<>(solField, fileName, 1000,mesh);
              fileName = "data_parallel" + util::to_string(k) + "0";
              dataCollection2.addPart(fileName + ".vts",times[k]);
              if (mesh) dataCollection2.addPart(fileName + "_mesh.vtp",times[k]);
            }
          }
          if (plot)
          {
            dataCollection2.save();
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
gsMultiPatch<T> FrustrumDomain(int n, int p, T R1, T R2, T h)
{
  // -------------------------------------------------------------------------
  // --------------------------Make beam geometry-----------------------------
  // -------------------------------------------------------------------------
  // n = number of uniform refinements over the height; n = 0, only top and bottom part

  int dim = 3; //physical dimension
  gsKnotVector<> kv0;
  kv0.initUniform(0,1,0,3,1);
  gsKnotVector<> kv1;
  kv1.initUniform(0,1,0,3,1);

  // Refine n times
  for(index_t i = 0; i< n; ++i)
      kv1.uniformRefine();

  // Make basis
  // gsTensorNurbsBasis<2,T> basis(kv0,kv1);

  // Initiate coefficient matrix
  index_t N = math::pow(2,n)+2;
  gsMatrix<> coefs(3*N,dim);
  gsMatrix<> tmp(3,3);
  T R,H;

  gsMatrix<> weights(3*N,1);
  for (index_t k=0; k!= N; k++)
  {
    R = k*(R2-R1)/(N-1) + R1;
    H = k*h/(N-1);
    tmp<< R,0,H,
          R,R,H,
          0,R,H;

    coefs.block(3*k,0,3,3) = tmp;

    weights.block(3*k,0,3,1) << 1,0.70711,1;
  }

  // Create gsGeometry-derived object for the patch
  gsTensorNurbs<2,real_t> shape(kv0,kv1,coefs,weights);

  gsMultiPatch<T> mp;
  mp.addPatch(shape);
  mp.addAutoBoundaries();

  // Elevate up to order p
  if (p>2)
  {
    for(index_t i = 2; i< p; ++i)
        mp.patch(0).degreeElevate();    // Elevate the degree
  }

  // // Refine n times
  // for(index_t i = 0; i< n; ++i)
  //     mp.patch(0).uniformRefine();

  return mp;
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
void writeStepOutput(const gsALMBase<T> * arcLength, const gsMultiPatch<T> & deformation, const std::string name, const gsMatrix<T> & points, const index_t extreme, const index_t kmax) // extreme: the column of point indices to compute the extreme over (default -1)
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
          << arcLength->solutionU().norm() << ",";
          for (index_t p=0; p!=points.cols(); p++)
          {
            file<< out(0,p) << ","
                << out(1,p) << ","
                << out(2,p) << ",";
          }

    file  << arcLength->solutionL() << ","
          << arcLength->indicator() << ","
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
          << arcLength->solutionU().norm() << ",";
          for (index_t p=0; p!=points.cols(); p++)
          {
            file<< out(0,p) << ","
                << out(1,p) << ","
                << std::max(abs(out2.col(p).maxCoeff()),abs(out2.col(p).minCoeff())) << ",";
          }

    file  << arcLength->solutionL() << ","
          << arcLength->indicator() << ","
          << "\n";
  }
  else
    GISMO_ERROR("Extremes setting unknown");

  file.close();
}
