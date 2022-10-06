/** @file gsThinShell_WrinklingPerturbed.cpp

    @brief Performs wrinkling simulations of different cases USING A PERTURBATION from a multipatch

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
#include <gsStructuralAnalysis/gsALMLoadControl.h>
#include <gsStructuralAnalysis/gsALMRiks.h>
#include <gsStructuralAnalysis/gsALMCrisfield.h>
#include <gsStructuralAnalysis/gsStaticNewton.h>

using namespace gismo;

template <class T>
gsMultiPatch<T> RectangularDomain(int n, int m, int p, int q, T L, T B, bool clamped = false, T offset = 0.1);
template <class T>
gsMultiPatch<T> RectangularDomain(int n, int p, T L, T B, bool clamped = false, T offset = 0.1);

template <class T>
gsMultiPatch<T> Rectangle(T L, T B);

template <class T>
gsMultiPatch<T> AnnularDomain(int n, int p, T R1, T R2);
template <class T>
gsMultiPatch<T> FrustrumDomain(int n, int p, T R1, T R2, T h);

template <class T>
void addClamping(gsMultiPatch<T> &mp, index_t patch, std::vector<boxSide> sides, T offset);

void writeToCSVfile(std::string name, gsMatrix<> matrix)
{
    std::ofstream file(name.c_str());
    for(int  i = 0; i < matrix.rows(); i++){
        for(int j = 0; j < matrix.cols(); j++){
           std::string str = std::to_string(matrix(i,j));
           if(j+1 == matrix.cols()){
               file<<str;
           }else{
               file<<str<<',';
           }
        }
        file<<'\n';
    }
  }

template <class T>
void initStepOutput( const std::string name, const gsMatrix<T> & points);

template <class T>
void writeStepOutput(const gsMatrix<T> & Uvec, const T lambda, const T indicator, const gsMultiPatch<T> & deformation, const std::string name, const gsMatrix<T> & points, const index_t extreme=-1, const index_t kmax=100);

void initSectionOutput( const std::string dirname, bool undeformed=false);

template <class T>
void writeSectionOutput(const gsMultiPatch<T> & mp, const std::string dirname, const index_t coordinate=0, const T coordVal=0.0, const index_t N=100, bool undeformed=false);

template<class T>
void writeBifurcation(std::string exportdir, gsMultiPatch<T> & U, gsMultiPatch<T> & DU, gsMultiPatch<T> & V, T & L, T & DL);

template<class T>
void readBifurcation(std::string importdir, gsMultiPatch<T> & U, gsMultiPatch<T> & DU, gsMultiPatch<T> & V, T & L, T & DL);

int main (int argc, char** argv)
{
    // Input options
    int numElevate  = 1;
    int numHref     = 1;
    int numElevateL = -1;
    int numHrefL    = -1;
    bool plot       = false;
    bool plotfiles  = false;
    bool stress       = false;
    bool membrane       = false;
    bool mesh  = false;
    bool quasiNewton = false;
    int quasiNewtonInt = -1;
    bool adaptive = false;
    int step = 10;
    int method = 2; // (0: Load control; 1: Riks' method; 2: Crisfield's method; 3: consistent crisfield method; 4: extended iterations)
    bool deformed = false;
    real_t perturbation = 0;

    real_t thickness = 25e-3;
    real_t E_modulus     = 1;
    real_t PoissonRatio = 0;
    real_t Density = 1e0;
    gsMultiPatch<> mp, mpBspline;
    real_t tau = 1e2;

    index_t Compressibility = 0;
    index_t material = 0;
    real_t Ratio = 7.0;
    bool composite = false;
    index_t impl = 1; // 1= analytical, 2= generalized, 3= spectral

    real_t aDim = 2.5;
    real_t bDim = 1.0;

    real_t relax = 1.0;

    int testCase = 1;

    int result = 0;

    bool write = false;
    bool writeG = false;
    bool writeP = false;


    bool crosssection = false;

    bool THB = false;


    index_t maxit = 20;

    // Arc length method options
    real_t dLini= 0; // General arc length
    real_t dL = 1e-2; // Ard length to find bifurcation
    real_t tol = 1e-6;
    real_t tolU = 1e-2;
    real_t tolF = 1e-2;

    std::string wn("data.csv");

    std::string importdir;

    std::string assemberOptionsFile("options/solver_options.xml");

    std::string fn;

    gsCmdLine cmd("Wrinkling analysis with thin shells using XML perturbation.");
    cmd.addString( "f", "file", "Input XML file for assembler options", assemberOptionsFile );
    cmd.addString( "d", "dir", "Input dir for start data", importdir );

    cmd.addInt("t", "testcase", "Test case: 0: clamped-clamped, 1: pinned-pinned, 2: clamped-free", testCase);

    cmd.addInt("r","hRefine", "Number of dyadic h-refinement (bisection) steps to perform before solving", numHref);
    cmd.addInt("e","degreeElevation", "Number of degree elevation steps to perform on the Geometry's basis before solving", numElevate);
    cmd.addInt("R","hRefine2", "Number of dyadic h-refinement (bisection) steps to perform before solving (secondary direction)", numHrefL);
    cmd.addInt("E","degreeElevation2", "Number of degree elevation steps to perform on the Geometry's basis before solving (secondary direction)", numElevateL);

    cmd.addInt( "M", "Material", "Material law",  material );
    cmd.addInt( "c", "Compressibility", "1: compressible, 0: incompressible",  Compressibility );
    cmd.addInt( "I", "Implementation", "Implementation: 1= analytical, 2= generalized, 3= spectral",  impl );
    cmd.addSwitch("composite", "Composite material", composite);

    cmd.addReal("T","hdim", "thickness of the plate", thickness);
    cmd.addReal("a","adim", "dimension a", aDim);
    cmd.addReal("b","bdim", "dimension b", bDim);

    cmd.addInt("m","Method", "Arc length method; 1: Crisfield's method; 2: RIks' method.", method);
    cmd.addReal("L","dL", "arc length", dL);
    cmd.addReal("l","dLini", "arc length at first step", dLini);
    cmd.addReal("A","relaxation", "Relaxation factor for arc length method", relax);

    cmd.addReal("P","perturbation", "perturbation factor", perturbation);

    cmd.addReal("F","factor", "factor for bifurcation perturbation", tau);
    cmd.addInt("q","QuasiNewtonInt","Use the Quasi Newton method every INT iterations",quasiNewtonInt);
    cmd.addInt("N", "maxsteps", "Maximum number of steps", step);

    cmd.addReal("U","tolU","displacement tolerance",tolU);

    cmd.addSwitch("adaptive", "Adaptive length ", adaptive);
    cmd.addSwitch("quasi", "Use the Quasi Newton method", quasiNewton);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addSwitch("plotfiles", "Write files for prostprocessing", plotfiles);
    cmd.addSwitch("mesh", "Plot mesh?", mesh);
    cmd.addSwitch("stress", "Plot stress in ParaView format", stress);
    cmd.addSwitch("write", "Write output to file", write);
    cmd.addSwitch("writeP", "Write perturbation", writeP);
    cmd.addSwitch("writeG", "Write refined geometry", writeG);
    cmd.addSwitch("cross", "Write cross-section to file", crosssection);
    cmd.addSwitch("membrane", "Use membrane model (no bending)", membrane);
    cmd.addSwitch("deformed", "plot on deformed shape", deformed);

    cmd.addSwitch("THB", "Use refinement", THB);

    cmd.addString("i","input", "Perturbation filename", fn);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo<<"-> Initializing parameters...";
    gsStopwatch time;

    gsFileData<> fd(assemberOptionsFile);
    gsOptionList opts;
    fd.getFirst<gsOptionList>(opts);

    if (dLini==0)
    {
      dLini= dL;
    }

    if (numHrefL==-1)
      numHrefL = numHref;
    if (numElevateL==-1)
      numElevateL = numElevate;

    if ((!Compressibility) && (material!=0))
      PoissonRatio = 0.5;
    else
      PoissonRatio = 0.499;

    // ![Material data]

    E_modulus = 3500;
    PoissonRatio = 0.31;
    gsDebug<<"E = "<<E_modulus<<"; nu = "<<PoissonRatio<<"\n";

    aDim = 380;
    bDim = 128;
    // ![Material data]

    // ![Read Geometry files]
    std::vector<boxSide> sides;
    if (testCase==2)
    {
    	sides.push_back(boundary::west);
    	sides.push_back(boundary::east);
    }

    mpBspline = Rectangle(aDim,    bDim   );

    for(index_t i = 0; i< numElevate; ++i)
      mpBspline.patch(0).degreeElevate();    // Elevate the degree

    // h-refine
    for(index_t i = 0; i< numHref; ++i)
      mpBspline.patch(0).uniformRefine();

    addClamping(mpBspline,0,sides, 1e-2);

    index_t N = mpBspline.patch(0).coefs().rows();
    mpBspline.patch(0).coefs().col(2) = gsMatrix<>::Random(N,1);


    mpBspline.patch(0).coefs().col(2) *= perturbation;

    // Cast all patches of the mp object to THB splines
    gsTHBSpline<2,real_t> thb;
    if (THB)
    {
      for (size_t k=0; k!=mpBspline.nPatches(); ++k)
      {
          gsTensorBSpline<2,real_t> *geo = dynamic_cast< gsTensorBSpline<2,real_t> * > (&mpBspline.patch(k));
          thb = gsTHBSpline<2,real_t>(*geo);
          mp.addPatch(thb);
      }

      gsMatrix<> refBoxes(2,2);
      if      (testCase==2 || testCase==3)
      {
        refBoxes.col(0) << 0.25,0.25;
        refBoxes.col(1) << 0.75,0.75;
      }
      else if (testCase==4 || testCase==5)
      {
        refBoxes.col(0) << 0.25,0.00;
        refBoxes.col(1) << 0.75,0.25;
      }
      else if (testCase==6 || testCase==7)
      {
        refBoxes.col(0) << 0.00,0.00;
        refBoxes.col(1) << 0.25,0.25;
      }

      int refExtension = 1;
      std::vector<index_t> elements = mp.patch(0).basis().asElements(refBoxes, refExtension);
      mp.patch(0).refineElements( elements );
    }
    else
      mp = mpBspline;

    gsInfo<<"finished ["<<time.stop()<<" s]\n";

    gsMultiBasis<> dbasis(mp);
    gsInfo<<"Basis (patch 0): "<< mp.patch(0).basis() << "\n";

    gsInfo<<"-> Initializing boundary conditions...";
    time.restart();

    // Boundary conditions
    gsBoundaryConditions<> BCs,BCs_ini;
    BCs.setGeoMap(mp);

    gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();

    // Initiate Surface forces
    std::string tx("0");
    std::string ty("0");
    std::string tz("0");

    gsVector<> tmp(3);
    gsVector<> neu(3);
    tmp << 0, 0, 0;

    gsConstantFunction<> displ(0.05,3);

    // Buckling coefficient
    real_t Load = 0;

    std::string output = "solution";
    std::string dirname = "ArcLengthResults";

    gsMatrix<> writePoints(2,2);
    writePoints.col(0)<< 1.0,0.5;
    writePoints.col(1)<< 1.0,1.0;
    index_t cross_coordinate = 0;
    real_t cross_val = 0.5;

    if (testCase == 1)
    {
      displ.setValue(0.05,3);
      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0 ,false,0);
      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0 ,false,1);
      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0 ,false,2);

      BCs.addCondition(boundary::north, condition_type::dirichlet, &displ, 0 ,false,1);
      BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0 ,false,2);

      BCs_ini = BCs;
      BCs_ini.addCondition(boundary::north, condition_type::dirichlet, 0, 0 ,false,0);

      BCs.addCondition(boundary::north, condition_type::collapsed, 0, 0 ,false,0);

      real_t Load = 1.0;
      gsVector<> point(2); point<< 1.0, 1.0 ;
      gsVector<> load (3); load << Load,0.0, 0.0;
      pLoads.addLoad(point, load, 0 );

      std::stringstream ss;
      ss<<perturbation;
      dirname = dirname + "/ShearSheet_Perturbed=" + ss.str() + "_r=" + std::to_string(numHref) + "_e=" + std::to_string(numElevate) + "_M=" + std::to_string(material) + "_c=" + std::to_string(Compressibility) + "_t=" + std::to_string(thickness);
      output =  "solution";
      wn = output + "data";
    }
    else if (testCase == 2)
    {
      displ.setValue(0.05,3);
      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0 ,false,0);
      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0 ,false,1);
      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0 ,false,2);

      BCs.addCondition(boundary::north, condition_type::dirichlet, &displ, 0 ,false,1);
      BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0 ,false,2);

      BCs.addCondition(boundary::east, condition_type::clamped, 0, 0 ,false,2);
      BCs.addCondition(boundary::west, condition_type::clamped, 0, 0 ,false,2);

      BCs_ini = BCs;
      BCs_ini.addCondition(boundary::north, condition_type::dirichlet, 0, 0 ,false,0);

      BCs.addCondition(boundary::north, condition_type::collapsed, 0, 0 ,false,0);

      real_t Load = 1.0;
      gsVector<> point(2); point<< 1.0, 1.0 ;
      gsVector<> load (3); load << Load,0.0, 0.0;
      pLoads.addLoad(point, load, 0 );

      std::stringstream ss;
      ss<<perturbation;
      dirname = dirname + "/ShearSheetRestrained_Perturbed=" + ss.str() + "_r=" + std::to_string(numHref) + "_e=" + std::to_string(numElevate) + "_M=" + std::to_string(material) + "_c=" + std::to_string(Compressibility);
      output =  "solution";
      wn = output + "data";
    }

    if (THB)
      dirname = dirname + "_THB";

    if (plot || writeG || writeP || write || plotfiles)
    {
      std::string commands = "mkdir -p " + dirname;
      const char *command = commands.c_str();
      system(command);
    }

    // plot geometry
    if (plot)
      gsWriteParaview(mp,dirname + "/" + "mp",1000,true);

    if (writeG || plotfiles)
    {
      gsWrite(mp,dirname + "/" + "geometry");
      gsInfo<<"Geometry written in: " + dirname + "/" + "geometry.xml\n";
    }

    gsInfo<<"finished ["<<time.stop()<<" s]\n";

    gsInfo<<"-> Initializing shell...";
    time.restart();

    gsFunctionExpr<> surfForce(tx,ty,tz,3);
    // Initialise solution object
    gsMultiPatch<> mp_def = mp;
    gsSparseSolver<>::LU solver;

    // Linear isotropic material model
    gsConstantFunction<> force(tmp,3);
    gsFunctionExpr<> t(std::to_string(thickness), 3);
    gsFunctionExpr<> E(std::to_string(E_modulus),3);
    gsFunctionExpr<> nu(std::to_string(PoissonRatio),3);
    gsFunctionExpr<> rho(std::to_string(Density),3);
    gsConstantFunction<> ratio(Ratio,3);

    real_t mu = E_modulus / (2 * (1 + PoissonRatio));
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
            materialMatrix = new gsMaterialMatrixComposite<3,real_t>(mp,Ts,Gs,Phis);
        }
        else
        {
            parameters.resize(2);
            options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
            options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
            materialMatrix = getMaterialMatrix<3,real_t>(mp,t,parameters,rho,options);
        }
    }
    else
    {
        options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",material);
        options.addSwitch("Compressibility","Compressibility: (false): Imcompressible | (true): Compressible",Compressibility);
        options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",impl);
        materialMatrix = getMaterialMatrix<3,real_t>(mp,t,parameters,rho,options);
    }

    gsThinShellAssemblerBase<real_t>* assembler;
    if(membrane)
        assembler = new gsThinShellAssembler<3, real_t, false>(mp,dbasis,BCs,force,materialMatrix);
    else
        assembler = new gsThinShellAssembler<3, real_t, true >(mp,dbasis,BCs,force,materialMatrix);


    // Construct assembler object
    assembler->setOptions(opts);
    assembler->setPointLoads(pLoads);

    typedef std::function<gsSparseMatrix<real_t> (gsVector<real_t> const &)>                                Jacobian_t;
    typedef std::function<gsVector<real_t> (gsVector<real_t> const &, real_t, gsVector<real_t> const &) >   ALResidual_t;
    typedef std::function<gsVector<real_t> (gsVector<real_t> const &) >                                     Residual_t;
    // Function for the Jacobian
    Jacobian_t Jacobian = [&assembler](gsVector<real_t> const &x)
    {
      gsMultiPatch<> mp_def;
      assembler->homogenizeDirichlet();
      assembler->constructSolution(x,mp_def);
      assembler->assembleMatrix(mp_def);
      gsSparseMatrix<real_t> m = assembler->matrix();
      // gsInfo<<"matrix = \n"<<m.toDense()<<"\n";
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

    // Function for the Residual
    Residual_t Residual = [&assembler](gsVector<real_t> const &x)
    {
      gsMultiPatch<> mp_def;
      assembler->homogenizeDirichlet();
      assembler->constructSolution(x,mp_def);
      assembler->assembleVector(mp_def);
      return assembler->rhs(); // - lam * force;
    };

    gsInfo<<"finished ["<<time.stop()<<" s]\n";

    gsInfo<<"-> Initializing assembly...";
    time.restart();


    if (importdir.empty())
    {
      // Assemble vector and matrix for first step static solve
      assembler->updateBCs(BCs_ini);
      assembler->assemble();
      gsVector<> vector = assembler->rhs();
      gsSparseMatrix<> matrix = assembler->matrix();

      // Run first static solve
      gsStaticNewton<real_t> staticSolver(matrix,vector,Jacobian,Residual);
      gsVector<> Uini = staticSolver.solveNonlinear();
      assembler->constructDisplacement(Uini,mp_def);
    }

    // Assemble vector and matrix for ALM
    assembler->updateBCs(BCs);
    // Force vector excl Dirichlet contributions
    assembler->homogenizeDirichlet();
    assembler->assemble();
    gsVector<> Force = assembler->rhs();

    gsInfo<<"finished ["<<time.stop()<<" s]\n";

    gsInfo<<"-> Initializing ALM...";
    time.restart();

    gsALMBase<real_t> * arcLength;
    if (method==0)
      arcLength = new gsALMLoadControl<real_t>(Jacobian, ALResidual, Force);
    else if (method==1)
      arcLength = new gsALMRiks<real_t>(Jacobian, ALResidual, Force);
    else if (method==2)
      arcLength = new gsALMCrisfield<real_t>(Jacobian, ALResidual, Force);
    else
      GISMO_ERROR("Method "<<method<<" unknown");

    // if (!membrane)
    // {
      arcLength->options().setInt("Solver",0); // LDLT solver
      arcLength->options().setInt("BifurcationMethod",0); // 0: determinant, 1: eigenvalue
    // }
    // else
    // {
      // arcLength->options().setInt("Solver",1); // CG solver
      // arcLength->options().setInt("BifurcationMethod",1); // 0: determinant, 1: eigenvalue
    // }

    arcLength->options().setReal("Length",dLini);
    if (method==2)
    {
      arcLength->options().setInt("AngleMethod",2); // 0: step, 1: iteration, 2: predictor
      arcLength->options().setReal("Scaling",0);
    }
    arcLength->options().setSwitch("AdaptiveLength",adaptive);
    arcLength->options().setInt("AdaptiveIterations",5);
    arcLength->options().setReal("Perturbation",tau);
    arcLength->options().setReal("Tol",tol);
    arcLength->options().setReal("TolU",tolU);
    arcLength->options().setReal("TolF",tolF);
    arcLength->options().setInt("MaxIter",maxit);
    arcLength->options().setSwitch("Verbose",true);
    arcLength->options().setReal("Relaxation",relax);
    arcLength->options().setInt("SingularPointFailure",0);
    if (quasiNewtonInt>0)
    {
      quasiNewton = true;
      arcLength->options().setInt("QuasiIterations",quasiNewtonInt);
    }
    arcLength->options().setSwitch("Quasi",quasiNewton);

    gsInfo<<arcLength->options();
    arcLength->applyOptions();
    arcLength->initialize(false); //  don't compute stability

    gsMultiPatch<> deformation = mp;

    // Make objects for previous solutions
    real_t Lold;
    gsMatrix<> Uold, Vold, Uinit;
    real_t indicator = 0.0;
    real_t indicatorOld = 0.0;
    index_t stability = 0.0;
    index_t stabilityOld = 0.0;

    gsInfo<<"finished ["<<time.stop()<<" s]\n";

    gsInfo<<"-> Setting start point...";
    time.restart();
    if (importdir.empty())
    {
      gsVector<> Uini = assembler->constructSolutionVector(mp_def);
      Lold = 0;
      Uold = Uini;
      Vold = gsMatrix<>::Zero(Uold.rows(),1);
      arcLength->setSolution(Uold,Lold);

      importdir =  dirname + "/" + "start" + "_U=" + std::to_string(Uold.norm()) + "_L=" + std::to_string(Lold);
      std::string commands = "mkdir -p " + importdir;
      const char *command = commands.c_str();
      system(command);
    }
    else
    {
      real_t DL;
      gsMatrix<> DUvec;

      gsMultiPatch<> U, DU, V;

      readBifurcation(importdir,U,DU,V,Lold,DL);

      gsDebugVar(Lold);
      Uold = assembler->constructSolutionVector(U);
      DUvec= assembler->constructSolutionVector(DU);
      Vold = assembler->constructSolutionVector(V);

      arcLength->setSolutionStep(1. / tau * Vold,0);
    }

    gsInfo<<"Starting branch from |U| = "<<(Uold).norm()
      <<"; lambda = "<<Lold<<"\n";
    Uinit = Uold + 1. / tau * Vold; // add perturbation
    arcLength->setSolution(Uinit,Lold);

    arcLength->computeStability(Uinit,true,-1e-1);
    indicatorOld = 0;
    indicator = arcLength->indicator();
    stabilityOld = 0;
    stability = arcLength->isStable();

    gsMatrix<> solVector;
    bool bisected = false;

    bool quit = false;

    gsParaviewCollection collection(importdir + "/" + output);
    if (write)
      initStepOutput(importdir + "/" + wn, writePoints);

    assembler->constructSolution(Uinit,mp_def);

    deformation = mp_def;
    deformation.patch(0).coefs() -= mp.patch(0).coefs();// assuming 1 patch here

    if (plot)
    {
      gsField<> solField;
      if (deformed)
        solField= gsField<>(mp_def,deformation);
      else
        solField= gsField<>(mp,deformation);

      std::string fileName = importdir + "/" + output + util::to_string(0);
      gsWriteParaview<>(solField, fileName, 5000,mesh);
      fileName = output + util::to_string(0) + "0";
      collection.addTimestep(fileName,0,".vts");
      if (mesh) collection.addTimestep(fileName,0,"_mesh.vtp");
    }

    if (write)
      writeStepOutput(Uinit, Lold, 0.0,deformation, importdir + "/" + wn, writePoints,1, 201);

    quit = false;

    real_t dL0 = dLini;

    gsInfo<<"finished ["<<time.stop()<<" s]\n";

    for (index_t k=1; k<step; k++)
    {
      gsInfo<<"Load step "<< k<<"\n";
      arcLength->step();

      if (!(arcLength->converged()))
      {
        gsInfo<<"Error: Loop terminated, arc length method did not converge.\n";
        dL0 = dL0 / 2.;
        GISMO_ASSERT(dL0 / dL > 1e-6,"Step size is becoming very small...");
        arcLength->setLength(dL0);
        arcLength->setSolution(Uold,Lold);
        bisected = true;
        k-=1;
        continue;
      }
      if (perturbation==0)
      {
        arcLength->computeStability(arcLength->solutionU(),true);

        stability = arcLength->stability();
        gsInfo<<"stability = "<<stability<<"; stability old = "<<stabilityOld<<"\n";
        if (stability != stabilityOld && stabilityOld!=0)
        {
          gsInfo<<"Bifurcation spotted!"<<"\n";
          arcLength->computeSingularPoint(1e-4, 25, Uold, Lold, 1e-7, 1e-2, false, true);

          Uold = arcLength->solutionU();
          Lold = arcLength->solutionL();

          //////////////////////////////////////////////////////////////////////////////////////////////////
          ///////////////////////////////////////// Export to file /////////////////////////////////////////
          //////////////////////////////////////////////////////////////////////////////////////////////////

          gsInfo<<"Bifurcation point (1) or limit point (0)? : "<<arcLength->testSingularPoint(arcLength->solutionU(), arcLength->solutionL(),1e-4,25,true )<<"\n";

          gsMultiPatch<> U, DU, V;
          gsMultiPatch<> mp_perturbation;
          gsField<> solField;
          real_t  L = arcLength->solutionL(),
                  DL = arcLength->solutionDL();

          assembler->constructDisplacement(arcLength->solutionU(),U);
          assembler->constructDisplacement(arcLength->solutionDU(),DU);
          assembler->constructDisplacement(arcLength->solutionV(),V);

          gsMatrix<> modes = arcLength->computeModes(arcLength->solutionU(),true);
          gsVector<> Vvec = gsVector<>::Zero(modes.rows());

          gsInfo<<"Found "<<modes.cols()<<" modes\n";

          std::string exportdir1 = dirname + "/" + "start" + "_U=" + std::to_string(arcLength->solutionU().norm()) + "_L=" + std::to_string(arcLength->solutionL()) + "_pos";
          writeBifurcation(exportdir1,U,DU,V,L,DL);
          mp_perturbation = V;
          mp_perturbation.patch(0).coefs() -= mp.patch(0).coefs();// assuming 1 patch here
          mp_perturbation.patch(0).coefs() *= 1./ mp_perturbation.patch(0).coefs().col(2).maxCoeff();
          solField= gsField<>(mp,mp_perturbation);
          gsWriteParaview(solField,exportdir1 + "/" +"bifurcation",10000,mesh);

          V.patch(0).coefs().col(2) = -V.patch(0).coefs().col(2);
          std::string exportdir2 = dirname + "/" + "start" + "_U=" + std::to_string(arcLength->solutionU().norm()) + "_L=" + std::to_string(arcLength->solutionL()) + "_neg";
          writeBifurcation(exportdir2,U,DU,V,L,DL);
          mp_perturbation = V;
          mp_perturbation.patch(0).coefs() -= mp.patch(0).coefs();// assuming 1 patch here
          mp_perturbation.patch(0).coefs() *= 1./ mp_perturbation.patch(0).coefs().col(2).maxCoeff();
          solField= gsField<>(mp,mp_perturbation);
          gsWriteParaview(solField,exportdir2 + "/" +"bifurcation",10000,mesh);

          if (modes.cols() > 1)
          {
            for (index_t k = 0; k!=modes.cols(); ++k)
              Vvec += modes.col(k);

            assembler->constructDisplacement(Vvec,V);

            for (index_t k = 0; k!=modes.cols(); ++k)
            {
              assembler->constructDisplacement(modes.col(k),V);
              mp_perturbation = V;
              mp_perturbation.patch(0).coefs() -= mp.patch(0).coefs();// assuming 1 patch here
              mp_perturbation.patch(0).coefs() *= 1./ mp_perturbation.patch(0).coefs().col(2).maxCoeff();
              solField= gsField<>(mp,mp_perturbation);
              gsWriteParaview(solField,exportdir1 + "/" +"bifurcation_mode" + std::to_string(k) + "_",10000,mesh);
            }



            for (index_t k = 0; k!=modes.cols(); ++k)
            {
              assembler->constructDisplacement(modes.col(k),V);
              mp_perturbation = V;
              mp_perturbation.patch(0).coefs() -= mp.patch(0).coefs();// assuming 1 patch here
              mp_perturbation.patch(0).coefs() *= 1./ mp_perturbation.patch(0).coefs().col(2).maxCoeff();
              solField= gsField<>(mp,mp_perturbation);
              gsWriteParaview(solField,exportdir2 + "/" +"bifurcation_mode" + std::to_string(k) + "_",10000,mesh);
            }
          }

          /////////////////////////////////////////////////////////////////////////////////////////////////

          gsMatrix<> Utmp = arcLength->solutionU();
          real_t Ltmp = arcLength->solutionL();
          assembler->constructSolution(Utmp,mp_def);

          deformation = mp_def;
          deformation.patch(0).coefs() -= mp.patch(0).coefs();// assuming 1 patch here

          if (plot)
          {
            gsField<> solField;
            if (deformed)
              solField= gsField<>(mp_def,deformation);
            else
              solField= gsField<>(mp,deformation);

            std::string fileName = importdir + "/" + output + util::to_string(k);
            gsWriteParaview<>(solField, fileName, 5000,mesh);
            fileName = output + util::to_string(k) + "0";
            collection.addTimestep(fileName,k,".vts");
            if (mesh) collection.addTimestep(fileName,k,"_mesh.vtp");
          }

          if (write)
            writeStepOutput(Utmp,Ltmp,static_cast<real_t>(indicator),deformation, importdir + "/" + wn, writePoints,1, 201);

          quit = false;
          stabilityOld = stability;
          stability = 1;

          arcLength->setSolution(Uold,Lold);
        }
        else
          stabilityOld = stability;
      }

      indicatorOld = indicator;
      indicator = arcLength->indicator();
      solVector = arcLength->solutionU();

      Uold = solVector;
      Lold = arcLength->solutionL();

      // Remove perturbation
      // if (k==2)
      // {
      //   solVector = - 1. / tau * Vold;
      //   Uold = solVector;
      //   arcLength->setSolution(Uold,Lold);
      // }

      assembler->constructSolution(solVector,mp_def);

      deformation = mp_def;
      deformation.patch(0).coefs() -= mp.patch(0).coefs();// assuming 1 patch here

      if (plot)
      {
        gsField<> solField;
        if (deformed)
          solField= gsField<>(mp_def,deformation);
        else
          solField= gsField<>(mp,deformation);

        std::string fileName = importdir + "/" + output + util::to_string(k);
        gsWriteParaview<>(solField, fileName, 5000,mesh);
        fileName = output + util::to_string(k) + "0";
        collection.addTimestep(fileName,k,".vts");
        if (mesh) collection.addTimestep(fileName,k,"_mesh.vtp");
      }

      if (write)
        writeStepOutput(Uold,Lold,static_cast<real_t>(indicator),deformation, importdir + "/" + wn, writePoints,1, 201);

      if (!bisected)
      {
        dL0 = dL;
        arcLength->setLength(dL0);
      }
      bisected = false;

      if (quit)
        break;

    }

    if (plot)
    {
      collection.save();
    }

  return result;
}

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
  for (index_t k = 0; k < len1; k++)
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
         geo->insertKnot(1 - std::min(offset, dknot1 / 2),1);
       else if (*it==boundary::south) // south
         geo->insertKnot(std::min(offset, dknot1 / 2),1);
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


template <class T>
gsMultiPatch<T> AnnularDomain(int n, int p, T R1, T R2)
{
  // -------------------------------------------------------------------------
  // --------------------------Make beam geometry-----------------------------
  // -------------------------------------------------------------------------
  int dim = 3; //physical dimension
  gsKnotVector<> kv0;
  kv0.initUniform(0,1,0,3,1);
  gsKnotVector<> kv1;
  kv1.initUniform(0,1,0,3,1);

  // Make basis
  // gsTensorNurbsBasis<2,T> basis(kv0,kv1);

  // Initiate coefficient matrix
  gsMatrix<> coefs(9,dim);

  coefs<<R1,0,0,
  (R1+R2)/2,0,0,
  R2,0,0,
  R1,R1,0,
  (R1+R2)/2,(R1+R2)/2,0,
  R2,R2,0,
  0,R1,0,
  0,(R1+R2)/2,0,
  0,R2,0;

  gsMatrix<> weights(9,1);
  weights<<1,1,1,
  0.707106781186548,0.707106781186548,0.707106781186548,
  1,1,1;

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

  // Refine n times
  for(index_t i = 0; i< n; ++i)
      mp.patch(0).uniformRefine();

  return mp;
}

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

  gsDebug<<kv1;

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
        for (index_t k=0; k < points.cols(); k++)
        {
          file<< "point "<<k<<" - x" << ","
              << "point "<<k<<" - y" << ","
              << "point "<<k<<" - z" << ",";
        }

  file  << "Lambda" << ","
        << "Indicator" << ","
        << "Branch" << "\n";
  file.close();

  gsInfo<<"Step results will be written in file: "<<name<<"\n";
}

template <class T>
void writeStepOutput(const gsMatrix<T> & Uvec, const T lambda, const T indicator, const gsMultiPatch<T> & deformation, const std::string name, const gsMatrix<T> & points,  const index_t extreme, const index_t kmax) // extreme: the column of point indices to compute the extreme over (default -1)
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
          << Uvec.norm() << ",";
          for (index_t p=0; p!=points.cols(); p++)
          {
            file<< out(0,p) << ","
                << out(1,p) << ","
                << out(2,p) << ",";
          }

    file  << lambda << ","
          << indicator << "\n";
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
          << Uvec.norm() << ",";
          for (index_t p=0; p!=points.cols(); p++)
          {
            T min = out2.col(p).minCoeff();
            T max = out2.col(p).maxCoeff();
            T ext = std::abs(min) > std::abs(max) ? std::abs(min) : max; // return absolute value
            file<< out(0,p) << ","
                << out(1,p) << ","
                << ext      << ",";
          }

    file  << lambda << ","
          << indicator << "\n";
  }
  else
    GISMO_ERROR("Extremes setting unknown");

  file.close();
}

void initSectionOutput(const std::string dirname, bool undeformed)
{
  std::ofstream file2, file3, file4;
  std::string wn2,wn3,wn4;

  if (! undeformed)
  {
    wn2 = dirname + "/" + "pointdataX.txt";
    wn3 = dirname + "/" + "pointdataY.txt";
    wn4 = dirname + "/" + "pointdataZ.txt";
  }
  else
  {
    wn2 = dirname + "/" + "pointdataX0.txt";
    wn3 = dirname + "/" + "pointdataY0.txt";
    wn4 = dirname + "/" + "pointdataZ0.txt";
  }

  file2.open(wn2,std::ofstream::out);
  file2.close();

  file3.open(wn3,std::ofstream::out);
  file3.close();

  file4.open(wn4,std::ofstream::out);
  file4.close();

  gsInfo<<"Cross-section results will be written in directory: "<<dirname<<"\n";
}

template <class T>
void writeSectionOutput(const gsMultiPatch<T> & mp, const std::string dirname, const index_t coordinate, const T coordVal, const index_t N, bool undeformed) // coordinate: the column which remains constant at coordVal
{
  gsMatrix<T> P(2,1);
  gsMatrix<T> tmp;
  P.setZero();
  P.at(coordinate) = coordVal;

  std::ofstream file2, file3, file4;
  std::string wn2,wn3,wn4;

  if (! undeformed)
  {
    wn2 = dirname + "/" + "pointdataX.txt";
    wn3 = dirname + "/" + "pointdataY.txt";
    wn4 = dirname + "/" + "pointdataZ.txt";
  }
  else
  {
    wn2 = dirname + "/" + "pointdataX0.txt";
    wn3 = dirname + "/" + "pointdataY0.txt";
    wn4 = dirname + "/" + "pointdataZ0.txt";
  }

  file2.open(wn2,std::ofstream::out | std::ofstream::app);
  file3.open(wn3,std::ofstream::out | std::ofstream::app);
  file4.open(wn4,std::ofstream::out | std::ofstream::app);


  gsMatrix<T> out(3,N); // evaluation points in the rows, output (per coordinate) in columns
    for (int k = 0; k != N; k ++)
    {
      P.at(1-coordinate) = 1.0*k/(N-1);

      mp.patch(0).eval_into(P,tmp);
      out.col(k) = tmp; // z coordinate

      std::string str2 = std::to_string(out(0,k));
      std::string str3 = std::to_string(out(1,k));
      std::string str4 = std::to_string(out(2,k));
      if(k+1 == N)
      {
          file2<<str2;
          file3<<str3;
          file4<<str4;
      }
      else{
          file2<<str2<<',';
          file3<<str3<<',';
          file4<<str4<<',';
      }
    }
    file2<<'\n';
    file2.close();
    file3<<'\n';
    file3.close();
    file4<<'\n';
    file4.close();
}

template<class T>
void writeBifurcation(std::string exportdir, gsMultiPatch<T> & U, gsMultiPatch<T> & DU, gsMultiPatch<T> & V, T & L, T & DL)
{
  //////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////// Export to file /////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////

  std::string commands = "mkdir -p " + exportdir;
  const char *command = commands.c_str();
  system(command);

  gsWrite(U ,exportdir + "/" + "U" + ".xml");
  gsWrite(DU,exportdir + "/" + "DU"+ ".xml");
  gsWrite(V ,exportdir + "/" + "V" + ".xml");

  std::ofstream file;
  file.open(exportdir + "/" + "L" + ".txt",std::ofstream::out);
  file<<L;
  file.close();

  file.open(exportdir + "/" + "DL"+ ".txt",std::ofstream::out);
  file<<DL;
  file.close();

  /////////////////////////////////////////////////////////////////////////////////////////////////

}

template<class T>
void readBifurcation(std::string importdir, gsMultiPatch<T> & U, gsMultiPatch<T> & DU, gsMultiPatch<T> & V, T & L, T & DL)
{
  //////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////// Export to file /////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////

  gsReadFile<T>(importdir + "/" + "U" + ".xml",U );
  gsReadFile<T>(importdir + "/" + "DU"+ ".xml",DU);
  gsReadFile<T>(importdir + "/" + "V" + ".xml",V );

  std::ifstream file;
  std::string string;

  file.open(importdir + "/" + "L" + ".txt");
  file>>string;
  L = std::stod(string);

  gsDebugVar(L);

  file.open(importdir + "/" + "DL"+ ".txt");
  file>>string;
  DL= std::stod(string);
  gsDebugVar(DL);

  /////////////////////////////////////////////////////////////////////////////////////////////////

}