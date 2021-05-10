/** @file benchmark_MaterialTest.cpp

    @brief Stretches a planar material with both ends clamped

    See fig 1 of Roohbakhshan and Sauer 2017
    Roohbakhshan, F., & Sauer, R. A. (2017). Efficient isogeometric thin shell formulations for soft biological materials. Biomechanics and Modeling in Mechanobiology, 16(5), 1569–1597. https://doi.org/10.1007/s10237-017-0906-6

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#include <gismo.h>

#include <gsKLShell/gsThinShellAssembler.h>
#include <gsKLShell/getMaterialMatrix.h>

// #include <gsThinShell/gsArcLengthIterator.h>
#include <gsStructuralAnalysis/gsArcLengthIterator.h>

using namespace gismo;

template <class T>
gsMultiPatch<T> Rectangle(T L, T B);

template <class T>
void initStepOutput( const std::string name, const gsMatrix<T> & points);

template <class T>
void writeStepOutput(const gsArcLengthIterator<T> & arcLength, const gsMultiPatch<T> & deformation, const std::string name, const gsMatrix<T> & points);

int main (int argc, char** argv)
{
    // Input options
    int numElevate  = -1;
    int numHref     = -1;
    bool plot       = false;
    bool stress       = false;
    bool quasiNewton = false;
    int quasiNewtonInt = -1;
    bool adaptive = false;
    int step = 10;
    int method = -1; // (0: Load control; 1: Riks' method; 2: Crisfield's method; 3: consistent crisfield method; 4: extended iterations)
    bool deformed = false;

    index_t Compressibility = 0;
    index_t material = 0;
    bool composite = false;
    index_t impl = 1; // 1= analytical, 2= generalized, 3= spectral

    real_t relax = 1.0;

    int testCase = 0;

    int result = 0;

    bool write = false;

    index_t maxit = 20;

    // Arc length method options
    real_t dL = -1; // General arc length
    real_t tol = 1e-6;
    real_t tolU = 1e-6;
    real_t tolF = 1e-3;

    std::string wn("data.csv");

    gsCmdLine cmd("Wrinkling analysis with thin shells.");

    cmd.addInt("t", "testcase", "Test case: 0: clamped-clamped, 1: pinned-pinned, 2: clamped-free", testCase);

    cmd.addInt("r","hRefine", "Number of dyadic h-refinement (bisection) steps to perform before solving", numHref);
    cmd.addInt("e","degreeElevation", "Number of degree elevation steps to perform on the Geometry's basis before solving", numElevate);

    cmd.addInt( "M", "Material", "Material law",  material );
    cmd.addInt( "c", "Compressibility", "1: compressible, 0: incompressible",  Compressibility );
    cmd.addInt( "I", "Implementation", "Implementation: 1= analytical, 2= generalized, 3= spectral",  impl );
    cmd.addSwitch("composite", "Composite material", composite);

    cmd.addInt("m","Method", "Arc length method; 1: Crisfield's method; 2: RIks' method.", method);
    cmd.addReal("L","dLb", "arc length", dL);
    cmd.addReal("A","relaxation", "Relaxation factor for arc length method", relax);

    cmd.addInt("q","QuasiNewtonInt","Use the Quasi Newton method every INT iterations",quasiNewtonInt);
    cmd.addInt("N", "maxsteps", "Maximum number of steps", step);

    cmd.addSwitch("adaptive", "Adaptive length ", adaptive);
    cmd.addSwitch("quasi", "Use the Quasi Newton method", quasiNewton);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addSwitch("stress", "Plot stress in ParaView format", stress);
    cmd.addSwitch("write", "Write output to file", write);
    cmd.addSwitch("deformed", "plot on deformed shape", deformed);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    real_t aDim,bDim;
    real_t thickness;
    real_t E_modulus;
    real_t PoissonRatio;
    real_t Density = 1e0;
    real_t Ratio;
    real_t mu;

    real_t A1,A2,A3;
    real_t M1,M2,M3;

    real_t Load;

    if (testCase==0)
    {
      if (dL == -1) { dL = 1e-1; }
      if (numHref == -1) { numHref = 3; }
      if (numElevate == -1) { numElevate = 2; }
      if (method == -1) { method = 2; }

      bDim = 0.14; aDim = 2*bDim;
      thickness = 0.14e-3;
      if ((!Compressibility) && (material!=0))
        PoissonRatio = 0.5;
      else
        PoissonRatio = 0.499;

      real_t C01,C10;
      if (material==3)
      {
        C10 = 6.21485502e4; // c1/2
        C01 = 15.8114570e4; // c2/2
        Ratio = C10/C01;
        mu = 2*(C01+C10);
      }
      else //if (material==1 || material==4)
      {
        C10 = 19.1010178e4;
        mu = 2*C10;
        if (material==14)
        {
          A1 = 1.1;
          A2 = -7.0;
          A3 = -3.0;
          M1 = 1.0*mu;
          M2 = -0.003*mu;
          M3 = -0.4*mu;
        }
      }
      E_modulus = 2*mu*(1+PoissonRatio);
      Load = 0.25e0;
    }
    /*
      Case 1: Constrained tension (see Roohbakhshan2017)
    */
    else if (testCase==1)
    {
      if (dL == -1) { dL = 1e-2; }
      if (numHref == -1) { numHref = 4; }
      if (numElevate == -1) { numElevate = 2; }
      if (method == -1) { method = 2; }

      aDim = 9e-3;
      bDim = 3e-3;
      thickness = 0.3e-3;

      mu = 10e3;
      if ((!Compressibility) && (material!=0))
        PoissonRatio = 0.5;
      else
        PoissonRatio = 0.45;

      if (material==1)
      {
        mu = 10e3;
        if (material==4)
        {
          A1 = 1.3;
          A2 = 5.0;
          A3 = -2.0;
          M1 = 6.3e5/4.225e5*mu;
          M2 = 0.012e5/4.225e5*mu;
          M3 = -0.1e5/4.225e5*mu;
        }
      }
      else if (material==3)
      {
        mu = 30e3;
      }

      E_modulus = 2*mu*(1+PoissonRatio);

      Ratio = 0.5;
      Load = 1e-1;
    }
    else if (testCase==2)
    {
      if (dL == -1) { dL = 11.8421052632; }
      if (numHref == -1) { numHref = 0; }
      if (numElevate == -1) { numElevate = 2; }
      if (method == -1) { method = 0; }

      thickness = 0.15;
      bDim = thickness / 1.9e-3;
      aDim = 2*bDim;

      if ((!Compressibility) && (material!=0))
        PoissonRatio = 0.5;
      else
        PoissonRatio = 0.45;

      E_modulus = 1.0;
      mu = E_modulus / (2 * (1 + PoissonRatio));

      A1 = 1.3;
      A2 = 5.0;
      A3 = -2.0;
      M1 = 6.3e5/4.225e5*mu;
      M2 = 0.012e5/4.225e5*mu;
      M3 = -0.1e5/4.225e5*mu;

      Ratio = 0.5;
      Load = 1e-1;
    }
    else
      GISMO_ERROR("TESTCASE UNKNOWN");

    gsMultiPatch<> mp,mp_def;

    mp = Rectangle(aDim   , bDim   );

    for(index_t i = 0; i< numElevate; ++i)
      mp.patch(0).degreeElevate();    // Elevate the degree

    // h-refine
    for(index_t i = 0; i< numHref; ++i)
      mp.patch(0).uniformRefine();

    mp_def = mp;
    gsInfo<<"alpha = "<<aDim/bDim<<"; beta = "<<bDim/thickness<<"\n";


    gsMultiBasis<> dbasis(mp);
    gsInfo<<"Basis (patch 0): "<< mp.patch(0).basis() << "\n";

    // Boundary conditions
    gsBoundaryConditions<> BCs;
    gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();

    std::string output = "solution";
    std::string dirname = "ArcLengthResults";

    gsMatrix<> writePoints(2,3);
    writePoints.col(0)<< 0.0,0.5;
    writePoints.col(1)<< 0.5,0.5;
    writePoints.col(2)<< 1.0,0.5;

    BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 2 - z
    BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 2 - z

    BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0 ,false,1);
    BCs.addCondition(boundary::east, condition_type::collapsed, 0, 0 ,false,0);

    gsVector<> point(2); point<< 1.0, 0.5 ;
    gsVector<> load (2); load << Load,0.0;
    pLoads.addLoad(point, load, 0 );

    dirname = dirname + "/MaterialTest";

    output =  "solution";
    wn = output + "data.txt";

    std::string commands = "mkdir -p " + dirname;
    const char *command = commands.c_str();
    system(command);

    // plot geometry
    if (plot)
      gsWriteParaview(mp,dirname + "/" + "mp",1000,true);
    if (write)
      initStepOutput(dirname + "/" + wn, writePoints);

    // Linear isotropic material model
    gsFunctionExpr<> force("0","0",2);
    gsConstantFunction<> t(thickness,2);
    gsConstantFunction<> E(E_modulus,2);
    gsConstantFunction<> nu(PoissonRatio,2);
    gsConstantFunction<> rho(Density,2);
    gsConstantFunction<> ratio(Ratio,2);

    mu = E_modulus / (2 * (1 + PoissonRatio));
    gsConstantFunction<> alpha1(A1,2);
    gsConstantFunction<> mu1(M1,2);
    gsConstantFunction<> alpha2(A2,2);
    gsConstantFunction<> mu2(M2,2);
    gsConstantFunction<> alpha3(A3,2);
    gsConstantFunction<> mu3(M3,2);

    index_t kmax = 1;

    std::vector<gsFunctionSet<> * > Gs(kmax);
    std::vector<gsFunctionSet<> * > Ts(kmax);
    std::vector<gsFunctionSet<> * > Phis(kmax);

    gsMatrix<> Gmat = gsCompositeMatrix(E_modulus,E_modulus,0.5 * E_modulus / (1+PoissonRatio),PoissonRatio,PoissonRatio);
    Gmat.resize(Gmat.rows()*Gmat.cols(),1);
    gsConstantFunction<> Gfun(Gmat,2);
    Gs[0] = &Gfun;

    gsConstantFunction<> phi;
    phi.setValue(0,2);

    Phis[0] = &phi;

    gsConstantFunction<> thicks(thickness/kmax,2);
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
            materialMatrix = new gsMaterialMatrixComposite<2,real_t>(mp,Ts,Gs,Phis);
        }
        else
        {
            parameters.resize(2);
            options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
            options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
            materialMatrix = getMaterialMatrix<2,real_t>(mp,t,parameters,rho,options);
        }
    }
    else
    {
        options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",material);
        options.addSwitch("Compressibility","Compressibility: (false): Imcompressible | (true): Compressible",Compressibility);
        options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",impl);
        materialMatrix = getMaterialMatrix<2,real_t>(mp,t,parameters,rho,options);
    }

    gsThinShellAssemblerBase<real_t>* assembler;
    assembler = new gsThinShellAssembler<2,real_t,false>(mp,dbasis,BCs,force,materialMatrix);

    // Construct assembler object
    assembler->setPointLoads(pLoads);

    gsStopwatch stopwatch;
    real_t time = 0.0;

    typedef std::function<gsSparseMatrix<real_t> (gsVector<real_t> const &)>                                Jacobian_t;
    typedef std::function<gsVector<real_t> (gsVector<real_t> const &, real_t, gsVector<real_t> const &) >   ALResidual_t;
    // Function for the Jacobian
    Jacobian_t Jacobian = [&time,&stopwatch,&assembler,&mp_def](gsVector<real_t> const &x)
    {
      stopwatch.restart();
      assembler->constructSolution(x,mp_def);
      assembler->assembleMatrix(mp_def);
      time += stopwatch.stop();

      gsSparseMatrix<real_t> m = assembler->matrix();
      return m;
    };
    // Function for the Residual
    ALResidual_t ALResidual = [&time,&stopwatch,&assembler,&mp_def](gsVector<real_t> const &x, real_t lam, gsVector<real_t> const &force)
    {
      stopwatch.restart();
      assembler->constructSolution(x,mp_def);
      assembler->assembleVector(mp_def);
      gsVector<real_t> Fint = -(assembler->rhs() - force);
      gsVector<real_t> result = Fint - lam * force;
      time += stopwatch.stop();
      return result; // - lam * force;
    };
    // Assemble linear system to obtain the force vector
    assembler->assemble();
    gsVector<> Force = assembler->rhs();


    gsArcLengthIterator<real_t> arcLength(Jacobian, ALResidual, Force);

    arcLength.options().setInt("Solver",0); // LDLT solver
    arcLength.options().setInt("BifurcationMethod",0); // 0: determinant, 1: eigenvalue
    arcLength.options().setInt("Method",method);
    arcLength.options().setReal("Length",dL);
    arcLength.options().setInt("AngleMethod",0); // 0: step, 1: iteration
    arcLength.options().setSwitch("AdaptiveLength",adaptive);
    arcLength.options().setInt("AdaptiveIterations",5);
    arcLength.options().setReal("Scaling",0.0);
    arcLength.options().setReal("Tol",tol);
    arcLength.options().setReal("TolU",tolU);
    arcLength.options().setReal("TolF",tolF);
    arcLength.options().setInt("MaxIter",maxit);
    arcLength.options().setSwitch("Verbose",true);
    arcLength.options().setReal("Relaxation",relax);
    if (quasiNewtonInt>0)
    {
      quasiNewton = true;
      arcLength.options().setInt("QuasiIterations",quasiNewtonInt);
    }
    arcLength.options().setSwitch("Quasi",quasiNewton);


    gsInfo<<arcLength.options();
    arcLength.applyOptions();
    arcLength.initialize();

    gsParaviewCollection collection(dirname + "/" + output);
    gsParaviewCollection PStretch(dirname + "/" + "stretch");
    gsParaviewCollection PStretch1dir(dirname + "/" + "stretch1dir");
    gsParaviewCollection PStretch2dir(dirname + "/" + "stretch2dir");
    gsParaviewCollection Smembrane(dirname + "/" + "membrane");
    gsParaviewCollection Smembrane_p(dirname + "/" + "membrane_p");
    gsMultiPatch<> deformation = mp;

    gsMatrix<> solVector;
    for (index_t k=0; k<step; k++)
    {
      gsInfo<<"Load step "<< k<<"\n";
      arcLength.step();

      if (!(arcLength.converged()))
        GISMO_ERROR("Error: Loop terminated, arc length method did not converge.\n");

      solVector = arcLength.solutionU();

      assembler->constructSolution(solVector,mp_def);

      deformation = mp_def;
      deformation.patch(0).coefs() -= mp.patch(0).coefs();// assuming 1 patch here

      gsInfo<<"Total ellapsed assembly time: "<<time<<" s\n";

      if (plot)
      {
        gsField<> solField;
        if (deformed)
          solField= gsField<>(mp_def,deformation);
        else
          solField= gsField<>(mp,deformation);

        std::string fileName = dirname + "/" + output + util::to_string(k);
        gsWriteParaview<>(solField, fileName, 1000,true);
        fileName = output + util::to_string(k) + "0";
        collection.addTimestep(fileName,k,".vts");
        collection.addTimestep(fileName,k,"_mesh.vtp");
      }
      if (stress)
      {
        gsField<> stretch, stretch1dir, stretch2dir, membraneStress, membraneStress_p;

        gsPiecewiseFunction<> stretches;
        assembler->constructStress(mp_def,stretches,stress_type::principal_stretch);
        if (deformed)
          stretch = gsField<>(mp_def,stretches, true);
        else
          stretch = gsField<>(mp,stretches, true);

        gsPiecewiseFunction<> stretch1dirs;
        assembler->constructStress(mp_def,stretch1dirs,stress_type::principal_stretch_dir1);
        if (deformed)
          stretch1dir = gsField<>(mp_def,stretch1dirs, true);
        else
          stretch1dir = gsField<>(mp,stretch1dirs, true);

        gsPiecewiseFunction<> stretch2dirs;
        assembler->constructStress(mp_def,stretch2dirs,stress_type::principal_stretch_dir2);
        if (deformed)
          stretch2dir = gsField<>(mp_def,stretch2dirs, true);
        else
          stretch2dir = gsField<>(mp,stretch2dirs, true);

        gsPiecewiseFunction<> membraneStresses;
        assembler->constructStress(mp_def,membraneStresses,stress_type::membrane);
        if (deformed)
          membraneStress = gsField<>(mp_def,membraneStresses,true);
        else
          membraneStress = gsField<>(mp,membraneStresses,true);

        gsPiecewiseFunction<> membraneStresses_p;
        assembler->constructStress(mp_def,membraneStresses_p,stress_type::principal_stress_membrane);
        if (deformed)
          membraneStress_p = gsField<>(mp_def,membraneStresses_p, true);
        else
          membraneStress_p = gsField<>(mp,membraneStresses_p, true);


        std::string fileName;
        fileName = dirname + "/" + "stretch" + util::to_string(k);
        gsWriteParaview( stretch, fileName, 5000);
        fileName = "stretch" + util::to_string(k) + "0";
        PStretch.addTimestep(fileName,k,".vts");

        fileName = dirname + "/" + "stretch1dir" + util::to_string(k);
        gsWriteParaview( stretch1dir, fileName, 5000);
        fileName = "stretch1dir" + util::to_string(k) + "0";
        PStretch1dir.addTimestep(fileName,k,".vts");

        fileName = dirname + "/" + "stretch2dir" + util::to_string(k);
        gsWriteParaview( stretch2dir, fileName, 5000);
        fileName = "stretch2dir" + util::to_string(k) + "0";
        PStretch2dir.addTimestep(fileName,k,".vts");

        fileName = dirname + "/" + "membrane" + util::to_string(k);
        gsWriteParaview( membraneStress, fileName, 5000);
        fileName = "membrane" + util::to_string(k) + "0";
        Smembrane.addTimestep(fileName,k,".vts");

        fileName = dirname + "/" + "membrane_p" + util::to_string(k);
        gsWriteParaview( membraneStress_p, fileName, 5000);
        fileName = "membrane_p" + util::to_string(k) + "0";
        Smembrane_p.addTimestep(fileName,k,".vts");

      }

      if (write)
        writeStepOutput(arcLength,deformation, dirname + "/" + wn, writePoints);
    }

    if (plot)
    {
      collection.save();
    }
    if (stress)
    {
      PStretch.save();
      PStretch1dir.save();
      PStretch2dir.save();
      Smembrane.save();
      Smembrane_p.save();
    }

  return result;
}


template <class T>
gsMultiPatch<T> Rectangle(T L, T B) //, int n, int m, std::vector<boxSide> sides, T offset)
{
  // -------------------------------------------------------------------------
  // --------------------------Make beam geometry-----------------------------
  // -------------------------------------------------------------------------
  int dim = 2; //physical dimension
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
void initStepOutput(const std::string name, const gsMatrix<T> & points)
{
  std::ofstream file;
  file.open(name,std::ofstream::out);
  file  << std::setprecision(6)
        << "Deformation norm" << ",";
        for (index_t k=0; k!=points.cols(); k++)
        {
          file<< "point "<<k<<" - x" << ","
              << "point "<<k<<" - y" << ",";
        }

  file  << "Lambda" << ","
        << "Indicator"
        << "\n";
  file.close();

  gsInfo<<"Step results will be written in file: "<<name<<"\n";
}

template <class T>
void writeStepOutput(const gsArcLengthIterator<T> & arcLength, const gsMultiPatch<T> & deformation, const std::string name, const gsMatrix<T> & points)
{
  gsMatrix<T> P(2,1), Q(2,1);
  gsMatrix<T> out(2,points.cols());
  gsMatrix<T> tmp;

  for (index_t p=0; p!=points.cols(); p++)
  {
    P<<points.col(p);
    deformation.patch(0).eval_into(P,tmp);
    out.col(p) = tmp;
  }

  std::ofstream file;
  file.open(name,std::ofstream::out | std::ofstream::app);
  file  << std::setprecision(12)
        << arcLength.solutionU().norm() << ",";
        for (index_t p=0; p!=points.cols(); p++)
        {
          file<< out(0,p) << ","
              << out(1,p) << ",";
        }

  file  << arcLength.solutionL() << ","
        << arcLength.indicator() << ","
        << "\n";

  file.close();
}

