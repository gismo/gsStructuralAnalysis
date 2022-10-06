/** @file benchmark_FrustrumALM.cpp

    @brief Benchmark for the collapsing frustrum with the Arc-Length Method

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

#include <gsKLShell/gsThinShellAssembler.h>
#include <gsKLShell/getMaterialMatrix.h>

#include <gsStructuralAnalysis/gsALMBase.h>
#include <gsStructuralAnalysis/gsALMLoadControl.h>
#include <gsStructuralAnalysis/gsALMRiks.h>
#include <gsStructuralAnalysis/gsALMCrisfield.h>

using namespace gismo;

template <class T>
gsMultiPatch<T> FrustrumDomain(int n, int p, T R1, T R2, T h);

template <class T>
void initStepOutput( const std::string name, const gsMatrix<T> & points);

template <class T>
void writeStepOutput(const gsALMBase<T> * arcLength, const gsMultiPatch<T> & deformation, const std::string name, const gsMatrix<T> & points, const index_t extreme=-1, const index_t kmax=100);

int main (int argc, char** argv)
{
    // Input options
    int numElevate    = 0;
    int numRefine     = 5;
    bool plot         = false;
    bool write         = false;
    bool stress       = false;
    bool quasiNewton  = false;
    int quasiNewtonInt= -1;
    bool adaptive     = false;
    int step          = 10;
    int method        = 2; // (0: Load control; 1: Riks' method; 2: Crisfield's method; 3: consistent crisfield method; 4: extended iterations)

    int testCase = 0;

    index_t material  = 3;
    bool composite = false;
    index_t impl = 1; // 1= analytical, 2= generalized, 3= spectral

    real_t relax      = 1.0;
    int result        = 0;
    index_t maxit     = 50;
    // Arc length method options
    real_t dL         = -1; // General arc length
    real_t tol        = 1e-6;
    real_t tolU       = 1e-6;
    real_t tolF       = 1e-3;

    std::string wn("data.csv");

    gsCmdLine cmd("Arc-length analysis of a collapsing frustrum.");

    cmd.addInt("t", "testcase", "Test case: 0: free; 1: restrained", testCase);
    cmd.addInt("r","hRefine", "Number of dyadic h-refinement (bisection) steps to perform before solving", numRefine);
    cmd.addInt("e","degreeElevation", "Number of degree elevation steps to perform on the Geometry's basis before solving", numElevate);
    cmd.addInt( "M", "Material", "Material law",  material );
    cmd.addInt( "I", "Implementation", "Implementation: 1= analytical, 2= generalized, 3= spectral",  impl );
    cmd.addSwitch("composite", "Composite material", composite);

    cmd.addInt("m","Method", "Arc length method; 1: Crisfield's method; 2: RIks' method.", method);
    cmd.addReal("L","dL", "arc length", dL);
    cmd.addReal("A","relaxation", "Relaxation factor for arc length method", relax);

    cmd.addInt("q","QuasiNewtonInt","Use the Quasi Newton method every INT iterations",quasiNewtonInt);
    cmd.addInt("N", "maxsteps", "Maximum number of steps", step);

    cmd.addSwitch("adaptive", "Adaptive length ", adaptive);
    cmd.addSwitch("quasi", "Use the Quasi Newton method", quasiNewton);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addSwitch("stress", "Plot stress in ParaView format", stress);
    cmd.addSwitch("write", "write to file", write);

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
    mp = FrustrumDomain(numRefine,numElevate+2,2.0,1.0,1.0);
    // Initialise solution object
    mp_def = mp;

    gsMultiBasis<> dbasis(mp);
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

    gsMatrix<> writePoints(2,3);
    writePoints.col(0)<<0.0,1.0;
    writePoints.col(1)<<0.5,1.0;
    writePoints.col(2)<<1.0,1.0;

    /*
      Case 0: Frustrum with constrained top boundary        (Complete cycle with -N 1200)
      Case 1: Frustrum with unconstrained top boundary      (Complete cycle with -N 850)
    */
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

      dirname = dirname + "/" + "Frustrum_-r=" + std::to_string(numRefine) + "-e" + std::to_string(numElevate) + "-M" + std::to_string(material) + "_solution";
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

      dirname = dirname + "/" + "Frustrum2_-r=" + std::to_string(numRefine) + "-e" + std::to_string(numElevate) + "-M" + std::to_string(material) + "_solution";
      if (dL == -1) { dL = 4e-2; }
    }
    else
      GISMO_ERROR("Test case" << testCase << " does not exist!");

    wn = output + "data.txt";

    std::string commands = "mkdir -p " + dirname;
    const char *command = commands.c_str();
    int systemRet = system(command);
    GISMO_ASSERT(systemRet!=-1,"Something went wrong with calling the system argument");

    // plot geometry
    if (plot)
      gsWriteParaview(mp,dirname + "/" + "mp",1000,true);

    if (write)
      initStepOutput(dirname + "/" + wn, writePoints);

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


    gsDebug<<arcLength->options();
    arcLength->applyOptions();
    arcLength->initialize();


    gsParaviewCollection collection(dirname + "/" + output);
    gsParaviewCollection Smembrane(dirname + "/" + "membrane");
    gsParaviewCollection Sflexural(dirname + "/" + "flexural");
    gsParaviewCollection Smembrane_p(dirname + "/" + "membrane_p");
    gsMultiPatch<> deformation = mp;

    gsMatrix<> solVector;
    real_t indicator = 0.0;
    arcLength->setIndicator(indicator); // RESET INDICATOR

    for (index_t k=0; k<step; k++)
    {
      gsInfo<<"Load step "<< k<<"\n";
      arcLength->step();

      if (!(arcLength->converged()))
        GISMO_ERROR("Loop terminated, arc length method did not converge.\n");

      indicator = arcLength->indicator();
      solVector = arcLength->solutionU();

      assembler->constructSolution(solVector,mp_def);

      deformation = mp_def;
      deformation.patch(0).coefs() -= mp.patch(0).coefs();// assuming 1 patch here

      gsInfo<<"Total ellapsed assembly time: "<<time<<" s\n";

      if (plot)
      {
        gsField<> solField;
        solField= gsField<>(mp,deformation);

        std::string fileName = dirname + "/" + output + util::to_string(k);
        gsWriteParaview<>(solField, fileName, 1000,true);
        fileName = output + util::to_string(k) + "0";
        collection.addTimestep(fileName,k,".vts");
        collection.addTimestep(fileName,k,"_mesh.vtp");
      }
      if (stress)
      {
        std::string fileName;

        gsField<> membraneStress, flexuralStress, membraneStress_p;

        gsPiecewiseFunction<> membraneStresses;
        assembler->constructStress(mp_def,membraneStresses,stress_type::membrane);
        membraneStress = gsField<>(mp,membraneStresses,true);

        fileName = dirname + "/" + "membrane" + util::to_string(k);
        gsWriteParaview( membraneStress, fileName, 1000);
        fileName = "membrane" + util::to_string(k) + "0";
        Smembrane.addTimestep(fileName,k,".vts");

        gsPiecewiseFunction<> flexuralStresses;
        assembler->constructStress(mp_def,flexuralStresses,stress_type::flexural);
        flexuralStress = gsField<>(mp,flexuralStresses, true);

        fileName = dirname + "/" + "flexural" + util::to_string(k);
        gsWriteParaview( flexuralStress, fileName, 1000);
        fileName = "flexural" + util::to_string(k) + "0";
        Sflexural.addTimestep(fileName,k,".vts");

        if (impl==3)
        {
          gsPiecewiseFunction<> membraneStresses_p;
          assembler->constructStress(mp_def,membraneStresses_p,stress_type::principal_stress_membrane);
          membraneStress_p = gsField<>(mp,membraneStresses_p, true);

          fileName = dirname + "/" + "membrane_p" + util::to_string(k);
          gsWriteParaview( membraneStress_p, fileName, 1000);
          fileName = "membrane_p" + util::to_string(k) + "0";
          Smembrane_p.addTimestep(fileName,k,".vts");
        }
      }
      if (write)
        writeStepOutput(arcLength,deformation, dirname + "/" + wn, writePoints,1, 201);

    }

    if (plot)
    {
      collection.save();
      Smembrane.save();
      Sflexural.save();
      Smembrane_p.save();
    }

    delete materialMatrix;
    delete assembler;
    delete arcLength;

  return result;
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
