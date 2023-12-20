/** @file benchmark_Roof.cpp

    @brief Collapse of a cylindrical roof, from

    Guo, Y., Do, H., & Ruess, M. (2019). Isogeometric stability analysis of thin shells:
    From simple geometries to engineering models.
    International Journal for Numerical Methods in Engineering.
    https://doi.org/10.1002/nme.6020

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
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMRiks.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMLoadControl.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMCrisfield.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMConsistentCrisfield.h>

using namespace gismo;

template <class T>
void initStepOutput( const std::string name, const gsMatrix<T> & points);

template <class T>
void writeStepOutput(const gsALMBase<T> * arcLength, const gsMultiPatch<T> & deformation, const std::string name, const gsMatrix<T> & points, const index_t extreme=-1, const index_t kmax=100);

#ifdef gsKLShell_ENABLED
int main (int argc, char** argv)
{
    // Input options
    int numElevate    = 1;
    int numHref       = 1;
    bool plot         = false;
    bool mesh         = false;
    bool stress       = false;
    bool membrane     = false;
    bool quasiNewton  = false;
    int quasiNewtonInt= -1;
    bool adaptive     = false;
    int step          = 10;
    int method        = 2; // (0: Load control; 1: Riks' method; 2: Crisfield's method; 3: consistent crisfield method)
    bool deformed     = false;
    bool symmetry = false;

    bool composite = false;

    real_t relax      = 1.0;

    int testCase      = 0;

    int result        = 0;

    bool write        = false;

    index_t maxit     = 20;

    // Arc length method options
    real_t dL        = 0.5; // Ard length to find bifurcation
    real_t tol        = 1e-6;
    real_t tolU       = 1e-6;
    real_t tolF       = 1e-3;

    std::string wn("data.csv");

    std::string assemberOptionsFile("options/solver_options.xml");

    gsCmdLine cmd("Arc-length analysis of a collapsing roof.");
    cmd.addString( "f", "file", "Input XML file for assembler options", assemberOptionsFile );

    cmd.addInt("t", "testcase", "Test case: 0: clamped-clamped, 1: pinned-pinned, 2: clamped-free", testCase);

    cmd.addInt("r","hRefine", "Number of dyadic h-refinement (bisection) steps to perform before solving", numHref);
    cmd.addInt("e","degreeElevation", "Number of degree elevation steps to perform on the Geometry's basis before solving", numElevate);
    cmd.addSwitch("composite", "Composite material", composite);

    cmd.addInt("m","Method", "Arc length method; 1: Crisfield's method; 2: RIks' method.", method);
    cmd.addReal("L","dL", "arc length", dL);
    cmd.addReal("A","relaxation", "Relaxation factor for arc length method", relax);

    cmd.addInt("q","QuasiNewtonInt","Use the Quasi Newton method every INT iterations",quasiNewtonInt);
    cmd.addInt("N", "maxsteps", "Maximum number of steps", step);

    cmd.addSwitch("adaptive", "Adaptive length ", adaptive);
    cmd.addSwitch("quasi", "Use the Quasi Newton method", quasiNewton);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addSwitch("mesh", "Plot mesh?", mesh);
    cmd.addSwitch("stress", "Plot stress in ParaView format", stress);
    cmd.addSwitch("membrane", "Use membrane model (no bending)", membrane);
    cmd.addSwitch("deformed", "plot on deformed shape", deformed);
    cmd.addSwitch("write", "write to file", write);
    cmd.addSwitch("symmetry", "Apply symmetry?", symmetry);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

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

    if (!symmetry)
        gsReadFile<>("surface/scordelis_lo_roof_shallow.xml", mp);
    else
        gsReadFile<>("surface/scordelis_lo_roof_shallow_symm.xml", mp);

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

    gsMatrix<> writePoints;
    if (!symmetry)
    {
      writePoints.resize(2,3);
      writePoints.col(0)<< 0.0,0.5;
      writePoints.col(1)<< 0.5,0.5;
      writePoints.col(2)<< 1.0,0.5;
    }
    else
    {
      writePoints.resize(2,2);
      writePoints.col(0)<< 0.0,0.0;
      writePoints.col(1)<< 1.0,0.0;
    }

    GISMO_ASSERT(mp.targetDim()==3,"Geometry must be surface (targetDim=3)!");
    if (!symmetry)
    {
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
    }
    else
    {
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

        BCs.addCondition(boundary::north, condition_type::clamped, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::north, condition_type::clamped, 0, 0, false, 2 ); // unknown 2 - z

        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 2 ); // unknown 2 - z
    }

    BCs.setGeoMap(mp);

    Load = -1e1;
    // Point loads
    gsVector<> point(2);
    gsVector<> load (3);
    if (!symmetry)
        point<< 0.5, 0.5 ;
    else
        point<< 1.0, 1.0 ;

    load << 0.0, 0.0, Load ;
    pLoads.addLoad(point, load, 0 );

    dirname = dirname + "/" +  "Roof_t="+ std::to_string(thickness) + "-r=" + std::to_string(numHref) + "-e" + std::to_string(numElevate) +"_solution";
    if (symmetry)
        dirname = dirname + "_symm";
    output =  "solution";
    wn = output + "data.txt";

    std::string commands = "mkdir -p " + dirname;
    const char *command = commands.c_str();
    int systemRet = system(command);
    GISMO_ASSERT(systemRet!=-1,"Something went wrong with calling the system argument");


    // plot geometry
    if (plot)
      gsWriteParaview(mp,dirname + "/" + "mp",1000,mesh);

    if (write)
      initStepOutput(dirname + "/" + wn, writePoints);

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
    // phi1.setValue(pi/2,3);
    // phi2.setValue(0,3);
    // phi3.setValue(pi/2,3);
    phi1.setValue(0,3);
    phi2.setValue(pi/2,3);
    phi3.setValue(0,3);

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

    std::vector<gsFunctionSet<>*> parameters;
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

    // (0: Load control; 1: Riks' method; 2: Crisfield's method; 3: consistent crisfield method)
    gsALMBase<real_t> * arcLength;
    if (method==0)
      arcLength = new gsALMLoadControl<real_t>(Jacobian, ALResidual, Force);
    else if (method==1)
      arcLength = new gsALMRiks<real_t>(Jacobian, ALResidual, Force);
    else if (method==2)
      arcLength = new gsALMCrisfield<real_t>(Jacobian, ALResidual, Force);
    else if (method==3)
      arcLength = new gsALMConsistentCrisfield<real_t>(Jacobian, ALResidual, Force);
    else
      GISMO_ERROR("Method unknown");

    arcLength->options().setString("Solver","SimplicialLDLT"); // LDLT solver
    arcLength->options().setInt("BifurcationMethod",0); // 0: determinant, 1: eigenvalue
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


    gsInfo<<arcLength->options();
    arcLength->applyOptions();
    arcLength->initialize();


    gsParaviewCollection collection(dirname + "/" + output);
    gsParaviewCollection Smembrane(dirname + "/" + "membrane");
    gsParaviewCollection Sflexural(dirname + "/" + "flexural");
    gsParaviewCollection Smembrane_p(dirname + "/" + "membrane_p");
    gsMultiPatch<> deformation = mp;

    // Make objects for previous solutions
    real_t Lold = 0;
    gsMatrix<> Uold = Force;
    Uold.setZero();

    gsMatrix<> solVector;
    real_t indicator = 0.0;
    arcLength->setIndicator(indicator); // RESET INDICATOR
    bool bisected = false;
    real_t dL0 = dL;

    gsVector<> pt(2);
    pt<<0.5,0.5;
    gsMatrix<> pt_result;
    gsMultiPatch<> mp_tmp;
    /////////////////////////////////////////////////////////////////////////////////////////////
    std::ofstream file;
    file.open("sol.txt",std::ofstream::out);
    assembler->constructSolution(Uold,mp_tmp);
    deformation.patch(0) = mp_tmp.patch(0);
    deformation.patch(0).coefs() -= mp.patch(0).coefs();// assuming 1 patch here
    deformation.patch(0).eval_into(pt,pt_result);
    file  << std::setprecision(20) << pt_result(2,0) << "," << Lold<<"\n";
    //////////////////////////////////////////////////////////////////////////////////////////////


    for (index_t k=0; k<step; k++)
    {
      gsInfo<<"Load step "<< k<<"\n";
        gsInfo<<"dL = "<<dL<<"\n";
      // assembler->constructSolution(solVector,solution);

      gsStatus status = arcLength->step();
      if (status==gsStatus::NotConverged || status==gsStatus::AssemblyError)
      {
        gsInfo<<"Error: Loop terminated, arc length method did not converge.\n";
        dL = dL / 2.;
        arcLength->setLength(dL);
        arcLength->setSolution(Uold,Lold);
        bisected = true;
        k -= 1;
        continue;
      }

      indicator = arcLength->indicator();

      solVector = arcLength->solutionU();
      Uold = solVector;
      Lold = arcLength->solutionL();
      assembler->constructSolution(solVector,mp_def);

      //////////////////////////////////////////////////////////////////////////////////////////////
      assembler->constructSolution(Uold,mp_tmp);
      deformation.patch(0) = mp_tmp.patch(0);
      deformation.patch(0).coefs() -= mp.patch(0).coefs();// assuming 1 patch here
      deformation.patch(0).eval_into(pt,pt_result);
      file  << std::setprecision(20) << pt_result(2,0) << "," << Lold<<"\n";
      //////////////////////////////////////////////////////////////////////////////////////////////



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
        gsWriteParaview<>(solField, fileName, 1000,mesh);
        fileName = output + util::to_string(k) + "0";
        collection.addPart(fileName + ".vts",k);
        if (mesh) collection.addPart(fileName + "_mesh.vtp",k);
      }
      if (stress)
      {
        std::string fileName;

        gsField<> membraneStress, flexuralStress, membraneStress_p;

        gsPiecewiseFunction<> membraneStresses;
        assembler->constructStress(mp_def,membraneStresses,stress_type::membrane);
        if (deformed)
          membraneStress = gsField<>(mp_def,membraneStresses,true);
        else
          membraneStress = gsField<>(mp,membraneStresses,true);

        fileName = dirname + "/" + "membrane" + util::to_string(k);
        gsWriteParaview( membraneStress, fileName, 1000);
        fileName = "membrane" + util::to_string(k) + "0";
        Smembrane.addPart(fileName + ".vts",k);

        gsPiecewiseFunction<> flexuralStresses;
        assembler->constructStress(mp_def,flexuralStresses,stress_type::flexural);
        if (deformed)
          flexuralStress = gsField<>(mp_def,flexuralStresses, true);
        else
          flexuralStress = gsField<>(mp,flexuralStresses, true);

        fileName = dirname + "/" + "flexural" + util::to_string(k);
        gsWriteParaview( flexuralStress, fileName, 1000);
        fileName = "flexural" + util::to_string(k) + "0";
        Sflexural.addPart(fileName + ".vts",k);
      }

      if (write)
        writeStepOutput(arcLength,deformation, dirname + "/" + wn, writePoints,1, 201);

      if (!bisected)
      {
        dL = dL0;
        arcLength->setLength(dL);
      }
      bisected = false;

    }

    if (plot)
    {
      collection.save();
      Smembrane.save();
      Sflexural.save();
      Smembrane_p.save();
    }

    //////////////////////////////////////////////////////////////////////////////////////////////
    file.close();
    //////////////////////////////////////////////////////////////////////////////////////////////


    delete arcLength;
    delete materialMatrix;
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
