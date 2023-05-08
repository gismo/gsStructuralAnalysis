/** @file benchmark_Balloon.cpp

    @brief Benchmark for the inflated pillow using the Arc-Length Method

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

int main (int argc, char** argv)
{
    // Input options
    int numElevate    = 0;
    int numRefine     = 1;
    bool plot         = false;
    bool stress       = false;
    bool quasiNewton  = false;
    int quasiNewtonInt= -1;
    bool adaptive     = false;
    int step          = 21;
    int method        = 2; // (0: Load control; 1: Riks' method; 2: Crisfield's method; 3: consistent crisfield method; 4: extended iterations)
    bool membrane = false;

    real_t relax      = 1.0;
    int result        = 0;
    index_t maxit     = 50;
    // Arc length method options
    real_t dL         = 5e0; // General arc length
    real_t tol        = 1e-6;
    real_t tolU       = 1e-6;
    real_t tolF       = 1e-3;

    std::string wn("data.csv");

    gsCmdLine cmd("Example for an inflating balloon.");

    cmd.addInt("r","hRefine", "Number of dyadic h-refinement (bisection) steps to perform before solving", numRefine);
    cmd.addInt("e","degreeElevation", "Number of degree elevation steps to perform on the Geometry's basis before solving", numElevate);
    cmd.addSwitch("membrane", "Membrane element", membrane);

    cmd.addInt("m","Method", "Arc length method; 1: Crisfield's method; 2: RIks' method.", method);
    cmd.addReal("L","dL", "arc length", dL);
    cmd.addReal("A","relaxation", "Relaxation factor for arc length method", relax);

    cmd.addInt("q","QuasiNewtonInt","Use the Quasi Newton method every INT iterations",quasiNewtonInt);
    cmd.addInt("N", "maxsteps", "Maximum number of steps", step);

    cmd.addSwitch("adaptive", "Adaptive length ", adaptive);
    cmd.addSwitch("quasi", "Use the Quasi Newton method", quasiNewton);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addSwitch("stress", "Plot stress in ParaView format", stress);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    /*
      Uniaxial tension of a square plate                            --- Validation settings: -L 1eX -l 1eX -M 14 -N 500 -r X -e X
      (bottom boundary fixed in Y, left boundary fixed in X, right boundary normal load)
    */
    real_t E_modulus;
    real_t PoissonRatio;
    real_t mu;
    real_t Ratio;
    real_t thickness = 0.1e-3;
    real_t Density    = 1e0;
    E_modulus = 2e9; // 2GPa
    PoissonRatio = 0.3;
    Ratio = 0;
    // }


    gsMultiPatch<> mp,mp_def;
    gsMatrix<> box(2,2);
    box.col(0)<<0,0;
    box.col(1)<<0.5,0.5;
    mp.addPatch(gsNurbsCreator<>::BSplineSquare(box));
    mp.embed(3);

    // p-refine
    for(index_t i = 0; i< numElevate; ++i)
      mp.patch(0).degreeElevate();    // Elevate the degree

    // h-refine
    for(index_t i = 0; i< numRefine; ++i)
      mp.patch(0).uniformRefine();

    gsInfo<<"mu = "<<E_modulus / (2 * (1 + PoissonRatio))<<"\n";

    gsMultiBasis<> dbasis(mp,true);
    gsInfo<<"Basis (patch 0): "<< mp.patch(0).basis() << "\n";
    mp_def = mp;

    // Boundary conditions
    gsBoundaryConditions<> BCs;
    BCs.setGeoMap(mp);

    BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2);
    BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 2);

    BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 );
    BCs.addCondition(boundary::north, condition_type::clamped, 0, 0, false, 0 );
    BCs.addCondition(boundary::north, condition_type::clamped, 0, 0, false, 2 );

    // Symmetry in y-direction:
    BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 );
    BCs.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 1 );
    BCs.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 2 );

    BCs.setGeoMap(mp);

    // Pressure
    real_t pressure = 5e3;

    std::string dirname = "ArcLengthResults";
    dirname = dirname + "/Pillow";
    std::string output =  "solution";

    std::string commands = "mkdir -p " + dirname;
    const char *command = commands.c_str();
    int systemRet = system(command);
    GISMO_ASSERT(systemRet!=-1,"Something went wrong with calling the system argument");


    // plot geometry
    if (plot)
      gsWriteParaview(mp,dirname + "/" + "mp",1000,true);


    // Linear isotropic material model
    gsFunctionExpr<> force("0","0","0",3);
    gsConstantFunction<> pressFun(pressure,3);
    gsConstantFunction<> t(thickness,3);
    gsConstantFunction<> E(E_modulus,3);
    gsConstantFunction<> nu(PoissonRatio,3);
    gsConstantFunction<> rho(Density,3);
    gsConstantFunction<> ratio(Ratio,3);

    gsMaterialMatrixBase<real_t>* materialMatrix;
    materialMatrix = new gsMaterialMatrixLinear<3,real_t>(mp,t,E,nu,rho);

    gsThinShellAssemblerBase<real_t>* assembler;
    if (membrane)
      assembler = new gsThinShellAssembler<3, real_t, true >(mp,dbasis,BCs,force,materialMatrix);
    else
      assembler = new gsThinShellAssembler<3, real_t, false >(mp,dbasis,BCs,force,materialMatrix);

    // Construct assembler object
    // assembler->setPressure(pressFun);

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
    ALResidual_t ALResidual = [&pressure,&time,&stopwatch,&assembler,&mp_def](gsVector<real_t> const &x, real_t lam, gsVector<real_t> const &force)
    {
      stopwatch.restart();
      assembler->constructSolution(x,mp_def);
      assembler->assembleVector(mp_def);
      gsVector<real_t> F = assembler->rhs();
      assembler->assemblePressureVector(pressure,mp_def);
      F += assembler->rhs();
      gsVector<real_t> Fint = -(F - force);
      gsVector<real_t> result = Fint - lam * force;
      time += stopwatch.stop();
      return result; // - lam * force;
    };
    // Assemble linear system to obtain the force vector
    assembler->assemble();
    gsVector<> Force = assembler->rhs();
    assembler->assemblePressureVector(pressure,mp_def);
    Force += assembler->rhs();

    gsALMBase<real_t> * arcLength;
    if (method==0)
      arcLength = new gsALMLoadControl<real_t>(Jacobian, ALResidual, Force);
    else if (method==1)
      arcLength = new gsALMRiks<real_t>(Jacobian, ALResidual, Force);
    else if (method==2)
      arcLength = new gsALMCrisfield<real_t>(Jacobian, ALResidual, Force);
    else
      GISMO_ERROR("Method "<<method<<" unknown");

    arcLength->options().setString("Solver","SimplicialLDLT");
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

    // Make objects for previous solutions
    gsMatrix<> Uold = Force;
    Uold.setZero();

    gsMatrix<> solVector;
    real_t indicator = 0.0;
    arcLength->setIndicator(indicator); // RESET INDICATOR

    gsMatrix<> lambdas(3,step);
    gsMatrix<> pressures(2,step);
    for (index_t k=0; k<step; k++)
    {
      gsInfo<<"Load step "<< k<<"\n";
      arcLength->step();

      if (!(arcLength->converged()))
        GISMO_ERROR("Loop terminated, arc length method did not converge.\n");

      solVector = arcLength->solutionU();
      Uold = solVector;

      assembler->constructSolution(solVector,mp_def);

      gsMatrix<> pts(2,1);
      pts<<0.0,1.0;

      lambdas.col(k) = assembler->computePrincipalStretches(pts,mp_def,0);
      pressures.col(k) << pressure*arcLength->solutionL(), pressure*arcLength->solutionL() * assembler->getArea(mp) / assembler->getArea(mp_def);

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

      }
    }

    gsInfo<<"Lambdas:\n"<<lambdas<<"\n";
    gsInfo<<"Pressures:\n"<<pressures<<"\n";

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
