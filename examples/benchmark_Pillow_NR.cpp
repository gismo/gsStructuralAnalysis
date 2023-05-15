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
#include <gsKLShell/gsMaterialMatrixTFT.h>

#include <gsStructuralAnalysis/gsStaticDR.h>
#include <gsStructuralAnalysis/gsStaticNewton.h>

using namespace gismo;

int main (int argc, char** argv)
{
    // Input options
    int numElevate    = 0;
    int numRefine     = 1;
    bool plot         = false;
    bool stress       = false;
    bool TFT       = false;
    bool quasiNewton  = false;
    int quasiNewtonInt= -1;
    bool adaptive     = false;
    int step          = 21;
    int method        = 2; // (0: Load control; 1: Riks' method; 2: Crisfield's method; 3: consistent crisfield method; 4: extended iterations)
    int verbose = 0;

    bool membrane = false;
    bool nonfollow     = false;

    int result        = 0;
    index_t maxIt     = 1e3;
    // Arc length method options
    real_t tol        = 1e-6;

    real_t alpha = 1.0;
    real_t damping = 0.1;

    real_t dt = 0.1;

    index_t testCase = 0;

    std::string wn("data.csv");

    gsCmdLine cmd("Example for an inflating balloon.");

    cmd.addInt("r","hRefine", "Number of dyadic h-refinement (bisection) steps to perform before solving", numRefine);
    cmd.addInt("e","degreeElevation", "Number of degree elevation steps to perform on the Geometry's basis before solving", numElevate);
    cmd.addSwitch("membrane", "Membrane element", membrane);
    cmd.addSwitch("nonfollow", "No follower load", nonfollow);

    cmd.addInt( "N", "maxit", "maxit",  maxIt );
    cmd.addInt( "t", "testCase", "testCase",  testCase );

    cmd.addReal( "d", "dt", "dt",  dt );
    cmd.addReal( "a", "alpha", "alpha",  alpha );
    cmd.addReal( "c", "damping", "damping",  damping );

    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addSwitch("stress", "Plot stress in ParaView format", stress);
    cmd.addSwitch("TFT", "Use Tension Field Theory", TFT);
    cmd.addInt("v","verbose", "0: no; 1: iteration output; 2: Full matrix and vector output", verbose);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    /*
      Uniaxial tension of a square plate                            --- Validation settings: -L 1eX -l 1eX -M 14 -N 500 -r X -e X
      (bottom boundary fixed in Y, left boundary fixed in X, right boundary normal load)
    */
    real_t E_modulus;
    real_t PoissonRatio;
    real_t mu;
    real_t Ratio;
    real_t thickness;
    real_t Density    = 1e0;

    gsMatrix<> box(2,2);
    if (testCase == 0)
    {
      //Jarasjarungkiat2009
      E_modulus = 588e6;
      PoissonRatio = 0.4;
      thickness = 0.0006;
      box.col(0)<<0,0;
      real_t AM = 0.6;
      real_t AB = math::sqrt(math::pow(AM,2)/2.);
      gsDebugVar(AB);
      box.col(1)<<math::sqrt(0.18),math::sqrt(0.18);
    }
    if (testCase == 1)
    {

      E_modulus = 2e9; // 2GPa
      PoissonRatio = 0.3;
      thickness = 0.1e-3;
      box.col(0)<<0,0;
      box.col(1)<<0.5,0.5;
    }

    gsMultiPatch<> mp,mp_def;
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

    gsVector<> neuDataX(3);
    neuDataX<<5e3,0,0;
    gsConstantFunction<> neuX(neuDataX,3);
    gsVector<> neuDataY(3);
    neuDataY<<0,-5e3,0;
    gsConstantFunction<> neuY(neuDataY,3);

    // Boundary conditions
    gsBoundaryConditions<> BCs;


    BCs.setGeoMap(mp);

    BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2);
    BCs.addCondition(boundary::south, condition_type::neumann, &neuY);
    BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 2);
    BCs.addCondition(boundary::east, condition_type::neumann, &neuX);

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
    // real_t pressure =0;

    std::string dirname = "ArcLengthResults";
    dirname = dirname + "/Pillow";
    std::string output =  "solution";

    std::string commands = "mkdir -p " + dirname;
    const char *command = commands.c_str();
    int systemRet = system(command);
    GISMO_ASSERT(systemRet!=-1,"Something went wrong with calling the system argument");


    // plot geometry
    if (plot)
      gsWriteParaview(mp,"mp",1000,true);

    // Linear isotropic material model
    gsVector<> forceVec(3);
    forceVec.setZero();

    gsConstantFunction<> force(forceVec,3);
    gsConstantFunction<> pressFun(pressure,3);
    gsConstantFunction<> t(thickness,3);
    gsConstantFunction<> E(E_modulus,3);
    gsConstantFunction<> nu(PoissonRatio,3);
    gsConstantFunction<> rho(Density,3);
    gsConstantFunction<> ratio(Ratio,3);

    gsMaterialMatrixBase<real_t>* materialMatrix = new gsMaterialMatrixLinear<3,real_t>(mp,t,E,nu,rho);
    gsMaterialMatrixTFT<3,real_t,false> * materialMatrixTFT = new gsMaterialMatrixTFT<3,real_t,false>(static_cast<gsMaterialMatrixBaseDim<3,real_t> * >(materialMatrix));
    materialMatrixTFT->options().setReal("SlackMultiplier",1e-6);    

    gsThinShellAssemblerBase<real_t>* assembler;
    if (membrane && TFT)
      assembler = new gsThinShellAssembler<3, real_t, true >(mp,dbasis,BCs,force,materialMatrixTFT);
    else if (membrane && !TFT)
      assembler = new gsThinShellAssembler<3, real_t, true >(mp,dbasis,BCs,force,materialMatrix);
    else if (!membrane && TFT)
      assembler = new gsThinShellAssembler<3, real_t, false >(mp,dbasis,BCs,force,materialMatrixTFT);
    else
      assembler = new gsThinShellAssembler<3, real_t, false >(mp,dbasis,BCs,force,materialMatrix);

    // Construct assembler object
    gsStopwatch stopwatch;
    real_t time = 0.0;

    typedef std::function<gsSparseMatrix<real_t> (gsVector<real_t> const &)>                                Jacobian_t;
    typedef std::function<gsVector<real_t> (gsVector<real_t> const &) >   Residual_t;
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
    Residual_t Residual = [&pressFun,&assembler,&mp_def](gsVector<real_t> const &x)
    {
      assembler->constructSolution(x,mp_def);
      assembler->assembleVector(mp_def);
      gsVector<real_t> result = assembler->rhs();
      // assembler->assemblePressureVector(pressFun,mp_def);
      // result += assembler->rhs();
      return result; // - lam * force;
    };
    // Assemble linear system to obtain the force vector
    assembler->assemble();
    gsSparseMatrix<> K = assembler->matrix();
    gsVector<> Fneu = assembler->rhs();
    real_t Fnorm = Fneu.norm();

    gsSparseMatrix<> diag(K.rows(),K.cols());
    diag.setIdentity();
    diag *= 0.14;
    K+=diag;

    gsVector<> solVector(assembler->numDofs());
    solVector.setZero();

    /*
     * Load stepping:
     * - First 3-5 steps: K+0.14*I
     * - First 20 steps: tension on free boundaries of 4e3 N/m decreasing to 0 N/m
     * - First 20 steps: follower pressure increasing from 0 Pa to 5000 Pa
     */
    real_t load_fac = 0;
    real_t dload_fac= 1./maxIt;
    real_t dload_fac0 = dload_fac;
    index_t k=0;
    while (load_fac <= 1.)
    {
      /*
       * Load stepping:
       * - First 3-5 steps: K+0.14*I
       * - First 20 steps: tension on free boundaries of 4e3 N/m decreasing to 0 N/m
       * - First 20 steps: follower pressure increasing from 0 Pa to 5000 Pa
       */
      neuDataX<<5e3*(1-load_fac),0,0;
      neuDataY<<0,-5e3*(1-load_fac),0;
      neuX.setValue(neuDataX,3);
      neuY.setValue(neuDataY,3);

      assembler->updateBCs(BCs);

      gsDebugVar(neuDataX);
      if (nonfollow)
      {
        forceVec<<0,0,pressure*load_fac;
        force.setValue(forceVec,3);
      }
      else
      {

        pressFun.setValue(pressure*load_fac,3);
        assembler->setPressure(pressFun);
      }


      gsInfo<<"Load step "<<k<<"; p = "<<pressure*load_fac<<"; load factor = "<<load_fac<<":\n";
      real_t tol = 1e-4;
      real_t resNorm = 1;

      assembler->assemble();
      K = assembler->matrix();
      gsVector<> F = assembler->rhs();
      assembler->assembleMass(true);
      gsVector<> M = assembler->rhs();

      bool bisected = false;

      maxIt = 1e4;
      gsStaticDR<real_t> DRM(M,F,Residual);
      gsOptionList DROptions = DRM.options();
      DROptions.setReal("damping",damping);
      DROptions.setReal("alpha",alpha);
      DROptions.setInt("maxIt",maxIt);
      DROptions.setReal("tol",1e-1);
      DROptions.setReal("tolE",1e-1);
      DROptions.setInt("verbose",verbose);
      DROptions.setInt("ResetIt",(index_t)(0.1*maxIt));
      DRM.setDisplacement(solVector);
      DRM.setOptions(DROptions);
      DRM.initialize();
      DRM.solve();
      if (!DRM.converged())
      {
        gsWarn<<"Load step "<<k<<" did not converge\n";
        GISMO_ASSERT(load_fac!=0,"load_fac is zero but no convergence on the first step. Try to increase the number of iterations");
        load_fac -= dload_fac;
        dload_fac /= 2;
        load_fac += dload_fac;
        bisected = true;
        continue;
      }

      gsVector<> updateVector = DRM.solution() - solVector;

      maxIt = 50;
      // gsVector<> updateVector(assembler->numDofs());
      // updateVector.setZero();
      gsStaticNewton<real_t> NWT(K,F,Jacobian,Residual);
      gsOptionList NWTOptions = NWT.options();
      NWTOptions.setInt("maxIt",maxIt);
      NWTOptions.setReal("tol",tol);
      NWTOptions.setInt("verbose",verbose);
      NWT.setOptions(NWTOptions);
      NWT.reset();
      NWT.setDisplacement(solVector);
      NWT.setUpdate(updateVector);
      NWT.solve();
      if (!NWT.converged())
      {
        gsWarn<<"Load step "<<k<<" did not converge\n";
        GISMO_ASSERT(load_fac!=0,"load_fac is zero but no convergence on the first step. Try to increase the number of iterations");
        load_fac -= dload_fac;
        dload_fac /= 2;
        load_fac += dload_fac;
        bisected = true;
        continue;
      }
      solVector = NWT.solution();

      mp_def = assembler->constructSolution(solVector);

      gsMultiPatch<> deformation = mp_def;
      for (size_t k = 0; k != mp_def.nPatches(); ++k)
          deformation.patch(k).coefs() -= mp.patch(k).coefs();


      // ! [Export visualization in ParaView]
      if (plot)
      {
          gsField<> solField(mp_def, deformation);
          gsInfo<<"Plotting in Paraview...\n";
          gsWriteParaview<>( solField, "solution", 1000, true);
          // ev.options().setSwitch("plot.elements", true);
          // ev.writeParaview( u_sol   , G, "solution");

          // gsFileManager::open("solution.pvd");


          assembler->constructSolution(solVector,mp_def);
          gsPiecewiseFunction<> TFes;
          assembler->constructStress(mp_def,TFes,stress_type::tension_field);
          gsField<> TF(mp_def,TFes, true);
          gsWriteParaview(TF,"tensionfield",5000);


          gsInfo <<"Maximum deformation coef: "
                 << deformation.patch(0).coefs().colwise().maxCoeff() <<".\n";
          gsInfo <<"Minimum deformation coef: "
                 << deformation.patch(0).coefs().colwise().minCoeff() <<".\n";
      }


      // for (index_t i=0; i!=10; i++)
      // {
      //   gsVector<> R = Residual(solVector);
      //   gsSparseMatrix<> K_NL = Jacobian(solVector);

      //   // if (k < 3)
      //   //   K_NL += diag;
      //   // gsDebugVar(K_NL.toDense());
      //   // gsDebugVar(R.transpose());

      //   solver.compute(K_NL);
      //   gsVector<> updateVector = solver.solve(R);
      //   // gsDebugVar(updateVector);
      //   // gsDebugVar(Fneu.norm());

      //   solVector+= updateVector;
      //   // gsDebugVar(solVector);
      //   resNorm = math::abs(R.norm()/Fnorm);
      //   gsInfo<<"\titeration "<<i<<": resNorm = "<<resNorm<<"; Unorm = "<<solVector.norm()<<"; dUnorm = "<<updateVector.norm()<<"\n";

      //   if (resNorm < tol)
      //     break;
      // }
      // if (resNorm > tol)
      // {
      //   gsWarn<<"Load step "<<k<<" did not converge\n";
      //   GISMO_ASSERT(load_fac!=0,"load_fac is zero but no convergence on the first step. Try to increase the number of iterations");
      //   load_fac -= dload_fac;
      //   dload_fac /= 2;
      //   load_fac += dload_fac;
      //   continue;
      // }

      if (bisected)
        bisected = false;
      else
        dload_fac = dload_fac0;
      
      k++;
      if (1 - load_fac < dload_fac && 1 - load_fac > 0)
        dload_fac = 1 - load_fac;
      load_fac += dload_fac;
    }



    // gsStaticDR<real_t> DRM(M,F,Residual);
    // gsOptionList DROptions = DRM.options();
    // DROptions.setReal("damping",damping);
    // DROptions.setReal("alpha",alpha);
    // DROptions.setInt("maxIt",maxIt);
    // DROptions.setReal("tol",1e-1);
    // DROptions.setReal("tolE",1e-1);
    // DROptions.setInt("verbose",verbose);
    // DRM.setOptions(DROptions);
    // DRM.initialize();
    // DRM.solve();

    // maxIt = 100;
    // gsStaticNewton<real_t> NWT(K,F,Jacobian,Residual);
    // gsOptionList NWTOptions = NWT.options();
    // NWTOptions.setInt("maxIt",maxIt);
    // NWTOptions.setReal("tol",1e-6);
    // NWTOptions.setInt("verbose",verbose);
    // NWT.setOptions(NWTOptions);
    // NWT.reset();
    // NWT.setUpdate(DRM.update());
    // gsDebugVar(DRM.update().norm());
    // gsDebugVar(NWT.update().norm());
    // NWT.solve();
    // gsMatrix<> solVector = NWT.solution();

    mp_def = assembler->constructSolution(solVector);

    gsMultiPatch<> deformation = mp_def;
    for (size_t k = 0; k != mp_def.nPatches(); ++k)
        deformation.patch(k).coefs() -= mp.patch(k).coefs();


    // ! [Export visualization in ParaView]
    if (plot)
    {
        gsField<> solField(mp_def, deformation);
        gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview<>( solField, "solution", 1000, true);
        // ev.options().setSwitch("plot.elements", true);
        // ev.writeParaview( u_sol   , G, "solution");

        // gsFileManager::open("solution.pvd");


        assembler->constructSolution(solVector,mp_def);
        gsPiecewiseFunction<> TFes;
        assembler->constructStress(mp_def,TFes,stress_type::tension_field);
        gsField<> TF(mp_def,TFes, true);
        gsWriteParaview(TF,"tensionfield",5000);


        gsInfo <<"Maximum deformation coef: "
               << deformation.patch(0).coefs().colwise().maxCoeff() <<".\n";
        gsInfo <<"Minimum deformation coef: "
               << deformation.patch(0).coefs().colwise().minCoeff() <<".\n";
    }
    if (stress)
    {
      gsPiecewiseFunction<> membraneStresses;
      assembler->constructStress(mp_def,membraneStresses,stress_type::membrane);
      gsField<> membraneStress(mp_def,membraneStresses, true);

      gsPiecewiseFunction<> flexuralStresses;
      assembler->constructStress(mp_def,flexuralStresses,stress_type::flexural);
      gsField<> flexuralStress(mp_def,flexuralStresses, true);

      gsPiecewiseFunction<> stretches;
      assembler->constructStress(mp_def,stretches,stress_type::principal_stretch);
      gsField<> Stretches(mp_def,stretches, true);

      gsPiecewiseFunction<> pstrain_m;
      assembler->constructStress(mp_def,pstrain_m,stress_type::principal_membrane_strain);
      gsField<> pstrainM(mp_def,pstrain_m, true);

      gsPiecewiseFunction<> pstrain_f;
      assembler->constructStress(mp_def,pstrain_f,stress_type::principal_flexural_strain);
      gsField<> pstrainF(mp_def,pstrain_f, true);

      gsPiecewiseFunction<> pstress_m;
      assembler->constructStress(mp_def,pstress_m,stress_type::principal_stress_membrane);
      gsField<> pstressM(mp_def,pstress_m, true);

      gsPiecewiseFunction<> pstress_f;
      assembler->constructStress(mp_def,pstress_f,stress_type::principal_stress_flexural);
      gsField<> pstressF(mp_def,pstress_f, true);

      gsPiecewiseFunction<> stretch1;
      assembler->constructStress(mp_def,stretch1,stress_type::principal_stretch_dir1);
      gsField<> stretchDir1(mp_def,stretch1, true);

      gsPiecewiseFunction<> stretch2;
      assembler->constructStress(mp_def,stretch2,stress_type::principal_stretch_dir2);
      gsField<> stretchDir2(mp_def,stretch2, true);

      gsPiecewiseFunction<> stretch3;
      assembler->constructStress(mp_def,stretch3,stress_type::principal_stretch_dir3);
      gsField<> stretchDir3(mp_def,stretch3, true);

      gsPiecewiseFunction<> VMStresses;
      assembler->constructStress(mp_def,VMStresses,stress_type::von_mises_membrane);
      gsField<> VMStress(mp_def,VMStresses, true);


      gsPiecewiseFunction<> TFes;
      assembler->constructStress(mp_def,TFes,stress_type::tension_field);
      gsField<> TF(mp_def,TFes, true);


      gsWriteParaview(membraneStress,"MembraneStress",5000);
      gsWriteParaview(VMStress,"MembraneStressVM",5000);
      gsWriteParaview(Stretches,"PrincipalStretch",5000);
      gsWriteParaview(pstrainM,"PrincipalMembraneStrain",5000);
      gsWriteParaview(pstrainF,"PrincipalFlexuralStrain",5000);
      gsWriteParaview(pstressM,"PrincipalMembraneStress",5000);
      gsWriteParaview(pstressF,"PrincipalFlexuralStress",5000);
      gsWriteParaview(stretchDir1,"PrincipalDirection1",5000);
      gsWriteParaview(stretchDir1,"PrincipalDirection1",5000);
      gsWriteParaview(stretchDir2,"PrincipalDirection2",5000);
      gsWriteParaview(stretchDir3,"PrincipalDirection3",5000);
      gsWriteParaview(TF,"tensionfield",5000);
    }
    gsInfo<<"Total ellapsed assembly time: \t\t"<<time<<" s\n";

    delete materialMatrix;
    delete materialMatrixTFT;
    delete assembler;

    return EXIT_SUCCESS;
}
