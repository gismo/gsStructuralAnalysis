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
#include <gsKLShell/gsMaterialMatrixEval.h>

#include <gsStructuralAnalysis/gsStaticDR.h>
#include <gsStructuralAnalysis/gsStaticNewton.h>
#include <gsStructuralAnalysis/gsStructuralAnalysisUtils.h>

using namespace gismo;

int main (int argc, char** argv)
{
    // Input options
    int numElevate    = 0;
    int numRefine     = 1;
    bool plot         = false;
    bool write        = false;
    bool stress       = false;
    bool quasiNewton  = false;
    int quasiNewtonInt= -1;
    bool adaptive     = false;
    int method        = 2; // (0: Load control; 1: Riks' method; 2: Crisfield's method; 3: consistent crisfield method; 4: extended iterations)
    int verbose = 0;

    bool membrane = false;
    bool nonfollow     = false;

    index_t Compressibility = 0;
    index_t material = 0;
    index_t impl = 1; // 1= analytical, 2= generalized, 3= spectral
    bool TFT = false;

    int result        = 0;
    index_t maxIt     = 1e3;
    // Arc length method options
    real_t tol        = 1e-3;

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

    cmd.addInt( "M", "Material", "Material law",  material );
    cmd.addInt( "C", "Compressibility", "1: compressible, 0: incompressible",  Compressibility );
    cmd.addInt( "I", "Implementation", "Implementation: 1= analytical, 2= generalized, 3= spectral",  impl );
    cmd.addSwitch( "TFT", "Use Tension-Field Theory",  TFT );

    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addSwitch("write", "Write data to file", write);
    cmd.addSwitch("stress", "Plot stress in ParaView format", stress);
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

    if (material==0) Compressibility = 0;

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

    std::string dirname = ".";
    dirname = dirname + "/Pillow_r" + std::to_string(numRefine) + "_e" + std::to_string(numElevate) + "-M" + std::to_string(material) + "-c" + std::to_string(Compressibility);
    if (TFT)
      dirname = dirname + "_TFT";
    std::string output =  "solution";

    gsFileManager::mkdir(dirname);

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

    gsMaterialMatrixBase<real_t>* materialMatrix;
    gsMaterialMatrixBase<real_t>* materialMatrixTFT;
    gsOptionList options;

    std::vector<gsFunction<>*> parameters;
    parameters.resize(2);
    parameters[0] = &E;
    parameters[1] = &nu;
    if      (material==0 && impl==1)
    {
      parameters.resize(2);
      options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
      options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
      materialMatrix = getMaterialMatrix<3,real_t>(mp,t,parameters,rho,options);
      materialMatrixTFT = new gsMaterialMatrixTFT<3,real_t,true>(static_cast<gsMaterialMatrixBaseDim<3,real_t> * >(materialMatrix));
      // dynamic_cast<gsMaterialMatrixTFT<3,real_t,true> *>(materialMatrixTFT)->updateDeformed(&mp_def);
    }
    else if (material==1 || material==2)
    {
      options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",material);
      options.addSwitch("Compressibility","Compressibility: (false): Imcompressible | (true): Compressible",Compressibility);
      options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",impl);
      materialMatrix = getMaterialMatrix<3,real_t>(mp,t,parameters,rho,options);
      materialMatrixTFT = new gsMaterialMatrixTFT<3,real_t,false>(static_cast<gsMaterialMatrixBaseDim<3,real_t> * >(materialMatrix));
      // dynamic_cast<gsMaterialMatrixTFT<3,real_t,false> *>(materialMatrixTFT)->updateDeformed(&mp_def);
    }
    else
      GISMO_ERROR("Other materials not implemented in this example");

    // gsMaterialMatrixBase<real_t>* materialMatrix = new gsMaterialMatrixLinear<3,real_t>(mp,t,E,nu,rho);
    // gsMaterialMatrixTFT<3,real_t,false> * materialMatrixTFT = new gsMaterialMatrixTFT<3,real_t,false>(static_cast<gsMaterialMatrixBaseDim<3,real_t> * >(materialMatrix));

    materialMatrixTFT->options().setReal("SlackMultiplier",1e-6);    
    // materialMatrixTFT->options().setSwitch("Explicit",true);    

    gsThinShellAssemblerBase<real_t>* assembler;
    if (membrane && TFT)
      assembler = new gsThinShellAssembler<3, real_t, true >(mp,dbasis,BCs,force,materialMatrixTFT);
    else if (membrane && !TFT)
      assembler = new gsThinShellAssembler<3, real_t, true >(mp,dbasis,BCs,force,materialMatrix);
    else if (!membrane && TFT)
      assembler = new gsThinShellAssembler<3, real_t, false >(mp,dbasis,BCs,force,materialMatrixTFT);
    else
      assembler = new gsThinShellAssembler<3, real_t, false >(mp,dbasis,BCs,force,materialMatrix);

    // Assemble linear system to obtain the force vector
    assembler->assemble();
    gsSparseMatrix<> K = assembler->matrix();
    gsVector<> Fneu = assembler->rhs();
    real_t Fnorm = Fneu.norm();
    gsDebugVar(assembler->numDofs());

    // gsSparseMatrix<> diag(assembler->numDofs(),assembler->numDofs());
    // diag.setIdentity();
    // diag *= 0.14;
    // Function for the Jacobian
    gsStructuralAnalysisOps<real_t>::Jacobian_t Jacobian = [&assembler,&mp_def](gsVector<real_t> const &x, gsSparseMatrix<real_t> & m)
    {
      // &diag,
      ThinShellAssemblerStatus status;
      assembler->constructSolution(x,mp_def);
      status = assembler->assembleMatrix(mp_def);
      m = assembler->matrix();
      // m += diag;
      return status == ThinShellAssemblerStatus::Success;
    };
    // Function for the Residual
    gsStructuralAnalysisOps<real_t>::Residual_t Residual = [&assembler,&mp_def](gsVector<real_t> const &x, gsVector<real_t> & result)
    {
      ThinShellAssemblerStatus status;
      assembler->constructSolution(x,mp_def);
      status = assembler->assembleVector(mp_def);
      result = assembler->rhs();
      return status == ThinShellAssemblerStatus::Success;
    };

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
    index_t step=0;
    bool bisected = false;

    gsMatrix<> writePoints(2,1);
    writePoints<<0,1;

    gsParaviewCollection TFcollection(dirname + "/" + "tensionfield");
    gsParaviewCollection solcollection(dirname + "/" + "solution");

    gsStructuralAnalysisOutput<real_t> writer(dirname + "/" + wn,writePoints);
    std::vector<std::string> pointheaders = {"x","y","z"};
    std::vector<std::string> otherheaders = {"load"};
    writer.init(pointheaders,otherheaders);
    while (load_fac <= 1.)
    {
      /*
       * Load stepping:
       * - First 3-5 steps: K+0.14*I
       * - First 20 steps: tension on free boundaries of 4e3 N/m decreasing to 0 N/m
       * - First 20 steps: follower pressure increasing from 0 Pa to 5000 Pa
       */
      if (true)
      {
        neuDataX<<5e3*(1-load_fac),0,0;
        neuDataY<<0,-5e3*(1-load_fac),0;
        neuX.setValue(neuDataX,3);
        neuY.setValue(neuDataY,3);        
      }
      else 
      {
        neuDataX<<0,0,0;
        neuDataY<<0,0,0;
        neuX.setValue(neuDataX,3);
        neuY.setValue(neuDataY,3);
      }


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

      gsInfo<<"Load step "<<step<<"; p = "<<pressure*load_fac<<"; load factor = "<<load_fac<<":\n";
      real_t resNorm = 1;

      // if (load_fac > 0.5)
      //   diag.setZero();

      assembler->assemble();
      K = assembler->matrix();
      // K += diag;
      gsVector<> F = assembler->rhs();
      assembler->assembleMass(true);
      gsVector<> M = assembler->rhs();

      maxIt = 1e3;
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
      // if (!DRM.converged())
      // {
      //   gsWarn<<"Load step "<<step<<" did not converge\n";
      //   GISMO_ASSERT(load_fac!=0,"load_fac is zero but no convergence on the first step. Try to increase the number of iterations");
      //   load_fac -= dload_fac;
      //   dload_fac /= 2;
      //   load_fac += dload_fac;
      //   bisected = true;
      //   continue;
      // }

      gsVector<> updateVector = DRM.solution() - solVector;

      maxIt = 1000;
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
        gsWarn<<"Load step "<<step<<" did not converge\n";
        GISMO_ASSERT(load_fac!=0,"load_fac is zero but no convergence on the first step. Try to increase the number of iterations");
        load_fac -= dload_fac;
        dload_fac /= 2;
        load_fac += dload_fac;
        bisected = true;
        continue;
      }
      solVector = NWT.solution();

      mp_def = assembler->constructSolution(solVector);
      // if (material==0)
      //   dynamic_cast<gsMaterialMatrixTFT<3,real_t,true> *>(materialMatrixTFT)->updateDeformed(&mp_def);
      // else if (material==1 || material==2)
      //   dynamic_cast<gsMaterialMatrixTFT<3,real_t,false> *>(materialMatrixTFT)->updateDeformed(&mp_def);
      // else GISMO_ERROR("Material unknown");

      gsMultiPatch<> deformation = mp_def;
      for (size_t k = 0; k != mp_def.nPatches(); ++k)
          deformation.patch(k).coefs() -= mp.patch(k).coefs();

      // ! [Export visualization in ParaView]
      if (plot)
      {
          std::string fileName;
          gsField<> solField(mp_def, deformation);
          gsInfo<<"Plotting in Paraview...\n";
          
          fileName = dirname + "/" + "solution" + util::to_string(step);
          gsWriteParaview<>( solField, fileName, 1000, true);
          solcollection.addPart(gsFileManager::getFilename(fileName) + "0.vts",step);

          // ev.options().setSwitch("plot.elements", true);
          // ev.writeParaview( u_sol   , G, "solution");

          // gsFileManager::open("solution.pvd");


          assembler->constructSolution(solVector,mp_def);
          gsPiecewiseFunction<> TFes;
          assembler->constructStress(mp_def,TFes,stress_type::tension_field);
          gsField<> TF(mp_def,TFes, true);

          fileName = dirname + "/" + "tensionfield" + util::to_string(step);
          gsWriteParaview(TF, fileName,5000);
          TFcollection.addPart(gsFileManager::getFilename(fileName) + "0.vts",step);


          gsInfo <<"Maximum deformation coef: "
                 << deformation.patch(0).coefs().colwise().maxCoeff() <<".\n";
          gsInfo <<"Minimum deformation coef: "
                 << deformation.patch(0).coefs().colwise().minCoeff() <<".\n";
      }
      if (write)
      {
        gsMatrix<> result = deformation.patch(0).eval(writePoints);
        gsVector<> data(1);
        data<<load_fac;
        writer.add(result,data);
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
      {
        bisected = false;
        dload_fac *= 2;
      }
      else
        dload_fac = dload_fac0;
      
      step++;
      if (1 - load_fac < dload_fac && 1 - load_fac > 1e-10)
        dload_fac = 1 - load_fac;
      load_fac += dload_fac;
    }

    solcollection.save();
    TFcollection.save();


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


      gsWriteParaview(membraneStress,dirname + "/" + "MembraneStress",5000);
      gsWriteParaview(VMStress,dirname + "/" + "MembraneStressVM",5000);
      gsWriteParaview(Stretches,dirname + "/" + "PrincipalStretch",5000);
      gsWriteParaview(pstrainM,dirname + "/" + "PrincipalMembraneStrain",5000);
      gsWriteParaview(pstrainF,dirname + "/" + "PrincipalFlexuralStrain",5000);
      gsWriteParaview(pstressM,dirname + "/" + "PrincipalMembraneStress",5000);
      gsWriteParaview(pstressF,dirname + "/" + "PrincipalFlexuralStress",5000);
      gsWriteParaview(stretchDir1,dirname + "/" + "PrincipalDirection1",5000);
      gsWriteParaview(stretchDir1,dirname + "/" + "PrincipalDirection1",5000);
      gsWriteParaview(stretchDir2,dirname + "/" + "PrincipalDirection2",5000);
      gsWriteParaview(stretchDir3,dirname + "/" + "PrincipalDirection3",5000);
      
    }

    delete materialMatrix;
    delete materialMatrixTFT;
    delete assembler;

    return EXIT_SUCCESS;
}
