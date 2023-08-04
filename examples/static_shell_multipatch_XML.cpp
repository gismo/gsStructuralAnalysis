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
#include <gsKLShell/gsFunctionSum.h>

#include <gsStructuralAnalysis/gsStaticDR.h>
#include <gsStructuralAnalysis/gsStaticNewton.h>
#include <gsStructuralAnalysis/gsStaticComposite.h>

#include <gsStructuralAnalysis/gsStructuralAnalysisUtils.h>

#include <gsUnstructuredSplines/src/gsSmoothInterfaces.h>

using namespace gismo;

int main (int argc, char** argv)
{
    // Input options
    int numElevate    = 0;
    int numRefine     = 1;
    bool plot         = false;
    bool write        = false;
    bool stress       = false;
    bool DR  = false;
    bool NR  = false;

    bool membrane = false;

    real_t perturb = 0;

    // Arc length method options

    std::string bvp;

    std::string dirname = ".";

    std::string wn("data.csv");

    gsCmdLine cmd("Example for an inflating balloon.");

    cmd.addInt("r","hRefine", "Number of dyadic h-refinement (bisection) steps to perform before solving", numRefine);
    cmd.addInt("e","degreeElevation", "Number of degree elevation steps to perform on the Geometry's basis before solving", numElevate);

    cmd.addReal( "P" , "perturb" , "Perturb the z coordinate of the refined geometry",  perturb );

    cmd.addSwitch("membrane", "Membrane element", membrane);

    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addSwitch("write", "Write data to file", write);
    cmd.addSwitch("stress", "Plot stress in ParaView format", stress);
    cmd.addSwitch("DR", "Use Dynamic Relaxation", DR);
    cmd.addSwitch("NR", "Use Newton Raphson", NR);

    cmd.addString("i","inputFile", "Input file", bvp);
    cmd.addString("o","outputDir", "Output directory", dirname);


    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // ![Initialize members]
    gsMultiPatch<> mp,mp_def;
    // Boundary conditions
    gsBoundaryConditions<> BCs;

    gsFunctionExpr<> forceFun;
    gsFunctionExpr<> pressFun;

    gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();
    gsMatrix<> points,loads;
    gsMatrix<index_t> pid_ploads;

    gsOptionList DROptions, NROptions, assemblerOptions;
    // ![Initialize members]

    // ![Read data]
    gsFileData<real_t> fd(bvp);
    gsInfo<<"Reading geometry (ID=0) ...";
    fd.getId(0,mp);
    gsInfo<<"Finished\n";

    // p-refine
    for (size_t p=0; p!=mp.nPatches(); p++)
    {
      for(index_t i = 0; i< numElevate; ++i)
        mp.patch(p).degreeElevate();    // Elevate the degree

      // h-refine
      for(index_t i = 0; i< numRefine; ++i)
        mp.patch(p).uniformRefine();
    }

    if (perturb!=0)
      for (size_t p=0; p!=mp.nPatches(); p++)
      {
        mp.patch(p).coefs().col(2).setRandom();
        mp.patch(p).coefs().col(2) *= perturb;
      }

    gsMultiBasis<> dbasis(mp,true);
    gsInfo<<"Basis (patch 0): "<< mp.patch(0).basis() << "\n";
    mp_def = mp;


    gsInfo<<"Looking for material matrices ...\n";
    gsMaterialMatrixContainer<real_t> materialMatrixContainer;
    gsMaterialMatrixBase<real_t> * materialMatrix;
    if (fd.hasAny<gsMaterialMatrixContainer<real_t>>())
    {
      gsInfo<<"Reading material matrix container (ID=11) ...";
      fd.getId(11,materialMatrixContainer);
    }
    else
    {
      gsInfo<<"Reading material matrix (ID=10) ...";
      materialMatrix = fd.getId<gsMaterialMatrixBase<real_t>>(10).release();
      for (size_t p = 0; p!=mp.nPatches(); p++)
        materialMatrixContainer.add(materialMatrix);
    }
    gsInfo<<"Finished\n";
    gsInfo<<"Finished\n";

    gsInfo<<"Reading boundary conditions (ID=20) ...";
    fd.getId(20,BCs);
    gsInfo<<"Finished\n";

    // Finalize BCs
    BCs.setGeoMap(mp);

    gsInfo<<"Reading force function (ID=21) ...";
    fd.getId(21,forceFun);
    gsInfo<<"Finished\n";
    bool pressure = false;
    if ( fd.hasId(22) )
    {
      gsInfo<<"Reading pressure function (ID=22) ...";
      pressure = true;
      fd.getId(22,pressFun);
      gsInfo<<"Finished\n";
    }
    // Point loads
    gsInfo<<"Reading point load locations (ID=30) ...";
    if ( fd.hasId(30) ) fd.getId(30,points);
    gsInfo<<"Finished\n";
    gsInfo<<"Reading point load vectors (ID=31) ...";
    if ( fd.hasId(31) ) fd.getId(31,loads);
    gsInfo<<"Finished\n";
    gsInfo<<"Reading point load patch indices (ID=32) ...";
    if ( fd.hasId(32) ) fd.getId(32,pid_ploads);
    gsInfo<<"Finished\n";

    if ( !fd.hasId(30) || !fd.hasId(31) || !fd.hasId(32) )
        pid_ploads = gsMatrix<index_t>::Zero(1,points.cols());

    for (index_t k =0; k!=points.cols(); k++)
        pLoads.addLoad(points.col(k), loads.col(k), pid_ploads.at(k) ); // in parametric domain!


    // Reference points
    gsMatrix<index_t> refPatches;
    gsMatrix<> refPoints, refPars, refValue; // todo: add refValue..
    gsInfo<<"Reading reference point locations (ID=50) ...";
    if ( fd.hasId(50) ) fd.getId(50,refPoints);
    gsInfo<<"Finished\n";
    gsInfo<<"Reading reference patches (ID=51) ...";
    if ( fd.hasId(51) ) fd.getId(51,refPatches);
    gsInfo<<"Finished\n";
    gsInfo<<"Reading reference values (ID=52) ...";
    if ( fd.hasId(52) ) fd.getId(52,refValue);
    gsInfo<<"Finished\n";

    if ( !fd.hasId(50) || !fd.hasId(51) || !fd.hasId(52) )
        refValue = gsMatrix<>::Zero(mp.geoDim(),refPoints.cols());
    GISMO_ENSURE(refPatches.cols()==refPoints.cols(),"Number of reference points and patches do not match");

    if (refPoints.rows()==2)
    {
        refPars = refPoints;
        gsInfo<<"Reference points are provided in parametric coordinates.\n";
    }
    else if (refPoints.rows()==3)
        gsInfo<<"Reference points are provided in physical coordinates.\n";
    else
        gsInfo<<"No reference points are provided.\n";

    gsInfo<<"Reading DR solver options (ID=90) ...";
    if ( fd.hasId(90) ) fd.getId(90,DROptions);
    gsInfo<<"Finished\n";
    gsInfo<<"Reading NR solver options (ID=91) ...";
    if ( fd.hasId(91) ) fd.getId(91,NROptions);
    gsInfo<<"Finished\n";
    gsInfo<<"Reading assembler options (ID=92) ...";
    if ( fd.hasId(92) ) fd.getId(92,assemblerOptions);
    gsInfo<<"Finished\n";

    // ![Read data]
    std::string output =  "solution";

    // Fix path
    dirname = gsFileManager::getCanonicRepresentation(dirname,true);
    char sep = gsFileManager::getNativePathSeparator();
    gsFileManager::mkdir(dirname);

    // plot geometry
    if (plot)
      gsWriteParaview(mp,"mp",1000,true);

    // materialMatrixTFT->options().setSwitch("Explicit",true);


    mp.computeTopology();

    // Make unstructured spline
    gsMultiPatch<> geom;
    gsMappedBasis<2,real_t> bb2;
    gsSparseMatrix<> global2local;
    gsSmoothInterfaces<2,real_t> smoothInterfaces(mp);
    smoothInterfaces.options().setSwitch("SharpCorners",false);
    smoothInterfaces.compute();
    smoothInterfaces.matrix_into(global2local);

    global2local = global2local.transpose();
    geom = smoothInterfaces.exportToPatches();
    dbasis = smoothInterfaces.localBasis();
    bb2.init(dbasis,global2local);

    // Make assembler
    gsThinShellAssemblerBase<real_t>* assembler;
    if (membrane)
      assembler = new gsThinShellAssembler<3, real_t, false >(mp,dbasis,BCs,forceFun,materialMatrixContainer);
    else
      assembler = new gsThinShellAssembler<3, real_t, true  >(mp,dbasis,BCs,forceFun,materialMatrixContainer);
    assembler->setOptions(assemblerOptions);
    assembler->options().setInt("Continuity",-1);
    assembler->setSpaceBasis(bb2);
    assembler->setPointLoads(pLoads);
    if (pressure)
      assembler->setPressure(pressFun);

    // Assemble linear system to obtain the force vector
    assembler->assemble();
    gsSparseMatrix<> K = assembler->matrix();
    gsVector<> F = assembler->rhs();
    assembler->assembleMass(true);
    gsVector<> M = assembler->rhs();

    gsStructuralAnalysisOps<real_t>::Jacobian_t Jacobian;
    gsStructuralAnalysisOps<real_t>::Residual_t Residual;
    // Function for the Jacobian
    Jacobian = [&assembler,&bb2,&geom](gsVector<real_t> const &x, gsSparseMatrix<real_t> & m)
    {
      ThinShellAssemblerStatus status;
      gsMatrix<real_t> solFull = assembler->fullSolutionVector(x);
      size_t d = geom.targetDim();
      GISMO_ASSERT(solFull.rows() % d==0,"Rows of the solution vector does not match the number of control points");
      solFull.resize(solFull.rows()/d,d);
      gsMappedSpline<2,real_t> mspline(bb2,solFull);
      gsFunctionSum<real_t> def(&geom,&mspline);
      assembler->assembleMatrix(def);
      status = assembler->assembleMatrix(def);
      m = assembler->matrix();
      return status == ThinShellAssemblerStatus::Success;
    };
    // Function for the Residual
    Residual = [&assembler,&bb2,&geom](gsVector<real_t> const &x, gsVector<real_t> & result)
    {
      ThinShellAssemblerStatus status;
      gsMatrix<real_t> solFull = assembler->fullSolutionVector(x);
      size_t d = geom.targetDim();
      GISMO_ASSERT(solFull.rows() % d==0,"Rows of the solution vector does not match the number of control points");
      solFull.resize(solFull.rows()/d,d);

      gsMappedSpline<2,real_t> mspline(bb2,solFull);
      gsFunctionSum<real_t> def(&geom,&mspline);

      status = assembler->assembleVector(def);
      result = assembler->rhs(); // assembler rhs - force = Finternal
      return status == ThinShellAssemblerStatus::Success;
    };

    gsStructuralAnalysisOutput<real_t> writer(dirname + sep + wn,refPoints);
    std::vector<std::string> pointheaders = {"x","y","z"};
    std::vector<std::string> otherheaders = {};
    writer.init(pointheaders,otherheaders);

    gsStaticNewton<real_t> LIN(K,F);
    LIN.setOptions(NROptions);

    gsStaticDR<real_t> DRM(M,F,Residual);
    DRM.setOptions(DROptions);

    gsStaticNewton<real_t> NRM(K,F,Jacobian,Residual);
    NRM.setOptions(NROptions);

    std::vector<gsStaticBase<real_t> *> solvers{&LIN};
    if (DR)
    {
      gsInfo<<"Using Dynamic Relaxation solver\n";
      solvers.push_back(&DRM);
    }
    if (NR)
    {
      gsInfo<<"Using Newton-Raphson solver\n";
      solvers.push_back(&NRM);
    }
    gsStaticComposite<real_t> solver(solvers);
    solver.initialize();
    solver.solve();
    GISMO_ASSERT(solver.converged(),"Solver failed");
    gsVector<> solVector = solver.solution();

    // ! [Export visualization in ParaView]
    if (plot)
    {

      /// Make a gsMappedSpline to represent the solution
      // 1. Get all the coefficients (including the ones from the eliminated BCs.)
      gsMatrix<real_t> solFull = assembler->fullSolutionVector(solVector);

      // 2. Reshape all the coefficients to a Nx3 matrix
      size_t d = geom.targetDim();
      GISMO_ASSERT(solFull.rows() % d==0,"Rows of the solution vector does not match the number of control points");
      solFull.resize(solFull.rows()/d,d);

      // 3. Make the mapped spline
      gsMappedSpline<2,real_t> mspline(bb2,solFull);

      // 4. Plot the mapped spline on the original geometry
      gsField<> solField(geom, mspline,true);

      std::string fileName;
      gsInfo<<"Plotting in Paraview...\n";

      fileName = dirname + sep + "solution";
      gsWriteParaview<>( solField, fileName, 10000, true);

      assembler->constructSolution(solVector,mp_def);
      gsPiecewiseFunction<> TFes;
      assembler->constructStress(mp_def,TFes,stress_type::tension_field);
      gsField<> TF(mp_def,TFes, true);

      fileName = dirname + sep + "tensionfield";
      gsWriteParaview(TF, fileName,5000);
    }
    if (write)
    {
      // gsMatrix<> result(mp.geoDim(),refPoints.cols());
      // for (index_t k=0; k!=refPoints.cols(); k++)
      //   result.col(k) = deformation.patch(refPatches.at(k)).eval(refPoints.col(k));
      // writer.add(result);
    }

    // ! [Export visualization in ParaView]
    if (stress)
    {
      gsField<> tensionField;
      gsPiecewiseFunction<> tensionFields;
      std::string fileName = dirname + "/" + "tensionfield";

      /// Make a gsMappedSpline to represent the solution
      // 1. Get all the coefficients (including the ones from the eliminated BCs.)
      gsMatrix<real_t> solFull = assembler->fullSolutionVector(solVector);

      // 2. Reshape all the coefficients to a Nx3 matrix
      size_t d = geom.targetDim();
      GISMO_ASSERT(solFull.rows() % d==0,"Rows of the solution vector does not match the number of control points");
      solFull.resize(solFull.rows()/d,d);

      // 3. Make the mapped spline
      gsMappedSpline<2,real_t> mspline(bb2,solFull);

      // 4. Create deformation spline
      gsFunctionSum<real_t> def(&geom,&mspline);

      // 5. Construct stress
      assembler->constructStress(def,tensionFields,stress_type::tension_field);
      gsWriteParaview(def,tensionFields,fileName,20000,"_");


      // gsPiecewiseFunction<> membraneStresses;
      // assembler->constructStress(mp_def,membraneStresses,stress_type::membrane);
      // gsField<> membraneStress(mp_def,membraneStresses, true);

      // gsPiecewiseFunction<> flexuralStresses;
      // assembler->constructStress(mp_def,flexuralStresses,stress_type::flexural);
      // gsField<> flexuralStress(mp_def,flexuralStresses, true);

      // gsPiecewiseFunction<> stretches;
      // assembler->constructStress(mp_def,stretches,stress_type::principal_stretch);
      // gsField<> Stretches(mp_def,stretches, true);

      // gsPiecewiseFunction<> pstrain_m;
      // assembler->constructStress(mp_def,pstrain_m,stress_type::principal_membrane_strain);
      // gsField<> pstrainM(mp_def,pstrain_m, true);

      // gsPiecewiseFunction<> pstrain_f;
      // assembler->constructStress(mp_def,pstrain_f,stress_type::principal_flexural_strain);
      // gsField<> pstrainF(mp_def,pstrain_f, true);

      // gsPiecewiseFunction<> pstress_m;
      // assembler->constructStress(mp_def,pstress_m,stress_type::principal_stress_membrane);
      // gsField<> pstressM(mp_def,pstress_m, true);

      // gsPiecewiseFunction<> pstress_f;
      // assembler->constructStress(mp_def,pstress_f,stress_type::principal_stress_flexural);
      // gsField<> pstressF(mp_def,pstress_f, true);

      // gsPiecewiseFunction<> stretch1;
      // assembler->constructStress(mp_def,stretch1,stress_type::principal_stretch_dir1);
      // gsField<> stretchDir1(mp_def,stretch1, true);

      // gsPiecewiseFunction<> stretch2;
      // assembler->constructStress(mp_def,stretch2,stress_type::principal_stretch_dir2);
      // gsField<> stretchDir2(mp_def,stretch2, true);

      // gsPiecewiseFunction<> stretch3;
      // assembler->constructStress(mp_def,stretch3,stress_type::principal_stretch_dir3);
      // gsField<> stretchDir3(mp_def,stretch3, true);

      // gsPiecewiseFunction<> VMStresses;
      // assembler->constructStress(mp_def,VMStresses,stress_type::von_mises_membrane);
      // gsField<> VMStress(mp_def,VMStresses, true);


      // gsWriteParaview(membraneStress,dirname + sep + "MembraneStress",5000);
      // gsWriteParaview(VMStress,dirname + sep + "MembraneStressVM",5000);
      // gsWriteParaview(Stretches,dirname + sep + "PrincipalStretch",5000);
      // gsWriteParaview(pstrainM,dirname + sep + "PrincipalMembraneStrain",5000);
      // gsWriteParaview(pstrainF,dirname + sep + "PrincipalFlexuralStrain",5000);
      // gsWriteParaview(pstressM,dirname + sep + "PrincipalMembraneStress",5000);
      // gsWriteParaview(pstressF,dirname + sep + "PrincipalFlexuralStress",5000);
      // gsWriteParaview(stretchDir1,dirname + sep + "PrincipalDirection1",5000);
      // gsWriteParaview(stretchDir1,dirname + sep + "PrincipalDirection1",5000);
      // gsWriteParaview(stretchDir2,dirname + sep + "PrincipalDirection2",5000);
      // gsWriteParaview(stretchDir3,dirname + sep + "PrincipalDirection3",5000);
      
    }

    delete assembler;

    return EXIT_SUCCESS;
}
