/** @file quasistatic_shell_XML.cpp

    @brief Blackbox solver for quasistatic shell analysis

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#include <gismo.h>

#ifdef gsKLShell_ENABLED
#include <gsKLShell/src/gsThinShellAssembler.h>
#include <gsKLShell/src/gsMaterialMatrixTFT.h>
#include <gsKLShell/src/gsFunctionSum.h>
#endif

#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicBase.h>
#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicExplicitEuler.h>
#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicImplicitEuler.h>
#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicRK4.h>
#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicNewmark.h>
#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicBathe.h>
#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicWilson.h>

#include <gsStructuralAnalysis/src/gsStructuralAnalysisTools/gsStructuralAnalysisUtils.h>

#ifdef gsUnstructuredSplines_ENABLED
#include <gsUnstructuredSplines/src/gsSmoothInterfaces.h>
#include <gsUnstructuredSplines/src/gsAlmostC1.h>
#include <gsUnstructuredSplines/src/gsDPatch.h>
#include <gsUnstructuredSplines/src/gsC1SurfSpline.h>
#endif


using namespace gismo;

#ifdef gsKLShell_ENABLED
int main (int argc, char** argv)
{
    // Input options
    int numElevate    = 0;
    int numRefine     = 1;
    bool plot         = false;
    bool write        = false;
    bool remap        = false;

    int step = 10;

    real_t dt = 1e-2;

    index_t method = 0;

    std::string bvp;

    std::string dirname = ".";

    std::string wn("data.csv");

    gsCmdLine cmd("Shell quasistatic solver.");

    cmd.addInt("r","hRefine", "Number of dyadic h-refinement (bisection) steps to perform before solving", numRefine);
    cmd.addInt("e","degreeElevation", "Number of degree elevation steps to perform on the Geometry's basis before solving", numElevate);
    cmd.addInt("m","method", "Smoothing method to use: 0: smoothInterfaces, 1: Almost C1, 2: D-Patch", method);

    cmd.addReal("d","dt","Time Step",dt);

    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addSwitch("write", "Write data to file", write);
    cmd.addSwitch("remap", "Remap the geometry", remap);

    cmd.addString("i","inputFile", "Input file", bvp);
    cmd.addString("o","outputDir", "Output directory", dirname);

    cmd.addInt("N", "maxsteps", "solMaximum number of steps", step);

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

    gsOptionList DynamicOptions, assemblerOptions;
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

    gsInfo<<"Reading assembler options (ID=92) ...";
    if ( fd.hasId(92) ) fd.getId(92,assemblerOptions);
    gsInfo<<"Finished\n";
    gsInfo<<"Reading ALM solver options (ID=96) ...";
    if ( fd.hasId(96) ) fd.getId(96,DynamicOptions);
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

    // Make unstructured spline
    gsMultiPatch<> geom;
    gsMappedBasis<2,real_t> bb2;
    gsSparseMatrix<> global2local;
    if (method==0)
    {
      gsSmoothInterfaces<2,real_t> smoothInterfaces(mp);
      smoothInterfaces.options().setSwitch("SharpCorners",false);
      smoothInterfaces.compute();
      smoothInterfaces.matrix_into(global2local);

      global2local = global2local.transpose();
      geom = smoothInterfaces.exportToPatches();
      dbasis = smoothInterfaces.localBasis();
      bb2.init(dbasis,global2local);
    }
    else if (method==1)
    {
      gsAlmostC1<2,real_t> almostC1(mp);
      almostC1.options().setSwitch("SharpCorners",true);
      almostC1.compute();
      almostC1.matrix_into(global2local);

      global2local = global2local.transpose();
      geom = almostC1.exportToPatches();
      dbasis = almostC1.localBasis();
      bb2.init(dbasis,global2local);

      if (remap)
      {
        gsDofMapper mapper(dbasis);
        mapper.finalize();
        gsMatrix<> coefs;
        // First project the geometry mp onto bb2 and make a mapped spline
        gsInfo<<"L2-Projection error of mp on bb2 = "<<gsL2Projection<real_t>::projectGeometry(dbasis,bb2,mp,coefs)<<"\n";
        coefs.resize(coefs.rows()/mp.geoDim(),mp.geoDim());
        gsMappedSpline<2,real_t> mspline;
        mspline.init(bb2,coefs);
        if (plot) gsWriteParaview( mspline, "mspline");

        // Then project onto dbasis so that geom represents the mapped geometry
        gsInfo<<"L2-Projection error of mspline on dbasis = "<<gsL2Projection<real_t>::projectGeometry(dbasis,mspline,coefs)<<"\n";
        coefs.resize(coefs.rows()/mp.geoDim(),mp.geoDim());

        index_t offset = 0;
        for (size_t p = 0; p != geom.nPatches(); p++)
        {
            geom.patch(p) = give(*dbasis.basis(p).makeGeometry((coefs.block(offset,0,mapper.patchSize(p),mp.geoDim()))));
            offset += mapper.patchSize(p);
        }
      }
    }
    else if (method==2)
    {
      gsDPatch<2,real_t> dpatch(mp);
      dpatch.options().setSwitch("SharpCorners",true);
      dpatch.compute();
      dpatch.matrix_into(global2local);

      global2local = global2local.transpose();
      geom = dpatch.exportToPatches();
      dbasis = dpatch.localBasis();
      bb2.init(dbasis,global2local);

      if (remap)
          for (size_t p=0; p!=geom.nPatches(); p++)
          {
              gsMatrix<> coefs;
              gsQuasiInterpolate<real_t>::localIntpl(dbasis.basis(p),mp.patch(p),coefs);
              geom.patch(p).coefs() = coefs;
              gsInfo<<"basis "<<p<<":\n"<<dbasis.basis(0)<<"\n";
          }
      gsWriteParaview(geom,"geom");
    }
    else if (method==3)
    {
      dbasis = gsMultiBasis<>(mp);
      gsC1SurfSpline<2,real_t> smoothC1(mp,dbasis);
      smoothC1.init();
      smoothC1.compute();

      global2local = smoothC1.getSystem();
      global2local = global2local.transpose();
      smoothC1.getMultiBasis(dbasis);
      bb2.init(dbasis,global2local);
      geom = mp;
    }
    else
      GISMO_ERROR("Method "<<method<<" unknown");

    // Finalize BCs
    BCs.setGeoMap(geom);

    gsThinShellAssemblerBase<real_t>* assembler;
    assembler = new gsThinShellAssembler<3, real_t, true  >(geom,dbasis,BCs,forceFun,materialMatrixContainer);
    assembler->setOptions(assemblerOptions);
    assembler->options().setInt("Continuity",-1);
    assembler->setSpaceBasis(bb2);
    assembler->setPointLoads(pLoads);
    if (pressure)
      assembler->setPressure(pressFun);

    // Assemble linear system to obtain the force vector
    assembler->assemble();
    gsVector<> Force = assembler->rhs();
    // Compute mass matrix (since it is constant over time)
    assembler->assembleMass();
    gsSparseMatrix<> M = assembler->massMatrix();

    // Function for the Jacobian
    gsStructuralAnalysisOps<real_t>::Jacobian_t Jacobian = [&assembler,&bb2,&geom](gsVector<real_t> const &x, gsSparseMatrix<real_t> & m)
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
    gsStructuralAnalysisOps<real_t>::TResidual_t TResidual = [&assembler,&bb2,&geom](gsVector<real_t> const &x, real_t time, gsVector<real_t> & result)
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

    gsSparseMatrix<> C = gsSparseMatrix<>(assembler->numDofs(),assembler->numDofs());
    gsStructuralAnalysisOps<real_t>::Damping_t Damping = [&C](const gsVector<real_t> &, gsSparseMatrix<real_t> & m) { m = C; return true; };
    gsStructuralAnalysisOps<real_t>::Mass_t    Mass    = [&M](                          gsSparseMatrix<real_t> & m) { m = M; return true; };

    gsStructuralAnalysisOutput<real_t> writer(dirname + sep + wn,refPoints);
    std::vector<std::string> pointheaders = {"x","y","z"};
    std::vector<std::string> otherheaders = {"time"};
    writer.init(pointheaders,otherheaders);

    gsDynamicBase<real_t> * dynamic;
    std::string timeMethod = DynamicOptions.askString("Method","Newmark");
    if      (timeMethod=="ExplicitEuler")
      dynamic = new gsDynamicExplicitEuler<real_t,true>(Mass,Damping,Jacobian,TResidual);
    else if (timeMethod=="ExplicitEuler")
      dynamic = new gsDynamicImplicitEuler<real_t,true>(Mass,Damping,Jacobian,TResidual);
    else if (timeMethod=="RK4")
      dynamic = new gsDynamicRK4<real_t,true>(Mass,Damping,Jacobian,TResidual);
    else if (timeMethod=="Newmark")
      dynamic = new gsDynamicNewmark<real_t,true>(Mass,Damping,Jacobian,TResidual);
    else if (timeMethod=="Bathe")
      dynamic = new gsDynamicBathe<real_t,true>(Mass,Damping,Jacobian,TResidual);
    else if (timeMethod=="Wilson")
      dynamic = new gsDynamicWilson<real_t,true>(Mass,Damping,Jacobian,TResidual);
    else
      GISMO_ERROR("Method "<<method<<" unknown");

    dynamic->options().setReal("DT",dt);
    dynamic->options().setReal("TolU",1e-3);
    dynamic->options().setSwitch("Verbose",true);

    size_t N = assembler->numDofs();
    gsVector<> U(N), V(N), A(N);
    U.setZero();
    V.setZero();
    A.setZero();

    gsParaviewCollection collection(dirname + sep + output);
    real_t time = 0;
    gsMultiPatch<> deformation = mp;

    for (index_t i=0; i<step; i++)
    {
      gsInfo<<"Time "<<time<<"\n";
      gsStatus status = dynamic->step(time,dt,U,V,A);
      GISMO_ASSERT(status == gsStatus::Success,"Time integrator did not succeed");
      time += dt;

      gsVector<> solVector = U;
      if (plot||write)
      {
        gsField<> solField;
        std::string fileName = dirname + sep + output + util::to_string(i);

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

        // Update the displacement vector
        gsMultiPatch<> solution = assembler->constructDisplacement(U);
        if (plot)
        {
          // 5. Plot the mapped spline on the original geometry
          gsPiecewiseFunction<> displacements;
          assembler->constructStress(def,displacements,stress_type::displacement);
          gsWriteParaview(def,displacements,fileName,1000,"_");

          for (size_t p = 0; p!=mp.nPatches(); p++)
          {
            fileName = output + util::to_string(i) + "_" + util::to_string(p);
            collection.addPart(fileName + ".vts",i);
          }
        }
        if (write)
        {
          gsMatrix<> result(mp.geoDim(),refPoints.cols());
          for (index_t k=0; k!=refPoints.cols(); k++)
            result.col(k) = def.piece(refPatches.at(k)).eval(refPoints.col(k));

          gsMatrix<> otherDataMatrix(1,1);
          otherDataMatrix<<time;
          writer.add(result,otherDataMatrix);
        }
      }


    }

    if (plot)
      collection.save();

    delete assembler;
    delete dynamic;

    return EXIT_SUCCESS;
}
#else//gsKLShell_ENABLED
int main(int argc, char *argv[])
{
    gsWarn<<"G+Smo is not compiled with the gsKLShell module.";
    return EXIT_FAILURE;
}
#endif
