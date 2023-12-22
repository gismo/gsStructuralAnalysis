/** @file benchmark_Balloon.cpp

    @brief Benchmark for the inflated pillow using the Arc-Length Method

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#include <gismo.h>

#include <gsKLShell/src/gsThinShellAssembler.h>
#include <gsKLShell/src/gsFunctionSum.h>

#include <gsStructuralAnalysis/src/gsEigenSolvers/gsBucklingSolver.h>

#include <gsStructuralAnalysis/src/gsStructuralAnalysisTools/gsStructuralAnalysisUtils.h>

#include <gsUnstructuredSplines/src/gsSmoothInterfaces.h>
#include <gsUnstructuredSplines/src/gsAlmostC1.h>
#include <gsUnstructuredSplines/src/gsDPatch.h>
#include <gsUnstructuredSplines/src/gsApproxC1Spline.h>
#include <gsUnstructuredSplines/src/gsC1SurfSpline.h>

using namespace gismo;

void writeToCSVfile(std::string name, gsMatrix<> matrix)
{
  std::ofstream file(name.c_str());
  for(int  i = 0; i < matrix.rows(); i++)
  {
    for(int j = 0; j < matrix.cols(); j++)
    {
       std::string str = std::to_string(matrix(i,j));
       if(j+1 == matrix.cols())
       {
           file<<std::setprecision(10)<<str;
       }
       else
       {
           file<<std::setprecision(10)<<str<<',';
       }
    }
    file<<'\n';
  }
}


int main (int argc, char** argv)
{
    // Input options
    int numElevate    = 0;
    int numRefine     = 1;
    bool plot         = false;
    bool write        = false;
    bool stress       = false;
    bool remap        = false;

    index_t degree = 3;
    index_t smoothness = 2;

    index_t method = 0;

    // Arc length method options
    std::string bvp;

    std::string dirname = "BucklingResults";

    index_t nmodes = 1;

    std::string wn("data.csv");

    gsCmdLine cmd("Example for an inflating balloon.");

    cmd.addInt("r","hRefine", "Number of dyadic h-refinement (bisection) steps to perform before solving", numRefine);
    cmd.addInt("e","degreeElevation", "Number of degree elevation steps to perform on the Geometry's basis before solving", numElevate);
    cmd.addInt("m","method", "Smoothing method to use: 0: smoothInterfaces, 1: Almost C1, 2: D-Patch", method);

    cmd.addInt( "p", "degree", "Set the polynomial degree of the basis.", degree );
    cmd.addInt( "s", "smoothness", "Set the smoothness of the basis.",  smoothness );

    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addSwitch("write", "Write data to file", write);
    cmd.addSwitch("stress", "Plot stress in ParaView format", stress);
    cmd.addSwitch("remap", "Remap the geometry", remap);

    cmd.addString("i","inputFile", "Input file", bvp);
    cmd.addString("o","outputDir", "Output directory", dirname);

    cmd.addInt( "N", "nmodes", "Number of modes",  nmodes );

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    GISMO_ENSURE(degree>smoothness,"Degree must be larger than the smoothness!");
    GISMO_ENSURE(smoothness>=0,"Degree must be larger than the smoothness!");
    if (method==3)
        GISMO_ENSURE(smoothness>=1 || smoothness <= degree-2,"Exact C1 method only works for smoothness <= p-2, but smoothness="<<smoothness<<" and p-2="<<degree-2);
    if (method==2 || method==3)
        GISMO_ENSURE(degree > 2,"Degree must be larger than 2 for the approx and exact C1 methods, but it is "<<degree);

    // ![Initialize members]
    gsMultiPatch<> mp,mp_def;
    // Boundary conditions
    gsBoundaryConditions<> BCs;

    gsFunctionExpr<> forceFun;
    gsFunctionExpr<> pressFun;

    gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();
    gsMatrix<> points,loads;
    gsMatrix<index_t> pid_ploads;

    gsOptionList BucklingOptions, assemblerOptions;
    // ![Initialize members]

    // ![Read data]
    gsFileData<real_t> fd(bvp);
    if (fd.hasAny<gsMultiPatch<real_t>>())
    {
      gsInfo<<"Reading geometry from "<<bvp<<" (ID=0) ...";
      fd.getId(0,mp);
    }
    else
    {
      std::string geomFileName = fd.getString(0);
      std::string path = gsFileManager::getPath(bvp) + gsFileManager::getNativePathSeparator();
      gsInfo<<"Reading geometry from "<<path + geomFileName<<" ...";
      gsReadFile<>(path + geomFileName,mp);
    }
    gsInfo<<"Finished\n";

    // p-refine
    for (size_t p=0; p!=mp.nPatches(); p++)
    {
      for(index_t i = 0; i< numElevate; ++i)
      {
        if (dynamic_cast<gsTensorNurbs<2,real_t> * >(&mp.patch(p)))
        {
          gsWarn<<"Degree elevation applied"<<"\n";
          mp.patch(p).degreeElevate();    // Elevate the degree
        }
        else
          mp.patch(p).degreeIncrease();    // Elevate the degree
      }

      // h-refine
      for(index_t i = 0; i< numRefine; ++i)
        mp.patch(p).uniformRefine();
    }

    gsMultiBasis<> dbasis(mp,true);
    gsInfo<<"Basis (patch 0): "<< mp.patch(0).basis() << "\n";
    if (gsTensorNurbsBasis<2,real_t> * basis = dynamic_cast<gsTensorNurbsBasis<2,real_t> * >(&mp.basis(0)))
      for (index_t dim = 0; dim!=2; dim++)
        gsDebug<<"dir "<<dim<<": "<<basis->knots(dim).asMatrix()<<"\n";
    
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
    gsInfo<<"Reading buckling solver options (ID=94) ...";
    if ( fd.hasId(94) ) fd.getId(94,BucklingOptions);
    gsInfo<<"Finished\n";
    // ![Read data]

    std::string output =  "modes";

    // Fix path
    dirname = gsFileManager::getCanonicRepresentation(dirname,true);
    char sep = gsFileManager::getNativePathSeparator();
    gsFileManager::mkdir(dirname);

    // plot geometry
    if (plot)
      gsWriteParaview(mp,"mp",1000,true);

    mp.computeTopology();

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
    }
    else if (method==3)
    {
        // The approx. C1 space
        gsApproxC1Spline<2,real_t> approxC1(mp,dbasis);
        // approxC1.options().setSwitch("info",info);
        // approxC1.options().setSwitch("plot",plot);
        approxC1.options().setSwitch("interpolation",true);
        approxC1.options().setInt("gluingDataDegree",-1);
        approxC1.options().setInt("gluingDataSmoothness",-1);
        approxC1.update(bb2);
        geom = mp;
    }
    else if (method==4)
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

    // Make assembler
    gsThinShellAssembler<3, real_t, true > assembler(geom,dbasis,BCs,forceFun,materialMatrixContainer);
    assembler.setOptions(assemblerOptions);
    assembler.options().setInt("Continuity",-1);
    assembler.setSpaceBasis(bb2);

    assembler.setPointLoads(pLoads);
    if (pressure)
      assembler.setPressure(pressFun);

    // Assemble linear system to obtain the force vector
    assembler.assemble();
    gsSparseMatrix<> K_L = assembler.matrix();
    gsVector<> F = assembler.rhs();

    // Function for the Jacobian
    gsStructuralAnalysisOps<real_t>::Jacobian_t Jacobian;
    Jacobian = [&assembler,&bb2,&geom](gsVector<real_t> const &x, gsSparseMatrix<real_t> & m)
    {
      ThinShellAssemblerStatus status;
      gsMatrix<real_t> solFull = assembler.fullSolutionVector(x);
      size_t d = geom.targetDim();
      GISMO_ASSERT(solFull.rows() % d==0,"Rows of the solution vector does not match the number of control points");
      solFull.resize(solFull.rows()/d,d);
      gsMappedSpline<2,real_t> mspline(bb2,solFull);
      gsFunctionSum<real_t> def(&geom,&mspline);
      assembler.assembleMatrix(def);
      status = assembler.assembleMatrix(def);
      m = assembler.matrix();
      return status == ThinShellAssemblerStatus::Success;
    };

    gsSparseSolver<>::CGDiagonal linearSolver;
    linearSolver.compute(K_L);
    gsVector<> solVector = linearSolver.solve(F);

    gsInfo<<"Assembling nonlinear stiffness matrix..."<<std::flush;
    gsMatrix<> solFull = assembler.fullSolutionVector(solVector);
    GISMO_ASSERT(solFull.rows() % 3==0,"Rows of the solution vector does not match the number of control points");
    solFull.resize(solFull.rows()/3,3);
    gsMappedSpline<2,real_t> mspline(bb2,solFull);
    gsFunctionSum<real_t> def(&geom,&mspline);

    gsField<> solField(geom, mspline,true);
    gsWriteParaview(solField,"LinearSolution");

    gsSparseMatrix<> K_NL;
    GISMO_ENSURE(Jacobian(solVector,K_NL),"Jacobian assembly failed.");
    K_NL -= K_L;

    gsBucklingSolver<real_t> solver(K_L,K_NL);
    solver.setOptions(BucklingOptions);

    gsStatus status = solver.computeSparse(nmodes);//,2,Spectra::SortRule::LargestMagn,Spectra::SortRule::SmallestMagn);

    GISMO_ENSURE(status == gsStatus::Success,"Buckling solver failed");

    gsMatrix<> values = solver.values();
    gsMatrix<> vectors = solver.vectors();

    gsInfo<< "First 10 eigenvalues:\n";
    for (index_t k = 0; k<10; k++)
        gsInfo<<"\t"<<std::setprecision(20)<<values.at(k)<<"\n";
    gsInfo<<"\n";

    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";
        gsMatrix<> modeShape;
        gsParaviewCollection collection(dirname + sep + "modes");

        int N = 1;
        // if (!first)
          N = nmodes;
        //    N = vectors.cols();
        for (index_t m=0; m<N; m++)
        {
            gsVector<> solVector = vectors.col(m).normalized();

            /// Make a gsMappedSpline to represent the solution
            // 1. Get all the coefficients (including the ones from the eliminated BCs.)
            gsMatrix<real_t> solFull = assembler.fullSolutionVector(solVector);
            gsMatrix<real_t> solZero = solFull;
            solZero.setZero();

            // 2. Reshape all the coefficients to a Nx3 matrix
            GISMO_ASSERT(solFull.rows() % 3==0,"Rows of the solution vector does not match the number of control points");
            solZero.resize(solZero.rows()/3,3);
            solFull.resize(solFull.rows()/3,3);

            // 3. Make the mapped spline
            gsMappedSpline<2,real_t> mspline(bb2,solFull);

            gsField<> solField(geom, mspline,true);

            std::string fileName = dirname + sep + "modes" + util::to_string(m) + "_";
            gsWriteParaview<>(solField, fileName, 1000,false);
            for (size_t p = 0; p!=geom.nPatches(); p++)
            {
                fileName = output + util::to_string(m);
                collection.addTimestep(fileName,p,m,".vts");
            }
        }
        collection.save();
    }
    if (write)
    {
        std::ofstream file;
        file.open(dirname + sep + "eigenvalues.csv",std::ofstream::out);
        for (index_t k=0; k!=values.size(); k++)
            file<<std::setprecision(12)<<values.at(k)<<"\n";

        file.close();
    }
    return EXIT_SUCCESS;
}
