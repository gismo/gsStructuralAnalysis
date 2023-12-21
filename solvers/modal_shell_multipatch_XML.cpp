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

#include <gsStructuralAnalysis/src/gsEigenSolvers/gsModalSolver.h>

#include <gsStructuralAnalysis/src/gsStructuralAnalysisTools/gsStructuralAnalysisUtils.h>

#include <gsUnstructuredSplines/src/gsSmoothInterfaces.h>
#include <gsUnstructuredSplines/src/gsAlmostC1.h>
#include <gsUnstructuredSplines/src/gsDPatch.h>

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

    index_t method = 0;

    // Arc length method options
    std::string bvp;

    std::string dirname = "ModalResults";

    index_t nmodes = 1;

    std::string wn("data.csv");

    gsCmdLine cmd("Example for an inflating balloon.");

    cmd.addInt("r","hRefine", "Number of dyadic h-refinement (bisection) steps to perform before solving", numRefine);
    cmd.addInt("e","degreeElevation", "Number of degree elevation steps to perform on the Geometry's basis before solving", numElevate);
    cmd.addInt("m","method", "Smoothing method to use: 0: smoothInterfaces, 1: Almost C1, 2: D-Patch", method);

    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addSwitch("write", "Write data to file", write);
    cmd.addSwitch("stress", "Plot stress in ParaView format", stress);

    cmd.addString("i","inputFile", "Input file", bvp);
    cmd.addString("o","outputDir", "Output directory", dirname);

    cmd.addInt( "N", "nmodes", "Number of modes",  nmodes );

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

    gsOptionList ModalOptions, assemblerOptions;
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
    gsInfo<<"Reading modal solver options (ID=95) ...";
    if ( fd.hasId(95) ) fd.getId(95,ModalOptions);
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
    }
    else
      GISMO_ERROR("Method "<<method<<" unknown");

    bb2.init(dbasis,global2local);
    // Finalize BCs
    BCs.setGeoMap(geom);

    // Make assembler
    gsThinShellAssembler<3, real_t, true > assembler(geom,dbasis,BCs,forceFun,materialMatrixContainer);
    assembler.setOptions(assemblerOptions);
    assembler.options().setInt("Continuity",-1);
    assembler.setSpaceBasis(bb2);

    // Assemble linear system to obtain the force vector
    assembler.assemble();
    gsSparseMatrix<> K = assembler.matrix();
    assembler.assembleMass();
    gsSparseMatrix<> M = assembler.matrix();

    gsModalSolver<real_t> solver(K,M);
    solver.setOptions(ModalOptions);

    // buckling.options().setInt("solver",2);
    // buckling.options().setInt("selectionRule",0);
    // buckling.options().setInt("sortRule",4);
    // buckling.options().setSwitch("verbose",true);
    // buckling.options().setInt("ncvFac",2);
    // buckling.options().setInt("shift",0.1);

    gsStatus status = solver.computeSparse(nmodes);//,2,Spectra::SortRule::LargestMagn,Spectra::SortRule::SmallestMagn);

    GISMO_ENSURE(status == gsStatus::Success,"Modal solver failed");

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
            // gsWriteParaview<>(solField, fileName, 1000,mesh);
            for (index_t p = 0; p!=geom.nPatches(); p++)
            {
                fileName = output + util::to_string(m);
                collection.addTimestep(fileName,p,m,".vts");
                // if (mesh)
                //     collection.addTimestep(fileName,p,m,"_mesh.vtp");
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
