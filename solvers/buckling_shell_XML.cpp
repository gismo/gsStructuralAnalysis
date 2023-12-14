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

#include <gsStructuralAnalysis/src/gsEigenSolvers/gsBucklingSolver.h>

#include <gsStructuralAnalysis/src/gsStructuralAnalysisTools/gsStructuralAnalysisUtils.h>

using namespace gismo;

int main (int argc, char** argv)
{
    // Input options
    int numElevate    = 0;
    int numRefine     = 1;
    bool plot         = false;
    bool write        = false;
    bool stress       = false;

    // Arc length method options
    std::string bvp;

    std::string dirname = "BucklingResults";

    index_t nmodes = 1;

    std::string wn("data.csv");

    gsCmdLine cmd("Example for an inflating balloon.");

    cmd.addInt("r","hRefine", "Number of dyadic h-refinement (bisection) steps to perform before solving", numRefine);
    cmd.addInt("e","degreeElevation", "Number of degree elevation steps to perform on the Geometry's basis before solving", numElevate);

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

    gsOptionList BucklingOptions, assemblerOptions;
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
    gsInfo<<"Reading buckling solver options (ID=91) ...";
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

    gsThinShellAssembler<3, real_t, true > assembler(mp,dbasis,BCs,forceFun,materialMatrixContainer);
    assembler.setOptions(assemblerOptions);
    if (mp.nPatches()>1)
    {
      assembler.addWeakC0(mp.topology().interfaces());
      assembler.addWeakC1(mp.topology().interfaces());
      assembler.initInterfaces();
    }
    assembler.setPointLoads(pLoads);
    if (pressure)
      assembler.setPressure(pressFun);

    // Assemble linear system to obtain the force vector
    assembler.assemble();
    gsSparseMatrix<> K_L = assembler.matrix();
    gsVector<> F = assembler.rhs();

    // Function for the Jacobian
    gsStructuralAnalysisOps<real_t>::Jacobian_t Jacobian = [&assembler,&mp_def](gsVector<real_t> const &x, gsSparseMatrix<real_t> & m)
    {
      // &diag,
      ThinShellAssemblerStatus status;
      assembler.constructSolution(x,mp_def);
      status = assembler.assembleMatrix(mp_def);
      m = assembler.matrix();
      // m += diag;
      return status == ThinShellAssemblerStatus::Success;
    };

    gsSparseSolver<>::CGDiagonal linearSolver;
    linearSolver.compute(K_L);
    gsVector<> solVector = linearSolver.solve(F);

    gsMultiPatch<> displacement;
    assembler.constructDisplacement(solVector,displacement);
    gsField<> field(mp,displacement);
    gsWriteParaview(field,"linearSolution");

    gsSparseMatrix<> K_NL;
    GISMO_ENSURE(Jacobian(solVector,K_NL),"Jacobian assembly failed.");

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
        gsMultiPatch<> deformation = mp;
        gsMatrix<> modeShape;
        gsParaviewCollection collection(dirname + sep + output);

        index_t N = nmodes;
        for (index_t m=0; m<N; m++)
        {

          // Compute solution based on eigenmode with number 'mode'
          modeShape = vectors.col(m);
          assembler.constructSolution(modeShape, deformation);

          // compute the deformation spline
          deformation.patch(0).coefs() -= mp.patch(0).coefs();// assuming 1 patch here

          // Normalize mode shape amplitude in z coordinate
          real_t maxAmpl = std::max(math::abs(deformation.patch(0).coefs().col(2).maxCoeff()),math::abs(deformation.patch(0).coefs().col(2).minCoeff()));
          if (maxAmpl!=0.0)
          {
            deformation.patch(0).coefs() = deformation.patch(0).coefs()/maxAmpl;
          }

          gsField<> solField(mp,deformation);
          std::string fileName = dirname + sep + output + util::to_string(m);
          gsWriteParaview<>(solField, fileName, 5000);
          fileName = output + util::to_string(m) + "0";
          collection.addPart(fileName + ".vts",m);
        }
        collection.save();
    }

    if (write)
    {
        int systemRet = system("mkdir -p BucklingResults");
        GISMO_ASSERT(systemRet!=-1,"Something went wrong with calling the system argument");
        std::string wnM = "BucklingResults/eigenvalues.txt";
        // writeToCSVfile(wnM,values);
    }
    return EXIT_SUCCESS;
}
