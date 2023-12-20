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
#include <gsKLShell/src/gsMaterialMatrixTFT.h>

#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMBase.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMLoadControl.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMRiks.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMCrisfield.h>

#include <gsStructuralAnalysis/src/gsStructuralAnalysisTools/gsStructuralAnalysisUtils.h>

using namespace gismo;

template <class T>
void writeStepOutput(const gsALMBase<T> * arcLength, const gsMultiPatch<T> & deformation, const std::string name, const gsMatrix<T> & points, const gsMatrix<index_t> & patches);

int main (int argc, char** argv)
{
    // Input options
    int numElevate    = 0;
    int numRefine     = 1;
    bool plot         = false;
    bool write        = false;
    bool stress       = false;

    bool SingularPoint = false;
    int step = 10;

    real_t perturb = 0;

    // Arc length method options

    std::string bvp;

    std::string dirname = ".";

    std::string wn("data.csv");

    gsCmdLine cmd("Example for an inflating balloon.");

    cmd.addInt("r","hRefine", "Number of dyadic h-refinement (bisection) steps to perform before solving", numRefine);
    cmd.addInt("e","degreeElevation", "Number of degree elevation steps to perform on the Geometry's basis before solving", numElevate);

    cmd.addReal( "P" , "perturb" , "Perturb the z coordinate of the refined geometry",  perturb );

    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addSwitch("write", "Write data to file", write);
    cmd.addSwitch("stress", "Plot stress in ParaView format", stress);
    cmd.addSwitch("bifurcation", "Compute singular points and bifurcation paths", SingularPoint);

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

    gsOptionList ALMOptions, assemblerOptions;
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

    gsInfo<<"Reading assembler options (ID=92) ...";
    if ( fd.hasId(92) ) fd.getId(92,assemblerOptions);
    gsInfo<<"Finished\n";
    gsInfo<<"Reading ALM solver options (ID=93) ...";
    if ( fd.hasId(93) ) fd.getId(93,ALMOptions);
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

    mp.computeTopology();
    gsThinShellAssemblerBase<real_t>* assembler;
    assembler = new gsThinShellAssembler<3, real_t, true  >(mp,dbasis,BCs,forceFun,materialMatrixContainer);

    if (mp.nPatches()>1)
    {
      assembler->addWeakC0(mp.topology().interfaces());
      assembler->addWeakC1(mp.topology().interfaces());
      assembler->initInterfaces();
    }
    assembler->setPointLoads(pLoads);
    if (pressure)
      assembler->setPressure(pressFun);

    // Assemble linear system to obtain the force vector
    assembler->assemble();
    gsVector<> Force = assembler->rhs();

    // Function for the Jacobian
    gsStructuralAnalysisOps<real_t>::Jacobian_t Jacobian = [&assembler,&mp_def](gsVector<real_t> const &x, gsSparseMatrix<real_t> & m)
    {
      ThinShellAssemblerStatus status;
      assembler->constructSolution(x,mp_def);
      status = assembler->assembleMatrix(mp_def);
      m = assembler->matrix();
      return status == ThinShellAssemblerStatus::Success;
    };
    // Function for the Residual
    gsStructuralAnalysisOps<real_t>::ALResidual_t ALResidual = [&assembler,&mp_def,&Force](gsVector<real_t> const &x, real_t lam, gsVector<real_t> & result)
    {
        ThinShellAssemblerStatus status;
        assembler->constructSolution(x,mp_def);
        status = assembler->assembleVector(mp_def);
        result = Force - lam * Force - assembler->rhs(); // assembler rhs - force = Finternal
        return status == ThinShellAssemblerStatus::Success;
    };

    gsStructuralAnalysisOutput<real_t> writer(dirname + sep + wn,refPoints);
    std::vector<std::string> pointheaders = {"x","y","z"};
    std::vector<std::string> otherheaders = {"L"};
    writer.init(pointheaders,otherheaders);

    gsALMBase<real_t> * arcLength;
    std::string method = ALMOptions.askString("Method","Crisfield");
    if (method=="LoadControl")
      arcLength = new gsALMLoadControl<real_t>(Jacobian, ALResidual, Force);
    else if (method=="Riks")
      arcLength = new gsALMRiks<real_t>(Jacobian, ALResidual, Force);
    else if (method=="Crisfield")
      arcLength = new gsALMCrisfield<real_t>(Jacobian, ALResidual, Force);
    else
      GISMO_ERROR("Method "<<method<<" unknown");

    arcLength->setOptions(ALMOptions);
    arcLength->applyOptions();
    arcLength->initialize();

    gsParaviewCollection collection(dirname + sep + output);
    gsParaviewCollection TensionFields(dirname + sep + "tensionfield");
    gsMultiPatch<> deformation = mp;

    // Make objects for previous solutions
    real_t Lold = 0;
    gsMatrix<> Uold(Force.rows(),1);
    Uold.setZero();

    real_t indicator = 0.0;
    arcLength->setIndicator(indicator); // RESET INDICATOR
    bool bisected = false;
    real_t dLb, dLb0;
    dLb = dLb0 = arcLength->getLength();
    for (index_t k=0; k<step; k++)
    {
      gsInfo<<"Load step "<< k<<"\n";
      arcLength->step();

      if (!(arcLength->converged()))
      {
        gsInfo<<"Error: Loop terminated, arc length method did not converge.\n";
        dLb = dLb / 2.;
        arcLength->setLength(dLb);
        arcLength->setSolution(Uold,Lold);
        k -= 1;
        continue;
      }

      arcLength->computeStability(false);
      if (arcLength->stabilityChange() && SingularPoint)
      {
          gsInfo<<"Bifurcation spotted!"<<"\n";
          arcLength->computeSingularPoint(Uold,Lold,false);
          dLb0 = dLb = ALMOptions.askReal("Length2",dLb0);
          arcLength->setLength(dLb);
      }

      indicator = arcLength->indicator();

      gsMatrix<> solVector;
      solVector = arcLength->solutionU();
      Uold = solVector;
      Lold = arcLength->solutionL();

      gsMappedSpline<2,real_t> mspline;
      if (plot)
      {
        gsField<> solField;
        std::string fileName = dirname + sep + output + util::to_string(k);
        assembler->constructSolution(solVector,mp_def);

        deformation = mp_def;
        for (size_t p=0; p!=mp_def.nPatches(); p++)
          deformation.patch(p).coefs() -= mp.patch(p).coefs();// assuming 1 patch here

        solField= gsField<>(mp_def,deformation);

        gsWriteParaview<>(solField, fileName, 1000,false,"_");
        // gsWriteParaview<>(solField, fileName, 1000,mesh,"_");

        for (size_t p = 0; p!=mp.nPatches(); p++)
        {
          fileName = output + util::to_string(k) + "_" + util::to_string(p);
          collection.addPart(fileName + ".vts",k);
          // if (mesh) collection.addPart(fileName + "_mesh.vtp",k);
        }
      }
      if (stress)
      {
        gsField<> tensionField;
        gsPiecewiseFunction<> tensionFields;
        std::string fileName;
        fileName = dirname + sep + "tensionfield" + util::to_string(k);

        assembler->constructStress(mp_def,tensionFields,stress_type::tension_field);
        tensionField = gsField<>(mp_def,tensionFields, true);

        gsWriteParaview( tensionField, fileName, 1000,false,"_");

        for (size_t p = 0; p!=mp.nPatches(); p++)
        {
          fileName = "tensionfield" + util::to_string(k) + "_" + util::to_string(p);
          TensionFields.addPart(fileName + ".vts",k);
          // if (mesh) TensionFields.addPart(fileName + "_mesh.vtp",k);
        }
      }

      if (write)
      {
        gsMultiPatch<> deformation;
        assembler->constructDisplacement(solVector,deformation);

        gsMatrix<> result(mp.geoDim(),refPoints.cols());
        for (index_t k=0; k!=refPoints.cols(); k++)
          result.col(k) = deformation.patch(refPatches.at(k)).eval(refPoints.col(k));

        gsVector<> otherdata(1);
        otherdata<<Lold;
        writer.add(result,otherdata);

        // writeStepOutput(arcLength,deformation, dirname + sep + wn, refPoints,refPatches);
      }

      if (!bisected)
      {
        dLb = dLb0;
        arcLength->setLength(dLb);
      }
      bisected = false;

    }

    if (plot)
    {
      collection.save();
    }
    if (stress)
    {
      TensionFields.save();
    }

    delete assembler;

    return EXIT_SUCCESS;
}

template <class T>
void writeStepOutput(const gsALMBase<T> * arcLength, const gsMultiPatch<T> & deformation, const std::string name, const gsMatrix<T> & points, const gsMatrix<index_t> & patches)
{
  gsMatrix<T> P(2,1), Q(2,1);
  gsMatrix<T> out(deformation.geoDim(),points.cols());
  gsMatrix<T> tmp;

  for (index_t p=0; p!=points.cols(); p++)
  {
    P<<points.col(p);
    deformation.patch(patches(0,p)).eval_into(P,tmp);
    out.col(p) = tmp;
  }

  std::ofstream file;
  file.open(name,std::ofstream::out | std::ofstream::app);
  file  << std::setprecision(6)
        << arcLength->solutionU().norm() << ",";
        for (index_t p=0; p!=points.cols(); p++)
        {
          file<< out(0,p) << ","
              << out(1,p) << ",";
        }

  file  << arcLength->solutionL() << ","
        << arcLength->indicator() << ","
        << "\n";
  file.close();
}