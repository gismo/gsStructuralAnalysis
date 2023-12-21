/** @file benchmark_Roof_DWR.cpp

    @brief Collapse of a cylindrical roof using DWR adaptivity

    Inspired from

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

#include <gsKLShell/src/gsThinShellAssembler.h>
#include <gsKLShell/src/gsThinShellAssemblerDWR.h>
#include <gsKLShell/src/gsThinShellDWRHelper.h>
#include <gsKLShell/src/getMaterialMatrix.h>

#include <gsAssembler/gsAdaptiveRefUtils.h>
#include <gsAssembler/gsAdaptiveMeshing.h>

#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMBase.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMRiks.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMLoadControl.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMCrisfield.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMConsistentCrisfield.h>

using namespace gismo;

template <class T>
void initStepOutput( const std::string name, const gsMatrix<T> & points);

template <class T>
void writeStepOutput(const T deformationNorm, const T duNorm, const T duOldNorm, const T L, const T indicator, const gsMultiPatch<T> & deformation, const T error, const index_t nDoFs, const std::string name, const gsMatrix<T> & points, const index_t extreme=-1, const index_t kmax=100);

template <class T>
void PlotResults(   index_t k,
                    gsThinShellAssemblerDWRBase<T> * assembler,
                    const gsMultiPatch<T> & mp, const gsMultiPatch<T> & mp_def,
                    bool plot, bool stress, bool write, bool mesh, bool deformed,
                    const std::string dirname, const std::string output,
                    gsParaviewCollection & collection,
                    gsParaviewCollection & Smembrane,
                    gsParaviewCollection & Sflexural,
                    gsParaviewCollection & Smembrane_p);

int main (int argc, char** argv)
{
    // Input options
    int numElevate    = 1;
    int numHref       = 2;
    bool plot         = false;
    bool plotError  = false;
    bool mesh         = false;
    bool stress       = false;
    bool membrane     = false;
    bool quasiNewton  = false;
    int quasiNewtonInt= -1;
    bool adaptive     = false;
    bool adaptiveMesh = false;
    bool admissible = true;
    int maxSteps          = 500;
    int method        = 2; // (0: Load control; 1: Riks' method; 2: Crisfield's method; 3: consistent crisfield method)
    bool deformed     = false;
    bool symmetry     = true;

    bool interior = true;

    bool composite = false;

    real_t relax      = 1.0;

    int testCase      = 0;

    int result        = 0;

    bool write        = false;

    index_t maxit     = 20;
    index_t maxRefIt = 5;

    // Arc length method options
    real_t dL        = 0.5; // Ard length to find bifurcation
    real_t tol        = 1e-6;
    real_t tolU       = 1e-6;
    real_t tolF       = 1e-3;

    real_t target   =1e-3;
    real_t bandwidth=1;

    index_t goal = 1;
    index_t component = 9;

    std::string wn("data.csv");

    std::string dirname = "ArcLengthResults";

    std::string assemberOptionsFile("options/solver_options.xml");
    std::string mesherOptionsFile("options/mesher_options.xml");

    gsCmdLine cmd("Arc-length analysis of a collapsing roof.");
    cmd.addString( "o", "assemblerOpt", "Input XML file for assembler options", assemberOptionsFile );
    cmd.addString( "O", "mesherOpt", "Input XML file for mesher options", mesherOptionsFile );

    cmd.addInt("t", "testcase", "Test case: 0: clamped-clamped, 1: pinned-pinned, 2: clamped-free", testCase);

    cmd.addInt("r","hRefine", "Number of dyadic h-refinement (bisection) steps to perform before solving", numHref);
    cmd.addInt("e","degreeElevation", "Number of degree elevation steps to perform on the Geometry's basis before solving", numElevate);
    cmd.addSwitch("composite", "Composite material", composite);

    cmd.addInt("m","Method", "Arc length method; 1: Crisfield's method; 2: RIks' method.", method);
    cmd.addReal("L","dL", "arc length", dL);
    cmd.addReal("A","relaxation", "Relaxation factor for arc length method", relax);

    cmd.addInt("q","QuasiNewtonInt","Use the Quasi Newton method every INT iterations",quasiNewtonInt);
    cmd.addInt("N", "maxsteps", "Maximum number of load steps", maxSteps);
    cmd.addInt("i", "maxRefIt", "Maximum number of refinement iterations per load step", maxRefIt);

    cmd.addString("U","output", "outputDirectory", dirname);

    cmd.addReal("T","target", "Refinement target error", target);
    cmd.addReal("B","band", "Refinement target error bandwidth", bandwidth);

    cmd.addInt( "g", "goal", "Goal function to use", goal );
    cmd.addInt( "C", "comp", "Component", component );

    cmd.addSwitch("adaptive", "Adaptive length ", adaptive);
    cmd.addSwitch("adaptMesh", "Adaptive mesh ", adaptiveMesh);
    // cmd.addSwitch("admissible", "Admissible refinement", admissible);
    cmd.addSwitch("quasi", "Use the Quasi Newton method", quasiNewton);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addSwitch("noInterior", "Error computation not on the interior", interior);
    cmd.addSwitch("plotError", "Plot error in ParaView format", plotError);
    cmd.addSwitch("mesh", "Plot mesh?", mesh);
    cmd.addSwitch("stress", "Plot stress in ParaView format", stress);
    cmd.addSwitch("write", "Write output to file", write);
    cmd.addSwitch("deformed", "plot on deformed shape", deformed);
    cmd.addSwitch("nosymmetry", "Turn off symmetry", symmetry);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsFileData<> fd(assemberOptionsFile);
    gsOptionList opts;
    fd.getFirst<gsOptionList>(opts);

    gsMultiPatch<> mp, mp0;

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

    mp0 = mp;

    for (index_t i = 0; i< numElevate; ++i)
        mp.patch(0).degreeElevate();    // Elevate the degree


    // Cast all patches of the mp object to THB splines
    if (adaptiveMesh)
    {
        gsMultiPatch<> mp_thb;
        gsTHBSpline<2,real_t> thb;
        for (index_t k=0; k!=mp.nPatches(); ++k)
        {
            if(gsTensorBSpline<2,real_t> *geo = dynamic_cast< gsTensorBSpline<2,real_t> * > (&mp.patch(k)))
            {
                thb = gsTHBSpline<2,real_t>(geo->basis().source(),geo->coefs());
                gsMatrix<> bbox = geo->support();
                for (index_t i = 0; i< numHref; ++i)
                    thb.refineElements(thb.basis().asElements(bbox));

                mp_thb.addPatch(thb);
            }
        }
        mp = mp_thb;
        // h-refine
    }
    else
    {
        for (index_t i = 0; i< numHref; ++i)
            mp.patch(0).uniformRefine();
    }

    gsMultiBasis<> basisL(mp);
    gsMultiBasis<> basisH(mp);
    basisH.degreeElevate(1);
    gsInfo<<"Basis (patch 0): "<< mp.patch(0).basis() << "\n";

    // Boundary conditions
    gsBoundaryConditions<> BCs;
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

    gsMatrix<> writePoints, midPoint(2,1);
    if (!symmetry)
    {
      writePoints.resize(2,3);
      writePoints.col(0)<< 0.0,0.5;
      writePoints.col(1)<< 0.5,0.5;
      writePoints.col(2)<< 1.0,0.5;
      midPoint.col(0)<<0.5,0.5;
    }
    else
    {
      writePoints.resize(2,2);
      writePoints.col(0)<< 0.0,1.0;
      writePoints.col(1)<< 1.0,1.0;
      midPoint.col(0)<<1.0,1.0;
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

    BCs.setGeoMap(mp0);

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

    dirname = dirname + "/" +  "Roof_t="+ std::to_string(thickness) + "-r=" + std::to_string(numHref) + "-e" + std::to_string(numElevate) + "-g" + std::to_string(goal) + "-C" + std::to_string(component);
    if (symmetry)
        dirname = dirname + "_symm";
    if (adaptiveMesh)
        dirname = dirname + "_adaptive";

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
    gsMultiPatch<> mp_def = mp0;

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
    phi1.setValue(pi/2,3);
    phi2.setValue(0,3);
    phi3.setValue(pi/2,3);

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

    std::vector<gsFunctionSet<real_t>*> parameters;
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

    gsThinShellAssemblerDWRBase<real_t>* assembler;
    assembler = new gsThinShellAssemblerDWR<3, real_t, true >(mp,basisL,basisH,BCs,force,materialMatrix);
    if (goal==1)
        assembler->setGoal(GoalFunction::Displacement,component);
    else if (goal==2)
        assembler->setGoal(GoalFunction::Stretch,component);
    else if (goal==3)
        assembler->setGoal(GoalFunction::MembraneStrain,component);
    else if (goal==4)
        assembler->setGoal(GoalFunction::PStrain,component);
    else if (goal==5)
        assembler->setGoal(GoalFunction::MembraneStress,component);
    else if (goal==6)
        assembler->setGoal(GoalFunction::PStress,component);
    else if (goal==7)
        assembler->setGoal(GoalFunction::MembraneForce,component);
    else if (goal==8)
        assembler->setGoal(GoalFunction::FlexuralStrain,component);
    else if (goal==9)
        assembler->setGoal(GoalFunction::FlexuralStress,component);
    else if (goal==10)
        assembler->setGoal(GoalFunction::FlexuralMoment,component);
    else
        GISMO_ERROR("Goal function unknown");

    // Construct assembler object
    assembler->setOptions(opts);
    assembler->setPointLoads(pLoads);

    gsStopwatch stopwatch;
    real_t time = 0.0;

    // Function for the Jacobian
    gsStructuralAnalysisOps<real_t>::Jacobian_t Jacobian = [&time,&stopwatch,&assembler,&mp_def](gsVector<real_t> const &x, gsSparseMatrix<real_t> & m)
    {
        stopwatch.restart();
        ThinShellAssemblerStatus status;
        assembler->constructSolutionL(x,mp_def);
        status = assembler->assembleMatrixL(mp_def);
        m = assembler->matrixL();
        time += stopwatch.stop();
        return status == ThinShellAssemblerStatus::Success;
    };

    // Function for the Residual
    gsStructuralAnalysisOps<real_t>::ALResidual_t ALResidual = [&pLoads,&time,&stopwatch,&assembler,&mp_def](gsVector<real_t> const &x, real_t lam, gsVector<real_t> & result)
    {
        gsPointLoads<real_t> pLoads_tmp = pLoads;
        pLoads_tmp[0].value *= lam;

        stopwatch.restart();
        ThinShellAssemblerStatus status;
        assembler->setPointLoads(pLoads_tmp);
        assembler->constructSolutionL(x,mp_def);
        status = assembler->assemblePrimalL(mp_def);
        result = -assembler->primalL(); // assembler rhs - force = Finternal
        time += stopwatch.stop();
        return status == ThinShellAssemblerStatus::Success;
    };

    // Assemble linear system to obtain the force vector
    assembler->setPointLoads(pLoads);
    assembler->assembleL();
    gsVector<> Force = assembler->primalL();

    gsParaviewCollection collection(dirname + "/" + output);
    gsParaviewCollection Smembrane(dirname + "/" + "membrane");
    gsParaviewCollection Sflexural(dirname + "/" + "flexural");
    gsParaviewCollection Smembrane_p(dirname + "/" + "membrane_p");

// Make objects for previous solutions
    real_t Lold = 0, deltaLold = 0;
    real_t L = 0, deltaL = 0;
    gsMatrix<> U(Force.size(),1), deltaU(Force.size(),1);
    U.setZero();
    deltaU.setZero();
    gsMatrix<> Uold(Force.size(),1), deltaUold(Force.size(),1);
    Uold.setZero();
    deltaUold.setZero();

    gsMatrix<> solVector;
    real_t indicator_prev = 0.0;
    real_t indicator = 0.0;
    bool bisected = false;
    bool unstable_prev = false;
    real_t dL0 = dL;

    gsFileData<> fd_mesher(mesherOptionsFile);
    gsOptionList mesherOpts;
    fd_mesher.getFirst<gsOptionList>(mesherOpts);
    gsAdaptiveMeshing<real_t> mesher;
    if (adaptiveMesh)
    {
        mesher = gsAdaptiveMeshing<real_t>(mp);
        mesher.options() = mesherOpts;
        mesher.getOptions();
    }

    gsHBoxContainer<2,real_t> markRef, markCrs;

    gsMultiPatch<> U_patch, deltaU_patch;
    gsMultiPatch<> Uold_patch, deltaUold_patch;

    assembler->constructMultiPatchL(Uold,Uold_patch);
    assembler->constructMultiPatchL(deltaUold,deltaUold_patch);

    std::vector<std::vector<std::pair<index_t,real_t>>> write_errors; // per load step, iteration, numDoFs, error
    std::vector<std::pair<index_t,real_t>> loadstep_errors;

    gsALMCrisfield<real_t> arcLength(Jacobian, ALResidual, Force);
    gsOptionList ALMoptions = arcLength.options();

#ifdef GISMO_WITH_PARDISO
    ALMoptions.setString("Solver","PardisoLU"); // LDLT solver
#else
    ALMoptions.setString("Solver","SimplicialLDLT"); // LDLT solver
#endif
    ALMoptions.setInt("BifurcationMethod",0); // 0: determinant, 1: eigenvalue
    ALMoptions.setReal("Length",dL);
    ALMoptions.setInt("AngleMethod",0); // 0: step, 1: iteration
    ALMoptions.setSwitch("AdaptiveLength",adaptive);
    ALMoptions.setInt("AdaptiveIterations",5);
    ALMoptions.setReal("Scaling",0.0);
    ALMoptions.setReal("Tol",tol);
    ALMoptions.setReal("TolU",tolU);
    ALMoptions.setReal("TolF",tolF);
    ALMoptions.setInt("MaxIter",maxit);
    ALMoptions.setSwitch("Verbose",true);
    ALMoptions.setReal("Relaxation",relax);
    if (quasiNewtonInt>0)
    {
      quasiNewton = true;
      ALMoptions.setInt("QuasiIterations",quasiNewtonInt);
    }
    ALMoptions.setSwitch("Quasi",quasiNewton);

    gsInfo<<ALMoptions;

    arcLength.options() = ALMoptions;
    arcLength.applyOptions();
    arcLength.initialize();

    gsThinShellDWRHelper<real_t> helper(assembler);
    typename gsBoxTopology::bContainer goalSides;
    //goalSides.push_back(patchSide(0,boundary::west));
    gsMatrix<> points;
    real_t error = std::numeric_limits<real_t>::max();
    index_t numDofs = assembler->numDofsL();

    real_t refTol = target / bandwidth; // refine if error is above
    real_t crsTol = target * bandwidth; // coarsen if error is below
    real_t Umid = 0;
    real_t Umidmax =  30;
    real_t Umidmin = -30;
    index_t k = 0;

    gsParaviewCollection errors(dirname + "/" + "errors");
    std::vector<real_t> elErrors;
    GISMO_ENSURE(refTol >= crsTol,"Refinement tolerance should be bigger than the coarsen tolerance");
    while (Umid < Umidmax && Umid > Umidmin && k < maxSteps)
    {
        loadstep_errors.clear();
        gsInfo<<"Load step "<< k<<"; \tUmid,L = "<<Umid<<","<<Lold<<"\tSystem size = "<<Uold.size()<<" x "<<Uold.size()<<"\n";

        index_t maxIt = 10;
        index_t it = 0;
        bool refined = true;
        bool coarsened = true;
        error = std::numeric_limits<real_t>::max();
        bool bandtest = (bandwidth==1) ? error > refTol : ((error < crsTol )|| (error >= refTol));
        gsMultiPatch<> mp_prev;
        while (bandtest && it < maxRefIt && (refined || coarsened))
        {
            gsInfo<<"Iteration "<<it<<"/"<<maxRefIt<<", refTol < prev error < crsTol : "<<refTol<<" < "<<error<<" < "<<crsTol<<"\n";
            gsInfo<<"New basis (L): \n"<<mp.basis(0)<<"\n";
            mp_prev = mp;

            assembler->setPointLoads(pLoads);
            assembler->assembleL();
            Force = assembler->primalL();
            Uold = assembler->constructSolutionVectorL(Uold_patch);
            deltaUold = assembler->constructSolutionVectorL(deltaUold_patch);

            gsALMCrisfield<real_t> arcLength(Jacobian, ALResidual, Force);
            arcLength.options() = ALMoptions;
            arcLength.applyOptions();
            arcLength.initialize();
            arcLength.setIndicator(indicator); // RESET INDICATOR
            arcLength.setSolution(Uold,Lold);
            arcLength.setSolutionStep(deltaUold,deltaLold);
            arcLength.setLength(dL);

            gsInfo<<"Starting from U.norm()="<<Uold.norm()<<", L="<<Lold<<"\n";
            arcLength.step();

            if (!(arcLength.converged()))
            {
              gsInfo<<"Error: Loop terminated, arc length method did not converge.\n";
              dL = dL / 2.;
              arcLength.setLength(dL);
              arcLength.setSolution(Uold,Lold);
              bisected = true;
              it -= 1;
              continue;
            }
            indicator = arcLength.indicator();
            gsInfo<<"indicator: (old = )"<<indicator_prev<<"; (new = )"<<indicator<<"\n";

            L = arcLength.solutionL();
            deltaL = arcLength.solutionDL();
            U = arcLength.solutionU();
            deltaU = arcLength.solutionDU();

            // Deformed geometry
            assembler->constructSolutionL(U,mp_def);
            // Deformation (primal)
            assembler->constructMultiPatchL(U,U_patch);
            // delta Deformation
            assembler->constructMultiPatchL(deltaU,deltaU_patch);

            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // ERROR ESTIMATION PART
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            helper.computeError(mp_def,U_patch,goalSides,points,interior,false);

            error = std::abs(helper.error());
            numDofs = assembler->numDofsL();

            gsInfo<<"Error = "<<error<<", numDofs = "<<numDofs<<"\n";
            loadstep_errors.push_back(std::make_pair(assembler->numDofsL(),error));

            elErrors = helper.sqErrors(true);

            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // ADAPTIVE MESHING PART
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            if (adaptiveMesh)
            {
                if (unstable_prev)
                {
                    unstable_prev = false;
                    break;
                }
                else
                {
                    if (error > refTol)
                    {
                        gsInfo<<"Load Step "<<k<<": Error is too big! Error = "<<error<<", refTol = "<<refTol<<"\n";
                        mesher.markRef_into(elErrors,markRef);
                        gsDebugVar(markRef);
                        gsInfo<<"Marked "<<markRef.totalSize()<<" elements for refinement\n";
                        refined = mesher.refine(markRef);
                    }
                    else if (error < refTol && error > crsTol)
                    {
                        gsInfo<<"Load Step "<<k<<": Error is within bounds. Error = "<<error<<", refTol = "<<refTol<<", crsTol = "<<crsTol<<"\n";
                        mesher.markRef_into(elErrors,markRef);
                        mesher.markCrs_into(elErrors,markRef,markCrs);
			            gsInfo<<"Marked "<<markRef.totalSize()<<" elements for refinement\n";
                        gsInfo<<"Marked "<<markCrs.totalSize()<<" elements for coarsening\n";
                        refined = mesher.refine(markRef);
                        coarsened = mesher.unrefine(markCrs);
                        gsInfo<<"No elements marked\n";
                    }
                    else if (error < crsTol)
                    {
                        //gsInfo<<"Error is too small!\n";
                        gsInfo<<"Load Step "<<k<<": Error is too small! Error = "<<error<<", crsTol = "<<crsTol<<"\n";
                        mesher.markCrs_into(elErrors,markCrs);
                        gsInfo<<"Marked "<<markCrs.totalSize()<<" elements for coarsening\n";
                        coarsened = mesher.unrefine(markCrs);
                    }
                    bandtest = (bandwidth==1) ? error > refTol : ((error < crsTol )|| (error >= refTol));

                    basisL = gsMultiBasis<>(mp);
                    basisH = basisL;
                    basisH.degreeElevate(1);

                    // Project the solution from old mesh to new mesh
                    gsMatrix<> coefs;

                    // Which of those are needed?

                    gsQuasiInterpolate<real_t>::localIntpl(basisL.basis(0), mp.patch(0), coefs);
                    mp.patch(0) = *basisL.basis(0).makeGeometry(give(coefs));

                    gsQuasiInterpolate<real_t>::localIntpl(basisL.basis(0), mp_def.patch(0), coefs);
                    mp_def.patch(0) = *basisL.basis(0).makeGeometry(give(coefs));

                    gsQuasiInterpolate<real_t>::localIntpl(basisL.basis(0), U_patch.patch(0), coefs);
                    U_patch.patch(0) = *basisL.basis(0).makeGeometry(give(coefs));

                    gsQuasiInterpolate<real_t>::localIntpl(basisL.basis(0), Uold_patch.patch(0), coefs);
                    Uold_patch.patch(0) = *basisL.basis(0).makeGeometry(give(coefs));

                    gsQuasiInterpolate<real_t>::localIntpl(basisL.basis(0), deltaUold_patch.patch(0), coefs);
                    deltaUold_patch.patch(0) = *basisL.basis(0).makeGeometry(give(coefs));

                    assembler->setBasisL(basisL);
                    assembler->setBasisH(basisH);
                    assembler->setUndeformed(mp);

                    mesher.rebuild();

                    // assembler->constructSolutionL(U,mp_def);
                    unstable_prev = false;

                }
                it++;
            }
            else
                break;
        }

        if (plotError)
        {
            for (size_t p=0; p!=mp_prev.nPatches(); p++)
            {
                gsElementErrorPlotter<real_t> err_eh(mp_prev.basis(p),elErrors);
                const gsField<> elemError_eh( mp_prev.patch(p), err_eh, true );
                std::string fileName = dirname + "/" + "error" + util::to_string(k);
                writeSinglePatchField<>(mp_prev.patch(p), err_eh, true, fileName + "_" + util::to_string(p), 1000);
                if (mesh)
                    writeSingleCompMesh<>(mp_prev.basis(p), mp_prev.patch(p),fileName + "_mesh" + "_" + util::to_string(p));
                fileName = "error" + util::to_string(k);
                errors.addTimestep(fileName,p,k,".vts");
                if (mesh)
                    errors.addTimestep(fileName + "_mesh",p,k,".vtp");
            }
        }

        deltaU_patch = U_patch;
        for (index_t p=0; p!=deltaU_patch.nPatches(); p++)
            deltaU_patch.patch(p).coefs() -= Uold_patch.patch(p).coefs();

        Umid = U_patch.patch(0).eval(midPoint)(2,0);

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        real_t deformationNorm  = assembler->deformationNorm(U_patch,mp);
        real_t duNorm           = assembler->deformationNorm(deltaU_patch,mp);
        real_t duOldNorm        = assembler->deformationNorm(deltaUold_patch,mp);

        PlotResults(k,assembler,mp,mp_def,plot,stress,write,mesh,deformed,dirname,output,
                    collection,Smembrane,Sflexural,Smembrane_p);

        if (write)
            writeStepOutput(deformationNorm,duNorm,duOldNorm,L,indicator,U_patch, error, numDofs, dirname + "/" + wn, writePoints);

        write_errors.push_back(loadstep_errors);

        // Update Uold
        Uold_patch = U_patch;
        deltaUold_patch = deltaU_patch;
        Lold = L;
        deltaLold = deltaL;

        indicator_prev = indicator;
        k++;
    }

    if (plotError)
    {
        errors.save();
    }

    if (plot)
    {
      collection.save();
    }
    if (stress)
    {
      Smembrane.save();
      Sflexural.save();
      Smembrane_p.save();
    }

    std::ofstream file;
    file.open(dirname + "/" + "errors.csv",std::ofstream::out);
    index_t loadstep=0;
    file<<"load_step,iteration,numDofs,error\n";
    for (std::vector<std::vector<std::pair<index_t,real_t>>>::const_iterator it = write_errors.begin(); it!=write_errors.end(); it++, loadstep++)
    {
        index_t iteration=0;
        for (std::vector<std::pair<index_t,real_t>>::const_iterator iit = it->begin(); iit!=it->end(); iit++, iteration++)
            file<<loadstep<<","<<iteration<<","<<iit->first<<","<<iit->second<<"\n";

    }
    file.close();

    delete materialMatrix;
    delete assembler;

    return result;
}

template <class T>
void initStepOutput(const std::string name, const gsMatrix<T> & points)
{
  std::ofstream file;
  file.open(name,std::ofstream::out);
  file  << std::setprecision(20)
        << "Deformation norm" << ","
        << "DU norm" << ","
        << "DU old norm" << ",";
        for (index_t k=0; k!=points.cols(); k++)
        {
          file<< "point "<<k<<" - x" << ","
              << "point "<<k<<" - y" << ","
              << "point "<<k<<" - z" << ",";
        }

  file  << "Lambda" << ","
        << "Indicator" << ","
        << "NumDofs" << ","
        << "Error"
        << "\n";
  file.close();

  gsInfo<<"Step results will be written in file: "<<name<<"\n";
}

template <class T>
void writeStepOutput(const T deformationNorm, const T duNorm, const T duOldNorm, const T L, const T indicator, const gsMultiPatch<T> & deformation, const T error, const index_t nDoFs, const std::string name, const gsMatrix<T> & points, const index_t extreme, const index_t kmax) // extreme: the column of point indices to compute the extreme over (default -1)
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
          << deformationNorm << ","
          << duNorm << ","
          << duOldNorm << ",";
          for (index_t p=0; p!=points.cols(); p++)
          {
            file<< out(0,p) << ","
                << out(1,p) << ","
                << out(2,p) << ",";
          }

    file  << L << ","
          << indicator << ","
          << nDoFs << ","
          << error << ","
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
          << deformationNorm << ","
          << duNorm << ","
          << duOldNorm << ",";
          for (index_t p=0; p!=points.cols(); p++)
          {
            file<< out(0,p) << ","
                << out(1,p) << ","
                << std::max(abs(out2.col(p).maxCoeff()),abs(out2.col(p).minCoeff())) << ",";
          }

    file  << L << ","
          << indicator << ","
          << nDoFs << ","
          << error << ","
          << "\n";
  }
  else
    GISMO_ERROR("Extremes setting unknown");

  file.close();
}

template <class T>
void PlotResults(   index_t k,
                    gsThinShellAssemblerDWRBase<T> * assembler,
                    const gsMultiPatch<T> & mp, const gsMultiPatch<T> & mp_def,
                    bool plot, bool stress, bool write, bool mesh, bool deformed,
                    const std::string dirname, const std::string output,
                    gsParaviewCollection & collection,
                    gsParaviewCollection & Smembrane,
                    gsParaviewCollection & Sflexural,
                    gsParaviewCollection & Smembrane_p)
{
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    gsMultiPatch<T> deformation = mp_def;

    deformation.patch(0).coefs() -= mp.patch(0).coefs();// assuming 1 patch here

    gsInfo<<"Total ellapsed assembly time: "<<time<<" s\n";

    if (plot)
    {
        gsField<T> solField;
        if (deformed)
          solField= gsField<>(mp_def,deformation);
        else
          solField= gsField<>(mp,deformation);

        std::string fileName = dirname + "/" + output + util::to_string(k);
        gsWriteParaview<T>(solField, fileName, 1000,mesh);
        fileName = output + util::to_string(k) + "0";
        collection.addTimestep(fileName,k,".vts");
        if (mesh) collection.addTimestep(fileName,k,"_mesh.vtp");
    }
    if (stress)
    {
        gsField<T> membraneStress, flexuralStress, membraneStress_p;

        gsPiecewiseFunction<T> membraneStresses;
        assembler->constructStress(mp_def,membraneStresses,stress_type::membrane);
        if (deformed)
          membraneStress = gsField<>(mp_def,membraneStresses,true);
        else
          membraneStress = gsField<>(mp,membraneStresses,true);

        gsPiecewiseFunction<T> flexuralStresses;
        assembler->constructStress(mp_def,flexuralStresses,stress_type::flexural);
        if (deformed)
          flexuralStress = gsField<>(mp_def,flexuralStresses, true);
        else
          flexuralStress = gsField<>(mp,flexuralStresses, true);

        gsPiecewiseFunction<T> membraneStresses_p;
        assembler->constructStress(mp_def,membraneStresses_p,stress_type::principal_stress);
        if (deformed)
          membraneStress_p = gsField<>(mp_def,membraneStresses_p, true);
        else
          membraneStress_p = gsField<>(mp,membraneStresses_p, true);

        std::string fileName;
        fileName = dirname + "/" + "membrane" + util::to_string(k);
        gsWriteParaview( membraneStress, fileName, 1000);
        fileName = "membrane" + util::to_string(k) + "0";
        Smembrane.addTimestep(fileName,k,".vts");

        fileName = dirname + "/" + "flexural" + util::to_string(k);
        gsWriteParaview( flexuralStress, fileName, 1000);
        fileName = "flexural" + util::to_string(k) + "0";
        Sflexural.addTimestep(fileName,k,".vts");

        fileName = dirname + "/" + "membrane_p" + util::to_string(k);
        gsWriteParaview( membraneStress_p, fileName, 1000);
        fileName = "membrane_p" + util::to_string(k) + "0";
        Smembrane_p.addTimestep(fileName,k,".vts");
    }
}
