
/** @file benchmark_cylinder_DC.cpp

    @brief Computes the wrinkling behaviour of a cylinder or annulus

    Twisting of an annulus (test 0) or a cylinder (test 1), Verhelst et al., CMAME 441 (2025) 117955,
    Secs. 6.3 and 6.4. The driven boundary (v = 1) is rotated by pi/2*D about the z axis (D = load factor)
    and displaced axially by a constant amount; the other boundary (v = 0) has all displacement components fixed (rotations are not constrained). The geometry is
    discretised with a multi-patch basis with smoothed (C1) interfaces (gsUnstructuredSplines).

    Material: incompressible Neo-Hookean, optionally wrapped by the tension-field theory
    (--TFT, gsMaterialMatrixTFT) or by the Zhang-Kiendl split (--tftsplit, gsMaterialMatrixTFTSplit);
    --tftlinear selects the linear flag of the wrapped law for both wrappers. --eta acts on --TFT as well as on
    --tftsplit; --etatangent (--TFT only) applies it to the tangent only.
    --membrane drops the bending stiffness.

    Output (directory ArcLengthResults/<case>_-r<r>-e<e>[_TFT|_Split][_membrane], or --out):
    - out.csv: load factor and the three virtual-work moments of the driven boundary
    - steps.csv: one row per solve attempt, see the column description in the source of main()
    - ParaView files and line_*.csv with --plot

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#include <gismo.h>

#ifdef gsKLShell_ENABLED
#include <gsKLShell/src/gsThinShellAssembler.h>
#include <gsKLShell/src/getMaterialMatrix.h>
#include <gsKLShell/src/gsMaterialMatrixTFT.h>
#include <gsKLShell/src/gsMaterialMatrixTFTSplit.h>
#include <gsKLShell/src/gsFunctionSum.h>
#endif

#ifdef gsUnstructuredSplines_ENABLED
#include <gsUnstructuredSplines/src/gsSmoothInterfaces.h>
#endif

#ifdef gsOptim_ENABLED
#include <gsOptim/gsOptim.h>
#endif

#include <gsStructuralAnalysis/src/gsStaticSolvers/gsStaticDR.h>
#include <gsStructuralAnalysis/src/gsStaticSolvers/gsStaticOpt.h>
#include <gsStructuralAnalysis/src/gsStaticSolvers/gsStaticNewton.h>
#include <gsStructuralAnalysis/src/gsStaticSolvers/gsStaticComposite.h>
#include <gsStructuralAnalysis/src/gsStaticSolvers/gsControlDisplacement.h>
using namespace gismo;

void writeToCSVfile(const std::string& name, const gsMatrix<>& matrix)
{
  std::ofstream file(name.c_str());
  for(int  i = 0; i < matrix.rows(); i++)
  {
    for(int j = 0; j < matrix.cols(); j++)
    {
       if(j+1 == matrix.cols())
       {
           file<<std::setprecision(10)<<matrix(i,j);
       }
       else
       {
           file<<std::setprecision(10)<<matrix(i,j)<<',';
       }
    }
    file<<'\n';
  }
}

/// Creates \a path and its missing parents (gsFileManager::mkdir is not recursive)
bool makeDirs(const std::string & path)
{
    bool ok = true;
    for (size_t pos = path.find('/',1); pos != std::string::npos; pos = path.find('/',pos+1))
        ok = gsFileManager::mkdir(path.substr(0,pos)) && ok;
    return gsFileManager::mkdir(path) && ok;
}

#ifdef gsKLShell_ENABLED

/// Wrinkling decorator of the base law \a base and the wrinkling-state oracle.
///
/// \a model: 0 no decorator, 1 gsMaterialMatrixTFT, 2 gsMaterialMatrixTFTSplit; \a linear is the linear flag
/// of the decorators. The decorators clone \a base, so the criterion (\a tensionField) must already be set on
/// \a base. All options are set here, before the caller hands the decorator to a gsMaterialMatrixContainer,
/// which clones it. \a exactTol, \a exactStagTol and \a exactMaxIt are the ExactTol, ExactStagTol and
/// ExactMaxIt options of the Split decorator (ignored for the other models).
///
/// The oracle is a Split decorator around \a base (not owned) with Eta = 0 and S2crit = 1e300: its
/// TensionField output is the pure criterion state (1 taut, 0 wrinkled, -1 slack) of the base law,
/// independent of the law the solver integrates.
template <bool linear>
void makeWrinkleMaterials(gsMaterialMatrixBase<real_t>::uPtr & wrapped,
                          gsMaterialMatrixBase<real_t>::uPtr & oracle,
                          gsMaterialMatrixBase<real_t> & base,
                          const int model, const real_t eta, const bool etaTangent, const bool exact,
                          const real_t slackMultiplier, const index_t tensionField,
                          const real_t exactTol, const real_t exactStagTol, const index_t exactMaxIt)
{
    if (model==1)
    {
        wrapped = memory::make_unique(new gsMaterialMatrixTFT<3,real_t,linear>(base));
        wrapped->options().setReal("SlackMultiplier",slackMultiplier);
        wrapped->options().setReal("Eta",eta);
        wrapped->options().setSwitch("EtaTangentOnly",etaTangent);
    }
    else if (model==2)
    {
        wrapped = memory::make_unique(new gsMaterialMatrixTFTSplit<3,real_t,linear>(base));
        gsOptionList o = wrapped->options();
        o.setReal("Eta",eta);
        o.setSwitch("Exact",exact);
        o.setReal("ExactTol",exactTol);
        o.setReal("ExactStagTol",exactStagTol);
        o.setInt("ExactMaxIt",exactMaxIt);
        o.setReal("SlackMultiplier",slackMultiplier);
        o.setInt("TensionField",tensionField);
        wrapped->setOptions(o);
    }

    oracle = memory::make_unique(new gsMaterialMatrixTFTSplit<3,real_t,linear>(&base));
    gsOptionList o = oracle->options();
    o.setReal("Eta",0);
    o.setSwitch("Exact",false);
    o.setReal("S2crit",1e300);
    o.setReal("SlackMultiplier",slackMultiplier);
    o.setInt("TensionField",tensionField);
    oracle->setOptions(o);
}

#ifdef gsUnstructuredSplines_ENABLED

int main (int argc, char** argv)
{
    // Input options
    int numElevate  = 2;
    int numHref     = 5;
    bool plot       = false;
    bool mesh = false;
    int step = 0;

    bool    TFT = false;
    bool    tftsplit = false;
    index_t tftlinear = 1;
    real_t  eta = 0;
    bool    etaTangent = false;
    bool    exact = false;
    real_t  exactTol = 1e-12;
    real_t  exactStagTol = 1e-9;
    index_t exactMaxIt = 50;
    real_t  slackmult = 0;
    index_t tensionField = 2;
    index_t maxbisect = 8;
    std::string outDir;
    bool membrane= false;
    real_t dL = 1;
    // Solvers
    // * Dynamic Relaxation
    index_t maxitDR = 1e4;
    real_t  alpha = 1.0;
    real_t  damping = 0.1;
    bool DR  = false;
    // * Optimizer
    index_t maxitOPT = 1000;
    bool OP  = false;
    // * Newton Raphson
    index_t maxitNR = 50;
    bool lineSearch = false;
    std::string linSolver("CGDiagonal");
    real_t  linTol = 0;
    index_t linMaxIt = 0;
    bool NR  = false;


    index_t testCase = 0;


    std::string wn("data.csv");

    std::string assemberOptionsFile("options/solver_options.xml");

    gsCmdLine cmd("Wrinkling analysis with thin shells.");
    // Input settings
    cmd.addString( "f", "file", "Input XML file for assembler options", assemberOptionsFile );

    // Mesh settings
    cmd.addInt("r","hRefine", "Number of dyadic h-refinement (bisection) steps to perform before solving", numHref);
    cmd.addInt("e","degreeElevation", "Number of degree elevation steps to perform on the Geometry's basis before solving", numElevate);

    // Material settings
    cmd.addSwitch("TFT", "Use TFT material matrix", TFT);
    cmd.addSwitch("tftsplit", "Use the Zhang-Kiendl split tension-field material (gsMaterialMatrixTFTSplit) around the base law", tftsplit);
    cmd.addInt("", "tftlinear", "Linear flag of the TFT / Split wrapper: (1): wrapped law linear in E (default) | (0): general law", tftlinear);
    cmd.addReal("", "eta", "TFT and Split: fraction of the unmodified stress kept in wrinkled/slack states (default 0)", eta);
    cmd.addSwitch("etatangent", "TFT: apply --eta to the tangent only (stress stays the pure tension field)", etaTangent);
    cmd.addSwitch("exact", "Split: exact (non-truncated) wrinkled stress; requires --tftlinear 0", exact);
    cmd.addReal("", "exacttol", "Split: ExactTol, convergence tolerance of the exact wrinkled-stress solve (default 1e-12)", exactTol);
    cmd.addReal("", "exactstagtol", "Split: ExactStagTol, stagnation tolerance of the exact solve; must be >= exacttol (default 1e-9)", exactStagTol);
    cmd.addInt("", "exactmaxit", "Split: ExactMaxIt, maximum iterations of the exact solve (default 50)", exactMaxIt);
    cmd.addReal("", "slackmult", "TFT and Split: SlackMultiplier of the slack state (default 0)", slackmult);
    cmd.addInt("", "tensionfield", "Wrinkling criterion: 0 strain, 1 stress, 2 mixed (default 2)", tensionField);
    cmd.addSwitch("membrane", "Use membrane material", membrane);

    // Test case
    cmd.addInt("t", "test", "Test case: (0): Annulus | (1): Cylinder", testCase);

    // Load/displ stepping settings
    cmd.addInt("N", "maxsteps", "Maximum number of accepted steps (0: walk to D = 1)", step);
    cmd.addInt("", "maxbisect", "Maximum number of consecutive halvings of the load step after a failed solve", maxbisect);
    cmd.addInt("", "maxitDR", "Maximum number of Dynamic Relaxation iterations", maxitDR);
    cmd.addReal( "L", "step", "load step size", dL);
    cmd.addReal( "a", "alpha", "alpha",  alpha );
    cmd.addReal( "c", "damping", "damping",  damping );
    cmd.addSwitch("DR", "Use Dynamic Relaxation", DR);
    cmd.addSwitch("OP", "Use Optimizer", OP);
    cmd.addSwitch("NR", "Use Newton Raphson", NR);
    cmd.addSwitch("linesearch", "Newton: backtracking line search on the residual norm (gsStaticNewton option LineSearch)", lineSearch);
    cmd.addString("", "solver", "Linear solver of the Newton iterations: LU, CGDiagonal, SimplicialLDLT, ... (default CGDiagonal)", linSolver);
    cmd.addReal("", "lintol", "Relative tolerance of the iterative linear solver (gsStaticBase option SolverTol; 0: solver default, machine epsilon)", linTol);
    cmd.addInt("", "linmaxit", "Iteration cap of the iterative linear solver (SolverMaxIt; 0: solver default, 2 x size)", linMaxIt);

    // Output settings
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addSwitch("mesh", "Plot mesh?", mesh);
    cmd.addString("", "out", "Output directory (default: ArcLengthResults/<case>_-r<r>-e<e>[_TFT|_Split][_membrane])", outDir);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    if (!(eta >= 0 && eta <= 1)) { gsWarn << "--eta must be in [0,1], got " << eta << "\n"; return EXIT_FAILURE; }
    if (!(linTol >= 0)) { gsWarn << "--lintol must be non-negative, got " << linTol << "\n"; return EXIT_FAILURE; }
    if (linMaxIt < 0) { gsWarn << "--linmaxit must be non-negative, got " << linMaxIt << "\n"; return EXIT_FAILURE; }
    if (!(exactStagTol >= exactTol && exactTol > 0 && exactMaxIt >= 1))
    {
        gsWarn << "Need exactstagtol >= exacttol > 0 and exactmaxit >= 1, got " << exactStagTol << ", " << exactTol << ", " << exactMaxIt << "\n";
        return EXIT_FAILURE;
    }

    GISMO_ENSURE(!(TFT && tftsplit),"--TFT and --tftsplit are mutually exclusive");
    GISMO_ENSURE(tensionField>=0 && tensionField<=2,"--tensionfield must be 0 (strain), 1 (stress) or 2 (mixed)");
    GISMO_ENSURE(tftlinear==0 || tftlinear==1,"--tftlinear must be 0 or 1");
    GISMO_ENSURE(!(exact && tftlinear==1),"--exact requires --tftlinear 0 (exact and linearised wrinkled stress coincide for a law linear in E)");
    GISMO_ENSURE(maxbisect>=0,"--maxbisect must be non-negative");
    if (!TFT && !tftsplit && eta!=0)
        gsWarn<<"--eta is ignored without a wrinkling model\n";
    if (!TFT && etaTangent)
        gsWarn<<"--etatangent only acts with --TFT and is ignored\n";
    if (!tftsplit && exact)
        gsWarn<<"--exact only acts with --tftsplit and is ignored\n";

    std::string output = "solution";
    std::string dirname = "ArcLengthResults";

    gsFileData<> fd(assemberOptionsFile);
    gsOptionList opts;
    fd.getFirst<gsOptionList>(opts);

    real_t thickness = 0.05e-3;
    real_t YoungsModulus = 1e9;
    real_t PoissonRatio = 0.5;
    real_t Density = 1e0;

    gsMultiPatch<> mp, mp_def;
    if (testCase==0)
    {
        mp = gsNurbsCreator<>::MultiPatchAnnulus(0.0625,0.25);
        mp.embed(3);
        // MultiPatchAnnulus leaves the third coefficient column of patches 1-3 uninitialised
        for (size_t p = 0; p!=mp.nPatches(); ++p)
            mp.patch(p).coefs().col(2).setZero();
    }
    else if (testCase==1)
        mp = gsNurbsCreator<>::MultiPatchCylinder(0.25,1.0);

    mp.computeTopology();

    for (size_t p = 0; p!=mp.nPatches(); ++p)
    {
        for(index_t i = 0; i< numElevate; ++i)
            mp.patch(p).degreeElevate();    // Elevate the degree
        for(index_t i = 0; i< numHref; ++i)
            mp.patch(p).uniformRefine();
    }

    mp_def = mp;

    // Make unstructured spline
    gsMultiPatch<> geom;
    gsMultiBasis<> dbasis;
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

    // Assign BCs
    gsFunctionExpr<> dx("sqrt(x^2+y^2)*(cos(atan2(y,x) + pi/2*t)- cos(atan2(y,x)))",3);
    gsFunctionExpr<> dy("sqrt(x^2+y^2)*(sin(atan2(y,x) + pi/2*t)- sin(atan2(y,x)))",3);
    real_t dz_val = 0;
    if (testCase==0)
        dz_val = 125e-3;
    else if (testCase==1)
        dz_val = 1.0;
    gsFunctionExpr<> dz(util::to_string(dz_val) + "*1",3);

    // Boundary conditions
    gsBoundaryConditions<> BCs;
    BCs.setGeoMap(geom);
    for (size_t p = 0; p!=geom.nPatches(); ++p)
    {
        BCs.addCondition(p,boundary::south, condition_type::dirichlet, 0  , 0, false, -1);
        BCs.addCondition(p,boundary::north, condition_type::dirichlet, &dx, 0, false, 0);
        BCs.addCondition(p,boundary::north, condition_type::dirichlet, &dy, 0, false, 1);
        BCs.addCondition(p,boundary::north, condition_type::dirichlet, &dz, 0, false, 2);
    }

    if (testCase==0)
        dirname = dirname + "/Annulus_-r" + std::to_string(numHref) + "-e" + std::to_string(numElevate);
    else if (testCase==1)
        dirname = dirname + "/Cylinder_-r" + std::to_string(numHref) + "-e" + std::to_string(numElevate);
    else
        GISMO_ERROR("Test case "<<testCase<<" not implemented");
    if (TFT)
        dirname = dirname + "_TFT";
    if (tftsplit)
        dirname = dirname + "_Split";
    if (membrane)
      	dirname = dirname + "_membrane";
    if (!outDir.empty())
        dirname = outDir;

    if (!makeDirs(dirname))
    {
        gsWarn<<"Could not create the output directory "<<dirname<<"\n";
        return EXIT_FAILURE;
    }

    output =  "solution";
    wn = output + "data.txt";

    // plot geometry
    if (plot)
      gsWriteParaview(mp,dirname + "/" + "mp",1000,true);

    // Linear isotropic material model
    gsFunctionExpr<> force("0","0","0",3);

    gsConstantFunction<> t(thickness,3);
    gsConstantFunction<> E(YoungsModulus,3);
    gsConstantFunction<> nu(PoissonRatio,3);
    gsConstantFunction<> rho(Density,3);

    std::vector<gsFunctionSet<>*> parameters;
    parameters.resize(2);
    parameters[0] = &E;
    parameters[1] = &nu;

    gsMaterialMatrixBase<real_t>::uPtr materialMatrix;
    gsMaterialMatrixBase<real_t>::uPtr materialMatrixTFT;
    gsMaterialMatrixBase<real_t>::uPtr oracle;
    gsThinShellAssemblerBase<real_t>* assembler;

    gsOptionList options;
    options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",1);
    options.addSwitch("Compressibility","Compressibility: (false): Imcompressible | (true): Compressible",false);
    options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
    materialMatrix = getMaterialMatrix<3,real_t>(mp,t,parameters,rho,options);
    // The wrinkling decorators clone the base law and read the criterion from the clone
    materialMatrix->options().setInt("TensionField",tensionField);
    const int materialModel = TFT ? 1 : (tftsplit ? 2 : 0);
    if (tftlinear)
        makeWrinkleMaterials<true>(materialMatrixTFT,oracle,*materialMatrix,materialModel,eta,etaTangent,exact,slackmult,tensionField,exactTol,exactStagTol,exactMaxIt);
    else
        makeWrinkleMaterials<false>(materialMatrixTFT,oracle,*materialMatrix,materialModel,eta,etaTangent,exact,slackmult,tensionField,exactTol,exactStagTol,exactMaxIt);
    gsMaterialMatrixContainer<real_t> materialMatrixContainer;
    if (TFT || tftsplit)
    {
        for (size_t p = 0; p!=mp.nPatches(); p++)
            materialMatrixContainer.add(*materialMatrixTFT);
        assembler = new gsThinShellAssembler<3, real_t, false>(geom,dbasis,BCs,force,materialMatrixContainer);
    }
    else if (membrane)
    {
        for (size_t p = 0; p!=mp.nPatches(); p++)
            materialMatrixContainer.add(*materialMatrix);
        assembler = new gsThinShellAssembler<3, real_t, false >(geom,dbasis,BCs,force,materialMatrixContainer);
    }
    else
    {
        for (size_t p = 0; p!=mp.nPatches(); p++)
            materialMatrixContainer.add(*materialMatrix);
        assembler = new gsThinShellAssembler<3, real_t, true >(geom,dbasis,BCs,force,materialMatrixContainer);
    }


    assembler->setOptions(opts);
    assembler->options().setInt("Continuity",-1);
    assembler->setSpaceBasis(bb2);

    // Assemble linear system to obtain the force vector
    assembler->assemble();
    gsVector<> F = assembler->rhs();
    gsSparseMatrix<> K = assembler->matrix();
    assembler->assembleMass(true);
    gsVector<> M = assembler->rhs();

    // Function for the Jacobian
    gsStructuralAnalysisOps<real_t>::Jacobian_t Jacobian = [&assembler,&bb2, &geom](gsVector<real_t> const &x, gsSparseMatrix<real_t> & m)
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
    gsStructuralAnalysisOps<real_t>::Residual_t Residual = [&assembler,&bb2, &geom](gsVector<real_t> const &x, gsVector<real_t> & result)
    {
        ThinShellAssemblerStatus status;
        gsMatrix<real_t> solFull = assembler->fullSolutionVector(x);
        size_t d = geom.targetDim();
        GISMO_ASSERT(solFull.rows() % d==0,"Rows of the solution vector does not match the number of control points");
        solFull.resize(solFull.rows()/d,d);

        gsMappedSpline<2,real_t> mspline(bb2,solFull);
        gsFunctionSum<real_t> def(&geom,&mspline);

        status = assembler->assembleVector(def);
        result = assembler->rhs();
        return status == ThinShellAssemblerStatus::Success;

    };

    gsParaviewCollection<> collection(dirname + "/" + output);
    gsParaviewCollection<> collectionTF(dirname + "/" + "tensionfield");
    gsMultiPatch<> deformation = mp;

    // Define solvers
    gsStaticDR<real_t> DynamicRelaxationSolver(M,F,Residual);
    gsOptionList DynamicRelaxationOptions = DynamicRelaxationSolver.options();
    DynamicRelaxationOptions.setReal("damping",damping);
    DynamicRelaxationOptions.setReal("alpha",alpha);
    DynamicRelaxationOptions.setInt("maxIt",maxitDR);
    DynamicRelaxationOptions.setReal("tolE",1e-2);
    DynamicRelaxationOptions.setReal("tolF",1e-4);
    DynamicRelaxationOptions.setReal("tolU",1e-4);
    DynamicRelaxationOptions.setInt("verbose",1);
    DynamicRelaxationSolver.setOptions(DynamicRelaxationOptions);
    DynamicRelaxationSolver.initialize();

#ifdef gsOptim_ENABLED
    gsStaticOpt<real_t,gsOptim<real_t>::LBFGS> OptimizerSolver(Residual,assembler->numDofs());
    gsOptionList OptimizerSolverOptions = OptimizerSolver.options();
    OptimizerSolverOptions.setInt("maxIt",maxitOPT);
    OptimizerSolverOptions.setReal("tolF",1e-5);
    OptimizerSolverOptions.setReal("tolU",1e-5);
    OptimizerSolverOptions.setInt("verbose",1);
    OptimizerSolver.setOptions(OptimizerSolverOptions);
    OptimizerSolver.initialize();
#endif

    gsStaticNewton<real_t> NewtonSolver(K,F,Jacobian,Residual);
    gsOptionList NewtonOptions = NewtonSolver.options();
    NewtonOptions.setInt("verbose",true);
    NewtonOptions.setInt("maxIt",maxitNR);
    NewtonOptions.setReal("tolU",1e-6);
    NewtonOptions.setReal("tolF",1e-8);
    NewtonOptions.setString("Solver",linSolver);
    NewtonOptions.setReal("SolverTol",linTol);
    NewtonOptions.setInt("SolverMaxIt",linMaxIt);
    NewtonOptions.setSwitch("LineSearch",lineSearch);
    NewtonSolver.setOptions(NewtonOptions);
    gsInfo<<"Linear solver: "<<NewtonSolver.options().getString("Solver")
          <<" SolverTol "<<NewtonSolver.options().getReal("SolverTol")
          <<" SolverMaxIt "<<NewtonSolver.options().getInt("SolverMaxIt");
    if (tftsplit)
        gsInfo<<" | Split exact options: ExactTol "<<materialMatrixTFT->options().getReal("ExactTol")
              <<" ExactStagTol "<<materialMatrixTFT->options().getReal("ExactStagTol")
              <<" ExactMaxIt "<<materialMatrixTFT->options().getInt("ExactMaxIt");
    gsInfo<<"\n";

    // The assembler integrates clones of the Split material (one per patch), not materialMatrixTFT:
    // the counters are summed over those clones.
    std::vector<gsMaterialMatrixTFTSplit<3,real_t,true>  *> splitLin;
    std::vector<gsMaterialMatrixTFTSplit<3,real_t,false> *> splitNL;
    if (tftsplit)
        for (size_t p = 0; p!=mp.nPatches(); p++)
        {
            gsMaterialMatrixBase<real_t> * m = assembler->material(p);
            if (tftlinear)
            {
                splitLin.push_back(dynamic_cast<gsMaterialMatrixTFTSplit<3,real_t,true> *>(m));
                GISMO_ENSURE(splitLin.back()!=NULL,"The assembler material of patch "<<p<<" is not a Split material");
            }
            else
            {
                splitNL.push_back(dynamic_cast<gsMaterialMatrixTFTSplit<3,real_t,false> *>(m));
                GISMO_ENSURE(splitNL.back()!=NULL,"The assembler material of patch "<<p<<" is not a Split material");
            }
        }
    auto resetSplitCounters = [&splitLin,&splitNL]()
    {
        for (auto s : splitLin) s->resetCounters();
        for (auto s : splitNL)  s->resetCounters();
    };
    // Order: thetaFallbacks, tftFailures, exactNonconverged, exactSaddle, exactWrongBranch, exactStagnated,
    // exactReclassifiedTaut, exactReclassifiedSlack, exactSkipped, exactUsed2x2; -1 without the Split material
    auto readSplitCounters = [&splitLin,&splitNL](index_t * c)
    {
        for (index_t i = 0; i!=10; ++i) c[i] = (splitLin.empty() && splitNL.empty()) ? -1 : 0;
        auto add = [&c](const auto * s)
        {
            c[0] += s->numThetaFallbacks();         c[1] += s->numTFTFailures();
            c[2] += s->numExactNonconverged();      c[3] += s->numExactSaddle();
            c[4] += s->numExactWrongBranch();       c[5] += s->numExactStagnated();
            c[6] += s->numExactReclassifiedTaut();  c[7] += s->numExactReclassifiedSlack();
            c[8] += s->numExactSkipped();           c[9] += s->numExactUsed2x2();
        };
        for (auto s : splitLin) add(s);
        for (auto s : splitNL)  add(s);
    };

    std::vector<gsStaticBase<real_t> *> solvers;
    if (DR)
    {
      gsInfo<<"Using Dynamic Relaxation solver\n";
      solvers.push_back(&DynamicRelaxationSolver);
    }
    if (OP)
    {
#ifdef gsOptim_ENABLED
      gsInfo<<"Using Optimizer solver\n";
      solvers.push_back(&OptimizerSolver);
#else
        gsWarn<<"G+Smo is not compiled with the gsOptim module, cannot use Optimizer solver.\n";
#endif
    }
    if (NR)
    {
      gsInfo<<"Using Newton-Raphson solver\n";
      solvers.push_back(&NewtonSolver);
    }

    gsStaticComposite<real_t> StaticSolver(solvers);
    StaticSolver.options().setInt("verbose",1);
    StaticSolver.initialize();

    gsMatrix<> U = gsMatrix<>::Zero(assembler->numDofs(),1);

    // Sample points of the wrinkling-state oracle: the Gauss-Legendre nodes (quA = 1, quB = 1, i.e. p+1 per direction per element) of every patch
    gsOptionList quOpt;
    quOpt.addReal("quA","",1.0);
    quOpt.addInt("quB","",1);
    std::vector<gsMatrix<> > sampleNodes(dbasis.nBases());
    index_t numSamples = 0;
    for (size_t p = 0; p!=dbasis.nBases(); ++p)
    {
        sampleNodes[p] = gsQuadrature::getAllNodes(dbasis.basis(p),quOpt);
        numSamples += sampleNodes[p].cols();
    }
    gsMatrix<> Zmid(1,1);
    Zmid << 0;

    std::string writeName = dirname + "/out.csv";
    std::ofstream file;
    file.open(writeName);
    file<<"Load factor, Moment1, Moment2, Moment3\n";
    file<<std::setprecision(10);
    file<<0<<","<<0<<"\n";
    file.close();

    /*
        steps.csv: one row per solve attempt. Columns:
          step        index of the load step (accepted step) the attempt belongs to
          attempt     running counter over all attempts
          D           load factor of the attempt
          theta       prescribed rotation of the driven boundary, pi/2*D [rad]
          uz          prescribed axial displacement of the driven boundary (independent of D) [m]
          converged   1 if the static solver reported gsStatus::Success, 0 otherwise (also 0 if it threw)
          accepted    1 if the attempt is kept (converged, or the bisection cap was reached with a returned
                      solution), 0 if the load step is halved or the run stops
          status_dr, status_nr   static_cast<int>(gsStatus) of the solver: 0 Success, 1 NotConverged,
                      2 AssemblyError, 3 SolverError, 4 NotStarted, 5 OtherError; -1 if the solver did not
                      run (disabled, or NR after a failed DR). If the solve threw, every enabled solver is
                      recorded as 3.
          it_dr, it_nr   iterations of the solver, -1 if it did not run or threw
          bisections  cumulative number of halvings of the load step in the run, including the one
                      triggered by this attempt
          n_points    number of wrinkling-state sample points (Gauss-Legendre nodes, p+1 per direction per element, of all patches)
          n_taut, n_wrinkled, n_slack   number of sample points in the taut, wrinkled and slack state of the
                      mid-surface, by the criterion --tensionfield applied to the base law; -1 for a
                      rejected attempt or if the evaluation failed
          n_theta_fallbacks, n_tft_failures, n_exact_nonconverged, n_exact_saddle, n_exact_wrongbranch,
          n_exact_stagnated, n_exact_taut, n_exact_slack, n_exact_skipped, n_exact_used2x2
                      the counters of gsMaterialMatrixTFTSplit (numThetaFallbacks, numTFTFailures,
                      numExactNonconverged, numExactSaddle, numExactWrongBranch, numExactStagnated,
                      numExactReclassifiedTaut, numExactReclassifiedSlack, numExactSkipped, numExactUsed2x2),
                      summed over the material clones of the patches, over this attempt's whole static solve
                      (every residual and Jacobian assembly). They count material evaluations, not points.
                      nonconverged + saddle + wrongbranch = numExactFallbacks() and taut + slack =
                      numExactReclassified(); stagnated and used2x2 are independent of the outcome counters.
                      -1 for a rejected attempt or a run without --tftsplit
          M1, M2, M3  the Moment1..3 of out.csv
          Mz          torque about the z axis [N m]: internal virtual work with the rigid-rotation field
                      v = (-y, x, 0) of the current configuration on the driven side
          Fz          axial reaction force [N]: internal virtual work with v = (0, 0, 1) on the driven
                      side; positive if the boundary has to pull to hold the elongation
          (M1..Fz are nan for a rejected attempt)
          time_s      wall time of the static solve [s]
    */
    std::ofstream stepsFile(dirname + "/steps.csv");
    GISMO_ENSURE(stepsFile.is_open(),"Could not open "<<dirname<<"/steps.csv");
    stepsFile<<"step,attempt,D,theta,uz,converged,accepted,status_dr,status_nr,it_dr,it_nr,bisections,n_points,n_taut,n_wrinkled,n_slack,n_theta_fallbacks,n_tft_failures,n_exact_nonconverged,n_exact_saddle,n_exact_wrongbranch,n_exact_stagnated,n_exact_taut,n_exact_slack,n_exact_skipped,n_exact_used2x2,M1,M2,M3,Mz,Fz,time_s\n";
    stepsFile<<std::setprecision(12);

    index_t k=0;
    index_t attempt=0;
    index_t bisections=0;
    index_t stepBisections=0;
    bool runFailed = false;
    real_t lastAccepted = 0; // last displacement whose solution was accepted
    gsStopwatch stopwatch;
    std::queue<real_t> Dcheckpoints;
    Dcheckpoints.push(0.7);
    Dcheckpoints.push(0.8);
    Dcheckpoints.push(0.9);
    Dcheckpoints.push(1.0);
    real_t dL0 = dL;
    dL = std::min(dL,Dcheckpoints.front());
    real_t D   = dL;
    while (true)
    {
        gsInfo<<"Displacement step "<<k<<"; D = "<<D<<"; dL = "<<dL<<"\n";
        dx.set_t(D);
        dy.set_t(D);
        dz.set_t(D);
        assembler->updateBCs(BCs);
        StaticSolver.setDisplacement(U);

        // The solvers signal a failed residual or Jacobian evaluation by throwing
        bool thrown = false;
        resetSplitCounters();
        stopwatch.restart();
        try { StaticSolver.solve(); }
        catch (int rv) { thrown = true; gsWarn<<"Static solver threw "<<rv<<"\n"; }
        catch (const std::exception & e) { thrown = true; gsWarn<<"Static solver threw: "<<e.what()<<"\n"; }
        catch (...) { thrown = true; gsWarn<<"Static solver threw an unknown exception\n"; }
        const real_t solveTime = stopwatch.stop();
        // Read before post-processing: the Moments below evaluate the same material clones.
        index_t splitCnt[10];
        readSplitCounters(splitCnt);
        ++attempt;

        const bool converged = !thrown && StaticSolver.status() == gsStatus::Success;
        int statusDR = -1, statusNR = -1;
        index_t itDR = -1, itNR = -1;
        if (thrown)
        {
            if (DR) statusDR = static_cast<int>(gsStatus::SolverError);
            if (NR) statusNR = static_cast<int>(gsStatus::SolverError);
        }
        else
        {
            if (DR)
            {
                statusDR = static_cast<int>(DynamicRelaxationSolver.status());
                itDR = DynamicRelaxationSolver.iterations();
            }
            if (NR && (!DR || DynamicRelaxationSolver.status() == gsStatus::Success))
            {
                statusNR = static_cast<int>(NewtonSolver.status());
                itNR = NewtonSolver.iterations();
            }
        }

        // A failed load step is halved at most maxbisect times in a row; after that a returned
        // (non-converged) solution is accepted only if it is finite and has the size of the free
        // DOFs: a solve that diverged (e.g. through an inverted element) leaves no usable solution.
        const bool halve = !converged && stepBisections < maxbisect;
        const bool usable = !thrown && StaticSolver.solution().rows() == assembler->numDofs()
                            && StaticSolver.solution().allFinite();
        const bool accepted = !halve && usable;
        if (halve)
        {
            ++bisections;
            ++stepBisections;
        }

        if (!accepted)
        {
            stepsFile<<k<<","<<attempt<<","<<D<<","<<EIGEN_PI/2*D<<","<<dz_val<<","
                     <<0<<","<<0<<","<<statusDR<<","<<statusNR<<","<<itDR<<","<<itNR<<","
                     <<bisections<<","<<numSamples<<",-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,nan,nan,nan,nan,nan,"<<solveTime<<std::endl;
            if (halve)
            {
                dL = dL/2;
                D = D - dL;
                continue;
            }
            gsWarn<<"The solve failed after "<<stepBisections<<" halvings of the load step, stopping at D = "<<D<<"\n";
            runFailed = true;
            break;
        }

        U = StaticSolver.solution();
        lastAccepted = D;

        // 1. Get all the coefficients (including the ones from the eliminated BCs.)
        // gsMatrix<real_t> solFull = assembler->fullSolutionVector(controlDC.solutionU());
        gsMatrix<real_t> solFull = assembler->fullSolutionVector(U);

        // 2. Reshape all the coefficients to a Nx3 matrix
        size_t d = geom.targetDim();
        GISMO_ASSERT(solFull.rows() % d==0,"Rows of the solution vector does not match the number of control points");
        solFull.resize(solFull.rows()/d,d);

        // 3. Make the mapped spline
        gsMappedSpline<2,real_t> mspline(bb2,solFull);

        // 4. Create deformation spline
        gsFunctionSum<real_t> def(&geom,&mspline);

        if (plot)
        {
            // 5. Plot the mapped spline on the original geometry
            gsPiecewiseFunction<> displacements;
            assembler->constructStress(def,displacements,stress_type::displacement);
            std::string fileName = dirname + "/" + output + util::to_string(k);
            gsWriteParaview(def,displacements,fileName,10000,"_");
            fileName = gsFileManager::getFilename(fileName);
            for (size_t p=0; p<geom.nPatches(); ++p)
                collection.addPart(fileName + "_" + util::to_string(p) + ".vts",k,"solution",p);

            // 5. Construct stress
            gsPiecewiseFunction<> tensionFields;
            assembler->constructStress(def,tensionFields,stress_type::tension_field);
            std::string fileNameTF = dirname + "/" + "tensionfield" + util::to_string(k);
            gsWriteParaview(def,tensionFields,fileNameTF,10000,"_");
            fileNameTF = gsFileManager::getFilename(fileNameTF);
            for (size_t p=0; p<geom.nPatches(); ++p)
                collectionTF.addPart(fileNameTF + "_" + util::to_string(p) + ".vts",k,"solution",p);


            index_t N = 100;
            gsMatrix<> coords(2,N), supp, result, tmp;
            index_t dir = 1;
            std::vector<real_t> par;
            if (testCase==0)
                par = {0.5,0.7,0.9};
            else if (testCase==1)
                par = {0.5,0.75};
            else
                GISMO_ERROR("Test case "<<testCase<<" not implemented");

            for (size_t p=0; p!=par.size(); p++)
            {
                result.resize(mp.geoDim(),mp.nPatches()*N);
                for (index_t patch=0; patch!=def.nPieces(); patch++)
                {
                    supp = def.piece(patch).support();
                    coords.row(dir).setConstant( par[p]*( supp(dir,1)-supp(dir,0) )+supp(dir,0) );
                    coords.row(1-dir).setLinSpaced(N,0,1);
                    coords.row(1-dir).array() *= ( supp(1-dir,1)-supp(1-dir,0) );
                    coords.row(1-dir).array() += supp(1-dir,0);
                    def.piece(patch).eval_into(coords,tmp);
                    result.block(0,patch*N,tmp.rows(),tmp.cols()) = tmp;
                }
                result.transposeInPlace();
                writeToCSVfile(dirname + "/" + "line_D=" + std::to_string(D) + "_u=" + util::to_string(par[p]) + "_dir=" + std::to_string(dir) + ".csv",result);
            }
        }

        /////////////////////////////////////////////////////////////////////////////////
        // MOMENT COMPUTATION ///////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////
        // See: https://comet-fenics.readthedocs.io/en/latest/demo/tips_and_tricks/computing_reactions.html

        // Virtual displacement fields v on the driven side (zero elsewhere), as component expressions
        // (x, y, z) of the reference position and the load factor t = D. Per field the components
        // that are prescribed are given; the others are free and therefore zero.
        //  0 (Moment1): v = (-y, x, 0), rigid rotation about z in the reference configuration
        //  1 (Moment2): v = (dx, dy, 0) at t = D, the prescribed in-plane displacement of the driven
        //               boundary, so Moment2 = sum_i R_i . u_i is the work of the in-plane reactions on
        //               that displacement. A copied gsFunctionExpr starts with t = 0 (gsFunctionExprPrivate
        //               copy constructor, src/gsCore/gsFunctionExpr.hpp), hence the explicit set_t.
        //  2 (Moment3): v = (r sin(theta + pi/2*D), r cos(theta + pi/2*D), 0)
        //  3 (Mz)     : v = (-r sin(theta + pi/2*D), r cos(theta + pi/2*D), 0) = (-y_cur, x_cur, 0), the rigid
        //               rotation about z of the CURRENT configuration (the boundary is at angle theta + phi,
        //               phi = pi/2*D). The torque is sum_i R_i . (e_z x x_cur,i), so this is the field that
        //               gives it. For a tangential reaction tau (-sin(theta+phi), cos(theta+phi)), field 0
        //               only gives tau r cos(phi), which vanishes at D = 1.
        //  4 (Fz)     : v = (0, 0, 1), the axial reaction
        // Component 0 is stored in Mz_y_vec and component 1 in Mz_x_vec (by the field's role in the
        // boundary conditions below).
        std::vector<gsFunctionExpr<>> Mz_x_vec(4);
        std::vector<gsFunctionExpr<>> Mz_y_vec(4);
        Mz_x_vec[0] = gsFunctionExpr<>("sqrt(x^2+y^2)*cos(atan2(y,x))",3);
        Mz_x_vec[1] = dy;
        Mz_x_vec[1].set_t(D);
        Mz_x_vec[2] = gsFunctionExpr<>("sqrt(x^2+y^2)*cos(atan2(y,x)+pi/2*t)",3);
        Mz_x_vec[2].set_t(D);
        Mz_x_vec[3] = gsFunctionExpr<>("sqrt(x^2+y^2)*cos(atan2(y,x)+pi/2*t)",3);
        Mz_x_vec[3].set_t(D);

        Mz_y_vec[0] = gsFunctionExpr<>("-sqrt(x^2+y^2)*sin(atan2(y,x))",3);
        Mz_y_vec[1] = dx;
        Mz_y_vec[1].set_t(D);
        Mz_y_vec[2] = gsFunctionExpr<>("sqrt(x^2+y^2)*sin(atan2(y,x)+pi/2*t)",3);
        Mz_y_vec[2].set_t(D);
        Mz_y_vec[3] = gsFunctionExpr<>("-sqrt(x^2+y^2)*sin(atan2(y,x)+pi/2*t)",3);
        Mz_y_vec[3].set_t(D);
        gsFunctionExpr<> unitZ("1",3);

        // Prescribed function per field and component (nullptr: component not prescribed)
        gsFunctionExpr<> * vField[5][3] = {
            {&Mz_y_vec[0], &Mz_x_vec[0], nullptr},
            {&Mz_y_vec[1], &Mz_x_vec[1], nullptr},
            {&Mz_y_vec[2], &Mz_x_vec[2], nullptr},
            {&Mz_y_vec[3], &Mz_x_vec[3], nullptr},
            {nullptr,      nullptr,      &unitZ}};

        std::vector<real_t> Moments(5);

        for (index_t f=0; f!=5; f++)
        {
            gsExprAssembler<> A(1,1);
            A.setIntegrationElements(dbasis);
            gsExprEvaluator<> evA(A);
            gsExprAssembler<>::space u = A.getSpace(dbasis,3);

            // Define dummy boundary conditions.
            // These conditions are used to define the test space on the boundary and compute its values there.
            gsBoundaryConditions<> bc_dummy;
            for (index_t c=0; c!=3; c++)
                if (vField[f][c])
                    for (size_t p=0; p!=mp.nPatches(); p++)
                        bc_dummy.addCondition(p,boundary::north, condition_type::dirichlet, vField[f][c], 0 ,false,c);

            bc_dummy.setGeoMap(mp);
            // L2 projection: the basis of the smoothed interfaces is not a tensor basis
            u.setup(bc_dummy,dirichlet::l2Projection,-1);

            // Set all coefficients corresponding to the space to zero, except for the boundary conditions.
            gsMultiPatch<> geo;
            for (size_t b=0; b!=dbasis.nBases(); b++)
            {
                gsMatrix<> coefs(dbasis.basis(b).size(),mp.geoDim());
                coefs.setZero();
                for (size_t d=0; d!=mp.geoDim(); d++)
                    for (index_t k=0; k!=dbasis.basis(b).size(); k++)
                    {
                        const index_t ii = u.mapper().index(k,b,d);
                        if (u.mapper().is_free(k,b,d))
                            coefs(k,d) = 0.0;
                        else if (u.mapper().is_boundary(k,b,d))
                            coefs(k,d) = u.fixedPart().at(u.mapper().global_to_bindex(ii));
                        else
                            GISMO_ERROR("WHAT?");
                    }
                geo.addPatch(*dbasis.basis(b).makeGeometry(give(coefs)));
            }

            // Cleanup the assembler
            A.initSystem();
            // Set the geometry map and the deformation map.
            gsExprAssembler<>::geometryMap ori   = A.getMap(mp);
            gsExprAssembler<>::geometryMap deff  = A.getMap(def);
            // Register the custom solution with coefficients defined above to the assembler.
            auto v = A.getCoeff(geo);
            // Add the stress tensors based on the materials defined for the shell solver
            gsMaterialMatrixIntegrate<real_t,MaterialOutput::VectorN> S0f(materialMatrixContainer,&mp,&def);
            gsMaterialMatrixIntegrate<real_t,MaterialOutput::VectorM> S1f(materialMatrixContainer,&mp,&def);
            auto S0  = A.getCoeff(S0f);
            auto S1  = A.getCoeff(S1f);
            // Helper matrix for flexural components (see gsThinShellAssembler::assembleMatrix)
            gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","2",2);
            auto m2 = A.getCoeff(mult2t);
            // Define the components used in the variational formulation, using the solution v as the deformation
            auto N  = S0.tr();
            auto dEm= flat( jac(deff).tr() * jac(v) ) ;
            auto M  = S1.tr();
            auto dEf= -( deriv2(v,sn(deff).normalized().tr() ) + deriv2(deff,var1(v,deff) ) ) * reshape(m2,3,3);
            // Integrate the variational formulation, this gives the moments.
            Moments[f] = evA.integral((
                                    N*dEm.tr()
                                    + M*dEf.tr()
                                    )*meas(ori)
                                    );
        }

        // Wrinkling states of the base law at the mid-surface sample points of the accepted state
        index_t nTaut = -1, nWrinkled = -1, nSlack = -1;
        try
        {
            gsMaterialMatrixEval<real_t,MaterialOutput::TensionField> stateEval(oracle.get(),&def,Zmid);
            nTaut = nWrinkled = nSlack = 0;
            for (size_t p = 0; p!=sampleNodes.size(); ++p)
            {
                gsMatrix<> state;
                stateEval.piece(p).eval_into(sampleNodes[p],state);
                for (index_t i = 0; i!=state.cols(); ++i)
                {
                    if (state(0,i) > 0.5)       ++nTaut;
                    else if (state(0,i) < -0.5) ++nSlack;
                    else                        ++nWrinkled;
                }
            }
        }
        catch (const std::exception & e)
        {
            gsWarn<<"Wrinkling state evaluation failed: "<<e.what()<<"\n";
            nTaut = nWrinkled = nSlack = -1;
        }

        stepsFile<<k<<","<<attempt<<","<<D<<","<<EIGEN_PI/2*D<<","<<dz_val<<","
                 <<(converged ? 1 : 0)<<","<<1<<","<<statusDR<<","<<statusNR<<","<<itDR<<","<<itNR<<","
                 <<bisections<<","<<numSamples<<","<<nTaut<<","<<nWrinkled<<","<<nSlack<<","
                 <<splitCnt[0]<<","<<splitCnt[1]<<","<<splitCnt[2]<<","<<splitCnt[3]<<","<<splitCnt[4]<<","
                 <<splitCnt[5]<<","<<splitCnt[6]<<","<<splitCnt[7]<<","<<splitCnt[8]<<","<<splitCnt[9]<<","
                 <<Moments[0]<<","<<Moments[1]<<","<<Moments[2]<<","<<Moments[3]<<","<<Moments[4]<<","
                 <<solveTime<<std::endl;

        file.open(writeName,std::ofstream::out | std::ofstream::app);
        file<< std::setprecision(10)
            << D << ","
            << Moments[0] << ","
            << Moments[1] << ","
            << Moments[2] << "\n";
        file.close();

        if (gsClose(D,Dcheckpoints.front(),1e-8))
        {
            Dcheckpoints.pop();
            if (Dcheckpoints.empty())
                break;
        }

        if (step > 0 && k+1 >= step)
        {
            gsInfo<<"Stopping after "<<step<<" load steps at D = "<<D<<"\n";
            break;
        }

        // GO TO THE NEXT STEPs
        dL = std::min(dL0,Dcheckpoints.front()-D);
        D = D + dL;
        k++;
        stepBisections = 0;
    }

    if (plot)
    {
        collection.save();
        collectionTF.save();
    }

    delete assembler;

    // The checkpoint queue is empty exactly when the final displacement was accepted
    const bool reached = !runFailed && Dcheckpoints.empty();
    if (!reached)
        gsInfo<<"FINAL LOAD NOT REACHED: last converged "<<lastAccepted<<" of 1\n";
    return reached ? EXIT_SUCCESS : EXIT_FAILURE;
}


#else //gsUnstructuredSplines_ENABLED
int main(int argc, char *argv[])
{
    gsWarn<<"G+Smo is not compiled with the gsUnstructuredSplines module, cannot use unstructured splines for this benchmark.\n";
    return EXIT_FAILURE;
}
#endif//gsUnstructuredSplines_ENABLED
#else//gsKLShell_ENABLED
int main(int argc, char *argv[])
{
    gsWarn<<"G+Smo is not compiled with the gsKLShell, gsUnstructuredSplines and gsOptim modules.\n";
    return EXIT_FAILURE;
}
#endif//gsKLShell_ENABLED
