/** @file benchmark_HangingBlanket.cpp

    @brief Hanging blanket: a square membrane under self-weight, supported at its four corners by springs

    Reproduces the hanging blanket of Zhang & Kiendl, JMPS 205 (2025) 106331 (Z25), Sec. 5.5, Fig. 16,
    Table 3, which repeats Zhang & Kiendl, CMAME 432 (2024) 117386 (Z24), Sec. 5.4, Fig. 21, Tables 3-4.

    Setup (SI units: m, N, Pa, kg)
    - Geometry: flat square [0,L]^2 in the plane z = 0, one tensor B-spline patch of degree \a -p with
      \a --nx x \a --nx elements (Greville control points, so the flat square is reproduced exactly).
      Parametric (xi,eta) maps to (x,y) = L(xi,eta).
    - Material: linear Saint Venant-Kirchhoff, E, nu, thickness t, surface density rho (mass per reference
      area); the wrinkling model \a --model is none, tft (gsMaterialMatrixTFT, criterion
      \a --tensionfield: 0 strain, 1 stress, 2 mixed; options \a --eta, \a --etatangent) or split
      (gsMaterialMatrixTFTSplit, \a --eta).
    - Load: self-weight, a dead load of -rho g lambda e_z per unit reference area [N/m^2], with the load
      factor lambda ramped from 0 to 1 in \a -N equal steps. The thickness does not enter the force.
    - Supports: all four corners have u_z = 0. In-plane, each corner component is rigid or carries a
      linear spring of stiffness \a --k [N/m] (the spring force k u_i acts against the displacement):

        | corner                       | u_x    | u_y    |
        |------------------------------|--------|--------|
        | bottom-left  (southwest)     | rigid  | rigid  |
        | top-left     (northwest)     | rigid  | spring |
        | top-right    (northeast), A  | spring | spring |
        | bottom-right (southeast)     | spring | rigid  |

      The springs are diagonal penalty terms added to the residual (-k u_i) and to the tangent and
      linear stiffness (+k) of the four spring degrees of freedom.
    - Points (Fig. 16; the origin of the driver is the bottom-left corner, which leaves displacements
      unchanged): M = centre, parametric (0.5,0.5); A = top-right corner, (1,1); B = midpoint of the
      bottom edge, (0.5,0).

    Solver strategy: the flat, stress-free start has no out-of-plane stiffness, so Newton cannot start
    there. The first load step (and every step with \a --DRall) is therefore solved by dynamic relaxation
    (gsStaticDR) followed by Newton (gsStaticNewton); later steps by Newton alone, warm-started from the
    last converged solution. With a wrinkling model, every step that uses DR starts with the same
    DR and Newton solve with the plain SvK material (the predictor); its solution is the starting point of
    DR and Newton on the wrinkling material. The zero strain of the flat start leaves the
    wrinkling state undefined and DR on the wrinkling material does not settle there. A failed attempt is retried with the load step halved and with DR followed by
    Newton, at most 20 times.

    The mass scaling alpha of the DR (time step 1, mass M/alpha) is chosen automatically with \a -a 0:
    central differences are stable for dt^2 lambda_max(M^-1 K)/alpha <= 4, and Gershgorin bounds
    lambda_max <= max_i (sum_j |K_ij|)/M_i. With the safety factor 2,
    alpha = 0.5 max_i (sum_j |K_lin,ij|)/M_i, where K_lin is the linear SvK stiffness of the reference
    configuration plus the spring stiffness and M the lumped mass. The maximum runs over the degrees of
    freedom without a spring; the lumped mass of the four spring degrees of freedom is raised to
    (sum_j |K_lin,ij|)/lambda_max,free (fictitious DR mass, no influence on the static solution), because
    the spring stiffness would otherwise dictate a very large alpha, hence very slow DR, for the whole mesh.
    Wrinkling only softens the base material and the geometric stiffness is small, so this bounds the
    stiffness of the run. \a -a > 0 overrides alpha.

    Source ambiguities (the defaults follow the printed values unless stated)
    - (i) Thickness: Z25 states 1.77 mm, Z24 states 1.177 mm. The reference rows of Z25 Table 3 and Z24
      Table 4 (nu = 0.3) are identical, so one of the two is a typo. The default is 1.77e-3 m;
      \a --thickness switches.
    - (ii) Spring: printed as 22.95 kN/m (penalty alpha = 1e-2 after Herrema et al. 2019). The corner
      reactions are of order 0.3-1 N (total weight 1.41 N), so k = 22950 N/m gives |u_x^A| of order
      1e-5 m, against 0.034 m in Table 3; k = 22.95 N/m is consistent with the table. The default is the
      printed value; \a --k switches.
    - (iii) Stress unit: the tables print sigma_1^M in "MPa", but Z24 Fig. 22 is labelled Pa and
      N/t of about 1 N/m / 1.77e-3 m is about 600 Pa, so the table values are Pa.
    - (iv) Density unit: Z25 writes 0.144 kg mm^-2, Z24 0.144 kg m^-2. kg/m^2 is used.

    Reference values of Z25 Table 3 (nu = 0.3, 25 x 25 elements), u_z^M, u_x^A, u_x^B [m], sigma_1^M [Pa]:

        | row                           | u_z^M   | u_x^A    | u_x^B    | sigma_1^M |
        |-------------------------------|---------|----------|----------|-----------|
        | Reference (Nakashino-Natori)  | -0.2833 | -0.03406 | -0.01703 | 642.66    |
        | Mixed                         | -0.2833 | -0.03406 | -0.01703 | 623.95    |
        | Strain                        | -0.2956 | -0.03290 | -0.01645 | 586.14    |
        | Stress                        | -0.2887 | -0.03365 | -0.01683 | 621.34    |

    Output (directory \a --out, default ./HangingBlanket_<model>)

    steps.csv: one row per attempted load step (std::setprecision(12)); columns
    - step:      index of the load step being solved (1-based)
    - attempt:   running index of the attempt
    - load:      load factor lambda of the attempt
    - converged: 1 if the solver chain converged, else 0
    - itDR, itNR: iterations of DR and Newton on the run material, -1 if that solver did not run (the SvK
      predictor is not counted)
    - time:      solve time [s]
    - uz_M:      z-displacement at M [m]
    - ux_A, uy_A: x- and y-displacement at A [m]
    - ux_B:      x-displacement at B [m]
    - S1_M, S2_M: principal mid-surface stresses at M, S1 >= S2 [Pa]. These are principal values of the
      2nd Piola-Kirchhoff stress (MaterialOutput::PStress), which differ from the Cauchy values at the
      order of the strain.
    The QoI columns are NaN for non-converged attempts.

    supports.csv: one row per converged load step; columns
    - step, load
    - for each spring degree of freedom, in the order NW_y, NE_x, NE_y, SE_x:
      u_* (displacement [m]), ku_* (spring force k u [N]), Rm_* (membrane residual f_ext - f_int at the
      converged solution without the spring term [N]) and relerr_* = |Rm - k u| / |k u|
    - rigid_max: max |u| over the 8 constrained corner components [m]

    Command line: all options have defaults; run without arguments for the default 25 x 25, p = 2 case.
    The exit code is 1 if the run stops before load factor 1 (more than 20 successive halvings), else 0.
    A DR solve that exceeds \a --maxitDR iterations counts as failed and triggers the halved-step retry.

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
#include <gsKLShell/src/gsMaterialMatrixEval.h>
#endif

#include <gsStructuralAnalysis/src/gsStaticSolvers/gsStaticDR.h>
#include <gsStructuralAnalysis/src/gsStaticSolvers/gsStaticNewton.h>
#include <gsStructuralAnalysis/src/gsStaticSolvers/gsStaticComposite.h>

#include <gsStructuralAnalysis/src/gsStructuralAnalysisTools/gsStructuralAnalysisUtils.h>

using namespace gismo;

#ifdef gsKLShell_ENABLED
namespace
{

typedef gsMaterialMatrixBase<real_t> MaterialBase;

/// Evaluates a material output of \a mat on \a def at the parametric points \a nodes and the
/// through-thickness coordinates \a Z (one column, normalised, +-0.5 = outer fibres, 0 = mid-surface).
/// Column j*nodes.cols()+k of the result belongs to fibre j and point k.
template <MaterialOutput out>
gsMatrix<> evalMaterial(MaterialBase * mat, const gsMultiPatch<> & def, const gsMatrix<> & nodes, const gsMatrix<> & Z)
{
    gsMaterialMatrixEval<real_t,out> ev(mat,&def,Z);
    gsMatrix<> result;
    ev.piece(0).eval_into(nodes,result);
    return result;
}

const real_t nanValue = std::numeric_limits<real_t>::quiet_NaN();

}// namespace

int main (int argc, char** argv)
{
    // Input options
    std::string modelName("none");
    real_t eta = 0;
    bool etaTangent = false;
    index_t tensionfield = 2;
    real_t slackmult = 1e-6;
    real_t E_modulus = 30000;
    real_t PoissonRatio = 0.3;
    real_t thickness = 1.77e-3;
    real_t rho = 0.144;
    real_t grav = 9.81;
    real_t L = 1.0;
    index_t nx = 25;
    index_t degree = 2;
    real_t kSpring = 22950;
    index_t loadSteps = 10;
    real_t alpha = 0;
    real_t damping = 0;
    real_t tolDR = 1e-1;
    real_t tolE = 1e-5;
    real_t tolNR = 1e-4;
    index_t maxitNR = 100;
    index_t maxitDR = 100000;
    std::string linSolver("LU");
    bool lineSearch = false;
    bool DRall = false;
    std::string outDir;
    bool plot = false;
    index_t verbose = 0;

    gsCmdLine cmd("Hanging blanket: square membrane under self-weight on corner springs (Zhang & Kiendl 2025, Sec. 5.5).");
    cmd.addString("","model","Wrinkling model: none, tft or split (default none)",modelName);
    cmd.addReal("","eta","TFT and Split: stress fraction retained in the wrinkled/slack state (default 0; needs --model tft or split)",eta);
    cmd.addSwitch("etatangent","TFT: apply --eta to the tangent only (stress stays the pure tension field)",etaTangent);
    cmd.addInt("","tensionfield","Wrinkling criterion: 0 strain, 1 stress, 2 mixed (default 2)",tensionfield);
    cmd.addReal("","slackmult","Multiplier of the stiffness in the slack state, TFT and Split (default 1e-6)",slackmult);
    cmd.addReal("E","young","Young's modulus [Pa] (default 30000)",E_modulus);
    cmd.addReal("","nu","Poisson's ratio (default 0.3)",PoissonRatio);
    cmd.addReal("","thickness","Thickness t [m] (default 1.77e-3)",thickness);
    cmd.addReal("","rho","Surface density, mass per reference area [kg/m^2] (default 0.144)",rho);
    cmd.addReal("g","gravity","Gravity [m/s^2] (default 9.81)",grav);
    cmd.addReal("","L","Side length [m] (default 1)",L);
    cmd.addInt("","nx","Elements per direction (default 25)",nx);
    cmd.addInt("p","degree","Spline degree (default 2)",degree);
    cmd.addReal("","k","Spring stiffness [N/m] (default 22950)",kSpring);
    cmd.addInt("N","loadSteps","Number of load steps (default 10)",loadSteps);
    cmd.addReal("a","alpha","DR mass scaling alpha; 0: automatic, see the file documentation (default 0)",alpha);
    cmd.addReal("c","damping","DR damping (default 0)",damping);
    cmd.addReal("","tolDR","DR residual tolerance (default 1e-1)",tolDR);
    cmd.addReal("","tolE","DR kinetic-energy tolerance (default 1e-5)",tolE);
    cmd.addReal("","tolNR","Newton tolerance (default 1e-4)",tolNR);
    cmd.addInt("","maxitNR","Maximum number of Newton iterations (default 100)",maxitNR);
    cmd.addInt("","maxitDR","Maximum number of DR iterations; a DR solve that exceeds it counts as failed and the load step is halved (default 100000)",maxitDR);
    cmd.addSwitch("linesearch", "Newton: backtracking line search on the residual norm (gsStaticNewton option LineSearch)", lineSearch);
    cmd.addString("","solver","Linear solver of the Newton iterations: LU, CGDiagonal, SimplicialLDLT, ... (default LU)",linSolver);
    cmd.addSwitch("DRall","DR followed by Newton on every load step (default: DR on the first step only)",DRall);
    cmd.addString("","out","Output directory (default ./HangingBlanket_<model>)",outDir);
    cmd.addSwitch("plot","Plot the deformed blanket in ParaView format",plot);
    cmd.addInt("v","verbose","Solver verbosity (default 0)",verbose);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    if (!(eta >= 0 && eta <= 1)) { gsWarn << "--eta must be in [0,1], got " << eta << "\n"; return EXIT_FAILURE; }

    int model = -1;
    if      (modelName=="none")  model = 0;
    else if (modelName=="tft")   model = 1;
    else if (modelName=="split") model = 2;
    else GISMO_ERROR("Unknown --model "<<modelName<<"; use none, tft or split");
    if (eta!=0 && model==0)
        gsWarn<<"--eta is ignored without a wrinkling model\n";
    if (etaTangent && model!=1)
        gsWarn<<"--etatangent only acts together with --model tft and is ignored\n";
    GISMO_ENSURE(nx>=1 && degree>=1 && loadSteps>=1,"--nx, -p and -N must be positive");
    GISMO_ENSURE(maxitDR>=2,"--maxitDR must be at least 2");

    const index_t d = 3;

    // Single patch with the Greville points of the spline basis as control points; the flat square is reproduced exactly
    gsKnotVector<> kx(0,1,nx-1,degree+1), ky(0,1,nx-1,degree+1);
    gsTensorBSplineBasis<2,real_t> tb(kx,ky);
    const gsMatrix<> anchors = tb.anchors();
    gsMatrix<> coefs(tb.size(),d);
    coefs.setZero();
    for (index_t i = 0; i!=coefs.rows(); ++i)
    {
        coefs(i,0) = L*anchors(0,i);
        coefs(i,1) = L*anchors(1,i);
    }
    gsMultiPatch<> mp, mp_def;
    mp.addPatch(gsTensorBSpline<2,real_t>(tb,coefs));
    mp_def = mp;
    gsMultiBasis<> dbasis(mp,true);
    gsInfo<<"Hanging blanket, model "<<modelName<<", basis (patch 0): "<<mp.patch(0).basis()<<"\n";

    // Boundary conditions: u_z = 0 at all four corners; in-plane the corners are rigid or sprung (see the file documentation)
    gsBoundaryConditions<> BCs;
    BCs.setGeoMap(mp);
    BCs.addCornerValue(boundary::southwest,0.0,0,0,2);
    BCs.addCornerValue(boundary::northwest,0.0,0,0,2);
    BCs.addCornerValue(boundary::northeast,0.0,0,0,2);
    BCs.addCornerValue(boundary::southeast,0.0,0,0,2);
    BCs.addCornerValue(boundary::southwest,0.0,0,0,0);
    BCs.addCornerValue(boundary::southwest,0.0,0,0,1);
    BCs.addCornerValue(boundary::northwest,0.0,0,0,0);
    BCs.addCornerValue(boundary::southeast,0.0,0,0,1);

    // Material: linear SvK base, optionally wrapped by a wrinkling model. The wrappers do not own the base, and
    // TFT and Split read the TensionField criterion from the base and from their own option list, respectively.
    gsConstantFunction<> t(thickness,d);
    gsConstantFunction<> E(E_modulus,d);
    gsConstantFunction<> nu(PoissonRatio,d);
    gsConstantFunction<> rhoVol(rho/thickness,d);

    std::vector<gsFunctionSet<>*> parameters;
    parameters.push_back(&E);
    parameters.push_back(&nu);

    gsOptionList options;
    options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
    options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);

    MaterialBase::uPtr base = getMaterialMatrix<3,real_t>(mp,t,parameters,rhoVol,options);
    base->options().setInt("TensionField",tensionfield);
    MaterialBase::uPtr wrapped;
    if (model==1)
    {
        wrapped = memory::make_unique(new gsMaterialMatrixTFT<3,real_t,true>(base.get()));
        wrapped->options().setReal("SlackMultiplier",slackmult);
        wrapped->options().setReal("Eta",eta);
        wrapped->options().setSwitch("EtaTangentOnly",etaTangent);
    }
    else if (model==2)
    {
        wrapped = memory::make_unique(new gsMaterialMatrixTFTSplit<3,real_t,true>(base.get()));
        gsOptionList o = wrapped->options();
        o.setReal("Eta",eta);
        o.setSwitch("Exact",false);
        o.setSwitch("Midsurface",false);
        o.setReal("S2crit",0);
        o.setReal("SlackMultiplier",slackmult);
        o.setInt("TensionField",tensionfield);
        wrapped->setOptions(o);
    }
    MaterialBase * runMaterial = wrapped ? wrapped.get() : base.get();

    // Self-weight: force per unit reference area [N/m^2], updated in every load step
    gsVector<> forceVec(d);
    forceVec.setZero();
    gsConstantFunction<> force(forceVec,d);

    std::unique_ptr<gsThinShellAssemblerBase<real_t> > assembler(new gsThinShellAssembler<3,real_t,false>(mp,dbasis,BCs,force,runMaterial));
    GISMO_ENSURE(assembler->material(0) == runMaterial,"The assembler does not use the expected material object");

    // Spring degrees of freedom (free-vector indices): the corner control points interpolate and the reference is
    // flat, so the free-vector entry is the corner displacement component
    const gsDofMapper mapper = assembler->getMapper();
    const gsBasis<> & b0 = dbasis.basis(0);
    const index_t cSW = b0.functionAtCorner(boxCorner(boundary::southwest));
    const index_t cSE = b0.functionAtCorner(boxCorner(boundary::southeast));
    const index_t cNW = b0.functionAtCorner(boxCorner(boundary::northwest));
    const index_t cNE = b0.functionAtCorner(boxCorner(boundary::northeast));
    const index_t springCorner[4] = {cNW,cNE,cNE,cSE};
    const index_t springComp[4]   = {1,0,1,0};
    const char * springName[4]    = {"NW_y","NE_x","NE_y","SE_x"};
    std::vector<index_t> springDofs(4);
    for (index_t s = 0; s!=4; ++s)
    {
        GISMO_ENSURE(mapper.is_free(springCorner[s],0,springComp[s]),"Spring DOF "<<springName[s]<<" is not free");
        springDofs[s] = mapper.index(springCorner[s],0,springComp[s]);
    }
    const index_t rigidCorner[8] = {cSW,cSW,cNW,cSE,cSW,cNW,cNE,cSE};
    const index_t rigidComp[8]   = {0,1,0,1,2,2,2,2};
    for (index_t s = 0; s!=8; ++s)
        GISMO_ENSURE(!mapper.is_free(rigidCorner[s],0,rigidComp[s]),"A constrained corner component is free");

    // Output
    std::string dirname = outDir.empty() ? "./HangingBlanket_" + modelName : outDir;
    for (size_t pos = dirname.find('/',1); pos != std::string::npos; pos = dirname.find('/',pos+1))
        gsFileManager::mkdir(dirname.substr(0,pos));
    GISMO_ENSURE(gsFileManager::mkdir(dirname),"Cannot create the output directory "<<dirname);

    // Assembly at the reference configuration (zero load) for the matrix structure
    gsStopwatch stopwatch;
    assembler->assemble();
    gsSparseMatrix<> K = assembler->matrix();
    for (index_t s = 0; s!=4; ++s)
        K.coeffRef(springDofs[s],springDofs[s]) += kSpring;
    gsVector<> F = assembler->rhs();
    assembler->assembleMass(true);
    gsVector<> M = assembler->rhs();
    gsInfo<<"Initial assembly: "<<stopwatch.stop()<<" s, "<<assembler->numDofs()<<" DOFs\n";

    // Gershgorin bound of the row stiffness over lumped mass of the linear SvK stiffness plus springs. The four stiff
    // spring rows would force a huge alpha for the whole mesh, so the lumped mass of those rows is raised (fictitious
    // DR mass, which does not change the static solution) to the ratio of the remaining rows, and alpha follows
    // the remaining rows: alpha = 0.5 max_i sum_j |K_lin,ij| / M_i.
    gsVector<> rowSum;
    gsThinShellAssembler<3,real_t,false> baseAssembler(mp,dbasis,BCs,force,base.get());
    {
        baseAssembler.assemble();
        gsSparseMatrix<> Klin = baseAssembler.matrix();
        for (index_t s = 0; s!=4; ++s)
            Klin.coeffRef(springDofs[s],springDofs[s]) += kSpring;
        rowSum = gsVector<>::Zero(Klin.rows());
        for (index_t j = 0; j!=Klin.outerSize(); ++j)
            for (gsSparseMatrix<>::InnerIterator it(Klin,j); it; ++it)
                rowSum(it.row()) += math::abs(it.value());
    }
    std::vector<bool> isSpring(rowSum.size(),false);
    for (index_t s = 0; s!=4; ++s)
        isSpring[springDofs[s]] = true;
    real_t lambdaFree = 0;
    for (index_t i = 0; i!=rowSum.size(); ++i)
        if (!isSpring[i])
            lambdaFree = std::max(lambdaFree,rowSum(i)/M(i));
    GISMO_ENSURE(lambdaFree>0,"Automatic DR mass scaling needs at least one free degree of freedom without a spring; use a finer mesh");
    auto balanceSpringMass = [&]()
    {
        for (index_t s = 0; s!=4; ++s)
            M(springDofs[s]) = std::max(M(springDofs[s]),rowSum(springDofs[s])/lambdaFree);
    };
    balanceSpringMass();
    if (alpha <= 0)
        alpha = 0.5*lambdaFree;
    gsInfo<<"DR mass scaling alpha = "<<alpha<<"\n";

    // Residual f_ext - f_int minus the spring forces; Jacobian plus the spring stiffness
    typedef gsStructuralAnalysisOps<real_t>::Jacobian_t Jacobian_t;
    typedef gsStructuralAnalysisOps<real_t>::Residual_t Residual_t;
    auto makeJacobian = [springDofs,kSpring](gsThinShellAssemblerBase<real_t> * a, gsMultiPatch<> * def)
    {
        return Jacobian_t([a,def,springDofs,kSpring](gsVector<real_t> const &x, gsSparseMatrix<real_t> & m)
        {
            a->constructSolution(x,*def);
            ThinShellAssemblerStatus status = a->assembleMatrix(*def);
            m = a->matrix();
            for (size_t s = 0; s!=springDofs.size(); ++s)
                m.coeffRef(springDofs[s],springDofs[s]) += kSpring;
            return status == ThinShellAssemblerStatus::Success;
        });
    };
    auto makeResidual = [springDofs,kSpring](gsThinShellAssemblerBase<real_t> * a, gsMultiPatch<> * def)
    {
        return Residual_t([a,def,springDofs,kSpring](gsVector<real_t> const &x, gsVector<real_t> & result)
        {
            a->constructSolution(x,*def);
            ThinShellAssemblerStatus status = a->assembleVector(*def);
            result = a->rhs();
            for (size_t s = 0; s!=springDofs.size(); ++s)
                result(springDofs[s]) -= kSpring*x(springDofs[s]);
            return status == ThinShellAssemblerStatus::Success;
        });
    };
    Jacobian_t Jacobian = makeJacobian(assembler.get(),&mp_def);
    Residual_t Residual = makeResidual(assembler.get(),&mp_def);
    gsMultiPatch<> mp_pred = mp;
    Jacobian_t JacobianPred = makeJacobian(&baseAssembler,&mp_pred);
    Residual_t ResidualPred = makeResidual(&baseAssembler,&mp_pred);

    gsVector<> solVector(assembler->numDofs());
    solVector.setZero();

    gsStaticDR<real_t> DRM(M,F,Residual);
    gsOptionList DROptions = DRM.options();
    DROptions.setReal("damping",damping);
    DROptions.setReal("alpha",alpha);
    DROptions.setInt("maxIt",maxitDR);
    DROptions.setReal("tol",tolDR);
    DROptions.setReal("tolE",tolE);
    DROptions.setInt("verbose",verbose);
    DROptions.setInt("ResetIt",(index_t)(100));
    DRM.setOptions(DROptions);

    gsStaticNewton<real_t> NWT(K,F,Jacobian,Residual);
    gsOptionList NWTOptions = NWT.options();
    NWTOptions.setInt("maxIt",maxitNR);
    NWTOptions.setReal("tol",tolNR);
    NWTOptions.setString("Solver",linSolver);
    NWTOptions.setSwitch("LineSearch",lineSearch);
    NWTOptions.setInt("verbose",verbose);
    NWT.setOptions(NWTOptions);

    std::vector<gsStaticBase<real_t> *> solvers;
    solvers.push_back(&DRM);
    solvers.push_back(&NWT);
    gsStaticComposite<real_t> solver(solvers);
    solver.options().setInt("verbose",verbose);

    std::vector<gsStaticBase<real_t> *> solversNROnly;
    solversNROnly.push_back(&NWT);
    gsStaticComposite<real_t> solverNROnly(solversNROnly);
    solverNROnly.options().setInt("verbose",verbose);

    // Predictor for the wrinkling models: the flat, stress-free start has zero membrane strain, where the wrinkling state
    // (taut, wrinkled, slack) is undefined and DR does not settle. The base SvK problem is solved first (DR, then Newton)
    // and its solution is the starting point of the wrinkling problem.
    gsSparseMatrix<> Kpred;
    gsStaticDR<real_t> DRP(M,F,ResidualPred);
    DRP.setOptions(DROptions);
    gsStaticNewton<real_t> NWTP(Kpred,F,JacobianPred,ResidualPred);
    NWTP.setOptions(NWTOptions);
    std::vector<gsStaticBase<real_t> *> predSolvers;
    predSolvers.push_back(&DRP);
    predSolvers.push_back(&NWTP);
    gsStaticComposite<real_t> solverPred(predSolvers);
    solverPred.options().setInt("verbose",verbose);

    // Points M, A, B and the four corners (SW, SE, NW, NE)
    gsMatrix<> pts(2,3);
    pts.col(0)<<0.5,0.5;
    pts.col(1)<<1,1;
    pts.col(2)<<0.5,0;
    gsMatrix<> cornerPts(2,4);
    cornerPts.col(0)<<0,0;
    cornerPts.col(1)<<1,0;
    cornerPts.col(2)<<0,1;
    cornerPts.col(3)<<1,1;
    const gsMatrix<> Z0 = gsMatrix<>::Zero(1,1);
    gsMatrix<> ptM(2,1);
    ptM<<0.5,0.5;

    std::ofstream stepsFile(dirname + "/steps.csv");
    stepsFile<<std::setprecision(12)<<"step,attempt,load,converged,itDR,itNR,time,uz_M,ux_A,uy_A,ux_B,S1_M,S2_M"<<std::endl;
    std::ofstream supportsFile(dirname + "/supports.csv");
    supportsFile<<std::setprecision(12)<<"step,load";
    for (index_t s = 0; s!=4; ++s)
        supportsFile<<",u_"<<springName[s]<<",ku_"<<springName[s]<<",Rm_"<<springName[s]<<",relerr_"<<springName[s];
    supportsFile<<",rigid_max"<<std::endl;

    gsParaviewCollection<> solcollection(dirname + "/solution");

    // Load stepping with bisection on failure (exit code 1 if the load factor 1 is not reached); DR is used on the first step, on every step with --DRall,
    // and on the retry of a failed attempt
    real_t load = 0;
    real_t dload = 1./loadSteps;
    const real_t dload0 = dload;
    index_t step = 1, attempt = 0, bisected = 0;
    bool retry = false, failed = false;
    real_t totalSolveTime = 0;
    gsStopwatch total;
    while (load < 1 - 1e-12)
    {
        const real_t target = std::min(1.0,load + dload);
        forceVec(2) = -rho*grav*target;
        force.setValue(forceVec,d);

        gsInfo<<"Load step "<<step<<": load factor = "<<target<<" (step "<<target-load<<")\n";
        assembler->assemble();
        K = assembler->matrix();
        for (index_t s = 0; s!=4; ++s)
            K.coeffRef(springDofs[s],springDofs[s]) += kSpring;
        F = assembler->rhs();
        assembler->assembleMass(true);
        M = assembler->rhs();
        balanceSpringMass();

        const bool useDR = step==1 || DRall || retry;
        gsVector<> startVector = solVector;
        if (useDR && model!=0)
        {
            baseAssembler.assemble();
            Kpred = baseAssembler.matrix();
            for (index_t s = 0; s!=4; ++s)
                Kpred.coeffRef(springDofs[s],springDofs[s]) += kSpring;
            solverPred.initialize();
            solverPred.setDisplacement(solVector);
            bool predThrew = false;
            stopwatch.restart();
            try { solverPred.solve(); } catch (int) { predThrew = true; }
            gsInfo<<"  SvK predictor: "<<(!predThrew && solverPred.converged() ? "converged" : "NOT converged")
                  <<", DR its "<<DRP.iterations()<<", solve time "<<stopwatch.stop()<<" s\n";
            if (!predThrew && solverPred.converged())
                startVector = solverPred.solution();
        }
        gsStaticComposite<real_t> & activeSolver = useDR ? solver : solverNROnly;
        activeSolver.initialize();
        activeSolver.setDisplacement(startVector);

        // gsStaticDR::solve does not catch its error codes (1: maximum iterations reached, 2: assembly error),
        // whereas gsStaticNewton::solve reports its failures through its status
        bool drThrew = false;
        stopwatch.restart();
        try
        {
            activeSolver.solve();
        }
        catch (int errorCode)
        {
            drThrew = true;
        }
        const real_t solveTime = stopwatch.stop();
        totalSolveTime += solveTime;
        const bool converged = !drThrew && activeSolver.converged();
        const bool ranNR = !useDR || DRM.converged();
        const index_t itDR = useDR ? DRM.iterations() : -1;
        const index_t itNR = ranNR ? NWT.iterations() : -1;
        gsInfo<<"  "<<(converged ? "converged" : "NOT converged")<<", DR its "<<itDR<<", NR its "<<itNR<<", solve time "<<solveTime<<" s\n";

        if (!converged)
        {
            stepsFile<<step<<","<<attempt++<<","<<target<<",0,"<<itDR<<","<<itNR<<","<<solveTime;
            for (index_t k = 0; k!=6; ++k)
                stepsFile<<","<<nanValue;
            stepsFile<<std::endl;

            if (++bisected > 20)
            {
                gsWarn<<"Simulation terminated because of too many failures\n";
                failed = true;
                break;
            }
            gsWarn<<"Load step "<<step<<" did not converge; halving the load step\n";
            dload = (target-load)/2;
            retry = true;
            continue;
        }

        solVector = activeSolver.solution();
        mp_def = assembler->constructSolution(solVector);
        gsMultiPatch<> deformation = mp_def;
        deformation.patch(0).coefs() -= mp.patch(0).coefs();

        const gsMatrix<> disp = deformation.patch(0).eval(pts);
        const gsMatrix<> SM = evalMaterial<MaterialOutput::PStress>(runMaterial,mp_def,ptM,Z0);
        const real_t S1 = std::max(SM(0,0),SM(1,0));
        const real_t S2 = std::min(SM(0,0),SM(1,0));
        const real_t uzM = disp(2,0), uxA = disp(0,1), uyA = disp(1,1), uxB = disp(0,2);
        stepsFile<<step<<","<<attempt++<<","<<target<<",1,"<<itDR<<","<<itNR<<","<<solveTime
                 <<","<<uzM<<","<<uxA<<","<<uyA<<","<<uxB<<","<<S1<<","<<S2<<std::endl;
        gsInfo<<"  uz_M = "<<uzM<<", ux_A = "<<uxA<<", uy_A = "<<uyA<<", ux_B = "<<uxB<<", S1_M = "<<S1<<", S2_M = "<<S2<<"\n";

        // Spring reactions: the membrane residual at the converged state must equal the spring force
        assembler->assembleVector(mp_def);
        const gsVector<> Rm = assembler->rhs();
        supportsFile<<step<<","<<target;
        for (index_t s = 0; s!=4; ++s)
        {
            const real_t u = solVector(springDofs[s]);
            const real_t ku = kSpring*u;
            supportsFile<<","<<u<<","<<ku<<","<<Rm(springDofs[s])<<","<<math::abs(Rm(springDofs[s])-ku)/math::abs(ku);
        }
        const gsMatrix<> cd = deformation.patch(0).eval(cornerPts);
        real_t rigidMax = 0;
        for (index_t s = 0; s!=8; ++s)
        {
            const index_t c = (rigidCorner[s]==cSW) ? 0 : (rigidCorner[s]==cSE) ? 1 : (rigidCorner[s]==cNW) ? 2 : 3;
            rigidMax = std::max(rigidMax,math::abs(cd(rigidComp[s],c)));
        }
        supportsFile<<","<<rigidMax<<std::endl;

        if (plot)
        {
            gsField<> solField(mp_def,deformation);
            const std::string fileName = dirname + "/solution" + util::to_string(step);
            gsWriteParaview<>(solField,fileName,1000,true);
            solcollection.addPart(gsFileManager::getFilename(fileName) + "0.vts",step);
        }

        load = target;
        step++;
        bisected = 0;
        retry = false;
        dload = dload0;
    }

    solcollection.save();
    gsInfo<<"Reached load factor "<<load<<" in "<<step-1<<" converged steps and "<<attempt<<" attempts; solve time "
          <<totalSolveTime<<" s, total time "<<total.stop()<<" s\nResults in "<<dirname<<"\n";

    const bool reached = !failed && load >= 1 - 1e-12;
    if (!reached)
        gsInfo<<"FINAL LOAD NOT REACHED: last converged "<<load<<" of 1\n";
    return reached ? EXIT_SUCCESS : EXIT_FAILURE;
}
#else//gsKLShell_ENABLED
int main(int argc, char *argv[])
{
    gsWarn<<"G+Smo is not compiled with the gsKLShell module.\n";
    return 0;
}
#endif
