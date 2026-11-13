/** @file benchmark_SteinHedgepeth.cpp

    @brief In-plane pure bending of a prestressed rectangular membrane (Stein & Hedgepeth 1961)
           with the tension-field and split wrinkling formulations under load control (Newton)

    Zhang & Kiendl (2025), Sec. 5.1: a rectangular membrane L x H = 2 x 1, thickness t = 0.01,
    E = 100, nu = 0.3, is loaded by the uniform axial prestress sigma0 = 5e-4 and, on top of it, by
    a bending moment M applied through a linearly varying edge stress. The units are consistent
    and dimensionless, as in the paper. The material is a linear (St. Venant-Kirchhoff) plane
    stress membrane, optionally wrapped by a wrinkling model (\a --model: \a none, \a tft for
    gsMaterialMatrixTFT, \a split for gsMaterialMatrixTFTSplit).

    Setup. The right half [0,1] x [0,1] of the membrane is modelled with symmetry at x = 0
    (u_x = 0 there). The mesh has 11 x 5 bi-quadratic elements (Nakashino & Natori 2005) and
    consists of two patches coupled strongly (C^0, option "Continuity" = 0 of the assembler; the
    default weak C^1 coupling does not act on planar problems):
    - patch 0: [0,10/11] x [0,1], 10 x 5 elements, with the wrinkling model;
    - patch 1: [10/11,1] x [0,1], 1 x 5 elements, the taut cap with the plain SvK material.
    The degrees of freedom are those of 98 control points with two components, minus the 7
    constrained u_x of the symmetry edge: 189.

    Loads (dead, force per unit undeformed edge length). With P = sigma0 t H the axial load and
    r = 2 M/(P H) the moment ratio, the equivalent edge stress is
    sigma_M(y) = 6 M/(t H^2) (2y/H - 1) = 3 r sigma0 (2y/H - 1). Fig. 5(b) of the paper is loaded by
    - top edge y = H (both patches):    (0, +sigma0 t lambda_P)
    - bottom edge y = 0 (both patches): (0, -sigma0 t lambda_P)
    - right edge x = L/2 (patch 1):     (sigma0 t lambda_P + 3 sigma0 t r (2y/H - 1), 0),
    i.e. (sigma_p + sigma_M(y)) t with sigma_p = sigma0 lambda_P. The load coordinate t of the
    expressions is lambda_P in phase 1 and 1 + r in phase 2.
    The vertical edge loads are self-equilibrated, so the y-translation of the half model is a null
    mode. It is removed with a linear spring of stiffness kSpring = 1e3 E t on the control point
    (0, H/2): E t is the membrane stiffness scale per degree of freedom, so the spring is stiff
    against the structure while its force stays ~0 (reported as F_spring); a larger value would only
    worsen the conditioning.

    Load schedule. Phase 1 raises lambda_P from 0 to 1 (prestress and P) in \a --presteps steps.
    Phase 2 raises r from 0 to \a --mmax in \a -N nominal steps at constant lambda_P = 1; the
    values r in {0.4, 0.6, 0.75, 0.825} below \a --mmax are hit exactly, as is \a --mmax. A step
    that does not converge is halved, at most \a --maxbisect times in a row; the run then stops
    and the exit code is 1.

    Options: \a --model, \a --eta (TFT and Split), \a --etatangent (TFT: \a --eta in the tangent only), \a --tensionfield (wrinkling criterion: 0 strain-based,
    1 stress-based, 2 mixed, the default), \a --slackmult, \a --nu, \a --presteps, \a -N (steps of
    phase 2), \a --mmax (final r, < 1), \a --tolNR, \a --maxitNR, \a --solver, \a --linesearch, \a --maxbisect,
    \a --out, \a -v.

    Output in \a <out> (default ./SteinHedgepeth_<model>):
    - steps.csv, one row per attempted step. Columns: step, attempt, phase (1 or 2), D (load
      coordinate of the phase), lambdaP, r, itNR (Newton iterations), converged, halvings
      (successive halvings before this attempt), Rx_west, Ry_west (boundaryForce of the symmetry
      edge), P_applied = lambdaP sigma0 t H (|Rx_west| should equal it), F_spring (spring force),
      hH_num (largest y/H of the profile points that are not taut, 0 if none), hH_ana (Eq. 40).
      The response columns are NaN for attempts that did not converge.
    - profile_r<r>.csv at the end of phase 1 (r = 0.000) and at every checkpoint of phase 2:
      101 points on x = 0 (y = eta H). Columns: y, yH = y/H, Sxx, Syy, Sxy (Cartesian PK2
      membrane stress of the run material at z = 0), S1, S2 (principal stresses, S1 >= S2), state
      (1 taut, 0 wrinkled, -1 slack, judged by a Split oracle with Eta = 0 around the base law),
      Sxx_s0 = Sxx/sigma0, Sxx_ana_s0 and hH_ana (analytic solution).

    Analytic solution (Stein & Hedgepeth 1961; Zhang & Kiendl 2025, Eqs. 40-41), with
    m = M/(P H) = r/2 and h the height of the wrinkled zone measured from y = 0:
    - h/H = 0 for m < 1/6 (r < 1/3), h/H = 3m - 1/2 = 1.5 r - 0.5 for 1/6 <= m < 1/2 (Eq. 40);
    - m < 1/6, no wrinkles: S_xx/sigma0 = 1 + 3 r (2y/H - 1);
    - m >= 1/6 (Eq. 41): S_xx/sigma0 = 0 for y/H <= h/H and 2 (y/H - h/H)/(1 - h/H)^2 above.
    Both branches give 2 y/H at r = 1/3. For r = 0.825: h/H = 0.7375, S_xx/sigma0(y = H) = 7.62.
    The analytic solution assumes a stress criterion and is matched exactly only for nu = 0 (Zhang
    & Kiendl 2024, Sec. 5.1); \a --nu 0 reproduces that setting.

    References:
    - Stein, M., Hedgepeth, J.M. (1961), Analysis of partly wrinkled membranes, NASA TN D-813.
    - Zhang, Y., Kiendl, J. (2025), J. Mech. Phys. Solids 205, 106331, Sec. 5.1.
    - Zhang, Y., Kiendl, J. (2024), Comput. Methods Appl. Mech. Engrg. 432, 117386, Sec. 5.1.
    - Nakashino, K., Natori, M.C. (2005), Efficient modification scheme of stiffness matrices to
      account for membrane wrinkling in structural analysis of inflatable structures.

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

#include <gsStructuralAnalysis/src/gsStaticSolvers/gsStaticNewton.h>

#include <gsStructuralAnalysis/src/gsStructuralAnalysisTools/gsStructuralAnalysisUtils.h>

using namespace gismo;

#ifdef gsKLShell_ENABLED
namespace
{

typedef gsMaterialMatrixBase<real_t> MaterialBase;

/// Evaluates a material output of \a mat on \a def at the parametric points \a nodes of patch 0 and
/// the through-thickness coordinates \a Z (one column, normalised, +-0.5 = outer fibres, 0 = mid-surface).
/// Column j*nodes.cols()+k of the result belongs to fibre j and point k.
template <MaterialOutput out>
gsMatrix<> evalMaterial(MaterialBase * mat, const gsMultiPatch<> & def, const gsMatrix<> & nodes, const gsMatrix<> & Z)
{
    gsMaterialMatrixEval<real_t,out> ev(mat,&def,Z);
    gsMatrix<> result;
    ev.piece(0).eval_into(nodes,result);
    return result;
}

/// Builds the wrinkling-state oracle around \a base (not owned): a Split decorator with Eta = 0, no exact
/// solve and S2crit = 1e300, so that the TensionField output is the criterion state; the oracle
/// evaluates the wrinkling criterion \a tensionField
template <short_t d, bool linear>
gsMaterialMatrixBase<real_t>::uPtr makeOracle(MaterialBase * base, real_t slackMultiplier, index_t tensionField)
{
    gsMaterialMatrixBase<real_t>::uPtr oracle =
        memory::make_unique(new gsMaterialMatrixTFTSplit<d,real_t,linear>(base));
    gsOptionList o = oracle->options();
    o.setReal("Eta",0);
    o.setSwitch("Exact",false);
    o.setReal("S2crit",1e300);
    o.setReal("SlackMultiplier",slackMultiplier);
    o.setInt("TensionField",tensionField);
    oracle->setOptions(o);
    return oracle;
}

/// Sets the options of a Split material, including its own wrinkling criterion \a tensionField
void setSplitOptions(MaterialBase * split, real_t eta, bool exact, bool midsurface, real_t s2crit, real_t slackMultiplier, index_t tensionField)
{
    gsOptionList o = split->options();
    o.setReal("Eta",eta);
    o.setSwitch("Exact",exact);
    o.setSwitch("Midsurface",midsurface);
    o.setReal("S2crit",s2crit);
    o.setReal("SlackMultiplier",slackMultiplier);
    o.setInt("TensionField",tensionField);
    split->setOptions(o);
}

/// Creates \a path and its missing parents
bool makeDirs(const std::string & path)
{
    bool ok = true;
    for (size_t pos = path.find('/',1); pos != std::string::npos; pos = path.find('/',pos+1))
        ok = gsFileManager::mkdir(path.substr(0,pos)) && ok;
    return gsFileManager::mkdir(path) && ok;
}

/// A phase of the load schedule: the load coordinate runs from \a from to \a to in \a nSteps nominal
/// steps and hits every value of \a checkpoints (and \a to) exactly
struct Phase
{
    real_t from, to;
    index_t nSteps;
    std::vector<real_t> checkpoints;
};

/// Height of the wrinkled zone h/H of the analytic solution (Eq. 40) at the moment ratio r = 2M/(PH)
real_t analyticHeight(real_t r)
{
    return (r < 1./3.) ? 0.0 : 1.5*r - 0.5;
}

/// Axial stress S_xx/sigma0 of the analytic solution at y/H = \a yH (Eq. 41, and the linear elastic
/// branch for r < 1/3)
real_t analyticSxx(real_t r, real_t yH)
{
    if (r < 1./3.)
        return 1 + 3*r*(2*yH - 1);
    const real_t h = analyticHeight(r);
    return (yH <= h) ? 0.0 : 2*(yH - h)/((1 - h)*(1 - h));
}

/// Stress and wrinkling state along the line x = 0, sampled at \a y/H = i/(n-1)
struct Profile
{
    std::vector<real_t> y, Sxx, Syy, Sxy, S1, S2, state;

    /// Largest y/H of the points that are not taut, 0 if there is none
    real_t wrinkleHeight(real_t H) const
    {
        real_t h = 0;
        for (size_t i = 0; i!=y.size(); ++i)
            if (state[i] != 1)
                h = std::max(h,y[i]/H);
        return h;
    }
};

/**
 * @brief Evaluates the profile of the run material \a run on \a def along x = 0 of patch 0.
 *
 * The material returns the PK2 stress in contravariant components of the reference covariant basis
 * a_i = dX/dxi_i. With the reference Jacobian J = [a_1 a_2], the Cartesian stress is
 * S = J [S^11 S^12; S^12 S^22] J^T. The principal stresses are invariants and are taken as returned.
 * Each point is evaluated in a separate call, because the multi-point evaluation of the wrapped
 * materials reuses the metric of the last point.
 */
Profile evalProfile(MaterialBase * run, MaterialBase * oracle, const gsMultiPatch<> & mp, const gsMultiPatch<> & def,
                    real_t H, index_t n)
{
    const gsMatrix<> Z0 = gsMatrix<>::Zero(1,1);
    Profile p;
    p.y.resize(n); p.Sxx.resize(n); p.Syy.resize(n); p.Sxy.resize(n);
    p.S1.resize(n); p.S2.resize(n); p.state.resize(n);
    gsMatrix<> pt(2,1);
    for (index_t i = 0; i!=n; ++i)
    {
        const real_t eta = static_cast<real_t>(i)/(n-1);
        pt<<0.0,eta;
        const gsMatrix<> S  = evalMaterial<MaterialOutput::Stress>(run,def,pt,Z0);
        const gsMatrix<> PS = evalMaterial<MaterialOutput::PStress>(run,def,pt,Z0);
        const gsMatrix<> TF = evalMaterial<MaterialOutput::TensionField>(oracle,def,pt,Z0);
        gsMatrix<> Scon(2,2);
        Scon<<S(0,0),S(2,0),S(2,0),S(1,0);
        const gsMatrix<> J = mp.patch(0).jacobian(pt);
        const gsMatrix<> Scart = J*Scon*J.transpose();
        p.y[i]   = eta*H;
        p.Sxx[i] = Scart(0,0);
        p.Syy[i] = Scart(1,1);
        p.Sxy[i] = Scart(0,1);
        p.S1[i]  = std::max(PS(0,0),PS(1,0));
        p.S2[i]  = std::min(PS(0,0),PS(1,0));
        p.state[i] = TF(0,0);
    }
    return p;
}

/// Writes \a p and the analytic solution at the moment ratio \a r as CSV. Returns the maximum of
/// |Sxx_s0 - Sxx_ana_s0|.
real_t writeProfile(const std::string & fileName, const Profile & p, real_t r, real_t sigma0, real_t H)
{
    std::ofstream file(fileName);
    file<<std::setprecision(12)<<"y,yH,Sxx,Syy,Sxy,S1,S2,state,Sxx_s0,Sxx_ana_s0,hH_ana"<<std::endl;
    real_t maxDev = 0;
    for (size_t i = 0; i!=p.y.size(); ++i)
    {
        const real_t yH = p.y[i]/H;
        const real_t Sxx_s0 = p.Sxx[i]/sigma0;
        const real_t Sxx_ana = analyticSxx(r,yH);
        maxDev = std::max(maxDev,std::abs(Sxx_s0 - Sxx_ana));
        file<<p.y[i]<<","<<yH<<","<<p.Sxx[i]<<","<<p.Syy[i]<<","<<p.Sxy[i]<<","<<p.S1[i]<<","<<p.S2[i]
            <<","<<p.state[i]<<","<<Sxx_s0<<","<<Sxx_ana<<","<<analyticHeight(r)<<"\n";
    }
    file.flush();
    return maxDev;
}

/// Numeric literal with full double precision, for use in function expressions
std::string literal(real_t v)
{
    std::ostringstream ss;
    ss<<std::setprecision(17)<<v;
    return ss.str();
}

}// namespace

int main (int argc, char** argv)
{
    // Geometry and material parameters of Zhang & Kiendl (2025), Sec. 5.1 (dimensionless, consistent units)
    const real_t Lhalf     = 1.0;    // half length L/2 of the modelled right half
    const real_t H         = 1.0;
    const real_t thickness = 0.01;
    const real_t E_modulus = 100.0;
    const real_t sigma0    = 5e-4;
    const index_t degree   = 2;
    const index_t nxMain   = 10;     // elements of the wrinkling patch in x
    const index_t ny       = 5;
    const index_t nPoints  = 101;    // sample points of a profile
    // Stiffness of the spring that fixes u_y at (0,H/2): E t is the membrane stiffness scale per DOF
    const real_t kSpring   = 1e3 * E_modulus * thickness;

    // Input options
    std::string modelName("split");
    real_t eta         = 0;
    bool etaTangent    = false;
    index_t tensionField = 2;
    real_t slackmult   = 1e-6;
    real_t PoissonRatio = 0.3;
    index_t presteps   = 1;
    index_t steps      = 33;
    real_t mmax        = 0.825;
    // The loads are ~sigma0*t = 5e-6, so the reference residual of a step is ~1e-7 and a relative
    // tolerance of 1e-8 asks for a residual below the round-off floor of the assembled residual
    // (~eps times the reaction scale); the criterion would then be met only by chance.
    real_t tolNR       = 1e-6;
    index_t maxitNR    = 50;
    std::string linSolver("LU");
    bool lineSearch = false;
    index_t maxBisect  = 8;
    std::string outDir;
    index_t verbose    = 0;

    gsCmdLine cmd("In-plane pure bending of a prestressed membrane (Stein & Hedgepeth), Zhang & Kiendl (2025) Sec. 5.1.");

    cmd.addString("","model","Wrinkling model: none, tft or split (default split)",modelName);
    cmd.addReal("","eta","TFT and Split: stress fraction retained in the wrinkled/slack state (default 0)",eta);
    cmd.addSwitch("etatangent","TFT: apply --eta to the tangent only (stress stays the pure tension field)",etaTangent);
    cmd.addInt("","tensionfield","Wrinkling criterion: 0 strain-based, 1 stress-based, 2 mixed (default 2)",tensionField);
    cmd.addReal("","slackmult","Multiplier of the stiffness in the slack state, TFT and Split (default 1e-6)",slackmult);
    cmd.addReal("","nu","Poisson ratio (default 0.3; 0 matches the analytic solution exactly)",PoissonRatio);
    cmd.addInt("","presteps","Number of load steps of phase 1, prestress and axial load (default 1)",presteps);
    cmd.addInt("N","steps","Nominal number of load steps of phase 2, from r = 0 to --mmax (default 33)",steps);
    cmd.addReal("","mmax","Final moment ratio r = 2M/(PH), below 1 (default 0.825)",mmax);
    cmd.addReal("","tolNR","Relative tolerance of the Newton solver, residual and update (default 1e-6)",tolNR);
    cmd.addInt("","maxitNR","Maximum number of Newton iterations (default 50)",maxitNR);
    cmd.addSwitch("linesearch", "Newton: backtracking line search on the residual norm (gsStaticNewton option LineSearch)", lineSearch);
    cmd.addString("","solver","Linear solver of the Newton iterations: LU, CGDiagonal, SimplicialLDLT, ... (default LU)",linSolver);
    cmd.addInt("","maxbisect","Maximum number of successive halvings of a load step before the run stops (default 8)",maxBisect);
    cmd.addString("","out","Output directory (default: ./SteinHedgepeth_<model>)",outDir);
    cmd.addInt("v","verbose","0: no; 1: Newton iteration output",verbose);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    if (!(eta >= 0 && eta <= 1)) { gsWarn << "--eta must be in [0,1], got " << eta << "\n"; return EXIT_FAILURE; }

    int model = 0;
    if      (modelName=="none")  model = 0;
    else if (modelName=="tft")   model = 1;
    else if (modelName=="split") model = 2;
    else
        GISMO_ERROR("Unknown model \""<<modelName<<"\"; use none, tft or split.");
    if (model==0 && eta!=0)
        gsWarn<<"--eta is ignored without a wrinkling model.\n";
    if (model!=1 && etaTangent)
        gsWarn<<"--etatangent only acts together with --model tft and is ignored.\n";
    if (tensionField<0 || tensionField>2)
        GISMO_ERROR("--tensionfield must be 0, 1 or 2, got "<<tensionField);
    if (presteps<1 || steps<1 || maxBisect<0)
        GISMO_ERROR("--presteps and --steps must be at least 1 and --maxbisect must not be negative.");
    if (mmax<=0)
        GISMO_ERROR("--mmax must be positive, got "<<mmax);
    if (mmax>=1)
        GISMO_ERROR("--mmax must be below 1: at r = 1 (M/(PH) = 1/2) the analytic wrinkled height reaches H, got "<<mmax);
    if (outDir.empty())
        outDir = "./SteinHedgepeth_" + modelName;
    if (!makeDirs(outDir))
    {
        gsWarn<<"Could not create the output directory "<<outDir<<"\n";
        return 1;
    }

    gsStopwatch totalTime, stageTime;
    stageTime.restart();

    // Two patches with the Greville points of the spline basis as control points (bilinear maps are
    // reproduced exactly): the wrinkling patch [0,10/11] x [0,1] and the taut cap [10/11,1] x [0,1]
    const real_t xInterface = Lhalf*nxMain/(nxMain+1.0);
    gsKnotVector<> kx0(0,1,nxMain-1,degree+1), kx1(0,1,0,degree+1), ky(0,1,ny-1,degree+1);
    gsTensorBSplineBasis<2,real_t> tb0(kx0,ky), tb1(kx1,ky);
    const gsMatrix<> anchors0 = tb0.anchors();
    const gsMatrix<> anchors1 = tb1.anchors();
    gsMatrix<> coefs0(tb0.size(),2), coefs1(tb1.size(),2);
    for (index_t i = 0; i!=coefs0.rows(); ++i)
    {
        coefs0(i,0) = xInterface*anchors0(0,i);
        coefs0(i,1) = H*anchors0(1,i);
    }
    for (index_t i = 0; i!=coefs1.rows(); ++i)
    {
        coefs1(i,0) = xInterface + (Lhalf-xInterface)*anchors1(0,i);
        coefs1(i,1) = H*anchors1(1,i);
    }
    gsMultiPatch<> mp, mp_def;
    mp.addPatch(gsTensorBSpline<2,real_t>(tb0,coefs0));
    mp.addPatch(gsTensorBSpline<2,real_t>(tb1,coefs1));
    mp.computeTopology();
    GISMO_ENSURE(mp.topology().interfaces().size()==1,"The two patches do not form a single conforming interface");
    gsMultiBasis<> dbasis(mp,true);
    mp_def = mp;

    // Control point on the symmetry edge at y = H/2, pinned in u_y
    const index_t kPin = tb0.index(0,3);
    GISMO_ENSURE(math::abs(anchors0(1,kPin)-0.5)<1e-14 && anchors0(0,kPin)==0,"The pinned control point is not (0,H/2)");

    // Loads: force per unit undeformed length, lambda_P = min(t,1) scales prestress and axial load,
    // r = max(t-1,0) is the moment ratio (Zhang & Kiendl (2025) Sec. 5.1, Fig. 5(b)):
    //   top    y = H: ( 0,  s0t lambda_P)
    //   bottom y = 0: ( 0, -s0t lambda_P)
    //   right  x = L/2: ( s0t lambda_P + 3 s0t r (2y/H - 1), 0 ),   s0t = sigma0 t,   H = 1
    const std::string s0t = literal(sigma0*thickness);
    gsFunctionExpr<> fTop("0",s0t+"*min(t,1)",2);
    gsFunctionExpr<> fBot("0","-"+s0t+"*min(t,1)",2);
    gsFunctionExpr<> fRight(s0t+"*min(t,1) + 3*"+s0t+"*max(t-1,0)*(2*y-1)","0",2);

    gsConstantFunction<> zeroFun(0.0,2);
    gsBoundaryConditions<> BCs;
    BCs.addCondition(0,boundary::west,condition_type::dirichlet,&zeroFun,0,false,0);
    BCs.addCondition(0,boundary::north,condition_type::neumann,&fTop,0,false,-1);
    BCs.addCondition(1,boundary::north,condition_type::neumann,&fTop,0,false,-1);
    BCs.addCondition(0,boundary::south,condition_type::neumann,&fBot,0,false,-1);
    BCs.addCondition(1,boundary::south,condition_type::neumann,&fBot,0,false,-1);
    BCs.addCondition(1,boundary::east,condition_type::neumann,&fRight,0,false,-1);
    BCs.setGeoMap(mp);

    // Materials: linear SvK on both patches; the wrinkling model wraps the base law of patch 0 only
    gsConstantFunction<> t(thickness,2), E(E_modulus,2), nu(PoissonRatio,2), rho(1.0,2);
    std::vector<gsFunctionSet<>*> parameters{&E,&nu};
    gsOptionList options;
    options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
    options.addSwitch("Compressibility","Compressibility: (false): Imcompressible | (true): Compressible",false);
    options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);

    gsMaterialMatrixBase<real_t>::uPtr base = getMaterialMatrix<2,real_t>(mp,t,parameters,rho,options);
    gsMaterialMatrixBase<real_t>::uPtr cap  = getMaterialMatrix<2,real_t>(mp,t,parameters,rho,options);
    // The TFT decorator clones the base law and takes its criterion from the clone, so the
    // option has to be set before anything wraps the base.
    base->options().setInt("TensionField",tensionField);
    gsMaterialMatrixBase<real_t>::uPtr oracle = makeOracle<2,true>(base.get(),slackmult,tensionField);
    gsMaterialMatrixBase<real_t>::uPtr wrapped;
    MaterialBase * run = base.get();
    if (model==1)
    {
        wrapped = memory::make_unique(new gsMaterialMatrixTFT<2,real_t,true>(*base));
        wrapped->options().setReal("SlackMultiplier",slackmult);
        wrapped->options().setReal("Eta",eta);
        wrapped->options().setSwitch("EtaTangentOnly",etaTangent);
        wrapped->options().setInt("TensionField",tensionField);
        run = wrapped.get();
    }
    else if (model==2)
    {
        wrapped = memory::make_unique(new gsMaterialMatrixTFTSplit<2,real_t,true>(*base));
        setSplitOptions(wrapped.get(),eta,false,false,0,slackmult,tensionField);
        run = wrapped.get();
    }

    // Pointer overload of add: the container does not own the materials
    gsMaterialMatrixContainer<real_t> materials;
    materials.add(run);
    materials.add(cap.get());

    gsConstantFunction<> force(gsVector<>::Zero(2),2);
    std::unique_ptr<gsThinShellAssembler<2,real_t,false> > assembler(
        new gsThinShellAssembler<2,real_t,false>(mp,dbasis,BCs,force,materials));
    // Strong C^0 coupling of the interface; the default weak C^1 coupling does not act in 2D
    assembler->options().setInt("Continuity",0);
    GISMO_ENSURE(assembler->material(0)==run && assembler->material(1)==cap.get(),"The assembler does not use the expected material objects");

    assembler->assemble();
    GISMO_ENSURE(assembler->numDofs()==189,"Unexpected number of degrees of freedom "<<assembler->numDofs()<<" (the interface is not glued?)");

    // The mapper is refreshed by the first assemble (re-initialisation of the space)
    const gsDofMapper mapper = assembler->getMapper();
    GISMO_ENSURE(mapper.is_free(kPin,0,1),"The pinned u_y degree of freedom is not free");
    const index_t iPin = mapper.index(kPin,0,1);

    gsSparseMatrix<> K = assembler->matrix();
    gsVector<> F = assembler->rhs();
    K.coeffRef(iPin,iPin) += kSpring;

    // Newton solves J du = R with R = F_ext - F_int; the spring adds kSpring u_pin to F_int
    gsStructuralAnalysisOps<real_t>::Jacobian_t Jacobian = [&assembler,&mp_def,iPin,kSpring](gsVector<real_t> const &x, gsSparseMatrix<real_t> & m)
    {
        assembler->constructSolution(x,mp_def);
        ThinShellAssemblerStatus status = assembler->assembleMatrix(mp_def);
        m = assembler->matrix();
        m.coeffRef(iPin,iPin) += kSpring;
        return status == ThinShellAssemblerStatus::Success;
    };
    gsStructuralAnalysisOps<real_t>::Residual_t Residual = [&assembler,&mp_def,iPin,kSpring](gsVector<real_t> const &x, gsVector<real_t> & result)
    {
        assembler->constructSolution(x,mp_def);
        ThinShellAssemblerStatus status = assembler->assembleVector(mp_def);
        result = assembler->rhs();
        result[iPin] -= kSpring*x[iPin];
        return status == ThinShellAssemblerStatus::Success;
    };

    gsStaticNewton<real_t> NWT(K,F,Jacobian,Residual);
    gsOptionList NWTOptions = NWT.options();
    NWTOptions.setInt("maxIt",maxitNR);
    NWTOptions.setReal("tol",tolNR);
    NWTOptions.setString("Solver",linSolver);
    NWTOptions.setSwitch("LineSearch",lineSearch);
    NWTOptions.setInt("verbose",verbose);
    NWT.setOptions(NWTOptions);

    gsInfo<<"Model "<<modelName<<", nu = "<<PoissonRatio<<", setup time: "<<stageTime.stop()
          <<" s; degrees of freedom: "<<assembler->numDofs()<<"; pinned dof "<<iPin<<"; kSpring = "<<kSpring<<"\n";

    // Load schedule. The load coordinate D of a phase is
    //   phase 1: lambda_P (prestress and axial load), global t = D
    //   phase 2: r = 2M/(PH), global t = 1 + D
    std::vector<Phase> phases;
    {
        Phase p1 = {0.0,1.0,presteps,std::vector<real_t>()};
        Phase p2 = {0.0,mmax,steps,std::vector<real_t>()};
        const real_t cps[] = {0.4,0.6,0.75,0.825};
        for (size_t i = 0; i!=4; ++i)
            if (cps[i] < mmax-1e-12)
                p2.checkpoints.push_back(cps[i]);
        phases.push_back(p1);
        phases.push_back(p2);
    }
    for (size_t ph = 0; ph!=phases.size(); ++ph)
        phases[ph].checkpoints.push_back(phases[ph].to);

    auto applyLoad = [&](size_t ph, real_t D)
    {
        const real_t tg = (ph==0) ? D : 1.0 + D;
        fTop.set_t(tg);
        fBot.set_t(tg);
        fRight.set_t(tg);
        assembler->updateBCs(BCs);
    };

    const real_t nan = std::numeric_limits<real_t>::quiet_NaN();

    std::ofstream stepsFile(outDir + "/steps.csv");
    if (!stepsFile)
    {
        gsWarn<<"Could not open "<<outDir<<"/steps.csv for writing\n";
        return 1;
    }
    stepsFile<<std::setprecision(12)
             <<"step,attempt,phase,D,lambdaP,r,itNR,converged,halvings,Rx_west,Ry_west,P_applied,F_spring,hH_num,hH_ana"
             <<std::endl;

    gsVector<> solVector = gsVector<>::Zero(assembler->numDofs());
    index_t step = 0, attempt = 0, bisected = 0;
    bool failed = false;
    // Last converged value of the load parameter D of the phase in which the run stands, and that phase
    real_t lastConverged = phases.front().from;
    size_t stopPhase = 0;
    gsStopwatch stopwatch;
    real_t solveTimeTotal = 0;

    for (size_t ph = 0; ph!=phases.size() && !failed; ++ph)
    {
        const Phase & P = phases[ph];
        lastConverged = P.from;
        stopPhase = ph;
        const real_t dL0 = (P.to-P.from)/P.nSteps;
        const real_t tolD = 1e-12*std::max<real_t>(1.0,std::abs(P.to));
        size_t cpIdx = 0;
        real_t Dold = P.from;
        real_t dL = std::min(dL0,P.checkpoints[cpIdx]-Dold);
        index_t halvings = 0;

        while (!failed)
        {
            // Snap to the checkpoint when the step reaches it
            const bool hitsCheckpoint = (dL >= P.checkpoints[cpIdx]-Dold-tolD);
            const real_t D = hitsCheckpoint ? P.checkpoints[cpIdx] : Dold + dL;
            const real_t lambdaP = (ph==0) ? D : 1.0;
            const real_t r = (ph==0) ? 0.0 : D;

            gsInfo<<"Load step "<<step<<" (attempt "<<attempt<<", phase "<<ph+1<<"): D = "<<D<<"; dD = "<<D-Dold<<"\n";
            applyLoad(ph,D);

            NWT.initialize();
            NWT.setDisplacement(solVector);
            stopwatch.restart();
            NWT.solve();
            const real_t solveTime = stopwatch.stop();
            solveTimeTotal += solveTime;
            const bool converged = NWT.converged();
            gsInfo<<"Solve time: "<<solveTime<<" s; Newton iterations: "<<NWT.iterations()<<"\n";

            real_t Rx = nan, Ry = nan, Fspring = nan, hNum = nan;
            const real_t hAna = analyticHeight(r);
            Profile prof;
            if (converged)
            {
                solVector = NWT.solution();
                mp_def = assembler->constructSolution(solVector);
                const gsMatrix<> bf = assembler->boundaryForce(mp_def,patchSide(0,boundary::west));
                Rx = bf(0,0);
                Ry = bf(1,0);
                Fspring = kSpring*solVector[iPin];
                prof = evalProfile(run,oracle.get(),mp,mp_def,H,nPoints);
                hNum = prof.wrinkleHeight(H);
            }

            stepsFile<<step<<","<<attempt<<","<<ph+1<<","<<D<<","<<lambdaP<<","<<r<<","<<NWT.iterations()<<","
                     <<(converged ? 1 : 0)<<","<<halvings<<","<<Rx<<","<<Ry<<","<<lambdaP*sigma0*thickness*H<<","
                     <<Fspring<<","<<hNum<<","<<hAna<<std::endl;
            ++attempt;

            if (!converged)
            {
                gsWarn<<"Load step "<<step<<" did not converge (status "<<static_cast<int>(NWT.status())<<")\n";
                if (halvings>=maxBisect)
                {
                    gsWarn<<"Simulation terminated: the load step was halved "<<maxBisect<<" times without convergence\n";
                    failed = true;
                    break;
                }
                dL = (D-Dold)/2;
                halvings++;
                bisected++;
                continue;
            }

            // Profile along x = 0 at the end of phase 1 (r = 0) and at every checkpoint of phase 2
            if (hitsCheckpoint)
            {
                std::ostringstream name;
                name<<outDir<<"/profile_r"<<std::fixed<<std::setprecision(3)<<r<<".csv";
                const real_t maxDev = writeProfile(name.str(),prof,r,sigma0,H);
                gsInfo<<"Checkpoint r = "<<r<<": hH_num = "<<hNum<<", hH_ana = "<<hAna
                      <<", max|Sxx_s0 - Sxx_ana_s0| = "<<maxDev<<"\n";
            }

            Dold = D;
            lastConverged = D;
            halvings = 0;
            step++;
            if (hitsCheckpoint)
            {
                cpIdx++;
                if (cpIdx==P.checkpoints.size())
                    break;
            }
            dL = std::min(dL0,P.checkpoints[cpIdx]-Dold);
        }
    }

    gsInfo<<"Total time: "<<totalTime.stop()<<" s (Newton solves: "<<solveTimeTotal<<" s); accepted steps: "
          <<step<<"; attempts: "<<attempt<<"; halvings: "<<bisected<<"\n";
    const bool reached = !failed;
    if (!reached)
        gsInfo<<"FINAL LOAD NOT REACHED: last converged "<<lastConverged<<" of "<<phases[stopPhase].to<<"\n";
    return reached ? EXIT_SUCCESS : EXIT_FAILURE;
}
#else//gsKLShell_ENABLED
int main(int argc, char *argv[])
{
    GISMO_UNUSED(argc); GISMO_UNUSED(argv);
    gsWarn<<"G+Smo is not compiled with the gsKLShell module.";
    return EXIT_FAILURE;
}
#endif
