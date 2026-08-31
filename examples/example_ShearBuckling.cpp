/** @file example_ShearBuckling.cpp

    @brief Linear (prescribed-displacement) buckling reference for the shear
    sheet, independent of the arc-length method (ALM).

    This driver poses and solves the *linear* buckling eigenproblem on the
    same sheet, material, thickness and boundary conditions as
    example_ShearExploration.cpp, and reports a critical load factor and a
    critical-mode half-wave count that are directly comparable to the ALM's
    own bifurcation load and traced-branch mode. It is a separate numerical
    method (a generalized eigenproblem, never the nonlinear arc-length path)
    so that the ALM's mode identification can be checked against an
    independent measurement instead of only against itself.

    -------------------------------------------------------------------------
    Formulation
    -------------------------------------------------------------------------
    The north edge is driven by a prescribed x-displacement lambda*d over a
    sheet of width bDim = 1, so lambda is numerically the engineering shear
    strain -- the same lambda the ALM reports. Under linear kinematics and a
    linear material, the response to lambda*d is exactly lambda*u_L(1): the
    pre-buckling membrane stress is linear in lambda, and so is the
    initial-stress (geometric) stiffness K_G. The buckling mode is a
    perturbation of the free (non-Dirichlet) degrees of freedom, so the
    classical bifurcation condition

        [ K_L + lambda K_G ] v = 0

    holds on exactly the reduced (Dirichlet-eliminated) system the shell
    assembler already produces -- only the way u_L(1) is obtained differs
    from the constant-force case: it is the solution of the LINEAR system at
    unit prescribed shear, K_L u_L(1) = assembler->rhs()|_{lambda=1,U_F=0},
    the standard Dirichlet-elimination right-hand side gsExprAssembler
    produces for an inhomogeneous BC with zero external force.

    Beware: this is NOT the "Force = -rhs()" vector example_ShearExploration.cpp
    builds for the ALM's own predictor -- gsALMBase applies that vector
    inside its own Newton-step bookkeeping, which carries an implicit sign
    flip relative to a standalone K^-1*rhs solve of the free-DoF equilibrium.
    Using the ALM's convention directly here would flip the sign of u_L(1)
    (and hence of the whole pre-buckling state); the residual-discrimination
    check below ("pre-buckling sign check") is what a caller must use to
    confirm which convention a given assembler/BC setup actually needs,
    rather than assume either one.

    K_G is defined the same way gsBucklingSolver defines its pencil's first
    matrix (K_G := K_NL(w) - K_L, gsBucklingSolver.h ~:170-172), evaluated as
    the directional derivative of the tangent stiffness at lambda=0:

        K_G(eps) = [ K_NL(eps*u_L(1)) - K_L ] / eps

    K_NL(w) - K_L also carries the initial-displacement term, which is
    quadratic in w, so K_G(eps) is contaminated at O(eps); an eps-sweep and a
    central-difference control (which cancels the even part of that
    contamination to O(eps^2)) both measure the size of that error rather
    than assume it away.

    -------------------------------------------------------------------------
    Why the eigenproblem is posed as K_G v = nu K_L v, not K_L v = mu K_G v
    -------------------------------------------------------------------------
    gsBucklingSolver itself poses K_G v = nu K_L v via
    GeneralizedSelfAdjointEigenSolver::compute(A,B), which Cholesky-factors
    B: only the symmetric positive definite K_L is ever factored, and
    gsBucklingSolver certifies that positive-definiteness before solving,
    reporting a non-Success status instead of a silently wrong answer when it
    does not hold. This driver poses the identical formulation independently
    of the library under test -- being a separate code path from
    gsBucklingSolver, not a call into it, is the entire evidential value of
    comparing the two -- by solving

        K_G v = nu K_L v ,   K_L symmetric positive definite,

    so nu < 0 <=> lambda = -1/nu > 0, and modes with negligible pre-stress
    coupling land at nu ~ 0, i.e. |lambda| -> infinity, harmlessly at the far
    end of the spectrum rather than as spurious near-zero eigenvalues of the
    K_L-first formulation.

    Complexity: the dense generalized eigensolve is O(n^3) in the number of
    free degrees of freedom (n ~ 300-1000 at the meshes below); do not run
    this driver past -r 4, where n grows into the low thousands.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gismo.h>

#ifdef gsKLShell_ENABLED
#include <gsKLShell/src/gsThinShellAssembler.h>
#include <gsKLShell/src/getMaterialMatrix.h>
#endif

using namespace gismo;

template <class T>
gsMultiPatch<T> Rectangle(T L, T B);

#ifdef gsKLShell_ENABLED

/// Infinity norm (max absolute row sum) of a sparse matrix, via its dense
/// form -- these systems are a few hundred to a few thousand DoFs, so
/// densifying once for a norm is negligible next to the O(n^3) eigensolve.
real_t matInfNorm(const gsSparseMatrix<real_t> & M)
{
    return M.toDense().cwiseAbs().rowwise().sum().maxCoeff();
}

/// Sign-change counter for a sampled displacement profile, referenced to
/// PHYSICAL ZERO rather than to the profile's own mean. A mean-referenced
/// counter is wrong for a one-sided lobe: sin(pi*s) is a single non-negative
/// hump, but subtracting its own mean manufactures two spurious crossings
/// near the domain ends where the samples dip below that mean. Referencing
/// zero instead reproduces the textbook count directly: sin(k*pi*s) has
/// exactly k-1 interior sign changes for every k, matching the half-wave
/// convention below. Samples with |u| at or below relTol of the profile's
/// own peak magnitude are treated as noise and dropped before counting, so a
/// near-zero station does not manufacture spurious crossings either.
index_t countSignChanges(const gsVector<real_t> & u, real_t relTol)
{
    const real_t scaleu = u.cwiseAbs().maxCoeff();
    if (scaleu == 0.0) return 0;

    std::vector<int> signs;
    signs.reserve(u.size());
    for (index_t i = 0; i != u.size(); ++i)
    {
        if (std::abs(u[i]) > relTol * scaleu)
            signs.push_back(u[i] > 0 ? 1 : -1);
    }
    if (signs.size() < 2) return 0;

    index_t n = 0;
    for (size_t i = 1; i != signs.size(); ++i)
        if (signs[i] != signs[i-1]) ++n;
    return n;
}

// A half-wave is one lobe between consecutive zeros (or between a zero and
// an end of the sampling line): n interior sign changes => n+1 half-waves;
// a profile that never changes sign carries 1 half-wave if it is genuine
// signal, 0 if it is flat. countSignChanges alone cannot distinguish "one
// smooth lobe" from "no signal at all" (both give 0 sign changes for the
// same reason a constant array does), so flatness is judged separately, on
// an ABSOLUTE scale: a station whose own peak magnitude is negligible next
// to the mode's global out-of-plane amplitude (refScale) carries no
// resolvable signal, independent of relTol. relTol governs only which
// *interior* wiggles of an already-significant profile count as a genuine
// sign change; it cannot by itself flag "no signal", because
// countSignChanges renormalizes by the profile's own peak -- literal
// floating-point noise, once it is nonzero at all, looks the same as a
// genuine oscillation to a purely self-normalized filter (verified in the
// self-test below: i.i.d. noise at any absolute scale produces ~n/2 sign
// changes under countSignChanges alone, not 0).
index_t halfWaveCount(const gsVector<real_t> & u, real_t relTol, real_t refScale)
{
    static const real_t flatFrac = 1e-6; // station peak / mode's global amplitude
    const real_t scaleu = u.cwiseAbs().maxCoeff();
    if (refScale > 0 && scaleu < flatFrac * refScale) return 0;
    if (scaleu == 0.0) return 0;
    return countSignChanges(u, relTol) + 1;
}

/// Samples the z-component of \a mp_mode along direction \a dir (0=u/length,
/// 1=v/width) at the fixed parametric station \a other, using \a nSamples
/// points on [0,1]. The undeformed sheet has z==0 everywhere, so this is the
/// mode's z-displacement profile directly. Never build the instrument on the
/// in-plane components -- they carry the Dirichlet lift, not the wrinkle.
gsVector<real_t> zProfile(const gsMultiPatch<real_t> & mp_mode, index_t dir,
                           real_t other, index_t nSamples)
{
    gsMatrix<real_t> uv(2, nSamples);
    gsVector<real_t> line = gsVector<real_t>::LinSpaced(nSamples, 0.0, 1.0);
    if (dir == 0) { uv.row(0) = line.transpose(); uv.row(1).setConstant(other); }
    else          { uv.row(1) = line.transpose(); uv.row(0).setConstant(other); }
    const gsMatrix<real_t> vals = mp_mode.patch(0).eval(uv);
    return vals.row(2).transpose();
}

/// One eigenpair's diagnostics: its physical load factor, the singularity
/// residual that shows it actually solves [K_L + lambda K_G] v = 0, its
/// out-of-plane/in-plane amplitude split, and the sign-change/half-wave
/// count of a single representative profile (v=0.5, along u) used for the
/// per-mode table. The detailed multi-station/multi-direction sweep is done
/// separately, only for the selected critical mode.
struct ModeInfo
{
    index_t idx;
    real_t nu, lambda, residual, vz, vxy;
    index_t signChanges, halfWaves;
    bool outOfPlane, residualOk, nuOk;
};

ModeInfo evalMode(index_t idx, const gsVector<real_t> & nuVec, const gsMatrix<real_t> & V,
                   const gsSparseMatrix<real_t> & K_L, const gsSparseMatrix<real_t> & K_G,
                   real_t nuFloor, gsThinShellAssemblerBase<real_t> * assembler,
                   gsBoundaryConditions<> & BCs, gsConstantFunction<> & displ,
                   const gsMultiPatch<> & mp, real_t halfWaveTol)
{
    ModeInfo info;
    info.idx    = idx;
    info.nu     = nuVec[idx];
    info.lambda = -1.0 / info.nu;

    const gsVector<real_t> v  = V.col(idx);
    const gsVector<real_t> Kv = K_L * v;
    const gsVector<real_t> r  = Kv + info.lambda * (K_G * v);
    info.residual = r.norm() / Kv.norm();

    // Zero Dirichlet lift for the mode shape: the shear BC must NOT be
    // injected into the eigenvector's constructed displacement field, or the
    // north-edge rigid shear would swamp the (much smaller) buckling
    // perturbation and make vz/vxy meaningless.
    gsMultiPatch<> mp_mode = mp;
    displ.setValue(0.0, 3);
    assembler->updateBCs(BCs);
    assembler->constructSolution(v, mp_mode);
    const gsMatrix<real_t> Vd = mp_mode.patch(0).coefs() - mp.patch(0).coefs();
    info.vz  = Vd.col(2).cwiseAbs().maxCoeff();
    info.vxy = std::max(Vd.col(0).cwiseAbs().maxCoeff(), Vd.col(1).cwiseAbs().maxCoeff());

    const gsVector<real_t> prof = zProfile(mp_mode, 0, 0.5, 101);
    info.signChanges = countSignChanges(prof, halfWaveTol);
    info.halfWaves   = halfWaveCount(prof, halfWaveTol, info.vz);

    info.outOfPlane = info.vz > info.vxy;
    info.residualOk = info.residual < 1e-8;
    info.nuOk       = std::abs(info.nu) > nuFloor;
    return info;
}

void printModeRow(const ModeInfo & m, const std::string & tag)
{
    gsInfo << std::setw(5) << std::left << m.idx
           << std::setw(16)<< std::left << m.lambda
           << std::setw(16)<< std::left << m.nu
           << std::setw(14)<< std::left << m.residual
           << std::setw(14)<< std::left << (m.vz / (m.vxy + 1e-30))
           << std::setw(10)<< std::left << m.signChanges
           << std::setw(10)<< std::left << m.halfWaves
           << tag << "\n";
}

/// Walks the spectrum from one end inward (positiveSide: idx=0,1,2,... where
/// nu is most negative, i.e. lambda is smallest positive; else idx=N-1,N-2,...
/// where nu is most positive, i.e. lambda is smallest-magnitude negative),
/// evaluating and printing every candidate examined. Selection stops at the
/// first candidate that passes all three filters (|nu|>nuFloor, singularity
/// residual<1e-8, out-of-plane dominance), but PRINTING continues until at
/// least \a minPrint candidates have been shown (or \a maxWalk is exhausted),
/// so a near-degenerate second or third mode stays visible next to the one
/// actually selected -- a table that stopped at the first pass could hide
/// exactly the near-degenerate case worth seeing. Every rejected candidate is
/// printed with the clause that rejected it; nothing examined is left silent.
index_t walkSelect(index_t nModes, bool positiveSide, const gsVector<real_t> & nuVec,
                    const gsMatrix<real_t> & V, const gsSparseMatrix<real_t> & K_L,
                    const gsSparseMatrix<real_t> & K_G, real_t nuFloor,
                    gsThinShellAssemblerBase<real_t> * assembler, gsBoundaryConditions<> & BCs,
                    gsConstantFunction<> & displ, const gsMultiPatch<> & mp,
                    real_t halfWaveTol, index_t maxWalk, index_t skipUpTo = -1,
                    index_t minPrint = 1)
{
    const index_t start = positiveSide ? 0 : nModes - 1;
    const index_t step  = positiveSide ? 1 : -1;
    index_t idx = start;
    index_t selected = -1;
    index_t printed  = 0;
    for (index_t k = 0; k < maxWalk && idx >= 0 && idx < nModes; ++k, idx += step)
    {
        if (positiveSide && idx <= skipUpTo) continue;
        if (!positiveSide && skipUpTo >= 0 && idx >= skipUpTo) continue;

        ModeInfo info = evalMode(idx, nuVec, V, K_L, K_G, nuFloor, assembler, BCs, displ, mp, halfWaveTol);
        const bool passes = info.nuOk && info.residualOk && info.outOfPlane;
        std::string reason;
        if      (!info.nuOk)       reason = "REJECTED: |nu| below floor (lambda -> infinity)";
        else if (!info.residualOk) reason = "REJECTED: singularity residual >= 1e-8";
        else if (!info.outOfPlane) reason = "REJECTED: in-plane dominated (vz <= vxy)";
        // Only the first passing candidate is the one this walk RETURNS
        // (`selected`); later passing rows are printed for the near-degeneracy
        // table (minPrint) but are not what the caller receives, so tagging
        // them SELECTED as well would claim multiple critical modes per side.
        else if (selected < 0)     reason = "SELECTED";
        else                       reason = "PASSED (not returned; walk continued for minPrint)";
        printModeRow(info, "  " + reason);
        ++printed;
        if (selected < 0 && passes)
            selected = idx;
        if (selected >= 0 && printed >= minPrint)
            break;
    }
    return selected;
}

int main (int argc, char** argv)
{
    // ---------------------- Input options (defaults run in seconds) --------
    index_t numElevate = 2;
    index_t numHref    = 3;
    real_t aDim        = 2.0;
    real_t bDim        = 1.0;
    real_t thickness   = 0.1;
    real_t E_modulus   = 70e3;
    real_t PoissonRatio= 0.3;
    real_t Density     = 2710e9;

    // Comparison to the ALM's own certified critical load; sentinel -1 => the
    // comparison is skipped rather than silently passed.
    real_t lambdaOracle = -1.0;
    // Agreement band on |lambda_lin - lambdaOracle|/lambdaOracle: linear
    // buckling assumes small-strain pre-buckling behaviour, and the ALM's own
    // crossing sits at ~13% shear strain, so a factor-of-2 discrepancy is the
    // expectation, not a defect.
    real_t lambdaBand   = 0.5;
    std::string epsListStr = "1.0,0.1,0.01";
    index_t profileSamples = 201;
    real_t halfWaveTol  = 1e-6;
    // Sentinel -1 => normal filtered selection; else force spectrum index
    // critMode (0-based, ascending nu) as "the" critical mode -- for driving
    // the out-of-plane-dominance check to a reproducible failure only.
    index_t critMode    = -1;
    bool plot           = false;

    gsCmdLine cmd("Independent linear-buckling reference for the shear sheet "
                  "(critical load + mode shape), for comparison against the ALM.");
    cmd.addInt   ("r","hRefine",        "Number of uniform h-refinement steps", numHref);
    cmd.addInt   ("e","degreeElevation","Number of degree elevation steps", numElevate);
    cmd.addReal  ("a","aDim",           "Sheet length (x)", aDim);
    cmd.addReal  ("b","bDim",           "Sheet width (y)", bDim);
    cmd.addReal  ("T","thickness",      "Sheet thickness", thickness);
    cmd.addReal  ("","lambdaOracle",    "ALM certified critical lambda for comparison "
                                         "(negative = skip the comparison)", lambdaOracle);
    cmd.addReal  ("","lambdaBand",      "Max |lambda_lin-lambdaOracle|/lambdaOracle accepted "
                                         "as agreement", lambdaBand);
    cmd.addString("","epsList",         "Comma-separated eps sweep for K_G(eps) (descending "
                                         "recommended; the smallest value is the primary estimate)",
                                         epsListStr);
    cmd.addInt   ("","profileSamples",  "Sampling resolution of the half-wave z-profile", profileSamples);
    cmd.addReal  ("","halfWaveTol",     "Relative tolerance of the sign-change/half-wave counter", halfWaveTol);
    cmd.addInt   ("","critMode",        "Force selection of spectrum index (0-based, ascending nu) "
                                         "as the critical mode, bypassing the selection filter "
                                         "(falsification use only)", critMode);
    cmd.addSwitch("plot", "Write the selected critical mode to ParaView format", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    GISMO_ENSURE(profileSamples >= 2, "--profileSamples must be >= 2, got "<<profileSamples);

    bool allOk = true;
    auto report = [&allOk](bool cond, const std::string & msg)
    {
        gsInfo << (cond ? "[ OK ] " : "[FAIL] ") << msg << "\n";
        if (!cond) allOk = false;
        return cond;
    };
    auto softReport = [](const std::string & msg)
    {
        gsInfo << "[FIND] " << msg << "\n";
    };

    std::vector<real_t> epsList;
    {
        std::stringstream ss(epsListStr);
        std::string tok;
        while (std::getline(ss, tok, ','))
        {
            real_t val;
            try
            {
                size_t consumed = 0;
                val = std::stod(tok, &consumed);
                GISMO_ENSURE(consumed == tok.size(), "trailing characters after the number");
            }
            catch (const std::exception &)
            {
                GISMO_ERROR("--epsList entry '"<<tok<<"' could not be parsed as a real number");
            }
            GISMO_ENSURE(val > 0.0, "--epsList entries must be strictly positive "
                                     "(K_G(eps) divides by eps; sign is applied internally)");
            epsList.push_back(val);
        }
    }
    GISMO_ENSURE(!epsList.empty(), "epsList must contain at least one value");
    const real_t epsPrimary = *std::min_element(epsList.begin(), epsList.end());

    gsInfo << "E = "<<E_modulus<<"; nu = "<<PoissonRatio<<"\n";
    gsInfo << "L = "<<aDim<<"; W = "<<bDim<<"; t = "<<thickness
           <<"; beta = W/t = "<<bDim/thickness<<"\n";
    gsInfo << "hRefine="<<numHref<<" degreeElevation="<<numElevate<<"\n";

    // ---------------------- Self-test (runs on EVERY invocation) -----------
    // sin(k*pi*s), k=1..4, must give the EXACT textbook counts (k-1 interior
    // sign changes, k half-waves each) -- not merely four distinct numbers,
    // which a wrong-but-monotonic counter could also produce. A constant
    // classifier fails every one of these four cases plus the noise case
    // below (it cannot reproduce four different values, and its one fixed
    // output can match at most one of {1,2,3,4}); a classifier that merely
    // orders inputs correctly without matching the convention (as a
    // mean-referenced counter does for a one-sided lobe) is caught by the
    // exact-value check, not just the distinctness check.
    {
        gsInfo << "\n-- half-wave counter self-test --\n";
        const index_t nSelf = 201;
        gsVector<real_t> s = gsVector<real_t>::LinSpaced(nSelf, 0.0, 1.0);
        std::vector<index_t> counts, signChanges;
        for (index_t k = 1; k <= 4; ++k)
        {
            gsVector<real_t> u(nSelf);
            for (index_t i = 0; i != nSelf; ++i) u[i] = std::sin(k*EIGEN_PI*s[i]);
            index_t sc = countSignChanges(u, halfWaveTol);
            index_t hw = halfWaveCount(u, halfWaveTol, 1.0);
            counts.push_back(hw);
            signChanges.push_back(sc);
            gsInfo << "  sin("<<k<<"*pi*s): sign changes="<<sc<<" half-waves="<<hw<<"\n";
        }
        bool exact = true;
        for (index_t k = 1; k <= 4; ++k)
            exact = exact && (signChanges[k-1] == k-1) && (counts[k-1] == k);
        report(exact, "self-test: sin(k*pi*s), k=1..4, give the exact expected counts "
                       "(k-1 sign changes, k half-waves each)");

        std::srand(12345);
        gsVector<real_t> noise(nSelf);
        for (index_t i = 0; i != nSelf; ++i)
            noise[i] = 1e-14 * (2.0*std::rand()/RAND_MAX - 1.0);
        index_t hwNoise = halfWaveCount(noise, halfWaveTol, 1.0);
        report(hwNoise == 0, "self-test: round-off-level noise (1e-14 relative to a unit reference "
                              "scale) reads as flat (0 half-waves)");
    }

    // ---------------------- Geometry (mirrors example_ShearExploration) ----
    gsStopwatch clock; clock.restart();

    gsMultiPatch<> mp = Rectangle(aDim,bDim);
    for (index_t i = 0; i < numElevate; ++i) mp.patch(0).degreeElevate();
    for (index_t i = 0; i < numHref;    ++i) mp.patch(0).uniformRefine();

    gsMultiBasis<> dbasis(mp);
    const index_t n0 = static_cast<index_t>(mp.patch(0).basis().component(0).size());
    const index_t n1 = static_cast<index_t>(mp.patch(0).basis().component(1).size());
    gsInfo << "Basis (patch 0): "<< mp.patch(0).basis()
           <<"  [resolvable half-wave bound: "<<n0<<" along u, "<<n1<<" along v]\n";

    // ---------------------- Boundary conditions (pure shear) ---------------
    // South edge fully clamped; north edge sheared in x (driven by lambda),
    // held in y and z. No geometric imperfection anywhere. lambda is the
    // prescribed north-edge x-displacement in physical units; over width
    // bDim=1 that makes lambda numerically equal to the engineering shear
    // strain, and directly comparable to the ALM's lambda without conversion.
    gsConstantFunction<> displ(0.0,3);
    gsConstantFunction<> displ_const(0.0,3);

    gsBoundaryConditions<> BCs;
    BCs.setGeoMap(mp);
    BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 0);
    BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1);
    BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2);
    BCs.addCondition(boundary::north, condition_type::dirichlet, &displ,       0, false, 0);
    BCs.addCondition(boundary::north, condition_type::dirichlet, &displ_const, 0, false, 1);
    BCs.addCondition(boundary::north, condition_type::dirichlet, 0,            0, false, 2);

    gsVector<> tmp(3); tmp.setZero();
    gsConstantFunction<> force(tmp,3);
    gsFunctionExpr<> t (std::to_string(thickness),   3);
    gsFunctionExpr<> E (std::to_string(E_modulus),   3);
    gsFunctionExpr<> nu(std::to_string(PoissonRatio),3);
    gsFunctionExpr<> rho(std::to_string(Density),    3);
    gsMaterialMatrixLinear<3,real_t> materialMatrix(mp,t,E,nu,rho);

    gsMultiPatch<> mp_def = mp;
    gsThinShellAssemblerBase<real_t>* assembler =
        new gsThinShellAssembler<3,real_t,true>(mp,dbasis,BCs,force,&materialMatrix);

    gsInfo<<"Setup / assembly preparation: "<<clock.stop()<<" s\n";

    // ---------------------- Pre-buckling state, K_L ------------------------
    // Assemble at unit prescribed shear from the undeformed configuration:
    // this is a LINEAR system (no constructSolution call precedes it), so
    // assembler->matrix() is K_L and assembler->rhs() is the standard
    // Dirichlet-elimination right-hand side gsExprAssembler forms for that
    // inhomogeneous BC (zero external force here, so it is entirely the
    // Dirichlet lift). Which SIGN of K_L^-1*rhs() is the physical free-DoF
    // equilibrium response is confirmed below, not assumed (see the file
    // header note on the ALM's differing internal convention).
    clock.restart();
    displ.setValue(1.0,3);
    assembler->updateBCs(BCs);
    assembler->assemble();
    gsSparseMatrix<> K_L = assembler->matrix();
    gsVector<>       rhs1 = assembler->rhs();
    gsInfo<<"K_L assembly: "<<clock.stop()<<" s ("<<K_L.rows()<<" free DoFs)\n";
    // At -r 0 -e 0 every control point of the bilinear patch sits on the
    // clamped south edge or the fully-prescribed north edge, so the free-DoF
    // system is 0x0; the dense eigensolvers below have no defined behaviour
    // on an empty matrix (Eigen's bounds asserts are compiled out in release
    // builds, so this reads as a bare crash rather than a diagnosable error
    // without the check below).
    GISMO_ENSURE(K_L.rows() > 0, "no free degrees of freedom at -r "<<numHref<<" -e "<<numElevate
                                  <<"; refine or elevate further");

    gsSparseSolver<real_t>::SimplicialLDLT KLsolver;
    KLsolver.compute(K_L);
    GISMO_ENSURE(KLsolver.info()==gsEigen::Success, "K_L factorization failed");
    gsVector<> uL1cand = KLsolver.solve(rhs1); // sign resolved by the measurement below

    // ---------------------- K_L soundness ------------------------------------
    gsInfo << "\n-- K_L soundness --\n";
    gsEigen::SelfAdjointEigenSolver<gsMatrix<real_t>::Base> saes;
    saes.compute(K_L.toDense(), gsEigen::EigenvaluesOnly);
    const real_t lamMin = saes.eigenvalues().minCoeff();
    const real_t lamMax = saes.eigenvalues().maxCoeff();
    gsInfo << "  lambda_min(K_L) = "<<lamMin<<"   lambda_max(K_L) = "<<lamMax
           <<"   ratio = "<<lamMax/lamMin<<"\n";
    report(lamMin > 0.0, "K_L is SPD (lambda_min(K_L) > 0; rules out rigid-body / spurious "
                          "zero modes in the free-DoF space)");

    {
        gsMultiPatch<> mp_zero = mp_def;
        gsVector<> zero(K_L.rows()); zero.setZero();
        displ.setValue(0.0,3);
        assembler->updateBCs(BCs);
        assembler->constructSolution(zero, mp_zero);
        assembler->assembleMatrix(mp_zero);
        gsSparseMatrix<> K0 = assembler->matrix();
        const real_t rel = matInfNorm(K0 - K_L) / matInfNorm(K_L);
        gsInfo << "  ||K_NL(0) - K_L||_inf / ||K_L||_inf = "<<rel<<"\n";
        if (rel > 1e-10) softReport("K_NL(0) and K_L disagree by relative infinity-norm "+std::to_string(rel));
    }

    // ---------------------- pre-buckling sign check --------------------------
    // The nonlinear residual at a small step discriminates the sign: the
    // correct sign leaves a residual that is higher order in eps (the free
    // DoFs satisfy the LINEAR equilibrium exactly, so only the O(eps^2)
    // geometric nonlinearity remains), the wrong sign leaves an O(eps)
    // mismatch of the LINEAR equilibrium itself.
    gsInfo << "\n-- pre-buckling sign check (eps=1e-3) --\n";
    real_t signU = 1.0;
    {
        const real_t epsSign = 1e-3;
        auto residNorm = [&](real_t sgn) -> real_t
        {
            gsMultiPatch<> mp_r = mp_def;
            displ.setValue(epsSign,3);
            assembler->updateBCs(BCs);
            assembler->constructSolution(sgn*epsSign*uL1cand, mp_r);
            assembler->assembleVector(mp_r);
            return assembler->rhs().norm();
        };
        const real_t rPlus  = residNorm(+1.0);
        const real_t rMinus = residNorm(-1.0);
        gsInfo << "  ||r(+eps*uL1)|| = "<<rPlus<<"   ||r(-eps*uL1)|| = "<<rMinus
               <<"   ratio(larger/smaller) = "<<std::max(rPlus,rMinus)/std::min(rPlus,rMinus)<<"\n";
        signU = (rPlus < rMinus) ? 1.0 : -1.0;
        gsInfo << "  winning sign: "<<(signU>0?"+":"-")
               <<" (expected discrimination ratio ~ 1/eps = "<<1.0/epsSign<<")\n";
        const real_t ratio = std::max(rPlus,rMinus)/std::min(rPlus,rMinus);
        if (!report(ratio > 0.1/epsSign, "the winning sign's residual is smaller by roughly the "
                                          "expected 1/eps order of magnitude"))
        {
            gsInfo << "  pre-buckling sign is NOT confirmed by measurement; stopping "
                      "rather than proceeding on an unverified sign.\n";
            delete assembler;
            return EXIT_FAILURE;
        }
    }
    const gsVector<> uL1 = signU * uL1cand;
    displ.setValue(0.0,3);
    assembler->updateBCs(BCs);

    // ---------------------- Geometric stiffness builders --------------------
    // K_NL(w) at prescribed shear `scale`, free-DoF response `w` (Dirichlet
    // part of the state is set by displ, NOT by w -- constructSolution takes
    // fixed DoFs from the BCs).
    auto buildKNL = [&](real_t scale, const gsVector<>& w) -> gsSparseMatrix<>
    {
        gsMultiPatch<> mp_w = mp_def;
        displ.setValue(scale,3);
        assembler->updateBCs(BCs);
        assembler->constructSolution(w, mp_w);
        assembler->assembleMatrix(mp_w);
        return assembler->matrix();
    };

    // ---------------------- Eigenproblem at one eps -------------------------
    // Solves K_G(eps) v = nu K_L v (K_L definite, factored; see file header),
    // returns ascending eigenvalues/eigenvectors and lambda = -1/nu.
    auto solveEigen = [&](const gsSparseMatrix<>& K_G, gsVector<real_t>& nuOut, gsMatrix<real_t>& Vout) -> bool
    {
        gsEigen::GeneralizedSelfAdjointEigenSolver<gsMatrix<real_t>::Base> es;
        es.compute(K_G.toDense(), K_L.toDense());
        nuOut = es.eigenvalues();
        Vout  = es.eigenvectors();
        return es.info()==gsEigen::Success;
    };

    // Runs the full selection (positive + negative side walks) at a given
    // (nu,V) spectrum, returns the selected positive-lambda mode's info, or
    // an idx==-1 sentinel if none passed.
    const real_t nuFloorRel = 1e-8;

    // ---------------------- eps sweep: convergence of the contamination ------
    gsInfo << "\n-- eps sweep: lambda_lin(eps), forward difference K_G(eps) --\n";
    gsInfo << std::setw(10)<<std::left<<"eps"
           << std::setw(16)<<std::left<<"lambda+(eps)"
           << std::setw(16)<<std::left<<"lambda-(eps)"
           << std::setw(16)<<std::left<<"asymmetry"<<"\n";
    std::vector<real_t> lambdaPlusSweep, lambdaMinusSweep;
    for (real_t eps : epsList)
    {
        gsSparseMatrix<> K_NLp = buildKNL(eps, eps*uL1);
        gsSparseMatrix<> K_G   = (K_NLp - K_L) / eps;
        gsVector<real_t> nuVec; gsMatrix<real_t> V;
        bool ok = solveEigen(K_G, nuVec, V);
        GISMO_ENSURE(ok, "eigensolve failed in eps sweep at eps="<<eps);
        const real_t nuFloor = nuFloorRel * nuVec.cwiseAbs().maxCoeff();

        std::ostringstream mute; // silence the per-mode walk table during the sweep
        std::streambuf* old = gsInfo.rdbuf(mute.rdbuf());
        index_t iPos = walkSelect(nuVec.size(), true,  nuVec, V, K_L, K_G, nuFloor,
                                   assembler, BCs, displ, mp, halfWaveTol, 30);
        index_t iNeg = walkSelect(nuVec.size(), false, nuVec, V, K_L, K_G, nuFloor,
                                   assembler, BCs, displ, mp, halfWaveTol, 30);
        gsInfo.rdbuf(old);

        const real_t lp = (iPos>=0) ? -1.0/nuVec[iPos] : std::numeric_limits<real_t>::quiet_NaN();
        const real_t ln = (iNeg>=0) ? -1.0/nuVec[iNeg] : std::numeric_limits<real_t>::quiet_NaN();
        lambdaPlusSweep.push_back(lp);
        lambdaMinusSweep.push_back(ln);
        const real_t asym = (iPos>=0 && iNeg>=0) ? (std::abs(ln)-std::abs(lp))/std::abs(lp)
                                                  : std::numeric_limits<real_t>::quiet_NaN();
        gsInfo << std::setw(10)<<std::left<<eps
               << std::setw(16)<<std::left<<lp
               << std::setw(16)<<std::left<<ln
               << std::setw(16)<<std::left<<asym<<"\n";
    }

    // ---------------------- primary eps: full analysis ----------------------
    gsInfo << "\n-- primary eps = "<<epsPrimary<<" : forward K_G --\n";
    gsSparseMatrix<> K_NLp = buildKNL(epsPrimary, epsPrimary*uL1);
    gsSparseMatrix<> K_G   = (K_NLp - K_L) / epsPrimary;

    gsSparseMatrix<> K_NLm = buildKNL(-epsPrimary, -epsPrimary*uL1);
    gsSparseMatrix<> K_Gcd = (K_NLp - K_NLm) / (2.0*epsPrimary);

    gsVector<real_t> nuVec; gsMatrix<real_t> V;
    bool eigOk = solveEigen(K_G, nuVec, V);
    report(eigOk, "GeneralizedSelfAdjointEigenSolver::info() == Success on K_G(eps_primary), K_L");
    const real_t nuFloor = nuFloorRel * nuVec.cwiseAbs().maxCoeff();
    const index_t nModes = nuVec.size();

    gsVector<real_t> nuVecCd; gsMatrix<real_t> Vcd;
    bool eigOkCd = solveEigen(K_Gcd, nuVecCd, Vcd);
    report(eigOkCd, "GeneralizedSelfAdjointEigenSolver::info() == Success on K_G_cd(eps_primary), K_L");

    gsInfo << "\n-- spectrum walk, positive-lambda side (idx=0,1,2,... ascending nu) --\n";
    gsInfo << std::setw(5) <<std::left<<"idx"
           << std::setw(16)<<std::left<<"lambda"
           << std::setw(16)<<std::left<<"nu"
           << std::setw(14)<<std::left<<"residual"
           << std::setw(14)<<std::left<<"vz/vxy"
           << std::setw(10)<<std::left<<"signChg"
           << std::setw(10)<<std::left<<"halfWave"<<"\n";
    if (critMode>=0)
        GISMO_ENSURE(critMode < nModes, "critMode ("<<critMode<<") is out of range [0,"<<nModes<<")");
    index_t iPos = (critMode>=0) ? critMode
                 : walkSelect(nModes, true, nuVec, V, K_L, K_G, nuFloor,
                              assembler, BCs, displ, mp, halfWaveTol, std::max<index_t>(5,20),
                              /*skipUpTo=*/-1, /*minPrint=*/5);
    GISMO_ENSURE(iPos>=0, "no positive-lambda mode passed the selection filter");
    ModeInfo mode1 = evalMode(iPos, nuVec, V, K_L, K_G, nuFloor, assembler, BCs, displ, mp, halfWaveTol);
    if (critMode>=0) printModeRow(mode1, "  FORCED (--critMode)");

    gsInfo << "\n-- spectrum walk, negative-lambda side (idx=N-1,N-2,... descending nu) --\n";
    gsInfo << std::setw(5) <<std::left<<"idx"
           << std::setw(16)<<std::left<<"lambda"
           << std::setw(16)<<std::left<<"nu"
           << std::setw(14)<<std::left<<"residual"
           << std::setw(14)<<std::left<<"vz/vxy"
           << std::setw(10)<<std::left<<"signChg"
           << std::setw(10)<<std::left<<"halfWave"<<"\n";
    index_t iNeg = walkSelect(nModes, false, nuVec, V, K_L, K_G, nuFloor,
                               assembler, BCs, displ, mp, halfWaveTol, std::max<index_t>(5,20),
                               /*skipUpTo=*/-1, /*minPrint=*/5);
    GISMO_ENSURE(iNeg>=0, "no negative-lambda mode passed the selection filter");
    ModeInfo modeNeg = evalMode(iNeg, nuVec, V, K_L, K_G, nuFloor, assembler, BCs, displ, mp, halfWaveTol);

    gsInfo << "\n-- selected critical mode (positive side) --\n";
    gsInfo << "  lambda_lin = "<<mode1.lambda<<"   nu = "<<mode1.nu
           <<"   residual = "<<mode1.residual<<"   vz/vxy = "<<mode1.vz/(mode1.vxy+1e-30)<<"\n";
    report(mode1.outOfPlane, "selected critical mode is out-of-plane dominated (vz > vxy)");
    report(mode1.residual < 1e-8, "selected critical mode singularity residual "
                                   "||(K_L + lambda K_G)v||/||K_L v|| < 1e-8");

    // ---------------------- symmetry of the +/- lambda spectrum --------------
    // The rectangle and its BCs are symmetric under x -> aDim-x, which maps
    // shear +lambda to -lambda, so the true linearized critical loads satisfy
    // |lambda+| = |lambda-| exactly; the measured asymmetry is therefore a
    // second, independent estimate of the eps-contamination in K_G(eps).
    gsInfo << "\n-- lambda spectrum symmetry (primary eps) --\n";
    const real_t asymPrimary = (std::abs(modeNeg.lambda)-std::abs(mode1.lambda))/std::abs(mode1.lambda);
    gsInfo << "  lambda+ = "<<mode1.lambda<<"   lambda- = "<<modeNeg.lambda
           <<"   (|lambda-|-|lambda+|)/|lambda+| = "<<asymPrimary<<"\n";

    // ---------------------- mode gap ratio ------------------------------------
    // Scale-free (a ratio of two loads under the SAME eps-contamination), so it
    // cancels much of the linearization error a single absolute load carries.
    gsInfo << "\n-- mode gap ratio (positive side, continuing the walk past the selected mode) --\n";
    index_t iPos2 = walkSelect(nModes, true, nuVec, V, K_L, K_G, nuFloor,
                                assembler, BCs, displ, mp, halfWaveTol,
                                std::max<index_t>(5,20), /*skipUpTo=*/iPos);
    real_t gapRatio = std::numeric_limits<real_t>::quiet_NaN();
    if (iPos2>=0)
    {
        ModeInfo mode2 = evalMode(iPos2, nuVec, V, K_L, K_G, nuFloor, assembler, BCs, displ, mp, halfWaveTol);
        gapRatio = mode2.lambda/mode1.lambda;
        gsInfo << "  lambda1 = "<<mode1.lambda<<"   lambda2 = "<<mode2.lambda
               <<"   lambda2/lambda1 = "<<gapRatio<<"\n";
    }
    else
        gsInfo << "  no second filter-passing positive-lambda mode found within the walk budget\n";

    // ---------------------- central-difference control -----------------------
    gsInfo << "\n-- central-difference control at eps_primary --\n";
    {
        const real_t nuFloorCd = nuFloorRel * nuVecCd.cwiseAbs().maxCoeff();
        index_t iPosCd = walkSelect(nuVecCd.size(), true, nuVecCd, Vcd, K_L, K_Gcd, nuFloorCd,
                                     assembler, BCs, displ, mp, halfWaveTol, std::max<index_t>(5,20));
        const real_t lambdaCd = (iPosCd>=0) ? -1.0/nuVecCd[iPosCd] : std::numeric_limits<real_t>::quiet_NaN();
        gsInfo << "  lambda_lin (forward, eps="<<epsPrimary<<")          = "<<mode1.lambda<<"\n";
        gsInfo << "  lambda_lin_cd (central diff, eps="<<epsPrimary<<")  = "<<lambdaCd<<"\n";
        gsInfo << "  |forward - cd| / cd = "<<std::abs(mode1.lambda-lambdaCd)/lambdaCd<<"\n";
    }

    // ---------------------- half-wave detail for the critical mode -----------
    gsInfo << "\n-- critical mode half-wave count, both directions, 3 stations, "
           <<profileSamples<<" samples --\n";
    gsInfo << "  (resolvable bound: "<<n0<<" along u, "<<n1<<" along v; a count at or near "
              "that bound is a discretization artifact)\n";
    {
        gsMultiPatch<> mp_mode = mp;
        displ.setValue(0.0,3); assembler->updateBCs(BCs);
        assembler->constructSolution(V.col(iPos), mp_mode);

        gsInfo << std::setw(10)<<std::left<<"direction"
               << std::setw(10)<<std::left<<"station"
               << std::setw(10)<<std::left<<"signChg"
               << std::setw(10)<<std::left<<"halfWave"<<"\n";
        real_t critHalfWaveProfileScale = mode1.vz;
        for (index_t dir = 0; dir < 2; ++dir)
            for (real_t station : {0.25, 0.5, 0.75})
            {
                gsVector<real_t> prof = zProfile(mp_mode, dir, station, profileSamples);
                index_t sc = countSignChanges(prof, halfWaveTol);
                index_t hw = halfWaveCount(prof, halfWaveTol, critHalfWaveProfileScale);
                gsInfo << std::setw(10)<<std::left<<(dir==0?"u":"v")
                       << std::setw(10)<<std::left<<station
                       << std::setw(10)<<std::left<<sc
                       << std::setw(10)<<std::left<<hw<<"\n";
                if (dir==0 && station==0.5)
                    report(hw > 0, "critical mode profile (u-direction, v=0.5) is not flat "
                                    "(half-wave count > 0)");
            }

        if (plot)
        {
            std::string dirname = "ShearBucklingResults";
            GISMO_ENSURE(gsFileManager::mkdir(dirname), "failed to create output directory "+dirname);
            gsMultiPatch<> deformation = mp_mode;
            deformation.patch(0).coefs() -= mp.patch(0).coefs();
            gsField<> modeField(mp, deformation);
            gsWriteParaview<>(modeField, dirname+"/criticalMode", 1000);
        }
    }

    // ---------------------- comparison to the ALM's certified crossing -------
    gsInfo << "\n-- comparison to ALM lambda* --\n";
    if (lambdaOracle >= 0.0)
    {
        const real_t relErr = (mode1.lambda - lambdaOracle)/lambdaOracle;
        gsInfo << "  lambda_lin = "<<mode1.lambda<<"   lambda* (ALM) = "<<lambdaOracle
               <<"   signed relative difference = "<<relErr<<"\n";
        report(std::abs(relErr) <= lambdaBand,
               "lambda_lin agrees with the ALM's lambda* within the pre-registered band "
               "("+std::to_string(lambdaBand)+")");
    }
    else
        gsInfo << "  [SKIP] --lambdaOracle not given (sentinel -1): comparison not evaluated\n";

    gsInfo << "\nfree DoFs = "<<K_L.rows()<<"\n";
    delete assembler;
    return allOk ? EXIT_SUCCESS : EXIT_FAILURE;
}
#else//gsKLShell_ENABLED
int main(int /*argc*/, char ** /*argv*/)
{
    gsWarn<<"G+Smo is not compiled with the gsKLShell module.";
    return EXIT_FAILURE;
}
#endif

template <class T>
gsMultiPatch<T> Rectangle(T L, T B)
{
  // Single-patch bi-linear rectangle [0,L] x [0,B] embedded in 3D (z=0).
  int dim = 3;
  gsKnotVector<> kv0; kv0.initUniform(0,1,0,2,1);
  gsKnotVector<> kv1; kv1.initUniform(0,1,0,2,1);

  gsTensorBSplineBasis<2,T> basis(kv0,kv1);

  gsMatrix<> coefs(basis.size(),dim);
  size_t len0 = basis.component(0).size();
  size_t len1 = basis.component(1).size();
  gsVector<> coefvec0(len0); coefvec0.setLinSpaced(len0,0.0,L);
  gsVector<> coefvec1(len1); coefvec1.setLinSpaced(len1,0.0,B);

  coefs.col(2).setZero();
  gsVector<> temp(len0); temp.setOnes();
  for (size_t k = 0; k < len1; k++)
  {
    coefs.col(0).segment(k*len0,len0) = coefvec0;
    coefs.col(1).segment(k*len0,len0) = temp*coefvec1.at(k);
  }

  gsTensorBSpline<2,T> shape(basis,coefs);
  gsMultiPatch<T> mp;
  mp.addPatch(shape);
  mp.addAutoBoundaries();
  return mp;
}
