/** @file benchmark_Pillow.cpp

    @brief Benchmark for the inflated pillow using DR and Newton-iterations

    Quarter model of an inflated square airbag (Verhelst et al. (2026), Fig. 3 and Table 1;
    Zhang et al. (2026), Table 1). Compressible Neo-Hookean membrane, E = 588 MPa, nu = 0.4,
    t = 0.1 mm, pressure 5000 Pa, full-square diagonal \a --diag (1.2 m: Verhelst et al. (2026);
    1.22 m: Zhang et al. (2026)).

    Two load stages (\a -N 1): load factor 0 applies the in-plane edge load of 5000 N/m on the free
    edges (pressure 0), load factor 1 applies the pressure (edge load 0).

    The membrane can be run without a wrinkling model, with the tension-field model
    (gsMaterialMatrixTFT, \a --TFT; options \a --eta and \a --etatangent), or with the split formulation (gsMaterialMatrixTFTSplit,
    \a --tftsplit; options \a --eta, \a --exact, \a --midsurface, \a --s2crit, \a --slackmult).

    Every attempted load step appends one row to \a <out>/steps.csv (see WrinkleMonitor), followed by
    the pressure, the displacements at the four corners and the two free-edge midpoints of the
    quarter, and the principal mid-surface stresses at M (columns documented in the
    steps.csv comment in main). \a --nx N replaces the \a -r refined bilinear square by a uniform N x N tensor B-spline
    patch of degree 1+\a -e (C^{p-1}).

    Solver strategies:
    - \a --DR: dynamic relaxation (\a -T, \a --tolE); \a --NR: Newton (\a --tolNR, \a --maxitNR,
      \a --solver, \a --linesearch). Both together run DR as predictor of Newton.
    - \a --warmstart: each load step starts from the last converged displacement instead of zero.
    - \a --drfirst: with \a --DR \a --NR \a --warmstart, the DR predictor is used on the first
      load step only (the flat, stress-free state has no out-of-plane stiffness, so Newton
      cannot start there); every later load step is solved by Newton alone.

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

/// Solver data of one attempted load step, written to the first columns of a steps.csv row
struct StepInfo
{
    index_t step, attempt;
    real_t load;
    bool converged;
    int statusNR, statusDR;   ///< static_cast<int>(gsStatus), -1 if the solver did not run
    index_t itNR, itDR;       ///< iterations, -1 if the solver is not enabled
    index_t bisections, failures, singular;
    bool usesNR;
    index_t splitCounters[4]; ///< theta fallbacks, TFT failures, exact fallbacks, exact reclassified; -1 without Split
    real_t time;
};

/// Statistics over the sample points on one deformed configuration
struct PointStats
{
    PointStats()
    :
    nPoints(-1), nTaut(-1), nWrinkled(-1), nSlack(-1), stateChanges(-1),
    ghMax(std::numeric_limits<real_t>::quiet_NaN()), ghP95(std::numeric_limits<real_t>::quiet_NaN()),
    n4SNonpos(-1), tftThetaMismatch(-1), tftThetaZero(-1), tftStateMismatch(-1)
    {}

    index_t nPoints, nTaut, nWrinkled, nSlack, stateChanges;
    real_t ghMax, ghP95;
    index_t n4SNonpos, tftThetaMismatch, tftThetaZero, tftStateMismatch;
};

/**
 * @brief Wrinkling diagnostics: the wrinkling state at the quadrature points of the assembler, and
 *        its change between Newton iterations and between load steps.
 *
 * The sample points are the Gauss points of the assembler (all elements) times the fibres z = +-1/2
 * (z = 0 for a mid-surface integrated material). The state is judged by an oracle, a
 * gsMaterialMatrixTFTSplit with Eta = 0, Exact = false and S2crit = 1e300 wrapped (without ownership)
 * around the base material. The oracle is never handed to the assembler, so the statistics are
 * comparable across models and its counters do not mix with those of the run material. State
 * convention of the TensionField output: 1 = taut, 0 = wrinkled, -1 = slack.
 *
 * Complexity: one evaluation of the oracle at all sample points per call, i.e. of the order of one
 * residual assembly.
 */
class WrinkleMonitor
{
public:
    /// \param oracle        wrinkling-state oracle (not owned)
    /// \param nodes         parametric sample points (dim x n)
    /// \param midsurface    sample z = 0 only
    /// \param level         0: no point statistics, 1: per-step statistics, 2: also per-iteration chatter
    WrinkleMonitor(MaterialBase * oracle, const gsMatrix<> & nodes, bool midsurface, index_t level)
    :
    m_oracle(oracle), m_nodes(nodes), m_Z(midsurface ? 1 : 2, 1), m_level(level),
    m_hasAccepted(false), m_hasPrev(false), m_jacCalls(0), m_chatter(0)
    {
        if (midsurface)
            m_Z<<0;
        else
            m_Z<<-0.5,0.5;
    }

    /// Clears the per-solve iteration observations; call before every solve
    void beginSolve()
    {
        m_hasPrev = false;
        m_jacCalls = 0;
        m_chatter = 0;
    }

    /// Number of iterates observed since \ref beginSolve (-1 if the hook is off)
    index_t jacobianCalls() const { return m_level >= 2 ? m_jacCalls : -1; }

    /// Number of state flips between consecutive iterates of the current solve (-1 if the hook is off)
    index_t chatter() const { return m_level >= 2 ? m_chatter : -1; }

    /// Compares the state on \a def with the previous iterate of the same solve. Call after the tangent
    /// matrix was assembled on \a def, once per Newton iteration.
    void observeIterate(const gsMultiPatch<> & def)
    {
        if (m_level < 2)
            return;
        try
        {
            const gsMatrix<> state = evalMaterial<MaterialOutput::TensionField>(m_oracle,def,m_nodes,m_Z);
            if (m_hasPrev && state.cols()==m_prev.cols())
                for (index_t i = 0; i!=state.cols(); ++i)
                    if (state(0,i)!=m_prev(0,i))
                        ++m_chatter;
            m_prev = state;
            m_hasPrev = true;
            ++m_jacCalls;
        }
        catch (const std::exception & e)
        {
            gsWarn<<"Wrinkle monitor: iterate observation failed: "<<e.what()<<"\n";
        }
        catch (...)
        {
            gsWarn<<"Wrinkle monitor: iterate observation failed\n";
        }
    }

    /**
     * @brief Statistics on the configuration \a def.
     *
     * @param def    deformed configuration
     * @param run    material that the assembler uses
     * @param model  0: no wrinkling model, 1: gsMaterialMatrixTFT, 2: gsMaterialMatrixTFTSplit
     *
     * The statistics of the oracle state (n_taut, n_wrinkled, n_slack, state_changes) are available for
     * every model; those of \a run need a wrinkling model. Nothing is evaluated on level 0.
     */
    PointStats evaluate(const gsMultiPatch<> & def, MaterialBase * run, int model)
    {
        PointStats st;
        m_candidate.resize(0,0);
        if (m_level < 1)
            return st;
        try
        {
            const gsMatrix<> state = evalMaterial<MaterialOutput::TensionField>(m_oracle,def,m_nodes,m_Z);
            st.nPoints = state.cols();
            st.nTaut = st.nWrinkled = st.nSlack = 0;
            for (index_t i = 0; i!=state.cols(); ++i)
            {
                if      (state(0,i)==1) ++st.nTaut;
                else if (state(0,i)==0) ++st.nWrinkled;
                else if (state(0,i)==-1) ++st.nSlack;
            }
            st.stateChanges = 0;
            if (m_hasAccepted && m_accepted.cols()==state.cols())
                for (index_t i = 0; i!=state.cols(); ++i)
                    if (state(0,i)!=m_accepted(0,i))
                        ++st.stateChanges;
            m_candidate = state;

            if (model==0)
                return st;

            const gsMatrix<> gh = evalMaterial<MaterialOutput::GammaHat>(run,def,m_nodes,m_Z);
            std::vector<real_t> finite;
            for (index_t i = 0; i!=gh.cols(); ++i)
                if (std::isfinite(gh(0,i)))
                    finite.push_back(gh(0,i));
            if (!finite.empty())
            {
                std::sort(finite.begin(),finite.end());
                st.ghMax = finite.back();
                // nearest-rank percentile
                st.ghP95 = finite[(size_t)std::ceil(0.95*finite.size()) - 1];
            }

            const gsMatrix<> theta = evalMaterial<MaterialOutput::Theta>(run,def,m_nodes,m_Z);
            const gsMatrix<> S = evalMaterial<MaterialOutput::Generic>(run,def,m_nodes,m_Z);

            // Tension along the wrinkles: n4.S = m1^2 S11 + m2^2 S22 + 2 m1 m2 S12, m = (-sin(theta),cos(theta)).
            // With Eta > 0 this is the blended stress that the assembler integrates (for TFT with --etatangent the
            // integrated stress is the unblended tension-field stress).
            st.n4SNonpos = 0;
            for (index_t i = 0; i!=theta.cols(); ++i)
            {
                if (!std::isfinite(theta(0,i)))
                    continue;
                const real_t m1 = -std::sin(theta(0,i)), m2 = std::cos(theta(0,i));
                if (m1*m1*S(0,i) + m2*m2*S(1,i) + 2*m1*m2*S(2,i) <= 0)
                    ++st.n4SNonpos;
            }

            if (model==1)
            {
                // The oracle's theta is the closed-form angle in [0,pi); TFT returns the root of its own
                // iteration, which is only defined modulo pi, and silently falls back to 0.
                const gsMatrix<> thetaCf = evalMaterial<MaterialOutput::Theta>(m_oracle,def,m_nodes,m_Z);
                st.tftThetaMismatch = st.tftThetaZero = st.tftStateMismatch = 0;
                for (index_t i = 0; i!=theta.cols(); ++i)
                {
                    const bool f1 = std::isfinite(theta(0,i)), f2 = std::isfinite(thetaCf(0,i));
                    if (f1 != f2)
                        ++st.tftStateMismatch;
                    else if (f1)
                    {
                        if (std::abs(std::sin(theta(0,i)-thetaCf(0,i))) > 1e-8)
                            ++st.tftThetaMismatch;
                        if (theta(0,i)==0 && std::abs(std::sin(thetaCf(0,i))) > 1e-8)
                            ++st.tftThetaZero;
                    }
                }
            }
        }
        catch (const std::exception & e)
        {
            gsWarn<<"Wrinkle monitor: point statistics failed: "<<e.what()<<"\n";
        }
        catch (...)
        {
            gsWarn<<"Wrinkle monitor: point statistics failed\n";
        }
        return st;
    }

    /// Makes the state of the last \ref evaluate call the reference of \c state_changes
    void accept()
    {
        if (m_candidate.cols() > 0)
        {
            m_accepted = m_candidate;
            m_hasAccepted = true;
        }
    }

    /// Header of steps.csv (without the case-specific columns).
    ///
    /// Columns, in order (counts are integers, -1 marks "not available"):
    /// - step:            index of the accepted load step the attempt belongs to (-)
    /// - attempt:         index of the attempt in the run; increments on every solve (-)
    /// - load:            load factor (-)
    /// - converged:       1 if the solver chain converged, else 0
    /// - status_nr, status_dr: static_cast<int>(gsStatus) of the Newton / DR solver (0 = Success, 1 = NotConverged,
    ///                    2 = AssemblyError, 3 = SolverError); -1 if the solver did not run in this attempt
    /// - it_nr, it_dr:    iterations of the Newton / DR solver in this attempt; -1 if the solver is not enabled,
    ///                    0 if it is enabled but did not run
    /// - bisections:      cumulative number of load step halvings
    /// - failures:        cumulative number of non-converged attempts
    /// - singular:        cumulative number of attempts that ended with gsStatus::SolverError (singular tangent)
    /// - jac_calls:       number of tangent matrix assemblies in this attempt (Newton iterations + 1); -1 if not observed
    /// - chatter:         number of oracle state flips between consecutive tangent assemblies of this attempt,
    ///                    summed over all pairs of consecutive assemblies (a point that flips k times counts k
    ///                    times, so the value can exceed n_points); -1 if not observed
    /// - n_points:        number of sample points (Gauss points times fibres); -1 without point statistics
    /// - n_taut, n_wrinkled, n_slack: number of points in each state of the oracle
    /// - state_changes:   number of points whose oracle state differs from that of the last accepted step
    /// - gh_max, gh_p95:  maximum and nearest-rank 95th percentile of the wrinkle strain gamma_hat
    ///                    (= gamma G^{ab} n_a n_b, dimensionless) of the run material over wrinkled points;
    ///                    NaN without wrinkling model or wrinkled points
    /// - n4S_nonpos:      number of wrinkled points where the stress along the wrinkles, n4.S, of the run
    ///                    material is <= 0; -1 without wrinkling model
    /// - tft_theta_mismatch, tft_theta_zero, tft_state_mismatch: TFT runs only (else -1). Points where the
    ///                    wrinkle angle of TFT differs from the closed-form angle modulo pi (tolerance 1e-8
    ///                    in the sine), where TFT returned exactly 0 against a nonzero closed-form angle, and
    ///                    where only one of the two is wrinkled
    /// - split_theta_fb, split_tft_fail, split_exact_fb, split_exact_recl: Split runs only (else -1). Counters
    ///                    of gsMaterialMatrixTFTSplit over this attempt's solve: closed-form angle rejected by the
    ///                    isotropy guard, TFT root finder also failed, exact solve not converged or at a saddle,
    ///                    exact solve reclassified the point
    /// - time_s:          wall time of the solve (s)
    static std::string header()
    {
        return "step,attempt,load,converged,status_nr,status_dr,it_nr,it_dr,bisections,failures,singular,"
               "jac_calls,chatter,n_points,n_taut,n_wrinkled,n_slack,state_changes,gh_max,gh_p95,n4S_nonpos,"
               "tft_theta_mismatch,tft_theta_zero,tft_state_mismatch,"
               "split_theta_fb,split_tft_fail,split_exact_fb,split_exact_recl,time_s";
    }

    /// Writes the columns of \ref header for one attempt, without trailing separator or newline
    void writeColumns(std::ostream & os, const StepInfo & info, const PointStats & st) const
    {
        os<<info.step<<","<<info.attempt<<","<<info.load<<","<<(info.converged ? 1 : 0)<<","
          <<info.statusNR<<","<<info.statusDR<<","<<info.itNR<<","<<info.itDR<<","
          <<info.bisections<<","<<info.failures<<","<<info.singular<<","
          <<(info.usesNR ? jacobianCalls() : -1)<<","<<(info.usesNR ? chatter() : -1)<<","
          <<st.nPoints<<","<<st.nTaut<<","<<st.nWrinkled<<","<<st.nSlack<<","<<st.stateChanges<<","
          <<st.ghMax<<","<<st.ghP95<<","<<st.n4SNonpos<<","
          <<st.tftThetaMismatch<<","<<st.tftThetaZero<<","<<st.tftStateMismatch<<","
          <<info.splitCounters[0]<<","<<info.splitCounters[1]<<","<<info.splitCounters[2]<<","<<info.splitCounters[3]<<","
          <<info.time;
    }

private:
    MaterialBase * m_oracle;
    gsMatrix<> m_nodes, m_Z;
    index_t m_level;
    gsMatrix<> m_accepted, m_candidate, m_prev;
    bool m_hasAccepted, m_hasPrev;
    index_t m_jacCalls, m_chatter;
};

/// Builds the wrinkling-state oracle around \a base (not owned): a Split decorator with Eta = 0, no exact
/// solve and S2crit = 1e300, so that the TensionField output is the criterion state
template <bool linear>
gsMaterialMatrixBase<real_t>::uPtr makeOracle(MaterialBase * base, real_t slackMultiplier)
{
    gsMaterialMatrixBase<real_t>::uPtr oracle =
        memory::make_unique(new gsMaterialMatrixTFTSplit<3,real_t,linear>(base));
    gsOptionList o = oracle->options();
    o.setReal("Eta",0);
    o.setSwitch("Exact",false);
    o.setReal("S2crit",1e300);
    o.setReal("SlackMultiplier",slackMultiplier);
    oracle->setOptions(o);
    return oracle;
}

/// Sets the options of a Split material
void setSplitOptions(MaterialBase * split, real_t eta, bool exact, bool midsurface, real_t s2crit, real_t slackMultiplier)
{
    gsOptionList o = split->options();
    o.setReal("Eta",eta);
    o.setSwitch("Exact",exact);
    o.setSwitch("Midsurface",midsurface);
    o.setReal("S2crit",s2crit);
    o.setReal("SlackMultiplier",slackMultiplier);
    split->setOptions(o);
}

}// namespace

int main (int argc, char** argv)
{
    // Input options
    int numElevate    = 0;
    int numRefine     = 1;
    int nx            = 0;
    bool plot         = false;
    bool write        = false;
    bool stress       = false;
    bool split        = false;
    bool smooth       = false;
    bool DR  = false;
    bool NR  = false;
    int verbose = 0;
    bool membrane = false;

    real_t tolDR = 1e-1;

    index_t Compressibility = 1;
    index_t material = 1;
    index_t impl = 1; // 1= analytical, 2= generalized, 3= spectral
    bool TFT = false;

    index_t loadSteps     = 1e0;
    // Arc length method options
    real_t tol        = 1e-4;

    real_t alpha = 1.0;
    real_t damping = 0;

    bool tftsplit = false;
    real_t eta = 0;
    bool etaTangent = false;
    bool exact = false;
    bool midsurface = false;
    real_t s2crit = 0;
    real_t slackmult = 1e-6;
    index_t maxitNR = 100;
    std::string linSolver("LU");
    bool lineSearch = false;
    real_t tolE = 1e-5;
    real_t diag = 1.2;
    real_t thicknessOpt = -1;
    bool warmstart = false;
    bool drfirst = false;
    std::string outDir;
    index_t monitorLevel = 2;

    index_t testCase = 0;
    real_t maxLoad = 5e3;

    std::string wn("data.csv");

    gsCmdLine cmd("Example for an inflating pillow.");

    cmd.addInt("r","hRefine", "Number of dyadic h-refinement (bisection) steps to perform before solving", numRefine);
    cmd.addInt("","nx", "Uniform N x N quarter mesh of degree 1+e (replaces -r); 0: BSplineSquare refined by -r (default 0)", nx);
    cmd.addInt("e","degreeElevation", "Number of degree elevation steps to perform on the Geometry's basis before solving", numElevate);

    cmd.addInt( "N", "maxit", "maxit",  loadSteps );
    cmd.addInt( "t", "testCase", "testCase",  testCase );

    cmd.addReal( "a", "alpha", "alpha",  alpha );
    cmd.addReal( "c", "damping", "damping",  damping );

    cmd.addInt( "M", "Material", "Material law",  material );
    cmd.addInt( "C", "Compressibility", "1: compressible, 0: incompressible",  Compressibility );
    cmd.addInt( "I", "Implementation", "Implementation: 1= analytical, 2= generalized, 3= spectral",  impl );
    cmd.addSwitch( "TFT", "Use Tension-Field Theory",  TFT );

    cmd.addReal( "L", "maxLoad", "Maximum load",  maxLoad );

    cmd.addReal( "T", "tolDR", "Tolerance for the DR solver",  tolDR );
    cmd.addReal( "", "tolE", "Tolerance on the relative kinetic energy of the DR solver (default 1e-5)",  tolE );
    cmd.addReal( "", "tolNR", "Relative tolerance of the Newton solver, residual and update (default 1e-4)",  tol );
    cmd.addInt( "", "maxitNR", "Maximum number of Newton iterations (default 100)",  maxitNR );
    cmd.addSwitch( "linesearch", "Newton: backtracking line search on the residual norm (gsStaticNewton option LineSearch)", lineSearch );
    cmd.addString( "", "solver", "Linear solver of the Newton iterations: LU, CGDiagonal, SimplicialLDLT, ... (default LU)",  linSolver );

    cmd.addSwitch( "tftsplit", "Use the split tension-field formulation (gsMaterialMatrixTFTSplit)",  tftsplit );
    cmd.addReal( "", "eta", "TFT and Split: stress fraction retained in the wrinkled/slack state (default 0; needs --TFT or --tftsplit)",  eta );
    cmd.addSwitch( "etatangent", "TFT: apply --eta to the tangent only (stress stays the pure tension field)",  etaTangent );
    cmd.addSwitch( "exact", "Split: exact wrinkled-state solve for nonlinear materials (needs --tftsplit)",  exact );
    cmd.addSwitch( "midsurface", "Split: evaluate the material at the mid-surface only (needs --tftsplit)",  midsurface );
    cmd.addReal( "", "s2crit", "Split: minor principal stress below which a point is displayed as wrinkled (default 0; needs --tftsplit)",  s2crit );
    cmd.addReal( "", "slackmult", "Multiplier of the stiffness in the slack state, TFT and Split (default 1e-6)",  slackmult );

    cmd.addReal( "", "diag", "Diagonal of the full square, testCase 0 [m] (default 1.2)",  diag );
    cmd.addReal( "", "thickness", "Thickness [m]; negative: value of the test case (default -1)",  thicknessOpt );
    cmd.addSwitch( "warmstart", "Start every load step from the last converged displacement (default: from zero)",  warmstart );
    cmd.addSwitch( "drfirst", "With --DR --NR --warmstart: DR on the first load step only, Newton alone afterwards",  drfirst );
    cmd.addString( "", "out", "Output directory (default: ./Pillow_<r<r>|nx<N>>_e<e>-M<M>-c<C>[_TFT|_TFTSplit])",  outDir );
    cmd.addInt( "", "monitor", "Point diagnostics in steps.csv: 0: none; 1: per step; 2: also per Newton iteration (default 2)",  monitorLevel );

    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addSwitch("write", "Write data to file", write);
    cmd.addSwitch("stress", "Plot stress in ParaView format", stress);
    cmd.addSwitch("DR", "Use Dynamic Relaxation", DR);
    cmd.addSwitch("NR", "Use Newton Raphson", NR);
    cmd.addSwitch("split", "Split to a multi-patch", split);
    cmd.addSwitch("smooth", "Use a smooth basis (maximum regularity)", smooth);
    cmd.addInt("v","verbose", "0: no; 1: iteration output; 2: Full matrix and vector output", verbose);
    cmd.addSwitch("membrane", "Membrane element", membrane);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    if (!(eta >= 0 && eta <= 1)) { gsWarn << "--eta must be in [0,1], got " << eta << "\n"; return EXIT_FAILURE; }
    if (!(slackmult >= 0)) { gsWarn << "--slackmult must be non-negative, got " << slackmult << "\n"; return EXIT_FAILURE; }
    GISMO_ASSERT(NR || DR,"At least a Newton Raphson or a Dynamic Relaxation solver needs to be enabled. Run with option NR or DR.");
    if (TFT && tftsplit)
      GISMO_ERROR("--TFT and --tftsplit select different wrinkling models; use only one of them.");
    if (!TFT && !tftsplit && eta!=0)
      gsWarn<<"--eta is ignored without a wrinkling model.\n";
    if (!TFT && etaTangent)
      gsWarn<<"--etatangent only acts together with --TFT and is ignored.\n";
    if (!tftsplit && (exact || midsurface || s2crit!=0))
      gsWarn<<"--exact, --midsurface and --s2crit only act together with --tftsplit and are ignored.\n";
    if (nx<0)
      GISMO_ERROR("--nx must be non-negative.");
    if (nx>0)
    {
      // -r defaults to 1, so only an explicitly given non-zero -r is worth a warning
      bool rGiven = false;
      for (int i = 1; i < argc; ++i)
      {
        const std::string arg(argv[i]);
        if (arg=="-r" || arg.compare(0,9,"--hRefine")==0 || (arg.size()>2 && arg.compare(0,2,"-r")==0 && std::isdigit(static_cast<unsigned char>(arg[2]))))
          rGiven = true;
      }
      if (rGiven && numRefine!=0)
        gsWarn<<"--nx given: -r is ignored.\n";
    }
    if (drfirst && !(DR && NR && warmstart))
      GISMO_ERROR("--drfirst needs --DR, --NR and --warmstart.");

    /*
      Uniaxial tension of a square plate                            --- Validation settings: -L 1eX -l 1eX -M 14 -N 500 -r X -e X
      (bottom boundary fixed in Y, left boundary fixed in X, right boundary normal load)
    */
    real_t E_modulus = 0;
    real_t PoissonRatio = 0;
    real_t Ratio = 0;
    real_t thickness = 0;
    real_t Density    = 1e0;

    if (material==0) Compressibility = 0;

    gsMatrix<> box(2,2);
    if      (testCase == 0)
    {
      //Jarasjarungkiat2009,Diaby2006
      E_modulus = 588e6;
      PoissonRatio = 0.4;
      // thickness = 0.0006; //Jarasjarungkiat2009
      thickness = 1e-4; //Diaby2006
      box.col(0)<<0,0;
      real_t AM = diag/2.;
      real_t AB = math::sqrt(math::pow(AM,2)/2.);
      box.col(1)<<AB,AB;
    }
    else if (testCase == 1)
    {
      E_modulus = 2e9; // 2GPa
      PoissonRatio = 0.3;
      thickness = 0.1e-3;
      box.col(0)<<0,0;
      box.col(1)<<0.5,0.5;
    }
    if (thicknessOpt > 0)
      thickness = thicknessOpt;

    gsMultiPatch<> mp,mp_def;
    if (nx > 0)
    {
      // Uniform nx x nx tensor B-spline of degree 1+numElevate and regularity C^{p-1}. Greville control
      // points reproduce the affine map of the box exactly.
      const index_t degree = 1 + numElevate;
      gsKnotVector<> kx(0,1,nx-1,degree+1), ky(0,1,nx-1,degree+1);
      gsTensorBSplineBasis<2,real_t> tb(kx,ky);
      const gsMatrix<> anchors = tb.anchors();
      gsMatrix<> coefs(tb.size(),2);
      for (index_t i = 0; i != tb.size(); ++i)
      {
        coefs(i,0) = box(0,0) + anchors(0,i)*(box(0,1)-box(0,0));
        coefs(i,1) = box(1,0) + anchors(1,i)*(box(1,1)-box(1,0));
      }
      mp.addPatch(gsTensorBSpline<2,real_t>(tb,coefs));
      mp.embed(3);
    }
    else
    {
      mp.addPatch(gsNurbsCreator<>::BSplineSquare(box));
      mp.embed(3);

      // p-refine
      for(index_t i = 0; i< numElevate; ++i)
        mp.patch(0).degreeElevate();    // Elevate the degree

      // h-refine
      for(index_t i = 0; i< numRefine; ++i)
        mp.patch(0).uniformRefine();
    }

    gsInfo<<"mu = "<<E_modulus / (2 * (1 + PoissonRatio))<<"\n";

    gsMultiBasis<> dbasis(mp,true);
    gsInfo<<"Basis (patch 0): "<< mp.patch(0).basis() << "\n";
    mp_def = mp;

    gsVector<> neuDataX(3);
    neuDataX<<maxLoad,0,0;
    gsConstantFunction<> neuX(neuDataX,3);
    gsVector<> neuDataY(3);
    neuDataY<<0,-maxLoad,0;
    gsConstantFunction<> neuY(neuDataY,3);

    // Boundary conditions
    gsBoundaryConditions<> BCs;


    BCs.setGeoMap(mp);

    BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2);
    BCs.addCondition(boundary::south, condition_type::neumann, &neuY);
    BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 2);
    BCs.addCondition(boundary::east, condition_type::neumann, &neuX);

    BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 );
    BCs.addCondition(boundary::north, condition_type::clamped, 0, 0, false, 0 );
    BCs.addCondition(boundary::north, condition_type::clamped, 0, 0, false, 2 );

    // Symmetry in y-direction:
    BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 );
    BCs.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 1 );
    BCs.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 2 );

    BCs.setGeoMap(mp);

    // Pressure
    real_t pressure = maxLoad;
    // real_t pressure =0;

    std::string dirname = ".";
    dirname = dirname + "/Pillow_" + (nx>0 ? "nx" + std::to_string(nx) : "r" + std::to_string(numRefine)) + "_e" + std::to_string(numElevate) + "-M" + std::to_string(material) + "-c" + std::to_string(Compressibility);
    if (TFT)
      dirname = dirname + "_TFT";
    if (tftsplit)
      dirname = dirname + "_TFTSplit";
    if (!outDir.empty())
      dirname = outDir;
    std::string output =  "solution";

    gsFileManager::mkdir(dirname);

    // plot geometry
    if (plot)
      gsWriteParaview(mp,"mp",1000,true);

    // Linear isotropic material model
    gsVector<> forceVec(3);
    forceVec.setZero();

    gsConstantFunction<> force(forceVec,3);
    gsConstantFunction<> pressFun(pressure,3);
    gsConstantFunction<> t(thickness,3);
    gsConstantFunction<> E(E_modulus,3);
    gsConstantFunction<> nu(PoissonRatio,3);
    gsConstantFunction<> rho(Density,3);
    gsConstantFunction<> ratio(Ratio,3);

    gsMaterialMatrixBase<real_t>::uPtr materialMatrix;
    gsMaterialMatrixBase<real_t>::uPtr materialMatrixTFT;
    gsMaterialMatrixBase<real_t>::uPtr materialMatrixTFTSplit;
    gsMaterialMatrixBase<real_t>::uPtr oracle;
    gsOptionList options;

    std::vector<gsFunctionSet<>*> parameters;
    parameters.resize(2);
    parameters[0] = &E;
    parameters[1] = &nu;
    if      (material==0 && impl==1)
    {
      parameters.resize(2);
      options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
      options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
      materialMatrix = getMaterialMatrix<3,real_t>(mp,t,parameters,rho,options);
      materialMatrixTFT = memory::make_unique(new gsMaterialMatrixTFT<3,real_t,true>(*materialMatrix));
      if (tftsplit)
        materialMatrixTFTSplit = memory::make_unique(new gsMaterialMatrixTFTSplit<3,real_t,true>(*materialMatrix));
      oracle = makeOracle<true>(materialMatrix.get(),slackmult);
      // materialMatrixTFT = new gsMaterialMatrixTFT<3,real_t,false>(materialMatri));
    }
    else if (material==1 || material==2)
    {
      options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",material);
      options.addSwitch("Compressibility","Compressibility: (false): Imcompressible | (true): Compressible",Compressibility);
      options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",impl);
      materialMatrix = getMaterialMatrix<3,real_t>(mp,t,parameters,rho,options);
      materialMatrixTFT = memory::make_unique(new gsMaterialMatrixTFT<3,real_t,false>(*materialMatrix));
      if (tftsplit)
        materialMatrixTFTSplit = memory::make_unique(new gsMaterialMatrixTFTSplit<3,real_t,false>(*materialMatrix));
      oracle = makeOracle<false>(materialMatrix.get(),slackmult);
      // dynamic_cast<gsMaterialMatrixTFT<3,real_t,false> *>(materialMatrixTFT)->updateDeformed(&mp_def);
    }
    else
      GISMO_ERROR("Other materials not implemented in this example");

    // gsMaterialMatrixBase<real_t>* materialMatrix = new gsMaterialMatrixLinear<3,real_t>(mp,t,E,nu,rho);
    // gsMaterialMatrixTFT<3,real_t,false> * materialMatrixTFT = new gsMaterialMatrixTFT<3,real_t,false>(static_cast<gsMaterialMatrixBaseDim<3,real_t> * >(materialMatrix));

    if (TFT)
    {
      materialMatrixTFT->options().setReal("SlackMultiplier",slackmult);
      materialMatrixTFT->options().setReal("Eta",eta);
      materialMatrixTFT->options().setSwitch("EtaTangentOnly",etaTangent);
    }
    if (tftsplit)
      setSplitOptions(materialMatrixTFTSplit.get(),eta,exact,midsurface,s2crit,slackmult);

    // Material that the assembler integrates, and its wrinkling model (0: none, 1: TFT, 2: Split)
    const bool wrapped = TFT || tftsplit;
    const int model = tftsplit ? 2 : (TFT ? 1 : 0);
    gsMaterialMatrixBase<real_t>::uPtr & materialMatrixWrapped = tftsplit ? materialMatrixTFTSplit : materialMatrixTFT;
    // materialMatrixTFT->options().setSwitch("Explicit",true);

    gsThinShellAssemblerBase<real_t>* assembler;
    if (!membrane && wrapped)
      assembler = new gsThinShellAssembler<3, real_t, true >(mp,dbasis,BCs,force,materialMatrixWrapped);
    else if (!membrane && !wrapped)
      assembler = new gsThinShellAssembler<3, real_t, true >(mp,dbasis,BCs,force,materialMatrix);
    else if (membrane && wrapped)
      assembler = new gsThinShellAssembler<3, real_t, false >(mp,dbasis,BCs,force,materialMatrixWrapped);
    else
      assembler = new gsThinShellAssembler<3, real_t, false >(mp,dbasis,BCs,force,materialMatrix);

    // The assembler uses the material objects as they are (it does not clone them), so the counters of the
    // Split material that is read below are those of the material that is integrated.
    gsMaterialMatrixBase<real_t> * runMaterial = wrapped ? materialMatrixWrapped.get() : materialMatrix.get();
    GISMO_ENSURE(assembler->material(0) == runMaterial,"The assembler does not use the expected material object");

    // Counters of the Split material
    gsMaterialMatrixTFTSplit<3,real_t,true>  * splitLin = (tftsplit && material==0) ? static_cast<gsMaterialMatrixTFTSplit<3,real_t,true> *>(runMaterial) : NULL;
    gsMaterialMatrixTFTSplit<3,real_t,false> * splitNL  = (tftsplit && material!=0) ? static_cast<gsMaterialMatrixTFTSplit<3,real_t,false> *>(runMaterial) : NULL;
    auto resetSplitCounters = [&splitLin,&splitNL]()
    {
      if (splitLin) splitLin->resetCounters();
      if (splitNL)  splitNL->resetCounters();
    };
    auto readSplitCounters = [&splitLin,&splitNL](index_t * c)
    {
      if (splitLin)
      {
        c[0] = splitLin->numThetaFallbacks(); c[1] = splitLin->numTFTFailures();
        c[2] = splitLin->numExactFallbacks(); c[3] = splitLin->numExactReclassified();
      }
      else if (splitNL)
      {
        c[0] = splitNL->numThetaFallbacks(); c[1] = splitNL->numTFTFailures();
        c[2] = splitNL->numExactFallbacks(); c[3] = splitNL->numExactReclassified();
      }
      else
        c[0] = c[1] = c[2] = c[3] = -1;
    };

    // Sample points of the diagnostics: the Gauss points of the assembler
    gsOptionList quOpt;
    quOpt.addReal("quA","",1.0);
    quOpt.addInt("quB","",1);
    const gsMatrix<> sampleNodes = gsQuadrature::getAllNodes(dbasis.basis(0),quOpt);
    WrinkleMonitor monitor(oracle.get(),sampleNodes,tftsplit && midsurface,monitorLevel);

    // Assemble linear system to obtain the force vector
    assembler->assemble();
    gsSparseMatrix<> K = assembler->matrix();
    gsVector<> Fneu = assembler->rhs();
    gsDebugVar(assembler->numDofs());

    // gsSparseMatrix<> diag(assembler->numDofs(),assembler->numDofs());
    // diag.setIdentity();
    // diag *= 0.14;
    // Function for the Jacobian
    gsStructuralAnalysisOps<real_t>::Jacobian_t Jacobian = [&assembler,&mp_def,&monitor](gsVector<real_t> const &x, gsSparseMatrix<real_t> & m)
    {
      // &diag,
      ThinShellAssemblerStatus status;
      assembler->constructSolution(x,mp_def);
      status = assembler->assembleMatrix(mp_def);
      monitor.observeIterate(mp_def);
      m = assembler->matrix();
      // m += diag;
      return status == ThinShellAssemblerStatus::Success;
    };
    // Function for the Residual
    gsStructuralAnalysisOps<real_t>::Residual_t Residual = [&assembler,&mp_def](gsVector<real_t> const &x, gsVector<real_t> & result)
    {
      ThinShellAssemblerStatus status;
      assembler->constructSolution(x,mp_def);
      status = assembler->assembleVector(mp_def);
      result = assembler->rhs();
      return status == ThinShellAssemblerStatus::Success;
    };

    gsVector<> solVector(assembler->numDofs());
    solVector.setZero();

    /*
     * Load stepping:
     * - First 3-5 steps: K+0.14*I
     * - First 20 steps: tension on free boundaries of 4e3 N/m decreasing to 0 N/m
     * - First 20 steps: follower pressure increasing from 0 Pa to 5000 Pa
     */
    real_t load_fac = 0;
    real_t dload_fac= 1./loadSteps;
    real_t dload_fac0 = dload_fac;
    index_t step=0;
    index_t bisected = 0;

    gsMatrix<> writePoints(2,3);
    writePoints.col(0)<<1,0;
    writePoints.col(1)<<0,0;
    writePoints.col(2)<<0,1;

    gsParaviewCollection<> TFcollection(dirname + "/" + "tensionfield");
    gsParaviewCollection<> WDcollection(dirname + "/" + "wrinkle_dir");
    gsParaviewCollection<> GHcollection(dirname + "/" + "wrinkle_gamma_hat");
    gsParaviewCollection<> solcollection(dirname + "/" + "solution");

    // Parametric points of the displacement columns of steps.csv, in column order (see the steps.csv comment below)
    gsMatrix<> qoiPoints(2,6);
    qoiPoints.col(0)<<0,0;
    qoiPoints.col(1)<<1,0;
    qoiPoints.col(2)<<0,1;
    qoiPoints.col(3)<<1,1;
    qoiPoints.col(4)<<0.5,0;
    qoiPoints.col(5)<<1,0.5;

    gsStructuralAnalysisOutput<real_t> writer(dirname + "/" + wn,writePoints);
    std::vector<std::string> pointheaders = {"x","y","z"};
    std::vector<std::string> otherheaders = {"load"};
    writer.init(pointheaders,otherheaders);

    gsVector<> F, M;
    gsStaticDR<real_t> DRM(M,F,Residual);
    gsOptionList DROptions = DRM.options();
    DROptions.setReal("damping",damping);
    DROptions.setReal("alpha",alpha);
    DROptions.setInt("maxIt",1e6);
    DROptions.setReal("tol",tolDR);
    DROptions.setReal("tolE",tolE);
    DROptions.setInt("verbose",verbose);
    DROptions.setInt("ResetIt",(index_t)(100));
    DRM.setOptions(DROptions);

    gsStaticNewton<real_t> NWT(K,F,Jacobian,Residual);
    gsOptionList NWTOptions = NWT.options();
    NWTOptions.setInt("maxIt",maxitNR);
    NWTOptions.setReal("tol",tol);
    NWTOptions.setString("Solver",linSolver);
    NWTOptions.setSwitch("LineSearch",lineSearch);
    NWTOptions.setInt("verbose",verbose);
    NWT.setOptions(NWTOptions);

    std::vector<gsStaticBase<real_t> *> solvers;
    if (DR)
      solvers.push_back(&DRM);
    if (NR)
      solvers.push_back(&NWT);
    gsStaticComposite<real_t> solver(solvers);
    solver.options().setInt("verbose",verbose);

    // Solver chain of the load steps after the first one when only the first load step uses DR
    std::vector<gsStaticBase<real_t> *> solversNROnly;
    solversNROnly.push_back(&NWT);
    gsStaticComposite<real_t> solverNROnly(solversNROnly);
    solverNROnly.options().setInt("verbose",verbose);

    // Diagnostics file: one row per attempted load step. Columns: see WrinkleMonitor::header, followed by
    //   pressure: applied pressure (Pa)
    //   uz_M:     z-displacement at the centre M of the full square, parametric point (0,1) (m)
    //   ux_A:     x-displacement at the free corner A, parametric point (1,0) (m)
    //   ux_c00,uy_c00,uz_c00: displacement at parametric point (0,0), where the west symmetry line meets the south free edge (m)
    //   ux_c10,uy_c10,uz_c10: displacement at parametric point (1,0), the free corner A (m)
    //   ux_c01,uy_c01,uz_c01: displacement at parametric point (0,1), the centre M (m)
    //   ux_c11,uy_c11,uz_c11: displacement at parametric point (1,1), where the north symmetry line meets the east free edge (m)
    //   ux_es,uy_es,uz_es:    displacement at parametric point (0.5,0), the midpoint of the south free edge (m)
    //   ux_ee,uy_ee,uz_ee:    displacement at parametric point (1,0.5), the midpoint of the east free edge (m)
    //   S1_M,S2_M: principal 2nd Piola-Kirchhoff stresses at the mid-surface at M, parametric point (0,1), S1 >= S2 (Pa)
    // The QoI columns are NaN for non-converged attempts.
    std::ofstream stepsFile(dirname + "/steps.csv");
    stepsFile<<std::setprecision(12)<<WrinkleMonitor::header()<<",pressure,uz_M,ux_A"
              <<",ux_c00,uy_c00,uz_c00,ux_c10,uy_c10,uz_c10,ux_c01,uy_c01,uz_c01,ux_c11,uy_c11,uz_c11"
              <<",ux_es,uy_es,uz_es,ux_ee,uy_ee,uz_ee,S1_M,S2_M"<<std::endl;
    index_t attempt = 0, failures = 0, singular = 0;
    real_t lastConverged = 0; // load factor of the last converged step
    gsStopwatch stopwatch;

    while (load_fac <= 1)
    {
      /*
       * Load stepping:
       * - First 3-5 steps: K+0.14*I
       * - First 20 steps: tension on free boundaries of 4e3 N/m decreasing to 0 N/m
       * - First 20 steps: follower pressure increasing from 0 Pa to 5000 Pa
       */
      if (true)
      {
        neuDataX<<maxLoad*(1-load_fac),0,0;
        neuDataY<<0,-maxLoad*(1-load_fac),0;
        neuX.setValue(neuDataX,3);
        neuY.setValue(neuDataY,3);
      }
      else
      {
        neuDataX<<0,0,0;
        neuDataY<<0,0,0;
        neuX.setValue(neuDataX,3);
        neuY.setValue(neuDataY,3);
      }


      assembler->updateBCs(BCs);
      pressFun.setValue(pressure*load_fac,3);
      assembler->setPressure(pressFun);

      gsInfo<<"Load step "<<step<<"; p = "<<pressure*load_fac<<"; load factor = "<<load_fac<<"; load factor step = "<<dload_fac<<":\n";
      assembler->assemble();
      K = assembler->matrix();
      // K += diag;
      F = assembler->rhs();
      assembler->assembleMass(true);
      M = assembler->rhs();

      // With --drfirst only the first load step (and its retries) uses DR
      const bool useDR = DR && !(drfirst && step>0);
      gsStaticComposite<real_t> & activeSolver = (drfirst && step>0) ? solverNROnly : solver;

      activeSolver.initialize();
      if (warmstart)
        activeSolver.setDisplacement(solVector);

      monitor.beginSolve();
      resetSplitCounters();
      stopwatch.restart();
      // gsStaticDR::solve does not catch its error codes (1: not converged, 2: assembly error, 3: solver error), so a failed DR
      // solve surfaces here as an int; it is recorded as a failed attempt like any other
      bool drThrew = false;
      int drStatus = -1;
      try
      {
        activeSolver.solve();
      }
      catch (int errorCode)
      {
        drThrew = true;
        drStatus = static_cast<int>(errorCode==1 ? gsStatus::NotConverged :
                                    errorCode==2 ? gsStatus::AssemblyError :
                                    errorCode==3 ? gsStatus::SolverError : gsStatus::OtherError);
      }
      const real_t solveTime = stopwatch.stop();
      gsInfo<<"Solve time: "<<solveTime<<" s\n";

      // Newton runs only if DR (when used) succeeded
      const bool ranNR = NR && !drThrew && (!useDR || DRM.converged());

      StepInfo info;
      info.step = step;
      info.attempt = attempt++;
      info.load = load_fac;
      info.converged = !drThrew && activeSolver.converged();
      info.statusDR = useDR ? (drThrew ? drStatus : static_cast<int>(DRM.status())) : -1;
      info.statusNR = ranNR ? static_cast<int>(NWT.status()) : -1;
      info.itDR = DR ? (useDR ? DRM.iterations() : 0) : -1;
      info.itNR = NR ? (ranNR ? NWT.iterations() : 0) : -1;
      info.usesNR = ranNR;
      info.time = solveTime;
      readSplitCounters(info.splitCounters);
      if (!info.converged)
        failures++;
      if (info.statusNR==static_cast<int>(gsStatus::SolverError) || info.statusDR==static_cast<int>(gsStatus::SolverError))
        singular++;
      info.bisections = bisected;
      info.failures = failures;
      info.singular = singular;

      if (info.converged)
      {
        solVector = activeSolver.solution();
        mp_def = assembler->constructSolution(solVector);
      }
      const PointStats pointStats = monitor.evaluate(mp_def,runMaterial,model);
      if (info.converged)
      {
        monitor.accept();
        lastConverged = load_fac;
      }

      if (!info.converged)
      {
        monitor.writeColumns(stepsFile,info,pointStats);
        stepsFile<<","<<pressure*load_fac<<","<<std::numeric_limits<real_t>::quiet_NaN()<<","<<std::numeric_limits<real_t>::quiet_NaN();
        for (index_t k = 0; k != 20; ++k)
          stepsFile<<","<<std::numeric_limits<real_t>::quiet_NaN();
        stepsFile<<std::endl;

        gsDebug<<"Number of failures: "<<bisected<<"\n";
        if (bisected > 20) // failed already 20 times
        {
          gsWarn<<"Simulation terminated because of too many failures...\n";
          break;
        }
        gsWarn<<"Load step "<<step<<" did not converge\n";
        GISMO_ASSERT(load_fac!=0,"load_fac is zero but no convergence on the first step. Try to increase the number of iterations");
        if (load_fac==0)
        {
          gsWarn<<"Simulation terminated because the first load step did not converge\n";
          break;
        }
        load_fac -= dload_fac;
        if (math::abs(load_fac-1) < 1e-3)
        {
          gsInfo<<"Simulation terminated because the load factor is close to 1\n";
          break;
        }
        dload_fac /= 2;
        load_fac += dload_fac;
        bisected++;
        continue;
      }

      // if (material==0)
      //   dynamic_cast<gsMaterialMatrixTFT<3,real_t,true> *>(materialMatrixTFT)->updateDeformed(&mp_def);
      // else if (material==1 || material==2)
      //   dynamic_cast<gsMaterialMatrixTFT<3,real_t,false> *>(materialMatrixTFT)->updateDeformed(&mp_def);
      // else GISMO_ERROR("Material unknown");

      gsMultiPatch<> deformation = mp_def;
      for (size_t k = 0; k != mp_def.nPatches(); ++k)
          deformation.patch(k).coefs() -= mp.patch(k).coefs();

      {
        const gsMatrix<> qoi = deformation.patch(0).eval(writePoints);
        monitor.writeColumns(stepsFile,info,pointStats);
        stepsFile<<","<<pressure*load_fac<<","<<qoi(2,2)<<","<<qoi(0,0);
        const gsMatrix<> qoiDisp = deformation.patch(0).eval(qoiPoints);
        for (index_t k = 0; k != qoiPoints.cols(); ++k)
          stepsFile<<","<<qoiDisp(0,k)<<","<<qoiDisp(1,k)<<","<<qoiDisp(2,k);

        // Principal stresses at M: one point per call, the order of the pair depends on the material
        const gsMatrix<> Z0 = gsMatrix<>::Zero(1,1);
        gsMatrix<> ptM(2,1);
        ptM<<0,1;
        const gsMatrix<> SM = evalMaterial<MaterialOutput::PStress>(runMaterial,mp_def,ptM,Z0);
        stepsFile<<","<<std::max(SM(0,0),SM(1,0))<<","<<std::min(SM(0,0),SM(1,0))<<std::endl;
      }

      // ! [Export visualization in ParaView]
      if (plot)
      {
          std::string fileName;
          gsField<> solField(mp_def, deformation);
          gsInfo<<"Plotting in Paraview...\n";

          fileName = dirname + "/" + "solution" + util::to_string(step);
          gsWriteParaview<>( solField, fileName, 1000, true);
          solcollection.addPart(gsFileManager::getFilename(fileName) + "0.vts",step);

          // ev.options().setSwitch("plot.elements", true);
          // ev.writeParaview( u_sol   , G, "solution");

          // gsFileManager::open("solution.pvd");


          assembler->constructSolution(solVector,mp_def);
          gsPiecewiseFunction<> TFes;
          assembler->constructStress(mp_def,TFes,stress_type::tension_field);
          gsField<> TF(mp_def,TFes, true);

          fileName = dirname + "/" + "tensionfield" + util::to_string(step);
          gsWriteParaview(TF, fileName,5000);
          TFcollection.addPart(gsFileManager::getFilename(fileName) + "0.vts",step);

          if (wrapped)
          {
            gsPiecewiseFunction<> WDes;
            assembler->constructStress(mp_def,WDes,stress_type::wrinkle_dir);
            gsField<> WD(mp_def,WDes, true);
            fileName = dirname + "/" + "wrinkle_dir" + util::to_string(step);
            gsWriteParaview(WD, fileName,5000);
            WDcollection.addPart(gsFileManager::getFilename(fileName) + "0.vts",step);

            gsPiecewiseFunction<> GHes;
            assembler->constructStress(mp_def,GHes,stress_type::wrinkle_gamma_hat);
            gsField<> GH(mp_def,GHes, true);
            fileName = dirname + "/" + "wrinkle_gamma_hat" + util::to_string(step);
            gsWriteParaview(GH, fileName,5000);
            GHcollection.addPart(gsFileManager::getFilename(fileName) + "0.vts",step);
          }


          gsInfo <<"Maximum deformation coef: "
                 << deformation.patch(0).coefs().colwise().maxCoeff() <<".\n";
          gsInfo <<"Minimum deformation coef: "
                 << deformation.patch(0).coefs().colwise().minCoeff() <<".\n";
      }
      if (write)
      {
        gsMatrix<> result = deformation.patch(0).eval(writePoints);
        gsVector<> data(1);
        data<<load_fac;
        writer.add(result,data);

        gsWrite(mp_def,dirname + "/" + "deformed");
      }

/*
      if (bisected!=0)
      {
        bisected = 0;
        dload_fac *= 2;
      }
      else
        dload_fac = dload_fac0;
*/
      dload_fac = dload_fac0;
      step++;
      if (1 - load_fac < dload_fac && 1 - load_fac > 1e-10)
        dload_fac = 1 - load_fac;
      load_fac += dload_fac;
      gsInfo<<"load factor = "<<load_fac<<"; load factor step = "<<dload_fac<<"\n";
    }

    solcollection.save();
    TFcollection.save();
    if (wrapped)
    {
      WDcollection.save();
      GHcollection.save();
    }

    mp_def = assembler->constructSolution(solVector);

    gsMultiPatch<> deformation = mp_def;
    for (size_t k = 0; k != mp_def.nPatches(); ++k)
        deformation.patch(k).coefs() -= mp.patch(k).coefs();


    // ! [Export visualization in ParaView]
    if (stress)
    {
      gsPiecewiseFunction<> membraneStresses;
      assembler->constructStress(mp_def,membraneStresses,stress_type::membrane);
      gsField<> membraneStress(mp_def,membraneStresses, true);

      gsPiecewiseFunction<> flexuralStresses;
      assembler->constructStress(mp_def,flexuralStresses,stress_type::flexural);
      gsField<> flexuralStress(mp_def,flexuralStresses, true);

      gsPiecewiseFunction<> stretches;
      assembler->constructStress(mp_def,stretches,stress_type::principal_stretch);
      gsField<> Stretches(mp_def,stretches, true);

      gsPiecewiseFunction<> pstrain_m;
      assembler->constructStress(mp_def,pstrain_m,stress_type::principal_membrane_strain);
      gsField<> pstrainM(mp_def,pstrain_m, true);

      gsPiecewiseFunction<> pstrain_f;
      assembler->constructStress(mp_def,pstrain_f,stress_type::principal_flexural_strain);
      gsField<> pstrainF(mp_def,pstrain_f, true);

      gsPiecewiseFunction<> pstress_m;
      assembler->constructStress(mp_def,pstress_m,stress_type::principal_stress_membrane);
      gsField<> pstressM(mp_def,pstress_m, true);

      gsPiecewiseFunction<> pstress_f;
      assembler->constructStress(mp_def,pstress_f,stress_type::principal_stress_flexural);
      gsField<> pstressF(mp_def,pstress_f, true);

      gsPiecewiseFunction<> stretch1;
      assembler->constructStress(mp_def,stretch1,stress_type::principal_stretch_dir1);
      gsField<> stretchDir1(mp_def,stretch1, true);

      gsPiecewiseFunction<> stretch2;
      assembler->constructStress(mp_def,stretch2,stress_type::principal_stretch_dir2);
      gsField<> stretchDir2(mp_def,stretch2, true);

      gsPiecewiseFunction<> stretch3;
      assembler->constructStress(mp_def,stretch3,stress_type::principal_stretch_dir3);
      gsField<> stretchDir3(mp_def,stretch3, true);

      gsPiecewiseFunction<> VMStresses;
      assembler->constructStress(mp_def,VMStresses,stress_type::von_mises_membrane);
      gsField<> VMStress(mp_def,VMStresses, true);


      gsWriteParaview(membraneStress,dirname + "/" + "MembraneStress",5000);
      gsWriteParaview(VMStress,dirname + "/" + "MembraneStressVM",5000);
      gsWriteParaview(Stretches,dirname + "/" + "PrincipalStretch",5000);
      gsWriteParaview(pstrainM,dirname + "/" + "PrincipalMembraneStrain",5000);
      gsWriteParaview(pstrainF,dirname + "/" + "PrincipalFlexuralStrain",5000);
      gsWriteParaview(pstressM,dirname + "/" + "PrincipalMembraneStress",5000);
      gsWriteParaview(pstressF,dirname + "/" + "PrincipalFlexuralStress",5000);
      gsWriteParaview(stretchDir1,dirname + "/" + "PrincipalDirection1",5000);
      gsWriteParaview(stretchDir1,dirname + "/" + "PrincipalDirection1",5000);
      gsWriteParaview(stretchDir2,dirname + "/" + "PrincipalDirection2",5000);
      gsWriteParaview(stretchDir3,dirname + "/" + "PrincipalDirection3",5000);

    }

    delete assembler;

    const bool reached = lastConverged >= 1 - 1e-12;
    if (!reached)
      gsInfo<<"FINAL LOAD NOT REACHED: last converged "<<lastConverged<<" of 1\n";
    return reached ? EXIT_SUCCESS : EXIT_FAILURE;
}
#else//gsKLShell_ENABLED
int main(int argc, char *argv[])
{
    gsWarn<<"G+Smo is not compiled with the gsKLShell module.";
    return EXIT_FAILURE;
}
#endif