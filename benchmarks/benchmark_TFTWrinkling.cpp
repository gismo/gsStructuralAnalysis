/** @file benchmark_TFTWrinkling.cpp

    @brief Benchmarks for membrane wrinkling with the tension-field and split formulations under
           displacement control (Newton)

    Test cases (\a -t):
    - 0: shear of a rectangular SvK membrane, plane (d = 2) membrane, 380 mm x 128 mm,
         t = 0.025 mm, E = 3500 MPa, nu = 0.31 (Zhang & Kiendl (2025), Sec. 5.2). The bottom edge is
         fixed. Phase 1 displaces the top edge vertically by 0.05 mm (\a --presteps steps), phase 2
         displaces it horizontally by \a --ux (\a -N steps) with the vertical displacement held.
         Whenever \a --ux >= 1.6, the schedule contains u_x = 1.6 mm exactly as a load step. Reaction
         force on the top edge and the principal stresses along the midline y = H/2 are reported.
    - 1: stretched incompressible Mooney-Rivlin sheet, quarter model, d = 2, L x W = 280 mm x 140 mm,
         t = 0.14 mm, SI units (Zhang & Kiendl (2026), Sec. 5.2). The grip (right edge) is displaced
         uniformly to the strain \a --eps. Reaction force on the grip and the width of the wrinkled
         zone along the vertical centre line are reported.
    - 2: twist of a trapezoidal incompressible Neo-Hookean membrane, d = 3 (Zhang & Kiendl (2026),
         Sec. 5.4). Top width 100 mm, bottom width 50 mm, height 100 mm, t = 5e-3 mm, E = 1e3 MPa. The
         top edge is clamped. Phase 1 displaces the bottom edge downwards by 50 mm (\a --presteps
         steps), phase 2 rotates it about the vertical central axis up to \a --theta (times pi,
         \a -N steps). The torque on the rotated edge is reported.

    Wrinkling model (\a --model): \a none, \a tft (gsMaterialMatrixTFT) or \a split
    (gsMaterialMatrixTFTSplit; options \a --eta, \a --exact, \a --midsurface, \a --s2crit,
    \a --slackmult).

    Every attempted load step appends one row to \a <out>/steps.csv (see WrinkleMonitor and the
    column list next to the header of the case). A step that does not converge is halved, at most
    \a --maxbisect times in a row; the run then stops and the exit code is 1. With \a --plot, the
    displacement and the wrinkling fields are written as ParaView collections.

    With \a --newtondiag (needs \a --monitor 2), the files \a newton_iter.csv, \a newton_flips.csv
    and \a converged_tangent.csv are written next to steps.csv (see NewtonDiag for the columns):
    per Newton iterate the residual, the update and the state flips of the sample points, and, for
    the attempts of the last phase with load >= \a --diagfrom, the spectrum of the symmetric part and
    the singular values of the tangent matrix. The diagnostics only observe the iteration. Without
    the flag no diagnostic file is written.

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
#include <gsKLShell/src/gsMaterialMatrixIntegrate.h>
#endif

#include <gsStructuralAnalysis/src/gsStaticSolvers/gsStaticNewton.h>

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

    /// State of the sample points at the iterate observed last in the current solve (convention of the
    /// TensionField output: 1 = taut, 0 = wrinkled, -1 = slack); nullptr before the first observation
    const gsMatrix<> * iterateState() const { return m_hasPrev ? &m_prev : nullptr; }

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
            // With Eta > 0 this is the blended stress that the assembler integrates.
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
    /// - load:            load coordinate of the active phase (case-specific, see the case columns)
    /// - converged:       1 if the solver chain converged, else 0
    /// - status_nr, status_dr: static_cast<int>(gsStatus) of the Newton / DR solver (0 = Success, 1 = NotConverged,
    ///                    3 = SolverError); -1 if the solver did not run in this attempt
    /// - it_nr, it_dr:    iterations of the Newton / DR solver in this attempt; -1 if the solver is not enabled,
    ///                    0 if it is enabled but did not run
    /// - bisections:      cumulative number of load step halvings
    /// - failures:        cumulative number of non-converged attempts
    /// - singular:        cumulative number of attempts that ended with gsStatus::SolverError (singular tangent)
    /// - jac_calls:       number of tangent matrix assemblies in this attempt (Newton iterations + 1); -1 if not observed
    /// - chatter:         total number of state flips of the sample points between consecutive tangent
    ///                    assemblies of this attempt, summed over all consecutive pairs (a point that flips k
    ///                    times counts k times, so the value can exceed n_points); -1 if not observed
    /// - n_points:        number of sample points (Gauss points times fibres); -1 without point statistics
    /// - n_taut, n_wrinkled, n_slack: number of points in each state of the oracle
    /// - state_changes:   number of points whose oracle state differs from that of the last accepted step
    /// - gh_max, gh_p95:  maximum and nearest-rank 95th percentile of the wrinkle strain gamma_hat
    ///                    (= gamma G^{ab} n_a n_b, dimensionless) of the run material over wrinkled points;
    ///                    NaN without wrinkling model or wrinkled points
    /// - n4S_nonpos:      number of wrinkled points where the stress along the wrinkles, n4.S, of the run
    ///                    material is <= 0 (with Eta > 0: of the blended stress that the assembler integrates);
    ///                    -1 without wrinkling model
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

/// Spectral quantities of one tangent matrix K = -dR/dx (R = the assembler's right-hand side).
/// Count fields are stored as reals so that "not computed" can be NaN in the CSV.
struct TangentStats
{
    TangentStats()
    {
        const real_t nan = std::numeric_limits<real_t>::quiet_NaN();
        asymRel = lamMin = lam2 = lam3 = lamMax = nNeg = nSmall = sigMin = sigMax = time = nan;
    }

    real_t asymRel;                 ///< ||K - K^T||_F / ||K||_F
    real_t lamMin, lam2, lam3, lamMax; ///< smallest three and largest eigenvalue of (K + K^T)/2
    real_t nNeg, nSmall;            ///< #{lambda < -1e-12 max|lambda|}, #{|lambda| <= 1e-12 max|lambda|}
    real_t sigMin, sigMax;          ///< extreme singular values of K
    real_t time;                    ///< wall time of the analysis (s)
};

/**
 * @brief Dense spectral analysis of the tangent \a m.
 *
 * The tangent of a tension-field material is not exactly symmetric where the state changes, so the
 * eigenvalues are taken of the symmetric part K_s = (K + K^T)/2 (the quadratic form x^T K x depends
 * on K_s only; a negative eigenvalue of K_s signals loss of positive definiteness). The singular
 * values of K itself are computed with a bidiagonal divide and conquer SVD, not from the eigenvalues
 * of K^T K, which would square the condition number.
 *
 * Complexity: O(n^3), n = number of free degrees of freedom.
 */
TangentStats analyseTangent(const gsSparseMatrix<real_t> & m)
{
    gsStopwatch sw;
    TangentStats s;
    const gsMatrix<> Kd = m.toDense();
    const gsMatrix<> Kt = Kd.transpose();
    const real_t normK = Kd.norm();
    s.asymRel = (Kd - Kt).norm() / normK;

    const gsMatrix<> Ks = 0.5*(Kd + Kt);
    gsMatrix<>::SelfAdjEigenSolver es(Ks,gsEigen::EigenvaluesOnly);
    const gsMatrix<> lam = es.eigenvalues();   // ascending
    const index_t n = lam.size();
    if (n > 0)
    {
        s.lamMin = lam(0);
        s.lam2   = lam(std::min<index_t>(1,n-1));
        s.lam3   = lam(std::min<index_t>(2,n-1));
        s.lamMax = lam(n-1);
        const real_t tol = 1e-12 * std::max(std::abs(lam(0)),std::abs(lam(n-1)));
        s.nNeg = s.nSmall = 0;
        for (index_t i = 0; i!=n; ++i)
        {
            if (lam(i) < -tol)
                s.nNeg += 1;
            if (std::abs(lam(i)) <= tol)
                s.nSmall += 1;
        }

        gsEigen::BDCSVD<gsMatrix<real_t>::Base> svd(Kd);
        const gsMatrix<> sv = svd.singularValues();   // descending
        s.sigMax = sv(0);
        s.sigMin = sv(sv.size()-1);
    }
    s.time = sw.stop();
    return s;
}

/**
 * @brief Per-Newton-iterate diagnostics of a load-step attempt: residual and update history, state
 *        flips of the sample points and, inside a window of the load, the spectrum of the tangent.
 *
 * The object only observes. It reads the iterate x, the residual vector, the tangent matrix and the
 * state of the WrinkleMonitor, and evaluates no material, so that the iteration is unchanged.
 *
 * Iterate k of an attempt is x_k = U0 + DeltaU_k (x_0 = U0). The residual at x_k is observed by the
 * residual callback, the tangent at x_k (all but the last iterate) by the Jacobian callback. A row is
 * finished when the next iterate arrives or the attempt ends.
 *
 * Files (header line first, std::setprecision(12), every line flushed; -1 or NaN = not available):
 * - newton_iter.csv, one row per iterate of every attempt:
 *   attempt,step,phase,load,in_window,k,res_norm,res_rel,du_rel,flips,flips_t2w,flips_w2t,flips_w2s,
 *   flips_s2w,flips_t2s,flips_s2t,pingpong,n_taut,n_wrinkled,n_slack, then the tangent columns
 *   asym_rel,lam_sym_min,lam_sym_2,lam_sym_3,lam_sym_max,n_neg_sym,n_small_sym,sig_min,sig_max,
 *   diag_time_s (see TangentStats; only inside the window and when the iterate has a tangent).
 *   - res_norm = ||R(x_k)||_2, res_rel = res_norm / ||R(x_0)||_2 (the norm 1 if ||R(x_0)|| = 0, as in
 *     gsStaticNewton), du_rel = ||x_k - x_{k-1}|| / ||x_k - x_0|| (NaN for k = 0).
 *   - flips = number of sample points whose state at x_k differs from that at x_{k-1}, split by
 *     transition (t = taut, w = wrinkled, s = slack); -1 for k = 0 and for the last iterate. The
 *     sum over an attempt equals the chatter of steps.csv.
 *   - pingpong = flipped points whose new state equals their state at x_{k-2}.
 *   - n_taut, n_wrinkled, n_slack = states at x_k (-1 for the last iterate).
 *   - k = -1 marks a tangent that was not evaluated at the open iterate (never expected).
 * - newton_flips.csv (attempts in the window): attempt,k,point,fibre,u,v,old,new, one line per flipped
 *   point; point = column of the state matrix, fibre = point / (nodes per fibre), (u,v) = its
 *   parametric position, old/new = state at x_{k-1} / x_k.
 * - converged_tangent.csv, one row per accepted step: step,phase,load,n_taut,n_wrinkled,n_slack and
 *   the tangent columns, with the tangent assembled on the converged configuration.
 *
 * Complexity: O(n^3) per in-window Jacobian call for the dense eigen/SVD (n = free dofs; 840 for the
 * shear case); O(number of sample points) per iterate otherwise.
 */
class NewtonDiag
{
public:
    /// \param outDir      directory of the output files
    /// \param nodes       parametric sample points (dim x n) of one fibre
    /// \param diagFrom    attempts of the last phase with load >= diagFrom are analysed in full
    NewtonDiag(const std::string & outDir, const gsMatrix<> & nodes, real_t diagFrom)
    :
    m_nodes(nodes), m_diagFrom(diagFrom), m_active(false), m_rowOpen(false), m_inWindow(false),
    m_havePrev(false), m_havePrev2(false), m_attempt(0), m_step(0), m_phase(0), m_nextK(0), m_load(0), m_resIni(1)
    {
        m_iter.open(outDir + "/newton_iter.csv");
        m_flips.open(outDir + "/newton_flips.csv");
        m_conv.open(outDir + "/converged_tangent.csv");
        if (!m_iter || !m_flips || !m_conv)
            GISMO_ERROR("Could not open the diagnostic files in "<<outDir);
        m_iter<<std::setprecision(12)
              <<"attempt,step,phase,load,in_window,k,res_norm,res_rel,du_rel,flips,flips_t2w,flips_w2t,"
                "flips_w2s,flips_s2w,flips_t2s,flips_s2t,pingpong,n_taut,n_wrinkled,n_slack,"
                "asym_rel,lam_sym_min,lam_sym_2,lam_sym_3,lam_sym_max,n_neg_sym,n_small_sym,sig_min,sig_max,"
                "diag_time_s"<<std::endl;
        m_flips<<std::setprecision(12)<<"attempt,k,point,fibre,u,v,old,new"<<std::endl;
        m_conv<<std::setprecision(12)
              <<"step,phase,load,n_taut,n_wrinkled,n_slack,"
                "asym_rel,lam_sym_min,lam_sym_2,lam_sym_3,lam_sym_max,n_neg_sym,n_small_sym,sig_min,sig_max,"
                "diag_time_s"<<std::endl;
    }

    /// Whether an attempt at load \a D of phase index \a ph (of \a nPhases) is analysed in full
    bool inWindow(size_t ph, size_t nPhases, real_t D) const
    {
        return ph + 1 == nPhases && D >= m_diagFrom - 1e-12;
    }

    /// Starts an attempt; call before the solve
    void beginAttempt(index_t attempt, index_t step, index_t phase, real_t load, bool inWindow)
    {
        m_attempt = attempt; m_step = step; m_phase = phase; m_load = load; m_inWindow = inWindow;
        m_active = true;
        m_rowOpen = false;
        m_nextK = 0;
        m_havePrev = m_havePrev2 = false;
    }

    /// Residual callback: observes R(x). The first call of an attempt is the iterate x_0.
    void observeResidual(const gsVector<real_t> & x, const gsVector<real_t> & result)
    {
        if (!m_active)
            return;
        const real_t nan = std::numeric_limits<real_t>::quiet_NaN();
        if (m_rowOpen)
            writeRow(m_row);
        m_row = Row();
        m_row.k = m_nextK++;
        m_row.resNorm = result.norm();
        if (m_row.k == 0)
        {
            m_x0 = x;
            m_resIni = (m_row.resNorm == 0) ? 1.0 : m_row.resNorm;
            m_row.duRel = nan;
        }
        else
        {
            const real_t DU = (x - m_x0).norm();
            m_row.duRel = (DU == 0) ? nan : (x - m_xk).norm() / DU;
        }
        m_row.resRel = m_row.resNorm / m_resIni;
        m_xk = x;
        m_rowOpen = true;
    }

    /// Jacobian callback: observes the tangent \a m at x. \a state is the oracle state at x
    /// (WrinkleMonitor::iterateState after the observation of this iterate).
    void observeJacobian(const gsVector<real_t> & x, const gsSparseMatrix<real_t> & m, const gsMatrix<> * state)
    {
        if (!m_active)
            return;
        if (!m_rowOpen || (x - m_xk).norm() != 0)
        {
            Row bad;
            bad.k = -1;
            writeRow(bad);
            return;
        }

        if (state != nullptr)
        {
            m_row.nTaut = m_row.nWrinkled = m_row.nSlack = 0;
            for (index_t i = 0; i!=state->cols(); ++i)
            {
                const real_t s = (*state)(0,i);
                if      (s==1)  ++m_row.nTaut;
                else if (s==0)  ++m_row.nWrinkled;
                else if (s==-1) ++m_row.nSlack;
            }

            if (m_havePrev && m_prev.cols()==state->cols())
            {
                m_row.flips = m_row.pingpong = 0;
                for (int j = 0; j!=6; ++j)
                    m_row.flipsByType[j] = 0;
                const index_t nNodes = m_nodes.cols();
                for (index_t i = 0; i!=state->cols(); ++i)
                {
                    const real_t sNew = (*state)(0,i), sOld = m_prev(0,i);
                    if (sNew == sOld)
                        continue;
                    ++m_row.flips;
                    const int type = flipType(sOld,sNew);
                    if (type >= 0)
                        ++m_row.flipsByType[type];
                    if (m_havePrev2 && m_prev2.cols()==state->cols() && m_prev2(0,i)==sNew)
                        ++m_row.pingpong;
                    if (m_inWindow)
                    {
                        const index_t node = i % nNodes;
                        m_flips<<m_attempt<<","<<m_row.k<<","<<i<<","<<i/nNodes<<","
                               <<m_nodes(0,node)<<","<<m_nodes(1,node)<<","
                               <<static_cast<int>(sOld)<<","<<static_cast<int>(sNew)<<std::endl;
                    }
                }
            }
            m_prev2 = m_prev; m_havePrev2 = m_havePrev;
            m_prev = *state; m_havePrev = true;
        }

        if (m_inWindow)
        {
            m_row.tangent = analyseTangent(m);
            m_row.hasTangent = true;
        }
    }

    /// Finishes the attempt; call right after the solve
    void endAttempt()
    {
        if (m_rowOpen)
            writeRow(m_row);
        m_rowOpen = false;
        m_active = false;
    }

    /// Writes one row of converged_tangent.csv for the tangent \a m on the converged configuration
    void writeConverged(index_t step, index_t phase, real_t load, const PointStats & st, const gsSparseMatrix<real_t> & m)
    {
        const TangentStats t = analyseTangent(m);
        m_conv<<step<<","<<phase<<","<<load<<","<<st.nTaut<<","<<st.nWrinkled<<","<<st.nSlack<<",";
        writeTangent(m_conv,t);
        m_conv<<std::endl;
    }

private:
    struct Row
    {
        Row() : k(-1), resNorm(std::numeric_limits<real_t>::quiet_NaN()), resRel(resNorm), duRel(resNorm),
                flips(-1), pingpong(-1), nTaut(-1), nWrinkled(-1), nSlack(-1), hasTangent(false)
        {
            for (int i = 0; i!=6; ++i)
                flipsByType[i] = -1;
        }

        index_t k;
        real_t resNorm, resRel, duRel;
        index_t flips, flipsByType[6], pingpong, nTaut, nWrinkled, nSlack;
        bool hasTangent;
        TangentStats tangent;
    };

    /// Index of the transition old -> new in the order t2w, w2t, w2s, s2w, t2s, s2t (-1: other values)
    static int flipType(real_t o, real_t n)
    {
        if (o==1  && n==0)  return 0;
        if (o==0  && n==1)  return 1;
        if (o==0  && n==-1) return 2;
        if (o==-1 && n==0)  return 3;
        if (o==1  && n==-1) return 4;
        if (o==-1 && n==1)  return 5;
        return -1;
    }

    static void writeTangent(std::ostream & os, const TangentStats & t)
    {
        os<<t.asymRel<<","<<t.lamMin<<","<<t.lam2<<","<<t.lam3<<","<<t.lamMax<<","<<t.nNeg<<","<<t.nSmall<<","
          <<t.sigMin<<","<<t.sigMax<<","<<t.time;
    }

    void writeRow(const Row & r)
    {
        m_iter<<m_attempt<<","<<m_step<<","<<m_phase<<","<<m_load<<","<<(m_inWindow ? 1 : 0)<<","<<r.k<<","
              <<r.resNorm<<","<<r.resRel<<","<<r.duRel<<","<<r.flips;
        for (int i = 0; i!=6; ++i)
            m_iter<<","<<r.flipsByType[i];
        m_iter<<","<<r.pingpong<<","<<r.nTaut<<","<<r.nWrinkled<<","<<r.nSlack<<",";
        writeTangent(m_iter,r.hasTangent ? r.tangent : TangentStats());
        m_iter<<std::endl;
    }

    std::ofstream m_iter, m_flips, m_conv;
    gsMatrix<> m_nodes;
    real_t m_diagFrom;
    bool m_active, m_rowOpen, m_inWindow, m_havePrev, m_havePrev2;
    index_t m_attempt, m_step, m_phase, m_nextK;
    real_t m_load, m_resIni;
    Row m_row;
    gsVector<real_t> m_x0, m_xk;
    gsMatrix<> m_prev, m_prev2;
};

/// Builds the wrinkling-state oracle around \a base (not owned): a Split decorator with Eta = 0, no exact
/// solve and S2crit = 1e300, so that the TensionField output is the criterion state
template <short_t d, bool linear>
gsMaterialMatrixBase<real_t>::uPtr makeOracle(MaterialBase * base, real_t slackMultiplier)
{
    gsMaterialMatrixBase<real_t>::uPtr oracle =
        memory::make_unique(new gsMaterialMatrixTFTSplit<d,real_t,linear>(base));
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

/// Material objects of a run and access to the counters of the Split material
struct RunMaterials
{
    gsMaterialMatrixBase<real_t>::uPtr base, wrapped, oracle;
    MaterialBase * run;                       ///< material that the assembler integrates
    int model;                                ///< 0: none, 1: TFT, 2: Split
    std::function<void()> resetCounters;
    std::function<void(index_t *)> readCounters; ///< theta fallbacks, TFT failures, exact fallbacks, exact reclassified
};

/// Creates the base material (\a options select the law), its wrinkling decorator and the oracle
template <short_t d, bool linear>
void buildMaterials(RunMaterials & m, const gsMultiPatch<> & mp, const gsFunctionSet<> & t, const gsFunctionSet<> & rho,
                    std::vector<gsFunctionSet<>*> & parameters, const gsOptionList & options, int model,
                    real_t eta, bool exact, bool midsurface, real_t s2crit, real_t slackMultiplier)
{
    m.base = getMaterialMatrix<d,real_t>(mp,t,parameters,rho,options);
    m.model = model;
    m.oracle = makeOracle<d,linear>(m.base.get(),slackMultiplier);
    m.resetCounters = [](){};
    m.readCounters = [](index_t * c){ c[0] = c[1] = c[2] = c[3] = -1; };
    if (model==0)
    {
        m.run = m.base.get();
    }
    else if (model==1)
    {
        m.wrapped = memory::make_unique(new gsMaterialMatrixTFT<d,real_t,linear>(*m.base));
        m.wrapped->options().setReal("SlackMultiplier",slackMultiplier);
        m.run = m.wrapped.get();
    }
    else
    {
        m.wrapped = memory::make_unique(new gsMaterialMatrixTFTSplit<d,real_t,linear>(*m.base));
        setSplitOptions(m.wrapped.get(),eta,exact,midsurface,s2crit,slackMultiplier);
        m.run = m.wrapped.get();
        gsMaterialMatrixTFTSplit<d,real_t,linear> * split = static_cast<gsMaterialMatrixTFTSplit<d,real_t,linear> *>(m.run);
        m.resetCounters = [split](){ split->resetCounters(); };
        m.readCounters = [split](index_t * c)
        {
            c[0] = split->numThetaFallbacks(); c[1] = split->numTFTFailures();
            c[2] = split->numExactFallbacks(); c[3] = split->numExactReclassified();
        };
    }
}

/// Creates \a path and its missing parents
bool makeDirs(const std::string & path)
{
    bool ok = true;
    for (size_t pos = path.find('/',1); pos != std::string::npos; pos = path.find('/',pos+1))
        ok = gsFileManager::mkdir(path.substr(0,pos)) && ok;
    return gsFileManager::mkdir(path) && ok;
}

/// Principal stresses of \a run on the line y = H/2 of the shear rectangle at \a n points with
/// xi in [0.25,0.75]; the major value is S1, the minor S2.
/// Each point is evaluated in a separate call, because the multi-point evaluation of the principal
/// stresses of the wrapped materials reuses the metric of the last point.
void midlineStress(MaterialBase * run, const gsMultiPatch<> & def, real_t length, index_t n,
                   std::vector<real_t> & x, std::vector<real_t> & S1, std::vector<real_t> & S2)
{
    const gsMatrix<> Z0 = gsMatrix<>::Zero(1,1);
    x.resize(n); S1.resize(n); S2.resize(n);
    gsMatrix<> pt(2,1);
    for (index_t i = 0; i!=n; ++i)
    {
        const real_t xi = 0.25 + 0.5*i/(n-1);
        pt<<xi,0.5;
        const gsMatrix<> S = evalMaterial<MaterialOutput::PStress>(run,def,pt,Z0);
        x[i]  = xi*length;
        S1[i] = std::max(S(0,0),S(1,0));
        S2[i] = std::min(S(0,0),S(1,0));
    }
}

/**
 * @brief Torque on the rotated edge of the trapezoid: the virtual work of the membrane stress resultants
 *        N on the rotation velocity field v = du/dtheta of the edge,
 *        T = int N : sym(F^T grad v) dA0.
 *
 * The velocity field v is the lift of the edge values (zero at all other control points), so that the
 * integral equals the sum of the internal nodal forces of the edge times v (N mm for mm / N units).
 */
real_t virtualWorkTorque(const gsMultiBasis<> & dbasis, const gsMultiPatch<> & mp, const gsMultiPatch<> & def,
                         const gsMaterialMatrixContainer<real_t> & materials,
                         gsFunctionExpr<> & vx, gsFunctionExpr<> & vz)
{
    gsExprAssembler<> A(1,1);
    A.setIntegrationElements(dbasis);
    gsExprEvaluator<> ev(A);
    gsExprAssembler<>::space u = A.getSpace(dbasis,3);

    gsBoundaryConditions<> bcDummy;
    bcDummy.addCondition(0,boundary::south,condition_type::dirichlet,&vx,0,false,0);
    bcDummy.addCondition(0,boundary::south,condition_type::dirichlet,&vz,0,false,2);
    bcDummy.setGeoMap(mp);
    u.setup(bcDummy,dirichlet::interpolation,-1);

    // Velocity field: the Dirichlet values on the edge, zero elsewhere
    gsMultiPatch<> geo;
    for (size_t b = 0; b!=dbasis.nBases(); b++)
    {
        gsMatrix<> coefs(dbasis.basis(b).size(),mp.geoDim());
        coefs.setZero();
        for (short_t c = 0; c!=mp.geoDim(); c++)
            for (index_t k = 0; k!=dbasis.basis(b).size(); k++)
            {
                const index_t ii = u.mapper().index(k,b,c);
                if (u.mapper().is_boundary(k,b,c))
                    coefs(k,c) = u.fixedPart().at(u.mapper().global_to_bindex(ii));
            }
        geo.addPatch(*dbasis.basis(b).makeGeometry(give(coefs)));
    }

    A.initSystem();
    gsExprAssembler<>::geometryMap ori  = A.getMap(mp);
    gsExprAssembler<>::geometryMap deff = A.getMap(def);
    auto v = A.getCoeff(geo);
    gsMaterialMatrixIntegrate<real_t,MaterialOutput::VectorN> S0f(materials,&mp,&def);
    auto S0 = A.getCoeff(S0f);
    auto N   = S0.tr();
    auto dEm = flat( jac(deff).tr() * jac(v) );
    return ev.integral( (N*dEm.tr())*meas(ori) );
}

/// Writes the sampled midline stresses as CSV (x, S1, S2)
void writeMidline(const std::string & fileName, const std::vector<real_t> & x, const std::vector<real_t> & S1, const std::vector<real_t> & S2)
{
    std::ofstream file(fileName);
    file<<std::setprecision(12)<<"x,S1,S2\n";
    for (size_t i = 0; i!=x.size(); ++i)
        file<<x[i]<<","<<S1[i]<<","<<S2[i]<<"\n";
}

/// A phase of the load schedule: the load coordinate runs from \a from to \a to in \a nSteps nominal
/// steps and hits every value of \a checkpoints (and \a to) exactly
struct Phase
{
    real_t from, to;
    index_t nSteps;
    std::vector<real_t> checkpoints;
};

/// Writes the material output \a out of \a mat on \a def (evaluated at z = 0) as a ParaView file.
/// Planar assemblers cannot construct the wrinkle fields through constructStress (which needs the surface
/// normal), so the fields of planar cases are evaluated from the material directly.
template <MaterialOutput out>
void writeMaterialField(MaterialBase * mat, const gsMultiPatch<> & def, const std::string & fileName)
{
    gsMaterialMatrixEval<real_t,out> ev(mat,&def,gsMatrix<>::Zero(1,1));
    gsField<> field(def,ev,true);
    gsWriteParaview(field,fileName,5000);
}

/// A ParaView collection of one field of the assembler
struct PlotField
{
    std::string name;
    stress_type::type type;
    MaterialOutput output;   ///< equivalent material output, used for planar cases
    std::unique_ptr<gsParaviewCollection<> > collection;
    bool failed;
};

}// namespace

int main (int argc, char** argv)
{
    // Input options
    index_t testCase   = 0;
    std::string modelName("split");
    real_t eta         = 0;
    bool exact         = false;
    real_t slackmult   = 1e-6;
    bool midsurface    = false;
    real_t s2crit      = 0;
    index_t degree     = 2;
    index_t nx         = -1;
    index_t ny         = -1;
    index_t steps      = -1;
    index_t presteps   = -1;
    real_t uxEnd       = 3.0;
    real_t epsEnd      = 0.3;
    real_t thetaEnd    = 0.6;
    real_t ratioValue  = 2.548;
    real_t tolNR       = 1e-8;
    index_t maxitNR    = 50;
    std::string linSolver("LU");
    index_t maxBisect  = 8;
    index_t monitorLevel = 2;
    bool newtonDiag    = false;
    real_t diagFrom    = 0;
    bool plot          = false;
    std::string outDir;
    index_t verbose    = 0;

    gsCmdLine cmd("Wrinkling benchmarks under displacement control: shear (0), stretched Mooney-Rivlin sheet (1), twisted trapezoid (2).");

    cmd.addInt("t","testCase","Test case: 0: shear rectangle (SvK, d=2); 1: stretched Mooney-Rivlin sheet (d=2); 2: twisted trapezoid (Neo-Hookean, d=3)",testCase);
    cmd.addString("","model","Wrinkling model: none, tft or split (default split)",modelName);
    cmd.addReal("","eta","Split: stress fraction retained in the wrinkled/slack state (default 0)",eta);
    cmd.addSwitch("exact","Split: exact wrinkled-state solve for nonlinear materials",exact);
    cmd.addReal("","slackmult","Multiplier of the stiffness in the slack state, TFT and Split (default 1e-6)",slackmult);
    cmd.addSwitch("midsurface","Split: evaluate the material at the mid-surface only",midsurface);
    cmd.addReal("","s2crit","Split: minor principal stress below which a point is displayed as wrinkled (default 0; pass negative values as --s2crit -4000)",s2crit);
    cmd.addInt("p","degree","Spline degree (default 2)",degree);
    cmd.addInt("","nx","Elements in the first direction (default: 40, 8, 16 for the test cases 0, 1, 2)",nx);
    cmd.addInt("","ny","Elements in the second direction (default: 10, 4, 16)",ny);
    cmd.addInt("N","steps","Nominal number of load steps of the main phase (default: 30, 30, 60)",steps);
    cmd.addInt("","presteps","Nominal number of load steps of the first phase (default: 1 and 10 for the test cases 0 and 2)",presteps);
    cmd.addReal("","ux","Test case 0: final horizontal displacement of the top edge [mm] (default 3.0)",uxEnd);
    cmd.addReal("","eps","Test case 1: final strain (default 0.3)",epsEnd);
    cmd.addReal("","theta","Test case 2: final rotation in units of pi (default 0.6)",thetaEnd);
    cmd.addReal("","ratio","Test case 1: Mooney-Rivlin ratio c1/c2 (default 2.548)",ratioValue);
    cmd.addReal("","tolNR","Relative tolerance of the Newton solver, residual and update (default 1e-8)",tolNR);
    cmd.addInt("","maxitNR","Maximum number of Newton iterations (default 50)",maxitNR);
    cmd.addString("","solver","Linear solver of the Newton iterations: LU, CGDiagonal, SimplicialLDLT, ... (default LU)",linSolver);
    cmd.addInt("","maxbisect","Maximum number of successive halvings of a load step before the run stops (default 8)",maxBisect);
    cmd.addInt("","monitor","Point diagnostics in steps.csv: 0: none; 1: per step; 2: also per Newton iteration (default 2)",monitorLevel);
    cmd.addSwitch("newtondiag","Write per-Newton-iterate diagnostics (newton_iter.csv, newton_flips.csv, converged_tangent.csv); needs --monitor 2",newtonDiag);
    cmd.addReal("","diagfrom","With --newtondiag: attempts of the last phase with load >= this value also get the tangent spectrum analysed (default 0)",diagFrom);
    cmd.addSwitch("plot","Write ParaView files and the midline stresses",plot);
    cmd.addString("","out","Output directory (default: ./TFTWrinkling_t<case>_<model>)",outDir);
    cmd.addInt("v","verbose","0: no; 1: Newton iteration output",verbose);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    if (testCase<0 || testCase>2)
        GISMO_ERROR("Test case "<<testCase<<" not implemented");
    int model = 0;
    if      (modelName=="none")  model = 0;
    else if (modelName=="tft")   model = 1;
    else if (modelName=="split") model = 2;
    else
        GISMO_ERROR("Unknown model \""<<modelName<<"\"; use none, tft or split.");
    if (model!=2 && (eta!=0 || exact || midsurface || s2crit!=0))
        gsWarn<<"--eta, --exact, --midsurface and --s2crit only act together with --model split and are ignored.\n";
    if (model==2 && exact && testCase==0)
        gsInfo<<"The shear case uses a linear (SvK) material, for which --exact has no effect.\n";
    if (degree<1)
        GISMO_ERROR("The degree must be at least 1.");
    if (newtonDiag && monitorLevel<2)
        GISMO_ERROR("--newtondiag needs --monitor 2: the state flips come from the monitor's iterate states.");
    if (nx<0) nx = (testCase==0 ? 40 : testCase==1 ? 8 : 16);
    if (ny<0) ny = (testCase==0 ? 10 : testCase==1 ? 4 : 16);
    if (steps<0) steps = (testCase==2 ? 60 : 30);
    if (presteps<0) presteps = (testCase==2 ? 10 : 1);
    if (nx<1 || ny<1 || steps<1 || presteps<1 || maxBisect<0)
        GISMO_ERROR("The numbers of elements and steps must be positive.");
    if (outDir.empty())
        outDir = "./TFTWrinkling_t" + std::to_string(testCase) + "_" + modelName;
    if (!makeDirs(outDir))
    {
        gsWarn<<"Could not create the output directory "<<outDir<<"\n";
        return 1;
    }

    gsStopwatch totalTime, stageTime;
    stageTime.restart();

    // Geometry and material parameters. Test case 0 and 1 are plane (d = 2), test case 2 is a membrane in space.
    const short_t d = (testCase==2) ? 3 : 2;
    real_t E_modulus = 0, PoissonRatio = 0, thickness = 0, ratio = 0;
    index_t materialId = 0;
    if (testCase==0)
    {
        E_modulus = 3500; PoissonRatio = 0.31; thickness = 0.025; // N/mm^2, mm
        materialId = 0;
    }
    else if (testCase==1)
    {
        E_modulus = 1.32e6; PoissonRatio = 0.5; thickness = 0.14e-3; ratio = ratioValue; // Pa, m
        materialId = 3;
        const real_t mu = E_modulus/(2*(1+PoissonRatio));
        const real_t c2 = mu/(ratio+1);
        const real_t c1 = ratio*c2;
        gsInfo<<"Mooney-Rivlin: mu = "<<mu<<" Pa, Ratio = "<<ratio<<", c1 = "<<c1<<" Pa, c2 = "<<c2<<" Pa\n";
    }
    else
    {
        E_modulus = 1e3; PoissonRatio = 0.5; thickness = 5e-3; // N/mm^2, mm
        materialId = 1;
    }

    // Single patch with the Greville points of the spline basis as control points; (bi)linear maps are reproduced exactly
    gsKnotVector<> kx(0,1,nx-1,degree+1), ky(0,1,ny-1,degree+1);
    gsTensorBSplineBasis<2,real_t> tb(kx,ky);
    const gsMatrix<> anchors = tb.anchors();
    gsMatrix<> coefs(tb.size(),d);
    coefs.setZero();
    for (index_t i = 0; i!=coefs.rows(); ++i)
    {
        const real_t xi = anchors(0,i), et = anchors(1,i);
        if      (testCase==0)
        {
            coefs(i,0) = 380.*xi;
            coefs(i,1) = 128.*et;
        }
        else if (testCase==1)
        {
            coefs(i,0) = 0.14*xi;
            coefs(i,1) = 0.07*et;
        }
        else
        {
            coefs(i,0) = (xi-0.5)*(50.+50.*et);
            coefs(i,1) = 100.*et - 50.;
        }
    }
    gsMultiPatch<> mp, mp_def;
    mp.addPatch(gsTensorBSpline<2,real_t>(tb,coefs));
    mp_def = mp;
    gsMultiBasis<> dbasis(mp,true);
    gsInfo<<"Test case "<<testCase<<", model "<<modelName<<", basis (patch 0): "<<mp.patch(0).basis()<<"\n";

    // Boundary conditions. The functions in the conditions are updated in every load step.
    gsConstantFunction<> uxFun(0.0,d), uyFun(0.0,d);
    // Test case 2: bottom edge. Phase 1 (parameter t = s): u = (0,-50 s,0). Phase 2 (t = theta/pi): rotation about
    // the vertical axis through x = 0 at the vertical displacement -50, u = (x (cos(theta)-1), -50, x sin(theta)).
    gsFunctionExpr<> p1y("-50*t",3);
    gsConstantFunction<> p2y(-50.0,3);
    gsFunctionExpr<> p2x("x*(cos(pi*t)-1)",3), p2z("x*sin(pi*t)",3);
    // Rotation velocity du/dtheta of the bottom edge
    gsFunctionExpr<> vx("-x*sin(pi*t)",3), vz("x*cos(pi*t)",3);

    gsBoundaryConditions<> BCs, BCs2;
    if (testCase==0)
    {
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 0);
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1);
        BCs.addCondition(boundary::north, condition_type::dirichlet, &uxFun, 0, false, 0);
        BCs.addCondition(boundary::north, condition_type::dirichlet, &uyFun, 0, false, 1);
    }
    else if (testCase==1)
    {
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0);
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1);
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 1);
        BCs.addCondition(boundary::east, condition_type::dirichlet, &uxFun, 0, false, 0);
    }
    else
    {
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, -1);
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 0);
        BCs.addCondition(boundary::south, condition_type::dirichlet, &p1y, 0, false, 1);
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2);

        BCs2.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, -1);
        BCs2.addCondition(boundary::south, condition_type::dirichlet, &p2x, 0, false, 0);
        BCs2.addCondition(boundary::south, condition_type::dirichlet, &p2y, 0, false, 1);
        BCs2.addCondition(boundary::south, condition_type::dirichlet, &p2z, 0, false, 2);
        BCs2.setGeoMap(mp);
    }
    BCs.setGeoMap(mp);

    // Materials
    gsConstantFunction<> t(thickness,d);
    gsConstantFunction<> E(E_modulus,d);
    gsConstantFunction<> nu(PoissonRatio,d);
    gsConstantFunction<> rho(1.0,d);
    gsConstantFunction<> ratioFun(ratio,d);

    std::vector<gsFunctionSet<>*> parameters;
    parameters.push_back(&E);
    parameters.push_back(&nu);
    if (materialId==3)
        parameters.push_back(&ratioFun);

    gsOptionList options;
    options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",materialId);
    options.addSwitch("Compressibility","Compressibility: (false): Imcompressible | (true): Compressible",false);
    options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);

    RunMaterials materials;
    if      (testCase==0)
        buildMaterials<2,true >(materials,mp,t,rho,parameters,options,model,eta,exact,midsurface,s2crit,slackmult);
    else if (testCase==1)
        buildMaterials<2,false>(materials,mp,t,rho,parameters,options,model,eta,exact,midsurface,s2crit,slackmult);
    else
        buildMaterials<3,false>(materials,mp,t,rho,parameters,options,model,eta,exact,midsurface,s2crit,slackmult);

    gsVector<> forceVec = gsVector<>::Zero(d);
    gsConstantFunction<> force(forceVec,d);

    // The assembler integrates the material objects as they are (it does not clone them)
    std::unique_ptr<gsThinShellAssemblerBase<real_t> > assembler;
    if (d==2)
        assembler.reset(new gsThinShellAssembler<2,real_t,false>(mp,dbasis,BCs,force,materials.run));
    else
        assembler.reset(new gsThinShellAssembler<3,real_t,false>(mp,dbasis,BCs,force,materials.run));
    GISMO_ENSURE(assembler->material(0) == materials.run,"The assembler does not use the expected material object");

    // Sample points of the diagnostics: the Gauss points of the assembler
    gsOptionList quOpt;
    quOpt.addReal("quA","",1.0);
    quOpt.addInt("quB","",1);
    const gsMatrix<> sampleNodes = gsQuadrature::getAllNodes(dbasis.basis(0),quOpt);
    WrinkleMonitor monitor(materials.oracle.get(),sampleNodes,model==2 && midsurface,monitorLevel);

    std::unique_ptr<NewtonDiag> diag;
    if (newtonDiag)
        diag.reset(new NewtonDiag(outDir,sampleNodes,diagFrom));

    assembler->assemble();
    gsSparseMatrix<> K = assembler->matrix();
    gsVector<> F = assembler->rhs();

    gsStructuralAnalysisOps<real_t>::Jacobian_t Jacobian = [&assembler,&mp_def,&monitor,&diag](gsVector<real_t> const &x, gsSparseMatrix<real_t> & m)
    {
        assembler->constructSolution(x,mp_def);
        ThinShellAssemblerStatus status = assembler->assembleMatrix(mp_def);
        monitor.observeIterate(mp_def);
        m = assembler->matrix();
        if (diag)
            diag->observeJacobian(x,m,monitor.iterateState());
        return status == ThinShellAssemblerStatus::Success;
    };
    gsStructuralAnalysisOps<real_t>::Residual_t Residual = [&assembler,&mp_def,&diag](gsVector<real_t> const &x, gsVector<real_t> & result)
    {
        assembler->constructSolution(x,mp_def);
        ThinShellAssemblerStatus status = assembler->assembleVector(mp_def);
        result = assembler->rhs();
        if (diag)
            diag->observeResidual(x,result);
        return status == ThinShellAssemblerStatus::Success;
    };

    gsStaticNewton<real_t> NWT(K,F,Jacobian,Residual);
    gsOptionList NWTOptions = NWT.options();
    NWTOptions.setInt("maxIt",maxitNR);
    NWTOptions.setReal("tol",tolNR);
    NWTOptions.setString("Solver",linSolver);
    NWTOptions.setInt("verbose",verbose);
    NWT.setOptions(NWTOptions);

    // The materials of the assembler, for the virtual work integral
    const gsMaterialMatrixContainer<real_t> materialContainer = assembler->materials();

    gsInfo<<"Setup time: "<<stageTime.stop()<<" s; degrees of freedom: "<<assembler->numDofs()<<"\n";

    // Load schedule. The load coordinate D of a phase is
    //   case 0: phase 1: u_y of the top edge [mm]; phase 2: u_x of the top edge [mm]
    //   case 1: the strain eps
    //   case 2: phase 1: the fraction s of the 50 mm downward displacement; phase 2: theta/pi
    std::vector<Phase> phases;
    if (testCase==0)
    {
        Phase p1 = {0.0,0.05,presteps,std::vector<real_t>()};
        Phase p2 = {0.0,uxEnd,steps,std::vector<real_t>()};
        if (uxEnd>1.6+1e-12)
            p2.checkpoints.push_back(1.6);
        phases.push_back(p1);
        phases.push_back(p2);
    }
    else if (testCase==1)
    {
        Phase p1 = {0.0,epsEnd,steps,std::vector<real_t>()};
        phases.push_back(p1);
    }
    else
    {
        Phase p1 = {0.0,1.0,presteps,std::vector<real_t>()};
        Phase p2 = {0.0,thetaEnd,steps,std::vector<real_t>()};
        phases.push_back(p1);
        phases.push_back(p2);
    }
    for (size_t ph = 0; ph!=phases.size(); ++ph)
        phases[ph].checkpoints.push_back(phases[ph].to);

    // Sets the prescribed values of phase \a ph at the load coordinate \a D
    auto applyLoad = [&](size_t ph, real_t D)
    {
        if (testCase==0)
        {
            uyFun.setValue(ph==0 ? D : 0.05,d);
            uxFun.setValue(ph==0 ? 0.0 : D,d);
            assembler->updateBCs(BCs);
        }
        else if (testCase==1)
        {
            uxFun.setValue(D*0.14,d);
            assembler->updateBCs(BCs);
        }
        else if (ph==0)
        {
            p1y.set_t(D);
            assembler->updateBCs(BCs);
        }
        else
        {
            p2x.set_t(D);
            p2z.set_t(D);
            vx.set_t(D);
            vz.set_t(D);
            assembler->updateBCs(BCs2);
        }
    };

    // Case-specific columns of steps.csv, after the columns of WrinkleMonitor::header:
    //   all cases:  phase (-)           1 or 2
    //   test case 0 (mm, N, MPa):
    //     ux, uy               prescribed displacement of the top edge (mm)
    //     Rx, Ry               reaction force on the top edge, components x and y (N)
    //     S1mid_mean, S1mid_max, S2mid_mean
    //                          mean and maximum of the major and mean of the minor principal stress of the run
    //                          material over 101 points of the midline y = H/2, x in [L/4, 3L/4] (MPa)
    //   test case 1 (m, N, Pa):
    //     eps                  prescribed strain (-)
    //     Rx_quarter           reaction force of the quarter model on the grip (N)
    //     Rx_full              reaction force of the full sheet on the grip, twice the quarter value (N)
    //     wwidth               width of the wrinkled zone along the vertical centre line (m): twice the length
    //                          of the centre line of the quarter on which the run material is wrinkled, sampled
    //                          at 201 points; -1 without wrinkling model
    //   test case 2 (mm, N, MPa):
    //     theta_pi             rotation of the bottom edge divided by pi (-)
    //     uy_bot               vertical displacement of the bottom edge (mm)
    //     Torque               torque on the rotated bottom edge (N mm); NaN in phase 1
    //     Fy_bot               vertical reaction force on the bottom edge (N)
    // Response columns are NaN for attempts that did not converge.
    std::vector<std::string> qoiNames;
    if      (testCase==0)
        qoiNames = {"ux","uy","Rx","Ry","S1mid_mean","S1mid_max","S2mid_mean"};
    else if (testCase==1)
        qoiNames = {"eps","Rx_quarter","Rx_full","wwidth"};
    else
        qoiNames = {"theta_pi","uy_bot","Torque","Fy_bot"};

    const real_t nan = std::numeric_limits<real_t>::quiet_NaN();
    const gsMatrix<> Z0 = gsMatrix<>::Zero(1,1);

    std::ofstream stepsFile(outDir + "/steps.csv");
    if (!stepsFile)
    {
        gsWarn<<"Could not open "<<outDir<<"/steps.csv for writing\n";
        return 1;
    }
    stepsFile<<std::setprecision(12)<<WrinkleMonitor::header()<<",phase";
    for (size_t i = 0; i!=qoiNames.size(); ++i)
        stepsFile<<","<<qoiNames[i];
    stepsFile<<std::endl;

    // ParaView output
    gsParaviewCollection<> solcollection(outDir + "/solution");
    std::vector<PlotField> plotFields;
    if (plot)
    {
        std::vector<PlotField> fieldList;
        if (model==0)
            fieldList.push_back(PlotField{"principal_stress_membrane",stress_type::principal_stress_membrane,MaterialOutput::PStress,nullptr,false});
        else
        {
            fieldList.push_back(PlotField{"tension_field",stress_type::tension_field,MaterialOutput::TensionField,nullptr,false});
            fieldList.push_back(PlotField{"wrinkle_theta",stress_type::wrinkle_theta,MaterialOutput::Theta,nullptr,false});
            fieldList.push_back(PlotField{"wrinkle_dir",stress_type::wrinkle_dir,MaterialOutput::WrinkleDir,nullptr,false});
            fieldList.push_back(PlotField{"wrinkle_intensity",stress_type::wrinkle_intensity,MaterialOutput::WrinkleIntensity,nullptr,false});
            fieldList.push_back(PlotField{"wrinkle_gamma_hat",stress_type::wrinkle_gamma_hat,MaterialOutput::GammaHat,nullptr,false});
        }
        for (size_t i = 0; i!=fieldList.size(); ++i)
        {
            PlotField pf;
            pf.name = fieldList[i].name;
            pf.type = fieldList[i].type;
            pf.output = fieldList[i].output;
            pf.collection.reset(new gsParaviewCollection<>(outDir + "/" + pf.name));
            pf.failed = false;
            plotFields.push_back(std::move(pf));
        }
    }

    // Initial state: undeformed
    gsVector<> solVector = gsVector<>::Zero(assembler->numDofs());
    index_t step = 0, attempt = 0, failures = 0, singular = 0, bisected = 0;
    bool failed = false;
    gsStopwatch stopwatch;

    for (size_t ph = 0; ph!=phases.size() && !failed; ++ph)
    {
        const Phase & P = phases[ph];
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

            gsInfo<<"Load step "<<step<<" (attempt "<<attempt<<", phase "<<ph+1<<"): D = "<<D<<"; dD = "<<D-Dold<<"\n";
            applyLoad(ph,D);

            // Initial guess: the previous converged state; the first step of the trapezoid starts from a stretched
            // flat state, since the unstressed flat membrane has no out-of-plane stiffness
            gsVector<> U0 = solVector;
            if (testCase==2 && step==0)
            {
                // constructSolutionVector takes the displacement coefficients
                gsMultiPatch<> dispPred = mp;
                dispPred.patch(0).coefs().setZero();
                for (index_t i = 0; i!=anchors.cols(); ++i)
                    dispPred.patch(0).coefs()(i,1) = -50.*D*(1.-anchors(1,i));
                U0 = assembler->constructSolutionVector(dispPred);
            }

            NWT.initialize();
            NWT.setDisplacement(U0);
            monitor.beginSolve();
            if (diag)
                diag->beginAttempt(attempt,step,ph+1,D,diag->inWindow(ph,phases.size(),D));
            materials.resetCounters();
            stopwatch.restart();
            NWT.solve();
            const real_t solveTime = stopwatch.stop();
            if (diag)
                diag->endAttempt();
            gsInfo<<"Solve time: "<<solveTime<<" s\n";

            StepInfo info;
            info.step = step;
            info.attempt = attempt++;
            info.load = D;
            info.converged = NWT.converged();
            info.statusNR = static_cast<int>(NWT.status());
            info.statusDR = -1;
            info.itNR = NWT.iterations();
            info.itDR = -1;
            info.usesNR = true;
            info.time = solveTime;
            // The counters have to be read before any further evaluation of the run material
            materials.readCounters(info.splitCounters);
            if (!info.converged)
                failures++;
            if (info.statusNR==static_cast<int>(gsStatus::SolverError))
                singular++;
            info.bisections = bisected;
            info.failures = failures;
            info.singular = singular;

            if (info.converged)
            {
                solVector = NWT.solution();
                mp_def = assembler->constructSolution(solVector);
            }
            const PointStats pointStats = monitor.evaluate(mp_def,materials.run,materials.model);
            if (info.converged)
                monitor.accept();

            // Prescribed quantities and responses of the case
            std::vector<real_t> qoi(qoiNames.size(),nan);
            std::vector<real_t> midX, midS1, midS2;
            if (testCase==0)
            {
                qoi[0] = ph==0 ? 0.0 : D;
                qoi[1] = ph==0 ? D : 0.05;
            }
            else if (testCase==1)
            {
                qoi[0] = D;
            }
            else
            {
                qoi[0] = ph==0 ? 0.0 : D;
                qoi[1] = ph==0 ? -50.*D : -50.;
            }

            if (info.converged)
            {
                if (testCase==0)
                {
                    const gsMatrix<> bf = assembler->boundaryForce(mp_def,patchSide(0,boundary::north));
                    qoi[2] = -bf(0,0);
                    qoi[3] = -bf(1,0);
                    midlineStress(materials.run,mp_def,380.,101,midX,midS1,midS2);
                    qoi[4] = std::accumulate(midS1.begin(),midS1.end(),0.0)/midS1.size();
                    qoi[5] = *std::max_element(midS1.begin(),midS1.end());
                    qoi[6] = std::accumulate(midS2.begin(),midS2.end(),0.0)/midS2.size();
                }
                else if (testCase==1)
                {
                    const gsMatrix<> bf = assembler->boundaryForce(mp_def,patchSide(0,boundary::east));
                    qoi[1] = -bf(0,0);
                    qoi[2] = 2*qoi[1];
                    qoi[3] = -1;
                    if (model!=0)
                    {
                        const index_t nLine = 201;
                        gsMatrix<> line(2,nLine);
                        line.row(0).setZero();
                        line.row(1).setLinSpaced(nLine,0,1);
                        const gsMatrix<> state = evalMaterial<MaterialOutput::TensionField>(materials.run,mp_def,line,Z0);
                        index_t nWrinkled = 0;
                        for (index_t i = 0; i!=state.cols(); ++i)
                            if (state(0,i)==0)
                                ++nWrinkled;
                        qoi[3] = 2*nWrinkled*0.07/(nLine-1);
                    }
                }
                else
                {
                    qoi[2] = nan;
                    if (ph==1)
                        qoi[2] = virtualWorkTorque(dbasis,mp,mp_def,materialContainer,vx,vz);
                    const gsMatrix<> bf = assembler->boundaryForce(mp_def,patchSide(0,boundary::south));
                    qoi[3] = -bf(1,0);
                }
            }

            monitor.writeColumns(stepsFile,info,pointStats);
            stepsFile<<","<<ph+1;
            for (size_t i = 0; i!=qoi.size(); ++i)
                stepsFile<<","<<qoi[i];
            stepsFile<<std::endl;

            if (diag && info.converged)
            {
                assembler->assembleMatrix(mp_def);
                diag->writeConverged(info.step,ph+1,D,pointStats,assembler->matrix());
            }

            if (!info.converged)
            {
                gsWarn<<"Load step "<<step<<" did not converge (status "<<info.statusNR<<")\n";
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

            // ! [Export visualization in ParaView]
            if (plot)
            {
                gsMultiPatch<> deformation = mp_def;
                deformation.patch(0).coefs() -= mp.patch(0).coefs();
                gsField<> solField(mp_def,deformation);
                std::string fileName = outDir + "/" + "solution" + util::to_string(step);
                gsWriteParaview<>(solField,fileName,1000,true);
                solcollection.addPart(gsFileManager::getFilename(fileName) + "0.vts",step);

                for (size_t i = 0; i!=plotFields.size(); ++i)
                {
                    PlotField & pf = plotFields[i];
                    if (pf.failed)
                        continue;
                    try
                    {
                        fileName = outDir + "/" + pf.name + util::to_string(step);
                        if (d==2)
                        {
                            switch (pf.output)
                            {
                                case MaterialOutput::PStress:           writeMaterialField<MaterialOutput::PStress>(materials.run,mp_def,fileName); break;
                                case MaterialOutput::TensionField:      writeMaterialField<MaterialOutput::TensionField>(materials.run,mp_def,fileName); break;
                                case MaterialOutput::Theta:             writeMaterialField<MaterialOutput::Theta>(materials.run,mp_def,fileName); break;
                                case MaterialOutput::WrinkleDir:        writeMaterialField<MaterialOutput::WrinkleDir>(materials.run,mp_def,fileName); break;
                                case MaterialOutput::WrinkleIntensity:  writeMaterialField<MaterialOutput::WrinkleIntensity>(materials.run,mp_def,fileName); break;
                                default:                                writeMaterialField<MaterialOutput::GammaHat>(materials.run,mp_def,fileName); break;
                            }
                        }
                        else
                        {
                            gsPiecewiseFunction<> pieces;
                            assembler->constructStress(mp_def,pieces,pf.type);
                            gsField<> field(mp_def,pieces,true);
                            gsWriteParaview(field,fileName,5000);
                        }
                        pf.collection->addPart(gsFileManager::getFilename(fileName) + "0.vts",step);
                    }
                    catch (...)
                    {
                        gsWarn<<"The field "<<pf.name<<" could not be constructed for this case and is skipped\n";
                        pf.failed = true;
                    }
                }
            }

            if (plot && testCase==0 && ph==1 && hitsCheckpoint)
                writeMidline(outDir + "/midline_" + util::to_string(step) + ".csv",midX,midS1,midS2);

            Dold = D;
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

    solcollection.save();
    for (size_t i = 0; i!=plotFields.size(); ++i)
        if (!plotFields[i].failed)
            plotFields[i].collection->save();

    gsInfo<<"Total time: "<<totalTime.stop()<<" s; accepted steps: "<<step<<"; attempts: "<<attempt<<"; halvings: "<<bisected<<"\n";
    return failed ? 1 : 0;
}
#else//gsKLShell_ENABLED
int main(int argc, char *argv[])
{
    GISMO_UNUSED(argc); GISMO_UNUSED(argv);
    gsWarn<<"G+Smo is not compiled with the gsKLShell module.";
    return EXIT_FAILURE;
}
#endif
