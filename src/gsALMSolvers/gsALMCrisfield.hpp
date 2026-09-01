/** @file gsALMCrisfield.hpp

    @brief Performs the arc length method to solve a nonlinear equation system.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#pragma once

#include <typeinfo>
#include <limits>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMHelper.h>

namespace gismo
{

template <class T>
void gsALMCrisfield<T>::defaultOptions()
{
    Base::defaultOptions();
    // THE DEFAULT IS Scaling = -1: the AUTOMATIC state-dependent phi, i.e. Lam & Morley 1992
    // eq. (11) A0 = p0'p0/lambda0^2, with their eq. (22) A0 = dq'dq at the origin exception --
    // see predictor(). p0 is a DISPLACEMENT, so A0 = |U|^2/Lambda^2 here. -1 is the ONLY
    // value for which m_phi_user is false, so any other value is used verbatim as a fixed psi.
    //
    // Scaling = 0 gives the CYLINDRICAL constraint |DeltaU|^2 = ds^2, i.e. A0 = 0: the load
    // term is dropped from the metric. It is a legitimate option and most in-tree call sites
    // pin it themselves, but it is NOT the shipped default -- see the literature note below.
    //
    // THE FAILING CONFIGURATION is Scaling = 0 with BorderedFallback = false: fixture F's
    // fundamental path (u2 == 0, A0 == 0) advances u1 by EXACTLY ds, so a step lands
    // head-on on the fold at u1 = 1.0 (det K = 0), the elimination's discriminant cancels
    // to exactly 0 (see (a) on the class doc), and crisfield_root_selection_never_retraces
    // / crisfield_fold_rounding_is_scale_invariant fail. Two independent ways out: the
    // shipped Scaling = -1 (the automatic phi avoids the head-on step), or
    // BorderedFallback = true (the bordered chart rounds the fold). Turning
    // BorderedFallback ON at the shipped Scaling = -1 reopens nothing. The Scaling = 0
    // precondition here is a fixed (state-independent) metric: under Scaling = -1 phi is
    // RECOMPUTED every step from the
    // state (see predictor()), so the constraint metric itself changes step to step and the
    // fixture's ds is only commensurate with the fold location under the fixed metric.
    //
    // Literature is split on cylindrical (Scaling = 0) vs. the automatic phi. Crisfield 1981
    // and Bellini & Chulya 1987 recommend cylindrical; against it: Lam & Morley 1992 (p.173,
    // "Many analysts set A0 = 0, BUT WE TAKE A0 = p0'p0/lambda0^2"), Crisfield 1983 ("a
    // variable scaling could be usefully developed for the load term"), Schweizerhof &
    // Wriggers 1986 ("none of the above-mentioned schemes can be recommended as best ... the
    // load term HAS TO BE PRESENT if stiffening structures are analyzed"), and Ritto-Correa &
    // Camotim 2008 (Sec. 5: "a cylindrical constraint is NOT SUITABLE to compute this kind of
    // 'unfolded' equilibrium paths"). No surveyed source endorses eq. (11) as a shipped
    // default either. Several in-tree drivers and XML configs pin Scaling explicitly (both
    // -1 and 0), so the default has real blast radius: do not flip it without a fresh,
    // two-sided literature case and the user's sign-off.
    m_options.addReal("Scaling","Set Scaling factor Phi (-1 = automatic, Lam & Morley 1992; 0 = cylindrical)",-1);
    m_options.addInt ("AngleMethod","Angle determination method: 0 = Previous step; 1 = Previous iteration",angmethod::Step);
    // BorderedFallback defaults OFF (opt-in). The MECHANISM is sound and stays in the
    // library, fully gated: the G-A test, and the adversarial re-derivation in
    // MATH-VERIFICATION-adversarial.md, re-derive every bordered-chart identity on exact
    // rationals. It is not the default because no in-tree driver pins it, so an ON default
    // would reach every consumer of this class, the large majority never measured with it on.
    // Turning it ON is a good choice for a caller that hits fold-step failures, but should
    // remain an informed opt-in until the real consumer population has been measured.
    //
    // DEPRECATED: superseded by the string option
    // BorderedMode below, which can also select the bordered chart as the PRIMARY corrector
    // (the references' own practice) rather than only as a fallback. Kept registered and
    // functional -- removing it would break existing options().setSwitch("BorderedFallback",...)
    // call sites -- and mapped onto BorderedMode="Fallback" by getOptions() with a
    // once-per-object warning; see getOptions() for the precedence rule.
    m_options.addSwitch("BorderedFallback","DEPRECATED, use BorderedMode=\"Fallback\" instead. Retry a NotConverged step with the bordered corrector chart (see gsALMCrisfield)",false);
    m_options.addString("BorderedSolver","Sparse solver for the bordered (n+1)x(n+1) system; must PIVOT (LU or QR)","LU");
    // Corrector chart selector -- see the class documentation and step(). "Off" (default) is
    // today's shipped behaviour, unchanged in every respect. "Fallback" attempts the
    // elimination first and retries with the bordered chart on failure (today's
    // BorderedFallback=true). "Primary" attempts the bordered chart FIRST and retries with the
    // elimination -- the references' own practice (AUTO and pde2path both use the bordered
    // (n+1) solve as THE corrector, not a fallback).
    m_options.addString("BorderedMode","Corrector chart: Off (elimination only), Fallback (elimination, bordered retry), Primary (bordered, elimination retry)","Off");
    // Seeded here, NOT in initMethods() and NOT in getOptions(): it must survive repeated
    // getOptions() calls (every gsALMBase::applyOptions() call re-runs getOptions()), so a
    // caller that calls applyOptions() three times still gets exactly one warning.
    m_borderedAliasWarned = false;
}

template <class T>
void gsALMCrisfield<T>::getOptions()
{
    Base::getOptions();
    m_phi                 = m_options.getReal("Scaling");
    m_phi_user = m_phi == -1 ? false : true;

    m_angleDetermine      = m_options.getInt ("AngleMethod");

    // Resolve BorderedMode (string, exact match, case sensitive). An unrecognised literal is
    // a hard error -- GISMO_ASSERT is compiled out under -O3 -DNDEBUG (see the class C8 note),
    // so this MUST be GISMO_ENSURE.
    const std::string borderedModeStr = m_options.getString("BorderedMode");
    if (borderedModeStr == "Off")
      m_borderedMode = borderedmode::Off;
    else if (borderedModeStr == "Fallback")
      m_borderedMode = borderedmode::Fallback;
    else if (borderedModeStr == "Primary")
      m_borderedMode = borderedmode::Primary;
    else
      GISMO_ENSURE(false,"gsALMCrisfield: unknown BorderedMode \""<<borderedModeStr<<"\"; "
                         "expected one of \"Off\", \"Fallback\", \"Primary\".");

    // Deprecated alias: BorderedFallback=true maps onto BorderedMode="Fallback", but ONLY
    // when BorderedMode itself was left at "Off" -- an explicitly non-"Off" BorderedMode
    // WINS over the alias, so a caller migrating to the new option is never silently
    // downgraded. Warn once per solver object (m_borderedAliasWarned, seeded in
    // defaultOptions() so it survives repeated getOptions() calls).
    if (m_borderedMode == borderedmode::Off && m_options.getSwitch("BorderedFallback"))
    {
      m_borderedMode = borderedmode::Fallback;
      if (!m_borderedAliasWarned)
      {
        gsWarn<<"gsALMCrisfield: option BorderedFallback is DEPRECATED; use BorderedMode=\"Fallback\". Mapping it for this run.\n";
        m_borderedAliasWarned = true;
      }
    }

    // The solver object is built lazily and then CACHED, so a mid-analysis change of the
    // option would otherwise be silently ignored forever. Release it when the name moves.
    const std::string borderedSolverName = m_options.getString("BorderedSolver");
    if (borderedSolverName != m_borderedSolverName)
      m_borderedSolver.reset();
    m_borderedSolverName  = borderedSolverName;
    // Never live outside a bordered attempt of step(); re-asserted here so that re-applying
    // options mid-analysis cannot leave the corrector in the bordered chart.
    m_useBordered         = false;
}

template <class T>
void gsALMCrisfield<T>::initMethods()
{
  m_numDof = m_forcing.size();
  m_DeltaU = m_U = gsVector<T>::Zero(m_numDof);
  m_DeltaL = m_L = 0.0;

  m_DeltaUold = gsVector<T>::Zero(m_numDof);
  m_DeltaLold = 0.0;

  // BORDERED-ONLY, POISONED on purpose -- see the member documentation. Not 1 and 0.
  m_deltaLt = m_deltaLbar = std::numeric_limits<T>::quiet_NaN();

  // Transversality guard, see the class documentation (h). The cosine is BORDERED-ONLY and
  // poisoned for the same reason as the two above; the counter is cumulative and starts at 0.
  m_chartCosine       = std::numeric_limits<T>::quiet_NaN();
  m_chartReborderings = 0;
  // ⚠ OURS AND UNSOURCED -- the two thresholds are documented at their declaration.
  m_chartCosTol       = 1e-6;
  m_chartRankTol      = 1e-10;
}

// ------------------------------------------------------------------------------------------------------------
// ---------------------------------------Crisfield's method---------------------------------------------------
// ------------------------------------------------------------------------------------------------------------

template <class T>
void gsALMCrisfield<T>::quasiNewtonPredictor()
{
  m_jacMat = computeJacobian();
  this->factorizeMatrix(m_jacMat);
  computeUt(); // rhs is state-independent for a dead load; with a forcing callback set it IS state-dependent (see gsALMBase::setForcingFunction)
  computeUbar(); // rhs contains residual and should be computed every time

}

template <class T>
void gsALMCrisfield<T>::quasiNewtonIteration()
{
  if (m_useBordered)
  {
    // The bordered chart never inverts K on its own, so K is assembled but NOT factorized
    // (factorizing it is what would throw 3 on an exactly singular tangent), and there is
    // no delta_u_t: its role is taken by (p_t,q_t), computed in iteration().
    m_jacMat = this->_assembleJacobianUnfactorized();
    return;
  }
  m_jacMat = computeJacobian();
  this->factorizeMatrix(m_jacMat);
  computeUt(); // rhs is state-independent for a dead load; with a forcing callback set it IS state-dependent (see gsALMBase::setForcingFunction)
}

template <class T>
void gsALMCrisfield<T>::iteration()
{
  if (m_useBordered)
    this->_borderedSolve(); // fills m_deltaUt/m_deltaLt = (p_t,q_t) and m_deltaUbar/m_deltaLbar = (p,q)
  else
    computeUbar(); // rhs contains residual and should be computed every time

  // Compute next solution
  m_eta = 1.0;

  T lamold = m_deltaL;
  computeLambdas();

  // Relaxation against oscillating load factor
  if (( (lamold*m_deltaL < 0) && (abs(m_deltaL) <= abs(lamold) ) ) && m_relax != 1.0 )
  {
     m_note += "\t relaxated solution!";
     m_deltaU = m_relax * (m_deltaL*m_deltaUt + m_eta*m_deltaUbar);
     m_deltaL = m_relax * m_deltaL;
  }

  m_DeltaU += m_deltaU;
  m_DeltaL += m_deltaL;

  if (m_angleDetermine == angmethod::Iteration)
  {
    m_DeltaUold = m_DeltaU;
    m_DeltaLold = m_DeltaL;
  }
}

template <class T>
void gsALMCrisfield<T>::initiateStep()
{
  m_DeltaU = m_deltaUbar = m_deltaUt = gsVector<T>::Zero(m_numDof);
  m_DeltaL = m_deltaL = 0.0;
  m_eta = 1.0;
  // Freeze the constraint-metric forcing at the predictor state (m_U,m_L): the ||f||^2
  // scalings (A0) and phi must define ONE fixed quadratic constraint for the whole
  // corrector, otherwise the arc-length surface moves under the root solve. A no-op
  // without a forcing callback.
  this->_freezeStepForcing();
}

template <class T>
void gsALMCrisfield<T>::predictor()
{
  m_jacMat = computeJacobian();
  this->factorizeMatrix(m_jacMat);
  // Predictor tangent at the Jacobian's linearization state (m_U,m_L).
  m_deltaUt = this->solveSystem(this->computeForcing(m_U,m_L));

  // Choose Solution
  if (m_DeltaUold.dot(m_DeltaUold) == 0 && m_DeltaLold*m_DeltaLold == 0) // no information about previous step.
  {
    // gsWarn<<"Different predictor!!!!!\n";
    m_note+= "predictor\t";
    // Constraint-metric scaling: frozen step forcing (see gsALMBase::stepForcing).
    // MUST be computed BEFORE the predictor length below, which is the positive root of
    // |DeltaU|^2 + A0*DeltaLambda^2 = ds^2 and therefore depends on A0 = phi^2 ||f||^2
    // (Ritto-Correa & Camotim 2008, eqs. (16)/(21): DeltaLambda_P = +-L/sqrt(a^Q.a^Q + psi^2)).
    // Specialising that root to A0 = |delta_u_t|^2 -- i.e. to the automatic phi set two
    // lines below it -- would violate the constraint by
    // sqrt((|delta_u_t|^2 + A0)/(2|delta_u_t|^2)) for any user Scaling.
    // This is the same denominator computeLambdaMU() uses.
    if (!m_phi_user)
      m_phi = math::pow( m_deltaUt.dot(m_deltaUt) / this->stepForcing().dot(this->stepForcing()),0.5);
    const T A0 = math::pow(m_phi,2) * this->stepForcing().dot(this->stepForcing());
    m_deltaL = m_arcLength / math::pow( m_deltaUt.dot(m_deltaUt) + A0 , 0.5 );
    m_deltaU = m_deltaUbar + m_deltaL*m_deltaUt;
    m_note += " phi=" + std::to_string(m_phi);

  }
  // {
  //   m_note+= "predictor\t";

  //   if (!m_phi_user)
  //     m_phi = math::pow( m_deltaUt.dot(m_deltaUt) / m_forcing.dot(m_forcing),0.5);

  //   // m_deltaL = m_arcLength * DL / math::sqrt( 2*( m_deltaUt.dot(m_deltaUt) + m_DeltaL*DL ) );
  //   m_deltaL = m_arcLength / math::sqrt( ( m_deltaUt.dot(m_deltaUt) + m_phi*m_phi ) );

  //   // m_deltaU = m_arcLength * m_deltaUt / math::sqrt( m_deltaUt.dot(m_deltaUt) + m_DeltaL*DL );
  //   m_deltaU = m_deltaL*m_deltaUt;


  //   m_note += " phi=" + std::to_string(m_phi);

  //   // m_DeltaUold = m_deltaU;
  //   // m_DeltaLold = m_deltaL;
  // }
  else // previous point is not in the origin
  {
    // Constraint-metric scaling: frozen step forcing (see gsALMBase::stepForcing).
    //
    // The automatic phi implements Lam & Morley 1992's TWO-formula prescription, in their
    // own order of precedence:
    //
    //   A0 = p0'p0 / lambda0^2      eq. (11)   -- normally, at the last converged point
    //   A0 = delta_q' delta_q       eq. (22)   -- at the exception, where (11) is undefined
    //
    // In our symbols A0 = phi^2 ||f||^2, so eq. (11) is phi = |U|/(|Lambda| ||f||) and
    // eq. (22) is phi = |delta_u_t|/||f|| -- the very formula the fresh-start branch above
    // and predictorGuess() already use. The two are ONE formula and its limit, not rivals:
    // the paper derives (22) from (11) by substituting p0 = lambda0 K^-1 q = lambda0 delta_q
    // (their eq. (21)). The branch split in this function is therefore the paper's own
    // structure and not an inconsistency.
    //
    // WHY eq. (22) MUST NOT BE USED AWAY FROM THE EXCEPTION. It makes A0 = |K^-1 f|^2,
    // which diverges at a LIMIT POINT -- Ritto-Correa & Camotim 2008, footnote 3, p.1356
    // ("the norm of t^Q ... tend[s] to infinity as one approaches a limit point where the
    // slope is exactly zero"), where psi^2 (= A0) is a fixed analysis PARAMETER, never a
    // state function. Measured on the 2-DOF fold fixture at ds = 0.05: with eq. (22)
    // applied everywhere, A0 grows 3.1 -> 8.5 -> 26.8 -> 75.6
    // over steps 20..26 where eq. (11) stays at 0.70 -> 0.85 -> 0.92; at the state step 27
    // starts from (u1 = 0.984759, K11 = 2-2*u1 = 0.0305) it is A0 = 1076 against 0.970, and
    // the corrector then does not round the fold at all -- both
    // crisfield_root_selection_never_retraces and crisfield_fold_rounding_is_scale_invariant
    // die at step 27, max u1 = 0.985 < 1.05, at all three load scales. Eq. (11) never
    // touches K, so it is finite exactly where eq. (22) is not.
    //
    // ***OUR EXTENSION, NOT THE PAPER'S***. Lam & Morley except only THE ORIGIN
    // {p0 = lambda0 = 0}, and there eq. (22) genuinely is eq. (11)'s removable-singularity
    // limit, because their eq. (21) assumes p0 is small. At a zero-load crossing at FINITE
    // displacement (lambda0 = 0, p0 != 0) that assumption fails: eq. (11) really does
    // diverge, and eq. (22) is NOT its limit. NO SURVEYED SOURCE COVERS THAT POINT. We
    // fall back to eq. (22) there anyway, by our own choice, on the reasoning that it is
    // the paper's answer at the only nearby case it does treat, that it is finite at a
    // regular point (det K != 0), and that it is a pure function of the restart state
    // (m_deltaUt is recomputed at (m_U,m_L) at the top of this call), so restart
    // determinism is preserved. Do not read the citations above as mandating it.
    //
    // KNOWN OPEN ISSUES, deliberately NOT guarded here:
    //  * lambda0 near zero but NOT zero. Eq. (11) stays finite but huge and the step
    //    degenerates, with no assertion firing. No source covers it, and an epsilon clamp
    //    is NOT the answer: MEASURED, clamping |Lambda| to 1e-12 gives phi ~ 2e12,
    //    mu ~ 2.5e-14 and |DeltaU| ~ 1.3e-14 -- a step that goes nowhere while satisfying
    //    the arc-length identity tautologically. Left open on purpose.
    //  * lambda0 = 0 coinciding with det K = 0, where BOTH formulas diverge. No surveyed
    //    source covers it; we warn rather than hand back a silent non-finite phi.
    if (!m_phi_user)
    {
      const T ff    = this->stepForcing().dot(this->stepForcing());
      const T denom = math::pow(m_L,2) * ff;
      const T sec   = m_U.dot(m_U) / denom;           // Lam & Morley 1992, eq. (11)
      if (denom > 0 && math::isfinite(sec))
        m_phi = math::pow(sec,0.5);
      else                                            // Lam & Morley 1992, eq. (22)
      {
        m_phi = math::pow( m_deltaUt.dot(m_deltaUt) / ff, 0.5);
        if (!math::isfinite(m_phi))
          gsWarn<<"gsALMCrisfield: automatic Scaling is undefined at this state "
                <<"(Lambda = "<<m_L<<", |delta_u_t| = "<<m_deltaUt.norm()<<"): both "
                <<"Lam & Morley 1992 eq. (11) and eq. (22) diverge. Set the Scaling "
                <<"option explicitly (0 = cylindrical) for this analysis.\n";
      }
    }
    m_note += " phi=" + std::to_string(m_phi);
    computeLambdaMU();
  }

  // Compute Temporary updates of DeltaL and DeltaU
  m_DeltaU += m_deltaU;
  m_DeltaL += m_deltaL;

  if (m_angleDetermine == angmethod::Iteration || m_angleDetermine == angmethod::Predictor)
  {
   m_DeltaUold = m_DeltaU;
   m_DeltaLold = m_DeltaL;
  }
}

template <class T>
void gsALMCrisfield<T>::predictorGuess()
{
  GISMO_ASSERT(m_Uguess.rows()!=0 && m_Uguess.cols()!=0,"Guess is empty");

  m_jacMat = computeJacobian();
  this->factorizeMatrix(m_jacMat);
  // Predictor tangent at the Jacobian's linearization state (m_U,m_L) - NOT at the guess,
  // which only supplies a secant direction and its sign below.
  m_deltaUt = this->solveSystem(this->computeForcing(m_U,m_L));
  // Constraint-metric scaling: frozen step forcing (see gsALMBase::stepForcing)
  if (!m_phi_user)
    m_phi = math::pow( m_deltaUt.dot(m_deltaUt) / this->stepForcing().dot(this->stepForcing()),0.5);
  m_note += " phi=" + std::to_string(m_phi);

  //
  m_DeltaUold = -(m_Uguess - m_U);
  m_DeltaLold = -(m_Lguess - m_L);

  // m_DeltaUold *= m_arcLength / math::sqrt( m_deltaU.dot(m_deltaU));
  // m_DeltaLold *= m_arcLength / math::sqrt( m_deltaU.dot(m_deltaU));

  computeLambdaMU();

  m_DeltaU = m_deltaU;
  m_DeltaL = m_deltaL;

  m_Uguess.resize(0);
}

template <class T>
void gsALMCrisfield<T>::iterationFinish()
{
  m_converged = true;
  m_Uprev = m_U;
  m_Lprev = m_L;
  m_U += m_DeltaU;
  m_L += m_DeltaL;
  if (m_angleDetermine == angmethod::Step)
  {
    m_DeltaUold = m_DeltaU;
    m_DeltaLold = m_DeltaL;
  }
}

// ------------------------------------------------------------------------------------------------------------
// ---------------------------------------Lambda computations--------------------------------------------------
// ------------------------------------------------------------------------------------------------------------

template <class T>
void gsALMCrisfield<T>::computeLambdasSimple() //Ritto-Corrêa et al. 2008
{
  // Constraint metric: frozen step forcing, so that the quadratic constraint solved for
  // delta_lambda is the SAME surface in every corrector iteration of this step.
  T A0 = math::pow(m_phi,2)* this->stepForcing().dot(this->stepForcing()); // see Lam & Morley 1992, eq. (11)

  m_a0 = m_deltaUt.dot(m_deltaUt) + A0;
  m_b0 = 2*( m_deltaUt.dot(m_DeltaU) + m_DeltaL * A0 );
  m_b1 = 2*( m_deltaUbar.dot(m_deltaUt) );
  m_c0 = m_DeltaU.dot(m_DeltaU) + m_DeltaL*m_DeltaL * A0 - math::pow(m_arcLength,2);
  m_c1 = 2*( m_DeltaU.dot(m_deltaUbar) );
  m_c2 = m_deltaUbar.dot(m_deltaUbar);

  /// Calculate the coefficients of the polynomial
  m_alpha1 = m_a0;
  m_alpha2 = m_b0 + m_eta*m_b1;
  m_alpha3 = m_c0 + m_eta*m_c1 + m_eta*m_eta*m_c2;

  m_discriminant = math::pow(m_alpha2 ,2) - 4 * m_alpha1 * m_alpha3;
  // m_note += "\t D = " + std::to_string(m_discriminant);
}

template <class T>
void gsALMCrisfield<T>::computeLambdasEta()
{
  m_alpha1 = m_a0;
  m_alpha2 = m_b0 + m_eta*m_b1;
  m_alpha3 = m_c0 + m_eta*m_c1 + m_eta*m_eta*m_c2;

  // Lam & Morley 1992
  // m_discriminant = math::pow(m_alpha2 ,2) - 4 * m_alpha1 * m_alpha3;
  // m_deltaLs[0] = (-m_alpha2 + math::sqrt(m_discriminant))/(2*m_alpha1);
  // m_deltaLs[1] = (-m_alpha2 - math::sqrt(m_discriminant))/(2*m_alpha1);

  // Zhou 1995
  m_deltaLs[0] = (-m_alpha2 )/(2*m_alpha1);
  m_deltaLs[1] = (-m_alpha2 )/(2*m_alpha1);
}

template <class T>
void gsALMCrisfield<T>::computeLambdasModified()
{
  m_alpha1 = m_b1*m_b1 - 4.0*m_a0*m_c2;
  m_alpha2 = 2.0*m_b0*m_b1 - 4.0*m_a0*m_c1;
  m_alpha3 = m_b0*m_b0 - 4.0*m_a0*m_c0;

  m_discriminant = math::pow(m_alpha2 ,2.0) - 4.0 * m_alpha1 * m_alpha3;

  gsVector<T> etas(2);
  etas.setZero();
  if (m_discriminant >= 0)
  {
    etas[0] = (-m_alpha2 + math::pow(m_discriminant,0.5))/(2.0*m_alpha1);
    etas[1] = (-m_alpha2 - math::pow(m_discriminant,0.5))/(2.0*m_alpha1);

    T eta1 = std::min(etas[0],etas[1]);
    T eta2 = std::max(etas[0],etas[1]);
    if (m_verbose) {gsInfo<<"eta 1 = "<<eta1<<"\t eta2 = "<<eta2<<"\n";}

    // Approach of Zhou 1995
    // m_eta = std::min(1.0,eta2);
    // if (m_eta <= 0)
    //   gsInfo<<"Warning: both etas are non-positive!\n";
    // if (m_eta <=0.5)
    // {
    //   gsInfo<<"Warning: eta is small; bisecting step length\n";
    //   m_arcLength=m_arcLength/2.;
    // }

    // Approach of Lam & Morley 1992
    T xi = 0.05*abs(eta2-eta1);
    if (eta2<1.0)
      m_eta = eta2-xi;
    else if ( (eta2 > 1.0) && (-m_alpha2/m_alpha1 < 1.0) )
      m_eta = eta2+xi;
    else if ( (eta1 < 1.0) && (-m_alpha2/m_alpha1 > 1.0) )
      m_eta = eta1-xi;
    else if ( eta1 > 1.0 )
      m_eta = eta1 + xi;

    if (eta2<1.0)
      m_note = m_note + " option 1";
    else if ( (eta2 > 1.0) && (-m_alpha2/m_alpha1 < 1.0) )
      m_note = m_note + " option 2";
    else if ( (eta1 < 1.0) && (-m_alpha2/m_alpha1 > 1.0) )
      m_note = m_note + " option 3";
    else if ( eta1 > 1.0 )
      m_note = m_note + " option 4";

    // m_eta = eta2;
  }
  else
  {
    gsInfo<<"Discriminant was negative in modified method\n";
  }
}
template <class T>
void gsALMCrisfield<T>::computeLambdasComplex()
{
  // Complex-root fallback (Lam & Morley 1992, eqs 13-17): a point is constructed ON the same
  // constraint surface, so it uses the frozen step forcing throughout. Both occurrences
  // in Lcr must come from that single source, otherwise Lcr is no longer a load factor
  // with respect to the reference load used by A0.
  const gsVector<T> & Fref = this->stepForcing();
  T A0 = math::pow(m_phi,2)* Fref.dot(Fref); // see Lam & Morley 1992, eq. (11)
  gsVector<T> DeltaUcr = m_DeltaU + m_deltaUbar;

  // Internal force AT Lam's trial point U + DeltaU_cr, recovered from the residual:
  // R(U,Lambda) = F_int(U) - Lambda f  =>  F_int(U) = R(U,Lambda) + Lambda f.
  // K_T*(U+DeltaU) would give the internal force of a LINEAR problem
  // only, and would moreover be taken at U+DeltaU rather than at U+DeltaU_cr; both errors
  // are O(1) on a nonlinear path and would corrupt DeltaLambda_cr (the constraint itself is
  // still met exactly, because mu rescales onto the sphere afterwards, so the defect is
  // invisible in |(DU,DL)|_A -- only the load factor of the fallback point is wrong).
  // One extra residual assembly per fallback; the fallback is a rare branch.
  const gsVector<T> R = this->computeResidual(m_U + DeltaUcr, m_L + m_DeltaL);
  const gsVector<T> Fint = R + (m_L + m_DeltaL) * Fref;
  T Lcr = Fint.dot(Fref)/Fref.dot(Fref);
  T DeltaLcr = Lcr - m_L;

  T arcLength_cr = math::pow( DeltaUcr.dot(DeltaUcr) + A0 * math::pow(DeltaLcr, 2.0) ,0.5);
  T mu = m_arcLength/arcLength_cr;

  m_deltaL = mu*DeltaLcr - m_DeltaL;
  m_deltaU = mu*DeltaUcr - m_DeltaU;
}

template <class T>
void gsALMCrisfield<T>::computeLambdas()
{
  m_deltaLs.setZero(2);
  if (m_useBordered)
  {
    this->_computeLambdasBordered();
    return;
  }
  computeLambdasSimple();
  if (m_discriminant >= 0)
  {
    m_eta = 1.0;
    m_deltaLs[0] = (-m_alpha2 + math::pow(m_discriminant,0.5))/(2*m_alpha1);
    m_deltaLs[1] = (-m_alpha2 - math::pow(m_discriminant,0.5))/(2*m_alpha1);
    computeLambdaDOT();
  }
  else
  {
    m_note += "\tC";
    // Compute eta
    computeLambdasModified();
    if ((m_discriminant >= 0) && m_eta > 0.05)
    {
      // recompute lambdas with new eta
      computeLambdasEta();
      // gsInfo<<"2: dL1 = "<<m_deltaLs[0]<<"\tdL2 = "<<m_deltaLs[1]<<"\t eta = "<<m_eta<<"\n";
      computeLambdaDOT();
      // gsInfo<<"2: dL1 = "<<m_deltaL<<"\t m_deltaU.norm = "<<m_deltaU.norm()<<"\t eta = "<<m_eta<<"\n";
      if (m_verbose) {gsInfo<<"Modified Complex Root Solve\n";}
    }
    else
    {
      // if the roots of the modified method are still complex, we use the following function (see Lam & Morley 1992, eqs 13-17)
      m_eta = 1.0;
      computeLambdasComplex();
      // gsInfo<<"3: dL1 = "<<m_deltaL<<"\t m_deltaU.norm = "<<m_deltaU.norm()<<"\t eta = "<<m_eta<<"\n";
      if (m_verbose) {gsInfo<<"Simplified Complex Root Solve\n";}
      // Note: no selection of roots is needed
    }
  }
}

template <class T>
void gsALMCrisfield<T>::computeLambdaDET()
{
    computeLambdas();

    if (sign(m_DeltaL + m_deltaLs[0]) == sign(m_detKT))
      m_deltaL = m_deltaLs[0];
    else
      m_deltaL = m_deltaLs[1];

    // Compute update of U (NOTE: m_eta=1.0)
    m_deltaU = m_deltaUbar + m_deltaL*m_deltaUt;

    // gsInfo<<"\t\t Choice based on DETERMINANT. Options:\n";
    // gsInfo<<"\t\t DeltaL = "<<m_DeltaL+m_deltaLs[0]<<" DeltaU.norm = "<<(m_DeltaU + m_deltaUbar + m_deltaLs[0]*m_deltaUt).norm()<<"\n";
    // gsInfo<<"\t\t DeltaL = "<<m_DeltaL+m_deltaLs[1]<<" DeltaU.norm = "<<(m_DeltaU + m_deltaUbar + m_deltaLs[1]*m_deltaUt).norm()<<"\n";
}

template <class T>
void gsALMCrisfield<T>::computeLambdaMU()
{
    // Constraint metric: frozen step forcing (see gsALMBase::stepForcing). Called from
    // the predictors, i.e. at the very state the freeze was taken at.
    T A0 = math::pow(m_phi,2)* this->stepForcing().dot(this->stepForcing()); // see Lam & Morley 1992, eq. (11)
    index_t dir = sign(m_DeltaUold.dot(m_deltaUt) + A0*m_DeltaLold); // Feng et al. 1995 with H = psi^2 * I
    T denum = ( math::pow( m_deltaUt.dot(m_deltaUt) + A0 ,0.5) ); // Feng et al. 1995 with H = psi^2 * I

    T mu;
    if (denum==0)
      mu = m_arcLength;
    else
      mu = m_arcLength / denum;

    m_deltaL = dir*mu;
    m_deltaU = m_deltaL*m_deltaUt;
}

template <class T>
void gsALMCrisfield<T>::computeLambdaDOT()
{
    gsVector<T> deltaU1, deltaU2;
    deltaU1 = m_eta*m_deltaUbar + m_deltaUt*m_deltaLs[0];
    deltaU2 = m_eta*m_deltaUbar + m_deltaUt*m_deltaLs[1];

    // ---------------------------------------------------------------------------------
    // Method by Ritto-Corea et al. 2008
    //
    // The load term MUST be weighted by the very A0 of the constraint whose roots are
    // being selected (see computeLambdasSimple() and computeLambdaMU() in this file).
    // Ritto-Correa & Camotim 2008 eq. (22) is t = Na_A . a_Q + w^2 Na_lambda, and their
    // eqs. (9)/(11) define w^2 through the constraint L^2 = Na.Na + w^2 Na_lambda^2 as
    // the lumped scaling "that renders the product dimensionally consistent". Matching
    // that constraint term by term against m_c0/m_a0 gives w^2 == A0 = phi^2 ||f||^2.
    // Using phi^2 = A0/||f||^2 directly would add a displacement^2 to a bare load factor
    // (dimensionally inhomogeneous) and, since the automatic phi makes A0 scale free,
    // would multiply the load term by ||f||^-2. At a small load scale that weight blows
    // up and the corrector selects the root that stalls the step -- the step then fails
    // outright rather than converging (the failure mode is NOT a silent retrace).
    // Frozen step forcing, as everywhere else in this metric (gsALMBase::stepForcing).
    const T A0 = math::pow(m_phi,2)* this->stepForcing().dot(this->stepForcing()); // see Lam & Morley 1992, eq. (11)
    T DOT1,DOT2;
    DOT1 = m_deltaLs[0]*(m_DeltaUold.dot(m_deltaUt) + A0*m_DeltaLold);
    DOT2 = m_deltaLs[1]*(m_DeltaUold.dot(m_deltaUt) + A0*m_DeltaLold);

    if (DOT1 > DOT2)
    {
      m_deltaL = m_deltaLs[0];
      m_deltaU = deltaU1;
    }
    else if (DOT1 < DOT2)
    {
      m_deltaL = m_deltaLs[1];
      m_deltaU = deltaU2;
    }
    else
    {
      m_deltaL = m_deltaLs[0];
      m_deltaU = deltaU1;
    }

    // ---------------------------------------------------------------------------------
    // // Method by Crisfield 1981
    // T DOT1,DOT2;
    // DOT1 = (m_DeltaUold+deltaU1).dot(m_DeltaUold);
    // DOT2 = (m_DeltaUold+deltaU2).dot(m_DeltaUold);

    // DOT1 = (m_DeltaU+deltaU1).dot(m_DeltaUold);
    // DOT2 = (m_DeltaU+deltaU2).dot(m_DeltaUold);

    // if ((DOT1 > DOT2) && DOT2 <= 0)
    // {
    //   m_deltaL = m_deltaLs[0];
    //   m_deltaU = deltaU1;
    // }
    // else if ((DOT1 < DOT2) && DOT1 <= 0)
    // {
    //   m_deltaL = m_deltaLs[1];
    //   m_deltaU = deltaU2;
    // }
    // else if ((DOT1 >=0) && (DOT2 >=0))
    // {
    //   T linsol = -m_alpha3/m_alpha2;
    //   T diff1 = abs(m_deltaLs[0]-linsol);
    //   T diff2 = abs(m_deltaLs[1]-linsol);
    //   m_note += "\t linear solution!";
    //   // m_note += "Linsol\t" + std::to_string(diff1) + "\t" + std::to_string(diff2) + "\t" + std::to_string(linsol) + "\t" + std::to_string(m_deltaLs[0]) + "\t" + std::to_string(m_deltaLs[1]) + "\n";
    //   if (diff1 > diff2)
    //   {
    //     m_deltaL = m_deltaLs[1];
    //     m_deltaU = deltaU2;
    //   }
    //   else
    //   {
    //     m_deltaL = m_deltaLs[0];
    //     m_deltaU = deltaU1;
    //   }
    // }
    // else
    // {
    //   m_deltaL = m_deltaLs[0];
    //   m_deltaU = deltaU1;
    // }

}

// ------------------------------------------------------------------------------------------------------------
// ---------------------------------------Bordered-solve fallback----------------------------------------------
// ------------------------------------------------------------------------------------------------------------
// Read the class documentation of gsALMCrisfield first: (a) says WHY (a cancellation, not a
// singular matrix), (b) says that this is an EXACT change of basis, (c) gives the coefficient
// map and (d) the qt<0 index-order trap.

template <class T>
gsSparseMatrix<T> gsALMCrisfield<T>::_assembleJacobianUnfactorized()
{
  // gsALMBase::_computeJacobian(), which every computeJacobian() overload delegates to,
  // calls factorizeMatrix() internally and would throw 3 on an exactly singular tangent --
  // the very tangent this chart exists to survive. Same shape as the factorizeShifted()
  // ladder in gsALMBase::_extendedSystemIteration(): m_djacobian direct, throw 2 on false.
  if (m_deltaU.rows() == 0)
    m_deltaU = gsVector<T>::Zero(m_DeltaU.rows());
  gsSparseMatrix<T> J;
  m_note += "J";
  if (!m_djacobian(m_U + m_DeltaU, m_deltaU, J))
    throw 2;
  return J;
}

template <class T>
void gsALMCrisfield<T>::_borderedSolve()
{
  // Constraint metric: frozen step forcing, exactly as computeLambdasSimple() uses it.
  const gsVector<T> & Fstep = this->stepForcing();
  const T A0 = math::pow(m_phi,2) * Fstep.dot(Fstep); // see Lam & Morley 1992, eq. (11)

  // f AT THE CURRENT ITERATE, i.e. the very right-hand side computeUt() would have used
  // (see gsALMBase::computeUt()): the bordered chart must border the same operator the
  // elimination would have inverted, otherwise it is not the same solution set. Copied,
  // because with a forcing callback the returned reference aliases m_forcingEval.
  const gsVector<T> f = this->computeForcing();

  // Border row and corner = the gradient of Crisfield's OWN quadratic constraint at the
  // current iterate, so B is the exact Jacobian of the system {R = 0, c = 0} the corrector
  // is solving, and "B nonsingular" is exactly "this Crisfield step is well posed". At the
  // first corrector iteration m_DeltaU is the predictor increment, which is nonzero, so the
  // border row is never identically zero.
  const gsVector<T> w     = 2.0*m_DeltaU;
  const T           gamma = 2.0*A0*m_DeltaL;

  if (!this->_borderedChartSolve(w,gamma,f))
  {
    // NEVER throw 3 here (see (e)): a SolverError is not retried by gsALMExploration and
    // would be recorded as a converged point. This is the case a SINGULAR B is DETECTED in
    // -- the guard below exists for the complementary one, where the factorization succeeds
    // and hands back a plausible-looking chart (see (h)).
    gsWarn<<"gsALMCrisfield: the bordered system could not be factorized by '"
          <<m_borderedSolverName<<"', or the chart came back non-finite (|p_t| = "
          <<m_deltaUt.norm()<<", q_t = "<<m_deltaLt<<"). This is a "
          <<"BRANCH point or a constraint surface tangent to the equilibrium curve, which "
          <<"the corrector cannot resolve; use computeSingularPoint() there. Reporting the "
          <<"step as not converged.\n";
    throw 1;
  }

  // ------------------------------------------------------------------------------------
  // (h) THE TRANSVERSALITY GUARD. Riks 1984 eq. (2.31): J = (D/lambda')(n^t.x') with
  // D = det K. Our b0 = w.delta_u_t + gamma IS his (n^t.x') and det B = det(K)*b0, so with
  // K nonsingular the chart degenerates exactly at b0 = 0 -- a TURNING POINT WITH RESPECT
  // TO THE BORDERING DIRECTION, i.e. a property of the chart and not of the path. The test
  // is written as (cos >= tol) rather than (cos < tol) so that a NaN cosine fires too.
  // ------------------------------------------------------------------------------------
  m_chartCosine = this->_chartTransversality(w,2.0*m_DeltaL,A0);
  if (m_chartCosine >= m_chartCosTol)
    return;                            // the ordinary path: not one further line executes

  // Riks's own split, which a bare |b0| test CANNOT make: (2.20a) is "G singular, G* full
  // rank" -- incompatible equations, a removable CHART artifact -- while (2.20b) is "a rank
  // deficiency of at least 1 for BOTH G and G*", which "indicates the existence of multiple
  // solutions x' of (2.18) and thus the manifestation of BIFURCATIONS of the solution
  // curve". His summary (2.29) supplies the discriminator: at a limit point J != 0 with
  // D = 0, so b0 != 0 there; b0 ~ 0 TOGETHER with a rank-deficient K is the bifurcation
  // signature.
  //
  // ORDER MATTERS, and not for tidiness. Riks's own remedy is the first half of the test --
  // "in principle, it [case (2.20a)] can always be avoided if we take another f*_{N+1} so
  // that n* != n" -- so the guard RE-BORDERS FIRST, with a GENERIC direction (Govaerts 2000,
  // Prop. 3.2.1: "a generic choice of B, C, D will do ... could be generated by a random
  // number generator"; his optimal, kernel-based choice is the recorded improvement, and it
  // would need gsALMBase::m_V, which is sized only inside _computeCriticalMode and is stale
  // or empty on this path). The re-bordered chart spans the SAME affine solution set: the top
  // block K deltaU - deltaLambda f = -r is untouched, so (p,q) is another point of it and
  // (p_t,q_t) another parameterisation of the same null line.
  //
  // ★ AND THE ACCEPTED POINT IS UNCHANGED FOR EVERY eta, NOT ONLY FOR eta = 1 -- which is the
  // case that matters, because a tangency is exactly where the discriminant is 0 and
  // _computeLambdasBordered() takes computeLambdasModified()'s eta != 1 branch. Write the new
  // base point as (p',q') = (p,q) + s*(p_t,q_t) (it must be, both being points of the same
  // affine set with the same null line). Then the chart point is
  //     eta*p' + t*p_t = eta*p + (t + eta*s)*p_t,
  // i.e. the constraint quadratic in t is the old one SHIFTED IN ITS VARIABLE by eta*s. A shift
  // leaves the discriminant invariant, so computeLambdasModified() selects the SAME eta and the
  // accepted (deltaU,deltaLambda) is identical. The exactness invariant of (c) -- the one test
  // crisfield_bordered_chart_equals_the_elimination_away_from_the_fold pins -- is therefore not
  // touched by re-bordering.
  //
  // It is ALSO what makes the rank
  // probe below trustworthy: that probe reconstructs K^-1 v by cancelling the bordering's own
  // contribution, and through a chart with cosine c that cancellation costs a factor 1/c --
  // i.e. everything, in the very chart that triggered the guard. MEASURED: probing
  // through the DEGENERATE chart returns rank = 1 on a tangent with det K = 1e-12, because
  // there z collapses onto delta_u_t and the estimate degenerates into |f|/|delta_u_t|, which
  // is blind to a null space the load does not excite. Through the re-bordered chart the same
  // state gives rank ~ 1e-12.
  const index_t n = m_numDof;
  gsVector<T> du(n);
  T           dl = 0;
  this->_genericBorderDirection(du,dl,math::sqrt(w.dot(w) + A0*4.0*m_DeltaL*m_DeltaL),A0);

  const bool ok   = this->_borderedChartSolve(du,A0*dl,f);
  const T    cos2 = ok ? this->_chartTransversality(du,dl,A0) : (T)0;
  if (!(cos2 >= m_chartCosTol))
  {
    // A generic bordering is transversal with probability one (Govaerts 2000, Prop. 3.2.1),
    // so this is evidence that the degeneracy is NOT removable by re-bordering: G* is rank
    // deficient too, which is Riks's (2.20b) reached through his own remedy. Report; do not
    // classify further, and do not continue on a chart known to be degenerate.
    gsWarn<<"gsALMCrisfield: the bordered chart is degenerate (transversality cosine "
          <<m_chartCosine<<" < "<<m_chartCosTol<<") and a GENERIC re-bordering does not "
          <<"remove it (cosine "<<cos2<<", factorization "<<(ok?"ok":"failed")<<"). Since a "
          <<"generic bordering is transversal with probability one (Govaerts 2000, Prop. "
          <<"3.2.1), the degeneracy is not a chart artifact; this is Riks 1984 case (2.20b), "
          <<"i.e. a candidate BIFURCATION. Use computeSingularPoint()/switchBranch() at this "
          <<"state. Reporting the step as not converged.\n";
    throw 1;
  }

  T rank = 0, modeCos = 0;
  this->_chartSingularityProbe(f,rank,modeCos);

  if (rank < m_chartRankTol)
  {
    // K is rank deficient as well. The second, INDEPENDENT signal is the phi^T P criterion
    // of Wriggers & Simo 1990 eqs. (9a)/(9b) (limit point iff phi^T P != 0, bifurcation iff
    // phi^T P = 0), the same statistic gsALMBase::_testSingularPoint thresholds with
    // SingularPointTestTol -- reused here rather than duplicated. AGREEMENT (a mode
    // orthogonal to the load) is Riks's (2.20b); DISAGREEMENT means the two sourced criteria
    // point opposite ways and nothing here is entitled to pick one. Neither outcome may
    // re-border: doing so at a genuine bifurcation walks straight past it.
    if (modeCos < m_SPTestTol)
      gsWarn<<"gsALMCrisfield: the bordered chart is degenerate (transversality cosine "
            <<m_chartCosine<<" < "<<m_chartCosTol<<") AND the tangent is numerically rank "
            <<"deficient (sigma_min/|K|_F ~ "<<rank<<"), with an approximate critical mode "
            <<"orthogonal to the load (|phi.f|/|f| = "<<modeCos<<" < SingularPointTestTol = "
            <<m_SPTestTol<<"). That is the BIFURCATION signature of Riks 1984 case (2.20b) "
            <<"and of Wriggers & Simo 1990 eq. (9b); the corrector cannot resolve it and it "
            <<"is NOT re-bordered, because re-bordering would continue straight past the "
            <<"bifurcation. Use computeSingularPoint()/switchBranch() at this state. "
            <<"Reporting the step as not converged.\n";
    else
      gsWarn<<"gsALMCrisfield: the bordered chart is degenerate (transversality cosine "
            <<m_chartCosine<<" < "<<m_chartCosTol<<") and the tangent is numerically rank "
            <<"deficient (sigma_min/|K|_F ~ "<<rank<<"), but the approximate critical mode "
            <<"is NOT orthogonal to the load (|phi.f|/|f| = "<<modeCos<<" >= "
            <<"SingularPointTestTol = "<<m_SPTestTol<<"): Riks 1984 (2.29) reads this state "
            <<"as case (2.20b) while Wriggers & Simo 1990 eq. (9a) reads it as a limit "
            <<"point. The two criteria DISAGREE, so this step is NOT classified and NOT "
            <<"re-bordered. Reporting the step as not converged.\n";
    throw 1;
  }

  // Confirmed (2.20a): b0 ~ 0 on a tangent that is NOT rank deficient, and the degeneracy
  // IS removable -- the generic re-bordering above already removed it. Keep that chart.
  gsWarn<<"gsALMCrisfield: the bordered chart was turning with respect to its own bordering "
        <<"direction (transversality cosine "<<m_chartCosine<<" < "<<m_chartCosTol<<", i.e. "
        <<"Riks 1984 case (2.20a): the constraint gradient is orthogonal to the tangent, so "
        <<"det B = det(K)*b0 vanishes with det K != 0). The tangent is well conditioned "
        <<"(sigma_min/|K|_F ~ "<<rank<<"), so the step has been RE-BORDERED with a generic "
        <<"direction (cosine "<<cos2<<") and continues on the same solution set.\n";
  m_note += "\treborder";
  m_chartCosine = cos2;
  ++m_chartReborderings;
}

template <class T>
bool gsALMCrisfield<T>::_borderedChartSolve(const gsVector<T> & w, const T gamma,
                                           const gsVector<T> & f)
{
  const index_t n = m_numDof;

  gsSparseMatrix<T> B(n+1,n+1);
  gsSparseEntries<T> se;
  se.reserve(m_jacMat.nonZeros() + 2*n + 1);
  for (index_t k = 0; k != m_jacMat.outerSize(); ++k)
    for (typename gsSparseMatrix<T>::InnerIterator it(m_jacMat,k); it; ++it)
      se.add(it.row(),it.col(),it.value());
  for (index_t i = 0; i != n; ++i) if (f[i] != (T)0) se.add(i,n,-f[i]);
  for (index_t i = 0; i != n; ++i) if (w[i] != (T)0) se.add(n,i, w[i]);
  if (gamma != (T)0) se.add(n,n,gamma);
  B.setFrom(se);
  // Eigen's SparseLU needs COMPRESSED COLUMN-MAJOR storage; gsSparseMatrix is ColMajor by
  // default (see the gsSparseMatrix class documentation: "_Option zero is ColMajor order").
  B.makeCompressed();

  // Own solver, never gsALMBase::m_solver -- see the m_borderedSolver documentation.
  if (!m_borderedSolver)
    m_borderedSolver = gsSparseSolver<T>::get(m_borderedSolverName);

  m_borderedSolver->compute(B);
  // The CALLER decides what a failure means -- see the declaration. A factorization that
  // fails is a DETECTED degeneracy, reported to the caller rather than warned-and-thrown
  // here; the case (h) exists for is the one where this returns Success.
  if (m_borderedSolver->info() != gsEigen::Success)
    return false;

  // The two bordered solves, in one call: B[p_t;q_t] = [0;1] and B[p;q] = [-r;0].
  gsMatrix<T> rhs(n+1,2);
  rhs.setZero();
  rhs(n,0) = 1.0;
  rhs.col(1).head(n) = -m_resVec;
  const gsMatrix<T> sol = m_borderedSolver->solve(rhs);

  m_deltaUt   = sol.col(0).head(n);   // p_t
  m_deltaLt   = sol(n,0);             // q_t
  m_deltaUbar = sol.col(1).head(n);   // p
  m_deltaLbar = sol(n,1);             // q

  return ( m_deltaUt.allFinite() && m_deltaUbar.allFinite() &&
           math::isfinite(m_deltaLt) && math::isfinite(m_deltaLbar) );
}

template <class T>
T gsALMCrisfield<T>::_chartTransversality(const gsVector<T> & du, const T dl,
                                          const T A0) const
{
  // |<d,x'>_A| / (||d||_A ||x'||_A), the cosine Riks's (n^t.x') = 0 sets to zero, in the
  // constraint metric <(a,alpha),(b,beta)>_A = a.b + A0*alpha*beta -- the metric in which
  // the corrector's own quadratic is written, hence the one that makes this dimensionless
  // (a Euclidean cosine would add a displacement to a bare load factor; see (h) and the C1
  // precedent at gsALMBase::_testSingularPoint).
  //
  // The numerator is IDENTICALLY 1: the last row of B reads w.p_t + gamma*q_t = 1, which is
  // <d,(p_t,q_t)>_A for the direction d whose metric dual (d_u, A0*d_lambda) is that row.
  // And (p_t,q_t) = q_t*(delta_u_t,1) whenever K is invertible, so
  //   1/(||d||_A ||(p_t,q_t)||_A) = |b0| / (||d||_A ||(delta_u_t,1)||_A),
  // i.e. this IS the normalised b0 -- computed without ever forming the cancelling
  // difference, and without a second copy of the b0 expression anywhere.
  //
  // ⚠ NOT computeLambdasSimple()'s m_b0 (which the elimination assigns at the top of that
  // function): on this path that member is never written by it, and _computeLambdasBordered()
  // writes 2*(p_t.DeltaU + A0*q_t*DeltaLambda) -- which is the SAME last row applied to
  // (p_t,q_t), i.e. identically 1, carrying no transversality information at all.
  //
  // At A0 = 0 (Scaling = 0, the exposed configuration) the form is only positive
  // SEMI-definite, but Cauchy-Schwarz still holds, so the value stays in [0,1] and reduces
  // to the plain cosine between DeltaU and delta_u_t.
  const T nd = math::sqrt( du.dot(du) + A0*dl*dl );
  const T nx = math::sqrt( m_deltaUt.dot(m_deltaUt) + A0*m_deltaLt*m_deltaLt );
  return (T)1 / (nd*nx);
}

template <class T>
void gsALMCrisfield<T>::_chartSingularityProbe(const gsVector<T> & f, T & rank, T & modeCos)
{
  // ONE extra solve against the bordered factorization that is already in hand; K itself is
  // never factorized (that is what would throw 3 on the tangent this chart exists to
  // survive) and gsALMBase::m_V is never touched. Complexity: one triangular solve pair,
  // O(nnz(B)) -- and it only ever runs on the rare degenerate-chart branch.
  const index_t n = m_numDof;
  gsMatrix<T> rhs(n+1,1);
  rhs.setZero();

  // A GENERIC right-hand side, from the fixed-seed LCG of gsALMBase::_computeCriticalMode
  // (Numerical Recipes constants; unsigned arithmetic is exactly modular and the map to
  // [-1/2,1/2) uses only powers of two, so this is bit-reproducible run to run and machine
  // to machine). Generic and not e0/ones for the reason finding m8 established there: a
  // structured probe can be exactly orthogonal to the very mode it must excite.
  unsigned seed = 2463534242u;
  for (index_t i = 0; i != n; ++i)
  {
    seed = 1664525u*seed + 1013904223u;
    rhs(i,0) = (T)( (index_t)(seed >> 8) ) / (T)16777216 - (T)0.5;
  }
  const T vn = rhs.col(0).head(n).norm();
  if (vn > (T)0) rhs.col(0).head(n) /= vn;

  const gsMatrix<T> sol  = m_borderedSolver->solve(rhs);
  const T           zeta = sol(n,0);

  // ★ THE ONE STEP OF INVERSE ITERATION ON K ITSELF, without ever factorizing K. The top
  // block of B gives K z - f*zeta = v, so the raw z is NOT K^-1 v: it carries the bordering's
  // own zeta*delta_u_t contribution, and in a degenerate chart that contribution DOMINATES
  // (zeta ~ 1/b0), which is what makes a naive |Kz|/|z| collapse to the load-direction
  // estimate |f|/|delta_u_t|. Subtracting the multiple of the chart tangent (p_t,q_t) that
  // zeroes the last component removes it EXACTLY:
  //   [z*;0] = [z;zeta] - (zeta/q_t)*[p_t;q_t]  =>  K z* = v + f*zeta - (zeta/q_t)*(f q_t) = v,
  // using K p_t = f q_t (the top block of the first bordered solve). So z* = K^-1 v exactly,
  // in exact arithmetic; in floating point the cancellation costs a factor 1/cos of the chart
  // this runs through, which is why the caller re-borders BEFORE probing.
  if (!math::isfinite(m_deltaLt) || m_deltaLt == (T)0)
  {
    rank = 0; modeCos = 0;                        // no usable chart tangent: assume the worst
    return;
  }
  const gsVector<T> z  = sol.col(0).head(n) - (zeta/m_deltaLt)*m_deltaUt;
  const T           zn = z.norm();

  // ||K z*||/||z*|| = ||v||/||z*|| = 1/||z*|| (v is a unit vector) is an upper bound on
  // sigma_min(K), and the ratio to ||K||_F makes it the dimensionless "how rank deficient is
  // the tangent" of (h). It sees a near-null direction of K even when the load never excites
  // it -- the branch-point configuration f in range(K), where delta_u_t = K^-1 f stays
  // bounded and any estimate taken along the load direction is blind.
  const T kn = m_jacMat.norm();                   // Frobenius
  if (!(zn > (T)0) || !math::isfinite(zn) || !(kn > (T)0))
  {
    rank = 0; modeCos = 0;                        // degenerate probe: assume the worst
    return;
  }
  rank = ((T)1/zn)/kn;

  // phi = z*/||z*|| is that approximate critical mode; |phi.f|/||f|| is the phi^T P statistic
  // of Wriggers & Simo 1990 eqs. (9a)/(9b), i.e. exactly what _testSingularPoint thresholds.
  // ||f|| = 0 takes the same convention as there (branch point).
  const T fn = f.norm();
  modeCos = (fn > (T)0) ? math::abs(z.dot(f))/(zn*fn) : (T)0;
}

template <class T>
void gsALMCrisfield<T>::_genericBorderDirection(gsVector<T> & du, T & dl, const T scaleA,
                                                const T A0) const
{
  // Govaerts 2000, Prop. 3.2.1: "a generic choice of B, C, D will do ... could be generated
  // by a random number generator". Same fixed-seed LCG as above and as
  // gsALMBase::_computeCriticalMode, with a DIFFERENT seed so that the direction is
  // independent of the probe vector; re-seeded on entry, so the re-bordering is a pure
  // function of the state and a restart reproduces it bit for bit.
  const index_t n = m_numDof;
  du.resize(n);
  unsigned seed = 1013904223u;                    // arbitrary but FIXED, != the probe's
  for (index_t i = 0; i != n; ++i)
  {
    seed = 1664525u*seed + 1013904223u;
    du[i] = (T)( (index_t)(seed >> 8) ) / (T)16777216 - (T)0.5;
  }
  seed = 1664525u*seed + 1013904223u;
  dl = (T)( (index_t)(seed >> 8) ) / (T)16777216 - (T)0.5;

  // Match the metric norm of the bordering it replaces, so that B keeps its row scale: the
  // chart is invariant under a rescaling of the row (it only reparameterises t), but the
  // CONDITIONING of B is not, and Govaerts's optimal construction is precisely a statement
  // about the norms of the bordering columns.
  const T nd = math::sqrt( du.dot(du) + A0*dl*dl );
  if (nd > (T)0 && scaleA > (T)0)
  {
    du *= scaleA/nd;
    dl *= scaleA/nd;
  }
}

template <class T>
void gsALMCrisfield<T>::_computeLambdasBordered()
{
  const gsVector<T> & Fstep = this->stepForcing();
  const T A0 = math::pow(m_phi,2)* Fstep.dot(Fstep); // see Lam & Morley 1992, eq. (11)

  // (c) THE COEFFICIENT MAP -- computeLambdasSimple()'s six coefficients, TERM FOR TERM in
  // the same order, under (delta_u_t,1) -> (p_t,q_t) and (u_bar,0) -> (p,q). Note that a0
  // and c0 keep their meaning: a0 = |p_t|^2 + A0 q_t^2 is the constraint metric applied to
  // the tangent of the solution set, c0 is untouched because it involves no chart vector.
  m_a0 = m_deltaUt.dot(m_deltaUt) + A0*m_deltaLt*m_deltaLt;
  m_b0 = 2*( m_deltaUt.dot(m_DeltaU)    + A0*m_deltaLt*m_DeltaL   );
  m_b1 = 2*( m_deltaUt.dot(m_deltaUbar) + A0*m_deltaLt*m_deltaLbar );
  m_c0 = m_DeltaU.dot(m_DeltaU) + m_DeltaL*m_DeltaL * A0 - math::pow(m_arcLength,2);
  m_c1 = 2*( m_DeltaU.dot(m_deltaUbar)  + A0*m_DeltaL *m_deltaLbar );
  m_c2 = m_deltaUbar.dot(m_deltaUbar)   + A0*m_deltaLbar*m_deltaLbar;

  /// Calculate the coefficients of the polynomial -- the same three lines as the
  /// m_alpha1/m_alpha2/m_alpha3 block of computeLambdasSimple() (and computeLambdasEta())
  m_alpha1 = m_a0;
  m_alpha2 = m_b0 + m_eta*m_b1;
  m_alpha3 = m_c0 + m_eta*m_c1 + m_eta*m_eta*m_c2;

  // The ONLY division on this path, and the only way a NaN could reach _step()'s residual
  // test SILENTLY instead of through throw 1 like every other failure here. It is
  // unreachable while B is nonsingular and f != 0: B[p_t;q_t] = [0;1] then forces
  // (p_t,q_t) != 0, and p_t = 0 would give q_t f = 0 hence q_t = 0. It CAN vanish in the
  // doubly degenerate case A0 = 0 AND f = 0, where the elimination divides by zero in
  // exactly the same way (its a0 = |delta_u_t|^2 + A0 = 0) -- pre-existing, but guarded
  // rather than inherited.
  if (!(m_alpha1 != (T)0))
  {
    gsWarn<<"gsALMCrisfield: the bordered constraint quadratic is degenerate (a0 = 0, i.e. "
          <<"a vanishing forcing together with Scaling = 0), so delta_lambda is not "
          <<"determined. Reporting the step as not converged.\n";
    throw 1;
  }

  m_discriminant = math::pow(m_alpha2 ,2) - 4 * m_alpha1 * m_alpha3;

  T t0, t1;
  if (m_discriminant >= 0)
  {
    m_eta = 1.0;
    t0 = (-m_alpha2 + math::pow(m_discriminant,0.5))/(2*m_alpha1);
    t1 = (-m_alpha2 - math::pow(m_discriminant,0.5))/(2*m_alpha1);
  }
  else
  {
    m_note += "\tC";
    // INVARIANT (c): D_t(eta) = q_t^2 D_lambda(eta) as POLYNOMIALS in eta, so
    // computeLambdasModified()'s three coefficients all scale by q_t^2 and its eta roots,
    // its -alpha2/alpha1 comparisons and its four option branches are bit-identically the
    // same decisions. It therefore needs -- and gets -- no change.
    computeLambdasModified();
    if ((m_discriminant >= 0) && m_eta > 0.05)
    {
      // computeLambdasEta()'s double root (Zhou 1995), in t. Invariant for the same reason:
      // -alpha2/(2 alpha1) mapped through delta_lambda = eta*q + t*q_t is the same
      // delta_lambda the elimination's computeLambdasEta() returns. Written out here only
      // because m_deltaLs must hold delta_lambda and not t (see below).
      m_alpha1 = m_a0;
      m_alpha2 = m_b0 + m_eta*m_b1;
      m_alpha3 = m_c0 + m_eta*m_c1 + m_eta*m_eta*m_c2;
      t0 = t1 = (-m_alpha2)/(2*m_alpha1);
      if (m_verbose) {gsInfo<<"Modified Complex Root Solve (bordered)\n";}
    }
    else
    {
      // The Lam & Morley 1992 eqs 13-17 complex-root fallback is NOT re-derived in the
      // bordered chart -- see (g). It COULD be: computeLambdasComplex()'s trial point
      // DeltaU + u_bar is the chart point (eta,delta_lambda) = (1,0), i.e. t = -q/q_t, so
      // DeltaU_cr = m_DeltaU + p - (q/q_t) p_t, which is bounded because q and q_t vanish
      // together. It is left out because it is untested: MEASURED on the
      // step-19 trace of fixture F at Scaling = 0, ds = 0.05 -- the trace this fallback
      // exists to rescue -- the branch does not fire at ANY of the three load scales. All
      // 99 corrector iterations of all three runs report an empty note delta and a
      // discriminant in {0, 1.00786}, both >= 0, so control never leaves the simple branch.
      // Nothing is lost there; anywhere else this is a documented limitation.
      gsWarn<<"gsALMCrisfield: the bordered corrector reached the complex-root branch of "
            <<"Lam & Morley 1992, which is not implemented in the bordered chart. "
            <<"Reporting the step as not converged.\n";
      throw 1;
    }
  }

  // m_deltaLs must keep holding DELTA_LAMBDA, not t: it is a gsALMBase member with existing
  // readers (computeLambdaDET, and the probes in the unit tests). Map back through the
  // affine chart change delta_lambda = eta*q + t*q_t.
  //
  // ⚠ (d) THE INDEX ORDER FLIPS WHEN q_t < 0. Both quadratics have a positive leading
  // coefficient, so t0 is the larger ROOT IN t and m_deltaLs[0] the larger in delta_lambda
  // only while q_t > 0. Harmless here because the selection below compares values and never
  // indices -- fatal for anyone who assumes m_deltaLs[0] came from t0.
  m_deltaLs[0] = m_eta*m_deltaLbar + t0*m_deltaLt;
  m_deltaLs[1] = m_eta*m_deltaLbar + t1*m_deltaLt;

  // Root selection: computeLambdaDOT()'s decision, computed on BOUNDED quantities. Its
  // DOT_i = m_deltaLs[i]*(m_DeltaUold.delta_u_t + A0*m_DeltaLold) contains the unbounded
  // delta_u_t, but only the COMPARISON matters and
  //   DOT1 - DOT2 = (dL0 - dL1)*(DeltaUold.u_t + A0 DeltaLold)
  //               = (t0 - t1)*q_t*(DeltaUold.u_t + A0 DeltaLold)
  //               = (t0 - t1)*(DeltaUold.p_t + A0*q_t*DeltaLold),
  // using dL0 - dL1 = (t0-t1) q_t and q_t*u_t = p_t. Both factors below stay O(1) as the
  // tangent degenerates, and the SIGN -- hence the decision -- is the elimination's.
  const T sel = (t0 - t1)*( m_DeltaUold.dot(m_deltaUt) + A0*m_deltaLt*m_DeltaLold );
  if (sel < 0)                        // DOT1 < DOT2
  {
    m_deltaL = m_deltaLs[1];
    m_deltaU = m_eta*m_deltaUbar + t1*m_deltaUt;
  }
  else                                // DOT1 > DOT2; the exact tie takes root 0, as
                                      // computeLambdaDOT()'s own DOT1 == DOT2 branch does
  {
    m_deltaL = m_deltaLs[0];
    m_deltaU = m_eta*m_deltaUbar + t0*m_deltaUt;
  }
}

template <class T>
gsStatus gsALMCrisfield<T>::step()
{
  // "Off" -- today's shipped behaviour, unchanged in every respect: no retry, no snapshot, and
  // m_useBordered is NOT touched here (getOptions() already re-asserts it false on every
  // applyOptions() call), which is what lets the forceBordered white-box idiom of the unit
  // tests keep working. Return BEFORE anything else happens, so a diff reader can see that no
  // new code executes on this path.
  if (m_borderedMode == borderedmode::Off)
    return Base::step();

  const bool primary = (m_borderedMode == borderedmode::Primary);

  // Attempt 1: "Primary" starts on the bordered chart (the references' own practice -- both
  // AUTO and pde2path use the bordered (n+1) solve as THE corrector, see the class-doc
  // heading); "Fallback" starts on the elimination, exactly like today.
  if (primary)
    m_note += "\tbordered";
  m_useBordered = primary;
  const gsStatus st = Base::step();     // attempt 1
  m_useBordered = false;

  if (st == gsStatus::Success)
    return st;

  // W3 -- the widened retry gate. Was `st != gsStatus::NotConverged`; a failed factorization
  // at a near-singular CORRECTOR iterate (quasiNewtonIteration()) is exactly the case the
  // bordered chart exists for (audit B.6: "a failed factorization at a near-singular tangent
  // -- the case the chart exists for -- never reaches the retry"), and NotConverged alone
  // never reached it: predictor() factorizes K on BOTH charts regardless of BorderedMode
  // (unchanged, out of scope, see the class C4 note), so a tangent that is already singular
  // at the step's SEED state is still not rescuable here; the corrector-iterate case is, and
  // that is the case that throws 3 -> SolverError instead of 1 -> NotConverged.
  if (!(st == gsStatus::NotConverged || st == gsStatus::SolverError))
    return st;

  // I5 -- a failed retry must leave no trace.
  const gsVector<T> DeltaUold_s     = m_DeltaUold;
  const T           DeltaLold_s     = m_DeltaLold;
  const T           arcLength_s     = m_arcLength;
  const T           arcLengthPrev_s = this->m_arcLength_prev;
  const bool        stepTaken_s     = this->m_stepTaken;
  const index_t     stabilityPrev_s = this->m_stabilityPrev;
  const T           phi_s           = m_phi;
  const gsVector<T> U_s             = m_U;
  const T           L_s             = m_L;
  // gsALMBase::_step() calls computeStability(false) inside EVERY corrector iteration
  // (gsALMBase.hpp:454), so a failed attempt 2 really does overwrite these, and the restore
  // is observable, not decorative. All FOUR must be carried: m_stability is derived from
  // m_indicator in the SAME call (gsALMBase.hpp:875-880, m_stability = this->stability()), so
  // restoring the other three without it would leave the pair mutually inconsistent -- worse
  // than merely stale.
  const index_t     negatives_s     = m_negatives;
  const T           indicator_s     = m_indicator;
  const gsVector<T> stabilityVec_s  = this->m_stabilityVec;
  const index_t     stability_s     = this->m_stability;

  // Attempt 2: the OTHER chart.
  m_note += primary ? "\telimination" : "\tbordered";
  m_useBordered = !primary;
  const gsStatus st2 = Base::step();    // Base::step() catches every throw itself
  m_useBordered = false;

  if (st2 == gsStatus::Success)
    return st2;

  // A failed _step() commits nothing -- gsALMBase::_step() calls iterationFinish() only
  // inside its converged branch -- but assert it rather than assume it.
  GISMO_ASSERT((m_U-U_s).norm()==(T)0 && m_L==L_s,
               "gsALMCrisfield: the failed retry (attempt 2) moved the converged state.");
  m_DeltaUold            = DeltaUold_s;
  m_DeltaLold            = DeltaLold_s;
  m_arcLength            = arcLength_s;
  this->m_arcLength_prev = arcLengthPrev_s;
  this->m_stepTaken      = stepTaken_s;
  this->m_stabilityPrev  = stabilityPrev_s;
  m_phi                  = phi_s;
  m_negatives            = negatives_s;
  m_indicator            = indicator_s;
  this->m_stabilityVec   = stabilityVec_s;
  this->m_stability      = stability_s;

  // Never report a status WORSE than attempt 1's: clamp to attempt 1's OWN status `st`
  // rather than a literal gsStatus::NotConverged. Under "Primary" this is load-bearing in
  // a direction the Fallback path never needed -- attempt 2 (the elimination) CAN emit
  // SolverError, while attempt 1 (the bordered chart) never does (class doc (e)) -- so
  // clamping to `st` keeps that asymmetry invisible to the caller. The clamp is
  // self-contained: the gate above means control reaches this line only when attempt 1
  // returned a definitive NotConverged or SolverError, so downgrading attempt 2's status
  // to attempt 1's own is correct on its own terms, independent of how any particular
  // caller (gsALMExploration, gsAPALM, ...) reacts to a non-Success status.
  return (this->m_status = st);
}

// ------------------------------------------------------------------------------------------------------------
// ---------------------------------------Output functions-----------------------------------------------------
// ------------------------------------------------------------------------------------------------------------

template <class T>
void gsALMCrisfield<T>::initOutput()
{
  gsInfo<<"\t";
  gsInfo<<std::setw(4)<<std::left<<"It.";
  gsInfo<<std::setw(17)<<std::left<<"Res. F";
  gsInfo<<std::setw(17)<<std::left<<"|dU|/|Du|";
  gsInfo<<std::setw(17)<<std::left<<"dL/DL";
  gsInfo<<std::setw(17)<<std::left<<"|U|";
  gsInfo<<std::setw(17)<<std::left<<"L";
  gsInfo<<std::setw(17)<<std::left<<"|DU|";
  gsInfo<<std::setw(17)<<std::left<<"DL";
  gsInfo<<std::setw(17)<<std::left<<"|dU|";
  gsInfo<<std::setw(17)<<std::left<<"dL";
  // "ds", not "ds²": the column below prints distance(DeltaU,DeltaL), a SQUARE ROOT --
  // kept identical to gsALMConsistentCrisfield::initOutput()'s header.
  gsInfo<<std::setw(17)<<std::left<<"ds";
  gsInfo<<std::setw(17)<<std::left<<"|dU|²";
  gsInfo<<std::setw(17)<<std::left<<"dL²";
  gsInfo<<std::setw(17)<<std::left<<"Dmin";
  gsInfo<<std::setw(17)<<std::left<<"note";
  gsInfo<<"\n";

  m_note = "";
}

template <class T>
void gsALMCrisfield<T>::stepOutput()
{
  // if (!m_quasiNewton)
  // {
    computeStability(false);
  // }
  // else
  //   m_indicator = 0;

  // Output only; uses the same (frozen) constraint metric as distance()
  T A0 = math::pow(m_phi,2)*this->stepForcing().dot(this->stepForcing());

  gsInfo<<"\t";
  gsInfo<<std::setw(4)<<std::left<<m_numIterations;
  gsInfo<<std::setw(17)<<std::left<<m_residueF;
  gsInfo<<std::setw(17)<<std::left<<m_residueU;
  gsInfo<<std::setw(17)<<std::left<<m_residueL;
  gsInfo<<std::setw(17)<<std::left<<(m_U+m_DeltaU).norm();
  gsInfo<<std::setw(17)<<std::left<<(m_L + m_DeltaL);
  gsInfo<<std::setw(17)<<std::left<<m_DeltaU.norm();
  gsInfo<<std::setw(17)<<std::left<<m_DeltaL;
  gsInfo<<std::setw(17)<<std::left<<m_deltaU.norm();
  gsInfo<<std::setw(17)<<std::left<<m_deltaL;
  gsInfo<<std::setw(17)<<std::left<<this->distance(m_DeltaU,m_DeltaL);//math::pow(m_DeltaU.dot(m_DeltaU) + A0*math::pow(m_DeltaL,2.0),0.5);
  gsInfo<<std::setw(17)<<std::left<<math::pow(m_DeltaU.norm(),2.0);
  gsInfo<<std::setw(17)<<std::left<<A0*math::pow(m_DeltaL,2.0);
  gsInfo<<std::setw(17)<<std::left<<m_indicator <<std::left << " (" <<std::left<< m_negatives<<std::left << ")";
  gsInfo<<std::setw(17)<<std::left<<m_note;
  gsInfo<<"\n";

  m_note = "";
}

} // namespace gismo