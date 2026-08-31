/** @file example_ShearExploration.cpp

    @brief Perturbation-free shear-wrinkling de-risk driver (Milestone M0.5).

    Shear-loads a thick, coarse-mesh Kirchhoff-Love shell sheet with the
    arc-length method (ALM), detects the first wrinkling bifurcation WITHOUT any
    geometric perturbation / imperfection, then SWITCHES ONTO and TRACES the
    wrinkled (post-buckling) branch -- still perturbation-free. This demonstrates
    that gsALMBase's singular-point detection (determinant method) + extended
    (bordered-Newton) singular-point solve + branch switch can follow ONE branch
    end-to-end, which is the prerequisite for the gsALMExploration orchestrator.

    That the followed branch is genuinely the WRINKLING mode is confirmed two ways:
    the critical mode solutionV() is z-dominated (out-of-plane), and the
    out-of-plane amplitude max|coefs.col(2)| GROWS monotonically along the branch
    while the flat fundamental path (z==0) would remain unstable.

    The shear load is applied as a Dirichlet (prescribed-displacement) BC on the
    north edge, scaled by the arc-length parameter lambda. Because gsALMBase
    re-evaluates the residual at the current lambda every iteration
    (gsALMBase::computeResidual evaluates m_residualFun(U+DeltaU, L+DeltaL)),
    the displacement-control residual (result = assembler->rhs()) is consistent
    with the ALM solver; the constant Force vector only sets predictor scaling
    and residual normalization, so it must be non-zero.

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

#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMBase.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMLoadControl.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMCrisfield.h>

using namespace gismo;

template <class T>
gsMultiPatch<T> Rectangle(T L, T B);

#ifdef gsKLShell_ENABLED
int main (int argc, char** argv)
{
    // ---------------------- Input options (defaults run in seconds) --------
    index_t numElevate = 2;      // -> cubic basis
    index_t numHref    = 3;      // coarse mesh (thick sheet needs no fine mesh)
    index_t step       = 24;     // number of arc-length steps (bifurcation ~step 12);
                                 // >= detection + postSteps so the post-switch break fires
    index_t maxit      = 25;     // max Newton iterations per step

    real_t aDim        = 2.0;    // sheet length (x, sheared)
    real_t bDim        = 1.0;    // sheet width  (y)
    // Genuinely thick sheet (beta = W/t = 10). A thin sheet (t <~ 0.05, beta >~ 20)
    // makes the tangent so ill-conditioned that the determinant bifurcation
    // indicator is numerical noise (spurious detection at the first step); the
    // thick sheet gives a clean, monotone determinant that crosses zero at the
    // physical shear-wrinkling load on this coarse mesh.
    real_t thickness   = 0.1;

    real_t E_modulus   = 70e3;
    real_t PoissonRatio= 0.3;
    real_t Density     = 2710e9;

    // Arc-length method options
    real_t dLb         = 1e-2;   // arc length = prescribed-shear increment / step
    real_t dLbSwitch   = 1e-2;   // arc length used AFTER the branch switch
    real_t tol         = 1e-6;
    real_t tolU        = 1e-6;
    real_t tolF        = 1e-3;
    // Branch-switch predictor scale (Perturbation): switchBranch() nudges
    // m_U += (|V|/tau)*V, so the nudge magnitude is ~1/tau. COUNTERINTUITIVELY for
    // this PERTURBATION-FREE symmetric problem a SMALL tau (LARGE nudge) is needed:
    // the flat post-buckling branch is an EXACT Newton root (z==0 is a genuine, if
    // unstable, equilibrium because there is no geometric imperfection), so a tiny
    // nudge (library default tau=1e3 -> nudge 1e-3, or tau=1e2 -> 1e-2) is pulled
    // straight back to the flat branch by the fixed-lambda corrector (z-amp stays
    // exactly 0). tau=10 (nudge 0.1 ~ thickness) lands in the WRINKLED branch's
    // basin; tau=1 gives the identical wrinkled equilibrium. So the task's
    // 1e2-1e4 sweep is INVERTED here -- see the report for the full matrix.
    real_t tau         = 1e1;
    // Singular-point (extended bordered-Newton) tolerances. The extended solve
    // seeded straight from the pre-crossing point DIVERGES to a NaN shell metric on
    // this coarse thick mesh; the bisection first stage (spTolB != 0) brackets the
    // singular point first and hands the extended solve a seed inside its Newton
    // radius -> convergence. spTolE=1e-6 is reachable (1e-10 is near round-off).
    real_t spTolE      = 1e-6;   // extended-system tolerance (SingularPointComputeTolE)
    real_t spTolB      = 1e-4;   // bisection first-stage tol (0 => disabled) -- REQUIRED here
    real_t spTolT      = -1.0;   // sentinel: unset (negative) => library SingularPointTestTol
                                 // default governs. The test below the option block is `>= 0`,
                                 // not `> 0`, because 0 is itself a meaningful, valid value for
                                 // this option (see gsALMBase.hpp's _testSingularPoint
                                 // unresolved-band guard) and must not be swallowed by the sentinel.
    index_t postSteps  = 10;     // post-switch steps to trace (needs -N >= detection+this
                                 // else the outer step budget ends the trace first)
    index_t bifMethod  = 0;      // 0: determinant, 1: eigenvalue
    index_t method     = 0;      // ALM: 0: LoadControl, 1: Crisfield

    bool switchBranch  = true;   // switch onto + trace the wrinkled branch (default ON)
    bool noSwitch      = false;  // --noSwitch => detection only
    bool plot          = false;

    gsCmdLine cmd("Perturbation-free shear-wrinkling ALM de-risk driver.");
    cmd.addInt ("r","hRefine", "Number of uniform h-refinement steps", numHref);
    cmd.addInt ("e","degreeElevation", "Number of degree elevation steps", numElevate);
    cmd.addInt ("N","maxsteps", "Maximum number of arc-length steps", step);
    cmd.addInt ("I","maxit", "Maximum number of Newton iterations per step", maxit);
    cmd.addInt ("B","bifMethod", "Bifurcation method: 0: determinant, 1: eigenvalue", bifMethod);
    cmd.addInt ("m","method", "ALM method: 0: LoadControl, 1: Crisfield", method);
    cmd.addReal("a","aDim", "Sheet length (x)", aDim);
    cmd.addReal("b","bDim", "Sheet width (y)", bDim);
    cmd.addReal("T","thickness", "Sheet thickness", thickness);
    cmd.addReal("L","dLb", "Initial arc length", dLb);
    cmd.addReal("l","dLbSwitch", "Arc length used after the branch switch", dLbSwitch);
    cmd.addReal("F","factor", "Branch-switch perturbation scale tau", tau);
    cmd.addReal("E","spTolE", "SingularPointComputeTolE (extended-system tol)", spTolE);
    cmd.addReal("D","spTolB", "SingularPointComputeTolB (bisection tol, 0=off)", spTolB);
    cmd.addReal("t","spTolT", "SingularPointTestTol (detection tol) (unset/negative: library SingularPointTestTol default governs)", spTolT);
    cmd.addInt ("P","postSteps", "Number of post-switch steps to trace", postSteps);
    cmd.addSwitch("noSwitch", "Disable branch switch (detection only)", noSwitch);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    switchBranch = !noSwitch;

    gsInfo<<"E = "<<E_modulus<<"; nu = "<<PoissonRatio<<"\n";
    gsInfo<<"L = "<<aDim<<"; W = "<<bDim<<"; t = "<<thickness
          <<"; beta = W/t = "<<bDim/thickness<<"\n";

    // ---------------------- Geometry (NO perturbation) ---------------------
    gsStopwatch clock;
    clock.restart();

    gsMultiPatch<> mp = Rectangle(aDim,bDim);

    for (index_t i = 0; i < numElevate; ++i)
      mp.patch(0).degreeElevate();
    for (index_t i = 0; i < numHref; ++i)
      mp.patch(0).uniformRefine();

    gsMultiBasis<> dbasis(mp);
    gsInfo<<"Basis (patch 0): "<< mp.patch(0).basis() << "\n";

    // ---------------------- Boundary conditions (pure shear) ---------------
    // South edge fully clamped; north edge sheared in x (driven by lambda),
    // fixed in y and z. Fixing y top and bottom sets up the diagonal tension
    // field that drives shear wrinkling. No geometric imperfection anywhere.
    gsConstantFunction<> displ(0.0,3);        // x-shear, set to lambda by ALResidual
    gsConstantFunction<> displ_const(0.0,3);  // north-edge y offset (0 => pure shear
                                              // from an equilibrium rest state)

    gsBoundaryConditions<> BCs;
    BCs.setGeoMap(mp);

    BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 0);
    BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1);
    BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2);

    BCs.addCondition(boundary::north, condition_type::dirichlet, &displ,       0, false, 0);
    BCs.addCondition(boundary::north, condition_type::dirichlet, &displ_const, 0, false, 1);
    BCs.addCondition(boundary::north, condition_type::dirichlet, 0,            0, false, 2);

    std::string dirname = "ShearExplorationResults";
    std::string output  = "solution";
    gsFileManager::mkdir(dirname);
    if (plot)
      gsWriteParaview(mp,dirname + "/" + "mp",1000,true);

    // ---------------------- Material (SvK, linear) -------------------------
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

    // ---------------------- ALM operators ----------------------------------
    gsStopwatch stopwatch;
    real_t time = 0.0;

    // homogenizeDirichlet() zeroes the (inhomogeneous) prescribed shear when
    // forming the tangent -- required for a consistent displacement-controlled
    // tangent (cf. example_ShearWrinkling's Jacobian).
    gsStructuralAnalysisOps<real_t>::Jacobian_t Jacobian =
      [&time,&stopwatch,&assembler,&mp_def](gsVector<real_t> const &x, gsSparseMatrix<real_t> & m)
    {
      stopwatch.restart();
      assembler->constructSolution(x,mp_def);
      ThinShellAssemblerStatus status = assembler->assembleMatrix(mp_def);
      m = assembler->matrix();
      time += stopwatch.stop();
      return status == ThinShellAssemblerStatus::Success;
    };
    // Dirichlet (displacement-control) residual: lambda drives the x-shear.
    gsStructuralAnalysisOps<real_t>::ALResidual_t ALResidual =
      [&time,&stopwatch,&displ,&BCs,&assembler,&mp_def](gsVector<real_t> const &x, real_t lambda, gsVector<real_t> & result)
    {
      stopwatch.restart();
      displ.setValue(lambda,3);
      assembler->updateBCs(BCs);
      assembler->constructSolution(x,mp_def);
      ThinShellAssemblerStatus status = assembler->assembleVector(mp_def);
      // gsALMBase uses the opposite residual sign convention from gsStaticNewton
      // (cf. gsThinShell_ArcLength: result = Force - lam*Force - rhs()).
      result = -assembler->rhs();
      time += stopwatch.stop();
      return status == ThinShellAssemblerStatus::Success;
    };

    // Reference forcing vector: assemble at UNIT prescribed shear so the ALM
    // predictor tangent deltaUt = K^-1 * Force equals the physical dU/dlambda
    // (top edge moved by 1). Assembling at a smaller value would rescale the
    // predictor and make the first lambda-step overshoot -> divergence.
    displ.setValue(1.0,3);
    assembler->updateBCs(BCs);
    assembler->assemble();
    // rhs() = residual = F_int - F_ext; at the undeformed unit-shear state this
    // equals -L (minus the Dirichlet lift per unit lambda). The ALM predictor
    // deltaUt = K^-1 * Force must point toward the sheared equilibrium, so use
    // Force = +L = -rhs().
    gsVector<> Force = -assembler->rhs();
    displ.setValue(0.0,3);
    assembler->updateBCs(BCs);

    // ---------------------- ALM construction -------------------------------
    clock.restart();
    gsALMBase<real_t> * arcLength;
    if (method==0)
      arcLength = new gsALMLoadControl<real_t>(Jacobian,ALResidual,Force);
    else
      arcLength = new gsALMCrisfield<real_t>(Jacobian,ALResidual,Force);

    arcLength->options().setString("Solver","SimplicialLDLT");
    arcLength->options().setInt   ("BifurcationMethod",bifMethod); // 0: determinant
    arcLength->options().setReal  ("Length",dLb);
    arcLength->options().setSwitch("AdaptiveLength",true);
    arcLength->options().setInt   ("AdaptiveIterations",5);
    arcLength->options().setReal  ("Perturbation",tau);
    if (method!=0) // Crisfield-specific options
    {
      arcLength->options().setInt ("AngleMethod",0);
      arcLength->options().setReal("Scaling",0.0);
    }
    arcLength->options().setReal  ("Tol",tol);
    arcLength->options().setReal  ("TolU",tolU);
    arcLength->options().setReal  ("TolF",tolF);
    arcLength->options().setInt   ("MaxIter",maxit);
    arcLength->options().setSwitch("Verbose",true);
    arcLength->options().setReal  ("SingularPointComputeTolB",spTolB);
    arcLength->options().setReal  ("SingularPointComputeTolE",spTolE);
    // Only override when the flag was actually given (local left
    // at its sentinel), so a bare invocation is governed by the library default set inside the
    // solver's own constructor -- see the sentinel note at the spTolT declaration.
    if (spTolT >= (real_t)0) arcLength->options().setReal("SingularPointTestTol",spTolT);

    arcLength->applyOptions();
    arcLength->initialize();
    gsInfo<<"ALM setup: "<<clock.stop()<<" s\n";

    // ---------------------- Stepping loop ----------------------------------
    gsParaviewCollection collection(dirname + "/" + output);
    gsMultiPatch<> deformation = mp;
    patchSide ps(0,boundary::north);

    real_t Lold = 0;
    gsVector<> Uold = Force; Uold.setZero();
    gsMatrix<> solVector;
    real_t indicator = 0.0;
    arcLength->setIndicator(indicator);
    bool bisected  = false;
    bool switched  = false;      // true once we have hopped onto the wrinkled branch
    real_t dLb0 = dLb;
    index_t nSingular = 0;
    index_t nPost = 0;           // number of post-switch steps traced

    // Out-of-plane (z) amplitude of the deformed sheet: undeformed z==0, so this is
    // just max|coefs.col(2)|. On the flat fundamental path it stays ~0; on the
    // wrinkled branch it GROWS -> the signature that we are following wrinkling.
    auto outOfPlaneAmp = [&mp](const gsMultiPatch<>& def) -> real_t
    {
      return (def.patch(0).coefs().col(2) - mp.patch(0).coefs().col(2))
                 .cwiseAbs().maxCoeff();
    };

    gsInfo<<"\n"
          <<std::setw(6) <<std::left<<"step"
          <<std::setw(16)<<std::left<<"lambda"
          <<std::setw(16)<<std::left<<"reaction"
          <<std::setw(16)<<std::left<<"|U|"
          <<std::setw(16)<<std::left<<"z-amp"
          <<std::setw(16)<<std::left<<"indicator"
          <<std::setw(10)<<std::left<<"stable"<<"\n";
    gsInfo<<std::string(96,'-')<<"\n";

    gsStopwatch loopclock; loopclock.restart();
    for (index_t k = 0; k < step; k++)
    {
      gsStatus status;
      try { status = arcLength->step(); }
      catch (...) { status = gsStatus::AssemblyError; } // NaN metric on a bad predictor

      if (status==gsStatus::NotConverged || status==gsStatus::AssemblyError)
      {
        dLb = dLb / 2.;
        GISMO_ENSURE(dLb/dLb0 > 1e-6, "Arc length is becoming vanishingly small.");
        arcLength->setLength(dLb);
        arcLength->setSolution(Uold,Lold);
        bisected = true;
        k -= 1;
        gsInfo<<"  (step did not converge; halving arc length to "<<dLb<<")\n";
        continue;
      }

      // Bifurcation detection (perturbation-free) + branch switch. Only switch ONCE:
      // this de-risk must follow a single branch end-to-end, not hop repeatedly.
      arcLength->computeStability(false);
      if (arcLength->stabilityChange() && !switched)
      {
        nSingular++;
        gsInfo<<">>> SINGULAR POINT DETECTED at step "<<k
              <<", lambda="<<arcLength->solutionL()<<"\n";
      }
      if (arcLength->stabilityChange() && switchBranch && !switched)
      {
        // Seed the extended (bordered-Newton) singular-point solve from the LAST
        // CONVERGED point BEFORE the crossing (Uold,Lold, still on the stable side)
        // -- not from the current post-crossing point, whose tangent is already
        // indefinite. This is far more likely inside the Newton convergence radius
        // (cf. gsThinShell_ArcLength.cpp, which also seeds from Uold,Lold).
        // switchBranch=true makes computeSingularPoint call switchBranch()
        // INTERNALLY once the extended solve converges -- do NOT call it again here
        // (that would double-nudge m_U).
        gsStatus spStatus = gsStatus::NotConverged;
        try {
          spStatus = arcLength->computeSingularPoint(Uold,Lold,
                                                     /*switchBranch=*/true,
                                                     /*jacobian=*/false,
                                                     /*testPoint=*/true);
        } catch (...) { spStatus = gsStatus::OtherError; }

        if (spStatus == gsStatus::Success)
        {
          switched = true;
          // Confirm the critical mode is out-of-plane (z-dominated) -> wrinkling,
          // not an in-plane mode. solutionV() is the mode vector m_V (free DOFs).
          // Map it with a ZERO Dirichlet boundary so the north-edge shear does not
          // inject a spurious in-plane component; the mode displacement is then
          // mp_def.coefs()-mp.coefs().
          {
            gsMultiPatch<> mp_mode = mp;
            displ.setValue(0.0,3); assembler->updateBCs(BCs);
            assembler->constructSolution(arcLength->solutionV(),mp_mode);
            gsMatrix<> Vd = mp_mode.patch(0).coefs() - mp.patch(0).coefs();
            real_t vz  = Vd.col(2).cwiseAbs().maxCoeff();
            real_t vxy = std::max(Vd.col(0).cwiseAbs().maxCoeff(),
                                  Vd.col(1).cwiseAbs().maxCoeff());
            gsInfo<<"    branch switch CONVERGED (status=Success). "
                  <<"critical mode |V_z|/|V_xy| = "<<vz/(vxy+1e-30)
                  <<(vz>vxy ? "  (z-dominated => wrinkling mode)\n"
                            : "  (in-plane dominated)\n");
          }
          gsInfo<<"    tracing wrinkled branch with tau="<<tau
                <<", dLbSwitch="<<dLbSwitch<<" ...\n";
          // Re-seed the stepper from the switched state and use the post-switch
          // arc length.
          Uold = arcLength->solutionU();
          Lold = arcLength->solutionL();
          arcLength->setSolution(Uold,Lold);
          dLb = dLb0 = dLbSwitch;
          arcLength->setLength(dLb);
          bisected = false;
          continue; // next step() advances ALONG the wrinkled branch
        }
        else
        {
          // Extended iterations did not converge: honest de-risk signal.
          gsInfo<<"    branch switch did NOT converge (status="
                <<(spStatus==gsStatus::NotConverged?"NotConverged":"error")
                <<"); detection succeeded. Try finer mesh / looser spTolE / "
                  "enable bisection (spTolB>0).\n";
          break;
        }
      }

      indicator = arcLength->indicator();
      solVector = arcLength->solutionU();
      Uold = solVector;
      Lold = arcLength->solutionL();
      assembler->constructSolution(solVector,mp_def);

      real_t reaction = assembler->boundaryForce(mp_def,ps)(0,0);
      real_t zamp     = outOfPlaneAmp(mp_def);

      gsInfo<<std::setw(6) <<std::left<<k
            <<std::setw(16)<<std::left<<Lold
            <<std::setw(16)<<std::left<<reaction
            <<std::setw(16)<<std::left<<solVector.norm()
            <<std::setw(16)<<std::left<<zamp
            <<std::setw(16)<<std::left<<indicator
            // Report stability from the determinant indicator sign -- the same
            // signal that drives stabilityChange() detection (isStable()/
            // m_stability is only refreshed by the pre-loop setIndicator call).
            <<std::setw(10)<<std::left<<(indicator > 0 ? "stable" : "UNSTABLE")
            <<(switched ? "  <-- wrinkled branch" : "")<<"\n";

      if (switched)
      {
        nPost++;
        if (nPost >= postSteps) { gsInfo<<"    ("<<nPost<<" post-switch steps traced)\n"; break; }
      }

      if (plot)
      {
        deformation = mp_def;
        deformation.patch(0).coefs() -= mp.patch(0).coefs();
        gsField<> solField(mp,deformation);
        std::string fileName = dirname + "/" + output + util::to_string(k);
        gsWriteParaview<>(solField, fileName, 1000);
        collection.addPart(output + util::to_string(k) + "0" + ".vts",k);
      }

      if (!bisected)
      {
        dLb = dLb0;
        arcLength->setLength(dLb);
      }
      bisected = false;
    }
    gsInfo<<std::string(96,'-')<<"\n";
    gsInfo<<"Total arc-length stepping time: "<<loopclock.stop()<<" s\n";
    gsInfo<<"Total elapsed assembly time:    "<<time<<" s\n";
    gsInfo<<"Singular points detected:       "<<nSingular<<"\n";
    gsInfo<<"Post-switch steps traced:       "<<nPost<<"\n";

    if (plot)
      collection.save();

    delete arcLength;
    delete assembler;

    return EXIT_SUCCESS;
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
