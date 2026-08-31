![GitHub commits since latest release](https://img.shields.io/github/commits-since/gismo/gsStructuralAnalysis/latest?color=008A00)
![GitHub commit activity](https://img.shields.io/github/commit-activity/m/gismo/gsStructuralAnalysis?color=008A00)


# gsStructuralAnalysis

Module for structural analysis with solids ([`gsElasticity`](https://github.com/gismo/gsElasticity/)) or Kirchhoff-Love shells ([`gsKLShell`](https://github.com/gismo/gsKLShell/src/)).

|CMake flags|```-DGISMO_OPTIONAL="<other submodules>;gsStructuralAnalysis;gsSpectra"```|
|--:|---|
|License|![GitHub License](https://img.shields.io/github/license/gismo/gismo?color=008A00)|
|DOI|[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15167527.svg)](https://doi.org/10.5281/zenodo.15167527)|
|OS support|Linux, Windows, macOS|
|Build status| [![ci](https://github.com/gismo/gsStructuralAnalysis/actions/workflows/ci.yml/badge.svg)](https://github.com/gismo/gsStructuralAnalysis/actions/workflows/ci.yml) |
|Developers/maintainers| [![Static Badge](https://img.shields.io/badge/@hverhelst-008A00)](https://github.com/hverhelst) [![Static Badge](https://img.shields.io/badge/@Crazy--Rich--Meghan-008A00)](https://github.com/Crazy-Rich-Meghan)|

#### Dependencies
`gsSpectra` via `-cmake . -DGISMO_OPTIONAL="<other submodules>;gsSpectra"`. The use of `gsSpectra` is not required, but strongly adviced. 

#### Installation
```
cd path/to/build/dir
cmake . -DGISMO_OPTIONAL="<other submodules>;gsStructuralAnalysis;gsSpectra"
make
```

***

#### Use of the `gsStructuralAnalysis` module
The `gsStructuralAnalysis` 	module provides the following analysis tools:
* `gsStaticAnalysis`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Requires (nonlinear) stiffness matrix and a right-hand side (residual for nonlinear). Simply solves Newton iterations.
* `gsModalSolver`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Solves the vibration problem to find eigenfrequencies and mode shapes given linear mass and stiffness matrices.
* `gsBucklingSolver`&nbsp;&nbsp;&nbsp;&nbsp;Solves the a buckling eigenvalue problem given a solution **u** from a linear analysis, the linear stiffness matrix and the jacobian given **u**.
* `gsALMBase`&nbsp;&nbsp;Used for nonlinear buckling analysis (i.e. *post buckling analysis*). It includes arc-length schemes, extended arc-length methods and branch-switching methods.
* `gsAPALM`&nbsp;&nbsp;Parallel implementation of the arc-length method
* `gsALMExploration`&nbsp;&nbsp;Automates a full landscape (bifurcation-diagram) exploration on top of any `gsALMBase` solver: traces curves, detects and classifies singular points, queues and de-duplicates branch-switch jobs, and collects the result in a `gsALMLandscape`.
* `gsALMLandscape`&nbsp;&nbsp;Container for an explored landscape: curves of (solution, load-factor) points carrying stability and bifurcation flags plus parent-curve connectivity, with CSV, ParaView and (in `gsHDF5`-enabled builds) HDF5 writers. Pure container, no solver logic.
* `gsTimeIntegrator`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Solves the (nonlinear) second-order structural dynamics problem.

All the tools in the `gsStructuralAnalysis` structural mass matrices, (linear/nonlinear) siffness matrices and forcing vectors/jacobians. The nonlinear modules typically work with jacobians and residuals of the following form (example using `gsThinShellAssembler`):
* Jacobian with solution **u**; K(**u**):
```
gsStructuralAnalysisOps<real_t>::Jacobian_t Jacobian = [&assembler,&mp_def](gsVector<real_t> const &x, gsSparseMatrix<real_t> & m)
{
    ThinShellAssemblerStatus status;
    assembler->constructSolution(x,mp_def);
    status = assembler->assembleMatrix(mp_def);
    m = assembler->matrix();
    return status == ThinShellAssemblerStatus::Success;
};
```
* Residual with solution **u**; R(**u**):
```
// Function for the Residual
gsStructuralAnalysisOps<real_t>::Residual_t Residual = [&assembler,&mp_def](gsVector<real_t> const &x, gsVector<real_t> & result)
{
    ThinShellAssemblerStatus status;
    assembler->constructSolution(x,mp_def);
    status = assembler->assembleVector(mp_def);
    result = assembler->rhs();
    return status == ThinShellAssemblerStatus::Success;
};

```
* Arc-Length method residual with solution **u**, load factor lambda and linear forcing vector **F**; R(u,\lambda,**F**):
```
gsStructuralAnalysisOps<real_t>::Residual_t Residual = [&assembler,&mp_def](gsVector<real_t> const &x, real_t lambda, gsVector<real_t> & result)
{
  ThinShellAssemblerStatus status;
  assembler.constructSolution(x,mp_def);
  assembler.assembleVector(mp_def);
  gsVector<T> Fint = -(assembler.rhs() - force);
  gsVector<T> result = Fint - lam * force;
  return status == ThinShellAssemblerStatus::Success;
};

```

Where the `std::function` types are the ones accepted by the gsStructuralAnalysis module. See the `struct` `gsStructuralAnalysisOps` in the file `gsStructuralAnalysisTools/gsStructuralAnalysisTypes`


#### Linear and nonlinear static analysis with gsStaticAnalysis
To use the `gsStaticAnalysis` class for a structural assembler (`gsElasticityAssembler` or `gsThinShellAssembler`), one simply performs the steps below.
##### Initialization of nonlinear solver
```
gsSparseMatrix<T>   matrix = any_assembler.function_for_StiffnessMatrix();
gsVector<T>         vector = any_assembler.function_for_rhs();
gsStaticNewton<T>   staticSolver(matrix,vector);

```
##### Initialization of nonlinear solver
```
gsSparseMatrix<T>   matrix = any_assembler.function_for_StiffnessMatrix();
gsVector<T>         vector = any_assembler.function_for_rhs();
Jacobian_t<T>       Jacobian = { your_jacobian };
Residual_t<T>       Residual = { your_residual };
gsStaticNewton<T>   staticSolver(matrix,vector,Jacobian,Residual); // see above documentation for definitions of Jacobian_t and Residual_t
```
##### General use
```
// get options
gsOptionList solverOptions = staticSolver.options();
// change some options
solverOptions.setInt("Verbose",1);
solverOptions.setInt("MaxIterations",10);
solverOptions.setReal("Tolerance",1e-6);
// set options
staticSolver.setOptions(solverOptions);

gsVector<T> solVector = staticSolver.solveNonlinear();
```

#### Linear buckling analysis with gsBucklingSolver
To use the `gsBucklingSolver` class for a structural assembler (`gsElasticityAssembler` or `gsThinShellAssembler`), one simply performs the following steps:
```
Jacobian_t<T>       Jacobian = { your_jacobian };
Residual_t<T>       Residual = { your_residual };
gsBucklingSolver<T> buckling(K_L,rhs,K_NL);

// computation using Eigen
buckling.compute();
// computation using gsSpectra for 10 buckling modes using a shift
buckling.computeSparse(shift,10);
// get results
gsMatrix<T> values = buckling.values();
gsMatrix<T> vectors = buckling.vectors();
```

#### Post-Buckling analysis using arc-length methods
The implementation includes the *Riks Method*, the *(Consistent) Crisfield Method* and a simple *Load Control Method*.

To use the `gsALMBase` class (here the derived `gsALMCrisfield`) for a structural assembler (`gsElasticityAssembler` or `gsThinShellAssembler`), one simply performs the following steps:
```
gsVector<T>         vector = any_assembler.function_for_rhs(); // this is the force of the linear system
Jacobian_t<T>       Jacobian = { your_jacobian };
ALResidual_t<T>     ALResidual = { your_arclenght_residual };
gsALMCrisfield<T> arclength(Jacobian, ALResidual, Force);

// example for setting options
arcLength.options().setInt("Method",method); // method 0: 1: 2: 3: 4:
arcLength.setLength(dL); // set arclength

arcLength.applyOptions();
arcLength.initialize();

for (index_t k=0; k<step; k++)
{
  gsInfo<<"Load step "<< k<<"\n";
  arcLength.step();
  arcLength.computeStability(quasiNewton);
  if (arcLength.stabilityChange())
  {
    gsInfo<<"Bifurcation spotted!"<<"\n";
    arcLength.computeSingularPoint(false);
    arcLength.switchBranch();
  }

  gsVector<T> solVector = arcLength.solutionU();
  T           LoadFactor = arcLength.solutionL();
}
```

##### Automating the loop with gsALMExploration
The loop above — step, stability check, singular-point computation, branch switch — is exactly what the `gsALMExploration` orchestrator automates. It takes a caller-configured `gsALMBase` solver plus one converged seed state and works off a job queue until the landscape is complete or `MaxCurves` is reached:
```
gsALMCrisfield<T>      arcLength(Jacobian, ALResidual, Force);
gsALMExploration<T>    explorer(&arcLength);
explorer.options().setInt   ("MaxCurves",8);
explorer.options().setInt   ("MaxPointsPerCurve",50);
explorer.options().setReal  ("Length",1e-2);       // arc length of a seed curve
explorer.options().setReal  ("SwitchLength",1e-2); // arc length of a branch curve
explorer.options().setString("OutputPrefix","results/landscape");
// U0,L0 must be a CONVERGED equilibrium; use (0,0) for an undeformed rest state.
// An overload takes several seeds and grows ONE landscape from all of them.
gsStatus status = explorer.solve(U0,L0);
const gsALMLandscape<T> & landscape = explorer.landscape();
```
Points to keep in mind when configuring it:
* **One job produces one curve, traced in BOTH arc-length directions.** The two halves of a seed (or of a branch) are assembled into a single curve object, running from one far end through the start state to the other, so they can never be duplicates of each other.
* **`MaxPointsPerCurve` is a per-sweep budget.** A curve traced in both directions may hold up to twice that many accepted points, plus the singular points stored on it, which the budget never counted.
* **Three independent comparison tolerances, with three different scopes.** `DedupTol` (default `1e-4`) de-duplicates branch *jobs*; `MoveTol` (default `1e-6`) is the no-progress guard that rejects a stalled step; `RetraceTol` (default `8e-2`) decides whether a *sweep* is retracing an already-stored curve. `RetraceTol` is deliberately the loosest: it compares samplings of two curves that never coincide, and it is not compared directly but scaled by `sqrt(StartSteps*SwitchLength/Length)`, because a genuinely emanating branch departs from its parent only as the square root of the distance travelled. A retracing sweep is rewound; the curve itself survives on the other sweep.
* **Export.** A non-empty `OutputPrefix` makes the explorer rewrite `<prefix>.csv` — and, in `gsHDF5`-enabled builds, `<prefix>.h5` — after every completed curve, as a crash checkpoint. `gsALMLandscape` also offers `writeCsv`, `writeParaview` and `saveHDF5`/`loadHDF5` directly.
* **Post-processing.** `examples/plot_landscape.py` renders a landscape CSV as a lambda-vs-||u|| equilibrium diagram (stable/unstable line styles, bifurcation and fold markers, parent-child links). It needs NumPy and Matplotlib; its HDF5-based modes (`--merge-loci` and the thesis-comparison metrics, which reclassify points by their solution profile) additionally need `h5py` and the `.h5` sibling of the CSV.

Complete drivers: `examples/example_BratuExploration.cpp` and `examples/example_ModifiedBratuExploration.cpp`. `examples/example_ShearExploration.cpp` runs the manual loop above on a shear-loaded Kirchhoff-Love sheet, i.e. it shows what the orchestrator automates.

#### Linear vibration analysis with gsModalAnalysis
To use the `gsModalAnalysis` class for a structural assembler (`gsElasticityAssembler` or `gsThinShellAssembler`), one simply performs the following steps:
```
gsSparseMatrix<T>   stif = any_assembler.function_for_StiffnessMatrix();
gsSparseMatrix<T>   mass = any_assembler.function_for_MassMatrix();
gsBucklingSolver<T> modal(stif,mass);

// computation using Eigen
modal.compute();
// computation using gsSpectra for 10 buckling modes using a shift
modal.computeSparse(shift,10);
// get results
gsMatrix<T> values = modal.values();
gsMatrix<T> vectors = modal.vectors();
```


---

## TODO — test-coverage gaps (recorded 2026-08-03)

- **`gsAPALM` has no regression oracle.** Every CSV oracle used by this module's gates exercises
  `gsALMExploration` on the Bratu drivers. No gate runs any `gsAPALM` driver, serial or MPI, so for
  any change touching `gsAPALM` a byte-identical oracle set is a *neutrality* check only — it
  proves the explorer still works and says nothing about the code being changed. Task 49's
  termination-guard fix therefore rests entirely on its own unit tests.
  Deferred deliberately; worth a driver + pinned CSV when someone has the appetite.
  ⚠ **Two of the seven oracles were re-based on 2026-08-05** and an older md5 raises a false alarm
  against them: both `example_ModifiedBratuExploration` oracles (with and without
  `--forcingCallback`) went from 92 to 72 CSV lines when the branch-job de-duplication was
  reworked and the two arc-length directions of a seed started tracing into one curve — the
  removed rows are exactly one retrace of the constant branch. Gate on the current values only.
  Note also what an identical md5 does *not* prove: `gsALMLandscape::writeCsv` sets no stream
  precision, so a match is exact in the integer columns and in the row/curve structure, but only
  about six significant figures in the load factor and the solution norm.
- **MPI runtime paths are unexercised.** Task 49's MPI guards are compile-verified only
  (`-DGISMO_WITH_MPI`, guarded branches confirmed active). No test runs them.
- **No mechanical guard against line-number citations in comments.** Comments citing `file.hpp:NNN`
  drift on every edit and silently start pointing at unrelated code. `src/gsALMSolvers/` was cleaned
  on 2026-08-03. Re-counted on 2026-08-10 with
  `grep -rEo "[A-Za-z_]+\.(h|hpp):[0-9]+" examples unittests benchmarks`: **18 such citations remain,
  all of them in `unittests/`** (17 in `gsALMSolvers_test.cpp`, one in `gsALMTestProblems.h`);
  `examples/` and `benchmarks/` carry none, so the earlier "21, spread over three directories" is
  superseded. Checking every cited target on that date, only **three** still land on what the citing
  comment describes: the `gsInfo` definition in `gsCore/gsDebug.h`, the corrector sequence of
  `gsALMBase::_step()`, and the `gsALMLoadControl` constructor signature. The rest have drifted —
  some by a few lines, several into a different function altogether: the citation described as "the
  scalar the corrector divides by" now lands on the commented-out legacy `gsALMRiks<T>::iteration()`
  above the live one; the one described as the halve-and-retry recovery of `gsALMExploration` now
  lands on `retraceThreshold()`; two describing `computeStability` and `_computeCriticalMode` land
  inside `computeSingularPoint`; others land on a blank line, a bare brace or a warning string.
  Most were invalidated *by the very edits that removed the citations from `src/`* — which is the
  argument for a mechanical gate rather than vigilance. A test grepping `\.(h|hpp):[0-9]` over
  `src/`, `examples/`, `unittests/` and `benchmarks/` would close the class. Symbol-anchor instead.
  (The count is a snapshot: repairs to individual citations were in flight when it was taken.)
