![GitHub commits since latest release](https://img.shields.io/github/commits-since/gismo/gsStructuralAnalysis/latest?color=008A00)
![GitHub commit activity](https://img.shields.io/github/commit-activity/m/gismo/gsStructuralAnalysis?color=008A00)


# gsStructuralAnalysis

Module for structural analysis with solids ([`gsElasticity`](https://github.com/gismo/gsElasticity/)) or Kirchhoff-Love shells ([`gsKLShell`](https://github.com/gismo/gsKLShell/src/)).

|CMake flags|```-DGISMO_OPTIONAL="<other submodules>;gsStructuralAnalysis;gsSpectra"```|
|--:|---|
|License|![GitHub License](https://img.shields.io/github/license/gismo/gismo?color=008A00)|
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

