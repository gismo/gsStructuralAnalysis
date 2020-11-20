# gsStructuralAnalysis

Module for structural analysis with solids ([`gsElasticity`](https://github.com/gismo/gsElasticity/)) or Kirchhoff-Love shells ([`gsKLShell`](https://github.com/gismo/gsKLShell/)).

|CMake flags|```-DGISMO_STRUCTURALANALYSIS=ON``` (default ```OFF```)|
|--:|---|
|License|[MPL 2.0](https://www.mozilla.org/en-US/MPL/2.0/)|
|OS support|Linux, Windows, macOS|
|Build status| [CDash](link) |
|Repository|[gismo/gismo](https://github.com/gismo/gismo)|
|Status|completed|
|Developer|[Hugo Verhelst](https://github.com/hverhelst)|
|Maintainer|[h.m.verhelst@tudelft.nl](mailto:h.m.verhelst@tudelft.nl)|
|Last checked|13-11-2020|

#### Dependencies
`gsSpectra` via `-cmake . -DGISMO_WITH_SPECTRA=ON`.

#### Installation
```
cd path/to/build/dir
cmake . -DGISMO_WITH_SPECTRA=ON -DGISMO_STRUCTURALANALYSIS=ON
make
```

***

#### Use of the `gsStructuralAnalysis` module
The `gsStructuralAnalysis` 	module provides the following analysis tools:
* `gsStaticAnalysis`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Requires (nonlinear) stiffness matrix and a right-hand side (residual for nonlinear). Simply solves Newton iterations.
* `gsModalSolver`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Solves the vibration problem to find eigenfrequencies and mode shapes given linear mass and stiffness matrices.
* `gsBucklingSolver`&nbsp;&nbsp;&nbsp;&nbsp;Solves the a buckling eigenvalue problem given a solution **u** from a linear analysis, the linear stiffness matrix and the jacobian given **u**.
* `gsArcLengthIterator`&nbsp;&nbsp;Used for nonlinear buckling analysis (i.e. *post buckling analysis*). It includes arc-length schemes, extended arc-length methods and branch-switching methods.
* `gsTimeIntegrator`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Solves the (nonlinear) second-order structural dynamics problem.

All the tools in the `gsStructuralAnalysis` structural mass matrices, (linear/nonlinear) siffness matrices and forcing vectors/jacobians. The nonlinear modules typically work with jacobians and residuals of the following form (example using `gsThinShellAssembler`):
* Jacobian with solution **u**; K(**u**):
```
std::function<gsSparseMatrix<real_t> (gsVector<real_t> const &)> Jacobian;
Jacobian = [&assembler,&mp_def](gsVector<real_t> const &x)
{
  assembler.constructSolution(x,mp_def);
  assembler.assembleMatrix(mp_def);
  return assembler.matrix();
};
```
* Residual with solution **u**; R(**u**):
```
std::function<gsVector<real_t> (gsVector<real_t> const &, real_t, gsVector<real_t> const &) > Residual;
Residual = [&assembler,&mp_def](gsVector<real_t> const &x)
{
  assembler.constructSolution(x,mp_def);
  assembler.assembleVector(mp_def);
  return assembler.rhs();
};
```
* Arc-Length method residual with solution **u**, load factor lambda and linear forcing vector **F**; R(u,\lambda,**F**):
```
std::function<gsVector<real_t> (gsVector<real_t> const &, real_t, gsVector<real_t> const &) > ALResidual;
ALResidual = [&time,&stopwatch,&assembler,&mp_def](gsVector<real_t> const &x, real_t lam, gsVector<real_t> const &force)
{
  assembler.constructSolution(x,mp_def);
  assembler.assembleVector(mp_def);
  gsVector<real_t> Fint = -(assembler.rhs() - force);
  gsVector<real_t> result = Fint - lam * force;
  return result;
};

```

Where the `std::function` types are the ones accepted by the gsStructuralAnalysis module


#### Linear and nonlinear static analysis with gsStaticAnalysis
To use the `gsStaticAnalysis` class for a structural assembler (`gsElasticityAssembler` or `gsThinShellAssembler`), one simply performs the following steps:
```

```
#### Linear buckling analysis with gsBucklingSolver
To use the `gsStaticAnalysis` class for a structural assembler (`gsElasticityAssembler` or `gsThinShellAssembler`), one simply performs the following steps:
```

```

#### Post-Buckling analysis using gsArcLengthIterations
The implementation includes the *Riks Method*, the *(Consistent) Crisfield Method*, the *Extended Iterations Method* and a simple *Load Control Method*.

To use the `gsStaticAnalysis` class for a structural assembler (`gsElasticityAssembler` or `gsThinShellAssembler`), one simply performs the following steps:
```

```
#### Linear vibration analysis with gsModalAnalysis
To use the `gsStaticAnalysis` class for a structural assembler (`gsElasticityAssembler` or `gsThinShellAssembler`), one simply performs the following steps:
```

```

