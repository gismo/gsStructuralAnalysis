/** @file gsElasticity_Static.cpp

    @brief Static simulations of a solid

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#include <gismo.h>

#include <gsElasticity/gsGeoUtils.h>
#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>

#include <gsStructuralAnalysis/gsStaticSolver.h>


using namespace gismo;

template <class T>
gsMultiPatch<T> BrickDomain(int n, int m, int o, int p, int q ,int r, T L, T B, T H);
template <class T>
gsMultiPatch<T> BrickDomain(int n, int p, T L, T B, T H);

int main (int argc, char** argv)
{
    // Input options
    int numElevate  = 0;
    int numHref     = 0;
    bool plot       = false;
    bool nonlinear  = false;

    real_t E_modulus     = 1;
    real_t PoissonRatio = 0;

    int testCase = 0;

    bool write = false;

    std::string wn("data.csv");
    std::string fn1,fn2,fn3;
    fn1 = "volumes/brick.xml";
    fn2 = "pde/elasticity_brick.xml";
    fn3 = "options/static_solver.xml";

    gsCmdLine cmd("Thin shell plate example.");
    cmd.addInt("r","hRefine",
       "Number of dyadic h-refinement (bisection) steps to perform before solving",
       numHref);
    cmd.addInt("t", "testcase",
        "Test case: 0: clamped-clamped, 1: pinned-pinned, 2: clamped-free",
       testCase);
    cmd.addString( "f", "GEOMfile", "Input XML Geometry file", fn1 );
    cmd.addString( "F", "PDEfile", "Input XML PDE file", fn2 );
    cmd.addInt("e","degreeElevation",
      "Number of degree elevation steps to perform on the Geometry's basis before solving",
      numElevate);
    cmd.addSwitch("nl", "Nonlinear elasticity (otherwise linear)", nonlinear);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addSwitch("write", "Write convergence data to file", write);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //! [Read input file]
    gsMultiPatch<> mp, mp_def, deformation;
    gsFileData<> fd;

    if (testCase == 1)
    {
        fn1 = "volumes/brick.xml";
        fn2 = "pde/elasticity_brick.xml";
        fn3 = "options/static_solver.xml";
    }
    gsReadFile<>(fn1, mp);

    // define basis
    gsMultiBasis<> dbasis(mp);

    gsReadFile<>(fn1, mp);
    fd.read(fn2);
    gsInfo << "Loaded file "<< fd.lastPath() <<"\n";

    // Boundary conditions
    gsBoundaryConditions<> bc;
    // fd.getId(20, bc); // id=2: boundary conditions
    fd.getFirst<gsBoundaryConditions<>>(bc);
    bc.setGeoMap(mp);
    gsInfo<<"Boundary conditions:\n"<< bc <<"\n";

    // Loads
    gsFunctionExpr<> bodyForce, pressure;
    fd.getFirst<gsFunctionExpr<>>(bodyForce);
    // fd.getId(21, bodyForce); // id=1: source function
    gsInfo<<"Body force function "<< bodyForce << "\n";

    // Material properties
    gsOptionList assemblerOptions;
    fd.getFirst<gsOptionList>(assemblerOptions);

    for (index_t i = 0; i < numElevate; ++i)
    {
        dbasis.degreeElevate();
        dbasis.uniformRefine();
    }
    for (index_t i = 0; i < numHref; ++i)
        dbasis.uniformRefine();

    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
    gsInfo<<mp<<"\n";
    gsInfo << dbasis.basis(0)<<"\n";

    if (plot)
        gsWriteParaview<>( mp_def    , "mp", 1000, true);

    fd.read(fn3);
    gsOptionList solverOptions;
    fd.getFirst<gsOptionList>(solverOptions);


    std::string output = "solutionElasticity";

    // plot geometry
    if (plot)
      gsWriteParaview(mp,"mp",1000,true);



    // -----------------------------------------------------------------------------------------------
    // --------------------------------------Solve Static problem-------------------------------------
    // -----------------------------------------------------------------------------------------------

    gsElasticityAssembler<real_t> assembler(mp,dbasis,bc,bodyForce);
    assembler.options().update(assemblerOptions,gsOptionList::addIfUnknown);

    // assembler.options().setReal("YoungsModulus",E_modulus);
    // assembler.options().setReal("PoissonsRatio",PoissonRatio);
    // assembler.options().setInt("MaterialLaw",material_law::hooke);

    gsStopwatch stopwatch,stopwatch2;
    // Define Matrices
    stopwatch.restart();
    stopwatch2.restart();
    real_t time = 0.0;
    real_t totaltime = 0.0;

    std::vector<gsMatrix<> > fixedDofs = assembler.allFixedDofs();
    typedef std::function<gsSparseMatrix<real_t> (gsVector<real_t> const &)>    Jacobian_t;
    typedef std::function<gsVector<real_t> (gsVector<real_t> const &) >         Residual_t;
    // Function for the Jacobian
    Jacobian_t Jacobian = [&time,&stopwatch,&assembler,&fixedDofs](gsVector<real_t> const &x)
    {
      stopwatch.restart();
      assembler.assemble(x,fixedDofs);
      time += stopwatch.stop();

      gsSparseMatrix<real_t> m = assembler.matrix();
      // gsInfo<<"matrix = \n"<<m.toDense()<<"\n";
      return m;
    };
    // Function for the Residual
    Residual_t Residual = [&time,&stopwatch,&assembler,&fixedDofs](gsVector<real_t> const &x)
    {
      stopwatch.restart();
      assembler.assemble(x,fixedDofs);
      gsVector<real_t> result = assembler.rhs();
      time += stopwatch.stop();
      return result; // - lam * force;
    };

    index_t materialLaw = assembler.options().getInt("MaterialLaw");
    assembler.options().setInt("MaterialLaw",material_law::hooke);
    // Assemble linear system to obtain the force vector
    assembler.assemble();
    time += stopwatch.stop();
    gsVector<> Force = assembler.rhs();

    assembler.options().setInt("MaterialLaw",materialLaw);

    //=============================================//
                  // Solving //
    //=============================================//

    gsSparseMatrix<> matrix = assembler.matrix();
    gsVector<> vector = assembler.rhs();

    // Configure Structural Analsysis module
    gsStaticSolver<real_t> staticSolver(matrix,vector,Jacobian,Residual);
    staticSolver.setOptions(solverOptions);

    gsInfo << "Solving...\n";
    // Solve linear problem
    gsVector<> solVector;
    solVector = staticSolver.solveLinear();
    if (nonlinear)
        solVector = staticSolver.solveNonlinear();

    totaltime += stopwatch2.stop();

    if (plot)
    {
      // solution to the nonlinear problem as an isogeometric displacement field
      gsMultiPatch<> displacement;
      assembler.constructSolution(solVector,assembler.allFixedDofs(),displacement);
      gsPiecewiseFunction<> stresses;
      assembler.constructCauchyStresses(displacement,stresses,stress_components::von_mises);

      // constructing an IGA field (geometry + solution)
      gsField<> displacementField(assembler.patches(),displacement);
      gsField<> stressField(assembler.patches(),stresses,true);
      // creating a container to plot all fields to one Paraview file
      std::map<std::string,const gsField<> *> fields;
      fields["Displacement"] = &displacementField;
      fields["von Mises"] = &stressField;
      gsWriteParaviewMultiPhysics(fields,"solution",1000,true);
    }

    return 1;
}


template <class T>
gsMultiPatch<T> BrickDomain(int n, int p, T L, T B, T H)
{
  int q = p;
  int r = p;

  int m = n;
  int o = n;
  gsMultiPatch<T> mp = BrickDomain(n, m, o, p, q, r, L, B, H);
  return mp;
}

template <class T>
gsMultiPatch<T> BrickDomain(int n, int m, int o, int p, int q ,int r, T L, T B, T H)
{
  // -------------------------------------------------------------------------
  // --------------------------Make beam geometry-----------------------------
  // -------------------------------------------------------------------------
  int dim = 3; //physical dimension
  gsKnotVector<> kv0;
  kv0.initUniform(0,1,0,p+1,1);
  gsKnotVector<> kv1;
  kv1.initUniform(0,1,0,q+1,1);
  gsKnotVector<> kv2;
  kv2.initUniform(0,1,0,r+1,1);

  for(index_t i = 0; i< n; ++i)
      kv0.uniformRefine();
  for(index_t i = 0; i< m; ++i)
      kv1.uniformRefine();
  for(index_t i = 0; i< o; ++i)
      kv2.uniformRefine();

  // Make basis
  gsTensorBSplineBasis<3,T> basis(kv0,kv1,kv2);

  // Initiate coefficient matrix
  gsMatrix<> coefs(basis.size(),dim);
  // Number of control points needed per component
  size_t len0 = basis.component(0).size();
  size_t len1 = basis.component(1).size();
  size_t len2 = basis.component(2).size();
  // Uniformly distribute control points per component
  gsVector<> coefvec0(len0);
  coefvec0.setLinSpaced(len0,0.0,L);
  gsVector<> coefvec1(basis.component(1).size());
  coefvec1.setLinSpaced(len1,0.0,B);
  gsVector<> coefvec2(basis.component(2).size());
  coefvec2.setLinSpaced(len2,0.0,H);

  // Define a matrix with ones
  gsVector<> temp(len0);
  temp.setOnes();
  for (index_t l = 0; l < len2; l++)
    {
        for (index_t k = 0; k < len1; k++)
        {
            index_t offset = l*len0*len1;
            // First column contains x-coordinates (length)
            coefs.col(0).segment(k*len0+offset,len0) = coefvec0;
            // Second column contains y-coordinates (width)
            coefs.col(1).segment(k*len0+offset,len0) = temp*coefvec1.at(k);

            coefs.col(2).segment(k*len0+offset,len0) = temp*coefvec2.at(l);
        }
    }
  // gsInfo<<"\n"<<coefs<<"\n";
  // Create gsGeometry-derived object for the patch
  gsTensorBSpline<3,real_t> shape(basis,coefs);

  gsMultiPatch<T> mp;
  mp.addPatch(shape);
  mp.addAutoBoundaries();

  return mp;
}
