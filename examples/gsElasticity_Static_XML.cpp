/** @file gsElasticity_Static_XML.cpp

    @brief Static simulations of a solid reading from an XML file

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#include <gismo.h>

#ifdef gsElasticity_ENABLED
#include <gsElasticity/gsGeoUtils.h>
#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>
#endif

#include <gsStructuralAnalysis/src/gsStructuralAnalysisTools/gsStructuralAnalysisTypes.h>

#include <gsStructuralAnalysis/src/gsStaticSolvers/gsStaticNewton.h>


using namespace gismo;

#ifdef gsElasticity_ENABLED
int main (int argc, char** argv)
{
    // Input options
    int numElevate  = 0;
    int numHref     = 0;
    bool plot       = false;
    bool nonlinear  = false;

    int testCase = 1;

    bool write = false;

    std::string wn("data.csv");
    std::string fn1,fn2,fn3;
    fn1 = "volumes/brick.xml";
    fn2 = "pde/elasticity_brick.xml";
    fn3 = "options/static_solver.xml";

    gsCmdLine cmd("Static analysis using gsElasticity.");
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
    // Function for the Jacobian
    gsStructuralAnalysisOps<real_t>::Jacobian_t Jacobian = [&time,&stopwatch,&assembler,&fixedDofs](gsVector<real_t> const &x, gsSparseMatrix<real_t> &m)
    {
      stopwatch.restart();
      assembler.assemble(x,fixedDofs);
      time += stopwatch.stop();

      m = assembler.matrix();
      // gsInfo<<"matrix = \n"<<m.toDense()<<"\n";
      return true;
    };
    // Function for the Residual
    gsStructuralAnalysisOps<real_t>::Residual_t Residual = [&time,&stopwatch,&assembler,&fixedDofs](gsVector<real_t> const &x, gsVector<real_t> & result)
    {
      stopwatch.restart();
      assembler.assemble(x,fixedDofs);
      result = assembler.rhs();
      time += stopwatch.stop();
      return true;
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
    gsStaticNewton<real_t> staticSolver(matrix,vector,Jacobian,Residual);
    staticSolver.setOptions(solverOptions);

    gsInfo << "Solving...\n";
    // Solve linear problem
    gsVector<> solVector;
    gsStatus status;
    if (!nonlinear)
        status = staticSolver.solveLinear(solVector);
    else
        status = staticSolver.solveNonlinear(solVector);
    GISMO_ENSURE(status==gsStatus::Success,"Newton solver failed");

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
      gsWriteParaviewMultiPhysics(fields,"solutionElasticity",1000,false);
    }

    return 1;
}
#else//gsElasticity_ENABLED
int main(int argc, char *argv[])
{
    gsWarn<<"G+Smo is not compiled with the gsElasticity module.";
    return EXIT_FAILURE;
}
#endif