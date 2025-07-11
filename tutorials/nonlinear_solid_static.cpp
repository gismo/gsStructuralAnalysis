/** @file gsKLShell/tutorials/linear_solid.cpp

    @brief Tutorial for assembling solid elasticity

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst
*/

#include <gismo.h>

//! [Includes]
#ifdef gsElasticity_ENABLED
#include <gsElasticity/src/gsElasticityAssembler.h>
#include <gsElasticity/src/gsWriteParaviewMultiPhysics.h>
#include <gsElasticity/src/gsGeoUtils.h>
#endif


#include <gsStructuralAnalysis/src/gsStaticSolvers/gsStaticNewton.h>
#include <gsStructuralAnalysis/src/gsStructuralAnalysisTools/gsStructuralAnalysisUtils.h>
//! [Includes]

using namespace gismo;

int main(int argc, char *argv[])
{
#ifdef gsElasticity_ENABLED
    //! [Parse command line]
    bool plot  = false;
    index_t numRefine  = 1;
    index_t numElevate = 1;

    gsCmdLine cmd("Linear shell tutorial.");
    cmd.addInt( "e", "degreeElevation","Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving",  numRefine );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Read geometry]
    // Initialize [ori]ginal and [def]ormed geometry
    gsMultiPatch<> ori;
    std::string fileName = "paraboloid_volume.xml";
    gsReadFile<>(fileName, ori);
    //! [Read geometry]

    //! [Initialize geometry]
    // p-refine
    if (numElevate!=0)
        ori.degreeElevate(numElevate);
    // h-refine (only in first two directions)
    for (int r =0; r < numRefine; ++r)
    {
        ori.uniformRefine(1,1,0);
        ori.uniformRefine(1,1,1);
    }

    // creating basis
    gsMultiBasis<> basis(ori);
    //! [Initialize geometry]

    //! [Set boundary conditions]
    // Define the boundary conditions object
    gsBoundaryConditions<> bc;
    // Set the geometry map for computation of the Dirichlet BCs
    bc.setGeoMap(ori);

    // Set the boundary conditions
    gsFunctionExpr<> surf_force("0","0","-1e5",3);
    for (index_t c=0; c!=3; c++)
    {
        bc.addCornerValue(boundary::southwestfront, 0.0, 0, c); // (corner,value, patch, unknown)
        bc.addCornerValue(boundary::southeastfront, 0.0, 0, c); // (corner,value, patch, unknown)
        bc.addCornerValue(boundary::northwestfront, 0.0, 0, c); // (corner,value, patch, unknown)
        bc.addCornerValue(boundary::northeastfront, 0.0, 0, c); // (corner,value, patch, unknown)
    }
    bc.addCondition(boundary::front, condition_type::neumann, &surf_force);
    bc.addCondition(boundary::front, condition_type::neumann, &surf_force);
    //! [Set boundary conditions]

    //! [Set surface force]
    // The surface force is defined in the physical space, i.e. 3D
    gsFunctionExpr<> body_force("0","0","0",3);
    //! [Set surface force]

    //! [Define the assembler]
    gsElasticityAssembler<real_t> assembler(ori,basis,bc,body_force);
    assembler.options().setReal("YoungsModulus",1.E9);
    assembler.options().setReal("PoissonsRatio",0.45);
    assembler.options().setInt("MaterialLaw",material_law::saint_venant_kirchhoff);
    assembler.options().setInt("DirichletValues",dirichlet::l2Projection);

    //! [Define nonlinear residual functions]
    std::vector<gsMatrix<> > fixedDofs = assembler.allFixedDofs();
    // Function for the Jacobian
    gsStructuralAnalysisOps<real_t>::Jacobian_t Jacobian = [&assembler,&fixedDofs](gsVector<real_t> const &x, gsSparseMatrix<real_t> & m)
    {
        assembler.assemble(x,fixedDofs);
        m = assembler.matrix();
        return true;
    };

    // Function for the Residual
    gsStructuralAnalysisOps<real_t>::Residual_t Residual = [&fixedDofs,&assembler](gsVector<real_t> const &x, gsVector<real_t> & v)
    {
        assembler.assemble(x,fixedDofs);
        v = assembler.rhs();
        return true;
    };
    //! [Define nonlinear residual functions]

    //! [Assemble linear part]
    assembler.assemble();
    gsSparseMatrix<> K = assembler.matrix();
    gsVector<> F = assembler.rhs();
    //! [Assemble linear part]

    //! [Set static solver]
    gsStaticNewton<real_t> solver(K,F,Jacobian,Residual);
    solver.options().setInt("verbose",1);
    solver.initialize();
    //! [Set static solver]

    //! [Solve nonlinear problem]
    gsInfo<<"Solving system with "<<assembler.numDofs()<<" DoFs\n";
    gsStatus status = solver.solve();
    GISMO_ASSERT(status==gsStatus::Success,"Solver failed");
    gsVector<> solVector = solver.solution();
    //! [Solve nonlinear problem]

    //! [Construct solution and deformed geometry]
    // constructing solution as an IGA function
    gsMultiPatch<> displ;
    assembler.constructSolution(solVector,assembler.allFixedDofs(),displ);
    gsMultiPatch<> def = ori;
    for (index_t p = 0; p < def.nPieces(); ++p)
        def.patch(p).coefs() += displ.patch(p).coefs();
    //! [Construct solution and deformed geometry]

    //! [Evaluate solution]
    // Evaluate the solution on a reference point (parametric coordinates)
    // The reference points are stored column-wise
    gsMatrix<> refPoint(3,1);
    refPoint<<0.5,0.5,1.0;
    // Compute the values
    gsVector<> physpoint = def.patch(0).eval(refPoint);
    gsVector<> refVals = displ.patch(0).eval(refPoint);
    // gsInfo << "Displacement at reference point: "<<numVal<<"\n";
    gsInfo  << "Displacement at reference point ("
            <<ori.patch(0).eval(refPoint).transpose()<<"): "
            <<": "<<refVals.transpose()<<"\n";
    //! [Evaluate solution]

    // ! [Export visualization in ParaView]
    if (plot)
    {
        // Plot the displacements on the deformed geometry
        gsField<> solField(def, displ);
        gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview<>( solField, "Deformation", 10000, true);

        // Plot the membrane Von-Mises stresses on the geometry
        gsPiecewiseFunction<> VMm;
        assembler.constructCauchyStresses(displ,VMm,stress_components::von_mises);
        gsField<> membraneStress(def,VMm, true);
        gsWriteParaview(membraneStress,"MembraneStress",10000);
    }
    // ! [Export visualization in ParaView]

    return EXIT_SUCCESS;
#else
    GISMO_ERROR("The tutorial needs to be compiled with gsElasticity enabled");
    return EXIT_FAILURE;
#endif

}// end main