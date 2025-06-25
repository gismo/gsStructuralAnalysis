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
#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>
#include <gsElasticity/gsGeoUtils.h>
#endif


#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMCrisfield.h>
#include <gsStructuralAnalysis/src/gsStructuralAnalysisTools/gsStructuralAnalysisUtils.h>
//! [Includes]

using namespace gismo;

#ifdef gsElasticity_ENABLED
int main(int argc, char *argv[])
{
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
    gsFunctionExpr<> surf_force("0","0","-1e10",3);
    for (index_t c=0; c!=3; c++)
    {
        bc.addCornerValue(boundary::southwestfront, 0.0, 0, c); // (corner,value, patch, unknown)
        bc.addCornerValue(boundary::southeastfront, 0.0, 0, c); // (corner,value, patch, unknown)
        bc.addCornerValue(boundary::northwestfront, 0.0, 0, c); // (corner,value, patch, unknown)
        bc.addCornerValue(boundary::northeastfront, 0.0, 0, c); // (corner,value, patch, unknown)
    }
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

    //! [Assemble linear part]
    assembler.assemble();
    gsSparseMatrix<> K = assembler.matrix();
    gsVector<> F = assembler.rhs();
    //! [Assemble linear part]

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
    gsStructuralAnalysisOps<real_t>::ALResidual_t ALResidual = [&fixedDofs,&assembler,&F](gsVector<real_t> const &x, real_t lam, gsVector<real_t> & v)
    {
        assembler.assemble(x,fixedDofs);
        v = F - lam * F - assembler.rhs(); // assembler rhs - force = Finternal
        return true;
    };
    //! [Define nonlinear residual functions]

    //! [Set ALM solver]
    gsInfo<<"Solving system with "<<assembler.numDofs()<<" DoFs\n";
    gsALMCrisfield<real_t> solver(Jacobian,ALResidual,F);
    solver.options().setSwitch("Verbose",true);
    solver.setLength(0.1);
    solver.applyOptions();
    solver.initialize();
    //! [Set ALM solver]

    //! [Solve nonlinear problem]
    index_t step = 50;
    gsParaviewCollection collection("Deformation");
    gsVector<> solVector;
    gsMultiPatch<> displ, def;
    for (index_t k=0; k<step; k++)
    {
        gsInfo<<"Load step "<< k<<"\n";
        solver.step();
        solVector = solver.solutionU();

        // constructing solution as an IGA function
        gsMultiPatch<> displ;
        assembler.constructSolution(solVector,assembler.allFixedDofs(),displ);
        gsMultiPatch<> def = ori;
        for (index_t p = 0; p < def.nPieces(); ++p)
            def.patch(p).coefs() += displ.patch(p).coefs();

        // ! [Export visualization in ParaView]
        if (plot)
        {
            // Plot the displacements on the deformed geometry
            gsField<> solField(def, displ);
            std::string outputName = "Deformation" + std::to_string(k) + "_";
            gsWriteParaview<>( solField, outputName, 1000, true);
            collection.addPart(outputName + "0.vts",k);
        }
    }
    //! [Solve nonlinear problem]

    // ![Save the paraview collection]
    if (plot)
        collection.save();
    // ![Save the paraview collection]

    return EXIT_SUCCESS;
#else
    GISMO_ERROR("The tutorial needs to be compiled with gsElasticity enabled");
    return EXIT_FAILED;
#endif

}// end main