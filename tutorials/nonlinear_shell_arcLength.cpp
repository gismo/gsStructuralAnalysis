/** @file gsStructuralAnalysis/tutorials/nonlinear_shell_arclength.cpp

    @brief Tutorial for assembling the Kirchhoff-Love shell

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst
*/

#include <gismo.h>

#ifdef gsKLShell_ENABLED
#include <gsKLShell/src/gsThinShellAssembler.h>
#include <gsKLShell/src/getMaterialMatrix.h>
#endif

#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMCrisfield.h>
#include <gsStructuralAnalysis/src/gsStructuralAnalysisTools/gsStructuralAnalysisUtils.h>


using namespace gismo;

// Choose among various shell examples, default = Thin Plate
int main(int argc, char *argv[])
{
#ifdef gsKLShell_ENABLED
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
    gsReadFile<>("surfaces/paraboloid.xml", ori);
    ori.computeTopology();
    //! [Read geometry]

    //! [Initialize geometry]
    // p-refine
    if (numElevate!=0)
        ori.degreeElevate(numElevate);
    // h-refine
    for (int r =0; r < numRefine; ++r)
        ori.uniformRefine();
    //! [Initialize geometry]

    //! [Construct basis]
    gsMultiBasis<> bases(ori);
    //! [Construct basis]

    //! [Set boundary conditions]
    // Define the boundary conditions object
    gsBoundaryConditions<> bc;
    // Set the geometry map for computation of the Dirichlet BCs
    bc.setGeoMap(ori);

    // Set the boundary conditions
    bc.addCornerValue(boundary::southwest, 0.0, 0, -1); // (corner,value, patch, unknown)
    bc.addCornerValue(boundary::southeast, 0.0, 0, -1); // (corner,value, patch, unknown)
    bc.addCornerValue(boundary::northwest, 0.0, 0, -1); // (corner,value, patch, unknown)
    bc.addCornerValue(boundary::northeast, 0.0, 0, -1); // (corner,value, patch, unknown)

    gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();
    gsVector<> point(2);
    gsVector<> load (3);
    point<< 0.5, 0.5 ; load << 0, 0.0, -1e4 ;
    pLoads.addLoad(point, load, 0 );
    //! [Set boundary conditions]

    //! [Set surface force]
    // The surface force is defined in the physical space, i.e. 3D
    gsFunctionExpr<> force("0","0","0",3);
    //! [Set surface force]


    //! [Define the material matrix class]
    // Define material parameters
    // The material parameters are defined in the physical domain as well (!)
    gsConstantFunction<> E (1.E9,3);
    gsConstantFunction<> nu(0.45,3);
    gsConstantFunction<> t (1E-2,3);

    // Define a linear material, see \ref gsMaterialMatrixLinear.h
    // The first parameter is the physical domain dimension

    // gsMaterialMatrixLinear<3,real_t> materialMatrix(ori,t,E,nu);

    std::vector<gsFunctionSet<>*> parameters{&E,&nu};
    gsOptionList options;
    options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",1);
    options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
    gsMaterialMatrixBase<real_t> * materialMatrix = getMaterialMatrix<3,real_t>(ori,t,parameters,options);
    //! [Define the material matrix class]

    //! [Define the assembler]
    // Define the assembler.
    // The first parameter is the physical domain dimension
    // The second parameter is the real type
    // The third parameter controls the bending stiffness
    gsThinShellAssembler<3, real_t, true > assembler(ori,bases,bc,force,materialMatrix);
    assembler.setPointLoads(pLoads);
    //! [Define the assembler]

    //! [Assemble linear part]
    assembler.assemble();
    gsSparseMatrix<> K = assembler.matrix();
    gsVector<> F = assembler.rhs();
    //! [Assemble linear part]

    //! [Define nonlinear jacobian]
    // Function for the Jacobian
    gsStructuralAnalysisOps<real_t>::Jacobian_t Jacobian = [&assembler](gsVector<real_t> const &x, gsSparseMatrix<real_t> & m)
    {
        gsMultiPatch<> def;
        assembler.constructSolution(x,def);
        ThinShellAssemblerStatus status = assembler.assembleMatrix(def);
        m = assembler.matrix();
        return status==ThinShellAssemblerStatus::Success;
    };
    //! [Define nonlinear jacobian]

    //! [Define nonlinear residual functions]
    // Function for the Residual
    gsStructuralAnalysisOps<real_t>::ALResidual_t ALResidual = [&assembler,&F](gsVector<real_t> const &x, real_t lam, gsVector<real_t> & v)
    {
        gsMultiPatch<> def;
        assembler.constructSolution(x,def);
        ThinShellAssemblerStatus status = assembler.assembleVector(def);
        v = F - lam * F - assembler.rhs(); // assembler rhs - force = Finternal
        return status == ThinShellAssemblerStatus::Success;
    };
    //! [Define nonlinear residual functions]

    //! [Set ALM solver]
    gsInfo<<"Solving system with "<<assembler.numDofs()<<" DoFs\n";
    gsALMCrisfield<real_t> solver(Jacobian,ALResidual,F);
    solver.options().setSwitch("Verbose",true);
    solver.setLength(0.01);
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

        displ = assembler.constructDisplacement(solVector);
        def = assembler.constructSolution(solVector);

        // Export the files for each load-step
        if (plot)
        {
            // Plot the displacements on the deformed geometry
            gsField<> solField(def, displ);
            std::string outputName = "Deformation" + std::to_string(k) + "_";
            gsWriteParaview<>( solField, outputName, 1000, true);
            collection.addPart(outputName + "0.vts",k);
        }
    }

    //! [Save the paraview collection]
    if (plot)
        collection.save();
    //! [Save the paraview collection]

    return EXIT_SUCCESS;
#else
    GISMO_ERROR("The tutorial needs to be compiled with gsKLShell enabled");
    return EXIT_FAILED;
#endif
}// end main