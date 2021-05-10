/** @file gsThinShell_Static_XML.cpp

    @brief Static simulations of a shell given a XML file

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#include <gismo.h>

#include <gsKLShell/gsThinShellAssembler.h>
#include <gsKLShell/getMaterialMatrix.h>

#include <gsElasticity/gsWriteParaviewMultiPhysics.h>

#include <gsStructuralAnalysis/gsStaticSolver.h>

//#include <gsThinShell/gsNewtonIterator.h>

using namespace gismo;

// Choose among various shell examples, default = Thin Plate
int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot  = false;
    bool stress= false;
    bool mesh= false;
    index_t numRefine  = 1;
    index_t numElevate = 1;

    index_t Compressibility = 0;
    index_t material = 0;
    bool composite = false;
    index_t impl = 1; // 1= analytical, 2= generalized, 3= spectral

    index_t testCase = -1;
    bool nonlinear = false;
    bool verbose = false;
    std::string fn1,fn2,fn3,fn4;
    fn1 = "planar/unitplate.xml";
    fn2 = "pde/kirchhoff_shell1.xml";
    fn3 = "options/solver_options.xml";
    bool membrane = false;

    gsCmdLine cmd("Static analysis for thin shells.");
    cmd.addString( "s", "file", "Input XML file for assembler options", fn3 );
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving",  numRefine );

    cmd.addInt( "M", "Material", "Material law",  material );
    cmd.addInt( "c", "Compressibility", "1: compressible, 0: incompressible",  Compressibility );
    cmd.addInt( "I", "Implementation", "Implementation: 1= analytical, 2= generalized, 3= spectral",  impl );
    cmd.addSwitch("composite", "Composite material", composite);

    cmd.addInt( "t", "testCase", "Define test case",  testCase );
    cmd.addString( "f", "GEOMfile", "Input XML Geometry file", fn1 );
    cmd.addString( "F", "PDEfile", "Input XML PDE file", fn2 );
    cmd.addString( "L", "LayupFile", "Layup file", fn4 );
    cmd.addSwitch("nl", "Solve nonlinear problem", nonlinear);
    cmd.addSwitch("verbose", "Full matrix and vector output", verbose);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("mesh", "Plot mesh", mesh);
    cmd.addSwitch("stress", "Create a ParaView visualization file with the stresses", stress);
    cmd.addSwitch("membrane", "Use membrane model (no bending)", membrane);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Read input file]
    gsMultiPatch<> mp, mp_def, deformation;
    gsFileData<> fd;

    if (testCase == 1)
    {
        fn1 = "planar/unitplate.xml";
        fn2 = "pde/kirchhoff_shell1.xml";
    }
    else if (testCase == 2)
    {
        fn1 = "surface/scordelis_lo_roof.xml";
        fn2 = "pde/kirchhoff_shell_scordelis.xml";
    }
    else if (testCase == 3)
    {
        fn1 = "surface/quarter_hemisphere.xml";
        fn2 = "pde/kirchhoff_shell_hemisphere.xml";
    }
    else if (testCase == 4)
    {
        fn1 = "surface/pinched_cylinder.xml";
        fn2 = "pde/kirchhoff_shell_pinchedCylinder.xml";
    }

    gsReadFile<>(fn1, mp);
    fd.read(fn2);
    gsInfo << "Loaded file "<< fd.lastPath() <<"\n";

    // Boundary conditions
    gsBoundaryConditions<> bc;
    fd.getId(20, bc); // id=2: boundary conditions
    bc.setGeoMap(mp);
    gsInfo<<"Boundary conditions:\n"<< bc <<"\n";

    // Loads
    gsFunctionExpr<> force, pressure;
    fd.getId(21, force); // id=1: source function
    gsInfo<<"Force function "<< force << "\n";
    // fd.getId(22, pressure); // id=1: source function ------- TO DO!
    // gsInfo<<"Pressure function "<< force << "\n";

    // Loads
    gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();
    gsMatrix<> points,loads;
    fd.getId(30,points);
    fd.getId(31,loads);
    for (index_t k =0; k!=points.cols(); k++)
        pLoads.addLoad(points.col(k), loads.col(k), 0 ); // in parametric domain!

    // Material properties
    gsFunctionExpr<> t,E,nu,rho;
    fd.getId(10,t);
    gsInfo<<"thickness =  "<< t << "\n";
    fd.getId(11,E);
    fd.getId(12,nu);
    fd.getId(13,rho);

    gsFunctionExpr<> ratio;
    if (material==3)
    {
        fd.getId(14,ratio);
    }
    gsFunctionExpr<> alpha1, mu1, alpha2, mu2;
    if (material==4)
    {
        fd.getId(15,alpha1);
        fd.getId(16,mu1);
        fd.getId(17,alpha1);
        fd.getId(18,mu1);
    }

    std::vector<gsFunction<>*> parameters;
    if (material==0 || material==1 || material==2)
    {
        parameters.resize(2);
        parameters[0] = &E;
        parameters[1] = &nu;
    }
    else if (material==3)
    {
        parameters.resize(3);
        parameters[0] = &E;
        parameters[1] = &nu;
        parameters[2] = &ratio;
    }
    else if (material==4)
    {
        parameters.resize(6);
        parameters[0] = &E;
        parameters[1] = &nu;
        parameters[2] = &mu1;
        parameters[3] = &alpha1;
        parameters[4] = &mu2;
        parameters[5] = &alpha2;
    }
    //! [Read input file]

    // p-refine
    if (numElevate!=0)
        mp.degreeElevate(numElevate);

    // h-refine
    for (int r =0; r < numRefine; ++r)
        mp.uniformRefine();

    gsOptionList solverOptions;
    fd.getId(90, solverOptions); // id=4: assembler options

    // set initial deformation to undeformed state
    mp_def = mp;

    // define basis
    gsMultiBasis<> dbasis(mp);

    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
    gsInfo<<mp_def<<"\n";
    gsInfo << dbasis.basis(0)<<"\n";

    if (plot)
        gsWriteParaview<>( mp_def    , "mp", 1000, mesh);


    // Set MaterialMatrix
    gsMaterialMatrixBase<real_t>* materialMatrix;

    gsOptionList options;

    std::vector<gsMatrix<> > Gmats;
    std::vector<gsFunctionSet<> * > Gs;
    std::vector<gsFunctionSet<> * > Ts;
    std::vector<gsFunctionSet<> * > Phis;

    std::vector<gsConstantFunction<> > Gfuns;
    std::vector<gsConstantFunction<> > Tfuns;
    std::vector<gsConstantFunction<> > Phifuns;


    gsMatrix<> E11,E22,G12,nu12,nu21,alpha,thick;
    if (composite)
    {
        GISMO_ENSURE(!fn4.empty(),"Layup file must be provided!");
        fd.read(fn4);
        GISMO_ENSURE(fd.count<gsMatrix<>>()==7,"Composites must have 7 parameters!");

        fd.getId(0,E11);
        fd.getId(1,E22);
        fd.getId(2,G12);
        fd.getId(3,nu12);
        fd.getId(4,nu21);
        fd.getId(5,alpha);
        fd.getId(6,thick);

        index_t nLayers = E11.rows();
        Gmats.resize(nLayers);
        Gs.resize(nLayers);
        Ts.resize(nLayers);
        Phis.resize(nLayers);

        Gfuns.resize(nLayers);
        Tfuns.resize(nLayers);
        Phifuns.resize(nLayers);

        for (index_t k=0; k!=nLayers; k++)
        {
            Gmats[k] = gsCompositeMatrix(E11(k,0),E22(k,0),G12(k,0),nu12(k,0),nu21(k,0));
            Gmats[k].resize(Gmats[k].rows()*Gmats[k].cols(),1);
            Gfuns[k] = gsConstantFunction<>(Gmats[k],3);

            Phifuns[k] = gsConstantFunction<>(alpha(k,0),3);
            Tfuns[k] = gsConstantFunction<>(thick(k,0),3);

            Gs[k] = &Gfuns[k];
            Ts[k] = &Tfuns[k];
            Phis[k] = &Phifuns[k];
        }


    }

    if      (material==0 && impl==1)
    {
        if (composite)
        {
            materialMatrix = new gsMaterialMatrixComposite<3,real_t>(mp,Ts,Gs,Phis);
        }
        else
        {
            parameters.resize(2);
            options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
            options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
            materialMatrix = getMaterialMatrix<3,real_t>(mp,t,parameters,rho,options);
        }
    }
    else
    {
        options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",material);
        options.addSwitch("Compressibility","Compressibility: (false): Imcompressible | (true): Compressible",Compressibility);
        options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",impl);
        materialMatrix = getMaterialMatrix<3,real_t>(mp,t,parameters,rho,options);
    }

    gsThinShellAssemblerBase<real_t>* assembler;
    if(membrane)
        assembler = new gsThinShellAssembler<3, real_t, false>(mp,dbasis,bc,force,materialMatrix);
    else
        assembler = new gsThinShellAssembler<3, real_t, true >(mp,dbasis,bc,force,materialMatrix);

    fd.read(fn3);
    gsOptionList opts;
    fd.getFirst<gsOptionList>(opts);


    // Construct assembler object
    assembler->setOptions(opts);
    assembler->setPointLoads(pLoads);

    // Define Matrices
    assembler->assemble();
    gsSparseMatrix<> matrix = assembler->matrix();
    gsVector<> vector = assembler->rhs();
    typedef std::function<gsSparseMatrix<real_t> (gsVector<real_t> const &)>    Jacobian_t;
    typedef std::function<gsVector<real_t> (gsVector<real_t> const &) >         Residual_t;
    // Function for the Jacobian
    Jacobian_t Jacobian = [&assembler,&mp_def](gsVector<real_t> const &x)
    {
      assembler->constructSolution(x,mp_def);
      assembler->assembleMatrix(mp_def);
      gsSparseMatrix<real_t> m = assembler->matrix();
      return m;
    };
    // Function for the Residual
    Residual_t Residual = [&assembler,&mp_def](gsVector<real_t> const &x)
    {
      assembler->constructSolution(x,mp_def);
      assembler->assembleVector(mp_def);
      return assembler->rhs();
    };

    // Configure Structural Analsysis module
    gsStaticSolver<real_t> staticSolver(matrix,vector,Jacobian,Residual);
    gsDebugVar(solverOptions);
    staticSolver.setOptions(solverOptions);

    // Solve linear problem
    gsVector<> solVector;
    solVector = staticSolver.solveLinear();
    if (nonlinear)
        solVector = staticSolver.solveNonlinear();

    mp_def = assembler->constructSolution(solVector);

    deformation = mp_def;
    for (size_t k = 0; k != mp_def.nPatches(); ++k)
        deformation.patch(k).coefs() -= mp.patch(k).coefs();

    // ! [Export visualization in ParaView]
    if (plot)
    {
        gsField<> solField(mp_def, deformation);
        gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview<>( solField, "solution", 1000, mesh);
        // ev.options().setSwitch("plot.elements", true);
        // ev.writeParaview( u_sol   , G, "solution");

        // gsFileManager::open("solution.pvd");

        gsInfo <<"Maximum deformation coef: "
               << deformation.patch(0).coefs().colwise().maxCoeff() <<".\n";
        gsInfo <<"Minimum deformation coef: "
               << deformation.patch(0).coefs().colwise().minCoeff() <<".\n";
    }
    if (stress)
    {

        gsPiecewiseFunction<> membraneStresses;
        assembler->constructStress(mp_def,membraneStresses,stress_type::membrane);
        gsField<> membraneStress(mp_def,membraneStresses, true);

        gsPiecewiseFunction<> flexuralStresses;
        assembler->constructStress(mp_def,flexuralStresses,stress_type::flexural);
        gsField<> flexuralStress(mp_def,flexuralStresses, true);

        gsPiecewiseFunction<> stretches;
        assembler->constructStress(mp_def,stretches,stress_type::principal_stretch);
        gsField<> Stretches(mp_def,stretches, true);

        // gsPiecewiseFunction<> membraneStresses_p;
        // assembler->constructStress(mp_def,membraneStresses_p,stress_type::principal_stress_membrane);
        // gsField<> membraneStress_p(mp_def,membraneStresses_p, true);

        // gsPiecewiseFunction<> flexuralStresses_p;
        // assembler->constructStress(mp_def,flexuralStresses_p,stress_type::principal_stress_flexural);
        // gsField<> flexuralStress_p(mp_def,flexuralStresses_p, true);

        gsPiecewiseFunction<> stretch1;
        assembler->constructStress(mp_def,stretch1,stress_type::principal_stretch_dir1);
        gsField<> stretchDir1(mp_def,stretch1, true);

        gsPiecewiseFunction<> stretch2;
        assembler->constructStress(mp_def,stretch2,stress_type::principal_stretch_dir2);
        gsField<> stretchDir2(mp_def,stretch2, true);

        gsPiecewiseFunction<> stretch3;
        assembler->constructStress(mp_def,stretch3,stress_type::principal_stretch_dir3);
        gsField<> stretchDir3(mp_def,stretch3, true);


        gsField<> solutionField(mp,deformation, true);



        #ifdef GISMO_ELASTICITY
        std::map<std::string,const gsField<> *> fields;
        fields["Deformation"] = &solutionField;
        fields["Membrane Stress"] = &membraneStress;
        fields["Flexural Stress"] = &flexuralStress;
        fields["Principal Stretch"] = &Stretches;
        // fields["Principal Membrane Stress"] = &membraneStress_p;
        // fields["Principal Flexural Stress"] = &flexuralStress_p;
        fields["Principal Direction 1"] = &stretchDir1;
        fields["Principal Direction 2"] = &stretchDir2;
        fields["Principal Direction 3"] = &stretchDir3;

        gsWriteParaviewMultiPhysics(fields,"stress",5000,true);
        #else
        gsWriteParaview(solutionField, "Deformation");
        gsWriteParaview(membraneStress, "MembraneStress");
        gsWriteParaview(flexuralStress, "FlexuralStress");
        gsWriteParaview(Stretches, "PrincipalStretch");
        gsWriteParaview(stretchDir1, "PrincipalDirection1");
        gsWriteParaview(stretchDir2, "PrincipalDirection2");
        gsWriteParaview(stretchDir3, "PrincipalDirection3");
        #endif
    }

    return EXIT_SUCCESS;

}// end main