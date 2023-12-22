/** @file gsThinShell_Static_XML.cpp

    @brief Static simulations of a shell given a XML file

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#include <gismo.h>

#ifdef gsKLShell_ENABLED
#include <gsKLShell/src/gsThinShellAssembler.h>
#include <gsKLShell/src/getMaterialMatrix.h>
#endif

#ifdef gsElasticity_ENABLED
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>
#endif

#include <gsStructuralAnalysis/src/gsStaticSolvers/gsStaticNewton.h>

//#include <gsThinShell/gsNewtonIterator.h>

using namespace gismo;

// Choose among various shell examples, default = Thin Plate
#ifdef gsKLShell_ENABLED
int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot  = false;
    bool stress= false;
    bool mesh= false;
    index_t numRefine  = 0;
    index_t numElevate = 0;

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
    bool membrane = false;

    gsCmdLine cmd("Static analysis for thin shells.");
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

    bool NURBS = false; // flag needed for benchmarking with NURBS geometries, compared to v21.6
    if (testCase == 1)
    {
        fn1 = "planar/unitplate.xml";
        fn2 = "pde/kirchhoff_shell1.xml";
    }
    else if (testCase == 2)
    {
        fn1 = "surface/scordelis_lo_roof.xml";
        fn2 = "pde/kirchhoff_shell_scordelis.xml";
        NURBS = true;
    }
    else if (testCase == 3)
    {
        fn1 = "surface/quarter_hemisphere.xml";
        fn2 = "pde/kirchhoff_shell_hemisphere.xml";
        NURBS = true;
    }
    else if (testCase == 4)
    {
        fn1 = "surface/pinched_cylinder.xml";
        fn2 = "pde/kirchhoff_shell_pinchedCylinder.xml";
        NURBS = true;
    }

    gsReadFile<>(fn1, mp);
    fd.read(fn2);
    gsInfo << "Loaded file "<< fd.lastPath() <<"\n";

    mp.clearTopology();
    mp.computeTopology();

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
    gsMatrix<index_t> pid_ploads;
    if ( fd.hasId(30) )
        fd.getId(30,points);
    if ( fd.hasId(31) )
        fd.getId(31,loads);

    if ( fd.hasId(32) )
        fd.getId(32,pid_ploads);
    else
        pid_ploads = gsMatrix<index_t>::Zero(1,points.cols());

    for (index_t k =0; k!=points.cols(); k++)
        pLoads.addLoad(points.col(k), loads.col(k), pid_ploads.at(k) ); // in parametric domain!

    gsInfo<<pLoads;

    // Reference points
    gsMatrix<index_t> refPatches;
    gsMatrix<> refPoints, refValue; // todo: add refValue..
    gsInfo<<"Reading reference point locations from "<<fn2<<" (ID=50) ...";
    if ( fd.hasId(50) )
        fd.getId(50,refPoints);
    gsInfo<<"Finished\n";
    gsInfo<<"Reading reference patches from "<<fn2<<" (ID=51) ...";
    if ( fd.hasId(51) )
        fd.getId(51,refPatches);
    gsInfo<<"Finished\n";
    gsInfo<<"Reading reference values from "<<fn2<<" (ID=52) ...";
    if ( fd.hasId(52) )
        fd.getId(52,refValue);
    else
        refValue = gsMatrix<>::Zero(mp.geoDim(),refPoints.cols());
    gsInfo<<"Finished\n";
    GISMO_ENSURE(refPatches.cols()==refPoints.cols(),"Number of reference points and patches do not match");

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

    std::vector<gsFunctionSet<>*> parameters;
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

    gsOptionList solverOptions, assemblerOptions;
    if (fd.hasId(90))
        fd.getId(90, solverOptions);
    if (fd.hasId(91))
        fd.getId(91, assemblerOptions);

    // set initial deformation to undeformed state
    mp_def = mp;

    // define basis
    gsMultiBasis<> dbasis(mp,NURBS);

    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() << "\n";
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

    // Construct assembler object
    assembler->setOptions(assemblerOptions);
    assembler->setPointLoads(pLoads);

    // Define Matrices
    assembler->assemble();
    const gsSparseMatrix<> & matrix = assembler->matrix();
    const gsVector<> & vector = assembler->rhs();

    // Function for the Jacobian
    gsStructuralAnalysisOps<real_t>::Jacobian_t Jacobian = [&assembler,&mp_def](gsVector<real_t> const &x, gsSparseMatrix<real_t> & m)
    {
      ThinShellAssemblerStatus status;
      assembler->constructSolution(x,mp_def);
      status = assembler->assembleMatrix(mp_def);
      m = assembler->matrix();
      return status == ThinShellAssemblerStatus::Success;
    };

    // Function for the Residual
    gsStructuralAnalysisOps<real_t>::Residual_t Residual = [&assembler,&mp_def](gsVector<real_t> const &x, gsVector<real_t> & result)
    {
      ThinShellAssemblerStatus status;
      assembler->constructSolution(x,mp_def);
      status = assembler->assembleVector(mp_def);
      result = assembler->rhs();
      return status == ThinShellAssemblerStatus::Success;
    };

    gsInfo << "Assemble done ("<<matrix.rows()<<"x"<<matrix.cols()<<")\n";
    // Configure Structural Analsysis module
    gsStaticNewton<real_t> staticSolver(matrix,vector,Jacobian,Residual);
    staticSolver.setOptions(solverOptions);

    // Solve linear problem
    gsInfo << "Solving linear system..\n";
    gsVector<> solVector;
    gsStatus status;
    status = staticSolver.solveLinear(solVector);
    GISMO_ENSURE(status==gsStatus::Success,"Newton solver failed");
    gsInfo << "Solving done.\n";
    if (nonlinear)
    {
        status = staticSolver.solveNonlinear(solVector);
        GISMO_ENSURE(status==gsStatus::Success,"Newton solver failed");
    }

    mp_def = assembler->constructSolution(solVector);
    gsInfo << "Solution constructed.\n";

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
    }

    // gsInfo <<"Maximum deformation coef: "
    //        << deformation.patch(0).coefs().colwise().maxCoeff() <<".\n";
    // gsInfo <<"Minimum deformation coef: "
    //        << deformation.patch(0).coefs().colwise().minCoeff() <<".\n";

    if (refPoints.cols()!=0)
    {
        gsMatrix<> refs(1,mp.geoDim()*refPoints.cols());
        for (index_t p=0; p!=refPoints.cols(); p++)
            refs.block(0,p*mp.geoDim(),1,mp.geoDim()) = deformation.piece(refPatches(0,p)).eval(refPoints.col(p)).transpose();

        gsInfo<<"Computed values\n";
        for (index_t p=0; p!=refPoints.cols(); ++p)
            gsInfo<<"x"<<std::to_string(p)<<"\ty"<<std::to_string(p)<<"\tz"<<std::to_string(p)<<"\t";
        gsInfo<<"\n";
        for (index_t p=0; p!=refPoints.cols(); ++p)
            gsInfo<<refs(0,mp.geoDim()*p)<<"\t"<<refs(0,mp.geoDim()*p+1)<<"\t"<<refs(0,mp.geoDim()*p+2)<<"\t";
        gsInfo<<"\n";

        gsInfo<<"Reference values\n"; // provided as mp.geoDim() x points.cols() matrix
        for (index_t p=0; p!=refValue.cols(); ++p)
            gsInfo<<"x"<<std::to_string(p)<<"\ty"<<std::to_string(p)<<"\tz"<<std::to_string(p)<<"\t";
        gsInfo<<"\n";
        for (index_t p=0; p!=refValue.cols(); ++p)
            for (index_t d=0; d!=mp.geoDim(); d++)
                gsInfo<<refValue(d,p)<<"\t";
        gsInfo<<"\n";
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

    delete materialMatrix;
    delete assembler;

    return EXIT_SUCCESS;

}// end main
#else//gsKLShell_ENABLED
int main(int argc, char *argv[])
{
    gsWarn<<"G+Smo is not compiled with the gsKLShell module.";
    return EXIT_FAILURE;
}
#endif
