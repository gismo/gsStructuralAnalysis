/** @file example_StaticSolvers.cpp

    @brief Static simulations of a shell

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

#include <gsStructuralAnalysis/src/gsStaticSolvers/gsStaticDR.h>
#include <gsStructuralAnalysis/src/gsStaticSolvers/gsStaticNewton.h>
#include <gsStructuralAnalysis/src/gsStaticSolvers/gsControlDisplacement.h>
using namespace gismo;

void writeToCSVfile(std::string name, gsMatrix<> matrix)
{
    std::ofstream file(name.c_str());
    for(int  i = 0; i < matrix.rows(); i++){
        for(int j = 0; j < matrix.cols(); j++){
           std::string str = std::to_string(matrix(i,j));
           if(j+1 == matrix.cols()){
               file<<str;
           }else{
               file<<str<<',';
           }
        }
        file<<'\n';
    }
}

// Choose among various shell examples, default = Thin Plate
#ifdef gsKLShell_ENABLED
int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot  = false;
    bool stress= false;
    index_t numRefine  = 1;
    index_t numElevate = 1;

    index_t maxIt = 1e3;

    bool nonlinear = false;
    int verbose = 0;
    std::string fn;

    index_t Compressibility = 0;
    index_t material = 0;
    index_t impl = 1; // 1= analytical, 2= generalized, 3= spectral

    real_t E_modulus = 1.0;
    real_t PoissonRatio = 0.0;
    real_t Density = 1.0;
    real_t thickness = 1.0;
    real_t Ratio = 7.0;

    real_t alpha = 1.0;
    real_t damping = 0.1;

    std::string assemberOptionsFile("options/solver_options.xml");

    real_t dt = 0.1;

    gsCmdLine cmd("Simple example showing the use of the static solvers.");
    cmd.addString( "f", "file", "Input XML file for assembler options", assemberOptionsFile );
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving",  numRefine );
    cmd.addInt( "N", "maxit", "maxit",  maxIt );

    cmd.addReal( "d", "dt", "dt",  dt );
    cmd.addReal( "a", "alpha", "alpha",  alpha );
    cmd.addReal( "c", "damping", "damping",  damping );
    cmd.addInt( "m", "Material", "Material law",  material );
    cmd.addSwitch("nl", "Solve nonlinear problem", nonlinear);
    cmd.addInt("v","verbose", "0: no; 1: iteration output; 2: Full matrix and vector output", verbose);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("stress", "Create a ParaView visualization file with the stresses", stress);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Read input file]
    gsFileData<> fd(assemberOptionsFile);
    gsOptionList opts;
    fd.getFirst<gsOptionList>(opts);

    gsMultiPatch<> mp;
    gsMultiPatch<> mp_def;

    E_modulus = 1;
    // thickness = 0.15;
    thickness = 1;
    if (!Compressibility)
      PoissonRatio = 0.499;
    else
      PoissonRatio = 0.45;

    E_modulus = 1;
    real_t bDim = 1;
    real_t aDim = 1;

    mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // degree
    mp.patch(0).coefs().col(0) *= aDim;
    mp.patch(0).coefs().col(1) *= bDim;
    mp.addAutoBoundaries();
    mp.embed(3);

    // p-refine
    if (numElevate!=0)
        mp.degreeElevate(numElevate);

    // h-refine
    for (int r =0; r < numRefine; ++r)
        mp.uniformRefine();

    mp_def = mp;
    gsWriteParaview<>( mp_def    , "mp", 1000, true);

    //! [Refinement]
    gsMultiBasis<> dbasis(mp);

    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
    gsInfo<<mp_def<<"\n";
    gsInfo << dbasis.basis(0)<<"\n";

    gsBoundaryConditions<> bc;
    bc.setGeoMap(mp);

    gsConstantFunction<> displ(0.0,3);

    gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();
    bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 2 - z
    bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 2 - z
    bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

    // bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
    // bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
    bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0 ,false,1);

    bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0 ,false,2);
    bc.addCondition(boundary::east, condition_type::dirichlet, &displ, 0 ,false,0);

    gsVector<> tmp(3);
    tmp<<0,0,0;

    //! [Refinement]
    gsConstantFunction<> force(tmp,3);
    gsFunctionExpr<> t(std::to_string(thickness), 3);
    gsFunctionExpr<> E(std::to_string(E_modulus),3);
    gsFunctionExpr<> nu(std::to_string(PoissonRatio),3);
    gsFunctionExpr<> rho(std::to_string(Density),3);
    gsConstantFunction<> ratio(Ratio,3);

    real_t mu = E_modulus / (2 * (1 + PoissonRatio));
    gsConstantFunction<> alpha1(1.3,3);
    gsConstantFunction<> mu1(6.3e5/4.225e5*mu,3);
    gsConstantFunction<> alpha2(5.0,3);
    gsConstantFunction<> mu2(0.012e5/4.225e5*mu,3);
    gsConstantFunction<> alpha3(-2.0,3);
    gsConstantFunction<> mu3(-0.1e5/4.225e5*mu,3);

    std::vector<gsFunctionSet<>*> parameters;
    if (material==0) // SvK & Composites
    {
        parameters.resize(2);
        parameters[0] = &E;
        parameters[1] = &nu;
    }
    else if (material==1 || material==2) // NH & NH_ext
    {
        parameters.resize(2);
        parameters[0] = &E;
        parameters[1] = &nu;
    }
    else if (material==3) // MR
    {
        parameters.resize(3);
        parameters[0] = &E;
        parameters[1] = &nu;
        parameters[2] = &ratio;
    }
    else if (material==4) // OG
    {
        parameters.resize(8);
        parameters[0] = &E;
        parameters[1] = &nu;
        parameters[2] = &mu1;
        parameters[3] = &alpha1;
        parameters[4] = &mu2;
        parameters[5] = &alpha2;
        parameters[6] = &mu3;
        parameters[7] = &alpha3;
    }

    gsMaterialMatrixBase<real_t>* materialMatrix;

    gsOptionList options;
    if      (material==0 && impl==1)
    {
        parameters.resize(2);
        options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
        options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
    }
    else
    {
        options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",material);
        options.addSwitch("Compressibility","Compressibility: (false): Imcompressible | (true): Compressible",Compressibility);
        options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",impl);
    }
    materialMatrix = getMaterialMatrix<3,real_t>(mp,t,parameters,rho,options);

    gsThinShellAssemblerBase<real_t>* assembler;
    assembler = new gsThinShellAssembler<3, real_t, true >(mp,dbasis,bc,force,materialMatrix);

    opts.addReal("WeakDirichlet","Penalty parameter weak dirichlet conditions",1e5);
    opts.addReal("WeakClamped","Penalty parameter weak clamped conditions",1e5);
    // Construct assembler object
    assembler->setOptions(opts);
    assembler->setPointLoads(pLoads);

    // Function for the Jacobian
    gsStructuralAnalysisOps<real_t>::Jacobian_t Jacobian = [&assembler,&mp_def](gsVector<real_t> const &x, gsSparseMatrix<real_t> & m)
    {
      ThinShellAssemblerStatus status;
      assembler->constructSolution(x,mp_def);
      status = assembler->assembleMatrix(mp_def);
      m = assembler->matrix();
      return status == ThinShellAssemblerStatus::Success;
    };

    gsStructuralAnalysisOps<real_t>::ALResidual_t ALResidual = [&displ,&bc,&assembler,&mp_def](gsVector<real_t> const &x, real_t lam, gsVector<real_t> & result)
    {
        ThinShellAssemblerStatus status;
        displ.setValue(lam,3);
        assembler->updateBCs(bc);
        assembler->constructSolution(x,mp_def);
        status = assembler->assembleVector(mp_def);
        result = assembler->rhs();
        return status == ThinShellAssemblerStatus::Success;
    };

    displ.setValue(1.0,3);
    assembler->updateBCs(bc);

    assembler->assemble();
    gsSparseMatrix<> K = assembler->matrix();
    gsVector<> F = assembler->rhs();
    assembler->assembleMass(true);
    gsVector<> M = assembler->rhs();


    gsStaticDR<real_t> DRM(M,F,ALResidual);
    gsOptionList DROptions = DRM.options();
    DROptions.setReal("damping",damping);
    DROptions.setReal("alpha",alpha);
    DROptions.setInt("maxIt",maxIt);
    DROptions.setReal("tol",1e-2);
    DROptions.setReal("tolE",1e-4);
    DROptions.setInt("verbose",verbose);
    DRM.setOptions(DROptions);


    gsStaticNewton<real_t> NWT(K,F,Jacobian,ALResidual);
    gsOptionList NWTOptions = NWT.options();
    NWTOptions.setInt("maxIt",maxIt);
    NWTOptions.setReal("tol",1e-6);
    NWTOptions.setInt("verbose",verbose);
    NWT.setOptions(NWTOptions);

    DRM.initialize();

    gsControlDisplacement<real_t> controlDR(&DRM);
    gsControlDisplacement<real_t> controlDC(&NWT);
    // control.step(0.5);
    // control.step(0.5);
    controlDR.step(1.0);

    NWT.reset();
    NWT.setUpdate(DRM.update());
    controlDC.step(1.0);

    gsVector<> displacements = controlDR.solutionU();

    mp_def = assembler->constructSolution(displacements);
    gsMultiPatch<> deformation = mp_def;
    for (size_t k = 0; k != mp_def.nPatches(); ++k)
        deformation.patch(k).coefs() -= mp.patch(k).coefs();

    // ! [Export visualization in ParaView]
    if (plot)
    {
        gsField<> solField(mp_def, deformation);
        gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview<>( solField, "solution", 1000, true);
        // ev.options().setSwitch("plot.elements", true);
        // ev.writeParaview( u_sol   , G, "solution");

        // gsFileManager::open("solution.pvd");

        gsInfo <<"Maximum deformation coef: "
               << deformation.patch(0).coefs().colwise().maxCoeff() <<".\n";
        gsInfo <<"Minimum deformation coef: "
               << deformation.patch(0).coefs().colwise().minCoeff() <<".\n";
    }

    delete assembler;
    delete materialMatrix;

    return EXIT_SUCCESS;

}// end main
#else//gsKLShell_ENABLED
int main(int argc, char *argv[])
{
    gsWarn<<"G+Smo is not compiled with the gsKLShell module.";
    return EXIT_FAILURE;
}
#endif
