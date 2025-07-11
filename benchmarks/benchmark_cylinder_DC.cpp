
/** @file benchmark_cylinder_DC.cpp

    @brief Computes the wrinkling behaviour of a cylinder or annulus

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
#include <gsKLShell/src/gsMaterialMatrixTFT.h>
#include <gsKLShell/src/gsFunctionSum.h>
#endif

#ifdef gsUnstructuredSplines_ENABLED
#include <gsUnstructuredSplines/src/gsSmoothInterfaces.h>
#endif

#ifdef gsHLBFGS_ENABLED
#include <gsHLBFGS/gsHLBFGS.h>
#endif

#ifdef gsOptim_ENABLED
#include <gsOptim/gsOptim.h>
#endif

#include <gsStructuralAnalysis/src/gsStaticSolvers/gsStaticDR.h>
#include <gsStructuralAnalysis/src/gsStaticSolvers/gsStaticOpt.h>
#include <gsStructuralAnalysis/src/gsStaticSolvers/gsStaticNewton.h>
#include <gsStructuralAnalysis/src/gsStaticSolvers/gsStaticComposite.h>
#include <gsStructuralAnalysis/src/gsStaticSolvers/gsControlDisplacement.h>
using namespace gismo;

void writeToCSVfile(const std::string& name, const gsMatrix<>& matrix)
{
  std::ofstream file(name.c_str());
  for(int  i = 0; i < matrix.rows(); i++)
  {
    for(int j = 0; j < matrix.cols(); j++)
    {
       if(j+1 == matrix.cols())
       {
           file<<std::setprecision(10)<<matrix(i,j);
       }
       else
       {
           file<<std::setprecision(10)<<matrix(i,j)<<',';
       }
    }
    file<<'\n';
  }
}

#ifdef gsKLShell_ENABLED

int main (int argc, char** argv)
{
    // Input options
    int numElevate  = 2;
    int numHref     = 5;
    bool plot       = false;
    bool mesh = false;
    int step = 10;

    bool    TFT = false;
    bool membrane= false;
    real_t dL = 1;
    // Solvers
    // * Dynamic Relaxation
    index_t maxitDR = 1e5;
    real_t  alpha = 1.0;
    real_t  damping = 0.1;
    bool DR  = false;
    // * Optimizer
    index_t maxitOPT = 1000;
    bool OP  = false;
    // * Newton Raphson
    index_t maxitNR = 50;
    bool NR  = false;


    index_t testCase = 0;


    std::string wn("data.csv");

    std::string assemberOptionsFile("options/solver_options.xml");

    gsCmdLine cmd("Wrinkling analysis with thin shells.");
    // Input settings
    cmd.addString( "f", "file", "Input XML file for assembler options", assemberOptionsFile );

    // Mesh settings
    cmd.addInt("r","hRefine", "Number of dyadic h-refinement (bisection) steps to perform before solving", numHref);
    cmd.addInt("e","degreeElevation", "Number of degree elevation steps to perform on the Geometry's basis before solving", numElevate);

    // Material settings
    cmd.addSwitch("TFT", "Use TFT material matrix", TFT);
    cmd.addSwitch("membrane", "Use membrane material", membrane);

    // Test case
    cmd.addInt("t", "test", "Test case: (0): Annulus | (1): Cylinder", testCase);

    // Load/displ stepping settings
    cmd.addInt("N", "maxsteps", "Maximum number of steps", step);
    cmd.addReal( "L", "step", "load step size", dL);
    cmd.addReal( "a", "alpha", "alpha",  alpha );
    cmd.addReal( "c", "damping", "damping",  damping );
    cmd.addSwitch("DR", "Use Dynamic Relaxation", DR);
    cmd.addSwitch("OP", "Use Optimizer", OP);
    cmd.addSwitch("NR", "Use Newton Raphson", NR);

    // Output settings
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addSwitch("mesh", "Plot mesh?", mesh);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    std::string output = "solution";
    std::string dirname = "ArcLengthResults";

    gsFileData<> fd(assemberOptionsFile);
    gsOptionList opts;
    fd.getFirst<gsOptionList>(opts);

    real_t thickness = 0.05e-3;
    real_t YoungsModulus = 1e9;
    real_t PoissonRatio = 0.5;
    real_t Density = 1e0;

    gsMultiPatch<> mp, mp_def;
    if (testCase==0)
    {
        mp = gsNurbsCreator<>::MultiPatchAnnulus(0.0625,0.25);
        mp.embed(3);
    }
    else if (testCase==1)
        mp = gsNurbsCreator<>::MultiPatchCylinder(0.25,1.0);

    mp.computeTopology();

    for (size_t p = 0; p!=mp.nPatches(); ++p)
    {
        for(index_t i = 0; i< numElevate; ++i)
            mp.patch(p).degreeElevate();    // Elevate the degree
        for(index_t i = 0; i< numHref; ++i)
            mp.patch(p).uniformRefine();
    }

    mp_def = mp;

    // Make unstructured spline
    gsMultiPatch<> geom;
    gsMultiBasis<> dbasis;
    gsMappedBasis<2,real_t> bb2;
    gsSparseMatrix<> global2local;
    gsSmoothInterfaces<2,real_t> smoothInterfaces(mp);
    smoothInterfaces.options().setSwitch("SharpCorners",false);
    smoothInterfaces.compute();
    smoothInterfaces.matrix_into(global2local);

    global2local = global2local.transpose();
    geom = smoothInterfaces.exportToPatches();
    dbasis = smoothInterfaces.localBasis();
    bb2.init(dbasis,global2local);

    // Assign BCs
    gsFunctionExpr<> dx("sqrt(x^2+y^2)*(cos(atan2(y,x) + pi/2*t)- cos(atan2(y,x)))",3);
    gsFunctionExpr<> dy("sqrt(x^2+y^2)*(sin(atan2(y,x) + pi/2*t)- sin(atan2(y,x)))",3);
    real_t dz_val = 0;
    if (testCase==0)
        dz_val = 125e-3;
    else if (testCase==1)
        dz_val = 1.0;
    gsFunctionExpr<> dz(util::to_string(dz_val) + "*1",3);

    // Boundary conditions
    gsBoundaryConditions<> BCs;
    BCs.setGeoMap(geom);
    for (size_t p = 0; p!=geom.nPatches(); ++p)
    {
        BCs.addCondition(p,boundary::south, condition_type::dirichlet, 0  , 0, false, -1);
        BCs.addCondition(p,boundary::north, condition_type::dirichlet, &dx, 0, false, 0);
        BCs.addCondition(p,boundary::north, condition_type::dirichlet, &dy, 0, false, 1);
        BCs.addCondition(p,boundary::north, condition_type::dirichlet, &dz, 0, false, 2);
    }

    if (testCase==0)
        dirname = dirname + "/Annulus_-r" + std::to_string(numHref) + "-e" + std::to_string(numElevate);
    else if (testCase==1)
        dirname = dirname + "/Cylinder_-r" + std::to_string(numHref) + "-e" + std::to_string(numElevate);
    else
        GISMO_ERROR("Test case "<<testCase<<" not implemented");
    if (TFT)
        dirname = dirname + "_TFT";
    if (membrane)
      	dirname = dirname + "_membrane";

    gsFileManager::mkdir(dirname);

    output =  "solution";
    wn = output + "data.txt";

    // plot geometry
    if (plot)
      gsWriteParaview(mp,dirname + "/" + "mp",1000,true);

    // Linear isotropic material model
    gsFunctionExpr<> force("0","0","0",3);

    gsConstantFunction<> t(thickness,3);
    gsConstantFunction<> E(YoungsModulus,3);
    gsConstantFunction<> nu(PoissonRatio,3);
    gsConstantFunction<> rho(Density,3);

    std::vector<gsFunctionSet<>*> parameters;
    parameters.resize(2);
    parameters[0] = &E;
    parameters[1] = &nu;

    gsMaterialMatrixBase<real_t>::uPtr materialMatrix;
    gsMaterialMatrixBase<real_t>::uPtr materialMatrixTFT;
    gsThinShellAssemblerBase<real_t>* assembler;

    gsOptionList options;
    options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",1);
    options.addSwitch("Compressibility","Compressibility: (false): Imcompressible | (true): Compressible",false);
    options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
    materialMatrix = getMaterialMatrix<3,real_t>(mp,t,parameters,rho,options);
    gsMaterialMatrixContainer<real_t> materialMatrixContainer;
    if (TFT)
    {
        materialMatrixTFT = memory::make_unique(new gsMaterialMatrixTFT<3,real_t,true>(*materialMatrix));
        for (size_t p = 0; p!=mp.nPatches(); p++)
            materialMatrixContainer.add(*materialMatrixTFT);
        assembler = new gsThinShellAssembler<3, real_t, false>(geom,dbasis,BCs,force,materialMatrixContainer);
    }
    else if (membrane)
    {
        for (size_t p = 0; p!=mp.nPatches(); p++)
            materialMatrixContainer.add(*materialMatrix);
        assembler = new gsThinShellAssembler<3, real_t, false >(geom,dbasis,BCs,force,materialMatrixContainer);
    }
    else
    {
        for (size_t p = 0; p!=mp.nPatches(); p++)
            materialMatrixContainer.add(*materialMatrix);
        assembler = new gsThinShellAssembler<3, real_t, true >(geom,dbasis,BCs,force,materialMatrixContainer);
    }


    assembler->setOptions(opts);
    assembler->options().setInt("Continuity",-1);
    assembler->setSpaceBasis(bb2);

    // Assemble linear system to obtain the force vector
    assembler->assemble();
    gsVector<> F = assembler->rhs();
    gsSparseMatrix<> K = assembler->matrix();
    assembler->assembleMass(true);
    gsVector<> M = assembler->rhs();

    // Function for the Jacobian
    gsStructuralAnalysisOps<real_t>::Jacobian_t Jacobian = [&assembler,&bb2, &geom](gsVector<real_t> const &x, gsSparseMatrix<real_t> & m)
    {
        ThinShellAssemblerStatus status;
        gsMatrix<real_t> solFull = assembler->fullSolutionVector(x);
        size_t d = geom.targetDim();
        GISMO_ASSERT(solFull.rows() % d==0,"Rows of the solution vector does not match the number of control points");
        solFull.resize(solFull.rows()/d,d);
        gsMappedSpline<2,real_t> mspline(bb2,solFull);
        gsFunctionSum<real_t> def(&geom,&mspline);
        assembler->assembleMatrix(def);
        status = assembler->assembleMatrix(def);
        m = assembler->matrix();
        return status == ThinShellAssemblerStatus::Success;
    };

    // Function for the Residual
    gsStructuralAnalysisOps<real_t>::Residual_t Residual = [&assembler,&bb2, &geom](gsVector<real_t> const &x, gsVector<real_t> & result)
    {
        ThinShellAssemblerStatus status;
        gsMatrix<real_t> solFull = assembler->fullSolutionVector(x);
        size_t d = geom.targetDim();
        GISMO_ASSERT(solFull.rows() % d==0,"Rows of the solution vector does not match the number of control points");
        solFull.resize(solFull.rows()/d,d);

        gsMappedSpline<2,real_t> mspline(bb2,solFull);
        gsFunctionSum<real_t> def(&geom,&mspline);

        status = assembler->assembleVector(def);
        result = assembler->rhs();
        return status == ThinShellAssemblerStatus::Success;

    };

    gsParaviewCollection collection(dirname + "/" + output);
    gsParaviewCollection collectionTF(dirname + "/" + "tensionfield");
    gsMultiPatch<> deformation = mp;

    // Define solvers
    gsStaticDR<real_t> DynamicRelaxationSolver(M,F,Residual);
    gsOptionList DynamicRelaxationOptions = DynamicRelaxationSolver.options();
    DynamicRelaxationOptions.setReal("damping",damping);
    DynamicRelaxationOptions.setReal("alpha",alpha);
    DynamicRelaxationOptions.setInt("maxIt",maxitDR);
    DynamicRelaxationOptions.setReal("tolE",1e-2);
    DynamicRelaxationOptions.setReal("tolF",1e-4);
    DynamicRelaxationOptions.setReal("tolU",1e-4);
    DynamicRelaxationOptions.setInt("verbose",1);
    DynamicRelaxationSolver.setOptions(DynamicRelaxationOptions);
    DynamicRelaxationSolver.initialize();

    gsStaticOpt<real_t,gsOptim<real_t>::LBFGS> OptimizerSolver(Residual,assembler->numDofs());
    gsOptionList OptimizerSolverOptions = OptimizerSolver.options();
    OptimizerSolverOptions.setInt("maxIt",maxitOPT);
    OptimizerSolverOptions.setReal("tolF",1e-5);
    OptimizerSolverOptions.setReal("tolU",1e-5);
    OptimizerSolverOptions.setInt("verbose",1);
    OptimizerSolver.setOptions(OptimizerSolverOptions);
    OptimizerSolver.initialize();


    gsStaticNewton<real_t> NewtonSolver(K,F,Jacobian,Residual);
    gsOptionList NewtonOptions = NewtonSolver.options();
    NewtonOptions.setInt("verbose",true);
    NewtonOptions.setInt("maxIt",maxitNR);
    NewtonOptions.setReal("tolU",1e-6);
    NewtonOptions.setReal("tolF",1e-8);
    NewtonSolver.setOptions(NewtonOptions);

    std::vector<gsStaticBase<real_t> *> solvers;
    if (DR)
    {
      gsInfo<<"Using Dynamic Relaxation solver\n";
      solvers.push_back(&DynamicRelaxationSolver);
    }
    if (OP)
    {
      gsInfo<<"Using Optimizer solver\n";
      solvers.push_back(&OptimizerSolver);
    }
    if (NR)
    {
      gsInfo<<"Using Newton-Raphson solver\n";
      solvers.push_back(&NewtonSolver);
    }

    gsStaticComposite<real_t> StaticSolver(solvers);
    StaticSolver.options().setInt("verbose",1);
    StaticSolver.initialize();

    gsMatrix<> U = gsMatrix<>::Zero(assembler->numDofs(),1);

    std::string writeName = dirname + "/out.csv";
    std::ofstream file;
    file.open(writeName);
    file<<"Load factor, Moment1, Moment2, Moment3\n";
    file<<std::setprecision(10);
    file<<0<<","<<0<<"\n";
    file.close();
    index_t k=0;
    std::queue<real_t> Dcheckpoints;
    Dcheckpoints.push(0.7);
    Dcheckpoints.push(0.8);
    Dcheckpoints.push(0.9);
    Dcheckpoints.push(1.0);
    real_t dL0 = dL;
    dL = std::min(dL,Dcheckpoints.front());
    real_t D   = dL;
    while (true)
    {
        gsInfo<<"Displacement step "<<k<<"; D = "<<D<<"; dL = "<<dL<<"\n";
        dx.set_t(D);
        dy.set_t(D);
        dz.set_t(D);
        assembler->updateBCs(BCs);
        StaticSolver.setDisplacement(U);
        StaticSolver.solve();

        if ((StaticSolver.status() != gsStatus::Success && dL/dL0 > math::pow(0.5,3)))
        {
            dL = dL/2;
            D = D - dL;
            continue;
        }

        U = StaticSolver.solution();

        // 1. Get all the coefficients (including the ones from the eliminated BCs.)
        // gsMatrix<real_t> solFull = assembler->fullSolutionVector(controlDC.solutionU());
        gsMatrix<real_t> solFull = assembler->fullSolutionVector(U);

        // 2. Reshape all the coefficients to a Nx3 matrix
        size_t d = geom.targetDim();
        GISMO_ASSERT(solFull.rows() % d==0,"Rows of the solution vector does not match the number of control points");
        solFull.resize(solFull.rows()/d,d);

        // 3. Make the mapped spline
        gsMappedSpline<2,real_t> mspline(bb2,solFull);

        // 4. Create deformation spline
        gsFunctionSum<real_t> def(&geom,&mspline);

        if (plot)
        {
            // 5. Plot the mapped spline on the original geometry
            gsPiecewiseFunction<> displacements;
            assembler->constructStress(def,displacements,stress_type::displacement);
            std::string fileName = dirname + "/" + output + util::to_string(k);
            gsWriteParaview(def,displacements,fileName,10000,"_");
            fileName = gsFileManager::getFilename(fileName);
            for (size_t p=0; p<geom.nPatches(); ++p)
                collection.addPart(fileName + "_" + util::to_string(p) + ".vts",k,"solution",p);

            // 5. Construct stress
            gsPiecewiseFunction<> tensionFields;
            assembler->constructStress(def,tensionFields,stress_type::tension_field);
            std::string fileNameTF = dirname + "/" + "tensionfield" + util::to_string(k);
            gsWriteParaview(def,tensionFields,fileNameTF,10000,"_");
            fileNameTF = gsFileManager::getFilename(fileNameTF);
            for (size_t p=0; p<geom.nPatches(); ++p)
                collectionTF.addPart(fileNameTF + "_" + util::to_string(p) + ".vts",k,"solution",p);


            index_t N = 100;
            gsMatrix<> coords(2,N), supp, result, tmp;
            index_t dir = 1;
            std::vector<real_t> par;
            if (testCase==0)
                par = {0.5,0.7,0.9};
            else if (testCase==1)
                par = {0.5,0.75};
            else
                GISMO_ERROR("Test case "<<testCase<<" not implemented");

            for (size_t p=0; p!=par.size(); p++)
            {
                result.resize(mp.geoDim(),mp.nPatches()*N);
                for (index_t patch=0; patch!=def.nPieces(); patch++)
                {
                    supp = def.piece(patch).support();
                    coords.row(dir).setConstant( par[p]*( supp(dir,1)-supp(dir,0) )+supp(dir,0) );
                    coords.row(1-dir).setLinSpaced(N,0,1);
                    coords.row(1-dir).array() *= ( supp(1-dir,1)-supp(1-dir,0) );
                    coords.row(1-dir).array() += supp(1-dir,0);
                    def.piece(patch).eval_into(coords,tmp);
                    result.block(0,patch*N,tmp.rows(),tmp.cols()) = tmp;
                }
                result.transposeInPlace();
                writeToCSVfile(dirname + "/" + "line_D=" + std::to_string(D) + "_u=" + util::to_string(par[p]) + "_dir=" + std::to_string(dir) + ".csv",result);
            }
        }

        /////////////////////////////////////////////////////////////////////////////////
        // MOMENT COMPUTATION ///////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////
        // See: https://comet-fenics.readthedocs.io/en/latest/demo/tips_and_tricks/computing_reactions.html

        // Define vectors representing the distance from the origin to the boundary for the moment around the x and y axes.
        std::vector<gsFunctionExpr<>> Mz_x_vec(3);
        std::vector<gsFunctionExpr<>> Mz_y_vec(3);
        Mz_x_vec[0] = gsFunctionExpr<>("sqrt(x^2+y^2)*cos(atan2(y,x))",3);
        Mz_x_vec[1] = dx;
        Mz_x_vec[2] = gsFunctionExpr<>("sqrt(x^2+y^2)*cos(atan2(y,x)+pi/2*t)",3);
        Mz_x_vec[2].set_t(D);

        Mz_y_vec[0] = gsFunctionExpr<>("-sqrt(x^2+y^2)*sin(atan2(y,x))",3);
        Mz_y_vec[1] = dy;
        Mz_y_vec[2] = gsFunctionExpr<>("sqrt(x^2+y^2)*sin(atan2(y,x)+pi/2*t)",3);
        Mz_y_vec[2].set_t(D);

        std::vector<real_t> Moments(3);

        for (index_t k=0; k!=Mz_x_vec.size(); k++)
        {
            gsExprAssembler<> A(1,1);
            A.setIntegrationElements(dbasis);
            gsExprEvaluator<> evA(A);
            gsExprAssembler<>::space u = A.getSpace(dbasis,3);

            // Define dummy boundary conditions.
            // These conditions are used to define the test space on the boundary and compute its values there.
            gsBoundaryConditions<> bc_dummy;
            bc_dummy.addCondition(0,boundary::north, condition_type::dirichlet, &Mz_y_vec[k], 0 ,false,0);
            bc_dummy.addCondition(1,boundary::north, condition_type::dirichlet, &Mz_y_vec[k], 0 ,false,0);
            bc_dummy.addCondition(2,boundary::north, condition_type::dirichlet, &Mz_y_vec[k], 0 ,false,0);
            bc_dummy.addCondition(3,boundary::north, condition_type::dirichlet, &Mz_y_vec[k], 0 ,false,0);

            bc_dummy.addCondition(0,boundary::north, condition_type::dirichlet, &Mz_x_vec[k], 0 ,false,1);
            bc_dummy.addCondition(1,boundary::north, condition_type::dirichlet, &Mz_x_vec[k], 0 ,false,1);
            bc_dummy.addCondition(2,boundary::north, condition_type::dirichlet, &Mz_x_vec[k], 0 ,false,1);
            bc_dummy.addCondition(3,boundary::north, condition_type::dirichlet, &Mz_x_vec[k], 0 ,false,1);

            bc_dummy.setGeoMap(mp);
            u.setup(bc_dummy,dirichlet::interpolation,-1);

            // Set all coefficients corresponding to the space to zero, except for the boundary conditions.
            gsMultiPatch<> geo;
            for (size_t b=0; b!=dbasis.nBases(); b++)
            {
                gsMatrix<> coefs(dbasis.basis(b).size(),mp.geoDim());
                coefs.setZero();
                for (size_t d=0; d!=mp.geoDim(); d++)
                    for (index_t k=0; k!=dbasis.basis(b).size(); k++)
                    {
                        const index_t ii = u.mapper().index(k,b,d);
                        if (u.mapper().is_free(k,b,d))
                            coefs(k,d) = 0.0;
                        else if (u.mapper().is_boundary(k,b,d))
                            coefs(k,d) = u.fixedPart().at(u.mapper().global_to_bindex(ii));
                        else
                            GISMO_ERROR("WHAT?");
                    }
                geo.addPatch(*dbasis.basis(b).makeGeometry(give(coefs)));
            }

            // Cleanup the assembler
            A.initSystem();
            // Set the geometry map and the deformation map.
            gsExprAssembler<>::geometryMap ori   = A.getMap(mp);
            gsExprAssembler<>::geometryMap deff  = A.getMap(def);
            // Register the custom solution with coefficients defined above to the assembler.
            auto v = A.getCoeff(geo);
            // Add the stress tensors based on the materials defined for the shell solver
            gsMaterialMatrixIntegrate<real_t,MaterialOutput::VectorN> S0f(materialMatrixContainer,&mp,&def);
            gsMaterialMatrixIntegrate<real_t,MaterialOutput::VectorM> S1f(materialMatrixContainer,&mp,&def);
            auto S0  = A.getCoeff(S0f);
            auto S1  = A.getCoeff(S1f);
            // Helper matrix for flexural components (see gsThinShellAssembler::assembleMatrix)
            gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","2",2);
            auto m2 = A.getCoeff(mult2t);
            // Define the components used in the variational formulation, using the solution v as the deformation
            auto N  = S0.tr();
            auto dEm= flat( jac(deff).tr() * jac(v) ) ;
            auto M  = S1.tr();
            auto dEf= -( deriv2(v,sn(deff).normalized().tr() ) + deriv2(deff,var1(v,deff) ) ) * reshape(m2,3,3);
            // Integrate the variational formulation, this gives the moments.
            Moments[k] = evA.integral((
                                    N*dEm.tr()
                                    + M*dEf.tr()
                                    )*meas(ori)
                                    );
        }

        file.open(writeName,std::ofstream::out | std::ofstream::app);
        file<< std::setprecision(10)
            << D << ","
            << Moments[0] << ","
            << Moments[1] << ","
            << Moments[2] << "\n";
        file.close();

        if (gsClose(D,Dcheckpoints.front(),1e-8))
        {
            Dcheckpoints.pop();
            if (Dcheckpoints.empty())
                break;
        }

        // GO TO THE NEXT STEPs
        dL = std::min(dL0,Dcheckpoints.front()-D);
        D = D + dL;
        k++;
    }

    if (plot)
    {
        collection.save();
        collectionTF.save();
    }

    delete assembler;
    return EXIT_SUCCESS;
}



#else//gsKLShell_ENABLED
int main(int argc, char *argv[])
{
    gsWarn<<"G+Smo is not compiled with the gsKLShell module.";
    return EXIT_FAILURE;
}
#endif
