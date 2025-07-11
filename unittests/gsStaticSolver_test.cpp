/** @file gsThinShellAssembler_test.cpp

    @brief Provides unittests for the gsThinShellAssembler class

    * Balloon: unit-test based on a hyperelastic balloon inflated by a follower pressure.
               This test allows to test the follower pressure as well as the hyperelastic material models (incompressible)

    * UAT:     unit-test based on a uni-axial tension test
               This test allows to test the hyperelastic material models (incompressible and compressible)

    * Modal:   unit-test based on a modal analysis
               This test allows to test the mass matrix and the compressible material model


    == BASIC REFERENCE ==
         - TEST(NAME_OF_TEST) { body_of_test }
         - TEST_FIXTURE(NAME_OF_FIXTURE,NAME_OF_TEST){ body_of_test }

    == CHECK MACRO REFERENCE ==
         - CHECK(EXPR);
         - CHECK_EQUAL(EXPECTED,ACTUAL);
         - CHECK_CLOSE(EXPECTED,ACTUAL,EPSILON);
         - CHECK_ARRAY_EQUAL(EXPECTED,ACTUAL,LENGTH);
         - CHECK_ARRAY_CLOSE(EXPECTED,ACTUAL,LENGTH,EPSILON);
         - CHECK_ARRAY2D_EQUAL(EXPECTED,ACTUAL,ROWCOUNT,COLCOUNT);
         - CHECK_ARRAY2D_CLOSE(EXPECTED,ACTUAL,ROWCOUNT,COLCOUNT,EPSILON);
         - CHECK_THROW(EXPR,EXCEPTION_TYPE_EXPECTED);

    == TIME CONSTRAINTS ==
         - UNITTEST_TIME_CONSTRAINT(TIME_IN_MILLISECONDS);
         - UNITTEST_TIME_CONSTRAINT_EXEMPT();

    == MORE INFO ==
         See: https://unittest-cpp.github.io/

    Author(s): H.M.Verhelst (2019 - ..., TU Delft)
 **/

#include "gismo_unittest.h"       // Brings in G+Smo and the UnitTest++ framework
#ifdef gsKLShell_ENABLED
#include <gsKLShell/gsKLShell.h>
#endif

#ifdef gsHLBFGS_ENABLED
#include <gsHLBFGS/gsHLBFGS.h>
#endif

#include <gsStructuralAnalysis/src/gsStaticSolvers/gsStaticDR.h>
#include <gsStructuralAnalysis/src/gsStaticSolvers/gsStaticNewton.h>
#include <gsStructuralAnalysis/src/gsStaticSolvers/gsStaticOpt.h>
#include <gsStructuralAnalysis/src/gsStaticSolvers/gsStaticComposite.h>

SUITE(gsStaticSolver_test)                 // The suite should have the same name as the file
{

#ifdef gsKLShell_ENABLED
    // solver: 0: newton, 1: DR, 2: OPT
    std::pair<real_t,real_t> UAT_analytical(const std::vector<index_t> solver, const index_t material=1, const index_t impl=1, const bool Compressibility=false);
    std::pair<real_t,real_t> UAT_numerical(const std::vector<index_t> solver, const index_t material=1, const index_t impl=1, const bool Compressibility=false);
    void UAT_CHECK(const std::vector<index_t> solver, const index_t material=1, const index_t impl=1, const bool Compressibility=false);
#endif

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef gsKLShell_ENABLED
    TEST(StaticSolver_UAT_NR)
    {
        std::vector<index_t> solver({0});
        UAT_CHECK(solver);
    }

    // TEST(StaticSolver_UAT_DR)
    // {
    //     std::vector<index_t> solver({1});
    //     UAT_CHECK(solver);
    // }
#ifdef gsHLBFGS_ENABLED
    TEST(StaticSolver_UAT_OP)
    {
        std::vector<index_t> solver({2});
        UAT_CHECK(solver);
    }
#endif

#ifdef gsHLBFGS_ENABLED
    TEST(StaticSolver_UAT_OP_NR)
    {
        std::vector<index_t> solver({2,0});
        UAT_CHECK(solver);
    }
#endif

    TEST(StaticSolver_UAT_DR_NR)
    {
        std::vector<index_t> solver({1,0});
        UAT_CHECK(solver);
    }

#ifdef gsHLBFGS_ENABLED
    TEST(StaticSolver_UAT_DR_OP)
    {
        std::vector<index_t> solver({1,2});
        UAT_CHECK(solver);
    }
#endif

    std::pair<real_t,real_t> UAT_numerical(const std::vector<index_t> solver, const index_t material, const index_t impl, const bool Compressibility)
    {
        //! [Parse command line]
        index_t numRefine  = 1;
        index_t numElevate = 1;

        real_t E_modulus = 1.0;
        real_t PoissonRatio;
        real_t Density = 1.0;
        real_t Ratio = 7.0;

        real_t mu = 1.5e6;
        real_t thickness = 0.001;

        real_t alpha1,alpha2,alpha3,mu1,mu2,mu3;
        alpha1 = 1.3;
        mu1    = 6.3e5/4.225e5*mu;
        alpha2 = 5.0;
        mu2    = 0.012e5/4.225e5*mu;
        alpha3 = -2.0;
        mu3    = -0.1e5/4.225e5*mu;

        if (!Compressibility)
          PoissonRatio = 0.5;
        else
          PoissonRatio = 0.45;

        E_modulus = 2*mu*(1+PoissonRatio);

        //! [Parse command line]

        //! [Read input file]
        gsMultiPatch<> mp, mp_def;

        mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // degree

        if (numElevate!=0)
            mp.degreeElevate(numElevate);

        // h-refine
        for (int r =0; r < numRefine; ++r)
            mp.uniformRefine();

        mp_def = mp;

        //! [Refinement]
        gsMultiBasis<> dbasis(mp);

        gsBoundaryConditions<> bc;
        bc.setGeoMap(mp);

        gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();

        real_t lambda = 2.0;
        gsConstantFunction<> displx(lambda-1.0,2);

        GISMO_ENSURE(mp.targetDim()==2,"Geometry must be planar (targetDim=2)!");
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 );

        bc.addCondition(boundary::east, condition_type::dirichlet, &displx, 0, false, 0 );

        bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 );

        //! [Refinement]

        // Linear isotropic material model
        gsVector<> tmp(2);
        tmp.setZero();
        gsConstantFunction<> force(tmp,2);
        gsFunctionExpr<> t(std::to_string(thickness),2);
        gsFunctionExpr<> E(std::to_string(E_modulus),2);
        gsFunctionExpr<> nu(std::to_string(PoissonRatio),2);
        gsFunctionExpr<> rho(std::to_string(Density),2);
        gsConstantFunction<> ratio(Ratio,2);

        gsConstantFunction<> alpha1fun(alpha1,2);
        gsConstantFunction<> mu1fun(mu1,2);
        gsConstantFunction<> alpha2fun(alpha2,2);
        gsConstantFunction<> mu2fun(mu2,2);
        gsConstantFunction<> alpha3fun(alpha3,2);
        gsConstantFunction<> mu3fun(mu3,2);

        std::vector<gsFunctionSet<>*> parameters(3);
        parameters[0] = &E;
        parameters[1] = &nu;
        parameters[2] = &ratio;
        gsMaterialMatrixBase<real_t>::uPtr materialMatrix;

        if (material==4)
        {
            parameters.resize(8);
            parameters[0] = &E;
            parameters[1] = &nu;
            parameters[2] = &mu1fun;
            parameters[3] = &alpha1fun;
            parameters[4] = &mu2fun;
            parameters[5] = &alpha2fun;
            parameters[6] = &mu3fun;
            parameters[7] = &alpha3fun;
        }

        gsOptionList options;
        if      (material==0)
        {
            GISMO_ERROR("This test is not available for SvK models");
        }
        else
        {
            options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",material);
            options.addSwitch("Compressibility","Compressibility: (false): Imcompressible | (true): Compressible",Compressibility);
            options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",impl);
            materialMatrix = getMaterialMatrix<2,real_t>(mp,t,parameters,rho,options);
        }

        gsThinShellAssemblerBase<real_t>* assembler;
        assembler = new gsThinShellAssembler<2, real_t, false >(mp,dbasis,bc,force,materialMatrix);

        assembler->setPointLoads(pLoads);

        // Function for the Jacobian
        typedef std::function<gsSparseMatrix<real_t> (gsVector<real_t> const &)>    Jacobian_t;
        typedef std::function<gsVector<real_t> (gsVector<real_t> const &) >         Residual_t;
        // Function for the Jacobian
        gsStructuralAnalysisOps<real_t>::Jacobian_t Jacobian = [&assembler,&mp_def](gsVector<real_t> const &x, gsSparseMatrix<real_t> & m)
        {
        ThinShellAssemblerStatus status;
        assembler->constructSolution(x,mp_def);
        status = assembler->assembleMatrix(mp_def);
        m = assembler->matrix();
        return status == ThinShellAssemblerStatus::Success;
        };

        gsStructuralAnalysisOps<real_t>::Residual_t Residual = [&assembler,&mp_def](gsVector<real_t> const &x, gsVector<real_t> & result)
        {
        ThinShellAssemblerStatus status;
        assembler->constructSolution(x,mp_def);
        status = assembler->assembleVector(mp_def);
        result = assembler->rhs();
        return status == ThinShellAssemblerStatus::Success;
        };

        // Define Matrices
        assembler->assemble();
        gsSparseMatrix<> K = assembler->matrix();
        gsVector<> F = assembler->rhs();
        assembler->assembleMass(true);
        gsVector<> M = assembler->rhs();

        gsStaticDR<real_t> DRM(M,F,Residual);
        gsOptionList DROptions = DRM.options();
        DROptions.setReal("damping",0.0);
        DROptions.setReal("alpha",1e9);
        DROptions.setInt("maxIt",1e6);
        DROptions.setReal("tolF",1e-8);
        DROptions.setReal("tolU",1e-6);
        DROptions.setReal("tolE",1e-4);
        DROptions.setInt("verbose",0);
        DRM.setOptions(DROptions);
        DRM.initialize();

#ifdef gsHLBFGS_ENABLED
        gsStaticOpt<real_t,gsHLBFGS<real_t>> OPT(Residual,assembler->numDofs());
#else
        gsStaticOpt<real_t,gsGradientDescent<real_t>> OPT(Residual,assembler->numDofs());
#endif
        gsOptionList OPTOptions = OPT.options();
        OPTOptions.setInt("maxIt",100);
        OPTOptions.setReal("tolF",1e-8);
        OPTOptions.setReal("tolU",1e-6);
        OPTOptions.setInt("verbose",0);
        OPT.setOptions(OPTOptions);
        OPT.initialize();

        gsStaticNewton<real_t> NWT(K,F,Jacobian,Residual);
        gsOptionList NWTOptions = NWT.options();
        NWTOptions.setInt("maxIt",20);
        NWTOptions.setReal("tolF",1e-8);
        NWTOptions.setReal("tolU",1e-6);
        NWTOptions.setInt("verbose",0);
        NWT.setOptions(NWTOptions);
        NWT.initialize();

        std::vector<gsStaticBase<real_t>*> solverVec(solver.size());
        for (size_t s=0; s!=solver.size(); s++)
            if (solver[s]==0)
                solverVec[s] = &NWT;
            else if (solver[s]==1)
                solverVec[s] = &DRM;
            else if (solver[s]==2)
                solverVec[s] = &OPT;

        gsStaticComposite<real_t> compositeSolver(solverVec);
        compositeSolver.initialize();
        compositeSolver.solve();
        CHECK(compositeSolver.status() == gsStatus::Success);
        mp_def = assembler->constructSolution(compositeSolver.solution());

        gsMultiPatch<> deformation = mp_def;
        for (size_t k = 0; k != mp_def.nPatches(); ++k)
            deformation.patch(k).coefs() -= mp.patch(k).coefs();

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Check solutions
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // NOTE: all analytical solutions for compressible materials are fixed for displ=1; (lambda=2)

        // Compute stretches (should be the same everywhere)
        // Ordering: lambda(0) < lambda(1); lambda(2) is ALWAYS the through-thickness stretch
        gsVector<> pt(2);
        pt<<1,0;
        gsMatrix<> lambdas = assembler->computePrincipalStretches(pt,mp_def,0);

        // Get the total force on the tension boundary
        patchSide ps(0,boundary::east);
        gsMatrix<> forceVector = assembler->boundaryForce(mp_def,ps);
        real_t sideForce = forceVector.sum();
        real_t S   = -sideForce / (thickness*lambdas(0)*lambdas(2));
        real_t L   = lambdas(0);

        std::pair<real_t,real_t> result;
        result.first = L;
        result.second = S;

        delete assembler;

        return result;
    }

    std::pair<real_t,real_t> UAT_analytical(const std::vector<index_t>, const index_t material, const index_t, const bool Compressibility)
    {
        real_t PoissonRatio;
        real_t Ratio = 7.0;

        real_t mu = 1.5e6;

        real_t alpha1,alpha2,alpha3,mu1,mu2,mu3;
        alpha1 = 1.3;
        mu1    = 6.3e5/4.225e5*mu;
        alpha2 = 5.0;
        mu2    = 0.012e5/4.225e5*mu;
        alpha3 = -2.0;
        mu3    = -0.1e5/4.225e5*mu;

        if (!Compressibility)
          PoissonRatio = 0.5;
        else
          PoissonRatio = 0.45;

        real_t lambda = 2.0;

        real_t San,J,K,Lan;
        if      (material==1 && Compressibility)
        {
            K = 2*mu*(1+PoissonRatio)/(3-6*PoissonRatio);
            J = 1.105598565;// specific for lambda==2!!
            San = lambda*(0.5*mu*(-(2*(math::pow(lambda,2)+2*J/lambda))/(3*math::pow(J,2./3.)*lambda)+2*lambda/math::pow(J,2./3.))+0.25*K*(2*math::pow(J,2)/lambda-2./lambda))/J;
            Lan = math::pow(J/lambda,0.5);
        }
        else if (material==1 && !Compressibility)
        {
            San = mu * (lambda*lambda - 1/lambda);
            Lan = math::pow(1./lambda,0.5);
        }
        else if (material==3 && Compressibility)
        {
            real_t c2 = 1.0 / (Ratio+1);
            real_t c1 = 1.0 - c2;
            K = 2*mu*(1+PoissonRatio)/(3-6*PoissonRatio);
            J = 1.099905842;// specific for lambda==2!!
            San = lambda*(0.5*c1*mu*(-(2*(math::pow(lambda,2)+2*J/lambda))/(3*math::pow(J,2./3.)*lambda)+2*lambda/math::pow(J,2./3.))+0.5*c2*mu*(-(4*(2*lambda*J+math::pow(J,2)/math::pow(lambda,2)))/(3*math::pow(J,4./3.)*lambda)+4/math::pow(J,1./3.))+0.25*K*(2*math::pow(J,2)/lambda-2/lambda))/J;
            Lan = math::pow(J/lambda,0.5);
        }
        else if (material==3 && !Compressibility)
        {
            real_t c2 = 1.0 / (Ratio+1);
            real_t c1 = 1.0 - c2;
            San =-mu*(c2*lambda*lambda+c2/lambda+c1)/lambda+lambda*(c1*lambda*mu+2*c2*mu);
            Lan = math::pow(1./lambda,0.5);
        }
        else if (material==4 && Compressibility)
        {
            K = 2*mu*(1+PoissonRatio)/(3-6*PoissonRatio);
            J = 1.088778638;// specific for lambda==2!!
            San = 1./J* (lambda *( mu1*(2*math::pow(lambda/math::pow(J,1./3.),alpha1)*alpha1/(3*lambda)-2*math::pow(math::pow(J/lambda,0.5)/math::pow(J,1./3.),alpha1)*alpha1/(3*lambda))/alpha1+mu2*(2*math::pow(lambda/math::pow(J,1./3.),alpha2)*alpha2/(3*lambda)-2*math::pow(math::pow(J/lambda,0.5)/math::pow(J,1./3.),alpha2)*alpha2/(3*lambda))/alpha2+mu3*(2*math::pow(lambda/math::pow(J,1./3.),alpha3)*alpha3/(3*lambda)-2*math::pow(math::pow(J/lambda,0.5)/math::pow(J,1./3.),alpha3)*alpha3/(3*lambda))/alpha3+0.25*K*(2*math::pow(J,2)/lambda-2/lambda) ) );
            Lan = math::pow(J/lambda,0.5);
        }
        else if (material==4 && !Compressibility)
        {
            San =-mu1*math::pow((1./lambda),0.5*alpha1)-mu2*math::pow((1./lambda),0.5*alpha2)-mu3*math::pow((1./lambda),0.5*alpha3)+mu1*math::pow(lambda,alpha1)+mu2*math::pow(lambda,alpha2)+mu3*math::pow(lambda,alpha3);
            Lan = math::pow(1./lambda,0.5);
        }
        else
            GISMO_ERROR("Material not treated");

        std::pair<real_t,real_t> result;
        result.first = Lan;
        result.second = San;
        return result;
    }

    void UAT_CHECK(const std::vector<index_t> solver, const index_t material, const index_t impl, const bool Compressibility)
    {
        if (material==4 && impl!=3)
          CHECK(true);

        real_t Lnum, Snum, Lana, Sana;
        std::tie(Lnum,Snum) = UAT_numerical(solver,material,impl,Compressibility);
        std::tie(Lana,Sana) = UAT_analytical(solver,material,impl,Compressibility);
        CHECK_CLOSE(std::abs(Lnum-Lana)/Lana,0,1e-7);
    }
#endif

}
