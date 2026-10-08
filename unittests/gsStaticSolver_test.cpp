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
#include <gsStructuralAnalysis/src/gsStaticSolvers/gsControlDisplacement.h>
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

    // Decoupled cubic springs: R(x) = f - k x - c x^3, K(x) = diag(k + 3 c x^2).
    // A solve that exhausts maxIt must leave solution() untouched, and the DR
    // maxIt = N must give N iterations and iterations() == N.
    TEST(StaticSolver_MaxIt_NoCommit)
    {
        const index_t n = 3;
        gsVector<> k(n), f(n), M(n);
        k << 1.0, 2.0, 4.0;
        f << 3.0, 2.0, 1.0;
        M.setConstant(1.0);
        const real_t c = 10.0;

        index_t nRes = 0;
        gsStructuralAnalysisOps<real_t>::Residual_t Residual = [&](const gsVector<> & x, gsVector<> & R)
        {
            ++nRes;
            R = f - k.cwiseProduct(x) - c * x.array().cube().matrix();
            return true;
        };
        gsStructuralAnalysisOps<real_t>::Jacobian_t Jacobian = [&](const gsVector<> & x, gsSparseMatrix<> & K)
        {
            K.resize(n,n);
            for (index_t i = 0; i != n; ++i)
                K.insert(i,i) = k(i) + 3 * c * x(i) * x(i);
            K.makeCompressed();
            return true;
        };

        // Dynamic relaxation: N iterations, N+1 residual evaluations (the start and one per iteration)
        {
            const index_t N = 3;
            gsStaticDR<real_t> DRM(M,f,Residual);
            gsOptionList opt = DRM.options();
            opt.setReal("damping",1.0);
            opt.setReal("alpha",2.0);
            opt.setInt("maxIt",N);
            opt.setReal("tolF",1e-14);
            opt.setReal("tolU",1e-14);
            opt.setReal("tolE",1e-14);
            opt.setInt("verbose",0);
            DRM.setOptions(opt);
            DRM.initialize();
            nRes = 0;
            gsStatus status = DRM.solve();
            CHECK(status == gsStatus::NotConverged);
            CHECK_EQUAL(N, DRM.iterations());
            CHECK_EQUAL(N+1, nRes);
            CHECK_EQUAL(0.0, DRM.solution().norm());
        }

        // Newton: N updates, N+1 residual evaluations; the last iterate is not committed
        {
            const index_t N = 2;
            gsSparseMatrix<> K0(n,n);
            gsStaticNewton<real_t> NWT(K0,f,Jacobian,Residual);
            gsOptionList opt = NWT.options();
            opt.setInt("maxIt",N);
            opt.setReal("tolF",1e-14);
            opt.setReal("tolU",1e-14);
            opt.setInt("verbose",0);
            NWT.setOptions(opt);
            NWT.initialize();
            nRes = 0;
            gsStatus status = NWT.solve();
            CHECK(status == gsStatus::NotConverged);
            CHECK_EQUAL(N-1, NWT.iterations());
            CHECK_EQUAL(N+1, nRes);
            CHECK_EQUAL(0.0, NWT.solution().norm());

            // With enough iterations the same solver converges and commits
            opt.setInt("maxIt",20);
            opt.setReal("tolF",1e-10);
            opt.setReal("tolU",1e-8);
            NWT.setOptions(opt);
            NWT.initialize();
            status = NWT.solve();
            CHECK(status == gsStatus::Success);
            gsVector<> R;
            Residual(NWT.solution(),R);
            CHECK(R.norm() < 1e-8);
        }
    }

    // A throwing residual after the linear-predictor branch leaves the displacements as they were
    TEST(StaticSolver_Newton_ThrowNoCommit)
    {
        const index_t n = 3;
        gsSparseMatrix<> K0(n,n);
        for (index_t i = 0; i != n; ++i)
            K0.insert(i,i) = 1.0;
        gsVector<> f(n);
        f << 1.0, 2.0, 3.0;
        index_t nRes = 0;
        gsStructuralAnalysisOps<real_t>::Residual_t Residual = [&](const gsVector<> & x, gsVector<> & R)
        {
            if (++nRes == 2)
                GISMO_ERROR("residual failure");
            R = f - x;
            return true;
        };
        gsStructuralAnalysisOps<real_t>::Jacobian_t Jacobian = [&](const gsVector<> &, gsSparseMatrix<> & K)
        {
            K = K0;
            return true;
        };
        gsStaticNewton<real_t> NWT(K0,f,Jacobian,Residual);
        NWT.initialize();
        gsVector<> U0(n);
        U0 << 0.1, 0.2, 0.3;
        NWT.setDisplacement(U0);
        NWT.setUpdate(gsVector<>());
        gsStatus status = NWT.solve();
        CHECK(status == gsStatus::OtherError);
        CHECK_EQUAL(0.0, (NWT.solution() - U0).norm());
    }

    // Coupled linear Dirichlet ramp through the Load-constraint preset: one update per step
    TEST(StaticSolver_ControlDisplacement_LoadRamp)
    {
        const index_t n = 20;
        gsSparseMatrix<> Kl(n,n);
        for (index_t i = 0; i != n; ++i)
        {
            Kl.insert(i,i) = 2.3 + 0.01 * i;
            if (i > 0)
            {
                Kl.insert(i,i-1) = -1.1;
                Kl.insert(i-1,i) = -1.1;
            }
        }
        Kl.makeCompressed();
        gsVector<> g(n);
        for (index_t i = 0; i != n; ++i)
            g(i) = 1.7 * math::sin(1.0 + i);

        const real_t steps[3] = {0.5, 1e-2, 1e-6};
        for (index_t s = 0; s != 3; ++s)
        {
            const real_t dL = steps[s];
            gsVectorStateSpace<real_t> S(gsVector<>::Zero(n));
            gsContinuationProblem<real_t> P;
            P.residual = [&](real_t lam, gsVector<> & R, gsEvalMode)
            { R = lam * g - Kl * S.working(); return true; };
            P.jacobian = [&](real_t, gsSparseMatrix<> & K) { K = Kl; return true; };
            P.loadDirection = [&](real_t, gsVector<> & d) { d = g; return true; };
            P.F = gsVector<>::Zero(n);

            gsControlDisplacement<real_t> CD(P,S);
            gsStatus status = CD.step(dL);
            CHECK(status == gsStatus::Success);
            CHECK_EQUAL(1, CD.driver().numIterations());
            CHECK_CLOSE(dL, CD.solutionL(), 1e-14);
            const gsVector<> r = Kl * CD.solutionU() - dL * g;
            CHECK(r.norm() / (dL * g.norm()) < 1e-12);

            // A second step continues from the committed state
            status = CD.step(dL);
            CHECK(status == gsStatus::Success);
            CHECK_CLOSE(2*dL, CD.solutionL(), 1e-14);
            const gsVector<> r2 = Kl * CD.solutionU() - 2 * dL * g;
            CHECK(r2.norm() / (2 * dL * g.norm()) < 1e-12);
        }
    }

    // Decoupled cubic springs R(x) = f - k x - c x^3 and their tangent
    struct StaticFacadeSprings
    {
        gsVector<> k, f, M;
        real_t c;
        gsStructuralAnalysisOps<real_t>::Residual_t Residual;
        gsStructuralAnalysisOps<real_t>::Jacobian_t Jacobian;

        StaticFacadeSprings() : k(3), f(3), M(gsVector<>::Ones(3)), c(10.0)
        {
            k << 1.0, 2.0, 4.0;
            f << 3.0, 2.0, 1.0;
            Residual = [this](const gsVector<> & x, gsVector<> & R)
            {
                R = f - k.cwiseProduct(x) - c * x.array().cube().matrix();
                return true;
            };
            Jacobian = [this](const gsVector<> & x, gsSparseMatrix<> & K)
            {
                K.resize(3,3);
                for (index_t i = 0; i != 3; ++i)
                    K.insert(i,i) = k(i) + 3 * c * x(i) * x(i);
                K.makeCompressed();
                return true;
            };
        }
    private:
        StaticFacadeSprings(const StaticFacadeSprings &);
        StaticFacadeSprings & operator=(const StaticFacadeSprings &);
    };

    // A solve that exhausts maxIt from a non-zero displacement leaves solution()
    // at that displacement while the discarded update is not zero, so a commit
    // on failure would be visible. The same solver then converges and commits.
    template<class S>
    void staticFacadeNewtonNoCommit()
    {
        StaticFacadeSprings sp;
        gsVector<> U0(3);
        U0 << 0.3, -0.2, 0.1;
        gsSparseMatrix<> K0(3,3);
        for (index_t i = 0; i != 3; ++i)
            K0.insert(i,i) = sp.k(i);
        S NWT(K0, sp.f, sp.Jacobian, sp.Residual);
        gsOptionList opt = NWT.options();
        opt.setInt("maxIt",2);
        opt.setReal("tolF",1e-14);
        opt.setReal("tolU",1e-14);
        opt.setInt("verbose",0);
        NWT.setOptions(opt);
        NWT.initialize();
        NWT.setDisplacement(U0);
        CHECK(NWT.solve() == gsStatus::NotConverged);
        CHECK_EQUAL(1, NWT.iterations());
        CHECK(NWT.update().norm() > 0);
        CHECK_EQUAL(U0.size(), NWT.solution().size());
        CHECK((NWT.solution() - U0).norm() == 0);

        const real_t discarded = NWT.update().norm();

        opt.setInt("maxIt",20);
        opt.setReal("tolF",1e-10);
        opt.setReal("tolU",1e-8);
        NWT.setOptions(opt);
        CHECK(NWT.solve() == gsStatus::Success);
        gsVector<> R;
        sp.Residual(NWT.solution(), R);
        CHECK(R.norm() < 1e-8);
        gsInfo << "StaticFacade Newton: discarded update " << discarded << "\n";
    }

    template<class S>
    void staticFacadeDRNoCommit()
    {
        StaticFacadeSprings sp;
        gsVector<> U0(3);
        U0 << 0.3, -0.2, 0.1;
        S DRM(sp.M, sp.f, sp.Residual);
        gsOptionList opt = DRM.options();
        opt.setReal("damping",1.0);
        opt.setReal("alpha",2.0);
        opt.setInt("maxIt",3);
        opt.setReal("tolF",1e-14);
        opt.setReal("tolU",1e-14);
        opt.setReal("tolE",1e-14);
        opt.setInt("verbose",0);
        DRM.setOptions(opt);
        DRM.initialize();
        DRM.setDisplacement(U0);
        CHECK(DRM.solve() == gsStatus::NotConverged);
        CHECK_EQUAL(3, DRM.iterations());
        CHECK(DRM.update().norm() > 0);
        CHECK_EQUAL(U0.size(), DRM.solution().size());
        CHECK((DRM.solution() - U0).norm() == 0);
        gsInfo << "StaticFacade DR: discarded update " << DRM.update().norm() << "\n";
    }

    TEST(StaticFacade_NoCommit_nonzero_state)
    {
        staticFacadeNewtonNoCommit<gsStaticNewton<real_t> >();
        staticFacadeDRNoCommit<gsStaticDR<real_t> >();
    }

    // Chain of three springs s_i(e) = k_i e + c e^3 between a fixed node 0 and a
    // node 3 displaced by lambda; the nodes 1 and 2 are free. A nonlinear
    // Dirichlet ramp: the load direction dR/dlambda depends on the state.
    struct StaticFacadeChain
    {
        gsVectorStateSpace<real_t> space;
        gsContinuationProblem<real_t> P;

        static real_t kk(index_t i) { return i == 0 ? 1.0 : (i == 1 ? 2.0 : 4.0); }
        static real_t cc() { return 20.0; }
        static real_t s (index_t i, real_t e) { return kk(i) * e + cc() * e * e * e; }
        static real_t ds(index_t i, real_t e) { return kk(i) + 3 * cc() * e * e; }

        StaticFacadeChain() : space(2)
        {
            P.F = gsVector<real_t>::Zero(2);
            P.residual = [this](real_t lam, gsVector<real_t> & R, gsEvalMode) -> bool
            {
                const gsVector<real_t> & x = space.working();
                R.resize(2);
                R << s(1, x[1] - x[0]) - s(0, x[0]), s(2, lam - x[1]) - s(1, x[1] - x[0]);
                return true;
            };
            P.jacobian = [this](real_t lam, gsSparseMatrix<real_t> & K) -> bool
            {
                const gsVector<real_t> & x = space.working();
                const real_t a = ds(0, x[0]), b = ds(1, x[1] - x[0]), d = ds(2, lam - x[1]);
                K.resize(2,2);
                K.insert(0,0) = a + b;  K.insert(0,1) = -b;
                K.insert(1,0) = -b;     K.insert(1,1) = b + d;
                K.makeCompressed();
                return true;
            };
            P.loadDirection = [this](real_t lam, gsVector<real_t> & d) -> bool
            {
                const gsVector<real_t> & x = space.working();
                d.resize(2);
                d << 0.0, ds(2, lam - x[1]);
                return true;
            };
        }
    private:
        StaticFacadeChain(const StaticFacadeChain &);
        StaticFacadeChain & operator=(const StaticFacadeChain &);
    };

    // The acceptance gate of gsControlDisplacement: at least one update,
    // relative residual and relative update below 1e-6, at most 25 updates
    static memory::shared_ptr< gsConvergenceTest<real_t> > staticFacadeGate(const memory::shared_ptr< gsStaticResidualTest<real_t> > & rt,
                                                                            const real_t tol, const index_t maxIt)
    {
        typedef gsConvergenceTest<real_t>::Ptr TestPtr;
        memory::shared_ptr< gsAllTest<real_t> > all(new gsAllTest<real_t>());
        all->add(TestPtr(new gsMinIterTest<real_t>(1))).add(rt).add(TestPtr(new gsStaticUpdateTest<real_t>(tol)));
        memory::shared_ptr< gsAnyTest<real_t> > any(new gsAnyTest<real_t>());
        any->add(all).add(TestPtr(new gsMaxIterTest<real_t>(maxIt)));
        return any;
    }

    // gsControlDisplacement on a problem is a gsContinuationDriver with a load
    // constraint, the gate above and the residual basis floored by |dL||dR/dlambda|.
    // \a beforeStep may change the object before each of its steps.
    static void staticFacadeControlDisplacement(const std::function<void(gsControlDisplacement<real_t> &)> & beforeStep
                                                = std::function<void(gsControlDisplacement<real_t> &)>())
    {
        StaticFacadeChain rigA, rigB;
        gsControlDisplacement<real_t> cd(rigA.P, rigA.space);

        gsContinuationDriver<real_t> drv(rigB.P, rigB.space);
        drv.setConstraint(gsPathConstraint<real_t>::Ptr(new gsLoadConstraint<real_t>()));
        memory::shared_ptr< gsStaticResidualTest<real_t> > rt(new gsStaticResidualTest<real_t>(1e-6));
        drv.setTest(staticFacadeGate(rt, 1e-6, 25));
        drv.setup();

        const real_t dL = 0.1;
        index_t maxIt = 0;
        for (index_t step = 0; step < 3; ++step)
        {
            if (beforeStep)
                beforeStep(cd);
            gsVector<real_t> d;
            CHECK(rigB.P.loadDirection(drv.lambda() + dL, d));
            rt->setBasisFloor(math::abs(dL) * d.norm());
            drv.controller().setLength(dL);
            const gsStatus sd = drv.step();
            const gsStatus sc = cd.step(dL);
            // the floor binds: without it the basis would be the residual at the predicted point
            CHECK(rt->basis() > rt->predictedResidual());

            CHECK(sc == gsStatus::Success);
            CHECK(sd == gsStatus::Success);
            CHECK_EQUAL(drv.lambda(), cd.solutionL());
            CHECK_EQUAL(drv.numIterations(), cd.driver().numIterations());
            CHECK_EQUAL(rigB.space.committed().size(), cd.solutionU().size());
            CHECK((cd.solutionU() - rigB.space.committed()).norm() == 0);
            maxIt = math::max(maxIt, cd.driver().numIterations());
        }
        CHECK(maxIt >= 2);

        // after three steps the springs carry one force
        const gsVector<real_t> & x = cd.solutionU();
        CHECK_CLOSE(0.3, cd.solutionL(), 1e-14);
        CHECK(x.norm() > 0);
        const real_t f1 = StaticFacadeChain::s(0, x[0]);
        const real_t f2 = StaticFacadeChain::s(1, x[1] - x[0]);
        const real_t f3 = StaticFacadeChain::s(2, cd.solutionL() - x[1]);
        CHECK(math::abs(f1 - f3) <= 1e-6 * math::abs(f3));
        CHECK(math::abs(f2 - f3) <= 1e-6 * math::abs(f3));

        // a step that does not converge leaves the load and the displacements as they were
        const gsVector<real_t> U3 = cd.solutionU();
        const real_t L3 = cd.solutionL();
        memory::shared_ptr< gsStaticResidualTest<real_t> > strict(new gsStaticResidualTest<real_t>(1e-14));
        cd.driver().setTest(staticFacadeGate(strict, 1e-14, 1));
        cd.driver().setup();
        CHECK(cd.step(dL) == gsStatus::NotConverged);
        CHECK(cd.driver().context().DeltaU.norm() > 0);
        CHECK_EQUAL(U3.size(), cd.solutionU().size());
        CHECK((cd.solutionU() - U3).norm() == 0);
        CHECK(cd.solutionL() == L3);
        gsInfo << "StaticFacade ControlDisplacement: up to " << maxIt << " updates per step, lambda " << L3 << "\n";
    }

    TEST(StaticFacade_ControlDisplacement_matches_driver)
    {
        staticFacadeControlDisplacement();
    }

}
