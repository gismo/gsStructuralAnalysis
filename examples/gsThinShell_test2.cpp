/** @file gsThinShell_test2.cpp

    @brief Example testing and debugging thin shell solver. Based on gsThinShell_test

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst
*/

#include <gismo.h>
#include <gismo_dev.h>

#include <gsThinShell2/gsThinShellAssembler.h>

//#include <gsThinShell/gsNewtonIterator.h>

using namespace gismo;


// Choose among various shell examples, default = Thin Plate
int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    index_t numRefine  = 1;
    index_t numElevate = 1;
    index_t testCase = 1;
    bool nonlinear = false;
    std::string fn("pde/poisson2d_bvp.xml");

    real_t E_modulus = 1.0;
    real_t PoissonRatio = 0.0;
    real_t thickness = 1.0;

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving",  numRefine );
    cmd.addInt( "t", "testCase", "Test case to run: 1 = unit square; 2 = Scordelis Lo Roof",  testCase );
    cmd.addString( "f", "file", "Input XML file", fn );
    cmd.addSwitch("nl", "Solve nonlinear problem", nonlinear);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Read input file]
    gsMultiPatch<> mp;
    gsMultiPatch<> mp_def;
    if (testCase==1)
    {
        // Unit square
        mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // degree
        mp.addAutoBoundaries();
        mp.embed(3);
        E_modulus = 1.0;
        thickness = 0.5;

    }
    else if (testCase == 2 || testCase == 3)
    {
        thickness = 0.125;
        E_modulus = 4.32E8;
        fn = "../extensions/unsupported/filedata/scordelis_lo_roof.xml";
        gsReadFile<>(fn, mp);
    }
    else if (testCase==9)
    {
        gsFileData<> fd(fn);
        gsInfo << "Loaded file "<< fd.lastPath() <<"\n";
        // Annulus
        fd.getId(0, mp); // id=0: Multipatch domain
        mp.embed(3);
        E_modulus = 1.0;
        thickness = 0.5;
    }
    //! [Read input file]

    // p-refine
    if (numElevate!=0)
        mp.degreeElevate(numElevate);

    // h-refine
    for (int r =0; r < numRefine; ++r)
        mp.uniformRefine();

    mp_def = mp;

    //! [Refinement]
    gsMultiBasis<> dbasis(mp);

    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
    gsInfo << dbasis.basis(0)<<"\n";

    gsBoundaryConditions<> bc;
    gsVector<> tmp(3);
    tmp << 0, 0, 0;

    gsVector<> neu(3);
    neu << 0, 0, 0;

    gsConstantFunction<> neuData(neu,3);
    if (testCase == 1)
    {
        for (index_t i=0; i!=3; ++i)
        {
            bc.addCondition(boundary::north, condition_type::dirichlet, 0, i ); // unknown 0 - x
            bc.addCondition(boundary::east, condition_type::dirichlet, 0, i ); // unknown 1 - y
            bc.addCondition(boundary::south, condition_type::dirichlet, 0, i ); // unknown 2 - z
            bc.addCondition(boundary::west, condition_type::dirichlet, 0, i ); // unknown 2 - z
        }
        tmp << 0,0,-1;
    }
    else if (testCase == 2)
    {
        // Diaphragm conditions
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 1 ); // unknown 1 - y
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 2 ); // unknown 2 - z

        // ORIGINAL
        // bc.addCornerValue(boundary::southwest, 0.0, 0, 0); // (corner,value, patch, unknown)

        bc.addCondition(boundary::east, condition_type::dirichlet, 0, 1 ); // unknown 1 - y
        bc.addCondition(boundary::east, condition_type::dirichlet, 0, 2 ); // unknown 2 - z

        // NOT ORIGINAL
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0 ); // unknown 1 - y
        bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0 ); // unknown 1 - y

        // Surface forces
        tmp << 0, 0, -90;
    }
    else if (testCase == 3)
    {
        neu << 0, 0, -90;
        neuData.setValue(neu,3);
        // Diaphragm conditions
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 1 ); // unknown 1 - y
        // bc.addCondition(boundary::west, condition_type::dirichlet, 0, 2 ); // unknown 2 - z
        bc.addCondition(boundary::west, condition_type::neumann, &neuData );

        // ORIGINAL
        // bc.addCornerValue(boundary::southwest, 0.0, 0, 0); // (corner,value, patch, unknown)

        bc.addCondition(boundary::east, condition_type::dirichlet, 0, 1 ); // unknown 1 - y
        bc.addCondition(boundary::east, condition_type::dirichlet, 0, 2 ); // unknown 2 - z

        // NOT ORIGINAL
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0 ); // unknown 1 - y
        bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0 ); // unknown 1 - y

        // Surface forces
        tmp << 0, 0, -90;


    }
    else if (testCase == 9)
    {
        for (index_t i=0; i!=3; ++i)
        {
            // patch 0
            bc.addCondition(0,boundary::north, condition_type::dirichlet, 0, i ); // unknown 0 - x
            // bc.addCondition(0,boundary::east, condition_type::dirichlet, 0, i ); // unknown 1 - y
            bc.addCondition(0,boundary::south, condition_type::dirichlet, 0, i ); // unknown 2 - z
            bc.addCondition(0,boundary::west, condition_type::dirichlet, 0, i ); // unknown 2 - z
            // patch 1
            bc.addCondition(1,boundary::north, condition_type::dirichlet, 0, i ); // unknown 0 - x
            bc.addCondition(1,boundary::east, condition_type::dirichlet, 0, i ); // unknown 1 - y
            bc.addCondition(1,boundary::south, condition_type::dirichlet, 0, i ); // unknown 2 - z
            // bc.addCondition(1,boundary::west, condition_type::dirichlet, 0, i ); // unknown 2 - z
        }
        tmp << 0,0,-1;
    }
    //! [Refinement]

    gsConstantFunction<> force(tmp,3);
    gsFunctionExpr<> t(std::to_string(thickness), 3);
    gsFunctionExpr<> E(std::to_string(E_modulus),3);
    gsFunctionExpr<> nu(std::to_string(PoissonRatio),3);

    gsThinShellAssembler assembler(mp,dbasis,bc,force,t,E,nu);


    gsSparseSolver<>::CGDiagonal solver;

    // gsInfo<< A.numDofs() <<"\n"<<std::flush;



    gsVector<> pt(2); pt.setConstant(0.25);
    // evaluateFunction(ev, u * ff * meas(G), pt); // evaluates an expression on a point

    // ! [Solve linear problem]

    // assemble system
    assembler.assemble();
    // solve system
    solver.compute( assembler.matrix() );
    gsVector<> solVector = solver.solve(assembler.rhs());

    // update deformed patch
    // gsMatrix<> cc;
    // for ( size_t k =0; k!=mp_def.nPatches(); ++k) // Deform the geometry
    // {
    //     // extract deformed geometry
    //     u_sol.extract(cc, k);
    //     mp_def.patch(k).coefs() += cc;  // defG points to mp_def, therefore updated
    // }

    /*Something with Dirichlet homogenization*/

    // ! [Solve linear problem]

    // ! [Solve nonlinear problem]

    // real_t residual = A.rhs().norm();
    // if (nonlinear)
    // {
    //     index_t itMax = 10;
    //     real_t tol = 1e-8;
    //     for (index_t it = 0; it != itMax; ++it)
    //     {
    //         A.initSystem();
    //         // assemble system
    //         A.assemble(
    //             (
    //             (tt.val()) * (E_m_der * reshape(mm,3,3) * E_m_der.tr() + E_m_der2)
    //             +
    //             (tt.val() * tt.val() * tt.val())/3.0 * (E_f_der * reshape(mm,3,3) * E_f_der.tr() -  E_f_der2)
    //             ) * meas(G)
    //             , u * ff * meas(G)
    //             -
    //             (
    //              (
    //                 (tt.val()) *(E_m * reshape(mm,3,3) * E_m_der.tr()) -
    //                 (tt.val() * tt.val() * tt.val())/3.0 * (E_f * reshape(mm,3,3) * E_f_der.tr())
    //              ) * meas(G)
    //             ).tr()
    //             );

    //         // A.assemble(tt.val() * tt.val() * tt.val() / 3.0 * E_f_der2);
    //         // solve system
    //         solver.compute( A.matrix() );
    //         solVector = solver.solve(A.rhs()); // this is the UPDATE
    //         residual = A.rhs().norm();

    //         gsInfo<<"Iteration: "<< it
    //                <<", residue: "<< residual
    //                <<", update norm: "<<solVector.norm()
    //                <<"\n";

    //         // update deformed patch
    //         gsMatrix<> cc;
    //         for ( size_t k =0; k!=mp_def.nPatches(); ++k) // Deform the geometry
    //         {
    //             // extract deformed geometry
    //             u_sol.extract(cc, k);
    //             mp_def.patch(k).coefs() += cc;  // defG points to mp_def, therefore updated


    //             // gsInfo<<"coefficients = "<<cc<<"\n";
    //             // gsInfo<<"solVector = "<<solVector.transpose()<<"\n";
    //         }


    //         if (residual < tol)
    //             break;
    //     }
    // }


    // ! [Solve nonlinear problem]

    // For Neumann (same for Dirichlet/Nitche) conditions
    // variable g_N = A.getBdrFunction();
    // A.assembleRhsBc(u * g_N.val() * nv(G).norm(), bc.neumannSides() );

    // Penalize the matrix? (we need values for the DoFs to be enforced..
    // function/call:  penalize_matrix(DoF_indices, DoF_values)
    // otherwise: should we tag the DoFs inside "u" ?

    // gsInfo<<"RHS rows = "<<A.rhs().rows()<<"\n";
    // gsInfo<<"RHS cols = "<<A.rhs().cols()<<"\n";
    // gsInfo<<"MAT rows = "<<A.matrix().rows()<<"\n";
    // gsInfo<<"MAT cols = "<<A.matrix().cols()<<"\n";

    // gsInfo<< A.rhs().transpose() <<"\n";
    // gsInfo<< A.matrix().toDense()<<"\n";

    // ADD BOUNDARY CONDITIONS! (clamped will be tricky..............)


    //! [Export visualization in ParaView]
    // if (plot)
    // {
    //     gsMultiPatch<> deformation = mp_def;
    //     for (index_t k = 0; k != mp_def.nPatches(); ++k)
    //         deformation.patch(0).coefs() -= mp.patch(0).coefs();

    //     gsField<> solField(mp, deformation);
    //     gsInfo<<"Plotting in Paraview...\n";
    //     gsWriteParaview<>( solField, "solution", 1000, true);

    //     // ev.options().setSwitch("plot.elements", true);
    //     // ev.writeParaview( u_sol   , G, "solution");

    //     // gsFileManager::open("solution.pvd");
    // }

    return EXIT_SUCCESS;

}// end main

template <class T>
void evaluateFunction(gsExprEvaluator<T> ev, auto expression, gsVector<T> pt)
{
    gsMatrix<T> evresult = ev.eval( expression,pt );
    gsInfo << "Eval on point ("<<pt.at(0)<<" , "<<pt.at(1)<<") :\n"<< evresult;
    gsInfo << "\nEnd ("<< evresult.rows()<< " x "<<evresult.cols()<<")\n";
};

/*
    to do:
    =  make function for construction of the solution given the space and the mp
*/



/*
template<class T>
void gsShellAssembler<T>::applyLoads()
{
    gsMatrix<T>        bVals;
    gsMatrix<unsigned> acts,globalActs;

    for (size_t i = 0; i< m_pLoads.numLoads(); ++i )
    {
        if ( m_pLoads[i].parametric )
        {
            m_bases.front().basis(m_pLoads[i].patch).active_into( m_pLoads[i].point, acts );
            m_bases.front().basis(m_pLoads[i].patch).eval_into  ( m_pLoads[i].point, bVals);
        }
        else
        {
            gsMatrix<> forcePoint;
            m_patches.patch(m_pLoads[i].patch).invertPoints(m_pLoads[i].point,forcePoint);
            u.source().piece(m_pLoads[i].patch).active_into( forcePoint, acts );
            u.source().piece(m_pLoads[i].patch).active_into( forcePoint, bVals);
        }

        // translate patch-local indices to global dof indices
        for (size_t j = 0; j< 3; ++j)
        {
            if (m_pLoads[i].value[j] != 0.0)
            {
                u.dofMappers[j].localToGlobal(acts, m_pLoads[i].patch, globalActs);

                for (index_t k=0; k < globalActs.rows(); ++k)
                {
                    if (int(globalActs(k,0)) < m_dofs)
                        m_rhs(globalActs(k,0), 0) += bVals(k,0) * m_pLoads[i].value[j];
                }
            }
        }
    }
}
*/
