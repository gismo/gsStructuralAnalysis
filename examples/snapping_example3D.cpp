/** @file bSplineCurve_example.cpp

    @brief Tutorial on gsBSpline class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <iostream>

#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMBase.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMLoadControl.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMRiks.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMCrisfield.h>

#ifdef gsElasticity_ENABLED
#include <gsElasticity/gsGeoUtils.h>
#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>

#include <gsElasticity/gsNeoHookLogMaterial.h>
#include <gsElasticity/gsLinearMaterial.h>
#endif

#include <gismo.h>

using namespace gismo;

#ifdef gsElasticity_ENABLED
template <class T>
std::vector<gsBSpline<T>> makeCurve(const T tw, const T tg, const T tb, const T ts, const T l, const T a, const std::string expr, const gsKnotVector<T> & kv1);

template <class T>
gsMultiPatch<T> makeTop(const T tw, const T tg, const T tb, const T ts, const T l, const T a, const std::string expr, const gsKnotVector<T> & kv1);
template <class T>
gsMultiPatch<T> makeTop(const T tw, const T tg, const T tb, const T ts, const T l, const T a, const std::vector<gsBSpline<T>> & curves);

template <class T>
gsMultiPatch<T> makeElement(const T tw, const T tg, const T tb, const T ts, const T l, const T a, const std::string expr, const gsKnotVector<T> & kv1);
template <class T>
gsMultiPatch<T> makeElement(const T tw, const T tg, const T tb, const T ts, const T l, const T a, const std::vector<gsBSpline<T>> & curves);

template <class T>
gsMultiPatch<T> makeBottom(const T tw, const T tg, const T tb, const T ts, const T l, const T a, const std::string expr, const gsKnotVector<T> & kv1);
template <class T>
gsMultiPatch<T> makeBottom(const T tw, const T tg, const T tb, const T ts, const T l, const T a, const std::vector<gsBSpline<T>> & curves);

int main(int argc, char *argv[])
{
    // Input options
    int numElevate  = 0;
    int numHref     = 0;
    bool plot       = false;
    bool write = false;

    // Arc length method options
    real_t dL = 0; // General arc length
    real_t dLb = 0.1; // Ard length to find bifurcation
    real_t tol = 1e-6;
    real_t tolU = 1e-6;
    real_t tolF = 1e-3;
    real_t relax = 1.0;
    real_t tau = 1e4;
    int method = 2; // (0: Load control; 1: Riks' method; 2: Crisfield's method; 3: consistent crisfield method; 4: extended iterations)
    bool SingularPoint = false;
    bool quasiNewton = false;
    int quasiNewtonInt = -1;
    bool adaptive = false;
    int step = 10;
    index_t maxit = 20;

    std::string wn("data.csv");

    std::string output = "solutionElasticity";

    // geometry options
    real_t start = 0; // starting knot
    real_t end = 1; // ending knot
    index_t interior = 4; // number of interior knots
    index_t multEnd = 3; // multiplicity at the two end knots

    index_t Nx = 1;
    index_t Ny = 0;

    real_t tw = 1.5e-3;
    real_t tg = 1.0e-3;
    real_t tb = 1.5e-3;
    real_t ts = 1.0e-3;
    real_t l  = 10.0e-3;
    real_t al = 0.1;
    real_t b = 3e-3;

    real_t Emax = 1.5;

    gsCmdLine cmd("");
    cmd.addInt("n","interior","Number of interior knots",interior);
    cmd.addInt("X","Nx","Number of element in x-direction",Nx);
    cmd.addInt("Y","Ny","Number of element in y-direction",Ny);

    cmd.addInt("r","hRefine",
     "Number of dyadic h-refinement (bisection) steps to perform before solving",
     numHref);
    cmd.addInt("e","degreeElevation",
      "Number of degree elevation steps to perform on the Geometry's basis before solving",
      numElevate);

    cmd.addReal("a","a/l", "Value of a/l", al);
    cmd.addReal("E","Emax", "Maximum strain", Emax);

    cmd.addInt("m","Method", "Arc length method; 1: Crisfield's method; 2: RIks' method.", method);
    cmd.addReal("L","dLb", "arc length", dLb);
    cmd.addReal("l","dL", "arc length after bifurcation", dL);
    cmd.addReal("A","relaxation", "Relaxation factor for arc length method", relax);

    cmd.addReal("P","perturbation", "perturbation factor", tau);

    cmd.addInt("q","QuasiNewtonInt","Use the Quasi Newton method every INT iterations",quasiNewtonInt);
    cmd.addInt("N", "maxsteps", "Maximum number of steps", step);

    cmd.addSwitch("adaptive", "Adaptive length ", adaptive);
    cmd.addSwitch("bifurcation", "Compute singular points and bifurcation paths", SingularPoint);
    cmd.addSwitch("quasi", "Use the Quasi Newton method", quasiNewton);

    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addSwitch("write", "Write convergence data to file", write);


    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    GISMO_ENSURE(Nx>0,"Number of elements in x must be larger than 0");
    GISMO_ENSURE(Ny>=0,"Number of elements in x must be positive");

    real_t a  = al * l / 2.;

    gsKnotVector<> kv1(start, end, interior, multEnd);
    std::string expr = "-" + std::to_string(a) + "*cos(2*x*pi)";

    gsInfo<<"-----------------------------------------------------------------------------------------------------\n";
    gsInfo<<"-----------------------------------Making the curve--------------------------------------------------\n";
    gsInfo<<"-----------------------------------------------------------------------------------------------------\n";
    std::vector<gsBSpline<>> curves = makeCurve(tw,tg,tb,ts,l,a, expr, kv1);

    gsInfo<<"-----------------------------------------------------------------------------------------------------\n";
    gsInfo<<"-----------------------------------Making the element------------------------------------------------\n";
    gsInfo<<"-----------------------------------------------------------------------------------------------------\n";
    gsMultiPatch<> element = makeElement(tw,tg,tb,ts,l,a,curves);

    gsWriteParaview(element,"element",1000,true);

    gsInfo<<"-----------------------------------------------------------------------------------------------------\n";
    gsInfo<<"-----------------------------------Making the bottom block-------------------------------------------\n";
    gsInfo<<"-----------------------------------------------------------------------------------------------------\n";
    gsMultiPatch<> bottom = makeBottom(tw,tg,tb,ts,l,a,curves);

    gsWriteParaview(bottom,"bottom",1000,true);

    gsInfo<<"-----------------------------------------------------------------------------------------------------\n";
    gsInfo<<"-----------------------------------Making the top block-------------------------------------------\n";
    gsInfo<<"-----------------------------------------------------------------------------------------------------\n";
    gsMultiPatch<> top = makeTop(tw,tg,tb,ts,l,a,curves);

    gsWriteParaview(top,"top",1000,true);

    gsInfo<<"-----------------------------------------------------------------------------------------------------\n";
    gsInfo<<"-----------------------------------Making the geometry--------------------------------------------\n";
    gsInfo<<"-----------------------------------------------------------------------------------------------------\n";
    gsMultiPatch<> mp;
    real_t dx = l;
    real_t dy = 2*tg+ts+tb;
    index_t topmid_ID;
    for (index_t kx = 0; kx!=Nx; kx++)
    {
        gsMultiPatch<> tmp;

        // bottom
        tmp = bottom;
        gsNurbsCreator<>::shift2D(tmp,kx*dx,0);
        for (size_t p = 0; p!=tmp.nPatches(); p ++)
        {
            gsTensorBSpline<3,real_t>::uPtr tmp2 = gsNurbsCreator<>::lift3D(dynamic_cast<gsTensorBSpline<2,real_t> &>(tmp.patch(p)),b);
            mp.addPatch(*tmp2,tmp.patch(p).label());
        }

        // mid
        for (index_t ky = 0; ky!=Ny; ky++)
        {
            tmp = element;
            gsNurbsCreator<>::shift2D(tmp,kx*dx,ky*dy);
            for (size_t p = 0; p!=tmp.nPatches(); p ++)
            {
                gsTensorBSpline<3,real_t>::uPtr tmp2 = gsNurbsCreator<>::lift3D(dynamic_cast<gsTensorBSpline<2,real_t> &>(tmp.patch(p)),b);
                mp.addPatch(*tmp2,tmp.patch(p).label());
            }
        }

        // top
        tmp = top;
        if (kx == math::floor(Nx/2.))
        {
            if (Nx % 2 ==0)
                topmid_ID = mp.nPatches()+8;
            else
                topmid_ID = mp.nPatches()+10;
        }
        gsNurbsCreator<>::shift2D(tmp,kx*dx,Ny*dy);
        for (size_t p = 0; p!=tmp.nPatches(); p ++)
        {
            gsTensorBSpline<3,real_t>::uPtr tmp2 = gsNurbsCreator<>::lift3D(dynamic_cast<gsTensorBSpline<2,real_t> &>(tmp.patch(p)),b);
            mp.addPatch(*tmp2,tmp.patch(p).label());
        }
    }
    mp.computeTopology();

    if (plot)
    {
        gsWriteParaview<>(mp,"mp",1000,true,false);
    }

    gsMultiPatch<> mp_def = mp;
    gsMultiPatch<> deformation;

    // define basis
    gsMultiBasis<> dbasis(mp);

    // Boundary conditions
    gsBoundaryConditions<> bc;

    gsVector<> neuData(3);
    neuData<<0,1,0;
    gsConstantFunction<> neu(neuData,3);

    gsMatrix<> bbox;
    mp.boundingBox(bbox);
    real_t band = 0.01 * (bbox(1,1) - bbox(1,0));
    real_t ymin = bbox(1,0);
    real_t ymax = bbox(1,1);
    for (typename gsBoxTopology::const_biterator bit = mp.bBegin(); bit!=mp.bEnd(); bit++)
    {
        // get the parametric center of the boundary
        gsMatrix<> pRange = mp.patch(bit->patch).boundary(bit->side())->support();
        real_t center = (pRange(0,1)-pRange(0,0)) / 2. + pRange(0,0);
        gsMatrix<> centerpt(2,1);
        centerpt<<center,0.5;
        gsMatrix<> eval = mp.patch(bit->patch).boundary(bit->side())->eval(centerpt);

        if ( (eval(1,0) <= ymin + band) && (eval(1,0) >= ymin - band))
        {
            bc.addCondition(bit->patch,bit->side(),condition_type::dirichlet,0,0);
            bc.addCondition(bit->patch,bit->side(),condition_type::dirichlet,0,1);
            bc.addCondition(bit->patch,bit->side(),condition_type::dirichlet,0,2);
        }
        else if ((eval(1,0) <= ymax + band) && (eval(1,0) >= ymax - band))
        {
            bc.addCondition(bit->patch,bit->side(),condition_type::dirichlet,0,0);
            bc.addCondition(bit->patch,bit->side(),condition_type::collapsed,0,1);
            bc.addCondition(bit->patch,bit->side(),condition_type::dirichlet,0,2);
            bc.addCondition(bit->patch,bit->side(),condition_type::neumann,&neu);
        }

    }

    // fd.getId(20, bc); // id=2: boundary conditions
    bc.setGeoMap(mp);
    gsInfo<<"Boundary conditions:\n"<< bc <<"\n";

    // Loads
    gsFunctionExpr<> bodyForce("0","0","0",3);
    gsInfo<<"Body force function "<< bodyForce << "\n";

    for (index_t i = 0; i < numElevate; ++i)
    {
        dbasis.degreeElevate();
        dbasis.uniformRefine(1,1,0);
        dbasis.uniformRefine(1,1,1);
    }
    for (index_t i = 0; i < numHref; ++i)
    {
        dbasis.uniformRefine(1,1,0);
        dbasis.uniformRefine(1,1,1);
    }

    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
    gsInfo<<mp<<"\n";

    // -----------------------------------------------------------------------------------------------
    // --------------------------------------Solve Static problem-------------------------------------
    // -----------------------------------------------------------------------------------------------

    real_t E = 78e6;
    real_t nu = 0.4;
    // gsLinearMaterial<real_t> material(E,nu,&mp,&mp_def);

    gsElasticityAssembler<real_t> assembler(mp,dbasis,bc,bodyForce);//,&material);
    gsInfo << "Initialized system with " << assembler.numDofs() << " dofs.\n";

    assembler.options().setReal("YoungsModulus",E);
    assembler.options().setReal("PoissonsRatio",nu);
    assembler.options().setInt("MaterialLaw",material_law::hooke);

    gsStopwatch stopwatch,stopwatch2;
    // Define Matrices
    stopwatch.restart();
    stopwatch2.restart();
    real_t time = 0.0;

    assembler.assemble();
    time += stopwatch.stop();
    gsVector<> Force = assembler.rhs();

    std::vector<gsMatrix<> > fixedDofs = assembler.allFixedDofs();
    // Function for the Jacobian
    gsStructuralAnalysisOps<real_t>::Jacobian_t Jacobian = [&time,&stopwatch,&assembler,&fixedDofs](gsVector<real_t> const &x, gsSparseMatrix<real_t> & m)
    {
        stopwatch.restart();
        assembler.assemble(x,fixedDofs);
        time += stopwatch.stop();

        m = assembler.matrix();
        // gsInfo<<"matrix = \n"<<m.toDense()<<"\n";
        return true;
    };
    // Function for the Residual
    gsStructuralAnalysisOps<real_t>::ALResidual_t ALResidual = [&time,&stopwatch,&assembler,&fixedDofs,&Force](gsVector<real_t> const &x, real_t lam, gsVector<real_t> & result)
    {
        stopwatch.restart();
        assembler.assemble(x,fixedDofs);
        result = Force - lam * Force - assembler.rhs(); // assembler rhs - force = Finternal
        // gsInfo<<"vector = \n"<<result.transpose()<<"\n";
        time += stopwatch.stop();
        return true;
    };

    assembler.options().setInt("MaterialLaw",material_law::neo_hooke_quad);
    
    //=============================================//
                  // Solving //
    //=============================================//

  gsALMBase<real_t> * arcLength;
  if (method==0)
    arcLength = new gsALMLoadControl<real_t>(Jacobian, ALResidual, Force);
else if (method==1)
    arcLength = new gsALMRiks<real_t>(Jacobian, ALResidual, Force);
else if (method==2)
    arcLength = new gsALMCrisfield<real_t>(Jacobian, ALResidual, Force);
else
    GISMO_ERROR("Method "<<method<<" unknown");

      arcLength->options().setString("Solver","SimplicialLDLT"); // LDLT solver
      arcLength->options().setInt("BifurcationMethod",0); // 0: determinant, 1: eigenvalue
      arcLength->options().setReal("Length",dLb);
      arcLength->options().setInt("AngleMethod",0); // 0: step, 1: iteration
      arcLength->options().setSwitch("AdaptiveLength",adaptive);
      arcLength->options().setInt("AdaptiveIterations",5);
      arcLength->options().setReal("Perturbation",tau);
      arcLength->options().setReal("Scaling",0.0);
      arcLength->options().setReal("Tol",tol);
      arcLength->options().setReal("TolU",tolU);
      arcLength->options().setReal("TolF",tolF);
      arcLength->options().setInt("MaxIter",maxit);
      arcLength->options().setSwitch("Verbose",true);
      arcLength->options().setReal("Relaxation",relax);
      if (quasiNewtonInt>0)
      {
        quasiNewton = true;
        arcLength->options().setInt("QuasiIterations",quasiNewtonInt);
    }
    arcLength->options().setSwitch("Quasi",quasiNewton);
    arcLength->options().setReal("SingularPointComputeTolB",1e-1);
    arcLength->options().setReal("SingularPointComputeTolE",1e-10);

    gsDebug<<arcLength->options();
    arcLength->applyOptions();
    arcLength->initialize();

    gsMatrix<> writePoints(3,1);
    gsVector<index_t> writePatches(1);
    writePatches<<topmid_ID;
    real_t xBoundary;
    if (Nx % 2 == 0)
    {
        gsMatrix<> supp = mp.patch(writePatches[0]).boundary(boundary::north)->support();
        xBoundary = supp(0,0);
    }
    else
    {
        gsMatrix<> northCenter = mp.patch(writePatches[0]).boundary(boundary::north)->parameterCenter();
        xBoundary = northCenter(0,0);
    }
    gsMatrix<> supp = mp.patch(writePatches[0]).support();
    writePoints.col(0)<<xBoundary,supp(1,1),0.5;

    gsMatrix<> evalPoints;
    mp.patch(topmid_ID).eval_into(writePoints,evalPoints);

    std::string dirname = "ArcLengthResults/snapping_3D_al=" + std::to_string(al);

    gsParaviewCollection collection(dirname + "/" + output);
    deformation = mp;

    // Make objects for previous solutions
    real_t Lold = 0;
    gsMatrix<> Uold = Force;
    Uold.setZero();

    gsMatrix<> solVector;
    real_t indicator = 0.0;
    arcLength->setIndicator(indicator); // RESET INDICATOR
    bool bisected = false;
    real_t dLb0 = dLb;

    if (write)
    {
        std::ofstream file;
        file.open(dirname + "/" + wn,std::ofstream::out);
        file    << std::setprecision(20)
        << "Deformation norm" << ",";
        for (index_t k=0; k!=writePatches.size(); k++)
        {
            file<< "point "<<k<<" - x" << ","
            << "point "<<k<<" - y" << ",";

        }
        file  << "Eps" << ",";
        file  << "Lambda" << ","
        << "Indicator"
        << "\n";
        file.close();
    }

    index_t  k = 0;
    real_t eps = 0;
    while (eps<=Emax && k < step)
    {

        gsInfo<<"Load step "<< k<<"\n";
        arcLength->step();

        if (!(arcLength->converged()))
        {
            gsInfo<<"Error: Loop terminated, arc length method did not converge.\n";
            dLb = dLb / 2.;
            arcLength->setLength(dLb);
            arcLength->setSolution(Uold,Lold);
            bisected = true;
            continue;
        }

        if (SingularPoint)
        {
            arcLength->computeStability(quasiNewton);
            if (arcLength->stabilityChange())
            {
                gsInfo<<"Bifurcation spotted!"<<"\n";
                arcLength->computeSingularPoint(false);
                arcLength->switchBranch();
                dLb0 = dLb = dL;
                arcLength->setLength(dLb);
            }
        }
        indicator = arcLength->indicator();

        solVector = arcLength->solutionU();
        Uold = solVector;
        Lold = arcLength->solutionL();
        assembler.constructSolution(solVector,fixedDofs,mp_def);
    // constructing an IGA field (geometry + solution)
        gsField<> solutionField(assembler.patches(),mp_def);
    // constructing stress tensor
        gsPiecewiseFunction<> stresses;
        assembler.constructCauchyStresses(mp_def,stresses,stress_components::von_mises);
        gsField<> stressField(assembler.patches(),stresses,true);

        gsInfo<<"Total ellapsed assembly time: "<<time<<" s\n";

        if (plot)
        {
            std::string fileName = dirname + "/" + output + util::to_string(k) + "_";
            // creating a container to plot all fields to one Paraview file
            std::map<std::string,const gsField<> *> fields;
            fields["Deformation"] = &solutionField;
            fields["Stress"] = &stressField;
            gsWriteParaviewMultiPhysics(fields,fileName,200,true);
            fileName = gsFileManager::getFilename(fileName);
            for (size_t p=0; p!=mp.nPatches(); p++)
            {
                collection.addPart(fileName + std::to_string(p) + ".vts",k,"",p);
                collection.addPart(fileName + std::to_string(p) + "_mesh.vtp",k,"",p);
            }
        }

        gsMatrix<> out_def = mp_def.patch(writePatches(0)).eval(writePoints.col(0));
        gsMatrix<> out_ori = mp    .patch(writePatches(0)).eval(writePoints.col(0));
        eps = out_def(1,0) / out_ori(1,0);
        gsDebugVar(eps);
        gsMatrix<> out = mp_def.patch(writePatches(0)).eval(writePoints.col(0));
        if (write)
        {

            std::ofstream file;
            file.open(dirname + "/" + wn,std::ofstream::out | std::ofstream::app);
            file    << std::setprecision(6)
                    << arcLength->solutionU().norm() << ",";
            file<< out(0,0) << ","
                << out(1,0) << ",";
            file<< eps<<",";
            file    << arcLength->solutionL() << "," << arcLength->indicator()<<"\n";
            file.close();
        }

        if (!bisected)
        {
          dLb = dLb0;
          arcLength->setLength(dLb);
      }
      bisected = false;
      k++;
    }

    if (plot)
    {
        collection.save();
    }

    delete arcLength;

    return 1;
}

template <class T>
std::vector<gsBSpline<T>> makeCurve(const T tw, const T tg, const T tb, const T ts, const T l, const T a, const std::string expr, const gsKnotVector<T> & kv1)
{

    gsFunctionExpr<T> fun(expr,2);
    gsKnotVector<T> kv2(0, 1, 0, 2);
    gsBSplineBasis<real_t> basis1(kv1);
    gsTensorBSplineBasis<2,real_t> basis2(kv1,kv2);

    gsMultiPatch<T> mp;
    mp.addPatch(gsNurbsCreator<T>::BSplineSquare());
    mp.patch(0).coefs()<<   0 , 0 , 1 , 0
    , 0 , 1 , 1 , 1 ;

    gsBoundaryConditions<T> bc;
    bc.addCondition(boundary::west,condition_type::dirichlet,&fun,0);
    bc.addCondition(boundary::east,condition_type::dirichlet,&fun,0);
    bc.setGeoMap(mp);

    gsMultiBasis<T> mb(basis2);
    gsMatrix<T> solution;

    gsExprAssembler<T> A(1,1);
    A.setIntegrationElements(mb);
    auto u = A.getSpace(mb,1);
    auto G = A.getMap(mp);
    auto f = A.getCoeff(fun,G);
    auto sol=A.getSolution(u,solution);

    u.setup(bc,dirichlet::l2Projection,0);
    A.initSystem();
    A.assemble(u*u.tr()*meas(G),u*f*meas(G));
    typename gsSparseSolver<T>::LU solver;
    solver.compute(A.matrix());
    solution = solver.solve(A.rhs());

    gsMatrix<T> spline_coefs;
    sol.extract(spline_coefs);

    gsTensorBSpline<2,real_t> spline(basis2,spline_coefs);
    gsWriteParaview<T>(spline,"spline",1000,true,true);

    gsTensorBSpline<2,real_t> left, midleft, mid, midright, right;
    real_t splitleft    = tw/2. / l;
    real_t splitmidleft = 0.5 - tw/2. / l;
    real_t splitmidright= 0.5 + tw/2. / l;
    real_t splitright   = 1 - tw/2. / l;

    spline.splitAt(0,splitleft,left,spline);
    spline.splitAt(0,splitmidleft,midleft,spline);
    spline.splitAt(0,splitmidright,mid,spline);
    spline.splitAt(0,splitright,midright,right);

    gsBSpline<T> left_crv, midleft_crv, mid_crv, midright_crv, right_crv;
    left.slice(1,0,left_crv);
    midleft.slice(1,0,midleft_crv);
    mid.slice(1,0,mid_crv);
    midright.slice(1,0,midright_crv);
    right.slice(1,0,right_crv);

    std::vector<gsBSpline<T>> result;
    result.push_back(left_crv);
    result.push_back(midleft_crv);
    result.push_back(mid_crv);
    result.push_back(midright_crv);
    result.push_back(right_crv);
    return result;
}

template <class T>
gsMultiPatch<T> makeTop(const T tw, const T tg, const T tb, const T ts, const T l, const T a, const std::string expr, const gsKnotVector<T> & kv1)
{
    std::vector<gsBSpline<T>> curves = makeCurve(tw,tg,tb,ts,l,a, expr, kv1);
    return makeTop(tw,tg,tb,ts,l,a,curves);
}

template <class T>
gsMultiPatch<T> makeTop(const T tw, const T tg, const T tb, const T ts, const T l, const T a, const std::vector<gsBSpline<T>> & curves)
{
    gsKnotVector<T> kv2(0, 1, 0, 2);

    gsBSpline<T> left_crv = curves[0];
    gsBSpline<T> midleft_crv = curves[1];
    gsBSpline<T> mid_crv = curves[2];
    gsBSpline<T> midright_crv = curves[3];
    gsBSpline<T> right_crv = curves[4];

    gsBSplineBasis<T> left_basis, midleft_basis, mid_basis, midright_basis, right_basis;
    left_basis = left_crv.basis();
    midleft_basis = midleft_crv.basis();
    mid_basis = mid_crv.basis();
    midright_basis = midright_crv.basis();
    right_basis = right_crv.basis();

    gsMatrix<T> left_coefs, midleft_coefs, mid_coefs, midright_coefs, right_coefs;
    gsMatrix<T> tmp_coefs;

    gsTensorBSplineBasis<2,real_t> left_tbasis(left_basis.knots(),kv2);
    gsTensorBSplineBasis<2,real_t> midleft_tbasis(midleft_basis.knots(),kv2);
    gsTensorBSplineBasis<2,real_t> mid_tbasis(mid_basis.knots(),kv2);
    gsTensorBSplineBasis<2,real_t> midright_tbasis(midright_basis.knots(),kv2);
    gsTensorBSplineBasis<2,real_t> right_tbasis(right_basis.knots(),kv2);

    gsMatrix<T> coefs_tmp, tmp, ones;
    gsMultiPatch<T> top;
    // Get the coefficients
    tmp_coefs.resize(left_crv.coefs().rows(),2);
    tmp_coefs.col(0) = left_basis.anchors().transpose() * l;
    tmp_coefs.col(1) = left_crv.coefs();
    left_coefs = tmp_coefs;

    tmp_coefs.resize(midleft_crv.coefs().rows(),2);
    tmp_coefs.col(0) = midleft_basis.anchors().transpose() * l;
    tmp_coefs.col(1) = midleft_crv.coefs();
    midleft_coefs = tmp_coefs;

    tmp_coefs.resize(mid_crv.coefs().rows(),2);
    tmp_coefs.col(0) = mid_basis.anchors().transpose() * l;
    tmp_coefs.col(1) = mid_crv.coefs();
    mid_coefs = tmp_coefs;

    tmp_coefs.resize(midright_crv.coefs().rows(),2);
    tmp_coefs.col(0) = midright_basis.anchors().transpose() * l;
    tmp_coefs.col(1) = midright_crv.coefs();
    midright_coefs = tmp_coefs;

    tmp_coefs.resize(right_crv.coefs().rows(),2);
    tmp_coefs.col(0) = right_basis.anchors().transpose() * l;
    tmp_coefs.col(1) = right_crv.coefs();
    right_coefs = tmp_coefs;

    // The cosine starts at a, but we want it to start at tg/2
    ones = gsMatrix<T>::Ones(left_coefs.rows(),1);
    left_coefs.col(1) += (tg/2.-a) * ones;
    ones = gsMatrix<T>::Ones(midleft_coefs.rows(),1);
    midleft_coefs.col(1) += (tg/2.-a) * ones;
    ones = gsMatrix<T>::Ones(mid_coefs.rows(),1);
    mid_coefs.col(1) += (tg/2.-a) * ones;
    ones = gsMatrix<T>::Ones(midright_coefs.rows(),1);
    midright_coefs.col(1) += (tg/2.-a) * ones;
    ones = gsMatrix<T>::Ones(right_coefs.rows(),1);
    right_coefs.col(1) += (tg/2.-a) * ones;

    // we start at the bottom block
    tmp = mid_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1) -= tg * ones;
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    top.addPatch(gsTensorBSpline<2,real_t>(mid_tbasis,coefs_tmp));

    // Now the bottom cosine
    tmp = left_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1) += ts * ones;
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    top.addPatch(gsTensorBSpline<2,real_t>(left_tbasis,coefs_tmp));

    tmp = midleft_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1) += ts * ones;
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    top.addPatch(gsTensorBSpline<2,real_t>(midleft_tbasis,coefs_tmp));

    tmp = mid_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1) += ts * ones;
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    top.addPatch(gsTensorBSpline<2,real_t>(mid_tbasis,coefs_tmp));

    tmp = midright_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1) += ts * ones;
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    top.addPatch(gsTensorBSpline<2,real_t>(midright_tbasis,coefs_tmp));

    tmp = right_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1) += ts * ones;
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    top.addPatch(gsTensorBSpline<2,real_t>(right_tbasis,coefs_tmp));

    // Lift all coefficient matrices with ts
    ones = gsMatrix<T>::Ones(left_coefs.rows(),1);
    left_coefs.col(1) += ts * ones;
    ones = gsMatrix<T>::Ones(midleft_coefs.rows(),1);
    midleft_coefs.col(1) += ts * ones;
    ones = gsMatrix<T>::Ones(mid_coefs.rows(),1);
    mid_coefs.col(1) += ts * ones;
    ones = gsMatrix<T>::Ones(midright_coefs.rows(),1);
    midright_coefs.col(1) += ts * ones;
    ones = gsMatrix<T>::Ones(right_coefs.rows(),1);
    right_coefs.col(1) += ts * ones;

    // mid-row
    tmp = left_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1) += tg*ones;
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    top.addPatch(gsTensorBSpline<2,real_t>(left_tbasis,coefs_tmp));

    tmp = right_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1) += tg*ones;
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    top.addPatch(gsTensorBSpline<2,real_t>(right_tbasis,coefs_tmp));

    // Lift all coefficient matrices with tg/2.
    ones = gsMatrix<T>::Ones(left_coefs.rows(),1);
    left_coefs.col(1) += tg * ones;
    ones = gsMatrix<T>::Ones(midleft_coefs.rows(),1);
    midleft_coefs.col(1) += tg * ones;
    ones = gsMatrix<T>::Ones(mid_coefs.rows(),1);
    mid_coefs.col(1) += tg * ones;
    ones = gsMatrix<T>::Ones(midright_coefs.rows(),1);
    midright_coefs.col(1) += tg * ones;
    ones = gsMatrix<T>::Ones(right_coefs.rows(),1);
    right_coefs.col(1) += tg * ones;

    // Now the top cosine
    tmp = left_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1).setConstant(tg/2+ts+tg+tb+(5*tg));
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    top.addPatch(gsTensorBSpline<2,real_t>(left_tbasis,coefs_tmp));

    tmp = midleft_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1).setConstant(tg/2+ts+tg+tb+(5*tg));
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    top.addPatch(gsTensorBSpline<2,real_t>(midleft_tbasis,coefs_tmp));

    tmp = mid_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1).setConstant(tg/2+ts+tg+tb+(5*tg));
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    top.addPatch(gsTensorBSpline<2,real_t>(mid_tbasis,coefs_tmp));

    tmp = midright_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1).setConstant(tg/2+ts+tg+tb+(5*tg));
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    top.addPatch(gsTensorBSpline<2,real_t>(midright_tbasis,coefs_tmp));

    tmp = right_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1).setConstant(tg/2+ts+tg+tb+(5*tg));
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    top.addPatch(gsTensorBSpline<2,real_t>(right_tbasis,coefs_tmp));

    for (index_t p=8; p!=13; p++)
        top.patch(p).uniformRefine(1,1,1);

    return top;
}

template <class T>
gsMultiPatch<T> makeElement(const T tw, const T tg, const T tb, const T ts, const T l, const T a, const std::string expr, const gsKnotVector<T> & kv1)
{
    std::vector<gsBSpline<T>> curves = makeCurve(tw,tg,tb,ts,l,a, expr, kv1);
    return makeElement(tw,tg,tb,ts,l,a,curves);
}

template <class T>
gsMultiPatch<T> makeElement(const T tw, const T tg, const T tb, const T ts, const T l, const T a, const std::vector<gsBSpline<T>> & curves)
{
    gsKnotVector<T> kv2(0, 1, 0, 2);

    gsBSpline<T> left_crv = curves[0];
    gsBSpline<T> midleft_crv = curves[1];
    gsBSpline<T> mid_crv = curves[2];
    gsBSpline<T> midright_crv = curves[3];
    gsBSpline<T> right_crv = curves[4];

    gsBSplineBasis<T> left_basis, midleft_basis, mid_basis, midright_basis, right_basis;
    left_basis = left_crv.basis();
    midleft_basis = midleft_crv.basis();
    mid_basis = mid_crv.basis();
    midright_basis = midright_crv.basis();
    right_basis = right_crv.basis();

    gsMatrix<T> left_coefs, midleft_coefs, mid_coefs, midright_coefs, right_coefs;
    gsMatrix<T> tmp_coefs;

    gsTensorBSplineBasis<2,real_t> left_tbasis(left_basis.knots(),kv2);
    gsTensorBSplineBasis<2,real_t> midleft_tbasis(midleft_basis.knots(),kv2);
    gsTensorBSplineBasis<2,real_t> mid_tbasis(mid_basis.knots(),kv2);
    gsTensorBSplineBasis<2,real_t> midright_tbasis(midright_basis.knots(),kv2);
    gsTensorBSplineBasis<2,real_t> right_tbasis(right_basis.knots(),kv2);

    gsMatrix<T> coefs_tmp, tmp, ones;
    gsMultiPatch<T> element;
    // Get the coefficients
    tmp_coefs.resize(left_crv.coefs().rows(),2);
    tmp_coefs.col(0) = left_basis.anchors().transpose() * l;
    tmp_coefs.col(1) = left_crv.coefs();
    left_coefs = tmp_coefs;

    tmp_coefs.resize(midleft_crv.coefs().rows(),2);
    tmp_coefs.col(0) = midleft_basis.anchors().transpose() * l;
    tmp_coefs.col(1) = midleft_crv.coefs();
    midleft_coefs = tmp_coefs;

    tmp_coefs.resize(mid_crv.coefs().rows(),2);
    tmp_coefs.col(0) = mid_basis.anchors().transpose() * l;
    tmp_coefs.col(1) = mid_crv.coefs();
    mid_coefs = tmp_coefs;

    tmp_coefs.resize(midright_crv.coefs().rows(),2);
    tmp_coefs.col(0) = midright_basis.anchors().transpose() * l;
    tmp_coefs.col(1) = midright_crv.coefs();
    midright_coefs = tmp_coefs;

    tmp_coefs.resize(right_crv.coefs().rows(),2);
    tmp_coefs.col(0) = right_basis.anchors().transpose() * l;
    tmp_coefs.col(1) = right_crv.coefs();
    right_coefs = tmp_coefs;

    // The cosine starts at a, but we want it to start at tg/2
    ones = gsMatrix<T>::Ones(left_coefs.rows(),1);
    left_coefs.col(1) += (tg/2.-a) * ones;
    ones = gsMatrix<T>::Ones(midleft_coefs.rows(),1);
    midleft_coefs.col(1) += (tg/2.-a) * ones;
    ones = gsMatrix<T>::Ones(mid_coefs.rows(),1);
    mid_coefs.col(1) += (tg/2.-a) * ones;
    ones = gsMatrix<T>::Ones(midright_coefs.rows(),1);
    midright_coefs.col(1) += (tg/2.-a) * ones;
    ones = gsMatrix<T>::Ones(right_coefs.rows(),1);
    right_coefs.col(1) += (tg/2.-a) * ones;

    // we start at the bottom two blocks
    tmp = mid_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1) -= tg*ones;
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    element.addPatch(gsTensorBSpline<2,real_t>(mid_tbasis,coefs_tmp));

    // Now the bottom cosine
    tmp = left_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1) += ts * ones;
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    element.addPatch(gsTensorBSpline<2,real_t>(left_tbasis,coefs_tmp));

    tmp = midleft_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1) += ts * ones;
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    element.addPatch(gsTensorBSpline<2,real_t>(midleft_tbasis,coefs_tmp));

    tmp = mid_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1) += ts * ones;
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    element.addPatch(gsTensorBSpline<2,real_t>(mid_tbasis,coefs_tmp));

    tmp = midright_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1) += ts * ones;
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    element.addPatch(gsTensorBSpline<2,real_t>(midright_tbasis,coefs_tmp));

    tmp = right_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1) += ts * ones;
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    element.addPatch(gsTensorBSpline<2,real_t>(right_tbasis,coefs_tmp));

    // Lift all coefficient matrices with ts
    ones = gsMatrix<T>::Ones(left_coefs.rows(),1);
    left_coefs.col(1) += ts * ones;
    ones = gsMatrix<T>::Ones(midleft_coefs.rows(),1);
    midleft_coefs.col(1) += ts * ones;
    ones = gsMatrix<T>::Ones(mid_coefs.rows(),1);
    mid_coefs.col(1) += ts * ones;
    ones = gsMatrix<T>::Ones(midright_coefs.rows(),1);
    midright_coefs.col(1) += ts * ones;
    ones = gsMatrix<T>::Ones(right_coefs.rows(),1);
    right_coefs.col(1) += ts * ones;

    // mid-row
    tmp = left_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1) += tg*ones;
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    element.addPatch(gsTensorBSpline<2,real_t>(left_tbasis,coefs_tmp));

    tmp = right_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1) += tg*ones;
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    element.addPatch(gsTensorBSpline<2,real_t>(right_tbasis,coefs_tmp));

    // Lift all coefficient matrices with tg.
    ones = gsMatrix<T>::Ones(left_coefs.rows(),1);
    left_coefs.col(1) += tg * ones;
    ones = gsMatrix<T>::Ones(midleft_coefs.rows(),1);
    midleft_coefs.col(1) += tg * ones;
    ones = gsMatrix<T>::Ones(mid_coefs.rows(),1);
    mid_coefs.col(1) += tg * ones;
    ones = gsMatrix<T>::Ones(midright_coefs.rows(),1);
    midright_coefs.col(1) += tg * ones;
    ones = gsMatrix<T>::Ones(right_coefs.rows(),1);
    right_coefs.col(1) += tg * ones;

    // Now the top cosine
    tmp = left_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1) += tb * ones;
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    element.addPatch(gsTensorBSpline<2,real_t>(left_tbasis,coefs_tmp));

    tmp = midleft_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1) += tb * ones;
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    element.addPatch(gsTensorBSpline<2,real_t>(midleft_tbasis,coefs_tmp));

    tmp = mid_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1) += tb * ones;
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    element.addPatch(gsTensorBSpline<2,real_t>(mid_tbasis,coefs_tmp));

    tmp = midright_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1) += tb * ones;
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    element.addPatch(gsTensorBSpline<2,real_t>(midright_tbasis,coefs_tmp));

    tmp = right_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1) += tb * ones;
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    element.addPatch(gsTensorBSpline<2,real_t>(right_tbasis,coefs_tmp));

    // Lift all coefficient matrices with tb
    ones = gsMatrix<T>::Ones(left_coefs.rows(),1);
    left_coefs.col(1) += tb * ones;
    ones = gsMatrix<T>::Ones(midleft_coefs.rows(),1);
    midleft_coefs.col(1) += tb * ones;
    ones = gsMatrix<T>::Ones(mid_coefs.rows(),1);
    mid_coefs.col(1) += tb * ones;
    ones = gsMatrix<T>::Ones(midright_coefs.rows(),1);
    midright_coefs.col(1) += tb * ones;
    ones = gsMatrix<T>::Ones(right_coefs.rows(),1);
    right_coefs.col(1) += tb * ones;
    return element;
}

template <class T>
gsMultiPatch<T> makeBottom(const T tw, const T tg, const T tb, const T ts, const T l, const T a, const std::string expr, const gsKnotVector<T> & kv1)
{
    std::vector<gsBSpline<T>> curves = makeCurve(tw,tg,tb,ts,l,a, expr, kv1);
    return makeBottom(tw,tg,tb,ts,l,a,curves);
}

template <class T>
gsMultiPatch<T> makeBottom(const T tw, const T tg, const T tb, const T ts, const T l, const T a, const std::vector<gsBSpline<T>> & curves)
{
    gsKnotVector<T> kv2(0, 1, 0, 2);

    gsBSpline<T> left_crv = curves[0];
    gsBSpline<T> midleft_crv = curves[1];
    gsBSpline<T> mid_crv = curves[2];
    gsBSpline<T> midright_crv = curves[3];
    gsBSpline<T> right_crv = curves[4];

    gsBSplineBasis<T> left_basis, midleft_basis, mid_basis, midright_basis, right_basis;
    left_basis = left_crv.basis();
    midleft_basis = midleft_crv.basis();
    mid_basis = mid_crv.basis();
    midright_basis = midright_crv.basis();
    right_basis = right_crv.basis();

    gsMatrix<T> left_coefs, midleft_coefs, mid_coefs, midright_coefs, right_coefs;
    gsMatrix<T> tmp_coefs;

    gsTensorBSplineBasis<2,real_t> left_tbasis(left_basis.knots(),kv2);
    gsTensorBSplineBasis<2,real_t> midleft_tbasis(midleft_basis.knots(),kv2);
    gsTensorBSplineBasis<2,real_t> mid_tbasis(mid_basis.knots(),kv2);
    gsTensorBSplineBasis<2,real_t> midright_tbasis(midright_basis.knots(),kv2);
    gsTensorBSplineBasis<2,real_t> right_tbasis(right_basis.knots(),kv2);

    gsMatrix<T> coefs_tmp, tmp, ones;
    gsMultiPatch<T> bottom;
    // Get the coefficients
    tmp_coefs.resize(left_crv.coefs().rows(),2);
    tmp_coefs.col(0) = left_basis.anchors().transpose() * l;
    tmp_coefs.col(1) = left_crv.coefs();
    left_coefs = tmp_coefs;

    tmp_coefs.resize(midleft_crv.coefs().rows(),2);
    tmp_coefs.col(0) = midleft_basis.anchors().transpose() * l;
    tmp_coefs.col(1) = midleft_crv.coefs();
    midleft_coefs = tmp_coefs;

    tmp_coefs.resize(mid_crv.coefs().rows(),2);
    tmp_coefs.col(0) = mid_basis.anchors().transpose() * l;
    tmp_coefs.col(1) = mid_crv.coefs();
    mid_coefs = tmp_coefs;

    tmp_coefs.resize(midright_crv.coefs().rows(),2);
    tmp_coefs.col(0) = midright_basis.anchors().transpose() * l;
    tmp_coefs.col(1) = midright_crv.coefs();
    midright_coefs = tmp_coefs;

    tmp_coefs.resize(right_crv.coefs().rows(),2);
    tmp_coefs.col(0) = right_basis.anchors().transpose() * l;
    tmp_coefs.col(1) = right_crv.coefs();
    right_coefs = tmp_coefs;

    // The cosine starts at a, but we want it to start at -tg/2
    ones = gsMatrix<T>::Ones(left_coefs.rows(),1);
    left_coefs.col(1) += (-tg/2.-a) * ones;
    ones = gsMatrix<T>::Ones(midleft_coefs.rows(),1);
    midleft_coefs.col(1) += (-tg/2.-a) * ones;
    ones = gsMatrix<T>::Ones(mid_coefs.rows(),1);
    mid_coefs.col(1) += (-tg/2.-a) * ones;
    ones = gsMatrix<T>::Ones(midright_coefs.rows(),1);
    midright_coefs.col(1) += (-tg/2.-a) * ones;
    ones = gsMatrix<T>::Ones(right_coefs.rows(),1);
    right_coefs.col(1) += (-tg/2.-a) * ones;

    // Now the bottom cosine
    tmp = left_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1).setConstant(-tg/2-2*a-(5*tg)); // the amplitude of the cosine is 2a so we have to cover at least 2a, plus a margin (here 5tg)
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    bottom.addPatch(gsTensorBSpline<2,real_t>(left_tbasis,coefs_tmp));

    tmp = midleft_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1).setConstant(-tg/2-2*a-(5*tg));
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    bottom.addPatch(gsTensorBSpline<2,real_t>(midleft_tbasis,coefs_tmp));

    tmp = mid_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1).setConstant(-tg/2-2*a-(5*tg));
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    bottom.addPatch(gsTensorBSpline<2,real_t>(mid_tbasis,coefs_tmp));

    tmp = midright_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1).setConstant(-tg/2-2*a-(5*tg));
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    bottom.addPatch(gsTensorBSpline<2,real_t>(midright_tbasis,coefs_tmp));

    tmp = right_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1).setConstant(-tg/2-2*a-(5*tg));
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    bottom.addPatch(gsTensorBSpline<2,real_t>(right_tbasis,coefs_tmp));

    for (index_t p=0; p!=5; p++)
        bottom.patch(p).uniformRefine(1,1,1);

    return bottom;
}
#else//gsElasticity_ENABLED

int main (int argc, char** argv)
{
  gsInfo<<"To run this example, compile G+Smo with gsElasticity\n";
  return EXIT_SUCCESS;
}

#endif


