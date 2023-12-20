/** @file bSplineCurve_example.cpp

    @brief Tutorial on gsBSpline class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <iostream>

#include <gsKLShell/src/gsThinShellAssembler.h>
#include <gsKLShell/src/gsMaterialMatrixBase.h>
#include <gsKLShell/src/gsMaterialMatrixLinear.h>
#include <gsKLShell/src/getMaterialMatrix.h>
#include <gsKLShell/src/gsMaterialMatrixIntegrate.h>
#include <gsKLShell/src/gsThinShellUtils.h>

#include <gsStructuralAnalysis/src/gsStaticSolvers/gsStaticNewton.h>
#include <gsStructuralAnalysis/src/gsStaticSolvers/gsStaticComposite.h>
#include <gsStructuralAnalysis/src/gsStaticSolvers/gsControlDisplacement.h>

#include <gismo.h>

using namespace gismo;

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
    int numHref     = 2;
    bool plot       = false;
    bool write = false;

    // Arc length method options
    real_t dL = 0; // General arc length
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
    real_t perturbation = 0;
    int step = 1000;
    index_t maxit = 50;

    std::string wn("data.csv");

    std::string output = "solutionElasticity";

    // geometry options
    real_t start = 0; // starting knot
    real_t end = 1; // ending knot
    index_t interior = 4; // number of interior knots
    index_t multEnd = 3; // multiplicity at the two end knots

    index_t Nx = 3;
    index_t Ny = 2;

    real_t tw = 1.5e-3;
    real_t tg = 1.0e-3;
    real_t tb = 1.5e-3;
    real_t ts = 1.0e-3;
    real_t l  = 10.0e-3;
    real_t al = 0.3;
    real_t b  = 3e-3;
    real_t h  = 2*tg+tb+ts;

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
    cmd.addReal("L","dLb", "arc length", dL);
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
    auto elementlabels = element.getBoxProperty<std::string>("label");

    gsWriteParaview(element,"element",1000,true);

    gsInfo<<"-----------------------------------------------------------------------------------------------------\n";
    gsInfo<<"-----------------------------------Making the bottom block-------------------------------------------\n";
    gsInfo<<"-----------------------------------------------------------------------------------------------------\n";
    gsMultiPatch<> bottom = makeBottom(tw,tg,tb,ts,l,a,curves);
    auto bottomlabels = bottom.getBoxProperty<std::string>("label");

    gsWriteParaview(bottom,"bottom",1000,true);

    gsInfo<<"-----------------------------------------------------------------------------------------------------\n";
    gsInfo<<"-----------------------------------Making the top block-------------------------------------------\n";
    gsInfo<<"-----------------------------------------------------------------------------------------------------\n";
    gsMultiPatch<> top = makeTop(tw,tg,tb,ts,l,a,curves);
    auto toplabels = top.getBoxProperty<std::string>("label");

    gsWriteParaview(top,"top",1000,true);

    gsInfo<<"-----------------------------------------------------------------------------------------------------\n";
    gsInfo<<"-----------------------------------Making the geometry--------------------------------------------\n";
    gsInfo<<"-----------------------------------------------------------------------------------------------------\n";
    gsMultiPatch<> mp;
    auto labels = mp.addBoxProperty("label",std::string());

    real_t dx = l;
    real_t dy = 2*tg+ts+tb;
    index_t topmid_ID;
    index_t pIndex;
    gsMultiPatch<> tmp;
    for (index_t kx = 0; kx!=Nx; kx++)
    {
        tmp = bottom;
        gsNurbsCreator<>::shift2D(tmp,kx*dx,0);
        for (auto patch = tmp.begin(); patch != tmp.end(); patch++)
        {
            pIndex = mp.addPatch(**patch);
            labels[pIndex] = bottomlabels[(*patch)->id()];
        }

        for (index_t ky = 0; ky!=Ny; ky++)
        {
            tmp = element;
            gsNurbsCreator<>::shift2D(tmp,kx*dx,ky*dy);
            for (auto patch = tmp.begin(); patch != tmp.end(); patch++)
            {
                pIndex = mp.addPatch(**patch);
                labels[pIndex] = elementlabels[(*patch)->id()];
            }
        }
        tmp = top;
        gsNurbsCreator<>::shift2D(tmp,kx*dx,Ny*dy);
        if (kx == math::floor(Nx/2.))
        {
            if (Nx % 2 ==0)
                topmid_ID = mp.nPatches()+8;
            else
                topmid_ID = mp.nPatches()+10;

            gsDebugVar(topmid_ID);
        }
        for (auto patch = tmp.begin(); patch != tmp.end(); patch++)
        {
            pIndex = mp.addPatch(**patch);
            labels[pIndex] = toplabels[(*patch)->id()];
        }
    }
    mp.computeTopology();

    // For different materials for snapping or bearing patches.
    index_t i=0;
    std::vector<size_t> snapIdx, bearIdx;
    for (auto label = labels.vector().begin(); label!=labels.vector().end(); label++, i++)
    {
        if      (*label=="snapping") snapIdx.push_back(i);
        else if (*label=="bearing")  bearIdx.push_back(i);
    }

    if (plot)
    {
        gsWriteParaview<>(mp,"mp",100,true,false);
    }

    gsMultiPatch<> mp_def = mp;
    gsMultiPatch<> deformation;

    // p-refine
    if (numElevate!=0)
        mp.degreeElevate(numElevate);

    // h-refine
    for (int r =0; r < numHref; ++r)
        mp.uniformRefine();

    gsInfo<<"Multi-patch info:\n";
    gsInfo<<"* Number of patches: "<<mp.nPatches()<<"\n";
    gsInfo<<"* First basis      : "<<mp.basis(0)<<"\n";

    // define basis
    gsMultiBasis<> dbasis(mp);

    // Boundary conditions
    gsBoundaryConditions<> bc;

    gsMatrix<> bbox;
    mp.boundingBox(bbox);
    real_t band = 0.01 * (bbox(1,1) - bbox(1,0));
    real_t ymin = bbox(1,0);
    real_t ymax = bbox(1,1);
    if (dL==0) dL = Ny*h / step * 2;
    gsConstantFunction<> displ(dL,2);
    std::vector<patchSide> topBdr;
    for (typename gsBoxTopology::const_biterator bit = mp.bBegin(); bit!=mp.bEnd(); bit++)
    {
        // get the parametric center of the boundary
        gsMatrix<> pRange = mp.patch(bit->patch).boundary(bit->side())->support();
        real_t center = (pRange(0,1)-pRange(0,0)) / 2. + pRange(0,0);
        gsMatrix<> centerpt(1,1);
        centerpt<<center;
        gsMatrix<> eval = mp.patch(bit->patch).boundary(bit->side())->eval(centerpt);

        if ( (eval(1,0) <= ymin + band) && (eval(1,0) >= ymin - band))
        {
            bc.addCondition(bit->patch,bit->side(),condition_type::dirichlet,0,0,false,0);
            bc.addCondition(bit->patch,bit->side(),condition_type::dirichlet,0,0,false,1);
        }
        else if ((eval(1,0) <= ymax + band) && (eval(1,0) >= ymax - band))
        {
            bc.addCondition(bit->patch,bit->side(),condition_type::dirichlet,0,0,false,0);
            bc.addCondition(bit->patch,bit->side(),condition_type::dirichlet,&displ,0,false,1);
            topBdr.push_back(*bit);
        }
    }
    bc.addCornerValue(boundary::corner::southwest,0.0,0,0,-1);

    // Ref points
    gsMatrix<> writePoints(2,1);
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
    writePoints.col(0)<<xBoundary,supp(1,1);

    gsMatrix<> evalPoints;
    mp.patch(topmid_ID).eval_into(writePoints,evalPoints);

    // fd.getId(20, bc); // id=2: boundary conditions
    bc.setGeoMap(mp);
    gsInfo<<"Boundary conditions:\n"<< bc <<"\n";

    // Loads
    gsFunctionExpr<> bodyForce("0","0",2);
    gsInfo<<"Body force function "<< bodyForce << "\n";

    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
    gsInfo<<mp<<"\n";

    // -----------------------------------------------------------------------------------------------
    // --------------------------------------Solve Static problem-------------------------------------
    // -----------------------------------------------------------------------------------------------

    // Loads
    gsFunctionExpr<> force("0","0",2);

    gsMaterialMatrixBase<real_t>* materialMatrix;
    gsOptionList options;
    std::vector<gsFunctionSet<>*> parameters(2);
    gsFunctionExpr<> t(std::to_string(b),2);
    gsFunctionExpr<> E(std::to_string(78e6),2);
    gsConstantFunction<> nu(0.4,2);
    gsFunctionExpr<> rho(std::to_string(1.0),2);
    parameters[0] = &E;
    parameters[1] = &nu;
    options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",1);
    options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
    options.addSwitch("Compressibility","Compressibility",true);
    materialMatrix = getMaterialMatrix<2,real_t>(mp,t,parameters,rho,options);

    gsThinShellAssembler<2, real_t, false> assembler(mp,dbasis,bc,force,materialMatrix);
    // Set the penalty parameter for the interface C1 continuity
    gsOptionList assemblerOptions = assembler.options();
    assemblerOptions.setInt("Continuity",0);
    assembler.setOptions(assemblerOptions);

    //! [Define jacobian and residual]
    // Function for the Jacobian
    gsStructuralAnalysisOps<real_t>::Jacobian_t Jacobian = [&assembler,&mp_def](gsVector<real_t> const &x, gsSparseMatrix<real_t> & m)
    {
        ThinShellAssemblerStatus status;
        assembler.constructSolution(x,mp_def);
        status = assembler.assembleMatrix(mp_def);
        m = assembler.matrix();
        return status == ThinShellAssemblerStatus::Success;
    };
    // Function for the Residual
    gsStructuralAnalysisOps<real_t>::ALResidual_t ALResidual = [&bc,&displ,&assembler,&mp_def](gsVector<real_t> const &x, real_t lam, gsVector<real_t> & result)
    {
        ThinShellAssemblerStatus status;
        displ.setValue(lam,2);
        assembler.updateBCs(bc);
        assembler.constructSolution(x,mp_def);
        status = assembler.assembleVector(mp_def);
        result = assembler.rhs(); 
        return status == ThinShellAssemblerStatus::Success;
    };
    //! [Define jacobian and residual]

    // Assemble linear system to obtain the force vector
    assembler.assemble();
    gsSparseMatrix<> matrix = assembler.matrix();
    gsVector<> vector = assembler.rhs();
    gsInfo << "Initialized system with " << assembler.numDofs() << " dofs.\n";

    //=============================================//
                  // Solving //
    //=============================================//

    gsStaticNewton<real_t> staticNR(matrix,vector,Jacobian,ALResidual);
    gsOptionList NROptions = staticNR.options();
    NROptions.setInt("verbose",1);
    NROptions.setInt("maxIt",maxit);
    NROptions.setReal("tolU",tolU);
    NROptions.setReal("tolF",tolF);
    staticNR.setOptions(NROptions);

    gsStaticComposite<real_t> staticSolver({&staticNR});
    staticSolver.initialize();

    gsControlDisplacement<real_t> control(&staticSolver);

    std::string dirname = "ArcLengthResults/snapping_DC_2D_" + std::to_string(Nx) + "x" + std::to_string(Ny+1) + "_al=" + std::to_string(al) + "_r=" + std::to_string(numHref) + "_e=" + std::to_string(numElevate);

    // Prepare and create directory with dirname
    dirname = gsFileManager::getCurrentPath() + dirname;
    GISMO_ENSURE(gsFileManager::mkdir(dirname),"Failed to create directory " + dirname);
    // Made directory

    gsParaviewCollection collection(dirname + "/" + output);
    deformation = mp;

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
        file  << "Sigma" << ",";
        file  << "Lambda" << ","
        << "Indicator"
        << "\n";
        file.close();
    }

    real_t dL0 = dL;
    int reset = 0;
    gsMultiPatch<> mp_def0 = mp_def;
    real_t indicator;

    index_t  k = 0;
    real_t D = 0, eps = 0, sig = 0;
    displ.setValue(D - dL,2);
    assembler.updateBCs(bc);
    gsMatrix<> solVector;
    real_t time = 0;
    while (eps<=Emax && k < step)
    {
        gsInfo<<"Load step "<<k<<"; D = "<<D<<"; dL = "<<dL<<"eps = "<<eps<<"\n";

	    gsStopwatch timer;
        gsStatus status = control.step(dL);
	    time += timer.stop();
        if (status==gsStatus::NotConverged || status==gsStatus::AssemblyError)
        {
            dL = dL/2;
            displ.setValue(D+dL,2);
            reset = 1;
            mp_def = mp_def0;
            gsInfo<<"Iterations did not converge\n";
            continue;
        }

        solVector = control.solutionU();
        indicator = staticNR.indicator();

        // Recover force on the north boundaries
        assembler.constructSolution(solVector,mp_def);
        gsVector<> F = assembler.boundaryForce(mp_def,topBdr);
        real_t Load = F[1];

        assembler.constructDisplacement(solVector,mp_def);
        if (plot)
        {
            std::string fileName = dirname + "/" + output + util::to_string(k) + "_";
            // creating a container to plot all fields to one Paraview file
            gsField<> solField(mp,mp_def);
            gsWriteParaview<>(solField, fileName, 50,true);
            fileName = gsFileManager::getFilename(fileName);
            for (size_t p=0; p!=mp.nPatches(); p++)
            {
                collection.addPart(fileName + std::to_string(p) + ".vts",k,"",p);
                collection.addPart(fileName + std::to_string(p) + "_mesh.vtp",k,"",p);
            }
        }

        gsMatrix<> out_def = mp_def.patch(writePatches(0)).eval(writePoints.col(0));
        gsDebugVar(out_def);
        eps = out_def(1,0) / (Ny*h);
        sig = Load / (Nx*b*l);
        gsDebugVar(eps);
        if (write)
        {

            std::ofstream file;
            file.open(dirname + "/" + wn,std::ofstream::out | std::ofstream::app);
            file    << std::setprecision(6)
            << control.solutionU().norm() << ",";
            file<< out_def(0,0) << ","
            << out_def(1,0) << ",";
            file<< eps<<",";
            file<< sig<<",";
            file    << control.solutionL() << "," << indicator<<"\n";
            file.close();
        }


        // if (reset!=1)
        // {
          dL = dL0;
        // }
        reset = 0;

        mp_def0 = mp_def;
        D += dL;
        k++;

        gsInfo<<"--------------------------------------------------------------------------------------------------------------\n";
    }

    if (plot)
    {
        collection.save();
    }

    if (write)
    {
        std::ofstream file;
        file.open(dirname + "/times.txt");
        file<<"solution time: "<<time<<" s\n";
        file.close();
    }

    delete materialMatrix;
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
    auto labels = top.addBoxProperty("label",std::string());
    index_t pIndex;
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
    pIndex = top.addPatch(gsTensorBSpline<2,real_t>(mid_tbasis,coefs_tmp));
    labels[pIndex] = "bearing";

    // Now the bottom cosine
    tmp = left_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1) += ts * ones;
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    pIndex = top.addPatch(gsTensorBSpline<2,real_t>(left_tbasis,coefs_tmp));
    labels[pIndex] = "snapping";

    tmp = midleft_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1) += ts * ones;
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    pIndex = top.addPatch(gsTensorBSpline<2,real_t>(midleft_tbasis,coefs_tmp));
    labels[pIndex] = "snapping";

    tmp = mid_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1) += ts * ones;
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    pIndex = top.addPatch(gsTensorBSpline<2,real_t>(mid_tbasis,coefs_tmp));
    labels[pIndex] = "snapping";

    tmp = midright_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1) += ts * ones;
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    pIndex = top.addPatch(gsTensorBSpline<2,real_t>(midright_tbasis,coefs_tmp));
    labels[pIndex] = "snapping";

    tmp = right_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1) += ts * ones;
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    pIndex = top.addPatch(gsTensorBSpline<2,real_t>(right_tbasis,coefs_tmp));
    labels[pIndex] = "snapping";

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
    pIndex = top.addPatch(gsTensorBSpline<2,real_t>(left_tbasis,coefs_tmp));
    labels[pIndex] = "bearing";

    tmp = right_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1) += tg*ones;
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    pIndex = top.addPatch(gsTensorBSpline<2,real_t>(right_tbasis,coefs_tmp));
    labels[pIndex] = "bearing";

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
    tmp.col(1).setConstant(tg/2+ts+tg+tb+(5*tg)); //last number is margin
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    pIndex = top.addPatch(gsTensorBSpline<2,real_t>(left_tbasis,coefs_tmp));
    labels[pIndex] = "bearing";

    tmp = midleft_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1).setConstant(tg/2+ts+tg+tb+(5*tg));
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    pIndex = top.addPatch(gsTensorBSpline<2,real_t>(midleft_tbasis,coefs_tmp));
    labels[pIndex] = "bearing";

    tmp = mid_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1).setConstant(tg/2+ts+tg+tb+(5*tg));
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    pIndex = top.addPatch(gsTensorBSpline<2,real_t>(mid_tbasis,coefs_tmp));
    labels[pIndex] = "bearing";

    tmp = midright_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1).setConstant(tg/2+ts+tg+tb+(5*tg));
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    pIndex = top.addPatch(gsTensorBSpline<2,real_t>(midright_tbasis,coefs_tmp));
    labels[pIndex] = "bearing";

    tmp = right_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1).setConstant(tg/2+ts+tg+tb+(5*tg));
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    pIndex = top.addPatch(gsTensorBSpline<2,real_t>(right_tbasis,coefs_tmp));
    labels[pIndex] = "bearing";

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
    auto labels = element.addBoxProperty("label",std::string());
    index_t pIndex;
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
    pIndex = element.addPatch(gsTensorBSpline<2,real_t>(mid_tbasis,coefs_tmp));
    labels[pIndex] = "bearing";

    // Now the bottom cosine
    tmp = left_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1) += ts * ones;
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    pIndex = element.addPatch(gsTensorBSpline<2,real_t>(left_tbasis,coefs_tmp));
    labels[pIndex] = "snapping";

    tmp = midleft_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1) += ts * ones;
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    pIndex = element.addPatch(gsTensorBSpline<2,real_t>(midleft_tbasis,coefs_tmp));
    labels[pIndex] = "snapping";

    tmp = mid_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1) += ts * ones;
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    pIndex = element.addPatch(gsTensorBSpline<2,real_t>(mid_tbasis,coefs_tmp));
    labels[pIndex] = "snapping";

    tmp = midright_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1) += ts * ones;
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    pIndex = element.addPatch(gsTensorBSpline<2,real_t>(midright_tbasis,coefs_tmp));
    labels[pIndex] = "snapping";

    tmp = right_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1) += ts * ones;
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    pIndex = element.addPatch(gsTensorBSpline<2,real_t>(right_tbasis,coefs_tmp));
    labels[pIndex] = "snapping";

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
    pIndex = element.addPatch(gsTensorBSpline<2,real_t>(left_tbasis,coefs_tmp));
    labels[pIndex] = "bearing";

    tmp = right_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1) += tg*ones;
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    pIndex = element.addPatch(gsTensorBSpline<2,real_t>(right_tbasis,coefs_tmp));
    labels[pIndex] = "bearing";

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
    pIndex = element.addPatch(gsTensorBSpline<2,real_t>(left_tbasis,coefs_tmp));
    labels[pIndex] = "bearing";

    tmp = midleft_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1) += tb * ones;
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    pIndex = element.addPatch(gsTensorBSpline<2,real_t>(midleft_tbasis,coefs_tmp));
    labels[pIndex] = "bearing";

    tmp = mid_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1) += tb * ones;
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    pIndex = element.addPatch(gsTensorBSpline<2,real_t>(mid_tbasis,coefs_tmp));
    labels[pIndex] = "bearing";

    tmp = midright_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1) += tb * ones;
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    pIndex = element.addPatch(gsTensorBSpline<2,real_t>(midright_tbasis,coefs_tmp));
    labels[pIndex] = "bearing";

    tmp = right_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1) += tb * ones;
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    pIndex = element.addPatch(gsTensorBSpline<2,real_t>(right_tbasis,coefs_tmp));
    labels[pIndex] = "bearing";

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
    auto labels = bottom.addBoxProperty("label",std::string());
    index_t pIndex;
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
    tmp.col(1).setConstant(-tg/2-2*a-(5*tg)); // the amplitude of the cosine is 2a so we have to cover at least 2a, plus a margin (here 5*tg)
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    pIndex = bottom.addPatch(gsTensorBSpline<2,real_t>(left_tbasis,coefs_tmp));
    labels[pIndex] = "bearing";

    tmp = midleft_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1).setConstant(-tg/2-2*a-(5*tg));
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    pIndex = bottom.addPatch(gsTensorBSpline<2,real_t>(midleft_tbasis,coefs_tmp));
    labels[pIndex] = "bearing";

    tmp = mid_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1).setConstant(-tg/2-2*a-(5*tg));
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    pIndex = bottom.addPatch(gsTensorBSpline<2,real_t>(mid_tbasis,coefs_tmp));
    labels[pIndex] = "bearing";

    tmp = midright_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1).setConstant(-tg/2-2*a-(5*tg));
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    pIndex = bottom.addPatch(gsTensorBSpline<2,real_t>(midright_tbasis,coefs_tmp));
    labels[pIndex] = "bearing";

    tmp = right_coefs;
    ones = gsMatrix<T>::Ones(tmp.rows(),1);
    coefs_tmp.resize(2*tmp.rows(),tmp.cols());
    coefs_tmp.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    tmp.col(1).setConstant(-tg/2-2*a-(5*tg));
    coefs_tmp.block(tmp.rows(),0,tmp.rows(),tmp.cols()) = tmp;
    pIndex = bottom.addPatch(gsTensorBSpline<2,real_t>(right_tbasis,coefs_tmp));
    labels[pIndex] = "bearing";

    for (index_t p=0; p!=5; p++)
        bottom.patch(p).uniformRefine(1,1,1);

    return bottom;
}


