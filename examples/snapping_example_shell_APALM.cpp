/** @file bSplineCurve_example.cpp

    @brief Tutorial on gsBSpline class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <gismo.h>

#ifdef gsKLShell_ENABLED
#include <gsKLShell/src/gsThinShellAssembler.h>
#include <gsKLShell/src/gsMaterialMatrixBase.h>
#include <gsKLShell/src/gsMaterialMatrixLinear.h>
#include <gsKLShell/src/getMaterialMatrix.h>
#endif
#include <gsStructuralAnalysis/src/gsStructuralAnalysisTools/gsStructuralAnalysisUtils.h>

#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMBase.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMLoadControl.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMRiks.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMCrisfield.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsAPALMData.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsAPALM.h>


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

#ifdef gsKLShell_ENABLED
template<class T>
class gsAPALMSnapping : public gsAPALM<T>
{

  using Base = gsAPALM<T>;
  typedef typename Base::solution_t solution_t;

public:
  gsAPALMSnapping(
                    const gsMpiComm & comm,
                    gsALMBase<T> * ALM,
                    const gsAPALMData<T,solution_t> & Data,
                    const gsThinShellAssemblerBase<T> * assembler,
                    std::string dirname,
                    const gsMatrix<T> & refPoints,
                    const gsVector<index_t> & refPatches,
                    const index_t Nx, const index_t Ny, const T l, const T h, const T b)
  :
  Base(ALM,Data,comm),
  m_assembler(assembler),
  m_dirname(dirname),
  m_refPoints(refPoints),
  m_refPatches(refPatches),
  m_Nx(Nx),
  m_Ny(Ny),
  m_l(l),
  m_h(h),
  m_b(b)
  {
    GISMO_ASSERT(refPoints.cols()==1,"Only one point will be used");
    GISMO_ASSERT(refPatches.size()==1,"Only one point will be used");
  }

  void parallelIntervalOutput(const std::vector<std::pair<gsVector<T>,T>> & stepSolutions, const std::vector<T> & stepTimes, index_t level, index_t ID)
  {
    GISMO_ASSERT(stepSolutions.size()==stepTimes.size(),"Solutions and times must have same size, but solutions.size() = "<<stepSolutions.size()<<" and times.size() = "<<stepTimes.size());

    gsStructuralAnalysisOutput<real_t> data(m_dirname + "/interval_"+std::to_string(ID)+".csv",m_refPoints);
    gsMultiPatch<T> deformation,mp_tmp, mp;
    deformation = mp = m_assembler->geometry();
    std::vector<std::string> pointheaders = {"u_x","u_y"};
    std::vector<std::string> otherheaders = {"U-norm","Eps","Sigma","Lambda","time","level"};

    data.init(pointheaders,otherheaders);

    for (size_t k=0; k!=stepSolutions.size(); k++)
    {
        m_assembler->constructDisplacement(stepSolutions[k].first,mp_tmp);
        gsField<T> solField(mp,deformation);

        gsMatrix<> pointResults = mp_tmp.patch(m_refPatches(0)).eval(m_refPoints.col(0));

        T eps = pointResults(1,0) / (m_Ny*m_h);
        T sig = stepSolutions[k].second / (m_Nx*m_b*m_l);

        gsVector<> otherData(6);
        otherData<<stepSolutions[k].first.norm(),eps,sig,stepSolutions[k].second,stepTimes[k],level;
        data.add(pointResults,otherData);
    }
  }

protected:
  const gsThinShellAssemblerBase<T> * m_assembler;

  const std::string m_dirname;

  const gsMatrix<T> m_refPoints;
  const gsVector<index_t> m_refPatches;

  const index_t m_Nx;
  const index_t m_Ny;

  const T m_l;
  const T m_h;
  const T m_b;
};


int main(int argc, char *argv[])
{
    // Input options
    index_t numElevate  = 0;
    index_t numHref     = 2;
    bool    plot       = false;
    bool    write = false;
    bool sequential   = false;

    // Arc length method options
    real_t  dL = 0.1; // Arc length
    real_t  tol = 1e-6;
    real_t  tolU = 1e-6;
    real_t  tolF = 1e-3;
    real_t  relax = 1.0;
    real_t  tau = 1e4;
    index_t SubIntervals  = 2;
    index_t method = 2; // (0: Load control; 1: Riks' method; 2: Crisfield's method; 3: consistent crisfield method; 4: extended iterations)
    bool    quasiNewton = false;
    index_t quasiNewtonInt = -1;
    bool    adaptive = false;
    real_t  perturbation = 0;
    index_t step = 1000;
    index_t maxit = 10;

    index_t verbose = 0;

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

    index_t maxLevel  = 2;

    // real_t Emax = 1.5;

    gsCmdLine cmd("");
    cmd.addInt("i","interior","Number of interior knots",interior);
    cmd.addInt("X","Nx","Number of element in x-direction",Nx);
    cmd.addInt("Y","Ny","Number of element in y-direction",Ny);

    cmd.addInt("r","hRefine",
       "Number of dyadic h-refinement (bisection) steps to perform before solving",
       numHref);
    cmd.addInt("e","degreeElevation",
      "Number of degree elevation steps to perform on the Geometry's basis before solving",
      numElevate);

    cmd.addReal("a","a/l", "Value of a/l", al);
    // cmd.addReal("E","Emax", "Maximum strain", Emax);

    cmd.addInt("m","Method", "Arc length method; 0: Load control; 1: Riks' method; 2: Crisfield's method", method);
    cmd.addReal("L","dL", "arc length", dL);
    cmd.addInt("l","level", "Max level", maxLevel);
    cmd.addReal("A","relaxation", "Relaxation factor for arc length method", relax);

    cmd.addReal("P","perturbation", "perturbation factor", tau);

    cmd.addInt("q","QuasiNewtonInt","Use the Quasi Newton method every INT iterations",quasiNewtonInt);
    cmd.addInt("N", "maxsteps", "Maximum number of steps", step);
    cmd.addInt("n", "SubIntervals", "Number of steps in subintervals", SubIntervals);

    cmd.addInt("v", "verbose", "verbose", verbose);

    cmd.addSwitch("adaptive", "Adaptive length ", adaptive);
    cmd.addSwitch("quasi", "Use the Quasi Newton method", quasiNewton);

    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addSwitch("write", "Write convergence data to file", write);
    cmd.addSwitch("sequential", "Solve sequential (serial -> parallel) instead of mixed", sequential);


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
    index_t pIndex;
    gsMultiPatch<> tmp;
    for (index_t kx = 0; kx!=Nx; kx++)
    {
        tmp = bottom;
        gsNurbsCreator<>::shift2D(tmp,kx*dx,0);
        for (auto patch = tmp.begin(); patch != tmp.end(); patch++)
            pIndex = mp.addPatch(**patch);

        for (index_t ky = 0; ky!=Ny; ky++)
        {
            tmp = element;
            gsNurbsCreator<>::shift2D(tmp,kx*dx,ky*dy);
            for (auto patch = tmp.begin(); patch != tmp.end(); patch++)
                pIndex = mp.addPatch(**patch);
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
            pIndex = mp.addPatch(**patch);
    }
    mp.computeTopology();

    if (plot)
    {
        gsWriteParaview<>(mp,"mp",1000,true,false);
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
            bc.addCondition(bit->patch,bit->side(),condition_type::collapsed,0,0,false,1);
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

    // point load
    gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();
    // point load
    gsVector<> point = writePoints.col(0);
    gsVector<> load (2);
    load << 0.0, 1.0;
    pLoads.addLoad(point, load, writePatches(0) );

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
    parameters[0] = &E;
    parameters[1] = &nu;
    options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",1);
    options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
    options.addSwitch("Compressibility","Compressibility",true);
    materialMatrix = getMaterialMatrix<2,real_t>(mp,t,parameters,options);

    gsThinShellAssembler<2, real_t, false> assembler(mp,dbasis,bc,force,materialMatrix);
    // Set the penalty parameter for the interface C1 continuity
    gsOptionList assemblerOptions = assembler.options();
    assemblerOptions.setInt("Continuity",0);
    assembler.setOptions(assemblerOptions);
    assembler.setPointLoads(pLoads);

    // Assemble linear system to obtain the force vector
    assembler.assemble();
    gsVector<> Force = assembler.rhs();
    gsDebugVar(Force.norm());
    gsInfo << "Initialized system with " << assembler.numDofs() << " dofs.\n";

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
    gsStructuralAnalysisOps<real_t>::ALResidual_t ALResidual = [&assembler,&mp_def,&Force](gsVector<real_t> const &x, real_t lam, gsVector<real_t> & result)
    {
        ThinShellAssemblerStatus status;
        assembler.constructSolution(x,mp_def);
        status = assembler.assembleVector(mp_def);
        result = Force - lam * Force - assembler.rhs(); // assembler rhs - force = Finternal
        return status == ThinShellAssemblerStatus::Success;
    };
    //! [Define jacobian and residual]

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
      arcLength->options().setReal("Length",dL);
      if (method==2)
      {
          arcLength->options().setInt("AngleMethod",0); // 0: step, 1: iteration
          arcLength->options().setReal("Scaling",0.0);
      }
      arcLength->options().setSwitch("AdaptiveLength",adaptive);
      arcLength->options().setInt("AdaptiveIterations",5);
      arcLength->options().setReal("Perturbation",tau);
      arcLength->options().setReal("Tol",tol);
      arcLength->options().setReal("TolU",tolU);
      arcLength->options().setReal("TolF",tolF);
      arcLength->options().setInt("MaxIter",maxit);
      arcLength->options().setSwitch("Verbose",verbose);
      arcLength->options().setReal("Relaxation",relax);
      if (quasiNewtonInt>0)
      {
        quasiNewton = true;
        arcLength->options().setInt("QuasiIterations",quasiNewtonInt);
    }
    arcLength->options().setSwitch("Quasi",quasiNewton);


    gsDebug<<arcLength->options();
    arcLength->applyOptions();
    arcLength->initialize();

    deformation = mp;

    // Make objects for previous solutions
    real_t Lold = 0;
    gsMatrix<> Uold = Force;
    Uold.setZero();

    typedef std::pair<gsVector<real_t>,real_t> solution_t;

    /*
      \a solutions is a container the for each level contains the solutions per point
      \a points contains all the points across levels in the format (level, U, lambda)
      \a refPoints is a container that contains (level, U, lambda) of the points from which a refinement should START in level+1
      \a errors is a container that contains the error[l][i] e_i at the ith point of level l
    */

    std::vector<solution_t> solutions;
    std::vector<real_t> times;
    std::vector<index_t> levels;

    index_t level = 0;
    gsInfo<<"------------------------------------------------------------------------------------\n";
    gsInfo<<"\t\t\tLevel "<<level<<" (dL = "<<dL<<") -- Coarse grid \n";
    gsInfo<<"------------------------------------------------------------------------------------\n";

    gsAPALMData<real_t,solution_t> apalmData;
    apalmData.options().setInt("MaxLevel",maxLevel);
    apalmData.options().setInt("Verbose",verbose);
    apalmData.options().setReal("Tolerance",1e-3);

    const gsMpi & mpi = gsMpi::init(argc, argv);
    gsMpiComm comm = mpi.worldComm();

    std::string dirname = "ArcLengthResults/snapping_2D_al=" + std::to_string(al) + "_r=" + std::to_string(numHref) + "_e=" + std::to_string(numElevate);
    std::string cores = "_ncores="+std::to_string(comm.size());
    dirname = dirname + cores;

    // Prepare and create directory with dirname
    dirname = gsFileManager::getCurrentPath() + dirname;
    GISMO_ENSURE(gsFileManager::mkdir(dirname),"Failed to create directory " + dirname);
    // Made directory
    gsParaviewCollection collection(dirname + "/" + output);

    gsAPALMSnapping<real_t> apalm(comm,arcLength,apalmData,&assembler,dirname,writePoints,writePatches,Nx,Ny,l,h,b);
    apalm.options().setSwitch("Verbose",(verbose>0));
    apalm.options().setInt("SubIntervals",SubIntervals);
    apalm.options().setInt("MaxIt",10000);
    apalm.initialize();

    real_t serialTime = 0, parallelTime = 0;
    if (!sequential)
    {
        parallelTime = mpi.wallTime();
        apalm.solve(step+1);
        parallelTime = mpi.wallTime() - parallelTime;
        if (apalm.isMain()) gsInfo<<"Parallel time = "<<parallelTime<<"\n";
    }
    else
    {
        serialTime = mpi.wallTime();
        apalm.serialSolve(step+1);
        serialTime = mpi.wallTime() - serialTime;
        if (apalm.isMain()) gsInfo<<"Serial time = "<<serialTime<<"\n";

        if (maxLevel!=0)
        {
            parallelTime = mpi.wallTime();
            apalm.parallelSolve();
            parallelTime = mpi.wallTime() - parallelTime;
            if (apalm.isMain()) gsInfo<<"Parallel time = "<<parallelTime<<"\n";
        }
    }

    if (apalm.isMain())
    {
        solutions = apalm.getFlatSolutions();
        times     = apalm.getFlatTimes();
        levels    = apalm.getFlatLevels();

        if (plot || write)
        {
            gsField<> solField;
            gsMultiPatch<> mp_tmp, displacement;
            std::vector<std::string> pointheaders = {"u_x","u_y"};
            std::vector<std::string> otherheaders = {"U-norm","Eps","Sigma","Lambda","time","level"};

            gsStructuralAnalysisOutput<real_t> data(dirname + "/data.csv",writePoints);
            if (write)
                data.init(pointheaders,otherheaders);

            for (size_t k=0; k!= solutions.size(); k++)
            {
                Lold = solutions[k].second;
                Uold = solutions[k].first;

                assembler.constructDisplacement(Uold,mp_tmp);
                solField = gsField<>(mp,mp_tmp);

                if (write)
                    if (writePoints.cols()!=0)
                    {
                        gsMatrix<> pointResults = mp_tmp.patch(writePatches(0)).eval(writePoints.col(0));

                        real_t eps = pointResults(1,0) / (Ny*h);
                        real_t sig = Lold / (Nx*b*l);

                        gsVector<> otherData(6);
                        otherData<<Uold.norm(),eps,sig,Lold,times[k],levels[k];
                        data.add(pointResults,otherData);
                    }

                if (plot)
                {
                    std::string fileName = dirname + "/" + output + util::to_string(k) + "_";
                    // creating a container to plot all fields to one Paraview file
                    gsWriteParaview<>(solField, fileName, 200,true);
                    fileName = gsFileManager::getFilename(fileName);
                    for (size_t p=0; p!=mp.nPatches(); p++)
                    {
                        collection.addPart(fileName + std::to_string(p) + ".vts",k,"",p);
                        collection.addPart(fileName + std::to_string(p) + "_mesh.vtp",k,"",p);
                    }
                }
            }
            if (plot)
            {
                collection.save();
            }
        }
    }
    if (apalm.isMain())
    {
        std::ofstream file;
        file.open(dirname + "/times.txt");
        if (sequential)
            file<<"serial   time: "<<serialTime<<" s\n";
        file<<"parallel time: "<<parallelTime<<" s\n";
        file.close();
    }

    delete arcLength;
    return EXIT_SUCCESS;
}
#else//gsKLShell_ENABLED
int main(int argc, char *argv[])
{
    gsWarn<<"G+Smo is not compiled with the gsKLShell module.";
    return EXIT_FAILURE;
}
#endif

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


