/** @file benchmark_Wrinkling.cpp

    @brief Computes the wrinkling behaviour of a thin sheet

    Fig 12 of:

    Verhelst, H. M., Möller, M., Den Besten, J. H., Mantzaflaris, A., & Kaminski, M. L. (2021).
    Stretch-Based Hyperelastic Material Formulations for Isogeometric Kirchhoff–Love Shells with Application to Wrinkling.
    Computer-Aided Design, 139, 103075. https://doi.org/10.1016/j.cad.2021.103075

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#include <gismo.h>

#include <gsKLShell/gsThinShellAssembler.h>
#include <gsKLShell/getMaterialMatrix.h>
#include <gsKLShell/gsMaterialMatrixTFT.h>
#include <gsKLShell/gsFunctionSum.h>

#include <gsStructuralAnalysis/gsALMBase.h>
#include <gsStructuralAnalysis/gsALMLoadControl.h>
#include <gsStructuralAnalysis/gsALMRiks.h>
#include <gsStructuralAnalysis/gsALMCrisfield.h>

#include <gsUnstructuredSplines/src/gsApproxC1Spline.h>
#include <gsUnstructuredSplines/src/gsAlmostC1.h>
#include <gsUnstructuredSplines/src/gsSmoothInterfaces.h>

using namespace gismo;

template <class T>
void initStepOutput( const std::string name, const gsMatrix<T> & points);

template <class T>
void writeStepOutput(const gsALMBase<T> * arcLength, const gsMultiPatch<T> & deformation, const std::string name, const gsMatrix<T> & points, const gsMatrix<index_t> & patches);

int main (int argc, char** argv)
{
    // Input options
    int numElevate  = 2;
    int numHref     = 5;
    bool plot       = false;
    bool smooth     = false;
    bool mesh = false;
    bool stress       = false;
    bool quasiNewton = false;
    int quasiNewtonInt = -1;
    bool adaptive = false;
    int step = 10;
    int method = 2; // (0: Load control; 1: Riks' method; 2: Crisfield's method; 3: consistent crisfield method; 4: extended iterations)
    bool symmetry = false;
    bool deformed = false;

    index_t Compressibility = 0;
    index_t material = 3;
    index_t impl = 1; // 1= analytical, 2= generalized, 3= spectral
    bool TFT = false;

    real_t relax = 1.0;

    int result = 0;

    bool write = false;
    bool writeG = false;

    index_t maxit = 50;

    // Arc length method options
    real_t dL = 1; // General arc length
    real_t tol = 1e-3;
    real_t tolU = 1e-3;
    real_t tolF = 1e-3;

    index_t testCase = 0;

    std::string wn("data.csv");

    std::string assemberOptionsFile("options/solver_options.xml");

    gsCmdLine cmd("Wrinkling analysis with thin shells.");
    cmd.addString( "f", "file", "Input XML file for assembler options", assemberOptionsFile );

    cmd.addInt("r","hRefine", "Number of dyadic h-refinement (bisection) steps to perform before solving", numHref);
    cmd.addInt("e","degreeElevation", "Number of degree elevation steps to perform on the Geometry's basis before solving", numElevate);

    cmd.addInt( "M", "Material", "Material law",  material );
    cmd.addInt( "c", "Compressibility", "1: compressible, 0: incompressible",  Compressibility );
    cmd.addInt( "I", "Implementation", "Implementation: 1= analytical, 2= generalized, 3= spectral",  impl );
    cmd.addSwitch( "TFT", "Use Tension-Field Theory",  TFT );

    cmd.addInt("m","Method", "Arc length method; 1: Crisfield's method; 2: RIks' method.", method);
    cmd.addReal("L","dL", "arc length", dL);
    cmd.addReal("A","relaxation", "Relaxation factor for arc length method", relax);

    cmd.addInt("q","QuasiNewtonInt","Use the Quasi Newton method every INT iterations",quasiNewtonInt);
    cmd.addInt("N", "maxsteps", "Maximum number of steps", step);

    cmd.addSwitch("smooth", "Use a smooth basis (maximum regularity)", smooth);

    cmd.addInt("t", "testcase", "testcase", testCase);
    cmd.addSwitch("adaptive", "Adaptive length ", adaptive);
    cmd.addSwitch("quasi", "Use the Quasi Newton method", quasiNewton);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addSwitch("mesh", "Plot mesh?", mesh);
    cmd.addSwitch("stress", "Plot stress in ParaView format", stress);
    cmd.addSwitch("write", "Write output to file", write);
    cmd.addSwitch("writeG", "Write refined geometry", writeG);
    cmd.addSwitch("symmetry", "Use symmetry boundary condition (different per problem)", symmetry);
    cmd.addSwitch("deformed", "plot on deformed shape", deformed);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsFileData<> fd(assemberOptionsFile);
    gsOptionList opts;
    fd.getFirst<gsOptionList>(opts);


    real_t aDim,bDim;
    real_t thickness, E_modulus, PoissonRatio, Density, Ratio, MU1, ALPHA1;
    real_t mu, C01,C10;
    if (testCase==0 || testCase==1 || testCase == 2)
    {
      thickness = 0.14e-3;
      E_modulus     = 1;
      PoissonRatio = 0;
      Density = 1e0;
      Ratio = 7.0;
      if ((!Compressibility) && (material!=0))
        PoissonRatio = 0.5;
      else
        PoissonRatio = 0.499;

      if (material==3)
      {
        C10 = 6.21485502e4; // c1/2
        C01 = 15.8114570e4; // c2/2
        Ratio = C10/C01;
        mu = 2*(C01+C10);
      }
      else
      {
        C10 = 19.1010178e4;
        mu = 2*C10;
      }
      E_modulus = 2*mu*(1+PoissonRatio);
      gsDebug<<"E = "<<E_modulus<<"; nu = "<<PoissonRatio<<"; mu = "<<mu<<"; ratio = "<<Ratio<<"\n";
    }
    else if (testCase==3)
    {
      thickness = 1.0e-3;
      E_modulus     = 1;
      PoissonRatio = 0.4999;
      Density = 1e0;
      Ratio = 0;
      MU1 = 749.18;
      ALPHA1 = 17.14;
      // material = 4;
      // impl = 3;
    }
    else
      GISMO_ERROR("Test case" << testCase<<" unknown.");

    /////////////////////////////////////////////////////////////////////////////////
    //////////////////////////Geometric data/////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////
    // Geometry
    gsMultiPatch<> mp,mp_def;
    //Probes
    gsMatrix<> writePoints;
    gsMatrix<index_t> writePatches;    
    // BCs
    gsBoundaryConditions<> BCs;
    gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();
    // Export
    std::string output = "solution";
    std::string dirname = "ArcLengthResults";

    if (testCase==0)
    {
      // Geometry
      bDim = 0.14; aDim = 2*bDim;
      mp.addPatch(gsNurbsCreator<>::BSplineRectangle(0,0,aDim/2., bDim/2.));
      gsInfo<<"alpha = "<<aDim/bDim<<"; beta = "<<bDim/thickness<<"\n";

      // Probes
      writePoints.resize(2,1);
      writePoints.col(0)<< 1.0,1.0;

      writePatches.resize(1,1);
      writePatches.row(0)<<0;

      // BCs
      BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0 ,false,0);

      BCs.addCondition(boundary::east, condition_type::collapsed, 0, 0 ,false,0);
      BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0 ,false,1);

      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 2 - z.

      real_t Load = 1e0;
      gsVector<> point(2); point<< 1.0, 0.5 ;
      gsVector<> load (2); load << Load,0.0;
      pLoads.addLoad(point, load, 0 );
      
      // Export
      dirname = dirname + "/QuarterSheet_-r" + std::to_string(numHref) + "-e" + std::to_string(numElevate) + "-M" + std::to_string(material) + "-c" + std::to_string(Compressibility) + "-alpha" + std::to_string(aDim/bDim) + "-beta" + std::to_string(bDim/thickness);
      output =  "solution";
      wn = output + "data.txt";

    }
    else if (testCase==1)
    {
      // Geometry
      bDim = 0.14; aDim = 2*bDim;
      gsInfo<<"alpha = "<<aDim/bDim<<"; beta = "<<bDim/thickness<<"\n";
      gsMultiPatch<> mp_tmp;
      mp_tmp.addPatch(gsNurbsCreator<>::BSplineRectangle(0,0,aDim/2., bDim/2.));
      mp = mp_tmp.uniformSplit();
      for (size_t p=0; p!=mp.nPatches(); p++)
      {
        gsTensorBSplineBasis<2,real_t> * basis;
          if ((basis = dynamic_cast<gsTensorBSplineBasis<2,real_t> *>(&mp.basis(p))))
          for (size_t d=0; d!=2; d++)
            basis->knots(d).transform(0,1);
      }

      gsWrite(mp,"mp");

      // Probes
      writePoints.resize(2,1);
      writePoints.col(0)<< 1.0,1.0;

      writePatches.resize(1,1);
      writePatches.row(0)<<3;

      BCs.addCondition(0,boundary::west, condition_type::dirichlet, 0, 0 ,false,0);
      BCs.addCondition(1,boundary::west, condition_type::dirichlet, 0, 0 ,false,0);

      BCs.addCondition(2,boundary::east, condition_type::collapsed, 0, 0 ,false,0);
      BCs.addCondition(3,boundary::east, condition_type::collapsed, 0, 0 ,false,0);
      BCs.addCondition(2,boundary::east, condition_type::dirichlet, 0, 0 ,false,1);
      BCs.addCondition(3,boundary::east, condition_type::dirichlet, 0, 0 ,false,1);
      BCs.addCoupled(2,boundary::east,3,boundary::east,2,0,0);

      BCs.addCondition(0,boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 2 - z.
      BCs.addCondition(2,boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 2 - z.

      real_t Load = 2*1e0;
      gsVector<> point(2); point<< 1.0, 1.0 ;
      gsVector<> load (2); load << Load,0.0;
      pLoads.addLoad(point, load, 3 );
      
      // Export
      dirname = dirname + "/QuarterSheet_-r" + std::to_string(numHref) + "-e" + std::to_string(numElevate) + "-M" + std::to_string(material) + "-c" + std::to_string(Compressibility) + "-alpha" + std::to_string(aDim/bDim) + "-beta" + std::to_string(bDim/thickness);
      output =  "solution";
      wn = output + "data.txt";

    }
    else if (testCase==2)
    {
      // Geometry
      /////////////////////////////////////////////////////////////////////////////////////
      gsMultiPatch<> mp_nurbs;
      real_t Radius = 0.25;
      bDim = 0.14; aDim = 2*bDim;
      GISMO_ASSERT(Radius<1,"Radius cannot be larger than 1");
      real_t Ri = bDim*Radius;
      real_t Ro = Ri+(bDim - Ri)/2; // middle between Ri an bDim
      gsKnotVector<> KVx (0,1,0,2) ;
      gsKnotVector<> KVy (0,1,0,3) ;
      gsMatrix<> C(6,2) ;
      C <<  Ri , 0  ,  Ro, 0
          , Ri , Ri ,  Ro, Ro
          , 0 ,  Ri ,  0 , Ro ;

      // Set weights
      gsMatrix<> ww(6,1) ;
      ww.setOnes();
      ww.at(2)= 0.707106781186548 ;
      ww.at(3)= 0.707106781186548 ;
      gsTensorNurbs<2,real_t> annulus(KVx,KVy,C,ww);
      annulus.swapDirections(0,1);
      std::vector<gsGeometry<> *> pieces = annulus.uniformSplit(0);

      gsTensorNurbs<2,real_t> bottom = static_cast<gsTensorNurbs<2,real_t> &>(*pieces[0]);
      gsTensorNurbs<2,real_t> top    = static_cast<gsTensorNurbs<2,real_t> &>(*pieces[1]);
      mp_nurbs.addPatch(bottom);
      mp_nurbs.addPatch(top);

      bottom.coefs().block(0,0,3,2)   = bottom.coefs().block(3,0,3,2);
      bottom.weights().block(0,0,3,1) = bottom.weights().block(3,0,3,1);
      bottom.coefs().block(3,0,3,1).setConstant(aDim);
      bottom.weights().block(3,0,3,1).setOnes();
      mp_nurbs.addPatch(bottom);

      top.coefs().block(0,0,3,2)   = top.coefs().block(3,0,3,2);
      top.weights().block(0,0,3,1) = top.weights().block(3,0,3,1);
      top.coefs().block(3,1,3,1).setConstant(bDim);
      top.weights().block(3,0,3,1).setOnes();
      mp_nurbs.addPatch(top);

      // we use C to store the coefs of the new patch
      C.resize(4,2);
      C.row(0) = top.coefs().row(0); // the first coefficient of top is the EV
      C.row(1)<<aDim,top.coefs()(0,1);
      C.row(2)<<top.coefs()(0,1),bDim; // the fourth coefficient of top is the second control point we need
      C.row(3)<<aDim,bDim; // the fourth coefficient of top is the second control point we need
      gsTensorBSpline<2,real_t> patch(KVx,KVx,C);
      mp_nurbs.addPatch(patch);

      // Construct mp and project mp_nurbs on it
      gsTensorBSpline<2,real_t> * tb;
      gsTensorNurbs<2,real_t> * tn;

      index_t degree = 0;
      for (size_t p=0; p!=mp_nurbs.nPatches(); p++)
      {
          tn = dynamic_cast<gsTensorNurbs<2,real_t> *>(&mp_nurbs.patch(p));
      }
      for (size_t p=0; p!=mp_nurbs.nPatches(); p++)
      {
          if ((tn = dynamic_cast<gsTensorNurbs<2,real_t> *>(&mp_nurbs.patch(p))))
              mp.addPatch(gsTensorBSpline<2,real_t>(tn->basis().knots(0),tn->basis().knots(1),tn->coefs()));
          else if ((tb = dynamic_cast<gsTensorBSpline<2,real_t> *>(&mp_nurbs.patch(p))))
              mp.addPatch(*tb);
          else
              GISMO_ERROR("Cannot construct multipatch");

          tb = dynamic_cast<gsTensorBSpline<2,real_t> *>(&mp.patch(p));
          if (tb->degree(0) > tb->degree(1))
              tb->degreeElevate(1,1);
          else if (tb->degree(0) < tb->degree(1))
              tb->degreeElevate(1,0);

          degree = math::max(degree,tb->basis().maxDegree());
      }

      for (size_t p=0; p!=mp_nurbs.nPatches(); p++)
      {
        tb = dynamic_cast<gsTensorBSpline<2,real_t> *>(&mp.patch(p));
        if (tb->basis().maxDegree() < degree)
            tb->degreeElevate(degree-tb->basis().maxDegree());
      }

      mp.computeTopology();

      // Probes
      writePoints.resize(2,1);
      writePoints.col(0)<< 1.0,1.0;

      writePatches.resize(1,1);
      writePatches.row(0)<<4;

      // BCs
      BCs.addCondition(2,boundary::north, condition_type::collapsed, 0, 0 ,false,0);
      BCs.addCondition(2,boundary::north, condition_type::dirichlet, 0, 0 ,false,1);

      BCs.addCondition(4,boundary::east, condition_type::collapsed, 0, 0 ,false,0);
      BCs.addCondition(4,boundary::east, condition_type::dirichlet, 0, 0 ,false,1);


      BCs.addCondition(1,boundary::east, condition_type::dirichlet, 0, 0 ,false,0);
      BCs.addCondition(3,boundary::east, condition_type::dirichlet, 0, 0 ,false,0);


      BCs.addCondition(0,boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 2 - z.
      BCs.addCondition(2,boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 2 - z.

      real_t Load = 1e0;
      gsVector<> point(2); point<< 1.0, 1.0 ;
      gsVector<> load (2); load << Load,0.0;
      pLoads.addLoad(point, load, 4 );
      
      // Export
      dirname = dirname + "/QuarterSheetHole_-r" + std::to_string(numHref) + "-e" + std::to_string(numElevate) + "-M" + std::to_string(material) + "-c" + std::to_string(Compressibility) + "-alpha" + std::to_string(aDim/bDim) + "-beta" + std::to_string(bDim/thickness);
      output =  "solution";
      wn = output + "data.txt";
    }
    else if (testCase==3)
    {
      // Geometry
      real_t H = 25e-3;
      real_t W  = 75e-3;
      bDim = 25e-3; aDim = 75e-3;
      mp.addPatch(gsNurbsCreator<>::BSplineRectangle(0,0,2./5.*W,H));
      mp.addPatch(gsNurbsCreator<>::BSplineRectangle(2./5.*W,0,3./5.*W,H));
      mp.addPatch(gsNurbsCreator<>::BSplineRectangle(3./5.*W,0,W,H));

      mp.degreeReduce(1);
      mp.computeTopology();

      gsWrite(mp,"mp");

      gsInfo<<"alpha = "<<aDim/bDim<<"; beta = "<<bDim/thickness<<"\n";

      // Probes
      writePoints.resize(2,1);
      writePoints.col(0)<< 0.5,0.0;

      writePatches.resize(1,1);
      writePatches.row(0)<<1;

      // BCs
      BCs.addCondition(0,boundary::west , condition_type::dirichlet, 0, 0 ,false, -1);
      BCs.addCondition(1,boundary::south, condition_type::collapsed, 0, 0, false, -1);
      BCs.addCondition(2,boundary::east , condition_type::dirichlet, 0, 0 ,false, -1);

      real_t Load = 1e0;
      gsVector<> point(2); point<< 0.5, 0.0 ;
      gsVector<> load (2); load << 0.0, -Load;
      pLoads.addLoad(point, load, 1 );

      // Export
      dirname = dirname + "/WebbBridge_-r" + std::to_string(numHref) + "-e" + std::to_string(numElevate) + "-M" + std::to_string(material) + "-c" + std::to_string(Compressibility) + "-alpha" + std::to_string(aDim/bDim) + "-beta" + std::to_string(bDim/thickness);
      output =  "solution";
      wn = output + "data.txt";
    }
    else
      GISMO_ERROR("Test case" << testCase<<" unknown.");

    mp.addAutoBoundaries();
    mp.computeTopology();

    gsDebugVar(mp);

    for(index_t i = 0; i< numElevate; ++i)
      mp.degreeElevate();    // Elevate the degree

    // h-refine
    for(index_t i = 0; i< numHref; ++i)
      mp.uniformRefine();

    mp_def = mp;
    BCs.setGeoMap(mp);

    gsMultiBasis<> dbasis(mp);
    gsInfo<<"Basis (patch 0): "<< mp.patch(0).basis() << "\n";

    gsFileManager::mkdir(dirname);

    // plot geometry
    if (plot)
      gsWriteParaview(mp,dirname + "/" + "mp",1000,true,false,"_");

    if (writeG)
    {
      gsWrite(mp,dirname + "/" + "geometry");
      gsInfo<<"Geometry written in: " + dirname + "/" + "geometry.xml\n";
    }

    if (write)
      initStepOutput(dirname + "/" + wn, writePoints);
    gsWriteParaview(dbasis.basis(0),"basis");

    gsWrite(mp,"mp");

    gsMappedBasis<2,real_t> bb2;
    gsMultiPatch<> geom = mp;
    if (smooth)
    {
      // // The approx. C1 space
      // gsApproxC1Spline<2,real_t> approxC1(mp,dbasis);
      // // approxC1.options().setSwitch("info",info);
      // // approxC1.options().setSwitch("plot",plot);
      // approxC1.options().setSwitch("interpolation",true);
      // approxC1.options().setInt("gluingDataDegree",-1);
      // approxC1.options().setInt("gluingDataSmoothness",-1);
      // approxC1.update(bb2);
      
      // gsSparseMatrix<> global2local;
      // gsAlmostC1<2,real_t> almostC1(geom);
      // almostC1.compute();
      // almostC1.matrix_into(global2local);

      // global2local = global2local.transpose();
      // geom = almostC1.exportToPatches();
      // dbasis = almostC1.localBasis();
      // bb2.init(dbasis,global2local);
      
      gsSparseMatrix<> global2local;
      gsSmoothInterfaces<2,real_t> smoothInterfaces(geom);
      smoothInterfaces.compute();
      smoothInterfaces.matrix_into(global2local);

      global2local = global2local.transpose();
      geom = smoothInterfaces.exportToPatches();
      dbasis = smoothInterfaces.localBasis();
      bb2.init(dbasis,global2local);
    }

    gsSparseSolver<>::LU solver;

    // Linear isotropic material model
    gsFunctionExpr<> force("0","0",2);
    gsConstantFunction<> t(thickness,2);
    gsConstantFunction<> E(E_modulus,2);
    gsConstantFunction<> nu(PoissonRatio,2);
    gsConstantFunction<> rho(Density,2);
    gsConstantFunction<> ratio(Ratio,2);

    mu = E_modulus / (2 * (1 + PoissonRatio));
    gsConstantFunction<> mu1(MU1,2);
    gsConstantFunction<> alpha1(ALPHA1,2);

    std::vector<gsFunction<>*> parameters;
    if (material==0) // SvK
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
      parameters.resize(4);
      parameters[0] = &E;
      parameters[1] = &nu;
      parameters[2] = &mu1;
      parameters[3] = &alpha1;
    }

    gsMaterialMatrixBase<real_t>* materialMatrix;

    gsOptionList options;
    if      (material==0 && impl==1)
    {
        parameters.resize(2);
        options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
        options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
        materialMatrix = getMaterialMatrix<2,real_t>(mp,t,parameters,rho,options);
    }
    else
    {
        options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",material);
        options.addSwitch("Compressibility","Compressibility: (false): Imcompressible | (true): Compressible",Compressibility);
        options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",impl);
        materialMatrix = getMaterialMatrix<2,real_t>(mp,t,parameters,rho,options);
    }

    gsMaterialMatrixTFT<2,real_t,false> * materialMatrixTFT = new gsMaterialMatrixTFT<2,real_t,false>(materialMatrix);
    materialMatrixTFT->options().setReal("SlackMultiplier",0);    
    // materialMatrixTFT->options().setSwitch("Explicit",true);    
    // materialMatrixTFT->updateDeformed(&mp_def);
    gsThinShellAssemblerBase<real_t>* assembler;
    if (TFT)
      assembler = new gsThinShellAssembler<2, real_t, false >(geom,dbasis,BCs,force,materialMatrixTFT);
    else
      assembler = new gsThinShellAssembler<2, real_t, false >(geom,dbasis,BCs,force,materialMatrix);


    // Construct assembler object
    if (smooth) assembler->options().setInt("Continuity",-1);
    else        assembler->options().setInt("Continuity",0);
    gsDebugVar(pLoads);
    assembler->setPointLoads(pLoads);
    if (smooth)
      assembler->setSpaceBasis(bb2);


    // Assemble linear system to obtain the force vector
    assembler->assemble();
    gsVector<> Force = assembler->rhs();

    // Function for the Jacobian
    gsStructuralAnalysisOps<real_t>::Jacobian_t Jacobian;
    gsStructuralAnalysisOps<real_t>::ALResidual_t ALResidual;
    if (smooth)
    {
      Jacobian = [&assembler,&bb2,&geom](gsVector<real_t> const &x, gsSparseMatrix<real_t> & m)
      {
        ThinShellAssemblerStatus status;
        gsMatrix<real_t> solFull = assembler->fullSolutionVector(x);
        GISMO_ASSERT(solFull.rows() % 2==0,"Rows of the solution vector does not match the number of control points");
        solFull.resize(solFull.rows()/2,2);
        gsMappedSpline<2,real_t> mspline(bb2,solFull);
        gsFunctionSum<real_t> def(&geom,&mspline);
        assembler->assembleMatrix(def);
        status = assembler->assembleMatrix(def);
        m = assembler->matrix();
        return status == ThinShellAssemblerStatus::Success;
      };
      // Function for the Residual
      ALResidual = [&assembler,&bb2,&geom,&Force](gsVector<real_t> const &x, real_t lam, gsVector<real_t> & result)
      {
        ThinShellAssemblerStatus status;
        gsMatrix<real_t> solFull = assembler->fullSolutionVector(x);
        GISMO_ASSERT(solFull.rows() % 2==0,"Rows of the solution vector does not match the number of control points");
        solFull.resize(solFull.rows()/2,2);

        gsMappedSpline<2,real_t> mspline(bb2,solFull);
        gsFunctionSum<real_t> def(&geom,&mspline);

        status = assembler->assembleVector(def);
        result = Force - lam * Force - assembler->rhs(); // assembler rhs - force = Finternal
        return status == ThinShellAssemblerStatus::Success;
      };
    }
    else
    {
      Jacobian = [&assembler,&mp_def](gsVector<real_t> const &x, gsSparseMatrix<real_t> & m)
      {
        ThinShellAssemblerStatus status;
        assembler->constructSolution(x,mp_def);
        status = assembler->assembleMatrix(mp_def);
        m = assembler->matrix();
        return status == ThinShellAssemblerStatus::Success;
      };
      // Function for the Residual
      ALResidual = [&assembler,&mp_def,&Force](gsVector<real_t> const &x, real_t lam, gsVector<real_t> & result)
      {
          ThinShellAssemblerStatus status;
          assembler->constructSolution(x,mp_def);
          status = assembler->assembleVector(mp_def);
          result = Force - lam * Force - assembler->rhs(); // assembler rhs - force = Finternal
          return status == ThinShellAssemblerStatus::Success;
      };
    }
    
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
    arcLength->options().setReal("Length",dL);
    if (method==2)
    {
      arcLength->options().setInt("AngleMethod",0); // 0: step, 1: iteration
      arcLength->options().setReal("Scaling",0.0);
    }
    arcLength->options().setSwitch("AdaptiveLength",adaptive);
    arcLength->options().setInt("AdaptiveIterations",5);
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


    gsInfo<<arcLength->options();
    arcLength->applyOptions();
    arcLength->initialize();

    gsParaviewCollection collection(dirname + "/" + output);
    gsParaviewCollection TensionFields(dirname + "/" + "tensionfield");
    gsMultiPatch<> deformation = mp;

    // Make objects for previous solutions
    real_t Lold = 0;
    gsMatrix<> Uold = Force;
    Uold.setZero();

    gsMatrix<> solVector;
    real_t indicator = 0.0;
    arcLength->setIndicator(indicator); // RESET INDICATOR
    real_t dL0 = dL;
    for (index_t k=0; k<step; k++)
    {
      gsInfo<<"Load step "<< k<<"\n";
      arcLength->step();

      if (!(arcLength->converged()))
      {
        gsInfo<<"Error: Loop terminated, arc length method did not converge.\n";
        dL = dL / 2.;
        arcLength->setLength(dL);
        arcLength->setSolution(Uold,Lold);
        k -= 1;
        continue;
      }

      indicator = arcLength->indicator();

      solVector = arcLength->solutionU();
      Uold = solVector;
      Lold = arcLength->solutionL();

      gsMappedSpline<2,real_t> mspline;
      if (plot)
      {
        gsField<> solField;
        std::string fileName = dirname + "/" + output + util::to_string(k);
        if (smooth)
        {
          /// Make a gsMappedSpline to represent the solution
          // 1. Get all the coefficients (including the ones from the eliminated BCs.)
          gsMatrix<real_t> solFull = assembler->fullSolutionVector(solVector);

          // 2. Reshape all the coefficients to a Nx3 matrix
          GISMO_ASSERT(solFull.rows() % 2==0,"Rows of the solution vector does not match the number of control points");
          solFull.resize(solFull.rows()/2,2);

          // 3. Make the mapped spline
          gsMappedSpline<2,real_t> mspline(bb2,solFull);

          // 4. Plot the mapped spline on the original geometry
          solField = gsField<>(geom, mspline,true);

          gsWriteParaview<>(solField, fileName, 1000,mesh,"_");
        }
        else
        {
          assembler->constructSolution(solVector,mp_def);

          deformation = mp_def;
          for (size_t p=0; p!=mp_def.nPatches(); p++)
            deformation.patch(p).coefs() -= mp.patch(p).coefs();// assuming 1 patch here

          if (deformed)
            solField= gsField<>(mp_def,deformation);
          else
            solField= gsField<>(mp,deformation);

          gsWriteParaview<>(solField, fileName, 1000,mesh,"_");
        }

        for (size_t p = 0; p!=mp.nPatches(); p++)
        {
          fileName = output + util::to_string(k) + "_" + util::to_string(p);
          collection.addPart(fileName + ".vts",k);
          if (mesh) collection.addPart(fileName + "_mesh.vtp",k);
        }
      }
      if (stress)
      {
        gsField<> tensionField;
        gsPiecewiseFunction<> tensionFields;
        std::string fileName;
        fileName = dirname + "/" + "tensionfield" + util::to_string(k);

        if (smooth)
        {
          /// Make a gsMappedSpline to represent the solution
          // 1. Get all the coefficients (including the ones from the eliminated BCs.)
          gsMatrix<real_t> solFull = assembler->fullSolutionVector(solVector);

          // 2. Reshape all the coefficients to a Nx3 matrix
          GISMO_ASSERT(solFull.rows() % 2==0,"Rows of the solution vector does not match the number of control points");
          solFull.resize(solFull.rows()/2,2);

          // 3. Make the mapped spline
          gsMappedSpline<2,real_t> mspline(bb2,solFull);

          // 4. Create deformation spline
          gsFunctionSum<real_t> def(&geom,&mspline);
          
          // 5. Construct stress
          assembler->constructStress(def,tensionFields,stress_type::tension_field);
          gsWriteParaview(def,tensionFields,fileName,1000,"_");
        }
        else
        {
          assembler->constructStress(mp_def,tensionFields,stress_type::tension_field);
          if (deformed)
            tensionField = gsField<>(mp_def,tensionFields, true);
          else
            tensionField = gsField<>(mp,tensionFields, true);

          gsWriteParaview( tensionField, fileName, 1000,false,"_");
        }

        for (size_t p = 0; p!=mp.nPatches(); p++)
        {
          fileName = "tensionfield" + util::to_string(k) + "_" + util::to_string(p);
          TensionFields.addPart(fileName + ".vts",k);
          if (mesh) TensionFields.addPart(fileName + "_mesh.vtp",k);
        }
      }

      // materialMatrixTFT->updateDeformed(&mp_def);


      if (write)
        writeStepOutput(arcLength,deformation, dirname + "/" + wn, writePoints,writePatches);

    }

    if (plot)
    {
      collection.save();
    }
    if (stress)
    {
      TensionFields.save();
    }

  delete materialMatrix;
  delete materialMatrixTFT;
  delete assembler;
  delete arcLength;

  return result;
}

template <class T>
void initStepOutput(const std::string name, const gsMatrix<T> & points)
{
  std::ofstream file;
  file.open(name,std::ofstream::out);
  file  << std::setprecision(20)
        << "Deformation norm" << ",";
        for (index_t k=0; k!=points.cols(); k++)
        {
          file<< "point "<<k<<" - x" << ","
              << "point "<<k<<" - y" << ",";
        }

  file  << "Lambda" << ","
        << "Indicator"
        << "\n";
  file.close();

  gsInfo<<"Step results will be written in file: "<<name<<"\n";
}

template <class T>
void writeStepOutput(const gsALMBase<T> * arcLength, const gsMultiPatch<T> & deformation, const std::string name, const gsMatrix<T> & points, const gsMatrix<index_t> & patches)
{
  gsMatrix<T> P(2,1), Q(2,1);
  gsMatrix<T> out(2,points.cols());
  gsMatrix<T> tmp;

  for (index_t p=0; p!=points.cols(); p++)
  {
    P<<points.col(p);
    deformation.patch(patches(0,p)).eval_into(P,tmp);
    out.col(p) = tmp;
  }

  std::ofstream file;
  file.open(name,std::ofstream::out | std::ofstream::app);
  file  << std::setprecision(6)
        << arcLength->solutionU().norm() << ",";
        for (index_t p=0; p!=points.cols(); p++)
        {
          file<< out(0,p) << ","
              << out(1,p) << ",";
        }

  file  << arcLength->solutionL() << ","
        << arcLength->indicator() << ","
        << "\n";
  file.close();
}

