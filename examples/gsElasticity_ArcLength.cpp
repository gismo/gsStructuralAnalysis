/** @file gsThinShell_ArcLength.cpp

    @brief Code for the arc-length method of a solid based on loads

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#include <gismo.h>

#include <gsStructuralAnalysis/gsArcLengthIterator.h>

#include <gsElasticity/gsGeoUtils.h>
#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>

using namespace gismo;

template <class T>
gsMultiPatch<T> BrickDomain(int n, int m, int o, int p, int q ,int r, T L, T B, T H);
template <class T>
gsMultiPatch<T> BrickDomain(int n, int p, T L, T B, T H);

template <class T>
void initStepOutput( const std::string name, const gsMatrix<T> & points);

template <class T>
void writeStepOutput(const gsArcLengthIterator<T> & arcLength, const gsMultiPatch<T> & deformation, const std::string name, const gsMatrix<T> & points, const index_t extreme=-1, const index_t kmax=100);

void initSectionOutput( const std::string dirname, bool undeformed=false);

template <class T>
void writeSectionOutput(const gsMultiPatch<T> & mp, const std::string dirname, const index_t coordinate=0, const T coordVal=0.0, const index_t N=100, bool undeformed=false);

int main (int argc, char** argv)
{
    // Input options
    int numElevate  = 1;
    int numHref     = 1;
    int numElevateL = -1;
    int numHrefL    = -1;
    bool plot       = false;
    bool stress       = false;
    bool membrane       = false;
    bool first  = false;
    bool SingularPoint = false;
    bool quasiNewton = false;
    int quasiNewtonInt = -1;
    bool adaptive = false;
    int step = 10;
    int method = 2; // (0: Load control; 1: Riks' method; 2: Crisfield's method; 3: consistent crisfield method; 4: extended iterations)
    bool symmetry = false;
    bool deformed = false;
    real_t perturbation = 0;

    real_t thickness = 1e-3;
    real_t E_modulus     = 1;
    real_t PoissonRatio = 0;
    real_t Density = 1e0;
    gsMultiPatch<> mp;
    real_t tau = 1e4;

    index_t Compressibility = 0;
    index_t material = 0;
    real_t Ratio = 7.0;

    real_t aDim = 2.5;
    real_t bDim = 1.0;
    real_t eta = 0;
    real_t Spring = 0;

    real_t relax = 1.0;

    int testCase = 0;

    int result = 0;

    bool write = false;

    bool THB = false;

    index_t maxit = 20;

    // Arc length method options
    real_t dL = 0; // General arc length
    real_t dLb = 0.1; // Ard length to find bifurcation
    real_t tol = 1e-6;
    real_t tolU = 1e-6;
    real_t tolF = 1e-3;

    std::string wn("data.csv");

    gsCmdLine cmd("Thin shell plate example.");

    cmd.addInt("t", "testcase", "Test case: 0: clamped-clamped, 1: pinned-pinned, 2: clamped-free", testCase);

    cmd.addInt("r","hRefine", "Number of dyadic h-refinement (bisection) steps to perform before solving", numHref);
    cmd.addInt("e","degreeElevation", "Number of degree elevation steps to perform on the Geometry's basis before solving", numElevate);
    cmd.addInt("R","hRefine2", "Number of dyadic h-refinement (bisection) steps to perform before solving (secondary direction)", numHrefL);
    cmd.addInt("E","degreeElevation2", "Number of degree elevation steps to perform on the Geometry's basis before solving (secondary direction)", numElevateL);
    cmd.addInt( "M", "Material", "Material law",  material );
    cmd.addReal("C", "MaterialRatio", "Material Ratio",  Ratio );
    cmd.addInt( "c", "Compressibility", "1: compressible, 0: incompressible",  Compressibility );

    cmd.addReal("T","hdim", "thickness of the plate", thickness);
    cmd.addReal("a","adim", "dimension a", aDim);
    cmd.addReal("b","bdim", "dimension b", bDim);

    cmd.addReal("S","spring", "Nondimensional Spring Stiffness (case 2 and 3 only!)", eta);

    cmd.addInt("m","Method", "Arc length method; 1: Crisfield's method; 2: RIks' method.", method);
    cmd.addReal("L","dLb", "arc length", dLb);
    cmd.addReal("l","dL", "arc length after bifurcation", dL);
    cmd.addReal("A","relaxation", "Relaxation factor for arc length method", relax);

    cmd.addReal("P","perturbation", "perturbation factor", perturbation);

    cmd.addReal("F","factor", "factor for bifurcation perturbation", tau);
    cmd.addInt("q","QuasiNewtonInt","Use the Quasi Newton method every INT iterations",quasiNewtonInt);
    cmd.addInt("N", "maxsteps", "Maximum number of steps", step);

    cmd.addSwitch("adaptive", "Adaptive length ", adaptive);
    cmd.addSwitch("bifurcation", "Compute singular points and bifurcation paths", SingularPoint);
    cmd.addSwitch("quasi", "Use the Quasi Newton method", quasiNewton);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addSwitch("stress", "Plot stress in ParaView format", stress);
    cmd.addSwitch("first", "Plot only first", first);
    cmd.addSwitch("write", "Write output to file", write);
    cmd.addSwitch("membrane", "Use membrane model (no bending)", membrane);
    cmd.addSwitch("symmetry", "Use symmetry boundary condition (different per problem)", symmetry);
    cmd.addSwitch("deformed", "plot on deformed shape", deformed);

    cmd.addSwitch("THB", "Use refinement", THB);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    if (dL==0)
    {
      dL = dLb;
    }

    real_t L,B,H;
    if (testCase==0)
    {
      std::string fn;
      fn = "volumes/heart3d_rotated.xml";
      gsReadFile<>(fn, mp);
    }
    else if (testCase==1 || testCase==2)
    {
        E_modulus = 1.;
        H = 0.01;
        PoissonRatio = 0.0;
        L = 1.0;
        B = 0.01;
        mp = BrickDomain(numHref, 0, 0, numElevate, 1, 1, L, B,H);
    }
    else
    {
        gsInfo<<"No geometry found\n";
    }

    gsMultiBasis<> dbasis(mp);

    if (testCase!=1 && testCase!=2)
    {
        for (index_t i = 0; i < numElevate; ++i)
        {
            dbasis.degreeElevate();
            dbasis.uniformRefine();
        }
        for (index_t i = 0; i < numHref; ++i)
            dbasis.uniformRefine();
    }

    gsInfo<<"Basis (patch 0): "<< dbasis.basis(0) << "\n";

    gsBoundaryConditions<> BCs;

    real_t displ = 1e0;
    gsFunctionExpr<> displ1("1",3);
    gsConstantFunction<> displ2(-displ,3);

    gsVector<> tmp(3);
    tmp << 0, 0, 0;
    gsConstantFunction<> neuData(tmp,3);
    gsConstantFunction<> neuData2(tmp,3);

    std::string dirname = "ArcLengthResults";
    std::string output = "solutionElasticity";

    if (testCase == 0)
    {
        tmp << 0, 0, 100;
        neuData.setValue(tmp,3);

        BCs.addCondition(0,boundary::front,condition_type::dirichlet,nullptr,0); // last number is a component (coordinate) number
        BCs.addCondition(0,boundary::front,condition_type::dirichlet,nullptr,1);
        BCs.addCondition(0,boundary::front,condition_type::dirichlet,nullptr,2);

        BCs.addCondition(0,boundary::back,condition_type::dirichlet,nullptr,0);
        BCs.addCondition(0,boundary::back,condition_type::dirichlet,nullptr,1);
        // BCs.addCondition(0,boundary::back,condition_type::dirichlet,&displ1,2);
        BCs.addCondition(0,boundary::back,condition_type::neumann,&neuData);

        BCs.addCondition(1,boundary::front,condition_type::dirichlet,nullptr,0); // last number is a component (coordinate) number
        BCs.addCondition(1,boundary::front,condition_type::dirichlet,nullptr,1);
        BCs.addCondition(1,boundary::front,condition_type::dirichlet,nullptr,2);

        BCs.addCondition(1,boundary::back,condition_type::dirichlet,nullptr,0);
        BCs.addCondition(1,boundary::back,condition_type::dirichlet,nullptr,1);
        // BCs.addCondition(1,boundary::back,condition_type::dirichlet,&displ1,2);
        BCs.addCondition(1,boundary::back,condition_type::neumann,&neuData);

        wn = output + "data.txt";
    }
    else if (testCase == 1)
    {
        real_t Load = 1e1;
        tmp << 0, 0, Load/(B*H);
        neuData.setValue(tmp,3);

        BCs.addCondition(0,boundary::west,condition_type::dirichlet,nullptr,0); // last number is a component (coordinate) number
        BCs.addCondition(0,boundary::west,condition_type::dirichlet,nullptr,1);
        BCs.addCondition(0,boundary::west,condition_type::dirichlet,nullptr,2);


        // BCs.addCondition(0,boundary::east,condition_type::dirichlet,nullptr,0);
        // BCs.addCondition(0,boundary::east,condition_type::dirichlet,nullptr,1);
        // vertical load
        BCs.addCondition(0,boundary::east,condition_type::neumann,&neuData);

        wn = output + "data.txt";
    }
    else if (testCase == 2)
    {
        real_t Load = 1e1;
        tmp << -Load/(B*H), 0, 0;
        gsInfo<<"Applied "<<tmp.transpose()<<" as load vector\n";
        neuData.setValue(tmp,3);

        BCs.addCondition(0,boundary::west,condition_type::dirichlet,nullptr,0); // last number is a component (coordinate) number
        BCs.addCondition(0,boundary::west,condition_type::dirichlet,nullptr,1);
        BCs.addCondition(0,boundary::west,condition_type::dirichlet,nullptr,2);

        BCs.addCondition(0,boundary::north,condition_type::dirichlet,nullptr,1);
        BCs.addCondition(0,boundary::south,condition_type::dirichlet,nullptr,1);

        // BCs.addCondition(0,boundary::east,condition_type::dirichlet,nullptr,0);
        // BCs.addCondition(0,boundary::east,condition_type::dirichlet,nullptr,1);
        // vertical load
        BCs.addCondition(0,boundary::east,condition_type::neumann,&neuData);

        wn = output + "data.txt";
        SingularPoint = true;
    }
    else if (testCase == 3)
    {
        real_t Load = 1e1;
        tmp << -Load/(B*H), 0, 0;
        neuData.setValue(tmp,3);

        BCs.addCondition(0,boundary::west,condition_type::dirichlet,nullptr,0); // last number is a component (coordinate) number
        BCs.addCondition(0,boundary::west,condition_type::dirichlet,nullptr,1);
        BCs.addCondition(0,boundary::west,condition_type::dirichlet,nullptr,2);

        BCs.addCondition(0,boundary::north,condition_type::dirichlet,nullptr,1);
        BCs.addCondition(0,boundary::south,condition_type::dirichlet,nullptr,1);

        // BCs.addCondition(0,boundary::east,condition_type::dirichlet,nullptr,0);
        // BCs.addCondition(0,boundary::east,condition_type::dirichlet,nullptr,1);
        // vertical load
        BCs.addCondition(0,boundary::east,condition_type::dirichlet,&displ1);

        wn = output + "data.txt";
        SingularPoint = true;
    }

    std::string commands = "mkdir -p " + dirname;
    const char *command = commands.c_str();
    system(command);

    // plot geometry
    if (plot)
      gsWriteParaview(mp,"mp",1000,true);

    gsMatrix<> writePoints(3,3);
    writePoints.col(0)<<0.5,0.5,0.5;
    writePoints.col(1)<<0,0.5,0.5;
    writePoints.col(2)<<1,0.5,0.5;

    if (write)
      initStepOutput(dirname + "/" + wn, writePoints);

    gsFunctionExpr<> g("0","0","0",3);

    gsElasticityAssembler<real_t> assembler(mp,dbasis,BCs,g);
    assembler.options().setReal("YoungsModulus",E_modulus);
    assembler.options().setReal("PoissonsRatio",PoissonRatio);
    assembler.options().setInt("MaterialLaw",material_law::hooke);
    gsInfo << "Initialized system with " << assembler.numDofs() << " dofs.\n";

    gsMultiPatch<> mp_def;

    gsStopwatch stopwatch;
    real_t time = 0.0;

    std::vector<gsMatrix<> > fixedDofs = assembler.allFixedDofs();
    typedef std::function<gsSparseMatrix<real_t> (gsVector<real_t> const &)>                                Jacobian_t;
    typedef std::function<gsVector<real_t> (gsVector<real_t> const &, real_t, gsVector<real_t> const &) >   ALResidual_t;
    // Function for the Jacobian
    Jacobian_t Jacobian = [&time,&stopwatch,&assembler,&fixedDofs](gsVector<real_t> const &x)
    {
      stopwatch.restart();
      assembler.assemble(x,fixedDofs);
      time += stopwatch.stop();

      gsSparseMatrix<real_t> m = assembler.matrix();
      // gsInfo<<"matrix = \n"<<m.toDense()<<"\n";
      return m;
    };
    // Function for the Residual
    ALResidual_t ALResidual = [&time,&stopwatch,&assembler,&fixedDofs](gsVector<real_t> const &x, real_t lam, gsVector<real_t> const &force)
    {
      stopwatch.restart();
      assembler.assemble(x,fixedDofs);
      gsVector<real_t> Fint = -(assembler.rhs() - force);
      gsVector<real_t> result = Fint - lam * force;
      time += stopwatch.stop();
      return result; // - lam * force;
    };
    // Assemble linear system to obtain the force vector
    assembler.assemble();
    gsVector<> Force = assembler.rhs();

    if (material==0)
      assembler.options().setInt("MaterialLaw",material_law::saint_venant_kirchhoff);
    else if (material==2)
      assembler.options().setInt("MaterialLaw",material_law::neo_hooke_ln);
    else if (material==22)
      assembler.options().setInt("MaterialLaw",material_law::neo_hooke_quad);

    gsArcLengthIterator<real_t> arcLength(Jacobian, ALResidual, Force);

    arcLength.options().setInt("Solver",0); // LDLT solver
    arcLength.options().setInt("BifurcationMethod",0); // 0: determinant, 1: eigenvalue
    arcLength.options().setInt("Method",method);
    arcLength.options().setReal("Length",dLb);
    arcLength.options().setInt("AngleMethod",0); // 0: step, 1: iteration
    arcLength.options().setSwitch("AdaptiveLength",adaptive);
    arcLength.options().setInt("AdaptiveIterations",5);
    arcLength.options().setReal("Perturbation",tau);
    arcLength.options().setReal("Scaling",0.0);
    arcLength.options().setReal("Tol",tol);
    arcLength.options().setReal("TolU",tolU);
    arcLength.options().setReal("TolF",tolF);
    arcLength.options().setInt("MaxIter",maxit);
    arcLength.options().setSwitch("Verbose",true);
    arcLength.options().setReal("Relaxation",relax);
    if (quasiNewtonInt>0)
    {
      quasiNewton = true;
      arcLength.options().setInt("QuasiIterations",quasiNewtonInt);
    }
    arcLength.options().setSwitch("Quasi",quasiNewton);


    gsDebug<<arcLength.options();
    arcLength.applyOptions();
    arcLength.initialize();


    gsParaviewCollection collection(dirname + "/" + output);
    gsMultiPatch<> deformation = mp;

    // Make objects for previous solutions
    real_t Lold = 0;
    gsMatrix<> Uold = Force;
    Uold.setZero();

    gsMatrix<> solVector;
    real_t indicator = 0.0;
    arcLength.setIndicator(indicator); // RESET INDICATOR
    bool bisected = false;
    real_t dLb0 = dLb;
    for (index_t k=0; k<step; k++)
    {
      gsInfo<<"Load step "<< k<<"\n";
      arcLength.step();

      if (!(arcLength.converged()))
      {
        gsInfo<<"Error: Loop terminated, arc length method did not converge.\n";
        dLb = dLb / 2.;
        arcLength.setLength(dLb);
        arcLength.setSolution(Uold,Lold);
        bisected = true;
        k -= 1;
        continue;
      }

      if (SingularPoint)
      {
        arcLength.computeStability(arcLength.solutionU(),quasiNewton);
        if (arcLength.stabilityChange())
        {
          gsInfo<<"Bifurcation spotted!"<<"\n";
          arcLength.computeSingularPoint(1e-4, 5, Uold, Lold, 1e-10, 1e-1, false);
          arcLength.switchBranch();
          dLb0 = dLb = dL;
          arcLength.setLength(dLb);
        }
      }
      indicator = arcLength.indicator();

      solVector = arcLength.solutionU();
      Uold = solVector;
      Lold = arcLength.solutionL();
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
        std::string fileName = dirname + "/" + output + util::to_string(k);
        // creating a container to plot all fields to one Paraview file
        std::map<std::string,const gsField<> *> fields;
        fields["Deformation"] = &solutionField;
        fields["Stress"] = &stressField;
        gsWriteParaviewMultiPhysics(fields,fileName,1000,true);
        fileName = output + util::to_string(k) + "0";
        collection.addTimestep(fileName,k,".vts");
        collection.addTimestep(fileName,k,"_mesh.vtp");
      }

      if (write)
        writeStepOutput(arcLength,deformation, dirname + "/" + wn, writePoints,1, 201);

      if (!bisected)
      {
        dLb = dLb0;
        arcLength.setLength(dLb);
      }
      bisected = false;

    }

    if (plot)
    {
      collection.save();
    }

  return result;
}

template <class T>
gsMultiPatch<T> BrickDomain(int n, int p, T L, T B, T H)
{
  int q = p;
  int r = p;

  int m = n;
  int o = n;
  gsMultiPatch<T> mp = BrickDomain(n, m, o, p, q, r, L, B, H);
  return mp;
}

template <class T>
gsMultiPatch<T> BrickDomain(int n, int m, int o, int p, int q ,int r, T L, T B, T H)
{
  // -------------------------------------------------------------------------
  // --------------------------Make beam geometry-----------------------------
  // -------------------------------------------------------------------------
  int dim = 3; //physical dimension
  gsKnotVector<> kv0;
  kv0.initUniform(0,1,0,p+1,1);
  gsKnotVector<> kv1;
  kv1.initUniform(0,1,0,q+1,1);
  gsKnotVector<> kv2;
  kv2.initUniform(0,1,0,r+1,1);

  for(index_t i = 0; i< n; ++i)
      kv0.uniformRefine();
  for(index_t i = 0; i< m; ++i)
      kv1.uniformRefine();
  for(index_t i = 0; i< o; ++i)
      kv2.uniformRefine();

  // Make basis
  gsTensorBSplineBasis<3,T> basis(kv0,kv1,kv2);

  // Initiate coefficient matrix
  gsMatrix<> coefs(basis.size(),dim);
  // Number of control points needed per component
  size_t len0 = basis.component(0).size();
  size_t len1 = basis.component(1).size();
  size_t len2 = basis.component(2).size();
  // Uniformly distribute control points per component
  gsVector<> coefvec0(len0);
  coefvec0.setLinSpaced(len0,0.0,L);
  gsVector<> coefvec1(basis.component(1).size());
  coefvec1.setLinSpaced(len1,0.0,B);
  gsVector<> coefvec2(basis.component(2).size());
  coefvec2.setLinSpaced(len2,0.0,H);

  // Define a matrix with ones
  gsVector<> temp(len0);
  temp.setOnes();
  for (index_t l = 0; l < len2; l++)
    {
        for (index_t k = 0; k < len1; k++)
        {
            index_t offset = l*len0*len1;
            // First column contains x-coordinates (length)
            coefs.col(0).segment(k*len0+offset,len0) = coefvec0;
            // Second column contains y-coordinates (width)
            coefs.col(1).segment(k*len0+offset,len0) = temp*coefvec1.at(k);

            coefs.col(2).segment(k*len0+offset,len0) = temp*coefvec2.at(l);
        }
    }
  // gsInfo<<"\n"<<coefs<<"\n";
  // Create gsGeometry-derived object for the patch
  gsTensorBSpline<3,real_t> shape(basis,coefs);

  gsMultiPatch<T> mp;
  mp.addPatch(shape);
  mp.addAutoBoundaries();

  return mp;
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
              << "point "<<k<<" - y" << ","
              << "point "<<k<<" - z" << ",";
        }

  file  << "Lambda" << ","
        << "Indicator"
        << "\n";
  file.close();

  gsInfo<<"Step results will be written in file: "<<name<<"\n";
}

template <class T>
void writeStepOutput(const gsArcLengthIterator<T> & arcLength, const gsMultiPatch<T> & deformation, const std::string name, const gsMatrix<T> & points, const index_t extreme, const index_t kmax) // extreme: the column of point indices to compute the extreme over (default -1)
{
  gsMatrix<T> P(2,1), Q(2,1);
  gsMatrix<T> out(3,points.cols());
  gsMatrix<T> tmp;

  for (index_t p=0; p!=points.cols(); p++)
  {
    P<<points.col(p);
    deformation.patch(0).eval_into(P,tmp);
    out.col(p) = tmp;
  }

  std::ofstream file;
  file.open(name,std::ofstream::out | std::ofstream::app);
  if (extreme==-1)
  {
    file  << std::setprecision(6)
          << arcLength.solutionU().norm() << ",";
          for (index_t p=0; p!=points.cols(); p++)
          {
            file<< out(0,p) << ","
                << out(1,p) << ","
                << out(2,p) << ",";
          }

    file  << arcLength.solutionL() << ","
          << arcLength.indicator() << ","
          << "\n";
  }
  else if (extreme==0 || extreme==1)
  {
    gsMatrix<T> out2(kmax,points.cols()); // evaluation points in the rows, output (per coordinate) in columns
    for (int p = 0; p != points.cols(); p ++)
    {
      Q.at(1-extreme) = points(1-extreme,p);
      for (int k = 0; k != kmax; k ++)
      {
        Q.at(extreme) = 1.0*k/(kmax-1);
        deformation.patch(0).eval_into(Q,tmp);
        out2(k,p) = tmp.at(2); // z coordinate
      }
    }

    file  << std::setprecision(6)
          << arcLength.solutionU().norm() << ",";
          for (index_t p=0; p!=points.cols(); p++)
          {
            file<< out(0,p) << ","
                << out(1,p) << ","
                << std::max(abs(out2.col(p).maxCoeff()),abs(out2.col(p).minCoeff())) << ",";
          }

    file  << arcLength.solutionL() << ","
          << arcLength.indicator() << ","
          << "\n";
  }
  else
    GISMO_ERROR("Extremes setting unknown");

  file.close();
}

void initSectionOutput(const std::string dirname, bool undeformed)
{
  std::ofstream file2, file3, file4;
  std::string wn2,wn3,wn4;

  if (! undeformed)
  {
    wn2 = dirname + "/" + "pointdataX.txt";
    wn3 = dirname + "/" + "pointdataY.txt";
    wn4 = dirname + "/" + "pointdataZ.txt";
  }
  else
  {
    wn2 = dirname + "/" + "pointdataX0.txt";
    wn3 = dirname + "/" + "pointdataY0.txt";
    wn4 = dirname + "/" + "pointdataZ0.txt";
  }

  file2.open(wn2,std::ofstream::out);
  file2.close();

  file3.open(wn3,std::ofstream::out);
  file3.close();

  file4.open(wn4,std::ofstream::out);
  file4.close();

  gsInfo<<"Cross-section results will be written in directory: "<<dirname<<"\n";
}

template <class T>
void writeSectionOutput(const gsMultiPatch<T> & mp, const std::string dirname, const index_t coordinate, const T coordVal, const index_t N, bool undeformed) // coordinate: the column which remains constant at coordVal
{
  gsMatrix<T> P(2,1);
  gsMatrix<T> tmp;
  P.setZero();
  P.at(coordinate) = coordVal;

  std::ofstream file2, file3, file4;
  std::string wn2,wn3,wn4;

  if (! undeformed)
  {
    wn2 = dirname + "/" + "pointdataX.txt";
    wn3 = dirname + "/" + "pointdataY.txt";
    wn4 = dirname + "/" + "pointdataZ.txt";
  }
  else
  {
    wn2 = dirname + "/" + "pointdataX0.txt";
    wn3 = dirname + "/" + "pointdataY0.txt";
    wn4 = dirname + "/" + "pointdataZ0.txt";
  }

  file2.open(wn2,std::ofstream::out | std::ofstream::app);
  file3.open(wn3,std::ofstream::out | std::ofstream::app);
  file4.open(wn4,std::ofstream::out | std::ofstream::app);


  gsMatrix<T> out(3,N); // evaluation points in the rows, output (per coordinate) in columns
    for (int k = 0; k != N; k ++)
    {
      P.at(1-coordinate) = 1.0*k/(N-1);

      mp.patch(0).eval_into(P,tmp);
      out.col(k) = tmp; // z coordinate

      std::string str2 = std::to_string(out(0,k));
      std::string str3 = std::to_string(out(1,k));
      std::string str4 = std::to_string(out(2,k));
      if(k+1 == N)
      {
          file2<<str2;
          file3<<str3;
          file4<<str4;
      }
      else{
          file2<<str2<<',';
          file3<<str3<<',';
          file4<<str4<<',';
      }
    }
    file2<<'\n';
    file2.close();
    file3<<'\n';
    file3.close();
    file4<<'\n';
    file4.close();
}
