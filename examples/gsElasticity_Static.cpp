/** @file gsElasticity_Static.cpp

    @brief Static simulations of a solid

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#include <gismo.h>

#include <gsElasticity/gsGeoUtils.h>
#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>

#include <gsStructuralAnalysis/gsStaticSolver.h>


using namespace gismo;

template <class T>
gsMultiPatch<T> BrickDomain(int n, int m, int o, int p, int q ,int r, T L, T B, T H);
template <class T>
gsMultiPatch<T> BrickDomain(int n, int p, T L, T B, T H);

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

int main (int argc, char** argv)
{
    // Input options
    int numElevate  = 1;
    int numHref     = 1;
    bool plot       = false;
    bool nonlinear  = false;

    index_t material = 0;

    int Dirichlet = 1; // 1 = elimination; 2 = nitsche

    real_t E_modulus     = 1;
    real_t PoissonRatio = 0;
    gsMultiPatch<> mp;

    int testCase = 0;

    bool write = false;

    std::string wn("data.csv");

    gsCmdLine cmd("Thin shell plate example.");
    cmd.addInt("r","hRefine",
       "Number of dyadic h-refinement (bisection) steps to perform before solving",
       numHref);
    cmd.addInt("t", "testcase",
        "Test case: 0: clamped-clamped, 1: pinned-pinned, 2: clamped-free",
       testCase);
    cmd.addInt("e","degreeElevation",
      "Number of degree elevation steps to perform on the Geometry's basis before solving",
      numElevate);
    cmd.addInt( "M", "Material", "Material law",  material );
    cmd.addInt("D","DirichletStrategy",
      "Dirichlet Strategy to be used; 1: elimination; 2: nitsche",
      Dirichlet);
    cmd.addSwitch("nl", "Nonlinear elasticity (otherwise linear)", nonlinear);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addSwitch("write", "Write convergence data to file", write);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    real_t L,B,H;
    gsMultiBasis<> dbasis;
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
        mp = BrickDomain(numHref, 0, 0, 1, 1, 1, L, B,H);
        dbasis = gsMultiBasis<>(mp);
        dbasis.degreeIncrease(numElevate,0);
    }
    else if (testCase==3)
    {
        E_modulus = 1e6;
        H = 0.005;
        PoissonRatio = 0.3;
        L = 1.0;
        B = 1.0;
        mp = BrickDomain(numHref, numHref, 0, 1, 1, 1, L, B,H);
        dbasis = gsMultiBasis<>(mp);
        dbasis.degreeIncrease(numElevate,0);
        dbasis.degreeIncrease(numElevate,1);
    }
    else if (testCase==4)
    {
        E_modulus = 1.;
        H = 0.14;
        PoissonRatio = 0.4999;
        L = 140;
        B = 70;
        mp = BrickDomain(numHref, numHref, 0, numElevate, numElevate, 2, L, B,H);
        dbasis = gsMultiBasis<>(mp);
        dbasis.degreeIncrease(numElevate,0);
        dbasis.degreeIncrease(numElevate,1);
    }
    else
    {
        gsInfo<<"No geometry found\n";
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

    gsVector<> bodyforce(3);
    bodyforce << 0, 0, 0;

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
    }
    else if (testCase == 3)
    {
        // tmp << 0, 0, -1;
        // gsInfo<<"Applied "<<tmp.transpose()<<" as load vector\n";
        // neuData.setValue(tmp,3);

        BCs.addCondition(0,boundary::north,condition_type::dirichlet,nullptr,0);
        BCs.addCondition(0,boundary::north,condition_type::dirichlet,nullptr,1);
        BCs.addCondition(0,boundary::north,condition_type::dirichlet,nullptr,2);

        BCs.addCondition(0,boundary::east,condition_type::dirichlet,nullptr,0);
        BCs.addCondition(0,boundary::east,condition_type::dirichlet,nullptr,1);
        BCs.addCondition(0,boundary::east,condition_type::dirichlet,nullptr,2);

        BCs.addCondition(0,boundary::south,condition_type::dirichlet,nullptr,0);
        BCs.addCondition(0,boundary::south,condition_type::dirichlet,nullptr,1);
        BCs.addCondition(0,boundary::south,condition_type::dirichlet,nullptr,2);

        BCs.addCondition(0,boundary::west,condition_type::dirichlet,nullptr,0);
        BCs.addCondition(0,boundary::west,condition_type::dirichlet,nullptr,1);
        BCs.addCondition(0,boundary::west,condition_type::dirichlet,nullptr,2);

        // BCs.addCondition(0,boundary::front,condition_type::neumann,&neuData);

        bodyforce<<0,0,-1/H;

        wn = output + "data.txt";
    }
    else if (testCase == 4)
    {
        real_t Load = 1e0;
        tmp << Load/(B*H), 0, 0;
        neuData.setValue(tmp,3);

        BCs.addCondition(0,boundary::west,condition_type::dirichlet,nullptr,0); // last number is a component (coordinate) number

        // BCs.addCondition(0,boundary::east,condition_type::collapsed,nullptr,0); // last number is a component (coordinate) number
        BCs.addCondition(0,boundary::east,condition_type::dirichlet,nullptr,1); // last number is a component (coordinate) number
        BCs.addCondition(0,boundary::east,condition_type::dirichlet,nullptr,2); // last number is a component (coordinate) number

        BCs.addCondition(0,boundary::south,condition_type::dirichlet,nullptr,1);
        BCs.addCondition(0,boundary::south,condition_type::dirichlet,nullptr,2);

        BCs.addCondition(0,boundary::east,condition_type::dirichlet,&displ1);

        wn = output + "data.txt";
    }

    // plot geometry
    if (plot)
      gsWriteParaview(mp,"mp",1000,true);

    gsConstantFunction<> g(bodyforce,3);


    // -----------------------------------------------------------------------------------------------
    // --------------------------------------Solve Static problem-------------------------------------
    // -----------------------------------------------------------------------------------------------

    gsElasticityAssembler<real_t> assembler(mp,dbasis,BCs,g);
    assembler.options().setReal("YoungsModulus",E_modulus);
    assembler.options().setReal("PoissonsRatio",PoissonRatio);
    assembler.options().setInt("MaterialLaw",material_law::hooke);
    gsInfo << "Initialized system with " << assembler.numDofs() << " dofs.\n";

    gsMultiPatch<> mp_def;

    gsStopwatch stopwatch,stopwatch2;
    // Define Matrices
    stopwatch.restart();
    stopwatch2.restart();
    real_t time = 0.0;
    real_t totaltime = 0.0;

    std::vector<gsMatrix<> > fixedDofs = assembler.allFixedDofs();
    typedef std::function<gsSparseMatrix<real_t> (gsVector<real_t> const &)>    Jacobian_t;
    typedef std::function<gsVector<real_t> (gsVector<real_t> const &) >         Residual_t;
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
    Residual_t Residual = [&time,&stopwatch,&assembler,&fixedDofs](gsVector<real_t> const &x)
    {
      stopwatch.restart();
      assembler.assemble(x,fixedDofs);
      gsVector<real_t> result = assembler.rhs();
      time += stopwatch.stop();
      return result; // - lam * force;
    };

    switch (Dirichlet)
    {
      case 1:
        assembler.options().setInt("DirichletStrategy", dirichlet::elimination);
        break;
      case 2:
        assembler.options().setInt("DirichletStrategy", dirichlet::nitsche);
        gsInfo<<"Nitsche!\n";
        break;
      default:
        gsInfo<<"no strategy known... using elimination\n";
        assembler.options().setInt("DirichletStrategy", dirichlet::elimination);
        break;
    }
    // assembler.options().setInt("Verbosity",newton_verbosity::all);


    gsInfo << "Initialized system with " << assembler.numDofs() << " dofs.\n";

    // Assemble linear system to obtain the force vector
    assembler.assemble();
    time += stopwatch.stop();
    gsVector<> Force = assembler.rhs();

    if (material==0)
      assembler.options().setInt("MaterialLaw",material_law::saint_venant_kirchhoff);
    else if (material==2)
      assembler.options().setInt("MaterialLaw",material_law::neo_hooke_ln);
    else if (material==22)
      assembler.options().setInt("MaterialLaw",material_law::neo_hooke_quad);

    //=============================================//
                  // Solving //
    //=============================================//

    gsSparseMatrix<> matrix = assembler.matrix();
    gsVector<> vector = assembler.rhs();

    // Configure Structural Analsysis module
    gsStaticSolver<real_t> staticSolver(matrix,vector,Jacobian,Residual);
    gsOptionList solverOptions = staticSolver.options();
    solverOptions.setInt("Verbose",1);
    solverOptions.setInt("MaxIterations",10);
    solverOptions.setReal("Tolerance",1e-6);
    staticSolver.setOptions(solverOptions);

    gsInfo << "Solving...\n";
    // Solve linear problem
    gsVector<> solVector;
    solVector = staticSolver.solveLinear();
    if (nonlinear)
        solVector = staticSolver.solveNonlinear();

    totaltime += stopwatch2.stop();

    if (plot)
    {
      // solution to the nonlinear problem as an isogeometric displacement field
      gsMultiPatch<> displacement;
      assembler.constructSolution(solVector,assembler.allFixedDofs(),displacement);
      gsPiecewiseFunction<> stresses;
      assembler.constructCauchyStresses(displacement,stresses,stress_components::von_mises);

      // constructing an IGA field (geometry + solution)
      gsField<> displacementField(assembler.patches(),displacement);
      gsField<> stressField(assembler.patches(),stresses,true);
      // creating a container to plot all fields to one Paraview file
      std::map<std::string,const gsField<> *> fields;
      fields["Displacement"] = &displacementField;
      fields["von Mises"] = &stressField;
      gsWriteParaviewMultiPhysics(fields,"solution",1000,true);
    }

    return 1;
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
