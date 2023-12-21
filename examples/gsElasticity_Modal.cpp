/** @file gsElasticity_Static.cpp

    @brief Static simulations of a solid

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#include <gismo.h>

#ifdef gsElasticity_ENABLED
#include <gsElasticity/gsGeoUtils.h>
#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsMassAssembler.h>
#endif

#include <gsStructuralAnalysis/src/gsEigenSolvers/gsModalSolver.h>


using namespace gismo;

#ifdef gsElasticity_ENABLED

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
    bool first       = false;
    bool sparse       = false;
    bool nonlinear  = false;

    index_t material = 0;

    int Dirichlet = 1; // 1 = elimination; 2 = nitsche

    real_t E_modulus     = 1e0;
    real_t PoissonRatio = 0;
    real_t Density = 1e0;
    gsMultiPatch<> mp;

    real_t shift = 0.0;

    int testCase = 1;

    bool write = false;

    std::string wn("data.csv");

    gsCmdLine cmd("Elasticity shell plate example.");
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
    cmd.addReal( "s", "shift", "shift",  shift );
    cmd.addSwitch("nl", "Nonlinear elasticity (otherwise linear)", nonlinear);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addSwitch("first", "Plot only the first mode", first);
    cmd.addSwitch("sparse", "Compute using sparse solver", sparse);
    cmd.addSwitch("write", "Write convergence data to file", write);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    real_t L,B,H;
    if (testCase==1)
    {
        E_modulus = 1.;
        H = 0.01;
        PoissonRatio = 0.3;
        L = 1.0;
        B = 1.0;
        mp = BrickDomain(numHref,numHref,0,numElevate,numElevate,1, L, B,H);
    }
    else
    {
        gsInfo<<"No geometry found\n";
    }

    gsMultiBasis<> dbasis(mp);

    // if (testCase!=1 && testCase!=2 && testCase!=3)
    // {
    //     for (index_t i = 0; i < numElevate; ++i)
    //     {
    //         dbasis.degreeElevate();
    //         dbasis.uniformRefine();
    //     }
    //     for (index_t i = 0; i < numHref; ++i)
    //         dbasis.uniformRefine();
    // }

    gsInfo<<"Basis (patch 0): "<< dbasis.basis(0) << "\n";

    gsBoundaryConditions<> BCs;

    if (testCase == 1)
    {
        BCs.addCondition(0,boundary::north,condition_type::dirichlet,nullptr,0); // last number is a component (coordinate) number
        BCs.addCondition(0,boundary::north,condition_type::dirichlet,nullptr,1);
        BCs.addCondition(0,boundary::north,condition_type::dirichlet,nullptr,2);

        BCs.addCondition(0,boundary::east,condition_type::dirichlet,nullptr,0); // last number is a component (coordinate) number
        BCs.addCondition(0,boundary::east,condition_type::dirichlet,nullptr,1);
        BCs.addCondition(0,boundary::east,condition_type::dirichlet,nullptr,2);

        BCs.addCondition(0,boundary::south,condition_type::dirichlet,nullptr,0); // last number is a component (coordinate) number
        BCs.addCondition(0,boundary::south,condition_type::dirichlet,nullptr,1);
        BCs.addCondition(0,boundary::south,condition_type::dirichlet,nullptr,2);

        BCs.addCondition(0,boundary::west,condition_type::dirichlet,nullptr,0); // last number is a component (coordinate) number
        BCs.addCondition(0,boundary::west,condition_type::dirichlet,nullptr,1);
        BCs.addCondition(0,boundary::west,condition_type::dirichlet,nullptr,2);
    }

    // plot geometry
    if (plot)
      gsWriteParaview(mp,"mp",1000,true);

    gsFunctionExpr<> g("0","0","0",3);


    // -----------------------------------------------------------------------------------------------
    // --------------------------------------Solve Static problem-------------------------------------
    // -----------------------------------------------------------------------------------------------

    gsElasticityAssembler<real_t> assembler(mp,dbasis,BCs,g);
    assembler.options().setReal("YoungsModulus",E_modulus);
    assembler.options().setReal("PoissonsRatio",PoissonRatio);
    assembler.options().setInt("MaterialLaw",material_law::hooke);
    gsInfo << "Initialized system with " << assembler.numDofs() << " dofs.\n";

    gsMultiPatch<> mp_def;

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

    gsMassAssembler<real_t> assemblerMass(mp,dbasis,BCs,g);
    assemblerMass.options() = assembler.options();
    assemblerMass.options().addReal("Density","Density of the material",1.);

    // Configure Structural Analsysis module
    assembler.assemble();
    gsSparseMatrix<> K =  assembler.matrix();
    assemblerMass.assemble();
    gsSparseMatrix<> M =  assemblerMass.matrix();

    gsModalSolver<real_t> modal(K,M);
    modal.verbose();

    if (!sparse)
      modal.compute();
    else
      modal.computeSparse(shift,10);


    gsMatrix<> values = modal.values();
    gsMatrix<> vectors = modal.vectors();

    gsInfo<< "First 10 eigenvalues:\n";
    for (index_t k = 0; k<10; k++)
        gsInfo<<"\t"<<std::setprecision(20)<<values.at(k)<<"\n";
    gsInfo<<"\n";

    if (plot)
    {
      gsInfo<<"Plotting in Paraview...\n";
      system("mkdir -p ModalResults");
      gsMatrix<> modeShape;
      gsMultiPatch<> displacement;
      gsParaviewCollection collection("ModalResults/modes_solid");

      int N = 1;
      if (!first)
        N = vectors.cols();
      for (index_t m=0; m<N; m++)
      {

        // Compute solution based on eigenmode with number 'mode'
        modeShape = modal.vector(m);//solver.solve( assembler->rhs() );

        // solution to the nonlinear problem as an isogeometric displacement field
        assembler.constructSolution(modeShape,assembler.allFixedDofs(),displacement);

        // Normalize mode shape amplitude in z coordinate
        real_t maxAmpl = std::max(math::abs(displacement.patch(0).coefs().col(2).maxCoeff()),math::abs(displacement.patch(0).coefs().col(2).minCoeff()));
        if (maxAmpl!=0.0)
        {
          displacement.patch(0).coefs() = displacement.patch(0).coefs()/maxAmpl;
        }

        // constructing an IGA field (geometry + solution)
        gsField<> displacementField(assembler.patches(),displacement);

        std::string fileName = "ModalResults/modes_solid" + util::to_string(m);
        gsWriteParaview(displacementField,fileName,5000,true);
        fileName = "modes_solid" + util::to_string(m) + "0";
        collection.addTimestep(fileName,m,".vts");
      }

      collection.save();
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

#else

int main (int argc, char** argv)
{
  gsInfo<<"To run this example, compile G+Smo with gsElasticity\n";
  return EXIT_SUCCESS;
}

#endif

