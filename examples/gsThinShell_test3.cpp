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
#include <gsThinShell2/gsMaterialMatrix.h>

#include <gsElasticity/gsWriteParaviewMultiPhysics.h>

#include <gsStructuralAnalysis/gsStaticSolver.h>

//#include <gsThinShell/gsNewtonIterator.h>

using namespace gismo;

template <class T>
gsMultiPatch<T> RectangularDomain(int n, int m, int p, int q, T L, T B);
template <class T>
gsMultiPatch<T> RectangularDomain(int n, int p, T L, T B);

template <class T>
gsMultiPatch<T> Rectangle(T L, T B);

template <class T>
gsMultiPatch<T> RectangularDomainVert(int n, int m, int p, int q, T L, T B);
template <class T>
gsMultiPatch<T> RectangularDomainVert(int n, int p, T L, T B);

template <class T>
gsMultiPatch<T> RectangularDomain90(int n, int m, int p, int q, T L, T B);
template <class T>
gsMultiPatch<T> RectangularDomain90(int n, int p, T L, T B);

template <class T>
gsMultiPatch<T> FrustrumDomain(int n, int p, T R1, T R2, T h);

// Choose among various shell examples, default = Thin Plate
int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot  = false;
    bool stress= false;
    index_t numRefine  = 1;
    index_t numElevate = 1;
    index_t Compressibility = 0;
    index_t material = 0;
    index_t testCase = -1;
    bool nonlinear = false;
    bool verbose = false;
    std::string fn1,fn2;
    fn1 = "planar/unitplate.xml";
    fn2 = "pde/shells/kirchhoff_shell1.xml";
    bool membrane = false;


    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving",  numRefine );
    cmd.addInt( "m", "Material", "Material law",  material );
    cmd.addInt( "c", "Compressibility", "1: compressible, 0: incompressible",  Compressibility );
    cmd.addInt( "t", "testCase", "Define test case",  testCase );
    cmd.addString( "f", "GEOMfile", "Input XML Geometry file", fn1 );
    cmd.addString( "F", "PDEfile", "Input XML PDE file", fn2 );
    cmd.addSwitch("nl", "Solve nonlinear problem", nonlinear);
    cmd.addSwitch("verbose", "Full matrix and vector output", verbose);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("stress", "Create a ParaView visualization file with the stresses", stress);
    cmd.addSwitch("membrane", "Use membrane model (no bending)", membrane);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Read input file]
    gsMultiPatch<> mp, mp_def, deformation;
    gsFileData<> fd;

    if (testCase == 1)
    {
        fn1 = "planar/unitplate.xml";
        fn2 = "pde/shells/kirchhoff_shell1.xml";
    }
    else if (testCase == 2)
    {
        fn1 = "scordelis_lo_roof.xml";
        fn2 = "pde/shells/kirchhoff_shell_scordelis.xml";
    }

    gsReadFile<>(fn1, mp);
    fd.read(fn2);
    gsInfo << "Loaded file "<< fd.lastPath() <<"\n";

    // Boundary conditions
    gsBoundaryConditions<> bc;
    fd.getId(20, bc); // id=2: boundary conditions
    bc.setGeoMap(mp);
    gsInfo<<"Boundary conditions:\n"<< bc <<"\n";

    // Loads
    gsFunctionExpr<> force, pressure;
    fd.getId(21, force); // id=1: source function
    gsInfo<<"Force function "<< force << "\n";
    // fd.getId(22, pressure); // id=1: source function ------- TO DO!
    // gsInfo<<"Pressure function "<< force << "\n";

    // Loads
    gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();
    gsMatrix<> points,loads;
    fd.getId(30,points);
    fd.getId(31,loads);
    for (index_t k =0; k!=points.cols(); k++)
        pLoads.addLoad(points.col(k), loads.col(k), 0 ); // in parametric domain!

    // Material properties
    gsFunctionExpr<> t,E,nu,rho;
    fd.getId(10,t);
    gsInfo<<"thickness =  "<< t << "\n";
    fd.getId(11,E);
    fd.getId(12,nu);
    fd.getId(13,rho);

    gsFunctionExpr<> ratio;
    if (material==3 || material==13 || material == 23)
    {
        fd.getId(14,ratio);
    }
    gsFunctionExpr<> alpha1, mu1, alpha2, mu2;
    if (material==14)
    {
        fd.getId(15,alpha1);
        fd.getId(16,mu1);
        fd.getId(17,alpha1);
        fd.getId(18,mu1);
    }

    std::vector<gsFunction<>*> parameters;
    if (material==2 || material==12 ||  material==22 || material==0)
    {
        parameters.resize(2);
        parameters[0] = &E;
        parameters[1] = &nu;
    }
    else if (material==3 || material==13 || material == 23)
    {
        parameters.resize(3);
        parameters[0] = &E;
        parameters[1] = &nu;
        parameters[2] = &ratio;
    }
    else if (material==14)
    {
        parameters.resize(6);
        parameters[0] = &E;
        parameters[1] = &nu;
        parameters[2] = &mu1;
        parameters[3] = &alpha1;
        parameters[4] = &mu2;
        parameters[5] = &alpha2;
    }
    //! [Read input file]

    // p-refine
    if (numElevate!=0)
        mp.degreeElevate(numElevate);

    // h-refine
    for (int r =0; r < numRefine; ++r)
        mp.uniformRefine();

    // set initial deformation to undeformed state
    mp_def = mp;

    // define basis
    gsMultiBasis<> dbasis(mp);

    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
    gsInfo<<mp_def<<"\n";
    gsInfo << dbasis.basis(0)<<"\n";

    if (plot)
        gsWriteParaview<>( mp_def    , "mp", 1000, true);


    // Set MaterialMatrix
    gsMaterialMatrix materialMatrix(mp,mp_def,t,parameters,rho);
    materialMatrix.options().setInt("MaterialLaw",material);
    materialMatrix.options().setInt("Compressibility",Compressibility);

    // Set Shell Assembler
    gsThinShellAssembler assembler(mp,dbasis,bc,force,materialMatrix);
    assembler.setPointLoads(pLoads);
    if (membrane)
        assembler.setMembrane();
    // if (pressure!= 0.0)
    //     assembler.setPressure(pressFun);

    // Define Matrices
    assembler.assemble();
    gsSparseMatrix<> matrix = assembler.matrix();
    gsVector<> vector = assembler.rhs();
    // Function for the Jacobian
    std::function<gsSparseMatrix<real_t> (gsVector<real_t> const &)> Jacobian;
    Jacobian = [&assembler,&mp_def](gsVector<real_t> const &x)
    {
      assembler.constructSolution(x,mp_def);
      assembler.assembleMatrix(mp_def);
      gsSparseMatrix<real_t> m = assembler.matrix();
      return m;
    };
    // Function for the Residual
    std::function<gsVector<real_t> (gsVector<real_t> const &) > Residual;
    Residual = [&assembler,&mp_def](gsVector<real_t> const &x)
    {
      assembler.constructSolution(x,mp_def);
      assembler.assembleVector(mp_def);
      return assembler.rhs();
    };

    // Configure Structural Analsysis module
    gsOptionList solverOptions;
    fd.getId(90, solverOptions); // id=4: assembler options
    gsStaticSolver<real_t> staticSolver(matrix,vector,Jacobian,Residual);
    staticSolver.setOptions(solverOptions);

    // Solve linear problem
    gsVector<> solVector;
    if (!nonlinear)
    {
        solVector = staticSolver.solveLinear();

        mp_def = assembler.constructSolution(solVector);
        deformation = mp_def;
        for (size_t k = 0; k != mp_def.nPatches(); ++k)
            deformation.patch(k).coefs() -= mp.patch(k).coefs();

        if (plot)
        {
            gsInfo <<"Maximum deformation coef: "
                   << deformation.patch(0).coefs().colwise().maxCoeff() <<".\n";
            gsInfo <<"Minimum deformation coef: "
                   << deformation.patch(0).coefs().colwise().minCoeff() <<".\n";


            gsField<> solField(mp_def, deformation);
            gsInfo<<"Plotting in Paraview...\n";
            gsWriteParaview<>( solField, "solution", 1000, true);

            gsWriteParaview<>( mp_def, "mp_def", 1000, true);

         }
        if (stress)
        {
            gsPiecewiseFunction<> membraneStresses;
            assembler.constructStress(mp_def,membraneStresses,stress_type::membrane);
            gsField<> membraneStress(mp_def,membraneStresses, true);

            gsPiecewiseFunction<> flexuralStresses;
            assembler.constructStress(mp_def,flexuralStresses,stress_type::flexural);
            gsField<> flexuralStress(mp_def,flexuralStresses, true);

            gsPiecewiseFunction<> stretches;
            assembler.constructStress(mp_def,stretches,stress_type::principal_stretch);
            gsField<> Stretches(mp_def,stretches, true);

            // gsPiecewiseFunction<> membraneStresses_p;
            // assembler.constructStress(mp_def,membraneStresses_p,stress_type::principal_stress_membrane);
            // gsField<> membraneStress_p(mp_def,membraneStresses_p, true);

            // gsPiecewiseFunction<> flexuralStresses_p;
            // assembler.constructStress(mp_def,flexuralStresses_p,stress_type::principal_stress_flexural);
            // gsField<> flexuralStress_p(mp_def,flexuralStresses_p, true);

            gsPiecewiseFunction<> stretch1;
            assembler.constructStress(mp_def,stretch1,stress_type::principal_stretch_dir1);
            gsField<> stretchDir1(mp_def,stretch1, true);

            gsPiecewiseFunction<> stretch2;
            assembler.constructStress(mp_def,stretch2,stress_type::principal_stretch_dir2);
            gsField<> stretchDir2(mp_def,stretch2, true);

            gsPiecewiseFunction<> stretch3;
            assembler.constructStress(mp_def,stretch3,stress_type::principal_stretch_dir3);
            gsField<> stretchDir3(mp_def,stretch3, true);


            gsField<> solutionField(mp,deformation, true);


            // gsField<> stressField = assembler.constructStress(mp_def,stress_type::membrane_strain);

            std::map<std::string,const gsField<> *> fields;
            fields["Deformation"] = &solutionField;
            fields["Membrane Stress"] = &membraneStress;
            fields["Flexural Stress"] = &flexuralStress;
            fields["Principal Stretch"] = &Stretches;
            // fields["Principal Membrane Stress"] = &membraneStress_p;
            // fields["Principal Flexural Stress"] = &flexuralStress_p;
            fields["Principal Direction 1"] = &stretchDir1;
            fields["Principal Direction 2"] = &stretchDir2;
            fields["Principal Direction 3"] = &stretchDir3;

            gsWriteParaviewMultiPhysics(fields,"stress",500,true);

        }
    }
    else
    {
        solVector = staticSolver.solveNonlinear();
        mp_def = assembler.constructSolution(solVector);

        deformation = mp_def;
        for (size_t k = 0; k != mp_def.nPatches(); ++k)
            deformation.patch(k).coefs() -= mp.patch(k).coefs();

        // ! [Export visualization in ParaView]
        if ( (plot) && (nonlinear) )
        {
            gsField<> solField(mp_def, deformation);
            gsInfo<<"Plotting in Paraview...\n";
            gsWriteParaview<>( solField, "solution", 1000, true);
            // ev.options().setSwitch("plot.elements", true);
            // ev.writeParaview( u_sol   , G, "solution");

            // gsFileManager::open("solution.pvd");

            gsInfo <<"Maximum deformation coef: "
                   << deformation.patch(0).coefs().colwise().maxCoeff() <<".\n";
            gsInfo <<"Minimum deformation coef: "
                   << deformation.patch(0).coefs().colwise().minCoeff() <<".\n";
        }
        if ((stress) && (nonlinear))
        {

            gsPiecewiseFunction<> membraneStresses;
            assembler.constructStress(mp_def,membraneStresses,stress_type::membrane);
            gsField<> membraneStress(mp_def,membraneStresses, true);

            gsPiecewiseFunction<> flexuralStresses;
            assembler.constructStress(mp_def,flexuralStresses,stress_type::flexural);
            gsField<> flexuralStress(mp_def,flexuralStresses, true);

            gsPiecewiseFunction<> stretches;
            assembler.constructStress(mp_def,stretches,stress_type::principal_stretch);
            gsField<> Stretches(mp_def,stretches, true);

            // gsPiecewiseFunction<> membraneStresses_p;
            // assembler.constructStress(mp_def,membraneStresses_p,stress_type::principal_stress_membrane);
            // gsField<> membraneStress_p(mp_def,membraneStresses_p, true);

            // gsPiecewiseFunction<> flexuralStresses_p;
            // assembler.constructStress(mp_def,flexuralStresses_p,stress_type::principal_stress_flexural);
            // gsField<> flexuralStress_p(mp_def,flexuralStresses_p, true);

            gsPiecewiseFunction<> stretch1;
            assembler.constructStress(mp_def,stretch1,stress_type::principal_stretch_dir1);
            gsField<> stretchDir1(mp_def,stretch1, true);

            gsPiecewiseFunction<> stretch2;
            assembler.constructStress(mp_def,stretch2,stress_type::principal_stretch_dir2);
            gsField<> stretchDir2(mp_def,stretch2, true);

            gsPiecewiseFunction<> stretch3;
            assembler.constructStress(mp_def,stretch3,stress_type::principal_stretch_dir3);
            gsField<> stretchDir3(mp_def,stretch3, true);


            gsField<> solutionField(mp,deformation, true);


            // gsField<> stressField = assembler.constructStress(mp_def,stress_type::membrane_strain);

            std::map<std::string,const gsField<> *> fields;
            fields["Deformation"] = &solutionField;
            fields["Membrane Stress"] = &membraneStress;
            fields["Flexural Stress"] = &flexuralStress;
            fields["Principal Stretch"] = &Stretches;
            // fields["Principal Membrane Stress"] = &membraneStress_p;
            // fields["Principal Flexural Stress"] = &flexuralStress_p;
            fields["Principal Direction 1"] = &stretchDir1;
            fields["Principal Direction 2"] = &stretchDir2;
            fields["Principal Direction 3"] = &stretchDir3;

            gsWriteParaviewMultiPhysics(fields,"stress",5000,true);
        }
    }

    return EXIT_SUCCESS;

}// end main

template <class T>
void evaluateFunction(gsExprEvaluator<T> ev, auto expression, gsVector<T> pt)
{
    gsMatrix<T> evresult = ev.eval( expression,pt );
    gsInfo << "Eval on point ("<<pt.at(0)<<" , "<<pt.at(1)<<") :\n"<< evresult;
    gsInfo << "\nEnd ("<< evresult.rows()<< " x "<<evresult.cols()<<")\n";
};

template <class T>
gsMultiPatch<T> RectangularDomain(int n, int p, T L, T B)
{
  gsMultiPatch<T> mp = RectangularDomain(n, n, p, p, L, B);
  return mp;
}

template <class T>
gsMultiPatch<T> RectangularDomain(int n, int m, int p, int q, T L, T B)
{
  // -------------------------------------------------------------------------
  // --------------------------Make beam geometry-----------------------------
  // -------------------------------------------------------------------------
  int dim = 3; //physical dimension
  gsKnotVector<> kv0;
  kv0.initUniform(0,1,0,p+1,1);
  gsKnotVector<> kv1;
  kv1.initUniform(0,1,0,q+1,1);

  for(index_t i = 0; i< n; ++i)
      kv0.uniformRefine();
  for(index_t i = 0; i< m; ++i)
      kv1.uniformRefine();

  // Make basis
  gsTensorBSplineBasis<2,T> basis(kv0,kv1);

  // Initiate coefficient matrix
  gsMatrix<> coefs(basis.size(),dim);
  // Number of control points needed per component
  size_t len0 = basis.component(0).size();
  size_t len1 = basis.component(1).size();
  gsVector<> coefvec0(len0);
  // Uniformly distribute control points per component
  coefvec0.setLinSpaced(len0,0.0,L);
  gsVector<> coefvec1(basis.component(1).size());
  coefvec1.setLinSpaced(len1,0.0,B);

  // Z coordinate is zero
  coefs.col(2).setZero();

  // Define a matrix with ones
  gsVector<> temp(len0);
  temp.setOnes();
  for (size_t k = 0; k < len1; k++)
  {
    // First column contains x-coordinates (length)
    coefs.col(0).segment(k*len0,len0) = coefvec0;
    // Second column contains y-coordinates (width)
    coefs.col(1).segment(k*len0,len0) = temp*coefvec1.at(k);
  }
  // Create gsGeometry-derived object for the patch
  gsTensorBSpline<2,real_t> shape(basis,coefs);

  gsMultiPatch<T> mp;
  mp.addPatch(shape);
  mp.addAutoBoundaries();

  return mp;
}

template <class T>
gsMultiPatch<T> Rectangle(T L, T B)
{
  // -------------------------------------------------------------------------
  // --------------------------Make beam geometry-----------------------------
  // -------------------------------------------------------------------------
  int dim = 3; //physical dimension
  gsKnotVector<> kv0;
  kv0.initUniform(0,1,0,2,1);
  gsKnotVector<> kv1;
  kv1.initUniform(0,1,0,2,1);

  // Make basis
  gsTensorBSplineBasis<2,T> basis(kv0,kv1);

  // Initiate coefficient matrix
  gsMatrix<> coefs(basis.size(),dim);
  // Number of control points needed per component
  size_t len0 = basis.component(0).size();
  size_t len1 = basis.component(1).size();
  gsVector<> coefvec0(len0);
  // Uniformly distribute control points per component
  coefvec0.setLinSpaced(len0,0.0,L);
  gsVector<> coefvec1(basis.component(1).size());
  coefvec1.setLinSpaced(len1,0.0,B);

  // Z coordinate is zero
  coefs.col(2).setZero();

  // Define a matrix with ones
  gsVector<> temp(len0);
  temp.setOnes();
  for (size_t k = 0; k < len1; k++)
  {
    // First column contains x-coordinates (length)
    coefs.col(0).segment(k*len0,len0) = coefvec0;
    // Second column contains y-coordinates (width)
    coefs.col(1).segment(k*len0,len0) = temp*coefvec1.at(k);
  }
  // Create gsGeometry-derived object for the patch
  gsTensorBSpline<2,real_t> shape(basis,coefs);

  gsMultiPatch<T> mp;
  mp.addPatch(shape);
  mp.addAutoBoundaries();

  return mp;
}


template <class T>
gsMultiPatch<T> RectangularDomainVert(int n, int p, T L, T B)
{
  gsMultiPatch<T> mp = RectangularDomainVert(n, n, p, p, L, B);
  return mp;
}

template <class T>
gsMultiPatch<T> RectangularDomainVert(int n, int m, int p, int q, T L, T B)
{
  // -------------------------------------------------------------------------
  // --------------------------Make beam geometry-----------------------------
  // -------------------------------------------------------------------------
  int dim = 3; //physical dimension
  gsKnotVector<> kv0;
  kv0.initUniform(0,1,0,p+1,1);
  gsKnotVector<> kv1;
  kv1.initUniform(0,1,0,q+1,1);

  for(index_t i = 0; i< n; ++i)
      kv0.uniformRefine();
  for(index_t i = 0; i< m; ++i)
      kv1.uniformRefine();

  // Make basis
  gsTensorBSplineBasis<2,T> basis(kv0,kv1);

  // Initiate coefficient matrix
  gsMatrix<> coefs(basis.size(),dim);
  // Number of control points needed per component
  size_t len0 = basis.component(0).size();
  size_t len1 = basis.component(1).size();
  gsVector<> coefvec0(len0);
  // Uniformly distribute control points per component
  coefvec0.setLinSpaced(len0,0.0,L);
  gsVector<> coefvec1(basis.component(1).size());
  coefvec1.setLinSpaced(len1,0.0,B);

  // Z coordinate is zero
  coefs.col(2).setZero();

  // Define a matrix with ones
  gsVector<> temp(len0);
  temp.setOnes();
  for (index_t k = 0; k < len1; k++)
  {
    // First column contains x-coordinates (length)
    coefs.col(0).segment(k*len0,len0) = coefvec0;
    // Second column contains z-coordinates (height)
    coefs.col(2).segment(k*len0,len0) = temp*coefvec1.at(k);
  }
  // Create gsGeometry-derived object for the patch
  gsTensorBSpline<2,real_t> shape(basis,coefs);

  gsMultiPatch<T> mp;
  mp.addPatch(shape);
  mp.addAutoBoundaries();

  return mp;
}

template <class T>
gsMultiPatch<T> RectangularDomain90(int n, int p, T L, T B)
{
  gsMultiPatch<T> mp = RectangularDomain90(n, n, p, p, L, B);
  return mp;
}

template <class T>
gsMultiPatch<T> RectangularDomain90(int n, int m, int p, int q, T L, T B)
{
  // -------------------------------------------------------------------------
  // --------------------------Make beam geometry-----------------------------
  // -------------------------------------------------------------------------
  int dim = 3; //physical dimension
  gsKnotVector<> kv0;
  kv0.initUniform(0,1,0,p+1,1);
  gsKnotVector<> kv1;
  kv1.initUniform(0,1,0,q+1,1);

  for(index_t i = 0; i< n; ++i)
      kv0.uniformRefine();
  for(index_t i = 0; i< m; ++i)
      kv1.uniformRefine();

  // Make basis
  gsTensorBSplineBasis<2,T> basis(kv0,kv1);

  // Initiate coefficient matrix
  gsMatrix<> coefs(basis.size(),dim);
  // Number of control points needed per component
  size_t len0 = basis.component(0).size();
  size_t len1 = basis.component(1).size();
  gsVector<> coefvec0(len0);
  // Uniformly distribute control points per component
  coefvec0.setLinSpaced(len0,0.0,L);
  gsVector<> coefvec1(basis.component(1).size());
  coefvec1.setLinSpaced(len1,0.0,B);

  // Z coordinate is zero
  coefs.col(2).setZero();

  // Define a matrix with ones
  gsVector<> temp(len0);
  temp.setOnes();
  for (index_t k = 0; k < len1; k++)
  {
    // First column contains x-coordinates (length)
    coefs.col(1).segment(k*len0,len0) = coefvec0;
    // Second column contains z-coordinates (height)
    coefs.col(0).segment(k*len0,len0) = temp*coefvec1.at(k);
  }
  // Create gsGeometry-derived object for the patch
  gsTensorBSpline<2,real_t> shape(basis,coefs);

  gsMultiPatch<T> mp;
  mp.addPatch(shape);
  mp.addAutoBoundaries();

  return mp;
}

template <class T>
gsMultiPatch<T> FrustrumDomain(int n, int p, T R1, T R2, T h)
{
  // -------------------------------------------------------------------------
  // --------------------------Make beam geometry-----------------------------
  // -------------------------------------------------------------------------
  // n = number of uniform refinements over the height; n = 0, only top and bottom part

  int dim = 3; //physical dimension
  gsKnotVector<> kv0;
  kv0.initUniform(0,1,0,3,1);
  gsKnotVector<> kv1;
  kv1.initUniform(0,1,0,3,1);

  // Refine n times
  for(index_t i = 0; i< n; ++i)
      kv1.uniformRefine();

  gsDebug<<kv1;

  // Make basis
  // gsTensorNurbsBasis<2,T> basis(kv0,kv1);

  // Initiate coefficient matrix
  index_t N = math::pow(2,n)+2;
  gsMatrix<> coefs(3*N,dim);
  gsMatrix<> tmp(3,3);
  T R,H;

  gsMatrix<> weights(3*N,1);
  for (index_t k=0; k!= N; k++)
  {
    R = k*(R2-R1)/(N-1) + R1;
    H = k*h/(N-1);
    tmp<< R,0,H,
          R,R,H,
          0,R,H;

    coefs.block(3*k,0,3,3) = tmp;

    weights.block(3*k,0,3,1) << 1,0.70711,1;
  }

  // Create gsGeometry-derived object for the patch
  gsTensorNurbs<2,real_t> shape(kv0,kv1,coefs,weights);

  gsMultiPatch<T> mp;
  mp.addPatch(shape);
  mp.addAutoBoundaries();

  // Elevate up to order p
  if (p>2)
  {
    for(index_t i = 2; i< p; ++i)
        mp.patch(0).degreeElevate();    // Elevate the degree
  }

  // // Refine n times
  // for(index_t i = 0; i< n; ++i)
  //     mp.patch(0).uniformRefine();

  return mp;
}


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
