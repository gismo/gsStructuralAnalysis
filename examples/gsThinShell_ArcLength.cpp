/** @file gsThinShell_ArcLength.cpp

    @brief Code for the arc-length method of a shell based on loads

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#include <gismo.h>

#include <gsKLShell/gsThinShellAssembler.h>
#include <gsKLShell/getMaterialMatrix.h>
// #include <gsThinShell/gsArcLengthIterator.h>
#include <gsStructuralAnalysis/gsArcLengthIterator.h>

using namespace gismo;

template <class T>
gsMultiPatch<T> RectangularDomain(int n, int m, int p, int q, T L, T B, bool clamped = false, T offset = 0.1);
template <class T>
gsMultiPatch<T> RectangularDomain(int n, int p, T L, T B, bool clamped = false, T offset = 0.1);

template <class T>
gsMultiPatch<T> Rectangle(T L, T B);

template <class T>
gsMultiPatch<T> AnnularDomain(int n, int p, T R1, T R2);
template <class T>
gsMultiPatch<T> FrustrumDomain(int n, int p, T R1, T R2, T h);

template <class T>
void addClamping(gsMultiPatch<T> &mp, index_t patch, std::vector<boxSide> sides, T offset);

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
    int numElevate    = 1;
    int numHref       = 1;
    bool plot         = false;
    bool stress       = false;
    bool membrane     = false;
    bool SingularPoint= false;
    bool quasiNewton  = false;
    int quasiNewtonInt= -1;
    bool adaptive     = false;
    int step          = 10;
    int method        = 2; // (0: Load control; 1: Riks' method; 2: Crisfield's method; 3: consistent crisfield method; 4: extended iterations)
    bool symmetry     = false;
    bool deformed     = false;

    real_t thickness  = 1e-3;
    real_t E_modulus  = 1;
    real_t PoissonRatio = 0;
    real_t Density    = 1e0;
    real_t tau        = 1e4;

    index_t Compressibility = 0;
    index_t material  = 0;
    real_t Ratio      = 7.0;
    bool composite = false;
    index_t impl = 1; // 1= analytical, 2= generalized, 3= spectral

    real_t eta        = 0;
    real_t Spring     = 0;

    real_t relax      = 1.0;

    int testCase      = 0;

    int result        = 0;

    bool write        = false;
    bool crosssection = false;

    index_t maxit     = 20;

    // Arc length method options
    real_t dL         = 0; // General arc length
    real_t dLb        = 0.5; // Ard length to find bifurcation
    real_t tol        = 1e-6;
    real_t tolU       = 1e-6;
    real_t tolF       = 1e-3;

    std::string wn("data.csv");

    std::string assemberOptionsFile("options/solver_options.xml");

    gsCmdLine cmd("Arc-length analysis for thin shells.");
    cmd.addString( "f", "file", "Input XML file for assembler options", assemberOptionsFile );

    cmd.addInt("t", "testcase", "Test case: 0: clamped-clamped, 1: pinned-pinned, 2: clamped-free", testCase);

    cmd.addInt("r","hRefine", "Number of dyadic h-refinement (bisection) steps to perform before solving", numHref);
    cmd.addInt("e","degreeElevation", "Number of degree elevation steps to perform on the Geometry's basis before solving", numElevate);
    cmd.addInt( "M", "Material", "Material law",  material );
    cmd.addReal("C", "MaterialRatio", "Material Ratio",  Ratio );
    cmd.addInt( "c", "Compressibility", "1: compressible, 0: incompressible",  Compressibility );
    cmd.addInt( "I", "Implementation", "Implementation: 1= analytical, 2= generalized, 3= spectral",  impl );
    cmd.addSwitch("composite", "Composite material", composite);

    cmd.addReal("T","hdim", "thickness of the plate", thickness);

    cmd.addReal("S","spring", "Nondimensional Spring Stiffness (case 2 and 3 only!)", eta);

    cmd.addInt("m","Method", "Arc length method; 1: Crisfield's method; 2: RIks' method.", method);
    cmd.addReal("L","dLb", "arc length", dLb);
    cmd.addReal("l","dL", "arc length after bifurcation", dL);
    cmd.addReal("A","relaxation", "Relaxation factor for arc length method", relax);

    cmd.addReal("F","factor", "factor for bifurcation perturbation", tau);
    cmd.addInt("q","QuasiNewtonInt","Use the Quasi Newton method every INT iterations",quasiNewtonInt);
    cmd.addInt("N", "maxsteps", "Maximum number of steps", step);

    cmd.addSwitch("adaptive", "Adaptive length ", adaptive);
    cmd.addSwitch("bifurcation", "Compute singular points and bifurcation paths", SingularPoint);
    cmd.addSwitch("quasi", "Use the Quasi Newton method", quasiNewton);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addSwitch("stress", "Plot stress in ParaView format", stress);
    cmd.addSwitch("cross", "Write cross-section to file", crosssection);
    cmd.addSwitch("membrane", "Use membrane model (no bending)", membrane);
    cmd.addSwitch("symmetry", "Use symmetry boundary condition (different per problem)", symmetry);
    cmd.addSwitch("deformed", "plot on deformed shape", deformed);
    cmd.addSwitch("write", "write to file", write);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsFileData<> fd(assemberOptionsFile);
    gsOptionList opts;
    fd.getFirst<gsOptionList>(opts);

    gsMultiPatch<> mp;
    real_t aDim;
    real_t bDim;

    if (dL==0)
    {
      dL = dLb;
    }

    /*
      Case 0: Simply supported beam under compressive load                  --- Validation settings: -L 1eX -l 1eX -M 14 -N 500 -r X -e X
      Case 1: Clamped beam under compressive load                           --- Validation settings: -L 1eX -l 1eX -M 14 -N 500 -r X -e X
    */
    if (testCase==0 || testCase==1)
    {
      E_modulus = 1e8;
      thickness = 0.005313292845913*2;
      PoissonRatio = 0.0;
      aDim = 1.0;
      bDim = 0.1;
      real_t EI = 1.0/12.0*(bDim*math::pow(thickness,3))*E_modulus;
      Spring = math::pow(eta/aDim,4)*EI/bDim;
      mp = RectangularDomain(numHref, 0, numElevate+2, 2, aDim, bDim);
      gsInfo<<"S = "<<Spring<<"; eta = "<<eta<<"\n";
    }
    /*
      Case 2: Clamped beam (left) under vertical end load                   --- Validation settings: -L 5e-1 -M 0 -N 10 -r 2 -e 1 (--plot --write -q 5)
                                                                                Fig 3b from: Pagani, A., & Carrera, E. (2018). Unified formulation of geometrically nonlinear refined beam theories. Mechanics of Advanced Materials and Structures, 25(1), 15–31. https://doi.org/10.1080/15376494.2016.1232458
      Case 3: Clamped beam (left) under horizontal compressive end load     --- Validation settings: -L 5e-5 -l 1e1 -M 0 -N 100 -r 3 -e 1
                                                                                Fig 5  from: Pagani, A., & Carrera, E. (2018). Unified formulation of geometrically nonlinear refined beam theories. Mechanics of Advanced Materials and Structures, 25(1), 15–31. https://doi.org/10.1080/15376494.2016.1232458
    */
    else if (testCase==2 || testCase==3)
    {
      E_modulus = 75e6;
      thickness = 0.01;
      PoissonRatio = 0.0;
      aDim = 1.0;
      bDim = 0.01;
      mp = RectangularDomain(numHref, 0, numElevate+2, 2, aDim, bDim);
    }
    /*
      Case 4: Uniaxial tension of a square plate                            --- Validation settings: -L 1eX -l 1eX -M 14 -N 500 -r X -e X
              (bottom boundary fixed in Y, left boundary fixed in X, right boundary normal load)
      Case 5: Biaxial tension of a square plate                             --- Validation settings: -L 1eX -l 1eX -M 14 -N 500 -r X -e X
              (bottom boundary fixed in Y, left boundary fixed in X, right and top boundaries normal load)
    */
    else if (testCase==4 || testCase==5)
    {
      aDim = 1.0;
      bDim = 1.0;
      real_t mu = 1.5e6;
      thickness = 0.001;
      if (!Compressibility)
        PoissonRatio = 0.5;
      else
        PoissonRatio = 0.45;
      E_modulus = 2*mu*(1+PoissonRatio);
      // PoissonRatio = 0;
      mp = RectangularDomain(numHref, numElevate+2, aDim, bDim);

      gsInfo<<"mu = "<<E_modulus / (2 * (1 + PoissonRatio))<<"\n";
    }
    /*
      Case 6: Constrained tension (see Roohbakhshan2017)
    */
    else if (testCase==6)
    {
      aDim = 10.0e-3;
      bDim = 10.0e-3;
      real_t mu = 10e3;
      thickness = 0.25e-3;
      if ((!Compressibility) && (material!=0))
        PoissonRatio = 0.5;
      else
        PoissonRatio = 0.45;

      if (material==2 || material==12)
        mu = 10e3;
      else if (material==3 || material==13)
        mu = 30e3;

      E_modulus = 2*mu*(1+PoissonRatio);

      Ratio = 0.5;

      mp = RectangularDomain(numHref, numElevate+2, aDim/2., bDim/2.);
    }
    /*
      Case 7: Constrained tension (Chopin2019)
      Chopin, J., Panaitescu, A., & Kudrolli, A. (2018). Corner singularities and shape of stretched elastic sheets. Physical Review E, 98(4), 043003. https://doi.org/10.1103/PhysRevE.98.043003Chopin, J., Panaitescu, A., & Kudrolli, A. (2018). Corner singularities and shape of stretched elastic sheets. Physical Review E, 98(4), 043003. https://doi.org/10.1103/PhysRevE.98.043003
    */
    else if (testCase == 7)
    {
      E_modulus = 1;
      thickness = 0.15;
      if (!Compressibility)
        PoissonRatio = 0.5;
      else
        PoissonRatio = 0.45;

      E_modulus = 1;

      bDim = thickness / 1.9e-3;
      aDim = 2*bDim;

      // Ratio = 2.5442834138486314;
      Ratio = 0.5;

      mp = Rectangle(aDim, bDim);

      for(index_t i = 0; i< numElevate; ++i)
        mp.patch(0).degreeElevate();    // Elevate the degree

      // h-refine
      for(index_t i = 0; i< numHref; ++i)
        mp.patch(0).uniformRefine();
    }
    /*
      Case 8: Balloon subject to increasing internal pressure               --- Validation settings: -L 1eX -l 1eX -M 14 -N 500 -r X -e X
    */
    else if (testCase == 8)
    {
        thickness = 0.1;
        real_t mu = 4.225e5;
        if (!Compressibility)
          PoissonRatio = 0.5;
        else
          PoissonRatio = 0.45;
        E_modulus = 2*mu*(1+PoissonRatio);
        gsReadFile<>("surface/eighth_sphere.xml", mp);

        for(index_t i = 0; i< numElevate; ++i)
          mp.patch(0).degreeElevate();    // Elevate the degree

        // h-refine
        for(index_t i = 0; i< numHref; ++i)
          mp.patch(0).uniformRefine();
    }
    /*
      Case 9: Frustrum with constrained top boundary                          --- Validation settings: -L 1eX -l 1eX -M 14 -N 500 -r X -e X
      Case 10: Frustrum with unconstrained top boundary                        --- Validation settings: -L 1eX -l 1eX -M 14 -N 500 -r X -e X
    */
    else if (testCase == 9 || testCase == 10)
    {
        thickness = 0.1;
        real_t mu = 4.225;
        PoissonRatio = 0.5;
        E_modulus = 2*mu*(1+PoissonRatio);
        // gsReadFile<>("quarter_frustrum.xml", mp);

        // R1 is radius on bottom, R2 is radius on top
        mp = FrustrumDomain(numHref,numElevate+2,2.0,1.0,1.0);
    }
    /*
      Half cylinder compressed from one side (see Kiendl2015)
    */
    else if (testCase == 11)
    {
        thickness = 2e-3;
        PoissonRatio = 0.4;
        E_modulus = 168e9; // GPa
        gsReadFile<>("surface/half_cylinder.xml", mp);
        Ratio = 4;

        for(index_t i = 0; i< numElevate; ++i)
          mp.patch(0).degreeElevate();    // Elevate the degree

        // h-refine
        for(index_t i = 0; i< numHref; ++i)
          mp.patch(0).uniformRefine();
    }
    /*
        Shallow scordelis loo roof
     */
    else if (testCase==12 || testCase==13 || testCase==14)
    {
      // thickness = 0.5*2.286;
      // E_modulus = 3102.75e2;
      // PoissonRatio = 0.3;

      if (testCase==12)
        thickness = 6.35;
      if (testCase==13)
        thickness = 12.7;
      if (testCase==14)
        thickness = 16.75;

      E_modulus = 3102.75;
      PoissonRatio = 0.3;

      gsReadFile<>("surface/scordelis_lo_roof_shallow.xml", mp);

      for(index_t i = 0; i< numElevate; ++i)
        mp.patch(0).degreeElevate();    // Elevate the degree

      // h-refine
      for(index_t i = 0; i< numHref; ++i)
        mp.patch(0).uniformRefine();
    }
    // /*
    //     Case 15: Lifted circular ring with slit
    //  */
    // else if (testCase==15)
    // {
    //   thickness = 3;
    //   E_modulus = 21e6;
    //   PoissonRatio = 0.0;
    //   gsReadFile<>("planar/circle_strip.xml", mp);
    //   mp.embed(3);

    //   for(index_t i = 0; i< numElevate; ++i)
    //     mp.patch(0).degreeElevate();    // Elevate the degree

    //   // h-refine
    //   for(index_t i = 0; i< numHref; ++i)
    //     mp.patch(0).uniformRefine();
    // }

    gsMultiBasis<> dbasis(mp);
    gsInfo<<"Basis (patch 0): "<< mp.patch(0).basis() << "\n";

    // Boundary conditions
    gsBoundaryConditions<> BCs;
    gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();

    // Initiate Surface forces
    std::string tx("0");
    std::string ty("0");
    std::string tz("0");

    gsVector<> tmp(mp.targetDim());
    gsVector<> neu(mp.targetDim());
    tmp.setZero();
    neu.setZero();
    gsConstantFunction<> neuData(neu,mp.targetDim());

    // Unscaled load
    real_t Load = 0;

    std::string output = "solution";
    std::string dirname = "ArcLengthResults";
    real_t pressure = 0.0;
    gsVector<> foundation(3);
    foundation<<0,0,Spring;

    gsMatrix<> writePoints(2,3);
    writePoints.col(0)<< 0.0,0.5;
    writePoints.col(1)<< 0.5,0.5;
    writePoints.col(2)<< 1.0,0.5;
    index_t cross_coordinate = -1;
    real_t cross_val = 0.0;

    if (testCase == 0)
    {
        GISMO_ASSERT(mp.targetDim()==3,"Geometry must be surface (targetDim=3)!");
        // Pinned-Pinned
        tmp << 1e-1, 0, 0;
        neuData.setValue(tmp,3);
        // // Clamped-Clamped
        BCs.addCondition(boundary::west, condition_type::neumann, &neuData ); // unknown 0 - x
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

        // dL =  1e-2;
        // dLb = 1e-4;
        // tol = 1e-3;

        dirname = dirname + "/Beam_pinned-pinned";
        output =  "solution";
        wn = output + "data.txt";
        SingularPoint = true;
    }
    else if (testCase == 1)
    {
        GISMO_ASSERT(mp.targetDim()==3,"Geometry must be surface (targetDim=3)!");
        real_t Area = bDim*thickness;
        real_t EA = E_modulus*Area;
        Load = EA*1e-6;
        tmp << Load, 0, 0;

        // tmp << 0, 0, Load;
        neuData.setValue(tmp,3);
        // // Clamped-Clamped
        BCs.addCondition(boundary::west, condition_type::neumann, &neuData ); // unknown 0 - x
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y

        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

        BCs.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 2 );
        BCs.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 2 );

        // dL =  1e-3;
        // dLb = 1e-3;

        dirname = dirname + "/Beam_clamped-clamped";
        output =  "solution";
        wn = output + "data.txt";
        SingularPoint = true;
    }
    else if (testCase == 2)
    {
        GISMO_ASSERT(mp.targetDim()==3,"Geometry must be surface (targetDim=3)!");
      BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
      BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
      BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

      BCs.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 2 );

      Load = 1e-4;
      gsVector<> point(2);
      gsVector<> load (3);
      point<< 1.0, 0.5 ;
      load << 0.0, 0.0, Load ;
      pLoads.addLoad(point, load, 0 );

      // dL =  1e-3;
      // dLb = 2e0;

      dirname = dirname + "/Beam_clamped-verticalLoad";
      output =  "solution";
      wn = output + "data.txt";
      SingularPoint = false;
    }
    else if (testCase == 3)
    {
      GISMO_ASSERT(mp.targetDim()==3,"Geometry must be surface (targetDim=3)!");
      BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
      BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
      BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

      BCs.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 2 );

      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
      BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y

      Load = 1e-1;

      // dL =  3e-0;
      // dLb = 0.8e-4;

      Load = 1e-1;
      gsVector<> point(2);
      gsVector<> load (3);
      point<< 1.0, 0.5 ;
      load << -Load, 0.0, 0.0 ;
      pLoads.addLoad(point, load, 0 );

      dirname = dirname + "/Beam_clamped-horizontalLoad";
      output =  "solution";
      wn = output + "data.txt";
      SingularPoint = true;
    }
    else if (testCase == 4) // Uniaxial tension; use with hyperelastic material model!
    {
      GISMO_ASSERT(mp.targetDim()==3,"Geometry must be surface (targetDim=3)!");
      BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
      BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

      BCs.addCondition(boundary::east, condition_type::collapsed, 0, 0, false, 0 ); // unknown 1 - y
      BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z.


      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 1 - y
      BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 1 - y

      Load = 1e0;
      gsVector<> point(2);
      gsVector<> load (3);
      point<< 1.0, 0.5 ;
      load << Load,0.0, 0.0;
      pLoads.addLoad(point, load, 0 );

      dirname = dirname + "/UniaxialTension";
      output =  "solution";
      wn = output + "data.txt";
      SingularPoint = true;
    }
    else if (testCase == 5) // Bi-axial tension; use with hyperelastic material model!
    {
      GISMO_ASSERT(mp.targetDim()==3,"Geometry must be surface (targetDim=3)!");
      BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
      BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

      BCs.addCondition(boundary::east, condition_type::collapsed, 0, 0, false, 0 ); // unknown 1 - y
      BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z.

      BCs.addCondition(boundary::north, condition_type::collapsed, 0, 0, false, 1 ); // unknown 1 - y
      BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z.

      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 1 - y

      real_t load_factor = 1;
      Load = 1e0;
      gsVector<> point(2);
      gsVector<> load (3);
      point<< 1.0, 0.5 ;
      load << Load,0.0, 0.0;
      pLoads.addLoad(point, load, 0 );

      point<< 0.5, 1.0 ;
      load << 0.0, Load/load_factor, 0.0;
      pLoads.addLoad(point, load, 0 );

      dirname = dirname + "/BiaxialTension";
      output =  "solution";
      wn = output + "data.txt";
      SingularPoint = true;
    }
    // else if (testCase == 6) //???
    // {
    //   GISMO_ASSERT(mp.targetDim()==3,"Geometry must be surface (targetDim=3)!");
    //   BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 1 - x
    //   BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - x
    //   BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 1 - x

    //   BCs.addCondition(boundary::east,  condition_type::dirichlet, 0, 0, false, 0 ); // unknown 1 - x
    //   BCs.addCondition(boundary::east,  condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - x
    //   BCs.addCondition(boundary::east,  condition_type::dirichlet, 0, 0, false, 2 ); // unknown 1 - x

    //   BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - x
    //   BCs.addCondition(boundary::south, condition_type::clamped,   0, 0, false, 2 ); // unknown 1 - x

    //   BCs.addCondition(boundary::west,  condition_type::dirichlet, 0, 0, false, 0 ); // unknown 1 - x
    //   BCs.addCondition(boundary::west,  condition_type::clamped,   0, 0, false, 2 ); // unknown 1 - x

    //   pressure = 1.0;

    //   dirname = dirname + "/" + "Case" + std::to_string(testCase);
    //   output =  "solution";
    //   wn = output + "data.txt";
    //   SingularPoint = false;

    //   writePoints.resize(2,3);
    //   writePoints.col(0)<< 0.0,0.5;
    //   writePoints.col(1)<< 0.5,0.5;
    //   writePoints.col(2)<< 1.0,0.5;
    // }
    else if (testCase == 6 || testCase == 7) // Uniaxial tension with fixed ends
    {
      GISMO_ASSERT(mp.targetDim()==3,"Geometry must be surface (targetDim=3)!");
      for (index_t i=0; i!=3; ++i)
      {
         BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, i ); // unknown 2 - z
      }
      BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

      BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0 ,false,1);
      BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0 ,false,2);
      BCs.addCondition(boundary::east, condition_type::collapsed, 0, 0 ,false,0);

      gsVector<> point(2); point<< 1.0, 0.5 ;
      gsVector<> load (3); load << 0.1, 0.0, 0.0 ;
      pLoads.addLoad(point, load, 0 );

      real_t alpha, beta;
      alpha = bDim/thickness;
      beta = aDim/bDim;
      dirname = dirname + "/" + "Tension_-r" + std::to_string(numHref) + "-e" + std::to_string(numElevate) + "-M" + std::to_string(material) + "-c" + std::to_string(Compressibility) + "-alpha" + std::to_string(alpha) + "-beta" + std::to_string(beta);
      output =  "solution";
      wn = output + "data.txt";

      cross_coordinate = 1;
      cross_val = 1.0;
    }
    else if (testCase == 8)
    {
      GISMO_ASSERT(mp.targetDim()==3,"Geometry must be surface (targetDim=3)!");
      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 2 - z
      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 2 - z

      BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z


      // Symmetry in x-direction:
      BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 );
      BCs.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 1 );
      BCs.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 2 );

      // Symmetry in y-direction:
      BCs.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 0 );
      BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 );
      BCs.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 2 );

      // Pressure
      pressure = 1e3;
      maxit = 50;

      dirname = dirname + "/" + "Balloon";
      output =  "solution";
      wn = output + "data.txt";
    }
    else if (testCase == 9)
    {
      GISMO_ASSERT(mp.targetDim()==3,"Geometry must be surface (targetDim=3)!");
      Load = -1;
      neu << 0, 0, Load;
      neuData.setValue(neu,3);

      BCs.addCondition(boundary::north, condition_type::neumann, &neuData );
      BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 2 - z
      BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 2 - z
      BCs.addCondition(boundary::north, condition_type::collapsed, 0, 0, false, 2 ); // unknown 1 - y

      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

      // Symmetry in x-direction:
      BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 );
      BCs.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 1 );
      BCs.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 2 );

      // Symmetry in y-direction:
      BCs.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 0 );
      BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 );
      BCs.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 2 );

      dirname = dirname + "/" + "Frustrum_-r=" + std::to_string(numHref) + "-e" + std::to_string(numElevate) + "-M" + std::to_string(material) + "_solution";
      output =  "solution";
      wn = output + "data.txt";

      writePoints.resize(2,3);
      writePoints.col(0)<<0.0,1.0;
      writePoints.col(1)<<0.5,1.0;
      writePoints.col(2)<<1.0,1.0;
    }
    else if (testCase == 10)
    {
      GISMO_ASSERT(mp.targetDim()==3,"Geometry must be surface (targetDim=3)!");
      Load = -1;
      neu << 0, 0, Load;
      neuData.setValue(neu,3);

      BCs.addCondition(boundary::north, condition_type::neumann, &neuData );
      // BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 2 - z
      // BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 2 - z
      BCs.addCondition(boundary::north, condition_type::collapsed, 0, 0, false, 2 ); // unknown 1 - y

      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

      // Symmetry in x-direction:
      BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 );
      BCs.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 1 );
      BCs.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 2 );

      // Symmetry in y-direction:
      BCs.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 0 );
      BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 );
      BCs.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 2 );

      dirname = dirname + "/" + "Frustrum2_-r=" + std::to_string(numHref) + "-e" + std::to_string(numElevate) + "-M" + std::to_string(material) + "_solution";
      output = "solution";
      wn = output + "data.txt";

      writePoints.resize(2,3);
      writePoints.col(0)<<0.0,1.0;
      writePoints.col(1)<<0.5,1.0;
      writePoints.col(2)<<1.0,1.0;
    }
    else if (testCase == 11)
    {
      GISMO_ASSERT(mp.targetDim()==3,"Geometry must be surface (targetDim=3)!");
      Load = -1;
      neu << 0, 0, Load;
      neuData.setValue(neu,3);

      BCs.addCondition(boundary::north, condition_type::neumann, &neuData );
      BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 );
      BCs.addCondition(boundary::north, condition_type::clamped, 0, 0, false, 2 );

      BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 );
      BCs.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 2 );

      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 0 );
      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 );
      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 );
      BCs.addCondition(boundary::south, condition_type::clamped, 0, 0, false, 2 );

      dirname = dirname + "/" + "Cylinder-r=" + std::to_string(numHref) + "-e" + std::to_string(numElevate) + "-M" + std::to_string(material) + "_solution";
      output =  "solution";
      wn = output + "data.txt";
      SingularPoint = false;

      writePoints.resize(2,3);
      writePoints.col(0)<<0.0,1.0;
      writePoints.col(1)<<0.5,1.0;
      writePoints.col(2)<<1.0,1.0;

      cross_coordinate = 0; // Constant on x-axis
      cross_val = 1.0; // parametric value x=1.0; this corresponds with the symmetry edge
    }
    else if (testCase == 12 || testCase == 13 || testCase == 14)
    {
      GISMO_ASSERT(mp.targetDim()==3,"Geometry must be surface (targetDim=3)!");
      // Diaphragm conditions
      BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
      BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
      BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
      // BCs.addCornerValue(boundary::southwest, 0.0, 0, 0); // (corner,value, patch, unknown)
      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

      Load = -1e1;
      // Point loads
      gsVector<> point(2);
      gsVector<> load (3);
      point<< 0.5, 0.5 ;
      load << 0.0, 0.0, Load ;
      pLoads.addLoad(point, load, 0 );

      dirname = dirname + "/" +  "Roof_t="+ std::to_string(thickness) + "-r=" + std::to_string(numHref) + "-e" + std::to_string(numElevate) +"_solution";
      output =  "solution";
      wn = output + "data.txt";
      SingularPoint = false;
    }
    // Needs a C1 description of a circle!
    // else if (testCase == 15)
    // {
    //   Load = 0.8;
    //   neu << 0, 0,Load;
    //   neuData.setValue(neu,3);
    //   BCs.addCondition(boundary::east, condition_type::neumann, &neuData );

    //   BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
    //   BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
    //   BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
    //   BCs.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 2 ); // unknown 2 - z

    //   dirname = dirname + "/Circle";
    //   output =  "solution";
    //   wn = output + "data.txt";
    //   SingularPoint = false;
    // }

    std::string commands = "mkdir -p " + dirname;
    const char *command = commands.c_str();
    system(command);

    // plot geometry
    if (plot)
      gsWriteParaview(mp,dirname + "/" + "mp",1000,true);

    if (write)
      initStepOutput(dirname + "/" + wn, writePoints);
    if (crosssection && cross_coordinate!=-1)
    {
      initSectionOutput(dirname,false); // write pointdataX.txt, pointdataY.txt, pointdataZ.txt
      initSectionOutput(dirname,true); // write pointdataX0.txt, pointdataY0.txt, pointdataZ0.txt
      writeSectionOutput(mp,dirname,cross_coordinate,cross_val,201,true);
    }
    else if (crosssection && cross_coordinate==-1)
    {
      gsInfo<<"No cross section can be exported if no coordinate is given...\n";
      crosssection=false;
    }


    gsConstantFunction<> pressFun(pressure,3);
    gsConstantFunction<> foundFun(foundation,3);
    // Initialise solution object
    gsMultiPatch<> mp_def = mp;

    // Linear isotropic material model
    gsConstantFunction<> force(tmp,3);
    gsFunctionExpr<> t(std::to_string(thickness), 3);
    gsFunctionExpr<> E(std::to_string(E_modulus),3);
    gsFunctionExpr<> nu(std::to_string(PoissonRatio),3);
    gsFunctionExpr<> rho(std::to_string(Density),3);
    gsConstantFunction<> ratio(Ratio,3);

    real_t mu = E_modulus / (2 * (1 + PoissonRatio));
    gsConstantFunction<> alpha1(1.3,3);
    gsConstantFunction<> mu1(6.3e5/4.225e5*mu,3);
    gsConstantFunction<> alpha2(5.0,3);
    gsConstantFunction<> mu2(0.012e5/4.225e5*mu,3);
    gsConstantFunction<> alpha3(-2.0,3);
    gsConstantFunction<> mu3(-0.1e5/4.225e5*mu,3);

    index_t kmax = 1;

    std::vector<gsFunctionSet<> * > Gs(kmax);
    std::vector<gsFunctionSet<> * > Ts(kmax);
    std::vector<gsFunctionSet<> * > Phis(kmax);

    gsMatrix<> Gmat = gsCompositeMatrix(E_modulus,E_modulus,0.5 * E_modulus / (1+PoissonRatio),PoissonRatio,PoissonRatio);
    Gmat.resize(Gmat.rows()*Gmat.cols(),1);
    gsConstantFunction<> Gfun(Gmat,3);
    Gs[0] = &Gfun;

    gsConstantFunction<> phi;
    phi.setValue(0,3);

    Phis[0] = &phi;

    gsConstantFunction<> thicks(thickness/kmax,3);
    Ts[0] = &thicks;

    std::vector<gsFunction<>*> parameters;
    if (material==0) // SvK & Composites
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
      parameters.resize(8);
      parameters[0] = &E;
      parameters[1] = &nu;
      parameters[2] = &mu1;
      parameters[3] = &alpha1;
      parameters[4] = &mu2;
      parameters[5] = &alpha2;
      parameters[6] = &mu3;
      parameters[7] = &alpha3;
    }

    gsMaterialMatrixBase<real_t>* materialMatrix;

    gsOptionList options;
    if      (material==0 && impl==1)
    {
        if (composite)
        {
            materialMatrix = new gsMaterialMatrixComposite<3,real_t>(mp,Ts,Gs,Phis);
        }
        else
        {
            options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
            options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
            materialMatrix = getMaterialMatrix<3,real_t>(mp,t,parameters,rho,options);
        }
    }
    else
    {
        options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",material);
        options.addSwitch("Compressibility","Compressibility: (false): Imcompressible | (true): Compressible",Compressibility);
        options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",impl);
        materialMatrix = getMaterialMatrix<3,real_t>(mp,t,parameters,rho,options);
    }

    gsThinShellAssemblerBase<real_t>* assembler;
    if(membrane)
        assembler = new gsThinShellAssembler<3, real_t, false>(mp,dbasis,BCs,force,materialMatrix);
    else
        assembler = new gsThinShellAssembler<3, real_t, true >(mp,dbasis,BCs,force,materialMatrix);


    // Construct assembler object
    assembler->setOptions(opts);
    assembler->setPointLoads(pLoads);
    if (pressure!= 0.0)
        assembler->setPressure(pressFun);
    if (Spring!= 0.0)
        assembler->setFoundation(foundFun);

    gsStopwatch stopwatch;
    real_t time = 0.0;

    typedef std::function<gsSparseMatrix<real_t> (gsVector<real_t> const &)>                                Jacobian_t;
    typedef std::function<gsVector<real_t> (gsVector<real_t> const &, real_t, gsVector<real_t> const &) >   ALResidual_t;
    // Function for the Jacobian
    Jacobian_t Jacobian = [&time,&stopwatch,&assembler,&mp_def](gsVector<real_t> const &x)
    {
      stopwatch.restart();
      assembler->constructSolution(x,mp_def);
      assembler->assembleMatrix(mp_def);
      time += stopwatch.stop();

      gsSparseMatrix<real_t> m = assembler->matrix();
      return m;
    };
    // Function for the Residual
    ALResidual_t ALResidual = [&time,&stopwatch,&assembler,&mp_def](gsVector<real_t> const &x, real_t lam, gsVector<real_t> const &force)
    {
      stopwatch.restart();
      assembler->constructSolution(x,mp_def);
      assembler->assembleVector(mp_def);
      gsVector<real_t> Fint = -(assembler->rhs() - force);
      gsVector<real_t> result = Fint - lam * force;
      time += stopwatch.stop();
      return result; // - lam * force;
    };
    // Assemble linear system to obtain the force vector
    assembler->assemble();
    gsVector<> Force = assembler->rhs();

    gsArcLengthIterator<real_t> arcLength(Jacobian, ALResidual, Force);

    if (!membrane)
    {
      arcLength.options().setInt("Solver",0); // LDLT solver
      arcLength.options().setInt("BifurcationMethod",0); // 0: determinant, 1: eigenvalue
    }
    else
    {
      arcLength.options().setInt("Solver",1); // CG solver
      arcLength.options().setInt("BifurcationMethod",1); // 0: determinant, 1: eigenvalue
    }

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


    gsInfo<<arcLength.options();
    arcLength.applyOptions();
    arcLength.initialize();


    gsParaviewCollection collection(dirname + "/" + output);
    gsParaviewCollection Smembrane(dirname + "/" + "membrane");
    gsParaviewCollection Sflexural(dirname + "/" + "flexural");
    gsParaviewCollection Smembrane_p(dirname + "/" + "membrane_p");
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
      // assembler->constructSolution(solVector,solution);
      arcLength.step();

      // gsInfo<<"m_U = "<<arcLength.solutionU()<<"\n";
      if (!(arcLength.converged()))
      {
        gsInfo<<"Error: Loop terminated, arc length method did not converge.\n";
        dLb = dLb / 2.;
        arcLength.setLength(dLb);
        arcLength.setSolution(Uold,Lold);
        bisected = true;
        k -= 1;
        continue;
        // if (plot)
        // {
        //   solVector = arcLength.solutionU();
        //   Uold = solVector;
        //   Lold = arcLength.solutionL();
        //   assembler->constructSolution(solVector,mp_def);

        //   deformation = mp_def;
        //   deformation.patch(0).coefs() -= mp.patch(0).coefs();// assuming 1 patch here

        //   gsField<> solField(mp,deformation);
        //   std::string fileName = dirname + "/" + output + util::to_string(k);
        //   gsWriteParaview<>(solField, fileName, 5000);
        //   fileName = output + util::to_string(k) + "0";
        //   collection.addTimestep(fileName,k,".vts");
        // }
        // break;
      }

      if (SingularPoint)
      {
        arcLength.computeStability(arcLength.solutionU(),quasiNewton);
        if (arcLength.stabilityChange())
        {
          gsInfo<<"Bifurcation spotted!"<<"\n";
          arcLength.computeSingularPoint(1e-4, 5, Uold, Lold, 1e-10, 0, false);
          arcLength.switchBranch();
          dLb0 = dLb = dL;
          arcLength.setLength(dLb);
        }
      }
      indicator = arcLength.indicator();

      solVector = arcLength.solutionU();
      Uold = solVector;
      Lold = arcLength.solutionL();
      assembler->constructSolution(solVector,mp_def);

      if (testCase==4 || testCase==8 || testCase==9)
      {
        gsMatrix<> pts(2,1);
        pts<<0.5,0.5;
        if (testCase==8 || testCase==9)
        {
          pts.resize(2,3);
          pts.col(0)<<0.0,1.0;
          pts.col(1)<<0.5,1.0;
          pts.col(2)<<1.0,1.0;
        }

        gsMatrix<> lambdas = assembler->computePrincipalStretches(pts,mp_def,0);
        std::streamsize ss = std::cout.precision();
        std::cout <<std::setprecision(20)
                  <<"lambdas = \n"<<lambdas<<"\n";
        std::cout<<std::setprecision(ss);

        real_t S = Lold / 1e-3 / lambdas(0) / lambdas(2);
        real_t San = mu * (math::pow(lambdas(1),2)-1/lambdas(1));
        gsDebugVar(S);
        gsDebugVar(San);
        gsDebugVar(abs(S-San));
      }

      deformation = mp_def;
      deformation.patch(0).coefs() -= mp.patch(0).coefs();// assuming 1 patch here

      // gsDebugVar(mp_def.patch(0).coefs());

      if (testCase==8)
      {
        std::streamsize ss = std::cout.precision();
        std::cout<<std::setprecision(20)
              <<"Pressures:\n"<<pressure*arcLength.solutionL()<<"\n"
                              <<pressure*arcLength.solutionL() * assembler->getArea(mp) / assembler->getArea(mp_def)<<"\n";
        std::cout<<std::setprecision(ss);
      }

      gsInfo<<"Total ellapsed assembly time: "<<time<<" s\n";

      if (plot)
      {
        gsField<> solField;
        if (deformed)
          solField= gsField<>(mp_def,deformation);
        else
          solField= gsField<>(mp,deformation);

        std::string fileName = dirname + "/" + output + util::to_string(k);
        gsWriteParaview<>(solField, fileName, 1000,true);
        fileName = output + util::to_string(k) + "0";
        collection.addTimestep(fileName,k,".vts");
        collection.addTimestep(fileName,k,"_mesh.vtp");
      }
      if (stress)
      {
        std::string fileName;

        gsField<> membraneStress, flexuralStress, membraneStress_p;

        gsPiecewiseFunction<> membraneStresses;
        assembler->constructStress(mp_def,membraneStresses,stress_type::membrane);
        if (deformed)
          membraneStress = gsField<>(mp_def,membraneStresses,true);
        else
          membraneStress = gsField<>(mp,membraneStresses,true);

        fileName = dirname + "/" + "membrane" + util::to_string(k);
        gsWriteParaview( membraneStress, fileName, 1000);
        fileName = "membrane" + util::to_string(k) + "0";
        Smembrane.addTimestep(fileName,k,".vts");

        gsPiecewiseFunction<> flexuralStresses;
        assembler->constructStress(mp_def,flexuralStresses,stress_type::flexural);
        if (deformed)
          flexuralStress = gsField<>(mp_def,flexuralStresses, true);
        else
          flexuralStress = gsField<>(mp,flexuralStresses, true);

        fileName = dirname + "/" + "flexural" + util::to_string(k);
        gsWriteParaview( flexuralStress, fileName, 1000);
        fileName = "flexural" + util::to_string(k) + "0";
        Sflexural.addTimestep(fileName,k,".vts");

        if (impl==3)
        {
          gsPiecewiseFunction<> membraneStresses_p;
          assembler->constructStress(mp_def,membraneStresses_p,stress_type::principal_stress_membrane);
          if (deformed)
            membraneStress_p = gsField<>(mp_def,membraneStresses_p, true);
          else
            membraneStress_p = gsField<>(mp,membraneStresses_p, true);

          fileName = dirname + "/" + "membrane_p" + util::to_string(k);
          gsWriteParaview( membraneStress_p, fileName, 1000);
          fileName = "membrane_p" + util::to_string(k) + "0";
          Smembrane_p.addTimestep(fileName,k,".vts");
        }

      }

      if (write)
        writeStepOutput(arcLength,deformation, dirname + "/" + wn, writePoints,1, 201);

      if (crosssection && cross_coordinate!=-1)
        writeSectionOutput(deformation,dirname,cross_coordinate,cross_val,201,false);

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
      Smembrane.save();
      Sflexural.save();
      Smembrane_p.save();
    }

  return result;
}

template <class T>
gsMultiPatch<T> RectangularDomain(int n, int p, T L, T B, bool clamped, T clampoffset)
{
  gsMultiPatch<T> mp = RectangularDomain(n, n, p, p, L, B, clamped, clampoffset);
  return mp;
}

template <class T>
gsMultiPatch<T> RectangularDomain(int n, int m, int p, int q, T L, T B, bool clamped, T clampoffset)
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

  if (clamped)
  {
    T knotval;
    knotval = kv0.uValue(1);
    kv0.insert(std::min(clampoffset,knotval/2.));

    knotval = kv0.uValue(kv0.uSize()-2);
    kv0.insert(std::max(1-clampoffset,knotval/2.));
  }

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
void addClamping(gsMultiPatch<T>& mp, index_t patch, std::vector<boxSide> sides, T offset) //, std::vector<boxSide> sides, T offset)
{

    gsTensorBSpline<2,T> *geo = dynamic_cast< gsTensorBSpline<2,real_t> * > (&mp.patch(patch));

    T dknot0 = geo->basis().component(0).knots().minIntervalLength();
    T dknot1 = geo->basis().component(1).knots().minIntervalLength();

    gsInfo<<"sides.size() = "<<sides.size()<<"\n";

    index_t k =0;


    for (std::vector<boxSide>::iterator it = sides.begin(); it != sides.end(); it++)
    {
        gsInfo<<"side = "<<(*it)<<"\n";

      if (*it==boundary::west || *it==boundary::east) // west or east
      {
        if (*it==boundary::east) // east, val = 1
          geo->insertKnot(1 - std::min(offset, dknot0 / 2),0);
        else if (*it==boundary::west) // west
          geo->insertKnot(std::min(offset, dknot0 / 2),0);
      }
      else if (*it==boundary::south || *it==boundary::north) // west or east
      {
       if (*it==boundary::north) // north
         geo->insertKnot(1 - std::min(offset, dknot0 / 2),1);
       else if (*it==boundary::south) // south
         geo->insertKnot(std::min(offset, dknot0 / 2),1);
      }
      else if (*it==boundary::none)
        gsWarn<<*it<<"\n";
      else
        GISMO_ERROR("Side unknown, side = " <<*it);

        k++;
gsInfo<<"k = "<<k<<"\n";
    }
}

template <class T>
gsMultiPatch<T> Rectangle(T L, T B) //, int n, int m, std::vector<boxSide> sides, T offset)
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
gsMultiPatch<T> AnnularDomain(int n, int p, T R1, T R2)
{
  // -------------------------------------------------------------------------
  // --------------------------Make beam geometry-----------------------------
  // -------------------------------------------------------------------------
  int dim = 3; //physical dimension
  gsKnotVector<> kv0;
  kv0.initUniform(0,1,0,3,1);
  gsKnotVector<> kv1;
  kv1.initUniform(0,1,0,3,1);

  // Make basis
  // gsTensorNurbsBasis<2,T> basis(kv0,kv1);

  // Initiate coefficient matrix
  gsMatrix<> coefs(9,dim);

  coefs<<R1,0,0,
  (R1+R2)/2,0,0,
  R2,0,0,
  R1,R1,0,
  (R1+R2)/2,(R1+R2)/2,0,
  R2,R2,0,
  0,R1,0,
  0,(R1+R2)/2,0,
  0,R2,0;

  gsMatrix<> weights(9,1);
  weights<<1,1,1,
  0.707106781186548,0.707106781186548,0.707106781186548,
  1,1,1;

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

  // Refine n times
  for(index_t i = 0; i< n; ++i)
      mp.patch(0).uniformRefine();

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
