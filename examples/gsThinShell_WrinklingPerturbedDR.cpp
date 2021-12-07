/** @file gsThinShell_WrinklingPerturbed.cpp

    @brief Performs wrinkling simulations of different cases USING A PERTURBATION from a multipatch

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#include <gismo.h>

#include <gsKLShell/gsThinShellAssembler.h>
#include <gsKLShell/getMaterialMatrix.h>

#include <gsStructuralAnalysis/gsDynamicRelaxationLC.h>
#include <gsStructuralAnalysis/gsStaticDR.h>
#include <gsStructuralAnalysis/gsControlDisplacement.h>

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
void writeStepOutput(const gsMatrix<T> & solution, const T &load, const gsMultiPatch<T> & deformation, const std::string name, const gsMatrix<T> & points, const index_t extreme=-1, const index_t kmax=100);

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
    bool energyPlot = false;
    bool stress       = false;
    bool membrane       = false;
    bool mesh = false;
    int step = 10;
    bool symmetry = false;
    bool deformed = false;
    real_t perturbation = 0;
    int verbose = 0;

    real_t thickness = 1e-3;
    real_t E_modulus     = 1;
    real_t PoissonRatio = 0;
    real_t Density = 1e0;
    gsMultiPatch<> mp, mpBspline;

    index_t Compressibility = 0;
    index_t material = 0;
    real_t Ratio = 7.0;
    bool composite = false;
    index_t impl = 1; // 1= analytical, 2= generalized, 3= spectral

    real_t aDim = 2.5;
    real_t bDim = 1.0;

    int testCase = 2;

    int result = 0;

    bool write = false;
    bool writeG = false;
    bool crosssection = false;

    bool THB = false;

    bool weak = false;

    index_t maxit = 2e6;

    real_t alpha = 2.0;
    real_t damping = 1.0;

    // Arc length method options
    real_t tolE = 1e-3;
    real_t tolF = 1e-5;

    std::string wn("data.csv");

    std::string assemberOptionsFile("options/solver_options.xml");

    std::string fn;

    gsCmdLine cmd("Wrinkling analysis with thin shells using XML perturbation.");
    cmd.addString( "f", "file", "Input XML file for assembler options", assemberOptionsFile );

    cmd.addInt("t", "testcase", "Test case: 0: clamped-clamped, 1: pinned-pinned, 2: clamped-free", testCase);

    cmd.addInt("r","hRefine", "Number of dyadic h-refinement (bisection) steps to perform before solving", numHref);
    cmd.addInt("e","degreeElevation", "Number of degree elevation steps to perform on the Geometry's basis before solving", numElevate);
    cmd.addInt("R","hRefine2", "Number of dyadic h-refinement (bisection) steps to perform before solving (secondary direction)", numHrefL);
    cmd.addInt("E","degreeElevation2", "Number of degree elevation steps to perform on the Geometry's basis before solving (secondary direction)", numElevateL);

    cmd.addInt( "M", "Material", "Material law",  material );
    cmd.addInt( "c", "Compressibility", "1: compressible, 0: incompressible",  Compressibility );
    cmd.addInt( "I", "Implementation", "Implementation: 1= analytical, 2= generalized, 3= spectral",  impl );
    cmd.addSwitch("composite", "Composite material", composite);

    cmd.addReal("T","hdim", "thickness of the plate", thickness);

    cmd.addReal( "a", "alpha", "alpha",  alpha );
    cmd.addReal( "D", "damping", "damping",  damping );
    cmd.addReal("P","perturbation", "perturbation factor", perturbation);

    cmd.addInt("N", "maxsteps", "Maximum number of steps", step);
    cmd.addInt("v","verbose", "0: no; 1: iteration output; 2: Full matrix and vector output", verbose);

    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addSwitch("mesh", "Plot mesh?", mesh);
    cmd.addSwitch("eplot", "Plot the energies", energyPlot);
    cmd.addSwitch("stress", "Plot stress in ParaView format", stress);
    cmd.addSwitch("write", "Write output to file", write);
    cmd.addSwitch("writeG", "Write refined geometry", writeG);
    cmd.addSwitch("cross", "Write cross-section to file", crosssection);
    cmd.addSwitch("membrane", "Use membrane model (no bending)", membrane);
    cmd.addSwitch("symmetry", "Use symmetry boundary condition (different per problem)", symmetry);
    cmd.addSwitch("deformed", "plot on deformed shape", deformed);
    cmd.addSwitch("weak", "Use weak clamping", weak);

    cmd.addSwitch("THB", "Use refinement", THB);

    cmd.addString("i","input", "Perturbation filename", fn);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsFileData<> fd(assemberOptionsFile);
    gsOptionList opts;
    fd.getFirst<gsOptionList>(opts);

    if (numHrefL==-1)
      numHrefL = numHref;
    if (numElevateL==-1)
      numElevateL = numElevate;

    if ((!Compressibility) && (material!=0))
      PoissonRatio = 0.5;
    else
      PoissonRatio = 0.499;

    real_t mu, C01,C10;

    // ![Material data]

    /*
        WRINKLING
        Fu  & PANAITESCU
        2   & 3           Full sheet
        4   & 5           Half sheet
        6   & 7           Quarter sheet
    */
    if (testCase==2 || testCase==3 || testCase==4 || testCase==5 || testCase==6 || testCase==7)
    {
      if (material==3||material==13||material==23)
      {
        if      (testCase==2 || testCase==4 || testCase==6)
        {
          C10 = (0.5-1/22.)*1e6;      // c1/2
          C01 = (1/22.)*1e6;          // c2/2
        }
        else if (testCase==3 || testCase==5 || testCase==7)
        {
          C10 = 6.21485502e4; // c1/2
          C01 = 15.8114570e4; // c2/2
        }
        Ratio = C10/C01;
        mu = 2*(C01+C10);
      }
      else
      {
        if      (testCase==2 || testCase==4 || testCase==6)
          C10 = (0.5)*1e6;
        else if (testCase==3 || testCase==5 || testCase==7)
          C10 = 19.1010178e4;

        mu = 2*C10;
      }
      E_modulus = 2*mu*(1+PoissonRatio);
      gsDebug<<"E = "<<E_modulus<<"; nu = "<<PoissonRatio<<"; mu = "<<mu<<"; ratio = "<<Ratio<<"\n";

      aDim = 0.28;
      bDim = 0.14;
      thickness = 0.14e-3;
    }
    // ![Material data]

    // ![Read Geometry files]
    if (fn.empty())
    {
      std::vector<boxSide> sides;
      if (testCase >= 2 && testCase<=7)
      {
      	sides.push_back(boundary::west);
      	sides.push_back(boundary::east);
      	if (symmetry && (testCase==4 || testCase==5 || testCase==6 || testCase==7))
          sides.push_back(boundary::south);
      }

      if        (testCase==2 || testCase==3)
        mpBspline = Rectangle(aDim,    bDim   );
      else if   (testCase==4 || testCase==5)
        mpBspline = Rectangle(aDim   , bDim/2.);
      else if   (testCase==6 || testCase==7)
        mpBspline = Rectangle(aDim/2., bDim/2.);

      for(index_t i = 0; i< numElevate; ++i)
        mpBspline.patch(0).degreeElevate();    // Elevate the degree

      // h-refine
      for(index_t i = 0; i< numHref; ++i)
        mpBspline.patch(0).uniformRefine();

      addClamping(mpBspline,0,sides, 1e-2);

      index_t N = mpBspline.patch(0).coefs().rows();
      mpBspline.patch(0).coefs().col(2) = gsMatrix<>::Random(N,1);

      // if (testCase==2)
      //   fn = "deformations/wrinklingFu_full.xml";
      // else if (testCase==3)
      //   fn = "deformations/wrinklingPanaitescu_full.xml";
      // else if (testCase==4)
      //   fn = "deformations/wrinklingFu_half.xml";
      // else if (testCase==5)
      //   fn = "deformations/wrinklingPanaitescu_half.xml";
      // else if (testCase==6)
      //   fn = "deformations/wrinklingFu_quarter.xml";
      // else if (testCase==7)
      //   fn = "deformations/wrinklingPanaitescu_quarter.xml";
      // else
      //   GISMO_ERROR("No filename provided..");
    }
    else
      gsReadFile<>(fn,mpBspline);

    mpBspline.patch(0).coefs().col(2) *= perturbation;

      // std::vector<boxSide> sides;
      // sides.push_back(boundary::west);
      // sides.push_back(boundary::east);
      // if (symmetry && (testCase==4 || testCase==5 || testCase==6 || testCase==7))
      //   sides.push_back(boundary::south);

      // if        (testCase==2 || testCase==3)
      //   mpBspline = Rectangle(aDim/2., bDim/2.);
      // else if   (testCase==4 || testCase==5)
      //   mpBspline = Rectangle(aDim   , bDim/2.);
      // else if   (testCase==6 || testCase==7)
      //   mpBspline = Rectangle(aDim   , bDim   );

      // addClamping(mpBspline,0,sides, 1e-2);

    // ![Read Geometry files]

    // Cast all patches of the mp object to THB splines
    gsTHBSpline<2,real_t> thb;
    if (THB)
    {
      for (size_t k=0; k!=mpBspline.nPatches(); ++k)
      {
          gsTensorBSpline<2,real_t> *geo = dynamic_cast< gsTensorBSpline<2,real_t> * > (&mpBspline.patch(k));
          thb = gsTHBSpline<2,real_t>(*geo);
          mp.addPatch(thb);
      }

      gsMatrix<> refBoxes(2,2);
      if      (testCase==2 || testCase==3)
      {
        refBoxes.col(0) << 0.25,0.25;
        refBoxes.col(1) << 0.75,0.75;
      }
      else if (testCase==4 || testCase==5)
      {
        refBoxes.col(0) << 0.25,0.00;
        refBoxes.col(1) << 0.75,0.25;
      }
      else if (testCase==6 || testCase==7)
      {
        refBoxes.col(0) << 0.00,0.00;
        refBoxes.col(1) << 0.25,0.25;
      }

      int refExtension = 1;
      std::vector<index_t> elements = mp.patch(0).basis().asElements(refBoxes, refExtension);
      mp.patch(0).refineElements( elements );
    }
    else
      mp = mpBspline;

    gsMultiBasis<> dbasis(mp);
    gsInfo<<"Basis (patch 0): "<< mp.patch(0).basis() << "\n";

    // Boundary conditions
    gsBoundaryConditions<> BCs;
    BCs.setGeoMap(mp);

    std::string output = "solution";
    std::string dirname = "DynamicRelaxationResults";

    gsMatrix<> writePoints(2,3);
    writePoints.col(0)<< 0.0,0.5;
    writePoints.col(1)<< 0.5,0.5;
    writePoints.col(2)<< 1.0,0.5;
    index_t cross_coordinate = -1;
    real_t cross_val = 0.0;
    gsConstantFunction<> displ(0.0,3);
    gsConstantFunction<> disply(0.0,3);

    real_t Displ;

    if (testCase == 2 || testCase == 3)
    {
      Displ = -aDim/2;
      BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0 ,false,0);
      BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0 ,false,1);
      BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0 ,false,2);

      BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0 ,false,1);
      BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0 ,false,2);
      BCs.addCondition(boundary::east, condition_type::dirichlet, &displ, 0 ,false,2);

      if (weak)
      {
        BCs.addCondition(boundary::east, condition_type::weak_clamped, 0, 0, false, 2);
        BCs.addCondition(boundary::west, condition_type::weak_clamped, 0, 0, false, 2);
      }
      else
      {
        BCs.addCondition(boundary::east, condition_type::clamped  , 0, 0, false,2);
        BCs.addCondition(boundary::west, condition_type::clamped  , 0, 0, false,2);
      }

      std::stringstream ss;
      ss<<perturbation;
      dirname = dirname + "/FullSheet_Perturbed=" + ss.str() + "_r=" + std::to_string(numHref) + "_e=" + std::to_string(numElevate) + "_M=" + std::to_string(material) + "_c=" + std::to_string(Compressibility);

      output =  "solution";
      wn = output + "data.txt";
      cross_coordinate = 0;
      cross_val = 0.5;
    }
    else if (testCase == 4 || testCase == 5)
    {
      Displ = -aDim/2;
      BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0 ,false,0);
      BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0 ,false,1);
      BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0 ,false,2);

      BCs.addCondition(boundary::east, condition_type::dirichlet, &displ, 0 ,false,0);
      BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0 ,false,1);
      BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0 ,false,2);

      if (weak)
      {
        BCs.addCondition(boundary::east, condition_type::weak_clamped, 0, 0, false, 2);
        BCs.addCondition(boundary::west, condition_type::weak_clamped, 0, 0, false, 2);
      }
      else
      {
        BCs.addCondition(boundary::east, condition_type::clamped  , 0, 0, false,2);
        BCs.addCondition(boundary::west, condition_type::clamped  , 0, 0, false,2);
      }

      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 2 - z.
      if (symmetry)
        BCs.addCondition(boundary::south, condition_type::clamped, 0, 0, false, 2 ); // unknown 2 - z.
      else
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z.

      std::stringstream ss;
      ss<<perturbation;
      dirname = dirname + "/HalfSheet_Perturbed=" + ss.str() + "_r=" + std::to_string(numHref) + "_e=" + std::to_string(numElevate) + "_M=" + std::to_string(material) + "_c=" + std::to_string(Compressibility);

      output =  "solution";
      wn = output + "data.txt";
      cross_coordinate = 0;
      cross_val = 0.5;
    }
    else if (testCase == 6 || testCase == 7)
    {
      Displ = aDim/4;
      BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0 ,false,0);

      BCs.addCondition(boundary::east, condition_type::dirichlet, &displ, 0 ,false,0);
      BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0 ,false,1);
      BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0 ,false,2);

      if (weak)
      {
        BCs.addCondition(boundary::east, condition_type::weak_clamped, 0, 0, false, 2);
        BCs.addCondition(boundary::west, condition_type::weak_clamped, 0, 0, false, 2);
      }
      else
      {
        BCs.addCondition(boundary::east, condition_type::clamped  , 0, 0, false,2);
        BCs.addCondition(boundary::west, condition_type::clamped  , 0, 0, false,2);
      }

      BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 2 - z.
      if (symmetry)
        BCs.addCondition(boundary::south, condition_type::clamped, 0, 0, false, 2 ); // unknown 2 - z.
      else
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z.

      std::stringstream ss;
      ss<<perturbation;
      dirname = dirname + "/QuarterSheet_Perturbed=" + ss.str() + "_r=" + std::to_string(numHref) + "_e=" + std::to_string(numElevate) + "_M=" + std::to_string(material) + "_c=" + std::to_string(Compressibility);

      output =  "solution";
      wn = output + "data.txt";
      cross_coordinate = 0;
      cross_val = 0.0;
    }

    if (THB)
      dirname = dirname + "_THB";
    if (symmetry)
      dirname = dirname + "_symmetryBC";
    if (weak)
      dirname = dirname + "_weak";

    std::string commands = "mkdir -p " + dirname;
    const char *command = commands.c_str();
    system(command);

    // plot geometry
    if (plot)
      gsWriteParaview(mp,dirname + "/" + "mp",1000,true);

    if (writeG)
    {
      gsWrite(mp,dirname + "/" + "geometry");
      gsInfo<<"Geometry written in: " + dirname + "/" + "geometry.xml\n";
    }

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


    // Initialise solution object
    gsMultiPatch<> mp_def = mp;
    gsSparseSolver<>::LU solver;

    // Linear isotropic material model
    gsFunctionExpr<> force("0","0","0",3);
    gsFunctionExpr<> t(std::to_string(thickness), 3);
    gsFunctionExpr<> E(std::to_string(E_modulus),3);
    gsFunctionExpr<> nu(std::to_string(PoissonRatio),3);
    gsFunctionExpr<> rho(std::to_string(Density),3);
    gsConstantFunction<> ratio(Ratio,3);

    mu = E_modulus / (2 * (1 + PoissonRatio));
    gsConstantFunction<> alpha1(1.3,3);
    gsConstantFunction<> mu1(6.3e5/4.225e5*mu,3);
    gsConstantFunction<> alpha2(5.0,3);
    gsConstantFunction<> mu2(0.012e5/4.225e5*mu,3);
    gsConstantFunction<> alpha3(-2.0,3);
    gsConstantFunction<> mu3(-0.1e5/4.225e5*mu,3);

    real_t pi = math::atan(1)*4;
    index_t kmax = 1;
    gsVector<> E11(kmax), E22(kmax), G12(kmax), nu12(kmax), nu21(kmax), thick(kmax), phi(kmax);
    E11.setZero(); E22.setZero(); G12.setZero(); nu12.setZero(); nu21.setZero(); thick.setZero(); phi.setZero();
    for (index_t k=0; k != kmax; ++k)
    {
        E11.at(k) = E22.at(k) = E_modulus;
        nu12.at(k) = nu21.at(k) = PoissonRatio;
        G12.at(k) = 0.5 * E_modulus / (1+PoissonRatio);
        thick.at(k) = thickness/kmax;
        phi.at(kmax) = k / kmax * pi/2.0;
    }

    gsConstantFunction<> E11fun(E11,3);
    gsConstantFunction<> E22fun(E22,3);
    gsConstantFunction<> G12fun(G12,3);
    gsConstantFunction<> nu12fun(nu12,3);
    gsConstantFunction<> nu21fun(nu21,3);
    gsConstantFunction<> thickfun(thick,3);
    gsConstantFunction<> phifun(phi,3);

    std::vector<gsFunction<>*> parameters;
    if (material==0) // SvK & Composites
    {
      if (composite)
      {
        parameters.resize(6);
        parameters[0] = &E11fun;
        parameters[1] = &E22fun;
        parameters[2] = &G12fun;
        parameters[3] = &nu12fun;
        parameters[4] = &nu21fun;
        parameters[5] = &phifun;
      }
      else
      {
        parameters.resize(2);
        parameters[0] = &E;
        parameters[1] = &nu;
      }
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
            options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
            options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",0);
            materialMatrix = getMaterialMatrix<3,real_t>(mp,thickfun,parameters,rho,options);
        }
        else
        {
            parameters.resize(2);
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

    // Function for the Jacobian
    typedef std::function<gsSparseMatrix<real_t> (gsVector<real_t> const &)>    Jacobian_t;
    typedef std::function<gsVector<real_t> (gsVector<real_t> const &, real_t) >   ALResidual_t;
    Jacobian_t Jacobian = [&assembler,&mp_def](gsVector<real_t> const &x)
    {
      assembler->constructSolution(x,mp_def);
      assembler->assembleMatrix(mp_def);
      gsSparseMatrix<real_t> m = assembler->matrix();
      return m;
    };
    // Function for the Residual
    ALResidual_t ALResidual = [&displ,&BCs,&assembler,&mp_def](gsVector<real_t> const &x, real_t lam)
    {
        displ.setValue(lam,3);
        assembler->updateBCs(BCs);
        assembler->constructSolution(x,mp_def);
        assembler->assembleVector(mp_def);
        return assembler->rhs(); // - lam * force;
    };

    displ.setValue(1.0,3);
    assembler->updateBCs(BCs);
    // Assemble linear system to obtain the force vector
    assembler->assemble();
    gsSparseMatrix<> K = assembler->matrix();
    gsVector<> F = assembler->rhs();
    assembler->assembleMass(true);
    gsVector<> M = assembler->rhs();

    gsStaticDR<real_t> DRM(M,F,ALResidual);
    gsOptionList DROptions = DRM.options();
    DROptions.setReal("damping",damping);
    DROptions.setReal("alpha",alpha);
    DROptions.setInt("maxIt",maxit);
    DROptions.setReal("tol",tolF);
    DROptions.setReal("tolE",tolE);
    DROptions.setInt("verbose",verbose);
    DRM.setOptions(DROptions);

    gsParaviewCollection collection(dirname + "/" + output);
    gsParaviewCollection Smembrane(dirname + "/" + "membrane");
    gsParaviewCollection Sflexural(dirname + "/" + "flexural");
    gsParaviewCollection Smembrane_p(dirname + "/" + "membrane_p");
    gsMultiPatch<> deformation = mp;

    gsMatrix<> solVector;
    DRM.initialize();

    gsControlDisplacement<real_t> control(&DRM);

    real_t Load;
    gsVector<real_t> energies;

    if (energyPlot)
    {
#ifdef GISMO_WITH_MATPLOTLIB
        plt::figure();
#endif
    }

    for (index_t k=0; k<step; k++)
    {
      Load = Displ/step;
      gsInfo<<"Load step "<< k<<", Load = "<<Load<<"\n";
      control.step(Load);
      gsInfo<<"Step finished in "<<DRM.iterations()<<" iterations";

      energies = DRM.relEnergies();

      if (energyPlot)
      {
#ifdef GISMO_WITH_MATPLOTLIB
        std::vector<real_t> xa(DRM.iterations()+1);
        std::iota(std::begin(xa), std::end(xa), 0);

        std::vector<real_t> E(energies.data(), energies.data() + energies.rows() * energies.cols());

        plt::clf();
        plt::title("Kinetic Energy at step " + std::to_string(k));
        plt::semilogy(xa, E);
        plt::show();
        //plt::save("./poisson2_example.png");
#endif
      }

      solVector = control.solutionU();
      mp_def = assembler->constructSolution(solVector);

      // assembler->constructSolution(solVector,solution);

      deformation = mp_def;
      deformation.patch(0).coefs() -= mp.patch(0).coefs();// assuming 1 patch here


      if (plot)
      {
        gsField<> solField;
        if (deformed)
          solField= gsField<>(mp_def,deformation);
        else
          solField= gsField<>(mp,deformation);

        std::string fileName = dirname + "/" + output + util::to_string(k);
        gsWriteParaview<>(solField, fileName, 1000,mesh);
        fileName = output + util::to_string(k) + "0";
        collection.addTimestep(fileName,k,".vts");
        if (mesh) collection.addTimestep(fileName,k,"_mesh.vtp");
      }
      if (stress)
      {
        gsField<> membraneStress, flexuralStress, membraneStress_p;

        gsPiecewiseFunction<> membraneStresses;
        assembler->constructStress(mp_def,membraneStresses,stress_type::membrane);
        if (deformed)
          membraneStress = gsField<>(mp_def,membraneStresses,true);
        else
          membraneStress = gsField<>(mp,membraneStresses,true);

        gsPiecewiseFunction<> flexuralStresses;
        assembler->constructStress(mp_def,flexuralStresses,stress_type::flexural);
        if (deformed)
          flexuralStress = gsField<>(mp_def,flexuralStresses, true);
        else
          flexuralStress = gsField<>(mp,flexuralStresses, true);

        gsPiecewiseFunction<> membraneStresses_p;
        assembler->constructStress(mp_def,membraneStresses_p,stress_type::principal_stress_membrane);
        if (deformed)
          membraneStress_p = gsField<>(mp_def,membraneStresses_p, true);
        else
          membraneStress_p = gsField<>(mp,membraneStresses_p, true);

        std::string fileName;
        fileName = dirname + "/" + "membrane" + util::to_string(k);
        gsWriteParaview( membraneStress, fileName, 1000);
        fileName = "membrane" + util::to_string(k) + "0";
        Smembrane.addTimestep(fileName,k,".vts");

        fileName = dirname + "/" + "flexural" + util::to_string(k);
        gsWriteParaview( flexuralStress, fileName, 1000);
        fileName = "flexural" + util::to_string(k) + "0";
        Sflexural.addTimestep(fileName,k,".vts");

        fileName = dirname + "/" + "membrane_p" + util::to_string(k);
        gsWriteParaview( membraneStress_p, fileName, 1000);
        fileName = "membrane_p" + util::to_string(k) + "0";
        Smembrane_p.addTimestep(fileName,k,".vts");
      }



      if (write)
        writeStepOutput(solVector,Load,deformation, dirname + "/" + wn, writePoints,1, 201);

      if (crosssection && cross_coordinate!=-1)
        writeSectionOutput(deformation,dirname,cross_coordinate,cross_val,201,false);

      // DRM.setDisplacement(solVector);

    }

    if (plot)
    {
      collection.save();
    }
    if (stress)
    {
      Smembrane.save();
      Sflexural.save();
      Smembrane_p.save();
    }

    if (energyPlot)
    {
      #ifdef GISMO_WITH_MATPLOTLIB
        Py_Finalize();
      #endif
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
         geo->insertKnot(1 - std::min(offset, dknot1 / 2),1);
       else if (*it==boundary::south) // south
         geo->insertKnot(std::min(offset, dknot1 / 2),1);
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
        for (index_t k=0; k < points.cols(); k++)
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
void writeStepOutput(const gsMatrix<T> & solution, const T & load, const gsMultiPatch<T> & deformation, const std::string name, const gsMatrix<T> & points, const index_t extreme, const index_t kmax) // extreme: the column of point indices to compute the extreme over (default -1)
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
          << solution.norm() << ",";
          for (index_t p=0; p!=points.cols(); p++)
          {
            file<< out(0,p) << ","
                << out(1,p) << ","
                << out(2,p) << ",";
          }

    file  << load << ","
          << 0.0 << ","
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
          << solution.norm() << ",";
          for (index_t p=0; p!=points.cols(); p++)
          {
            file<< out(0,p) << ","
                << out(1,p) << ","
                << std::max(abs(out2.col(p).maxCoeff()),abs(out2.col(p).minCoeff())) << ",";
          }

    file  << load << ","
          << 0.0 << ","
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
