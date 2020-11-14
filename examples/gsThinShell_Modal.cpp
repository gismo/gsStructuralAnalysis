/** @file gsThinShell_Buckling.cpp

    @brief Example to compute eigenvalues and eigenmodes of a buckled shell

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gismo.h>
#include <gismo_dev.h>

#include <gsIO/gsMatrixToFile.h>
#include <gsThinShell2/gsThinShellAssembler.h>

#include <gsStructuralAnalysis/gsModalSolver.h>

using namespace gismo;

void writeToCSVfile(std::string name, gsMatrix<> matrix)
{
    std::ofstream file(name.c_str());
    for(int  i = 0; i < matrix.rows(); i++){
        for(int j = 0; j < matrix.cols(); j++){
           std::string str = std::to_string(matrix(i,j));
           if(j+1 == matrix.cols()){
               file<<std::setprecision(10)<<str;
           }else{
               file<<std::setprecision(10)<<str<<',';
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
    bool sparse     = false;
    bool nonlinear  = false;
    bool first  = false;
    int mode = 0;

    real_t thickness     = 1;
    real_t width = 1; // Width of the strip is equal to 1.
    real_t length = 10; // Length of the strip is equal to 10.
    real_t Area = thickness*width;

    real_t E_modulus     = 1e0;
    real_t PoissonRatio = 0;
    real_t Density = 1e0;
    gsMultiPatch<> mp;

    real_t shift = 0.0;

    int testCase = 0;

    int result = 0;

    bool write = false;

    gsCmdLine cmd("Thin shell plate example.");
    cmd.addInt("r","hRefine",
               "Number of dyadic h-refinement (bisection) steps to perform before solving",
               numHref);
    cmd.addInt("t", "testcase",
                "Test case: 0: Beam - pinned-pinned, 1: Beam - fixed-fixed, 2: beam - fixed-free, 3: plate - fully pinned, 4: plate - fully fixed, 5: circle - fully pinned, 6: 5: circle - fully fixed",
               testCase);
    cmd.addInt("m", "nmode",
               "Mode shape number, starting from 0",
              mode);
    cmd.addInt("e","degreeElevation",
               "Number of degree elevation steps to perform on the Geometry's basis before solving",
               numElevate);
    cmd.addReal("s","shift", "eigenvalue shift", shift);
    cmd.addSwitch("nl", "Nonlinear elasticity (otherwise linear)", nonlinear);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addSwitch("first", "Plot only first", first);
    cmd.addSwitch("write", "Write convergence data to file", write);
    cmd.addSwitch("sparse", "Use sparse solver", sparse);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    std::string fn("planar/strip.xml");


    real_t EA,EI,r,D;
    if (testCase==0 || testCase==1 || testCase==2)
    {
        fn = "planar/strip.xml";
        EA = E_modulus*Area;
        EI = 1.0/12.0*(width*math::pow(thickness,3))*E_modulus;
        r = math::sqrt(EI/EA);
        gsInfo<<"EI = "<<EI<<"; EA = "<<EA<<"; r = "<<r<<"\n";
    }
    else if (testCase==3 || testCase==4)
    {
      fn = "planar/unitplate.xml";
    }

    else if (testCase==5 || testCase==6)
    {
      fn = "planar/unitcircle.xml";
    }

    gsReadFile<>(fn, mp);

    for(index_t i = 0; i< numElevate; ++i)
        mp.patch(0).degreeElevate();    // Elevate the degree

    // h-refine
    for(index_t i = 0; i< numHref; ++i)
        mp.patch(0).uniformRefine();

    gsMultiBasis<> dbasis(mp);


    gsInfo<<"Basis (patch 0): "<< mp.patch(0).basis() << "\n";

    gsWriteParaview<>(mp, "mp", 500,true,false);

    // Boundary conditions
    std::vector< std::pair<patchSide,int> > clamped;
    gsBoundaryConditions<> BCs;
    gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();

    // Initiate Surface forces
    std::string tx("0");
    std::string ty("0");
    std::string tz("0");

    gsVector<> tmp(3);
    tmp << 0, 0, 0;

    std::vector<real_t> omegas;

    if (testCase == 0)
    {
        // Beam
        // Pinned-Pinned [no loads]
            // Left
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
            // Right
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

        // Restrain sides
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y

        for (index_t n=0; n!=10; n++)
          omegas.push_back(math::pow((n*3.1415926535/length),2)*math::pow(EI/(Density*Area),0.5));

    }
    else if (testCase == 1)
    {
        // Beam
        // Fixed-Fixed
            // Left
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
        BCs.addCondition(boundary::west, condition_type::clamped,0,0,false,2);
            // Right
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
        BCs.addCondition(boundary::east, condition_type::clamped,0,0,false,2);

        // Restrain sides
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y

        omegas.push_back(math::pow((4.73004/length),2)*math::pow(EI/(Density*Area),0.5));

    }
    else if (testCase == 2)
    {
        // Beam
        // Fixed-Free
            // Left
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
        BCs.addCondition(boundary::west, condition_type::clamped,0,0,false,2);

        // Restrain sides
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        // Restrain tip in horizontal and side direction
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y

        omegas.push_back(3.52/(2*3.1415926535)*math::pow((1/length),2)*math::pow(EI/(Density*Area),0.5));
    }
    else if (testCase == 3)
    {
        thickness = 0.01;
        PoissonRatio = 0.3;
        // Plate
        // Pinned-Pinned-Pinned-Pinned
            // Left
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
            // Right
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
            // Top
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
            // Bottom
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z


        real_t D = E_modulus*math::pow(thickness,3)/(12*(1-math::pow(PoissonRatio,2)));
        gsInfo<<"D = "<<D<<"\n";
        for (index_t m=1; m!=10; m++)
          for (index_t n=1; n!=10; n++)
            omegas.push_back((math::pow(m/1.0,2)+math::pow(n/1.0,2))*math::pow(3.1415926535,2)*math::sqrt(D / (Density * thickness)));

        std::sort(omegas.begin(),omegas.end());
    }
    else if (testCase == 4)
    {
        thickness = 0.01;
        PoissonRatio = 0.3;
        // Plate
        // Pinned-Pinned-Pinned-Pinned
            // Left
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
        BCs.addCondition(boundary::west, condition_type::clamped,0,0,false,2);

            // Right
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
        BCs.addCondition(boundary::east, condition_type::clamped,0,0,false,2);

            // Top
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
        BCs.addCondition(boundary::north, condition_type::clamped,0,0,false,2);

            // Bottom
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
        BCs.addCondition(boundary::south, condition_type::clamped,0,0,false,2);


        // real_t D = E_modulus*math::pow(thickness,3)/(12*(1-math::pow(PoissonRatio,2)));
        // gsInfo<<"D = "<<D<<"\n";
        // for (index_t m=0; m!=10; m++)
        //   for (index_t n=0; n!=10; n++)
        //     omegas.push_back((math::pow(m/1.0,2)+math::pow(n/1.0,2))*math::pow(3.1415926535,2)*math::sqrt(D / (Density * thickness)));

        // std::sort(omegas.begin(),omegas.end());

        omegas.push_back(0);
    }
    else if (testCase == 5)
    {
        thickness = 0.01;
        PoissonRatio = 0.3;
        // Circle
        // Pinned-Pinned-Pinned-Pinned
            // Left
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

            // Right
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

            // Top
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

            // Bottom
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z


        omegas.push_back(0);
    }
    else if (testCase == 6)
    {
        thickness = 0.01;
        PoissonRatio = 0.3;
        // Circle
        // Pinned-Pinned-Pinned-Pinned
            // Left
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
        BCs.addCondition(boundary::west, condition_type::clamped,0,0,false,2);

            // Right
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
        BCs.addCondition(boundary::east, condition_type::clamped,0,0,false,2);

            // Top
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
        BCs.addCondition(boundary::north, condition_type::clamped,0,0,false,2);

            // Bottom
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
        BCs.addCondition(boundary::south, condition_type::clamped,0,0,false,2);

        real_t D = E_modulus*math::pow(thickness,3)/(12*(1-math::pow(PoissonRatio,2)));
        gsInfo<<"D = "<<D<<"\n";
        gsVector<> gammas(10);
        gammas<<3.1962206158252, 4.6108999, 5.9056782, 6.3064370, 7.1435310, 7.7992738, 8.3466059, 9.1968826, 9.4394991, 9.5257014;
        for (index_t n=0; n!=gammas.size(); n++)
          omegas.push_back(math::pow(math::pow(gammas[n],4)*D/(Density*thickness),0.5));
    }


    gsFunctionExpr<> surfForce(tx,ty,tz,3);
    // Initialise solution object
    gsMultiPatch<> mp_def = mp;
    gsSparseSolver<>::LU solver;

    // Linear isotropic material model
    gsConstantFunction<> force(tmp,3);
    gsFunctionExpr<> t(std::to_string(thickness), 3);
    gsFunctionExpr<> E(std::to_string(E_modulus),3);
    gsFunctionExpr<> nu(std::to_string(PoissonRatio),3);
    gsFunctionExpr<> rho(std::to_string(Density),3);

    std::vector<gsFunction<>*> parameters(2);
    parameters[0] = &E;
    parameters[1] = &nu;
    gsMaterialMatrix materialMat(mp,mp_def,t,parameters,rho);

    gsThinShellAssembler assembler(mp,dbasis,BCs,surfForce,materialMat);
    assembler.setPointLoads(pLoads);

    // Initialise solution object
    gsMultiPatch<> solution = mp;

    assembler.assemble();
    gsSparseMatrix<> K =  assembler.matrix();
    assembler.assembleMass();
    gsSparseMatrix<> M =  assembler.matrix();

    gsModalSolver<real_t> modal(K,M);
    modal.verbose();

    if (!sparse)
      modal.compute();
    else
      modal.computeSparse(shift,10);


    gsMatrix<> values = modal.values();
    gsMatrix<> vectors = modal.vectors();


    gsInfo<<"First eigenfrequency: "<<"\t[Analytical]: "<< omegas[0]<<"\t[Numerical]: "<<math::sqrt(values.at(0))<<"\n";


    gsInfo<< "First 10 eigenvalues:\n";
    for (index_t k = 0; k<10; k++)
        gsInfo<<"\t"<<std::setprecision(20)<<values.at(k)<<"\n";
    gsInfo<<"\n";

    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";
        system("mkdir -p ModalResults");
        gsMultiPatch<> deformation = solution;
        gsMatrix<> modeShape;
        gsParaviewCollection collection("ModalResults/modes");

        int N = 1;
        if (!first)
          N = vectors.cols();
        for (index_t m=0; m<N; m++)
        {

          // Compute solution based on eigenmode with number 'mode'
          modeShape = modal.vector(m);//solver.solve( assembler.rhs() );
          assembler.constructSolution(modeShape, solution);

          // compute the deformation spline
          deformation = solution;
          deformation.patch(0).coefs() -= mp.patch(0).coefs();// assuming 1 patch here

          // Normalize mode shape amplitude in z coordinate
          real_t maxAmpl = std::max(math::abs(deformation.patch(0).coefs().col(2).maxCoeff()),math::abs(deformation.patch(0).coefs().col(2).minCoeff()));
          if (maxAmpl!=0.0)
          {
            deformation.patch(0).coefs() = deformation.patch(0).coefs()/maxAmpl;
          }

          gsField<> solField(mp,deformation);
          std::string fileName = "ModalResults/modes" + util::to_string(m);
          gsWriteParaview<>(solField, fileName, 5000);
          fileName = "modes" + util::to_string(m) + "0";
          collection.addTimestep(fileName,m,".vts");
        }
        collection.save();
    }

    if (write)
    {
        system("mkdir -p ModalResults");
        std::string wnM = "ModalResults/eigenvalues.txt";
        writeToCSVfile(wnM,values);
    }

    return result;
}