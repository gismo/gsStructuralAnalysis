/** @file DRMethoid.cpp

    @brief Static simulations of a shell

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#include <gismo.h>

#include <gsKLShell/gsThinShellAssembler.h>
#include <gsKLShell/getMaterialMatrix.h>

#include <gsElasticity/gsWriteParaviewMultiPhysics.h>

#include <gsStructuralAnalysis/gsStaticSolver.h>
#include <gsStructuralAnalysis/gsTimeIntegrator.h>

#include <gsSpectra/gsSpectra.h>

//#include <gsThinShell/gsNewtonIterator.h>

using namespace gismo;

template <class T>
gsMultiPatch<T> Rectangle(T L, T B);

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
class gsDynamicRelaxationALM
{
public:
    /// Constructor given the matrices and timestep
    gsDynamicRelaxationALM( const gsSparseMatrix<T> & K,
                            const gsVector<T> & F,
                            const std::function < gsVector<T> ( gsVector<T> const & , T, gsVector<T> const &) > &Residual
                        )
    :
    m_stif(K),
    m_forcing(F),
    m_residualFun(Residual)
    {
        m_dt = 1.0;
        m_dofs = m_stif.rows();
        m_V = m_U = m_DeltaU = m_deltaU = m_DeltaUold = gsVector<T>(m_dofs);
        m_V.setZero();
        m_U.setZero();
        m_DeltaU.setZero();
        m_deltaU.setZero();
        m_DeltaUold.setZero();

        m_L = m_DeltaL = m_deltaL = m_DeltaLold = 0.0;
        m_phi = 1.0;

        defaultOptions();
    }

    int sign(T val)
    {
        return (T(0) < val) - (val < T(0));
    }

    gsVector<T> displacements() const {return m_U+m_DeltaU;}
    gsVector<T> velocities() const {return m_V;}

    void setDisplacement(gsVector<T> displacement) { m_U = displacement; }

    T kineticEnergy() const { return m_Ek; }

    void init()
    {
        getOptions();
        m_mass = m_alpha * m_dt*m_dt / 2. * m_stif;
        m_massInv = m_mass.toDense().inverse();
        m_damp = m_c * m_mass;
    }

    void predictor()
    {
        m_DeltaU.setZero();
        m_V.setZero();
        m_DeltaL = m_deltaL = 0.0;

        m_deltaUt = 0.5*m_dt*m_dt*m_massInv*m_forcing;

        if (m_DeltaUold.dot(m_DeltaUold) == 0 && m_DeltaLold*m_DeltaLold == 0) // no information about previous step.
        {
            T DL = 1.;
            m_deltaL = m_arcLength * DL / math::sqrt( ( m_deltaUt.dot(m_deltaUt) + m_DeltaL*DL ) );

            m_deltaU = m_deltaL*m_deltaUt;

            m_DeltaUold = m_deltaU;
            m_DeltaLold = m_deltaL;
        }
        else
        {
            m_phi = math::pow(m_U.dot(m_U)/( math::pow(m_L,2) * m_forcing.dot(m_forcing) ),0.5);
            T A0 = math::pow(m_phi,2)* m_forcing.dot(m_forcing); // see Lam et al. 1991
            int dir = sign(m_DeltaUold.dot(m_deltaUt) + A0*m_DeltaLold); // Feng et al. 1995 with H = \Psi^2
            T denum = ( math::pow( m_deltaUt.dot(m_deltaUt) + A0 ,0.5) ); // Feng et al. 1995 with H = \Psi^2

            T mu;
            if (denum==0)
              mu = m_arcLength;
            else
              mu = m_arcLength / denum;

            m_deltaL = dir*mu;
            m_deltaU = m_deltaL*m_deltaUt;
        }
        m_V = m_deltaU/m_dt;

        m_DeltaU += m_deltaU;
        m_DeltaL += m_deltaL;
    }

    void step()
    {
        m_R = m_residualFun(m_U+m_DeltaU,m_L+m_DeltaL,m_forcing) - m_damp*m_V;
        m_V += m_dt * m_massInv * m_R;                    // Velocities at t+dt/2
        m_deltaUbar = m_dt* (m_massInv * m_R + m_V);               // Velocities at t+dt

        T A0 = math::pow(m_phi,2)* m_forcing.dot(m_forcing); // see Lam et al. 1991,

        T a0 = m_deltaUt.dot(m_deltaUt) + A0;
        T b0 = 2*( m_deltaUt.dot(m_DeltaU) + m_DeltaL * A0 );
        T b1 = 2*( m_deltaUbar.dot(m_deltaUt) );
        T c0 = m_DeltaU.dot(m_DeltaU) + m_DeltaL*m_DeltaL * A0 - math::pow(m_arcLength,2);
        T c1 = 2*( m_DeltaU.dot(m_deltaUbar) );
        T c2 = m_deltaUbar.dot(m_deltaUbar);

        /// Calculate the coefficients of the polynomial
        T alpha1 = a0;
        T alpha2 = b0 + b1;
        T alpha3 = c0 + c1 + c2;

        T discriminant = math::pow(alpha2 ,2) - 4 * alpha1 * alpha3;

        gsVector<T> deltaLs(2);
        deltaLs.setZero();
        if (discriminant >= 0)
        {
          deltaLs[0] = (-alpha2 + math::pow(discriminant,0.5))/(2*alpha1);
          deltaLs[1] = (-alpha2 - math::pow(discriminant,0.5))/(2*alpha1);

          gsVector<T> deltaU1, deltaU2;
          deltaU1 = m_deltaUbar + m_deltaUt*deltaLs[0];
          deltaU2 = m_deltaUbar + m_deltaUt*deltaLs[1];

          // ---------------------------------------------------------------------------------
          // Method by Ritto-Corea et al. 2008
          T DOT1,DOT2;
          DOT1 = deltaLs[0]*(m_DeltaUold.dot(m_deltaUt) + math::pow(m_phi,2)*m_DeltaLold);
          DOT2 = deltaLs[1]*(m_DeltaUold.dot(m_deltaUt) + math::pow(m_phi,2)*m_DeltaLold);

          if (DOT1 > DOT2)
          {
            m_deltaL = deltaLs[0];
            m_deltaU = deltaU1;
          }
          else if (DOT1 < DOT2)
          {
            m_deltaL = deltaLs[1];
            m_deltaU = deltaU2;
          }
          else
          {
            m_deltaL = deltaLs[0];
            m_deltaU = deltaU1;
          }
        }
        else


        m_DeltaU += m_deltaU;
        m_DeltaL += m_deltaL;


        m_Ek_prev = m_Ek;
        m_R = m_residualFun(m_U+m_DeltaU,m_L+m_DeltaL,m_forcing) - m_damp*m_V;
        m_V += m_dt * m_massInv * m_R;                    // Velocities at t+dt/2
        m_deltaU = m_dt * m_V;               // Velocities at t+dt
        m_DeltaU += m_deltaU;               // Velocities at t+dt
        m_Ek = m_V.transpose() * m_mass * m_V;
    }

    // void reset()
    // {
    //     start();
    // }

    void peak()
    {
        m_R = m_residualFun(m_U+m_DeltaU,m_L+m_DeltaL,m_forcing)- m_damp*m_V;
        m_deltaU = - 1.5 * m_dt * m_V + m_dt*m_dt / 2. * m_massInv * m_R;
        m_DeltaU += m_deltaU;

        m_V.setZero();
        // m_V = 0.5 * m_dt * m_massInv * m_R; // Velocities at dt/2
    }

    void finish()
    {
        m_U += m_DeltaU;
        m_L += m_DeltaL;
        m_DeltaUold = m_DeltaU;
        m_DeltaLold = m_DeltaL;
    }

    void defaultOptions()
    {
        m_options.addReal("damping","damping factor",1.0);
        m_options.addReal("alpha","mass coefficient",2.0);
        m_options.addReal("length","arc-length",1.0);
    }

    gsOptionList options() const {return m_options;}

    void getOptions()
    {
        m_c = m_options.getReal("damping");
        m_alpha = m_options.getReal("alpha");
        m_arcLength = m_options.getReal("length");
    }

    void setOptions(gsOptionList & options) {m_options.update(options,gsOptionList::addIfUnknown); }

    T residualNorm() const { return m_R.norm(); }

protected:
    const gsSparseMatrix<T> m_stif;
    const gsVector<T> m_forcing;
    const std::function<gsVector<T> ( gsVector<T> const &, T, gsVector<T> const &) > m_residualFun;
    gsVector<T> m_U, m_V;
    gsVector<T> m_DeltaU, m_deltaU, m_deltaUt, m_deltaUbar;
    gsVector<T> m_DeltaUold;
    T m_DeltaLold;
    gsVector<T> m_R;
    gsSparseMatrix<T> m_mass, m_damp;
    gsMatrix<T> m_massInv;
    index_t m_dofs;
    T m_dt, m_alpha, m_c;
    gsOptionList m_options;
    T m_Ek, m_Ek_prev;
    T m_L, m_DeltaL, m_deltaL;
    T m_arcLength;
    T m_phi;
};


// Choose among various shell examples, default = Thin Plate
int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot  = false;
    bool stress= false;
    index_t numRefine  = 1;
    index_t numElevate = 1;
    index_t testCase = 1;

    bool nonlinear = false;
    int verbose = 0;
    std::string fn;
    bool membrane = false;
    bool weak = false;


    index_t Compressibility = 0;
    index_t material = 0;
    bool composite = false;
    index_t impl = 1; // 1= analytical, 2= generalized, 3= spectral

    real_t E_modulus = 1.0;
    real_t PoissonRatio = 0.0;
    real_t Density = 1.0;
    real_t thickness = 1.0;
    real_t Ratio = 7.0;

    real_t alpha = 1.0;
    real_t damping = 0.1;

    std::string assemberOptionsFile("options/solver_options.xml");

    real_t dt = 0.1;

    gsCmdLine cmd("Static analysis for thin shells.");
    cmd.addString( "f", "file", "Input XML file for assembler options", assemberOptionsFile );
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving",  numRefine );

    cmd.addReal( "d", "dt", "dt",  dt );
    cmd.addReal( "a", "alpha", "alpha",  alpha );
    cmd.addReal( "c", "damping", "damping",  damping );
    cmd.addInt( "m", "Material", "Material law",  material );
    cmd.addSwitch("nl", "Solve nonlinear problem", nonlinear);
    cmd.addInt("v","verbose", "0: no; 1: iteration output; 2: Full matrix and vector output", verbose);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("stress", "Create a ParaView visualization file with the stresses", stress);
    cmd.addSwitch("membrane", "Use membrane model (no bending)", membrane);
    cmd.addSwitch("weak", "Impose boundary conditions weakly", weak);



    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    // std::string methodName;
    // if (method==1)
    //   methodName = "ExplEuler";
    // else if (method==2)
    //   methodName = "ImplEuler";
    // else if (method==3)
    //   methodName = "Newmark";
    // else if (method==4)
    //   methodName = "Bathe";
    // else if (method==5)
    //   methodName = "CentralDiff";
    // else if (method==6)
    //   methodName = "RK4";


    //! [Read input file]
    gsFileData<> fd(assemberOptionsFile);
    gsOptionList opts;
    fd.getFirst<gsOptionList>(opts);

    gsMultiPatch<> mp;
    gsMultiPatch<> mp_def;

    E_modulus = 1;
    // thickness = 0.15;
    thickness = 1;
    if (!Compressibility)
      PoissonRatio = 0.499;
    else
      PoissonRatio = 0.45;

    E_modulus = 1;

    real_t bDim = 1;
    real_t aDim = 1;

    mp = Rectangle(aDim, bDim);
    // p-refine
    if (numElevate!=0)
        mp.degreeElevate(numElevate);

    // h-refine
    for (int r =0; r < numRefine; ++r)
        mp.uniformRefine();

    mp_def = mp;
    gsWriteParaview<>( mp_def    , "mp", 1000, true);


    //! [Refinement]
    gsMultiBasis<> dbasis(mp);

    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
    gsInfo<<mp_def<<"\n";
    gsInfo << dbasis.basis(0)<<"\n";

    gsBoundaryConditions<> bc;
    bc.setGeoMap(mp);

    gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();
    for (index_t i=0; i!=3; ++i)
    {
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, i ); // unknown 2 - z
    }
    bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
    bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

    bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0 ,false,1);
    bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0 ,false,2);
    bc.addCondition(boundary::east, condition_type::collapsed, 0, 0 ,false,0);

    gsVector<> point(2); point<< 1.0, 0.5 ;
    gsVector<> load (3); load << 0.125, 0.0, 0.0 ;
    pLoads.addLoad(point, load, 0 );

    gsVector<> tmp(3);
    tmp<<0,0,0;

    //! [Refinement]
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

    real_t pi = math::atan(1)*4;
    index_t kmax = 2;
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

    // gsConstantFunction<> E11fun(E11,3);
    // gsConstantFunction<> E22fun(E22,3);
    // gsConstantFunction<> G12fun(G12,3);
    // gsConstantFunction<> nu12fun(nu12,3);
    // gsConstantFunction<> nu21fun(nu21,3);
    // gsConstantFunction<> thickfun(thick,3);
    // gsConstantFunction<> phifun(phi,3);
    gsFunctionExpr<> E11fun(std::to_string(E_modulus),std::to_string(E_modulus),3);
    gsFunctionExpr<> E22fun(std::to_string(E_modulus),std::to_string(E_modulus),3);
    gsFunctionExpr<> G12fun(std::to_string(0.5 * E_modulus / (1+PoissonRatio)),std::to_string(0.5 * E_modulus / (1+PoissonRatio)),3);
    gsFunctionExpr<> nu12fun(std::to_string(PoissonRatio),std::to_string(PoissonRatio),3);
    gsFunctionExpr<> nu21fun(std::to_string(PoissonRatio),std::to_string(PoissonRatio),3);
    gsFunctionExpr<> thickfun(std::to_string(thickness/kmax),std::to_string(thickness/kmax), 3);
    gsFunctionExpr<> phifun("0","0", 3);

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
    assembler = new gsThinShellAssembler<3, real_t, true >(mp,dbasis,bc,force,materialMatrix);

    opts.addReal("WeakDirichlet","Penalty parameter weak dirichlet conditions",1e5);
    opts.addReal("WeakClamped","Penalty parameter weak clamped conditions",1e5);
    // Construct assembler object
    assembler->setOptions(opts);
    assembler->setPointLoads(pLoads);

    // Function for the Jacobian
    typedef std::function<gsSparseMatrix<real_t> (gsVector<real_t> const &)>    Jacobian_t;
    typedef std::function<gsVector<real_t> (gsVector<real_t> const &) >         Residual_t;
    typedef std::function<gsVector<real_t> (gsVector<real_t> const &, real_t, gsVector<real_t> const &) >   ALMResidual_t;
    Jacobian_t Jacobian = [&assembler,&mp_def](gsVector<real_t> const &x)
    {
      assembler->constructSolution(x,mp_def);
      assembler->assembleMatrix(mp_def);
      gsSparseMatrix<real_t> m = assembler->matrix();
      return m;
    };
    // Function for the Residual
    Residual_t Residual = [&assembler,&mp_def](gsVector<real_t> const &x)
    {
      assembler->constructSolution(x,mp_def);
      assembler->assembleVector(mp_def);
      return assembler->rhs();
    };

    // Function for the Residual
    ALMResidual_t ALMResidual = [&assembler,&mp_def](gsVector<real_t> const &x, real_t lam, gsVector<real_t> const &force)
    {
        assembler->constructSolution(x,mp_def);
        assembler->assembleVector(mp_def);
        gsVector<real_t> Fint = -(assembler->rhs() - force);
        gsVector<real_t> result = Fint - lam * force;
        return result; // - lam * force;
    };

    assembler->assemble();
    gsSparseMatrix<> K = assembler->matrix();
    gsVector<> F = assembler->rhs();

    gsParaviewCollection collection("incr_solution");

    gsDynamicRelaxationALM<real_t> DRM(K,F,ALMResidual);
    gsOptionList DROptions = DRM.options();
    DROptions.setReal("damping",damping);
    DROptions.setReal("alpha",alpha);
    DRM.setOptions(DROptions);
    DRM.init();


    index_t steps = 1;
    for (index_t k=0; k!=steps; k++)
    {
        gsVector<> displacements;
        DRM.predictor();
        real_t Ek = DRM.kineticEnergy();
        real_t EkOld = Ek;
        for (index_t i=1; i<1e2; i++)
        {
            displacements = DRM.displacements();
            mp_def = assembler->constructSolution(displacements);
            gsMultiPatch<> deformation = mp_def;
            deformation.patch(0).coefs() -= mp.patch(0).coefs();// assuming 1 patch here
            gsField<> solField(mp,deformation);
            std::string fileName = "incr_solution" + util::to_string(i);
            gsWriteParaview<>(solField, fileName, 500);
            fileName = "incr_solution" + util::to_string(i) + "0";
            collection.addTimestep(fileName,i,".vts");

            gsInfo<<"Step "<<i<<"\t";
            DRM.step();
            if (damping==0)
            {
                if (EkOld > DRM.kineticEnergy())
                    DRM.peak();
                EkOld = DRM.kineticEnergy();
            }
            gsInfo<<"Res norm = "<<DRM.residualNorm()/F.norm()<<"\t Kin energy = "<<DRM.kineticEnergy()/Ek<<"\n";

            if (DRM.residualNorm()/F.norm() < 1e-10 && DRM.kineticEnergy()/Ek < 1e-3)
            {
                DRM.finish();
                break;
            }


        }
        displacements = DRM.displacements();
        DRM.setDisplacement(displacements);
    }



    collection.save();

    gsMultiPatch<> deformation = mp_def;
    for (size_t k = 0; k != mp_def.nPatches(); ++k)
        deformation.patch(k).coefs() -= mp.patch(k).coefs();

    // ! [Export visualization in ParaView]
    if (plot)
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

    return EXIT_SUCCESS;

}// end main


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
