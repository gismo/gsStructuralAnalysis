/** @file gsTimeIntegrator.h

    @brief Provides temporal solvers for structural analysis problems


    TO DO:
    * [V] Fix Implicit Euler Nonlinear
    * [ ] Fix Implicit Euler BlockOp
    * [V] Fix Explicit Euler with nonlinearities (NOT VERIFIED)
    * [V] Fix linear Newmark and Bathe
    * [ ] Fix Implicit/Explicit Linear euler

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/


#pragma once


#include <gsSolver/gsBlockOp.h>
#include <gsSolver/gsMatrixOp.h>

namespace gismo
{

// enum class TMethod : short_t
// {
//     ExplEuler = 0;
//     ImplEuler = 1,
//     Newmark = 2,
//     Bathe = 3
// }


/**
    @brief Performs Newton iterations to solve a nonlinear equation system.

    \tparam T coefficient type

    \ingroup ThinShell
*/
template <class T>
class gsTimeIntegrator
{
public:

    /// Constructor given the matrices and timestep
    gsTimeIntegrator(   gsSparseMatrix<T> &Mass,
                        gsSparseMatrix<T> &Damp,
                        gsSparseMatrix<T> &Stif,
                        std::function < gsMatrix<T> ( T) > &Force,
                        T dt
                        )
    : m_mass(Mass),
      m_damp(Damp),
      m_stif(Stif),
      m_forceFun(Force),
      m_dt(dt)
    {
        m_dofs = m_mass.cols();
        m_method = "ExplEuler";
        m_NL = false;
        this->initializeCoefficients();
        this->initiate();
    }

    /// Constructor given the matrices and timestep
    gsTimeIntegrator(   gsSparseMatrix<T> &Mass,
                        gsSparseMatrix<T> &Stif,
                        std::function < gsMatrix<T> ( T) > &Force,
                        T dt
                        )
    : m_mass(Mass),
      m_stif(Stif),
      m_forceFun(Force),
      m_dt(dt)
    {
        m_dofs = m_mass.cols();
        m_damp = gsSparseMatrix<T>(m_dofs,m_dofs);
        m_method = "ExplEuler";
        m_NL = false;
        this->initializeCoefficients();
        this->initiate();
    }

    /// Constructor given the matrices and timestep
    gsTimeIntegrator(   gsSparseMatrix<T> &Mass,
                        gsSparseMatrix<T> &Damp,
                        std::function < gsSparseMatrix<T> ( gsMatrix<T> const & ) > &Jacobian,
                        std::function < gsMatrix<T> ( gsMatrix<T> const &, T) > &Residual,
                        T dt
                        )
    : m_mass(Mass),
      m_damp(Damp),
      m_jacobian(Jacobian),
      m_residualFun(Residual),
      m_dt(dt)
    {
        m_dofs = m_mass.cols();
        m_method = "Newmark";
        m_NL = true;
        this->initializeCoefficients();
        this->initiate();
    }

    /// Constructor given the matrices and timestep
    gsTimeIntegrator(   gsSparseMatrix<T> &Mass,
                        std::function<gsSparseMatrix<T> (gsMatrix<real_t> const & )> Jacobian,
                        std::function<gsMatrix<T> (gsMatrix<T> const &, T) > Residual,
                        T dt
                        )
    : m_mass(Mass),
      m_jacobian(Jacobian),
      m_residualFun(Residual),
      m_dt(dt)
    {
        m_dofs = m_mass.cols();
        m_damp = gsSparseMatrix<T>(m_dofs,m_dofs);
        m_method = "Newmark";
        m_NL = true;
        this->initializeCoefficients();
        this->initiate();
    }

public:

    void step();

    void constructSystem();

public:

    /// Tells if the Newton method converged
    bool converged() const {return m_converged;}

    /// Returns the number of Newton iterations performed
    index_t numIterations() const { return m_numIterations;}

    /// Returns the tolerance value used
    T tolerance() const {return m_tolerance;}

    /// Returns the error after solving the nonlinear system
    T residue()   const {return m_residue;}

    /// Set the maximum number of Newton iterations allowed
    void setMaxIterations(index_t nIter) {m_maxIterations = nIter;}

    /// Set the tolerance for convergence
    void setTolerance(T tol) {m_tolerance = tol;}

    void setTime(T time) {m_t = time; }
    T currentTime() const {return m_t; }

    /// set mass matrix
    void setMassMatrix(gsMatrix<T>& Mass) {m_mass = Mass; }
    void setDampingMatrix(gsMatrix<T>& Damp) {m_damp = Damp; }
    void setStiffnessMatrix(gsMatrix<T>& Stif) {m_stif = Stif; }
    void setJacobian(std::function < gsSparseMatrix<T> ( gsMatrix<T> const & ) > &Jacobian) {m_jacobian = Jacobian; }

    // set solutions
    void setDisplacement(gsMatrix<T>& displ);
    void setVelocity(gsMatrix<T>& velo);
    void setAcceleration(gsMatrix<T>& accel);

    // set time integration method
    void setMethod(std::string method);

    /// verbose
    void verbose() {m_verbose = true; }

    void quasiNewton() {m_quasiNewton = true; }

    const gsMatrix<T>& displacements() const { return m_displacements;}
    const gsMatrix<T>& velocities() const { return m_velocities;}
    const gsMatrix<T>& accelerations() const { return m_accelerations;}

    void constructSolution();

    void resetIteration();

    /// set post-processing routine
    // void setPostProcessing( std::function

protected:

    /// Linear solver employed
    gsSparseSolver<>::LU  m_solver;
    //gsSparseSolver<>::BiCGSTABDiagonal solver;
    //gsSparseSolver<>::QR  solver;

protected:

    size_t m_dofs;
    T m_dt;
    T m_t;
    bool m_NL;

    bool m_verbose;

    bool m_quasiNewton;

    bool m_converged;

    int m_maxIterations;
    int m_numIterations;

    T m_tolerance;


protected:

    gsSparseMatrix<T> m_mass;
    gsSparseMatrix<T> m_damp;
    gsSparseMatrix<T> m_stif;
    std::function < gsSparseMatrix<T> ( gsMatrix<T> const & ) > m_jacobian;
    std::function < gsMatrix<T> ( gsMatrix<T> const &, T ) > m_residualFun;
    std::function < gsMatrix<T> ( T) > m_forceFun;
    gsSparseMatrix<T> m_jacMat;
    gsMatrix<T> m_resVec;

    gsMatrix<T> m_massInv;
    gsMatrix<T> m_forceVec;

    gsSparseMatrix<T> m_sysmat;
    gsSparseMatrix<T> m_syseye;
    typename gsBlockOp<T>::Ptr m_sysmatBlock;
    gsMatrix<T> m_sysvec;

    // For Euler-type methods
    gsMatrix<T> m_sol;
    gsMatrix<T> m_solOld;
    gsMatrix<T> m_dsol;

    // For Newmark-type of methods
    gsMatrix<T> m_uNew;
    gsMatrix<T> m_vNew;
    gsMatrix<T> m_aNew;

    // gsFunction<T> m_dini;
    // gsFunction<T> m_vini;

    gsMultiPatch<T> m_solution;

    std::string m_method;

    void initializeCoefficients();
    void initializeSolution();

    T m_updateNorm;
    T m_residue;

    // void constructSolution();
    gsMatrix<T> m_displacements;
    gsMatrix<T> m_velocities;
    gsMatrix<T> m_accelerations;

    void constructSystemExplEuler();
    void constructSystemImplEuler();
    void constructSystemBlock();
    void initiate();
    void stepExplEuler();
    void stepExplEulerNL();
    void stepImplEuler();
    void stepImplEulerNL();
    void stepNewmark();
    void stepNewmarkNL();
    void stepBathe();
    void stepBatheNL();
};


} // namespace gismo

//////////////////////////////////////////////////
//////////////////////////////////////////////////

namespace gismo
{

    template <class T>
    void gsTimeIntegrator<T>::initializeCoefficients()
    {
        m_verbose = false;
        m_quasiNewton = false;
        m_numIterations = 0;
        m_maxIterations = 100;
        m_tolerance = 1e-12;
        m_converged = false;
        m_t = 0;
        m_residue = 1.0;
        m_updateNorm = 1.0;
    }

    template <class T>
    void gsTimeIntegrator<T>::initializeSolution()
    {
        // Defines homogenous solutions as initial solution. Can be overwritten with initialDisplacement() or initialVelocity() later
        if (m_method=="ExplEuler" || m_method=="ImplEuler")
        {
            m_sol = gsMatrix<T>::Zero(2*m_dofs,1);
        }
        else if ((m_method=="Newmark") || (m_method=="Bathe"))
        {
            m_uNew = gsMatrix<T>::Zero(m_dofs,1);
            m_vNew = gsMatrix<T>::Zero(m_dofs,1);
            m_aNew = gsMatrix<T>::Zero(m_dofs,1);
        }
        else
        {
            gsInfo<<"Method unknown...";
        }
    }

    template <class T>
    void gsTimeIntegrator<T>::setDisplacement(gsMatrix<T>& displ)
    {
        if (m_method=="ExplEuler" || m_method=="ImplEuler")
            m_sol.block(0,0,m_dofs,1) = displ;
        else if ((m_method=="Newmark") || (m_method=="Bathe"))
            m_uNew = displ;
    }

    template <class T>
    void gsTimeIntegrator<T>::setVelocity(gsMatrix<T>& velo)
    {
        if (m_method=="ExplEuler" || m_method=="ImplEuler")
            m_sol.block(m_dofs,0,m_dofs,1) = velo;
        else if ((m_method=="Newmark") || (m_method=="Bathe"))
            m_vNew = velo;
    }

    template <class T>
    void gsTimeIntegrator<T>::setAcceleration(gsMatrix<T>& accel)
    {
        if (m_method=="ExplEuler" || m_method=="ImplEuler")
            gsWarn<<"Initial acceleration not implemented for method "<<m_method<<"\n";
        else if ((m_method=="Newmark") || (m_method=="Bathe"))
            m_aNew = accel;
    }

    template <class T>
    void gsTimeIntegrator<T>::resetIteration()
    {
        m_numIterations = 0;
        m_residue = 1.0;
        m_updateNorm = 1.0;
    }

    template <class T>
    void gsTimeIntegrator<T>::constructSystemExplEuler()
    {
            m_sysmat = gsSparseMatrix<T>(2*m_dofs,2*m_dofs);
            m_sysmat.setZero();

            m_sysvec.setZero(2*m_dofs,1);

            m_massInv = m_mass.toDense().inverse();
            // Explicit Euler Method
            gsSparseMatrix<T>Imat(m_dofs,m_dofs);
            Imat.setIdentity();

            gsMatrix<T> sub1 = -m_massInv*m_stif;
            gsMatrix<T> sub2 = -m_massInv*m_damp;

            // top-right
            for (index_t i=0; i < m_dofs; i++)
                m_sysmat.insert(i,i+m_dofs) = 1.0;
            // bottom-left
            for (index_t i=0; i < m_dofs; i++)
            {
              for (index_t j=0; j < m_dofs; j++)
              {
                m_sysmat.insert(i+m_dofs,j) = sub1(i,j);
              }
            }
            // bottom-right
            for (index_t i=0; i < m_dofs; i++)
            {
              for (index_t j=0; j < m_dofs; j++)
              {
                m_sysmat.insert(i+m_dofs,j+m_dofs) = sub2(i,j);
              }
            }

            m_sysvec.setZero(2*m_dofs,1);

            m_sysmat *= m_dt;
    }

    template <class T>
    void gsTimeIntegrator<T>::constructSystemImplEuler()
    {
            constructSystemExplEuler();

            gsSparseMatrix<T> eye(2*m_dofs,2*m_dofs);
            eye.setIdentity();
            m_sysmat = eye-m_sysmat;

            // gsSparseMatrix<T> eye(2*m_dofs,2*m_dofs);
            // for (index_t i=0; i < 2*m_dofs; i++)
            //     m_sysmat(i,i+m_dofs) = 1.0 - m_sysmat(i,i);

    }

    template <class T>
    void gsTimeIntegrator<T>::constructSystemBlock()
    {
            m_massInv = m_mass.toDense().inverse();

            m_sysmatBlock = gsBlockOp<>::make(2,2);
            m_sysvec = gsMatrix<T>::Zero(2*m_dofs,1);

            gsSparseMatrix<> eye(m_dofs,m_dofs);
            eye.setIdentity();
            // m_sysmatBlock->addOperator(0,0,gsIdentityOp<real_t>::make(m_dofs) );
            // m_sysmatBlock->addOperator(1,1,gsIdentityOp<real_t>::make(m_dofs) );
            // m_sysmatBlock->addOperator(0,0,makeMatrixOp(eye) );
            // m_sysmatBlock->addOperator(1,1,makeMatrixOp(eye) );
            m_sysmatBlock->addOperator(0,0,gsIdentityOp<real_t>::make(m_dofs) );
            m_sysmatBlock->addOperator(0,1,makeMatrixOp(-m_dt*eye) );
            m_sysmatBlock->addOperator(1,0,makeMatrixOp(m_dt*m_massInv*m_stif) );
            m_sysmatBlock->addOperator(1,1,makeMatrixOp(eye + m_dt*m_massInv*m_damp) );

            gsMatrix<T> uNew = m_sol.topRows(m_dofs);
            gsMatrix<T> uOld = m_solOld.topRows(m_dofs);
            gsMatrix<T> vNew = m_sol.bottomRows(m_dofs);
            gsMatrix<T> vOld = m_solOld.bottomRows(m_dofs);

            m_sysvec.topRows(m_dofs) = uNew - uOld - m_dt*vNew;
            m_sysvec.bottomRows(m_dofs) = vNew - vOld + m_dt*(m_massInv*(-m_resVec));
    }

    // template <class T>
    // void gsTimeIntegrator<T>::initiate(const gsBoundaryConditions<T> BCs, );
    // {
    //     m_forcing.bottomRows(m_dofs) = force;
    //     gsFunctionExpr<> initialShape(dinix,diniy,diniz,3);
    //
    //     // Construct assembler object for initial solution. Hence, thickness and
    //     // density are chosen such that M-matrix corresponds with the one needed for
    //     // the L2-projection
    //     real_t unitThickness = 0.5; // HALF thickness
    //     real_t unitDensity = 1;
    //     gsShellAssembler<real_t> assemblerInitial(mp, 1.0, 0.5, 1.0,
    //                                        0.0, BCs, initialShape, clamped, pLoads);
    //
    //     // TO DO: Check of ik ook gewoon een expression assembler kan gebruiken.
    // }

    // template <class T>
    // void gsTimeIntegrator<T>::solution();
    // {
    //   assembler.constructSolution(m_sol,m_solution);
    // }

    template <class T>
    void gsTimeIntegrator<T>::setMethod(std::string string)
    {
        if ( (string != "ExplEuler") &&
             (string != "ImplEuler") &&
             (string != "Newmark") &&
             (string != "Bathe") )
            {
                gsInfo<<"Method Unknown... Process terminating\n";
            }
        else
        {
            m_method = string;
            initializeSolution();
        }
    }

    template <class T>
    void gsTimeIntegrator<T>::initiate()
    {
        // BUILD SYSTEM
        if (m_method=="ExplEuler")
        {
            initializeSolution();
            if (!m_NL)
                constructSystemExplEuler();
        }
        else if (m_method=="ImplEuler")
        {
            initializeSolution();
            if (!m_NL)
                constructSystemImplEuler();
        }
        else if (m_method=="Newmark")
        {
            initializeSolution();
        }
        else if (m_method=="Bathe")
        {
            initializeSolution();
        }
        else
        {
            gsInfo<<"Method unknown...";
        }

      // case "ImplEulerBlockOp":
      // {

      // }
    }

    template <class T>
    void gsTimeIntegrator<T>::stepExplEuler()
    {
        // ONLY FOR LINEAR SYSTEM!
        m_forceVec = m_forceFun(m_t);
        m_sysvec.bottomRows(m_dofs) = m_dt*m_massInv*m_forceVec;

        m_sol += m_sysmat * m_sol + m_sysvec;
        m_t += m_dt;
    }

    template <class T>
    void gsTimeIntegrator<T>::stepExplEulerNL()
    {
        m_solOld = m_sol;
        gsMatrix<T> uOld = m_solOld.topRows(m_dofs);
        gsMatrix<T> vOld = m_solOld.bottomRows(m_dofs);

        m_resVec = m_residualFun(uOld,m_t);

        m_sysmat = gsSparseMatrix<T>(2*m_dofs,2*m_dofs);
        m_sysmat.setZero();
        m_sysvec.setZero(2*m_dofs,1);

        // Generate system matrix [I , 0; 0 , M]
        // top-left
        for (index_t i=0; i < m_dofs; i++)
            m_sysmat.insert(i,i) = 1.0;
        // bottom-right
        for (index_t i=0; i < m_dofs; i++)
        {
          for (index_t j=0; j < m_dofs; j++)
          {
            m_sysmat.insert(i+m_dofs,j+m_dofs) = m_mass(i,j);
          }
        }

        m_sysvec.topRows(m_dofs) = vOld;
        m_sysvec.bottomRows(m_dofs) = -m_damp*vOld - m_resVec;

        m_solver.compute(m_sysmat);
        m_dsol = m_solver.solve(m_dt*m_sysvec);

        gsInfo<<m_dsol.topRows(m_dofs).transpose()<<"\n";

        m_sol += m_dsol;

        if (m_verbose)
        {
            gsInfo  <<"Explicit Euler Iteration Finished."
                    << m_dsol.norm()
                    <<"\n";
        }

        m_t += m_dt;
    }

    template <class T>
    void gsTimeIntegrator<T>::stepImplEuler()
    {
        m_forceVec = m_forceFun(m_t);
        m_sysvec.bottomRows(m_dofs) = m_dt*m_massInv*m_forceVec;

        m_solver.compute(m_sysmat);
        m_sol = m_solver.solve(m_sysvec+m_sol);
        m_t += m_dt;
    }

    template <class T>
    void gsTimeIntegrator<T>::stepImplEulerNL()
    {
        gsMatrix<T> solOld = m_sol;
        gsMatrix<T> uOld = solOld.topRows(m_dofs);
        gsMatrix<T> vOld = solOld.bottomRows(m_dofs);
        gsMatrix<T> uNew,vNew;

        m_sysmat = gsSparseMatrix<T>(2*m_dofs,2*m_dofs);
        m_sysmat.setZero();

        m_sysvec.setZero(2*m_dofs,1);

        gsSparseMatrix<T> eyeN(m_dofs,m_dofs);
        eyeN.setIdentity();
        gsSparseMatrix<T> eye2N(2*m_dofs,2*m_dofs);
        eye2N.setIdentity();

        while ( (m_updateNorm>m_tolerance || m_residue>m_tolerance) && m_numIterations<=m_maxIterations )
        {
            uNew = m_sol.topRows(m_dofs);
            vNew = m_sol.bottomRows(m_dofs);

            if ( (!m_quasiNewton) || (m_numIterations==0) )
                m_jacMat = m_jacobian(uNew);

            m_resVec = m_residualFun(uNew,m_t);

            m_sysvec.topRows(m_dofs) = uNew - uOld - m_dt*vNew;
            m_sysvec.bottomRows(m_dofs) = m_mass*(vNew - vOld) + m_dt*((-m_resVec));

            gsMatrix<T> sub1 = m_dt*m_jacMat;
            gsMatrix<T> sub2 = m_dt*m_damp + m_mass*eyeN;

            m_sysmat.setZero();
            // top-left
            for (index_t i=0; i < m_dofs; i++)
                m_sysmat.insert(i,i) = 1.0;
            // top-right
            for (index_t i=0; i < m_dofs; i++)
                m_sysmat.insert(i,i+m_dofs) = -m_dt*1.0;
            // bottom-left
            for (index_t i=0; i < m_dofs; i++)
            {
              for (index_t j=0; j < m_dofs; j++)
              {
                m_sysmat.insert(i+m_dofs,j) = sub1(i,j);
              }
            }
            // bottom-right
            for (index_t i=0; i < m_dofs; i++)
            {
              for (index_t j=0; j < m_dofs; j++)
              {
                m_sysmat.insert(i+m_dofs,j+m_dofs) = sub2(i,j);
              }
            }

            m_residue= m_sysvec.norm();
            m_solver.compute(m_sysmat);
            m_dsol = m_solver.solve(-m_sysvec);

            m_sol += m_dsol;
            m_updateNorm = m_dsol.norm();

            if (m_verbose)
            {
                gsInfo  <<"\tIteration: "   << m_numIterations
                        <<", residue: "     << m_residue
                        <<", update norm: " << m_updateNorm
                        <<"\n";
            }
            m_numIterations++;
        }

        m_converged = true;
        if ((m_numIterations==m_maxIterations+1) && m_verbose)
        {
            gsInfo<<"\tWARNING: Maximum number of iterations reached"<<"\n";
            m_converged = false;
        }

        m_t += m_dt;

        resetIteration();
    }

    // template <class T>
    // void gsTimeIntegrator<T>::stepImplEulerNLOp()
    // {
    //     m_solOld = m_sol;
    //     m_dsol = gsMatrix<T>::Zero(2*m_dofs,1);

    //     gsBlockOp<>::Ptr Amat=gsBlockOp<>::make(2,2);

    //     gsSparseMatrix<T> eye(m_dofs,m_dofs);
    //     eye.setIdentity();

    //     // m_sysmatBlock=gsBlockOp<>::make(2,2);
    //     m_sysvec = gsMatrix<T>::Zero(2*m_dofs,1);

    //     gsGMRes<> gmres(Amat);
    //     gmres.setMaxIterations(100);
    //     gmres.setTolerance(1e-8);

    //     while ( (m_updateNorm>m_tolerance || m_residue>m_tolerance) && m_numIterations<=m_maxIterations )
    //     {
    //         m_uNew = m_sol.topRows(m_dofs);

    //         m_jacMat = m_jacobian(m_uNew);
    //         m_resVec = m_residualFun(m_uNew,m_t);

    //         m_sysvec = gsMatrix<T>::Zero(2*m_dofs,1);

    //         gsSparseMatrix<> eye(m_dofs,m_dofs);
    //         eye.setIdentity();
    //         // m_sysmatBlock->addOperator(0,0,gsIdentityOp<real_t>::make(m_dofs) );
    //         // m_sysmatBlock->addOperator(1,1,gsIdentityOp<real_t>::make(m_dofs) );
    //         // m_sysmatBlock->addOperator(0,0,makeMatrixOp(eye) );
    //         // m_sysmatBlock->addOperator(1,1,makeMatrixOp(eye) );
    //         Amat->addOperator(0,0,gsIdentityOp<real_t>::make(m_dofs) );
    //         Amat->addOperator(0,1,makeMatrixOp(-m_dt*eye) );
    //         Amat->addOperator(1,0,makeMatrixOp(m_dt*m_massInv*m_stif) );
    //         Amat->addOperator(1,1,makeMatrixOp(eye + m_dt*m_massInv*m_damp) );

    //         gsMatrix<T> uNew = m_sol.topRows(m_dofs);
    //         gsMatrix<T> uOld = m_solOld.topRows(m_dofs);
    //         gsMatrix<T> vNew = m_sol.bottomRows(m_dofs);
    //         gsMatrix<T> vOld = m_solOld.bottomRows(m_dofs);

    //         m_sysvec.topRows(m_dofs) = uNew - uOld - m_dt*vNew;
    //         m_sysvec.bottomRows(m_dofs) = vNew - vOld + m_dt*(m_massInv*(-m_resVec));


    //         m_residue= m_resVec.norm();
    //         gmres.solve(-m_sysvec,m_dsol);

    //         gsInfo<<m_sysvec.transpose()<<"\n";
    //         gsInfo<<m_dsol.transpose()<<"\n";

    //         // m_sol += m_dsol;
    //         // m_updateNorm = m_dsol.norm();

    //         // if (m_verbose)
    //         // {
    //         //     gsInfo  <<"\tIteration: "   << m_numIterations
    //         //             <<", residue: "     << m_residue
    //         //             <<", update norm: " << m_updateNorm
    //         //             <<"\n";
    //         // }
    //         m_numIterations++;
    //     }

    //     m_converged = true;
    //     if ((m_numIterations==m_maxIterations+1) && m_verbose)
    //     {
    //         gsInfo<<"\tWARNING: Maximum number of iterations reached"<<"\n";
    //         m_converged = false;
    //     }

    //     m_uNew = m_sol.topRows(m_dofs);
    //     m_vNew = m_sol.bottomRows(m_dofs);
    //     m_t += m_dt;
    // }



    template <class T>
    void gsTimeIntegrator<T>::stepNewmark()
    {
        T gamma = 1;
        m_t += gamma*m_dt;

        gsMatrix<T> uOld,vOld,aOld;
        uOld = m_uNew;
        vOld = m_vNew;
        aOld = m_aNew;

        m_forceVec = m_forceFun(m_t);

        gsSparseMatrix<> lhs = m_stif + 4/(math::pow(gamma*m_dt,2))*m_mass + 2/(gamma*m_dt)*m_damp;
        gsMatrix<> rhs = m_forceVec + m_mass*(4/(math::pow(gamma*m_dt,2))*(uOld)+4/(gamma*m_dt)*vOld+aOld) + m_damp*(2/(gamma*m_dt)*(uOld)+vOld);

        m_solver.compute(lhs);
        m_sol = m_solver.solve(rhs);

        m_uNew = m_sol;

        m_vNew = 2/(gamma*m_dt)*(m_uNew-uOld) - vOld;
        m_aNew = (m_uNew-uOld-vOld*m_dt)*4/(math::pow(gamma*m_dt,2)) - aOld;
    }

    template <class T>
    void gsTimeIntegrator<T>::stepNewmarkNL()
    {
        T gamma = 1;
        m_t += gamma*m_dt;

        gsMatrix<T> uOld,vOld,aOld;
        uOld = m_uNew;
        vOld = m_vNew;
        aOld = m_aNew;

        while ( (m_updateNorm>m_tolerance || m_residue>m_tolerance) && m_numIterations<=m_maxIterations )
        {
            if ( (!m_quasiNewton) || (m_numIterations==0) )
                m_jacMat = m_jacobian(m_uNew);

            m_resVec = m_residualFun(m_uNew,m_t);

            gsSparseMatrix<> lhs = m_jacMat + 4/(math::pow(gamma*m_dt,2))*m_mass + 2/(gamma*m_dt)*m_damp;
            gsMatrix<> rhs = m_resVec - m_mass*(4/(math::pow(gamma*m_dt,2))*(m_uNew-uOld)-4/(gamma*m_dt)*vOld-aOld) - m_damp*(2/(gamma*m_dt)*(m_uNew-uOld)-vOld);
            m_solver.compute(lhs);
            m_dsol = m_solver.solve(rhs);
            m_updateNorm = m_dsol.norm();
            m_uNew += m_dsol;
            m_residue = rhs.norm();

            if (m_verbose)
            {
                gsInfo  <<"\tIteration: "   << m_numIterations
                        <<", residue: "     << m_residue
                        <<", update norm: " << m_updateNorm
                        <<"\n";
            }
            m_numIterations++;
        }
        // Construct velocities and accelerations
        m_vNew = 2/(gamma*m_dt)*(m_uNew-uOld) - vOld;
        m_aNew = (m_uNew-uOld-vOld*m_dt)*4/(math::pow(gamma*m_dt,2)) - aOld;

        m_converged = true;
        if ((m_numIterations==m_maxIterations+1) && m_verbose)
        {
            gsInfo<<"\tWARNING: Maximum number of iterations reached"<<"\n";
            m_converged = false;
        }

        resetIteration();
    }

    template <class T>
    void gsTimeIntegrator<T>::stepBathe()
    {
        T gamma = 0.5;
        m_t += gamma*m_dt;

        gsMatrix<T> uOld,vOld,aOld;
        gsMatrix<T> uStep,vStep,aStep;
        uOld = m_uNew;
        vOld = m_vNew;
        aOld = m_aNew;

        m_forceVec = m_forceFun(m_t);

        gsSparseMatrix<> lhs = m_stif + 4/(math::pow(gamma*m_dt,2))*m_mass + 2/(gamma*m_dt)*m_damp;
        gsMatrix<> rhs = m_forceVec + m_mass*(4/(math::pow(gamma*m_dt,2))*(uOld)+4/(gamma*m_dt)*vOld+aOld) + m_damp*(2/(gamma*m_dt)*(uOld)+vOld);;
        m_solver.compute(lhs);
        m_sol = m_solver.solve(rhs);
        m_uNew = m_sol;

        m_vNew = 2/(gamma*m_dt)*(m_uNew-uOld) - vOld;
        m_aNew = (m_uNew-uOld-vOld*m_dt)*4/(math::pow(gamma*m_dt,2)) - aOld;

        uStep = m_uNew;
        vStep = m_vNew;
        aStep = m_aNew;

        T c1 = (1-gamma)/(gamma*m_dt);
        T c2 = (-1)/((1-gamma)*gamma*m_dt);
        T c3 = (2-gamma)/((1-gamma)*m_dt);

        m_t += (1-gamma)*m_dt;

        lhs = m_stif + c3*c3*m_mass + c3*m_damp;
        rhs = m_forceVec - m_mass*(c1*vOld+c2*vStep+c1*c3*uOld+c3*c2*uStep) - m_damp*(c2*uStep+c1*uOld);
        m_solver.compute(lhs);
        m_sol = m_solver.solve(rhs);
        m_uNew = m_sol;
        m_residue = rhs.norm();

        m_vNew = c1*uOld + c2*uStep + c3*m_uNew;
        m_aNew = c1*vOld + c2*vStep + c3*m_vNew;

    }

    template <class T>
    void gsTimeIntegrator<T>::stepBatheNL()
    {
        T gamma = 0.5;
        m_t += gamma*m_dt;

        gsMatrix<T> uOld,vOld,aOld;
        gsMatrix<T> uStep,vStep,aStep;
        uOld = m_uNew;
        vOld = m_vNew;
        aOld = m_aNew;

        if (m_verbose)
        {
            gsInfo  << "\t Stage 1:\n";
        }
        while ( (m_updateNorm>m_tolerance || m_residue>m_tolerance) && m_numIterations<=m_maxIterations )
        {
            if ( (!m_quasiNewton) || (m_numIterations==0) )
                m_jacMat = m_jacobian(m_uNew);

            m_resVec = m_residualFun(m_uNew,m_t);

            gsSparseMatrix<> lhs = m_jacMat + 4/(math::pow(gamma*m_dt,2))*m_mass + 2/(gamma*m_dt)*m_damp;
            gsMatrix<> rhs = m_resVec - m_mass*(4/(math::pow(gamma*m_dt,2))*(m_uNew-uOld)-4/(gamma*m_dt)*vOld-aOld) - m_damp*(2/(gamma*m_dt)*(m_uNew-uOld)-vOld);
            m_solver.compute(lhs);
            m_dsol = m_solver.solve(rhs);
            m_updateNorm = m_dsol.norm();
            m_uNew += m_dsol;
            m_residue = rhs.norm();

            if (m_verbose)
            {
                gsInfo  <<"\tIteration: "   << m_numIterations
                        <<", residue: "     << m_residue
                        <<", update norm: " << m_updateNorm
                        <<"\n";
            }
            m_numIterations++;
        }
        // Construct velocities and accelerations
        m_vNew = 2/(gamma*m_dt)*(m_uNew-uOld) - vOld;
        m_aNew = (m_uNew-uOld-vOld*m_dt)*4/(math::pow(gamma*m_dt,2)) - aOld;
        uStep = m_uNew;
        vStep = m_vNew;
        aStep = m_aNew;

        T c1 = (1-gamma)/(gamma*m_dt);
        T c2 = (-1)/((1-gamma)*gamma*m_dt);
        T c3 = (2-gamma)/((1-gamma)*m_dt);

        if (m_verbose)
        {
            gsInfo  << "\tStage 2:\n";
        }

        m_t += (1-gamma)*m_dt;
        resetIteration();

        while ( (m_updateNorm>1e-6 || m_residue>1e-6) && m_numIterations<=m_maxIterations )
        {
            if ( (!m_quasiNewton) || (m_numIterations==0) )
                m_jacMat = m_jacobian(m_uNew);

            m_resVec = m_residualFun(m_uNew,m_t);

            gsSparseMatrix<T> lhs = m_jacMat + c3*c3*m_mass + c3*m_damp;
            gsMatrix<T> rhs = m_resVec - m_mass*(c1*vOld+c2*vStep+c1*c3*uOld+c3*c2*uStep+c3*c3*m_uNew) - m_damp*(c1*uOld+c2*uStep+c3*m_uNew);
            m_solver.compute(lhs);
            m_dsol = m_solver.solve(rhs);
            m_updateNorm = m_dsol.norm();
            m_uNew += m_dsol;
            m_residue = rhs.norm();

            if (m_verbose)
            {
                gsInfo  <<"\tIteration: "   << m_numIterations
                        <<", residue: "     << m_residue
                        <<", update norm: " << m_updateNorm
                        <<"\n";
            }
            m_numIterations++;
        }

        m_vNew = c1*uOld + c2*uStep + c3*m_uNew;
        m_aNew = c1*vOld + c2*vStep + c3*m_vNew;

        m_converged = true;
        if ((m_numIterations==m_maxIterations+1) && m_verbose)
        {
            gsInfo<<"\tWARNING: Maximum number of iterations reached"<<"\n";
            m_converged = false;
        }

        resetIteration();
    }

    template <class T>
    void gsTimeIntegrator<T>::step()
    {
        if (m_verbose)
            gsInfo  << "time = "<<m_t<<"\n";

        // The actual step
        if (m_method=="ExplEuler")
        {
            if (m_NL)
                stepExplEulerNL();
            else
                stepExplEuler();
        }
        else if (m_method=="ImplEuler")
        {
            if (m_NL)
                stepImplEulerNL();
            else
                stepImplEuler();
        }
        else if (m_method=="Newmark")
        {
            if (m_NL)
                stepNewmarkNL();
            else
                stepNewmark();
        }
        else if (m_method=="Bathe")
        {
            if (m_NL)
                stepBatheNL();
            else
                stepBathe();
        }
        else
        {
            gsInfo<<"Method unknown, terminating..\n";
        }
        m_numIterations = 0;
    }

    template <class T>
    void gsTimeIntegrator<T>::constructSolution()
    {
        if ((m_method=="ExplEuler") || (m_method=="ImplEuler"))
        {
            m_displacements = m_sol.topRows(m_dofs);
            m_velocities = m_sol.bottomRows(m_dofs);
        }
        else if ((m_method=="Newmark") || (m_method=="Bathe"))
        {
            m_displacements = m_uNew;
            m_velocities = m_vNew;
            m_accelerations = m_aNew;
        }
    }

} // namespace gismo
