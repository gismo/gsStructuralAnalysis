/** @file gsTimeIntegrator.hpp

    @brief Provides temporal solvers for structural analysis problems

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#pragma once

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
        m_first = true;
    }

    template <class T>
    void gsTimeIntegrator<T>::initializeSolution()
    {
        // Defines homogenous solutions as initial solution. Can be overwritten with initialDisplacement() or initialVelocity() later
        if (m_method=="ExplEuler" || m_method=="ImplEuler" || m_method=="RK4")
        {
            m_sol = gsMatrix<T>::Zero(2*m_dofs,1);
        }
        else if ((m_method=="Newmark") || (m_method=="Bathe") || (m_method=="CentralDiff"))
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
        if (m_method=="ExplEuler" || m_method=="ImplEuler" || m_method=="RK4")
            m_sol.block(0,0,m_dofs,1) = displ;
        else if ((m_method=="Newmark") || (m_method=="Bathe") || (m_method=="CentralDiff"))
            m_uNew = displ;
    }

    template <class T>
    void gsTimeIntegrator<T>::setVelocity(gsMatrix<T>& velo)
    {
        if (m_method=="ExplEuler" || m_method=="ImplEuler" || m_method=="RK4")
            m_sol.block(m_dofs,0,m_dofs,1) = velo;
        else if ((m_method=="Newmark") || (m_method=="Bathe") || (m_method=="CentralDiff"))
            m_vNew = velo;
    }

    template <class T>
    void gsTimeIntegrator<T>::setAcceleration(gsMatrix<T>& accel)
    {
        if (m_method=="ExplEuler" || m_method=="ImplEuler" || m_method=="RK4")
            gsWarn<<"Initial acceleration not implemented for method "<<m_method<<"\n";
        else if ((m_method=="Newmark") || (m_method=="Bathe") || (m_method=="CentralDiff"))
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
            for (size_t i=0; i < m_dofs; i++)
                m_sysmat.insert(i,i+m_dofs) = 1.0;
            // bottom-left
            for (size_t i=0; i < m_dofs; i++)
            {
              for (size_t j=0; j < m_dofs; j++)
              {
                m_sysmat.insert(i+m_dofs,j) = sub1(i,j);
              }
            }
            // bottom-right
            for (size_t i=0; i < m_dofs; i++)
            {
              for (size_t j=0; j < m_dofs; j++)
              {
                m_sysmat.insert(i+m_dofs,j+m_dofs) = sub2(i,j);
              }
            }

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
             (string != "Bathe") &&
             (string !="CentralDiff") &&
             (string !="RK4")            )
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
        else if (m_method=="CentralDiff")
        {
            initializeSolution();
        }
        else if (m_method=="RK4")
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
        m_solOld = m_sol;
        gsMatrix<T> uOld = m_solOld.topRows(m_dofs);
        gsMatrix<T> vOld = m_solOld.bottomRows(m_dofs);
        m_forceVec = m_forceFun(m_t);
        m_massInv = m_mass.toDense().inverse();

        m_sol.topRows(m_dofs) += m_dt * vOld;
        m_sol.bottomRows(m_dofs) += m_dt * m_massInv * ( m_forceVec - m_stif * uOld - m_damp * vOld);

        m_t += m_dt;
    }

    template <class T>
    void gsTimeIntegrator<T>::stepExplEulerNL()
    {
        m_solOld = m_sol;
        gsMatrix<T> uOld = m_solOld.topRows(m_dofs);
        gsMatrix<T> vOld = m_solOld.bottomRows(m_dofs);
        m_resVec = m_residualFun(uOld,m_t);
        m_massInv = m_mass.toDense().inverse();

        m_sol.topRows(m_dofs) += m_dt * vOld;
        m_sol.bottomRows(m_dofs) += m_dt * m_massInv * ( m_resVec - m_damp * vOld);

        m_t += m_dt;
    }

    template <class T>
    void gsTimeIntegrator<T>::stepImplEuler()
    {
        this->constructSystemImplEuler();
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
            m_sysvec.bottomRows(m_dofs) = m_mass*(vNew - vOld) + m_dt * m_damp * vNew + m_dt*((-m_resVec));

            gsMatrix<T> sub1 = m_dt*m_jacMat;
            gsMatrix<T> sub2 = m_dt*m_damp + m_mass*eyeN;

            m_sysmat.setZero();
            // top-left
            for (size_t i=0; i < m_dofs; i++)
                m_sysmat.insert(i,i) = 1.0;
            // top-right
            for (size_t i=0; i < m_dofs; i++)
                m_sysmat.insert(i,i+m_dofs) = -m_dt*1.0;
            // bottom-left
            for (size_t i=0; i < m_dofs; i++)
            {
              for (size_t j=0; j < m_dofs; j++)
              {
                m_sysmat.insert(i+m_dofs,j) = sub1(i,j);
              }
            }
            // bottom-right
            for (size_t i=0; i < m_dofs; i++)
            {
              for (size_t j=0; j < m_dofs; j++)
              {
                m_sysmat.insert(i+m_dofs,j+m_dofs) = sub2(i,j);
              }
            }

            m_residue = m_sysvec.norm();
            m_solver.compute(m_sysmat);
            m_dsol = m_solver.solve(-m_sysvec);

            m_sol += m_dsol;
            m_updateNorm = m_dsol.norm() / m_sol.norm();

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

    template <class T>
    void gsTimeIntegrator<T>::stepImplEulerNLOp()
    {
        gsMatrix<T> solOld = m_sol;
        gsMatrix<T> uOld = solOld.topRows(m_dofs);
        gsMatrix<T> vOld = solOld.bottomRows(m_dofs);
        gsMatrix<T> uNew,vNew;

        gsBlockOp<>::Ptr Amat=gsBlockOp<>::make(2,2);

        gsSparseMatrix<T> eye(m_dofs,m_dofs);
        eye.setIdentity();

        // m_sysmatBlock=gsBlockOp<>::make(2,2);
        m_sysvec = gsMatrix<T>::Zero(2*m_dofs,1);

        gsGMRes<> gmres(Amat);
        gmres.setMaxIterations(100);
        gmres.setTolerance(1e-8);
        while ( (m_updateNorm>m_tolerance || m_residue>m_tolerance) && m_numIterations<=m_maxIterations )
        {
            uNew = m_sol.topRows(m_dofs);
            vNew = m_sol.bottomRows(m_dofs);

            if ( (!m_quasiNewton) || (m_numIterations==0) )
                m_jacMat = m_jacobian(uNew);

            m_resVec = m_residualFun(uNew,m_t);

            m_sysvec.topRows(m_dofs) = uNew - uOld - m_dt*vNew;
            m_sysvec.bottomRows(m_dofs) = m_mass*(vNew - vOld) + m_dt * m_damp * vNew + m_dt*(-m_resVec);

            m_jacMat = m_jacobian(uNew);

            Amat->addOperator(0,0,gsIdentityOp<T>::make(m_dofs) );
            Amat->addOperator(0,1,makeMatrixOp(-m_dt*eye) );
            Amat->addOperator(1,0,makeMatrixOp(m_dt*m_jacMat) );
            Amat->addOperator(1,1,makeMatrixOp(m_mass + m_dt*m_damp) );

            m_residue= m_sysvec.norm();
            gmres.solve(-m_sysvec,m_dsol);

            m_sol += m_dsol;
            m_updateNorm = m_dsol.norm() / m_sol.norm();

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


        gsMatrix<> rhs;
        gsSparseMatrix<> lhs;

        m_resVec = m_residualFun(m_uNew,m_t);
        rhs = m_resVec - m_mass*(4/(math::pow(gamma*m_dt,2))*(m_uNew-uOld)-4/(gamma*m_dt)*vOld-aOld) - m_damp*(2/(gamma*m_dt)*(m_uNew-uOld)-vOld);
        T res0 = rhs.norm();

        while ( (m_updateNorm>m_tolerance || m_residue/res0>m_tolerance) && m_numIterations<=m_maxIterations )
        {
            if ( (!m_quasiNewton) || (m_numIterations==0) )
                m_jacMat = m_jacobian(m_uNew);

            lhs = m_jacMat + 4/(math::pow(gamma*m_dt,2))*m_mass + 2/(gamma*m_dt)*m_damp;
            m_solver.compute(lhs);
            m_dsol = m_solver.solve(rhs);
            m_updateNorm = m_dsol.norm()/m_uNew.norm();
            m_uNew += m_dsol;

            m_resVec = m_residualFun(m_uNew,m_t);
            rhs = m_resVec - m_mass*(4/(math::pow(gamma*m_dt,2))*(m_uNew-uOld)-4/(gamma*m_dt)*vOld-aOld) - m_damp*(2/(gamma*m_dt)*(m_uNew-uOld)-vOld);
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
            gsInfo  << "\tStage 1:\n";
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
            m_uNew += m_dsol;
            m_updateNorm = m_dsol.norm() / m_uNew.norm();
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

        while ( (m_updateNorm>m_tolerance || m_residue>m_tolerance) && m_numIterations<=m_maxIterations )
        {
            if ( (!m_quasiNewton) || (m_numIterations==0) )
                m_jacMat = m_jacobian(m_uNew);

            m_resVec = m_residualFun(m_uNew,m_t);

            gsSparseMatrix<T> lhs = m_jacMat + c3*c3*m_mass + c3*m_damp;
            gsMatrix<T> rhs = m_resVec - m_mass*(c1*vOld+c2*vStep+c1*c3*uOld+c3*c2*uStep+c3*c3*m_uNew) - m_damp*(c1*uOld+c2*uStep+c3*m_uNew);
            m_solver.compute(lhs);
            m_dsol = m_solver.solve(rhs);
            m_uNew += m_dsol;
            m_updateNorm = m_dsol.norm() / m_uNew.norm();
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
    void gsTimeIntegrator<T>::stepCentralDiff()
    {
        gsMatrix<T> uOld,vOld,aOld;
        m_massInv = m_mass.toDense().inverse();

        uOld = m_uNew; // u_n
        vOld = m_vNew; // v_n+1/2
        aOld = m_aNew; // a_n
        m_forceVec = m_forceFun(m_t);

        if (m_first)
        {
            aOld = m_massInv * (m_forceVec - m_stif * uOld - m_damp * vOld);
            vOld = vOld + m_dt / 2 * aOld; // v_1/2 = v_0 + dt/2*a_0
            m_first = false;
        }

        m_uNew = uOld + m_dt * vOld;
        m_aNew = m_massInv * (m_forceVec - m_stif * m_uNew - m_damp * vOld);
        m_vNew = vOld + m_dt * m_aNew; // v_n+1/2 = v_n-1/2 + dt*a_n

        m_t+=m_dt;
    }

    template <class T>
    void gsTimeIntegrator<T>::stepCentralDiffNL()
    {
        gsMatrix<T> uOld,vOld,aOld;
        m_massInv = m_mass.toDense().inverse();

        uOld = m_uNew; // u_n
        vOld = m_vNew; // v_n+1/2
        aOld = m_aNew; // a_n

        if (m_first)
        {
            m_resVec = m_residualFun(uOld,m_t);
            aOld = m_massInv * (m_resVec - m_damp * vOld);
            vOld = vOld + m_dt / 2 * aOld; // v_1/2 = v_0 + dt/2*a_0
            m_first = false;
        }

        m_uNew = uOld + m_dt * vOld;
        m_t+=m_dt;
        m_resVec = m_residualFun(uOld,m_t);
        m_aNew = m_massInv * (m_resVec - m_damp * vOld);
        m_vNew = vOld + m_dt * m_aNew; // v_n+1/2 = v_n-1/2 + dt*a_n

    }

    template <class T>
    void gsTimeIntegrator<T>::stepRK4()
    {
        m_massInv = m_mass.toDense().inverse();
        gsVector<T> k1(2*m_dofs), k2(2*m_dofs), k3(2*m_dofs), k4(2*m_dofs);
        gsVector<T> uTmp, vTmp;
        m_solOld = m_sol;

        gsVector<T> uOld = m_solOld.topRows(m_dofs);
        gsVector<T> vOld = m_solOld.bottomRows(m_dofs);

        m_resVec = m_forceFun(m_t) - m_stif * uOld;
        k1.topRows(m_dofs) = vOld;
        k1.bottomRows(m_dofs) = m_massInv * ( m_resVec - m_damp * vOld);

        uTmp = uOld + m_dt/2. * k1.topRows(m_dofs);
        vTmp = vOld + m_dt/2. * k1.bottomRows(m_dofs);
        m_resVec = m_forceFun(m_t + m_dt/2.) - m_stif * uTmp;
        k2.topRows(m_dofs) = vTmp;
        k2.bottomRows(m_dofs) = m_massInv * ( m_resVec - m_damp * vTmp);

        uTmp = uOld + m_dt/2. * k2.topRows(m_dofs);
        vTmp = vOld + m_dt/2. * k2.bottomRows(m_dofs);
        m_resVec = m_forceFun(m_t + m_dt/2.) - m_stif * uTmp;
        k3.topRows(m_dofs) = vTmp;
        k3.bottomRows(m_dofs) = m_massInv * ( m_resVec - m_damp * vTmp);

        uTmp = uOld + m_dt * k3.topRows(m_dofs);
        vTmp = vOld + m_dt * k3.bottomRows(m_dofs);
        m_resVec = m_forceFun(m_t + m_dt) - m_stif * uTmp;
        k4.topRows(m_dofs) = vTmp;
        k4.bottomRows(m_dofs) = m_massInv * ( m_resVec - m_damp * vTmp);

        m_sol += 1./6. * m_dt *  ( k1 + 2.*k2 + 2.*k3 + k4 );
        m_t += m_dt;
    }

    template <class T>
    void gsTimeIntegrator<T>::stepRK4NL()
    {
        m_massInv = m_mass.toDense().inverse();
        gsVector<T> k1(2*m_dofs), k2(2*m_dofs), k3(2*m_dofs), k4(2*m_dofs);
        gsVector<T> uTmp, vTmp;
        m_solOld = m_sol;

        gsVector<T> uOld = m_solOld.topRows(m_dofs);
        gsVector<T> vOld = m_solOld.bottomRows(m_dofs);

        m_resVec = m_residualFun(uOld,m_t);
        k1.topRows(m_dofs) = vOld;
        k1.bottomRows(m_dofs) = m_massInv * ( m_resVec - m_damp * vOld);

        uTmp = uOld + m_dt/2. * k1.topRows(m_dofs);
        vTmp = vOld + m_dt/2. * k1.bottomRows(m_dofs);
        m_resVec = m_residualFun(uTmp,m_t + m_dt/2.);
        k2.topRows(m_dofs) = vTmp;
        k2.bottomRows(m_dofs) = m_massInv * ( m_resVec - m_damp * vTmp);

        uTmp = uOld + m_dt/2. * k2.topRows(m_dofs);
        vTmp = vOld + m_dt/2. * k2.bottomRows(m_dofs);
        m_resVec = m_residualFun(uTmp,m_t + m_dt/2.);
        k3.topRows(m_dofs) = vTmp;
        k3.bottomRows(m_dofs) = m_massInv * ( m_resVec - m_damp * vTmp);

        uTmp = uOld + m_dt * k3.topRows(m_dofs);
        vTmp = vOld + m_dt * k3.bottomRows(m_dofs);
        m_resVec = m_residualFun(uTmp,m_t + m_dt);
        k4.topRows(m_dofs) = vTmp;
        k4.bottomRows(m_dofs) = m_massInv * ( m_resVec - m_damp * vTmp);

        m_sol += 1./6. * m_dt * ( k1 + 2.*k2 + 2.*k3 + k4 );

        m_t += m_dt;
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
                stepImplEulerNLOp();
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
        else if (m_method=="CentralDiff")
        {
            if (m_NL)
                stepCentralDiffNL();
            else
                stepCentralDiff();
        }
        else if (m_method=="RK4")
        {
            if (m_NL)
                stepRK4NL();
            else
                stepRK4();
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
        if ((m_method=="ExplEuler") || (m_method=="ImplEuler") || (m_method=="RK4"))
        {
            m_displacements = m_sol.topRows(m_dofs);
            m_velocities = m_sol.bottomRows(m_dofs);
        }
        else if ((m_method=="Newmark") || (m_method=="Bathe") || (m_method=="CentralDiff"))
        {
            m_displacements = m_uNew;
            m_velocities = m_vNew;
            m_accelerations = m_aNew;
        }
    }

} // namespace gismo