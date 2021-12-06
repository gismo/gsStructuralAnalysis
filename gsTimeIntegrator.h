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

    void initializeCoefficients();
    void initializeSolution();
    void constructSystemExplEuler();
    void constructSystemImplEuler();
    void constructSystemBlock();
    void initiate();
    void stepExplEuler();
    void stepExplEulerNL();
    void stepImplEuler();
    void stepImplEulerNL();
    void stepImplEulerNLOp();
    void stepNewmark();
    void stepNewmarkNL();
    void stepBathe();
    void stepBatheNL();
    void stepCentralDiff();
    void stepCentralDiffNL();
    void stepRK4();
    void stepRK4NL();

protected:

    /// Linear solver employed
    gsSparseSolver<>::LU  m_solver;
    //gsSparseSolver<>::BiCGSTABDiagonal solver;
    //gsSparseSolver<>::QR  solver;

protected:

    gsSparseMatrix<T> m_mass;
    gsSparseMatrix<T> m_damp;
    gsSparseMatrix<T> m_stif;
    std::function < gsSparseMatrix<T> ( gsMatrix<T> const & ) > m_jacobian;
    std::function < gsMatrix<T> ( gsMatrix<T> const &, T ) > m_residualFun;
    std::function < gsMatrix<T> ( T) > m_forceFun;
    gsSparseMatrix<T> m_jacMat;
    gsMatrix<T> m_resVec;

    size_t m_dofs;
    T m_dt;
    T m_t;
    bool m_NL;

    bool m_verbose;

    bool m_quasiNewton;

    bool m_converged;

    int m_maxIterations;
    int m_numIterations;

    bool m_first;

    T m_tolerance;

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

    T m_updateNorm;
    T m_residue;

    // void constructSolution();
    gsMatrix<T> m_displacements;
    gsMatrix<T> m_velocities;
    gsMatrix<T> m_accelerations;

};


} // namespace gismo

//////////////////////////////////////////////////
//////////////////////////////////////////////////



#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsTimeIntegrator.hpp)
#endif