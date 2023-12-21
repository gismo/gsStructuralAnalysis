/** @file gsTimeIntegrator.h

    @brief Provides temporal solvers for structural analysis problems

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/


#pragma once

#include <gsSolver/gsBlockOp.h>
#include <gsStructuralAnalysis/src/gsStructuralAnalysisTools/gsStructuralAnalysisTypes.h>

namespace gismo
{
/**
    @brief Provides temporal solvers for structural analysis problems

    \tparam T coefficient type

    \ingroup gsTimeIntegrator
*/
template <class T>
class gsTimeIntegrator
{

    typedef typename gsStructuralAnalysisOps<T>::TForce_t    TForce_t;
    typedef typename gsStructuralAnalysisOps<T>::TResidual_t TResidual_t;
    typedef typename gsStructuralAnalysisOps<T>::Jacobian_t  Jacobian_t;
    typedef typename gsStructuralAnalysisOps<T>::dJacobian_t dJacobian_t;

public:


    /**
     * @brief      Constructor
     *
     * @param[in]  Mass   The mass matrix
     * @param[in]  Damp   The damping matrix
     * @param[in]  Stif   The linear stiffness matrix
     * @param[in]  Force  The time-dependent force
     * @param[in]  dt     The time step
     */
    gsTimeIntegrator(   const gsSparseMatrix<T> &Mass,
                        const gsSparseMatrix<T> &Damp,
                        const gsSparseMatrix<T> &Stif,
                        const TForce_t &Force,
                        const T dt
                        )
    : m_mass(Mass),
      m_damp(Damp),
      m_stif(Stif),
      m_jacobian(nullptr),
      m_djacobian(nullptr),
      m_forceFun(Force),
      m_dt(dt)
    {
        m_dofs = m_mass.cols();
        m_method = "ExplEuler";
        m_NL = false;
        this->initializeCoefficients();
        this->initiate();
    }

    /**
     * @brief      Constructor
     *
     * @param[in]  Mass   The mass matrix
     * @param[in]  Stif   The linear stiffness matrix
     * @param[in]  Force  The time-dependent force
     * @param[in]  dt     The time step
     */
    gsTimeIntegrator(   const gsSparseMatrix<T> &Mass,
                        const gsSparseMatrix<T> &Stif,
                        const TForce_t &Force,
                        const T dt
                        )
    : m_mass(Mass),
      m_stif(Stif),
      m_jacobian(nullptr),
      m_djacobian(nullptr),
      m_forceFun(Force),
      m_dt(dt)
    {
        m_dofs = m_mass.cols();
        m_damp = gsSparseMatrix<T>(m_dofs,m_dofs);
        m_method = "ExplEuler";
        m_NL = false;
        this->initializeCoefficients();
        this->initiate();
        m_solver = gsSparseSolver<T>::get( "LU" ); // todo: make proper options
    }

private:

    /**
     * @brief      Constructor
     *
     * @param[in]  Mass      The mass matrix
     * @param[in]  Damp      The damping matrix
     * @param[in]  Residual  The residual vector
     * @param[in]  dt        The time step
     */
    gsTimeIntegrator(   const gsSparseMatrix<T> &Mass,
                        const gsSparseMatrix<T> &Damp,
                        const TResidual_t &Residual,
                        const T dt
                        )
    : m_mass(Mass),
      m_damp(Damp),
      m_jacobian(nullptr),
      m_residualFun(Residual),
      m_dt(dt)
    {
        m_dofs = m_mass.cols();
        m_method = "Newmark";
        m_NL = true;
        this->initializeCoefficients();
        this->initiate();
        m_solver = gsSparseSolver<T>::get( "LU" ); // todo: make proper options
    }

    /**
     * @brief      Constructor
     *
     * @param[in]  Mass      The mass matrix
     * @param[in]  Residual  The residual vector
     * @param[in]  dt        The time step
     */
    gsTimeIntegrator(   const gsSparseMatrix<T> &Mass,
                        const TResidual_t &Residual,
                        const T dt
                        )
    : m_mass(Mass),
      m_jacobian(nullptr),
      m_residualFun(Residual),
      m_dt(dt)
    {
        m_dofs = m_mass.cols();
        m_method = "Newmark";
        m_NL = true;
        m_damp = gsSparseMatrix<T>(m_dofs,m_dofs);
        this->initializeCoefficients();
        this->initiate();
        m_solver = gsSparseSolver<T>::get( "LU" ); // todo: make proper options
    }


public:

    /**
     * @brief      Constructor
     *
     * @param[in]  Mass      The mass matrix
     * @param[in]  Damp      The damping matrix
     * @param[in]  Jacobian  The Jacobian matrix
     * @param[in]  Residual  The residual vector
     * @param[in]  dt        The time step
     */
    gsTimeIntegrator(   const gsSparseMatrix<T> &Mass,
                        const gsSparseMatrix<T> &Damp,
                        const Jacobian_t  &Jacobian,
                        const TResidual_t &Residual,
                        const T dt
                        )
    : gsTimeIntegrator(Mass,Damp,Residual,dt)
    {
        m_jacobian = Jacobian;
        m_djacobian = [this](gsVector<T> const & x, gsVector<T> const & dx, gsSparseMatrix<T> & m) -> bool
        {
            return m_jacobian(x,m);
        };
        m_solver = gsSparseSolver<T>::get( "LU" ); // todo: make proper options
    }


    /**
     * @brief      Constructor
     *
     * @param[in]  Mass      The mass matrix
     * @param[in]  Jacobian  The Jacobian matrix
     * @param[in]  Residual  The residual vector
     * @param[in]  dt        The time step
     */
    gsTimeIntegrator(   const gsSparseMatrix<T> &Mass,
                        const Jacobian_t &Jacobian,
                        const TResidual_t &Residual,
                        const T dt
                        )
    : gsTimeIntegrator(Mass,Residual,dt)
    {
        m_jacobian = Jacobian;
        m_djacobian = [this](gsVector<T> const & x, gsVector<T> const & dx, gsSparseMatrix<T> & m) -> bool
        {
            return m_jacobian(x,m);
        };
        m_solver = gsSparseSolver<T>::get( "LU" ); // todo: make proper options
    }


    /**
     * @brief      Constructor
     *
     * @param[in]  Mass      The mass matrix
     * @param[in]  Damp      The damping matrix
     * @param[in]  dJacobian The Jacobian matrix taking the solution as first argument and the update as second
     * @param[in]  Residual  The residual vector
     * @param[in]  dt        The time step
     */
    gsTimeIntegrator(   const gsSparseMatrix<T> &Mass,
                        const gsSparseMatrix<T> &Damp,
                        const dJacobian_t &dJacobian,
                        const TResidual_t &Residual,
                        const T dt
                        )
    : gsTimeIntegrator(Mass,Damp,Residual,dt)
    {
      m_djacobian = dJacobian;
        m_solver = gsSparseSolver<T>::get( "LU" ); // todo: make proper options
    }


    /**
     * @brief      Constructor
     *
     * @param[in]  Mass      The mass matrix
     * @param[in]  dJacobian The Jacobian matrix taking the solution as first argument and the update as second
     * @param[in]  Residual  The residual vector
     * @param[in]  dt        The time step
     */
    gsTimeIntegrator(   const gsSparseMatrix<T> &Mass,
                        const dJacobian_t &dJacobian,
                        const TResidual_t &Residual,
                        const T dt
                        )
    : gsTimeIntegrator(Mass,Residual,dt)
    {
        m_djacobian = dJacobian;
        m_damp = gsSparseMatrix<T>(m_dofs,m_dofs);
        m_solver = gsSparseSolver<T>::get( "LU" ); // todo: make proper options
    }


public:
    /// Performs a step
    gsStatus step();

    /// Constructs the system
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

    /// Sets the time
    void setTime(T time) {m_t = time; }

    /// Sets the time
    void setTimeStep(T dt) {m_dt = dt; }
    /// Gets the time
    T currentTime() const {return m_t; }

    /// Set mass matrix
    void setMassMatrix(const gsSparseMatrix<T>& Mass)
    {
        m_mass = Mass;
        m_dofs = m_mass.cols();
    }
    /// Set damping matrix
    void setDampingMatrix(const gsSparseMatrix<T>& Damp)
    {
        m_damp = Damp;
        m_dofs = m_damp.cols();
    }
    /// Reset damping matrix
    void resetDampingMatrix()
    {
        m_damp = gsSparseMatrix<T>(m_dofs,m_dofs);
    }
    /// Set stiffness matrix
    void setStiffnessMatrix(const gsSparseMatrix<T>& Stif)
    {
        m_stif = Stif;
        m_dofs = m_mass.cols();
    }
    /// Set Jacobian matrix
    void setJacobian(Jacobian_t &Jacobian)
    {
        m_jacobian = Jacobian;
        m_djacobian = [this](gsVector<T> const & x, gsVector<T> const & dx, gsSparseMatrix<T> & m)
        {
            return m_jacobian(x,m);
        };
        m_dofs = m_mass.cols();
    }
    /// Set the Jacobian matrix taking the solution in the first argument and the update as second
    void setJacobian(dJacobian_t &dJacobian)
    {
        m_djacobian = dJacobian;
    }
    /// Resets the solution
    void resetSolution();

    /// Set the displacement
    void setDisplacement(gsMatrix<T>& displ);
    /// Set the velocity
    void setVelocity(gsMatrix<T>& velo);
    /// Set the acceleration
    void setAcceleration(gsMatrix<T>& accel);

    /// Set time integration method
    void setMethod(std::string method);

    /// Toggle verbosity
    void verbose() {m_verbose = true; }

    /// Toggle quasi newton method
    void quasiNewton() {m_quasiNewton = true; }

    /// Return the displacements
    const gsMatrix<T>& displacements() const { return m_displacements;}
    /// Return the velocities
    const gsMatrix<T>& velocities() const { return m_velocities;}
    /// Return the accelerations
    const gsMatrix<T>& accelerations() const { return m_accelerations;}

    /// Constructs the solution
    void constructSolution();

    /// Resets the iterations
    void resetIteration();

protected:

    /// Compute the time-dependent force
    gsVector<T>         _computeForce(const T & t);
    
    /// Compute the time-dependent residual
    gsVector<T>         _computeResidual(const gsVector<T> & U, const T & t);
    
    /// Compute the time-indepentent jacobian. TODO: Make time-dependent
    gsSparseMatrix<T>   _computeJacobian(const gsVector<T> & U, const gsVector<T> & dU);

    /// Factorize the matrix \a M
    void _factorizeMatrix(const gsSparseMatrix<T> & M);

    /// Solve the system with right-hand side \a F
    gsVector<T> _solveSystem(const gsVector<T> & F);


    /// Initializes some parameters for the class
    void initializeCoefficients();
    /// Initializes the solution
    void initializeSolution();
    /// Constructs the system for the explicit Euler method
    void constructSystemExplEuler();
    /// Constructs the system for the implicit Euler method
    void constructSystemImplEuler();
    /// Constructs the block-system for the implicit Euler method
    void constructSystemBlock();
    /// Initialize the solution
    void initiate();

    void _step();
    /// Perform a step with the explicit euler method
    void stepExplEuler();
    /// Perform a step with the nonlinear explicit euler method
    void stepExplEulerNL();
    /// Perform a step with the implicit euler method
    void stepImplEuler();
    /// Perform a step with the nonlinear implicit euler method
    void stepImplEulerNL();
    /// Perform a step with the nonlinear implicit euler method using block operations
    void stepImplEulerNLOp();
    /// Perform a step with the Newmark method
    void stepNewmark();
    /// Perform a step with the nonlinear Newmark method
    void stepNewmarkNL();
    /// Perform a step with the Bathe method
    void stepBathe();
    /// Perform a step with the nonlinear Bathe method
    void stepBatheNL();
    /// Perform a step with the Central difference method
    void stepCentralDiff();
    /// Perform a step with the nonlinear Central difference method
    void stepCentralDiffNL();
    /// Perform a step with the Runge-Kutta 4 method
    void stepRK4();
    /// Perform a step with the nonlinear Runge-Kutta 4 method
    void stepRK4NL();

protected:

    /// Linear solver employed
    mutable typename gsSparseSolver<T>::uPtr m_solver;

protected:

    gsSparseMatrix<T> m_mass;
    gsSparseMatrix<T> m_damp;
    gsSparseMatrix<T> m_stif;
    Jacobian_t  m_jacobian;
    dJacobian_t m_djacobian;
    TResidual_t m_residualFun;
    TForce_t    m_forceFun;
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

    std::string m_method;

    T m_updateNorm;
    T m_residue;

    // void constructSolution();
    gsMatrix<T> m_displacements;
    gsMatrix<T> m_velocities;
    gsMatrix<T> m_accelerations;

    gsStatus m_status;
};


} // namespace gismo

//////////////////////////////////////////////////
//////////////////////////////////////////////////



#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsTimeIntegrator.hpp)
#endif
