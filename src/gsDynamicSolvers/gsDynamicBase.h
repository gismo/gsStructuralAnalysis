 /** @file gsDynamicBase.h

    @brief Base class to perform time integration of second-order structural dynamics systems

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)

    TODO (June 2023):
    *    Change inputs to const references!
*/

#pragma once
#include <gsCore/gsLinearAlgebra.h>

#include <gsIO/gsOptionList.h>
#include <gsStructuralAnalysis/src/gsStructuralAnalysisTools/gsStructuralAnalysisTypes.h>

namespace gismo
{

/**
    @brief Performs the arc length method to solve a nonlinear system of equations.

    \tparam T coefficient type

    \ingroup gsStructuralAnalysis
*/
template <class T>
class gsDynamicBase
{
protected:

    typedef typename gsStructuralAnalysisOps<T>::Force_t     Force_t;
    typedef typename gsStructuralAnalysisOps<T>::TForce_t    TForce_t;
    typedef typename gsStructuralAnalysisOps<T>::Residual_t  Residual_t;
    typedef typename gsStructuralAnalysisOps<T>::TResidual_t TResidual_t;
    typedef typename gsStructuralAnalysisOps<T>::Mass_t      Mass_t;
    typedef typename gsStructuralAnalysisOps<T>::TMass_t     TMass_t;
    typedef typename gsStructuralAnalysisOps<T>::Damping_t   Damping_t;
    typedef typename gsStructuralAnalysisOps<T>::TDamping_t  TDamping_t;
    typedef typename gsStructuralAnalysisOps<T>::Stiffness_t Stiffness_t;
    typedef typename gsStructuralAnalysisOps<T>::Jacobian_t  Jacobian_t;
    typedef typename gsStructuralAnalysisOps<T>::TJacobian_t TJacobian_t;

public:

    virtual ~gsDynamicBase() {};

    /// Constructor
    gsDynamicBase(
                    const Mass_t        & Mass,
                    const Damping_t     & Damping,
                    const Stiffness_t   & Stiffness,
                    const Force_t       & Force
                )
    :
    m_mass(Mass),
    m_damping(Damping),
    m_stiffness(Stiffness),
    m_force(Force),
    m_nonlinear(false)
    {
        m_Tmass       = [this](                       const T time, gsSparseMatrix<T> & result) -> bool {return m_mass(result);};
        m_Tdamping    = [this](gsVector<T> const & x, const T time, gsSparseMatrix<T> & result) -> bool {return m_damping(x,result);};
        m_Tjacobian   = [this](gsVector<T> const & x, const T time, gsSparseMatrix<T> & result) -> bool {return m_stiffness(result);};
        m_Tforce      = [this](                       const T time, gsVector<T>       & result) -> bool {return m_force(result);};
        m_Tresidual   = [this](gsVector<T> const & x, const T time, gsVector<T>       & result) -> bool {GISMO_ERROR("time-dependent residual not available");};
        _init();
    }

    /// Constructor
    gsDynamicBase(
                    const Mass_t        & Mass,
                    const Damping_t     & Damping,
                    const Stiffness_t   & Stiffness,
                    const TForce_t      & TForce
                )
    :
    m_mass(Mass),
    m_damping(Damping),
    m_stiffness(Stiffness),
    m_Tforce(TForce),
    m_nonlinear(false)
    {
        m_Tmass       = [this](                       const T time, gsSparseMatrix<T> & result) -> bool {return m_mass(result);};
        m_Tdamping    = [this](gsVector<T> const & x, const T time, gsSparseMatrix<T> & result) -> bool {return m_damping(x,result);};
        m_Tjacobian   = [this](gsVector<T> const & x, const T time, gsSparseMatrix<T> & result) -> bool {return m_stiffness(result);};
        m_Tresidual   = [this](gsVector<T> const & x, const T time, gsVector<T>       & result) -> bool {GISMO_ERROR("time-dependent residual not available");};
        _init();
    }

    /// Constructor
    gsDynamicBase(
                    const Mass_t        & Mass,
                    const Damping_t     & Damping,
                    const Jacobian_t    & Jacobian,
                    const Residual_t    & Residual
                )
    :
    m_mass(Mass),
    m_damping(Damping),
    m_jacobian(Jacobian),
    m_residual(Residual),
    m_nonlinear(true)
    {
        m_Tmass       = [this](                       const T time, gsSparseMatrix<T> & result) -> bool {return m_mass(result);};
        m_Tdamping    = [this](gsVector<T> const & x, const T time, gsSparseMatrix<T> & result) -> bool {return m_damping(x,result);};
        m_Tjacobian   = [this](gsVector<T> const & x, const T time, gsSparseMatrix<T> & result) -> bool {return m_jacobian(x,result);};
        m_Tforce      = [this](                       const T time, gsVector<T>       & result) -> bool {GISMO_ERROR("time-dependent force not available");};
        m_Tresidual   = [this](gsVector<T> const & x, const T time, gsVector<T>       & result) -> bool {return m_residual(x,result);};
        _init();
    }

    /// Constructor
    gsDynamicBase(
                    const Mass_t        & Mass,
                    const Damping_t     & Damping,
                    const TJacobian_t   & TJacobian,
                    const TResidual_t   & TResidual
                )
    :
    m_mass(Mass),
    m_damping(Damping),
    m_Tjacobian(TJacobian),
    m_Tresidual(TResidual),
    m_nonlinear(true)
    {
        m_Tmass       = [this](                       const T time, gsSparseMatrix<T> & result) -> bool {return m_mass(result);};
        m_Tdamping    = [this](gsVector<T> const & x, const T time, gsSparseMatrix<T> & result) -> bool {return m_damping(x,result);};
        _init();
    }

    /// Constructor
    gsDynamicBase(
                    const TMass_t       & TMass,
                    const TDamping_t    & TDamping,
                    const TJacobian_t   & TJacobian,
                    const TResidual_t   & TResidual
                )
    :
    m_Tmass(TMass),
    m_Tdamping(TDamping),
    m_Tjacobian(TJacobian),
    m_Tresidual(TResidual),
    m_nonlinear(true)
    {
        _init();
    }

    gsDynamicBase() {_init();}

protected:
    void _init()
    {
        // initialize variables
        m_numIterations = 0;
        this->defaultOptions();
        m_converged = false;
        m_initialized = false;

        m_status = gsStatus::NotStarted;
    }

// General functions
public:

    virtual gsStatus status() { return m_status; }

    // Returns the number of DoFs
    // virtual index_t numDofs() {return m_forcing.size();}

    // Returns the current time step
    virtual T getTimeStep() const {return m_dt; }

    /// Perform one arc-length step
    virtual gsStatus step(T dt)
    {
        return this->_step(dt);
    }

    virtual gsStatus step()
    {
        return this->step(m_dt);
    }

    /// Initialize the arc-length method, computes the stability of the initial configuration if \a stability is true
    virtual void initialize()
    {
        m_initialized = true;
        this->getOptions();
    }

    /// Set time step to \a dt
    virtual void setTimeStep(T dt)
    {
      // m_options.setReal("Length",length);
      // m_arcLength = m_arcLength_prev = m_arcLength_ori = m_options.getReal("Length");
    }

    // Output
    /// True if the Arc Length method converged
    virtual bool converged() const {return m_converged;}

    /// Returns the number of Newton iterations performed
    virtual index_t numIterations() const { return m_numIterations;}

    virtual const gsVector<T> & displacements() const {return this->solutionU();}
    virtual const gsVector<T> & velocities()    const {return this->solutionV();}
    virtual const gsVector<T> & accelerations() const {return this->solutionA();}

    virtual const gsVector<T> & solutionU() const {return m_U;}
    virtual const gsVector<T> & solutionV() const {return m_V;}
    virtual const gsVector<T> & solutionA() const {return m_A;}

    virtual const T & time() {return m_time;}

    virtual void setU(const gsVector<T> & U) {m_U = U;}
    virtual void setV(const gsVector<T> & V) {m_V = V;}
    virtual void setA(const gsVector<T> & A) {m_A = A;}

    virtual void setDisplacements(const gsVector<T> & U) {this->setU(U);}
    virtual void setVelocities(const gsVector<T> & V)    {this->setV(V);}
    virtual void setAccelerations(const gsVector<T> & A) {this->setA(A);}

    /// Access the options
    virtual gsOptionList & options() {return m_options;};

    /// Set the options to \a options
    virtual void setOptions(gsOptionList options) {m_options.update(options,gsOptionList::addIfUnknown); this->getOptions(); };

    /// Return the options into \a options
    virtual const void options_into(gsOptionList options) {options = m_options;};

    /// Apply the options
    virtual void applyOptions() {this->getOptions(); }

// ------------------------------------------------------------------------------------------------------------
// ---------------------------------------Computations---------------------------------------------------------
// ------------------------------------------------------------------------------------------------------------
protected:
    /// Set default options
    virtual void defaultOptions();

    /// Apply options
    virtual void getOptions();

    /// Compute the residual
    virtual void _computeForce(const T time, gsVector<T> & F)
    {
        if (!m_Tforce(time,F))
            throw 2;
    }

    /// Compute the residual
    virtual void _computeResidual(const gsVector<T> & U, const T time, gsVector<T> & R)
    {
        if (!m_Tresidual(U,time,R))
            throw 2;
    }

    /// Compute the mass matrix
    virtual void _computeMass(const T time, gsSparseMatrix<T> & M)
    {
        if (!m_Tmass(time,M))
            throw 2;
    }

    /// Compute the damping matrix
    virtual void _computeDamping(const gsVector<T> & U, const T time, gsSparseMatrix<T> & C)
    {
        if (!m_Tdamping(U,time,C))
            throw 2;
    }

    /// Compute the Jacobian matrix
    virtual void _computeJacobian(const gsVector<T> & U, const T time, gsSparseMatrix<T> & K)
    {
        if (!m_Tjacobian(U,time,K))
            throw 2;
    }

// Purely virtual functions
protected:
    /// Initialize the ALM
    virtual gsStatus _step(const T dt) = 0;

protected:

    // Number of degrees of freedom
    index_t m_numDof;

    Mass_t      m_mass;
    TMass_t     m_Tmass;

    Damping_t   m_damping;
    TDamping_t  m_Tdamping;

    Stiffness_t m_stiffness;

    Jacobian_t  m_jacobian;
    TJacobian_t m_Tjacobian;

    Force_t     m_force;
    TForce_t    m_Tforce;

    Residual_t  m_residual;
    TResidual_t m_Tresidual;

    bool m_nonlinear;

    mutable typename gsSparseSolver<T>::uPtr m_solver; // Cholesky by default

protected:

    mutable gsOptionList m_options;

protected:


    /// Number of Arc Length iterations performed
    index_t m_numIterations;

    /// Maximum number of Arc Length iterations allowed
    index_t m_maxIterations;

    /// Number of desired iterations
    index_t m_desiredIterations;

    /// Time step
    T m_dt;

    /// Time
    T m_time;

    /// Tolerance value to decide convergence
    T m_tolerance;

    /// Tolerance value to decide convergence - Force criterion
    T m_toleranceF;

    /// Tolerance value to decide convergence - Displacement criterion
    T m_toleranceU;

    bool m_verbose;
    bool m_initialized;

    bool m_quasiNewton;
    index_t m_quasiNewtonInterval;

    gsStatus m_status;

protected:

    /// Convergence result
    bool m_converged;

    /// Force residuum
    T m_residualNorm;

    /// Update norm
    T m_updateNorm;

protected:

    // /// Current update
    // gsVector<T> m_DeltaU;
    // gsVector<T> m_DeltaV;
    // gsVector<T> m_DeltaA;

    // /// Iterative update
    // gsVector<T> m_deltaU;
    // gsVector<T> m_deltaV;
    // gsVector<T> m_deltaA;

    // Previous update
    gsVector<T> m_DeltaUold;
    gsVector<T> m_DeltaVold;
    gsVector<T> m_DeltaAold;

    /// Displacement vector (present, at previously converged point)
    gsVector<T> m_U, m_Uold;
    gsVector<T> m_V, m_Vold;
    gsVector<T> m_A, m_Aold;

    /// Time
    T m_t, m_tprev;
};


} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsDynamicBase.hpp)
#endif
