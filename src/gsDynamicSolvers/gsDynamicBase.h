 /** @file gsDynamicBase.h

    @brief Base class to perform time integration of second-order structural dynamics systems

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
3
    TODO (June 2023):
    *    Change inputs to const references!
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsIO/gsOptionList.h>
#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicMassInverse.h>
#include <gsStructuralAnalysis/src/gsStructuralAnalysisTools/gsStructuralAnalysisTypes.h>

namespace gismo
{

/**
    @brief Base class for solution strategies of nonlinear second-order dynamic 
    systems.

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

    /// @brief Callback that is evaluated before each timestep.
    typedef std::function<bool(const T /* time */, const T /* dtime */, 
        const gsOptionList & /* solver options */)> StepCallback_t;
    typedef gsDynamicMassInverse<T> MassInverse_t;
    typedef typename MassInverse_t::Ptr MassInversePtr;
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
    m_force(Force)
    {
        m_Tmass       = [this](                           const T /*time*/, gsSparseMatrix<T> & result)     -> bool {return m_mass(result);};
        m_Tdamping    = [this](gsVector<T> const & x,     const T /*time*/, gsSparseMatrix<T> & result)     -> bool {return m_damping(x,result);};
        m_Tjacobian   = [this](gsVector<T> const & /*x*/, const T /*time*/, gsSparseMatrix<T> & result)     -> bool {return m_stiffness(result);};
        m_Tforce      = [this](                           const T /*time*/, gsVector<T>       & result)     -> bool {return m_force(result);};
        m_Tresidual   = [    ](gsVector<T> const & /*x*/, const T /*time*/, gsVector<T>       & /*result*/) -> bool {GISMO_ERROR("time-dependent residual not available");};
        _init(false);
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
    m_Tforce(TForce)
    {
        m_Tmass       = [this](                           const T /*time*/, gsSparseMatrix<T> & result)     -> bool {return m_mass(result);};
        m_Tdamping    = [this](gsVector<T> const & x,     const T /*time*/, gsSparseMatrix<T> & result)     -> bool {return m_damping(x,result);};
        m_Tjacobian   = [this](gsVector<T> const & /*x*/, const T /*time*/, gsSparseMatrix<T> & result)     -> bool {return m_stiffness(result);};
        m_Tresidual   = [    ](gsVector<T> const & /*x*/, const T /*time*/, gsVector<T>       & /*result*/) -> bool {GISMO_ERROR("time-dependent residual not available");};
        _init(false);
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
    m_residual(Residual)
    {
        m_Tmass       = [this](                       const T /*time*/, gsSparseMatrix<T> & result)     -> bool {return m_mass(result);};
        m_Tdamping    = [this](gsVector<T> const & x, const T /*time*/, gsSparseMatrix<T> & result)     -> bool {return m_damping(x,result);};
        m_Tjacobian   = [this](gsVector<T> const & x, const T /*time*/, gsSparseMatrix<T> & result)     -> bool {return m_jacobian(x,result);};
        m_Tforce      = [    ](                       const T /*time*/, gsVector<T>       & /*result*/) -> bool {GISMO_ERROR("time-dependent force not available");};
        m_Tresidual   = [this](gsVector<T> const & x, const T /*time*/, gsVector<T>       & result)     -> bool {return m_residual(x,result);};
        _init(false);
    }

    /// Constructor
    gsDynamicBase(
                    const Mass_t        & Mass,
                    const Damping_t     & Damping,
                    const Jacobian_t   & Jacobian,
                    const TResidual_t   & TResidual
                )
    :
    m_mass(Mass),
    m_damping(Damping),
    m_jacobian(Jacobian),
    m_Tresidual(TResidual)
    {
        m_Tmass       = [this](                       const T /*time*/, gsSparseMatrix<T> & result) -> bool {return m_mass(result);};
        m_Tjacobian   = [this](gsVector<T> const & x, const T /*time*/, gsSparseMatrix<T> & result) -> bool {return m_jacobian(x,result);};
        m_Tdamping    = [this](gsVector<T> const & x, const T /*time*/, gsSparseMatrix<T> & result) -> bool {return m_damping(x,result);};
        _init(false);
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
    m_Tresidual(TResidual)
    {
        m_Tmass       = [this](                       const T /*time*/, gsSparseMatrix<T> & result) -> bool {return m_mass(result);};
        m_Tdamping    = [this](gsVector<T> const & x, const T /*time*/, gsSparseMatrix<T> & result) -> bool {return m_damping(x,result);};
        _init(false);
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
    m_Tresidual(TResidual)
    {
        _init(true);
    }

    gsDynamicBase() {_init(true);}

protected:
    void _init(const bool mass_timedependent)
    {
        m_time = 0;
        // initialize variables
        m_numIterations = -1;
        defaultOptions();

        m_status = gsStatus::NotStarted;

        m_massIsTimeDependent = mass_timedependent;
        m_massInverse = MassInversePtr(new gsDynamicExplicitMassInverse<T>);
        m_massInverseComputed = false;
    }

// General functions
public:

    virtual gsStatus status() { return m_status; }

    // Returns the number of DoFs
    // virtual index_t numDofs() {return m_forcing.size();}

    // Returns the current time step
    virtual T getTimeStep() const {return m_options.getReal("DT"); }

    /// Perform one time step
    virtual gsStatus step(T dt)
    {
        if (m_preStepCallback && !m_preStepCallback(m_time,dt,m_options))
        {
            m_status = gsStatus::OtherError;
            return m_status;
        }

        gsStatus status = this->_step(m_time,dt,m_U,m_V,m_A);
        m_time += dt;
        m_status = status;
        return status;
    }

    virtual gsStatus step()
    {
        return this->step(m_options.getReal("DT"));
    }

    virtual gsStatus step(const T t, const T dt, gsVector<T> & U, gsVector<T> & V, gsVector<T> & A) const
    {
        return _step(t,dt,U,V,A);
    }

    /// Set time step to \a dt
    virtual void setTimeStep(T dt)
    {
        m_options.setReal("DT",dt);
    }

    // Output
    /// True if the method converged
    virtual bool converged() const {return m_status==gsStatus::Success;}

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

    /// @brief Set a function that is executed when non-const step() is called before
    ///        the actual step is performed. 
    virtual void setPreStepCallback(const StepCallback_t & callback)
    {
        m_preStepCallback = callback;
    }

    /// @brief Clear the pre-step callback.
    virtual void clearPreStepCallback()
    {
        m_preStepCallback = StepCallback_t();
    }

    /// @brief Set a custom inverse mass function object.
    virtual void setMassInverse(const MassInversePtr & massInverse)
    {
        GISMO_ENSURE(massInverse, "Cannot set an empty inverse mass function object.");
        m_massInverse = massInverse;
        m_massInverseComputed = false;
    }

    /// @brief Restore the default explicit inverse mass function object.
    virtual void clearMassInverse()
    {
        m_massInverse = MassInversePtr(new gsDynamicExplicitMassInverse<T>);
        m_massInverseComputed = false;
    }

    /// Access the options
    virtual gsOptionList & options() {return m_options;};

    /// Set the options to \a options
    virtual void setOptions(gsOptionList options) {m_options.update(options,gsOptionList::addIfUnknown); };

    /// Return the options into \a options
    virtual void options_into(gsOptionList options) {options = m_options;};

    /// Return the number of degrees of freedom
    virtual index_t numDofs() {return m_numDofs; }
// ------------------------------------------------------------------------------------------------------------
// ---------------------------------------Computations---------------------------------------------------------
// ------------------------------------------------------------------------------------------------------------
protected:
    /// Set default options
    virtual void defaultOptions();

    /// Compute the residual
    virtual void _computeForce(const T time, gsVector<T> & F) const
    {
        if (!m_Tforce(time,F))
            throw 2;
    }

    /// Compute the residual
    virtual void _computeResidual(const gsVector<T> & U, const T time, gsVector<T> & R) const
    {
        if (!m_Tresidual(U,time,R))
            throw 2;
    }

    /// Compute the mass matrix
    virtual void _computeMass(const T time, gsSparseMatrix<T> & M) const
    {
        if (!m_Tmass(time,M))
            throw 2;
    }

    /// Return true if the inverse mass representation needs a fresh mass matrix.
    virtual bool _massInverseNeedsUpdate() const
    {
        return m_massIsTimeDependent || !m_massInverseComputed;
    }

    /// Compute or update the inverse mass representation.
    virtual void _computeMassInverse(const gsSparseMatrix<T> & M) const
    {
        if (this->_massInverseNeedsUpdate())
        {
            GISMO_ENSURE(M.rows()!=0 && M.cols()!=0, "Cannot compute inverse mass from an empty mass matrix.");
            m_massInverse->compute(M);
            m_massInverseComputed = true;
        }
    }

    /// Apply the inverse mass matrix without necessarily forming it explicitly.
    virtual void _applyMassInverse(const T time, const gsVector<T> & rhs, gsVector<T> & result) const
    {
        if (this->_massInverseNeedsUpdate())
        {
            gsSparseMatrix<T> M;
            this->_computeMass(time,M);
            this->_computeMassInverse(M);
        }
        m_massInverse->apply(rhs,result);
    }

    /// Compute the damping matrix
    virtual void _computeDamping(const gsVector<T> & U, const T time, gsSparseMatrix<T> & C) const
    {
        if (!m_Tdamping(U,time,C))
            throw 2;
    }

    /// Compute the Jacobian matrix
    virtual void _computeJacobian(const gsVector<T> & U, const T time, gsSparseMatrix<T> & K) const
    {
        if (!m_Tjacobian(U,time,K))
            throw 2;
    }

// Purely virtual functions
protected:
    /// Initialize the ALM
    virtual gsStatus _step(const T t, const T dt, gsVector<T> & U, gsVector<T> & V, gsVector<T> & A) const = 0;

protected:

    // Number of degrees of freedom
    index_t m_numDofs;

    Mass_t      m_mass;
    TMass_t     m_Tmass;
    bool m_massIsTimeDependent;
    mutable MassInversePtr m_massInverse;
    mutable bool m_massInverseComputed;

    Damping_t   m_damping;
    TDamping_t  m_Tdamping;

    Stiffness_t m_stiffness;

    Jacobian_t  m_jacobian;
    TJacobian_t m_Tjacobian;

    Force_t     m_force;
    TForce_t    m_Tforce;

    Residual_t  m_residual;
    TResidual_t m_Tresidual;

    StepCallback_t m_preStepCallback;

    mutable typename gsSparseSolver<T>::uPtr m_solver; // Cholesky by default

protected:

    mutable gsOptionList m_options;

protected:


    /// Number of iterations performed
    mutable index_t m_numIterations;

    // /// Maximum number of iterations allowed
    // index_t m_maxIterations;

    // /// Number of desired iterations
    // index_t m_desiredIterations;

    // /// Time step
    // T m_dt;

    /// Time
    T m_time;

    // /// Tolerance value to decide convergence
    // T m_tolerance;

    // /// Tolerance value to decide convergence - Force criterion
    // T m_toleranceF;

    // /// Tolerance value to decide convergence - Displacement criterion
    // T m_toleranceU;

    // bool m_verbose;
    bool m_initialized;

    // bool m_quasiNewton;
    // index_t m_quasiNewtonInterval;

    gsStatus m_status;

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
