 /** @file gsStaticBase.h

    @brief Base class for static solvers

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#include <typeinfo>

#include <gsCore/gsLinearAlgebra.h>
#ifdef gsSpectra_ENABLED
#include <gsSpectra/gsSpectra.h>
#endif
#include <gsIO/gsOptionList.h>
#include <gsStructuralAnalysis/src/gsStructuralAnalysisTools/gsStructuralAnalysisTypes.h>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsStabilityAnalysis.h>

#pragma once


namespace gismo
{

/**
    @brief Base class for static solvers

    \tparam T coefficient type

    \ingroup gsStaticBase
*/
template <class T>
class gsStaticBase
{
protected:

    typedef typename gsStructuralAnalysisOps<T>::Residual_t      Residual_t;
    typedef typename gsStructuralAnalysisOps<T>::ALResidual_t    ALResidual_t;
    typedef typename gsStructuralAnalysisOps<T>::Jacobian_t      Jacobian_t;
    typedef typename gsStructuralAnalysisOps<T>::dJacobian_t     dJacobian_t;

public:

    virtual ~gsStaticBase() {};

    /// Solve. On any status other than Success, solution() is unchanged.
    virtual gsStatus solve() = 0;

    /// See \ref gsStaticBase
    virtual void initialize()
    {
        this->reset();
        this->getOptions();
    }

    /// Initialize output
    virtual void initOutput() {};
    /// Stepwise output
    virtual void stepOutput(index_t /*k*/) {};

    /// Get default options
    virtual void defaultOptions()
    {
        m_options.addReal("tol","Relative Tolerance",1e-6);
        m_options.addReal("tolF","Residual relative tolerance",-1);
        m_options.addReal("tolU","Solution step relative tolerance",-1);
        m_options.addInt("maxIt","Maximum number of iterations",25);
        m_options.addInt("verbose","Verbose output",0);
        m_options.addInt ("BifurcationMethod","Bifurcation Identification based on: 0: Determinant;  1: Eigenvalue",stabmethod::Eigenvalue);
        m_options.addString("Solver","Sparse linear solver", "SimplicialLDLT");
    }

    /// Apply the options
    virtual void getOptions()
    {
        m_tolF = m_options.getReal("tolF")!=-1 ? m_options.getReal("tolF") : m_options.getReal("tol");
        m_tolU = m_options.getReal("tolU")!=-1 ? m_options.getReal("tolU") : m_options.getReal("tol");
        m_maxIterations = m_options.getInt("maxIt");
        m_verbose = m_options.getInt("verbose");
        m_stabilityMethod     = m_options.getInt ("BifurcationMethod");
        if (m_stabilityMethod!=stabmethod::Determinant && m_stabilityMethod!=stabmethod::Eigenvalue)
            GISMO_ERROR("Unknown BifurcationMethod "<<m_stabilityMethod<<"; use 0 (determinant) or 1 (eigenvalue).");
        m_solver = gsSparseSolver<T>::get( m_options.askString("Solver","SimplicialLDLT") );
    }

    /// Set the options from \a options
    virtual void setOptions(gsOptionList & options) {m_options.update(options,gsOptionList::addIfUnknown); }

    /// Get options
    virtual const gsOptionList & options() const {return m_options;}
    virtual       gsOptionList & options()       {return m_options;}

    /// Access the solution; it is advanced only by a solve() that returns Success
    virtual gsVector<T> solution() const {return m_U;}

    /// Access the update
    virtual gsVector<T> update() const {return m_DeltaU;}

    /// Set the displacement
    virtual void setDisplacement(const gsVector<T> & displacement)
    {
        m_start = true;
        m_U = displacement;
    }

    /// Set the load
    virtual void setLoad(const T L) { m_L = L;}

    /// Set the displacement and the load
    virtual void setSolution(const gsVector<T> & displacement, const T L)
    {
        this->setDisplacement(displacement);
        this->setLoad(L);
    }
    virtual void setUpdate(const gsVector<T> & update)
    {
        m_headstart = true;
        m_DeltaU = update;
    }

    /// Returns the number of iterations
    virtual index_t iterations() const { return m_numIterations; }

    /// Returns whether the solver converged or not
    virtual bool converged() const { return m_status==gsStatus::Success; }

    /// Returns the status
    virtual gsStatus status() const { return m_status; }

    /// Returns the number of DoFs of the system
    virtual index_t numDofs() { return m_dofs; }

    /// Reset the stored solution
    virtual void reset()
    {
        m_U.setZero(m_dofs);
        m_DeltaU.setZero(m_dofs);
        m_deltaU.setZero(m_dofs);
        m_R.setZero(m_dofs);
        m_L = m_DeltaL = m_deltaL = 0.0;
        m_headstart = false;
    }

    /// Returns the stability indicator
    virtual T indicator(const gsSparseMatrix<T> & jacMat, T shift = -1e-2)
    {
       _computeStability(jacMat, shift);
       return m_indicator;
    }

    /// Returns the stability vector
    virtual gsVector<T> stabilityVec(const gsSparseMatrix<T> & jacMat, T shift = -1e-2)
    {
       _computeStability(jacMat, shift);
       return m_stabilityVec;
    }

protected:
    /// Computes the stability vector using the determinant of the Jacobian.
    /// If \a m_solver is a SimplicialLDLT, it is factorized and its pivots are used;
    /// otherwise the pivots of a local factorization are used.
    virtual bool _computeStabilityDet(const gsSparseMatrix<T> & jacMat)
    {
        if ( auto * s = dynamic_cast<typename gsSparseSolver<T>::SimplicialLDLT*>(m_solver.get()) )
        {
            m_solver->compute(jacMat);
            // If 1: matrix is not SPD
            if (m_solver->info()!=gsEigen::ComputationInfo::Success)
            {
                gsInfo<<"Solver error with code "<<m_solver->info()<<". See Eigen documentation on ComputationInfo \n"
                         <<gsEigen::ComputationInfo::Success<<": Success"<<"\n"
                         <<gsEigen::ComputationInfo::NumericalIssue<<": NumericalIssue"<<"\n"
                         <<gsEigen::ComputationInfo::NoConvergence<<": NoConvergence"<<"\n"
                         <<gsEigen::ComputationInfo::InvalidInput<<": InvalidInput"<<"\n";
                return false;
            }
            m_stabilityVec = s->vectorD();
            return m_stabilityAnalysis.classify(m_stabilityVec);
        }

        m_stabilityAnalysis.setMethod(stabmethod::Determinant);
        return m_stabilityAnalysis.compute(jacMat, T(0), m_stabilityVec);
    }

    /// Computes the stability vector using the eigenvalues of the Jacobian, optionally applying a shift
    virtual bool _computeStabilityEig(const gsSparseMatrix<T> & jacMat, T shift)
    {
        m_stabilityAnalysis.setMethod(stabmethod::Eigenvalue);
        bool ok = m_stabilityAnalysis.compute(jacMat, shift, m_stabilityVec);
        if (ok)
            m_indicator = m_stabilityAnalysis.indicator();
        return ok;
    }

    /// Computes the stability of the Jacobian, optionally applying a shift (if provided)
    virtual bool _computeStability   (const gsSparseMatrix<T> & jacMat, T shift)
    {
        bool success = false;
        if (m_stabilityMethod == stabmethod::Determinant)
            success = _computeStabilityDet(jacMat);
        else if (m_stabilityMethod == stabmethod::Eigenvalue)
            success = _computeStabilityEig(jacMat, shift);
        else
          gsInfo<<"bifurcation method unknown!";

        if (success)
            m_indicator = m_stabilityAnalysis.indicator();
        return success;
    }

protected:
    mutable gsOptionList m_options;

    T           m_L, m_DeltaL, m_deltaL;
    gsVector<T> m_U, m_DeltaU, m_deltaU;
    gsVector<T> m_R;
    T           m_residual, m_residualIni, m_residualOld;

    index_t m_numIterations;
    index_t m_maxIterations;
    bool m_start;
    index_t m_headstart;

    index_t m_verbose;

    index_t m_dofs;

    T m_tolF,m_tolU;

    T m_indicator;

    gsVector<T> m_stabilityVec;

    gsStabilityAnalysis<T> m_stabilityAnalysis;

    mutable typename gsSparseSolver<T>::uPtr  m_solver;   // Cholesky by default

    index_t m_stabilityMethod;

    struct stabmethod
    {
        enum type
        {
            Determinant = 0,
            Eigenvalue  = 1,
        };
    };

    gsStatus m_status;

};


} // namespace gismo


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsStaticBase.hpp)
#endif
