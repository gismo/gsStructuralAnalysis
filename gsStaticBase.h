 /** @file gsStaticBase.h

    @brief Base class for static solvers

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#include <typeinfo>
#include <gsSpectra/gsSpectra.h>
#pragma once


namespace gismo
{

/**
    @brief Base class for static solvers

    \tparam T coefficient type

    \ingroup ThinShell
*/
template <class T>
class gsStaticBase
{
public:

    virtual ~gsStaticBase() {};

    /// Solve
    virtual void solve() = 0;

    /// Initialize output
    virtual void initOutput() = 0;
    /// Stepwise output
    virtual void stepOutput(index_t k) = 0;

    /// Get default options
    virtual void defaultOptions()
    {
        m_options.addReal("tol","Relative Tolerance",1e-6);
        m_options.addReal("tolF","(Force) Residual tolerance",-1);
        m_options.addReal("tolU","(Force) Residual tolerance",-1);
        m_options.addInt("maxIt","Maximum number of iterations",25);
        m_options.addInt("verbose","Verbose output",0);
        m_options.addInt ("BifurcationMethod","Bifurcation Identification based on: 0: Determinant;  1: Eigenvalue",stabmethod::Eigenvalue);
    }

    /// Apply the options
    virtual void getOptions()
    {
        m_tolF = m_options.getReal("tolF")!=-1 ? m_options.getReal("tolF") : m_options.getReal("tol");
        m_tolU = m_options.getReal("tolU")!=-1 ? m_options.getReal("tolU") : m_options.getReal("tol");
        m_maxIterations = m_options.getInt("maxIt");
        m_verbose = m_options.getInt("verbose");
        m_stabilityMethod     = m_options.getInt ("BifurcationMethod");
    }

    /// Set the options from \a options
    void setOptions(gsOptionList & options) {m_options.update(options,gsOptionList::addIfUnknown); }

    /// Get options
    gsOptionList options() const {return m_options;}

    /// Access the solution
    gsVector<T> solution() const {return m_U;}

    /// Access the update
    gsVector<T> update() const {return m_DeltaU;}

    /// Set the displacement
    void setDisplacement(const gsVector<T> displacement)
    {
        m_start = true;
        m_U = displacement;
    }

    /// Set the load
    void setLoad(const T L) { m_L = L;}

    /// Set the displacement and the load
    void setSolution(const gsVector<T> displacement, const T L)
    {
        this->setDisplacement(displacement);
        this->setLoad(L);
    }
    virtual void setUpdate(const gsVector<T> update)
    {
        m_headstart = true;
        m_DeltaU = update;
    }

    /// Returns the number of iterations
    index_t iterations() const { return m_numIterations; }

    /// Returns whether the solver converged or not
    bool converged() const { return m_converged; }

    /// Returns the number of DoFs of the system
    index_t numDofs() { return m_dofs; }

    /// Reset the stored solution
    void reset()
    {
        m_U.setZero();
        m_DeltaU.setZero();
        m_deltaU.setZero();
        m_L = m_DeltaL = m_deltaL = 0.0;
        m_headstart = false;
    }

    /// Returns the stability indicator
    T indicator(const gsSparseMatrix<T> & jacMat, T shift = -1e-2)
    {
       _computeStability(jacMat, shift);
       return m_indicator;
    }

    /// Returns the stability vector
    gsVector<T> stabilityVec(const gsSparseMatrix<T> & jacMat, T shift = -1e-2)
    {
       _computeStability(jacMat, shift);
       return m_stabilityVec;
    }

protected:
    /// Computes the stability vector using the determinant of the Jacobian
    virtual void _computeStabilityDet(const gsSparseMatrix<T> & jacMat)
    {
        m_LDLTsolver.compute(jacMat);
        // If 1: matrix is not SPD
        GISMO_ASSERT(m_LDLTsolver.info()==Eigen::ComputationInfo::Success,"Solver error with code "<<m_LDLTsolver.info()<<". See Eigen documentation on ComputationInfo \n"
                                                                  <<Eigen::ComputationInfo::Success<<": Success"<<"\n"
                                                                  <<Eigen::ComputationInfo::NumericalIssue<<": NumericalIssue"<<"\n"
                                                                  <<Eigen::ComputationInfo::NoConvergence<<": NoConvergence"<<"\n"
                                                                  <<Eigen::ComputationInfo::InvalidInput<<": InvalidInput"<<"\n");

        m_stabilityVec = m_LDLTsolver.vectorD();
    }

    /// Computes the stability vector using the eigenvalues of the Jacobian, optionally applying a shift
    virtual void _computeStabilityEig(const gsSparseMatrix<T> & jacMat, T shift)
    {
#ifdef GISMO_WITH_SPECTRA
        index_t number = std::min(static_cast<index_t>(std::floor(jacMat.cols()/5.)),10);
        /*
        // Without shift!
        // This one can sometimes not converge, because spectra is better at finding large values.
          gsSpectraSymSolver<gsSparseMatrix<T>> es(jacMat,number,5*number);
          es.init();
          es.compute(Spectra::SortRule::SmallestAlge,1000,1e-6,Spectra::SortRule::SmallestAlge);
          GISMO_ASSERT(es.info()==Spectra::CompInfo::Successful,"Spectra did not converge!"); // Reason for not converging can be due to the value of ncv (last input in the class member), which is too low.
        */

        // With shift!
        // This one converges easier. However, a shift must be provided!
        gsSpectraSymShiftSolver<gsSparseMatrix<T>> es(jacMat,number,5*number,shift);
        es.init();
        es.compute(Spectra::SortRule::LargestAlge,1000,1e-6,Spectra::SortRule::SmallestAlge);
        GISMO_ENSURE(es.info()==Spectra::CompInfo::Successful,"Spectra did not converge!"); // Reason for not converging can be due to the value of ncv (last input in the class member), which is too low.

        // if (es.info()==Spectra::CompInfo::NotComputed)
        // if (es.info()==Spectra::CompInfo::NotConverging)
        // if (es.info()==Spectra::CompInfo::NumericalIssue)
        // Eigen::SelfAdjointEigenSolver< gsMatrix<T> > es(jacMat);
        m_stabilityVec = es.eigenvalues();
#else
        Eigen::SelfAdjointEigenSolver<gsMatrix<T>> es2(jacMat);
        m_stabilityVec = es2.eigenvalues();
#endif

        m_indicator = m_stabilityVec.colwise().minCoeff()[0]; // This is required since D does not necessarily have one column.
    }

    /// Computes the stability of the Jacobian, optionally applying a shift (if provided)
    virtual void _computeStability   (const gsSparseMatrix<T> & jacMat, T shift)
    {
        if (m_stabilityMethod == stabmethod::Determinant)
        {
            _computeStabilityDet(jacMat);
        }
        else if (m_stabilityMethod == stabmethod::Eigenvalue)
        {
            _computeStabilityEig(jacMat, shift);
        }
        else
          gsInfo<<"bifurcation method unknown!";

        m_indicator = m_stabilityVec.colwise().minCoeff()[0]; // This is required since D does not necessarily have one column
    }

protected:
    mutable gsOptionList m_options;

    T           m_L, m_DeltaL, m_deltaL;
    gsVector<T> m_U, m_DeltaU, m_deltaU;
    gsVector<T> m_R;
    T           m_residual, m_residualIni, m_residualOld;

    index_t m_numIterations;
    index_t m_maxIterations;
    bool m_converged;
    bool m_start;
    index_t m_headstart;

    index_t m_verbose;

    index_t m_dofs;

    T m_tolF,m_tolU;

    T m_indicator;

    gsVector<T> m_stabilityVec;

    mutable gsSparseSolver<>::SimplicialLDLT  m_LDLTsolver;   // Cholesky

    index_t m_stabilityMethod;

    struct stabmethod
    {
        enum type
        {
            Determinant = 0,
            Eigenvalue  = 1,
        };
    };

};


} // namespace gismo


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsStaticBase.hpp)
#endif