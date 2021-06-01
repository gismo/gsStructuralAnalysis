 /** @file gsStaticSolver.hpp

    @brief Performs linear modal analysis given a matrix or functions of a matrix

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#include <typeinfo>
#pragma once

namespace gismo
{

template <class T>
gsVector<T> gsStaticSolver<T>::solveLinear()
{
    this->getOptions();
    if (m_verbose==2)
    {
        gsInfo<<"Matrix: \n"<<m_linear.toDense()<<"\n";
        gsInfo<<"Vector: \n"<<m_force<<"\n";
    }
    factorizeMatrix( m_linear );
    m_solVec = solveSystem(m_force);
    m_converged = true;

    return m_solVec;
};

template <class T>
gsVector<T> gsStaticSolver<T>::solveNonlinear()
{
    m_converged = false;
    this->getOptions();
    if (m_solVec.rows()==0)
        m_solVec = this->solveLinear();

    gsVector<T> resVec;
    T residual,residual0,residualOld;
    if (!m_headstart) // no headstart
    {
        m_DeltaU.setZero(m_solVec.rows());// ->>>>>>> THIS GOES WRONG WITH THE NORMS?!
        resVec = m_residual(m_solVec);
        if (residual==0) residual=1;
        residual0 = residualOld = resVec.norm();
    }
    else
    {
        resVec = m_residual(m_solVec + m_DeltaU);
        residual = resVec.norm();
        if (residual==0) residual=1;
        residualOld = residual;
        residual0 = m_residual(m_solVec).norm();
        if (residual0==0) residual0=1;
    }

    // gsDebugVar(resVec.transpose());

    m_headstart = false;

    m_deltaU.setZero(m_solVec.rows());

    gsSparseMatrix<T> jacMat;

    if (m_verbose>0)
    {
        gsInfo<<"\t";
        gsInfo<<std::setw(4)<<std::left<<"It.";
        gsInfo<<std::setw(17)<<std::left<<"|R|";
        gsInfo<<std::setw(17)<<std::left<<"|R|/|R0|";
        gsInfo<<std::setw(17)<<std::left<<"|dU|";
        gsInfo<<std::setw(17)<<std::left<<"|dU|/|DU|";
        gsInfo<<std::setw(17)<<std::left<<"|dU|/|U|";
        gsInfo<<std::setw(17)<<std::left<<"log(Ri/R0):";
        gsInfo<<std::setw(17)<<std::left<<"log(Ri+1/R0)";
        gsInfo<<"\n";
    }

    for (m_iterations = 0; m_iterations != m_maxIterations; ++m_iterations)
    {
        jacMat = m_nonlinear(m_solVec+m_DeltaU);
        if (m_verbose==2)
        {
            gsInfo<<"Matrix: \n"<<jacMat.toDense()<<"\n";
            gsInfo<<"Vector: \n"<<resVec<<"\n";
        }
        factorizeMatrix(jacMat);
        m_deltaU = solveSystem(resVec); // this is the UPDATE
        m_DeltaU += m_relax * m_deltaU;

        resVec = m_residual(m_solVec+m_DeltaU);
        residual = resVec.norm();

        if (m_verbose>0)
        {
            gsInfo<<"\t";
            gsInfo<<std::setw(4)<<std::left<<m_iterations;
            gsInfo<<std::setw(17)<<std::left<<residual;
            gsInfo<<std::setw(17)<<std::left<<residual/residual0;
            gsInfo<<std::setw(17)<<std::left<<m_relax * m_deltaU.norm();
            gsInfo<<std::setw(17)<<std::left<<m_relax * m_deltaU.norm()/m_DeltaU.norm();
            gsInfo<<std::setw(17)<<std::left<<m_relax * m_deltaU.norm()/m_solVec.norm();
            gsInfo<<std::setw(17)<<std::left<<math::log10(residualOld/residual0);
            gsInfo<<std::setw(17)<<std::left<<math::log10(residual/residual0);
            gsInfo<<"\n";
        }

        residualOld = residual;

        if (m_relax * m_deltaU.norm()/m_solVec.norm()  < m_toleranceU && residual/residual0 < m_toleranceF)
        {
            m_converged = true;
            m_solVec+=m_DeltaU;
            break;
        }
        else if (m_iterations+1 == m_maxIterations)
        {
            m_converged = false;
            gsWarn<<"Maximum iterations reached!\n";
        }


            // ADD DIRICHLET HOMOGENIZE
    }
    return m_solVec;
};

template <class T>
void gsStaticSolver<T>::_computeStability(const gsVector<T> x) const
{
    gsSparseMatrix<T> jacMat = m_nonlinear(x);
    // gsInfo<<"x = \n"<<x.transpose()<<"\n";
    if (m_bifurcationMethod == bifmethod::Determinant)
    {
      factorizeMatrix(jacMat);
      m_stabilityVec = m_LDLTsolver.vectorD();
    }
    else if (m_bifurcationMethod == bifmethod::Eigenvalue)
    {
      #ifdef GISMO_WITH_SPECTRA
      index_t number = std::min(static_cast<index_t>(std::floor(jacMat.cols()/3.)),10);
      gsSpectraSymSolver<gsSparseMatrix<T>> es(jacMat,number,5*number);
      es.init();
      es.compute(Spectra::SortRule::SmallestAlge,1000,1e-6,Spectra::SortRule::SmallestAlge);
      GISMO_ASSERT(es.info()==Spectra::CompInfo::Successful,"Spectra did not converge!"); // Reason for not converging can be due to the value of ncv (last input in the class member), which is too low.
      // if (es.info()==Spectra::CompInfo::NotComputed)
      // if (es.info()==Spectra::CompInfo::NotConverging)
      // if (es.info()==Spectra::CompInfo::NumericalIssue)
      // Eigen::SelfAdjointEigenSolver< gsMatrix<T> > es(m_jacMat);
      m_stabilityVec = es.eigenvalues();
      #else
      Eigen::SelfAdjointEigenSolver<gsMatrix<T>> es2(jacMat);
      m_stabilityVec = es2.eigenvalues();
      #endif
    }
    else
      gsInfo<<"bifurcation method unknown!";

    m_indicator = m_stabilityVec.colwise().minCoeff()[0]; // This is required since D does not necessarily have one column
}

template <class T>
void gsStaticSolver<T>::defaultOptions()
{
    m_options.addInt("Verbose","Verbose output",verbose::iterations);
    m_options.addInt("MaxIterations","Maximum number of iterations",25);
    m_options.addReal("Tolerance","Relative Tolerance Force",1e-6);
    m_options.addReal("ToleranceF","Relative Tolerance Force",-1);
    m_options.addReal("ToleranceU","Relative Tolerance Displacements",-1);
    m_options.addReal("Relaxation","Relaxation parameter",1);
    m_options.addInt ("BifurcationMethod","Bifurcation Identification based on: 0: Determinant;  1: Eigenvalue",bifmethod::Eigenvalue);
    m_options.addInt ("Solver","Linear solver: 0 = LDLT; 1 = CG",solver::CG); // The CG solver is robust for membrane models, where zero-blocks in the matrix might occur.

};

template <class T>
void gsStaticSolver<T>::getOptions() const
{
    m_maxIterations = m_options.getInt("MaxIterations");
    m_toleranceF = m_options.getReal("ToleranceF")!=-1 ? m_options.getReal("ToleranceF") : m_options.getReal("Tolerance");
    m_toleranceU = m_options.getReal("ToleranceU")!=-1 ? m_options.getReal("ToleranceU") : m_options.getReal("Tolerance");
    m_verbose = m_options.getInt("Verbose");
    m_relax = m_options.getReal("Relaxation");

    m_bifurcationMethod   = m_options.getInt ("BifurcationMethod");
    m_solverType          = m_options.getInt ("Solver");
    if (m_solverType!=solver::LDLT && m_bifurcationMethod==bifmethod::Determinant)
    {
      gsWarn<<"Determinant method cannot be used with solvers other than LDLT. Bifurcation method will be set to 'Eigenvalue'.\n";
      m_bifurcationMethod = bifmethod::Eigenvalue;
    }
};

template <class T>
void gsStaticSolver<T>::factorizeMatrix(const gsSparseMatrix<T> & M) const
{
  if (m_solverType==solver::LDLT)
  {
    m_LDLTsolver.compute(M);
    // If 1: matrix is not SPD
    GISMO_ASSERT(m_LDLTsolver.info()==Eigen::ComputationInfo::Success,"Solver error with code "<<m_LDLTsolver.info()<<". See Eigen documentation on ComputationInfo \n"
                                                                <<Eigen::ComputationInfo::Success<<": Success"<<"\n"
                                                                <<Eigen::ComputationInfo::NumericalIssue<<": NumericalIssue"<<"\n"
                                                                <<Eigen::ComputationInfo::NoConvergence<<": NoConvergence"<<"\n"
                                                                <<Eigen::ComputationInfo::InvalidInput<<": InvalidInput"<<"\n");

  }
  else if (m_solverType==solver::CG)
  {
    m_CGsolver.compute(M);

    GISMO_ASSERT(m_CGsolver.info()==Eigen::ComputationInfo::Success,"Solver error with code "<<m_CGsolver.info()<<". See Eigen documentation on ComputationInfo \n"
                                                                <<Eigen::ComputationInfo::Success<<": Success"<<"\n"
                                                                <<Eigen::ComputationInfo::NumericalIssue<<": NumericalIssue"<<"\n"
                                                                <<Eigen::ComputationInfo::NoConvergence<<": NoConvergence"<<"\n"
                                                                <<Eigen::ComputationInfo::InvalidInput<<": InvalidInput"<<"\n");

  }
  else
    GISMO_ERROR("Solver type "<<m_solverType<<" unknown.");

}

template <class T>
gsVector<T> gsStaticSolver<T>::solveSystem(const gsVector<T> & F)
{
  if (m_solverType==solver::LDLT)
    return m_LDLTsolver.solve(F);
  else if (m_solverType==solver::CG)
    return m_CGsolver.solve(F);
  else
    GISMO_ERROR("Solver type "<<m_solverType<<" unknown.");
}

} // namespace gismo
