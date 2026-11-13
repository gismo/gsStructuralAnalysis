 /** @file gsStaticNewton.hpp

    @brief Performs linear modal analysis given a matrix or functions of a matrix

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
void gsStaticNewton<T>::defaultOptions()
{
    Base::defaultOptions();
    m_options.setString("Solver","CGDiagonal"); // The CG solver is robust for membrane models, where zero-blocks in the matrix might occur.
    m_options.addReal("Relaxation","Relaxation parameter",1);
    m_options.addSwitch("LineSearch","Backtracking (Armijo) line search on the residual norm",false);
    m_options.addReal("LineSearchC","Line search: sufficient decrease constant c, 0 < c < 1",1e-4);
    m_options.addReal("LineSearchRho","Line search: step reduction factor rho, 0 < rho < 1",0.5);
    m_options.addInt("LineSearchMaxRed","Line search: maximum number of step reductions",8);
};

template <class T>
void gsStaticNewton<T>::getOptions()
{
    Base::getOptions();

    // if (m_solverType!=solver::LDLT && m_stabilityMethod==stabmethod::Determinant)
    // {
    //   gsWarn<<"Determinant method cannot be used with solvers other than LDLT. Bifurcation method will be set to 'Eigenvalue'.\n";
    //   m_stabilityMethod = stabmethod::Eigenvalue;
    // }

    m_relax = m_options.getReal("Relaxation");

    m_lineSearch = m_options.getSwitch("LineSearch");
    m_lsC        = m_options.getReal("LineSearchC");
    m_lsRho      = m_options.getReal("LineSearchRho");
    m_lsMaxRed   = m_options.getInt("LineSearchMaxRed");
    if (m_lineSearch)
    {
        GISMO_ENSURE(m_lsC>0 && m_lsC<1,"LineSearchC must be in (0,1), got "<<m_lsC);
        GISMO_ENSURE(m_lsRho>0 && m_lsRho<1,"LineSearchRho must be in (0,1), got "<<m_lsRho);
        GISMO_ENSURE(m_lsMaxRed>=0,"LineSearchMaxRed must be non-negative, got "<<m_lsMaxRed);
    }
};

template <class T>
void gsStaticNewton<T>::initOutput()
{
    gsInfo<<"\t";
    gsInfo<<std::setw(4)<<std::left<<"It.";
    gsInfo<<std::setw(17)<<std::left<<"|R|";
    gsInfo<<std::setw(17)<<std::left<<"|R|/|R0|";
    gsInfo<<std::setw(17)<<std::left<<"|DU|";
    gsInfo<<std::setw(17)<<std::left<<"|dU|";
    gsInfo<<std::setw(17)<<std::left<<"|dU|/|DU|";
    gsInfo<<std::setw(17)<<std::left<<"|dU|/|U+DU|";
    gsInfo<<std::setw(17)<<std::left<<"log(Ri/R0):";
    gsInfo<<std::setw(17)<<std::left<<"log(Ri+1/R0)";
    if (m_lineSearch)
        gsInfo<<std::setw(17)<<std::left<<"alpha";
    gsInfo<<"\n";
}

template <class T>
void gsStaticNewton<T>::stepOutput(index_t k)
{
    gsInfo<<"\t";
    gsInfo<<std::setw(4)<<std::left<<k;
    gsInfo<<std::setw(17)<<std::left<<m_residual;
    gsInfo<<std::setw(17)<<std::left<<m_residual/m_residualIni;
    gsInfo<<std::setw(17)<<std::left<<m_DeltaU.norm();
    gsInfo<<std::setw(17)<<std::left<<m_relax * m_deltaU.norm();
    gsInfo<<std::setw(17)<<std::left<<m_relax * m_deltaU.norm()/m_DeltaU.norm();
    gsInfo<<std::setw(17)<<std::left<<m_relax * m_deltaU.norm()/(m_U+m_DeltaU).norm();
    gsInfo<<std::setw(17)<<std::left<<math::log10(m_residualOld/m_residualIni);
    gsInfo<<std::setw(17)<<std::left<<math::log10(m_residual/m_residualIni);
    if (m_lineSearch)
        gsInfo<<std::setw(17)<<std::left<<m_alpha;
    gsInfo<<"\n";

}

template <class T>
gsStatus gsStaticNewton<T>::solveLinear()
{
    try
    {
        _solveLinear();
        m_status = gsStatus::Success;
    }
    catch (int errorCode)
    {
        if      (errorCode==1)
            m_status = gsStatus::NotConverged;
        else if (errorCode==2)
            m_status = gsStatus::AssemblyError;
        else if (errorCode==3)
            m_status = gsStatus::SolverError;
        else
            m_status = gsStatus::OtherError;
    }
    catch (...)
    {
        m_status = gsStatus::OtherError;
    }
    return m_status;
};

template <class T>
gsVector<T> gsStaticNewton<T>::_solveLinear()
{
    this->getOptions();
    if (m_verbose==2)
    {
        gsInfo<<"Matrix: \n"<<m_linear.toDense()<<"\n";
        gsInfo<<"Vector: \n"<<m_force<<"\n";
    }
    _factorizeMatrix( m_linear );
    m_U = _solveSystem(m_force);
    return m_U;
};

template <class T>
gsStatus gsStaticNewton<T>::solveNonlinear()
{
    try
    {
        _solveNonlinear();
        m_status = gsStatus::Success;
    }
    catch (int errorCode)
    {
        if      (errorCode==1)
            m_status = gsStatus::NotConverged;
        else if (errorCode==2)
            m_status = gsStatus::AssemblyError;
        else if (errorCode==3)
            m_status = gsStatus::SolverError;
        else
            m_status = gsStatus::OtherError;
    }
    catch (...)
    {
        m_status = gsStatus::OtherError;
    }
    return m_status;
}


template <class T>
gsVector<T> gsStaticNewton<T>::_solveNonlinear()
{
    this->getOptions();
    // m_start: true -> m_U given
    // m_headstart: true -> m_DeltaU given

    if (m_DeltaU.norm()==0 && m_DeltaU.rows()==0) ///
    {
        m_deltaU = m_DeltaU = this->_solveLinear();
        m_U.setZero(); // Needed because linear solve modifies m_U.
        m_headstart = true; // due to this, the relative residual is based on the solution of the linear solve
    }
    _start();

    gsSparseMatrix<T> jacMat;

    if (m_verbose>0) { initOutput(); }

    for (m_numIterations = 0; m_numIterations != m_maxIterations; ++m_numIterations)
    {
        jacMat = this->_computeJacobian(m_U+m_DeltaU,m_deltaU);
        if (m_verbose==2)
        {
            gsInfo<<"Matrix: \n"<<jacMat.toDense()<<"\n";
            gsInfo<<"Vector: \n"<<m_R<<"\n";
        }

        _factorizeMatrix(jacMat);
        m_deltaU = _solveSystem(m_R); // this is the UPDATE
        if (!m_lineSearch)
        {
            m_DeltaU += m_relax * m_deltaU;

            m_R = this->_computeResidual(m_U+m_DeltaU);
        }
        else
        {
            // Backtracking on the residual norm: accept the first alpha = rho^j,
            // j = 0..maxRed, with finite trial residual and
            //   ||R(u + alpha du)|| <= (1 - c alpha) ||R(u)||   (Armijo)
            // At most maxRed extra residual assemblies per iteration.
            const T R0n = m_R.norm();
            // At convergence level the decrease test is dominated by round-off.
            const bool converged = (R0n==0) || (R0n/m_residualIni < m_tolF);
            gsVector<T> Rt;
            bool accepted = false;
            T alpha = 1;
            for (index_t j = 0; j <= m_lsMaxRed; ++j, alpha *= m_lsRho)
            {
                bool failed = false;
                try { Rt = this->_computeResidual(m_U + m_DeltaU + alpha * m_relax * m_deltaU); }
                catch (...) { failed = true; }

                if (failed)
                {
                    if (m_verbose==2) gsInfo<<"\t\tline search: alpha = "<<alpha<<" residual evaluation threw\n";
                }
                else if (!Rt.allFinite())
                {
                    if (m_verbose==2) gsInfo<<"\t\tline search: alpha = "<<alpha<<" non-finite residual\n";
                }
                else if (converged || Rt.norm() <= (1 - m_lsC*alpha)*R0n)
                {
                    accepted = true;
                    break;
                }
                else if (m_verbose==2)
                    gsInfo<<"\t\tline search: alpha = "<<alpha<<" |R_trial| = "<<Rt.norm()<<" > "<<(1 - m_lsC*alpha)*R0n<<"\n";
            }

            if (!accepted)
            {
                m_U += m_DeltaU;
                gsInfo<<"Line search failed: no step length in [rho^maxRed, 1] gives sufficient decrease. Solution did not converge\n";
                throw 1;
            }

            m_alpha = alpha;
            m_deltaU *= alpha;
            m_DeltaU += m_relax * m_deltaU;
            m_R = give(Rt);
        }
        m_residual = m_R.norm();

        if (m_verbose>0) { stepOutput(m_numIterations); }

        m_residualOld = m_residual;

        if (m_relax * m_deltaU.norm()/m_DeltaU.norm()  < m_tolU && m_residual/m_residualIni < m_tolF)
        {
            m_U+=m_DeltaU;
            break;
        }
        else if (m_numIterations+1 == m_maxIterations)
        {
            m_U+=m_DeltaU;
            gsInfo<<"Maximum iterations reached. Solution did not converge\n";
            throw 1;
        }
    }
    return m_U;
};

template <class T>
gsVector<T> gsStaticNewton<T>::_computeResidual(const gsVector<T> & U)
{
  gsVector<T> resVec;
  if (!m_residualFun(U, resVec))
    throw 2;
  return resVec;
}

template <class T>
gsSparseMatrix<T> gsStaticNewton<T>::_computeJacobian(const gsVector<T> & U, const gsVector<T> & deltaU)
{
  // Compute Jacobian
  gsSparseMatrix<T> m;
  if (!m_dnonlinear(U,deltaU,m))
    throw 2;
  return m;
}

template <class T>
void gsStaticNewton<T>::_factorizeMatrix(const gsSparseMatrix<T> & jacMat) const
{
    m_solver->compute(jacMat);
    if (m_solver->info()!=gsEigen::ComputationInfo::Success)
    {
      gsInfo<<"Solver error with code "<<m_solver->info()<<". See Eigen documentation on ComputationInfo \n"
                                                                    <<gsEigen::ComputationInfo::Success<<": Success"<<"\n"
                                                                    <<gsEigen::ComputationInfo::NumericalIssue<<": NumericalIssue"<<"\n"
                                                                    <<gsEigen::ComputationInfo::NoConvergence<<": NoConvergence"<<"\n"
                                                                    <<gsEigen::ComputationInfo::InvalidInput<<": InvalidInput"<<"\n";
      throw 3;
    }
}

template <class T>
gsVector<T> gsStaticNewton<T>::_solveSystem(const gsVector<T> & F)
{
    try
    {
      return m_solver->solve(F);
    }
    catch (...)
    {
      throw 3;
    }
}

template <class T>
void gsStaticNewton<T>::reset()
{
    m_dofs = m_force.rows();
    // resets m_U, m_DeltaU, m_deltaU, m_R, m_L, m_DeltaL, m_deltaL and m_headstart
    Base::reset();
}

template <class T>
void gsStaticNewton<T>::_init()
{
    this->reset();
    if( m_dnonlinear==nullptr || m_residualFun==nullptr)
        m_NL=false;
    else
        m_NL = true;

    m_stabilityMethod = 0;
    m_start = false;

    m_dofs = m_linear.rows();
    if (m_dofs==0)
        gsWarn<<"The number of degrees of freedom is equal to zero. This can lead to bad initialization.\n";

    m_residual = m_residualIni = m_residualOld = 0;
    m_alpha = 1;

    defaultOptions();

    m_status = gsStatus::NotStarted;
}

template <class T>
void gsStaticNewton<T>::_start()
{
    // Define residual measures:
    // residual    = current residual
    // residual0   = residual on m_U
    // residualOld = residual in previous step

    if (!m_headstart) // no headstart
    {
        // We can reset the update to ensure we properly restart
        if (m_dofs==0) gsWarn<<"The number of degrees of freedom is equal to zero. This can lead to bad initialization.\n";
        m_DeltaU.setZero(m_dofs);
        // Compute current residual and its norm
        m_R = this->_computeResidual(m_U);
        m_residual = m_R.norm();
        // If the residual is 0 (e.g. with purely displacment loading), we set it to 1 to allow divisions
        if (m_residual==0) m_residual=1;
        // All residual norms are equal
        m_residualIni = m_residualOld = m_residual;
    }
    else
    {
        // If we have a headstart, we need to compute Residual0 on the solution m_U
        // Compute current residual and its norm
        m_R = this->_computeResidual(m_U + m_DeltaU);
        m_residual = m_R.norm();
        // If the residual is 0 (e.g. with purely displacment loading), we set it to 1 to allow divisions
        if (m_residual==0) m_residual=1;
        // The previous step residual is the same as the residual
        m_residualOld = m_residual;
        // Residual0 is the residual without m_DeltaU
        m_residualIni = this->_computeResidual(m_U).norm();
        // If the residual is 0 (e.g. with purely displacment loading), we set it to 1 to allow divisions
        if (m_residualIni==0) m_residualIni=1;

        // Now we can reset the headstart
        m_headstart = false;
    }

    // Reset incremental update
    m_deltaU.setZero(m_dofs);

}

} // namespace gismo
