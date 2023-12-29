/** @file gsDynamicImplicitEuler.hpp

    @brief Performs the arc length method to solve a nonlinear equation system.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#pragma once

#include <gsSolver/gsGMRes.h>
#include <gsSolver/gsMatrixOp.h>
#include <gsSolver/gsBlockOp.h>

namespace gismo
{

template <class T>
template <bool _nonlinear>
typename std::enable_if<(_nonlinear==false), gsStatus>::type
gsDynamicImplicitEuler<T>::_step_impl(const T dt)
{
  m_Uold = m_U;
  m_Vold = m_V;

  index_t N = m_U.rows();
  gsVector<T> sol(2*N);
  sol.topRows(N) = m_U;
  sol.bottomRows(N) = m_V;

  typename gsBlockOp<>::Ptr Amat;

  Amat=gsBlockOp<>::make(2,2);
  gsGMRes<T> gmres(Amat);

  gsSparseMatrix<T> eye(N,N);
  eye.setIdentity();

  gsVector<T> F;
  gsSparseMatrix<T> M, C, K;

  this->_computeMass(m_time,M);
  this->_computeForce(m_time,F);
  this->_computeDamping(m_U,m_time,C);
  this->_computeJacobian(m_U,m_time,K);

  gsSparseSolver<>::LU solver(M);
  gsMatrix<T> MinvK = solver.solve(K);
  gsMatrix<T> MinvC = solver.solve(C);

  gsDebugVar(MinvK);
  gsDebugVar(MinvC);


  // top-left
  Amat->addOperator(0,0,makeMatrixOp( eye) );
  // top-right
  Amat->addOperator(0,1,makeMatrixOp( m_dt*eye  ) );
  // bottom-left
  Amat->addOperator(1,0,makeMatrixOp( m_dt*MinvK ) );
  // bottom-right
  Amat->addOperator(1,1,makeMatrixOp( eye+m_dt*MinvC ) );

  gsVector<T> MinvF = solver.solve(F);
  gsVector<T> rhs(2*N);
  rhs.setZero();

  gsDebugVar(m_dt*MinvF);
  rhs.bottomRows(N) = m_dt*MinvF;

  gsDebugVar(M.toDense().inverse()*F);
  gsDebugVar(MinvF);

  gsDebugVar(rhs);
  gsDebugVar(M.toDense());

/*
  P = M^-1 F
  MP = F
*/

  gsMatrix<T> tmpsol;
  gmres.solve(rhs+sol,tmpsol);

  gsDebugVar(tmpsol);

  m_U = tmpsol.topRows(N);
  m_V = tmpsol.bottomRows(N);
  m_time += m_dt;
  return gsStatus::Success;
}


template <class T>
template <bool _nonlinear>
typename std::enable_if<(_nonlinear==true), gsStatus>::type
gsDynamicImplicitEuler<T>::_step_impl(const T dt)
{
  m_Uold = m_U;
  m_Vold = m_V;

  index_t N = m_U.rows();
  gsVector<T> sol(2*N);
  sol.topRows(N) = m_U;
  sol.bottomRows(N) = m_V;

  gsMatrix<T> dsol;

  typename gsBlockOp<>::Ptr Amat;

  Amat=gsBlockOp<>::make(2,2);
  gsGMRes<T> gmres(Amat);

  gsSparseMatrix<T> eye(N,N);
  eye.setIdentity();

  gsVector<T> rhs(2*N);
  rhs.setZero();

  gsVector<T> R;
  gsSparseMatrix<T> M, C, K;

  this->_computeMass(m_time,M);
  this->_computeResidual(m_U,m_time,R);

  m_updateNorm   = 10*m_tolerance;
  m_residualNorm  = R.norm();
  T residualNorm0 = m_residualNorm;

  this->_initOutput();
  for (m_numIterations = 1; m_numIterations <= m_maxIterations; ++m_numIterations)
  {
      // Quasi newton
      // Quasi newton
      // Quasi newton
      // Quasi newton

      this->_computeDamping(m_U,m_time,C);
      this->_computeJacobian(m_U,m_time,K);

      Amat->addOperator(0,0,gsIdentityOp<T>::make(N) );
      Amat->addOperator(0,1,makeMatrixOp(-m_dt*eye) );
      Amat->addOperator(1,0,makeMatrixOp(m_dt*K) );
      Amat->addOperator(1,1,makeMatrixOp(M + m_dt*C) );

      rhs.topRows(m_U.rows()) = m_U - m_Uold + m_dt * m_V;
      rhs.bottomRows(m_V.rows()) = M*(m_Vold - m_V) + m_dt * C * m_V + m_dt*(-R);

      gmres.solve(-rhs,dsol);
      sol += dsol;
      m_updateNorm = dsol.norm() / sol.norm();

      m_U = sol.topRows(N);
      m_V = sol.bottomRows(N);

      this->_computeResidual(m_U,m_time,R);
      m_residualNorm = R.norm() / residualNorm0;

      this->_stepOutput(m_numIterations,m_residualNorm,m_updateNorm);

      if ( (m_updateNorm>m_tolerance && m_residualNorm>m_tolerance) )
      {
          m_converged = true;
          break;
      }
  }

  if (!m_converged)
  {
    gsInfo<<"maximum iterations reached. Solution did not converge\n";
    return gsStatus::NotConverged;
  }

  m_time += m_dt;
  return gsStatus::Success;

}

/*
template <class T>
gsStatus gsDynamicImplicitEuler<T>::_step(const T dt)
{
  m_Uold = m_U;
  m_Vold = m_V;

  index_t N = m_U.rows();
  gsVector<T> sol(2*N);
  sol.topRows(N) = m_U;
  sol.bottomRows(N) = m_V;

  gsMatrix<T> dsol;

  typename gsBlockOp<>::Ptr Amat;

  Amat=gsBlockOp<>::make(2,2);
  gsGMRes<T> gmres(Amat);

  gsSparseMatrix<T> eye(N,N);
  eye.setIdentity();

  gsVector<T> rhs(2*N);
  rhs.setZero();

  gsVector<T> R;
  gsSparseMatrix<T> M, C, K;

  this->_computeMass(m_time,M);
  this->_computeResidual(m_U,m_time,R);

  m_updateNorm   = 10*m_tolerance;
  m_residualNorm  = R.norm();
  T residualNorm0 = m_residualNorm;

  this->_initOutput();
  for (m_numIterations = 1; m_numIterations <= m_maxIterations; ++m_numIterations)
  {
      // Quasi newton
      // Quasi newton
      // Quasi newton
      // Quasi newton

      this->_computeDamping(m_U,m_time,C);
      this->_computeJacobian(m_U,m_time,K);

      Amat->addOperator(0,0,gsIdentityOp<T>::make(N) );
      Amat->addOperator(0,1,makeMatrixOp(-m_dt*eye) );
      Amat->addOperator(1,0,makeMatrixOp(m_dt*K) );
      Amat->addOperator(1,1,makeMatrixOp(M + m_dt*C) );

      rhs.topRows(m_U.rows()) = m_U - m_Uold + m_dt * m_V;
      rhs.bottomRows(m_V.rows()) = M*(m_Vold - m_V) + m_dt * C * m_V + m_dt*(-R);

      gmres.solve(-rhs,dsol);
      sol += dsol;
      m_updateNorm = dsol.norm() / sol.norm();

      m_U = sol.topRows(N);
      m_V = sol.bottomRows(N);

      this->_computeResidual(m_U,m_time,R);
      m_residualNorm = R.norm() / residualNorm0;

      this->_stepOutput(m_numIterations,m_residualNorm,m_updateNorm);

      if ( (m_updateNorm>m_tolerance && m_residualNorm>m_tolerance) )
      {
          m_converged = true;
          break;
      }
  }

  if (!m_converged)
  {
    gsInfo<<"maximum iterations reached. Solution did not converge\n";
    return gsStatus::NotConverged;
  }

  m_time += m_dt;
  return gsStatus::Success;

}
*/

template <class T>
void gsDynamicImplicitEuler<T>::_initOutput()
{
  gsInfo<<"\t";
  gsInfo<<std::setw(4)<<std::left<<"It.";
  gsInfo<<std::setw(17)<<std::left<<"|R|/|R0|";
  gsInfo<<std::setw(17)<<std::left<<"|dU|/|U0|";
}

template <class T>
void gsDynamicImplicitEuler<T>::_stepOutput(const index_t it, const index_t resnorm, const index_t updatenorm)
{
  gsInfo<<"\t";
  gsInfo<<std::setw(4)<<std::left<<it;
  gsInfo<<std::setw(17)<<std::left<<resnorm;
  gsInfo<<std::setw(17)<<std::left<<updatenorm;
}

} // namespace gismo