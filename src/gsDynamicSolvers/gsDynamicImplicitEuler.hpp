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

template <class T, bool _NL>
template <bool _nonlinear>
typename std::enable_if<(_nonlinear==false), gsStatus>::type
gsDynamicImplicitEuler<T,_NL>::_step_impl(const T t, const T dt, gsVector<T> & U, gsVector<T> & V, gsVector<T> & A) const
{
  gsVector<T> Uold = U;
  gsVector<T> Vold = V;
  gsVector<T> Aold = A;

  index_t N = U.rows();
  gsVector<T> sol(2*N);
  sol.topRows(N) = U;
  sol.bottomRows(N) = V;

  typename gsBlockOp<>::Ptr Amat;

  Amat=gsBlockOp<>::make(2,2);
  gsGMRes<T> gmres(Amat); 

  gsSparseMatrix<T> eye(N,N);
  eye.setIdentity();

  gsVector<T> F;
  gsSparseMatrix<T> M, C, K;
  gsSparseMatrix<T> Minv;

  // Computed at t=t0+dt
  this->_computeMass(t+dt,M);
  this->_computeForce(t+dt,F);
  this->_computeDamping(U,t+dt,C);
  this->_computeJacobian(U,t,K);

  // top-left
  Amat->addOperator(0,0,gsIdentityOp<T>::make(N) );
  // top-right
  Amat->addOperator(0,1,makeMatrixOp( -dt*eye  ) );
  // bottom-left
  Amat->addOperator(1,0,makeMatrixOp( dt*K ) );
  // bottom-right
  Amat->addOperator(1,1,makeMatrixOp( M + dt*C ) );

  gsVector<T> rhs(2*N);
  rhs.topRows(N) = U;
  rhs.bottomRows(N) = dt*F + M*V;

  gsMatrix<T> tmpsol;
  gmres.solve(rhs,tmpsol);

  U = tmpsol.topRows(N);
  V = tmpsol.bottomRows(N);
  this->_initOutput();
  this->_stepOutput(0,sol.norm(),0.);

  return gsStatus::Success;
}


template <class T, bool _NL>
template <bool _nonlinear>
typename std::enable_if<(_nonlinear==true), gsStatus>::type
gsDynamicImplicitEuler<T,_NL>::_step_impl(const T t, const T dt, gsVector<T> & U, gsVector<T> & V, gsVector<T> & A) const
{
  gsVector<T> Uold = U;
  gsVector<T> Vold = V;
  gsVector<T> Aold = A;

  index_t N = U.rows();
  gsVector<T> sol(2*N);
  sol.topRows(N) = U;
  sol.bottomRows(N) = V;

  gsMatrix<T> dsol;

  typename gsBlockOp<>::Ptr Amat;

  Amat=gsBlockOp<>::make(2,2);
  gsGMRes<T> gmres(Amat);

  gsSparseMatrix<T> eye(N,N);
  eye.setIdentity();

  gsVector<T> R;
  gsSparseMatrix<T> M, C, K;

  // Computed at t=t0+dt
  this->_computeMass(t+dt,M);
  this->_computeDamping(U,t+dt,C);
  this->_computeResidual(U,t+dt,R);

  gsVector<T> rhs(2*N);
  rhs.topRows(N) = - dt * V; // same as below, but U-Uold=0
  rhs.bottomRows(N) = dt * C * V + dt*(-R); // same as below, but V-Vold=0

  T tolU = m_options.getReal("TolU");
  T tolF = m_options.getReal("TolF");
  T updateNorm   = 10.0*tolU;
  T residualNorm  = rhs.norm();
  T residualNorm0 = (residualNorm!=0) ? residualNorm : 1;
  T solnorm, dsolnorm;
  this->_initOutput();
  for (index_t numIterations = 0; numIterations < m_options.getInt("MaxIter"); ++numIterations)
  {
      // TODO: Quasi newton
    if ((!m_options.getSwitch("Quasi")) || ((numIterations==0) || (numIterations % m_options.getInt("QuasiIterations") == 0)))
    {
      // Computed at t=t0+dt
      this->_computeDamping(U,t+dt,C);
      this->_computeJacobian(U,t+dt,K);
    }
      
    Amat->addOperator(0,0,gsIdentityOp<T>::make(N) );
    Amat->addOperator(0,1,makeMatrixOp(-dt*eye) );
    Amat->addOperator(1,0,makeMatrixOp(dt*K) );
    Amat->addOperator(1,1,makeMatrixOp(M + dt*C) );

    rhs.topRows(N) = U - Uold - dt * V;
    rhs.bottomRows(N) = M*(V - Vold) + dt * C * V + dt*(-R);

    gmres.solve(-rhs,dsol);
    sol += dsol;

    solnorm = sol.norm();
    dsolnorm = dsol.norm();
    updateNorm = (solnorm != 0) ? dsolnorm / solnorm : dsolnorm;

    U = sol.topRows(N);
    V = sol.bottomRows(N);

    this->_computeResidual(U,t,R);
    residualNorm = rhs.norm() / residualNorm0;

    this->_stepOutput(numIterations,residualNorm,updateNorm);

    if ( (updateNorm<tolU && residualNorm<tolF) )
    {
        return gsStatus::Success;
    }
  }

  gsInfo<<"maximum iterations reached. Solution did not converge\n";
  return gsStatus::NotConverged;
}

template <class T, bool _NL>
gsStatus gsDynamicImplicitEuler<T,_NL>::_step(const T t, const T dt,
                                        gsVector<T> & U, gsVector<T> & V,
                                        gsVector<T> & A) const
{
    gsStatus status = gsStatus::NotStarted;
    status = _step_impl<_NL>(t,dt,U,V,A);
    return status;
}

template <class T, bool _NL>
void gsDynamicImplicitEuler<T,_NL>::_initOutput() const
{
  if (m_options.getSwitch("Verbose"))
  {
    gsInfo<<"\t";
    gsInfo<<std::setw(4)<<std::left<<"It.";
    gsInfo<<std::setw(17)<<std::left<<"|R|/|R0|";
    gsInfo<<std::setw(17)<<std::left<<"|dU|/|U0|"<<"\n";
  }
}

template <class T, bool _NL>
void gsDynamicImplicitEuler<T,_NL>::_stepOutput(const index_t it, const T resnorm, const T updatenorm) const
{
  if (m_options.getSwitch("Verbose"))
  {
    gsInfo<<"\t";
    gsInfo<<std::setw(4)<<std::left<<it;
    gsInfo<<std::setw(17)<<std::left<<resnorm;
    gsInfo<<std::setw(17)<<std::left<<updatenorm<<"\n";
  }
}

} // namespace gismo