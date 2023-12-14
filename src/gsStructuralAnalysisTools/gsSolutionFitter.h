 /** @file gsSolutionFitter.h

    @brief XXX

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#include <gsModeling/gsFitting.h>
#include <gsIO/gsOptionList.h>

#pragma once

namespace gismo
{

template<class T>
class gsSolutionFitter
{
public:
  gsSolutionFitter  ( const std::vector<gsVector<T>> & data,
                      const gsVector<T> & times,
                      const index_t deg = 2)
  :
  m_deg(deg)
  {
    this->_defaultOptions();

    GISMO_ENSURE(data.size() != 0,"No data provided!");
    GISMO_ENSURE(data.size() == (size_t)times.size(),"Solution coefs and times should have the same size");
    // GISMO_ENSURE((size_t)m_ptimes.size() == (size_t)m_times.size(),"(Parametric)times should have the same size");

    for (index_t k=0; k!=times.size(); k++)
      m_data.insert({times[k],data[k]});
  }

  // gsSolutionFitter  ( const std::vector<gsVector<T>> & data,
  //                     const gsVector<T> & times,
  //                     const gsVector<T> & ptimes,
  //                     const index_t deg = 2)
  // :
  // {
  //   std::vector<gsVector<T>> tmp = data;
  //   GISMO_ENSURE(data.size() == (size_t)times.rows(),"Solution coefs and times should have the same size");
  //   for (size_t k = 0; k!=data.size(); k++)
  //   {
  //     tmp[k] = conservativeResize(tmp[k].rows()+1);
  //     tmp[k].tail() = times[k];
  //   }

  //   gsSolutionFitter<T>(tmp,ptimes,deg);
  // }

  void setDegree(index_t deg) {m_deg = deg;}

  void addDataPoint(gsVector<T> & point, T ptime)//, index_t continuity )
  {
    m_data.insert({ptime,point});
    gsDebugVar(m_data.size());
  }

  gsOptionList& options() {return m_options;}

protected:

  void _defaultOptions()
  {
    m_options.addInt("Parameterization","Type of parameterization to be used: 0: uniform, 1: chord-length, 2: centripetal",1);
  }

  gsBSplineBasis<T> _basis(const gsKnotVector<T> &kv, index_t nsteps)
  {
    gsBSplineBasis<T> tbasis(kv);
    return tbasis;
  }

public:

  const gsBSpline<T> & fit() const {return m_fit; }

  gsVector<T> slice(T xi)
  {
    gsVector<T> pt(1);
    pt<<xi;
    return m_fit.eval(pt);
  }

  T nearestParam(const gsMatrix<T> & solution, gsMatrix<T> & params, index_t nTrialPoints = 25, index_t nIterations = 100)
  {
      gsMatrix<T> supp = m_fit.support();
      gsVector<T> start = supp.col(0), end = supp.col(1);
      gsMatrix<T> trialPoints = uniformPointGrid(start, end, nTrialPoints);
      gsMatrix<T> surfVal, surfDeriv, surfDeriv2;
      gsMatrix<T> closestParam;
      gsMatrix<T> u, du;
      T tol = 1e-8;
      T err = 1;
      T closestSqDist(10e100);

      for(int idxTrial = 0; idxTrial < nTrialPoints; idxTrial++)
      {
        u = trialPoints.col(idxTrial);
        // apply Newton's method to the function
        // f(u) = ||p - x(u)||^2
        // where x(u) is the parametrisation of the curve in space.
        // (todo - also check the distances at the endpoints of the curve?
        // although this is for splitting the curve and we would never want
        // to split at the endpoint.)
        for(int iteration = 0; iteration < nIterations; iteration++)
        {
            m_fit.eval_into(u, surfVal);
            m_fit.jacobian_into(u, surfDeriv);
            m_fit.deriv2_into(u, surfDeriv2);

            // evaluate derivative of f
            gsMatrix<T> sqDistDeriv = -2 * (solution - surfVal).transpose() *
                    surfDeriv;
            GISMO_ASSERT(sqDistDeriv.rows() == 1 && sqDistDeriv.cols() == 1, "Derivative should be 1x1");

            gsMatrix<T> sqDistDeriv2 = 2*surfDeriv.transpose() * surfDeriv - 2*(solution - surfVal).transpose() * surfDeriv2;
            GISMO_ASSERT(sqDistDeriv2.rows() == 1 && sqDistDeriv2.cols() == 1, "Second derivative should be 1x1");

            du = sqDistDeriv / sqDistDeriv2(0, 0);
            u -= du;

            if (du.norm()/u.norm() < tol)
              break;

            u(0, 0) = (u(0, 0) < supp(0, 0))? supp(0, 0): ((u(0, 0) > supp(0, 1)? supp(0, 1): u(0, 0)));
        }
        // compute sqDist for the point found by the last iteration, and compare against the best seen so far
        m_fit.eval_into(u, surfVal);
        T sqDist = (solution - surfVal).squaredNorm();
        if(idxTrial == 0 || sqDist < closestSqDist)
        {
            closestParam = u;
            closestSqDist = sqDist;
        }
    }
    return closestParam.value();
  }

  void compute()
  {
    GISMO_ASSERT(m_data.size()>(m_deg),"Data needs to have more than "<<m_deg<<" points, but only "<<m_data.size()<<" are stored!");
    index_t nsteps = m_data.size();

    gsDebugVar(nsteps);

    // gsMatrix<T> ptimes(1,nsteps);
    std::vector<T> ptimes(nsteps);
    gsMatrix<T> rhs(nsteps,m_data.at(0).rows());
    index_t k=0;
    for(typename std::map<T,gsVector<T>>::iterator it = m_data.begin(); it != m_data.end(); ++it, ++k)
    {
      // ptimes(0,k) = it->first;
      ptimes.at(k) = it->first;
      rhs.row(k) = it->second;
    }

    // Prepare fitting basis
    std::vector<T> u_bar(ptimes.size());
    u_bar.front() = 0;
    u_bar.back()  = 1;

    gsKnotVector<T> kv;
    index_t param = m_options.getInt("Parameterization");
    if (param==0)
    {
      kv = gsKnotVector<T>(gsAsVector<T>(ptimes).minCoeff(),gsAsVector<T>(ptimes).maxCoeff(),nsteps-(m_deg+1),m_deg+1, 1);
    }
    else if (param==1 || param==2)
    {
      // Chord-length parameterization
      T d = 0.;
      for (typename std::vector<T>::iterator pit = std::next(ptimes.begin()); pit!=ptimes.end(); pit++)
        d += math::pow(*pit - *std::prev(pit),1./param);

      typename std::vector<T>::iterator pit = std::next(ptimes.begin());
      for (typename std::vector<T>::iterator it=std::next(u_bar.begin()); it!=std::prev(u_bar.end()); it++, pit++)
        *it = *std::prev(it) + 1. / d * math::pow(*pit - *std::prev(pit),1./param);

      std::vector<T> u(nsteps+m_deg+1);
      for (index_t k=0; k!=m_deg+1; k++)
      {
        u.at(k) = 0;
        u.at(u.size()-k-1) = 1;
      }
      for (index_t j=1; j!=nsteps-m_deg; j++)
      {
        u.at(j+m_deg) = 0;
        for (index_t i=j; i!=j+m_deg; i++)
          u.at(j+m_deg) += u_bar.at(i);

        u.at(j+m_deg) /= m_deg;
      }

      kv = gsKnotVector<T>(u,m_deg);
      kv.affineTransformTo(ptimes.front(),ptimes.back());
    }
    else
      GISMO_ERROR("Parameterization unknown");


    m_basis = gsBSplineBasis<T>(kv);

    // get the Greville Abcissae (anchors)
    gsMatrix<T> anchors = m_basis.anchors();

    // Get the collocation matrix at the anchors
    gsSparseMatrix<T> C = m_basis.collocationMatrix(anchors);

    // gsSparseSolver<>::LU solver;
    // solver.compute(C);

    // m_coefs = solver.solve(rhs);
    // m_fit = gsBSpline<T>(m_basis,give(m_coefs));

    gsFitting<T> fitter(gsAsMatrix<T>(ptimes,1,ptimes.size()),rhs.transpose(),m_basis);
    gsSparseMatrix<T> lhs(nsteps,nsteps);
    lhs.setIdentity();
    fitter.setConstraints(C,rhs);
    fitter.compute(0.0);

    m_fit = gsBSpline<T>(m_basis,give(fitter.result()->coefs()));
  }

protected:
  std::map<T,gsVector<T>> m_data;

  gsVector<T> m_ptimes;
  index_t m_deg;

  mutable gsBSplineBasis<T> m_basis;
  gsMatrix<T> m_coefs;

  mutable gsBSpline<T> m_fit;

  gsOptionList m_options;

};


} // namespace gismo

// #ifndef GISMO_BUILD_LIB
// #include GISMO_HPP_HEADER(gsSpaceTimeHierarchy.hpp)
// #endif
