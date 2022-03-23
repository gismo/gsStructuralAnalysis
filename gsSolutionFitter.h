 /** @file gsSolutionFitter.h

    @brief XXX

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#pragma once

namespace gismo
{

template<class T>
class gsSolutionFitter
{
public:
  gsSolutionFitter  ( const std::vector<gsMatrix<T>> & solutionCoefs,
                      const gsVector<T> & times,
                      const gsVector<T> & ptimes,
                      const index_t deg = 2)
  :
  m_data(solutionCoefs),
  m_times(times),
  m_ptimes(ptimes),
  m_deg(deg)
  {
    GISMO_ENSURE(m_data.size() != 0,"No data provided!");
    GISMO_ENSURE(m_data.size() == (size_t)m_times.size(),"Solution coefs and times should have the same size");
    GISMO_ENSURE((size_t)m_ptimes.size() == (size_t)m_times.size(),"(Parametric)times should have the same size");
    m_targetDim = m_data.at(0).cols();
  }

  void setDegree(index_t deg) {m_deg = deg;}

  void addDataPoint(gsMatrix<T> & solution, T time, T ptime, index_t continuity )
  {
    index_t N = m_ptimes.rows();
    m_ptimes.conservativeResize(N+1);
    m_ptimes.row(N) = ptime;
    gsDebugVar(m_ptimes);
    std::sort(begin(m_ptimes),end(m_ptimes));
    gsDebugVar(m_ptimes);
  }

  const gsVector<T> & ptimes() const { return m_ptimes; }
  const gsVector<T> &  times() const { return  m_times; }
  const std::vector<gsMatrix<T>> & data() const { return  m_data; }

protected:

  gsBSplineBasis<T> _basis(const gsKnotVector<T> &kv, index_t nsteps)
  {
    gsBSplineBasis<T> tbasis(kv);
    return tbasis;
  }

public:

  std::pair<T,gsMatrix<T> > slice(T xi)
  {
    T load;
    gsVector<T> pt(1);
    pt<<xi;
    gsMatrix<T> res = m_fit.eval(pt);
    load = res(res.rows()-1,0); // the last entry is the load factor
    res.conservativeResize(res.rows()-1,1);
    index_t bsize = m_data.at(0).rows();
    res.resize(bsize,m_targetDim);
    return std::make_pair(load,res);
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
    GISMO_ASSERT(m_data.size()==(size_t)m_times.rows(),"Number of time and solution steps should match! "<<m_data.size()<<"!="<<m_times.rows());
    // GISMO_ASSERT(m_data.at(0).cols() == dim,"Is the dimension correct?"<<m_data.at(0).cols()<<"!="<<dim);
    index_t nsteps = m_times.rows();
    index_t bsize = m_data.at(0).rows();

    // Prepare fitting basis
    gsKnotVector<> kv(m_ptimes.minCoeff(),m_ptimes.maxCoeff(),nsteps-(m_deg+1),m_deg+1);
    m_basis = gsBSplineBasis<T>(kv);

    gsDebugVar(m_basis.size());
    gsDebugVar(kv);

    gsMatrix<T> rhs(m_times.size(),m_targetDim*bsize+1);
    gsVector<T> ones; ones.setOnes(bsize);

    for (index_t lam = 0; lam!=nsteps; ++lam)
    {
      rhs.block(lam,0,1,m_targetDim * bsize) = m_data.at(lam).reshape(1,m_targetDim * bsize);
      rhs(lam,m_targetDim*bsize) = m_times.at(lam);
    }

    // get the Greville Abcissae (anchors)
    gsMatrix<T> anchors = m_basis.anchors();

    // Get the collocation matrix at the anchors
    gsSparseMatrix<T> C = m_basis.collocationMatrix(anchors);

    gsSparseSolver<>::LU solver;
    solver.compute(C);

    m_coefs = solver.solve(rhs);

    m_fit = gsBSpline<T>(m_basis,give(m_coefs));
  }

protected:
  index_t m_targetDim;
  std::vector<gsMatrix<T>> m_data;
  gsVector<T> m_times;
  gsVector<T> m_ptimes;
  index_t m_deg;

  mutable gsBSplineBasis<T> m_basis;
  gsMatrix<T> m_coefs;

  mutable gsBSpline<T> m_fit;

};


} // namespace gismo

// #ifndef GISMO_BUILD_LIB
// #include GISMO_HPP_HEADER(gsSpaceTimeHierarchy.hpp)
// #endif
