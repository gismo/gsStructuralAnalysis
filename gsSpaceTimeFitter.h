 /** @file gsSpaceTimeHierarchy.h

    @brief Performs the arc length method to solve a nonlinear equation system.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#pragma once

namespace gismo
{

template<size_t domainDim, size_t targetDim, class T>
class gsSpaceTimeFitter
{
    typedef typename gsTensorBSpline<domainDim+1,real_t>::BoundaryGeometryType slice_t;


public:
  gsSpaceTimeFitter  ( const std::vector<gsMatrix<T>> & solutionCoefs,
                      const gsVector<T> & times,
                      const gsVector<T> & ptimes,
                      const gsMultiBasis<T> & spatialBasis,
                      const index_t deg = 2)
  :
  m_data(solutionCoefs),
  m_times(times),
  m_ptimes(ptimes),
  m_bases(spatialBasis),
  m_deg(deg)
  {
    GISMO_ENSURE(m_data.size() != 0,"No data provided!");
    GISMO_ENSURE(m_data.size() == m_times.size(),"Solution coefs and times should have the same size");
    GISMO_ENSURE(m_ptimes.size() == m_times.size(),"(Parametric)times should have the same size");
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

protected:

  gsTensorBSplineBasis<domainDim+1,T> _basis(const gsKnotVector<T> &kv, index_t nsteps)
  {
    return _basis_impl<domainDim+1>(kv,nsteps);
  }

  template<index_t _domainDim>
  typename std::enable_if<_domainDim==2, gsTensorBSplineBasis<_domainDim,T>>::type
  _basis_impl(const gsKnotVector<T> &kv, index_t nsteps)
  {
    gsTensorBSplineBasis<_domainDim,T> tbasis(
                                            static_cast<gsBSplineBasis<T> *>(&m_bases.basis(0).component(1))->knots(),
                                            kv
                                            );
    return tbasis;
  }

  template<index_t _domainDim>
  typename std::enable_if<_domainDim==3, gsTensorBSplineBasis<_domainDim,T>>::type
  _basis_impl(const gsKnotVector<T> &kv, index_t nsteps)
  {
    gsTensorBSplineBasis<_domainDim,T> tbasis(
                                            static_cast<gsBSplineBasis<T> *>(&m_bases.basis(0).component(0))->knots(),
                                            static_cast<gsBSplineBasis<T> *>(&m_bases.basis(0).component(1))->knots(),
                                            kv
                                            );
    return tbasis;
  }

  template<index_t _domainDim>
  typename std::enable_if<_domainDim==4, gsTensorBSplineBasis<_domainDim,T>>::type
  _basis_impl(const gsKnotVector<T> &kv, index_t nsteps)
  {
    gsTensorBSplineBasis<_domainDim,T> tbasis(
                                            static_cast<gsBSplineBasis<T> *>(&m_bases.basis(0).component(0))->knots(),
                                            static_cast<gsBSplineBasis<T> *>(&m_bases.basis(0).component(1))->knots(),
                                            static_cast<gsBSplineBasis<T> *>(&m_bases.basis(0).component(2))->knots(),
                                            kv
                                            );
    return tbasis;
  }

public:

  std::pair<T,gsGeometry<T> *> slice(T xi)
  {
    slice_t res;
    m_fit.slice(domainDim,xi,res);
    gsGeometry<T> * geom = res.clone().release();
    T load = geom->coefs()(0,domainDim+1);
    geom->embed(targetDim);
    return std::make_pair(load,geom);
  }

  void compute()
  {
    GISMO_ASSERT(m_data.size()==m_times.rows(),"Number of time and solution steps should match! "<<m_data.size()<<"!="<<m_times.rows());
    // GISMO_ASSERT(m_data.at(0).cols() == dim,"Is the dimension correct?"<<m_data.at(0).cols()<<"!="<<dim);
    index_t nsteps = m_times.rows();
    index_t bsize = m_data.at(0).rows();

    // Prepare fitting basis
    gsKnotVector<> kv(m_ptimes.minCoeff(),m_ptimes.maxCoeff(),nsteps-(m_deg+1),m_deg+1);
    gsBSplineBasis<T> lbasis(kv);

    //////// TO DO:
    //////// - Make multi-patch compatible
    //////// - Include dimension d
    //////// - Make compatible with other basis types (dynamic casts)

    m_basis = _basis(kv,nsteps);

    gsMatrix<> rhs(m_times.size(),(targetDim+1)*bsize);
    gsVector<> ones; ones.setOnes(bsize);

    for (index_t lam = 0; lam!=nsteps; ++lam)
    {
      rhs.block(lam,0,1,targetDim * bsize) = m_data.at(lam).reshape(1,targetDim * bsize);
      rhs.block(lam,targetDim*bsize,1,bsize) = m_times.at(lam) * ones.transpose();
    }

    // get the Greville Abcissae (anchors)
    gsMatrix<> anchors = lbasis.anchors();

    // Get the collocation matrix at the anchors
    gsSparseMatrix<> C;
    lbasis.collocationMatrix(anchors,C);

    gsSparseSolver<>::LU solver;
    solver.compute(C);

    gsMatrix<> sol;
    m_coefs.resize((nsteps)*bsize,targetDim+1);

    sol = solver.solve(rhs);

    for (index_t lam = 0; lam!=nsteps; ++lam)
    {
      gsMatrix<> tmp = sol.block(lam,0,1,targetDim * bsize);
      m_coefs.block(lam * bsize,0,bsize,targetDim) = tmp.reshape(bsize,targetDim);
      m_coefs.block(lam * bsize,targetDim,bsize,1) = sol.block(lam,targetDim*bsize,1,bsize).transpose();
    }

    // gsTensorBSpline<3,T> tspline = tbasis.makeGeometry(give(coefs)).release();
    // gsTensorBSpline<dim,T> tspline(m_basis,give(m_coefs));




    m_fit = gsTensorBSpline<domainDim+1,T>(m_basis,give(m_coefs));
  }

protected:
  std::vector<gsMatrix<T>> m_data;
  gsVector<T> m_times;
  gsVector<T> m_ptimes;
  gsMultiBasis<T> m_bases;
  index_t m_deg;

  mutable gsTensorBSplineBasis<domainDim+1,T> m_basis;
  gsMatrix<T> m_coefs;

  mutable gsTensorBSpline<domainDim+1,T> m_fit;

};


} // namespace gismo

// #ifndef GISMO_BUILD_LIB
// #include GISMO_HPP_HEADER(gsSpaceTimeHierarchy.hpp)
// #endif
