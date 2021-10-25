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

template<index_t dim, class T>
class gsSpaceTimeFitter
{
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

  gsTensorBSplineBasis<dim,T> _basis(const gsKnotVector<T> &kv, index_t nsteps)
  {
    return _basis_impl<dim>(kv,nsteps);
  }

  template<index_t _dim>
  typename std::enable_if<_dim==1, gsTensorBSplineBasis<_dim,T>>::type
  _basis_impl(const gsKnotVector<T> &kv, index_t nsteps)
  {
    gsTensorBSplineBasis<_dim,T> tbasis(
                                            kv
                                            );
    return tbasis;
  }

  template<index_t _dim>
  typename std::enable_if<_dim==2, gsTensorBSplineBasis<_dim,T>>::type
  _basis_impl(const gsKnotVector<T> &kv, index_t nsteps)
  {
    gsTensorBSplineBasis<_dim,T> tbasis(
                                            static_cast<gsBSplineBasis<T> *>(&m_bases.basis(0).component(1))->knots(),
                                            kv
                                            );
    return tbasis;
  }

  template<index_t _dim>
  typename std::enable_if<_dim==3, gsTensorBSplineBasis<_dim,T>>::type
  _basis_impl(const gsKnotVector<T> &kv, index_t nsteps)
  {
    gsTensorBSplineBasis<_dim,T> tbasis(
                                            static_cast<gsBSplineBasis<T> *>(&m_bases.basis(0).component(0))->knots(),
                                            static_cast<gsBSplineBasis<T> *>(&m_bases.basis(0).component(1))->knots(),
                                            kv
                                            );
    return tbasis;
  }

  template<index_t _dim>
  typename std::enable_if<_dim==4, gsTensorBSplineBasis<_dim,T>>::type
  _basis_impl(const gsKnotVector<T> &kv, index_t nsteps)
  {
    gsTensorBSplineBasis<_dim,T> tbasis(
                                            static_cast<gsBSplineBasis<T> *>(&m_bases.basis(0).component(0))->knots(),
                                            static_cast<gsBSplineBasis<T> *>(&m_bases.basis(0).component(1))->knots(),
                                            static_cast<gsBSplineBasis<T> *>(&m_bases.basis(0).component(2))->knots(),
                                            kv
                                            );
    return tbasis;
  }

public:

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

    // for (index_t p = 0; p!=dbasis.nBases(); ++p)
    // {
      // gsTensorBSplineBasis<dim,T> tbasis(
      //                                         static_cast<gsBSplineBasis<T> *>(&m_bases.basis(0).component(0))->knots(),
      //                                         static_cast<gsBSplineBasis<T> *>(&m_bases.basis(0).component(1))->knots(),
      //                                         kv
      //                                         );
    // }

    m_basis = _basis(kv,nsteps);

    gsMatrix<> rhs(m_times.size(),(dim+1)*bsize);
    gsVector<> ones; ones.setOnes(bsize);

    for (index_t lam = 0; lam!=nsteps; ++lam)
    {
      rhs.block(lam,0,1,dim * bsize) = m_data.at(lam).reshape(1,dim * bsize);
      rhs.block(lam,dim*bsize,1,bsize) = m_times.at(lam) * ones.transpose();
    }

    // get the Greville Abcissae (anchors)
    gsMatrix<> anchors = lbasis.anchors();

    // Get the collocation matrix at the anchors
    gsSparseMatrix<> C;
    lbasis.collocationMatrix(anchors,C);

    gsSparseSolver<>::LU solver;
    solver.compute(C);

    gsMatrix<> sol;
    m_coefs.resize((nsteps)*bsize,dim+1);

    sol = solver.solve(rhs);

    for (index_t lam = 0; lam!=nsteps; ++lam)
    {
      gsMatrix<> tmp = sol.block(lam,0,1,dim * bsize);
      m_coefs.block(lam * bsize,0,bsize,dim) = tmp.reshape(bsize,dim);
      m_coefs.block(lam * bsize,dim,bsize,1) = sol.block(lam,dim*bsize,1,bsize).transpose();
    }

    // gsTensorBSpline<3,T> tspline = tbasis.makeGeometry(give(coefs)).release();
    // gsTensorBSpline<dim,T> tspline(tbasis,give(coefs));




    // return tbasis.makeGeometry(give(coefs));
  }

protected:
  std::vector<gsMatrix<T>> m_data;
  gsVector<T> m_times;
  gsVector<T> m_ptimes;
  gsMultiBasis<T> m_bases;
  index_t m_deg;

  mutable gsTensorBSplineBasis<dim,T> m_basis;
  gsMatrix<T> m_coefs;

};


} // namespace gismo

// #ifndef GISMO_BUILD_LIB
// #include GISMO_HPP_HEADER(gsSpaceTimeHierarchy.hpp)
// #endif
