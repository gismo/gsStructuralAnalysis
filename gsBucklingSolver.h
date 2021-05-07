 /** @file gsBucklingSolver.h

    @brief Performs linear buckling analysis given a matrix or functions of a matrix

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
    @brief Performs the arc length method to solve a nonlinear equation system.

    \tparam T coefficient type

    \ingroup ThinShell
*/
template <class T, Spectra::GEigsMode GEigsMode = Spectra::GEigsMode::Cholesky>
class gsBucklingSolver
{
protected:
    // typedef typename std::vector<std::pair<T,gsMatrix<T>> > modes_t;

public:

  /// Constructor giving access to the gsShellAssembler object to create a linear system per iteration
  gsBucklingSolver(     gsSparseMatrix<T> &linear,
                        gsVector<T> &rhs,
                        std::function < gsSparseMatrix<T> ( gsVector<T> const & ) > &nonlinear,
                        T scaling = 1.0) :
    m_A(linear),
    m_rhs(rhs),
    m_nonlinearFun(nonlinear),
    m_scaling(scaling)
  {
    m_verbose = false;
    this->initializeMatrix();
  }

  /// Constructor giving access to the gsShellAssembler object to create a linear system per iteration
  gsBucklingSolver(     gsSparseMatrix<T> &linear,
                        gsSparseMatrix<T> &nonlinear ) :
    m_A(linear)
  {
    m_verbose = false;
    m_B = nonlinear-m_A;
  }
public:

    void verbose() {m_verbose=true; };

    void compute(T shift = 0.0);
    void computeSparse( T shift = 0.0,
                        index_t number = 10,
                        index_t ncvFac = 3,
                        Spectra::SortRule selectionRule = Spectra::SortRule::SmallestMagn,
                        Spectra::SortRule sortRule = Spectra::SortRule::SmallestMagn)
    {computeSparse_impl<GEigsMode>(shift,number,ncvFac,selectionRule,sortRule);};

    void computePower();

    gsMatrix<T> values() const { return m_values; };
    T value(int k) const { return m_values.at(k); };

    gsMatrix<T> vectors() const { return m_vectors; };
    gsMatrix<T> vector(int k) const { return m_vectors.col(k); };

    std::vector<std::pair<T,gsMatrix<T>> > mode(int k) const {return makeMode(k); }

private:
    template<Spectra::GEigsMode _GEigsMode>
    typename std::enable_if<_GEigsMode==Spectra::GEigsMode::Cholesky ||
                            _GEigsMode==Spectra::GEigsMode::RegularInverse
                            ,
                            void>::type computeSparse_impl(T shift, index_t number, index_t ncvFac, Spectra::SortRule selectionRule, Spectra::SortRule sortRule);

    template<Spectra::GEigsMode _GEigsMode>
    typename std::enable_if<_GEigsMode==Spectra::GEigsMode::ShiftInvert ||
                            _GEigsMode==Spectra::GEigsMode::Buckling ||
                            _GEigsMode==Spectra::GEigsMode::Cayley
                            ,
                            void>::type computeSparse_impl(T shift, index_t number, index_t ncvFac, Spectra::SortRule selectionRule, Spectra::SortRule sortRule);


protected:

    const gsSparseMatrix<T> m_A;
    gsSparseMatrix<T> m_B;
    const gsVector<T> m_rhs;
    const std::function < gsSparseMatrix<T> ( gsVector<T> const & ) > m_nonlinearFun;
    T m_scaling;

    /// Linear solver employed
    gsSparseSolver<>::SimplicialLDLT  m_solver;
    Eigen::GeneralizedSelfAdjointEigenSolver< gsMatrix<real_t>::Base >  m_eigSolver;

    gsSparseMatrix<T> m_diff;
    gsVector<T> m_solVec;
    gsMatrix<T> m_values,m_vectors;

    bool m_verbose;

    index_t m_num;

protected:

    void initializeMatrix();
    std::vector<std::pair<T,gsMatrix<T>> > makeMode(int k) const;

};


} // namespace gismo


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsBucklingSolver.hpp)
#endif