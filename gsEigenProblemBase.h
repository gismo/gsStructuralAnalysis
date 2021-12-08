 /** @file gsEigenProblemBase.h

    @brief Base class for buckling and modal analyses

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
    @brief Base class for buckling and modal analyses

    \tparam T coefficient type

    \ingroup gsStructuralAnalysis
*/
template <class T, Spectra::GEigsMode GEigsMode>
class gsEigenProblemBase
{

public:

    ~gsEigenProblemBase() {};

public:

    virtual void verbose() {m_verbose=true; };

    virtual void compute();
    virtual void compute(T shift);

    virtual void computeSparse(T shift = 0.0, index_t number = 10, index_t ncvFac = 3, Spectra::SortRule selectionRule = Spectra::SortRule::SmallestMagn, Spectra::SortRule sortRule = Spectra::SortRule::SmallestMagn)
    {computeSparse_impl<GEigsMode>(shift,number,ncvFac,selectionRule,sortRule);};

    virtual void computePower();

    virtual gsMatrix<T> values() const { return m_values; };
    virtual T value(int k) const { return m_values.at(k); };

    virtual gsMatrix<T> vectors() const { return m_vectors; };
    virtual gsMatrix<T> vector(int k) const { return m_vectors.col(k); };

    virtual std::vector<std::pair<T,gsMatrix<T>> > mode(int k) const {return makeMode(k); }



protected:

    virtual std::vector<std::pair<T,gsMatrix<T>> > makeMode(int k) const;

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

    gsSparseMatrix<T> m_A;
    gsSparseMatrix<T> m_B;

    Eigen::GeneralizedSelfAdjointEigenSolver< gsMatrix<real_t>::Base >  m_eigSolver;

    gsMatrix<T> m_values,m_vectors;

    bool m_verbose;

    index_t m_num;

};


} // namespace gismo


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsEigenProblemBase.hpp)
#endif