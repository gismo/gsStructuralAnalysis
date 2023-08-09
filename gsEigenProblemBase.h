 /** @file gsEigenProblemBase.h

    @brief Base class for buckling and modal analyses

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#include <typeinfo>

#include <gsCore/gsLinearAlgebra.h>
#ifdef GISMO_WITH_SPECTRA
#include <gsSpectra/gsSpectra.h>
#endif
#include <gsIO/gsOptionList.h>
#include <gsStructuralAnalysis/gsStructuralAnalysisTools.h>

#pragma once


namespace gismo
{

/**
    @brief Base class for buckling and modal analyses

    \tparam T coefficient type

    \ingroup gsStructuralAnalysis
*/
template <class T>
class gsEigenProblemBase
{

public:

    gsEigenProblemBase() 
    { 
        m_options = defaultOptions(); 
        m_status  = gsStatus::NotStarted;
    }

    ~gsEigenProblemBase() {};

public:

    virtual gsOptionList defaultOptions()
    {
        gsOptionList options;
        options.addSwitch("verbose","Verbose output",false);
        options.addInt("solver","Spectra solver to be used (see Spectra documentation):"
                                    "0: Cholesky (default)"
                                    "1: RegularInverse"
                                    "2: ShiftInvert"
                                    "3: Buckling"
                                    "4: Cayley",0);

        options.addInt("selectionRule","Selection rule to be used (see Spectra documentation):"
                                        "0: LargestMagn"
                                        "1: LargestReal"
                                        "2: LargestImag"
                                        "3: LargestAlge"
                                        "4: SmallestMagn"
                                        "5: SmallestReal"
                                        "6: SmallestImag"
                                        "7: SmallestAlge"
                                        "8: BothEnds",4);

        options.addInt("sortRule","Sort rule to be used (see Spectra documentation):"
                                        "0: LargestMagn"
                                        "1: LargestReal"
                                        "2: LargestImag"
                                        "3: LargestAlge"
                                        "4: SmallestMagn"
                                        "5: SmallestReal"
                                        "6: SmallestImag"
                                        "7: SmallestAlge"
                                        "8: BothEnds",4);

        options.addInt("ncvFac","Factor for Spectra's ncv number. Ncv = ncvFac * numEigenvalues",3);
        return options;
    }

    gsOptionList & options() {return m_options; };

    virtual gsStatus compute();
    virtual gsStatus compute(const T shift);

    virtual gsStatus computeSparse(const T shift = 0.0, const index_t number = 10);

    virtual gsStatus computePower();

    virtual gsMatrix<T> values() const { return m_values; };
    virtual T value(int k) const { return m_values.at(k); };

    virtual gsMatrix<T> vectors() const { return m_vectors; };
    virtual gsMatrix<T> vector(int k) const { return m_vectors.col(k); };

    virtual std::vector<std::pair<T,gsMatrix<T>> > mode(int k) const {return makeMode(k); }

protected:

    virtual std::vector<std::pair<T,gsMatrix<T>> > makeMode(int k) const;

private:
    #ifdef GISMO_WITH_SPECTRA
    template<Spectra::GEigsMode _GEigsMode>
    typename std::enable_if<_GEigsMode==Spectra::GEigsMode::Cholesky ||
                            _GEigsMode==Spectra::GEigsMode::RegularInverse
                            ,
                            gsStatus>::type computeSparse_impl(T shift, index_t number);
    #endif

    #ifdef GISMO_WITH_SPECTRA
    template<Spectra::GEigsMode _GEigsMode>
    typename std::enable_if<_GEigsMode==Spectra::GEigsMode::ShiftInvert ||
                            _GEigsMode==Spectra::GEigsMode::Buckling ||
                            _GEigsMode==Spectra::GEigsMode::Cayley
                            ,
                            gsStatus>::type computeSparse_impl(T shift, index_t number);
    #endif

protected:

    gsSparseMatrix<T> m_A;
    gsSparseMatrix<T> m_B;

    gsOptionList m_options;

    gsEigen::GeneralizedSelfAdjointEigenSolver< typename gsMatrix<T>::Base >  m_eigSolver;

    gsMatrix<T> m_values,m_vectors;

    index_t m_num;

    gsStatus m_status;
};


} // namespace gismo


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsEigenProblemBase.hpp)
#endif
