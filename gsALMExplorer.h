 /** @file gsArcLengthIterator.h

    @brief Performs the arc length method to solve a nonlinear equation system.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#pragma once
#include <gsStructuralAnalysis/gsALMBase.h>

namespace gismo
{

/**
    @brief Performs the arc length method to solve a nonlinear equation system.

    \tparam T coefficient type

    \ingroup gsStructuralAnalysis
*/
template <class T>
class gsALMExplorer
{
public:

    virtual ~gsALMExplorer() {};

    /// Constructor giving access to the gsShellAssembler object to create a linear system per iteration
    gsALMExplorer(  const gsALMBase<T> * ALM )
    : m_BaseALM(ALM)
    {
        // m_DummyALM = m_BaseALM->clone();
        // gsDebugVar(m_DummyALM->options());
    }

// General functions
public:

    void defaultOptions();
    void setOptions(gsOptionList opt) { m_options.update(opt,gsOptionList::addIfUnknown); } // gsOptionList opt
    void getOptions() const;

protected:
    const gsALMBase<T> * m_BaseALM;
    typename gsALMBase<T>::uPtr m_DummyALM;
    gsOptionList m_ALMOptions;
    gsOptionList m_options;
};


} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsALMExplorer.hpp)
#endif