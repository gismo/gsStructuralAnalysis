 /** @file gsContinuationBase.h

    @brief Base class for continuation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#pragma once

namespace gismo
{

/**
    @brief Base class for simple continuation schemes

    \tparam T coefficient type

    \ingroup gsStructuralAnalysis
*/
template <class T>
class gsContinuationBase
{
public:
    virtual ~gsContinuationBase() {};
    gsContinuationBase()
    {
        defaultOptions();
    }


    virtual gsStatus step(T dL) = 0;
    virtual gsVector<T> & solutionU() = 0;
    virtual T solutionL() = 0;

public:
    void defaultOptions()
    {
        m_options.addInt("numSteps","number of steps to be taken",10);
        m_options.addReal("dL","Step size",1.);
    }

    gsOptionList options() { return m_options; }

protected:
    gsOptionList m_options;

};

} // namespace gismo
