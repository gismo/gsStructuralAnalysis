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
    @brief Performs the arc length method to solve a nonlinear equation system.

    \tparam T coefficient type

    \ingroup gsStructuralAnalysis
*/
template <class T>
class gsContinuationBase
{
public:
    virtual ~gsContinuationBase() {};

    virtual void step(T dL) = 0;
    virtual gsVector<T> solutionU() = 0;
    virtual T solutionL() = 0;

};

} // namespace gismo
