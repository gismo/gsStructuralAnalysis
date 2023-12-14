/** @file gsALMHelper.hpp

    @brief Provides helper functions for the ALM

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#pragma once

#include <typeinfo>

using namespace gismo;

// Miscelaneous functions
/// sign function
template <class T>
index_t sign(T val)
{
    return (T(0) < val) - (val < T(0));
}

/// sort vector
template <class T>
index_t countNegatives(gsVector<T> vec)
{
  index_t count = 0;
  index_t N = vec.cols();
  index_t M = vec.rows();
  for(index_t i = 0; i < M; i++)
        for(index_t j = 0; j < N; j++)
        {
            if( vec(i,j) < 0 )
                count += 1;
        }
    return count;
}