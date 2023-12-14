 /** @file gsAPALMBase.h

    @brief Performs the arc length method to solve a nonlinear equation system.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsIO/gsOptionList.h>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMBase.h>

namespace gismo
{

template<class T, class solution_t>
class gsAPALMBase
{

public:

  virtual ~gsAPALMBase() { }

  /**
   * @brief      Constructs a new instance.
   *
   * @param[in]  ALM   The alms
   */
  gsAPALMBase(      gsALMBase<T> * ALM,
              const gsAPALMData<T,solution_t> & Data);

  gsAPALMBase() {};

private:

  virtual void _defaultOptions();

  virtual void _getOptions();

public:
  virtual void initialize();

  virtual void serialSolve(index_t Nsteps = 10);
  virtual void parallelSolve();

  virtual void serialStepOutput(const gsVector<T> & U, const T & L) {};
  virtual void parallelStepOutput(const gsVector<T> & U, const T & L) {};

  gsOptionList & options() { return m_options; }

  std::vector<solution_t>   getSolutions()  { return m_solutions; }
  std::vector<T>            getTimes()      { return m_times; }
  gsAPALMData<T,solution_t> getHierarchy()  { return m_data; }

protected:

  std::vector<solution_t> m_solutions;
  std::vector<T> m_times;

  gsALMBase<T> * m_ALM;
  gsAPALMData<T,solution_t> m_data;

  bool m_verbose;

  gsOptionList m_options;
};

}

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsAPALMBase.hpp)
#endif
