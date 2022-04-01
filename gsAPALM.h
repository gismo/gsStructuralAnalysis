 /** @file gsAPALM.h

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
#include <gsStructuralAnalysis/gsALMBase.h>

namespace gismo
{

template<class T, class solution_t>
class gsAPALM
{

public:

  virtual ~gsAPALM() { }

  /**
   * @brief      Constructs a new instance.
   *
   * @param[in]  ALM   The alms
   */
  gsAPALM(      gsALMBase<T> * ALM,
              const gsAPALMData<T,solution_t> & Data);

  gsAPALM() {};

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

  const gsAPALMData<T,solution_t> & getHierarchy() const { return m_data; }
  gsAPALMData<T,solution_t> getHierarchy()  { return m_data; }

  const std::vector<solution_t>   & getFlatSolutions() const { return m_solutions; }
  const std::vector<T>            & getFlatTimes()     const { return m_times; }
  const std::vector<index_t>      & getFlatLevels()    const { return m_levels; }

  std::vector<solution_t>   getFlatSolutions()  { return m_solutions; }
  std::vector<T>            getFlatTimes()      { return m_times; }
  std::vector<index_t>      getFlatLevels()     { return m_levels; }

  const std::vector<solution_t *>   & getSolutions(index_t level) const { return m_lvlSolutions[level]; }
  const std::vector<T *>            & getTimes(index_t level)     const { return m_lvlTimes[level]; }

  std::vector<solution_t>   getSolutions(index_t level)
  {
    std::vector<solution_t> result;
    for (typename std::vector<solution_t *>::iterator it=m_lvlSolutions[level].begin(); it!=m_lvlSolutions[level].end(); it++)
      result.push_back(**it);
    return result;
  }

  std::vector<T>   getTimes(index_t level)
  {
    std::vector<T> result;
    for (typename std::vector<T *>::iterator it=m_lvlTimes[level].begin(); it!=m_lvlTimes[level].end(); it++)
      result.push_back(**it);
    return result;
  }

  std::vector<std::vector<solution_t>> getSolutionsPerLevel()
  {
    std::vector<std::vector<solution_t>> result(m_lvlSolutions.size());
    for (index_t l=0; l!=m_lvlSolutions.size(); l++)
      for (typename std::vector<solution_t *>::iterator it=m_lvlSolutions[l].begin(); it!=m_lvlSolutions[l].end(); it++)
        result[l].push_back(**it);
    return result;
  }

  std::vector<std::vector<T>> getTimesPerLevel()
  {
    std::vector<std::vector<T>> result(m_lvlTimes.size());
    for (index_t l=0; l!=m_lvlSolutions.size(); l++)
      for (typename std::vector<T *>::iterator it=m_lvlTimes[l].begin(); it!=m_lvlTimes[l].end(); it++)
        result[l].push_back(**it);
    return result;
  }

protected:

  std::vector<solution_t> m_solutions;
  std::vector<T> m_times;
  std::vector<index_t> m_levels;

  std::vector<std::vector<solution_t * > >  m_lvlSolutions;
  std::vector<std::vector<T * > >           m_lvlTimes;

  gsALMBase<T> * m_ALM;
  gsAPALMData<T,solution_t> m_data;

  bool m_verbose;
  index_t m_subIntervals;

  gsOptionList m_options;
};

}

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsAPALM.hpp)
#endif
