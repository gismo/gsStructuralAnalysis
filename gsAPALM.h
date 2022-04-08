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
#include <gsStructuralAnalysis/gsAPALMDataContainer.h>

namespace gismo
{

template<class T>
class gsAPALM
{

public:
  typedef typename std::pair<gsVector<T>,T> solution_t;

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

  void _initStart(const gsVector<T> & Ustart, const T & Lstart, const T & dL);
  void _initStart(const std::vector<gsVector<T>> & Ustart, const std::vector<T> & Lstart, const std::vector<T> & dLs);

public:
  virtual void initialize();

  virtual void serialSolve(index_t Nsteps = 10);
  virtual void parallelSolve();

  virtual void serialStepOutput(const std::pair<gsVector<T>,T> & pair, const T & time, index_t step) {};
  virtual void parallelStepOutput(const std::pair<gsVector<T>,T> & pair, const T & time, index_t step) {};
  virtual void parallelIntervalOutput(const std::vector<std::pair<gsVector<T>,T>> & stepSolutions, const std::vector<T> & stepTimes, index_t level, index_t ID) {};

  gsOptionList & options() { return m_options; }

  const gsAPALMDataContainer<T,solution_t> & getHierarchy() const { return m_data; }
  gsAPALMDataContainer<T,solution_t> getHierarchy()  { return m_data; }

  const std::vector<solution_t>   & getFlatSolutions(index_t branch = 0) const { return m_solutions[branch]; }
  const std::vector<T>            & getFlatTimes(index_t branch = 0)     const { return m_times[branch]; }
  const std::vector<index_t>      & getFlatLevels(index_t branch = 0)    const { return m_levels[branch]; }

  std::vector<solution_t>   getFlatSolutions(index_t branch = 0)  { return m_solutions[branch]; }
  std::vector<T>            getFlatTimes(index_t branch = 0)      { return m_times[branch]; }
  std::vector<index_t>      getFlatLevels(index_t branch = 0)     { return m_levels[branch]; }

  const std::vector<solution_t *>   & getSolutions(index_t level, index_t branch = 0) const { return m_lvlSolutions[branch][level]; }
  const std::vector<T *>            & getTimes(index_t level, index_t branch = 0)     const { return m_lvlTimes[branch][level]; }

  std::vector<solution_t>   getSolutions(index_t level, index_t branch = 0)
  {
    std::vector<solution_t> result;
    for (typename std::vector<solution_t *>::iterator it=m_lvlSolutions[branch][level].begin(); it!=m_lvlSolutions[branch][level].end(); it++)
      result.push_back(**it);
    return result;
  }

  std::vector<T>   getTimes(index_t level, index_t branch = 0)
  {
    std::vector<T> result;
    for (typename std::vector<T *>::iterator it=m_lvlTimes[branch][level].begin(); it!=m_lvlTimes[branch][level].end(); it++)
      result.push_back(**it);
    return result;
  }

  std::vector<std::vector<solution_t>> getSolutionsPerLevel(index_t branch = 0)
  {
    std::vector<std::vector<solution_t>> result(m_lvlSolutions[branch].size());
    for (index_t l=0; l!=m_lvlSolutions[branch].size(); l++)
      for (typename std::vector<solution_t *>::iterator it=m_lvlSolutions[branch][l].begin(); it!=m_lvlSolutions[branch][l].end(); it++)
        result[l].push_back(**it);
    return result;
  }

  std::vector<std::vector<T>> getTimesPerLevel(index_t branch = 0)
  {
    std::vector<std::vector<T>> result(m_lvlTimes[branch].size());
    for (index_t l=0; l!=m_lvlSolutions[branch].size(); l++)
      for (typename std::vector<T *>::iterator it=m_lvlTimes[branch][l].begin(); it!=m_lvlTimes[branch][l].end(); it++)
        result[l].push_back(**it);
    return result;
  }

protected:

  std::vector<std::vector<solution_t>>  m_solutions;
  std::vector<std::vector<T>>           m_times;
  std::vector<std::vector<index_t>>     m_levels;

  std::queue<std::tuple<solution_t,T,bool>> m_starts; // solution, step length, true if start was a bifurcation

  std::vector<std::vector<std::vector<solution_t * > > >  m_lvlSolutions;
  std::vector<std::vector<std::vector<T * > > >           m_lvlTimes;

  gsALMBase<T> * m_ALM;
  gsAPALMData<T,solution_t> m_dataEmpty;
  gsAPALMDataContainer<T,solution_t> m_data;

  bool m_verbose;
  index_t m_subIntervals;

  gsOptionList m_options;

  bool m_singularPoint;
  T m_branchLengthMult;
  T m_bifLengthMult;
};

}

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsAPALM.hpp)
#endif
