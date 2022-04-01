 /** @file gsAPALMData.h

    @brief Performs the arc length method to solve a nonlinear equation system.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsNurbs/gsKnotVector.h>
#include <gsIO/gsOptionList.h>
#include <gsHSplines/gsKdNode.h>
#include <queue>

namespace gismo
{

template<class T, class solution_t >
class gsAPALMData
{

public:

  ~gsAPALMData() { }

  /**
   * @brief      Initializes the hierarchy given the first level of solutions
   *
   * @param[in]  times      The times (MONOTONICALLY INCREASING!)
   * @param[in]  solutions  The solutions
   */
  gsAPALMData(const std::vector<T> & times, const  std::vector<solution_t> & solutions);

  /**
   * @brief      Initializes the hierarchy without solutions
   *
   * @param[in]  times      The times (MONOTONICALLY INCREASING!)
   */
  gsAPALMData(const std::vector<T> & times);

  gsAPALMData();

private:

  void _defaultOptions();

  void _applyOptions();

public:

  void setData(const std::vector<T> & times, const  std::vector<solution_t> & solutions);

  void init(const std::vector<T> & times, const  std::vector<solution_t> & solutions);

  void init();

  gsOptionList & options() { return m_options; }

  //         ID,      tstart, tend, dt, start     , guess
  std::tuple<index_t, T     , T   , T , solution_t, solution_t> pop();

  bool getReferenceByTime(T time, solution_t & result);

  bool getReferenceByPar(T xi, solution_t & result);

  bool getReferenceByID(index_t ID, solution_t & result);

  /**
   * @brief      { function_description }
   *
   * @param[in]  ID         { parameter_description }
   * @param      distances  The distances
   * @param      solutions  The solutions
   * @param      upperError The error at the end of the domain
   * @param      lowerError The error in the first (known) part of the domain. By default this is 0, meaning that the lower part will never be refined
   */
  void submit(index_t ID, const std::vector<T> & distances,  std::vector<solution_t> solutions, const T & upperError, const T & lowerError = 0);

  void finishJob(index_t ID);

  void printActiveJobs();

  void removeJob(index_t ID);

  bool empty();

  std::pair<std::vector<T>,std::vector<solution_t>> getFlatSolution();

  void print();

  void printQueue();

  void printKnots();

  size_t nActive() { return m_tmp.size(); }
  size_t nWaiting() { return m_queue.size(); }

protected:
  std::pair<index_t,T> _activeJob(index_t ID);

  void _buildMap();


protected:
  index_t m_points;
  index_t m_maxLevel;
  index_t m_verbose;
  bool m_initialized;

  T m_tolerance;

  // solution map, stores the solution per parametric value
  std::map<T,solution_t>       m_solutions;
  std::map<T,solution_t * >    m_guesses;

  // map that maps GIVEN a parametric value TO a time
  std::map<T,T> m_tmap;
  // map that maps GIVEN a time TO a parametric value
  std::map<T,T> m_ximap;
  // map that maps GIVEN a parametric value TO a level
  std::map<T,index_t> m_lmap;

  // Parametric domain
  gsKnotVector<T> m_xi;
  // Temporal domain
  gsKnotVector<T> m_t;

  // Stores [xi_i,xi_i+1]
  std::queue<std::pair<T,T>>          m_queue;

  // Stores ID, [xi_i,xi_i+1]
  // map that maps GIVEN an ID TO a parametric interval
  std::map<index_t,std::pair<T,T>>    m_jobs;

  mutable index_t                           m_ID = (0);

  gsOptionList m_options;

  std::map<index_t,std::vector<T>> m_tmp;

};

}

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsAPALMData.hpp)
#endif
