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
#include <deque>

namespace gismo
{

template<class T, class solution_t >
class gsAPALMData
{

public:

  ~gsAPALMData()
  {

  }

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

  void addStartPoint(const T & time, const solution_t & solution, bool priority = false);

  void appendPoint(bool priority = false);

  void appendData(const T & time, const solution_t & solution, bool priority=false);

  void setData(const std::vector<T> & times, const std::vector<solution_t> & solutions);

  void init(const std::vector<T> & times, const std::vector<solution_t> & solutions);

  void init();

  gsOptionList & options() { return m_options; }

  //         ID,      dt, start     , prev
  std::tuple<index_t, T , solution_t, solution_t> pop();

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

  T jobStartTime(index_t ID);
  T jobStartPar(index_t ID);

  std::pair<T,T> jobTimes(index_t ID);
  std::pair<T,T> jobPars(index_t ID);
  index_t jobLevel(index_t ID);

  void printActiveJobs();

  bool empty();

  std::tuple<std::vector<T>,std::vector<solution_t>,std::vector<index_t>> getFlatSolution(index_t level=-1);

  void print();

  void printQueue();

  void printKnots();

  size_t nActive() { return m_jobs.size(); }
  size_t nWaiting() { return m_queue.size(); }

  size_t maxLevel() { return m_maxLevel; }

  void setLength(T dt) { m_dt = dt; }
  T getLength() { return m_dt; }

protected:
  std::tuple<T,T,index_t> _activeJob(index_t ID);

  void _buildMap();


protected:
  index_t m_points;
  index_t m_maxLevel;
  index_t m_verbose;
  bool m_initialized;

  // Default arc-length
  T m_dt;

  T m_tolerance;

  // solution map, stores the solution per parametric value
  std::map<T,std::shared_ptr<solution_t>> m_solutions;
  std::map<T,std::shared_ptr<solution_t>> m_prevs;
  // map that maps GIVEN a parametric value TO a level
  std::map<T,index_t>          m_levels;

  // map that maps GIVEN a parametric value TO a time
  std::map<T,T> m_tmap;
  // map that maps GIVEN a time TO a parametric value
  std::map<T,T> m_ximap;

  // Parametric domain
  gsKnotVector<T> m_xi;
  // Temporal domain
  gsKnotVector<T> m_t;

  // Stores [xi_i,xi_i+1,level]
  std::deque<std::tuple<T,T,index_t>>          m_queue;

  // Stores ID, [xi_i,xi_i+1,level]
  // map that maps GIVEN an ID TO a parametric interval
  std::map<index_t,std::tuple<T,T,index_t>>    m_jobs;

  mutable index_t                           m_ID = (0);

  gsOptionList m_options;

};

}

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsAPALMData.hpp)
#endif
