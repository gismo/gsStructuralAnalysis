 /** @file gsAdaptiveSpaceTime.h

    @brief Performs the arc length method to solve a nonlinear equation system.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#pragma once

#include <gsHSplines/gsKdNode.h>

namespace gismo
{

template<class T, class solution_t >
class gsAdaptiveSpaceTime
{

public:

  ~gsAdaptiveSpaceTime()
  {

  }

  /**
   * @brief      Initializes the hierarchy given the first level of solutions
   *
   * @param[in]  times      The times (MONOTONICALLY INCREASING!)
   * @param[in]  solutions  The solutions
   */
  gsAdaptiveSpaceTime(const std::vector<T> & times, const  std::vector<solution_t> & solutions)
  :
  m_initialized(false)
  {
    GISMO_ASSERT(times.size()==solutions.size(),"Sizes must agree");

    // Initialize the knot vectors
    m_t = gsKnotVector<T>(times,1,0);
    m_t.degreeDecrease(1);

    m_xi = m_t;
    m_xi.transform(0,1);

    for (size_t k=0; k!=m_xi.size(); ++k)
      m_solutions.insert({m_xi.at(k),solutions.at(k)});

    _defaultOptions();
  }

  /**
   * @brief      Initializes the hierarchy without solutions
   *
   * @param[in]  times      The times (MONOTONICALLY INCREASING!)
   */
  gsAdaptiveSpaceTime(const std::vector<T> & times)
  :
  m_initialized(false)
  {
    // Initialize the knot vectors
    m_t = gsKnotVector<T>(times,1,0);
    m_t.degreeDecrease(1);

    m_xi = m_t;
    m_xi.transform(0,1);

    _defaultOptions();
  }

  // gsAdaptiveSpaceTime(const T tmin, const T tmax, const index_t N = 10)
  // :
  // m_initialized(false)
  // {
  //   std::vector<T>
  // }


  gsAdaptiveSpaceTime() {};

private:

  void _defaultOptions()
  {
    m_maxLevel = 5;
    m_tolerance = 1e-1;
    m_split = false;
    m_options.addInt("MaxLevel","Sets the maximum level for hierarchical refinement",m_maxLevel);
    m_options.addInt("Split","Adds the jobs as sub-intervals, so that references with intermediate points can be computed. If this option is false (default) you have to perform TWO steps per item"
                      ,m_split);
    m_options.addReal("Tolerance","Relative tolerance",m_tolerance);
  }

  void _applyOptions()
  {
    m_maxLevel = m_options.getInt("MaxLevel");
    m_split = m_options.getInt("Split");
    m_tolerance = m_options.getReal("Tolerance");
  }

public:

  void init()
  {
    this->_applyOptions();
    // Clean queue
    std::queue<std::pair<T,T>> queue;
    m_queue.swap(queue);

    T low,upp;
    for (typename gsKnotVector<T>::const_iterator it = m_xi.begin(); it != std::prev(m_xi.end()); )
    {
      low = *it;
      it++;
      upp = *it;

      gsDebugVar(*it);

      m_queue.push({low,upp});
    }

    this->_buildMap();

    m_initialized = true;
  }

  gsOptionList & options() { return m_options; }

  //         ID,      tstart, tend, dt, start     , guess
  std::tuple<index_t, T     , T   , T , solution_t, solution_t> pop()
  {
    GISMO_ASSERT(m_initialized,"Structure is not initialized");
    GISMO_ASSERT(!m_queue.empty(),"The queue is empty! Something went wrong.");

    T xilow = m_queue.front().first;
    T xiupp = m_queue.front().second;
    m_queue.pop();
    T tlow = m_tmap[xilow];
    T tupp = m_tmap[xiupp];

    T dt  = (tupp-tlow);

    GISMO_ASSERT(m_solutions.count(xilow)!=0,"Cannot find start point at tstart = "<<tlow<<"\n");
    GISMO_ASSERT(m_solutions.count(xiupp)!=0,"Cannot find end point at tend = "<<tupp<<"\n");

    solution_t start = m_solutions[xilow];
    solution_t guess = m_solutions[xiupp];

    // add the job to the active jobs
    m_jobs[m_ID++] = std::make_pair(xilow,xiupp);
    gsDebug<<"Added active job (ID="<<m_ID-1<<") on interval = ["<<xilow<<","<<xiupp<<"] = ["<<tlow<<","<<tupp<<"]\n";

    return std::make_tuple(m_ID-1,tlow, tupp, dt, start, guess);
  }

  bool getReferenceByTime(T time, solution_t & result)
  {
    if (m_solutions.count(m_ximap[time])!=0)
    {
      result =  m_solutions[m_ximap[time]];
      return true;
    }
    else
      return false;
  }

  bool getReferenceByPar(T xi, solution_t & result)
  {
    if (m_solutions.count(xi)!=0)
    {
      result =  m_solutions[xi];
      return true;
    }
    else
      return false;
  }

  bool getReferenceByID(index_t ID, solution_t & result)
  {
    if (m_solutions.count(m_jobs[ID].second)!=0)
    {
      result =  m_solutions[m_jobs[ID].second];
      return true;
    }
    else
      return false;
  }


  // /**
  //  * @brief      Submits job \a ID to the solution database.
  //  *
  //  * The solution for job\a ID is added to the database
  //  *
  //  * @param[in]  ID        ID of the job
  //  * @param[in]  solution  The solution for the specific ID
  //  */
  // void submit(index_t ID, std::vector<solution_t> solutions)
  // {
  //   GISMO_ASSERT(m_initialized,"Structure is not initialized");
  //   GISMO_ASSERT(solutions.size()==m_steps,"The size of solutions ("<<solutions.size()<<") does not match the number of intermediate steps ("<<m_steps<<")");
  //   T low, upp;

  //   std::tie(low,upp) = m_jobs[ID];
  //   T ds = upp-low;
  //   for (index_t k=0; k!=m_steps; k++)
  //     m_solutions[low + (k+1)*ds] = solutions[k];

  //   gsDebug<<"Jobs submitted at t = ["<<low<","<<upp<<"]"<<"\n";
  // }

  /**
   * @brief      { function_description }
   *
   * @param[in]  ID         { parameter_description }
   * @param      distances  The distances
   * @param      ptimes     The parametric times on which the solutions are computed ()
   * @param      solutions  The solutions
   */
  void submit(index_t ID, const std::vector<T> & distances, const std::vector<T> & times, const std::vector<solution_t> & solutions)
  {
    GISMO_ASSERT(m_initialized,"Structure is not initialized");
    GISMO_ASSERT(distances.size()==solutions.size()+1,"You must provide one more distance than solutions");
    GISMO_ASSERT(times.size()==solutions.size(),"number of parametric times must be equal to the number of solutions");

    T tlow, tupp, xilow, xiupp, dxi, Dt, dt;
    std::vector<T> t,xi;

    // Compute interval (xi and t)
    std::tie(xilow,xiupp) = m_jobs[ID];
    tlow = m_tmap[xilow];
    tupp = m_tmap[xiupp];

    // Compute interval distance
    dxi = xiupp - xilow;
    Dt = tupp - tlow; // original time step

    // Get time points
    t.resize(distances.size()+1);
    t.at(0) = tlow;
    for (size_t k=0; k!=t.size()-1; k++)
      t.at(k+1) = t.at(k) + distances.at(k);

    T dt_tmp = 0;
    for(typename std::vector<T>::const_iterator it = distances.begin(); it != std::prev(distances.end()); ++it)
        dt_tmp += *it;
    GISMO_ASSERT((Dt-dt_tmp)/Dt<1e-12,"Total distance of the computed intervals should be equal to the original interval length ("<<dt_tmp<<"!="<<Dt<<")");

    // Get actual time step and the error
    dt = t.back()-tlow;

    // Fill and scale xi vector
    xi.resize(distances.size()+1);
    xi.front() = xilow;
    xi.back() = xiupp;
    for (size_t k=1; k!=xi.size()-1; k++)
        xi.at(k) = xilow + dxi * ( t.at(k) - t.front() ) / dt;

    for (index_t k = 0; k!= xi.size(); k++)
      gsDebug<<"xi = "<<xi.at(k)<<"; t = "<<t.at(k)<<"\n";

    // Transform the time vector
    m_t.addConstant(tupp,dt-Dt);

    for (size_t k=1; k!=t.size()-1; k++)
    {
        m_t.insert(t.at(k));
        gsDebug<<"Will submit "<<t.at(k)<<" in t vector\n";
        m_xi.insert(xi.at(k));
    }

    this->_buildMap();

    for (size_t k=1; k!=xi.size()-1; k++) // add only INTERIOR solutions
    {
      gsInfo<<"Added a solution on time = "<<m_tmap[xi.at(k)]<<" (parametric time = "<<xi.at(k)<<")\n";
      m_solutions[xi.at(k)] = solutions.at(k-1);
    }

    m_tmp[ID] = xi;
  }

  void addJobs(index_t ID)
  {
    GISMO_ASSERT(m_tmp[ID].size()>2,"Job for ID "<<ID<<" has not been submitted yet, or no data within the interval is known");

    T tlow, tupp, xilow, xiupp, dt, tcurr, tprev, dt_last;

    xilow = m_tmp[ID].front();
    xiupp = m_tmp[ID].back();
    tlow = m_tmap[xilow];
    tupp = m_tmap[xiupp];

    index_t N = m_tmp[ID].size()-1;
    dt = tupp - tlow;
    dt_last = tupp - m_tmap[m_tmp[ID].at(N-1)];

    std::vector<T> distances(N);
    tprev = m_tmap[m_tmp[ID].at(0)];
    for (index_t k = 1; k!=N+1; k++)
    {
      tcurr = m_tmap[m_tmp[ID].at(k)];
      distances.at(k-1) = (tcurr - tprev) / dt;
      tprev = tcurr;
    }

    gsDebugVar(gsAsVector(distances));

    if (!(distances.back() < m_tolerance))
    {
      for (index_t k = 1; k!=m_tmp[ID].size(); k++)
      {
        m_queue.push(std::make_pair(m_tmp[ID].at(k-1),m_tmp[ID].at(k)));
        gsInfo<<"Interval ["<<m_tmp[ID].at(k-1)<<","<<m_tmp[ID].at(k)<<"] added to queue\n";
      }
    }
  }

  void finishJob(index_t ID)
  {
    // GISMO_ASSERT(m_tmp[ID].size()==m_steps,"Too early to finish; I only have "<<m_tmp[ID].size()<<" points for ID "<<ID<<" but I need "<<m_steps<<" points");
    m_tmp.erase(ID);
  }

  void printActiveJobs()
  {
    // GISMO_ASSERT(m_tmp[ID].size()==m_steps,"Too early to finish; I only have "<<m_tmp[ID].size()<<" points for ID "<<ID<<" but I need "<<m_steps<<" points");
    for (typename std::map<index_t,std::vector<T>>::const_iterator it = m_tmp.cbegin(); it!= m_tmp.cend(); it++)
    {
      std::string tmp = it->second.size()==0 ? "unfinished " : "finished ";
      gsInfo<<"ID = "<<it->first<<": "<<tmp<<"\n";
    }
  }



  // void submit(index_t ID, std::vector<T> & reltime, std::vector<solution_t> & solutions)
  // {
  //   GISMO_ASSERT(m_initialized,"Structure is not initialized");
  //   GISMO_ASSERT(solutions.size()==m_steps,"The size of solutions ("<<solutions.size()<<") does not match the number of intermediate steps ("<<m_steps<<")");
  //   GISMO_ASSERT(solutions.size()==reltime,"The size of solutions ("<<solutions.size()<<") does not match the number of relative times ("<<reltime.size()<<")");

  //   T low, upp;

  //   std::tie(low,upp) = m_jobs[ID];

  //   gsDebug<<"Job submitted at t = "<<time<<"\n";
  // }


  /**
   * @brief       Adds a job based on the job \a ID in the tree.
   *              Example: If current job is {0.2,4} on level l,
   *              the job that is added is XXXX
   *
   * @param[in]  ID    { parameter_description }
   */
  // void addJob(index_t ID)
  // {
  //   GISMO_ASSERT(m_initialized,"Structure is not initialized");
  //   index_t level;
  //   T time;
  //   std::tie(level,time) = m_jobs[ID];

  //   if (level+1 > m_maxLevel)
  //   {
  //     gsWarn<<"Max level reached!";
  //     return;
  //   }

  //   gsInfo<<"ID = "<<ID<<" level = "<<level<<" time = "<<time<<"\n";

  //   // Collects the time interval in the tree on \a level and \a time
  //   gsVector<T> low, upp;
  //   std::tie(low,upp) = _interval(m_tree[level][time]);

  //   // Splits the current time interval [low,upp].
  //   // m_tree[level][time]->split();
  //   m_tree[level][time]->split(0,(upp[0]+low[0])/2);

  //   // Add the new intervals to the tree
  //   m_tree[level+1][m_tree[level][time]->left ->lowCorner()[0]] = m_tree[level][time]->left;
  //   m_tree[level+1][m_tree[level][time]->right->lowCorner()[0]] = m_tree[level][time]->right;

  //   // Adds two new jobs on level \a l+1, for both new intervals. The jobs start at tje lower corner (i.e. the start of the intervals)
  //   m_queue.push({level+1,m_tree[level][time]->left ->lowCorner()[0]});
  //   m_queue.push({level+1,m_tree[level][time]->right->lowCorner()[0]});
  // }

  // void addJob(T low, T upp)
  void addJob(index_t ID)
  {
    // Find all sub-intervals in the interval of job ID
    T low, upp;
    std::tie(low,upp) = m_jobs[ID];

    // for (index_t)

    // m_queue.push({low,upp});
  }

  void removeJob(index_t ID)
  {
    // Remove ID
    gsDebug<<"Erasing finished job on level = "<<m_jobs[ID].first<<"; time = "<<m_jobs[ID].second<<"\n";
    m_jobs.erase(ID);
    gsDebug<<"Current active jobs\n";
    for (typename std::map<index_t,std::pair<index_t,T>>::const_iterator it = m_jobs.begin(); it!=m_jobs.end(); ++it)
      gsDebug<<"ID = "<<it->first<<"; level = "<<it->second.first<<"; time = "<<it->second.second<<"\n";
  }

  // index_t currentLevel(index_t ID)
  // {
  //   return _activeJob(ID).first;
  // }

  bool empty()
  {
    return m_queue.empty();
  }

  std::pair<std::vector<T>,std::vector<solution_t>> getFlatSolution()
  {
    std::vector<T> times;
    std::vector<solution_t> solutions;

    for (typename std::map<T,solution_t>::const_iterator it = m_solutions.begin(); it!=m_solutions.end(); ++it)
    {
      times.push_back(m_tmap[it->first]);
      solutions.push_back(it->second);
    }

    return std::make_pair(times,solutions);
  }

  void print()
  {
    for (typename std::map<T,solution_t>::const_iterator it = m_solutions.begin(); it!=m_solutions.end(); ++it)
    {
      gsInfo<<"\t";
      gsInfo<<"time = "<<it->first<<"\tsol = "<<it->second<<"\n";
    }
  }

  void printQueue()
  {
    std::queue<std::pair<T,T>> queue = m_queue;
    gsInfo<<"Queue has "<<queue.size()<<" elements\n";
    while (!queue.empty())
    {
      gsInfo<<"Start time "<<queue.front().first<<"; end time = "<<queue.front().second<<"\n";
      queue.pop();
    }
    gsInfo<<"Queue is empty\n";

  }

  void printKnots()
  {
    gsInfo<<"xi = \n"<<m_xi.asMatrix()<<"\n";
    gsInfo<<"t = \n"<<m_t.asMatrix()<<"\n";
  }

  // void printActive()
  // {
  //   gsVector<T,1> low, upp;
  //   gsInfo<<"The following "<<m_jobs.size()<<" jobs are active\n";
  //   for (index_t k=0; k!=m_jobs.size(); k++)
  //   {
  //     gsInfo<<"Job "<<k<<" has level"
  //     gsInfo<<"Level "<<queue.front().first<<"; time "<<queue.front().second<<"\n";
  //     queue.pop();
  //   }
  //   gsInfo<<"Queue is empty\n";

  // }

protected:
  std::pair<index_t,T> _activeJob(index_t ID)
  {
    if (m_jobs.count(ID) == 1)
      return m_jobs[ID];
    else if (m_jobs.count(ID) > 1)
    {
      gsWarn<<"Multiple jobs with ID "<<ID<<" exist! Did something go wrong?\n";
      return m_jobs[ID];
    }
    else
      GISMO_ERROR("No job found with ID "<<ID<<". Did something go wrong?");
  }

  void _buildMap()
  {
    GISMO_ASSERT(m_xi.size()==m_t.size(),"Time and parametric time don't have equal length");
    m_tmap.clear();
    m_ximap.clear();
    for (size_t k=0; k!=m_xi.size(); k++)
    {
        m_tmap.insert({m_xi[k],m_t[k]});
        m_ximap.insert({m_t[k],m_xi[k]});
    }
  }


protected:
  index_t m_points;
  index_t m_maxLevel;
  bool m_split;
  bool m_initialized;

  T m_tolerance;

  // solution map, stores the solution per parametric value
  std::map<T,solution_t>       m_solutions;

  // map that maps GIVEN a parametric value TO a time
  std::map<T,T> m_tmap, m_ximap;

  // Parametric domain
  gsKnotVector<T> m_xi;
  // Temporal domain
  gsKnotVector<T> m_t;

  // Stores [xi_i,xi_i+1]
  std::queue<std::pair<T,T>>          m_queue;

  // Stores ID, [xi_i,xi_i+1]
  std::map<index_t,std::pair<T,T>>    m_jobs;

  mutable index_t                           m_ID = (0);

  gsOptionList m_options;

  std::map<index_t,std::vector<T>> m_tmp;

};

}


// #ifndef GISMO_BUILD_LIB
// #include GISMO_HPP_HEADER(gsAdaptiveSpaceTime.hpp)
// #endif
