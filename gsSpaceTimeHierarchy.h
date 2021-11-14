 /** @file gsSpaceTimeHierarchy.h

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
class gsSpaceTimeHierarchy
{

  typedef kdnode<1, T> node_t;

  typedef typename node_t::point point_t;

public:

  ~gsSpaceTimeHierarchy()
  {
    for (size_t l=0; l!=m_tree.size(); ++l)
      for (typename std::map<T,node_t * >::const_iterator it = m_tree[l].begin(); it != m_tree[l].end(); ++it)
      {
        gsDebug<<"delete m_tree["<<l<<"]["<<it->first<<"])\n";
        delete it->second;
      }
  }

  /**
   * @brief      { function_description }
   *
   * @param[in]  times      The times (MONOTONICALLY INCREASING!)
   * @param[in]  solutions  The solutions
   * @param[in]  maxLevels  The maximum levels
   */
  gsSpaceTimeHierarchy(std::vector<T> times, std::vector<solution_t> solutions)
  :
  m_points(times.size()),
  m_initialized(false)
  {
    GISMO_ASSERT(times.size()==solutions.size(),"Sizes must agree");

    m_solutions.push_back(std::map<T,solution_t>());
    m_solutions.push_back(std::map<T,solution_t>());
    for (size_t k=0; k!=times.size(); ++k)
    {
      m_solutions[0].insert({times.at(k),solutions.at(k)});

    }

    _defaultOptions();
  }

  gsSpaceTimeHierarchy() {};

private:

  void _defaultOptions()
  {
    m_maxLevel = 5;
    m_split = false;
    m_options.addInt("MaxLevel","Sets the maximum level for hierarchical refinement",m_maxLevel);
    m_options.addInt("Split","Adds the jobs as sub-intervals, so that references with intermediate points can be computed. If this option is false (default) you have to perform TWO steps per item"
                      ,m_split);
  }

  void _applyOptions()
  {
    m_maxLevel = m_options.getInt("MaxLevel");
    m_split = m_options.getInt("Split");
    m_tree.resize(m_maxLevel+1);
  }

  std::pair<point_t,point_t> _interval(node_t * node)
  {
    std::pair<point_t,point_t> interval;

    if (node->isLeaf())
    {
      interval.first  = node->lowCorner();
      interval.second = node->uppCorner();
      return interval;
    }
    else
    {
      node_t * left  = node->left;
      node_t * right = node->right;
      while (!left->isLeaf() || !right->isLeaf())
      {
        if (!left->isLeaf())
          left  = left ->left;
        if (!right->isLeaf())
          right = right->right;
      }
      interval.first  = left->lowCorner();
      interval.second = right->uppCorner();
      return interval;
    }
  }

  // void _initLevel(index_t level)
  // {
  //   GISMO_ASSERT(m_ptimes.size()==m_indices.size(),"times and indices have different level");
  //   if (m_ptimes.size() < level+1 && level < m_maxLevel)
  //   {
  //     gsVector<T> time;
  //     time.setLinSpaced(m_points*math::pow(2,level),0,1);
  //     m_ptime = std::vector(time.data(),time.data() + time.rows()*time.cols());
  //     gsVector<index_t> index;
  //     index.setLinSpaced(m_points*math::pow(2,level),0,m_points*math::pow(2,level)-1);
  //     m_index = std::vector(index.data(),index.data() + index.rows()*index.cols());

  //     m_ptimes.push_back(std::map<index_t,T>());
  //     m_indices.push_back(std::map<T,index_t>());
  //     GISMO_ASSERT(m_ptime.size()==m_index.size(),"Sizes are not the same!");
  //     for (index_t k=0; k!=m_ptime.size(); k++)
  //     {
  //       m_ptimes[level].insert({m_index[k],m_ptime[k]});
  //       m_indices[level].insert({m_ptime[k],m_index[k]});
  //     }

  //     // m_solutions.push_back(std::map<index_t,solution_t>());
  //     // m_times.push_back(std::map<index_t,T>());
  //   }
  // }
public:

  void init()
  {
    this->_applyOptions();
    // Clean tree and queue
    m_tree.clear();
    std::queue<std::pair<index_t,T>> queue;
    m_queue.swap(queue);

    m_tree.resize(m_maxLevel+1);
    gsVector<T,1> low,upp;
    index_t k=0;
    for (typename std::map<T,solution_t>::const_iterator it = m_solutions[0].begin(); it != std::prev(m_solutions[0].end()); it++, k++)
    {
      /// Add parents to tree
      low<<it->first;
      upp<<std::next(it)->first;
      m_tree[0][it->first] = new node_t(low,upp);
    }

    for (typename std::map<T,node_t * >::const_iterator it = m_tree[0].begin(); it != m_tree[0].end(); ++it)
    {
      gsDebugVar(m_split);
      if (m_split)
      {
        low = it->second->lowCorner();
        upp = it->second->uppCorner();
        it->second->split(0,(upp[0]+low[0])/2);
        m_tree[1][it->second->left ->lowCorner()[0]] = it->second->left;
        m_tree[1][it->second->right->lowCorner()[0]] = it->second->right;

        /// Add childs to queue
        m_queue.push({1,it->second->left ->lowCorner()[0]});
        m_queue.push({1,it->second->right->lowCorner()[0]});
      }
      else
      {
        low = it->second->lowCorner();
        upp = it->second->uppCorner();

        m_tree[1][it->second->lowCorner()[0]] = it->second;
        /// Add childs to queue
        m_queue.push({1,it->second->lowCorner()[0]});
      }
    }

    /// Add childs to queue
    // for (typename std::map<T,node_t * >::const_iterator it = m_tree[1].begin(); it != m_tree[1].end(); ++it)
    //   m_queue.push({1,it->first});


    m_initialized = true;
  }

  gsOptionList & options() { return m_options; }

  std::tuple<index_t, T,    T,  solution_t, solution_t> pop()
  //         ID,   time, dt, start,      guess
  {
    GISMO_ASSERT(m_initialized,"Structure is not initialized");
    GISMO_ASSERT(!m_queue.empty(),"The queue is empty! Something went wrong.");

    T time = m_queue.front().second;
    index_t level = m_queue.front().first;
    m_queue.pop();

    GISMO_ASSERT(m_tree[level].count(time)>0,"Node at level "<<level<<" with time "<<time<<"not found");

    // GISMO_ASSERT(!m_tree[level][time]->isRoot(),"Node cannot be a parent!");

    gsVector<T> low, upp;
    std::tie(low,upp) = _interval(m_tree[level][time]);
    T tstart = low[0];
    T tend   = upp[0];
    T dt     = tend-tstart;

    solution_t start, guess;
    if (m_tree[level][time]->isLeftChild())
    {
      std::tie(low,upp) = _interval(m_tree[level][time]->parent);
      tend = upp(0);

      // Find best available start point
      index_t startLevel = level;
      while (startLevel > 0 && m_solutions[startLevel-1].count(tstart)==0)
        startLevel -=1;

      GISMO_ASSERT(m_solutions[startLevel-1].count(tstart)!=0,"Cannot find start point at tstart = "<<tstart<<" in level "<<startLevel-1);
      start = m_solutions[startLevel-1][tstart];

      // Find best available guess point
      index_t guessLevel = level;
      while (guessLevel > 0 && m_solutions[guessLevel-1].count(tend)==0)
        guessLevel -=1;

      GISMO_ASSERT(m_solutions[guessLevel-1].count(tend)!=0,"Cannot find end point at tend = "<<tend<<" in level "<<guessLevel-1);
      guess = m_solutions[guessLevel-1][tend];
    }
    else if (m_tree[level][time]->isRightChild())
    {
      std::tie(low,upp) = _interval(m_tree[level][time]->parent);
      tend = upp(0);
      GISMO_ASSERT(m_solutions[level].count(tstart)!=0,"Cannot find start point at tstart = "<<tstart<<" in level "<<level);
      start = m_solutions[level][tstart];
      GISMO_ASSERT(m_solutions[level-1].count(tend)!=0,"Cannot find end point at tend = "<<tend<<" in level "<<level-1);
      guess = m_solutions[level-1][tend];
    }
    else if (level>0)// m_tree[level][time]->isRoot()
    {
      std::tie(low,upp) = _interval(m_tree[level][time]);
      tend = upp(0);

      // Find best available start point
      index_t startLevel = level;
      while (startLevel > 0 && m_solutions[startLevel-1].count(tstart)==0)
        startLevel -=1;

      GISMO_ASSERT(m_solutions[startLevel-1].count(tstart)!=0,"Cannot find start point at tstart = "<<tstart<<" in level "<<startLevel-1);
      start = m_solutions[startLevel-1][tstart];

      // Find best available guess point
      index_t guessLevel = level;
      while (guessLevel > 0 && m_solutions[guessLevel-1].count(tend)==0)
        guessLevel -=1;

      GISMO_ASSERT(m_solutions[guessLevel-1].count(tend)!=0,"Cannot find end point at tend = "<<tend<<" in level "<<guessLevel-1);
      guess = m_solutions[guessLevel-1][tend];
    }
    else
      GISMO_ERROR("Node is a root!");

    // add the job to the active jobs
    m_jobs[m_ID++] = std::make_pair(level,time);
    gsDebug<<"Added active job (ID="<<m_ID-1<<") on level = "<<m_jobs[m_ID-1].first<<"; time = "<<m_jobs[m_ID-1].second<<"\n";

    return std::make_tuple(m_ID-1,tstart, dt, start, guess);
  }

  bool getReference(index_t ID, solution_t & result)
  {
    index_t level;
    T time;
    std::tie(level,time) = m_jobs[ID];
    gsVector<T> low, upp;
    std::tie(low,upp) = _interval(m_tree[level][time]);

    if (m_solutions[level-1].count(upp[0])!=0)
    {
      result =  m_solutions[level-1][upp[0]];
      return true;
    }
    else
      return false;
  }

  /**
   * @brief      Submits job \a ID to the solution database.
   *
   * The solution for job\a ID is added to the database
   *
   * @param[in]  ID        ID of the job
   * @param[in]  solution  The solution for the specific ID
   */
  void submit(index_t ID, solution_t solution)
  {
    GISMO_ASSERT(m_initialized,"Structure is not initialized");
    index_t level;
    T time;
    std::tie(level,time) = m_jobs[ID];
    node_t * tmp = m_tree[level][time];

    gsVector<T> low, upp;
    std::tie(low,upp) = _interval(tmp);
    T tend   = upp[0];

    if (m_solutions.size() < level + 1)
      m_solutions.resize(level+1);

    m_solutions[level][tend] = solution;
    gsDebug<<"Job submitted at t = "<<tend<<"; level "<<level<<"\n";
  }

  void submitLeft(index_t ID, solution_t solution)
  {
    GISMO_ASSERT(m_initialized,"Structure is not initialized");
    GISMO_ASSERT(!m_split,"This function is only useful when Split is false");
    gsVector<T> low, upp;
    index_t level;
    T time;
    std::tie(level,time) = m_jobs[ID];
    node_t * tmp = m_tree[level][time];

    std::tie(low,upp) = _interval(tmp);

    T tend   = (upp[0]+low[0])/2;

    if (m_solutions.size() < level + 1)
      m_solutions.resize(level+1);

    m_solutions[level][tend] = solution;
    gsDebug<<"Job submitted at t = "<<tend<<"; level "<<level<<"\n";
  }


  void submitRight(index_t ID, solution_t solution)
  {
    GISMO_ASSERT(m_initialized,"Structure is not initialized");
    GISMO_ASSERT(!m_split,"This function is only useful when Split is false");
    gsVector<T> low, upp;
    index_t level;
    T time;
    std::tie(level,time) = m_jobs[ID];
    node_t * tmp = m_tree[level][time];

    std::tie(low,upp) = _interval(tmp);

    T tend   = upp[0];

    if (m_solutions.size() < level + 1)
      m_solutions.resize(level+1);

    m_solutions[level][tend] = solution;
    gsDebug<<"Job submitted at t = "<<tend<<"; level "<<level<<"\n";
  }

  /**
   * @brief       Adds a job based on the job \a ID in the tree.
   *              Example: If current job is {0.2,4} on level l,
   *              the job that is added is XXXX
   *
   * @param[in]  ID    { parameter_description }
   */
  void addJob(index_t ID)
  {
    GISMO_ASSERT(m_initialized,"Structure is not initialized");
    index_t level;
    T time;
    std::tie(level,time) = m_jobs[ID];

    if (level+1 > m_maxLevel)
    {
      gsWarn<<"Max level reached!";
      return;
    }

    gsInfo<<"ID = "<<ID<<" level = "<<level<<" time = "<<time<<"\n";

    // Collects the time interval in the tree on \a level and \a time
    gsVector<T> low, upp;
    std::tie(low,upp) = _interval(m_tree[level][time]);

    // Splits the current time interval [low,upp].
    // m_tree[level][time]->split();
    m_tree[level][time]->split(0,(upp[0]+low[0])/2);

    // Add the new intervals to the tree
    m_tree[level+1][m_tree[level][time]->left ->lowCorner()[0]] = m_tree[level][time]->left;
    m_tree[level+1][m_tree[level][time]->right->lowCorner()[0]] = m_tree[level][time]->right;

    // Adds two new jobs on level \a l+1, for both new intervals. The jobs start at tje lower corner (i.e. the start of the intervals)
    m_queue.push({level+1,m_tree[level][time]->left ->lowCorner()[0]});
    m_queue.push({level+1,m_tree[level][time]->right->lowCorner()[0]});
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

  index_t currentLevel(index_t ID)
  {
    return _activeJob(ID).first;
  }

  bool empty()
  {
    return m_queue.empty();
  }

  std::pair<std::vector<T>,std::vector<solution_t>> getFlatSolution()
  {
    std::vector<T> times;
    std::vector<solution_t> solutions;

    std::map<T,solution_t> map;
    for (size_t l = 0; l!=m_solutions.size(); ++l)
      for (typename std::map<T,solution_t>::const_iterator it = m_solutions[l].begin(); it!=m_solutions[l].end(); ++it)
        map[it->first] = it->second; // overwrites if duplicate

    for (typename std::map<T,solution_t>::const_iterator it = map.begin(); it!=map.end(); ++it)
    {
      times.push_back(it->first);
      solutions.push_back(it->second);
    }

    return std::make_pair(times,solutions);
  }

  void print()
  {
    for (size_t l = 0; l!=m_solutions.size(); ++l)
    {
      gsInfo<<"level = "<<l<<":\n";
      for (typename std::map<T,solution_t>::const_iterator it = m_solutions[l].begin(); it!=m_solutions[l].end(); ++it)
      {
        gsInfo<<"\t";
        gsInfo<<"time = "<<it->first<<"\tsol = "<<it->second<<"\n";
      }
    }
  }

  void printTree()
  {
    gsVector<T,1> low, upp;
    for (size_t l=0; l!=m_tree.size(); ++l)
      for (typename std::map<T,node_t * >::const_iterator it = m_tree[l].begin(); it != m_tree[l].end(); ++it)
      {
        std::tie(low,upp) = _interval(it->second);
        if (l!=0)
          gsDebug<<"interval = ["<<low[0]<<","<<upp[0]<<"], level = "<<l<<"\n";
        else
          gsDebug<<"interval = ["<<low[0]<<","<<upp[0]<<"], level = "<<l<<"\n";
      }
  }

  void printQueue()
  {
    gsVector<T,1> low, upp;
    std::queue<std::pair<index_t,T>> queue = m_queue;
    gsInfo<<"Queue has "<<queue.size()<<" elements\n";
    while (!queue.empty())
    {
      gsInfo<<"Level "<<queue.front().first<<"; time "<<queue.front().second<<"\n";
      queue.pop();
    }
    gsInfo<<"Queue is empty\n";

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


protected:
  index_t m_points;
  index_t m_maxLevel;
  bool m_split;
  bool m_initialized;

  std::vector<std::map<T,solution_t>>       m_solutions;

  std::vector<std::map<T,node_t *>>         m_tree;

  std::queue<std::pair<index_t,T>>          m_queue;

  std::map<index_t,std::pair<index_t,T>>    m_jobs;
  mutable index_t                           m_ID = (0);

  gsOptionList m_options;

};

}

// template<class T, class solution_t >
// class gsSpaceTimeHierarchy
// {
// public:

//   typedef kdnode<1, T> node_t;

//   typedef typename node_t::point point_t;

// public:

//   ~gsSpaceTimeHierarchy()
//   {
//     for (size_t l=0; l!=m_tree.size(); ++l)
//       for (typename std::map<T,node_t * >::const_iterator it = m_tree[l].begin(); it != m_tree[l].end(); ++it)
//       {
//         gsDebug<<"delete m_tree["<<l<<"]["<<it->first<<"])\n";
//         delete it->second;
//       }
//   }

//   /**
//    * @brief      { function_description }
//    *
//    * @param[in]  times      The times (MONOTONICALLY INCREASING!)
//    * @param[in]  solutions  The solutions
//    * @param[in]  maxLevels  The maximum levels
//    */
//   gsSpaceTimeHierarchy(std::vector<T> times, std::vector<solution_t> solutions, index_t maxLevels = 5)
//   :
//   m_points(times.size()),
//   m_maxLevel(maxLevels),
//   m_ID(0)
//   {
//     GISMO_ASSERT(times.size()==solutions.size(),"Sizes must agree");

//     m_solutions.push_back(std::map<T,solution_t>());
//     m_solutions.push_back(std::map<T,solution_t>());
//     for (size_t k=0; k!=times.size(); ++k)
//     {
//       m_solutions[0].insert({times.at(k),solutions.at(k)});

//     }

//     _init();

//   }

// private:
//   void _init();

//   std::pair<point_t,point_t> _interval(node_t * node);

//   void _defaultOptions()
//   {
//     m_maxLevel = 5;
//     m_options.addInt("MaxLevel","Sets the maximum level for hierarchical refinement",m_maxLevel);
//   }

//   void _applyOptions()
//   {
//     m_maxLevel = m_options.getInt("MaxLevel");
//     m_tree.resize(m_maxLevel+1);

//   }

// public:

//   gsOptionList & options() { return m_options; }

//   std::tuple<index_t, T,    T,  solution_t, solution_t> pop();

//   bool getReference(index_t ID, solution_t & result);


//   /**
//    * @brief      Submits job \a ID to the solution database.
//    *
//    * The solution for job\a ID is added to the database and the job ID is removed from the
//    *
//    * @param[in]  ID        ID of the job
//    * @param[in]  solution  The solution for the specific ID
//    */
//   void submit(index_t ID, solution_t solution);


//   /**
//    * @brief       Adds a job based on the job \a ID in the tree.
//    *              Example: If current job is {0.2,4} on level l,
//    *              the job that is added is XXXX
//    *
//    * @param[in]  ID    { parameter_description }
//    */
//   void addJob(index_t ID);

//   index_t currentLevel(index_t ID);

//   bool empty()
//   {
//     return m_queue.empty();
//   }

//   std::pair<std::vector<T>,std::vector<solution_t>> getFlatSolution();

//   void print();

//   void printTree();

//   void printQueue();

// protected:
//   std::pair<index_t,T> _activeJob(index_t ID);


// protected:
//   index_t m_points;
//   index_t m_maxLevel;

//   std::vector<std::map<T,solution_t>>       m_solutions;

//   std::vector<std::map<T,node_t *>>         m_tree;

//   std::queue<std::pair<index_t,T>>          m_queue;

//   std::map<index_t,std::pair<index_t,T>>    m_jobs;
//   mutable index_t                           m_ID;

//   gsOptionList m_options;

// };

// } // namespace gismo

// #ifndef GISMO_BUILD_LIB
// #include GISMO_HPP_HEADER(gsSpaceTimeHierarchy.hpp)
// #endif
