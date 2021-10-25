/** @file gsSpaceTimeHierarchy.hpp

    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

namespace gismo
{

  template<class T, class solution_t >
  void gsSpaceTimeHierarchy<T,solution_t>::_init()
  {
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

        m_tree[1][it->second->lowCorner()[0]] = it->second->left;
        /// Add childs to queue
        m_queue.push({1,it->second->lowCorner()[0]});
      }
    }

    /// Add childs to queue
    // for (typename std::map<T,node_t * >::const_iterator it = m_tree[1].begin(); it != m_tree[1].end(); ++it)
    //   m_queue.push({1,it->first});
  }

  template<class T, class solution_t >
  std::pair<point_t,point_t> gsSpaceTimeHierarchy<T,solution_t>::_interval(node_t * node)
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

  template<class T, class solution_t >
  std::tuple<index_t, T,    T,  solution_t, solution_t>
  gsSpaceTimeHierarchy<T,solution_t>::pop()
{
    this->_applyOptions();

    GISMO_ASSERT(!m_queue.empty(),"The queue is empty! Something went wrong.");

    T time = m_queue.front().second;
    index_t level = m_queue.front().first;

    GISMO_ASSERT(m_tree[level].count(time)>0,"Node at level "<<level<<" with time "<<time<<"not found");

    GISMO_ASSERT(!m_tree[level][time]->isRoot(),"Node cannot be a parent!");

    gsVector<T> low, upp;
    std::tie(low,upp) = _interval(m_tree[level][time]);
    T tstart = low[0];
    T tend   = upp[0];
    T dt     = tend-tstart;

    gsDebugVar(level);

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

      // Find best available start point
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
    else
      GISMO_ERROR("Node is not a child!");

    // add the job to the active jobs
    gsDebugVar(level);
    gsDebugVar(time);
    m_jobs[m_ID++] = std::make_pair(level,time);
    gsDebug<<"Added active job (ID="<<m_ID-1<<") on level = "<<m_jobs[m_ID-1].first<<"; time = "<<m_jobs[m_ID].second<<"\n";

    return std::make_tuple(m_ID-1,tstart, dt, start, guess);
  }

  template<class T, class solution_t >
  bool gsSpaceTimeHierarchy<T,solution_t>::getReference(index_t ID, solution_t & result)
  {
    index_t level;
    T time;
    std::tie(level,time) = m_jobs[ID];
    gsDebugVar(level);
    gsDebugVar(time);

    gsDebugVar(m_queue.front().second);
    gsDebugVar(m_queue.front().first);

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

  template<class T, class solution_t >
  void gsSpaceTimeHierarchy<T,solution_t>::submit(index_t ID, solution_t solution)
  {
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

    // Remove ID
    gsDebug<<"Erasing finished job on level = "<<m_jobs[m_ID].first<<"; time = "<<m_jobs[m_ID].second<<"\n";
    m_jobs.erase(ID);
    gsDebug<<"Current active jobs\n";
    for (typename std::map<index_t,std::pair<index_t,T>>::const_iterator it = m_jobs.begin(); it!=m_jobs.end(); ++it)
      gsDebug<<"level = "<<it->first<<"; time = "<<it->second<<"\n";
  }

  template<class T, class solution_t >
  void gsSpaceTimeHierarchy<T,solution_t>::addJob(index_t ID)
  {
    index_t level;
    T time;
    std::tie(level,time) = m_jobs[ID];

    if (level+1 > m_maxLevel)
    {
      gsWarn<<"Max level reached!";
      return;
    }

    // Collects the time interval in the tree on \a level and \a time
    gsVector<T> low, upp;
    std::tie(low,upp) = _interval(m_tree[level][time]);

    // Splits the current time interval [low,upp].
    m_tree[level][time]->split(0,(upp[0]+low[0])/2);

    // Add the new intervals to the tree
    m_tree[level+1][m_tree[level][time]->left ->lowCorner()[0]] = m_tree[level][time]->left;
    m_tree[level+1][m_tree[level][time]->right->lowCorner()[0]] = m_tree[level][time]->right;

    gsDebugVar(m_tree[level][time]->left ->lowCorner()[0]);
    gsDebugVar(m_tree[level][time]->right->lowCorner()[0]);

    // Adds two new jobs on level \a l+1, for both new intervals. The jobs start at tje lower corner (i.e. the start of the intervals)
    m_queue.push({level+1,m_tree[level][time]->left ->lowCorner()[0]});
    m_queue.push({level+1,m_tree[level][time]->right->lowCorner()[0]});
  }

  template<class T, class solution_t >
  index_t gsSpaceTimeHierarchy<T,solution_t>::currentLevel(index_t ID)
  {
    return _activeJob(ID).first;
  }

  template<class T, class solution_t >
  std::pair<std::vector<T>,std::vector<solution_t>> gsSpaceTimeHierarchy<T,solution_t>::getFlatSolution()
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

  template<class T, class solution_t >
  void gsSpaceTimeHierarchy<T,solution_t>::print()
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

  template<class T, class solution_t >
  void gsSpaceTimeHierarchy<T,solution_t>::printTree()
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


  template<class T, class solution_t >
  void gsSpaceTimeHierarchy<T,solution_t>::printQueue()
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


  template<class T, class solution_t >
  std::pair<index_t,T> gsSpaceTimeHierarchy<T,solution_t>::_activeJob(index_t ID)
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


} // namespace gismo