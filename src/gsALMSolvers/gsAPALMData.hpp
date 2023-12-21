/** @file gsAPALMData.hpp

    @brief XXXX

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#pragma once

namespace gismo
{

template <class T, class solution_t >
gsAPALMData<T,solution_t>::gsAPALMData( const std::vector<T> & times,
                                  const  std::vector<solution_t> & solutions)
:
m_initialized(false),
m_dt(1.0)
{
  GISMO_ASSERT(times.size()==solutions.size(),"Sizes must agree");

  setData(times,solutions);

  _defaultOptions();
}

template <class T, class solution_t >
gsAPALMData<T,solution_t>::gsAPALMData(const std::vector<T> & times)
:
m_initialized(false),
m_dt(1.0)
{
  // Initialize the knot vectors
  m_t = gsKnotVector<T>(times);
  m_xi = m_t;
  m_xi.transform(0,1);

  _defaultOptions();
}

template <class T, class solution_t >
gsAPALMData<T,solution_t>::gsAPALMData()
:
m_initialized(false),
m_dt(1.0)
{
  _defaultOptions();
  m_t = gsKnotVector<T>(0);
}

template <class T, class solution_t >
void gsAPALMData<T,solution_t>::_defaultOptions()
{
  m_maxLevel = 5;
  m_tolerance = 1e-1;
  m_verbose = 0;
  m_options.addInt("MaxLevel","Sets the maximum level for hierarchical refinement",m_maxLevel);
  m_options.addReal("Tolerance","Relative tolerance",m_tolerance);
  m_options.addInt("Verbose","Verbosity; 0=none, 1=minimal, 2=full",m_verbose);
}

template <class T, class solution_t >
void gsAPALMData<T,solution_t>::_applyOptions()
{
  m_maxLevel = m_options.getInt("MaxLevel");
  m_tolerance = m_options.getReal("Tolerance");
  m_verbose = m_options.getInt("Verbose");
}


template <class T, class solution_t >
void gsAPALMData<T,solution_t>::init()
{
  this->_applyOptions();
  // Clean queue
  std::deque<std::tuple<T,T,index_t>> queue;
  m_queue.swap(queue);

  T low,upp;
  for (typename gsKnotVector<T>::const_iterator it = m_xi.begin(); it != std::prev(m_xi.end()); )
  {
    low = *it;
    it++;
    upp = *it;
    m_queue.push_back({low,upp,1});
  }

  this->_buildMap();

  m_initialized = true;
}

template <class T, class solution_t >
void gsAPALMData<T,solution_t>::addStartPoint(const T & time, const solution_t & solution, bool priority)
{
  GISMO_ASSERT(m_t.size()==m_xi.size(),"Sizes must be the same!");
  GISMO_ASSERT(m_t.size()==0,"Start point can only be added in a new branch");
  this->_applyOptions();
  T xi  = 0;
  m_t.insert(time);
  m_xi.insert(xi);

  auto sol = std::make_shared<solution_t>(solution);
  m_solutions.insert({xi,sol});
  m_levels.insert({xi,0});

  auto prev = std::make_shared<solution_t>(*(m_solutions[xi]).get());
  m_prevs  .insert({xi,prev});

  // same start and end coordinate, will be recognized by the pop function
  if (priority)
    m_queue.push_front({xi,xi,0});
  else
    m_queue.push_back({xi,xi,0});
  m_initialized = true;
}

template <class T, class solution_t >
void gsAPALMData<T,solution_t>::appendPoint(bool priority)
{
  // same start and end coordinate, will be recognized by the pop function
  if (priority)
    m_queue.push_front({m_xi.last(),m_xi.last(),0});
  else
    m_queue.push_back({m_xi.last(),m_xi.last(),0});
  m_initialized = true;
}

template <class T, class solution_t >
void gsAPALMData<T,solution_t>::appendData(const T & time, const solution_t & solution, bool priority)
{
  // Initialize the knot vectors
  GISMO_ASSERT(m_t.size()==m_xi.size(),"Sizes must be the same!");
  GISMO_ASSERT(m_t.size()==0 || time>m_t.last(),"Time must be bigger than last stored time!");
  T xi;
  if (m_t.size()==0)
    xi = 0;
  else if (m_t.size()==1)
    xi = 1;
  else
    xi = (time-m_t.last()) * (m_xi.last()-m_xi.first()) / (m_t.last()-m_t.first()) + m_xi.last();

  auto sol = std::make_shared<solution_t>(solution);
  m_solutions.insert({xi,sol});
  m_levels.insert({xi,0});

  // add a previous, if exists
  if (m_t.size()==0)
  {
    auto prev = std::make_shared<solution_t>(*(m_solutions[xi]).get());
    m_prevs  .insert({xi,prev});
  }
  else
  {
    auto prev = std::make_shared<solution_t>(*(m_solutions[m_xi.last()]).get());
    m_prevs  .insert({xi,prev});
    if (priority)
      m_queue.push_front({m_xi.last(),xi,1});
    else
      m_queue.push_back({m_xi.last(),xi,1});
    m_initialized = true;
  }


  m_t.insert(time);
  m_xi.insert(xi);

  this->_buildMap();
}

template <class T, class solution_t >
void gsAPALMData<T,solution_t>::setData(const std::vector<T> & times,const  std::vector<solution_t> & solutions)
{
  // Initialize the knot vectors
  m_t = gsKnotVector<T>(times);
  m_xi = m_t;
  m_xi.transform(0,1);

  auto sol0 = std::make_shared<solution_t>(solutions.at(0));
  m_solutions.insert({m_xi.at(0),sol0});
  m_levels.insert({m_xi.at(0),0});
  auto prev0 = std::make_shared<solution_t>(*(m_solutions[m_xi.at(0)]).get());
  m_prevs  .insert({m_xi.at(0),prev0});
  for (size_t k=1; k!=m_xi.size(); ++k)
  {
    auto solk = std::make_shared<solution_t>(solutions.at(k));
    m_solutions.insert({m_xi.at(k),solk});
    m_levels.insert({m_xi.at(k),0});
    auto prevk = std::make_shared<solution_t>(*(m_solutions[m_xi.at(k-1)]).get());
    m_prevs .insert({m_xi.at(k),prevk});
  }
}

// template <class T, class solution_t >
// void gsAPALMData<T,solution_t>::initEmptyQueue(const index_t & N)
// {
//   setData(times,solutions);
//   init();
// }

template <class T, class solution_t >
void gsAPALMData<T,solution_t>::init(const std::vector<T> & times,const  std::vector<solution_t> & solutions)
{
  setData(times,solutions);
  init();
}

template <class T, class solution_t >
std::tuple<index_t, T     , solution_t, solution_t> gsAPALMData<T,solution_t>::pop()
{
  GISMO_ASSERT(m_initialized,"Structure is not initialized");
  GISMO_ASSERT(!m_queue.empty(),"The queue is empty! Something went wrong.");

  T xilow = std::get<0>(m_queue.front());
  T xiupp = std::get<1>(m_queue.front());
  index_t level = std::get<2>(m_queue.front());
  m_queue.pop_front();

  T tlow = m_tmap[xilow];
  T tupp = m_tmap[xiupp];

  T dt;
  if (xilow==xiupp && level==0)
    dt = m_dt;
  else
    dt = (tupp-tlow);

  GISMO_ASSERT(m_solutions.count(xilow)!=0,"Cannot find start point at tstart = "<<tlow<<"\n");
  GISMO_ASSERT(m_prevs  .count(xilow)!=0,"Cannot find previous solution for tstart = "<<tlow<<"\n");

  solution_t start = *(m_solutions[xilow].get());

  solution_t prev = *(m_prevs[xilow].get());

  // add the job to the active jobs
  m_jobs[m_ID++] = std::make_tuple(xilow,xiupp,level);
  if (m_verbose==2) gsInfo<<"Got active job (ID="<<m_ID-1<<") on interval = ["<<xilow<<","<<xiupp<<"] = ["<<tlow<<","<<tupp<<"] with level "<<level<<"\n";

  for (typename std::map<index_t,std::tuple<T,T,index_t>>::const_iterator it = m_jobs.begin(); it!=m_jobs.end(); it++)
  {
    std::tie(xilow,xiupp,level) = it->second;
    gsDebug<<"job "<<it->first<<" on ["<<xilow<<","<<xiupp<<"] = ["<<tlow<<","<<tupp<<"] with level "<<level<<"\n";
  }

  return std::make_tuple(m_ID-1, dt, start, prev);
}

template <class T, class solution_t >
bool gsAPALMData<T,solution_t>::getReferenceByTime(T time, solution_t & result)
{
  if (m_solutions.count(m_ximap[time])!=0)
  {
    result =  *(m_solutions[m_ximap[time]].get());
    return true;
  }
  else
    return false;
}

template <class T, class solution_t >
bool gsAPALMData<T,solution_t>::getReferenceByPar(T xi, solution_t & result)
{
  if (m_solutions.count(xi)!=0)
  {
    result =  *(m_solutions[xi].get());
    return true;
  }
  else
    return false;
}

template <class T, class solution_t >
bool gsAPALMData<T,solution_t>::getReferenceByID(index_t ID, solution_t & result)
{
  if (m_solutions.count(std::get<1>(m_jobs[ID]))!=0)
  {
    result =  *(m_solutions[std::get<1>(m_jobs[ID])].get());
    return true;
  }
  else
    return false;
}

template <class T, class solution_t >
void gsAPALMData<T,solution_t>::submit(index_t ID, const std::vector<T> & distances, std::vector<solution_t> solutions, const T & upperDistance, const T & lowerDistance)
{
  GISMO_ASSERT(m_initialized,"Structure is not initialized");
  GISMO_ASSERT(distances.size()==solutions.size()+1,"You must provide one more distance than solutions");

  //  <d(xilow,xin)                                            >
  //  <d(xilow,xi1)> <d(xi1,xi2) > <d(xi2,Dxin-1)> <d(xin-1,xin)> <d(xin,xiupp)>
  // |--------------|-------------|..............|--------------|--------------|
  // xilow         xi1           xi2           xin-1           xin           xiupp

  // All distances d(.,.) are measured in the solution space (not the parametric space!).
  // D(.,.) is the sum of intervals
  // By the triangle inequality, d(xilow,xin)\leq D(xilow,xin)

  // N: number of interior points
  // Xi:        [xilow,xi1,...,xiupp] (size N+2)
  // distances: [d(xilow,xi1),d(xi1,xi2),...,d(xin-1,xin),d(xin,xiupp)] (size N+1)
  // solutions: [u(xi1),...,u(xin)] (size N)
  // lowerDistance = d(xilow,xin)
  // upperDistance = d(xilow,xiupp)

  //////////////////////////////////////////////////////////////////////////////
  // Rescale intervals
  //////////////////////////////////////////////////////////////////////////////
  T tlow, tupp, xilow, xiupp, dxi, Dt, dt;
  std::vector<T> t,xi;
  index_t level;

  // Compute interval (xi and t)
  std::tie(xilow,xiupp,level) = m_jobs[ID];
  tlow = m_tmap[xilow];
  tupp = m_tmap[xiupp];

  // Compute interval distance
  dxi = xiupp - xilow;
  // Get the original distance
  Dt = tupp - tlow; // original time step

  // Get time points
  t.resize(distances.size()+1);
  t.at(0) = tlow;
  for (size_t k=0; k!=distances.size(); k++)
    t.at(k+1) = t.at(k) + distances.at(k);

  // check if the total distance matches the original time step
  T dt_tmp = std::accumulate(distances.begin(), std::prev(distances.end()), 0.0);
  GISMO_ASSERT((Dt-dt_tmp)/Dt<1e-12,"Total distance of the computed intervals should be equal to the original interval length ("<<dt_tmp<<"!="<<Dt<<")");

  // get the total travelled distance
  dt = t.back()-tlow;

  // Compute the lower and upper errors
  // The lowerError is the surplus distance over the first intervals
  // The upperError is the distance between the last computed point and the referemce, minus the lowerError
  T totalError = dt-upperDistance;
  T lowerError = Dt-lowerDistance;
  T upperError = totalError-lowerError;

  // T lowerError = Dt-lowerDistance;
  // T upperError = upperDistance-Dt;
  // T totalError = upperError+lowerError;
  // gsDebugVar(totalError);
  // gsDebugVar(lowerError);
  // gsDebugVar(upperError);

  // Fill and scale xi vector
  xi.resize(solutions.size()+2);
  xi.front() = xilow;
  xi.back() = xiupp;
  for (size_t k=1; k!=xi.size()-1; k++)
    xi.at(k) = xilow + dxi * ( t.at(k) - t.front() ) / dt;

  // If the upper error is small enough, we don't store the solution
  // Also in this case, the lower error is likely to be small as well

  if (m_verbose==1)
  {
    gsInfo<<"Relative errors:\n";
    gsInfo<<"\tlowerError/dt = "<<lowerError/dt<<"\n";
    gsInfo<<"\tupperError/dt = "<<upperError/dt<<"\n";
  }
  // if (totalError/dt < m_tolerance)
  if (false)
  {
    GISMO_ASSERT(lowerError/dt<m_tolerance,"Lower error is larger than total error. How is this possible?");
    t .erase(t .end()-2);
    xi.erase(xi.end()-2);
    solutions.erase(solutions.end()-1);

    // IS THIS CORRECT????
  }
  else
  {
    size_t kmin = 1;
    size_t kmax = xi.size();
    if (lowerError/Dt<m_tolerance)
      kmin = xi.size()-1;
    if (upperError/Dt<m_tolerance)
      kmax = xi.size()-1;
    for (size_t k = kmin; k!=kmax; k++)
    {
      if (level < m_maxLevel)
      {
        m_queue.push_back(std::make_tuple(xi.at(k-1),xi.at(k),level+1));
        if (m_verbose==1) gsInfo<<"Interval ["<<xi.at(k-1)<<","<<xi.at(k)<<"] on level "<<level+1<<" added to queue\n";
      }
      else
      {
        if (m_verbose==1) gsInfo<<"Interval ["<<xi.at(k-1)<<","<<xi.at(k)<<"] on level "<<level+1<<" NOT added to queue (max level reached)\n";
      }

    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // Push the data
  //////////////////////////////////////////////////////////////////////////////

  // Transform the time vector
  m_t.addConstant(tupp,dt-Dt);

  for (size_t k=1; k!=t.size()-1; k++)
  {
      m_t.insert(t.at(k));
      m_xi.insert(xi.at(k));
  }

  this->_buildMap();

  for (size_t k=1; k!=xi.size()-1; k++) // add only INTERIOR solutions
  {
    if (m_verbose==2) gsInfo<<"Added a solution on time = "<<m_tmap[xi.at(k)]<<" (parametric time = "<<xi.at(k)<<")\n";
    m_solutions[xi.at(k)] = std::make_shared<solution_t>(solutions.at(k-1));
    m_levels[xi.at(k)] = level;
  }

  // The previous solution for the first computed point (xi[1]) is the solution on xi[0]
  for (size_t k=1; k!=xi.size()-1; k++) // add only INTERIOR solutions
  {
    if (m_verbose==2) gsInfo<<"Added a previous solution on time = "<<m_tmap[xi.at(k)]<<" (parametric time = "<<xi.at(k)<<")\n";
    m_prevs[xi.at(k)] = std::make_shared<solution_t>(*(m_solutions[xi.at(k-1)].get()));
  }
}

template <class T, class solution_t >
T gsAPALMData<T,solution_t>::jobStartTime(index_t ID)
{
  return m_tmap[std::get<0>(m_jobs[ID])];
}

template <class T, class solution_t >
T gsAPALMData<T,solution_t>::jobStartPar(index_t ID)
{
  return std::get<0>(m_jobs[ID]);
}

template <class T, class solution_t >
std::pair<T,T> gsAPALMData<T,solution_t>::jobTimes(index_t ID)
{
  T xilow,xiupp;
  std::tie(xilow,xiupp,std::ignore) = m_jobs[ID];
  return std::make_pair(m_tmap[xilow],m_tmap[xiupp]);

}

template <class T, class solution_t >
std::pair<T,T> gsAPALMData<T,solution_t>::jobPars(index_t ID)
{
  T xilow,xiupp;
  std::tie(xilow,xiupp,std::ignore) = m_jobs[ID];
  return std::make_pair(xilow,xiupp);
}

template <class T, class solution_t >
index_t gsAPALMData<T,solution_t>::jobLevel(index_t ID)
{
  return std::get<2>(m_jobs[ID]);
}

template <class T, class solution_t >
void gsAPALMData<T,solution_t>::printActiveJobs()
{
  for (typename std::map<index_t,std::tuple<T,T,index_t>>::const_iterator it = m_jobs.cbegin(); it!= m_jobs.cend(); it++)
  {
    gsInfo<<"ID = "<<it->first<<": xilow = "<<std::get<0>(it->second)<<"; xiupp = "<<std::get<1>(it->second)<<"; level = "<<std::get<2>(it->second)<<"\n";
  }
}

template <class T, class solution_t >
void gsAPALMData<T,solution_t>::finishJob(index_t ID)
{
  // Remove ID
  if (m_verbose==2) gsInfo<<"Erasing finished job on level = "<<std::get<0>(m_jobs[ID])
                          <<"; time = "<<std::get<1>(m_jobs[ID])<<"\n";
  m_jobs.erase(ID);
  if (m_verbose==2)
  {
    gsDebug<<"Current active jobs\n";
    for (typename std::map<index_t,std::tuple<T,T,index_t>>::const_iterator it = m_jobs.begin(); it!=m_jobs.end(); ++it)
      gsDebug<<"ID = "<<it->first<<"; xilow = "<<std::get<0>(it->second)<<"; xiupp = "<<std::get<1>(it->second)<<"; level = "<<std::get<2>(it->second)<<"\n";
  }
}

template <class T, class solution_t >
bool gsAPALMData<T,solution_t>::empty()
{
  return m_queue.empty();
}

template <class T, class solution_t >
std::tuple<std::vector<T>,std::vector<solution_t>,std::vector<index_t>> gsAPALMData<T,solution_t>::getFlatSolution(index_t level)
{
  std::vector<T> times;
  std::vector<solution_t> solutions;
  std::vector<index_t> levels;

  for (typename std::map<T,std::shared_ptr<solution_t>>::const_iterator
        it = m_solutions.begin();
        it!=m_solutions.end();
        ++it)
  {
    if (m_levels[it->first]==level || level==-1)
    {
      times.push_back(m_tmap[it->first]);
      solutions.push_back(*(it->second.get()));
      levels.push_back(m_levels[it->first]);
    }
  }
  return std::make_tuple(times,solutions,levels);
}

template <class T, class solution_t >
void gsAPALMData<T,solution_t>::print()
{
  for (typename std::map<T,std::shared_ptr<solution_t>>::const_iterator
        it = m_solutions.begin();
        it!=m_solutions.end();
        ++it)
  {
    gsInfo<<"\t";
    gsInfo<<"time = "<<it->first<<"\n";
  }
}

template <class T, class solution_t >
void gsAPALMData<T,solution_t>::printQueue()
{
  std::deque<std::tuple<T,T,index_t>> queue = m_queue;
  gsInfo<<"Queue has "<<queue.size()<<" elements\n";
  while (!queue.empty())
  {
    gsInfo<<"Start time "<<std::get<0>(queue.front())
          <<"; end time = "<<std::get<1>(queue.front())
          <<"; level = "<<std::get<2>(queue.front())<<"\n";
    queue.pop_front();
  }
  gsInfo<<"Queue is empty\n";
}

template <class T, class solution_t >
void gsAPALMData<T,solution_t>::printKnots()
{
  gsInfo<<"xi = \n"<<m_xi.asMatrix()<<"\n";
  gsInfo<<"t = \n"<<m_t.asMatrix()<<"\n";
}

template <class T, class solution_t >
std::tuple<T,T,index_t> gsAPALMData<T,solution_t>::_activeJob(index_t ID)
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

template <class T, class solution_t >
void gsAPALMData<T,solution_t>::_buildMap()
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

} // namespace gismo