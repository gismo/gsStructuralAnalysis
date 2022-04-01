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
m_initialized(false)
{
  GISMO_ASSERT(times.size()==solutions.size(),"Sizes must agree");

  setData(times,solutions);

  _defaultOptions();
}

template <class T, class solution_t >
gsAPALMData<T,solution_t>::gsAPALMData(const std::vector<T> & times)
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

template <class T, class solution_t >
gsAPALMData<T,solution_t>::gsAPALMData()
:
m_initialized(false)
{
  _defaultOptions();
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
  std::queue<std::pair<T,T>> queue;
  m_queue.swap(queue);

  T low,upp;
  for (typename gsKnotVector<T>::const_iterator it = m_xi.begin(); it != std::prev(m_xi.end()); )
  {
    low = *it;
    it++;
    upp = *it;
    m_queue.push({low,upp});
  }

  this->_buildMap();

  m_initialized = true;
}

template <class T, class solution_t >
void gsAPALMData<T,solution_t>::setData(const std::vector<T> & times,const  std::vector<solution_t> & solutions)
{
  // Initialize the knot vectors
  m_t = gsKnotVector<T>(times,1,0);
  m_t.degreeDecrease(1);

  m_xi = m_t;
  m_xi.transform(0,1);

  m_solutions.insert({m_xi.at(0),solutions.at(0)});
  m_guesses  .insert({m_xi.at(0),&(m_solutions[m_xi.at(0)])});
  for (size_t k=1; k!=m_xi.size(); ++k)
  {
    m_solutions.insert({m_xi.at(k),solutions.at(k)});
    m_guesses  .insert({m_xi.at(k),&(m_solutions[m_xi.at(k)])});
  }
}

template <class T, class solution_t >
void gsAPALMData<T,solution_t>::init(const std::vector<T> & times,const  std::vector<solution_t> & solutions)
{
  setData(times,solutions);
  init();
}

template <class T, class solution_t >
std::tuple<index_t, T     , T   , T , solution_t, solution_t> gsAPALMData<T,solution_t>::pop()
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
  GISMO_ASSERT(m_guesses  .count(xilow)!=0,"Cannot find end point at tend = "<<tupp<<"\n");

  solution_t start = m_solutions[xilow];
  solution_t guess = *(m_guesses[xilow]);

  // add the job to the active jobs
  m_jobs[m_ID++] = std::make_pair(xilow,xiupp);
  if (m_verbose==2) gsInfo<<"Got active job (ID="<<m_ID-1<<") on interval = ["<<xilow<<","<<xiupp<<"] = ["<<tlow<<","<<tupp<<"]\n";

  return std::make_tuple(m_ID-1,tlow, tupp, dt, start, guess);
}

template <class T, class solution_t >
bool gsAPALMData<T,solution_t>::getReferenceByTime(T time, solution_t & result)
{
  if (m_solutions.count(m_ximap[time])!=0)
  {
    result =  m_solutions[m_ximap[time]];
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
    result =  m_solutions[xi];
    return true;
  }
  else
    return false;
}

template <class T, class solution_t >
bool gsAPALMData<T,solution_t>::getReferenceByID(index_t ID, solution_t & result)
{
  if (m_solutions.count(m_jobs[ID].second)!=0)
  {
    result =  m_solutions[m_jobs[ID].second];
    return true;
  }
  else
    return false;
}

template <class T, class solution_t >
void gsAPALMData<T,solution_t>::submit(index_t ID, const std::vector<T> & distances, std::vector<solution_t> solutions, const T & upperError, const T & lowerError)
{
  GISMO_ASSERT(m_initialized,"Structure is not initialized");
  GISMO_ASSERT(distances.size()==solutions.size(),"You must provide one distance per solution");

  //  <d(xilow,xin)                                            >
  //  <d(xilow,xi1)> <d(xi1,xi2) > <d(xi2,Dxin-1)> <d(xin-1,xin)> <d(xin,xiupp)>
  // |--------------|-------------|..............|--------------|--------------|
  // xilow        xi1          xi2          xin-1        xin          xiupp

  // All distances d(.,.) are measured in the solution space (not the parametric space!).
  // D(.,.) is the sum of intervals
  // By the triangle inequality, d(xilow,xin)\leq D(xilow,xin)

  // Xi:        [xilow,xi1,...,xiupp]
  // distances: [d(xilow,xi1),d(xi1,xi2),...,d(xin-1,xin)]
  // solutions: [u(xi1),...,u(xin)]
  // lowerError = d(xilow,xin)
  // upperError = d(xin,xiupp)

  //////////////////////////////////////////////////////////////////////////////
  // Rescale intervals
  //////////////////////////////////////////////////////////////////////////////
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
  t.resize(solutions.size()+1);
  t.at(0) = tlow;
  for (size_t k=0; k!=solutions.size(); k++)
    t.at(k+1) = t.at(k) + distances.at(k);
  t.push_back(t.back()+upperError);

  T dt_tmp = 0;
  for(typename std::vector<T>::const_iterator it = distances.begin(); it != distances.end(); ++it)
      dt_tmp += *it;
  GISMO_ASSERT((Dt-dt_tmp)/Dt<1e-12,"Total distance of the computed intervals should be equal to the original interval length ("<<dt_tmp<<"!="<<Dt<<")");

  // Get actual time step and the error
  dt = t.back()-tlow;

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
    gsInfo<<"\tlowerError/Dt = "<<lowerError/Dt<<"\n";
    gsInfo<<"\tupperError/Dt = "<<upperError/Dt<<"\n";
  }
  if (upperError/Dt < m_tolerance)
  {
    GISMO_ASSERT(lowerError/Dt<m_tolerance,"Lower error is big, but upper error is small. How is this possible?");
    t .erase(t .end()-2);
    xi.erase(xi.end()-2);
    solutions.erase(solutions.end()-1);

    // IS THIS CORRECT????
  }
  else
  {
    size_t kstart = xi.size()-1;
    if (lowerError/Dt<m_tolerance)
      kstart = 1;
    for (size_t k = kstart; k!=xi.size(); k++)
    {
      m_queue.push(std::make_pair(xi.at(k-1),xi.at(k)));
      if (m_verbose==1) gsInfo<<"Interval ["<<xi.at(k-1)<<","<<xi.at(k)<<"] added to queue\n";
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
    m_solutions[xi.at(k)] = solutions.at(k-1);
  }

  // The guess for the first computed point (xi[1]) is the solution on xi[0]
  m_guesses[xi.at(1)] = &(m_solutions[xi.at(0)]);
  for (size_t k=2; k!=xi.size(); k++) // add only INTERIOR solutions
  {
    if (m_verbose==2) gsInfo<<"Added a guess on time = "<<m_tmap[xi.at(k)]<<" (parametric time = "<<xi.at(k)<<")\n";
    m_guesses[xi.at(k)] = &(m_solutions[xi.at(k-1)]);
  }
}

template <class T, class solution_t >
void gsAPALMData<T,solution_t>::finishJob(index_t ID)
{
  // GISMO_ASSERT(m_tmp[ID].size()==m_steps,"Too early to finish; I only have "<<m_tmp[ID].size()<<" points for ID "<<ID<<" but I need "<<m_steps<<" points");
  m_tmp.erase(ID);
}

template <class T, class solution_t >
void gsAPALMData<T,solution_t>::printActiveJobs()
{
  // GISMO_ASSERT(m_tmp[ID].size()==m_steps,"Too early to finish; I only have "<<m_tmp[ID].size()<<" points for ID "<<ID<<" but I need "<<m_steps<<" points");
  for (typename std::map<index_t,std::vector<T>>::const_iterator it = m_tmp.cbegin(); it!= m_tmp.cend(); it++)
  {
    std::string tmp = it->second.size()==0 ? "unfinished " : "finished ";
    gsInfo<<"ID = "<<it->first<<": "<<tmp<<"\n";
  }
}

template <class T, class solution_t >
void gsAPALMData<T,solution_t>::removeJob(index_t ID)
{
  // Remove ID
  if (m_verbose==2) gsInfo<<"Erasing finished job on level = "<<m_jobs[ID].first<<"; time = "<<m_jobs[ID].second<<"\n";
  m_jobs.erase(ID);
  if (m_verbose==2)
  {
    gsDebug<<"Current active jobs\n";
    for (typename std::map<index_t,std::pair<T,T>>::const_iterator it = m_jobs.begin(); it!=m_jobs.end(); ++it)
      gsDebug<<"ID = "<<it->first<<"; level = "<<it->second.first<<"; time = "<<it->second.second<<"\n";
  }
}

template <class T, class solution_t >
bool gsAPALMData<T,solution_t>::empty()
{
  return m_queue.empty();
}

template <class T, class solution_t >
std::pair<std::vector<T>,std::vector<solution_t>> gsAPALMData<T,solution_t>::getFlatSolution()
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

template <class T, class solution_t >
void gsAPALMData<T,solution_t>::print()
{
  for (typename std::map<T,solution_t>::const_iterator it = m_solutions.begin(); it!=m_solutions.end(); ++it)
  {
    gsInfo<<"\t";
    gsInfo<<"time = "<<it->first<<"\n";
  }
}

template <class T, class solution_t >
void gsAPALMData<T,solution_t>::printQueue()
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

template <class T, class solution_t >
void gsAPALMData<T,solution_t>::printKnots()
{
  gsInfo<<"xi = \n"<<m_xi.asMatrix()<<"\n";
  gsInfo<<"t = \n"<<m_t.asMatrix()<<"\n";
}

template <class T, class solution_t >
std::pair<index_t,T> gsAPALMData<T,solution_t>::_activeJob(index_t ID)
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