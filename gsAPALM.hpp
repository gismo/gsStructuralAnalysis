/** @file gsAPALM.hpp

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

template <class T>
gsAPALM<T>::gsAPALM(       gsALMBase<T> * ALM,
                                        const gsAPALMData<T,solution_t> & Data)
:
m_ALM(ALM),
m_dataEmpty(Data)
{
  GISMO_ASSERT(m_dataEmpty.empty(),"gsAPALMData must be empty; it will be used to define the options!");
  this->_defaultOptions();
  this->_getOptions();

#ifdef GISMO_WITH_MPI
  gsInfo << "Gismo was compiled with MPI support.\n";

  // Initialize the MPI environment
  // Get the world communicator
  const gsMpi & mpi = gsMpi::init();
  m_comm = mpi.worldComm();
  // MPI_Request req;

  //Get size and rank of the processor
  m_proc_count = m_comm.size();
  m_rank = m_comm.rank();
#else
  m_proc_count = 1;
  m_rank = 0;
#endif
}

template <class T>
void gsAPALM<T>::_defaultOptions()
{
  m_options.addInt("MaxIt","Maximum number of steps",1000);
  m_options.addSwitch("Verbose","Verbosity",false);
  m_options.addInt("SubIntervals","Number of subintervals",2);
  m_options.addSwitch("SingularPoint","Enable singular point detection",false);
  m_options.addReal("BranchLengthMultiplier","Multiplier for the length of branches after bifurcation",1);
  m_options.addReal("BifLengthMultiplier","Multiplier for the length of the first interval after detection of a bifurcation point",0.1);
}

template <class T>
void gsAPALM<T>::_getOptions()
{
  m_maxIterations = m_options.getInt("MaxIt");
  m_verbose = m_options.getSwitch("Verbose");
  m_subIntervals = m_options.getInt("SubIntervals");
  m_singularPoint = m_options.getSwitch("SingularPoint");
  m_branchLengthMult = m_options.getReal("BranchLengthMultiplier");
  m_bifLengthMult = m_options.getReal("BifLengthMultiplier");
}

template <class T>
void gsAPALM<T>::initialize()
{
  this->_getOptions();

  // Initialize solutions
  if (m_starts.size()==0)
  {
    T L0;
    gsVector<T> U0(m_ALM->numDofs());
    U0.setZero();
    L0 = 0.0;
    this->_initStart(U0,L0,m_ALM->getLength());
  }
}

template <class T>
void gsAPALM<T>::_initStart(const gsVector<T> & Ustart, const T & Lstart, const T & dL)
{
  // Initialize solutions
  std::vector<gsVector<T>>  Ustarts; Ustarts.push_back(Ustart);
  std::vector<T>            Lstarts; Lstarts.push_back(Lstart);
  std::vector<T>            dLs;     dLs    .push_back(dL);
  this->_initStart(Ustarts,Lstarts,dLs);
}

template <class T>
void gsAPALM<T>::_initStart(const std::vector<gsVector<T>> & Ustarts, const std::vector<T> & Lstarts, const std::vector<T> & dLs)
{
  GISMO_ASSERT(Ustarts.size()==Lstarts.size(),"Size mismatch in start solutions");
  GISMO_ASSERT(Ustarts.size()==dLs.size(),"Size mismatch in start solutions");

  for (size_t k=0; k!=Ustarts.size(); k++)
    m_starts.push(std::make_tuple(std::make_pair(Ustarts[k],Lstarts[k]),dLs[k],false));
}

template <class T>
void gsAPALM<T>::serialSolve(index_t Nsteps)
{
#ifdef GISMO_WITH_MPI
  if (m_rank!=0)
    m_comm.barrier();
  else
  {
#endif
  GISMO_ASSERT(m_starts.size()>0,"No start point is created. Call initialize first?");

  T Lold, L0;
  gsMatrix<T> Uold, U0;
  T dL, dL0;
  bool bisected, finished, diverged;
  std::vector<solution_t>   solutions;
  std::vector<T>            times;
  std::vector<index_t>      levels;
  while (!m_starts.empty())
  {
    finished = false;
    diverged = false;
    solutions.clear();
    times.clear();
    levels.clear();
    // Initialize solutions
    U0 = Uold = std::get<0>(m_starts.front()).first;
    L0 = Lold = std::get<0>(m_starts.front()).second;
    dL = dL0  = std::get<1>(m_starts.front());
    bisected  = std::get<2>(m_starts.front());
    m_starts.pop();
    T s = 0;

    // If start point comes after a bisection, then multiply it by the corresponding multiplier.
    if (bisected)
      dL *= m_bifLengthMult;
    else
    {
      // Store initial solution
      solutions.push_back({U0,L0});
      times.push_back(s);
      levels.push_back(0);
    }

    m_ALM->setLength(dL);
    m_ALM->setSolution(U0,L0);
    m_ALM->setIndicator(0.0); // RESET INDICATOR
    index_t k=1;

    while (k<Nsteps && !finished)
    {
      if (m_verbose) gsMPIInfo(m_rank)<<"Load step "<< k<<"\t"<<"dL = "<<m_ALM->getLength()<<"; curve time = "<<s<<"\n";

      // Set a step
      m_ALM->step();
      // If not converged, bisect the arc-length
      if ((diverged = !(m_ALM->converged())))
      {
        if (m_verbose) gsMPIInfo(m_rank)<<"Error: Loop terminated, arc length method did not converge.\n";
        dL = m_ALM->reduceLength();
        m_ALM->setSolution(Uold,Lold);
        // bisected = true;
        continue;
      }

      if (m_singularPoint)
      {
        m_ALM->computeStability(m_ALM->solutionU(),m_ALM->options().getSwitch("Quasi"));
        if (m_ALM->stabilityChange())
        {
          gsMPIInfo(m_rank)<<"Bifurcation spotted!"<<"\n";
          m_ALM->computeSingularPoint(1e-4, 5, Uold, Lold, 1e-10, 0, false);
          finished = true;

          gsVector<T> DeltaU = m_ALM->solutionU()-Uold;
          T           DeltaL = m_ALM->solutionL()-Lold;

          dL = m_ALM->distance(DeltaU,DeltaL);
        }
      }

      s += dL;

      T lambda = m_ALM->solutionL();
      std::pair<gsVector<T>,T> pair = std::make_pair(m_ALM->solutionU(),lambda);
      solutions.push_back(pair);
      times.push_back(s);
      levels.push_back(0);
      this->serialStepOutput(pair,s,k);

      Uold = m_ALM->solutionU();
      Lold = m_ALM->solutionL();

      if (finished)
      {
        m_ALM->switchBranch();
        m_starts.push(std::make_tuple(std::make_pair(m_ALM->solutionU(),m_ALM->solutionL()),dL0*m_branchLengthMult,true));
      }

      if (!diverged)
      {
        if (k==1)
        {
          m_ALM->setLength(dL0);
          dL = dL0; // set the length to the length that we want after the first iteration
        }
        else
          dL = m_ALM->resetLength();
      }

      diverged = false;
      k++;
    } // end of steps

    gsAPALMData<T,solution_t> data = m_dataEmpty;
    data.setData(times,solutions);
    data.init();
    m_data.add(data);

    m_solutions.push_back(solutions);
    m_times.push_back(times);
    m_levels.push_back(levels);

  } // end of start points

#ifdef GISMO_WITH_MPI
    m_comm.barrier();
  }
#endif
}


// NOTE: This does not make new branches!
template <class T>
void gsAPALM<T>::parallelSolve()
{
#ifdef GISMO_WITH_MPI
  if (m_comm.size()==1)
    this->parallelSolve_impl<false>();
  else
    this->parallelSolve_impl<true>();
#else
    this->parallelSolve_impl<false>();
#endif
}

template <class T>
template <bool _hasWorkers>
typename std::enable_if< _hasWorkers, void>::type
gsAPALM<T>::parallelSolve_impl()
{
#ifdef GISMO_WITH_MPI
  GISMO_ASSERT(m_comm.size()>1,"Something went wrong. This implementation only works when >1 processes are availbale, but nprocesses = "<<m_comm.size());

  bool stop = false; // stop signal
  solution_t reference;
  index_t ID;
  index_t it = 0;
  std::vector<T> distances;
  std::vector<solution_t> stepSolutions;
  T lowerDistance, upperDistance;

  std::tuple<index_t, T     , solution_t, solution_t, solution_t> dataEntry;

  //----------------------------------------------------------------------------------------
  //MAIN------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------
  if (m_rank==0)
  {
    index_t branch;
    index_t njobs = 0;

    gsMPIInfo(m_rank)<<"Adding workers...\n";
    for (index_t w = 1; w!=m_proc_count; w++)
      m_workers.push(w);

    while (!m_data.empty() && it < m_maxIterations && !m_workers.empty())
    {
      gsMPIInfo(m_rank)<<"Available workers:\t";
      std::queue<index_t> workers_copy = m_workers;
      while (!workers_copy.empty())
      {
          gsMPIInfo(m_rank)<<workers_copy.front()<<"\t";
          workers_copy.pop();
      }
      gsMPIInfo(m_rank)<<"\n";


      branch = m_data.getFirstNonEmptyBranch();
      gsMPIInfo(m_rank)<<"There are "<<m_data.branch(branch).nActive()<<" active jobs and "<<m_data.branch(branch).nWaiting()<<" jobs in the queue of branch "<<branch<<"\n";

      dataEntry = m_data.branch(branch).pop();
      ID = std::get<0>(dataEntry);
      bool success = m_data.branch(branch).getReferenceByID(ID,reference);
      GISMO_ASSERT(success,"Reference not found");

      std::pair<T,T> dataInterval = m_data.branch(branch).jobTimes(ID);
      T tstart = 0;
      T tend = 0;
      std::tie(tstart,tend) = dataInterval;

      this->_sendMainToWorker(m_workers.front(),stop);
      this->_sendMainToWorker(m_workers.front(),
                              branch,
                              dataEntry,
                              m_data.branch(branch).jobTimes(ID),
                              m_data.branch(branch).jobLevel(ID),
                              reference
                              );

      m_workers.pop();
      it++;
      njobs++;
    }

    while (njobs > 0)
    {
      gsMPIInfo(m_rank)<<njobs<<" job(s) running\n";

      MPI_Status status;

      index_t source = MPI_ANY_SOURCE;
      this->_recvWorkerToMain(source,
                              branch,
                              ID,
                              distances,
                              stepSolutions,
                              upperDistance,
                              lowerDistance,
                              &status);

      m_data.branch(branch).submit(ID,distances,stepSolutions,upperDistance,lowerDistance);
      m_data.branch(branch).finishJob(ID);

      // Remove job
      njobs--;
      // Add worker to pool
      m_workers.push(source);

      gsMPIInfo(m_rank)<<"Available workers:\t";
      std::queue<index_t> workers_copy = m_workers;
      while (!workers_copy.empty())
      {
          gsInfo<<workers_copy.front()<<"\t";
          workers_copy.pop();
      }
      gsInfo<<"\n";

      // As long as the queue is not empty, the iterations are lower than the max number and the pool of workers is not empty
      while (!m_data.empty() && it < m_maxIterations && !m_workers.empty())
      {
        // SAME WHILE LOOP AS ABOVE!!!!!!!!
        branch = m_data.getFirstNonEmptyBranch();
        gsMPIInfo(m_rank)<<"There are "<<m_data.branch(branch).nActive()<<" active jobs and "<<m_data.branch(branch).nWaiting()<<" jobs in the queue of branch "<<branch<<"\n";

        dataEntry = m_data.branch(branch).pop();
        ID = std::get<0>(dataEntry);
        bool success = m_data.branch(branch).getReferenceByID(ID,reference);
        GISMO_ASSERT(success,"Reference not found");
        this->_sendMainToWorker(m_workers.front(),stop);
        this->_sendMainToWorker(m_workers.front(),
                                branch,
                                dataEntry,
                                m_data.branch(branch).jobTimes(ID),
                                m_data.branch(branch).jobLevel(ID),
                                reference
                                );
        m_workers.pop();

        it++;
        njobs++;
      }
    }
    stop = true;
    this->_sendMainToAll(stop);

    this->_parallelSolve_finalize();
  }
  else
  {
    while (!stop) // until loop breaks due to stop signal
    {
      this->_recvMainToWorker(0,stop);
      if (stop)
        break;

      std::pair<T,T> dataInterval;
      index_t dataLevel, branch;

      this->_recvMainToWorker(0,
                              branch,
                              dataEntry,
                              dataInterval,
                              dataLevel,
                              reference
                              );

      this->_parallelSolve_worker(dataEntry,
                                  dataInterval,
                                  dataLevel,
                                  reference,
                                  distances,
                                  stepSolutions,
                                  upperDistance,
                                  lowerDistance);

      this->_sendWorkerToMain(0,
                              branch,
                              std::get<0>(dataEntry),
                              distances,
                              stepSolutions,
                              upperDistance,
                              lowerDistance);
    }
  }

#else
  GISMO_NO_IMPLEMENTATION;
#endif

    // #ifdef GISMO_WITH_MPI
    //   gsMPIInfo(m_rank)<<"[MPI Process "<<m_rank<<" of "<<m_comm.size()<<"] Hi!\n";
    // #endif
}

template <class T>
template <bool _hasWorkers>
typename std::enable_if<!_hasWorkers, void>::type
gsAPALM<T>::parallelSolve_impl()
{
  solution_t reference;
  index_t ID;
  index_t it = 0;
  std::vector<T> distances;
  std::vector<solution_t> stepSolutions;
  T lowerDistance, upperDistance;
  index_t branch;

  m_data.print();
  std::tuple<index_t, T     , solution_t, solution_t, solution_t> dataEntry;
  while (!m_data.empty() && it < m_maxIterations)
  {
    branch = m_data.getFirstNonEmptyBranch();
    gsMPIInfo(m_rank)<<"There are "<<m_data.branch(branch).nActive()<<" active jobs and "<<m_data.branch(branch).nWaiting()<<" jobs in the queue of branch "<<branch<<"\n";

    dataEntry = m_data.branch(branch).pop();
    ID = std::get<0>(dataEntry);
    bool success = m_data.branch(branch).getReferenceByID(ID,reference);
    GISMO_ASSERT(success,"Reference not found");
    this->_parallelSolve_worker(dataEntry,
                                m_data.branch(branch).jobTimes(ID),
                                m_data.branch(branch).jobLevel(ID),
                                reference,
                                distances,
                                stepSolutions,
                                upperDistance,
                                lowerDistance
                                );

    m_data.branch(branch).submit(ID,distances,stepSolutions,upperDistance,lowerDistance);
    m_data.branch(branch).finishJob(ID);

    it++;
  }
  this->_parallelSolve_finalize();
}

template <class T>
void gsAPALM<T>::_parallelSolve_finalize()
{
  std::vector<solution_t>   solutions;
  std::vector<T>            times;
  std::vector<index_t>      levels;
  index_t nBranches = m_data.nBranches();
  m_solutions.resize(nBranches);
  m_times.resize(nBranches);
  m_levels.resize(nBranches);
  m_lvlSolutions.resize(nBranches);
  m_lvlTimes.resize(nBranches);
  for (index_t b=0; b!=nBranches; b++)
  {
    std::tie(times,solutions,levels) = m_data.branch(b).getFlatSolution();
    m_times.at(b)     = times;
    m_solutions.at(b) = solutions;
    m_levels.at(b)    = levels;

    index_t maxLevel = m_data.branch(b).maxLevel();
    m_lvlSolutions[b].resize(maxLevel+1);
    m_lvlTimes[b].resize(maxLevel+1);
    for (size_t k=0; k!=m_solutions.size(); k++)
    {
      GISMO_ASSERT(m_lvlSolutions[b].size() > (size_t)m_levels[b][k],"level mismatch, maxLevel = " << maxLevel << "level = "<<m_levels[b][k]);
      GISMO_ASSERT(m_lvlTimes[b].size() > (size_t)m_levels[b][k],"level mismatch, maxLevel = " << maxLevel << "level = "<<m_levels[b][k]);
      m_lvlSolutions[b][m_levels[b][k]].push_back(&(m_solutions[b][k]));
      m_lvlTimes[b][m_levels[b][k]].push_back(&(m_times[b][k]));
    }
  }
}

// NOTE: This does not make new branches!
template <class T>
void gsAPALM<T>::_parallelSolve_worker( const std::tuple<index_t, T     , solution_t, solution_t, solution_t> & dataEntry,
                                        const std::pair<T,T> &  dataInterval,
                                        const index_t &         dataLevel,
                                        const solution_t &      dataReference,
                                        std::vector<T> &        distances,
                                        std::vector<solution_t>&stepSolutions,
                                        T &                     upperDistance,
                                        T &                     lowerDistance )
{
  solution_t start, prev, next, reference;
  index_t ID;
  T tstart = 0;
  T tend = 0;
  T dL, dL0;
  gsVector<T> DeltaU;
  T DeltaL;
  index_t Nintervals = m_subIntervals;
  bool bisected = false;
  T dL_rem = 0;

  T Lprev, Lnext, Lold;
  gsMatrix<T> Uprev, Unext,Uold;

  /// Worker
  stepSolutions.resize(Nintervals);
  std::vector<T> stepTimes(Nintervals);
  distances.resize(Nintervals+1);

  std::tie(ID,dL0,start,prev,next) = dataEntry;
  std::tie(Uold,Lold) = start;
  std::tie(Uprev,Lprev) = prev;
  std::tie(Unext,Lnext) = next;
  std::tie(tstart,tend) = dataInterval;

  gsMPIDebug(m_rank)<<"tstart = "<<tstart<<" , "<<"tend = "<<tend<<"\n";

  gsMatrix<T> Uori = Uold;
  T Lori = Lold;

  dL0 = dL0 / Nintervals;
  dL  = dL0;

  m_ALM->setLength(dL);
  m_ALM->setSolution(Uold,Lold);
  // m_ALM->resetStep();
  m_ALM->setPrevious(Uprev,Lprev);
  gsMPIDebug(m_rank)<<"Start - ||u|| = "<<Uold.norm()<<", L = "<<Lold<<"\n";
  gsMPIDebug(m_rank)<<"Prev  - ||u|| = "<<Uprev.norm()<<", L = "<<Lprev<<"\n";
  gsMPIDebug(m_rank)<<"Next  - ||u|| = "<<Unext.norm()<<", L = "<<Lnext<<"\n";

  gsVector<T> tmpU = Uold-Uprev;

  T s = 0;
  T time = tstart;

  gsMPIInfo(m_rank)<<"Starting with ID "<<ID<<" from (|U|,L) = ("<<Uold.norm()<<","<<Lold<<"), curve time = "<<time<<"\n";
  for (index_t k = 0; k!=Nintervals; k++)
  {
    gsMPIDebug(m_rank)<<"Interval "<<k+1<<" of "<<Nintervals<<"\n";
    gsMPIDebug(m_rank)<<"Start - ||u|| = "<<Uold.norm()<<", L = "<<Lold<<"\n";
    m_ALM->step();
    if (!(m_ALM->converged()))
    {
      gsMPIInfo(m_rank)<<"Error: Loop terminated, arc length method did not converge.\n";
      dL = dL / 2.;
      dL_rem += dL; // add the remainder of the interval to dL_rem
      m_ALM->setLength(dL);
      m_ALM->setSolution(Uold,Lold);
      bisected = true;
      k -= 1;
      continue;
    }
    GISMO_ENSURE(m_ALM->converged(),"Loop terminated, arc length method did not converge.\n");

    std::pair<gsVector<T>,T> pair = std::make_pair(m_ALM->solutionU(),m_ALM->solutionL());
    stepSolutions.at(k) = pair;
    DeltaU = m_ALM->solutionU() - Uold;
    DeltaL = m_ALM->solutionL() - Lold;

    distances.at(k) = m_ALM->distance(DeltaU,DeltaL);

    s += distances.at(k);
    time += distances.at(k);
    stepTimes.at(k) = time;

    Uold = m_ALM->solutionU();
    Lold = m_ALM->solutionL();

    if (!bisected) // if bisected = false
      dL = dL0;
    else
    {
      // if the current interval has finished, but was refined before.
      // The next interval should have the remaining length.
      // Also, Nintervals should increase
      //
      dL = dL_rem;
      Nintervals++;
      stepSolutions.resize(Nintervals);
      stepTimes.resize(Nintervals);
      distances.resize(Nintervals+1);
    }

    this->parallelStepOutput(pair,time,k);

    m_ALM->setLength(dL);
    dL_rem = 0;
    bisected = false;
  }

  gsMPIDebug(m_rank)<<"Ref   - ||u|| = "<<dataReference.first.norm()<<", L = "<<dataReference.second<<"\n";

  GISMO_ASSERT((dataReference.first.normalized()).dot(m_ALM->solutionU().normalized())>-0.8,
                "Reference is almost in the opposite direction of the solution. Are branches mirrored? result:" << (reference.first.normalized()).dot(m_ALM->solutionU().normalized()));

  DeltaU = dataReference.first - m_ALM->solutionU();
  DeltaL = dataReference.second - m_ALM->solutionL();
  distances.back() = m_ALM->distance(DeltaU,DeltaL);

  DeltaU = dataReference.first - Uori;
  DeltaL = dataReference.second - Lori;
  upperDistance = m_ALM->distance(DeltaU,DeltaL);

  DeltaU = stepSolutions.back().first - Uori;
  DeltaL = stepSolutions.back().second - Lori;
  lowerDistance = m_ALM->distance(DeltaU,DeltaL);

  std::vector<solution_t> stepSolutionsExport(Nintervals+2);
  std::vector<T> stepTimesExport(Nintervals+2);
  // Export parallel interval output
  // Solutions
  std::pair<gsVector<T>,T> front = std::make_pair(start.first,start.second);
  stepSolutionsExport.front() = front;
  std::pair<gsVector<T>,T> back  = std::make_pair(dataReference.first,dataReference.second);
  stepSolutionsExport.back() = back;
  // Times
  stepTimesExport.front() = tstart;
  stepTimesExport.back() = tend;

  // Temporarily swap the interval solutions and times
  for (index_t k=0; k!=Nintervals; k++)
  {
    std::swap(stepSolutions.at(k),stepSolutionsExport.at(k+1));
    std::swap(stepTimes.at(k),stepTimesExport.at(k+1));
  }

  this->parallelIntervalOutput(stepSolutionsExport,stepTimesExport,dataLevel,ID);

  // And swap them back
  for (index_t k=0; k!=Nintervals; k++)
  {
    std::swap(stepSolutionsExport.at(k+1),stepSolutions.at(k));
    std::swap(stepTimesExport.at(k+1),stepTimes.at(k));
  }
}


// -----------------------------------------------------------------------------------------------------
// MPI functions
// -----------------------------------------------------------------------------------------------------
#ifdef GISMO_WITH_MPI

template <class T>
void gsAPALM<T>::_sendMainToWorker( const index_t &         workerID,
                                    const index_t &         branch,
                                    const std::tuple<index_t, T     , solution_t, solution_t, solution_t> & dataEntry,
                                    const std::pair<T,T> &  dataInterval,
                                    const index_t &         dataLevel,
                                    const solution_t &      dataReference )
{
  gsMPIInfo(m_rank)<<"Sending data from "<<m_rank<<" to "<<workerID<<"\n";

  index_t     ID;
  T           dL0;
  solution_t  start, prev, next;
  gsVector<T> sol;
  T           load;

  T           tstart, tend;
  index_t     level;

  index_t     vectorSize;

  index_t tag = 0;

  std::tie(ID,dL0,start,prev,next) = dataEntry;
  std::tie(tstart,tend) = dataInterval;
  level = dataLevel;

  MPI_Request req;
  m_comm.isend(&ID        ,1,workerID,&req,tag++);
  m_comm.isend(&branch    ,1,workerID,&req,tag++);
  m_comm.isend(&dL0       ,1,workerID,&req,tag++);
  m_comm.isend(&tstart    ,1,workerID,&req,tag++);
  m_comm.isend(&tend      ,1,workerID,&req,tag++);
  m_comm.isend(&dataLevel ,1,workerID,&req,tag++);

  std::tie(sol,load) = start;
  vectorSize = sol.size();
  m_comm.isend(&vectorSize,1,workerID,        &req,tag++);
  m_comm.isend(sol.data(),sol.size(),workerID,&req,tag++);
  m_comm.isend(&load,1,workerID,              &req,tag++);

  std::tie(sol,load) = prev;
  vectorSize = sol.size();
  m_comm.isend(&vectorSize,1,workerID,        &req,tag++);
  m_comm.isend(sol.data(),sol.size(),workerID,&req,tag++);
  m_comm.isend(&load,1,workerID,              &req,tag++);

  std::tie(sol,load) = next;
  vectorSize = sol.size();
  m_comm.isend(&vectorSize,1,workerID,        &req,tag++);
  m_comm.isend(sol.data(),sol.size(),workerID,&req,tag++);
  m_comm.isend(&load,1,workerID,              &req,tag++);

  std::tie(sol,load) = dataReference;
  vectorSize = sol.size();
  m_comm.isend(&vectorSize,1,workerID,        &req,tag++);
  m_comm.isend(sol.data(),sol.size(),workerID,&req,tag++);
  m_comm.isend(&load,1,workerID,              &req,tag++);
}

template <class T>
void gsAPALM<T>::_sendMainToWorker( const index_t &   workerID,
                                    const bool &      stop )
{
  // gsMPIInfo(m_rank)<<"Sending stop signal from "<<m_rank<<" to "<<workerID<<"\n";

  MPI_Request req;
  m_comm.isend(&stop,1,workerID,&req,0);
}

template <class T>
void gsAPALM<T>::_sendMainToAll(  const bool &      stop )
{
  for (index_t w = 1; w!=m_proc_count; w++)
    this->_sendMainToWorker(w,stop);
}

template <class T>
void gsAPALM<T>::_recvMainToWorker( const index_t &   sourceID,
                                          index_t &         branch,
                                          std::tuple<index_t, T     , solution_t, solution_t, solution_t> & dataEntry,
                                          std::pair<T,T> &  dataInterval,
                                          index_t &         dataLevel,
                                          solution_t &      dataReference,
                                          MPI_Status *      status)
{
  gsMPIInfo(m_rank)<<"Receiving data on "<<m_rank<<" from "<<sourceID<<"\n";

  index_t     ID;
  T           dL0;
  solution_t  start, prev, next;
  gsVector<T> sol;
  T           load;

  T           tstart, tend;
  index_t     level;

  index_t     vectorSize;

  index_t tag = 0;

  level = dataLevel;

  m_comm.recv(&ID        ,1,sourceID,tag++, status);
  m_comm.recv(&branch    ,1,sourceID,tag++, status);
  m_comm.recv(&dL0       ,1,sourceID,tag++, status);

  m_comm.recv(&tstart    ,1,sourceID,tag++, status);
  m_comm.recv(&tend      ,1,sourceID,tag++, status);

  dataInterval = std::make_pair(tstart,tend);
  m_comm.recv(&dataLevel ,1,sourceID,tag++, status);


  m_comm.recv(&vectorSize,1,        sourceID,tag++, status);
  sol.resize(vectorSize);
  m_comm.recv(sol.data(),sol.size(),sourceID,tag++, status);
  m_comm.recv(&load,1,              sourceID,tag++, status);
  start = std::make_pair(sol,load);

  m_comm.recv(&vectorSize,1,        sourceID,tag++, status);
  sol.resize(vectorSize);
  m_comm.recv(sol.data(),sol.size(),sourceID,tag++, status);
  m_comm.recv(&load,1,              sourceID,tag++, status);
  prev = std::make_pair(sol,load);

  m_comm.recv(&vectorSize,1,        sourceID,tag++, status);
  sol.resize(vectorSize);
  m_comm.recv(sol.data(),sol.size(),sourceID,tag++, status);
  m_comm.recv(&load,1,              sourceID,tag++, status);
  next = std::make_pair(sol,load);

  m_comm.recv(&vectorSize,1,        sourceID,tag++, status);
  sol.resize(vectorSize);
  m_comm.recv(sol.data(),sol.size(),sourceID,tag++, status);
  m_comm.recv(&load,1,              sourceID,tag++, status);
  dataReference = std::make_pair(sol,load);

  dataEntry = std::make_tuple(ID,dL0,start,prev,next);
}

template <class T>
void gsAPALM<T>::_recvMainToWorker(  const  index_t &   sourceID,
                                            bool &      stop,
                                            MPI_Status *status   )
{
  // gsMPIInfo(m_rank)<<"Receiving stop signal on "<<m_rank<<" from "<<sourceID<<"\n";

  m_comm.recv(&stop,1,sourceID,0,status);
}

template <class T>
void gsAPALM<T>::_sendWorkerToMain( const index_t &                   mainID,
                                    const index_t &                   branch,
                                    const index_t &                   jobID,
                                    const std::vector<T> &            distances,
                                    const std::vector<solution_t> &   stepSolutions,
                                    const T &                         upperDistance,
                                    const T &                         lowerDistance )
{
  gsMPIInfo(m_rank)<<"Sending data from "<<m_rank<<" to "<<mainID<<"\n";

  index_t tag = 0;
  index_t size = stepSolutions.size();
  index_t vectorSize;

  MPI_Request req;
  m_comm.isend(&branch        ,1,mainID,&req,tag++);
  m_comm.isend(&jobID         ,1,mainID,&req,tag++);
  m_comm.isend(&upperDistance ,1,mainID,&req,tag++);
  m_comm.isend(&lowerDistance ,1,mainID,&req,tag++);
  m_comm.isend(&size          ,1,mainID,&req,tag++);

  for (index_t k=0; k!=size; k++)
  {
    m_comm.isend(&(distances.at(k))              ,1            ,mainID,&req,tag++);
    vectorSize = stepSolutions.at(k).first.size();
    m_comm.isend(&vectorSize                     ,1            ,mainID,&req,tag++);
    m_comm.isend(stepSolutions.at(k).first.data(),vectorSize   ,mainID,&req,tag++);
    m_comm.isend(&(stepSolutions.at(k).second)   ,1            ,mainID,&req,tag++);
  }
  m_comm.isend(&(distances.at(size))             ,1            ,mainID,&req,tag++);
}

template <class T>
void gsAPALM<T>::_recvWorkerToMain( index_t &                   sourceID,
                                    index_t &                   branch,
                                    index_t &                   jobID,
                                    std::vector<T> &            distances,
                                    std::vector<solution_t> &   stepSolutions,
                                    T &                         upperDistance,
                                    T &                         lowerDistance,
                                    MPI_Status *                status)
{

  index_t tag = 0;
  index_t size;
  index_t vectorSize;

  m_comm.recv(&branch         ,1,sourceID,tag++, status);
  gsMPIInfo(m_rank)<<"Receiving data on "<<m_rank<<" from "<<sourceID<<"\n";
  if (sourceID==MPI_ANY_SOURCE)  sourceID = status->MPI_SOURCE; // Overwrite source ID to be unique
  m_comm.recv(&jobID          ,1,sourceID,tag++, status);
  m_comm.recv(&upperDistance  ,1,sourceID,tag++, status);
  m_comm.recv(&lowerDistance  ,1,sourceID,tag++, status);
  m_comm.recv(&size           ,1,sourceID,tag++, status);

  stepSolutions.resize(size);
  distances.resize(size+1);

  gsVector<T> sol;
  T load;
  for (index_t k=0; k!=size; k++)
  {
    m_comm.recv(&(distances.at(k))  ,1            ,sourceID,tag++, status);
    m_comm.recv(&vectorSize         ,1            ,sourceID,tag++, status);
    sol.resize(vectorSize);
    m_comm.recv(sol.data()          ,vectorSize   ,sourceID,tag++, status);
    m_comm.recv(&load               ,1            ,sourceID,tag++, status);

    stepSolutions.at(k) = std::make_pair(sol,load);
  }
  m_comm.recv(&(distances.at(size)) ,1            ,sourceID,tag++, status);
}

#endif

// -----------------------------------------------------------------------------------------------------
} // namespace gismo