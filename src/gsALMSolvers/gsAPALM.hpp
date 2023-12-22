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
#ifdef GISMO_WITH_MPI
gsAPALM<T>::gsAPALM(  gsALMBase<T> * ALM,
                      const gsAPALMData<T,solution_t> & Data,
                      const gsMpiComm & comm                 )
#else
gsAPALM<T>::gsAPALM(  gsALMBase<T> * ALM,
                      const gsAPALMData<T,solution_t> & Data,
                      const gsSerialComm & comm              )
#endif
:
m_ALM(ALM),
m_dataEmpty(Data),
m_comm(comm)
{
  GISMO_ASSERT(m_dataEmpty.empty(),"gsAPALMData must be empty; it will be used to define the options!");
  this->_defaultOptions();
  this->_getOptions();

  //Get size and rank of the processor
  m_proc_count = m_comm.size();
  m_rank = m_comm.rank();
}

template <class T>
gsAPALM<T>::gsAPALM(  gsALMBase<T> * ALM,
                      const gsAPALMData<T,solution_t> & Data)
:
m_ALM(ALM),
m_dataEmpty(Data),
m_comm(gsSerialComm())
{
  GISMO_ASSERT(m_dataEmpty.empty(),"gsAPALMData must be empty; it will be used to define the options!");
  this->_defaultOptions();
  this->_getOptions();

  m_proc_count = 1;
  m_rank = 0;
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

  T Lold, Lprev, L0;
  gsMatrix<T> Uold, Uprev, U0;
  T dL, dL0;
  bool bisected, finished;
  std::vector<solution_t>   solutions;
  std::vector<T>            times;
  std::vector<index_t>      levels;
  while (!m_starts.empty())
  {
    finished = false;
    solutions.clear();
    times.clear();
    levels.clear();
    // Initialize solutions
    U0 = Uold = Uprev = std::get<0>(m_starts.front()).first;
    L0 = Lold = Lprev = std::get<0>(m_starts.front()).second;
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

    index_t k=1;
    std::vector<solution_t> tmpSolutions;
    while (k<Nsteps && !finished)
    {
      if (m_verbose) gsMPIInfo(m_rank)<<"Load step "<< k<<"\t"<<"dL = "<<m_ALM->getLength()<<"; curve time = "<<s<<"\n";
      solution_t start = std::make_pair(Uold,Lold);
      solution_t prev = std::make_pair(Uprev,Lprev);

      std::tuple<index_t, T, solution_t, solution_t> dataEntry = std::make_tuple(k,dL,start,prev);

      this->_initiation(dataEntry,s,dL,tmpSolutions,finished);

      // Update curve length
      s += dL;

      // Update solutions
      solutions.push_back(tmpSolutions[0]);
      times.push_back(s);
      levels.push_back(0);

      // Update previous and old solution
      Uprev = Uold;
      Lprev = Lold;
      Uold = tmpSolutions[0].first;
      Lold = tmpSolutions[0].second;

      // Store new branch as start point
      if (finished)
        m_starts.push(std::make_tuple(tmpSolutions[1],dL0*m_branchLengthMult,true));

      dL = dL0; // set the length to the length that we want after the first iteration

      // Update step counter
      k++;
    } // end of steps

    // Store data of the branch
    gsAPALMData<T,solution_t> data = m_dataEmpty; // create new data set
    data.setData(times,solutions); // initialize the dataset with the newly computed data
    data.init(); // initialize the dataset
    m_data.add(data); // add the dataset to the tree

    // Store the solutions, times and levels
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
void gsAPALM<T>::solve(index_t Nsteps)
{
#ifdef GISMO_WITH_MPI
  if (m_comm.size()==1)
    this->_solve_impl<false>(Nsteps);
  else
    this->_solve_impl<true>(Nsteps);
#else
    this->_solve_impl<false>(Nsteps);
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
  index_t branch;

  std::vector<T> distances;
  std::vector<solution_t> stepSolutions;
  T lowerDistance, upperDistance;
  std::pair<T,T> dataInterval;

  index_t dataLevel;

  std::tuple<index_t, T     , solution_t, solution_t> dataEntry;

  //----------------------------------------------------------------------------------------
  //MAIN------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------
  if (m_rank==0)
  {
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
      dataLevel = m_data.branch(branch).jobLevel(ID);
      bool success = m_data.branch(branch).getReferenceByID(ID,reference);
      GISMO_ASSERT(success,"Reference not found");

      this->_sendMainToWorker(m_workers.front(),stop);
      this->_sendMainToWorker(m_workers.front(),
                              branch,
                              ID,
                              dataLevel);
      this->_sendMainToWorker(m_workers.front(),
                              dataEntry,
                              m_data.branch(branch).jobTimes(ID),
                              reference
                              );

      m_workers.pop();
      it++;
      njobs++;
    }

    while (njobs > 0)
    {
      gsMPIInfo(m_rank)<<njobs<<" job(s) running\n";

      index_t source = MPI_ANY_SOURCE;
      this->_recvWorkerToMain(source,
                              branch,
                              ID);
      this->_recvWorkerToMain(source,
                              distances,
                              stepSolutions,
                              upperDistance,
                              lowerDistance);

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
        dataLevel = m_data.branch(branch).jobLevel(ID);
        bool success = m_data.branch(branch).getReferenceByID(ID,reference);
        GISMO_ASSERT(success,"Reference not found");
        this->_sendMainToWorker(m_workers.front(),stop);
        this->_sendMainToWorker(m_workers.front(),
                                branch,
                                ID,
                                dataLevel);
        this->_sendMainToWorker(m_workers.front(),
                                dataEntry,
                                m_data.branch(branch).jobTimes(ID),
                                reference
                                );
        m_workers.pop();

        it++;
        njobs++;
      }
    }
    stop = true;
    this->_sendMainToAll(stop);

    this->_finalize();
  }
  else
  {
    while (!stop) // until loop breaks due to stop signal
    {
      this->_recvMainToWorker(0,stop);
      if (stop)
        break;

      this->_recvMainToWorker(0,
                              branch,
                              ID,
                              dataLevel);
      this->_recvMainToWorker(0,
                              dataEntry,
                              dataInterval,
                              reference
                              );

      this->_correction(dataEntry,
                        dataInterval,
                        dataLevel,
                        reference,
                        distances,
                        stepSolutions,
                        upperDistance,
                        lowerDistance);

      this->_sendWorkerToMain(0,
                              branch,
                              ID);
      this->_sendWorkerToMain(0,
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
  index_t dataLevel;

  m_data.print();
  std::tuple<index_t, T     , solution_t, solution_t> dataEntry;
  while (!m_data.empty() && it < m_maxIterations)
  {
    branch = m_data.getFirstNonEmptyBranch();
    gsMPIInfo(m_rank)<<"There are "<<m_data.branch(branch).nActive()<<" active jobs and "<<m_data.branch(branch).nWaiting()<<" jobs in the queue of branch "<<branch<<"\n";

    dataEntry = m_data.branch(branch).pop();
    ID = std::get<0>(dataEntry);
    dataLevel = m_data.branch(branch).jobLevel(ID);
    bool success = m_data.branch(branch).getReferenceByID(ID,reference);
    GISMO_ASSERT(success,"Reference not found");
    this->_correction(dataEntry,
                      m_data.branch(branch).jobTimes(ID),
                      dataLevel,
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
  this->_finalize();
}

template <class T>
void gsAPALM<T>::_finalize()
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

template <class T>
template <bool _hasWorkers>
typename std::enable_if< _hasWorkers, void>::type
gsAPALM<T>::_solve_impl(index_t Nsteps)
{
#ifdef GISMO_WITH_MPI
  GISMO_ASSERT(m_comm.size()>1,"Something went wrong. This implementation only works when >1 processes are availbale, but nprocesses = "<<m_comm.size());

  bool stop = false; // stop signal
  solution_t reference;
  index_t ID;
  index_t it = 0;
  index_t branch;
  index_t dataLevel;
  std::tuple<index_t, T     , solution_t, solution_t> dataEntry;

  // Correction
  std::pair<T,T> dataInterval;
  std::vector<T> distances;
  std::vector<solution_t> stepSolutions;
  T lowerDistance, upperDistance;

  // Initiation
  T tstart;
  T distance;
  std::vector<solution_t> solutions;
  bool bifurcation;


  //----------------------------------------------------------------------------------------
  //MAIN------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------
  if (m_rank==0)
  {
    index_t njobs = 0;

    // K keeps track of the iterations per branch
    std::vector<index_t> K(1);
    K[0] = 1;

    // Fill the start interval list
    T L0;
    gsMatrix<T> U0;
    T dL;
    while (!m_starts.empty())
    {
      // Initialize solutions
      U0 = std::get<0>(m_starts.front()).first;
      L0 = std::get<0>(m_starts.front()).second;
      dL = std::get<1>(m_starts.front());
      m_starts.pop();

      gsAPALMData<T,solution_t> data = m_dataEmpty;
      branch = m_data.add(data);
      m_data.branch(branch).addStartPoint(T(0),std::make_pair(U0,L0));
      m_data.branch(branch).setLength(dL); // sets the default length that is used when started from a Point
    }

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
      dataLevel = m_data.branch(branch).jobLevel(ID);
      // Initialization intervals start at level 0
      if (dataLevel==0)
      {
        this->_sendMainToWorker(m_workers.front(),stop);
        this->_sendMainToWorker(m_workers.front(),
                                branch,
                                ID,
                                dataLevel);
        this->_sendMainToWorker(m_workers.front(),
                                dataEntry,
                                m_data.branch(branch).jobStartTime(ID)
                                );
      }
      // Correction intervals have level > 0
      else
      {
        bool success = m_data.branch(branch).getReferenceByID(ID,reference);
        GISMO_ASSERT(success,"Reference not found");

        this->_sendMainToWorker(m_workers.front(),stop);
        this->_sendMainToWorker(m_workers.front(),
                                branch,
                                ID,
                                dataLevel);
        this->_sendMainToWorker(m_workers.front(),
                                dataEntry,
                                m_data.branch(branch).jobTimes(ID),
                                reference
                                );
      }
      m_workers.pop();
      it++;
      njobs++;
    }

    while (njobs > 0)
    {
      gsMPIInfo(m_rank)<<njobs<<" job(s) running\n";
      index_t source = MPI_ANY_SOURCE;
      // Receive meta-data
      this->_recvWorkerToMain(source,
                              branch,
                              ID);

      // Initialization intervals start at level 0
      dataLevel = m_data.branch(branch).jobLevel(ID);
      if (dataLevel == 0)
      {
        this->_recvWorkerToMain(source,
                                distance,
                                solutions,
                                bifurcation);

        T tstart = m_data.branch(branch).jobStartTime(ID);
        m_data.branch(branch).appendData(tstart+distance,solutions[0],false);
        m_data.branch(branch).finishJob(ID);

        if (K[branch]++ < Nsteps-1) // add a new point
        {
          if (!bifurcation)
          {
            GISMO_ASSERT(solutions.size()==1,"There must be one solution, but solutions.size() = "<<solutions.size());
            m_data.branch(branch).appendPoint(true);
            // sets the default length that is used when started from a Point.
            // After bifurcation, it's the original one times the branch length multiplier times the bifurcation length multiplier
            if (branch!=0 && ID==0)
              m_data.branch(branch).setLength(m_data.branch(branch).getLength()/m_bifLengthMult);
          }
          else
          {
            GISMO_ASSERT(solutions.size()==2,"There must be two solutions!");
            gsAPALMData<T,solution_t> data = m_dataEmpty;
            branch = m_data.add(data);
            K.push_back(1);
            m_data.branch(branch).addStartPoint(T(0),solutions[1],true);
            // Sets the default length that is used when started from a Point.
            // After bifurcation, it's the original one times the branch length multiplier times the bifurcation length multiplier
            m_data.branch(branch).setLength(m_ALM->getLength()*m_branchLengthMult*m_bifLengthMult);
          }
        }
      }
      // Correction intervals have level > 0
      else
      {
        this->_recvWorkerToMain(source,
                                distances,
                                stepSolutions,
                                upperDistance,
                                lowerDistance);

        m_data.branch(branch).submit(ID,distances,stepSolutions,upperDistance,lowerDistance);
        m_data.branch(branch).finishJob(ID);
      }

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
        // // SAME WHILE LOOP AS ABOVE!!!!!!!!
        branch = m_data.getFirstNonEmptyBranch();
        gsMPIInfo(m_rank)<<"There are "<<m_data.branch(branch).nActive()<<" active jobs and "<<m_data.branch(branch).nWaiting()<<" jobs in the queue of branch "<<branch<<"\n";

        dataEntry = m_data.branch(branch).pop();
        ID = std::get<0>(dataEntry);

        // Initialization intervals start at level 0
        dataLevel = m_data.branch(branch).jobLevel(ID);
        if (dataLevel==0)
        {
          this->_sendMainToWorker(m_workers.front(),stop);
          this->_sendMainToWorker(m_workers.front(),
                                  branch,
                                  ID,
                                  dataLevel);

          this->_sendMainToWorker(m_workers.front(),
                                  dataEntry,
                                  m_data.branch(branch).jobStartTime(ID)
                                  );
        }
        // Correction intervals have level > 0
        else
        {
          bool success = m_data.branch(branch).getReferenceByID(ID,reference);
          GISMO_ASSERT(success,"Reference not found");

          this->_sendMainToWorker(m_workers.front(),stop);
          this->_sendMainToWorker(m_workers.front(),
                                  branch,
                                  ID,
                                  dataLevel);
          this->_sendMainToWorker(m_workers.front(),
                                  dataEntry,
                                  m_data.branch(branch).jobTimes(ID),
                                  reference
                                  );
        }
        m_workers.pop();
        it++;
        njobs++;
      }
    }
    stop = true;
    this->_sendMainToAll(stop);

    this->_finalize();
  }
  else
  {
    while (!stop) // until loop breaks due to stop signal
    {
      this->_recvMainToWorker(0,stop);
      if (stop)
        break;

      // Receive meta-data
      this->_recvMainToWorker(0,
                              branch,
                              ID,
                              dataLevel);

      // Initialization intervals start at level 0
      if (dataLevel==0)
      {
        this->_recvMainToWorker(0,
                                dataEntry,
                                tstart
                                );

        this->_initiation(dataEntry,
                          tstart,
                          distance,
                          solutions,
                          bifurcation
                          );

        this->_sendWorkerToMain(0,
                                branch,
                                ID);
        this->_sendWorkerToMain(0,
                                distance,
                                solutions,
                                bifurcation);
      }
      // Correction intervals have level > 0
      else
      {
        this->_recvMainToWorker(0,
                                dataEntry,
                                dataInterval,
                                reference
                                );

        this->_correction(dataEntry,
                          dataInterval,
                          dataLevel,
                          reference,
                          distances,
                          stepSolutions,
                          upperDistance,
                          lowerDistance);

        this->_sendWorkerToMain(0,
                                branch,
                                ID);
        this->_sendWorkerToMain(0,
                                distances,
                                stepSolutions,
                                upperDistance,
                                lowerDistance);
      }
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
gsAPALM<T>::_solve_impl(index_t Nsteps)
{
  index_t ID;
  index_t it = 0;
  std::vector<index_t> K(1);
  K[0] = 1;

  index_t branch;
  index_t dataLevel;

  T L0;
  gsMatrix<T> U0;
  T dL;
  while (!m_starts.empty())
  {
    // Initialize solutions
    U0 = std::get<0>(m_starts.front()).first;
    L0 = std::get<0>(m_starts.front()).second;
    dL = std::get<1>(m_starts.front());
    m_starts.pop();

    gsAPALMData<T,solution_t> data = m_dataEmpty;
    branch = m_data.add(data);
    m_data.branch(branch).addStartPoint(T(0),std::make_pair(U0,L0));
    m_data.branch(branch).setLength(dL); // sets the default length that is used when started from a Point
  }

  m_data.print();
  std::tuple<index_t, T     , solution_t, solution_t> dataEntry;
  while (!m_data.empty() && it < m_maxIterations)
  {
    branch = m_data.getFirstNonEmptyBranch();
    gsMPIInfo(m_rank)<<"There are "<<m_data.branch(branch).nActive()<<" active jobs and "<<m_data.branch(branch).nWaiting()<<" jobs in the queue of branch "<<branch<<"\n";

    dataEntry = m_data.branch(branch).pop();
    ID = std::get<0>(dataEntry);

    dataLevel = m_data.branch(branch).jobLevel(ID);
    if (dataLevel==0)
    {
      if (m_verbose) gsMPIInfo(m_rank)<<"Initialization\n";
      T startTime = m_data.branch(branch).jobStartTime(ID);
      T distance;
      std::vector<solution_t> solutions;
      bool bifurcation;
      this->_initiation(dataEntry,
                        startTime,
                        distance,
                        solutions,
                        bifurcation
                        );

      m_data.branch(branch).appendData(startTime+distance,solutions[0],false);
      m_data.branch(branch).finishJob(ID);
      // if no bifurcation is hit, we mark the end of the interval as the start point

      if (K[branch]++ < Nsteps - 1) // add a new point
      {
        if (!bifurcation)
        {
          GISMO_ASSERT(solutions.size()==1,"There must be one solution!");
          m_data.branch(branch).appendPoint(true);
          // sets the default length that is used when started from a Point.
          // After bifurcation, it's the original one times the branch length multiplier times the bifurcation length multiplier
          if (branch!=0 && ID==0)
            m_data.branch(branch).setLength(m_data.branch(branch).getLength()/m_bifLengthMult);
        }
        else
        {
          GISMO_ASSERT(solutions.size()==2,"There must be two solutions!");
          gsAPALMData<T,solution_t> data = m_dataEmpty;
          branch = m_data.add(data);
          K.push_back(1);
          m_data.branch(branch).addStartPoint(T(0),solutions[1],true); // use the second solution that is given in solutions
          // Sets the default length that is used when started from a Point.
          // After bifurcation, it's the original one times the branch length multiplier times the bifurcation length multiplier
          m_data.branch(branch).setLength(m_ALM->getLength()*m_branchLengthMult*m_bifLengthMult);
        }
      }
    }
    else
    {
      if (m_verbose) gsMPIInfo(m_rank)<<"Correction\n";
      solution_t reference;
      std::vector<T> distances;
      std::vector<solution_t> stepSolutions;
      T lowerDistance, upperDistance;
      bool success = m_data.branch(branch).getReferenceByID(ID,reference);
      GISMO_ASSERT(success,"Reference not found");
      this->_correction(dataEntry,
                        m_data.branch(branch).jobTimes(ID),
                        dataLevel,
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
  }
  this->_finalize();
}

// NOTE: This does not make new branches!
template <class T>
void gsAPALM<T>::_initiation( const std::tuple<index_t, T     , solution_t, solution_t> & dataEntry,
                              const T &               startTime,
                              T &                     distance,
                              std::vector<solution_t>&solutions,
                              bool &                  bifurcation )
{
  solution_t start, prev;
  index_t ID;
  T tstart = startTime;
  T dL, dL0;
  gsVector<T> DeltaU;
  T DeltaL;

  T Lprev, Lold;
  gsMatrix<T> Uprev, Uold;

  /// Worker
  solutions.clear();
  std::tie(ID,dL,start,prev) = dataEntry;
  std::tie(Uold,Lold) = start;
  std::tie(Uprev,Lprev) = prev;
  dL0 = dL;

  m_ALM->setLength(dL);
  m_ALM->setSolution(Uold,Lold);
  // m_ALM->resetStep();
  m_ALM->setPrevious(Uprev,Lprev);
  gsMPIDebug(m_rank)<<"Start - ||u|| = "<<Uold.norm()<<", L = "<<Lold<<"\n";
  gsMPIDebug(m_rank)<<"Prev  - ||u|| = "<<Uprev.norm()<<", L = "<<Lprev<<"\n";

  bool diverged = true;
  bifurcation = false;
  while (diverged)
  {
    gsMPIInfo(m_rank)<<"Starting with ID "<<ID<<" from (|U|,L) = ("<<Uold.norm()<<","<<Lold<<"), curve time = "<<tstart<<", arc-length = "<<dL<<"\n";
    // Set a step
    gsStatus status = m_ALM->step();
    diverged = (status!=gsStatus::Success);
    if (status==gsStatus::NotConverged || status==gsStatus::AssemblyError)
    {
      if (m_verbose) gsMPIInfo(m_rank)<<"Error: Loop terminated, arc length method did not converge.\n";
      dL = dL / 2;
      m_ALM->setLength(dL);
      m_ALM->setSolution(Uold,Lold);
      continue;
    }
  }
  dL = dL0;
  if (m_singularPoint)
  {
    m_ALM->computeStability(m_ALM->options().getSwitch("Quasi"));
    if (m_ALM->stabilityChange())
    {
      gsMPIInfo(m_rank)<<"Bifurcation spotted!"<<"\n";
      if (m_ALM->isBifurcation(false))
      {
        m_ALM->computeSingularPoint(Uold,Lold,false,false,false);
        bifurcation = true;
      }
    }
  }

  DeltaU = m_ALM->solutionU() - Uold;
  DeltaL = m_ALM->solutionL() - Lold;
  distance = m_ALM->distance(DeltaU,DeltaL);
  solutions.push_back(std::make_pair(m_ALM->solutionU(),m_ALM->solutionL()));

  this->serialStepOutput(solutions[0],tstart,ID);

  if (bifurcation)
  {
    m_ALM->switchBranch();
    // add the extra point after bifurcation to the export
    solutions.push_back(std::make_pair(m_ALM->solutionU(),m_ALM->solutionL()));
  }
}

// NOTE: This does not make new branches!
template <class T>
void gsAPALM<T>::_correction( const std::tuple<index_t, T     , solution_t, solution_t> & dataEntry,
                              const std::pair<T,T> &  dataInterval,
                              const index_t &         dataLevel,
                              const solution_t &      dataReference,
                              std::vector<T> &        distances,
                              std::vector<solution_t>&stepSolutions,
                              T &                     upperDistance,
                              T &                     lowerDistance )
{
  solution_t start, prev, reference;
  index_t ID;
  T tstart = 0;
  T tend = 0;
  T dL, dL0;
  gsVector<T> DeltaU;
  T DeltaL;
  index_t Nintervals = m_subIntervals;
  bool bisected = false;
  T dL_rem = 0;

  T Lprev, Lold;
  gsMatrix<T> Uprev, Uold;

  /// Worker
  stepSolutions.resize(Nintervals);
  std::vector<T> stepTimes(Nintervals);
  distances.resize(Nintervals+1);

  std::tie(ID,dL0,start,prev) = dataEntry;
  std::tie(Uold,Lold) = start;
  std::tie(Uprev,Lprev) = prev;
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

  T s = 0;
  T time = tstart;

  gsMPIInfo(m_rank)<<"Starting with ID "<<ID<<" from (|U|,L) = ("<<Uold.norm()<<","<<Lold<<"), curve time = "<<time<<"\n";
  for (index_t k = 0; k!=Nintervals; k++)
  {
    gsMPIDebug(m_rank)<<"Interval "<<k+1<<" of "<<Nintervals<<"\n";
    gsMPIDebug(m_rank)<<"Start - ||u|| = "<<Uold.norm()<<", L = "<<Lold<<"\n";

    gsStatus status = m_ALM->step();
    if (status==gsStatus::NotConverged || status==gsStatus::AssemblyError)
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

  GISMO_ASSERT(dataReference.first.normalized().size()==m_ALM->solutionU().normalized().size(),"Reference solution and current solution have different size! ref.size() = "<<dataReference.first.normalized().size()<<"; sol.size() = "<<m_ALM->solutionU().normalized().size());
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
void gsAPALM<T>::_sendMainToWorker( const index_t & workerID,
                                    const index_t & branch,
                                    const index_t & jobID,
                                    const index_t & dataLevel)
{
  gsMPIInfo(m_rank)<<"Sending data from "<<m_rank<<" to "<<workerID<<"\n";
  MPI_Request req[3];
  m_comm.isend(&branch        ,1,workerID,&req[0] ,0 );
  m_comm.isend(&jobID         ,1,workerID,&req[1] ,1 );
  m_comm.isend(&dataLevel     ,1,workerID,&req[2] ,2 );
  MPI_Waitall( 3, req, MPI_STATUSES_IGNORE );
}

template <class T>
void gsAPALM<T>::_sendMainToWorker( const index_t &         workerID,
                                    const std::tuple<index_t, T     , solution_t, solution_t> & dataEntry,
                                    const std::pair<T,T> &  dataInterval,
                                    const solution_t &      dataReference )
{
  gsMPIInfo(m_rank)<<"Sending data from "<<m_rank<<" to "<<workerID<<"\n";

  index_t     ID;
  T           dL0;
  solution_t  start,    prev;
  gsVector<T> startU,   prevU,    refU;
  index_t     startSize,prevSize, refSize;
  T           startL,   prevL,    refL;

  T           tstart, tend;

  std::tie(ID,dL0,start,prev) = dataEntry;
  std::tie(tstart,tend) = dataInterval;

  std::tie(startU,startL) = start;
  startSize = startU.size();
  std::tie(prevU,prevL) = prev;
  prevSize  = prevU.size();
  std::tie(refU,refL) = dataReference;
  refSize   = refU.size();

  // Put all data 0-14 in a struct and send a single struct
  MPI_Request req[10];
  m_comm.isend(&ID        ,1,workerID,&req[0] ,0 );
  m_comm.isend(&dL0       ,1,workerID,&req[1] ,1 );
  m_comm.isend(&tstart    ,1,workerID,&req[2] ,2 );
  m_comm.isend(&tend      ,1,workerID,&req[3] ,3 );

  m_comm.isend(&startSize ,1,workerID,&req[4] ,4 );
  m_comm.isend(&startL    ,1,workerID,&req[5] ,5 );

  m_comm.isend(&prevSize  ,1,workerID,&req[6] ,6 );
  m_comm.isend(&prevL     ,1,workerID,&req[7] ,7 );

  m_comm.isend(&refSize   ,1,workerID,&req[8] ,8 );
  m_comm.isend(&refL      ,1,workerID,&req[9] ,9 );

  MPI_Waitall( 10, req     , MPI_STATUSES_IGNORE );
  MPI_Request req_data[3];
  m_comm.isend(startU.data(), startSize,workerID,&req_data[0],10);
  m_comm.isend(prevU.data(),  prevSize, workerID,&req_data[1],11);
  m_comm.isend(refU.data(),   refSize,  workerID,&req_data[2],12);

  MPI_Waitall( 3 , req_data, MPI_STATUSES_IGNORE );
}

template <class T>
void gsAPALM<T>::_sendMainToWorker( const index_t &         workerID,
                                    const std::tuple<index_t, T     , solution_t, solution_t> & dataEntry,
                                    const T       &         startTime)
{
  gsMPIInfo(m_rank)<<"Sending data from "<<m_rank<<" to "<<workerID<<"\n";

  index_t     ID;
  T           dL0;
  solution_t  start,    prev;
  gsVector<T> startU,   prevU;
  index_t     startSize,prevSize;
  T           startL,   prevL;

  T           tstart;

  std::tie(ID,dL0,start,prev) = dataEntry;
  tstart = startTime;

  std::tie(startU,startL) = start;
  startSize = startU.size();
  std::tie(prevU,prevL) = prev;
  prevSize  = prevU.size();

  // Put all data 0-10 in a struct and send a single struct
  MPI_Request req[7];
  m_comm.isend(&ID        ,1,workerID,&req[0] ,0 );
  m_comm.isend(&dL0       ,1,workerID,&req[1] ,1 );
  m_comm.isend(&tstart    ,1,workerID,&req[2] ,2 );

  m_comm.isend(&startSize ,1,workerID,&req[3] ,3 );
  m_comm.isend(&startL    ,1,workerID,&req[4] ,4 );

  m_comm.isend(&prevSize  ,1,workerID,&req[5] ,5 );
  m_comm.isend(&prevL     ,1,workerID,&req[6] ,6 );
  MPI_Waitall( 7, req     , MPI_STATUSES_IGNORE );

  MPI_Request req_data[2];
  m_comm.isend(startU.data(), startSize,workerID,&req_data[0],7);
  m_comm.isend(prevU.data(),  prevSize, workerID,&req_data[1],8);

  MPI_Waitall( 2, req_data, MPI_STATUSES_IGNORE );
}

template <class T>
void gsAPALM<T>::_sendMainToWorker( const index_t &   workerID,
                                    const bool &      stop )
{
  // NOTE: when we handle the stop signal with isend/irecv, we get a deadlock
  // TODO: handle stop with broadcast
  gsMPIInfo(m_rank)<<"Sending stop signal from "<<m_rank<<" to "<<workerID<<"\n";
  m_comm.send(&stop,1,workerID,99);
}

template <class T>
void gsAPALM<T>::_sendMainToAll(  const bool &      stop )
{
  for (index_t w = 1; w!=m_proc_count; w++)
    this->_sendMainToWorker(w,stop);
}

template <class T>
void gsAPALM<T>::_recvMainToWorker( const index_t &   sourceID,
                                          index_t &   branch,
                                          index_t &   jobID,
                                          index_t &   dataLevel)
{
  gsMPIInfo(m_rank)<<"Receiving data on "<<m_rank<<" from "<<sourceID<<"\n";
  MPI_Request req[3];
  m_comm.irecv(&branch    ,1,sourceID,&req[0] ,0 );
  m_comm.irecv(&jobID     ,1,sourceID,&req[1] ,1 );
  m_comm.irecv(&dataLevel ,1,sourceID,&req[2] ,2 );
  MPI_Waitall( 3, req, MPI_STATUSES_IGNORE );
}


template <class T>
void gsAPALM<T>::_recvMainToWorker( const index_t &         sourceID,
                                          std::tuple<index_t, T     , solution_t, solution_t> & dataEntry,
                                          std::pair<T,T> &  dataInterval,
                                          solution_t &      dataReference)
{
  gsMPIInfo(m_rank)<<"Receiving data on "<<m_rank<<" from "<<sourceID<<"\n";

  index_t     ID;
  T           dL0;
  solution_t  start,    prev;
  gsVector<T> startU,   prevU,    refU;
  index_t     startSize,prevSize, refSize;
  T           startL,   prevL,    refL;

  T           tstart, tend;

  // Put all data 0-14 in a struct and recv a single struct
  MPI_Request req[10];
  m_comm.irecv(&ID        ,1,sourceID,&req[0] ,0 );
  m_comm.irecv(&dL0       ,1,sourceID,&req[1] ,1 );
  m_comm.irecv(&tstart    ,1,sourceID,&req[2] ,2 );
  m_comm.irecv(&tend      ,1,sourceID,&req[3] ,3 );

  m_comm.irecv(&startSize ,1,sourceID,&req[4] ,4 );
  m_comm.irecv(&startL    ,1,sourceID,&req[5] ,5 );

  m_comm.irecv(&prevSize  ,1,sourceID,&req[6] ,6 );
  m_comm.irecv(&prevL     ,1,sourceID,&req[7] ,7 );

  m_comm.irecv(&refSize   ,1,sourceID,&req[8] ,8 );
  m_comm.irecv(&refL      ,1,sourceID,&req[9] ,9 );

  MPI_Waitall( 10, req, MPI_STATUSES_IGNORE );
  dataInterval = std::make_pair(tstart,tend);

  startU.resize(startSize);
  prevU .resize(prevSize);
  refU  .resize(refSize);

  MPI_Request req_data[3];
  m_comm.irecv(startU.data(), startSize,sourceID,&req_data[0],10);
  m_comm.irecv(prevU.data(),  prevSize, sourceID,&req_data[1],11);
  m_comm.irecv(refU.data(),   refSize,  sourceID,&req_data[2],12);

  MPI_Waitall( 3, req_data, MPI_STATUSES_IGNORE );

  start         = std::make_pair(startU,startL);
  prev          = std::make_pair(prevU ,prevL);
  dataReference = std::make_pair(refU  ,refL);

  dataEntry = std::make_tuple(ID,dL0,start,prev);
}

template <class T>
void gsAPALM<T>::_recvMainToWorker( const index_t &   sourceID,
                                          std::tuple<index_t, T     , solution_t, solution_t> & dataEntry,
                                          T &               startTime)
{
  gsMPIInfo(m_rank)<<"Receiving data on "<<m_rank<<" from "<<sourceID<<"\n";

  index_t     ID;
  T           dL0;
  solution_t  start,    prev;
  gsVector<T> startU,   prevU;
  index_t     startSize,prevSize;
  T           startL,   prevL;

  // Put all data 0-14 in a struct and recv a single struct
  MPI_Request req[7];
  m_comm.irecv(&ID        ,1,sourceID,&req[0] ,0 );
  m_comm.irecv(&dL0       ,1,sourceID,&req[1] ,1 );
  m_comm.irecv(&startTime ,1,sourceID,&req[2] ,2 );

  m_comm.irecv(&startSize ,1,sourceID,&req[3] ,3 );
  m_comm.irecv(&startL    ,1,sourceID,&req[4] ,4 );

  m_comm.irecv(&prevSize  ,1,sourceID,&req[5] ,5 );
  m_comm.irecv(&prevL     ,1,sourceID,&req[6] ,6 );

  MPI_Waitall( 7, req, MPI_STATUSES_IGNORE );

  startU.resize(startSize);
  prevU .resize(prevSize);

  MPI_Request req_data[2];
  m_comm.irecv(startU.data(), startU.size(),sourceID,&req_data[0],7);
  m_comm.irecv(prevU.data(),  prevU.size(), sourceID,&req_data[1],8);

  MPI_Waitall( 2, req_data, MPI_STATUSES_IGNORE );

  start         = std::make_pair(startU,startL);
  prev          = std::make_pair(prevU ,prevL);

  dataEntry = std::make_tuple(ID,dL0,start,prev);
}

template <class T>
void gsAPALM<T>::_recvMainToWorker(  const  index_t &   sourceID,
                                            bool &      stop)
{
  // when we handle the stop signal with isend/irecv, we get a deadlock
  m_comm.recv(&stop,1,sourceID,99);
}

// Sends metadata
template <class T>
void gsAPALM<T>::_sendWorkerToMain( const index_t &                   mainID,
                                    const index_t &                   branch,
                                    const index_t &                   jobID)
{
  gsMPIInfo(m_rank)<<"Sending data from "<<m_rank<<" to "<<mainID<<"\n";

  index_t tag = 0;

  MPI_Request req_meta[2];
  m_comm.isend(&branch        ,1, mainID, &req_meta[0], tag++);
  m_comm.isend(&jobID         ,1, mainID, &req_meta[1], tag++);
  MPI_Waitall( 2, req_meta, MPI_STATUSES_IGNORE );
}

// Receives metadata
template <class T>
void gsAPALM<T>::_recvWorkerToMain( index_t &                   sourceID,
                                    index_t &                   branch,
                                    index_t &                   jobID)
{

  index_t tag = 0;

  MPI_Status  status;
  MPI_Request req_meta[2];
  m_comm.irecv(&branch         ,1, sourceID, &req_meta[0], tag++);
  MPI_Wait ( &req_meta[0], &status );
  if (sourceID==MPI_ANY_SOURCE)  sourceID = status.MPI_SOURCE; // Overwrite source ID to be unique

  // Request metadata
  m_comm.irecv(&jobID          ,1, sourceID, &req_meta[1], tag++);
  MPI_Waitall( 2, req_meta, MPI_STATUSES_IGNORE );
}

template <class T>
void gsAPALM<T>::_sendWorkerToMain( const index_t &                   mainID,
                                    const std::vector<T> &            distances,
                                    const std::vector<solution_t> &   stepSolutions,
                                    const T &                         upperDistance,
                                    const T &                         lowerDistance )
{
  gsMPIInfo(m_rank)<<"Sending data from "<<m_rank<<" to "<<mainID<<"\n";

  index_t tag = 0;
  index_t size = stepSolutions.size();
  index_t vectorSize;
  T load;

  MPI_Request req_meta[2];
  m_comm.isend(&upperDistance ,1, mainID, &req_meta[0], tag++);
  m_comm.isend(&lowerDistance ,1, mainID, &req_meta[1], tag++);

  MPI_Request req_solSize;
  m_comm.isend(&size          ,1, mainID, &req_solSize, tag++);

  MPI_Request req_vecSizes[size];
  MPI_Request req_loads[size];
  MPI_Request req_distances[size+1];
  for (index_t k=0; k!=size; k++)
  {
    m_comm.isend(&distances.at(k)   ,1, mainID, &req_distances[k]   , tag++);
    load      = stepSolutions.at(k).second;
    m_comm.isend(&load              ,1, mainID, &req_loads[k]       , tag++);
    vectorSize= stepSolutions.at(k).first.size();
    m_comm.isend(&vectorSize        ,1, mainID, &req_vecSizes[k]    , tag++);
  }
  m_comm.isend(&(distances.at(size)),1, mainID, &req_distances[size], tag++);

  MPI_Request req_data[size];
  for (index_t k=0; k!=size; k++)
    m_comm.isend( stepSolutions.at(k).first.data(),
                  stepSolutions.at(k).first.size(), mainID, &req_data[k], tag++);

  MPI_Wait(&req_solSize,MPI_STATUS_IGNORE);
  MPI_Waitall( size, req_vecSizes, MPI_STATUSES_IGNORE );
  MPI_Waitall( size, req_data, MPI_STATUSES_IGNORE );
  MPI_Waitall( size, req_loads, MPI_STATUSES_IGNORE );
  MPI_Waitall( size+1, req_distances, MPI_STATUSES_IGNORE );
  MPI_Waitall( 2, req_meta, MPI_STATUSES_IGNORE );
}

template <class T>
void gsAPALM<T>::_sendWorkerToMain( const index_t &                   mainID,
                                    const T &                         distance,
                                    const std::vector<solution_t> &   solutions,
                                    const bool &                      bifurcation)
{
  gsMPIInfo(m_rank)<<"Sending data from "<<m_rank<<" to "<<mainID<<"\n";

  index_t tag = 0;
  index_t size = solutions.size();
  index_t vectorSize;
  T load;

  MPI_Request req_meta[2];
  m_comm.isend(&distance ,1, mainID, &req_meta[0], tag++);
  m_comm.isend(&bifurcation ,1, mainID, &req_meta[1], tag++);

  MPI_Request req_solSize;
  m_comm.isend(&size          ,1, mainID, &req_solSize, tag++);

  MPI_Request req_vecSizes[size];
  MPI_Request req_loads[size];
  for (index_t k=0; k!=size; k++)
  {
    load      = solutions.at(k).second;
    m_comm.isend(&load              ,1, mainID, &req_loads[k]       , tag++);
    vectorSize= solutions.at(k).first.size();
    m_comm.isend(&vectorSize        ,1, mainID, &req_vecSizes[k]    , tag++);
  }

  MPI_Request req_data[size];
  for (index_t k=0; k!=size; k++)
    m_comm.isend( solutions.at(k).first.data(),
                  solutions.at(k).first.size(), mainID, &req_data[k], tag++);

  MPI_Wait(&req_solSize,MPI_STATUS_IGNORE);
  MPI_Waitall( size, req_vecSizes, MPI_STATUSES_IGNORE );
  MPI_Waitall( size, req_data, MPI_STATUSES_IGNORE );
  MPI_Waitall( size, req_loads, MPI_STATUSES_IGNORE );
  MPI_Waitall( 2, req_meta, MPI_STATUSES_IGNORE );
}

template <class T>
void gsAPALM<T>::_recvWorkerToMain( index_t &                   sourceID,
                                    std::vector<T> &            distances,
                                    std::vector<solution_t> &   stepSolutions,
                                    T &                         upperDistance,
                                    T &                         lowerDistance)
{

  index_t tag = 0;
  index_t size;

  MPI_Status  status;
  MPI_Request req_meta[2];
  m_comm.irecv(&upperDistance  ,1, sourceID, &req_meta[0], tag++);
  MPI_Wait ( &req_meta[0], &status );
  if (sourceID==MPI_ANY_SOURCE)  sourceID = status.MPI_SOURCE; // Overwrite source ID to be unique
  gsMPIInfo(m_rank)<<"Receiving data on "<<m_rank<<" from "<<sourceID<<"\n";

  // Request metadata
  m_comm.irecv(&lowerDistance  ,1, sourceID, &req_meta[1], tag++);

  // Request the number of solution intervals
  MPI_Request req_solSize;
  m_comm.irecv(&size           ,1, sourceID, &req_solSize, tag++);

  MPI_Wait(&req_solSize,MPI_STATUS_IGNORE);

  // temporary objects to store solutions and the sizes of the solution vectors
  std::vector<gsVector<T>> Us(size);
  std::vector<T          > Ls(size);
  std::vector<index_t    > vectorSizes(size);
  // final containers for distances and step solutions
  distances.resize(size+1);
  stepSolutions.resize(size);

  // Collect distances, solution vector sizes and loads
  MPI_Request req_vecSizes[size];
  MPI_Request req_loads[size];
  MPI_Request req_distances[size+1];
  for (index_t k=0; k!=size; k++)
  {
    m_comm.irecv(&distances.at(k)   ,1, sourceID, &req_distances[k]   , tag++);
    m_comm.irecv(&Ls.at(k)          ,1, sourceID, &req_loads[k]       , tag++);
    m_comm.irecv(&vectorSizes.at(k) ,1, sourceID, &req_vecSizes[k]    , tag++);
  }
  m_comm.irecv(&(distances.at(size)),1, sourceID, &req_distances[size], tag++);

  // When all solution vector sizes are collected, resize the solution vectors
  MPI_Waitall( size, req_vecSizes, MPI_STATUSES_IGNORE );
  for (index_t k=0; k!=size; k++)
    Us[k].resize(vectorSizes.at(k));

  // Collect the solution data
  MPI_Request req_data[size];
  for (index_t k=0; k!=size; k++)
    m_comm.irecv(Us.at(k).data()    ,vectorSizes.at(k)  , sourceID, &req_data[k], tag++);

  // When the solution data and the loads are collected, put them in the stepSolutions container.
  MPI_Waitall( size, req_data, MPI_STATUSES_IGNORE );
  MPI_Waitall( size, req_loads, MPI_STATUSES_IGNORE );
  for (index_t k=0; k!=size; k++)
    stepSolutions.at(k) = std::make_pair(Us.at(k),Ls.at(k));

  // Wait for the distances to finish
  MPI_Waitall( size+1, req_distances, MPI_STATUSES_IGNORE );
  // Wait for the meta data to finish
  MPI_Waitall( 2, req_meta, MPI_STATUSES_IGNORE );
}

template <class T>
void gsAPALM<T>::_recvWorkerToMain( index_t &                 sourceID,
                                    T &                       distance,
                                    std::vector<solution_t> & solutions,
                                    bool &                    bifurcation)
{
  index_t tag = 0;
  index_t size;

  MPI_Status  status;
  MPI_Request req_meta[2];
  m_comm.irecv(&distance  ,1, sourceID, &req_meta[0], tag++);
  MPI_Wait ( &req_meta[0], &status );
  if (sourceID==MPI_ANY_SOURCE)  sourceID = status.MPI_SOURCE; // Overwrite source ID to be unique
  gsMPIInfo(m_rank)<<"Receiving data on "<<m_rank<<" from "<<sourceID<<"\n";

  // Request metadata
  m_comm.irecv(&bifurcation  ,1, sourceID, &req_meta[1], tag++);

  // Request the number of solution intervals
  MPI_Request req_solSize;
  m_comm.irecv(&size           ,1, sourceID, &req_solSize, tag++);

  MPI_Wait(&req_solSize,MPI_STATUS_IGNORE);

  // temporary objects to store solutions and the sizes of the solution vectors
  std::vector<gsVector<T>> Us(size);
  std::vector<T          > Ls(size);
  std::vector<index_t    > vectorSizes(size);
  // final containers for distances and step solutions
  solutions.resize(size);

  // Collect distances, solution vector sizes and loads
  MPI_Request req_vecSizes[size];
  MPI_Request req_loads[size];
  for (index_t k=0; k!=size; k++)
  {
    m_comm.irecv(&Ls.at(k)          ,1, sourceID, &req_loads[k]       , tag++);
    m_comm.irecv(&vectorSizes.at(k) ,1, sourceID, &req_vecSizes[k]    , tag++);
  }

  // When all solution vector sizes are collected, resize the solution vectors
  MPI_Waitall( size, req_vecSizes, MPI_STATUSES_IGNORE );
  for (index_t k=0; k!=size; k++)
    Us[k].resize(vectorSizes.at(k));

  // Collect the solution data
  MPI_Request req_data[size];
  for (index_t k=0; k!=size; k++)
    m_comm.irecv(Us.at(k).data()    ,vectorSizes.at(k)  , sourceID, &req_data[k], tag++);

  // When the solution data and the loads are collected, put them in the solutions container.
  MPI_Waitall( size, req_data, MPI_STATUSES_IGNORE );
  MPI_Waitall( size, req_loads, MPI_STATUSES_IGNORE );
  for (index_t k=0; k!=size; k++)
    solutions.at(k) = std::make_pair(Us.at(k),Ls.at(k));

  // Wait for the meta data to finish
  MPI_Waitall( 2, req_meta, MPI_STATUSES_IGNORE );
}


#endif

// -----------------------------------------------------------------------------------------------------
} // namespace gismo
