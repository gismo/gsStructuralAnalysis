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

#ifdef GISMO_WITH_MPI
#include <gsMpi/gsMpi.h>
#endif

namespace gismo
{

#define gsMPIInfo(rank) gsInfo<<"[MPI process "<<rank<<"]: "
#define gsMPIDebug(rank) gsDebug<<"[MPI process "<<rank<<"]: "

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

  gsAPALM()
#ifdef GISMO_WITH_MPI
  :m_mpi(gsMpi::init())
#endif
  {};

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
    for (size_t l=0; l!=m_lvlSolutions[branch].size(); l++)
      for (typename std::vector<solution_t *>::iterator it=m_lvlSolutions[branch][l].begin(); it!=m_lvlSolutions[branch][l].end(); it++)
        result[l].push_back(**it);
    return result;
  }

  std::vector<std::vector<T>> getTimesPerLevel(index_t branch = 0)
  {
    std::vector<std::vector<T>> result(m_lvlTimes[branch].size());
    for (size_t l=0; l!=m_lvlSolutions[branch].size(); l++)
      for (typename std::vector<T *>::iterator it=m_lvlTimes[branch][l].begin(); it!=m_lvlTimes[branch][l].end(); it++)
        result[l].push_back(**it);
    return result;
  }

#ifdef GISMO_WITH_MPI
  bool isMain() {return (m_rank==0); }
  index_t rank() {return (m_rank); }
  const gsMpi & mpi() { return m_mpi; }
  real_t wallTime() { return m_mpi.wallTime(); }
#else
  bool isMain() {return true; }
  index_t rank() {return 0; }
#endif

protected:

#ifdef GISMO_WITH_MPI
  // For data
  void _sendMainToWorker( const index_t &         workerID,
                          const index_t &         branch,
                          const std::tuple<index_t, T     , solution_t, solution_t, solution_t> & dataEntry,
                          const std::pair<T,T> &  dataInterval,
                          const index_t &         dataLevel,
                          const solution_t &      dataReference );

  // For stop signal
  void _sendMainToWorker( const index_t &   workerID,
                          const bool &      stop      );
  void _sendMainToAll(    const bool &      stop      );

  // For data
  void _recvMainToWorker( const index_t &         sourceID,
                                index_t &         branch,
                                std::tuple<index_t, T     , solution_t, solution_t, solution_t> & dataEntry,
                                std::pair<T,T> &  dataInterval,
                                index_t &         dataLevel,
                                solution_t &      dataReference,
                                MPI_Status *      status = NULL);

  // For stop signal
  void _recvMainToWorker(   const index_t &   sourceID,
                                  bool &      stop,
                                  MPI_Status *status = NULL );

  // For data
  void _sendWorkerToMain( const index_t &                   mainID,
                          const index_t &                   branch,
                          const index_t &                   jobID,
                          const std::vector<T> &            distances,
                          const std::vector<solution_t> &   stepSolutions,
                          const T &                         upperDistance,
                          const T &                         lowerDistance );

  // For data
  void _recvWorkerToMain( index_t &                   sourceID, // source ID will change to MPI
                          index_t &                   branch,
                          index_t &                   jobID,
                          std::vector<T> &            distances,
                          std::vector<solution_t>&    stepSolutions,
                          T &                         upperDistance,
                          T &                         lowerDistance,
                          MPI_Status *                status = NULL);

  gsMpiComm & comm() { return m_comm; }
#endif

  template <bool _hasWorkers>
  typename std::enable_if< _hasWorkers, void>::type
  parallelSolve_impl();

  template <bool _hasWorkers>
  typename std::enable_if<!_hasWorkers, void>::type
  parallelSolve_impl();

  void _parallelSolve_worker(   const std::tuple<index_t, T     , solution_t, solution_t, solution_t> & dataEntry,
                                        const std::pair<T,T> &  dataInterval,
                                        const index_t &         dataLevel,
                                        const solution_t &      dataReference,
                                        std::vector<T> &        distances,
                                        std::vector<solution_t>&stepSolutions,
                                        T &                     upperDistance,
                                        T &                     lowerDistance   );
  void _parallelSolve_main(     gsAPALMDataContainer<T,solution_t> & data       );

  void _parallelSolve_finalize();

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

  index_t m_maxIterations;

  // Conditional compilation
#ifdef GISMO_WITH_MPI
  const gsMpi & m_mpi;
  mutable gsMpiComm m_comm ;
  std::queue<index_t> m_workers;
#endif

  index_t m_proc_count, m_rank;
};

}

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsAPALM.hpp)
#endif
