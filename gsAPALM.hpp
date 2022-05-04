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
  GISMO_ASSERT(m_starts.size()>0,"No start point is created. Call initialize first?");

  T Lold, L0;
  gsMatrix<T> Uold, U0;
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
      if (m_verbose) gsInfo<<"Load step "<< k<<"\t"<<"dL = "<<dL<<"; curve time = "<<s<<"\n";

      // Set a step
      m_ALM->step();
      // If not converged, bisect the arc-length
      if (!(m_ALM->converged()))
      {
        if (m_verbose) gsInfo<<"Error: Loop terminated, arc length method did not converge.\n";
        m_ALM->reduceLength();
        m_ALM->setSolution(Uold,Lold);
        bisected = true;
        continue;
      }

      if (m_singularPoint)
      {
        m_ALM->computeStability(m_ALM->solutionU(),m_ALM->options().getSwitch("Quasi"));
        if (m_ALM->stabilityChange())
        {
          gsInfo<<"Bifurcation spotted!"<<"\n";
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

      if (!bisected)
      {
        m_ALM->resetLength();
      }
      bisected = false;
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

}

// NOTE: This does not make new branches!
template <class T>
void gsAPALM<T>::parallelSolve()
{
  solution_t start, prev, next, reference;
  index_t ID;
  T tstart = 0;
  T tend = 0;
  T dL, dL0;
  index_t it = 0;
  gsVector<T> DeltaU;
  T DeltaL;
  index_t Nintervals = m_subIntervals;
  std::vector<solution_t> stepSolutions;
  std::vector<T> stepTimes;
  std::vector<T> distances;
  T lowerDistance, upperDistance;
  bool bisected = false;
  T dL_rem = 0;

  T Lprev, Lnext, Lold;
  gsMatrix<T> Uprev, Unext,Uold;
  index_t level;
  index_t branch;
  m_data.print();
  while (!m_data.empty() && it < m_maxIterations)
  {
    branch = m_data.getFirstNonEmptyBranch();
    Nintervals = m_subIntervals;
    gsInfo<<"There are "<<m_data.branch(branch).nActive()<<" active jobs and "<<m_data.branch(branch).nWaiting()<<" jobs in the queue of branch "<<branch<<"\n";
    stepSolutions.resize(Nintervals);
    stepTimes.resize(Nintervals);
    distances.resize(Nintervals+1);

    std::tie(ID,dL0,start,prev,next) = m_data.branch(branch).pop();
    std::tie(Uold,Lold) = start;
    std::tie(Uprev,Lprev) = prev;
    std::tie(Unext,Lnext) = next;
    std::tie(tstart,tend) = m_data.branch(branch).jobTimes(ID);
    level                 = m_data.branch(branch).jobLevel(ID);

    gsMatrix<T> Uori = Uold;
    T Lori = Lold;

    dL0 = dL0 / Nintervals;
    dL  = dL0;

    m_ALM->setLength(dL);
    m_ALM->setSolution(Uold,Lold);
    // m_ALM->resetStep();
    m_ALM->setPrevious(Uprev,Lprev);
    gsDebug<<"Start - ||u|| = "<<Uold.norm()<<", L = "<<Lold<<"\n";
    gsDebug<<"Prev  - ||u|| = "<<Uprev.norm()<<", L = "<<Lprev<<"\n";
    gsDebug<<"Next  - ||u|| = "<<Unext.norm()<<", L = "<<Lnext<<"\n";

    gsVector<T> tmpU = Uold-Uprev;
    T tmpL = Lold-Lprev;


    T s = 0;
    T time = tstart;

    gsInfo<<"Starting with ID "<<ID<<" from (|U|,L) = ("<<Uold.norm()<<","<<Lold<<"), curve time = "<<time<<"\n";
    for (index_t k = 0; k!=Nintervals; k++)
    {
      gsDebug<<"Interval "<<k+1<<" of "<<Nintervals<<"\n";
      gsDebug<<"Start - ||u|| = "<<Uold.norm()<<", L = "<<Lold<<"\n";
      m_ALM->step();
      if (!(m_ALM->converged()))
      {
        gsInfo<<"Error: Loop terminated, arc length method did not converge.\n";
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

    bool success = m_data.branch(branch).getReferenceByID(ID,reference);
    GISMO_ASSERT(success,"Reference not found");

    GISMO_ASSERT((reference.first.normalized()).dot(m_ALM->solutionU().normalized())>-0.8,
                  "Reference is almost in the opposite direction of the solution. Are branches mirrored? result:" << (reference.first.normalized()).dot(m_ALM->solutionU().normalized()));

    DeltaU = reference.first - m_ALM->solutionU();
    DeltaL = reference.second - m_ALM->solutionL();
    distances.back() = m_ALM->distance(DeltaU,DeltaL);

    DeltaU = reference.first - Uori;
    DeltaL = reference.second - Lori;
    upperDistance = m_ALM->distance(DeltaU,DeltaL);

    DeltaU = stepSolutions.back().first - Uori;
    DeltaL = stepSolutions.back().second - Lori;
    lowerDistance = m_ALM->distance(DeltaU,DeltaL);

    if (upperDistance > 0.1)

  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////

    m_data.branch(branch).submit(ID,distances,stepSolutions,upperDistance,lowerDistance);
    m_data.branch(branch).finishJob(ID);

    this->parallelIntervalOutput(stepSolutions,stepTimes,level,ID);

    it++;
  }

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

} // namespace gismo