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
m_data(Data)
{
  this->_defaultOptions();
  this->_getOptions();
}

template <class T>
void gsAPALM<T>::_defaultOptions()
{
  m_options.addSwitch("Verbose","Verbosity",false);
  m_options.addInt("SubIntervals","Number of subintervals",2);
  m_options.addSwitch("SingularPoint","Enable singular point detection",false);
  m_options.addReal("LengthMultiplier","Multiplier for the length after detection of a bifurcation point",1);
}

template <class T>
void gsAPALM<T>::_getOptions()
{
  m_verbose = m_options.getSwitch("Verbose");
  m_subIntervals = m_options.getInt("SubIntervals");
  m_singularPoint = m_options.getSwitch("SingularPoint");
  m_lengthMult = m_options.getReal("LengthMultiplier");

}

template <class T>
void gsAPALM<T>::initialize()
{
  this->_getOptions();
}

template <class T>
void gsAPALM<T>::serialSolve(index_t Nsteps)
{
  bool bisected = false;
  T s = 0;
  T dL = m_ALM->getLength();
  T dL0 = dL;

  // Initialize solutions
  T Lold, L0;
  gsMatrix<T> Uold, U0;
  Uold.setZero(m_ALM->numDofs(),1);
  U0.setZero(m_ALM->numDofs(),1);
  L0 = Lold = 0.0;

  // Store initial solution
  m_solutions.push_back({U0,L0});
  m_times.push_back(s);
  m_levels.push_back(0);
  m_ALM->setSolution(U0,L0);

  for (index_t k=1; k!=Nsteps; k++)
  {
    s+=dL;
    if (m_verbose) gsInfo<<"Load step "<< k<<"\t"<<"dL = "<<dL<<"; curve time = "<<s<<"\n";

    // Set a step
    m_ALM->step();
    // If not converged, bisect the arc-length
    if (!(m_ALM->converged()))
    {
      s -= dL;
      if (m_verbose) gsInfo<<"Error: Loop terminated, arc length method did not converge.\n";
      dL = dL / 2.;
      m_ALM->setLength(dL);
      m_ALM->setSolution(Uold,Lold);
      bisected = true;
      k -= 1;
      continue;
    }

    Uold = m_ALM->solutionU();
    Lold = m_ALM->solutionL();

    T lambda = m_ALM->solutionL();
    std::pair<gsVector<T>,T> pair = std::make_pair(m_ALM->solutionU(),lambda);
    m_solutions.push_back(pair);
    m_times.push_back(s);
    m_levels.push_back(0);

    this->serialStepOutput(pair,s,k);

    if (!bisected)
    {
      dL = dL0;
      m_ALM->setLength(dL);
    }
    bisected = false;

  } // end of steps

  m_data.setData(m_times,m_solutions);
  m_data.init();
  m_data.printQueue();
}

template <class T>
void gsAPALM<T>::parallelSolve()
{
  solution_t start, guess, reference;
  index_t ID;
  T tstart = 0;
  T tend = 0;
  T dt, dt0;
  index_t it = 0;
  index_t itmax = 100;
  T TOL = 1e-2;
  gsVector<T> DeltaU;
  T DeltaL;
  index_t Nintervals = m_subIntervals;
  std::vector<solution_t> stepSolutions;
  std::vector<T> stepTimes;
  std::vector<T> distances;
  T lowerError, upperError;
  bool bisected = false;
  T dt_rem = 0;

  T Lguess, Lold, L0;
  gsMatrix<T> Uguess,Uold, U0;
  index_t level;

  while (!m_data.empty() && it < itmax)
  {
    Nintervals = m_subIntervals;
    gsInfo<<"There are "<<m_data.nActive()<<" active jobs and "<<m_data.nWaiting()<<" jobs in the queue\n";
    stepSolutions.resize(Nintervals);
    stepTimes.resize(Nintervals);
    distances.resize(Nintervals+1);

    std::tie(ID,tstart,tend,dt0,start,guess,level) = m_data.pop();
    std::tie(Uold,Lold) = start;
    std::tie(Uguess,Lguess) = guess;

    gsMatrix<T> Uori = Uold;
    T Lori = Lold;

    dt0 = dt0 / Nintervals;
    dt = dt0;

    m_ALM->setLength(dt);
    m_ALM->setSolution(Uold,Lold);
    // m_ALM->resetStep();
    m_ALM->setPrevious(Uguess,Lguess);
    gsDebug<<"Start - ||u|| = "<<Uold.norm()<<", L = "<<Lold<<"\n";
    gsDebug<<"Guess - ||u|| = "<<Uguess.norm()<<", L = "<<Lguess<<"\n";

    gsVector<T> tmpU = Uold-Uguess;
    T tmpL = Lold-Lguess;


    T s = 0;

    gsInfo<<"Starting with ID "<<ID<<" from (|U|,L) = ("<<Uold.norm()<<","<<Lold<<"), curve time = "<<tstart<<"\n";
    for (index_t k = 0; k!=Nintervals; k++)
    {
      gsDebug<<"Interval "<<k+1<<" of "<<Nintervals<<"\n";
      gsDebug<<"Start - ||u|| = "<<Uold.norm()<<", L = "<<Lold<<"\n";
      m_ALM->step();
      if (!(m_ALM->converged()))
      {
        gsInfo<<"Error: Loop terminated, arc length method did not converge.\n";
        dt = dt / 2.;
        dt_rem += dt; // add the remainder of the interval to dt_rem
        m_ALM->setLength(dt);
        m_ALM->setSolution(Uold,Lold);
        bisected = true;
        k -= 1;
        continue;
      }
      GISMO_ENSURE(m_ALM->converged(),"Loop terminated, arc length method did not converge.\n");

      std::pair<gsVector<T>,T> pair = std::make_pair(m_ALM->solutionU(),m_ALM->solutionL());
      T time = tstart + dt;
      stepSolutions.at(k) = pair;
      stepTimes.at(k) = time;
      DeltaU = m_ALM->solutionU() - Uold;
      DeltaL = m_ALM->solutionL() - Lold;

      distances.at(k) = m_ALM->distance(DeltaU,DeltaL);

      s += distances.at(k);

      Uold = m_ALM->solutionU();
      Lold = m_ALM->solutionL();

      if (!bisected) // if bisected = false
        dt = dt0;
      else
      {
        // if the current interval has finished, but was refined before.
        // The next interval should have the remaining length.
        // Also, Nintervals should increase
        //
        dt = dt_rem;
        Nintervals++;
        stepSolutions.resize(Nintervals);
        stepTimes.resize(Nintervals);
        distances.resize(Nintervals+1);
      }

      this->parallelStepOutput(pair,time,k);

      m_ALM->setLength(dt);
      dt_rem = 0;
      bisected = false;
    }

    bool success = m_data.getReferenceByID(ID,reference);
    GISMO_ASSERT(success,"Reference not found");

    DeltaU = reference.first - m_ALM->solutionU();
    DeltaL = reference.second - m_ALM->solutionL();
    upperError = distances.back() = m_ALM->distance(DeltaU,DeltaL);

    DeltaU = reference.first - Uori;
    DeltaL = reference.second - Lori;
    T upperDistance = m_ALM->distance(DeltaU,DeltaL);

    DeltaU = stepSolutions.back().first - Uori;
    DeltaL = stepSolutions.back().second - Lori;
    T lowerDistance = m_ALM->distance(DeltaU,DeltaL);

    lowerError = s - lowerDistance;

  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////

    m_data.submit(ID,distances,stepSolutions,upperDistance,lowerDistance);
    m_data.finishJob(ID);

    this->parallelIntervalOutput(stepSolutions,stepTimes,level,ID);

    it++;
  }

  std::tie(m_times,m_solutions,m_levels) = m_data.getFlatSolution();
  m_lvlSolutions.resize(m_data.maxLevel()+1);
  m_lvlTimes.resize(m_data.maxLevel()+1);
  for (index_t k=0; k!=m_solutions.size(); k++)
  {
    GISMO_ASSERT(m_lvlSolutions.size() > (size_t)m_levels[k],"level mismatch, maxLevel = " << m_data.maxLevel() << "level = "<<m_levels[k]);
    GISMO_ASSERT(m_lvlTimes.size() > (size_t)m_levels[k],"level mismatch, maxLevel = " << m_data.maxLevel() << "level = "<<m_levels[k]);
    m_lvlSolutions[m_levels[k]].push_back(&(m_solutions[k]));
    m_lvlTimes[m_levels[k]].push_back(&(m_times[k]));
  }
}

} // namespace gismo