/** @file gsAPALMBase.hpp

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
gsAPALMBase<T,solution_t>::gsAPALMBase(       gsALMBase<T> * ALM,
                                        const gsAPALMData<T,solution_t> & Data)
:
m_ALM(ALM),
m_data(Data)
{
  this->_defaultOptions();
}

template <class T, class solution_t >
void gsAPALMBase<T,solution_t>::_defaultOptions()
{
  m_options.addSwitch("Verbose","Verbosity",false);
}

template <class T, class solution_t >
void gsAPALMBase<T,solution_t>::_getOptions()
{
  m_verbose = m_options.getSwitch("verbose");
}

template <class T, class solution_t >
void gsAPALMBase<T,solution_t>::initialize()
{
  this->_getOptions();
}

template <class T, class solution_t >
void gsAPALMBase<T,solution_t>::serialSolve(index_t Nsteps)
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
  m_ALM->setSolution(U0,L0);

  for (index_t k=0; k!=Nsteps; k++)
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
    m_solutions.push_back({m_ALM->solutionU(),lambda});
    m_times.push_back(s);

    this->serialStepOutput(m_ALM->solutionU(),lambda);

    // // add data point to fit
    // index_t N = m_ALM->solutionU().rows();
    // gsVector<T> sol(N+1);
    // sol.head(N) = m_ALM->solutionU();
    // sol.at(N) = lambda;
    // cfitter.addDataPoint(sol,s);
    // if (k>deg_z+1)
    //   cfitter.compute();
    //


    if (!bisected)
    {
      dL = dL0;
      m_ALM->setLength(dL);
    }
    bisected = false;

  } // end of steps
}

template <class T, class solution_t >
void gsAPALMBase<T,solution_t>::parallelSolve()
{
  solution_t start, guess, reference;
  index_t ID;
  real_t tstart = 0;
  real_t tend = 0;
  real_t dt, dt0;
  index_t it = 0;
  index_t itmax = 100;
  real_t TOL = 1e-2;
  gsVector<> DeltaU;
  real_t DeltaL;
  index_t Nintervals;
  std::vector<solution_t> stepSolutions;
  std::vector<real_t> stepTimes;
  std::vector<real_t> distances;
  bool bisected = false;
  real_t dt_rem = 0;

  T Lguess, Lold, L0;
  gsMatrix<T> Uguess,Uold, U0;

  while (!m_data.empty() && it < itmax)
  {
    Nintervals = 3;
    stepSolutions.resize(Nintervals);
    stepTimes.resize(Nintervals);
    distances.resize(Nintervals+1);

    std::tie(ID,tstart,tend,dt0,start,guess) = m_data.pop();
    std::tie(Uold,Lold) = start;
    std::tie(Uguess,Lguess) = guess;

    gsMatrix<> Uori = Uold;
    real_t Lori = Lold;

    dt0 = dt0 / Nintervals;
    dt = dt0;

    m_ALM->setLength(dt);
    m_ALM->setSolution(Uold,Lold);
    // m_ALM->resetStep();
    m_ALM->setPrevious(Uguess,Lguess);
    gsDebug<<"Start - ||u|| = "<<Uold.norm()<<", L = "<<Lold<<"\n";
    gsDebug<<"Guess - ||u|| = "<<Uguess.norm()<<", L = "<<Lguess<<"\n";

    gsVector<> tmpU = Uold-Uguess;
    real_t tmpL = Lold-Lguess;


    real_t s = 0;

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

      stepSolutions.at(k) = std::make_pair(m_ALM->solutionU(),m_ALM->solutionL());
      stepTimes.at(k) = tstart + dt;
      DeltaU = m_ALM->solutionU() - Uold;
      DeltaL = m_ALM->solutionL() - Lold;

      distances.at(k) = m_ALM->distance(DeltaU,DeltaL);
      gsDebugVar(m_ALM->distance(DeltaU,DeltaL));

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

      m_ALM->setLength(dt);
      dt_rem = 0;
      bisected = false;
    }

    gsDebugVar(s);

    bool success = m_data.getReferenceByID(ID,reference);
    GISMO_ASSERT(success,"Reference not found");

    DeltaU = reference.first - m_ALM->solutionU();
    DeltaL = reference.second - m_ALM->solutionL();
    distances.back() = m_ALM->distance(DeltaU,DeltaL);
    // gsDebugVar(m_ALM->distance(DeltaU,DeltaL));

    // DeltaU = m_ALM->solutionU() - Uori;
    // DeltaL = m_ALM->solutionL() - Lori;
    // gsDebugVar(m_ALM->distance(DeltaU,DeltaL));

    s += distances.back();
    // gsDebugVar(s);

  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////

    m_data.submit(ID,distances,stepTimes,stepSolutions);
    m_data.addJobs(ID);
    m_data.finishJob(ID);

    it++;
  }

  m_data.printQueue();
  std::tie(m_times,m_solutions) = m_data.getFlatSolution();

  m_data.printKnots();

  gsDebugVar(gsAsVector(m_times));

}

} // namespace gismo