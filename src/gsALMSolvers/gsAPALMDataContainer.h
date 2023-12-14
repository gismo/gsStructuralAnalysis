 /** @file gsAPALMDataContainer.h

    @brief Performs the arc length method to solve a nonlinear equation system.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsNurbs/gsKnotVector.h>
#include <gsIO/gsOptionList.h>
#include <gsHSplines/gsKdNode.h>
#include <queue>

namespace gismo
{

template<class T, class solution_t >
class gsAPALMDataContainer
{

public:

  ~gsAPALMDataContainer() { }

  gsAPALMDataContainer() { }

  gsAPALMDataContainer(const gsAPALMData<T,solution_t> & data)
  {
    this->add(data);
  }

  bool empty()
  {
    bool result = true;
    for (typename std::vector<gsAPALMData<T,solution_t>>::iterator it=m_container.begin(); it!=m_container.end(); it++)
      result &= it->empty();
    return result;
  }

  index_t nBranches()
  {
    return m_container.size();
  }

  // std::tuple<std::vector<T>,std::vector<solution_t>,std::vector<index_t>> getFlatSolution(index_t level=-1);

  index_t add(const gsAPALMData<T,solution_t> & data)
  {
    m_container.push_back(data);
    return m_container.size()-1;
  }

  gsAPALMData<T,solution_t> & branch(index_t k)
  {
    return m_container.at(k);
  }

  index_t getFirstNonEmptyBranch()
  {
    index_t k=0;
    while (m_container[k].empty() && (size_t)k < m_container.size())
      k++;
    return k;
  }

  void print()
  {
    index_t k=0;
    for (typename std::vector<gsAPALMData<T,solution_t>>::iterator it=m_container.begin(); it!=m_container.end(); it++, k++)
    {
      gsInfo<<"-------------------------------------------\n";
      gsInfo<<"Branch "<<k<<"\n";
      it->print();
    }
  }

  void printQueue()
  {
    index_t k=0;
    for (typename std::vector<gsAPALMData<T,solution_t>>::iterator it=m_container.begin(); it!=m_container.end(); it++, k++)
    {
      gsInfo<<"-------------------------------------------\n";
      gsInfo<<"Branch "<<k<<"\n";
      it->printQueue();
    }
  }

  void printKnots()
  {
    index_t k=0;
    for (typename std::vector<gsAPALMData<T,solution_t>>::iterator it=m_container.begin(); it!=m_container.end(); it++, k++)
    {
      gsInfo<<"-------------------------------------------\n";
      gsInfo<<"Branch "<<k<<"\n";
      it->printKnots();
    }
  }

  size_t nActive()
  {
    index_t k=0;
    for (typename std::vector<gsAPALMData<T,solution_t>>::iterator it=m_container.begin(); it!=m_container.end(); it++)
      k+= it->nActive();
    return k;
  }
  size_t nWaiting()
  {
    index_t k=0;
    for (typename std::vector<gsAPALMData<T,solution_t>>::iterator it=m_container.begin(); it!=m_container.end(); it++)
      k+= it->nWaiting();
    return k;
  }

  size_t maxLevel()
  {
    size_t k=0;
    for (typename std::vector<gsAPALMData<T,solution_t>>::iterator it=m_container.begin(); it!=m_container.end(); it++)
      k = it->maxLevel() > k ? it->maxLevel() : k;
    return k;
  }

protected:
  std::vector<gsAPALMData<T,solution_t>> m_container;

};

}

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsAPALMDataContainer.hpp)
#endif
