 /** @file gsArcLengthIterator.h

    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#pragma once
#include <gsStructuralAnalysis/gsSpaceTimeHierarchy.h>
#include <gsStructuralAnalysis/gsALMBase.h>

namespace gismo
{

/**
    @brief

    \tparam T coefficient type

    \ingroup gsStructuralAnalysis
*/
template <class T>
class gsAdaptiveALM
{
    typedef std::pair<gsVector<T>,T> solution_t;

public:

    ~gsAdaptiveALM() {};

    ///
    gsAdaptiveALM( gsALMBase<T> * ALM )
    :
    m_ALM(ALM),
    m_hierarchy()
    {
        m_ALM->initialize();
        m_verbose = true;
    }

    ///
    gsAdaptiveALM( gsALMBase<T> * ALM, const std::vector<T> times, const std::vector<solution_t> solutions )
    :
    m_ALM(ALM),
    m_times(times),
    m_solutions(solutions),
    m_hierarchy(times,solutions)
    {
        m_hierarchy.init();
        m_verbose = true;
    }

// General functions
public:
    void step()
    {
        if (!m_hierarchy.empty())
        {
            // FOR NON-SPLIT!
            std::tie(m_ID,m_s,m_ds,m_start,m_guess) = m_hierarchy.pop();
            std::tie(m_Uold,m_Lold) = m_start;
            std::tie(m_Uguess,m_Lguess) = m_guess;

            m_ALM->setLength(m_ds/2);
            m_ALM->setSolution(m_Uold,m_Lold);
            m_ALM->resetStep();
            m_ALM->setInitialGuess(m_Uguess,m_Lguess);

            if (m_verbose) gsInfo<<"Starting with ID "<<m_ID<<" from (lvl,|U|,L) = ("<<m_hierarchy.currentLevel(m_ID)<<","<<m_Uold.norm()<<","<<m_Lold<<"), curve time = "<<m_s<<"\n";

            m_ALM->step();
            m_hierarchy.submitLeft(m_ID,std::make_pair(m_ALM->solutionU(),m_ALM->solutionL()));
            m_ALM->step();
            m_hierarchy.submitRight(m_ID,std::make_pair(m_ALM->solutionU(),m_ALM->solutionL()));

            bool success = m_hierarchy.getReference(m_ID,m_reference);
            if (success)
            {
              std::tie(m_Uref,m_Lref) = m_reference;
            }
            else
                GISMO_ERROR("Reference not found; something went wrong.");

            m_DeltaU = m_Uref - m_ALM->solutionU();
            m_DeltaL = m_Lref - m_ALM->solutionL();

            m_error = m_ALM->distance(m_DeltaU,m_DeltaL) / (m_ds);


            if (m_error > m_tolerance)
            {
              gsDebug<<"ERROR > TOL\n";
              m_hierarchy.addJob(m_ID);
            }

            m_hierarchy.removeJob(m_ID);
        }
        else
        {
            gsWarn<<"Hierarchy is empty\n";
        }
    }

    void solve()
    {
        m_it = 0;
        m_maxIt = 100;
        while (!m_hierarchy.empty() && m_it < m_maxIt)
        {
            this->step();
        }
    }

    std::pair<std::vector<T>,std::vector<solution_t>> getSolution()
    {
        return m_hierarchy.getFlatSolution();
    }

    void _init()
    {

    }

    void initialize(index_t N = 10)
    {
        // Initialize solution
        m_U.setZero(m_ALM->numDofs(),1);
        m_L = 0;

        m_solutions.push_back({m_U,m_L});
        m_times.push_back(m_s);

        m_s = 0;
        m_ds = m_ALM->getLength();
        for (index_t k=0; k<N; k++)
        {
            m_s += m_ds;

            if (m_verbose) gsInfo<<"Load step "<< k<<"\t"<<"dL = "<<m_ds<<"; curve time = "<<m_s<<"\n";
            // assembler->constructSolution(solVector,solution);
            m_ALM->step();

            if (!(m_ALM->converged()))
              GISMO_ERROR("Loop terminated, arc length method did not converge.\n");

            m_L = m_ALM->solutionL();
            m_solutions.push_back({m_ALM->solutionU(),m_L});
            m_times.push_back(m_s);
        }

        _makeHierarchy();
    }

    void _makeHierarchy()
    {
        m_hierarchy = gsSpaceTimeHierarchy<T,solution_t> (m_times,m_solutions);
        m_hierarchy.init();
    }

    void setHierarchy(const gsSpaceTimeHierarchy<T,solution_t> & hierarchy)
    {
        m_hierarchy = hierarchy;
        m_hierarchy.init();
    }

    void printTree() { m_hierarchy.printTree(); }
    void printQueue() { m_hierarchy.printQueue(); }


    void computeFit()
    {

    }


protected:
    gsALMBase<T> * m_ALM;
    gsOptionList m_Hoptions;
    gsOptionList m_ALMoptions;

    gsSpaceTimeHierarchy<T,solution_t> m_hierarchy;

    bool m_verbose;

    std::vector<T> m_times;
    std::vector<solution_t> m_solutions;

    gsVector<T> m_U, m_Uold,m_Uguess,m_Uref,m_DeltaU;
    T m_L, m_Lold, m_Lguess,m_Lref,m_DeltaL;

    index_t m_maxIt, m_ID, m_it;
    T m_tolerance, m_s, m_ds;
    T m_error;

    solution_t m_start, m_guess, m_reference;

};


} // namespace gismo

// #ifndef GISMO_BUILD_LIB
// #include GISMO_HPP_HEADER(gsAdaptiveALM.hpp)
// #endif