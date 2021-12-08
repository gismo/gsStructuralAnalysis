 /** @file gsControlDisplacement.h

    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#include <typeinfo>
#include <gsStructuralAnalysis/gsContinuationBase.h>
#pragma once


namespace gismo
{

/**
    @brief

    \tparam T

    \ingroup gsStructuralAnalysis
*/
template <class T>
class gsControlDisplacement : public gsContinuationBase<T>
{
protected:

public:

    gsControlDisplacement(gsStaticBase<T> * solver) :
    m_solver(solver),
    first(true)
    {
        m_U = gsVector<T>::Zero(m_solver->numDofs());
        m_L = 0;
    }

    void step(T dL)
    {
        m_solver->setLoad(m_L + dL);
        if (!first)
            m_solver->setDisplacement(m_U);
        else first = false;
        m_solver->solve();

        if (m_solver->converged())
        {
            m_L += dL;
            m_U += m_solver->update();
        }
    }

    gsVector<T> solutionU()
    {
        return m_U;
    }

    T solutionL()
    {
        return m_L;
    }

    void reset()
    {
        m_solver->reset();
    }

protected:
    mutable gsStaticBase<T> * m_solver;
    T m_L;
    bool first;

    gsVector<T> m_U;


};


} // namespace gismo
