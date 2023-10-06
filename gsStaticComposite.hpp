 /** @file gsStaticComposite.hpp

    @brief Performs linear modal analysis given a matrix or functions of a matrix

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#include <typeinfo>
#pragma once

namespace gismo
{

template <class T>
gsStatus gsStaticComposite<T>::solve()
{
    this->getOptions();
    m_numIterations = 0;
    m_status = gsStatus::Success;
    for (size_t k=0; k!=m_solvers.size() && m_status==gsStatus::Success; k++)
    {
        if (m_verbose > 0)
            gsInfo<<"Solver "<<k<<" out of "<<m_solvers.size()<<"\n";
        
        if (k>0)
        {
            m_solvers[k]->setUpdate(m_solvers[k-1]->update());
        }

        m_status = m_solvers[k]->solve();
        m_numIterations += m_solvers[k]->iterations();
        m_U = m_solvers[k]->solution();
        m_DeltaU = m_solvers[k]->update();
    }
    return m_status;
};

template <class T>
void gsStaticComposite<T>::initialize()
{
    for (typename std::vector<gsStaticBase<T> *>::iterator solver = m_solvers.begin(); solver!=m_solvers.end(); solver++)
        (*solver)->initialize();
};

template <class T>
void gsStaticComposite<T>::defaultOptions()
{
    // Commented because it is assumed that options are handled in each individual solver
    // for (typename std::vector<gsStaticBase<T> *>::iterator solver = m_solvers.begin(); solver!=m_solvers.end(); solver++)
    //     (*solver)->defaultOptions();

    m_options.addInt("verbose","Verbose output",0);
};

template <class T>
void gsStaticComposite<T>::getOptions()
{
    // Commented because it is assumed that options are handled in each individual solver
    // for (typename std::vector<gsStaticBase<T> *>::iterator solver = m_solvers.begin(); solver!=m_solvers.end(); solver++)
    //     (*solver)->getOptions();

    m_verbose = m_options.getInt("verbose");
};

template <class T>
void gsStaticComposite<T>::setOptions(gsOptionList & options)
{
    gsWarn<<"setOptions cannot be used on a gsStaticComposite solver. Call setOptions on each solver individually\n";
};

template <class T>
void gsStaticComposite<T>::setDisplacement(const gsVector<T> & displacement)
{
    for (typename std::vector<gsStaticBase<T> *>::iterator solver = m_solvers.begin(); solver!=m_solvers.end(); solver++)
        (*solver)->setDisplacement(displacement);
};

template <class T>
void gsStaticComposite<T>::setLoad(const T L)
{
    for (typename std::vector<gsStaticBase<T> *>::iterator solver = m_solvers.begin(); solver!=m_solvers.end(); solver++)
        (*solver)->setLoad(L);
};

template <class T>
void gsStaticComposite<T>::setSolution(const gsVector<T> & displacement, const T L)
{
    for (typename std::vector<gsStaticBase<T> *>::iterator solver = m_solvers.begin(); solver!=m_solvers.end(); solver++)
        (*solver)->setSolution(displacement,L);
};

template <class T>
void gsStaticComposite<T>::setUpdate(const gsVector<T> & update)
{
    for (typename std::vector<gsStaticBase<T> *>::iterator solver = m_solvers.begin(); solver!=m_solvers.end(); solver++)
        (*solver)->setUpdate(update);
};

template <class T>
void gsStaticComposite<T>::reset()
{
    for (typename std::vector<gsStaticBase<T> *>::iterator solver = m_solvers.begin(); solver!=m_solvers.end(); solver++)
        (*solver)->reset();
};


} // namespace gismo
