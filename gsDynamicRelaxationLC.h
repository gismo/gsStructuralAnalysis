/** @file .cpp

    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#include <gismo.h>

namespace gismo
{

template <class T>
class gsDynamicRelaxationLC
{
public:
    /// Constructor given the matrices and timestep
    gsDynamicRelaxationLC(  const gsVector<T> & M, // lumped
                            const gsVector<T> & F,
                            const std::function < gsVector<T> ( gsVector<T> const & , T, gsVector<T> const &) > &Residual
                        )
    :
    m_mass(M),
    m_forcing(F),
    m_residualFun(Residual)
    {
        m_dt = 1.0;
        m_dofs = m_mass.rows();
        m_v = m_U = m_DeltaU = m_deltaU = gsVector<T>(m_dofs);
        m_v.setZero();
        m_U.setZero();
        m_DeltaU.setZero();
        m_deltaU.setZero();
        m_L = m_DeltaL = 0.0;

        defaultOptions();
        m_massInv = M.array().inverse();
    }

    gsVector<T> displacements() const {return m_U;}
    gsVector<T> increment() const {return m_DeltaU;}
    T load() const {return m_L;}
    gsVector<T> velocities() const {return m_v;}

    void setDisplacement(gsVector<T> displacement) { m_U = displacement; }

    T kineticEnergy() const { return m_Ek; }

    index_t iterations() const { return m_iterations; }

    gsVector<T> energies() const { return gsAsVector<T>(m_Eks); }
    gsVector<T> relEnergies() const { return gsAsVector<T>(m_Eks)/m_Ek0; }

    void init()
    {
        getOptions();
        m_massInv *= 1./m_alpha;
        m_damp = m_c * m_mass;
    }

    void predictor(T loadStep)
    {
        m_DeltaL = loadStep;
        predictor();
    }

    void predictor()
    {
        start();
        m_Ek = m_v.transpose() * m_mass.cwiseProduct(m_v);
    }

    void iteration(T loadStep)
    {
        m_DeltaL = loadStep;
        iteration();
    }

    void iteration()
    {
        m_Ek_prev = m_Ek;

        m_R = m_residualFun(m_U+m_DeltaU,m_L+m_DeltaL,m_forcing) - m_damp.cwiseProduct(m_v);
//----------------------------------------------------------------------------------
        m_v += m_dt * m_massInv.cwiseProduct(m_R);                    // Velocities at t+dt/2
        m_deltaU = m_dt * m_v;               // Velocities at t+dt
//----------------------------------------------------------------------------------
        // m_deltaU = m_deltaU + m_dt*m_dt*m_massInv*m_R;
        // m_v = m_deltaU/m_dt;
//----------------------------------------------------------------------------------

        m_DeltaU += m_deltaU;               // Velocities at t+dt
        m_Ek = m_v.transpose() * m_mass.cwiseProduct(m_v);
    }

    void step(T loadStep)
    {
        m_Eks.clear();
        m_Eks.reserve(m_maxIt);

        if (m_verbose) printHeader();
        m_DeltaL = loadStep;
        predictor();
        m_Ek0 = m_Ek;
        m_Eks.push_back(m_Ek);

        m_R0 = m_R;
        if (m_verbose) stepInfo(0);
        for (m_iterations=1; m_iterations!=m_maxIt; m_iterations++)
        {
            iteration();
            if ((m_c==0 && m_Ek_prev > m_Ek))// || (m_Ek/m_Ek_prev > 1/m_tolE && m_Ek_prev!=0))
                peak();

            if (m_verbose) stepInfo(m_iterations);
            m_Eks.push_back(m_Ek);

            if (m_R.norm()/m_R0.norm() < m_tolF && m_Ek/m_Ek0 < m_tolE)
            {
                m_L += m_DeltaL;
                m_U += m_DeltaU;
                break;
            }
            if (m_iterations==m_maxIt-1)
            {
                m_L += m_DeltaL;
                m_U += m_DeltaU;
                gsWarn<<"Maximum iterations reached!\n";
            }
        }
        gsInfo<<"\n";
    }

    void stepBack()
    {
        m_L -= m_DeltaL;
        m_U -= m_DeltaU;
    }

    void reset()
    {
        start();
    }

    void peak(T load)
    {
        m_L = load;
        peak();
    }

    void peak()
    {
        m_R = m_residualFun(m_U+m_DeltaU,m_L+m_DeltaL,m_forcing)- m_damp.cwiseProduct(m_v);
        m_deltaU = - 1.5 * m_dt * m_v + m_dt*m_dt / 2. * m_massInv.cwiseProduct(m_R);
        m_DeltaU += m_deltaU;

        m_v.setZero();
        m_Ek = 0.0;
        // m_v = 0.5 * m_dt * m_massInv * m_R; // Velocities at dt/2
    }

    void start()
    {
        m_DeltaU.setZero();
        m_v.setZero();

//----------------------------------------------------------------------------------
        m_R = m_residualFun(m_U+m_DeltaU,m_L+m_DeltaL,m_forcing);
        m_deltaU = 0.5*m_dt*m_dt*m_massInv.cwiseProduct(m_R);
        m_v = m_deltaU/m_dt;
//----------------------------------------------------------------------------------
        // m_R = m_residualFun(m_U+m_DeltaU,load,m_forcing)- m_damp*m_v;
        // m_v = 0.5 * m_dt * m_massInv * m_R; // Velocities at dt/2
        // m_deltaU = m_dt * m_v;               // Velocities at t+dt
//----------------------------------------------------------------------------------

        m_DeltaU += m_deltaU;
    }

    void defaultOptions()
    {
        m_options.addReal("damping","damping factor",1.0);
        m_options.addReal("alpha","mass coefficient",2.0);
        m_options.addInt("maxIt","maximum number of iterations",1e2);
        m_options.addReal("tolF","(Force) Residual tolerance",1e-6);
        m_options.addReal("tolE","Kinetic energy tolerance",1e-12);
        m_options.addInt("verbose","Verbose output",0);
    }

    gsOptionList options() const {return m_options;}

    void getOptions()
    {
        m_c = m_options.getReal("damping");
        m_alpha = m_options.getReal("alpha");
        m_maxIt = m_options.getInt("maxIt");
        m_tolF = m_options.getReal("tolF");
        m_tolE = m_options.getReal("tolE");
        m_verbose = m_options.getInt("verbose");
    }

    void setOptions(gsOptionList & options) {m_options.update(options,gsOptionList::addIfUnknown); }

    T residualNorm() const { return m_R.norm(); }

    void printHeader()
    {
        gsInfo  <<std::setw(8)<<"It."
                <<std::setw(16)<<"|R|"
                <<std::setw(16)<<"Ek"
                <<std::setw(16)<<"|u|"
                <<std::setw(16)<<"|Δu|"
                <<std::setw(16)<<"|δu|"
                <<std::setw(16)<<"λ"
                <<std::setw(16)<<"Δλ"
                <<std::setw(16)<<"ΔL"
                <<std::setw(16)<<"|Δv|"
                <<"\n";
    }

    void stepInfo(index_t k)
    {
        gsInfo  <<std::setw(8)<<k
                <<std::setw(16)<<m_R.norm()/m_R0.norm()
                <<std::setw(16)<<m_Ek/m_Ek0
                <<std::setw(16)<<(m_U+m_DeltaU).norm()
                <<std::setw(16)<<m_DeltaU.norm()
                <<std::setw(16)<<m_deltaU.norm()
                <<std::setw(16)<<m_L+m_DeltaL
                <<std::setw(16)<<m_DeltaL
                <<std::setw(16)<<m_DeltaU.norm()+m_DeltaL*m_DeltaL
                <<std::setw(16)<<m_v.norm()
                <<"\n";
    }

protected:
    const gsVector<T> m_mass;
    const gsVector<T> m_forcing;
    const std::function<gsVector<T> ( gsVector<T> const &, T, gsVector<T> const &) > m_residualFun;
    gsVector<T> m_U, m_v;
    gsVector<T> m_DeltaU, m_deltaU;
    gsVector<T> m_R, m_R0;
    gsVector<T> m_massInv, m_damp;
    index_t m_dofs;

    bool m_verbose;

    T m_dt, m_alpha, m_c;
    gsOptionList m_options;
    T m_Ek, m_Ek_prev, m_Ek0;
    T m_L, m_DeltaL;

    index_t m_maxIt, m_iterations;
    T m_tolF, m_tolE;

    mutable std::vector<T> m_Eks;
};

} //namespace