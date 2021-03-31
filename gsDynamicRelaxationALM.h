/** @file .cpp

    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#include <gismo.h>
using namespace gismo;

template <class T>
class gsDynamicRelaxationALM
{
public:
    /// Constructor given the matrices and timestep
    gsDynamicRelaxationALM( const gsSparseMatrix<T> & K,
                            const gsVector<T> & F,
                            const std::function < gsVector<T> ( gsVector<T> const & , T, gsVector<T> const &) > &Residual
                        )
    :
    m_stif(K),
    m_forcing(F),
    m_residualFun(Residual)
    {
        m_dt = 1.0;
        m_dofs = m_stif.rows();
        m_V = m_U = m_DeltaU = m_deltaU = m_DeltaUold = gsVector<T>(m_dofs);
        m_V.setZero();
        m_U.setZero();
        m_DeltaU.setZero();
        m_deltaU.setZero();
        m_DeltaUold.setZero();

        m_L = m_DeltaL = m_deltaL = m_DeltaLold = 0.0;
        m_phi = 0.0;

        defaultOptions();
    }

    int sign(T val)
    {
        return (T(0) < val) - (val < T(0));
    }

    gsVector<T> displacements() const {return m_U+m_DeltaU;}
    gsVector<T> velocities() const {return m_V;}

    void setDisplacement(gsVector<T> displacement) { m_U = displacement; }

    T kineticEnergy() const { return m_Ek; }

    index_t iterations() const { return m_iterations; }

    gsVector<T> energies() const { return gsAsVector<T>(m_Eks); }
    gsVector<T> relEnergies() const { return gsAsVector<T>(m_Eks)/m_Ek0; }

    void init()
    {
        getOptions();
        m_mass = m_alpha * m_dt*m_dt / 2. * m_stif;
        m_massInv = m_mass.toDense().inverse();
        m_damp = m_c * m_mass;
    }

    void step(T length)
    {
        m_Eks.clear();
        m_Eks.reserve(m_maxIt);

        printHeader();
        m_arcLength = length;
        this->predictor();
        m_Ek0 = m_Ek;
        m_Eks.push_back(m_Ek);
        stepInfo(0);
        for (index_t m_iterations=1; m_iterations!=m_maxIt; m_iterations++)
        {
            iteration();
            if (m_c==0 && m_Ek_prev > m_Ek)
                peak();

        	m_Eks.push_back(m_Ek);

            stepInfo(m_iterations);

            if (m_R.norm()/m_forcing.norm() < m_tolF && m_Ek/m_Ek0 < m_tolE)
            {
                finish();
                break;
            }
            if (m_iterations==m_maxIt-1)
            {
                finish();
                gsWarn<<"Maximum iterations reached!\n";
            }
        }
        gsInfo<<"\n";
    }

    void predictor()
    {
        // initiateStep
        m_V = m_DeltaU = m_deltaU = m_deltaUbar = m_deltaUt = gsVector<T>::Zero(m_dofs);
        m_DeltaL = m_deltaL = 0.0;

        // compute deltaUt
        m_R = m_forcing;
        m_V = 0.5*m_dt*m_massInv*m_R;
        m_Ek = m_V.transpose() * m_mass * m_V;
        m_deltaUt = m_V/m_dt;

        // Predictor
        if (m_DeltaUold.dot(m_DeltaUold) == 0 && m_DeltaLold*m_DeltaLold == 0) // no information about previous step.
        {
            T DL = 1.;
            m_deltaL = m_arcLength * DL / math::sqrt( ( m_deltaUt.dot(m_deltaUt) + m_DeltaL*DL ) );

            m_deltaU = m_deltaL*m_deltaUt;

            // m_phi = math::pow( m_deltaUt.dot(m_deltaUt) / m_forcing.dot(m_forcing),0.5);

            m_DeltaUold = m_deltaU;
            m_DeltaLold = m_deltaL;
        }
        else
        {
            m_phi = math::pow(m_U.dot(m_U)/( math::pow(m_L,2) * m_forcing.dot(m_forcing) ),0.5);
            m_phi = 0.0;
            T A0 = math::pow(m_phi,2)* m_forcing.dot(m_forcing); // see Lam et al. 1991
            int dir = sign(m_DeltaUold.dot(m_deltaUt) + A0*m_DeltaLold); // Feng et al. 1995 with H = \Psi^2
            T denum = ( math::pow( m_deltaUt.dot(m_deltaUt) + A0 ,0.5) ); // Feng et al. 1995 with H = \Psi^2

            T mu;
            if (denum==0)
              mu = m_arcLength;
            else
              mu = m_arcLength / denum;

            m_deltaL = dir*mu;
            m_deltaU = m_deltaL*m_deltaUt;
        }

        m_DeltaU += m_deltaU;
        m_DeltaL += m_deltaL;
    }

    void iteration()
    {
        m_Ek_prev = m_Ek;

        m_R = m_residualFun(m_U+m_DeltaU,m_L+m_DeltaL,m_forcing) - m_damp*m_V;
        m_V += m_dt * m_massInv * m_R;                   // Velocities at t+dt/2

        m_Ek = m_V.transpose() * m_mass * m_V;

        m_deltaUbar = m_dt*m_V;               // Velocities at t+dt

        T A0 = math::pow(m_phi,2)* m_forcing.dot(m_forcing); // see Lam et al. 1991,

        T a0 = m_deltaUt.dot(m_deltaUt) + A0;
        T b0 = 2*( m_deltaUt.dot(m_DeltaU) + m_DeltaL * A0 );
        T b1 = 2*( m_deltaUbar.dot(m_deltaUt) );
        T c0 = m_DeltaU.dot(m_DeltaU) + m_DeltaL*m_DeltaL * A0 - math::pow(m_arcLength,2);
        T c1 = 2*( m_DeltaU.dot(m_deltaUbar) );
        T c2 = m_deltaUbar.dot(m_deltaUbar);

        /// Calculate the coefficients of the polynomial
        T alpha1 = a0;
        T alpha2 = b0 + b1;
        T alpha3 = c0 + c1 + c2;

        T discriminant = math::pow(alpha2 ,2) - 4 * alpha1 * alpha3;

        gsVector<T> deltaLs(2);
        deltaLs.setZero();
        if (discriminant >= 0)
        {
          deltaLs[0] = (-alpha2 + math::pow(discriminant,0.5))/(2*alpha1);
          deltaLs[1] = (-alpha2 - math::pow(discriminant,0.5))/(2*alpha1);

          gsVector<T> deltaU1, deltaU2;
          deltaU1 = m_deltaUbar + m_deltaUt*deltaLs[0];
          deltaU2 = m_deltaUbar + m_deltaUt*deltaLs[1];

          // ---------------------------------------------------------------------------------
          // Method by Ritto-Corea et al. 2008
          T DOT1,DOT2;
          DOT1 = deltaLs[0]*(m_DeltaUold.dot(m_deltaUt) + math::pow(m_phi,2)*m_DeltaLold);
          DOT2 = deltaLs[1]*(m_DeltaUold.dot(m_deltaUt) + math::pow(m_phi,2)*m_DeltaLold);

          if (DOT1 > DOT2)
          {
            m_deltaL = deltaLs[0];
            m_deltaU = deltaU1;
          }
          else if (DOT1 < DOT2)
          {
            m_deltaL = deltaLs[1];
            m_deltaU = deltaU2;
          }
          else
          {
            m_deltaL = deltaLs[0];
            m_deltaU = deltaU1;
          }
        }
        else
            GISMO_ERROR("Discriminant is negative!");


        m_DeltaU += m_deltaU;
        m_DeltaL += m_deltaL;
    }

    // void reset()
    // {
    //     start();
    // }

    void peak()
    {
        m_R = m_residualFun(m_U+m_DeltaU,m_L+m_DeltaL,m_forcing)- m_damp*m_V;
        m_deltaU = - 1.5 * m_dt * m_V + m_dt*m_dt / 2. * m_massInv * m_R;
        m_DeltaU += m_deltaU;

        m_V.setZero();
        // m_V = 0.5 * m_dt * m_massInv * m_R; // Velocities at dt/2
    }

    void finish()
    {
        m_U += m_DeltaU;
        m_L += m_DeltaL;
        m_DeltaUold = m_DeltaU;
        m_DeltaLold = m_DeltaL;
    }

    void defaultOptions()
    {
        m_options.addReal("damping","damping factor",1.0);
        m_options.addReal("alpha","mass coefficient",2.0);
        m_options.addReal("length","arc-length",1.0);
        m_options.addInt("maxIt","maximum number of iterations",1e2);
        m_options.addReal("tolF","(Force) Residual tolerance",1e-6);
        m_options.addReal("tolE","Kinetic energy tolerance",1e-12);
    }

    gsOptionList options() const {return m_options;}

    void getOptions()
    {
        m_c = m_options.getReal("damping");
        m_alpha = m_options.getReal("alpha");
        m_arcLength = m_options.getReal("length");
        m_maxIt = m_options.getInt("maxIt");
        m_tolF = m_options.getReal("tolF");
        m_tolE = m_options.getReal("tolE");
    }

    void setOptions(gsOptionList & options) {m_options.update(options,gsOptionList::addIfUnknown); }

    T residualNorm() const { return m_R.norm(); }

    void printHeader()
    {
        gsInfo  <<std::setw(4)<<std::left<<"It."
                <<std::setw(16)<<std::left<<"|R|"
                <<std::setw(16)<<std::left<<"Ek"
                <<std::setw(16)<<std::left<<"|u|"
                <<std::setw(16)<<std::left<<"|Du|" // Δu
                <<std::setw(16)<<std::left<<"|du|" // δu
                <<std::setw(16)<<std::left<<"L" // λ
                <<std::setw(16)<<std::left<<"DL" // Δλ
                <<std::setw(16)<<std::left<<"dL" // δλ
                <<std::setw(16)<<std::left<<"Dl"
                <<std::setw(16)<<std::left<<"Dlu"
                <<std::setw(16)<<std::left<<"Dlλ"
                <<"\n";
    }

    void stepInfo(index_t k)
    {
        T A0 = math::pow(m_phi,2)*m_forcing.dot(m_forcing);

        gsInfo  <<std::setw(4)<<std::left<<k
                <<std::setw(16)<<std::left<<m_R.norm()/m_forcing.norm()
                <<std::setw(16)<<std::left<<m_Ek/m_Ek0
                <<std::setw(16)<<std::left<<(m_U+m_DeltaU).norm()
                <<std::setw(16)<<std::left<<m_DeltaU.norm()
                <<std::setw(16)<<std::left<<m_deltaU.norm()
                <<std::setw(16)<<std::left<<m_L+m_DeltaL
                <<std::setw(16)<<std::left<<m_DeltaL
                <<std::setw(16)<<std::left<<m_deltaL
                <<std::setw(16)<<std::left<<m_arcLength //math::pow(m_DeltaU.dot(m_DeltaU) + A0*math::pow(m_DeltaL,2.0),0.5);
                <<std::setw(16)<<std::left<<math::pow(m_DeltaU.norm(),2.0)
                <<std::setw(16)<<std::left<<A0*math::pow(m_DeltaL,2.0)
                <<"\n";
    }

protected:
    const gsSparseMatrix<T> m_stif;
    const gsVector<T> m_forcing;
    const std::function<gsVector<T> ( gsVector<T> const &, T, gsVector<T> const &) > m_residualFun;
    gsVector<T> m_U, m_V;
    gsVector<T> m_DeltaU, m_deltaU, m_deltaUt, m_deltaUbar;
    gsVector<T> m_DeltaUold;
    T m_DeltaLold;
    gsVector<T> m_R;
    gsSparseMatrix<T> m_mass, m_damp;
    gsMatrix<T> m_massInv;
    index_t m_dofs;
    T m_dt, m_alpha, m_c;
    gsOptionList m_options;
    T m_Ek, m_Ek_prev, m_Ek0;
    T m_L, m_DeltaL, m_deltaL;
    T m_arcLength;
    T m_phi;
    index_t m_maxIt, m_iterations;
    T m_tolF, m_tolE;

    std::vector<T> m_Eks;
};