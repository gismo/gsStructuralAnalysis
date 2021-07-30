/** @file .cpp

    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#include <gismo.h>
#include <gsStructuralAnalysis/gsStaticBase.h>

namespace gismo
{

template <class T>
class gsStaticDR : public gsStaticBase<T>
{
protected:

    typedef gsStaticBase<T> Base;

    typedef std::function<gsVector<T>(gsVector<T> const &)   >  Residual_t;
    typedef std::function<gsVector<T>(gsVector<T> const &, T)>  ALResidual_t;

public:
    /// Constructor given the matrices and timestep
    gsStaticDR( const gsVector<T> & M, // lumped
                const gsVector<T> & F,
                const Residual_t &Residual
               )
    :
    m_mass(M),
    m_forcing(F),
    m_residualFun(Residual),
    m_ALresidualFun(nullptr)
    {
        this->_init();
    }

    /// Constructor given the matrices and timestep
    gsStaticDR( const gsVector<T>   & M, // lumped
                const gsVector<T>   & F,
                const ALResidual_t  & ALResidual
               )
    :
    m_mass(M),
    m_forcing(F),
    m_ALresidualFun(ALResidual)
    {
        m_L = 1.0;
        m_residualFun = [this](gsVector<T> const & x)
        {
            return m_ALresidualFun(x,m_L);
        };
        this->_init();
    }

/// gsStaticBase base functions
public:
    void solve();
    void initialize();

    void initOutput();
    void stepOutput(index_t k);

    void defaultOptions();
    void getOptions();

    void setOptions(gsOptionList & options) {m_options.update(options,gsOptionList::addIfUnknown); }

/// Class-specific
public:
    T kineticEnergy() const { return m_Ek; }
    gsVector<T> energies() const { return gsAsVector<T>(m_Eks); }
    gsVector<T> relEnergies() const { return gsAsVector<T>(m_Eks)/m_Ek0; }
    gsVector<T> velocities() const {return m_v;}

protected:

    void _iteration();
    void _peak();
    void _start();
    void _init();

public:
    //// ????
    void _stepBack()
    {
        m_U -= m_DeltaU;
    }

    //// ????
    void _reset()
    {
        _start();
    }




    T residualNorm() const { return m_R.norm(); }



protected:
    const gsVector<T> m_mass;
    const gsVector<T> m_forcing;
    Residual_t m_residualFun;
    const ALResidual_t m_ALresidualFun;

    // Solution
    using Base::m_U;
    using Base::m_DeltaU;
    using Base::m_deltaU;

    using Base::m_L;
    using Base::m_DeltaL;
    using Base::m_deltaL;


    // Iterations
    using Base::m_numIterations;
    using Base::m_maxIterations;

    // Residuals
    using Base::m_R;

    // Tolerances
    using Base::m_converged;
    using Base::m_tolF;
    using Base::m_tolU;
    T m_tolE;

    // Residual norms
    using Base::m_residual;
    using Base::m_residualIni;
    using Base::m_residualOld;

    // Options
    using Base::m_options;
    using Base::m_verbose;

    // DoFs
    using Base::m_dofs;

    // Headstart
    using Base::m_headstart;



    // Velocities
    gsVector<T> m_v;                // (class-specific)
    gsVector<T> m_massInv, m_damp;  // (class-specific)



    // Coefficients
    T m_dt, m_alpha, m_c;           // (class-specific)

    // Kinetic energy
    T m_Ek, m_Ek_prev, m_Ek0;       // (class-specific)
    mutable std::vector<T> m_Eks;   // (class-specific)
};

} //namespace


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsStaticDR.hpp)
#endif