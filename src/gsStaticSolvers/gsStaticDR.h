/** @file gsStaticDR.cpp

    @brief Static solver using the Dynamic Relaxation method

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#include <gsStructuralAnalysis/src/gsStaticSolvers/gsStaticBase.h>

namespace gismo
{

/**
 * @brief Static solver using the Dynamic Relaxation method
 *
 * @tparam     T     coefficient type
 *
 * \ingroup gsStaticBase
 */
template <class T>
class gsStaticDR : public gsStaticBase<T>
{
protected:

    typedef gsStaticBase<T> Base;

    typedef typename Base::Residual_t    Residual_t;
    typedef typename Base::ALResidual_t  ALResidual_t;

public:

    /**
     * @brief      Constructor
     *
     * @param[in]  M         Lumped mass matrix (as vector)
     * @param[in]  F         External forcing vector
     * @param[in]  Residual  The residual
     */
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


    /**
     * @brief      Constructs a new instance.
     *
     * @param[in]  M         Lumped mass matrix (as vector)
     * @param[in]  F         External forcing vector
     * @param[in]  Residual  The residual as arc-length object
     */
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
        m_residualFun = [this](gsVector<T> const & x, gsVector<T> & result) -> bool
        {
            return m_ALresidualFun(x,m_L,result);
        };
        this->_init();
    }

/// gsStaticBase base functions
public:
    /// See \ref gsStaticBase
    gsStatus solve() override;

    /// See \ref gsStaticBase
    void initialize() override;

    /// See \ref gsStaticBase
    void initOutput() override;
    
    /// See \ref gsStaticBase
    void stepOutput(index_t k) override;

    /// See \ref gsStaticBase
    void defaultOptions() override;

    /// See \ref gsStaticBase
    void reset() override;

    /// See \ref gsStaticBase
    void getOptions() override;

public:
    /// Returns the kinetic energy
    T kineticEnergy() const { return m_Ek; }
    /// Returns the kinetic energy in all steps
    gsVector<T> energies() const { return gsAsVector<T>(m_Eks); }
    /// Returns the kinetic energy relative to the first iteration
    gsVector<T> relEnergies() const { return gsAsVector<T>(m_Eks)/m_Ek0; }
    /// Returns the velocity
    gsVector<T> velocities() const {return m_v;}

protected:
    /// See \ref solve()
    void _solve();
    /// Performs an iteration
    void _iteration();
    /// Identifies a peak
    void _peak();
    /// Starts the method
    void _start();
    /// Initializes the method
    void _init();

    gsVector<T> _computeResidual(const gsVector<T> & U);

public:
    //// Perform a step back
    void _stepBack()
    {
        m_U -= m_DeltaU;
    }

    /// Start over again
    void _reset()
    {
        _start();
    }

    /// Return the residual norm
    T residualNorm() const { return m_R.norm(); }

protected:
    const gsVector<T> & m_mass;
    const gsVector<T> & m_forcing;
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
    index_t m_resetIterations;

    // Residuals
    using Base::m_R;

    // Tolerances
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

    // Solver status
    using Base::m_status;

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
