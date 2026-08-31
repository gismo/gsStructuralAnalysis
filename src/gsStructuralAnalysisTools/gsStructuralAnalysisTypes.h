 /** @file gsStructuralAnalysisTypes.h

    @brief Provides a status object and typedefs

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#pragma once

#include <functional>
#include <gsCore/gsLinearAlgebra.h>

namespace gismo
{

enum struct gsStatus
{
    Success,         ///< Successful
    NotConverged,    ///< Step did not converge
    AssemblyError,   ///< Assembly problem in step
    SolverError,     ///< Assembly problem in step
    NotStarted,      ///< ALM has not started yet
    OtherError       ///< Other error
};

// ALTERNATIVE IMPLEMENTATION USING FUNCTORS
// template<typename T>
// struct DynamicFunctor
// {
//     Force();

//     gsVector<T>& operator(T time = 0) {};

//     virtual index_t rows() = 0;
//     virtual index_t cols() = 0;
// };


// template<typename T>
// struct DynamicForce : public DynamicFunctor<T>
// {
//     DynamicForce(... ) {}


// };

/**
 * @brief      Operators for the gsStructuralAnalysis module
 *
 * @tparam     T     Double type
 */
template<class T>
struct gsStructuralAnalysisOps
{
    /// Residual energy
    typedef std::function < bool ( gsVector<T> const &, T &)>     Energy_t;

    /// Force
    typedef std::function < bool (           gsVector<T> & )>     Force_t;
    /// Time-dependent force
    typedef std::function < bool ( const T,  gsVector<T> & )>     TForce_t;

    /// Residual, Fint-Fext
    typedef std::function < bool ( gsVector<T> const &,     gsVector<T> & )>    Residual_t;
    /// Arc-Length Residual, Fint-lambda*Fext
    typedef std::function < bool ( gsVector<T> const &, const T,  gsVector<T> & )>    ALResidual_t;
    /// Arc-Length Force, i.e. the load derivative \f$ f(U,\Lambda) = -\partial R/\partial\Lambda \f$
    /// evaluated at the state \f$(U,\Lambda)\f$. For a dead load this is the constant
    /// external force vector \f$F_{ext}\f$; for state-dependent (follower, configuration-
    /// dependent) loads it varies with \f$U\f$ and/or \f$\Lambda\f$. See
    /// gsALMBase::setForcingFunction for the role this plays in the arc-length correctors.
    typedef std::function < bool ( gsVector<T> const &, const T,  gsVector<T> & )>    ALForce_t;
    /// Time-dependent Residual Fint(t)-Fext(t)
    typedef std::function < bool ( gsVector<T> const &, const T,  gsVector<T> & )>    TResidual_t;
    
    /// Mass matrix
    typedef std::function < bool (                                              gsSparseMatrix<T> & ) > Mass_t;
    /// Time-dependent mass matrix
    typedef std::function < bool (                      const T,                gsSparseMatrix<T> & ) > TMass_t;
    /// Damping matrix
    typedef std::function < bool ( gsVector<T> const &,                         gsSparseMatrix<T> & ) > Damping_t;
    /// Time-dependent Damping matrix
    typedef std::function < bool ( gsVector<T> const &, const T,                gsSparseMatrix<T> & ) > TDamping_t;

    /// Stiffness matrix
    typedef std::function < bool (                                              gsSparseMatrix<T> & ) > Stiffness_t;
    /// Jacobian
    typedef std::function < bool ( gsVector<T> const &,                         gsSparseMatrix<T> & ) > Jacobian_t;
    /// Jacobian
    typedef std::function < bool ( gsVector<T> const &, const T,                gsSparseMatrix<T> & ) > TJacobian_t;
    /// Jacobian with solution update as argument
    typedef std::function < bool ( gsVector<T> const &, gsVector<T> const &,    gsSparseMatrix<T> & ) > dJacobian_t;
};

}
