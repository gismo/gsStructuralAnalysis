 /** @file gsStructuralAnalysisTypes.h

    @brief Provides a status object and typedefs

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

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

    template<class T>
    struct gsStructuralAnalysisOps
    {
        /// Force
        typedef std::function < bool (           gsVector<T> & )>     Force_t;
        /// Time-dependent force
        typedef std::function < bool ( const T,  gsVector<T> & )>     TForce_t;

        /// Residual, Fint-Fext
        typedef std::function < bool ( gsVector<T> const &,     gsVector<T> & )>    Residual_t;
        /// Arc-Length Residual, Fint-lambda*Fext
        typedef std::function < bool ( gsVector<T> const &, const T,  gsVector<T> & )>    ALResidual_t;
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
} // namespace gismo
