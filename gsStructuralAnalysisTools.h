 /** @file gsStructuralAnalysisTools.h

    @brief Provides a status object and typedefs

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#pragma once

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
        typedef std::function < bool ( T,  gsVector<T> & )>    TForce_t;

        // Residual, Fint-Fext
        typedef std::function < bool ( gsVector<T> const &,     gsVector<T> & )>    Residual_t;
        // Arc-Length Residual, Fint-lambda*Fext
        typedef std::function < bool ( gsVector<T> const &, T,  gsVector<T> & )>    ALResidual_t;
        // Time-dependent Residual Fint(t)-Fext(t)
        typedef std::function < bool ( gsVector<T> const &, T,  gsVector<T> & )>    TResidual_t;
        
        // Jacobian
        typedef std::function < bool ( gsVector<T> const &,                         gsSparseMatrix<T> & ) > Jacobian_t;
        // Jacobian with solution update as argument
        typedef std::function < bool ( gsVector<T> const &, gsVector<T> const &,    gsSparseMatrix<T> & ) > dJacobian_t;
    };
} // namespace gismo
