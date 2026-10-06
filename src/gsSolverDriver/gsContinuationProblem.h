 /** @file gsContinuationProblem.h

    @brief The problem definition of the solver driver: residual, tangent
    stiffness and optional load direction as closures.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#pragma once
#include <functional>

#include <gsCore/gsLinearAlgebra.h>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsSolverDriverTypes.h>

namespace gismo
{

/**
   @brief Closures that define a nonlinear problem
   \f$R(U,\lambda) = \lambda F_{ext} - F_{int}(U) = 0\f$ for the solver driver.

   Every op is evaluated at the client's current working state, which is
   reached only through the gsStateSpace; no \f$U\f$ argument exists. Every op
   returns false on failure instead of throwing; the driver maps a false to a
   gsSolverReport with status AssemblyError and failedOp Residual, Jacobian or
   LoadDirection.

   Sign convention: the residual is \f$R = \lambda F - F_{int}\f$ and the
   Jacobian is the tangent stiffness \f$K = \partial F_{int}/\partial U =
   -\partial R/\partial U\f$, so a Newton correction solves
   \f$K\,\delta U = R\f$. This is the sign of the Jacobian_t of the gsALM* and
   gsStaticNewton classes. The gsALMBase facade negates the user's
   ALResidual_t (\f$F_{int} - \lambda F_{ext}\f$); gsStaticNewton and
   dynamic-relaxation facades pass theirs (\f$\lambda F_{ext} - F_{int}\f$)
   through unchanged.

   @tparam T coefficient type

   \ingroup gsStructuralAnalysis
*/
template<class T>
struct gsContinuationProblem
{
    /// Writes \f$R = \lambda F_{ext} - F_{int}(\text{working})\f$ of tangent
    /// size into @a R. Applies the Dirichlet ramp for @a lambda itself (there
    /// is no setLoadFactor). @a mode is passed explicitly by the driver
    /// (std::function cannot carry a default argument); with
    /// gsEvalMode::Probe the op must not record state for a later commit.
    typedef std::function<bool(T lambda, gsVector<T>& R, gsEvalMode mode)> Residual_t;

    /// Writes the tangent stiffness \f$K = \partial F_{int}/\partial U =
    /// -\partial R/\partial U\f$, n x n, linearised consistently with
    /// gsStateSpace::apply, so that a Newton correction solves
    /// \f$K\,\delta = R\f$.
    typedef std::function<bool(T lambda, gsSparseMatrix<T>& K)> Jacobian_t;

    /// Writes \f$\partial R/\partial\lambda\f$ of tangent size into @a dRdl.
    typedef std::function<bool(T lambda, gsVector<T>& dRdl)> LoadDirection_t;

    Residual_t      residual;       ///< Required
    Jacobian_t      jacobian;       ///< Required
    LoadDirection_t loadDirection;  ///< Optional, for follower loads or lambda-dependent Dirichlet data; empty means \f$\partial R/\partial\lambda = F\f$
    gsVector<T>     F;              ///< Reference external force, tangent size
};

} // namespace gismo
