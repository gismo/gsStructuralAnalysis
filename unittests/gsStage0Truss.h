/** @file gsStage0Truss.h

    @brief Two-dof von Mises truss used as a dependency-free reference problem
    for the gsStructuralAnalysis path-following solvers.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#pragma once

#include <gsStructuralAnalysis/src/gsStructuralAnalysisTools/gsStructuralAnalysisTypes.h>

namespace gismo
{

/**
   @brief Shallow symmetric two-bar truss (rise h, Green-Lagrange strain) with
   a linear spring k between the apex and the loaded degree of freedom.

   Unknowns are U = (w, v): w is the apex deflection (positive downward) and
   v the displacement of the loaded dof, which is connected to the apex by the
   spring. The load is F = (0, 1). With the strain energy EA L0 eps^2 and
   eps = (w^2 - 2 h w) / (2 L0^2), the bar force is
   \f[ f(w) = c\, w (w-h)(w-2h), \qquad c = EA/L_0^3, \f]
   and with the residual convention R = Fint - lambda F
   \f[ R(U,\lambda) = \begin{pmatrix} f(w) + k (w - v) \\ k (v - w) - \lambda \end{pmatrix}, \qquad
       K(U) = \begin{pmatrix} f'(w) + k & -k \\ -k & k \end{pmatrix}, \f]
   with \f$ f'(w) = c (3 w^2 - 6 h w + 2 h^2) \f$.

   The equilibrium path is known in closed form, \f$ \lambda(w) = f(w) \f$ and
   \f$ v(w) = w + f(w)/k \f$. The limit point is at
   \f$ w^* = h (1 - 1/\sqrt 3) \f$ with
   \f$ \lambda_{\max} = 2 c h^3 / (3 \sqrt 3) \f$. For \f$ k > c h^2 \ge \max(-f') \f$
   there is no snap-back in v: K is singular only at the limit point(s), and the
   first pivot f'(w) + k is positive everywhere. The path returns to
   \f$ \lambda = \lambda_{\max} \f$ on the stable far branch at
   \f$ w = h (1 + 2/\sqrt 3) \f$.

   Default parameters (h = 1, c = 3 sqrt(3)/2, k = 5) give lambda_max = 1.
*/
template<class T>
class gsStage0Truss
{
public:
    /// Constructor: bar stiffness \a c = EA/L0^3, rise \a h, spring stiffness \a k
    explicit gsStage0Truss(T c = T(1.5)*math::sqrt(T(3)), T h = 1, T k = 5)
    : m_c(c), m_h(h), m_k(k)
    { }

    /// Bar force f(w)
    T fint(T w) const { return m_c * w * (w - m_h) * (w - 2*m_h); }

    /// Tangent stiffness f'(w) of the bar
    T dfint(T w) const { return m_c * (3*w*w - 6*m_h*w + 2*m_h*m_h); }

    /// Load factor on the equilibrium path at apex deflection \a w
    T lambdaOf(T w) const { return fint(w); }

    /// Equilibrium point U(w) = (w, w + f(w)/k)
    gsVector<T> pathPoint(T w) const
    {
        gsVector<T> U(2);
        U << w, w + fint(w)/m_k;
        return U;
    }

    /// Apex deflection of the limit point
    T wLimit() const { return m_h * (1 - 1/math::sqrt(T(3))); }

    /// Load factor of the limit point
    T lambdaMax() const { return m_c * 2 * m_h*m_h*m_h / (3*math::sqrt(T(3))); }

    /// Reference load vector (0,1)
    gsVector<T> force() const
    {
        gsVector<T> F(2);
        F << 0, 1;
        return F;
    }

    /// Residual R = Fint(U) - lambda F; always assemblable
    bool residual(const gsVector<T>& U, T lambda, gsVector<T>& R) const
    {
        const T w = U[0], v = U[1];
        R.resize(2);
        R[0] = fint(w) + m_k*(w - v);
        R[1] = m_k*(v - w) - lambda;
        return true;
    }

    /// Tangent stiffness matrix K(U); resizes \a K to 2x2
    bool jacobianMatrix(const gsVector<T>& U, gsSparseMatrix<T>& K) const
    {
        K.resize(2,2);
        K.setZero();
        K.insert(0,0) = dfint(U[0]) + m_k;
        K.insert(0,1) = -m_k;
        K.insert(1,0) = -m_k;
        K.insert(1,1) = m_k;
        return true;
    }

    /// Residual functor; captures a copy of the truss, so it outlives it
    typename gsStructuralAnalysisOps<T>::ALResidual_t alResidual() const
    {
        const gsStage0Truss<T> self(*this);
        return [self](gsVector<T> const& U, const T L, gsVector<T>& R)
        { return self.residual(U, L, R); };
    }

    /// Jacobian functor; captures a copy of the truss, so it outlives it
    typename gsStructuralAnalysisOps<T>::Jacobian_t jacobian() const
    {
        const gsStage0Truss<T> self(*this);
        return [self](gsVector<T> const& U, gsSparseMatrix<T>& K)
        { return self.jacobianMatrix(U, K); };
    }

private:
    T m_c, m_h, m_k;
};

} // namespace gismo
