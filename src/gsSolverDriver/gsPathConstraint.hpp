 /** @file gsPathConstraint.hpp

    @brief Implementation of the path constraints.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#pragma once
#include <gsStructuralAnalysis/src/gsSolverDriver/gsPathConstraint.h>

namespace gismo
{

// ---------------------------------------------------------------- Fixed

template<class T>
bool gsFixedConstraint<T>::predictorNeedsTangent(const gsStepContext<T> &) const
{
    return false;
}

template<class T>
bool gsFixedConstraint<T>::predict(gsStepContext<T> & ctx, const gsVector<T> &,
                                   const gsStateSpace<T> & space)
{
    const index_t n = space.tangentSize();
    ctx.DeltaU.setZero(n);
    ctx.deltaU.setZero(n);
    ctx.DeltaL = ctx.deltaL = 0;
    return true;
}

template<class T>
bool gsFixedConstraint<T>::correct(gsStepContext<T> & ctx, const gsVector<T> & dUbar,
                                   const gsVector<T> &)
{
    if (dUbar.size() != ctx.DeltaU.size())
        return false;
    ctx.deltaU = dUbar;
    ctx.deltaL = 0;
    return true;
}

template<class T>
T gsFixedConstraint<T>::residual(const gsStepContext<T> &) const
{
    return T(0);
}

template<class T>
T gsFixedConstraint<T>::distance(const gsVector<T> & DeltaU, T) const
{
    return DeltaU.norm();
}

template<class T>
bool gsFixedConstraint<T>::stepIsLoadIncrement() const
{
    return false;
}

template<class T>
std::string gsFixedConstraint<T>::name() const
{
    return "Fixed";
}

// ----------------------------------------------------------------- Load

template<class T>
bool gsLoadConstraint<T>::predictorNeedsTangent(const gsStepContext<T> &) const
{
    return true;
}

template<class T>
bool gsLoadConstraint<T>::predict(gsStepContext<T> & ctx, const gsVector<T> & dUt,
                                  const gsStateSpace<T> & space)
{
    if (dUt.size() != space.tangentSize())
        return false;
    ctx.DeltaL = ctx.deltaL = ctx.stepSize;
    m_target = ctx.DeltaL;
    ctx.deltaU = ctx.deltaL*dUt;
    ctx.DeltaU = ctx.deltaU;
    return true;
}

template<class T>
bool gsLoadConstraint<T>::correct(gsStepContext<T> & ctx, const gsVector<T> & dUbar,
                                  const gsVector<T> &)
{
    if (dUbar.size() != ctx.DeltaU.size())
        return false;
    ctx.deltaU = dUbar;
    ctx.deltaL = 0;
    return true;
}

template<class T>
T gsLoadConstraint<T>::residual(const gsStepContext<T> & ctx) const
{
    const T d = ctx.DeltaL - m_target;
    return (m_target == 0) ? d : d / math::abs(m_target);
}

template<class T>
T gsLoadConstraint<T>::distance(const gsVector<T> &, T DeltaL) const
{
    return math::abs(DeltaL);
}

template<class T>
bool gsLoadConstraint<T>::stepIsLoadIncrement() const
{
    return true;
}

template<class T>
std::string gsLoadConstraint<T>::name() const
{
    return "Load";
}

// --------------------------------------------------------- Displacement

template<class T>
gsDisplacementConstraint<T>::gsDisplacementConstraint(index_t dof)
: m_dof(dof)
{
    GISMO_ENSURE(dof >= 0, "gsDisplacementConstraint: the dof index must be non-negative, got " << dof);
}

template<class T>
void gsDisplacementConstraint<T>::setDof(index_t dof)
{
    GISMO_ENSURE(dof >= 0, "gsDisplacementConstraint::setDof: the dof index must be non-negative, got " << dof);
    m_dof = dof;
}

template<class T>
bool gsDisplacementConstraint<T>::predictorNeedsTangent(const gsStepContext<T> &) const
{
    return true;
}

template<class T>
bool gsDisplacementConstraint<T>::predict(gsStepContext<T> & ctx, const gsVector<T> & dUt,
                                          const gsStateSpace<T> & space)
{
    const index_t n = space.tangentSize();
    if (m_dof >= n || dUt.size() != n || dUt(m_dof) == 0)
        return false;
    ctx.deltaL = ctx.stepSize / dUt(m_dof);
    ctx.deltaU = ctx.deltaL*dUt;
    ctx.DeltaL = ctx.deltaL;
    ctx.DeltaU = ctx.deltaU;
    return true;
}

template<class T>
bool gsDisplacementConstraint<T>::correct(gsStepContext<T> & ctx, const gsVector<T> & dUbar,
                                          const gsVector<T> & dUt)
{
    if (dUbar.size() != ctx.DeltaU.size() || dUt.size() != ctx.DeltaU.size()
        || m_dof >= ctx.DeltaU.size() || dUt(m_dof) == 0)
        return false;
    const T g = residual(ctx);
    ctx.deltaL = -( g + dUbar(m_dof) ) / dUt(m_dof);
    ctx.deltaU = dUbar + ctx.deltaL*dUt;
    return true;
}

template<class T>
T gsDisplacementConstraint<T>::residual(const gsStepContext<T> & ctx) const
{
    return m_dof < ctx.DeltaU.size() ? ctx.DeltaU(m_dof) - ctx.stepSize
                                     : std::numeric_limits<T>::quiet_NaN();
}

template<class T>
T gsDisplacementConstraint<T>::distance(const gsVector<T> & DeltaU, T) const
{
    return m_dof < DeltaU.size() ? math::abs(DeltaU(m_dof))
                                 : std::numeric_limits<T>::quiet_NaN();
}

template<class T>
bool gsDisplacementConstraint<T>::stepIsLoadIncrement() const
{
    return false;
}

template<class T>
std::string gsDisplacementConstraint<T>::name() const
{
    return "Displacement";
}

// ----------------------------------------------------------------- Riks

template<class T>
T gsRiksConstraint<T>::scaling(index_t n)
{
    return n == 1 ? T(0.99999999999) : (n > 1 ? T(1)/static_cast<T>(n) : T(0));
}

template<class T>
bool gsRiksConstraint<T>::firstStep(const gsStepContext<T> & ctx) const
{
    const T tol = 1e-10;
    return !ctx.prev.valid
        || ctx.prev.DeltaU.size() != ctx.F.size()
        || ( ctx.prev.DeltaU.norm() < tol && ctx.prev.DeltaL*ctx.prev.DeltaL < tol );
}

template<class T>
bool gsRiksConstraint<T>::predictorNeedsTangent(const gsStepContext<T> & ctx) const
{
    return firstStep(ctx);
}

template<class T>
bool gsRiksConstraint<T>::predict(gsStepContext<T> & ctx, const gsVector<T> & dUt,
                                  const gsStateSpace<T> & space)
{
    const index_t n = space.tangentSize();
    const T phi = scaling(n);
    const T s = ctx.stepSize;

    if (firstStep(ctx))
    {
        if (dUt.size() != n)
            return false;
        T DL = 1.;
        const T nrm = math::sqrt( phi * dUt.dot(dUt) + (1.0 - phi) * DL * DL );
        ctx.deltaL = DL / nrm;
        ctx.deltaU = dUt * (DL / nrm);
    }
    else
    {
        const T sPrev = ctx.prev.stepSize;
        ctx.deltaL = 1./sPrev*(ctx.prev.DeltaL);
        ctx.deltaU = 1./sPrev*(ctx.prev.DeltaU);
    }

    ctx.deltaL *= s;
    ctx.deltaU *= s;
    ctx.DeltaL = ctx.deltaL;
    ctx.DeltaU = ctx.deltaU;
    return true;
}

template<class T>
bool gsRiksConstraint<T>::correct(gsStepContext<T> & ctx, const gsVector<T> & dUbar,
                                  const gsVector<T> & dUt)
{
    if (dUbar.size() != ctx.DeltaU.size() || dUt.size() != ctx.DeltaU.size())
        return false;
    const T phi = scaling(ctx.DeltaU.size());
    const T r = residual(ctx);
    ctx.deltaL = - ( r + 2*phi*(ctx.DeltaU).dot(dUbar) ) / ( 2*(1-phi)*(ctx.DeltaL) + 2*phi*(ctx.DeltaU).dot(dUt) );
    ctx.deltaU = dUbar + ctx.deltaL*dUt;
    return true;
}

template<class T>
T gsRiksConstraint<T>::residual(const gsStepContext<T> & ctx) const
{
    const T phi = scaling(ctx.DeltaU.size());
    const T s = ctx.stepSize;
    const T s2 = s*s;
    return phi*ctx.DeltaU.dot(ctx.DeltaU) + (1-phi)*ctx.DeltaL*ctx.DeltaL - s2;
}

template<class T>
T gsRiksConstraint<T>::distance(const gsVector<T> & DeltaU, T DeltaL) const
{
    const T phi = scaling(DeltaU.size());
    return math::pow(phi * math::pow(DeltaU.norm(),2.0) + (1.0-phi) * math::pow(DeltaL,2.0),0.5);
}

template<class T>
bool gsRiksConstraint<T>::stepIsLoadIncrement() const
{
    return false;
}

template<class T>
std::string gsRiksConstraint<T>::name() const
{
    return "Riks";
}

} // namespace gismo
