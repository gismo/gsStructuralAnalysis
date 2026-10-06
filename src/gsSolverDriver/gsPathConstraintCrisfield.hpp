 /** @file gsPathConstraintCrisfield.hpp

    @brief Implementation of the Crisfield path constraints.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#pragma once
#include <algorithm>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsPathConstraintCrisfield.h>

namespace gismo
{

// ------------------------------------------------------------- Crisfield

template<class T>
void gsCrisfieldConstraint<T>::setAngleMethod(index_t m)
{
    GISMO_ENSURE(m >= 0 && m <= 2, "gsCrisfieldConstraint::setAngleMethod: the angle method must be 0, 1 or 2, got " << m);
    m_angle = m;
}

template<class T>
void gsCrisfieldConstraint<T>::resetCounters()
{
    m_last = Root::None;
    m_numModified = 0;
    m_numSimplified = 0;
    m_numRelaxed = 0;
}

template<class T>
bool gsCrisfieldConstraint<T>::firstStep(const gsStepContext<T> & ctx) const
{
    return !ctx.prev.valid
        || ctx.prev.DeltaU.size() != ctx.F.size()
        || ( ctx.prev.DeltaU.dot(ctx.prev.DeltaU) == 0 && ctx.prev.DeltaL*ctx.prev.DeltaL == 0 );
}

template<class T>
bool gsCrisfieldConstraint<T>::predictorNeedsTangent(const gsStepContext<T> &) const
{
    return true;
}

template<class T>
bool gsCrisfieldConstraint<T>::predict(gsStepContext<T> & ctx, const gsVector<T> & dUt,
                                       const gsStateSpace<T> & space)
{
    const index_t n = space.tangentSize();
    if (dUt.size() != n || ctx.F.size() != n)
        return false;
    const T s = ctx.stepSize;
    m_FF = ctx.F.dot(ctx.F);

    if (firstStep(ctx))
    {
        // With the automatic scaling phi = |u_t|/|F|:
        // |DeltaU|^2 + phi^2 |F|^2 DeltaL^2 = s^2/2 + s^2/2 = s^2.
        ctx.deltaL = s / math::pow(2*( dUt.dot(dUt) ) , 0.5);
        ctx.deltaU = ctx.deltaL*dUt;
        if (!m_phiUser)
            m_phi = math::pow( dUt.dot(dUt) / ctx.F.dot(ctx.F),0.5);
    }
    else
    {
        if (!m_phiUser)
        {
            // phi = |U| / (|lambda| |F|) at the committed state; unguarded, so
            // lambda = 0 gives an infinite or NaN value.
            gsVector<T> U;
            if (!space.magnitude(U) || U.size() == 0)
            {
                if (!space.exportWorking(U))
                    return false;
            }
            if (U.size() != n)
                return false;
            m_phi = math::pow(U.dot(U)/( math::pow(ctx.lambda,2) * ctx.F.dot(ctx.F) ),0.5);
        }

        // MU predictor, previous increment (oldU, oldL).
        const bool usePred = (m_angle == angle::Predictor) && m_DeltaUpred.size() == n;
        const gsVector<T> & oldU = usePred ? m_DeltaUpred : ctx.prev.DeltaU;
        const T oldL = usePred ? m_DeltaLpred : ctx.prev.DeltaL;

        T A0 = math::pow(m_phi,2)* ctx.F.dot(ctx.F);
        index_t dir = math::getSign(oldU.dot(dUt) + A0*oldL);
        T denum = ( math::pow( dUt.dot(dUt) + A0 ,0.5) );
        T mu = (denum==0) ? s : s / denum;
        ctx.deltaL = dir*mu;
        ctx.deltaU = ctx.deltaL*dUt;
    }

    ctx.DeltaL = ctx.deltaL;
    ctx.DeltaU = ctx.deltaU;

    if (m_angle == angle::Predictor)
    {
        m_DeltaUpred = ctx.DeltaU;
        m_DeltaLpred = ctx.DeltaL;
    }
    return true;
}

template<class T>
bool gsCrisfieldConstraint<T>::correct(gsStepContext<T> & ctx, const gsVector<T> & dUbar,
                                       const gsVector<T> & dUt)
{
    const index_t n = ctx.DeltaU.size();
    if (dUbar.size() != n || dUt.size() != n || ctx.F.size() != n || ctx.R.size() != n)
        return false;
    const T s = ctx.stepSize;

    // Previous increment (oldU, oldL) of the angle criterion
    const gsVector<T> zero = gsVector<T>::Zero(n);
    const bool usePred = (m_angle == angle::Predictor) && m_DeltaUpred.size() == n;
    const bool useIter = (m_angle == angle::Iteration);
    const bool useStep = !useIter && !usePred && !firstStep(ctx);
    const gsVector<T> & oldU = useIter ? ctx.DeltaU : (usePred ? m_DeltaUpred : (useStep ? ctx.prev.DeltaU : zero));
    const T oldL = useIter ? ctx.DeltaL : (usePred ? m_DeltaLpred : (useStep ? ctx.prev.DeltaL : T(0)));

    const T lamold = ctx.deltaL;
    T eta = 1;

    // The quadratic a0 dl^2 + (b0 + eta b1) dl + (c0 + eta c1 + eta^2 c2) = 0
    T A0 = math::pow(m_phi,2)* ctx.F.dot(ctx.F);
    const T a0 = dUt.dot(dUt) + A0;
    const T b0 = 2*( dUt.dot(ctx.DeltaU) + ctx.DeltaL * A0 );
    const T b1 = 2*( dUbar.dot(dUt) );
    const T c0 = ctx.DeltaU.dot(ctx.DeltaU) + ctx.DeltaL*ctx.DeltaL * A0 - math::pow(s,2);
    const T c1 = 2*( ctx.DeltaU.dot(dUbar) );
    const T c2 = dUbar.dot(dUbar);

    T alpha1 = a0;
    T alpha2 = b0 + eta*b1;
    T alpha3 = c0 + eta*c1 + eta*eta*c2;
    T disc = math::pow(alpha2 ,2) - 4 * alpha1 * alpha3;

    // Root choice by the angle criterion of Ritto-Correa and Camotim (2008): the
    // root with the larger inner product with the previous increment. The load
    // term is weighted by phi^2, not by psi = phi^2 |F|^2 as in the MU predictor.
    auto chooseRoot = [&](const T dl0, const T dl1)
    {
        const gsVector<T> deltaU1 = eta*dUbar + dUt*dl0;
        const gsVector<T> deltaU2 = eta*dUbar + dUt*dl1;
        const T DOT1 = dl0*(oldU.dot(dUt) + math::pow(m_phi,2)*oldL);
        const T DOT2 = dl1*(oldU.dot(dUt) + math::pow(m_phi,2)*oldL);
        if (DOT1 > DOT2)
        {
            ctx.deltaL = dl0;
            ctx.deltaU = deltaU1;
        }
        else if (DOT1 < DOT2)
        {
            ctx.deltaL = dl1;
            ctx.deltaU = deltaU2;
        }
        else
        {
            ctx.deltaL = dl0;
            ctx.deltaU = deltaU1;
        }
    };

    if (disc >= 0)
    {
        eta = 1;
        const T dl0 = (-alpha2 + math::pow(disc,0.5))/(2*alpha1);
        const T dl1 = (-alpha2 - math::pow(disc,0.5))/(2*alpha1);
        chooseRoot(dl0, dl1);
        m_last = Root::Real;
    }
    else
    {
        // Lam and Morley (1992): eta that makes the discriminant of the quadratic vanish
        alpha1 = b1*b1 - 4.0*a0*c2;
        alpha2 = 2.0*b0*b1 - 4.0*a0*c1;
        alpha3 = b0*b0 - 4.0*a0*c0;
        disc = math::pow(alpha2 ,2.0) - 4.0 * alpha1 * alpha3;

        if (disc >= 0)
        {
            const T e0 = (-alpha2 + math::pow(disc,0.5))/(2.0*alpha1);
            const T e1 = (-alpha2 - math::pow(disc,0.5))/(2.0*alpha1);
            const T eta1 = std::min(e0,e1);
            const T eta2 = std::max(e0,e1);
            const T xi = 0.05*math::abs(eta2-eta1);
            if (eta2<1.0)
                eta = eta2-xi;
            else if ( (eta2 > 1.0) && (-alpha2/alpha1 < 1.0) )
                eta = eta2+xi;
            else if ( (eta1 < 1.0) && (-alpha2/alpha1 > 1.0) )
                eta = eta1-xi;
            else if ( eta1 > 1.0 )
                eta = eta1 + xi;
        }

        if ((disc >= 0) && eta > 0.05)
        {
            // double root of the quadratic at the modified eta
            alpha1 = a0;
            alpha2 = b0 + eta*b1;
            alpha3 = c0 + eta*c1 + eta*eta*c2;
            const T dl0 = (-alpha2)/(2*alpha1);
            const T dl1 = (-alpha2)/(2*alpha1);
            chooseRoot(dl0, dl1);
            m_last = Root::Modified;
            ++m_numModified;
        }
        else
        {
            // Rescale the Newton point onto the sphere; the internal force follows
            // from the residual R = (lambda+DeltaL) F - F_int at the current iterate.
            eta = 1;
            const gsVector<T> DeltaUcr = ctx.DeltaU + dUbar;
            const gsVector<T> Fint = (ctx.lambda + ctx.DeltaL)*ctx.F - ctx.R;
            const T Lcr = Fint.dot(ctx.F)/ctx.F.dot(ctx.F);
            const T DeltaLcr = Lcr - ctx.lambda;
            const T arcLength_cr = math::pow( DeltaUcr.dot(DeltaUcr) + A0 * math::pow(DeltaLcr, 2.0) ,0.5);
            const T mu = s/arcLength_cr;
            ctx.deltaL = mu*DeltaLcr - ctx.DeltaL;
            ctx.deltaU = mu*DeltaUcr - ctx.DeltaU;
            m_last = Root::Simplified;
            ++m_numSimplified;
        }
    }

    // Relaxation against an oscillating load factor: sign change without growth.
    // The caller scales deltaU and deltaL by ctx.relaxation.
    if ( ( (lamold*ctx.deltaL < 0) && (math::abs(ctx.deltaL) <= math::abs(lamold)) ) && m_relax != 1.0 )
    {
        ctx.deltaU = ctx.deltaL*dUt + eta*dUbar;
        ctx.relaxation = m_relax;
        ++m_numRelaxed;
    }
    return true;
}

template<class T>
T gsCrisfieldConstraint<T>::residual(const gsStepContext<T> & ctx) const
{
    const T s2 = ctx.stepSize*ctx.stepSize;
    return ctx.DeltaU.dot(ctx.DeltaU) + m_phi*m_phi*ctx.F.dot(ctx.F)*ctx.DeltaL*ctx.DeltaL - s2;
}

template<class T>
T gsCrisfieldConstraint<T>::distance(const gsVector<T> & DeltaU, T DeltaL) const
{
    T A0 = math::pow(m_phi,2)*m_FF;
    return math::pow(DeltaU.dot(DeltaU) + A0*math::pow(DeltaL,2.0),0.5);
}

template<class T>
bool gsCrisfieldConstraint<T>::stepIsLoadIncrement() const
{
    return false;
}

template<class T>
std::string gsCrisfieldConstraint<T>::name() const
{
    return "Crisfield";
}

// ----------------------------------------------------- ConsistentCrisfield

template<class T>
bool gsConsistentCrisfieldConstraint<T>::firstStep(const gsStepContext<T> & ctx) const
{
    const T tol = 1e-10;
    return !ctx.prev.valid
        || ctx.prev.DeltaU.size() != ctx.F.size()
        || ( ctx.prev.DeltaU.norm() < tol && ctx.prev.DeltaL*ctx.prev.DeltaL < tol );
}

template<class T>
bool gsConsistentCrisfieldConstraint<T>::predictorNeedsTangent(const gsStepContext<T> & ctx) const
{
    return firstStep(ctx) || (!m_phiUser && !m_phi0Set);
}

template<class T>
bool gsConsistentCrisfieldConstraint<T>::predict(gsStepContext<T> & ctx, const gsVector<T> & dUt,
                                                 const gsStateSpace<T> & space)
{
    const index_t n = space.tangentSize();
    if (ctx.F.size() != n)
        return false;
    const T s = ctx.stepSize;
    m_FF = ctx.F.dot(ctx.F);

    if (firstStep(ctx))
    {
        if (dUt.size() != n)
            return false;
        const T nrm = math::sqrt( dUt.dot(dUt) );
        ctx.deltaU = dUt / nrm;
        ctx.deltaL = T(1) / nrm;
        if (!m_phiUser && !m_phi0Set)
        {
            m_phi0 = math::pow( dUt.dot(dUt) / ctx.F.dot(ctx.F),0.5);
            m_phi0Set = true;
        }
    }
    else
    {
        const T sPrev = ctx.prev.stepSize;
        ctx.deltaL = 1./sPrev*(ctx.prev.DeltaL);
        ctx.deltaU = 1./sPrev*(ctx.prev.DeltaU);
        if (!m_phiUser && !m_phi0Set)
        {
            if (dUt.size() != n)
                return false;
            m_phi0 = math::pow( dUt.dot(dUt) / ctx.F.dot(ctx.F),0.5);
            m_phi0Set = true;
        }
    }

    if (!m_phiUser)
        m_phi = m_phi0;

    ctx.deltaL *= s;
    ctx.deltaU *= s;
    ctx.DeltaL = ctx.deltaL;
    ctx.DeltaU = ctx.deltaU;
    return true;
}

template<class T>
bool gsConsistentCrisfieldConstraint<T>::correct(gsStepContext<T> & ctx, const gsVector<T> & dUbar,
                                                 const gsVector<T> & dUt)
{
    const index_t n = ctx.DeltaU.size();
    if (dUbar.size() != n || dUt.size() != n || ctx.F.size() != n)
        return false;

    T cres = ctx.DeltaU.dot(ctx.DeltaU) + m_phi*m_phi*ctx.DeltaL*ctx.DeltaL * ctx.F.dot(ctx.F) - ctx.stepSize*ctx.stepSize;
    T num = cres + (2*ctx.DeltaU).dot(dUbar);
    T denum = (2*ctx.DeltaU).dot(dUt) + m_phi*m_phi*2*ctx.DeltaL*ctx.F.dot(ctx.F);
    ctx.deltaL = - num / denum;
    ctx.deltaU = ctx.deltaL * dUt + dUbar;
    return true;
}

template<class T>
T gsConsistentCrisfieldConstraint<T>::residual(const gsStepContext<T> & ctx) const
{
    const T s2 = ctx.stepSize*ctx.stepSize;
    return ctx.DeltaU.dot(ctx.DeltaU) + m_phi*m_phi*ctx.DeltaL*ctx.DeltaL*ctx.F.dot(ctx.F) - s2;
}

template<class T>
T gsConsistentCrisfieldConstraint<T>::distance(const gsVector<T> & DeltaU, T DeltaL) const
{
    return math::sqrt( DeltaU.dot(DeltaU) + m_phi*m_phi*DeltaL*DeltaL*m_FF );
}

template<class T>
bool gsConsistentCrisfieldConstraint<T>::stepIsLoadIncrement() const
{
    return false;
}

template<class T>
std::string gsConsistentCrisfieldConstraint<T>::name() const
{
    return "ConsistentCrisfield";
}

} // namespace gismo
