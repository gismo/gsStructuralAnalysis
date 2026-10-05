/** @file gsStage0Oscillator.h

    @brief Dependency-free nonlinear 2-dof oscillator used as a dynamics regression problem.

    Chain: wall --spring 1-- mass m1 --spring 2-- mass m2, with a harmonic load on m2.
    Every spring is hardening, f(e) = k e + k3 e^3, with extensions e1 = u1 and e2 = u2 - u1.
    The tangent f'(e) = k + 3 k3 e^2 is positive, so K(u) is symmetric positive definite.

    Sign convention of the residual op: R(u,t) = Fext(t) - Fint(u). This is the convention
    gsDynamicNewmark expects; the solver iterates on M A + C V = R.

    Author(s): H.M.Verhelst
 **/

#pragma once

#include <cmath>
#include <gsCore/gsLinearAlgebra.h>
#include <gsStructuralAnalysis/src/gsStructuralAnalysisTools/gsStructuralAnalysisTypes.h>

namespace gismo
{

/**
   @brief Nonlinear 2-dof spring-mass chain with a harmonic load.

   All vectors are 2-vectors and all matrices are 2x2 sparse matrices.
   The ops returned by the accessors capture `this`, hence the object is neither
   copyable nor movable and must outlive every solver built from its ops.
   The call counters live here (not in the lambdas) because solvers copy the
   std::function objects.
 */
template <class T>
class gsStage0Oscillator
{
public:
    /// Mass of the first body
    T m1;
    /// Mass of the second body
    T m2;
    /// Linear spring stiffness
    T k;
    /// Cubic spring stiffness
    T k3;
    /// Damping coefficient, C = diag(c,c)
    T c;
    /// Load amplitude
    T F0;
    /// Load angular frequency
    T omega;

    /// Number of Jacobian-op calls since the last resetCounters()
    index_t jacobianCalls;
    /// Number of residual-op calls since the last resetCounters()
    index_t residualCalls;

    gsStage0Oscillator()
    : m1(1), m2(1), k(1), k3(10), c(T(0.05)), F0(1), omega(2),
      jacobianCalls(0), residualCalls(0)
    {}

    gsStage0Oscillator(const gsStage0Oscillator &) = delete;
    gsStage0Oscillator(gsStage0Oscillator &&) = delete;
    gsStage0Oscillator & operator=(const gsStage0Oscillator &) = delete;
    gsStage0Oscillator & operator=(gsStage0Oscillator &&) = delete;

    /// Sets both call counters to zero
    void resetCounters() { jacobianCalls = 0; residualCalls = 0; }

    /// Spring force for extension \a e
    T spring(const T e) const { return k*e + k3*e*e*e; }
    /// Spring tangent stiffness for extension \a e
    T springTangent(const T e) const { return k + 3*k3*e*e; }

    /// Internal force Fint(u) = [f(e1) - f(e2), f(e2)], a 2-vector
    gsVector<T> Fint(const gsVector<T> & u) const
    {
        const T f1 = spring(u(0)), f2 = spring(u(1) - u(0));
        gsVector<T> F(2);
        F << f1 - f2, f2;
        return F;
    }

    /// Tangent stiffness K(u) = dFint/du, a 2x2 sparse matrix
    gsSparseMatrix<T> tangent(const gsVector<T> & u) const
    {
        const T t1 = springTangent(u(0)), t2 = springTangent(u(1) - u(0));
        gsSparseMatrix<T> K(2,2);
        K.insert(0,0) = t1 + t2; K.insert(0,1) = -t2;
        K.insert(1,0) = -t2;     K.insert(1,1) = t2;
        K.makeCompressed();
        return K;
    }

    /// Constant mass matrix diag(m1,m2), a 2x2 sparse matrix
    gsSparseMatrix<T> mass() const
    {
        gsSparseMatrix<T> M(2,2);
        M.insert(0,0) = m1; M.insert(1,1) = m2;
        M.makeCompressed();
        return M;
    }

    /// Constant damping matrix diag(c,c), always sized 2x2
    gsSparseMatrix<T> damping() const
    {
        gsSparseMatrix<T> C(2,2);
        C.insert(0,0) = c; C.insert(1,1) = c;
        C.makeCompressed();
        return C;
    }

    /// External force F(t) = [0, F0 sin(omega t)], a 2-vector
    gsVector<T> force(const T t) const
    {
        gsVector<T> F(2);
        F << 0, F0*std::sin(omega*t);
        return F;
    }

    /// Residual Fext(t) - Fint(u), a 2-vector
    gsVector<T> residual(const gsVector<T> & u, const T t) const
    {
        return force(t) - Fint(u);
    }

    /// Acceleration consistent with the equation of motion, M^{-1} (F(t0) - Fint(U0) - C V0)
    gsVector<T> initialAcceleration(const gsVector<T> & U0, const gsVector<T> & V0, const T t0) const
    {
        gsVector<T> r = residual(U0,t0) - c*V0;
        r(0) /= m1;
        r(1) /= m2;
        return r;
    }

    /// Mass op, writes M
    typename gsStructuralAnalysisOps<T>::Mass_t massOp()
    {
        return [this](gsSparseMatrix<T> & M) { M = this->mass(); return true; };
    }

    /// Damping op, writes C (independent of the state)
    typename gsStructuralAnalysisOps<T>::Damping_t dampingOp()
    {
        return [this](const gsVector<T> &, gsSparseMatrix<T> & C) { C = this->damping(); return true; };
    }

    /// Jacobian op, writes K(u) and counts the call
    typename gsStructuralAnalysisOps<T>::Jacobian_t jacobianOp()
    {
        return [this](const gsVector<T> & u, gsSparseMatrix<T> & K)
        { ++this->jacobianCalls; K = this->tangent(u); return true; };
    }

    /// Time-dependent residual op, writes R = Fext(t) - Fint(u) and counts the call
    typename gsStructuralAnalysisOps<T>::TResidual_t residualOp()
    {
        return [this](const gsVector<T> & u, const T t, gsVector<T> & R)
        { ++this->residualCalls; R = this->residual(u,t); return true; };
    }

    /// Time-dependent mass op, writes M
    typename gsStructuralAnalysisOps<T>::TMass_t tMassOp()
    {
        return [this](const T, gsSparseMatrix<T> & M) { M = this->mass(); return true; };
    }

    /// Time-dependent damping op, writes C
    typename gsStructuralAnalysisOps<T>::TDamping_t tDampingOp()
    {
        return [this](const gsVector<T> &, const T, gsSparseMatrix<T> & C) { C = this->damping(); return true; };
    }

    /// Time-dependent Jacobian op, writes K(u) and counts the call
    typename gsStructuralAnalysisOps<T>::TJacobian_t tJacobianOp()
    {
        return [this](const gsVector<T> & u, const T, gsSparseMatrix<T> & K)
        { ++this->jacobianCalls; K = this->tangent(u); return true; };
    }
};

} // namespace gismo
