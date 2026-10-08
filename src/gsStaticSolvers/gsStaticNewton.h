 /** @file gsStaticNewton.h

    @brief Static solver using a newton method

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#include <gsStructuralAnalysis/src/gsStaticSolvers/gsStaticBase.h>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsStateSpace.h>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsContinuationProblem.h>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsConvergenceTest.h>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsPathConstraint.h>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsContinuationDriver.h>
#pragma once


namespace gismo
{

/**
    @brief Relative residual test of the static Newton solver:
    \f$\|R\|/\text{basis} < \text{tol}\f$, tested strictly.

    The basis is the norm of the residual at the predicted point (the
    committed state, or the committed state plus the headstart increment),
    unless an external basis was armed with setBasis() for the next step
    attempt. A floor armed with setBasisFloor() bounds the basis from below:
    basis = max(basis, floor). A zero basis is replaced by 1, so that a purely
    displacement-driven problem, whose residual is zero at the start, still
    has a finite ratio.

    Theory of the floor: at an equilibrium \f$R(U_n,\lambda)=0\f$, so
    \f$R(U_n,\lambda+\Delta\lambda)\approx\Delta\lambda\,\partial R/\partial\lambda\f$
    is, to first order, the initial residual of a fixed-load Newton solve. When
    the predictor is exact (linear problems) the residual at the predicted
    point is round-off and must not serve as the basis.

    Cost of test(): O(n).

    \tparam T coefficient type

    \ingroup gsStructuralAnalysis
*/
template <class T>
class gsStaticResidualTest : public gsConvergenceTest<T>
{
public:
    /// @param[in] tol tolerance of the relative residual
    explicit gsStaticResidualTest(T tol)
    : m_tol(tol), m_armed(false), m_floorArmed(false), m_ext(T(1)), m_floor(T(0)), m_basis(T(1)), m_pred(T(1)),
      m_value(std::numeric_limits<T>::quiet_NaN())
    {}

    /// Arms @a basis as the basis of the next initStep(), which disarms it again
    void setBasis(T basis) { m_ext = basis; m_armed = true; }

    /// Arms @a floor as a lower bound of the basis of the next initStep(), which disarms it again
    void setBasisFloor(T floor) { m_floor = floor; m_floorArmed = true; }

    /// The basis in use; valid after initStep()
    T basis() const { return m_basis; }

    /// Norm of the residual at the predicted point (1 if it is zero); valid after initStep()
    T predictedResidual() const { return m_pred; }

    void initStep(const gsStepContext<T> & ctx, const gsStateSpace<T> &)
    {
        m_pred = ctx.R.norm();
        if (m_pred == 0) m_pred = T(1);
        m_basis = m_armed ? m_ext : m_pred;
        if (m_floorArmed && m_floor > m_basis) m_basis = m_floor;
        if (m_basis == 0) m_basis = T(1);
        m_armed = m_floorArmed = false;
    }

    gsConvergenceStatus test(const gsStepContext<T> & ctx, const gsStateSpace<T> &)
    {
        m_value = ctx.R.norm() / m_basis;
        return (m_value < m_tol) ? gsConvergenceStatus::Converged : gsConvergenceStatus::Unconverged;
    }

    T value() const { return m_value; }
    std::string name() const { return "StaticRelResidual"; }

private:
    T m_tol;
    bool m_armed, m_floorArmed;
    T m_ext, m_floor, m_basis, m_pred, m_value;
};

/**
    @brief Relative update test of the static Newton solver:
    \f$\|\delta U\|/\|\Delta U\| < \text{tol}\f$, tested strictly.

    \f$\delta U\f$ is the latest applied (relaxed) update and \f$\Delta U\f$
    the running increment of the step including that update. The basis grows
    with the iteration, unlike gsRelUpdateTest, whose basis is the predictor
    increment and is zero for a fixed-load step. A non-finite ratio is never
    Converged.

    Cost of test(): O(n).

    \tparam T coefficient type

    \ingroup gsStructuralAnalysis
*/
template <class T>
class gsStaticUpdateTest : public gsConvergenceTest<T>
{
public:
    /// @param[in] tol tolerance of the relative update
    explicit gsStaticUpdateTest(T tol)
    : m_tol(tol), m_value(std::numeric_limits<T>::quiet_NaN())
    {}

    gsConvergenceStatus test(const gsStepContext<T> & ctx, const gsStateSpace<T> &)
    {
        m_value = ctx.deltaU.norm() / ctx.DeltaU.norm();
        return (m_value < m_tol) ? gsConvergenceStatus::Converged : gsConvergenceStatus::Unconverged;
    }

    T value() const { return m_value; }
    std::string name() const { return "StaticRelUpdate"; }

private:
    T m_tol, m_value;
};

/**
    @brief Fixed-load constraint whose predictor is a given increment
    \f$\Delta U\f$ instead of zero (the headstart of a static solve).

    \tparam T coefficient type

    \ingroup gsStructuralAnalysis
*/
template <class T>
class gsStaticHeadstartConstraint : public gsFixedConstraint<T>
{
public:
    /// @param[in] DeltaU the predicted increment
    explicit gsStaticHeadstartConstraint(const gsVector<T> & DeltaU)
    : m_DeltaU(DeltaU)
    {}

    /// DeltaU = deltaU = the given increment, DeltaL = 0; false iff its size differs from space.tangentSize()
    bool predict(gsStepContext<T> & ctx, const gsVector<T> &,
                 const gsStateSpace<T> & space)
    {
        if (m_DeltaU.size() != space.tangentSize())
            return false;
        ctx.DeltaU = ctx.deltaU = m_DeltaU;
        ctx.DeltaL = ctx.deltaL = 0;
        return true;
    }

private:
    gsVector<T> m_DeltaU;
};

/**
    @brief Static solver using a newton method

    \tparam T coefficient type

    \ingroup gsStructuralAnalysis
*/
template <class T>
class gsStaticNewton : public gsStaticBase<T>
{
protected:

    typedef gsStaticBase<T> Base;

    typedef typename Base::Residual_t    Residual_t;
    typedef typename Base::ALResidual_t  ALResidual_t;
    typedef typename Base::Jacobian_t    Jacobian_t;
    typedef typename Base::dJacobian_t   dJacobian_t;

public:

    /**
     * @brief      Constructor
     *
     * @param[in]  linear  The linear stiffness matrix
     * @param[in]  force   The external force
     */
    gsStaticNewton(     const gsSparseMatrix<T> &linear,
                        const gsVector<T> &force        )
    :
        m_linear(linear),
        m_force(force),
        m_nonlinear(nullptr),
        m_residualFun(nullptr)
    {
        this->_init();
    }

    /**
     * @brief      Constructor
     *
     * @param[in]  linear     The linear stiffness matrix
     * @param[in]  force      The external force
     * @param[in]  nonlinear  The Jacobian
     * @param[in]  residual   The residual
     */
    gsStaticNewton( const gsSparseMatrix<T> &linear,
                    const gsVector<T> &force,
                    const Jacobian_t &nonlinear,
                    const Residual_t &residual      )
    :
        m_linear(linear),
        m_force(force),
        m_nonlinear(nonlinear),
        m_dnonlinear(nullptr),
        m_residualFun(residual),
        m_ALresidualFun(nullptr)
    {
        m_dnonlinear = [this](gsVector<T> const & x, gsVector<T> const & /*dx*/, gsSparseMatrix<T> & m) -> bool
        {
            return m_nonlinear(x,m);
        };

        this->_init();
    }

    /**
     * @brief      Constructor
     *
     * @param[in]  linear     The linear stiffness matrix
     * @param[in]  force      The external force
     * @param[in]  nonlinear  The Jacobian
     * @param[in]  residual   The residual as arc-length object
     */
    gsStaticNewton( const gsSparseMatrix<T> &linear,
                    const gsVector<T>  &force,
                    const Jacobian_t   &nonlinear,
                    const ALResidual_t &ALResidual  )
    :
        m_linear(linear),
        m_force(force),
        m_nonlinear(nonlinear),
        m_dnonlinear(nullptr),
        m_residualFun(nullptr),
        m_ALresidualFun(ALResidual)
    {
        m_L = 1.0;
        m_residualFun = [this](gsVector<T> const & x, gsVector<T> & result) -> bool
        {
            return m_ALresidualFun(x,m_L,result);
        };

        m_dnonlinear = [this](gsVector<T> const & x, gsVector<T> const & /*dx*/, gsSparseMatrix<T> & m) -> bool
        {
            return m_nonlinear(x,m);
        };

        this->_init();
    }

    /**
     * @brief      { function_description }
     *
     * @param[in]  linear      The linear matrix
     * @param[in]  force       The force vector
     * @param[in]  dnonlinear  The jacobian taking the solution x and the iterative update dx
     * @param[in]  residual    The residual function
     */
    gsStaticNewton(   const gsSparseMatrix<T> &linear,
                      const gsVector<T> &force,
                      const dJacobian_t &dnonlinear,
                      const Residual_t &residual
                      )
    :
        m_linear(linear),
        m_force(force),
        m_nonlinear(nullptr),
        m_dnonlinear(dnonlinear),
        m_residualFun(residual),
        m_ALresidualFun(nullptr)
    {
        this->_init();
    }

public:

    /**
        @brief Solves the linear or the nonlinear problem, see \ref gsStaticBase.

        The nonlinear solve is a Newton iteration at the fixed load
        \f$\lambda\f$ = the load set by setLoad(), run by a
        gsContinuationDriver with a gsFixedConstraint. Each update solves
        \f$K\,\delta U = R\f$ with \f$R\f$ the residual of the user, which is
        handed to the driver as it is. The solve is accepted after at least
        one update when both
        \f$\|R\|/\|R_0\| < \texttt{tolF}\f$ and
        \f$\|\delta U\|/\|\Delta U\| < \texttt{tolU}\f$ hold, with
        \f$R_0\f$ the residual at the start and \f$\Delta U\f$ the running
        increment. At most <tt>maxIt</tt> updates are applied.

        On Success, solution() is advanced by update(). On any other status,
        in particular NotConverged after <tt>maxIt</tt> updates, solution()
        is unchanged.

        @return the status, see gsStatus
    */
    gsStatus solve() override
    {
        if (m_NL)
            return solveNonlinear();
        else
            return solveLinear();
    }
    /// Perform a linear solve
    gsStatus solveLinear();
    gsStatus solveLinear(gsVector<T> & solution)
    {
        gsStatus status = this->solveLinear();
        solution = m_U;
        return status;
    }
    /// Perform a nonlinearg solve
    gsStatus solveNonlinear(gsVector<T> & solution)
    {
        gsStatus status = this->solveNonlinear();
        solution = m_U;
        return status;
    }
    gsStatus solveNonlinear();

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

    using Base::indicator;
    using Base::stabilityVec;

    /// See \ref gsStaticBase
    T indicator()
    {
        gsSparseMatrix<T> m;
        GISMO_ENSURE(m_nonlinear(m_U, m),"Assembly failed");
        return indicator(m);
    }

    /// See \ref gsStaticBase
    gsVector<T> stabilityVec()
    {
        gsSparseMatrix<T> m;
        GISMO_ENSURE(m_nonlinear(m_U, m),"Assembly failed");
        return stabilityVec(m);
    }

protected:

    /// Perform a linear solve
    gsVector<T> _solveLinear();
    /// Perform a nonlinear solve on the continuation driver
    gsStatus _solveNonlinear();

    using Base::_computeStability;

    /// Initializes the method
    void _init();
    /// Factorizes the \a jacMat
    void _factorizeMatrix(const gsSparseMatrix<T> & jacMat) const;
    /// Solves the system with RHS \a F and LHS the Jacobian
    gsVector<T> _solveSystem(const gsVector<T> & F);

    gsVector<T> _computeResidual(const gsVector<T> & U);

protected:

    const gsSparseMatrix<T> & m_linear;
    const gsVector<T> & m_force;
    const Jacobian_t      m_nonlinear;
          dJacobian_t     m_dnonlinear;
          Residual_t      m_residualFun;
    const ALResidual_t    m_ALresidualFun;

    // The nonlinear solve is rebuilt per solve: options and m_U may change in between.
    // Declaration order is destruction order reversed: the driver refers to the problem and the space.
    gsContinuationProblem<T> m_problem;
    memory::unique_ptr< gsVectorStateSpace<T> >      m_space;
    memory::unique_ptr< gsContinuationDriver<T> >    m_driver;
    memory::shared_ptr< gsStaticResidualTest<T> >    m_resTest;

    using Base::m_R;

    using Base::m_U;
    using Base::m_DeltaU;
    using Base::m_deltaU;

    using Base::m_L;

    using Base::m_stabilityVec;

    using Base::m_verbose;

    bool m_NL;

    using Base::m_options;

    using Base::m_tolF;
    using Base::m_tolU;

    mutable T m_relax;

    using Base::m_maxIterations;
    using Base::m_numIterations;

    // Residual norms
    using Base::m_residual;
    using Base::m_residualIni;
    using Base::m_residualOld;

    using Base::m_start;
    using Base::m_headstart;

    /// Linear solver employed
    using Base::m_solver;   // Cholesky by default

    using Base::m_stabilityMethod;

    /// Indicator for bifurcation
    using Base::m_indicator;

    using Base::m_dofs;

    // Solver status
    using Base::m_status;
};


} // namespace gismo


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsStaticNewton.hpp)
#endif
