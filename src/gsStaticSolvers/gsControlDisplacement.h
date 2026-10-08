 /** @file gsControlDisplacement.h

    @brief Simple class for displacement control given a static solver

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#include <typeinfo>
#include <gsStructuralAnalysis/src/gsStaticSolvers/gsStaticNewton.h>
#include <gsStructuralAnalysis/src/gsStaticSolvers/gsContinuationBase.h>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsContinuationDriver.h>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsPathConstraint.h>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsContinuationProblem.h>
#include <gsStructuralAnalysis/src/gsSolverDriver/gsStateSpace.h>

#pragma once


namespace gismo
{

/**
    @brief Simple class for displacement control

    Two ways of construction, with the same public interface:
    - Given a static solver: each step() sets the load \f$\lambda\f$ of the
      solver, solves, and adds the update of the solver to the displacements.
    - Given a continuation problem and its state space: a preset over the
      driver's gsLoadConstraint, where \f$\lambda\f$ ramps the Dirichlet data.

    \tparam T coefficient type

    \ingroup gsStructuralAnalysis
*/
template <class T>
class gsControlDisplacement : public gsContinuationBase<T>
{
protected:

public:

    /// Constructor given a static solver
    gsControlDisplacement(gsStaticBase<T> * solver) :
    m_solver(solver),
    m_space(nullptr),
    m_problem(nullptr),
    first(true),
    m_status(gsStatus::NotStarted)
    {
        m_U = gsVector<T>::Zero(m_solver->numDofs());
        m_L = 0;
    }

    /**
        @brief Constructor given a continuation problem and its state space

        The step is a gsContinuationDriver step with a gsLoadConstraint:
        \f$\Delta\lambda\f$ is the step passed to step(), and the predictor
        is \f$\Delta U = \Delta\lambda\,K^{-1}\partial R/\partial\lambda\f$.
        The acceptance rule is that of gsStaticNewton with the defaults of
        gsStaticBase: at least one update, \f$\|R\|/\|R_0\| < 10^{-6}\f$ and
        \f$\|\delta U\|/\|\Delta U\| < 10^{-6}\f$ with at most 25 updates, where
        \f$\|R_0\|\f$ is bounded below by \f$|\Delta\lambda|\,\|\partial R/\partial\lambda\|\f$
        (see gsStaticResidualTest). The
        tangent solver is the default of gsTangentSolver. Both can be changed
        through driver().

        \f$\lambda\f$ is the parameter of the Dirichlet ramp and is applied by
        <tt>P.residual(lambda, ...)</tt>. <tt>P.loadDirection</tt> should return
        \f$\partial R/\partial\lambda\f$ of the ramp: without it the predictor
        uses <tt>P.F</tt> instead, and for <tt>P.F = 0</tt> the predictor is
        zero. The Load predictor costs one extra Jacobian assembly and
        factorization per step.

        @param[in] P the problem; it must outlive this object
        @param[in] S the state space; it must outlive this object
    */
    gsControlDisplacement(gsContinuationProblem<T> & P, gsStateSpace<T> & S) :
    m_solver(nullptr),
    m_space(&S),
    m_problem(&P),
    m_driver(new gsContinuationDriver<T>(P, S)),
    first(true),
    m_status(gsStatus::NotStarted)
    {
        typedef typename gsConvergenceTest<T>::Ptr TestPtr;
        m_driver->setConstraint(typename gsPathConstraint<T>::Ptr(new gsLoadConstraint<T>()));

        // The LoadNorm basis of the arc-length gate is zero for pure displacement control (F = 0)
        memory::shared_ptr< gsAllTest<T> > all(new gsAllTest<T>());
        m_resTest.reset(new gsStaticResidualTest<T>(T(1e-6)));
        all->add(TestPtr(new gsMinIterTest<T>(1)))
            .add(m_resTest)
            .add(TestPtr(new gsStaticUpdateTest<T>(T(1e-6))));
        memory::shared_ptr< gsAnyTest<T> > any(new gsAnyTest<T>());
        any->add(all).add(TestPtr(new gsMaxIterTest<T>(25)));
        m_driver->setTest(any);
        m_driver->setup();

        if (!S.exportWorking(m_U))
            m_U.resize(0);
        m_L = 0;
    }

    /// The driver of the constructor given a problem; throws (GISMO_ENSURE) if this object was constructed with a static solver
    gsContinuationDriver<T> & driver()
    {
        GISMO_ENSURE(m_driver, "gsControlDisplacement::driver: constructed with a static solver, not with a problem");
        return *m_driver;
    }

    /// Performs a step of size \a dL in the load and returns its status. On failure, the load and the displacements are unchanged.
    gsStatus step(T dL)
    {
        if (m_driver)
        {
            // R(U_n, L+dL) ~ dL * dR/dL at equilibrium, so |dL*dR/dL| is, to first order, the
            // legacy initial residual; |R| at the Load-predicted point is round-off for an
            // exact (linear) predictor and cannot be the basis of the relative residual.
            gsVector<T> d;
            if (m_problem->loadDirection)
            {
                if (!m_problem->loadDirection(m_L + dL, d))
                    return m_status = gsStatus::AssemblyError;
            }
            else
                d = m_problem->F;
            m_resTest->setBasisFloor(math::abs(dL) * d.norm());

            m_driver->controller().setLength(dL);
            m_status = m_driver->step();
            if (m_status == gsStatus::Success)
            {
                m_L = m_driver->lambda();
                m_space->exportWorking(m_U);
            }
            return m_status;
        }

        m_solver->setLoad(m_L + dL);
        if (!first)
            m_solver->setDisplacement(m_U);
        else first = false;
        m_status = m_solver->solve();

        if (m_status == gsStatus::Success)
        {
            m_L += dL;
            m_U += m_solver->update();
        }
        return m_status;
    }

    /// Return the displacements
    gsVector<T> & solutionU()
    {
        return m_U;
    }

    /// Return the load
    T solutionL()
    {
        return m_L;
    }

    /// Reset the solver, or the step history of the driver
    void reset()
    {
        if (m_driver)
            m_driver->resetHistory();
        else
            m_solver->reset();
    }

    /**
        @brief Sets the load and the displacements to zero and resets the solver.

        With a problem and a state space, the state space is set to zero only
        if it is copyable(); otherwise it keeps its state.
    */
    void setZero()
    {
        if (m_driver)
        {
            m_L = 0;
            m_driver->setLambda(0);
            m_driver->resetHistory();
            m_U.setZero();
            if (m_space->copyable())
                m_space->importCommitted(m_U);
            return;
        }

        m_solver->reset();
        m_U = gsVector<T>::Zero(m_solver->numDofs());
        m_L = 0;
        m_solver->setLoad(m_L);
    }

protected:
    mutable gsStaticBase<T> * m_solver;
    gsStateSpace<T> * m_space;
    gsContinuationProblem<T> * m_problem;
    memory::shared_ptr< gsStaticResidualTest<T> > m_resTest;
    memory::unique_ptr< gsContinuationDriver<T> > m_driver;
    T m_L;
    bool first;

    gsVector<T> m_U;

    gsStatus m_status;
};


} // namespace gismo
