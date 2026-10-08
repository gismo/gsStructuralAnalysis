 /** @file gsStaticNewton.hpp

    @brief Performs linear modal analysis given a matrix or functions of a matrix

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#pragma once

namespace gismo
{

template <class T>
void gsStaticNewton<T>::defaultOptions()
{
    Base::defaultOptions();
    m_options.setString("Solver","CGDiagonal"); // The CG solver is robust for membrane models, where zero-blocks in the matrix might occur.
    m_options.addReal("Relaxation","Relaxation parameter",1);
};

template <class T>
void gsStaticNewton<T>::getOptions()
{
    Base::getOptions();

    // if (m_solverType!=solver::LDLT && m_stabilityMethod==stabmethod::Determinant)
    // {
    //   gsWarn<<"Determinant method cannot be used with solvers other than LDLT. Bifurcation method will be set to 'Eigenvalue'.\n";
    //   m_stabilityMethod = stabmethod::Eigenvalue;
    // }

    m_relax = m_options.getReal("Relaxation");
};

template <class T>
void gsStaticNewton<T>::initOutput()
{
    gsInfo<<"\t";
    gsInfo<<std::setw(4)<<std::left<<"It.";
    gsInfo<<std::setw(17)<<std::left<<"|R|";
    gsInfo<<std::setw(17)<<std::left<<"|R|/|R0|";
    gsInfo<<std::setw(17)<<std::left<<"|DU|";
    gsInfo<<std::setw(17)<<std::left<<"|dU|";
    gsInfo<<std::setw(17)<<std::left<<"|dU|/|DU|";
    gsInfo<<std::setw(17)<<std::left<<"|dU|/|U+DU|";
    gsInfo<<std::setw(17)<<std::left<<"log(Ri/R0):";
    gsInfo<<std::setw(17)<<std::left<<"log(Ri+1/R0)";
    gsInfo<<"\n";
}

template <class T>
void gsStaticNewton<T>::stepOutput(index_t k)
{
    gsInfo<<"\t";
    gsInfo<<std::setw(4)<<std::left<<k;
    gsInfo<<std::setw(17)<<std::left<<m_residual;
    gsInfo<<std::setw(17)<<std::left<<m_residual/m_residualIni;
    gsInfo<<std::setw(17)<<std::left<<m_DeltaU.norm();
    gsInfo<<std::setw(17)<<std::left<<m_relax * m_deltaU.norm();
    gsInfo<<std::setw(17)<<std::left<<m_relax * m_deltaU.norm()/m_DeltaU.norm();
    gsInfo<<std::setw(17)<<std::left<<m_relax * m_deltaU.norm()/(m_U+m_DeltaU).norm();
    gsInfo<<std::setw(17)<<std::left<<math::log10(m_residualOld/m_residualIni);
    gsInfo<<std::setw(17)<<std::left<<math::log10(m_residual/m_residualIni);
    gsInfo<<"\n";

}

template <class T>
gsStatus gsStaticNewton<T>::solveLinear()
{
    try
    {
        _solveLinear();
        m_status = gsStatus::Success;
    }
    catch (int errorCode)
    {
        if      (errorCode==1)
            m_status = gsStatus::NotConverged;
        else if (errorCode==2)
            m_status = gsStatus::AssemblyError;
        else if (errorCode==3)
            m_status = gsStatus::SolverError;
        else
            m_status = gsStatus::OtherError;
    }
    catch (...)
    {
        m_status = gsStatus::OtherError;
    }
    return m_status;
};

template <class T>
gsVector<T> gsStaticNewton<T>::_solveLinear()
{
    this->getOptions();
    if (m_verbose==2)
    {
        gsInfo<<"Matrix: \n"<<m_linear.toDense()<<"\n";
        gsInfo<<"Vector: \n"<<m_force<<"\n";
    }
    _factorizeMatrix( m_linear );
    m_U = _solveSystem(m_force);
    return m_U;
};

template <class T>
gsStatus gsStaticNewton<T>::solveNonlinear()
{
    // The single try of the nonlinear solve: user operations (assemblers) and the
    // driver setup may throw; the corrector loop itself never does.
    // Whatever the failure, the displacements are those of the call.
    const gsVector<T> U0 = m_U;
    try
    {
        m_status = _solveNonlinear();
    }
    catch (int errorCode)
    {
        if      (errorCode==1)
            m_status = gsStatus::NotConverged;
        else if (errorCode==2)
            m_status = gsStatus::AssemblyError;
        else if (errorCode==3)
            m_status = gsStatus::SolverError;
        else
            m_status = gsStatus::OtherError;
    }
    catch (...)
    {
        m_status = gsStatus::OtherError;
    }
    if (m_status != gsStatus::Success)
        m_U = U0;
    return m_status;
}

/*
   Newton iteration at fixed load as a continuation step with a Fixed constraint:
   U <- U + dU with K dU = R, R = Fext - Fint(U) (the user residual, evaluated at
   the load m_L by the residual functor), K = dFint/dU.

   Acceptance: at least one update, |R|/|R0| < tolF and |dU|/|DU| < tolU with DU
   the running increment of the step (not the predictor, which is zero). With a
   headstart, the predictor is the given DU and R0 is the residual at U.
*/
template <class T>
gsStatus gsStaticNewton<T>::_solveNonlinear()
{
    this->getOptions();
    // m_start: true -> m_U given
    // m_headstart: true -> m_DeltaU given

    if (m_DeltaU.norm()==0 && m_DeltaU.rows()==0)
    {
        m_deltaU = m_DeltaU = this->_solveLinear();
        m_U.setZero(); // Needed because linear solve modifies m_U.
        m_headstart = true; // due to this, the relative residual is based on the solution of the linear solve
    }

    const index_t n = m_U.size();
    m_space.reset(new gsVectorStateSpace<T>(m_U));

    // The user residual Fext − Fint is passed to the driver unchanged (no negation) because the driver solves K δ = R.
    m_problem.residual = [this](T, gsVector<T> & R, gsEvalMode) -> bool
    {
        return m_residualFun(m_space->working(), R);
    };
    m_problem.jacobian = [this](T, gsSparseMatrix<T> & K) -> bool
    {
        const gsStepContext<T> & ctx = m_driver->context();
        // The Jacobian of the first update sees a zero iterative update, later ones the unrelaxed one.
        gsVector<T> dx;
        if (ctx.iteration == 0)
            dx.setZero(m_space->working().size());
        else
        {
            dx = ctx.deltaU;
            if (ctx.relaxation != T(1)) dx /= ctx.relaxation;
        }
        if (!m_dnonlinear(m_space->working(), dx, K))
            return false;
        if (m_verbose==2)
        {
            gsInfo<<"Matrix: \n"<<K.toDense()<<"\n";
            gsInfo<<"Vector: \n"<<ctx.R<<"\n";
        }
        return true;
    };
    m_problem.loadDirection = nullptr;
    m_problem.F = gsVector<T>::Zero(n); // read only by the per-update load tangent solve, trivial for a zero right-hand side

    m_driver.reset(new gsContinuationDriver<T>(m_problem, *m_space));
    if (m_headstart)
        m_driver->setConstraint(typename gsPathConstraint<T>::Ptr(new gsStaticHeadstartConstraint<T>(m_DeltaU)));
    else
        m_driver->setConstraint(typename gsPathConstraint<T>::Ptr(new gsFixedConstraint<T>()));

    typedef typename gsConvergenceTest<T>::Ptr TestPtr;
    m_resTest.reset(new gsStaticResidualTest<T>(m_tolF));
    memory::shared_ptr< gsAllTest<T> > all(new gsAllTest<T>());
    all->add(TestPtr(new gsMinIterTest<T>(1))).add(m_resTest).add(TestPtr(new gsStaticUpdateTest<T>(m_tolU)));
    memory::shared_ptr< gsAnyTest<T> > any(new gsAnyTest<T>());
    any->add(all).add(TestPtr(new gsMaxIterTest<T>(m_maxIterations)));
    m_driver->setTest(any);

    m_driver->corrector().setRelaxation(m_relax);
    m_driver->linearSolver().setSolver(m_options.askString("Solver","SimplicialLDLT"), true);
    m_driver->setLambda(m_L);
    m_driver->setup();

    if (m_headstart)
    {
        // The relative residual of a headstart is based on the residual at U, not at U + DeltaU.
        gsVector<T> R0;
        if (!m_residualFun(m_U, R0))
            return gsStatus::AssemblyError;
        m_resTest->setBasis(R0.norm());
        m_headstart = false;
    }

    if (m_verbose>0)
    {
        initOutput();
        m_driver->addObserver([this](const gsStepEvent<T> & ev)
        {
            if (ev.kind != gsStepEvent<T>::IterationDone)
                return;
            const gsStepContext<T> & ctx = ev.ctx;
            m_residualIni = m_resTest->basis();
            if (ctx.iteration == 1)
                m_residualOld = m_resTest->predictedResidual();
            m_residual = ctx.R.norm();
            m_DeltaU = ctx.DeltaU;
            m_deltaU = ctx.deltaU;
            if (ctx.relaxation != T(1)) m_deltaU /= ctx.relaxation; // stepOutput applies m_relax
            this->stepOutput(ctx.iteration-1);
            m_residualOld = m_residual;
        });
    }

    const gsStatus st = m_driver->step();

    const gsStepContext<T> & ctx = m_driver->context();
    m_DeltaU = ctx.DeltaU;
    m_R = ctx.R;
    m_residual = m_R.norm();

    const index_t k = m_driver->numIterations();
    const gsDriverOp op = m_driver->report().failedOp;
    if (st==gsStatus::Success || st==gsStatus::NotConverged
        || (st==gsStatus::AssemblyError && op==gsDriverOp::Residual && k>0))
        m_numIterations = math::max(k-1,index_t(0));
    else
        m_numIterations = k;

    if (st==gsStatus::Success)
        m_U = m_space->committed();
    else if (st==gsStatus::NotConverged && op==gsDriverOp::None)
        gsInfo<<"Maximum iterations reached. Solution did not converge\n";
    return st;
}

template <class T>
gsVector<T> gsStaticNewton<T>::_computeResidual(const gsVector<T> & U)
{
  gsVector<T> resVec;
  if (!m_residualFun(U, resVec))
    throw 2;
  return resVec;
}

template <class T>
void gsStaticNewton<T>::_factorizeMatrix(const gsSparseMatrix<T> & jacMat) const
{
    m_solver->compute(jacMat);
    if (m_solver->info()!=gsEigen::ComputationInfo::Success)
    {
      gsInfo<<"Solver error with code "<<m_solver->info()<<". See Eigen documentation on ComputationInfo \n"
                                                                    <<gsEigen::ComputationInfo::Success<<": Success"<<"\n"
                                                                    <<gsEigen::ComputationInfo::NumericalIssue<<": NumericalIssue"<<"\n"
                                                                    <<gsEigen::ComputationInfo::NoConvergence<<": NoConvergence"<<"\n"
                                                                    <<gsEigen::ComputationInfo::InvalidInput<<": InvalidInput"<<"\n";
      throw 3;
    }
}

template <class T>
gsVector<T> gsStaticNewton<T>::_solveSystem(const gsVector<T> & F)
{
    try
    {
      return m_solver->solve(F);
    }
    catch (...)
    {
      throw 3;
    }
}

template <class T>
void gsStaticNewton<T>::reset()
{
    m_dofs = m_force.rows();
    // resets m_U, m_DeltaU, m_deltaU, m_R, m_L, m_DeltaL, m_deltaL and m_headstart
    Base::reset();
}

template <class T>
void gsStaticNewton<T>::_init()
{
    this->reset();
    if( m_dnonlinear==nullptr || m_residualFun==nullptr)
        m_NL=false;
    else
        m_NL = true;

    m_stabilityMethod = 0;
    m_start = false;

    m_dofs = m_linear.rows();
    if (m_dofs==0)
        gsWarn<<"The number of degrees of freedom is equal to zero. This can lead to bad initialization.\n";

    m_residual = m_residualIni = m_residualOld = 0;

    defaultOptions();

    m_status = gsStatus::NotStarted;
}

} // namespace gismo
