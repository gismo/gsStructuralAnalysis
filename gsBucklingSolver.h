 /** @file gsBucklingSolver.h

    @brief Performs linear buckling analysis given a matrix or functions of a matrix

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#include <gsStructuralAnalysis/gsEigenProblemBase.h>

#pragma once

namespace gismo
{

/**
    @brief Performs the arc length method to solve a nonlinear equation system.

    \tparam T coefficient type

    \ingroup ThinShell
*/
template <class T, Spectra::GEigsMode GEigsMode = Spectra::GEigsMode::Cholesky>
class gsBucklingSolver : public gsEigenProblemBase<T,GEigsMode>
{
protected:

    typedef gsEigenProblemBase<T,GEigsMode> Base;

public:

    /// Constructor giving access to the gsShellAssembler object to create a linear system per iteration
    gsBucklingSolver(     gsSparseMatrix<T> &linear,
                        gsVector<T> &rhs,
                        std::function < gsSparseMatrix<T> ( gsVector<T> const & ) > &nonlinear,
                        T scaling = 1.0) :
        m_rhs(rhs),
        m_nonlinearFun(nonlinear),
        m_scaling(scaling)
    {
        m_A = linear;
        m_verbose = false;
        this->initializeMatrix();
    }

    /// Constructor giving access to the gsShellAssembler object to create a linear system per iteration
    gsBucklingSolver(   gsSparseMatrix<T> &linear,
                        gsSparseMatrix<T> &nonlinear )
    {
        m_A = linear;
        m_verbose = false;
        m_B = nonlinear-m_A;
    }

protected:

    void initializeMatrix()
    {
        if (m_verbose) { gsInfo<<"Computing matrices" ; }
        m_solver.compute(m_A);
        if (m_verbose) { gsInfo<<"." ; }
        m_solVec = m_solver.solve(m_scaling*m_rhs);
        if (m_verbose) { gsInfo<<"." ; }
        m_B = m_nonlinearFun(m_solVec)-m_A;
        if (m_verbose) { gsInfo<<"." ; }
        if (m_verbose) { gsInfo<<"Finished\n" ; }
    }

protected:

    using Base::m_A;
    const gsVector<T> m_rhs;
    const std::function < gsSparseMatrix<T> ( gsVector<T> const & ) > m_nonlinearFun;
    T m_scaling;
    using Base::m_B;


    /// Linear solver employed
    gsSparseSolver<>::SimplicialLDLT  m_solver;
    gsVector<> m_solVec;

    using Base::m_verbose;

};

} // namespace gismo
