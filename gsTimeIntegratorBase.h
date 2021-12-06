// /** @file gsTimeIntegratorBase.h

//     @brief Provides temporal solvers for structural analysis problems


//     TO DO:
//     * [V] Fix Implicit Euler Nonlinear
//     * [ ] Fix Implicit Euler BlockOp
//     * [V] Fix Explicit Euler with nonlinearities (NOT VERIFIED)
//     * [V] Fix linear Newmark and Bathe
//     * [ ] Fix Implicit/Explicit Linear euler

//     This file is part of the G+Smo library.

//     This Source Code Form is subject to the terms of the Mozilla Public
//     License, v. 2.0. If a copy of the MPL was not distributed with this
//     file, You can obtain one at http://mozilla.org/MPL/2.0/.

//     Author(s): H.M. Verhelst (2019-..., TU Delft)
// */


// #pragma once


// #include <gsSolver/gsBlockOp.h>
// #include <gsSolver/gsMatrixOp.h>

// namespace gismo
// {

// enum class TMethod : short_t
// {
//     ExplEuler = 0;
//     ImplEuler = 1,
//     Newmark = 2,
//     Bathe = 3,
//     RK4 = 4,
// }

// /**
//     @brief  ...

//     \tparam T coefficient type

//     \ingroup ThinShell
// */
// template <class T>
// class gsTimeIntegratorBase
// {
// public:

//     void step();

//     /// Tells if the Newton method converged
//     bool converged() const {return m_converged;}

//     void setTime(T time) {m_t = time; }
//     T currentTime() const {return m_t; }

//     /// set mass matrix
//     void setMassMatrix(gsMatrix<T>& Mass) {m_mass = Mass; }
//     void setDampingMatrix(gsMatrix<T>& Damp) {m_damp = Damp; }
//     void setStiffnessMatrix(gsMatrix<T>& Stif) {m_stif = Stif; }
//     void setJacobian(std::function < gsSparseMatrix<T> ( gsMatrix<T> const & ) > &Jacobian) {m_jacobian = Jacobian; }

//     // set solutions
//     void setDisplacement(gsMatrix<T>& displ);
//     void setVelocity(gsMatrix<T>& velo);
//     void setAcceleration(gsMatrix<T>& accel);

//     const gsMatrix<T>& displacements() const { return m_displacements;}
//     const gsMatrix<T>& velocities() const { return m_velocities;}
//     const gsMatrix<T>& accelerations() const { return m_accelerations;}

//     void constructSolution();

//     void resetIteration();

//     gsOptionList & options() {return m_options;}

//     void applyOptions() {this->getOptions(); }
//     void resetOptions() {this->defaultOptions(); this->getOptions(); }

//     void printSettings() { gsInfo<<m_options<<"\n"; }

// protected:
//     index_t m_maxIterations;
//     T m_tolerance;
//     bool m_verbose;
//     bool m_quasiNewton;

//     T m_t;

//     gsSparseMatrix<T> m_mass;
//     gsSparseMatrix<T> m_damp;
//     gsSparseMatrix<T> m_stif;
//     std::function < gsSparseMatrix<T> ( gsMatrix<T> const & ) > m_jacobian;
//     std::function < gsMatrix<T> ( gsMatrix<T> const &, T ) > m_residualFun;
//     std::function < gsMatrix<T> ( T) > m_forceFun;
//     gsSparseMatrix<T> m_jacMat;
//     gsMatrix<T> m_resVec;
// }

// } //namespace

// #ifndef GISMO_BUILD_LIB
// #include GISMO_HPP_HEADER(gsTimeIntegratorBase.hpp)
// #endif
