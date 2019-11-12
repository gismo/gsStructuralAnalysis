 /** @file gsArcLengthIterator.h

    @brief Performs the arc length method to solve a nonlinear equation system.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <typeinfo>
#pragma once


namespace gismo
{

/**
    @brief Performs the arc length method to solve a nonlinear equation system.

    \tparam T coefficient type

    \ingroup ThinShell
*/
template <class T>
class gsArcLengthIterator
{
public:

    /// Constructor giving access to the gsShellAssembler object to create a linear system per iteration
    gsArcLengthIterator(  std::function < gsSparseMatrix<T> ( gsVector<T> const & ) > &Jacobian,
                          std::function < gsVector<T> ( gsVector<T> const &, T, gsVector<T> const & ) > &Residual,
                          gsVector<T> &Force )
    : m_jacobian(Jacobian),
      m_residualFun(Residual),
      m_forcing(Force),
      m_numIterations(0),
      m_maxIterations(100),
      m_arcLength(1e-2),
      m_arcLength_prev(1e-2),
      m_tolerance(1e-6),
      m_toleranceF(1e-3),
      m_toleranceU(1e-6),
      m_converged(false),
      m_tau(1e3),
      m_verbose(false)
    {
      // gsInfo<<".";
      m_bifurcationMethod = "Eigenvalue";
      m_bifurcation = false;
      m_angleDetermine = 0; // Angle determination method: 0: determine based on previous load step. 1: determine based on previous iteration ONLY FOR CRISFIELD'S METHOD

      m_adaptiveLength = false;
      m_quasiNewton = false;
      m_quasiNewtonInterval = -1;

      m_relax = 1.0;

      m_method = "Crisfield";

      // initialize errors
      m_basisResidualF = 0.0;
      m_basisResidualU = 0.0;

      this -> initialize();
      // m_deltaLs = gsMatrix<T>::Zero(2,1);
    }


public:

    void initialize();

    void step();

    void computeSingularPoint();
    void computeSingularPoint(T tol);

    void computeSingularPoint(const gsVector<T> U, T L);
    void computeSingularPoint(const gsVector<T> U, T L, T tol);

    void checkBifurcation();

    void switchBranch();

    void bisectionSolve();
    void bisectionSolve(T tol);

    void bisectionSolve(const gsVector<T> U, T L);
    void bisectionSolve(const gsVector<T> U, T L, T tol);


public:

    // Control parameters

      /// Set Arclength method
      void setMethod(std::string method) {m_method = method; this -> initialize();}

      /// Set the maximum number of Newton iterations allowed
      void setMaxIterations(index_t nIter) {m_maxIterations = nIter;}

      /// Set the tolerance for convergence
      void setTolerance(T tol) {m_tolerance = tol;}
      void setToleranceF(T tol) {m_toleranceF = tol;}
      void setToleranceU(T tol) {m_toleranceU = tol;}

      /// Set mode shape parameter tau
      void setTau(T tau) {m_tau = tau;}

      /// Set arc length
      void setLength(T length) {m_arcLength_prev = m_arcLength; m_arcLength = length;}
      void setLength(T length, bool adaptive)
      {
        m_arcLength_prev = m_arcLength;
        m_arcLength = length;
        m_adaptiveLength = adaptive;
        m_desiredIterations = 10;
      } // number of desired iterations defaults to 10
      void setLength(T length, bool adaptive, index_t iterations)
      {
        m_arcLength_prev = m_arcLength;
        m_arcLength = length;
        m_adaptiveLength = adaptive;
        m_desiredIterations = iterations;
      }
      void setLength(T length, index_t iterations)
      {
        m_arcLength_prev = m_arcLength;
        m_arcLength = length;
        m_adaptiveLength = true;
        m_desiredIterations = iterations;
      }

      // Use quasi-newton or not
      void quasiNewton() {m_quasiNewton = true; };
      void quasiNewton(int interval) {m_quasiNewton = true; m_quasiNewtonInterval = interval;}
    // Methods and procedures
      /// Set method for singular point detection
      void setBifurcationMethod(std::string method) {m_bifurcationMethod = method;}

      /// Set AngleDetermination
      void setAngleDeterminationMethod(int method) {m_angleDetermine = method;}

    // Toggles
      // Verbose output
      void verbose() {m_verbose = true;}

    // Output
      /// Tells if the Arc Length method converged
      bool converged() const {return m_converged;}

      /// Returns the number of Newton iterations performed
      index_t numIterations() const { return m_numIterations;}

      /// Returns the tolerance value used
      T tolerance() const {return m_tolerance;}

      /// Returns the error after solving the nonlinear system
      T residue()   const {return m_residue;}

      /// Returns the value of the Determinant or other indicator
      T indicator() const {return m_indicator;}

      /// Return the solution vector and factor
      const gsVector<T> & solutionU() const {return  m_U;}
      const gsVector<T> & solutionDU() const {return  m_DeltaU;}
      T  solutionL() const {return  m_L;}
      T  solutionDL() const {return  m_DeltaL;}
      const gsVector<T> & solutionV() const {return  m_V;}

      // Returns if solution passed a bifurcation point
      bool isBifurcation() const {return m_bifurcation;}

      T determinant() const {return m_jacobian(m_U).toDense().determinant();}

    // Miscelaneous
      void resetStep() {m_DeltaUold.setZero();}

      // Set initial guess for solution
      void setInitialGuess(const gsVector<T> guess) {m_U = guess;}

      void setSolution(const gsVector<T> U, T L) {m_L = L; m_U = U; m_DeltaUold.setZero(); }

      void setRelaxation(T relaxation) {m_relax = relaxation; }

      void printSettings();

protected:

    void computeLambdas();
    void computeLambdasSimple();
    void computeLambdasModified();
    void computeLambdasComplex();
    void computeLambdasEta();

    void computeLambdaDET();

    void computeLambdaDOT();

    void computeLambdaMU();

    void extendedSystemIteration();

    int bisectionObjectiveFunction(const gsVector<T> x);

    T bisectionTerminationFunction(const gsVector<T> x);

    void computeResidual();
    void computeResidualNorms();
    void computeJacobian();
    void computeLength();
    void computeUt();
    void computeUbar();

    void initOutput();
    void initOutputRiks();
    void initOutputCrisfield();
    void initOutputExtended();

    void stepOutput();
    void stepOutputRiks();
    void stepOutputCrisfield();
    void stepOutputExtended();

    void iteration();
    void iterationRiks();
    void iterationConsistentCrisfield();
    void iterationExplicitIterations();
    void iterationCrisfield();

    void initiateStep();
    void initiateStepRiks();
    void initiateStepConsistentCrisfield();
    void initiateStepExplicitIterations();
    void initiateStepCrisfield();

    void predictor();
    void predictorRiks();
    void predictorConsistentCrisfield();
    void predictorExplicitIterations();
    void predictorCrisfield();

    void iterationFinish();
    void iterationFinishRiks();
    void iterationFinishConsistentCrisfield();
    void iterationFinishExplicitIterations();
    void iterationFinishCrisfield();

protected:

    // Number of degrees of freedom
    index_t m_numDof;

    std::function < gsSparseMatrix<T> ( gsVector<T> const & ) > m_jacobian;
    std::function < gsMatrix<T> ( gsVector<T> const &, T, gsVector<T> const & ) > m_residualFun;
    std::function < gsMatrix<T> ( gsVector<T> const &, T, gsVector<T> const & ) > m_residualFunMod;
    gsVector<T> m_forcing;

    /// Linear solver employed
    gsSparseSolver<>::SimplicialLDLT  m_solver;

protected:

    /// Number of Arc Length iterations performed
    index_t m_numIterations;

    /// Maximum number of Arc Length iterations allowed
    index_t m_maxIterations;

    /// Number of desired iterations
    index_t m_desiredIterations;

    /// Length of the step in the u,f plane
    T m_arcLength;
    T m_arcLength_prev;
    bool m_adaptiveLength;

    /// Tuning parameter
    T m_phi;

    /// Tolerance value to decide convergence
    T m_tolerance;

    /// Tolerance value to decide convergence - Force criterion
    T m_toleranceF;

    /// Tolerance value to decide convergence - Displacement criterion
    T m_toleranceU;

    bool m_verbose;

    bool m_quasiNewton;
    int m_quasiNewtonInterval;

    std::string note;

protected:

    /// Convergence result
    bool m_converged;

    /// Force residuum
    T m_residue;

    /// Force residuum
    T m_residueF;
    T m_basisResidualF;

    /// Displacement residuum
    T m_residueU;
    T m_basisResidualU;

    /// Load residuum
    T m_residueL;

    /// Singular point
    T m_residueKTPhi;
    T m_basisResidualKTPhi;

    /// Indicator for bifurcation
    T m_indicator;

    /// Relaxation factor
    T m_relax;

protected:

    // Arc length method (either "Crisfield" or "Riks")
    std::string m_method;

    // Previous update
    gsVector<T> m_DeltaUold;
    real_t m_DeltaLold;
    /// Displacement vector (present, at previously converged point)
    gsVector<T> m_U, m_Uprev;
    /// Update of displacement vector
    gsVector<T> m_DeltaU;
    /// u_bar
    gsVector<T> m_deltaUbar;
    /// u_t
    gsVector<T> m_deltaUt;
    /// Update of update of displacement vector
    gsVector<T> m_deltaU;

    /// Lambda (present, at previously converged point)
    T m_L, m_Lprev;
    /// Update of lambdaGeneralizedSelfAdjointEigenSolver
    T m_DeltaL;
    /// Update of update of lambda
    T m_deltaL;
    /// Vector with lambda updates
    gsVector<T> m_deltaLs;

    /// Jacobian matrix
    gsSparseMatrix<T> m_jacMat;
    T m_detKT;

    /// Value of residual function
    gsVector<T> m_resVec;

    // eigenvector
    gsVector<T> m_V;
    // step eigenvector
    gsVector<T> m_DeltaV;
    gsVector<T> m_DeltaVbar;
    gsVector<T> m_DeltaVt;

    // Boolean if previous step contained was bifurcation point
    bool m_bifurcation;
    // Method to check if a point is a bifurcation point
    std::string m_bifurcationMethod;

    // Angle determination method: 0: determine based on previous load step. 1: determine based on previous iteration
    int m_angleDetermine;

    // Branch switch parameter
    T m_tau;

    // MODIFIED ARC LENGTH METHOD
    /// factor (modified arc length method)
    T m_eta;

    // discriminant
    T m_discriminant;

    T m_alpha1;
    T m_alpha2;
    T m_alpha3;

    T m_a0;
    T m_b0,m_b1;
    T m_c0,m_c1,m_c2;

    T m_tangentLength;
};


} // namespace gismo

//////////////////////////////////////////////////
//////////////////////////////////////////////////

namespace gismo
{

// Miscelaneous functions
/// sign function
template <class T>
int sign(T val)
{
    return (T(0) < val) - (val < T(0));
}

/// Modulus function
template <class T>
int mod(T x, T y)
{
    T m = x - floor(x / y) * y;
    m *= sign(y);
    return m;
}

/// sort vector
template <class T>
int countNegatives(gsVector<T> vec)
{
  int count = 0;
  index_t N = vec.cols();
  index_t M = vec.rows();
  for(int i = 0; i < M; i++)
    {
        for(int j = 0; j < N; j++)
        {
            if( vec(i,j) < 0 )
                count += 1;
        }
    }
    return count;
}

template <class T>
void gsArcLengthIterator<T>::initialize()
{
  m_numDof = m_forcing.size();
  m_U = gsVector<T>::Zero(m_numDof);
  m_L = 0.0;

  if ((m_method == "Crisfield"))
  {
    m_DeltaUold = gsVector<T>::Zero(m_numDof);
    m_DeltaLold = 0.0;
  }
  if ((m_method == "Riks") || (m_method == "ConsistentCrisfield") || (m_method == "ExplicitIterations"))
  {
    m_Uprev = gsVector<T>::Zero(m_numDof);
    m_Lprev = 0.0;
  }
}

template <class T>
void gsArcLengthIterator<T>::computeLength()
{
  // m_desiredIterations := optimal number of iterations
  T fac = m_desiredIterations / m_numIterations;
  if (fac < 0.5)
    fac = 0.5;
  else if (fac > 2.0)
    fac = 2.0;

  m_arcLength_prev = m_arcLength;
  m_arcLength = m_arcLength*fac;
}

template <class T>
void gsArcLengthIterator<T>::computeResidual()
{
  m_resVec = m_residualFun(m_U + m_DeltaU, m_L + m_DeltaL, m_forcing);
}

template <class T>
void gsArcLengthIterator<T>::computeResidualNorms()
{
  if (m_numIterations ==0 ) // then define the residual
  {
    m_basisResidualF = m_resVec.norm();
    m_basisResidualU = m_deltaU.norm();
    m_residueF = m_basisResidualF;
    m_residueU = m_basisResidualU;
    m_residueL = m_DeltaL;
  }
  else
  {
    m_residueF = m_resVec.norm() / m_basisResidualF;
    m_residueU = m_deltaU.norm() / m_basisResidualU;
    m_residueL = m_deltaL / m_DeltaL;  
    // gsInfo<<"basis residue ="<<m_basisResidualF<<"\t residue ="<<m_residueF<<"\t resvec norm ="<<m_resVec.norm()<<"\tCTD"<<"\n";
  }
}

template <class T>
void gsArcLengthIterator<T>::computeJacobian()
{
  // Compute Jacobian
  m_jacMat = m_jacobian(m_U + m_DeltaU);
  m_solver.compute(m_jacMat);
  note += "J";
}

template <class T>
void gsArcLengthIterator<T>::iteration()
{
  if (m_method=="Riks")
    iterationRiks();
  else if (m_method=="ConsistentCrisfield")
    iterationConsistentCrisfield();
  else if (m_method=="ExplicitIterations")
    iterationExplicitIterations();
  else if (m_method=="Crisfield")
    iterationCrisfield();
  else
  {
    gsInfo<<"Error: Method unknown...\n Terminating process...\n";
    std::terminate();
  }
}

template <class T> void gsArcLengthIterator<T>::computeUbar() { m_deltaUbar = m_solver.solve(-m_resVec); }
template <class T> void gsArcLengthIterator<T>::computeUt()   { m_deltaUt = m_solver.solve(m_forcing); }

template <class T>
void gsArcLengthIterator<T>::initiateStep()
{
  m_converged = false;
  m_numIterations = 0;
  if (m_method=="Riks")
    initiateStepRiks();
  else if (m_method=="ConsistentCrisfield")
    initiateStepConsistentCrisfield();
  else if (m_method=="ExplicitIterations")
    initiateStepExplicitIterations();
  else if (m_method=="Crisfield")
    initiateStepCrisfield();
  else
  {
    gsInfo<<"Error: Method unknown...\n Terminating process...\n";
    std::terminate();
  }
}

template <class T>
void gsArcLengthIterator<T>::predictor()
{
  if (m_method=="Riks")
    predictorRiks();
  else if (m_method=="ConsistentCrisfield")
    predictorConsistentCrisfield();
  else if (m_method=="ExplicitIterations")
    predictorExplicitIterations();
  else if (m_method=="Crisfield")
    predictorCrisfield();
  else
  {
    gsInfo<<"Error: Method unknown...\n Terminating process...\n";
    std::terminate();
  }
}

template <class T>
void gsArcLengthIterator<T>::iterationFinish()
{
  if (m_method=="Riks")
    iterationFinishRiks();
  else if (m_method=="ConsistentCrisfield")
    iterationFinishConsistentCrisfield();
  else if (m_method=="ExplicitIterations")
    iterationFinishExplicitIterations();
  else if (m_method=="Crisfield")
    iterationFinishCrisfield();
  else
  {
    gsInfo<<"Error: Method unknown...\n Terminating process...\n";
    std::terminate();
  }
}

template <class T>
void gsArcLengthIterator<T>::step()
{

  initiateStep();
  computeJacobian();
  predictor();
  computeResidual();
  computeResidualNorms();

  if (m_verbose)
     stepOutput();

  if (m_quasiNewton)
  {
    computeJacobian();
    computeUt(); // rhs does not depend on solution
    computeUbar(); // rhs contains residual and should be computed every time
  }

  for (m_numIterations = 1; m_numIterations < m_maxIterations; ++m_numIterations)
  {
    if ( (!m_quasiNewton) || ( ( m_quasiNewtonInterval>0 ) && ( mod(m_numIterations,m_quasiNewtonInterval) < 1e-10 ) ) )
    {
      computeJacobian();
      computeUt(); // rhs does not depend on solution
    }

    computeUbar(); // rhs contains residual and should be computed every time
    iteration();
    computeResidual();
    computeResidualNorms();

    if (m_verbose)
       stepOutput();
    else
    {
      gsInfo<<"Residual: "<<m_residue<<"\n";
    }

    // Termination criteria
    if ( m_residueF < m_toleranceF && m_residueU < m_toleranceU )
    // if ( m_residue < m_tolerance )
    {
      iterationFinish();
      // Change arc length
      if (m_adaptiveLength)
        computeLength();
      else
        m_arcLength_prev = m_arcLength;
      break;
    }
  }
}

// ------------------------------------------------------------------------------------------------------------
// ---------------------------------------Riks's method--------------------------------------------------------
// ------------------------------------------------------------------------------------------------------------

// template <class T>
// void gsArcLengthIterator<T>::iterationRiks()
// {
//   m_deltaUbar = m_solver.solve(-m_resVec);
//   m_deltaUt = m_solver.solve(m_forcing); // note the minus!!

//   m_deltaL = - ( (m_DeltaU).dot(m_deltaUbar) ) / ( (m_DeltaL)*m_forcing.dot(m_forcing) + (m_DeltaU).dot(m_deltaUt) );
//   m_deltaU = m_deltaUbar + m_deltaL*m_deltaUt;

//   m_DeltaU += m_deltaU;
//   m_DeltaL += m_deltaL;
// }

template <class T>
void gsArcLengthIterator<T>::iterationRiks()
{
  // Residual function
  // r = m_phi * m_DeltaU.dot(m_DeltaU)  + (1.0-m_phi) * m_DeltaL*m_DeltaL - m_arcLength*m_arcLength;
  T r = m_phi*(m_U + m_DeltaU - m_U).dot(m_U + m_DeltaU - m_U) + (1-m_phi)*math::pow(m_L + m_DeltaL - m_L,2.0) - m_arcLength*m_arcLength;

  m_deltaL = - ( r + 2*m_phi*(m_DeltaU).dot(m_deltaUbar) ) / ( 2*(1-m_phi)*(m_DeltaL) + 2*m_phi*(m_DeltaU).dot(m_deltaUt) );
  m_deltaU = m_deltaUbar + m_deltaL*m_deltaUt;

  m_DeltaU += m_deltaU;
  m_DeltaL += m_deltaL;
}

template <class T>
void gsArcLengthIterator<T>::initiateStepRiks()
{
  // (m_U,m_L) is the present solution (iteratively updated)
  // (m_Uprev,m_Lprev) is the previously converged solution before (m_Lprev,m_Uprev)

  // Reset step
  m_DeltaU = m_deltaU =  gsVector<T>::Zero(m_numDof);
  m_DeltaL = m_deltaL = 0.0;

  // Initiate verbose
  if (m_verbose)
    initOutputRiks();
}

template <class T>
void gsArcLengthIterator<T>::predictorRiks()
{
  // Define scaling
  if (m_numDof ==1)
    m_phi = 0.99999999999;
  else
    m_phi = 1./m_numDof;

  // Check if the solution on start and prev are similar.
  // Then compute predictor of the method
  T tol = 1e-10;

  if ( ((m_U-m_Uprev).norm() < tol) && ((m_L - m_Lprev) * (m_L - m_Lprev) < tol ) )
  {
    note+= "predictor\t";
    T DL = 1.;
    m_deltaUt = m_solver.solve(m_forcing);
    m_deltaU = m_deltaUt / math::sqrt( m_deltaUt.dot(m_deltaUt) + m_DeltaL*DL );
    m_deltaL = DL / math::sqrt( m_deltaUt.dot(m_deltaUt) + m_DeltaL*DL );
  }
  else
  {
    m_deltaL = 1./m_arcLength_prev*(m_L - m_Lprev);
    m_deltaU = 1./m_arcLength_prev*(m_U - m_Uprev);
  }

  // Update iterative step
  m_deltaL *= m_arcLength;
  m_deltaU *= m_arcLength;

  // Update load step
  m_DeltaU += m_deltaU;
  m_DeltaL += m_deltaL;
}

// template <class T>
// void gsArcLengthIterator<T>::predictorRiks()
// {
//   // Check if the solution on start and prev are similar.
//   // Then compute predictor of the method
//   T tol = 1e-10;

//   if ( ((m_U-m_Uprev).norm() < tol) && ((m_L - m_Lprev) * (m_L - m_Lprev) < tol ) )
//   {
//     note+= "predictor\t";
//     T DL = 1.;
//     computeJacobian();
//     m_deltaUt = m_solver.solve(m_forcing);
//     m_deltaU = m_deltaUt / math::sqrt( m_deltaUt.dot(m_deltaUt) + m_DeltaL*DL );
//     m_deltaL = DL / math::sqrt( m_deltaUt.dot(m_deltaUt) + m_DeltaL*DL );
//   }
//   else
//   {
//     m_deltaL = 1./m_arcLength*(m_L - m_Lprev);
//     m_deltaU = 1./m_arcLength*(m_U - m_Uprev);
//   }

//   // Update iterative step
//   m_deltaL *= m_arcLength;
//   m_deltaU *= m_arcLength;

//   // Update load step
//   m_DeltaU += m_deltaU;
//   m_DeltaL += m_deltaL;
// }

template <class T>
void gsArcLengthIterator<T>::iterationFinishRiks()
{
  m_converged = true;
  m_Uprev = m_U;
  m_Lprev = m_L;
  m_U += m_DeltaU;
  m_L += m_DeltaL;
}

// ------------------------------------------------------------------------------------------------------------
// ---------------------------------------Consistent Crisfield's method----------------------------------------
// ------------------------------------------------------------------------------------------------------------

template <class T>
void gsArcLengthIterator<T>::iterationConsistentCrisfield() // see Carrera1994 eqs 20 and 21
{
  T cres = m_DeltaU.dot(m_DeltaU) + m_phi*m_phi*m_DeltaL*m_DeltaL * m_forcing.dot(m_forcing) - m_arcLength*m_arcLength;
  T num = cres + (2*m_DeltaU).dot(m_deltaUbar);
  T denum = (2*m_DeltaU).dot(m_deltaUt) + m_phi*m_phi*2*m_DeltaL*m_forcing.dot(m_forcing);
  m_deltaL = - num / denum;
  m_deltaU = m_deltaL * m_deltaUt + m_deltaUbar;

  m_DeltaU += m_deltaU;
  m_DeltaL += m_deltaL;
}

template <class T>
void gsArcLengthIterator<T>::initiateStepConsistentCrisfield()
{
  // (m_U,m_L) is the present solution (iteratively updated)
  // (m_Uprev,m_Lprev) is the previously converged solution before (m_Lprev,m_Uprev)

  // Reset step
  m_DeltaU = m_deltaU =  gsVector<T>::Zero(m_numDof);
  m_DeltaL = m_deltaL = 0.0;

  // Initiate verbose
  if (m_verbose)
    initOutputRiks();
}

template <class T>
void gsArcLengthIterator<T>::predictorConsistentCrisfield()
{
  // Check if the solution on start and prev are similar.
  // Then compute predictor of the method
  T tol = 1e-10;

  if ( ((m_U-m_Uprev).norm() < tol) && ((m_L - m_Lprev) * (m_L - m_Lprev) < tol ) )
  {
    note+= "predictor\t";
    T DL = 1.;
    m_deltaUt = m_solver.solve(m_forcing);
    m_deltaU = m_deltaUt / math::sqrt( m_deltaUt.dot(m_deltaUt) + m_DeltaL*DL );
    m_deltaL = DL / math::sqrt( m_deltaUt.dot(m_deltaUt) + m_DeltaL*DL );

    // m_phi = 0;
    m_phi = 1;
    // m_phi = math::pow( m_deltaUt.dot(m_deltaUt) / m_forcing.dot(m_forcing),0.5);
  }
  else
  {
    m_deltaL = 1./m_arcLength_prev*(m_L - m_Lprev);
    m_deltaU = 1./m_arcLength_prev*(m_U - m_Uprev);

    // m_phi = 0;
    m_phi = 1;
    // m_phi = math::pow( m_deltaUt.dot(m_deltaUt) / m_forcing.dot(m_forcing),0.5);
  }

  // Update iterative step
  m_deltaL *= m_arcLength;
  m_deltaU *= m_arcLength;

  // Update load step
  m_DeltaU += m_deltaU;
  m_DeltaL += m_deltaL;
}

template <class T>
void gsArcLengthIterator<T>::iterationFinishConsistentCrisfield()
{
  m_converged = true;
  m_Uprev = m_U;
  m_Lprev = m_L;
  m_U += m_DeltaU;
  m_L += m_DeltaL;
  m_DeltaUold = m_DeltaU;
  m_DeltaLold = m_DeltaL;
}

// ------------------------------------------------------------------------------------------------------------
// ---------------------------------------ExplicitIterations Method--------------------------------------------
// ------------------------------------------------------------------------------------------------------------

template <class T>
void gsArcLengthIterator<T>::iterationExplicitIterations() // see Forde 1987 page 627
{
  // first do the Riks-Wempner-Tangent Method (Carrera1994 section 4.1)
  T num = (m_DeltaU).dot(m_deltaUbar);
  T denum = (m_DeltaU).dot(m_deltaUt) + m_phi*m_phi*m_DeltaL*m_forcing.dot(m_forcing);
  m_deltaL = - num / denum;
  m_deltaU = m_deltaL * m_deltaUt + m_deltaUbar;

  // find the length of the tangent in the potential configuration
  m_tangentLength = math::sqrt( m_tangentLength * m_tangentLength + m_deltaU.dot(m_deltaU) + m_phi*m_phi*m_deltaL*m_deltaL );

  // residual of tangents
  T res = ( m_arcLength*m_arcLength / ( m_tangentLength*m_tangentLength ) ) * ( m_tangentLength - m_arcLength );

  m_deltaL = - ( res + m_DeltaU.dot(m_deltaUbar) ) / (m_DeltaU).dot(m_deltaUt) + m_phi*m_phi*m_DeltaL*m_forcing.dot(m_forcing);
  m_deltaU = m_deltaL * m_deltaUt + m_deltaUbar;

  m_DeltaU += m_deltaU;
  m_DeltaL += m_deltaL;
}

template <class T>
void gsArcLengthIterator<T>::initiateStepExplicitIterations()
{
  // (m_U,m_L) is the present solution (iteratively updated)
  // (m_Uprev,m_Lprev) is the previously converged solution before (m_Lprev,m_Uprev)

  // Reset step
  m_DeltaU = m_deltaU =  gsVector<T>::Zero(m_numDof);
  m_DeltaL = m_deltaL = 0.0;

  // Initiate verbose
  if (m_verbose)
    initOutputRiks();
}

template <class T>
void gsArcLengthIterator<T>::predictorExplicitIterations()
{
  // Check if the solution on start and prev are similar.
  // Then compute predictor of the method
  T tol = 1e-10;

  if ( ((m_U-m_Uprev).norm() < tol) && ((m_L - m_Lprev) * (m_L - m_Lprev) < tol ) )
  {
    note+= "predictor\t";
    T DL = 1.;
    m_deltaUt = m_solver.solve(m_forcing);
    m_deltaU = m_deltaUt / math::sqrt( m_deltaUt.dot(m_deltaUt) + m_DeltaL*DL );
    m_deltaL = DL / math::sqrt( m_deltaUt.dot(m_deltaUt) + m_DeltaL*DL );

    // m_phi = 0;
    m_phi = 1;
    // m_phi = math::pow( m_deltaUt.dot(m_deltaUt) / m_forcing.dot(m_forcing),0.5);
  }
  else
  {
    m_deltaL = 1./m_arcLength_prev*(m_L - m_Lprev);
    m_deltaU = 1./m_arcLength_prev*(m_U - m_Uprev);

    // m_phi = 0;
    m_phi = 1;
    // m_phi = math::pow( m_deltaUt.dot(m_deltaUt) / m_forcing.dot(m_forcing),0.5);
  }

  // Update iterative step
  m_deltaL *= m_arcLength;
  m_deltaU *= m_arcLength;

  // Update load step
  m_DeltaU += m_deltaU;
  m_DeltaL += m_deltaL;

  m_tangentLength = math::sqrt( m_DeltaU.dot(m_DeltaU) + m_phi*m_phi*m_DeltaL*m_DeltaL );
}

template <class T>
void gsArcLengthIterator<T>::iterationFinishExplicitIterations()
{
  m_converged = true;
  m_Uprev = m_U;
  m_Lprev = m_L;
  m_U += m_DeltaU;
  m_L += m_DeltaL;
  m_DeltaUold = m_DeltaU;
  m_DeltaLold = m_DeltaL;
}

// ------------------------------------------------------------------------------------------------------------
// ---------------------------------------Crisfield's method---------------------------------------------------
// ------------------------------------------------------------------------------------------------------------
template <class T>
void gsArcLengthIterator<T>::iterationCrisfield()
{
  m_eta = 1.0;

  T lamold = m_deltaL;
  computeLambdas();

  // Relaxation against oscillating load factor
  if (( (lamold*m_deltaL < 0) && (abs(m_deltaL) <= abs(lamold) ) ) && m_relax != 1.0 )
  {
     note += "\t relaxated solution!";
     m_deltaU = m_relax * (m_deltaL*m_deltaUt + m_eta*m_deltaUbar);
     m_deltaL = m_relax * m_deltaL;
  }

  m_DeltaU += m_deltaU;
  m_DeltaL += m_deltaL;

  if (m_angleDetermine == 1)
  {
    m_DeltaUold = m_DeltaU;
    m_DeltaLold = m_DeltaL;
  }
}

template <class T>
void gsArcLengthIterator<T>::initiateStepCrisfield()
{
  if (m_verbose)
    initOutputCrisfield();

  m_DeltaU = m_deltaUbar = m_deltaUt = gsVector<T>::Zero(m_numDof);
  m_DeltaL = m_deltaL = 0.0;
  m_eta = 1.0;
}

template <class T>
void gsArcLengthIterator<T>::predictorCrisfield()
{
  m_deltaUt = m_solver.solve(m_forcing);

  // Choose Solution
  if (m_DeltaUold.dot(m_DeltaUold) == 0 && m_DeltaLold*m_DeltaLold == 0) // no information about previous step.
  // {
  //  m_deltaL = m_arcLength / math::pow(2*( m_deltaUt.dot(m_deltaUt) ) , 0.5);
  //  m_deltaU = m_deltaUbar + m_deltaL*m_deltaUt;
  //  m_phi = math::pow( m_deltaUt.dot(m_deltaUt) / m_forcing.dot(m_forcing),0.5);

  // }
  {
    note+= "predictor\t";
    T DL = 1.;
    m_deltaL = m_arcLength * DL / math::sqrt( 2*( m_deltaUt.dot(m_deltaUt) + m_DeltaL*DL ) );
    // m_deltaU = m_arcLength * m_deltaUt / math::sqrt( m_deltaUt.dot(m_deltaUt) + m_DeltaL*DL );
    m_deltaU = m_deltaL*m_deltaUt;

    // m_phi = 0;
    // m_phi = 1;
    m_phi = math::pow( m_deltaUt.dot(m_deltaUt) / m_forcing.dot(m_forcing),0.5);

    m_DeltaUold = m_deltaU;
    m_DeltaLold = m_deltaL;
  }
  else // previous point is not in the origin
  {
    // m_phi = 0;
    // m_phi = 1;
    m_phi = math::pow(m_U.dot(m_U)/( math::pow(m_L,2) * m_forcing.dot(m_forcing) ),0.5);
    note += " phi=" + std::to_string(m_phi);
    computeLambdaMU();
  }

  // Compute Temporary updates of DeltaL and DeltaU
  m_DeltaU += m_deltaU;
  m_DeltaL += m_deltaL;

  if (m_angleDetermine == 1)
  {
   m_DeltaUold = m_DeltaU;
   m_DeltaLold = m_DeltaL;
  }
}

template <class T>
void gsArcLengthIterator<T>::iterationFinishCrisfield()
{
  m_converged = true;
  m_U += m_DeltaU;
  m_L += m_DeltaL;
  if (m_angleDetermine == 0)
  {
    m_DeltaUold = m_DeltaU;
    m_DeltaLold = m_DeltaL;
  }
}

template <class T>
void gsArcLengthIterator<T>::computeLambdasSimple() //Ritto-Corrêa et al. 2008
{
  T A0 = math::pow(m_phi,2)* m_forcing.dot(m_forcing); // see Lam et al. 1991,

  m_a0 = m_deltaUt.dot(m_deltaUt) + A0;
  m_b0 = 2*( m_deltaUt.dot(m_DeltaU) + m_DeltaL * A0 );
  m_b1 = 2*( m_deltaUbar.dot(m_deltaUt) );
  m_c0 = m_DeltaU.dot(m_DeltaU) + m_DeltaL*m_DeltaL * A0 - math::pow(m_arcLength,2);
  m_c1 = 2*( m_DeltaU.dot(m_deltaUbar) );
  m_c2 = m_deltaUbar.dot(m_deltaUbar);

  /// Calculate the coefficients of the polynomial
  m_alpha1 = m_a0;
  m_alpha2 = m_b0 + m_eta*m_b1;
  m_alpha3 = m_c0 + m_eta*m_c1 + m_eta*m_eta*m_c2;

  m_discriminant = math::pow(m_alpha2 ,2) - 4 * m_alpha1 * m_alpha3;
}

template <class T>
void gsArcLengthIterator<T>::computeLambdasEta()
{
  m_alpha1 = m_a0;
  m_alpha2 = m_b0 + m_eta*m_b1;
  m_alpha3 = m_c0 + m_eta*m_c1 + m_eta*m_eta*m_c2;

  // Lam 1992
  m_discriminant = math::pow(m_alpha2 ,2) - 4 * m_alpha1 * m_alpha3;
  m_deltaLs[0] = (-m_alpha2 + math::sqrt(m_discriminant))/(2*m_alpha1);
  m_deltaLs[1] = (-m_alpha2 - math::sqrt(m_discriminant))/(2*m_alpha1);

  // Zhou 1995
  // m_deltaLs[0] = (-m_alpha2 )/(2*m_alpha1);
  // m_deltaLs[1] = (-m_alpha2 )/(2*m_alpha1);
}

template <class T>
void gsArcLengthIterator<T>::computeLambdasModified()
{
  m_alpha1 = m_b1*m_b1 - 4*m_a0*m_c2;
  m_alpha2 = 2*m_b0*m_b1 - 4*m_a0*m_c1;
  m_alpha3 = m_b0*m_b0 - 4*m_a0*m_c0;

  m_discriminant = math::pow(m_alpha2 ,2) - 4 * m_alpha1 * m_alpha3;

  gsVector<T> etas(2);
  etas.setZero();
  if (m_discriminant >= 0)
    {
      etas[0] = (-m_alpha2 + math::pow(m_discriminant,0.5))/(2*m_alpha1);
      etas[1] = (-m_alpha2 - math::pow(m_discriminant,0.5))/(2*m_alpha1);

      T eta1 = std::min(etas[0],etas[1]);
      T eta2 = std::max(etas[0],etas[1]);
      gsInfo<<"eta 1 = "<<eta1<<"\t eta2 = "<<eta2<<"\n";

      // Approach of Zhou 1995
      // m_eta = std::min(1.0,eta2);
      // if (m_eta <= 0)
      //   gsInfo<<"Warning: both etas are non-positive!\n";
      // if (m_eta <=0.5)
      // {
      //   gsInfo<<"Warning: eta is small; bisecting step length\n";
      //   m_arcLength=m_arcLength/2.;
      // }

      // Approach of Lam 1992
      T xi = 0.05*abs(eta2-eta1);
      if (eta2<1.0)
        m_eta = eta2-xi;
      else if ( (eta2 > 1.0) && (-m_alpha2/m_alpha1 < 1.0) )
        m_eta = eta2+xi;
      else if ( (eta1 < 1.0) && (-m_alpha2/m_alpha1 > 1.0) )
        m_eta = eta1-xi;
      else if ( eta1 > 1.0 )
        m_eta = eta1 + xi;

      if (eta2<1.0)
        note = note + " option 1";
      else if ( (eta2 > 1.0) && (-m_alpha2/m_alpha1 < 1.0) )
        note = note + " option 2";
      else if ( (eta1 < 1.0) && (-m_alpha2/m_alpha1 > 1.0) )
        note = note + " option 3";
      else if ( eta1 > 1.0 )
        note = note + " option 4";

      // m_eta = eta2;
    }
  else
  {
    gsInfo<<"Discriminant was negative in modified method\n";
  }
}
template <class T>
void gsArcLengthIterator<T>::computeLambdasComplex()
{
  T A0 = math::pow(m_phi,2)* m_forcing.dot(m_forcing); // see Lam et al. 1991
  gsVector<T> DeltaUcr = m_DeltaU + m_deltaUbar;

  // Compute internal loads from residual and
  // gsVector<T> R = m_residualFun(m_U + DeltaUcr, m_L + m_DeltaL , m_forcing);
  // gsVector<T> Fint = R + (m_L + m_DeltaL) * m_forcing;
  gsVector<T> Fint = m_jacMat*(m_U+m_DeltaU);
  T Lcr = Fint.dot(m_forcing)/m_forcing.dot(m_forcing);
  T DeltaLcr = Lcr - m_L;

  T arcLength_cr = math::pow( DeltaUcr.dot(DeltaUcr) + A0 * math::pow(DeltaLcr, 2.0) ,0.5);
  T mu = m_arcLength/arcLength_cr;

  m_deltaL = mu*DeltaLcr - m_DeltaL;
  m_deltaU = mu*DeltaUcr - m_DeltaU;
}

template <class T>
void gsArcLengthIterator<T>::computeLambdas()
{
  m_deltaLs.setZero(2);
  computeLambdasSimple();
  if (m_discriminant >= 0)
  {
    m_eta = 1.0;
    m_deltaLs[0] = (-m_alpha2 + math::pow(m_discriminant,0.5))/(2*m_alpha1);
    m_deltaLs[1] = (-m_alpha2 - math::pow(m_discriminant,0.5))/(2*m_alpha1);
    computeLambdaDOT();
  }
  else
  {
    note += "\tC";
    // Compute eta
    computeLambdasModified();
    if ((m_discriminant >= 0) && m_eta > 0.002)
    {
      // recompute lambdas with new eta
      computeLambdasEta();
      // gsInfo<<"2: dL1 = "<<m_deltaLs[0]<<"\tdL2 = "<<m_deltaLs[1]<<"\t eta = "<<m_eta<<"\n";
      computeLambdaDOT();
      // gsInfo<<"2: dL1 = "<<m_deltaL<<"\t m_deltaU.norm = "<<m_deltaU.norm()<<"\t eta = "<<m_eta<<"\n";
      gsInfo<<"Modified Complex Root Solve\n";
    }
    else
    {
      // if the roots of the modified method are still complex, we use the following function (see Lam 1992, eq 13-17)
      m_eta = 1.0;
      computeLambdasComplex();
      // gsInfo<<"3: dL1 = "<<m_deltaL<<"\t m_deltaU.norm = "<<m_deltaU.norm()<<"\t eta = "<<m_eta<<"\n";
      gsInfo<<"Simplified Complex Root Solve\n";
      // Note: no selection of roots is needed
    }
  }
}

template <class T>
void gsArcLengthIterator<T>::computeLambdaDET()
{
    computeLambdas();

    if (sign(m_DeltaL + m_deltaLs[0]) == sign(m_detKT))
      m_deltaL = m_deltaLs[0];
    else
      m_deltaL = m_deltaLs[1];

    // Compute update of U (NOTE: m_eta=1.0)
    m_deltaU = m_deltaUbar + m_deltaL*m_deltaUt;

    // gsInfo<<"\t\t Choice based on DETERMINANT. Options:\n";
    // gsInfo<<"\t\t DeltaL = "<<m_DeltaL+m_deltaLs[0]<<" DeltaU.norm = "<<(m_DeltaU + m_deltaUbar + m_deltaLs[0]*m_deltaUt).norm()<<"\n";
    // gsInfo<<"\t\t DeltaL = "<<m_DeltaL+m_deltaLs[1]<<" DeltaU.norm = "<<(m_DeltaU + m_deltaUbar + m_deltaLs[1]*m_deltaUt).norm()<<"\n";
}

template <class T>
void gsArcLengthIterator<T>::computeLambdaMU()
{
    T A0 = math::pow(m_phi,2)* m_forcing.dot(m_forcing); // see Lam et al. 1991
    int dir = sign(m_DeltaUold.dot(m_deltaUt) + A0*m_DeltaLold); // Feng et al. 1995 with H = \Psi^2
    T denum = ( math::pow( m_deltaUt.dot(m_deltaUt) + A0 ,0.5) ); // Feng et al. 1995 with H = \Psi^2

    T mu;
    if (denum==0)
      mu = m_arcLength;
    else
      mu = m_arcLength / denum;

    m_deltaL = dir*mu;
    m_deltaU = m_deltaL*m_deltaUt;
}

template <class T>
void gsArcLengthIterator<T>::computeLambdaDOT()
{
    gsVector<T> deltaU1, deltaU2;
    deltaU1 = m_eta*m_deltaUbar + m_deltaUt*m_deltaLs[0];
    deltaU2 = m_eta*m_deltaUbar + m_deltaUt*m_deltaLs[1];

    // ---------------------------------------------------------------------------------
    // Method by Ritto-Corea et al. 2008
    T DOT1,DOT2;
    DOT1 = m_deltaLs[0]*(m_DeltaUold.dot(m_deltaUt) + math::pow(m_phi,2)*m_DeltaLold);
    DOT2 = m_deltaLs[1]*(m_DeltaUold.dot(m_deltaUt) + math::pow(m_phi,2)*m_DeltaLold);

    if (DOT1 > DOT2)
    {
      m_deltaL = m_deltaLs[0];
      m_deltaU = deltaU1;
    }
    else if (DOT1 < DOT2)
    {
      m_deltaL = m_deltaLs[1];
      m_deltaU = deltaU2;
    }
    else
    {
      m_deltaL = m_deltaLs[0];
      m_deltaU = deltaU1;
    }

    // ---------------------------------------------------------------------------------
    // // Method by Crisfield 1981
    // T DOT1,DOT2;
    // DOT1 = (m_DeltaUold+deltaU1).dot(m_DeltaUold);
    // DOT2 = (m_DeltaUold+deltaU2).dot(m_DeltaUold);

    // DOT1 = (m_DeltaU+deltaU1).dot(m_DeltaUold);
    // DOT2 = (m_DeltaU+deltaU2).dot(m_DeltaUold);

    // if ((DOT1 > DOT2) && DOT2 <= 0)
    // {
    //   m_deltaL = m_deltaLs[0];
    //   m_deltaU = deltaU1;
    // }
    // else if ((DOT1 < DOT2) && DOT1 <= 0)
    // {
    //   m_deltaL = m_deltaLs[1];
    //   m_deltaU = deltaU2;
    // }
    // else if ((DOT1 >=0) && (DOT2 >=0))
    // {
    //   T linsol = -m_alpha3/m_alpha2;
    //   T diff1 = abs(m_deltaLs[0]-linsol);
    //   T diff2 = abs(m_deltaLs[1]-linsol);
    //   note += "\t linear solution!";
    //   // note += "Linsol\t" + std::to_string(diff1) + "\t" + std::to_string(diff2) + "\t" + std::to_string(linsol) + "\t" + std::to_string(m_deltaLs[0]) + "\t" + std::to_string(m_deltaLs[1]) + "\n";
    //   if (diff1 > diff2)
    //   {
    //     m_deltaL = m_deltaLs[1];
    //     m_deltaU = deltaU2;
    //   }
    //   else
    //   {
    //     m_deltaL = m_deltaLs[0];
    //     m_deltaU = deltaU1;
    //   }
    // }
    // else
    // {
    //   m_deltaL = m_deltaLs[0];
    //   m_deltaU = deltaU1;
    // }

}

// ------------------------------------------------------------------------------------------------------------
// ---------------------------------------Singular point methods-----------------------------------------------
// ------------------------------------------------------------------------------------------------------------
template <class T>
void gsArcLengthIterator<T>::computeSingularPoint()
{
  computeSingularPoint(m_tolerance);
}

template <class T>
void gsArcLengthIterator<T>::computeSingularPoint(T tol)
{
  // First, approximate the eigenvector of the Jacobian by a few arc length iterations
  // Initiate m_V and m_DeltaVDET
  gsInfo<<"Starting with U.norm = "<<m_U.norm()<<" and L = "<<m_L<<"\n";

  m_V = gsVector<T>::Ones(m_numDof);
  m_V.normalize();
  m_DeltaV = gsVector<T>::Zero(m_numDof);
  m_DeltaU.setZero();
  m_DeltaL = 0.0;
  m_converged = false;
  m_solver.compute(m_jacMat);
  for (index_t k = 0; k<5; k++)
  {
    m_V = m_solver.solve(m_V);
    m_V.normalize();
  }

  m_basisResidualKTPhi = (m_jacobian(m_U).toDense()*m_V).norm();

  real_t tol2 = 1e-2;
  gsInfo<<"Extended Iterations:\t";
  // If a singular point is detected, the product of phi with the RHS (m_forcing) should be approximately zero. If so, start procedure to find singular point
  T test = (abs(m_V.dot(m_forcing)));
  if (test < tol2)
  {
    gsInfo<<"Singular point detected"<<"Singular point test returns: "<<test<<"\n";
    if (m_verbose)
      initOutputExtended();
    for (m_numIterations = 1; m_numIterations < m_maxIterations; ++m_numIterations)
    {
      extendedSystemIteration();
      m_U = m_U + m_DeltaU;
      m_L = m_L + m_DeltaL;
      m_V = m_V + m_DeltaV;
      // m_V.normalize();

      // m_resVec = m_residualFun(m_U, m_L, m_forcing);
      // m_residue = m_resVec.norm() / ( m_L * m_forcing.norm() );
      // m_residue = (m_jacobian(m_U).toDense()*m_V).norm() / refError;
      m_residueKTPhi = (m_jacobian(m_U).toDense()*m_V).norm(); // /m_basisResidualKTPhi;
      m_resVec = m_residualFun(m_U,m_L,m_forcing);
      m_residueF = m_resVec.norm();
      m_residueU = m_DeltaU.norm();
      m_residueL = m_DeltaL;
      if (m_verbose)
        stepOutputExtended();

      // termination criteria
      // if ( m_residueF < m_toleranceF && m_residueU < m_toleranceU &&  m_residueKTPhi< tol )
      if ( m_residueKTPhi< tol )
      {
          m_converged = true;
          gsInfo<<"Iterations finished. U.norm() = "<<m_U.norm()<<"\t L = "<<m_L<<"\n";
          break;
      }
      if ( (m_numIterations == m_maxIterations-1) && (!m_converged) )
      {
        gsInfo<<"Warning: Extended iterations did not converge! \n";
        gsInfo<<"Iterations finished. U.norm() = "<<m_U.norm()<<"\t L = "<<m_L<<"\n";
      }
    }
  }
  else
  {
    gsInfo<<"Limit point detected"<<"Singular point test returns: "<<test<<"\n";
  }

}

template <class T>
void gsArcLengthIterator<T>::computeSingularPoint(const gsVector<T> U, T L)
{
  computeSingularPoint(U, L, m_tolerance);
}

template <class T>
void gsArcLengthIterator<T>::computeSingularPoint(const gsVector<T> U, T L, T tol)
{
  m_U = U;
  m_L = L;
  computeSingularPoint(tol);
}

template <class T>
void gsArcLengthIterator<T>::extendedSystemIteration()
{
  m_resVec = m_residualFun(m_U, m_L, m_forcing);
  m_jacMat = m_jacobian(m_U);

  m_deltaUt = m_solver.solve(m_forcing); // DeltaV1
  m_deltaUbar = m_solver.solve(-m_resVec); // DeltaV2

  real_t eps = 1e-12;
  gsSparseMatrix<T> jacMatEps = m_jacobian(m_U + eps*m_V);
  gsVector<T> h1 = 1/eps * ( jacMatEps * m_deltaUt ) - 1/eps * m_forcing;
  gsVector<T> h2 = m_jacMat * m_V + 1/eps * ( jacMatEps * m_deltaUbar + m_resVec );

  m_solver.compute(m_jacMat);
  m_DeltaVt = m_solver.solve(-h1); // DeltaV1
  m_DeltaVbar = m_solver.solve(-h2); // DeltaV2

  m_DeltaL = -( ( m_V/m_V.norm() ).dot(m_DeltaVbar)  + m_V.norm() - 1) / ( m_V/m_V.norm() ).dot( m_DeltaVt );

  m_DeltaU = m_DeltaL * m_deltaUt + m_deltaUbar;
  m_DeltaV = m_DeltaL * m_DeltaVt + m_DeltaVbar;
}

template <class T>
int gsArcLengthIterator<T>::bisectionObjectiveFunction(const gsVector<T> x)
{
  gsVector<T> D;
  if (m_bifurcationMethod == "Determinant")
  {
    gsSparseMatrix<T> jac = m_jacobian(x);
    m_solver.compute(jac);
    D = m_solver.vectorD();
  }
  else // if (m_bifurcationMethod == "Determinant")
  {
    m_jacMat = m_jacobian(x);
    Eigen::SelfAdjointEigenSolver< gsMatrix<T> > es(m_jacMat);
    // gsInfo<<es.eigenvalues();
    D = es.eigenvalues();
  }
  return countNegatives(D);
}

template <class T>
T gsArcLengthIterator<T>::bisectionTerminationFunction(const gsVector<T> x)
{
  gsVector<T> D;
  gsVector<T> min;

  if (m_bifurcationMethod == "Determinant")
  {
    gsSparseMatrix<T> jac = m_jacobian(x);
    m_solver.compute(jac);
    D = m_solver.vectorD();
    min = D.colwise().minCoeff(); // This is required since D does not necessarily have one column.
    return min[0];
  }
  else // if (m_bifurcationMethod == "Determinant")
  {
    m_jacMat = m_jacobian(x);
    Eigen::SelfAdjointEigenSolver< gsMatrix<T> > es(m_jacMat);
    D = es.eigenvalues();
    min = D.colwise().minCoeff(); // This is required since D does not necessarily have one column.
    return min[0];
  }
}

template <class T>
void gsArcLengthIterator<T>::bisectionSolve(T tol)
{
  gsVector<T> U_old  = m_U;
  T L_old = m_L;
  // Store original arc length
  real_t dL = m_arcLength;
  bool adaptiveBool = m_adaptiveLength;
  m_adaptiveLength = false;

  T referenceError = bisectionTerminationFunction(m_U);

  // Store termination function value of initial solution U
  // T termU =  bisectionTerminationFunction(m_U);
  for (int k = 1; k < m_maxIterations; ++k)
  {
    // Reset start point
    m_U = U_old;
    m_L = L_old;
    m_arcLength = m_arcLength/2.0;

    // Make an arc length step; m_U and U_old and m_L and L_old are different
    step();

    // Objective functions on previous point
    int fa = bisectionObjectiveFunction(U_old);
    // Objective functions on new point
    int fb = bisectionObjectiveFunction(m_U);

    if (fa==fb)
    {
      U_old = m_U;
      L_old = m_L;
    }

    // termination criteria
    // relative error
    T term = bisectionTerminationFunction(m_U);
    if (m_verbose)
    {
      gsInfo<<"\t bisection iteration "<<k<<"\t arc length = "<<m_arcLength<<"\t relative error = "<< abs(term/referenceError)<<"\t"<<"obj.value = "<<term<<"\n";
    }
    // if (  m_arcLength < tol && abs(term) < tol )
    if (  abs(term/referenceError) < tol )
    {
        m_converged = true;
        break;
    }

  }
  // Reset arc length
  m_arcLength = dL;
  m_adaptiveLength = adaptiveBool;

  // Compute eigenvector
  m_V = gsVector<T>::Ones(m_numDof);
  m_V.normalize();
  gsVector<T> Vold = gsVector<T>::Zero(m_numDof);
  for (index_t k = 0; k<5; k++)
  {
    m_solver.compute(m_jacMat);
    m_V = m_solver.solve(m_V);
    m_V.normalize();
    if (  (m_V-Vold).norm() < m_tolerance )
    {
        m_converged = true;
        break;
    }
    Vold = m_V;
  }
}

template <class T>
void gsArcLengthIterator<T>::bisectionSolve()
{
  // First step: compute objective functions
  int fa = bisectionObjectiveFunction(m_U);
    // arc length step
    step();
  int fb = bisectionObjectiveFunction(m_U);

  if (fa == fb)
  {
    gsInfo<<"Improper choice of initial point and step"<<"\n";
    std::terminate();
  }
  else
  {
    bisectionSolve(m_tolerance);
  }
}

template <class T>
void gsArcLengthIterator<T>::bisectionSolve(const gsVector<T> U, T L, T tol)
{
  // Initialise: load U and L in memory positions m_U and m_L. Also, define temporary objects
  m_U = U;
  m_L = L;
  bisectionSolve(tol);
}

template <class T>
void gsArcLengthIterator<T>::bisectionSolve(const gsVector<T> U, T L)
{
  bisectionSolve(U,L,m_tolerance);
}


template <class T>
void gsArcLengthIterator<T>::checkBifurcation()
{
  // EXPENSIVE!!
  if (m_bifurcationMethod=="Eigenvalue")
  {
    m_jacMat = m_jacobian(m_U);
    Eigen::SelfAdjointEigenSolver< gsMatrix<T> > es(m_jacMat);
    // gsInfo<<es.eigenvalues();
    gsVector<T> eigenvalues = es.eigenvalues();
    if (eigenvalues[0] > 0)
    {
      m_bifurcation = false;
    }
    else
    {
      m_bifurcation = true;
    }
    m_indicator = eigenvalues.colwise().minCoeff()[0];
    gsInfo<<"\tMinimum eigenvalue : "<<m_indicator<<"\n";
  }
  else if (m_bifurcationMethod=="Determinant")
  {
    m_jacMat = m_jacobian(m_U);
    m_solver.compute(m_jacMat);
    gsVector<T> D = m_solver.vectorD();
    int N = countNegatives(D);
    if (N!=0)
    {
      m_bifurcation = true;
    }
    else
    {
      m_bifurcation = false;
    }
    m_indicator = D.colwise().minCoeff()[0];
    gsInfo<<"\tMinimum diagonal value : "<<m_indicator<<"\n";
  }
}

template <class T>
void gsArcLengthIterator<T>::switchBranch()
{
  m_V.normalize();
  real_t lenPhi = m_V.norm();
  real_t xi = lenPhi/m_tau;
  // gsInfo<<xi<<"\n";
  m_DeltaU = xi*m_V;
  m_U = m_U + m_DeltaU;

  // m_DeltaLold = 0.0;
  // m_DeltaUold.setZero();
  m_DeltaUold = m_DeltaU;
}

// ------------------------------------------------------------------------------------------------------------
// ---------------------------------------Output functions-----------------------------------------------------
// ------------------------------------------------------------------------------------------------------------

template <class T>
void gsArcLengthIterator<T>::initOutput()
{
  if ( (m_method == "Riks") || (m_method == "ConsistentCrisfield") || (m_method == "ExplicitIterations") )
    initOutputRiks();
  else if (m_method == "Crisfield")
    initOutputCrisfield();
  else
    gsInfo<<"Init output: Method unknown \n";
}


template <class T>
void gsArcLengthIterator<T>::initOutputRiks()
{
  gsInfo<<"\t";
  gsInfo<<std::setw(12)<<std::left<<"Iteration";
  gsInfo<<std::setw(17)<<std::left<<"ResidualF";
  gsInfo<<std::setw(17)<<std::left<<"ResidualU";
  gsInfo<<std::setw(17)<<std::left<<"ResidualL";
  gsInfo<<std::setw(17)<<std::left<<"U.norm";
  gsInfo<<std::setw(17)<<std::left<<"L";
  gsInfo<<std::setw(17)<<std::left<<"DU.norm";
  gsInfo<<std::setw(17)<<std::left<<"DL";
  gsInfo<<std::setw(17)<<std::left<<"dU.norm";
  gsInfo<<std::setw(17)<<std::left<<"dL";
  gsInfo<<std::setw(17)<<std::left<<"ds";
  gsInfo<<std::setw(17)<<std::left<<"dsU";
  gsInfo<<std::setw(17)<<std::left<<"dsL";
  gsInfo<<std::setw(17)<<std::left<<"Dmin";
  gsInfo<<std::setw(17)<<std::left<<"note";
  gsInfo<<"\n";

  note = "";
}

template <class T>
void gsArcLengthIterator<T>::initOutputCrisfield()
{
  gsInfo<<"\t";
  gsInfo<<std::setw(12)<<std::left<<"Iteration";
  gsInfo<<std::setw(17)<<std::left<<"ResidualF";
  gsInfo<<std::setw(17)<<std::left<<"ResidualU";
  gsInfo<<std::setw(17)<<std::left<<"ResidualL";
  gsInfo<<std::setw(17)<<std::left<<"U.norm";
  gsInfo<<std::setw(17)<<std::left<<"L";
  gsInfo<<std::setw(17)<<std::left<<"DU.norm";
  gsInfo<<std::setw(17)<<std::left<<"DL";
  gsInfo<<std::setw(17)<<std::left<<"dU.norm";
  gsInfo<<std::setw(17)<<std::left<<"dL";
  gsInfo<<std::setw(17)<<std::left<<"ds";
  gsInfo<<std::setw(17)<<std::left<<"dsU";
  gsInfo<<std::setw(17)<<std::left<<"dsL";
  gsInfo<<std::setw(17)<<std::left<<"Dmin";
  gsInfo<<std::setw(17)<<std::left<<"note";
  gsInfo<<"\n";

  note = "";
}

template <class T>
void gsArcLengthIterator<T>::initOutputExtended()
{
  gsInfo<<"\t";
  gsInfo<<std::setw(12)<<std::left<<"Iteration";
  gsInfo<<std::setw(17)<<std::left<<"ResidualF";
  gsInfo<<std::setw(17)<<std::left<<"ResidualU";
  gsInfo<<std::setw(17)<<std::left<<"ResidualL";
  gsInfo<<std::setw(17)<<std::left<<"K_T * phi";
  gsInfo<<std::setw(17)<<std::left<<"U.norm";
  gsInfo<<std::setw(17)<<std::left<<"phi.norm";
  gsInfo<<std::setw(17)<<std::left<<"L";
  gsInfo<<std::setw(17)<<std::left<<"DU.norm";
  gsInfo<<std::setw(17)<<std::left<<"Dphi.norm";
  gsInfo<<std::setw(17)<<std::left<<"DL";
  gsInfo<<std::setw(17)<<std::left<<"Dmin";
  gsInfo<<std::setw(17)<<std::left<<"note";
  gsInfo<<"\n";

  note = "";
}

template <class T>
void gsArcLengthIterator<T>::stepOutput()
{
  if ( (m_method == "Riks") || (m_method == "ConsistentCrisfield") || (m_method == "ExplicitIterations") )
    stepOutputRiks();
  else if (m_method == "Crisfield")
    stepOutputCrisfield();
  else
    gsInfo<<"Step output: Method unknown \n";
}

template <class T>
void gsArcLengthIterator<T>::stepOutputRiks()
{
  if (!m_quasiNewton)
  {
    gsVector<T> D = m_solver.vectorD();
    m_indicator = D.colwise().minCoeff()[0];
  }
  else
    m_indicator = 0;

  gsInfo<<"\t";
  gsInfo<<std::setw(12)<<std::left<<m_numIterations;
  gsInfo<<std::setw(17)<<std::left<<m_residueF;
  gsInfo<<std::setw(17)<<std::left<<m_residueU;
  gsInfo<<std::setw(17)<<std::left<<m_residueL;
  gsInfo<<std::setw(17)<<std::left<<(m_U+m_DeltaU).norm();
  gsInfo<<std::setw(17)<<std::left<<(m_L + m_DeltaL);
  gsInfo<<std::setw(17)<<std::left<<m_DeltaU.norm();
  gsInfo<<std::setw(17)<<std::left<<m_DeltaL;
  gsInfo<<std::setw(17)<<std::left<<m_deltaU.norm();
  gsInfo<<std::setw(17)<<std::left<<m_deltaL;
  gsInfo<<std::setw(17)<<std::left<<m_phi * math::pow(m_DeltaU.norm(),2.0) + (1.0-m_phi) * math::pow(m_DeltaL,2.0);
  gsInfo<<std::setw(17)<<std::left<<m_phi * math::pow(m_DeltaU.norm(),2.0);
  gsInfo<<std::setw(17)<<std::left<<(1-m_phi) * math::pow(m_DeltaL,2.0);
  gsInfo<<std::setw(17)<<std::left<<m_indicator;
  gsInfo<<std::setw(17)<<std::left<<note;
  gsInfo<<"\n";

  note = "";
}

template <class T>
void gsArcLengthIterator<T>::stepOutputCrisfield()
{
  if (!m_quasiNewton)
  {
    gsVector<T> D = m_solver.vectorD();
    m_indicator = D.colwise().minCoeff()[0];
  }
  else
    m_indicator = 0;

  T A0 = math::pow(m_phi,2)*m_forcing.dot(m_forcing);

  gsInfo<<"\t";
  gsInfo<<std::setw(12)<<std::left<<m_numIterations;
  gsInfo<<std::setw(17)<<std::left<<m_residueF;
  gsInfo<<std::setw(17)<<std::left<<m_residueU;
  gsInfo<<std::setw(17)<<std::left<<m_residueL;
  gsInfo<<std::setw(17)<<std::left<<(m_U+m_DeltaU).norm();
  gsInfo<<std::setw(17)<<std::left<<(m_L + m_DeltaL);
  gsInfo<<std::setw(17)<<std::left<<m_DeltaU.norm();
  gsInfo<<std::setw(17)<<std::left<<m_DeltaL;
  gsInfo<<std::setw(17)<<std::left<<m_deltaU.norm();
  gsInfo<<std::setw(17)<<std::left<<m_deltaL;
  gsInfo<<std::setw(17)<<std::left<<m_arcLength; //math::pow(m_DeltaU.dot(m_DeltaU) + A0*math::pow(m_DeltaL,2.0),0.5);
  gsInfo<<std::setw(17)<<std::left<<math::pow(m_DeltaU.norm(),2.0);
  gsInfo<<std::setw(17)<<std::left<<A0*math::pow(m_DeltaL,2.0);
  gsInfo<<std::setw(17)<<std::left<<m_indicator;
  gsInfo<<std::setw(17)<<std::left<<note;
  gsInfo<<"\n";

  note = "";
}

template <class T>
void gsArcLengthIterator<T>::stepOutputExtended()
{
  gsInfo<<"\t";
  gsInfo<<std::setw(12)<<std::left<<m_numIterations;
  gsInfo<<std::setw(17)<<std::left<<m_residueF;
  gsInfo<<std::setw(17)<<std::left<<m_residueU;
  gsInfo<<std::setw(17)<<std::left<<m_residueL;
  gsInfo<<std::setw(17)<<std::left<<m_residueKTPhi;
  gsInfo<<std::setw(17)<<std::left<<(m_U).norm();
  gsInfo<<std::setw(17)<<std::left<<(m_V).norm();
  gsInfo<<std::setw(17)<<std::left<<(m_L);
  gsInfo<<std::setw(17)<<std::left<<m_DeltaU.norm();
  gsInfo<<std::setw(17)<<std::left<<m_DeltaV.norm();
  gsInfo<<std::setw(17)<<std::left<<m_DeltaL;
  gsInfo<<std::setw(17)<<std::left<<bisectionTerminationFunction(m_U);
  gsInfo<<std::setw(17)<<std::left<<note;
  gsInfo<<"\n";

  note = "";
}

template <class T>
void gsArcLengthIterator<T>::printSettings()
{
  gsInfo<<"-------------------------------------------------------------------------------------\n";
  gsInfo<<"--------------------------------gsArcLengthIterator class----------------------------\n";
  gsInfo<<"-------------------------------------------------------------------------------------\n";
  gsInfo<<"Arc-Length Method:\t\t"<<m_method<<"\n";
  gsInfo<<"Arc-Length Distance:\t\t"<<m_arcLength<<"\n";
  gsInfo<<"Adaptive Length:\t\t"<<m_adaptiveLength<<"\n";
  gsInfo<<"Max. no. iterations:\t\t"<<m_maxIterations<<"\n";
  gsInfo<<"Quasi Newton Method:\t\t"<<m_quasiNewton<<"\n";
  gsInfo<<"Quasi Newton Internal:\t\t"<<m_quasiNewtonInterval<<"\n";
  gsInfo<<"Singular Point Identification:\t"<<m_bifurcationMethod<<"\n";
  gsInfo<<"Branching Perturbation:\t\t"<<m_tau<<"\n";
  gsInfo<<"Relaxation Factor:\t\t"<<m_relax<<"\n";
  gsInfo<<"-------------------------------------------------------------------------------------\n";
}

} // namespace gismo
