/** @file gsALMCrisfield.hpp

    @brief Performs the arc length method to solve a nonlinear equation system.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-..., TU Delft)
*/

#pragma once

#include <typeinfo>
#include <gsStructuralAnalysis/src/gsALMSolvers/gsALMHelper.h>

namespace gismo
{

template <class T>
void gsALMCrisfield<T>::defaultOptions()
{
    Base::defaultOptions();
    m_options.addReal("Scaling","Set Scaling factor Phi",-1);
    m_options.addInt ("AngleMethod","Angle determination method: 0 = Previous step; 1 = Previous iteration",angmethod::Step);
}

template <class T>
void gsALMCrisfield<T>::getOptions()
{
    Base::getOptions();
    m_phi                 = m_options.getReal("Scaling");
    m_phi_user = m_phi == -1 ? false : true;

    m_angleDetermine      = m_options.getInt ("AngleMethod");
}

template <class T>
void gsALMCrisfield<T>::initMethods()
{
  m_numDof = m_forcing.size();
  m_DeltaU = m_U = gsVector<T>::Zero(m_numDof);
  m_DeltaL = m_L = 0.0;

  m_DeltaUold = gsVector<T>::Zero(m_numDof);
  m_DeltaLold = 0.0;
}

// ------------------------------------------------------------------------------------------------------------
// ---------------------------------------Crisfield's method---------------------------------------------------
// ------------------------------------------------------------------------------------------------------------

template <class T>
void gsALMCrisfield<T>::quasiNewtonPredictor()
{
  m_jacMat = computeJacobian();
  this->factorizeMatrix(m_jacMat);
  computeUt(); // rhs does not depend on solution
  computeUbar(); // rhs contains residual and should be computed every time

}

template <class T>
void gsALMCrisfield<T>::quasiNewtonIteration()
{
  m_jacMat = computeJacobian();
  this->factorizeMatrix(m_jacMat);
  computeUt(); // rhs does not depend on solution
}

template <class T>
void gsALMCrisfield<T>::iteration()
{
  computeUbar(); // rhs contains residual and should be computed every time

  // Compute next solution
  m_eta = 1.0;

  T lamold = m_deltaL;
  computeLambdas();

  // Relaxation against oscillating load factor
  if (( (lamold*m_deltaL < 0) && (abs(m_deltaL) <= abs(lamold) ) ) && m_relax != 1.0 )
  {
     m_note += "\t relaxated solution!";
     m_deltaU = m_relax * (m_deltaL*m_deltaUt + m_eta*m_deltaUbar);
     m_deltaL = m_relax * m_deltaL;
  }

  m_DeltaU += m_deltaU;
  m_DeltaL += m_deltaL;

  if (m_angleDetermine == angmethod::Iteration)
  {
    m_DeltaUold = m_DeltaU;
    m_DeltaLold = m_DeltaL;
  }
}

template <class T>
void gsALMCrisfield<T>::initiateStep()
{
  m_DeltaU = m_deltaUbar = m_deltaUt = gsVector<T>::Zero(m_numDof);
  m_DeltaL = m_deltaL = 0.0;
  m_eta = 1.0;
}

template <class T>
void gsALMCrisfield<T>::predictor()
{
  m_jacMat = computeJacobian();
  this->factorizeMatrix(m_jacMat);
  m_deltaUt = this->solveSystem(m_forcing);

  // Choose Solution
  if (m_DeltaUold.dot(m_DeltaUold) == 0 && m_DeltaLold*m_DeltaLold == 0) // no information about previous step.
  {
    // gsWarn<<"Different predictor!!!!!\n";
    m_note+= "predictor\t";
    m_deltaL = m_arcLength / math::pow(2*( m_deltaUt.dot(m_deltaUt) ) , 0.5);
    m_deltaU = m_deltaUbar + m_deltaL*m_deltaUt;
    if (!m_phi_user)
      m_phi = math::pow( m_deltaUt.dot(m_deltaUt) / m_forcing.dot(m_forcing),0.5);
    m_note += " phi=" + std::to_string(m_phi);

  }
  // {
  //   m_note+= "predictor\t";

  //   if (!m_phi_user)
  //     m_phi = math::pow( m_deltaUt.dot(m_deltaUt) / m_forcing.dot(m_forcing),0.5);

  //   // m_deltaL = m_arcLength * DL / math::sqrt( 2*( m_deltaUt.dot(m_deltaUt) + m_DeltaL*DL ) );
  //   m_deltaL = m_arcLength / math::sqrt( ( m_deltaUt.dot(m_deltaUt) + m_phi*m_phi ) );

  //   // m_deltaU = m_arcLength * m_deltaUt / math::sqrt( m_deltaUt.dot(m_deltaUt) + m_DeltaL*DL );
  //   m_deltaU = m_deltaL*m_deltaUt;


  //   m_note += " phi=" + std::to_string(m_phi);

  //   // m_DeltaUold = m_deltaU;
  //   // m_DeltaLold = m_deltaL;
  // }
  else // previous point is not in the origin
  {
    if (!m_phi_user)
      m_phi = math::pow(m_U.dot(m_U)/( math::pow(m_L,2) * m_forcing.dot(m_forcing) ),0.5);
    m_note += " phi=" + std::to_string(m_phi);
    computeLambdaMU();
  }

  // Compute Temporary updates of DeltaL and DeltaU
  m_DeltaU += m_deltaU;
  m_DeltaL += m_deltaL;

  if (m_angleDetermine == angmethod::Iteration || m_angleDetermine == angmethod::Predictor)
  {
   m_DeltaUold = m_DeltaU;
   m_DeltaLold = m_DeltaL;
  }
}

template <class T>
void gsALMCrisfield<T>::predictorGuess()
{
  GISMO_ASSERT(m_Uguess.rows()!=0 && m_Uguess.cols()!=0,"Guess is empty");

  m_jacMat = computeJacobian();
  this->factorizeMatrix(m_jacMat);
  m_deltaUt = this->solveSystem(m_forcing);
  if (!m_phi_user)
    m_phi = math::pow( m_deltaUt.dot(m_deltaUt) / m_forcing.dot(m_forcing),0.5);
  m_note += " phi=" + std::to_string(m_phi);

  //
  m_DeltaUold = -(m_Uguess - m_U);
  m_DeltaLold = -(m_Lguess - m_L);

  // m_DeltaUold *= m_arcLength / math::sqrt( m_deltaU.dot(m_deltaU));
  // m_DeltaLold *= m_arcLength / math::sqrt( m_deltaU.dot(m_deltaU));

  computeLambdaMU();

  m_DeltaU = m_deltaU;
  m_DeltaL = m_deltaL;

  m_Uguess.resize(0);
}

template <class T>
void gsALMCrisfield<T>::iterationFinish()
{
  m_converged = true;
  m_Uprev = m_U;
  m_Lprev = m_L;
  m_U += m_DeltaU;
  m_L += m_DeltaL;
  if (m_angleDetermine == angmethod::Step)
  {
    m_DeltaUold = m_DeltaU;
    m_DeltaLold = m_DeltaL;
  }
}

// ------------------------------------------------------------------------------------------------------------
// ---------------------------------------Lambda computations--------------------------------------------------
// ------------------------------------------------------------------------------------------------------------

template <class T>
void gsALMCrisfield<T>::computeLambdasSimple() //Ritto-Corrêa et al. 2008
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
  // m_note += "\t D = " + std::to_string(m_discriminant);
}

template <class T>
void gsALMCrisfield<T>::computeLambdasEta()
{
  m_alpha1 = m_a0;
  m_alpha2 = m_b0 + m_eta*m_b1;
  m_alpha3 = m_c0 + m_eta*m_c1 + m_eta*m_eta*m_c2;

  // Lam 1992
  // m_discriminant = math::pow(m_alpha2 ,2) - 4 * m_alpha1 * m_alpha3;
  // m_deltaLs[0] = (-m_alpha2 + math::sqrt(m_discriminant))/(2*m_alpha1);
  // m_deltaLs[1] = (-m_alpha2 - math::sqrt(m_discriminant))/(2*m_alpha1);

  // Zhou 1995
  m_deltaLs[0] = (-m_alpha2 )/(2*m_alpha1);
  m_deltaLs[1] = (-m_alpha2 )/(2*m_alpha1);
}

template <class T>
void gsALMCrisfield<T>::computeLambdasModified()
{
  m_alpha1 = m_b1*m_b1 - 4.0*m_a0*m_c2;
  m_alpha2 = 2.0*m_b0*m_b1 - 4.0*m_a0*m_c1;
  m_alpha3 = m_b0*m_b0 - 4.0*m_a0*m_c0;

  m_discriminant = math::pow(m_alpha2 ,2.0) - 4.0 * m_alpha1 * m_alpha3;

  gsVector<T> etas(2);
  etas.setZero();
  if (m_discriminant >= 0)
  {
    etas[0] = (-m_alpha2 + math::pow(m_discriminant,0.5))/(2.0*m_alpha1);
    etas[1] = (-m_alpha2 - math::pow(m_discriminant,0.5))/(2.0*m_alpha1);

    T eta1 = std::min(etas[0],etas[1]);
    T eta2 = std::max(etas[0],etas[1]);
    if (m_verbose) {gsInfo<<"eta 1 = "<<eta1<<"\t eta2 = "<<eta2<<"\n";}

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
      m_note = m_note + " option 1";
    else if ( (eta2 > 1.0) && (-m_alpha2/m_alpha1 < 1.0) )
      m_note = m_note + " option 2";
    else if ( (eta1 < 1.0) && (-m_alpha2/m_alpha1 > 1.0) )
      m_note = m_note + " option 3";
    else if ( eta1 > 1.0 )
      m_note = m_note + " option 4";

    // m_eta = eta2;
  }
  else
  {
    gsInfo<<"Discriminant was negative in modified method\n";
  }
}
template <class T>
void gsALMCrisfield<T>::computeLambdasComplex()
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
void gsALMCrisfield<T>::computeLambdas()
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
    m_note += "\tC";
    // Compute eta
    computeLambdasModified();
    gsDebugVar(m_discriminant);
    if ((m_discriminant >= 0) && m_eta > 0.05)
    {
      // recompute lambdas with new eta
      computeLambdasEta();
      // gsInfo<<"2: dL1 = "<<m_deltaLs[0]<<"\tdL2 = "<<m_deltaLs[1]<<"\t eta = "<<m_eta<<"\n";
      computeLambdaDOT();
      // gsInfo<<"2: dL1 = "<<m_deltaL<<"\t m_deltaU.norm = "<<m_deltaU.norm()<<"\t eta = "<<m_eta<<"\n";
      if (m_verbose) {gsInfo<<"Modified Complex Root Solve\n";}
    }
    else
    {
      // if the roots of the modified method are still complex, we use the following function (see Lam 1992, eq 13-17)
      m_eta = 1.0;
      computeLambdasComplex();
      // gsInfo<<"3: dL1 = "<<m_deltaL<<"\t m_deltaU.norm = "<<m_deltaU.norm()<<"\t eta = "<<m_eta<<"\n";
      if (m_verbose) {gsInfo<<"Simplified Complex Root Solve\n";}
      // Note: no selection of roots is needed
    }
  }
}

template <class T>
void gsALMCrisfield<T>::computeLambdaDET()
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
void gsALMCrisfield<T>::computeLambdaMU()
{
    T A0 = math::pow(m_phi,2)* m_forcing.dot(m_forcing); // see Lam et al. 1991
    index_t dir = sign(m_DeltaUold.dot(m_deltaUt) + A0*m_DeltaLold); // Feng et al. 1995 with H = \Psi^2
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
void gsALMCrisfield<T>::computeLambdaDOT()
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
    //   m_note += "\t linear solution!";
    //   // m_note += "Linsol\t" + std::to_string(diff1) + "\t" + std::to_string(diff2) + "\t" + std::to_string(linsol) + "\t" + std::to_string(m_deltaLs[0]) + "\t" + std::to_string(m_deltaLs[1]) + "\n";
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
// ---------------------------------------Output functions-----------------------------------------------------
// ------------------------------------------------------------------------------------------------------------

template <class T>
void gsALMCrisfield<T>::initOutput()
{
  gsInfo<<"\t";
  gsInfo<<std::setw(4)<<std::left<<"It.";
  gsInfo<<std::setw(17)<<std::left<<"Res. F";
  gsInfo<<std::setw(17)<<std::left<<"|dU|/|Du|";
  gsInfo<<std::setw(17)<<std::left<<"dL/DL";
  gsInfo<<std::setw(17)<<std::left<<"|U|";
  gsInfo<<std::setw(17)<<std::left<<"L";
  gsInfo<<std::setw(17)<<std::left<<"|DU|";
  gsInfo<<std::setw(17)<<std::left<<"DL";
  gsInfo<<std::setw(17)<<std::left<<"|dU|";
  gsInfo<<std::setw(17)<<std::left<<"dL";
  gsInfo<<std::setw(17)<<std::left<<"ds²";
  gsInfo<<std::setw(17)<<std::left<<"|dU|²";
  gsInfo<<std::setw(17)<<std::left<<"dL²";
  gsInfo<<std::setw(17)<<std::left<<"Dmin";
  gsInfo<<std::setw(17)<<std::left<<"note";
  gsInfo<<"\n";

  m_note = "";
}

template <class T>
void gsALMCrisfield<T>::stepOutput()
{
  // if (!m_quasiNewton)
  // {
    computeStability(false);
  // }
  // else
  //   m_indicator = 0;

  T A0 = math::pow(m_phi,2)*m_forcing.dot(m_forcing);

  gsInfo<<"\t";
  gsInfo<<std::setw(4)<<std::left<<m_numIterations;
  gsInfo<<std::setw(17)<<std::left<<m_residueF;
  gsInfo<<std::setw(17)<<std::left<<m_residueU;
  gsInfo<<std::setw(17)<<std::left<<m_residueL;
  gsInfo<<std::setw(17)<<std::left<<(m_U+m_DeltaU).norm();
  gsInfo<<std::setw(17)<<std::left<<(m_L + m_DeltaL);
  gsInfo<<std::setw(17)<<std::left<<m_DeltaU.norm();
  gsInfo<<std::setw(17)<<std::left<<m_DeltaL;
  gsInfo<<std::setw(17)<<std::left<<m_deltaU.norm();
  gsInfo<<std::setw(17)<<std::left<<m_deltaL;
  gsInfo<<std::setw(17)<<std::left<<this->distance(m_DeltaU,m_DeltaL);//math::pow(m_DeltaU.dot(m_DeltaU) + A0*math::pow(m_DeltaL,2.0),0.5);
  gsInfo<<std::setw(17)<<std::left<<math::pow(m_DeltaU.norm(),2.0);
  gsInfo<<std::setw(17)<<std::left<<A0*math::pow(m_DeltaL,2.0);
  gsInfo<<std::setw(17)<<std::left<<m_indicator <<std::left << " (" <<std::left<< m_negatives<<std::left << ")";
  gsInfo<<std::setw(17)<<std::left<<m_note;
  gsInfo<<"\n";

  m_note = "";
}

} // namespace gismo