/*!
 * \file flow_sources.cpp
 * \brief Implementation of numerics classes for integration
 *        of source terms in fluid flow problems.
 * \author F. Palacios, T. Economon
 * \version 7.3.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../../../include/numerics/flow/flow_sources.hpp"
#include "../../../../Common/include/toolboxes/geometry_toolbox.hpp"

CSourceBase_Flow::CSourceBase_Flow(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config) :
                                   CNumerics(val_nDim, val_nVar, config) {
  residual = new su2double [nVar]();
  jacobian = new su2double* [nVar];
  for(unsigned short iVar = 0; iVar < nVar; ++iVar)
    jacobian[iVar] = new su2double [nVar]();
}

CSourceBase_Flow::~CSourceBase_Flow() {
  delete [] residual;
  if(jacobian) {
    for(unsigned short iVar = 0; iVar < nVar; ++iVar)
      delete [] jacobian[iVar];
    delete [] jacobian;
  }
}

CSourceAxisymmetric_Flow::CSourceAxisymmetric_Flow(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config) :
                          CSourceBase_Flow(val_nDim, val_nVar, config) {

  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;

  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  viscous = config->GetViscous();
  rans = (config->GetKind_Turb_Model() != TURB_MODEL::NONE);

}

CNumerics::ResidualType<> CSourceAxisymmetric_Flow::ComputeResidual(const CConfig* config) {

  su2double Pressure_i, Enthalpy_i, Velocity_i, sq_vel;
  unsigned short iDim, iVar, jVar;

  if (Coord_i[1] > EPS) {

    yinv = 1.0/Coord_i[1];

    sq_vel = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Velocity_i = U_i[iDim+1] / U_i[0];
      sq_vel += Velocity_i *Velocity_i;
    }

    Pressure_i = (Gamma-1.0)*U_i[0]*(U_i[nDim+1]/U_i[0]-0.5*sq_vel);
    Enthalpy_i = (U_i[nDim+1] + Pressure_i) / U_i[0];

    residual[0] = yinv*Volume*U_i[2];
    residual[1] = yinv*Volume*U_i[1]*U_i[2]/U_i[0];
    residual[2] = yinv*Volume*(U_i[2]*U_i[2]/U_i[0]);
    residual[3] = yinv*Volume*Enthalpy_i*U_i[2];

    /*--- Inviscid component of the source term. ---*/

    if (implicit) {
      jacobian[0][0] = 0.0;
      jacobian[0][1] = 0.0;
      jacobian[0][2] = 1.0;
      jacobian[0][3] = 0.0;

      jacobian[1][0] = -U_i[1]*U_i[2]/(U_i[0]*U_i[0]);
      jacobian[1][1] = U_i[2]/U_i[0];
      jacobian[1][2] = U_i[1]/U_i[0];
      jacobian[1][3] = 0.0;

      jacobian[2][0] = -U_i[2]*U_i[2]/(U_i[0]*U_i[0]);
      jacobian[2][1] = 0.0;
      jacobian[2][2] = 2*U_i[2]/U_i[0];
      jacobian[2][3] = 0.0;

      jacobian[3][0] = -Gamma*U_i[2]*U_i[3]/(U_i[0]*U_i[0]) + (Gamma-1)*U_i[2]*(U_i[1]*U_i[1]+U_i[2]*U_i[2])/(U_i[0]*U_i[0]*U_i[0]);
      jacobian[3][1] = -(Gamma-1)*U_i[2]*U_i[1]/(U_i[0]*U_i[0]);
      jacobian[3][2] = Gamma*U_i[3]/U_i[0] - 1/2*(Gamma-1)*( (U_i[1]*U_i[1]+U_i[2]*U_i[2])/(U_i[0]*U_i[0]) + 2*U_i[2]*U_i[2]/(U_i[0]*U_i[0]) );
      jacobian[3][3] = Gamma*U_i[2]/U_i[0];

      for (iVar=0; iVar < nVar; iVar++)
        for (jVar=0; jVar < nVar; jVar++)
          jacobian[iVar][jVar] *= yinv*Volume;

    }

    /*--- Add the viscous terms if necessary. ---*/

    if (viscous) ResidualDiffusion();

  }

  else {

    for (iVar=0; iVar < nVar; iVar++)
      residual[iVar] = 0.0;

    if (implicit) {
      for (iVar=0; iVar < nVar; iVar++) {
        for (jVar=0; jVar < nVar; jVar++)
          jacobian[iVar][jVar] = 0.0;
      }
    }

  }

  return ResidualType<>(residual, jacobian, nullptr);
}

void CSourceAxisymmetric_Flow::ResidualDiffusion(){

  if (!rans){ turb_ke_i = 0.0; }

  su2double laminar_viscosity_i    = V_i[nDim+5];
  su2double eddy_viscosity_i       = V_i[nDim+6];
  su2double thermal_conductivity_i = V_i[nDim+7];
  su2double heat_capacity_cp_i     = V_i[nDim+8];

  su2double total_viscosity_i = laminar_viscosity_i + eddy_viscosity_i;
  su2double total_conductivity_i = thermal_conductivity_i + heat_capacity_cp_i*eddy_viscosity_i/Prandtl_Turb;

  su2double u = U_i[1]/U_i[0];
  su2double v = U_i[2]/U_i[0];

  residual[0] -= 0.0;
  residual[1] -= Volume*(yinv*total_viscosity_i*(PrimVar_Grad_i[1][1]+PrimVar_Grad_i[2][0])
                         -TWO3*AuxVar_Grad_i[0][0]);
  residual[2] -= Volume*(yinv*total_viscosity_i*2*(PrimVar_Grad_i[2][1]-v*yinv)
                         -TWO3*AuxVar_Grad_i[0][1]);
  residual[3] -= Volume*(yinv*(total_viscosity_i*(u*(PrimVar_Grad_i[2][0]+PrimVar_Grad_i[1][1])
                                                 +v*TWO3*(2*PrimVar_Grad_i[2][1]-PrimVar_Grad_i[1][0]
                                                 -v*yinv+U_i[0]*turb_ke_i))
                                                 +total_conductivity_i*PrimVar_Grad_i[0][1])
                                                 -TWO3*(AuxVar_Grad_i[1][1]+AuxVar_Grad_i[2][0]));
}


CNumerics::ResidualType<> CSourceGeneralAxisymmetric_Flow::ComputeResidual(const CConfig* config) {
  unsigned short iVar, jVar;

  if (Coord_i[1] > EPS) {

    yinv = 1.0/Coord_i[1];

    su2double Density_i = U_i[0];
    su2double Velocity1_i = U_i[1]/U_i[0];
    su2double Velocity2_i = U_i[2]/U_i[0];
    su2double Energy_i = U_i[3]/U_i[0];

    su2double Pressure_i = V_j[3];
    su2double Enthalpy_i = Energy_i + Pressure_i/Density_i;

    /*--- Inviscid component of the source term. ---*/

    residual[0] = yinv*Volume*U_i[2];
    residual[1] = yinv*Volume*U_i[1]*Velocity2_i;
    residual[2] = yinv*Volume*U_i[2]*Velocity2_i;
    residual[3] = yinv*Volume*U_i[2]*Enthalpy_i;

    if (implicit) {

      su2double dPdrho_e_i = S_i[0];
      su2double dPde_rho_i = S_i[1];

      jacobian[0][0] = 0.0;
      jacobian[0][1] = 0.0;
      jacobian[0][2] = 1.0;
      jacobian[0][3] = 0.0;

      jacobian[1][0] = -Velocity1_i*Velocity2_i;
      jacobian[1][1] = Velocity2_i;
      jacobian[1][2] = Velocity1_i;
      jacobian[1][3] = 0.0;

      jacobian[2][0] = -Velocity2_i*Velocity2_i;
      jacobian[2][1] = 0.0;
      jacobian[2][2] = 2*Velocity2_i;
      jacobian[2][3] = 0.0;

      jacobian[3][0] = Velocity2_i*(dPdrho_e_i + dPde_rho_i/Density_i*(Velocity1_i*Velocity1_i
                                                                        + Velocity2_i*Velocity2_i
                                                                        - Energy_i) - Enthalpy_i);
      jacobian[3][1] = -Velocity1_i*Velocity2_i/Density_i *dPde_rho_i;
      jacobian[3][2] = Enthalpy_i - Velocity2_i*Velocity2_i/Density_i *dPde_rho_i;
      jacobian[3][3] = Velocity2_i + Velocity2_i/Density_i *dPde_rho_i;

      for (iVar=0; iVar < nVar; iVar++)
        for (jVar=0; jVar < nVar; jVar++)
          jacobian[iVar][jVar] *= yinv*Volume;

    }

    /*--- Add the viscous terms if necessary. ---*/

    if (viscous) ResidualDiffusion();

  }

  else {

    for (iVar=0; iVar < nVar; iVar++)
      residual[iVar] = 0.0;

    if (implicit) {
      for (iVar=0; iVar < nVar; iVar++) {
        for (jVar=0; jVar < nVar; jVar++)
          jacobian[iVar][jVar] = 0.0;
      }
    }

  }

  return ResidualType<>(residual, jacobian, nullptr);
}

CSourceBodyForce::CSourceBodyForce(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config) :
                  CSourceBase_Flow(val_nDim, val_nVar, config) {

  /*--- Store the pointer to the constant body force vector. ---*/

  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    Body_Force_Vector[iDim] = config->GetBody_Force_Vector()[iDim];

}

CNumerics::ResidualType<> CSourceBodyForce::ComputeResidual(const CConfig* config) {

  unsigned short iDim;
  su2double Force_Ref = config->GetForce_Ref();

  /*--- Zero the continuity contribution ---*/

  residual[0] = 0.0;

  /*--- Momentum contribution ---*/

  for (iDim = 0; iDim < nDim; iDim++)
    residual[iDim+1] = -Volume * U_i[0] * Body_Force_Vector[iDim] / Force_Ref;

  /*--- Energy contribution ---*/

  residual[nDim+1] = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    residual[nDim+1] += -Volume * U_i[iDim+1] * Body_Force_Vector[iDim] / Force_Ref;

  return ResidualType<>(residual, jacobian, nullptr);
}

CSourceGravity::CSourceGravity(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config) :
                CSourceBase_Flow(val_nDim, val_nVar, config) {
                  Force_Ref = config->GetForce_Ref();
                }

CNumerics::ResidualType<> CSourceGravity::ComputeResidual(const CConfig* config) {

  unsigned short iVar;

  for (iVar = 0; iVar < nVar; iVar++)
    residual[iVar] = 0.0;

  /*--- Evaluate the source term  ---*/
  residual[nDim] = Volume * U_i[0] * STANDARD_GRAVITY / Force_Ref;

  return ResidualType<>(residual, jacobian, nullptr);
}

CSourceRotatingFrame_Flow::CSourceRotatingFrame_Flow(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config) :
                           CSourceBase_Flow(val_nDim, val_nVar, config) {

  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;

}

CNumerics::ResidualType<> CSourceRotatingFrame_Flow::ComputeResidual(const CConfig* config) {

  unsigned short iDim, iVar, jVar;
  su2double Omega[MAXNDIM] = {0}, Momentum[MAXNDIM] = {0};

  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

  /*--- Retrieve the angular velocity vector from config. ---*/

  for (iDim = 0; iDim < 3; iDim++){
    Omega[iDim] = config->GetRotation_Rate(iDim)/config->GetOmega_Ref();
  }

  /*--- Get the momentum vector at the current node. ---*/

  for (iDim = 0; iDim < nDim; iDim++)
    Momentum[iDim] = U_i[iDim+1];

  /*--- Calculate rotating frame source term as ( Omega X Rho-U ) ---*/

  if (nDim == 2) {
    residual[0] = 0.0;
    residual[1] = (Omega[1]*Momentum[2] - Omega[2]*Momentum[1])*Volume;
    residual[2] = (Omega[2]*Momentum[0] - Omega[0]*Momentum[2])*Volume;
    residual[3] = 0.0;
  } else {
    residual[0] = 0.0;
    residual[1] = (Omega[1]*Momentum[2] - Omega[2]*Momentum[1])*Volume;
    residual[2] = (Omega[2]*Momentum[0] - Omega[0]*Momentum[2])*Volume;
    residual[3] = (Omega[0]*Momentum[1] - Omega[1]*Momentum[0])*Volume;
    residual[4] = 0.0;
  }

  /*--- Calculate the source term Jacobian ---*/

  if (implicit) {
    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++)
        jacobian[iVar][jVar] = 0.0;
    if (nDim == 2) {
      jacobian[1][2] = -Omega[2]*Volume;
      jacobian[2][1] =  Omega[2]*Volume;
    } else {
      jacobian[1][2] = -Omega[2]*Volume;
      jacobian[1][3] =  Omega[1]*Volume;
      jacobian[2][1] =  Omega[2]*Volume;
      jacobian[2][3] = -Omega[0]*Volume;
      jacobian[3][1] = -Omega[1]*Volume;
      jacobian[3][2] =  Omega[0]*Volume;
    }
  }

  return ResidualType<>(residual, jacobian, nullptr);
}

