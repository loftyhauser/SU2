/*!
 * \file CSolver.cpp
 * \brief Main subroutines for CSolver class.
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


#include "../../include/solvers/CSolver.hpp"
#include "../../include/gradients/computeGradientsGreenGauss.hpp"
#include "../../include/gradients/computeGradientsLeastSquares.hpp"
#include "../../include/limiters/computeLimiters.hpp"
#include "../../../Common/include/toolboxes/printing_toolbox.hpp"
#include "../../../Common/include/toolboxes/C1DInterpolation.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"
#include "../../../Common/include/toolboxes/CLinearPartitioner.hpp"
#include "../../../Common/include/adt/CADTPointsOnlyClass.hpp"
#include "../../include/CMarkerProfileReaderFVM.hpp"


CSolver::CSolver(LINEAR_SOLVER_MODE linear_solver_mode) : System(linear_solver_mode) {

  rank = SU2_MPI::GetRank();
  size = SU2_MPI::GetSize();

  adjoint = false;

  /*--- Set the multigrid level to the finest grid. This can be
        overwritten in the constructors of the derived classes. ---*/
  MGLevel = MESH_0;

  /*--- Array initialization ---*/

  OutputHeadingNames = nullptr;
  Residual           = nullptr;
  Residual_i         = nullptr;
  Residual_j         = nullptr;
  Solution           = nullptr;
  Solution_i         = nullptr;
  Solution_j         = nullptr;
  Vector             = nullptr;
  Vector_i           = nullptr;
  Vector_j           = nullptr;
  Res_Conv           = nullptr;
  Res_Visc           = nullptr;
  Res_Sour           = nullptr;
  Res_Conv_i         = nullptr;
  Res_Visc_i         = nullptr;
  Res_Conv_j         = nullptr;
  Res_Visc_j         = nullptr;
  Jacobian_i         = nullptr;
  Jacobian_j         = nullptr;
  Jacobian_ii        = nullptr;
  Jacobian_ij        = nullptr;
  Jacobian_ji        = nullptr;
  Jacobian_jj        = nullptr;
  Restart_Vars       = nullptr;
  Restart_Data       = nullptr;
  base_nodes         = nullptr;
  nOutputVariables   = 0;
  ResLinSolver       = 0.0;

  /*--- Variable initialization to avoid valgrid warnings when not used. ---*/

  IterLinSolver = 0;

  /*--- Flags for the periodic BC communications. ---*/

  rotate_periodic   = false;
  implicit_periodic = false;

  /*--- Containers to store the markers. ---*/
  nMarker = 0;

  /*--- Flags for the dynamic grid (rigid movement or unsteady deformation). ---*/
  dynamic_grid = false;

  /*--- Auxiliary data needed for CFL adaption. ---*/

  Old_Func = 0;
  New_Func = 0;
  NonLinRes_Counter = 0;

  nPrimVarGrad = 0;
  nPrimVar     = 0;

}

CSolver::~CSolver(void) {

  unsigned short iVar;

  /*--- Public variables, may be accessible outside ---*/

  delete [] OutputHeadingNames;

  /*--- Private ---*/

  delete [] Residual;
  delete [] Residual_i;
  delete [] Residual_j;
  delete [] Solution;
  delete [] Solution_i;
  delete [] Solution_j;
  delete [] Vector;
  delete [] Vector_i;
  delete [] Vector_j;
  delete [] Res_Conv;
  delete [] Res_Visc;
  delete [] Res_Sour;
  delete [] Res_Conv_i;
  delete [] Res_Visc_i;
  delete [] Res_Visc_j;

  if (Jacobian_i != nullptr) {
    for (iVar = 0; iVar < nVar; iVar++)
      delete [] Jacobian_i[iVar];
    delete [] Jacobian_i;
  }

  if (Jacobian_j != nullptr) {
    for (iVar = 0; iVar < nVar; iVar++)
      delete [] Jacobian_j[iVar];
    delete [] Jacobian_j;
  }

  if (Jacobian_ii != nullptr) {
    for (iVar = 0; iVar < nVar; iVar++)
      delete [] Jacobian_ii[iVar];
    delete [] Jacobian_ii;
  }

  if (Jacobian_ij != nullptr) {
    for (iVar = 0; iVar < nVar; iVar++)
      delete [] Jacobian_ij[iVar];
    delete [] Jacobian_ij;
  }

  if (Jacobian_ji != nullptr) {
    for (iVar = 0; iVar < nVar; iVar++)
      delete [] Jacobian_ji[iVar];
    delete [] Jacobian_ji;
  }

  if (Jacobian_jj != nullptr) {
    for (iVar = 0; iVar < nVar; iVar++)
      delete [] Jacobian_jj[iVar];
    delete [] Jacobian_jj;
  }

  delete [] Restart_Vars;
  delete [] Restart_Data;

}

void CSolver::GetPeriodicCommCountAndType(const CConfig* config,
                                          unsigned short commType,
                                          unsigned short &COUNT_PER_POINT,
                                          unsigned short &MPI_TYPE,
                                          unsigned short &ICOUNT,
                                          unsigned short &JCOUNT) const {
  switch (commType) {
    case PERIODIC_VOLUME:
      COUNT_PER_POINT  = 1;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      break;
    case PERIODIC_NEIGHBORS:
      COUNT_PER_POINT  = 1;
      MPI_TYPE         = COMM_TYPE_UNSIGNED_SHORT;
      break;
    case PERIODIC_RESIDUAL:
      COUNT_PER_POINT  = nVar + nVar*nVar + 1;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      break;
    case PERIODIC_IMPLICIT:
      COUNT_PER_POINT  = nVar;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      break;
    case PERIODIC_LAPLACIAN:
      COUNT_PER_POINT  = nVar;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      break;
    case PERIODIC_MAX_EIG:
      COUNT_PER_POINT  = 1;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      break;
    case PERIODIC_SENSOR:
      COUNT_PER_POINT  = 2;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      break;
    case PERIODIC_SOL_GG:
    case PERIODIC_SOL_GG_R:
      COUNT_PER_POINT  = nVar*nDim;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      ICOUNT           = nVar;
      JCOUNT           = nDim;
      break;
    case PERIODIC_PRIM_GG:
    case PERIODIC_PRIM_GG_R:
      COUNT_PER_POINT  = nPrimVarGrad*nDim;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      ICOUNT           = nPrimVarGrad;
      JCOUNT           = nDim;
      break;
    case PERIODIC_SOL_LS:
    case PERIODIC_SOL_ULS:
    case PERIODIC_SOL_LS_R:
    case PERIODIC_SOL_ULS_R:
      COUNT_PER_POINT  = nDim*nDim + nVar*nDim;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      ICOUNT           = nVar;
      JCOUNT           = nDim;
      break;
    case PERIODIC_PRIM_LS:
    case PERIODIC_PRIM_ULS:
    case PERIODIC_PRIM_LS_R:
    case PERIODIC_PRIM_ULS_R:
      COUNT_PER_POINT  = nDim*nDim + nPrimVarGrad*nDim;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      ICOUNT           = nPrimVarGrad;
      JCOUNT           = nDim;
      break;
    case PERIODIC_LIM_PRIM_1:
      COUNT_PER_POINT  = nPrimVarGrad*2;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      ICOUNT           = nPrimVarGrad;
      break;
    case PERIODIC_LIM_PRIM_2:
      COUNT_PER_POINT  = nPrimVarGrad;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      ICOUNT           = nPrimVarGrad;
      break;
    case PERIODIC_LIM_SOL_1:
      COUNT_PER_POINT  = nVar*2;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      ICOUNT           = nVar;
      break;
    case PERIODIC_LIM_SOL_2:
      COUNT_PER_POINT  = nVar;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      ICOUNT           = nVar;
      break;
    default:
      SU2_MPI::Error("Unrecognized quantity for periodic communication.",
                     CURRENT_FUNCTION);
      break;
  }
}

namespace PeriodicCommHelpers {
  CVectorOfMatrix& selectGradient(CVariable* nodes, unsigned short commType) {
    switch(commType) {
      case PERIODIC_PRIM_GG:
      case PERIODIC_PRIM_LS:
      case PERIODIC_PRIM_ULS:
        return nodes->GetGradient_Primitive();
        break;
      case PERIODIC_SOL_GG:
      case PERIODIC_SOL_LS:
      case PERIODIC_SOL_ULS:
        return nodes->GetGradient();
        break;
      default:
        return nodes->GetGradient_Reconstruction();
        break;
    }
  }

  const su2activematrix& selectField(CVariable* nodes, unsigned short commType) {
    switch(commType) {
      case PERIODIC_PRIM_GG:
      case PERIODIC_PRIM_LS:
      case PERIODIC_PRIM_ULS:
      case PERIODIC_PRIM_GG_R:
      case PERIODIC_PRIM_LS_R:
      case PERIODIC_PRIM_ULS_R:
      case PERIODIC_LIM_PRIM_1:
      case PERIODIC_LIM_PRIM_2:
        return nodes->GetPrimitive();
        break;
      default:
        return nodes->GetSolution();
        break;
    }
  }

  su2activematrix& selectLimiter(CVariable* nodes, unsigned short commType) {
    switch(commType) {
      case PERIODIC_LIM_PRIM_1:
      case PERIODIC_LIM_PRIM_2:
        return nodes->GetLimiter_Primitive();
        break;
      default:
        return nodes->GetLimiter();
        break;
    }
  }
}

void CSolver::InitiatePeriodicComms(CGeometry *geometry,
                                    const CConfig *config,
                                    unsigned short val_periodic_index,
                                    unsigned short commType) {

  /*--- Check for dummy communication. ---*/

  if (commType == PERIODIC_NONE) return;

  /*--- Local variables ---*/

  bool boundary_i, boundary_j;
  bool weighted = true;

  unsigned short iVar, jVar, iDim;
  unsigned short nNeighbor       = 0;
  unsigned short COUNT_PER_POINT = 0;
  unsigned short MPI_TYPE        = 0;
  unsigned short ICOUNT          = nVar;
  unsigned short JCOUNT          = nVar;

  int iMessage, iSend, nSend;

  unsigned long iPoint, msg_offset, buf_offset, iPeriodic;

  su2double *Diff      = new su2double[nVar];
  su2double *Und_Lapl  = new su2double[nVar];
  su2double *Sol_Min   = new su2double[nPrimVarGrad];
  su2double *Sol_Max   = new su2double[nPrimVarGrad];
  su2double *rotPrim_i = new su2double[nPrimVar];
  su2double *rotPrim_j = new su2double[nPrimVar];

  su2double Sensor_i = 0.0, Sensor_j = 0.0, Pressure_i, Pressure_j;
  const su2double *Coord_i, *Coord_j;
  su2double r11, r12, r13, r22, r23_a, r23_b, r33, weight;
  const su2double *center, *angles, *trans;
  su2double rotMatrix2D[2][2] = {{1.0,0.0},{0.0,1.0}};
  su2double rotMatrix3D[3][3] = {{1.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};
  su2double rotCoord_i[3] = {0.0}, rotCoord_j[3] = {0.0};
  su2double translation[3] = {0.0}, distance[3] = {0.0};
  const su2double zeros[3] = {0.0};
  su2activematrix Cvector;

  auto Rotate = [&](const su2double* origin, const su2double* direction, su2double* rotated) {
    if(nDim==2) GeometryToolbox::Rotate(rotMatrix2D, origin, direction, rotated);
    else GeometryToolbox::Rotate(rotMatrix3D, origin, direction, rotated);
  };

  string Marker_Tag;

  /*--- Set the size of the data packet and type depending on quantity. ---*/

  GetPeriodicCommCountAndType(config, commType, COUNT_PER_POINT, MPI_TYPE, ICOUNT, JCOUNT);

  /*--- Allocate buffers for matrices that need rotation. ---*/

  su2activematrix jacBlock(ICOUNT,JCOUNT);
  su2activematrix rotBlock(ICOUNT,JCOUNT);

  /*--- Check to make sure we have created a large enough buffer
   for these comms during preprocessing. It will be reallocated whenever
   we find a larger count per point than currently exists. After the
   first cycle of comms, this should be inactive. ---*/

  geometry->AllocatePeriodicComms(COUNT_PER_POINT);

  /*--- Set some local pointers to make access simpler. ---*/

  su2double *bufDSend = geometry->bufD_PeriodicSend;

  unsigned short *bufSSend = geometry->bufS_PeriodicSend;

  /*--- Handle the different types of gradient and limiter. ---*/

  auto& gradient = PeriodicCommHelpers::selectGradient(base_nodes, commType);
  auto& limiter = PeriodicCommHelpers::selectLimiter(base_nodes, commType);
  auto& field = PeriodicCommHelpers::selectField(base_nodes, commType);

  /*--- Load the specified quantity from the solver into the generic
   communication buffer in the geometry class. ---*/

  if (geometry->nPeriodicSend > 0) {

    /*--- Post all non-blocking recvs first before sends. ---*/

    geometry->PostPeriodicRecvs(geometry, config, MPI_TYPE, COUNT_PER_POINT);

    for (iMessage = 0; iMessage < geometry->nPeriodicSend; iMessage++) {

      /*--- Get the offset in the buffer for the start of this message. ---*/

      msg_offset = geometry->nPoint_PeriodicSend[iMessage];

      /*--- Get the number of periodic points we need to
       communicate on the current periodic marker. ---*/

      nSend = (geometry->nPoint_PeriodicSend[iMessage+1] -
               geometry->nPoint_PeriodicSend[iMessage]);

      SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
      for (iSend = 0; iSend < nSend; iSend++) {

        /*--- Get the local index for this communicated data. We need
         both the node and periodic face index (for rotations). ---*/

        iPoint    = geometry->Local_Point_PeriodicSend[msg_offset  + iSend];
        iPeriodic = geometry->Local_Marker_PeriodicSend[msg_offset + iSend];

        /*--- Retrieve the supplied periodic information. ---*/

        Marker_Tag = config->GetMarker_All_TagBound(iPeriodic);
        center     = config->GetPeriodicRotCenter(Marker_Tag);
        angles     = config->GetPeriodicRotAngles(Marker_Tag);
        trans      = config->GetPeriodicTranslation(Marker_Tag);

        /*--- Store (center+trans) as it is constant and will be added. ---*/

        translation[0] = center[0] + trans[0];
        translation[1] = center[1] + trans[1];
        translation[2] = center[2] + trans[2];

        /*--- Store angles separately for clarity. Compute sines/cosines. ---*/

        su2double Theta = angles[0];
        su2double Phi = angles[1];
        su2double Psi = angles[2];

        /*--- Compute the rotation matrix. Note that the implicit
         ordering is rotation about the x-axis, y-axis, then z-axis. ---*/

        if (nDim==2) {
          GeometryToolbox::RotationMatrix(Psi, rotMatrix2D);
        } else {
          GeometryToolbox::RotationMatrix(Theta, Phi, Psi, rotMatrix3D);
        }

        /*--- Compute the offset in the recv buffer for this point. ---*/

        buf_offset = (msg_offset + iSend)*COUNT_PER_POINT;

        /*--- Load the send buffers depending on the particular value
         that has been requested for communication. ---*/

        switch (commType) {

          case PERIODIC_VOLUME:

            /*--- Load the volume of the current periodic CV so that
             we can accumulate the total control volume size on all
             periodic faces. ---*/

            bufDSend[buf_offset] = geometry->nodes->GetVolume(iPoint) +
            geometry->nodes->GetPeriodicVolume(iPoint);

            break;

          case PERIODIC_NEIGHBORS:

            nNeighbor = 0;
            for (auto jPoint : geometry->nodes->GetPoints(iPoint)) {

              /*--- Check if this neighbor lies on the periodic face so
               that we avoid double counting neighbors on both sides. If
               not, increment the count of neighbors for the donor. ---*/

              if (!geometry->nodes->GetPeriodicBoundary(jPoint))
                nNeighbor++;
            }

            /*--- Store the number of neighbors in bufffer. ---*/

            bufSSend[buf_offset] = nNeighbor;

            break;

          case PERIODIC_RESIDUAL:

            /*--- Communicate the residual from our partial control
             volume to the other side of the periodic face. ---*/

            for (iVar = 0; iVar < nVar; iVar++) {
              bufDSend[buf_offset+iVar] = LinSysRes(iPoint, iVar);
            }

            /*--- Rotate the momentum components of the residual array. ---*/

            if (rotate_periodic) {
              Rotate(zeros, &LinSysRes(iPoint,1), &bufDSend[buf_offset+1]);
            }
            buf_offset += nVar;

            /*--- Load the time step for the current point. ---*/

            bufDSend[buf_offset] = base_nodes->GetDelta_Time(iPoint);
            buf_offset++;

            /*--- For implicit calculations, we will communicate the
             contributions to the Jacobian block diagonal, i.e., the
             impact of the point upon itself, J_ii. ---*/

            if (implicit_periodic) {

              for (iVar = 0; iVar < nVar; iVar++) {
                for (jVar = 0; jVar < nVar; jVar++) {
                  jacBlock[iVar][jVar] = Jacobian.GetBlock(iPoint, iPoint, iVar, jVar);
                }
              }

              /*--- Rotate the momentum columns of the Jacobian. ---*/

              if (rotate_periodic) {
                for (iVar = 0; iVar < nVar; iVar++) {
                  if (nDim == 2) {
                    jacBlock[1][iVar] = (rotMatrix2D[0][0]*Jacobian.GetBlock(iPoint, iPoint, 1, iVar) +
                                         rotMatrix2D[0][1]*Jacobian.GetBlock(iPoint, iPoint, 2, iVar));
                    jacBlock[2][iVar] = (rotMatrix2D[1][0]*Jacobian.GetBlock(iPoint, iPoint, 1, iVar) +
                                         rotMatrix2D[1][1]*Jacobian.GetBlock(iPoint, iPoint, 2, iVar));
                  } else {

                    jacBlock[1][iVar] = (rotMatrix3D[0][0]*Jacobian.GetBlock(iPoint, iPoint, 1, iVar) +
                                         rotMatrix3D[0][1]*Jacobian.GetBlock(iPoint, iPoint, 2, iVar) +
                                         rotMatrix3D[0][2]*Jacobian.GetBlock(iPoint, iPoint, 3, iVar));
                    jacBlock[2][iVar] = (rotMatrix3D[1][0]*Jacobian.GetBlock(iPoint, iPoint, 1, iVar) +
                                         rotMatrix3D[1][1]*Jacobian.GetBlock(iPoint, iPoint, 2, iVar) +
                                         rotMatrix3D[1][2]*Jacobian.GetBlock(iPoint, iPoint, 3, iVar));
                    jacBlock[3][iVar] = (rotMatrix3D[2][0]*Jacobian.GetBlock(iPoint, iPoint, 1, iVar) +
                                         rotMatrix3D[2][1]*Jacobian.GetBlock(iPoint, iPoint, 2, iVar) +
                                         rotMatrix3D[2][2]*Jacobian.GetBlock(iPoint, iPoint, 3, iVar));
                  }
                }
              }

              /*--- Load the Jacobian terms into the buffer for sending. ---*/

              for (iVar = 0; iVar < nVar; iVar++) {
                for (jVar = 0; jVar < nVar; jVar++) {
                  bufDSend[buf_offset] = jacBlock[iVar][jVar];
                  buf_offset++;
                }
              }
            }

            break;

          case PERIODIC_IMPLICIT:

            /*--- Communicate the solution from our master set of periodic
             nodes (from the linear solver perspective) to the passive
             periodic nodes on the matching face. This is done at the
             end of the iteration to synchronize the solution after the
             linear solve. ---*/

            for (iVar = 0; iVar < nVar; iVar++) {
              bufDSend[buf_offset+iVar] = base_nodes->GetSolution(iPoint, iVar);
            }

            /*--- Rotate the momentum components of the solution array. ---*/

            if (rotate_periodic) {
              Rotate(zeros, &base_nodes->GetSolution(iPoint)[1], &bufDSend[buf_offset+1]);
            }

            break;

          case PERIODIC_LAPLACIAN:

            /*--- For JST, the undivided Laplacian must be computed
             consistently by using the complete control volume info
             from both sides of the periodic face. ---*/

            for (iVar = 0; iVar < nVar; iVar++)
              Und_Lapl[iVar] = 0.0;

            for (auto jPoint : geometry->nodes->GetPoints(iPoint)) {

              /*--- Avoid periodic boundary points so that we do not
               duplicate edges on both sides of the periodic BC. ---*/

              if (!geometry->nodes->GetPeriodicBoundary(jPoint)) {

                /*--- Solution differences ---*/

                for (iVar = 0; iVar < nVar; iVar++)
                Diff[iVar] = (base_nodes->GetSolution(iPoint, iVar) -
                              base_nodes->GetSolution(jPoint,iVar));

                boundary_i = geometry->nodes->GetPhysicalBoundary(iPoint);
                boundary_j = geometry->nodes->GetPhysicalBoundary(jPoint);

                /*--- Both points inside the domain, or both in the boundary ---*/
                /*--- iPoint inside the domain, jPoint on the boundary ---*/

                if (!(boundary_i && !boundary_j)) {
                  if (geometry->nodes->GetDomain(iPoint)){
                    for (iVar = 0; iVar< nVar; iVar++)
                    Und_Lapl[iVar] -= Diff[iVar];
                  }
                }
              }
            }

            /*--- Store the components to be communicated in the buffer. ---*/

            for (iVar = 0; iVar < nVar; iVar++)
              bufDSend[buf_offset+iVar] = Und_Lapl[iVar];

            /*--- Rotate the momentum components of the Laplacian. ---*/

            if (rotate_periodic) {
              Rotate(zeros, &Und_Lapl[1], &bufDSend[buf_offset+1]);
            }

            break;

          case PERIODIC_MAX_EIG:

            /*--- Simple summation of eig calc on both periodic faces. ---*/

            bufDSend[buf_offset] = base_nodes->GetLambda(iPoint);

            break;

          case PERIODIC_SENSOR:

            /*--- For the centered schemes, the sensor must be computed
             consistently using info from the entire control volume
             on both sides of the periodic face. ---*/

            Sensor_i = 0.0; Sensor_j = 0.0;
            for (auto jPoint : geometry->nodes->GetPoints(iPoint)) {

              /*--- Avoid halos and boundary points so that we don't
               duplicate edges on both sides of the periodic BC. ---*/

              if (!geometry->nodes->GetPeriodicBoundary(jPoint)) {

                /*--- Use density instead of pressure for incomp. flows. ---*/

                  Pressure_i = base_nodes->GetPressure(iPoint);
                  Pressure_j = base_nodes->GetPressure(jPoint);

                boundary_i = geometry->nodes->GetPhysicalBoundary(iPoint);
                boundary_j = geometry->nodes->GetPhysicalBoundary(jPoint);

                /*--- Both points inside domain, or both on boundary ---*/
                /*--- iPoint inside the domain, jPoint on the boundary ---*/

                if (!(boundary_i && !boundary_j)) {
                  if (geometry->nodes->GetDomain(iPoint)) {
                    Sensor_i += (Pressure_j - Pressure_i);
                    Sensor_j += (Pressure_i + Pressure_j);
                  }
                }

              }
            }

            /*--- Store the sensor increments to buffer. After summing
             all contributions, these will be divided. ---*/

            bufDSend[buf_offset] = Sensor_i;
            buf_offset++;
            bufDSend[buf_offset] = Sensor_j;

            break;

          case PERIODIC_SOL_GG:
          case PERIODIC_SOL_GG_R:
          case PERIODIC_PRIM_GG:
          case PERIODIC_PRIM_GG_R:

            /*--- Access and rotate the partial G-G gradient. These will be
             summed on both sides of the periodic faces before dividing
             by the volume to complete the Green-Gauss gradient calc. ---*/

            for (iVar = 0; iVar < ICOUNT; iVar++) {
              for (iDim = 0; iDim < nDim; iDim++) {
                jacBlock[iVar][iDim] = gradient(iPoint, iVar, iDim);
              }
            }

            /*--- Rotate the gradients in x,y,z space for all variables. ---*/

            for (iVar = 0; iVar < ICOUNT; iVar++) {
              Rotate(zeros, jacBlock[iVar], rotBlock[iVar]);
            }

            /*--- Rotate the vector components of the solution. ---*/

            if (rotate_periodic) {
              for (iDim = 0; iDim < nDim; iDim++) {
                su2double d_diDim[3] = {0.0};
                for (iVar = 1; iVar < 1+nDim; ++iVar) {
                  d_diDim[iVar-1] = rotBlock(iVar, iDim);
                }
                su2double rotated[3] = {0.0};
                Rotate(zeros, d_diDim, rotated);
                for (iVar = 1; iVar < 1+nDim; ++iVar) {
                  rotBlock(iVar, iDim) = rotated[iVar-1];
                }
              }
            }

            /*--- Store the partial gradient in the buffer. ---*/

            for (iVar = 0; iVar < ICOUNT; iVar++) {
              for (iDim = 0; iDim < nDim; iDim++) {
                bufDSend[buf_offset+iVar*nDim+iDim] = rotBlock[iVar][iDim];
              }
            }

            break;

          case PERIODIC_SOL_LS: case PERIODIC_SOL_ULS:
          case PERIODIC_SOL_LS_R: case PERIODIC_SOL_ULS_R:
          case PERIODIC_PRIM_LS: case PERIODIC_PRIM_ULS:
          case PERIODIC_PRIM_LS_R: case PERIODIC_PRIM_ULS_R:

            /*--- For L-S gradient calculations with rotational periodicity,
             we will need to rotate the x,y,z components. To make the process
             easier, we choose to rotate the initial periodic point and their
             neighbor points into their location on the donor marker before
             computing the terms that we need to communicate. ---*/

            /*--- Set a flag for unweighted or weighted least-squares. ---*/

            switch(commType) {
              case PERIODIC_SOL_ULS:
              case PERIODIC_SOL_ULS_R:
              case PERIODIC_PRIM_ULS:
              case PERIODIC_PRIM_ULS_R:
                weighted = false;
                break;
              default:
                weighted = true;
                break;
            }

            /*--- Get coordinates for the current point. ---*/

            Coord_i = geometry->nodes->GetCoord(iPoint);

            /*--- Get the position vector from rotation center to point. ---*/

            GeometryToolbox::Distance(nDim, Coord_i, center, distance);

            /*--- Compute transformed point coordinates. ---*/

            Rotate(translation, distance, rotCoord_i);

            /*--- Get conservative solution and rotate if necessary. ---*/

            for (iVar = 0; iVar < ICOUNT; iVar++)
              rotPrim_i[iVar] = field(iPoint, iVar);

            if (rotate_periodic) {
              Rotate(zeros, &field(iPoint,1), &rotPrim_i[1]);
            }

            /*--- Inizialization of variables ---*/

            Cvector.resize(ICOUNT,nDim) = su2double(0.0);

            r11 = 0.0;   r12 = 0.0;   r22 = 0.0;
            r13 = 0.0; r23_a = 0.0; r23_b = 0.0;  r33 = 0.0;

            for (auto jPoint : geometry->nodes->GetPoints(iPoint)) {

              /*--- Avoid periodic boundary points so that we do not
               duplicate edges on both sides of the periodic BC. ---*/

              if (!geometry->nodes->GetPeriodicBoundary(jPoint)) {

                /*--- Get coordinates for the neighbor point. ---*/

                Coord_j = geometry->nodes->GetCoord(jPoint);

                /*--- Get the position vector from rotation center. ---*/

                GeometryToolbox::Distance(nDim, Coord_j, center, distance);

                /*--- Compute transformed point coordinates. ---*/

                Rotate(translation, distance, rotCoord_j);

                /*--- Get conservative solution and rotate if necessary. ---*/

                for (iVar = 0; iVar < ICOUNT; iVar++)
                  rotPrim_j[iVar] = field(jPoint,iVar);

                if (rotate_periodic) {
                  Rotate(zeros, &field(jPoint,1), &rotPrim_j[1]);
                }

                if (weighted) {
                  weight = GeometryToolbox::SquaredDistance(nDim, rotCoord_j, rotCoord_i);
                } else {
                  weight = 1.0;
                }

                /*--- Sumations for entries of upper triangular matrix R ---*/

                if (weight != 0.0) {

                  r11 += ((rotCoord_j[0]-rotCoord_i[0])*
                          (rotCoord_j[0]-rotCoord_i[0])/weight);
                  r12 += ((rotCoord_j[0]-rotCoord_i[0])*
                          (rotCoord_j[1]-rotCoord_i[1])/weight);
                  r22 += ((rotCoord_j[1]-rotCoord_i[1])*
                          (rotCoord_j[1]-rotCoord_i[1])/weight);

                  if (nDim == 3) {
                    r13   += ((rotCoord_j[0]-rotCoord_i[0])*
                              (rotCoord_j[2]-rotCoord_i[2])/weight);
                    r23_a += ((rotCoord_j[1]-rotCoord_i[1])*
                              (rotCoord_j[2]-rotCoord_i[2])/weight);
                    r23_b += ((rotCoord_j[0]-rotCoord_i[0])*
                              (rotCoord_j[2]-rotCoord_i[2])/weight);
                    r33   += ((rotCoord_j[2]-rotCoord_i[2])*
                              (rotCoord_j[2]-rotCoord_i[2])/weight);
                  }

                  /*--- Entries of c:= transpose(A)*b ---*/

                  for (iVar = 0; iVar < ICOUNT; iVar++)
                  for (iDim = 0; iDim < nDim; iDim++)
                  Cvector(iVar,iDim) += ((rotCoord_j[iDim]-rotCoord_i[iDim])*
                                          (rotPrim_j[iVar]-rotPrim_i[iVar])/weight);

                }
              }
            }

            /*--- We store and communicate the increments for the matching
             upper triangular matrix (weights) and the r.h.s. vector.
             These will be accumulated before completing the L-S gradient
             calculation for each periodic point. ---*/

            if (nDim == 2) {
              bufDSend[buf_offset] = r11;   buf_offset++;
              bufDSend[buf_offset] = r12;   buf_offset++;
              bufDSend[buf_offset] = 0.0;   buf_offset++;
              bufDSend[buf_offset] = r22;   buf_offset++;
            }
            if (nDim == 3) {
              bufDSend[buf_offset] = r11;   buf_offset++;
              bufDSend[buf_offset] = r12;   buf_offset++;
              bufDSend[buf_offset] = r13;   buf_offset++;

              bufDSend[buf_offset] = 0.0;   buf_offset++;
              bufDSend[buf_offset] = r22;   buf_offset++;
              bufDSend[buf_offset] = r23_a; buf_offset++;

              bufDSend[buf_offset] = 0.0;   buf_offset++;
              bufDSend[buf_offset] = r23_b; buf_offset++;
              bufDSend[buf_offset] = r33;   buf_offset++;
            }

            for (iVar = 0; iVar < ICOUNT; iVar++) {
              for (iDim = 0; iDim < nDim; iDim++) {
                bufDSend[buf_offset] = Cvector(iVar,iDim);
                buf_offset++;
              }
            }

            break;

          case PERIODIC_LIM_PRIM_1:
          case PERIODIC_LIM_SOL_1:

            /*--- The first phase of the periodic limiter calculation
             ensures that the proper min and max of the solution are found
             among all nodes adjacent to periodic faces. ---*/

            /*--- We send the min and max over "our" neighbours. ---*/

            for (iVar = 0; iVar < ICOUNT; iVar++) {
              Sol_Min[iVar] = base_nodes->GetSolution_Min()(iPoint, iVar);
              Sol_Max[iVar] = base_nodes->GetSolution_Max()(iPoint, iVar);
            }

            for (auto jPoint : geometry->nodes->GetPoints(iPoint)) {
              for (iVar = 0; iVar < ICOUNT; iVar++) {
                Sol_Min[iVar] = min(Sol_Min[iVar], field(jPoint, iVar));
                Sol_Max[iVar] = max(Sol_Max[iVar], field(jPoint, iVar));
              }
            }

            for (iVar = 0; iVar < ICOUNT; iVar++) {
              bufDSend[buf_offset+iVar]        = Sol_Min[iVar];
              bufDSend[buf_offset+ICOUNT+iVar] = Sol_Max[iVar];
            }

            /*--- Rotate the momentum components of the min/max. ---*/

            if (rotate_periodic) {
              Rotate(zeros, &Sol_Min[1], &bufDSend[buf_offset+1]);
              Rotate(zeros, &Sol_Max[1], &bufDSend[buf_offset+ICOUNT+1]);
            }

            break;

          case PERIODIC_LIM_PRIM_2:
          case PERIODIC_LIM_SOL_2:

            /*--- The second phase of the periodic limiter calculation
             ensures that the correct minimum value of the limiter is
             found for a node on a periodic face and stores it. ---*/

            for (iVar = 0; iVar < ICOUNT; iVar++) {
              bufDSend[buf_offset+iVar] = limiter(iPoint, iVar);
            }

            if (rotate_periodic) {
              Rotate(zeros, &limiter(iPoint,1), &bufDSend[buf_offset+1]);
            }

            break;

          default:
            SU2_MPI::Error("Unrecognized quantity for periodic communication.",
                           CURRENT_FUNCTION);
            break;
        }
      }
      END_SU2_OMP_FOR

      /*--- Launch the point-to-point MPI send for this message. ---*/

      geometry->PostPeriodicSends(geometry, config, MPI_TYPE, COUNT_PER_POINT, iMessage);

    }
  }

  delete [] Diff;
  delete [] Und_Lapl;
  delete [] Sol_Min;
  delete [] Sol_Max;
  delete [] rotPrim_i;
  delete [] rotPrim_j;

}

void CSolver::CompletePeriodicComms(CGeometry *geometry,
                                    const CConfig *config,
                                    unsigned short val_periodic_index,
                                    unsigned short commType) {

  /*--- Check for dummy communication. ---*/

  if (commType == PERIODIC_NONE) return;

  /*--- Set the size of the data packet and type depending on quantity. ---*/

  unsigned short COUNT_PER_POINT = 0, MPI_TYPE = 0, ICOUNT = 0, JCOUNT = 0;
  GetPeriodicCommCountAndType(config, commType, COUNT_PER_POINT, MPI_TYPE, ICOUNT, JCOUNT);

  /*--- Local variables ---*/

  unsigned short nPeriodic = config->GetnMarker_Periodic();
  unsigned short iDim, jDim, iVar, jVar, iPeriodic, nNeighbor;

  unsigned long iPoint, iRecv, nRecv, msg_offset, buf_offset, total_index;

  int source, iMessage, jRecv;

  /*--- Status is global so all threads can see the result of Waitany. ---*/
  static SU2_MPI::Status status;

  su2double *Diff = new su2double[nVar];

  su2double Time_Step, Volume;

  su2double **Jacobian_i = nullptr;
  if ((commType == PERIODIC_RESIDUAL) && implicit_periodic) {
    Jacobian_i = new su2double* [nVar];
    for (iVar = 0; iVar < nVar; iVar++)
      Jacobian_i[iVar] = new su2double [nVar];
  }

  /*--- Set some local pointers to make access simpler. ---*/

  const su2double *bufDRecv = geometry->bufD_PeriodicRecv;

  const unsigned short *bufSRecv = geometry->bufS_PeriodicRecv;

  /*--- Handle the different types of gradient and limiter. ---*/

  auto& gradient = PeriodicCommHelpers::selectGradient(base_nodes, commType);
  auto& limiter = PeriodicCommHelpers::selectLimiter(base_nodes, commType);

  /*--- Store the data that was communicated into the appropriate
   location within the local class data structures. ---*/

  if (geometry->nPeriodicRecv > 0) {

    for (iMessage = 0; iMessage < geometry->nPeriodicRecv; iMessage++) {

      /*--- For efficiency, recv the messages dynamically based on
       the order they arrive. ---*/

#ifdef HAVE_MPI
      /*--- Once we have recv'd a message, get the source rank. ---*/
      int ind;
      SU2_OMP_MASTER
      SU2_MPI::Waitany(geometry->nPeriodicRecv,
                       geometry->req_PeriodicRecv,
                       &ind, &status);
      END_SU2_OMP_MASTER
      SU2_OMP_BARRIER
      source = status.MPI_SOURCE;
#else
      /*--- For serial calculations, we know the rank. ---*/
      source = rank;
      SU2_OMP_BARRIER
#endif

      /*--- We know the offsets based on the source rank. ---*/

      jRecv = geometry->PeriodicRecv2Neighbor[source];

      /*--- Get the offset in the buffer for the start of this message. ---*/

      msg_offset = geometry->nPoint_PeriodicRecv[jRecv];

      /*--- Get the number of packets to be received in this message. ---*/

      nRecv = (geometry->nPoint_PeriodicRecv[jRecv+1] -
               geometry->nPoint_PeriodicRecv[jRecv]);

      SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
      for (iRecv = 0; iRecv < nRecv; iRecv++) {

        /*--- Get the local index for this communicated data. ---*/

        iPoint    = geometry->Local_Point_PeriodicRecv[msg_offset  + iRecv];
        iPeriodic = geometry->Local_Marker_PeriodicRecv[msg_offset + iRecv];

        /*--- While all periodic face data was accumulated, we only store
         the values for the current pair of periodic faces. This is slightly
         inefficient when we have multiple pairs of periodic faces, but
         it simplifies the communications. ---*/

        if ((iPeriodic == val_periodic_index) ||
            (iPeriodic == val_periodic_index + nPeriodic/2)) {

          /*--- Compute the offset in the recv buffer for this point. ---*/

          buf_offset = (msg_offset + iRecv)*COUNT_PER_POINT;

          /*--- Store the data correctly depending on the quantity. ---*/

          switch (commType) {

            case PERIODIC_VOLUME:

              /*--- The periodic points need to keep track of their
               total volume spread across the periodic faces. ---*/

              Volume = (bufDRecv[buf_offset] +
                        geometry->nodes->GetPeriodicVolume(iPoint));
              geometry->nodes->SetPeriodicVolume(iPoint, Volume);

              break;

            case PERIODIC_NEIGHBORS:

              /*--- Store the extra neighbors on the periodic face. ---*/

              nNeighbor = (geometry->nodes->GetnNeighbor(iPoint) +
                           bufSRecv[buf_offset]);
              geometry->nodes->SetnNeighbor(iPoint, nNeighbor);

              break;

            case PERIODIC_RESIDUAL:

              /*--- Add contributions to total residual. ---*/

              LinSysRes.AddBlock(iPoint, &bufDRecv[buf_offset]);
              buf_offset += nVar;

              /*--- Check the computed time step against the donor
               value and keep the minimum in order to be conservative. ---*/

              Time_Step = base_nodes->GetDelta_Time(iPoint);
              if (bufDRecv[buf_offset] < Time_Step)
                base_nodes->SetDelta_Time(iPoint,bufDRecv[buf_offset]);
              buf_offset++;

              /*--- For implicit integration, we choose the first
               periodic face of each pair to be the master/owner of
               the solution for the linear system while fixing the
               solution at the matching face during the solve. Here,
               we remove the Jacobian and residual contributions from
               the passive face such that it does not participate in
               the linear solve. ---*/

              if (implicit_periodic) {

                for (iVar = 0; iVar < nVar; iVar++) {
                  for (jVar = 0; jVar < nVar; jVar++) {
                    Jacobian_i[iVar][jVar] = bufDRecv[buf_offset];
                    buf_offset++;
                  }
                }

                Jacobian.AddBlock2Diag(iPoint, Jacobian_i);

                if (iPeriodic == val_periodic_index + nPeriodic/2) {
                  for (iVar = 0; iVar < nVar; iVar++) {
                    LinSysRes(iPoint, iVar) = 0.0;
                    total_index = iPoint*nVar+iVar;
                    Jacobian.DeleteValsRowi(total_index);
                  }
                }

              }

              break;

            case PERIODIC_IMPLICIT:

              /*--- For implicit integration, we choose the first
               periodic face of each pair to be the master/owner of
               the solution for the linear system while fixing the
               solution at the matching face during the solve. Here,
               we are updating the solution at the passive nodes
               using the new solution from the master. ---*/

              if ((implicit_periodic) &&
                  (iPeriodic == val_periodic_index + nPeriodic/2)) {

                /*--- Directly set the solution on the passive periodic
                 face that is provided from the master. ---*/

                for (iVar = 0; iVar < nVar; iVar++) {
                  base_nodes->SetSolution(iPoint, iVar, bufDRecv[buf_offset]);
                  base_nodes->SetSolution_Old(iPoint, iVar, bufDRecv[buf_offset]);
                  buf_offset++;
                }

              }

              break;

            case PERIODIC_LAPLACIAN:

              /*--- Adjust the undivided Laplacian. The accumulation was
               with a subtraction before communicating, so now just add. ---*/

              for (iVar = 0; iVar < nVar; iVar++)
                base_nodes->AddUnd_Lapl(iPoint, iVar, bufDRecv[buf_offset+iVar]);

              break;

            case PERIODIC_MAX_EIG:

              /*--- Simple accumulation of the max eig on periodic faces. ---*/

              base_nodes->AddLambda(iPoint,bufDRecv[buf_offset]);

              break;

            case PERIODIC_SENSOR:

              /*--- Simple accumulation of the sensors on periodic faces. ---*/

              iPoint_UndLapl[iPoint] += bufDRecv[buf_offset]; buf_offset++;
              jPoint_UndLapl[iPoint] += bufDRecv[buf_offset];

              break;

            case PERIODIC_SOL_GG:
            case PERIODIC_SOL_GG_R:
            case PERIODIC_PRIM_GG:
            case PERIODIC_PRIM_GG_R:

              /*--- For G-G, we accumulate partial gradients then compute
               the final value using the entire volume of the periodic cell. ---*/

              for (iVar = 0; iVar < ICOUNT; iVar++)
                for (iDim = 0; iDim < nDim; iDim++)
                  gradient(iPoint, iVar, iDim) += bufDRecv[buf_offset+iVar*nDim+iDim];

              break;

            case PERIODIC_SOL_LS: case PERIODIC_SOL_ULS:
            case PERIODIC_SOL_LS_R: case PERIODIC_SOL_ULS_R:
            case PERIODIC_PRIM_LS: case PERIODIC_PRIM_ULS:
            case PERIODIC_PRIM_LS_R: case PERIODIC_PRIM_ULS_R:

              /*--- For L-S, we build the upper triangular matrix and the
               r.h.s. vector by accumulating from all periodic partial
               control volumes. ---*/

              for (iDim = 0; iDim < nDim; iDim++) {
                for (jDim = 0; jDim < nDim; jDim++) {
                  base_nodes->AddRmatrix(iPoint, iDim,jDim,bufDRecv[buf_offset]);
                  buf_offset++;
                }
              }
              for (iVar = 0; iVar < ICOUNT; iVar++) {
                for (iDim = 0; iDim < nDim; iDim++) {
                  gradient(iPoint, iVar, iDim) += bufDRecv[buf_offset];
                  buf_offset++;
                }
              }

              break;

            case PERIODIC_LIM_PRIM_1:
            case PERIODIC_LIM_SOL_1:

              /*--- Update solution min/max with min/max between "us" and
               the periodic match plus its neighbors, computation will need to
               be concluded on "our" side to account for "our" neighbors. ---*/

              for (iVar = 0; iVar < ICOUNT; iVar++) {

                /*--- Solution minimum. ---*/

                su2double Solution_Min = min(base_nodes->GetSolution_Min()(iPoint, iVar),
                                             bufDRecv[buf_offset+iVar]);
                base_nodes->GetSolution_Min()(iPoint, iVar) = Solution_Min;

                /*--- Solution maximum. ---*/

                su2double Solution_Max = max(base_nodes->GetSolution_Max()(iPoint, iVar),
                                             bufDRecv[buf_offset+ICOUNT+iVar]);
                base_nodes->GetSolution_Max()(iPoint, iVar) = Solution_Max;
              }

              break;

            case PERIODIC_LIM_PRIM_2:
            case PERIODIC_LIM_SOL_2:

              /*--- Check the min values found on the matching periodic
               faces for the limiter, and store the proper min value. ---*/

              for (iVar = 0; iVar < ICOUNT; iVar++)
                limiter(iPoint, iVar) = min(limiter(iPoint, iVar), bufDRecv[buf_offset+iVar]);

              break;

            default:

              SU2_MPI::Error("Unrecognized quantity for periodic communication.",
                             CURRENT_FUNCTION);
              break;

          }
        }
      }
      END_SU2_OMP_FOR
    }

    /*--- Verify that all non-blocking point-to-point sends have finished.
     Note that this should be satisfied, as we have received all of the
     data in the loop above at this point. ---*/

#ifdef HAVE_MPI
    SU2_OMP_MASTER
    SU2_MPI::Waitall(geometry->nPeriodicSend,
                     geometry->req_PeriodicSend,
                     MPI_STATUS_IGNORE);
    END_SU2_OMP_MASTER
#endif
    SU2_OMP_BARRIER
  }

  delete [] Diff;

  if (Jacobian_i)
    for (iVar = 0; iVar < nVar; iVar++)
      delete [] Jacobian_i[iVar];
  delete [] Jacobian_i;

}

void CSolver::GetCommCountAndType(const CConfig* config,
                                  unsigned short commType,
                                  unsigned short &COUNT_PER_POINT,
                                  unsigned short &MPI_TYPE) const {
  switch (commType) {
    case SOLUTION:
    case SOLUTION_OLD:
    case UNDIVIDED_LAPLACIAN:
    case SOLUTION_LIMITER:
      COUNT_PER_POINT  = nVar;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      break;
    case MAX_EIGENVALUE:
    case SENSOR:
      COUNT_PER_POINT  = 1;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      break;
    case SOLUTION_GRADIENT:
    case SOLUTION_GRAD_REC:
      COUNT_PER_POINT  = nVar*nDim;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      break;
    case PRIMITIVE_GRADIENT:
    case PRIMITIVE_GRAD_REC:
      COUNT_PER_POINT  = nPrimVarGrad*nDim;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      break;
    case PRIMITIVE_LIMITER:
      COUNT_PER_POINT  = nPrimVarGrad;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      break;
    case SOLUTION_EDDY:
      COUNT_PER_POINT  = nVar+1;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      break;
    case SOLUTION_FEA:
      if (config->GetTime_Domain())
        COUNT_PER_POINT  = nVar*3;
      else
        COUNT_PER_POINT  = nVar;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      break;
    case AUXVAR_GRADIENT:
      COUNT_PER_POINT  = nDim*base_nodes->GetnAuxVar();
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      break;
    case MESH_DISPLACEMENTS:
      COUNT_PER_POINT  = nDim;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      break;
    case SOLUTION_TIME_N:
      COUNT_PER_POINT  = nVar;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      break;
    case SOLUTION_TIME_N1:
      COUNT_PER_POINT  = nVar;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      break;
    default:
      SU2_MPI::Error("Unrecognized quantity for point-to-point MPI comms.",
                     CURRENT_FUNCTION);
      break;
  }
}

namespace CommHelpers {
  CVectorOfMatrix& selectGradient(CVariable* nodes, unsigned short commType) {
    switch(commType) {
      case SOLUTION_GRAD_REC: return nodes->GetGradient_Reconstruction();
      case PRIMITIVE_GRADIENT: return nodes->GetGradient_Primitive();
      case PRIMITIVE_GRAD_REC: return nodes->GetGradient_Reconstruction();
      case AUXVAR_GRADIENT: return nodes->GetAuxVarGradient();
      default: return nodes->GetGradient();
    }
  }

  su2activematrix& selectLimiter(CVariable* nodes, unsigned short commType) {
    if (commType == PRIMITIVE_LIMITER) return nodes->GetLimiter_Primitive();
    return nodes->GetLimiter();
  }
}

void CSolver::InitiateComms(CGeometry *geometry,
                            const CConfig *config,
                            unsigned short commType) {

  /*--- Local variables ---*/

  unsigned short iVar, iDim;
  unsigned short COUNT_PER_POINT = 0;
  unsigned short MPI_TYPE        = 0;

  unsigned long iPoint, msg_offset, buf_offset;

  int iMessage, iSend, nSend;

  /*--- Set the size of the data packet and type depending on quantity. ---*/

  GetCommCountAndType(config, commType, COUNT_PER_POINT, MPI_TYPE);

  /*--- Check to make sure we have created a large enough buffer
   for these comms during preprocessing. This is only for the su2double
   buffer. It will be reallocated whenever we find a larger count
   per point. After the first cycle of comms, this should be inactive. ---*/

  geometry->AllocateP2PComms(COUNT_PER_POINT);

  /*--- Set some local pointers to make access simpler. ---*/

  su2double *bufDSend = geometry->bufD_P2PSend;

  /*--- Handle the different types of gradient and limiter. ---*/

  const auto nVarGrad = COUNT_PER_POINT / nDim;
  auto& gradient = CommHelpers::selectGradient(base_nodes, commType);
  auto& limiter = CommHelpers::selectLimiter(base_nodes, commType);

  /*--- Load the specified quantity from the solver into the generic
   communication buffer in the geometry class. ---*/

  if (geometry->nP2PSend > 0) {

    /*--- Post all non-blocking recvs first before sends. ---*/

    geometry->PostP2PRecvs(geometry, config, MPI_TYPE, COUNT_PER_POINT, false);

    for (iMessage = 0; iMessage < geometry->nP2PSend; iMessage++) {

      /*--- Get the offset in the buffer for the start of this message. ---*/

      msg_offset = geometry->nPoint_P2PSend[iMessage];

      /*--- Total count can include multiple pieces of data per element. ---*/

      nSend = (geometry->nPoint_P2PSend[iMessage+1] -
               geometry->nPoint_P2PSend[iMessage]);

      SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
      for (iSend = 0; iSend < nSend; iSend++) {

        /*--- Get the local index for this communicated data. ---*/

        iPoint = geometry->Local_Point_P2PSend[msg_offset + iSend];

        /*--- Compute the offset in the recv buffer for this point. ---*/

        buf_offset = (msg_offset + iSend)*COUNT_PER_POINT;

        switch (commType) {
          case SOLUTION:
            for (iVar = 0; iVar < nVar; iVar++)
              bufDSend[buf_offset+iVar] = base_nodes->GetSolution(iPoint, iVar);
            break;
          case SOLUTION_OLD:
            for (iVar = 0; iVar < nVar; iVar++)
              bufDSend[buf_offset+iVar] = base_nodes->GetSolution_Old(iPoint, iVar);
            break;
          case SOLUTION_EDDY:
            for (iVar = 0; iVar < nVar; iVar++)
              bufDSend[buf_offset+iVar] = base_nodes->GetSolution(iPoint, iVar);
            bufDSend[buf_offset+nVar]   = base_nodes->GetmuT(iPoint);
            break;
          case UNDIVIDED_LAPLACIAN:
            for (iVar = 0; iVar < nVar; iVar++)
              bufDSend[buf_offset+iVar] = base_nodes->GetUndivided_Laplacian(iPoint, iVar);
            break;
          case SOLUTION_LIMITER:
          case PRIMITIVE_LIMITER:
            for (iVar = 0; iVar < COUNT_PER_POINT; iVar++)
              bufDSend[buf_offset+iVar] = limiter(iPoint, iVar);
            break;
          case MAX_EIGENVALUE:
            bufDSend[buf_offset] = base_nodes->GetLambda(iPoint);
            break;
          case SENSOR:
            bufDSend[buf_offset] = base_nodes->GetSensor(iPoint);
            break;
          case SOLUTION_GRADIENT:
          case PRIMITIVE_GRADIENT:
          case SOLUTION_GRAD_REC:
          case PRIMITIVE_GRAD_REC:
          case AUXVAR_GRADIENT:
            for (iVar = 0; iVar < nVarGrad; iVar++)
              for (iDim = 0; iDim < nDim; iDim++)
                bufDSend[buf_offset+iVar*nDim+iDim] = gradient(iPoint, iVar, iDim);
            break;
          case SOLUTION_FEA:
            for (iVar = 0; iVar < nVar; iVar++) {
              bufDSend[buf_offset+iVar] = base_nodes->GetSolution(iPoint, iVar);
              if (config->GetTime_Domain()) {
                bufDSend[buf_offset+nVar+iVar]   = base_nodes->GetSolution_Vel(iPoint, iVar);
                bufDSend[buf_offset+nVar*2+iVar] = base_nodes->GetSolution_Accel(iPoint, iVar);
              }
            }
            break;
          case MESH_DISPLACEMENTS:
            for (iDim = 0; iDim < nDim; iDim++)
              bufDSend[buf_offset+iDim] = base_nodes->GetBound_Disp(iPoint, iDim);
            break;
          case SOLUTION_TIME_N:
            for (iVar = 0; iVar < nVar; iVar++)
              bufDSend[buf_offset+iVar] = base_nodes->GetSolution_time_n(iPoint, iVar);
            break;
          case SOLUTION_TIME_N1:
            for (iVar = 0; iVar < nVar; iVar++)
              bufDSend[buf_offset+iVar] = base_nodes->GetSolution_time_n1(iPoint, iVar);
            break;
          default:
            SU2_MPI::Error("Unrecognized quantity for point-to-point MPI comms.",
                           CURRENT_FUNCTION);
            break;
        }
      }
      END_SU2_OMP_FOR

      /*--- Launch the point-to-point MPI send for this message. ---*/

      geometry->PostP2PSends(geometry, config, MPI_TYPE, COUNT_PER_POINT, iMessage, false);

    }
  }

}

void CSolver::CompleteComms(CGeometry *geometry,
                            const CConfig *config,
                            unsigned short commType) {

  /*--- Local variables ---*/

  unsigned short iDim, iVar;
  unsigned long iPoint, iRecv, nRecv, msg_offset, buf_offset;
  unsigned short COUNT_PER_POINT = 0;
  unsigned short MPI_TYPE = 0;

  int ind, source, iMessage, jRecv;

  /*--- Global status so all threads can see the result of Waitany. ---*/
  static SU2_MPI::Status status;

  /*--- Set the size of the data packet and type depending on quantity. ---*/

  GetCommCountAndType(config, commType, COUNT_PER_POINT, MPI_TYPE);

  /*--- Set some local pointers to make access simpler. ---*/

  const su2double *bufDRecv = geometry->bufD_P2PRecv;

  /*--- Handle the different types of gradient and limiter. ---*/

  const auto nVarGrad = COUNT_PER_POINT / nDim;
  auto& gradient = CommHelpers::selectGradient(base_nodes, commType);
  auto& limiter = CommHelpers::selectLimiter(base_nodes, commType);

  /*--- Store the data that was communicated into the appropriate
   location within the local class data structures. ---*/

  if (geometry->nP2PRecv > 0) {

    for (iMessage = 0; iMessage < geometry->nP2PRecv; iMessage++) {

      /*--- For efficiency, recv the messages dynamically based on
       the order they arrive. ---*/

      SU2_OMP_MASTER
      SU2_MPI::Waitany(geometry->nP2PRecv, geometry->req_P2PRecv, &ind, &status);
      END_SU2_OMP_MASTER
      SU2_OMP_BARRIER

      /*--- Once we have recv'd a message, get the source rank. ---*/

      source = status.MPI_SOURCE;

      /*--- We know the offsets based on the source rank. ---*/

      jRecv = geometry->P2PRecv2Neighbor[source];

      /*--- Get the offset in the buffer for the start of this message. ---*/

      msg_offset = geometry->nPoint_P2PRecv[jRecv];

      /*--- Get the number of packets to be received in this message. ---*/

      nRecv = (geometry->nPoint_P2PRecv[jRecv+1] -
               geometry->nPoint_P2PRecv[jRecv]);

      SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
      for (iRecv = 0; iRecv < nRecv; iRecv++) {

        /*--- Get the local index for this communicated data. ---*/

        iPoint = geometry->Local_Point_P2PRecv[msg_offset + iRecv];

        /*--- Compute the offset in the recv buffer for this point. ---*/

        buf_offset = (msg_offset + iRecv)*COUNT_PER_POINT;

        /*--- Store the data correctly depending on the quantity. ---*/

        switch (commType) {
          case SOLUTION:
            for (iVar = 0; iVar < nVar; iVar++)
              base_nodes->SetSolution(iPoint, iVar, bufDRecv[buf_offset+iVar]);
            break;
          case SOLUTION_OLD:
            for (iVar = 0; iVar < nVar; iVar++)
              base_nodes->SetSolution_Old(iPoint, iVar, bufDRecv[buf_offset+iVar]);
            break;
          case SOLUTION_EDDY:
            for (iVar = 0; iVar < nVar; iVar++)
              base_nodes->SetSolution(iPoint, iVar, bufDRecv[buf_offset+iVar]);
            base_nodes->SetmuT(iPoint,bufDRecv[buf_offset+nVar]);
            break;
          case UNDIVIDED_LAPLACIAN:
            for (iVar = 0; iVar < nVar; iVar++)
              base_nodes->SetUnd_Lapl(iPoint, iVar, bufDRecv[buf_offset+iVar]);
            break;
          case SOLUTION_LIMITER:
          case PRIMITIVE_LIMITER:
            for (iVar = 0; iVar < COUNT_PER_POINT; iVar++)
              limiter(iPoint,iVar) = bufDRecv[buf_offset+iVar];
            break;
          case MAX_EIGENVALUE:
            base_nodes->SetLambda(iPoint,bufDRecv[buf_offset]);
            break;
          case SENSOR:
            base_nodes->SetSensor(iPoint,bufDRecv[buf_offset]);
            break;
          case SOLUTION_GRADIENT:
          case PRIMITIVE_GRADIENT:
          case SOLUTION_GRAD_REC:
          case PRIMITIVE_GRAD_REC:
          case AUXVAR_GRADIENT:
            for (iVar = 0; iVar < nVarGrad; iVar++)
              for (iDim = 0; iDim < nDim; iDim++)
                gradient(iPoint,iVar,iDim) = bufDRecv[buf_offset+iVar*nDim+iDim];
            break;
          case SOLUTION_FEA:
            for (iVar = 0; iVar < nVar; iVar++) {
              base_nodes->SetSolution(iPoint, iVar, bufDRecv[buf_offset+iVar]);
              if (config->GetTime_Domain()) {
                base_nodes->SetSolution_Vel(iPoint, iVar, bufDRecv[buf_offset+nVar+iVar]);
                base_nodes->SetSolution_Accel(iPoint, iVar, bufDRecv[buf_offset+nVar*2+iVar]);
              }
            }
            break;
          case MESH_DISPLACEMENTS:
            for (iDim = 0; iDim < nDim; iDim++)
              base_nodes->SetBound_Disp(iPoint, iDim, bufDRecv[buf_offset+iDim]);
            break;
          case SOLUTION_TIME_N:
            for (iVar = 0; iVar < nVar; iVar++)
              base_nodes->Set_Solution_time_n(iPoint, iVar, bufDRecv[buf_offset+iVar]);
            break;
          case SOLUTION_TIME_N1:
            for (iVar = 0; iVar < nVar; iVar++)
              base_nodes->Set_Solution_time_n1(iPoint, iVar, bufDRecv[buf_offset+iVar]);
            break;
          default:
            SU2_MPI::Error("Unrecognized quantity for point-to-point MPI comms.",
                           CURRENT_FUNCTION);
            break;
        }
      }
      END_SU2_OMP_FOR
    }

    /*--- Verify that all non-blocking point-to-point sends have finished.
     Note that this should be satisfied, as we have received all of the
     data in the loop above at this point. ---*/

#ifdef HAVE_MPI
    SU2_OMP_MASTER
    SU2_MPI::Waitall(geometry->nP2PSend, geometry->req_P2PSend, MPI_STATUS_IGNORE);
    END_SU2_OMP_MASTER
#endif
    SU2_OMP_BARRIER
  }

}

void CSolver::ResetCFLAdapt() {
  NonLinRes_Series.clear();
  Old_Func = 0;
  New_Func = 0;
  NonLinRes_Counter = 0;
}


void CSolver::AdaptCFLNumber(CGeometry **geometry,
                             CSolver   ***solver_container,
                             CConfig   *config) {

  /* Adapt the CFL number on all multigrid levels using an
   exponential progression with under-relaxation approach. */

  vector<su2double> MGFactor(config->GetnMGLevels()+1,1.0);
  const su2double CFLFactorDecrease = config->GetCFL_AdaptParam(0);
  const su2double CFLFactorIncrease = config->GetCFL_AdaptParam(1);
  const su2double CFLMin            = config->GetCFL_AdaptParam(2);
  const su2double CFLMax            = config->GetCFL_AdaptParam(3);
  const su2double acceptableLinTol  = config->GetCFL_AdaptParam(4);
  const bool fullComms              = (config->GetComm_Level() == COMM_FULL);

  /* Number of iterations considered to check for stagnation. */
  const auto Res_Count = min(100ul, config->GetnInner_Iter()-1);

  static bool reduceCFL, resetCFL, canIncrease;

  for (unsigned short iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {

    /* Store the mean flow, and turbulence solvers more clearly. */

    CSolver *solverFlow = solver_container[iMesh][FLOW_SOL];
    CSolver *solverTurb = solver_container[iMesh][TURB_SOL];

    /* Compute the reduction factor for CFLs on the coarse levels. */

    if (iMesh == MESH_0) {
      MGFactor[iMesh] = 1.0;
    } else {
      const su2double CFLRatio = config->GetCFL(iMesh)/config->GetCFL(iMesh-1);
      MGFactor[iMesh] = MGFactor[iMesh-1]*CFLRatio;
    }

    /* Check whether we achieved the requested reduction in the linear
     solver residual within the specified number of linear iterations. */

    su2double linResTurb = 0.0;
    if ((iMesh == MESH_0) && solverTurb) linResTurb = solverTurb->GetResLinSolver();

    /* Max linear residual between flow and turbulence. */
    const su2double linRes = max(solverFlow->GetResLinSolver(), linResTurb);

    /* Tolerance limited to an acceptable value. */
    const su2double linTol = max(acceptableLinTol, config->GetLinear_Solver_Error());

    /* Check that we are meeting our nonlinear residual reduction target
     over time so that we do not get stuck in limit cycles, this is done
     on the fine grid and applied to all others. */

    SU2_OMP_MASTER
    { /* Only the master thread updates the shared variables. */

    /* Check if we should decrease or if we can increase, the 20% is to avoid flip-flopping. */
    resetCFL = linRes > 0.99;
    reduceCFL = linRes > 1.2*linTol;
    canIncrease = linRes < linTol;

    if ((iMesh == MESH_0) && (Res_Count > 0)) {
      Old_Func = New_Func;
      if (NonLinRes_Series.empty()) NonLinRes_Series.resize(Res_Count,0.0);

      /* Sum the RMS residuals for all equations. */

      New_Func = 0.0;
      for (unsigned short iVar = 0; iVar < solverFlow->GetnVar(); iVar++) {
        New_Func += log10(solverFlow->GetRes_RMS(iVar));
      }
      if ((iMesh == MESH_0) && solverTurb) {
        for (unsigned short iVar = 0; iVar < solverTurb->GetnVar(); iVar++) {
          New_Func += log10(solverTurb->GetRes_RMS(iVar));
        }
      }

      /* Compute the difference in the nonlinear residuals between the
       current and previous iterations, taking care with very low initial
       residuals (due to initialization). */

      if ((config->GetInnerIter() == 1) && (New_Func - Old_Func > 10)) {
        Old_Func = New_Func;
      }
      NonLinRes_Series[NonLinRes_Counter] = New_Func - Old_Func;

      /* Increment the counter, if we hit the max size, then start over. */

      NonLinRes_Counter++;
      if (NonLinRes_Counter == Res_Count) NonLinRes_Counter = 0;

      /* Detect flip-flop convergence to reduce CFL and large increases
       to reset to minimum value, in that case clear the history. */

      if (config->GetInnerIter() >= Res_Count) {
        unsigned long signChanges = 0;
        su2double totalChange = 0.0;
        auto prev = NonLinRes_Series.front();
        for (auto val : NonLinRes_Series) {
          totalChange += val;
          signChanges += (prev > 0) ^ (val > 0);
          prev = val;
        }
        reduceCFL |= (signChanges > Res_Count/4) && (totalChange > -0.5);

        if (totalChange > 2.0) { // orders of magnitude
          resetCFL = true;
          NonLinRes_Counter = 0;
          for (auto& val : NonLinRes_Series) val = 0.0;
        }
      }
    }
    } /* End SU2_OMP_MASTER, now all threads update the CFL number. */
    END_SU2_OMP_MASTER
    SU2_OMP_BARRIER

    /* Loop over all points on this grid and apply CFL adaption. */

    su2double myCFLMin = 1e30, myCFLMax = 0.0, myCFLSum = 0.0;

    SU2_OMP_MASTER
    if ((iMesh == MESH_0) && fullComms) {
      Min_CFL_Local = 1e30;
      Max_CFL_Local = 0.0;
      Avg_CFL_Local = 0.0;
    }
    END_SU2_OMP_MASTER

    SU2_OMP_FOR_STAT(roundUpDiv(geometry[iMesh]->GetnPointDomain(),omp_get_max_threads()))
    for (unsigned long iPoint = 0; iPoint < geometry[iMesh]->GetnPointDomain(); iPoint++) {

      /* Get the current local flow CFL number at this point. */

      su2double CFL = solverFlow->GetNodes()->GetLocalCFL(iPoint);

      /* Get the current under-relaxation parameters that were computed
       during the previous nonlinear update. If we have a turbulence model,
       take the minimum under-relaxation parameter between the mean flow
       and turbulence systems. */

      su2double underRelaxationFlow = solverFlow->GetNodes()->GetUnderRelaxation(iPoint);
      su2double underRelaxationTurb = 1.0;
      if ((iMesh == MESH_0) && solverTurb)
        underRelaxationTurb = solverTurb->GetNodes()->GetUnderRelaxation(iPoint);
      const su2double underRelaxation = min(underRelaxationFlow,underRelaxationTurb);

      /* If we apply a small under-relaxation parameter for stability,
       then we should reduce the CFL before the next iteration. If we
       are able to add the entire nonlinear update (under-relaxation = 1)
       then we schedule an increase the CFL number for the next iteration. */

      su2double CFLFactor = 1.0;
      if (underRelaxation < 0.1 || reduceCFL) {
        CFLFactor = CFLFactorDecrease;
      } else if ((underRelaxation >= 0.1 && underRelaxation < 1.0) || !canIncrease) {
        CFLFactor = 1.0;
      } else {
        CFLFactor = CFLFactorIncrease;
      }

      /* Check if we are hitting the min or max and adjust. */

      if (CFL*CFLFactor <= CFLMin) {
        CFL       = CFLMin;
        CFLFactor = MGFactor[iMesh];
      } else if (CFL*CFLFactor >= CFLMax) {
        CFL       = CFLMax;
        CFLFactor = MGFactor[iMesh];
      }

      /* If we detect a stalled nonlinear residual, then force the CFL
       for all points to the minimum temporarily to restart the ramp. */

      if (resetCFL) {
        CFL       = CFLMin;
        CFLFactor = MGFactor[iMesh];
      }

      /* Apply the adjustment to the CFL and store local values. */

      CFL *= CFLFactor;
      solverFlow->GetNodes()->SetLocalCFL(iPoint, CFL);
      if ((iMesh == MESH_0) && solverTurb) {
        solverTurb->GetNodes()->SetLocalCFL(iPoint, CFL);
      }

      /* Store min and max CFL for reporting on the fine grid. */

      if ((iMesh == MESH_0) && fullComms) {
        myCFLMin = min(CFL,myCFLMin);
        myCFLMax = max(CFL,myCFLMax);
        myCFLSum += CFL;
      }

    }
    END_SU2_OMP_FOR

    /* Reduce the min/max/avg local CFL numbers. */

    if ((iMesh == MESH_0) && fullComms) {
      SU2_OMP_CRITICAL
      { /* OpenMP reduction. */
        Min_CFL_Local = min(Min_CFL_Local,myCFLMin);
        Max_CFL_Local = max(Max_CFL_Local,myCFLMax);
        Avg_CFL_Local += myCFLSum;
      }
      END_SU2_OMP_CRITICAL
      SU2_OMP_BARRIER

      SU2_OMP_MASTER
      { /* MPI reduction. */
        myCFLMin = Min_CFL_Local; myCFLMax = Max_CFL_Local; myCFLSum = Avg_CFL_Local;
        SU2_MPI::Allreduce(&myCFLMin, &Min_CFL_Local, 1, MPI_DOUBLE, MPI_MIN, SU2_MPI::GetComm());
        SU2_MPI::Allreduce(&myCFLMax, &Max_CFL_Local, 1, MPI_DOUBLE, MPI_MAX, SU2_MPI::GetComm());
        SU2_MPI::Allreduce(&myCFLSum, &Avg_CFL_Local, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
        Avg_CFL_Local /= su2double(geometry[iMesh]->GetGlobal_nPointDomain());
      }
      END_SU2_OMP_MASTER
      SU2_OMP_BARRIER
    }

  }

}

void CSolver::SetResidual_RMS(const CGeometry *geometry, const CConfig *config) {

  if (geometry->GetMGLevel() != MESH_0) return;

  SU2_OMP_MASTER {

  /*--- Set the L2 Norm residual in all the processors. ---*/

  vector<su2double> rbuf_res(nVar);
  unsigned long Global_nPointDomain = 0;

  if (config->GetComm_Level() == COMM_FULL) {

    SU2_MPI::Allreduce(Residual_RMS.data(), rbuf_res.data(), nVar, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
    Global_nPointDomain = geometry->GetGlobal_nPointDomain();
  }
  else {
    /*--- Reduced MPI comms have been requested. Use a local residual only. ---*/

    for (unsigned short iVar = 0; iVar < nVar; iVar++) rbuf_res[iVar] = Residual_RMS[iVar];
    Global_nPointDomain = geometry->GetnPointDomain();
  }

  for (unsigned short iVar = 0; iVar < nVar; iVar++) {

    if (std::isnan(SU2_TYPE::GetValue(rbuf_res[iVar]))) {
      SU2_MPI::Error("SU2 has diverged (NaN detected).", CURRENT_FUNCTION);
    }

    Residual_RMS[iVar] = max(EPS*EPS, sqrt(rbuf_res[iVar]/Global_nPointDomain));

    if (log10(GetRes_RMS(iVar)) > 20.0) {
      SU2_MPI::Error("SU2 has diverged (Residual > 10^20 detected).", CURRENT_FUNCTION);
    }
  }

  /*--- Set the Maximum residual in all the processors. ---*/

  if (config->GetComm_Level() == COMM_FULL) {

    const unsigned long nProcessor = size;

    su2activematrix rbuf_residual(nProcessor,nVar);
    su2matrix<unsigned long> rbuf_point(nProcessor,nVar);
    su2activematrix rbuf_coord(nProcessor*nVar, nDim);

    SU2_MPI::Allgather(Residual_Max.data(), nVar, MPI_DOUBLE, rbuf_residual.data(), nVar, MPI_DOUBLE, SU2_MPI::GetComm());
    SU2_MPI::Allgather(Point_Max.data(), nVar, MPI_UNSIGNED_LONG, rbuf_point.data(), nVar, MPI_UNSIGNED_LONG, SU2_MPI::GetComm());
    SU2_MPI::Allgather(Point_Max_Coord.data(), nVar*nDim, MPI_DOUBLE, rbuf_coord.data(), nVar*nDim, MPI_DOUBLE, SU2_MPI::GetComm());

    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
      for (auto iProcessor = 0ul; iProcessor < nProcessor; iProcessor++) {
        AddRes_Max(iVar, rbuf_residual(iProcessor,iVar), rbuf_point(iProcessor,iVar), rbuf_coord[iProcessor*nVar+iVar]);
      }
    }
  }

  }
  END_SU2_OMP_MASTER
  SU2_OMP_BARRIER
}

void CSolver::SetResidual_BGS(const CGeometry *geometry, const CConfig *config) {

  if (geometry->GetMGLevel() != MESH_0) return;

  SU2_OMP_MASTER {

  /*--- Set the L2 Norm residual in all the processors. ---*/

  vector<su2double> rbuf_res(nVar);

  SU2_MPI::Allreduce(Residual_BGS.data(), rbuf_res.data(), nVar, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
  const auto Global_nPointDomain = geometry->GetGlobal_nPointDomain();

  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    Residual_BGS[iVar] = max(EPS*EPS, sqrt(rbuf_res[iVar]/Global_nPointDomain));
  }

  if (config->GetComm_Level() == COMM_FULL) {

    /*--- Set the Maximum residual in all the processors. ---*/

    const unsigned long nProcessor = size;

    su2activematrix rbuf_residual(nProcessor,nVar);
    su2matrix<unsigned long> rbuf_point(nProcessor,nVar);
    su2activematrix rbuf_coord(nProcessor*nVar, nDim);

    SU2_MPI::Allgather(Residual_Max_BGS.data(), nVar, MPI_DOUBLE, rbuf_residual.data(), nVar, MPI_DOUBLE, SU2_MPI::GetComm());
    SU2_MPI::Allgather(Point_Max_BGS.data(), nVar, MPI_UNSIGNED_LONG, rbuf_point.data(), nVar, MPI_UNSIGNED_LONG, SU2_MPI::GetComm());
    SU2_MPI::Allgather(Point_Max_Coord_BGS.data(), nVar*nDim, MPI_DOUBLE, rbuf_coord.data(), nVar*nDim, MPI_DOUBLE, SU2_MPI::GetComm());

    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
      for (auto iProcessor = 0ul; iProcessor < nProcessor; iProcessor++) {
        AddRes_Max_BGS(iVar, rbuf_residual(iProcessor,iVar), rbuf_point(iProcessor,iVar), rbuf_coord[iProcessor*nVar+iVar]);
      }
    }
  }

  }
  END_SU2_OMP_MASTER
  SU2_OMP_BARRIER
}

void CSolver::SetRotatingFrame_GCL(CGeometry *geometry, const CConfig *config) {

  /*--- Loop interior points ---*/

  SU2_OMP_FOR_STAT(roundUpDiv(nPointDomain,2*omp_get_max_threads()))
  for (auto iPoint = 0ul; iPoint < nPointDomain; ++iPoint) {

    const su2double* GridVel_i = geometry->nodes->GetGridVel(iPoint);
    const su2double* Solution_i = base_nodes->GetSolution(iPoint);

    for (auto iNeigh = 0u; iNeigh < geometry->nodes->GetnPoint(iPoint); iNeigh++) {

      const auto iEdge = geometry->nodes->GetEdge(iPoint, iNeigh);
      const su2double* Normal = geometry->edges->GetNormal(iEdge);

      const auto jPoint = geometry->nodes->GetPoint(iPoint, iNeigh);
      const su2double* GridVel_j = geometry->nodes->GetGridVel(jPoint);

      /*--- Determine whether to consider the normal outward or inward. ---*/
      su2double dir = (iPoint < jPoint)? 0.5 : -0.5;

      su2double Flux = 0.0;
      for (auto iDim = 0u; iDim < nDim; iDim++)
        Flux += dir*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];

      for (auto iVar = 0u; iVar < nVar; iVar++)
        LinSysRes(iPoint,iVar) += Flux * Solution_i[iVar];
    }
  }
  END_SU2_OMP_FOR

  /*--- Loop boundary edges ---*/

  for (auto iMarker = 0u; iMarker < geometry->GetnMarker(); iMarker++) {
    if ((config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY) &&
        (config->GetMarker_All_KindBC(iMarker) != NEARFIELD_BOUNDARY) &&
        (config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY)) {

      SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
      for (auto iVertex = 0u; iVertex < geometry->GetnVertex(iMarker); iVertex++) {

        const auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

        /*--- Grid Velocity at each edge point ---*/

        const su2double* GridVel = geometry->nodes->GetGridVel(iPoint);

        /*--- Summed normal components ---*/

        const su2double* Normal = geometry->vertex[iMarker][iVertex]->GetNormal();

        su2double Flux = GeometryToolbox::DotProduct(nDim, Normal, GridVel);

        for (auto iVar = 0u; iVar < nVar; iVar++)
          LinSysRes(iPoint,iVar) -= Flux * base_nodes->GetSolution(iPoint,iVar);
      }
      END_SU2_OMP_FOR
    }
  }

}

void CSolver::SetAuxVar_Gradient_GG(CGeometry *geometry, const CConfig *config) {

  const auto& solution = base_nodes->GetAuxVar();
  auto& gradient = base_nodes->GetAuxVarGradient();

  computeGradientsGreenGauss(this, AUXVAR_GRADIENT, PERIODIC_NONE, *geometry,
                             *config, solution, 0, base_nodes->GetnAuxVar(), gradient);
}

void CSolver::SetAuxVar_Gradient_LS(CGeometry *geometry, const CConfig *config) {

  bool weighted = true;
  const auto& solution = base_nodes->GetAuxVar();
  auto& gradient = base_nodes->GetAuxVarGradient();
  auto& rmatrix  = base_nodes->GetRmatrix();

  computeGradientsLeastSquares(this, AUXVAR_GRADIENT, PERIODIC_NONE, *geometry, *config,
                               weighted, solution, 0, base_nodes->GetnAuxVar(), gradient, rmatrix);
}

void CSolver::SetSolution_Gradient_GG(CGeometry *geometry, const CConfig *config, bool reconstruction) {

  const auto& solution = base_nodes->GetSolution();
  auto& gradient = reconstruction? base_nodes->GetGradient_Reconstruction() : base_nodes->GetGradient();
  const auto comm = reconstruction? SOLUTION_GRAD_REC : SOLUTION_GRADIENT;
  const auto commPer = reconstruction? PERIODIC_SOL_GG_R : PERIODIC_SOL_GG;

  computeGradientsGreenGauss(this, comm, commPer, *geometry, *config, solution, 0, nVar, gradient);
}

void CSolver::SetSolution_Gradient_LS(CGeometry *geometry, const CConfig *config, bool reconstruction) {

  /*--- Set a flag for unweighted or weighted least-squares. ---*/
  bool weighted;
  PERIODIC_QUANTITIES commPer;

  if (reconstruction) {
    weighted = (config->GetKind_Gradient_Method_Recon() == WEIGHTED_LEAST_SQUARES);
    commPer = weighted? PERIODIC_SOL_LS_R : PERIODIC_SOL_ULS_R;
  }
  else {
    weighted = (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES);
    commPer = weighted? PERIODIC_SOL_LS : PERIODIC_SOL_ULS;
  }

  const auto& solution = base_nodes->GetSolution();
  auto& rmatrix = base_nodes->GetRmatrix();
  auto& gradient = reconstruction? base_nodes->GetGradient_Reconstruction() : base_nodes->GetGradient();
  const auto comm = reconstruction? SOLUTION_GRAD_REC : SOLUTION_GRADIENT;

  computeGradientsLeastSquares(this, comm, commPer, *geometry, *config, weighted, solution, 0, nVar, gradient, rmatrix);
}

void CSolver::SetUndivided_Laplacian(CGeometry *geometry, const CConfig *config) {

  /*--- Loop domain points. ---*/

  SU2_OMP_FOR_DYN(256)
  for (unsigned long iPoint = 0; iPoint < nPointDomain; ++iPoint) {

    const bool boundary_i = geometry->nodes->GetPhysicalBoundary(iPoint);

    /*--- Initialize. ---*/
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
      base_nodes->SetUnd_Lapl(iPoint, iVar, 0.0);

    /*--- Loop over the neighbors of point i. ---*/
    for (auto jPoint : geometry->nodes->GetPoints(iPoint)) {

      bool boundary_j = geometry->nodes->GetPhysicalBoundary(jPoint);

      /*--- If iPoint is boundary it only takes contributions from other boundary points. ---*/
      if (boundary_i && !boundary_j) continue;

      /*--- Add solution differences, with correction for compressible flows which use the enthalpy. ---*/

      for (unsigned short iVar = 0; iVar < nVar; iVar++) {
        su2double delta = base_nodes->GetSolution(jPoint,iVar)-base_nodes->GetSolution(iPoint,iVar);
        base_nodes->AddUnd_Lapl(iPoint, iVar, delta);
      }
    }
  }
  END_SU2_OMP_FOR

  /*--- Correct the Laplacian across any periodic boundaries. ---*/

  for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
    InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_LAPLACIAN);
    CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_LAPLACIAN);
  }

  /*--- MPI parallelization ---*/

  InitiateComms(geometry, config, UNDIVIDED_LAPLACIAN);
  CompleteComms(geometry, config, UNDIVIDED_LAPLACIAN);

}

void CSolver::Add_External_To_Solution() {
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {
    base_nodes->AddSolution(iPoint, base_nodes->Get_External(iPoint));
  }

  base_nodes->Add_ExternalExtra_To_SolutionExtra();
}

void CSolver::Add_Solution_To_External() {
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {
    base_nodes->Add_External(iPoint, base_nodes->GetSolution(iPoint));
  }

  base_nodes->Set_ExternalExtra_To_SolutionExtra();
}

void CSolver::Update_Cross_Term(CConfig *config, su2passivematrix &cross_term) {

  /*--- This method is for discrete adjoint solvers and it is used in multi-physics
   *    contexts, "cross_term" is the old value, the new one is in "Solution".
   *    We update "cross_term" and the sum of all cross terms (in "External")
   *    with a fraction of the difference between new and old.
   *    When "alpha" is 1, i.e. no relaxation, we effectively subtract the old
   *    value and add the new one to the total ("External"). ---*/

  vector<su2double> solution(nVar);
  passivedouble alpha = SU2_TYPE::GetValue(config->GetAitkenStatRelax());

  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
      passivedouble
      new_val = SU2_TYPE::GetValue(base_nodes->GetSolution(iPoint,iVar)),
      delta = alpha * (new_val - cross_term(iPoint,iVar));
      /*--- Update cross term. ---*/
      cross_term(iPoint,iVar) += delta;
      solution[iVar] = delta;
    }
    /*--- Update the sum of all cross-terms. ---*/
    base_nodes->Add_External(iPoint, solution.data());
  }
}

void CSolver::SetGridVel_Gradient(CGeometry *geometry, const CConfig *config) {

  /// TODO: No comms needed for this gradient? The Rmatrix should be allocated somewhere.

  const auto& gridVel = geometry->nodes->GetGridVel();
  auto& gridVelGrad = geometry->nodes->GetGridVel_Grad();
  auto rmatrix = CVectorOfMatrix(nPoint,nDim,nDim);

  computeGradientsLeastSquares(nullptr, GRID_VELOCITY, PERIODIC_NONE, *geometry, *config,
                               true, gridVel, 0, nDim, gridVelGrad, rmatrix);
}

void CSolver::SetSolution_Limiter(CGeometry *geometry, const CConfig *config) {

  const auto kindLimiter = config->GetKind_SlopeLimit();
  const auto& solution = base_nodes->GetSolution();
  const auto& gradient = base_nodes->GetGradient_Reconstruction();
  auto& solMin = base_nodes->GetSolution_Min();
  auto& solMax = base_nodes->GetSolution_Max();
  auto& limiter = base_nodes->GetLimiter();

  computeLimiters(kindLimiter, this, SOLUTION_LIMITER, PERIODIC_LIM_SOL_1, PERIODIC_LIM_SOL_2,
                  *geometry, *config, 0, nVar, solution, gradient, solMin, solMax, limiter);
}

void CSolver::Gauss_Elimination(su2double** A, su2double* rhs, unsigned short nVar) {

  short iVar, jVar, kVar;
  su2double weight, aux;

  if (nVar == 1)
    rhs[0] /= A[0][0];
  else {

    /*--- Transform system in Upper Matrix ---*/

    for (iVar = 1; iVar < (short)nVar; iVar++) {
      for (jVar = 0; jVar < iVar; jVar++) {
        weight = A[iVar][jVar]/A[jVar][jVar];
        for (kVar = jVar; kVar < (short)nVar; kVar++)
          A[iVar][kVar] -= weight*A[jVar][kVar];
        rhs[iVar] -= weight*rhs[jVar];
      }
    }

    /*--- Backwards substitution ---*/

    rhs[nVar-1] = rhs[nVar-1]/A[nVar-1][nVar-1];
    for (iVar = (short)nVar-2; iVar >= 0; iVar--) {
      aux = 0;
      for (jVar = iVar+1; jVar < (short)nVar; jVar++)
        aux += A[iVar][jVar]*rhs[jVar];
      rhs[iVar] = (rhs[iVar]-aux)/A[iVar][iVar];
      if (iVar == 0) break;
    }
  }

}

void CSolver::Restart_OldGeometry(CGeometry *geometry, CConfig *config) {

  SU2_OMP_MASTER {

  /*--- This function is intended for dual time simulations ---*/

  int Unst_RestartIter;
  ifstream restart_file_n;

  string filename = config->GetSolution_FileName();
  string filename_n;

  /*--- Auxiliary vector for storing the coordinates ---*/
  su2double Coord[3] = {0.0};

  /*--- Variables for reading the restart files ---*/
  string text_line;
  long iPoint_Local;
  unsigned long iPoint_Global_Local = 0, iPoint_Global = 0;

  /*--- First, we load the restart file for time n ---*/

  /*-------------------------------------------------------------------------------------------*/

  /*--- Modify file name for an unsteady restart ---*/
  if (config->GetRestart()) Unst_RestartIter = SU2_TYPE::Int(config->GetRestart_Iter())-1;
  else Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_AdjointIter())-1;
  filename_n = config->GetFilename(filename, ".csv", Unst_RestartIter);

  /*--- Open the restart file, throw an error if this fails. ---*/

  restart_file_n.open(filename_n.data(), ios::in);
  if (restart_file_n.fail()) {
    SU2_MPI::Error(string("There is no flow restart file ") + filename_n, CURRENT_FUNCTION);
  }

  /*--- First, set all indices to a negative value by default, and Global n indices to 0 ---*/
  iPoint_Global_Local = 0; iPoint_Global = 0;

  /*--- Read all lines in the restart file ---*/
  /*--- The first line is the header ---*/

  getline (restart_file_n, text_line);

  for (iPoint_Global = 0; iPoint_Global < geometry->GetGlobal_nPointDomain(); iPoint_Global++ ) {

    getline (restart_file_n, text_line);

    vector<string> point_line = PrintingToolbox::split(text_line, ',');

    /*--- Retrieve local index. If this node from the restart file lives
     on the current processor, we will load and instantiate the vars. ---*/

    iPoint_Local = geometry->GetGlobal_to_Local_Point(iPoint_Global);

    if (iPoint_Local > -1) {

      Coord[0] = PrintingToolbox::stod(point_line[1]);
      Coord[1] = PrintingToolbox::stod(point_line[2]);
      if (nDim == 3){
        Coord[2] = PrintingToolbox::stod(point_line[3]);
      }
      geometry->nodes->SetCoord_n(iPoint_Local, Coord);

      iPoint_Global_Local++;
    }
  }

  /*--- Detect a wrong solution file ---*/

  if (iPoint_Global_Local < geometry->GetnPointDomain()) {
    SU2_MPI::Error(string("The solution file ") + filename + string(" doesn't match with the mesh file!\n") +
                   string("It could be empty lines at the end of the file."), CURRENT_FUNCTION);
  }

  /*--- Close the restart file ---*/

  restart_file_n.close();

  /*-------------------------------------------------------------------------------------------*/
  /*-------------------------------------------------------------------------------------------*/

  /*--- Now, we load the restart file for time n-1, if the simulation is 2nd Order ---*/

  if (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND) {

    ifstream restart_file_n1;
    string filename_n1;

    /*--- Modify file name for an unsteady restart ---*/
    if (config->GetRestart()) Unst_RestartIter = SU2_TYPE::Int(config->GetRestart_Iter())-2;
    else Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_AdjointIter())-2;
    filename_n1 = config->GetFilename(filename, ".csv", Unst_RestartIter);

    /*--- Open the restart file, throw an error if this fails. ---*/

    restart_file_n1.open(filename_n1.data(), ios::in);
    if (restart_file_n1.fail()) {
        SU2_MPI::Error(string("There is no flow restart file ") + filename_n1, CURRENT_FUNCTION);

    }

    /*--- First, set all indices to a negative value by default, and Global n indices to 0 ---*/
    iPoint_Global_Local = 0; iPoint_Global = 0;

    /*--- Read all lines in the restart file ---*/
    /*--- The first line is the header ---*/

    getline (restart_file_n1, text_line);

    for (iPoint_Global = 0; iPoint_Global < geometry->GetGlobal_nPointDomain(); iPoint_Global++ ) {

      getline (restart_file_n1, text_line);

      vector<string> point_line = PrintingToolbox::split(text_line, ',');

      /*--- Retrieve local index. If this node from the restart file lives
       on the current processor, we will load and instantiate the vars. ---*/

      iPoint_Local = geometry->GetGlobal_to_Local_Point(iPoint_Global);

      if (iPoint_Local > -1) {

        Coord[0] = PrintingToolbox::stod(point_line[1]);
        Coord[1] = PrintingToolbox::stod(point_line[2]);
        if (nDim == 3){
          Coord[2] = PrintingToolbox::stod(point_line[3]);
        }

        geometry->nodes->SetCoord_n1(iPoint_Local, Coord);

        iPoint_Global_Local++;
      }

    }

    /*--- Detect a wrong solution file ---*/

    if (iPoint_Global_Local < geometry->GetnPointDomain()) {
      SU2_MPI::Error(string("The solution file ") + filename + string(" doesn't match with the mesh file!\n") +
                     string("It could be empty lines at the end of the file."), CURRENT_FUNCTION);
    }

    /*--- Close the restart file ---*/

    restart_file_n1.close();

  }

  }
  END_SU2_OMP_MASTER
  SU2_OMP_BARRIER

  /*--- It's necessary to communicate this information ---*/

  geometry->InitiateComms(geometry, config, COORDINATES_OLD);
  geometry->CompleteComms(geometry, config, COORDINATES_OLD);

}

void CSolver::Read_SU2_Restart_ASCII(CGeometry *geometry, const CConfig *config, string val_filename) {

  ifstream restart_file;
  string text_line, Tag;
  unsigned short iVar;
  long iPoint_Local = 0; unsigned long iPoint_Global = 0;
  int counter = 0;
  fields.clear();

  Restart_Vars = new int[5];

  string error_string = "Note: ASCII restart files must be in CSV format since v7.0.\n"
                        "Check https://su2code.github.io/docs/Guide-to-v7 for more information.";

  /*--- First, check that this is not a binary restart file. ---*/

  char fname[100];
  val_filename += ".csv";
  strcpy(fname, val_filename.c_str());
  int magic_number;

#ifndef HAVE_MPI

  /*--- Serial binary input. ---*/

  FILE *fhw;
  fhw = fopen(fname,"rb");
  size_t ret;

  /*--- Error check for opening the file. ---*/

  if (!fhw) {
    SU2_MPI::Error(string("Unable to open SU2 restart file ") + fname, CURRENT_FUNCTION);
  }

  /*--- Attempt to read the first int, which should be our magic number. ---*/

  ret = fread(&magic_number, sizeof(int), 1, fhw);
  if (ret != 1) {
    SU2_MPI::Error("Error reading restart file.", CURRENT_FUNCTION);
  }

  /*--- Check that this is an SU2 binary file. SU2 binary files
   have the hex representation of "SU2" as the first int in the file. ---*/

  if (magic_number == 535532) {
    SU2_MPI::Error(string("File ") + string(fname) + string(" is a binary SU2 restart file, expected ASCII.\n") +
                   string("SU2 reads/writes binary restart files by default.\n") +
                   string("Note that backward compatibility for ASCII restart files is\n") +
                   string("possible with the READ_BINARY_RESTART option."), CURRENT_FUNCTION);
  }

  fclose(fhw);

#else

  /*--- Parallel binary input using MPI I/O. ---*/

  MPI_File fhw;
  int ierr;

  /*--- All ranks open the file using MPI. ---*/

  ierr = MPI_File_open(SU2_MPI::GetComm(), fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fhw);

  /*--- Error check opening the file. ---*/

  if (ierr) {
    SU2_MPI::Error(string("SU2 ASCII restart file ") + string(fname) + string(" not found.\n") + error_string,
                   CURRENT_FUNCTION);
  }

  /*--- Have the master attempt to read the magic number. ---*/

  if (rank == MASTER_NODE)
    MPI_File_read(fhw, &magic_number, 1, MPI_INT, MPI_STATUS_IGNORE);

  /*--- Broadcast the number of variables to all procs and store clearly. ---*/

  SU2_MPI::Bcast(&magic_number, 1, MPI_INT, MASTER_NODE, SU2_MPI::GetComm());

  /*--- Check that this is an SU2 binary file. SU2 binary files
   have the hex representation of "SU2" as the first int in the file. ---*/

  if (magic_number == 535532) {
    SU2_MPI::Error(string("File ") + string(fname) + string(" is a binary SU2 restart file, expected ASCII.\n") +
                   string("SU2 reads/writes binary restart files by default.\n") +
                   string("Note that backward compatibility for ASCII restart files is\n") +
                   string("possible with the READ_BINARY_RESTART option."), CURRENT_FUNCTION);
  }

  MPI_File_close(&fhw);

#endif

  /*--- Open the restart file ---*/

  restart_file.open(val_filename.data(), ios::in);

  /*--- In case there is no restart file ---*/

  if (restart_file.fail()) {
    SU2_MPI::Error(string("SU2 ASCII restart file ") + string(fname) + string(" not found.\n") + error_string,
                   CURRENT_FUNCTION);
  }

  /*--- Identify the number of fields (and names) in the restart file ---*/

  getline (restart_file, text_line);

  char delimiter = ',';
  fields = PrintingToolbox::split(text_line, delimiter);

  if (fields.size() <= 1) {
    SU2_MPI::Error(string("Restart file does not seem to be a CSV file.\n") + error_string, CURRENT_FUNCTION);
  }

  for (unsigned short iField = 0; iField < fields.size(); iField++){
    PrintingToolbox::trim(fields[iField]);
  }

  /*--- Set the number of variables, one per field in the
   restart file (without including the PointID) ---*/

  Restart_Vars[1] = (int)fields.size() - 1;

  /*--- Allocate memory for the restart data. ---*/

  Restart_Data = new passivedouble[Restart_Vars[1]*geometry->GetnPointDomain()];

  /*--- Read all lines in the restart file and extract data. ---*/

  for (iPoint_Global = 0; iPoint_Global < geometry->GetGlobal_nPointDomain(); iPoint_Global++) {

    if (!getline (restart_file, text_line)) break;

    /*--- Retrieve local index. If this node from the restart file lives
     on the current processor, we will load and instantiate the vars. ---*/

    iPoint_Local = geometry->GetGlobal_to_Local_Point(iPoint_Global);

    if (iPoint_Local > -1) {

      vector<string> point_line = PrintingToolbox::split(text_line, delimiter);

      /*--- Store the solution (starting with node coordinates) --*/

      for (iVar = 0; iVar < Restart_Vars[1]; iVar++)
        Restart_Data[counter*Restart_Vars[1] + iVar] = SU2_TYPE::GetValue(PrintingToolbox::stod(point_line[iVar+1]));

      /*--- Increment our local point counter. ---*/

      counter++;

    }
  }

  if (iPoint_Global != geometry->GetGlobal_nPointDomain())
    SU2_MPI::Error("The solution file does not match the mesh, currently only binary files can be interpolated.",
                   CURRENT_FUNCTION);

}

void CSolver::Read_SU2_Restart_Binary(CGeometry *geometry, const CConfig *config, string val_filename) {

  char str_buf[CGNS_STRING_SIZE], fname[100];
  val_filename += ".dat";
  strcpy(fname, val_filename.c_str());
  const int nRestart_Vars = 5;
  Restart_Vars = new int[nRestart_Vars];
  fields.clear();

#ifndef HAVE_MPI

  /*--- Serial binary input. ---*/

  FILE *fhw;
  fhw = fopen(fname,"rb");
  size_t ret;

  /*--- Error check for opening the file. ---*/

  if (!fhw) {
    SU2_MPI::Error(string("Unable to open SU2 restart file ") + string(fname), CURRENT_FUNCTION);
  }

  /*--- First, read the number of variables and points. ---*/

  ret = fread(Restart_Vars, sizeof(int), nRestart_Vars, fhw);
  if (ret != (unsigned long)nRestart_Vars) {
    SU2_MPI::Error("Error reading restart file.", CURRENT_FUNCTION);
  }

  /*--- Check that this is an SU2 binary file. SU2 binary files
   have the hex representation of "SU2" as the first int in the file. ---*/

  if (Restart_Vars[0] != 535532) {
    SU2_MPI::Error(string("File ") + string(fname) + string(" is not a binary SU2 restart file.\n") +
                   string("SU2 reads/writes binary restart files by default.\n") +
                   string("Note that backward compatibility for ASCII restart files is\n") +
                   string("possible with the READ_BINARY_RESTART option."), CURRENT_FUNCTION);
  }

  /*--- Store the number of fields and points to be read for clarity. ---*/

  const unsigned long nFields = Restart_Vars[1];
  const unsigned long nPointFile = Restart_Vars[2];

  /*--- Read the variable names from the file. Note that we are adopting a
   fixed length of 33 for the string length to match with CGNS. This is
   needed for when we read the strings later. We pad the beginning of the
   variable string vector with the Point_ID tag that wasn't written. ---*/

  fields.push_back("Point_ID");
  for (auto iVar = 0u; iVar < nFields; iVar++) {
    ret = fread(str_buf, sizeof(char), CGNS_STRING_SIZE, fhw);
    if (ret != (unsigned long)CGNS_STRING_SIZE) {
      SU2_MPI::Error("Error reading restart file.", CURRENT_FUNCTION);
    }
    fields.push_back(str_buf);
  }

  /*--- For now, create a temp 1D buffer to read the data from file. ---*/

  Restart_Data = new passivedouble[nFields*nPointFile];

  /*--- Read in the data for the restart at all local points. ---*/

  ret = fread(Restart_Data, sizeof(passivedouble), nFields*nPointFile, fhw);
  if (ret != nFields*nPointFile) {
    SU2_MPI::Error("Error reading restart file.", CURRENT_FUNCTION);
  }

  /*--- Close the file. ---*/

  fclose(fhw);

#else

  /*--- Parallel binary input using MPI I/O. ---*/

  MPI_File fhw;
  SU2_MPI::Status status;
  MPI_Datatype etype, filetype;
  MPI_Offset disp;

  /*--- All ranks open the file using MPI. ---*/

  int ierr = MPI_File_open(SU2_MPI::GetComm(), fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fhw);

  if (ierr) SU2_MPI::Error(string("Unable to open SU2 restart file ") + string(fname), CURRENT_FUNCTION);

  /*--- First, read the number of variables and points (i.e., cols and rows),
   which we will need in order to read the file later. Also, read the
   variable string names here. Only the master rank reads the header. ---*/

  if (rank == MASTER_NODE)
    MPI_File_read(fhw, Restart_Vars, nRestart_Vars, MPI_INT, MPI_STATUS_IGNORE);

  /*--- Broadcast the number of variables to all procs and store clearly. ---*/

  SU2_MPI::Bcast(Restart_Vars, nRestart_Vars, MPI_INT, MASTER_NODE, SU2_MPI::GetComm());

  /*--- Check that this is an SU2 binary file. SU2 binary files
   have the hex representation of "SU2" as the first int in the file. ---*/

  if (Restart_Vars[0] != 535532) {
    SU2_MPI::Error(string("File ") + string(fname) + string(" is not a binary SU2 restart file.\n") +
                   string("SU2 reads/writes binary restart files by default.\n") +
                   string("Note that backward compatibility for ASCII restart files is\n") +
                   string("possible with the READ_BINARY_RESTART option."), CURRENT_FUNCTION);
  }

  /*--- Store the number of fields and points to be read for clarity. ---*/

  const unsigned long nFields = Restart_Vars[1];
  const unsigned long nPointFile = Restart_Vars[2];

  /*--- Read the variable names from the file. Note that we are adopting a
   fixed length of 33 for the string length to match with CGNS. This is
   needed for when we read the strings later. ---*/

  char *mpi_str_buf = new char[nFields*CGNS_STRING_SIZE];
  if (rank == MASTER_NODE) {
    disp = nRestart_Vars*sizeof(int);
    MPI_File_read_at(fhw, disp, mpi_str_buf, nFields*CGNS_STRING_SIZE,
                     MPI_CHAR, MPI_STATUS_IGNORE);
  }

  /*--- Broadcast the string names of the variables. ---*/

  SU2_MPI::Bcast(mpi_str_buf, nFields*CGNS_STRING_SIZE, MPI_CHAR,
                 MASTER_NODE, SU2_MPI::GetComm());

  /*--- Now parse the string names and load into the config class in case
   we need them for writing visualization files (SU2_SOL). ---*/

  fields.push_back("Point_ID");
  for (auto iVar = 0u; iVar < nFields; iVar++) {
    const auto index = iVar*CGNS_STRING_SIZE;
    string field_buf("\"");
    for (int iChar = 0; iChar < CGNS_STRING_SIZE; iChar++) {
      str_buf[iChar] = mpi_str_buf[index + iChar];
    }
    field_buf.append(str_buf);
    field_buf.append("\"");
    fields.push_back(field_buf.c_str());
  }

  /*--- Free string buffer memory. ---*/

  delete [] mpi_str_buf;

  /*--- We're writing only su2doubles in the data portion of the file. ---*/

  etype = MPI_DOUBLE;

  /*--- We need to ignore the 4 ints describing the nVar_Restart and nPoints,
   along with the string names of the variables. ---*/

  disp = nRestart_Vars*sizeof(int) + CGNS_STRING_SIZE*nFields*sizeof(char);

  /*--- Define a derived datatype for this rank's set of non-contiguous data
   that will be placed in the restart. Here, we are collecting each one of the
   points which are distributed throughout the file in blocks of nVar_Restart data. ---*/

  int nBlock;
  int *blocklen = nullptr;
  MPI_Aint *displace = nullptr;

  if (nPointFile == geometry->GetGlobal_nPointDomain() ||
      config->GetKind_SU2() == SU2_COMPONENT::SU2_SOL) {
    /*--- No interpolation, each rank reads the indices it needs. ---*/
    nBlock = geometry->GetnPointDomain();

    blocklen = new int[nBlock];
    displace = new MPI_Aint[nBlock];
    int counter = 0;
    for (auto iPoint_Global = 0ul; iPoint_Global < geometry->GetGlobal_nPointDomain(); ++iPoint_Global) {
      if (geometry->GetGlobal_to_Local_Point(iPoint_Global) > -1) {
        blocklen[counter] = nFields;
        displace[counter] = iPoint_Global*nFields*sizeof(passivedouble);
        counter++;
      }
    }
  }
  else {
    /*--- Interpolation required, read large blocks of data. ---*/
    nBlock = 1;

    blocklen = new int[nBlock];
    displace = new MPI_Aint[nBlock];

    const auto partitioner = CLinearPartitioner(nPointFile,0);

    blocklen[0] = nFields*partitioner.GetSizeOnRank(rank);
    displace[0] = nFields*partitioner.GetFirstIndexOnRank(rank)*sizeof(passivedouble);;
  }

  MPI_Type_create_hindexed(nBlock, blocklen, displace, MPI_DOUBLE, &filetype);
  MPI_Type_commit(&filetype);

  /*--- Set the view for the MPI file write, i.e., describe the location in
   the file that this rank "sees" for writing its piece of the restart file. ---*/

  MPI_File_set_view(fhw, disp, etype, filetype, (char*)"native", MPI_INFO_NULL);

  /*--- For now, create a temp 1D buffer to read the data from file. ---*/

  const int bufSize = nBlock*blocklen[0];
  Restart_Data = new passivedouble[bufSize];

  /*--- Collective call for all ranks to read from their view simultaneously. ---*/

  MPI_File_read_all(fhw, Restart_Data, bufSize, MPI_DOUBLE, &status);

  /*--- All ranks close the file after writing. ---*/

  MPI_File_close(&fhw);

  /*--- Free the derived datatype and release temp memory. ---*/

  MPI_Type_free(&filetype);

  delete [] blocklen;
  delete [] displace;

#endif

  if (nPointFile != geometry->GetGlobal_nPointDomain() &&
      config->GetKind_SU2() != SU2_COMPONENT::SU2_SOL) {
    InterpolateRestartData(geometry, config);
  }
}

void CSolver::InterpolateRestartData(const CGeometry *geometry, const CConfig *config) {

  if (geometry->GetGlobal_nPointDomain() == 0) return;

  if (size != SINGLE_NODE && size % 2)
    SU2_MPI::Error("Number of ranks must be multiple of 2.", CURRENT_FUNCTION);

  /* Challenges:
   *  - Do not use too much memory by gathering the restart data in all ranks.
   *  - Do not repeat too many computations in all ranks.
   * Solution?:
   *  - Build a local ADT for the domain points (not the restart points).
   *  - Find the closest target point for each donor, which does not match all targets.
   *  - "Diffuse" the data to neighbor points.
   *  Complexity is approx. Nlt + (Nlt + Nd) log(Nlt) where Nlt is the LOCAL number
   *  of target points and Nd the TOTAL number of donors. */

  const unsigned long nFields = Restart_Vars[1];
  const unsigned long nPointFile = Restart_Vars[2];
  const auto t0 = SU2_MPI::Wtime();
  auto nRecurse = 0;

  if (rank == MASTER_NODE) {
    cout << "\nThe number of points in the restart file (" << nPointFile << ") does not match "
            "the mesh (" << geometry->GetGlobal_nPointDomain() << ").\n"
            "A recursive nearest neighbor interpolation will be performed." << endl;
  }

  su2activematrix localVars(nPointDomain, nFields);
  localVars = su2double(0.0);
  {
  su2vector<uint8_t> isMapped(nPoint);
  isMapped = false;

  /*--- ADT of local target points. ---*/
  {
  const auto& coord = geometry->nodes->GetCoord();
  vector<unsigned long> index(nPointDomain);
  iota(index.begin(), index.end(), 0ul);

  CADTPointsOnlyClass adt(nDim, nPointDomain, coord.data(), index.data(), false);
  vector<unsigned long>().swap(index);

  /*--- Copy local donor restart data, which will circulate over all ranks. ---*/

  const auto partitioner = CLinearPartitioner(nPointFile,0);

  unsigned long nPointDonorMax = 0;
  for (int i=0; i<size; ++i)
    nPointDonorMax = max(nPointDonorMax, partitioner.GetSizeOnRank(i));

  su2activematrix sendBuf(nPointDonorMax, nFields);

  for (auto iPoint = 0ul; iPoint < nPointDonorMax; ++iPoint) {
    const auto iPointDonor = min(iPoint,partitioner.GetSizeOnRank(rank)-1ul);
    for (auto iVar = 0ul; iVar < nFields; ++iVar)
      sendBuf(iPoint,iVar) = Restart_Data[iPointDonor*nFields+iVar];
  }

  delete [] Restart_Data;
  Restart_Data = nullptr;

  /*--- Make room to receive donor data from other ranks, and to map it to target points. ---*/

  su2activematrix donorVars(nPointDonorMax, nFields);
  vector<su2double> donorDist(nPointDomain, 1e12);

  /*--- Circle over all ranks. ---*/

  const int dst = (rank+1) % size; // send to next
  const int src = (rank-1+size) % size; // receive from prev.
  const int count = sendBuf.size();

  for (int iStep = 0; iStep < size; ++iStep) {

    swap(sendBuf, donorVars);

    if (iStep) {
      /*--- Odd ranks send and then receive, and vice versa. ---*/
      if (rank%2) SU2_MPI::Send(sendBuf.data(), count, MPI_DOUBLE, dst, 0, SU2_MPI::GetComm());
      else SU2_MPI::Recv(donorVars.data(), count, MPI_DOUBLE, src, 0, SU2_MPI::GetComm(), MPI_STATUS_IGNORE);

      if (rank%2==0) SU2_MPI::Send(sendBuf.data(), count, MPI_DOUBLE, dst, 0, SU2_MPI::GetComm());
      else SU2_MPI::Recv(donorVars.data(), count, MPI_DOUBLE, src, 0, SU2_MPI::GetComm(), MPI_STATUS_IGNORE);
    }

    /*--- Find the closest target for each donor. ---*/

    vector<su2double> targetDist(donorVars.rows());
    vector<unsigned long> iTarget(donorVars.rows());

    SU2_OMP_PARALLEL_(for schedule(dynamic,4*OMP_MIN_SIZE))
    for (auto iDonor = 0ul; iDonor < donorVars.rows(); ++iDonor) {
      int r=0;
      adt.DetermineNearestNode(donorVars[iDonor], targetDist[iDonor], iTarget[iDonor], r);
    }
    END_SU2_OMP_PARALLEL

    /*--- Keep the closest donor for each target (this is separate for OpenMP). ---*/

    for (auto iDonor = 0ul; iDonor < donorVars.rows(); ++iDonor) {
      const auto iPoint = iTarget[iDonor];
      const auto dist = targetDist[iDonor];

      if (dist < donorDist[iPoint]) {
        donorDist[iPoint] = dist;
        isMapped[iPoint] = true;
        for (auto iVar = 0ul; iVar < donorVars.cols(); ++iVar)
          localVars(iPoint,iVar) = donorVars(iDonor,iVar);
      }
    }
  }
  } // everything goes out of scope except "localVars" and "isMapped"

  /*--- Recursively diffuse the nearest neighbor data. ---*/

  auto nDonor = isMapped;
  bool done = false;

  SU2_OMP_PARALLEL
  while (!done) {
    SU2_OMP_FOR_DYN(roundUpDiv(nPointDomain,2*omp_get_num_threads()))
    for (auto iPoint = 0ul; iPoint < nPointDomain; ++iPoint) {
      /*--- Do not change points that are already interpolated. ---*/
      if (isMapped[iPoint]) continue;

      /*--- Boundaries to boundaries and domain to domain. ---*/
      const bool boundary_i = geometry->nodes->GetSolidBoundary(iPoint);

      for (const auto jPoint : geometry->nodes->GetPoints(iPoint)) {
        if (!isMapped[jPoint]) continue;
        if (boundary_i != geometry->nodes->GetSolidBoundary(jPoint)) continue;

        nDonor[iPoint]++;

        for (auto iVar = 0ul; iVar < localVars.cols(); ++iVar)
          localVars(iPoint,iVar) += localVars(jPoint,iVar);
      }

      if (nDonor[iPoint] > 0) {
        for (auto iVar = 0ul; iVar < localVars.cols(); ++iVar)
          localVars(iPoint,iVar) /= nDonor[iPoint];
        nDonor[iPoint] = true;
      }
    }
    END_SU2_OMP_FOR

    /*--- Repeat while all points are not mapped. ---*/

    SU2_OMP_MASTER {
      done = true;
      ++nRecurse;
    }
    END_SU2_OMP_MASTER

    bool myDone = true;

    SU2_OMP_FOR_STAT(16*OMP_MIN_SIZE)
    for (auto iPoint = 0ul; iPoint < nPointDomain; ++iPoint) {
      isMapped[iPoint] = nDonor[iPoint];
      myDone &= nDonor[iPoint];
    }
    END_SU2_OMP_FOR

    SU2_OMP_ATOMIC
    done &= myDone;

    SU2_OMP_BARRIER
  }
  END_SU2_OMP_PARALLEL

  } // everything goes out of scope except "localVars"

  /*--- Move to Restart_Data in ascending order of global index, which is how a matching restart would have been read. ---*/

  Restart_Data = new passivedouble[nPointDomain*nFields];
  Restart_Vars[2] = nPointDomain;

  int counter = 0;
  for (auto iPoint_Global = 0ul; iPoint_Global < geometry->GetGlobal_nPointDomain(); ++iPoint_Global) {
    const auto iPoint = geometry->GetGlobal_to_Local_Point(iPoint_Global);
    if (iPoint >= 0) {
      for (auto iVar = 0ul; iVar < nFields; ++iVar)
        Restart_Data[counter*nFields+iVar] = SU2_TYPE::GetValue(localVars(iPoint,iVar));
      counter++;
    }
  }

  if (rank == MASTER_NODE) {
    cout << "Number of recursions: " << nRecurse << ".\n"
            "Elapsed time: " << SU2_MPI::Wtime()-t0 << "s.\n" << endl;
  }
}

void CSolver::Read_SU2_Restart_Metadata(CGeometry *geometry, CConfig *config, bool adjoint, string val_filename) const {

  su2double AoA_ = config->GetAoA();
  su2double AoS_ = config->GetAoS();
  su2double BCThrust_ = config->GetInitial_BCThrust();
  su2double dCD_dCL_ = config->GetdCD_dCL();
  su2double dCMx_dCL_ = config->GetdCMx_dCL();
  su2double dCMy_dCL_ = config->GetdCMy_dCL();
  su2double dCMz_dCL_ = config->GetdCMz_dCL();
  su2double SPPressureDrop_ = config->GetStreamwise_Periodic_PressureDrop();
  string::size_type position;
  unsigned long InnerIter_ = 0;
  ifstream restart_file;

  /*--- Carry on with ASCII metadata reading. ---*/

  restart_file.open(val_filename.data(), ios::in);
  if (restart_file.fail()) {
    if (rank == MASTER_NODE) {
      cout << " Warning: There is no restart file (" << val_filename.data() << ")."<< endl;
      cout << " Computation will continue without updating metadata parameters." << endl;
    }
  }
  else {

    string text_line;

    /*--- Space for extra info (if any) ---*/

    while (getline (restart_file, text_line)) {

      /*--- External iteration ---*/

      position = text_line.find ("ITER=",0);
      if (position != string::npos) {
        // TODO: 'ITER=' has 5 chars, not 9!
        text_line.erase (0,9); InnerIter_ = atoi(text_line.c_str());
      }

      /*--- Angle of attack ---*/

      position = text_line.find ("AOA=",0);
      if (position != string::npos) {
        text_line.erase (0,4); AoA_ = atof(text_line.c_str());
      }

      /*--- Sideslip angle ---*/

      position = text_line.find ("SIDESLIP_ANGLE=",0);
      if (position != string::npos) {
        text_line.erase (0,15); AoS_ = atof(text_line.c_str());
      }

      /*--- BCThrust angle ---*/

      position = text_line.find ("INITIAL_BCTHRUST=",0);
      if (position != string::npos) {
        text_line.erase (0,17); BCThrust_ = atof(text_line.c_str());
      }

      /*--- dCD_dCL coefficient ---*/

      position = text_line.find ("DCD_DCL_VALUE=",0);
      if (position != string::npos) {
        text_line.erase (0,14); dCD_dCL_ = atof(text_line.c_str());
      }

      /*--- dCMx_dCL coefficient ---*/

      position = text_line.find ("DCMX_DCL_VALUE=",0);
      if (position != string::npos) {
        text_line.erase (0,15); dCMx_dCL_ = atof(text_line.c_str());
      }

      /*--- dCMy_dCL coefficient ---*/

      position = text_line.find ("DCMY_DCL_VALUE=",0);
      if (position != string::npos) {
        text_line.erase (0,15); dCMy_dCL_ = atof(text_line.c_str());
      }

      /*--- dCMz_dCL coefficient ---*/

      position = text_line.find ("DCMZ_DCL_VALUE=",0);
      if (position != string::npos) {
        text_line.erase (0,15); dCMz_dCL_ = atof(text_line.c_str());
      }

      /*--- Streamwise periodic pressure drop for prescribed massflow cases. ---*/

      position = text_line.find ("STREAMWISE_PERIODIC_PRESSURE_DROP=",0);
      if (position != string::npos) {
        // Erase the name from the line, 'STREAMWISE_PERIODIC_PRESSURE_DROP=' has 34 chars.
        text_line.erase (0,34); SPPressureDrop_ = atof(text_line.c_str());
      }

    }

    /*--- Close the restart meta file. ---*/

    restart_file.close();

  }


  /*--- Load the metadata. ---*/

  /*--- Angle of attack ---*/

  if (config->GetDiscard_InFiles() == false) {
    if ((config->GetAoA() != AoA_) && (rank == MASTER_NODE)) {
      cout.precision(6);
      cout <<"WARNING: AoA in the solution file (" << AoA_ << " deg.) +" << endl;
      cout << "         AoA offset in mesh file (" << config->GetAoA_Offset() << " deg.) = " << AoA_ + config->GetAoA_Offset() << " deg." << endl;
    }
    config->SetAoA(AoA_ + config->GetAoA_Offset());
  }

  else {
    if ((config->GetAoA() != AoA_) && (rank == MASTER_NODE))
      cout <<"WARNING: Discarding the AoA in the solution file." << endl;
  }

  /*--- Sideslip angle ---*/

  if (config->GetDiscard_InFiles() == false) {
    if ((config->GetAoS() != AoS_) && (rank == MASTER_NODE)) {
      cout.precision(6);
      cout <<"WARNING: AoS in the solution file (" << AoS_ << " deg.) +" << endl;
      cout << "         AoS offset in mesh file (" << config->GetAoS_Offset() << " deg.) = " << AoS_ + config->GetAoS_Offset() << " deg." << endl;
    }
    config->SetAoS(AoS_ + config->GetAoS_Offset());
  }
  else {
    if ((config->GetAoS() != AoS_) && (rank == MASTER_NODE))
      cout <<"WARNING: Discarding the AoS in the solution file." << endl;
  }

  /*--- BCThrust ---*/

  if (config->GetDiscard_InFiles() == false) {
    if ((config->GetInitial_BCThrust() != BCThrust_) && (rank == MASTER_NODE))
      cout <<"WARNING: SU2 will use the initial BC Thrust provided in the solution file: " << BCThrust_ << " lbs." << endl;
    config->SetInitial_BCThrust(BCThrust_);
  }
  else {
    if ((config->GetInitial_BCThrust() != BCThrust_) && (rank == MASTER_NODE))
      cout <<"WARNING: Discarding the BC Thrust in the solution file." << endl;
  }


  if (config->GetDiscard_InFiles() == false) {

    if ((config->GetdCD_dCL() != dCD_dCL_) && (rank == MASTER_NODE))
      cout <<"WARNING: SU2 will use the dCD/dCL provided in the direct solution file: " << dCD_dCL_ << "." << endl;
    config->SetdCD_dCL(dCD_dCL_);

    if ((config->GetdCMx_dCL() != dCMx_dCL_) && (rank == MASTER_NODE))
      cout <<"WARNING: SU2 will use the dCMx/dCL provided in the direct solution file: " << dCMx_dCL_ << "." << endl;
    config->SetdCMx_dCL(dCMx_dCL_);

    if ((config->GetdCMy_dCL() != dCMy_dCL_) && (rank == MASTER_NODE))
      cout <<"WARNING: SU2 will use the dCMy/dCL provided in the direct solution file: " << dCMy_dCL_ << "." << endl;
    config->SetdCMy_dCL(dCMy_dCL_);

    if ((config->GetdCMz_dCL() != dCMz_dCL_) && (rank == MASTER_NODE))
      cout <<"WARNING: SU2 will use the dCMz/dCL provided in the direct solution file: " << dCMz_dCL_ << "." << endl;
    config->SetdCMz_dCL(dCMz_dCL_);

  }

  else {

    if ((config->GetdCD_dCL() != dCD_dCL_) && (rank == MASTER_NODE))
      cout <<"WARNING: Discarding the dCD/dCL in the direct solution file." << endl;

    if ((config->GetdCMx_dCL() != dCMx_dCL_) && (rank == MASTER_NODE))
      cout <<"WARNING: Discarding the dCMx/dCL in the direct solution file." << endl;

    if ((config->GetdCMy_dCL() != dCMy_dCL_) && (rank == MASTER_NODE))
      cout <<"WARNING: Discarding the dCMy/dCL in the direct solution file." << endl;

    if ((config->GetdCMz_dCL() != dCMz_dCL_) && (rank == MASTER_NODE))
      cout <<"WARNING: Discarding the dCMz/dCL in the direct solution file." << endl;

  }

  if (config->GetDiscard_InFiles() == false) {
    if ((config->GetStreamwise_Periodic_PressureDrop() != SPPressureDrop_) && (rank == MASTER_NODE))
      cout <<"WARNING: SU2 will use the STREAMWISE_PERIODIC_PRESSURE_DROP provided in the direct solution file: " << std::setprecision(16) << SPPressureDrop_ << endl;
    config->SetStreamwise_Periodic_PressureDrop(SPPressureDrop_);
  }
  else {
    if ((config->GetStreamwise_Periodic_PressureDrop() != SPPressureDrop_) && (rank == MASTER_NODE))
      cout <<"WARNING: Discarding the STREAMWISE_PERIODIC_PRESSURE_DROP in the direct solution file." << endl;
  }

  /*--- External iteration ---*/

  if ((config->GetDiscard_InFiles() == false) && (!adjoint || (adjoint && config->GetRestart())))
    config->SetExtIter_OffSet(InnerIter_);

}

void CSolver::ComputeVertexTractions(CGeometry *geometry, const CConfig *config){

  /*--- Compute the constant factor to dimensionalize pressure and shear stress. ---*/
  const su2double *Velocity_ND, *Velocity_Real;
  su2double Density_ND,  Density_Real, Velocity2_Real, Velocity2_ND;
  su2double factor;

  unsigned short iDim;

  // Check whether the problem is viscous
  bool viscous_flow = config->GetViscous();

  // Parameters for the calculations
  su2double Pn = 0.0;
  su2double auxForce[3] = {1.0, 0.0, 0.0};

  unsigned short iMarker;
  unsigned long iVertex, iPoint;
  const su2double* iNormal;

  su2double Pressure_Inf = config->GetPressure_FreeStreamND();

  Velocity_Real = config->GetVelocity_FreeStream();
  Density_Real  = config->GetDensity_FreeStream();

  Velocity_ND = config->GetVelocity_FreeStreamND();
  Density_ND  = config->GetDensity_FreeStreamND();

  Velocity2_Real = GeometryToolbox::SquaredNorm(nDim, Velocity_Real);
  Velocity2_ND   = GeometryToolbox::SquaredNorm(nDim, Velocity_ND);

  factor = Density_Real * Velocity2_Real / ( Density_ND * Velocity2_ND );

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

    /*--- If this is defined as a wall ---*/
    if (!config->GetSolid_Wall(iMarker)) continue;

    // Loop over the vertices
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {

      // Recover the point index
      iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
      // Get the normal at the vertex: this normal goes inside the fluid domain.
      iNormal = geometry->vertex[iMarker][iVertex]->GetNormal();

      /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
      if (geometry->nodes->GetDomain(iPoint)) {

        // Retrieve the values of pressure
        Pn = base_nodes->GetPressure(iPoint);

        // Calculate tn in the fluid nodes for the inviscid term --> Units of force (non-dimensional).
        for (iDim = 0; iDim < nDim; iDim++)
          auxForce[iDim] = -(Pn-Pressure_Inf)*iNormal[iDim];

        // Calculate tn in the fluid nodes for the viscous term
        if (viscous_flow) {
          su2double Viscosity = base_nodes->GetLaminarViscosity(iPoint);
          su2double Tau[3][3];
          CNumerics::ComputeStressTensor(nDim, Tau, base_nodes->GetVelocityGradient(iPoint), Viscosity);
          for (iDim = 0; iDim < nDim; iDim++) {
            auxForce[iDim] += GeometryToolbox::DotProduct(nDim, Tau[iDim], iNormal);
          }
        }

        // Redimensionalize the forces
        for (iDim = 0; iDim < nDim; iDim++) {
          VertexTraction[iMarker][iVertex][iDim] = factor * auxForce[iDim];
        }
      }
      else{
        for (iDim = 0; iDim < nDim; iDim++) {
          VertexTraction[iMarker][iVertex][iDim] = 0.0;
        }
      }
    }
  }

}

void CSolver::RegisterVertexTractions(CGeometry *geometry, const CConfig *config){

  unsigned short iMarker, iDim;
  unsigned long iVertex, iPoint;

  /*--- Loop over all the markers ---*/
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

    /*--- If this is defined as a wall ---*/
    if (!config->GetSolid_Wall(iMarker)) continue;

    /*--- Loop over the vertices ---*/
    SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {

      /*--- Recover the point index ---*/
      iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

      /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
      if (!geometry->nodes->GetDomain(iPoint)) continue;

      /*--- Register the vertex traction as output ---*/
      for (iDim = 0; iDim < nDim; iDim++) {
        AD::RegisterOutput(VertexTraction[iMarker][iVertex][iDim]);
      }
    }
    END_SU2_OMP_FOR
  }

}

void CSolver::SetVertexTractionsAdjoint(CGeometry *geometry, const CConfig *config){

  unsigned short iMarker, iDim;
  unsigned long iVertex, iPoint;

  /*--- Loop over all the markers ---*/
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

    /*--- If this is defined as a wall ---*/
    if (!config->GetSolid_Wall(iMarker)) continue;

    /*--- Loop over the vertices ---*/
    SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {

      /*--- Recover the point index ---*/
      iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

      /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
      if (!geometry->nodes->GetDomain(iPoint)) continue;

      /*--- Set the adjoint of the vertex traction from the value received ---*/
      for (iDim = 0; iDim < nDim; iDim++) {
        SU2_TYPE::SetDerivative(VertexTraction[iMarker][iVertex][iDim],
                                SU2_TYPE::GetValue(VertexTractionAdjoint[iMarker][iVertex][iDim]));
      }
    }
    END_SU2_OMP_FOR
  }

}

void CSolver::ComputeResidual_Multizone(const CGeometry *geometry, const CConfig *config){

  SU2_OMP_PARALLEL {

  /*--- Set Residuals to zero ---*/
  SU2_OMP_MASTER
  for (unsigned short iVar = 0; iVar < nVar; iVar++){
    Residual_BGS[iVar] = 0.0;
    Residual_Max_BGS[iVar] = 0.0;
  }
  END_SU2_OMP_MASTER

  vector<su2double> resMax(nVar,0.0), resRMS(nVar,0.0);
  vector<const su2double*> coordMax(nVar,nullptr);
  vector<unsigned long> idxMax(nVar,0);

  /*--- Set the residuals and BGSSolution_k to solution for next multizone outer iteration. ---*/
  SU2_OMP_FOR_STAT(roundUpDiv(nPoint,2*omp_get_num_threads()))
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {
    const su2double domain = (iPoint < nPointDomain);
    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
      const su2double Res = (base_nodes->Get_BGSSolution(iPoint,iVar) - base_nodes->Get_BGSSolution_k(iPoint,iVar))*domain;

      /*--- Update residual information for current thread. ---*/
      resRMS[iVar] += Res*Res;
      if (fabs(Res) > resMax[iVar]) {
        resMax[iVar] = fabs(Res);
        idxMax[iVar] = iPoint;
        coordMax[iVar] = geometry->nodes->GetCoord(iPoint);
      }
    }
  }
  END_SU2_OMP_FOR

  /*--- Reduce residual information over all threads in this rank. ---*/
  SU2_OMP_CRITICAL
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    Residual_BGS[iVar] += resRMS[iVar];
    AddRes_Max_BGS(iVar, resMax[iVar], geometry->nodes->GetGlobalIndex(idxMax[iVar]), coordMax[iVar]);
  }
  END_SU2_OMP_CRITICAL
  SU2_OMP_BARRIER

  SetResidual_BGS(geometry, config);

  }
  END_SU2_OMP_PARALLEL
}

void CSolver::BasicLoadRestart(CGeometry *geometry, const CConfig *config, const string& filename, unsigned long skipVars) {

  /*--- Read and store the restart metadata. ---*/

//  Read_SU2_Restart_Metadata(geometry[MESH_0], config, true, filename);

  /*--- Read the restart data from either an ASCII or binary SU2 file. ---*/

  if (config->GetRead_Binary_Restart()) {
    Read_SU2_Restart_Binary(geometry, config, filename);
  } else {
    Read_SU2_Restart_ASCII(geometry, config, filename);
  }

  /*--- Load data from the restart into correct containers. ---*/

  unsigned long iPoint_Global_Local = 0;

  for (auto iPoint_Global = 0ul; iPoint_Global < geometry->GetGlobal_nPointDomain(); iPoint_Global++ ) {

    /*--- Retrieve local index. If this node from the restart file lives
     on the current processor, we will load and instantiate the vars. ---*/

    const auto iPoint_Local = geometry->GetGlobal_to_Local_Point(iPoint_Global);

    if (iPoint_Local > -1) {

      /*--- We need to store this point's data, so jump to the correct
       offset in the buffer of data from the restart file and load it. ---*/

      const auto index = iPoint_Global_Local*Restart_Vars[1] + skipVars;

      for (auto iVar = 0u; iVar < nVar; iVar++) {
        base_nodes->SetSolution(iPoint_Local, iVar, Restart_Data[index+iVar]);
      }

      iPoint_Global_Local++;
    }

  }

  /*--- Delete the class memory that is used to load the restart. ---*/

  delete [] Restart_Vars;  Restart_Vars = nullptr;
  delete [] Restart_Data;  Restart_Data = nullptr;

  /*--- Detect a wrong solution file ---*/

  if (iPoint_Global_Local != nPointDomain) {
    SU2_MPI::Error(string("The solution file ") + filename + string(" doesn't match with the mesh file!\n") +
                   string("It could be empty lines at the end of the file."), CURRENT_FUNCTION);
  }
}

