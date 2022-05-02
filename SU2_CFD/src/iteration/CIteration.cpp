/*!
 * \file iteration_structure.cpp
 * \brief Main subroutines used by SU2_CFD
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

#include "../../include/iteration/CIteration.hpp"

#include "../../include/output/COutput.hpp"
#include "../../include/solvers/CFEASolver.hpp"

void CIteration::SetGrid_Movement(CGeometry** geometry, CSurfaceMovement* surface_movement,
                                  CVolumetricMovement* grid_movement, CSolver*** solver, CConfig* config,
                                  unsigned long IntIter, unsigned long TimeIter) {
  unsigned short Kind_Grid_Movement = config->GetKind_GridMovement();

  unsigned short val_iZone = config->GetiZone();

  /*--- Perform mesh movement depending on specified type ---*/
  switch (Kind_Grid_Movement) {
    case RIGID_MOTION:

      if (rank == MASTER_NODE) cout << endl << " Performing rigid mesh transformation." << endl;

      /*--- Move each node in the volume mesh using the specified type
       of rigid mesh motion. These routines also compute analytic grid
       velocities for the fine mesh. ---*/

      grid_movement->Rigid_Translation(geometry[MESH_0], config, val_iZone, TimeIter);
      grid_movement->Rigid_Plunging(geometry[MESH_0], config, val_iZone, TimeIter);
      grid_movement->Rigid_Pitching(geometry[MESH_0], config, val_iZone, TimeIter);
      grid_movement->Rigid_Rotation(geometry[MESH_0], config, val_iZone, TimeIter);

      /*--- Update the multigrid structure after moving the finest grid,
       including computing the grid velocities on the coarser levels. ---*/

      grid_movement->UpdateMultiGrid(geometry, config);

      break;

      /*--- Already initialized in the static mesh movement routine at driver level. ---*/
    case STEADY_TRANSLATION:
    case ROTATING_FRAME:
      break;
  }

}

void CIteration::SetMesh_Deformation(CGeometry** geometry, CSolver** solver, CNumerics*** numerics, CConfig* config,
                                     RECORDING kind_recording) {
  if (!config->GetDeform_Mesh()) return;

  /*--- Perform the elasticity mesh movement ---*/

  bool wasActive = false;
  if ((kind_recording != RECORDING::MESH_DEFORM) && !config->GetMultizone_Problem()) {
    /*--- In a primal run, AD::TapeActive returns a false ---*/
    /*--- In any other recordings, the tape is passive during the deformation. ---*/
    wasActive = AD::BeginPassive();
  }

  /*--- Set the stiffness of each element mesh into the mesh numerics ---*/

  solver[MESH_SOL]->SetMesh_Stiffness(geometry, numerics[MESH_SOL], config);

  /*--- Deform the volume grid around the new boundary locations ---*/

  solver[MESH_SOL]->DeformMesh(geometry, numerics[MESH_SOL], config);

  /*--- Continue recording. ---*/
  AD::EndPassive(wasActive);
}

void CIteration::Output(COutput* output, CGeometry**** geometry, CSolver***** solver, CConfig** config,
                        unsigned long InnerIter, bool StopCalc, unsigned short val_iZone, unsigned short val_iInst) {
  output->SetResult_Files(geometry[val_iZone][INST_0][MESH_0], config[val_iZone], solver[val_iZone][INST_0][MESH_0],
                          InnerIter);
}
