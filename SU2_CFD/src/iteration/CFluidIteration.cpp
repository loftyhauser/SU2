/*!
 * \file CFluidIteration.cpp
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

#include "../../include/iteration/CFluidIteration.hpp"
#include "../../include/output/COutput.hpp"

void CFluidIteration::Preprocess(COutput* output, CIntegration**** integration, CGeometry**** geometry,
                                 CSolver***** solver, CNumerics****** numerics, CConfig** config,
                                 CSurfaceMovement** surface_movement, CVolumetricMovement*** grid_movement,
                                 CFreeFormDefBox*** FFDBox, unsigned short val_iZone, unsigned short val_iInst) {
  unsigned long TimeIter = config[val_iZone]->GetTimeIter();

  bool fsi = config[val_iZone]->GetFSI_Simulation();
  unsigned long OuterIter = config[val_iZone]->GetOuterIter();

  /*--- Set the initial condition for FSI problems with subiterations ---*/
  /*--- This is done only in the first block subiteration.---*/
  /*--- From then on, the solver reuses the partially converged solution obtained in the previous subiteration ---*/
  if (fsi && !config[val_iZone]->GetDiscrete_Adjoint() && (OuterIter == 0)) {
    solver[val_iZone][val_iInst][MESH_0][FLOW_SOL]->SetInitialCondition(
        geometry[val_iZone][val_iInst], solver[val_iZone][val_iInst], config[val_iZone], TimeIter);
  }

}

void CFluidIteration::Iterate(COutput* output, CIntegration**** integration, CGeometry**** geometry,
                              CSolver***** solver, CNumerics****** numerics, CConfig** config,
                              CSurfaceMovement** surface_movement, CVolumetricMovement*** grid_movement,
                              CFreeFormDefBox*** FFDBox, unsigned short val_iZone, unsigned short val_iInst) {

  const bool frozen_visc = (config[val_iZone]->GetContinuous_Adjoint() && config[val_iZone]->GetFrozen_Visc_Cont()) ||
                           (config[val_iZone]->GetDiscrete_Adjoint() && config[val_iZone]->GetFrozen_Visc_Disc());
  const bool disc_adj = (config[val_iZone]->GetDiscrete_Adjoint());

  /* --- Setting up iteration values depending on if this is a
   steady or an unsteady simulation */

  /*--- Update global parameters ---*/

  MAIN_SOLVER main_solver = MAIN_SOLVER::NONE;

  switch (config[val_iZone]->GetKind_Solver()) {
    case MAIN_SOLVER::EULER:
    case MAIN_SOLVER::DISC_ADJ_EULER:
    case MAIN_SOLVER::INC_EULER:
    case MAIN_SOLVER::DISC_ADJ_INC_EULER:
      main_solver = MAIN_SOLVER::EULER;
      break;

    case MAIN_SOLVER::NAVIER_STOKES:
    case MAIN_SOLVER::DISC_ADJ_NAVIER_STOKES:
    case MAIN_SOLVER::INC_NAVIER_STOKES:
    case MAIN_SOLVER::DISC_ADJ_INC_NAVIER_STOKES:
      main_solver = MAIN_SOLVER::NAVIER_STOKES;
      break;

    case MAIN_SOLVER::RANS:
    case MAIN_SOLVER::DISC_ADJ_RANS:
    case MAIN_SOLVER::INC_RANS:
    case MAIN_SOLVER::DISC_ADJ_INC_RANS:
      main_solver = MAIN_SOLVER::RANS;
      break;

    default:
      break;
  }
  config[val_iZone]->SetGlobalParam(main_solver, RUNTIME_FLOW_SYS);

  /*--- Solve the Euler, Navier-Stokes or Reynolds-averaged Navier-Stokes (RANS) equations (one iteration) ---*/

  integration[val_iZone][val_iInst][FLOW_SOL]->MultiGrid_Iteration(geometry, solver, numerics, config, RUNTIME_FLOW_SYS,
                                                                   val_iZone, val_iInst);

  /*--- If the flow integration is not fully coupled, run the various single grid integrations. ---*/

  if ((main_solver == MAIN_SOLVER::RANS) && !frozen_visc) {
    /*--- Solve the turbulence model ---*/

    config[val_iZone]->SetGlobalParam(main_solver, RUNTIME_TURB_SYS);
    integration[val_iZone][val_iInst][TURB_SOL]->SingleGrid_Iteration(geometry, solver, numerics, config,
                                                                      RUNTIME_TURB_SYS, val_iZone, val_iInst);

    /*--- Solve transition model ---*/

    if (config[val_iZone]->GetKind_Trans_Model() == TURB_TRANS_MODEL::LM) {
      config[val_iZone]->SetGlobalParam(main_solver, RUNTIME_TRANS_SYS);
      integration[val_iZone][val_iInst][TRANS_SOL]->SingleGrid_Iteration(geometry, solver, numerics, config,
                                                                         RUNTIME_TRANS_SYS, val_iZone, val_iInst);
    }
  }

  /*--- Adapt the CFL number using an exponential progression with under-relaxation approach. ---*/

  if ((config[val_iZone]->GetCFL_Adapt() == YES) && (!disc_adj)) {
    SU2_OMP_PARALLEL
    solver[val_iZone][val_iInst][MESH_0][FLOW_SOL]->AdaptCFLNumber(geometry[val_iZone][val_iInst],
                                                                   solver[val_iZone][val_iInst], config[val_iZone]);
    END_SU2_OMP_PARALLEL
  }

}

void CFluidIteration::Update(COutput* output, CIntegration**** integration, CGeometry**** geometry, CSolver***** solver,
                             CNumerics****** numerics, CConfig** config, CSurfaceMovement** surface_movement,
                             CVolumetricMovement*** grid_movement, CFreeFormDefBox*** FFDBox, unsigned short val_iZone,
                             unsigned short val_iInst) {
  unsigned short iMesh;

  /*--- Dual time stepping strategy ---*/

  if ((config[val_iZone]->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST) ||
      (config[val_iZone]->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND)) {
    /*--- Update dual time solver on all mesh levels ---*/

    for (iMesh = 0; iMesh <= config[val_iZone]->GetnMGLevels(); iMesh++) {
      integration[val_iZone][val_iInst][FLOW_SOL]->SetDualTime_Solver(geometry[val_iZone][val_iInst][iMesh],
                                                                      solver[val_iZone][val_iInst][iMesh][FLOW_SOL],
                                                                      config[val_iZone], iMesh);

      integration[val_iZone][val_iInst][FLOW_SOL]->SetDualTime_Geometry(geometry[val_iZone][val_iInst][iMesh],
                                                                        solver[val_iZone][val_iInst][iMesh][MESH_SOL],
                                                                        config[val_iZone], iMesh);

      integration[val_iZone][val_iInst][FLOW_SOL]->SetConvergence(false);
    }

    /*--- Update dual time solver for the turbulence model ---*/

    if ((config[val_iZone]->GetKind_Solver() == MAIN_SOLVER::RANS) || (config[val_iZone]->GetKind_Solver() == MAIN_SOLVER::DISC_ADJ_RANS) ||
        (config[val_iZone]->GetKind_Solver() == MAIN_SOLVER::INC_RANS) ||
        (config[val_iZone]->GetKind_Solver() == MAIN_SOLVER::DISC_ADJ_INC_RANS)) {
      integration[val_iZone][val_iInst][TURB_SOL]->SetDualTime_Solver(geometry[val_iZone][val_iInst][MESH_0],
                                                                      solver[val_iZone][val_iInst][MESH_0][TURB_SOL],
                                                                      config[val_iZone], MESH_0);
      integration[val_iZone][val_iInst][TURB_SOL]->SetConvergence(false);
    }

    /*--- Update dual time solver for the transition model ---*/

    if (config[val_iZone]->GetKind_Trans_Model() == TURB_TRANS_MODEL::LM) {
      integration[val_iZone][val_iInst][TRANS_SOL]->SetDualTime_Solver(geometry[val_iZone][val_iInst][MESH_0],
                                                                       solver[val_iZone][val_iInst][MESH_0][TRANS_SOL],
                                                                       config[val_iZone], MESH_0);
      integration[val_iZone][val_iInst][TRANS_SOL]->SetConvergence(false);
    }

  }
}

bool CFluidIteration::Monitor(COutput* output, CIntegration**** integration, CGeometry**** geometry,
                              CSolver***** solver, CNumerics****** numerics, CConfig** config,
                              CSurfaceMovement** surface_movement, CVolumetricMovement*** grid_movement,
                              CFreeFormDefBox*** FFDBox, unsigned short val_iZone, unsigned short val_iInst) {
  bool StopCalc = false;

  StopTime = SU2_MPI::Wtime();

  UsedTime = StopTime - StartTime;

  if (config[val_iZone]->GetMultizone_Problem() || config[val_iZone]->GetSinglezone_Driver()) {
    output->SetHistory_Output(geometry[val_iZone][INST_0][MESH_0], solver[val_iZone][INST_0][MESH_0], config[val_iZone],
                              config[val_iZone]->GetTimeIter(), config[val_iZone]->GetOuterIter(),
                              config[val_iZone]->GetInnerIter());
  }

  /*--- If convergence was reached --*/
  StopCalc = output->GetConvergence();

  /* --- Checking convergence of Fixed CL mode to target CL, and perform finite differencing if needed  --*/

  if (config[val_iZone]->GetFixed_CL_Mode()) {
    StopCalc = MonitorFixed_CL(output, geometry[val_iZone][INST_0][MESH_0], solver[val_iZone][INST_0][MESH_0],
                               config[val_iZone]);
  }

  return StopCalc;
}

void CFluidIteration::Postprocess(COutput* output, CIntegration**** integration, CGeometry**** geometry,
                                  CSolver***** solver, CNumerics****** numerics, CConfig** config,
                                  CSurfaceMovement** surface_movement, CVolumetricMovement*** grid_movement,
                                  CFreeFormDefBox*** FFDBox, unsigned short val_iZone, unsigned short val_iInst) {

  /*--- Temporary: enable only for single-zone driver. This should be removed eventually when generalized. ---*/
  if (config[val_iZone]->GetSinglezone_Driver()) {

    /*--- Compute the tractions at the vertices ---*/
    solver[val_iZone][val_iInst][MESH_0][FLOW_SOL]->ComputeVertexTractions(geometry[val_iZone][val_iInst][MESH_0],
                                                                           config[val_iZone]);
  }
}

void CFluidIteration::Solve(COutput* output, CIntegration**** integration, CGeometry**** geometry, CSolver***** solver,
                            CNumerics****** numerics, CConfig** config, CSurfaceMovement** surface_movement,
                            CVolumetricMovement*** grid_movement, CFreeFormDefBox*** FFDBox, unsigned short val_iZone,
                            unsigned short val_iInst) {
  /*--- Boolean to determine if we are running a static or dynamic case ---*/
  bool steady = !config[val_iZone]->GetTime_Domain();

  unsigned long Inner_Iter, nInner_Iter = config[val_iZone]->GetnInner_Iter();
  bool StopCalc = false;

  /*--- Synchronization point before a single solver iteration.
        Compute the wall clock time required. ---*/

  StartTime = SU2_MPI::Wtime();

  /*--- Preprocess the solver ---*/
  Preprocess(output, integration, geometry, solver, numerics, config, surface_movement, grid_movement, FFDBox,
             val_iZone, INST_0);

  /*--- For steady-state flow simulations, we need to loop over ExtIter for the number of time steps ---*/
  /*--- However, ExtIter is the number of FSI iterations, so nIntIter is used in this case ---*/

  for (Inner_Iter = 0; Inner_Iter < nInner_Iter; Inner_Iter++) {
    config[val_iZone]->SetInnerIter(Inner_Iter);

    /*--- Run a single iteration of the solver ---*/
    Iterate(output, integration, geometry, solver, numerics, config, surface_movement, grid_movement, FFDBox, val_iZone,
            INST_0);

    /*--- Monitor the pseudo-time ---*/
    StopCalc = Monitor(output, integration, geometry, solver, numerics, config, surface_movement, grid_movement, FFDBox,
                       val_iZone, INST_0);

    /*--- Output files at intermediate iterations if the problem is single zone ---*/

    if (singlezone && steady) {
      Output(output, geometry, solver, config, Inner_Iter, StopCalc, val_iZone, val_iInst);
    }

    /*--- If the iteration has converged, break the loop ---*/
    if (StopCalc) break;
  }

  if (multizone && steady) {
    Output(output, geometry, solver, config, config[val_iZone]->GetOuterIter(), StopCalc, val_iZone, val_iInst);

    /*--- Set the convergence to false (to make sure outer subiterations converge) ---*/

      integration[val_iZone][INST_0][FLOW_SOL]->SetConvergence(false);
  }
}

void CFluidIteration::InitializeVortexDistribution(unsigned long& nVortex, vector<su2double>& x0, vector<su2double>& y0,
                                                   vector<su2double>& vort_strength, vector<su2double>& r_core) {
  /*--- Read in Vortex Distribution ---*/
  std::string line;
  std::ifstream file;
  su2double x_temp, y_temp, vort_strength_temp, r_core_temp;
  file.open("vortex_distribution.txt");
  /*--- In case there is no vortex file ---*/
  if (file.fail()) {
    SU2_MPI::Error("There is no vortex data file!!", CURRENT_FUNCTION);
  }

  // Ignore line containing the header
  getline(file, line);
  // Read in the information of the vortices (xloc, yloc, lambda(strength), eta(size, gradient))
  while (file.good()) {
    getline(file, line);
    std::stringstream ss(line);
    if (line.size() != 0) {  // ignore blank lines if they exist.
      ss >> x_temp;
      ss >> y_temp;
      ss >> vort_strength_temp;
      ss >> r_core_temp;
      x0.push_back(x_temp);
      y0.push_back(y_temp);
      vort_strength.push_back(vort_strength_temp);
      r_core.push_back(r_core_temp);
    }
  }
  file.close();
  // number of vortices
  nVortex = x0.size();
}

bool CFluidIteration::MonitorFixed_CL(COutput *output, CGeometry *geometry, CSolver **solver, CConfig *config) {

  CSolver* flow_solver= solver[FLOW_SOL];

  bool fixed_cl_convergence = flow_solver->FixedCL_Convergence(config, output->GetConvergence());

  /* --- If Fixed CL mode has ended and Finite Differencing has started: --- */

  if (flow_solver->GetStart_AoA_FD() && flow_solver->GetIter_Update_AoA() == config->GetInnerIter()){

    /* --- Print convergence history and volume files since fixed CL mode has converged--- */
    if (rank == MASTER_NODE) output->PrintConvergenceSummary();

    output->SetResult_Files(geometry, config, solver,
                            config->GetInnerIter(), true);

    /* --- Set finite difference mode in config (disables output) --- */
    config->SetFinite_Difference_Mode(true);
  }

  /* --- Set convergence based on fixed CL convergence  --- */
  return fixed_cl_convergence;
}

