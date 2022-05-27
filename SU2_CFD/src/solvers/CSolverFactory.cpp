/*!
 * \file CSolverFactory.cpp
 * \brief Main subroutines for CSolverFactoryclass.
 * \author T. Albring
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
#include "../../include/solvers/CSolverFactory.hpp"
#include "../../include/solvers/CEulerSolver.hpp"
#include "../../include/solvers/CTemplateSolver.hpp"

map<const CSolver*, SolverMetaData> CSolverFactory::allocatedSolvers;

CSolver** CSolverFactory::CreateSolverContainer(MAIN_SOLVER kindMainSolver, CConfig *config, CGeometry *geometry, int iMGLevel){

  CSolver** solver;

  solver = new CSolver*[MAX_SOLS]();

  switch (kindMainSolver) {
    case MAIN_SOLVER::TEMPLATE_SOLVER:
      solver[FLOW_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::TEMPLATE, solver, geometry, config, iMGLevel);
      break;
    case MAIN_SOLVER::EULER:
      solver[FLOW_SOL] = CreateSubSolver(SUB_SOLVER_TYPE::EULER, solver, geometry, config, iMGLevel);
      break;
     default:
      solver = nullptr;
  }

  return solver;

}

CSolver* CSolverFactory::CreateSubSolver(SUB_SOLVER_TYPE kindSolver, CSolver **solver, CGeometry *geometry, CConfig *config, int iMGLevel){

  CSolver *genericSolver = nullptr;

  SolverMetaData metaData;

  metaData.solverType = kindSolver;

  switch (kindSolver) {
    case SUB_SOLVER_TYPE::EULER:
      genericSolver = CreateFlowSolver(kindSolver, solver, geometry, config, iMGLevel);
      if (!config->GetNewtonKrylov() || config->GetDiscrete_Adjoint() || config->GetContinuous_Adjoint())
        metaData.integrationType = INTEGRATION_TYPE::MULTIGRID;
      else
        metaData.integrationType = INTEGRATION_TYPE::NEWTON;
      break;
    case SUB_SOLVER_TYPE::TEMPLATE:
      genericSolver = new CTemplateSolver(geometry, config);
      metaData.integrationType = INTEGRATION_TYPE::SINGLEGRID;
      break;
    default:
      SU2_MPI::Error("No proper allocation found for requested sub solver", CURRENT_FUNCTION);
      break;
  }

  if (genericSolver != nullptr)
    allocatedSolvers[genericSolver] = metaData;

  return genericSolver;

}

CSolver* CSolverFactory::CreateFlowSolver(SUB_SOLVER_TYPE kindFlowSolver, CSolver **solver,  CGeometry *geometry, CConfig *config, int iMGLevel){

  CSolver *flowSolver = nullptr;

  switch (kindFlowSolver) {
    case SUB_SOLVER_TYPE::EULER:
      flowSolver = new CEulerSolver(geometry, config, iMGLevel);
      flowSolver->Preprocessing(geometry, solver, config, iMGLevel, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
      break;
    default:
      SU2_MPI::Error("Flow solver not found", CURRENT_FUNCTION);
      break;
  }

  return flowSolver;
}
