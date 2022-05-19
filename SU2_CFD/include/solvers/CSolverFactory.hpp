/*!
 * \file CSolverFactory.hpp
 * \brief Headers of the CSolverFactory class
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
#pragma once

#include "../../../Common/include/option_structure.hpp"

/*!
 * \brief Enum of different sub solvers the main solver can use. There is not a 1-to-1 correspondence between the actual classes
 * and the types of sub solvers, as one class can be used for several sub solvers.
 */
enum class SUB_SOLVER_TYPE {
  CONT_ADJ_EULER,          /*!< \brief Continuous Adjoint Euler solver  */
  CONT_ADJ_NAVIER_STOKES,  /*!< \brief Continuous Adjoint Navier Stokes solver  */
  CONT_ADJ_TURB,           /*!< \brief Continuous Adjoint Turbulent solver  */
  BASELINE,                /*!< \brief Baseline solver  */
  TEMPLATE,                /*!< \brief Template solver  */
  DISC_ADJ_FLOW,           /*!< \brief Discrete adjoint flow solver */
  DISC_ADJ_TURB,           /*!< \brief Discrete adjoint turbulence solver */
  DISC_ADJ_HEAT,           /*!< \brief Discrete adjoint heat solver */
  EULER,                   /*!< \brief Compressible Euler solver */
  NAVIER_STOKES,           /*!< \brief Compressible Navier-Stokes solver */
  TRANSITION,              /*!< \brief Transition model solver*/
  TURB_SA,                 /*!< \brief SA turbulence model solver */
  TURB_SST,                /*!< \brief SST turbulence model solver */
  TURB,                    /*!< \brief Turbulence model solver */
  NONE
};

enum class INTEGRATION_TYPE{
  MULTIGRID,
  NEWTON,
  SINGLEGRID,
  DEFAULT,
  NONE
};

struct SolverMetaData{
  SUB_SOLVER_TYPE  solverType        = SUB_SOLVER_TYPE::NONE;
  INTEGRATION_TYPE integrationType   = INTEGRATION_TYPE::NONE;
};

class CSolver;
class CGeometry;
class CConfig;

class CSolverFactory {

private:

  static std::map<const CSolver*, SolverMetaData> allocatedSolvers;

  /*!
   * \brief Create a flow solver
   * \param[in] kindFlowSolver - Kind of flow solver
   * \param[in] solver         - The solver container
   * \param[in] geometry       - The geometry definition
   * \param[in] config         - The configuration
   * \param[in] iMGLevel       - The multigrid level
   * \return                   - A pointer to the allocated flow solver
   */
  static CSolver* CreateFlowSolver(SUB_SOLVER_TYPE kindFlowSolver, CSolver **solver, CGeometry *geometry, CConfig *config, int iMGLevel);

  /*!
   * \brief Generic routine to create a solver
   * \param[in] kindSolver    - Kind of solver
   * \param[in] solver        - The solver container
   * \param[in] geometry      - The geometry definition
   * \param[in] config        - The configuration
   * \param[in] iMGLevel      - The multigrid level
   * \return                  - A pointer to the allocated solver
   */
  static CSolver* CreateSubSolver(SUB_SOLVER_TYPE kindSolver, CSolver **solver, CGeometry *geometry, CConfig *config, int iMGLevel);

public:

  /*!
   * \brief Deleted constructor to avoid creating instances of this class
   */
  CSolverFactory() = delete;

  /*!
   * \brief Create the solver container by allocating the primary solver
   * and secondary solvers like heat solver, turbulent solver etc
   * \param[in] kindSolver    - The kind of primary solver
   * \param[in] config        - The configuration
   * \param[in] geometry      - The geometry definition
   * \param[in] iMGLevel      - The multigrid level
   * \return                  - Pointer to the allocated solver array
   */
  static CSolver** CreateSolverContainer(MAIN_SOLVER kindSolver, CConfig *config, CGeometry *geometry, int iMGLevel);


  /*!
   * \brief Return a sub solver object that contains information about the solver allocated at a specific memory address
   * \param[in] solver - Address of the solver
   * \return sub solver info struct.
   */
  static SolverMetaData GetSolverMeta(const CSolver* solver) { return allocatedSolvers.at(solver); }

  /*!
   * \brief Clear the solver meta data
   */
  static void ClearSolverMeta() { allocatedSolvers.clear(); }

};
