/*!
 * \file CIterationFactory.cpp
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

#include "../../include/iteration/CIterationFactory.hpp"
#include "../../include/iteration/CIteration.hpp"
#include "../../include/iteration/CFluidIteration.hpp"

CIteration* CIterationFactory::CreateIteration(MAIN_SOLVER kindSolver, const CConfig* config){

  CIteration *iteration = nullptr;

  const auto rank = SU2_MPI::GetRank();

  /*--- Loop over all zones and instantiate the physics iteration. ---*/

  switch (kindSolver) {

    case MAIN_SOLVER::EULER: case MAIN_SOLVER::NAVIER_STOKES: case MAIN_SOLVER::RANS:
        if (rank == MASTER_NODE)
          cout << "Euler/Navier-Stokes/RANS fluid iteration." << endl;
        iteration = new CFluidIteration(config);
      break;

    case MAIN_SOLVER::NONE: case MAIN_SOLVER::TEMPLATE_SOLVER: case MAIN_SOLVER::MULTIPHYSICS:
      SU2_MPI::Error("No iteration found for specified solver.", CURRENT_FUNCTION);
      break;
  }

  return iteration;
}
