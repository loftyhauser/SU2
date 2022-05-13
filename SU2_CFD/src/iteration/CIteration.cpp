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

void CIteration::Output(COutput* output, CGeometry**** geometry, CSolver***** solver, CConfig** config,
                        unsigned long InnerIter, bool StopCalc, unsigned short val_iZone, unsigned short val_iInst) {
  output->SetResult_Files(geometry[val_iZone][INST_0][MESH_0], config[val_iZone], solver[val_iZone][INST_0][MESH_0],
                          InnerIter);
}
