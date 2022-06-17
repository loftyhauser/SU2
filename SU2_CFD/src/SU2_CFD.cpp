/*!
 * \file SU2_CFD.cpp
 * \brief Main file of the SU2 Computational Fluid Dynamics code
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

#include "../include/SU2_CFD.hpp"

/* Include file, needed for the runtime NaN catching. You also have to include feenableexcept(...) below. */
//#include <fenv.h>

using namespace std;

int main(int argc, char *argv[]) {

  char config_file_name[MAX_STRING_SIZE];
  std::string filename = "default.cfg";

  /*--- Command line parsing ---*/

  CLI::App app{"SU2 v7.3.1 \"Blackbird\", The Open-Source CFD Code"};
  app.add_option("configfile", filename, "A config file.")->check(CLI::ExistingFile);

  CLI11_PARSE(app, argc, argv)

  /*--- Create a pointer to the main SU2 Driver ---*/

  CDriver* driver = nullptr;

  /*--- Load in the number of zones and spatial dimensions in the mesh file (If no config
   file is specified, default.cfg is used) ---*/
  strcpy(config_file_name, filename.c_str());

  /*--- Read the name and format of the input mesh file to get from the mesh
   file the number of zones and dimensions from the numerical grid (required
   for variables allocation). ---*/

  const CConfig config(config_file_name, SU2_COMPONENT::SU2_CFD);

  /*--- Generic single zone problem: instantiate the single zone driver class. ---*/

  driver = new CSinglezoneDriver(config_file_name, 1);

  /*--- Launch the main external loop of the solver. ---*/

  driver->StartSolver();

  /*--- Postprocess all the containers, close history file, exit SU2. ---*/

  driver->Postprocessing();

  delete driver;

  return EXIT_SUCCESS;

}
