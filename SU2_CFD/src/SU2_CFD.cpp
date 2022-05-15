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
  bool dry_run = false;
  int num_threads = omp_get_max_threads();
  bool use_thread_mult = false;
  std::string filename = "default.cfg";

  /*--- Command line parsing ---*/

  CLI::App app{"SU2 v7.3.1 \"Blackbird\", The Open-Source CFD Code"};
  app.add_flag("-d,--dryrun", dry_run, "Enable dry run mode.\n"
                                       "Only execute preprocessing steps using a dummy geometry.");
  app.add_option("-t,--threads", num_threads, "Number of OpenMP threads per MPI rank.");
  app.add_flag("--thread_multiple", use_thread_mult, "Request MPI_THREAD_MULTIPLE thread support.");
  app.add_option("configfile", filename, "A config file.")->check(CLI::ExistingFile);

  CLI11_PARSE(app, argc, argv)

  /*--- MPI initialization, and buffer setting ---*/

  SU2_MPI::Init(&argc, &argv);
  SU2_MPI::Comm MPICommunicator = SU2_MPI::GetComm();

  /*--- Uncomment the following line if runtime NaN catching is desired. ---*/
  // feenableexcept(FE_INVALID | FE_OVERFLOW | FE_DIVBYZERO );

  /*--- Create a pointer to the main SU2 Driver ---*/

  CDriver* driver = nullptr;

  /*--- Load in the number of zones and spatial dimensions in the mesh file (If no config
   file is specified, default.cfg is used) ---*/
  strcpy(config_file_name, filename.c_str());

  /*--- Read the name and format of the input mesh file to get from the mesh
   file the number of zones and dimensions from the numerical grid (required
   for variables allocation). ---*/

  const CConfig config(config_file_name, SU2_COMPONENT::SU2_CFD);
  const unsigned short nZone = config.GetnZone();

  /*--- First, given the basic information about the number of zones and the
   solver types from the config, instantiate the appropriate driver for the problem
   and perform all the preprocessing. ---*/

  if (dry_run) {

    /*--- Dry Run. ---*/
    driver = new CDummyDriver(config_file_name, nZone, MPICommunicator);

  }
  else {

    /*--- Generic single zone problem: instantiate the single zone driver class. ---*/
    if (nZone != 1)
      SU2_MPI::Error("The required solver doesn't support multizone simulations", CURRENT_FUNCTION);
    else {
      driver = new CSinglezoneDriver(config_file_name, nZone, MPICommunicator);
    }

  }

  /*--- Launch the main external loop of the solver. ---*/

  driver->StartSolver();

  /*--- Postprocess all the containers, close history file, exit SU2. ---*/

  driver->Postprocessing();

  delete driver;

  /*--- Finalize MPI parallelization. ---*/
  SU2_MPI::Finalize();

  return EXIT_SUCCESS;

}
