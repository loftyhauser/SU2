/*!
 * \file CDriver.cpp
 * \brief The main subroutines for driving single or multi-zone problems.
 * \author T. Economon, H. Kline, R. Sanchez, F. Palacios
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

#include "../../include/drivers/CDriver.hpp"
#include "../../include/definition_structure.hpp"

#include "../../../Common/include/geometry/CDummyGeometry.hpp"
#include "../../../Common/include/geometry/CPhysicalGeometry.hpp"
#include "../../../Common/include/geometry/CMultiGridGeometry.hpp"

#include "../../include/solvers/CSolverFactory.hpp"

#include "../../include/output/COutputFactory.hpp"
#include "../../include/output/COutput.hpp"

#include "../../include/variables/CEulerVariable.hpp"

#include "../../include/numerics/template.hpp"
#include "../../include/numerics/flow/convection/roe.hpp"
#include "../../include/numerics/flow/convection/fvs.hpp"
#include "../../include/numerics/flow/convection/cusp.hpp"
#include "../../include/numerics/flow/convection/hllc.hpp"
#include "../../include/numerics/flow/convection/ausm_slau.hpp"

#include "../../include/integration/CIntegrationFactory.hpp"

#include "../../include/iteration/CIterationFactory.hpp"

#include "../../../Common/include/parallelization/omp_structure.hpp"

#include <cassert>

#include <fenv.h>

CDriver::CDriver(char* confFile, unsigned short val_nZone, SU2_Comm MPICommunicator, bool dummy_geo) :
  config_file_name(confFile), StartTime(0.0), StopTime(0.0), UsedTime(0.0),
  TimeIter(0), nZone(val_nZone), StopCalc(false), dry_run(dummy_geo) {

  SU2_MPI::SetComm(MPICommunicator);

  rank = SU2_MPI::GetRank();
  size = SU2_MPI::GetSize();

  /*--- Start timer to track preprocessing for benchmarking. ---*/

  StartTime = SU2_MPI::Wtime();

  /*--- Initialize containers with null --- */

  SetContainers_Null();

  /*--- Preprocessing of the config files. ---*/

  Input_Preprocessing(config_container, driver_config);

  /*--- Retrieve dimension from mesh file ---*/

  nDim = CConfig::GetnDim(config_container[ZONE_0]->GetMesh_FileName(),
                          config_container[ZONE_0]->GetMesh_FileFormat());

  /*--- Output preprocessing ---*/

  Output_Preprocessing(config_container, driver_config, output_container, driver_output);


  for (iZone = 0; iZone < nZone; iZone++) {

    /*--- Read the number of instances for each zone ---*/

    nInst[iZone] = 1;

    geometry_container[iZone]    = new CGeometry**    [nInst[iZone]] ();
    iteration_container[iZone]   = new CIteration*    [nInst[iZone]] ();
    solver_container[iZone]      = new CSolver***     [nInst[iZone]] ();
    integration_container[iZone] = new CIntegration** [nInst[iZone]] ();
    numerics_container[iZone]    = new CNumerics****  [nInst[iZone]] ();

    for (iInst = 0; iInst < nInst[iZone]; iInst++) {

      config_container[iZone]->SetiInst(iInst);

      /*--- Preprocessing of the geometry for all zones. In this routine, the edge-
       based data structure is constructed, i.e. node and cell neighbors are
       identified and linked, face areas and volumes of the dual mesh cells are
       computed, and the multigrid levels are created using an agglomeration procedure. ---*/

      Geometrical_Preprocessing(config_container[iZone], geometry_container[iZone][iInst], dry_run);

    }
  }

  /*--- Before we proceed with the zone loop we have to compute the wall distances.
     * This computation depends on all zones at once. ---*/
  for (iZone = 0; iZone < nZone; iZone++) {

    for (iInst = 0; iInst < nInst[iZone]; iInst++){

      /*--- Definition of the solver class: solver_container[#ZONES][#INSTANCES][#MG_GRIDS][#EQ_SYSTEMS].
       The solver classes are specific to a particular set of governing equations,
       and they contain the subroutines with instructions for computing each spatial
       term of the PDE, i.e. loops over the edges to compute convective and viscous
       fluxes, loops over the nodes to compute source terms, and routines for
       imposing various boundary condition type for the PDE. ---*/

      Solver_Preprocessing(config_container[iZone], geometry_container[iZone][iInst], solver_container[iZone][iInst]);

      /*--- Definition of the numerical method class:
       numerics_container[#ZONES][#INSTANCES][#MG_GRIDS][#EQ_SYSTEMS][#EQ_TERMS].
       The numerics class contains the implementation of the numerical methods for
       evaluating convective or viscous fluxes between any two nodes in the edge-based
       data structure (centered, upwind, galerkin), as well as any source terms
       (piecewise constant reconstruction) evaluated in each dual mesh volume. ---*/

      Numerics_Preprocessing(config_container[iZone], geometry_container[iZone][iInst],
                             solver_container[iZone][iInst], numerics_container[iZone][iInst]);

      /*--- Definition of the integration class: integration_container[#ZONES][#INSTANCES][#EQ_SYSTEMS].
       The integration class orchestrates the execution of the spatial integration
       subroutines contained in the solver class (including multigrid) for computing
       the residual at each node, R(U) and then integrates the equations to a
       steady state or time-accurately. ---*/

      Integration_Preprocessing(config_container[iZone], solver_container[iZone][iInst][MESH_0],
                                integration_container[iZone][iInst]);

      /*--- Instantiate the type of physics iteration to be executed within each zone. For
       example, one can execute the same physics across multiple zones (mixing plane),
       different physics in different zones (fluid-structure interaction), or couple multiple
       systems tightly within a single zone by creating a new iteration class (e.g., RANS). ---*/

      Iteration_Preprocessing(config_container[iZone], iteration_container[iZone][iInst]);

    }

  }

  PythonInterface_Preprocessing(config_container, geometry_container, solver_container);


  /*--- Preprocessing time is reported now, but not included in the next compute portion. ---*/

  StopTime = SU2_MPI::Wtime();

  /*--- Compute/print the total time for performance benchmarking. ---*/

  UsedTime = StopTime-StartTime;
  UsedTimePreproc    = UsedTime;
  UsedTimeCompute    = 0.0;
  UsedTimeOutput     = 0.0;
  IterCount          = 0;
  OutputCount        = 0;
  MDOFs              = 0.0;
  MDOFsDomain        = 0.0;
  Mpoints            = 0.0;
  MpointsDomain      = 0.0;
  for (iZone = 0; iZone < nZone; iZone++) {
    Mpoints       += geometry_container[iZone][INST_0][MESH_0]->GetGlobal_nPoint()/(1.0e6);
    MpointsDomain += geometry_container[iZone][INST_0][MESH_0]->GetGlobal_nPointDomain()/(1.0e6);
    MDOFs         += DOFsPerPoint*geometry_container[iZone][INST_0][MESH_0]->GetGlobal_nPoint()/(1.0e6);
    MDOFsDomain   += DOFsPerPoint*geometry_container[iZone][INST_0][MESH_0]->GetGlobal_nPointDomain()/(1.0e6);
  }

  /*--- Reset timer for compute/output performance benchmarking. ---*/

  StopTime = SU2_MPI::Wtime();

  /*--- Compute/print the total time for performance benchmarking. ---*/

  UsedTime = StopTime-StartTime;
  UsedTimePreproc = UsedTime;

  /*--- Reset timer for compute performance benchmarking. ---*/

  StartTime = SU2_MPI::Wtime();

}

void CDriver::SetContainers_Null(){

  /*--- Create pointers to all of the classes that may be used throughout
   the SU2_CFD code. In general, the pointers are instantiated down a
   hierarchy over all zones, multigrid levels, equation sets, and equation
   terms as described in the comments below. ---*/

  ConvHist_file                  = nullptr;
  iteration_container            = nullptr;
  output_container               = nullptr;
  integration_container          = nullptr;
  geometry_container             = nullptr;
  solver_container               = nullptr;
  numerics_container             = nullptr;
  config_container               = nullptr;
  interface_types                = nullptr;
  nInst                          = nullptr;

  /*--- Definition and of the containers for all possible zones. ---*/

  iteration_container            = new CIteration**[nZone] ();
  solver_container               = new CSolver****[nZone] ();
  integration_container          = new CIntegration***[nZone] ();
  numerics_container             = new CNumerics*****[nZone] ();
  config_container               = new CConfig*[nZone] ();
  geometry_container             = new CGeometry***[nZone] ();
  interface_types                = new unsigned short*[nZone] ();
  output_container               = new COutput*[nZone] ();
  nInst                          = new unsigned short[nZone] ();
  driver_config                  = nullptr;
  driver_output                  = nullptr;

  for (iZone = 0; iZone < nZone; iZone++) {
    interface_types[iZone] = new unsigned short[nZone];
    nInst[iZone] = 1;
  }

  strcpy(runtime_file_name, "runtime.dat");

}


void CDriver::Postprocessing() {

  const bool wrt_perf = config_container[ZONE_0]->GetWrt_Performance();

    /*--- Output some information to the console. ---*/

  if (rank == MASTER_NODE) {

    /*--- Print out the number of non-physical points and reconstructions ---*/

    if (config_container[ZONE_0]->GetNonphysical_Points() > 0)
      cout << "Warning: there are " << config_container[ZONE_0]->GetNonphysical_Points() << " non-physical points in the solution." << endl;
    if (config_container[ZONE_0]->GetNonphysical_Reconstr() > 0)
      cout << "Warning: " << config_container[ZONE_0]->GetNonphysical_Reconstr() << " reconstructed states for upwinding are non-physical." << endl;
  }

  if (rank == MASTER_NODE)
    cout << endl <<"------------------------- Solver Postprocessing -------------------------" << endl;

  for (iZone = 0; iZone < nZone; iZone++) {
    for (iInst = 0; iInst < nInst[iZone]; iInst++){
      Numerics_Postprocessing(numerics_container[iZone], solver_container[iZone][iInst],
          geometry_container[iZone][iInst], config_container[iZone], iInst);
    }
    delete [] numerics_container[iZone];
  }
  delete [] numerics_container;
  if (rank == MASTER_NODE) cout << "Deleted CNumerics container." << endl;

  for (iZone = 0; iZone < nZone; iZone++) {
    for (iInst = 0; iInst < nInst[iZone]; iInst++){
      Integration_Postprocessing(integration_container[iZone],
          geometry_container[iZone][iInst],
          config_container[iZone],
          iInst);
    }
    delete [] integration_container[iZone];
  }
  delete [] integration_container;
  if (rank == MASTER_NODE) cout << "Deleted CIntegration container." << endl;

  for (iZone = 0; iZone < nZone; iZone++) {
    for (iInst = 0; iInst < nInst[iZone]; iInst++){
      Solver_Postprocessing(solver_container[iZone],
          geometry_container[iZone][iInst],
          config_container[iZone],
          iInst);
    }
    delete [] solver_container[iZone];
  }
  delete [] solver_container;
  if (rank == MASTER_NODE) cout << "Deleted CSolver container." << endl;

  for (iZone = 0; iZone < nZone; iZone++) {
    for (iInst = 0; iInst < nInst[iZone]; iInst++)
      delete iteration_container[iZone][iInst];
    delete [] iteration_container[iZone];
  }
  delete [] iteration_container;
  if (rank == MASTER_NODE) cout << "Deleted CIteration container." << endl;

  if (interface_types != nullptr) {
    for (iZone = 0; iZone < nZone; iZone++)
      delete [] interface_types[iZone];
    delete [] interface_types;
  }

  for (iZone = 0; iZone < nZone; iZone++) {
    if (geometry_container[iZone] != nullptr) {
      for (iInst = 0; iInst < nInst[iZone]; iInst++){
        for (unsigned short iMGlevel = 0; iMGlevel < config_container[iZone]->GetnMGLevels()+1; iMGlevel++)
          delete geometry_container[iZone][iInst][iMGlevel];
        delete [] geometry_container[iZone][iInst];
      }
      delete [] geometry_container[iZone];
    }
  }
  delete [] geometry_container;
  if (rank == MASTER_NODE) cout << "Deleted CGeometry container." << endl;

  /*--- Deallocate config container ---*/
  if (config_container!= nullptr) {
    for (iZone = 0; iZone < nZone; iZone++)
      delete config_container[iZone];
    delete [] config_container;
  }
  delete driver_config;
  if (rank == MASTER_NODE) cout << "Deleted CConfig container." << endl;

  delete [] nInst;
  if (rank == MASTER_NODE) cout << "Deleted nInst container." << endl;

  /*--- Deallocate output container ---*/

  if (output_container!= nullptr) {
    for (iZone = 0; iZone < nZone; iZone++)
      delete output_container[iZone];
    delete [] output_container;
  }

  delete driver_output;

  if (rank == MASTER_NODE) cout << "Deleted COutput class." << endl;

  if (rank == MASTER_NODE) cout << "-------------------------------------------------------------------------" << endl;


  /*--- Stop the timer and output the final performance summary. ---*/

  StopTime = SU2_MPI::Wtime();

  UsedTime = StopTime-StartTime;
  UsedTimeCompute += UsedTime;

  if ((rank == MASTER_NODE) && (wrt_perf)) {
    su2double TotalTime = UsedTimePreproc + UsedTimeCompute + UsedTimeOutput;
    cout.precision(6);
    cout << endl << endl <<"-------------------------- Performance Summary --------------------------" << endl;
    cout << "Simulation totals:" << endl;
    cout << setw(25) << "Wall-clock time (hrs):" << setw(12) << (TotalTime)/(60.0*60.0) << " | ";
    cout << setw(20) << "Core-hrs:" << setw(12) << size*TotalTime/(60.0*60.0) << endl;
    cout << setw(25) << "Cores:" << setw(12) << size << " | ";
    cout << setw(20) << "DOFs/point:" << setw(12) << DOFsPerPoint << endl;
    cout << setw(25) << "Points/core:" << setw(12) << 1.0e6*MpointsDomain/size << " | ";
    cout << setw(20) << "Ghost points/core:" << setw(12) << 1.0e6*(Mpoints-MpointsDomain)/size << endl;
    cout << setw(25) << "Ghost/Owned Point Ratio:" << setw(12) << (Mpoints-MpointsDomain)/MpointsDomain << " | " << endl;
    cout << endl;
    cout << "Preprocessing phase:" << endl;
    cout << setw(25) << "Preproc. Time (s):"  << setw(12)<< UsedTimePreproc << " | ";
    cout << setw(20) << "Preproc. Time (%):" << setw(12)<< ((UsedTimePreproc * 100.0) / (TotalTime)) << endl;
    cout << endl;
    cout << "Compute phase:" << endl;
    cout << setw(25) << "Compute Time (s):"  << setw(12)<< UsedTimeCompute << " | ";
    cout << setw(20) << "Compute Time (%):" << setw(12)<< ((UsedTimeCompute * 100.0) / (TotalTime)) << endl;
    cout << setw(25) << "Iteration count:"  << setw(12)<< IterCount << " | ";
    if (IterCount != 0) {
      cout << setw(20) << "Avg. s/iter:" << setw(12)<< UsedTimeCompute/IterCount << endl;
      cout << setw(25) << "Core-s/iter/Mpoints:" << setw(12)<< size*UsedTimeCompute/IterCount/Mpoints << " | ";
      cout << setw(20) << "Mpoints/s:" << setw(12)<< Mpoints*IterCount/UsedTimeCompute << endl;
    } else cout << endl;
    cout << endl;
    cout << "Output phase:" << endl;
    cout << setw(25) << "Output Time (s):"  << setw(12)<< UsedTimeOutput << " | ";
    cout << setw(20) << "Output Time (%):" << setw(12)<< ((UsedTimeOutput * 100.0) / (TotalTime)) << endl;
    cout << setw(25) << "Output count:" << setw(12)<< OutputCount << " | ";
    if (OutputCount != 0) {
      cout << setw(20)<< "Avg. s/output:" << setw(12)<< UsedTimeOutput/OutputCount << endl;
      if (BandwidthSum > 0) {
        cout << setw(25)<< "Restart Aggr. BW (MB/s):" << setw(12)<< BandwidthSum/OutputCount << " | ";
        cout << setw(20)<< "MB/s/core:" << setw(12)<< BandwidthSum/OutputCount/size << endl;
      }
    } else cout << endl;
    cout << "-------------------------------------------------------------------------" << endl;
    cout << endl;
  }

  /*--- Exit the solver cleanly ---*/

  if (rank == MASTER_NODE)
    cout << endl <<"------------------------- Exit Success (SU2_CFD) ------------------------" << endl << endl;

}


void CDriver::Input_Preprocessing(CConfig **&config, CConfig *&driver_config) {

  char zone_file_name[MAX_STRING_SIZE];

  /*--- Initialize the configuration of the driver ---*/

  driver_config = new CConfig(config_file_name, SU2_COMPONENT::SU2_CFD, false);

  for (iZone = 0; iZone < nZone; iZone++) {

    if (rank == MASTER_NODE){
      cout  << endl << "Parsing config file for zone " << iZone << endl;
    }
    /*--- Definition of the configuration option class for all zones. In this
     constructor, the input configuration file is parsed and all options are
     read and stored. ---*/

    if (driver_config->GetnConfigFiles() > 0){

      strcpy(zone_file_name, driver_config->GetConfigFilename(iZone).c_str());
      config[iZone] = new CConfig(driver_config, zone_file_name, SU2_COMPONENT::SU2_CFD, iZone, nZone, true);
    }
    else{
      config[iZone] = new CConfig(driver_config, config_file_name, SU2_COMPONENT::SU2_CFD, iZone, nZone, true);
    }

    /*--- Set the MPI communicator ---*/

    config[iZone]->SetMPICommunicator(SU2_MPI::GetComm());
  }


}

void CDriver::Geometrical_Preprocessing(CConfig* config, CGeometry **&geometry, bool dummy){

  if (!dummy){
    if (rank == MASTER_NODE)
      cout << endl <<"------------------- Geometry Preprocessing ( Zone " << config->GetiZone() <<" ) -------------------" << endl;

    Geometrical_Preprocessing_FVM(config, geometry);
  } else {
    if (rank == MASTER_NODE)
      cout << endl <<"-------------------------- Using Dummy Geometry -------------------------" << endl;

    unsigned short iMGlevel;

    geometry = new CGeometry*[config->GetnMGLevels()+1] ();

      for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
        geometry[iMGlevel] = new CDummyGeometry(config);
      }

    nDim = geometry[MESH_0]->GetnDim();
  }

  /*--- Computation of positive surface area in the z-plane which is used for
     the calculation of force coefficient (non-dimensionalization). ---*/

  geometry[MESH_0]->SetPositive_ZArea(config);

  /*--- Set the actuator disk boundary conditions, if necessary. ---*/

  for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
    geometry[iMesh]->MatchActuator_Disk(config);
  }

  /*--- If we have any periodic markers in this calculation, we must
       match the periodic points found on both sides of the periodic BC.
       Note that the current implementation requires a 1-to-1 matching of
       periodic points on the pair of periodic faces after the translation
       or rotation is taken into account. ---*/

  if ((config->GetnMarker_Periodic() != 0)) {
    for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {

      /*--- Note that we loop over pairs of periodic markers individually
           so that repeated nodes on adjacent periodic faces are properly
           accounted for in multiple places. ---*/

      for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
        geometry[iMesh]->MatchPeriodic(config, iPeriodic);
      }

      /*--- For Streamwise Periodic flow, find a unique reference node on the dedicated inlet marker. ---*/
      if (config->GetKind_Streamwise_Periodic() != ENUM_STREAMWISE_PERIODIC::NONE)
        geometry[iMesh]->FindUniqueNode_PeriodicBound(config);

      /*--- Initialize the communication framework for the periodic BCs. ---*/
      geometry[iMesh]->PreprocessPeriodicComms(geometry[iMesh], config);

    }
  }

  /*--- If activated by the compile directive, perform a partition analysis. ---*/
#if PARTITION
  if (!dummy){
    else Partition_Analysis(geometry[MESH_0], config);
  }
#endif

  /*--- Check if Euler & Symmetry markers are straight/plane. This information
        is used in the Euler & Symmetry boundary routines. ---*/
  if((config_container[iZone]->GetnMarker_Euler() != 0 ||
     config_container[iZone]->GetnMarker_SymWall() != 0)) {

    if (rank == MASTER_NODE)
      cout << "Checking if Euler & Symmetry markers are straight/plane:" << endl;

    for (iMesh = 0; iMesh <= config_container[iZone]->GetnMGLevels(); iMesh++)
      geometry_container[iZone][iInst][iMesh]->ComputeSurf_Straightness(config_container[iZone], (iMesh==MESH_0) );

  }

}

void CDriver::Geometrical_Preprocessing_FVM(CConfig *config, CGeometry **&geometry) {

  unsigned short iZone = config->GetiZone(), iMGlevel;
  unsigned short requestedMGlevels = config->GetnMGLevels();

  /*--- Definition of the geometry class to store the primal grid in the partitioning process.
   *    All ranks process the grid and call ParMETIS for partitioning ---*/

  CGeometry *geometry_aux = new CPhysicalGeometry(config, iZone, nZone);

  /*--- Set the dimension --- */

  nDim = geometry_aux->GetnDim();

  /*--- Allocate the memory of the current domain, and divide the grid
     between the ranks. ---*/

  geometry = new CGeometry *[config->GetnMGLevels()+1] ();

  /*--- Build the grid data structures using the ParMETIS coloring. ---*/

  geometry[MESH_0] = new CPhysicalGeometry(geometry_aux, config);

  /*--- Deallocate the memory of geometry_aux and solver_aux ---*/

  delete geometry_aux;

  /*--- Add the Send/Receive boundaries ---*/
  geometry[MESH_0]->SetSendReceive(config);

  /*--- Add the Send/Receive boundaries ---*/
  geometry[MESH_0]->SetBoundaries(config);

  /*--- Compute elements surrounding points, points surrounding points ---*/

  if (rank == MASTER_NODE) cout << "Setting point connectivity." << endl;
  geometry[MESH_0]->SetPoint_Connectivity();

  /*--- Renumbering points using Reverse Cuthill McKee ordering ---*/

  if (rank == MASTER_NODE) cout << "Renumbering points (Reverse Cuthill McKee Ordering)." << endl;
  geometry[MESH_0]->SetRCM_Ordering(config);

  /*--- recompute elements surrounding points, points surrounding points ---*/

  if (rank == MASTER_NODE) cout << "Recomputing point connectivity." << endl;
  geometry[MESH_0]->SetPoint_Connectivity();

  /*--- Compute elements surrounding elements ---*/

  if (rank == MASTER_NODE) cout << "Setting element connectivity." << endl;
  geometry[MESH_0]->SetElement_Connectivity();

  /*--- Check the orientation before computing geometrical quantities ---*/

  geometry[MESH_0]->SetBoundVolume();
  if (config->GetReorientElements()) {
    if (rank == MASTER_NODE) cout << "Checking the numerical grid orientation." << endl;
    geometry[MESH_0]->Check_IntElem_Orientation(config);
    geometry[MESH_0]->Check_BoundElem_Orientation(config);
  }

  /*--- Create the edge structure ---*/

  if (rank == MASTER_NODE) cout << "Identifying edges and vertices." << endl;
  geometry[MESH_0]->SetEdges();
  geometry[MESH_0]->SetVertex(config);

  /*--- Create the control volume structures ---*/

  if (rank == MASTER_NODE) cout << "Setting the control volume structure." << endl;
  SU2_OMP_PARALLEL {
    geometry[MESH_0]->SetControlVolume(config, ALLOCATE);
    geometry[MESH_0]->SetBoundControlVolume(config, ALLOCATE);
  }
  END_SU2_OMP_PARALLEL

  /*--- Visualize a dual control volume if requested ---*/

  if ((config->GetVisualize_CV() >= 0) &&
      (config->GetVisualize_CV() < (long)geometry[MESH_0]->GetGlobal_nPointDomain()))
    geometry[MESH_0]->VisualizeControlVolume(config);

  /*--- Identify closest normal neighbor ---*/

  if (rank == MASTER_NODE) cout << "Searching for the closest normal neighbors to the surfaces." << endl;
  geometry[MESH_0]->FindNormal_Neighbor(config);

  /*--- Store the global to local mapping. ---*/

  if (rank == MASTER_NODE) cout << "Storing a mapping from global to local point index." << endl;
  geometry[MESH_0]->SetGlobal_to_Local_Point();

  /*--- Compute the surface curvature ---*/

    if (rank == MASTER_NODE) cout << "Compute the surface curvature." << endl;
    geometry[MESH_0]->ComputeSurf_Curvature(config);

  /*--- Compute the global surface areas for all markers. ---*/

  geometry[MESH_0]->ComputeSurfaceAreaCfgFile(config);

  /*--- Check for periodicity and disable MG if necessary. ---*/

  if (rank == MASTER_NODE) cout << "Checking for periodicity." << endl;
  geometry[MESH_0]->Check_Periodicity(config);

  /*--- Compute mesh quality statistics on the fine grid. ---*/

    if (rank == MASTER_NODE)
      cout << "Computing mesh quality statistics for the dual control volumes." << endl;
    geometry[MESH_0]->ComputeMeshQualityStatistics(config);

  geometry[MESH_0]->SetMGLevel(MESH_0);
  if ((config->GetnMGLevels() != 0) && (rank == MASTER_NODE))
    cout << "Setting the multigrid structure." << endl;

  /*--- Loop over all the new grid ---*/

  for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {

    /*--- Create main agglomeration structure ---*/

    geometry[iMGlevel] = new CMultiGridGeometry(geometry[iMGlevel-1], config, iMGlevel);

    /*--- Compute points surrounding points. ---*/

    geometry[iMGlevel]->SetPoint_Connectivity(geometry[iMGlevel-1]);

    /*--- Create the edge structure ---*/

    geometry[iMGlevel]->SetEdges();
    geometry[iMGlevel]->SetVertex(geometry[iMGlevel-1], config);

    /*--- Create the control volume structures ---*/

    geometry[iMGlevel]->SetControlVolume(geometry[iMGlevel-1], ALLOCATE);
    geometry[iMGlevel]->SetBoundControlVolume(geometry[iMGlevel-1], ALLOCATE);
    geometry[iMGlevel]->SetCoord(geometry[iMGlevel-1]);

    /*--- Find closest neighbor to a surface point ---*/

    geometry[iMGlevel]->FindNormal_Neighbor(config);

    /*--- Store our multigrid index. ---*/

    geometry[iMGlevel]->SetMGLevel(iMGlevel);

    /*--- Protect against the situation that we were not able to complete
       the agglomeration for this level, i.e., there weren't enough points.
       We need to check if we changed the total number of levels and delete
       the incomplete CMultiGridGeometry object. ---*/

    if (config->GetnMGLevels() != requestedMGlevels) {
      delete geometry[iMGlevel];
      geometry[iMGlevel] = nullptr;
      break;
    }

  }

  if (config->GetWrt_MultiGrid()) geometry[MESH_0]->ColorMGLevels(config->GetnMGLevels(), geometry);

  /*--- For unsteady simulations, initialize the grid volumes
   and coordinates for previous solutions. Loop over all zones/grids ---*/


  /*--- Create the data structure for MPI point-to-point communications. ---*/

  for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
    geometry[iMGlevel]->PreprocessP2PComms(geometry[iMGlevel], config);


  /*--- Perform a few preprocessing routines and communications. ---*/

  for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {

    /*--- Compute the max length. ---*/

      if ((rank == MASTER_NODE) && (iMGlevel == MESH_0))
        cout << "Finding max control volume width." << endl;
      geometry[iMGlevel]->SetMaxLength(config);

    /*--- Communicate the number of neighbors. This is needed for
         some centered schemes and for multigrid in parallel. ---*/

    if ((rank == MASTER_NODE) && (size > SINGLE_NODE) && (iMGlevel == MESH_0))
      cout << "Communicating number of neighbors." << endl;
    geometry[iMGlevel]->InitiateComms(geometry[iMGlevel], config, NEIGHBORS);
    geometry[iMGlevel]->CompleteComms(geometry[iMGlevel], config, NEIGHBORS);
  }

}

void CDriver::Solver_Preprocessing(CConfig* config, CGeometry** geometry, CSolver ***&solver) {

  MAIN_SOLVER kindSolver = config->GetKind_Solver();

  if (rank == MASTER_NODE)
    cout << endl <<"-------------------- Solver Preprocessing ( Zone " << config->GetiZone() <<" ) --------------------" << endl;

  solver = new CSolver**[config->GetnMGLevels()+1] ();

  for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++){
    solver[iMesh] = CSolverFactory::CreateSolverContainer(kindSolver, config, geometry[iMesh], iMesh);
  }

  /*--- Count the number of DOFs per solution point. ---*/

  DOFsPerPoint = 0;

  for (unsigned int iSol = 0; iSol < MAX_SOLS; iSol++)
    if (solver[MESH_0][iSol]) DOFsPerPoint += solver[MESH_0][iSol]->GetnVar();

  /*--- Restart solvers, for FSI the geometry cannot be updated because the interpolation classes
   * should always use the undeformed mesh (otherwise the results would not be repeatable). ---*/

  if (!fsi) Solver_Restart(solver, geometry, config, true);

  /*--- Set up any necessary inlet profiles ---*/

  Inlet_Preprocessing(solver, geometry, config);

}

void CDriver::Inlet_Preprocessing(CSolver ***solver, CGeometry **geometry,
                                  CConfig *config) const {

  /*--- Adjust iteration number for unsteady restarts. ---*/

  const bool adjoint = config->GetDiscrete_Adjoint() || config->GetContinuous_Adjoint();

  int val_iter = 0;

  if (config->GetTime_Domain()) {
    val_iter = adjoint? config->GetUnst_AdjointIter() : config->GetRestart_Iter();
    val_iter -= 1;
    if (!adjoint && config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND)
      val_iter -= 1;
    if (!adjoint && !config->GetRestart()) val_iter = 0;
  }

  /*--- Load inlet profile files for any of the active solver containers.
   Note that these routines fill the fine grid data structures for the markers
   and restrict values down to all coarser MG levels. ---*/


    /*--- Uniform inlets or python-customized inlets ---*/

    /* --- Initialize quantities for inlet boundary
     * This routine does not check if the python wrapper is being used to
     * set custom boundary conditions.  This is intentional; the
     * default values for python custom BCs are initialized with the default
     * values specified in the config (avoiding non physical values) --- */

    for (unsigned short iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
      for(unsigned short iMarker=0; iMarker < config->GetnMarker_All(); iMarker++) {
        if (solver[iMesh][FLOW_SOL]) solver[iMesh][FLOW_SOL]->SetUniformInlet(config, iMarker);
        if (solver[iMesh][TURB_SOL]) solver[iMesh][TURB_SOL]->SetUniformInlet(config, iMarker);
      }
    }


}

void CDriver::Solver_Restart(CSolver ***solver, CGeometry **geometry,
                             CConfig *config, bool update_geo) {

  /*--- Check for restarts and use the LoadRestart() routines. ---*/

  const bool restart = config->GetRestart();
  const bool restart_flow = config->GetRestart_Flow();

  /*--- Adjust iteration number for unsteady restarts. ---*/

  int val_iter = 0;

  const bool adjoint = (config->GetDiscrete_Adjoint() || config->GetContinuous_Adjoint());
  const bool time_domain = config->GetTime_Domain();
  const bool dt_step_2nd = (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND) &&
                           !adjoint && time_domain;

  if (time_domain) {
    if (adjoint) val_iter = config->GetUnst_AdjointIter() - 1;
    else val_iter = config->GetRestart_Iter() - 1 - dt_step_2nd;
  }

  /*--- Restart direct solvers. ---*/

  if (restart || restart_flow) {
    for (auto iSol = 0u; iSol < MAX_SOLS; ++iSol) {
      auto sol = solver[MESH_0][iSol];
      if (sol && !sol->GetAdjoint()) {
        /*--- Note that the mesh solver always loads the most recent file (and not -2). ---*/
        SU2_OMP_PARALLEL_(if(sol->GetHasHybridParallel()))
        sol->LoadRestart(geometry, solver, config, val_iter + (iSol==MESH_SOL && dt_step_2nd), update_geo);
        END_SU2_OMP_PARALLEL
      }
    }
  }

  /*--- Restart adjoint solvers. ---*/

  if (restart) {
    if ((config->GetKind_Solver() == MAIN_SOLVER::TEMPLATE_SOLVER) && !config->GetFrozen_Visc_Cont()) {
      SU2_MPI::Error("A restart capability has not been implemented yet for this solver.\n"
                     "Please set RESTART_SOL= NO and try again.", CURRENT_FUNCTION);
    }

    for (auto iSol = 0u; iSol < MAX_SOLS; ++iSol) {
      auto sol = solver[MESH_0][iSol];
      if (sol && sol->GetAdjoint())
        sol->LoadRestart(geometry, solver, config, val_iter, update_geo);
    }
  }

}

void CDriver::Solver_Postprocessing(CSolver ****solver, CGeometry **geometry,
                                    CConfig *config, unsigned short val_iInst) {

  for (int iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
    for (unsigned int iSol = 0; iSol < MAX_SOLS; iSol++){
      delete solver[val_iInst][iMGlevel][iSol];
    }
    delete [] solver[val_iInst][iMGlevel];
  }
  delete [] solver[val_iInst];

  CSolverFactory::ClearSolverMeta();

}

void CDriver::Integration_Preprocessing(CConfig *config, CSolver **solver, CIntegration **&integration) const {

  if (rank == MASTER_NODE)
    cout << endl <<"----------------- Integration Preprocessing ( Zone " << config->GetiZone() <<" ) ------------------" << endl;

  MAIN_SOLVER kindMainSolver = config->GetKind_Solver();

  integration = CIntegrationFactory::CreateIntegrationContainer(kindMainSolver, solver);

}

void CDriver::Integration_Postprocessing(CIntegration ***integration, CGeometry **geometry, CConfig *config, unsigned short val_iInst) {

  for (unsigned int iSol = 0; iSol < MAX_SOLS; iSol++){
    delete integration[val_iInst][iSol];
  }

  delete [] integration[val_iInst];

}

void CDriver::Numerics_Preprocessing(CConfig *config, CGeometry **geometry, CSolver ***solver, CNumerics ****&numerics) const {

  if (rank == MASTER_NODE)
    cout << endl <<"------------------- Numerics Preprocessing ( Zone " << config->GetiZone() <<" ) -------------------" << endl;

  unsigned short iMGlevel, iSol,

  nVar_Template         = 0,
  nVar_Flow             = 0;

  numerics = new CNumerics***[config->GetnMGLevels()+1] ();

  bool compressible = false;
  bool ideal_gas = (config->GetKind_FluidModel() == STANDARD_AIR) || (config->GetKind_FluidModel() == IDEAL_GAS);
  bool roe_low_dissipation = (config->GetKind_RoeLowDiss() != NO_ROELOWDISS);

  /*--- Initialize some useful booleans ---*/
  bool euler, ns, turbulent, adj_euler, adj_ns, adj_turb;
  bool template_solver;

  euler = ns = turbulent = adj_euler = adj_ns = adj_turb = false;
  template_solver = false;

  /*--- Assign booleans ---*/
  switch (config->GetKind_Solver()) {
    case MAIN_SOLVER::TEMPLATE_SOLVER:
      template_solver = true; break;

    case MAIN_SOLVER::EULER:
      euler = compressible = true; break;

    default:
      break;

  }

  /*--- Number of variables for the template ---*/

  if (template_solver) nVar_Flow = solver[MESH_0][FLOW_SOL]->GetnVar();

  /*--- Number of variables for direct problem ---*/

  if (euler)        nVar_Flow = solver[MESH_0][FLOW_SOL]->GetnVar();
  if (ns)           nVar_Flow = solver[MESH_0][FLOW_SOL]->GetnVar();

  /*--- Definition of the Class for the numerical method:
    numerics_container[INSTANCE_LEVEL][MESH_LEVEL][EQUATION][EQ_TERM] ---*/

  for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
    numerics[iMGlevel] = new CNumerics** [MAX_SOLS];
    for (iSol = 0; iSol < MAX_SOLS; iSol++)
      numerics[iMGlevel][iSol] = new CNumerics* [MAX_TERMS*omp_get_max_threads()]();
  }

  /*--- Instantiate one numerics object per thread for each required term. ---*/

  for (int thread = 0; thread < omp_get_max_threads(); ++thread)
  {
  const int offset = thread * MAX_TERMS;

  const int conv_term = CONV_TERM + offset;
  const int visc_term = VISC_TERM + offset;

  const int source_first_term = SOURCE_FIRST_TERM + offset;
  const int source_second_term = SOURCE_SECOND_TERM + offset;

  const int conv_bound_term = CONV_BOUND_TERM + offset;

  /*--- Solver definition for the template problem ---*/
  if (template_solver) {

    /*--- Definition of the convective scheme for each equation and mesh level ---*/
    switch (config->GetKind_ConvNumScheme_Template()) {
      case SPACE_CENTERED : case SPACE_UPWIND :
        for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
          numerics[iMGlevel][TEMPLATE_SOL][conv_term] = new CConvective_Template(nDim, nVar_Template, config);
        break;
      default:
        SU2_MPI::Error("Convective scheme not implemented (template_solver).", CURRENT_FUNCTION);
        break;
    }

    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
      numerics[iMGlevel][TEMPLATE_SOL][visc_term] = new CViscous_Template(nDim, nVar_Template, config);

    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
      numerics[iMGlevel][TEMPLATE_SOL][source_first_term] = new CSource_Template(nDim, nVar_Template, config);

    /*--- Definition of the boundary condition method ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      numerics[iMGlevel][TEMPLATE_SOL][conv_bound_term] = new CConvective_Template(nDim, nVar_Template, config);
    }

  }

  /*--- Solver definition for the Potential, Euler, Navier-Stokes problems ---*/
  if ((euler) || (ns)) {

    /*--- Definition of the convective scheme for each equation and mesh level ---*/
    switch (config->GetKind_ConvNumScheme_Flow()) {
      case NO_CONVECTIVE :
        SU2_MPI::Error("Config file is missing the CONV_NUM_METHOD_FLOW option.", CURRENT_FUNCTION);
        break;

      case SPACE_CENTERED :
        if (compressible) {
          /*--- "conv_term" is not instantiated as all compressible centered schemes are vectorized. ---*/

          /*--- Definition of the boundary condition method ---*/
          for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwRoe_Flow(nDim, nVar_Flow, config, false);

        }
        break;
      case SPACE_UPWIND :
        if (compressible) {
          /*--- Compressible flow ---*/
          switch (config->GetKind_Upwind_Flow()) {
            case ROE:
              if (ideal_gas) {

                for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                  numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwRoe_Flow(nDim, nVar_Flow, config, roe_low_dissipation);
                  numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwRoe_Flow(nDim, nVar_Flow, config, false);
                }
              } else {

                for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                  numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwGeneralRoe_Flow(nDim, nVar_Flow, config);
                  numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwGeneralRoe_Flow(nDim, nVar_Flow, config);
                }
              }
              break;

            case AUSM:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwAUSM_Flow(nDim, nVar_Flow, config);
                numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwAUSM_Flow(nDim, nVar_Flow, config);
              }
              break;

            case AUSMPLUSUP:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwAUSMPLUSUP_Flow(nDim, nVar_Flow, config);
                numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwAUSMPLUSUP_Flow(nDim, nVar_Flow, config);
              }
              break;

            case AUSMPLUSUP2:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwAUSMPLUSUP2_Flow(nDim, nVar_Flow, config);
                numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwAUSMPLUSUP2_Flow(nDim, nVar_Flow, config);
              }
              break;

            case TURKEL:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwTurkel_Flow(nDim, nVar_Flow, config);
                numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwTurkel_Flow(nDim, nVar_Flow, config);
              }
              break;

            case L2ROE:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwL2Roe_Flow(nDim, nVar_Flow, config);
                numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwL2Roe_Flow(nDim, nVar_Flow, config);
              }
              break;
            case LMROE:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwLMRoe_Flow(nDim, nVar_Flow, config);
                numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwLMRoe_Flow(nDim, nVar_Flow, config);
              }
              break;

            case SLAU:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwSLAU_Flow(nDim, nVar_Flow, config, roe_low_dissipation);
                numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwSLAU_Flow(nDim, nVar_Flow, config, false);
              }
              break;

            case SLAU2:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwSLAU2_Flow(nDim, nVar_Flow, config, roe_low_dissipation);
                numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwSLAU2_Flow(nDim, nVar_Flow, config, false);
              }
              break;

            case HLLC:
              if (ideal_gas) {
                for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                  numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwHLLC_Flow(nDim, nVar_Flow, config);
                  numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwHLLC_Flow(nDim, nVar_Flow, config);
                }
              }
              else {
                for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                  numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwGeneralHLLC_Flow(nDim, nVar_Flow, config);
                  numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwGeneralHLLC_Flow(nDim, nVar_Flow, config);
                }
              }
              break;

            case MSW:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwMSW_Flow(nDim, nVar_Flow, config);
                numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwMSW_Flow(nDim, nVar_Flow, config);
              }
              break;

            case CUSP:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwCUSP_Flow(nDim, nVar_Flow, config);
                numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwCUSP_Flow(nDim, nVar_Flow, config);
              }
              break;

            default:
              SU2_MPI::Error("Invalid upwind scheme or not implemented.", CURRENT_FUNCTION);
              break;
          }

        }
        break;

      default:
        SU2_MPI::Error("Invalid convective scheme for the Euler / Navier-Stokes equations.", CURRENT_FUNCTION);
        break;
    }

    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {

        numerics[iMGlevel][FLOW_SOL][source_first_term] = new CSourceNothing(nDim, nVar_Flow, config);

      /*--- At the moment it is necessary to have the RHT equation in order to have a volumetric heat source. ---*/
        numerics[iMGlevel][FLOW_SOL][source_second_term] = new CSourceNothing(nDim, nVar_Flow, config);
    }

  }

  } // end "per-thread" allocation loop

}

void CDriver::Numerics_Postprocessing(CNumerics *****numerics, CSolver***, CGeometry**,
                                      CConfig *config, unsigned short val_iInst) {

  for (unsigned short iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {

    for (unsigned int iSol = 0; iSol < MAX_SOLS; iSol++) {

      for (unsigned int iTerm = 0; iTerm < MAX_TERMS*omp_get_max_threads(); iTerm++) {

        delete numerics[val_iInst][iMGlevel][iSol][iTerm];
      }
      delete [] numerics[val_iInst][iMGlevel][iSol];
    }
    delete[] numerics[val_iInst][iMGlevel];
  }
  delete[] numerics[val_iInst];

}

void CDriver::Iteration_Preprocessing(CConfig* config, CIteration *&iteration) const {

  if (rank == MASTER_NODE)
    cout << endl <<"------------------- Iteration Preprocessing ( Zone " << config->GetiZone() <<" ) ------------------" << endl;

  iteration = CIterationFactory::CreateIteration(config->GetKind_Solver(), config);

}

void CDriver::Output_Preprocessing(CConfig **config, CConfig *driver_config, COutput **&output, COutput *&driver_output){

  /*--- Definition of the output class (one for each zone). The output class
   manages the writing of all restart, volume solution, surface solution,
   surface comma-separated value, and convergence history files (both in serial
   and in parallel). ---*/

  for (iZone = 0; iZone < nZone; iZone++){

    if (rank == MASTER_NODE)
      cout << endl <<"-------------------- Output Preprocessing ( Zone " << iZone <<" ) --------------------" << endl;

    MAIN_SOLVER kindSolver = config[iZone]->GetKind_Solver();

    output[iZone] = COutputFactory::CreateOutput(kindSolver, config[iZone], nDim);

    /*--- If dry-run is used, do not open/overwrite history file. ---*/
    output[iZone]->PreprocessHistoryOutput(config[iZone], !dry_run);

    output[iZone]->PreprocessVolumeOutput(config[iZone]);

  }

  /*--- Check for an unsteady restart. Update ExtIter if necessary. ---*/
  if (config_container[ZONE_0]->GetTime_Domain() && config_container[ZONE_0]->GetRestart())
    TimeIter = config_container[ZONE_0]->GetRestart_Iter();

}


CDriver::~CDriver(void) {}

void CDriver::Print_DirectResidual(RECORDING kind_recording) {

  if (!(rank == MASTER_NODE && kind_recording == RECORDING::SOLUTION_VARIABLES)) return;

  const bool multizone = config_container[ZONE_0]->GetMultizone_Problem();

  /*--- Helper lambda func to return lenghty [iVar][iZone] string.  ---*/
  auto iVar_iZone2string = [&](unsigned short ivar, unsigned short izone) {
    if (multizone)
      return "[" + std::to_string(ivar) + "][" + std::to_string(izone) + "]";
    else
      return "[" + std::to_string(ivar) + "]";
  };

  /*--- Print residuals in the first iteration ---*/

  const unsigned short fieldWidth = 15;
  PrintingToolbox::CTablePrinter RMSTable(&std::cout);
  RMSTable.SetPrecision(config_container[ZONE_0]->GetOutput_Precision());

  /*--- The CTablePrinter requires two sweeps:
    *--- 0. Add the colum names (addVals=0=false) plus CTablePrinter.PrintHeader()
    *--- 1. Add the RMS-residual values (addVals=1=true) plus CTablePrinter.PrintFooter() ---*/
  for (int addVals = 0; addVals < 2; addVals++) {

    for (unsigned short iZone = 0; iZone < nZone; iZone++) {

      auto solvers = solver_container[iZone][INST_0][MESH_0];
      auto configs = config_container[iZone];

      /*--- Note: the FEM-Flow solvers are availalbe for disc. adjoint runs only for SingleZone. ---*/
      if (configs->GetFluidProblem()) {

        for (unsigned short iVar = 0; iVar < solvers[FLOW_SOL]->GetnVar(); iVar++) {
          if (!addVals)
            RMSTable.AddColumn("rms_Flow" + iVar_iZone2string(iVar, iZone), fieldWidth);
          else
            RMSTable << log10(solvers[FLOW_SOL]->GetRes_RMS(iVar));
        }

      }
        SU2_MPI::Error("Invalid KindSolver for SingleZone-Driver.", CURRENT_FUNCTION);
    } // loop iZone

    if (!addVals) RMSTable.PrintHeader();
    else RMSTable.PrintFooter();

  } // for addVals

  cout << "\n-------------------------------------------------------------------------\n" << endl;

}

CFluidDriver::CFluidDriver(char* confFile, unsigned short val_nZone, SU2_Comm MPICommunicator) : CDriver(confFile, val_nZone, MPICommunicator, false) {
  Max_Iter = config_container[ZONE_0]->GetnInner_Iter();
}

CFluidDriver::~CFluidDriver(void) { }

void CFluidDriver::StartSolver(){

  /*--- Main external loop of the solver. Within this loop, each iteration ---*/

  if (rank == MASTER_NODE){
    cout << endl <<"------------------------------ Begin Solver -----------------------------" << endl;
  }

  unsigned long Iter = 0;
  while ( Iter < Max_Iter ) {

    /*--- Perform some external iteration preprocessing. ---*/

    Preprocess(Iter);

    /*--- Run a single iteration of the problem (fluid, elasticity, heat, ...). ---*/

    Run();

    /*--- Update the solution for dual time stepping strategy ---*/

    Update();

    /*--- Terminate the simulation if only the Jacobian must be computed. ---*/
    if (config_container[ZONE_0]->GetJacobian_Spatial_Discretization_Only()) break;

    /*--- Monitor the computations after each iteration. ---*/

    Monitor(Iter);

    /*--- Output the solution in files. ---*/

    Output(Iter);

    /*--- If the convergence criteria has been met, terminate the simulation. ---*/

    if (StopCalc) break;

    Iter++;

  }
}


void CFluidDriver::Preprocess(unsigned long Iter) {

  /*--- Set the value of the external iteration and physical time. ---*/

  for (iZone = 0; iZone < nZone; iZone++) {
    config_container[iZone]->SetInnerIter(Iter);
    if (config_container[iZone]->GetTime_Marching() != TIME_MARCHING::STEADY)
      config_container[iZone]->SetPhysicalTime(static_cast<su2double>(Iter)*config_container[iZone]->GetDelta_UnstTimeND());
    else
      config_container[iZone]->SetPhysicalTime(0.0);
  }

  /*--- Set the initial condition for EULER/N-S/RANS and for a non FSI simulation ---*/

  if(!fsi) {
    for (iZone = 0; iZone < nZone; iZone++) {
      if (config_container[iZone]->GetFluidProblem()) {
        for (iInst = 0; iInst < nInst[iZone]; iInst++) {
          solver_container[iZone][iInst][MESH_0][FLOW_SOL]->SetInitialCondition(geometry_container[iZone][INST_0], solver_container[iZone][iInst], config_container[iZone], Iter);
        }
      }
    }
  }
}

void CFluidDriver::Run() {

  unsigned short iZone, checkConvergence;
  unsigned long IntIter, nIntIter;
  bool unsteady;

  /*--- Run a single iteration of a multi-zone problem by looping over all
   zones and executing the iterations. Note that data transers between zones
   and other intermediate procedures may be required. ---*/

  unsteady = (config_container[MESH_0]->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST) ||
             (config_container[MESH_0]->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND);

  /*--- Zone preprocessing ---*/

  for (iZone = 0; iZone < nZone; iZone++)
    iteration_container[iZone][INST_0]->Preprocess(output_container[iZone], integration_container, geometry_container, solver_container, numerics_container, config_container, iZone, INST_0);

  /*--- Begin Unsteady pseudo-time stepping internal loop, if not unsteady it does only one step --*/

  if (unsteady)
    nIntIter = config_container[MESH_0]->GetnInner_Iter();
  else
    nIntIter = 1;

  for (IntIter = 0; IntIter < nIntIter; IntIter++) {

    /*--- For each zone runs one single iteration ---*/

    for (iZone = 0; iZone < nZone; iZone++) {
      config_container[iZone]->SetInnerIter(IntIter);
      iteration_container[iZone][INST_0]->Iterate(output_container[iZone], integration_container, geometry_container, solver_container, numerics_container,
                                                  config_container, iZone, INST_0);
    }

    /*--- Check convergence in each zone --*/

    checkConvergence = 0;
    for (iZone = 0; iZone < nZone; iZone++)
    checkConvergence += (int) integration_container[iZone][INST_0][FLOW_SOL]->GetConvergence();

    /*--- If convergence was reached in every zone --*/

  if (checkConvergence == nZone) break;
  }

}

void CFluidDriver::Update() {

  for(iZone = 0; iZone < nZone; iZone++)
    iteration_container[iZone][INST_0]->Update(output_container[iZone], integration_container, geometry_container,
         solver_container, numerics_container, config_container,
         iZone, INST_0);
}

bool CFluidDriver::Monitor(unsigned long ExtIter) {

  /*--- Synchronization point after a single solver iteration. Compute the
   wall clock time required. ---*/

  StopTime = SU2_MPI::Wtime();

  IterCount++;
  UsedTime = (StopTime - StartTime) + UsedTimeCompute;

  /*--- Check if there is any change in the runtime parameters ---*/

  CConfig *runtime = nullptr;
  strcpy(runtime_file_name, "runtime.dat");
  runtime = new CConfig(runtime_file_name, config_container[ZONE_0]);
  runtime->SetTimeIter(ExtIter);
  delete runtime;

  /*--- Check whether the current simulation has reached the specified
   convergence criteria, and set StopCalc to true, if so. ---*/

  switch (config_container[ZONE_0]->GetKind_Solver()) {
    case MAIN_SOLVER::EULER: 
      StopCalc = integration_container[ZONE_0][INST_0][FLOW_SOL]->GetConvergence(); break;
    default:
      break;
  }

  /*--- Set StopCalc to true if max. number of iterations has been reached ---*/

  StopCalc = StopCalc || (ExtIter == Max_Iter - 1);

  return StopCalc;

}


void CFluidDriver::Output(unsigned long InnerIter) {

  for (iZone = 0; iZone < nZone; iZone++) {
    const auto inst = config_container[iZone]->GetiInst();

    for (iInst = 0; iInst < nInst[iZone]; ++iInst) {
      config_container[iZone]->SetiInst(iInst);
      output_container[iZone]->SetResult_Files(geometry_container[iZone][iInst][MESH_0],
                                               config_container[iZone],
                                               solver_container[iZone][iInst][MESH_0],
                                               InnerIter, StopCalc);
    }
    config_container[iZone]->SetiInst(inst);
  }

}

