﻿/*!
 * \file driver_structure.hpp
 * \brief Headers of the main subroutines for driving single or multi-zone problems.
 *        The subroutines and functions are in the <i>driver_structure.cpp</i> file.
 * \author T. Economon, H. Kline, R. Sanchez
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

#include "../../../Common/include/parallelization/mpi_structure.hpp"

#include "../integration/CIntegration.hpp"
#include "../solvers/CSolver.hpp"

#include "../../../Common/include/geometry/CGeometry.hpp"

using namespace std;

class CIteration;
class COutput;

/*!
 * \class CDriver
 * \brief Parent class for driving an iteration of a single or multi-zone problem.
 * \author T. Economon
 */
class CDriver {
protected:
  char* config_file_name;                       /*!< \brief Configuration file name of the problem.*/
  char runtime_file_name[MAX_STRING_SIZE];
  su2double StartTime,                          /*!< \brief Start point of the timer for performance benchmarking.*/
            StopTime,                           /*!< \brief Stop point of the timer for performance benchmarking.*/
            UsedTimePreproc,                    /*!< \brief Elapsed time between Start and Stop point of the timer for tracking preprocessing phase.*/
            UsedTimeCompute,                    /*!< \brief Elapsed time between Start and Stop point of the timer for tracking compute phase.*/
            UsedTimeOutput,                     /*!< \brief Elapsed time between Start and Stop point of the timer for tracking output phase.*/
            UsedTime;                           /*!< \brief Elapsed time between Start and Stop point of the timer.*/
  su2double BandwidthSum = 0.0;                 /*!< \brief Aggregate value of the bandwidth for writing restarts (to be average later).*/
  unsigned long IterCount,                      /*!< \brief Iteration count stored for performance benchmarking.*/
  OutputCount;                                  /*!< \brief Output count stored for performance benchmarking.*/
  unsigned long DOFsPerPoint;                   /*!< \brief Number of unknowns at each vertex, i.e., number of equations solved. */
  su2double Mpoints;                            /*!< \brief Total number of grid points in millions in the calculation (including ghost points).*/
  su2double MpointsDomain;                      /*!< \brief Total number of grid points in millions in the calculation (excluding ghost points).*/
  unsigned long TimeIter;                       /*!< \brief External iteration.*/
  ofstream **ConvHist_file;                     /*!< \brief Convergence history file.*/
  unsigned short iMesh,                         /*!< \brief Iterator on mesh levels.*/
                iZone,                          /*!< \brief Iterator on zones.*/
                nZone,                          /*!< \brief Total number of zones in the problem. */
                nDim,                           /*!< \brief Number of dimensions.*/
                iInst,                          /*!< \brief Iterator on instance levels.*/
                *nInst;                         /*!< \brief Total number of instances in the problem (per zone). */
  bool StopCalc;                                /*!< \brief Stop computation flag.*/
  CIteration ***iteration_container;            /*!< \brief Container vector with all the iteration methods. */
  COutput **output_container;                   /*!< \brief Pointer to the COutput class. */
  CIntegration ****integration_container;       /*!< \brief Container vector with all the integration methods. */
  CGeometry ****geometry_container;             /*!< \brief Geometrical definition of the problem. */
  CSolver *****solver_container;                /*!< \brief Container vector with all the solutions. */
  CNumerics ******numerics_container;           /*!< \brief Description of the numerical method (the way in which the equations are solved). */
  CConfig **config_container;                   /*!< \brief Definition of the particular problem. */
  CConfig *driver_config;                       /*!< \brief Definition of the driver configuration. */
  COutput *driver_output;                       /*!< \brief Definition of the driver output. */
  bool dry_run;                                 /*!< \brief Flag if SU2_CFD was started as dry-run via "SU2_CFD -d <config>.cfg" */

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] confFile - Configuration file name.
   * \param[in] val_nZone - Total number of zones.
   * \param[in] val_nDim - Number of dimensions.
   * \param[in] MPICommunicator - MPI communicator for SU2.
   */
  CDriver(char* confFile,
          unsigned short val_nZone);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CDriver(void);

  /*!
   * \brief A virtual member.
   */
  virtual void Run() { };

protected:

  /*!
   * \brief Init_Containers
   */
  void SetContainers_Null();

  /*!
   * \brief Read in the config and mesh files.
   */
  void Input_Preprocessing(CConfig **&config, CConfig *&driver_config);

  /*!
   * \brief Construction of the edge-based data structure and the multigrid structure.
   */
  void Geometrical_Preprocessing(CConfig *config, CGeometry **&geometry, bool dummy);

  /*!
   * \brief Geometrical_Preprocessing_FVM
   */
  void Geometrical_Preprocessing_FVM(CConfig *config, CGeometry **&geometry);

  /*!
   * \brief Definition of the physics iteration class or within a single zone.
   * \param[in] iteration_container - Pointer to the iteration container to be instantiated.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iZone - Index of the zone.
   */
  void Iteration_Preprocessing(CConfig *config, CIteration *&iteration) const;

  /*!
   * \brief Definition and allocation of all solution classes.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Solver_Preprocessing(CConfig *config, CGeometry **geometry, CSolver ***&solver);

  /*!
   * \brief Restart of the solvers from the restart files.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Solver_Restart(CSolver ***solver, CGeometry **geometry, CConfig *config, bool update_geo);

  /*!
   * \brief Definition and allocation of all solution classes.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Solver_Postprocessing(CSolver ****solver, CGeometry **geometry, CConfig *config, unsigned short val_iInst);

  /*!
   * \brief Definition and allocation of all integration classes.
   * \param[in] config - Definition of the particular problem.
   * \param[in] solver - Container vector with all the solutions.
   * \param[out] integration - Container vector with all the integration methods.
   */
  void Integration_Preprocessing(CConfig *config, CSolver **solver, CIntegration **&integration) const;

  /*!
   * \brief Definition and allocation of all integration classes.
   * \param[in] integration_container - Container vector with all the integration methods.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Integration_Postprocessing(CIntegration ***integration, CGeometry **geometry, CConfig *config, unsigned short val_iInst);

  /*!
   * \brief Definition and allocation of all solver classes.
   * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void Numerics_Preprocessing(CConfig *config, CGeometry **geometry, CSolver ***solver, CNumerics ****&numerics) const;

  /*!
   * \brief Definition and allocation of all solver classes.
   * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Numerics_Postprocessing(CNumerics *****numerics, CSolver ***solver, CGeometry **geometry, CConfig *config, unsigned short val_iInst);

  /*!
   * \brief Initialize Python interface functionalities
   */
  void PythonInterface_Preprocessing(CConfig** config, CGeometry**** geometry, CSolver***** solver);

  /*!
   * \brief Preprocess the output container.
   */
  void Output_Preprocessing(CConfig **config, CConfig *driver_config, COutput **&output_container, COutput *&driver_output);

  /*!
   * \brief Initiate value for static mesh movement such as the gridVel for the ROTATING frame.
   */
  void StaticMesh_Preprocessing(const CConfig *config, CGeometry **geometry);

  /*!
   * \brief A virtual member to run a Block Gauss-Seidel iteration in multizone problems.
   */
  virtual void Run_GaussSeidel(){}

  /*!
   * \brief A virtual member to run a Block-Jacobi iteration in multizone problems.
   */
  virtual void Run_Jacobi(){}

  /*!
   * \brief A virtual member.
   */
  virtual void Update() {}

  /*!
   * \brief Print out the direct residuals.
   * \param[in] kind_recording - Type of recording (full list in ENUM_RECORDING, option_structure.hpp)
   */
  void Print_DirectResidual(RECORDING kind_recording);

public:

  /*!
   * \brief Launch the computation for all zones and all physics.
   */
  virtual void StartSolver() {}

  /*!
   * \brief Deallocation routine
   */
  void Postprocessing();

  /*!
   * \brief A virtual member.
   */
  virtual void ResetConvergence();

  /*!
   * \brief Perform some pre-processing before an iteration of the physics.
   */
  virtual void Preprocess(unsigned long TimeIter){ }

  /*!
   * \brief Monitor the computation.
   */
  virtual bool Monitor(unsigned long TimeIter){ return false; }

  /*!
   * \brief Output the solution in solution file.
   */
  virtual void Output(unsigned long TimeIter){ }

  /*!
   * \brief Process the boundary conditions and update the multigrid structure.
   */
  void BoundaryConditionsUpdate();

  /*!
   * \brief Get the total drag.
   * \return Total drag.
   */
  passivedouble Get_Drag() const;

  /*!
   * \brief Get the total lift.
   * \return Total lift.
   */
  passivedouble Get_Lift() const;

  /*!
   * \brief Get the total x moment.
   * \return Total x moment.
   */
  passivedouble Get_Mx() const;

  /*!
   * \brief Get the total y moment.
   * \return Total y moment.
   */
  passivedouble Get_My() const;

  /*!
   * \brief Get the total z moment.
   * \return Total z moment.
   */
  passivedouble Get_Mz() const;

  /*!
   * \brief Get the total drag coefficient.
   * \return Total drag coefficient.
   */
  passivedouble Get_DragCoeff() const;

  /*!
   * \brief Get the total lift coefficient.
   * \return Total lift coefficient.
   */
  passivedouble Get_LiftCoeff() const;

  /*!
   * \brief Get the number of vertices (halo nodes included) from a specified marker.
   * \param[in] iMarker -  Marker identifier.
   * \return Number of vertices.
   */
  unsigned long GetNumberVertices(unsigned short iMarker) const;

  /*!
   * \brief Get the number of halo vertices from a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \return Number of vertices.
   */
  unsigned long GetNumberHaloVertices(unsigned short iMarker) const;

  /*!
   * \brief Check if a vertex is physical or not (halo node) on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return True if the specified vertex is a halo node.
   */
  bool IsAHaloNode(unsigned short iMarker, unsigned long iVertex) const;

  /*!
   * \brief Get the number of external iterations.
   * \return Number of external iterations.
   */
  unsigned long GetnTimeIter() const;

  /*!
   * \brief Get the current external iteration.
   * \return Current external iteration.
   */
  unsigned long GetTime_Iter() const;

  /*!
   * \brief Get the unsteady time step.
   * \return Unsteady time step.
   */
  passivedouble GetUnsteady_TimeStep() const;

  /*!
   * \brief Get the name of the output file for the surface.
   * \return File name for the surface output.
   */
  string GetSurfaceFileName() const;

  /*!
   * \brief Get the global index of a vertex on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return Vertex global index.
   */
  unsigned long GetVertexGlobalIndex(unsigned short iMarker, unsigned long iVertex) const;

  /*!
   * \brief Get the temperature at a vertex on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return Temperature of the vertex.
   */
  passivedouble GetVertexTemperature(unsigned short iMarker, unsigned long iVertex) const;

  /*!
   * \brief Set the temperature of a vertex on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \param[in] val_WallTemp - Value of the temperature.
   */
  void SetVertexTemperature(unsigned short iMarker, unsigned long iVertex, passivedouble val_WallTemp);

  /*!
   * \brief Get the normal (vector) at a vertex on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \param[in] unitNormal - Bool to normalise the vector.
   * \return Normal (vector) at the vertex.
   */
  vector<passivedouble> GetVertexNormal(unsigned short iMarker, unsigned long iVertex, bool unitNormal = false) const;

  /*!
   * \brief Get the unit normal (vector) at a vertex on a specified marker.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \return Unit normal (vector) at the vertex.
   */
  inline vector<passivedouble> GetVertexUnitNormal(unsigned short iMarker, unsigned long iVertex) const {
    return GetVertexNormal(iMarker, iVertex, true);
  }

  /*!
   * \brief Get all the boundary markers tags.
   * \return List of boundary markers tags.
   */
  vector<string> GetAllBoundaryMarkersTag() const;

  /*!
   * \brief Get all the boundary markers tags with their associated indices.
   * \return List of boundary markers tags with their indices.
   */
  map<string, int> GetAllBoundaryMarkers() const;

  /*!
   * \brief Get all the boundary markers tags with their associated types.
   * \return List of boundary markers tags with their types.
   */
  map<string, string> GetAllBoundaryMarkersType() const;

  /*!
   * \brief Get the sensitivity of the flow loads for the structural solver.
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   * \param[in] LoadX - Value of the load in the direction X.
   * \param[in] LoadX - Value of the load in the direction Y.
   * \param[in] LoadX - Value of the load in the direction Z.
   */
  vector<passivedouble> GetFlowLoad_Sensitivity(unsigned short iMarker, unsigned long iVertex) const;

  /*!
   * \brief Get the flow load (from the extra step - the repeated methods should be unified once the postprocessing
   * strategy is in place).
   * \param[in] iMarker - Marker identifier.
   * \param[in] iVertex - Vertex identifier.
   */
  vector<passivedouble> GetFlowLoad(unsigned short iMarker, unsigned long iVertex) const;

  /*!
   * \brief Sum the number of primal or adjoint variables for all solvers in a given zone.
   * \param[in] iZone - Index of the zone.
   * \param[in] adjoint - True to consider adjoint solvers instead of primal.
   * \return Total number of solution variables.
   */
  unsigned short GetTotalNumberOfVariables(unsigned short iZone, bool adjoint) const {
    unsigned short nVar = 0;
    for (auto iSol = 0u; iSol < MAX_SOLS; iSol++) {
      auto solver = solver_container[iZone][INST_0][MESH_0][iSol];
      if (solver && (solver->GetAdjoint() == adjoint)) nVar += solver->GetnVar();
    }
    return nVar;
  }

  /*!
   * \brief Set the solution of all solvers (adjoint or primal) in a zone.
   * \param[in] iZone - Index of the zone.
   * \param[in] adjoint - True to consider adjoint solvers instead of primal.
   * \param[in] solution - Solution object with interface (iPoint,iVar).
   * \tparam Old - If true set "old solutions" instead.
   */
  template<class Container, bool Old = false>
  void SetAllSolutions(unsigned short iZone, bool adjoint, const Container& solution) {
    const auto nPoint = geometry_container[iZone][INST_0][MESH_0]->GetnPoint();
    for (auto iSol = 0u, offset = 0u; iSol < MAX_SOLS; ++iSol) {
      auto solver = solver_container[iZone][INST_0][MESH_0][iSol];
      if (!(solver && (solver->GetAdjoint() == adjoint))) continue;
      for (auto iPoint = 0ul; iPoint < nPoint; ++iPoint)
        for (auto iVar = 0ul; iVar < solver->GetnVar(); ++iVar)
          if (!Old) solver->GetNodes()->SetSolution(iPoint, iVar, solution(iPoint,offset+iVar));
          else solver->GetNodes()->SetSolution_Old(iPoint, iVar, solution(iPoint,offset+iVar));
      offset += solver->GetnVar();
    }
  }

  /*!
   * \brief Set the "old solution" of all solvers (adjoint or primal) in a zone.
   */
  template<class Container>
  void SetAllSolutionsOld(unsigned short iZone, bool adjoint, const Container& solution) {
    SetAllSolutions<Container,true>(iZone, adjoint, solution);
  }

  /*!
   * \brief Get the solution of all solvers (adjoint or primal) in a zone.
   * \param[in] iZone - Index of the zone.
   * \param[in] adjoint - True to consider adjoint solvers instead of primal.
   * \param[out] solution - Solution object with interface (iPoint,iVar).
   */
  template<class Container>
  void GetAllSolutions(unsigned short iZone, bool adjoint, Container& solution) const {
    const auto nPoint = geometry_container[iZone][INST_0][MESH_0]->GetnPoint();
    for (auto iSol = 0u, offset = 0u; iSol < MAX_SOLS; ++iSol) {
      auto solver = solver_container[iZone][INST_0][MESH_0][iSol];
      if (!(solver && (solver->GetAdjoint() == adjoint))) continue;
      const auto& sol = solver->GetNodes()->GetSolution();
      for (auto iPoint = 0ul; iPoint < nPoint; ++iPoint)
        for (auto iVar = 0ul; iVar < solver->GetnVar(); ++iVar)
          solution(iPoint,offset+iVar) = SU2_TYPE::GetValue(sol(iPoint,iVar));
      offset += solver->GetnVar();
    }
  }

};

/*!
 * \class CFluidDriver
 * \brief Class for driving an iteration of the physics within multiple zones.
 * \author T. Economon, G. Gori
 */
class CFluidDriver : public CDriver {

protected:
   unsigned long Max_Iter;

protected:

  /*!
   * \brief Constructor of the class.
   * \param[in] confFile - Configuration file name.
   * \param[in] val_nZone - Total number of zones.
   * \param[in] val_nDim - Number of dimensions.
   * \param[in] MPICommunicator - MPI communicator for SU2.
   */
  CFluidDriver(char* confFile,
               unsigned short val_nZone);

public:
  /*!
   * \brief Destructor of the class.
   */
  ~CFluidDriver(void) override;

  /*!
   * \brief Launch the computation for all zones and all physics.
   */
  void StartSolver() override;

  /*!
   * \brief Run a single iteration of the physics within multiple zones.
   */
  void Run() override;

  /*!
   * \brief Update the dual-time solution within multiple zones.
   */
  void Update() override;

  /*!
   * \brief Output the solution in solution file.
   */
  void Output(unsigned long InnerIter) override;

  /*!
   * \brief Monitor the computation.
   */
  bool Monitor(unsigned long ExtIter) override;

  /*!
   * \brief Perform some pre-processing before an iteration of the physics.
   */
  void Preprocess(unsigned long Iter) override;

  /*!
   * \brief Transfer data among different zones (multiple zone).
   */
  void Transfer_Data(unsigned short donorZone, unsigned short targetZone);

};

