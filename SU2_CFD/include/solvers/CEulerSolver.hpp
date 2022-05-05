/*!
 * \file CEulerSolver.hpp
 * \brief Headers of the CEulerSolver class
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

#pragma once

#include "CFVMFlowSolverBase.hpp"
#include "../variables/CEulerVariable.hpp"

/*!
 * \class CEulerSolver
 * \brief Class for compressible inviscid flow problems, serves as base for Navier-Stokes/RANS.
 * \author F. Palacios
 */
class CEulerSolver : public CFVMFlowSolverBase<CEulerVariable, ENUM_REGIME::COMPRESSIBLE> {
protected:
  using BaseClass = CFVMFlowSolverBase<CEulerVariable, ENUM_REGIME::COMPRESSIBLE>;

  su2double
  Prandtl_Lam = 0.0,    /*!< \brief Laminar Prandtl number. */
  Prandtl_Turb = 0.0;   /*!< \brief Turbulent Prandtl number. */

  su2double AllBound_CEquivArea_Inv=0.0; /*!< \brief equivalent area coefficient (inviscid contribution) for all the boundaries. */
  vector<su2double> CEquivArea_Mnt;      /*!< \brief Equivalent area (inviscid contribution) for each boundary. */
  vector<su2double> CEquivArea_Inv;      /*!< \brief Equivalent area (inviscid contribution) for each boundary. */

  vector<su2double> Inflow_MassFlow;     /*!< \brief Mass flow rate for each boundary. */
  vector<su2double> Exhaust_MassFlow;    /*!< \brief Mass flow rate for each boundary. */
  vector<su2double> Inflow_Pressure;     /*!< \brief Fan face pressure for each boundary. */
  vector<su2double> Inflow_Mach;         /*!< \brief Fan face mach number for each boundary. */
  vector<su2double> Inflow_Area;         /*!< \brief Boundary total area. */
  vector<su2double> Exhaust_Area;        /*!< \brief Boundary total area. */
  vector<su2double> Exhaust_Pressure;    /*!< \brief Fan face pressure for each boundary. */
  vector<su2double> Exhaust_Temperature; /*!< \brief Fan face mach number for each boundary. */
  su2double
  Inflow_MassFlow_Total = 0.0,   /*!< \brief Mass flow rate for each boundary. */
  Exhaust_MassFlow_Total = 0.0,  /*!< \brief Mass flow rate for each boundary. */
  Inflow_Pressure_Total = 0.0,   /*!< \brief Fan face pressure for each boundary. */
  Inflow_Mach_Total = 0.0,       /*!< \brief Fan face mach number for each boundary. */
  InverseDesign = 0.0;           /*!< \brief Inverse design functional for each boundary. */
  vector<vector<unsigned long> > DonorGlobalIndex;  /*!< \brief Value of the donor global index. */
  vector<su2activematrix> DonorPrimVar;       /*!< \brief Value of the donor variables at each boundary. */
  vector<vector<su2double> > ActDisk_DeltaP;  /*!< \brief Value of the Delta P. */
  vector<vector<su2double> > ActDisk_DeltaT;  /*!< \brief Value of the Delta T. */

  su2activevector
  ActDisk_R;         /*!< \brief Value of the actuator disk Radius. */
  su2activematrix
  ActDisk_C,         /*!< \brief Value of the actuator disk Center. */
  ActDisk_Axis;      /*!< \brief Value of the actuator disk Axis. */
  vector<vector<su2double> > ActDisk_Fa; /*!< \brief Value of the actuator disk Axial Force per Unit Area. */
  vector<vector<su2double> > ActDisk_Fx; /*!< \brief Value of the actuator disk X component of the radial and tangential forces per Unit Area resultant. */
  vector<vector<su2double> > ActDisk_Fy; /*!< \brief Value of the actuator disk Y component of the radial and tangential forces per Unit Area resultant. */
  vector<vector<su2double> > ActDisk_Fz; /*!< \brief Value of the actuator disk Z component of the radial and tangential forces per Unit Area resultant. */

  su2double
  Total_CL_Prev = 0.0,        /*!< \brief Total lift coefficient for all the boundaries (fixed lift mode). */
  Total_SolidCD = 0.0,        /*!< \brief Total drag coefficient for all the boundaries. */
  Total_CD_Prev = 0.0,        /*!< \brief Total drag coefficient for all the boundaries (fixed lift mode). */
  Total_NetThrust = 0.0,      /*!< \brief Total drag coefficient for all the boundaries. */
  Total_Power = 0.0,          /*!< \brief Total drag coefficient for all the boundaries. */
  Total_ReverseFlow = 0.0,    /*!< \brief Total drag coefficient for all the boundaries. */
  Total_IDC = 0.0,            /*!< \brief Total IDC coefficient for all the boundaries. */
  Total_IDC_Mach = 0.0,       /*!< \brief Total IDC coefficient for all the boundaries. */
  Total_IDR = 0.0,            /*!< \brief Total IDC coefficient for all the boundaries. */
  Total_DC60 = 0.0,           /*!< \brief Total IDC coefficient for all the boundaries. */
  Total_MFR = 0.0,            /*!< \brief Total Mass Flow Ratio for all the boundaries. */
  Total_Prop_Eff = 0.0,       /*!< \brief Total Mass Flow Ratio for all the boundaries. */
  Total_ByPassProp_Eff = 0.0, /*!< \brief Total Mass Flow Ratio for all the boundaries. */
  Total_Adiab_Eff = 0.0,      /*!< \brief Total Mass Flow Ratio for all the boundaries. */
  Total_Poly_Eff = 0.0,       /*!< \brief Total Mass Flow Ratio for all the boundaries. */
  Total_CMx_Prev = 0.0,       /*!< \brief Total drag coefficient for all the boundaries (fixed lift mode). */
  Total_CMy_Prev = 0.0,       /*!< \brief Total drag coefficient for all the boundaries (fixed lift mode). */
  Total_CMz_Prev = 0.0,       /*!< \brief Total drag coefficient for all the boundaries (fixed lift mode). */
  Total_AeroCD = 0.0,         /*!< \brief Total aero drag coefficient for all the boundaries. */
  Total_CEquivArea = 0.0;     /*!< \brief Total Equivalent Area coefficient for all the boundaries. */

  su2double AoA_Prev,  /*!< \brief Old value of the angle of attack (monitored). */
  AoA_inc;
  bool Start_AoA_FD = false,  /*!< \brief Boolean for start of finite differencing for FixedCL mode */
  End_AoA_FD = false,         /*!< \brief Boolean for end of finite differencing for FixedCL mode */
  Update_AoA = false;         /*!< \brief Boolean to signal Angle of Attack Update */
  unsigned long Iter_Update_AoA = 0; /*!< \brief Iteration at which AoA was updated last */
  su2double dCL_dAlpha;              /*!< \brief Value of dCL_dAlpha used to control CL in fixed CL mode */
  unsigned long BCThrust_Counter;

  vector<CFluidModel*> FluidModel;   /*!< \brief fluid model used in the solver. */

  /*!
   * \brief Preprocessing actions common to the Euler and NS solvers.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] Output - boolean to determine whether to print output.
   */
  void CommonPreprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh,
                           unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output);

  /*!
   * \brief Update the AoA and freestream velocity at the farfield.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - current mesh level for the multigrid.
   * \param[in] Output - boolean to determine whether to print output.
   */
  void SetFarfield_AoA(CGeometry *geometry, CSolver **solver_container,
                       CConfig *config, unsigned short iMesh, bool Output);

  /*!
   * \brief Compute Ducros Sensor for Roe Dissipation.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void SetUpwind_Ducros_Sensor(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Compute the Fan face Mach number.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solution - Container vector with all the solutions.
   */
  void GetPower_Properties(CGeometry *geometry, CConfig *config,
                           unsigned short iMesh, bool Output);

  /*!
   * \brief Parallelization of Undivided Laplacian.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Set_MPI_ActDisk(CSolver **solver_container, CGeometry *geometry, CConfig *config);

  /*!
   * \brief Update the AoA and freestream velocity at the farfield.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - current mesh level for the multigrid.
   * \param[in] Output - boolean to determine whether to print output.
   */
  void SetActDisk_BCThrust(CGeometry *geometry, CSolver **solver_container,
                           CConfig *config, unsigned short iMesh, bool Output);

  /*!
   * \brief Read the actuator disk input file for the VARIABLE_LOAD type.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - current mesh level for the multigrid.
   * \param[in] Output - boolean to determine whether to print output.
   */
  void ReadActDisk_InputFile(CGeometry *geometry, CSolver **solver_container,
                           CConfig *config, unsigned short iMesh, bool Output);

  /*!
   * \brief Compute the max eigenvalue.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetMax_Eigenvalue(CGeometry *geometry, const CConfig *config);

  /*!
   * \brief Compute the undivided laplacian for the solution.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetUndivided_Laplacian(CGeometry *geometry, const CConfig *config);

  /*!
   * \brief Compute the dissipation sensor for centered schemes.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetCentered_Dissipation_Sensor(CGeometry *geometry, const CConfig *config);

  /*!
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition.
   * \param[in] config - Definition of the particular problem.
   */
  inline virtual void SetRoe_Dissipation(CGeometry *geometry, CConfig *config) { }

  /*!
   * \brief Compute the velocity^2, SoundSpeed, Pressure, Enthalpy, Viscosity.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \return - The number of non-physical points.
   */
  virtual unsigned long SetPrimitive_Variables(CSolver **solver_container,
                                               const CConfig *config);

  /*!
   * \brief Set gradients of coefficients for fixed CL mode
   * \param[in] config - Definition of the particular problem.
   */
  void SetCoefficient_Gradients(CConfig *config) const;

  /*!
   * \brief Instantiate a SIMD numerics object.
   * \param[in] solvers - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void InstantiateEdgeNumerics(const CSolver* const* solvers, const CConfig* config) final;

  /*!
   * \brief Set the solver nondimensionalization.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  void SetNondimensionalization(CConfig *config, unsigned short iMesh);

  /*!
   * \brief Set reference values for pressure, forces, etc.
   */
  void SetReferenceValues(const CConfig& config) final;

public:
  CEulerSolver() = delete;

  /*!
   * \brief Main constructor of this class.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Grid level.
   * \param[in] navier_stokes - True when the constructor is called by the derived class CNSSolver.
   */
  CEulerSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh, const bool navier_stokes = false);

  /*!
   * \brief Destructor of the class.
   */
  ~CEulerSolver(void) override;

  /*!
   * \brief Compute the pressure at the infinity.
   * \return Value of the pressure at the infinity.
   */
  inline CFluidModel* GetFluidModel(void) const final { return FluidModel[omp_get_thread_num()]; }

  /*!
   * \brief Compute the time step for solving the Euler equations.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] Iteration - Value of the current iteration.
   */
  void SetTime_Step(CGeometry *geometry,
                    CSolver **solver_container,
                    CConfig *config,
                    unsigned short iMesh,
                    unsigned long Iteration) final;

  /*!
   * \brief Compute the spatial integration using a centered scheme.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics_container - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   */
  void Centered_Residual(CGeometry *geometry,
                         CSolver **solver_container,
                         CNumerics **numerics_container,
                         CConfig *config,
                         unsigned short iMesh,
                         unsigned short iRKStep) final;

  /*!
   * \brief Compute the spatial integration using a upwind scheme.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics_container - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  void Upwind_Residual(CGeometry *geometry,
                       CSolver **solver_container,
                       CNumerics **numerics_container,
                       CConfig *config,
                       unsigned short iMesh) final;

  /*!
   * \brief Recompute the extrapolated quantities, after MUSCL reconstruction,
   *        in a more thermodynamically consistent way.
   * \note This method is static to improve the chances of it being used in a
   *       thread-safe manner.
   * \param[in,out] fluidModel - The fluid model.
   * \param[in] nDim - Number of physical dimensions.
   * \param[in,out] primitive - Primitive variables.
   * \param[out] secondary - Secondary variables.
   */
  static void ComputeConsistentExtrapolation(CFluidModel *fluidModel,
                                             unsigned short nDim,
                                             su2double *primitive,
                                             su2double *secondary);

  /*!
   * \brief Apply low Mach number correction to the primitives at two points,
   *        usually connected by an edge.
   * \note This method is static to improve the chances of it being used in a
   *       thread-safe manner.
   * \param[in,out] fluidModel - The fluid model.
   * \param[in] nDim - Number of physical dimensions.
   * \param[in,out] primitive_i - Primitive variables at point i.
   * \param[in,out] primitive_j - Primitive variables at point j.
   */
  static void LowMachPrimitiveCorrection(CFluidModel *fluidModel,
                                         unsigned short nDim,
                                         su2double *primitive_i,
                                         su2double *primitive_j);

  /*!
   * \brief Source term integration.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics_container - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  void Source_Residual(CGeometry *geometry,
                       CSolver **solver_container,
                       CNumerics **numerics_container,
                       CConfig *config,
                       unsigned short iMesh) override;

  /*!
   * \brief Source term integration.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  void Source_Template(CGeometry *geometry,
                       CSolver **solver_container,
                       CNumerics *numerics,
                       CConfig *config,
                       unsigned short iMesh) final;

  /*!
   * \brief Compute primitive variables and their gradients.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] Output - boolean to determine whether to print output.
   */
  void Preprocessing(CGeometry *geometry,
                     CSolver **solver_container,
                     CConfig *config,
                     unsigned short iMesh,
                     unsigned short iRKStep,
                     unsigned short RunTime_EqSystem,
                     bool Output) override;

  /*!
   * \brief Compute the preconditioner for convergence acceleration by Roe-Turkel method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iPoint - Index of the grid point.
   * \param[in] delta - Volume over delta t.
   * \param[in,out] preconditioner - The preconditioner matrix, must be allocated outside.
   */
  void SetPreconditioner(const CConfig *config, unsigned long iPoint,
                         su2double delta, su2activematrix& preconditioner) const;

  /*!
   * \author H. Kline
   * \brief Compute weighted-sum "combo" objective output
   * \param[in] config - Definition of the particular problem.
   * \param[in] solver - Container vector with all the solutions.
   */
  void Evaluate_ObjFunc(const CConfig *config, CSolver **solver) override;

  /*!
   * \brief Impose the far-field boundary condition using characteristics.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Far_Field(CGeometry *geometry,
                    CSolver **solver_container,
                    CNumerics *conv_numerics,
                    CNumerics *visc_numerics,
                    CConfig *config,
                    unsigned short val_marker) final;

  /*!
   * \brief Impose the engine inflow boundary condition.
    * \param[in] geometry - Geometrical definition of the problem.
    * \param[in] solver_container - Container vector with all the solutions.
    * \param[in] conv_numerics - Description of the numerical method.
    * \param[in] visc_numerics - Description of the numerical method.
    * \param[in] config - Definition of the particular problem.
    * \param[in] val_marker - Surface marker where the boundary condition is applied.
    */
  void BC_ActDisk_Inlet(CGeometry *geometry,
                        CSolver **solver_container,
                        CNumerics *conv_numerics,
                        CNumerics *visc_numerics,
                        CConfig *config,
                        unsigned short val_marker) final;

  /*!
   * \brief Impose the engine exhaust boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_ActDisk_Outlet(CGeometry *geometry,
                         CSolver **solver_container,
                         CNumerics *conv_numerics,
                         CNumerics *visc_numerics,
                         CConfig *config,
                         unsigned short val_marker) final;

  /*!
   * \brief Impose an actuator disk inlet boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   * \param[in] val_inlet_surface - Boolean for whether val_marker is an inlet
   */
  void BC_ActDisk(CGeometry *geometry,
                  CSolver **solver_container,
                  CNumerics *conv_numerics,
                  CNumerics *visc_numerics,
                  CConfig *config,
                  unsigned short val_marker,
                  bool val_inlet_surface) final;

  /*!
   * \brief Impose an actuator disk with variable load boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   * \param[in] val_inlet_surface - Boolean for whether val_marker is an inlet
   */
  void BC_ActDisk_VariableLoad(CGeometry *geometry,
                               CSolver **solver_container,
                               CNumerics *conv_numerics,
                               CNumerics *visc_numerics,
                               CConfig *config,
                               unsigned short val_marker,
                               bool val_inlet_surface);

  /*!
   * \author: G.Gori, S.Vitale, M.Pini, A.Guardone, P.Colonna
   *
   * \brief Impose the boundary condition using characteristic recostruction.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Riemann(CGeometry *geometry,
                  CSolver **solver_container,
                  CNumerics *conv_numerics,
                  CNumerics *visc_numerics,
                  CConfig *config,
                  unsigned short val_marker) final;

  /*!
   * \brief Impose a subsonic inlet boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Inlet(CGeometry *geometry,
                CSolver **solver_container,
                CNumerics *conv_numerics,
                CNumerics *visc_numerics,
                CConfig *config,
                unsigned short val_marker) final;

  /*!
   * \brief Impose a supersonic inlet boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Supersonic_Inlet(CGeometry *geometry,
                           CSolver **solver_container,
                           CNumerics *conv_numerics,
                           CNumerics *visc_numerics,
                           CConfig *config,
                           unsigned short val_marker) final;

  /*!
   * \brief Impose a supersonic outlet boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Supersonic_Outlet(CGeometry *geometry,
                            CSolver **solver_container,
                            CNumerics *conv_numerics,
                            CNumerics *visc_numerics,
                            CConfig *config,
                            unsigned short val_marker) final;

  /*!
   * \brief Impose the outlet boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Outlet(CGeometry *geometry,
                 CSolver **solver_container,
                 CNumerics *conv_numerics,
                 CNumerics *visc_numerics,
                 CConfig *config,
                 unsigned short val_marker) final;

  /*!
   * \brief Impose the nacelle inflow boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Engine_Inflow(CGeometry *geometry,
                        CSolver **solver_container,
                        CNumerics *conv_numerics,
                        CNumerics *visc_numerics,
                        CConfig *config,
                        unsigned short val_marker) final;

  /*!
   * \brief Impose the ancelle exhaust boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Engine_Exhaust(CGeometry *geometry,
                         CSolver **solver_container,
                         CNumerics *conv_numerics,
                         CNumerics *visc_numerics,
                         CConfig *config,
                         unsigned short val_marker) final;

  /*!
   * \brief Update the solution using a Runge-Kutta scheme.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   */
  void ExplicitRK_Iteration(CGeometry *geometry,
                            CSolver **solver_container,
                            CConfig *config,
                            unsigned short iRKStep) final;

  /*!
   * \brief Update the solution using the classical fourth-order Runge-Kutta scheme.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   */
  void ClassicalRK4_Iteration(CGeometry *geometry,
                              CSolver **solver_container,
                              CConfig *config,
                              unsigned short iRKStep) final;

  /*!
   * \brief Check for convergence of the Fixed CL mode to the target CL
   * \param[in] config - Definition of the particular problem.
   * \param[in] convergence - boolean for whether the solution is converged
   * \return boolean for whether the Fixed CL mode is converged to target CL
   */
  bool FixedCL_Convergence(CConfig *config, bool convergence) final;

  /*!
   * \brief Checking whether fixed CL mode in finite-differencing mode
   * \return boolean for whether the Fixed CL mode is currently in finite-differencing mode
   */
  inline bool GetStart_AoA_FD(void) const final { return Start_AoA_FD; }

  /*!
   * \brief Checking whether fixed CL mode in finite-differencing mode
   * \return boolean for whether the Fixed CL mode is currently in finite-differencing mode
   */
  inline bool GetEnd_AoA_FD(void) const final { return End_AoA_FD; }

    /*!
   * \brief Get the iteration of the last AoA update (Fixed CL Mode)
   * \return value for the last iteration that the AoA was updated
   */
  inline unsigned long GetIter_Update_AoA(void) const final { return Iter_Update_AoA; }

  /*!
   * \brief Get the AoA before the most recent update
   * \return value of the AoA before most recent update
   */
  inline su2double GetPrevious_AoA(void) const final { return AoA_Prev; }

  /*!
   * \brief Get the CL Driver's control command
   * \return value of CL Driver control command (AoA_inc)
   */
  inline su2double GetAoA_inc(void) const final { return AoA_inc; }

  /*!
   * \brief Update the solution using the explicit Euler scheme.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void ExplicitEuler_Iteration(CGeometry *geometry,
                               CSolver **solver_container,
                               CConfig *config) final;

  /*!
   * \brief Prepare an implicit iteration.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void PrepareImplicitIteration(CGeometry *geometry, CSolver**, CConfig *config) final;

  /*!
   * \brief Complete an implicit iteration.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void CompleteImplicitIteration(CGeometry *geometry, CSolver**, CConfig *config) final;

  /*!
   * \brief Provide the mass flow rate.
   * \param val_marker Surface where the coeficient is going to be computed.
   * \return Value of the mass flow rate on the surface <i>val_marker</i>.
   */
  inline su2double GetInflow_MassFlow(unsigned short val_marker) const final { return Inflow_MassFlow[val_marker]; }

  /*!
   * \brief Provide the mass flow rate.
   * \param val_marker Surface where the coeficient is going to be computed.
   * \return Value of the mass flow rate on the surface <i>val_marker</i>.
   */
  inline su2double GetExhaust_MassFlow(unsigned short val_marker) const final { return Exhaust_MassFlow[val_marker]; }

  /*!
   * \brief Provide the mass flow rate.
   * \param val_marker Surface where the coeficient is going to be computed.
   * \return Value of the fan face pressure on the surface <i>val_marker</i>.
   */
  inline su2double GetInflow_Pressure(unsigned short val_marker) const final { return Inflow_Pressure[val_marker]; }

  /*!
   * \brief Provide the mass flow rate.
   * \param val_marker Surface where the coeficient is going to be computed.
   * \return Value of the fan face mach on the surface <i>val_marker</i>.
   */
  inline su2double GetInflow_Mach(unsigned short val_marker) const final { return Inflow_Mach[val_marker]; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional aero CD.
   * \return Value of the Aero CD coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_AeroCD() const final { return Total_AeroCD; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional Near-Field pressure coefficient.
   * \return Value of the NearField pressure coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_CNearFieldOF() const final { return Total_CNearFieldOF; }

  /*!
   * \brief Set the value of the Aero drag.
   * \param[in] val_cequivarea - Value of the aero drag.
   */
  inline void SetTotal_AeroCD(su2double val_aerocd) final { Total_AeroCD = val_aerocd; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_NetThrust() const final { return Total_NetThrust; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_Power() const final { return Total_Power; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_SolidCD() const final { return Total_SolidCD; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_ReverseFlow() const final { return Total_ReverseFlow; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_MFR() const final { return Total_MFR; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_Prop_Eff() const final { return Total_Prop_Eff; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_ByPassProp_Eff() const final { return Total_ByPassProp_Eff; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_Adiab_Eff() const final { return Total_Adiab_Eff; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_Poly_Eff() const final { return Total_Poly_Eff; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_IDC() const final { return Total_IDC; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_IDC_Mach() const final { return Total_IDC_Mach; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_IDR() const final { return Total_IDR; }

  /*!
   * \brief Provide the total (inviscid + viscous) non dimensional drag coefficient.
   * \return Value of the drag coefficient (inviscid + viscous contribution).
   */
  inline su2double GetTotal_DC60() const final { return Total_DC60; }

  /*!
   * \brief Store the total (inviscid + viscous) non dimensional drag coefficient.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  inline void SetTotal_NetThrust(su2double val_Total_NetThrust) final { Total_NetThrust = val_Total_NetThrust; }

  /*!
   * \brief Store the total (inviscid + viscous) non dimensional drag coefficient.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  inline void SetTotal_Power(su2double val_Total_Power) final { Total_Power = val_Total_Power; }

  /*!
   * \brief Store the total (inviscid + viscous) non dimensional drag coefficient.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  inline void SetTotal_SolidCD(su2double val_Total_SolidCD) final { Total_SolidCD = val_Total_SolidCD; }

  /*!
   * \brief Store the total (inviscid + viscous) non dimensional drag coefficient.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  inline void SetTotal_ReverseFlow(su2double val_Total_ReverseFlow) final { Total_ReverseFlow = val_Total_ReverseFlow; }

  /*!
   * \brief Store the total (inviscid + viscous) non dimensional drag coefficient.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  inline void SetTotal_MFR(su2double val_Total_MFR) final { Total_MFR = val_Total_MFR; }

  /*!
   * \brief Store the total (inviscid + viscous) non dimensional drag coefficient.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  inline void SetTotal_Prop_Eff(su2double val_Total_Prop_Eff) final { Total_Prop_Eff = val_Total_Prop_Eff; }

  /*!
   * \brief Store the total (inviscid + viscous) non dimensional drag coefficient.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  inline void SetTotal_ByPassProp_Eff(su2double val_Total_ByPassProp_Eff) final { Total_ByPassProp_Eff = val_Total_ByPassProp_Eff; }

  /*!
   * \brief Store the total (inviscid + viscous) non dimensional drag coefficient.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  inline void SetTotal_Adiab_Eff(su2double val_Total_Adiab_Eff) final { Total_Adiab_Eff = val_Total_Adiab_Eff; }

  /*!
   * \brief Store the total (inviscid + viscous) non dimensional drag coefficient.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  inline void SetTotal_Poly_Eff(su2double val_Total_Poly_Eff) final { Total_Poly_Eff = val_Total_Poly_Eff; }

  /*!
   * \brief Store the total (inviscid + viscous) non dimensional drag coefficient.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  inline void SetTotal_IDC(su2double val_Total_IDC) final { Total_IDC = val_Total_IDC; }

  /*!
   * \brief Store the total (inviscid + viscous) non dimensional drag coefficient.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  inline void SetTotal_IDC_Mach(su2double val_Total_IDC_Mach) final { Total_IDC_Mach = val_Total_IDC_Mach; }

  /*!
   * \brief Store the total (inviscid + viscous) non dimensional drag coefficient.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  inline void SetTotal_IDR(su2double val_Total_IDR) final { Total_IDR = val_Total_IDR; }

  /*!
   * \brief Store the total (inviscid + viscous) non dimensional drag coefficient.
   * \param[in] val_Total_CD - Value of the total drag coefficient.
   */
  inline void SetTotal_DC60(su2double val_Total_DC60) final { Total_DC60 = val_Total_DC60; }

  /*!
   * \brief Value of the characteristic global index at the boundaries.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  inline unsigned long GetDonorGlobalIndex(unsigned short val_marker,
                                           unsigned long val_vertex) const final {
    return DonorGlobalIndex[val_marker][val_vertex];
  }

  /*!
   * \brief Value of the characteristic global index at the boundaries.
   * \param[in] val_marker - Surface marker where the coefficient is computed.
   * \param[in] val_vertex - Vertex of the marker <i>val_marker</i> where the coefficient is evaluated.
   * \return Value of the pressure coefficient.
   */
  inline void SetDonorGlobalIndex(unsigned short val_marker,
                                  unsigned long val_vertex,
                                  unsigned long val_index) final {
    DonorGlobalIndex[val_marker][val_vertex] = val_index;
  }

  /*!
   * \brief Update the multi-grid structure for the customized boundary conditions
   * \param geometry_container - Geometrical definition.
   * \param config - Definition of the particular problem.
   */
  void UpdateCustomBoundaryConditions(CGeometry **geometry_container, CConfig *config) final;

  /*!
   * \brief Set the initial condition for the Euler Equations.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] ExtIter - External iteration.
   */
  void SetInitialCondition(CGeometry **geometry,
                           CSolver ***solver_container,
                           CConfig *config,
                           unsigned long TimeIter) final;

  /*!
   * \brief Set the solution using the Freestream values.
   * \param[in] config - Definition of the particular problem.
   */
  void SetFreeStream_Solution(const CConfig *config) final;

  /*!
   * \brief Print verification error to screen.
   * \param[in] config - Definition of the particular problem.
   */
  void PrintVerificationError(const CConfig* config) const final;

  /*!
   * \brief The Euler and NS solvers support MPI+OpenMP.
   */
  inline bool GetHasHybridParallel() const final { return true; }

};
