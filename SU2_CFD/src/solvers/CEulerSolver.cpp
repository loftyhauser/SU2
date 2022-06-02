/*!
 * \file CEulerSolver.cpp
 * \brief Main subrotuines for solving Finite-Volume Euler flow problems.
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

#include "../../include/solvers/CEulerSolver.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"
#include "../../../Common/include/toolboxes/printing_toolbox.hpp"
#include "../../include/fluid/CIdealGas.hpp"
#include "../../include/numerics_simd/CNumericsSIMD.hpp"
#include "../../include/solvers/CFVMFlowSolverBase.inl"

CEulerSolver::CEulerSolver(CGeometry *geometry, CConfig *config,
                           unsigned short iMesh, const bool navier_stokes) :
  CFVMFlowSolverBase<CEulerVariable, ENUM_REGIME::COMPRESSIBLE>(*geometry, *config) {

  /*--- Based on the navier_stokes boolean, determine if this constructor is
   *    being called by itself, or by its derived class CNSSolver. ---*/
  string description;
  unsigned short nSecVar;
  if (navier_stokes) {
    description = "Navier-Stokes";
    nSecVar = 8;
  }
  else {
    description = "Euler";
    nSecVar = 2;
  }

  const auto nZone = geometry->GetnZone();
  const bool restart = (config->GetRestart() || config->GetRestart_Flow());
  const bool dual_time = (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST) ||
                         (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND);
  const bool time_stepping = (config->GetTime_Marching() == TIME_MARCHING::TIME_STEPPING);

  int Unst_RestartIter = 0;
  unsigned long iPoint, iMarker, counter_local = 0, counter_global = 0;
  unsigned short iDim, nLineLets;
  su2double StaticEnergy, Density, Velocity2, Pressure, Temperature;

  /*--- Store the multigrid level. ---*/
  MGLevel = iMesh;

  /*--- Check for a restart file to evaluate if there is a change in the angle of attack
   before computing all the non-dimesional quantities. ---*/

  if (!(!restart || (iMesh != MESH_0) || nZone > 1) && config->GetFixed_CL_Mode()) {

    /*--- Modify file name for a dual-time unsteady restart ---*/

    if (dual_time) {
      if (adjoint) Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_AdjointIter())-1;
      else if (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST)
        Unst_RestartIter = SU2_TYPE::Int(config->GetRestart_Iter())-1;
      else Unst_RestartIter = SU2_TYPE::Int(config->GetRestart_Iter())-2;
    }

    /*--- Modify file name for a time stepping unsteady restart ---*/

    if (time_stepping) {
      if (adjoint) Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_AdjointIter())-1;
      else Unst_RestartIter = SU2_TYPE::Int(config->GetRestart_Iter())-1;
    }

    /*--- Read and store the restart metadata. ---*/

    string filename_ = "flow";
    filename_ = config->GetFilename(filename_, ".meta", Unst_RestartIter);
    Read_SU2_Restart_Metadata(geometry, config, adjoint, filename_);

  }

  /*--- Set the gamma value ---*/

  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;

  /*--- Define geometry constants in the solver structure
   Compressible flow, primitive variables (T, vx, vy, vz, P, rho, h, c, lamMu, EddyMu, ThCond, Cp).
   ---*/

  nDim = geometry->GetnDim();

  nVar = nDim+2;
  nPrimVar = nDim+9; nPrimVarGrad = nDim+4;
  nSecondaryVar = nSecVar; nSecondaryVarGrad = 2;

  /*--- Initialize nVarGrad for deallocation ---*/

  nVarGrad = nPrimVarGrad;

  nMarker      = config->GetnMarker_All();
  nPoint       = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();

  /*--- Store the number of vertices on each marker for deallocation later ---*/

  nVertex.resize(nMarker);
  for (iMarker = 0; iMarker < nMarker; iMarker++)
    nVertex[iMarker] = geometry->nVertex[iMarker];

  /*--- Perform the non-dimensionalization for the flow equations using the
   specified reference values. ---*/

  SetNondimensionalization(config, iMesh);

  /*--- Allocate base class members. ---*/

  Allocate(*config);

  /*--- MPI + OpenMP initialization. ---*/

  HybridParallelInitialization(*config, *geometry);

  /*--- Jacobians and vector structures for implicit computations ---*/

  if (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT) {

    if (rank == MASTER_NODE)
      cout << "Initialize Jacobian structure (" << description << "). MG level: " << iMesh <<"." << endl;

    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config, ReducerStrategy);

    if (config->GetKind_Linear_Solver_Prec() == LINELET) {
      nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
      if (rank == MASTER_NODE)
        cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
    }
  }
  else {
    if (rank == MASTER_NODE)
      cout << "Explicit scheme. No Jacobian structure (" << description << "). MG level: " << iMesh <<"." << endl;
  }

  /*--- Store the value of the primitive variables + 2 turb variables at the boundaries,
   used for IO with a donor cell ---*/

  AllocVectorOfMatrices(nVertex, nPrimVar, DonorPrimVar);

  /*--- Store the value of the characteristic primitive variables index at the boundaries ---*/

  DonorGlobalIndex.resize(nMarker);
  for (iMarker = 0; iMarker < nMarker; iMarker++)
    DonorGlobalIndex[iMarker].resize(nVertex[iMarker],0);

  /*--- Supersonic coefficients ---*/

  CEquivArea_Inv.resize(nMarker);

  /*--- Engine simulation ---*/

  Inflow_MassFlow.resize(nMarker);
  Inflow_Pressure.resize(nMarker);
  Inflow_Mach.resize(nMarker);
  Inflow_Area.resize(nMarker);

  Exhaust_Temperature.resize(nMarker);
  Exhaust_MassFlow.resize(nMarker);
  Exhaust_Pressure.resize(nMarker);
  Exhaust_Area.resize(nMarker);

  /*--- Read farfield conditions from config ---*/

  Temperature_Inf = config->GetTemperature_FreeStreamND();
  Velocity_Inf = config->GetVelocity_FreeStreamND();
  Pressure_Inf = config->GetPressure_FreeStreamND();
  Density_Inf = config->GetDensity_FreeStreamND();
  Energy_Inf = config->GetEnergy_FreeStreamND();
  Mach_Inf = config->GetMach();

  SetReferenceValues(*config);

  /*--- Initialize fan face pressure, fan face mach number, and mass flow rate ---*/

  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    Inflow_MassFlow[iMarker]     = 0.0;
    Inflow_Mach[iMarker]         = Mach_Inf;
    Inflow_Pressure[iMarker]     = Pressure_Inf;
    Inflow_Area[iMarker]         = 0.0;

    Exhaust_MassFlow[iMarker]    = 0.0;
    Exhaust_Temperature[iMarker] = Temperature_Inf;
    Exhaust_Pressure[iMarker]    = Pressure_Inf;
    Exhaust_Area[iMarker]        = 0.0;
  }

  /*--- Initialize the solution to the far-field state everywhere. ---*/

  nodes = new CEulerVariable(Density_Inf, Velocity_Inf, Energy_Inf, nPoint, nDim, nVar, config);
  SetBaseClassPointerToNodes();

  /*--- Check that the initial solution is physical, report any non-physical nodes ---*/

  counter_local = 0;

  for (iPoint = 0; iPoint < nPoint; iPoint++) {

    Density = nodes->GetDensity(iPoint);

    Velocity2 = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      Velocity2 += pow(nodes->GetSolution(iPoint,iDim+1)/Density,2);

    StaticEnergy= nodes->GetEnergy(iPoint) - 0.5*Velocity2;

    GetFluidModel()->SetTDState_rhoe(Density, StaticEnergy);
    Pressure= GetFluidModel()->GetPressure();
    Temperature= GetFluidModel()->GetTemperature();

    /*--- Use the values at the infinity ---*/

    su2double Solution[MAXNVAR] = {0.0};
    if ((Pressure < 0.0) || (Density < 0.0) || (Temperature < 0.0)) {
      Solution[0] = Density_Inf;
      for (iDim = 0; iDim < nDim; iDim++)
        Solution[iDim+1] = Velocity_Inf[iDim]*Density_Inf;
      Solution[nDim+1] = Energy_Inf*Density_Inf;
      nodes->SetSolution(iPoint,Solution);
      nodes->SetSolution_Old(iPoint,Solution);
      counter_local++;
    }

  }

  /*--- Warning message about non-physical points ---*/

  if (config->GetComm_Level() == COMM_FULL) {

    SU2_MPI::Reduce(&counter_local, &counter_global, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, SU2_MPI::GetComm());

    if ((rank == MASTER_NODE) && (counter_global != 0))
      cout << "Warning. The original solution contains " << counter_global << " points that are not physical." << endl;
  }

  /*--- Initial comms. ---*/

  CommunicateInitialState(geometry, config);

  /*--- Add the solver name (max 8 characters). ---*/
  SolverName = "C.FLOW";

  /*--- Finally, check that the static arrays will be large enough (keep this
   *    check at the bottom to make sure we consider the "final" values). ---*/
  if((nDim > MAXNDIM) || (nPrimVar > MAXNVAR) || (nSecondaryVar > MAXNVAR))
    SU2_MPI::Error("Oops! The CEulerSolver static array sizes are not large enough.",CURRENT_FUNCTION);
}

CEulerSolver::~CEulerSolver(void) {

  for(auto& model : FluidModel) delete model;
}

void CEulerSolver::InstantiateEdgeNumerics(const CSolver* const* solver_container, const CConfig* config) {

  SU2_OMP_BARRIER
  SU2_OMP_MASTER {

  if (config->Low_Mach_Correction())
    SU2_MPI::Error("Low-Mach correction is not supported with vectorization.", CURRENT_FUNCTION);

  edgeNumerics = CNumericsSIMD::CreateNumerics(*config, nDim, MGLevel);

  if (!edgeNumerics)
    SU2_MPI::Error("The numerical scheme + gas model in use do not "
                   "support vectorization.", CURRENT_FUNCTION);

  }
  END_SU2_OMP_MASTER
  SU2_OMP_BARRIER
}

void CEulerSolver::SetNondimensionalization(CConfig *config, unsigned short iMesh) {

  su2double Temperature_FreeStream = 0.0, Mach2Vel_FreeStream = 0.0, ModVel_FreeStream = 0.0,
  Energy_FreeStream = 0.0, ModVel_FreeStreamND = 0.0, 
  Density_FreeStream = 0.0, Pressure_FreeStream = 0.0, 
  Length_Ref = 0.0, Density_Ref = 0.0, Pressure_Ref = 0.0, Velocity_Ref = 0.0,
  Temperature_Ref = 0.0, Time_Ref = 0.0, Omega_Ref = 0.0, Force_Ref = 0.0,
  Gas_Constant_Ref = 0.0, Viscosity_Ref = 0.0, Conductivity_Ref = 0.0, Energy_Ref= 0.0,
  Froude = 0.0, Pressure_FreeStreamND = 0.0, Density_FreeStreamND = 0.0,
  Temperature_FreeStreamND = 0.0, Gas_ConstantND = 0.0,
  Velocity_FreeStreamND[3] = {0.0, 0.0, 0.0},
  Energy_FreeStreamND = 0.0,
  Total_UnstTimeND = 0.0, Delta_UnstTimeND = 0.0, Heat_Flux_Ref = 0.0;

  unsigned short iDim;

  /*--- Local variables ---*/

  su2double Alpha         = config->GetAoA()*PI_NUMBER/180.0;
  su2double Beta          = config->GetAoS()*PI_NUMBER/180.0;
  su2double Mach          = config->GetMach();
  bool unsteady           = (config->GetTime_Marching() != TIME_MARCHING::STEADY);
  bool gravity            = config->GetGravityForce();
  bool free_stream_temp   = (config->GetKind_FreeStreamOption() == FREESTREAM_OPTION::TEMPERATURE_FS);

  /*--- Compressible non dimensionalization ---*/

  /*--- Compute the Free Stream velocity, using the Mach number ---*/

  Pressure_FreeStream = config->GetPressure_FreeStream();
  Density_FreeStream  = config->GetDensity_FreeStream();
  Temperature_FreeStream = config->GetTemperature_FreeStream();

  CFluidModel* auxFluidModel = nullptr;

  switch (config->GetKind_FluidModel()) {

    case STANDARD_AIR:

      switch (config->GetSystemMeasurements()) {
        case SI: config->SetGas_Constant(287.058); break;
        case US: config->SetGas_Constant(1716.49); break;
      }

      auxFluidModel = new CIdealGas(1.4, config->GetGas_Constant());

      break;

    case IDEAL_GAS:

      auxFluidModel = new CIdealGas(Gamma, config->GetGas_Constant());
      break;

    default:
      SU2_MPI::Error("Unknown fluid model.", CURRENT_FUNCTION);
      break;
  }

  if (free_stream_temp) {
    auxFluidModel->SetTDState_PT(Pressure_FreeStream, Temperature_FreeStream);
    Density_FreeStream = auxFluidModel->GetDensity();
    config->SetDensity_FreeStream(Density_FreeStream);
  }
  else {
    auxFluidModel->SetTDState_Prho(Pressure_FreeStream, Density_FreeStream );
    Temperature_FreeStream = auxFluidModel->GetTemperature();
    config->SetTemperature_FreeStream(Temperature_FreeStream);
  }

  Mach2Vel_FreeStream = auxFluidModel->GetSoundSpeed();

  /*--- Compute the Free Stream velocity, using the Mach number ---*/

  if (nDim == 2) {
    config->GetVelocity_FreeStream()[0] = cos(Alpha)*Mach*Mach2Vel_FreeStream;
    config->GetVelocity_FreeStream()[1] = sin(Alpha)*Mach*Mach2Vel_FreeStream;
  }
  if (nDim == 3) {
    config->GetVelocity_FreeStream()[0] = cos(Alpha)*cos(Beta)*Mach*Mach2Vel_FreeStream;
    config->GetVelocity_FreeStream()[1] = sin(Beta)*Mach*Mach2Vel_FreeStream;
    config->GetVelocity_FreeStream()[2] = sin(Alpha)*cos(Beta)*Mach*Mach2Vel_FreeStream;
  }

  /*--- Compute the modulus of the free stream velocity ---*/

  ModVel_FreeStream = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    ModVel_FreeStream += config->GetVelocity_FreeStream()[iDim]*config->GetVelocity_FreeStream()[iDim];
  ModVel_FreeStream = sqrt(ModVel_FreeStream); config->SetModVel_FreeStream(ModVel_FreeStream);



    /*--- For inviscid flow, energy is calculated from the specified
       FreeStream quantities using the proper gas law. ---*/

  Energy_FreeStream = auxFluidModel->GetStaticEnergy() + 0.5*ModVel_FreeStream*ModVel_FreeStream;

  /*-- Compute the freestream energy. ---*/

  config->SetEnergy_FreeStream(Energy_FreeStream);

  /*--- Compute non dimensional quantities. By definition,
     Lref is one because we have converted the grid to meters. ---*/

  if (config->GetRef_NonDim() == DIMENSIONAL) {
    Pressure_Ref      = 1.0;
    Density_Ref       = 1.0;
    Temperature_Ref   = 1.0;
  }
  else if (config->GetRef_NonDim() == FREESTREAM_PRESS_EQ_ONE) {
    Pressure_Ref      = Pressure_FreeStream;     // Pressure_FreeStream = 1.0
    Density_Ref       = Density_FreeStream;      // Density_FreeStream = 1.0
    Temperature_Ref   = Temperature_FreeStream;  // Temperature_FreeStream = 1.0
  }
  else if (config->GetRef_NonDim() == FREESTREAM_VEL_EQ_MACH) {
    Pressure_Ref      = Gamma*Pressure_FreeStream; // Pressure_FreeStream = 1.0/Gamma
    Density_Ref       = Density_FreeStream;        // Density_FreeStream = 1.0
    Temperature_Ref   = Temperature_FreeStream;    // Temp_FreeStream = 1.0
  }
  else if (config->GetRef_NonDim() == FREESTREAM_VEL_EQ_ONE) {
    Pressure_Ref      = Mach*Mach*Gamma*Pressure_FreeStream; // Pressure_FreeStream = 1.0/(Gamma*(M_inf)^2)
    Density_Ref       = Density_FreeStream;        // Density_FreeStream = 1.0
    Temperature_Ref   = Temperature_FreeStream;    // Temp_FreeStream = 1.0
  }
  config->SetPressure_Ref(Pressure_Ref);
  config->SetDensity_Ref(Density_Ref);
  config->SetTemperature_Ref(Temperature_Ref);

  Length_Ref        = 1.0;                                                         config->SetLength_Ref(Length_Ref);
  Velocity_Ref      = sqrt(config->GetPressure_Ref()/config->GetDensity_Ref());    config->SetVelocity_Ref(Velocity_Ref);
  Time_Ref          = Length_Ref/Velocity_Ref;                                     config->SetTime_Ref(Time_Ref);
  Omega_Ref         = Velocity_Ref/Length_Ref;                                     config->SetOmega_Ref(Omega_Ref);
  Force_Ref         = config->GetDensity_Ref()*Velocity_Ref*Velocity_Ref*Length_Ref*Length_Ref; config->SetForce_Ref(Force_Ref);
  Heat_Flux_Ref     = Density_Ref*Velocity_Ref*Velocity_Ref*Velocity_Ref;           config->SetHeat_Flux_Ref(Heat_Flux_Ref);
  Gas_Constant_Ref  = Velocity_Ref*Velocity_Ref/config->GetTemperature_Ref();      config->SetGas_Constant_Ref(Gas_Constant_Ref);
  Viscosity_Ref     = config->GetDensity_Ref()*Velocity_Ref*Length_Ref;            config->SetViscosity_Ref(Viscosity_Ref);
  Conductivity_Ref  = Viscosity_Ref*Gas_Constant_Ref;                              config->SetConductivity_Ref(Conductivity_Ref);
  Froude            = ModVel_FreeStream/sqrt(STANDARD_GRAVITY*Length_Ref);         config->SetFroude(Froude);

  /*--- Divide by reference values, to compute the non-dimensional free-stream values ---*/

  Pressure_FreeStreamND = Pressure_FreeStream/config->GetPressure_Ref();  config->SetPressure_FreeStreamND(Pressure_FreeStreamND);
  Density_FreeStreamND  = Density_FreeStream/config->GetDensity_Ref();    config->SetDensity_FreeStreamND(Density_FreeStreamND);

  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_FreeStreamND[iDim] = config->GetVelocity_FreeStream()[iDim]/Velocity_Ref; config->SetVelocity_FreeStreamND(Velocity_FreeStreamND[iDim], iDim);
  }

  Temperature_FreeStreamND = Temperature_FreeStream/config->GetTemperature_Ref(); config->SetTemperature_FreeStreamND(Temperature_FreeStreamND);

  Gas_ConstantND = config->GetGas_Constant()/Gas_Constant_Ref;    config->SetGas_ConstantND(Gas_ConstantND);


  ModVel_FreeStreamND = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) ModVel_FreeStreamND += Velocity_FreeStreamND[iDim]*Velocity_FreeStreamND[iDim];
  ModVel_FreeStreamND    = sqrt(ModVel_FreeStreamND); config->SetModVel_FreeStreamND(ModVel_FreeStreamND);

  /*--- Initialize the dimensionless Fluid Model that will be used to solve the dimensionless problem ---*/

  /*--- Auxilary (dimensional) FluidModel no longer needed. ---*/
  delete auxFluidModel;

  /*--- Create one final fluid model object per OpenMP thread to be able to use them in parallel.
   *    GetFluidModel() should be used to automatically access the "right" object of each thread. ---*/

  assert(FluidModel.empty() && "Potential memory leak!");
  FluidModel.resize(omp_get_max_threads());

  SU2_OMP_PARALLEL
  {
    const int thread = omp_get_thread_num();

    switch (config->GetKind_FluidModel()) {

      case STANDARD_AIR:
        FluidModel[thread] = new CIdealGas(1.4, Gas_ConstantND);
        break;

      case IDEAL_GAS:
        FluidModel[thread] = new CIdealGas(Gamma, Gas_ConstantND);
        break;

    }

    GetFluidModel()->SetEnergy_Prho(Pressure_FreeStreamND, Density_FreeStreamND);

  }
  END_SU2_OMP_PARALLEL

  Energy_FreeStreamND = GetFluidModel()->GetStaticEnergy() + 0.5*ModVel_FreeStreamND*ModVel_FreeStreamND;

  config->SetEnergy_FreeStreamND(Energy_FreeStreamND);

  Energy_Ref = Energy_FreeStream/Energy_FreeStreamND; config->SetEnergy_Ref(Energy_Ref);

  Total_UnstTimeND = config->GetTotal_UnstTime() / Time_Ref;    config->SetTotal_UnstTimeND(Total_UnstTimeND);
  Delta_UnstTimeND = config->GetDelta_UnstTime() / Time_Ref;    config->SetDelta_UnstTimeND(Delta_UnstTimeND);

  /*--- Write output to the console if this is the master node and first domain ---*/

  if ((rank == MASTER_NODE) && (MGLevel == MESH_0)) {

    cout.precision(6);

    cout << "Inviscid flow: Computing density based on free-stream" << endl;
    cout << "temperature and pressure using the ideal gas law." << endl;

    cout << "Force coefficients computed using free-stream values." << endl;

    stringstream NonDimTableOut, ModelTableOut;
    stringstream Unit;

    cout << endl;
    PrintingToolbox::CTablePrinter ModelTable(&ModelTableOut);
    ModelTableOut <<"-- Models:"<< endl;

    ModelTable.AddColumn("Viscosity Model", 25);
    ModelTable.AddColumn("Conductivity Model", 26);
    ModelTable.AddColumn("Fluid Model", 25);
    ModelTable.SetAlign(PrintingToolbox::CTablePrinter::RIGHT);
    ModelTable.PrintHeader();

    PrintingToolbox::CTablePrinter NonDimTable(&NonDimTableOut);
    NonDimTable.AddColumn("Name", 22);
    NonDimTable.AddColumn("Dim. value", 14);
    NonDimTable.AddColumn("Ref. value", 14);
    NonDimTable.AddColumn("Unit", 10);
    NonDimTable.AddColumn("Non-dim. value", 14);
    NonDimTable.SetAlign(PrintingToolbox::CTablePrinter::RIGHT);

    NonDimTableOut <<"-- Fluid properties:"<< endl;

    NonDimTable.PrintHeader();

    ModelTable << "-" << "-";

    if      (config->GetSystemMeasurements() == SI) Unit << "N.m/kg.K";
    else if (config->GetSystemMeasurements() == US) Unit << "lbf.ft/slug.R";
    NonDimTable << "Gas Constant" << config->GetGas_Constant() << config->GetGas_Constant_Ref() << Unit.str() << config->GetGas_ConstantND();
    Unit.str("");
    if      (config->GetSystemMeasurements() == SI) Unit << "N.m/kg.K";
    else if (config->GetSystemMeasurements() == US) Unit << "lbf.ft/slug.R";
    NonDimTable << "Spec. Heat Ratio" << "-" << "-" << "-" << Gamma;
    Unit.str("");

    switch(config->GetKind_FluidModel()){
    case STANDARD_AIR:
      ModelTable << "STANDARD_AIR";
      break;
    case IDEAL_GAS:
      ModelTable << "IDEAL_GAS";
      break;
    }

    NonDimTable.PrintFooter();

    NonDimTableOut <<"-- Initial and free-stream conditions:"<< endl;

    NonDimTable.PrintHeader();

    if      (config->GetSystemMeasurements() == SI) Unit << "Pa";
    else if (config->GetSystemMeasurements() == US) Unit << "psf";
    NonDimTable << "Static Pressure" << config->GetPressure_FreeStream() << config->GetPressure_Ref() << Unit.str() << config->GetPressure_FreeStreamND();
    Unit.str("");
    if      (config->GetSystemMeasurements() == SI) Unit << "kg/m^3";
    else if (config->GetSystemMeasurements() == US) Unit << "slug/ft^3";
    NonDimTable << "Density" << config->GetDensity_FreeStream() << config->GetDensity_Ref() << Unit.str() << config->GetDensity_FreeStreamND();
    Unit.str("");
    if      (config->GetSystemMeasurements() == SI) Unit << "K";
    else if (config->GetSystemMeasurements() == US) Unit << "R";
    NonDimTable << "Temperature" << config->GetTemperature_FreeStream() << config->GetTemperature_Ref() << Unit.str() << config->GetTemperature_FreeStreamND();
    Unit.str("");
    if      (config->GetSystemMeasurements() == SI) Unit << "m^2/s^2";
    else if (config->GetSystemMeasurements() == US) Unit << "ft^2/s^2";
    NonDimTable << "Total Energy" << config->GetEnergy_FreeStream() << config->GetEnergy_Ref() << Unit.str() << config->GetEnergy_FreeStreamND();
    Unit.str("");
    if      (config->GetSystemMeasurements() == SI) Unit << "m/s";
    else if (config->GetSystemMeasurements() == US) Unit << "ft/s";
    NonDimTable << "Velocity-X" << config->GetVelocity_FreeStream()[0] << config->GetVelocity_Ref() << Unit.str() << config->GetVelocity_FreeStreamND()[0];
    NonDimTable << "Velocity-Y" << config->GetVelocity_FreeStream()[1] << config->GetVelocity_Ref() << Unit.str() << config->GetVelocity_FreeStreamND()[1];
    if (nDim == 3){
      NonDimTable << "Velocity-Z" << config->GetVelocity_FreeStream()[2] << config->GetVelocity_Ref() << Unit.str() << config->GetVelocity_FreeStreamND()[2];
    }
    NonDimTable << "Velocity Magnitude" << config->GetModVel_FreeStream() << config->GetVelocity_Ref() << Unit.str() << config->GetModVel_FreeStreamND();
    Unit.str("");

    NonDimTable.PrintFooter();
    NonDimTable << "Mach Number" << "-" << "-" << "-" << config->GetMach();
    if (gravity) {
      NonDimTable << "Froude Number" << "-" << "-" << "-" << Froude;
      NonDimTable << "Wave Length"   << "-" << "-" << "-" << 2.0*PI_NUMBER*Froude*Froude;
    }
    NonDimTable.PrintFooter();
    ModelTable.PrintFooter();

    if (unsteady){
      NonDimTableOut << "-- Unsteady conditions" << endl;
      NonDimTable.PrintHeader();
      NonDimTable << "Total Time" << config->GetMax_Time() << config->GetTime_Ref() << "s" << config->GetMax_Time()/config->GetTime_Ref();
      Unit.str("");
      NonDimTable << "Time Step" << config->GetTime_Step() << config->GetTime_Ref() << "s" << config->GetDelta_UnstTimeND();
      Unit.str("");
      NonDimTable.PrintFooter();
    }

    cout << ModelTableOut.str();
    cout << NonDimTableOut.str();

  }

}

void CEulerSolver::SetReferenceValues(const CConfig& config) {

  /*--- Evaluate reference values for non-dimensionalization. For dynamic meshes,
   use the motion Mach number as a reference value for computing the force coefficients.
   Otherwise, use the freestream values, which is the standard convention. ---*/

  su2double RefVel2;

  RefVel2 = GeometryToolbox::SquaredNorm(nDim, Velocity_Inf);

  DynamicPressureRef = 0.5 * Density_Inf * RefVel2;
  AeroCoeffForceRef =  DynamicPressureRef * config.GetRefArea();

}

void CEulerSolver::SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long TimeIter) {

  const bool restart = (config->GetRestart() || config->GetRestart_Flow());
  const bool SubsonicEngine = config->GetSubsonicEngine();

  /*--- Use default implementation, then add solver-specifics. ---*/

  BaseClass::SetInitialCondition(geometry, solver_container, config, TimeIter);

  /*--- Set subsonic initial condition for engine intakes at iteration 0 ---*/

  if (!SubsonicEngine || (TimeIter != 0) || restart) return;

  /*--- Start OpenMP parallel region. ---*/

  SU2_OMP_PARALLEL {

    unsigned long iPoint;
    unsigned short iMesh, iDim;
    su2double X0[MAXNDIM] = {0.0}, X1[MAXNDIM] = {0.0}, X2[MAXNDIM] = {0.0},
    X1_X0[MAXNDIM] = {0.0}, X2_X0[MAXNDIM] = {0.0}, X2_X1[MAXNDIM] = {0.0},
    CP[MAXNDIM] = {0.0}, Distance, DotCheck, Radius;

    su2double Velocity_Cyl[MAXNDIM] = {0.0}, Velocity_CylND[MAXNDIM] = {0.0},
    Density_Cyl, Density_CylND, Pressure_CylND, ModVel_CylND, Energy_CylND;
    const su2double *Coord, *SubsonicEngine_Cyl, *SubsonicEngine_Values;

    SubsonicEngine_Values = config->GetSubsonicEngine_Values();
    su2double Mach_Cyl        = SubsonicEngine_Values[0];
    su2double Alpha_Cyl       = SubsonicEngine_Values[1];
    su2double Beta_Cyl        = SubsonicEngine_Values[2];
    su2double Pressure_Cyl    = SubsonicEngine_Values[3];
    su2double Temperature_Cyl = SubsonicEngine_Values[4];

    su2double Alpha = Alpha_Cyl*PI_NUMBER/180.0;
    su2double Beta  = Beta_Cyl*PI_NUMBER/180.0;

    su2double Gamma_Minus_One = Gamma - 1.0;
    su2double Gas_Constant = config->GetGas_Constant();

    su2double Mach2Vel_Cyl = sqrt(Gamma*Gas_Constant*Temperature_Cyl);

    for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {

      auto FlowNodes = solver_container[iMesh][FLOW_SOL]->GetNodes();

      SU2_OMP_FOR_STAT(omp_chunk_size)
      for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {

        Velocity_Cyl[0] = cos(Alpha)*cos(Beta)*Mach_Cyl*Mach2Vel_Cyl;
        Velocity_Cyl[1] = sin(Beta)*Mach_Cyl*Mach2Vel_Cyl;
        Velocity_Cyl[2] = sin(Alpha)*cos(Beta)*Mach_Cyl*Mach2Vel_Cyl;

        Density_Cyl = Pressure_Cyl/(Gas_Constant*Temperature_Cyl);

        Density_CylND  = Density_Cyl/config->GetDensity_Ref();
        Pressure_CylND = Pressure_Cyl/config->GetPressure_Ref();

        for (iDim = 0; iDim < nDim; iDim++) {
          Velocity_CylND[iDim] = Velocity_Cyl[iDim]/config->GetVelocity_Ref();
        }

        ModVel_CylND = GeometryToolbox::Norm(nDim, Velocity_CylND);

        Energy_CylND = Pressure_CylND/(Density_CylND*Gamma_Minus_One)+0.5*ModVel_CylND*ModVel_CylND;

        Coord = geometry[iMesh]->nodes->GetCoord(iPoint);

        SubsonicEngine_Cyl = config->GetSubsonicEngine_Cyl();

        X0[0] = Coord[0];               X0[1] = Coord[1];               if (nDim==3) X0[2] = Coord[2];
        X1[0] = SubsonicEngine_Cyl[0];  X1[1] = SubsonicEngine_Cyl[1];  X1[2] = SubsonicEngine_Cyl[2];
        X2[0] = SubsonicEngine_Cyl[3];  X2[1] = SubsonicEngine_Cyl[4];  X2[2] = SubsonicEngine_Cyl[5];
        Radius = SubsonicEngine_Cyl[6];

        GeometryToolbox::Distance(3, X1, X2, X2_X1);
        GeometryToolbox::Distance(3, X0, X1, X1_X0);
        GeometryToolbox::Distance(3, X0, X2, X2_X0);

        GeometryToolbox::CrossProduct(X2_X1, X1_X0, CP);

        Distance = sqrt(GeometryToolbox::SquaredNorm(3,CP) / GeometryToolbox::SquaredNorm(3,X2_X1));

        DotCheck = -GeometryToolbox::DotProduct(3, X1_X0, X2_X1);
        if (DotCheck < 0.0) Distance = GeometryToolbox::Norm(3, X1_X0);

        DotCheck = GeometryToolbox::DotProduct(3, X2_X0, X2_X1);
        if (DotCheck < 0.0) Distance = GeometryToolbox::Norm(3, X2_X0);

        if (Distance < Radius) {
          FlowNodes->SetSolution(iPoint, 0, Density_CylND);
          for (iDim = 0; iDim < nDim; iDim++)
            FlowNodes->SetSolution(iPoint, iDim+1, Density_CylND*Velocity_CylND[iDim]);
          FlowNodes->SetSolution(iPoint, nVar-1, Density_CylND*Energy_CylND);
        }

      }
      END_SU2_OMP_FOR

      FlowNodes->Set_OldSolution();

    }

  }
  END_SU2_OMP_PARALLEL

}

void CEulerSolver::CommonPreprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh,
                                       unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {

  bool cont_adjoint     = config->GetContinuous_Adjoint();
  bool disc_adjoint     = config->GetDiscrete_Adjoint();
  bool implicit         = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  bool center           = (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED);
  bool center_jst       = (config->GetKind_Centered_Flow() == JST) && (iMesh == MESH_0);
  bool center_jst_ke    = (config->GetKind_Centered_Flow() == JST_KE) && (iMesh == MESH_0);
  bool center_jst_mat   = (config->GetKind_Centered_Flow() == JST_MAT) && (iMesh == MESH_0);
  bool fixed_cl         = config->GetFixed_CL_Mode();
  unsigned short kind_row_dissipation = config->GetKind_RoeLowDiss();
  bool roe_low_dissipation  = (kind_row_dissipation != NO_ROELOWDISS) &&
                              (config->GetKind_Upwind_Flow() == ROE ||
                               config->GetKind_Upwind_Flow() == SLAU ||
                               config->GetKind_Upwind_Flow() == SLAU2);

  /*--- Set the primitive variables ---*/

  ompMasterAssignBarrier(ErrorCounter, 0);

  SU2_OMP_ATOMIC
  ErrorCounter += SetPrimitive_Variables(solver_container, config);
  SU2_OMP_BARRIER

  SU2_OMP_MASTER { /*--- Ops that are not OpenMP parallel go in this block. ---*/

    if ((iMesh == MESH_0) && (config->GetComm_Level() == COMM_FULL)) {
      unsigned long tmp = ErrorCounter;
      SU2_MPI::Allreduce(&tmp, &ErrorCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
      config->SetNonphysical_Points(ErrorCounter);
    }

    /*--- Update the angle of attack at the far-field for fixed CL calculations (only direct problem). ---*/

    if (fixed_cl && !disc_adjoint && !cont_adjoint) {
      SetFarfield_AoA(geometry, solver_container, config, iMesh, Output);
    }

  }
  END_SU2_OMP_MASTER
  SU2_OMP_BARRIER

  /*--- Artificial dissipation ---*/

  if (center && !Output) {
    if (!center_jst_mat) SetMax_Eigenvalue(geometry, config);
    if (center_jst || center_jst_ke || center_jst_mat) {
      SetCentered_Dissipation_Sensor(geometry, config);
      if (!center_jst_ke) SetUndivided_Laplacian(geometry, config);
    }
  }

  /*--- Roe Low Dissipation Sensor ---*/

  if (roe_low_dissipation) {
    SetRoe_Dissipation(geometry, config);
    if (kind_row_dissipation == FD_DUCROS || kind_row_dissipation == NTS_DUCROS){
      SetUpwind_Ducros_Sensor(geometry, config);
    }
  }

  /*--- Initialize the Jacobian matrix and residual, not needed for the reducer strategy
   *    as we set blocks (including diagonal ones) and completely overwrite. ---*/

  if(!ReducerStrategy && !Output) {
    LinSysRes.SetValZero();
    if (implicit) Jacobian.SetValZero();
    else {SU2_OMP_BARRIER} // because of "nowait" in LinSysRes
  }

}

void CEulerSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh,
                                 unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
  const auto InnerIter = config->GetInnerIter();
  const bool muscl = config->GetMUSCL_Flow() && (iMesh == MESH_0);
  const bool center = (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED);
  const bool limiter = (config->GetKind_SlopeLimit_Flow() != LIMITER::NONE) && (InnerIter <= config->GetLimiterIter());
  const bool van_albada = (config->GetKind_SlopeLimit_Flow() == LIMITER::VAN_ALBADA_EDGE);

  /*--- Common preprocessing steps. ---*/

  CommonPreprocessing(geometry, solver_container, config, iMesh, iRKStep, RunTime_EqSystem, Output);

  /*--- Upwind second order reconstruction ---*/

  if (!Output && muscl && !center) {

    /*--- Gradient computation for MUSCL reconstruction. ---*/

    switch (config->GetKind_Gradient_Method_Recon()) {
      case GREEN_GAUSS:
        SetPrimitive_Gradient_GG(geometry, config, true); break;
      case LEAST_SQUARES:
      case WEIGHTED_LEAST_SQUARES:
        SetPrimitive_Gradient_LS(geometry, config, true); break;
      default: break;
    }

    /*--- Limiter computation ---*/

    if (limiter && !van_albada) SetPrimitive_Limiter(geometry, config);
  }
}

unsigned long CEulerSolver::SetPrimitive_Variables(CSolver **solver_container, const CConfig *config) {

  /*--- Number of non-physical points, local to the thread, needs
   *    further reduction if function is called in parallel ---*/
  unsigned long nonPhysicalPoints = 0;


  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint ++) {

    /*--- Compressible flow, primitive variables nDim+9, (T, vx, vy, vz, P, rho, h, c, lamMu, eddyMu, ThCond, Cp) ---*/

    bool physical = nodes->SetPrimVar(iPoint, GetFluidModel());
    nodes->SetSecondaryVar(iPoint, GetFluidModel());

    /* Check for non-realizable states for reporting. */

    if (!physical) nonPhysicalPoints++;
  }
  END_SU2_OMP_FOR


  return nonPhysicalPoints;
}

void CEulerSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                unsigned short iMesh, unsigned long Iteration) {

  /*--- Define an object to compute the speed of sound. ---*/
  struct SoundSpeed {
    FORCEINLINE su2double operator() (const CEulerVariable& nodes, unsigned long iPoint, unsigned long jPoint) const {
      return 0.5 * (nodes.GetSoundSpeed(iPoint) + nodes.GetSoundSpeed(jPoint));
    }

    FORCEINLINE su2double operator() (const CEulerVariable& nodes, unsigned long iPoint) const {
      return nodes.GetSoundSpeed(iPoint);
    }

  } soundSpeed;

  /*--- Now instantiate the generic implementation with the two functors above. ---*/

  SetTime_Step_impl(soundSpeed, geometry, solver_container, config, iMesh, Iteration);

}

void CEulerSolver::Centered_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics_container,
                                     CConfig *config, unsigned short iMesh, unsigned short iRKStep) {

  EdgeFluxResidual(geometry, solver_container, config);
}

void CEulerSolver::Upwind_Residual(CGeometry *geometry, CSolver **solver_container,
                                   CNumerics **numerics_container, CConfig *config, unsigned short iMesh) {

  if (config->GetUseVectorization()) {
    EdgeFluxResidual(geometry, solver_container, config);
    return;
  }

  const bool implicit         = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  const bool ideal_gas        = (config->GetKind_FluidModel() == STANDARD_AIR) ||
                                (config->GetKind_FluidModel() == IDEAL_GAS);

  const bool roe_turkel       = (config->GetKind_Upwind_Flow() == TURKEL);
  const bool low_mach_corr    = config->Low_Mach_Correction();
  const auto kind_dissipation = config->GetKind_RoeLowDiss();

  const bool muscl            = (config->GetMUSCL_Flow() && (iMesh == MESH_0));
  const bool limiter          = (config->GetKind_SlopeLimit_Flow() != LIMITER::NONE);
  const bool van_albada       = (config->GetKind_SlopeLimit_Flow() == LIMITER::VAN_ALBADA_EDGE);

  /*--- Non-physical counter. ---*/
  unsigned long counter_local = 0;
  SU2_OMP_MASTER
  ErrorCounter = 0;
  END_SU2_OMP_MASTER

  /*--- Pick one numerics object per thread. ---*/
  CNumerics* numerics = numerics_container[CONV_TERM + omp_get_thread_num()*MAX_TERMS];

  /*--- Static arrays of MUSCL-reconstructed primitives and secondaries (thread safety). ---*/
  su2double Primitive_i[MAXNVAR] = {0.0}, Primitive_j[MAXNVAR] = {0.0};
  su2double Secondary_i[MAXNVAR] = {0.0}, Secondary_j[MAXNVAR] = {0.0};

  /*--- For hybrid parallel AD, pause preaccumulation if there is shared reading of
  * variables, otherwise switch to the faster adjoint evaluation mode. ---*/

  /*--- Loop over edge colors. ---*/
  for (auto color : EdgeColoring)
  {
  /*--- Chunk size is at least OMP_MIN_SIZE and a multiple of the color group size. ---*/
  SU2_OMP_FOR_DYN(nextMultiple(OMP_MIN_SIZE, color.groupSize))
  for(auto k = 0ul; k < color.size; ++k) {

    auto iEdge = color.indices[k];

    unsigned short iDim, iVar;

    /*--- Points in edge and normal vectors ---*/

    auto iPoint = geometry->edges->GetNode(iEdge,0);
    auto jPoint = geometry->edges->GetNode(iEdge,1);

    numerics->SetNormal(geometry->edges->GetNormal(iEdge));

    auto Coord_i = geometry->nodes->GetCoord(iPoint);
    auto Coord_j = geometry->nodes->GetCoord(jPoint);

    /*--- Roe Turkel preconditioning ---*/

    if (roe_turkel) {
      numerics->SetVelocity2_Inf(GeometryToolbox::SquaredNorm(nDim, config->GetVelocity_FreeStream()));
    }

    /*--- Get primitive and secondary variables ---*/

    auto V_i = nodes->GetPrimitive(iPoint); auto V_j = nodes->GetPrimitive(jPoint);
    auto S_i = nodes->GetSecondary(iPoint); auto S_j = nodes->GetSecondary(jPoint);

    /*--- Set them with or without high order reconstruction using MUSCL strategy. ---*/

    if (!muscl) {

      numerics->SetPrimitive(V_i, V_j);
      numerics->SetSecondary(S_i, S_j);

    }
    else {
      /*--- Reconstruction ---*/

      su2double Vector_ij[MAXNDIM] = {0.0};
      for (iDim = 0; iDim < nDim; iDim++) {
        Vector_ij[iDim] = 0.5*(Coord_j[iDim] - Coord_i[iDim]);
      }

      auto Gradient_i = nodes->GetGradient_Reconstruction(iPoint);
      auto Gradient_j = nodes->GetGradient_Reconstruction(jPoint);

      for (iVar = 0; iVar < nPrimVarGrad; iVar++) {

        su2double Project_Grad_i = 0.0;
        su2double Project_Grad_j = 0.0;

        for (iDim = 0; iDim < nDim; iDim++) {
          Project_Grad_i += Vector_ij[iDim]*Gradient_i[iVar][iDim];
          Project_Grad_j -= Vector_ij[iDim]*Gradient_j[iVar][iDim];
        }

        su2double lim_i = 1.0;
        su2double lim_j = 1.0;

        if (van_albada) {
          su2double V_ij = V_j[iVar] - V_i[iVar];
          lim_i = LimiterHelpers<>::vanAlbadaFunction(Project_Grad_i, V_ij, EPS);
          lim_j = LimiterHelpers<>::vanAlbadaFunction(-Project_Grad_j, V_ij, EPS);
        }
        else if (limiter) {
          lim_i = nodes->GetLimiter_Primitive(iPoint, iVar);
          lim_j = nodes->GetLimiter_Primitive(jPoint, iVar);
        }

        Primitive_i[iVar] = V_i[iVar] + lim_i * Project_Grad_i;
        Primitive_j[iVar] = V_j[iVar] + lim_j * Project_Grad_j;

      }

      /*--- Recompute the reconstructed quantities in a thermodynamically consistent way. ---*/

      if (!ideal_gas || low_mach_corr) {
        ComputeConsistentExtrapolation(GetFluidModel(), nDim, Primitive_i, Secondary_i);
        ComputeConsistentExtrapolation(GetFluidModel(), nDim, Primitive_j, Secondary_j);
      }

      /*--- Low-Mach number correction. ---*/

      if (low_mach_corr) {
        LowMachPrimitiveCorrection(GetFluidModel(), nDim, Primitive_i, Primitive_j);
      }

      /*--- Check for non-physical solutions after reconstruction. If found, use the
       cell-average value of the solution. This is a locally 1st order approximation,
       which is typically only active during the start-up of a calculation. ---*/

      bool neg_pres_or_rho_i = (Primitive_i[prim_idx.Pressure()] < 0.0) || (Primitive_i[prim_idx.Density()] < 0.0);
      bool neg_pres_or_rho_j = (Primitive_j[prim_idx.Pressure()] < 0.0) || (Primitive_j[prim_idx.Density()] < 0.0);

      su2double R = sqrt(fabs(Primitive_j[prim_idx.Density()]/Primitive_i[prim_idx.Density()]));
      su2double sq_vel = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        su2double RoeVelocity = (R * Primitive_j[iDim + prim_idx.Velocity()] +
                                 Primitive_i[iDim + prim_idx.Velocity()]) / (R+1);
        sq_vel += pow(RoeVelocity, 2);
      }
      su2double RoeEnthalpy = (R * Primitive_j[prim_idx.Enthalpy()] + Primitive_i[prim_idx.Enthalpy()]) / (R+1);

      bool neg_sound_speed = ((Gamma-1)*(RoeEnthalpy-0.5*sq_vel) < 0.0);

      bool bad_i = neg_sound_speed || neg_pres_or_rho_i;
      bool bad_j = neg_sound_speed || neg_pres_or_rho_j;

      nodes->SetNon_Physical(iPoint, bad_i);
      nodes->SetNon_Physical(jPoint, bad_j);

      /*--- Get updated state, in case the point recovered after the set. ---*/
      bad_i = nodes->GetNon_Physical(iPoint);
      bad_j = nodes->GetNon_Physical(jPoint);

      counter_local += bad_i+bad_j;

      numerics->SetPrimitive(bad_i? V_i : Primitive_i,  bad_j? V_j : Primitive_j);
      numerics->SetSecondary(bad_i? S_i : Secondary_i,  bad_j? S_j : Secondary_j);

    }

    /*--- Roe Low Dissipation Scheme ---*/

    if (kind_dissipation != NO_ROELOWDISS) {

      numerics->SetDissipation(nodes->GetRoe_Dissipation(iPoint),
                               nodes->GetRoe_Dissipation(jPoint));

      if (kind_dissipation == FD_DUCROS || kind_dissipation == NTS_DUCROS){
        numerics->SetSensor(nodes->GetSensor(iPoint),
                            nodes->GetSensor(jPoint));
      }
      if (kind_dissipation == NTS || kind_dissipation == NTS_DUCROS){
        numerics->SetCoord(Coord_i, Coord_j);
      }
    }

    /*--- Compute the residual ---*/

    auto residual = numerics->ComputeResidual(config);

    /*--- Set the final value of the Roe dissipation coefficient ---*/

    if ((kind_dissipation != NO_ROELOWDISS) && (MGLevel != MESH_0)) {
      nodes->SetRoe_Dissipation(iPoint,numerics->GetDissipation());
      nodes->SetRoe_Dissipation(jPoint,numerics->GetDissipation());
    }

    /*--- Update residual value ---*/

    if (ReducerStrategy) {
      EdgeFluxes.SetBlock(iEdge, residual);
      if (implicit)
        Jacobian.SetBlocks(iEdge, residual.jacobian_i, residual.jacobian_j);
    }
    else {
      LinSysRes.AddBlock(iPoint, residual);
      LinSysRes.SubtractBlock(jPoint, residual);

      /*--- Set implicit computation ---*/
      if (implicit)
        Jacobian.UpdateBlocks(iEdge, iPoint, jPoint, residual.jacobian_i, residual.jacobian_j);
    }

  }
  END_SU2_OMP_FOR
  } // end color loop

  if (ReducerStrategy) {
    SumEdgeFluxes(geometry);
    if (implicit)
      Jacobian.SetDiagonalAsColumnSum();
  }

  /*--- Warning message about non-physical reconstructions. ---*/

  if ((iMesh == MESH_0) && (config->GetComm_Level() == COMM_FULL)) {
    /*--- Add counter results for all threads. ---*/
    SU2_OMP_ATOMIC
    ErrorCounter += counter_local;
    SU2_OMP_BARRIER

    /*--- Add counter results for all ranks. ---*/
    SU2_OMP_MASTER {
      counter_local = ErrorCounter;
      SU2_MPI::Reduce(&counter_local, &ErrorCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, SU2_MPI::GetComm());
      config->SetNonphysical_Reconstr(ErrorCounter);
    }
    END_SU2_OMP_MASTER
    SU2_OMP_BARRIER
  }

}

void CEulerSolver::ComputeConsistentExtrapolation(CFluidModel *fluidModel, unsigned short nDim,
                                                  su2double *primitive, su2double *secondary) {
  const CEulerVariable::CIndices<unsigned short> prim_idx(nDim, 0);
  const su2double density = primitive[prim_idx.Density()];
  const su2double pressure = primitive[prim_idx.Pressure()];
  const su2double velocity2 = GeometryToolbox::SquaredNorm(nDim, &primitive[prim_idx.Velocity()]);

  fluidModel->SetTDState_Prho(pressure, density);

  primitive[prim_idx.Temperature()] = fluidModel->GetTemperature();
  primitive[prim_idx.Enthalpy()] = fluidModel->GetStaticEnergy() + pressure / density + 0.5*velocity2;
  primitive[prim_idx.SoundSpeed()] = fluidModel->GetSoundSpeed();
  secondary[0] = fluidModel->GetdPdrho_e();
  secondary[1] = fluidModel->GetdPde_rho();

}

void CEulerSolver::LowMachPrimitiveCorrection(CFluidModel *fluidModel, unsigned short nDim,
                                              su2double *primitive_i, su2double *primitive_j) {
  unsigned short iDim;

  su2double velocity2_i = 0.0;
  su2double velocity2_j = 0.0;

  for (iDim = 0; iDim < nDim; iDim++) {
    velocity2_i += pow(primitive_i[iDim+1], 2);
    velocity2_j += pow(primitive_j[iDim+1], 2);
  }
  su2double mach_i = sqrt(velocity2_i)/primitive_i[nDim+4];
  su2double mach_j = sqrt(velocity2_j)/primitive_j[nDim+4];

  su2double z = min(max(mach_i,mach_j),1.0);
  velocity2_i = 0.0;
  velocity2_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    su2double vel_i_corr = ( primitive_i[iDim+1] + primitive_j[iDim+1] )/2.0
                     + z * ( primitive_i[iDim+1] - primitive_j[iDim+1] )/2.0;
    su2double vel_j_corr = ( primitive_i[iDim+1] + primitive_j[iDim+1] )/2.0
                     + z * ( primitive_j[iDim+1] - primitive_i[iDim+1] )/2.0;

    velocity2_i += pow(vel_i_corr, 2);
    velocity2_j += pow(vel_j_corr, 2);

    primitive_i[iDim+1] = vel_i_corr;
    primitive_j[iDim+1] = vel_j_corr;
  }

  fluidModel->SetEnergy_Prho(primitive_i[nDim+1], primitive_i[nDim+2]);
  primitive_i[nDim+3]= fluidModel->GetStaticEnergy() + primitive_i[nDim+1]/primitive_i[nDim+2] + 0.5*velocity2_i;

  fluidModel->SetEnergy_Prho(primitive_j[nDim+1], primitive_j[nDim+2]);
  primitive_j[nDim+3]= fluidModel->GetStaticEnergy() + primitive_j[nDim+1]/primitive_j[nDim+2] + 0.5*velocity2_j;

}

void CEulerSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container,
                                   CNumerics **numerics_container, CConfig *config, unsigned short iMesh) {

  const bool implicit         = config->GetKind_TimeIntScheme() == EULER_IMPLICIT;
  const bool rotating_frame   = config->GetRotating_Frame();
  const bool gravity          = (config->GetGravityForce() == YES);
  const bool body_force       = config->GetBody_Force();

  /*--- Pick one numerics object per thread. ---*/
  CNumerics* numerics = numerics_container[SOURCE_FIRST_TERM + omp_get_thread_num()*MAX_TERMS];

  unsigned long iPoint;

  if (body_force) {

    /*--- Loop over all points ---*/
    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

      /*--- Load the conservative variables ---*/
      numerics->SetConservative(nodes->GetSolution(iPoint),
                                nodes->GetSolution(iPoint));

      /*--- Load the volume of the dual mesh cell ---*/
      numerics->SetVolume(geometry->nodes->GetVolume(iPoint));

      /*--- Compute the rotating frame source residual ---*/
      auto residual = numerics->ComputeResidual(config);

      /*--- Add the source residual to the total ---*/
      LinSysRes.AddBlock(iPoint, residual);

    }
    END_SU2_OMP_FOR
  }

  if (rotating_frame) {

    /*--- Include the residual contribution from GCL due to the static
     mesh movement that is set for rotating frame. ---*/

    SetRotatingFrame_GCL(geometry, config);

    /*--- Loop over all points ---*/
    SU2_OMP_FOR_DYN(omp_chunk_size)
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

      /*--- Load the conservative variables ---*/
      numerics->SetConservative(nodes->GetSolution(iPoint),
                                nodes->GetSolution(iPoint));

      /*--- Load the volume of the dual mesh cell ---*/
      numerics->SetVolume(geometry->nodes->GetVolume(iPoint));

      /*--- Compute the rotating frame source residual ---*/
      auto residual = numerics->ComputeResidual(config);

      /*--- Add the source residual to the total ---*/
      LinSysRes.AddBlock(iPoint, residual);

      /*--- Add the implicit Jacobian contribution ---*/
      if (implicit) Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);

    }
    END_SU2_OMP_FOR
  }

  if (gravity) {

    /*--- loop over points ---*/
    SU2_OMP_FOR_DYN(omp_chunk_size)
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

      /*--- Set solution  ---*/
      numerics->SetConservative(nodes->GetSolution(iPoint), nodes->GetSolution(iPoint));

      /*--- Set control volume ---*/
      numerics->SetVolume(geometry->nodes->GetVolume(iPoint));

      /*--- Compute Source term Residual ---*/
      auto residual = numerics->ComputeResidual(config);

      /*--- Add Residual ---*/
      LinSysRes.AddBlock(iPoint, residual);

    }
    END_SU2_OMP_FOR

  }

}

void CEulerSolver::Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                   CConfig *config, unsigned short iMesh) {

  /* This method should be used to call any new source terms for a particular problem*/
  /* This method calls the new child class in CNumerics, where the new source term should be implemented.  */

  /* Next we describe how to get access to some important quanties for this method */
  /* Access to all points in the current geometric mesh by saying: nPointDomain */
  /* Get the vector of conservative variables at some point iPoint = nodes->GetSolution(iPoint) */
  /* Get the volume (or area in 2D) associated with iPoint = nodes->GetVolume(iPoint) */
  /* Get the vector of geometric coordinates of point iPoint = nodes->GetCoord(iPoint) */

}

void CEulerSolver::SetMax_Eigenvalue(CGeometry *geometry, const CConfig *config) {

  /*--- Define an object to compute the speed of sound. ---*/
  struct SoundSpeed {
    FORCEINLINE su2double operator() (const CEulerVariable& nodes, unsigned long iPoint, unsigned long jPoint) const {
      return 0.5 * (nodes.GetSoundSpeed(iPoint) + nodes.GetSoundSpeed(jPoint));
    }

    FORCEINLINE su2double operator() (const CEulerVariable& nodes, unsigned long iPoint) const {
      return nodes.GetSoundSpeed(iPoint);
    }

  } soundSpeed;

  /*--- Instantiate generic implementation. ---*/

  SetMax_Eigenvalue_impl(soundSpeed, geometry, config);

}

void CEulerSolver::SetUndivided_Laplacian(CGeometry *geometry, const CConfig *config) {

  /*--- Loop domain points. ---*/

  SU2_OMP_FOR_DYN(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPointDomain; ++iPoint) {

    const bool boundary_i = geometry->nodes->GetPhysicalBoundary(iPoint);
    const su2double Pressure_i = nodes->GetPressure(iPoint);

    /*--- Initialize. ---*/
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
      nodes->SetUnd_Lapl(iPoint, iVar, 0.0);

    /*--- Loop over the neighbors of point i. ---*/
    for (auto jPoint : geometry->nodes->GetPoints(iPoint)) {

      bool boundary_j = geometry->nodes->GetPhysicalBoundary(jPoint);

      /*--- If iPoint is boundary it only takes contributions from other boundary points. ---*/
      if (boundary_i && !boundary_j) continue;

      /*--- Add solution differences, with correction for compressible flows which use the enthalpy. ---*/

      for (unsigned short iVar = 0; iVar < nVar; iVar++)
        nodes->AddUnd_Lapl(iPoint, iVar, nodes->GetSolution(jPoint,iVar)-nodes->GetSolution(iPoint,iVar));

      su2double Pressure_j = nodes->GetPressure(jPoint);
      nodes->AddUnd_Lapl(iPoint, nVar-1, Pressure_j-Pressure_i);
    }
  }
  END_SU2_OMP_FOR

  /*--- Correct the Laplacian across any periodic boundaries. ---*/

  for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
    InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_LAPLACIAN);
    CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_LAPLACIAN);
  }

  /*--- MPI parallelization ---*/

  InitiateComms(geometry, config, UNDIVIDED_LAPLACIAN);
  CompleteComms(geometry, config, UNDIVIDED_LAPLACIAN);

}

void CEulerSolver::SetCentered_Dissipation_Sensor(CGeometry *geometry, const CConfig *config) {

  /*--- Define an object for the sensor variable, pressure. ---*/
  struct SensVar {
    FORCEINLINE su2double operator() (const CEulerVariable& nodes, unsigned long iPoint) const {
      return nodes.GetPressure(iPoint);
    }
  } sensVar;

  /*--- Instantiate generic implementation. ---*/
  SetCentered_Dissipation_Sensor_impl(sensVar, geometry, config);
}

void CEulerSolver::SetUpwind_Ducros_Sensor(CGeometry *geometry, CConfig *config){

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {

    /*---- Ducros sensor for iPoint and its neighbor points to avoid lower dissipation near shocks. ---*/

    su2double Ducros_i = 0.0;
    const auto nNeigh = geometry->nodes->GetnPoint(iPoint);

    for (unsigned short iNeigh = 0; iNeigh <= nNeigh; iNeigh++) {

      auto jPoint = iPoint; // when iNeigh == nNeigh
      if (iNeigh < nNeigh) jPoint = geometry->nodes->GetPoint(iPoint, iNeigh);

      /*---- Dilatation for jPoint ---*/

      su2double uixi=0.0;
      for(unsigned short iDim = 0; iDim < nDim; iDim++){
        uixi += nodes->GetGradient_Primitive(jPoint,iDim+1, iDim);
      }

      /*--- Compute norm of vorticity ---*/

      const su2double* Vorticity = nodes->GetVorticity(jPoint);
      su2double Omega = 0.0;
      for (unsigned short iDim = 0; iDim < nDim; iDim++) {
        Omega += pow(Vorticity[iDim], 2);
      }
      Omega = sqrt(Omega);

      su2double Ducros_j = 0.0;

      if (config->GetKind_RoeLowDiss() == FD_DUCROS) {
        Ducros_j = -uixi / (fabs(uixi) + Omega + 1e-20);
      }
      else if (config->GetKind_RoeLowDiss() == NTS_DUCROS) {
        Ducros_j = pow(uixi,2.0) /(pow(uixi,2.0)+ pow(Omega,2.0) + 1e-20);
      }
      Ducros_i = max(Ducros_i, Ducros_j);
    }

    nodes->SetSensor(iPoint, Ducros_i);
  }
  END_SU2_OMP_FOR

  InitiateComms(geometry, config, SENSOR);
  CompleteComms(geometry, config, SENSOR);

}

void CEulerSolver::ExplicitRK_Iteration(CGeometry *geometry, CSolver **solver_container,
                                        CConfig *config, unsigned short iRKStep) {

  Explicit_Iteration<RUNGE_KUTTA_EXPLICIT>(geometry, solver_container, config, iRKStep);
}

void CEulerSolver::ClassicalRK4_Iteration(CGeometry *geometry, CSolver **solver_container,
                                        CConfig *config, unsigned short iRKStep) {

  Explicit_Iteration<CLASSICAL_RK4_EXPLICIT>(geometry, solver_container, config, iRKStep);
}

void CEulerSolver::ExplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {

  Explicit_Iteration<EULER_EXPLICIT>(geometry, solver_container, config, 0);
}

void CEulerSolver::PrepareImplicitIteration(CGeometry *geometry, CSolver**, CConfig *config) {

  struct LowMachPrec {
    const CEulerSolver* solver;
    const bool active;
    su2activematrix matrix;

    LowMachPrec(const CEulerSolver* s, bool a, unsigned short nVar) : solver(s), active(a) {
      if (active) matrix.resize(nVar,nVar);
    }

    FORCEINLINE const su2activematrix& operator() (const CConfig* config, unsigned long iPoint, su2double delta) {
      solver->SetPreconditioner(config, iPoint, delta, matrix);
      return matrix;
    }

  } precond(this, config->Low_Mach_Preconditioning() || (config->GetKind_Upwind_Flow() == TURKEL), nVar);

  PrepareImplicitIteration_impl(precond, geometry, config);
}

void CEulerSolver::CompleteImplicitIteration(CGeometry *geometry, CSolver**, CConfig *config) {

  CompleteImplicitIteration_impl<true>(geometry, config);
}

void CEulerSolver::SetPreconditioner(const CConfig *config, unsigned long iPoint,
                                     su2double delta, su2activematrix& preconditioner) const {

  unsigned short iDim, jDim, iVar, jVar;
  su2double local_Mach, rho, enthalpy, soundspeed, sq_vel;
  su2double *U_i = nullptr;
  su2double Beta_max = config->GetmaxTurkelBeta();
  su2double Mach_infty2, Mach_lim2, aux, parameter;

  /*--- Variables to calculate the preconditioner parameter Beta ---*/
  local_Mach = sqrt(nodes->GetVelocity2(iPoint))/nodes->GetSoundSpeed(iPoint);

  /*--- Weiss and Smith Preconditioning---*/
  Mach_infty2 = pow(config->GetMach(),2.0);
  Mach_lim2 = pow(0.00001,2.0);
  aux = max(pow(local_Mach,2.0),Mach_lim2);
  parameter = min(1.0, max(aux,Beta_max*Mach_infty2));

  U_i = nodes->GetSolution(iPoint);

  rho = U_i[0];
  enthalpy = nodes->GetEnthalpy(iPoint);
  soundspeed = nodes->GetSoundSpeed(iPoint);
  sq_vel = nodes->GetVelocity2(iPoint);

  /*---Calculating the inverse of the preconditioning matrix that multiplies the time derivative  */
  preconditioner[0][0] = 0.5*sq_vel;
  preconditioner[0][nVar-1] = 1.0;
  for (iDim = 0; iDim < nDim; iDim ++)
    preconditioner[0][1+iDim] = -1.0*U_i[iDim+1]/rho;

  for (iDim = 0; iDim < nDim; iDim ++) {
    preconditioner[iDim+1][0] = 0.5*sq_vel*U_i[iDim+1]/rho;
    preconditioner[iDim+1][nVar-1] = U_i[iDim+1]/rho;
    for (jDim = 0; jDim < nDim; jDim ++) {
      preconditioner[iDim+1][1+jDim] = -1.0*U_i[jDim+1]/rho*U_i[iDim+1]/rho;
    }
  }

  preconditioner[nVar-1][0] = 0.5*sq_vel*enthalpy;
  preconditioner[nVar-1][nVar-1] = enthalpy;
  for (iDim = 0; iDim < nDim; iDim ++)
    preconditioner[nVar-1][1+iDim] = -1.0*U_i[iDim+1]/rho*enthalpy;


  for (iVar = 0; iVar < nVar; iVar ++ ) {
    for (jVar = 0; jVar < nVar; jVar ++ ) {
      preconditioner[iVar][jVar] = (parameter - 1.0) * ((Gamma-1.0)/(soundspeed*soundspeed))*preconditioner[iVar][jVar];
      if (iVar == jVar)
        preconditioner[iVar][iVar] += 1.0;

      preconditioner[iVar][jVar] *= delta;
    }
  }

}

void CEulerSolver::SetFarfield_AoA(CGeometry *geometry, CSolver **solver_container,
                                   CConfig *config, unsigned short iMesh, bool Output) {

  const auto InnerIter = config->GetInnerIter();
  const su2double AoS = config->GetAoS()*PI_NUMBER/180.0;

  /* --- Initialize values at first iteration --- */

  if (InnerIter == 0) {
    Total_CD_Prev = 0.0;
    Total_CL_Prev = 0.0;
    Total_CMx_Prev = 0.0;
    Total_CMy_Prev = 0.0;
    Total_CMz_Prev = 0.0;
    AoA_Prev = config->GetAoA();
    dCL_dAlpha = config->GetdCL_dAlpha();
    AoA_inc = 0.0;
  }

  /*--- Retrieve the AoA (degrees) ---*/

  su2double AoA = config->GetAoA();

  /* --- Set new AoA if needed --- */

  if (fabs(AoA_inc) > 0.0 && Output) {

    /* --- Update *_Prev values with current coefficients --- */

    SetCoefficient_Gradients(config);

    Total_CD_Prev = TotalCoeff.CD;
    Total_CL_Prev = TotalCoeff.CL;
    Total_CMx_Prev = TotalCoeff.CMx;
    Total_CMy_Prev = TotalCoeff.CMy;
    Total_CMz_Prev = TotalCoeff.CMz;
    AoA_Prev = AoA;

    /*--- Compute a new value for AoA on the fine mesh only (degrees)---*/

    if (iMesh == MESH_0) {
      AoA = AoA + AoA_inc;
      config->SetAoA(AoA);
    }

    AoA *= PI_NUMBER/180.0;

    /*--- Update the freestream velocity vector at the farfield
     * Compute the new freestream velocity with the updated AoA,
     * "Velocity_Inf" is shared with config. ---*/

    const su2double Vel_Infty_Mag = GeometryToolbox::Norm(nDim, Velocity_Inf);

    if (nDim == 2) {
      Velocity_Inf[0] = cos(AoA)*Vel_Infty_Mag;
      Velocity_Inf[1] = sin(AoA)*Vel_Infty_Mag;
    }
    else {
      Velocity_Inf[0] = cos(AoA)*cos(AoS)*Vel_Infty_Mag;
      Velocity_Inf[1] = sin(AoS)*Vel_Infty_Mag;
      Velocity_Inf[2] = sin(AoA)*cos(AoS)*Vel_Infty_Mag;
    }
  }
}

bool CEulerSolver::FixedCL_Convergence(CConfig* config, bool convergence) {
  const su2double Target_CL = config->GetTarget_CL();
  const auto curr_iter = config->GetInnerIter();
  const auto Iter_dCL_dAlpha = config->GetIter_dCL_dAlpha();
  bool fixed_cl_conv = false;
  AoA_inc = 0.0;

  /*--- if in Fixed CL mode, before finite differencing --- */

  if (!Start_AoA_FD){
    if (convergence){

      /* --- C_L and solution are converged, start finite differencing --- */

      if (fabs(TotalCoeff.CL-Target_CL) < (config->GetCauchy_Eps()/2)) {

        /* --- If no finite differencing required --- */

        if (Iter_dCL_dAlpha == 0){
          fixed_cl_conv = true;
          return fixed_cl_conv;
        }

        /* --- Else, set up finite differencing routine ---*/

        Iter_Update_AoA = curr_iter;
        Start_AoA_FD = true;
        fixed_cl_conv = false;
        AoA_inc = 0.001;
      }

      /* --- C_L is not converged to target value and some iterations
          have passed since last update, so update AoA --- */

      else if ((curr_iter - Iter_Update_AoA) > config->GetStartConv_Iter()){
        Iter_Update_AoA = curr_iter;
        fixed_cl_conv = false;
        if (fabs(TotalCoeff.CL-Target_CL) > (config->GetCauchy_Eps()/2)) {
          AoA_inc = (1.0/dCL_dAlpha)*(Target_CL - TotalCoeff.CL);
        }
      }
    }

    /* --- If the iteration limit between AoA updates is met, so update AoA --- */

    else if ((curr_iter - Iter_Update_AoA) == config->GetUpdate_AoA_Iter_Limit()) {
      Iter_Update_AoA = curr_iter;
      fixed_cl_conv = false;
      if (fabs(TotalCoeff.CL-Target_CL) > (config->GetCauchy_Eps()/2)) {
        AoA_inc = (1.0/dCL_dAlpha)*(Target_CL - TotalCoeff.CL);
      }
    }

    /* --- If the total iteration limit is reached, start finite differencing --- */

    if (curr_iter == config->GetnInner_Iter() - Iter_dCL_dAlpha){
      if (Iter_dCL_dAlpha == 0){
        End_AoA_FD = true;
      }
      Iter_Update_AoA = curr_iter;
      Start_AoA_FD = true;
      fixed_cl_conv = false;
      AoA_inc = 0.001;
    }
  }

  /* --- If Finite Difference Mode has ended, end simulation --- */

  if (End_AoA_FD){
    //fixed_cl_conv = true;
    return true;
  }

  /* --- If starting Finite Difference Mode --- */

  if (Start_AoA_FD){

    /* --- Disable history writing --- */

    config->SetHistory_Wrt_Freq(2, 0);

    /* --- End Finite Difference Mode if iteration limit is reached, so simualtion is converged --- */

    End_AoA_FD = ((curr_iter - Iter_Update_AoA - 2) == Iter_dCL_dAlpha ||
      curr_iter == config->GetnInner_Iter()- 2 );

    if (convergence && (curr_iter - Iter_Update_AoA) > config->GetStartConv_Iter())
      End_AoA_FD = true;

    /* --- If Finite Difference mode is ending, reset AoA and calculate Coefficient Gradients --- */

    if (End_AoA_FD){
      SetCoefficient_Gradients(config);
      config->SetAoA(AoA_Prev);
    }
  }

  return fixed_cl_conv;

}

void CEulerSolver::SetCoefficient_Gradients(CConfig *config) const{

  const su2double AoA = config->GetAoA();

  if (AoA == AoA_Prev) return;

  /*--- Calculate gradients of coefficients w.r.t. CL ---*/

  const su2double dCL = TotalCoeff.CL - Total_CL_Prev;
  const su2double dCL_dAlpha_ = dCL / (AoA - AoA_Prev);
  const su2double dCD_dCL_ = (TotalCoeff.CD-Total_CD_Prev) / dCL;
  const su2double dCMx_dCL_ = (TotalCoeff.CMx-Total_CMx_Prev) / dCL;
  const su2double dCMy_dCL_ = (TotalCoeff.CMy-Total_CMy_Prev) / dCL;
  const su2double dCMz_dCL_ = (TotalCoeff.CMz-Total_CMz_Prev) / dCL;

  /*--- Set the value of the  dOF/dCL in the config file ---*/

  config->SetdCD_dCL(dCD_dCL_);
  config->SetdCMx_dCL(dCMx_dCL_);
  config->SetdCMy_dCL(dCMy_dCL_);
  config->SetdCMz_dCL(dCMz_dCL_);
  config->SetdCL_dAlpha(dCL_dAlpha_);
}

void CEulerSolver::UpdateCustomBoundaryConditions(CGeometry **geometry_container, CConfig *config){

  unsigned short nMGlevel;
  unsigned long iMarker;

  // TODO: Update the fluid boundary conditions for MG
  nMGlevel = config->GetnMGLevels();
  if (nMGlevel > 1) {
    for (iMarker=0; iMarker < nMarker; iMarker++) {
      bool isCustomizable = config->GetMarker_All_PyCustom(iMarker);
      bool isInlet = (config->GetMarker_All_KindBC(iMarker) == INLET_FLOW);
      if (isCustomizable && isInlet)
        SU2_MPI::Error("Custom inlet BCs are not currently compatible with multigrid.", CURRENT_FUNCTION);
    }
  }
}

void CEulerSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  unsigned short iDim;
  unsigned long iVertex, iPoint;

  su2double Area, UnitNormal[MAXNDIM] = {0.0};
  su2double Density, Pressure, Energy,  Velocity[MAXNDIM] = {0.0};
  su2double Density_Bound, Pressure_Bound, Vel_Bound[MAXNDIM] = {0.0};
  su2double Density_Infty, Pressure_Infty, Vel_Infty[MAXNDIM] = {0.0};
  su2double SoundSpeed, Entropy, Velocity2, Vn;
  su2double SoundSpeed_Bound, Entropy_Bound, Vel2_Bound, Vn_Bound;
  su2double SoundSpeed_Infty, Entropy_Infty, Vel2_Infty, Vn_Infty, Qn_Infty;
  su2double RiemannPlus, RiemannMinus;
  su2double *V_infty, *V_domain;

  su2double Gas_Constant     = config->GetGas_ConstantND();

  bool implicit       = config->GetKind_TimeIntScheme() == EULER_IMPLICIT;

  su2double *Normal = new su2double[nDim];

  /*--- Loop over all the vertices on this boundary marker ---*/

  SU2_OMP_FOR_DYN(OMP_MIN_SIZE)
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Allocate the value at the infinity ---*/
    V_infty = GetCharacPrimVar(val_marker, iVertex);

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/

    if (geometry->nodes->GetDomain(iPoint)) {

      /*--- Index of the closest interior node ---*/


      /*--- Normal vector for this vertex (negate for outward convention) ---*/

      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);

      /*--- Retrieve solution at the farfield boundary node ---*/
      V_domain = nodes->GetPrimitive(iPoint);

      /*--- Construct solution state at infinity for compressible flow by
         using Riemann invariants, and then impose a weak boundary condition
         by computing the flux using this new state for U. See CFD texts by
         Hirsch or Blazek for more detail. Adapted from an original
         implementation in the Stanford University multi-block (SUmb) solver
         in the routine bcFarfield.f90 written by Edwin van der Weide,
         last modified 06-12-2005. First, compute the unit normal at the
         boundary nodes. ---*/

      Area = GeometryToolbox::Norm(nDim, Normal);
      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = Normal[iDim]/Area;

      /*--- Store primitive variables (density, velocities, velocity squared,
         energy, pressure, and sound speed) at the boundary node, and set some
         other quantities for clarity. Project the current flow velocity vector
         at this boundary node into the local normal direction, i.e. compute
         v_bound.n.  ---*/

      Density_Bound = V_domain[nDim+2];
      Vel2_Bound = 0.0; Vn_Bound = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Vel_Bound[iDim] = V_domain[iDim+1];
        Vel2_Bound     += Vel_Bound[iDim]*Vel_Bound[iDim];
        Vn_Bound       += Vel_Bound[iDim]*UnitNormal[iDim];
      }
      Pressure_Bound   = nodes->GetPressure(iPoint);
      SoundSpeed_Bound = sqrt(Gamma*Pressure_Bound/Density_Bound);
      Entropy_Bound    = pow(Density_Bound, Gamma)/Pressure_Bound;

      /*--- Store the primitive variable state for the freestream. Project
         the freestream velocity vector into the local normal direction,
         i.e. compute v_infty.n. ---*/

      Density_Infty = GetDensity_Inf();
      Vel2_Infty = 0.0; Vn_Infty = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Vel_Infty[iDim] = GetVelocity_Inf(iDim);
        Vel2_Infty     += Vel_Infty[iDim]*Vel_Infty[iDim];
        Vn_Infty       += Vel_Infty[iDim]*UnitNormal[iDim];
      }
      Pressure_Infty   = GetPressure_Inf();
      SoundSpeed_Infty = sqrt(Gamma*Pressure_Infty/Density_Infty);
      Entropy_Infty    = pow(Density_Infty, Gamma)/Pressure_Infty;

      /*--- Adjust the normal freestream velocity for grid movement ---*/

      Qn_Infty = Vn_Infty;

      /*--- Compute acoustic Riemann invariants: R = u.n +/- 2c/(gamma-1).
         These correspond with the eigenvalues (u+c) and (u-c), respectively,
         which represent the acoustic waves. Positive characteristics are
         incoming, and a physical boundary condition is imposed (freestream
         state). This occurs when either (u.n+c) > 0 or (u.n-c) > 0. Negative
         characteristics are leaving the domain, and numerical boundary
         conditions are required by extrapolating from the interior state
         using the Riemann invariants. This occurs when (u.n+c) < 0 or
         (u.n-c) < 0. Note that grid movement is taken into account when
         checking the sign of the eigenvalue. ---*/

      /*--- Check whether (u.n+c) is greater or less than zero ---*/

      if (Qn_Infty > -SoundSpeed_Infty) {
        /*--- Subsonic inflow or outflow ---*/
        RiemannPlus = Vn_Bound + 2.0*SoundSpeed_Bound/Gamma_Minus_One;
      } else {
        /*--- Supersonic inflow ---*/
        RiemannPlus = Vn_Infty + 2.0*SoundSpeed_Infty/Gamma_Minus_One;
      }

      /*--- Check whether (u.n-c) is greater or less than zero ---*/

      if (Qn_Infty > SoundSpeed_Infty) {
        /*--- Supersonic outflow ---*/
        RiemannMinus = Vn_Bound - 2.0*SoundSpeed_Bound/Gamma_Minus_One;
      } else {
        /*--- Subsonic outflow ---*/
        RiemannMinus = Vn_Infty - 2.0*SoundSpeed_Infty/Gamma_Minus_One;
      }

      /*--- Compute a new value for the local normal velocity and speed of
         sound from the Riemann invariants. ---*/

      Vn = 0.5 * (RiemannPlus + RiemannMinus);
      SoundSpeed = 0.25 * (RiemannPlus - RiemannMinus)*Gamma_Minus_One;

      /*--- Construct the primitive variable state at the boundary for
         computing the flux for the weak boundary condition. The values
         that we choose to construct the solution (boundary or freestream)
         depend on whether we are at an inflow or outflow. At an outflow, we
         choose boundary information (at most one characteristic is incoming),
         while at an inflow, we choose infinity values (at most one
         characteristic is outgoing). ---*/

      if (Qn_Infty > 0.0)   {
        /*--- Outflow conditions ---*/
        for (iDim = 0; iDim < nDim; iDim++)
          Velocity[iDim] = Vel_Bound[iDim] + (Vn-Vn_Bound)*UnitNormal[iDim];
        Entropy = Entropy_Bound;
      } else  {
        /*--- Inflow conditions ---*/
        for (iDim = 0; iDim < nDim; iDim++)
          Velocity[iDim] = Vel_Infty[iDim] + (Vn-Vn_Infty)*UnitNormal[iDim];
        Entropy = Entropy_Infty;
      }

      /*--- Recompute the primitive variables. ---*/

      Density = pow(Entropy*SoundSpeed*SoundSpeed/Gamma,1.0/Gamma_Minus_One);
      Velocity2 = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Velocity2 += Velocity[iDim]*Velocity[iDim];
      }
      Pressure = Density*SoundSpeed*SoundSpeed/Gamma;
      Energy   = Pressure/(Gamma_Minus_One*Density) + 0.5*Velocity2;

      /*--- Store new primitive state for computing the flux. ---*/

      V_infty[0] = Pressure/(Gas_Constant*Density);
      for (iDim = 0; iDim < nDim; iDim++)
        V_infty[iDim+1] = Velocity[iDim];
      V_infty[nDim+1] = Pressure;
      V_infty[nDim+2] = Density;
      V_infty[nDim+3] = Energy + Pressure/Density;



      /*--- Set various quantities in the numerics class ---*/

      conv_numerics->SetPrimitive(V_domain, V_infty);

      /*--- Compute the convective residual using an upwind scheme ---*/

      auto residual = conv_numerics->ComputeResidual(config);

      /*--- Update residual value ---*/

      LinSysRes.AddBlock(iPoint, residual);

      /*--- Convective Jacobian contribution for implicit integration ---*/

      if (implicit)
        Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);

    }
  }
  END_SU2_OMP_FOR

  /*--- Free locally allocated memory ---*/
  delete [] Normal;

}

void CEulerSolver::BC_Riemann(CGeometry *geometry, CSolver **solver_container,
                              CNumerics *conv_numerics, CNumerics *visc_numerics,
                              CConfig *config, unsigned short val_marker) {

  unsigned short iDim, iVar, jVar, kVar;
  unsigned long iVertex, iPoint;
  const su2double *Flow_Dir, *Mach;
  su2double P_Total, T_Total, P_static, T_static, Rho_static, Area, UnitNormal[MAXNDIM];
  su2double *Velocity_b, Velocity2_b, Enthalpy_b, Energy_b, StaticEnergy_b, Density_b, Kappa_b, Chi_b, Pressure_b;
  su2double *Velocity_e, Velocity2_e, VelMag_e, Enthalpy_e, Entropy_e, Energy_e = 0.0, StaticEnthalpy_e, StaticEnergy_e, Density_e = 0.0, Pressure_e;
  su2double *Velocity_i, Velocity2_i, Enthalpy_i, Energy_i, StaticEnergy_i, Density_i, Kappa_i, Chi_i, Pressure_i, SoundSpeed_i;
  su2double ProjVelocity_i;
  su2double **P_Tensor, **invP_Tensor, *Lambda_i, **Jacobian_b, **Jacobian_i, **DubDu, *dw, *u_e, *u_i, *u_b;
  su2double *Residual;
  su2double *S_boundary;

  bool implicit             = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  string Marker_Tag         = config->GetMarker_All_TagBound(val_marker);
  bool gravity = (config->GetGravityForce());

  su2double *Normal, *FlowDirMix, TangVelocity, NormalVelocity;
  Normal = new su2double[nDim];

  Velocity_i = new su2double[nDim];
  Velocity_b = new su2double[nDim];
  Velocity_e = new su2double[nDim];
  FlowDirMix = new su2double[nDim];
  Lambda_i = new su2double[nVar];
  u_i = new su2double[nVar];
  u_e = new su2double[nVar];
  u_b = new su2double[nVar];
  dw = new su2double[nVar];

  Residual = new su2double[nVar];

  S_boundary = new su2double[8];

  P_Tensor = new su2double*[nVar];
  invP_Tensor = new su2double*[nVar];
  Jacobian_i = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++)
  {
    P_Tensor[iVar] = new su2double[nVar];
    invP_Tensor[iVar] = new su2double[nVar];
    Jacobian_i[iVar] = new su2double[nVar];
  }

  /*--- Loop over all the vertices on this boundary marker ---*/

  SU2_OMP_FOR_DYN(OMP_MIN_SIZE)
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    if (geometry->nodes->GetDomain(iPoint)) {


      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);

      Area = GeometryToolbox::Norm(nDim, Normal);

      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = Normal[iDim]/Area;

      /*--- Compute the internal state u_i ---*/
      Velocity2_i = 0;
      for (iDim=0; iDim < nDim; iDim++)
      {
        Velocity_i[iDim] = nodes->GetVelocity(iPoint,iDim);
        Velocity2_i += Velocity_i[iDim]*Velocity_i[iDim];
      }


      Density_i = nodes->GetDensity(iPoint);

      Energy_i = nodes->GetEnergy(iPoint);
      StaticEnergy_i = Energy_i - 0.5*Velocity2_i;

      GetFluidModel()->SetTDState_rhoe(Density_i, StaticEnergy_i);

      Pressure_i = GetFluidModel()->GetPressure();
      Enthalpy_i = Energy_i + Pressure_i/Density_i;

      SoundSpeed_i = GetFluidModel()->GetSoundSpeed();

      Kappa_i = GetFluidModel()->GetdPde_rho() / Density_i;
      Chi_i = GetFluidModel()->GetdPdrho_e() - Kappa_i * StaticEnergy_i;

      ProjVelocity_i = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        ProjVelocity_i += Velocity_i[iDim]*UnitNormal[iDim];

      /*--- Build the external state u_e from boundary data and internal node ---*/

      switch(config->GetKind_Data_Riemann(Marker_Tag))
      {
      case TOTAL_CONDITIONS_PT:

        /*--- Retrieve the specified total conditions for this boundary. ---*/
        if (gravity) P_Total = config->GetRiemann_Var1(Marker_Tag) - geometry->nodes->GetCoord(iPoint, nDim-1)*STANDARD_GRAVITY;/// check in which case is true (only freesurface?)
        else P_Total  = config->GetRiemann_Var1(Marker_Tag);
        T_Total  = config->GetRiemann_Var2(Marker_Tag);
        Flow_Dir = config->GetRiemann_FlowDir(Marker_Tag);

        /*--- Non-dim. the inputs if necessary. ---*/
        P_Total /= config->GetPressure_Ref();
        T_Total /= config->GetTemperature_Ref();

        /* --- Computes the total state --- */
        GetFluidModel()->SetTDState_PT(P_Total, T_Total);
        Enthalpy_e = GetFluidModel()->GetStaticEnergy()+ GetFluidModel()->GetPressure()/GetFluidModel()->GetDensity();
        Entropy_e = GetFluidModel()->GetEntropy();

        /* --- Compute the boundary state u_e --- */
        Velocity2_e = Velocity2_i;
        if (nDim == 2){
          NormalVelocity= -sqrt(Velocity2_e)*Flow_Dir[0];
          TangVelocity= -sqrt(Velocity2_e)*Flow_Dir[1];
          Velocity_e[0]= UnitNormal[0]*NormalVelocity - UnitNormal[1]*TangVelocity;
          Velocity_e[1]= UnitNormal[1]*NormalVelocity + UnitNormal[0]*TangVelocity;
        }else{
          for (iDim = 0; iDim < nDim; iDim++)
            Velocity_e[iDim] = sqrt(Velocity2_e)*Flow_Dir[iDim];
        }
        StaticEnthalpy_e = Enthalpy_e - 0.5 * Velocity2_e;
        GetFluidModel()->SetTDState_hs(StaticEnthalpy_e, Entropy_e);
        Density_e = GetFluidModel()->GetDensity();
        StaticEnergy_e = GetFluidModel()->GetStaticEnergy();
        Energy_e = StaticEnergy_e + 0.5 * Velocity2_e;
        break;

      case STATIC_SUPERSONIC_INFLOW_PT:

        /*--- Retrieve the specified total conditions for this boundary. ---*/
        if (gravity) P_static = config->GetRiemann_Var1(Marker_Tag) - geometry->nodes->GetCoord(iPoint, nDim-1)*STANDARD_GRAVITY;/// check in which case is true (only freesurface?)
        else P_static  = config->GetRiemann_Var1(Marker_Tag);
        T_static  = config->GetRiemann_Var2(Marker_Tag);
        Mach = config->GetRiemann_FlowDir(Marker_Tag);

        /*--- Non-dim. the inputs if necessary. ---*/
        P_static /= config->GetPressure_Ref();
        T_static /= config->GetTemperature_Ref();

        /* --- Computes the total state --- */
        GetFluidModel()->SetTDState_PT(P_static, T_static);

        /* --- Compute the boundary state u_e --- */
        Velocity2_e = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          Velocity_e[iDim] = Mach[iDim]*GetFluidModel()->GetSoundSpeed();
          Velocity2_e += Velocity_e[iDim]*Velocity_e[iDim];
        }
        Density_e = GetFluidModel()->GetDensity();
        StaticEnergy_e = GetFluidModel()->GetStaticEnergy();
        Energy_e = StaticEnergy_e + 0.5 * Velocity2_e;
        break;

      case STATIC_SUPERSONIC_INFLOW_PD:

        /*--- Retrieve the specified total conditions for this boundary. ---*/

        if (gravity) P_static = config->GetRiemann_Var1(Marker_Tag) - geometry->nodes->GetCoord(iPoint, nDim-1)*STANDARD_GRAVITY;/// check in which case is true (only freesurface?)
        else P_static  = config->GetRiemann_Var1(Marker_Tag);
        Rho_static  = config->GetRiemann_Var2(Marker_Tag);
        Mach = config->GetRiemann_FlowDir(Marker_Tag);

        /*--- Non-dim. the inputs if necessary. ---*/
        P_static /= config->GetPressure_Ref();
        Rho_static /= config->GetDensity_Ref();

        /* --- Computes the total state --- */
        GetFluidModel()->SetTDState_Prho(P_static, Rho_static);

        /* --- Compute the boundary state u_e --- */
        Velocity2_e = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          Velocity_e[iDim] = Mach[iDim]*GetFluidModel()->GetSoundSpeed();
          Velocity2_e += Velocity_e[iDim]*Velocity_e[iDim];
        }
        Density_e = GetFluidModel()->GetDensity();
        StaticEnergy_e = GetFluidModel()->GetStaticEnergy();
        Energy_e = StaticEnergy_e + 0.5 * Velocity2_e;
        break;

      case DENSITY_VELOCITY:

        /*--- Retrieve the specified density and velocity magnitude ---*/
        Density_e  = config->GetRiemann_Var1(Marker_Tag);
        VelMag_e   = config->GetRiemann_Var2(Marker_Tag);
        Flow_Dir = config->GetRiemann_FlowDir(Marker_Tag);

        /*--- Non-dim. the inputs if necessary. ---*/
        Density_e /= config->GetDensity_Ref();
        VelMag_e /= config->GetVelocity_Ref();

        for (iDim = 0; iDim < nDim; iDim++)
          Velocity_e[iDim] = VelMag_e*Flow_Dir[iDim];
        Energy_e = Energy_i;
        break;

      case STATIC_PRESSURE:

        /*--- Retrieve the static pressure for this boundary. ---*/
        Pressure_e = config->GetRiemann_Var1(Marker_Tag);
        Pressure_e /= config->GetPressure_Ref();
        Density_e = Density_i;

        /* --- Compute the boundary state u_e --- */
        GetFluidModel()->SetTDState_Prho(Pressure_e, Density_e);
        Velocity2_e = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          Velocity_e[iDim] = Velocity_i[iDim];
          Velocity2_e += Velocity_e[iDim]*Velocity_e[iDim];
        }
        Energy_e = GetFluidModel()->GetStaticEnergy() + 0.5*Velocity2_e;
        break;

      default:
        SU2_MPI::Error("Invalid Riemann input!", CURRENT_FUNCTION);
        break;
      }

      /*--- Compute P (matrix of right eigenvectors) ---*/
      conv_numerics->GetPMatrix(&Density_i, Velocity_i, &SoundSpeed_i, &Enthalpy_i, &Chi_i, &Kappa_i, UnitNormal, P_Tensor);

      /*--- Compute inverse P (matrix of left eigenvectors)---*/
      conv_numerics->GetPMatrix_inv(invP_Tensor, &Density_i, Velocity_i, &SoundSpeed_i, &Chi_i, &Kappa_i, UnitNormal);

      /*--- Flow eigenvalues ---*/
      for (iDim = 0; iDim < nDim; iDim++)
        Lambda_i[iDim] = ProjVelocity_i;
      Lambda_i[nVar-2] = ProjVelocity_i + SoundSpeed_i;
      Lambda_i[nVar-1] = ProjVelocity_i - SoundSpeed_i;

      /*--- Compute the boundary state u_e ---*/
      u_e[0] = Density_e;
      for (iDim = 0; iDim < nDim; iDim++)
        u_e[iDim+1] = Velocity_e[iDim]*Density_e;
      u_e[nVar-1] = Energy_e*Density_e;

      /*--- Compute the boundary state u_i ---*/
      u_i[0] = Density_i;
      for (iDim = 0; iDim < nDim; iDim++)
        u_i[iDim+1] = Velocity_i[iDim]*Density_i;
      u_i[nVar-1] = Energy_i*Density_i;

      /*--- Compute the characteristic jumps ---*/
      for (iVar = 0; iVar < nVar; iVar++)
      {
        dw[iVar] = 0;
        for (jVar = 0; jVar < nVar; jVar++)
          dw[iVar] += invP_Tensor[iVar][jVar] * (u_e[jVar] - u_i[jVar]);

      }

      /*--- Compute the boundary state u_b using characteristics ---*/
      for (iVar = 0; iVar < nVar; iVar++)
      {
        u_b[iVar] = u_i[iVar];

        for (jVar = 0; jVar < nVar; jVar++)
        {
          if (Lambda_i[jVar] < 0)
          {
            u_b[iVar] += P_Tensor[iVar][jVar]*dw[jVar];

          }
        }
      }


      /*--- Compute the thermodynamic state in u_b ---*/
      Density_b = u_b[0];
      Velocity2_b = 0;
      for (iDim = 0; iDim < nDim; iDim++)
      {
        Velocity_b[iDim] = u_b[iDim+1]/Density_b;
        Velocity2_b += Velocity_b[iDim]*Velocity_b[iDim];
      }
      Energy_b = u_b[nVar-1]/Density_b;
      StaticEnergy_b = Energy_b - 0.5*Velocity2_b;
      GetFluidModel()->SetTDState_rhoe(Density_b, StaticEnergy_b);
      Pressure_b = GetFluidModel()->GetPressure();
      Enthalpy_b = Energy_b + Pressure_b/Density_b;
      Kappa_b = GetFluidModel()->GetdPde_rho() / Density_b;
      Chi_b = GetFluidModel()->GetdPdrho_e() - Kappa_b * StaticEnergy_b;

      /*--- Compute the residuals ---*/
      conv_numerics->GetInviscidProjFlux(&Density_b, Velocity_b, &Pressure_b, &Enthalpy_b, Normal, Residual);

      if (implicit) {

        Jacobian_b = new su2double*[nVar];
        DubDu = new su2double*[nVar];
        for (iVar = 0; iVar < nVar; iVar++)
        {
          Jacobian_b[iVar] = new su2double[nVar];
          DubDu[iVar] = new su2double[nVar];
        }

        /*--- Initialize DubDu to unit matrix---*/

        for (iVar = 0; iVar < nVar; iVar++)
        {
          for (jVar = 0; jVar < nVar; jVar++)
            DubDu[iVar][jVar]= 0;

          DubDu[iVar][iVar]= 1;
        }

        /*--- Compute DubDu -= RNL---*/
        for (iVar=0; iVar<nVar; iVar++)
        {
          for (jVar=0; jVar<nVar; jVar++)
          {
            for (kVar=0; kVar<nVar; kVar++)
            {
              if (Lambda_i[kVar]<0)
                DubDu[iVar][jVar] -= P_Tensor[iVar][kVar] * invP_Tensor[kVar][jVar];
            }
          }
        }

        /*--- Compute flux Jacobian in state b ---*/
        conv_numerics->GetInviscidProjJac(Velocity_b, &Enthalpy_b, &Chi_b, &Kappa_b, Normal, 1.0, Jacobian_b);


        /*--- initiate Jacobian_i to zero matrix ---*/
        for (iVar=0; iVar<nVar; iVar++)
          for (jVar=0; jVar<nVar; jVar++)
            Jacobian_i[iVar][jVar] = 0.0;

        /*--- Compute numerical flux Jacobian at node i ---*/
        for (iVar=0; iVar<nVar; iVar++) {
          for (jVar=0; jVar<nVar; jVar++) {
            for (kVar=0; kVar<nVar; kVar++) {
              Jacobian_i[iVar][jVar] += Jacobian_b[iVar][kVar] * DubDu[kVar][jVar];
            }
          }
        }

        for (iVar = 0; iVar < nVar; iVar++) {
          delete [] Jacobian_b[iVar];
          delete [] DubDu[iVar];
        }
        delete [] Jacobian_b;
        delete [] DubDu;
      }

      /*--- Update residual value ---*/
      LinSysRes.AddBlock(iPoint, Residual);

      /*--- Jacobian contribution for implicit integration ---*/
      if (implicit)
        Jacobian.AddBlock2Diag(iPoint, Jacobian_i);

    }
  }
  END_SU2_OMP_FOR

  /*--- Free locally allocated memory ---*/
  delete [] Normal;
  delete [] Velocity_e;
  delete [] Velocity_b;
  delete [] Velocity_i;
  delete [] FlowDirMix;

  delete [] S_boundary;
  delete [] Lambda_i;
  delete [] u_i;
  delete [] u_e;
  delete [] u_b;
  delete [] dw;

  delete [] Residual;

  for (iVar = 0; iVar < nVar; iVar++)
  {
    delete [] P_Tensor[iVar];
    delete [] invP_Tensor[iVar];
    delete [] Jacobian_i[iVar];
  }
  delete [] P_Tensor;
  delete [] invP_Tensor;
  delete [] Jacobian_i;

}

void CEulerSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container,
                            CNumerics *conv_numerics, CNumerics *visc_numerics,
                            CConfig *config, unsigned short val_marker) {
  unsigned short iDim;
  unsigned long iVertex, iPoint;
  su2double P_Total, T_Total, Velocity[3], Velocity2, H_Total, Temperature, Riemann,
  Pressure, Density, Energy, *Flow_Dir, Mach2, SoundSpeed2, SoundSpeed_Total2, Vel_Mag,
  alpha, aa, bb, cc, dd, Area, UnitNormal[3];
  su2double *V_inlet, *V_domain;

  bool implicit             = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  su2double Two_Gamma_M1       = 2.0/Gamma_Minus_One;
  su2double Gas_Constant       = config->GetGas_ConstantND();
  INLET_TYPE Kind_Inlet= config->GetKind_Inlet();
  string Marker_Tag         = config->GetMarker_All_TagBound(val_marker);
  su2double *Normal = new su2double[nDim];

  /*--- Loop over all the vertices on this boundary marker ---*/

  SU2_OMP_FOR_DYN(OMP_MIN_SIZE)
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    /*--- Allocate the value at the inlet ---*/

    V_inlet = GetCharacPrimVar(val_marker, iVertex);

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/

    if (geometry->nodes->GetDomain(iPoint)) {

      /*--- Normal vector for this vertex (negate for outward convention) ---*/

      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);

      Area = GeometryToolbox::Norm(nDim, Normal);
      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = Normal[iDim]/Area;

      /*--- Retrieve solution at this boundary node ---*/

      V_domain = nodes->GetPrimitive(iPoint);

      /*--- Build the fictitious intlet state based on characteristics ---*/


      /*--- Subsonic inflow: there is one outgoing characteristic (u-c),
         therefore we can specify all but one state variable at the inlet.
         The outgoing Riemann invariant provides the final piece of info.
         Adapted from an original implementation in the Stanford University
         multi-block (SUmb) solver in the routine bcSubsonicInflow.f90
         written by Edwin van der Weide, last modified 04-20-2009. ---*/

      switch (Kind_Inlet) {

        /*--- Total properties have been specified at the inlet. ---*/

        case INLET_TYPE::TOTAL_CONDITIONS:

          /*--- Retrieve the specified total conditions for this inlet. ---*/

          P_Total  = Inlet_Ptotal[val_marker][iVertex];
          T_Total  = Inlet_Ttotal[val_marker][iVertex];
          Flow_Dir = Inlet_FlowDir[val_marker][iVertex];

          /*--- Non-dim. the inputs if necessary. ---*/

          P_Total /= config->GetPressure_Ref();
          T_Total /= config->GetTemperature_Ref();

          /*--- Store primitives and set some variables for clarity. ---*/

          Density = V_domain[nDim+2];
          Velocity2 = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            Velocity[iDim] = V_domain[iDim+1];
            Velocity2 += Velocity[iDim]*Velocity[iDim];
          }
          Energy      = V_domain[nDim+3] - V_domain[nDim+1]/V_domain[nDim+2];
          Pressure    = V_domain[nDim+1];
          H_Total     = (Gamma*Gas_Constant/Gamma_Minus_One)*T_Total;
          SoundSpeed2 = Gamma*Pressure/Density;

          /*--- Compute the acoustic Riemann invariant that is extrapolated
             from the domain interior. ---*/

          Riemann   = 2.0*sqrt(SoundSpeed2)/Gamma_Minus_One;
          for (iDim = 0; iDim < nDim; iDim++)
            Riemann += Velocity[iDim]*UnitNormal[iDim];

          /*--- Total speed of sound ---*/

          SoundSpeed_Total2 = Gamma_Minus_One*(H_Total - (Energy + Pressure/Density)+0.5*Velocity2) + SoundSpeed2;

          /*--- Dot product of normal and flow direction. This should
             be negative due to outward facing boundary normal convention. ---*/

          alpha = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            alpha += UnitNormal[iDim]*Flow_Dir[iDim];

          /*--- Coefficients in the quadratic equation for the velocity ---*/

          aa =  1.0 + 0.5*Gamma_Minus_One*alpha*alpha;
          bb = -1.0*Gamma_Minus_One*alpha*Riemann;
          cc =  0.5*Gamma_Minus_One*Riemann*Riemann
              -2.0*SoundSpeed_Total2/Gamma_Minus_One;

          /*--- Solve quadratic equation for velocity magnitude. Value must
             be positive, so the choice of root is clear. ---*/

          dd = bb*bb - 4.0*aa*cc;
          dd = sqrt(max(0.0, dd));
          Vel_Mag   = (-bb + dd)/(2.0*aa);
          Vel_Mag   = max(0.0, Vel_Mag);
          Velocity2 = Vel_Mag*Vel_Mag;

          /*--- Compute speed of sound from total speed of sound eqn. ---*/

          SoundSpeed2 = SoundSpeed_Total2 - 0.5*Gamma_Minus_One*Velocity2;

          /*--- Mach squared (cut between 0-1), use to adapt velocity ---*/

          Mach2 = Velocity2/SoundSpeed2;
          Mach2 = min(1.0, Mach2);
          Velocity2   = Mach2*SoundSpeed2;
          Vel_Mag     = sqrt(Velocity2);
          SoundSpeed2 = SoundSpeed_Total2 - 0.5*Gamma_Minus_One*Velocity2;

          /*--- Compute new velocity vector at the inlet ---*/

          for (iDim = 0; iDim < nDim; iDim++)
            Velocity[iDim] = Vel_Mag*Flow_Dir[iDim];

          /*--- Static temperature from the speed of sound relation ---*/

          Temperature = SoundSpeed2/(Gamma*Gas_Constant);

          /*--- Static pressure using isentropic relation at a point ---*/

          Pressure = P_Total*pow((Temperature/T_Total), Gamma/Gamma_Minus_One);

          /*--- Density at the inlet from the gas law ---*/

          Density = Pressure/(Gas_Constant*Temperature);

          /*--- Using pressure, density, & velocity, compute the energy ---*/

          Energy = Pressure/(Density*Gamma_Minus_One) + 0.5*Velocity2;

          /*--- Primitive variables, using the derived quantities ---*/

          V_inlet[0] = Temperature;
          for (iDim = 0; iDim < nDim; iDim++)
            V_inlet[iDim+1] = Velocity[iDim];
          V_inlet[nDim+1] = Pressure;
          V_inlet[nDim+2] = Density;
          V_inlet[nDim+3] = Energy + Pressure/Density;

          break;

          /*--- Mass flow has been specified at the inlet. ---*/

        case INLET_TYPE::MASS_FLOW:

          /*--- Retrieve the specified mass flow for the inlet. ---*/

          Density  = Inlet_Ttotal[val_marker][iVertex];
          Vel_Mag  = Inlet_Ptotal[val_marker][iVertex];
          Flow_Dir = Inlet_FlowDir[val_marker][iVertex];

          /*--- Non-dim. the inputs if necessary. ---*/

          Density /= config->GetDensity_Ref();
          Vel_Mag /= config->GetVelocity_Ref();

          /*--- Get primitives from current inlet state. ---*/

          for (iDim = 0; iDim < nDim; iDim++)
            Velocity[iDim] = nodes->GetVelocity(iPoint,iDim);
          Pressure    = nodes->GetPressure(iPoint);
          SoundSpeed2 = Gamma*Pressure/V_domain[nDim+2];

          /*--- Compute the acoustic Riemann invariant that is extrapolated
             from the domain interior. ---*/

          Riemann = Two_Gamma_M1*sqrt(SoundSpeed2);
          for (iDim = 0; iDim < nDim; iDim++)
            Riemann += Velocity[iDim]*UnitNormal[iDim];

          /*--- Speed of sound squared for fictitious inlet state ---*/

          SoundSpeed2 = Riemann;
          for (iDim = 0; iDim < nDim; iDim++)
            SoundSpeed2 -= Vel_Mag*Flow_Dir[iDim]*UnitNormal[iDim];

          SoundSpeed2 = max(0.0,0.5*Gamma_Minus_One*SoundSpeed2);
          SoundSpeed2 = SoundSpeed2*SoundSpeed2;

          /*--- Pressure for the fictitious inlet state ---*/

          Pressure = SoundSpeed2*Density/Gamma;

          /*--- Energy for the fictitious inlet state ---*/

          Energy = Pressure/(Density*Gamma_Minus_One) + 0.5*Vel_Mag*Vel_Mag;

          /*--- Primitive variables, using the derived quantities ---*/

          V_inlet[0] = Pressure / ( Gas_Constant * Density);
          for (iDim = 0; iDim < nDim; iDim++)
            V_inlet[iDim+1] = Vel_Mag*Flow_Dir[iDim];
          V_inlet[nDim+1] = Pressure;
          V_inlet[nDim+2] = Density;
          V_inlet[nDim+3] = Energy + Pressure/Density;

          break;

        default:
          SU2_MPI::Error("Unsupported INLET_TYPE.", CURRENT_FUNCTION);
          break;
      }

      /*--- Set various quantities in the solver class ---*/

      conv_numerics->SetPrimitive(V_domain, V_inlet);


      /*--- Compute the residual using an upwind scheme ---*/

      auto residual = conv_numerics->ComputeResidual(config);

      /*--- Update residual value ---*/

      LinSysRes.AddBlock(iPoint, residual);

      /*--- Jacobian contribution for implicit integration ---*/

      if (implicit)
        Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);

    }
  }
  END_SU2_OMP_FOR

  /*--- Free locally allocated memory ---*/

  delete [] Normal;

}

void CEulerSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container,
                             CNumerics *conv_numerics, CNumerics *visc_numerics,
                             CConfig *config, unsigned short val_marker) {
  unsigned short iVar, iDim;
  unsigned long iVertex, iPoint;
  su2double Pressure, P_Exit, Velocity[3],
  Velocity2, Entropy, Density, Energy, Riemann, Vn, SoundSpeed, Mach_Exit, Vn_Exit,
  Area, UnitNormal[3];
  su2double *V_outlet, *V_domain;

  bool implicit           = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  su2double Gas_Constant     = config->GetGas_ConstantND();
  string Marker_Tag       = config->GetMarker_All_TagBound(val_marker);
  bool gravity = (config->GetGravityForce());
  su2double *Normal = new su2double[nDim];

  /*--- Loop over all the vertices on this boundary marker ---*/

  SU2_OMP_FOR_DYN(OMP_MIN_SIZE)
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    /*--- Allocate the value at the outlet ---*/
    V_outlet = GetCharacPrimVar(val_marker, iVertex);

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    if (geometry->nodes->GetDomain(iPoint)) {

      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);

      Area = GeometryToolbox::Norm(nDim, Normal);
      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = Normal[iDim]/Area;

      /*--- Current solution at this boundary node ---*/
      V_domain = nodes->GetPrimitive(iPoint);

      /*--- Build the fictitious inlet state based on characteristics ---*/

      /*--- Retrieve the specified back pressure for this outlet. ---*/
      if (gravity) P_Exit = config->GetOutlet_Pressure(Marker_Tag) - geometry->nodes->GetCoord(iPoint, nDim-1)*STANDARD_GRAVITY;
      else P_Exit = config->GetOutlet_Pressure(Marker_Tag);

      /*--- Non-dim. the inputs if necessary. ---*/
      P_Exit = P_Exit/config->GetPressure_Ref();

      /*--- Check whether the flow is supersonic at the exit. The type
         of boundary update depends on this. ---*/
      Density = V_domain[nDim+2];
      Velocity2 = 0.0; Vn = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Velocity[iDim] = V_domain[iDim+1];
        Velocity2 += Velocity[iDim]*Velocity[iDim];
        Vn += Velocity[iDim]*UnitNormal[iDim];
      }
      Pressure   = V_domain[nDim+1];
      SoundSpeed = sqrt(Gamma*Pressure/Density);
      Mach_Exit  = sqrt(Velocity2)/SoundSpeed;

      if (Mach_Exit >= 1.0) {

        /*--- Supersonic exit flow: there are no incoming characteristics,
           so no boundary condition is necessary. Set outlet state to current
           state so that upwinding handles the direction of propagation. ---*/
        for (iVar = 0; iVar < nPrimVar; iVar++) V_outlet[iVar] = V_domain[iVar];

      } else {

        /*--- Subsonic exit flow: there is one incoming characteristic,
           therefore one variable can be specified (back pressure) and is used
           to update the conservative variables. Compute the entropy and the
           acoustic Riemann variable. These invariants, as well as the
           tangential velocity components, are extrapolated. Adapted from an
           original implementation in the Stanford University multi-block
           (SUmb) solver in the routine bcSubsonicOutflow.f90 by Edwin van
           der Weide, last modified 09-10-2007. ---*/

        Entropy = Pressure*pow(1.0/Density, Gamma);
        Riemann = Vn + 2.0*SoundSpeed/Gamma_Minus_One;

        /*--- Compute the new fictious state at the outlet ---*/
        Density    = pow(P_Exit/Entropy,1.0/Gamma);
        Pressure   = P_Exit;
        SoundSpeed = sqrt(Gamma*P_Exit/Density);
        Vn_Exit    = Riemann - 2.0*SoundSpeed/Gamma_Minus_One;
        Velocity2  = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          Velocity[iDim] = Velocity[iDim] + (Vn_Exit-Vn)*UnitNormal[iDim];
          Velocity2 += Velocity[iDim]*Velocity[iDim];
        }
        Energy = P_Exit/(Density*Gamma_Minus_One) + 0.5*Velocity2;

        /*--- Conservative variables, using the derived quantities ---*/
        V_outlet[0] = Pressure / ( Gas_Constant * Density);
        for (iDim = 0; iDim < nDim; iDim++)
          V_outlet[iDim+1] = Velocity[iDim];
        V_outlet[nDim+1] = Pressure;
        V_outlet[nDim+2] = Density;
        V_outlet[nDim+3] = Energy + Pressure/Density;

      }

      /*--- Set various quantities in the solver class ---*/
      conv_numerics->SetPrimitive(V_domain, V_outlet);


      /*--- Compute the residual using an upwind scheme ---*/

      auto residual = conv_numerics->ComputeResidual(config);

      /*--- Add Residuals and Jacobians ---*/

      LinSysRes.AddBlock(iPoint, residual);
      if (implicit)
        Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);

    }
  }
  END_SU2_OMP_FOR

  /*--- Free locally allocated memory ---*/
  delete [] Normal;

}

void CEulerSolver::BC_Supersonic_Inlet(CGeometry *geometry, CSolver **solver_container,
                                       CNumerics *conv_numerics, CNumerics *visc_numerics,
                                       CConfig *config, unsigned short val_marker) {
  unsigned short iDim;
  unsigned long iVertex, iPoint;
  su2double *V_inlet, *V_domain;

  su2double Density, Energy, Velocity2;
  su2double Gas_Constant = config->GetGas_ConstantND();

  bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  su2double *Normal = new su2double[nDim];
  su2double *Velocity = new su2double[nDim];

  /*--- Supersonic inlet flow: there are no outgoing characteristics,
   so all flow variables can be imposed at the inlet.
   First, retrieve the specified values for the primitive variables. ---*/

  auto Temperature = config->GetInlet_Temperature(Marker_Tag);
  auto Pressure    = config->GetInlet_Pressure(Marker_Tag);
  auto Vel         = config->GetInlet_Velocity(Marker_Tag);

  /*--- Non-dim. the inputs if necessary. ---*/

  Temperature /= config->GetTemperature_Ref();
  Pressure    /= config->GetPressure_Ref();
  for (iDim = 0; iDim < nDim; iDim++)
    Velocity[iDim] = Vel[iDim] / config->GetVelocity_Ref();

  /*--- Density at the inlet from the gas law ---*/

  Density = Pressure/(Gas_Constant*Temperature);

  /*--- Compute the energy from the specified state ---*/

  Velocity2 = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Velocity2 += Velocity[iDim]*Velocity[iDim];
  Energy = Pressure/(Density*Gamma_Minus_One)+0.5*Velocity2;

  /*--- Loop over all the vertices on this boundary marker ---*/

  SU2_OMP_FOR_DYN(OMP_MIN_SIZE)
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    /*--- Allocate the value at the outlet ---*/

    V_inlet = GetCharacPrimVar(val_marker, iVertex);

    /*--- Primitive variables, using the derived quantities ---*/

    V_inlet[0] = Temperature;
    for (iDim = 0; iDim < nDim; iDim++)
      V_inlet[iDim+1] = Velocity[iDim];
    V_inlet[nDim+1] = Pressure;
    V_inlet[nDim+2] = Density;
    V_inlet[nDim+3] = Energy + Pressure/Density;

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/

    if (geometry->nodes->GetDomain(iPoint)) {

      /*--- Current solution at this boundary node ---*/

      V_domain = nodes->GetPrimitive(iPoint);

      /*--- Normal vector for this vertex (negate for outward convention) ---*/

      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];

      /*--- Set various quantities in the solver class ---*/

      conv_numerics->SetNormal(Normal);
      conv_numerics->SetPrimitive(V_domain, V_inlet);


      /*--- Compute the residual using an upwind scheme ---*/

      auto residual = conv_numerics->ComputeResidual(config);

      LinSysRes.AddBlock(iPoint, residual);

      /*--- Jacobian contribution for implicit integration ---*/

      if (implicit)
        Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);

    }
  }
  END_SU2_OMP_FOR

  /*--- Free locally allocated memory ---*/

  delete [] Normal;
  delete [] Velocity;

}

void CEulerSolver::BC_Supersonic_Outlet(CGeometry *geometry, CSolver **solver_container,
                                        CNumerics *conv_numerics, CNumerics *visc_numerics,
                                        CConfig *config, unsigned short val_marker) {
  unsigned short iDim;
  unsigned long iVertex, iPoint;
  su2double *V_outlet, *V_domain;

  bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  su2double *Normal = new su2double[nDim];

  /*--- Supersonic outlet flow: there are no ingoing characteristics,
   so all flow variables can should be interpolated from the domain. ---*/

  /*--- Loop over all the vertices on this boundary marker ---*/

  SU2_OMP_FOR_DYN(OMP_MIN_SIZE)
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/

    if (geometry->nodes->GetDomain(iPoint)) {

      /*--- Current solution at this boundary node ---*/

      V_domain = nodes->GetPrimitive(iPoint);

      /*--- Allocate the value at the outlet ---*/

      V_outlet = GetCharacPrimVar(val_marker, iVertex);

      /*--- Primitive variables, using the derived quantities ---*/

      V_outlet[0] = V_domain[0];
      for (iDim = 0; iDim < nDim; iDim++)
        V_outlet[iDim+1] = V_domain[iDim+1];
      V_outlet[nDim+1] = V_domain[nDim+1];
      V_outlet[nDim+2] = V_domain[nDim+2];
      V_outlet[nDim+3] = V_domain[nDim+3];

      /*--- Current solution at this boundary node ---*/

      V_domain = nodes->GetPrimitive(iPoint);

      /*--- Normal vector for this vertex (negate for outward convention) ---*/

      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];

      /*--- Set various quantities in the solver class ---*/

      conv_numerics->SetNormal(Normal);
      conv_numerics->SetPrimitive(V_domain, V_outlet);


      /*--- Compute the residual using an upwind scheme ---*/

      auto residual = conv_numerics->ComputeResidual(config);

      LinSysRes.AddBlock(iPoint, residual);

      /*--- Jacobian contribution for implicit integration ---*/

      if (implicit)
        Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);

    }
  }
  END_SU2_OMP_FOR

  /*--- Free locally allocated memory ---*/

  delete [] Normal;

}

void CEulerSolver::BC_Engine_Inflow(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  unsigned short iDim;
  unsigned long iVertex, iPoint;
  su2double Pressure, Inflow_Pressure = 0.0, Velocity[3], Velocity2, Entropy, Target_Inflow_MassFlow = 0.0, Target_Inflow_Mach = 0.0, Density, Energy,
  Riemann, Area, UnitNormal[3], Vn, SoundSpeed, Vn_Exit, Inflow_Pressure_inc, Inflow_Pressure_old, Inflow_Mach_old, Inflow_MassFlow_old;
  su2double *V_inflow, *V_domain;

  su2double DampingFactor = config->GetDamp_Engine_Inflow();
  bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  unsigned short Kind_Engine_Inflow = config->GetKind_Engine_Inflow();
  su2double Gas_Constant = config->GetGas_ConstantND();
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  su2double Baseline_Press = 0.75 * config->GetPressure_FreeStreamND();
  bool Engine_HalfModel = config->GetEngine_HalfModel();

  su2double *Normal = new su2double[nDim];


  if (Kind_Engine_Inflow == FAN_FACE_MACH) {

    /*--- Retrieve the specified target fan face mach at the nacelle. ---*/

    Target_Inflow_Mach = config->GetEngineInflow_Target(Marker_Tag);

    /*--- Retrieve the old fan face pressure and mach number in the nacelle (this has been computed in a preprocessing). ---*/

    Inflow_Pressure_old = config->GetInflow_Pressure(Marker_Tag);  // Note that has been computed by the code (non-dimensional).
    Inflow_Mach_old = config->GetInflow_Mach(Marker_Tag);

    /*--- Compute the pressure increment (note that increasing pressure decreases flow speed) ---*/

    Inflow_Pressure_inc = - (1.0 - (Inflow_Mach_old/Target_Inflow_Mach)) * Baseline_Press;

    /*--- Estimate the new fan face pressure ---*/

    Inflow_Pressure = (1.0 - DampingFactor)*Inflow_Pressure_old + DampingFactor * (Inflow_Pressure_old + Inflow_Pressure_inc);

  }

  if (Kind_Engine_Inflow == FAN_FACE_MDOT) {

    /*--- Retrieve the specified target mass flow (non-dimensional) at the nacelle. ---*/

    Target_Inflow_MassFlow = config->GetEngineInflow_Target(Marker_Tag) / (config->GetDensity_Ref() * config->GetVelocity_Ref());

    if (config->GetSystemMeasurements() == US) Target_Inflow_MassFlow /= 32.174;

    if (Engine_HalfModel) Target_Inflow_MassFlow /= 2.0;

    /*--- Retrieve the old fan face pressure and mach number in the nacelle (this has been computed in a preprocessing). ---*/

    Inflow_Pressure_old = config->GetInflow_Pressure(Marker_Tag);  // Note that has been computed by the code (non-dimensional).
    Inflow_MassFlow_old = config->GetInflow_MassFlow(Marker_Tag);  // same here... it is a non dimensional value

    /*--- Compute the pressure increment (note that increasing pressure decreases flow speed) ---*/

    Inflow_Pressure_inc = - (1.0 - (Inflow_MassFlow_old/Target_Inflow_MassFlow)) * Baseline_Press;

    /*--- Estimate the new fan face pressure ---*/

    Inflow_Pressure = (1.0 - DampingFactor)*Inflow_Pressure_old + DampingFactor * (Inflow_Pressure_old + Inflow_Pressure_inc);

  }

  /*--- No iterative scheme if we provide the static pressure ---*/

  if (Kind_Engine_Inflow == FAN_FACE_PRESSURE) {

    /*--- Retrieve the specified pressure (non-dimensional) at the nacelle. ---*/

    Inflow_Pressure = config->GetEngineInflow_Target(Marker_Tag) / config->GetPressure_Ref();

  }


  /*--- Loop over all the vertices on this boundary marker ---*/

  SU2_OMP_FOR_DYN(OMP_MIN_SIZE)
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    /*--- Allocate the value at the outlet ---*/

    V_inflow = GetCharacPrimVar(val_marker, iVertex);

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/

    if (geometry->nodes->GetDomain(iPoint)) {

      /*--- Normal vector for this vertex (negate for outward convention) ---*/

      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];

      Area = GeometryToolbox::Norm(nDim, Normal);
      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = Normal[iDim]/Area;

      /*--- Current solution at this boundary node ---*/

      V_domain = nodes->GetPrimitive(iPoint);

      /*--- Subsonic nacelle inflow: there is one incoming characteristic,
       therefore one variable can be specified (back pressure) and is used
       to update the conservative variables.

       Compute the entropy and the acoustic variable. These
       riemann invariants, as well as the tangential velocity components,
       are extrapolated. ---*/

      Density = V_domain[nDim+2];
      Velocity2 = 0.0; Vn = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Velocity[iDim] = V_domain[iDim+1];
        Velocity2 += Velocity[iDim]*Velocity[iDim];
        Vn += Velocity[iDim]*UnitNormal[iDim];
      }
      Pressure   = V_domain[nDim+1];
      SoundSpeed = sqrt(Gamma*Pressure/Density);
      Entropy = Pressure*pow(1.0/Density, Gamma);
      Riemann = Vn + 2.0*SoundSpeed/Gamma_Minus_One;

      /*--- Compute the new fictious state at the outlet ---*/

      Density    = pow(Inflow_Pressure/Entropy,1.0/Gamma);
      Pressure   = Inflow_Pressure;
      SoundSpeed = sqrt(Gamma*Inflow_Pressure/Density);
      Vn_Exit    = Riemann - 2.0*SoundSpeed/Gamma_Minus_One;
      Velocity2  = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Velocity[iDim] = Velocity[iDim] + (Vn_Exit-Vn)*UnitNormal[iDim];
        Velocity2 += Velocity[iDim]*Velocity[iDim];
      }

      Energy = Inflow_Pressure/(Density*Gamma_Minus_One) + 0.5*Velocity2;

      /*--- Conservative variables, using the derived quantities ---*/

      V_inflow[0] = Pressure / ( Gas_Constant * Density);
      for (iDim = 0; iDim < nDim; iDim++)
        V_inflow[iDim+1] = Velocity[iDim];
      V_inflow[nDim+1] = Pressure;
      V_inflow[nDim+2] = Density;
      V_inflow[nDim+3] = Energy + Pressure/Density;
      V_inflow[nDim+4] = SoundSpeed;

      /*--- Set various quantities in the solver class ---*/

      conv_numerics->SetNormal(Normal);
      conv_numerics->SetPrimitive(V_domain, V_inflow);

      /*--- Set grid movement ---*/


      /*--- Compute the residual using an upwind scheme ---*/

      auto residual = conv_numerics->ComputeResidual(config);

      LinSysRes.AddBlock(iPoint, residual);

      /*--- Jacobian contribution for implicit integration ---*/

      if (implicit)
        Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);

    }
  }
  END_SU2_OMP_FOR

  delete [] Normal;

}

void CEulerSolver::BC_Engine_Exhaust(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  unsigned short iDim;
  unsigned long iVertex, iPoint;
  su2double Exhaust_Pressure, Exhaust_Temperature, Velocity[3], Velocity2, H_Exhaust, Temperature, Riemann, Area, UnitNormal[3], Pressure, Density, Energy, Mach2, SoundSpeed2, SoundSpeed_Exhaust2, Vel_Mag, alpha, aa, bb, cc, dd, Flow_Dir[3];
  su2double *V_exhaust, *V_domain, Target_Exhaust_Pressure, Exhaust_Pressure_old, Exhaust_Pressure_inc;

  su2double Gas_Constant = config->GetGas_ConstantND();
  bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  su2double DampingFactor = config->GetDamp_Engine_Exhaust();
  su2double Baseline_Press = 0.75 * config->GetPressure_FreeStreamND();

  su2double *Normal = new su2double[nDim];

  /*--- Retrieve the specified exhaust pressure in the engine (non-dimensional). ---*/

  Target_Exhaust_Pressure = config->GetExhaust_Pressure_Target(Marker_Tag) / config->GetPressure_Ref();

  /*--- Retrieve the old exhaust pressure in the engine exhaust (this has been computed in a preprocessing). ---*/

  Exhaust_Pressure_old = config->GetExhaust_Pressure(Marker_Tag);

  /*--- Compute the Pressure increment ---*/

  Exhaust_Pressure_inc = (1.0 - (Exhaust_Pressure_old/Target_Exhaust_Pressure)) * Baseline_Press;

  /*--- Estimate the new exhaust pressure ---*/

  Exhaust_Pressure = (1.0 - DampingFactor) * Exhaust_Pressure_old + DampingFactor * (Exhaust_Pressure_old + Exhaust_Pressure_inc);

  /*--- The temperature is given (no iteration is required) ---*/

  Exhaust_Temperature  = config->GetExhaust_Temperature_Target(Marker_Tag);
  Exhaust_Temperature /= config->GetTemperature_Ref();

  /*--- The pressure is given (no iteration is required) ---*/

  Exhaust_Pressure  = config->GetExhaust_Pressure_Target(Marker_Tag);
  Exhaust_Pressure /= config->GetPressure_Ref();

  /*--- Loop over all the vertices on this boundary marker ---*/

  SU2_OMP_FOR_DYN(OMP_MIN_SIZE)
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    /*--- Allocate the value at the exhaust ---*/

    V_exhaust = GetCharacPrimVar(val_marker, iVertex);

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/

    if (geometry->nodes->GetDomain(iPoint)) {

      /*--- Normal vector for this vertex (negate for outward convention) ---*/

      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];

      Area = GeometryToolbox::Norm(nDim, Normal);
      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = Normal[iDim]/Area;

      /*--- Current solution at this boundary node ---*/

      V_domain = nodes->GetPrimitive(iPoint);

      /*--- Subsonic inflow: there is one outgoing characteristic (u-c),
       therefore we can specify all but one state variable at the inlet.
       The outgoing Riemann invariant provides the final piece of info. ---*/

      /*--- Store primitives and set some variables for clarity. ---*/

      Density = V_domain[nDim+2];
      Velocity2 = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Velocity[iDim] = V_domain[iDim+1];
        Velocity2 += Velocity[iDim]*Velocity[iDim];
      }
      Energy      = V_domain[nDim+3] - V_domain[nDim+1]/V_domain[nDim+2];
      Pressure    = V_domain[nDim+1];
      H_Exhaust   = (Gamma*Gas_Constant/Gamma_Minus_One)*Exhaust_Temperature;
      SoundSpeed2 = Gamma*Pressure/Density;

      /*--- Compute the acoustic Riemann invariant that is extrapolated
       from the domain interior. ---*/

      Riemann   = 2.0*sqrt(SoundSpeed2)/Gamma_Minus_One;
      for (iDim = 0; iDim < nDim; iDim++)
        Riemann += Velocity[iDim]*UnitNormal[iDim];

      /*--- Total speed of sound ---*/

      SoundSpeed_Exhaust2 = Gamma_Minus_One*(H_Exhaust - (Energy + Pressure/Density)+0.5*Velocity2) + SoundSpeed2;

      /*--- The flow direction is defined by the surface normal ---*/

      for (iDim = 0; iDim < nDim; iDim++)
        Flow_Dir[iDim] = -UnitNormal[iDim];

      /*--- Dot product of normal and flow direction. This should
       be negative due to outward facing boundary normal convention. ---*/

      alpha = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        alpha += UnitNormal[iDim]*Flow_Dir[iDim];

      /*--- Coefficients in the quadratic equation for the velocity ---*/

      aa =  1.0 + 0.5*Gamma_Minus_One*alpha*alpha;
      bb = -1.0*Gamma_Minus_One*alpha*Riemann;
      cc =  0.5*Gamma_Minus_One*Riemann*Riemann - 2.0*SoundSpeed_Exhaust2/Gamma_Minus_One;

      /*--- Solve quadratic equation for velocity magnitude. Value must
       be positive, so the choice of root is clear. ---*/

      dd      = bb*bb - 4.0*aa*cc;
      dd      = sqrt(max(0.0, dd));
      Vel_Mag = (-bb + dd)/(2.0*aa);

      if (Vel_Mag >= 0.0) {

        Velocity2 = Vel_Mag*Vel_Mag;

        /*--- Compute speed of sound from total speed of sound eqn. ---*/

        SoundSpeed2 = SoundSpeed_Exhaust2 - 0.5*Gamma_Minus_One*Velocity2;
        Mach2       = Velocity2/SoundSpeed2;
        Velocity2   = Mach2*SoundSpeed2;
        Vel_Mag     = sqrt(Velocity2);
        SoundSpeed2 = SoundSpeed_Exhaust2 - 0.5*Gamma_Minus_One*Velocity2;

        /*--- Compute new velocity vector at the inlet ---*/

        for (iDim = 0; iDim < nDim; iDim++)
          Velocity[iDim] = Vel_Mag*Flow_Dir[iDim];

        /*--- Static temperature from the speed of sound relation ---*/

        Temperature = SoundSpeed2/(Gamma*Gas_Constant);

        /*--- Static pressure using isentropic relation at a point ---*/

        Pressure = Exhaust_Pressure*pow((Temperature/Exhaust_Temperature), Gamma/Gamma_Minus_One);

        /*--- Density at the exhaust from the gas law ---*/

        Density = Pressure/(Gas_Constant*Temperature);

        /*--- Using pressure, density, & velocity, compute the energy ---*/

        Energy = Pressure/(Density*Gamma_Minus_One) + 0.5*Velocity2;

        /*--- Primitive variables, using the derived quantities ---*/

        V_exhaust[0] = Temperature;
        for (iDim = 0; iDim < nDim; iDim++)
          V_exhaust[iDim+1] = Velocity[iDim];
        V_exhaust[nDim+1] = Pressure;
        V_exhaust[nDim+2] = Density;
        V_exhaust[nDim+3] = Energy + Pressure/Density;
        V_exhaust[nDim+4] = sqrt(SoundSpeed2);

      }
      /*--- The flow goes in the wrong direction ---*/

      else {

        V_exhaust[0] = V_domain[0];
        for (iDim = 0; iDim < nDim; iDim++)
          V_exhaust[iDim+1] = V_domain[iDim+1];
        V_exhaust[nDim+1] = V_domain[nDim+1];
        V_exhaust[nDim+2] = V_domain[nDim+2];
        V_exhaust[nDim+3] = V_domain[nDim+3];
        V_exhaust[nDim+4] = V_domain[nDim+4];

      }

      /*--- Set various quantities in the solver class ---*/

      conv_numerics->SetNormal(Normal);
      conv_numerics->SetPrimitive(V_domain, V_exhaust);

      /*--- Set grid movement ---*/


      /*--- Compute the residual using an upwind scheme ---*/

      auto residual = conv_numerics->ComputeResidual(config);

      LinSysRes.AddBlock(iPoint, residual);

      /*--- Jacobian contribution for implicit integration ---*/

      if (implicit)
        Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);

    }
  }
  END_SU2_OMP_FOR

  delete [] Normal;

}

void CEulerSolver::SetFreeStream_Solution(const CConfig *config) {

  unsigned long iPoint;
  unsigned short iDim;

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    nodes->SetSolution(iPoint,0, Density_Inf);
    for (iDim = 0; iDim < nDim; iDim++) {
      nodes->SetSolution(iPoint,iDim+1, Density_Inf*Velocity_Inf[iDim]);
    }
    nodes->SetSolution(iPoint,nVar-1, Density_Inf*Energy_Inf);
  }
  END_SU2_OMP_FOR
}
