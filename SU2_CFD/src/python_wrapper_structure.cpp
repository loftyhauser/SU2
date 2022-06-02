/*!
 * \file python_wrapper_structure.cpp
 * \brief Driver subroutines that are used by the Python wrapper. Those routines are usually called from an external Python environment.
 * \author D. Thomas
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


#include "../include/drivers/CDriver.hpp"
#include "../include/drivers/CSinglezoneDriver.hpp"
#include "../../Common/include/toolboxes/geometry_toolbox.hpp"

void CDriver::PythonInterface_Preprocessing(CConfig **config, CGeometry ****geometry, CSolver *****solver){

  int rank = MASTER_NODE;
  SU2_MPI::Comm_rank(SU2_MPI::GetComm(), &rank);

  /* --- Initialize boundary conditions customization, this is achieve through the Python wrapper --- */
  for(iZone=0; iZone < nZone; iZone++){

    if (config[iZone]->GetnMarker_PyCustom() > 0){

      if (rank == MASTER_NODE) cout << endl << "----------------- Python Interface Preprocessing ( Zone "<< iZone <<" ) -----------------" << endl;

      if (rank == MASTER_NODE) cout << "Setting customized boundary conditions for zone " << iZone << endl;
      for (iMesh = 0; iMesh <= config[iZone]->GetnMGLevels(); iMesh++) {
        geometry[iZone][INST_0][iMesh]->SetCustomBoundary(config[iZone]);
      }
      geometry[iZone][INST_0][MESH_0]->UpdateCustomBoundaryConditions(geometry[iZone][INST_0], config[iZone]);

      if (config[iZone]->GetKind_Solver() == MAIN_SOLVER::EULER) {

        solver[iZone][INST_0][MESH_0][FLOW_SOL]->UpdateCustomBoundaryConditions(geometry[iZone][INST_0], config[iZone]);
      }
    }
  }

}

/////////////////////////////////////////////////////////////////////////////
/* Functions related to the global performance indices (Lift, Drag, ecc..) */
/////////////////////////////////////////////////////////////////////////////

passivedouble CDriver::Get_Drag() const {

  unsigned short val_iZone = ZONE_0;
  unsigned short FinestMesh = config_container[val_iZone]->GetFinestMesh();
  su2double CDrag, factor, val_Drag;

  /*--- Calculate drag force based on drag coefficient ---*/
  factor = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetAeroCoeffsReferenceForce();
  CDrag = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CD();

  val_Drag = CDrag*factor;

  return SU2_TYPE::GetValue(val_Drag);
}

passivedouble CDriver::Get_Lift() const {

  unsigned short val_iZone = ZONE_0;
  unsigned short FinestMesh = config_container[val_iZone]->GetFinestMesh();
  su2double CLift, factor, val_Lift;

  /*--- Calculate drag force based on drag coefficient ---*/
  factor = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetAeroCoeffsReferenceForce();
  CLift = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CL();

  val_Lift = CLift*factor;

  return SU2_TYPE::GetValue(val_Lift);
}

passivedouble CDriver::Get_Mx() const {

  unsigned short val_iZone = ZONE_0;
  unsigned short FinestMesh = config_container[val_iZone]->GetFinestMesh();
  su2double CMx, RefLengthCoeff, factor, val_Mx;

  RefLengthCoeff = config_container[val_iZone]->GetRefLength();

  /*--- Calculate moment around x-axis based on coefficients ---*/
  factor = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetAeroCoeffsReferenceForce();
  CMx = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CMx();

  val_Mx = CMx*factor*RefLengthCoeff;

  return SU2_TYPE::GetValue(val_Mx);
}

passivedouble CDriver::Get_My() const {

  unsigned short val_iZone = ZONE_0;
  unsigned short FinestMesh = config_container[val_iZone]->GetFinestMesh();
  su2double CMy, RefLengthCoeff, factor, val_My;

  RefLengthCoeff = config_container[val_iZone]->GetRefLength();

  /*--- Calculate moment around x-axis based on coefficients ---*/
  factor = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetAeroCoeffsReferenceForce();
  CMy = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CMy();

  val_My = CMy*factor*RefLengthCoeff;

  return SU2_TYPE::GetValue(val_My);
}

passivedouble CDriver::Get_Mz() const {

  unsigned short val_iZone = ZONE_0;
  unsigned short FinestMesh = config_container[val_iZone]->GetFinestMesh();
  su2double CMz, RefLengthCoeff, factor, val_Mz;

  RefLengthCoeff = config_container[val_iZone]->GetRefLength();

  /*--- Calculate moment around z-axis based on coefficients ---*/
  factor = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetAeroCoeffsReferenceForce();
  CMz = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CMz();

  val_Mz = CMz*factor*RefLengthCoeff;

  return SU2_TYPE::GetValue(val_Mz);
}

passivedouble CDriver::Get_DragCoeff() const {

  unsigned short val_iZone = ZONE_0;
  unsigned short FinestMesh = config_container[val_iZone]->GetFinestMesh();
  su2double CDrag;

  CDrag = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CD();

  return SU2_TYPE::GetValue(CDrag);
}

passivedouble CDriver::Get_LiftCoeff() const {

  unsigned short val_iZone = ZONE_0;
  unsigned short FinestMesh = config_container[val_iZone]->GetFinestMesh();
  su2double CLift;

  CLift = solver_container[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CL();

  return SU2_TYPE::GetValue(CLift);
}

/////////////////////////////////////////////////////////////////////////////
/* Functions to obtain information from the geometry/mesh                  */
/////////////////////////////////////////////////////////////////////////////

unsigned long CDriver::GetNumberVertices(unsigned short iMarker) const {

  return geometry_container[ZONE_0][INST_0][MESH_0]->nVertex[iMarker];

}

unsigned long CDriver::GetNumberHaloVertices(unsigned short iMarker) const {

  unsigned long nHaloVertices, iVertex, iPoint;

  nHaloVertices = 0;
  for(iVertex = 0; iVertex < geometry_container[ZONE_0][INST_0][MESH_0]->nVertex[iMarker]; iVertex++){
    iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
    if(!(geometry_container[ZONE_0][INST_0][MESH_0]->nodes->GetDomain(iPoint))) nHaloVertices += 1;
  }

  return nHaloVertices;

}

unsigned long CDriver::GetVertexGlobalIndex(unsigned short iMarker, unsigned long iVertex) const {

  unsigned long iPoint, GlobalIndex;

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
  GlobalIndex = geometry_container[ZONE_0][INST_0][MESH_0]->nodes->GetGlobalIndex(iPoint);

  return GlobalIndex;

}

bool CDriver::IsAHaloNode(unsigned short iMarker, unsigned long iVertex) const {

  unsigned long iPoint;

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
  if(geometry_container[ZONE_0][INST_0][MESH_0]->nodes->GetDomain(iPoint)) return false;
  else return true;

}

vector<passivedouble> CDriver::GetVertexNormal(unsigned short iMarker, unsigned long iVertex, bool unitNormal) const {

  su2double *Normal;
  su2double Area;
  vector<su2double> ret_Normal(3, 0.0);
  vector<passivedouble> ret_Normal_passive(3, 0.0);

  Normal = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNormal();

  if (!unitNormal) {

    ret_Normal_passive[0] = SU2_TYPE::GetValue(Normal[0]);
    ret_Normal_passive[1] = SU2_TYPE::GetValue(Normal[1]);
    if(nDim>2) ret_Normal_passive[2] = SU2_TYPE::GetValue(Normal[2]);

    return ret_Normal_passive;
  }

  Area = GeometryToolbox::Norm(nDim, Normal);

  ret_Normal[0] = Normal[0]/Area;
  ret_Normal[1] = Normal[1]/Area;
  if(nDim>2) ret_Normal[2] = Normal[2]/Area;

  ret_Normal_passive[0] = SU2_TYPE::GetValue(ret_Normal[0]);
  ret_Normal_passive[1] = SU2_TYPE::GetValue(ret_Normal[1]);
  ret_Normal_passive[2] = SU2_TYPE::GetValue(ret_Normal[2]);

  return ret_Normal_passive;
}

//////////////////////////////////////////////////////////////////////////////////
/* Functions to obtain global parameters from SU2 (time steps, delta t, ecc...) */
//////////////////////////////////////////////////////////////////////////////////

unsigned long CDriver::GetnTimeIter() const {

  return config_container[ZONE_0]->GetnTime_Iter();
}

unsigned long CDriver::GetTime_Iter() const{

  return TimeIter;
}

passivedouble CDriver::GetUnsteady_TimeStep() const {

  return SU2_TYPE::GetValue(config_container[ZONE_0]->GetTime_Step());
}

string CDriver::GetSurfaceFileName() const {

  return config_container[ZONE_0]->GetSurfCoeff_FileName();
}

///////////////////////////////////////////////////////////////////////////////
/* Functions related to CHT solver                                           */
///////////////////////////////////////////////////////////////////////////////

passivedouble CDriver::GetVertexTemperature(unsigned short iMarker, unsigned long iVertex) const {

  unsigned long iPoint;
  su2double vertexWallTemp(0.0);

  bool compressible = (config_container[ZONE_0]->GetKind_Regime() == ENUM_REGIME::COMPRESSIBLE);

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();

  if(geometry_container[ZONE_0][INST_0][MESH_0]->nodes->GetDomain(iPoint) && compressible){
    vertexWallTemp = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetNodes()->GetTemperature(iPoint);
  }

  return SU2_TYPE::GetValue(vertexWallTemp);

}

void CDriver::SetVertexTemperature(unsigned short iMarker, unsigned long iVertex, passivedouble val_WallTemp){

  geometry_container[ZONE_0][INST_0][MESH_0]->SetCustomBoundaryTemperature(iMarker, iVertex, val_WallTemp);
}

vector<passivedouble> CDriver::GetVertexHeatFluxes(unsigned short iMarker, unsigned long iVertex) const {

  unsigned long iPoint;
  unsigned short iDim;
  su2double Prandtl_Lam  = config_container[ZONE_0]->GetPrandtl_Lam();
  su2double Gas_Constant = config_container[ZONE_0]->GetGas_ConstantND();
  su2double Gamma = config_container[ZONE_0]->GetGamma();
  su2double Gamma_Minus_One = Gamma - 1.0;
  su2double Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
  su2double laminar_viscosity, thermal_conductivity;
  vector<su2double> GradT (3,0.0);
  vector<su2double> HeatFlux (3,0.0);
  vector<passivedouble> HeatFluxPassive (3,0.0);

  bool compressible = (config_container[ZONE_0]->GetKind_Regime() == ENUM_REGIME::COMPRESSIBLE);

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();

  if(compressible){
    laminar_viscosity    = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetNodes()->GetLaminarViscosity(iPoint);
    thermal_conductivity = Cp * (laminar_viscosity/Prandtl_Lam);
    for(iDim=0; iDim < nDim; iDim++){
      GradT[iDim] = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetNodes()->GetGradient_Primitive(iPoint, 0, iDim);
      HeatFlux[iDim] = -thermal_conductivity*GradT[iDim];
    }
  }

  HeatFluxPassive[0] = SU2_TYPE::GetValue(HeatFlux[0]);
  HeatFluxPassive[1] = SU2_TYPE::GetValue(HeatFlux[1]);
  HeatFluxPassive[2] = SU2_TYPE::GetValue(HeatFlux[2]);

  return HeatFluxPassive;
}

passivedouble CDriver::GetVertexNormalHeatFlux(unsigned short iMarker, unsigned long iVertex) const{

  unsigned long iPoint;
  unsigned short iDim;
  su2double vertexWallHeatFlux;
  su2double Prandtl_Lam  = config_container[ZONE_0]->GetPrandtl_Lam();
  su2double Gas_Constant = config_container[ZONE_0]->GetGas_ConstantND();
  su2double Gamma = config_container[ZONE_0]->GetGamma();
  su2double Gamma_Minus_One = Gamma - 1.0;
  su2double Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
  su2double Area;
  su2double laminar_viscosity, thermal_conductivity, dTdn;
  su2double *Normal, GradT[3] = {0.0,0.0,0.0}, UnitNormal[3] = {0.0,0.0,0.0};

  bool compressible = (config_container[ZONE_0]->GetKind_Regime() == ENUM_REGIME::COMPRESSIBLE);

  vertexWallHeatFlux = 0.0;
  dTdn = 0.0;

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();

  if(geometry_container[ZONE_0][INST_0][MESH_0]->nodes->GetDomain(iPoint) && compressible){
    Normal = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNormal();

    Area = GeometryToolbox::Norm(nDim, Normal);

    for (iDim = 0; iDim < nDim; iDim++)
      UnitNormal[iDim] = Normal[iDim]/Area;

    laminar_viscosity    = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetNodes()->GetLaminarViscosity(iPoint);
    thermal_conductivity = Cp * (laminar_viscosity/Prandtl_Lam);
    /*Compute wall heat flux (normal to the wall) based on computed temperature gradient*/
    for(iDim=0; iDim < nDim; iDim++){
      GradT[iDim] = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetNodes()->GetGradient_Primitive(iPoint, 0, iDim);
      dTdn += GradT[iDim]*UnitNormal[iDim];
    }

    vertexWallHeatFlux = -thermal_conductivity*dTdn;
  }

  return SU2_TYPE::GetValue(vertexWallHeatFlux);
}

void CDriver::SetVertexNormalHeatFlux(unsigned short iMarker, unsigned long iVertex, passivedouble val_WallHeatFlux){

  geometry_container[ZONE_0][INST_0][MESH_0]->SetCustomBoundaryHeatFlux(iMarker, iVertex, val_WallHeatFlux);
}

passivedouble CDriver::GetThermalConductivity(unsigned short iMarker, unsigned long iVertex) const {

  unsigned long iPoint;
  su2double Prandtl_Lam  = config_container[ZONE_0]->GetPrandtl_Lam();
  su2double Gas_Constant = config_container[ZONE_0]->GetGas_ConstantND();
  su2double Gamma = config_container[ZONE_0]->GetGamma();
  su2double Gamma_Minus_One = Gamma - 1.0;
  su2double Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
  su2double laminar_viscosity, thermal_conductivity;

  iPoint = geometry_container[ZONE_0][INST_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
  laminar_viscosity    = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetNodes()->GetLaminarViscosity(iPoint);
  thermal_conductivity = Cp * (laminar_viscosity/Prandtl_Lam);

  return SU2_TYPE::GetValue(thermal_conductivity);

}

////////////////////////////////////////////////////////////////////////////////
/* Functions related to the management of markers                             */
////////////////////////////////////////////////////////////////////////////////

vector<string> CDriver::GetAllBoundaryMarkersTag() const {

  vector<string> boundariesTagList;
  unsigned short iMarker,nBoundariesMarkers;
  string Marker_Tag;

  nBoundariesMarkers = config_container[ZONE_0]->GetnMarker_All();
  boundariesTagList.resize(nBoundariesMarkers);

  for(iMarker=0; iMarker < nBoundariesMarkers; iMarker++){
    Marker_Tag = config_container[ZONE_0]->GetMarker_All_TagBound(iMarker);
    boundariesTagList[iMarker] = Marker_Tag;
  }

  return boundariesTagList;
}

vector<string> CDriver::GetAllDeformMeshMarkersTag() const {

  vector<string> interfaceBoundariesTagList;
  unsigned short iMarker, nBoundariesMarker;
  string Marker_Tag;

  nBoundariesMarker = config_container[ZONE_0]->GetnMarker_Deform_Mesh();
  interfaceBoundariesTagList.resize(nBoundariesMarker);

  for(iMarker=0; iMarker < nBoundariesMarker; iMarker++){
    Marker_Tag = config_container[ZONE_0]->GetMarker_Deform_Mesh_TagBound(iMarker);
    interfaceBoundariesTagList[iMarker] = Marker_Tag;
  }

  return interfaceBoundariesTagList;
}

vector<string> CDriver::GetAllCHTMarkersTag() const {

  vector<string> CHTBoundariesTagList;
  unsigned short iMarker, nBoundariesMarker;
  string Marker_Tag;

  nBoundariesMarker = config_container[ZONE_0]->GetnMarker_All();

  //The CHT markers can be identified as the markers that are customizable with a BC type HEAT_FLUX or ISOTHERMAL.
  for(iMarker=0; iMarker<nBoundariesMarker; iMarker++){
    if((config_container[ZONE_0]->GetMarker_All_KindBC(iMarker) == HEAT_FLUX || config_container[ZONE_0]->GetMarker_All_KindBC(iMarker) == ISOTHERMAL) && config_container[ZONE_0]->GetMarker_All_PyCustom(iMarker)){
      Marker_Tag = config_container[ZONE_0]->GetMarker_All_TagBound(iMarker);
      CHTBoundariesTagList.push_back(Marker_Tag);
    }
  }

  return CHTBoundariesTagList;
}

vector<string> CDriver::GetAllInletMarkersTag() const {

  vector<string> BoundariesTagList;
  unsigned short iMarker, nBoundariesMarker;
  string Marker_Tag;

  nBoundariesMarker = config_container[ZONE_0]->GetnMarker_All();

  for(iMarker=0; iMarker<nBoundariesMarker; iMarker++){
    bool isCustomizable = config_container[ZONE_0]->GetMarker_All_PyCustom(iMarker);
    bool isInlet = (config_container[ZONE_0]->GetMarker_All_KindBC(iMarker) == INLET_FLOW);
    if(isCustomizable && isInlet) {
      Marker_Tag = config_container[ZONE_0]->GetMarker_All_TagBound(iMarker);
      BoundariesTagList.push_back(Marker_Tag);
    }
  }

  return BoundariesTagList;
}

map<string, int> CDriver::GetAllBoundaryMarkers() const {

  map<string, int>  allBoundariesMap;
  unsigned short iMarker, nBoundaryMarkers;
  string Marker_Tag;

  nBoundaryMarkers = config_container[ZONE_0]->GetnMarker_All();

  for(iMarker=0; iMarker < nBoundaryMarkers; iMarker++){
    Marker_Tag = config_container[ZONE_0]->GetMarker_All_TagBound(iMarker);
    allBoundariesMap[Marker_Tag] = iMarker;
  }

  return allBoundariesMap;
}

map<string, string> CDriver::GetAllBoundaryMarkersType() const {

  map<string, string> allBoundariesTypeMap;
  unsigned short iMarker, KindBC;
  string Marker_Tag, Marker_Type;

  for(iMarker=0; iMarker < config_container[ZONE_0]->GetnMarker_All(); iMarker++){
    Marker_Tag = config_container[ZONE_0]->GetMarker_All_TagBound(iMarker);
    KindBC = config_container[ZONE_0]->GetMarker_All_KindBC(iMarker);
    switch(KindBC){
      case EULER_WALL:
        Marker_Type = "EULER_WALL";
        break;
      case FAR_FIELD:
        Marker_Type = "FARFIELD";
        break;
      case ISOTHERMAL:
        Marker_Type = "ISOTHERMAL";
        break;
      case HEAT_FLUX:
        Marker_Type = "HEATFLUX";
        break;
      case INLET_FLOW:
        Marker_Type = "INLET_FLOW";
        break;
      case OUTLET_FLOW:
        Marker_Type = "OUTLET_FLOW";
        break;
      case SYMMETRY_PLANE:
        Marker_Type = "SYMMETRY";
        break;
      case SEND_RECEIVE:
        Marker_Type = "SEND_RECEIVE";
        break;
      default:
        Marker_Type = "UNKNOWN_TYPE";
    }
    allBoundariesTypeMap[Marker_Tag] = Marker_Type;
  }

  return allBoundariesTypeMap;
}

void CDriver::SetInlet_Angle(unsigned short iMarker, passivedouble alpha){

  su2double alpha_rad = alpha * PI_NUMBER/180.0;

  unsigned long iVertex;

  for (iVertex = 0; iVertex < geometry_container[ZONE_0][INST_0][MESH_0]->nVertex[iMarker]; iVertex++){
    solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->SetInlet_FlowDir(iMarker, iVertex, 0, cos(alpha_rad));
    solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->SetInlet_FlowDir(iMarker, iVertex, 1, sin(alpha_rad));
  }

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Functions related to simulation control, high level functions (reset convergence, set initial mesh, ecc...) */
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CDriver::ResetConvergence() {

  for(iZone = 0; iZone < nZone; iZone++) {
    switch (config_container[iZone]->GetKind_Solver()) {

    case MAIN_SOLVER::EULER:
      integration_container[iZone][INST_0][FLOW_SOL]->SetConvergence(false);
      break;

    default:
      break;
    }
  }

}

void CDriver::BoundaryConditionsUpdate(){

  int rank = MASTER_NODE;
  unsigned short iZone;

  SU2_MPI::Comm_rank(SU2_MPI::GetComm(), &rank);

  if(rank == MASTER_NODE) cout << "Updating boundary conditions." << endl;
  for(iZone = 0; iZone < nZone; iZone++){
    geometry_container[iZone][INST_0][MESH_0]->UpdateCustomBoundaryConditions(geometry_container[iZone][INST_0], config_container[iZone]);
  }
}

////////////////////////////////////////////////////////////////////////////////
/* Functions related to adjoint simulations                                   */
////////////////////////////////////////////////////////////////////////////////

void CDriver::SetFlowLoad_Adjoint(unsigned short iMarker, unsigned long iVertex, passivedouble val_AdjointX,
                                  passivedouble val_AdjointY, passivedouble val_AdjointZ) {

  CSolver *solver = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL];
  CGeometry *geometry = geometry_container[ZONE_0][INST_0][MESH_0];

  solver->StoreVertexTractionsAdjoint(iMarker, iVertex, 0, val_AdjointX);
  solver->StoreVertexTractionsAdjoint(iMarker, iVertex, 1, val_AdjointY);
  if (geometry->GetnDim() == 3)
    solver->StoreVertexTractionsAdjoint(iMarker, iVertex, 2, val_AdjointZ);

}

////////////////////////////////////////////////////////////////////////////////
/* Functions related to flow loads */
////////////////////////////////////////////////////////////////////////////////

vector<passivedouble> CDriver::GetFlowLoad(unsigned short iMarker, unsigned long iVertex) const {

  vector<su2double> FlowLoad(3, 0.0);
  vector<passivedouble> FlowLoad_passive(3, 0.0);

  CSolver *solver = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL];
  CGeometry *geometry = geometry_container[ZONE_0][INST_0][MESH_0];

  if (config_container[ZONE_0]->GetSolid_Wall(iMarker)) {
    FlowLoad[0] = solver->GetVertexTractions(iMarker, iVertex, 0);
    FlowLoad[1] = solver->GetVertexTractions(iMarker, iVertex, 1);
    if (geometry->GetnDim() == 3)
      FlowLoad[2] = solver->GetVertexTractions(iMarker, iVertex, 2);
    else
      FlowLoad[2] = 0.0;
  }

  FlowLoad_passive[0] = SU2_TYPE::GetValue(FlowLoad[0]);
  FlowLoad_passive[1] = SU2_TYPE::GetValue(FlowLoad[1]);
  FlowLoad_passive[2] = SU2_TYPE::GetValue(FlowLoad[2]);

  return FlowLoad_passive;

}
