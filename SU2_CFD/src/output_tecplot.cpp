/*!
 * \file output_tecplot.cpp
 * \brief Main subroutines for output solver information.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.7
 *
 * Stanford University Unstructured (SU2).
 * Copyright (C) 2012-2013 Aerospace Design Laboratory (ADL).
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

#include "../include/output_structure.hpp"

string AssembleVariableNames(CGeometry *geometry, CConfig *config, unsigned short nVar_Consv, unsigned short *NVar);

void COutput::SetTecplot_ASCII(CConfig *config, CGeometry *geometry, unsigned short val_iZone, unsigned short val_nZone, bool surf_sol) {
    
	/*--- Local variables and initialization ---*/
	unsigned short iDim, iVar, nDim = geometry->GetnDim();
	unsigned short Kind_Solver = config->GetKind_Solver();
    
	unsigned long iPoint, iElem, iNode;
	unsigned long iExtIter = config->GetExtIter();
  unsigned long *LocalIndex = NULL;
  bool *SurfacePoint = NULL;
  
	bool grid_movement  = config->GetGrid_Movement();
	bool adjoint = config->GetAdjoint();
    
	char cstr[200], buffer[50];
	string filename;
  
	/*--- Write file name with extension ---*/
  if (surf_sol) {
    if (adjoint)
      filename = config->GetSurfAdjCoeff_FileName();
    else
      filename = config->GetSurfFlowCoeff_FileName();
  }
  else {
    if (adjoint)
      filename = config->GetAdj_FileName();
    else
      filename = config->GetFlow_FileName();
  }
  
	if (Kind_Solver == WAVE_EQUATION)
		filename = config->GetWave_FileName().c_str();
  
	if ((Kind_Solver == WAVE_EQUATION) && (Kind_Solver == ADJ_AEROACOUSTIC_EULER))
		filename = config->GetAdjWave_FileName().c_str();
  
	if (Kind_Solver == ELECTRIC_POTENTIAL)
		filename = config->GetStructure_FileName().c_str();
  
	if (Kind_Solver == PLASMA_EULER) {
		if (val_iZone == 0) Kind_Solver = PLASMA_EULER;
		if (val_iZone == 1) Kind_Solver = ELECTRIC_POTENTIAL;
	}
	if (Kind_Solver == PLASMA_NAVIER_STOKES) {
		if (val_iZone == 0) Kind_Solver = PLASMA_NAVIER_STOKES;
		if (val_iZone == 1) Kind_Solver = ELECTRIC_POTENTIAL;
	}
    
	strcpy (cstr, filename.c_str());
	if (Kind_Solver == ELECTRIC_POTENTIAL) strcpy (cstr, config->GetStructure_FileName().c_str());
    
	/*--- Special cases where a number needs to be appended to the file name. ---*/
	if ((Kind_Solver == EULER || Kind_Solver == NAVIER_STOKES || Kind_Solver == RANS) &&
        (val_nZone > 1) && (config->GetUnsteady_Simulation() != TIME_SPECTRAL)) {
		sprintf (buffer, "_%d", int(val_iZone));
		strcat(cstr,buffer);
	}
    
	/*--- Special cases where a number needs to be appended to the file name. ---*/
	if (((Kind_Solver == ADJ_EULER) || (Kind_Solver == ADJ_NAVIER_STOKES) || (Kind_Solver == ADJ_RANS)) &&
        (val_nZone > 1) && (config->GetUnsteady_Simulation() != TIME_SPECTRAL)) {
		sprintf (buffer, "_%d", int(val_iZone));
		strcat(cstr,buffer);
	}
    
    /*--- Special cases where a number needs to be appended to the file name. ---*/
	if (((Kind_Solver == ELECTRIC_POTENTIAL) &&( (Kind_Solver == PLASMA_EULER) || (Kind_Solver == PLASMA_NAVIER_STOKES)) )
        && config->GetUnsteady_Simulation()) {
		sprintf (buffer, "_%d", int(iExtIter));
		strcat(cstr,buffer);
	}
    
	if (config->GetUnsteady_Simulation() == TIME_SPECTRAL) {

		/*--- SU2_SOL requires different names. It is only called for parallel cases. ---*/
		if (config->GetKind_SU2() == SU2_SOL) {
			val_iZone = iExtIter;
		}
		if (int(val_iZone) < 10) sprintf (buffer, "_0000%d.dat", int(val_iZone));
		if ((int(val_iZone) >= 10) && (int(val_iZone) < 100)) sprintf (buffer, "_000%d.dat", int(val_iZone));
		if ((int(val_iZone) >= 100) && (int(val_iZone) < 1000)) sprintf (buffer, "_00%d.dat", int(val_iZone));
		if ((int(val_iZone) >= 1000) && (int(val_iZone) < 10000)) sprintf (buffer, "_0%d.dat", int(val_iZone));
		if (int(val_iZone) >= 10000) sprintf (buffer, "_%d.dat", int(val_iZone));
        
	} else if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
		if (int(iExtIter) < 10) sprintf (buffer, "_0000%d.dat", int(iExtIter));
		if ((int(iExtIter) >= 10) && (int(iExtIter) < 100)) sprintf (buffer, "_000%d.dat", int(iExtIter));
		if ((int(iExtIter) >= 100) && (int(iExtIter) < 1000)) sprintf (buffer, "_00%d.dat", int(iExtIter));
		if ((int(iExtIter) >= 1000) && (int(iExtIter) < 10000)) sprintf (buffer, "_0%d.dat", int(iExtIter));
		if (int(iExtIter) >= 10000) sprintf (buffer, "_%d.dat", int(iExtIter));
	} else {
		sprintf (buffer, ".dat");
	}
    
	strcat(cstr,buffer);
    
	/*--- Open Tecplot ASCII file and write the header. ---*/
	ofstream Tecplot_File;
	Tecplot_File.open(cstr, ios::out);
  Tecplot_File.precision(6);
  if (surf_sol) Tecplot_File << "TITLE = \"Visualization of the surface solution\"" << endl;
  else Tecplot_File << "TITLE = \"Visualization of the volumetric solution\"" << endl;

	/*--- Prepare the variable lists. ---*/
  
  /*--- Write the list of the fields in the restart file.
   Without including the PointID---*/
  if (config->GetKind_SU2() == SU2_SOL) {
    
    /*--- If SU2_SOL called this routine, we already have a set of output
     variables with the appropriate string tags stored in the config class. ---*/
    Tecplot_File << "VARIABLES = ";
    nVar_Total = config->fields.size() - 1;
    for (unsigned short iField = 1; iField < config->fields.size(); iField++) {
      Tecplot_File << config->fields[iField];
    }
    Tecplot_File << endl;
    
  } else {
    
    if (nDim == 2) {
      Tecplot_File << "VARIABLES = \"x\",\"y\"";
    } else {
      Tecplot_File << "VARIABLES = \"x\",\"y\",\"z\"";
    }
    
    /*--- Add names for conservative and residual variables ---*/
    for (iVar = 0; iVar < nVar_Consv; iVar++) {
      Tecplot_File << ",\"Conservative_" << iVar+1 << "\"";
    }
    if (config->GetWrt_Residuals()) {
      for (iVar = 0; iVar < nVar_Consv; iVar++) {
        Tecplot_File << ",\"Residual_" << iVar+1 << "\"";
      }
    }
    
    /*--- Add names for any extra variables (this will need to be adjusted). ---*/
    if (grid_movement) {
      if (nDim == 2) {
        Tecplot_File << ",\"Grid_Velx\",\"Grid_Vely\"";
      } else {
        Tecplot_File << ",\"Grid_Velx\",\"Grid_Vely\",\"Grid_Velz\"";
      }
    }
    
    if (config->GetKind_Regime() == FREESURFACE) {
      Tecplot_File << ",\"Density\"";
    }
    
    if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
      Tecplot_File << ",\"Pressure\",\"Pressure_Coefficient\",\"Mach\"";
    }
    
    if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
      Tecplot_File << ", \"Temperature\", \"Laminar_Viscosity\", \"Skin_Friction_Coefficient\", \"Heat_Transfer\", \"Y_Plus\"";
    }
    
    if (Kind_Solver == RANS) {
      Tecplot_File << ", \"Eddy_Viscosity\"";
    }
    
    if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
      Tecplot_File << ", \"Sharp_Edge_Dist\"";
    }
    
    if ((Kind_Solver == PLASMA_EULER) || (Kind_Solver == PLASMA_NAVIER_STOKES)) {
      unsigned short iSpecies;
      for (iSpecies = 0; iSpecies < config->GetnSpecies(); iSpecies++)
        Tecplot_File << ",\"Pressure_" << iSpecies << "\"";
      for (iSpecies = 0; iSpecies < config->GetnSpecies(); iSpecies++)
        Tecplot_File << ",\"Temperature_" << iSpecies << "\"";
      for (iSpecies = 0; iSpecies < config->GetnDiatomics(); iSpecies++)
        Tecplot_File << ",\"TemperatureVib_" << iSpecies << "\"";
      for (iSpecies = 0; iSpecies < config->GetnSpecies(); iSpecies++)
        Tecplot_File << ",\"Mach_" << iSpecies << "\"";
    }
    
    if (Kind_Solver == PLASMA_NAVIER_STOKES) {
      unsigned short iSpecies;
      for (iSpecies = 0; iSpecies < config->GetnSpecies(); iSpecies++)
        Tecplot_File << ",\"LaminaryViscosity_" << iSpecies << "\"";
      
      if ( Kind_Solver == PLASMA_NAVIER_STOKES  && (config->GetMagnetic_Force() == YES) && (geometry->GetnDim() == 3)) {
        for (iDim = 0; iDim < nDim; iDim++)
          Tecplot_File << ",\"Magnet_Field" << iDim << "\"";
      }
    }
    
    if (Kind_Solver == ELECTRIC_POTENTIAL) {
      unsigned short iDim;
      for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
        Tecplot_File << ",\"ElectricField_" << iDim+1 << "\"";
    }
    
    if ((Kind_Solver == ADJ_EULER) || (Kind_Solver == ADJ_NAVIER_STOKES) || (Kind_Solver == ADJ_RANS) || (Kind_Solver == ADJ_PLASMA_EULER) || (Kind_Solver == ADJ_PLASMA_NAVIER_STOKES)) {
      Tecplot_File << ", \"Surface_Sensitivity\", \"Solution_Sensor\"";
    }
    
    if (config->GetExtraOutput()) {
      for (iVar = 0; iVar < nVar_Extra; iVar++) {
        Tecplot_File << ", \"ExtraOutput_" << iVar+1<<"\"";
      }
    }
    
    Tecplot_File << endl;
    
  }
  
  /*--- If it's a surface output, print only the points 
   that are in the element list, change the numbering ---*/
  
  if (surf_sol) {
        
    LocalIndex = new unsigned long [nGlobal_Poin+1];
    SurfacePoint = new bool [nGlobal_Poin+1];

    for (iPoint = 0; iPoint < nGlobal_Poin+1; iPoint++) SurfacePoint[iPoint] = false;

    for(iElem = 0; iElem < nGlobal_Line; iElem++) {
      iNode = iElem*N_POINTS_LINE;
      SurfacePoint[Conn_Line[iNode+0]] = true;
      SurfacePoint[Conn_Line[iNode+1]] = true;
    }
    for(iElem = 0; iElem < nGlobal_BoundTria; iElem++) {
      iNode = iElem*N_POINTS_TRIANGLE;
      SurfacePoint[Conn_BoundTria[iNode+0]] = true;
      SurfacePoint[Conn_BoundTria[iNode+1]] = true;
      SurfacePoint[Conn_BoundTria[iNode+2]] = true;
    }
    for(iElem = 0; iElem < nGlobal_BoundQuad; iElem++) {
      iNode = iElem*N_POINTS_QUADRILATERAL;
      SurfacePoint[Conn_BoundQuad[iNode+0]] = true;
      SurfacePoint[Conn_BoundQuad[iNode+1]] = true;
      SurfacePoint[Conn_BoundQuad[iNode+2]] = true;
      SurfacePoint[Conn_BoundQuad[iNode+3]] = true;
    }
    
    nSurf_Poin = 0;
    for (iPoint = 0; iPoint < nGlobal_Poin+1; iPoint++) {
      LocalIndex[iPoint] = 0;
      if (SurfacePoint[iPoint]) { nSurf_Poin++; LocalIndex[iPoint] = nSurf_Poin; }
    }
    
  }
  
  /*--- Write the header ---*/
  Tecplot_File << "ZONE ";
	if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
		Tecplot_File << "STRANDID="<<int(iExtIter+1)<<", SOLUTIONTIME="<<config->GetDelta_UnstTime()*iExtIter<<", ";
	} else if (config->GetUnsteady_Simulation() == TIME_SPECTRAL) {
		/*--- Compute period of oscillation & compute time interval using nTimeInstances ---*/
		double period = config->GetTimeSpectral_Period();
		double deltaT = period/(double)(config->GetnTimeInstances());
		Tecplot_File << "STRANDID="<<int(iExtIter+1)<<", SOLUTIONTIME="<<deltaT*iExtIter<<", ";
	}

	if (nDim == 2) {
		if (surf_sol) Tecplot_File << "NODES= "<< nSurf_Poin <<", ELEMENTS= "<< nSurf_Elem <<", DATAPACKING=POINT, ZONETYPE=FELINESEG"<< endl;
		else Tecplot_File << "NODES= "<< nGlobal_Poin <<", ELEMENTS= "<< nGlobal_Elem <<", DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL"<< endl;
	} else {
		if (surf_sol) Tecplot_File << "NODES= "<< nSurf_Poin<<", ELEMENTS= "<< nSurf_Elem <<", DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL"<< endl;
		else Tecplot_File << "NODES= "<< nGlobal_Poin <<", ELEMENTS= "<< nGlobal_Elem <<", DATAPACKING=POINT, ZONETYPE=FEBRICK"<< endl;
	}
  
	/*--- Write surface and volumetric solution data. ---*/
  
  for (iPoint = 0; iPoint < nGlobal_Poin; iPoint++) {
    
    if (surf_sol) {
      
      if (LocalIndex[iPoint+1] != 0) {
        
        /*--- Write the node coordinates ---*/
        if (config->GetKind_SU2() != SU2_SOL) {
          for(iDim = 0; iDim < nDim; iDim++)
            Tecplot_File << scientific << Coords[iDim][iPoint] << "\t";
        }
        
        /*--- Loop over the vars/residuals and write the values to file ---*/
        for (iVar = 0; iVar < nVar_Total; iVar++)
          Tecplot_File << scientific << Data[iVar][iPoint] << "\t";
        
        Tecplot_File << endl;

      }
      
    } else {
      
      /*--- Write the node coordinates ---*/
      if (config->GetKind_SU2() != SU2_SOL) {
        for(iDim = 0; iDim < nDim; iDim++)
          Tecplot_File << scientific << Coords[iDim][iPoint] << "\t";
      }
      
      /*--- Loop over the vars/residuals and write the values to file ---*/
      for (iVar = 0; iVar < nVar_Total; iVar++)
        Tecplot_File << scientific << Data[iVar][iPoint] << "\t";
      
      Tecplot_File << endl;
      
    }

  }
  

	/*--- Write connectivity data. ---*/
  if (surf_sol) {
    
    iNode = 0;
    for(iElem = 0; iElem < nGlobal_Line; iElem++) {
      iNode = iElem*N_POINTS_LINE;
      Tecplot_File << LocalIndex[Conn_Line[iNode+0]] << "\t";
      Tecplot_File << LocalIndex[Conn_Line[iNode+1]] << "\n";
      
    }
    
    iNode = 0;
    for(iElem = 0; iElem < nGlobal_BoundTria; iElem++) {
      iNode = iElem*N_POINTS_TRIANGLE;
      Tecplot_File << LocalIndex[Conn_BoundTria[iNode+0]] << "\t";
      Tecplot_File << LocalIndex[Conn_BoundTria[iNode+1]] << "\t";
      Tecplot_File << LocalIndex[Conn_BoundTria[iNode+2]] << "\t";
      Tecplot_File << LocalIndex[Conn_BoundTria[iNode+2]] << "\n";
    }
    
    iNode = 0;
    for(iElem = 0; iElem < nGlobal_BoundQuad; iElem++) {
      iNode = iElem*N_POINTS_QUADRILATERAL;
      Tecplot_File << LocalIndex[Conn_BoundQuad[iNode+0]] << "\t";
      Tecplot_File << LocalIndex[Conn_BoundQuad[iNode+1]] << "\t";
      Tecplot_File << LocalIndex[Conn_BoundQuad[iNode+2]] << "\t";
      Tecplot_File << LocalIndex[Conn_BoundQuad[iNode+3]] << "\n";
    }
    
  }
  else {
    
    iNode = 0;
    for(iElem = 0; iElem < nGlobal_Tria; iElem++) {
      iNode = iElem*N_POINTS_TRIANGLE;
      Tecplot_File << Conn_Tria[iNode+0] << "\t";
      Tecplot_File << Conn_Tria[iNode+1] << "\t";
      Tecplot_File << Conn_Tria[iNode+2] << "\t";
      Tecplot_File << Conn_Tria[iNode+2] << "\n";
    }
    
    iNode = 0;
    for(iElem = 0; iElem < nGlobal_Quad; iElem++) {
      iNode = iElem*N_POINTS_QUADRILATERAL;
      Tecplot_File << Conn_Quad[iNode+0] << "\t";
      Tecplot_File << Conn_Quad[iNode+1] << "\t";
      Tecplot_File << Conn_Quad[iNode+2] << "\t";
      Tecplot_File << Conn_Quad[iNode+3] << "\n";
    }
    
    iNode = 0;
    for(iElem = 0; iElem < nGlobal_Tetr; iElem++) {
      iNode = iElem*N_POINTS_TETRAHEDRON;
      Tecplot_File << Conn_Tetr[iNode+0] << "\t" << Conn_Tetr[iNode+1] << "\t";
      Tecplot_File << Conn_Tetr[iNode+2] << "\t" << Conn_Tetr[iNode+2] << "\t";
      Tecplot_File << Conn_Tetr[iNode+3] << "\t" << Conn_Tetr[iNode+3] << "\t";
      Tecplot_File << Conn_Tetr[iNode+3] << "\t" << Conn_Tetr[iNode+3] << "\n";
    }
    
    iNode = 0;
    for(iElem = 0; iElem < nGlobal_Hexa; iElem++) {
      iNode = iElem*N_POINTS_HEXAHEDRON;
      Tecplot_File << Conn_Hexa[iNode+0] << "\t" << Conn_Hexa[iNode+1] << "\t";
      Tecplot_File << Conn_Hexa[iNode+2] << "\t" << Conn_Hexa[iNode+3] << "\t";
      Tecplot_File << Conn_Hexa[iNode+4] << "\t" << Conn_Hexa[iNode+5] << "\t";
      Tecplot_File << Conn_Hexa[iNode+6] << "\t" << Conn_Hexa[iNode+7] << "\n";
    }
    
    iNode = 0;
    for(iElem = 0; iElem < nGlobal_Wedg; iElem++) {
      iNode = iElem*N_POINTS_WEDGE;
      Tecplot_File << Conn_Wedg[iNode+0] << "\t" << Conn_Wedg[iNode+1] << "\t";
      Tecplot_File << Conn_Wedg[iNode+1] << "\t" << Conn_Wedg[iNode+2] << "\t";
      Tecplot_File << Conn_Wedg[iNode+3] << "\t" << Conn_Wedg[iNode+4] << "\t";
      Tecplot_File << Conn_Wedg[iNode+4] << "\t" << Conn_Wedg[iNode+5] << "\n";
    }
    
    iNode = 0;
    for(iElem = 0; iElem < nGlobal_Pyra; iElem++) {
      iNode = iElem*N_POINTS_PYRAMID;
      Tecplot_File << Conn_Pyra[iNode+0] << "\t" << Conn_Pyra[iNode+1] << "\t";
      Tecplot_File << Conn_Pyra[iNode+2] << "\t" << Conn_Pyra[iNode+3] << "\t";
      Tecplot_File << Conn_Pyra[iNode+4] << "\t" << Conn_Pyra[iNode+4] << "\t";
      Tecplot_File << Conn_Pyra[iNode+4] << "\t" << Conn_Pyra[iNode+4] << "\n";
    }
  }
    
	Tecplot_File.close();
  
  if (surf_sol) delete [] LocalIndex;

}

string AssembleVariableNames(CGeometry *geometry, CConfig *config, unsigned short nVar_Consv, unsigned short *NVar) {
  
  /*--- Local variables ---*/
  stringstream variables; variables.str(string());
	unsigned short iVar;
	*NVar = 0;
  unsigned short iDim, nDim = geometry->GetnDim();
  unsigned short Kind_Solver  = config->GetKind_Solver();
  bool grid_movement = config->GetGrid_Movement();
  bool Wrt_Unsteady = config->GetWrt_Unsteady();
  
  
	/*--- Write the basic variable header based on the particular solution ----*/
  
  /*--- Write the list of the fields in the restart file.
   Without including the PointID---*/
  if (config->GetKind_SU2() == SU2_SOL) {
    
    /*--- If SU2_SOL called this routine, we already have a set of output
     variables with the appropriate string tags stored in the config class. ---*/
    
    /*--- Set the number of variables to be written. Subtract off an index for
     the PointID as well as each coordinate (x,y,z). ---*/
    string varname;
    
    if (Wrt_Unsteady && grid_movement) {
      
      *NVar = config->fields.size()-1;
      for (unsigned short iField = 1; iField < config->fields.size(); iField++) {
        varname = config->fields[iField];
        varname.erase (varname.begin(), varname.begin()+1);
        varname.erase (varname.end()-2, varname.end());
        variables << varname << " ";
      }
    } else {
      
      *NVar = config->fields.size()-1-nDim;
      for (unsigned short iField = 1+nDim; iField < config->fields.size(); iField++) {
        varname = config->fields[iField];
        varname.erase (varname.begin(), varname.begin()+1);
        varname.erase (varname.end()-2, varname.end());
        variables << varname << " ";
      }
    }

  } else {
    
    if (Wrt_Unsteady && grid_movement) {
      if (nDim == 2) {
        variables << "x y "; *NVar += 2;
      } else {
        variables << "x y z "; *NVar += 3;
      }
    }
    
    for (iVar = 0; iVar < nVar_Consv; iVar++) {
      variables << "Conservative_" << iVar+1<<" "; *NVar += 1;
    }
    if (config->GetWrt_Residuals()) {
      for (iVar = 0; iVar < nVar_Consv; iVar++) {
        variables << "Residual_" << iVar+1<<" "; *NVar += 1;
      }
    }
    
    /*--- Add names for any extra variables (this will need to be adjusted). ---*/
    if (grid_movement) {
      if (nDim == 2) {
        variables << "Grid_Velx Grid_Vely "; *NVar += 2;
      } else {
        variables << "Grid_Velx Grid_Vely Grid_Velz "; *NVar += 3;
      }
    }
    
    if (config->GetKind_Regime() == FREESURFACE) {
      variables << "Density ";
      *NVar += 1;
    }
    
    if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
      variables << "Pressure Pressure_Coefficient Mach ";
      *NVar += 3;
    }
    
    if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
      variables << "Temperature Laminar_Viscosity Skin_Friction_Coefficient Heat_Transfer Y_Plus ";
      *NVar += 5;
    }
    
    if (Kind_Solver == RANS) {
      variables << "Eddy_Viscosity ";
      *NVar += 1;
    }
    
    if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
      variables << "Sharp_Edge_Dist ";
      *NVar += 1;
    }
    
    if ((Kind_Solver == PLASMA_EULER) || (Kind_Solver == PLASMA_NAVIER_STOKES)) {
      unsigned short iSpecies;
      for (iSpecies = 0; iSpecies < config->GetnSpecies(); iSpecies++) {
        variables << "Pressure_" << iSpecies << " ";
        *NVar += 1;
      }
      for (iSpecies = 0; iSpecies < config->GetnSpecies(); iSpecies++) {
        variables << "Temperature_" << iSpecies << " ";
        *NVar += 1;
      }
      for (iSpecies = 0; iSpecies < config->GetnDiatomics(); iSpecies++) {
        variables << "TemperatureVib_" << iSpecies << " ";
        *NVar += 1;
      }
      for (iSpecies = 0; iSpecies < config->GetnSpecies(); iSpecies++) {
        variables << "Mach_" << iSpecies << " ";
        *NVar += 1;
      }
    }
    
    if (Kind_Solver == PLASMA_NAVIER_STOKES) {
      unsigned short iSpecies;
      for (iSpecies = 0; iSpecies < config->GetnSpecies(); iSpecies++) {
        variables << "LaminaryViscosity_" << iSpecies << " ";
        *NVar += 1;
      }
      
      if ( Kind_Solver == PLASMA_NAVIER_STOKES  && (config->GetMagnetic_Force() == YES) && (geometry->GetnDim() == 3)) {
        for (iDim = 0; iDim < nDim; iDim++) {
          variables << "Magnet_Field" << iDim << " ";
          *NVar += 1;
        }
      }
    }
    
    if (Kind_Solver == ELECTRIC_POTENTIAL) {
      for (iDim = 0; iDim < geometry->GetnDim(); iDim++) {
        variables << "ElectricField_" << iDim+1 << " ";
        *NVar += 1;
      }
    }
    
    if ((Kind_Solver == ADJ_EULER) || (Kind_Solver == ADJ_NAVIER_STOKES) || (Kind_Solver == ADJ_RANS) || (Kind_Solver == ADJ_PLASMA_EULER) || (Kind_Solver == ADJ_PLASMA_NAVIER_STOKES)) {
      variables << "Surface_Sensitivity Solution_Sensor ";
      *NVar += 2;
    }
  }

	return variables.str();
  
}
