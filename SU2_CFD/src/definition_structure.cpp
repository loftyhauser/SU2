/*!
 * \file definition_structure.cpp
 * \brief Main subroutines used by SU2_CFD.
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

#include "../include/definition_structure.hpp"


unsigned short GetnZone(string val_mesh_filename, unsigned short val_format, CConfig *config) {
	string text_line, Marker_Tag;
	ifstream mesh_file;
	short nZone = 1;
	bool isFound = false;
	char cstr[200];
	string::size_type position;
	int rank = MASTER_NODE;
    
	/*--- Search the mesh file for the 'NZONE' keyword. ---*/
	switch (val_format) {
        case SU2:
            
            /*--- Open grid file ---*/
            strcpy (cstr, val_mesh_filename.c_str());
            mesh_file.open(cstr, ios::in);
            if (mesh_file.fail()) {
                cout << "cstr=" << cstr << endl;
                cout << "There is no geometry file (GetnZone))!" << endl;
                cout << "Press any key to exit..." << endl;
                cin.get();
                exit(1);
            }
            
            /*--- Open the SU2 mesh file ---*/
            while (getline (mesh_file,text_line)) {
                
                /*--- Search for the "NZONE" keyword to see if there are multiple Zones ---*/
                position = text_line.find ("NZONE=",0);
                if (position != string::npos) {
                    text_line.erase (0,6); nZone = atoi(text_line.c_str()); isFound = true;
                    if (rank == MASTER_NODE) {
                        //					if (nZone == 1) cout << "SU2 mesh file format with a single zone." << endl;
                        //					else if (nZone >  1) cout << "SU2 mesh file format with " << nZone << " zones." << endl;
                        //					else
                        if (nZone <= 0) {
                            cout << "Error: Number of mesh zones is less than 1 !!!" << endl;
                            cout << "Press any key to exit..." << endl;
                            cin.get();
                            exit(1);
                        }
                    }
                }
            }
            /*--- If the "NZONE" keyword was not found, assume this is an ordinary
             simulation on a single Zone ---*/
            if (!isFound) {
                nZone = 1;
                //			if (rank == MASTER_NODE) cout << "SU2 mesh file format with a single zone." << endl;
            }
            break;
            
        case NETCDF_ASCII:
            
            nZone = 1;
            //		if (rank == MASTER_NODE) cout << "NETCDF mesh file format with a single zone." << endl;
            break;
            
	}
    
	/*--- For time spectral integration, nZones = nTimeInstances. ---*/
	if (config->GetUnsteady_Simulation() == TIME_SPECTRAL) {
		nZone = config->GetnTimeInstances();
	}
    
	return (unsigned short) nZone;
}

unsigned short GetnDim(string val_mesh_filename, unsigned short val_format) {
    
	string text_line, Marker_Tag;
	ifstream mesh_file;
	short nDim = 3;
	bool isFound = false;
	char cstr[200];
	string::size_type position;
    
	switch (val_format) {
        case SU2:
            
            /*--- Open grid file ---*/
            strcpy (cstr, val_mesh_filename.c_str());
            mesh_file.open(cstr, ios::in);
            
            /*--- Read SU2 mesh file ---*/
            while (getline (mesh_file,text_line)) {
                /*--- Search for the "NDIM" keyword to see if there are multiple Zones ---*/
                position = text_line.find ("NDIME=",0);
                if (position != string::npos) {
                    text_line.erase (0,6); nDim = atoi(text_line.c_str()); isFound = true;
                }
            }
            break;
            
        case NETCDF_ASCII:
            nDim = 3;
            break;
	}
	return (unsigned short) nDim;
}



void Geometrical_Preprocessing(CGeometry ***geometry, CConfig **config, unsigned short val_nZone) {
    
    
	unsigned short iMGlevel, iZone;
	unsigned long iPoint;
	int rank = MASTER_NODE;
    
	for (iZone = 0; iZone < val_nZone; iZone++) {
        
		/*--- Compute elements surrounding points, points surrounding points,
         and elements surrounding elements ---*/
		if (rank == MASTER_NODE) cout << "Setting local point and element connectivity." <<endl;
		geometry[iZone][MESH_0]->SetEsuP();
		geometry[iZone][MESH_0]->SetPsuP();
		geometry[iZone][MESH_0]->SetEsuE();
        
		/*--- Check the orientation before computing geometrical quantities ---*/
		if (rank == MASTER_NODE) cout << "Checking the numerical grid orientation." <<endl;
		geometry[iZone][MESH_0]->SetBoundVolume();
		geometry[iZone][MESH_0]->Check_Orientation(config[iZone]);
        
		/*--- Create the edge structure ---*/
		if (rank == MASTER_NODE) cout << "Identifying edges and vertices." <<endl;
		geometry[iZone][MESH_0]->SetEdges();
		geometry[iZone][MESH_0]->SetVertex(config[iZone]);
        
		/*--- Compute center of gravity ---*/
		if (rank == MASTER_NODE) cout << "Computing centers of gravity." << endl;
		geometry[iZone][MESH_0]->SetCG();
        
		/*--- Create the control volume structures ---*/
		if (rank == MASTER_NODE) cout << "Setting the control volume structure." << endl;
		geometry[iZone][MESH_0]->SetControlVolume(config[iZone], ALLOCATE);
		geometry[iZone][MESH_0]->SetBoundControlVolume(config[iZone], ALLOCATE);
        
		/*--- Identify closest normal neighbor ---*/
		if (rank == MASTER_NODE) cout << "Searching for the closest normal neighbors to the surfaces." << endl;
		geometry[iZone][MESH_0]->FindNormal_Neighbor(config[iZone]);
        
        /*--- Compute the surface curvature ---*/
		if (rank == MASTER_NODE) cout << "Compute the surface curvature." << endl;
        geometry[iZone][MESH_0]->ComputeSurf_Curvature(config[iZone]);
        
		if ((config[iZone]->GetMGLevels() != 0) && (rank == MASTER_NODE))
			cout << "Setting the multigrid structure." <<endl;
        
	}
    
	/*--- Loop over all the new grid ---*/
	for (iMGlevel = 1; iMGlevel <= config[ZONE_0]->GetMGLevels(); iMGlevel++) {
        
		/*--- Loop over all zones at each grid level. ---*/
		for (iZone = 0; iZone < val_nZone; iZone++) {
            
			/*--- Create main agglomeration structure (including MPI calls) ---*/
			geometry[iZone][iMGlevel] = new CMultiGridGeometry(geometry, config, iMGlevel, iZone);
            
			/*--- Compute points surrounding points. ---*/
			geometry[iZone][iMGlevel]->SetPsuP(geometry[iZone][iMGlevel-1]);
            
			/*--- Create the edge structure ---*/
			geometry[iZone][iMGlevel]->SetEdges();
			geometry[iZone][iMGlevel]->SetVertex(geometry[iZone][iMGlevel-1], config[iZone]);
            
			/*--- Create the control volume structures ---*/
			geometry[iZone][iMGlevel]->SetControlVolume(config[iZone],geometry[iZone][iMGlevel-1], ALLOCATE);
			geometry[iZone][iMGlevel]->SetBoundControlVolume(config[iZone],geometry[iZone][iMGlevel-1], ALLOCATE);
			geometry[iZone][iMGlevel]->SetCoord(geometry[iZone][iMGlevel-1]);
            
			/*--- Find closest neighbor to a surface point ---*/
			geometry[iZone][iMGlevel]->FindNormal_Neighbor(config[iZone]);
            
		}
	}
    
	/*--- For unsteady simulations, initialize the grid volumes
     and coordinates for previous solutions. Loop over all zones/grids. Is this
     the best place for this? ---*/
	for (iZone = 0; iZone < val_nZone; iZone++) {
		if (config[iZone]->GetUnsteady_Simulation() && config[iZone]->GetGrid_Movement()) {
			for (iMGlevel = 0; iMGlevel <= config[iZone]->GetMGLevels(); iMGlevel++) {
				for (iPoint = 0; iPoint < geometry[iZone][iMGlevel]->GetnPoint(); iPoint++) {
					geometry[iZone][iMGlevel]->node[iPoint]->SetVolume_n();
					geometry[iZone][iMGlevel]->node[iPoint]->SetVolume_nM1();
                    
					geometry[iZone][iMGlevel]->node[iPoint]->SetCoord_n();
					geometry[iZone][iMGlevel]->node[iPoint]->SetCoord_n1();
				}
			}
		}
	}
    
}

void Solver_Preprocessing(CSolver ***solver_container, CGeometry **geometry, CConfig *config, unsigned short iZone) {
    
	unsigned short iMGlevel;
	bool euler, ns, plasma_euler, plasma_ns,
	plasma_monatomic, plasma_diatomic,
	lin_euler, lin_ns, turbulent, lin_turb, electric, wave, fea,
	spalart_allmaras, menter_sst, template_solver, transition;
    
    
	/*--- Initialize some useful booleans ---*/
	euler = false;		ns = false; turbulent = false;	electric = false;	plasma_monatomic = false;
	plasma_diatomic = false; plasma_euler = false; plasma_ns = false; transition = false;
	wave = false; fea = false; spalart_allmaras = false;
	lin_euler = false;	lin_ns = false;			lin_turb = false;	menter_sst = false;
	template_solver = false;
    
	/*--- Assign booleans ---*/
	switch (config->GetKind_Solver()) {
        case TEMPLATE_SOLVER: template_solver = true; break;
        case EULER : euler = true; break;
        case NAVIER_STOKES: ns = true; break;
        case RANS : ns = true; turbulent = true; if (config->GetKind_Trans_Model() == LM) transition = true; break;
        case FLUID_STRUCTURE_EULER: euler = true; fea = true; break;
        case FLUID_STRUCTURE_NAVIER_STOKES: ns = true; fea = true; break;
        case FLUID_STRUCTURE_RANS: ns = true; turbulent = true; fea = true; break;
        case AEROACOUSTIC_NAVIER_STOKES: ns = true; wave = true; break;
        case AEROACOUSTIC_RANS: ns = true; turbulent = true; wave = true; break;
        case ELECTRIC_POTENTIAL: electric = true; break;
        case WAVE_EQUATION: wave = true; break;
        case LIN_EULER: euler = true; lin_euler = true; break;
            
            /*--- Specify by zone for the aeroacoustic problem ---*/
        case AEROACOUSTIC_EULER:
            if (iZone == ZONE_0) {
                euler = true;
            } else if (iZone == ZONE_1) {
                wave = true;
            }
            break;
        case PLASMA_EULER:
            if (iZone == ZONE_0) {
                plasma_euler = true;
            } else if (iZone == ZONE_1) {
                electric = true;
            }
            break;
        case PLASMA_NAVIER_STOKES:
            if (iZone == ZONE_0) {
                plasma_ns = true;
            } else if (iZone == ZONE_1) {
                electric = true;
            }
            break;
	}
	/*--- Assign turbulence model booleans --- */
	if (turbulent)
		switch (config->GetKind_Turb_Model()){
            case SA: spalart_allmaras = true; break;
            case SST: menter_sst = true; break;
            default: cout << "Specified turbulence model unavailable or none selected" << endl; cin.get(); break;
		}
    
	if (plasma_euler || plasma_ns) {
		switch (config->GetKind_GasModel()){
            case AIR7: plasma_diatomic = true; break;
            case O2: plasma_diatomic = true; break;
            case N2: plasma_diatomic = true; break;
            case AIR5: plasma_diatomic = true; break;
            case ARGON: plasma_monatomic = true; break;
            case ARGON_SID: plasma_diatomic = true; break;
            default: cout << "Specified plasma model unavailable or none selected" << endl; cin.get(); break;
		}
	}
    
	/*--- Definition of the Class for the solution: solver_container[DOMAIN][MESH_LEVEL][EQUATION]. Note that euler, ns
     and potential are incompatible, they use the same position in sol container ---*/
	for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
        
		/*--- Allocate solution for a template problem ---*/
		if (template_solver) {
			solver_container[iMGlevel][TEMPLATE_SOL] = new CTemplateSolver(geometry[iMGlevel], config);
		}
        
		/*--- Allocate solution for direct problem, and run the preprocessing and postprocessing ---*/
		if (euler) {
			solver_container[iMGlevel][FLOW_SOL] = new CEulerSolver(geometry[iMGlevel], config, iMGlevel);
		}
		if (ns) {
			solver_container[iMGlevel][FLOW_SOL] = new CNSSolver(geometry[iMGlevel], config, iMGlevel);
		}
        
	}
    
}


void Integration_Preprocessing(CIntegration **integration_container, CGeometry **geometry, CConfig *config, unsigned short iZone) {
    
	bool euler, ns, plasma_euler, plasma_ns, plasma_monatomic, plasma_diatomic, 
	lin_euler, lin_ns, turbulent, lin_turb, electric, wave, fea, spalart_allmaras, menter_sst, template_solver, transition;
    
	/*--- Initialize some useful booleans ---*/
	euler = false;		ns = false; turbulent = false;	electric = false;	plasma_monatomic = false;
	plasma_diatomic = false; plasma_euler = false; plasma_ns = false;
	wave = false; fea = false; spalart_allmaras = false;
	lin_euler = false;	lin_ns = false;			lin_turb = false;	menter_sst = false; transition = false;
	template_solver = false;
    
	/*--- Assign booleans ---*/
	switch (config->GetKind_Solver()) {
        case TEMPLATE_SOLVER: template_solver = true; break;
        case EULER : euler = true; break;
        case NAVIER_STOKES: ns = true; break;
        case RANS : ns = true; turbulent = true; if (config->GetKind_Trans_Model() == LM) transition = true; break;
        case FLUID_STRUCTURE_EULER: euler = true; fea = true; break;
        case FLUID_STRUCTURE_NAVIER_STOKES: ns = true; fea = true; break;
        case FLUID_STRUCTURE_RANS: ns = true; turbulent = true; fea = true; break;
        case AEROACOUSTIC_NAVIER_STOKES: ns = true; wave = true; break;
        case AEROACOUSTIC_RANS: ns = true; turbulent = true; wave = true; break;
        case ELECTRIC_POTENTIAL: electric = true; break;
        case WAVE_EQUATION: wave = true; break;
        case LIN_EULER: euler = true; lin_euler = true; break;
            
            /*--- Specify by zone for the aeroacoustic problem ---*/
        case AEROACOUSTIC_EULER:
            if (iZone == ZONE_0) {
                euler = true;
            } else if (iZone == ZONE_1) {
                wave = true;
            }
            break;
        case PLASMA_EULER:
            if (iZone == ZONE_0) {
                plasma_euler = true;
            } else if (iZone == ZONE_1) {
                electric = true;
            }
            break;
        case PLASMA_NAVIER_STOKES:
            if (iZone == ZONE_0) {
                plasma_ns = true;
            } else if (iZone == ZONE_1) {
                electric = true;
            }
            break;
	}
    
	/*--- Assign turbulence model booleans --- */
	if (turbulent) {
		switch (config->GetKind_Turb_Model()) {
            case SA: spalart_allmaras = true; break;
            case SST: menter_sst = true; break;
            default: cout << "Specified turbulence model unavailable or none selected" << endl; cin.get(); break;
		}
    }
    
	if (plasma_euler || plasma_ns) {
		switch (config->GetKind_GasModel()) {
            case AIR7: plasma_diatomic = true; break;
            case O2: plasma_diatomic = true; break;
            case N2: plasma_diatomic = true; break;
            case AIR5: plasma_diatomic = true; break;
            case ARGON: plasma_monatomic = true; break;
            case AIR21: plasma_diatomic = true; break;
            case ARGON_SID: plasma_diatomic = true; break;
            default: cout << "Specified plasma model unavailable or none selected" << endl; cin.get(); break;
		}
	}
    
	/*--- Allocate solution for a template problem ---*/
	if (template_solver) integration_container[TEMPLATE_SOL] = new CSingleGridIntegration(config);
    
	/*--- Allocate solution for direct problem ---*/
	if (euler) integration_container[FLOW_SOL] = new CMultiGridIntegration(config);
	if (ns) integration_container[FLOW_SOL] = new CMultiGridIntegration(config);
	if (turbulent) integration_container[TURB_SOL] = new CSingleGridIntegration(config);
	if (transition) integration_container[TRANS_SOL] = new CSingleGridIntegration(config);
	if (electric) integration_container[ELEC_SOL] = new CPotentialIntegration(config);
	if (plasma_euler) integration_container[PLASMA_SOL] = new CMultiGridIntegration(config);
	if (plasma_ns) integration_container[PLASMA_SOL] = new CMultiGridIntegration(config);
	if (wave) integration_container[WAVE_SOL] = new CSingleGridIntegration(config);
	if (fea) integration_container[FEA_SOL] = new CSingleGridIntegration(config);
    
	/*--- Allocate solution for linear problem (at the moment we use the same scheme as the adjoint problem) ---*/
	if (lin_euler) integration_container[LINFLOW_SOL] = new CMultiGridIntegration(config);
	if (lin_ns) { cout <<"Equation not implemented." << endl; cin.get(); }
    
}


void Numerics_Preprocessing(CNumerics ****numerics_container, CSolver ***solver_container, CGeometry **geometry,
                            CConfig *config, unsigned short iZone) {
    
	unsigned short iMGlevel, iSol, nDim, nVar_Template = 0, nVar_Flow = 0, nVar_Trans = 0, nVar_Adj_Flow = 0, nVar_Plasma = 0,
    nVar_Turb = 0, nVar_Adj_Turb = 0, nVar_Elec = 0, nVar_FEA = 0, nVar_Wave = 0, nVar_Lin_Flow = 0,
    nVar_Adj_Plasma = 0, nSpecies = 0, nDiatomics = 0, nMonatomics = 0;
    
	double *constants = NULL;
    
	bool euler, ns, plasma_euler, plasma_ns, plasma_monatomic, plasma_diatomic,
	adj_plasma_euler, adj_plasma_ns, adj_euler, lin_euler, adj_ns, lin_ns, turbulent,
	adj_turb, lin_turb, electric, wave, fea, spalart_allmaras, menter_sst, template_solver, transition;
    
    bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
	bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
	bool freesurface = (config->GetKind_Regime() == FREESURFACE);
    
	/*--- Initialize some useful booleans ---*/
	euler          = false;   ns               = false;   turbulent        = false;
	electric       = false;   plasma_monatomic = false;   plasma_diatomic  = false;  plasma_euler     = false;
	plasma_ns      = false;   adj_euler        = false;	  adj_ns           = false;	 adj_turb         = false;
	wave           = false;   fea              = false;   spalart_allmaras = false;
	lin_euler      = false;   lin_ns           = false;   lin_turb         = false;	 menter_sst       = false;   adj_plasma_euler = false;
	adj_plasma_ns  = false;   transition       = false;   template_solver  = false;
    
	/*--- Assign booleans ---*/
	switch (config->GetKind_Solver()) {
        case TEMPLATE_SOLVER: template_solver = true; break;
        case EULER : euler = true; break;
        case NAVIER_STOKES: ns = true; break;
        case RANS : ns = true; turbulent = true; if (config->GetKind_Trans_Model() == LM) transition = true; break;
        case FLUID_STRUCTURE_EULER: euler = true; fea = true; break;
        case FLUID_STRUCTURE_NAVIER_STOKES: ns = true; fea = true; break;
        case FLUID_STRUCTURE_RANS: ns = true; turbulent = true; fea = true; break;
        case AEROACOUSTIC_NAVIER_STOKES: ns = true; wave = true; break;
        case AEROACOUSTIC_RANS: ns = true; turbulent = true; wave = true; break;
        case ELECTRIC_POTENTIAL: electric = true; break;
        case WAVE_EQUATION: wave = true; break;
        case ADJ_EULER : euler = true; adj_euler = true; break;
        case ADJ_NAVIER_STOKES : ns = true; turbulent = (config->GetKind_Turb_Model() != NONE); adj_ns = true; break;
        case ADJ_RANS : ns = true; turbulent = true; adj_ns = true; adj_turb = (!config->GetFrozen_Visc()); break;
        case ADJ_PLASMA_EULER : plasma_euler = true; adj_plasma_euler = true; break;
        case ADJ_PLASMA_NAVIER_STOKES : plasma_ns = true; adj_plasma_ns = true; break;
        case LIN_EULER: euler = true; lin_euler = true; break;
            
            /*--- Specify by zone for the aeroacoustic problem ---*/
        case AEROACOUSTIC_EULER:
            if (iZone == ZONE_0) {
                euler = true;
            } else if (iZone == ZONE_1) {
                wave = true;
            }
            break;
        case ADJ_AEROACOUSTIC_EULER:
            if (iZone == ZONE_0) {
                euler = true; adj_euler = true;
            } else if (iZone == ZONE_1) {
                wave = true;
            }
            break;
        case PLASMA_EULER:
            if (iZone == ZONE_0) {
                plasma_euler = true;
            } else if (iZone == ZONE_1) {
                electric = true;
            }
            break;
        case PLASMA_NAVIER_STOKES:
            if (iZone == ZONE_0) {
                plasma_ns = true;
            } else if (iZone == ZONE_1) {
                electric = true;
            }
            break;
	}
    
	/*--- Assign turbulence model booleans --- */
	if (turbulent)
		switch (config->GetKind_Turb_Model()){
            case SA: spalart_allmaras = true; break;
            case SST: menter_sst = true; constants = solver_container[MESH_0][TURB_SOL]->GetConstants(); break;
            default: cout << "Specified turbulence model unavailable or none selected" << endl; cin.get(); break;
		}
    
	/*--- Number of variables for the template ---*/
	if (template_solver) nVar_Flow = solver_container[MESH_0][FLOW_SOL]->GetnVar();
    
	/*--- Number of variables for direct problem ---*/
	if (euler)				nVar_Flow = solver_container[MESH_0][FLOW_SOL]->GetnVar();
	if (ns)	nVar_Flow = solver_container[MESH_0][FLOW_SOL]->GetnVar();
	if (turbulent)		nVar_Turb = solver_container[MESH_0][TURB_SOL]->GetnVar();
	if (transition)		nVar_Trans = solver_container[MESH_0][TRANS_SOL]->GetnVar();
    
	/*--- Number of variables for the linear problem ---*/
	if (lin_euler)	nVar_Lin_Flow = solver_container[MESH_0][LINFLOW_SOL]->GetnVar();
    
	/*--- Number of dimensions ---*/
	nDim = geometry[MESH_0]->GetnDim();
    
	/*--- Definition of the Class for the numerical method: numerics_container[MESH_LEVEL][EQUATION][EQ_TERM] ---*/
	for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
		numerics_container[iMGlevel] = new CNumerics** [MAX_SOLS];
		for (iSol = 0; iSol < MAX_SOLS; iSol++)
			numerics_container[iMGlevel][iSol] = new CNumerics* [MAX_TERMS];
	}
    
	/*--- Solver definition for the template problem ---*/
	if (template_solver) {
        
		/*--- Definition of the convective scheme for each equation and mesh level ---*/
		switch (config->GetKind_ConvNumScheme_Template()) {
            case SPACE_CENTERED : case SPACE_UPWIND :
                for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
                    numerics_container[iMGlevel][TEMPLATE_SOL][CONV_TERM] = new CConvective_Template(nDim, nVar_Template, config);
                break;
            default : cout << "Convective scheme not implemented (template_solver)." << endl; cin.get(); break;
		}
        
		/*--- Definition of the viscous scheme for each equation and mesh level ---*/
		switch (config->GetKind_ViscNumScheme_Template()) {
            case AVG_GRAD : case AVG_GRAD_CORRECTED : case GALERKIN :
                for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
                    numerics_container[iMGlevel][TEMPLATE_SOL][VISC_TERM] = new CViscous_Template(nDim, nVar_Template, config);
                break;
            default : cout << "Viscous scheme not implemented." << endl; cin.get(); break;
		}
        
		/*--- Definition of the source term integration scheme for each equation and mesh level ---*/
		switch (config->GetKind_SourNumScheme_Template()) {
            case PIECEWISE_CONSTANT :
                for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
                    numerics_container[iMGlevel][TEMPLATE_SOL][SOURCE_FIRST_TERM] = new CSource_Template(nDim, nVar_Template, config);
                break;
            default : cout << "Source term not implemented." << endl; cin.get(); break;
		}
        
		/*--- Definition of the boundary condition method ---*/
		for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
			numerics_container[iMGlevel][TEMPLATE_SOL][CONV_BOUND_TERM] = new CConvective_Template(nDim, nVar_Template, config);
		}
        
	}
    
	/*--- Solver definition for the Potential, Euler, Navier-Stokes problems ---*/
	if ((euler) || (ns)) {
        
		/*--- Definition of the convective scheme for each equation and mesh level ---*/
		switch (config->GetKind_ConvNumScheme_Flow()) {
            case NO_CONVECTIVE :
                cout << "No convective scheme." << endl; cin.get();
                break;
                
            case SPACE_CENTERED :
                if (compressible) {
                    /*--- Compressible flow ---*/
                    switch (config->GetKind_Centered_Flow()) {
                        case NO_CENTERED : cout << "No centered scheme." << endl; break;
                        case LAX : numerics_container[MESH_0][FLOW_SOL][CONV_TERM] = new CCentLax_Flow(nDim,nVar_Flow, config); break;
                        case JST : numerics_container[MESH_0][FLOW_SOL][CONV_TERM] = new CCentJST_Flow(nDim,nVar_Flow, config); break;
                        default : cout << "Centered scheme not implemented." << endl; cin.get(); break;
                    }
                    
                    if (!config->GetLowFidelitySim()) {
                        for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++)
                            numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CCentLax_Flow(nDim, nVar_Flow, config);
                    }
                    else {
                        numerics_container[MESH_1][FLOW_SOL][CONV_TERM] = new CCentJST_Flow(nDim, nVar_Flow, config);
                        for (iMGlevel = 2; iMGlevel <= config->GetMGLevels(); iMGlevel++)
                            numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CCentLax_Flow(nDim, nVar_Flow, config);
                    }
                    
                    /*--- Definition of the boundary condition method ---*/
                    for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
                        numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwRoe_Flow(nDim, nVar_Flow, config);
                    
                }
                if (incompressible) {
                    /*--- Incompressible flow, use artificial compressibility method ---*/
                    switch (config->GetKind_Centered_Flow()) {
                        case NO_CENTERED : cout << "No centered scheme." << endl; break;
                        case LAX : numerics_container[MESH_0][FLOW_SOL][CONV_TERM] = new CCentLaxArtComp_Flow(nDim, nVar_Flow, config); break;
                        case JST : numerics_container[MESH_0][FLOW_SOL][CONV_TERM] = new CCentJSTArtComp_Flow(nDim, nVar_Flow, config); break;
                        default : cout << "Centered scheme not implemented." << endl; cin.get(); break;
                    }
                    for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++)
                        numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CCentLaxArtComp_Flow(nDim,nVar_Flow, config);
                    
                    /*--- Definition of the boundary condition method ---*/
                    for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
                        numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwRoeArtComp_Flow(nDim, nVar_Flow, config);
                    
                }
                if (freesurface) {
                    /*--- FreeSurface flow, use artificial compressibility method ---*/
                    cout << "Centered scheme not implemented." << endl; cin.get();
                }
                break;
            case SPACE_UPWIND :
                if (compressible) {
                    /*--- Compressible flow ---*/
                    switch (config->GetKind_Upwind_Flow()) {
                        case NO_UPWIND : cout << "No upwind scheme." << endl; break;
                        case ROE_1ST : case ROE_2ND :
                            for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
                                numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwRoe_Flow(nDim, nVar_Flow, config);
                                numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwRoe_Flow(nDim, nVar_Flow, config);
                            }
                            break;
                            
                        case AUSM_1ST : case AUSM_2ND :
                            for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
                                numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwAUSM_Flow(nDim, nVar_Flow, config);
                                numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwAUSM_Flow(nDim, nVar_Flow, config);
                            }
                            break;
                            
                        case ROE_TURKEL_1ST : case ROE_TURKEL_2ND :
                            for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
                                numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwRoe_Turkel_Flow(nDim, nVar_Flow, config);
                                numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwRoe_Turkel_Flow(nDim, nVar_Flow, config);
                            }
                            break;
                            
                        case HLLC_1ST : case HLLC_2ND :
                            for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
                                numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwHLLC_Flow(nDim, nVar_Flow, config);
                                numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwHLLC_Flow(nDim, nVar_Flow, config);
                            }
                            break;
                            
                        default : cout << "Upwind scheme not implemented." << endl; cin.get(); break;
                    }
                    
                }
                if (incompressible) {
                    /*--- Incompressible flow, use artificial compressibility method ---*/
                    switch (config->GetKind_Upwind_Flow()) {
                        case NO_UPWIND : cout << "No upwind scheme." << endl; break;
                        case ROE_1ST : case ROE_2ND :
                            for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
                                numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwRoeArtComp_Flow(nDim, nVar_Flow, config);
                                numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwRoeArtComp_Flow(nDim, nVar_Flow, config);
                            }
                            break;
                        default : cout << "Upwind scheme not implemented." << endl; cin.get(); break;
                    }
                }
                if (freesurface) {
                    /*--- Incompressible flow, use artificial compressibility method ---*/
                    switch (config->GetKind_Upwind_Flow()) {
                        case NO_UPWIND : cout << "No upwind scheme." << endl; break;
                        case ROE_1ST : case ROE_2ND :
                            for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
                                numerics_container[iMGlevel][FLOW_SOL][CONV_TERM] = new CUpwRoeArtComp_Flow_FreeSurface(nDim, nVar_Flow, config);
                                numerics_container[iMGlevel][FLOW_SOL][CONV_BOUND_TERM] = new CUpwRoeArtComp_Flow_FreeSurface(nDim, nVar_Flow, config);
                            }
                            break;
                        default : cout << "Upwind scheme not implemented." << endl; cin.get(); break;
                    }
                }
                
                break;
                
            default :
                cout << "Convective scheme not implemented (euler and ns)." << endl; cin.get();
                break;
		}
        
		/*--- Definition of the viscous scheme for each equation and mesh level ---*/
		switch (config->GetKind_ViscNumScheme_Flow()) {
            case NONE :
                break;
            case AVG_GRAD :
                if (compressible) {
                    /*--- Compressible flow ---*/
                    for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
                        numerics_container[iMGlevel][FLOW_SOL][VISC_TERM] = new CAvgGrad_Flow(nDim, nVar_Flow, config);
                        numerics_container[iMGlevel][FLOW_SOL][VISC_BOUND_TERM] = new CAvgGrad_Flow(nDim, nVar_Flow, config);
                    }
                }
                if (incompressible) {
                    /*--- Incompressible flow, use artificial compressibility method ---*/
                    for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
                        numerics_container[iMGlevel][FLOW_SOL][VISC_TERM] = new CAvgGradArtComp_Flow(nDim, nVar_Flow, config);
                        numerics_container[iMGlevel][FLOW_SOL][VISC_BOUND_TERM] = new CAvgGradArtComp_Flow(nDim, nVar_Flow, config);
                    }
                }
                if (freesurface) {
                    /*--- Freesurface flow, use artificial compressibility method ---*/
                    for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
                        numerics_container[iMGlevel][FLOW_SOL][VISC_TERM] = new CAvgGradArtComp_Flow(nDim, nVar_Flow, config);
                        numerics_container[iMGlevel][FLOW_SOL][VISC_BOUND_TERM] = new CAvgGradArtComp_Flow(nDim, nVar_Flow, config);
                    }
                }
                break;
            case AVG_GRAD_CORRECTED :
                if (compressible) {
                    /*--- Compressible flow ---*/
                    numerics_container[MESH_0][FLOW_SOL][VISC_TERM] = new CAvgGradCorrected_Flow(nDim, nVar_Flow, config);
                    for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++)
                        numerics_container[iMGlevel][FLOW_SOL][VISC_TERM] = new CAvgGrad_Flow(nDim, nVar_Flow, config);
                    
                    /*--- Definition of the boundary condition method ---*/
                    for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
                        numerics_container[iMGlevel][FLOW_SOL][VISC_BOUND_TERM] = new CAvgGrad_Flow(nDim, nVar_Flow, config);
                }
                if (incompressible) {
                    /*--- Incompressible flow, use artificial compressibility method ---*/
                    numerics_container[MESH_0][FLOW_SOL][VISC_TERM] = new CAvgGradCorrectedArtComp_Flow(nDim, nVar_Flow, config);
                    for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++)
                        numerics_container[iMGlevel][FLOW_SOL][VISC_TERM] = new CAvgGradArtComp_Flow(nDim, nVar_Flow, config);
                    
                    /*--- Definition of the boundary condition method ---*/
                    for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
                        numerics_container[iMGlevel][FLOW_SOL][VISC_BOUND_TERM] = new CAvgGradArtComp_Flow(nDim, nVar_Flow, config);
                }
                if (freesurface) {
                    /*--- Freesurface flow, use artificial compressibility method ---*/
                    numerics_container[MESH_0][FLOW_SOL][VISC_TERM] = new CAvgGradCorrectedArtComp_Flow(nDim, nVar_Flow, config);
                    for (iMGlevel = 1; iMGlevel <= config->GetMGLevels(); iMGlevel++)
                        numerics_container[iMGlevel][FLOW_SOL][VISC_TERM] = new CAvgGradArtComp_Flow(nDim, nVar_Flow, config);
                    
                    /*--- Definition of the boundary condition method ---*/
                    for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++)
                        numerics_container[iMGlevel][FLOW_SOL][VISC_BOUND_TERM] = new CAvgGradArtComp_Flow(nDim, nVar_Flow, config);
                }
                break;
            case GALERKIN :
                cout << "Galerkin viscous scheme not implemented." << endl; cin.get(); exit(1);
                break;
            default :
                cout << "Numerical viscous scheme not recognized." << endl; cin.get(); exit(1);
                break;
		}
        
		/*--- Definition of the source term integration scheme for each equation and mesh level ---*/
		switch (config->GetKind_SourNumScheme_Flow()) {
            case NONE :
                break;
            case PIECEWISE_CONSTANT :
                
                for (iMGlevel = 0; iMGlevel <= config->GetMGLevels(); iMGlevel++) {
                    
                    if (config->GetRotating_Frame() == YES)
                        numerics_container[iMGlevel][FLOW_SOL][SOURCE_FIRST_TERM] = new CSourceRotatingFrame_Flow(nDim, nVar_Flow, config);
                    else if (config->GetAxisymmetric() == YES)
                        numerics_container[iMGlevel][FLOW_SOL][SOURCE_FIRST_TERM] = new CSourceAxisymmetric_Flow(nDim,nVar_Flow, config);
                    else if (config->GetGravityForce() == YES)
                        numerics_container[iMGlevel][FLOW_SOL][SOURCE_FIRST_TERM] = new CSourceGravity(nDim, nVar_Flow, config);
                    else if (config->GetJouleHeating() == YES)
                        numerics_container[iMGlevel][FLOW_SOL][SOURCE_FIRST_TERM] = new CSource_JouleHeating(nDim, nVar_Flow, config);
                    else if (config->GetWind_Gust() == YES)
                        numerics_container[iMGlevel][FLOW_SOL][SOURCE_FIRST_TERM] = new CSourceWindGust(nDim, nVar_Flow, config);
                    else
                        numerics_container[iMGlevel][FLOW_SOL][SOURCE_FIRST_TERM] = new CSourceNothing(nDim, nVar_Flow, config);
                    
                    numerics_container[iMGlevel][FLOW_SOL][SOURCE_SECOND_TERM] = new CSourceNothing(nDim, nVar_Flow, config);
                }
                
                break;
            default :
                cout << "Source term not implemented." << endl; cin.get();
                break;
		}
        
	}
    
}
