/*!
 * \file iteration_structure.cpp
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

#include "../include/iteration_structure.hpp"

void MeanFlowIteration(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
                       CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container) {
  
	double Physical_dt, Physical_t;
	unsigned short iMesh; // Index for multi-grid level
  unsigned short iZone; // Index for zone of the mesh
	bool time_spectral = (config_container[ZONE_0]->GetUnsteady_Simulation() == TIME_SPECTRAL);
	unsigned short nZone = geometry_container[ZONE_0][MESH_0]->GetnZone();
	if (time_spectral) nZone = config_container[ZONE_0]->GetnTimeInstances();
	bool relative_motion = config_container[ZONE_0]->GetRelative_Motion();
  unsigned long IntIter = 0; config_container[ZONE_0]->SetIntIter(IntIter);
  unsigned long ExtIter = config_container[ZONE_0]->GetExtIter();
  
  /*--- Initial set up for unsteady problems with dynamic meshes. ---*/
	for (iZone = 0; iZone < nZone; iZone++) {
        
        
        /*--- Apply a Wind Gust ---*/  //Will need to think about this if I do multiple zones and combined with Aeroelastic
        
        if (config_container[ZONE_0]->GetWind_Gust()) {
            SetWind_GustField(config_container[iZone],geometry_container[iZone],solver_container[iZone]);
        }
        
	}
  
	/*--- If any relative motion between zones was found, perform a search
   and interpolation for any sliding interfaces before the next timestep. ---*/
	if (relative_motion) {
		SetSliding_Interfaces(geometry_container, solver_container, config_container, nZone);
  }
  
	for (iZone = 0; iZone < nZone; iZone++) {
    
		/*--- Set the value of the internal iteration ---*/
		IntIter = ExtIter;
		if ((config_container[iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
				(config_container[iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) IntIter = 0;
    
		/*--- Set the initial condition ---*/
		solver_container[iZone][MESH_0][FLOW_SOL]->SetInitialCondition(geometry_container[iZone], solver_container[iZone], config_container[iZone], ExtIter);
    
		/*--- Update global parameters ---*/
		if (config_container[iZone]->GetKind_Solver() == EULER)         config_container[iZone]->SetGlobalParam(EULER, RUNTIME_FLOW_SYS, ExtIter);
		if (config_container[iZone]->GetKind_Solver() == NAVIER_STOKES) config_container[iZone]->SetGlobalParam(NAVIER_STOKES, RUNTIME_FLOW_SYS, ExtIter);
		if (config_container[iZone]->GetKind_Solver() == RANS)          config_container[iZone]->SetGlobalParam(RANS, RUNTIME_FLOW_SYS, ExtIter);
    
		/*--- Solve the Euler, Navier-Stokes or Reynolds-averaged Navier-Stokes (RANS) equations (one iteration) ---*/
		integration_container[iZone][FLOW_SOL]->MultiGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                config_container, RUNTIME_FLOW_SYS, IntIter, iZone);
        
		if (config_container[iZone]->GetKind_Solver() == RANS) {
      
      /*--- Solve the turbulence model ---*/
			config_container[iZone]->SetGlobalParam(RANS, RUNTIME_TURB_SYS, ExtIter);
			integration_container[iZone][TURB_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                   config_container, RUNTIME_TURB_SYS, IntIter, iZone);
      
			/*--- Solve transition model ---*/
			if (config_container[iZone]->GetKind_Trans_Model() == LM) {
				config_container[iZone]->SetGlobalParam(RANS, RUNTIME_TRANS_SYS, ExtIter);
				integration_container[iZone][TRANS_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,config_container, RUNTIME_TRANS_SYS, IntIter, iZone);
			}
		}
    
		/*--- Compute & store time-spectral source terms across all zones ---*/
		if (time_spectral)
			SetTimeSpectral(geometry_container, solver_container, config_container, nZone, (iZone+1)%nZone);
    
	}
  
	/*--- Dual time stepping strategy ---*/
	if ((config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
      (config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {
    
		for(IntIter = 1; IntIter < config_container[ZONE_0]->GetUnst_nIntIter(); IntIter++) {
      
      /*--- Write the convergence history (only screen output) ---*/
      output->SetConvergence_History(NULL, geometry_container, solver_container, config_container, integration_container, true, 0, ZONE_0);
      
      /*--- Set the value of the internal iteration ---*/
      config_container[ZONE_0]->SetIntIter(IntIter);
      
			/*--- All zones must be advanced and coupled with each pseudo timestep. ---*/
			for (iZone = 0; iZone < nZone; iZone++) {
        
				/*--- Pseudo-timestepping for the Euler, Navier-Stokes or Reynolds-averaged Navier-Stokes equations ---*/
				if (config_container[iZone]->GetKind_Solver() == EULER) config_container[iZone]->SetGlobalParam(EULER, RUNTIME_FLOW_SYS, ExtIter);
				if (config_container[iZone]->GetKind_Solver() == NAVIER_STOKES) config_container[iZone]->SetGlobalParam(NAVIER_STOKES, RUNTIME_FLOW_SYS, ExtIter);
				if (config_container[iZone]->GetKind_Solver() == RANS) config_container[iZone]->SetGlobalParam(RANS, RUNTIME_FLOW_SYS, ExtIter);
				integration_container[iZone][FLOW_SOL]->MultiGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                    config_container, RUNTIME_FLOW_SYS, IntIter, iZone);
        
				/*--- Pseudo-timestepping the turbulence model ---*/
				if (config_container[iZone]->GetKind_Solver() == RANS) {
					/*--- Turbulent model solution ---*/
					config_container[iZone]->SetGlobalParam(RANS, RUNTIME_TURB_SYS, ExtIter);
					integration_container[iZone][TURB_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                       config_container, RUNTIME_TURB_SYS, IntIter, iZone);
					if (config_container[iZone]->GetKind_Trans_Model() == LM) {
						/*--- Transition model solution ---*/
						config_container[iZone]->SetGlobalParam(RANS, RUNTIME_TRANS_SYS, ExtIter);
						integration_container[iZone][TRANS_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                          config_container, RUNTIME_TRANS_SYS, IntIter, iZone);
					}
				}
        
        
			}
      
      if (integration_container[ZONE_0][FLOW_SOL]->GetConvergence()) break;
      
		}
    
		for (iZone = 0; iZone < nZone; iZone++) {
			/*--- Update dual time solver on all mesh levels ---*/
			for (iMesh = 0; iMesh <= config_container[iZone]->GetMGLevels(); iMesh++) {
				integration_container[iZone][FLOW_SOL]->SetDualTime_Solver(geometry_container[iZone][iMesh], solver_container[iZone][iMesh][FLOW_SOL], config_container[iZone]);
				integration_container[iZone][FLOW_SOL]->SetConvergence(false);
			}
      
			/*--- Update dual time solver for the turbulence model ---*/
			if (config_container[iZone]->GetKind_Solver() == RANS) {
				integration_container[iZone][TURB_SOL]->SetDualTime_Solver(geometry_container[iZone][MESH_0], solver_container[iZone][MESH_0][TURB_SOL], config_container[iZone]);
				integration_container[iZone][TURB_SOL]->SetConvergence(false);
			}
      
			Physical_dt = config_container[iZone]->GetDelta_UnstTime();
			Physical_t  = (ExtIter+1)*Physical_dt;
			if (Physical_t >=  config_container[iZone]->GetTotal_UnstTime())
				integration_container[iZone][FLOW_SOL]->SetConvergence(true);
		}
	}
  
}

void SetWind_GustField(CConfig *config_container, CGeometry **geometry_container, CSolver ***solver_container) {
    // The gust is imposed on the flow field via the grid velocities. This method called the Field Velocity Method is described in the
    // NASA TM–2012-217771 - Development, Verification and Use of Gust Modeling in the NASA Computational Fluid Dynamics Code FUN3D
    // the desired gust is prescribed as the negative of the grid velocity.
    
    // If a source term is included to account for the gust field, the method is described by Jones et al. as the Split Velocity Method in
    // Simulation of Airfoil Gust Responses Using Prescribed Velocities.
    // In this routine the gust derivatives needed for the source term are calculated. The source term itself is implemented in the class CSourceWindGust
    
    int rank = MASTER_NODE;
    
    /*--- Gust Parameters from config ---*/
    unsigned short Gust_Type = config_container->GetGust_Type();
    double x0 = config_container->GetGust_Begin_Loc();        // Location at which the gust begins.
    double L = config_container->GetGust_WaveLength();        // Gust size
    double tbegin = config_container->GetGust_Begin_Time();   // Physical time at which the gust begins.
    double gust_amp = config_container->GetGust_Ampl();       // Gust amplitude
    double n = config_container->GetGust_Periods();           // Number of gust periods
    unsigned short GustDir = config_container->GetGust_Dir(); // Gust direction
    
    /*--- Check to make sure gust lenght is not zero or negative ---*/
    if (L <= 0.0) {
        cout << "ERROR: The gust length needs to be positive" << endl;
        cout << "Press any key to exit..." << endl;
        cin.get();
        exit(1);
    }

    /*--- Variables needed to compute the gust ---*/
    unsigned short iDim;
    unsigned short nDim = geometry_container[MESH_0]->GetnDim();
    unsigned long iPoint;

    double x, x_gust, gust, dgust_dx, dgust_dy, dgust_dt;
    double *GridVel;
    
    double Physical_dt = config_container->GetDelta_UnstTime();
    unsigned long ExtIter = config_container->GetExtIter();
    double Physical_t = (ExtIter+1)*Physical_dt;
    
    double Uinf = solver_container[MESH_0][FLOW_SOL]->GetVelocity_Inf(0); // Assumption gust moves at infinity velocity
    
    /*--- Loop over all multigrid levels ---*/
    unsigned short iMGlevel, nMGlevel = config_container->GetMGLevels();
    
    for (iMGlevel = 0; iMGlevel <= nMGlevel; iMGlevel++) { //<= ?
    
        /*--- Loop over each node in the volume mesh ---*/
        for (iPoint = 0; iPoint < geometry_container[iMGlevel]->GetnPoint(); iPoint++) {

            //initialize the gust and derivatives to zero everywhere
            gust = 0;
            dgust_dx = 0;
            dgust_dy = 0; // This always stays equal to zero. All the gust implemented don't vary in the y-dir.
            dgust_dt = 0;
            
            // Reset the Grid Velocity to zero if there is no grid movement
                for(iDim = 0; iDim < nDim; iDim++)
                    geometry_container[iMGlevel]->node[iPoint]->SetGridVel(iDim, 0.0);
            
            // Begin applying the gust
            if (Physical_t > tbegin) {
                x = geometry_container[iMGlevel]->node[iPoint]->GetCoord()[0]; // x-location of the node
                // Gust coordinate
                x_gust = (x - x0 - Uinf*(Physical_t-tbegin))/L;
                // Check if we are in the region where the gust is active
                if (x_gust > 0 && x_gust < n) {
                    
                    /*--- Calculate the specified gust ---*/
                    switch (Gust_Type) {
                        
                    case TOP_HAT:
                        gust = gust_amp;
                        // Still need to put the gust derivatives. Think about this.
                        break;
                    
                    case SINE:
                        gust = gust_amp*(sin(2*PI_NUMBER*x_gust));
                        
                        // Gust derivatives
                        dgust_dx = gust_amp*2*PI_NUMBER*(cos(2*PI_NUMBER*x_gust))/L;
                        dgust_dy = 0;
                        dgust_dt = gust_amp*2*PI_NUMBER*(cos(2*PI_NUMBER*x_gust))*(-Uinf)/L;
                        break;
                            
                    case ONE_M_COSINE:
                        gust = 0.5*gust_amp*(1-cos(2*PI_NUMBER*x_gust));
                        
                        // Gust derivatives
                        dgust_dx = 0.5*gust_amp*2*PI_NUMBER*(sin(2*PI_NUMBER*x_gust))/L;
                        dgust_dy = 0;
                        dgust_dt = 0.5*gust_amp*2*PI_NUMBER*(sin(2*PI_NUMBER*x_gust))*(-Uinf)/L;
                        break;
                            
                    case NONE: default:
                        
                        /*--- There is no wind gust specified. ---*/
                        if (rank == MASTER_NODE)
                            cout << "No wind gust specified." << endl;

                    }
                }
            }

            /*--- Set the Wind Gust, Wind Gust Derivatives and the Grid Velocities ---*/
            double NewGridVel[2] = {0,0};
            double Gust[2] = {0,0};
            Gust[GustDir] = gust;
            double GustDer[3] = {0,0,0};
            GustDer[0] = dgust_dx;
            GustDer[1] = dgust_dy;
            GustDer[2] = dgust_dt;
            
            solver_container[iMGlevel][FLOW_SOL]->node[iPoint]->SetWindGust(Gust);
            solver_container[iMGlevel][FLOW_SOL]->node[iPoint]->SetWindGustDer(GustDer);
            
            GridVel = geometry_container[iMGlevel]->node[iPoint]->GetGridVel();
            NewGridVel[GustDir] = GridVel[GustDir] - gust;
            
            // Store new grid velocity
            for(iDim = 0; iDim < nDim; iDim++) {
                geometry_container[iMGlevel]->node[iPoint]->SetGridVel(iDim, NewGridVel[iDim]);
            }
        }
    }
}

void SetTimeSpectral(CGeometry ***geometry_container, CSolver ****solver_container,
		CConfig **config_container, unsigned short nZone, unsigned short iZone) {

	int rank = MASTER_NODE;

	/*--- Local variables and initialization ---*/
	unsigned short iVar, kZone, jZone, iMGlevel;
	unsigned short nVar = solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetnVar();
	unsigned long iPoint;
	bool implicit = (config_container[ZONE_0]->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	bool adjoint = (config_container[ZONE_0]->GetAdjoint());
	if (adjoint) {
		implicit = (config_container[ZONE_0]->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
	}

	/*--- Retrieve values from the config file ---*/
	double *U      	= new double[nVar];
	double *U_old  	= new double[nVar];
	double *Psi	   	= new double[nVar];
	double *Psi_old	= new double[nVar];
	double *Source = new double[nVar];
	double deltaU, deltaPsi;
	//	double Omega[3], Center[3], Ampl[3];
	//	double Lref   = config_container[ZONE_0]->GetLength_Ref();
	//	double current_time;
	//	double alpha_dot[3];
	//  double DEG2RAD = PI_NUMBER/180.0;
	//	double GridVel[3];

//	/*--- rotational velocity ---*/
//	Omega[0]  = config_container[ZONE_0]->GetPitching_Omega_X(ZONE_0)/config_container[ZONE_0]->GetOmega_Ref();
//	Omega[1]  = config_container[ZONE_0]->GetPitching_Omega_Y(ZONE_0)/config_container[ZONE_0]->GetOmega_Ref();
//	Omega[2]  = config_container[ZONE_0]->GetPitching_Omega_Z(ZONE_0)/config_container[ZONE_0]->GetOmega_Ref();
//	double Omega_mag = sqrt(pow(Omega[0],2)+pow(Omega[1],2)+pow(Omega[2],2));
//
//	/*--- amplitude of motion ---*/
//	Ampl[0]   = config_container[ZONE_0]->GetPitching_Ampl_X(ZONE_0)*DEG2RAD;
//	Ampl[1]   = config_container[ZONE_0]->GetPitching_Ampl_Y(ZONE_0)*DEG2RAD;
//	Ampl[2]   = config_container[ZONE_0]->GetPitching_Ampl_Z(ZONE_0)*DEG2RAD;

	/*--- Compute period of oscillation & compute time interval using nTimeInstances ---*/
	double period = config_container[ZONE_0]->GetTimeSpectral_Period();
	//	double deltaT = period/(double)(config_container[ZONE_0]->GetnTimeInstances());

//	/*--- For now, hard code the origin for the pitching airfoil. ---*/
//	Center[0] = config_container[ZONE_0]->GetMotion_Origin_X(ZONE_0);
//	Center[1] = config_container[ZONE_0]->GetMotion_Origin_Y(ZONE_0);
//	Center[2] = config_container[ZONE_0]->GetMotion_Origin_Z(ZONE_0);

	//	/*--- Compute the Mesh Velocities ---*/
	//	if (config_container[ZONE_0]->GetExtIter() == 0) {
	//
	////
	////		/*--- Loop over all grid levels ---*/
	////		for (iMGlevel = 0; iMGlevel <= config_container[ZONE_0]->GetMGLevels(); iMGlevel++) {
	////
	////			/*--- Loop over each zone ---*/
	////			for (iZone = 0; iZone < nZone; iZone++) {
	////
	////				/*--- Time at current instance, i.e. at the current zone ---*/
	////				current_time = (static_cast<double>(iZone))*deltaT;
	////
	////				/*--- Angular velocity at the given time instance ---*/
	////				alpha_dot[0] = -Omega[0]*Ampl[0]*cos(Omega[0]*current_time);
	////				alpha_dot[1] = -Omega[1]*Ampl[1]*cos(Omega[1]*current_time);
	////				alpha_dot[2] = -Omega[2]*Ampl[2]*cos(Omega[2]*current_time);
	////
	////				/*--- Loop over each node in the volume mesh ---*/
	////				for (iPoint = 0; iPoint < geometry_container[iZone][iMGlevel]->GetnPoint(); iPoint++) {
	////
	////					/*--- Coordinates of the current point ---*/
	////					Coord = geometry_container[iZone][iMGlevel]->node[iPoint]->GetCoord();
	////
	////					/*--- Calculate non-dim. position from rotation center ---*/
	////					double r[3] = {0,0,0};
	////					for (iDim = 0; iDim < nDim; iDim++)
	////						r[iDim] = (Coord[iDim]-Center[iDim])/Lref;
	////					if (nDim == 2) r[nDim] = 0.0;
	////
	////					/*--- Cross Product of angular velocity and distance from center ---*/
	////					GridVel[0] = alpha_dot[1]*r[2] - alpha_dot[2]*r[1];
	////					GridVel[1] = alpha_dot[2]*r[0] - alpha_dot[0]*r[2];
	////					GridVel[2] = alpha_dot[0]*r[1] - alpha_dot[1]*r[0];
	////
	////					// PRINT STATEMENTS FOR DEBUGGING!
	////					if (rank == MASTER_NODE && iMGlevel == 0 && iZone == 0) {
	////						if (iPoint == 440) {
	////							cout << "Analytical Mesh Velocity = (" <<  GridVel[0] << ", " << GridVel[1] << ", " << GridVel[2] << ")" << endl;
	////						}
	////					}
	////					/*--- Set Grid Velocity for the point in the given zone ---*/
	////					for(iDim = 0; iDim < nDim; iDim++) {
	////
	////						/*--- Store grid velocity for this point ---*/
	////						geometry_container[iZone][iMGlevel]->node[iPoint]->SetGridVel(iDim, GridVel[iDim]);
	////
	////					}
	////				}
	////			}
	////		}
	//
	//		/*--- by fourier interpolation ---*/
	//
	//	}




	/*--- allocate dynamic memory for D ---*/
	double **D = new double*[nZone];
	for (kZone = 0; kZone < nZone; kZone++) {
		D[kZone] = new double[nZone];
	}

	/*--- Build the time-spectral operator matrix ---*/
	for (kZone = 0; kZone < nZone; kZone++) {
		for (jZone = 0; jZone < nZone; jZone++) {

			/*--- For an even number of time instances ---*/
			if (nZone%2 == 0) {
				if (kZone == jZone) {
					D[kZone][jZone] = 0.0;
				} else {
					D[kZone][jZone] = (PI_NUMBER/period)*pow(-1.0,(kZone-jZone))*(1/tan(PI_NUMBER*(kZone-jZone)/nZone));
				}
			} else {
				/*--- For an odd number of time instances ---*/
				if (kZone == jZone) {
					D[kZone][jZone] = 0.0;
				} else {
					D[kZone][jZone] = (PI_NUMBER/period)*pow(-1.0,(kZone-jZone))*(1/sin(PI_NUMBER*(kZone-jZone)/nZone));
				}
			}
		}
	}

	/*--- Compute various source terms for explicit direct, implicit direct, and adjoint problems ---*/
	/*--- Loop over all grid levels ---*/
	for (iMGlevel = 0; iMGlevel <= config_container[ZONE_0]->GetMGLevels(); iMGlevel++) {

		/*--- Loop over each node in the volume mesh ---*/
		for (iPoint = 0; iPoint < geometry_container[ZONE_0][iMGlevel]->GetnPoint(); iPoint++) {

			for (iVar = 0; iVar < nVar; iVar++) {
				Source[iVar] = 0.0;
			}

			/*--- Step across the columns ---*/
			for (jZone = 0; jZone < nZone; jZone++) {

				/*--- Retrieve solution at this node in current zone ---*/
				for (iVar = 0; iVar < nVar; iVar++) {

					if (!adjoint) {

						U[iVar] = solver_container[jZone][iMGlevel][FLOW_SOL]->node[iPoint]->GetSolution(iVar);
						Source[iVar] += U[iVar]*D[iZone][jZone];

						if (implicit) {
							U_old[iVar] = solver_container[jZone][iMGlevel][FLOW_SOL]->node[iPoint]->GetSolution_Old(iVar);
							deltaU = U[iVar] - U_old[iVar];
							Source[iVar] += deltaU*D[iZone][jZone];
						}

					} else {

						Psi[iVar] = solver_container[jZone][iMGlevel][ADJFLOW_SOL]->node[iPoint]->GetSolution(iVar);
						Source[iVar] += Psi[iVar]*D[jZone][iZone];

						if (implicit) {
							Psi_old[iVar] = solver_container[jZone][iMGlevel][ADJFLOW_SOL]->node[iPoint]->GetSolution_Old(iVar);
							deltaPsi = Psi[iVar] - Psi_old[iVar];
							Source[iVar] += deltaPsi*D[jZone][iZone];
						}
					}
				}


				/*--- Store sources for current row ---*/
				for (iVar = 0; iVar < nVar; iVar++) {
					if (!adjoint) {
						solver_container[iZone][iMGlevel][FLOW_SOL]->node[iPoint]->SetTimeSpectral_Source(iVar,Source[iVar]);
					} else {
						solver_container[iZone][iMGlevel][ADJFLOW_SOL]->node[iPoint]->SetTimeSpectral_Source(iVar,Source[iVar]);
					}
				}

			}
		}
	}



	//	/*--- Loop over all grid levels ---*/
	//	for (iMGlevel = 0; iMGlevel <= config_container[ZONE_0]->GetMGLevels(); iMGlevel++) {
	//
	//		/*--- Loop over each node in the volume mesh ---*/
	//		for (iPoint = 0; iPoint < geometry_container[ZONE_0][iMGlevel]->GetnPoint(); iPoint++) {
	//
	//			for (iZone = 0; iZone < nZone; iZone++) {
	//				for (iVar = 0; iVar < nVar; iVar++) Source[iVar] = 0.0;
	//				for (jZone = 0; jZone < nZone; jZone++) {
	//
	//					/*--- Retrieve solution at this node in current zone ---*/
	//					for (iVar = 0; iVar < nVar; iVar++) {
	//						U[iVar] = solver_container[jZone][iMGlevel][FLOW_SOL]->node[iPoint]->GetSolution(iVar);
	//						Source[iVar] += U[iVar]*D[iZone][jZone];
	//					}
	//				}
	//				/*--- Store sources for current iZone ---*/
	//				for (iVar = 0; iVar < nVar; iVar++)
	//					solver_container[iZone][iMGlevel][FLOW_SOL]->node[iPoint]->SetTimeSpectral_Source(iVar,Source[iVar]);
	//			}
	//		}
	//	}

	/*--- Source term for a turbulence model ---*/
	if (config_container[ZONE_0]->GetKind_Solver() == RANS) {

		/*--- Extra variables needed if we have a turbulence model. ---*/
		unsigned short nVar_Turb = solver_container[ZONE_0][MESH_0][TURB_SOL]->GetnVar();
		double *U_Turb      = new double[nVar_Turb];
		double *Source_Turb = new double[nVar_Turb];

		/*--- Loop over only the finest mesh level (turbulence is always solved
     on the original grid only). ---*/
		for (iPoint = 0; iPoint < geometry_container[ZONE_0][MESH_0]->GetnPoint(); iPoint++) {
			//for (iZone = 0; iZone < nZone; iZone++) {
			for (iVar = 0; iVar < nVar_Turb; iVar++) Source_Turb[iVar] = 0.0;
			for (jZone = 0; jZone < nZone; jZone++) {

				/*--- Retrieve solution at this node in current zone ---*/
				for (iVar = 0; iVar < nVar_Turb; iVar++) {
					U_Turb[iVar] = solver_container[jZone][MESH_0][TURB_SOL]->node[iPoint]->GetSolution(iVar);
					Source_Turb[iVar] += U_Turb[iVar]*D[iZone][jZone];
				}
			}
			/*--- Store sources for current iZone ---*/
			for (iVar = 0; iVar < nVar_Turb; iVar++)
				solver_container[iZone][MESH_0][TURB_SOL]->node[iPoint]->SetTimeSpectral_Source(iVar,Source_Turb[iVar]);
			//}
		}

		delete [] U_Turb;
		delete [] Source_Turb;
  }

  /*--- delete dynamic memory for D ---*/
  for (kZone = 0; kZone < nZone; kZone++) {
	  delete [] D[kZone];
  }
  delete [] D;
  delete [] U;
  delete [] U_old;
  delete [] Psi;
  delete [] Psi_old;
  delete [] Source;

  /*--- Write file with force coefficients ---*/
	ofstream TS_Flow_file;
  ofstream mean_TS_Flow_file;

  /*--- MPI Send/Recv buffers ---*/
  double *sbuf_force = NULL,  *rbuf_force = NULL;

  /*--- Other variables ---*/
  unsigned short nVar_Force = 8;
  unsigned long current_iter = config_container[ZONE_0]->GetExtIter();

  /*--- Allocate memory for send buffer ---*/
  sbuf_force = new double[nVar_Force];

  double *averages = new double[nVar_Force];
  for (iVar = 0; iVar < nVar_Force; iVar++)
	  averages[iVar] = 0;

  /*--- Allocate memory for receive buffer ---*/
  if (rank == MASTER_NODE) {
	  rbuf_force = new double[nVar_Force];

	  TS_Flow_file.precision(15);
	  TS_Flow_file.open("TS_force_coefficients.csv", ios::out);
	  TS_Flow_file <<  "\"time_instance\",\"lift_coeff\",\"drag_coeff\",\"moment_coeff_x\",\"moment_coeff_y\",\"moment_coeff_z\"" << endl;

	  mean_TS_Flow_file.precision(15);
	  if (current_iter == 0 && iZone == 1) {
		  mean_TS_Flow_file.open("history_TS_forces.plt", ios::trunc);
		  mean_TS_Flow_file << "TITLE = \"SU2 TIME-SPECTRAL SIMULATION\"" << endl;
		  mean_TS_Flow_file <<  "VARIABLES = \"Iteration\",\"CLift\",\"CDrag\",\"CMx\",\"CMy\",\"CMz\",\"CT\",\"CQ\",\"CMerit\"" << endl;
		  mean_TS_Flow_file << "ZONE T= \"Average Convergence History\"" << endl;
	  }
	  else
		  mean_TS_Flow_file.open("history_TS_forces.plt", ios::out | ios::app);
  }

  for (kZone = 0; kZone < nZone; kZone++) {

	  /*--- Flow solution coefficients (parallel) ---*/
	  sbuf_force[0] = solver_container[kZone][MESH_0][FLOW_SOL]->GetTotal_CLift();
	  sbuf_force[1] = solver_container[kZone][MESH_0][FLOW_SOL]->GetTotal_CDrag();
	  sbuf_force[2] = solver_container[kZone][MESH_0][FLOW_SOL]->GetTotal_CMx();
	  sbuf_force[3] = solver_container[kZone][MESH_0][FLOW_SOL]->GetTotal_CMy();
	  sbuf_force[4] = solver_container[kZone][MESH_0][FLOW_SOL]->GetTotal_CMz();
	  sbuf_force[5] = solver_container[kZone][MESH_0][FLOW_SOL]->GetTotal_CT();
	  sbuf_force[6] = solver_container[kZone][MESH_0][FLOW_SOL]->GetTotal_CQ();
	  sbuf_force[7] = solver_container[kZone][MESH_0][FLOW_SOL]->GetTotal_CMerit();

	  for (iVar = 0; iVar < nVar_Force; iVar++) {
		  rbuf_force[iVar] = sbuf_force[iVar];
	  }

	  if (rank == MASTER_NODE) {
		  TS_Flow_file << kZone << ", ";
		  for (iVar = 0; iVar < nVar_Force; iVar++)
			  TS_Flow_file << rbuf_force[iVar] << ", ";
		  TS_Flow_file << endl;

		  /*--- Increment the total contributions from each zone, dividing by nZone as you go ---*/
		  for (iVar = 0; iVar < nVar_Force; iVar++) {
			  averages[iVar] += (1.0/double(nZone))*rbuf_force[iVar];
		  }
	  }
  }


  if (rank == MASTER_NODE && iZone == ZONE_0) {

	  mean_TS_Flow_file << current_iter << ", ";
	  for (iVar = 0; iVar < nVar_Force; iVar++) {
		  mean_TS_Flow_file << averages[iVar];
		  if (iVar < nVar_Force-1)
			  mean_TS_Flow_file << ", ";
	  }
	  mean_TS_Flow_file << endl;
  }

  if (rank == MASTER_NODE) {
	  TS_Flow_file.close();
	  mean_TS_Flow_file.close();
	  delete [] rbuf_force;
  }

  delete [] sbuf_force;
  delete [] averages;

}

void SetSliding_Interfaces(CGeometry ***geometry_container, CSolver ****solver_container,
		CConfig **config_container, unsigned short nZone) {

	unsigned short nDim = geometry_container[ZONE_0][MESH_0]->GetnDim();
	unsigned short iMarker, iZone, donorZone;
	unsigned long iVertex, iPoint;
	double *Coord_i, N_0, N_1, N_2;
	unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0;
	double *Coord_0 = NULL, *Coord_1= NULL, *Coord_2= NULL;
	double a[4], b[4], c[4], Area;
	double eps = 1e-10;

	cout << endl << "Searching and interpolating across sliding interfaces." << endl;

	/*--- Check all zones for a sliding interface. ---*/
	for (iZone = 0; iZone < nZone; iZone++) {

		/*--- Check all markers for any SEND/RECV boundaries that might be sliding. ---*/
		for (iMarker = 0; iMarker < config_container[iZone]->GetnMarker_All(); iMarker++) {
			if (config_container[iZone]->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) {

				short SendRecv = config_container[iZone]->GetMarker_All_SendRecv(iMarker);
				/*--- We are only interested in receive information  ---*/
				if (SendRecv < 0) {

					/*--- Loop over all points on this boundary. ---*/
					for (iVertex = 0; iVertex < geometry_container[iZone][MESH_0]->GetnVertex(iMarker); iVertex++) {

						/*--- If the donor is not in a different zone, do nothing (periodic). ---*/
						donorZone = geometry_container[iZone][MESH_0]->vertex[iMarker][iVertex]->GetMatching_Zone();

						if (donorZone != iZone) {
							/*--- Get the index and coordinates of the point for which we need to interpolate a solution. ---*/
							iPoint  = geometry_container[iZone][MESH_0]->vertex[iMarker][iVertex]->GetNode();
							Coord_i = geometry_container[iZone][MESH_0]->node[iPoint]->GetCoord();

							/*--- Loop over the donor marker in the donor zone and search
               the neighboring elements to find the element that contains the
               point to be interpolated. This simply MUST be improved before
               going to any large cases: for now we are doing a brute force
               search, and there are much more efficient alternatives. ---*/

							for (iElem = 0; iElem < geometry_container[donorZone][MESH_0]->GetnElem(); iElem++) {

								/*--- Compute the basis functions for triangular elements. ---*/
								if (nDim == 2) {

									Point_0 = geometry_container[donorZone][MESH_0]->elem[iElem]->GetNode(0);
									Point_1 = geometry_container[donorZone][MESH_0]->elem[iElem]->GetNode(1);
									Point_2 = geometry_container[donorZone][MESH_0]->elem[iElem]->GetNode(2);

									Coord_0 = geometry_container[donorZone][MESH_0]->node[Point_0]->GetCoord();
									Coord_1 = geometry_container[donorZone][MESH_0]->node[Point_1]->GetCoord();
									Coord_2 = geometry_container[donorZone][MESH_0]->node[Point_2]->GetCoord();

									for (unsigned short iDim = 0; iDim < nDim; iDim++) {
										a[iDim] = Coord_0[iDim]-Coord_2[iDim];
										b[iDim] = Coord_1[iDim]-Coord_2[iDim];
									}

									/*--- Norm of the normal component of area, area = 1/2*cross(a,b) ---*/
											Area = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);

									a[0] = 0.5 * (Coord_1[0]*Coord_2[1]-Coord_2[0]*Coord_1[1]) / Area;
									a[1] = 0.5 * (Coord_2[0]*Coord_0[1]-Coord_0[0]*Coord_2[1]) / Area;
									a[2] = 0.5 * (Coord_0[0]*Coord_1[1]-Coord_1[0]*Coord_0[1]) / Area;

									b[0] = 0.5 * (Coord_1[1]-Coord_2[1]) / Area;
									b[1] = 0.5 * (Coord_2[1]-Coord_0[1]) / Area;
									b[2] = 0.5 * (Coord_0[1]-Coord_1[1]) / Area;

									c[0] = 0.5 * (Coord_2[0]-Coord_1[0]) / Area;
									c[1] = 0.5 * (Coord_0[0]-Coord_2[0]) / Area;
									c[2] = 0.5 * (Coord_1[0]-Coord_0[0]) / Area;

									/*--- Basis functions ---*/
									N_0 = a[0] + b[0]*Coord_i[0] + c[0]*Coord_i[1];
									N_1 = a[1] + b[1]*Coord_i[0] + c[1]*Coord_i[1];
									N_2 = a[2] + b[2]*Coord_i[0] + c[2]*Coord_i[1];

									/*--- If 0 <= N <= 1 for all three basis functions, then
                   iPoint is found within this element. Store the element index
                   and the value for the basis functions to be used during the
                   send receive along the sliding interfaces. ---*/
									if ((N_0 >= -eps && N_0 <= 1.0+eps) &&
											(N_1 >= -eps && N_1 <= 1.0+eps) &&
											(N_2 >= -eps && N_2 <= 1.0+eps)) {

										/*--- Store this element as the donor for the SEND/RECV ---*/
										geometry_container[iZone][MESH_0]->vertex[iMarker][iVertex]->SetDonorElem(iElem);
										geometry_container[iZone][MESH_0]->vertex[iMarker][iVertex]->SetBasisFunction(0, N_0);
										geometry_container[iZone][MESH_0]->vertex[iMarker][iVertex]->SetBasisFunction(1, N_1);
										geometry_container[iZone][MESH_0]->vertex[iMarker][iVertex]->SetBasisFunction(2, N_2);

									}

								} else {
									cout << "!!! Error: Sliding mesh interfaces not yet supported in 3-D. !!!" << endl;
									cout << "Press any key to exit..." << endl;
									cin.get();
									exit(1);
								}
							}
						}
					}
				}
			}
		}
	}
}

void SetTimeSpectral_Velocities(CGeometry ***geometry_container,
		CConfig **config_container, unsigned short nZone) {

	unsigned short iZone, jDegree, iDim, iMGlevel;
	unsigned short nDim = geometry_container[ZONE_0][MESH_0]->GetnDim();

	double angular_interval = 2.0*PI_NUMBER/(double)(nZone);
	double *Coord;
	unsigned long iPoint;


	/*--- Compute period of oscillation & compute time interval using nTimeInstances ---*/
	double period = config_container[ZONE_0]->GetTimeSpectral_Period();
	double deltaT = period/(double)(config_container[ZONE_0]->GetnTimeInstances());

	/*--- allocate dynamic memory for angular positions (these are the abscissas) ---*/
	double *angular_positions = new double [nZone];
	for (iZone = 0; iZone < nZone; iZone++) {
		angular_positions[iZone] = iZone*angular_interval;
	}

	/*--- find the highest-degree trigonometric polynomial allowed by the Nyquist criterion---*/
	double high_degree = (nZone-1)/2.0;
	int highest_degree = (int)(high_degree);

	/*--- allocate dynamic memory for a given point's coordinates ---*/
	double **coords = new double *[nZone];
	for (iZone = 0; iZone < nZone; iZone++) {
		coords[iZone] = new double [nDim];
	}

	/*--- allocate dynamic memory for vectors of Fourier coefficients ---*/
	double *a_coeffs = new double [highest_degree+1];
	double *b_coeffs = new double [highest_degree+1];

	/*--- allocate dynamic memory for the interpolated positions and velocities ---*/
	double *fitted_coords = new double [nZone];
	double *fitted_velocities = new double [nZone];

	/*--- Loop over all grid levels ---*/
	for (iMGlevel = 0; iMGlevel <= config_container[ZONE_0]->GetMGLevels(); iMGlevel++) {

		/*--- Loop over each node in the volume mesh ---*/
		for (iPoint = 0; iPoint < geometry_container[ZONE_0][iMGlevel]->GetnPoint(); iPoint++) {

			/*--- Populate the 2D coords array with the
			coordinates of a given mesh point across
			the time instances (i.e. the Zones) ---*/
			/*--- Loop over each zone ---*/
			for (iZone = 0; iZone < nZone; iZone++) {
				/*--- get the coordinates of the given point ---*/
				Coord = geometry_container[iZone][iMGlevel]->node[iPoint]->GetCoord();
				for (iDim = 0; iDim < nDim; iDim++) {
					/*--- add them to the appropriate place in the 2D coords array ---*/
					coords[iZone][iDim] = Coord[iDim];
				}
			}

			/*--- Loop over each Dimension ---*/
			for (iDim = 0; iDim < nDim; iDim++) {

				/*--- compute the Fourier coefficients ---*/
				for (jDegree = 0; jDegree < highest_degree+1; jDegree++) {
					a_coeffs[jDegree] = 0;
					b_coeffs[jDegree] = 0;
					for (iZone = 0; iZone < nZone; iZone++) {
						a_coeffs[jDegree] = a_coeffs[jDegree] + (2.0/(double)nZone)*cos(jDegree*angular_positions[iZone])*coords[iZone][iDim];
						b_coeffs[jDegree] = b_coeffs[jDegree] + (2.0/(double)nZone)*sin(jDegree*angular_positions[iZone])*coords[iZone][iDim];
					}
				}

				/*--- find the interpolation of the coordinates and its derivative (the velocities) ---*/
				for (iZone = 0; iZone < nZone; iZone++) {
					fitted_coords[iZone] = a_coeffs[0]/2.0;
					fitted_velocities[iZone] = 0.0;
					for (jDegree = 1; jDegree < highest_degree+1; jDegree++) {
						fitted_coords[iZone] = fitted_coords[iZone] + a_coeffs[jDegree]*cos(jDegree*angular_positions[iZone]) + b_coeffs[jDegree]*sin(jDegree*angular_positions[iZone]);
						fitted_velocities[iZone] = fitted_velocities[iZone] + (angular_interval/deltaT)*jDegree*(b_coeffs[jDegree]*cos(jDegree*angular_positions[iZone]) - a_coeffs[jDegree]*sin(jDegree*angular_positions[iZone]));
					}
				}

				/*--- Store grid velocity for this point, at this given dimension, across the Zones ---*/
				for (iZone = 0; iZone < nZone; iZone++) {
					geometry_container[iZone][iMGlevel]->node[iPoint]->SetGridVel(iDim, fitted_velocities[iZone]);
				}



			}
		}
	}

	/*--- delete dynamic memory for the abscissas, coefficients, et cetera ---*/
	delete [] angular_positions;
	delete [] a_coeffs;
	delete [] b_coeffs;
	delete [] fitted_coords;
	delete [] fitted_velocities;
	for (iZone = 0; iZone < nDim; iZone++) {
		delete [] coords[iZone];
	}
	delete [] coords;

}
