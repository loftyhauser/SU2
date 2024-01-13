/*!
 * \file solution_direct_levelset.cpp
 * \brief Main subrotuines for solving the level set problem.
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

#include "../include/solver_structure.hpp"

CLevelSetSolver::CLevelSetSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CSolver() {
	unsigned short iVar, iDim;
	unsigned long iPoint, index;
	double dull_val, levelset = 0.0, XCoord = 0.0, YCoord = 0.0, ZCoord = 0.0;
	string text_line;
  
	bool restart = (config->GetRestart() || config->GetRestart_Flow());
  bool rans = ((config->GetKind_Solver() == RANS) || (config->GetKind_Solver() == ADJ_RANS));
  
  int rank = MASTER_NODE;
  
	/*--- Define geometry constans in the solver structure ---*/
	nDim = geometry->GetnDim();
  nPoint = geometry->GetnPoint();
	nPointDomain = geometry->GetnPointDomain();
  
	node = new CVariable*[nPoint];
	
	/*--- Dimension of the problem ---*/
	nVar = 1;
	
  /*--- Single grid simulation ---*/
	if (iMesh == MESH_0) {
    
    /*--- Define some auxiliar vector related with the residual ---*/
    Residual      = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]      = 0.0;
    Residual_RMS  = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 0.0;
    Residual_Max  = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 0.0;
    Point_Max  = new unsigned long[nVar]; for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar]  = 0;
    Residual_i    = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_i[iVar]    = 0.0;
    Residual_j    = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_j[iVar]    = 0.0;
    
    /*--- Define some auxiliar vector related with the solution ---*/
    Solution    = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution[iVar]    = 0.0;
    Solution_i  = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution_i[iVar]  = 0.0;
    Solution_j  = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution_j[iVar]  = 0.0;
    
    /*--- Define some auxiliar vector related with the geometry ---*/
    Vector    = new double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector[iDim]    = 0.0;
    Vector_i  = new double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_i[iDim]  = 0.0;
    Vector_j  = new double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_j[iDim]  = 0.0;
    
    /*--- Define some auxiliar vector related with the flow solution ---*/
    FlowSolution_i = new double[nDim+1]; for (iVar = 0; iVar < nDim+1; iVar++) FlowSolution_i[iVar]  = 0.0;
    FlowSolution_j = new double[nDim+1]; for (iVar = 0; iVar < nDim+1; iVar++) FlowSolution_j[iVar]  = 0.0;
    
    LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
    LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
    
    /*--- Jacobians and vector structures for implicit computations ---*/
    if (config->GetKind_TimeIntScheme_LevelSet() == EULER_IMPLICIT) {
      
      /*--- Block auxiliar Jacobians ---*/
      Jacobian_i = new double* [nVar];
      Jacobian_j = new double* [nVar];
      for (iVar = 0; iVar < nVar; iVar++) {
        Jacobian_i[iVar] = new double [nVar];
        Jacobian_j[iVar] = new double [nVar];
      }
      
      Jacobian_MeanFlow_i = new double* [nDim+1];
      Jacobian_MeanFlow_j = new double* [nDim+1];
      for (iVar = 0; iVar < nDim+1; iVar++) {
        Jacobian_MeanFlow_i[iVar] = new double [nDim+1];
        Jacobian_MeanFlow_j[iVar] = new double [nDim+1];
      }
      
      /*--- Initialization of the structure of the whole Jacobian ---*/
      if (rank == MASTER_NODE) cout << "Initialize jacobian structure (Level Set). MG level: 0." << endl;
      Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, geometry);
    }
    
    /*--- Computation of gradients by least squares ---*/
    if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
      
      /*--- S matrix := inv(R)*traspose(inv(R)) ---*/
      Smatrix = new double* [nDim];
      for (iDim = 0; iDim < nDim; iDim++)
        Smatrix[iDim] = new double [nDim];
      /*--- c vector := transpose(WA)*(Wb) ---*/
      cvector = new double* [nVar];
      for (iVar = 0; iVar < nVar; iVar++)
        cvector[iVar] = new double [nDim];
    }
    
  }
	
	/*--- Restart the solution from file information ---*/
	if (!restart || geometry->GetFinestMGLevel() == false) {
    
    /*--- Restart the solution from uniform distance ---*/
		for (iPoint = 0; iPoint < nPoint; iPoint++) {
      XCoord = geometry->node[iPoint]->GetCoord(0);
      YCoord = geometry->node[iPoint]->GetCoord(1);
      if (nDim == 2) levelset = YCoord - config->GetFreeSurface_Zero();
      else {
        ZCoord = geometry->node[iPoint]->GetCoord(2);
        levelset = ZCoord - config->GetFreeSurface_Zero();
      }
			node[iPoint] = new CLevelSetVariable(levelset, nDim, nVar, config);
		}
    
	}
	else {
    
		/*--- Restart the solution from file information ---*/
		ifstream restart_file;
		string filename = config->GetSolution_FlowFileName();
    
		/*--- Append time step for unsteady restart ---*/
		if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
			char buffer[50];
			unsigned long flowIter = config->GetnExtIter() - 1;
      unsigned short lastindex = filename.find_last_of(".");
      filename = filename.substr(0, lastindex);
      
			if ((int(flowIter) >= 0) && (int(flowIter) < 10)) sprintf (buffer, "_0000%d.dat", int(flowIter));
			if ((int(flowIter) >= 10) && (int(flowIter) < 100)) sprintf (buffer, "_000%d.dat", int(flowIter));
			if ((int(flowIter) >= 100) && (int(flowIter) < 1000)) sprintf (buffer, "_00%d.dat", int(flowIter));
			if ((int(flowIter) >= 1000) && (int(flowIter) < 10000)) sprintf (buffer, "_0%d.dat", int(flowIter));
			if (int(flowIter) >= 10000) sprintf (buffer, "_%d.dat", int(flowIter));
			string UnstExt = string(buffer);
			filename.append(UnstExt);
		}
		restart_file.open(filename.data(), ios::in);
    
    /*--- In case there is no restart file ---*/
		if (restart_file.fail()) {
			cout << "There is no level set restart file!!" << endl;
			cout << "Press any key to exit..." << endl;
			cin.get();
			exit(1);
		}
    
    /*--- In case this is a parallel simulation, we need to perform the
     Global2Local index transformation first. ---*/
    long *Global2Local;
    Global2Local = new long[geometry->GetGlobal_nPointDomain()];
    /*--- First, set all indices to a negative value by default ---*/
    for(iPoint = 0; iPoint < geometry->GetGlobal_nPointDomain(); iPoint++) {
      Global2Local[iPoint] = -1;
    }
    /*--- Now fill array with the transform values only for local points ---*/
    for(iPoint = 0; iPoint < nPointDomain; iPoint++) {
      Global2Local[geometry->node[iPoint]->GetGlobalIndex()] = iPoint;
    }
    
    /*--- Read all lines in the restart file ---*/
    long iPoint_Local; unsigned long iPoint_Global = 0;
    
    /*--- The first line is the header ---*/
    getline (restart_file, text_line);
    
    while (getline (restart_file, text_line)) {
			istringstream point_line(text_line);
      
      /*--- Retrieve local index. If this node from the restart file lives
       on a different processor, the value of iPoint_Local will be -1.
       Otherwise, the local index for this node on the current processor
       will be returned and used to instantiate the vars. ---*/
      iPoint_Local = Global2Local[iPoint_Global];
      if (iPoint_Local >= 0) {
        if (!rans) {
          if (nDim == 2) point_line >> index >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> Solution[0];
          if (nDim == 3) point_line >> index >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> Solution[0];
        }
        else {
          if (nDim == 2) point_line >> index >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> Solution[0];
          if (nDim == 3) point_line >> index >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> Solution[0];
        }
        node[iPoint_Local] = new CLevelSetVariable(Solution[0], nDim, nVar, config);
      }
      iPoint_Global++;
    }
    
    /*--- Instantiate the variable class with an arbitrary solution
     at any halo/periodic nodes. The initial solution can be arbitrary,
     because a send/recv is performed immediately in the solver. ---*/
    for(iPoint = nPointDomain; iPoint < nPoint; iPoint++) {
      node[iPoint] = new CLevelSetVariable(Solution[0], nDim, nVar, config);
    }
    
		/*--- Close the restart file ---*/
		restart_file.close();
    
    /*--- Free memory needed for the transformation ---*/
    delete [] Global2Local;
	}
  
  /*--- MPI solution ---*/
  Set_MPI_Solution(geometry, config);
  
}

CLevelSetSolver::~CLevelSetSolver(void) { }

void CLevelSetSolver::Set_MPI_Solution(CGeometry *geometry, CConfig *config) {
	unsigned short iVar, iMarker, MarkerS, MarkerR;
	unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector, nBufferS_Scalar, nBufferR_Scalar;
	double *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;
	int send_to, receive_from;
  
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
		if ((config->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
			
			MarkerS = iMarker;  MarkerR = iMarker+1;
      
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
			receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
			
			nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
			nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;
      nBufferS_Scalar = nVertexS;             nBufferR_Scalar = nVertexR;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_U = new double [nBufferR_Vector];
      Buffer_Send_U = new double[nBufferS_Vector];
      
      /*--- Copy the solution that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution(iVar);
      }
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Receive_U[iVar*nVertexR+iVertex] = Buffer_Send_U[iVar*nVertexR+iVertex];
      }
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_U;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        
        /*--- Copy conservative variables. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetSolution(iVar, Buffer_Receive_U[iVar*nVertexR+iVertex]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_U;
      
    }
    
	}
  
}

void CLevelSetSolver::Set_MPI_Solution_Old(CGeometry *geometry, CConfig *config) {
	unsigned short iVar, iMarker, MarkerS, MarkerR;
	unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
	double *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;
	int send_to, receive_from;
  
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
		if ((config->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
			
			MarkerS = iMarker;  MarkerR = iMarker+1;
      
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
			receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
			
			nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
			nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_U = new double [nBufferR_Vector];
      Buffer_Send_U = new double[nBufferS_Vector];
      
      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution_Old(iVar);
      }
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Receive_U[iVar*nVertexR+iVertex] = Buffer_Send_U[iVar*nVertexR+iVertex];
      }
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_U;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        
        /*--- Copy transformed conserved variables back into buffer. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetSolution_Old(iVar, Buffer_Receive_U[iVar*nVertexR+iVertex]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_U;
      
    }
    
	}
}

void CLevelSetSolver::Set_MPI_Solution_Limiter(CGeometry *geometry, CConfig *config) {
	unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
	unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
	double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_Limit = NULL, *Buffer_Send_Limit = NULL;
	int send_to, receive_from;
  
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
		if ((config->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
			
			MarkerS = iMarker;  MarkerR = iMarker+1;
      
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
			receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
			
			nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
			nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Limit = new double [nBufferR_Vector];
      Buffer_Send_Limit = new double[nBufferS_Vector];
      
      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_Limit[iVar*nVertexS+iVertex] = node[iPoint]->GetLimiter(iVar);
      }
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Receive_Limit[iVar*nVertexR+iVertex] = Buffer_Send_Limit[iVar*nVertexR+iVertex];
      }
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_Limit;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        iPeriodic_Index = geometry->vertex[MarkerR][iVertex]->GetRotation_Type();
        
        /*--- Retrieve the supplied periodic information. ---*/
        angles = config->GetPeriodicRotation(iPeriodic_Index);
        
        /*--- Store angles separately for clarity. ---*/
        theta    = angles[0];   phi    = angles[1];     psi    = angles[2];
        cosTheta = cos(theta);  cosPhi = cos(phi);      cosPsi = cos(psi);
        sinTheta = sin(theta);  sinPhi = sin(phi);      sinPsi = sin(psi);
        
        /*--- Compute the rotation matrix. Note that the implicit
         ordering is rotation about the x-axis, y-axis,
         then z-axis. Note that this is the transpose of the matrix
         used during the preprocessing stage. ---*/
        rotMatrix[0][0] = cosPhi*cosPsi;    rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;     rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
        rotMatrix[0][1] = cosPhi*sinPsi;    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;     rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
        rotMatrix[0][2] = -sinPhi;          rotMatrix[1][2] = sinTheta*cosPhi;                              rotMatrix[2][2] = cosTheta*cosPhi;
        
        /*--- Copy conserved variables before performing transformation. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          Solution[iVar] = Buffer_Receive_Limit[iVar*nVertexR+iVertex];
        
        /*--- Rotate the momentum components. ---*/
        if (nDim == 2) {
          Solution[1] = rotMatrix[0][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_Limit[2*nVertexR+iVertex];
          Solution[2] = rotMatrix[1][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_Limit[2*nVertexR+iVertex];
        }
        else {
          Solution[1] = rotMatrix[0][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_Limit[2*nVertexR+iVertex] +
          rotMatrix[0][2]*Buffer_Receive_Limit[3*nVertexR+iVertex];
          Solution[2] = rotMatrix[1][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_Limit[2*nVertexR+iVertex] +
          rotMatrix[1][2]*Buffer_Receive_Limit[3*nVertexR+iVertex];
          Solution[3] = rotMatrix[2][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[2][1]*Buffer_Receive_Limit[2*nVertexR+iVertex] +
          rotMatrix[2][2]*Buffer_Receive_Limit[3*nVertexR+iVertex];
        }
        
        /*--- Copy transformed conserved variables back into buffer. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetLimiter(iVar, Solution[iVar]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_Limit;
      
    }
    
	}
}

void CLevelSetSolver::Set_MPI_Solution_Gradient(CGeometry *geometry, CConfig *config) {
	unsigned short iVar, iDim, iMarker, iPeriodic_Index, MarkerS, MarkerR;
	unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
	double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_Gradient = NULL, *Buffer_Send_Gradient = NULL;
	int send_to, receive_from;
  
  double **Gradient = new double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Gradient[iVar] = new double[nDim];
  
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
		if ((config->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
			
			MarkerS = iMarker;  MarkerR = iMarker+1;
      
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
			receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
			
			nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
			nBufferS_Vector = nVertexS*nVar*nDim;        nBufferR_Vector = nVertexR*nVar*nDim;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Gradient = new double [nBufferR_Vector];
      Buffer_Send_Gradient = new double[nBufferS_Vector];
      
      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Buffer_Send_Gradient[iDim*nVar*nVertexS+iVar*nVertexS+iVertex] = node[iPoint]->GetGradient(iVar, iDim);
      }
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Buffer_Receive_Gradient[iDim*nVar*nVertexR+iVar*nVertexR+iVertex] = Buffer_Send_Gradient[iDim*nVar*nVertexR+iVar*nVertexR+iVertex];
      }
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_Gradient;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        iPeriodic_Index = geometry->vertex[MarkerR][iVertex]->GetRotation_Type();
        
        /*--- Retrieve the supplied periodic information. ---*/
        angles = config->GetPeriodicRotation(iPeriodic_Index);
        
        /*--- Store angles separately for clarity. ---*/
        theta    = angles[0];   phi    = angles[1];     psi    = angles[2];
        cosTheta = cos(theta);  cosPhi = cos(phi);      cosPsi = cos(psi);
        sinTheta = sin(theta);  sinPhi = sin(phi);      sinPsi = sin(psi);
        
        /*--- Compute the rotation matrix. Note that the implicit
         ordering is rotation about the x-axis, y-axis,
         then z-axis. Note that this is the transpose of the matrix
         used during the preprocessing stage. ---*/
        rotMatrix[0][0] = cosPhi*cosPsi;    rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;     rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
        rotMatrix[0][1] = cosPhi*sinPsi;    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;     rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
        rotMatrix[0][2] = -sinPhi;          rotMatrix[1][2] = sinTheta*cosPhi;                              rotMatrix[2][2] = cosTheta*cosPhi;
        
        /*--- Copy conserved variables before performing transformation. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Gradient[iVar][iDim] = Buffer_Receive_Gradient[iDim*nVar*nVertexR+iVar*nVertexR+iVertex];
        
        /*--- Need to rotate the gradients for all conserved variables. ---*/
        for (iVar = 0; iVar < nVar; iVar++) {
          if (nDim == 2) {
            Gradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex];
          }
          else {
            Gradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][2]*Buffer_Receive_Gradient[2*nVar*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][2]*Buffer_Receive_Gradient[2*nVar*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][2] = rotMatrix[2][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[2][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[2][2]*Buffer_Receive_Gradient[2*nVar*nVertexR+iVar*nVertexR+iVertex];
          }
        }
        
        /*--- Store the received information ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            node[iPoint]->SetGradient(iVar, iDim, Gradient[iVar][iDim]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_Gradient;
      
    }
    
	}
  
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Gradient[iVar];
  delete [] Gradient;
  
}


void CLevelSetSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem) {
  
	unsigned long iPoint;
	
	bool implicit = (config->GetKind_TimeIntScheme_LevelSet() == EULER_IMPLICIT);
	bool upwind_2nd = (config->GetKind_Upwind_LevelSet() == SCALAR_UPWIND_2ND);
  bool limiter = (config->GetKind_SlopeLimit_LevelSet() != NONE);
  
	for (iPoint = 0; iPoint < nPoint; iPoint ++) {
    
    /*--- Initialize the convective residual vector ---*/
		LinSysRes.SetBlock_Zero(iPoint);

	}
  
  /*--- Upwind second order reconstruction ---*/
  if (upwind_2nd) {
		if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);
		if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);
    
    /*--- Limiter computation ---*/
    if (limiter) SetSolution_Limiter(geometry, config);
	}
  
  /*--- Implicit part ---*/
	if (implicit) Jacobian.SetValZero();
  
}

void CLevelSetSolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                     unsigned short iMesh) {
  bool output, reevaluation;
  
  /*--- Compute level set function using the distance to the free surface ---*/
  if (config->GetIntIter() == 0) output = true;
  else output = false;
  
  if (config->GetIntIter() % config->GetFreeSurface_Reevaluation() == 0) reevaluation = true;
  else reevaluation = false;
  
  SetFreeSurface_Distance(geometry, config, reevaluation, output);
    
}

void CLevelSetSolver::Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short iMesh) {
	double *LevelSet_i, *LevelSet_j, *Limiter_i = NULL, *Limiter_j = NULL, **Gradient_i, **Gradient_j, Project_Grad_i, Project_Grad_j;
	unsigned long iEdge, iPoint, jPoint;
	unsigned short iDim;
  
	bool implicit = (config->GetKind_TimeIntScheme_LevelSet() == EULER_IMPLICIT);
	bool high_order_diss = (config->GetKind_Upwind_LevelSet() == SCALAR_UPWIND_2ND);
  bool limiter = (config->GetKind_SlopeLimit_LevelSet() != NONE);
  
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
		/*--- Points in edge and normal vectors ---*/
		iPoint = geometry->edge[iEdge]->GetNode(0); jPoint = geometry->edge[iEdge]->GetNode(1);
		numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
		
		/*--- Primitive variables w/o reconstruction at time n (density, vel)---*/
    numerics->SetPrimitive(solver_container[FLOW_SOL]->node[iPoint]->GetPrimVar(), solver_container[FLOW_SOL]->node[jPoint]->GetPrimVar());
    
		/*--- Level Set variables w/o reconstruction ---*/
		LevelSet_i = node[iPoint]->GetSolution(); LevelSet_j = node[jPoint]->GetSolution();
    
		if (high_order_diss) {
      
			for (iDim = 0; iDim < nDim; iDim++) {
				Vector_i[iDim] = 0.5*(geometry->node[jPoint]->GetCoord(iDim) - geometry->node[iPoint]->GetCoord(iDim));
				Vector_j[iDim] = 0.5*(geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim));
			}
			
			/*--- Level Set variables using gradient reconstruction ---*/
			Gradient_i = node[iPoint]->GetGradient(); Gradient_j = node[jPoint]->GetGradient();
      if (limiter) { Limiter_i = node[iPoint]->GetLimiter(); Limiter_j = node[jPoint]->GetLimiter(); }
      
      Project_Grad_i = 0.0; Project_Grad_j = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Project_Grad_i += Vector_i[iDim]*Gradient_i[0][iDim];
        Project_Grad_j += Vector_j[iDim]*Gradient_j[0][iDim];
      }
      if (limiter) {
        Solution_i[0] = LevelSet_i[0] + Project_Grad_i*Limiter_i[0];
        Solution_j[0] = LevelSet_j[0] + Project_Grad_j*Limiter_j[0];
      }
      else {
        Solution_i[0] = LevelSet_i[0] + Project_Grad_i;
        Solution_j[0] = LevelSet_j[0] + Project_Grad_j;
      }
      
      numerics->SetLevelSetVar(Solution_i, Solution_j);
      
		} else {
      
      /*--- Set conservative variables without reconstruction ---*/
      numerics->SetLevelSetVar(LevelSet_i, LevelSet_j);
            
    }
		
		/*--- Add and subtract Residual ---*/
		numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, Jacobian_MeanFlow_i, Jacobian_MeanFlow_j, config);
		
		LinSysRes.AddBlock(iPoint, Residual);
		LinSysRes.SubtractBlock(jPoint, Residual);
		
		/*--- Implicit part ---*/
		if (implicit) {
			Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
			Jacobian.AddBlock(iPoint, jPoint, Jacobian_j);
			Jacobian.SubtractBlock(jPoint, iPoint, Jacobian_i);
			Jacobian.SubtractBlock(jPoint, jPoint, Jacobian_j);
		}
	}
}

void CLevelSetSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics, CConfig *config, unsigned short iMesh) {
	unsigned long iPoint;
	double Vol, x_o, x_od, x, z, levelset, DampingFactor;
  
	double factor = config->GetFreeSurface_Damping_Length();
	bool implicit = (config->GetKind_TimeIntScheme_LevelSet() == EULER_IMPLICIT);
  
	x_o = config->GetFreeSurface_Outlet();
	x_od = x_o - factor*2.0*PI_NUMBER*config->GetFroude()*config->GetFroude();
  
	for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
		
		Vol = geometry->node[iPoint]->GetVolume();
		x = geometry->node[iPoint]->GetCoord()[0];
		z = geometry->node[iPoint]->GetCoord()[nDim-1]-config->GetFreeSurface_Zero();
		levelset = node[iPoint]->GetSolution(0);
    
		DampingFactor = 0.0;
		if (x >= x_od)
			DampingFactor = config->GetFreeSurface_Damping_Coeff()*pow((x-x_od)/(x_o-x_od), 2.0);
    
		Residual[0] = Vol*(levelset-z)*DampingFactor;
		Jacobian_i[0][0] = Vol*DampingFactor;
		
		LinSysRes.AddBlock(iPoint, Residual);
		if (implicit)
			Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
	}
  
}

void CLevelSetSolver::BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short val_marker) {
  /*--- Convective fluxes across euler wall are equal to zero. ---*/
}

void CLevelSetSolver::BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                                   CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex, Point_Normal;
	unsigned short iVar, iDim;
	double *V_domain, *V_wall, *V_mirror, *LevelSet_domain, *LevelSet_wall, *LevelSet_mirror;
	
	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	
	V_domain = new double[nDim+1];
	V_wall = new double[nDim+1];
  V_mirror = new double[nDim+1];
  
	LevelSet_domain = new double[1];
	LevelSet_wall = new double[1];
  LevelSet_mirror = new double[1];
  
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
    
		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {
      
			for (iVar = 0; iVar < nDim+1; iVar++) {
				V_domain[iVar] = solver_container[FLOW_SOL]->node[Point_Normal]->GetPrimVar(iVar);
				V_wall[iVar] = solver_container[FLOW_SOL]->node[iPoint]->GetPrimVar(iVar);
        V_mirror[iVar] = 2.0*V_wall[iVar] - V_domain[iVar];
			}
      
      LevelSet_domain[0] = node[iPoint]->GetSolution(0);
			LevelSet_wall[0] = node[iPoint]->GetSolution(0);
      LevelSet_mirror[0] = 2.0*LevelSet_wall[0] - LevelSet_domain[0];
      
			/*--- Set the normal vector ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
			for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = -Vector[iDim];
			conv_numerics->SetNormal(Vector);
			
			conv_numerics->SetPrimitive(V_wall, V_wall);
      
			conv_numerics->SetLevelSetVar(LevelSet_domain, LevelSet_wall);
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, Jacobian_MeanFlow_i, Jacobian_MeanFlow_j, config);
			
			LinSysRes.AddBlock(iPoint, Residual);
			
			if (implicit) {
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      }
			
		}
	}
  
  delete[] V_domain;
  delete[] V_wall;
  delete[] V_mirror;
  delete[] LevelSet_domain;
  delete[] LevelSet_wall;
  delete[] LevelSet_mirror;
}

void CLevelSetSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  /*--- Convective fluxes across euler wall are equal to zero. ---*/
}

void CLevelSetSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex;
	
	bool implicit = (config->GetKind_TimeIntScheme_LevelSet() == EULER_IMPLICIT);
	
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		
		/*--- Set the solution to the original value ---*/
		Solution[0] = node[iPoint]->GetSolution(0);
		
		node[iPoint]->SetSolution_Old(Solution);
		LinSysRes.SetBlock_Zero(iPoint);
		
		/*--- Includes 1 in the diagonal ---*/
		if (implicit)
			Jacobian.DeleteValsRowi(iPoint);
	}
}

void CLevelSetSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex, Point_Normal;
	unsigned short iVar, iDim;
	double *V_domain, *V_wall, *V_mirror, *LevelSet_domain, *LevelSet_wall, *LevelSet_mirror;
	
	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	
	V_domain = new double[nDim+1];
	V_wall = new double[nDim+1];
  V_mirror = new double[nDim+1];
  
	LevelSet_domain = new double[1];
	LevelSet_wall = new double[1];
  LevelSet_mirror = new double[1];
  
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
    
		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {
      
			for (iVar = 0; iVar < nDim+1; iVar++) {
				V_domain[iVar] = solver_container[FLOW_SOL]->node[Point_Normal]->GetPrimVar(iVar);
				V_wall[iVar] = solver_container[FLOW_SOL]->node[iPoint]->GetPrimVar(iVar);
        V_mirror[iVar] = 2.0*V_wall[iVar] - V_domain[iVar];
			}
      
      LevelSet_domain[0] = node[iPoint]->GetSolution(0);
			LevelSet_wall[0] = node[iPoint]->GetSolution(0);
      LevelSet_mirror[0] = 2.0*LevelSet_wall[0] - LevelSet_domain[0];
      
			/*--- Set the normal vector ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
			for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = -Vector[iDim];
			conv_numerics->SetNormal(Vector);
			
			conv_numerics->SetPrimitive(V_wall, V_mirror);
      
			conv_numerics->SetLevelSetVar(LevelSet_domain, LevelSet_mirror);
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, Jacobian_MeanFlow_i, Jacobian_MeanFlow_j, config);
			
			LinSysRes.AddBlock(iPoint, Residual);
			
			if (implicit) {
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      }
			
		}
	}
  
  delete[] V_domain;
  delete[] V_wall;
  delete[] V_mirror;
  delete[] LevelSet_domain;
  delete[] LevelSet_wall;
  delete[] LevelSet_mirror;
  
}

void CLevelSetSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex;
	
	bool implicit = (config->GetKind_TimeIntScheme_LevelSet() == EULER_IMPLICIT);
	
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		
		/*--- Set the solution to the original value ---*/
		Solution[0] = node[iPoint]->GetSolution(0);
		
		node[iPoint]->SetSolution_Old(Solution);
		LinSysRes.SetBlock_Zero(iPoint);
		
		/*--- Includes 1 in the diagonal ---*/
		if (implicit)
			Jacobian.DeleteValsRowi(iPoint);
	}
  
}

void CLevelSetSolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
	unsigned long iPoint;
	double Delta = 0.0, Vol;
	
	/*--- Set maximum residual to zero ---*/
	SetRes_RMS(0, 0.0);
  SetRes_Max(0, 0.0, 0);
  
	/*--- Build implicit system ---*/
	for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
		
		/*--- Read the volume ---*/
		Vol = geometry->node[iPoint]->GetVolume();
		
		/*--- Modify matrix diagonal to assure diagonal dominance ---*/
		Delta = Vol / (config->GetLevelSet_CFLRedCoeff()*solver_container[FLOW_SOL]->node[iPoint]->GetDelta_Time());
    
		Jacobian.AddVal2Diag(iPoint,Delta);
    
		/*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/
		LinSysRes[iPoint] = -LinSysRes[iPoint];
		LinSysSol[iPoint] = 0.0;
		AddRes_RMS(0, LinSysRes[iPoint]*LinSysRes[iPoint]);
    AddRes_Max(0, fabs(LinSysRes[iPoint]), geometry->node[iPoint]->GetGlobalIndex());
	}
	
  /*--- Initialize residual and solution at the ghost points ---*/
  for (iPoint = nPointDomain; iPoint < nPoint; iPoint++) {
    LinSysRes[iPoint] = 0.0;
    LinSysSol[iPoint] = 0.0;
  }
	
	/*--- Solve the linear system (Krylov subspace methods) ---*/
  CMatrixVectorProduct* mat_vec = new CSysMatrixVectorProduct(Jacobian, geometry, config);
  
  CPreconditioner* precond = NULL;
  if (config->GetKind_Linear_Solver_Prec() == JACOBI) {
    Jacobian.BuildJacobiPreconditioner();
    precond = new CJacobiPreconditioner(Jacobian, geometry, config);
  }
  else if (config->GetKind_Linear_Solver_Prec() == LU_SGS) {
    precond = new CLU_SGSPreconditioner(Jacobian, geometry, config);
  }
  else if (config->GetKind_Linear_Solver_Prec() == LINELET) {
    Jacobian.BuildJacobiPreconditioner();
    Jacobian.BuildLineletPreconditioner(geometry, config);
    precond = new CLineletPreconditioner(Jacobian, geometry, config);
  }
  
  CSysSolve system;
  if (config->GetKind_Linear_Solver() == BCGSTAB)
    system.BCGSTAB(LinSysRes, LinSysSol, *mat_vec, *precond, config->GetLinear_Solver_Error(),
                   config->GetLinear_Solver_Iter(), false);
  else if (config->GetKind_Linear_Solver() == FGMRES)
    system.FGMRES(LinSysRes, LinSysSol, *mat_vec, *precond, config->GetLinear_Solver_Error(),
                  config->GetLinear_Solver_Iter(), false);
  
  delete mat_vec;
  delete precond;
	
	/*--- Update solution (system written in terms of increments), be careful with the update of the
	 scalar equations which includes the density ---*/
	for (iPoint = 0; iPoint < nPointDomain; iPoint++)
		node[iPoint]->AddSolution(0, config->GetLinear_Solver_Relax()*LinSysSol[iPoint]);
  
  /*--- MPI solution ---*/
  Set_MPI_Solution(geometry, config);
  
  /*--- Compute the root mean square residual ---*/
  SetResidual_RMS(geometry, config);
  
}

void CLevelSetSolver::SetFreeSurface_Distance(CGeometry *geometry, CConfig *config, bool Initialization, bool WriteLevelSet) {
	double *coord = NULL, dist2, *iCoord = NULL, *jCoord = NULL, *U_i = NULL, *U_j = NULL,
  **Coord_LevelSet = NULL, *xCoord = NULL, *yCoord = NULL, *zCoord = NULL, auxCoordx, auxCoordy,
  auxCoordz, FreeSurface, volume, LevelSetDiff_Squared, LevelSetDiff, dist, Min_dist;
	unsigned short iDim;
	unsigned long iPoint, jPoint, iVertex, jVertex, nVertex_LevelSet, iEdge;
	ifstream index_file;
	ofstream LevelSet_file;
	string text_line;
	int rank = MASTER_NODE;
	char cstr[200], buffer[50];
  
	unsigned short nDim = geometry->GetnDim();
  unsigned long iExtIter = config->GetExtIter();
  
  /*--- Identification of the 0 level set points and coordinates ---*/
  nVertex_LevelSet = 0;
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    iPoint = geometry->edge[iEdge]->GetNode(0); U_i = node[iPoint]->GetSolution();
    jPoint = geometry->edge[iEdge]->GetNode(1); U_j = node[jPoint]->GetSolution();
    if (U_i[0]*U_j[0] < 0.0) nVertex_LevelSet ++;
  }
  
  /*--- Allocate vector of boundary coordinates ---*/
  Coord_LevelSet = new double* [nVertex_LevelSet];
  for (iVertex = 0; iVertex < nVertex_LevelSet; iVertex++)
    Coord_LevelSet[iVertex] = new double [nDim];
  
  /*--- Get coordinates of the points of the surface ---*/
  nVertex_LevelSet = 0;
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    iPoint = geometry->edge[iEdge]->GetNode(0); U_i = node[iPoint]->GetSolution(); iCoord = geometry->node[iPoint]->GetCoord();
    jPoint = geometry->edge[iEdge]->GetNode(1); U_j = node[jPoint]->GetSolution(); jCoord = geometry->node[jPoint]->GetCoord();
    if (U_i[0]*U_j[0] < 0.0) {
      for (iDim = 0; iDim < nDim; iDim++)
        Coord_LevelSet[nVertex_LevelSet][iDim] = iCoord[iDim]-U_i[0]*(jCoord[iDim]-iCoord[iDim])/(U_j[0]-U_i[0]);
      nVertex_LevelSet++;
    }
  }
  
  if (Initialization) {
    
    /*--- Get coordinates of the points and compute distances to the surface ---*/
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      coord = geometry->node[iPoint]->GetCoord();
      
      /*--- Compute the min distance ---*/
      Min_dist = 1E20;
      for (iVertex = 0; iVertex < nVertex_LevelSet; iVertex++) {
        
        dist2 = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          dist2 += (coord[iDim]-Coord_LevelSet[iVertex][iDim])*(coord[iDim]-Coord_LevelSet[iVertex][iDim]);
        dist = sqrt(dist2);
        if (dist < Min_dist) { Min_dist = dist; }
        
      }
      
      /*--- Compute the sign using the current solution ---*/
      double NumberSign = 1.0;
      if (node[iPoint]->GetSolution(0) != 0.0)
        NumberSign = node[iPoint]->GetSolution(0)/fabs(node[iPoint]->GetSolution(0));
      
      /*--- Store the value of the primitive variable ---*/
      node[iPoint]->SetPrimVar(0, Min_dist*NumberSign);
    }
    
  }
  
  if (WriteLevelSet) {
    
    /*--- Order the arrays (x Coordinate, y Coordinate, z Coordiante) ---*/
    for (iVertex = 0; iVertex < nVertex_LevelSet; iVertex++) {
      for (jVertex = 0; jVertex < nVertex_LevelSet - 1 - iVertex; jVertex++) {
        if (Coord_LevelSet[jVertex][0] > Coord_LevelSet[jVertex+1][0]) {
          auxCoordx = Coord_LevelSet[jVertex][0]; Coord_LevelSet[jVertex][0] = Coord_LevelSet[jVertex+1][0]; Coord_LevelSet[jVertex+1][0] = auxCoordx;
          auxCoordy = Coord_LevelSet[jVertex][1]; Coord_LevelSet[jVertex][1] = Coord_LevelSet[jVertex+1][1]; Coord_LevelSet[jVertex+1][1] = auxCoordy;
          if (nDim == 3) { auxCoordz = Coord_LevelSet[jVertex][2]; Coord_LevelSet[jVertex][2] = Coord_LevelSet[jVertex+1][2]; Coord_LevelSet[jVertex+1][2] = auxCoordz; }
        }
      }
    }
    
    /*--- Get coordinates of the points and compute distances to the surface ---*/
    FreeSurface = 0.0;
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      coord = geometry->node[iPoint]->GetCoord();
      volume = geometry->node[iPoint]->GetVolume();
      LevelSetDiff_Squared = 0.0; LevelSetDiff = 0.0;
      
      if ((coord[0] > 1.0) && (coord[0] < 4.0)){
        LevelSetDiff = (node[iPoint]->GetSolution()[0] - coord[nDim-1]);
        LevelSetDiff_Squared = LevelSetDiff*LevelSetDiff;
        FreeSurface += 0.5*LevelSetDiff_Squared*volume;
      }
      else {
        LevelSetDiff = 0.0;
      }
      
      node[iPoint]->SetDiffLevelSet(LevelSetDiff);
      
    }
    
    if ((rank == MASTER_NODE) && (iExtIter % config->GetWrt_Sol_Freq_DualTime() == 0)) {
      
      /*--- Write the Level Set distribution, the target level set---*/
      LevelSet_file.precision(15);
      
      /*--- Write file name with extension ---*/
      strcpy (cstr, "LevelSet");
      if (config->GetUnsteady_Simulation()){
        if ((int(iExtIter) >= 0) && (int(iExtIter) < 10)) sprintf (buffer, "_0000%d.dat", int(iExtIter));
        if ((int(iExtIter) >= 10) && (int(iExtIter) < 100)) sprintf (buffer, "_000%d.dat", int(iExtIter));
        if ((int(iExtIter) >= 100) && (int(iExtIter) < 1000)) sprintf (buffer, "_00%d.dat", int(iExtIter));
        if ((int(iExtIter) >= 1000) && (int(iExtIter) < 10000)) sprintf (buffer, "_0%d.dat", int(iExtIter));
        if (int(iExtIter) >= 10000) sprintf (buffer, "_%d.dat", int(iExtIter));
      }
      else {
        sprintf (buffer, ".dat");
      }
      
      strcat(cstr,buffer);
      
      LevelSet_file.open(cstr, ios::out);
      LevelSet_file << "TITLE = \"SU2 Free surface simulation\"" << endl;
      if (nDim == 2) LevelSet_file << "VARIABLES = \"x coord\",\"y coord\"" << endl;
      if (nDim == 3) LevelSet_file << "VARIABLES = \"x coord\",\"y coord\",\"z coord\"" << endl;
      LevelSet_file << "ZONE T= \"Free Surface\"" << endl;
      
      for (iVertex = 0; iVertex < nVertex_LevelSet; iVertex++) {
        if (nDim == 2) LevelSet_file << scientific << Coord_LevelSet[iVertex][0] << ", " << Coord_LevelSet[iVertex][1] << endl;
        if (nDim == 3) LevelSet_file << scientific << Coord_LevelSet[iVertex][0] << ", " << Coord_LevelSet[iVertex][1] << ", " << Coord_LevelSet[iVertex][2] << endl;
      }
      LevelSet_file.close();
      
    }
    
    /*--- Store the value of the free surface coefficient ---*/
    SetTotal_CFreeSurface(FreeSurface);
        
    delete [] xCoord;
    delete [] yCoord;
    if (nDim == 3) delete [] zCoord;
    
  }
  
  /*--- Deallocate vector of boundary coordinates ---*/
  for (iVertex = 0; iVertex < nVertex_LevelSet; iVertex++)
    delete Coord_LevelSet[iVertex];
  delete [] Coord_LevelSet;
  
}

void CLevelSetSolver::SetResidual_DualTime(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iRKStep, unsigned short iMesh, unsigned short RunTime_EqSystem) {
	unsigned long iPoint;
	double *U_time_nM1, *U_time_n, *U_time_nP1, Volume_nM1, Volume_n, Volume_nP1, TimeStep;
	
	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	bool Grid_Movement = config->GetGrid_Movement();
  
	/*--- loop over points ---*/
	for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
		
		/*--- Solution at time n-1, n and n+1 ---*/
		U_time_nM1  = node[iPoint]->GetSolution_time_n1();
		U_time_n    = node[iPoint]->GetSolution_time_n();
		U_time_nP1  = node[iPoint]->GetSolution();
    
		/*--- Volume at time n-1 and n ---*/
		if (Grid_Movement) {
			Volume_nM1 = geometry->node[iPoint]->GetVolume_nM1();
			Volume_n = geometry->node[iPoint]->GetVolume_n();
			Volume_nP1 = geometry->node[iPoint]->GetVolume();
		}
		else {
			Volume_nM1 = geometry->node[iPoint]->GetVolume();
			Volume_n = geometry->node[iPoint]->GetVolume();
			Volume_nP1 = geometry->node[iPoint]->GetVolume();
		}
		
		/*--- Time Step ---*/
		TimeStep = config->GetDelta_UnstTimeND();
		
		/*--- Compute Residual ---*/
    if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
      Residual[0] = ( U_time_nP1[0]*Volume_nP1 - U_time_n[0]*Volume_n ) / TimeStep;
    if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
      Residual[0] = ( 3.0*U_time_nP1[0]*Volume_nP1 - 4.0*U_time_n[0]*Volume_n
                     +  1.0*U_time_nM1[0]*Volume_nM1 ) / (2.0*TimeStep);
    
		/*--- Add Residual ---*/
		LinSysRes.AddBlock(iPoint, Residual);
		
		if (implicit) {
      Jacobian_i[0][0] = 0.0;
      if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
        Jacobian_i[0][0] = Volume_nP1 / TimeStep;
      if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
        Jacobian_i[0][0] = (Volume_nP1*3.0)/(2.0*TimeStep);
			Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
		}
    
	}
}
