/*!
 * \file geometry_structure.cpp
 * \brief Main subroutines for creating the primal grid and multigrid structure.
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

#include "../include/geometry_structure.hpp"

CGeometry::CGeometry(void) {
    
	nEdge = 0;
    nPoint = 0;
	nElem = 0;
    
	nElem_Bound = NULL;
	Tag_to_Marker = NULL;
	elem = NULL;
	face = NULL;
	bound = NULL;
	node = NULL;
	edge = NULL;
	vertex = NULL;
	nVertex = NULL;
	newBound = NULL;
	nNewElem_Bound = NULL;
    Marker_All_SendRecv = NULL;
    
    //	PeriodicPoint[MAX_NUMBER_PERIODIC][2].clear();
    //	PeriodicElem[MAX_NUMBER_PERIODIC].clear();
    //	OldBoundaryElems[MAX_NUMBER_MARKER].clear();
    //  SendTransfLocal[MAX_NUMBER_DOMAIN].clear();
    //  ReceivedTransfLocal[MAX_NUMBER_DOMAIN].clear();
    //	SendDomainLocal[MAX_NUMBER_DOMAIN].clear();
    //	ReceivedDomainLocal[MAX_NUMBER_DOMAIN].clear();
    //	XCoordList.clear();
    
    //	Xcoord_plane.clear();
    //	Ycoord_plane.clear();
    //	Zcoord_plane.clear();
    //	FaceArea_plane.clear();
    //	Plane_points.clear();
    
}

CGeometry::~CGeometry(void) {
    unsigned long iElem, iElem_Bound, iPoint, iFace, iVertex, iEdge;
    unsigned short iMarker;
    
    if (elem != NULL) {
        for (iElem = 0; iElem < nElem; iElem++)
            if (elem[iElem] != NULL) delete elem[iElem];
        delete[] elem;
    }
    
    if (bound != NULL) {
        for (iMarker = 0; iMarker < nMarker; iMarker++) {
            for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
                if (bound[iMarker][iElem_Bound] != NULL) delete bound[iMarker][iElem_Bound];
            }
        }
        delete[] bound;
    }
    
    if (face != NULL) {
        for (iFace = 0; iFace < nFace; iFace ++)
            if (face[iFace] != NULL) delete face[iFace];
        delete[] face;
    }
    
    if (node != NULL) {
        for (iPoint = 0; iPoint < nPoint; iPoint ++)
            if (node[iPoint] != NULL) delete node[iPoint];
        delete[] node;
    }
    
    if (edge != NULL) {
        for (iEdge = 0; iEdge < nEdge; iEdge ++)
            if (edge[iPoint] != NULL) delete edge[iEdge];
        delete[] edge;
    }
    
    if (vertex != NULL)  {
        for (iMarker = 0; iMarker < nMarker; iMarker++) {
            for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
                if (vertex[iMarker][iVertex] != NULL) delete vertex[iMarker][iVertex];
            }
        }
        delete[] vertex;
    }
    
    if (newBound != NULL) {
        for (iMarker = 0; iMarker < nMarker; iMarker++) {
            for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
                if (newBound[iMarker][iElem_Bound] != NULL) delete newBound[iMarker][iElem_Bound];
            }
        }
        delete[] newBound;
    }
    
	if (nElem_Bound != NULL) delete[] nElem_Bound;
	if (nVertex != NULL) delete[] nVertex;
	if (nNewElem_Bound != NULL) delete[] nNewElem_Bound;
    if (Marker_All_SendRecv != NULL) delete[] Marker_All_SendRecv;
	if (Tag_to_Marker != NULL) delete[] Tag_to_Marker;
    
    //	PeriodicPoint[MAX_NUMBER_PERIODIC][2].~vector();
    //	PeriodicElem[MAX_NUMBER_PERIODIC].~vector();
    //	OldBoundaryElems[MAX_NUMBER_MARKER].~vector();
    //  SendTransfLocal[MAX_NUMBER_DOMAIN].~vector();
    //  ReceivedTransfLocal[MAX_NUMBER_DOMAIN].~vector();
    //	SendDomainLocal[MAX_NUMBER_DOMAIN].~vector();
    //	ReceivedDomainLocal[MAX_NUMBER_DOMAIN].~vector();
    //	XCoordList.~vector();
    
    //	Xcoord_plane.~vector()
    //	Ycoord_plane.~vector()
    //	Zcoord_plane.~vector()
    //	FaceArea_plane.~vector()
    //	Plane_points.~vector()
    
}

double CGeometry::Point2Plane_Distance(double *Coord, double *iCoord, double *jCoord, double *kCoord) {
	double CrossProduct[3], iVector[3], jVector[3], distance, modulus;
	unsigned short iDim;
    
	for (iDim = 0; iDim < 3; iDim ++) {
		iVector[iDim] = jCoord[iDim] - iCoord[iDim];
		jVector[iDim] = kCoord[iDim] - iCoord[iDim];
	}
    
	CrossProduct[0] = iVector[1]*jVector[2] - iVector[2]*jVector[1];
	CrossProduct[1] = iVector[2]*jVector[0] - iVector[0]*jVector[2];
	CrossProduct[2] = iVector[0]*jVector[1] - iVector[1]*jVector[0];
    
	modulus = sqrt(CrossProduct[0]*CrossProduct[0]+CrossProduct[1]*CrossProduct[1]+CrossProduct[2]*CrossProduct[2]);
    
	distance = 0.0;
	for (iDim = 0; iDim < 3; iDim ++)
		distance += CrossProduct[iDim]*(Coord[iDim]-iCoord[iDim]);
	distance /= modulus;
    
	return distance;
    
}

long CGeometry::FindEdge(unsigned long first_point, unsigned long second_point) {
	unsigned long iPoint = 0;
	unsigned short iNode;
	for (iNode = 0; iNode < node[first_point]->GetnPoint(); iNode++) {
		iPoint = node[first_point]->GetPoint(iNode);
		if (iPoint == second_point) break;
	}
    
	if (iPoint == second_point) return node[first_point]->GetEdge(iNode);
	else {
        cout << "\n\n   !!! Error !!!\n" << endl;
		cout <<"Can't find the edge that connects "<< first_point <<" and "<< second_point <<"."<< endl;
		exit(1);
		return -1;
	}
}

bool CGeometry::CheckEdge(unsigned long first_point, unsigned long second_point) {
	unsigned long iPoint = 0;
	unsigned short iNode;
	for (iNode = 0; iNode < node[first_point]->GetnPoint(); iNode++) {
		iPoint = node[first_point]->GetPoint(iNode);
		if (iPoint == second_point) break;
	}
    
	if (iPoint == second_point) return true;
	else return false;
    
}

void CGeometry::SetEdges(void) {
	unsigned long iPoint, jPoint;
    long iEdge;
	unsigned short jNode, iNode;
	long TestEdge = 0;
    
	nEdge = 0;
	for(iPoint = 0; iPoint < nPoint; iPoint++)
		for(iNode = 0; iNode < node[iPoint]->GetnPoint(); iNode++) {
			jPoint = node[iPoint]->GetPoint(iNode);
			for(jNode = 0; jNode < node[jPoint]->GetnPoint(); jNode++)
				if (node[jPoint]->GetPoint(jNode) == iPoint) {
					TestEdge = node[jPoint]->GetEdge(jNode);
					break;
				}
			if (TestEdge == -1) {
				node[iPoint]->SetEdge(nEdge, iNode);
				node[jPoint]->SetEdge(nEdge, jNode);
				nEdge++;
			}
		}
    
	edge = new CEdge*[nEdge];
    
	for(iPoint = 0; iPoint < nPoint; iPoint++)
		for(iNode = 0; iNode < node[iPoint]->GetnPoint(); iNode++) {
			jPoint = node[iPoint]->GetPoint(iNode);
			iEdge = FindEdge(iPoint, jPoint);
			if (iPoint < jPoint) edge[iEdge] = new CEdge(iPoint, jPoint, nDim);
		}
}

void CGeometry::SetFaces(void) {
	//	unsigned long iPoint, jPoint, iFace;
	//	unsigned short jNode, iNode;
	//	long TestFace = 0;
	//
	//	nFace = 0;
	//	for(iPoint = 0; iPoint < nPoint; iPoint++)
	//		for(iNode = 0; iNode < node[iPoint]->GetnPoint(); iNode++) {
	//			jPoint = node[iPoint]->GetPoint(iNode);
	//			for(jNode = 0; jNode < node[jPoint]->GetnPoint(); jNode++)
	//				if (node[jPoint]->GetPoint(jNode) == iPoint) {
	//					TestFace = node[jPoint]->GetFace(jNode);
	//					break;
	//				}
	//			if (TestFace == -1) {
	//				node[iPoint]->SetFace(nFace, iNode);
	//				node[jPoint]->SetFace(nFace, jNode);
	//				nFace++;
	//			}
	//		}
	//
	//	face = new CFace*[nFace];
	//
	//	for(iPoint = 0; iPoint < nPoint; iPoint++)
	//		for(iNode = 0; iNode < node[iPoint]->GetnPoint(); iNode++) {
	//			jPoint = node[iPoint]->GetPoint(iNode);
	//			iFace = FindFace(iPoint, jPoint);
	//			if (iPoint < jPoint) face[iFace] = new CFace(iPoint,jPoint,nDim);
	//		}
}

void CGeometry::TestGeometry(void) {
    
	ofstream para_file;
    
	para_file.open("test_geometry.dat", ios::out);
    
	double *Normal = new double[nDim];
    
	for(unsigned long iEdge = 0; iEdge < nEdge; iEdge++) {
		para_file << "Edge index: " << iEdge << endl;
		para_file << "   Point index: " << edge[iEdge]->GetNode(0) << "\t" << edge[iEdge]->GetNode(1) << endl;
		edge[iEdge]->GetNormal(Normal);
		para_file << "      Face normal : ";
		for(unsigned short iDim = 0; iDim < nDim; iDim++)
			para_file << Normal[iDim] << "\t";
		para_file << endl;
	}
    
	para_file << endl;
	para_file << endl;
	para_file << endl;
	para_file << endl;
    
	for(unsigned short iMarker =0; iMarker < nMarker; iMarker++) {
		para_file << "Marker index: " << iMarker << endl;
		for(unsigned long iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
			para_file << "   Vertex index: " << iVertex << endl;
			para_file << "      Point index: " << vertex[iMarker][iVertex]->GetNode() << endl;
			para_file << "      Point coordinates : ";
			for(unsigned short iDim = 0; iDim < nDim; iDim++) {
				para_file << node[vertex[iMarker][iVertex]->GetNode()]->GetCoord(iDim) << "\t";}
			para_file << endl;
			vertex[iMarker][iVertex]->GetNormal(Normal);
			para_file << "         Face normal : ";
			for(unsigned short iDim = 0; iDim < nDim; iDim++)
				para_file << Normal[iDim] << "\t";
			para_file << endl;
		}
	}
    
}

void CGeometry::SetSpline(vector<double> &x, vector<double> &y, unsigned long n, double yp1, double ypn, vector<double> &y2) {
	unsigned long i, k;
	double p, qn, sig, un, *u;
    
	u = new double [n];
    
	if (yp1 > 0.99e30)			// The lower boundary condition is set either to be "nat
		y2[0]=u[0]=0.0;			  // -ural"
	else {									// or else to have a specified first derivative.
		y2[0] = -0.5;
		u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
	}
    
	for (i=2; i<=n-1; i++) {									//  This is the decomposition loop of the tridiagonal al-
		sig=(x[i-1]-x[i-2])/(x[i]-x[i-2]);		//	gorithm. y2 and u are used for tem-
		p=sig*y2[i-2]+2.0;										//	porary storage of the decomposed
		y2[i-1]=(sig-1.0)/p;										//	factors.
		u[i-1]=(y[i]-y[i-1])/(x[i]-x[i-1]) - (y[i-1]-y[i-2])/(x[i-1]-x[i-2]);
		u[i-1]=(6.0*u[i-1]/(x[i]-x[i-2])-sig*u[i-2])/p;
	}
    
	if (ypn > 0.99e30)						// The upper boundary condition is set either to be
		qn=un=0.0;									// "natural"
	else {												// or else to have a specified first derivative.
		qn=0.5;
		un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
	}
	y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
	for (k=n-1; k>=1; k--)					// This is the backsubstitution loop of the tridiagonal
		y2[k-1]=y2[k-1]*y2[k]+u[k-1];	  // algorithm.
    
	delete[] u;
    
}

double CGeometry::GetSpline(vector<double>&xa, vector<double>&ya, vector<double>&y2a, unsigned long n, double x) {
	unsigned long klo, khi, k;
	double h, b, a, y;
    
    if (x < xa[0]) x = xa[0];       // Clip max and min values
    if (x > xa[n-1]) x = xa[n-1];
    
	klo = 1;										// We will find the right place in the table by means of
	khi = n;										// bisection. This is optimal if sequential calls to this
	while (khi-klo > 1) {			// routine are at random values of x. If sequential calls
		k = (khi+klo) >> 1;				// are in order, and closely spaced, one would do better
		if (xa[k-1] > x) khi = k;		// to store previous values of klo and khi and test if
		else klo=k;							// they remain appropriate on the next call.
	}								// klo and khi now bracket the input value of x
	h = xa[khi-1] - xa[klo-1];
	if (h == 0.0) cout << "Bad xa input to routine splint" << endl;	// The xa’s must be distinct.
	a = (xa[khi-1]-x)/h;
	b = (x-xa[klo-1])/h;				// Cubic spline polynomial is now evaluated.
	y = a*ya[klo-1]+b*ya[khi-1]+((a*a*a-a)*y2a[klo-1]+(b*b*b-b)*y2a[khi-1])*(h*h)/6.0;
    
	return y;
}

unsigned short CGeometry::ComputeSegmentPlane_Intersection(double *Segment_P0, double *Segment_P1, double *Plane_P0, double *Plane_Normal, double *Intersection) {
    double u[3], v[3], Denominator, Numerator, Aux;
    unsigned short iDim;
    
    for (iDim = 0; iDim < 3; iDim++) {
        u[iDim] = Segment_P1[iDim] - Segment_P0[iDim];
        v[iDim] = Plane_P0[iDim] - Segment_P0[iDim];
    }
    
    Numerator = Plane_Normal[0]*v[0] + Plane_Normal[1]*v[1] + Plane_Normal[2]*v[2];
    Denominator = Plane_Normal[0]*u[0] + Plane_Normal[1]*u[1] + Plane_Normal[2]*u[2];
    
    if (fabs(Denominator) <= 0.0) return 0; // No intersection.
    
    Aux = Numerator / Denominator;
    
    if (Aux < 0.0 || Aux > 1.0) return 0; // No intersection.
    
    for (iDim = 0; iDim < 3; iDim++)
        Intersection[iDim] = Segment_P0[iDim] + Aux * u[iDim];
    
    
    /*--- Check that the intersection is in the segment ---*/
    for (iDim = 0; iDim < 3; iDim++) {
        u[iDim] = Segment_P0[iDim] - Intersection[iDim];
        v[iDim] = Segment_P1[iDim] - Intersection[iDim];
    }
    
    Denominator = Plane_Normal[0]*u[0] + Plane_Normal[1]*u[1] + Plane_Normal[2]*u[2];
    Numerator = Plane_Normal[0]*v[0] + Plane_Normal[1]*v[1] + Plane_Normal[2]*v[2];
    
    Aux = Numerator * Denominator;
    
    if (Aux > 0.0) return 3; // Intersection outside the segment.
    else return 1;
    
}

CPhysicalGeometry::CPhysicalGeometry() : CGeometry() {}

CPhysicalGeometry::CPhysicalGeometry(CConfig *config, unsigned short val_iZone, unsigned short val_nZone) : CGeometry() {
    
    /*--- Local variables and initialization ---*/
	string text_line, Marker_Tag;
	ifstream mesh_file;
	unsigned short iNode_Surface, iMarker;
	unsigned long Point_Surface, iElem_Surface;
	double Conversion_Factor = 1.0;
	int rank = MASTER_NODE;
	nZone = val_nZone;
    
    string val_mesh_filename = config->GetMesh_FileName();
    unsigned short val_format = config->GetMesh_FileFormat();
    
    /*--- Initialize counters for local/global points & elements ---*/
    
	if (rank == MASTER_NODE)
		cout << endl <<"---------------------- Read grid file information -----------------------" << endl;
    
	switch (val_format) {
        case SU2:
            Read_SU2_Format(config, val_mesh_filename, val_iZone, val_nZone);
            break;
        case NETCDF_ASCII:
            Read_NETCDF_Format(config, val_mesh_filename, val_iZone, val_nZone);
            break;
        default:
            cout << "Unrecognized mesh format specified!!" << endl;
            cout << "Press any key to exit..." << endl;
            cin.get();
            exit(1);
            break;
	}
    
    if (config->GetKind_SU2() == SU2_CFD) Conversion_Factor = config->GetConversion_Factor();
    else Conversion_Factor = 1.0;
    
	/*--- Loop over the surface element to set the boundaries ---*/
	for (iMarker = 0; iMarker < nMarker; iMarker++)
		for (iElem_Surface = 0; iElem_Surface < nElem_Bound[iMarker]; iElem_Surface++)
			for (iNode_Surface = 0; iNode_Surface < bound[iMarker][iElem_Surface]->GetnNodes(); iNode_Surface++) {
				Point_Surface = bound[iMarker][iElem_Surface]->GetNode(iNode_Surface);
				node[Point_Surface]->SetBoundary(nMarker);
                if (config->GetMarker_All_Boundary(iMarker) != SEND_RECEIVE &&
                    config->GetMarker_All_Boundary(iMarker) != INTERFACE_BOUNDARY &&
                    config->GetMarker_All_Boundary(iMarker) != NEARFIELD_BOUNDARY &&
                    config->GetMarker_All_Boundary(iMarker) != PERIODIC_BOUNDARY)
                    node[Point_Surface]->SetPhysicalBoundary(true);
            }
    
	/*--- Write a new copy of the grid in meters if requested ---*/
	if (config->GetKind_SU2() == SU2_CFD)
		if (config->GetWrite_Converted_Mesh()) {
			SetMeshFile(config,config->GetMesh_Out_FileName());
			cout.precision(4);
			cout << "Converted mesh by a factor of " << Conversion_Factor << endl;
			cout << "  and wrote to the output file: " << config->GetMesh_Out_FileName() << endl;
		}
    
}

CPhysicalGeometry::~CPhysicalGeometry(void) {
    
}

void CPhysicalGeometry::Read_SU2_Format(CConfig *config, string val_mesh_filename, unsigned short val_iZone, unsigned short val_nZone) {
    
    /*--- Local variables and initialization ---*/
	string text_line, Marker_Tag;
	ifstream mesh_file;
	unsigned short VTK_Type, iMarker, iChar, iCount = 0;
	unsigned long iElem_Bound = 0, iPoint = 0, ielem_div = 0, ielem = 0, *Local2Global = NULL, vnodes_edge[2], vnodes_triangle[3], vnodes_quad[4], vnodes_tetra[4], vnodes_hexa[8], vnodes_wedge[6], vnodes_pyramid[5], dummyLong, GlobalIndex, iElem;
	char cstr[200];
	double Coord_2D[2], Coord_3D[3], Conversion_Factor = 1.0, dummyDouble;
	string::size_type position;
	bool time_spectral = (config->GetUnsteady_Simulation() == TIME_SPECTRAL);
	int rank = MASTER_NODE, size = SINGLE_NODE;
	bool domain_flag = false;
	bool found_transform = false;
	nZone = val_nZone;
    
    /*--- Initialize counters for local/global points & elements ---*/
    FinestMGLevel = true;
	Global_nPoint = 0; Global_nPointDomain = 0; Global_nElem = 0;
    nelem_edge     = 0; Global_nelem_edge     = 0;
    nelem_triangle = 0; Global_nelem_triangle = 0;
    nelem_quad     = 0; Global_nelem_quad     = 0;
    nelem_tetra    = 0; Global_nelem_tetra    = 0;
    nelem_hexa     = 0; Global_nelem_hexa     = 0;
    nelem_wedge    = 0; Global_nelem_wedge    = 0;
    nelem_pyramid  = 0; Global_nelem_pyramid  = 0;
    
    /*--- Open grid file ---*/
    strcpy (cstr, val_mesh_filename.c_str());
    mesh_file.open(cstr, ios::in);
    
    /*--- Check the grid ---*/
    if (mesh_file.fail()) {
        cout << "There is no geometry file (CPhysicalGeometry)!! " << cstr << endl;
        cout << "Press any key to exit..." << endl;
        cin.get();
        exit(1);
    }
    
    /*--- If divided grid, we need the global index to
     perform the right element division, this is just a hack in the future we
     should read first the coordinates ---*/
    if (config->GetDivide_Element()) {
        
        unsigned long nDim_, nElem_, nPoint_;
        
        while (getline (mesh_file, text_line)) {
            
            position = text_line.find ("NDIME=",0);
            if (position != string::npos) {
                text_line.erase (0,6); nDim_ = atoi(text_line.c_str());
            }
            
            position = text_line.find ("NELEM=",0);
            if (position != string::npos) {
                text_line.erase (0,6); nElem_ = atoi(text_line.c_str());
                for (iElem = 0; iElem < nElem_; iElem ++)
                    getline(mesh_file, text_line);
            }
            
            position = text_line.find ("NPOIN=",0);
            if (position != string::npos) {
                text_line.erase (0,6);
                
                /*--- Now read the number of points and ghost points. ---*/
                stringstream  stream_line(text_line);
                stream_line >> nPoint_;
                
                Local2Global = new unsigned long [nPoint_];
                
                for (iPoint = 0; iPoint < nPoint_; iPoint ++) {
                    getline(mesh_file, text_line);
                    istringstream point_line(text_line);
                    if (size == 1) { Local2Global[iPoint] = iPoint; }
                    else {
                        point_line >> dummyDouble; point_line >> dummyDouble; if (nDim_ == 3) point_line >> dummyDouble;
                        point_line >> dummyLong; point_line >> Local2Global[iPoint];
                    }
                }
                
            }
        }
        
        /*--- Close and open again the grid file ---*/
        mesh_file.close();
        strcpy (cstr, val_mesh_filename.c_str());
        mesh_file.open(cstr, ios::in);
        
    }
    
    /*--- If more than one, find the zone in the mesh file ---*/
    if (val_nZone > 1 || time_spectral) {
        if (time_spectral) {
            if (rank == MASTER_NODE) cout << "Reading time spectral instance " << val_iZone << ":" << endl;
        } else {
            while (getline (mesh_file,text_line)) {
                /*--- Search for the current domain ---*/
                position = text_line.find ("IZONE=",0);
                if (position != string::npos) {
                    text_line.erase (0,6);
                    unsigned short jDomain = atoi(text_line.c_str());
                    if (jDomain == val_iZone) {
                        if (rank == MASTER_NODE) cout << "Reading zone " << val_iZone << ":" << endl;
                        break;
                    }
                }
            }
        }
    }
    
    /*--- Read grid file with format SU2 ---*/
    while (getline (mesh_file,text_line)) {
        
        /*--- Read the dimension of the problem ---*/
        position = text_line.find ("NDIME=",0);
        if (position != string::npos) {
            if (domain_flag == false) {
                text_line.erase (0,6); nDim = atoi(text_line.c_str());
                if (rank == MASTER_NODE) {
                    if (nDim == 2) cout << "Two dimensional problem." << endl;
                    if (nDim == 3) cout << "Three dimensional problem." << endl;
                }
                domain_flag = true;
            } else {
                break;
            }
        }
        
        /*--- Read the information about inner elements ---*/
        position = text_line.find ("NELEM=",0);
        if (position != string::npos) {
            text_line.erase (0,6); nElem = atoi(text_line.c_str());
            
            Global_nElem = nElem;
            if (rank == MASTER_NODE)
                cout << Global_nElem << " interior elements. " << endl;
            
            /*--- Allocate space for elements ---*/
            if (!config->GetDivide_Element()) elem = new CPrimalGrid*[nElem];
            else {
                if (nDim == 2) elem = new CPrimalGrid*[2*nElem];
                if (nDim == 3) elem = new CPrimalGrid*[6*nElem];
            }
            
            unsigned short IndirectionPrism[6][6], CT_FromVTK_Prism[6], IndirectionHexa[9][9];
            unsigned short temp, zero, one, two, three, four, five, six, seven, eight, iNode, smallestNode = 0, lookupindex;
            unsigned long smallest;
            
            /*--- Indirection matrix for dividing prisms into tets, using vtk format, and conversion table for local numbering of prisms in vtk format ---*/
            CT_FromVTK_Prism[0] = 1;  CT_FromVTK_Prism[1] = 3;  CT_FromVTK_Prism[2] = 2;  CT_FromVTK_Prism[3] = 4;  CT_FromVTK_Prism[4] = 6;  CT_FromVTK_Prism[5] = 5;
            
            IndirectionPrism[0][0] = 0; IndirectionPrism[0][1] = 2; IndirectionPrism[0][2] = 1; IndirectionPrism[0][3] = 3; IndirectionPrism[0][4] = 5; IndirectionPrism[0][5] = 4;
            IndirectionPrism[1][0] = 2; IndirectionPrism[1][1] = 1; IndirectionPrism[1][2] = 0; IndirectionPrism[1][3] = 5; IndirectionPrism[1][4] = 4; IndirectionPrism[1][5] = 3;
            IndirectionPrism[2][0] = 1; IndirectionPrism[2][1] = 0; IndirectionPrism[2][2] = 2; IndirectionPrism[2][3] = 4; IndirectionPrism[2][4] = 3; IndirectionPrism[2][5] = 5;
            IndirectionPrism[3][0] = 3; IndirectionPrism[3][1] = 4; IndirectionPrism[3][2] = 5; IndirectionPrism[3][3] = 0; IndirectionPrism[3][4] = 1; IndirectionPrism[3][5] = 2;
            IndirectionPrism[4][0] = 5; IndirectionPrism[4][1] = 3; IndirectionPrism[4][2] = 4; IndirectionPrism[4][3] = 2; IndirectionPrism[4][4] = 0; IndirectionPrism[4][5] = 1;
            IndirectionPrism[5][0] = 4; IndirectionPrism[5][1] = 5; IndirectionPrism[5][2] = 3; IndirectionPrism[5][3] = 1; IndirectionPrism[5][4] = 2; IndirectionPrism[5][5] = 0;
            
            /*--- Indirection matrix for dividing hexahedron into tets ---*/
            IndirectionHexa[1][1] = 1; IndirectionHexa[1][2] = 2; IndirectionHexa[1][3] = 3; IndirectionHexa[1][4] = 4; IndirectionHexa[1][5] = 5; IndirectionHexa[1][6] = 6; IndirectionHexa[1][7] = 7; IndirectionHexa[1][8] = 8;
            IndirectionHexa[2][1] = 2; IndirectionHexa[2][2] = 1; IndirectionHexa[2][3] = 5; IndirectionHexa[2][4] = 6; IndirectionHexa[2][5] = 3; IndirectionHexa[2][6] = 4; IndirectionHexa[2][7] = 8; IndirectionHexa[2][8] = 7;
            IndirectionHexa[3][1] = 3; IndirectionHexa[3][2] = 2; IndirectionHexa[3][3] = 6; IndirectionHexa[3][4] = 7; IndirectionHexa[3][5] = 4; IndirectionHexa[3][6] = 1; IndirectionHexa[3][7] = 5; IndirectionHexa[3][8] = 8;
            IndirectionHexa[4][1] = 4; IndirectionHexa[4][2] = 1; IndirectionHexa[4][3] = 2; IndirectionHexa[4][4] = 3; IndirectionHexa[4][5] = 8; IndirectionHexa[4][6] = 5; IndirectionHexa[4][7] = 6; IndirectionHexa[4][8] = 7;
            IndirectionHexa[5][1] = 5; IndirectionHexa[5][2] = 1; IndirectionHexa[5][3] = 4; IndirectionHexa[5][4] = 8; IndirectionHexa[5][5] = 6; IndirectionHexa[5][6] = 2; IndirectionHexa[5][7] = 3; IndirectionHexa[5][8] = 7;
            IndirectionHexa[6][1] = 6; IndirectionHexa[6][2] = 2; IndirectionHexa[6][3] = 1; IndirectionHexa[6][4] = 5; IndirectionHexa[6][5] = 7; IndirectionHexa[6][6] = 3; IndirectionHexa[6][7] = 4; IndirectionHexa[6][8] = 8;
            IndirectionHexa[7][1] = 7; IndirectionHexa[7][2] = 3; IndirectionHexa[7][3] = 2; IndirectionHexa[7][4] = 6; IndirectionHexa[7][5] = 8; IndirectionHexa[7][6] = 4; IndirectionHexa[7][7] = 1; IndirectionHexa[7][8] = 5;
            IndirectionHexa[8][1] = 8; IndirectionHexa[8][2] = 4; IndirectionHexa[8][3] = 3; IndirectionHexa[8][4] = 7; IndirectionHexa[8][5] = 5; IndirectionHexa[8][6] = 1; IndirectionHexa[8][7] = 2; IndirectionHexa[8][8] = 6;
            
            
            /*--- Loop over all the volumetric elements ---*/
            while (ielem_div < nElem) {
                getline(mesh_file,text_line);
                istringstream elem_line(text_line);
                
                elem_line >> VTK_Type;
                
                switch(VTK_Type) {
                    case TRIANGLE:
                        
                        elem_line >> vnodes_triangle[0]; elem_line >> vnodes_triangle[1]; elem_line >> vnodes_triangle[2];
                        elem[ielem] = new CTriangle(vnodes_triangle[0], vnodes_triangle[1], vnodes_triangle[2], 2);
                        ielem_div++; ielem++; nelem_triangle++;
                        break;
                        
                    case RECTANGLE:
                        
                        elem_line >> vnodes_quad[0]; elem_line >> vnodes_quad[1]; elem_line >> vnodes_quad[2]; elem_line >> vnodes_quad[3];
                        if (!config->GetDivide_Element()) {
                            elem[ielem] = new CRectangle(vnodes_quad[0], vnodes_quad[1], vnodes_quad[2], vnodes_quad[3], 2);
                            ielem++; nelem_quad++; }
                        else {
                            
                            unsigned long index1, index2;
                            bool division1 = false;
                            index1 = min(Local2Global[vnodes_quad[0]], Local2Global[vnodes_quad[2]]);
                            index2 = min(Local2Global[vnodes_quad[1]], Local2Global[vnodes_quad[3]]);
                            if (index1 < index2) division1 = true;
                            
                            if (division1) {
                                elem[ielem] = new CTriangle(vnodes_quad[0], vnodes_quad[2], vnodes_quad[1], 2);
                                ielem++; nelem_triangle++;
                                elem[ielem] = new CTriangle(vnodes_quad[3], vnodes_quad[2], vnodes_quad[0], 2);
                                ielem++; nelem_triangle++;
                            }
                            else {
                                elem[ielem] = new CTriangle(vnodes_quad[2], vnodes_quad[1], vnodes_quad[3], 2);
                                ielem++; nelem_triangle++;
                                elem[ielem] = new CTriangle(vnodes_quad[1], vnodes_quad[0], vnodes_quad[3], 2);
                                ielem++; nelem_triangle++;
                            }
                            
                        }
                        
                        ielem_div++;
                        break;
                        
                    case TETRAHEDRON:
                        
                        elem_line >> vnodes_tetra[0]; elem_line >> vnodes_tetra[1]; elem_line >> vnodes_tetra[2]; elem_line >> vnodes_tetra[3];
                        elem[ielem] = new CTetrahedron(vnodes_tetra[0], vnodes_tetra[1], vnodes_tetra[2], vnodes_tetra[3]);
                        ielem_div++; ielem++; nelem_tetra++;
                        break;
                        
                    case HEXAHEDRON:
                        
                        elem_line >> vnodes_hexa[0]; elem_line >> vnodes_hexa[1]; elem_line >> vnodes_hexa[2];
                        elem_line >> vnodes_hexa[3]; elem_line >> vnodes_hexa[4]; elem_line >> vnodes_hexa[5];
                        elem_line >> vnodes_hexa[6]; elem_line >> vnodes_hexa[7];
                        
                        if (!config->GetDivide_Element()) {
                            elem[ielem] = new CHexahedron(vnodes_hexa[0], vnodes_hexa[1], vnodes_hexa[2], vnodes_hexa[3],
                                                          vnodes_hexa[4], vnodes_hexa[5], vnodes_hexa[6], vnodes_hexa[7]);
                            ielem++; nelem_hexa++;
                        }
                        else {
                            
                            smallest = Local2Global[vnodes_hexa[0]]; smallestNode = 0;
                            for (iNode = 1; iNode < 8; iNode ++) {
                                if ( smallest > Local2Global[vnodes_hexa[iNode]]) {
                                    smallest = Local2Global[vnodes_hexa[iNode]];
                                    smallestNode = iNode;
                                }
                            }
                            
                            one  = IndirectionHexa[smallestNode+1][1] - 1;
                            two  = IndirectionHexa[smallestNode+1][2] - 1;
                            three  = IndirectionHexa[smallestNode+1][3] - 1;
                            four = IndirectionHexa[smallestNode+1][4] - 1;
                            five = IndirectionHexa[smallestNode+1][5] - 1;
                            six = IndirectionHexa[smallestNode+1][6] - 1;
                            seven = IndirectionHexa[smallestNode+1][7] - 1;
                            eight = IndirectionHexa[smallestNode+1][8] - 1;
                            
                            unsigned long index1, index2;
                            unsigned short code1 = 0, code2 = 0, code3 = 0;
                            
                            index1 = min(Local2Global[vnodes_hexa[two]], Local2Global[vnodes_hexa[seven]]);
                            index2 = min(Local2Global[vnodes_hexa[three]], Local2Global[vnodes_hexa[six]]);
                            if (index1 < index2) code1 = 1;
                            
                            index1 = min(Local2Global[vnodes_hexa[four]], Local2Global[vnodes_hexa[seven]]);
                            index2 = min(Local2Global[vnodes_hexa[three]], Local2Global[vnodes_hexa[eight]]);
                            if (index1 < index2) code2 = 1;
                            
                            index1 = min(Local2Global[vnodes_hexa[five]], Local2Global[vnodes_hexa[seven]]);
                            index2 = min(Local2Global[vnodes_hexa[six]], Local2Global[vnodes_hexa[eight]]);
                            if (index1 < index2) code3 = 1;
                            
                            /*--- Rotation of 120 degrees ---*/
                            if ((!code1 && !code2 && code3) || (code1 && code2 && !code3)) {
                                temp = two; two = five; five = four; four = temp;
                                temp = six; six = eight; eight = three; three = temp;
                            }
                            
                            /*--- Rotation of 240 degrees ---*/
                            if ((!code1 && code2 && !code3) || (code1 && !code2 && code3)) {
                                temp = two; two = four; four = five; five = temp;
                                temp = six; six = three; three = eight; eight = temp;
                            }
                            
                            if ((code1 + code2 + code3) == 0) {
                                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[two], vnodes_hexa[three], vnodes_hexa[six]);
                                ielem++; nelem_tetra++;
                                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[three], vnodes_hexa[eight], vnodes_hexa[six]);
                                ielem++; nelem_tetra++;
                                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[three], vnodes_hexa[four], vnodes_hexa[eight]);
                                ielem++; nelem_tetra++;
                                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[six], vnodes_hexa[eight], vnodes_hexa[five]);
                                ielem++; nelem_tetra++;
                                elem[ielem] = new CTetrahedron(vnodes_hexa[three], vnodes_hexa[eight], vnodes_hexa[six], vnodes_hexa[seven]);
                                ielem++; nelem_tetra++;
                            }
                            if ((code1 + code2 + code3) == 1) {
                                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[six], vnodes_hexa[eight], vnodes_hexa[five]);
                                ielem++; nelem_tetra++;
                                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[two], vnodes_hexa[eight], vnodes_hexa[six]);
                                ielem++; nelem_tetra++;
                                elem[ielem] = new CTetrahedron(vnodes_hexa[two], vnodes_hexa[seven], vnodes_hexa[eight], vnodes_hexa[six]);
                                ielem++; nelem_tetra++;
                                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[eight], vnodes_hexa[three], vnodes_hexa[four]);
                                ielem++; nelem_tetra++;
                                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[eight], vnodes_hexa[two], vnodes_hexa[three]);
                                ielem++; nelem_tetra++;
                                elem[ielem] = new CTetrahedron(vnodes_hexa[two], vnodes_hexa[eight], vnodes_hexa[seven], vnodes_hexa[three]);
                                ielem++; nelem_tetra++;
                                
                                //                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[six], vnodes_hexa[eight], vnodes_hexa[five]);
                                //                ielem++; nelem_tetra++;
                                //                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[two], vnodes_hexa[seven], vnodes_hexa[six]);
                                //                ielem++; nelem_tetra++;
                                //                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[seven], vnodes_hexa[eight], vnodes_hexa[six]);
                                //                ielem++; nelem_tetra++;
                                //                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[eight], vnodes_hexa[three], vnodes_hexa[four]);
                                //                ielem++; nelem_tetra++;
                                //                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[eight], vnodes_hexa[seven], vnodes_hexa[three]);
                                //                ielem++; nelem_tetra++;
                                //                elem[ielem] = new CTetrahedron(vnodes_hexa[two], vnodes_hexa[one], vnodes_hexa[seven], vnodes_hexa[three]);
                                //                ielem++; nelem_tetra++;
                                
                            }
                            if ((code1 + code2 + code3) == 2) {
                                //                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[five], vnodes_hexa[six], vnodes_hexa[seven]);
                                //                ielem++; nelem_tetra++;
                                //                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[four], vnodes_hexa[eight], vnodes_hexa[seven]);
                                //                ielem++; nelem_tetra++;
                                //                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[eight], vnodes_hexa[five], vnodes_hexa[seven]);
                                //                ielem++; nelem_tetra++;
                                //                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[two], vnodes_hexa[three], vnodes_hexa[six]);
                                //                ielem++; nelem_tetra++;
                                //                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[four], vnodes_hexa[seven], vnodes_hexa[three]);
                                //                ielem++; nelem_tetra++;
                                //                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[seven], vnodes_hexa[six], vnodes_hexa[three]);
                                //                ielem++; nelem_tetra++;
                                
                                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[three], vnodes_hexa[four], vnodes_hexa[seven]);
                                ielem++; nelem_tetra++;
                                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[five], vnodes_hexa[seven], vnodes_hexa[eight]);
                                ielem++; nelem_tetra++;
                                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[seven], vnodes_hexa[four], vnodes_hexa[eight]);
                                ielem++; nelem_tetra++;
                                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[two], vnodes_hexa[three], vnodes_hexa[six]);
                                ielem++; nelem_tetra++;
                                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[seven], vnodes_hexa[five], vnodes_hexa[six]);
                                ielem++; nelem_tetra++;
                                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[three], vnodes_hexa[seven], vnodes_hexa[six]);
                                ielem++; nelem_tetra++;
                            }
                            if ((code1 + code2 + code3) == 3) {
                                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[three], vnodes_hexa[four], vnodes_hexa[seven]);
                                ielem++; nelem_tetra++;
                                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[four], vnodes_hexa[eight], vnodes_hexa[seven]);
                                ielem++; nelem_tetra++;
                                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[eight], vnodes_hexa[five], vnodes_hexa[seven]);
                                ielem++; nelem_tetra++;
                                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[six], vnodes_hexa[seven], vnodes_hexa[five]);
                                ielem++; nelem_tetra++;
                                elem[ielem] = new CTetrahedron(vnodes_hexa[two], vnodes_hexa[six], vnodes_hexa[seven], vnodes_hexa[one]);
                                ielem++; nelem_tetra++;
                                elem[ielem] = new CTetrahedron(vnodes_hexa[two], vnodes_hexa[seven], vnodes_hexa[three], vnodes_hexa[one]);
                                ielem++; nelem_tetra++;
                                
                                //                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[two], vnodes_hexa[seven], vnodes_hexa[six]);
                                //                ielem++; nelem_tetra++;
                                //                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[seven], vnodes_hexa[eight], vnodes_hexa[five]);
                                //                ielem++; nelem_tetra++;
                                //                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[six], vnodes_hexa[seven], vnodes_hexa[five]);
                                //                ielem++; nelem_tetra++;
                                //                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[seven], vnodes_hexa[two], vnodes_hexa[three]);
                                //                ielem++; nelem_tetra++;
                                //                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[eight], vnodes_hexa[seven], vnodes_hexa[four]);
                                //                ielem++; nelem_tetra++;
                                //                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[seven], vnodes_hexa[three], vnodes_hexa[four]);
                                //                ielem++; nelem_tetra++;
                                
                                //                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[two], vnodes_hexa[seven], vnodes_hexa[six]);
                                //                ielem++; nelem_tetra++;
                                //                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[two], vnodes_hexa[three], vnodes_hexa[seven]);
                                //                ielem++; nelem_tetra++;
                                //                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[three], vnodes_hexa[four], vnodes_hexa[seven]);
                                //                ielem++; nelem_tetra++;
                                //                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[six], vnodes_hexa[seven], vnodes_hexa[five]);
                                //                ielem++; nelem_tetra++;
                                //                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[seven], vnodes_hexa[four], vnodes_hexa[eight]);
                                //                ielem++; nelem_tetra++;
                                //                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[seven], vnodes_hexa[eight], vnodes_hexa[five]);
                                //                ielem++; nelem_tetra++;
                            }
                            
                        }
                        ielem_div++;
                        break;
                        
                    case WEDGE:
                        
                        elem_line >> vnodes_wedge[0]; elem_line >> vnodes_wedge[1]; elem_line >> vnodes_wedge[2];
                        elem_line >> vnodes_wedge[3]; elem_line >> vnodes_wedge[4]; elem_line >> vnodes_wedge[5];
                        
                        if (!config->GetDivide_Element()) {
                            
                            elem[ielem] = new CWedge(vnodes_wedge[0],vnodes_wedge[1],vnodes_wedge[2],vnodes_wedge[3],vnodes_wedge[4],vnodes_wedge[5]);
                            ielem++; nelem_wedge++;
                            
                        }
                        else {
                            
                            smallest = Local2Global[vnodes_wedge[0]]; smallestNode = 0;
                            for (iNode = 1; iNode < 6; iNode ++) {
                                if ( smallest > Local2Global[vnodes_wedge[iNode]]) {
                                    smallest = Local2Global[vnodes_wedge[iNode]];
                                    smallestNode = iNode;
                                }
                            }
                            
                            lookupindex = (CT_FromVTK_Prism[smallestNode] - 1);
                            zero  = IndirectionPrism[lookupindex][0];
                            one  = IndirectionPrism[lookupindex][1];
                            two  = IndirectionPrism[lookupindex][2];
                            three = IndirectionPrism[lookupindex][3];
                            four = IndirectionPrism[lookupindex][4];
                            five = IndirectionPrism[lookupindex][5];
                            
                            unsigned long index1, index2;
                            bool division = false;
                            index1 = min(Local2Global[vnodes_wedge[one]], Local2Global[vnodes_wedge[five]]);
                            index2 = min(Local2Global[vnodes_wedge[two]], Local2Global[vnodes_wedge[four]]);
                            if (index1 < index2) division = true;
                            
                            if (division) {
                                elem[ielem] = new CTetrahedron(vnodes_wedge[zero], vnodes_wedge[one], vnodes_wedge[two], vnodes_wedge[five]);
                                ielem++; nelem_tetra++;
                                elem[ielem] = new CTetrahedron(vnodes_wedge[zero], vnodes_wedge[one], vnodes_wedge[five], vnodes_wedge[four]);
                                ielem++; nelem_tetra++;
                                elem[ielem] = new CTetrahedron(vnodes_wedge[zero], vnodes_wedge[four], vnodes_wedge[five], vnodes_wedge[three]);
                                ielem++; nelem_tetra++;
                            }
                            else {
                                elem[ielem] = new CTetrahedron(vnodes_wedge[zero], vnodes_wedge[one], vnodes_wedge[two], vnodes_wedge[four]);
                                ielem++; nelem_tetra++;
                                elem[ielem] = new CTetrahedron(vnodes_wedge[zero], vnodes_wedge[four], vnodes_wedge[two], vnodes_wedge[five]);
                                ielem++; nelem_tetra++;
                                elem[ielem] = new CTetrahedron(vnodes_wedge[zero], vnodes_wedge[four], vnodes_wedge[five], vnodes_wedge[three]);
                                ielem++; nelem_tetra++;
                            }
                        }
                        
                        ielem_div++;
                        break;
                    case PYRAMID:
                        
                        elem_line >> vnodes_pyramid[0]; elem_line >> vnodes_pyramid[1]; elem_line >> vnodes_pyramid[2];
                        elem_line >> vnodes_pyramid[3]; elem_line >> vnodes_pyramid[4];
                        
                        if (!config->GetDivide_Element()) {
                            
                            elem[ielem] = new CPyramid(vnodes_pyramid[0],vnodes_pyramid[1],vnodes_pyramid[2],vnodes_pyramid[3],vnodes_pyramid[4]);
                            ielem++; nelem_pyramid++;
                            
                        }
                        else {
                            
                            unsigned long index1, index2;
                            bool division = false;
                            index1 = min(Local2Global[vnodes_pyramid[0]], Local2Global[vnodes_pyramid[2]]);
                            index2 = min(Local2Global[vnodes_pyramid[1]], Local2Global[vnodes_pyramid[3]]);
                            if (index1 < index2) division = true;
                            
                            if (division) {
                                elem[ielem] = new CTetrahedron(vnodes_pyramid[0], vnodes_pyramid[1], vnodes_pyramid[2], vnodes_pyramid[4]);
                                ielem++; nelem_tetra++;
                                elem[ielem] = new CTetrahedron(vnodes_pyramid[0], vnodes_pyramid[2], vnodes_pyramid[3], vnodes_pyramid[4]);
                                ielem++; nelem_tetra++;
                            }
                            else {
                                elem[ielem] = new CTetrahedron(vnodes_pyramid[1], vnodes_pyramid[2], vnodes_pyramid[3], vnodes_pyramid[4]);
                                ielem++; nelem_tetra++;
                                elem[ielem] = new CTetrahedron(vnodes_pyramid[1], vnodes_pyramid[3], vnodes_pyramid[0], vnodes_pyramid[4]);
                                ielem++; nelem_tetra++;
                            }
                            
                        }
                        
                        ielem_div++;
                        break;
                }
            }
            
            if (config->GetDivide_Element()) nElem = nelem_triangle + nelem_quad + nelem_tetra + nelem_hexa + nelem_wedge + nelem_pyramid;
            
            /*--- Communicate the number of each element type to all processors. ---*/
            Global_nelem_triangle = nelem_triangle;
            Global_nelem_quad     = nelem_quad;
            Global_nelem_tetra    = nelem_tetra;
            Global_nelem_hexa     = nelem_hexa;
            Global_nelem_wedge    = nelem_wedge;
            Global_nelem_pyramid  = nelem_pyramid;
            
            /*--- Print information about the elements to the console ---*/
            if (rank == MASTER_NODE) {
                if (Global_nelem_triangle > 0)
                    cout << Global_nelem_triangle << " triangles." << endl;
                if (Global_nelem_quad > 0)
                    cout << Global_nelem_quad << " quadrilaterals." << endl;
                if (Global_nelem_tetra > 0)
                    cout << Global_nelem_tetra << " tetrahedra." << endl;
                if (Global_nelem_hexa > 0)
                    cout << Global_nelem_hexa << " hexahedra." << endl;
                if (Global_nelem_wedge > 0)
                    cout << Global_nelem_wedge << " prisms." << endl;
                if (Global_nelem_pyramid > 0)
                    cout << Global_nelem_pyramid << " pyramids." << endl;
                if ((size > SINGLE_NODE) && (config->GetKind_SU2() != SU2_DDC))
                    cout << "Element totals include halo cells." << endl;
            }
        }
        
        /*--- Read number of points ---*/
        position = text_line.find ("NPOIN=",0);
        if (position != string::npos) {
            text_line.erase (0,6);
            
            /*--- Check for ghost points. ---*/
            stringstream test_line(text_line);
            while (test_line >> dummyLong)
                iCount++;
            
            /*--- Now read and store the number of points and possible ghost points. ---*/
            stringstream  stream_line(text_line);
            if (iCount == 2) {
                stream_line >> nPoint;
                stream_line >> nPointDomain;
                
                /*--- Set some important point information for parallel simulations. ---*/
                Global_nPoint = nPoint;
                Global_nPointDomain = nPointDomain;
                if (rank == MASTER_NODE)
                    cout << Global_nPointDomain << " points, and " << Global_nPoint-Global_nPointDomain << " ghost points." << endl;
            }
            else if (iCount == 1) {
                stream_line >> nPoint;
                nPointDomain = nPoint;
                Global_nPointDomain = nPoint;
                Global_nPoint = nPoint;
                if (rank == MASTER_NODE) cout << nPoint << " points." << endl;
            }
            else {
                cout << "NPOIN improperly specified!!" << endl;
                cout << "Press any key to exit..." << endl;
                cin.get();
                exit(1);
            }
            
            /*--- Retrieve grid conversion factor. The conversion is only
             applied for SU2_CFD. All other SU2 components leave the mesh
             as is. ---*/
            if (config->GetKind_SU2() == SU2_CFD)
                Conversion_Factor = config->GetConversion_Factor();
            else
                Conversion_Factor = 1.0;
            
            node = new CPoint*[nPoint];
            iPoint = 0;
            while (iPoint < nPoint) {
                getline(mesh_file,text_line);
                istringstream point_line(text_line);
                switch(nDim) {
                    case 2:
                        GlobalIndex = iPoint;
                        point_line >> Coord_2D[0]; point_line >> Coord_2D[1];
                        node[iPoint] = new CPoint(Conversion_Factor*Coord_2D[0], Conversion_Factor*Coord_2D[1], GlobalIndex, config);
                        iPoint++; break;
                    case 3:
                        GlobalIndex = iPoint;
                        point_line >> Coord_3D[0]; point_line >> Coord_3D[1]; point_line >> Coord_3D[2];
                        node[iPoint] = new CPoint(Conversion_Factor*Coord_3D[0], Conversion_Factor*Coord_3D[1], Conversion_Factor*Coord_3D[2], GlobalIndex, config);
                        iPoint++; break;
                }
            }
        }
        
        
        /*--- Read number of markers ---*/
        position = text_line.find ("NMARK=",0);
        if (position != string::npos) {
            text_line.erase (0,6); nMarker = atoi(text_line.c_str());
            if (size == 1) cout << nMarker << " surface markers." << endl;
            config->SetnMarker_All(nMarker);
            bound = new CPrimalGrid**[nMarker];
            nElem_Bound = new unsigned long [nMarker];
            Tag_to_Marker = new string [MAX_INDEX_VALUE];
            
            for (iMarker = 0 ; iMarker < nMarker; iMarker++) {
                getline (mesh_file,text_line);
                text_line.erase (0,11);
                string::size_type position;
                for (iChar = 0; iChar < 20; iChar++) {
                    position = text_line.find( " ", 0 );
                    if(position != string::npos) text_line.erase (position,1);
                    position = text_line.find( "\r", 0 );
                    if(position != string::npos) text_line.erase (position,1);
                    position = text_line.find( "\n", 0 );
                    if(position != string::npos) text_line.erase (position,1);
                }
                Marker_Tag = text_line.c_str();
                
                /*--- Physical boundaries definition ---*/
                if (Marker_Tag != "SEND_RECEIVE") {
                    getline (mesh_file,text_line);
                    text_line.erase (0,13); nElem_Bound[iMarker] = atoi(text_line.c_str());
                    if (size == 1)
                        cout << nElem_Bound[iMarker]  << " boundary elements in index "<< iMarker <<" (Marker = " <<Marker_Tag<< ")." << endl;
                    
                    
                    /*--- Allocate space for elements ---*/
                    if (!config->GetDivide_Element()) bound[iMarker] = new CPrimalGrid* [nElem_Bound[iMarker]];
                    else {
                        if (nDim == 2) bound[iMarker] = new CPrimalGrid* [2*nElem_Bound[iMarker]];;
                        if (nDim == 3) bound[iMarker] = new CPrimalGrid* [2*nElem_Bound[iMarker]];;
                    }
                    
                    nelem_edge_bound = 0; nelem_triangle_bound = 0; nelem_quad_bound = 0; ielem = 0;
                    for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
                        getline(mesh_file,text_line);
                        istringstream bound_line(text_line);
                        bound_line >> VTK_Type;
                        switch(VTK_Type) {
                            case LINE:
                                
                                if (nDim == 3) {
                                    cout << "Please remove line boundary conditions from the mesh file!" << endl;
                                    cout << "Press any key to exit..." << endl;
                                    cin.get();
                                    exit(1);
                                }
                                
                                bound_line >> vnodes_edge[0]; bound_line >> vnodes_edge[1];
                                bound[iMarker][ielem] = new CLine(vnodes_edge[0],vnodes_edge[1],2);
                                ielem++; nelem_edge_bound++; break;
                                
                            case TRIANGLE:
                                bound_line >> vnodes_triangle[0]; bound_line >> vnodes_triangle[1]; bound_line >> vnodes_triangle[2];
                                bound[iMarker][ielem] = new CTriangle(vnodes_triangle[0],vnodes_triangle[1],vnodes_triangle[2],3);
                                ielem++; nelem_triangle_bound++; break;
                                
                            case RECTANGLE:
                                
                                bound_line >> vnodes_quad[0]; bound_line >> vnodes_quad[1]; bound_line >> vnodes_quad[2]; bound_line >> vnodes_quad[3];
                                
                                if (!config->GetDivide_Element()) {
                                    
                                    bound[iMarker][ielem] = new CRectangle(vnodes_quad[0],vnodes_quad[1],vnodes_quad[2],vnodes_quad[3],3);
                                    ielem++; nelem_quad_bound++;
                                    
                                }
                                else {
                                    
                                    unsigned long index1, index2;
                                    bool division1 = false;
                                    index1 = min(Local2Global[vnodes_quad[0]], Local2Global[vnodes_quad[2]]);
                                    index2 = min(Local2Global[vnodes_quad[1]], Local2Global[vnodes_quad[3]]);
                                    if (index1 < index2) division1 = true;
                                    
                                    if (division1) {
                                        bound[iMarker][ielem] = new CTriangle(vnodes_quad[0], vnodes_quad[2], vnodes_quad[1], 3);
                                        ielem++; nelem_triangle_bound++;
                                        bound[iMarker][ielem] = new CTriangle(vnodes_quad[3], vnodes_quad[2], vnodes_quad[0], 3);
                                        ielem++; nelem_triangle_bound++;
                                    }
                                    else {
                                        bound[iMarker][ielem] = new CTriangle(vnodes_quad[2], vnodes_quad[1], vnodes_quad[3], 3);
                                        ielem++; nelem_triangle_bound++;
                                        bound[iMarker][ielem] = new CTriangle(vnodes_quad[1], vnodes_quad[0], vnodes_quad[3], 3);
                                        ielem++; nelem_triangle_bound++;
                                    }
                                    
                                }
                                
                                break;
                                
                                
                        }
                    }
                    if (config->GetDivide_Element()) nElem_Bound[iMarker] = nelem_edge_bound + nelem_triangle_bound + nelem_quad_bound;
                    
                    /*--- Update config information storing the boundary information in the right place ---*/
                    Tag_to_Marker[config->GetMarker_Config_Tag(Marker_Tag)] = Marker_Tag;
                    config->SetMarker_All_Tag(iMarker, Marker_Tag);
                    config->SetMarker_All_Boundary(iMarker, config->GetMarker_Config_Boundary(Marker_Tag));
                    config->SetMarker_All_Monitoring(iMarker, config->GetMarker_Config_Monitoring(Marker_Tag));
                    config->SetMarker_All_Designing(iMarker, config->GetMarker_Config_Designing(Marker_Tag));
                    config->SetMarker_All_Plotting(iMarker, config->GetMarker_Config_Plotting(Marker_Tag));
                    config->SetMarker_All_DV(iMarker, config->GetMarker_Config_DV(Marker_Tag));
                    config->SetMarker_All_Moving(iMarker, config->GetMarker_Config_Moving(Marker_Tag));
                    config->SetMarker_All_PerBound(iMarker, config->GetMarker_Config_PerBound(Marker_Tag));
                    config->SetMarker_All_Sliding(iMarker, config->GetMarker_Config_Sliding(Marker_Tag));
                    config->SetMarker_All_SendRecv(iMarker, NONE);
                    
                }
                
                /*--- Send-Receive boundaries definition ---*/
                else {
                    unsigned long nelem_vertex = 0, vnodes_vertex;
                    unsigned short transform, matching_zone;
                    getline (mesh_file,text_line);
                    text_line.erase (0,13); nElem_Bound[iMarker] = atoi(text_line.c_str());
                    bound[iMarker] = new CPrimalGrid* [nElem_Bound[iMarker]];
                    
                    nelem_vertex = 0; ielem = 0;
                    getline (mesh_file,text_line); text_line.erase (0,8);
                    config->SetMarker_All_Boundary(iMarker, SEND_RECEIVE);
                    config->SetMarker_All_SendRecv(iMarker, atoi(text_line.c_str()));
                    
                    for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
                        getline(mesh_file,text_line);
                        istringstream bound_line(text_line);
                        bound_line >> VTK_Type; bound_line >> vnodes_vertex; bound_line >> transform;
                        
                        if (val_nZone > 1) bound_line >> matching_zone;
                        bound[iMarker][ielem] = new CVertexMPI(vnodes_vertex, nDim);
                        bound[iMarker][ielem]->SetRotation_Type(transform);
                        if (val_nZone > 1) bound[iMarker][ielem]->SetMatching_Zone(matching_zone);
                        ielem++; nelem_vertex++;
                        if (config->GetMarker_All_SendRecv(iMarker) < 0)
                            node[vnodes_vertex]->SetDomain(false);
                    }
                    
                }
                
            }
        }
        /*--- Read periodic transformation info (center, rotation, translation) ---*/
        position = text_line.find ("NPERIODIC=",0);
        if (position != string::npos) {
            unsigned short nPeriodic, iPeriodic, iIndex;
            
            /*--- Set bool signifying that periodic transormations were found ---*/
            found_transform = true;
            
            /*--- Read and store the number of transformations. ---*/
            text_line.erase (0,10); nPeriodic = atoi(text_line.c_str());
            if (rank == MASTER_NODE) {
                if (nPeriodic - 1 != 0)
                    cout << nPeriodic - 1 << " periodic transformations." << endl;
            }
            config->SetnPeriodicIndex(nPeriodic);
            
            /*--- Store center, rotation, & translation in that order for each. ---*/
            for (iPeriodic = 0; iPeriodic < nPeriodic; iPeriodic++) {
                getline (mesh_file,text_line);
                position = text_line.find ("PERIODIC_INDEX=",0);
                if (position != string::npos) {
                    text_line.erase (0,15); iIndex = atoi(text_line.c_str());
                    if (iIndex != iPeriodic) {
                        cout << "PERIODIC_INDEX out of order in SU2 file!!" << endl;
                        cout << "Press any key to exit..." << endl;
                        cin.get();
                        exit(1);
                    }
                }
                double* center    = new double[3];
                double* rotation  = new double[3];
                double* translate = new double[3];
                getline (mesh_file,text_line);
                istringstream cent(text_line);
                cent >> center[0]; cent >> center[1]; cent >> center[2];
                config->SetPeriodicCenter(iPeriodic, center);
                getline (mesh_file,text_line);
                istringstream rot(text_line);
                rot >> rotation[0]; rot >> rotation[1]; rot >> rotation[2];
                config->SetPeriodicRotation(iPeriodic, rotation);
                getline (mesh_file,text_line);
                istringstream tran(text_line);
                tran >> translate[0]; tran >> translate[1]; tran >> translate[2];
                config->SetPeriodicTranslate(iPeriodic, translate);
            }
            
        }
        
    }
    
    /*--- If no periodic transormations were found, store default zeros ---*/
    if (!found_transform) {
        unsigned short nPeriodic = 1, iPeriodic = 0;
        config->SetnPeriodicIndex(nPeriodic);
        double* center    = new double[3];
        double* rotation  = new double[3];
        double* translate = new double[3];
        for (unsigned short iDim = 0; iDim < 3; iDim++) {
            center[iDim] = 0.0; rotation[iDim] = 0.0; translate[iDim] = 0.0;
        }
        config->SetPeriodicCenter(iPeriodic, center);
        config->SetPeriodicRotation(iPeriodic, rotation);
        config->SetPeriodicTranslate(iPeriodic, translate);
    }
    
    /*--- Close the input file ---*/
    mesh_file.close();
    
    if (config->GetDivide_Element()) {
        if (Local2Global != NULL) delete [] Local2Global;
    }
    
}

void CPhysicalGeometry::Read_NETCDF_Format(CConfig *config, string val_mesh_filename, unsigned short val_iZone, unsigned short val_nZone) {
    
    /*--- Local variables and initialization ---*/
	string text_line, Marker_Tag;
	ifstream mesh_file;
	unsigned short VTK_Type, iMarker;
	unsigned long ielem = 0,
    vnodes_triangle[3], vnodes_quad[4], vnodes_tetra[4], vnodes_hexa[8], vnodes_wedge[6];
	char cstr[200];
	string::size_type position;
	nZone = val_nZone;
    
    /*--- Initialize counters for local/global points & elements ---*/
    
    FinestMGLevel = true;
	Global_nPoint = 0; Global_nPointDomain = 0; Global_nElem = 0;
    nelem_edge     = 0; Global_nelem_edge     = 0;
    nelem_triangle = 0; Global_nelem_triangle = 0;
    nelem_quad     = 0; Global_nelem_quad     = 0;
    nelem_tetra    = 0; Global_nelem_tetra    = 0;
    nelem_hexa     = 0; Global_nelem_hexa     = 0;
    nelem_wedge    = 0; Global_nelem_wedge    = 0;
    nelem_pyramid  = 0; Global_nelem_pyramid  = 0;
    
    /*--- Throw error if not in serial mode. ---*/
    
    unsigned short Marker_Index, marker, icommas, iDim;
    unsigned long nmarker, ielem_triangle, ielem_hexa, ncoord, iSurfElem, ielem_wedge,
    ielem_quad, *marker_list, **surf_elem, ielem_surface, *surf_marker, nSurfElem;
    double coord;
    string::size_type position_;
    bool stop, add;
    
    ielem_surface = 0; nSurfElem = 0;
    surf_marker = NULL;
    surf_elem = NULL;
    
    
    nDim = 3; cout << "Three dimensional problem." << endl;
    
    /*--- Open grid file ---*/
    strcpy (cstr, val_mesh_filename.c_str());
    mesh_file.open(cstr, ios::in);
    if (mesh_file.fail()) {
        cout << "There is no geometry file (CPhysicalGeometry)!!" << endl;
        cout << "Press any key to exit..." << endl;
        cin.get();
        exit(1);
        
    }
    
    while (getline (mesh_file, text_line)) {
        
        position = text_line.find ("no_of_elements = ",0);
        if (position != string::npos) {
            text_line.erase (0,17); nElem = atoi(text_line.c_str());
            cout << nElem << " inner elements to store." << endl;
            elem = new CPrimalGrid*[nElem]; }
        
        position = text_line.find ("no_of_surfaceelements = ",0);
        if (position != string::npos) {
            text_line.erase (0,24); nSurfElem = atoi(text_line.c_str());
            cout << nSurfElem << " surface elements to store." << endl;
            surf_elem = new unsigned long* [nSurfElem];
            for (ielem_surface = 0; ielem_surface < nSurfElem; ielem_surface++)
                surf_elem[ielem_surface] = new unsigned long [5];
            ielem_surface = 0;
        }
        
        position = text_line.find ("no_of_points = ",0);
        if (position != string::npos) {
            text_line.erase (0,15); nPoint = atoi(text_line.c_str()); nPointDomain = nPoint;
            cout << nPoint << " points to store." << endl;
            node = new CPoint*[nPoint]; }
        
        position = text_line.find ("no_of_tetraeders = ",0);
        if (position != string::npos) {
            text_line.erase (0,19); nelem_tetra = atoi(text_line.c_str());
            cout << nelem_tetra << " tetraeders elements to store." << endl; }
        
        position = text_line.find ("no_of_prisms = ",0);
        if (position != string::npos) {
            text_line.erase (0,15); nelem_wedge = atoi(text_line.c_str());
            cout << nelem_wedge << " prims elements to store." << endl; }
        
        position = text_line.find ("no_of_hexaeders = ",0);
        if (position != string::npos) {
            text_line.erase (0,18); nelem_hexa = atoi(text_line.c_str());
            cout << nelem_hexa << " hexaeders elements to store." << endl; }
        
        position = text_line.find ("no_of_surfacetriangles = ",0);
        if (position != string::npos) {
            text_line.erase (0,25); nelem_triangle = atoi(text_line.c_str());
            cout << nelem_triangle << " surface triangle elements to store." << endl; }
        
        position = text_line.find ("no_of_surfacequadrilaterals = ",0);
        if (position != string::npos) {
            text_line.erase (0,30); nelem_quad = atoi(text_line.c_str());
            cout << nelem_quad << " surface quadrilaterals elements to store." << endl; }
        
        position = text_line.find ("points_of_tetraeders =",0);
        if (position != string::npos) {
            for (unsigned long ielem_tetra = 0; ielem_tetra < nelem_tetra; ielem_tetra++) {
                getline(mesh_file,text_line);
                for (unsigned short icommas = 0; icommas < 15; icommas++) {
                    position_ = text_line.find( ",", 0 ); if(position_!=string::npos) text_line.erase (position_,1);
                    position_ = text_line.find( ";", 0 ); if(position_!=string::npos) text_line.erase (position_,1);
                }
                istringstream elem_line(text_line);
                VTK_Type = TETRAHEDRON;
                elem_line >> vnodes_tetra[0]; elem_line >> vnodes_tetra[1]; elem_line >> vnodes_tetra[2]; elem_line >> vnodes_tetra[3];
                elem[ielem] = new CTetrahedron(vnodes_tetra[0],vnodes_tetra[1],vnodes_tetra[2],vnodes_tetra[3]);
                ielem++;
            }
            cout << "finish tetrahedron element reading" << endl;
        }
        
        position = text_line.find ("points_of_prisms =",0);
        if (position != string::npos) {
            for (ielem_wedge = 0; ielem_wedge < nelem_wedge; ielem_wedge++) {
                getline(mesh_file,text_line);
                for (icommas = 0; icommas < 15; icommas++) {
                    position_ = text_line.find( ",", 0 ); if(position_!=string::npos) text_line.erase (position_,1);
                    position_ = text_line.find( ";", 0 ); if(position_!=string::npos) text_line.erase (position_,1);
                }
                istringstream elem_line(text_line);
                VTK_Type = WEDGE;
                elem_line >> vnodes_wedge[0]; elem_line >> vnodes_wedge[1]; elem_line >> vnodes_wedge[2]; elem_line >> vnodes_wedge[3]; elem_line >> vnodes_wedge[4]; elem_line >> vnodes_wedge[5];
                elem[ielem] = new CWedge(vnodes_wedge[0],vnodes_wedge[1],vnodes_wedge[2],vnodes_wedge[3],vnodes_wedge[4],vnodes_wedge[5]);
                ielem++;
            }
            cout << "finish prims element reading" << endl;
        }
        
        position = text_line.find ("points_of_hexaeders =",0);
        if (position != string::npos) {
            for (ielem_hexa = 0; ielem_hexa < nelem_hexa; ielem_hexa++) {
                getline(mesh_file,text_line);
                for (icommas = 0; icommas < 15; icommas++) {
                    position_ = text_line.find( ",", 0 ); if(position_!=string::npos) text_line.erase (position_,1);
                    position_ = text_line.find( ";", 0 ); if(position_!=string::npos) text_line.erase (position_,1);
                }
                istringstream elem_line(text_line);
                VTK_Type = HEXAHEDRON;
                elem_line >> vnodes_hexa[0]; elem_line >> vnodes_hexa[1]; elem_line >> vnodes_hexa[2]; elem_line >> vnodes_hexa[3];
                elem_line >> vnodes_hexa[4]; elem_line >> vnodes_hexa[5]; elem_line >> vnodes_hexa[6]; elem_line >> vnodes_hexa[7];
                elem[ielem] = new CHexahedron(vnodes_hexa[0],vnodes_hexa[1],vnodes_hexa[2],vnodes_hexa[3],vnodes_hexa[4],vnodes_hexa[5],vnodes_hexa[6],vnodes_hexa[7]);
                ielem++;
            }
            cout << "finish hexaeders element reading" << endl;
        }
        
        position = text_line.find ("points_of_surfacetriangles =",0);
        if (position != string::npos) {
            for (ielem_triangle = 0; ielem_triangle < nelem_triangle; ielem_triangle++) {
                getline(mesh_file,text_line);
                for (icommas = 0; icommas < 15; icommas++) {
                    position_ = text_line.find( ",", 0 ); if(position_!=string::npos) text_line.erase (position_,1);
                    position_ = text_line.find( ";", 0 ); if(position_!=string::npos) text_line.erase (position_,1);
                }
                istringstream elem_line(text_line);
                elem_line >> vnodes_triangle[0]; elem_line >> vnodes_triangle[1]; elem_line >> vnodes_triangle[2];
                surf_elem[ielem_surface][0]= 3;
                surf_elem[ielem_surface][1]= vnodes_triangle[0];
                surf_elem[ielem_surface][2]= vnodes_triangle[1];
                surf_elem[ielem_surface][3]= vnodes_triangle[2];
                ielem_surface++;
            }
            cout << "finish surface triangles element reading" << endl;
        }
        
        position = text_line.find ("points_of_surfacequadrilaterals =",0);
        if (position != string::npos) {
            for (ielem_quad = 0; ielem_quad < nelem_quad; ielem_quad++) {
                getline(mesh_file,text_line);
                for (icommas = 0; icommas < 15; icommas++) {
                    position_ = text_line.find( ",", 0 ); if(position_!=string::npos) text_line.erase (position_,1);
                    position_ = text_line.find( ";", 0 ); if(position_!=string::npos) text_line.erase (position_,1);
                }
                istringstream elem_line(text_line);
                elem_line >> vnodes_quad[0]; elem_line >> vnodes_quad[1]; elem_line >> vnodes_quad[2]; elem_line >> vnodes_quad[3];
                surf_elem[ielem_surface][0]= 4;
                surf_elem[ielem_surface][1]= vnodes_quad[0];
                surf_elem[ielem_surface][2]= vnodes_quad[1];
                surf_elem[ielem_surface][3]= vnodes_quad[2];
                surf_elem[ielem_surface][4]= vnodes_quad[3];
                ielem_surface++;
            }
            cout << "finish surface quadrilaterals element reading" << endl;
        }
        
        position = text_line.find ("boundarymarker_of_surfaces =",0);
        if (position != string::npos) {
            nmarker=0;
            stop = false;
            surf_marker = new unsigned long [nelem_triangle + nelem_quad];
            
            text_line.erase (0,29);
            for (icommas = 0; icommas < 50; icommas++) {
                position_ = text_line.find( ",", 0 );
                if(position_!=string::npos) text_line.erase (position_,1);
            }
            
            stringstream  point_line(text_line);
            while (point_line >> marker,!point_line.eof()) {
                surf_marker[nmarker] = marker;
                nmarker++; }
            
            while (!stop) {
                getline(mesh_file,text_line);
                for (icommas = 0; icommas < 50; icommas++) {
                    position_ = text_line.find( ",", 0 );
                    if(position_!=string::npos) text_line.erase (position_,1);
                    position_ = text_line.find( ";", 0 );
                    if(position_!=string::npos) text_line.erase (position_,1);
                }
                stringstream  point_line(text_line);
                while (point_line>> marker,!point_line.eof()) {
                    surf_marker[nmarker] = marker;
                    if (nmarker == nSurfElem-1) {stop = true; break;}
                    nmarker++;
                }
            }
        }
        
        for (iDim = 0; iDim < nDim; iDim++) {
            ncoord = 0; stop = false;
            if (iDim == 0) position = text_line.find ("points_xc = ",0);
            if (iDim == 1) position = text_line.find ("points_yc = ",0);
            if (iDim == 2) position = text_line.find ("points_zc = ",0);
            
            if (position != string::npos) {
                text_line.erase (0,12);
                for (icommas = 0; icommas < 50; icommas++) {
                    position_ = text_line.find( ",", 0 );
                    if(position_!=string::npos) text_line.erase (position_,1);
                }
                stringstream  point_line(text_line);
                while (point_line>> coord,!point_line.eof()) {
                    if (iDim==0) node[ncoord] = new CPoint(coord, 0.0, 0.0, ncoord, config);
                    if (iDim==1) node[ncoord]->SetCoord(1, coord);
                    if (iDim==2) node[ncoord]->SetCoord(2, coord);
                    ncoord++; }
                while (!stop) {
                    getline(mesh_file,text_line);
                    for (icommas = 0; icommas < 50; icommas++) {
                        position_ = text_line.find( ",", 0 );
                        if(position_!=string::npos) text_line.erase (position_,1);
                        position_ = text_line.find( ";", 0 );
                        if(position_!=string::npos) text_line.erase (position_,1);
                    }
                    stringstream  point_line(text_line);
                    while (point_line>> coord,!point_line.eof()) {
                        if (iDim==0) node[ncoord] = new CPoint(coord, 0.0, 0.0, ncoord, config);
                        if (iDim==1) node[ncoord]->SetCoord(1, coord);
                        if (iDim==2) node[ncoord]->SetCoord(2, coord);
                        if (ncoord == nPoint-1) {stop = true; break;}
                        ncoord++;
                    }
                }
                if (iDim==0) cout << "finish point xc reading" << endl;
                if (iDim==1) cout << "finish point yc reading" << endl;
                if (iDim==2) cout << "finish point zc reading" << endl;
            }
        }
    }
    
    
    /*--- Create a list with all the markers ---*/
    marker_list = new unsigned long [MAX_NUMBER_MARKER];
    marker_list[0] = surf_marker[0]; nMarker = 1;
    for (iSurfElem = 0; iSurfElem < nSurfElem; iSurfElem++) {
        add = true;
        for (iMarker = 0; iMarker < nMarker; iMarker++)
            if (marker_list[iMarker] == surf_marker[iSurfElem]) {
                add = false; break; }
        if (add) {
            marker_list[nMarker] = surf_marker[iSurfElem];
            nMarker++;
        }
    }
    
    nElem_Bound = new unsigned long [nMarker];
    
    /*--- Compute the number of element per marker ---*/
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
        nElem_Bound[iMarker] = 0;
        for (iSurfElem = 0; iSurfElem < nSurfElem; iSurfElem++)
            if (surf_marker[iSurfElem] == marker_list[iMarker])
                nElem_Bound[iMarker]++;
    }
    
    
    /*--- Realate the marker index with the position in the array of markers ---*/
    unsigned short *Index_to_Marker;
    Index_to_Marker = new unsigned short [MAX_INDEX_VALUE];
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
        Marker_Index = marker_list[iMarker];
        Index_to_Marker[Marker_Index] = iMarker;
    }
    
    bound = new CPrimalGrid**[nMarker];
    
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
        Marker_Index = marker_list[iMarker];
        bound[iMarker] = new CPrimalGrid* [nElem_Bound[iMarker]];
        ielem_triangle = 0; ielem_quad = 0;
        for (iSurfElem = 0; iSurfElem < nSurfElem; iSurfElem++)
            if (surf_marker[iSurfElem] == Marker_Index) {
                if (surf_elem[iSurfElem][0] == 3) {
                    vnodes_triangle[0] = surf_elem[iSurfElem][1]; vnodes_triangle[1] = surf_elem[iSurfElem][2]; vnodes_triangle[2] = surf_elem[iSurfElem][3];
                    bound[iMarker][ielem_triangle+ielem_quad] = new CTriangle(vnodes_triangle[0],vnodes_triangle[1],vnodes_triangle[2],3);
                    ielem_triangle ++;
                }
                if (surf_elem[iSurfElem][0] == 4) {
                    vnodes_quad[0] = surf_elem[iSurfElem][1]; vnodes_quad[1] = surf_elem[iSurfElem][2]; vnodes_quad[2] = surf_elem[iSurfElem][3]; vnodes_quad[3] = surf_elem[iSurfElem][4];
                    bound[iMarker][ielem_triangle+ielem_quad] = new CRectangle(vnodes_quad[0],vnodes_quad[1],vnodes_quad[2],vnodes_quad[3],3);
                    ielem_quad ++;
                }
            }
    }
    
    
    
    /*--- Update config information storing the boundary information in the right place ---*/
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
        
        stringstream out;
        out << marker_list[iMarker];
        Marker_Tag = out.str();
        
        Tag_to_Marker = new string [MAX_INDEX_VALUE];
        Tag_to_Marker[config->GetMarker_Config_Tag(Marker_Tag)] = Marker_Tag;
        config->SetMarker_All_Tag(iMarker, Marker_Tag);
        config->SetMarker_All_Boundary(iMarker, config->GetMarker_Config_Boundary(Marker_Tag));
        config->SetMarker_All_Monitoring(iMarker, config->GetMarker_Config_Monitoring(Marker_Tag));
        config->SetMarker_All_Designing(iMarker, config->GetMarker_Config_Designing(Marker_Tag));
        config->SetMarker_All_Plotting(iMarker, config->GetMarker_Config_Plotting(Marker_Tag));
        config->SetMarker_All_DV(iMarker, config->GetMarker_Config_DV(Marker_Tag));
        config->SetMarker_All_Moving(iMarker, config->GetMarker_Config_Moving(Marker_Tag));
        config->SetMarker_All_SendRecv(iMarker, NONE);
    }
    
}

void CPhysicalGeometry::Check_Orientation(CConfig *config) {
	unsigned long Point_1, Point_2, Point_3, Point_4, Point_5, Point_6,
	iElem, Point_1_Surface, Point_2_Surface, Point_3_Surface, Point_4_Surface,
	iElem_Domain, Point_Domain = 0, Point_Surface, iElem_Surface;
	double test_1, test_2, test_3, test_4, *Coord_1, *Coord_2, *Coord_3, *Coord_4,
	*Coord_5, *Coord_6, a[3], b[3], c[3], n[3], test;
	unsigned short iDim, iMarker, iNode_Domain, iNode_Surface;
	bool find;
    
	/*--- Loop over all the elements ---*/
	for (iElem = 0; iElem < nElem; iElem++) {
        
		/*--- 2D grid, triangle case ---*/
		if (elem[iElem]->GetVTK_Type() == TRIANGLE) {
            
			Point_1 = elem[iElem]->GetNode(0); Coord_1 = node[Point_1]->GetCoord();
			Point_2 = elem[iElem]->GetNode(1); Coord_2 = node[Point_2]->GetCoord();
			Point_3 = elem[iElem]->GetNode(2); Coord_3 = node[Point_3]->GetCoord();
            
			for(iDim = 0; iDim < nDim; iDim++) {
				a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
				b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]); }
			test = a[0]*b[1]-b[0]*a[1];
            
			if (test < 0.0) elem[iElem]->Change_Orientation();
		}
        
		/*--- 2D grid, rectangle case ---*/
		if (elem[iElem]->GetVTK_Type() == RECTANGLE) {
            
			Point_1 = elem[iElem]->GetNode(0); Coord_1 = node[Point_1]->GetCoord();
			Point_2 = elem[iElem]->GetNode(1); Coord_2 = node[Point_2]->GetCoord();
			Point_3 = elem[iElem]->GetNode(2); Coord_3 = node[Point_3]->GetCoord();
			Point_4 = elem[iElem]->GetNode(3); Coord_4 = node[Point_4]->GetCoord();
            
			for(iDim = 0; iDim < nDim; iDim++) {
				a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
				b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]); }
			test_1 = a[0]*b[1]-b[0]*a[1];
            
			for(iDim = 0; iDim < nDim; iDim++) {
				a[iDim] = 0.5*(Coord_3[iDim]-Coord_2[iDim]);
				b[iDim] = 0.5*(Coord_4[iDim]-Coord_2[iDim]); }
			test_2 = a[0]*b[1]-b[0]*a[1];
            
			for(iDim = 0; iDim < nDim; iDim++) {
				a[iDim] = 0.5*(Coord_4[iDim]-Coord_3[iDim]);
				b[iDim] = 0.5*(Coord_1[iDim]-Coord_3[iDim]); }
			test_3 = a[0]*b[1]-b[0]*a[1];
            
			for(iDim = 0; iDim < nDim; iDim++) {
				a[iDim] = 0.5*(Coord_1[iDim]-Coord_4[iDim]);
				b[iDim] = 0.5*(Coord_3[iDim]-Coord_4[iDim]); }
			test_4 = a[0]*b[1]-b[0]*a[1];
            
			if ((test_1 < 0.0) && (test_2 < 0.0) && (test_3 < 0.0) && (test_4 < 0.0))
				elem[iElem]->Change_Orientation();
		}
        
		/*--- 3D grid, tetrahedron case ---*/
		if (elem[iElem]->GetVTK_Type() == TETRAHEDRON) {
            
			Point_1 = elem[iElem]->GetNode(0); Coord_1 = node[Point_1]->GetCoord();
			Point_2 = elem[iElem]->GetNode(1); Coord_2 = node[Point_2]->GetCoord();
			Point_3 = elem[iElem]->GetNode(2); Coord_3 = node[Point_3]->GetCoord();
			Point_4 = elem[iElem]->GetNode(3); Coord_4 = node[Point_4]->GetCoord();
            
			for(iDim = 0; iDim < nDim; iDim++) {
				a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
				b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]);
				c[iDim] = Coord_4[iDim]-Coord_1[iDim]; }
			n[0] = a[1]*b[2]-b[1]*a[2];
			n[1] = -(a[0]*b[2]-b[0]*a[2]);
			n[2] = a[0]*b[1]-b[0]*a[1];
            
			test = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
			if (test < 0.0) elem[iElem]->Change_Orientation();
            
		}
        
		/*--- 3D grid, wedge case ---*/
		if (elem[iElem]->GetVTK_Type() == WEDGE) {
            
			Point_1 = elem[iElem]->GetNode(0); Coord_1 = node[Point_1]->GetCoord();
			Point_2 = elem[iElem]->GetNode(1); Coord_2 = node[Point_2]->GetCoord();
			Point_3 = elem[iElem]->GetNode(2); Coord_3 = node[Point_3]->GetCoord();
			Point_4 = elem[iElem]->GetNode(3); Coord_4 = node[Point_4]->GetCoord();
			Point_5 = elem[iElem]->GetNode(4); Coord_5 = node[Point_5]->GetCoord();
			Point_6 = elem[iElem]->GetNode(5); Coord_6 = node[Point_6]->GetCoord();
            
			for(iDim = 0; iDim < nDim; iDim++) {
				a[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]);
				b[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
				c[iDim] = (Coord_4[iDim]-Coord_1[iDim])+
                (Coord_5[iDim]-Coord_2[iDim])+
                (Coord_6[iDim]-Coord_3[iDim]); }
            
			/*--- The normal vector should point to the interior of the element ---*/
			n[0] = a[1]*b[2]-b[1]*a[2];
			n[1] = -(a[0]*b[2]-b[0]*a[2]);
			n[2] = a[0]*b[1]-b[0]*a[1];
            
			test_1 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
            
			for(iDim = 0; iDim < nDim; iDim++) {
				a[iDim] = 0.5*(Coord_5[iDim]-Coord_4[iDim]);
				b[iDim] = 0.5*(Coord_6[iDim]-Coord_4[iDim]);
				c[iDim] = (Coord_1[iDim]-Coord_4[iDim])+
                (Coord_2[iDim]-Coord_5[iDim])+
                (Coord_3[iDim]-Coord_6[iDim]); }
            
			/*--- The normal vector should point to the interior of the element ---*/
			n[0] = a[1]*b[2]-b[1]*a[2];
			n[1] = -(a[0]*b[2]-b[0]*a[2]);
			n[2] = a[0]*b[1]-b[0]*a[1];
            
			test_2 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
            
			if ((test_1 < 0.0) || (test_2 < 0.0))
				elem[iElem]->Change_Orientation();
            
		}
        
		if (elem[iElem]->GetVTK_Type() == HEXAHEDRON) {
            
			Point_1 = elem[iElem]->GetNode(0); Coord_1 = node[Point_1]->GetCoord();
			Point_2 = elem[iElem]->GetNode(1); Coord_2 = node[Point_2]->GetCoord();
			Point_3 = elem[iElem]->GetNode(2); Coord_3 = node[Point_3]->GetCoord();
			Point_4 = elem[iElem]->GetNode(5); Coord_4 = node[Point_4]->GetCoord();
            
			for(iDim = 0; iDim < nDim; iDim++) {
				a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
				b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]);
				c[iDim] = Coord_4[iDim]-Coord_1[iDim]; }
			n[0] = a[1]*b[2]-b[1]*a[2];
			n[1] = -(a[0]*b[2]-b[0]*a[2]);
			n[2] = a[0]*b[1]-b[0]*a[1];
            
			test_1 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
            
			Point_1 = elem[iElem]->GetNode(2); Coord_1 = node[Point_1]->GetCoord();
			Point_2 = elem[iElem]->GetNode(3); Coord_2 = node[Point_2]->GetCoord();
			Point_3 = elem[iElem]->GetNode(0); Coord_3 = node[Point_3]->GetCoord();
			Point_4 = elem[iElem]->GetNode(7); Coord_4 = node[Point_4]->GetCoord();
            
			for(iDim = 0; iDim < nDim; iDim++) {
				a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
				b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]);
				c[iDim] = Coord_4[iDim]-Coord_1[iDim]; }
			n[0] = a[1]*b[2]-b[1]*a[2];
			n[1] = -(a[0]*b[2]-b[0]*a[2]);
			n[2] = a[0]*b[1]-b[0]*a[1];
            
			test_2 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
            
			Point_1 = elem[iElem]->GetNode(1); Coord_1 = node[Point_1]->GetCoord();
			Point_2 = elem[iElem]->GetNode(2); Coord_2 = node[Point_2]->GetCoord();
			Point_3 = elem[iElem]->GetNode(3); Coord_3 = node[Point_3]->GetCoord();
			Point_4 = elem[iElem]->GetNode(6); Coord_4 = node[Point_4]->GetCoord();
            
			for(iDim = 0; iDim < nDim; iDim++) {
				a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
				b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]);
				c[iDim] = Coord_4[iDim]-Coord_1[iDim]; }
			n[0] = a[1]*b[2]-b[1]*a[2];
			n[1] = -(a[0]*b[2]-b[0]*a[2]);
			n[2] = a[0]*b[1]-b[0]*a[1];
            
			test_3 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
            
			Point_1 = elem[iElem]->GetNode(3); Coord_1 = node[Point_1]->GetCoord();
			Point_2 = elem[iElem]->GetNode(0); Coord_2 = node[Point_2]->GetCoord();
			Point_3 = elem[iElem]->GetNode(1); Coord_3 = node[Point_3]->GetCoord();
			Point_4 = elem[iElem]->GetNode(4); Coord_4 = node[Point_4]->GetCoord();
            
			for(iDim = 0; iDim < nDim; iDim++) {
				a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
				b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]);
				c[iDim] = Coord_4[iDim]-Coord_1[iDim]; }
			n[0] = a[1]*b[2]-b[1]*a[2];
			n[1] = -(a[0]*b[2]-b[0]*a[2]);
			n[2] = a[0]*b[1]-b[0]*a[1];
            
			test_4 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
            
			if ((test_1 < 0.0) || (test_2 < 0.0) || (test_3 < 0.0)
                || (test_4 < 0.0)) elem[iElem]->Change_Orientation();
            
		}
        
		if (elem[iElem]->GetVTK_Type() == PYRAMID) {
            
			Point_1 = elem[iElem]->GetNode(0); Coord_1 = node[Point_1]->GetCoord();
			Point_2 = elem[iElem]->GetNode(1); Coord_2 = node[Point_2]->GetCoord();
			Point_3 = elem[iElem]->GetNode(2); Coord_3 = node[Point_3]->GetCoord();
			Point_4 = elem[iElem]->GetNode(4); Coord_4 = node[Point_4]->GetCoord();
            
			for(iDim = 0; iDim < nDim; iDim++) {
				a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
				b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]);
				c[iDim] = Coord_4[iDim]-Coord_1[iDim]; }
			n[0] = a[1]*b[2]-b[1]*a[2];
			n[1] = -(a[0]*b[2]-b[0]*a[2]);
			n[2] = a[0]*b[1]-b[0]*a[1];
            
			test_1 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
            
			Point_1 = elem[iElem]->GetNode(2); Coord_1 = node[Point_1]->GetCoord();
			Point_2 = elem[iElem]->GetNode(3); Coord_2 = node[Point_2]->GetCoord();
			Point_3 = elem[iElem]->GetNode(0); Coord_3 = node[Point_3]->GetCoord();
			Point_4 = elem[iElem]->GetNode(4); Coord_4 = node[Point_4]->GetCoord();
            
			for(iDim = 0; iDim < nDim; iDim++) {
				a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
				b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]);
				c[iDim] = Coord_4[iDim]-Coord_1[iDim]; }
			n[0] = a[1]*b[2]-b[1]*a[2];
			n[1] = -(a[0]*b[2]-b[0]*a[2]);
			n[2] = a[0]*b[1]-b[0]*a[1];
            
			test_2 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
            
			if ((test_1 < 0.0) || (test_2 < 0.0))
				elem[iElem]->Change_Orientation();
            
		}
        
	}
    
	/*--- Surface elements ---*/
	for (iMarker = 0; iMarker < nMarker; iMarker++)
		for (iElem_Surface = 0; iElem_Surface < nElem_Bound[iMarker]; iElem_Surface++) {
            
			iElem_Domain = bound[iMarker][iElem_Surface]->GetDomainElement();
			for (iNode_Domain = 0; iNode_Domain < elem[iElem_Domain]->GetnNodes(); iNode_Domain++) {
				Point_Domain = elem[iElem_Domain]->GetNode(iNode_Domain);
				find = false;
				for (iNode_Surface = 0; iNode_Surface < bound[iMarker][iElem_Surface]->GetnNodes(); iNode_Surface++) {
					Point_Surface = bound[iMarker][iElem_Surface]->GetNode(iNode_Surface);
					if (Point_Surface == Point_Domain) {find = true; break;}
				}
				if (!find) break;
			}
            
			/*--- 2D grid, line case ---*/
			if (bound[iMarker][iElem_Surface]->GetVTK_Type() == LINE) {
                
				Point_1_Surface = bound[iMarker][iElem_Surface]->GetNode(0); Coord_1 = node[Point_1_Surface]->GetCoord();
				Point_2_Surface = bound[iMarker][iElem_Surface]->GetNode(1); Coord_2 = node[Point_2_Surface]->GetCoord();
				Coord_3 = node[Point_Domain]->GetCoord();
                
				for(iDim = 0; iDim < nDim; iDim++) {
					a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
					b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]); }
				test = a[0]*b[1]-b[0]*a[1];
                
				if (test < 0.0) bound[iMarker][iElem_Surface]->Change_Orientation();
			}
            
			/*--- 3D grid, triangle case ---*/
			if (bound[iMarker][iElem_Surface]->GetVTK_Type() == TRIANGLE) {
                
				Point_1_Surface = bound[iMarker][iElem_Surface]->GetNode(0); Coord_1 = node[Point_1_Surface]->GetCoord();
				Point_2_Surface = bound[iMarker][iElem_Surface]->GetNode(1); Coord_2 = node[Point_2_Surface]->GetCoord();
				Point_3_Surface = bound[iMarker][iElem_Surface]->GetNode(2); Coord_3 = node[Point_3_Surface]->GetCoord();
				Coord_4 = node[Point_Domain]->GetCoord();
                
				for(iDim = 0; iDim < nDim; iDim++) {
					a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
					b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]);
					c[iDim] = Coord_4[iDim]-Coord_1[iDim]; }
				n[0] = a[1]*b[2]-b[1]*a[2];
				n[1] = -(a[0]*b[2]-b[0]*a[2]);
				n[2] = a[0]*b[1]-b[0]*a[1];
                
				test = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
				if (test < 0.0) bound[iMarker][iElem_Surface]->Change_Orientation();
			}
            
			if (bound[iMarker][iElem_Surface]->GetVTK_Type() == RECTANGLE) {
                
				Point_1_Surface = bound[iMarker][iElem_Surface]->GetNode(0); Coord_1 = node[Point_1_Surface]->GetCoord();
				Point_2_Surface = bound[iMarker][iElem_Surface]->GetNode(1); Coord_2 = node[Point_2_Surface]->GetCoord();
				Point_3_Surface = bound[iMarker][iElem_Surface]->GetNode(2); Coord_3 = node[Point_3_Surface]->GetCoord();
				Point_4_Surface = bound[iMarker][iElem_Surface]->GetNode(3); Coord_4 = node[Point_4_Surface]->GetCoord();
				Coord_5 = node[Point_Domain]->GetCoord();
                
				for(iDim = 0; iDim < nDim; iDim++) {
					a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
					b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]);
					c[iDim] = Coord_5[iDim]-Coord_1[iDim]; }
				n[0] = a[1]*b[2]-b[1]*a[2];
				n[1] = -(a[0]*b[2]-b[0]*a[2]);
				n[2] = a[0]*b[1]-b[0]*a[1];
				test_1 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
                
				for(iDim = 0; iDim < nDim; iDim++) {
					a[iDim] = 0.5*(Coord_3[iDim]-Coord_2[iDim]);
					b[iDim] = 0.5*(Coord_4[iDim]-Coord_2[iDim]);
					c[iDim] = Coord_5[iDim]-Coord_2[iDim]; }
				n[0] = a[1]*b[2]-b[1]*a[2];
				n[1] = -(a[0]*b[2]-b[0]*a[2]);
				n[2] = a[0]*b[1]-b[0]*a[1];
				test_2 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
                
				for(iDim = 0; iDim < nDim; iDim++) {
					a[iDim] = 0.5*(Coord_4[iDim]-Coord_3[iDim]);
					b[iDim] = 0.5*(Coord_1[iDim]-Coord_3[iDim]);
					c[iDim] = Coord_5[iDim]-Coord_3[iDim]; }
				n[0] = a[1]*b[2]-b[1]*a[2];
				n[1] = -(a[0]*b[2]-b[0]*a[2]);
				n[2] = a[0]*b[1]-b[0]*a[1];
				test_3 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
                
				for(iDim = 0; iDim < nDim; iDim++) {
					a[iDim] = 0.5*(Coord_1[iDim]-Coord_4[iDim]);
					b[iDim] = 0.5*(Coord_3[iDim]-Coord_4[iDim]);
					c[iDim] = Coord_5[iDim]-Coord_4[iDim]; }
				n[0] = a[1]*b[2]-b[1]*a[2];
				n[1] = -(a[0]*b[2]-b[0]*a[2]);
				n[2] = a[0]*b[1]-b[0]*a[1];
				test_4 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
                
				if ((test_1 < 0.0) && (test_2 < 0.0) && (test_3 < 0.0) && (test_4 < 0.0))
					bound[iMarker][iElem_Surface]->Change_Orientation();
			}
		}
}

void CPhysicalGeometry::SetEsuP(void) {
	unsigned long iPoint, iElem;
	unsigned short iNode;
    
	/*--- Loop over all the elements ---*/
	for(iElem = 0; iElem < nElem; iElem++)
    /*--- Loop over all the nodes of an element ---*/
		for(iNode = 0; iNode < elem[iElem]->GetnNodes(); iNode++) {
			iPoint = elem[iElem]->GetNode(iNode);
			/*--- Store the element into the point ---*/
			node[iPoint]->SetElem(iElem);
		}
}

void CPhysicalGeometry::ComputeWall_Distance(CConfig *config) {
    
	double *coord, dist2, dist;
	unsigned short iDim, iMarker;
	unsigned long iPoint, iVertex, nVertex_SolidWall;
    
    int rank = MASTER_NODE;
	if (rank == MASTER_NODE)
		cout << "Computing wall distances." << endl;
    
	/*--- Compute the total number of nodes on no-slip boundaries ---*/
    
	nVertex_SolidWall = 0;
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if ((config->GetMarker_All_Boundary(iMarker) == HEAT_FLUX) ||
            (config->GetMarker_All_Boundary(iMarker) == ISOTHERMAL))
			nVertex_SolidWall += GetnVertex(iMarker);
    
	/*--- Allocate an array to hold boundary node coordinates ---*/
    
	double **Coord_bound;
	Coord_bound = new double* [nVertex_SolidWall];
	for (iVertex = 0; iVertex < nVertex_SolidWall; iVertex++)
		Coord_bound[iVertex] = new double [nDim];
    
	/*--- Retrieve and store the coordinates of the no-slip boundary nodes ---*/
    
	nVertex_SolidWall = 0;
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
		if ((config->GetMarker_All_Boundary(iMarker) == HEAT_FLUX) ||
            (config->GetMarker_All_Boundary(iMarker) == ISOTHERMAL))
			for (iVertex = 0; iVertex < GetnVertex(iMarker); iVertex++) {
				iPoint = vertex[iMarker][iVertex]->GetNode();
				for (iDim = 0; iDim < nDim; iDim++)
					Coord_bound[nVertex_SolidWall][iDim] = node[iPoint]->GetCoord(iDim);
				nVertex_SolidWall++;
			}
    }
    
	/*--- Loop over all interior mesh nodes and compute the distances to each
     of the no-slip boundary nodes. Store the minimum distance to the wall for
     each interior mesh node. ---*/
    
	for (iPoint = 0; iPoint < GetnPoint(); iPoint++) {
		coord = node[iPoint]->GetCoord();
		dist = 1E20;
		for (iVertex = 0; iVertex < nVertex_SolidWall; iVertex++) {
			dist2 = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				dist2 += (coord[iDim]-Coord_bound[iVertex][iDim])
                *(coord[iDim]-Coord_bound[iVertex][iDim]);
			if (dist2 < dist) dist = dist2;
		}
		node[iPoint]->SetWall_Distance(sqrt(dist));
	}
    
	/*--- Deallocate the vector of boundary coordinates. ---*/
    
	for (iVertex = 0; iVertex < nVertex_SolidWall; iVertex++)
		delete[] Coord_bound[iVertex];
	delete[] Coord_bound;
    
}

void CPhysicalGeometry::SetPositive_ZArea(CConfig *config) {
	unsigned short iMarker, Boundary, Monitoring;
	unsigned long iVertex, iPoint;
	double *Normal, PositiveZArea;
	int rank = MASTER_NODE;
    
	PositiveZArea = 0.0;
	for (iMarker = 0; iMarker < nMarker; iMarker++) {
		Boundary = config->GetMarker_All_Boundary(iMarker);
		Monitoring = config->GetMarker_All_Monitoring(iMarker);
        
		if (((Boundary == EULER_WALL) || (Boundary == HEAT_FLUX) ||
             (Boundary == ISOTHERMAL)) && (Monitoring == YES))
			for(iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
				iPoint = vertex[iMarker][iVertex]->GetNode();
				if (node[iPoint]->GetDomain()) {
					Normal = vertex[iMarker][iVertex]->GetNormal();
					if (Normal[nDim-1] < 0) PositiveZArea -= Normal[nDim-1];
				}
			}
	}
    
	if (config->GetRefAreaCoeff() == 0.0)
		config->SetRefAreaCoeff(PositiveZArea);
    
	if (rank == MASTER_NODE) {
		if (nDim == 2) cout << "Area projection in the y-plane = "<< PositiveZArea << "." << endl;
		else cout << "Area projection in the z-plane = "<< PositiveZArea << "." << endl;
	}
    
}

void CPhysicalGeometry::SetPsuP(void) {
    
	unsigned short Node_Neighbor, iNode, iNeighbor;
	unsigned long jElem, Point_Neighbor, iPoint, iElem;
    
	/*--- Loop over all the points ---*/
	for(iPoint = 0; iPoint < nPoint; iPoint++)
    /*--- Loop over all elements shared by the point ---*/
		for(iElem = 0; iElem < node[iPoint]->GetnElem(); iElem++) {
			jElem = node[iPoint]->GetElem(iElem);
			/*--- If we find the point iPoint in the surronding element ---*/
			for(iNode = 0; iNode < elem[jElem]->GetnNodes(); iNode++)
				if (elem[jElem]->GetNode(iNode) == iPoint)
                /*--- Localize the local index of the neighbor of iPoint in the element ---*/
					for(iNeighbor = 0; iNeighbor < elem[jElem]->GetnNeighbor_Nodes(iNode); iNeighbor++) {
						Node_Neighbor = elem[jElem]->GetNeighbor_Nodes(iNode,iNeighbor);
						Point_Neighbor = elem[jElem]->GetNode(Node_Neighbor);
						/*--- Store the point into the point ---*/
						node[iPoint]->SetPoint(Point_Neighbor);
					}
		}
    
	/*--- Set the number of neighbors variable, this is
	 important for JST and multigrid in parallel ---*/
	for(iPoint = 0; iPoint < nPoint; iPoint++)
		node[iPoint]->SetnNeighbor(node[iPoint]->GetnPoint());
    
}

void CPhysicalGeometry::SetPsuP_FEA(void) {
    
	unsigned short Node_Neighbor, iNode, iNeighbor;
	unsigned long jElem, Point_Neighbor, iPoint, jPoint, iElem;
    
	/*--- Loop over all the points ---*/
	for(iPoint = 0; iPoint < nPoint; iPoint++) {
        
        /*--- Loop over all elements shared by the point ---*/
		for(iElem = 0; iElem < node[iPoint]->GetnElem(); iElem++) {
			jElem = node[iPoint]->GetElem(iElem);
            
			/*--- If we find the point iPoint in the surrounding element ---*/
			for(iNode = 0; iNode < elem[jElem]->GetnNodes(); iNode++)
				if (elem[jElem]->GetNode(iNode) == iPoint)
                /*--- Localize the local index of the neighbor of iPoint in the element ---*/
					for(iNeighbor = 0; iNeighbor < elem[jElem]->GetnNeighbor_Nodes(iNode); iNeighbor++) {
						Node_Neighbor = elem[jElem]->GetNeighbor_Nodes(iNode,iNeighbor);
						Point_Neighbor = elem[jElem]->GetNode(Node_Neighbor);
						/*--- Store the point into the point ---*/
						node[iPoint]->SetPoint(Point_Neighbor);
					}
		}
    }
    
    /*--- For grid deformation using the linear elasticity equations,
     we will cut each element into either triangles (2-D) or tetrahedra (3-D)
     because we already have these shape functions implemented. We only do
     this internally however, because we want the deformed mesh to retain
     the original element connectivity. Therefore, we add the new edges
     in this routine manually for these divisions so that the global stiffness
     matrix is constructed correctly. ---*/
    
    for(iElem = 0; iElem < nElem; iElem++) {
        
        /*--- Divide quads into 2 triangles ---*/
        
        if (elem[iElem]->GetVTK_Type() == RECTANGLE) {
            
            iPoint = elem[iElem]->GetNode(0);
            jPoint = elem[iElem]->GetNode(2);
            
            node[iPoint]->SetPoint(jPoint);
            node[jPoint]->SetPoint(iPoint);
            
        }
        
        /*--- Divide hexehedra into 5 tetrahedra ---*/
        
        if (elem[iElem]->GetVTK_Type() == HEXAHEDRON) {
            
            /*--- Cut each of the 6 quad faces of the hex cell ---*/
            
            iPoint = elem[iElem]->GetNode(0);
            jPoint = elem[iElem]->GetNode(2);
            node[iPoint]->SetPoint(jPoint);
            node[jPoint]->SetPoint(iPoint);
            
            iPoint = elem[iElem]->GetNode(2);
            jPoint = elem[iElem]->GetNode(7);
            node[iPoint]->SetPoint(jPoint);
            node[jPoint]->SetPoint(iPoint);
            
            iPoint = elem[iElem]->GetNode(5);
            jPoint = elem[iElem]->GetNode(7);
            node[iPoint]->SetPoint(jPoint);
            node[jPoint]->SetPoint(iPoint);
            
            iPoint = elem[iElem]->GetNode(0);
            jPoint = elem[iElem]->GetNode(5);
            node[iPoint]->SetPoint(jPoint);
            node[jPoint]->SetPoint(iPoint);
            
            iPoint = elem[iElem]->GetNode(0);
            jPoint = elem[iElem]->GetNode(7);
            node[iPoint]->SetPoint(jPoint);
            node[jPoint]->SetPoint(iPoint);
            
            iPoint = elem[iElem]->GetNode(2);
            jPoint = elem[iElem]->GetNode(5);
            node[iPoint]->SetPoint(jPoint);
            node[jPoint]->SetPoint(iPoint);
            
        }
        
        /*--- Divide prisms into 3 tetrahedra ---*/
        
        if (elem[iElem]->GetVTK_Type() == WEDGE) {
            
            /*--- Cut each of the 3 quad faces of the prism ---*/
            
            iPoint = elem[iElem]->GetNode(0);
            jPoint = elem[iElem]->GetNode(4);
            node[iPoint]->SetPoint(jPoint);
            node[jPoint]->SetPoint(iPoint);
            
            iPoint = elem[iElem]->GetNode(2);
            jPoint = elem[iElem]->GetNode(4);
            node[iPoint]->SetPoint(jPoint);
            node[jPoint]->SetPoint(iPoint);
            
            iPoint = elem[iElem]->GetNode(0);
            jPoint = elem[iElem]->GetNode(5);
            node[iPoint]->SetPoint(jPoint);
            node[jPoint]->SetPoint(iPoint);
            
        }
        
        /*--- Divide pyramids into 2 tetrahedra ---*/
        
        if (elem[iElem]->GetVTK_Type() == PYRAMID) {
            
            /*--- Cut the single quad face of the pyramid ---*/
            
            iPoint = elem[iElem]->GetNode(0);
            jPoint = elem[iElem]->GetNode(2);
            node[iPoint]->SetPoint(jPoint);
            node[jPoint]->SetPoint(iPoint);
            
        }
        
    }
    
	/*--- Set the number of neighbors variable, this is
	 important for JST and multigrid in parallel ---*/
	for(iPoint = 0; iPoint < nPoint; iPoint++)
		node[iPoint]->SetnNeighbor(node[iPoint]->GetnPoint());
    
}

void CPhysicalGeometry::SetEsuE(void) {
	unsigned short first_elem_face, second_elem_face, iFace, iNode, jElem;
	unsigned long face_point, Test_Elem, iElem;
    
	/*--- Loop over all the elements, faces and nodes ---*/
	for(iElem = 0; iElem < nElem; iElem++)
		for (iFace = 0; iFace < elem[iElem]->GetnFaces(); iFace++)
			for (iNode = 0; iNode < elem[iElem]->GetnNodesFace(iFace); iNode++) {
				face_point = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace,iNode));
				/*--- Loop over all elements sharing the face point ---*/
				for(jElem = 0; jElem < node[face_point]->GetnElem(); jElem++) {
					Test_Elem = node[face_point]->GetElem(jElem);
					/*--- If it is a new element in this face ---*/
					if ((elem[iElem]->GetNeighbor_Elements(iFace) == -1) && (iElem < Test_Elem) &&
						(FindFace(iElem, Test_Elem, first_elem_face, second_elem_face))) {
                        /*--- Localice which faces are sharing both elements ---*/
                        elem[iElem]->SetNeighbor_Elements(Test_Elem,first_elem_face);
                        /*--- Store the element for both elements ---*/
                        elem[Test_Elem]->SetNeighbor_Elements(iElem,second_elem_face);
                        
                    }
                }
            }
}

void CPhysicalGeometry::SetBoundVolume(void) {
	unsigned short cont, iMarker, iElem, iNode_Domain, iNode_Surface;
	unsigned long Point_Domain, Point_Surface, Point, iElem_Surface, iElem_Domain;
	bool CheckVol;
    
	for (iMarker = 0; iMarker < nMarker; iMarker++)
		for (iElem_Surface = 0; iElem_Surface < nElem_Bound[iMarker]; iElem_Surface++) {
            
			/*--- Choose and arbitrary point from the surface --*/
			Point = bound[iMarker][iElem_Surface]->GetNode(0);
			CheckVol = false;
            
			for (iElem = 0; iElem < node[Point]->GetnElem(); iElem++) {
				/*--- Look for elements surronding that point --*/
				cont = 0; iElem_Domain = node[Point]->GetElem(iElem);
				for (iNode_Domain = 0; iNode_Domain < elem[iElem_Domain]->GetnNodes(); iNode_Domain++) {
					Point_Domain = elem[iElem_Domain]->GetNode(iNode_Domain);
					for (iNode_Surface = 0; iNode_Surface < bound[iMarker][iElem_Surface]->GetnNodes(); iNode_Surface++) {
						Point_Surface = bound[iMarker][iElem_Surface]->GetNode(iNode_Surface);
						if (Point_Surface == Point_Domain) cont++;
						if (cont == bound[iMarker][iElem_Surface]->GetnNodes()) break;
					}
					if (cont == bound[iMarker][iElem_Surface]->GetnNodes()) break;
				}
                
				if (cont == bound[iMarker][iElem_Surface]->GetnNodes()) {
					bound[iMarker][iElem_Surface]->SetDomainElement(iElem_Domain);
					CheckVol = true;
					break;
				}
			}
			if (!CheckVol) {
				cout << "The surface element ("<< iMarker <<", "<< iElem_Surface << ") doesn't have an associated volume element." << endl;
				cout << "Press any key to exit..." << endl;
				cin.get();
				exit(1);
			}
		}
}

void CPhysicalGeometry::SetVertex(CConfig *config) {
	unsigned long  iPoint, iVertex, iElem;
	unsigned short iMarker, iNode;
    
	/*--- Initialize the Vertex vector for each node of the grid ---*/
	for (iPoint = 0; iPoint < nPoint; iPoint++)
		for (iMarker = 0; iMarker < nMarker; iMarker++)
			node[iPoint]->SetVertex(-1,iMarker);
    
	/*--- Create and compute the vector with the number of vertex per marker ---*/
	nVertex = new unsigned long [nMarker];
	for (iMarker = 0; iMarker < nMarker; iMarker++) {
		/*--- Initialize the number of Bound Vertex for each Marker ---*/
		nVertex[iMarker] = 0;
		for (iElem = 0; iElem < nElem_Bound[iMarker]; iElem++)
			for(iNode = 0; iNode < bound[iMarker][iElem]->GetnNodes(); iNode++) {
				iPoint = bound[iMarker][iElem]->GetNode(iNode);
				/*--- Set the vertex in the node information ---*/
				if ((node[iPoint]->GetVertex(iMarker) == -1) || (config->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE)) {
					iVertex = nVertex[iMarker];
					node[iPoint]->SetVertex(nVertex[iMarker],iMarker);
					nVertex[iMarker]++;
				}
			}
	}
    
	/*--- Initialize the Vertex vector for each node, the previous result is deleted ---*/
	for (iPoint = 0; iPoint < nPoint; iPoint++)
		for (iMarker = 0; iMarker < nMarker; iMarker++)
			node[iPoint]->SetVertex(-1,iMarker);
    
	/*--- Create the bound vertex structure, note that the order
	 is the same as in the input file, this is important for Send/Receive part ---*/
	vertex = new CVertex**[nMarker];
	for (iMarker = 0; iMarker < nMarker; iMarker++) {
		vertex[iMarker] = new CVertex* [nVertex[iMarker]];
		nVertex[iMarker] = 0;
        
		/*--- Initialize the number of Bound Vertex for each Marker ---*/
		for (iElem = 0; iElem < nElem_Bound[iMarker]; iElem++)
			for(iNode = 0; iNode < bound[iMarker][iElem]->GetnNodes(); iNode++) {
				iPoint = bound[iMarker][iElem]->GetNode(iNode);
				/*--- Set the vertex in the node information ---*/
				if ((node[iPoint]->GetVertex(iMarker) == -1) || (config->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE)){
					iVertex = nVertex[iMarker];
					vertex[iMarker][iVertex] = new CVertex(iPoint, nDim);
                    
					if (config->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) {
						vertex[iMarker][iVertex]->SetRotation_Type(bound[iMarker][iElem]->GetRotation_Type());
						vertex[iMarker][iVertex]->SetMatching_Zone(bound[iMarker][iElem]->GetMatching_Zone());
					}
					node[iPoint]->SetVertex(nVertex[iMarker],iMarker);
					nVertex[iMarker]++;
				}
			}
	}
}

void CPhysicalGeometry::SetCG(void) {
	unsigned short nNode, iDim, iMarker, iNode;
	unsigned long elem_poin, edge_poin, iElem, iEdge;
	double **Coord;
    
	/*--- Compute the center of gravity for elements ---*/
	for(iElem = 0; iElem<nElem; iElem++) {
		nNode = elem[iElem]->GetnNodes();
		Coord = new double* [nNode];
		/*--- Store the coordinates for all the element nodes ---*/
		for (iNode = 0; iNode < nNode; iNode++) {
			elem_poin = elem[iElem]->GetNode(iNode);
			Coord[iNode] = new double [nDim];
			for (iDim = 0; iDim < nDim; iDim++)
				Coord[iNode][iDim]=node[elem_poin]->GetCoord(iDim);
		}
		/*--- Compute the element CG coordinates ---*/
		elem[iElem]->SetCG(Coord);
        
		for (iNode = 0; iNode < nNode; iNode++)
			if (Coord[iNode] != NULL) delete[] Coord[iNode];
		if (Coord != NULL) delete[] Coord;
	}
    
	/*--- Center of gravity for face elements ---*/
	for(iMarker = 0; iMarker < nMarker; iMarker++)
		for(iElem = 0; iElem < nElem_Bound[iMarker]; iElem++) {
			nNode = bound[iMarker][iElem]->GetnNodes();
			Coord = new double* [nNode];
			/*--- Store the coordinates for all the element nodes ---*/
			for (iNode = 0; iNode < nNode; iNode++) {
				elem_poin = bound[iMarker][iElem]->GetNode(iNode);
				Coord[iNode] = new double [nDim];
				for (iDim = 0; iDim < nDim; iDim++)
					Coord[iNode][iDim]=node[elem_poin]->GetCoord(iDim);
			}
			/*--- Compute the element CG coordinates ---*/
			bound[iMarker][iElem]->SetCG(Coord);
			for (iNode=0; iNode < nNode; iNode++)
				if (Coord[iNode] != NULL) delete[] Coord[iNode];
			if (Coord != NULL) delete[] Coord;
		}
    
	/*--- Center of gravity for edges ---*/
	for (iEdge = 0; iEdge < nEdge; iEdge++) {
		nNode = edge[iEdge]->GetnNodes();
		Coord = new double* [nNode];
		/*--- Store the coordinates for all the element nodes ---*/
		for (iNode = 0; iNode < nNode; iNode++) {
			edge_poin=edge[iEdge]->GetNode(iNode);
			Coord[iNode] = new double [nDim];
			for (iDim = 0; iDim<nDim; iDim++)
				Coord[iNode][iDim]=node[edge_poin]->GetCoord(iDim);
		}
		/*--- Compute the edge CG coordinates ---*/
		edge[iEdge]->SetCG(Coord);
        
		for (iNode=0; iNode < nNode; iNode++)
			if (Coord[iNode] != NULL) delete[] Coord[iNode];
		if (Coord != NULL) delete[] Coord;
	}
}

void CPhysicalGeometry::SetBoundControlVolume(CConfig *config, unsigned short action) {
	unsigned short Neighbor_Node, iMarker, iNode, iNeighbor_Nodes, iDim;
	unsigned long Neighbor_Point, iVertex, iPoint, iElem;
    long iEdge;
    double Area, *NormalFace = NULL;
    
	/*--- Update values of faces of the edge ---*/
	if (action != ALLOCATE)
		for (iMarker = 0; iMarker < nMarker; iMarker++)
			for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
				vertex[iMarker][iVertex]->SetZeroValues();
    
	double *Coord_Edge_CG = new double [nDim];
	double *Coord_Elem_CG = new double [nDim];
	double *Coord_Vertex = new double [nDim];
    
	/*--- Loop over all the markers ---*/
	for (iMarker = 0; iMarker < nMarker; iMarker++)
    /*--- Loop over all the boundary elements ---*/
		for (iElem = 0; iElem < nElem_Bound[iMarker]; iElem++)
        /*--- Loop over all the nodes of the boundary ---*/
			for(iNode = 0; iNode < bound[iMarker][iElem]->GetnNodes(); iNode++) {
				iPoint = bound[iMarker][iElem]->GetNode(iNode);
				iVertex = node[iPoint]->GetVertex(iMarker);
				/*--- Loop over the neighbor nodes, there is a face for each one ---*/
				for(iNeighbor_Nodes = 0; iNeighbor_Nodes < bound[iMarker][iElem]->GetnNeighbor_Nodes(iNode); iNeighbor_Nodes++) {
					Neighbor_Node = bound[iMarker][iElem]->GetNeighbor_Nodes(iNode,iNeighbor_Nodes);
					Neighbor_Point = bound[iMarker][iElem]->GetNode(Neighbor_Node);
					/*--- Shared edge by the Neighbor Point and the point ---*/
					iEdge = FindEdge(iPoint, Neighbor_Point);
					for (iDim = 0; iDim < nDim; iDim++) {
						Coord_Edge_CG[iDim] = edge[iEdge]->GetCG(iDim);
						Coord_Elem_CG[iDim] = bound[iMarker][iElem]->GetCG(iDim);
						Coord_Vertex[iDim] = node[iPoint]->GetCoord(iDim);
					}
					switch (nDim) {
                        case 2:
                            /*--- Store the 2D face ---*/
                            if (iNode == 0) vertex[iMarker][iVertex]->SetNodes_Coord(Coord_Elem_CG, Coord_Vertex, config);
                            if (iNode == 1) vertex[iMarker][iVertex]->SetNodes_Coord(Coord_Vertex, Coord_Elem_CG, config);
                            break;
                        case 3:
                            /*--- Store the 3D face ---*/
                            if (iNeighbor_Nodes == 0) vertex[iMarker][iVertex]->SetNodes_Coord(Coord_Elem_CG, Coord_Edge_CG, Coord_Vertex, config);
                            if (iNeighbor_Nodes == 1) vertex[iMarker][iVertex]->SetNodes_Coord(Coord_Edge_CG, Coord_Elem_CG, Coord_Vertex, config);
                            break;
					}
				}
			}
    
	delete[] Coord_Edge_CG;
	delete[] Coord_Elem_CG;
	delete[] Coord_Vertex;
    
    /*--- Check if there is a normal with null area ---*/
	for (iMarker = 0; iMarker < nMarker; iMarker ++)
		for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
            NormalFace = vertex[iMarker][iVertex]->GetNormal();
            Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += NormalFace[iDim]*NormalFace[iDim];
            Area = sqrt(Area);
            if (Area == 0.0) for (iDim = 0; iDim < nDim; iDim++) NormalFace[iDim] = EPS*EPS;
        }
    
}

void CPhysicalGeometry::MatchNearField(CConfig *config) {
	double epsilon = 1e-1;
    
    unsigned short nMarker_NearfieldBound = config->GetnMarker_NearFieldBound();
    
    if (nMarker_NearfieldBound != 0) {
        
        unsigned short iMarker, jMarker;
        unsigned long iVertex, iPoint, jVertex, jPoint = 0, pPoint = 0;
        double *Coord_i, *Coord_j, dist = 0.0, mindist, maxdist;
        
        cout << "Set Near-Field boundary conditions. " <<endl;
        
        maxdist = 0.0;
        for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
            if (config->GetMarker_All_Boundary(iMarker) == NEARFIELD_BOUNDARY) {
                for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
                    iPoint = vertex[iMarker][iVertex]->GetNode();
                    Coord_i = node[iPoint]->GetCoord();
                    
                    mindist = 1e10;
                    for (jMarker = 0; jMarker < config->GetnMarker_All(); jMarker++)
                        if ((config->GetMarker_All_Boundary(jMarker) == NEARFIELD_BOUNDARY) && (iMarker != jMarker))
                            for (jVertex = 0; jVertex < nVertex[jMarker]; jVertex++) {
                                jPoint = vertex[jMarker][jVertex]->GetNode();
                                Coord_j = node[jPoint]->GetCoord();
                                if (nDim == 2) dist = sqrt(pow(Coord_j[0]-Coord_i[0],2.0) + pow(Coord_j[1]-Coord_i[1],2.0));
                                if (nDim == 3) dist = sqrt(pow(Coord_j[0]-Coord_i[0],2.0) + pow(Coord_j[1]-Coord_i[1],2.0) + pow(Coord_j[2]-Coord_i[2],2.0));
                                if (dist < mindist) { mindist = dist; pPoint = jPoint; }
                            }
                    maxdist = max(maxdist, mindist);
                    vertex[iMarker][iVertex]->SetDonorPoint(pPoint);
                    
                    if (mindist > epsilon) {
                        cout.precision(10);
                        cout << endl;
                        cout << "   Bad match for point " << iPoint << ".\tNearest";
                        cout << " donor distance: " << scientific << mindist << ".";
                        vertex[iMarker][iVertex]->SetDonorPoint(iPoint);
                        maxdist = min(maxdist, 0.0);
                    }
                }
                cout <<"The max distance between points is: " << maxdist <<"."<< endl;
            }
        
    }
    
}

void CPhysicalGeometry::MatchInterface(CConfig *config) {
	double epsilon = 1.5e-1;
    
    unsigned short nMarker_InterfaceBound = config->GetnMarker_InterfaceBound();
    
    if (nMarker_InterfaceBound != 0) {
        
        unsigned short iMarker, jMarker;
        unsigned long iVertex, iPoint, jVertex, jPoint = 0, pPoint = 0;
        double *Coord_i, *Coord_j, dist = 0.0, mindist, maxdist;
        
        cout << "Set Interface boundary conditions." << endl;
        
        maxdist = 0.0;
        for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
            if (config->GetMarker_All_Boundary(iMarker) == INTERFACE_BOUNDARY) {
                for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
                    iPoint = vertex[iMarker][iVertex]->GetNode();
                    Coord_i = node[iPoint]->GetCoord();
                    
                    mindist = 1E6;
                    for (jMarker = 0; jMarker < config->GetnMarker_All(); jMarker++)
                        if ((config->GetMarker_All_Boundary(jMarker) == INTERFACE_BOUNDARY) && (iMarker != jMarker))
                            for (jVertex = 0; jVertex < nVertex[jMarker]; jVertex++) {
                                jPoint = vertex[jMarker][jVertex]->GetNode();
                                Coord_j = node[jPoint]->GetCoord();
                                if (nDim == 2) dist = sqrt(pow(Coord_j[0]-Coord_i[0],2.0) + pow(Coord_j[1]-Coord_i[1],2.0));
                                if (nDim == 3) dist = sqrt(pow(Coord_j[0]-Coord_i[0],2.0) + pow(Coord_j[1]-Coord_i[1],2.0) + pow(Coord_j[2]-Coord_i[2],2.0));
                                if (dist < mindist) {mindist = dist; pPoint = jPoint;}
                            }
                    maxdist = max(maxdist, mindist);
                    vertex[iMarker][iVertex]->SetDonorPoint(pPoint);
                    
                    if (mindist > epsilon) {
                        cout.precision(10);
                        cout << endl;
                        cout << "   Bad match for point " << iPoint << ".\tNearest";
                        cout << " donor distance: " << scientific << mindist << ".";
                        vertex[iMarker][iVertex]->SetDonorPoint(iPoint);
                        maxdist = min(maxdist, 0.0);
                    }
                    
                }
                cout <<"The max distance between points is: " << maxdist <<"."<< endl;
            }
        
    }
}

void CPhysicalGeometry::MatchZone(CConfig *config, CGeometry *geometry_donor, CConfig *config_donor,
                                  unsigned short val_iZone, unsigned short val_nZone) {
    
	unsigned short iMarker, jMarker;
	unsigned long iVertex, iPoint, jVertex, jPoint = 0, pPoint = 0;
	double *Coord_i, *Coord_j, dist = 0.0, mindist, maxdist;
    
	if (val_iZone == ZONE_0) cout << "Set zone boundary conditions (if any)." << endl;
    
	maxdist = 0.0;
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
		if (config->GetMarker_All_Boundary(iMarker) != SLIDING_INTERFACE) {
			for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
				iPoint = vertex[iMarker][iVertex]->GetNode();
				Coord_i = node[iPoint]->GetCoord();
                
				mindist = 1E6;
				for (jMarker = 0; jMarker < config_donor->GetnMarker_All(); jMarker++)
					for (jVertex = 0; jVertex < geometry_donor->GetnVertex(jMarker); jVertex++) {
						jPoint = geometry_donor->vertex[jMarker][jVertex]->GetNode();
						Coord_j = geometry_donor->node[jPoint]->GetCoord();
						if (nDim == 2) dist = sqrt(pow(Coord_j[0]-Coord_i[0],2.0) + pow(Coord_j[1]-Coord_i[1],2.0));
						if (nDim == 3) dist = sqrt(pow(Coord_j[0]-Coord_i[0],2.0) + pow(Coord_j[1]-Coord_i[1],2.0) + pow(Coord_j[2]-Coord_i[2],2.0));
						if (dist < mindist) { mindist = dist; pPoint = jPoint; }
					}
                
				maxdist = max(maxdist, mindist);
				vertex[iMarker][iVertex]->SetDonorPoint(pPoint);
                
			}
		}
	}
    
}


void CPhysicalGeometry::SetControlVolume(CConfig *config, unsigned short action) {
	unsigned long face_iPoint = 0, face_jPoint = 0, iPoint, iElem;
    long iEdge;
	unsigned short nEdgesFace = 1, iFace, iEdgesFace, iDim;
	double *Coord_Edge_CG, *Coord_FaceElem_CG, *Coord_Elem_CG, *Coord_FaceiPoint, *Coord_FacejPoint, Area,
	Volume, DomainVolume, my_DomainVolume, *NormalFace = NULL;
	bool change_face_orientation;
    
	int rank = MASTER_NODE;
    
	/*--- Update values of faces of the edge ---*/
	if (action != ALLOCATE) {
		for(iEdge = 0; iEdge < nEdge; iEdge++)
			edge[iEdge]->SetZeroValues();
		for(iPoint = 0; iPoint < nPoint; iPoint++)
			node[iPoint]->SetVolume (0.0);
	}
    
	Coord_Edge_CG = new double [nDim];
	Coord_FaceElem_CG = new double [nDim];
	Coord_Elem_CG = new double [nDim];
	Coord_FaceiPoint = new double [nDim];
	Coord_FacejPoint = new double [nDim];
    
	my_DomainVolume = 0.0;
	for(iElem = 0; iElem < nElem; iElem++)
		for (iFace = 0; iFace < elem[iElem]->GetnFaces(); iFace++) {
            
			/*--- In 2D all the faces have only one edge ---*/
			if (nDim == 2) nEdgesFace = 1;
			/*--- In 3D the number of edges per face is the same as the number of point per face ---*/
			if (nDim == 3) nEdgesFace = elem[iElem]->GetnNodesFace(iFace);
            
			/*-- Loop over the edges of a face ---*/
			for (iEdgesFace = 0; iEdgesFace < nEdgesFace; iEdgesFace++) {
                
				/*--- In 2D only one edge (two points) per edge ---*/
				if (nDim == 2) {
					face_iPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace,0));
					face_jPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace,1));
				}
                
				/*--- In 3D there are several edges in each face ---*/
				if (nDim == 3) {
					face_iPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace,iEdgesFace));
					if (iEdgesFace != nEdgesFace-1)
						face_jPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace,iEdgesFace+1));
					else
						face_jPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace,0));
				}
                
				/*--- We define a direction (from the smalest index to the greatest) --*/
				change_face_orientation = false;
				if (face_iPoint > face_jPoint) change_face_orientation = true;
				iEdge = FindEdge(face_iPoint, face_jPoint);
                
				for (iDim = 0; iDim < nDim; iDim++) {
					Coord_Edge_CG[iDim] = edge[iEdge]->GetCG(iDim);
					Coord_Elem_CG[iDim] = elem[iElem]->GetCG(iDim);
					Coord_FaceElem_CG[iDim] = elem[iElem]->GetFaceCG(iFace,iDim);
					Coord_FaceiPoint[iDim] = node[face_iPoint]->GetCoord(iDim);
					Coord_FacejPoint[iDim] = node[face_jPoint]->GetCoord(iDim);
				}
                
				switch (nDim) {
                    case 2:
                        /*--- Two dimensional problem ---*/
                        if (change_face_orientation) edge[iEdge]->SetNodes_Coord(Coord_Elem_CG, Coord_Edge_CG, config);
                        else edge[iEdge]->SetNodes_Coord(Coord_Edge_CG,Coord_Elem_CG, config);
                        Area = edge[iEdge]->GetVolume(Coord_FaceiPoint,Coord_Edge_CG,Coord_Elem_CG);
                        node[face_iPoint]->AddVolume(Area); my_DomainVolume +=Area;
                        Area = edge[iEdge]->GetVolume(Coord_FacejPoint,Coord_Edge_CG,Coord_Elem_CG);
                        node[face_jPoint]->AddVolume(Area); my_DomainVolume +=Area;
                        break;
                    case 3:
                        /*--- Three dimensional problem ---*/
                        if (change_face_orientation) edge[iEdge]->SetNodes_Coord(Coord_FaceElem_CG,Coord_Edge_CG,Coord_Elem_CG, config);
                        else edge[iEdge]->SetNodes_Coord(Coord_Edge_CG,Coord_FaceElem_CG,Coord_Elem_CG, config);
                        Volume = edge[iEdge]->GetVolume(Coord_FaceiPoint,Coord_Edge_CG,Coord_FaceElem_CG, Coord_Elem_CG);
                        node[face_iPoint]->AddVolume(Volume); my_DomainVolume +=Volume;
                        Volume = edge[iEdge]->GetVolume(Coord_FacejPoint,Coord_Edge_CG,Coord_FaceElem_CG, Coord_Elem_CG);
                        node[face_jPoint]->AddVolume(Volume); my_DomainVolume +=Volume;
                        break;
				}
			}
		}
    
    /*--- Check if there is a normal with null area ---*/
    for (iEdge = 0; iEdge < nEdge; iEdge++) {
        NormalFace = edge[iEdge]->GetNormal();
        Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += NormalFace[iDim]*NormalFace[iDim];
        Area = sqrt(Area);
        if (Area == 0.0) for (iDim = 0; iDim < nDim; iDim++) NormalFace[iDim] = EPS*EPS;
    }
    
	//	/*--- Set the volume for the iterations n and n-1 (dual time stteping with grid movement) ---*/
	//	if (config->GetUnsteady_Simulation() != NO) {
	//		for (iPoint = 0; iPoint < nPoint; iPoint++) {
	//			node[iPoint]->SetVolume_n();
	//			node[iPoint]->SetVolume_nM1();
	//		}
	//	}
    
	DomainVolume = my_DomainVolume;
    
	if ((rank == MASTER_NODE) && (action == ALLOCATE)) {
		if (nDim == 2) cout <<"Area of the computational grid: "<< DomainVolume <<"."<<endl;
		if (nDim == 3) cout <<"Volume of the computational grid: "<< DomainVolume <<"."<<endl;
	}
    
	config->SetDomainVolume(DomainVolume);
    
	delete[] Coord_Edge_CG;
	delete[] Coord_FaceElem_CG;
	delete[] Coord_Elem_CG;
	delete[] Coord_FaceiPoint;
	delete[] Coord_FacejPoint;
}

void CPhysicalGeometry::SetMeshFile (CConfig *config, string val_mesh_out_filename) {
	unsigned long iElem, iPoint, iElem_Bound;
	unsigned short iMarker, iNodes, iDim;
	unsigned short iPeriodic, nPeriodic = 0;
	ofstream output_file;
	string Grid_Marker;
	char *cstr;
	double *center, *angles, *transl;
    
	cstr = new char [val_mesh_out_filename.size()+1];
	strcpy (cstr, val_mesh_out_filename.c_str());
    
	/*--- Open .su2 grid file ---*/
	output_file.precision(15);
	output_file.open(cstr, ios::out);
    
	/*--- Write dimension, number of elements and number of points ---*/
	output_file << "NDIME= " << nDim << endl;
	output_file << "NELEM= " << nElem << endl;
	for (iElem = 0; iElem < nElem; iElem++) {
		output_file << elem[iElem]->GetVTK_Type();
		for (iNodes = 0; iNodes < elem[iElem]->GetnNodes(); iNodes++)
			output_file << "\t" << elem[iElem]->GetNode(iNodes);
		output_file << "\t"<<iElem<<endl;
	}
    
	/*--- Write the node coordinates ---*/
	output_file << "NPOIN= " << nPoint << "\t" << nPointDomain << endl;
	output_file.precision(15);
	for (iPoint = 0; iPoint < nPoint; iPoint++) {
		for (iDim = 0; iDim < nDim; iDim++)
			output_file << scientific << "\t" << node[iPoint]->GetCoord(iDim) ;
		output_file << "\t" << iPoint << endl;
        
	}
    
	/*--- Loop through and write the boundary info ---*/
	output_file << "NMARK= " << nMarker << endl;
	for (iMarker = 0; iMarker < nMarker; iMarker++) {
        
		/*--- Ignore SEND_RECEIVE for the moment ---*/
		if (bound[iMarker][0]->GetVTK_Type() != VERTEX) {
            
			Grid_Marker = config->GetMarker_All_Tag(iMarker);
			output_file << "MARKER_TAG= " << Grid_Marker <<endl;
			output_file << "MARKER_ELEMS= " << nElem_Bound[iMarker]<< endl;
            
			if (nDim == 2) {
				for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
					output_file << bound[iMarker][iElem_Bound]->GetVTK_Type() << "\t" ;
					for (iNodes = 0; iNodes < bound[iMarker][iElem_Bound]->GetnNodes(); iNodes++)
						output_file << bound[iMarker][iElem_Bound]->GetNode(iNodes) << "\t" ;
					output_file	<< iElem_Bound << endl;
				}
			}
            
			if (nDim == 3) {
				for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
					output_file << bound[iMarker][iElem_Bound]->GetVTK_Type() << "\t" ;
					for (iNodes = 0; iNodes < bound[iMarker][iElem_Bound]->GetnNodes(); iNodes++)
						output_file << bound[iMarker][iElem_Bound]->GetNode(iNodes) << "\t" ;
					output_file	<< iElem_Bound << endl;
				}
			}
            
		} else if (bound[iMarker][0]->GetVTK_Type() == VERTEX) {
			output_file << "MARKER_TAG= SEND_RECEIVE" << endl;
			output_file << "MARKER_ELEMS= " << nElem_Bound[iMarker]<< endl;
			if (config->GetMarker_All_SendRecv(iMarker) > 0) output_file << "SEND_TO= " << config->GetMarker_All_SendRecv(iMarker) << endl;
			if (config->GetMarker_All_SendRecv(iMarker) < 0) output_file << "SEND_TO= " << config->GetMarker_All_SendRecv(iMarker) << endl;
            
			for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
				output_file << bound[iMarker][iElem_Bound]->GetVTK_Type() << "\t" <<
                bound[iMarker][iElem_Bound]->GetNode(0) << "\t" <<
                bound[iMarker][iElem_Bound]->GetRotation_Type() << "\t" <<
                bound[iMarker][iElem_Bound]->GetMatching_Zone()<< endl;
			}
            
		}
	}
    
	/*--- Get the total number of periodic transformations ---*/
	nPeriodic = config->GetnPeriodicIndex();
	output_file << "NPERIODIC= " << nPeriodic << endl;
    
	/*--- From iPeriodic obtain the iMarker ---*/
	for (iPeriodic = 0; iPeriodic < nPeriodic; iPeriodic++) {
        
		/*--- Retrieve the supplied periodic information. ---*/
		center = config->GetPeriodicCenter(iPeriodic);
		angles = config->GetPeriodicRotation(iPeriodic);
		transl = config->GetPeriodicTranslate(iPeriodic);
        
		output_file << "PERIODIC_INDEX= " << iPeriodic << endl;
		output_file << center[0] << "\t" << center[1] << "\t" << center[2] << endl;
		output_file << angles[0] << "\t" << angles[1] << "\t" << angles[2] << endl;
		output_file << transl[0] << "\t" << transl[1] << "\t" << transl[2] << endl;
        
	}
    
    
	output_file.close();
}

void CPhysicalGeometry::SetMeshFile(CConfig *config, string val_mesh_out_filename, string val_mesh_in_filename) {
	unsigned long iElem, iPoint, iElem_Bound, nElem_, nElem_Bound_, vnodes_edge[2], vnodes_triangle[3], vnodes_quad[4], vnodes_tetra[4], vnodes_hexa[8], vnodes_wedge[6], vnodes_pyramid[5], vnodes_vertex;
	unsigned short iMarker, iDim, iChar, iPeriodic, nPeriodic = 0, VTK_Type, nDim_, nMarker_, transform, matching_zone = 0;
    char *cstr;
	double *center, *angles, *transl;
    long SendRecv;
	ofstream output_file;
    ifstream input_file;
	string Grid_Marker, text_line, Marker_Tag;
    string::size_type position;
    
	/*--- Open output_file .su2 grid file ---*/
    cstr = new char [val_mesh_out_filename.size()+1];
	strcpy (cstr, val_mesh_out_filename.c_str());
	output_file.precision(15);
	output_file.open(cstr, ios::out);
    
    /*--- Open input_file .su2 grid file ---*/
    cstr = new char [val_mesh_in_filename.size()+1];
	strcpy (cstr, val_mesh_in_filename.c_str());
    input_file.open(cstr, ios::out);
    
    /*--- Read grid file with format SU2 ---*/
    while (getline (input_file, text_line)) {
        
        /*--- Read the dimension of the problem ---*/
        position = text_line.find ("NDIME=",0);
        if (position != string::npos) {
            text_line.erase (0,6); nDim_ = atoi(text_line.c_str());
            output_file << "NDIME= " << nDim_ << endl;
        }
        
        /*--- Read the information about inner elements ---*/
        position = text_line.find ("NELEM=",0);
        if (position != string::npos) {
            text_line.erase (0,6); nElem_ = atoi(text_line.c_str());
            output_file << "NELEM= " << nElem_ << endl;
            
            
            /*--- Loop over all the volumetric elements ---*/
            for (iElem = 0; iElem < nElem_;  iElem++) {
                getline(input_file, text_line);
                istringstream elem_line(text_line);
                
                elem_line >> VTK_Type;
                output_file << VTK_Type;
                
                switch(VTK_Type) {
                    case TRIANGLE:
                        elem_line >> vnodes_triangle[0]; elem_line >> vnodes_triangle[1]; elem_line >> vnodes_triangle[2];
                        output_file << "\t" << vnodes_triangle[0] << "\t" << vnodes_triangle[1] << "\t" << vnodes_triangle[2] << endl;
                        break;
                    case RECTANGLE:
                        elem_line >> vnodes_quad[0]; elem_line >> vnodes_quad[1]; elem_line >> vnodes_quad[2]; elem_line >> vnodes_quad[3];
                        output_file << "\t" << vnodes_quad[0] << "\t" << vnodes_quad[1] << "\t" << vnodes_quad[2] << "\t" << vnodes_quad[3] << endl;
                        break;
                    case TETRAHEDRON:
                        elem_line >> vnodes_tetra[0]; elem_line >> vnodes_tetra[1]; elem_line >> vnodes_tetra[2]; elem_line >> vnodes_tetra[3];
                        output_file << "\t" << vnodes_tetra[0] << "\t" << vnodes_tetra[1] << "\t" << vnodes_tetra[2] << "\t" << vnodes_tetra[3] << endl;
                        break;
                    case HEXAHEDRON:
                        elem_line >> vnodes_hexa[0]; elem_line >> vnodes_hexa[1]; elem_line >> vnodes_hexa[2];
                        elem_line >> vnodes_hexa[3]; elem_line >> vnodes_hexa[4]; elem_line >> vnodes_hexa[5];
                        elem_line >> vnodes_hexa[6]; elem_line >> vnodes_hexa[7];
                        output_file << "\t" << vnodes_hexa[0] << "\t" << vnodes_hexa[1] << "\t" << vnodes_hexa[2] << "\t" << vnodes_hexa[3] << "\t" << vnodes_hexa[4] << "\t" << vnodes_hexa[5] << "\t" << vnodes_hexa[6] << "\t" << vnodes_hexa[7] << endl;
                        break;
                    case WEDGE:
                        elem_line >> vnodes_wedge[0]; elem_line >> vnodes_wedge[1]; elem_line >> vnodes_wedge[2];
                        elem_line >> vnodes_wedge[3]; elem_line >> vnodes_wedge[4]; elem_line >> vnodes_wedge[5];
                        output_file << "\t" << vnodes_wedge[0] << "\t" << vnodes_wedge[1] << "\t" << vnodes_wedge[2] << "\t" << vnodes_wedge[3] << "\t" << vnodes_wedge[4] << "\t" << vnodes_wedge[5] << endl;
                        break;
                    case PYRAMID:
                        elem_line >> vnodes_pyramid[0]; elem_line >> vnodes_pyramid[1]; elem_line >> vnodes_pyramid[2];
                        elem_line >> vnodes_pyramid[3]; elem_line >> vnodes_pyramid[4];
                        output_file << "\t" << vnodes_pyramid[0] << "\t" << vnodes_pyramid[1] << "\t" << vnodes_pyramid[2] << "\t" << vnodes_pyramid[3] << "\t" << vnodes_pyramid[4] << endl;
                        break;
                }
            }
        }
        
        /*--- Coordinates ---*/
        position = text_line.find ("NPOIN=",0);
        if (position != string::npos) {
            
            /*--- Skip the lines about the points ---*/
            for (iPoint = 0; iPoint < nPoint;  iPoint++) {
                getline(input_file, text_line);
            }
            
            /*--- Add the new coordinates ---*/
            output_file << "NPOIN= " << nPoint << "\t" << nPointDomain << endl;
            for (iPoint = 0; iPoint < nPoint; iPoint++) {
                for (iDim = 0; iDim < nDim; iDim++)
                    output_file << scientific << node[iPoint]->GetCoord(iDim) << "\t";
                output_file << iPoint << endl;
            }
            
        }
        
        /*--- Write the physical boundaries ---*/
        position = text_line.find ("NMARK=",0);
        if (position != string::npos) {
            
            text_line.erase (0,6); nMarker_ = atoi(text_line.c_str());
            output_file << "NMARK= " << nMarker_ << endl;
            
            for (iMarker = 0 ; iMarker < nMarker_; iMarker++) {
                
                getline (input_file,text_line);
                text_line.erase (0,11);
                string::size_type position;
                for (iChar = 0; iChar < 20; iChar++) {
                    position = text_line.find( " ", 0 );
                    if(position != string::npos) text_line.erase (position,1);
                    position = text_line.find( "\r", 0 );
                    if(position != string::npos) text_line.erase (position,1);
                    position = text_line.find( "\n", 0 );
                    if(position != string::npos) text_line.erase (position,1);
                }
                Marker_Tag = text_line.c_str();
                
                /*--- Standart physical boundary ---*/
                if (Marker_Tag != "SEND_RECEIVE") {
                    
                    getline (input_file, text_line);
                    
                    text_line.erase (0,13); nElem_Bound_ = atoi(text_line.c_str());
                    output_file << "MARKER_TAG= " << Marker_Tag << endl;
                    output_file << "MARKER_ELEMS= " << nElem_Bound_<< endl;
                    
                    for (iElem_Bound = 0; iElem_Bound < nElem_Bound_; iElem_Bound++) {
                        
                        getline(input_file, text_line);
                        istringstream bound_line(text_line);
                        
                        bound_line >> VTK_Type;
                        output_file << VTK_Type;
                        
                        switch(VTK_Type) {
                            case LINE:
                                bound_line >> vnodes_edge[0]; bound_line >> vnodes_edge[1];
                                output_file << "\t" << vnodes_edge[0] << "\t" << vnodes_edge[1] << endl;
                                break;
                            case TRIANGLE:
                                bound_line >> vnodes_triangle[0]; bound_line >> vnodes_triangle[1]; bound_line >> vnodes_triangle[2];
                                output_file << "\t" << vnodes_triangle[0] << "\t" << vnodes_triangle[1] << "\t" << vnodes_triangle[2] << endl;
                                break;
                            case RECTANGLE:
                                bound_line >> vnodes_quad[0]; bound_line >> vnodes_quad[1]; bound_line >> vnodes_quad[2]; bound_line >> vnodes_quad[3];
                                output_file << "\t" << vnodes_quad[0] << "\t" << vnodes_quad[1] << "\t" << vnodes_quad[2] << "\t" << vnodes_quad[3] << endl;
                                break;
                        }
                    }
                    
                }
                
                /*--- Send-Receive boundaries definition ---*/
                else {
                    output_file << "MARKER_TAG= SEND_RECEIVE" << endl;
                    getline (input_file,text_line);
                    text_line.erase (0,13); nElem_Bound_ = atoi(text_line.c_str());
                    output_file << "MARKER_ELEMS= " << nElem_Bound_ << endl;
                    getline (input_file, text_line); text_line.erase (0,8);
                    SendRecv = atoi(text_line.c_str());
                    output_file << "SEND_TO= " << SendRecv << endl;
                    
                    for (iElem_Bound = 0; iElem_Bound < nElem_Bound_; iElem_Bound++) {
                        getline(input_file,text_line);
                        istringstream bound_line(text_line);
                        bound_line >> VTK_Type; bound_line >> vnodes_vertex; bound_line >> transform;
                        output_file << VTK_Type << "\t" << vnodes_vertex << "\t" << transform << "\t" << matching_zone << endl;
                    }
                }
                
            }
        }
    }
    
    
	/*--- Get the total number of periodic transformations ---*/
	nPeriodic = config->GetnPeriodicIndex();
	output_file << "NPERIODIC= " << nPeriodic << endl;
    
	/*--- From iPeriodic obtain the iMarker ---*/
	for (iPeriodic = 0; iPeriodic < nPeriodic; iPeriodic++) {
        
		/*--- Retrieve the supplied periodic information. ---*/
		center = config->GetPeriodicCenter(iPeriodic);
		angles = config->GetPeriodicRotation(iPeriodic);
		transl = config->GetPeriodicTranslate(iPeriodic);
        
		output_file << "PERIODIC_INDEX= " << iPeriodic << endl;
		output_file << center[0] << "\t" << center[1] << "\t" << center[2] << endl;
		output_file << angles[0] << "\t" << angles[1] << "\t" << angles[2] << endl;
		output_file << transl[0] << "\t" << transl[1] << "\t" << transl[2] << endl;
        
	}
    
    input_file.close();
	output_file.close();
    
}

void CPhysicalGeometry::SetTecPlot(char mesh_filename[200]) {
	unsigned long iElem, iPoint;
	unsigned short iDim;
	ofstream Tecplot_File;
    
	Tecplot_File.open(mesh_filename, ios::out);
	Tecplot_File << "TITLE = \"Visualization of the volumetric grid\"" << endl;
    
	if (nDim == 2) {
		Tecplot_File << "VARIABLES = \"x\",\"y\" " << endl;
		Tecplot_File << "ZONE NODES= "<< nPoint <<", ELEMENTS= "<< nElem <<", DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL"<< endl;
	}
	if (nDim == 3) {
		Tecplot_File << "VARIABLES = \"x\",\"y\",\"z\" " << endl;
		Tecplot_File << "ZONE NODES= "<< nPoint <<", ELEMENTS= "<< nElem <<", DATAPACKING=POINT, ZONETYPE=FEBRICK"<< endl;
	}
    
	for(iPoint = 0; iPoint < nPoint; iPoint++) {
		for(iDim = 0; iDim < nDim; iDim++)
			Tecplot_File << scientific << node[iPoint]->GetCoord(iDim) << "\t";
		Tecplot_File << "\n";
	}
    
	for(iElem = 0; iElem < nElem; iElem++) {
		if (elem[iElem]->GetVTK_Type() == TRIANGLE) {
			Tecplot_File <<
            elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
            elem[iElem]->GetNode(2)+1 <<" "<< elem[iElem]->GetNode(2)+1 << endl;
		}
		if (elem[iElem]->GetVTK_Type() == RECTANGLE) {
			Tecplot_File <<
            elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
            elem[iElem]->GetNode(2)+1 <<" "<< elem[iElem]->GetNode(3)+1 << endl;
		}
		if (elem[iElem]->GetVTK_Type() == TETRAHEDRON) {
			Tecplot_File <<
            elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
            elem[iElem]->GetNode(2)+1 <<" "<< elem[iElem]->GetNode(2)+1 <<" "<<
            elem[iElem]->GetNode(3)+1 <<" "<< elem[iElem]->GetNode(3)+1 <<" "<<
            elem[iElem]->GetNode(3)+1 <<" "<< elem[iElem]->GetNode(3)+1 << endl;
		}
		if (elem[iElem]->GetVTK_Type() == HEXAHEDRON) {
			Tecplot_File <<
            elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
            elem[iElem]->GetNode(2)+1 <<" "<< elem[iElem]->GetNode(3)+1 <<" "<<
            elem[iElem]->GetNode(4)+1 <<" "<< elem[iElem]->GetNode(5)+1 <<" "<<
            elem[iElem]->GetNode(6)+1 <<" "<< elem[iElem]->GetNode(7)+1 << endl;
		}
		if (elem[iElem]->GetVTK_Type() == PYRAMID) {
			Tecplot_File <<
            elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
            elem[iElem]->GetNode(2)+1 <<" "<< elem[iElem]->GetNode(3)+1 <<" "<<
            elem[iElem]->GetNode(4)+1 <<" "<< elem[iElem]->GetNode(4)+1 <<" "<<
            elem[iElem]->GetNode(4)+1 <<" "<< elem[iElem]->GetNode(4)+1 << endl;
		}
		if (elem[iElem]->GetVTK_Type() == WEDGE) {
			Tecplot_File <<
            elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
            elem[iElem]->GetNode(1)+1 <<" "<< elem[iElem]->GetNode(2)+1 <<" "<<
            elem[iElem]->GetNode(3)+1 <<" "<< elem[iElem]->GetNode(4)+1 <<" "<<
            elem[iElem]->GetNode(4)+1 <<" "<< elem[iElem]->GetNode(5)+1 << endl;
		}
	}
    
	Tecplot_File.close();
}

void CPhysicalGeometry::SetCoord_Smoothing (unsigned short val_nSmooth, double val_smooth_coeff, CConfig *config) {
	unsigned short iSmooth, nneigh, iMarker;
	double *Coord_Old, *Coord_Sum, *Coord, *Coord_i, *Coord_j, Position_Plane = 0.0;
	unsigned long iEdge, iPoint, jPoint, iVertex;
	double eps = 1E-6;
	bool NearField = false;
    
	Coord = new double [nDim];
    
	for (iPoint = 0; iPoint < GetnPoint(); iPoint++) {
		double *Coord = node[iPoint]->GetCoord();
		node[iPoint]->SetCoord_Old(Coord);
	}
    
	/*--- Jacobi iterations ---*/
	for (iSmooth = 0; iSmooth < val_nSmooth; iSmooth++) {
        
		for (iPoint = 0; iPoint < nPoint; iPoint++)
			node[iPoint]->SetCoord_SumZero();
        
        
		/*--- Loop over Interior edges ---*/
		for(iEdge = 0; iEdge < nEdge; iEdge++) {
			iPoint = edge[iEdge]->GetNode(0);
			Coord_i = node[iPoint]->GetCoord();
            
			jPoint = edge[iEdge]->GetNode(1);
			Coord_j = node[jPoint]->GetCoord();
            
			/*--- Accumulate nearest neighbor Coord to Res_sum for each variable ---*/
			node[iPoint]->AddCoord_Sum(Coord_j);
			node[jPoint]->AddCoord_Sum(Coord_i);
            
		}
        
		/*--- Loop over all mesh points (Update Coords with averaged sum) ---*/
		for (iPoint = 0; iPoint < nPoint; iPoint++) {
			nneigh = node[iPoint]->GetnPoint();
			Coord_Sum = node[iPoint]->GetCoord_Sum();
			Coord_Old = node[iPoint]->GetCoord_Old();
            
			if (nDim == 2) {
				Coord[0] =(Coord_Old[0] + val_smooth_coeff*Coord_Sum[0]) /(1.0 + val_smooth_coeff*double(nneigh));
				Coord[1] =(Coord_Old[1] + val_smooth_coeff*Coord_Sum[1]) /(1.0 + val_smooth_coeff*double(nneigh));
				if ((NearField) && ((Coord_Old[1] > Position_Plane-eps) && (Coord_Old[1] < Position_Plane+eps)))
					Coord[1] = Coord_Old[1];
			}
            
			if (nDim == 3) {
				Coord[0] =(Coord_Old[0] + val_smooth_coeff*Coord_Sum[0]) /(1.0 + val_smooth_coeff*double(nneigh));
				Coord[1] =(Coord_Old[1] + val_smooth_coeff*Coord_Sum[1]) /(1.0 + val_smooth_coeff*double(nneigh));
				Coord[2] =(Coord_Old[2] + val_smooth_coeff*Coord_Sum[2]) /(1.0 + val_smooth_coeff*double(nneigh));
				if ((NearField) && ((Coord_Old[2] > Position_Plane-eps) && (Coord_Old[2] < Position_Plane+eps)))
					Coord[2] = Coord_Old[2];
			}
            
			node[iPoint]->SetCoord(Coord);
		}
        
		/*--- Copy boundary values ---*/
		for (iMarker = 0; iMarker < nMarker; iMarker++)
			for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
				iPoint = vertex[iMarker][iVertex]->GetNode();
				Coord_Old = node[iPoint]->GetCoord_Old();
				node[iPoint]->SetCoord(Coord_Old);
			}
	}
    
	delete[] Coord;
}

bool CPhysicalGeometry::FindFace(unsigned long first_elem, unsigned long second_elem, unsigned short &face_first_elem,
                                 unsigned short &face_second_elem) {
    
	/*--- Find repeated nodes between two elements to identify the common face ---*/
	unsigned long iPoint = 0, jPoint = 0;
	unsigned short face_node, iFace, iNode, jNode, kNode, nNodesFace;
    vector<unsigned long> CommonPoints, PointFaceFirst, PointFaceSecond;
    vector<unsigned long>::iterator IterPoint;
    pair<vector <unsigned long>::iterator, vector <unsigned long>::iterator> mypair;
	bool face_first_found = false, face_second_found =false;
    
	if (first_elem == second_elem) return 0;
    
	kNode = 0;
	for (iNode = 0; iNode < elem[first_elem]->GetnNodes(); iNode++) {
		iPoint = elem[first_elem]->GetNode(iNode);
		for (jNode = 0; jNode < elem[second_elem]->GetnNodes(); jNode++) {
			jPoint = elem[second_elem]->GetNode(jNode);
			if (iPoint == jPoint) {
                CommonPoints.push_back(iPoint);
				break;
            }
		}
	}
    
	/*--- Sort point in face and check that the list is unique ---*/
    sort( CommonPoints.begin(), CommonPoints.end());
    IterPoint = unique( CommonPoints.begin(), CommonPoints.end());
    CommonPoints.resize( distance(CommonPoints.begin(), IterPoint) );
    
	/*--- Search the secuence in the first element ---*/
	for (iFace = 0; iFace < elem[first_elem]->GetnFaces(); iFace++) {
		nNodesFace = elem[first_elem]->GetnNodesFace(iFace);
		for (iNode = 0; iNode < nNodesFace; iNode++) {
			face_node = elem[first_elem]->GetFaces(iFace, iNode);
            PointFaceFirst.push_back(elem[first_elem]->GetNode(face_node));
		}
        
		/*--- Sort face_poin to perform comparison ---*/
        sort( PointFaceFirst.begin(), PointFaceFirst.end());
        
		/*--- List comparison ---*/
        mypair = mismatch (PointFaceFirst.begin(), PointFaceFirst.end(), CommonPoints.begin());
		if (mypair.first == PointFaceFirst.end()) {
            face_first_elem = iFace;
            face_first_found = true;
            break;
        }
        
        PointFaceFirst.erase (PointFaceFirst.begin(),PointFaceFirst.end());
	}
    
	/*--- Search the secuence in the second element ---*/
	for (iFace = 0; iFace < elem[second_elem]->GetnFaces(); iFace++) {
		nNodesFace = elem[second_elem]->GetnNodesFace(iFace);
		for (iNode = 0; iNode < nNodesFace; iNode++) {
			face_node = elem[second_elem]->GetFaces(iFace,iNode);
            PointFaceSecond.push_back(elem[second_elem]->GetNode(face_node));
		}
        
		/*--- Sort face_poin to perform comparison ---*/
        sort( PointFaceSecond.begin(), PointFaceSecond.end());
        
		/*--- List comparison ---*/
        mypair = mismatch (PointFaceSecond.begin(), PointFaceSecond.end(), CommonPoints.begin());
        if (mypair.first == PointFaceSecond.end()) {
            face_second_elem = iFace;
            face_second_found = true;
            break;
        }
        
        PointFaceSecond.erase (PointFaceSecond.begin(),PointFaceSecond.end());
	}
    
	if (face_first_found && face_second_found) return true;
	else return false;
    
}

void CPhysicalGeometry::SetBoundTecPlot (CConfig *config, char mesh_filename[200]) {
	ofstream Tecplot_File;
	unsigned long iPoint, Total_nElem_Bound, iElem, *PointSurface = NULL, nPointSurface = 0;
	unsigned short Coord_i, iMarker;
    
	/*--- It is important to do a renumering to don't add points
	 that do not belong to the surfaces ---*/
	PointSurface = new unsigned long[nPoint];
	for (iPoint = 0; iPoint < nPoint; iPoint++)
		if (node[iPoint]->GetBoundary()) {
			PointSurface[iPoint] = nPointSurface;
			nPointSurface++;
		}
    
	/*--- Compute the total number of elements ---*/
	Total_nElem_Bound = 0;
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
		if (config->GetMarker_All_Plotting(iMarker) == YES) {
			Total_nElem_Bound += nElem_Bound[iMarker];
		}
    }
    
	/*--- Open the tecplot file ---*/
	Tecplot_File.open(mesh_filename, ios::out);
	Tecplot_File << "TITLE = \"Visualization of the surface grid\"" << endl;
    
    if (Total_nElem_Bound != 0) {
        
        /*--- Write the header of the file ---*/
        if (nDim == 2) {
            Tecplot_File << "VARIABLES = \"x\",\"y\" " << endl;
            Tecplot_File << "ZONE NODES= "<< nPointSurface <<", ELEMENTS= "<< Total_nElem_Bound <<", DATAPACKING=POINT, ZONETYPE=FELINESEG"<< endl;
        }
        if (nDim == 3) {
            Tecplot_File << "VARIABLES = \"x\",\"y\",\"z\" " << endl;
            Tecplot_File << "ZONE NODES= "<< nPointSurface <<", ELEMENTS= "<< Total_nElem_Bound <<", DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL"<< endl;
        }
        
        /*--- Only write the coordiantes of the points that are on the surfaces ---*/
        if (nDim == 3) {
            for(iPoint = 0; iPoint < nPoint; iPoint++)
                if (node[iPoint]->GetBoundary()) {
                    for(Coord_i = 0; Coord_i < nDim-1; Coord_i++)
                        Tecplot_File << node[iPoint]->GetCoord(Coord_i) << " ";
                    Tecplot_File << node[iPoint]->GetCoord(nDim-1) << "\n";
                }
        }
        else {
            for(iPoint = 0; iPoint < nPoint; iPoint++)
                if (node[iPoint]->GetBoundary()){
                    for(Coord_i = 0; Coord_i < nDim; Coord_i++)
                        Tecplot_File << node[iPoint]->GetCoord(Coord_i) << " ";
                    Tecplot_File << "\n";
                }
        }
        
        /*--- Write the cells using the new numbering ---*/
        for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
            if (config->GetMarker_All_Plotting(iMarker) == YES)
                for(iElem = 0; iElem < nElem_Bound[iMarker]; iElem++) {
                    if (nDim == 2) {
                        Tecplot_File << PointSurface[bound[iMarker][iElem]->GetNode(0)]+1 << " "
                        << PointSurface[bound[iMarker][iElem]->GetNode(1)]+1 << endl;
                    }
                    if (nDim == 3) {
                        if (bound[iMarker][iElem]->GetnNodes() == 3) {
                            Tecplot_File << PointSurface[bound[iMarker][iElem]->GetNode(0)]+1 << " "
                            << PointSurface[bound[iMarker][iElem]->GetNode(1)]+1 << " "
                            << PointSurface[bound[iMarker][iElem]->GetNode(2)]+1 << " "
                            << PointSurface[bound[iMarker][iElem]->GetNode(2)]+1 << endl;
                        }
                        if (bound[iMarker][iElem]->GetnNodes() == 4) {
                            Tecplot_File << PointSurface[bound[iMarker][iElem]->GetNode(0)]+1 << " "
                            << PointSurface[bound[iMarker][iElem]->GetNode(1)]+1 << " "
                            << PointSurface[bound[iMarker][iElem]->GetNode(2)]+1 << " "
                            << PointSurface[bound[iMarker][iElem]->GetNode(3)]+1 << endl;
                        }
                    }
                }
    }
    else {
        /*--- No elements in the surface ---*/
        if (nDim == 2) {
            Tecplot_File << "VARIABLES = \"x\",\"y\" " << endl;
            Tecplot_File << "ZONE NODES= 1, ELEMENTS= 1, DATAPACKING=POINT, ZONETYPE=FELINESEG"<< endl;
            Tecplot_File << "0.0 0.0"<< endl;
            Tecplot_File << "1 1"<< endl;
        }
        if (nDim == 3) {
            Tecplot_File << "VARIABLES = \"x\",\"y\",\"z\" " << endl;
            Tecplot_File << "ZONE NODES= 1, ELEMENTS= 1, DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL"<< endl;
            Tecplot_File << "0.0 0.0 0.0"<< endl;
            Tecplot_File << "1 1 1 1"<< endl;
        }
    }
    
    
	/*--- Dealocate memory and close the file ---*/
	delete[] PointSurface;
	Tecplot_File.close();
}

void CPhysicalGeometry::SetBoundSTL (CConfig *config, char mesh_filename[200]) {
	ofstream STL_File;
	unsigned long this_node, iNode, nNode, iElem;
	unsigned short iDim, iMarker;
	double p[3], u[3], v[3], n[3], a;
    
	/*---	STL format:
     solid NAME
     ...
     facet normal 0.00 0.00 1.00
     outer loop
     vertex  2.00  2.00  0.00
     vertex -1.00  1.00  0.00
     vertex  0.00 -1.00  0.00
     endloop
     endfacet
     ...
     end solid
     --- */
    
	/*--- Open the STL file ---*/
	STL_File.open(mesh_filename, ios::out);
    
	/*--- Write the header of the file ---*/
	STL_File << "solid surface_mesh" << endl;
    
	/*--- Write facets of surface markers ---*/
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if (config->GetMarker_All_Plotting(iMarker) == YES)
			for(iElem = 0; iElem < nElem_Bound[iMarker]; iElem++) {
                
				/*--- number of nodes for this elemnt ---*/
				nNode = bound[iMarker][iElem]->GetnNodes();
                
				/*--- Calculate Normal Vector ---*/
				for (iDim=0; iDim<nDim; iDim++){
					p[0] = node[bound[iMarker][iElem]->GetNode(0)]      ->GetCoord(iDim);
					p[1] = node[bound[iMarker][iElem]->GetNode(1)]      ->GetCoord(iDim);
					p[2] = node[bound[iMarker][iElem]->GetNode(nNode-1)]->GetCoord(iDim);
					/*cout << p[0] <<endl;
                     cout << p[1] <<endl;
                     cout << p[2] <<endl;*/
					u[iDim] = p[1]-p[0];
					v[iDim] = p[2]-p[0];
				}
                
				n[0] = u[1]*v[2]-u[2]*v[1];
				n[1] = u[2]*v[0]-u[0]*v[2];
				n[2] = u[0]*v[1]-u[1]*v[0];
				a = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
				/*cout << n[0] <<endl;
                 cout << n[1] <<endl;
                 cout << n[2] <<endl;
                 cout << a << endl;*/
                
				/*--- Print normal vector ---*/
				STL_File << "  facet normal ";
				for (iDim=0; iDim<nDim; iDim++){
					STL_File << n[iDim]/a << " ";
				}
				STL_File << endl;
                
				/*--- STL Facet Loop --*/
				STL_File << "    outer loop" << endl;
                
				/*--- Print Nodes for Facet ---*/
				for(iNode=0; iNode<nNode; iNode++) {
					this_node = bound[iMarker][iElem]->GetNode(iNode);
					STL_File << "      vertex ";
					for (iDim = 0; iDim < nDim; iDim++)
						STL_File << node[this_node]->GetCoord(iDim) << " ";
					if (nDim==2)
						STL_File << 0.0 << " ";
					STL_File <<  endl;
				}
				STL_File << "    endloop" << endl;
				STL_File << "  endfacet" << endl;
			}
    
	/*--- Done with Surface Mesh ---*/
	STL_File << "endsolid" << endl;
    
	/*--- Close the file ---*/
	STL_File.close();
}

void CPhysicalGeometry::SetColorGrid(CConfig *config) {
    
}

void CPhysicalGeometry::GetQualityStatistics(double *statistics) {
	unsigned long jPoint, Point_2, Point_3, iElem;
	double *Coord_j, *Coord_2, *Coord_3;
	unsigned short iDim;
    
	statistics[0] = 1e06;
	statistics[1] = 0;
    
	/*--- Loop interior edges ---*/
	for (iElem = 0; iElem < this->GetnElem(); iElem++) {
        
		if ((this->GetnDim() == 2) && (elem[iElem]->GetVTK_Type() == TRIANGLE)) {
            
			jPoint = elem[iElem]->GetNode(0); Coord_j = node[jPoint]->GetCoord();
			Point_2 = elem[iElem]->GetNode(1); Coord_2 = node[Point_2]->GetCoord();
			Point_3 = elem[iElem]->GetNode(2); Coord_3 = node[Point_3]->GetCoord();
            
			/*--- Compute sides of the triangle ---*/
			double a = 0, b = 0, c = 0;
			for (iDim = 0; iDim < nDim; iDim++) {
				a += (Coord_2[iDim]-Coord_j[iDim])*(Coord_2[iDim]-Coord_j[iDim]);
				b += (Coord_3[iDim]-Coord_j[iDim])*(Coord_3[iDim]-Coord_j[iDim]);
				c += (Coord_3[iDim]-Coord_2[iDim])*(Coord_3[iDim]-Coord_2[iDim]);
			}
			a = sqrt(a); b = sqrt(b); c = sqrt(c);
            
			/*--- Compute semiperimeter (s) and area ---*/
			double s = 0.5*(a + b + c);
			double Area = sqrt(s*(s-a)*(s-b)*(s-c));
            
			/*--- Compute radius of the circumcircle (R) and of the incircle (r) ---*/
			double R = (a*b*c) / (4.0*Area);
			double r = Area / s;
			double roR = r / R;
            
			/*--- Update statistics ---*/
			if (roR < statistics[0])
				statistics[0] = roR;
			statistics[1] += roR;
            
		}
	}
	statistics[1] /= this->GetnElem();
    
}

void CPhysicalGeometry::SetRotationalVelocity(CConfig *config) {
    
	unsigned long iPoint;
	double RotVel[3], Distance[3], *Coord, Center[3], Omega[3], L_Ref;
    
    int rank = MASTER_NODE;
    
    /*--- Center of rotation & angular velocity vector from config ---*/
    
    Center[0] = config->GetMotion_Origin_X(ZONE_0);
    Center[1] = config->GetMotion_Origin_Y(ZONE_0);
    Center[2] = config->GetMotion_Origin_Z(ZONE_0);
    Omega[0]  = config->GetRotation_Rate_X(ZONE_0)/config->GetOmega_Ref();
    Omega[1]  = config->GetRotation_Rate_Y(ZONE_0)/config->GetOmega_Ref();
    Omega[2]  = config->GetRotation_Rate_Z(ZONE_0)/config->GetOmega_Ref();
    L_Ref     = config->GetLength_Ref();
    
    /*--- Print some information to the console ---*/
    
    if (rank == MASTER_NODE) {
        cout << " Rotational origin (x,y,z): ( " << Center[0] << ", " << Center[1];
        cout << ", " << Center[2] << " )" << endl;
        cout << " Angular velocity about x, y, z axes: ( " << Omega[0] << ", ";
        cout << Omega[1] << ", " << Omega[2] << " ) rad/s" << endl;
    }
    
	/*--- Loop over all nodes and set the rotational velocity ---*/
    
	for (iPoint = 0; iPoint < nPoint; iPoint++) {
        
		/*--- Get the coordinates of the current node ---*/
        
		Coord = node[iPoint]->GetCoord();
        
		/*--- Calculate the non-dim. distance from the rotation center ---*/
        
		Distance[0] = (Coord[0]-Center[0])/L_Ref;
		Distance[1] = (Coord[1]-Center[1])/L_Ref;
		Distance[2] = (Coord[2]-Center[2])/L_Ref;
        
		/*--- Calculate the angular velocity as omega X r ---*/
        
		RotVel[0] = Omega[1]*(Distance[2]) - Omega[2]*(Distance[1]);
		RotVel[1] = Omega[2]*(Distance[0]) - Omega[0]*(Distance[2]);
		RotVel[2] = Omega[0]*(Distance[1]) - Omega[1]*(Distance[0]);
        
        /*--- Store the grid velocity at this node ---*/
        
		node[iPoint]->SetGridVel(RotVel);
        
	}
    
}

void CPhysicalGeometry::SetGridVelocity(CConfig *config, unsigned long iter) {
    
	/*--- Local variables ---*/
    
	double *Coord_nP1 = NULL, *Coord_n = NULL, *Coord_nM1 = NULL;
    double TimeStep, GridVel = 0.0;
	unsigned long iPoint;
	unsigned short iDim;
    
	/*--- Compute the velocity of each node in the volume mesh ---*/
    
	for (iPoint = 0; iPoint < GetnPoint(); iPoint++) {
        
		/*--- Coordinates of the current point at n+1, n, & n-1 time levels ---*/
        
		Coord_nM1 = node[iPoint]->GetCoord_n1();
		Coord_n   = node[iPoint]->GetCoord_n();
		Coord_nP1 = node[iPoint]->GetCoord();
        
		/*--- Unsteady time step ---*/
        
		TimeStep = config->GetDelta_UnstTimeND();
        
		/*--- Compute mesh velocity with 1st or 2nd-order approximation ---*/
        
		for(iDim = 0; iDim < nDim; iDim++) {
			if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
				GridVel = ( Coord_nP1[iDim] - Coord_n[iDim] ) / TimeStep;
			if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
				GridVel = ( 3.0*Coord_nP1[iDim] - 4.0*Coord_n[iDim]
                           + 1.0*Coord_nM1[iDim] ) / (2.0*TimeStep);
            
			/*--- Store grid velocity for this point ---*/
            
			node[iPoint]->SetGridVel(iDim, GridVel);
		}
	}
    
}

void CPhysicalGeometry::Set_MPI_GridVel(CConfig *config) {
	unsigned short iVar, iMarker, iPeriodic_Index;
	unsigned long iVertex, iPoint, nVert, nBuffer_Vector;
	double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi,
    psi, cosPsi, sinPsi, *newGridVel = NULL, *Buffer_Receive_GridVel = NULL, *GridVel;
	short SendRecv;
	int send_to, receive_from;
    
	newGridVel = new double[nDim];
    
	/*--- Send-Receive boundary conditions ---*/
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
		if (config->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) {
			SendRecv = config->GetMarker_All_SendRecv(iMarker);
			nVert = nVertex[iMarker];
			nBuffer_Vector = nVert*nDim;
			send_to = SendRecv-1;
			receive_from = abs(SendRecv)-1;
            
			/*--- Receive information  ---*/
			if (SendRecv < 0) {
                Buffer_Receive_GridVel = new double [nBuffer_Vector];
                
				/*--- Receive information without MPI ---*/
				for (iVertex = 0; iVertex < nVert; iVertex++) {
                    iPoint = vertex[iMarker][iVertex]->GetNode();
                    GridVel = node[iPoint]->GetGridVel();
                    for (iVar = 0; iVar < nDim; iVar++)
                        Buffer_Receive_GridVel[iVar*nVert+iVertex] = GridVel[iVar];
                }
                
				/*--- Do the coordinate transformation ---*/
				for (iVertex = 0; iVertex < nVert; iVertex++) {
                    
					/*--- Find point and its type of transformation ---*/
					iPoint = vertex[iMarker][iVertex]->GetNode();
					iPeriodic_Index = vertex[iMarker][iVertex]->GetRotation_Type();
                    
					/*--- Retrieve the supplied periodic information. ---*/
					angles = config->GetPeriodicRotation(iPeriodic_Index);
                    
					/*--- Store angles separately for clarity. ---*/
					theta    = angles[0];   phi    = angles[1]; psi    = angles[2];
					cosTheta = cos(theta);  cosPhi = cos(phi);  cosPsi = cos(psi);
					sinTheta = sin(theta);  sinPhi = sin(phi);  sinPsi = sin(psi);
                    
					/*--- Compute the rotation matrix. Note that the implicit
					 ordering is rotation about the x-axis, y-axis,
					 then z-axis. Note that this is the transpose of the matrix
					 used during the preprocessing stage. ---*/
					rotMatrix[0][0] = cosPhi*cosPsi; rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi; rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
					rotMatrix[0][1] = cosPhi*sinPsi; rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi; rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
					rotMatrix[0][2] = -sinPhi; rotMatrix[1][2] = sinTheta*cosPhi; rotMatrix[2][2] = cosTheta*cosPhi;
                    
					/*--- Copy conserved variables before performing transformation. ---*/
					for (iVar = 0; iVar < nDim; iVar++)
						newGridVel[iVar] = Buffer_Receive_GridVel[iVar*nVert+iVertex];
                    
					/*--- Rotate the momentum components. ---*/
					if (nDim == 2) {
						newGridVel[1] = rotMatrix[0][0]*Buffer_Receive_GridVel[1*nVert+iVertex] + rotMatrix[0][1]*Buffer_Receive_GridVel[2*nVert+iVertex];
						newGridVel[2] = rotMatrix[1][0]*Buffer_Receive_GridVel[1*nVert+iVertex] + rotMatrix[1][1]*Buffer_Receive_GridVel[2*nVert+iVertex];
					}
					else {
						newGridVel[1] = rotMatrix[0][0]*Buffer_Receive_GridVel[1*nVert+iVertex] + rotMatrix[0][1]*Buffer_Receive_GridVel[2*nVert+iVertex] + rotMatrix[0][2]*Buffer_Receive_GridVel[3*nVert+iVertex];
						newGridVel[2] = rotMatrix[1][0]*Buffer_Receive_GridVel[1*nVert+iVertex] + rotMatrix[1][1]*Buffer_Receive_GridVel[2*nVert+iVertex] + rotMatrix[1][2]*Buffer_Receive_GridVel[3*nVert+iVertex];
						newGridVel[3] = rotMatrix[2][0]*Buffer_Receive_GridVel[1*nVert+iVertex] + rotMatrix[2][1]*Buffer_Receive_GridVel[2*nVert+iVertex] + rotMatrix[2][2]*Buffer_Receive_GridVel[3*nVert+iVertex];
					}
                    
					/*--- Copy transformed conserved variables back into buffer. ---*/
					for (iVar = 0; iVar < nDim; iVar++)
						Buffer_Receive_GridVel[iVar*nVert+iVertex] = newGridVel[iVar];
                    
                    for (iVar = 0; iVar < nDim; iVar++)
                        node[iPoint]->SetGridVel(iVar, Buffer_Receive_GridVel[iVar*nVert+iVertex]);
                    
				}
                delete [] Buffer_Receive_GridVel;
			}
		}
	}
	delete [] newGridVel;
    
}

void CPhysicalGeometry::SetPeriodicBoundary(CConfig *config) {
    
	unsigned short iMarker, jMarker, kMarker = 0, iPeriodic, iDim, nPeriodic = 0, VTK_Type;
	unsigned long iNode, iIndex, iVertex, iPoint, iElem, kElem;
	unsigned long jElem, kPoint = 0, jVertex = 0, jPoint = 0, pPoint = 0, nPointPeriodic, newNodes[4];
	vector<unsigned long>::iterator IterElem, IterPoint[MAX_NUMBER_PERIODIC][2], IterNewElem[MAX_NUMBER_MARKER];
	double *center, *angles, rotMatrix[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}},
    translation[3], *trans, theta, phi, psi, cosTheta, sinTheta, cosPhi, sinPhi, cosPsi, sinPsi,
    dx, dy, dz, rotCoord[3], epsilon = 1e-10, mindist = 1e6, *Coord_i, *Coord_j, dist = 0.0;
	bool isBadMatch = false;
    
    /*--- It only create the mirror structure for the second boundary ---*/
    bool CreateMirror[10];
    CreateMirror[1] = false;
    CreateMirror[2] = true;
    
	/*--- Send an initial message to the console. ---*/
	cout << "Setting the periodic boundary conditions." <<endl;
    
	/*--- Loop through each marker to find any periodic boundaries. ---*/
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if (config->GetMarker_All_Boundary(iMarker) == PERIODIC_BOUNDARY) {
            
			/*--- Evaluate the number of periodic boundary conditions defined
             in the geometry file ---*/
			nPeriodic++;
            
			/*--- Get marker index of the periodic donor boundary. ---*/
			jMarker = config->GetMarker_Periodic_Donor(config->GetMarker_All_Tag(iMarker));
            
			/*--- Write some info to the console. ---*/
			cout << "Checking " << config->GetMarker_All_Tag(iMarker);
			cout << " boundary against periodic donor, " << config->GetMarker_All_Tag(jMarker) << ". ";
            
			/*--- Retrieve the supplied periodic information. ---*/
			center = config->GetPeriodicRotCenter(config->GetMarker_All_Tag(iMarker));
			angles = config->GetPeriodicRotAngles(config->GetMarker_All_Tag(iMarker));
			trans  = config->GetPeriodicTranslation(config->GetMarker_All_Tag(iMarker));
            
			/*--- Store (center+trans) as it is constant and will be added on. ---*/
			translation[0] = center[0] + trans[0];
			translation[1] = center[1] + trans[1];
			translation[2] = center[2] + trans[2];
            
			/*--- Store angles separately for clarity. Compute sines/cosines. ---*/
			theta = angles[0];
			phi   = angles[1];
			psi   = angles[2];
            
			cosTheta = cos(theta);  cosPhi = cos(phi);  cosPsi = cos(psi);
			sinTheta = sin(theta);  sinPhi = sin(phi);  sinPsi = sin(psi);
            
			/*--- Compute the rotation matrix. Note that the implicit
             ordering is rotation about the x-axis, y-axis, then z-axis. ---*/
			rotMatrix[0][0] = cosPhi*cosPsi;
			rotMatrix[1][0] = cosPhi*sinPsi;
			rotMatrix[2][0] = -sinPhi;
            
			rotMatrix[0][1] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
			rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
			rotMatrix[2][1] = sinTheta*cosPhi;
            
			rotMatrix[0][2] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
			rotMatrix[1][2] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
			rotMatrix[2][2] = cosTheta*cosPhi;
            
			/*--- Loop through all vertices and find/set the periodic point. ---*/
			for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
                
				/*--- Retrieve node information for this boundary point. ---*/
				iPoint  = vertex[iMarker][iVertex]->GetNode();
				Coord_i = node[iPoint]->GetCoord();
                
				/*--- Get the position vector from rot center to point. ---*/
				dx = Coord_i[0] - center[0];
				dy = Coord_i[1] - center[1];
				if (nDim == 3) {
					dz = Coord_i[2] - center[2];
				} else {
					dz = 0.0;
				}
                
				/*--- Compute transformed point coordinates. ---*/
				rotCoord[0] = rotMatrix[0][0]*dx
                + rotMatrix[0][1]*dy
                + rotMatrix[0][2]*dz + translation[0];
                
				rotCoord[1] = rotMatrix[1][0]*dx
                + rotMatrix[1][1]*dy
                + rotMatrix[1][2]*dz + translation[1];
                
				rotCoord[2] = rotMatrix[2][0]*dx
                + rotMatrix[2][1]*dy
                + rotMatrix[2][2]*dz + translation[2];
                
				/*--- Perform a search to find the closest donor point. ---*/
				mindist = 1e10;
				for (jVertex = 0; jVertex < nVertex[jMarker]; jVertex++) {
                    
					/*--- Retrieve information for this jPoint. ---*/
					jPoint = vertex[jMarker][jVertex]->GetNode();
					Coord_j = node[jPoint]->GetCoord();
                    
					/*--- Check the distance between the computed periodic
					 location and this jPoint. ---*/
					dist = 0.0;
					for (iDim = 0; iDim < nDim; iDim++){
						dist += (Coord_j[iDim]-rotCoord[iDim])*(Coord_j[iDim]-rotCoord[iDim]);
					}
					dist = sqrt(dist);
                    
					/*---  Store vertex information if this is the closest
					 point found thus far. ---*/
					if (dist < mindist) { mindist = dist; pPoint = jPoint; }
				}
                
				/*--- Set the periodic point for this iPoint. ---*/
				vertex[iMarker][iVertex]->SetDonorPoint(pPoint);
                
				/*--- Print warning if the nearest point was not within
                 the specified tolerance. Computation will continue. ---*/
				if (mindist > epsilon) {
					isBadMatch = true;
					cout.precision(10);
					cout << endl;
					cout << "   Bad match for point " << iPoint << ".\tNearest";
					cout << " donor distance: " << scientific << mindist << ".";
				}
			}
            
			/*--- Print final warning when finding bad matches. ---*/
			if (isBadMatch) {
				cout << endl;
				cout << "\n !!! Warning !!!" << endl;
				cout << "Bad matches found. Computation will continue, but be cautious.\n";
			}
			cout << endl;
			isBadMatch = false;
            
		}
    
	/*--- Create a vector to identify the points that belong to each periodic boundary condition ---*/
	bool *PeriodicBC = new bool [nPoint];
	for (iPoint = 0; iPoint < nPoint; iPoint++) PeriodicBC[iPoint] = false;
    
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if (config->GetMarker_All_Boundary(iMarker) == PERIODIC_BOUNDARY)
			for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
				iPoint = vertex[iMarker][iVertex]->GetNode();
				PeriodicBC[iPoint] = true;
			}
    
	/*--- Determine the new points that must be added to each periodic boundary,
     note that only one of the boundaries require the extra data ---*/
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
		if (config->GetMarker_All_Boundary(iMarker) == PERIODIC_BOUNDARY) {
			iPeriodic = config->GetMarker_All_PerBound(iMarker);
            
			/*--- An integer identify the periodic boundary condition --*/
            for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
                
                /*--- iPoint is the original point on the surface and jPoint is the
                 equivalent point in the other periodic surface ---*/
                iPoint = vertex[iMarker][iVertex]->GetNode();
                jPoint = vertex[iMarker][iVertex]->GetDonorPoint();
                
                /*--- First the case in which it is necessary to create a mirror set of elements ---*/
                if (CreateMirror[iPeriodic]) {
                    /*--- Now we must determine the neighbor points (including indirect ones) to the periodic points
                     and store all the information (in this list we do not include the points
                     that already belong to the periodic boundary), we also add the elements that
                     share a point with the periodic boundary condition ---*/
                    for (iIndex = 0; iIndex < node[jPoint]->GetnElem(); iIndex++) {
                        iElem = node[jPoint]->GetElem(iIndex);
                        PeriodicElem[iPeriodic].push_back(iElem);
                        for (unsigned short iNode = 0; iNode <	elem[iElem]->GetnNodes(); iNode ++) {
                            kPoint = elem[iElem]->GetNode(iNode);
                            if (!PeriodicBC[kPoint]) PeriodicPoint[iPeriodic][0].push_back(kPoint);
                        }
                    }
                }
                /*--- Second the case where no new element is added, neither points ---*/
                else {
                    PeriodicPoint[iPeriodic][0].push_back(jPoint);
                    PeriodicPoint[iPeriodic][1].push_back(iPoint);
                }
            }
		}
    }
    
	/*--- Sort the points that must be sended and delete repeated points ---*/
	for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++) {
        if (CreateMirror[iPeriodic]) {
            sort( PeriodicPoint[iPeriodic][0].begin(), PeriodicPoint[iPeriodic][0].end());
            IterPoint[iPeriodic][0] = unique( PeriodicPoint[iPeriodic][0].begin(), PeriodicPoint[iPeriodic][0].end());
            PeriodicPoint[iPeriodic][0].resize( IterPoint[iPeriodic][0] - PeriodicPoint[iPeriodic][0].begin() );
        }
	}
    
	/*--- Create a list of the points that receive the values (only the new points) ---*/
	nPointPeriodic = nPoint;
	for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++) {
        if (CreateMirror[iPeriodic]) {
            for (iPoint = 0; iPoint < PeriodicPoint[iPeriodic][0].size(); iPoint++) {
                PeriodicPoint[iPeriodic][1].push_back(nPointPeriodic);
                nPointPeriodic++;
            }
        }
    }
    
	/*--- Sort the elements that must be replicated in the periodic boundary
	 and delete the repeated elements ---*/
	for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++) {
        if (CreateMirror[iPeriodic]) {
            sort( PeriodicElem[iPeriodic].begin(), PeriodicElem[iPeriodic].end());
            IterElem = unique( PeriodicElem[iPeriodic].begin(), PeriodicElem[iPeriodic].end());
            PeriodicElem[iPeriodic].resize( IterElem - PeriodicElem[iPeriodic].begin() );
        }
	}
    
	/*--- Check all SEND points to see if they also lie on another boundary. ---*/
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
		for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
            
			/*--- iPoint is a node that lies on the current marker. ---*/
			iPoint = vertex[iMarker][iVertex]->GetNode();
            
			/*--- Search through SEND points to check for iPoint. ---*/
			for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++) {
                if (CreateMirror[iPeriodic]) {
                    
                    /*--- jPoint is the SEND point. ---*/
                    for (iElem = 0; iElem < PeriodicPoint[iPeriodic][0].size(); iElem++) {
                        jPoint = PeriodicPoint[iPeriodic][0][iElem];
                        
                        /*--- If the two match, then jPoint lies on this boundary.
                         However, we are concerned with the new points, so we
                         will store kPoint instead. ---*/
                        if (iPoint == jPoint) {
                            kPoint = PeriodicPoint[iPeriodic][1][iElem];
                            
                            /*--- We also want the type of boundary element that this point
                             was within, so that we know what type of element to add
                             built from the new points. ---*/
                            bool isJPoint, isPeriodic;
                            for(jElem = 0; jElem < nElem_Bound[iMarker]; jElem++) {
                                isJPoint = false; isPeriodic = false;
                                for (iNode = 0; iNode < bound[iMarker][jElem]->GetnNodes(); iNode++) {
                                    if (bound[iMarker][jElem]->GetNode(iNode) == jPoint) isJPoint = true;
                                    if (PeriodicBC[bound[iMarker][jElem]->GetNode(iNode)]) isPeriodic = true;
                                }
                                
                                /*--- If both points were found, store this element. ---*/
                                if (isJPoint && isPeriodic) {
                                    OldBoundaryElems[iMarker].push_back(jElem);
                                }
                                
                            }
                            
                        }
                    }
                }
            }
		}
	}
    
	/*--- Sort the elements that must be added and remove duplicates. ---*/
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
		sort( OldBoundaryElems[iMarker].begin(), OldBoundaryElems[iMarker].end());
		IterNewElem[iMarker] = unique( OldBoundaryElems[iMarker].begin(), OldBoundaryElems[iMarker].end());
		OldBoundaryElems[iMarker].resize( IterNewElem[iMarker] - OldBoundaryElems[iMarker].begin() );
	}
    
	/*--- Create the new boundary elements. Points making up these new
     elements must either be SEND/RECEIVE or periodic points. ---*/
	nNewElem_Bound = new unsigned long[nMarker];
	newBound = new CPrimalGrid**[nMarker];
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
		nNewElem_Bound[iMarker] = OldBoundaryElems[iMarker].size();
		newBound[iMarker]       = new CPrimalGrid*[nNewElem_Bound[iMarker]];
        
		/*--- Loop through all new elements to be added. ---*/
		for (iElem = 0; iElem < nNewElem_Bound[iMarker]; iElem++) {
			jElem = OldBoundaryElems[iMarker][iElem];
            
			/*--- Loop through all nodes of this element. ---*/
			for (iNode = 0; iNode < bound[iMarker][jElem]->GetnNodes(); iNode++) {
				pPoint = bound[iMarker][jElem]->GetNode(iNode);
                
				/*--- Check if this node is a send point. If so, the corresponding
                 receive point will be used in making the new boundary element. ---*/
				for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++) {
                    for (kElem = 0; kElem < PeriodicPoint[iPeriodic][0].size(); kElem++) {
                        if (pPoint == PeriodicPoint[iPeriodic][0][kElem]) newNodes[iNode] = PeriodicPoint[iPeriodic][1][kElem];
                    }
				}
                
				/*--- Check if this node is a periodic point. If so, the corresponding
                 periodic point will be used in making the new boundary element. ---*/
				if (PeriodicBC[pPoint]) {
                    
					/*--- Find the corresponding periodic point. ---*/
					for (jMarker = 0; jMarker < config->GetnMarker_All(); jMarker++) {
						if (config->GetMarker_All_Boundary(jMarker) == PERIODIC_BOUNDARY) {
							for (iVertex = 0; iVertex < nVertex[jMarker]; iVertex++) {
								if (pPoint == vertex[jMarker][iVertex]->GetNode()) {kMarker = jMarker; jVertex = iVertex;}
							}
						}
					}
					newNodes[iNode] = vertex[kMarker][jVertex]->GetDonorPoint();
				}
			}
            
			/*--- Now instantiate the new element. ---*/
			VTK_Type = bound[iMarker][jElem]->GetVTK_Type();
			switch(VTK_Type) {
                case LINE:
                    newBound[iMarker][iElem] = new CLine(newNodes[0],newNodes[1],2);
                    break;
                case TRIANGLE:
                    newBound[iMarker][iElem] = new CTriangle(newNodes[0],newNodes[1],newNodes[2],3);
                    break;
                case RECTANGLE:
                    newBound[iMarker][iElem] = new CRectangle(newNodes[0],newNodes[1],newNodes[2],newNodes[3],3);
                    break;
			}
            
		}
	}
    
	delete[] PeriodicBC;
    
}

void CPhysicalGeometry::ComputeSurf_Curvature(CConfig *config) {
	unsigned short iMarker, iNeigh_Point, iDim, iNode, iNeighbor_Nodes, Neighbor_Node;
	unsigned long Neighbor_Point, iVertex, iPoint, jPoint,iElem_Bound, iEdge, nLocalVertex, nGlobalVertex , MaxLocalVertex , *Buffer_Send_nVertex, *Buffer_Receive_nVertex, nBuffer, TotalnPointDomain;
    int iProcessor, nProcessor;
    vector<unsigned long> Point_NeighborList, Elem_NeighborList, Point_Triangle, Point_Edge, Point_Critical;
    vector<unsigned long>::iterator it;
    double U[3], V[3], W[3], Length_U, Length_V, Length_W, CosValue, Angle_Value, *K, *Angle_Defect, *Area_Vertex, *Angle_Alpha, *Angle_Beta, **NormalMeanK, MeanK, GaussK, MaxPrinK, MinPrinK, cot_alpha, cot_beta, delta, X1, X2, X3, Y1, Y2, Y3, radius, *Buffer_Send_Coord, *Buffer_Receive_Coord, *Coord, Dist, MinDist, MaxK, MinK, SigmaK;
    bool *Check_Edge;
    
	int rank = MASTER_NODE;
    
    /*--- Allocate surface curvature ---*/
    K = new double [nPoint];
    
    
	if (nDim == 2) {
        
		/*--- Loop over all the markers ---*/
		for (iMarker = 0; iMarker < nMarker; iMarker++) {
            
            if (config->GetMarker_All_Boundary(iMarker) != SEND_RECEIVE) {
                
                /*--- Loop through all marker vertices again, this time also
                 finding the neighbors of each node.---*/
                for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
                    iPoint  = vertex[iMarker][iVertex]->GetNode();
                    
                    if (node[iPoint]->GetDomain()) {
                        
                        /*--- Loop through neighbors. In 2-D, there should be 2 nodes on either
                         side of this vertex that lie on the same surface. ---*/
                        Point_Edge.clear();
                        
                        for (iNeigh_Point = 0; iNeigh_Point < node[iPoint]->GetnPoint(); iNeigh_Point++) {
                            Neighbor_Point = node[iPoint]->GetPoint(iNeigh_Point);
                            
                            /*--- Check if this neighbor lies on the surface. If so,
                             add to the list of neighbors. ---*/
                            if (node[Neighbor_Point]->GetPhysicalBoundary()) {
                                Point_Edge.push_back(Neighbor_Point);
                            }
                            
                        }
                        
                        if (Point_Edge.size() == 2) {
                            
                            /*--- Compute the curvature using three points ---*/
                            X1 = node[iPoint]->GetCoord(0);
                            X2 = node[Point_Edge[0]]->GetCoord(0);
                            X3 = node[Point_Edge[1]]->GetCoord(0);
                            Y1 = node[iPoint]->GetCoord(1);
                            Y2 = node[Point_Edge[0]]->GetCoord(1);
                            Y3 = node[Point_Edge[1]]->GetCoord(1);
                            
                            radius = sqrt(((X2-X1)*(X2-X1) + (Y2-Y1)*(Y2-Y1))*
                                          ((X2-X3)*(X2-X3) + (Y2-Y3)*(Y2-Y3))*
                                          ((X3-X1)*(X3-X1) + (Y3-Y1)*(Y3-Y1)))/
                            (2.0*fabs(X1*Y2+X2*Y3+X3*Y1-X1*Y3-X2*Y1-X3*Y2));
                            
                            K[iPoint] = 1.0/radius;
                            
                        }
                        
                    }
                    
                }
                
            }
            
		}
        
	}
    
    else {
        
        Angle_Defect = new double [nPoint];
        Area_Vertex = new double [nPoint];
        for (iPoint = 0; iPoint < nPoint; iPoint++) {
            Angle_Defect[iPoint] = 2*PI_NUMBER;
            Area_Vertex[iPoint] = 0.0;
        }
        
        Angle_Alpha = new double [nEdge];
        Angle_Beta = new double [nEdge];
        Check_Edge = new bool [nEdge];
        for (iEdge = 0; iEdge < nEdge; iEdge++) {
            Angle_Alpha[iEdge] = 0.0;
            Angle_Beta[iEdge] = 0.0;
            Check_Edge[iEdge] = true;
        }
        
        NormalMeanK = new double *[nPoint];
        for (iPoint = 0; iPoint < nPoint; iPoint++) {
            NormalMeanK[iPoint] = new double [nDim];
            for (iDim = 0; iDim < nDim; iDim++) {
                NormalMeanK[iPoint][iDim] = 0.0;
            }
        }
        
        /*--- Loop over all the markers ---*/
        for (iMarker = 0; iMarker < nMarker; iMarker++) {
            
            if (config->GetMarker_All_Boundary(iMarker) != SEND_RECEIVE) {
                
                /*--- Loop over all the boundary elements ---*/
                for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
                    
                    /*--- Only triangles ---*/
                    if (bound[iMarker][iElem_Bound]->GetVTK_Type() == TRIANGLE) {
                        
                        /*--- Loop over all the nodes of the boundary element ---*/
                        for(iNode = 0; iNode < bound[iMarker][iElem_Bound]->GetnNodes(); iNode++) {
                            
                            iPoint = bound[iMarker][iElem_Bound]->GetNode(iNode);
                            
                            Point_Triangle.clear();
                            
                            for(iNeighbor_Nodes = 0; iNeighbor_Nodes < bound[iMarker][iElem_Bound]->GetnNeighbor_Nodes(iNode); iNeighbor_Nodes++) {
                                Neighbor_Node = bound[iMarker][iElem_Bound]->GetNeighbor_Nodes(iNode, iNeighbor_Nodes);
                                Neighbor_Point = bound[iMarker][iElem_Bound]->GetNode(Neighbor_Node);
                                Point_Triangle.push_back(Neighbor_Point);
                            }
                            
                            iEdge = FindEdge(Point_Triangle[0], Point_Triangle[1]);
                            
                            for (iDim = 0; iDim < nDim; iDim++) {
                                U[iDim] = node[Point_Triangle[0]]->GetCoord(iDim) - node[iPoint]->GetCoord(iDim);
                                V[iDim] = node[Point_Triangle[1]]->GetCoord(iDim) - node[iPoint]->GetCoord(iDim);
                            }
                            
                            W[0] = 0.5*(U[1]*V[2]-U[2]*V[1]); W[1] = -0.5*(U[0]*V[2]-U[2]*V[0]); W[2] = 0.5*(U[0]*V[1]-U[1]*V[0]);
                            
                            Length_U = 0.0, Length_V = 0.0, Length_W = 0.0, CosValue = 0.0;
                            for (iDim = 0; iDim < nDim; iDim++) { Length_U += U[iDim]*U[iDim]; Length_V += V[iDim]*V[iDim]; Length_W += W[iDim]*W[iDim]; }
                            Length_U = sqrt(Length_U); Length_V = sqrt(Length_V); Length_W = sqrt(Length_W);
                            for (iDim = 0; iDim < nDim; iDim++) { U[iDim] /= Length_U; V[iDim] /= Length_V; CosValue += U[iDim]*V[iDim]; }
                            if (CosValue >= 1.0) CosValue = 1.0;
                            if (CosValue <= -1.0) CosValue = -1.0;
                            
                            Angle_Value = acos(CosValue);
                            Area_Vertex[iPoint] += Length_W;
                            Angle_Defect[iPoint] -= Angle_Value;
                            if (Angle_Alpha[iEdge] == 0.0) Angle_Alpha[iEdge] = Angle_Value;
                            else Angle_Beta[iEdge] = Angle_Value;
                            
                        }
                    }
                }
            }
        }
        
        /*--- Compute mean curvature ---*/
        for (iMarker = 0; iMarker < nMarker; iMarker++) {
            if (config->GetMarker_All_Boundary(iMarker) != SEND_RECEIVE) {
                for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
                    if (bound[iMarker][iElem_Bound]->GetVTK_Type() == TRIANGLE) {
                        for(iNode = 0; iNode < bound[iMarker][iElem_Bound]->GetnNodes(); iNode++) {
                            iPoint = bound[iMarker][iElem_Bound]->GetNode(iNode);
                            
                            for(iNeighbor_Nodes = 0; iNeighbor_Nodes < bound[iMarker][iElem_Bound]->GetnNeighbor_Nodes(iNode); iNeighbor_Nodes++) {
                                Neighbor_Node = bound[iMarker][iElem_Bound]->GetNeighbor_Nodes(iNode, iNeighbor_Nodes);
                                jPoint = bound[iMarker][iElem_Bound]->GetNode(Neighbor_Node);
                                
                                iEdge = FindEdge(iPoint, jPoint);
                                
                                if (Check_Edge[iEdge]) {
                                    
                                    Check_Edge[iEdge] = false;
                                    
                                    cot_alpha = 1.0/tan(Angle_Alpha[iEdge]);
                                    cot_beta = 1.0/tan(Angle_Beta[iEdge]);
                                    
                                    /*--- iPoint, and jPoint ---*/
                                    for (iDim = 0; iDim < nDim; iDim++) {
                                        NormalMeanK[iPoint][iDim] += 3.0 * (cot_alpha + cot_beta) * (node[iPoint]->GetCoord(iDim) - node[jPoint]->GetCoord(iDim)) / Area_Vertex[iPoint];
                                        NormalMeanK[jPoint][iDim] += 3.0 * (cot_alpha + cot_beta) * (node[jPoint]->GetCoord(iDim) - node[iPoint]->GetCoord(iDim)) / Area_Vertex[jPoint];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        
        /*--- Compute Gauss, mean, max and min principal curvature,
         and set the list of critical points ---*/
        
        for (iMarker = 0; iMarker < nMarker; iMarker++) {
            if (config->GetMarker_All_Boundary(iMarker) != SEND_RECEIVE) {
                for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
                    iPoint  = vertex[iMarker][iVertex]->GetNode();
                    
                    if (node[iPoint]->GetDomain()) {
                        
                        if (Area_Vertex[iPoint] != 0.0) GaussK = 3.0*Angle_Defect[iPoint]/Area_Vertex[iPoint];
                        else GaussK = 0.0;
                        
                        MeanK = 0.0;
                        for (iDim = 0; iDim < nDim; iDim++)
                            MeanK += NormalMeanK[iPoint][iDim]*NormalMeanK[iPoint][iDim];
                        MeanK = sqrt(MeanK);
                        
                        delta = max((MeanK*MeanK - GaussK), 0.0);
                        
                        MaxPrinK = MeanK + sqrt(delta);
                        MinPrinK = MeanK - sqrt(delta);
                        
                        /*--- Store the curvature value ---*/
                        K[iPoint] = MaxPrinK;
                        
                    }
                    
                }
            }
        }
        
        delete [] Angle_Defect;
        delete [] Area_Vertex;
        delete [] Angle_Alpha;
        delete [] Angle_Beta;
        delete [] Check_Edge;
        
        for (iPoint = 0; iPoint < nPoint; iPoint++)
            delete NormalMeanK[iPoint];
        delete [] NormalMeanK;
        
    }
    
    /*--- Sharp edge detection is based in the statistical
     distribution of the curvature ---*/
    
    MaxK = K[0]; MinK = K[0]; MeanK = 0.0; TotalnPointDomain = 0;
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
        if (config->GetMarker_All_Boundary(iMarker) != SEND_RECEIVE) {
            for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
                iPoint  = vertex[iMarker][iVertex]->GetNode();
                if (node[iPoint]->GetDomain()) {
                    MaxK = max(MaxK, fabs(K[iPoint]));
                    MinK = min(MinK, fabs(K[iPoint]));
                    MeanK += fabs(K[iPoint]);
                    TotalnPointDomain++;
                }
            }
        }
    }
    
    /*--- Compute the mean ---*/
    MeanK /= double(TotalnPointDomain);
    
    /*--- Compute the standard deviation ---*/
    SigmaK = 0.0;
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
        if (config->GetMarker_All_Boundary(iMarker) != SEND_RECEIVE) {
            for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
                iPoint  = vertex[iMarker][iVertex]->GetNode();
                if (node[iPoint]->GetDomain()) {
                    SigmaK += (fabs(K[iPoint]) - MeanK) * (fabs(K[iPoint]) - MeanK);
                }
            }
        }
    }
    
    SigmaK = sqrt(SigmaK/double(TotalnPointDomain));
    
    if (rank == MASTER_NODE)
        cout << "Max K: " << MaxK << ". Mean K: " << MeanK << ". Standard deviation K: " << SigmaK << "." <<endl;
    
    Point_Critical.clear();
    
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
        if (config->GetMarker_All_Boundary(iMarker) != SEND_RECEIVE) {
            for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
                iPoint  = vertex[iMarker][iVertex]->GetNode();
                if (node[iPoint]->GetDomain()) {
                    if (fabs(K[iPoint]) > MeanK + config->GetRefSharpEdges()*SigmaK) {
                        Point_Critical.push_back(iPoint);
                    }
                }
            }
        }
    }
    
    /*--- Variables and buffers needed for MPI ---*/
    
	nProcessor = 1;
    
	nLocalVertex = 0, nGlobalVertex = 0, MaxLocalVertex = 0;
	Buffer_Send_nVertex    = new unsigned long [1];
	Buffer_Receive_nVertex = new unsigned long [nProcessor];
    
    /*--- Count the total number of critical edge nodes. ---*/
    nLocalVertex = Point_Critical.size();
    Buffer_Send_nVertex[0] = nLocalVertex;
    
    /*--- Communicate to all processors the total number of critical edge nodes. ---*/
    MaxLocalVertex = nLocalVertex;
    nGlobalVertex = nLocalVertex;
    Buffer_Receive_nVertex[0] = nLocalVertex;
    
    
    /*--- Create and initialize to zero some buffers to hold the coordinates
     of the boundary nodes that are communicated from each partition (all-to-all). ---*/
    
	Buffer_Send_Coord     = new double [MaxLocalVertex*nDim];
	Buffer_Receive_Coord  = new double [nProcessor*MaxLocalVertex*nDim];
    nBuffer               = MaxLocalVertex*nDim;
    
	for (iVertex = 0; iVertex < MaxLocalVertex; iVertex++) {
		for (iDim = 0; iDim < nDim; iDim++) {
			Buffer_Send_Coord[iVertex*nDim+iDim] = 0.0;
        }
    }
    
    /*--- Retrieve and store the coordinates of the sharp edges boundary nodes on
     the local partition and broadcast them to all partitions. ---*/
    
    for (iVertex = 0; iVertex < Point_Critical.size(); iVertex++) {
        iPoint = Point_Critical[iVertex];
        for (iDim = 0; iDim < nDim; iDim++)
            Buffer_Send_Coord[iVertex*nDim+iDim] = node[iPoint]->GetCoord(iDim);
    }
    
    for (iVertex = 0; iVertex < Point_Critical.size(); iVertex++) {
        for (iDim = 0; iDim < nDim; iDim++) {
            Buffer_Receive_Coord[iVertex*nDim+iDim] = Buffer_Send_Coord[iVertex*nDim+iDim];
        }
    }
    
    /*--- Loop over all interior mesh nodes on the local partition and compute
     the distances to each of the no-slip boundary nodes in the entire mesh.
     Store the minimum distance to the wall for each interior mesh node. ---*/
    
	for (iPoint = 0; iPoint < GetnPoint(); iPoint++) {
		Coord = node[iPoint]->GetCoord();
        
        MinDist = 1E20;
		for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
			for (iVertex = 0; iVertex < Buffer_Receive_nVertex[iProcessor]; iVertex++) {
				Dist = 0.0;
				for (iDim = 0; iDim < nDim; iDim++) {
					Dist += (Coord[iDim]-Buffer_Receive_Coord[(iProcessor*MaxLocalVertex+iVertex)*nDim+iDim])*
					(Coord[iDim]-Buffer_Receive_Coord[(iProcessor*MaxLocalVertex+iVertex)*nDim+iDim]);
                }
                Dist = sqrt(Dist);
				if (Dist < MinDist) MinDist = Dist;
			}
        }
        node[iPoint]->SetSharpEdge_Distance(MinDist);
	}
    
    /*--- Deallocate Max curvature ---*/
    delete[] K;
    
    /*--- Deallocate the buffers needed for the MPI communication. ---*/
	delete[] Buffer_Send_Coord;
	delete[] Buffer_Receive_Coord;
	delete[] Buffer_Send_nVertex;
	delete[] Buffer_Receive_nVertex;
    
}

void CPhysicalGeometry::FindNormal_Neighbor(CConfig *config) {
    double cos_max, scalar_prod, norm_vect, norm_Normal, cos_alpha, diff_coord, *Normal;
    unsigned long Point_Normal, jPoint;
    unsigned short iNeigh, iMarker, iDim;
	unsigned long iPoint, iVertex;
    
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
        
		if (config->GetMarker_All_Boundary(iMarker) != SEND_RECEIVE &&
            config->GetMarker_All_Boundary(iMarker) != INTERFACE_BOUNDARY &&
            config->GetMarker_All_Boundary(iMarker) != NEARFIELD_BOUNDARY ) {
			
			for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
                
				iPoint = vertex[iMarker][iVertex]->GetNode();
                Normal = vertex[iMarker][iVertex]->GetNormal();
                
                /*--- Compute closest normal neighbor, note that the normal are oriented inwards ---*/
                Point_Normal = 0; cos_max = -1.0;
                for (iNeigh = 0; iNeigh < node[iPoint]->GetnPoint(); iNeigh++) {
                    jPoint = node[iPoint]->GetPoint(iNeigh);
                    scalar_prod = 0.0; norm_vect = 0.0; norm_Normal = 0.0;
                    for(iDim = 0; iDim < nDim; iDim++) {
                        diff_coord = node[jPoint]->GetCoord(iDim)-node[iPoint]->GetCoord(iDim);
                        scalar_prod += diff_coord*Normal[iDim];
                        norm_vect += diff_coord*diff_coord;
                        norm_Normal += Normal[iDim]*Normal[iDim];
                    }
                    norm_vect = sqrt(norm_vect);
                    norm_Normal = sqrt(norm_Normal);
                    cos_alpha = scalar_prod/(norm_vect*norm_Normal);
                    
                    /*--- Get maximum cosine ---*/
                    if (cos_alpha >= cos_max) {
                        Point_Normal = jPoint;
                        cos_max = cos_alpha;
                    }
                }
                vertex[iMarker][iVertex]->SetNormal_Neighbor(Point_Normal);
			}
        }
    }
}

void CPhysicalGeometry::SetGeometryPlanes(CConfig *config) {
    
	bool loop_on;
	unsigned short iMarker = 0;
	double auxXCoord, auxYCoord, auxZCoord,	*Face_Normal = NULL, auxArea, *Xcoord = NULL, *Ycoord = NULL, *Zcoord = NULL, *FaceArea = NULL;
	unsigned long jVertex, iVertex,ixCoord, iPoint, iVertex_Wall, nVertex_Wall = 0;
    
	/*--- Compute the total number of points on the near-field ---*/
	nVertex_Wall = 0;
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if ((config->GetMarker_All_Boundary(iMarker) == HEAT_FLUX) ||
            (config->GetMarker_All_Boundary(iMarker) == ISOTHERMAL) ||
            (config->GetMarker_All_Boundary(iMarker) == EULER_WALL))
			nVertex_Wall += nVertex[iMarker];
    
    
	/*--- Create an array with all the coordinates, points, pressures, face area,
	 equivalent area, and nearfield weight ---*/
	Xcoord = new double[nVertex_Wall];
	Ycoord = new double[nVertex_Wall];
	if (nDim == 3)	Zcoord = new double[nVertex_Wall];
	FaceArea = new double[nVertex_Wall];
    
	/*--- Copy the boundary information to an array ---*/
	iVertex_Wall = 0;
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if ((config->GetMarker_All_Boundary(iMarker) == HEAT_FLUX) ||
            (config->GetMarker_All_Boundary(iMarker) == ISOTHERMAL) ||
            (config->GetMarker_All_Boundary(iMarker) == EULER_WALL))
			for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
				iPoint = vertex[iMarker][iVertex]->GetNode();
				Xcoord[iVertex_Wall] = node[iPoint]->GetCoord(0);
				Ycoord[iVertex_Wall] = node[iPoint]->GetCoord(1);
				if (nDim==3) Zcoord[iVertex_Wall] = node[iPoint]->GetCoord(2);
				Face_Normal = vertex[iMarker][iVertex]->GetNormal();
				FaceArea[iVertex_Wall] = fabs(Face_Normal[nDim-1]);
				iVertex_Wall ++;
			}
    
    
	//vector<double> XCoordList;
	vector<double>::iterator IterXCoordList;
    
	for (iVertex = 0; iVertex < nVertex_Wall; iVertex++)
		XCoordList.push_back(Xcoord[iVertex]);
    
	sort( XCoordList.begin(), XCoordList.end());
	IterXCoordList = unique( XCoordList.begin(), XCoordList.end());
	XCoordList.resize( IterXCoordList - XCoordList.begin() );
    
	/*--- Create vectors and distribute the values among the different PhiAngle queues ---*/
	Xcoord_plane.resize(XCoordList.size());
	Ycoord_plane.resize(XCoordList.size());
	if (nDim==3) Zcoord_plane.resize(XCoordList.size());
	FaceArea_plane.resize(XCoordList.size());
	Plane_points.resize(XCoordList.size());
    
    
	double dist_ratio;
	unsigned long iCoord;
    
	/*--- Distribute the values among the different PhiAngles ---*/
	for (iPoint = 0; iPoint < nPoint; iPoint++) {
		if (node[iPoint]->GetDomain()) {
			loop_on = true;
			for (ixCoord = 0; ixCoord < XCoordList.size()-1 && loop_on; ixCoord++) {
				dist_ratio = (node[iPoint]->GetCoord(0) - XCoordList[ixCoord])/(XCoordList[ixCoord+1]- XCoordList[ixCoord]);
				if (dist_ratio >= 0 && dist_ratio <= 1.0) {
					if (dist_ratio <= 0.5) iCoord = ixCoord;
					else iCoord = ixCoord+1;
					Xcoord_plane[iCoord].push_back(node[iPoint]->GetCoord(0) );
					Ycoord_plane[iCoord].push_back(node[iPoint]->GetCoord(1) );
					if (nDim==3) Zcoord_plane[iCoord].push_back(node[iPoint]->GetCoord(2) );
					FaceArea_plane[iCoord].push_back(node[iPoint]->GetVolume());   ///// CHECK AREA CALCULATION
					Plane_points[iCoord].push_back(iPoint );
					loop_on = false;
				}
			}
		}
	}
    
	unsigned long auxPoint;
	/*--- Order the arrays in ascending values of y ---*/
	for (ixCoord = 0; ixCoord < XCoordList.size(); ixCoord++)
		for (iVertex = 0; iVertex < Xcoord_plane[ixCoord].size(); iVertex++)
			for (jVertex = 0; jVertex < Xcoord_plane[ixCoord].size() - 1 - iVertex; jVertex++)
				if (Ycoord_plane[ixCoord][jVertex] > Ycoord_plane[ixCoord][jVertex+1]) {
					auxXCoord = Xcoord_plane[ixCoord][jVertex]; Xcoord_plane[ixCoord][jVertex] = Xcoord_plane[ixCoord][jVertex+1]; Xcoord_plane[ixCoord][jVertex+1] = auxXCoord;
					auxYCoord = Ycoord_plane[ixCoord][jVertex]; Ycoord_plane[ixCoord][jVertex] = Ycoord_plane[ixCoord][jVertex+1]; Ycoord_plane[ixCoord][jVertex+1] = auxYCoord;
					auxPoint = Plane_points[ixCoord][jVertex]; Plane_points[ixCoord][jVertex] = Plane_points[ixCoord][jVertex+1]; Plane_points[ixCoord][jVertex+1] = auxPoint;
					if (nDim==3) {
						auxZCoord = Zcoord_plane[ixCoord][jVertex]; Zcoord_plane[ixCoord][jVertex] = Zcoord_plane[ixCoord][jVertex+1]; Zcoord_plane[ixCoord][jVertex+1] = auxZCoord;
					}
					auxArea = FaceArea_plane[ixCoord][jVertex]; FaceArea_plane[ixCoord][jVertex] = FaceArea_plane[ixCoord][jVertex+1]; FaceArea_plane[ixCoord][jVertex+1] = auxArea;
				}
    
	/*--- Delete structures ---*/
	delete[] Xcoord; delete[] Ycoord;
	if (nDim==3) delete[] Zcoord;
	delete[] FaceArea;
}

void CPhysicalGeometry::ComputeAirfoil_Section(double *Plane_P0, double *Plane_Normal, unsigned short iSection, CConfig *config,
                                               vector<double> &Xcoord_Airfoil, vector<double> &Ycoord_Airfoil,
                                               vector<double> &Zcoord_Airfoil, vector<unsigned long> &point1_Airfoil,
                                               vector<unsigned long> &point2_Airfoil, bool original_surface) {
    unsigned short iMarker, iNode, jNode, iDim, intersect;
    long MinDist_Point, MinDistAngle_Point;
    unsigned long iPoint, jPoint, iElem, Trailing_Point, Airfoil_Point, iVertex, jVertex, n;
    double Segment_P0[3] = {0.0, 0.0, 0.0}, Segment_P1[3] = {0.0, 0.0, 0.0}, Intersection[3] = {0.0, 0.0, 0.0}, Trailing_Coord, MinDist_Value, MinDistAngle_Value, Dist_Value,
    Airfoil_Tangent[3] = {0.0, 0.0, 0.0}, Segment[3] = {0.0, 0.0, 0.0}, Length, Angle_Value, Normal[3], Tangent[3], BiNormal[3], auxXCoord,
    auxYCoord, auxZCoord, auxpoint1, auxpoint2, zp1, zpn, Camber_Line, MaxAngle = 15, CosValue;
    vector<double> Xcoord, Ycoord, Zcoord, Z2coord, Xcoord_Normal, Ycoord_Normal, Zcoord_Normal, Xcoord_Camber, Ycoord_Camber, Zcoord_Camber;
    vector<unsigned long> Duplicate, point1, point2;
    vector<unsigned long>::iterator it;
    int rank = MASTER_NODE;
    double **Coord_Variation;
    bool Boundary, Monitoring;
    
  	// clean the vector, just to be safe
  	Xcoord_Airfoil.clear();
	Ycoord_Airfoil.clear();
	Zcoord_Airfoil.clear();
	point1_Airfoil.clear();
	point2_Airfoil.clear();
    
    /*--- Set the right plane in 2D (note the change in Y-Z plane) ---*/
    if (nDim == 2) {
        iSection = 0;
        Plane_P0[0] = 0.0;      Plane_P0[1] = 0.0;      Plane_P0[2] = 0.0;
        Plane_Normal[0] = 0.0;  Plane_Normal[1] = 1.0;  Plane_Normal[2] = 0.0;
    }
    
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
        for (iMarker = 0; iMarker < nMarker; iMarker++) {
            Boundary   = config->GetMarker_All_Boundary(iMarker);
            Monitoring = config->GetMarker_All_Monitoring(iMarker);
            
            if ((Boundary == EULER_WALL) || (Boundary == HEAT_FLUX) ||
                (Boundary == ISOTHERMAL) || (Boundary == NEARFIELD_BOUNDARY)) {
                
                if (Monitoring) {
                    
                    for (iElem = 0; iElem < nElem_Bound[iMarker]; iElem++) {
                        for(iNode = 0; iNode < bound[iMarker][iElem]->GetnNodes(); iNode++) {
                            iPoint = bound[iMarker][iElem]->GetNode(iNode);
                            for(jNode = 0; jNode < bound[iMarker][iElem]->GetnNodes(); jNode++) {
                                jPoint = bound[iMarker][iElem]->GetNode(jNode);
                                
                                if (jPoint > iPoint) {
                                    
                                    Segment_P0[0] = 0.0;
                                    Segment_P0[1] = 0.0;
                                    Segment_P0[2] = 0.0;
                                    
                                    Segment_P1[0] = 0.0;
                                    Segment_P1[1] = 0.0;
                                    Segment_P1[2] = 0.0;
                                    
                                    for (iDim = 0; iDim < nDim; iDim++) {
                                        Segment_P0[iDim] = node[iPoint]->GetCoord(iDim);
                                        Segment_P1[iDim] = node[jPoint]->GetCoord(iDim);
                                        // get pressure as well
                                    }
                                    
                                    /*--- In 2D add the points directly (note the change between Y and Z coordinate) ---*/
                                    if (nDim == 2) {
                                        Xcoord.push_back(Segment_P0[0]);  Xcoord.push_back(Segment_P1[0]);
                                        Ycoord.push_back(Segment_P0[2]);  Ycoord.push_back(Segment_P1[2]);
                                        Zcoord.push_back(Segment_P0[1]);  Zcoord.push_back(Segment_P1[1]);
                                    }
                                    /*--- In 3D compute the intersection ---*/
                                    else if (nDim == 3) {
                                        intersect = ComputeSegmentPlane_Intersection(Segment_P0, Segment_P1, Plane_P0, Plane_Normal, Intersection);
                                        if (intersect == 1) {
                                            Xcoord.push_back(Intersection[0]);
                                            Ycoord.push_back(Intersection[1]);
                                            Zcoord.push_back(Intersection[2]);
                                            point1.push_back(iPoint);			// let point1 be i
                                            point2.push_back(jPoint);			// let point2 be j
                                            //interpolate the pressures at the endpoints.
                                        }
                                    }
                                    
                                }
                                
                            }
                        }
                    }
                }
            }
        }
    }
    
    if (rank == MASTER_NODE) {
        
        /*--- Create a list with the duplicated points ---*/
        for (iVertex = 0; iVertex < Xcoord.size()-1; iVertex++) {
            for (jVertex = iVertex+1; jVertex < Xcoord.size(); jVertex++) {
                Segment[0] = Xcoord[jVertex] - Xcoord[iVertex];
                Segment[1] = Ycoord[jVertex] - Ycoord[iVertex];
                Segment[2] = Zcoord[jVertex] - Zcoord[iVertex];
                Dist_Value = sqrt(pow(Segment[0], 2.0) + pow(Segment[1], 2.0) + pow(Segment[2], 2.0));
                if (Dist_Value < 1E-6) {
                    Duplicate.push_back (jVertex);
                }
            }
        }
        
        sort(Duplicate.begin(), Duplicate.end());
        it = unique(Duplicate.begin(), Duplicate.end());
        Duplicate.resize(it - Duplicate.begin());
        
        /*--- Remove duplicated points (starting from the back) ---*/
        for (iVertex = Duplicate.size(); iVertex > 0; iVertex--) {
            Xcoord.erase (Xcoord.begin() + Duplicate[iVertex-1]);
            Ycoord.erase (Ycoord.begin() + Duplicate[iVertex-1]);
            Zcoord.erase (Zcoord.begin() + Duplicate[iVertex-1]);
            //also erase the points from i and j
            point1.erase(point1.begin() + Duplicate.at(iVertex-1));
            point2.erase(point2.begin() + Duplicate.at(iVertex-1));
        }
        
        /*--- Find the trailing edge ---*/
        Trailing_Point = 0; Trailing_Coord = Xcoord[0];
        for (iVertex = 1; iVertex < Xcoord.size(); iVertex++) {
            if (Xcoord[iVertex] > Trailing_Coord) {
                Trailing_Point = iVertex; Trailing_Coord = Xcoord[iVertex];
            }
        }
        
        /*--- Add the trailing edge to the list, and remove from the original list ---*/
        Xcoord_Airfoil.push_back(Xcoord[Trailing_Point]);
        Ycoord_Airfoil.push_back(Ycoord[Trailing_Point]);
        Zcoord_Airfoil.push_back(Zcoord[Trailing_Point]);
        
        Xcoord.erase (Xcoord.begin() + Trailing_Point);
        Ycoord.erase (Ycoord.begin() + Trailing_Point);
        Zcoord.erase (Zcoord.begin() + Trailing_Point);
        
        point1_Airfoil.push_back(point1.at(Trailing_Point));
        point2_Airfoil.push_back(point2.at(Trailing_Point));
        
        point1.erase(point1.begin() + Trailing_Point);
        point2.erase(point2.begin() + Trailing_Point);
        
        /*--- Find the next point using the right hand side rule ---*/
        MinDist_Value = 1E6;
        for (iVertex = 0; iVertex < Xcoord.size(); iVertex++) {
            Segment[0] = Xcoord[iVertex] - Xcoord_Airfoil[0];
            Segment[1] = Ycoord[iVertex] - Ycoord_Airfoil[0];
            Segment[2] = Zcoord[iVertex] - Zcoord_Airfoil[0];
            Dist_Value = sqrt(pow(Segment[0], 2.0) + pow(Segment[1], 2.0) + pow(Segment[2], 2.0));
            Segment[0] /= Dist_Value; Segment[1] /= Dist_Value; Segment[2] /= Dist_Value;
            
            if ((Dist_Value < MinDist_Value) && (Segment[2] > 0.0)) { MinDist_Point = iVertex; MinDist_Value = Dist_Value; }
        }
        
        Xcoord_Airfoil.push_back(Xcoord[MinDist_Point]);
        Ycoord_Airfoil.push_back(Ycoord[MinDist_Point]);
        Zcoord_Airfoil.push_back(Zcoord[MinDist_Point]);
        
        Xcoord.erase (Xcoord.begin() + MinDist_Point);
        Ycoord.erase (Ycoord.begin() + MinDist_Point);
        Zcoord.erase (Zcoord.begin() + MinDist_Point);
        
        point1_Airfoil.push_back(point1.at(MinDist_Point));
        point2_Airfoil.push_back(point2.at(MinDist_Point));
        
        point1.erase(point1.begin() + MinDist_Point);
        point2.erase(point2.begin() + MinDist_Point);
        
        
        /*--- Algorithm for the rest of the points ---*/
        do {
            
            /*--- Last added point in the list ---*/
            Airfoil_Point = Xcoord_Airfoil.size() - 1;
            
            /*--- Compute the slope of the curve ---*/
            Airfoil_Tangent[0] = Xcoord_Airfoil[Airfoil_Point] - Xcoord_Airfoil[Airfoil_Point-1];
            Airfoil_Tangent[1] = Ycoord_Airfoil[Airfoil_Point] - Ycoord_Airfoil[Airfoil_Point-1];
            Airfoil_Tangent[2] = Zcoord_Airfoil[Airfoil_Point] - Zcoord_Airfoil[Airfoil_Point-1];
            Length = sqrt(pow(Airfoil_Tangent[0], 2.0) + pow(Airfoil_Tangent[1], 2.0) + pow(Airfoil_Tangent[2], 2.0));
            Airfoil_Tangent[0] /= Length; Airfoil_Tangent[1] /= Length; Airfoil_Tangent[2] /= Length;
            
            /*--- Find the closest point with the right slope ---*/
            MinDist_Value = 1E6; MinDistAngle_Value = 180;
            MinDist_Point = -1; MinDistAngle_Point = -1;
            for (iVertex = 0; iVertex < Xcoord.size(); iVertex++) {
                
                Segment[0] = Xcoord[iVertex] - Xcoord_Airfoil[Airfoil_Point];
                Segment[1] = Ycoord[iVertex] - Ycoord_Airfoil[Airfoil_Point];
                Segment[2] = Zcoord[iVertex] - Zcoord_Airfoil[Airfoil_Point];
                
                /*--- Compute the distance to each point ---*/
                Dist_Value = sqrt(pow(Segment[0], 2.0) + pow(Segment[1], 2.0) + pow(Segment[2], 2.0));
                
                /*--- Compute the angle of the point ---*/
                Segment[0] /= Dist_Value; Segment[1] /= Dist_Value; Segment[2] /= Dist_Value;
                
                /*--- Clip the value of the cosine, this is important due to the round errors ---*/
                CosValue = Airfoil_Tangent[0]*Segment[0] + Airfoil_Tangent[1]*Segment[1] + Airfoil_Tangent[2]*Segment[2];
                if (CosValue >= 1.0) CosValue = 1.0;
                if (CosValue <= -1.0) CosValue = -1.0;
                
                Angle_Value = acos(CosValue) * 180 / PI_NUMBER;
                
                if (Dist_Value < MinDist_Value) { MinDist_Point = iVertex; MinDist_Value = Dist_Value; }
                if ((Dist_Value < MinDistAngle_Value) && (Angle_Value < MaxAngle)) {MinDistAngle_Point = iVertex; MinDistAngle_Value = Dist_Value;}
                
            }
            
            if ( MinDistAngle_Point != -1) MinDist_Point = MinDistAngle_Point;
            
            /*--- Add and remove the min distance to the list ---*/
            Xcoord_Airfoil.push_back(Xcoord[MinDist_Point]);
            Ycoord_Airfoil.push_back(Ycoord[MinDist_Point]);
            Zcoord_Airfoil.push_back(Zcoord[MinDist_Point]);
            
            Xcoord.erase(Xcoord.begin() + MinDist_Point);
            Ycoord.erase(Ycoord.begin() + MinDist_Point);
            Zcoord.erase(Zcoord.begin() + MinDist_Point);
            
            point1_Airfoil.push_back(point1.at(MinDist_Point));
            point2_Airfoil.push_back(point2.at(MinDist_Point));
            
            point1.erase(point1.begin() + MinDist_Point);
            point2.erase(point2.begin() + MinDist_Point);
            
        } while (Xcoord.size() != 0);
        
        /*--- Clean the vector before using them again for storing the upper or the lower side ---*/
        Xcoord.clear();
        Ycoord.clear();
        Zcoord.clear();
        
        point1.clear();
        point2.clear();
        
        
        ////cout << endl << "About to start the tracing algorithm!" << endl;
        //
        //  	/*--- THIS IS THE START OF THE AIRFOIL-TRACING ALGORITHM
        //		(it requires vectors of the X, Y, and Z coordinates
        //		of the airfoil section, with no points repeated.) ---*/
        //	unsigned long iNode, trailing_loc;
        //	unsigned long nNode = Xcoord_Airfoil.size();
        //	double max_x;
        //	vector<double> X_dummy(nNode), Y_dummy(nNode), Z_dummy(nNode);
        //
        ////cout << endl << "finished first initializations!" << endl;
        //
        //	/*--- find the trailing-edge point
        //		(we're assuming that the airfoil is facing
        //		left and that it isn't sitting upside down.
        //		If either of those happens to be true, this
        //		algorithm will not work.) ---*/
        //	max_x = -1.0;
        //	for (iNode = 0; iNode < nNode; iNode++) {
        //		if (Xcoord_Airfoil.at(iNode) >= max_x) {
        //			max_x = Xcoord_Airfoil.at(iNode);
        //			trailing_loc = iNode;
        //		}
        //	}
        //
        ////cout << endl << "found the trailing edge!" << endl;
        //	/*--- move the trailing-edge point to the
        //		top of the coordinates lists ---*/
        //	//X_dummy.clear();
        //	//Y_dummy.clear();
        //	//Z_dummy.clear();
        //
        ////cout << endl << "cleared the vectors!" << endl;
        //
        //	X_dummy.at(0) = Xcoord_Airfoil.at(trailing_loc);
        //	Y_dummy.at(0) = Ycoord_Airfoil.at(trailing_loc);
        //	Z_dummy.at(0) = Zcoord_Airfoil.at(trailing_loc);
        //
        //	//X_dummy.push_back(Xcoord_Airfoil.at(trailing_loc));
        //	//Y_dummy.push_back(Ycoord_Airfoil.at(trailing_loc));
        //	//Z_dummy.push_back(Zcoord_Airfoil.at(trailing_loc));
        //
        ////cout << endl << "assigned trailing edge to beginning." << endl;
        //	if (trailing_loc != 0) {
        //		for (iNode = 0; iNode <= trailing_loc-1; iNode++) {
        //		X_dummy.at(iNode+1) = Xcoord_Airfoil.at(iNode);
        //		Y_dummy.at(iNode+1) = Ycoord_Airfoil.at(iNode);
        //		Z_dummy.at(iNode+1) = Zcoord_Airfoil.at(iNode);
        //
        ////cout << endl << "nNode = " << nNode << endl;
        ////cout << endl << "trailing_loc = " << trailing_loc << endl;
        //
        //	//		X_dummy.push_back(Xcoord_Airfoil.at(iNode));
        //	//		Y_dummy.push_back(Ycoord_Airfoil.at(iNode));
        //	//		Z_dummy.push_back(Zcoord_Airfoil.at(iNode));
        //		}
        //	}
        //
        ////cout << endl << "wrote points to just before the trailing_loc" << endl;
        //
        //	for (iNode = trailing_loc+1; iNode < nNode; iNode++) {
        //		X_dummy.at(iNode) = Xcoord_Airfoil.at(iNode);
        //		Y_dummy.at(iNode) = Ycoord_Airfoil.at(iNode);
        //		Z_dummy.at(iNode) = Zcoord_Airfoil.at(iNode);
        //
        //	//	X_dummy.push_back(Xcoord_Airfoil.at(iNode));
        //	//	Y_dummy.push_back(Ycoord_Airfoil.at(iNode));
        //	//	Z_dummy.push_back(Zcoord_Airfoil.at(iNode));
        //	}
        //
        ////cout << endl << "wrote points from the trailing loc to the end of the vector" << endl;
        //
        //	Xcoord_Airfoil.swap(X_dummy);
        //	Ycoord_Airfoil.swap(Y_dummy);
        //	Zcoord_Airfoil.swap(Z_dummy);
        //
        //	X_dummy.clear();
        //	Y_dummy.clear();
        //	Z_dummy.clear();
        //
        ////cout << endl << "moved the trailing edge to the beginning!" << endl;
        //
        //	/*--- name the trailing-edge point, for convenience ---*/
        //	double trailing_point [3];
        //	trailing_point[0] = Xcoord_Airfoil.at(0);
        //	trailing_point[1] = Ycoord_Airfoil.at(0);
        //	trailing_point[2] = Zcoord_Airfoil.at(0);
        //
        //	/*--- find the leading-edge point ---*/
        //	double maxDistance, xDist, yDist, zDist, distance;
        //	double leading_point [3];
        //	maxDistance = 0.0;
        //	for (iNode = 1; iNode < nNode; iNode++) {
        //
        //		/*--- compute the distance of the each
        //			point from the trailing edge ---*/
        //		xDist = trailing_point[0] - Xcoord_Airfoil.at(iNode);
        //		yDist = trailing_point[1] - Ycoord_Airfoil.at(iNode);
        //		zDist = trailing_point[2] - Zcoord_Airfoil.at(iNode);
        //		distance = sqrt(pow(xDist,2.0) + pow(yDist,2.0) + pow(zDist,2.0));
        //
        //		/*--- assume the leading-edge point
        //			is the one that is farthest
        //			from the trailing edge ---*/
        //		if (distance > maxDistance) {
        //			maxDistance = distance;
        //			leading_point[0] = Xcoord_Airfoil.at(iNode);
        //			leading_point[1] = Ycoord_Airfoil.at(iNode);
        //			leading_point[2] = Zcoord_Airfoil.at(iNode);
        //		}
        //	}
        //
        ////cout << endl << "found the leading edge!" << endl;
        //
        //	/*--- find the chord length ---*/
        //	double chord;
        //	xDist = trailing_point[0] - leading_point[0];
        //	yDist = trailing_point[1] - leading_point[1];
        //	zDist = trailing_point[2] - leading_point[2];
        //	chord = sqrt(pow(xDist,2.0) + pow(yDist,2.0) + pow(zDist,2.0));
        //
        //	/*--- define a vector pointing from the
        //		trailing edge to the leading edge ---*/
        //	double chord_line [3];
        //	for (iDim = 0; iDim < nDim; iDim++) {
        //		chord_line[iDim] = leading_point[iDim]-trailing_point[iDim];
        //	}
        //
        //	/*--- define a point that lies "above" the
        //		leading edge, in the plane of the
        //		airfoil. (Here, we need to hope
        //		that the Z-axis points "upward.") ---*/
        //	double above_point [3];
        //	above_point[0] = leading_point[0];
        //	above_point[1] = leading_point[1];
        //	above_point[2] = leading_point[2] + 1;
        //
        //	/*--- define a vector point from the trailing
        //		edge to the newly created "point above" ---*/
        //	double above_line [3];
        //	for (iDim = 0; iDim < nDim; iDim++) {
        //		above_line[iDim] = above_point[iDim] - trailing_point[iDim];
        //	}
        //
        //	/*--- find the vector that points "into the page,"
        //		i.e. if the airfoil is facing  left, this
        //		vector points away from you and lies
        //		normal to the section. ---*/
        //	double section_normal [3];
        //	section_normal[0] = chord_line[1]*above_line[2] - chord_line[2]*above_line[1];
        //	section_normal[1] = chord_line[2]*above_line[0] - chord_line[0]*above_line[2];
        //	section_normal[2] = chord_line[0]*above_line[1] - chord_line[1]*above_line[0];
        //
        //	/*--- beginning at the trailing edge,
        //		connect the points, in order. ---*/
        //	double startingPossibles, magnitude, outward_vec_mag, a_dot_b, cos_angle, weight_dist, weight_ang;
        //	double *dist, **vec, *vec_mag, *angles;
        //	double tangent [3], unit_tangent [3], unit_vec [3], inward_vec [3], outward_vec [3], penalty [2];
        //	vector<double> minDist, list, ang_list;
        //	vector<long> min_loc;
        //	unsigned long jNode, next_loc;
        //	unsigned long min_ang_loc [2];
        //	unsigned short  nPossibles, k, jPoss;
        //	unsigned short ang_index [2];
        //
        //	dist = new double [nNode];
        //
        ////cout << endl << "about to start the loop over the points!" << endl;
        //
        //	for (iNode = 0; iNode < nNode-1; iNode++) {
        //
        //		/*--- as we move around the airfoil, and as more
        //			points are connected, we don't need to
        //			test as many "possibles" at each station.
        //			here, we will begin with four and, by the
        //			end, linearly decrease the number to one. ---*/
        //		startingPossibles = 5.9;
        //		nPossibles = int(floor(((1.0-startingPossibles)/double(nNode))*double(iNode) + startingPossibles));
        //
        //		/*--- initialize the minimum-distance
        //			and corresponding index arrays ---*/
        //		minDist.clear();
        //		min_loc.clear();
        //
        //		minDist.assign(nPossibles, pow(chord,2.0));
        //		min_loc.assign(nPossibles, 0);
        //
        ////cout << endl << "initialized minDist and min_loc" << endl;
        //
        //		/*--- intialize all distances from the
        //			current point equal to zero ---*/
        //		for (jNode = 0; jNode < nNode; jNode++) {
        //			dist[jNode] = 0;
        //		}
        //
        //		/*--- of the points left to be connected,
        //			find the nPossible closest ones ---*/
        //		for (jNode = iNode+1; jNode < nNode; jNode++) {
        //
        //			/*--- compute th distance from the current point ---*/
        //			xDist = Xcoord_Airfoil.at(iNode) - Xcoord_Airfoil.at(jNode);
        //			yDist = Ycoord_Airfoil.at(iNode) - Ycoord_Airfoil.at(jNode);
        //			zDist = Zcoord_Airfoil.at(iNode) - Zcoord_Airfoil.at(jNode);
        //			dist[jNode] = sqrt(pow(xDist,2.0) + pow(yDist,2.0) + pow(zDist,2.0));
        //
        //			/*--- let the "list" be the nPossible minimum distances
        //				and also the newly computed distance ---*/
        //			list = minDist;
        //
        //			list.push_back(dist[jNode]);
        //
        //			/*--- sort the list, in ascending order ---*/
        //			sort(list.begin(),list.end());
        //
        //			/*--- let the new minimum-distance vector hold the
        //
        //				first nPossible entries of the list ---*/
        //			minDist.assign(list.begin(),list.begin()+nPossibles);
        //		}
        //
        ////cout << endl << "found the nPossible closest points" << endl;
        //
        //		/*--- find the index of the four closest nodes.
        //			(N.B. There might be a more elegant way
        //			to do this using unordered_sets. That
        //			being said, this should work, since the
        //			grids will be unstructured and the chance
        //			that two separate points have the same
        //			distance is slim.) ---*/
        //		for (jNode = iNode+1; jNode < nNode; jNode++){
        //
        //			/*--- check against the nPossible distances ---*/
        //			for(k = 0; k < nPossibles; k++) {
        //
        //				/*--- if the distances match,
        //					record the index ---*/
        //				if (minDist.at(k) == dist[jNode]) {
        //					min_loc.at(k) = jNode;
        //				}
        //			}
        //		}
        //
        ////cout << endl << "found the indices of the nPossible closest points" << endl;
        //
        //		/*--- move on to comparing angles. ---*/
        //		if (nPossibles == 1) {
        //
        //			/*--- if we only have one "possible" point,
        //				then there is nothing to compare ---*/
        //			 next_loc = min_loc.at(0);
        //
        //		}
        //		else {
        //			/*--- define vectors pointing from the current
        //				point to each of the nPossible points ---*/
        //			vec = new double * [nPossibles];
        //			for (k = 0; k < nPossibles; k ++) {
        //				vec[k] = new double [3];
        //			}
        //			vec_mag = new double [nPossibles];
        //
        //			for (jPoss = 0; jPoss < nPossibles; jPoss++) {
        //				vec[jPoss][0] = Xcoord_Airfoil.at(min_loc.at(jPoss)) - Xcoord_Airfoil.at(iNode);
        //				vec[jPoss][1] = Ycoord_Airfoil.at(min_loc.at(jPoss)) - Ycoord_Airfoil.at(iNode);
        //				vec[jPoss][2] = Zcoord_Airfoil.at(min_loc.at(jPoss)) - Zcoord_Airfoil.at(iNode);
        //				vec_mag[jPoss] = sqrt(pow(vec[jPoss][0],2.0) + pow(vec[jPoss][1],2.0) + pow(vec[jPoss][2],2.0));
        //			}
        //
        //			/*-- find the vector that points "inward." if there is a
        //				previous point, then use it to find a unit tangent
        //				vector along the surface. use this tangent as a
        //				component of the "inward" vector, along with the
        //				nPossible directions. If we're at the trailing edge,
        //				where there is not previous point, then just use
        //				the "mean" of the nPossible directions. ---*/
        //			if (iNode > 1) {
        //
        //				/*-- define the tangent vector ---*/
        //				tangent[0] = Xcoord_Airfoil.at(iNode) - Xcoord_Airfoil.at(iNode-1);
        //				tangent[1] = Ycoord_Airfoil.at(iNode) - Ycoord_Airfoil.at(iNode-1);
        //				tangent[2] = Zcoord_Airfoil.at(iNode) - Zcoord_Airfoil.at(iNode-1);
        //				magnitude = sqrt(pow(tangent[0],2.0) + pow(tangent[1],2.0) + pow(tangent[2],2.0));
        //
        //				/*--- define the unit tangent vector ---*/
        //				for (iDim = 0; iDim < nDim; iDim++) {
        //					unit_tangent[iDim] = tangent[iDim]/magnitude;
        //				}
        //
        //				/*--- initialize the "inward"-pointing vector ---*/
        //				for (iDim = 0; iDim < nDim; iDim++) {
        //					inward_vec[iDim] = unit_tangent[iDim];
        //				}
        //			}
        //			else {
        //				/*--- if there is no tangent vector to include ---*/
        //				for (iDim = 0; iDim < nDim; iDim ++) {
        //					inward_vec[iDim] = 0.0;
        //				}
        //			}
        //
        //			/*--- add the nPossible directions as components
        //				of the "inward"-pointing vector ---*/
        //			for (jPoss = 0; jPoss < nPossibles; jPoss++){
        //
        //				for (iDim = 0; iDim < nDim; iDim++) {
        //
        //					/*--- compute the i-th component of the unit vector
        //						pointing in the j-th "possible" direction ---*/
        //					unit_vec[iDim] = vec[jPoss][iDim]/vec_mag[jPoss];
        //
        //					/*--- add it to the definition of
        //						the "inward" vector ---*/
        //					inward_vec[iDim] += unit_vec[iDim];
        //				}
        //			}
        //
        ////cout << endl << "set the direction of the \"inward\" vector" << endl;
        //
        //			/*--- rotate the "inward" vector clockwise 90 degrees
        //				to get the "outward"-facing vector ---*/
        //			outward_vec[0] = section_normal[1]*inward_vec[2] - section_normal[2]*inward_vec[1];
        //			outward_vec[1] = section_normal[2]*inward_vec[0] - section_normal[0]*inward_vec[2];
        //			outward_vec[2] = section_normal[0]*inward_vec[1] - section_normal[1]*inward_vec[0];
        //			outward_vec_mag = sqrt(pow(outward_vec[0],2.0) + pow(outward_vec[1],2.0) + pow(outward_vec[2],2.0));
        //
        //			/*--- compute the angles between the "outward" vector and the nPossible "possible" points ---*/
        //			angles = new double [nPossibles];
        //			for (jPoss = 0; jPoss < nPossibles; jPoss++) {
        //
        //				/*--- dot the "outward" vector with the j-th possible direction ---*/
        //				a_dot_b = outward_vec[0]*vec[jPoss][0] + outward_vec[1]*vec[jPoss][1] + outward_vec[2]*vec[jPoss][2];
        //				/*--- compute the cosine of the angle between them ---*/
        //				cos_angle = a_dot_b/(outward_vec_mag*vec_mag[jPoss]);
        //
        //				/*--- store the actual angle ---*/
        //				angles[jPoss] = acos(cos_angle);
        //			}
        //
        //			/*--- copy the nPossible angles to the angle-list ---*/
        //			ang_list.assign(angles,angles+nPossibles);
        //
        //			/*--- sort the angle-list ---*/
        //			sort(ang_list.begin(),ang_list.end());
        //
        //			/*--- find the indices corresponding to the
        //				two smallest angles in the angle-list ---*/
        //			for (jPoss = 0; jPoss < 2; jPoss++) {
        //
        //				/*--- check against the
        //					nPossible angles ---*/
        //				for (k = 0; k < nPossibles; k++) {
        //
        //					/*--- if there is a match,
        //						note the index ---*/
        //					if (ang_list.at(jPoss) == angles[k]) {
        //						min_ang_loc[jPoss] = min_loc.at(k);
        //						ang_index[jPoss] = k;
        //					}
        //				}
        //			}
        //
        //			/*--- assign nondimensional penalties for long distanes and large angles
        //				(distances are nondimensionalized by chord length, while angles
        //				are nondimensionalized by 90 degrees) ---*/
        //			weight_dist = 2.0;
        //			weight_ang = 1.0;
        //			for (jPoss = 0; jPoss < 2; jPoss++) {
        //				penalty[jPoss] = weight_dist*(dist[min_ang_loc[jPoss]]/chord) + weight_ang*(angles[ang_index[jPoss]]/(3.14159/2));
        //			}
        //
        //			/*--- of the two points with the smallest angle,
        //				choose the one with the smaller penalty ---*/
        //			if (penalty[1] > penalty[0]) {
        //				next_loc = min_ang_loc[0];
        //			}
        //			else {
        //				next_loc = min_ang_loc[1];
        //			}
        //
        //			/*--- delete dynamically allocated memory
        //				within the else statement ---*/
        //			for (k = 0; k < nPossibles; k++) {
        //				delete [] vec[k];
        //			}
        //			delete [] vec;
        //			delete [] vec_mag;
        //			delete [] angles;
        //		}
        //
        //		/*--- re-initialize the dummy coordinate
        //			vectors to a length of nNode ---*/
        //		X_dummy.assign(nNode,0.0);
        //		Y_dummy.assign(nNode,0.0);
        //		Z_dummy.assign(nNode,0.0);
        //
        //		/*--- move the newly found "next point" to the
        //			appropriate location in the coordiates vectors ---*/
        //		for (jNode = 0; jNode <= iNode; jNode++) {
        //			X_dummy.at(jNode) = Xcoord_Airfoil.at(jNode);
        //			Y_dummy.at(jNode) = Ycoord_Airfoil.at(jNode);
        //			Z_dummy.at(jNode) = Zcoord_Airfoil.at(jNode);
        //		}
        //
        //		X_dummy.at(iNode+1) = Xcoord_Airfoil.at(next_loc);
        //		Y_dummy.at(iNode+1) = Ycoord_Airfoil.at(next_loc);
        //		Z_dummy.at(iNode+1) = Zcoord_Airfoil.at(next_loc);
        //
        //		for (jNode = iNode+1; jNode <= next_loc-1; jNode++) {
        //			X_dummy.at(jNode+1) = Xcoord_Airfoil.at(jNode);
        //			Y_dummy.at(jNode+1) = Ycoord_Airfoil.at(jNode);
        //			Z_dummy.at(jNode+1) = Zcoord_Airfoil.at(jNode);
        //		}
        //
        //		for (jNode = next_loc+1; jNode < nNode; jNode++) {
        //			X_dummy.at(jNode) = Xcoord_Airfoil.at(jNode);
        //			Y_dummy.at(jNode) = Ycoord_Airfoil.at(jNode);
        //			Z_dummy.at(jNode) = Zcoord_Airfoil.at(jNode);
        //		}
        //
        //		Xcoord_Airfoil.swap(X_dummy);
        //		Ycoord_Airfoil.swap(Y_dummy);
        //		Zcoord_Airfoil.swap(Z_dummy);
        //
        //		X_dummy.clear();
        //		Y_dummy.clear();
        //		Z_dummy.clear();
        //
        //	}
        //
        //	delete [] dist;
        
        
    	/*--- Write the output file (tecplot format) ---*/
    	if (original_surface == true) {
      		
            ofstream Tecplot_File;
      		if (iSection == 0) Tecplot_File.open("Airfoil_Sections.plt", ios::out);
      		else Tecplot_File.open("Airfoil_Sections.plt", ios::app);
            
      		if (iSection == 0) {
        		Tecplot_File << "TITLE = \"Wing airfoil sections\"" << endl;
        		Tecplot_File << "VARIABLES = \"X\",\"Y\",\"Z\"" << endl;
      		}
            
      		Tecplot_File << "ZONE T=\"SECTION_"<< (iSection+1) << "\", NODES= "<< Xcoord_Airfoil.size() << ", ELEMENTS= " << Xcoord_Airfoil.size()-1 << ", DATAPACKING= POINT, ZONETYPE= FELINESEG" << endl;
            
      		/*--- Coordinates ---*/
      		for (iVertex = 0; iVertex < Xcoord_Airfoil.size(); iVertex++) {
        		Tecplot_File << Xcoord_Airfoil[iVertex] <<" "<< Ycoord_Airfoil[iVertex] <<" "<< Zcoord_Airfoil[iVertex] << endl;
      		}
      		
            /*--- Conectivity ---*/
      		for (iVertex = 1; iVertex < Xcoord_Airfoil.size(); iVertex++) {
        		Tecplot_File << iVertex << "\t" << iVertex+1 << "\n";
      		}
            
      		Tecplot_File.close();
    	}
        
    }
    
}

CMultiGridGeometry::CMultiGridGeometry(CGeometry ***geometry, CConfig **config_container, unsigned short iMesh, unsigned short iZone) : CGeometry() {
    
	/*--- CGeometry & CConfig pointers to the fine grid level for clarity. We may
     need access to the other zones in the mesh for zone/sliding boundaries. ---*/
	CGeometry *fine_grid = geometry[iZone][iMesh-1];
	CConfig *config = config_container[iZone];
    
	/*--- Local variables ---*/
	unsigned long iPoint, Index_CoarseCV, CVPoint, iElem, iVertex, jPoint;
	bool agglomerate_seed = true, agglomerate_CV = true;
	unsigned short nChildren, iNode, counter, iMarker, jMarker, Marker_Boundary;
	short marker_seed;
    
	unsigned short priority;
    
	/*--- Set the boolean to indicate that this is a coarse multigrid level. ---*/
	FinestMGLevel = false;
    
	unsigned short max_children = config->GetMaxChildren();
    
    
	/*--- Create a queue system to deo the agglomeration ---*/
	/*--- 1st) More than two markers ---> Vertices (never agglomerate)                             ---*/
	/*--- 2nd) Two markers ---> Edges (agglomerate if same BC, never agglomerate if different BC)  ---*/
	/*--- 3rd) One marker ---> Surface (always agglomarate)                                        ---*/
	/*--- 4th) No marker ---> Internal Volume (always agglomarate)                                 ---*/
    
	/*--- Set a marker to indicate indirect agglomeration ---*/
	if (iMesh == MESH_1) {
		for (iPoint = 0; iPoint < fine_grid->GetnPoint(); iPoint ++)
			fine_grid->node[iPoint]->SetAgglomerate_Indirect(false);
		for (iElem = 0; iElem < fine_grid->GetnElem(); iElem++)
			if ((fine_grid->elem[iElem]->GetVTK_Type() == HEXAHEDRON) ||
                (fine_grid->elem[iElem]->GetVTK_Type() == RECTANGLE))
				for (iNode = 0; iNode < fine_grid->elem[iElem]->GetnNodes(); iNode++) {
					iPoint = fine_grid->elem[iElem]->GetNode(iNode);
					fine_grid->node[iPoint]->SetAgglomerate_Indirect(true);
				}
	}
    
	/*--- Write the number of dimensions of the coarse grid ---*/
	nDim = fine_grid->GetnDim();
    
	/*--- Create the coarse grid structure using as baseline the fine grid ---*/
	CMultiGridQueue MGQueue_InnerCV(fine_grid->GetnPoint());
    
	node = new CPoint*[fine_grid->GetnPoint()];
	for (iPoint = 0; iPoint < fine_grid->GetnPoint(); iPoint ++) {
		/*--- Create node structure ---*/
		node[iPoint] = new CPoint(nDim, iPoint, config);
		/*--- Set the indirect agglomeration to false ---*/
		node[iPoint]->SetAgglomerate_Indirect(false);
	}
    
	Index_CoarseCV = 0;
    
	/*--- The first step is the boundary agglomeration. ---*/
	for (iMarker = 0; iMarker < fine_grid->GetnMarker(); iMarker++) {
		Marker_Boundary = config->GetMarker_All_Boundary(iMarker);
        
		for (iVertex = 0; iVertex < fine_grid->GetnVertex(iMarker); iVertex++) {
			iPoint = fine_grid->vertex[iMarker][iVertex]->GetNode();
            
			/*--- If the element has not being previously agglomerated and it belongs to the physical domain,
			 then the agglomeration is studied ---*/
			if ((fine_grid->node[iPoint]->GetAgglomerate() == false) &&
                (fine_grid->node[iPoint]->GetDomain()) &&
                (GeometricalCheck(iPoint, fine_grid, config))) {
                
				nChildren = 1;
                
				/*--- We set an index for the parent control volume ---*/
				fine_grid->node[iPoint]->SetParent_CV(Index_CoarseCV);
                
				/*--- We add the seed point (child) to the parent control volume ---*/
				node[Index_CoarseCV]->SetChildren_CV(0,iPoint);
				agglomerate_seed = true; counter = 0; marker_seed = iMarker;
                
				/*--- For a particular point in the fine grid we save all the markers that are in that point ---*/
				unsigned short copy_marker[MAX_NUMBER_MARKER];
				for (jMarker = 0; jMarker < fine_grid->GetnMarker(); jMarker ++)
					if (fine_grid->node[iPoint]->GetVertex(jMarker) != -1) {
						copy_marker[counter] = jMarker;
						counter++;
					}
                
				/*--- To aglomerate a vertex it must have only one physical bc!!
				 This can be improved. ---*/
                
				/*--- If there is only a marker, it is a good candidate for agglomeration ---*/
				if (counter == 1) agglomerate_seed = true;
                
				/*--- If there are two markers, we will aglomerate if one of the marker is SEND_RECEIVE ---*/
				if (counter == 2) {
					if ((config->GetMarker_All_Boundary(copy_marker[0]) == SEND_RECEIVE) ||
                        (config->GetMarker_All_Boundary(copy_marker[1]) == SEND_RECEIVE)) agglomerate_seed = true;
					else agglomerate_seed = false;
				}
                
				/*--- If there are more than 2 markers on the , the aglomeration will be discarted ---*/
				if (counter > 2) agglomerate_seed = false;
                
				/*--- If the seed can be agglomerated, we try to agglomerate more points ---*/
				if (agglomerate_seed) {
                    
					/*--- Now we do a sweep over all the nodes that surround the seed point ---*/
					for (iNode = 0; iNode <	fine_grid->node[iPoint]->GetnPoint(); iNode ++) {
                        
						agglomerate_CV = false;
						CVPoint = fine_grid->node[iPoint]->GetPoint(iNode);
                        
						/*--- Determine if the CVPoint can be agglomerated ---*/
						agglomerate_CV = SetBoundAgglomeration(CVPoint, marker_seed, fine_grid, config);
                        
						/*--- The new point can be agglomerated ---*/
						if (agglomerate_CV && (nChildren < max_children))  {
                            
							/*--- We set the value of the parent ---*/
							fine_grid->node[CVPoint]->SetParent_CV(Index_CoarseCV);
                            
							/*--- We set the value of the child ---*/
							node[Index_CoarseCV]->SetChildren_CV(nChildren,CVPoint);
							nChildren++;
						}
					}
                    
					vector<unsigned long> Suitable_Indirect_Neighbors;
					if (fine_grid->node[iPoint]->GetAgglomerate_Indirect())
						SetSuitableNeighbors(&Suitable_Indirect_Neighbors, iPoint, Index_CoarseCV, fine_grid);
                    
					/*--- Now we do a sweep over all the indirect nodes that can be added ---*/
					for (iNode = 0; iNode <	Suitable_Indirect_Neighbors.size(); iNode ++) {
						agglomerate_CV = false;
						CVPoint = Suitable_Indirect_Neighbors[iNode];
                        
						/*--- Determine if the CVPoint can be agglomerated ---*/
						agglomerate_CV = SetBoundAgglomeration(CVPoint, marker_seed, fine_grid, config);
                        
						/*--- The new point can be agglomerated ---*/
						if (agglomerate_CV && (nChildren < max_children))  {
                            
							/*--- We set the value of the parent ---*/
							fine_grid->node[CVPoint]->SetParent_CV(Index_CoarseCV);
                            
							/*--- We set the indirect agglomeration information ---*/
							if (fine_grid->node[CVPoint]->GetAgglomerate_Indirect())
								node[Index_CoarseCV]->SetAgglomerate_Indirect(true);
                            
							/*--- We set the value of the child ---*/
							node[Index_CoarseCV]->SetChildren_CV(nChildren,CVPoint);
							nChildren++;
						}
					}
				}
                
				/*--- Update the number of child of the control volume ---*/
				node[Index_CoarseCV]->SetnChildren_CV(nChildren);
				Index_CoarseCV++;
			}
		}
	}
    
	/*--- Agglomerate all the nodes that have more than one physical boundary condition,
	 Maybe here we can add the posibility of merging the vertex that have the same number,
	 and kind  of markers---*/
	for (iMarker = 0; iMarker < fine_grid->GetnMarker(); iMarker++)
		for (iVertex = 0; iVertex < fine_grid->GetnVertex(iMarker); iVertex++) {
			iPoint = fine_grid->vertex[iMarker][iVertex]->GetNode();
			if ((fine_grid->node[iPoint]->GetAgglomerate() == false) &&
                (fine_grid->node[iPoint]->GetDomain())) {
				fine_grid->node[iPoint]->SetParent_CV(Index_CoarseCV);
				node[Index_CoarseCV]->SetChildren_CV(0,iPoint);
				node[Index_CoarseCV]->SetnChildren_CV(1);
				Index_CoarseCV++;
			}
		}
    
	/*--- Update the queue with the results from the boundary agglomeration ---*/
	for (iPoint = 0; iPoint < fine_grid->GetnPoint(); iPoint ++) {
		/*--- The CV has been agglomerated, remove form the list ---*/
		if (fine_grid->node[iPoint]->GetAgglomerate() == true) {
			MGQueue_InnerCV.RemoveCV(iPoint);
		}
		else {
			/*--- Count the number of agglomerated neighbors, and modify the queue ---*/
			priority = 0;
			for (iNode = 0; iNode <	fine_grid->node[iPoint]->GetnPoint(); iNode ++) {
				jPoint = fine_grid->node[iPoint]->GetPoint(iNode);
				if (fine_grid->node[jPoint]->GetAgglomerate() == true) priority++;
			}
			MGQueue_InnerCV.MoveCV(iPoint, priority);
		}
	}
    
	/*--- Agglomerate the domain nodes ---*/
	//		for (iPoint = 0; iPoint < fine_grid->GetnPoint(); iPoint ++) {
	unsigned long iteration = 0;
	while (!MGQueue_InnerCV.EmptyQueue() && (iteration < fine_grid->GetnPoint())) {
		iPoint = MGQueue_InnerCV.NextCV();
		iteration ++;
        
		/*--- If the element has not being previously agglomerated, belongs to the physical domain,
		 and satisfies several geometrical criteria then the seed CV is acepted for agglomeration ---*/
		if ((fine_grid->node[iPoint]->GetAgglomerate() == false) &&
            (fine_grid->node[iPoint]->GetDomain()) &&
            (GeometricalCheck(iPoint, fine_grid, config))) {
            
			nChildren = 1;
            
			/*--- We set an index for the parent control volume ---*/
			fine_grid->node[iPoint]->SetParent_CV(Index_CoarseCV);
            
			/*--- We add the seed point (child) to the parent control volume ---*/
			node[Index_CoarseCV]->SetChildren_CV(0, iPoint);
            
			/*--- Update the queue with the seed point (remove the seed and
			 increase the priority of the neighbors) ---*/
			MGQueue_InnerCV.Update(iPoint, fine_grid);
            
			/*--- Now we do a sweep over all the nodes that surround the seed point ---*/
			for (iNode = 0; iNode <	fine_grid->node[iPoint]->GetnPoint(); iNode ++) {
                
				CVPoint = fine_grid->node[iPoint]->GetPoint(iNode);
                
				/*--- Determine if the CVPoint can be agglomerated ---*/
				if ((fine_grid->node[CVPoint]->GetAgglomerate() == false) &&
                    (fine_grid->node[CVPoint]->GetDomain()) &&
                    (GeometricalCheck(CVPoint, fine_grid, config))) {
                    
					/*--- The new point can be agglomerated, note that the applicability of max_children depend on the seed ---*/
					if (nChildren < max_children)  {
                        
						/*--- We set the value of the parent ---*/
						fine_grid->node[CVPoint]->SetParent_CV(Index_CoarseCV);
                        
						/*--- We set the value of the child ---*/
						node[Index_CoarseCV]->SetChildren_CV(nChildren, CVPoint);
						nChildren++;
                        
						/*--- Update the queue with the new control volume (remove the CV and
						 increase the priority of the neighbors) ---*/
						MGQueue_InnerCV.Update(CVPoint, fine_grid);
                        
					}
				}
                
			}
            
			/*--- Subrotuine to identify the indirect neighbors ---*/
			vector<unsigned long> Suitable_Indirect_Neighbors;
			if (fine_grid->node[iPoint]->GetAgglomerate_Indirect())
				SetSuitableNeighbors(&Suitable_Indirect_Neighbors, iPoint, Index_CoarseCV, fine_grid);
            
			/*--- Now we do a sweep over all the indirect nodes that can be added ---*/
			for (iNode = 0; iNode <	Suitable_Indirect_Neighbors.size(); iNode ++) {
                
				agglomerate_CV = false;
				CVPoint = Suitable_Indirect_Neighbors[iNode];
                
				/*--- Determine if the CVPoint can be agglomerated ---*/
				if ((fine_grid->node[CVPoint]->GetAgglomerate() == false) && (fine_grid->node[CVPoint]->GetDomain()))
					agglomerate_CV = true;
                
				/*--- The new point can be agglomerated ---*/
				if ((agglomerate_CV) && (nChildren < max_children))  {
                    
					/*--- We set the value of the parent ---*/
					fine_grid->node[CVPoint]->SetParent_CV(Index_CoarseCV);
                    
					/*--- We set the indirect agglomeration information ---*/
					if (fine_grid->node[CVPoint]->GetAgglomerate_Indirect())
						node[Index_CoarseCV]->SetAgglomerate_Indirect(true);
                    
					/*--- We set the value of the child ---*/
					node[Index_CoarseCV]->SetChildren_CV(nChildren, CVPoint);
					nChildren++;
                    
					/*--- Update the queue with the new control volume (remove the CV and
					 increase the priority of the neighbors) ---*/
					MGQueue_InnerCV.Update(CVPoint, fine_grid);
                    
				}
			}
            
			/*--- Update the number of control of childrens ---*/
			node[Index_CoarseCV]->SetnChildren_CV(nChildren);
			Index_CoarseCV++;
		}
		else {
			/*--- The seed point can not be agglomerated because of size, domain, streching, etc.
			 move the point to the lowest priority ---*/
			MGQueue_InnerCV.MoveCV(iPoint, -1);
		}
        
	}
    
	/*--- Add all the elements that have not being agglomerated, in the previous stage ---*/
	for (iPoint = 0; iPoint < fine_grid->GetnPoint(); iPoint ++) {
		if ((fine_grid->node[iPoint]->GetAgglomerate() == false) && (fine_grid->node[iPoint]->GetDomain())) {
			nChildren = 1;
			fine_grid->node[iPoint]->SetParent_CV(Index_CoarseCV);
			if (fine_grid->node[iPoint]->GetAgglomerate_Indirect())
				node[Index_CoarseCV]->SetAgglomerate_Indirect(true);
			node[Index_CoarseCV]->SetChildren_CV(0,iPoint);
			node[Index_CoarseCV]->SetnChildren_CV(nChildren);
			Index_CoarseCV++;
		}
	}
    
    
	/*--- Dealing with MPI parallelization, the objective is that the received nodes must be agglomerated
	 in the same way as the donor nodes. ---*/
    
	nPointDomain = Index_CoarseCV;
    
	unsigned long jVertex;
	short Send_Recv;
    
	/*--- Send the node agglomeration information of the donor
	 (parent and children), Sending only occurs with MPI ---*/
    
	/*--- Receive the donor agglomeration (parent and children) ---*/
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
		Marker_Boundary = config->GetMarker_All_Boundary(iMarker);
		Send_Recv = config->GetMarker_All_SendRecv(iMarker);
        
		if ((Marker_Boundary == SEND_RECEIVE) && (Send_Recv < 0)) {
            
			int receive_from = abs(Send_Recv)-1;
			unsigned long nBuffer = fine_grid->nVertex[iMarker];
            
			unsigned long *Buffer_Receive_Parent = new unsigned long [nBuffer];
			unsigned long *Buffer_Receive_Children = new unsigned long [nBuffer];
            
			/*--- Retrieve the donor information from the matching marker ---*/
			unsigned short donorZone   = 0;
			unsigned short donorMarker = 1;
			unsigned long  donorPoint;
            
			/*--- Get the information from the donor directly. This is a serial
			 computation with access to all nodes. Note that there is an
			 implicit ordering in the list. ---*/
			for (iVertex = 0; iVertex < fine_grid->nVertex[iMarker]; iVertex++) {
                
				/*--- Check the donor zone in case there is more than one in the mesh. ---*/
				donorZone = fine_grid->vertex[iMarker][iVertex]->GetMatching_Zone();
                
				/*--- For now, search for donor marker for every receive point. Probably
                 a more efficient way to do this in the future. ---*/
				for (unsigned short iMark = 0; iMark < config_container[donorZone]->GetnMarker_All(); iMark++){
					if (config_container[donorZone]->GetMarker_All_SendRecv(iMark)-1 == receive_from) donorMarker = iMark;
				}
                
				donorPoint = geometry[donorZone][iMesh-1]->vertex[donorMarker][iVertex]->GetNode();
				Buffer_Receive_Children[iVertex] = donorPoint;
				Buffer_Receive_Parent[iVertex] = geometry[donorZone][iMesh-1]->node[donorPoint]->GetParent_CV();
			}
            
			/*--- Create a list of the parent nodes without repeated parents ---*/
			vector<unsigned long>::iterator it;
			vector<unsigned long> Aux_Parent;
            
			for (iVertex = 0; iVertex < fine_grid->nVertex[iMarker]; iVertex++)
				Aux_Parent.push_back (Buffer_Receive_Parent[iVertex]);
            
			sort(Aux_Parent.begin(), Aux_Parent.end());
			it = unique(Aux_Parent.begin(), Aux_Parent.end());
			Aux_Parent.resize(it - Aux_Parent.begin());
            
			unsigned long *Parent_Remote = new unsigned long[fine_grid->nVertex[iMarker]];
			unsigned long *Children_Remote = new unsigned long[fine_grid->nVertex[iMarker]];
			unsigned long *Parent_Local = new unsigned long[fine_grid->nVertex[iMarker]];
			unsigned long *Children_Local = new unsigned long[fine_grid->nVertex[iMarker]];
            
			/*--- Create the local vector and remote for the parents and the children ---*/
			for (iVertex = 0; iVertex < fine_grid->nVertex[iMarker]; iVertex++) {
				Parent_Remote[iVertex] = Buffer_Receive_Parent[iVertex];
                
				/*--- We use the same sorting as in the donor domain ---*/
				for (jVertex = 0; jVertex < Aux_Parent.size(); jVertex++) {
					if (Parent_Remote[iVertex] == Aux_Parent[jVertex]) {
						Parent_Local[iVertex] = jVertex + Index_CoarseCV;
						break;
					}
				}
                
				Children_Remote[iVertex] = Buffer_Receive_Children[iVertex];
				Children_Local[iVertex] = fine_grid->vertex[iMarker][iVertex]->GetNode();
                
			}
            
			Index_CoarseCV += Aux_Parent.size();
            
			unsigned short *nChildren_ = new unsigned short [Index_CoarseCV];
			for (unsigned long iParent = 0; iParent < Index_CoarseCV; iParent++)
				nChildren_[iParent] = 0;
            
			/*--- Create the final structure ---*/
			for (iVertex = 0; iVertex < fine_grid->nVertex[iMarker]; iVertex++) {
				/*--- Be careful, it is possible that a node change the agglomeration configuration, the priority
				 is always, when receive the information ---*/
				fine_grid->node[Children_Local[iVertex]]->SetParent_CV(Parent_Local[iVertex]);
				node[Parent_Local[iVertex]]->SetChildren_CV(nChildren_[Parent_Local[iVertex]],Children_Local[iVertex]);
				nChildren_[Parent_Local[iVertex]]++;
				node[Parent_Local[iVertex]]->SetnChildren_CV(nChildren_[Parent_Local[iVertex]]);
				node[Parent_Local[iVertex]]->SetDomain(false);
			}
            
			delete[] nChildren_;
			delete[] Buffer_Receive_Parent;
			delete[] Buffer_Receive_Children;
			delete[] Parent_Remote;
			delete[] Children_Remote;
			delete[] Parent_Local;
			delete[] Children_Local;
		}
	}
    
	nPoint = Index_CoarseCV;
    
	cout << "CVs of the MG level: " << nPoint << ". Agglom. rate 1/" << double(fine_grid->GetnPoint())/double(nPoint) <<". MG level: "<< iMesh <<"."<< endl;
}

CMultiGridGeometry::~CMultiGridGeometry(void) {
    
}

bool CMultiGridGeometry::SetBoundAgglomeration(unsigned long CVPoint, short marker_seed, CGeometry *fine_grid, CConfig *config) {
    
	bool agglomerate_CV = false;
	unsigned short counter, jMarker;
    
	/*--- Basic condition, the element has not being previously agglomerated and, it belong to the domain ---*/
	if ((fine_grid->node[CVPoint]->GetAgglomerate() == false) &&
        (fine_grid->node[CVPoint]->GetDomain()) &&
        (GeometricalCheck(CVPoint, fine_grid, config))) {
        
		/*--- If the element belong to the boundary, we must be careful ---*/
		if (fine_grid->node[CVPoint]->GetBoundary()) {
            
			/*--- Identify the markers of the vertex that we whant to agglomerate ---*/
			counter = 0;
			unsigned short copy_marker[MAX_NUMBER_MARKER];
			for (jMarker = 0; jMarker < fine_grid->GetnMarker(); jMarker ++)
				if (fine_grid->node[CVPoint]->GetVertex(jMarker) != -1) {
					copy_marker[counter] = jMarker;
					counter++;
				}
            
			/*--- The basic condition is that the aglomerated vertex must have the same physical marker,
			 but eventually a send-receive condition ---*/
            
			/*--- Only one marker in the vertex that is going to be aglomerated ---*/
			if (counter == 1) {
                
				/*--- We agglomerate if there is only a marker and is the same marker as the seed marker ---*/
				if (copy_marker[0] == marker_seed)
					agglomerate_CV = true;
                
				/*--- If there is only a marker, but the marker is the SEND_RECEIVE ---*/
				if (config->GetMarker_All_Boundary(copy_marker[0]) == SEND_RECEIVE)
					agglomerate_CV = true;
			}
            
			/*--- If there are two markers in the vertex that is going to be aglomerated ---*/
			if (counter == 2) {
                
				/*--- First we verify that the seed is a physical boundary ---*/
				if (config->GetMarker_All_Boundary(marker_seed) != SEND_RECEIVE) {
                    
					/*--- Then we check that one of the marker is equal to the seed marker, and the other is send/receive ---*/
					if (((copy_marker[0] == marker_seed) && (config->GetMarker_All_Boundary(copy_marker[1]) == SEND_RECEIVE)) ||
                        ((config->GetMarker_All_Boundary(copy_marker[0]) == SEND_RECEIVE) && (copy_marker[1] == marker_seed)))
						agglomerate_CV = true;
				}
			}
            
		}
        
		/*--- If the element belong to the domain, it is allways aglomerated ---*/
		else {
			agglomerate_CV = true;
		}
	}
    
	return agglomerate_CV;
}


bool CMultiGridGeometry::GeometricalCheck(unsigned long iPoint, CGeometry *fine_grid, CConfig *config) {
    
	double max_dimension = 1.0/config->GetMaxDimension();
    
	/*--- Evaluate the total size of the element ---*/
	bool Volume = true;
	double ratio = pow(fine_grid->node[iPoint]->GetVolume(), 1.0/double(nDim))*max_dimension;
	double limit = pow(config->GetDomainVolume(), 1.0/double(nDim));
	if ( ratio > limit ) Volume = false;
    
	/*--- Evaluate the stretching of the element ---*/
	bool Stretching = true;
	/*	unsigned short iNode, iDim;
     unsigned long jPoint;
     double *Coord_i = fine_grid->node[iPoint]->GetCoord();
     double max_dist = 0.0 ; double min_dist = 1E20;
     for (iNode = 0; iNode <	fine_grid->node[iPoint]->GetnPoint(); iNode ++) {
     jPoint = fine_grid->node[iPoint]->GetPoint(iNode);
     double *Coord_j = fine_grid->node[jPoint]->GetCoord();
     double distance = 0.0;
     for (iDim = 0; iDim < nDim; iDim++)
     distance += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
     distance = sqrt(distance);
     max_dist = max(distance, max_dist);
     min_dist = min(distance, min_dist);
     }
     if ( max_dist/min_dist > 100.0 ) Stretching = false;*/
    
	return (Stretching && Volume);
    
}

void CMultiGridGeometry::SetSuitableNeighbors(vector<unsigned long> *Suitable_Indirect_Neighbors, unsigned long iPoint,
                                              unsigned long Index_CoarseCV, CGeometry *fine_grid) {
    
	unsigned long jPoint, kPoint,  iOriginNeighbor, lPoint;
	unsigned short iNode, jNode, iNeighbor, jNeighbor, kNode;
	bool SecondNeighborSeed, check_1, ThirdNeighborSeed;
    
	/*--- Create a list with the first neighbors, including the seed ---*/
	vector<unsigned long> First_Neighbor_Points;
	First_Neighbor_Points.push_back(iPoint);
	for (iNode = 0; iNode <	fine_grid->node[iPoint]->GetnPoint(); iNode ++) {
		jPoint = fine_grid->node[iPoint]->GetPoint(iNode);
		First_Neighbor_Points.push_back(jPoint);
	}
    
	/*--- Create a list with the second neighbors, without first, and seed neighbors ---*/
	vector<unsigned long> Second_Neighbor_Points, Second_Origin_Points, Suitable_Second_Neighbors;
    
	for (iNode = 0; iNode <	fine_grid->node[iPoint]->GetnPoint(); iNode ++) {
		jPoint = fine_grid->node[iPoint]->GetPoint(iNode);
		for (jNode = 0; jNode <	fine_grid->node[jPoint]->GetnPoint(); jNode ++) {
			kPoint = fine_grid->node[jPoint]->GetPoint(jNode);
            
			/*--- Check that the second neighbor do not belong to the first neighbor or the seed ---*/
			SecondNeighborSeed = true;
			for (iNeighbor = 0; iNeighbor <	First_Neighbor_Points.size(); iNeighbor ++)
				if (kPoint == First_Neighbor_Points[iNeighbor]) {
					SecondNeighborSeed = false;
					break;
				}
            
			if (SecondNeighborSeed) {
				Second_Neighbor_Points.push_back(kPoint);
				Second_Origin_Points.push_back(jPoint);
			}
            
		}
	}
    
	/*---  Identify those second neighbors that are repeated (candidate to be added) ---*/
	for (iNeighbor = 0; iNeighbor <	Second_Neighbor_Points.size(); iNeighbor ++)
		for (jNeighbor = 0; jNeighbor <	Second_Neighbor_Points.size(); jNeighbor ++)
            
        /*--- Repeated second neighbor with different origin ---*/
			if ((Second_Neighbor_Points[iNeighbor] == Second_Neighbor_Points[jNeighbor]) &&
                (Second_Origin_Points[iNeighbor] != Second_Origin_Points[jNeighbor]) &&
                (iNeighbor < jNeighbor)) {
                
				/*--- Check that the origin nodes are not neighbor ---*/
				check_1 = true;
				for (iNode = 0; iNode <	fine_grid->node[Second_Origin_Points[iNeighbor]]->GetnPoint(); iNode ++) {
					iOriginNeighbor = fine_grid->node[Second_Origin_Points[iNeighbor]]->GetPoint(iNode);
					if (iOriginNeighbor == Second_Origin_Points[jNeighbor]) {
						check_1 = false;
						break;
					}
				}
                
				if (check_1) {
					Suitable_Indirect_Neighbors->push_back(Second_Neighbor_Points[iNeighbor]);
                    
					/*--- Create alist with the suitable second neighbor, that we will use
					 to compute the third neighbors --*/
					Suitable_Second_Neighbors.push_back(Second_Neighbor_Points[iNeighbor]);
				}
			}
    
	vector<unsigned long>::iterator it;
    
	/*--- Remove repeated from the suitable second neighbors ---*/
	sort(Suitable_Second_Neighbors.begin(), Suitable_Second_Neighbors.end());
	it = unique(Suitable_Second_Neighbors.begin(), Suitable_Second_Neighbors.end());
	Suitable_Second_Neighbors.resize(it - Suitable_Second_Neighbors.begin());
    
	/*--- Remove repeated from first neighbors ---*/
	sort(First_Neighbor_Points.begin(), First_Neighbor_Points.end());
	it = unique(First_Neighbor_Points.begin(), First_Neighbor_Points.end());
	First_Neighbor_Points.resize(it - First_Neighbor_Points.begin());
    
	/*--- Create a list with the third neighbors, without first, second, and seed neighbors ---*/
	vector<unsigned long> Third_Neighbor_Points, Third_Origin_Points;
    
	for (jNode = 0; jNode <	Suitable_Second_Neighbors.size(); jNode ++) {
		kPoint = Suitable_Second_Neighbors[jNode];
        
		for (kNode = 0; kNode <	fine_grid->node[kPoint]->GetnPoint(); kNode ++) {
			lPoint = fine_grid->node[kPoint]->GetPoint(kNode);
            
			/*--- Check that the third neighbor do not belong to the first neighbors or the seed ---*/
			ThirdNeighborSeed = true;
            
			for (iNeighbor = 0; iNeighbor <	First_Neighbor_Points.size(); iNeighbor ++)
				if (lPoint == First_Neighbor_Points[iNeighbor]) {
					ThirdNeighborSeed = false;
					break;
				}
            
			/*--- Check that the third neighbor do not belong to the second neighbors ---*/
			for (iNeighbor = 0; iNeighbor <	Suitable_Second_Neighbors.size(); iNeighbor ++)
				if (lPoint == Suitable_Second_Neighbors[iNeighbor]) {
					ThirdNeighborSeed = false;
					break;
				}
            
			if (ThirdNeighborSeed) {
				Third_Neighbor_Points.push_back(lPoint);
				Third_Origin_Points.push_back(kPoint);
			}
            
		}
	}
    
	/*---  Identify those third neighbors that are repeated (candidate to be added) ---*/
	for (iNeighbor = 0; iNeighbor <	Third_Neighbor_Points.size(); iNeighbor ++)
		for (jNeighbor = 0; jNeighbor <	Third_Neighbor_Points.size(); jNeighbor ++)
            
        /*--- Repeated second neighbor with different origin ---*/
			if ((Third_Neighbor_Points[iNeighbor] == Third_Neighbor_Points[jNeighbor]) &&
                (Third_Origin_Points[iNeighbor] != Third_Origin_Points[jNeighbor]) &&
                (iNeighbor < jNeighbor)) {
                
				/*--- Check that the origin nodes are not neighbor ---*/
				check_1 = true;
				for (iNode = 0; iNode <	fine_grid->node[Third_Origin_Points[iNeighbor]]->GetnPoint(); iNode ++) {
					iOriginNeighbor = fine_grid->node[Third_Origin_Points[iNeighbor]]->GetPoint(iNode);
					if (iOriginNeighbor == Third_Origin_Points[jNeighbor]) {
						check_1 = false;
						break;
					}
				}
                
				if (check_1)
					Suitable_Indirect_Neighbors->push_back(Third_Neighbor_Points[iNeighbor]);
                
			}
    
	/*--- Remove repeated from Suitable Indirect Neighbors List ---*/
	sort(Suitable_Indirect_Neighbors->begin(), Suitable_Indirect_Neighbors->end());
	it = unique(Suitable_Indirect_Neighbors->begin(), Suitable_Indirect_Neighbors->end());
	Suitable_Indirect_Neighbors->resize(it - Suitable_Indirect_Neighbors->begin());
    
}



void CMultiGridGeometry::SetPsuP(CGeometry *fine_grid) {
	unsigned long iFinePoint, iFinePoint_Neighbor, iParent, iCoarsePoint;
	unsigned short iChildren, iNode;
    
	/*--- Set the point suronfding a point ---*/
	for (iCoarsePoint = 0; iCoarsePoint < nPoint; iCoarsePoint ++)
		for (iChildren = 0; iChildren <  node[iCoarsePoint]->GetnChildren_CV(); iChildren ++) {
			iFinePoint = node[iCoarsePoint]->GetChildren_CV(iChildren);
			for (iNode = 0; iNode < fine_grid->node[iFinePoint]->GetnPoint(); iNode ++) {
				iFinePoint_Neighbor = fine_grid->node[iFinePoint]->GetPoint(iNode);
				iParent = fine_grid->node[iFinePoint_Neighbor]->GetParent_CV();
				if (iParent != iCoarsePoint) node[iCoarsePoint]->SetPoint(iParent);
			}
		}
    
	/*--- Set the number of neighbors variable, this is
	 important for JST and multigrid in parallel ---*/
	for (iCoarsePoint = 0; iCoarsePoint < nPoint; iCoarsePoint ++)
		node[iCoarsePoint]->SetnNeighbor(node[iCoarsePoint]->GetnPoint());
    
}

void CMultiGridGeometry::SetVertex(CGeometry *fine_grid, CConfig *config) {
	unsigned long  iVertex, iFinePoint, iCoarsePoint;
	unsigned short iMarker, iMarker_Tag, iChildren;
    
	nMarker = fine_grid->GetnMarker();
    
	/*--- If any children node belong to the boundary then the entire control
	 volume will belong to the boundary ---*/
	for (iCoarsePoint = 0; iCoarsePoint < nPoint; iCoarsePoint ++)
		for (iChildren = 0; iChildren <	node[iCoarsePoint]->GetnChildren_CV(); iChildren ++) {
			iFinePoint = node[iCoarsePoint]->GetChildren_CV(iChildren);
			if (fine_grid->node[iFinePoint]->GetBoundary()) {
				node[iCoarsePoint]->SetBoundary(nMarker);
				break;
			}
		}
    
	vertex = new CVertex**[nMarker];
	nVertex = new unsigned long [nMarker];
    
	Tag_to_Marker = new string [MAX_INDEX_VALUE];
	for (iMarker_Tag = 0; iMarker_Tag < MAX_INDEX_VALUE; iMarker_Tag++)
		Tag_to_Marker[iMarker_Tag] = fine_grid->GetMarker_Tag(iMarker_Tag);
    
	/*--- Compute the number of vertices to do the dimensionalization ---*/
	for (iMarker = 0; iMarker < nMarker; iMarker++) nVertex[iMarker] = 0;
    
    
	for (iCoarsePoint = 0; iCoarsePoint < nPoint; iCoarsePoint ++) {
		if (node[iCoarsePoint]->GetBoundary()) {
			for (iChildren = 0; iChildren <	node[iCoarsePoint]->GetnChildren_CV(); iChildren ++) {
				iFinePoint = node[iCoarsePoint]->GetChildren_CV(iChildren);
				for (iMarker = 0; iMarker < nMarker; iMarker ++) {
					if ((fine_grid->node[iFinePoint]->GetVertex(iMarker) != -1) && (node[iCoarsePoint]->GetVertex(iMarker) == -1)) {
						iVertex = nVertex[iMarker];
						node[iCoarsePoint]->SetVertex(iVertex,iMarker);
						nVertex[iMarker]++;
					}
				}
			}
		}
	}
    
	for (iMarker = 0; iMarker < nMarker; iMarker++) {
		vertex[iMarker] = new CVertex* [fine_grid->GetnVertex(iMarker)+1];
		nVertex[iMarker] = 0;
	}
    
	for (iCoarsePoint = 0; iCoarsePoint < nPoint; iCoarsePoint ++)
		if (node[iCoarsePoint]->GetBoundary())
			for (iMarker = 0; iMarker < nMarker; iMarker ++)
				node[iCoarsePoint]->SetVertex(-1,iMarker);
    
	for (iMarker = 0; iMarker < nMarker; iMarker++) nVertex[iMarker] = 0;
    
	for (iCoarsePoint = 0; iCoarsePoint < nPoint; iCoarsePoint ++)
		if (node[iCoarsePoint]->GetBoundary())
			for (iChildren = 0; iChildren <	node[iCoarsePoint]->GetnChildren_CV(); iChildren ++) {
				iFinePoint = node[iCoarsePoint]->GetChildren_CV(iChildren);
				for (iMarker = 0; iMarker < fine_grid->GetnMarker(); iMarker ++) {
					if ((fine_grid->node[iFinePoint]->GetVertex(iMarker) != -1) && (node[iCoarsePoint]->GetVertex(iMarker) == -1)) {
						iVertex = nVertex[iMarker];
						vertex[iMarker][iVertex] = new CVertex(iCoarsePoint, nDim);
						node[iCoarsePoint]->SetVertex(iVertex,iMarker);
                        
						/*--- Set the transformation to apply ---*/
						unsigned long ChildVertex = fine_grid->node[iFinePoint]->GetVertex(iMarker);
						unsigned short RotationKind = fine_grid->vertex[iMarker][ChildVertex]->GetRotation_Type();
						unsigned short MatchingZone = fine_grid->vertex[iMarker][ChildVertex]->GetMatching_Zone();
						vertex[iMarker][iVertex]->SetRotation_Type(RotationKind);
						vertex[iMarker][iVertex]->SetMatching_Zone(MatchingZone);
						nVertex[iMarker]++;
					}
				}
			}
}

void CMultiGridGeometry::MatchNearField(CConfig *config) {
    
	unsigned short iMarker;
	unsigned long iVertex, iPoint;
    
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if (config->GetMarker_All_Boundary(iMarker) == NEARFIELD_BOUNDARY)
			for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
				iPoint = vertex[iMarker][iVertex]->GetNode();
				vertex[iMarker][iVertex]->SetDonorPoint(iPoint);
			}
    
}

void CMultiGridGeometry::MatchInterface(CConfig *config) {
    
	unsigned short iMarker;
	unsigned long iVertex, iPoint;
    
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if (config->GetMarker_All_Boundary(iMarker) == INTERFACE_BOUNDARY)
			for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
				iPoint = vertex[iMarker][iVertex]->GetNode();
				vertex[iMarker][iVertex]->SetDonorPoint(iPoint);
			}
    
}


void CMultiGridGeometry::SetControlVolume(CConfig *config, CGeometry *fine_grid, unsigned short action) {
    
	unsigned long iFinePoint,iFinePoint_Neighbor, iCoarsePoint, iEdge, iParent, iPoint, jPoint;
    long FineEdge, CoarseEdge;
	unsigned short iChildren, iNode, iDim;
	bool change_face_orientation;
	double *Normal, Coarse_Volume, Area, *NormalFace = NULL;
	Normal = new double [nDim];
    
	/*--- Compute the area of the coarse volume ---*/
	for (iCoarsePoint = 0; iCoarsePoint < nPoint; iCoarsePoint ++) {
		node[iCoarsePoint]->SetVolume(0.0);
		Coarse_Volume = 0.0;
		for (iChildren = 0; iChildren < node[iCoarsePoint]->GetnChildren_CV(); iChildren ++) {
			iFinePoint = node[iCoarsePoint]->GetChildren_CV(iChildren);
			Coarse_Volume += fine_grid->node[iFinePoint]->GetVolume();
		}
		node[iCoarsePoint]->SetVolume(Coarse_Volume);
	}
    
	/*--- Update or not the values of faces at the edge ---*/
	if (action != ALLOCATE) {
		for(iEdge=0; iEdge < nEdge; iEdge++)
			edge[iEdge]->SetZeroValues();
	}
    
	for (iCoarsePoint = 0; iCoarsePoint < nPoint; iCoarsePoint ++)
		for (iChildren = 0; iChildren < node[iCoarsePoint]->GetnChildren_CV(); iChildren ++) {
			iFinePoint = node[iCoarsePoint]->GetChildren_CV(iChildren);
            
			for (iNode = 0; iNode < fine_grid->node[iFinePoint]->GetnPoint(); iNode ++) {
				iFinePoint_Neighbor = fine_grid->node[iFinePoint]->GetPoint(iNode);
				iParent = fine_grid->node[iFinePoint_Neighbor]->GetParent_CV();
				if ((iParent != iCoarsePoint) && (iParent < iCoarsePoint)) {
                    
					FineEdge = fine_grid->FindEdge(iFinePoint, iFinePoint_Neighbor);
                    
					change_face_orientation = false;
					if (iFinePoint < iFinePoint_Neighbor) change_face_orientation = true;
                    
					CoarseEdge = FindEdge(iParent, iCoarsePoint);
                    
					fine_grid->edge[FineEdge]->GetNormal(Normal);
                    
					if (change_face_orientation) {
						for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
						edge[CoarseEdge]->AddNormal(Normal);
					}
					else {
						edge[CoarseEdge]->AddNormal(Normal);
					}
				}
			}
		}
	delete[] Normal;
    
	/*--- Check if there isn't any element with only one neighbor...
	 a CV that is inside another CV ---*/
	for (iPoint = 0; iPoint < nPoint; iPoint ++) {
		if (node[iPoint]->GetnPoint() == 1) {
			jPoint = node[iPoint]->GetPoint(0);
			node[jPoint]->AddVolume(node[iPoint]->GetVolume());
		}
	}
    
    /*--- Check if there is a normal with null area ---*/
    for (iEdge = 0; iEdge < nEdge; iEdge++) {
        NormalFace = edge[iEdge]->GetNormal();
        Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += NormalFace[iDim]*NormalFace[iDim];
        Area = sqrt(Area);
        if (Area == 0.0) for (iDim = 0; iDim < nDim; iDim++) NormalFace[iDim] = EPS*EPS;
    }
    
}

void CMultiGridGeometry::SetBoundControlVolume(CConfig *config, CGeometry *fine_grid, unsigned short action) {
	unsigned long iCoarsePoint, iFinePoint, FineVertex, iVertex;
	unsigned short iMarker, iChildren, iDim;
	double *Normal, Area, *NormalFace = NULL;
    
	Normal = new double [nDim];
    
	if (action != ALLOCATE) {
		for (iMarker = 0; iMarker < nMarker; iMarker++)
			for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
				vertex[iMarker][iVertex]->SetZeroValues();
	}
    
	for (iMarker = 0; iMarker < nMarker; iMarker ++)
		for(iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
			iCoarsePoint = vertex[iMarker][iVertex]->GetNode();
			for (iChildren = 0; iChildren < node[iCoarsePoint]->GetnChildren_CV(); iChildren ++) {
				iFinePoint = node[iCoarsePoint]->GetChildren_CV(iChildren);
				if (fine_grid->node[iFinePoint]->GetVertex(iMarker)!=-1) {
					FineVertex = fine_grid->node[iFinePoint]->GetVertex(iMarker);
					fine_grid->vertex[iMarker][FineVertex]->GetNormal(Normal);
					vertex[iMarker][iVertex]->AddNormal(Normal);
				}
			}
		}
    
	delete[] Normal;
    
    /*--- Check if there is a normal with null area ---*/
	for (iMarker = 0; iMarker < nMarker; iMarker ++)
		for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
            NormalFace = vertex[iMarker][iVertex]->GetNormal();
            Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += NormalFace[iDim]*NormalFace[iDim];
            Area = sqrt(Area);
            if (Area == 0.0) for (iDim = 0; iDim < nDim; iDim++) NormalFace[iDim] = EPS*EPS;
        }
    
}

void CMultiGridGeometry::SetCoord(CGeometry *geometry) {
	unsigned long Point_Fine, Point_Coarse;
	unsigned short iChildren, iDim;
	double Area_Parent, Area_Children;
	double *Coordinates_Fine, *Coordinates;
	Coordinates = new double[nDim];
    
	for (Point_Coarse = 0; Point_Coarse < GetnPoint(); Point_Coarse++) {
		Area_Parent = node[Point_Coarse]->GetVolume();
		for (iDim = 0; iDim < nDim; iDim++) Coordinates[iDim] = 0.0;
		for (iChildren = 0; iChildren < node[Point_Coarse]->GetnChildren_CV(); iChildren++) {
			Point_Fine = node[Point_Coarse]->GetChildren_CV(iChildren);
			Area_Children = geometry->node[Point_Fine]->GetVolume();
			Coordinates_Fine = geometry->node[Point_Fine]->GetCoord();
			for (iDim = 0; iDim < nDim; iDim++)
				Coordinates[iDim] += Coordinates_Fine[iDim]*Area_Children/Area_Parent;
		}
		for (iDim = 0; iDim < nDim; iDim++)
			node[Point_Coarse]->SetCoord(iDim,Coordinates[iDim]);
	}
	delete[] Coordinates;
}

void CMultiGridGeometry::SetRotationalVelocity(CConfig *config) {
    
	unsigned long iPoint_Coarse;
	double RotVel[3], Distance[3], *Coord, Center[3], Omega[3], L_Ref;
    
    /*--- Center of rotation & angular velocity vector from config. ---*/
    
    Center[0] = config->GetMotion_Origin_X(ZONE_0);
    Center[1] = config->GetMotion_Origin_Y(ZONE_0);
    Center[2] = config->GetMotion_Origin_Z(ZONE_0);
    Omega[0]  = config->GetRotation_Rate_X(ZONE_0)/config->GetOmega_Ref();
    Omega[1]  = config->GetRotation_Rate_Y(ZONE_0)/config->GetOmega_Ref();
    Omega[2]  = config->GetRotation_Rate_Z(ZONE_0)/config->GetOmega_Ref();
    L_Ref     = config->GetLength_Ref();
    
	/*--- Loop over all nodes and set the rotational velocity. ---*/
    
	for (iPoint_Coarse = 0; iPoint_Coarse < GetnPoint(); iPoint_Coarse++) {
        
		/*--- Get the coordinates of the current node ---*/
        
		Coord = node[iPoint_Coarse]->GetCoord();
        
		/*--- Calculate the non-dim. distance from the rotation center ---*/
        
		Distance[0] = (Coord[0]-Center[0])/L_Ref;
		Distance[1] = (Coord[1]-Center[1])/L_Ref;
		Distance[2] = (Coord[2]-Center[2])/L_Ref;
        
		/*--- Calculate the angular velocity as omega X r ---*/
        
		RotVel[0] = Omega[1]*(Distance[2]) - Omega[2]*(Distance[1]);
		RotVel[1] = Omega[2]*(Distance[0]) - Omega[0]*(Distance[2]);
		RotVel[2] = Omega[0]*(Distance[1]) - Omega[1]*(Distance[0]);
        
        /*--- Store the grid velocity at this node ---*/
        
		node[iPoint_Coarse]->SetGridVel(RotVel);
        
	}
    
}

void CMultiGridGeometry::SetGridVelocity(CConfig *config, unsigned long iter) {
    
	/*--- Local variables ---*/
    
	double *Coord_nP1 = NULL, *Coord_n = NULL, *Coord_nM1 = NULL;
    double TimeStep, GridVel = 0.0;
	unsigned long Point_Coarse;
	unsigned short iDim;
    
	/*--- Compute the velocity of each node in the volume mesh ---*/
    
	for (Point_Coarse = 0; Point_Coarse < GetnPoint(); Point_Coarse++) {
        
		/*--- Coordinates of the current point at n+1, n, & n-1 time levels ---*/
        
		Coord_nM1 = node[Point_Coarse]->GetCoord_n1();
		Coord_n   = node[Point_Coarse]->GetCoord_n();
		Coord_nP1 = node[Point_Coarse]->GetCoord();
        
		/*--- Unsteady time step ---*/
        
		TimeStep = config->GetDelta_UnstTimeND();
        
		/*--- Compute mesh velocity with 1st or 2nd-order approximation ---*/
        
		for(iDim = 0; iDim < nDim; iDim++) {
			if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
				GridVel = ( Coord_nP1[iDim] - Coord_n[iDim] ) / TimeStep;
			if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
				GridVel = ( 3.0*Coord_nP1[iDim] - 4.0*Coord_n[iDim]
                           +  1.0*Coord_nM1[iDim] ) / (2.0*TimeStep);
            
			/*--- Store grid velocity for this point ---*/
            
			node[Point_Coarse]->SetGridVel(iDim,GridVel);
            
		}
	}
}

void CMultiGridGeometry::SetRestricted_GridVelocity(CGeometry *fine_mesh, CConfig *config) {
    
	/*--- Local variables ---*/
	unsigned short iDim, iChild;
	unsigned long Point_Coarse, Point_Fine;
	double Area_Parent, Area_Child, Grid_Vel[3], *Grid_Vel_Fine;
    
	/*--- Loop over all coarse mesh points ---*/
	for (Point_Coarse = 0; Point_Coarse < GetnPoint(); Point_Coarse++) {
		Area_Parent = node[Point_Coarse]->GetVolume();
        
		/*--- Zero out the grid velocity ---*/
		for (iDim = 0; iDim < nDim; iDim++)
			Grid_Vel[iDim] = 0.0;
        
		/*--- Loop over all of the children for this coarse CV and compute
         a grid velocity based on the values in the child CVs (fine mesh). ---*/
		for (iChild = 0; iChild < node[Point_Coarse]->GetnChildren_CV(); iChild++) {
			Point_Fine    = node[Point_Coarse]->GetChildren_CV(iChild);
			Area_Child    = fine_mesh->node[Point_Fine]->GetVolume();
			Grid_Vel_Fine = fine_mesh->node[Point_Fine]->GetGridVel();
			for (iDim = 0; iDim < nDim; iDim++)
				Grid_Vel[iDim] += Grid_Vel_Fine[iDim]*Area_Child/Area_Parent;
		}
        
		/*--- Set the grid velocity for this coarse node. ---*/
		for (iDim = 0; iDim < nDim; iDim++)
			node[Point_Coarse]->SetGridVel(iDim, Grid_Vel[iDim]);
	}
}


void CMultiGridGeometry::FindNormal_Neighbor(CConfig *config) {
    
	unsigned short iMarker, iDim;
	unsigned long iPoint, iVertex;
    
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
        
		if (config->GetMarker_All_Boundary(iMarker) != SEND_RECEIVE &&
            config->GetMarker_All_Boundary(iMarker) != INTERFACE_BOUNDARY &&
            config->GetMarker_All_Boundary(iMarker) != NEARFIELD_BOUNDARY ) {
            
			for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
                
				iPoint = vertex[iMarker][iVertex]->GetNode();
                
				/*--- If the node belong to the domain ---*/
				if (node[iPoint]->GetDomain()) {
                    
					/*--- Compute closest normal neighbor ---*/
					double cos_max, scalar_prod, norm_vect, norm_Normal, cos_alpha, diff_coord;
					unsigned long Point_Normal = 0, jPoint;
					unsigned short iNeigh;
					double *Normal = vertex[iMarker][iVertex]->GetNormal();
					cos_max = -1.0;
					for (iNeigh = 0; iNeigh < node[iPoint]->GetnPoint(); iNeigh++) {
						jPoint = node[iPoint]->GetPoint(iNeigh);
						scalar_prod = 0.0; norm_vect = 0.0; norm_Normal = 0.0;
						for(iDim = 0; iDim < nDim; iDim++) {
							diff_coord = node[jPoint]->GetCoord(iDim)-node[iPoint]->GetCoord(iDim);
							scalar_prod += diff_coord*Normal[iDim];
							norm_vect += diff_coord*diff_coord;
							norm_Normal += Normal[iDim]*Normal[iDim];
						}
						norm_vect = sqrt(norm_vect);
						norm_Normal = sqrt(norm_Normal);
						cos_alpha = scalar_prod/(norm_vect*norm_Normal);
                        
						/*--- Get maximum cosine (not minimum because normals are oriented inwards) ---*/
						if (cos_alpha >= cos_max) {
							Point_Normal = jPoint;
							cos_max = cos_alpha;
						}
					}
					vertex[iMarker][iVertex]->SetNormal_Neighbor(Point_Normal);
				}
			}
        }
    }
}


void CMultiGridGeometry::SetGeometryPlanes(CConfig *config) {
	bool loop_on;
	unsigned short iMarker = 0;
	double auxXCoord, auxYCoord, auxZCoord,	*Face_Normal = NULL, auxArea, *Xcoord = NULL, *Ycoord = NULL, *Zcoord = NULL, *FaceArea = NULL;
	unsigned long jVertex, iVertex,ixCoord, iPoint, iVertex_Wall, nVertex_Wall = 0;
    
	/*--- Compute the total number of points on the near-field ---*/
	nVertex_Wall = 0;
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if ((config->GetMarker_All_Boundary(iMarker) == HEAT_FLUX) ||
            (config->GetMarker_All_Boundary(iMarker) == ISOTHERMAL) ||
            (config->GetMarker_All_Boundary(iMarker) == EULER_WALL))
			nVertex_Wall += nVertex[iMarker];
    
    
	/*--- Create an array with all the coordinates, points, pressures, face area,
	 equivalent area, and nearfield weight ---*/
	Xcoord = new double[nVertex_Wall];
	Ycoord = new double[nVertex_Wall];
	if (nDim == 3)	Zcoord = new double[nVertex_Wall];
	FaceArea = new double[nVertex_Wall];
    
	/*--- Copy the boundary information to an array ---*/
	iVertex_Wall = 0;
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if ((config->GetMarker_All_Boundary(iMarker) == HEAT_FLUX) ||
            (config->GetMarker_All_Boundary(iMarker) == ISOTHERMAL) ||
            (config->GetMarker_All_Boundary(iMarker) == EULER_WALL))
			for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
				iPoint = vertex[iMarker][iVertex]->GetNode();
				Xcoord[iVertex_Wall] = node[iPoint]->GetCoord(0);
				Ycoord[iVertex_Wall] = node[iPoint]->GetCoord(1);
				if (nDim==3) Zcoord[iVertex_Wall] = node[iPoint]->GetCoord(2);
				Face_Normal = vertex[iMarker][iVertex]->GetNormal();
				FaceArea[iVertex_Wall] = fabs(Face_Normal[nDim-1]);
				iVertex_Wall ++;
			}
    
    
	//vector<double> XCoordList;
	vector<double>::iterator IterXCoordList;
    
	for (iVertex = 0; iVertex < nVertex_Wall; iVertex++)
		XCoordList.push_back(Xcoord[iVertex]);
    
	sort( XCoordList.begin(), XCoordList.end());
	IterXCoordList = unique( XCoordList.begin(), XCoordList.end());
	XCoordList.resize( IterXCoordList - XCoordList.begin() );
    
	/*--- Create vectors and distribute the values among the different PhiAngle queues ---*/
	Xcoord_plane.resize(XCoordList.size());
	Ycoord_plane.resize(XCoordList.size());
	if (nDim==3) Zcoord_plane.resize(XCoordList.size());
	FaceArea_plane.resize(XCoordList.size());
	Plane_points.resize(XCoordList.size());
    
    
	double dist_ratio;
	unsigned long iCoord;
    
	/*--- Distribute the values among the different PhiAngles ---*/
	for (iPoint = 0; iPoint < nPoint; iPoint++) {
		if (node[iPoint]->GetDomain()) {
			loop_on = true;
			for (ixCoord = 0; ixCoord < XCoordList.size()-1 && loop_on; ixCoord++) {
				dist_ratio = (node[iPoint]->GetCoord(0) - XCoordList[ixCoord])/(XCoordList[ixCoord+1]- XCoordList[ixCoord]);
				if (dist_ratio >= 0 && dist_ratio <= 1.0) {
					if (dist_ratio <= 0.5) iCoord = ixCoord;
					else iCoord = ixCoord+1;
					Xcoord_plane[iCoord].push_back(node[iPoint]->GetCoord(0) );
					Ycoord_plane[iCoord].push_back(node[iPoint]->GetCoord(1) );
					if (nDim==3) Zcoord_plane[iCoord].push_back(node[iPoint]->GetCoord(2) );
					FaceArea_plane[iCoord].push_back(node[iPoint]->GetVolume());   ///// CHECK AREA CALCULATION
					Plane_points[iCoord].push_back(iPoint );
					loop_on = false;
				}
			}
		}
	}
    
	unsigned long auxPoint;
	/*--- Order the arrays in ascending values of y ---*/
	for (ixCoord = 0; ixCoord < XCoordList.size(); ixCoord++)
		for (iVertex = 0; iVertex < Xcoord_plane[ixCoord].size(); iVertex++)
			for (jVertex = 0; jVertex < Xcoord_plane[ixCoord].size() - 1 - iVertex; jVertex++)
				if (Ycoord_plane[ixCoord][jVertex] > Ycoord_plane[ixCoord][jVertex+1]) {
					auxXCoord = Xcoord_plane[ixCoord][jVertex]; Xcoord_plane[ixCoord][jVertex] = Xcoord_plane[ixCoord][jVertex+1]; Xcoord_plane[ixCoord][jVertex+1] = auxXCoord;
					auxYCoord = Ycoord_plane[ixCoord][jVertex]; Ycoord_plane[ixCoord][jVertex] = Ycoord_plane[ixCoord][jVertex+1]; Ycoord_plane[ixCoord][jVertex+1] = auxYCoord;
					auxPoint = Plane_points[ixCoord][jVertex]; Plane_points[ixCoord][jVertex] = Plane_points[ixCoord][jVertex+1]; Plane_points[ixCoord][jVertex+1] = auxPoint;
					if (nDim==3) {
						auxZCoord = Zcoord_plane[ixCoord][jVertex]; Zcoord_plane[ixCoord][jVertex] = Zcoord_plane[ixCoord][jVertex+1]; Zcoord_plane[ixCoord][jVertex+1] = auxZCoord;
					}
					auxArea = FaceArea_plane[ixCoord][jVertex]; FaceArea_plane[ixCoord][jVertex] = FaceArea_plane[ixCoord][jVertex+1]; FaceArea_plane[ixCoord][jVertex+1] = auxArea;
				}
    
	/*--- Delete structures ---*/
	delete[] Xcoord; delete[] Ycoord;
	if (nDim==3) delete[] Zcoord;
	delete[] FaceArea;
}


CBoundaryGeometry::CBoundaryGeometry(CConfig *config, string val_mesh_filename, unsigned short val_format) : CGeometry() {
    
	string text_line;
	ifstream mesh_file;
	unsigned short iNode_Surface, VTK_Type, iMarker, iChar, iCount = 0, val_iZone = 1, val_nZone = 1;
	unsigned long Point_Surface, iElem_Surface, iElem_Bound = 0, iPoint = 0, iElem = 0, ielem = 0,
    nelem_edge = 0, nelem_triangle = 0, nelem_quad = 0,
    vnodes_edge[2], vnodes_triangle[3],
    vnodes_quad[4], dummy, GlobalIndex;
	string Marker_Tag;
	char cstr[200];
	int rank = MASTER_NODE, size = SINGLE_NODE;
	bool domain_flag = false;
	bool found_transform = false;
	nZone = val_nZone;
	double Coord_2D[2], Coord_3D[3];
	string::size_type position;
    
	Global_nPointDomain = 0;
	FinestMGLevel = true;
    
	if(rank == MASTER_NODE)
		cout << endl <<"---------------------- Read grid file information -----------------------" << endl;
    
    
	/*--- Determine whether there are multiplze zones first ---*/
	strcpy (cstr, val_mesh_filename.c_str());
	mesh_file.open(cstr, ios::in);
	if (mesh_file.fail()) {
		cout << "There is no geometry file (CBoundaryGeometry)!" << endl;
		cout << "Press any key to exit..." << endl;
		cin.get();
		exit(1);
	}
    
	/*--- If more than one, find the domain in the mesh file ---*/
	if (val_nZone > 1) {
		while (getline (mesh_file,text_line)) {
			/*--- Search for the current domain ---*/
			position = text_line.find ("IZONE=",0);
			if (position != string::npos) {
				text_line.erase (0,6);
				unsigned short jDomain = atoi(text_line.c_str());
				if (jDomain == val_iZone) {
					if (rank == MASTER_NODE) cout << "Reading zone " << val_iZone << ":" << endl;
					break;
				}
			}
		}
	}
    
	/*--- Read grid file with format SU2 ---*/
	while (getline (mesh_file,text_line)) {
        
		/*--- Read the dimension of the problem ---*/
		position = text_line.find ("NDIME=",0);
		if (position != string::npos) {
			if (domain_flag == false) {
				text_line.erase (0,6); nDim = atoi(text_line.c_str());
				if (rank == MASTER_NODE) {
					if (nDim == 2) cout << "Two dimensional problem." << endl;
					if (nDim == 3) cout << "Three dimensional problem." << endl;
				}
				domain_flag = true;
			} else {
				break;
			}
		}
        
		/*--- Read the information about inner elements ---*/
		position = text_line.find ("NELEM=",0);
		if (position != string::npos) {
			text_line.erase (0,6); nElem = atoi(text_line.c_str());
			if (size == 1)
				cout << nElem << " interior elements. ";
            
			while (iElem < nElem) {
				getline(mesh_file,text_line);
				iElem++;
			}
		}
        
		/*--- Read number of points ---*/
		position = text_line.find ("NPOIN=",0);
		if (position != string::npos) {
			text_line.erase (0,6);
            
			/*--- Check for ghost points. ---*/
			stringstream test_line(text_line);
			while (test_line >> dummy)
				iCount++;
            
			/*--- Now read and store the number of points and possible ghost points. ---*/
			stringstream  stream_line(text_line);
			if (iCount == 2) {
				stream_line >> nPoint;
				stream_line >> nPointDomain;
				if (size == 1)
					cout << nPoint << " points, and " << nPoint-nPointDomain << " ghost points." << endl;
                
				/*--- Set some important point information for parallel simulations. ---*/
				Global_nPointDomain = nPointDomain;
			}
			else if (iCount == 1) {
				stream_line >> nPoint;
				nPointDomain = nPoint;
				Global_nPointDomain = nPoint;
				if (rank == MASTER_NODE) cout << nPoint << " points." << endl;
			}
			else {
				cout << "NPOIN improperly specified!!" << endl;
				cout << "Press any key to exit..." << endl;
				cin.get();
				exit(1);
			}
            
			node = new CPoint*[nPoint];
			while (iPoint < nPoint) {
				getline(mesh_file,text_line);
				istringstream point_line(text_line);
				switch(nDim) {
                    case 2:
                        GlobalIndex = iPoint;
                        point_line >> Coord_2D[0]; point_line >> Coord_2D[1];
                        node[iPoint] = new CPoint(Coord_2D[0], Coord_2D[1], GlobalIndex, config);
                        iPoint++; break;
                    case 3:
                        GlobalIndex = iPoint;
                        point_line >> Coord_3D[0]; point_line >> Coord_3D[1]; point_line >> Coord_3D[2];
                        node[iPoint] = new CPoint(Coord_3D[0], Coord_3D[1], Coord_3D[2], GlobalIndex, config);
                        iPoint++; break;
				}
			}
		}
        
		/*--- Read number of markers ---*/
		position = text_line.find ("NMARK=",0);
		if (position != string::npos) {
			text_line.erase (0,6); nMarker = atoi(text_line.c_str());
			if (size == 1) cout << nMarker << " surface markers." << endl;
			config->SetnMarker_All(nMarker);
			bound = new CPrimalGrid**[nMarker];
			nElem_Bound = new unsigned long [nMarker];
			Tag_to_Marker = new string [MAX_INDEX_VALUE];
            
			for (iMarker = 0 ; iMarker < nMarker; iMarker++) {
				getline (mesh_file,text_line);
				text_line.erase (0,11);
				string::size_type position;
				for (iChar = 0; iChar < 20; iChar++) {
					position = text_line.find( " ", 0 );
					if(position != string::npos) text_line.erase (position,1);
					position = text_line.find( "\r", 0 );
					if(position != string::npos) text_line.erase (position,1);
					position = text_line.find( "\n", 0 );
					if(position != string::npos) text_line.erase (position,1);
				}
				Marker_Tag = text_line.c_str();
                
				/*--- Physical boundaries definition ---*/
				if (Marker_Tag != "SEND_RECEIVE") {
					getline (mesh_file,text_line);
					text_line.erase (0,13); nElem_Bound[iMarker] = atoi(text_line.c_str());
					if (size == 1)
						cout << nElem_Bound[iMarker]  << " boundary elements in index "<< iMarker <<" (Marker = " <<Marker_Tag<< ")." << endl;
					bound[iMarker] = new CPrimalGrid* [nElem_Bound[iMarker]];
                    
					nelem_edge = 0; nelem_triangle = 0; nelem_quad = 0; ielem = 0;
					for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
						getline(mesh_file,text_line);
						istringstream bound_line(text_line);
						bound_line >> VTK_Type;
						switch(VTK_Type) {
                            case LINE:
                                bound_line >> vnodes_edge[0]; bound_line >> vnodes_edge[1];
                                bound[iMarker][ielem] = new CLine(vnodes_edge[0],vnodes_edge[1],2);
                                ielem++; nelem_edge++; break;
                            case TRIANGLE:
                                bound_line >> vnodes_triangle[0]; bound_line >> vnodes_triangle[1]; bound_line >> vnodes_triangle[2];
                                bound[iMarker][ielem] = new CTriangle(vnodes_triangle[0],vnodes_triangle[1],vnodes_triangle[2],3);
                                ielem++; nelem_triangle++; break;
                            case RECTANGLE:
                                bound_line >> vnodes_quad[0]; bound_line >> vnodes_quad[1]; bound_line >> vnodes_quad[2]; bound_line >> vnodes_quad[3];
                                bound[iMarker][ielem] = new CRectangle(vnodes_quad[0],vnodes_quad[1],vnodes_quad[2],vnodes_quad[3],3);
                                ielem++; nelem_quad++; break;
						}
					}
                    
					/*--- Update config information storing the boundary information in the right place ---*/
					Tag_to_Marker[config->GetMarker_Config_Tag(Marker_Tag)] = Marker_Tag;
					config->SetMarker_All_Tag(iMarker, Marker_Tag);
					config->SetMarker_All_Boundary(iMarker, config->GetMarker_Config_Boundary(Marker_Tag));
					config->SetMarker_All_Monitoring(iMarker, config->GetMarker_Config_Monitoring(Marker_Tag));
					config->SetMarker_All_Designing(iMarker, config->GetMarker_Config_Designing(Marker_Tag));
					config->SetMarker_All_Plotting(iMarker, config->GetMarker_Config_Plotting(Marker_Tag));
					config->SetMarker_All_DV(iMarker, config->GetMarker_Config_DV(Marker_Tag));
                    config->SetMarker_All_Moving(iMarker, config->GetMarker_Config_Moving(Marker_Tag));
					config->SetMarker_All_PerBound(iMarker, config->GetMarker_Config_PerBound(Marker_Tag));
					config->SetMarker_All_SendRecv(iMarker, NONE);
                    
				}
                
				/*--- Send-Receive boundaries definition ---*/
				else {
					unsigned long nelem_vertex = 0, vnodes_vertex;
					unsigned short transform;
					getline (mesh_file,text_line);
					text_line.erase (0,13); nElem_Bound[iMarker] = atoi(text_line.c_str());
					bound[iMarker] = new CPrimalGrid* [nElem_Bound[iMarker]];
                    
					nelem_vertex = 0; ielem = 0;
					getline (mesh_file,text_line); text_line.erase (0,8);
					config->SetMarker_All_Boundary(iMarker, SEND_RECEIVE);
					config->SetMarker_All_SendRecv(iMarker, atoi(text_line.c_str()));
                    
					for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
						getline(mesh_file,text_line);
						istringstream bound_line(text_line);
						bound_line >> VTK_Type; bound_line >> vnodes_vertex; bound_line >> transform;
                        
						bound[iMarker][ielem] = new CVertexMPI(vnodes_vertex, nDim);
						bound[iMarker][ielem]->SetRotation_Type(transform);
						ielem++; nelem_vertex++;
						if (config->GetMarker_All_SendRecv(iMarker) < 0)
							node[vnodes_vertex]->SetDomain(false);
                        
					}
                    
				}
                
			}
		}
        
		/*--- Read periodic transformation info (center, rotation, translation) ---*/
		position = text_line.find ("NPERIODIC=",0);
		if (position != string::npos) {
			unsigned short nPeriodic, iPeriodic, iIndex;
            
			/*--- Set bool signifying that periodic transormations were found ---*/
			found_transform = true;
            
			/*--- Read and store the number of transformations. ---*/
			text_line.erase (0,10); nPeriodic = atoi(text_line.c_str());
			if (rank == MASTER_NODE) {
				if (nPeriodic - 1 != 0)
					cout << nPeriodic - 1 << " periodic transformations." << endl;
			}
			config->SetnPeriodicIndex(nPeriodic);
            
			/*--- Store center, rotation, & translation in that order for each. ---*/
			for (iPeriodic = 0; iPeriodic < nPeriodic; iPeriodic++) {
				getline (mesh_file,text_line);
				position = text_line.find ("PERIODIC_INDEX=",0);
				if (position != string::npos) {
					text_line.erase (0,15); iIndex = atoi(text_line.c_str());
					if (iIndex != iPeriodic) {
						cout << "PERIODIC_INDEX out of order in SU2 file!!" << endl;
						cout << "Press any key to exit..." << endl;
						cin.get();
						exit(1);
					}
				}
				double* center    = new double[3];
				double* rotation  = new double[3];
				double* translate = new double[3];
				getline (mesh_file,text_line);
				istringstream cent(text_line);
				cent >> center[0]; cent >> center[1]; cent >> center[2];
				config->SetPeriodicCenter(iPeriodic, center);
				getline (mesh_file,text_line);
				istringstream rot(text_line);
				rot >> rotation[0]; rot >> rotation[1]; rot >> rotation[2];
				config->SetPeriodicRotation(iPeriodic, rotation);
				getline (mesh_file,text_line);
				istringstream tran(text_line);
				tran >> translate[0]; tran >> translate[1]; tran >> translate[2];
				config->SetPeriodicTranslate(iPeriodic, translate);
			}
            
		}
        
	}
    
	/*--- If no periodic transormations were found, store default zeros ---*/
	if (!found_transform) {
		unsigned short nPeriodic = 1, iPeriodic = 0;
		config->SetnPeriodicIndex(nPeriodic);
		double* center    = new double[3];
		double* rotation  = new double[3];
		double* translate = new double[3];
		for (unsigned short iDim = 0; iDim < 3; iDim++) {
			center[iDim] = 0.0; rotation[iDim] = 0.0; translate[iDim] = 0.0;
		}
		config->SetPeriodicCenter(iPeriodic, center);
		config->SetPeriodicRotation(iPeriodic, rotation);
		config->SetPeriodicTranslate(iPeriodic, translate);
	}
    
	/*--- Close the input file ---*/
	mesh_file.close();
    
	/*--- Loop over the surface element to set the boundaries ---*/
	for (iMarker = 0; iMarker < nMarker; iMarker++)
		for (iElem_Surface = 0; iElem_Surface < nElem_Bound[iMarker]; iElem_Surface++)
			for (iNode_Surface = 0; iNode_Surface < bound[iMarker][iElem_Surface]->GetnNodes(); iNode_Surface++) {
				Point_Surface = bound[iMarker][iElem_Surface]->GetNode(iNode_Surface);
				node[Point_Surface]->SetBoundary(nMarker);
                if (config->GetMarker_All_Boundary(iMarker) != SEND_RECEIVE &&
                    config->GetMarker_All_Boundary(iMarker) != INTERFACE_BOUNDARY &&
                    config->GetMarker_All_Boundary(iMarker) != NEARFIELD_BOUNDARY &&
                    config->GetMarker_All_Boundary(iMarker) != PERIODIC_BOUNDARY)
                    node[Point_Surface]->SetPhysicalBoundary(true);
			}
    
}


CBoundaryGeometry::~CBoundaryGeometry(void) {
    
}

void CBoundaryGeometry::SetVertex(void) {
	unsigned long  iPoint, iVertex, iElem;
	unsigned short iMarker, iNode;
    
	/*--- Initialize the Vertex vector for each node of the grid ---*/
	for (iPoint = 0; iPoint < nPoint; iPoint++)
		for (iMarker = 0; iMarker < nMarker; iMarker++)
			node[iPoint]->SetVertex(-1,iMarker);
    
	/*--- Create and compute the vector with the number of vertex per marker ---*/
	nVertex = new unsigned long [nMarker];
	for (iMarker = 0; iMarker < nMarker; iMarker++) {
		/*--- Initialize the number of Bound Vertex for each Marker ---*/
		nVertex[iMarker] = 0;
		for (iElem = 0; iElem < nElem_Bound[iMarker]; iElem++)
			for(iNode = 0; iNode < bound[iMarker][iElem]->GetnNodes(); iNode++) {
				iPoint = bound[iMarker][iElem]->GetNode(iNode);
				/*--- Set the vertex in the node information ---*/
				if (node[iPoint]->GetVertex(iMarker) == -1) {
					iVertex = nVertex[iMarker];
					node[iPoint]->SetVertex(nVertex[iMarker],iMarker);
					nVertex[iMarker]++;
				}
			}
	}
    
	/*--- Initialize the Vertex vector for each node, the previous result is deleted ---*/
	for (iPoint = 0; iPoint < nPoint; iPoint++)
		for (iMarker = 0; iMarker < nMarker; iMarker++)
			node[iPoint]->SetVertex(-1,iMarker);
    
	/*--- Create the bound vertex structure, note that the order
	 is the same as in the input file, this is important for Send/Receive part ---*/
	vertex = new CVertex**[nMarker];
	for (iMarker = 0; iMarker < nMarker; iMarker++) {
		vertex[iMarker] = new CVertex* [nVertex[iMarker]];
		nVertex[iMarker] = 0;
		/*--- Initialize the number of Bound Vertex for each Marker ---*/
		for (iElem = 0; iElem < nElem_Bound[iMarker]; iElem++)
			for(iNode = 0; iNode < bound[iMarker][iElem]->GetnNodes(); iNode++) {
				iPoint = bound[iMarker][iElem]->GetNode(iNode);
				/*--- Set the vertex in the node information ---*/
				if (node[iPoint]->GetVertex(iMarker) == -1) {
					iVertex = nVertex[iMarker];
					vertex[iMarker][iVertex] = new CVertex(iPoint, nDim);
					node[iPoint]->SetVertex(nVertex[iMarker],iMarker);
					nVertex[iMarker]++;
				}
			}
	}
}

void CBoundaryGeometry::SetBoundControlVolume(CConfig *config, unsigned short action) {
	unsigned short Neighbor_Node, iMarker, iNode, iNeighbor_Nodes, iDim;
	unsigned long Neighbor_Point, iVertex, iPoint, iElem;
	double **Coord;
	unsigned short nNode;
	unsigned long elem_poin;
    
	/*--- Center of gravity for face elements ---*/
	for(iMarker = 0; iMarker < nMarker; iMarker++)
		for(iElem = 0; iElem < nElem_Bound[iMarker]; iElem++) {
			nNode = bound[iMarker][iElem]->GetnNodes();
			Coord = new double* [nNode];
			/*--- Store the coordinates for all the element nodes ---*/
			for (iNode = 0; iNode < nNode; iNode++) {
				elem_poin = bound[iMarker][iElem]->GetNode(iNode);
				Coord[iNode] = new double [nDim];
				for (iDim = 0; iDim < nDim; iDim++)
					Coord[iNode][iDim] = node[elem_poin]->GetCoord(iDim);
			}
			/*--- Compute the element CG coordinates ---*/
			bound[iMarker][iElem]->SetCG(Coord);
			for (iNode=0; iNode < nNode; iNode++)
				if (Coord[iNode] != NULL) delete[] Coord[iNode];
			if (Coord != NULL) delete[] Coord;
		}
    
	double *Coord_Edge_CG = new double [nDim];
	double *Coord_Elem_CG = new double [nDim];
	double *Coord_Vertex = new double [nDim];
    
	/*--- Loop over all the markers ---*/
	for (iMarker = 0; iMarker < nMarker; iMarker++)
    /*--- Loop over all the boundary elements ---*/
		for (iElem = 0; iElem < nElem_Bound[iMarker]; iElem++)
        /*--- Loop over all the nodes of the boundary ---*/
			for(iNode = 0; iNode < bound[iMarker][iElem]->GetnNodes(); iNode++) {
				iPoint = bound[iMarker][iElem]->GetNode(iNode);
				iVertex = node[iPoint]->GetVertex(iMarker);
				/*--- Loop over the neighbor nodes, there is a face for each one ---*/
				for(iNeighbor_Nodes = 0; iNeighbor_Nodes < bound[iMarker][iElem]->GetnNeighbor_Nodes(iNode); iNeighbor_Nodes++) {
					Neighbor_Node  = bound[iMarker][iElem]->GetNeighbor_Nodes(iNode,iNeighbor_Nodes);
					Neighbor_Point = bound[iMarker][iElem]->GetNode(Neighbor_Node);
					/*--- Shared edge by the Neighbor Point and the point ---*/
					for (iDim = 0; iDim < nDim; iDim++) {
						Coord_Edge_CG[iDim] = 0.5*(node[iPoint]->GetCoord(iDim) + node[Neighbor_Point]->GetCoord(iDim));
						Coord_Elem_CG[iDim] = bound[iMarker][iElem]->GetCG(iDim);
						Coord_Vertex[iDim]  = node[iPoint]->GetCoord(iDim);
					}
                    
					vertex[iMarker][iVertex]->SetCoord(Coord_Vertex);
                    
					switch (nDim) {
                        case 2:
                            /*--- Store the 2D face (ojo hay cambio de sentido para ajustarse al sentido del contorno de nodo 0 al 1) ---*/
                            if (iNode == 0) vertex[iMarker][iVertex]->SetNodes_Coord(Coord_Elem_CG, Coord_Vertex, config);
                            if (iNode == 1) vertex[iMarker][iVertex]->SetNodes_Coord(Coord_Vertex, Coord_Elem_CG, config);
                            break;
                        case 3:
                            /*--- Store the 3D face (ojo hay cambio de sentido para ajustarse al sentido del contorno de nodo 0 al 1) ---*/
                            if (iNeighbor_Nodes == 0) vertex[iMarker][iVertex]->SetNodes_Coord(Coord_Elem_CG, Coord_Edge_CG, Coord_Vertex, config);
                            if (iNeighbor_Nodes == 1) vertex[iMarker][iVertex]->SetNodes_Coord(Coord_Edge_CG, Coord_Elem_CG, Coord_Vertex, config);
					}
				}
			}
    
	delete[] Coord_Edge_CG;
	delete[] Coord_Elem_CG;
	delete[] Coord_Vertex;
}

void CBoundaryGeometry::SetBoundSensitivity(CConfig *config) {
	unsigned short iMarker, icommas;
	unsigned long iVertex, iPoint, (*Point2Vertex)[2], nPointLocal = 0, nPointGlobal = 0;
	double Sensitivity;
	bool *PointInDomain;
	int rank = MASTER_NODE;
	int size = SINGLE_NODE;
    
	nPointLocal = nPoint;
	nPointGlobal = nPointLocal;
    
	Point2Vertex = new unsigned long[nPointGlobal][2];
	PointInDomain = new bool[nPointGlobal];
    
	for (iPoint = 0; iPoint < nPointGlobal; iPoint ++)
		PointInDomain[iPoint] = false;
    
	for (iMarker = 0; iMarker < nMarker; iMarker++)
		if (config->GetMarker_All_DV(iMarker) == YES)
			for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
                
				/*--- The sensitivity file uses the global numbering ---*/
				iPoint = vertex[iMarker][iVertex]->GetNode();
				if (vertex[iMarker][iVertex]->GetNode() < GetnPointDomain()) {
					Point2Vertex[iPoint][0] = iMarker;
					Point2Vertex[iPoint][1] = iVertex;
					PointInDomain[iPoint] = true;
					vertex[iMarker][iVertex]->SetAuxVar(0.0);
				}
			}
    
	/*--- Time-average any unsteady surface sensitivities ---*/
	unsigned long iExtIter, nExtIter;
	double delta_T, total_T;
	if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
		nExtIter = config->GetnExtIter();
		delta_T  = config->GetDelta_UnstTimeND();
		total_T  = (double)nExtIter*delta_T;
	} else if (config->GetUnsteady_Simulation() == TIME_SPECTRAL) {
        
		/*--- Compute period of oscillation & compute time interval using nTimeInstances ---*/
		double period = config->GetTimeSpectral_Period();
		nExtIter  = config->GetnTimeInstances();
		delta_T   = period/(double)nExtIter;
		total_T   = period;
        
	} else {
		nExtIter = 1;
		delta_T  = 1.0;
		total_T  = 1.0;
	}
    
	for (iExtIter = 0; iExtIter < nExtIter; iExtIter++) {
        
		/*--- Prepare to read surface sensitivity files (CSV) ---*/
		string text_line;
		ifstream Surface_file;
		char buffer[50];
		char cstr[200];
		string surfadj_filename = config->GetSurfAdjCoeff_FileName();
        
		/*--- Remove the domain number from the surface csv filename ---*/
		if (size > SINGLE_NODE) {
			if ((rank+1 >= 0) && (rank+1 < 10)) surfadj_filename.erase (surfadj_filename.end()-2, surfadj_filename.end());
			if ((rank+1 >= 10) && (rank+1 < 100)) surfadj_filename.erase (surfadj_filename.end()-3, surfadj_filename.end());
			if ((rank+1 >= 100) && (rank+1 < 1000)) surfadj_filename.erase (surfadj_filename.end()-4, surfadj_filename.end());
			if ((rank+1 >= 1000) && (rank+1 < 10000)) surfadj_filename.erase (surfadj_filename.end()-5, surfadj_filename.end());
		}
		strcpy (cstr, surfadj_filename.c_str());
        
		/*--- Write file name with extension if unsteady or steady ---*/
		if ((config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) ||
            (config->GetUnsteady_Simulation() == TIME_SPECTRAL)) {
			if ((int(iExtIter) >= 0)    && (int(iExtIter) < 10))    sprintf (buffer, "_0000%d.csv", int(iExtIter));
			if ((int(iExtIter) >= 10)   && (int(iExtIter) < 100))   sprintf (buffer, "_000%d.csv",  int(iExtIter));
			if ((int(iExtIter) >= 100)  && (int(iExtIter) < 1000))  sprintf (buffer, "_00%d.csv",   int(iExtIter));
			if ((int(iExtIter) >= 1000) && (int(iExtIter) < 10000)) sprintf (buffer, "_0%d.csv",    int(iExtIter));
			if  (int(iExtIter) >= 10000) sprintf (buffer, "_%d.csv", int(iExtIter));
		}
		else
			sprintf (buffer, ".csv");
        
		strcat (cstr, buffer);
        
		/*--- Read the sensitivity file ---*/
		string::size_type position;
        
		Surface_file.open(cstr, ios::in);
		getline(Surface_file,text_line);
        
		while (getline(Surface_file,text_line)) {
			for (icommas = 0; icommas < 50; icommas++) {
				position = text_line.find( ",", 0 );
				if(position!=string::npos) text_line.erase (position,1);
			}
			stringstream  point_line(text_line);
			point_line >> iPoint >> Sensitivity;
            
			if (PointInDomain[iPoint]) {
                
				/*--- Find the vertex for the Point and Marker ---*/
				iMarker = Point2Vertex[iPoint][0];
				iVertex = Point2Vertex[iPoint][1];
                
				/*--- Increment the auxiliary variable with the contribution of
                 this unsteady timestep. For steady problems, this reduces to
                 a single sensitivity value multiplied by 1.0. ---*/
				vertex[iMarker][iVertex]->AddAuxVar(Sensitivity*(delta_T/total_T));
			}
            
		}
		Surface_file.close();
	}
    
	delete[] Point2Vertex;
}

void CBoundaryGeometry::ComputeAirfoil_Section(double *Plane_P0, double *Plane_Normal, unsigned short iSection, CConfig *config,
                                               vector<double> &Xcoord_Airfoil, vector<double> &Ycoord_Airfoil,
                                               vector<double> &Zcoord_Airfoil, bool original_surface) {
  
	unsigned short iMarker, iNode, jNode, iDim, intersect;
	long MinDist_Point, MinDistAngle_Point;
	unsigned long iPoint, jPoint, iElem, Trailing_Point, Airfoil_Point, iVertex, jVertex, n;
	double Segment_P0[3] = {0.0, 0.0, 0.0}, Segment_P1[3] = {0.0, 0.0, 0.0}, Intersection[3] = {0.0, 0.0, 0.0}, Trailing_Coord, MinDist_Value, MinDistAngle_Value, Dist_Value,
  Airfoil_Tangent[3] = {0.0, 0.0, 0.0}, Segment[3] = {0.0, 0.0, 0.0}, Length, Angle_Value, Normal[3], Tangent[3], BiNormal[3], auxXCoord,
  auxYCoord, auxZCoord, zp1, zpn, Camber_Line, MaxAngle = 25, *VarCoord = NULL, CosValue;
	vector<double> Xcoord, Ycoord, Zcoord, Z2coord, Xcoord_Normal, Ycoord_Normal, Zcoord_Normal, Xcoord_Camber, Ycoord_Camber, Zcoord_Camber;
	vector<unsigned long> Duplicate;
	vector<unsigned long>::iterator it;
	int rank = MASTER_NODE;
	double **Coord_Variation;
  
	Xcoord_Airfoil.clear();
	Ycoord_Airfoil.clear();
	Zcoord_Airfoil.clear();
  
	/*--- Set the right plane in 2D (note the change in Y-Z plane) ---*/
	if (nDim == 2) {
		iSection = 0;
		Plane_P0[0] = 0.0;      Plane_P0[1] = 0.0;      Plane_P0[2] = 0.0;
		Plane_Normal[0] = 0.0;  Plane_Normal[1] = 1.0;  Plane_Normal[2] = 0.0;
	}
  
	/*--- the grid variation is stored using a vertices information,
   we should go from vertex to points ---*/
	if (original_surface == false) {
    
		Coord_Variation = new double *[nPoint];
		for (iPoint = 0; iPoint < nPoint; iPoint++)
			Coord_Variation[iPoint] = new double [nDim];
    
		for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
			if (config->GetMarker_All_Monitoring(iMarker) == YES) {
				for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
					VarCoord = vertex[iMarker][iVertex]->GetVarCoord();
					iPoint = vertex[iMarker][iVertex]->GetNode();
					for (iDim = 0; iDim < nDim; iDim++)
						Coord_Variation[iPoint][iDim] = VarCoord[iDim];
				}
			}
		}
    
	}
  
	for (iMarker = 0; iMarker < nMarker; iMarker++) {
		if (config->GetMarker_All_Monitoring(iMarker) == YES) {
			for (iElem = 0; iElem < nElem_Bound[iMarker]; iElem++) {
				for(iNode = 0; iNode < bound[iMarker][iElem]->GetnNodes(); iNode++) {
					iPoint = bound[iMarker][iElem]->GetNode(iNode);
					for(jNode = 0; jNode < bound[iMarker][iElem]->GetnNodes(); jNode++) {
						jPoint = bound[iMarker][iElem]->GetNode(jNode);
            
						if (jPoint > iPoint) {
              
							Segment_P0[0] = 0.0;  Segment_P0[1] = 0.0;  Segment_P0[2] = 0.0;
							Segment_P1[0] = 0.0;  Segment_P1[1] = 0.0;  Segment_P1[2] = 0.0;
              
							for (iDim = 0; iDim < nDim; iDim++) {
								if (original_surface == true) {
									Segment_P0[iDim] = node[iPoint]->GetCoord(iDim);
									Segment_P1[iDim] = node[jPoint]->GetCoord(iDim);
								}
								else {
									Segment_P0[iDim] = node[iPoint]->GetCoord(iDim) + Coord_Variation[iPoint][iDim];
									Segment_P1[iDim] = node[jPoint]->GetCoord(iDim) + Coord_Variation[jPoint][iDim];
								}
							}
              
							/*--- In 2D add the points directly (note the change between Y and Z coordinate) ---*/
							if (nDim == 2) {
								Xcoord.push_back(Segment_P0[0]);  Xcoord.push_back(Segment_P1[0]);
								Ycoord.push_back(Segment_P0[2]);  Ycoord.push_back(Segment_P1[2]);
								Zcoord.push_back(Segment_P0[1]);  Zcoord.push_back(Segment_P1[1]);
							}
							/*--- In 3D compute the intersection ---*/
							else if (nDim == 3) {
								intersect = ComputeSegmentPlane_Intersection(Segment_P0, Segment_P1, Plane_P0, Plane_Normal, Intersection);
								if (intersect == 1) {
									Xcoord.push_back(Intersection[0]);
									Ycoord.push_back(Intersection[1]);
									Zcoord.push_back(Intersection[2]);
								}
							}
              
						}
            
					}
				}
			}
		}
	}
  
	if (original_surface == false) {
    
		for (iPoint = 0; iPoint < nPoint; iPoint++)
			delete [] Coord_Variation[iPoint];
		delete [] Coord_Variation;
    
	}
  
  
  if ((rank == MASTER_NODE) && (Xcoord.size() != 0)) {
    
    /*--- Create a list with the duplicated points ---*/
    for (iVertex = 0; iVertex < Xcoord.size()-1; iVertex++) {
      for (jVertex = iVertex+1; jVertex < Xcoord.size(); jVertex++) {
        Segment[0] = Xcoord[jVertex] - Xcoord[iVertex];
        Segment[1] = Ycoord[jVertex] - Ycoord[iVertex];
        Segment[2] = Zcoord[jVertex] - Zcoord[iVertex];
        Dist_Value = sqrt(pow(Segment[0], 2.0) + pow(Segment[1], 2.0) + pow(Segment[2], 2.0));
        if (Dist_Value < 1E-6) {
          Duplicate.push_back (jVertex);
        }
      }
    }
    
    sort(Duplicate.begin(), Duplicate.end());
    it = unique(Duplicate.begin(), Duplicate.end());
    Duplicate.resize(it - Duplicate.begin());
    
    /*--- Remove duplicated points (starting from the back) ---*/
    for (iVertex = Duplicate.size(); iVertex > 0; iVertex--) {
      Xcoord.erase (Xcoord.begin() + Duplicate[iVertex-1]);
      Ycoord.erase (Ycoord.begin() + Duplicate[iVertex-1]);
      Zcoord.erase (Zcoord.begin() + Duplicate[iVertex-1]);
    }
    
    /*--- Find the trailing edge ---*/
    Trailing_Point = 0; Trailing_Coord = Xcoord[0];
    for (iVertex = 1; iVertex < Xcoord.size(); iVertex++) {
      if (Xcoord[iVertex] > Trailing_Coord) {
        Trailing_Point = iVertex; Trailing_Coord = Xcoord[iVertex];
      }
    }
    
    /*--- Add the trailing edge to the list, and remove from the original list ---*/
    Xcoord_Airfoil.push_back(Xcoord[Trailing_Point]); Ycoord_Airfoil.push_back(Ycoord[Trailing_Point]); Zcoord_Airfoil.push_back(Zcoord[Trailing_Point]);
    Xcoord.erase (Xcoord.begin() + Trailing_Point); Ycoord.erase (Ycoord.begin() + Trailing_Point); Zcoord.erase (Zcoord.begin() + Trailing_Point);
    
    /*--- Find the next point using the right hand side rule ---*/
    MinDist_Value = 1E6;
    for (iVertex = 0; iVertex < Xcoord.size(); iVertex++) {
      Segment[0] = Xcoord[iVertex] - Xcoord_Airfoil[0];
      Segment[1] = Ycoord[iVertex] - Ycoord_Airfoil[0];
      Segment[2] = Zcoord[iVertex] - Zcoord_Airfoil[0];
      Dist_Value = sqrt(pow(Segment[0], 2.0) + pow(Segment[1], 2.0) + pow(Segment[2], 2.0));
      Segment[0] /= Dist_Value; Segment[1] /= Dist_Value; Segment[2] /= Dist_Value;
      
      if ((Dist_Value < MinDist_Value) && (Segment[2] > 0.0)) { MinDist_Point = iVertex; MinDist_Value = Dist_Value; }
    }
    
		Xcoord_Airfoil.push_back(Xcoord[MinDist_Point]);  Ycoord_Airfoil.push_back(Ycoord[MinDist_Point]);  Zcoord_Airfoil.push_back(Zcoord[MinDist_Point]);
		Xcoord.erase (Xcoord.begin() + MinDist_Point);    Ycoord.erase (Ycoord.begin() + MinDist_Point);    Zcoord.erase (Zcoord.begin() + MinDist_Point);
    
		/*--- Algorithm for the rest of the points ---*/
		do {
      
			/*--- Last added point in the list ---*/
			Airfoil_Point = Xcoord_Airfoil.size() - 1;
      
			/*--- Compute the slope of the curve ---*/
			Airfoil_Tangent[0] = Xcoord_Airfoil[Airfoil_Point] - Xcoord_Airfoil[Airfoil_Point-1];
			Airfoil_Tangent[1] = Ycoord_Airfoil[Airfoil_Point] - Ycoord_Airfoil[Airfoil_Point-1];
			Airfoil_Tangent[2] = Zcoord_Airfoil[Airfoil_Point] - Zcoord_Airfoil[Airfoil_Point-1];
			Length = sqrt(pow(Airfoil_Tangent[0], 2.0) + pow(Airfoil_Tangent[1], 2.0) + pow(Airfoil_Tangent[2], 2.0));
			Airfoil_Tangent[0] /= Length; Airfoil_Tangent[1] /= Length; Airfoil_Tangent[2] /= Length;
      
			/*--- Find the closest point with the right slope ---*/
			MinDist_Value = 1E6; MinDistAngle_Value = 180;
			MinDist_Point = -1; MinDistAngle_Point = -1;
			for (iVertex = 0; iVertex < Xcoord.size(); iVertex++) {
        
				Segment[0] = Xcoord[iVertex] - Xcoord_Airfoil[Airfoil_Point];
				Segment[1] = Ycoord[iVertex] - Ycoord_Airfoil[Airfoil_Point];
				Segment[2] = Zcoord[iVertex] - Zcoord_Airfoil[Airfoil_Point];
        
				/*--- Compute the distance to each point ---*/
				Dist_Value = sqrt(pow(Segment[0], 2.0) + pow(Segment[1], 2.0) + pow(Segment[2], 2.0));
        
				/*--- Compute the angle of the point ---*/
				Segment[0] /= Dist_Value; Segment[1] /= Dist_Value; Segment[2] /= Dist_Value;
        
				/*--- Clip the value of the cosine, this is important due to the round errors ---*/
				CosValue = Airfoil_Tangent[0]*Segment[0] + Airfoil_Tangent[1]*Segment[1] + Airfoil_Tangent[2]*Segment[2];
				if (CosValue >= 1.0) CosValue = 1.0;
				if (CosValue <= -1.0) CosValue = -1.0;
        
				Angle_Value = acos(CosValue) * 180 / PI_NUMBER;
        
				if (Dist_Value < MinDist_Value) { MinDist_Point = iVertex; MinDist_Value = Dist_Value; }
				if ((Dist_Value < MinDistAngle_Value) && (Angle_Value < MaxAngle)) {MinDistAngle_Point = iVertex; MinDistAngle_Value = Dist_Value;}
        
			}
      
			if ( MinDistAngle_Point != -1) MinDist_Point = MinDistAngle_Point;
      
			/*--- Add and remove the min distance to the list ---*/
			Xcoord_Airfoil.push_back(Xcoord[MinDist_Point]);  Ycoord_Airfoil.push_back(Ycoord[MinDist_Point]);  Zcoord_Airfoil.push_back(Zcoord[MinDist_Point]);
			Xcoord.erase(Xcoord.begin() + MinDist_Point);    Ycoord.erase(Ycoord.begin() + MinDist_Point);    Zcoord.erase(Zcoord.begin() + MinDist_Point);
      
		} while (Xcoord.size() != 0);
    
		/*--- Clean the vector before using them again for storing the upper or the lower side ---*/
		Xcoord.clear(); Ycoord.clear(); Zcoord.clear();
    
		/*--- Identify upper and lower side using the value of the normal vector --*/
		for (iVertex = 1; iVertex < Xcoord_Airfoil.size(); iVertex++) {
			Tangent[0] = Xcoord_Airfoil[iVertex] - Xcoord_Airfoil[iVertex-1];
			Tangent[1] = Ycoord_Airfoil[iVertex] - Ycoord_Airfoil[iVertex-1];
			Tangent[2] = Zcoord_Airfoil[iVertex] - Zcoord_Airfoil[iVertex-1];
			Length = sqrt(pow(Tangent[0], 2.0) + pow(Tangent[1], 2.0) + pow(Tangent[2], 2.0));
			Tangent[0] /= Length; Tangent[1] /= Length; Tangent[2] /= Length;
      
			BiNormal[0] = Plane_Normal[0];
			BiNormal[1] = Plane_Normal[1];
			BiNormal[2] = Plane_Normal[2];
			Length = sqrt(pow(BiNormal[0], 2.0) + pow(BiNormal[1], 2.0) + pow(BiNormal[2], 2.0));
			BiNormal[0] /= Length; BiNormal[1] /= Length; BiNormal[2] /= Length;
      
			Normal[0] = Tangent[1]*BiNormal[2] - Tangent[2]*BiNormal[1];
			Normal[1] = Tangent[2]*BiNormal[0] - Tangent[0]*BiNormal[2];
			Normal[2] = Tangent[0]*BiNormal[1] - Tangent[1]*BiNormal[0];
      
			Xcoord_Normal.push_back(Normal[0]); Ycoord_Normal.push_back(Normal[1]); Zcoord_Normal.push_back(Normal[2]);
      
			if (Normal[2] >= 0.0) {
				Xcoord.push_back(Xcoord_Airfoil[iVertex]);
				Ycoord.push_back(Ycoord_Airfoil[iVertex]);
				Zcoord.push_back(Zcoord_Airfoil[iVertex]);
			}
      
		}
    
		n = Xcoord.size();
    
		/*--- Order the arrays using the X component ---*/
		for (iVertex = 0; iVertex < Xcoord.size(); iVertex++) {
			for (jVertex = 0; jVertex < Xcoord.size() - 1 - iVertex; jVertex++) {
				if (Xcoord[jVertex] > Xcoord[jVertex+1]) {
					auxXCoord = Xcoord[jVertex]; Xcoord[jVertex] = Xcoord[jVertex+1]; Xcoord[jVertex+1] = auxXCoord;
					auxYCoord = Ycoord[jVertex]; Ycoord[jVertex] = Ycoord[jVertex+1]; Ycoord[jVertex+1] = auxYCoord;
					auxZCoord = Zcoord[jVertex]; Zcoord[jVertex] = Zcoord[jVertex+1]; Zcoord[jVertex+1] = auxZCoord;
				}
			}
		}
    
		zp1=(Zcoord[1]-Zcoord[0]) / (Xcoord[1]-Xcoord[0]);
		zpn=(Zcoord[n-1]-Zcoord[n-2]) / (Xcoord[n-1]-Xcoord[n-2]);
		Z2coord.resize(n+1);
		SetSpline(Xcoord, Zcoord, n, zp1, zpn, Z2coord);
    
		/*--- Compute the camber--*/
		for (iVertex = 0; iVertex < Xcoord_Airfoil.size(); iVertex++) {
			if (Zcoord_Normal[iVertex] < 0.0) {
				Camber_Line = 0.5 * (Zcoord_Airfoil[iVertex] + GetSpline(Xcoord, Zcoord, Z2coord, n, Xcoord_Airfoil[iVertex]));
				Xcoord_Camber.push_back(Xcoord_Airfoil[iVertex]);
				Ycoord_Camber.push_back(Ycoord_Airfoil[iVertex]);
				Zcoord_Camber.push_back(Camber_Line);
			}
		}
    
		/*--- Write the output file (tecplot format) ---*/
		if (original_surface == true) {
			ofstream Tecplot_File;
			if (iSection == 0) Tecplot_File.open("Airfoil_Sections.plt", ios::out);
			else Tecplot_File.open("Airfoil_Sections.plt", ios::app);
      
			if (iSection == 0) {
				Tecplot_File << "TITLE = \"Wing airfoil sections\"" << endl;
				Tecplot_File << "VARIABLES = \"X\",\"Y\",\"Z\"" << endl;
			}
      
			Tecplot_File << "ZONE T=\"SECTION_"<< (iSection+1) << "\", NODES= "<< Xcoord_Airfoil.size() << ", ELEMENTS= " << Xcoord_Airfoil.size()-1 << ", DATAPACKING= POINT, ZONETYPE= FELINESEG" << endl;
      
			/*--- Coordinates ---*/
			for (iVertex = 0; iVertex < Xcoord_Airfoil.size(); iVertex++) {
				Tecplot_File << Xcoord_Airfoil[iVertex] <<" "<< Ycoord_Airfoil[iVertex] <<" "<< Zcoord_Airfoil[iVertex] << endl;
			}
			/*--- Conectivity ---*/
			for (iVertex = 1; iVertex < Xcoord_Airfoil.size(); iVertex++) {
				Tecplot_File << iVertex << "\t" << iVertex+1 << "\n";
			}
      
			Tecplot_File.close();
		}
    
	}
  
}

double CBoundaryGeometry::Compute_MaxThickness(double *Plane_P0, double *Plane_Normal, unsigned short iSection, vector<double> &Xcoord_Airfoil, vector<double> &Ycoord_Airfoil, vector<double> &Zcoord_Airfoil, bool original_surface) {
	unsigned long iVertex, jVertex, n, Trailing_Point, LeadingPoint;
	double Normal[3], Tangent[3], BiNormal[3], auxXCoord, auxYCoord, auxZCoord, zp1, zpn, MaxThickness_Value = 0, MaxThickness_Location, Thickness, Length, Xcoord_Trailing, Ycoord_Trailing, Zcoord_Trailing, ValCos, ValSin, XValue, ZValue, MaxDistance, Distance, AoA;
	vector<double> Xcoord, Ycoord, Zcoord, Z2coord, Xcoord_Normal, Ycoord_Normal, Zcoord_Normal, Xcoord_Airfoil_, Ycoord_Airfoil_, Zcoord_Airfoil_;
    
    /*--- Find the leading and trailing edges and compute the angle of attack ---*/
    MaxDistance = 0.0; Trailing_Point = 0;
    for (iVertex = 1; iVertex < Xcoord_Airfoil.size(); iVertex++) {
        
        Distance = sqrt(pow(Xcoord_Airfoil[iVertex] - Xcoord_Airfoil[Trailing_Point], 2.0) +
                        pow(Ycoord_Airfoil[iVertex] - Ycoord_Airfoil[Trailing_Point], 2.0) +
                        pow(Zcoord_Airfoil[iVertex] - Zcoord_Airfoil[Trailing_Point], 2.0));
        
        if (MaxDistance < Distance) { MaxDistance = Distance; LeadingPoint = iVertex; }
    }
    
    AoA = atan((Zcoord_Airfoil[LeadingPoint] - Zcoord_Airfoil[Trailing_Point]) / (Xcoord_Airfoil[Trailing_Point] - Xcoord_Airfoil[LeadingPoint]))*180/PI_NUMBER;
    
    /*--- Translate to the origin ---*/
    Xcoord_Trailing = Xcoord_Airfoil[0];
    Ycoord_Trailing = Ycoord_Airfoil[0];
    Zcoord_Trailing = Zcoord_Airfoil[0];
    
    for (iVertex = 0; iVertex < Xcoord_Airfoil.size(); iVertex++) {
        Xcoord_Airfoil_.push_back(Xcoord_Airfoil[iVertex] - Xcoord_Trailing);
        Ycoord_Airfoil_.push_back(Ycoord_Airfoil[iVertex] - Ycoord_Trailing);
        Zcoord_Airfoil_.push_back(Zcoord_Airfoil[iVertex] - Zcoord_Trailing);
    }
    
    /*--- Rotate the airfoil ---*/
    ValCos = cos(AoA*PI_NUMBER/180.0);
    ValSin = sin(AoA*PI_NUMBER/180.0);
    
    for (iVertex = 0; iVertex < Xcoord_Airfoil.size(); iVertex++) {
        XValue = Xcoord_Airfoil_[iVertex];
        ZValue = Zcoord_Airfoil_[iVertex];
        
        Xcoord_Airfoil_[iVertex] = XValue*ValCos - ZValue*ValSin;
        Zcoord_Airfoil_[iVertex] = ZValue*ValCos + XValue*ValSin;
        
    }
    
    /*--- Identify upper and lower side, and store the value of the normal --*/
    for (iVertex = 1; iVertex < Xcoord_Airfoil_.size(); iVertex++) {
        Tangent[0] = Xcoord_Airfoil_[iVertex] - Xcoord_Airfoil_[iVertex-1];
        Tangent[1] = Ycoord_Airfoil_[iVertex] - Ycoord_Airfoil_[iVertex-1];
        Tangent[2] = Zcoord_Airfoil_[iVertex] - Zcoord_Airfoil_[iVertex-1];
        Length = sqrt(pow(Tangent[0], 2.0) + pow(Tangent[1], 2.0) + pow(Tangent[2], 2.0));
        Tangent[0] /= Length; Tangent[1] /= Length; Tangent[2] /= Length;
        
        BiNormal[0] = Plane_Normal[0];
        BiNormal[1] = Plane_Normal[1];
        BiNormal[2] = Plane_Normal[2];
        Length = sqrt(pow(BiNormal[0], 2.0) + pow(BiNormal[1], 2.0) + pow(BiNormal[2], 2.0));
        BiNormal[0] /= Length; BiNormal[1] /= Length; BiNormal[2] /= Length;
        
        Normal[0] = Tangent[1]*BiNormal[2] - Tangent[2]*BiNormal[1];
        Normal[1] = Tangent[2]*BiNormal[0] - Tangent[0]*BiNormal[2];
        Normal[2] = Tangent[0]*BiNormal[1] - Tangent[1]*BiNormal[0];
        
        Xcoord_Normal.push_back(Normal[0]); Ycoord_Normal.push_back(Normal[1]); Zcoord_Normal.push_back(Normal[2]);
        
        if (Normal[2] >= 0.0) {
            Xcoord.push_back(Xcoord_Airfoil_[iVertex]);
            Ycoord.push_back(Ycoord_Airfoil_[iVertex]);
            Zcoord.push_back(Zcoord_Airfoil_[iVertex]);
        }
        
    }
    
    /*--- Order the arrays using the X component ---*/
    for (iVertex = 0; iVertex < Xcoord.size(); iVertex++) {
        for (jVertex = 0; jVertex < Xcoord.size() - 1 - iVertex; jVertex++) {
            if (Xcoord[jVertex] > Xcoord[jVertex+1]) {
                auxXCoord = Xcoord[jVertex]; Xcoord[jVertex] = Xcoord[jVertex+1]; Xcoord[jVertex+1] = auxXCoord;
                auxYCoord = Ycoord[jVertex]; Ycoord[jVertex] = Ycoord[jVertex+1]; Ycoord[jVertex+1] = auxYCoord;
                auxZCoord = Zcoord[jVertex]; Zcoord[jVertex] = Zcoord[jVertex+1]; Zcoord[jVertex+1] = auxZCoord;
            }
        }
    }
    
    n = Xcoord.size();
    zp1 = (Zcoord[1]-Zcoord[0])/(Xcoord[1]-Xcoord[0]);
    zpn = (Zcoord[n-1]-Zcoord[n-2])/(Xcoord[n-1]-Xcoord[n-2]);
    Z2coord.resize(n+1);
    SetSpline(Xcoord, Zcoord, n, zp1, zpn, Z2coord);
    
    /*--- Compute the max thickness --*/
    MaxThickness_Value = 0.0; MaxThickness_Location = 0.0;
    for (iVertex = 0; iVertex < Xcoord_Airfoil_.size(); iVertex++) {
        if (Zcoord_Normal[iVertex] < 0.0) {
            Thickness = fabs(Zcoord_Airfoil_[iVertex]) + fabs(GetSpline(Xcoord, Zcoord, Z2coord, n, Xcoord_Airfoil_[iVertex]));
            if (Thickness > MaxThickness_Value) { MaxThickness_Value = Thickness; MaxThickness_Location = Xcoord_Airfoil_[iVertex]; }
        }
    }
    
    return MaxThickness_Value;
    
}

double CBoundaryGeometry::Compute_AoA(double *Plane_P0, double *Plane_Normal, unsigned short iSection, vector<double> &Xcoord_Airfoil, vector<double> &Ycoord_Airfoil, vector<double> &Zcoord_Airfoil, bool original_surface) {
	unsigned long iVertex, Trailing_Point, LeadingPoint;
	double MaxDistance, Distance, AoA = 0.0;
    
    /*--- Find the leading and trailing edges and compute the angle of attack ---*/
    MaxDistance = 0.0; Trailing_Point = 0;
    for (iVertex = 1; iVertex < Xcoord_Airfoil.size(); iVertex++) {
        
        Distance = sqrt(pow(Xcoord_Airfoil[iVertex] - Xcoord_Airfoil[Trailing_Point], 2.0) +
                        pow(Ycoord_Airfoil[iVertex] - Ycoord_Airfoil[Trailing_Point], 2.0) +
                        pow(Zcoord_Airfoil[iVertex] - Zcoord_Airfoil[Trailing_Point], 2.0));
        
        if (MaxDistance < Distance) { MaxDistance = Distance; LeadingPoint = iVertex; }
    }
    
    AoA = atan((Zcoord_Airfoil[LeadingPoint] - Zcoord_Airfoil[Trailing_Point]) / (Xcoord_Airfoil[Trailing_Point] - Xcoord_Airfoil[LeadingPoint]))*180/PI_NUMBER;
    
    return AoA;
    
}

double CBoundaryGeometry::Compute_Chord(double *Plane_P0, double *Plane_Normal, unsigned short iSection, vector<double> &Xcoord_Airfoil, vector<double> &Ycoord_Airfoil, vector<double> &Zcoord_Airfoil, bool original_surface) {
	unsigned long iVertex, Trailing_Point, LeadingPoint;
	double MaxDistance, Distance, Chord = 0.0;
    
    /*--- Find the leading and trailing edges and compute the angle of attack ---*/
    MaxDistance = 0.0; Trailing_Point = 0;
    for (iVertex = 1; iVertex < Xcoord_Airfoil.size(); iVertex++) {
        
        Distance = sqrt(pow(Xcoord_Airfoil[iVertex] - Xcoord_Airfoil[Trailing_Point], 2.0) +
                        pow(Ycoord_Airfoil[iVertex] - Ycoord_Airfoil[Trailing_Point], 2.0) +
                        pow(Zcoord_Airfoil[iVertex] - Zcoord_Airfoil[Trailing_Point], 2.0));
        
        if (MaxDistance < Distance) { MaxDistance = Distance; LeadingPoint = iVertex; }
    }
    
    Chord = MaxDistance;
    
    return Chord;
    
}

double CBoundaryGeometry::Compute_Thickness(double *Plane_P0, double *Plane_Normal, unsigned short iSection, double Location, vector<double> &Xcoord_Airfoil, vector<double> &Ycoord_Airfoil, vector<double> &Zcoord_Airfoil, bool original_surface) {
	unsigned long iVertex, jVertex, n_Upper, n_Lower, Trailing_Point, LeadingPoint;
	double Thickness_Location, Normal[3], Tangent[3], BiNormal[3], auxXCoord, auxYCoord, auxZCoord, Thickness_Value = 0.0, Length, Xcoord_Trailing, Ycoord_Trailing, Zcoord_Trailing, ValCos, ValSin, XValue, ZValue, zp1, zpn, Chord, MaxDistance, Distance, AoA;
	vector<double> Xcoord_Upper, Ycoord_Upper, Zcoord_Upper, Z2coord_Upper, Xcoord_Lower, Ycoord_Lower, Zcoord_Lower, Z2coord_Lower, Z2coord, Xcoord_Normal, Ycoord_Normal, Zcoord_Normal, Xcoord_Airfoil_, Ycoord_Airfoil_, Zcoord_Airfoil_;
    
    /*--- Find the leading and trailing edges and compute the angle of attack ---*/
    MaxDistance = 0.0; Trailing_Point = 0;
    for (iVertex = 1; iVertex < Xcoord_Airfoil.size(); iVertex++) {
        
        Distance = sqrt(pow(Xcoord_Airfoil[iVertex] - Xcoord_Airfoil[Trailing_Point], 2.0) +
                        pow(Ycoord_Airfoil[iVertex] - Ycoord_Airfoil[Trailing_Point], 2.0) +
                        pow(Zcoord_Airfoil[iVertex] - Zcoord_Airfoil[Trailing_Point], 2.0));
        
        if (MaxDistance < Distance) { MaxDistance = Distance; LeadingPoint = iVertex; }
    }
    
    AoA = atan((Zcoord_Airfoil[LeadingPoint] - Zcoord_Airfoil[Trailing_Point]) / (Xcoord_Airfoil[Trailing_Point] - Xcoord_Airfoil[LeadingPoint]))*180/PI_NUMBER;
    Chord = MaxDistance;
    
    /*--- Translate to the origin ---*/
    Xcoord_Trailing = Xcoord_Airfoil[0];
    Ycoord_Trailing = Ycoord_Airfoil[0];
    Zcoord_Trailing = Zcoord_Airfoil[0];
    
    for (iVertex = 0; iVertex < Xcoord_Airfoil.size(); iVertex++) {
        Xcoord_Airfoil_.push_back(Xcoord_Airfoil[iVertex] - Xcoord_Trailing);
        Ycoord_Airfoil_.push_back(Ycoord_Airfoil[iVertex] - Ycoord_Trailing);
        Zcoord_Airfoil_.push_back(Zcoord_Airfoil[iVertex] - Zcoord_Trailing);
    }
    
    /*--- Rotate the airfoil ---*/
    ValCos = cos(AoA*PI_NUMBER/180.0);
    ValSin = sin(AoA*PI_NUMBER/180.0);
    
    for (iVertex = 0; iVertex < Xcoord_Airfoil.size(); iVertex++) {
        XValue = Xcoord_Airfoil_[iVertex];
        ZValue = Zcoord_Airfoil_[iVertex];
        
        Xcoord_Airfoil_[iVertex] = XValue*ValCos - ZValue*ValSin;
        Zcoord_Airfoil_[iVertex] = ZValue*ValCos + XValue*ValSin;
    }
    
    /*--- Identify upper and lower side, and store the value of the normal --*/
    for (iVertex = 1; iVertex < Xcoord_Airfoil_.size(); iVertex++) {
        Tangent[0] = Xcoord_Airfoil_[iVertex] - Xcoord_Airfoil_[iVertex-1];
        Tangent[1] = Ycoord_Airfoil_[iVertex] - Ycoord_Airfoil_[iVertex-1];
        Tangent[2] = Zcoord_Airfoil_[iVertex] - Zcoord_Airfoil_[iVertex-1];
        Length = sqrt(pow(Tangent[0], 2.0) + pow(Tangent[1], 2.0) + pow(Tangent[2], 2.0));
        Tangent[0] /= Length; Tangent[1] /= Length; Tangent[2] /= Length;
        
        BiNormal[0] = Plane_Normal[0];
        BiNormal[1] = Plane_Normal[1];
        BiNormal[2] = Plane_Normal[2];
        Length = sqrt(pow(BiNormal[0], 2.0) + pow(BiNormal[1], 2.0) + pow(BiNormal[2], 2.0));
        BiNormal[0] /= Length; BiNormal[1] /= Length; BiNormal[2] /= Length;
        
        Normal[0] = Tangent[1]*BiNormal[2] - Tangent[2]*BiNormal[1];
        Normal[1] = Tangent[2]*BiNormal[0] - Tangent[0]*BiNormal[2];
        Normal[2] = Tangent[0]*BiNormal[1] - Tangent[1]*BiNormal[0];
        
        Xcoord_Normal.push_back(Normal[0]); Ycoord_Normal.push_back(Normal[1]); Zcoord_Normal.push_back(Normal[2]);
        
        if (Normal[2] >= 0.0) {
            Xcoord_Upper.push_back(Xcoord_Airfoil_[iVertex]);
            Ycoord_Upper.push_back(Ycoord_Airfoil_[iVertex]);
            Zcoord_Upper.push_back(Zcoord_Airfoil_[iVertex]);
        }
        else {
            Xcoord_Lower.push_back(Xcoord_Airfoil_[iVertex]);
            Ycoord_Lower.push_back(Ycoord_Airfoil_[iVertex]);
            Zcoord_Lower.push_back(Zcoord_Airfoil_[iVertex]);
        }
        
    }
    
    /*--- Order the arrays using the X component ---*/
    for (iVertex = 0; iVertex < Xcoord_Upper.size(); iVertex++) {
        for (jVertex = 0; jVertex < Xcoord_Upper.size() - 1 - iVertex; jVertex++) {
            if (Xcoord_Upper[jVertex] > Xcoord_Upper[jVertex+1]) {
                auxXCoord = Xcoord_Upper[jVertex]; Xcoord_Upper[jVertex] = Xcoord_Upper[jVertex+1]; Xcoord_Upper[jVertex+1] = auxXCoord;
                auxYCoord = Ycoord_Upper[jVertex]; Ycoord_Upper[jVertex] = Ycoord_Upper[jVertex+1]; Ycoord_Upper[jVertex+1] = auxYCoord;
                auxZCoord = Zcoord_Upper[jVertex]; Zcoord_Upper[jVertex] = Zcoord_Upper[jVertex+1]; Zcoord_Upper[jVertex+1] = auxZCoord;
            }
        }
    }
    
    /*--- Order the arrays using the X component ---*/
    for (iVertex = 0; iVertex < Xcoord_Lower.size(); iVertex++) {
        for (jVertex = 0; jVertex < Xcoord_Lower.size() - 1 - iVertex; jVertex++) {
            if (Xcoord_Lower[jVertex] > Xcoord_Lower[jVertex+1]) {
                auxXCoord = Xcoord_Lower[jVertex]; Xcoord_Lower[jVertex] = Xcoord_Lower[jVertex+1]; Xcoord_Lower[jVertex+1] = auxXCoord;
                auxYCoord = Ycoord_Lower[jVertex]; Ycoord_Lower[jVertex] = Ycoord_Lower[jVertex+1]; Ycoord_Lower[jVertex+1] = auxYCoord;
                auxZCoord = Zcoord_Lower[jVertex]; Zcoord_Lower[jVertex] = Zcoord_Lower[jVertex+1]; Zcoord_Lower[jVertex+1] = auxZCoord;
            }
        }
    }
    
    n_Upper = Xcoord_Upper.size();
    zp1 = (Xcoord_Upper[1]-Xcoord_Upper[0])/(Xcoord_Upper[1]-Xcoord_Upper[0]);
    zpn = (Xcoord_Upper[n_Upper-1]-Xcoord_Upper[n_Upper-2])/(Xcoord_Upper[n_Upper-1]-Xcoord_Upper[n_Upper-2]);
    Z2coord_Upper.resize(n_Upper+1);
    SetSpline(Xcoord_Upper, Zcoord_Upper, n_Upper, zp1, zpn, Z2coord_Upper);
    
    n_Lower = Xcoord_Lower.size();
    zp1 = (Xcoord_Lower[1]-Xcoord_Lower[0])/(Xcoord_Lower[1]-Xcoord_Lower[0]);
    zpn = (Xcoord_Lower[n_Lower-1]-Xcoord_Lower[n_Lower-2])/(Xcoord_Lower[n_Lower-1]-Xcoord_Lower[n_Lower-2]);
    Z2coord_Lower.resize(n_Lower+1);
    SetSpline(Xcoord_Lower, Zcoord_Lower, n_Lower, zp1, zpn, Z2coord_Lower);
    
    /*--- Compute thickness location ---*/
    Thickness_Location = - Chord*(1.0-Location);
    
    Thickness_Value = fabs(GetSpline(Xcoord_Upper, Zcoord_Upper, Z2coord_Upper, n_Upper, Thickness_Location)) + fabs(GetSpline(Xcoord_Lower, Zcoord_Lower, Z2coord_Lower, n_Lower, Thickness_Location));
    
    return Thickness_Value;
    
}

double CBoundaryGeometry::Compute_Area(double *Plane_P0, double *Plane_Normal, unsigned short iSection, vector<double> &Xcoord_Airfoil, vector<double> &Ycoord_Airfoil, vector<double> &Zcoord_Airfoil, bool original_surface) {
	unsigned long iVertex, jVertex;
	double Normal[3], Tangent[3], BiNormal[3], auxXCoord, auxYCoord, auxZCoord, Area_Value = 0.0, Length, Xcoord_Trailing, Ycoord_Trailing, Zcoord_Trailing, ValCos, ValSin, XValue, ZValue;
	vector<double> Xcoord_Upper, Ycoord_Upper, Zcoord_Upper, Xcoord_Lower, Ycoord_Lower, Zcoord_Lower, Z2coord, Xcoord_Normal, Ycoord_Normal, Zcoord_Normal, Xcoord_Airfoil_, Ycoord_Airfoil_, Zcoord_Airfoil_;
    unsigned long Trailing_Point, LeadingPoint;
	double MaxDistance, Distance, AoA;
    
    /*--- Find the leading and trailing edges and compute the angle of attack ---*/
    MaxDistance = 0.0; Trailing_Point = 0;
    for (iVertex = 1; iVertex < Xcoord_Airfoil.size(); iVertex++) {
        
        Distance = sqrt(pow(Xcoord_Airfoil[iVertex] - Xcoord_Airfoil[Trailing_Point], 2.0) +
                        pow(Ycoord_Airfoil[iVertex] - Ycoord_Airfoil[Trailing_Point], 2.0) +
                        pow(Zcoord_Airfoil[iVertex] - Zcoord_Airfoil[Trailing_Point], 2.0));
        
        if (MaxDistance < Distance) { MaxDistance = Distance; LeadingPoint = iVertex; }
    }
    
    AoA = atan((Zcoord_Airfoil[LeadingPoint] - Zcoord_Airfoil[Trailing_Point]) / (Xcoord_Airfoil[Trailing_Point] - Xcoord_Airfoil[LeadingPoint]))*180/PI_NUMBER;
    
    /*--- Translate to the origin ---*/
    Xcoord_Trailing = Xcoord_Airfoil[0];
    Ycoord_Trailing = Ycoord_Airfoil[0];
    Zcoord_Trailing = Zcoord_Airfoil[0];
    
    for (iVertex = 0; iVertex < Xcoord_Airfoil.size(); iVertex++) {
        Xcoord_Airfoil_.push_back(Xcoord_Airfoil[iVertex] - Xcoord_Trailing);
        Ycoord_Airfoil_.push_back(Ycoord_Airfoil[iVertex] - Ycoord_Trailing);
        Zcoord_Airfoil_.push_back(Zcoord_Airfoil[iVertex] - Zcoord_Trailing);
    }
    
    /*--- Rotate the airfoil ---*/
    ValCos = cos(AoA*PI_NUMBER/180.0);
    ValSin = sin(AoA*PI_NUMBER/180.0);
    
    for (iVertex = 0; iVertex < Xcoord_Airfoil.size(); iVertex++) {
        XValue = Xcoord_Airfoil_[iVertex];
        ZValue = Zcoord_Airfoil_[iVertex];
        
        Xcoord_Airfoil_[iVertex] = XValue*ValCos - ZValue*ValSin;
        Zcoord_Airfoil_[iVertex] = ZValue*ValCos + XValue*ValSin;
        
    }
    
    /*--- Identify upper and lower side, and store the value of the normal --*/
    for (iVertex = 1; iVertex < Xcoord_Airfoil_.size(); iVertex++) {
        Tangent[0] = Xcoord_Airfoil_[iVertex] - Xcoord_Airfoil_[iVertex-1];
        Tangent[1] = Ycoord_Airfoil_[iVertex] - Ycoord_Airfoil_[iVertex-1];
        Tangent[2] = Zcoord_Airfoil_[iVertex] - Zcoord_Airfoil_[iVertex-1];
        Length = sqrt(pow(Tangent[0], 2.0) + pow(Tangent[1], 2.0) + pow(Tangent[2], 2.0));
        Tangent[0] /= Length; Tangent[1] /= Length; Tangent[2] /= Length;
        
        BiNormal[0] = Plane_Normal[0];
        BiNormal[1] = Plane_Normal[1];
        BiNormal[2] = Plane_Normal[2];
        Length = sqrt(pow(BiNormal[0], 2.0) + pow(BiNormal[1], 2.0) + pow(BiNormal[2], 2.0));
        BiNormal[0] /= Length; BiNormal[1] /= Length; BiNormal[2] /= Length;
        
        Normal[0] = Tangent[1]*BiNormal[2] - Tangent[2]*BiNormal[1];
        Normal[1] = Tangent[2]*BiNormal[0] - Tangent[0]*BiNormal[2];
        Normal[2] = Tangent[0]*BiNormal[1] - Tangent[1]*BiNormal[0];
        
        Xcoord_Normal.push_back(Normal[0]); Ycoord_Normal.push_back(Normal[1]); Zcoord_Normal.push_back(Normal[2]);
        
        if (Normal[2] >= 0.0) {
            Xcoord_Upper.push_back(Xcoord_Airfoil_[iVertex]);
            Ycoord_Upper.push_back(Ycoord_Airfoil_[iVertex]);
            Zcoord_Upper.push_back(Zcoord_Airfoil_[iVertex]);
        }
        else {
            Xcoord_Lower.push_back(Xcoord_Airfoil_[iVertex]);
            Ycoord_Lower.push_back(Ycoord_Airfoil_[iVertex]);
            Zcoord_Lower.push_back(Zcoord_Airfoil_[iVertex]);
        }
        
    }
    
    /*--- Order the arrays using the X component ---*/
    for (iVertex = 0; iVertex < Xcoord_Upper.size(); iVertex++) {
        for (jVertex = 0; jVertex < Xcoord_Upper.size() - 1 - iVertex; jVertex++) {
            if (Xcoord_Upper[jVertex] > Xcoord_Upper[jVertex+1]) {
                auxXCoord = Xcoord_Upper[jVertex]; Xcoord_Upper[jVertex] = Xcoord_Upper[jVertex+1]; Xcoord_Upper[jVertex+1] = auxXCoord;
                auxYCoord = Ycoord_Upper[jVertex]; Ycoord_Upper[jVertex] = Ycoord_Upper[jVertex+1]; Ycoord_Upper[jVertex+1] = auxYCoord;
                auxZCoord = Zcoord_Upper[jVertex]; Zcoord_Upper[jVertex] = Zcoord_Upper[jVertex+1]; Zcoord_Upper[jVertex+1] = auxZCoord;
            }
        }
    }
    
    /*--- Order the arrays using the X component ---*/
    for (iVertex = 0; iVertex < Xcoord_Lower.size(); iVertex++) {
        for (jVertex = 0; jVertex < Xcoord_Lower.size() - 1 - iVertex; jVertex++) {
            if (Xcoord_Lower[jVertex] > Xcoord_Lower[jVertex+1]) {
                auxXCoord = Xcoord_Lower[jVertex]; Xcoord_Lower[jVertex] = Xcoord_Lower[jVertex+1]; Xcoord_Lower[jVertex+1] = auxXCoord;
                auxYCoord = Ycoord_Lower[jVertex]; Ycoord_Lower[jVertex] = Ycoord_Lower[jVertex+1]; Ycoord_Lower[jVertex+1] = auxYCoord;
                auxZCoord = Zcoord_Lower[jVertex]; Zcoord_Lower[jVertex] = Zcoord_Lower[jVertex+1]; Zcoord_Lower[jVertex+1] = auxZCoord;
            }
        }
    }
    
    /*--- Compute total area ---*/
    Area_Value = 0.0;
    for (iVertex = 0; iVertex < Xcoord_Upper.size()-1; iVertex++)
        Area_Value += fabs((Xcoord_Upper[iVertex+1] - Xcoord_Upper[iVertex]) * 0.5*(Zcoord_Upper[iVertex+1] + Zcoord_Upper[iVertex]));
    for (iVertex = 0; iVertex < Xcoord_Lower.size()-1; iVertex++)
        Area_Value += fabs((Xcoord_Lower[iVertex+1] - Xcoord_Lower[iVertex]) * 0.5*(Zcoord_Lower[iVertex+1] + Zcoord_Lower[iVertex]));
    
    return Area_Value;
    
}

CDomainGeometry::CDomainGeometry(CGeometry *geometry, CConfig *config) {
    
}

CDomainGeometry::~CDomainGeometry(void) {
    
    delete[] Global_to_Local_Point;
	delete[] Local_to_Global_Point;
    delete[] Global_to_Local_Marker;
	delete[] Local_to_Global_Marker;
    
}

void CDomainGeometry::SetDomainSerial(CGeometry *geometry, CConfig *config, unsigned short val_domain) {
	unsigned long iElemDomain, iPointDomain, iPointGhost, iPointReal, iPointPeriodic,
	nVertexDomain[MAX_NUMBER_MARKER], iPoint, iElem, iVertex,
	nelem_edge = 0, nelem_triangle = 0, nelem_quad = 0, nelem_tetra = 0,
	nelem_hexa = 0, nelem_wedge = 0, nelem_pyramid = 0, nPointPeriodic;
	long vnodes_global[8], vnodes_local[8], jPoint;
	unsigned short iNode, iDim, iMarker, nMarkerDomain, nDomain, iDomain, jDomain, jNode;
	double coord[3];
    
    vector<unsigned long> SendDomain_Periodic[MAX_NUMBER_DOMAIN];
    vector<unsigned long> ReceivedDomain_Periodic[MAX_NUMBER_DOMAIN];
    
    vector<unsigned short> SendDomain_PeriodicTrans[MAX_NUMBER_DOMAIN];
    vector<unsigned short> ReceivedDomain_PeriodicTrans[MAX_NUMBER_DOMAIN];
    
    unsigned long ReceptorColor, DonorColor;
    unsigned short jMarker;
    unsigned long Transformation;
    
	vector<unsigned long>::iterator it;
    bool *ElemIn = new bool [geometry->GetnElem()];
	bool *MarkerIn = new bool [geometry->GetnMarker()];
	bool **VertexIn = new bool* [geometry->GetnMarker()];
	for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
		VertexIn[iMarker] = new bool [geometry->GetnElem_Bound(iMarker)];
    
    /*--- Create a copy of the markers ---*/
    Marker_All_SendRecv = new short[100];
    for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
        Marker_All_SendRecv[iMarker] = config->GetMarker_All_SendRecv(iMarker);
    
	nDomain = config->GetnDomain();
	nDim = geometry->GetnDim();
    
	/*--- Auxiliar vector to change the numbering ---*/
	Global_to_Local_Point =  new long[geometry->GetnPoint()];
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
		Global_to_Local_Point[iPoint] = -1;
    
	/*--- Loop over the original grid to perform the dimensionalizaton of the new vectors ---*/
	nElem = 0; nPoint = 0; nPointGhost = 0; nPointDomain = 0; nPointPeriodic = 0;
    
	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
        
		/*--- Check if the element belong to the domain ---*/
		ElemIn[iElem] = false;
		for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++) {
			iPoint = geometry->elem[iElem]->GetNode(iNode);
			if ( geometry->node[iPoint]->GetColor() == val_domain ) {
				ElemIn[iElem] = true; break;
			}
		}
        
		/*--- If an element belong to the domain (at least one point belong has the
		 same color as the domain)---*/
		if (ElemIn[iElem]) {
			for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++) {
				iPoint = geometry->elem[iElem]->GetNode(iNode);
				if (Global_to_Local_Point[iPoint] == -1) {
					Global_to_Local_Point[iPoint] = 1;
                    
					nPoint++;
                    
					if ( geometry->node[iPoint]->GetColor() != val_domain ) nPointGhost++;
					else {
                        /*--- Note that the periodic BC (receive) are also ghost cell ---*/
                        if (iPoint > geometry->GetnPointDomain() - 1) {
                            nPointGhost++;
                            nPointPeriodic++;
                        }
                        else nPointDomain++;
                    }
                    
				}
			}
			nElem++;
		}
	}
    
	/*--- Auxiliar vector to store the local to global index ---*/
	Local_to_Global_Point =  new unsigned long[nPoint];
    
	/*--- Dimensionate the number of elements and nodes of the domain ---*/
	elem = new CPrimalGrid*[nElem];
	node = new CPoint*[nPoint];
    
	/*--- Reset auxiliar vector to change the numbering ---*/
	iElemDomain = 0;
    iPointDomain = 0;
    iPointReal = 0;
    iPointPeriodic = nPointDomain;
    iPointGhost = nPointDomain + nPointPeriodic;
    
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
		Global_to_Local_Point[iPoint] = -1;
    
	/*--- Loop over the original grid to create the point and element structures ---*/
	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
        
		if (ElemIn[iElem]) {
            
			for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++) {
                
				iPoint = geometry->elem[iElem]->GetNode(iNode);
                
				if (Global_to_Local_Point[iPoint] == -1) {
                    
                    if ( geometry->node[iPoint]->GetColor() == val_domain ) {
                        if ( iPoint > geometry->GetnPointDomain() - 1) iPointDomain = iPointPeriodic;
                        else iPointDomain = iPointReal;
                    }
					else iPointDomain = iPointGhost;
                    
					Global_to_Local_Point[iPoint] = iPointDomain;
                    
					for (iDim = 0; iDim < nDim; iDim++)
						coord[iDim] = geometry->node[iPoint]->GetCoord(iDim);
                    
					Local_to_Global_Point[iPointDomain] = iPoint;
					
					/*--- Create the point, and save the color information ---*/
					if ( nDim == 2 ) node[iPointDomain] = new CPoint(coord[0], coord[1], iPoint, config);
					if ( nDim == 3 ) node[iPointDomain] = new CPoint(coord[0], coord[1], coord[2], iPoint, config);
                    node[iPointDomain]->SetColor(geometry->node[iPoint]->GetColor());
                    
                    /*--- Compute the index for the PEriodic and extra Ghost points. ---*/
                    if ( geometry->node[iPoint]->GetColor() == val_domain ) {
                        if ( iPoint > geometry->GetnPointDomain() - 1) iPointPeriodic++;
                        else iPointReal++;
                    }
                    else iPointGhost++;
                    
				}
                
				vnodes_local[iNode] = Global_to_Local_Point[iPoint];
                
			}
            
            /*--- Create the elements ---*/
			switch(geometry->elem[iElem]->GetVTK_Type()) {
                case TRIANGLE:
                    elem[iElemDomain] = new CTriangle(vnodes_local[0], vnodes_local[1], vnodes_local[2], 2);
                    nelem_triangle++;
                    break;
                case RECTANGLE:
                    elem[iElemDomain] = new CRectangle(vnodes_local[0], vnodes_local[1], vnodes_local[2],
                                                       vnodes_local[3], 2);
                    nelem_quad++;
                    break;
                case TETRAHEDRON:
                    elem[iElemDomain] = new CTetrahedron(vnodes_local[0], vnodes_local[1], vnodes_local[2],
                                                         vnodes_local[3]);
                    nelem_tetra++;
                    break;
                case HEXAHEDRON:
                    elem[iElemDomain] = new CHexahedron(vnodes_local[0], vnodes_local[1], vnodes_local[2],
                                                        vnodes_local[3], vnodes_local[4], vnodes_local[5],
                                                        vnodes_local[6], vnodes_local[7]);
                    nelem_hexa++;
                    break;
                case WEDGE:
                    elem[iElemDomain] = new CWedge(vnodes_local[0], vnodes_local[1], vnodes_local[2],
                                                   vnodes_local[3], vnodes_local[4], vnodes_local[5]);
                    nelem_wedge++;
                    break;
                case PYRAMID:
                    elem[iElemDomain] = new CPyramid(vnodes_local[0], vnodes_local[1], vnodes_local[2],
                                                     vnodes_local[3], vnodes_local[4]);
                    nelem_pyramid++;
                    break;
			}
			iElemDomain++;
		}
	}
    
	SetnElem(nElem); SetnPoint(nPoint);
    
	/*--- Dimensionalization with physical boundaries ---*/
	nMarkerDomain = 0;
	for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
        if (config->GetMarker_All_Boundary(iMarker) != SEND_RECEIVE) {
            MarkerIn[iMarker] = false; nVertexDomain[nMarkerDomain] = 0;
            for (iVertex = 0; iVertex < geometry->GetnElem_Bound(iMarker); iVertex++) {
                VertexIn[iMarker][iVertex] = false;
                for (iNode = 0; iNode < geometry->bound[iMarker][iVertex]->GetnNodes(); iNode++) {
                    vnodes_global[iNode] = geometry->bound[iMarker][iVertex]->GetNode(iNode);
                    vnodes_local[iNode] = Global_to_Local_Point[vnodes_global[iNode]];
                    if (geometry->node[vnodes_global[iNode]]->GetColor() == val_domain ) VertexIn[iMarker][iVertex] = true;
                }
                if (VertexIn[iMarker][iVertex]) { nVertexDomain[nMarkerDomain] ++;  MarkerIn[iMarker] = true; }
            }
            if (MarkerIn[iMarker]) { nMarkerDomain++; }
        }
    }
    
    /*--- Evaluate the number of already existing periodic boundary conditions ---*/
    for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
        
        if (config->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) {
            
            for (iVertex = 0; iVertex < geometry->GetnElem_Bound(iMarker); iVertex++) {
                iPoint = geometry->bound[iMarker][iVertex]->GetNode(0);
                Transformation = geometry->bound[iMarker][iVertex]->GetRotation_Type();
                
                if (val_domain == geometry->node[iPoint]->GetColor()) {
                    
                    /*--- If the information is going to be sended, find the
                     domain of the receptor ---*/
                    if (Marker_All_SendRecv[iMarker] > 0) {
                        
                        /*--- Identify the color of the receptor ---*/
                        for (jMarker = 0; jMarker < geometry->GetnMarker(); jMarker++) {
                            if (Marker_All_SendRecv[jMarker] < 0) {
                                jPoint = geometry->bound[jMarker][iVertex]->GetNode(0);
                                ReceptorColor = geometry->node[jPoint]->GetColor();
                            }
                        }
                        
                        /*--- For each color of the receptor we will han an extra marker (+) ---*/
                        SendDomain_Periodic[ReceptorColor].push_back(Global_to_Local_Point[iPoint]);
                        SendDomain_PeriodicTrans[ReceptorColor].push_back(Transformation);
                    }
                    
                    
                    /*--- If the information is goint to be received,
                     find the domain if the donor ---*/
                    if (Marker_All_SendRecv[iMarker] < 0) {
                        
                        /*--- Identify the color of the donor ---*/
                        for (jMarker = 0; jMarker < geometry->GetnMarker(); jMarker++) {
                            if (Marker_All_SendRecv[jMarker] > 0) {
                                jPoint = geometry->bound[jMarker][iVertex]->GetNode(0);
                                DonorColor = geometry->node[jPoint]->GetColor();
                            }
                        }
                        
                        /*--- For each color of the donor we will han an extra marker (-) ---*/
                        ReceivedDomain_Periodic[DonorColor].push_back(Global_to_Local_Point[iPoint]);
                        ReceivedDomain_PeriodicTrans[DonorColor].push_back(Transformation);
                        
                    }
                }
            }
        }
    }
    
    for (iDomain = 0; iDomain < nDomain; iDomain++) {
        
        /*--- Add the new periodic markers to the domain ---*/
        if (SendDomain_Periodic[iDomain].size() != 0) {
            nVertexDomain[nMarkerDomain] = SendDomain_Periodic[iDomain].size();
            nMarkerDomain++;
        }
        if (ReceivedDomain_Periodic[iDomain].size() != 0) {
            nVertexDomain[nMarkerDomain] = ReceivedDomain_Periodic[iDomain].size();
            nMarkerDomain++;
        }
        
    }
    
	/*--- Loop over the all the points of the element
	 to find the points with different colours, and create the send/received list ---*/
	for (iElem = 0; iElem < nElem; iElem++) {
		for (iNode = 0; iNode < elem[iElem]->GetnNodes(); iNode++) {
			iPoint = elem[iElem]->GetNode(iNode);
            iDomain = node[iPoint]->GetColor();
            
            if (iDomain == val_domain) {
                for(jNode = 0; jNode < elem[iElem]->GetnNodes(); jNode++) {
                    jPoint = elem[iElem]->GetNode(jNode);
                    jDomain = node[jPoint]->GetColor();
                    
                    /*--- If different color and connected by an edge, then we add them to the list ---*/
                    if (iDomain != jDomain) {
                        
                        /*--- We send from iDomain to jDomain the value of iPoint, we save the
                         global value becuase we need to sort the lists ---*/
                        SendDomainLocal[jDomain].push_back(Local_to_Global_Point[iPoint]);
                        /*--- We send from jDomain to iDomain the value of jPoint, we save the
                         global value becuase we need to sort the lists ---*/
                        ReceivedDomainLocal[jDomain].push_back(Local_to_Global_Point[jPoint]);
                        
                    }
                }
            }
		}
	}
	
	/*--- Sort the points that must be sended and delete repeated points, note
     that the sortering should be done with the global point (not the local) ---*/
	for (iDomain = 0; iDomain < nDomain; iDomain++) {
		sort( SendDomainLocal[iDomain].begin(), SendDomainLocal[iDomain].end());
		it = unique( SendDomainLocal[iDomain].begin(), SendDomainLocal[iDomain].end());
		SendDomainLocal[iDomain].resize( it - SendDomainLocal[iDomain].begin() );
	}
	
	/*--- Sort the points that must be received and delete repeated points, note
     that the sortering should be done with the global point (not the local) ---*/
	for (iDomain = 0; iDomain < nDomain; iDomain++) {
		sort( ReceivedDomainLocal[iDomain].begin(), ReceivedDomainLocal[iDomain].end());
		it = unique( ReceivedDomainLocal[iDomain].begin(), ReceivedDomainLocal[iDomain].end());
		ReceivedDomainLocal[iDomain].resize( it - ReceivedDomainLocal[iDomain].begin() );
	}
    
    /*--- Add the new MPI send receive boundaries, reset the transformation, and save the local value ---*/
	for (iDomain = 0; iDomain < nDomain; iDomain++) {
		if (SendDomainLocal[iDomain].size() != 0) {
			nVertexDomain[nMarkerDomain] = SendDomainLocal[iDomain].size();
            for (iVertex = 0; iVertex < nVertexDomain[nMarkerDomain]; iVertex++) {
                SendDomainLocal[iDomain][iVertex] = Global_to_Local_Point[SendDomainLocal[iDomain][iVertex]];
                SendTransfLocal[iDomain].push_back(0);
            }
			nMarkerDomain++;
		}
    }
    
    /*--- Add the new MPI receive boundaries, reset the transformation, and save the local value ---*/
    for (iDomain = 0; iDomain < nDomain; iDomain++) {
 		if (ReceivedDomainLocal[iDomain].size() != 0) {
			nVertexDomain[nMarkerDomain] = ReceivedDomainLocal[iDomain].size();
            for (iVertex = 0; iVertex < nVertexDomain[nMarkerDomain]; iVertex++) {
                ReceivedDomainLocal[iDomain][iVertex] = Global_to_Local_Point[ReceivedDomainLocal[iDomain][iVertex]];
                ReceivedTransfLocal[iDomain].push_back(0);
            }
			nMarkerDomain++;
		}
    }
    
	SetnMarker(nMarkerDomain);
	nElem_Bound = new unsigned long [nMarkerDomain];
	Local_to_Global_Marker = new unsigned short [nMarkerDomain];
	Global_to_Local_Marker = new unsigned short [geometry->GetnMarker()];
    
	for (iMarker = 0; iMarker < nMarkerDomain; iMarker++)
		SetnElem_Bound(iMarker, nVertexDomain[iMarker]);
    
	bound = new CPrimalGrid**[GetnMarker()];
	for (iMarker = 0; iMarker < GetnMarker(); iMarker++)
		bound[iMarker] = new CPrimalGrid* [GetnElem_Bound(iMarker)];
    
	/*--- Loop over the original grid to create the boundaries (hre we need to add the already existing send receive, periodic bc?) ---*/
	nMarkerDomain = 0;
	for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
        if (config->GetMarker_All_Boundary(iMarker) != SEND_RECEIVE) {
            
            /*--- If the marker is in the domain ---*/
            if (MarkerIn[iMarker]) {
                
                nelem_edge = 0; nelem_triangle = 0; nelem_quad = 0;
                nVertexDomain[nMarkerDomain] = 0;
                
                for (iVertex = 0; iVertex < geometry->GetnElem_Bound(iMarker); iVertex++) {
                    
                    /*--- If the vertex is in the domain ---*/
                    if (VertexIn[iMarker][iVertex]) {
                        
                        /*--- Read the points in the local domain ---*/
                        for (iNode = 0; iNode < geometry->bound[iMarker][iVertex]->GetnNodes(); iNode++) {
                            vnodes_global[iNode] = geometry->bound[iMarker][iVertex]->GetNode(iNode);
                            vnodes_local[iNode] = Global_to_Local_Point[vnodes_global[iNode]];
                        }
                        
                        /*--- Create the data structure for the boundaries ---*/
                        switch(geometry->bound[iMarker][iVertex]->GetVTK_Type()) {
                            case LINE:
                                bound[nMarkerDomain][nVertexDomain[nMarkerDomain]] = new CLine(vnodes_local[0],vnodes_local[1],2);
                                nelem_edge++;
                                break;
                            case TRIANGLE:
                                bound[nMarkerDomain][nVertexDomain[nMarkerDomain]] = new CTriangle(vnodes_local[0],vnodes_local[1],
                                                                                                   vnodes_local[2],3);
                                nelem_triangle++;
                                break;
                            case RECTANGLE:
                                bound[nMarkerDomain][nVertexDomain[nMarkerDomain]] = new CRectangle(vnodes_local[0],vnodes_local[1],
                                                                                                    vnodes_local[2],vnodes_local[3],3);
                                nelem_quad++;
                                break;
                        }
                        nVertexDomain[nMarkerDomain]++;
                        
                    }
                }
                
                /*--- add the marker and update some structures ---*/
                Local_to_Global_Marker[nMarkerDomain] = iMarker;
                Global_to_Local_Marker[iMarker] = nMarkerDomain;
                nMarkerDomain++;
                
            }
        }
	}
	
    /*--- Add the periodic BC ---*/
    for (iDomain = 0; iDomain < nDomain; iDomain++) {
        
        /*--- Add the new periodic markers to the domain ---*/
        if (SendDomain_Periodic[iDomain].size() != 0) {
            nVertexDomain[nMarkerDomain] = 0;
            for (iVertex = 0; iVertex < SendDomain_Periodic[iDomain].size(); iVertex++) {
                bound[nMarkerDomain][iVertex] = new CVertexMPI(SendDomain_Periodic[iDomain][iVertex], 3);
                bound[nMarkerDomain][iVertex]->SetRotation_Type(SendDomain_PeriodicTrans[iDomain][iVertex]);
                nVertexDomain[nMarkerDomain]++;
            }
            Marker_All_SendRecv[nMarkerDomain] = iDomain+1;
            nMarkerDomain++;
        }
        if (ReceivedDomain_Periodic[iDomain].size() != 0) {
            nVertexDomain[nMarkerDomain] = 0;
            for (iVertex = 0; iVertex < ReceivedDomain_Periodic[iDomain].size(); iVertex++) {
                bound[nMarkerDomain][iVertex] = new CVertexMPI(ReceivedDomain_Periodic[iDomain][iVertex], 3);
                bound[nMarkerDomain][iVertex]->SetRotation_Type(ReceivedDomain_PeriodicTrans[iDomain][iVertex]);
                nVertexDomain[nMarkerDomain]++;
            }
            Marker_All_SendRecv[nMarkerDomain] = -(iDomain+1);
            nMarkerDomain++;
        }
    }
    
	cout << "Domain "<< val_domain + 1 << ": " << nPoint << " points (" << nPointGhost << " ghost points";
	if (nPointPeriodic == 0) cout <<")." << endl;
	else cout <<" including " << nPointPeriodic << " periodic points)." << endl;
    
	delete[] ElemIn;
	delete[] MarkerIn;
	for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
		delete VertexIn[iMarker];
	delete[] VertexIn;
    
}

void CDomainGeometry::SetSendReceive(CConfig *config) {
    
}

void CDomainGeometry::SetMeshFile(CConfig *config, string val_mesh_out_filename) {
	unsigned long iElem, iPoint, iElem_Bound;
	unsigned short iMarker, iNodes, iDim, iPeriodic, nPeriodic = 0;
	double *center, *angles, *transl;
	ofstream output_file;
	string Grid_Marker;
    
    /*--- Open the output file ---*/
	char *cstr = new char [val_mesh_out_filename.size()+1];
	strcpy (cstr, val_mesh_out_filename.c_str());
	output_file.open(cstr, ios::out);
    
	output_file << "NDIME= " << nDim << endl;
	output_file << "NELEM= " << nElem << endl;
    
	for (iElem = 0; iElem < nElem; iElem++) {
		output_file << elem[iElem]->GetVTK_Type();
		for (iNodes = 0; iNodes < elem[iElem]->GetnNodes(); iNodes++)
			output_file << "\t" << elem[iElem]->GetNode(iNodes);
		output_file << "\t"<<iElem<<endl;
	}
    
	output_file << "NPOIN= " << nPoint << "\t" << nPointDomain <<endl;
	output_file.precision(15);
	for (iPoint = 0; iPoint < nPoint; iPoint++) {
		for (iDim = 0; iDim < nDim; iDim++)
			output_file << scientific << "\t" << node[iPoint]->GetCoord(iDim) ;
		output_file << "\t" << iPoint << "\t" << Local_to_Global_Point[iPoint] << endl;
	}
	
    output_file << "NMARK= " << nMarker << endl;
	for (iMarker = 0; iMarker < nMarker; iMarker++) {
		if (bound[iMarker][0]->GetVTK_Type() != VERTEX) {
			Grid_Marker = config->GetMarker_All_Tag(Local_to_Global_Marker[iMarker]);
			output_file << "MARKER_TAG= " << Grid_Marker <<endl;
			output_file << "MARKER_ELEMS= " << nElem_Bound[iMarker]<< endl;
			
			for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
				output_file << bound[iMarker][iElem_Bound]->GetVTK_Type() << "\t" ;
				for (iNodes = 0; iNodes < bound[iMarker][iElem_Bound]->GetnNodes()-1; iNodes++)
					output_file << bound[iMarker][iElem_Bound]->GetNode(iNodes) << "\t" ;
				iNodes = bound[iMarker][iElem_Bound]->GetnNodes()-1;
				output_file << bound[iMarker][iElem_Bound]->GetNode(iNodes) << endl;
			}
        }
    }
    
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
		if (bound[iMarker][0]->GetVTK_Type() == VERTEX) {
            output_file << "MARKER_TAG= SEND_RECEIVE" << endl;
            output_file << "MARKER_ELEMS= " << nElem_Bound[iMarker]<< endl;
            if (Marker_All_SendRecv[iMarker] > 0) output_file << "SEND_TO= " << Marker_All_SendRecv[iMarker] << endl;
            if (Marker_All_SendRecv[iMarker] < 0) output_file << "SEND_TO= " << Marker_All_SendRecv[iMarker] << endl;
            
            for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
                output_file << bound[iMarker][iElem_Bound]->GetVTK_Type() << "\t" ;
                output_file << bound[iMarker][iElem_Bound]->GetNode(0) << "\t";
                output_file << bound[iMarker][iElem_Bound]->GetRotation_Type() << "\t";
                output_file << bound[iMarker][iElem_Bound]->GetMatching_Zone() << endl;
            }
        }
    }
    
	/*--- Get the total number of periodic transformations ---*/
	nPeriodic = config->GetnPeriodicIndex();
	output_file << "NPERIODIC= " << nPeriodic << endl;
    
	/*--- From iPeriodic obtain the iMarker ---*/
	for (iPeriodic = 0; iPeriodic < nPeriodic; iPeriodic++) {
        
		/*--- Retrieve the supplied periodic information. ---*/
		center = config->GetPeriodicCenter(iPeriodic);
		angles = config->GetPeriodicRotation(iPeriodic);
		transl = config->GetPeriodicTranslate(iPeriodic);
        
		output_file << "PERIODIC_INDEX= " << iPeriodic << endl;
		output_file << center[0] << "\t" << center[1] << "\t" << center[2] << endl;
		output_file << angles[0] << "\t" << angles[1] << "\t" << angles[2] << endl;
		output_file << transl[0] << "\t" << transl[1] << "\t" << transl[2] << endl;
        
	}
    
	output_file.close();
    
    delete[] cstr;
}

void CDomainGeometry::SetTecPlot(char mesh_filename[200]) {
	unsigned long iElem, iPoint;
	unsigned short iDim;
	ofstream Tecplot_File;
    
	Tecplot_File.open(mesh_filename, ios::out);
	Tecplot_File << "TITLE = \"Visualization of the volumetric grid\"" << endl;
    
	if (nDim == 2) {
		Tecplot_File << "VARIABLES = \"x\",\"y\" " << endl;
		Tecplot_File << "ZONE NODES= "<< nPoint <<", ELEMENTS= "<< nElem <<", DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL"<< endl;
	}
	if (nDim == 3) {
		Tecplot_File << "VARIABLES = \"x\",\"y\",\"z\" " << endl;
		Tecplot_File << "ZONE NODES= "<< nPoint <<", ELEMENTS= "<< nElem <<", DATAPACKING=POINT, ZONETYPE=FEBRICK"<< endl;
	}
    
	for(iPoint = 0; iPoint < nPoint; iPoint++) {
		for(iDim = 0; iDim < nDim; iDim++)
			Tecplot_File << scientific << node[iPoint]->GetCoord(iDim) << "\t";
		Tecplot_File << "\n";
	}
    
	for(iElem = 0; iElem < nElem; iElem++) {
		if (elem[iElem]->GetVTK_Type() == TRIANGLE) {
			Tecplot_File <<
            elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
            elem[iElem]->GetNode(2)+1 <<" "<< elem[iElem]->GetNode(2)+1 << endl;
		}
		if (elem[iElem]->GetVTK_Type() == RECTANGLE) {
			Tecplot_File <<
            elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
            elem[iElem]->GetNode(2)+1 <<" "<< elem[iElem]->GetNode(3)+1 << endl;
		}
		if (elem[iElem]->GetVTK_Type() == TETRAHEDRON) {
			Tecplot_File <<
            elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
            elem[iElem]->GetNode(2)+1 <<" "<< elem[iElem]->GetNode(2)+1 <<" "<<
            elem[iElem]->GetNode(3)+1 <<" "<< elem[iElem]->GetNode(3)+1 <<" "<<
            elem[iElem]->GetNode(3)+1 <<" "<< elem[iElem]->GetNode(3)+1 << endl;
		}
		if (elem[iElem]->GetVTK_Type() == HEXAHEDRON) {
			Tecplot_File <<
            elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
            elem[iElem]->GetNode(2)+1 <<" "<< elem[iElem]->GetNode(3)+1 <<" "<<
            elem[iElem]->GetNode(4)+1 <<" "<< elem[iElem]->GetNode(5)+1 <<" "<<
            elem[iElem]->GetNode(6)+1 <<" "<< elem[iElem]->GetNode(7)+1 << endl;
		}
		if (elem[iElem]->GetVTK_Type() == PYRAMID) {
			Tecplot_File <<
            elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
            elem[iElem]->GetNode(2)+1 <<" "<< elem[iElem]->GetNode(3)+1 <<" "<<
            elem[iElem]->GetNode(4)+1 <<" "<< elem[iElem]->GetNode(4)+1 <<" "<<
            elem[iElem]->GetNode(4)+1 <<" "<< elem[iElem]->GetNode(4)+1 << endl;
		}
		if (elem[iElem]->GetVTK_Type() == WEDGE) {
			Tecplot_File <<
            elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
            elem[iElem]->GetNode(1)+1 <<" "<< elem[iElem]->GetNode(2)+1 <<" "<<
            elem[iElem]->GetNode(3)+1 <<" "<< elem[iElem]->GetNode(4)+1 <<" "<<
            elem[iElem]->GetNode(4)+1 <<" "<< elem[iElem]->GetNode(5)+1 << endl;
		}
	}
    
	Tecplot_File.close();
    
}

CPeriodicGeometry::CPeriodicGeometry(CGeometry *geometry, CConfig *config) {
	unsigned long nElem_new, nPoint_new, jPoint, iPoint, iElem, jElem, iVertex, 
	nelem_triangle = 0, nelem_quad = 0, nelem_tetra = 0, nelem_hexa = 0, nelem_wedge = 0, 
	nelem_pyramid = 0, iIndex, newElementsBound = 0;
	unsigned short  iMarker, nPeriodic = 0, iPeriodic;
	double *center, *angles, rotMatrix[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}, 
    translation[3], *trans, theta, phi, psi, cosTheta, sinTheta, cosPhi, sinPhi, cosPsi, sinPsi,
    dx, dy, dz, rotCoord[3], *Coord_i;
    
    /*--- It only create the mirror structure for the second boundary ---*/
    bool CreateMirror[10];
    CreateMirror[1] = false;
    CreateMirror[2] = true;
    
	/*--- Compute the number of periodic bc on the geometry ---*/
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if (config->GetMarker_All_Boundary(iMarker) == PERIODIC_BOUNDARY)
			nPeriodic++;
    
	/*--- Write the number of dimensions of the problem ---*/
	nDim = geometry->GetnDim();
    
	/*--- Copy the new boundary element information from the geometry class.
     Be careful, as these are pointers to vectors/objects. ---*/
	nNewElem_BoundPer = geometry->nNewElem_Bound;
	newBoundPer       = geometry->newBound;
    
	/*--- Count the number of new boundary elements. ---*/
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		newElementsBound += nNewElem_BoundPer[iMarker];
    
	/*--- Loop over the original grid to perform the dimensionalizaton of the new vectors ---*/
	nElem_new = 0; nPoint_new = 0;
	for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++) {
        if (CreateMirror[iPeriodic]) {
            nElem_new += geometry->PeriodicElem[iPeriodic].size();
            nPoint_new += geometry->PeriodicPoint[iPeriodic][0].size();
        }
    }
    
	cout << "Number of new points: " << nPoint_new << "." << endl;
	cout << "Number of new interior elements: " << nElem_new << "." << endl;
	cout << "Number of new boundary elements added to preexisting markers: " << newElementsBound << "." << endl;
    
	/*--- Create a copy of the original grid ---*/
	elem = new CPrimalGrid*[geometry->GetnElem() + nElem_new];
	for (iElem = 0; iElem < geometry->GetnElem(); iElem ++) {
		switch(geometry->elem[iElem]->GetVTK_Type()) {
            case TRIANGLE:
                elem[iElem] = new CTriangle(geometry->elem[iElem]->GetNode(0),
                                            geometry->elem[iElem]->GetNode(1),
                                            geometry->elem[iElem]->GetNode(2), 2);
                nelem_triangle++;
                break;
                
            case RECTANGLE:
                elem[iElem] = new CRectangle(geometry->elem[iElem]->GetNode(0),
                                             geometry->elem[iElem]->GetNode(1),
                                             geometry->elem[iElem]->GetNode(2),
                                             geometry->elem[iElem]->GetNode(3), 2);
                nelem_quad++;
                break;
                
            case TETRAHEDRON:
                elem[iElem] = new CTetrahedron(geometry->elem[iElem]->GetNode(0),
                                               geometry->elem[iElem]->GetNode(1),
                                               geometry->elem[iElem]->GetNode(2),
                                               geometry->elem[iElem]->GetNode(3));
                nelem_tetra++;
                break;
                
            case HEXAHEDRON:
                elem[iElem] = new CHexahedron(geometry->elem[iElem]->GetNode(0),
                                              geometry->elem[iElem]->GetNode(1),
                                              geometry->elem[iElem]->GetNode(2),
                                              geometry->elem[iElem]->GetNode(3),
                                              geometry->elem[iElem]->GetNode(4),
                                              geometry->elem[iElem]->GetNode(5),
                                              geometry->elem[iElem]->GetNode(6),
                                              geometry->elem[iElem]->GetNode(7));
                nelem_hexa++;
                break;
                
            case WEDGE:
                elem[iElem] = new CWedge(geometry->elem[iElem]->GetNode(0),
                                         geometry->elem[iElem]->GetNode(1),
                                         geometry->elem[iElem]->GetNode(2),
                                         geometry->elem[iElem]->GetNode(3),
                                         geometry->elem[iElem]->GetNode(4),
                                         geometry->elem[iElem]->GetNode(5));
                nelem_wedge++;
                break;
                
            case PYRAMID:
                elem[iElem] = new CPyramid(geometry->elem[iElem]->GetNode(0),
                                           geometry->elem[iElem]->GetNode(1),
                                           geometry->elem[iElem]->GetNode(2),
                                           geometry->elem[iElem]->GetNode(3),
                                           geometry->elem[iElem]->GetNode(4));
                nelem_pyramid++;
                break;
                
		}
	}
    
	/*--- Create a list with all the points and the new index ---*/
	unsigned long *Index = new unsigned long [geometry->GetnPoint()];
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) Index[iPoint] = 0;
    
	for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++) {
        if (CreateMirror[iPeriodic]) {
            for (iIndex = 0; iIndex < geometry->PeriodicPoint[iPeriodic][0].size(); iIndex++) {
                iPoint =  geometry->PeriodicPoint[iPeriodic][0][iIndex];
                Index[iPoint] = geometry->PeriodicPoint[iPeriodic][1][iIndex];
            }
        }
    }
    
	for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
		if (config->GetMarker_All_Boundary(iMarker) == PERIODIC_BOUNDARY)
			for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				jPoint = geometry->vertex[iMarker][iVertex]->GetDonorPoint();
				Index[iPoint] = jPoint;
			}
    
	/*--- Add the new elements due to the periodic boundary condtion ---*/
	iElem = geometry->GetnElem();
    
	for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++) {
        if (CreateMirror[iPeriodic]) {
            for (iIndex = 0; iIndex < geometry->PeriodicElem[iPeriodic].size(); iIndex++) {
                jElem = geometry->PeriodicElem[iPeriodic][iIndex];
                
                switch(geometry->elem[jElem]->GetVTK_Type()) {
                    case TRIANGLE:
                        elem[iElem] = new CTriangle(Index[geometry->elem[jElem]->GetNode(0)],
                                                    Index[geometry->elem[jElem]->GetNode(1)],
                                                    Index[geometry->elem[jElem]->GetNode(2)], 2);
                        iElem++; nelem_triangle++;
                        break;
                        
                    case RECTANGLE:
                        elem[iElem] = new CRectangle(Index[geometry->elem[jElem]->GetNode(0)],
                                                     Index[geometry->elem[jElem]->GetNode(1)],
                                                     Index[geometry->elem[jElem]->GetNode(2)],
                                                     Index[geometry->elem[jElem]->GetNode(3)], 2);
                        iElem++; nelem_quad++;
                        break;
                        
                    case TETRAHEDRON:
                        elem[iElem] = new CTetrahedron(Index[geometry->elem[jElem]->GetNode(0)],
                                                       Index[geometry->elem[jElem]->GetNode(1)],
                                                       Index[geometry->elem[jElem]->GetNode(2)],
                                                       Index[geometry->elem[jElem]->GetNode(3)]);
                        iElem++; nelem_tetra++;
                        break;
                        
                    case HEXAHEDRON:
                        elem[iElem] = new CHexahedron(Index[geometry->elem[jElem]->GetNode(0)],
                                                      Index[geometry->elem[jElem]->GetNode(1)],
                                                      Index[geometry->elem[jElem]->GetNode(2)],
                                                      Index[geometry->elem[jElem]->GetNode(3)],
                                                      Index[geometry->elem[jElem]->GetNode(4)],
                                                      Index[geometry->elem[jElem]->GetNode(5)],
                                                      Index[geometry->elem[jElem]->GetNode(6)],
                                                      Index[geometry->elem[jElem]->GetNode(7)]);
                        iElem++; nelem_hexa++;
                        break;
                        
                    case WEDGE:
                        elem[iElem] = new CWedge(Index[geometry->elem[jElem]->GetNode(0)],
                                                 Index[geometry->elem[jElem]->GetNode(1)],
                                                 Index[geometry->elem[jElem]->GetNode(2)],
                                                 Index[geometry->elem[jElem]->GetNode(3)],
                                                 Index[geometry->elem[jElem]->GetNode(4)],
                                                 Index[geometry->elem[jElem]->GetNode(5)]);
                        iElem++; nelem_wedge++;
                        break;
                        
                    case PYRAMID:
                        elem[iElem] = new CPyramid(Index[geometry->elem[jElem]->GetNode(0)],
                                                   Index[geometry->elem[jElem]->GetNode(1)],
                                                   Index[geometry->elem[jElem]->GetNode(2)],
                                                   Index[geometry->elem[jElem]->GetNode(3)],
                                                   Index[geometry->elem[jElem]->GetNode(4)]);
                        iElem++; nelem_pyramid++;
                        break;
                        
                }
            }
        }
    }
    
	nElem = geometry->GetnElem() + nElem_new;
    
	/*--- Add the old points ---*/
	node = new CPoint*[geometry->GetnPoint() + nPoint_new];
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
		if (geometry->GetnDim() == 2)
			node[iPoint] = new CPoint(geometry->node[iPoint]->GetCoord(0),
                                      geometry->node[iPoint]->GetCoord(1), iPoint, config);
		if (geometry->GetnDim() == 3)
			node[iPoint] = new CPoint(geometry->node[iPoint]->GetCoord(0),
                                      geometry->node[iPoint]->GetCoord(1),
                                      geometry->node[iPoint]->GetCoord(2), iPoint, config);
	}
    
	/*--- Add the new points due to the periodic boundary condtion (only in the mirror part) ---*/
	for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++) {
        if (CreateMirror[iPeriodic]) {
            for (iIndex = 0; iIndex < geometry->PeriodicPoint[iPeriodic][0].size(); iIndex++) {
                
                /*--- From iPeriodic obtain the iMarker ---*/
                for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
                    if (iPeriodic == config->GetMarker_All_PerBound(iMarker)) break;
                
                /*--- Retrieve the supplied periodic information. ---*/
                center = config->GetPeriodicRotCenter(config->GetMarker_All_Tag(iMarker));
                angles = config->GetPeriodicRotAngles(config->GetMarker_All_Tag(iMarker));
                trans  = config->GetPeriodicTranslation(config->GetMarker_All_Tag(iMarker));
                
                /*--- Store center - trans as it is constant and will be added on.
                 Note the subtraction, as this is the inverse translation. ---*/
                translation[0] = center[0] - trans[0];
                translation[1] = center[1] - trans[1];
                translation[2] = center[2] - trans[2];
                
                /*--- Store angles separately for clarity. Compute sines/cosines.
                 Note the negative sign, as this is the inverse rotation. ---*/
                theta = -angles[0];
                phi   = -angles[1];
                psi   = -angles[2];
                
                cosTheta = cos(theta);  cosPhi = cos(phi);  cosPsi = cos(psi);
                sinTheta = sin(theta);  sinPhi = sin(phi);  sinPsi = sin(psi);
                
                /*--- Compute the rotation matrix. Note that the implicit
                 ordering is rotation about the x-axis, y-axis, then z-axis. ---*/
                rotMatrix[0][0] = cosPhi*cosPsi;
                rotMatrix[1][0] = cosPhi*sinPsi;
                rotMatrix[2][0] = -sinPhi;
                
                rotMatrix[0][1] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
                rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
                rotMatrix[2][1] = sinTheta*cosPhi;
                
                rotMatrix[0][2] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
                rotMatrix[1][2] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
                rotMatrix[2][2] = cosTheta*cosPhi;
                
                /*--- Retrieve node information for this boundary point. ---*/
                iPoint = geometry->PeriodicPoint[iPeriodic][0][iIndex];
                jPoint = geometry->PeriodicPoint[iPeriodic][1][iIndex];
                Coord_i = geometry->node[iPoint]->GetCoord();
                
                /*--- Get the position vector from rot center to point. ---*/
                dx = Coord_i[0] - center[0];
                dy = Coord_i[1] - center[1];
                if (nDim == 3) {
                    dz = Coord_i[2] - center[2];
                } else {
                    dz = 0.0;
                }
                
                /*--- Compute transformed point coordinates. ---*/
                rotCoord[0] = rotMatrix[0][0]*dx + rotMatrix[0][1]*dy + rotMatrix[0][2]*dz + translation[0];
                rotCoord[1] = rotMatrix[1][0]*dx + rotMatrix[1][1]*dy + rotMatrix[1][2]*dz + translation[1];
                rotCoord[2] = rotMatrix[2][0]*dx + rotMatrix[2][1]*dy + rotMatrix[2][2]*dz + translation[2];
                
                /*--- Save the new points with the new coordinates. ---*/
                if (geometry->GetnDim() == 2)
                    node[jPoint] = new CPoint(rotCoord[0], rotCoord[1], jPoint, config);
                if (geometry->GetnDim() == 3)
                    node[jPoint] = new CPoint(rotCoord[0], rotCoord[1], rotCoord[2], jPoint, config);
                
            }
        }
    }
    
	nPoint = geometry->GetnPoint() + nPoint_new;
    
	/*--- Add the old boundary, reserving space for two new bc (send/recive periodic bc) ---*/
	nMarker = geometry->GetnMarker() + 2;
	nElem_Bound = new unsigned long [nMarker];
	bound = new CPrimalGrid**[nMarker];	
	Tag_to_Marker = new string [MAX_INDEX_VALUE];
	config->SetnMarker_All(nMarker);
    
	/*--- Copy the olf boundary ---*/
	for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
        
		bound[iMarker] = new CPrimalGrid* [geometry->GetnElem_Bound(iMarker)];
        
		for (iVertex = 0; iVertex < geometry->GetnElem_Bound(iMarker); iVertex++) {
			if (geometry->bound[iMarker][iVertex]->GetVTK_Type() == LINE)
				bound[iMarker][iVertex] = new CLine(geometry->bound[iMarker][iVertex]->GetNode(0),
                                                    geometry->bound[iMarker][iVertex]->GetNode(1), 2);
			if (geometry->bound[iMarker][iVertex]->GetVTK_Type() == TRIANGLE)
				bound[iMarker][iVertex] = new CTriangle(geometry->bound[iMarker][iVertex]->GetNode(0),
                                                        geometry->bound[iMarker][iVertex]->GetNode(1),
                                                        geometry->bound[iMarker][iVertex]->GetNode(2), 3);
			if (geometry->bound[iMarker][iVertex]->GetVTK_Type() == RECTANGLE)
				bound[iMarker][iVertex] = new CRectangle(geometry->bound[iMarker][iVertex]->GetNode(0),
                                                         geometry->bound[iMarker][iVertex]->GetNode(1),
                                                         geometry->bound[iMarker][iVertex]->GetNode(2),
                                                         geometry->bound[iMarker][iVertex]->GetNode(3), 3);
		}
        
		nElem_Bound[iMarker] = geometry->GetnElem_Bound(iMarker);
		Tag_to_Marker[iMarker] = geometry->GetMarker_Tag(iMarker);
        
	}
    
	delete[] Index;
    
}

CPeriodicGeometry::~CPeriodicGeometry(void) {
    unsigned long iElem_Bound;
    unsigned short iMarker;
    
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
        for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
            if (newBoundPer[iMarker][iElem_Bound] != NULL) delete newBoundPer[iMarker][iElem_Bound];
        }
    }
    if (newBoundPer != NULL) delete[] newBoundPer;
    
    if (nNewElem_BoundPer != NULL) delete[] nNewElem_BoundPer;
    
}

void CPeriodicGeometry::SetPeriodicBoundary(CGeometry *geometry, CConfig *config) {
	unsigned short iMarker, iPeriodic, nPeriodic = 0, iMarkerSend, iMarkerReceive;
	unsigned long iVertex, Counter_Send = 0, Counter_Receive = 0, iIndex;
    
	/*--- Compute the number of periodic bc on the geometry ---*/
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if (config->GetMarker_All_Boundary(iMarker) == PERIODIC_BOUNDARY)
			nPeriodic++;
    
	/*--- First compute the Send/Receive boundaries, count the number of points ---*/
	Counter_Send = 0; 	Counter_Receive = 0;
	for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++) {
		if (geometry->PeriodicPoint[iPeriodic][0].size() != 0) 
			Counter_Send += geometry->PeriodicPoint[iPeriodic][0].size();
		if (geometry->PeriodicPoint[iPeriodic][1].size() != 0) 
			Counter_Receive += geometry->PeriodicPoint[iPeriodic][1].size();		
	}
    
	/*--- Adimensionalization of the new boundaries ---*/
	iMarkerSend = nMarker - 2; iMarkerReceive = nMarker - 1;
	config->SetMarker_All_SendRecv(iMarkerSend,1);
	config->SetMarker_All_SendRecv(iMarkerReceive,-1);
	nElem_Bound[iMarkerSend] = Counter_Send; 
	nElem_Bound[iMarkerReceive] = Counter_Receive;
	bound[iMarkerSend] = new CPrimalGrid* [Counter_Send];
	bound[iMarkerReceive] = new CPrimalGrid* [Counter_Receive];
    
	/*--- First we do the send ---*/
	iVertex = 0;
	for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++)
		if (geometry->PeriodicPoint[iPeriodic][0].size() != 0)
			for (iIndex = 0; iIndex < geometry->PeriodicPoint[iPeriodic][0].size(); iIndex++) {
				bound[iMarkerSend][iVertex] = new CVertexMPI(geometry->PeriodicPoint[iPeriodic][0][iIndex], 3);
				bound[iMarkerSend][iVertex]->SetRotation_Type(iPeriodic);
				iVertex++;
			}
    
	/*--- Second we do the receive ---*/
	iVertex = 0;
	for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++)
		if (geometry->PeriodicPoint[iPeriodic][1].size() != 0)
			for (iIndex = 0; iIndex < geometry->PeriodicPoint[iPeriodic][1].size(); iIndex++) {
				bound[iMarkerReceive][iVertex] = new CVertexMPI(geometry->PeriodicPoint[iPeriodic][1][iIndex], 3);
				bound[iMarkerReceive][iVertex]->SetRotation_Type(iPeriodic);
				iVertex++;
			}
    
}

void CPeriodicGeometry::SetMeshFile(CGeometry *geometry, CConfig *config, string val_mesh_out_filename) {
	unsigned long iElem, iPoint, iElem_Bound, GhostPoints;
	unsigned short iMarker, iNodes, iDim;
	unsigned short iMarkerReceive, iPeriodic, nPeriodic = 0;
	ofstream output_file;
	string Grid_Marker;
	char *cstr;
	double *center, *angles, *transl;
    
	cstr = new char [val_mesh_out_filename.size()+1];
	strcpy (cstr, val_mesh_out_filename.c_str());
    
	/*--- Open .su2 grid file ---*/
	output_file.precision(15);
	output_file.open(cstr, ios::out);
    
    /*--- Ghost points, look at the nodes in the send receive ---*/
	iMarkerReceive = nMarker - 1;
	GhostPoints = nElem_Bound[iMarkerReceive];
    
    /*--- Change the numbering to guarantee that the all the receive 
     points are at the end of the file ---*/
    unsigned long OldnPoint = geometry->GetnPoint();
    unsigned long NewSort[nPoint];
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
        NewSort[iPoint] = iPoint;
    }
    
    unsigned long Index = OldnPoint-1;    
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
		if (bound[iMarker][0]->GetVTK_Type() == VERTEX) {
			if (config->GetMarker_All_SendRecv(iMarker) < 0) {
                for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
                    if (bound[iMarker][iElem_Bound]->GetNode(0) < geometry->GetnPoint()) {
                        NewSort[bound[iMarker][iElem_Bound]->GetNode(0)] = Index;
                        NewSort[Index] = bound[iMarker][iElem_Bound]->GetNode(0);
                        Index--;
                    }
                }
            }
		}
	}
    
    
	/*--- Write dimension, number of elements and number of points ---*/
	output_file << "NDIME= " << nDim << endl;
	output_file << "NELEM= " << nElem << endl;
	for (iElem = 0; iElem < nElem; iElem++) {
		output_file << elem[iElem]->GetVTK_Type();
		for (iNodes = 0; iNodes < elem[iElem]->GetnNodes(); iNodes++)
			output_file << "\t" << NewSort[elem[iElem]->GetNode(iNodes)];
		output_file << "\t"<<iElem<<endl;	
	}
    
	output_file << "NPOIN= " << nPoint << "\t" << nPoint - GhostPoints << endl;
	for (iPoint = 0; iPoint < nPoint; iPoint++) {
		for (iDim = 0; iDim < nDim; iDim++)
			output_file << scientific << "\t" << node[NewSort[iPoint]]->GetCoord(iDim) ;
		output_file << "\t" << iPoint <<endl;
	}
    
	output_file << "NMARK= " << nMarker << endl;
	for (iMarker = 0; iMarker < nMarker; iMarker++) {
		if (bound[iMarker][0]->GetVTK_Type() != VERTEX) {
            
			Grid_Marker = config->GetMarker_All_Tag(iMarker);
			output_file << "MARKER_TAG= " << Grid_Marker <<endl;
			output_file << "MARKER_ELEMS= " << nElem_Bound[iMarker] + nNewElem_BoundPer[iMarker] << endl;
            
			for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
				output_file << bound[iMarker][iElem_Bound]->GetVTK_Type() << "\t" ;
				for (iNodes = 0; iNodes < bound[iMarker][iElem_Bound]->GetnNodes()-1; iNodes++)
					output_file << NewSort[bound[iMarker][iElem_Bound]->GetNode(iNodes)] << "\t" ;
				iNodes = bound[iMarker][iElem_Bound]->GetnNodes()-1;
				output_file << NewSort[bound[iMarker][iElem_Bound]->GetNode(iNodes)] << endl;
			}
            
			/*--- Write any new elements at the end of the list. ---*/
			if (nNewElem_BoundPer[iMarker] > 0) {
				for (iElem_Bound = 0; iElem_Bound < nNewElem_BoundPer[iMarker]; iElem_Bound++) {
					output_file << newBoundPer[iMarker][iElem_Bound]->GetVTK_Type() << "\t" ;
					for (iNodes = 0; iNodes < newBoundPer[iMarker][iElem_Bound]->GetnNodes()-1; iNodes++)
						output_file << NewSort[newBoundPer[iMarker][iElem_Bound]->GetNode(iNodes)] << "\t" ;
					iNodes = newBoundPer[iMarker][iElem_Bound]->GetnNodes()-1;
					output_file << NewSort[newBoundPer[iMarker][iElem_Bound]->GetNode(iNodes)] << endl;
				}
			}
            
		}
        
		if (bound[iMarker][0]->GetVTK_Type() == VERTEX) {
			output_file << "MARKER_TAG= SEND_RECEIVE" << endl;
			output_file << "MARKER_ELEMS= " << nElem_Bound[iMarker]<< endl;
			if (config->GetMarker_All_SendRecv(iMarker) > 0) output_file << "SEND_TO= " << config->GetMarker_All_SendRecv(iMarker) << endl;
			if (config->GetMarker_All_SendRecv(iMarker) < 0) output_file << "SEND_TO= " << config->GetMarker_All_SendRecv(iMarker) << endl;
            
			for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
				output_file << bound[iMarker][iElem_Bound]->GetVTK_Type() << "\t" << 
                NewSort[bound[iMarker][iElem_Bound]->GetNode(0)] << "\t" <<
                bound[iMarker][iElem_Bound]->GetRotation_Type() << "\t" <<
                bound[iMarker][iElem_Bound]->GetMatching_Zone() << endl;
			}
		}
	}
    
	/*--- Compute the number of periodic bc on the geometry ---*/
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if (config->GetMarker_All_Boundary(iMarker) == PERIODIC_BOUNDARY)
			nPeriodic++;
    
	output_file << "NPERIODIC= " << nPeriodic + 1 << endl;
    
	/*--- Periodic 0 correspond with no movement of the surface ---*/
	output_file << "PERIODIC_INDEX= 0" << endl;
	output_file << "0.000000000000000e+00" << "\t" << "0.000000000000000e+00" << "\t" << "0.000000000000000e+00" << endl;
	output_file << "0.000000000000000e+00" << "\t" << "0.000000000000000e+00" << "\t" << "0.000000000000000e+00" << endl;
	output_file << "0.000000000000000e+00" << "\t" << "0.000000000000000e+00" << "\t" << "0.000000000000000e+00" << endl;
    
	/*--- From iPeriodic obtain the iMarker ---*/
	for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++) {
		for (iMarker = 0; iMarker < nMarker; iMarker++)
			if (iPeriodic == config->GetMarker_All_PerBound(iMarker)) break;
        
		/*--- Retrieve the supplied periodic information. ---*/
		center = config->GetPeriodicRotCenter(config->GetMarker_All_Tag(iMarker));
		angles = config->GetPeriodicRotAngles(config->GetMarker_All_Tag(iMarker));
		transl = config->GetPeriodicTranslation(config->GetMarker_All_Tag(iMarker));
        
		output_file << "PERIODIC_INDEX= " << iPeriodic << endl;
		output_file << center[0] << "\t" << center[1] << "\t" << center[2] << endl;
		output_file << angles[0] << "\t" << angles[1] << "\t" << angles[2] << endl;
		output_file << transl[0] << "\t" << transl[1] << "\t" << transl[2] << endl;
        
	}
    
    
	output_file.close();
    
}

void CPeriodicGeometry::SetTecPlot(char mesh_filename[200]) {
	unsigned long iElem, iPoint;
	unsigned short iDim;
	ofstream Tecplot_File;
    
	Tecplot_File.open(mesh_filename, ios::out);
	Tecplot_File << "TITLE = \"Visualization of the volumetric grid\"" << endl;
    
	if (nDim == 2) {
		Tecplot_File << "VARIABLES = \"x\",\"y\" " << endl;
		Tecplot_File << "ZONE NODES= "<< nPoint <<", ELEMENTS= "<< nElem <<", DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL"<< endl;
	}
	if (nDim == 3) {
		Tecplot_File << "VARIABLES = \"x\",\"y\",\"z\" " << endl;	
		Tecplot_File << "ZONE NODES= "<< nPoint <<", ELEMENTS= "<< nElem <<", DATAPACKING=POINT, ZONETYPE=FEBRICK"<< endl;
	}
    
	for(iPoint = 0; iPoint < nPoint; iPoint++) {
		for(iDim = 0; iDim < nDim; iDim++)
			Tecplot_File << scientific << node[iPoint]->GetCoord(iDim) << "\t";
		Tecplot_File << "\n";
	}	 
    
	for(iElem = 0; iElem < nElem; iElem++) {
		if (elem[iElem]->GetVTK_Type() == TRIANGLE) {
			Tecplot_File <<
            elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
            elem[iElem]->GetNode(2)+1 <<" "<< elem[iElem]->GetNode(2)+1 << endl;
		}
		if (elem[iElem]->GetVTK_Type() == RECTANGLE) {
			Tecplot_File <<
            elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
            elem[iElem]->GetNode(2)+1 <<" "<< elem[iElem]->GetNode(3)+1 << endl;
		}
		if (elem[iElem]->GetVTK_Type() == TETRAHEDRON) {
			Tecplot_File <<
            elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
            elem[iElem]->GetNode(2)+1 <<" "<< elem[iElem]->GetNode(2)+1 <<" "<<
            elem[iElem]->GetNode(3)+1 <<" "<< elem[iElem]->GetNode(3)+1 <<" "<<
            elem[iElem]->GetNode(3)+1 <<" "<< elem[iElem]->GetNode(3)+1 << endl;
		}
		if (elem[iElem]->GetVTK_Type() == HEXAHEDRON) {
			Tecplot_File <<
            elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
            elem[iElem]->GetNode(2)+1 <<" "<< elem[iElem]->GetNode(3)+1 <<" "<<
            elem[iElem]->GetNode(4)+1 <<" "<< elem[iElem]->GetNode(5)+1 <<" "<<
            elem[iElem]->GetNode(6)+1 <<" "<< elem[iElem]->GetNode(7)+1 << endl;
		}
		if (elem[iElem]->GetVTK_Type() == PYRAMID) {
			Tecplot_File <<
            elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
            elem[iElem]->GetNode(2)+1 <<" "<< elem[iElem]->GetNode(3)+1 <<" "<<
            elem[iElem]->GetNode(4)+1 <<" "<< elem[iElem]->GetNode(4)+1 <<" "<<
            elem[iElem]->GetNode(4)+1 <<" "<< elem[iElem]->GetNode(4)+1 << endl;
		}
		if (elem[iElem]->GetVTK_Type() == WEDGE) {
			Tecplot_File <<
            elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
            elem[iElem]->GetNode(1)+1 <<" "<< elem[iElem]->GetNode(2)+1 <<" "<<
            elem[iElem]->GetNode(3)+1 <<" "<< elem[iElem]->GetNode(4)+1 <<" "<<
            elem[iElem]->GetNode(4)+1 <<" "<< elem[iElem]->GetNode(5)+1 << endl;
		}
	}
    
	Tecplot_File.close();
}

CMultiGridQueue::CMultiGridQueue(unsigned long val_npoint) {
	unsigned long iPoint;
    
	nPoint = val_npoint;
	Priority = new short[nPoint];
	RightCV = new bool[nPoint];
    
	QueueCV.resize(1); 
    
	/*--- Queue initialization with all the points in the finer grid ---*/
	for (iPoint = 0; iPoint < nPoint; iPoint ++) {
		QueueCV[0].push_back(iPoint);
		Priority[iPoint] = 0;
		RightCV[iPoint] = true;
	}
    
}

CMultiGridQueue::~CMultiGridQueue(void) {
    
	delete[] Priority;
	delete[] RightCV;
    
}

void CMultiGridQueue::AddCV(unsigned long val_new_point, unsigned short val_number_neighbors) {
    
	unsigned short Max_Neighbors = QueueCV.size()-1;
    
	/*--- Basic check ---*/
	if (val_new_point > nPoint) {
		cout << "The index of the CV is greater than the size of the priority list." << endl;
		exit(0);
	}
    
	/*--- Resize the list ---*/
	if (val_number_neighbors > Max_Neighbors)
		QueueCV.resize(val_number_neighbors+1);
    
	/*--- Find the point in the queue ---*/
	bool InQueue = false;
	if (Priority[val_new_point] == val_number_neighbors) InQueue = true;
    
	if (!InQueue) {
		/*--- Add the control volume, and update the priority list ---*/
		QueueCV[val_number_neighbors].push_back(val_new_point);
		Priority[val_new_point] = val_number_neighbors;
	}
    
}

void CMultiGridQueue::RemoveCV(unsigned long val_remove_point) {
	unsigned short iPoint;
	bool check;
    
	/*--- Basic check ---*/
	if (val_remove_point > nPoint) {
		cout << "The index of the CV is greater than the size of the priority list." << endl;
		exit(0);
	}
    
	/*--- Find priority of the Control Volume ---*/
	short Number_Neighbors = Priority[val_remove_point];
	if (Number_Neighbors == -1) {
		cout << "The CV "<< val_remove_point <<" is not in the priority list. (RemoveCV)" << endl;
		exit(0);
	}
    
	/*--- Find the point in the queue ---*/
	vector<unsigned long>::iterator ItQueue = find(QueueCV[Number_Neighbors].begin(), 
                                                   QueueCV[Number_Neighbors].end(),
                                                   val_remove_point);
	if( ItQueue != QueueCV[Number_Neighbors].end() ) QueueCV[Number_Neighbors].erase(ItQueue);
    
	Priority[val_remove_point] = -1;
    
	/*--- Check that the size of the queue is the right one ---*/
	unsigned short Size_QueueCV = 0;
	check = false;
	for (iPoint = 0; iPoint < QueueCV.size(); iPoint ++)
		if (QueueCV[iPoint].size() != 0) { Size_QueueCV = iPoint; check = true;}
    
	/*--- Resize the queue, if check = false, the queue is empty, at least 
	 we need one element in the queue ---*/
	if (check) QueueCV.resize(Size_QueueCV+1);
	else QueueCV.resize(1);
    
}

void CMultiGridQueue::MoveCV(unsigned long val_move_point, short val_number_neighbors) {
	unsigned short Priority;
    
	if (val_number_neighbors < 0) {
		val_number_neighbors = 0;
		RightCV[val_move_point] = false;
	}
	else {
		Priority = val_number_neighbors;
		RightCV[val_move_point] = true;
	}
    
	/*--- Remove the control volume ---*/
	RemoveCV(val_move_point);
    
	/*--- Add a new control volume ---*/
	AddCV(val_move_point, val_number_neighbors);
    
}

void CMultiGridQueue::IncrPriorityCV(unsigned long val_incr_point) {
    
	/*--- Find the priority list ---*/
	short Number_Neighbors = Priority[val_incr_point];
	if (Number_Neighbors == -1) {
		cout << "The CV "<< val_incr_point <<" is not in the priority list. (IncrPriorityCV)" << endl;
		exit(0);
	}	
    
	/*--- Remove the control volume ---*/
	RemoveCV(val_incr_point);
    
	/*--- Increase the priority ---*/
	AddCV(val_incr_point, Number_Neighbors+1);
    
}

void CMultiGridQueue::RedPriorityCV(unsigned long val_red_point) {
    
	/*--- Find the priority list ---*/
	short Number_Neighbors = Priority[val_red_point];
	if (Number_Neighbors == -1) {
		cout << "The CV "<< val_red_point <<" is not in the priority list. (RedPriorityCV)" << endl;
		exit(0);
	}	
    
	if (Number_Neighbors != 0) {
        
		/*--- Remove the control volume ---*/
		RemoveCV(val_red_point);
        
		/*--- Increase the priority ---*/
		AddCV(val_red_point, Number_Neighbors-1);
        
	}
    
}

void CMultiGridQueue::VisualizeQueue(void) {
	unsigned short iPoint;
	unsigned long jPoint;
    
	cout << endl;
	for (iPoint = 0; iPoint < QueueCV.size(); iPoint ++) {
		cout << "Number of neighbors " << iPoint <<": ";
		for (jPoint = 0; jPoint < QueueCV[iPoint].size(); jPoint ++) {
			cout << QueueCV[iPoint][jPoint] << " ";
		}
		cout << endl;
	}
    
}

void CMultiGridQueue::VisualizePriority(void) {
	unsigned long iPoint;
    
	for (iPoint = 0; iPoint < nPoint; iPoint ++)
		cout << "Control Volume: " << iPoint <<" Priority: " << Priority[iPoint] << endl;
    
}

long CMultiGridQueue::NextCV(void) {
	if (QueueCV.size() != 0) return QueueCV[QueueCV.size()-1][0];
	else return -1;
}

bool CMultiGridQueue::EmptyQueue(void) {
	unsigned short iPoint;
    
	/*--- In case there is only the no agglomerated elements, 
	 check if they can be agglomerated or we have already finished ---*/
	bool check = true;
    
	if ( QueueCV.size() == 1 ) {
		for (iPoint = 0; iPoint < QueueCV[0].size(); iPoint ++) {
			if (RightCV[QueueCV[0][iPoint]]) { check = false; break; }
		}
	}
	else {	
		for (iPoint = 1; iPoint < QueueCV.size(); iPoint ++)
			if (QueueCV[iPoint].size() != 0) { check = false; break;}
	}
    
	return check;
}

unsigned long CMultiGridQueue::TotalCV(void) {
	unsigned short iPoint;
	unsigned long TotalCV;
    
	TotalCV = 0;
	for (iPoint = 0; iPoint < QueueCV.size(); iPoint ++)
		if (QueueCV[iPoint].size() != 0) { TotalCV += QueueCV[iPoint].size(); }
    
	return TotalCV;
}

void CMultiGridQueue::Update(unsigned long iPoint, CGeometry *fine_grid) {
	unsigned short iNode;
	unsigned long jPoint;
    
	RemoveCV(iPoint);
	for (iNode = 0; iNode <	fine_grid->node[iPoint]->GetnPoint(); iNode ++) {
		jPoint = fine_grid->node[iPoint]->GetPoint(iNode);
		if (fine_grid->node[jPoint]->GetAgglomerate() == false) 
			IncrPriorityCV(jPoint);
	}
    
}
