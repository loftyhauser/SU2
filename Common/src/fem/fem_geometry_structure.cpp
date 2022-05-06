/*!
 * \file fem_geometry_structure.cpp
 * \brief Functions for creating the primal grid for the FEM solver.
 * \author E. van der Weide
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

#include "../../include/fem/fem_geometry_structure.hpp"
#include "../../include/geometry/primal_grid/CPrimalGridFEM.hpp"
#include "../../include/geometry/primal_grid/CPrimalGridBoundFEM.hpp"
#include "../../include/adt/CADTElemClass.hpp"
#include "../../include/adt/CADTPointsOnlyClass.hpp"

/* Prototypes for Lapack functions, if MKL or LAPACK is used. */
#if defined (HAVE_MKL) || defined(HAVE_LAPACK)
extern "C" void dpotrf_(char *, int*, passivedouble*, int*, int*);
extern "C" void dpotri_(char *, int*, passivedouble*, int*, int*);
#endif

bool CLong3T::operator<(const CLong3T &other) const {
  if(long0 != other.long0) return (long0 < other.long0);
  if(long1 != other.long1) return (long1 < other.long1);
  if(long2 != other.long2) return (long2 < other.long2);

  return false;
}

CReorderElements::CReorderElements(const unsigned long  val_GlobalElemID,
                                   const unsigned short val_TimeLevel,
                                   const bool           val_CommSolution,
                                   const unsigned short val_VTK_Type,
                                   const unsigned short val_nPolySol,
                                   const bool           val_JacConstant) {

  /* Copy the global elment ID, time level and whether or not this element
     must be communicated. */
  globalElemID = val_GlobalElemID;
  timeLevel    = val_TimeLevel;
  commSolution = val_CommSolution;

  /* Create the element type used in this class, which stores information of
     the VTK type, the polynomial degree of the solution and whether or not the
     Jacobian of the transformation is constant. As it is possible that the
     polynomial degree of the solution is zero, this convention is different
     from the convention used in the SU2 grid file. */
  elemType = val_VTK_Type + 100*val_nPolySol;
  if( !val_JacConstant ) elemType += 50;
}

bool CReorderElements::operator< (const CReorderElements &other) const {

  /* Elements with the lowest time level are stored first. */
  if(timeLevel != other.timeLevel) return timeLevel < other.timeLevel;

  /* Next comparison is whether or not the element must communicate its
     solution data to other ranks. Elements which do not need to do this
     are stored first. */
  if(commSolution != other.commSolution) return other.commSolution;

  /* Elements of the same element type must be stored as contiguously as
     possible to allow for the simultaneous treatment of elements in the
     matrix multiplications. */
  if(elemType != other.elemType) return elemType < other.elemType;

  /* The final comparison is based on the global element ID. */
  return globalElemID < other.globalElemID;
}

bool CSortFaces::operator()(const CFaceOfElement &f0,
                            const CFaceOfElement &f1) {

  /*--- Comparison in case both faces are boundary faces. ---*/
  if(f0.faceIndicator >= 0 && f1.faceIndicator >= 0) {

    /* Both faces are boundary faces. The first comparison is the boundary
       marker, which is stored in faceIndicator. */
    if(f0.faceIndicator != f1.faceIndicator) return f0.faceIndicator < f1.faceIndicator;

    /* Both faces belong to the same boundary marker. The second comparison is
       based on the on the local volume ID's of the adjacent elements. As the
       volumes are sorted according to the time levels for time accurate
       local stepping, there is no need to do this check seperately here. */
    unsigned long ind0 = f0.elemID0 < nVolElemTot ? f0.elemID0 : f0.elemID1;
    unsigned long ind1 = f1.elemID0 < nVolElemTot ? f1.elemID0 : f1.elemID1;

    return ind0 < ind1;
  }

  /*--- Comparison in case both faces are internal faces. ---*/
  if(f0.faceIndicator == -1 && f1.faceIndicator == -1) {

    /* Both faces are internal faces. First determine the minimum and maximum
       ID of its adjacent elements.  */
    unsigned long elemIDMin0 = min(f0.elemID0, f0.elemID1);
    unsigned long elemIDMax0 = max(f0.elemID0, f0.elemID1);

    unsigned long elemIDMin1 = min(f1.elemID0, f1.elemID1);
    unsigned long elemIDMax1 = max(f1.elemID0, f1.elemID1);

    /* Determine the situation. */
    if(elemIDMax0 < nVolElemTot && elemIDMax1 < nVolElemTot) {

      /* Both faces are matching internal faces. Determine whether or not these
         faces are local faces, i.e. faces between locally owned elements. */
      const bool face0IsLocal = elemIDMax0 < nVolElemOwned;
      const bool face1IsLocal = elemIDMax1 < nVolElemOwned;

      /* Check if both faces have the same status, i.e. either local or
         not local. */
      if(face0IsLocal == face1IsLocal) {

        /* Both faces are either local or not local. Determine the time level
           of the faces, which is the minimum value of the adjacent volume
           elements. */
        const unsigned short timeLevel0 = min(volElem[elemIDMin0].timeLevel,
                                              volElem[elemIDMax0].timeLevel);
        const unsigned short timeLevel1 = min(volElem[elemIDMin1].timeLevel,
                                              volElem[elemIDMax1].timeLevel);

        /* Internal faces with the same status are first sorted according to
           their time level. Faces with the smallest time level are numbered
           first. Note this is only relevant for time accurate local time
           stepping. */
        if(timeLevel0 != timeLevel1) return timeLevel0 < timeLevel1;

        /* The faces belong to the same time level. They are sorted according
           to their element ID's in order to increase cache performance. */
        if(elemIDMin0 != elemIDMin1) return elemIDMin0 < elemIDMin1;
        return elemIDMax0 < elemIDMax1;
      }
      else {

        /* One face is a local face and the other is not. Make sure that
           the local faces are numbered first. */
        if( face0IsLocal ) return true;
        else               return false;
      }
    }
    else if(elemIDMax0 >= nVolElemTot && elemIDMax1 >= nVolElemTot) {

      /* Both faces are non-matching internal faces. Sort them according to
         their relevant element ID. The time level is not taken into account
         yet, because non-matching faces are not possible at the moment with
         time accurate local time stepping. */
      return elemIDMin0 < elemIDMin1;
    }
    else {

      /* One face is a matching internal face and the other face is a
         non-matching internal face. Make sure that the non-matching face
         is numbered after the matching face. This is accomplished by comparing
         the maximum element ID's. */
      return elemIDMax0 < elemIDMax1;
    }
  }

  /*--- One face is a boundary face and the other face is an internal face.
        Make sure that the boundary face is numbered first. This can be
        accomplished by using the > operator for faceIndicator. ---*/
  return f0.faceIndicator > f1.faceIndicator;
}

bool CSortBoundaryFaces::operator()(const CSurfaceElementFEM &f0,
                                    const CSurfaceElementFEM &f1) {

  /* First sorting criterion is the index of the standard element. The
     boundary faces should be sorted per standard element. Note that the
     time level is not taken into account here, because it is assumed that
     the surface elements to be sorted belong to one time level. */
  if(f0.indStandardElement != f1.indStandardElement)
    return f0.indStandardElement < f1.indStandardElement;

  /* The standard elements are the same. The second criterion is the
     corresponding volume IDs of the surface elements. */
  return f0.volElemID < f1.volElemID;
}

bool CPointFEM::operator< (const CPointFEM &other) const {
  if(periodIndexToDonor != other.periodIndexToDonor)
    return periodIndexToDonor < other.periodIndexToDonor;
  return globalID < other.globalID;
 }

bool CPointFEM::operator==(const CPointFEM &other) const {
 return (globalID           == other.globalID &&
         periodIndexToDonor == other.periodIndexToDonor);
}

void CVolumeElementFEM::GetCornerPointsAllFaces(unsigned short &numFaces,
                                                unsigned short nPointsPerFace[],
                                                unsigned long  faceConn[6][4]) {

  /*--- Get the corner connectivities of the faces, local to the element. ---*/
  CPrimalGridFEM::GetLocalCornerPointsAllFaces(VTK_Type, nPolyGrid, nDOFsGrid,
                                               numFaces, nPointsPerFace, faceConn);

  /*--- Convert the local values of faceConn to global values. ---*/
  for(unsigned short i=0; i<numFaces; ++i) {
    for(unsigned short j=0; j<nPointsPerFace[i]; ++j) {
      unsigned long nn = faceConn[i][j];
      faceConn[i][j] = nodeIDsGrid[nn];
    }
  }
}

bool CInternalFaceElementFEM::operator<(const CInternalFaceElementFEM &other) const {

  /* First comparison is the standard element, such that elements with
     the same same standard elements are grouped together. */
  if(indStandardElement != other.indStandardElement)
    return indStandardElement < other.indStandardElement;

  /* The second comparison is the element ID on side 0. */
  if(elemID0 != other.elemID0) return elemID0 < other.elemID0;

  /* The final comparison is the element ID on side 1. */
  return elemID1 < other.elemID1;
}

void CSurfaceElementFEM::GetCornerPointsFace(unsigned short &nPointsPerFace,
                                             unsigned long  faceConn[]) {

  /*--- Get the corner connectivities of the face, local to the element. ---*/
  CPrimalGridBoundFEM::GetLocalCornerPointsFace(VTK_Type, nPolyGrid, nDOFsGrid,
                                                nPointsPerFace, faceConn);

  /*--- Convert the local values of faceConn to global values. ---*/
  for(unsigned short j=0; j<nPointsPerFace; ++j) {
    unsigned long nn = faceConn[j];
    faceConn[j] = nodeIDsGrid[nn];
  }
}

CMeshFEM::CMeshFEM(CGeometry *geometry, CConfig *config) {

  /*--- Allocate the memory for blasFunctions. ---*/
  blasFunctions = new CBlasStructure;

  /*--- The new FEM mesh class has the same problem dimension/zone. ---*/
  nDim  = geometry->GetnDim();
  nZone = geometry->GetnZone();

  /*--- Determine the number of variables stored per DOF.              ---*/
  /*--- This is specifically taylored towards perfect gas computations ---*/
  /*--- and must be generalized later. As this number is also needed   ---*/
  /*--- in the solver, it is worthwhile to create a separate function  ---*/
  /*--- for this purpose.                                              ---*/

  unsigned short nVar = nDim + 1;

  /*--- Determine a mapping from the global point ID to the local index
        of the points.            ---*/
  map<unsigned long,unsigned long> globalPointIDToLocalInd;
  for(unsigned long i=0; i<geometry->GetnPoint(); ++i)
    globalPointIDToLocalInd[geometry->nodes->GetGlobalIndex(i)] = i;

  /*----------------------------------------------------------------------------*/
  /*--- Step 1: Communicate the elements and the boundary elements to the    ---*/
  /*---         ranks where they will be stored during the computation.      ---*/
  /*----------------------------------------------------------------------------*/

  /*--- Determine the ranks to which I have to send my elements. ---*/
  vector<int> sendToRank(size, 0);

  for(unsigned long i=0; i<geometry->GetnElem(); ++i) {
    sendToRank[geometry->elem[i]->GetColor()] = 1;
  }

  map<int,int> rankToIndCommBuf;
  for(int i=0; i<size; ++i) {
    if( sendToRank[i] ) {
      int ind = (int)rankToIndCommBuf.size();
      rankToIndCommBuf[i] = ind;
    }
  }

  /*--- Definition of the communication buffers, used to send the element data
        to the correct ranks.                ---*/
  int nRankSend = (int)rankToIndCommBuf.size();
  vector<vector<short> >     shortSendBuf(nRankSend,  vector<short>(0));
  vector<vector<long>  >     longSendBuf(nRankSend,   vector<long>(0));
  vector<vector<su2double> > doubleSendBuf(nRankSend, vector<su2double>(0));

  /*--- The first element of longSendBuf will contain the number of elements, which
        are stored in the communication buffers. Initialize this value to 0. ---*/
  for(int i=0; i<nRankSend; ++i) longSendBuf[i].push_back(0);

  /*--- Determine the number of ranks, from which this rank will receive elements. ---*/
  int nRankRecv = nRankSend;

#ifdef HAVE_MPI
  vector<int> sizeRecv(size, 1);

  SU2_MPI::Reduce_scatter(sendToRank.data(), &nRankRecv, sizeRecv.data(),
                          MPI_INT, MPI_SUM, SU2_MPI::GetComm());
#endif

  /*--- Loop over the local elements to fill the communication buffers with element data. ---*/
  for(unsigned long i=0; i<geometry->GetnElem(); ++i) {
    int ind = (int)geometry->elem[i]->GetColor();
    map<int,int>::const_iterator MI = rankToIndCommBuf.find(ind);
    ind = MI->second;

    ++longSendBuf[ind][0];   /* The number of elements in the buffers must be incremented. */

    shortSendBuf[ind].push_back(geometry->elem[i]->GetVTK_Type());
    shortSendBuf[ind].push_back(geometry->elem[i]->GetNPolyGrid());
    shortSendBuf[ind].push_back(geometry->elem[i]->GetNPolySol());
    shortSendBuf[ind].push_back(geometry->elem[i]->GetNDOFsGrid());
    shortSendBuf[ind].push_back(geometry->elem[i]->GetNDOFsSol());
    shortSendBuf[ind].push_back(geometry->elem[i]->GetnFaces());
    shortSendBuf[ind].push_back(geometry->elem[i]->GetTimeLevel());
    shortSendBuf[ind].push_back( (short) geometry->elem[i]->GetJacobianConsideredConstant());

    longSendBuf[ind].push_back(geometry->elem[i]->GetGlobalElemID());
    longSendBuf[ind].push_back(geometry->elem[i]->GetGlobalOffsetDOFsSol());

    for(unsigned short j=0; j<geometry->elem[i]->GetNDOFsGrid(); ++j)
      longSendBuf[ind].push_back(geometry->elem[i]->GetNode(j));

    for(unsigned short j=0; j<geometry->elem[i]->GetnFaces(); ++j)
      longSendBuf[ind].push_back(geometry->elem[i]->GetNeighbor_Elements(j));

    for(unsigned short j=0; j<geometry->elem[i]->GetnFaces(); ++j) {
      shortSendBuf[ind].push_back(geometry->elem[i]->GetPeriodicIndex(j));
      shortSendBuf[ind].push_back( (short) geometry->elem[i]->GetJacobianConstantFace(j));
      shortSendBuf[ind].push_back( (short) geometry->elem[i]->GetOwnerFace(j));
    }

    doubleSendBuf[ind].push_back(geometry->elem[i]->GetLengthScale());
  }

  /*--- Determine for each rank to which I have to send elements the data of
        the corresponding nodes.   ---*/
  for(int i=0; i<nRankSend; ++i) {

    /*--- Determine the vector with node IDs in the connectivity
          of the elements for this rank.   ---*/
    vector<long> nodeIDs;

    unsigned long indL = 3;
    unsigned long indS = 3;
    for(long j=0; j<longSendBuf[i][0]; ++j) {
      short nDOFsGrid = shortSendBuf[i][indS], nFaces = shortSendBuf[i][indS+2];
      indS += 3*nFaces+8;

      for(short k=0; k<nDOFsGrid; ++k, ++indL)
        nodeIDs.push_back(longSendBuf[i][indL]);
      indL += nFaces+2;
    }

    /*--- Sort nodeIDs in increasing order and remove the double entities. ---*/
    sort(nodeIDs.begin(), nodeIDs.end());
    vector<long>::iterator lastNodeID = unique(nodeIDs.begin(), nodeIDs.end());
    nodeIDs.erase(lastNodeID, nodeIDs.end());

    /*--- Add the number of node IDs and the node IDs itself to longSendBuf[i]. ---*/
    longSendBuf[i].push_back(nodeIDs.size());
    longSendBuf[i].insert(longSendBuf[i].end(), nodeIDs.begin(), nodeIDs.end());

    /*--- Copy the coordinates to doubleSendBuf. ---*/
    for(unsigned long j=0; j<nodeIDs.size(); ++j) {
      map<unsigned long,unsigned long>::const_iterator LMI;
      LMI = globalPointIDToLocalInd.find(nodeIDs[j]);

      if(LMI == globalPointIDToLocalInd.end())
        SU2_MPI::Error("Entry not found in map", CURRENT_FUNCTION);

      unsigned long ind = LMI->second;
      for(unsigned short l=0; l<nDim; ++l)
        doubleSendBuf[i].push_back(geometry->nodes->GetCoord(ind, l));
    }
  }

  /*--- Loop over the boundaries to send the boundary data to the appropriate rank. ---*/
  nMarker = geometry->GetnMarker();
  for(unsigned short iMarker=0; iMarker<nMarker; ++iMarker) {

    /* Store the current indices in the longSendBuf, which are used to store the
       number of boundary elements sent to this rank. Initialize this value to 0. */
    vector<long> indLongBuf(nRankSend);
    for(int i=0; i<nRankSend; ++i) {
      indLongBuf[i] = longSendBuf[i].size();
      longSendBuf[i].push_back(0);
    }

    /* Loop over the local boundary elements in geometry for this marker. */
    for(unsigned long i=0; i<geometry->GetnElem_Bound(iMarker); ++i) {

      /* Determine the local ID of the corresponding domain element. */
      unsigned long elemID = geometry->bound[iMarker][i]->GetDomainElement()
                           - geometry->beg_node[rank];

      /* Determine to which rank this boundary element must be sent.
         That is the same as its corresponding domain element.
         Update the corresponding index in longSendBuf. */
      int ind = (int)geometry->elem[elemID]->GetColor();
      map<int,int>::const_iterator MI = rankToIndCommBuf.find(ind);
      ind = MI->second;

      ++longSendBuf[ind][indLongBuf[ind]];

      /* Get the donor information for the wall function treatment. */
      const unsigned short nDonors = geometry->bound[iMarker][i]->GetNDonorsWallFunctions();
      const unsigned long  *donors = geometry->bound[iMarker][i]->GetDonorsWallFunctions();

      /* Store the data for this boundary element in the communication buffers. */
      shortSendBuf[ind].push_back(geometry->bound[iMarker][i]->GetVTK_Type());
      shortSendBuf[ind].push_back(geometry->bound[iMarker][i]->GetNPolyGrid());
      shortSendBuf[ind].push_back(geometry->bound[iMarker][i]->GetNDOFsGrid());
      shortSendBuf[ind].push_back(nDonors);

      longSendBuf[ind].push_back(geometry->bound[iMarker][i]->GetDomainElement());
      longSendBuf[ind].push_back(geometry->bound[iMarker][i]->GetGlobalElemID());

      for(unsigned short j=0; j<geometry->bound[iMarker][i]->GetNDOFsGrid(); ++j)
        longSendBuf[ind].push_back(geometry->bound[iMarker][i]->GetNode(j));

      for(unsigned short j=0; j<nDonors; ++j)
        longSendBuf[ind].push_back(donors[j]);
    }
  }

  /*--- Definition of the communication buffers, used to receive
        the element data from the other correct ranks.        ---*/
  vector<vector<short> >     shortRecvBuf(nRankRecv,  vector<short>(0));
  vector<vector<long>  >     longRecvBuf(nRankRecv,   vector<long>(0));
  vector<vector<su2double> > doubleRecvBuf(nRankRecv, vector<su2double>(0));

  /*--- Communicate the data to the correct ranks. Make a distinction
        between parallel and sequential mode.    ---*/
  map<int,int>::const_iterator MI;

#ifdef HAVE_MPI

  /*--- Parallel mode. Send all the data using non-blocking sends. ---*/
  vector<SU2_MPI::Request> commReqs(3*nRankSend);
  MI = rankToIndCommBuf.begin();

  for(int i=0; i<nRankSend; ++i, ++MI) {

    int dest = MI->first;
    SU2_MPI::Isend(shortSendBuf[i].data(), shortSendBuf[i].size(), MPI_SHORT,
                   dest, dest, SU2_MPI::GetComm(), &commReqs[3*i]);
    SU2_MPI::Isend(longSendBuf[i].data(), longSendBuf[i].size(), MPI_LONG,
                   dest, dest+1, SU2_MPI::GetComm(), &commReqs[3*i+1]);
    SU2_MPI::Isend(doubleSendBuf[i].data(), doubleSendBuf[i].size(), MPI_DOUBLE,
                   dest, dest+2, SU2_MPI::GetComm(), &commReqs[3*i+2]);
  }

  /* Loop over the number of ranks from which I receive data. */
  for(int i=0; i<nRankRecv; ++i) {

    /* Block until a message with shorts arrives from any processor.
       Determine the source and the size of the message.   */
    SU2_MPI::Status status;
    SU2_MPI::Probe(MPI_ANY_SOURCE, rank, SU2_MPI::GetComm(), &status);
    int source = status.MPI_SOURCE;

    int sizeMess;
    SU2_MPI::Get_count(&status, MPI_SHORT, &sizeMess);

    /* Allocate the memory for the short receive buffer and receive the message. */
    shortRecvBuf[i].resize(sizeMess);
    SU2_MPI::Recv(shortRecvBuf[i].data(), sizeMess, MPI_SHORT,
                  source, rank, SU2_MPI::GetComm(), &status);

    /* Block until the corresponding message with longs arrives, determine
       its size, allocate the memory and receive the message. */
    SU2_MPI::Probe(source, rank+1, SU2_MPI::GetComm(), &status);
    SU2_MPI::Get_count(&status, MPI_LONG, &sizeMess);
    longRecvBuf[i].resize(sizeMess);

    SU2_MPI::Recv(longRecvBuf[i].data(), sizeMess, MPI_LONG,
                  source, rank+1, SU2_MPI::GetComm(), &status);

    /* Idem for the message with doubles. */
    SU2_MPI::Probe(source, rank+2, SU2_MPI::GetComm(), &status);
    SU2_MPI::Get_count(&status, MPI_DOUBLE, &sizeMess);
    doubleRecvBuf[i].resize(sizeMess);

    SU2_MPI::Recv(doubleRecvBuf[i].data(), sizeMess, MPI_DOUBLE,
                  source, rank+2, SU2_MPI::GetComm(), &status);
  }

  /* Complete the non-blocking sends. */
  SU2_MPI::Waitall(3*nRankSend, commReqs.data(), MPI_STATUSES_IGNORE);

  /* Wild cards have been used in the communication,
     so synchronize the ranks to avoid problems.    */
  SU2_MPI::Barrier(SU2_MPI::GetComm());

#else

  /*--- Sequential mode. Simply copy the buffers. ---*/
  shortRecvBuf[0]  = shortSendBuf[0];
  longRecvBuf[0]   = longSendBuf[0];
  doubleRecvBuf[0] = doubleSendBuf[0];

#endif

  /*--- Release the memory of the send buffers. To make sure that all
        the memory is deleted, the swap function is used. ---*/
  for(int i=0; i<nRankSend; ++i) {
    vector<short>().swap(shortSendBuf[i]);
    vector<long>().swap(longSendBuf[i]);
    vector<su2double>().swap(doubleSendBuf[i]);
  }

  /*--- Allocate the memory for the number of elements for every boundary
        marker and initialize them to zero.     ---*/
  nElem_Bound = new unsigned long[nMarker];
  for(unsigned short i=0; i<nMarker; ++i)
    nElem_Bound[i] = 0;

  /*--- Determine the global element ID's of the elements stored on this rank.
        Sort them in increasing order, such that an easy search can be done.
        In the same loop determine the upper bound for the local nodes (without
        halos) and the number of boundary elements for every marker. ---*/
  nElem = nPoint = 0;
  for(int i=0; i<nRankRecv; ++i) nElem += longRecvBuf[i][0];

  vector<unsigned long> globalElemID;
  globalElemID.reserve(nElem);

  for(int i=0; i<nRankRecv; ++i) {
    unsigned long indL = 1, indS = 0;
    for(long j=0; j<longRecvBuf[i][0]; ++j) {
      globalElemID.push_back(longRecvBuf[i][indL]);

      unsigned short nDOFsGrid = shortRecvBuf[i][indS+3];
      unsigned short nFaces    = shortRecvBuf[i][indS+5];
      indS += 3*nFaces + 8;
      indL += nDOFsGrid + nFaces + 2;
    }

    long nNodesThisRank = longRecvBuf[i][indL];
    nPoint += nNodesThisRank;
    indL   += nNodesThisRank+1;

    for(unsigned iMarker=0; iMarker<nMarker; ++iMarker) {
      long nBoundElemThisRank = longRecvBuf[i][indL]; ++indL;
      nElem_Bound[iMarker] += nBoundElemThisRank;

      for(long j=0; j<nBoundElemThisRank; ++j) {
        short nDOFsBoundElem      = shortRecvBuf[i][indS+2];
        short nDonorsWallFunction = shortRecvBuf[i][indS+3];
        indS += 4;
        indL += nDOFsBoundElem + nDonorsWallFunction + 2;
      }
    }
  }

  sort(globalElemID.begin(), globalElemID.end());

  /*--- Determine the global element ID's of the halo elements. A vector of
        CUnsignedLong2T is used for this purpose, such that a possible periodic
        transformation can be taken into account. Neighbors with a periodic
        transformation will always become a halo element, even if the element
        is stored on this rank. Furthermore a vector of CReorderElements
        is created for the owned elements to be able to reorder the elements,
        see step 2 below. ---*/
  vector<CUnsignedLong2T>  haloElements;
  vector<CReorderElements> ownedElements;
  unsigned short maxTimeLevelLoc = 0;

  for(int i=0; i<nRankRecv; ++i) {
    unsigned long indL = 1, indS = 0;
    for(long j=0; j<longRecvBuf[i][0]; ++j) {
      unsigned long  globalID    = longRecvBuf[i][indL];
      unsigned short VTK_Type    = shortRecvBuf[i][indS];
      unsigned short nPolySol    = shortRecvBuf[i][indS+2];
      unsigned short nDOFsGrid   = shortRecvBuf[i][indS+3];
      unsigned short nFaces      = shortRecvBuf[i][indS+5];
      unsigned short timeLevel   = shortRecvBuf[i][indS+6];
      bool           JacConstant = (bool) shortRecvBuf[i][indS+7];
      bool           commSol     = false;

      /* If the integration rule for constant and non-constant Jacobian
         elements is the same, simply set JacConstant to false, such that
         the elements with constant and non-constant Jacobians are
         considered the same. */
      if( JacConstant ) {
        const unsigned short orderExactStraight =
          (unsigned short) ceil(nPolySol*config->GetQuadrature_Factor_Straight());
        const unsigned short orderExactCurved =
          (unsigned short) ceil(nPolySol*config->GetQuadrature_Factor_Curved());
        if(orderExactStraight == orderExactCurved) JacConstant = false;
      }

      /* Update the local value of the maximum time level. */
      maxTimeLevelLoc = max(maxTimeLevelLoc, timeLevel);

      /*--- Loop over the faces of this element to determine the halo elements
            and the information needed to reorder the owned elements. ---*/
      indS += 8;
      indL += nDOFsGrid + 2;
      for(unsigned short k=0; k<nFaces; ++k, indS+=3, ++indL) {

        if(longRecvBuf[i][indL] != -1) {  /* -1 indicates a boundary face. */

          /* Check if the neighbor of the face is also an owned element.
             Per definition an element for which a periodic transformation
             is carried out, is a halo element, even if the parent element
             is stored locally. */
          bool neighborIsInternal = false;
          if(shortRecvBuf[i][indS] == -1) { /* -1 indicates no periodic transformation. */
            neighborIsInternal = binary_search(globalElemID.begin(),
                                               globalElemID.end(),
                                               longRecvBuf[i][indL]);
          }

          /* Check if this neighbor is not internal and if the element owns the face. */
          if( !neighborIsInternal ) {
            if( shortRecvBuf[i][indS+2] ) {

              /* The face is owned by this element. As the neighboring element
                 is not owned, this implies that a halo element must be created. */
              haloElements.push_back(CUnsignedLong2T(longRecvBuf[i][indL],
                                                     shortRecvBuf[i][indS]+1));  /* The +1, because haloElements */
            }                                                                    /* are unsigned longs.          */
            else {

              /* The face is not owned by this element and therefore it is owned
                 by the neighboring element on a different rank. Consequently the
                 solution of this element must be communicated. */
              commSol = true;
            }
          }
        }
      }

      /* Store the required data for this element in ownedElements. */
      ownedElements.push_back(CReorderElements(globalID, timeLevel, commSol,
                                               VTK_Type, nPolySol, JacConstant));
    }

    /* Skip the part with the node numbers. */
    long nNodesThisRank = longRecvBuf[i][indL];
    indL += nNodesThisRank+1;

    /* Loop over the boundary markers. */
    for(unsigned iMarker=0; iMarker<nMarker; ++iMarker) {
      long nBoundElemThisRank = longRecvBuf[i][indL]; ++indL;

      /* Loop over the boundary elements coming from this rank. */
      for(long j=0; j<nBoundElemThisRank; ++j) {
        short nDOFsBoundElem       = shortRecvBuf[i][indS+2];
        short nDonorsWallFunctions = shortRecvBuf[i][indS+3];
        indS += 4;
        indL += nDOFsBoundElem + 2;

        /* Loop over the donors for the wall functions. */
        for(short k=0; k<nDonorsWallFunctions; ++k, ++indL) {

          /* Check for an external donor. If external, store it in
             haloElements with no periodic transformation. */
          if( !binary_search(globalElemID.begin(), globalElemID.end(),
                             longRecvBuf[i][indL]) )
            haloElements.push_back(CUnsignedLong2T(longRecvBuf[i][indL], 0));
        }
      }
    }
  }

  /* Sort the halo elements in increasing order and remove the double entities. */
  sort(haloElements.begin(), haloElements.end());
  vector<CUnsignedLong2T>::iterator lastHalo = unique(haloElements.begin(), haloElements.end());
  haloElements.erase(lastHalo, haloElements.end());

  /* Determine the maximum global time level and possibly reset the number
     of time levels in config. */
  unsigned short maxTimeLevelGlob = maxTimeLevelLoc;

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&maxTimeLevelLoc, &maxTimeLevelGlob,
                     1, MPI_UNSIGNED_SHORT, MPI_MAX, SU2_MPI::GetComm());
#endif

  const unsigned short nTimeLevels = maxTimeLevelGlob+1;
  config->SetnLevels_TimeAccurateLTS(nTimeLevels);

  /*----------------------------------------------------------------------------*/
  /*--- Step 2: Find out on which rank the halo elements are stored.         ---*/
  /*----------------------------------------------------------------------------*/

  /* Determine the number of owned elements and the total number of
     elements stored on this rank. */
  nVolElemOwned = globalElemID.size();
  nVolElemTot   = nVolElemOwned + haloElements.size();

  /* Determine the map from the global element ID to the current storage
     sequence of ownedElements. */
  map<unsigned long, unsigned long> mapGlobalElemIDToInd;
  for(unsigned long i=0; i<nVolElemOwned; ++i)
    mapGlobalElemIDToInd[ownedElements[i].GetGlobalElemID()] = i;

  /* Determine the number of elements per rank of the originally partitioned grid
     stored in cumulative storage format. */
  vector<unsigned long> nElemPerRankOr(size+1);

  for(int i=0; i<size; ++i) nElemPerRankOr[i] = geometry->beg_node[i];
  nElemPerRankOr[size] = geometry->end_node[size-1];

  /* Determine to which ranks I have to send messages to find out the information
     of the halos stored on this rank. */
  sendToRank.assign(size, 0);

  for(unsigned long i=0; i<haloElements.size(); ++i) {

    /* Determine the rank where this halo element was originally stored. */
    vector<unsigned long>::iterator low;
    low = lower_bound(nElemPerRankOr.begin(), nElemPerRankOr.end(),
                      haloElements[i].long0);
    unsigned long rankHalo = low - nElemPerRankOr.begin();
    if(*low > haloElements[i].long0) --rankHalo;

    sendToRank[rankHalo] = 1;
  }

  rankToIndCommBuf.clear();
  for(int i=0; i<size; ++i) {
    if( sendToRank[i] ) {
      int ind = (int)rankToIndCommBuf.size();
      rankToIndCommBuf[i] = ind;
    }
  }

  /* Resize the first index of the long send buffers for the communication of
     the halo data.        */
  nRankSend = (int)rankToIndCommBuf.size();
  longSendBuf.resize(nRankSend);

  /* Determine the number of ranks, from which this rank will receive elements. */
  nRankRecv = nRankSend;

#ifdef HAVE_MPI
  SU2_MPI::Reduce_scatter(sendToRank.data(), &nRankRecv, sizeRecv.data(),
                          MPI_INT, MPI_SUM, SU2_MPI::GetComm());
#endif

  /* Loop over the local halo elements to fill the communication buffers. */
  for(unsigned long i=0; i<haloElements.size(); ++i) {

    /* Determine the rank where this halo element was originally stored. */
    vector<unsigned long>::iterator low;
    low = lower_bound(nElemPerRankOr.begin(), nElemPerRankOr.end(),
                      haloElements[i].long0);
    unsigned long ind = low - nElemPerRankOr.begin();
    if(*low > haloElements[i].long0) --ind;

    /* Convert this rank to the index in the send buffer. */
    MI = rankToIndCommBuf.find((int)ind);
    ind = MI->second;

    /* Store the global element ID and the periodic index in the long buffer.
       The subtraction of 1 is there to obtain the correct periodic index.
       In haloElements a +1 is added, because this variable is of unsigned long,
       which cannot handle negative numbers. */
    long perIndex = haloElements[i].long1 -1;

    longSendBuf[ind].push_back(haloElements[i].long0);
    longSendBuf[ind].push_back(perIndex);
  }

  /* Define a second set of long receive buffers, because the information from
     the first set is still needed later on. Also define the vector to
     store the ranks from which the message came. */
  vector<vector<long> > longSecondRecvBuf(nRankRecv, vector<long>(0));
  vector<int> sourceRank(nRankRecv);

  /*--- Communicate the data to the correct ranks. Make a distinction
        between parallel and sequential mode.    ---*/

#ifdef HAVE_MPI

  /* Parallel mode. Send all the data using non-blocking sends. */
  commReqs.resize(nRankSend);
  MI = rankToIndCommBuf.begin();

  for(int i=0; i<nRankSend; ++i, ++MI) {
    int dest = MI->first;
    SU2_MPI::Isend(longSendBuf[i].data(), longSendBuf[i].size(), MPI_LONG,
                   dest, dest, SU2_MPI::GetComm(), &commReqs[i]);
  }

  /* Loop over the number of ranks from which I receive data. */
  for(int i=0; i<nRankRecv; ++i) {

    /* Block until a message with longs arrives from any processor.
       Determine the source and the size of the message and receive it. */
    SU2_MPI::Status status;
    SU2_MPI::Probe(MPI_ANY_SOURCE, rank, SU2_MPI::GetComm(), &status);
    sourceRank[i] = status.MPI_SOURCE;

    int sizeMess;
    SU2_MPI::Get_count(&status, MPI_LONG, &sizeMess);

    longSecondRecvBuf[i].resize(sizeMess);
    SU2_MPI::Recv(longSecondRecvBuf[i].data(), sizeMess, MPI_LONG,
                  sourceRank[i], rank, SU2_MPI::GetComm(), &status);
  }

  /* Complete the non-blocking sends. */
  SU2_MPI::Waitall(nRankSend, commReqs.data(), MPI_STATUSES_IGNORE);

#else

  /*--- Sequential mode. Simply copy the buffer, if present at all. ---*/
  for(int i=0; i<nRankRecv; ++i)
    longSecondRecvBuf[i] = longSendBuf[i];

#endif

  /*--- Release the memory of the send buffers. To make sure that all the memory
        is deleted, the swap function is used. Afterwards resize the first index
        of the send buffers to nRankRecv, because this number of messages must
        be sent back to the sending ranks with halo information. ---*/
  for(int i=0; i<nRankSend; ++i) {
    vector<long>().swap(longSendBuf[i]);
  }

  longSendBuf.resize(nRankRecv);

#ifdef HAVE_MPI
  /* Resize the vector of the communication requests to the number of messages
     to be sent by this rank. Only in parallel node. */
  commReqs.resize(nRankRecv);
#endif

  /*--- Loop over the receive buffers to fill and send the send buffers again. ---*/
  for(int i=0; i<nRankRecv; ++i) {

    /* Determine the number of elements present in longSecondRecvBuf[i] and
       reserve the memory for the send buffer. */
    const long nElemBuf = longSecondRecvBuf[i].size()/2;
    longSendBuf[i].reserve(3*nElemBuf);

    /* Loop over the elements stored in the receive buffer. */
    for(long j=0; j<nElemBuf; ++j) {

      /* Get the global element ID and periodic index from the receive buffer. */
      const long globalID = longSecondRecvBuf[i][2*j];
      const long perInd   = longSecondRecvBuf[i][2*j+1];

      /* Determine the local index of the element in the original partitioning.
         Check if the index is valid. */
      const long localID = globalID - geometry->beg_node[rank];
      if(localID < 0 || localID >= (long) geometry->nPointLinear[rank]) {
        ostringstream message;
        message << localID << " " << geometry->nPointLinear[rank] << endl;
        message << "Invalid local element ID";
        SU2_MPI::Error(message.str(), CURRENT_FUNCTION);
      }

      /* Determine which rank owns this element and store everything in the
         send buffer. */
      longSendBuf[i].push_back(globalID);
      longSendBuf[i].push_back(perInd);
      longSendBuf[i].push_back(geometry->elem[localID]->GetColor());
    }

     /* Release the memory of this receive buffer. */
    vector<long>().swap(longSecondRecvBuf[i]);

    /*--- Send the send buffer back to the calling rank.
          Only in parallel mode of course.     ---*/
#ifdef HAVE_MPI
    int dest = sourceRank[i];
    SU2_MPI::Isend(longSendBuf[i].data(), longSendBuf[i].size(), MPI_LONG,
                   dest, dest+1, SU2_MPI::GetComm(), &commReqs[i]);
#endif
  }

  /*--- Resize the first index of the receive buffers to nRankSend, such that
        the requested halo information can be received.     ---*/
  longSecondRecvBuf.resize(nRankSend);

  /*--- Receive the communication data from the correct ranks. Make a distinction
        between parallel and sequential mode.    ---*/
#ifdef HAVE_MPI

  /* Parallel mode. Loop over the number of ranks from which I receive data
     in the return communication, i.e. nRankSend. */
  for(int i=0; i<nRankSend; ++i) {

    /* Block until a message with longs arrives from any processor.
       Determine the source and the size of the message.   */
    SU2_MPI::Status status;
    SU2_MPI::Probe(MPI_ANY_SOURCE, rank+1, SU2_MPI::GetComm(), &status);
    int source = status.MPI_SOURCE;

    int sizeMess;
    SU2_MPI::Get_count(&status, MPI_LONG, &sizeMess);

    /* Allocate the memory for the long receive buffer and receive the message. */
    longSecondRecvBuf[i].resize(sizeMess);
    SU2_MPI::Recv(longSecondRecvBuf[i].data(), sizeMess, MPI_LONG,
                  source, rank+1, SU2_MPI::GetComm(), &status);
  }

  /* Complete the non-blocking sends and synchronize the ranks, because
     wild cards have been used. */
  SU2_MPI::Waitall(nRankRecv, commReqs.data(), MPI_STATUSES_IGNORE);
  SU2_MPI::Barrier(SU2_MPI::GetComm());

#else

  /*--- Sequential mode. Simply copy the buffer, if present at all. ---*/
  for(int i=0; i<nRankRecv; ++i)
    longSecondRecvBuf[i] = longSendBuf[i];

#endif

  /* Release the memory of the send buffers. To make sure that all
     the memory is deleted, the swap function is used. */
  for(int i=0; i<nRankRecv; ++i)
    vector<long>().swap(longSendBuf[i]);

  /* Copy the data from the receive buffers into a class of CLong3T, such that
     it can be sorted in increasing order. Note that the rank of the element
     is stored first, followed by its global ID and last the periodic index. */
  vector<CLong3T> haloData;
  for(int i=0; i<nRankSend; ++i) {
    const long nElemBuf = longSecondRecvBuf[i].size()/3;

    for(long j=0; j<nElemBuf; ++j) {
      const long j3 = 3*j;
      haloData.push_back(CLong3T(longSecondRecvBuf[i][j3+2], longSecondRecvBuf[i][j3],
                                 longSecondRecvBuf[i][j3+1]));
    }

    /* Release the memory of this receive buffer. */
    vector<long>().swap(longSecondRecvBuf[i]);
  }

  /* Sort halo data in increasing order. */
  sort(haloData.begin(), haloData.end());

  /* Determine the number of halo elements per rank in cumulative storage.
     The first element of this vector is nVolElemOwned, such that this vector
     contains the starting position in the vector volElem. */
  vector<unsigned long> nHaloElemPerRank(size+1, 0);
  for(unsigned long i=0; i<haloData.size(); ++i)
    ++nHaloElemPerRank[haloData[i].long0+1];

  nHaloElemPerRank[0] = nVolElemOwned;
  for(int i=0; i<size; ++i)
    nHaloElemPerRank[i+1] += nHaloElemPerRank[i];

  if(nHaloElemPerRank[size] != nVolElemTot)
    SU2_MPI::Error("Inconsistency in total number of volume elements",
                   CURRENT_FUNCTION);

  /* Determine the number of ranks to which I have to send data in this cycle. */
  sendToRank.assign(size, 0);
  rankToIndCommBuf.clear();
  for(int i=0; i<size; ++i) {
    if(nHaloElemPerRank[i+1] > nHaloElemPerRank[i]) {
      sendToRank[i] = 1;
      int ind = (int)rankToIndCommBuf.size();
      rankToIndCommBuf[i] = ind;
    }
  }

  nRankSend = (int)rankToIndCommBuf.size();

  /* Store the value of nRankSend for later use. */
  const int nRankSendHaloInfo = nRankSend;

  /* Determine the number of ranks, from which this rank will receive elements. */
  nRankRecv = nRankSend;

#ifdef HAVE_MPI
  SU2_MPI::Reduce_scatter(sendToRank.data(), &nRankRecv, sizeRecv.data(),
                          MPI_INT, MPI_SUM, SU2_MPI::GetComm());
#endif

  /* Copy the data to be sent to the send buffers. */
  longSendBuf.resize(nRankSend);
  MI = rankToIndCommBuf.begin();

  for(int i=0; i<nRankSend; ++i, ++MI) {
    int dest = MI->first;
    for(unsigned long j=nHaloElemPerRank[dest]; j<nHaloElemPerRank[dest+1]; ++j) {
      const unsigned long jj = j - nVolElemOwned;
      longSendBuf[i].push_back(haloData[jj].long1);
      longSendBuf[i].push_back(haloData[jj].long2);
    }
  }

  /* Resize the first index of the long receive buffer. */
  longSecondRecvBuf.resize(nRankRecv);

  /*--- Communicate the data to the correct ranks. Make a distinction
        between parallel and sequential mode.    ---*/
#ifdef HAVE_MPI

  /* Parallel mode. Send all the data using non-blocking sends. */
  commReqs.resize(nRankSend);
  MI = rankToIndCommBuf.begin();

  for(int i=0; i<nRankSend; ++i, ++MI) {
    int dest = MI->first;
    SU2_MPI::Isend(longSendBuf[i].data(), longSendBuf[i].size(), MPI_LONG,
                   dest, dest, SU2_MPI::GetComm(), &commReqs[i]);
  }

  /* Resize the vector to store the ranks from which the message came. */
  sourceRank.resize(nRankRecv);

  /* Loop over the number of ranks from which I receive data. */
  for(int i=0; i<nRankRecv; ++i) {

    /* Block until a message with longs arrives from any processor.
       Determine the source and the size of the message and receive it. */
    SU2_MPI::Status status;
    SU2_MPI::Probe(MPI_ANY_SOURCE, rank, SU2_MPI::GetComm(), &status);
    sourceRank[i] = status.MPI_SOURCE;

    int sizeMess;
    SU2_MPI::Get_count(&status, MPI_LONG, &sizeMess);

    longSecondRecvBuf[i].resize(sizeMess);
    SU2_MPI::Recv(longSecondRecvBuf[i].data(), sizeMess, MPI_LONG,
                  sourceRank[i], rank, SU2_MPI::GetComm(), &status);
  }

  /* Complete the non-blocking sends and synchronize the ranks,
     because wild cards have been used. */
  SU2_MPI::Waitall(nRankSend, commReqs.data(), MPI_STATUSES_IGNORE);
  SU2_MPI::Barrier(SU2_MPI::GetComm());

#else

  /*--- Sequential mode. Simply copy the buffer, if present at all. ---*/
  for(int i=0; i<nRankSend; ++i)
    longSecondRecvBuf[i] = longSendBuf[i];

#endif

  /* Release the memory of the send buffers. To make sure that all the memory
     is deleted, the swap function is used. */
  for(int i=0; i<nRankSend; ++i) {
    vector<long>().swap(longSendBuf[i]);
  }

  /*--- Loop over the receive buffers to flag the locally owned elements for
        communication. Although the face information has already been used to
        do this when ownedElements are constructed, elements that are donors
        for the wall function treatment and are not direct neighbors may
        have been missed. ---*/
  for(int i=0; i<nRankRecv; ++i) {

    const unsigned long nElemBuf = longSecondRecvBuf[i].size()/2;
    for(unsigned long j=0; j<nElemBuf; ++j) {
      const unsigned long elemID = longSecondRecvBuf[i][2*j];

      map<unsigned long, unsigned long>::iterator MMI = mapGlobalElemIDToInd.find(elemID);
      if(MMI == mapGlobalElemIDToInd.end())
        SU2_MPI::Error("Entry not found in mapGlobalElemIDToInd", CURRENT_FUNCTION);

      ownedElements[MMI->second].SetCommSolution(true);
    }
  }

  /*----------------------------------------------------------------------------*/
  /*--- Step 3: Determine the numbering of the owned elements. The following ---*/
  /*---         criteria are used for the owned elements.                    ---*/
  /*---         - Time level of the element: elements with the smallest time ---*/
  /*---           level are number first, etc.                               ---*/
  /*---         - For each time level the elements that do not need to send  ---*/
  /*---           their data are numbered first, followed by the elements    ---*/
  /*---           that must send their data to other ranks. Note that not    ---*/
  /*---           sending the solution does not mean that the residual can   ---*/
  /*---           be built without communication. It is possible that a face ---*/
  /*---           is owned by a local element, but it is adjacent to an      ---*/
  /*---           element owned by a different rank. In that case the data   ---*/
  /*---           from the neighboring element is communicated and stored in ---*/
  /*---           a halo element. However, the residual of these internal    ---*/
  /*---           elements do not receive a contribution computed on a       ---*/
  /*---           different rank.                                            ---*/
  /*---         - Some of the elements that do not need to send their data   ---*/
  /*---           to other ranks may still be put in that part to increase   ---*/
  /*---           the efficiency when treating element simultaneously, see   ---*/
  /*---           next point.                                                ---*/
  /*---         - Elements of the same type are put contiguously in memory   ---*/
  /*---           as much as possible to facilitate the simultaneous         ---*/
  /*---           treatment in the matrix multiplications to increase        ---*/
  /*---           performance.                                               ---*/
  /*---         - A reverse Cuthill McKee renumbering takes place to obtain  ---*/
  /*---           better cache performance for the face residuals.           ---*/
  /*----------------------------------------------------------------------------*/

  /* Determine the number of elements/faces that are treated simultaneously */
  /* in the matrix products to obtain good gemm performance.                */
  const unsigned short nElemSimul = config->GetSizeMatMulPadding()/nVar;

  /* Determine the number of different element types present. */
  map<unsigned short, unsigned short> mapElemTypeToInd;
  for(vector<CReorderElements>::iterator OEI =ownedElements.begin();
                                         OEI!=ownedElements.end(); ++OEI) {
    const unsigned short elType = OEI->GetElemType();

    if(mapElemTypeToInd.find(elType) == mapElemTypeToInd.end()) {
      const unsigned short ind = mapElemTypeToInd.size();
      mapElemTypeToInd[elType] = ind;
    }
  }

  /* Determine the number of elements per element type per time level.
     Make a distinction between internal elements (elements that are not
     communicated) and elements that are communicated. Later on these vectors
     will be put in cumulative storage format, which explains the +1 for
     the second index. */
  vector<vector<unsigned long> > nInternalElem(nTimeLevels,
                                               vector<unsigned long>(mapElemTypeToInd.size()+1));
  vector<vector<unsigned long> > nCommElem(nTimeLevels,
                                           vector<unsigned long>(mapElemTypeToInd.size()+1));
  for(unsigned short i=0; i<nTimeLevels; ++i)
  {
    for(unsigned long j=0; j<=mapElemTypeToInd.size(); ++j) {
      nInternalElem[i][j] = 0;
      nCommElem[i][j] = 0;
    }
  }

  for(vector<CReorderElements>::iterator OEI =ownedElements.begin();
                                         OEI!=ownedElements.end(); ++OEI) {
    const unsigned short elType = OEI->GetElemType();
    map<unsigned short, unsigned short>::const_iterator MI = mapElemTypeToInd.find(elType);
    unsigned short ind = MI->second +1;

    if( OEI->GetCommSolution() )
      ++nCommElem[OEI->GetTimeLevel()][ind];
    else
      ++nInternalElem[OEI->GetTimeLevel()][ind];
  }

  /* Loop again over the owned elements and check if elements, which do not
     have to send their solution, should be flagged as such in order to improve
     the gemm performance. */
  for(vector<CReorderElements>::iterator OEI =ownedElements.begin();
                                         OEI!=ownedElements.end(); ++OEI) {

    /* Check for an internal element, i.e. an element for which the solution
       does not need to be communicated. */
    if( !OEI->GetCommSolution() ) {

      /* Determine the time level and the index of the element type. */
      const unsigned short tLev = OEI->GetTimeLevel();

      const unsigned short elType = OEI->GetElemType();
      map<unsigned short, unsigned short>::const_iterator MI = mapElemTypeToInd.find(elType);
      const unsigned short ind = MI->second +1;

      /* Determine whether or not this element must be moved from internal to
         comm elements to improve performance. Change the appropriate data
         when this must happen. */
      if(nInternalElem[tLev][ind]%nElemSimul && nCommElem[tLev][ind]%nElemSimul) {
        OEI->SetCommSolution(true);
        --nInternalElem[tLev][ind];
        ++nCommElem[tLev][ind];
      }
    }
  }

  /* Check whether there are enough elements in every partition for every
     time level to guarantee a good gemm performance. */
  unsigned long nFullChunks = 0, nPartialChunks = 0;
  for(unsigned short tLev=0; tLev<nTimeLevels; ++tLev) {
    for(unsigned long j=1; j<nInternalElem[tLev].size(); ++j) {
      nFullChunks += nInternalElem[tLev][j]/nElemSimul;
      if(nInternalElem[tLev][j]%nElemSimul) ++nPartialChunks;
    }

    for(unsigned long j=1; j<nCommElem[tLev].size(); ++j) {
      nFullChunks += nCommElem[tLev][j]/nElemSimul;
      if(nCommElem[tLev][j]%nElemSimul) ++nPartialChunks;
    }
  }

  unsigned long tooManyPartChunksLoc = 0;
  if(5*nPartialChunks >= (nPartialChunks+nFullChunks)) tooManyPartChunksLoc = 1;

  /* Determine the number of ranks which contain too many partial chunks.
     The result only needs to be known on the master node. */
  unsigned long nRanksTooManyPartChunks = tooManyPartChunksLoc;
#ifdef HAVE_MPI
  SU2_MPI::Reduce(&tooManyPartChunksLoc, &nRanksTooManyPartChunks, 1,
                  MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, SU2_MPI::GetComm());
#endif

  if((rank == MASTER_NODE) && (nRanksTooManyPartChunks != 0) && (size > 1)) {
    cout << endl << "                 WARNING" << endl;
    cout << "There are " << nRanksTooManyPartChunks << " partitions for which "
         << " the simultaneous treatment of volume elements is not optimal." << endl;
    cout << "This could lead to a significant load imbalance." << endl;
  }

  /* Put nInternalElem and nCommElem in cumulative storage format. */
  for(unsigned short tLev=0; tLev<nTimeLevels; ++tLev) {
    if( tLev ) nInternalElem[tLev][0] = nCommElem[tLev-1].back();
    for(unsigned long j=1; j<nInternalElem[tLev].size(); ++j)
      nInternalElem[tLev][j] += nInternalElem[tLev][j-1];

    nCommElem[tLev][0] = nInternalElem[tLev].back();
    for(unsigned long j=1; j<nCommElem[tLev].size(); ++j)
      nCommElem[tLev][j] += nCommElem[tLev][j-1];
  }

  /* In the solver only the number of elements per time level is needed.
     The number of owned elements, nVolElemOwnedPerTimeLevel, is stored in
     cumulative storage format and can be retrieved easily from nCommElem.
     The number of internal elements per time level is stored normally and
     can be retrieved from nInternalElem. */
  nVolElemOwnedPerTimeLevel.resize(nTimeLevels+1);
  nVolElemInternalPerTimeLevel.resize(nTimeLevels);

  nVolElemOwnedPerTimeLevel[0] = 0;
  for(unsigned short tLev=0; tLev<nTimeLevels; ++tLev) {
    nVolElemOwnedPerTimeLevel[tLev+1]  = nCommElem[tLev].back();
    nVolElemInternalPerTimeLevel[tLev] = nInternalElem[tLev].back()
                                       - nInternalElem[tLev][0];
  }

  /* Sort the elements of ownedElements in increasing order. */
  sort(ownedElements.begin(), ownedElements.end());

  /* At the moment the elements are sorted per time level and element type.  */
  /* This is in principle enough to make to code work with the simultaneous  */
  /* treatment for volume elements. Below follows a renumbering within the   */
  /* time level and element type (Reverse Cuthill-McKee), such that a better */
  /* cache performance is obtained for the surface elements. For high order  */
  /* elements this effect is probably marginal.                              */

  /* Determine the map from the global element ID to the current storage
     sequence of ownedElements. */
  mapGlobalElemIDToInd.clear();
  for(unsigned long i=0; i<nVolElemOwned; ++i)
    mapGlobalElemIDToInd[ownedElements[i].GetGlobalElemID()] = i;

  /*--- Create the graph of local elements. The halo elements are ignored. ---*/
  vector<vector<unsigned long> > neighElem(nVolElemOwned, vector<unsigned long>(0));

  nRankRecv = (int) longRecvBuf.size();
  for(int i=0; i<nRankRecv; ++i) {
    unsigned long indL = 1, indS = 0;
    for(long j=0; j<longRecvBuf[i][0]; ++j) {
      unsigned long  globalID  = longRecvBuf[i][indL];
      unsigned short nDOFsGrid = shortRecvBuf[i][indS+3];
      unsigned short nFaces    = shortRecvBuf[i][indS+5];

      map<unsigned long, unsigned long>::iterator MMI = mapGlobalElemIDToInd.find(globalID);
      unsigned long ind = MMI->second;

      indS += 8;
      indL += nDOFsGrid + 2;
      for(unsigned short k=0; k<nFaces; ++k, indS+=3, ++indL) {
        if((longRecvBuf[i][indL] != -1) && (shortRecvBuf[i][indS] == -1)) { // Check for internal owned node.

          MMI = mapGlobalElemIDToInd.find(longRecvBuf[i][indL]);
          if(MMI != mapGlobalElemIDToInd.end()) neighElem[ind].push_back(MMI->second);
        }
      }
    }
  }

  /* Sort the neighbors of each element in increasing order. */
  for(unsigned long i=0; i<nVolElemOwned; ++i)
    sort(neighElem[i].begin(), neighElem[i].end());

  /* Define the vector, which contains the new numbering of the owned elements
     w.r.t. to the numbering currently stored in ownedElements. Note that signed
     longs are used for this purpose, because the vector is initialized with -1
     to indicate that no new number has been assigned yet. */
  vector<long> oldElemToNewElem(nVolElemOwned, -1);

  /*--- While loop to carry out the renumbering. A while loop is present,
        because the local partition may not be contiguous. ---*/
  unsigned long nElemRenumbered = 0;
  while (nElemRenumbered < nVolElemOwned) {

    /* Determine the first element in the list that has not been renumbered. */
    unsigned long indBeg;
    for(indBeg=0; indBeg<nVolElemOwned; ++indBeg)
      if(oldElemToNewElem[indBeg] == -1) break;

    /* Determine the time level and element type of the element indBeg and
       determine the element range in which the starting element for the
       renumbering must be sought. */
    unsigned short timeLevel = ownedElements[indBeg].GetTimeLevel();
    unsigned short elType    = ownedElements[indBeg].GetElemType();

    map<unsigned short, unsigned short>::const_iterator MI = mapElemTypeToInd.find(elType);
    unsigned short ind = MI->second;

    unsigned long indEnd;
    if( ownedElements[indBeg].GetCommSolution() )
      indEnd = nCommElem[timeLevel][ind+1];
    else
      indEnd = nInternalElem[timeLevel][ind+1];

    /* Determine the element in the range [indBeg,indEnd) with the least number
       of neighbors that has not been renumbered yet. This is the starting
       element for the current renumbering round. */
    for(unsigned long i=(indBeg+1); i<indEnd; ++i) {
      if((oldElemToNewElem[i] == -1) &&
         (neighElem[i].size() < neighElem[indBeg].size())) indBeg = i;
    }

    /* Start of the Reverse Cuthil McKee renumbering. */
    vector<unsigned long> frontElements(1, indBeg);
    while( frontElements.size() ) {

      /* Vector, which stores the front for the next round. */
      vector<unsigned long> frontElementsNew;

      /* Loop over the elements of the current front. */
      for(unsigned long i=0; i<frontElements.size(); ++i) {

        /* Carry out the renumbering for this element. */
        const unsigned long iFront = frontElements[i];

        timeLevel = ownedElements[iFront].GetTimeLevel();
        elType    = ownedElements[iFront].GetElemType();
        MI  = mapElemTypeToInd.find(elType);
        ind = MI->second;

        if( ownedElements[iFront].GetCommSolution() )
          oldElemToNewElem[iFront] = nCommElem[timeLevel][ind]++;
        else
          oldElemToNewElem[iFront] = nInternalElem[timeLevel][ind]++;

        /* Store the neighbors that have not been renumbered yet in the front
           for the next round. Set its index to -2 to indicate that the element
           is already on the new front. */
        for(unsigned long j=0; j<neighElem[iFront].size(); ++j) {
          if(oldElemToNewElem[neighElem[iFront][j]] == -1) {
            frontElementsNew.push_back(neighElem[iFront][j]);
            oldElemToNewElem[neighElem[iFront][j]] = -2;
          }
        }
      }

      /* Update the counter nElemRenumbered. */
      nElemRenumbered += frontElements.size();

      /* Sort frontElementsNew in increasing order. */
      sort(frontElementsNew.begin(), frontElementsNew.end());

      /* Store the new front elements in frontElements for the next round. */
      frontElements = frontElementsNew;
    }
  }

  if(nElemRenumbered != nVolElemOwned)
    SU2_MPI::Error("Something went wrong in the renumbering", CURRENT_FUNCTION);

  /* Determine the final mapping from the global element number to the local
     entry for the owned elements. First clear mapGlobalElemIDToInd before
     it can be used to store its correct content. */
  mapGlobalElemIDToInd.clear();
  for(unsigned long i=0; i<nVolElemOwned; ++i)
    mapGlobalElemIDToInd[ownedElements[i].GetGlobalElemID()] = oldElemToNewElem[i];

  /*----------------------------------------------------------------------------*/
  /*--- Step 4: Store the elements, nodes and boundary elements in the data  ---*/
  /*---         structures used by the FEM solver.                           ---*/
  /*----------------------------------------------------------------------------*/

  /*--- Check in parallel mode for empty partitions. If present, print a warning.
        The solver is capable of handling empty partitions, but it may not be
        efficient. ---*/
#ifdef HAVE_MPI
  unsigned long thisPartitionEmpty = nVolElemOwned ? 0 : 1;
  unsigned long nEmptyPartitions = 0;

  SU2_MPI::Reduce(&thisPartitionEmpty, &nEmptyPartitions, 1,
                  MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, SU2_MPI::GetComm());

  if(rank == MASTER_NODE && nEmptyPartitions) {
    cout << endl << "         WARNING" << endl;
    cout << "There are " << nEmptyPartitions << " empty partitions present." << endl;
    cout << "SU2 is able to handle this, but it may be inefficient." << endl << endl;
  }
#endif

  /*--- Allocate the memory for the volume elements, the nodes
        and the surface elements of the boundaries.    ---*/
  volElem.resize(nVolElemTot);
  meshPoints.reserve(nPoint);

  boundaries.resize(nMarker);
  for(unsigned short iMarker=0; iMarker<nMarker; ++iMarker) {
    boundaries[iMarker].markerTag        = config->GetMarker_All_TagBound(iMarker);
    boundaries[iMarker].periodicBoundary = config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY;
    boundaries[iMarker].surfElem.reserve(nElem_Bound[iMarker]);
  }

  /*--- Copy the data from the communication buffers. ---*/
  for(int i=0; i<nRankRecv; ++i) {

    /* The data for the volume elements. Loop over these elements in the buffer. */
    unsigned long indL = 1, indS = 0, indD = 0;
    for(long j=0; j<longRecvBuf[i][0]; ++j) {

      /* Determine the location in volElem where this data must be stored. */
      unsigned long elemID = longRecvBuf[i][indL++];
      map<unsigned long, unsigned long>::iterator MMI = mapGlobalElemIDToInd.find(elemID);
      unsigned long ind = MMI->second;

      /* Store the data. */
      volElem[ind].elemIsOwned        = true;
      volElem[ind].rankOriginal       = rank;
      volElem[ind].periodIndexToDonor = -1;

      volElem[ind].VTK_Type  = shortRecvBuf[i][indS++];
      volElem[ind].nPolyGrid = shortRecvBuf[i][indS++];
      volElem[ind].nPolySol  = shortRecvBuf[i][indS++];
      volElem[ind].nDOFsGrid = shortRecvBuf[i][indS++];
      volElem[ind].nDOFsSol  = shortRecvBuf[i][indS++];
      volElem[ind].nFaces    = shortRecvBuf[i][indS++];
      volElem[ind].timeLevel = shortRecvBuf[i][indS++];

      volElem[ind].JacIsConsideredConstant = (bool) shortRecvBuf[i][indS++];

      volElem[ind].elemIDGlobal        = elemID;
      volElem[ind].offsetDOFsSolGlobal = longRecvBuf[i][indL++];

      volElem[ind].nodeIDsGrid.resize(volElem[ind].nDOFsGrid);
      volElem[ind].JacFacesIsConsideredConstant.resize(volElem[ind].nFaces);
      volElem[ind].ElementOwnsFaces.resize(volElem[ind].nFaces);

      for(unsigned short k=0; k<volElem[ind].nDOFsGrid; ++k)
        volElem[ind].nodeIDsGrid[k] = longRecvBuf[i][indL++];

      for(unsigned short k=0; k<volElem[ind].nFaces; ++k) {
        long neighBorID = longRecvBuf[i][indL++];

        ++indS; // At this location the periodic index of the face is stored in
                // shortRecvBuf, which is not stored in volElem.
        volElem[ind].JacFacesIsConsideredConstant[k] = (bool) shortRecvBuf[i][indS++];
        volElem[ind].ElementOwnsFaces[k]             = (bool) shortRecvBuf[i][indS++];

        if(neighBorID == -1)
          volElem[ind].ElementOwnsFaces[k] = true;  // Boundary faces are always owned.
      }

      volElem[ind].lenScale = doubleRecvBuf[i][indD++];
    }

    /* The data for the nodes. Loop over these nodes in the buffer and store
       them in meshPoints. */
    unsigned long nNodesThisRank = longRecvBuf[i][indL++];
    for(unsigned long j=0; j<nNodesThisRank; ++j) {
      CPointFEM thisPoint;
      thisPoint.globalID = longRecvBuf[i][indL++];
      thisPoint.periodIndexToDonor = -1;
      for(unsigned short k=0; k<nDim; ++k)
        thisPoint.coor[k] = doubleRecvBuf[i][indD++];

      meshPoints.push_back(thisPoint);
    }

    /* The data for the boundary markers. Loop over them. */
    for(unsigned short iMarker=0; iMarker<nMarker; ++iMarker) {

      unsigned long nElemThisRank = longRecvBuf[i][indL++];
      for(unsigned long j=0; j<nElemThisRank; ++j) {
        CSurfaceElementFEM thisSurfElem;

        thisSurfElem.VTK_Type  = shortRecvBuf[i][indS++];
        thisSurfElem.nPolyGrid = shortRecvBuf[i][indS++];
        thisSurfElem.nDOFsGrid = shortRecvBuf[i][indS++];
        const short nDonors    = shortRecvBuf[i][indS++];

        thisSurfElem.volElemID         = longRecvBuf[i][indL++];
        thisSurfElem.boundElemIDGlobal = longRecvBuf[i][indL++];

        thisSurfElem.nodeIDsGrid.resize(thisSurfElem.nDOFsGrid);
        for(unsigned short k=0; k<thisSurfElem.nDOFsGrid; ++k)
          thisSurfElem.nodeIDsGrid[k] = longRecvBuf[i][indL++];

        indL += nDonors;

        /* Convert the global volume element ID to the local one.
           It is essential to do this before the sorting. */
        map<unsigned long, unsigned long>::iterator MMI;
        MMI = mapGlobalElemIDToInd.find(thisSurfElem.volElemID);
        thisSurfElem.volElemID = MMI->second;

        /* Store the surface element in the data structure for this boundary. */
        boundaries[iMarker].surfElem.push_back(thisSurfElem);
      }
    }
  }

  /* Sort meshPoints in increasing order and remove the double entities. */
  sort(meshPoints.begin(), meshPoints.end());
  vector<CPointFEM>::iterator lastPoint = unique(meshPoints.begin(), meshPoints.end());
  meshPoints.erase(lastPoint, meshPoints.end());

  /* Clear the contents of the map globalPointIDToLocalInd and fill
     it with the information present in meshPoints. */
  globalPointIDToLocalInd.clear();
  for(unsigned long i=0; i<meshPoints.size(); ++i)
    globalPointIDToLocalInd[meshPoints[i].globalID] = i;

  /*--- All the data from the receive buffers has been copied in the local
        data structures. Release the memory of the receive buffers. To make
        sure that all the memory is deleted, the swap function is used. ---*/
  for(int i=0; i<nRankRecv; ++i) {
    vector<short>().swap(shortRecvBuf[i]);
    vector<long>().swap(longRecvBuf[i]);
    vector<su2double>().swap(doubleRecvBuf[i]);
  }

  /*--- Sort the surface elements of the boundaries in increasing order. ---*/
  for(unsigned short iMarker=0; iMarker<nMarker; ++iMarker)
    sort(boundaries[iMarker].surfElem.begin(), boundaries[iMarker].surfElem.end());

  /*----------------------------------------------------------------------------*/
  /*--- Step 5: Obtain the information of the halo elements, which are       ---*/
  /*---         sorted per time level and afterwards per rank, where the     ---*/
  /*---         sequence is determined by the numbering on the sending rank. ---*/
  /*----------------------------------------------------------------------------*/

  /*--- Resize the first index of the send buffers to nRankRecv, because this
        number of messages must be sent back to the sending ranks with halo
        information. ---*/
  nRankRecv = (int) longSecondRecvBuf.size();
  shortSendBuf.resize(nRankRecv);
  longSendBuf.resize(nRankRecv);
  doubleSendBuf.resize(nRankRecv);

#ifdef HAVE_MPI
  /* Resize the vector of the communication requests to the number of messages
     to be sent by this rank. Only in parallel node. */
  commReqs.resize(3*nRankRecv);
#endif

  /*--- Loop over the receive buffers to fill the send buffers again. ---*/
  for(int i=0; i<nRankRecv; ++i) {

    /* Loop over the elements in this receive buffer to determine the local
       index on this rank. Note that also the periodic index must be stored,
       hence use an CUnsignedLong2T for this purpose. As -1 cannot be stored
       for an unsigned long a 1 is added to the periodic transformation. */
    const unsigned long nElemBuf = longSecondRecvBuf[i].size()/2;
    vector<CUnsignedLong2T> elemBuf(nElemBuf);

    for(unsigned long j=0; j<nElemBuf; ++j) {
      const unsigned long j2 = 2*j;

      const unsigned long elemID = longSecondRecvBuf[i][j2];
      map<unsigned long, unsigned long>::iterator MMI = mapGlobalElemIDToInd.find(elemID);
      if(MMI == mapGlobalElemIDToInd.end())
        SU2_MPI::Error("Entry not found in mapGlobalElemIDToInd", CURRENT_FUNCTION);

      elemBuf[j].long0 = MMI->second;
      elemBuf[j].long1 = longSecondRecvBuf[i][j2+1] + 1;
    }

    /* Release the memory of the long receive buffer via the swap function
       and sort elemBuf in increasing order. */
    vector<long>().swap(longSecondRecvBuf[i]);

    sort(elemBuf.begin(), elemBuf.end());

    /* Store the number of elements in the first element of the long send buffer. */
    longSendBuf[i].push_back(nElemBuf);

    /* Vector with node IDs that must be returned to this calling rank.
       Note that also the periodic index must be stored, hence use an
       CUnsignedLong2T for this purpose. */
    vector<CUnsignedLong2T> nodeIDs;

    /* Loop over the elements to fill the send buffers. */
    for(unsigned long j=0; j<nElemBuf; ++j) {

      /* Store the global element ID in the long buffer,
         the periodic index in the short send buffer and
         the length scale in the double send buffer. */
      const unsigned long indV = elemBuf[j].long0;
      longSendBuf[i].push_back(volElem[indV].elemIDGlobal);

      const short perIndex = (short) elemBuf[j].long1 -1; // Note the -1.
      shortSendBuf[i].push_back(perIndex);

      doubleSendBuf[i].push_back(volElem[indV].lenScale);

      /* Store the other relevant information of this element in the short
         and long communication buffers. Also store the node IDs and the
         periodic transformation in nodeIDs. */
      shortSendBuf[i].push_back(volElem[indV].VTK_Type);
      shortSendBuf[i].push_back(volElem[indV].nPolyGrid);
      shortSendBuf[i].push_back(volElem[indV].nPolySol);
      shortSendBuf[i].push_back(volElem[indV].nDOFsGrid);
      shortSendBuf[i].push_back(volElem[indV].nDOFsSol);
      shortSendBuf[i].push_back(volElem[indV].nFaces);
      shortSendBuf[i].push_back(volElem[indV].timeLevel);

      for(unsigned short k=0; k<volElem[indV].nDOFsGrid; ++k) {
        longSendBuf[i].push_back(volElem[indV].nodeIDsGrid[k]);
        nodeIDs.push_back(CUnsignedLong2T(volElem[indV].nodeIDsGrid[k],
                                          elemBuf[j].long1));
      }

      for(unsigned short k=0; k<volElem[indV].nFaces; ++k)
        shortSendBuf[i].push_back((short) volElem[indV].JacFacesIsConsideredConstant[k]);
    }

    /* Sort nodeIDs in increasing order and remove the double entities. */
    sort(nodeIDs.begin(), nodeIDs.end());
    vector<CUnsignedLong2T>::iterator lastNodeID = unique(nodeIDs.begin(), nodeIDs.end());
    nodeIDs.erase(lastNodeID, nodeIDs.end());

    /* Add the number of node IDs and the node IDs itself to longSendBuf[i]
       and the periodix index to shortSendBuf. Note again the -1 for the
       periodic index, because an unsigned long cannot represent -1, the
       value for the periodic index when no peridicity is present. */
    longSendBuf[i].push_back(nodeIDs.size());
    for(unsigned long j=0; j<nodeIDs.size(); ++j) {
      longSendBuf[i].push_back(nodeIDs[j].long0);
      shortSendBuf[i].push_back( (short) nodeIDs[j].long1-1);
    }

    /*--- Copy the coordinates to doubleSendBuf. ---*/
    for(unsigned long j=0; j<nodeIDs.size(); ++j) {
      map<unsigned long,unsigned long>::const_iterator LMI;
      LMI = globalPointIDToLocalInd.find(nodeIDs[j].long0);

      if(LMI == globalPointIDToLocalInd.end())
        SU2_MPI::Error("Entry not found in map", CURRENT_FUNCTION);

      unsigned long ind = LMI->second;
      for(unsigned short l=0; l<nDim; ++l)
        doubleSendBuf[i].push_back(meshPoints[ind].coor[l]);
    }

    /*--- Send the communication buffers back to the calling rank.
          Only in parallel mode of course.     ---*/
#ifdef HAVE_MPI
    int dest = sourceRank[i];
    SU2_MPI::Isend(shortSendBuf[i].data(), shortSendBuf[i].size(), MPI_SHORT,
                   dest, dest+1, SU2_MPI::GetComm(), &commReqs[3*i]);
    SU2_MPI::Isend(longSendBuf[i].data(), longSendBuf[i].size(), MPI_LONG,
                   dest, dest+2, SU2_MPI::GetComm(), &commReqs[3*i+1]);
    SU2_MPI::Isend(doubleSendBuf[i].data(), doubleSendBuf[i].size(), MPI_DOUBLE,
                   dest, dest+3, SU2_MPI::GetComm(), &commReqs[3*i+2]);
#endif
  }

  /*--- Resize the first index of the receive buffers to nRankSendHaloInfo,
        such that the requested halo information can be received.     ---*/
  nRankSend = nRankSendHaloInfo;
  shortRecvBuf.resize(nRankSend);
  longRecvBuf.resize(nRankSend);
  doubleRecvBuf.resize(nRankSend);

  /* Resize the vector to store the ranks from which the message came. */
  sourceRank.resize(nRankSend);

  /*--- Receive the communication data from the correct ranks. Make a distinction
        between parallel and sequential mode.    ---*/
#ifdef HAVE_MPI

  /* Parallel mode. Loop over the number of ranks from which I receive data
     in the return communication, i.e. nRankSend. */
  for(int i=0; i<nRankSend; ++i) {

    /* Block until a message with shorts arrives from any processor.
       Determine the source and the size of the message.   */
    SU2_MPI::Status status;
    SU2_MPI::Probe(MPI_ANY_SOURCE, rank+1, SU2_MPI::GetComm(), &status);
    sourceRank[i] = status.MPI_SOURCE;

    int sizeMess;
    SU2_MPI::Get_count(&status, MPI_SHORT, &sizeMess);

    /* Allocate the memory for the short receive buffer and receive the message. */
    shortRecvBuf[i].resize(sizeMess);
    SU2_MPI::Recv(shortRecvBuf[i].data(), sizeMess, MPI_SHORT,
                  sourceRank[i], rank+1, SU2_MPI::GetComm(), &status);

    /* Block until the corresponding message with longs arrives, determine
       its size, allocate the memory and receive the message. */
    SU2_MPI::Probe(sourceRank[i], rank+2, SU2_MPI::GetComm(), &status);
    SU2_MPI::Get_count(&status, MPI_LONG, &sizeMess);
    longRecvBuf[i].resize(sizeMess);

    SU2_MPI::Recv(longRecvBuf[i].data(), sizeMess, MPI_LONG,
                  sourceRank[i], rank+2, SU2_MPI::GetComm(), &status);

    /* Idem for the message with doubles. */
    SU2_MPI::Probe(sourceRank[i], rank+3, SU2_MPI::GetComm(), &status);
    SU2_MPI::Get_count(&status, MPI_DOUBLE, &sizeMess);
    doubleRecvBuf[i].resize(sizeMess);

    SU2_MPI::Recv(doubleRecvBuf[i].data(), sizeMess, MPI_DOUBLE,
                  sourceRank[i], rank+3, SU2_MPI::GetComm(), &status);
  }

  /* Complete the non-blocking sends. */
  SU2_MPI::Waitall(3*nRankRecv, commReqs.data(), MPI_STATUSES_IGNORE);

  /* Wild cards have been used in the communication,
     so synchronize the ranks to avoid problems.    */
  SU2_MPI::Barrier(SU2_MPI::GetComm());

#else

  /* Sequential mode. Simply copy the buffers. Note that nRankSend is at most 1. */
  for(int i=0; i<nRankSend; ++i) {
    sourceRank[i]    = i;
    shortRecvBuf[i]  = shortSendBuf[i];
    longRecvBuf[i]   = longSendBuf[i];
    doubleRecvBuf[i] = doubleSendBuf[i];
  }

#endif

  /*--- Release the memory of the send buffers. To make sure that all
        the memory is deleted, the swap function is used. ---*/
  for(int i=0; i<nRankRecv; ++i) {
    vector<short>().swap(shortSendBuf[i]);
    vector<long>().swap(longSendBuf[i]);
    vector<su2double>().swap(doubleSendBuf[i]);
  }

  /*----------------------------------------------------------------------------*/
  /*--- Step 6: Build the layer of halo elements from the information in the ---*/
  /*---         receive buffers shortRecvBuf, longRecvBuf and doubleRecvBuf. ---*/
  /*----------------------------------------------------------------------------*/

  /*--- The halo elements must be sorted first based on the time level, followed
        by the rank and finally the index in the receive buffers (which represents
        the sequence on the sending rank. This sorting can be accomplished by
        using a vector of CLong3T. The contents of this vector is build below. ---*/
  vector<CLong3T> haloElemInfo;
  haloElemInfo.reserve(nVolElemTot - nVolElemOwned);

  for(int i=0; i<nRankSend; ++i) {
    unsigned long indS = 0;

    for(long j=0; j<longRecvBuf[i][0]; ++j) {
      const unsigned short nFaces = shortRecvBuf[i][indS+6];
      haloElemInfo.push_back(CLong3T(shortRecvBuf[i][indS+7], sourceRank[i], j));
      indS += nFaces + 8;
    }
  }

  sort(haloElemInfo.begin(), haloElemInfo.end());

  /*--- Loop over the receive buffers to store the information of the
        halo elements and the halo points. ---*/
  vector<CPointFEM> haloPoints;
  for(int i=0; i<nRankSend; ++i) {

    /* Initialization of the indices in the communication buffers. */
    unsigned long indL = 1, indS = 0, indD = 0;

    /* Loop over the halo elements received from this rank. */
    for(long j=0; j<longRecvBuf[i][0]; ++j) {

      /* Create an object of CLong3T and search for its position
         in haloElemInfo to determine the index in volElem, where
         this element is stored. */
      CLong3T thisElem(shortRecvBuf[i][indS+7], sourceRank[i], j);
      vector<CLong3T>::iterator low;
      low = lower_bound(haloElemInfo.begin(), haloElemInfo.end(), thisElem);

      unsigned long indV = low - haloElemInfo.begin();
      indV += nVolElemOwned;

      /* Retrieve the data from the communication buffers. */
      volElem[indV].elemIDGlobal = longRecvBuf[i][indL++];
      volElem[indV].rankOriginal = sourceRank[i];

      volElem[indV].periodIndexToDonor = shortRecvBuf[i][indS++];
      volElem[indV].VTK_Type           = shortRecvBuf[i][indS++];
      volElem[indV].nPolyGrid          = shortRecvBuf[i][indS++];
      volElem[indV].nPolySol           = shortRecvBuf[i][indS++];
      volElem[indV].nDOFsGrid          = shortRecvBuf[i][indS++];
      volElem[indV].nDOFsSol           = shortRecvBuf[i][indS++];
      volElem[indV].nFaces             = shortRecvBuf[i][indS++];
      volElem[indV].timeLevel          = shortRecvBuf[i][indS++];

      volElem[indV].nodeIDsGrid.resize(volElem[indV].nDOFsGrid);
      for(unsigned short k=0; k<volElem[indV].nDOFsGrid; ++k)
        volElem[indV].nodeIDsGrid[k] = longRecvBuf[i][indL++];

      volElem[indV].JacFacesIsConsideredConstant.resize(volElem[indV].nFaces);
      for(unsigned short k=0; k<volElem[indV].nFaces; ++k)
        volElem[indV].JacFacesIsConsideredConstant[k] = (bool) shortRecvBuf[i][indS++];

      /* Give the member variables that are not obtained via communication their
         values. Some of these variables are not used for halo elements. */
      volElem[indV].elemIsOwned             = false;
      volElem[indV].JacIsConsideredConstant = false;
      volElem[indV].offsetDOFsSolGlobal     = ULONG_MAX;

      /* Halo elements do not own a face per definition. */
      volElem[indV].ElementOwnsFaces.assign(volElem[indV].nFaces, false);

      /* Get the length scale from the double receive buffer.*/
      volElem[indV].lenScale = doubleRecvBuf[i][indD++];
    }

    /* Store the information of the points in haloPoints. */
    const long nPointsThisRank = longRecvBuf[i][indL++];
    for(long j=0; j<nPointsThisRank; ++j) {
      CPointFEM thisPoint;
      thisPoint.globalID           = longRecvBuf[i][indL++];
      thisPoint.periodIndexToDonor = shortRecvBuf[i][indS++];
      for(unsigned short l=0; l<nDim; ++l)
        thisPoint.coor[l] = doubleRecvBuf[i][indD++];

      haloPoints.push_back(thisPoint);
    }

    /* The communication buffers from this rank are not needed anymore.
       Delete them using the swap function. */
    vector<short>().swap(shortRecvBuf[i]);
    vector<long>().swap(longRecvBuf[i]);
    vector<su2double>().swap(doubleRecvBuf[i]);
  }

  /* Remove the duplicate entries from haloPoints. */
  sort(haloPoints.begin(), haloPoints.end());
  lastPoint = unique(haloPoints.begin(), haloPoints.end());
  haloPoints.erase(lastPoint, haloPoints.end());

  /* Initialization of some variables to sort the halo points. */
  Global_nPoint = geometry->GetGlobal_nPoint();
  unsigned long InvalidPointID = Global_nPoint + 10;
  short         InvalidPerInd  = SHRT_MAX;

  /*--- Search for the nonperiodic halo points in the local points to see
        if these points are already stored on this rank. If this is the
        case invalidate this halo and decrease the number of halo points.
        Afterwards remove the invalid halos from the vector.       ---*/
  unsigned long nHaloPoints = haloPoints.size();
  for(unsigned long i=0; i<haloPoints.size(); ++i) {
    if(haloPoints[i].periodIndexToDonor != -1) break;  // Test for nonperiodic.

    if( binary_search(meshPoints.begin(), meshPoints.end(), haloPoints[i]) ) {
      haloPoints[i].globalID           = InvalidPointID;
      haloPoints[i].periodIndexToDonor = InvalidPerInd;
      --nHaloPoints;
    }
  }

  sort(haloPoints.begin(), haloPoints.end());
  haloPoints.resize(nHaloPoints);

  /* Increase the capacity of meshPoints, such that the halo points can be
     stored in there as well. Note that in case periodic points are present
     this is an upper bound. Add the non-periodic halo points to meshPoints. */
  meshPoints.reserve(meshPoints.size() + nHaloPoints);

  for(unsigned long i=0; i<haloPoints.size(); ++i) {
    if(haloPoints[i].periodIndexToDonor != -1) break;  // Test for nonperiodic.

    meshPoints.push_back(haloPoints[i]);
  }

  /* Create a map from the global point ID and periodic index to the local
     index in the vector meshPoints. First store the points already present
     in meshPoints. */
  map<CUnsignedLong2T, unsigned long> mapGlobalPointIDToInd;
  for(unsigned long i=0; i<meshPoints.size(); ++i) {
    CUnsignedLong2T globIndAndPer;
    globIndAndPer.long0 = meshPoints[i].globalID;
    globIndAndPer.long1 = meshPoints[i].periodIndexToDonor+1;  // Note the +1 again.

    mapGlobalPointIDToInd[globIndAndPer] = i;
  }

  /*--- Convert the global indices in the boundary connectivities to local ones.
        Note that the volume ID's already contain the local number. ---*/
  for(unsigned short iMarker=0; iMarker<nMarker; ++iMarker) {
    for(unsigned long i=0; i<boundaries[iMarker].surfElem.size(); ++i) {

      /* Convert the global node ID's to local values. Note that for these node
         ID's no periodic transformation can be present. */
      for(unsigned short j=0; j<boundaries[iMarker].surfElem[i].nDOFsGrid; ++j) {
        CUnsignedLong2T searchItem(boundaries[iMarker].surfElem[i].nodeIDsGrid[j], 0);
        map<CUnsignedLong2T, unsigned long>::const_iterator LLMI;
        LLMI = mapGlobalPointIDToInd.find(searchItem);
        boundaries[iMarker].surfElem[i].nodeIDsGrid[j] = LLMI->second;
      }
    }
  }

  /*--- The only halo points that must be added to meshPoints are the periodic
        halo points. It must be checked whether or not the periodic points in
        haloPoints match with points in meshPoints. This is done below. ---*/
  for(unsigned long iLow=0; iLow<haloPoints.size(); ) {

    /* Determine the upper index for this periodic transformation. */
    unsigned long iUpp;
    for(iUpp=iLow+1; iUpp<haloPoints.size(); ++iUpp)
      if(haloPoints[iUpp].periodIndexToDonor != haloPoints[iLow].periodIndexToDonor) break;

    /* Check for a true periodic index. */
    short perIndex = haloPoints[iLow].periodIndexToDonor;
    if(perIndex != -1) {

      /* Easier storage of the surface elements. */
      vector<CSurfaceElementFEM> &surfElem = boundaries[perIndex].surfElem;

      /*--- In the loop below the coordinates of the points of this local
            periodic boundary as well as a matching tolerance are determined.
            A vector of point ID's is also created, which is needed later on
            when it is checked whether or not a matching point is already
            stored in meshPoints. ---*/
      vector<long>          indInPoints(meshPoints.size(), -1);
      vector<unsigned long> IDsPoints;
      vector<su2double>     coordPoints;
      vector<su2double>     tolPoints;

      for(unsigned long j=0; j<surfElem.size(); ++j) {

        /* Determine the tolerance for equal points, which is a small value
           times the length scale of the adjacent volume element. */
        const su2double tolElem = 1.e-2*volElem[surfElem[j].volElemID].lenScale;

        /* Loop over the nodes of this surface grid and update the points on
           this periodic boundary. */
        for(unsigned short k=0; k<surfElem[j].nDOFsGrid; ++k) {
          unsigned long nn = surfElem[j].nodeIDsGrid[k];

          if(indInPoints[nn] == -1) {

            /* Point is not stored yet in pointsBoundary. Do so now. */
            indInPoints[nn] = IDsPoints.size();
            IDsPoints.push_back(nn);
            tolPoints.push_back(tolElem);

            for(unsigned short l=0; l<nDim; ++l)
              coordPoints.push_back(meshPoints[nn].coor[l]);
          }
          else {

            /* Point is already stored. Update the tolerance. */
            nn = indInPoints[nn];
            tolPoints[nn] = min(tolPoints[nn], tolElem);
          }
        }
      }

      /* Create a local ADT of the points on the periodic boundary. */
      CADTPointsOnlyClass periodicADT(nDim, IDsPoints.size(), coordPoints.data(),
                                      IDsPoints.data(), false);

      /* Get the data for the periodic transformation to the donor. */
      auto center = config->GetPeriodicRotCenter(config->GetMarker_All_TagBound(perIndex));
      auto angles = config->GetPeriodicRotAngles(config->GetMarker_All_TagBound(perIndex));
      auto trans  = config->GetPeriodicTranslation(config->GetMarker_All_TagBound(perIndex));

      /*--- Compute the rotation matrix and translation vector for the
            transformation from the donor. This is the transpose of the
            transformation to the donor. ---*/

      /* Store (center-trans) as it is constant and will be added on. */
      su2double translation[] = {center[0] - trans[0],
                                 center[1] - trans[1],
                                 center[2] - trans[2]};

      /* Store angles separately for clarity. Compute sines/cosines. */
      su2double theta = angles[0];
      su2double phi   = angles[1];
      su2double psi   = angles[2];

      su2double cosTheta = cos(theta), cosPhi = cos(phi), cosPsi = cos(psi);
      su2double sinTheta = sin(theta), sinPhi = sin(phi), sinPsi = sin(psi);

      /* Compute the rotation matrix. Note that the implicit
         ordering is rotation about the x-axis, y-axis, then z-axis. */
      su2double rotMatrix[3][3];
      rotMatrix[0][0] =  cosPhi*cosPsi;
      rotMatrix[0][1] =  cosPhi*sinPsi;
      rotMatrix[0][2] = -sinPhi;

      rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
      rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
      rotMatrix[1][2] = sinTheta*cosPhi;

      rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
      rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
      rotMatrix[2][2] = cosTheta*cosPhi;

      /* Loop over the halo points for this periodic transformation. */
      for(unsigned long i=iLow; i<iUpp; ++i) {

        /* Apply the periodic transformation to the coordinates
           stored in this halo point. */
        su2double dx =             haloPoints[i].coor[0] - center[0];
        su2double dy =             haloPoints[i].coor[1] - center[1];
        su2double dz = nDim == 3 ? haloPoints[i].coor[2] - center[2] : su2double(0.0);

        haloPoints[i].coor[0] = rotMatrix[0][0]*dx + rotMatrix[0][1]*dy
                              + rotMatrix[0][2]*dz + translation[0];
        haloPoints[i].coor[1] = rotMatrix[1][0]*dx + rotMatrix[1][1]*dy
                              + rotMatrix[1][2]*dz + translation[1];
        haloPoints[i].coor[2] = rotMatrix[2][0]*dx + rotMatrix[2][1]*dy
                              + rotMatrix[2][2]*dz + translation[2];

        /* Search for the nearest coordinate in the ADT. */
        su2double dist;
        unsigned long pointID;
        int rankID;

        periodicADT.DetermineNearestNode(haloPoints[i].coor, dist,
                                         pointID, rankID);

        /* Check whether the distance is less equal to the tolerance for
           a matching point. */
        const unsigned long nn = indInPoints[pointID];
        if(dist <= tolPoints[nn]) {

          /* The distance to the nearest point is less than the tolerance,
             hence this periodically transformed point is present on the
             boundary. Store it as such in the map mapGlobalPointIDToInd. */
          CUnsignedLong2T globIndAndPer;
          globIndAndPer.long0 = haloPoints[i].globalID;
          globIndAndPer.long1 = haloPoints[i].periodIndexToDonor+1;  // Note the +1 again.

          mapGlobalPointIDToInd[globIndAndPer] = pointID;
        }
        else {

          /* The distance to the nearest point is larger than the tolerance,
             hence this periodically transformed point is not present yet on
             this rank. Store it in the mapping to the local points and
             create it in meshPoints. */
          CUnsignedLong2T globIndAndPer;
          globIndAndPer.long0 = haloPoints[i].globalID;
          globIndAndPer.long1 = haloPoints[i].periodIndexToDonor+1;  // Note the +1 again.

          mapGlobalPointIDToInd[globIndAndPer] = meshPoints.size();
          meshPoints.push_back(haloPoints[i]);
        }
      }
    }

    /* Set iLow to iUpp for the next periodic transformation. */
    iLow = iUpp;
  }

  /*--- Convert the global node numbering in the elements to a local numbering and
        determine the value of factTimeLevel. This is the number of local time steps
        of the element relative to the largest time step of an element in the mesh.
        This value can only differ from 1 when time accurate local time stepping is
        used. ---*/
  for(unsigned long i=0; i<nVolElemTot; ++i) {
    for(unsigned short j=0; j<volElem[i].nDOFsGrid; ++j) {
      CUnsignedLong2T searchItem(volElem[i].nodeIDsGrid[j],
                                 volElem[i].periodIndexToDonor+1); // Again the +1.
      map<CUnsignedLong2T, unsigned long>::const_iterator LLMI;
      LLMI = mapGlobalPointIDToInd.find(searchItem);
      volElem[i].nodeIDsGrid[j] = LLMI->second;
    }

    const unsigned short diffTimeLevel = nTimeLevels - 1 - volElem[i].timeLevel;
    volElem[i].factTimeLevel = pow(2, diffTimeLevel);
  }

  /* Determine the number of halo elements per time level in cumulative
     storage format. */
  nVolElemHaloPerTimeLevel.assign(nTimeLevels+1, 0);
  for(unsigned long i=nVolElemOwned; i<nVolElemTot; ++i)
    ++nVolElemHaloPerTimeLevel[volElem[i].timeLevel+1];

  nVolElemHaloPerTimeLevel[0] = nVolElemOwned;
  for(unsigned short i=0; i<nTimeLevels; ++i)
    nVolElemHaloPerTimeLevel[i+1] += nVolElemHaloPerTimeLevel[i];
}

void CMeshFEM::ComputeGradientsCoordinatesFace(const unsigned short nIntegration,
                                               const unsigned short nDOFs,
                                               const su2double      *matDerBasisInt,
                                               const unsigned long  *DOFs,
                                               su2double            *derivCoor,
                                               CConfig              *config) {

  /* Allocate the memory to store the values of dxdr, dydr, etc. */
  vector<su2double> helpDxdrVec(nIntegration*nDim*nDim);
  su2double *dxdrVec = helpDxdrVec.data();

  /* Determine the gradients of the Cartesian coordinates w.r.t. the
     parametric coordinates. */
  ComputeGradientsCoorWRTParam(nIntegration, nDOFs, matDerBasisInt, DOFs,
                               dxdrVec, config);

  /* Make a distinction between 2D and 3D to compute the derivatives drdx,
     drdy, etc. */
  switch( nDim ) {
    case 2: {
      /* 2D computation. Store the offset between the r and s derivatives. */
      const unsigned short off = 2*nIntegration;

      /* Loop over the integration points. */
      unsigned short ii = 0;
      for(unsigned short j=0; j<nIntegration; ++j) {

        /* Retrieve the values of dxdr, dydr, dxds and dyds from dxdrVec
           in this integration point. */
        const unsigned short jx = 2*j; const unsigned short jy = jx+1;
        const su2double dxdr = dxdrVec[jx],     dydr = dxdrVec[jy];
        const su2double dxds = dxdrVec[jx+off], dyds = dxdrVec[jy+off];

        /* Compute the inverse relations drdx, drdy, dsdx, dsdy. */
        const su2double Jinv = 1.0/(dxdr*dyds - dxds*dydr);

        derivCoor[ii++] =  dyds*Jinv;  // drdx
        derivCoor[ii++] = -dxds*Jinv;  // drdy
        derivCoor[ii++] = -dydr*Jinv;  // dsdx
        derivCoor[ii++] =  dxdr*Jinv;  // dsdy
      }
      break;
    }

    case 3: {
      /* 3D computation. Store the offset between the r and s and r and t derivatives. */
      const unsigned short offS = 3*nIntegration, offT = 6*nIntegration;

      /* Loop over the integration points. */
      unsigned short ii = 0;
      for(unsigned short j=0; j<nIntegration; ++j) {

        /* Retrieve the values of dxdr, dydr, dzdr, dxds, dyds, dzds, dxdt, dydt
           and dzdt from dxdrVec in this integration point. */
        const unsigned short jx = 3*j; const unsigned short jy = jx+1, jz = jx+2;
        const su2double dxdr = dxdrVec[jx],      dydr = dxdrVec[jy],      dzdr = dxdrVec[jz];
        const su2double dxds = dxdrVec[jx+offS], dyds = dxdrVec[jy+offS], dzds = dxdrVec[jz+offS];
        const su2double dxdt = dxdrVec[jx+offT], dydt = dxdrVec[jy+offT], dzdt = dxdrVec[jz+offT];

        /* Compute the inverse relations drdx, drdy, drdz, dsdx, dsdy, dsdz,
           dtdx, dtdy, dtdz. */
        const su2double Jinv = 1.0/(dxdr*(dyds*dzdt - dzds*dydt)
                             -      dxds*(dydr*dzdt - dzdr*dydt)
                             +      dxdt*(dydr*dzds - dzdr*dyds));

        derivCoor[ii++] = (dyds*dzdt - dzds*dydt)*Jinv;  // drdx
        derivCoor[ii++] = (dzds*dxdt - dxds*dzdt)*Jinv;  // drdy
        derivCoor[ii++] = (dxds*dydt - dyds*dxdt)*Jinv;  // drdz

        derivCoor[ii++] = (dzdr*dydt - dydr*dzdt)*Jinv;  // dsdx
        derivCoor[ii++] = (dxdr*dzdt - dzdr*dxdt)*Jinv;  // dsdy
        derivCoor[ii++] = (dydr*dxdt - dxdr*dydt)*Jinv;  // dsdz

        derivCoor[ii++] = (dydr*dzds - dzdr*dyds)*Jinv;  // dtdx
        derivCoor[ii++] = (dzdr*dxds - dxdr*dzds)*Jinv;  // dtdy
        derivCoor[ii++] = (dxdr*dyds - dydr*dxds)*Jinv;  // dtdz
      }

      break;
    }
  }
}

void CMeshFEM::ComputeGradientsCoorWRTParam(const unsigned short nIntegration,
                                            const unsigned short nDOFs,
                                            const su2double      *matDerBasisInt,
                                            const unsigned long  *DOFs,
                                            su2double            *derivCoor,
                                            CConfig              *config) {

  /* Allocate the memory to store the coordinates as right hand side. */
  vector<su2double> vecRHS(nDOFs*nDim);

  /* Loop over the grid DOFs of the element and copy the coordinates in
     vecRHS in row major order. */
  unsigned long ic = 0;
  for(unsigned short j=0; j<nDOFs; ++j) {
    for(unsigned short k=0; k<nDim; ++k, ++ic)
      vecRHS[ic] = meshPoints[DOFs[j]].coor[k];
  }

  /* Carry out the matrix matrix product The last argument is NULL, such
     that this gemm call is ignored in the profiling. Replace by config if
     if should be included. */
  blasFunctions->gemm(nDim*nIntegration, nDim, nDOFs, matDerBasisInt,
                      vecRHS.data(), derivCoor, nullptr);
}

void CMeshFEM::ComputeNormalsFace(const unsigned short nIntegration,
                                  const unsigned short nDOFs,
                                  const su2double      *dr,
                                  const su2double      *ds,
                                  const unsigned long  *DOFs,
                                  su2double            *normals) {

  /* Initialize the counter ii to 0. ii is the index in normals where the
     information is stored. */
  unsigned int ii = 0;

  /* Make a distinction between 2D and 3D. */
  switch( nDim ) {
    case 2: {
      /* 2D computation. Loop over the integration points of the face. */
      for(unsigned short j=0; j<nIntegration; ++j) {

        /*--- Loop over the number of DOFs of the face to compute dxdr
              and dydr. ---*/
        const su2double *drr = &dr[j*nDOFs];
        su2double dxdr = 0.0, dydr = 0.0;
        for(unsigned short k=0; k<nDOFs; ++k) {
          dxdr += drr[k]*meshPoints[DOFs[k]].coor[0];
          dydr += drr[k]*meshPoints[DOFs[k]].coor[1];
        }

        /* Determine the length of the tangential vector (dxdr, dydr), which
           is also the length of the corresponding normal vector. Also compute
           the inverse of the length. Make sure that a division by zero is
           avoided, although this is most likely never active. */
        const su2double lenNorm    = sqrt(dxdr*dxdr + dydr*dydr);
        const su2double invLenNorm = lenNorm < su2double(1.e-35) ? su2double(1.e+35) : 1.0/lenNorm;

        /* Store the corresponding unit normal vector and its length. The
           direction of the normal vector is such that it is outward pointing
           for the element on side 0 of the face. */
        normals[ii++] =  dydr*invLenNorm;
        normals[ii++] = -dxdr*invLenNorm;
        normals[ii++] =  lenNorm;
      }

      break;
    }

    case 3: {
      /* 3D computation. Loop over the integration points of the face. */
      for(unsigned short j=0; j<nIntegration; ++j) {

        /*--- Loop over the number of DOFs of the face to compute dxdr,
              dxds, dydr, dyds, dzdr and dzds. ---*/
        const su2double *drr = &dr[j*nDOFs], *dss = &ds[j*nDOFs];
        su2double dxdr = 0.0, dydr = 0.0, dzdr = 0.0,
                  dxds = 0.0, dyds = 0.0, dzds = 0.0;
        for(unsigned short k=0; k<nDOFs; ++k) {
          dxdr += drr[k]*meshPoints[DOFs[k]].coor[0];
          dydr += drr[k]*meshPoints[DOFs[k]].coor[1];
          dzdr += drr[k]*meshPoints[DOFs[k]].coor[2];

          dxds += dss[k]*meshPoints[DOFs[k]].coor[0];
          dyds += dss[k]*meshPoints[DOFs[k]].coor[1];
          dzds += dss[k]*meshPoints[DOFs[k]].coor[2];
        }

        /* Compute the vector product dxdr X dxds, where x is the coordinate
           vector (x,y,z). Compute the length of this vector, which is an area,
           as well as the inverse.  Make sure that a division by zero is
           avoided, although this is most likely never active. */
        const su2double nx = dydr*dzds - dyds*dzdr;
        const su2double ny = dxds*dzdr - dxdr*dzds;
        const su2double nz = dxdr*dyds - dxds*dydr;

        const su2double lenNorm    = sqrt(nx*nx + ny*ny + nz*nz);
        const su2double invLenNorm = lenNorm < su2double(1.e-35) ? su2double(1.e+35) : 1.0/lenNorm;

        /* Store the components of the unit normal as well as its length in the
           normals. Note that the current direction of the normal is pointing
           into the direction of the element on side 0 of the face. However,
           in the actual computation of the integral over the faces, it is
           assumed that the vector points in the opposite direction.
           Hence the normal vector must be reversed. */
        normals[ii++] = -nx*invLenNorm;
        normals[ii++] = -ny*invLenNorm;
        normals[ii++] = -nz*invLenNorm;
        normals[ii++] =  lenNorm;
      }
      break;
    }
  }
}

void CMeshFEM::MetricTermsBoundaryFaces(CBoundaryFEM *boundary,
                                        CConfig      *config) {

  /*--- Loop over the boundary faces stored on this rank. ---*/
  for(unsigned long i=0; i<boundary->surfElem.size(); ++i) {

    /*--------------------------------------------------------------------------*/
    /*--- Step 1: Allocate the memory for the face metric terms.             ---*/
    /*---         Depending on the case, not all of this memory is needed.   ---*/
    /*---         - Unit normals + area (nDim+1 per integration point)       ---*/
    /*---         - drdx, dsdx, etc. (nDim*nDim per integration point)       ---*/
    /*---         - grid velocities in the integration points.               ---*/
    /*--------------------------------------------------------------------------*/

    /* Determine the corresponding standard face element and get the
       relevant information from it. */
    const unsigned short ind  = boundary->surfElem[i].indStandardElement;
    const unsigned short nInt = standardBoundaryFacesSol[ind].GetNIntegration();

    /*--- Allocate the several metric terms. ---*/
    boundary->surfElem[i].metricNormalsFace.resize(nInt*(nDim+1));

    /* Allocate the memory for the grid velocities and initialize them to the
       default value of zero. */
    boundary->surfElem[i].gridVelocities.assign(nInt*nDim, 0.0);

    /*--------------------------------------------------------------------------*/
    /*--- Step 2: Determine the actual metric data in the integration points ---*/
    /*---         of the faces.                                              ---*/
    /*--------------------------------------------------------------------------*/

    /* Call the function ComputeNormalsFace to compute the unit normals and
       its corresponding area in the integration points. */
    unsigned short nDOFs = standardBoundaryFacesGrid[ind].GetNDOFsFace();
    const su2double *dr  = standardBoundaryFacesGrid[ind].GetDrBasisFaceIntegration();
    const su2double *ds  = standardBoundaryFacesGrid[ind].GetDsBasisFaceIntegration();

    ComputeNormalsFace(nInt, nDOFs, dr, ds, boundary->surfElem[i].DOFsGridFace.data(),
                       boundary->surfElem[i].metricNormalsFace.data());

    /* Compute the derivatives of the parametric coordinates w.r.t. the
       Cartesian coordinates, i.e. drdx, drdy, etc. in the integration points
       of the face, if needed. */
  }
}

void CMeshFEM::SetPositive_ZArea(CConfig *config) {

  /*---------------------------------------------------------------------------*/
  /*--- Step 1: Determine the local contribution to the positive z area.    ---*/
  /*---------------------------------------------------------------------------*/

  /* Loop over the boundary markers. */
  su2double PositiveZArea = 0.0;
  for(unsigned short iMarker=0; iMarker<nMarker; ++iMarker) {

    /* Determine whether or not this boundary contributes. */
    const unsigned short Boundary   = config->GetMarker_All_KindBC(iMarker);
    const unsigned short Monitoring = config->GetMarker_All_Monitoring(iMarker);

    if( ((Boundary == EULER_WALL)              ||
         (Boundary == HEAT_FLUX)               ||
         (Boundary == ISOTHERMAL)              ||
         (Boundary == LOAD_BOUNDARY)           ||
         (Boundary == DISPLACEMENT_BOUNDARY)) && (Monitoring == YES) ) {

      /* Easier storage of the surface elements for this marker. */
      const vector<CSurfaceElementFEM> &surfElem = boundaries[iMarker].surfElem;

      /* Loop over the surface elements. */
      for(unsigned long i=0; i<surfElem.size(); ++i) {

        /* Determine the number of integration points and their weights via
           the corresponding standard element. */
        const unsigned short ind     = surfElem[i].indStandardElement;
        const unsigned short nInt    = standardBoundaryFacesGrid[ind].GetNIntegration();
        const su2double     *weights = standardBoundaryFacesGrid[ind].GetWeightsIntegration();

        /* Loop over the integration points for this element and update PositiveZArea
           if the normal has a negative z-component. In that case it will give a
           positive contribution to PositiveZArea as this must take the normal pointing
           into the element into account. */
        for(unsigned short j=0; j<nInt; ++j) {

          /* Store the normal data for this integration point a bit easier and update
             PositiveZArea, if needed. */
          const su2double *normal = &surfElem[i].metricNormalsFace[j*(nDim+1)];
          if(normal[nDim-1] < 0.0)
            PositiveZArea -= weights[j]*normal[nDim-1]*normal[nDim];
        }
      }
    }
  }

  /*---------------------------------------------------------------------------*/
  /*--- Step 2: Perform an Allreduce such that the global value of          ---*/
  /*---         PositiveZArea is known on all ranks.                        ---*/
  /*---------------------------------------------------------------------------*/

#ifdef HAVE_MPI
  su2double locArea = PositiveZArea;
  SU2_MPI::Allreduce(&locArea, &PositiveZArea, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
#endif

  /*---------------------------------------------------------------------------*/
  /*--- Step 3: Set the reference area, if this was not specified and write ---*/
  /*---         a message about the projected area if I am the master.      ---*/
  /*---------------------------------------------------------------------------*/

  if (config->GetRefArea() == 0.0)
    config->SetRefArea(PositiveZArea);

  if (rank == MASTER_NODE) {
    if (nDim == 2) cout << "Area projection in the y-plane = "<< PositiveZArea << "." << endl;
    else           cout << "Area projection in the z-plane = "<< PositiveZArea << "." << endl;
  }
}
