/*!
 * \file mpi_structure.hpp
 * \brief Headers of the mpi interface for generalized datatypes.
 *        The subroutines and functions are in the <i>mpi_structure.cpp</i> file.
 * \author T. Albring
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

#ifdef HAVE_MPI
#include <map>
#include "mpi.h"
#endif

#include <stdlib.h>
#include "../basic_types/datatype_structure.hpp"
#ifndef _MSC_VER
#include <unistd.h>
#else
#include <io.h>
#endif

#include "omp_structure.hpp"

/* Depending on the compiler, define the correct macro to get the current function name */

#if defined(__GNUC__) || (defined(__ICC) && (__ICC >= 600))
#define CURRENT_FUNCTION __PRETTY_FUNCTION__
#elif defined(__FUNCSIG__)
#define CURRENT_FUNCTION __FUNCSIG__
#elif (defined(__INTEL_COMPILER) && (__INTEL_COMPILER >= 600))
#define CURRENT_FUNCTION __FUNCTION__
#elif defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 199901)
#define CURRENT_FUNCTION __func__
#elif defined(__cplusplus) && (__cplusplus >= 201103)
#define CURRENT_FUNCTION __func__
#else
#define CURRENT_FUNCTION "(unknown)"
#endif

#ifdef HAVE_MPI

/*--- Depending on the datatype used, the correct MPI wrapper class is defined.
 * For the default (double type) case this results in using the normal MPI routines. ---*/
class CBaseMPIWrapper;
typedef CBaseMPIWrapper SU2_MPI;

/*!
 * \class CMPIWrapper
 * \brief Class for defining the MPI wrapper routines; this class features as a base class for
 * MPI interfaces for non-primitive dataypes e.g. used by AD, complex etc.
 */
class CBaseMPIWrapper {
 public:
  typedef MPI_Request Request;
  typedef MPI_Status Status;
  typedef MPI_Datatype Datatype;
  typedef MPI_Op Op;
  typedef MPI_Comm Comm;
  typedef MPI_Win Win;

 protected:
  static int Rank, Size, MinRankError;
  static Comm currentComm;
  static bool winMinRankErrorInUse;
  static Win winMinRankError;

 public:
  static void Error(std::string ErrorMsg, std::string FunctionName);

  static inline int GetRank() { return Rank; }

  static inline int GetSize() { return Size; }

  static inline void SetComm(Comm newComm) {
    currentComm = newComm;
    MPI_Comm_rank(currentComm, &Rank);
    MPI_Comm_size(currentComm, &Size);

    if (winMinRankErrorInUse) MPI_Win_free(&winMinRankError);
    MinRankError = Size;
    MPI_Win_create(&MinRankError, sizeof(int), sizeof(int), MPI_INFO_NULL, currentComm, &winMinRankError);
    winMinRankErrorInUse = true;
  }

  static inline Comm GetComm() { return currentComm; }

  static inline void Init(int* argc, char*** argv) {
    MPI_Init(argc, argv);
    MPI_Comm_rank(currentComm, &Rank);
    MPI_Comm_size(currentComm, &Size);

    MinRankError = Size;
    MPI_Win_create(&MinRankError, sizeof(int), sizeof(int), MPI_INFO_NULL, currentComm, &winMinRankError);
    winMinRankErrorInUse = true;
  }

  static inline void Init_thread(int* argc, char*** argv, int required, int* provided) {
    MPI_Init_thread(argc, argv, required, provided);
    MPI_Comm_rank(currentComm, &Rank);
    MPI_Comm_size(currentComm, &Size);

    MinRankError = Size;
    MPI_Win_create(&MinRankError, sizeof(int), sizeof(int), MPI_INFO_NULL, currentComm, &winMinRankError);
    winMinRankErrorInUse = true;
  }

  static inline void Comm_rank(Comm comm, int* rank) { MPI_Comm_rank(comm, rank); }

  static inline void Comm_size(Comm comm, int* size) { MPI_Comm_size(comm, size); }

  static inline void Finalize() {
    if (winMinRankErrorInUse) MPI_Win_free(&winMinRankError);
    MPI_Finalize();
  }

  static inline void Barrier(Comm comm) { MPI_Barrier(comm); }

  static inline void Abort(Comm comm, int error) { MPI_Abort(comm, error); }

  static inline void Get_count(const Status* status, Datatype datatype, int* count) {
    MPI_Get_count(status, datatype, count);
  }

  static inline void Isend(const void* buf, int count, Datatype datatype, int dest, int tag, Comm comm,
                           Request* request) {
    MPI_Isend(buf, count, datatype, dest, tag, comm, request);
  }

  static inline void Irecv(void* buf, int count, Datatype datatype, int dest, int tag, Comm comm, Request* request) {
    MPI_Irecv(buf, count, datatype, dest, tag, comm, request);
  }

  static inline void Wait(Request* request, Status* status) { MPI_Wait(request, status); }

  static inline int Request_free(Request *request) { return MPI_Request_free(request); }

  static inline void Testall(int count, Request* array_of_requests, int* flag, Status* array_of_statuses) {
    MPI_Testall(count, array_of_requests, flag, array_of_statuses);
  }

  static inline void Waitall(int nrequests, Request* request, Status* status) {
    MPI_Waitall(nrequests, request, status);
  }

  static inline void Probe(int source, int tag, Comm comm, Status* status) { MPI_Probe(source, tag, comm, status); }

  static inline void Send(const void* buf, int count, Datatype datatype, int dest, int tag, Comm comm) {
    MPI_Send(buf, count, datatype, dest, tag, comm);
  }

  static inline void Recv(void* buf, int count, Datatype datatype, int dest, int tag, Comm comm, Status* status) {
    MPI_Recv(buf, count, datatype, dest, tag, comm, status);
  }

  static inline void Bcast(void* buf, int count, Datatype datatype, int root, Comm comm) {
    MPI_Bcast(buf, count, datatype, root, comm);
  }

  static inline void Reduce(const void* sendbuf, void* recvbuf, int count, Datatype datatype, Op op, int root,
                            Comm comm) {
    MPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm);
  }

  static inline void Allreduce(const void* sendbuf, void* recvbuf, int count, Datatype datatype, Op op, Comm comm) {
    MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
  }

  static inline void Gather(const void* sendbuf, int sendcnt, Datatype sendtype, void* recvbuf, int recvcnt,
                            Datatype recvtype, int root, Comm comm) {
    MPI_Gather(sendbuf, sendcnt, sendtype, recvbuf, recvcnt, recvtype, root, comm);
  }

  static inline void Scatter(const void* sendbuf, int sendcnt, Datatype sendtype, void* recvbuf, int recvcnt,
                             Datatype recvtype, int root, Comm comm) {
    MPI_Scatter(sendbuf, sendcnt, sendtype, recvbuf, recvcnt, recvtype, root, comm);
  }

  static inline void Allgather(const void* sendbuf, int sendcnt, Datatype sendtype, void* recvbuf, int recvcnt,
                               Datatype recvtype, Comm comm) {
    MPI_Allgather(sendbuf, sendcnt, sendtype, recvbuf, recvcnt, recvtype, comm);
  }

  static inline void Allgatherv(const void* sendbuf, int sendcount, Datatype sendtype, void* recvbuf,
                                const int* recvcounts, const int* displs, Datatype recvtype, Comm comm) {
    MPI_Allgatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, comm);
  }

  static inline void Alltoall(const void* sendbuf, int sendcount, Datatype sendtype, void* recvbuf, int recvcount,
                              Datatype recvtype, Comm comm) {
    MPI_Alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);
  }

  static inline void Alltoallv(const void* sendbuf, const int* sendcounts, const int* sdispls, Datatype sendtype,
                               void* recvbuf, const int* recvcounts, const int* recvdispls, Datatype recvtype,
                               Comm comm) {
    MPI_Alltoallv(sendbuf, sendcounts, sdispls, sendtype, recvbuf, recvcounts, recvdispls, recvtype, comm);
  }

  static inline void Sendrecv(const void* sendbuf, int sendcnt, Datatype sendtype, int dest, int sendtag, void* recvbuf,
                              int recvcnt, Datatype recvtype, int source, int recvtag, Comm comm, Status* status) {
    MPI_Sendrecv(sendbuf, sendcnt, sendtype, dest, sendtag, recvbuf, recvcnt, recvtype, source, recvtag, comm, status);
  }

  static inline void Reduce_scatter(const void* sendbuf, void* recvbuf, const int* recvcounts, Datatype datatype, Op op,
                                    Comm comm) {
    MPI_Reduce_scatter(sendbuf, recvbuf, recvcounts, datatype, op, comm);
  }

  static inline void Waitany(int nrequests, Request* request, int* index, Status* status) {
    MPI_Waitany(nrequests, request, index, status);
  }

  static inline passivedouble Wtime(void) { return MPI_Wtime(); }
};

typedef MPI_Comm SU2_Comm;


#else  // HAVE_MPI

#define MPI_UNSIGNED_LONG 1
#define MPI_LONG 2
#define MPI_UNSIGNED_SHORT 3
#define MPI_FLOAT 4
#define MPI_DOUBLE 4
#define MPI_ANY_SOURCE 5
#define MPI_SUM 6
#define MPI_CHAR 7
#define MPI_SHORT 8
#define MPI_MIN 9
#define MPI_MAX 10
#define MPI_INT 11
#define MPI_PROD 12
#define MPI_STATUS_IGNORE nullptr

/*!
 * \class CMPIWrapper
 * \brief Version for when there is no MPI.
 */
class CBaseMPIWrapper {
 public:
  typedef int Comm;
  typedef int Datatype;
  typedef int Request;
  typedef int Op;

  struct Status {
    int MPI_TAG;
    int MPI_SOURCE;
    Status() : MPI_TAG(0), MPI_SOURCE(0) {}
  };

 private:
  static int Rank, Size;
  static Comm currentComm;

  static void CopyData(const void* sendbuf, void* recvbuf, int size, Datatype datatype, int recvshift=0, int sendshift=0);

 public:
  static void Error(std::string ErrorMsg, std::string FunctionName);

  static inline int GetRank() { return Rank; }

  static inline int GetSize() { return Size; }

  static inline void SetComm(Comm newComm) { currentComm = newComm; }

  static inline Comm GetComm() { return currentComm; }

  static inline void Init(int* argc, char*** argv) {}

  static inline void Init_thread(int* argc, char*** argv, int required, int* provided) { *provided = required; }

  static inline void Barrier(Comm comm) {}

  static inline void Abort(Comm comm, int error) { exit(EXIT_FAILURE); }

  static inline void Comm_rank(Comm comm, int* rank) { *rank = 0; }

  static inline void Comm_size(Comm comm, int* size) { *size = 1; }

  static inline void Finalize() {}

  static inline void Isend(const void* buf, int count, Datatype datatype, int dest, int tag, Comm comm,
                           Request* request) {}

  static inline void Irecv(void* buf, int count, Datatype datatype, int source, int tag, Comm comm, Request* request) {}

  static inline void Wait(Request* request, Status* status) {}

  static inline int Request_free(Request *request) { return 0; }

  static inline void Waitall(int nrequests, Request* request, Status* status) {}

  static inline void Waitany(int nrequests, Request* request, int* index, Status* status) {}

  static inline void Send(const void* buf, int count, Datatype datatype, int dest, int tag, Comm comm) {}

  static inline void Recv(void* buf, int count, Datatype datatype, int dest, int tag, Comm comm, Status* status) {}

  static inline void Bcast(void* buf, int count, Datatype datatype, int root, Comm comm) {}

  static inline void Reduce(const void* sendbuf, void* recvbuf, int count, Datatype datatype, Op op, int root,
                            Comm comm) {
    CopyData(sendbuf, recvbuf, count, datatype);
  }

  static inline void Allreduce(const void* sendbuf, void* recvbuf, int count, Datatype datatype, Op op, Comm comm) {
    CopyData(sendbuf, recvbuf, count, datatype);
  }

  static inline void Gather(const void* sendbuf, int sendcnt, Datatype sendtype, void* recvbuf, int recvcnt,
                            Datatype recvtype, int root, Comm comm) {
    CopyData(sendbuf, recvbuf, sendcnt, sendtype);
  }

  static inline void Scatter(const void* sendbuf, int sendcnt, Datatype sendtype, void* recvbuf, int recvcnt,
                             Datatype recvtype, int root, Comm comm) {
    CopyData(sendbuf, recvbuf, sendcnt, sendtype);
  }

  static inline void Allgatherv(const void* sendbuf, int sendcnt, Datatype sendtype, void* recvbuf, const int* recvcnt,
                                const int* displs, Datatype recvtype, Comm comm) {
    CopyData(sendbuf, recvbuf, sendcnt, sendtype, displs[0]);
  }

  static inline void Allgather(const void* sendbuf, int sendcnt, Datatype sendtype, void* recvbuf, int recvcnt,
                               Datatype recvtype, Comm comm) {
    CopyData(sendbuf, recvbuf, sendcnt, sendtype);
  }

  static inline void Sendrecv(const void* sendbuf, int sendcnt, Datatype sendtype, int dest, int sendtag, void* recvbuf,
                              int recvcnt, Datatype recvtype, int source, int recvtag, Comm comm, Status* status) {
    CopyData(sendbuf, recvbuf, sendcnt, sendtype);
  }

  static inline void Reduce_scatter(const void* sendbuf, void* recvbuf, const int* recvcounts, Datatype datatype, Op op,
                                    Comm comm) {
    CopyData(sendbuf, recvbuf, recvcounts[0], datatype);
  }

  static inline void Alltoall(const void* sendbuf, int sendcount, Datatype sendtype, void* recvbuf, int recvcount,
                              Datatype recvtype, Comm comm) {
    CopyData(sendbuf, recvbuf, recvcount, sendtype);
  }

  static inline void Alltoallv(const void* sendbuf, const int* sendcounts, const int* sdispls, Datatype sendtype,
                               void* recvbuf, const int* recvcounts, const int* recvdispls, Datatype recvtype,
                               Comm comm) {
    CopyData(sendbuf, recvbuf, recvcounts[0], recvtype, recvdispls[0], sdispls[0]);
  }

  static inline void Probe(int source, int tag, Comm comm, Status* status) {}

  static inline passivedouble Wtime(void) { return omp_get_wtime(); }
};
typedef int SU2_Comm;
typedef CBaseMPIWrapper SU2_MPI;

#endif

/*--- Select the appropriate MPI wrapper based on datatype, to use in templated classes. ---*/
template <class T>
struct SelectMPIWrapper {
  typedef SU2_MPI W;
};
