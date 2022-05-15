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

/*--- Select the appropriate MPI wrapper based on datatype, to use in templated classes. ---*/
template <class T>
struct SelectMPIWrapper {
  typedef SU2_MPI W;
};
