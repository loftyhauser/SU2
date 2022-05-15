/*!
 * \file mpi_structure.cpp
 * \brief Main subroutines for the mpi structures.
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

#include "mpi_structure.hpp"


/* Initialise the MPI Communicator Rank and Size */
int CBaseMPIWrapper::Rank = 0;
int CBaseMPIWrapper::Size = 1;

/* Set the default MPI Communicator */
CBaseMPIWrapper::Comm CBaseMPIWrapper::currentComm = 0;  // dummy value

void CBaseMPIWrapper::Error(std::string ErrorMsg, std::string FunctionName){
  if (Rank == 0){
    std::cout << std::endl << std::endl;
    std::cout << "Error in \"" << FunctionName << "\": " << std::endl;
    std::cout <<  "-------------------------------------------------------------------------" << std::endl;
    std::cout << ErrorMsg << std::endl;
    std::cout <<  "------------------------------ Error Exit -------------------------------" << std::endl;
    std::cout << std::endl << std::endl;
  }
  Abort(currentComm, 0);
}

void CBaseMPIWrapper::CopyData(const void* sendbuf, void* recvbuf, int size, Datatype datatype, int recvshift, int sendshift) {
  switch (datatype) {
    case MPI_DOUBLE:
      for (int i = 0; i < size; i++) {
        static_cast<su2double*>(recvbuf)[i+recvshift] = static_cast<const su2double*>(sendbuf)[i+sendshift];
      }
      break;
    case MPI_UNSIGNED_LONG:
      for (int i = 0; i < size; i++) {
        static_cast<unsigned long*>(recvbuf)[i+recvshift] = static_cast<const unsigned long*>(sendbuf)[i+sendshift];
      }
      break;
    case MPI_LONG:
      for (int i = 0; i < size; i++) {
        static_cast<long*>(recvbuf)[i+recvshift] = static_cast<const long*>(sendbuf)[i+sendshift];
      }
      break;
    case MPI_UNSIGNED_SHORT:
      for (int i = 0; i < size; i++) {
        static_cast<unsigned short*>(recvbuf)[i+recvshift] = static_cast<const unsigned short*>(sendbuf)[i+sendshift];
      }
      break;
    case MPI_CHAR:
      for (int i = 0; i < size; i++) {
        static_cast<char*>(recvbuf)[i+recvshift] = static_cast<const char*>(sendbuf)[i+sendshift];
      }
      break;
    case MPI_SHORT:
      for (int i = 0; i < size; i++) {
        static_cast<short*>(recvbuf)[i+recvshift] = static_cast<const short*>(sendbuf)[i+sendshift];
      }
      break;
    case MPI_INT:
      for (int i = 0; i < size; i++) {
        static_cast<int*>(recvbuf)[i+recvshift] = static_cast<const int*>(sendbuf)[i+sendshift];
      }
      break;
    default:
      Error("Unknown type", CURRENT_FUNCTION);
      break;
  };
}
