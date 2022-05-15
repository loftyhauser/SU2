/*!
 * \file CParallelFileWriter.cpp
 * \brief Filewriter base class.
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

#include "../../../include/output/filewriter/CFileWriter.hpp"

CFileWriter::CFileWriter(CParallelDataSorter *valDataSorter, string valFileExt):
  fileExt(valFileExt),
  dataSorter(valDataSorter){

  rank = SU2_MPI::GetRank();
  size = SU2_MPI::GetSize();

  fileSize = 0.0;
  bandwidth = 0.0;

}

CFileWriter::CFileWriter(string valFileExt):
  fileExt(valFileExt){

  rank = SU2_MPI::GetRank();
  size = SU2_MPI::GetSize();

  fileSize = 0.0;
  bandwidth = 0.0;

}

CFileWriter::~CFileWriter(){

}

bool CFileWriter::WriteMPIBinaryDataAll(const void *data, unsigned long sizeInBytes,
                                        unsigned long totalSizeInBytes, unsigned long offsetInBytes){


  startTime = SU2_MPI::Wtime();

  unsigned long bytesWritten;

  /*--- Write binary data ---*/

  bytesWritten = fwrite(data, sizeof(char), sizeInBytes, fhw);
  fileSize += bytesWritten;

  stopTime = SU2_MPI::Wtime();

  usedTime += stopTime - startTime;

  return (bytesWritten == sizeInBytes);

}

bool CFileWriter::WriteMPIBinaryData(const void *data, unsigned long sizeInBytes, unsigned short processor){

  startTime = SU2_MPI::Wtime();

  unsigned long bytesWritten = sizeInBytes;

  /*--- Write the total size in bytes at the beginning of the binary data blob ---*/

  bytesWritten = fwrite(data, sizeof(char), sizeInBytes, fhw);

  stopTime = SU2_MPI::Wtime();

  usedTime += stopTime - startTime;

  return (bytesWritten == sizeInBytes);

}

bool CFileWriter::WriteMPIString(const string &str, unsigned short processor){

  startTime = SU2_MPI::Wtime();

  unsigned long bytesWritten;
  bytesWritten = fwrite(str.c_str(), sizeof(char), str.size(), fhw);

  fileSize += bytesWritten;

  stopTime = SU2_MPI::Wtime();

  usedTime += stopTime - startTime;

  return (bytesWritten == str.size()*sizeof(char));

}

bool CFileWriter::OpenMPIFile(string val_filename){

  /*--- We append the pre-defined suffix (extension) to the filename (prefix) ---*/
  val_filename.append(fileExt);

  fhw = fopen(val_filename.c_str(), "wb");
  /*--- Error check for opening the file. ---*/

  if (!fhw) {
    SU2_MPI::Error(string("Unable to open file ") +
                   val_filename, CURRENT_FUNCTION);
  }

  fileSize = 0.0;
  usedTime = 0;

  return true;
}

bool CFileWriter::CloseMPIFile(){

  fclose(fhw);

  /*--- Communicate the total file size for the restart ---*/

  su2double my_fileSize = fileSize;
  SU2_MPI::Allreduce(&my_fileSize, &fileSize, 1,
                     MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());

  /*--- Compute and store the bandwidth ---*/

  bandwidth = fileSize/(1.0e6)/usedTime;

  return true;
}

