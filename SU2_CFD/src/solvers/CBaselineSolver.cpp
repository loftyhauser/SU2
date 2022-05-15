/*!
 * \file CBaselineSolver.cpp
 * \brief Main subroutines for CBaselineSolver class.
 * \author F. Palacios, T. Economon
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


#include "../../include/solvers/CBaselineSolver.hpp"
#include "../../../Common/include/toolboxes/printing_toolbox.hpp"

CBaselineSolver::CBaselineSolver(void) : CSolver() { }

CBaselineSolver::CBaselineSolver(CGeometry *geometry, CConfig *config) {

  nPoint = geometry->GetnPoint();

  /*--- Define geometry constants in the solver structure ---*/

  nDim = geometry->GetnDim();

  /*--- Routines to access the number of variables and string names. ---*/

  SetOutputVariables(geometry, config);

  /*--- Initialize a zero solution and instantiate the CVariable class. ---*/

  Solution = new su2double[nVar];
  for (unsigned short iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;

  nodes = new CBaselineVariable(nPoint, nVar, config);
  SetBaseClassPointerToNodes();

}

CBaselineSolver::CBaselineSolver(CGeometry *geometry, CConfig *config, unsigned short val_nvar, vector<string> field_names) {

  /*--- Define geometry constants in the solver structure ---*/

  nPoint = geometry->GetnPoint();
  nDim = geometry->GetnDim();
  nVar = val_nvar;
  fields = field_names;

  /*--- Allocate the node variables ---*/

  nodes = new CBaselineVariable(nPoint, nVar, config);
  SetBaseClassPointerToNodes();


}

void CBaselineSolver::SetOutputVariables(CGeometry *geometry, CConfig *config) {

  /*--- Open the restart file and extract the nVar and field names. ---*/

  string Tag, text_line;

  ifstream restart_file;
  string filename;

  /*--- Retrieve filename from config ---*/

  if (config->GetContinuous_Adjoint() || config->GetDiscrete_Adjoint()) {
    filename = config->GetSolution_AdjFileName();
    filename = config->GetObjFunc_Extension(filename);
  } else {
    filename = config->GetSolution_FileName();
  }


  /*--- Read only the number of variables in the restart file. ---*/

  if (config->GetRead_Binary_Restart()) {

    /*--- Multizone problems require the number of the zone to be appended. ---*/

    filename = config->GetFilename(filename, ".dat", config->GetTimeIter());

    char fname[100];
    strcpy(fname, filename.c_str());
    int nVar_Buf = 5;
    int var_buf[5];

    /*--- Serial binary input. ---*/

    FILE *fhw;
    fhw = fopen(fname,"rb");
    size_t ret;

    /*--- Error check for opening the file. ---*/

    if (!fhw) {
      SU2_MPI::Error(string("Unable to open SU2 restart file ") + string(fname), CURRENT_FUNCTION);
    }

    /*--- First, read the number of variables and points. ---*/

    ret = fread(var_buf, sizeof(int), nVar_Buf, fhw);
    if (ret != (unsigned long)nVar_Buf) {
      SU2_MPI::Error("Error reading restart file.", CURRENT_FUNCTION);
    }

    /*--- Check that this is an SU2 binary file. SU2 binary files
     have the hex representation of "SU2" as the first int in the file. ---*/

    if (var_buf[0] != 535532) {
      SU2_MPI::Error(string("File ") + string(fname) + string(" is not a binary SU2 restart file.\n") +
                     string("SU2 reads/writes binary restart files by default.\n") +
                     string("Note that backward compatibility for ASCII restart files is\n") +
                     string("possible with the READ_BINARY_RESTART option."), CURRENT_FUNCTION);
    }

    /*--- Close the file. ---*/

    fclose(fhw);

    /*--- Set the number of variables, one per field in the
     restart file (without including the PointID) ---*/

    nVar = var_buf[1];
  } else {

    /*--- Multizone problems require the number of the zone to be appended. ---*/

    filename = config->GetFilename(filename, ".csv", config->GetTimeIter());

    /*--- First, check that this is not a binary restart file. ---*/

    char fname[100];
    strcpy(fname, filename.c_str());
    int magic_number;


    /*--- Serial binary input. ---*/

    FILE *fhw;
    fhw = fopen(fname,"rb");
    size_t ret;

    /*--- Error check for opening the file. ---*/

    if (!fhw) {
      SU2_MPI::Error(string("Unable to open SU2 restart file ") + string(fname), CURRENT_FUNCTION);
    }

    /*--- Attempt to read the first int, which should be our magic number. ---*/

    ret = fread(&magic_number, sizeof(int), 1, fhw);
    if (ret != 1) {
      SU2_MPI::Error("Error reading restart file.", CURRENT_FUNCTION);
    }

    /*--- Check that this is an SU2 binary file. SU2 binary files
     have the hex representation of "SU2" as the first int in the file. ---*/

    if (magic_number == 535532) {
      SU2_MPI::Error(string("File ") + string(fname) + string(" is a binary SU2 restart file, expected ASCII.\n") +
                     string("SU2 reads/writes binary restart files by default.\n") +
                     string("Note that backward compatibility for ASCII restart files is\n") +
                     string("possible with the READ_BINARY_RESTART option."), CURRENT_FUNCTION);
    }

    fclose(fhw);


    /*--- Open the restart file ---*/

    restart_file.open(filename.data(), ios::in);

    /*--- In case there is no restart file ---*/

    if (restart_file.fail()) {
      SU2_MPI::Error(string("SU2 solution file ") + filename + string(" not found"), CURRENT_FUNCTION);
    }

    /*--- Identify the number of fields (and names) in the restart file ---*/

    getline (restart_file, text_line);

    fields = PrintingToolbox::split(text_line, ',');

    for (unsigned short iField = 0; iField < fields.size(); iField++){
      PrintingToolbox::trim(fields[iField]);
    }

    /*--- Close the file (the solution date is read later). ---*/

    restart_file.close();

    /*--- Set the number of variables, one per field in the
     restart file (without including the PointID) ---*/

    nVar = fields.size() - 1;

  }

}

void CBaselineSolver::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo) {

  /*--- Restart the solution from file information ---*/

  string filename;
  unsigned long index;
  unsigned short iDim, iVar;
  bool adjoint = ( config->GetContinuous_Adjoint() || config->GetDiscrete_Adjoint() );
  unsigned short iInst = config->GetiInst();
  bool steady_restart = config->GetSteadyRestart();
  TURB_MODEL turb_model = config->GetKind_Turb_Model();

  su2double *Coord = new su2double [nDim];
  for (iDim = 0; iDim < nDim; iDim++)
    Coord[iDim] = 0.0;

  /*--- Skip coordinates ---*/

  unsigned short skipVars = geometry[iInst]->GetnDim();

  /*--- Retrieve filename from config ---*/

  if (adjoint) {
    filename = config->GetSolution_AdjFileName();
    filename = config->GetObjFunc_Extension(filename);
  } else {
    filename = config->GetSolution_FileName();
  }

  filename = config->GetFilename(filename, "", val_iter);

  /*--- Output the file name to the console. ---*/

  if (rank == MASTER_NODE)
    cout << "Reading and storing the solution from " << filename
    << "." << endl;

  /*--- Read the restart data from either an ASCII or binary SU2 file. ---*/

  if (config->GetRead_Binary_Restart()) {
    Read_SU2_Restart_Binary(geometry[iInst], config, filename);
  } else {
    Read_SU2_Restart_ASCII(geometry[iInst], config, filename);
  }

  int counter = 0;
  long iPoint_Local = 0; unsigned long iPoint_Global = 0;

  /*--- Load data from the restart into correct containers. ---*/

  for (iPoint_Global = 0; iPoint_Global < geometry[iInst]->GetGlobal_nPointDomain(); iPoint_Global++ ) {

    /*--- Retrieve local index. If this node from the restart file lives
     on the current processor, we will load and instantiate the vars. ---*/

    iPoint_Local = geometry[iInst]->GetGlobal_to_Local_Point(iPoint_Global);

    if (iPoint_Local > -1) {

      /*--- We need to store this point's data, so jump to the correct
       offset in the buffer of data from the restart file and load it. ---*/

      index = counter*Restart_Vars[1];
      for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = Restart_Data[index+iVar];
      nodes->SetSolution(iPoint_Local,Solution);


      /*--- Increment the overall counter for how many points have been loaded. ---*/
      counter++;
    }

  }

  /*--- MPI solution ---*/

  InitiateComms(geometry[iInst], config, SOLUTION);
  CompleteComms(geometry[iInst], config, SOLUTION);


  delete [] Coord;

  /*--- Delete the class memory that is used to load the restart. ---*/

  delete [] Restart_Vars;
  delete [] Restart_Data;
  Restart_Vars = nullptr; Restart_Data = nullptr;

}

void CBaselineSolver::LoadRestart_FSI(CGeometry *geometry, CConfig *config, int val_iter) {

  /*--- Restart the solution from file information ---*/
  string filename;
  unsigned long index;
  unsigned short iVar;
  bool adjoint = (config->GetContinuous_Adjoint() || config->GetDiscrete_Adjoint());

  /*--- Retrieve filename from config ---*/
  if (adjoint) {
    filename = config->GetSolution_AdjFileName();
    filename = config->GetObjFunc_Extension(filename);
  } else {
    filename = config->GetSolution_FileName();
  }

  /*--- Multizone problems require the number of the zone to be appended. ---*/

  filename = config->GetFilename(filename, "", val_iter);

  /*--- Output the file name to the console. ---*/

  if (rank == MASTER_NODE)
    cout << "Reading and storing the solution from " << filename
    << "." << endl;

  /*--- Read the restart data from either an ASCII or binary SU2 file. ---*/

  if (config->GetRead_Binary_Restart()) {
    Read_SU2_Restart_Binary(geometry, config, filename);
  } else {
    Read_SU2_Restart_ASCII(geometry, config, filename);
  }

  unsigned short nVar_Local = Restart_Vars[1];
  su2double *Solution_Local = new su2double[nVar_Local];

  int counter = 0;
  long iPoint_Local = 0; unsigned long iPoint_Global = 0;

  /*--- Load data from the restart into correct containers. ---*/

  for (iPoint_Global = 0; iPoint_Global < geometry->GetGlobal_nPointDomain(); iPoint_Global++ ) {

    /*--- Retrieve local index. If this node from the restart file lives
     on the current processor, we will load and instantiate the vars. ---*/

    iPoint_Local = geometry->GetGlobal_to_Local_Point(iPoint_Global);

    if (iPoint_Local > -1) {

      /*--- We need to store this point's data, so jump to the correct
       offset in the buffer of data from the restart file and load it. ---*/

      index = counter*Restart_Vars[1];
      for (iVar = 0; iVar < nVar_Local; iVar++) Solution[iVar] = Restart_Data[index+iVar];
      nodes->SetSolution(iPoint_Local,Solution);

      /*--- Increment the overall counter for how many points have been loaded. ---*/

      counter++;

    }

  }

  delete [] Solution_Local;

}

CBaselineSolver::~CBaselineSolver(void) {
  delete nodes;
}
