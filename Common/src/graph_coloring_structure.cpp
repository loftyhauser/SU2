/*!
 * \file graph_coloring_structure.cpp
 * \brief Functions used to carry out the coloring of a given graph.
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

#include "../include/graph_coloring_structure.hpp"

/* Constructor. Nothing to be done. */
CGraphColoringStructure::CGraphColoringStructure(void) {}

/* Destructor. Nothing to be done. */
CGraphColoringStructure::~CGraphColoringStructure(void) {}

/* Function, which determines the colors for the vertices of the given graph. */
void CGraphColoringStructure::GraphVertexColoring(
                         CConfig                              *config,
                         const vector<unsigned long>          &nVerticesPerRank,
                         const vector<vector<unsigned long> > &entriesVertices,
                         int                                  &nGlobalColors,
                         vector<int>                          &colorLocalVertices) {

  /* Determine the number of ranks and the current rank. */
  int nRank  = 1;
  int myRank = 0;

  /*--- Determine the algorithm to use for the graph coloring. ---*/
  switch( config->GetKind_Matrix_Coloring() ) {

    case GREEDY_COLORING: {

      /* Greedy algorithm, which is implemented sequentially.
         Make a distinction between the master rank and the other ranks. */
      if(myRank == 0) {

        /*--------------------------------------------------------------------*/
        /*             Master node, which does all the work.                  */
        /* Step 1: Create the global vector for the graph by gathering all the*/
        /*         data from the other ranks. This is done, because the greedy*/
        /*         algorithm used here is a sequential algorithm. Moreover the*/
        /*         algorithm requires the neighbors of neighbors, which would */
        /*         require a rather cumbersome communication step anyway.     */
        /**************************************************************************/

        /* Define the global vector and copy my data in it. */
        vector<vector<unsigned long> > entriesVert(nVerticesPerRank[nRank],
                                                   vector<unsigned long>(0));

        for(unsigned long i=nVerticesPerRank[0]; i<nVerticesPerRank[1]; ++i)
          entriesVert[i] = entriesVertices[i];


        /**********************************************************************/
        /* Step 2: The greedy algorithm to determine the vertex colors.       */
        /**********************************************************************/

        /* Allocate the memory to store the color of the vertices. */
        vector<int> colorVertices(nVerticesPerRank[nRank]);

        /* Allocate and reserve the memory for vectors that are used
           in the greedy algorithm. */
        vector<unsigned long> flagColorStored(nVerticesPerRank[nRank],
                                              nVerticesPerRank[nRank]);
        vector<int> colorNeighbors;
        colorNeighbors.reserve(1000);

        /* Loop over the vertices of the graph. */
        for(unsigned long i=0; i<nVerticesPerRank[nRank]; ++i) {

          /* Make sure that colorNeighbors is empty. */
          colorNeighbors.resize(0);

          /* Loop over the entries of this vertex. */
          for(unsigned long j=0; j<entriesVert[i].size(); ++j) {
            const unsigned long jj = entriesVert[i][j];

            /* Add the color of jj if jj is less than i and if its color
               is not stored yet. */
            if(jj < i) {
              const int cJJ = colorVertices[jj];
              if(flagColorStored[cJJ] != i) {
                flagColorStored[cJJ] = i;
                colorNeighbors.push_back(cJJ);
              }
            }

            /* Loop over the entries of vertex jj. */
            for(unsigned long k=0; k<entriesVert[jj].size(); ++k) {
              const unsigned long kk = entriesVert[jj][k];

              /* Add the color of kk if kk is less than i and if its color
                 is not stored yet. */
              if(kk < i) {
                const int cKK = colorVertices[kk];
                if(flagColorStored[cKK] != i) {
                  flagColorStored[cKK] = i;
                  colorNeighbors.push_back(cKK);
                }
              }
            }
          }

          /* Sort colorNeighbors in increasing order. */
          sort(colorNeighbors.begin(), colorNeighbors.end());

          /* Determine the color of this vertex. */
          int cVert;
          for(cVert=0; cVert<(int) colorNeighbors.size(); ++cVert) {
            if(cVert != colorNeighbors[cVert]) break;
          }

          colorVertices[i] = cVert;
        }

        /* Check if the coloring is done correctly. */
        flagColorStored.assign(nVerticesPerRank[nRank], nVerticesPerRank[nRank]);
        for(unsigned long i=0; i<entriesVert.size(); ++i) {
          for(unsigned long j=0; j<entriesVert[i].size(); ++j) {
            const int cJJ = colorVertices[entriesVert[i][j]];
            if(flagColorStored[cJJ] == i) {
              cout << "In function GraphVertexColoring: Color " << cJJ
                   << " appears more than once in the stencil of vertex "
                   << i << "." << endl;
              exit(1);
            }
            flagColorStored[cJJ] = i;
          }
        }

        /**********************************************************************/
        /* Step 3: Store the coloring information in the local vectors again. */
        /**********************************************************************/

        /* Store the data of the root rank. */
        colorLocalVertices.resize(nVerticesPerRank[1]);
        memcpy(colorLocalVertices.data(), colorVertices.data(),
               nVerticesPerRank[1]*sizeof(int));

        /* Send the color of the vertices to the other ranks. Only in parallel
           mode. Use blocking sends, because deadlock cannot occur. */
      }
      break;
    }

    /*-----------------------------------------------------------------------*/

    case NATURAL_COLORING: {

      /* Natural coloring, i.e. every DOF gets its own color. This is very
        inefficient and should only be used for debugging. */
      unsigned long nLocalVert = entriesVertices.size();
      colorLocalVertices.resize(nLocalVert);

      for(unsigned long i=0; i<nLocalVert; ++i)
        colorLocalVertices[i] = (int)(i + nVerticesPerRank[myRank]);

      break;
    }
  }

  /* Determine the total (global) number of colors. */
  int nLocalColors = 0;
  for(unsigned long i=0; i<colorLocalVertices.size(); ++i)
    nLocalColors = max(nLocalColors, colorLocalVertices[i]);
  ++nLocalColors;

  nGlobalColors = nLocalColors;
}
