// Routine to adjust a SYMMETRIC matrix which is not positive definite
//  to a matrix which will be positive definite
//  via a truncated eigen-values  back composition

// This code written by:
//
// Arnost Komarek
//
// Dept. of Probability and Mathematical Statistics
// Charles University
// Sokolovska 83
// CZ - 186 75, Praha 8
// the Czech Republic
//
// komarek@karlin.mff.cuni.cz
//
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation; either version 2 of the
// License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
//

#include "createPosDef.h"

// INPUT:
//    H ...... a SYMMETRIC matrix stored in an array in a column-major order (as in S)
//             (will be destroyed on the exit)
//
//    n ...... dimension of the H matrix, it is supposed to be n x n
//
//    eps .... the value to be used in the back composition of the matrix
//             instead of the eigen values lower than eps
//             (it's usually a small positive number)

// OUTPUT:
//    H ...... positive definite matrix obtained from H

// RETURN:
//    0 if no problems with a computation of eigen values and vectors
//    IERR != 0 returned by 'rsCPP' routine

int
createPosDef(double * H, const int n, const double eps)
{
   int i, j, k, l;
   int ierr = 0;
   int * ierrP = new int;
   int * nP = new int;
   *nP = n;
   int * matzP = new int;
   *matzP = 1;                             // I want to compute both eigen vectors and values

   double * lambdas = new double[n];      // array to store eigen values (in ascending order)
   double * vectors = new double[n*n];    // array to store eigen vectors

   // Compute eigen values and eigen vectors
   rsCPP(nP, nP, H, lambdas,  matzP, vectors, ierrP);

   if (*ierrP != 0)
      ierr = *ierrP;
   else{
      i = n - 1;
      // Back composition using sufficiently large eigen values
      while (lambdas[i] > eps && i > -1){
         for (j = 0; j < n*n; j++){
            k = j % n;              // row
            l = j / n;              // column
            H[j] += lambdas[i] * vectors[i*n + k] * vectors[i*n + l];
         }
         i--;
      }
      // Replacement of small and negative eigen values by eps and back comnposition
      for (i; i > -1; i--){
         for (j = 0; j < n*n; j++){
            k = j % n;              // row
            l = j / n;              // column
            H[j] += eps * vectors[i*n + k] * vectors[i*n + l];
         }
      }
   }

   delete ierrP;
   delete nP;
   delete matzP;
   delete [] lambdas;
   delete [] vectors;

   return ierr;
}

