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
#ifndef AK_MATRIX_CPP
#define AK_MATRIX_CPP

#include "AKmatrix.h"

#include "Scythe_Error.h"
#include "Scythe_Util.h"
#include "Scythe_Stat.h"
#include <algorithm>
#include <cmath>
#include <set>

namespace SCYTHE {

// sort the elements of matrix A
// return the ordering in the matrix ordering
// return the inverse ordering in the matrix invordering
// ordering must be of the same dimension as A

// i.e. ordering[i]    = rank(A[i])
//      invordering[i] = j          <=> rank(A[j]) = i
//            (j is the index of the ith element)

//    -- it sorts actually the array where the elements of A are stored
//       i.e. it uses row major order
template <class T>
Matrix<T> sortOrder(const Matrix<T> & A, Matrix<int> & ordering, Matrix<int> & invordering)
{
  if (ordering.rows() != A.rows() || ordering.cols() != A.cols())
       throw scythe_dimension_error (__FILE__, __AK_PRETTY_FUNCTION__,
                                    __LINE__, "A and ordering of different dimension");
  if (invordering.rows() != A.rows() || invordering.cols() != A.cols())
       throw scythe_dimension_error (__FILE__, __AK_PRETTY_FUNCTION__,
                                     __LINE__, "A and ordering of different dimension");

   Matrix<T> sorted = A;
   stable_sort(sorted.begin(), sorted.end());
   for (int i = 0; i < A.size(); i++){
      ordering[i] = 0;
      for (int j = 0; j < A.size(); j++){
          if (A[j] < A[i]){
              ordering[i]++;
          }
          else
             if (A[j] == A[i] && j < i){
                ordering[i]++;
             }
      }
      invordering[ordering[i]] = i;
   }

   return sorted;
}

}   // namespace SCYTHE

#endif

