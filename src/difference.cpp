// Stuff for working with differences

// This code written by:
//
// Arnost Komarek
//
// Biostatistical Centre
// Katholieke Universiteit Leuven
// Kapucijnenvoer 35
// B - 3000, Leuven
// Belgium
//
// arnost.komarek@med.kuleuven.ac.be
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

#include "difference.h"

using namespace SCYTHE;

// Return D matrix ((g - k) x g) to compute
// k-th ordered differences of vector A (g x 1) as
// Delta{k} = D*A
// INPUT:
//   g    - length of vector A
//   k    - order of the difference I want
//   !!!  g > k  !!!
// OUTPUT:
//   D(k) - matrix to compute the difference
// THEORY:
//   Delta{k+1}A = Delta{k}A[1:(g-1)] - Delta{k}A[0:(g-2)]
//               = D(k;g-1)*A[1:(g-1)] - D(k;g-1)*A[0:(g-2)]
//               = [D(k;g-1)*H(g) - D(k;g-1)*G(g)] * A
//     where H(g)*A = A[1:(g-1)] ==> H(g) is matrix (g-1) x g with I(g-1) in lower right corner
//           G(g)*A = A[0:(g-2)] ==> G(g) is matrix (g-1) x g with I(g-1) in upper left corner
Matrix<double>
D_operator(const int g, const int k)
{
   Matrix<double> D;
   //int gg = g;

   if (k == 0){
      D = eye<double>(g);
      return D;
   }
   else{
     Matrix<double> nullcol(g-1, 1, true, 0);
     Matrix<double> ident = eye<double>(g-1);
     Matrix<double> H = cbind(nullcol, ident);
     Matrix<double> G = cbind(ident, nullcol);
     D = D_operator(g-1, k-1) * (H - G);
     return D;
   }
}

// Compute k-th ordered differences of elements of vector A
// INPUT:
//   A    - matrix g x 1 (column vector)
//   k    - order of the difference I want
//   !!!  g > k  !!!
// OUTPUT:
//   B    - matrix (g-k) x 1 with ordered differences (I return always column vector)
//           B[0] = Delta{k} A[k]
//           B[1] = Delta{k} A[k+1]
//           ...
//           B[g-k-1] = Delta{k} A[g-1]
Matrix<double>
Diff(const Matrix<double> & A, const int k)
{

   int m = A.rows();
   int n = A.cols();
   int dimA = m * n;
   //int dimB = dimA - k;

   Matrix<double> B;

   if (k == 0){
     B = A;
     return B;
   }
   else{
     Matrix<double> Aleft = A(1, 0, dimA-1, 0);
     Matrix<double> Aright = A(0, 0, dimA-2, 0);
     B = Diff(Aleft, k-1) - Diff(Aright, k-1);
     return B;
   }
}




