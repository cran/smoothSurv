// Helping functions for smoothSurvReg
// ====================================

// Various functions to convert a's to c's and vice versa

// Functions to compute various derivatives
// of c's w.r.t. a's
// of a's w.r.t. d's etc.

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

#include "convertCAD.h"

using namespace SCYTHE;

// ===========================================================================
int
A_to_C(const Matrix<double> & A, Matrix<double> & C)
{
   if (A.size() != C.size()){
      REprintf("Incorrect dimensions of the input vectors in A_to_C");
      return 99;
   }

   Matrix<double> exp_a = exp(A);
   double sum_exp_a = sum(exp_a);

   if (sum_exp_a > FLT_MAX - 1) return 1;

   C = (1/sum_exp_a) * exp_a;
   return 0;
}


// ===========================================================================
// Function to compute derivatives of non-zero exp(a)'s w.r.t. exp(d)'s

//  * there are g - 1 non-zero a's
//  * there are g - 3 d's

//  * g - 3 a's are equal to d's
//  * 2 a's are a function of the rest
//  * 1 a is equal to zero

// INPUT: knots ....... vector of knots
//        sd0 ......... standard deviation of a basis G-spline
//        lastThree ... vector with indeces of the three a's which are to be computed from the rest
//             a[lastThree[0]] = 0
//             a[lastThree[1]] = first function(d's)
//             a[lastThree[2]] = second function(d's)
//        all ......... do I want the full  matrix or only two columns w.r.t.
//                      to the two a's

// OUTPUT: Omega.....
//               if all == true, matrix (g - 2) x g (there is one zero column)
//                  all == false, matrix (g - 2) x 2
//                         (the first row is always an intercept)

// RETURN: 0 if no problems
//         >0 if there are some problems

// if all == true then
//      exp(a) = t(Omega[-0,]) * exp(d) + t(Omega[0,])
// if all == false then
//      exp(a1, a2) = t(Omega[-0,]) * exp(d) + t(Omega[0,])
//          where a1 and a2 are the two a's
int
deriv_expAD(const Matrix<double> & knots, double sd0, const Matrix<int> & lastThree, Matrix<double> & Omega, bool all)
{
   int i, j, k, kk;

   int g = knots.size();
   if (g < 4){
      REprintf("Too short input 'knots' vector in deriv_expAD");
      return 1;
   }

   if (lastThree.size() != 3){
      REprintf("Incorrect 'lastThree' parameter in deriv_expAD");
      return 2;
   }

   int isInRange = 0;
   for (i = 0; i < 3; i++){
       for (k = 0; k < g; k++){
           if (lastThree[i] == k) isInRange++;
       }
   }
   if (isInRange != 3){
      REprintf("Incorrect 'lastThree' parameter in deriv_expAD");
      return 3;
   }

   if (lastThree[0] == lastThree[1] || lastThree[0] == lastThree[2] || lastThree[1] == lastThree[2]){
      REprintf("Incorrect 'lastThree' parameter in deriv_expAD");
      return 4;
   }

   if (sd0 >= 1 || sd0 <= 0){
      REprintf("Incorrect 'sd0' parameter in deriv_expAD");
      return 5;
   }

   int omegaCols = all ? g : 2;
   int omegaRows = g - 2;

   if (Omega.cols() != omegaCols || Omega.rows() != omegaRows){
      REprintf("Incorrect dimension of the matrix 'Omega' in deriv_expAD");
      return 6;
   }

   int whichZero = lastThree[0];
   int l1 = lastThree[1];
   int l2 = lastThree[2];
   double s02 = sd0 * sd0;

   if (fabs(knots[l1]) < 1e-4){
      REprintf("Zero reference knot in deriv_expAD");
      return 7;
   }

   double kn2_kn1 = knots[l2] - knots[l1];
//   if (fabs(kn2_kn1) < 1e-4){
//      REprintf("Too close reference knots in deriv_expAD");
//      return 8;
//   }

   double jsmm = 1 - s02 + knots[l1]*knots[l2];
//   if (fabs(jsmm) < 1e-4){
//      REprintf("Badly conditioned reference knots in deriv_expAD");
//      return 9;
//   }

// Knots with removed the two ones corresponding to the two special a's
   Matrix<double> knotsMin2 = Matrix<double>(g - 2, 1, false);
   for (i = 0, k = 0; k < g - 2; i++){
       if (i == l1 || i == l2) continue;
       knotsMin2[k] = knots[i];
       k++;
   }

// Index of the zero a in the shorter (by 2) a's sequence
   int minl1l2 = l1 < l2 ? l1 : l2;
   int maxl1l2 = l1 >= l2 ? l1 : l2;
   int whichZero2 = (whichZero < minl1l2 ? whichZero :
		      (whichZero < maxl1l2 ? whichZero - 1 : whichZero - 2));


// Compute the two columns of the resulting matrix corresponding to the two a's
   Matrix<double> vec2b = -(1/kn2_kn1) * (knotsMin2 - knots[l1]);
   Matrix<double> vec2c = (1/jsmm) * (1 - s02 + knots[l1] * knotsMin2);
   Matrix<double> vec2 = (vec2b & vec2c);

   double vec1a = knots[l2]/knots[l1];
   Matrix<double> vec1d = (1/knots[l1]) * knotsMin2;
   Matrix<double> vec1 = -vec1a * vec2 - vec1d;

   double int1 = vec1[whichZero2];
   double int2 = vec2[whichZero2];

   Matrix<double> slope1 = Matrix<double>(g - 3, 1, false);
   Matrix<double> slope2 = Matrix<double>(g - 3, 1, false);
   for (i = 0, k = 0; k < g - 3; i++){
       if (i == whichZero2) continue;
       slope1[k] = vec1[i];
       slope2[k] = vec2[i];
       k++;
   }

   Matrix<double> intercept = Matrix<double>(1, 2, false);
   intercept[0] = int1;
   intercept[1] = int2;
   Matrix<double> slope = cbind(slope1, slope2);

   Omega = rbind(intercept, slope);

   if (all){
      Matrix<double> leftMat = rbind(Matrix<double>(1, g - 3, true, 0.0), eye<double>(g - 3));
      Matrix<double> zeroCol = Matrix<double>(g - 2, 1, true, 0.0);
      zeroCol[0] = 1.0;
      for (j = 0, k = 0, kk = 0; j < g; j++){
         if (j == whichZero){
            Omega = cbind(Omega, zeroCol);
         }
         else{
            if (j == l1 || j == l2){
               Omega = cbind(Omega, Omega(0, k, g - 3, k));
               k++;
            }
            else{
               Omega = cbind(Omega, leftMat(0, kk, g - 3, kk));
               kk++;
            }
         }
      }
// the first two columns are now obscure -> remove them
      Omega = Omega(0, 2, g - 3, g + 1);
   }

   return 0;
}


// =================================================================================
// Function to compute derivatives of the two a's w.r.t. the remaining non-zero a's

// INPUT: exp_Acoef2 ........ exp of the two a's
//        exp_Dcoef ......... exp of the remaining (g - 3) non-zero a's
//        OmegaExpA2 ........ matrix (g - 3) x 2 with derivatives  d(exp(Acoef2))/d(exp(Dcoef))

// OUTPUT: dA2dD ............ matrix (g - 3) x 2
void
compute_dA2dD(Matrix<double> & dA2dD,
              const Matrix<double> & exp_Acoef2,
              const Matrix<double> & exp_Dcoef,
              const Matrix<double> & OmegaExpA2)
{
      dA2dD = cbind((1/exp_Acoef2[0]) * exp_Dcoef, (1/exp_Acoef2[1]) * exp_Dcoef);
      dA2dD = dA2dD & OmegaExpA2;
}


// ========================================================================================
// Function to compute second derivatives of the two a's w.r.t. the remaining non-zero a's

// It's going to be used (properly) inside penalLogLik function
// I do not do any checks for correctness of the input paramaters!!!

// INPUT: dA2dD ......... matrix (g - 3) x 2 with d(A_1, A_2)/dD
//        nD ............ number of non-zero remaining a's (g - 3)

// OUTPUT: ddAdDD1 ...... matrix (g - 3) x (g - 3) with d^2 A_1/dD dD^T
//         ddAdDD2 ...... matrix (g - 3) x (g - 3) with d^2 A_2/dD dD^T
void
compute_ddA2dDD(Matrix<double> & ddAdDD1,
                Matrix<double> & ddAdDD2,
                const Matrix<double> & dA2dD,
                const int nD)
{
   int i;

   ddAdDD1 = (-1) * dA2dD(0, 0, nD - 1, 0) * t(dA2dD(0, 0, nD - 1, 0));
   for (i = 0; i < nD; i++) ddAdDD1(i, i) += dA2dD(i, 0);

   ddAdDD2 = (-1) * dA2dD(0, 1, nD - 1, 1) * t(dA2dD(0, 1, nD - 1, 1));
   for (i = 0; i < nD; i++) ddAdDD2(i, i) += dA2dD(i, 1);
}


// ===========================================================================
// Function to compute derivatives of c's w.r.t. non-zero a's

// It's going to be used (properly) inside penalLogLik function
// I do not do any checks for correctness of the input paramaters!!!

// INPUT: Ccoef ......... c coefficients (g x 1)
//        tCcoef ........ its transposition
//        lastThreeA .... indeces of the zero a and of the two a's
//        restA ......... indeces of the remaining a's

// OUTPUT: dCdA2 ........ matrix 2 x g with dC/d(A_1, A_2)
//         dCdAg3 ....... matrix (g - 3) x g with dC/dAg3
void
compute_dCdA(Matrix<double> & dCdA2,
             Matrix<double> & dCdAg3,
             const Matrix<double> & Ccoef,
             const Matrix<double> & tCcoef,
             const Matrix<int> & lastThreeA,
             const Matrix<int> & restA,
             const int nD)
{
   int i;

   int nSplines = Ccoef.size();

// start with derivatives of c's w.r.t. all a'c (including the zero one)
   Matrix<double> dCdA = (-1) * Ccoef * tCcoef;    // matrix g x g
   for (i = 0; i < nSplines; i++) dCdA(i, i) += Ccoef[i];

// extract derivative of c's w.r.t. the two a's
   dCdA2 = dCdA(lastThreeA[1], 0, lastThreeA[1], nSplines - 1);
   dCdA2 = rbind(dCdA2, dCdA(lastThreeA[2], 0, lastThreeA[2], nSplines - 1));

// extract derivatives of c's w.r.t. the remaining g - 3 a's
   dCdAg3 = dCdA(restA[0], 0, restA[0], nSplines - 1);
   for (i = 1; i < nD; i++){
      dCdAg3 = rbind(dCdAg3, dCdA(restA[i], 0, restA[i], nSplines - 1));
   }

}


// =======================================================================
// Function to compute dC/dD

// It's going to be used (properly) inside penalLogLik function
// I do not do any checks for correctness of the input paramaters!!!

// INPUT: dCdA2 ......... matrix 2 x g with dC/d(A_1, A_2)
//        dCdAg3 ........ matrix (g - 3) x g with dC/dAg3
//        dA2dD ......... matrix (g - 3) x 2 with d(A_1, A_2)/dD
//        Ccoef ......... vactor (g x 1) with c coefficients
//        tCcoef ......... vactor (1 x g) with transposed c coefficients
//        restA ......... indeces of non-zero a's in the whole sequence
//        nD ............ number of remaining non-zero a's (g - 3) or (g - 1)
//        useD .......... which parametrization is used?

// OUTPUT: dCdD .... matrix (g - 3) x g
//                          = dC/dD
//                          =  dA_1/dD * dC/dA_1 +
//                           + dA_2/dD * dC/dA_2 +
//                           + dC/dAg3
//                 or matrix (g - 1) x g
//                          = dC/dA (if !useD)
void
compute_dCdD(Matrix<double> & dCdD,
             const Matrix<double> & dCdA2,
             const Matrix<double> & dCdAg3,
             const Matrix<double> & dA2dD,
             const Matrix<double> & Ccoef,
             const Matrix<double> & tCcoef,
             const Matrix<int> & restA,
             const int nD,
             const bool useD)
{
    int i;
    int nSplines = useD ? (nD + 3) : (nD + 1);

    if (useD){
       dCdD = dA2dD(0, 0, nD - 1, 0) * dCdA2(0, 0, 0, nSplines - 1)
             + dA2dD(0, 1, nD - 1, 1) * dCdA2(1, 0, 1, nSplines - 1)
             + dCdAg3;
    }
    else{
       Matrix<double> dCdA = (-1) * Ccoef * tCcoef;    // matrix g x g
       for (i = 0; i < nSplines; i++) dCdA(i, i) += Ccoef[i];
    // extract derivatives of c's w.r.t. the (g - 1) a's
       dCdD = dCdA(restA[0], 0, restA[0], nSplines - 1);
       for (i = 1; i < nD; i++)
          dCdD = rbind(dCdD, dCdA(restA[i], 0, restA[i], nSplines - 1));
    }


}


// ===========================================================================
// Function to compute the second derivatives of c's w.r.t. all non-zero a's

// One has to allocate properly enough space for the resulting matrices!!!
// I.e. ddCdAA must be a pointer to g matrices, each of them of size (g - 1) x (g - 1)
//   where g = length(Ccoef)

// It's going to be used (properly) inside penalLogLik function
// I do not do any checks for correctness of the input paramaters!!!

// INPUT: Ccoef ....... c coefficients
//        posZeroA .... index of the a which is zero

// OUTPUT: ddCdAA ..... an array of g matrices, each of them of type (g - 1) x (g - 1)
void
compute_ddCdAA(Matrix<double> * ddCdAA,
               const Matrix<double> & Ccoef,
               const int posZeroA)
{
   int i, j, k;

   int nSplines = Ccoef.size();
   int nA = nSplines - 1;


// ddCdAA[posZeroA] (it does not have "special" row and column)
// -----------------------------------------------------------
   for (i = 0; i < posZeroA; i++){
      ddCdAA[posZeroA](i, i) = -Ccoef[posZeroA] * Ccoef[i] * (1 - 2*Ccoef[i]);
      for (j = i + 1; j < posZeroA; j++){
          ddCdAA[posZeroA](i, j) = 2 * Ccoef[i] * Ccoef[j] * Ccoef[posZeroA];
          ddCdAA[posZeroA](j, i) = ddCdAA[posZeroA](i, j);
      }
      for (j = posZeroA; j < nA; j++){
          ddCdAA[posZeroA](i, j) = 2 * Ccoef[i] * Ccoef[j+1] * Ccoef[posZeroA];
          ddCdAA[posZeroA](j, i) = ddCdAA[posZeroA](i, j);
      }
   }
   for (i = posZeroA; i < nA; i++){
      ddCdAA[posZeroA](i, i) = -Ccoef[posZeroA] * Ccoef[i+1] * (1 - 2*Ccoef[i+1]);
      for (j = i + 1; j < nA; j++){
          ddCdAA[posZeroA](i, j) = 2 * Ccoef[i+1] * Ccoef[j+1] * Ccoef[posZeroA];
          ddCdAA[posZeroA](j, i) = ddCdAA[posZeroA](i, j);
      }
   }


// ddCdAA[k], k = 0, ..., posZeroA-1
// ---------------------------------
   for (k = 0; k < posZeroA; k++){
   // rows 0, ..., k - 1 (+ symmetry)
      for (i = 0; i < k; i++){
         ddCdAA[k](i, i) = -Ccoef[k] * Ccoef[i] * (1 - 2*Ccoef[i]);
         for (j = i + 1; j < posZeroA; j++){
             if (j == k){
                ddCdAA[k](i, k) = -Ccoef[i] * Ccoef[k] * (1 - 2*Ccoef[k]);
                ddCdAA[k](k, i) = ddCdAA[k](i, k);
             }
             else{
                ddCdAA[k](i, j) = 2 * Ccoef[i] * Ccoef[j] * Ccoef[k];
                ddCdAA[k](j, i) = ddCdAA[k](i, j);
             }
         }
         for (j = posZeroA; j < nA; j++){
             ddCdAA[k](i, j) = 2 * Ccoef[i] * Ccoef[j+1] * Ccoef[k];
             ddCdAA[k](j, i) = ddCdAA[k](i, j);
         }
      }
   // row k (special) (+ symmetry)
      ddCdAA[k](k, k) = Ccoef[k] * (1 - Ccoef[k]) * (1 - 2*Ccoef[k]);
      for (j = k + 1; j < posZeroA; j++){
          ddCdAA[k](k, j) = -Ccoef[j] * Ccoef[k] * (1 - 2*Ccoef[k]);
          ddCdAA[k](j, k) = ddCdAA[k](k, j);
      }
      for (j = posZeroA; j < nA; j++){
          ddCdAA[k](k, j) = -Ccoef[j+1] * Ccoef[k] * (1 - 2*Ccoef[k]);
          ddCdAA[k](j, k) = ddCdAA[k](k, j);
      }
   // rows k + 1, ..., posZeroA - 1 (+ symmetry)
      for (i = k + 1; i < posZeroA; i++){
         ddCdAA[k](i, i) = -Ccoef[k] * Ccoef[i] * (1 - 2*Ccoef[i]);
         for (j = i + 1; j < posZeroA; j++){
             ddCdAA[k](i, j) = 2 * Ccoef[i] * Ccoef[j] * Ccoef[k];
             ddCdAA[k](j, i) = ddCdAA[k](i, j);
         }
         for (j = posZeroA; j < nA; j++){
             ddCdAA[k](i, j) = 2 * Ccoef[i] * Ccoef[j+1] * Ccoef[k];
             ddCdAA[k](j, i) = ddCdAA[k](i, j);
         }
      }
   // rows posZeroA, ..., nA - 1 (+ symmetry)
      for (i = posZeroA; i < nA; i++){
          ddCdAA[k](i, i) = -Ccoef[k] * Ccoef[i+1] * (1 - 2*Ccoef[i+1]);
          for (j = i + 1; j < nA; j++){
              ddCdAA[k](i, j) = 2 * Ccoef[i+1] * Ccoef[j+1] * Ccoef[k];
              ddCdAA[k](j, i) = ddCdAA[k](i, j);
          }
      }
   }


// ddCdAA[k], k = posZeroA + 1, ..., nA ( = g-1)
// ---------------------------------------------
   for (k = posZeroA + 1; k < nSplines; k++){
   // rows 0, ..., posZeroA - 1 (+ symmetry)
      for (i = 0; i < posZeroA; i++){
         ddCdAA[k](i, i) = -Ccoef[k] * Ccoef[i] * (1 - 2*Ccoef[i]);
         for (j = i + 1; j < posZeroA; j++){
             ddCdAA[k](i, j) = 2 * Ccoef[i] * Ccoef[j] * Ccoef[k];
             ddCdAA[k](j, i) = ddCdAA[k](i, j);
         }
         for (j = posZeroA; j < nA; j++){
            if (j == k - 1){
               ddCdAA[k](i, k-1) = -Ccoef[i] * Ccoef[k] * (1 - 2*Ccoef[k]);
               ddCdAA[k](k-1, i) = ddCdAA[k](i, k-1);
            }
            else{
               ddCdAA[k](i, j) = 2 * Ccoef[i] * Ccoef[j+1] * Ccoef[k];
               ddCdAA[k](j, i) = ddCdAA[k](i, j);
            }
         }
      }
   // rows posZeroA, ..., k - 2 (+ symmetry)
      for (i = posZeroA; i < k - 1; i++){
         ddCdAA[k](i, i) = -Ccoef[k] * Ccoef[i+1] * (1 - 2*Ccoef[i+1]);
         for (j = i + 1; j < nA; j++){
            if (j == k - 1){
               ddCdAA[k](i, k-1) = -Ccoef[i+1] * Ccoef[k] * (1 - 2*Ccoef[k]);
               ddCdAA[k](k-1, i) = ddCdAA[k](i, k-1);
            }
            else{
               ddCdAA[k](i, j) = 2 * Ccoef[i+1] * Ccoef[j+1] * Ccoef[k];
               ddCdAA[k](j, i) = ddCdAA[k](i, j);
            }
         }
      }
   // row k - 1 (special) (+ symmetry)
      ddCdAA[k](k-1, k-1) = Ccoef[k] * (1 - Ccoef[k]) * (1 - 2*Ccoef[k]);
      for (j = k; j < nA; j++){
         ddCdAA[k](k-1, j) = -Ccoef[j+1] * Ccoef[k] * (1 - 2*Ccoef[k]);
         ddCdAA[k](j, k-1) = ddCdAA[k](k-1, j);
      }
   // rows k, ..., nA - 1 (+ symmetry)
      for (i = k; i < nA; i++){
         ddCdAA[k](i, i) = -Ccoef[k] * Ccoef[i+1] * (1 - 2*Ccoef[i+1]);
         for (j = i + 1; j < nA; j++){
            ddCdAA[k](i, j) = 2 * Ccoef[i+1] * Ccoef[j+1] * Ccoef[k];
            ddCdAA[k](j, i) = ddCdAA[k](i, j);
         }
      }
   }

}


// ==========================================================================================
// Function to compute ddCdDD

// It's going to be used (properly) inside penalLogLik function
// I do not do any checks for correctness of the input paramaters!!!

// INPUT: dCdA2 ......... matrix 2 x g with dC/d(A_1, A_2)
//        ddCdAA ........ array of g matrices, each of them of size (g - 1) x (g - 1)
//                          the jth of them should be equal to d^2 C_j/dA dA^T
//                          where A's are non-zero a's
//        dA2dD ......... matrix (g - 3) x 2 with d(A_1, A_2)/dD
//        ddAdDD1 ....... matrix (g - 3) x (g - 3) with ddA_1/dD dD^T
//        ddAdDD2 ....... matrix (g - 3) x (g - 3) with ddA_2/dD dD^T
//        lastTwoA ...... indeces of the two a's (in a shorter sequence where
//                          zero a is omitted)
//        restA ......... indeces of the remaining a's ()in a shorter sequence where
//                          zero a is omitted)

// OUTPUT: ddCdDD ........ array of g matrices (g - 3) x (g - 3)
//                             the kth one is equal to
//                             ddC/dD dD^T =
//                               dA_1/dD * d^2 C_k/dA_1 dA_1 * dA_1/dD^T +
//                             + dA_1/dD * d^2 C_k/dA_1 dA_2 * dA_2/dD^T +
//                             + dA_1/dD * d^2 C_k/dA_1 dAg3^T +
//                             + dA_2/dD * d^2 C_k/dA_2 dA_1 * dA_1/dD^T +
//                             + dA_2/dD * d^2 C_k/dA_2 dA_2 * dA_2/dD^T +
//                             + dA_2/dD * d^2 C_k/dA_2 dAg3^T +
//                             +           d^2 C_k/dAg3 dA_1^T * dA_1/dD^T +
//                             +           d^2 C_k/dAg3 dA_2^T * dA_2/dD^T +
//                             +           d^2 C_k/dAg3 dAg3^T
//                             + ddA_1/dD dD^T * dC_k/dA_1 +
//                             + ddA_2/dD dD^T * dC_k/dA_2 +
void
compute_ddCdDD(Matrix<double> * ddCdDD,
               const Matrix<double> & dCdA2,
               const Matrix<double> * ddCdAA,
               const Matrix<double> & dA2dD,
               const Matrix<double> & ddAdDD1,
               const Matrix<double> & ddAdDD2,
               const Matrix<int> & lastTwoA,
               const Matrix<int> & restA)
{

   int nSplines = lastTwoA.size() + restA.size() + 1;
   int nD = nSplines - 3;

   for (int k = 0; k < nSplines; k++){
      for (int i = 0; i < nD; i++){
      // diagonal
         ddCdDD[k](i, i) =
               2 * (dA2dD(i, 0) * ddCdAA[k](lastTwoA[0], lastTwoA[1]) * dA2dD(i, 1) +
                    dA2dD(i, 0) * ddCdAA[k](lastTwoA[0], restA[i]) +
                    dA2dD(i, 1) * ddCdAA[k](lastTwoA[1], restA[i])) +
               dA2dD(i, 0) * ddCdAA[k](lastTwoA[0], lastTwoA[0]) * dA2dD(i, 0) +
               dA2dD(i, 1) * ddCdAA[k](lastTwoA[1], lastTwoA[1]) * dA2dD(i, 1) +
               ddCdAA[k](restA[i], restA[i]) +
               dCdA2(0, k) * ddAdDD1(i, i) +
               dCdA2(1, k) * ddAdDD2(i, i);
      // rest
         for (int j = i + 1; j < nD; j++){
            ddCdDD[k](i, j) =
                 dA2dD(i, 0) * ddCdAA[k](lastTwoA[0], lastTwoA[0]) * dA2dD(j, 0) +
                 dA2dD(i, 0) * ddCdAA[k](lastTwoA[0], lastTwoA[1]) * dA2dD(j, 1) +
                 dA2dD(i, 0) * ddCdAA[k](lastTwoA[0], restA[j]) +
                 dA2dD(i, 1) * ddCdAA[k](lastTwoA[1], lastTwoA[0]) * dA2dD(j, 0) +
                 dA2dD(i, 1) * ddCdAA[k](lastTwoA[1], lastTwoA[1]) * dA2dD(j, 1) +
                 dA2dD(i, 1) * ddCdAA[k](lastTwoA[1], restA[j]) +
                 ddCdAA[k](restA[i], lastTwoA[0]) * dA2dD(j, 0) +
                 ddCdAA[k](restA[i], lastTwoA[1]) * dA2dD(j, 1) +
                 ddCdAA[k](restA[i], restA[j]) +
                 dCdA2(0, k) * ddAdDD1(i, j) +
                 dCdA2(1, k) * ddAdDD2(i, j);
            ddCdDD[k](j, i) = ddCdDD[k](i, j);
         }
      }
   }
}

