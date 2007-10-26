// Grabbed from FORTRAN routines by Arnost Komarek in May 2003

// SET OF FORTRAN ROUTINES FROM EISPACK to compute eigen vectors and values
// REWRITTEN TO C++
// ========================================================================
// * only routines for real symmetric matrices
// * when only eigen values and not eigen vectors are asked
//   the function rsCPP may not work correctly
//   there is still probably a small error in the subfunction 'tqlratCPP'

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

#include "eispackFORT.h"

extern "C"{

// ------------------------------------------------------------------
//
//     THIS SUBROUTINE CALLS THE RECOMMENDED SEQUENCE OF
//     SUBROUTINES FROM THE EIGENSYSTEM SUBROUTINE PACKAGE (EISPACK)
//     TO FIND THE EIGENVALUES AND EIGENVECTORS (IF DESIRED)
//     OF A REAL SYMMETRIC MATRIX.
//
//     ON INPUT
//
//        NM  MUST BE SET TO THE ROW DIMENSION OF THE TWO-DIMENSIONAL
//        ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
//        DIMENSION STATEMENT.
//
//        N  IS THE ORDER OF THE MATRIX  A.
//
//        A  CONTAINS THE REAL SYMMETRIC MATRIX.
//
//        MATZ  IS AN INTEGER VARIABLE SET EQUAL TO ZERO IF
//        ONLY EIGENVALUES ARE DESIRED.  OTHERWISE IT IS SET TO
//        ANY NON-ZERO INTEGER FOR BOTH EIGENVALUES AND EIGENVECTORS.
//
//     ON OUTPUT
//
//        W  CONTAINS THE EIGENVALUES IN ASCENDING ORDER.
//
//        Z  CONTAINS THE EIGENVECTORS IF MATZ IS NOT ZERO.
//
//        IERR  IS AN INTEGER OUTPUT VARIABLE SET EQUAL TO AN ERROR
//           COMPLETION CODE DESCRIBED IN THE DOCUMENTATION FOR TQLRAT
//           AND TQL2.  THE NORMAL COMPLETION CODE IS ZERO.
//
//        FV1  AND  FV2  ARE TEMPORARY STORAGE ARRAYS.
//
//     THIS VERSION DATED AUGUST 1983.
//
//     ------------------------------------------------------------------
//
void
rsCPP(int* NM, int* N, double* A, double* W, int* MATZ, double* Z, int* IERR)
{

//    dimensions: A(NM,N), W(N), Z(NM,N)

      if (*N  >  *NM){
        *IERR = 10 * (*N);
        return;
      }

      double * FV1 = new double[*N];
      double * FV2 = new double[*N];

      if (*MATZ == 0){
//     .......... FIND EIGENVALUES ONLY ..........
         tred1CPP(*NM, *N, A, W, FV1, FV2);
         tqlratCPP(*N, W, FV2, IERR);
      }
      else{
//     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........
         tred2CPP(*NM, *N, A, W, FV1, Z);
         tql2CPP(*NM, *N, W, FV1, Z, IERR);
      }

      delete [] FV1;
      delete [] FV2;
      return;
}

// ------------------------------------------------------------------
//
//     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED2,
//     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
//     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
//
//     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX TO A
//     SYMMETRIC TRIDIAGONAL MATRIX USING AND ACCUMULATING
//     ORTHOGONAL SIMILARITY TRANSFORMATIONS.
//
//     i.e. A = Z * DE * t(Z)
//          where DE is a tridiagonal matrix
//          and Z is the orthogonal transformation matrix
//
//     ON INPUT
//
//        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
//          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
//          DIMENSION STATEMENT.
//
//        N IS THE ORDER OF THE MATRIX.
//
//        A CONTAINS THE REAL SYMMETRIC INPUT MATRIX.  ONLY THE
//          LOWER TRIANGLE OF THE MATRIX NEED BE SUPPLIED.
//
//     ON OUTPUT
//
//        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX.
//
//        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
//          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO.
//
//        Z CONTAINS THE ORTHOGONAL TRANSFORMATION MATRIX
//          PRODUCED IN THE REDUCTION.
//
//        A AND Z MAY COINCIDE.  IF DISTINCT, A IS UNALTERED.
//
//     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
//     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
//
//     THIS VERSION DATED AUGUST 1983.
//
//     ------------------------------------------------------------------
void
tred2CPP(int NM, int N, double* A, double* D, double* E, double *Z)
{

      int I, J, K, L, II, JP1;
// dimensions: A(NM,N), D(N), E(N), Z(NM,N)
      double F, G, H, HH, SCALE;
      bool helpbool;

      for (I = 1; I <= N; I++){
         for (J = I; J <= N; J++){
           Z[(I-1)*NM + J-1] = A[(I-1)*NM + J-1];
           D[I-1] = A[(I-1)*NM + N-1];
         }
     }

     if (N == 1){
        for (I = 1; I <= N; I++){
           D[I-1] = Z[(I-1)*NM + N-1];
           Z[(I-1)*NM + N-1] = 0.0;
        }
        Z[(N-1)*NM + N-1] = 1.0;
        E[0] = 0.0;
        return;
     }

//     .......... FOR I=N STEP -1 UNTIL 2 DO -- ..........
     for (II = 2; II <= N; II++){   // DO 300
         I = N + 2 - II;
         L = I - 1;
         H = 0.0;
         SCALE = 0.0;
         helpbool = true;
         if (L >= 2){

//     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........
            for (K = 1; K <= L; K++)
               SCALE +=  fabs(D[K-1]);

            if (SCALE != 0.0){
               helpbool = false;
               for (K = 1; K <= L; K++){
                  D[K-1] /= SCALE;
                  H += D[K-1] * D[K-1];
               }

               F = D[L-1];
               G = (F < 0 ? sqrt(H) : -sqrt(H));
               E[I-1] = SCALE * G;
               H -= F * G;
               D[L-1] = F - G;

//     .......... FORM A*U ..........
               for (J = 1; J <= L; J++)
                  E[J-1] = 0.0;

               for (J = 1; J <= L; J++){
                  F = D[J-1];
                  Z[(I-1)*NM + J-1] = F;
                  G = E[J-1] + Z[(J-1)*NM + J-1] * F;
                  JP1 = J + 1;
                  if (L >= JP1){
                     for (K = JP1; K <= L; K++){
                        G += Z[(J-1)*NM + K-1] * D[K-1];
                        E[K-1] +=  Z[(J-1)*NM + K-1] * F;
                     }
                  }
                  E[J-1] = G;
               }

//     .......... FORM P ..........
               F = 0.0;
               for (J = 1; J <= L; J++){
                  E[J-1] /=  H;
                  F += E[J-1] * D[J-1];
               }

               HH = F / (H + H);

//     .......... FORM Q ..........
               for (J = 1; J <= L; J++) E[J-1] -= HH * D[J-1];

//     .......... FORM REDUCED A ..........
               for (J = 1; J <= L; J++){
                  F = D[J-1];
                  G = E[J-1];

                  for (K = J; K <= L; K++)
                     Z[(J-1)*NM + K-1] -= (F * E[K-1] + G * D[K-1]);

                  D[J-1] = Z[(J-1)*NM + L-1];
                  Z[(J-1)*NM + I-1] = 0.0;
               }

            }   // end of if (SCALE != 0.0)
         }      // end of if (L >= 2)

         if (helpbool){    // do this if either (L < 2) or (SCALE == 0.0)
            E[I-1] = D[L-1];

            for (J = 1; J <= L; J++){
               D[J-1] = Z[(J-1)*NM + L-1];
               Z[(J-1)*NM + I-1] = 0.0;
               Z[(I-1)*NM + J-1] = 0.0;
            }
         }

         D[I-1] = H;

     }   // end of DO 300

//     .......... ACCUMULATION OF TRANSFORMATION MATRICES ..........
     for (I = 2; I <= N; I++){   // DO 500
         L = I - 1;
         Z[(L-1)*NM + N-1] = Z[(L-1)*NM + L-1];
         Z[(L-1)*NM + L-1] = 1.0;
         H = D[I-1];
         if (H != 0.0){
            for (K = 1; K <= L; K++)
               D[K-1] = Z[(I-1)*NM + K-1] / H;
            for (J = 1; J <= L; J++){
               G = 0.0;
               for (K = 1; K <= L; K++)
                  G += Z[(I-1)*NM + K-1] * Z[(J-1)*NM + K-1];
               for (K = 1; K <= L; K++)
                  Z[(J-1)*NM + K-1] -= G * D[K-1];
            }
         }

         for (K = 1; K <= L; K++)     // LABEL 380
            Z[(I-1)*NM + K-1] = 0.0;
     }   // END of DO 500

     for (I = 1; I <= N; I++){
        D[I-1] = Z[(I-1)*NM + N-1];
        Z[(I-1)*NM + N-1] = 0.0;
     }
     Z[(N-1)*NM + N-1] = 1.0;
     E[0] = 0.0;
     return;
}

// ------------------------------------------------------------------
//
//
//     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQL2,
//     NUM. MATH. 11, 293-306(1968) BY BOWDLER, MARTIN, REINSCH, AND
//     WILKINSON.
//     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).
//
//     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS
//     OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE QL METHOD.
//     THE EIGENVECTORS OF A FULL SYMMETRIC MATRIX CAN ALSO
//     BE FOUND IF  TRED2  HAS BEEN USED TO REDUCE THIS
//     FULL MATRIX TO TRIDIAGONAL FORM.
//
//     ON INPUT
//
//        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
//          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
//          DIMENSION STATEMENT.
//
//        N IS THE ORDER OF THE MATRIX.
//
//        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
//
//        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
//          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.
//
//        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
//          REDUCTION BY  TRED2, IF PERFORMED.  IF THE EIGENVECTORS
//          OF THE TRIDIAGONAL MATRIX ARE DESIRED, Z MUST CONTAIN
//          THE IDENTITY MATRIX.
//
//      ON OUTPUT
//
//       D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
//          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT
//          UNORDERED FOR INDICES 1,2,...,IERR-1.
//
//        E HAS BEEN DESTROYED.
//
//        Z CONTAINS ORTHONORMAL EIGENVECTORS OF THE SYMMETRIC
//          TRIDIAGONAL (OR FULL) MATRIX.  IF AN ERROR EXIT IS MADE,
//          Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED
//          EIGENVALUES.
//
//        IERR IS SET TO
//          ZERO       FOR NORMAL RETURN,
//          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
//                     DETERMINED AFTER 30 ITERATIONS.
//
//     CALLS PYTHAG FOR  DSQRT(A*A + B*B) .
//
//     THIS VERSION DATED AUGUST 1983.
//
//     ------------------------------------------------------------------
//
//     unnecessary initialization of C3 and S2 to keep g77 -Wall happy
//
void
tql2CPP(int NM, int N, double* D, double* E, double* Z, int* IERR)
{
      int I, J, K, L, M, II, L1, L2, MML;
      double C, C2, C3, DL1, EL1, F, G, H, P, R, S, S2, TST1, TST2;

// dimensions: D(N), E(N), Z(NM,N)

      C3 = 0.0;
      S2 = 0.0;

      *IERR = 0;
      if (N == 1) return;

      for (I = 2; I <= N; I++)
          E[I-2] = E[I-1];

      F = 0.0;
      TST1 = 0.0;
      E[N-1] = 0.0;

      for (L = 1; L <= N; L++){   // DO 240
         J = 0;
         H = fabs(D[L-1]) + fabs(E[L-1]);
         if (TST1 < H) TST1 = H;

//     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........
         TST2 = 0.0;
         for (M = L; M <= N; M++){
            TST2 = TST1 + fabs(E[M-1]);
            if (TST2 == TST1) break;
//    .......... E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
//               THROUGH THE BOTTOM OF THE LOOP ..........
         }
         if (M == L){
            D[L-1] += F;
            continue;
         }

         while(TST2 > TST1 || J == 0){     // this is satisfied for the first time
            if (J == 30){    // 130
//     .......... SET ERROR -- NO CONVERGENCE TO AN
//                EIGENVALUE AFTER 30 ITERATIONS ..........
              *IERR = L;
              return;
            }
            J++;

//     .......... FORM SHIFT ..........
            L1 = L + 1;
            L2 = L1 + 1;
            G = D[L-1];
            P = (D[L1-1] - G) / (2.0 * E[L-1]);
            R = pythagCPP(P, 1.0);
            D[L-1] = E[L-1] / (P + (P < 0 ? -R : R));
            D[L1-1] = E[L-1] * (P + (P < 0 ? -R : R));
            DL1 = D[L1-1];
            H = G - D[L-1];
            if (L2 <= N){
               for (I = L2; I <= N; I++) D[I-1] -= H;
            }
            F += H;

//     .......... QL TRANSFORMATION ..........
            P = D[M-1];
            C = 1.0;
            C2 = C;
            EL1 = E[L1-1];
            S = 0.0;
            MML = M - L;

//     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
            for (II = 1; II <= MML; II++){   // DO 200
               C3 = C2;
               C2 = C;
               S2 = S;
               I = M - II;
               G = C * E[I-1];
               H = C * P;
               R = pythagCPP(P, E[I-1]);
               E[I] = S * R;
               S = E[I-1] / R;
               C = P / R;
               P = C * D[I-1] - S * G;
               D[I] = H + S * (C * G + S * D[I-1]);

//     .......... FORM VECTOR ..........
               for (K = 1; K <= N; K++){
                  H = Z[I * NM + K-1];
                  Z[I * NM + K-1] = S * Z[(I-1) * NM + K-1] + C * H;
                  Z[(I-1) * NM + K-1] = C * Z[(I-1) * NM + K-1] - S * H;
               }
            }      // end of DO 200

            P = -S * S2 * C3 * EL1 * E[L-1] / DL1;
            E[L-1] = S * P;
            D[L-1] = C * P;
            TST2 = TST1 + fabs(E[L-1]);
         }    // end of while

         D[L-1] += F;
      }    // end of DO 240


//     .......... ORDER EIGENVALUES AND EIGENVECTORS ..........
      for (II = 2; II <= N; II++){   // DO 300
         I = II - 1;
         K = I;
         P = D[I-1];

         for (J = II; J <= N; J++){
            if (D[J-1] >= P) continue;
            K = J;
            P = D[J-1];
         }

         if (K == I) continue;
         D[K-1] = D[I-1];
         D[I-1] = P;

         for (J = 1; J <= N; J++){
            P = Z[(I-1)*NM + J-1];
            Z[(I-1)*NM + J-1] = Z[(K-1)*NM + J-1];
            Z[(K-1)*NM + J-1] = P;
         }
      }
      return;
}

//     ------------------------------------------------------------------
//
//     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED1,
//     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
//     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
//
//     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX
//     TO A SYMMETRIC TRIDIAGONAL MATRIX USING
//     ORTHOGONAL SIMILARITY TRANSFORMATIONS.
//
//     ON INPUT
//
//        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
//          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
//          DIMENSION STATEMENT.
//
//        N IS THE ORDER OF THE MATRIX.
//
//        A CONTAINS THE REAL SYMMETRIC INPUT MATRIX.  ONLY THE
//          LOWER TRIANGLE OF THE MATRIX NEED BE SUPPLIED.
//
//     ON OUTPUT
//
//        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL TRANS-
//          FORMATIONS USED IN THE REDUCTION IN ITS STRICT LOWER
//          TRIANGLE.  THE FULL UPPER TRIANGLE OF A IS UNALTERED.
//
//        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX.
//
//        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
//          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO.
//
//        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
//          E2 MAY COINCIDE WITH E IF THE SQUARES ARE NOT NEEDED.
//
//     THIS VERSION DATED AUGUST 1983.
//
//   ------------------------------------------------------------------
void
tred1CPP(int NM, int N, double* A, double* D, double* E, double* E2)
{

      int I, J, K, L, II, JP1;
      double F, G, H, SCALE;

// dimensions: A(NM,N), D(N), E(N), E2(N)

      for (I = 1; I <= N; I++){
         D[I-1] = A[(I-1)*NM + N-1];
         A[(I-1)*NM + N-1] = A[(I-1)*NM + I-1];
     }

//     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........
     for (II = 1; II <= N; II++){   // DO 300
         I = N + 1 - II;
         L = I - 1;
         H = 0.0;
         SCALE = 0.0;
         if (L >= 1){
            for (K = 1; K <= L; K++) SCALE += fabs(D[K-1]);

//     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........

            if (SCALE == 0.0){
               for (J = 1; J <= L; J++){
                  D[J-1] = A[(J-1)*NM + L-1];
                  A[(J-1)*NM + L-1] = A[(J-1)*NM + I-1];
                  A[(J-1)*NM + I-1] = 0.0;
               }
               E[I-1] = 0.0;
               E2[I-1] = 0.0;
               continue;
            }
         }    // end of if (L >= 1)
         else{  // i.e. if L < 1
            E[I-1] = 0.0;
            E2[I-1] = 0.0;
            continue;
         }

         // the following is done if L >= 1 && SCALE != 0
         for (K = 1; K <= L; K++){
            D[K-1] /=  SCALE;
            H += D[K-1] * D[K-1];
         }

         E2[I-1] = SCALE * SCALE * H;
         F = D[L-1];
         G = (F < 0 ? sqrt(H) : -sqrt(H));
         E[I-1] = SCALE * G;
         H -= F * G;
         D[L-1] = F - G;

         if (L == 1){
            for (J = 1; J <= L; J++){
               F = D[J-1];
               D[J-1] = A[(J-1)*NM + L-1];
               A[(J-1)*NM + L-1] = A[(J-1)*NM + I-1];
               A[(J-1)*NM + I-1] = F * SCALE;
            }
            continue;
         }

//     .......... FORM A*U ..........
         for (J = 1; J <= L; J++) E[J-1] = 0.0;

         for (J = 1; J <= L; J++){   // DO 240
            F = D[J-1];
            G = E[J-1] + A[(J-1)*NM + J-1] * F;
            JP1 = J + 1;
            if (L >= JP1)
               for (K = JP1; K <= L; K++){
                  G += A[(J-1)*NM + K-1] * D[K-1];
                  E[K-1] += A[(J-1)*NM + K-1] * F;
               }
            E[J-1] = G;
         }   // end of DO 240

//     .......... FORM P ..........
         F = 0.0;

         for (J = 1; J <= L; J++){
            E[J-1] /= H;
            F += E[J-1] * D[J-1];
         }

         H = F / (H + H);

//     .......... FORM Q ..........
         for (J = 1; J <= L; J++) E[J-1] -= H * D[J-1];

//     .......... FORM REDUCED A ..........
         for (J = 1; J <= L; J++){   // DO 280
            F = D[J-1];
            G = E[J-1];
            for (K = J; K <= L; K++)
                A[(J-1)*NM + K-1] -= (F * E[K-1] + G * D[K-1]);
         }   // end of DO 280

         for (J = 1; J <= L; J++){
            F = D[J-1];
            D[J-1] = A[(J-1)*NM + L-1];
            A[(J-1)*NM + L-1] = A[(J-1)*NM + I-1];
            A[(J-1)*NM + I-1] = F * SCALE;
         }

     }   // end of DO 300
     return;
}


//     ------------------------------------------------------------------
//
//     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQLRAT,
//     ALGORITHM 464, COMM. ACM 16, 689(1973) BY REINSCH.
//
//     THIS SUBROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC
//     TRIDIAGONAL MATRIX BY THE RATIONAL QL METHOD.
//
//     ON INPUT
//
//        N IS THE ORDER OF THE MATRIX.
//
//        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
//
//        E2 CONTAINS THE SQUARES OF THE SUBDIAGONAL ELEMENTS OF THE
//          INPUT MATRIX IN ITS LAST N-1 POSITIONS.  E2(1) IS ARBITRARY.
//
//     ON OUTPUT
//
//        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
//          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND
//          ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE
//          THE SMALLEST EIGENVALUES.
//
//        E2 HAS BEEN DESTROYED.
//
//        IERR IS SET TO
//          ZERO       FOR NORMAL RETURN,
//          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
//                     DETERMINED AFTER 30 ITERATIONS.
//
//     CALLS PYTHAG FOR  DSQRT(A*A + B*B) .
//
//     THIS VERSION DATED AUGUST 1983.
//
//     ------------------------------------------------------------------
//
//     unnecessary initialization of B and C to keep g77 -Wall happy
//
//     ------------------------------------------------------------------
void
tqlratCPP(int N, double* D, double* E2, int* IERR)
{

      int I, J, L, M, II, L1, MML;
      double B, C, F, G, H, P, R, S, T;
      bool iterate;

// dimensions: D(N), E2(N)

      B = 0.0;
      C = 0.0;

      *IERR = 0;
      if (N == 1) return;

      for (I = 2; I <= N; I++)
         E2[I-2] = E2[I-1];

      F = 0.0;
      T = 0.0;
      E2[N-1] = 0.0;

      for (L = 1; L <= N; L++){    // DO 290
         J = 0;
         H = fabs(D[L-1]) + sqrt(E2[L-1]);
         if (T <= H){
            T = H;
            B = epslonCPP(T);
            C = B * B;
         }

//     .......... LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT ..........
         for (M = L; M <= N; M++){  // DO 110
            if (E2[M-1] <= C){
                  break;
            }
//     .......... E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
//                THROUGH THE BOTTOM OF THE LOOP ..........
         }   // end of DO 110

         iterate = true;
         if (M != L){    // else goto 210; */
            while (iterate){
               if (J == 30){
//     .......... SET ERROR -- NO CONVERGENCE TO AN
//                EIGENVALUE AFTER 30 ITERATIONS ..........
                  *IERR = L;
                  return;
               }
               J++;

//     .......... FORM SHIFT ..........
               L1 = L + 1;
               S = sqrt(E2[L-1]);
               G = D[L-1];
               P = (D[L1-1] - G) / (2.0 * S);
               R = pythagCPP(P, 1.0);
               D[L-1] = S / (P + (P < 0 ? -R : R));
               H = G - D[L-1];

               for (I = L1; I <= N; I++)
                  D[I-1] -= H;

               F += H;

//     .......... RATIONAL QL TRANSFORMATION ..........
               G = D[M-1];
               if (G == 0.0) G = B;
               H = G;
               S = 0.0;
               MML = M - L;

//     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
               for (II = 1; II <= MML; II++){
                  I = M - II;
                  P = G * H;
                  R = P + E2[I-1];
                  E2[I] = S * R;
                  S = E2[I-1] / R;
                  D[I] = H + S * (H + D[I-1]);
                  G = D[I-1] - E2[I-1] / G;
                  if (G == 0.0) G = B;
                  H = G * P / R;
               }

               E2[L-1] = S * G;
               D[L-1] = H;

//     .......... GUARD AGAINST UNDERFLOW IN CONVERGENCE TEST ..........
               if (H != 0.0)
                  if (fabs(E2[L-1]) > fabs(C/H)){
                     E2[L-1] *= H;
                     if (E2[L-1] == 0.0) iterate = false;
                  }
                  else iterate = false;
               else iterate = false;
            }  // end of while (iterate)
         }    // end of if (M != L)  */

         P = D[L-1] + F;

//     .......... ORDER EIGENVALUES ..........
         if (L != 1){   // else goto 250

//     .......... FOR I=L STEP -1 UNTIL 2 DO -- ..........
            for (II = 2; II <= L; II++){
               I = L + 2 - II;
               if (P >= D[I-2]) break;
               D[I-1] = D[I-2];
            }
            if (P < D[I-2]) I = 1;
         }
         else{   // i.e. if L == 1
            I = 1;
         }
         D[I-1] = P;

      }    // end of DO 290
      return;
}


//     ------------------------------------------------------------------
//
//     FINDS DSQRT(A*A + B*B) WITHOUT OVERFLOW OR DESTRUCTIVE UNDERFLOW
//
//     ------------------------------------------------------------------
double
pythagCPP(double A, double B)
{
      double P, R, S, T, U;

      P = (fabs(A) > fabs(B) ? fabs(A) : fabs(B));
      if (P == 0.0) return 0.0;
      R = (fabs(A) < fabs(B) ? fabs(A) : fabs(B))/P;
      R *= R;

      T = 4.0 + R;
      while (T != 4.0){
         S = R / T;
         U = 1.0 + 2.0 * S;
         P = U * P;
         R *= (S / U) * (S / U);
         T = 4.0 + R;
      }

      return P;
}


//     ------------------------------------------------------------------
//
//     ESTIMATE UNIT ROUNDOFF IN QUANTITIES OF SIZE X.
//
//     THIS PROGRAM SHOULD FUNCTION PROPERLY ON ALL SYSTEMS
//     SATISFYING THE FOLLOWING TWO ASSUMPTIONS,
//        1.  THE BASE USED IN REPRESENTING FLOATING POINT
//            NUMBERS IS NOT A POWER OF THREE.
//        2.  THE QUANTITY  A  IN STATEMENT 10 IS REPRESENTED TO
//            THE ACCURACY USED IN FLOATING POINT VARIABLES
//            THAT ARE STORED IN MEMORY.
//     THE STATEMENT NUMBER 10 AND THE GO TO 10 ARE INTENDED TO
//     FORCE OPTIMIZING COMPILERS TO GENERATE CODE SATISFYING
//     ASSUMPTION 2.
//     UNDER THESE ASSUMPTIONS, IT SHOULD BE TRUE THAT,
//            A  IS NOT EXACTLY EQUAL TO FOUR-THIRDS,
//            B  HAS A ZERO FOR ITS LAST BIT OR DIGIT,
//            C  IS NOT EXACTLY EQUAL TO ONE,
//            EPS  MEASURES THE SEPARATION OF 1.0 FROM
//                 THE NEXT LARGER FLOATING POINT NUMBER.
//     THE DEVELOPERS OF EISPACK WOULD APPRECIATE BEING INFORMED
//     ABOUT ANY SYSTEMS WHERE THESE ASSUMPTIONS DO NOT HOLD.
//
//     THIS VERSION DATED 4/6/83.
//
//     ------------------------------------------------------------------
double
epslonCPP(double X)
{
      double A, B, C, EPS;

      A = 4.0 / 3.0;
      EPS = 0.0;
      while (EPS == 0.0){
         B = A - 1.0;
         C = B + B + B;
         EPS = fabs(C - 1.0);
      }

      return (EPS * fabs(X));
}

}     // end of extern "C"


