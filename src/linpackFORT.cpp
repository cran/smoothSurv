// Grabbed from FORTRAN routines by Arnost Komarek in May 2003

// SET OF FORTRAN ROUTINES FROM LINPACK USED IN 'quadprog' R library
// REWRITTEN TO C++
// ==================================================================

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

#include "linpackFORT.h"

// =============================================================
//
//     forms the dot product of two vectors.
//     uses unrolled loops for increments equal to one.
//     jack dongarra, linpack, 3/11/78.
//     modified 12/3/93, array(1) declarations changed to array(*)
//     rewritten to C++ 5/5/03
//     n ............ how many components are to be multiplied and summed
//     dx, dy ....... vectors
//     incx, incy ... increments

// arrays: dx, dy
// =============================================================
double
ddotCPP(const int n, double* dx, const int incx, double* dy, const int incy)
{
   double dtemp;
   int i, ix, iy, m, mp1;

   double ddot = 0.0;
   dtemp = 0.0;
   if(n <= 0) return ddot;
   if(!(incx == 1 && incy == 1)){
//
//        code for unequal increments or equal increments
//          not equal to 1
//
      ix = 1;
      iy = 1;
      if(incx < 0) ix = (-n+1)*incx + 1;
      if(incy < 0) iy = (-n+1)*incy + 1;
      for (i = 1; i <= n; i++){
        dtemp += dx[ix-1]*dy[iy-1];
        ix += incx;
        iy += incy;
      }
      ddot = dtemp;
      return ddot;
   }
   else{
//
//        code for both increments equal to 1
//
//
//        clean-up loop
//
      m = n % 5;
      if(m != 0)
         for (i = 1; i <= m; i++)
           dtemp += dx[i-1]*dy[i-1];
      if(n >= 5){
         mp1 = m + 1;
         for (i = mp1; i <= n; i+=5)
           dtemp += dx[i - 1]*dy[i - 1] + dx[i]*dy[i] +
                dx[i + 1]*dy[i + 1] + dx[i + 2]*dy[i + 2] + dx[i + 3]*dy[i + 3];
      }
      ddot = dtemp;
      return ddot;
   }
}

// =============================================================
//
//     constant times a vector plus a vector.
//     uses unrolled loops for increments equal to one.
//     jack dongarra, linpack, 3/11/78.
//     modified 12/3/93, array(1) declarations changed to array(*)
//     rewritten to C++ 5/5/03
//     Input
//          n ....... how many components are to be multiplied and summed
//          da ...... scalar
//          dx ...... vector
//          incx .... increment for the vector dx
//          dy ...... second vector
//          incy .... increment for the second vector
//
//     Output
//          dy
// =============================================================
void
daxpyCPP(const int n, const double da, double* dx, const int incx, double* dy, const int incy)
{
      int i, ix, iy, m, mp1;

      if (n <= 0) return;
      if (da == 0.0) return;
      if(!(incx == 1 && incy ==1)){
//
//        code for unequal increments or equal increments
//          not equal to 1
//
         ix = 1;
         iy = 1;
         if(incx < 0) ix = (-n+1)*incx + 1;
         if(incy < 0) iy = (-n+1)*incy + 1;
         for (i = 1; i <= n; i++){
           dy[iy - 1] += da*dx[ix - 1];
           ix += incx;
           iy += incy;
         }
         return;
      }
      else{
//
//        code for both increments equal to 1
//
//
//        clean-up loop
//
         m = n % 4;
         if(m != 0)
            for (i = 1; i <= m; i++)
               dy[i - 1] += da*dx[i - 1];
         if(n < 4) return;
         mp1 = m + 1;
         for (i = mp1; i<= n; i+=4){
           dy[i - 1] += da*dx[i - 1];
           dy[i] += da*dx[i];
           dy[i + 1] += da*dx[i + 1];
           dy[i + 2] += da*dx[i + 2];
         }
         return;
      }
}


// =============================================================
//
//     scales a vector by a constant.
//     uses unrolled loops for increment equal to one.
//     jack dongarra, linpack, 3/11/78.
//     modified 3/93 to return if incx .le. 0.
//     modified 12/3/93, array(1) declarations changed to array(*)
//     rewritten to C++ 5/5/03
//
// =============================================================
void
dscalCPP(const int n, const double da, double* dx, const int incx)
{
      int i, m, mp1, nincx;

      if(n <= 0 || incx <= 0) return;
      if(incx != 1){
//
//        code for increment not equal to 1
//
         nincx = n*incx;
         for (i = 1; i <= nincx; i+=incx)
            dx[i - 1] *= da;
         return;
      }
      else{
//
//        code for increment equal to 1
//
//
//        clean-up loop
//

         m = n % 5;
         if(m != 0 )
            for (i = 1; i <= m; i++)
              dx[i - 1] *= da;
         if(n < 5) return;
         mp1 = m + 1;
         for (i = mp1; i<= n; i+=5){
           dx[i - 1] *= da;
           dx[i]  *= da;
           dx[i + 1] *= da;
           dx[i + 2] *= da;
           dx[i + 3] *= da;
         }
         return;
      }
}


// =============================================================
//
//  Modified 2002-05-20 for R to add a tolerance of positive definiteness.
//
//
//     dpofa factors a double precision symmetric positive definite
//     matrix.
//
//     dpofa is usually called by dpoco, but it can be called
//     directly with a saving in time if  rcond  is not needed.
//     (time for dpoco) = (1 + 18/n)*(time for dpofa) .
//
//     on entry
//
//        a       double precision(lda, n)
//                the symmetric matrix to be factored.  only the
//                diagonal and upper triangle are used.
//
//        lda     integer
//                the leading dimension of the array  a .
//
//        n       integer
//                the order of the matrix  a .
//
//        eps     double
//                tolerance for the Cholesky decomposition
//
//     on return
//
//        a       an upper triangular matrix  r  so that  a = trans(r)*r
//                where  trans(r)  is the transpose.
//                the strict lower triangle is unaltered.
//                if  info .ne. 0 , the factorization is not complete.
//
//        info    integer
//                = 0  for normal return.
//                = k  signals an error condition.  the leading minor
//                     of order  k  is not positive definite.
//
//     linpack.  this version dated 08/14/78 .
//     cleve moler, university of new mexico, argonne national lab.
//
//     subroutines and functions
//
//     blas ddot
//     fortran dsqrt
//
// a ..... array(lda * n)
//     --> matrix stored according to column major order
//
// =============================================================
void
dpofaCPP(double* a, const int lda, const int n, int* info_dpofa, const double eps)
{
//
//     internal variables
//
      double t;
      double s;
      double absa;
      int j, jm1, k;

//     begin block with ...exits to 40
//
//
      for (j = 1; j <= n; j++){
            *info_dpofa = j;
            s = 0.0;
            jm1 = j - 1;
            if (jm1 >= 1){
               for (k = 1; k <= jm1; k++){
                  t = a[(j-1)*lda+k-1] - ddotCPP(k-1, a + (k-1)*lda, 1, a + (j-1)*lda, 1);
                  t = t / a[(k-1)*lda + k - 1];
                  a[(j - 1)*lda + k - 1] = t;
                  s += t*t;
               }
            }
            s = a[(j - 1)*lda + j - 1] - s;
//     ......exit
//            if (s .le. 0.0d0) go to 40
            absa = (a[(j-1)*lda + j - 1] >= 0) ? a[(j-1)*lda + j - 1] : (-a[(j-1)*lda + j - 1]);
            if (s <= eps * absa){
                return;
            }
            a[(j-1)*lda + j - 1] = sqrt(s);
      }
      *info_dpofa = 0;
      return;
}


// =============================================================
//
//     dpori computes the inverse of the Cholesky factor of a
//     double precision symmetric positive definite matrix
//     using the factors computed by dpofa.
//
//     modification of dpodi by BaT 05/11/95
//
//     on entry
//
//        a       double precision(lda, n)
//                the output  a  from dpofa
//
//        lda     integer
//                the leading dimension of the array  a .
//
//        n       integer
//                the order of the matrix  a .
//
//     on return
//
//        a       if dpofa was used to factor  a  then
//                dpodi produces the upper half of inverse(a) .
//                elements of  a  below the diagonal are unchanged.
//
//     error condition
//
//        a division by zero will occur if the input factor contains
//        a zero on the diagonal and the inverse is requested.
//        it will not occur if the subroutines are called correctly
//        and if dpoco or dpofa has set iinfo .eq. 0 .
//
//     linpack.  this version dated 08/14/78 .
//     cleve moler, university of new mexico, argonne national lab.
//     modified by Berwin A. Turlach 05/11/95
//
//     subroutines and functions
//
//     blas daxpy,dscal
//     fortran mod
//
// =============================================================
void
dporiCPP(double* a, const int lda, const int n)
{
      double t;
      int j, k, kp1;
//
//     compute inverse(r)
//
      for (k = 1; k <= n; k++){
         a[(k-1)*lda + k - 1] = 1.0/a[(k-1)*lda + k - 1];
         t = -a[(k-1)*lda + k - 1];
         dscalCPP(k-1, t, a + (k-1)*lda, 1);
         kp1 = k + 1;
         if (n < kp1) continue;
         for (j = kp1; j <= n; j++){
            t = a[(j-1)*lda + k - 1];
            a[(j-1)*lda + k - 1] = 0.0;
            daxpyCPP(k, t, a + (k-1)*lda, 1, a + (j-1)*lda, 1);
         }
      }
      return;
}


// =============================================================
//
//     dposl solves the double precision symmetric positive definite
//     system a * x = b
//     using the factors computed by dpoco or dpofa.
//
//     on entry
//
//        a       double precision(lda, n)
//                the output from dpoco or dpofa.
//
//        lda     integer
//                the leading dimension of the array  a .
//
//        n       integer
//                the order of the matrix  a .
//
//        b       double precision(n)
//                the right hand side vector.
//
//     on return
//
//        b       the solution vector  x .
//
//     error condition
//
//        a division by zero will occur if the input factor contains
//        a zero on the diagonal.  technically this indicates
//        singularity but it is usually caused by improper subroutine
//        arguments.  it will not occur if the subroutines are called
//        correctly and  iinfo .eq. 0 .
//
//     to compute  inverse(a) * c  where  c  is a matrix
//     with  p  columns
//           call dpoco(a,lda,n,rcond,z,iinfo)
//           if (rcond is too small .or. iinfo .ne. 0) go to ...
//           do 10 j = 1, p
//              call dposl(a,lda,n,c(1,j))
//        10 continue
//
//     linpack.  this version dated 08/14/78 .
//     cleve moler, university of new mexico, argonne national lab.
//
//     subroutines and functions
//
//     blas daxpy,ddot
//
// =============================================================
void
dposlCPP(double* a, const int lda, const int n, double* b)
{
//
//     internal variables
//
      double t;
      int k, kb;
//
//     solve trans(r)*y = b
//
      for (k = 1; k <= n; k++){
         t = ddotCPP(k - 1, a + (k - 1)*lda, 1, b, 1);
         b[k - 1] = (b[k - 1] - t)/a[(k - 1)*lda + k - 1];
      }

//
//     solve r*x = y
//
      for (kb = 1; kb <= n; kb++){
         k = n + 1 - kb;
         b[k - 1] = b[k - 1]/a[(k - 1)*lda + k - 1];
         t = -b[k - 1];
         daxpyCPP(k - 1, t, a + (k-1)*lda, 1, b, 1);
      }
      return;
}
