// Grabbed from FORTRAN routine  in May 2003
// =========================================
//  FORTRAN routine qpgen1 written by Berwin A. Turlach <berwin@alphasun.anu.edu.au> in 1995
//     as a part of 'quadprog' R library.

// I kept several labels in C++ code, sorry...
//  label50, label55, label700, label797, label798, label799

// All matrices must be supplied in arrays where elements are stored in
//   column major order (as in S).

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

// ========================================================================
//  this routine uses the Goldfarb/Idnani algorithm to solve the
//  following minimization problem:
//        minimize  -d^T x + 1/2 *  x^T D x
//        where   A1^T x  = b1
//                A2^T x >= b2
//
//  the matrix D is assumed to be positive definite.  Especially,
//  w.l.o.g. D is assumed to be symmetric.
//
//  Input parameter:
// ------------------
//  dmat   nxn matrix, the matrix D from above (dp)
//         *** WILL BE DESTROYED ON EXIT ***
//         The user has two possibilities:
//         a) Give D (ierr=0), in this case we use routines from LINPACK
//            to decompose D.
//         b) To get the algorithm started we need R^-1, where D=R^TR.
//            So if it is cheaper to calculate R^-1 in another way (D may
//            be a band matrix) then with the general routine, the user
//            may pass R^{-1}.  Indicated by ierr not equal to zero.
//  dvec   nx1 vector, the vector d from above (dp)
//         *** WILL BE DESTROYED ON EXIT ***
//         contains on exit the solution to the initial, i.e.,
//         unconstrained problem
//  fddmat scalar, the leading dimension of the matrix dmat
//  n      the dimension of dmat and dvec (int)
//  amat   lxq matrix (dp)
//         *** ENTRIES CORRESPONDING TO EQUALITY CONSTRAINTS MAY HAVE
//             CHANGED SIGNES ON EXIT ***
//  iamat  (l+1)xq matrix (int)
//         these two matrices store the matrix A in compact form. the format
//         is: [ A=(A1 A2)^T ]
//           iamat(1,i) is the number of non-zero elements in column i of A
//           iamat(k,i) for k>=2, is equal to j if the (k-1)-th non-zero
//                      element in column i of A is A(i,j)
//            amat(k,i) for k>=1, is equal to the k-th non-zero element
//                      in column i of A.
//
//  bvec   qx1 vector, the vector of constants b in the constraints (dp)
//         [ b = (b1^T b2^T)^T ]
//         *** ENTRIES CORRESPONDING TO EQUALITY CONSTRAINTS MAY HAVE
//             CHANGED SIGNES ON EXIT ***
//  fdamat the first dimension of amat as declared in the calling program.
//         fdamat <= n (and iamat must have fdamat+1 as first dimension)
//         if fdamat == 0 then ONLY unconstrained solution is reported
//  q      integer, the number of constraints.
//  meq    integer, the number of equality constraints, 0 <= meq <= q.
//  ierr   integer, code for the status of the matrix D:
//            ierr =  0, we have to decompose D
//            ierr != 0, D is already decomposed into D=R^TR and we were
//                       given R^{-1}.
//  eps    tolerance for Cholesky decomposition (argument added by AK)
//
//  Output parameter:
// -------------------
//  sol   nx1 the final solution (x in the notation above)
//  crval scalar, the value of the criterion at the minimum
//  iact  qx1 vector, the constraints which are active in the final
//        fit (int)
//  nact  scalar, the number of constraints active in the final fit (int)
//  iter  2x1 vector, first component gives the number of "main"
//        iterations, the second one says how many constraints were
//        deleted after they became active
//  ierr  integer, error code on exit, if
//           ierr = 0, no problems
//           ierr = 1, the minimization problem has no solution
//           ierr = 2, problems with decomposing D, in this case sol
//                     contains garbage!!
//
//  Working space:
// ----------------
//  work  vector with length at least 2*n+r*(r+5)/2 + 2*q +1
//        where r=min(n,q)
//
#include "solve.QP.compact.h"

extern "C"{

using namespace std;

void
qpgen1CPP (double* dmat, double* dvec, int* fddmat, int* n,
           double* sol, double* crval,
           double* amat, int* iamat, double* bvec,
           int* fdamat, int* q, int* meq,
           int* iact, int* nact, int* iter, double* work, int* ierr,
           double* eps)
{
   int i, j, l, l1,
       it1, iwzv, iwrv, iwrm, iwsv, iwuv, nvl,
       r, iwnbv;
   int * iinfoqp = new int;
   double temp, sum, t1, tt, gc, gs, nu;
   bool t1inf, t2min;

   r = (*n < *q) ? *n : *q;
   l = 2*(*n) + (r*(r+5))/2 + 2*(*q) + 1;

//
// store the initial dvec to calculate below the unconstrained minima of
// the critical value.
//
   for (i = 1; i <= *n; i++) work[i-1] = dvec[i-1];
   for (i = *n+1; i <= l; i++) work[i-1] = 0.0;
   for (i = 1; i <= *q; i++) iact[i-1] = 0;

//
// get the initial solution
//
   if (*ierr == 0){
      dpofaCPP(dmat, *fddmat, *n, iinfoqp, *eps);
      if(*iinfoqp != 0){
         *ierr = 2;
         delete iinfoqp;
         return;
      }
      delete iinfoqp;
      dposlCPP(dmat, *fddmat, *n, dvec);
      dporiCPP(dmat, *fddmat, *n);
   }
   else{   // else 118
//
// Matrix D is already factorized, so we have to multiply d first with
// R^-T and then with R^-1.  R^-1 is stored in the upper half of the
// array dmat.
//
      delete iinfoqp;
      for (j = 1; j <= *n; j++){
         sol[j-1] = 0.0;
         for (i = 1; i <= j; i++)
             sol[j-1] += dmat[(j-1)*(*n) + i-1]*dvec[i-1];
      }    // 20
      for (j = 1; j <= *n; j++){
         dvec[j-1] = 0.0;
         for (i = j; i <= *n; i++)
             dvec[j-1] += dmat[(i-1)*(*n) + j-1]*sol[i-1];
      }   // 22
   }    // end of else 118

//
// set lower triangular of dmat to zero, store dvec in sol and
// calculate value of the criterion at unconstrained minima
//
   *crval = 0.0;
   for (j = 1; j <= *n; j++){
       sol[j-1] = dvec[j-1];
       *crval += work[j-1]*sol[j-1];
       work[j-1] = 0.0;
       for (i = j + 1; i <= *n; i++)
           dmat[(j-1)*(*n) + i-1] = 0.0;
   }    // 30
   *crval = -(*crval)/2.0;
   *ierr = 0;

//
// leave if only unconstrained solution is required
//
   if (*fdamat == 0) return;

//
// calculate some constants, i.e., from which index on the different
// quantities are stored in the work matrix
//
      iwzv  = *n;
      iwrv  = iwzv + (*n);
      iwuv  = iwrv + r;
      iwrm  = iwuv + r + 1;
      iwsv  = iwrm + (r*(r+1))/2;
      iwnbv = iwsv + *q;

//
// calculate the norm of each column of the A matrix
//
// AK: work[iwnbv : iwnbv + (*q) - 1] then contains norms of
//       each column of the matrix A
   for (i = 1; i <= *q; i++){
       sum = 0.0;
       for (j = 1; j <= iamat[(i-1)*(*fdamat+1)]; j++)
           sum += amat[(i-1)*(*fdamat) + j-1]*amat[(i-1)*(*fdamat) + j-1];
       work[iwnbv + i - 1] = sqrt(sum);
   }    // 51
   *nact = 0;
   iter[0] = 0;
   iter[1] = 0;
   label50:

//
// start a new iteration
//
   (iter[0])++;

//
// calculate all constraints and check which are still violated
// for the equality constraints we have to check whether the normal
// vector has to be negated (as well as bvec in that case)
//
// AK: work[iwsv : iwsv+(*q)-1] then contains
//       a recent value of the constraint minus right side
//       for inequality const., if this is non-negative then
//       the constraint is satisfied
   l = iwsv;
   for (i = 1; i <= *q; i++){
       l++;
       sum = -bvec[i-1];
       for (j = 1; j <= iamat[(i-1)*(*fdamat+1)]; j++)
           sum += amat[(i-1)*(*fdamat) + j-1]*sol[iamat[(i-1)*(*fdamat+1)+j]-1];
       if (i > *meq)
          work[l-1] = sum;
       else{
          work[l-1] = -abs(sum);
          if (sum > 0.0){
              for (j = 1; j <= iamat[(i-1)*(*fdamat+1)]; j++)
                  amat[(i-1)*(*fdamat)+j-1] = -amat[(i-1)*(*fdamat)+j-1];
              bvec[i-1] = -bvec[i-1];
          }
       }
   }    // 60

//
// as safeguard against rounding errors set already active constraints
// explicitly to zero
//
   for (i = 1; i <= *nact; i++)
       work[iwsv + iact[i-1] - 1] = 0.0;

//
// we weight each violation by the number of non-zero elements in the
// corresponding row of A. then we choose the violated constraint which
// has maximal absolute value, i.e., the minimum.
// by obvious commenting and uncommenting we can choose the strategy to
//
// take always the first constraint which is violated. ;-)
//
   nvl = 0;
   temp = 0.0;
   for (i = 1; i <= *q; i++){
       if (work[iwsv+i-1] < temp*work[iwnbv+i-1]){
          nvl = i;
          temp = work[iwsv+i-1]/work[iwnbv+i-1];
       }
//        if (work[iwsv+i-1] < 0.0){
//           nvl = i;
//           break;
//        }
   }      // 71

   if (nvl == 0) return;

//
// calculate d=J^Tn^+ where n^+ is the normal vector of the violated
// constraint. J is stored in dmat in this implementation!!
// if we drop a constraint, we have to jump back here.
//
   label55:
   for (i = 1; i <= *n; i++){
       sum = 0.0;
       for (j = 1; j <= iamat[(nvl-1)*(*fdamat + 1)]; j++)
             sum += dmat[(i-1)*(*n)+iamat[(nvl-1)*(*fdamat+1)+j]-1] * amat[(nvl-1)*(*fdamat)+j-1];
       work[i-1] = sum;
   }

//
// Now calculate z = J_2 d_2
//
   l1 = iwzv;
   for (i = 1; i <= *n; i++)
       work[l1+i-1] = 0.0;
   for (j = *nact+1; j <= *n; j++)
       for (i = 1; i <= *n; i++)
            work[l1+i-1] += dmat[(j-1)*(*n)+i-1] * work[j-1];

//
// and r = R^{-1} d_1, check also if r has positive elements (among the
// entries corresponding to inequalities constraints).
//
   t1inf = true;
   it1 = *nact;
   for (i = *nact; i >= 1; i--){
       sum = work[i-1];
       l  = iwrm + (i*(i+3))/2;
       l1 = l - i;
       for (j = i + 1; j <= *nact; j++){
          sum -= work[l-1] * work[iwrm+j-1];
       }   // 96
       sum = sum/work[l1-1];
       if (iact[i-1] <= *meq) break;
       if (sum <= 0.0) break;
       t1inf = false;
       it1 = i;
   }   // 95

//
// if r has positive elements, find the partial step length t1, which is
// the maximum step in dual space without violating dual feasibility.
// it1  stores in which component t1, the min of u/r, occurs.
//
   t1 = 0.0;
   if (!t1inf){
         t1   = work[iwuv+it1-1]/work[iwrv+it1-1];
         for (i = 1; i <= *nact; i++){
            if (iact[i-1] <=  *meq) break;
            if (work[iwrv+i-1] <= 0.0) break;
            temp = work[iwuv+i-1]/work[iwrv+i-1];
            if (temp < t1){
               t1   = temp;
               it1  = i;
            }
         }  // 100
   }

//
// test if the z vector is equal to zero
//
   sum = 0.0;
   for (i = iwzv+1; i <= iwzv+(*n); i++)
     sum += work[i-1]*work[i-1];
   temp = 1000.0;
   sum  += temp;

   if (temp == sum){
//
// No step in pmrimal space such that the new constraint becomes
// feasible. Take step in dual space and drop a constant.
//
         if (t1inf){
//
// No step in dual space possible either, problem is not solvable
//
            *ierr = 1;
            return;
         }
         else{
//
// we take a partial step in dual space and drop constraint it1,
// that is, we drop the it1-th active constraint.
// then we continue at step 2(a) (marked by label 55)
//
            for (i = 1; i <= *nact; i++)
               work[iwuv+i-1] = work[iwuv+i-1] - t1*work[iwrv+i-1];
            work[iwuv+(*nact)+1-1] = work[iwuv+(*nact)+1-1] + t1;
            goto label700;
         }
   }
   else{   // else 317

//
// compute full step length t2, minimum step in primal space such that/
// the constraint becomes feasible.
// keep sum (which is z^Tn^+) to update crval below!
//
         sum = 0.0;
         for (i = 1; i <= iamat[(nvl-1)*(*fdamat+1)]; i++)
            sum += work[iwzv+iamat[(nvl-1)*(*fdamat+1)+i]-1]*amat[(nvl-1)*(*fdamat)+i-1];
         tt = -work[iwsv+nvl-1]/sum;
         t2min = true;
         if (!t1inf){
            if (t1 < tt){
               tt    = t1;
               t2min = false;
            }
         }

//
// take step in primal and dual space
//
         for (i = 1; i <= *n; i++)
            sol[i-1] += tt*work[iwzv+i-1];
         *crval += tt*sum*(tt/2.0 + work[iwuv+(*nact)+1-1]);
         for (i = 1; i <= *nact; i++)
            work[iwuv+i-1] -= tt*work[iwrv+i-1];
         work[iwuv+(*nact)+1-1] += tt;

//
// if it was a full step, then we check wheter further constraints are
// violated otherwise we can drop the current constraint and iterate once
// more
         if(t2min){   //   if 348
//
// we took a full step. Thus add constraint nvl to the list of active
// constraints and update J and R
//
            (*nact)++;
            iact[*nact-1] = nvl;

//
// to update R we have to put the first nact-1 components of the d vector
// into column (nact) of R
//
            l = iwrm + ((*nact-1)*(*nact))/2 + 1;
            for (i = 1; i <= *nact-1; i++){
               work[l-1] = work[i-1];
               l++;
            }
//
// if now nact=n, then we just have to add the last element to the new
// row of R.
// Otherwise we use Givens transformations to turn the vector d(nact:n)
// into a multiple of the first unit vector. That multiple goes into the
// last element of the new row of R and J is accordingly updated by the
// Givens transformations.
//
            if (*nact == *n)
               work[l-1] = work[*n-1];
            else{   // else 374
               for (i = *n; i >= *nact+1; i--){   // do 160
//
// we have to find the Givens rotation which will reduce the element
// (l1) of d to zero.
// if it is already zero we don't have to do anything, except of
// decreasing l1
//
                    if (work[i-1] == 0.0) break;  // leave the for loop
                    gc = (abs(work[i-2]) > abs(work[i-1])) ? abs(work[i-2]) : abs(work[i-1]);
                    gs = (abs(work[i-2]) <= abs(work[i-1])) ? abs(work[i-2]) : abs(work[i-1]);
                    temp = (work[i-2] < 0) ? -gc*sqrt(1+gs*gs/(gc*gc)) : gc*sqrt(1+gs*gs/(gc*gc));
                    gc   = work[i-2]/temp;
                    gs   = work[i-1]/temp;
//
// The Givens rotation is done with the matrix (gc gs, gs -gc).
// If gc is one, then element (i) of d is zero compared with element
// (l1-1). Hence we don't have to do anything.
// If gc is zero, then we just have to switch column (i) and column (i-1)
// of J. Since we only switch columns in J, we have to be careful how we
// update d depending on the sign of gs.
// Otherwise we have to apply the Givens rotation to these columns.
// The i-1 element of d has to be updated to temp.
//
                    if (gc == 1.0) break;   // leave the for loop
                    if (gc == 0.0){
                       work[i-2] = gs * temp;
                       for (j = 1; j <= *n; j++){
                          temp                   = dmat[(i-2)*(*n) + j-1];
                          dmat[(i-2)*(*n) + j-1] = dmat[(i-1)*(*n) + j-1];
                          dmat[(i-1)*(*n) + j-1] = temp;
                       }
                    }
                    else{
                       work[i-2] = temp;
                       nu = gs/(1.0+gc);
                       for (j = 1; j <= *n; j++){
                          temp                   = gc*dmat[(i-2)*(*n) + j-1] + gs*dmat[(i-1)*(*n) + j-1];
                          dmat[(i-1)*(*n) + j-1] = nu*(dmat[(i-2)*(*n) + j-1]+temp) - dmat[(i-1)*(*n) + j-1];
                          dmat[(i-2)*(*n) + j-1] = temp;
                       }
                    }
               }     // 160
//
// l is still pointing to element (nact,nact) of the matrix R.
// So store d(nact) in R(nact,nact)
               work[l-1] = work[*nact-1];
            }   // end of else 374
         }       // end of if 348
         else{   // else from if 348

//
// we took a partial step in dual space. Thus drop constraint it1,
// that is, we drop the it1-th active constraint.
// then we continue at step 2(a) (marked by label 55)
// but since the fit changed, we have to recalculate now "how much"
// the fit violates the chosen constraint now.
//
            sum = -bvec[nvl-1];
            for (j = 1; j <= iamat[(nvl-1)*(*fdamat+1)]; j++)
               sum += sol[iamat[(nvl-1)*(*fdamat+1)+j]-1]*amat[(nvl-1)*(*fdamat)+j-1];
            if(nvl > *meq)
               work[iwsv+nvl-1] = sum;
            else{
               work[iwsv+nvl-1] = -abs(sum);
               if(sum > 0.0){
                  for (j = 1; j <= iamat[(nvl-1)*(*fdamat+1)]; j++)
                     amat[(nvl-1)*(*fdamat)+j-1] = -amat[(nvl-1)*(*fdamat)+j-1];
                  bvec[i-1] = -bvec[i-1];
               }
            }
            goto label700;
         }       // end of else from if 348
   }   // end of else 317
   goto label50;    // go to the next iteration

//
// Drop constraint it1
//
   label700:
//
// if it1 = nact it is only necessary to update the vector u and nact
//
   if (it1 == *nact) goto label799;
//
// After updating one row of R (column of J) we will also come back here
//

   label797:
//
// we have to find the Givens rotation which will reduce the element
// (it1+1,it1+1) of R to zero.
// if it is already zero we don't have to do anything except of updating
// u, iact, and shifting column (it1+1) of R to column (it1)
// l  will point to element (1,it1+1) of R
// l1 will point to element (it1+1,it1+1) of R
//
   l  = iwrm + (it1*(it1+1))/2 + 1;
   l1 = l+it1;
   if (work[l1-1] == 0.0) goto label798;
   gc = (abs(work[l1-2]) > abs(work[l1-1])) ? abs(work[l1-2]) : abs(work[l1-1]);
   gs = (abs(work[l1-2]) <= abs(work[l1-1])) ? abs(work[l1-2]) : abs(work[l1-1]);
   temp = (work[l1-2] < 0) ? -gc*sqrt(1+gs*gs/(gc*gc)) : gc*sqrt(1+gs*gs/(gc*gc));
   gc   = work[l1-2]/temp;
   gs   = work[l1-1]/temp;

//
// The Givens rotatin is done with the matrix (gc gs, gs -gc).
// If gc is one, then element (it1+1,it1+1) of R is zero compared with
// element (it1,it1+1). Hence we don't have to do anything.
// if gc is zero, then we just have to switch row (it1) and row (it1+1)
// of R and column (it1) and column (it1+1) of J. Since we swithc rows in
// R and columns in J, we can ignore the sign of gs.
// Otherwise we have to apply the Givens rotation to these rows/columns.
//
   if (gc == 1.0) goto label798;
   if (gc == 0.0){
      for (i = it1+1; i <= *nact; i++){
         temp       = work[l1-2];
         work[l1-2] = work[l1-1];
         work[l1-1]   = temp;
         l1 = l1+i;
      }
      for (i = 1; i <= *n; i++){
         temp          = dmat[(it1-1)*(*n)+i-1];
         dmat[(it1-1)*(*n)+i-1] = dmat[it1*(*n)+i-1];
         dmat[it1*(*n)+i-1] = temp;
      }
   }
   else{
      nu = gs/(1.0+gc);
      for (i = it1+1; i <= *nact; i++){
         temp         = gc*work[l1-2] + gs*work[l1-1];
         work[l1-1]   = nu*(work[l1-2]+temp) - work[l1-1];
         work[l1-2]   = temp;
         l1 = l1+i;
      }
      for (i = 1; i <= *n; i++){
         temp                   = gc*dmat[(it1-1)*(*n)+i-1] + gs*dmat[it1*(*n)+i-1];
         dmat[it1*(*n)+i-1]     = nu*(dmat[(it1-1)*(*n)+i-1]+temp) - dmat[it1*(*n)+i-1];
         dmat[(it1-1)*(*n)+i-1] = temp;
      }
   }

//
// shift column (it1+1) of R to column (it1) (that is, the first it1
// elements). The posit1on of element (1,it1+1) of R was calculated above
// and stored in l.
//
   label798:
   l1 = l-it1;
   for (i = 1; i <= it1; i++){
      work[l1-1]=work[l-1];
      l++;
      l1++;
   }
//
// update vector u and iact as necessary
// Continue with updating the matrices J and R
//
   work[iwuv+it1-1] = work[iwuv+it1];
   iact[it1-1]      = iact[it1];
   it1++;
   if (it1 < *nact) goto label797;
   label799:
   work[iwuv+*nact-1]   = work[iwuv+*nact];
   work[iwuv+*nact]     = 0.0;
   iact[*nact-1]        = 0;
   (*nact)--;
   (iter[1])++;
   goto label55;

   return;
}

}    // end of extern "C"




