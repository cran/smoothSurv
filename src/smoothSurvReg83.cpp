// Main C++ function to perform estimation in linear model
//  with possibly censored observations
//  with smoothed error distribution
// ================================================================
// version8.3: constrained maximization of penalized log-likleihood
//             with some improvements compared to version 8.2

//     * using 'a' parametrization

//     * penalty on a's
//       ==> it appeared to be a better choice than a penalty on c's

//     * equality constraints
//                   E(epsilon) = 0
//                   var(epsilon) = 1
//       are explicitely solved
//       ==> only g-3 a's are estimated
//           two a's are a function of the rest and one a is equal to zero
//       HOWEVER: one has to be careful since each of these two a's
//                is a logarithm of some expression and this expression
//                can be generally zero or even negative!!!
//       HOPE: expression will not be non-positive if the two a's are
//             chosen in the area where there are data, i.e. if
//             the two a's does not correspond to zero c's

//     * in comments etc. I denote the (g-3) a's which are really estimated
//       as d_1, ..., d_{g-3}

//     * lambda is multiplied by n (the sample size)
//
// ALGORITHM:
//     * conditional maximization
//         ==> I believe it's more stable than full maximization
//     * with conditional steps (regression -> a's -> regression -> a's etc.)
//     * step-halving to find a step-size
//     * at the end of the conditional maximization I do at least one more step
//       of the full maximization
//          --> a) to be sure I converge
//          --> b) to compute full Hessian matrix

// =====================================================================================

// start working on this version: July 3, 2003

// -------------------------------------------
//

// REMARKS:
//   * this version assumes that there is at least intercept
//     included in the model

//   * either scale or 'a' coefficients can be fixed and not estimated

//   * error distribution assumes epsilon ~ sum_{j=1}^g c_j(a) N(mu_j, sigma0^2)
//     where c_j = exp(a_j)/sum(exp(a_l))

//   * model: log(T) = Y = beta0 + beta'x + sigma*epsilon

//   * if a's (c's) are estimated there must be at least 4 knots


//   nSplinesP     - on INPUT: number of G-spline at the beginnig
//                   on OUTPUT: number of G-splines after removal of knots with zero G-spline coeff.
//
//   knotsP        - on both INPUT and OUTPUT: original knots
//
//   knotRemovedP  - on OUTPUT: vector of zeros and ones indicating which knots have been removed

//
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

#include "smoothSurvReg83.h"

extern "C"{

using namespace std;
using namespace SCYTHE;

// GLOBAL VARIABLES
// *****************
   int n,                        // number of observations
       nTheta,                   // number of all parameters to be really estimated
       nThetaSq,                 // nTheta^2
       nBeta,                    // number of columns in X matrix (number of beta parameters to be estimated)
       nScale,                   // number of scale parameters to be estimated (either 1 or 0)
       nRegres,                  // nBeta + nScale
       nD,                       // number of a coefficients to be estimated (either g-3 or 0)
       nSplines,                 // number of splines (after removal of some c's it's smaller)
       debug;

   double lambda,                // tuning parameter for the penalty term
          sigmaZero,             // standard deviation of one basis "spline"
          invsigmaZero,          // inversion of sigmaZero
          logsigmaZero;          // log(sigmaZero)

   bool estScale,                // true if scale is to be estimated
        estA,                    // true if a's are to be estimated
        info;                    // do I want to print some information during the iteration process?

   double correctLik;             // correction to likelihood due to log transformation of the response

// 'double' variables modified by other procedures
   double logLikelihood,          // corrected value of unpenalized log-likelihood
          penalty;                // penalty term

// Response and covariates matrices
   Matrix<int> statusMat;           // status vector (1=exact, 0=right censored, 3=interval censored)
   Matrix<double> resp1Mat,         // response (exact, right censored or lower limit)  n x 1
                  resp2Mat,         // upper limit of response for interval censored data  n x 1
                  XMat,             // design matrix (covariates)  n x nBeta
                  offset;           // offset vector n x 1

// Matrices and vectors used in Newton-Raphson steps
//   matrices UMat, HMat, IMat, Ccoef are modified by other procedures
   Matrix<double> UMat,          // score vector
                  UMatUP,        // zero vector ((nBeta+nScale) x 1)
                  HMat,          // minus Hessian matrix (minus second derivative) for the Lagrange function
                  IMat,          // minus Hessian matrix for UNpenalized log-likelihood
                  GMat,          // minus Hessian matrix for penalty term
                  GMatLEFT,      // building block for GMat
                  GMatUP,        // building block for GMat
                  UMatRegres,    // score w.r.t. beta and log(scale) only
                  HMatRegres,    // minus second deriv. w.r.t. beta and log(scale) only
                  UMatD,         // score w.r.t. (c_1, ..., c_{g-3})' only
                  HMatD;         // minus second deriv. w.r.t. (c_1, ..., c_{g-3})' only

// A and C coefficients
   Matrix<double> Acoef,         // (a_1, ..., a_g)'
                  Ccoef;         // (c_1, ..., c_g)'


// Lagrangians (for optimization with equality constraints)
   Matrix<double> XiMat;         // vector 2 x 1

// Equality constraints and their derivatives
   Matrix<double> minCon;           // minus constraints
   Matrix<double> dConInRows;
   Matrix<double> dCon0, dCon1;
   Matrix<double> ddCon0, ddCon1;

// Splines related vectors and matrices
   Matrix<double> knotsMat,       // knots
                  knotsMatSq,     // squared knots
                  OmegaExpA2,     // matrix with d(exp(A1, A2))/d(exp(D)) (g-3) x 2
                  tOmegaExpA2,    // its transposition                    2 x (g-3)
                  OmegaZeroExpA2, // intercept to compute the two a's     2 x 1
                         // ==> a[the two] = tOmegaExpA2 * d + OmegaZeroExpA2
                  minLamDDMat,    // matrix -lambda*D'*D (useful for derivatives) (matrix g x g)
                  LamDDg1,        // matrix (g - 1) x (g - 1)
                  minLamDD2,      // matrix 2 x 2
                  minLamDD2g3,    // matrix 2 x (g - 3)           blocks of the previous matrix
                  minLamDDg3;     // matrix (g - 3) x (g - 3)

// Various derivatives
   Matrix<double> * ddCdAA;
   Matrix<double> * ddCdDD;
   Matrix<double> dA2dD,
                  ddAdDD1,
                  ddAdDD2,
                  dCdA2,
                  dCdAg3,
                  dCdD;

// Matrices with possibly fixed some parameters
   Matrix<double> GammaGlobal;

// Matrix indicating which three c's are expressed as a function of the remaining ones
   Matrix<int> lastThreeA,     // indeces of the zero a and of the two a's
               restA,          // indeces of (g - 3) a's which are real parameters
               lastTwoAshort,  // indeces of the two a's in the sequence of a's with removed zero
               restAshort;     // indeces of (g - 3) a's in the sequence of a's with removed zero

// Help matrices
   Matrix<double> ZeroMat(1, 1, true, 0.0);       // matrix 1 x 1 with 0
   Matrix<double> OneMat(1, 1, true, 1.0);        // matrix 1 x 1 with 1


// MAIN FUNCTION
// **************
void
smoothSurvReg83(int *nP,   int *nyP,   int *nBetaP,   int *nSplinesP,
               double *XP,   double *YP,   double *offsetP,
               double *knotsP,   double *sigmaZeroP,
               int *lastThreeAP,
               int *estScaleP,   int *estAP,
               double *BetaP,    double *GammaP,   double *AcoefP,   double *CcoefP,
               double *penalloglikP,
               double *loglikP,   double *correctLikP,   double *penaltyP,
               double *HP,   double *IP,   double *GP,
               double *UP,   double *dCdDP,
               double *HaP,   double *IaP,   double *GaP,   double *dConP,
               double *lambdaP,   int *difforderP,
               int *maxiterP,   int *firstiterP,
               double *epsP,   double *tolCholP,   double *tolEigenP,
               int *maxhalfP,
               int *infoP,   int *debugP,
               int *failP,   int *nonPosDefHP)
{
// LONG LASTING VARIABLES
// ======================
     int i, j, k, kk;              // helping integers
     Matrix<double> ThetaMat,
                    newThetaMat;


// ASSIGN PROPER VALUES TO GLOBAL VARIABLES
// =========================================
     n = *nP;
     nBeta = *nBetaP;
     nSplines = *nSplinesP;
     estScale = (*estScaleP != 0);
     estA = (*estAP != 0);
     estA = (nSplines < 4) ? false : estA;      // do not estimate a's if there are too few a's
     nScale = estScale ? 1 : 0;
     nRegres = nBeta + nScale;
     if (estA) nD = (*maxiterP > 0) ? (nSplines - 1) : (nSplines - 3);    // I start with optimization with two eq. constraints
     else      nD = 0;
     nTheta = nBeta + nScale + nD;
     nThetaSq = nTheta * nTheta;

     debug = *debugP;
     info = (*infoP != 0);

     lambda = n*(*lambdaP);
     sigmaZero = *sigmaZeroP;
     invsigmaZero = (sigmaZero > 0) ? (1/sigmaZero) : FLT_MAX;
     logsigmaZero = (sigmaZero > 0) ? log(sigmaZero) : (-FLT_MAX);

  // Status vector and response vector/matrix.
     int *statusarray = new int[n];
     double *resp1array = new double[n];
     double *resp2array = new double[n];
     if(*nyP == 2){                               // only exact/right censored observations
       for(i = 0; i < n; i++){
          statusarray[i] = int(YP[n + i]);
          resp1array[i] = YP[i];
       }
       statusMat = Matrix<int>(n, 1, statusarray);                 // vector (n x 1) with status (0 = right censored, 1 = event)
       resp1Mat = Matrix<double>(n, 1, resp1array);                // vector (n x 1) with response (exact or right censored)
     }
     else{       // ny is 3                     // also interval censored observations
       for(i = 0; i < n; i++){
          statusarray[i] = int(YP[2*n + i]);
          resp1array[i] = YP[i];
          resp2array[i] = YP[n + i];
       }
       statusMat = Matrix<int>(n, 1, statusarray);                // vector (n x 1) with status (0 = right censored, 1 = event, 3 = interval censored).
       resp1Mat = Matrix<double>(n, 1, resp1array);               // vector (n x 1) with response (exact or censored or lower limit of intervals)
       resp2Mat = Matrix<double>(n, 1, resp2array);               // vector (n x 1) with upper limits of interval, arbitrary value for exact/right c. cases
     }
     delete [] statusarray;
     delete [] resp1array;
     delete [] resp2array;

  // Offset vector
     offset = Matrix<double>(n, 1, offsetP);

  // X matrix with covariates.
     XMat = Matrix<double>(nBeta, n, XP);
     XMat = t(XMat);                            // transposition since Scythe uses row-major order when filling the matrix

  // Correction to log-likelihood
     correctLik = *correctLikP;

  // Knots and related matrices
     knotsMat = Matrix<double>(nSplines, 1, knotsP);  // vector with knots (expectations of narrow normals)
     knotsMatSq = knotsMat & knotsMat;  // vector with squared knots

  // Indicators of the a coefficients expressed as a function of the rest
     lastThreeA = Matrix<int>(3, 1, true, 0);
     lastTwoAshort = Matrix<int>(2, 1, true, 0);
     if (estA){
        lastThreeA = Matrix<int>(3, 1, lastThreeAP);
        restA = Matrix<int>(nD, 1, false);

        if (*maxiterP > 0){
     // if maxiter > 0, restA ........... indeces of all non-zero a's
     //                 lastThreeA[0] ... index of the zero a
           for (i = 0, k = 0; k < nD; i++){
               if (i != lastThreeA[0]){
                   restA[k] = i;
                   k++;
               }
           }
        }

        else{
     // if *maxiterP == 0 I compute only derivatives and these with 'd' parametrization
     //   in that case, nD has a value nSplines - 3
           lastTwoAshort[0] = (lastThreeA[0] < lastThreeA[1] ? lastThreeA[1] - 1 : lastThreeA[1]);
           lastTwoAshort[1] = (lastThreeA[0] < lastThreeA[2] ? lastThreeA[2] - 1 : lastThreeA[2]);

           for (i = 0, k = 0; k < nD; i++){
               if (i != lastThreeA[0] && i != lastThreeA[1] && i != lastThreeA[2]){
                   restA[k] = i;
                   k++;
               }
           }

           restAshort = Matrix<int>(nD, 1, false);
           for (i = 0, k = 0; k < nD; i++){
               if (i != lastTwoAshort[0] && i != lastTwoAshort[1]){
                   restAshort[k] = i;
                   k++;
               }
           }
        }
     }


  // D matrix to compute ordered differences, DD matrix to compute squared ordered diff
  //   and matrices related to penalized log-likelihood
     Matrix<double> DMat, DDMat;
     if (estA){
         DMat = D_operator(nSplines, *difforderP);             // matrix (g-k) x g
         DDMat = t(DMat) * DMat;
         minLamDDMat = -lambda * DDMat;                        // matrix g x g

     // minus block of minLamDDMat corresponding to non-zero a's
         if (*maxiterP > 0){   // nD is nSplines - 1
            LamDDg1 = Matrix<double>(nSplines - 1, nSplines - 1, false);
            for (i = 0; i < nSplines - 1; i++)
                for (j = 0; j < nSplines - 1; j++)
                    LamDDg1(i, j) = -minLamDDMat(restA[i], restA[j]);
         }

         else{                // nD is nSplines - 3
            OmegaExpA2 = Matrix<double>(nSplines - 2, 2, false);
            if (deriv_expAD(knotsMat, sigmaZero, lastThreeA, OmegaExpA2, false) != 0){
               if (info) REprintf("      STOP: Untractable choice of the reference knots.\n");
               *failP = 101;
               return;
            }
            OmegaZeroExpA2 = t(OmegaExpA2(0, 0, 0, 1));
            OmegaExpA2 = OmegaExpA2(1, 0, nD, 1);
            tOmegaExpA2 = t(OmegaExpA2);

           // block of minLamDDMat corresponding to the two special a's
            minLamDD2 = Matrix<double>(2, 2, false);
            minLamDD2(0, 0) = minLamDDMat(lastThreeA[1], lastThreeA[1]);
            minLamDD2(0, 1) = minLamDDMat(lastThreeA[1], lastThreeA[2]);
            minLamDD2(1, 0) = minLamDDMat(lastThreeA[2], lastThreeA[1]);
            minLamDD2(1, 1) = minLamDDMat(lastThreeA[2], lastThreeA[2]);

           // block of minLamDDMat corresponding to the cross-product of two a's and d's
            minLamDD2g3 = Matrix<double>(2, nD, false);
            for (i = 0; i < 2; i++)
               for (j = 0; j < nD; j++)
                   minLamDD2g3(i, j) = minLamDDMat(lastThreeA[i + 1], restA[j]);

           // block of minLamDDMat corresponding to d's
            minLamDDg3 = Matrix<double>(nD, nD, false);
            for (i = 0; i < nD; i++)
               for (j = 0; j < nD; j++)
                   minLamDDg3(i, j) = minLamDDMat(restA[i], restA[j]);
        }
     }


  // Building blocks for score vector (null part for unpenalized parameters)
  // and for GMat (also null for unpenalized parameters)
     UMatUP = Matrix<double>(nRegres, 1, true, 0.0);
     if (estA){
        GMatLEFT = Matrix<double>(nTheta, nRegres, true, 0.0);
        GMatUP = Matrix<double>(nRegres, nD, true, 0.0);
     }


// INITIALIZE remaining GLOBAL VARIABLES (usually BY ZEROS)
// ========================================================
  // it is useful for returning zero matrices if something is not estimated
     penalty = 0;
     logLikelihood = 0;
     UMat = Matrix<double>(nTheta, 1, true, 0.0);
     HMat = Matrix<double>(nTheta, nTheta, true, 0.0);
     GMat = Matrix<double>(nTheta, nTheta, true, 0.0);
     IMat = Matrix<double>(nTheta, nTheta, true, 0.0);
     if (*maxiterP > 0){
        UMatRegres = Matrix<double>(nRegres, 1, true, 0.0);
        HMatRegres = Matrix<double>(nRegres, nRegres, true, 0.0);
        if (estA){
           UMatD = Matrix<double>(nD, 1, true, 0.0);
           HMatD = Matrix<double>(nD, nD, true, 0.0);
           minCon = Matrix<double>(2, 1, true, 0.0);
           dConInRows = Matrix<double>(2, nSplines - 1, true, 0.0);
           dCon0 = Matrix<double>(nSplines - 1, 1, true, 0.0);
           dCon1 = Matrix<double>(nSplines - 1, 1, true, 0.0);
           ddCon0 = Matrix<double>(nSplines - 1, nSplines - 1, true, 0.0);
           ddCon1 = Matrix<double>(nSplines - 1, nSplines - 1, true, 0.0);
           XiMat = Matrix<double>(2, 1, true, 0.0);
        }
     }


// INITIAL ESTIMATES and POSSIBLY ASSIGN VALUES TO ADDITIONAL GLOBAL VARIABLES
// ===========================================================================
     ThetaMat = Matrix<double>(nTheta, 1, true, 0.0);
     for (i = 0; i < nBeta; i++) ThetaMat[i] = BetaP[i];
     if (estScale) ThetaMat[i] = *GammaP;
     else          GammaGlobal = Matrix<double>(1, 1, GammaP);
     if (estA){
         j = i + nScale;
         if (*maxiterP > 0){
            for (k = 0; k < nSplines; k++){
                if (k == lastThreeA[0]) continue;
                ThetaMat[j] = AcoefP[k];
                j++;
            }
         }
         else{    // no iterations => only derivatives w.r.t. g - 3 a's
            for (k = 0; k < nSplines; k++){
                if (k == lastThreeA[0] || k == lastThreeA[1] || k == lastThreeA[2]) continue;
                ThetaMat[j] = AcoefP[k];
                j++;
            }
         }
     }
     Acoef = Matrix<double>(nSplines, 1, AcoefP);
     Ccoef = Matrix<double>(nSplines, 1, false);
     if (A_to_C(Acoef, Ccoef) != 0){
         if (info) REprintf("      STOP: Untractable initial 'a' coefficients.\n");
         *failP = 99;
         return;
     }
     newThetaMat = ThetaMat;


// INITIALIZE GLOBAL DERIVATIVE MATRICES
// ======================================
     ddCdDD = new Matrix<double> [nSplines];
     ddCdAA = new Matrix<double> [nSplines];
     if (estA){
        if (*maxiterP > 0){
           for (i = 0; i < nSplines; i++) ddCdDD[i] = Matrix<double>(nD, nD, true, 0.0);
           dCdD = Matrix<double>(nD , nSplines, true, 0.0);
        }
        else{
           for (i = 0; i < nSplines; i++){
              ddCdAA[i] = Matrix<double>(nSplines - 1, nSplines - 1, true, 0.0);
              ddCdDD[i] = Matrix<double>(nD, nD, true, 0.0);
           }
           dA2dD = Matrix<double>(nD, 2, true, 0.0);
           ddAdDD1 = Matrix<double>(nD, nD, true, 0.0);
           ddAdDD2 = Matrix<double>(nD, nD, true, 0.0);
           dCdA2 = Matrix<double>(2, nSplines, true, 0.0);
           dCdAg3 = Matrix<double>(nD , nSplines, true, 0.0);
           dCdD = Matrix<double>(nD , nSplines, true, 0.0);
        }
     }

// 0th ITERATION
// ==============
  // if maxiter == 0, compute all derivatives using 'd' parametrization
  //    maxiter > 0, compute only value of the penalized log-likelihood in initials,
  //        to check them, use 'a' parametrization
   *penalloglikP = penalLogLik(ThetaMat, 3, 2*(*maxiterP == 0), (*maxiterP == 0));

   if (*penalloglikP < -FLT_MAX + 11){
      if (*penalloglikP < -FLT_MAX + 1){
         if (info) REprintf("      STOP: Initials off the probability scale.\n");
         *failP = 100;
      }
      else{
         if (info) REprintf("      STOP: log(negative) in initials.\n");
         *failP = 102;
      }
      delete [] ddCdAA;
      delete [] ddCdDD;
      return;
   }
   else{
      if (info){
         Rprintf("\nIter.: %d, ", *firstiterP);
         Rprintf("Lp = %g = %g + (%g)\n",
               *penalloglikP, logLikelihood, penalty);
         if (*maxiterP > 0 && estA)
                Rprintf("                                                           Constraints = %g,  %g\n", -minCon[0], -minCon[1]);
      }
   }


// REMAINING ITERATIONS
// =====================
   Matrix<double> useHMat;

   int iteration = *maxiterP;
   double newpenalloglik = 0.0;

   if (*maxiterP > 0){
      *failP = 0;

 // FIRST, OPTIMIZATION WITH CONSTRAINTS (conditional regres -> c's etc.)
 // ======================================================================

 // More stuff for quadratic programming
      Matrix<int> AindA;
      if (estA){
         AindA = Matrix<int>(2, nSplines, true, nSplines - 1);
         for (j = 1; j < nSplines; j++){
             AindA(0, j) = j;
             AindA(1, j) = j;
         }
      }

      int nCon = 2;
      double crvalQP;
      Matrix<int> iterQP = Matrix<int>(2, 1, false);
      j = (nTheta + 2 < nCon) ? nTheta + 2: nCon;
      k = int(0.5*j*(j+5));
      i = 2*(nTheta + 2) + k + 2*nCon + 1;      // ***** length of the working array
      Matrix<double> workQP = Matrix<double>(i, 1, false);
      double dummydouble = 0.0;
      int zeroint = 0;
      int dummyint = 0;
      Matrix<int> iactQPA = Matrix<int>(nCon, 1, true, 0);
      int nactQPA = nCon;
      int ierrQP;

  // Variables used in the loop over the iterations
      bool stopRegres, stopA;
      bool break_for_loop = false;
      double objectfun, newobjectfun;
      double relatdiff = 0.0;
      double normU;
      double maxXi;
      int halfstep = 0;
      Matrix<double> keepHMat;
      Matrix<double> useUMat;
      Matrix<double> NRstep = Matrix<double>(nTheta, 1, false);
      Matrix<double> NRstepRegres = Matrix<double>(nRegres, 1, false);
      Matrix<double> NRstepA = Matrix<double>(nD, 1, false);
      Matrix<double> NRzeroRegres = Matrix<double>(nRegres, 1, true, 0.0);
      Matrix<double> NRzeroA = Matrix<double>(nD, 1, true, 0.0);

      for (iteration = *firstiterP + 1; iteration <= *maxiterP; iteration++){
         stopRegres = false;
         stopA = estA ? false : true;

  // 1. UPDATE a's
  // ===============
         if (estA){
         // 1a) derivatives
            maxXi = (fabs(XiMat[0]) < fabs(XiMat[1])) ? fabs(XiMat[1]) : fabs(XiMat[0]);
            *penalloglikP = penalLogLik(ThetaMat, 2, 2, false);
            objectfun = *penalloglikP - maxXi * (fabs(minCon[0]) + fabs(minCon[1]));
            useHMat = HMatD;     // I need original Hessian to make it positive definite if it is not
            useUMat = UMatD;     // I need original score to update lagrangians

         // 1b) quadratic program (if necessary modify the second derivative matrix to have it pos. def.)
            kk = -1;
            do{
               ierrQP = 0;
               keepHMat = useHMat;   // I need really used (adjusted) Hessian to update lagrangians
               qpgen1CPP(useHMat.getArray(), useUMat.getArray(), &nD, &nD,
                   NRstepA.getArray(), &crvalQP,
                   dConInRows.getArray(), AindA.getArray(), minCon.getArray(),
                   &nD, &nCon, &nCon,
                   iactQPA.getArray(), &nactQPA,
                   iterQP.getArray(), workQP.getArray(),
                   &ierrQP, tolCholP);

                   switch (ierrQP){
                      case 0:         // OK
                      break;
                      case 1:         // no solution to the quadratic program
                        if (info) REprintf("      STOP: No solution to the quadratic program.\n");
                        *failP = 51;
                        break_for_loop = true;
                      break;
                      case 2:         // problems to decompose the second derivative matrix
                         kk++;
//                         if (info && kk == 0) REprintf("      WARNING: Minus Hessian was not positive definite.\n");
                         useHMat = HMatD;
                         k = createPosDef(useHMat.getArray(), nD, *tolEigenP * pow(10.0, kk));
                         if (k != 0){
                            if (info) REprintf("      STOP: Problems with eigen values decomposition.\n");
                            *failP = 53;
                            break_for_loop = true;
                            ierrQP = 1;   // to leave the 'do' sequence
                         }
                      break;
                   }
            } while (ierrQP == 2);
            if (break_for_loop) break;

         // 1c) find the optimal Lagrange multipliers from the quadratic program (update them)
            k = findLagrangeQP(XiMat, UMatD, keepHMat, NRstepA, dConInRows);
            if (k > 0){     // this should never happen...
               if (info) REprintf("      STOP: No optimal lagrange multipliers for the quadratic program.\n");
               *failP = 52;
               break;
            }

         // 1d) find the length of the step that increases the objective function
         //    (use all the time old lagrangians!!!)
            NRstep = rbind(NRzeroRegres, NRstepA);
            newThetaMat = ThetaMat + NRstep;
            for (halfstep = 0; halfstep < *maxhalfP; halfstep++){
               newpenalloglik = penalLogLik(newThetaMat, 2, 0*(halfstep == 0), false);
               newobjectfun = newpenalloglik - maxXi * (fabs(minCon[0]) + fabs(minCon[1]));
               relatdiff = fabs(1 - (objectfun/newobjectfun));
               if (relatdiff <= *epsP || newobjectfun >= objectfun) break;
               newThetaMat = 0.5*(ThetaMat + newThetaMat);
            }
            if (halfstep == *maxhalfP){
               if (info) REprintf("      STOP: Conditionally not converging.\n");
               *failP = 54;
               break;
            }

          // 1e) check the convergence and update the estimates
            normU = EucNorm(minCon);
            if (halfstep == 0 && relatdiff <= *epsP && normU <= *epsP) stopA = true;

            if (info){
                Rprintf("Iter.: %d(a coef),  ", iteration);
                Rprintf("Lp = %g = %g + (%g),  Rel.diff. = %g,  %d halfsteps\n",
                     newpenalloglik, logLikelihood, penalty, relatdiff, halfstep);
            if (*maxiterP > 0 && estA)
                Rprintf("                                                           Constraints = %g,  %g\n", -minCon[0], -minCon[1]);
            }

            ThetaMat = newThetaMat;

         }  // end of the a's update

  // 2. UPDATE REGRESSION PART
  // ==========================
      // 2a) derivatives
         *penalloglikP = penalLogLik(ThetaMat, 1, 2, false);
         useHMat = HMatRegres;   // I need original Hessian to make it positive definite if it is not

      // 2b) Cholesky decomposition and compute H^{-1}U
      //      (if necessary modify the second derivative matrix to have it pos. def.)
         kk = -1;
         do{
            ierrQP = 0;
            // if HMat is positive def. UMatRegres contains NR step
            //                          NRstepRegres contains the same NR step
            qpgen1CPP(useHMat.getArray(), UMatRegres.getArray(), &nRegres, &nRegres,
                NRstepRegres.getArray(), &crvalQP,
                &dummydouble, &dummyint, &dummydouble,
                &zeroint, &zeroint, &zeroint,
                &dummyint, &dummyint,
                iterQP.getArray(), workQP.getArray(),
                &ierrQP, tolCholP);

                switch (ierrQP){
                   case 0:         // OK
                   break;
                   case 2:         // problems to decompose the second derivative matrix
                      kk++;
//                      if (info && kk == 0) REprintf("      WARNING: Minus Hessian was not positive definite.\n");
                      useHMat = HMatRegres;
                      k = createPosDef(useHMat.getArray(), nRegres, *tolEigenP * pow(10.0, kk));
                      if (k != 0){
                         if (info) REprintf("      STOP: Problems with eigen values decomposition.\n");
                         *failP = 63;
                         break_for_loop = true;
                         ierrQP = 1;   // to leave the 'do' sequence
                      }
                   break;
               }
         } while (ierrQP == 2);
         if (break_for_loop) break;

      // 2d) find the length of the step that increases the objective function
      //    (use all the time old lagrangians!!!)
         NRstep = rbind(NRstepRegres, NRzeroA);
         newThetaMat = ThetaMat + NRstep;
         for (halfstep = 0; halfstep < *maxhalfP; halfstep++){
            newpenalloglik = penalLogLik(newThetaMat, 1, 0*(halfstep == 0), false);
            relatdiff = fabs(1 - (*penalloglikP/newpenalloglik));
            if (relatdiff <= *epsP || newpenalloglik >= *penalloglikP) break;
            newThetaMat = 0.5*(ThetaMat + newThetaMat);
         }
         if (halfstep == *maxhalfP){
            if (info) REprintf("      STOP: Conditionally not converging.\n");
            *failP = 64;
            break;
         }

       // 2e) check the convergence and update the estimates
         if (halfstep == 0 && relatdiff <= *epsP) stopRegres = true;

         if (info){
             Rprintf("Iter.: %d(regres),  ", iteration);
             Rprintf("Lp = %g = %g + (%g),  Rel.diff. = %g,  %d halfsteps\n",
                  newpenalloglik, logLikelihood, penalty, relatdiff, halfstep);
         }

         ThetaMat = newThetaMat;

      // end of the regression update

  // 3. CHECK THE CONDITIONAL CONVERGENCE
  // =====================================
         if (stopA && stopRegres){
            if (info) REprintf("      Conditional convergence reached.\n");
            break;
         }

      }    // end of the loop over the iterations (for the first part of optimization)

 // SECOND, FULL OPTIMIZATION WITHOUT CONSTRAINTS (only if estA)
 // ============================================================
    // Stop if problems in convergence of a's -> no output
      if (*failP > 50 && *failP < 60){
          *failP = 103;
          delete [] ddCdAA;
          delete [] ddCdDD;
          return;
      }

      if (estA){
         *penalloglikP = penalLogLik(ThetaMat, 3, 0, false);   // this will set Acoef

         nD = nSplines - 3;
         nTheta = nRegres + nD;
         nThetaSq = nTheta * nTheta;


      // Find which knots should be the reference ones
         Matrix<int> Aord = Matrix<int>(nSplines, 1, false);
         Matrix<int> Ainvord = Matrix<int>(nSplines, 1, false);
         Matrix<double> Asorted = sortOrder(Acoef, Aord, Ainvord);

         lastThreeA[1] = (fabs(knotsMat[Ainvord[nSplines - 1]]) > 2e-1) ?
                         Ainvord[nSplines - 1] : Ainvord[nSplines - 2];
         lastThreeA[2] = -1;
         for (i = nSplines - 1; i >= 0; i--){
             if (Ainvord[i] == lastThreeA[1]) continue;
             if (fabs(1 - sigmaZero * sigmaZero + knotsMat[lastThreeA[1]] * knotsMat[Ainvord[i]]) > 1e-4
                    && fabs(knotsMat[Ainvord[i]] - knotsMat[lastThreeA[1]]) > 1e-4){
                 lastThreeA[2] = Ainvord[i];
                 break;
             }
         }
         if (lastThreeA[2] < 0){
             lastThreeA[2] = (lastThreeA[1] == Ainvord[nSplines - 1]) ?
                             Ainvord[nSplines - 2] : Ainvord[nSplines - 3];
         }

         for (i = nSplines - 1; i >= 0; i--){
             if (Ainvord[i] == lastThreeA[1] || Ainvord[i] == lastThreeA[2]) continue;
             lastThreeA[0] = Ainvord[i];
             break;
         }


      // Recalculate indeces of the reference knots
         lastTwoAshort[0] = (lastThreeA[0] < lastThreeA[1] ? lastThreeA[1] - 1 : lastThreeA[1]);
         lastTwoAshort[1] = (lastThreeA[0] < lastThreeA[2] ? lastThreeA[2] - 1 : lastThreeA[2]);

         restA = Matrix<int>(nD, 1, false);
         for (i = 0, k = 0; k < nD; i++){
             if (i != lastThreeA[0] && i != lastThreeA[1] && i != lastThreeA[2]){
                 restA[k] = i;
                 k++;
             }
         }

         restAshort = Matrix<int>(nD, 1, false);
         for (i = 0, k = 0; k < nD; i++){
             if (i != lastTwoAshort[0] && i != lastTwoAshort[1]){
                 restAshort[k] = i;
                 k++;
             }
         }


      // Recalculate fixed derivative matrices
         OmegaExpA2 = Matrix<double>(nSplines - 2, 2, false);
         if (deriv_expAD(knotsMat, sigmaZero, lastThreeA, OmegaExpA2, false) != 0){
            if (info) REprintf("      STOP: I cannot find the reference knots.\n");
            *failP = 75;
            // change number of parameters back to +2 higher value
            nD = nSplines - 1;
            nTheta = nRegres + nD;
            nThetaSq = nTheta * nTheta;
         }
         else{
            OmegaZeroExpA2 = t(OmegaExpA2(0, 0, 0, 1));
            OmegaExpA2 = OmegaExpA2(1, 0, nD, 1);
            tOmegaExpA2 = t(OmegaExpA2);

        // block of minLamDDMat corresponding to the two special a's
            minLamDD2 = Matrix<double>(2, 2, false);
            minLamDD2(0, 0) = minLamDDMat(lastThreeA[1], lastThreeA[1]);
            minLamDD2(0, 1) = minLamDDMat(lastThreeA[1], lastThreeA[2]);
            minLamDD2(1, 0) = minLamDDMat(lastThreeA[2], lastThreeA[1]);
            minLamDD2(1, 1) = minLamDDMat(lastThreeA[2], lastThreeA[2]);

        // block of minLamDDMat corresponding to the cross-product of two a's and d's
            minLamDD2g3 = Matrix<double>(2, nD, false);
            for (i = 0; i < 2; i++)
               for (j = 0; j < nD; j++)
                   minLamDD2g3(i, j) = minLamDDMat(lastThreeA[i + 1], restA[j]);

        // block of minLamDDMat corresponding to d's
            minLamDDg3 = Matrix<double>(nD, nD, false);
            for (i = 0; i < nD; i++)
               for (j = 0; j < nD; j++)
                   minLamDDg3(i, j) = minLamDDMat(restA[i], restA[j]);


      // Re-initialize global variables which will be used in the following (often a dimension changed)
            GMatLEFT = Matrix<double>(nTheta, nRegres, true, 0.0);
            GMatUP = Matrix<double>(nRegres, nD, true, 0.0);

            UMat = Matrix<double>(nTheta, 1, true, 0.0);
            HMat = Matrix<double>(nTheta, nTheta, true, 0.0);
            GMat = Matrix<double>(nTheta, nTheta, true, 0.0);
            IMat = Matrix<double>(nTheta, nTheta, true, 0.0);


      // Initialize derivative matrices (also change of the dimension)
            for (i = 0; i < nSplines; i++){
               ddCdAA[i] = Matrix<double>(nSplines - 1, nSplines - 1, true, 0.0);
               ddCdDD[i] = Matrix<double>(nD, nD, true, 0.0);
            }
            dA2dD = Matrix<double>(nD, 2, true, 0.0);
            ddAdDD1 = Matrix<double>(nD, nD, true, 0.0);
            ddAdDD2 = Matrix<double>(nD, nD, true, 0.0);
            dCdA2 = Matrix<double>(2, nSplines, true, 0.0);
            dCdAg3 = Matrix<double>(nD , nSplines, true, 0.0);
            dCdD = Matrix<double>(nD , nSplines, true, 0.0);


      // Remove reference a's from the vector of parameters
      //   and subtract from all a's the value which is now assigned to a[lastThreeA[0]]
            dummydouble = Acoef[lastThreeA[0]];
            Acoef -= dummydouble;

            newThetaMat = Matrix<double>(nTheta, 1, true, 0.0);
            for (i = 0; i < nBeta; i++) newThetaMat[i] = ThetaMat[i];
            if (estScale) newThetaMat[i] = ThetaMat[i];
            if (estA){
               j = nBeta + nScale;
               for (k = 0; k < nSplines; k++){
                   if (k == lastThreeA[0] || k == lastThreeA[1] || k == lastThreeA[2]) continue;
                   newThetaMat[j] = Acoef[k];
                   j++;
               }
            }
            ThetaMat = newThetaMat;

         // Continue with iterations (full optimization)
            break_for_loop = false;
            NRstep = Matrix<double>(nTheta, 1, false);
            double normU;

            bool firstFull = true;
            *failP = 0;

            for (++iteration; iteration <= *maxiterP; iteration++){
               stopA = false;

  // 4. UPDATE THE WHOLE PARAMETERS VECTOR
  // ======================================
      // 4a) derivatives
               *penalloglikP = penalLogLik(ThetaMat, 3, 2, true);
               normU = EucNorm(UMat);

               if (firstFull){
                  if (info){
                      Rprintf("Iter.: %d,  ", iteration - 1);
                      Rprintf("Lp = %g = %g + (%g),  Norm(U) = %g\n",
                           newpenalloglik, logLikelihood, penalty, normU);
                  }
                  firstFull = false;
               }
               useHMat = HMat;       // I need original Hessian to make it positive definite if it is not

      // 4b) Cholesky decomposition and compute H^{-1}U
      //      (if necessary modify the second derivative matrix to have it pos. def.)
               kk = -1;
               do{
                  ierrQP = 0;
                  // if HMat is positive def. UMatRegres contains NR step
                  //                          NRstepRegres contains the same NR step
                  qpgen1CPP(useHMat.getArray(), UMat.getArray(), &nTheta, &nTheta,
                      NRstep.getArray(), &crvalQP,
                      &dummydouble, &dummyint, &dummydouble,
                      &zeroint, &zeroint, &zeroint,
                      &dummyint, &dummyint,
                      iterQP.getArray(), workQP.getArray(),
                      &ierrQP, tolCholP);

                      switch (ierrQP){
                         case 0:         // OK
                         break;
                         case 2:         // problems to decompose the second derivative matrix
                            kk++;
//                            if (info && kk == 0) REprintf("      WARNING: Minus Hessian was not positive definite.\n");
                            useHMat = HMat;
                            k = createPosDef(useHMat.getArray(), nTheta, *tolEigenP * pow(10.0, kk));
                            if (k != 0){
                               if (info) REprintf("      STOP: Problems with eigen values decomposition.\n");
                               *failP = 73;
                               break_for_loop = true;
                               ierrQP = 1;   // to leave the 'do' sequence
                            }
                         break;
                     }
               } while (ierrQP == 2);
               if (break_for_loop) break;

      // 4d) find the length of the step that increases the objective function
      //    (use all the time old lagrangians!!!)
               newThetaMat = ThetaMat + NRstep;
               for (halfstep = 0; halfstep < *maxhalfP; halfstep++){
                  newpenalloglik = penalLogLik(newThetaMat, 3, 1*(halfstep == 0), true);
                    relatdiff = fabs(1 - (*penalloglikP/newpenalloglik));
                  if (relatdiff <= *epsP || newpenalloglik >= *penalloglikP) break;
//                  if (halfstep == 0) normU = EucNorm(UMat);
//                  if (newpenalloglik >= *penalloglikP) break;
                  newThetaMat = 0.5*(ThetaMat + newThetaMat);
               }
               if (halfstep == *maxhalfP){
                  if (info) REprintf("      STOP: Not converging.\n");
                  *failP = 74;
                  break;
               }

       // 4e) check the convergence and update the estimates
//               if (halfstep == 0 && normU <= 10*nTheta*(*epsP)) stopA = true;
               if (halfstep == 0 && relatdiff <= *epsP) stopA = true;

               if (info){
                   Rprintf("Iter.: %d,  ", iteration);
//                   Rprintf("Lp = %g = %g + (%g),  Norm(U) = %g,  %d halfsteps\n",
//                        newpenalloglik, logLikelihood, penalty, normU, halfstep);
                   Rprintf("Lp = %g = %g + (%g),  Rel.diff. = %g,  %d halfsteps\n",
                        newpenalloglik, logLikelihood, penalty, relatdiff, halfstep);
               }

               ThetaMat = newThetaMat;

  // 5. CHECK THE FULL CONVERGENCE
  // =====================================
               if (stopA){
                  if (info) REprintf("      Full convergence reached.\n");
                  break;
               }

            }   // end of all iterations

         }   // end of else (it was possible to find reference knots)

      }   // end of the second part of optimization (if (estA))

   }   // end of if (*maxiterP > 0)


   if (iteration > *maxiterP) *failP = 16;
   else                       *maxiterP = iteration;

   if (*failP == 75){      // I could not find the reference knots -> no 'd' parametrization
      *penalloglikP = penalLogLik(ThetaMat, 3, 0, false);
      *nonPosDefHP = 1;
   }
   else{
      *penalloglikP = penalLogLik(ThetaMat, 3, 2, true);

      // test whether H is positive definite
      useHMat = HMat;
      dpofaCPP(useHMat.getArray(), nTheta, nTheta, nonPosDefHP, *tolCholP);

      for (i = 0; i < nTheta; i++) UP[i] = UMat[i];

      for (i = 0; i < nThetaSq; i++){
         HP[i] = HMat[i];
         GP[i] = GMat[i];
         IP[i] = IMat[i];
      }

      for (i = 0; i < 3; i++) lastThreeAP[i] = lastThreeA[i] + 1;   // give S indeces (i.e. starting from 1)

      dCdD = t(dCdD);
      for (i = 0; i < nSplines * nD; i++) dCdDP[i] = dCdD[i];
   }    // end of if (*failP == 75)

   *loglikP = logLikelihood;
   *penaltyP = penalty;

   for (i = 0; i < nBeta; i++) BetaP[i] = ThetaMat[i];
   if (estScale) *GammaP = ThetaMat[nBeta];

   for (i = 0; i < nSplines; i++){
      AcoefP[i] = Acoef[i];
      CcoefP[i] = Ccoef[i];
   }

// Compute matrices used to get df
//  - back only 'a' parametrization
   if (estA){

   // Re-initialize some global variables
      nD = nSplines - 1;
      nTheta = nRegres + nD;
      nThetaSq = nTheta * nTheta;

      restA = Matrix<int>(nD, 1, false);
      for (i = 0, k = 0; k < nD; i++){
         if (i != lastThreeA[0]){
             restA[k] = i;
             k++;
         }
      }

      LamDDg1 = Matrix<double>(nSplines - 1, nSplines - 1, false);
      for (i = 0; i < nSplines - 1; i++)
          for (j = 0; j < nSplines - 1; j++)
              LamDDg1(i, j) = -minLamDDMat(restA[i], restA[j]);

      for (i = 0; i < nSplines; i++) ddCdDD[i] = Matrix<double>(nD, nD, true, 0.0);
      dCdD = Matrix<double>(nD , nSplines, true, 0.0);

      minCon = Matrix<double>(2, 1, true, 0.0);
      dConInRows = Matrix<double>(2, nSplines - 1, true, 0.0);
      dCon0 = Matrix<double>(nSplines - 1, 1, true, 0.0);
      dCon1 = Matrix<double>(nSplines - 1, 1, true, 0.0);
      ddCon0 = Matrix<double>(nSplines - 1, nSplines - 1, true, 0.0);
      ddCon1 = Matrix<double>(nSplines - 1, nSplines - 1, true, 0.0);
      XiMat = Matrix<double>(2, 1, true, 0.0);

      UMatD = Matrix<double>(nD, 1, true, 0.0);
      HMatD = Matrix<double>(nD, nD, true, 0.0);

   // Add the two a's back to the vector of parameters
      newThetaMat = Matrix<double>(nTheta, 1, true, 0.0);
      for (i = 0; i < nBeta; i++) newThetaMat[i] = ThetaMat[i];
      if (estScale) newThetaMat[i] = ThetaMat[i];
      if (estA){
         j = nBeta + nScale;
         for (k = 0; k < nSplines; k++){
             if (k == lastThreeA[0]) continue;
             newThetaMat[j] = Acoef[k];
             j++;
         }
      }

   // Compute the derivatives and store them to be returned
   //   (store also derivatives of constraints)
      // !!! Zero multipliers are entering the following function
      // --> no impact of constraints on derivative matrices
      newpenalloglik = penalLogLik(newThetaMat, 2, 2, false);
      if (info) Rprintf("      Final constraints - mean: %g,  variance: %g\n", -minCon[0], -minCon[1]);

   // Find optimal lagrange multipliers
      k = findLagrange(XiMat, UMatD, dConInRows);
      if (k > 0){     // this should never happen...
         if (info) REprintf("      STOP: No optimal lagrange multipliers for the constraint optimization.\n");
         *failP = 150;
      }
      else{
         if (info) Rprintf("      Lagrange multipliers: %g,  %g\n\n", XiMat[0], XiMat[1]);
      }

      if (debug > 0){
         Matrix<double> LagrScore = UMatD - XiMat[0] * dCon0 - XiMat[1] * dCon1;
         cout << "Unadjusted 'a' score = " << t(UMatD).toString() << "\n";
         cout << "Lagrange 'a' score = " << t(LagrScore).toString() << "\n";
      }

      for (i = 0; i < nD * nD; i++){
         HaP[i] = HMatD[i] + XiMat[0] * ddCon0[i] + XiMat[1] * ddCon1[i];
         GaP[i] = LamDDg1[i];
         IaP[i] = HaP[i] - GaP[i];
      }
      for (i = 0; i < nSplines - 1; i++) dConP[i] = dCon0[i];
      for (i = 0; i < nSplines - 1; i++) dConP[i + nSplines - 1] = dCon1[i];
   }

   delete [] ddCdAA;
   delete [] ddCdDD;

   return;

}     // end of smoothSurvReg83 function

}     // end of extern "C"



