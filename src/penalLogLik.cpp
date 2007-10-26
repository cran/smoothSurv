// SUBROUTINES FOR smoothSurvReg84 function
// ========================================

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

// Function to compute a current value of
//    the penalized log-likelihood function
//    and to compute derivatives to be able
//    to update the estimate
// =================================================================

// INPUT:
//     Theta ......... current value of all parameters (beta, log(scale) pars., d_1, ..., d_{g-3})
//
//     what .......... what types of derivatives are to be computed
//               1 = only w.r.t. beta and log(scale) pars.
//               2 = only w.r.t. (d_1, ..., d_{g-3})'
//                              or (a_1, ..., a_{g-1})'
//               3 = joint
//                     if what == 3 I do allow estA == false
//
//     derivOrder ... which quantities do I want to compute
//               0 = only value
//               1 = also first derivatives
//               2 = also second derivatives
//     OnlyValue  - true if no derivatives are desired to be computed
//
//     useD ......... true if 'd' parametrization is used (i.e. if two a's
//                     are expressed as the function of the (g - 3) a's)
//                        false if 'a' parametrization is used (i.e. there are
//                        (g - 1) a's to be estimated)
//        !!! when switching between useD and !useD one has to adjust dimensions of
//            all appropriate matrices !!!

// RETURN:
//     value of the penalized log-likelihood function evaluated in INPUT
//     * it returns -FLT_MAX if off the probability scale
//     * it returns -FLT_MAX + 10 if log(negative) occurs in computation of a's
//     * it returns -FLT_MAX + 20 if sum(exp(a)) is infinity

// SET GLOBAL VARIABLES INSIDE THE FUNCTION:
//     penalty
//
//     logLikelihood
//
//     if what == 1 then further UMatRegres
//                               HMatRegres
//     if what == 2 then further UMatD
//                               HMatD
//     if what == 3 then further UMat
//                               HMat
//                               IMat

#include "penalLogLik.h"

using namespace std;
using namespace SCYTHE;

// GLOBAL VARIABLES DEFINED in smoothSurvReg82.cpp
// ==============================================
   extern int n,                        // number of observations
              nTheta,                   // number of all parameters to be really estimated
              nThetaSq,                 // nTheta^2
              nBeta,                    // number of columns in X matrix (number of beta parameters to be estimated)
              nGamma,                   // number of columns in Z matrix (number of s parameters to be estimated)
              nScale,                   // number of scale parameters to be estimated (either 1 or 0)
              nRegres,                  // nBeta + nScale
              nD,                       // number of a coefficients to be estimated (either g-3 or g-1 or 0)
              nSplines,                 // number of G-splines
              debug;

   extern double lambda,                // tuning parameter for the penalty term
                 sigmaZero,             // standard deviation of one basis "spline"
                 sigmaZeroInv,          // inversion of the standard deviation of one basis "spline" (added on 15/08/2007)
                 invsigmaZero,          // inversion of sigmaZero
                 logsigmaZero;          // log(sigmaZero)

  extern Matrix<double> oneMat;         // 1x1 matrix containing 1 (added on 15/08/2007)


   extern bool estScale,                // true if scale is to be estimated
               estA,                    // true if a's are to be estimated
               info;                    // do I want to print some information during the iteration process?

   extern double correctLik,             // correction to likelihood due to log transformation of the response
                 logLikelihood,          // corrected value of unpenalized log-likelihood
                 penalty;                // penalty term

  // Response and covariates matrices
   extern Matrix<int> statusMat;           // status vector (1=exact, 0=right censored, 3=interval censored)
   extern Matrix<double> resp1Mat,         // response (exact, right censored or lower limit)  n x 1
                         resp2Mat,         // upper limit of response for interval censored data  n x 1
                         XMat,             // design matrix (covariates)  n x nBeta
                         ZMat,             // design matrix (covariates) for log(scale) n x nGamma
                         offset;           // offset vector n x 1

  // Matrices and vectors used in Newton-Raphson steps
   extern Matrix<double> UMat,        // score vector
                         UMatUP,      // zero vector ((nBeta+nScale) x 1)
                         HMat,        // minus Hessian matrix (minus second derivative) for the penalized log-likelihood
                         IMat,        // minus Hessian matrix for UNpenalized log-likelihood
                         GMat,        // minus Hessian matrix for penalty term
                         GMatLEFT,      // building block for GMat
                         GMatUP,        // building block for GMat
                         UMatRegres,
                         HMatRegres,
                         UMatD,
                         HMatD;
// A and C coefficients
   extern Matrix<double> Acoef,         // (a_1, ..., a_g)'
                         Ccoef;         // (c_1, ..., c_g)'

// Lagrangians (for optimization with equality constraints)
   extern Matrix<double> XiMat;         // vector 2 x 1


  // Equality constraints and their derivatives
   extern Matrix<double> minCon;
   extern Matrix <double> dConInRows;
   extern Matrix<double> dCon0, dCon1;
   extern Matrix<double> ddCon0, ddCon1;


  // Splines related vectors and matrices
   extern Matrix<double> knotsMat,       // knots
                         knotsMatSq,     // squared knots
                         OmegaExpA2,     // matrix with d(exp(A1, A2))/d(exp(D)) (g-3) x 2
                         tOmegaExpA2,    // its transposition                    2 x (g-3)
                         OmegaZeroExpA2, // intercept to compute the two a's     2 x 1
                                // ==> a[the two] = tOmegaExpA2 * d + OmegaZeroExpA2
                         minLamDDMat,    // matrix -lambda*D'*D (matrix g x g)
                         LamDDg1,        // matrix (g - 1) x (g - 1)
                         minLamDD2,      // matrix 2 x 2
                         minLamDD2g3,    // matrix 2 x (g - 3)
                         minLamDDg3;     // matrix (g - 3) x (g - 3)

  // Various derivatives
   extern Matrix<double> * ddCdAA;
   extern Matrix<double> * ddCdDD;
   extern Matrix<double> dA2dD,
                         ddAdDD1,
                         ddAdDD2,
                         dCdA2,
                         dCdAg3,
                         dCdD;

  // Matrices with possibly fixed some parameters
   extern Matrix<double> GammaGlobal;

  // Matrix indicating which three c's are expressed as a function of the remaining ones
   extern Matrix<int> lastThreeA,     // indeces of the zero a and of the two a's
                      restA,          // indeces of (g - 3) a's which are real parameters
                      lastTwoAshort,  // indeces of the two a's in the sequence of a's with removed zero
                      restAshort;     // indeces of (g - 3) a's in the sequence of a's with removed zero


// REMARK: nD is always a number of a's which are really estimated,
//         i.e. = nSplines - 3 if useD
//              = nSplines - 1 if !useD

// =========
// FUNCTION
// =========
double
penalLogLik(const Matrix<double> & Theta,
            const int what,
            const int derivOrder,
            const bool useD)
{
   int i, j, k;           // variables for 'for' loop
   int person;            // variable for 'for' loop


// EXTRACT BETA, LOG(SCALE) PARS., D'S FROM THE VECTOR OF PARAMETERS
// ==================================================================
   Matrix<double> Beta, Gamma, Dcoef;
   if (n > 0){
     Beta = Theta(0, 0, nBeta - 1, 0);
     if (estScale) Gamma = Theta(nBeta, 0, nBeta + nGamma - 1, 0);
     else          Gamma = GammaGlobal;
   }
   if (estA){
      Dcoef = Theta(nBeta+nScale, 0, nBeta+nScale+nD-1, 0);
   }


// SPLINE COEFFICIENTS (A and C COEFFICIENTS)
// ===========================================
   Matrix<double> exp_Acoef_two;
   Matrix<double> exp_Dcoef;
   Matrix<double> Acoef_two;

   if (estA){
       exp_Dcoef = exp(Dcoef);
       if (useD){
          exp_Acoef_two = tOmegaExpA2 * exp_Dcoef + OmegaZeroExpA2;
          if (exp_Acoef_two[0] <= 0.0 || exp_Acoef_two[1] <= 0.0)
              return (-FLT_MAX + 10);
          Acoef_two = log(exp_Acoef_two);
          Acoef[lastThreeA[1]] = Acoef_two[0];
          Acoef[lastThreeA[2]] = Acoef_two[1];
       }

       for (i = 0; i < nD; i++) Acoef[restA[i]] = Dcoef[i];
       Acoef[lastThreeA[0]] = 0.0;

       if (A_to_C(Acoef, Ccoef) != 0) return (-FLT_MAX + 20);

   }      // end of if (estA)

   Matrix<double> tCcoef = t(Ccoef);
   Matrix<double> tAcoef = t(Acoef);


// VARIOUS DERIVATIVES OF C'S, A'S and D'S
// ===========================================
   if (estA && (derivOrder > 0) && (what != 1)){

      if (useD){
         compute_dA2dD(dA2dD, exp_Acoef_two, exp_Dcoef, OmegaExpA2);

         if (derivOrder > 1){
            compute_ddA2dDD(ddAdDD1, ddAdDD2, dA2dD, nD);
            compute_dCdA(dCdA2, dCdAg3, Ccoef, tCcoef, lastThreeA, restA, nD);
         }
      }
      compute_dCdD(dCdD, dCdA2, dCdAg3, dA2dD, Ccoef, tCcoef, restA, nD, useD);

      if (derivOrder > 1){
         if (useD){
            compute_ddCdAA(ddCdAA, Ccoef, lastThreeA[0]);
            compute_ddCdDD(ddCdDD, dCdA2, ddCdAA, dA2dD, ddAdDD1, ddAdDD2, lastTwoAshort, restAshort);
         }
         else{  // store ddCdAA in ddCdDD
            compute_ddCdAA(ddCdDD, Ccoef, lastThreeA[0]);
         }
      }

   }


// SET GLOBALS TO BE CHANGED TO ZEROS (mostly) (also GMat!!!)
// ==========================================================
   // First derivatives
   if (derivOrder > 0){
      switch (what){
         case 1:
            UMatRegres = Matrix<double>(nBeta+nScale, 1, true, 0.0);
         break;
         case 2:
            UMatD = Matrix<double>(nD, 1, true, 0.0);
         break;
         case 3:
            UMat = Matrix<double>(nTheta, 1, true, 0.0);
         break;
      }
   }

   // Second derivatives
   if (derivOrder > 1){
      switch (what){
         case 1:
            HMatRegres = Matrix<double>(nBeta+nScale, nBeta+nScale, true, 0.0);
         break;
         case 2:
            HMatD = Matrix<double>(nD, nD, true, 0.0);
         break;
         case 3:
            IMat = Matrix<double>(nTheta, nTheta, true, 0.0);
            HMat = Matrix<double>(nTheta, nTheta, true, 0.0);
            GMat = Matrix<double>(nTheta, nTheta, true, 0.0);
         break;
      }
   }

   // Values
   logLikelihood = (n > 0) ? correctLik : 0;                        // not zero but the correction factor (if there are some data)
   penalty = 0;


// **********************************************************************************
// PENALTY TERM AND ITS DERIVATIVES
// ==================================
   Matrix<double> dQD = Matrix<double>(nD, 1, false);  // d(penalty)/dD;  (g - 3) x 1 or (g - 1) x 1
   Matrix<double> minLamDDA;                           // -lambda * DD * a(all);  g x 1
   Matrix<double> minLamDDA2;                          // -(lambda * DD * a(all))[the two]
   Matrix<double> minLamDDAg3;                         // -(lambda * DD * a(all))[g-3 comp.]
   Matrix<double> Gaa;                                 // minus second derivative of the penalty w.r.t. d's

   if (estA){
      minLamDDA = minLamDDMat * Acoef;
      penalty = 0.5 * (tAcoef * minLamDDA)[0];
   }

 // Derivatives of the penalty w.r.t. d's
   if ((derivOrder > 0) && estA && (what != 1)){
      if (useD){
         minLamDDA2 = Matrix<double>(2, 1, false);        // 2 x 1
         minLamDDAg3 = Matrix<double>(nD, 1, false);      // (g - 3) x 1
         minLamDDA2[0] = minLamDDA[lastThreeA[1]];
         minLamDDA2[1] = minLamDDA[lastThreeA[2]];
         for (i = 0; i < nD; i++) minLamDDAg3[i] = minLamDDA[restA[i]];
         dQD = dA2dD * minLamDDA2 + minLamDDAg3;
      }
      else{
         for (i = 0; i < nD; i++) dQD[i] = minLamDDA[restA[i]];
      }


   // d^2(penalty)/dD dD^T, nD x nD = (g - 3) x (g - 3) matrix
   // or (g - 1) x (g - 1) matrix
   // ---------------------------------------------------------
      if (derivOrder > 1){
         if (useD){
            Gaa = (dA2dD * minLamDD2 * t(dA2dD)) + minLamDDg3 +
                     (dA2dD * minLamDD2g3) + (t(minLamDD2g3) * t(dA2dD));
            Gaa += minLamDDA2[0] * ddAdDD1 + minLamDDA2[1] * ddAdDD2;
            Gaa *= (-1);
         }
         else{
            Gaa = LamDDg1;
         }
      }

   // Add appropriate values to UMat or UMatD and to GMat
      if (what == 3){
	 if (n > 0){ 
            UMat = rbind(UMatUP, dQD);
            if (derivOrder > 1){
               GMat = cbind(GMatLEFT, rbind(GMatUP, Gaa));
               HMat = GMat;
            }
	 }
         else{
           UMat = dQD;
           if (derivOrder > 1) HMat = Gaa;
         }
      }
      else{   // what == 2
         UMatD = dQD;
         if (derivOrder > 1) HMatD = Gaa;
      }
   }         // end of if (derivOrder > 0 && estA && (what != 1))


// **********************************************************************************
// UNPENALIZED LOG-LIKELIHOOD AND ITS DERIVATIVES
// ==============================================

   if (n > 0){

   // REGRESSION VARIABLES
   // =====================
      Matrix<double> scale, logscale;
      if (estScale) logscale = ZMat * Gamma;
      else          logscale = Matrix<double>(n, 1, true, Gamma[0]);
      scale = exp(logscale);
      Matrix<double> eta = XMat * Beta + offset;       // linear predictor  (n x 1)
      Matrix<double> w1Mat = (resp1Mat - eta)/scale;   // censored residuals (or lower limits of residuals in interval censored case) (n x 1)
      Matrix<double> ss0 = oneMat / (scale*sigmaZero);
      Matrix<double> s2s02 = ss0 & ss0;
      double s0 = invsigmaZero;
      double s02 = s0 * s0;
      Matrix<double> ss02 = (oneMat / scale) * s02;

    // Variables used in the loop over persons,
    //   variables used for interval censored observations have "subscript" 2
      double loglperson = 0.0;          // personal contribution to log(likelihood)
      Matrix<double> xvec;              // personal regression covariates
      Matrix<double> zvec;              // personal covariates for log(sigma)

      double w1, w2;                        // censored residuals w1 = lower limit = (y1 - eta)/sigma
      Matrix<double> wm1, wm2,              // personal (w1[i] - knot[j])/sigma0 (g x 1)
                     phiwm1, phiwm2,        // personal phi(wm1) (g x 1)
                     phiwmDot1, phiwmDot2,  // personal wm1 & phi(wm1) (g x 1)
                     phiwmDotDot1,          // personal (wm1*wm1 - 1) & phi(wm1) (g x 1)
                     Swm1, Swm2,            // personal S(wm1) or (F(wm1) for left censored obs.) (g x 1)
                     deltaSwm;              // personal S(wm1) - S(wm2) for interval censored observations
      double fw1, fw2,                      // personal sigma0 * f(w1)
             fwDot1, fwDot2,                // personal -sigma0^2 * f'(w1)
             fwDotDot1,                     // personal sigma0^3 * f''(w1)
             Sw1;                           // personal S(w1) or F(wm1) for left censored observation
      double deltaS;                        // personal S(w1) - S(w2) for interval censored obs.
      double dratio,                        // ratio involved in db, dg
             ddratio,                       // ratio involved in ddbb, ddgg, ddbg
             invDensity;                    // inversion of either density or survivor or PDF or deltaS
      Matrix<double> vec_bg, vec_bg2;

    // Components to store derivatives
    // d** and D** stuff changes over persons
    // U** and I** stuff sums over persons
      double db, dg, ddbb, ddgg, ddbg;
      double semidg = 0.0;
      Matrix<double> Db, Da,
                     DDba, DDga, DDaa1, DDaa2,
                     Ubeta, Ugamma, Ibb, Ibg, Igg,
                     Ua, Iaa, Iba, Iga;

    // Initialize by zeros these which will be summed up.
    //   Some of them may not be used (if what != 3 or derivOrder < 2)
    //   but initialization does not take almost any time and memory.
      if (derivOrder > 0){
         Ubeta = Matrix<double>(nBeta, 1, true, 0.0);
         Ibb = Matrix<double>(nBeta, nBeta, true, 0.0);
         if (estScale){
            Ugamma = Matrix<double>(nGamma, 1, true, 0.0);
            Igg = Matrix<double>(nGamma, nGamma, true, 0.0);
            Ibg = Matrix<double>(nBeta, nGamma, true, 0.0);
         }
         if (estA){
            Ua = Matrix<double>(nD, 1, true, 0.0);
            Iaa = Matrix<double>(nD, nD, true, 0.0);
            Iba = Matrix<double>(nBeta, nD, true, 0.0);
         }
         if (estScale && estA){
            Iga = Matrix<double>(nGamma, nD, true, 0.0);
         }
      }

    // Helping matrices
    // Matrix ww1(i, j) = w1[i]
      Matrix<double> ww1 = w1Mat;
      for(j = 1; j < nSplines; j++){        // create a matrix with g columns equal to w1, matrix (n x g)
         ww1 = cbind(ww1, w1Mat);
      }

    // Matrix mm(i, j) = knot[j]
      Matrix<double> mm = t(knotsMat);
      for(i = 1; i < n; i++){               // create a matrix with n rows equal to knots, matrix (n x g)
        mm = rbind(mm, t(knotsMat));
      }

    // Following quantities are used for all possible censoring patterns
      Matrix<double> wm1Mat = sigmaZeroInv*(ww1 - mm);        // matrix (n x g), wm1Mat(i, j) = (w1[i] - mu[j]) / sigma0

    // Loop over persons to compute derivatives and likelihood contributions
      for (person = 0; person < n; person++) {
         xvec = t(XMat(person, 0, person, nBeta - 1));                 // covariate values for given person, vector (nBeta x 1)
         if (estScale) zvec = t(ZMat(person, 0, person, nGamma - 1));
         w1 = w1Mat[person];                                           // residual
         wm1 = t(wm1Mat(person, 0, person, nSplines - 1));             // (w1 - knots) / sigmaZero   (nknots x 1)

         switch(statusMat[person]) {
  
            // == exact observation ==
            // **************************
            case 1:
                phiwm1 = fnorm(wm1);
                fw1 = (tCcoef * phiwm1)[0];
                if (fw1 <= 0){
                    logLikelihood = -FLT_MAX;
                    return -FLT_MAX;
                }
//                if (fw1 <= 0) fw1 = ZERO;
                loglperson = -logscale[0] + log(fw1) - logsigmaZero;

                if (derivOrder > 0){
                   if ((what == 1) || (what == 3)){
                      phiwmDot1 = wm1 & phiwm1;
                      fwDot1 = (tCcoef * phiwmDot1)[0];
                      dratio = fwDot1 / fw1;

                      db = ss0[person] * dratio;
                      Db = db * xvec;
                      if (estScale){
                         semidg = s0 * w1 * dratio;
                         dg = -1 + semidg;
                      }

                      if (derivOrder > 1){
                         phiwmDotDot1 = ((wm1 & wm1) - 1.0) & phiwm1;
                         fwDotDot1 = (tCcoef * phiwmDotDot1)[0];
                         ddratio = fwDotDot1 / fw1;

                         ddbb = s2s02[person] * ddratio - (db * db);
                         if (estScale){
                            ddgg = s02 * w1 * w1 * ddratio - semidg * (1 + semidg);
                            ddbg = ss02[person] * w1 * ddratio - db * (1 + semidg);
                         }
                      }
                   }
 
                   if (((what == 2) || (what == 3)) && estA){
                      invDensity = 1 / fw1;
                      Da = invDensity * dCdD * phiwm1;

                      if (derivOrder > 1){
                         DDaa1 = -Da * t(Da);
                         DDaa2 = phiwm1[0] * ddCdDD[0];
                         for (k = 1; k < nSplines; k++)
                            DDaa2 += phiwm1[k] * ddCdDD[k];
                         DDaa2 *= invDensity;
                      }
                   }

                   if ((what == 3) && estA && (derivOrder > 1)){
                      vec_bg = invDensity * t(dCdD * phiwmDot1);
                      DDba = -Db * t(Da) + ss0[person] * xvec * vec_bg;
                      if (estScale)
                         DDga = -semidg * t(Da) + s0 * w1 * vec_bg;
                   }
                }
            break;

            // === right censored observation ==
            // ***********************************
            case 0:
                Swm1 = Snorm(wm1);
                Sw1 = (tCcoef * Swm1)[0];
                if (Sw1 <= 0){
                    logLikelihood = -FLT_MAX;
                    return -FLT_MAX;
                }
//                if (Sw1 <= 0) Sw1 = ZERO;
                loglperson = log(Sw1);

                if (derivOrder > 0){
                   if ((what == 1) || (what == 3)){
                      phiwm1 = fnorm(wm1);
                      fw1 = (tCcoef * phiwm1)[0];
                      dratio = fw1 / Sw1;

                      db = ss0[person] * dratio;
                      Db = db * xvec;
                      if (estScale){
                         dg = s0 * w1 * dratio;
                      }

                      if (derivOrder > 1){
                         phiwmDot1 = wm1 & phiwm1;
                         fwDot1 = (tCcoef * phiwmDot1)[0];
                         ddratio = fwDot1 / Sw1;

                         ddbb = s2s02[person] * ddratio - (db * db);
                         if (estScale){
                            ddgg = s02 * w1 * w1 * ddratio - dg * (1 + dg);
                            ddbg = ss02[person] * w1 * ddratio - db * (1 + dg);
                         }
                      }
                   }

                   if (((what == 2) || (what == 3)) && estA){
                      invDensity = 1 / Sw1;
                      Da = invDensity * dCdD * Swm1;

                      if (derivOrder > 1){
                         DDaa1 = -Da * t(Da);
                         DDaa2 = Swm1[0] * ddCdDD[0];
                         for (k = 1; k < nSplines; k++)
                            DDaa2 += Swm1[k] * ddCdDD[k];
                         DDaa2 *= invDensity;
                      }
                   }

                   if ((what == 3) && estA && (derivOrder > 1)){
                      vec_bg = invDensity * t(dCdD * phiwm1);
                      DDba = -Db * t(Da) + ss0[person] * xvec * vec_bg;
                      if (estScale)
                         DDga = -dg * t(Da) + s0 * w1 * vec_bg;
                   }
                }
            break;

            // === left censored observation ==
            // ***********************************
            case 2:
                Swm1 = Fnorm(wm1);
                Sw1 = (tCcoef * Swm1)[0];

                if (Sw1 <= 0){
                   logLikelihood = -FLT_MAX;
                   return -FLT_MAX;
                }
//                if (Sw1 <= 0) Sw1 = ZERO;
                loglperson = log(Sw1);
  
                if (derivOrder > 0){
                   if ((what == 1) || (what == 3)){
                      phiwm1 = fnorm(wm1);
                      fw1 = (tCcoef * phiwm1)[0];
                      dratio = -(fw1 / Sw1);

                      db = ss0[person] * dratio;
                      Db = db * xvec;
                      if (estScale){
                         dg = s0 * w1 * dratio;
                      }

                      if (derivOrder > 1){
                         phiwmDot1 = wm1 & phiwm1;
                         fwDot1 = (tCcoef * phiwmDot1)[0];
                         ddratio = -(fwDot1 / Sw1);

                         ddbb = s2s02[person] * ddratio - (db * db);
                         if (estScale){
                            ddgg = s02 * w1 * w1 * ddratio - dg * (1 + dg);
                            ddbg = ss02[person] * w1 * ddratio - db * (1 + dg);
                         }
                      }
                   }
  
                   if (((what == 2) || (what == 3)) && estA){
                      invDensity = 1 / Sw1;
                      Da = invDensity * dCdD * Swm1;

                      if (derivOrder > 1){
                         DDaa1 = -Da * t(Da);
                         DDaa2 = Swm1[0] * ddCdDD[0];
                         for (k = 1; k < nSplines; k++)
                            DDaa2 += Swm1[k] * ddCdDD[k];
                         DDaa2 *= invDensity;
                      }
                   }

                   if ((what == 3) && estA && (derivOrder > 1)){
                      vec_bg = invDensity * t(dCdD * phiwm1);
                      DDba = -Db * t(Da) - ss0[person] * xvec * vec_bg;
                      if (estScale)
                         DDga = -dg * t(Da) - s0 * w1 * vec_bg;
                   }
                }
            break;

            // === interval censored observation ==
            // ***********************************
            case 3:
                Swm1 = Snorm(wm1);
                w2 = (resp2Mat[person] - eta[person])/scale[person]; // upper censored residual
                wm2 = sigmaZeroInv*(w2 - knotsMat);                  // vector (g x 1)
                Swm2 = Snorm(wm2);
                deltaSwm = Swm1 - Swm2;                              // vector (g x 1)
                deltaS = (tCcoef * deltaSwm)[0];
                if (deltaS <= 0){
                   logLikelihood = -FLT_MAX;
                   return -FLT_MAX;
                }
//                if (deltaS <= 0) deltaS = ZERO;
                loglperson = log(deltaS);

                if (derivOrder > 0){
                   if ((what == 1) || (what == 3)){
                      phiwm1 = fnorm(wm1);
                      fw1 = (tCcoef * phiwm1)[0];
                      phiwm2 = fnormZero(wm2);                            // vector (g x 1)
                      fw2 = (tCcoef * phiwm2)[0];

                      db = ss0[person] * ((fw1 -fw2)/deltaS);
                      Db = db * xvec;
                      if (estScale){
                         dg = s0 * ((w1*fw1 - w2*fw2)/deltaS);
                      }

                      if (derivOrder > 1){
                         phiwmDot1 = wm1 & phiwm1;
                         fwDot1 = (tCcoef * phiwmDot1)[0];
                         phiwmDot2 = wm2 & phiwm2;                           // vector (g x 1)
                         fwDot2 = (tCcoef * phiwmDot2)[0];

                         ddbb = s2s02[person] * ((fwDot1-fwDot2)/deltaS) - (db * db);
                         if (estScale){
                            ddgg = s02 * ((w1*w1*fwDot1 - w2*w2*fwDot2)/deltaS) - dg * (1 + dg);
                            ddbg = ss02[person] * ((w1*fwDot1 - w2*fwDot2)/deltaS) - db * (1 + dg);
                         }
                      }
                   }

                   if (((what == 2) || (what == 3)) && estA){
                      invDensity = 1 / deltaS;
                      Da = invDensity * dCdD * deltaSwm;
   
                      if (derivOrder > 1){
                         DDaa1 = -Da * t(Da);
                         DDaa2 = deltaSwm[0] * ddCdDD[0];
                         for (k = 1; k < nSplines; k++)
                            DDaa2 += deltaSwm[k] * ddCdDD[k];
                         DDaa2 *= invDensity;
                      }
                   }

                   if ((what == 3) && estA && (derivOrder > 1)){
                      vec_bg = invDensity * t(dCdD * phiwm1);
                      vec_bg2 = invDensity * t(dCdD * phiwm2);
                      DDba = -Db * t(Da) + ss0[person] * xvec * (vec_bg - vec_bg2);
                      if (estScale)
                         DDga = -dg * t(Da) + s0 * (w1*vec_bg - w2*vec_bg2);
                   }
                }
            break;

         }   // end of switch(statusMat[person])


         // add the value
         logLikelihood += loglperson;

         // add the first derivatives
         if (derivOrder > 0){
            if ((what == 1) || (what == 3)){
               Ubeta += Db;
               if (estScale) Ugamma += dg * zvec;
            }
            if (((what == 2) || (what == 3)) && estA) Ua += Da;
         }

         // add the second derivatives
         if (derivOrder > 1){
            if ((what == 1) || (what == 3)){
               Ibb -= ddbb * xvec * t(xvec);
               if (estScale){
                  Igg -= ddgg * zvec * t(zvec);
                  Ibg -= ddbg * xvec * t(zvec);
               }
            }
            if (((what == 2) || (what == 3)) && estA){
               Iaa -= (DDaa1 + DDaa2);
            }
            if ((what == 3) && estA){
               Iba -= DDba;
               if (estScale)
                  Iga -= zvec * DDga;
            }
         }

      }   // end of the loop over persons


     // Add appropriate values to UMat, create IMat
     //   and put everything together

      // first derivatives
      if (derivOrder > 0){
         switch (what){
            case 1:
               if (estScale) UMatRegres = rbind(Ubeta, Ugamma);
               else          UMatRegres = Ubeta;                  // both a and scale are fixed
            break;

            case 2:
               UMatD += Ua;
            break;

            case 3:
               if (estA && estScale) UMat += rbind(rbind(Ubeta, Ugamma), Ua);
               else
                  if (estA) UMat += rbind(Ubeta, Ua);              // scale is fixed
                  else
                     if (estScale) UMat += rbind(Ubeta, Ugamma);   // a is fixed
                     else          UMat += Ubeta;                 // both a and scale are fixed
            break;
         }  // end of switch
      }     // end of if (derivOrder > 0)

      // second derivatives
      if (derivOrder > 1){
         Matrix<double> IMatRow1, IMatRow2, IMatRow3;
         switch (what){
            case 1:
               if (estScale){
                  IMatRow1 = cbind(Ibb, Ibg);
                  IMatRow2 = cbind(t(Ibg), Igg);
                  HMatRegres = rbind(IMatRow1, IMatRow2);
               }
               else{           // both a and scale are fixed
                  HMatRegres = Ibb;
               }
            break;

            case 2:
               HMatD += Iaa;
            break;

            case 3:
               if (estA && estScale){
                  IMatRow1 = cbind(cbind(Ibb, Ibg), Iba);
                  IMatRow2 = cbind(cbind(t(Ibg), Igg), Iga);
                  IMatRow3 = cbind(cbind(t(Iba), t(Iga)), Iaa);
                  IMat = rbind(rbind(IMatRow1, IMatRow2), IMatRow3);
               }
               else
                  if (estA){     // scale is fixed
                     IMatRow1 = cbind(Ibb, Iba);
                     IMatRow3 = cbind(t(Iba), Iaa);
                     IMat = rbind(IMatRow1, IMatRow3);
                  }
                  else
                     if (estScale){   // a is fixed
                        IMatRow1 = cbind(Ibb, Ibg);
                        IMatRow2 = cbind(t(Ibg), Igg);
                        IMat = rbind(IMatRow1, IMatRow2);
                     }
                     else{           // both a and scale are fixed
                        IMat = Ibb;
                     }
               HMat += IMat;
            break;
         }  // end of switch
      }     // end of if (derivOrder > 1)
   }        // end of if (n > 0)


// **********************************************************************************
// EQUALITY CONSTRAINTS AND THEIR DERIVATIVES (if needed)
// =========================================================

   if ((what != 1) && (estA) && (!useD)){
      minCon[0] = -(tCcoef * knotsMat)[0];
      minCon[1] = -(tCcoef * knotsMatSq)[0] - sigmaZero * sigmaZero + 1.0;

      if (derivOrder > 0){
         dCon0 = dCdD * knotsMat;
         dCon1 = dCdD * knotsMatSq;
         dConInRows = rbind(t(dCon0), t(dCon1));

         if (derivOrder > 1){
            ddCon0 = knotsMat[0] * ddCdDD[0];
            ddCon1 = knotsMatSq[0] * ddCdDD[0];
            for (k = 1; k < nSplines; k++){
               ddCon0 += knotsMat[k] * ddCdDD[k];
               ddCon1 += knotsMatSq[k] * ddCdDD[k];
            }

         // add combination of second derivative of constr. to HMat
            Matrix<double> ddConComb = XiMat[0] * ddCon0 + XiMat[1] * ddCon1;
            if (what == 3){
              if (n > 0) HMat += cbind(GMatLEFT, rbind(GMatUP, ddConComb));
              else       HMat += ddConComb;
            }
            else         HMatD += ddConComb;              // what == 2
         }
      }
   }

   double penal_log_lik = logLikelihood + penalty;

   return penal_log_lik;

}    // end of penalLogLik function


// ========================================================================================
// Function used to find optimal values of the Lagrangians according to
// the solution to the quadratic program
//
// REMARK: At this moment implemented only for TWO EQUALITY CONSTRAINTS

// Quadratic function is assumed to be
// q(d) = U^T * d - 0.5 * d^T * H * d

// Lagrange function is assumed to have a form
// L(d, xi) = q(d) - sum(xi * constraint),
// Optimal multipliers are thus a solution to
// sum(xi * d(constraint)) = U - H * d

// !!! I do not do any checks of the input arguments
//     It is assumed to be used (properly) inside the smoothSurvReg function

// INPUT:  U ...... vector (p x 1) of the linear part of the quadratic function
//         H ...... symmetric matrix (p x p) from the quadratic part of the
//                  quadratic function
//         dvec ... optimal solution to the quadratic program (p x 1)
//         dCon ... matrix (q x p) with derivatives of the constraints in ROWS

// OUTPUT: Xi ..... vector (q x 1) with resulting Lagrangians

// RETURN: 0 if no problems
//         >0 if there are problems
//         99 if there is no solution for optimal Lagrangians
//            (this should theoretically never happen...)
int
findLagrangeQP(Matrix<double> & Xi,
               const Matrix<double> & U, const Matrix<double> & H,
               const Matrix<double> & dvec,
               const Matrix<double> & dCon)
{
//   const double ZERO = 1e-8;

   int p = U.size();
   int q = dCon.rows();
   if (H.rows() != p || H.cols() != p || dCon.cols() != p ||
          Xi.size() != q || dvec.size() != p || dvec.cols() != 1 || p < q){
       REprintf("Error in findLagrangeQP function");
       return 1;
   }

   if (q < 2 || q > 2){
       REprintf("Error in findLagrangeQP function");
       return 1;
   }

   Matrix<double> Ugeneral = U - H * dvec;

   int toReturn = findLagrange(Xi, Ugeneral, dCon);

   return(toReturn);
}


// ========================================================================================
// Function used to find optimal values of the Lagrangians according to
// the solution of general optimization with TWO EQUALITY constraints
//
// REMARK: At this moment implemented only for TWO EQUALITY CONSTRAINTS

// q(d) ..... objective function which was maximized

// Lagrange function is assumed to have a form
// L(d, xi) = q(d) - sum(xi * constraint),
// Optimal multipliers are thus a solution to
// sum(xi * d(constraint)) = U, where U = dq/dd

// !!! I do not do any checks of the input arguments
//     It is assumed to be used (properly) inside the smoothSurvReg function

// INPUT:  U ...... vector (p x 1) of the first derivatives of the objective function
//         dCon ... matrix (q x p) with derivatives of the constraints in ROWS

// OUTPUT: Xi ..... vector (q x 1) with resulting Lagrangians

// RETURN: 0 if no problems
//         >0 if there are problems
//         99 if there is no solution for optimal Lagrangians
//            (this should theoretically never happen...)
int
findLagrange(Matrix<double> & Xi,
             const Matrix<double> & U,
             const Matrix<double> & dCon)
{
   const double ZERO = 1e-8;

   int p = U.size();
   int q = dCon.rows();
   if (dCon.cols() != p || Xi.size() != q || p < q){
       REprintf("Error in findLagrange function");
       return 1;
   }

   if (q < 2 || q > 2){
       REprintf("Error in findLagrange function");
       return 1;
   }

   // Sort absolute values of derivatives at each row
   Matrix<double> absdCon = Matrix<double>(q, p, false);
   for (int i = 0; i < absdCon.size(); i++) absdCon[i] = fabs(dCon[i]);
   Matrix<int> orderdCon;
   Matrix<int> invorderdCon;
   Matrix<double> sortabsdCon;

   Matrix<int> orderRow = Matrix<int>(1, p, false);
   Matrix<int> invorderRow = Matrix<int>(1, p, false);
   Matrix<double> sortabsRow;

   for(int i = 0; i < q; i++){
      sortabsRow = sortOrder(absdCon(i, 0, i, p-1), orderRow, invorderRow);
      if (i == 0){
         orderdCon = orderRow;
         invorderdCon = invorderRow;
         sortabsdCon = sortabsRow;
      }
      else{
         orderdCon = rbind(orderdCon, orderRow);
         invorderdCon = rbind(invorderdCon, invorderRow);
         sortabsdCon = rbind(sortabsdCon, sortabsRow);
      }
   }

   // Find non-zero elements in dCon
   Matrix<bool> nonZero = Matrix<bool>(q, p, false);
   Matrix<int> nRowNonZero = Matrix<int>(q, 1, true, 0);
   Matrix<int> lastNonZero = Matrix<int>(q, 1, true, -1);
   for (int i = 0; i < q; i++)
      for (int j = 0; j < p; j++){
         nonZero(i, j) = fabs(dCon(i, j)) > ZERO;
         if (nonZero(i, j)){
            nRowNonZero[i]++;
            lastNonZero[i] = j;
         }
      }

   // If some row contains only zero's then the corresponding
   // Lagrangian may be whatever (e.g. zero)
   int q2 = q;    // number of constraints with non-zero derivatives
   for (int i = 0; i < q; i++)
      if (nRowNonZero[i] == 0){
         Xi[i] = 0.0;
         q2--;
      }

   // All constraints have zero derivatives, I return also zero lagrangians
   if (q2 == 0) return 0;

   // One non-zero d(constraint)
   int icon;      // index of the constraint which is non-zero
   if (q2 == 1){
      icon = (nRowNonZero[0] > 0) ? 0 : 1;
      int r = lastNonZero[icon];
      double rightS = U[r];
      Xi[icon] = rightS / dCon(icon, r);
      return 0;
   }

   // Both constraints have non-zero d(constraint)
   int icon0, icon1;
   if (q2 == 2){
      icon0 = 0;
      icon1 = 1;
      int r0 = invorderdCon(icon0, p - 1);
      int r1 = invorderdCon(icon1, p - 1);

   // Try to change r1 such that it's different from r0
      if (r0 == r1){
         for (int j = p - 1; j >= 0; j--){
            if (invorderdCon(icon1, j) == r0) continue;
            if (nonZero(icon1, invorderdCon(icon1, j))){
               r1 = invorderdCon(icon1, j);
               break;
            }
         }
      }

  // If still r1 = r0, try to change r0 such that it's different from r1
      if (r0 == r1){
         for (int j = p - 1; j >= 0; j--){
            if (invorderdCon(icon0, j) == r1) continue;
            if (nonZero(icon0, invorderdCon(icon0, j))){
               r0 = invorderdCon(icon0, j);
               break;
            }
         }
      }


   // If still r1 = r0, there is only one non-zero element at each row of dCon
   //  if possible, give equal lagrangians
      if (r0 == r1){
         double rightS = U[r0];
         double coef = dCon(icon0, r0) + dCon(icon1, r1);
         if (fabs(coef) > ZERO){   // equal lagrangians are for sure possible
            Xi[icon0] = rightS / coef;
            Xi[icon1] = Xi[icon0];
         }
         else{                     // give opposite lagrangians
            Xi[icon0] = rightS / (2 * dCon(icon0, r0));
            Xi[icon1] = -Xi[icon0];
         }
         return 0;
      }
   // r1 != r0 => I can solve 2 x 2 system of equation with non-zero diagonal
      else{
         double rightS0 = U[r0];
         double rightS1 = U[r1];
         double coef1 = dCon(icon1, r1) - ((dCon(icon0, r1) * dCon(icon1, r0))/dCon(icon0, r0));
         double right1 = rightS1 - ((dCon(icon0, r1) * rightS0)/dCon(icon0, r0));
         if (fabs(coef1) > ZERO){
            Xi[icon1] = right1 / coef1;
            Xi[icon0] = (1/dCon(icon0, r0)) * (rightS0 - dCon(icon1, r0) * Xi[icon1]);
            return 0;
         }
         else{     // singular system
            if (fabs(right1) < ZERO){
               Xi[icon1] = 0.0;
               Xi[icon0] = rightS0 / dCon(icon0, r0);
               return 0;
            }
            else{  // inconsistent system (this should theoretically never happen)
               return 99;
            }
         }
      }
   }
   return 0;
}


// ========================================================================================
// Function used to split a vector vec according to the logical vector cond
//   vecFalse has to have a length of length(cond) - sum(cond)
//   vecTrue has to have a length of sum(cond)
//   and the user has to know their length beforehand.

// INPUT:  vec ...... vector to be splitted
//         cond ..... logical vector of the same size as vec indicating how to split it

// OUTPUT: vecFalse ..... components of the vector vec corresponding to 'false' components of
//                        the vector cond
//         vecTrue ..... components of the vector vec corresponding to 'true' components of
//                       the vector cond

// RETURN: 0           if no problems
//         1, 2, 3, 4  if there were some problems to split the vector vec

int
splitVec(const Matrix<double> & vec, Matrix<double> & vecFalse, Matrix<double> & vecTrue, const Matrix<bool> cond)
{
   int n = vec.size();    // this must be also a size of a vector 'cond'
   if (n != cond.size()){
       REprintf("Error in a splitVec function");
       return 1;
   }

   int i, j, k;
   for (i = 0, j = 0, k = 0; i < n; i++){
       if (cond[i]){
          if (j < vecTrue.size()){
              vecTrue[j] = vec[i];
              j++;
          }
          else{
              REprintf("Too short 'vecTrue' vector in a splitVec function");
              return 2;
          }
       }
       else{
          if (k < vecFalse.size()){
              vecFalse[k] = vec[i];
              k++;
          }
          else{
              REprintf("Too short 'vecFalse' vector in a splitVec function");
              return 3;
          }
       }
   }

   if (j != vecTrue.size() || k != vecFalse.size()){
       REprintf("Too long 'vecTrue' or 'vecFalse' vector in a splitVec function");
       return 4;
   }

   return 0;

}   // end of splitVec function

