// Helping functions for smoothSurvReg
// ====================================

// Various functions to convert a's to c's
// and vice versa

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

#ifndef CONVERT_C_A_D_H
#define CONVERT_C_A_D_H

#include <R.h>
#include <Rmath.h>

#include <cmath>
#include <cfloat>

#include "Scythe_Matrix.h"
#include "Scythe_Math.h"
#include "Scythe_Stat.h"
#include "Scythe_LA.h"

using namespace SCYTHE;

int
A_to_C(const Matrix<double> & A,
       Matrix<double> & C);

int
deriv_expAD(const Matrix<double> & knots,
            double sd0,
            const Matrix<int> & lastThree,
            Matrix<double> & Omega,
            bool all = true);

void
compute_dA2dD(Matrix<double> & dA2dD,
              const Matrix<double> & exp_Acoef2,
              const Matrix<double> & exp_Dcoef,
              const Matrix<double> & OmegaExpA2);

void
compute_ddA2dDD(Matrix<double> & ddAdDD1,
                Matrix<double> & ddAdDD2,
                const Matrix<double> & dA2dD,
                const int nD);

void
compute_dCdA(Matrix<double> & dCdA2,
             Matrix<double> & dCdAg3,
             const Matrix<double> & Ccoef,
             const Matrix<double> & tCcoef,
             const Matrix<int> & lastThreeA,
             const Matrix<int> & restA,
             const int nD);

void
compute_dCdD(Matrix<double> & dCdD,
             const Matrix<double> & dCdA2,
             const Matrix<double> & dCdAg3,
             const Matrix<double> & dA2dD,
             const Matrix<double> & Ccoef,
             const Matrix<double> & tCcoef,
             const Matrix<int> & restA,
             const int nD,
             const bool useD);

void
compute_ddCdAA(Matrix<double> * ddCdAA,
               const Matrix<double> & Ccoef,  
               const int posZeroA);

void
compute_ddCdDD(Matrix<double> * ddCdDD,
               const Matrix<double> & dCdA2,
               const Matrix<double> * ddCdAA,
               const Matrix<double> & dA2dD,
               const Matrix<double> & ddAdDD1,
               const Matrix<double> & ddAdDD2,
               const Matrix<int> & lastTwoA,
               const Matrix<int> & restA);

#endif
