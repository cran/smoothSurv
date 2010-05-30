// Gaussian density and survival function

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
// 30/05/2010:  erf, erfc changed to erfAK, erfcAK
//              becaues of some conflicts with /usr/include/bits/mathcalls.h
//              reported by CRAN checks
//

#ifndef GAUSS_H
#define GAUSS_H

#include "Scythe_Matrix.h"
#include <R.h>
#include <Rmath.h>

const double SPI = 2.506628274631001;                     // sqrt(2*pi)
const double ROOT_2 = 1.414213562373095048801688724210;   // sqrt(2)
const double ZERO = 1e-50;

using namespace SCYTHE;

// Return 2 * Phi(x * sqrt(2)) - 1
inline double erfAK(double x)
{
  return 2*pnorm(x*ROOT_2, 0, 1, 1, 0) - 1;
}

// Return 2 * Phi(-x * sqrt(2))
inline double erfcAK(double x)
{
  return 2*pnorm(-x*ROOT_2, 0, 1, 1, 0);
}

double fnorm(double);
Matrix<double> fnorm(Matrix<double> &);
double Snorm(double);
Matrix<double> Snorm(Matrix<double> &);
double Fnorm(double);
Matrix<double> Fnorm(Matrix<double> &);

double fnormZero(double);
Matrix<double> fnormZero(Matrix<double> &);
double SnormZero(double);
Matrix<double> SnormZero(Matrix<double> &);
double FnormZero(double);
Matrix<double> FnormZero(Matrix<double> &);

#endif
