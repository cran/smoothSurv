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

#include "gauss.h"

using namespace SCYTHE;

// Density of standard normal distribution (scalar version).
double
fnorm(double z)
{
  double density = exp(-z*z/2)/SPI;
  return density;
}

// Density of standard normal distribution (matrix version).
Matrix<double>
fnorm(Matrix<double> &z)
{
  int nr = z.rows();
  int nc = z.cols();
  int nrc = nr * nc;
  double *densf = new double[nrc];
  for(int i = 0; i < nrc; i++){
     densf[i] = fnorm(z[i]);
  }
  Matrix<double> Densf(nr, nc, densf);
  delete [] densf;
  return Densf;
}

// Survivor function of standard normal distribution (scalar version).
double
Snorm(double z)
{
     double survf;
     if (z > 0) {
       survf =  erfcAK(z/ROOT_2) /2;
     }
     else {
       survf = (1 + erfAK(-z/ROOT_2))/2;
     }

     return survf;
}

// Survivor function of standard normal distribution (matrix version).
Matrix<double>
Snorm(Matrix<double> &z)
{
  int nr = z.rows();
  int nc = z.cols();
  int nrc = nr * nc;
  double *survf = new double[nrc];
  for(int i = 0; i < nrc; i++){
     survf[i] = Snorm(z[i]);
  }
  Matrix<double> Survf(nr, nc, survf);
  delete [] survf;
  return Survf;
}

// CDF of standard normal distribution (scalar version).
double
Fnorm(double z)
{
     double cdff;
     if (z > 0) {
       cdff = (1 + erfAK(z/ROOT_2))/2;
     }
     else {
       cdff = erfcAK(-z/ROOT_2) /2;
     }

     return cdff;
}

// CDF of standard normal distribution (matrix version).
Matrix<double>
Fnorm(Matrix<double> &z)
{
  int nr = z.rows();
  int nc = z.cols();
  int nrc = nr * nc;
  double *cdff = new double[nrc];
  for(int i = 0; i < nrc; i++){
     cdff[i] = Fnorm(z[i]);
  }
  Matrix<double> Cdff(nr, nc, cdff);
  delete [] cdff;
  return Cdff;
}

// VERSIONS OF ABOVE MENTIONED FUNCTIONS WHICH CHANGE COMPUTER ZEROS INTO
// ZERO DEFINED BY ME
// ======================================================================
// Density of standard normal distribution (scalar version).
double
fnormZero(double z)
{
  double density = exp(-z*z/2)/SPI;
  if (density < ZERO) density = ZERO;
  return density;
}

// Density of standard normal distribution (matrix version).
Matrix<double>
fnormZero(Matrix<double> &z)
{
  int nr = z.rows();
  int nc = z.cols();
  int nrc = nr * nc;
  double *densf = new double[nrc];
  for(int i = 0; i < nrc; i++){
     densf[i] = fnormZero(z[i]);
  }
  Matrix<double> Densf(nr, nc, densf);
  delete [] densf;
  return Densf;
}

// Survivor function of standard normal distribution (scalar version).
double
SnormZero(double z)
{
     double survf;
     if (z > 0) {
       survf =  erfcAK(z/ROOT_2) /2;
     }
     else {
       survf = (1 + erfAK(-z/ROOT_2))/2;
     }
     if (survf < ZERO) survf = ZERO;
     return survf;
}

// Survivor function of standard normal distribution (matrix version).
Matrix<double>
SnormZero(Matrix<double> &z)
{
  int nr = z.rows();
  int nc = z.cols();
  int nrc = nr * nc;
  double *survf = new double[nrc];
  for(int i = 0; i < nrc; i++){
     survf[i] = SnormZero(z[i]);
  }
  Matrix<double> Survf(nr, nc, survf);
  delete [] survf;
  return Survf;
}

// CDF of standard normal distribution (scalar version).
double
FnormZero(double z)
{
     double cdff;
     if (z > 0) {
       cdff = (1 + erfAK(z/ROOT_2))/2;
     }
     else {
       cdff = erfcAK(-z/ROOT_2) /2;
     }
     if (cdff < ZERO) cdff = ZERO;
     return cdff;
}

// CDF of standard normal distribution (matrix version).
Matrix<double>
FnormZero(Matrix<double> &z)
{
  int nr = z.rows();
  int nc = z.cols();
  int nrc = nr * nc;
  double *cdff = new double[nrc];
  for(int i = 0; i < nrc; i++){
     cdff[i] = FnormZero(z[i]);
  }
  Matrix<double> Cdff(nr, nc, cdff);
  delete [] cdff;
  return Cdff;
}



