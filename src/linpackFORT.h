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

#ifndef LINPACK_FORT_H_
#define LINPACK_FORT_H_

#include <R.h>
#include <Rmath.h>
#include <iostream>

double
ddotCPP(const int n, double* dx, const int incx, double* dy, const int incy);

void
daxpyCPP(const int n, const double da, double* dx, const int incx, double* dy, const int incy);

void
dscalCPP(const int n, const double da, double* dx, const int incx);

void
dpofaCPP(double* a, const int lda, const int n, int* iinfo, const double eps = 1e-14);

void
dporiCPP(double* a, const int lda, const int n);

void
dposlCPP(double* a, const int lda, const int n, double* b);

#endif

