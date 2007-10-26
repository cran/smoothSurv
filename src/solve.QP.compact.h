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

#ifndef SOLVE_QP_COMPACT_H_
#define SOLVE_QP_COMPACT_H_

#include <R.h>
#include <Rmath.h>

#include <iostream>
#include <cmath>

#include "linpackFORT.h"

extern "C" {

void
qpgen1CPP (double* dmat, double* dvec, int* fddmat, int* n,
           double* sol, double* crval,
           double* amat, int* iamat, double* bvec,
           int* fdamat, int* q, int* meq,
           int* iact, int* nact, int* iter, double* work, int* ierr, double* eps);

}

#endif

