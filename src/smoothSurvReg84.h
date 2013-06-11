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

#ifndef smooth_Surv84_S_H
#define smooth_Surv84_S_H

#include <R.h>
#include <Rmath.h>

#include <iostream>
#include <cmath>
#include <cfloat>

#include "penalLogLik.h"
#include "gauss.h"
#include "difference.h"
#include "convertCAD.h"
#include "solve.QP.compact.h"
#include "linpackFORT.h"
#include "createPosDef.h"
#include "eispackFORT.h"
#include "AKmatrix.h"

#include "Scythe_Matrix.h"
#include "Scythe_Math.h"
#include "Scythe_Stat.h"
#include "Scythe_IDE.h"
#include "Scythe_LA.h"

const int non_conv_flag = 1000;          // flag for non-convergence

using namespace SCYTHE;

// Return Euclidian norm of the vector
inline 
double EucNorm(Matrix<double> & A)
{
  double Asq = sum(A & A);
  if (Asq > FLT_MAX - 1) return FLT_MAX;
  else                   return sqrt(Asq);
}

#endif
