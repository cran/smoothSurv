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

#ifndef PENAL_LOG_LIK_H_
#define PENAL_LOG_LIK_H_

#include <R.h>
#include <Rmath.h>

#include <iostream>
#include <cmath>
#include <cfloat>

#include "gauss.h"
#include "convertCAD.h"

#include "AKmatrix.cpp"
#include "Scythe_Matrix.h"
#include "Scythe_Math.h"
#include "Scythe_Stat.h"
#include "Scythe_IDE.h"
#include "Scythe_LA.h"

using namespace SCYTHE;

double
penalLogLik(const Matrix<double> & Theta,
            const int what = 3,
            const int derivOrder = 2,
            const bool useD = true);

int
findLagrangeQP(Matrix<double> & Xi,
               const Matrix<double> & U,
               const Matrix<double> & H,
               const Matrix<double> & dvec,
               const Matrix<double> & dCon);

int
findLagrange(Matrix<double> & Xi,
             const Matrix<double> & U,
             const Matrix<double> & dCon);

int
splitVec(const Matrix<double> & vec,
         Matrix<double> & vecFalse,
         Matrix<double> & vecTrue,
         const Matrix<bool> cond);

#endif
