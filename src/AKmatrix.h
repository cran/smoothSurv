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
#ifndef AK_MATRIX_H
#define AK_MATRIX_H

#include "Scythe_Matrix.h"

namespace SCYTHE {

template <class T>
Matrix<T> sortOrder(const Matrix<T> &, Matrix<int> &, Matrix<int> &);

} // end namespace SCYTHE
#if defined (__GNUG__) || defined (__MWERKS__) || defined (_MSC_VER) || \
    defined (EXPLICIT_TEMPLATE_INSTANTIATION)
  // Necessary for template instantiation with some compilers.
# include "AKmatrix.cpp"
#endif  /* EXPLICIT_TEMPLATE_INSTANTIATION */

#endif
