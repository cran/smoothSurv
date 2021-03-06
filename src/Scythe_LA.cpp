/* Scythe_LA.cpp
 *
 * This file provides implementations of functions which deal with
 * linear algebra in the Scythe Statistical Library.
 * 
 * Scythe C++ Library
 * Copyright (C) Kevin M. Quinn, Andrew D. Martin,
 * and Daniel B. Pemstein
 *
 * This code written by:
 *
 * Kevin Quinn
 * Assistant Professor
 * Dept. of Political Science and
 * Center for Statistics and Social Sciences
 * Box 354322
 * University of Washington
 * Seattle, WA 98195-4322
 * quinn@stat.washington.edu
 *
 * Andrew D. Martin
 * Assistant Professor
 * Dept. of Political Science
 * Campus Box 1063
 * Washington University
 * St. Louis, MO 63130
 * admartin@artsci.wustl.edu
 * 
 * Daniel B. Pemstein
 * dbpemste@artsci.wustl.edu
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 * USA
 */

#ifndef SCYTHE_LA_CC
#define SCYTHE_LA_CC

#include <cmath>
#include <algorithm>
#include <set>
#include "Scythe_Error.h"
#include "Scythe_Util.h"
#include "Scythe_LA.h"
#include "Scythe_Stat.h"


namespace SCYTHE {
	/* Compute the transpose of a matrix. Kept for back-compatibility*/
	template <class T>
	Matrix<T>
	t (const Matrix<T> &old_matrix)
	{
		return (! old_matrix);
	}
	
	/* Create a matrix of ones from the given dimensions 
	 * Note:  call is of from ones<double>(4,3) or ones<int>(5,8)
	 */
	template <class T>
	Matrix<T>
	ones (const int& rows, const int& cols)
	{
  	if (rows < 1 || cols < 1) {
			throw scythe_dimension_error(__FILE__, __AK_PRETTY_FUNCTION__,
					__LINE__, std::string("Improper row (") & rows
					& ") or column (" & cols & ") dimension");
  	}
		
		return Matrix<T> (rows, cols, true, (T) 1);
	}
	
	/* Create a k x k identity matrix
	 * Note:  class is of form eye<double>(4) or eye<int>(7
	 */
	template <class T>
	Matrix<T>
	eye (const int &k)
	{
		Matrix<T> temp(k, k, false);
		for (int i = 0; i < temp.rows(); ++i) {
			for (int j = 0; j < temp.cols(); ++j) {
				if (i == j)
					temp(i,j) = (T) 1.0;
				else
					temp(i,j) = (T) 0.0;
			}
		}
		
  	return temp;
	}
	
	/* Create a k x 1 vector-additive sequence matrix */
	template <class T>
	Matrix<T>
	seqa (T start, const T& incr, const int& size)
	{
		Matrix<T> temp (size, 1, false);
		for (int i = 0; i < size; ++i) {
			temp[i] = start;
			start += incr;
		}
 	 
		return temp;
	}
	
	/* Uses the STL sort to sort a Matrix in ascending row-major order */
	template <class T>
	Matrix<T>
	sort (Matrix<T> A) {
		sort(A.begin(), A.end());
		return A;
	}
		
	template <class T>
	Matrix<T>
	sortc (Matrix<T> A)
	{
	  for (typename SCYTHE::Matrix<T>::col_major_Iterator it = A.beginc(); it < A.endc(); it.next_vec())
		sort(it, it + A.rows());

  	  return A;
	}
	
	/* Column bind two matrices */
	template <class T>
	Matrix<T>
	cbind (const Matrix<T> &A, const Matrix<T> &B)
	{
		if (A.rows() != B.rows()) {
			throw scythe_conformation_error(__FILE__, __AK_PRETTY_FUNCTION__,
					__LINE__, "Matrices have different number of rows");
		}
		
		Matrix<T> C(A.rows(), A.cols() + B.cols(), false);

		typename SCYTHE::Matrix<T>::col_major_Iterator write = C.beginc();

        	for (typename SCYTHE::Matrix<T>::const_col_major_Iterator read = A.beginc(); read < A.endc(); ++read)
			*(write++) = *read;
	
		for (typename SCYTHE::Matrix<T>::const_col_major_Iterator read = B.beginc(); read < B.endc(); ++read)
			*(write++) = *read;
	
		return C;
	}

	
	/* Row bind two matrices: kept for backwards compatibility */
	template <class T>
	Matrix<T>
	rbind (const Matrix<T> &A, const Matrix<T> &B)
	{
		if (A.cols() != B.cols()) {
			throw scythe_conformation_error(__FILE__, __AK_PRETTY_FUNCTION__,
					__LINE__, "Matrices have different number of rows");
		}
		
		Matrix<T> C(A.rows() + B.rows(), A.cols(), false);
		typename SCYTHE::Matrix<T>::row_major_Iterator write = C.begin();
	
		for (typename SCYTHE::Matrix<T>::const_row_major_Iterator read = A.begin(); read < A.end(); ++read)
			*(write++) = *read;
	
		for (typename SCYTHE::Matrix<T>::const_row_major_Iterator read = B.begin(); read < B.end(); ++read)
			*(write++) = *read;
	
		return C;
	}
	
	/* Calculates the order of each element in a Matrix */
	template <class T>
	Matrix<int>
	order(const Matrix<T> &A){
		if (! A.isColVector()) {
			throw scythe_dimension_error (__FILE__, __AK_PRETTY_FUNCTION__,
					__LINE__, "A not a column vector");
		}
		Matrix<int> temp(A.rows(), 1, false);
		for (int i = 0;  i < A.rows(); ++i) {
			temp[i] = sumc(A << A[i])[0];
		}
		
		return temp;
	}
	
	/* Selects all the rows of Matrix A for which binary column vector e
 	* has an element equal to 1
 	*/
	template <class T>
	Matrix<T>
	selif(const Matrix<T> &A, const Matrix<bool> &e)
	{
		if (A.rows() != e.rows()) {
			throw scythe_dimension_error (__FILE__, __AK_PRETTY_FUNCTION__,
					__LINE__, "A and e have different number of rows");
		}
	
		if (! e.isColVector()) {
			throw scythe_dimension_error (__FILE__, __AK_PRETTY_FUNCTION__,
					__LINE__, "e not a column vector");
		}
 	 
		// See how many rows are true
		int N = std::accumulate(e.begin(), e.end(), (int) 0);
 	 
		// declare and form output Matrix
		Matrix<T> temp(N, A.cols(), false);
		int cnt = 0;
		for (int i = 0; i < e.size(); ++i) {
			if (e[i])  {
				copy(A.vec(i), A.vec(i + 1), temp.vec(cnt++));
			}
		}
	
		return temp;
	}
	
	/* Find unique elements in a matrix and return a sorted row vector */
	template <class T>
	Matrix<T>
	unique(const Matrix<T> &A)
	{
		std::set<T> u(A.begin(), A.end());
		Matrix<T> temp(1, u.size(), false);
		
		copy(u.begin(), u.end(), temp.begin());

		return temp;
	}

	/* Reshape a matrix */
	template <class T>
	Matrix<T>
	reshape(const Matrix<T> &A, const int &r, const int &c) 
	{
		if (A.size() != r * c)
			throw scythe_invalid_arg(__FILE__, __AK_PRETTY_FUNCTION__, __LINE__,
					std::string("Input dimensions (") & r & "," & c & ") not" &
					" consistent with size of input matrix (" & A.size() & ")");
	
		Matrix<T> temp(r, c, A.getArray());
  	return temp;
	}
	
	/* Make vector out of unique elements of a symmetric Matrix.  
 	* NOTE: DOES NOT CHECK FOR SYMMETRY!!!
 	*/
	template <class T>
	Matrix<T>
	vech(const Matrix<T> &A)
	{
		if (! A.isSquare()) {
			throw scythe_dimension_error(__FILE__, __AK_PRETTY_FUNCTION__,
					__LINE__, "Matrix not square");
  	}
		Matrix<T> temp ((int) (0.5 * (A.size() - A.rows())) + A.rows(), 1,
				false);
		typename SCYTHE::Matrix<T>::row_major_Iterator iter = temp.begin();
		
		for (int i = 0; i < A.rows(); ++i)
			iter = copy(A.vecc(i) + i, A.vecc(i + 1), iter);
		
 	 
		return temp;
	}
	
	/* Expand xpnd(A) == B from A = vech(B) */
	template <class T>
	Matrix<T>
	xpnd(const Matrix<T> &A)
	{
	  double newrowsize_d = -.5 + .5 * ::sqrt(1 + 8 * A.size());      
	  if (std::fmod(newrowsize_d, 1.0) != 0.0)
			throw scythe_invalid_arg(__FILE__, __AK_PRETTY_FUNCTION__, __LINE__,
					"Can't turn input vector into a square matrix");
		
		int newrowsize = (int) newrowsize_d;
		Matrix<T> temp(newrowsize, newrowsize, false);
		int cnt = 0;
	
		for (int i = 0; i < newrowsize; ++i) {
			for (int j = i; j < newrowsize; ++j)
				temp(i, j) = temp(j, i) = A[cnt++];
		}
	
		return temp;
	}
	
	/* Get the diagonal of a Matrix. */
	template <class T>
	Matrix<T>
	diag(const Matrix<T> &A)
	{
		if (A.rows() != A.cols())
			throw scythe_dimension_error(__FILE__, __AK_PRETTY_FUNCTION__,
					__LINE__, "Matrix not square");
	
		Matrix<T> temp(A.rows(), 1, false);
		for (int i = 0; i < A.rows(); ++i)
			temp[i] = A(i, i);
 	 
		return temp;
	}

	template <class T>
	Matrix<T>
	gaxpy (const Matrix<T> &A, const Matrix<T> &B, const Matrix<T> &C)
	{
		Matrix<T> temp;

		if (A.isScalar() && B.rows() == C.rows() && B.cols() == C.cols()) {
			// Case 1: 1 x 1  *  n x k  +  n x k
			temp = Matrix<T> (B.rows(), B.cols(), false);
			
			for (int i = 0; i  < B.size(); ++i)
				temp[i] = A[0] * B[i] + C[i];
			
		} else if (B.isScalar() && A.rows() == C.rows() && 
				A.cols() == C.cols()) {
			// Case 2: m x n  *  1 x 1  +  m x n
			temp = Matrix<T> (A.rows(), A.cols(), false);
			
			for (int i = 0; i  < B.size(); ++i)
				temp[i] = A[i] * B[0] + C[i];

		} else if (A.cols() == B.rows() && A.rows() == C.rows() &&
				B.cols() == C.cols()) {
			// Case 3: m x n  *  n x k  +  m x n
			temp = Matrix<T> (A.rows(), B.cols(), false);
			
			for (int i = 0; i < A.rows(); ++i) {
				for (int j = 0; j < B.cols(); ++j) {
					temp[i * B.cols() + j] = C[i * B.cols() +j];
					for (int k = 0; k < B.rows(); ++k)
						temp[i * B.cols() + j] += A[i * A.cols() + k] * 
							B[k * B.cols() + j];
				}
			}

		} else {
			throw scythe_conformation_error(__FILE__, __AK_PRETTY_FUNCTION__,
					__LINE__, std::string("Expects (m x n  *  1 x 1  +  m x n)") &
						"or (1 x 1  *  n x k  +  n x k) or (m x n  *  n x k  +" &
						"  m x k");
		}

		return temp;
	}

	/* Fast calculation of A'A */
	template <class T>
	Matrix<T>
	crossprod (const Matrix<T> &A)
	{
		Matrix<T> temp(A.cols(), A.cols(), false);
	
		for (int i = 0; i < A.cols(); ++i) {
			for (int j = i; j < A.cols(); ++j) {
				temp(i,j) = T (0);
				for (int k = 0; k < A.rows(); ++k)
					temp(j,i) = temp(i,j) += A(k,i) * A(k,j);
			}
		}
	
		return temp;
	}

} // end namespace SCYTHE
#endif /* SCYTHE_LA_CC */
