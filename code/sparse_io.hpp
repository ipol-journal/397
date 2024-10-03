/*
  Copyright (c) 2022 Samy Blusseau <samy.blusseau@minesparis.psl.eu>

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as
  published by the Free Software Foundation, either version 3 of the
  License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public
  License (GNU AGPL) along with this program. If not, see
  <http://www.gnu.org/licenses/>.
*/

#ifndef SPARSE_IO_HPP
#define SPARSE_IO_HPP
#include "mysparsematrix.h"

using namespace std;

size_t count_sparse_matrix_non_zeros(MySparseMatrix<double> &S,
				     size_t n_col);

void sparse_matrix_to_arrays(MySparseMatrix<double> &S, 
			     size_t* I, size_t* J, double* W,
			     size_t n_col);

#endif // SPARSE_IO_HPP
