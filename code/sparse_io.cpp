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

#include "sparse_io.hpp"

/// Count the number of non-zero coefficients in the sparse matrix S
/// containing n_col colmuns
size_t count_sparse_matrix_non_zeros(MySparseMatrix<double> &S, size_t n_col)
{
  size_t n_non_zeros = 0;
  size_t j;
  std::map< unsigned int , double > col_j;
  for(j=0;j<n_col;j++)
    {
      col_j = S.get_column(j);
      for( std::map< unsigned int , double >::const_iterator it = col_j.begin() ; it != col_j.end() ; ++it )
	n_non_zeros++;
    }
  return n_non_zeros;
}

/// Translates the sparse matrix S counting n_col columns, into three
/// arrays I, J, W such that S[I[l], J[l]] = W[l]
void sparse_matrix_to_arrays(MySparseMatrix<double> &S, size_t* I, size_t* J, double* W, size_t n_col)
{
  size_t index = -1;
  std::map< unsigned int , double > col_j;
  size_t j;
  for(j=0;j<n_col;j++)
    {
      col_j = S.get_column(j);
      for( std::map< unsigned int , double >::const_iterator it = col_j.begin() ; it != col_j.end() ; ++it ) 
	{
	  index++;
	  J[index] = j;
	  I[index] = it->first;
	  W[index] = it->second;
	}
    }
}
