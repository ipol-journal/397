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


#include "image2D.h"
#include "mysparsematrix.h"

/// Computes the dilation as max-plus matrix-to-vector product
///            im_out[j] = max_i {im_in[i] + M[i,j]}
/// where the matrix M is supposed to be symmetric and is
/// represented by the three arrays I, J, W:
///          M[I[l], J[l]] = W[l] and zero for other indices.
int dilate_max_plus_symmetric(size_t* I, size_t* J, double* W,
			      size_t n_lines,
			      const Image2D<float> &imIn,
			      Image2D<float> &imOut);

/// Computes iterations of the symmetric max-plus dilation.
int dilate_max_plus_symmetric_iterated(size_t* I, size_t* J, double* W,
				      int k,
				      size_t n_lines,
				      const Image2D<float> &imIn,
				      Image2D<float> &imOut);

/// Computes the erosion as min-plus matrix-to-vector product
///            im_out[j] = min_i {im_in[i] - M[j,i]}
/// where the matrix M is supposed to be symmetric and is
/// represented by the three arrays I, J, W:
///          M[I[l], J[l]] = W[l] and zero for other indices.
int erode_max_plus_symmetric(size_t* I, size_t* J, double* W,
			     size_t n_lines,
			     const Image2D<float> &imIn,
			     Image2D<float> &imOut);


/// Computes iterations of the symmetric min-plus erosion.
int erode_max_plus_symmetric_iterated(size_t* I, size_t* J, double* W,
				      int k,
				      size_t n_lines,
				      const Image2D<float> &imIn,
				      Image2D<float> &imOut);


/// Computes the opening based on adjoint iterated symmetric max-plus
/// dilations and min-plus erosions
int open_max_plus_symmetric(size_t* I, size_t* J, double* Weights,
			    int it,
			    size_t n_non_zeros,
			    const Image2D<float> &imIn,
			    Image2D<float> &imOut);

/// Computes the closing based on adjoint iterated symmetric max-plus
/// dilations and min-plus erosions
int close_max_plus_symmetric(size_t* I, size_t* J, double* Weights,
			     int it,
			     size_t n_non_zeros,
			     const Image2D<float> &imIn,
			     Image2D<float> &imOut);
