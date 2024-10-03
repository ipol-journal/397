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

#ifndef COCIRCULAIRITY_GRAPHS_HPP
#define COCIRCULAIRITY_GRAPHS_HPP

#include "structure_tensors.hpp"
#include "mysparsematrix.h"

/// Stores a list of pixel coordinates in row_list and col_list. The
/// corresponding set of pixels is composed of two rectangles,
/// respectively bounded by [i_min1, i_max1]x[j_min1, j_max1] and
/// [i_min2, i_max2]x[j_min2, j_max2]
void pixel_list(size_t* row_list, size_t* col_list,
		size_t i_min1, size_t i_min2,
		size_t i_max1, size_t i_max2,
		size_t j_min1, size_t j_min2,
		size_t j_max1, size_t j_max2);

/// Fill in the sparse, symmetric matrix S representing the
/// co-circularity graph. Compared to the paper, here we actually
/// compute the exponential of the weights (they are non negative and
/// belong to [0,1]), and the log will be applied when applying
/// morphological max-plus/mn-plus operators.
///
/// S : sparse, symmetric matrix to be filled in
///
/// p : integer ruling the size of the (2p+1)x(2p+1) square
/// neighbourhood
///
/// W : width of the underlying image
///
/// H : height of the underlying image
///
/// eigvect2_x : array of size W*H containing the x-coordinate of the
/// pricipal direction of the structure tensor for each pixel
///
/// eigvect2_y : array of size W*H containing the y-coordinate of the
/// pricipal direction of the structure tensor for each pixel
///
/// thresh_cocirc : cos(alpha), with alpha the angular precision on
/// co-circularity (only useful for binary weighted graph, i.e. if
/// bin=1)
///
/// thresh_cone : cos(beta), with beta the angle which determines the
/// conic condition.
///
/// bin : 1 for binary-weighted graph, 0 otherwise.
void compute_adjacency_matrix(MySparseMatrix<double> &S, 
			      size_t p, size_t W, size_t H, 
			      const double* eigvect2_x,
			      const double* eigvect2_y,
			      double thresh_cocirc,
			      double thresh_cone,
			      int bin);


/// Determine whether point 1 (x1, y1) with direction (v1x, v1y) and
/// point 2 (x2, y2) with direction (v2x, v2y) are co-circular up to
/// the angular precision thresh_cocirc (binary case), or to which
/// degree they are co-circular (non-binary case) and if the conic
/// constraint is satisfied. The answer resp is in [0, 1] and we will
/// take 255*log(resp) as weight between the two corresponding
/// vertices.
double test_adjacency(double x1, double y1, double v1x, double v1y,
		      double x2, double y2, double v2x, double v2y, 
		      double thresh_cocirc, double thresh_cone,
		      int bin);

/// From the input image im, compute the co-circularity graph by
/// filling in the sparse, symmetric matrix S representing the
/// co-circularity graph; it also returns the anisotropy image
/// corresponding to parameter sigma.
///
/// im : input image
///
/// S : sparse, symmetric matrix to be filled in
///
/// im_anisotropy : anisotropy image for the chosen sigma; to be
/// filled in
///
/// sigma : smoothing parameter
///
/// p : integer ruling the size of the (2p+1)x(2p+1) square
/// neighbourhood
///
/// thresh_cocirc : cos(alpha), with alpha the angular precision on
/// co-circularity (only useful for binary weighted graph, i.e. if
/// mode_str=1)
///
/// thresh_cone : cos(beta), with beta the angle which determines the
/// conic condition.
///
/// mode_str : bin for binary-weighted graph, nbin otherwise.
void image_to_cocircularity_graph(const Image2D<float> &im, MySparseMatrix<double> &S, 
				  Image2D<float> &im_anisotropy,
				  double sigma, size_t p, 
				  double thresh_cocirc, double thresh_cone,
				  string mode_str);


#endif // COCIRCULAIRITY_GRAPHS_HPP
