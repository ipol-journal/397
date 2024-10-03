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

#ifndef STRUCTURE_TENSORS_HPP
#define STRUCTURE_TENSORS_HPP
#include "image2D.h"

/// Linear gradient of a single channel image im
/// gx and gy are approximations of the partial derivatives
/// For the pixel at raw i and column j:
///         gx(i,j) = ( im(i, j+1) - im(i, j-1) )/2
///         gy(i,j) = ( im(i+1, j) - im(i-1, j) )/2
/// when this makes sense. 
/// On the image's edges (i=0, i=H-1, j=0, j=W-1), the
/// approximation is then just the difference betwen two
/// adjacent pixels, for example
///         gx(i,0) = im(i, 1) - im(i, 0)
void lin_grad_gray(const Image2D<float> &im, double* gx, double* gy);

/// Compute the rank 1 structure tensors T0 from the gradient of the
/// image represented by its x and y components, contained in the two
/// arrays gx and gy
///
/// The returned array of symmetrical structure tensors T0 is
/// represented by three arrays of coefficients a, b, and c, such
/// that      
///            | a  c |
///       T0 = | c  b |
///
void structure_tensors0(double* gx, double* gy,
		       size_t H, size_t W,
		       double* a, double* b, double* c);

/// Define the 1D Gaussian kernel with variance sigma^2, cropped in a
/// 2*radius+1 large window centered on zero.
///
/// The kernel is NOT normalized and does NOT sum to one.
void gaussian_kernel(double* kernel, int radius, double sigma);

/// 2D Convolution with separable kernel
void convol(const double* imIn, double* kernel,
	    size_t H, size_t W, int radius,
	    double* imOut);


/// Gaussian smoothing of a 2D grayscale image
void gaussian_filter(const double* imIn, double sigma, size_t H, size_t W,
		     double* imOut);

/// Gaussian smoothing of structure tensors T_sigma, component by
/// component. The regularized symmetric structure tensors, one for
/// each pixel, are represented by the three arrays a_reg, b_reg,
/// c_reg, each one of length W*H
///
///                 | a_reg  c_reg |
///       T_sigma = | c_reg  b_reg |
///
/// T_sigma is computed from the input zero-scale structure tensors
///
///            | a  c |
///       T0 = | c  b |
///
/// represented by the three arays a, b and c.
void smooth_structure_tensors(const double* a, const double* b, const double* c,
			      double sigma, size_t H, size_t W,
			      double* a_reg,double* b_reg, double* c_reg);

/// Diagonalizes the 2 x 2 symmetric positive definite matrices, one
/// for each pixel, represented by the three arrays a, b and c.
///           
///           | a  c |
///       T = | c  b |
///
/// Returns arrays of eigenvalues (lambda1 and lambda2), eigenvectors
/// (eigvect2_x and eigvect2_y) and the anisotropy values
/// (anisotropy)
///
void diagonalize_structure_tensors(const double* a,
				   const double* b,
				   const double* c,
				   size_t W,
				   size_t H,
				   double* lambda1,
				   double* lambda2,
				   double* eigvect2_x,
				   double* eigvect2_y,
				   double* anisotropy);

#endif // STRUCTURE_TENSORS_HPP
