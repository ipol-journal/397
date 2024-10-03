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

#include <cmath>
#include "structure_tensors.hpp"

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
void lin_grad_gray(const Image2D<float> &im, double* gx, double* gy)
{ 
  size_t W = im.width();
  size_t H = im.height();
  
  size_t i, j;
  // Compute gx
  for (i=0; i<H; i++)
    {
      // Case j=0
      gx[i*W] = im[i*W+1] - im[i*W];
      // Case j=W-1
      gx[i*W + W-1] = im[i*W + W-1] - im[i*W + W-2];
      for (j=1;j<W-1;j++)
	gx[i*W + j] =  0.5*(im[i*W + j+1] - im[i*W + j-1]);
    }
  // Compute gy
  for (j=0; j<W; j++)
    {
      // Case i=0
      gy[j] = im[W+j] - im[j];
      // Case i=H-1
      gy[(H-1)*W + j] = im[(H-1)*W + j] - im[(H-2)*W + j];
      for (i=1; i<H-1; i++)
	gy[i*W + j] =  0.5*(im[(i+1)*W + j] - im[(i-1)*W + j]);
    }
}


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
			double* a, double* b, double* c)
{
  size_t i, j;
  double gx_ij, gy_ij;
  for (i=0; i<H; i++)
    for (j=0; j<W; j++)
      {
	gx_ij = gx[i*W + j];
	gy_ij = gy[i*W + j];
	a[i*W + j] = gx_ij*gx_ij;
	b[i*W + j] = gy_ij*gy_ij;
	c[i*W + j] = gx_ij*gy_ij;
      }
}


/// Define the 1D Gaussian kernel with variance sigma^2, cropped in a
/// 2*radius+1 large window centered on zero.
///
/// The kernel is NOT normalized and does NOT sum to one.
void gaussian_kernel(double* kernel, int radius, double sigma)
{ 
  for (int i = -radius; i <= radius; i++)
    kernel[i + radius] = exp(-(pow(i / sigma, 2) / 2));
}


/// 2D Convolution with separable kernel
void convol(const double* imIn, double* kernel,
	    size_t H, size_t W, int radius,
	    double* imOut)
{ 
  double* imOutX = new double[H*W];

  /*
   * convolution in X
   */
  for (off_t y = 0; y < H; y++) {
    for (off_t x = 0; x < W; x++) {
      off_t i0    = y * W + x;
      double sumV = 0.;
      double sumK = 0.;
      
      for (off_t i = -radius; i <= radius; i++) {
	if ((x + i < 0) || (x + i > W - 1))
	  continue;
	double valK = kernel[i + radius];
	sumV += imIn[i0 + i] * valK;
	sumK += valK;
      }
      imOutX[i0] = sumV / sumK;
    }
  }
  /*
   * convolution in Y
   */
  off_t stride = W;
  for (off_t y = 0; y < H; y++) {
    for (off_t x = 0; x < W; x++) {
      off_t i0    =  y * W + x;
      double sumV = 0.;
      double sumK = 0.;
      
      for (off_t i = -radius; i <= radius; i++) {
	if ((y + i < 0) || (y + i > H - 1))
	  continue;
	double valK = kernel[i + radius];
	sumV += imOutX[i0 + i * stride] * valK;
	sumK += valK;
      }
      imOut[i0] = sumV / sumK;
    }
  }
  
  delete[] imOutX;
}


/// Gaussian smoothing of a 2D grayscale image
void gaussian_filter(const double* imIn, double sigma, size_t H, size_t W,
		     double* imOut){
  
  int radius = ceil(2*sigma);
  double* kernel = new double[2*radius+1];
  gaussian_kernel(kernel, radius, sigma);
  convol(imIn, kernel, H, W, radius, imOut);

  delete[] kernel;
}


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
			      double* a_reg,double* b_reg, double* c_reg)
{ 
  gaussian_filter(a, sigma, H, W, a_reg);
  gaussian_filter(b, sigma, H, W, b_reg);
  gaussian_filter(c, sigma, H, W, c_reg);
}


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
				   double* anisotropy)
{
  size_t i, j;
  double delta, trace, x1, x2, norm_eig_vect;
  for (i=0; i<H; i++)
    for (j=0; j<W; j++)
      {
	delta = (a[i*W + j]-b[i*W + j])*(a[i*W + j]-b[i*W + j]) + 4*c[i*W + j]*c[i*W + j];
	trace = a[i*W + j]+b[i*W + j];
	lambda1[i*W + j] = 0.5*(trace + sqrt(delta));
	lambda2[i*W + j] = 0.5*(trace - sqrt(delta));
	if(trace==0)
	  anisotropy[i*W + j] = 0;
	else
	  anisotropy[i*W + j] = 1-2*lambda2[i*W + j]/trace;

	x1 = 2*c[i*W + j];
	x2 = b[i*W + j] - a[i*W + j] - sqrt(delta);
	norm_eig_vect = sqrt(x1*x1 + x2*x2);
	if(norm_eig_vect > 0)
	  {
	    eigvect2_x[i*W + j] = x1/norm_eig_vect;
	    eigvect2_y[i*W + j] = x2/norm_eig_vect;
	  }
	else
	  {
	    eigvect2_x[i*W + j] = 0;
	    eigvect2_y[i*W + j] = 0;
	  }
      }
}
