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

#include "cocircularity_graphs.hpp"

/// Stores a list of pixel coordinates in row_list and col_list. The
/// corresponding set of pixels is composed of two rectangles,
/// respectively bounded by [i_min1, i_max1]x[j_min1, j_max1] and
/// [i_min2, i_max2]x[j_min2, j_max2]
void pixel_list(size_t* row_list, size_t* col_list,
		size_t i_min1, size_t i_min2,
		size_t i_max1, size_t i_max2,
		size_t j_min1, size_t j_min2,
		size_t j_max1, size_t j_max2)
{
  size_t k, l, index;
  index=0;
  // First rectangle
  for(k=i_min1; k<=i_max1; ++k)
    for(l=j_min1; l<=j_max1;++l)
      {
	row_list[index]=k;
	col_list[index]=l;
	index++;
      }
  // Second rectangle
  for(k=i_min2; k<=i_max2; ++k)
    for(l=j_min2; l<=j_max2;++l)
      {
	row_list[index]=k;
	col_list[index]=l;
	index++;
      }
}

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
			      int bin)
{ 
  clock_t start_time = clock();
  size_t i, j, i_min1, i_max1, j_min1, j_max1, i_min2, i_max2, j_min2, j_max2, k, l;
  double x1, x2, y1, y2, v1x, v1y, v2x, v2y;
  double e_1 = 1/exp(1);// below this value, the weight will not
			// contribute to operators so it is not
			// included in the sparse matrix S. Indeed, we
			// will take w_{ij} = -255*log(S_{ij}) as
			// weight for 8-bits images, and w_{ij} <=
			// -255 is equivalent to w_{ij} = -\infty.

  // Browse all pixels
  for (i=0; i<H; i++)
    {
      // Define the first dimension of the squared neighbourhood
      // Takes into account the image boundaries
      i_min1 = i+1;
      i_max1 = ((i+p < H) ? i+p : H-1);
      i_min2 = ((i > p) ? i-p : 0);
      i_max2 = i;
      for (j=0; j<W; ++j)
	{
	  // Define the second dimension of the squared neighbourhood
	  // Takes into account the image boundaries
	  j_min1 = j;
	  j_max1 = ((j+p < W) ? j+p : W-1);
	  j_min2 = j+1;
	  j_max2 = j_max1;
	  // Coordinates of the central pixel and principal direction
	  // of its smoothed structure tensor
	  x1 = j;
	  y1 = i;
	  v1x = eigvect2_x[i*W+j];
	  v1y = eigvect2_y[i*W+j];

	  // Define the list of neighbour pixels to test
	  size_t h1 = ((i_max1-i_min1+1 > 0) ? i_max1-i_min1+1 : 0);
	  size_t w1 = ((j_max1-j_min1+1 > 0) ? j_max1-j_min1+1 : 0);
	  size_t h2 = ((i_max2-i_min2+1 > 0) ? i_max2-i_min2+1 : 0);
	  size_t w2 = ((j_max2-j_min2+1 > 0) ? j_max2-j_min2+1 : 0);
	  size_t n_pixels1 = h1*w1;
	  size_t n_pixels2 = h2*w2;
	  size_t n_pixels = n_pixels1+n_pixels2;
	  size_t* row_list = new size_t[n_pixels];
	  size_t* col_list = new size_t[n_pixels];
	  pixel_list(row_list, col_list,
	  	     i_min1, i_min2, i_max1, i_max2,
	  	     j_min1, j_min2, j_max1, j_max2);
	  size_t index;

	  for(index=0; index<n_pixels; index++)
	    {
	      l = col_list[index];
	      k = row_list[index];
	      // Coordinates of a neighbour pixel and principal
	      // direction of its smoothed structure tensor
	      x2 = l;
	      y2 = k;
	      v2x = eigvect2_x[k*W+l];
	      v2y = eigvect2_y[k*W+l];
	      // Determine if the central pixel and the neighbour
	      // pixel should be adjacent in the co-circularity graph
	      // and with which weight (if non-binary case)
	      double resp = test_adjacency(x1, y1, v1x, v1y,
	      				   x2, y2, v2x, v2y,
	      				   thresh_cocirc, thresh_cone, bin);
	      
	      if(resp>e_1)
		S(i*W+j,k*W+l) = resp;
	    }
	  delete[] row_list;
	  delete[] col_list;
	}
    }
  clock_t end_time = clock();
  double elapsed_secs = double(end_time - start_time) / CLOCKS_PER_SEC;
  std::cout << "Time elapsed: "<< elapsed_secs << " s." << std::endl;
}

/// Determine whether point 1 (x1, y1) with direction (v1x, v1y) and
/// point 2 (x2, y2) with direction (v2x, v2y) are co-circular up to
/// the angular precision thresh_cocirc (binary case), or to which
/// degree they are co-circular (non-binary case) and if the conic
/// constraint is satisfied. The answer resp is in [0, 1] and we will
/// take 255*log(resp) as weight between the two corresponding
/// vertices when computing morphological operators (see
/// dilate_max_plus_symmetric() and erode_max_plus_symmetric() in
/// max_plus_operators.cpp).
double test_adjacency(double x1, double y1, double v1x, double v1y,
		      double x2, double y2, double v2x, double v2y, 
		      double thresh_cocirc, double thresh_cone,
		      int bin)
{
  double dx, dy, n, v1_dot_d;
  double resp = 0.;
  int conic_constraint;
  double cocircularity=0.;

  // First test conic constraint
  dx = x2-x1;
  dy = y2-y1;
  n = sqrt(dx*dx + dy*dy);
  if(n>0)
    {
      dx = dx/n; 
      dy = dy/n; 
    }
  v1_dot_d = v1x*dx + v1y*dy;
  conic_constraint = (abs(v1_dot_d)>=thresh_cone);

  // Then test co-circularity if conic constraint is fullfilled
  if(conic_constraint)
    {
      double v0x, v0y;
      v0x = 2*v1_dot_d*dx - v1x;
      v0y = 2*v1_dot_d*dy - v1y;
      cocircularity = abs(v0x*v2x + v0y*v2y);
      if(bin)
      	resp=(cocircularity>=thresh_cocirc);
      else
      	resp = cocircularity;
    }
  return resp;
}

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
				  string mode_str)
{  
  size_t W = im.width();
  size_t H = im.height();
  size_t len = W*H;

  double* gx, *gy;
  gx = new double[len];
  gy = new double[len];
  lin_grad_gray(im, gx, gy);

  // Coeffs of the rank 1 structure tensor T0
  //           
  //            | a  c |
  //       T0 = | c  b |
  double* a, *b, *c;
  a = new double[len];
  b = new double[len];
  c = new double[len];
  // Compute the rank 1 structure tensor T0
  structure_tensors0(gx, gy, H, W, a, b, c);

  // Coeffs of the rank > 1 smoothed structure tensor T_sigma
  //           
  //                 | a_reg  c_reg |
  //       T_sigma = | c_reg  b_reg |
  double* a_reg, *b_reg, *c_reg;
  a_reg = new double[len];
  b_reg = new double[len];
  c_reg = new double[len];
  smooth_structure_tensors(a, b, c, sigma, H, W, a_reg, b_reg, c_reg);


  // Diagonalize structre tensor of scale sigma
  // Compute main eigenvector
  double* lambda1, *lambda2, *eigvect2_x, *eigvect2_y, *anisotropy;
  lambda1 = new double[len]; // Largest eigenvalue
  lambda2 = new double[len];// smallest eigenvalue
  eigvect2_x = new double[len];// x component of 2nd eigenvector
  eigvect2_y = new double[len];// y component of 2nd eigenvector
  anisotropy = new double[len];// anisotropy defined as 1 - 2*lambda2/(lambda1+lambda2)
  diagonalize_structure_tensors(a_reg, b_reg, c_reg, W, H, 
				lambda1, lambda2,
				eigvect2_x, eigvect2_y,
				anisotropy);
  
  // Build image of the tensors' anisotropy
  size_t i, j;
  for (i=0; i<H; i++)
    for (j=0; j<W; j++)
      im_anisotropy[i*W + j] = 255*anisotropy[i*W + j];

  // Fill in the sparse matrix S representing the weighted
  // co-circularity graph
  string bin_str ="bin";
  int bin;
  if(mode_str == bin_str)
    bin=1;
  else
    bin=0;
  compute_adjacency_matrix(S, p, W, H,  eigvect2_x, eigvect2_y, thresh_cocirc, thresh_cone, bin);
  
  delete[] gx;
  delete[] gy;
  delete[] a;
  delete[] a_reg;
  delete[] b;
  delete[] b_reg;
  delete[] c;
  delete[] c_reg;
  delete[] lambda1;
  delete[] lambda2;
  delete[] eigvect2_x;
  delete[] eigvect2_y;
  delete[] anisotropy;
}
