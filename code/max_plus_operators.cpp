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

#include "max_plus_operators.hpp"

/// Computes the dilation as max-plus matrix-to-vector product
///            im_out[j] = max_i {im_in[i] + M[i,j]}
/// where the matrix M is supposed to be symmetric and is
/// represented by the three arrays I, J, W:
///          M[I[l], J[l]] = W[l] and zero for other indices.
int dilate_max_plus_symmetric(size_t* I, size_t* J, double* W,
			      size_t n_lines,
			      const Image2D<float> &imIn,
			      Image2D<float> &imOut)
{
  size_t l, i, j;
  double w;
  int change=0;
  // copy(imIn, imOut);
  for(int k=0;k<imIn.size();k++){
    imOut[k]= imIn[k];
  }
  for(l=0;l<n_lines;l++)
    {
      i = I[l];
      j = J[l];
      w = 255*log(W[l]);
      if(imOut[i] < imIn[j]+w)
	{
	  imOut[i] = imIn[j]+w;
	  change=1;
	}
      if(imOut[j] < imIn[i]+w)
	{
	  imOut[j] = imIn[i]+w;
	  change=1;
	}
    }
  return change;
}


/// Computes iterations of the symmetric max-plus dilation.
int dilate_max_plus_symmetric_iterated(size_t* I, size_t* J, double* W,
				       int k,
				       size_t n_lines,
				       const Image2D<float> &imIn,
				       Image2D<float> &imOut)
{
  Image2D<float> imAux(imIn);
  int it, change;
  // copy(imIn, imAux);
  for(int i=0;i<imIn.size();i++){
    imAux[i]= imIn[i];
  }
  it=0;
  change=1;
  while(it<k && change)
    {
      it++;
      std::cout << "Iterated dilations: dilation "<<it<< std::endl;
      change=dilate_max_plus_symmetric(I, J, W, n_lines, imAux, imOut);
      if(change)
	{
	  // copy(imOut, imAux);
	  for(int i=0;i<imIn.size();i++){
	    imAux[i]= imOut[i];
	  }
	}
      else
	std::cout << "Iterated dilations: no change after iteration "<<it-1<< std::endl;
    }
  if(change)
    return k;
  else
    return it-1;
}


/// Computes the erosion as min-plus matrix-to-vector product
///            im_out[j] = min_i {im_in[i] - M[j,i]}
/// where the matrix M is supposed to be symmetric and is
/// represented by the three arrays I, J, W:
///          M[I[l], J[l]] = W[l] and zero for other indices.
int erode_max_plus_symmetric(size_t* I, size_t* J, double* W,
			     size_t n_lines,
			     const Image2D<float> &imIn,
			     Image2D<float> &imOut)
{
  size_t l, i, j;
  double w;
  int change=0;
  // copy(imIn, imOut);
  for(int k=0;k<imIn.size();k++)
    imOut[k]= imIn[k];
  
  for(l=0;l<n_lines;l++)
    {
      i = I[l];
      j = J[l];
      w = 255*log(W[l]);
      if(imOut[j]+w> imIn[i])
	{
	  imOut[j] = imIn[i]-w;
	  change=1;
	}
      if(imOut[i]+w> imIn[j])
	{
	  imOut[i] = imIn[j]-w;
	  change=1;
	}
    }
  return change;
}


/// Computes iterations of the symmetric min-plus erosion.
int erode_max_plus_symmetric_iterated(size_t* I, size_t* J, double* W,
				      int k,
				      size_t n_lines,
				      const Image2D<float> &imIn,
				      Image2D<float> &imOut)
{ 
  Image2D<float> imAux(imIn);
  int it, change;
  // copy(imIn, imAux);
  for(int i=0;i<imIn.size();i++)
    imAux[i]= imIn[i];
  it=0;
  change=1;
  while(it<k && change)
    {
      it++;
      std::cout << "Iterated erosions: erosion "<<it<< std::endl;
      change = erode_max_plus_symmetric(I, J, W, n_lines, imAux, imOut);
      if(change)
	{
	  // copy(imOut, imAux);
	  for(int i=0;i<imIn.size();i++)
	    imAux[i]= imOut[i];
	}
      else
	std::cout << "Iterated erosions: no change after iteration "<<it-1<< std::endl;
    }
  if(change)
    return k;
  else
    return it-1;
}


/// Computes the opening based on adjoint iterated symmetric max-plus
/// dilations and min-plus erosions
int open_max_plus_symmetric(size_t* I, size_t* J, double* Weights,
			    int it,
			    size_t n_non_zeros,
			    const Image2D<float> &imIn,
			    Image2D<float> &imOut)
{
  size_t W = imIn.width();
  size_t H = imIn.height();
  Image2D<float> im_erod(W,H);
  int kmax;
  kmax = erode_max_plus_symmetric_iterated(I, J, Weights, it, n_non_zeros, imIn, im_erod);
  if(kmax==it)
    dilate_max_plus_symmetric_iterated(I, J, Weights, kmax, n_non_zeros, im_erod, imOut);
  else
    {
      std::cout << "No need to dilate."<< std::endl;
      for(int i=0; i<imIn.size();i++)
	imOut[i]=im_erod[i];
    }
  return 1;
}

/// Computes the closing based on adjoint iterated symmetric max-plus
/// dilations and min-plus erosions
int close_max_plus_symmetric(size_t* I, size_t* J, double* Weights,
			     int it,
			     size_t n_non_zeros,
			     const Image2D<float> &imIn,
			     Image2D<float> &imOut)
{
  size_t W = imIn.width();
  size_t H = imIn.height();
  Image2D<float> im_dil(W,H);
  int kmax;
  kmax = dilate_max_plus_symmetric_iterated(I, J, Weights, it, n_non_zeros, imIn, im_dil);
  if(kmax==it)
    erode_max_plus_symmetric_iterated(I, J, Weights, kmax, n_non_zeros, im_dil, imOut);
  else
    {
      std::cout << "No need to erode."<< std::endl;
      for(int i=0; i<imIn.size();i++)
  	imOut[i]=im_dil[i];
    }
  return 1;
}
