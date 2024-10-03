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

#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "structure_tensors.hpp"
#include "image2D.h"
#include "cocircularity_graphs.hpp"
#include "sparse_io.hpp"
#include "mysparsematrix.h"
#include "max_plus_operators.hpp"

int main(int argc, char *argv[])
{
  if(argc!=12){
    std::cout << "" << std::endl;
    std::cout << "Wrong number of arguments: 11 expected." << std::endl;
    std::cout << "Typical use:" << std::endl;
    std::cout << "" << std::endl;
    std::cout << "anisop operator input_image output_image anisotropy_image sigma p ang_prec cone_ang n_it apply_to mode" << std::endl;
    std::cout << "" << std::endl;
    std::cout << "where:" << std::endl;
    std::cout << "" << std::endl;
    std::cout << "- operator can be one of the following strings: dilation, erosion, opening, closing, asfop (for an ASF starting with an opening), asfclos (for an ASF starting with a closing)" << std::endl;
    std::cout << "- input_image is the path to the input image" << std::endl;
    std::cout << "- output_image is the path to the output image" << std::endl;
    std::cout << "- anisotropy_image is the path to the anisotropy image (output by the program)" << std::endl;
    std::cout << "- sigma is the scale parameter as described in the paper (can be an interger or a float)" << std::endl;
    std::cout << "- p is the parameter which determines the size of the square neighbourhood, as described in the paper (must be an integer)" << std::endl;
    std::cout << "- ang_prec is the angular precision alpha for the binary-weighted co-circularity graph, as described in the paper, except here it is expressed in degrees and should be an integer between 0 and 90." << std::endl;
    std::cout << "- cone_ang is the angle beta ruling the conic condition in the computation of the co-circularity graph, as described in the paper, except here it is expressed in degrees and should be an integer between 0 and 90." << std::endl;
    std::cout << "- n_it is the number of iterations, or more precisely the parameter called k in the paper. It is an integer." << std::endl;
    std::cout << "- apply_to can be one of the two following strings: input or anisotropy. it determines to whch image to apply the operator, the original image or its corresponding anisotropy image, given the parameter sigma." << std::endl;
    std::cout << "- mode can be one of the two following strings: bin or nbin, depending on whether you want to compute the binary-weighted or non binary-weighted graph." << std::endl;
    
    return 0;
  }
  string oper = argv[1];
  string inputname = argv[2];
  string outname = argv[3];
  string anisotropy_outname = argv[4];
  double sigma = atof(argv[5]);
  size_t p = atoi(argv[6]);
  double angular_prec=atof(argv[7]);
  double cone_angle=atof(argv[8]);
  int it=atoi(argv[9]);
  string apply_to = argv[10];
  string mode_str = argv[11];
  
  std::cout << "Running anisotropic operator with the following parameters:" << std::endl;
  std::cout << "- Input image: " << inputname << std::endl;
  std::cout << "- Output image: " << outname << std::endl;
  std::cout << "- Anisotropy image: " << anisotropy_outname << std::endl;
  std::cout << "- Operator: " << oper << std::endl;
  std::cout << "- Sigma: " << sigma << std::endl;
  std::cout << "- Neighbourhood size p: " << p << std::endl;
  std::cout << "- Angular precision on co-circularity: " << angular_prec << std::endl;
  std::cout << "- Cone angle to avoid ladders: " << cone_angle << std::endl;
  std::cout << "- Number of iterations: " << it << std::endl;
  std::cout << "- Apply to: " << apply_to << std::endl;
  std::cout << "- Binary or non-binary weights: " << mode_str << std::endl;

  /// READ THE IMAGE AND ITS NAME
  imageRGB ImIn(argv[2]);
  if(ImIn.size()==0){
    printf("problems reading image %s \n",argv[1]);
  }
  size_t W = ImIn.width();
  size_t H = ImIn.height();

  /// CONVERT RGB IMAGE TO GRAYSCALE
  Image2D<float> im(W,H);
  for(int k=0;k<ImIn.size();k++){
    im[k]= (ImIn[k].R_+ImIn[k].G_+ImIn[k].B_)/3.0;
  }

  /// BUILD THE CO-CIRCULARITY GRAPH, REPRESENTED BY A SPARSE MATRIX
  std::cout << "Building the co-circularity graph..." << std::endl;
  double thresh_cocirc, thresh_cone;
  thresh_cocirc = cos(M_PI*angular_prec/180); //default: cos(M_PI/20);
  thresh_cone = cos(M_PI*cone_angle/180);//default: cos(M_PI/6);
  size_t n_col = W*H;
  MySparseMatrix<double> S(n_col);
  Image2D<float> im_anisotropy(W,H);
  image_to_cocircularity_graph(im, S, im_anisotropy, sigma, p, thresh_cocirc, thresh_cone, mode_str);
  std::cout << "DONE." << std::endl;

  size_t n_non_zeros;
  n_non_zeros = count_sparse_matrix_non_zeros(S, n_col);
  std::cout << "The sparse matrix contains "<< n_non_zeros << " non zero entries ("<< (float) 100*n_non_zeros/(W*H*W*H) << "%)." << std::endl;

  /// REPRESENTS THE CONTENT OF THE MATRIX BY THE INDICES AND VALUES
  /// OF NON-ZERO ENTRIES
  size_t* I, *J;
  I = new size_t[n_non_zeros];
  J = new size_t[n_non_zeros];
  double* Weights;
  Weights = new double[n_non_zeros];
  sparse_matrix_to_arrays(S, I, J, Weights, n_col);


  /// INSTANCIATE AND FILL THE IMAGE TO WHICH THE OPERATOR WILL BE
  /// APPLIED (INPUT OR ANISOTROPY)
  Image2D<float> im_in(W,H);
  if(apply_to == "anisotropy"){
    for(int k=0;k<ImIn.size();k++){
      im_in[k]= im_anisotropy[k];
    }
  }
  else{
    for(int k=0;k<ImIn.size();k++)
      im_in[k]= im[k];
  }

  /// APPLY OPERATOR
  Image2D<float> im_res(W,H);
  if(oper == "erosion")
    erode_max_plus_symmetric_iterated(I, J, Weights, it, n_non_zeros, im_in, im_res);
  if(oper == "dilation")
    dilate_max_plus_symmetric_iterated(I, J, Weights, it, n_non_zeros, im_in, im_res);
  if(oper == "opening")
    open_max_plus_symmetric(I, J, Weights, it, n_non_zeros, im_in, im_res);
  if(oper == "closing")
    close_max_plus_symmetric(I, J, Weights, it, n_non_zeros, im_in, im_res);
  if(oper == "asfop"){
    // ASF starting with an opening
    for(size_t n_it=1;n_it<it+1;n_it++){
      open_max_plus_symmetric(I, J, Weights, n_it, n_non_zeros, im_in, im_res);
      for(int k=0;k<ImIn.size();k++)
	im_in[k]= im_res[k];
      close_max_plus_symmetric(I, J, Weights, n_it, n_non_zeros, im_in, im_res);
      for(int k=0;k<ImIn.size();k++)
	im_in[k]= im_res[k];
    }
  }
  if(oper == "asfclos"){
    // ASF starting with a closing
    for(size_t n_it=1;n_it<it+1;n_it++){
      close_max_plus_symmetric(I, J, Weights, n_it, n_non_zeros, im_in, im_res);
      for(int k=0;k<ImIn.size();k++)
	im_in[k]= im_res[k];
      open_max_plus_symmetric(I, J, Weights, n_it, n_non_zeros, im_in, im_res);
      for(int k=0;k<ImIn.size();k++)
	im_in[k]= im_res[k];
    }
  }

  
  // INSTANCIATE AND WRITE OUTPUT IMAGE
  imageRGB im_out(W,H);
  for (size_t i=0; i<H; i++)
    for (size_t j=0; j<W; j++)
      im_out[i*W + j] = im_res[i*W + j];
  int nlen = outname.length();
  char name[nlen+1];
  strncpy(name,outname.c_str(), nlen);
  name[nlen]='\0';
  im_out.write(name);

  // INSTANCIATE AND WRITE ANISOTROPY IMAGE
  imageRGB im_anis(W,H);
  for (size_t i=0; i<H; i++)
    for (size_t j=0; j<W; j++)
      im_anis[i*W + j] = im_anisotropy[i*W + j];
  int nlen2 = anisotropy_outname.length();
  char name2[nlen2+1];
  strncpy(name2,anisotropy_outname.c_str(), nlen2);
  name2[nlen2]='\0';
  im_anis.write(name2);

  delete[] I;
  delete[] J;
  delete[] Weights;

  return 0;
}
