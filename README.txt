This directory contains code related to the paper

"Adaptive anisotropic morphological filtering based on co-circularity
of local orientations"

by

Samy Blusseau, Santiago Velasco-Forero, Jesús Angulo, Isabelle Bloch.

published in the IPOL special Issue on Reproducible Research in
Mathematical Morphology and Discrete Geometry (TC18 Special Issue).

It contains the corresponding C/C++ code. The main function is in the
  anisop.cpp file.


**********************
Compiling the code
**********************

To compile the code, the following sequence of five commands should be sufficient
(provided cmake is installed):

cd code
mkdir build
cd build
cmake ..
make


**********************
Running the code
**********************

Once compiled, you can run the program as follows:

./build/anisop operator input_image output_image anisotropy_image
sigma p ang_prec cone_ang n_it apply_to mode

where

- operator can be one of the following strings: dilation, erosion,
  opening, closing, asfop (for an ASF starting with an opening),
  asfclos (for an ASF starting with a closing)

- input_image is the path to the input image

- output_image is the path to the output image

- anisotropy_image is the path to the anisotropy image (output by the
  program)

- sigma is the scale parameter as described in the paper (can be an
  interger or a float)

- p is the parameter which determines the size of the square
  neighbourhood, as described in the paper (must be an integer)

- ang_prec is the angular precision alpha for the binary-weighted
  co-circularity graph, as described in the paper, except here it is
  expressed in degrees and should be an integer between 0 and 90.

- cone_ang is the angle beta ruling the conic condition in the
  computation of the co-circularity graph, as described in the paper,
  except here it is expressed in degrees and should be an integer
  between 0 and 90.

- n_it is the number of iterations, or more precisely the parameter
  called "k" in the paper. It is an integer.

- apply_to can be one of the two following strings: input or
  anisotropy. it determines to whch image to apply the operator, the
  original image or its corresponding anisotropy image, given the
  parameter sigma.

- mode can be one of the two following strings: bin or nbin, depending
  on whether you want to compute the binary-weighted or non
  binary-weighted graph.



**********************
Examples
**********************

You can test the code on one of the six test images with appropriate
parameters, using the test.sh script. To do so you should run (in the
code directory)

./test.sh IMNAME

where IMNAME must be one of the six image names: y, k, retine,
eye_crop, autumn or klok.
