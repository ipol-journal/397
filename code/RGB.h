/*
 * Copyright (c) 2013-2014, Luis Alvarez <lalvarez@dis.ulpgc.es>
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under, at your option, the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version, or
 * the terms of the simplified BSD license.
 *
 * You should have received a copy of these licenses along this
 * program. If not, see <http://www.gnu.org/licenses/> and
 * <http://www.opensource.org/licenses/bsd-license.html>.
 */

#ifndef _RGB_H
#define _RGB_H
#ifndef COLOR_H_
#define COLOR_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <math.h>

#define color_type float

#include "lodepng.h"

using namespace std;

typedef struct {
   unsigned short int type;                 /* Magic identifier            */
   unsigned int size;                       /* File size in bytes          */
   unsigned short int reserved1, reserved2;
   unsigned int offset;                     /* Offset to image data, bytes */
} HEADER;
typedef struct {
   unsigned int size;               /* Header size in bytes      */
   int width,height;                /* Width and height of image */
   unsigned short int planes;       /* Number of colour planes   */
   unsigned short int bits;         /* Bits per pixel            */
   unsigned int compression;        /* Compression type          */
   unsigned int imagesize;          /* Image size in bytes       */
   int xresolution,yresolution;     /* Pixels per meter          */
   unsigned int ncolours;           /* Number of colours         */
   unsigned int importantcolours;   /* Important colours         */
} INFOHEADER;


class color{
  public :
  color_type R_ /** red color channel */;
  color_type G_ /** green color channel */;
  color_type B_ /** blue color channel */;
  color(unsigned char  R,unsigned char G,unsigned char B){R_=(color_type) R; G_=(color_type)G; B_=(color_type)B;}
  color(color_type c_){R_=c_; G_=c_; B_=c_;}
  color(){};
};


class imageRGB{

  int width_;
  int height_;
  int size_;
  vector<color> RGB_;

  public:
  inline int width() const {return width_;}
  inline int height() const {return height_;}
  inline int size() const {return width_*height_;}
  inline color &operator[](int k){
    return(RGB_[k]);
  }
  inline color operator[](int k)const {return(RGB_[k]);}
  void equal(color c){for(int k_=0;k_<(int) RGB_.size();k_++) RGB_[k_]=c;}


  imageRGB(const char *fileName /** INPUT FILE NAME */ );
  int write(const char *name);
  imageRGB(){width_=height_=size_=0;}
  imageRGB(int width,int height){
    width_=width; height_=height; size_=width*height; RGB_.resize(size_);
  };
  imageRGB(int width,int height,color c){
    width_=width; height_=height; size_=width*height; RGB_.resize(size_);
    for(int k=0;k<size_;k++){RGB_[k]=c;}
  };
};

#endif /* !_COLOR_H */
#endif /* !_RGB_H */







