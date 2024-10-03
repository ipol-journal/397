/*
  Copyright (c) 2020 Luis Alvarez <lalvarez@ulpgc.es> and Yuchen He <yhe306@gatech.edu>

  This program is free software: you can redistribute it and/or modify it under
  the terms of the GNU Affero General Public License as published by the Free
  Software Foundation, either version 3 of the License, or (at your option) any
  later version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
  details.

  You should have received a copy of the GNU Affero General Public License along
  with this program. If not, see <http://www.gnu.org/licenses/>.


  NOTE: The original code by Alvarez and He was slightly modified here:

  - to rename classes with names starting with uppercase letters:
  image2D -> Image2D, point2d -> Point2d.

  - to rename methods with names starting with lowercase letters:
    Neighborhood -> neighborhood, NeighborhoodNormal ->
    neighborhoodNormal, Jacc -> jacc, DSC -> dsc, Bpn -> bpn.
*/

/**
 * \file mainSkeleton2D.cpp
 * \brief Methods to manipulate 2D images and implemente the flux-ordered thinning algorithm
 * \author Yuchen He and Luis Alvarez
 */

#ifndef image2D_H
#define image2D_H
#include <algorithm>
#include <vector>
#include <iostream>
#include <istream>
#include <ostream>
#include <typeinfo>
#include <string>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <fstream>
#include <cstring>
#include <limits>
#include <queue>
#include <set>


#ifdef AMI_OMP_H
  #include <omp.h>
#endif

#include "RGB.h"

#define IMAGE2D_DEBUG

using namespace std;

/**
 * \class  Point2d
 * \brief class  to store 2D points and basic methods
 */
class Point2d{
  public :
  float x /** point x coordinate */;
  float y /** point y coordinate */;
  Point2d():x( (float) 0), y( (float) 0){};
  ~Point2d(){};
  Point2d(const float xx , const float yy){x=xx; y=yy;}
  Point2d(const float scalar){x=y=scalar;}
  Point2d & operator=(const Point2d &p){ x=p.x; y=p.y; return *this;}
  Point2d & operator=(const float scalar){ x=scalar; y=scalar; return *this;}
  Point2d (const Point2d &p){x=p.x; y=p.y;}
  Point2d operator+(const Point2d &p)const {return Point2d(x+p.x,y+p.y);}
  Point2d operator-(const Point2d &p)const {return Point2d(x-p.x,y-p.y);}
  Point2d operator*(const float &a)const {return Point2d(a*x,a*y);}
  Point2d operator/(const float &a)const {return Point2d(x/a,y/a);}
  float operator*(const Point2d &p)const {return ( (float) x*p.x+y*p.y);}
  void print(){ std::cout << "point2d : (" << x << "," << y << ")" << std::endl;}
};


/**
 * \class  Image2D
 * \brief class  to store 3D image2Ds and basic methods
 */
template <class  T>
class Image2D{
  std::vector <T> Image2D_ /** std vector to allocate Image2D */;
  int width_ /** image width */;
  int height_ /** image height */;

  public:

  /// CONSTRUCTORS - DESTRUCTOR METHODS
  Image2D();
  Image2D(const int width,const int height,const T &a);
  Image2D(const int width,const int height);
  Image2D(const Image2D<T> &another_Image2D);
  ~Image2D();
  void clear(){Image2D_.clear(); Image2D_.shrink_to_fit();width_=height_=0;};
  //-----------------------------------------------------------------------------------

  /// class ELEMENT ACCESS METHODS
  inline int width() const {return width_;}
  inline int height() const {return height_;}
  inline int size() const {return Image2D_.size();}
  inline T& operator[](const int &i);
  inline T& operator()(const int &x,const int &y);
  inline double operator()(const double &x,const double &y,int boundary_type=0) const;
  inline const T& operator[](const int &i) const;
  //-----------------------------------------------------------------------------------


  /// BASIC OPERATION METHODS
  /// convert image to RGB. The image is scaled according to  min->0 and max->255
  imageRGB convert2RGB(const T &min,const T &max) const ;
  /// computation of max and min of the image
  T min()const;
  T max()const;
  /// GRADIENT COMPUTATION
  void gradient(Image2D<float> &Ix, Image2D<float> &Iy);

  ///Define 8-neighbors
  vector<int> neighborhood();

  ///compute neighbor normals
  void neighborhoodNormal(vector<float> &Nx,vector<float> &Ny);

  /// Jaccard index
  double jacc(Image2D<T> &I);

  /// DSC similarity coefficient
  double dsc(Image2D<T> &I);

  /// Bpn bias estimator
  double bpn(Image2D<T> &I);

};


/// Jaccard index
template <class  T>
double Image2D<T>::jacc(Image2D<T> &I){

  int size0=0,size1=0;
  for(size_t k=0;k<Image2D_.size();k++){
    if(Image2D_[k]==0.0 && I[k]==0.0)
      size0++;
    if(Image2D_[k]==0.0 || I[k]==0.0)
      size1++;
  }
  return (double) size0/size1;
}


/// DSC similarity coefficient
template <class  T>
double Image2D<T>::dsc(Image2D<T> &I){
  int size0=0,size1=0,size2=0;
  for(size_t k=0;k<Image2D_.size();k++){
    if(Image2D_[k]==0.0){
      size1++;
      if(I[k]==0.0){
        size2++;
        size0++;
      }
      continue;
    }
    if(I[k]==0.0) size2++;
  }
  return (double) 2.*size0/(size1+size2);
}

/// Bpn bias estimator
/// ROY
template <class  T>
double Image2D<T>::bpn(Image2D<T> &I){

  int size0=0,size1=0,size2=0;
  for(size_t k=0;k<Image2D_.size();k++){
    if(Image2D_[k]==0.0){
      if(I[k]==0.0){
        size0++;
      }
      else{
        size2++;
      }
    }else{ /* Image2D != 0.0 */
      if(I[k]==0.0){
        size1++;
      } 
    }
  }
  return (double) (size1-size2)/size0;
}


template <class  T>
Image2D<T>::Image2D(const Image2D<T> &another_Image2D){
  width_ = another_Image2D.width();
  height_ = another_Image2D.height();
  for(int i = 0; i < another_Image2D.size(); i++){
    Image2D_.push_back(another_Image2D[i]);
  }
}

/**
 * \brief Define 8-neighbors
*/
template <class  T>
vector<int> Image2D<T>::neighborhood(){
  vector<int> N(4);
  int k=0;
  N[k++]=1; //0
  N[k++]=-1; //1
  N[k++]=width_; //2
  N[k++]=-width_; //3
  N.push_back(1+width_); //4
  N.push_back(1-width_); //5
  N.push_back(-1+width_); //6
  N.push_back(-1-width_); //7
  return N;
}

//////////////////////////////////////////////////////////////////////////////////////////////
/**
 * \brief compute neighbor normals
*/
template <class  T>
void Image2D<T>::neighborhoodNormal(vector<float> &Nx,vector<float> &Ny){
  Nx=vector<float>(8);
  Ny=vector<float>(8);
  float c=0.707106781186547;
  int k=0;
  Nx[k]=1.; Ny[k++]=0.;
  Nx[k]=-1.; Ny[k++]=0.;
  Nx[k]=0.; Ny[k++]=1.;
  Nx[k]=0.; Ny[k++]=-1.;
  Nx[k]=c; Ny[k++]=c;
  Nx[k]=c; Ny[k++]=-c;
  Nx[k]=-c; Ny[k++]=c;
  Nx[k]=-c; Ny[k++]=-c;
}




//////////////////////////////////////////////////////////////////////////////////////////////
/// computation of min of the image
template <class  T>
T Image2D<T>::min()const {
  T Imin=Image2D_[0];
  for(size_t k=1;k<Image2D_.size();k++){
    if(Imin>Image2D_[k]) Imin=Image2D_[k];
  }
  return Imin;
}

/// computation of max of the image
template <class  T>
T Image2D<T>::max()const {
  T Imax=Image2D_[0];
  for(size_t k=1;k<Image2D_.size();k++){
    if(Imax<Image2D_[k]) Imax=Image2D_[k];
  }
  return Imax;
}
//////////////////////////////////////////////////////////////////////////////////////////////
/**
 * \fn template <class  T> const T& Image2D<T>::operator[](const int &i) const
 * \brief operator [] to acces image value
*/
template <class  T>
const T& Image2D<T>::operator[](const int &i) const
{
  #ifdef IMAGE2D_DEBUG
    if(i>=(int) (Image2D_.size())){
      printf("image<T>: bounds error vector access\n");
      printf("image size()=%ld index to acces the vector =%d\n",Image2D_.size(),i);
      int j; scanf("%d",&j);
#ifndef CMCH
      exit(0);
#else
      throw "Ask me, Carmelo Cuenca";
#endif
    }
  #endif
  return Image2D_[i];
}

//////////////////////////////////////////////////////////////////////////////////////////////
/**
 * \fn template <class  T> T& Image2D<T>::operator[](const int &i)
 * \brief operator [] to acces image value
*/
template <class  T>
T& Image2D<T>::operator[](const int &i)
{
  #ifdef IMAGE2D_DEBUG
    if(i<0 || i>=(int) (Image2D_.size())){
      printf("image<T>: bounds error vector access\n");
      printf("image size()=%ld index to accces the vector =%d\n",Image2D_.size(),i);
      int j; scanf("%d",&j);
      #ifndef CMCH
            exit(0);
      #else
            throw "Ask me, Carmelo Cuenca";
      #endif
    }
  #endif
  return  (Image2D_.at(i));
}

//////////////////////////////////////////////////////////////////////////////////////////////
/**
 * \fn template <class  T> T& image<T>::operator()(const int &x,const int &y,const int &channel)
 * \brief operator () to acces image value
*/
template <class  T>
T& Image2D<T>::operator()(const int &x,const int &y){
  #ifdef IMAGE2D_DEBUG
    unsigned int value=x+y*width_;
    if(value>=(Image2D_.size()) || value<0){
      printf("image<T>: bounds error vector access\n");
      printf("image size()=%d index to accces the vector =%d\n",
                                        Image2D_.size(),value);
      int j; scanf("%d",&j);
      #ifndef CMCH
            exit(0);
      #else
            throw "Ask me, Carmelo Cuenca";
      #endif
    }
    return Image2D_[value];
  #endif
  return Image2D_[x+y*width_];
}


//////////////////////////////////////////////////////////////////////////////////////////////
/**
 * \fn template <class  T> T& image<T>::operator()(const double &x,const double &y,const double &channel)
 * \brief operator () to acces image value using bilinear interpolation
*/
template <class  T>
double Image2D<T>::operator()(const double &x,const double &y,const int type) const{
    int ix=(int) x;
    int iy=(int) y;
    if(type==1 && (ix>=width_ || iy>=height_ || ix<-1|| iy<-1)) return (-3000.);
    else{
      if(ix>=width_) ix=width_-1;
      else if(ix<-1) ix=0;
      if(iy>=height_) iy=height_-1;
      else if(iy<-1) iy=0;
    }
    double dx=x-ix;
    double dy=y-iy;
    int ix1=(ix<(width_-1))?ix+1:ix;
    int iy1=(iy<(height_-1))?iy+1:iy;
    if(ix==-1) ix=0;
    if(iy==-1) iy=0;
    return(
        (1-dx)*(1-dy)*(*this)[ix+iy*width_]+
        (dx)*(1-dy)*(*this)[ix1+iy*width_]+
        (dx)*(dy)*(*this)[ix1+iy1*width_]+
        (1-dx)*(dy)*(*this)[ix+iy1*width_]
    );
}


//////////////////////////////////////////////////////////////////////////////////////////////
/**
 * \brief basic constructor
*/
template <class  T>
Image2D<T>::~Image2D(){width_=0; height_=0;}


//////////////////////////////////////////////////////////////////////////////////////////////
/**
 * \brief basic constructor
*/
template <class  T>
Image2D<T>::Image2D(const int width,const int height,const T &a)
{
  width_=width;
  height_=height;
  int size=width*height;

  Image2D_.resize(size);
  for(int k=0;k<size;k++){
    Image2D_[k]=a;
  }

}


//////////////////////////////////////////////////////////////////////////////////////////////
/**
 * \brief basic constructor
*/
template <class  T>
Image2D<T>::Image2D(const int width,const int height)
{
  width_=width;
  height_=height;
  Image2D_.resize(width*height);
}


//////////////////////////////////////////////////////////////////////////////////////////////
/**
 * \brief basic destructor
*/
template <class  T>
Image2D<T>::Image2D(){}


template <class  T>
imageRGB Image2D<T>::convert2RGB(const T &min,const T &max) const{
  if(min>=max) return imageRGB();

  int height_new=height_;
  int width_new=width_;

  imageRGB IRGB(width_new,height_new);
  double scale=255./(max-min);
  for(int x=0;x<width_new;x++){
    for(int y=0;y<height_new;y++){
      int znew=y*width_new+x;
      float tempf=scale*((float) Image2D_[y*width_+x]-min);
      unsigned char temp=tempf<0.?0:(tempf>255.?255:(unsigned char) tempf);
      IRGB[znew]=color(temp,temp,temp);
    }
  }

  return(IRGB);
}






//-----------------------------------------------------------------------------------

/// FOTA ALGORITHM


//////////////////////////////////////////////////////////////////////////////////////////////
/**
 * \brief Gradient computation
*/
/// GRADIENT COMPUTATION
template <class  T>
void Image2D<T>::gradient (Image2D<float> &Ix, Image2D<float> &Iy)
{
  Ix=Image2D<float>(width_,height_,0.);
  Iy=Image2D<float>(width_,height_,0.);

  double coef1,coef2,c1,d1;
  coef1=sqrt((double) 2.);
  coef2=0.25*(2.-coef1);
  coef1=0.5*(coef1-1);

  for(int y=width_;y<(int) Image2D_.size()-width_ ;y+=width_){
    int k_end=y+width_;
    for(int k=y+1; k<k_end-1; k++){
      c1=Image2D_[k+width_+1]-Image2D_[k-width_-1];
      d1=Image2D_[k-width_+1]-Image2D_[k+width_-1];
      Iy[k]=(float)(coef1*((float) Image2D_[k+width_]-(float)Image2D_[k-width_])+coef2*(c1-d1));
      Ix[k]=-(float)(-(coef1*((float) Image2D_[k+1]- (float) Image2D_[k-1])+coef2*(c1+d1)));
    }
  }

  for(int k=0;k<width_;k++){
    Ix[k]=Ix[k+width_];
    Iy[k]=Iy[k+width_];
    Ix[Image2D_.size()-1-k]=Ix[Image2D_.size()-1-k-width_];
    Iy[Image2D_.size()-1-k]=Iy[Image2D_.size()-1-k-width_];
  }
  for(int k=0;k<height_;k++){
    Ix[k*width_]=Ix[k*width_+1];
    Iy[k*width_]=Iy[k*width_+1];
    Ix[k*width_+width_-1]=Ix[k*width_+width_-2];
    Iy[k*width_+width_-1]=Iy[k*width_+width_-2];
  }
}






//-----------------------------------------------------------------------------------



#endif
