#include "RGB.h"
#include <stdint.h>
#include <iostream>     // std::cout
#include <fstream>      // std::ifstream

imageRGB::imageRGB(const char *name /** INPUT FILE NAME */ ){

    width_=height_=size_=0;

   /// WE CHECK THE EXTENSION OF THE FILE
    int N=strlen(name);
    if(N<4) return;
    if( (name[N-3]=='p' || name[N-3]=='P') && (name[N-2]=='n' || name[N-2]=='N') && (name[N-1]=='g' || name[N-1]=='G') ){
      /// PNG FORMAT FILE
      vector<unsigned char> image;
      unsigned int width0,height0;
      unsigned int error = lodepng::decode(image, width0, height0,name);
      if(error){
        printf("Problems reading image %s\n",name);
        return;
      }
      width_=width0;
      height_=height0;
      size_=width0*height0;
      RGB_.resize(size_);

      if(image.size()!=(unsigned int) 4*size_){
        printf("Problems with dimension of image %s\n",name);
        printf("size_=%d image.size()/4=%d\n",size_,(int) image.size()/4);
        return;
      }

      unsigned int cont=0;
      for(int y=height_-1;y>=0;y--){
        for(int x=0;x<width_;x++){
          unsigned int k=y*width_+x;
          RGB_[k].B_=image[cont++];
          RGB_[k].G_=image[cont++];
          RGB_[k].R_=image[cont++];
          cont++;
        }
      }
      return;
   }

   /// RGB FORMAT FILE

   HEADER header;
   INFOHEADER infoheader;
   FILE *fptr;
   long size,m,cont;
   unsigned char *image;


   /* Open file */
   if ((fptr = fopen(name,"rb")) == NULL) {
      fprintf(stderr,"Unable to open RGB file \"%s\"\n",name);
      return;
   }

   /* Read and check the header */
   if (fread(&header.type,sizeof(unsigned short),1,fptr) != 1) {
       fprintf(stderr,"Image read failed\n");
       return;
   }

   /* Read and check the header */
   if (header.type != (int) 'M'*256+'B') {
     printf("Be careful, the file is not a RGB image file (magic number not correct)\n");
   }

   if (fread(&header.size,sizeof(unsigned int),1,fptr) != 1) {
               fprintf(stderr,"Image read failed\n");
               return;
   }

   if (fread(&header.reserved1,sizeof(unsigned short),1,fptr) != 1) {
               fprintf(stderr,"Image read failed\n");
               return;
   }

   if (fread(&header.reserved2,sizeof(unsigned short),1,fptr) != 1) {
               fprintf(stderr,"Image read failed\n");
               return;
   }
   if (fread(&header.offset,sizeof(unsigned int),1,fptr) != 1) {
               fprintf(stderr,"Image read failed\n");
               return;
    }

   /* Read and check the information header */
   if (fread(&infoheader,sizeof(INFOHEADER),1,fptr) != 1) {
      fprintf(stderr,"Failed to read RGB info header\n");
      return;
   }

   /* WE CHECK THE 24 bits */
   if(infoheader.bits!=24){
	   fprintf(stderr,"Be carefull, function only implemented for 24 bits colour format\n");
   }
   /* WE CHECK COMPRESSION */
   if(infoheader.compression!=0){
	   fprintf(stderr,"Sorry, function only implemented image without compression \n");
       return;
   }
   /* WE CHECK THE 24 bits */
   if(infoheader.ncolours!=0){
	   fprintf(stderr,"Be carefull, function only implemented for 24 bits colour format\n");
   }
   /* WE CHECK THE 24 bits */
   if(infoheader.importantcolours!=0){
	   fprintf(stderr,"Be carefull, function only implemented for 24 bits colour format\n");
   }


   /* WE FIT DIMENSIONS AND ALLOCATE MEMORY FOR THE IMAGE */
   width_=infoheader.width; height_=infoheader.height;
   size=width_*height_;
   RGB_.resize(size);

  /* WE STORE THE IMAGE INFORMATION IN image */
   image=(unsigned char*)malloc(3*size*sizeof(unsigned char));

   /* We read the image */
   if (fread(image,size*3,1,fptr) != 1) {
     fprintf(stderr,"Image read failed\n");
     free(image);
     return;
   }

   cont=0;
   for(m=0;m<size;m++){
       RGB_[m].R_=image[cont];
       RGB_[m].G_=image[cont+1];
       RGB_[m].B_=image[cont+2]; cont+=3;
   }
   fclose(fptr);

}


int imageRGB::write(const char *name){

    /// WE CHECK THE EXTENSION OF THE FILE
    int N=strlen(name);
    if(N<4) return -1;
    if( (name[N-3]=='p' || name[N-3]=='P') && (name[N-2]=='n' || name[N-2]=='N') && (name[N-1]=='g' || name[N-1]=='G') ){
      /// PNG FORMAT FILE
      vector<unsigned char> image((unsigned int) width_ * height_ * 4);
      unsigned int cont=0;
      for(int y=height_-1;y>=0;y--){
        for(int x=0;x<width_;x++){
          unsigned int k=y*width_+x;
          image[cont++]=(unsigned int) RGB_[k].B_<0?0:RGB_[k].B_>255?255:(round(RGB_[k].B_));
          image[cont++]=(unsigned int) RGB_[k].G_<0?0:RGB_[k].G_>255?255:(round(RGB_[k].G_));
          image[cont++]=(unsigned int) RGB_[k].R_<0?0:RGB_[k].R_>255?255:(round(RGB_[k].R_));
          image[cont++]=255;
        }
      }
      unsigned error = lodepng::encode(name, image, (unsigned int) width_, (unsigned int) height_);

      if(error){
        printf("Problems writing image %s\n",name);
        return-1;
      }

      return 0;

    }



   /// RGB FORMAT FILE
   int width=4*(width_/4);
   int height=4*(height_/4);
   if(width==0 || height==0) return -1;
   if(width!=width_ || height!=height_){
     imageRGB I2(width,height);
     for(int x=0;x<width;x++){
       for(int y=0;y<height;y++){
         I2[y*width+x]=RGB_[y*width_+x];
       }
     }
     return I2.write(name);
   }

   HEADER header;
   INFOHEADER infoheader;
   FILE *fptr;
   long cont,size,m;
   unsigned char *image;

   /* Open file */
   if ((fptr = fopen(name,"wb")) == NULL) {
      fprintf(stderr,"Unable to open RGB file \"%s\"\n",name);
      return(-1);
   }

   /* Header Definition */
   header.type = (unsigned short int) 'M'*256+'B';
   header.size = (unsigned int) 3*width_*height_+54;
   header.reserved1 = 0;
   header.reserved2 = 0;
   header.offset = (unsigned int) 54;


  /* We write the header */
   if (fwrite(&header.type,sizeof(unsigned short),1,fptr) != 1) {
      fprintf(stderr,"Failed to write RGB info header\n");
      return(-1);
   }
   if (fwrite(&header.size,sizeof(unsigned int),1,fptr) != 1) {
      fprintf(stderr,"Failed to write RGB info header\n");
      return(-1);
   }
   if (fwrite(&header.reserved1,sizeof(unsigned short),1,fptr) != 1) {
      fprintf(stderr,"Failed to write RGB info header\n");
      return(-1);
   }
   if (fwrite(&header.reserved2,sizeof(unsigned short),1,fptr) != 1) {
      fprintf(stderr,"Failed to write RGB info header\n");
      return(-1);
   }
   if (fwrite(&header.offset,sizeof(unsigned int),1,fptr) != 1) {
      fprintf(stderr,"Failed to write RGB info header\n");
      return(-1);
   }


   /* INFOHEADER DEFINITION */
   infoheader.size = 40;
   infoheader.width = width_;
   infoheader.height = height_;
   infoheader.planes = (unsigned short int) 1;
   infoheader.bits = (unsigned short int) 24;
   infoheader.compression = (unsigned int) 0;
   infoheader.imagesize=3*width_*height_;
   infoheader.xresolution=2834;
   infoheader.yresolution=2834;
   infoheader.ncolours = (unsigned int) 0 ;
   infoheader.importantcolours = (unsigned int) 0;



   /* We write the infoheaderinformation */
   if (fwrite(&infoheader,sizeof(INFOHEADER),1,fptr) != 1) {
      fprintf(stderr,"Failed to write RGB info header\n");
      return(-1);
   }

   size=width_*height_;

   /* WE STORE THE IMAGE INFORMATION IN image */
   image=(unsigned char*)malloc((size+3)*sizeof(unsigned char));

   long lim[4]={0,size/3,2*size/3,size};

   for(int k=0;k<3;k++){
     cont=0;
     for(m=lim[k];m<lim[k+1];m++){
       image[cont]=RGB_[m].R_;
       image[cont+1]=RGB_[m].G_;
       image[cont+2]=RGB_[m].B_;
       cont+=3;
     }

     /* We write the image */
       if (fwrite(image,(int) (3*(lim[k+1]-lim[k])),1,fptr) != 1) {
       fprintf(stderr,"Image write failed\n");
       free(image);
       return(-1);
     }
   }

   fclose(fptr);
   free(image);
   return(0);

}








