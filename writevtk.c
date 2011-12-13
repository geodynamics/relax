/*************************************************************
*  export vectors and tensors in big-endian mixed ascii/binary
*  vtk format for Paraview.
*
*  sylvain barbot 10/27/11 - original form
*************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "config.h"

// check data alignment
#ifdef FFTW3
#define ALIGN_DATA 1
#else
#ifdef SGI_FFT
#define ALIGN_DATA 1
#else
#ifdef IMKL_FFT
#define ALIGN_DATA 1
#endif
#endif
#endif

int fix(int n){
  if(n<0) return n-1;return n+1;
}

float swap(float d){ 
	float a; 
	unsigned char *dst = (unsigned char *)&a; 
	unsigned char *src = (unsigned char *)&d; 
	dst[0] = src[3];
	dst[1] = src[2];
	dst[2] = src[1];
	dst[3] = src[0]; 
	return a; 
} 

// test endianness of the machine
unsigned char isbigendian(){

  typedef union{
        int i;
        char c[4];
  } u;
  u temp;

  temp.i = 0x12345678;
  
  switch(temp.c[0]) {
     case 0x12:
        return 1u; // big endian
     case 0x78:
        return 0u; // little endian
     default:
        fprintf(stderr,"invalid result for endianness test.\n");
        fprintf(stderr,"temp %x %x %x %x.\n",temp.c[0],temp.c[1],temp.c[2],temp.c[3]);
        fprintf(stderr,"temp %x \n",temp.c[0]);
        return 2u;
  }
}

/*************************************************************
*  subroutine ExportVTK_Vectors_Legacy
*  creates a .vtk file in the VTK Legacy binary format with 
*  structured points containing a vector field.
*
*  sylvain barbot 10/27/11 - original form
*************************************************************/
void exportvtk_vectors_legacy_(u1,u2,u3,sx1,sx2,sx3,dx1,dx2,dx3,j1,j2,j3,filename,title,name)
  float *u1, *u2, *u3;    /* data array for output */
  int *sx1, *sx2, *sx3;   /* number of points      */
  double*dx1, *dx2, *dx3; /* sampling distance     */ 
  int *j1, *j2, *j3;      /* subsampling rate      */
  char *filename;         /* output file name      */
  char *title;            /* output file name      */
  char *name;             /* output file name      */
{

  FILE * funit;
  float buffer[3];
  int i1,i2,k1,k2,k3,index;
  unsigned char endian;

  funit=fopen(filename,"wb");
  if (NULL==funit){
     fprintf(stderr,"could not open file %s for vtk output\n",filename);
     fprintf(stderr,"exiting.\n");
     return;
  }

  // find endianness
  endian=isbigendian();

  // writing header of file
  fprintf(funit,"# vtk DataFile Version 3.0\n");
  fprintf(funit,"%s\n",title);
  fprintf(funit,"BINARY\n");
  fprintf(funit,"DATASET STRUCTURED_POINTS\n");

  // structured points grid
  fprintf(funit,"DIMENSIONS %i %i %i\n",(*sx1)/(*j1),(*sx2)/(*j2),(*sx3)/(*j3)); 
  fprintf(funit,"ORIGIN %f %f %f\n",-(*dx1)*((*sx1)/2),-(*dx2)*((*sx2)/2),0.0); 
  fprintf(funit,"SPACING %f %f %f\n",(*dx1)*(*j1),(*dx2)*(*j2),(*dx3)*(*j3)); 

  // data header for this grid
  fprintf(funit,"POINT_DATA %i\n",((*sx1)/(*j1))*((*sx2)/(*j2))*((*sx3)/(*j3)));

  // data array
  fprintf(funit,"VECTORS %s float\n",name);

  // data values
  for (k3=0; k3<(*sx3); k3=k3+(*j3)){
     for (k2=-(*sx2)/2; k2<(*sx2)/2; k2=k2+(*j2)){
        i2=((*sx2)+fix(k2)) % (*sx2);

        for (k1=-(*sx1)/2; k1<(*sx1)/2; k1+=(*j1)){
           i1=((*sx1)+fix(k1)) % (*sx1);

#ifdef ALIGN_DATA
           index=i1+(i2+k3*(*sx2))*((*sx1)+2);
#else
           index=i1+(i2+k3*(*sx2))*(*sx1);
#endif

           // convert to big endian if necessary
           buffer[0]=(1u==endian)?u1[index]:swap(u1[index]);
           buffer[1]=(1u==endian)?u2[index]:swap(u2[index]);
           buffer[2]=(1u==endian)?u3[index]:swap(u3[index]);

           fwrite(buffer,12,1,funit);
        }
     }
  }

  // close binary file
  fclose(funit);

}


/*************************************************************
*  subroutine ExportVTK_tensors_Legacy
*  creates a .vtk file in the VTK Legacy binary format with 
*  structured points containing a vector field.
*
*  sylvain barbot 10/28/11 - original form
*************************************************************/
void exportvtk_tensors_legacy_(sig,sx1,sx2,sx3,dx1,dx2,dx3,j1,j2,j3,filename,title,name)
  float *sig;             /* data array for tensor output */
  int *sx1, *sx2, *sx3;   /* number of points             */
  double*dx1, *dx2, *dx3; /* sampling distance            */ 
  int *j1, *j2, *j3;      /* subsampling rate             */
  char *filename;         /* output file name             */
  char *title;            /* output file name             */
  char *name;             /* output file name             */
{

  FILE * funit;
  float buffer[9];
  int i1,i2,k1,k2,k3,index;
  unsigned char endian;
#define DOF 6

  funit=fopen(filename,"wb");
  if (NULL==funit){
     fprintf(stderr,"could not open file %s for vtk output\n",filename);
     fprintf(stderr,"exiting.\n");
     return;
  }

  // find endianness
  endian=isbigendian();

  // writing header of file
  fprintf(funit,"# vtk DataFile Version 3.0\n");
  fprintf(funit,"%s\n",title);
  fprintf(funit,"BINARY\n");
  fprintf(funit,"DATASET STRUCTURED_POINTS\n");

  // structured points grid
  fprintf(funit,"DIMENSIONS %i %i %i\n",(*sx1)/(*j1),(*sx2)/(*j2),(*sx3)/(*j3)); 
  fprintf(funit,"ORIGIN %f %f %f\n",-(*dx1)*((*sx1)/2),-(*dx2)*((*sx2)/2),0.0); 
  fprintf(funit,"SPACING %f %f %f\n",(*dx1)*(*j1),(*dx2)*(*j2),(*dx3)*(*j3)); 

  // data header for this grid
  fprintf(funit,"POINT_DATA %i\n",((*sx1)/(*j1))*((*sx2)/(*j2))*((*sx3)/(*j3)));

  // data array
  fprintf(funit,"TENSORS %s float\n",name);

  // data values
  for (k3=0; k3<(*sx3); k3=k3+(*j3)){
     for (k2=-(*sx2)/2; k2<(*sx2)/2; k2=k2+(*j2)){
        i2=((*sx2)+fix(k2)) % (*sx2);

        for (k1=-(*sx1)/2; k1<(*sx1)/2; k1+=(*j1)){
           i1=((*sx1)+fix(k1)) % (*sx1);

           // index of first stress component
           index=(i1+(i2+k3*(*sx2))*(*sx1))*DOF;

           // convert to big endian if necessary
           buffer[0]=(1u==endian)?sig[index+0]:swap(sig[index+0]);
           buffer[1]=(1u==endian)?sig[index+1]:swap(sig[index+1]);
           buffer[2]=(1u==endian)?sig[index+2]:swap(sig[index+2]);
           buffer[3]=(1u==endian)?sig[index+3]:swap(sig[index+3]);
           buffer[4]=(1u==endian)?sig[index+4]:swap(sig[index+4]);
           buffer[5]=(1u==endian)?sig[index+5]:swap(sig[index+5]);
           buffer[6]=(1u==endian)?sig[index+6]:swap(sig[index+6]);
           buffer[7]=(1u==endian)?sig[index+7]:swap(sig[index+7]);
           buffer[8]=(1u==endian)?sig[index+8]:swap(sig[index+8]);

           // write buffer to disk
           fwrite(buffer,36,1,funit);
        }
     }
  }

  // close binary file
  fclose(funit);

}

