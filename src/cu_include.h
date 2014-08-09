/*-----------------------------------------------------------------------
 * ! Copyright 2013 Sylvain Barbot
 * !
 * ! This file is part of RELAX
 * !
 * ! RELAX is free software: you can redistribute it and/or modify
 * ! it under the terms of the GNU General Public License as published by
 * ! the Free Software Foundation, either version 3 of the License, or
 * ! (at your option) any later version.
 * !
 * ! RELAX is distributed in the hope that it will be useful,
 * ! but WITHOUT ANY WARRANTY; without even the implied warranty of
 * ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * ! GNU General Public License for more details.
 * !
 * ! You should have received a copy of the GNU General Public License
 * ! along with RELAX.  If not, see <http://www.gnu.org/licenses/>.
 * !
 * ! \author Sagar Masuti 
 * !----------------------------------------------------------------------*/

#include "config.h"

#ifndef _CUDA_INCLUDE_
#define _CUDA_INCLUDE_

#include <cuda_runtime.h>

#ifdef USING_CUDA

/* --------------------------------  Definitions ---------------------------------- */


#define FFT_FORWARD -1 
#define FFT_INVERSE  1 
#define NUM_STREAMS 4

/* Reduction */
//#define REDUCTION 1


/* Constants */
#define PI2 6.28318530717958623199592693708837032318

/* Can enable this to print the memory available and required at certain stages before 
   the actual allocation is done */

/* #define GPU_MEMORY_LOG 1 */

/* Prints the some debug messages at certain stages to identify the possible bug at 
   a later point */

/* #define PRINT_DEBUG_INFO 1 */


#define CHECK_CUDA_ERROR(sFunction, Label)      if (cudaSuccess != cuError)                                             \
                                                {                                                                       \
                                                        printf ("%s : Failed  reason is : %s Line number : %d\n",       \
                                                        sFunction, cudaGetErrorString(cuError), __LINE__) ;             \
                                                        goto Label ;                                                    \
                                                }


#define CUDA_FREE_MEM(var)              if (NULL != (var))              \
                                        {                               \
                                                cudaFree ((var)) ;      \
                                                (var) = NULL ;          \
                                        }                               \


#define CHECK_ERROR(sFunction)   cuError = cudaGetLastError () ;        \
                                 if (cudaSuccess != cuError)            \
                                 {                                      \
                                    printf ("%s : There was an error executing  \
                                    and the error is %s\n", sFunction, cudaGetErrorString (cuError)) ; \
                                        return cuError;                 \
                                }                            

/* -------------------------------------------------------------------------------- */




/* ------------------------------ Enumeration types ------------------------------- */

/**
  * This is the enumeration for dimension. 
  *
  */

typedef enum _e_dimension
{
        E_INVALID_DIMENSION = 0,
        E_ONE_DIMENSION,
        E_TWO_DIMENSION,
        E_THREE_DIMENSION,
} E_DIMENSION ;


typedef enum _e_tensor_field
{
        E_INVALID_TENSOR_FIELD = 0,
        E_TENSOR_TAU_TAU,
        E_TENSOR_SIG_TAU,
        E_TENSOR_MOM_MOM,
        E_TENSOR_SIG_MOM,
        E_TENSOR_TAU_MOM
}E_TENSOR_FIELD ;


typedef enum _e_tensor_amp_type
{
        E_INVALID_TENSOR_AMP_TYPE=0,
        E_TENSOR_AMP_MOMENT,
        E_TENSOR_AMP_TAU
}E_TENSOR_AMP_TYPE ;

typedef enum _e_tensor_type
{
        E_INVALID_TENSOR_TYPE=0,
	E_TENSOR_SIG,
	E_TENSOR_MOMENT,
	E_TENSOR_TAU
}E_TENSOR_TYPE ;

typedef enum _e_type
{
	E_INVALID_TYPE=0,
        E_TYPE_U,
        E_TYPE_V
}E_TYPE ;

/* -------------------------------------------------------------------------------- */

__constant__  double DEV_CONST_FIR_1[1] = { 5.000e-01 } ;
__constant__  double DEV_CONST_FIR_7[7] = { 8.77856e-01,
                                -2.81913e-01,
                                +6.22696e-02,
                                +2.82441e-02,
                                -5.09029e-02,
                                +4.20471e-02,
                                -1.59409e-02} ;

__constant__  double DEV_CONST_FIR_14[14] = { 9.739464097198434e-01,
                                 -4.492955962260918e-01,
                                  2.606661503992121e-01,
                                 -1.590778397098753e-01,
                                  9.524605395168785e-02,
                                 -5.279001022321913e-02,
                                  2.452656124714124e-02,
                                 -6.434920307760272e-03,
                                 -4.122947453390886e-03,
                                  9.245789328795669e-03,
                                 -1.060146500976655e-02,
                                  9.786847569837574e-03,
                                 -9.114943973080788e-03,
                                  4.398360884720647e-03} ;



const double CONST_FIR_1[1] = { 5.000e-01 } ;
const double CONST_FIR_7[7] = { 8.77856e-01,
                                -2.81913e-01,
                                +6.22696e-02,
                                +2.82441e-02,
                                -5.09029e-02,
                                +4.20471e-02,
                                -1.59409e-02} ;

const double CONST_FIR_14[14] = { 9.739464097198434e-01,
                                 -4.492955962260918e-01,
                                  2.606661503992121e-01,
                                 -1.590778397098753e-01,
                                  9.524605395168785e-02,
                                 -5.279001022321913e-02,
                                  2.452656124714124e-02,
                                 -6.434920307760272e-03,
                                 -4.122947453390886e-03,
                                  9.245789328795669e-03,
                                 -1.060146500976655e-02,
                                  9.786847569837574e-03,
                                 -9.114943973080788e-03,
                                  4.398360884720647e-03} ;


/* -------------------------------------------------------------------------------- */


/* TENSOR equivalent of fortran code. */
typedef struct _st_tensor
{
        float s11 ;
        float s12 ;
        float s13 ;
        float s22 ;
        float s23 ;
        float s33 ;
} ST_TENSOR ;

typedef struct _st_layer
{
/*      double  dZ ;
        double  dGammaDot0 ;
        double  dStressExp ;
        double  dCohesion ;
        double  dFriction ;*/

        double  z ;
        double  gammadot0 ;
        double  stressexponent ;
        double  cohesion ;
        double  friction ;

} ST_LAYER ;

typedef struct _st_weak
{
/*      double dDgammaDot0 ;
        double dX; 
        double dY ;
        double dZ ;
        double dWidth ;
        double dLength ;
        double dThickness ; 
        double dStrike ;
        double dDip ;*/

        double dgammadot0 ;
        double x;
        double y ;
        double z ;
        double width ;
        double length ;
        double thickness ;
        double strike ;
        double dip ;

}ST_WEAK ;

cudaError_t memcpyUsingStreams (float           *fDest,
                                float           *fSrc,
                                int             iBytes,  
                                cudaMemcpyKind   eDirection,
                                cudaStream_t    *pCuStream) ;


#ifdef PAPI_PROF

extern "C" void papistartprofiling_ (char pcProfName[]) ;
extern "C" void papiendprofiling_ (char    pcName[]) ;

#endif

template <class T>
__global__ void scaling (T              *pCompData,
                         T              fScale,
                         int            iNx,
                         int            iNy,
                         int            iNz) ; 



template <class T>              
__global__ void scaling (T              *pCompData,
                         T              fScale,
                         int            iNx,
                         int            iNy,
                         int            iNz)
{                               
        int             iX ;            
        int             iY ;            
        int             iZ ;    
        unsigned int    iIdx ;  
                                
        iX = threadIdx.x ;              
        iY = blockIdx.y ;               
        iZ = blockIdx.x ;               
                                        
        if ((iX < iNx) && (iY < iNy) && (iZ < iNz))
        {
                iIdx = (((iZ * iNy) + iY) * iNx) + iX ;
                pCompData[iIdx] = pCompData[iIdx] * fScale ;
        }
}
void copyfftmemory (float        *fData1,
                    float        *fData2,
               float        *fData3,
               float        *fData4,
               float        *fData5,
               float        *fData6,
                    int         iSize,
                    int         iSize2,
                    int         iDirection) ;

int destroyPlanForFFT () ;

int createPlanForFFT(int iSx1, 
                                int iSx2,
                                int iSx3) ;

#endif /* USING_CUDA */

#endif /* _CUDA_INCLUDE_ */



