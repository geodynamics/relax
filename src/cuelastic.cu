/*-----------------------------------------------------------------------
! Copyright 2013-2016 Sylvain Barbot
!
! This file is part of RELAX
!
! RELAX is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! RELAX is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with RELAX.  If not, see <http://www.gnu.org/licenses/>.
!
! \author Sagar Masuti 
!----------------------------------------------------------------------*/

/* Main functions called from host are : 
    1) cuinit
    2) custressupdatewrapper
    3) cuviscouseigen
    4) cufrictioneigenstress
    5) cubodyforceswrapper
    6) cudeinit
    7) cusource

  Other functions called from host but not so important/significant
    1) curesetvectors
    2) cutensormemset
    3) copytau
    4) cufieldrep 
    5) cutensorfieldadd
    6) cutensoramp
    7) cucopytraction
    8) cufieldadd
*/


/* This include contains the common macros and definitions */
#include "cuinclude.h"
#include "thrust/extrema.h"
#include <cufft.h>
#include <stdio.h>

/* Switch for enabling/disabling the usage of GPU */
#ifdef USING_CUDA

/* Thrust is used for finding min/max element */
#include <thrust/device_vector.h>

/* Number of points in the filter. Possible values are : 1, 7 and 14 */
#define FILTER_SIZE 1

/*      #define ENABLE_REG_BLOCKING         */ //Dont enable this. 
/*      #define PRINT_DEBUG_INFO        */  
/*      #define STRESS_SHARED_MEM       */ //Dont enable this.

/* Some macros */
#define PI  3.141592653589793115997963468544185161 

#define DEG2RAD 0.01745329251994329547437168059786927

#define MAX_NUM(a,b) (((a) > (b)) ? (a) : (b))

#define MAX3(a,b,c) (MAX_NUM(a, b) > c ? MAX_NUM(a, b) : c)

#define MIN(a, b)  (((a) < (b)) ? (a) : (b))

#define DSIGN(a,b) (((b) > 0) ? (a) : -(a)) 


/* -------------------------------------- global variables ------------------------------------- */ 

            /* Device Variables */
float       *gpV1 = NULL ;          /* Device Pointer. No dereferencing in host. */
float       *gpV2 = NULL ;          /* Device Pointer. No dereferencing in host. */     
float       *gpV3 = NULL ;          /* Device Pointer. No dereferencing in host. */
float       *gpU1 = NULL ;          /* Device Pointer. No dereferencing in host. */
float       *gpU2 = NULL ;          /* Device Pointer. No dereferencing in host. */
float       *gpU3 = NULL ;          /* Device Pointer. No dereferencing in host. */
float       *gpGammadot0 = NULL ;   /* Device Pointer. No dereferencing in host. */
float       *pfDevTract1 = NULL ;   /* Device Pointer. No dereferencing in host. */
float       *pfDevTract2 = NULL ;   /* Device Pointer. No dereferencing in host. */
float       *pfDevTract3 = NULL ;   /* Device Pointer. No dereferencing in host. */
ST_TENSOR   *pstSig = NULL ;        /* Device Pointer. No dereferencing in host. */
ST_LAYER    *pstStruct = NULL ;     /* Device Pointer. No dereferencing in host. */
ST_TENSOR   *pstMoment = NULL ;     /* Device Pointer. No dereferencing in host. */
ST_TENSOR   *pstTau = NULL ;        /* Device Pointer. No dereferencing in host. */
ST_TENSOR_LAYER   *pstPrestress = NULL ;  /* Device Pointer. No dereferencing in host. */
ST_TENSOR   *pstEpsilonik = NULL ;
ST_TENSOR   *pstEpsilonikdot = NULL ;
            /* Host Variables */
int         ihSx1 ;             /* Contains sx1 value*/
int         ihSx2 ;             /* Contains sx2 value*/
int         ihSx3 ;             /* Contains sx3 value*/
int         iLen1 ;             /* Contains the number points in the filter. */
int         iLen2 ;             /* Contains the number points in the filter. */
int         iLen3 ;             /* Contains the number points in the filter. */

            /* Device Constants */
__constant__ double constdKer1 [14] ;
__constant__ double constdKer2 [14] ;
__constant__ double constdKer3 [14] ;

/* --------------------------------------------------------------------------------------------- */



/* -------------------------------------- Funtion declaration ---------------------------------- */ 

int cuOptimalFilter (int *, int *, int *, int, int, int, double, double, double) ;

void cuFreeCudaMemory() ;

int checkMemRequirement(int, int, int) ;

void copygammadot0 (int, int, int, float *, int *) ;

bool cuispresent (void *) ;

__host__ __device__ void cutDot (ST_TENSOR *, double *, double *) ;

__host__ __device__ double cuSum (double *, double *) ;

__host__ __device__ void cuMulSub (double, double *, double *, double *) ;

__host__ __device__ double mycuSinh (double ) ;

__host__ __device__ double cuTensorTrace (ST_TENSOR *) ;

__host__ __device__ void cuIsotrpicStressStrain (ST_TENSOR *, double, double) ;

__host__ __device__ void cuTensorOperate (ST_TENSOR *, void *, char) ;

__host__ __device__ void cuTensorMemset (ST_TENSOR *) ;

__host__ __device__ double cuGauss (double, double) ;

__host__ __device__ double cuOmega (double, double) ;

__host__ __device__ double cuGaussp (double, double ) ;

__host__ __device__ double cuOmegap (double, double) ;

__host__ __device__ void cuTensorDyadProd (ST_TENSOR *, double *, double *) ;

__host__ __device__ double cuTensorNorm (ST_TENSOR *) ;

__host__ __device__ double cuDgGammaDotNot (ST_WEAK *, int, double, double, double, double) ;

__host__ __device__ void cuTensorDeviatoric(ST_TENSOR *, ST_TENSOR *) ;

__host__ __device__ void cuTensorDecompose (ST_TENSOR *, double *, ST_TENSOR *) ;

__host__ __device__ void cuShiftedCoordinates (int, int, int, int, int, int, double, double, 
                                               double, double *, double *, double *) ;

__device__ __host__ void print_tensor (ST_TENSOR *) ;

__device__ void culocalstrain_fir2 (ST_TENSOR *, int, int, int, int, int, int, float *, float *,
                                    float *, int, int, int) ;

__device__ void cuLocalStrain_ani (ST_TENSOR  *, int, int, double, int, int, int, int, int, 
                                   float *, float *, float *, int, int, int) ;

__device__ void cuLocalDivergence_ani (ST_TENSOR *, int, int, double, int, int, int, int, int, 
                                       double *, double *, double *, int, int, int) ;

__device__ void cuLocalDivergence_fir (ST_TENSOR *, int, int, int, int, int, int, double *, 
                                       double *, double *, int, int, int) ;

__global__ void cuTensorAmpKernel (ST_TENSOR *, double *, int, int, int) ;

__global__ void cuTensorFieldKernel (ST_TENSOR *, ST_TENSOR *, float, float, int, int, int) ;

__global__ void cuFieldAddKernel (float *, float *, float, float, int, int, int) ;

__global__ void cuEquivalentTraction (float *, float *, float *, ST_TENSOR *, int, int) ;

__global__ void cuLocalDivergenceKernel (int, int, int, int, int, double, int, int, int, float *, 
                                         float *, float *, ST_TENSOR *) ;

__global__ void cuStressUpdateKernel (int, int, int, double, double, int, int, int, float *, 
                                      float *, float *, ST_TENSOR *) ;

__global__ void cuEquiBodyKernel (ST_TENSOR *, int, int, int, int, int, int, float *, float *, 
                                  float *) ;

__global__ void cuBuildGammadotKernel (int, int, int, double, double, double, 
                                       double, ST_WEAK *, int, float *) ;

__global__ void cuViscousEigenKernel (ST_LAYER *, ST_TENSOR *, ST_TENSOR *, ST_TENSOR_LAYER *,
                                      float *, double, int, int, int, double, double, 
                                      double, float *, float *, bool, bool, bool) ;

__global__ void cuTransientEigenKernel (ST_LAYER *, ST_TENSOR *, ST_TENSOR *, float *,
                                        ST_TENSOR *, ST_TENSOR *, double, int,
                                        int, int, double, double, double, int, float *,
                                        bool) ;

__global__ void cuLocalStressStrainKernel (int, int, int, int, int, double, double, double, int, 
                                           int, int, float *, float *, float *, ST_TENSOR *) ;

__global__ void cuFrictionStress (double, double, double, double, double, double, double, double,
                                  double, double, double, double, double, int, int, int, double, 
                                  double, double, double, double, double, double, double, double, 
                                  double, int, double, float *, ST_TENSOR *, ST_TENSOR *, 
                                  ST_LAYER *) ;

__global__ void cuSourceTractionKernel (int, int, int, double, double, double, double, double, 
                                        double, double, double, double, double, double, double, 
                                        double, double, double, double, double, double, double, 
                                        double, double, double, double, float *, float *, float *) ;

__global__ void cuSourceForceKernel (int, int, int, double, double, double, double, double, double,
                                     double, double, double, double, double, double, double, double,
                                     double, double, double, double, double, double, double, double,
                                     double, double, double, float *, float *, float *) ;

/* --------------------------------------------------------------------------------------------- */


/* ----------------------------------------- Intermediate functions ---------------------------- */

/*
!-----------------------------------------------------------------
  !> StressUpdate
  !! computes the 3-d stress tensor sigma_ij' from the current
  !! deformation field. Strain is the second order tensor
  !!
  !!  \f[ \epsilon_{ij} = \frac{1}{2} ( u_{i,j} + u_{j,i} ) \f]
  !!
  !! The displacement derivatives are approximated numerically by the
  !! application of a differentiator space-domain finite impulse
  !! response filter. Coefficients of the filter can be obtained with
  !! the MATLAB command line
  !!
  !!  \f[ \sigma' = - C' : E \f]
  !!
  !! or in indicial notation
  !!
  !!
  !!  \f[ \sigma_{ij}' = -\lambda'*\delta_{ij}*\epsilon_{kk} - 2*\mu'*\epsilon_{ij}\f]
  !!
  !! where C' is the heterogeneous elastic moduli tensor and lambda'
  !! and mu' are the inhomogeneous lame parameters
  !!
  !!  \f[ C' = C(x) - C_0 \f]
  !!
  !! For isotropic materials
  !!
  !!  \f[ \mu'(x) = \mu(x) - \mu_0 \f]
  !!  \f[ \lambda'(x) = \lambda(x) - \lambda_0 \f]
  !!
  !! Optionally, the surface traction sigma_i3 can be sampled.
  !! 
  !-----------------------------------------------------------------
*/

/**
 *  
 *
 * @param   dLambda[in]     Lame's first parameter
 * @param   dMu[in]         shear modulus or Lame's second parameter. 
 * @param   dDx1[in]        Sampling size in x1(north) direction.
 * @param   dDx2[in]        Sampling size in x2(east) direction.
 * @param   dDx3[in]        Sampling size in x3(down) direction.
 * @param   iSx1[in]        The size of array in x1 direction.
 * @param   iSx2[in]        The size of array in x2 direction.
 * @param   iSx3[in]        The size of array in x3 direction.
 * @param   pstHostSig[in]  Host pointer to the sigma. ** legacy code. Not used. reserved for future **
 * @param   fData?[in]      Host pointer to the data(i.e., v? or u?) ** reserved for future ** 
 * @param   pfInput?[in]    Device pointer to the data( either v? or u?)  
*/

extern "C" void custressupdate_ (double     dLambda,
                                 double     dMu,
                                 double     dDx1,
                                 double     dDx2,
                                 double     dDx3,
                                 int        iSx1,
                                 int        iSx2,
                                 int        iSx3,
                                 ST_TENSOR  *pstHostSig,
                                 float      *fData1,
                                 float      *fData2,
                                 float      *fData3,
                                 float      *pfInput1,
                                 float      *pfInput2,
                                 float      *pfInput3)
{
    int          iInd3 = 0 ;
    int          iInd3p = 0 ;
    int          iInd3m = 0 ;

    double       dPx3 = 0.0;
    cudaError_t  cuError = cudaSuccess ;
    dim3         dimGrid (iSx2, 1, 1) ;
    dim3         dimBlock (iSx1, 1, 1) ;
    dim3         dimGrid1 (iSx3, iSx2, 1) ;
    dim3         dimBlock1 (iSx1, 1, 1) ;

#ifdef PAPI_PROF        
    char        cTimerName[17] = "stress          " ;
#endif

#ifdef PAPI_PROF        
    papistartprofiling_(cTimerName) ;
#endif

    for (iInd3 = 0 ; iInd3 < iSx3 ; iInd3++)
    {
        if ((iInd3 >= iLen3) && (iInd3 < (iSx3-iLen3)))
        {
            continue ;
        }

        if (iInd3 == 0)
        {
            dPx3 = dDx3 ;
            iInd3p = 1 ;
            iInd3m = 0 ;
        }
        else
        {
            if (iInd3 == iSx3-1)
            {
                dPx3 = dDx3 ;
                iInd3p = iSx3-1 ;
                iInd3m = iSx3-2 ;
            }
            else
            {
                dPx3 = dDx3*2.0 ;
                iInd3m = iInd3-1 ;
                iInd3p = iInd3+1 ;
            }
        }
        cuLocalStressStrainKernel<<<dimGrid, dimBlock>>> (iInd3m, iInd3p, iLen1, iLen2, iInd3, dPx3,
                                                          dLambda, dMu, iSx1, iSx2, iSx3, pfInput1, 
                                                          pfInput2, pfInput3, pstSig) ;

        cuError = cudaGetLastError () ;
        CHECK_CUDA_ERROR ("stressupdate kernel launch failed 1\n", STRESS_UPDATE_EXIT_WITH_FREE)

    }
    if (cudaDeviceSynchronize() != cudaSuccess)
    {
        printf ("sync failed 1\n") ;
    }
    
    cuStressUpdateKernel <<<dimGrid1, dimBlock1>>>(iLen1, iLen2, iLen3, dLambda, dMu,
                                                   iSx1, iSx2, iSx3, pfInput1, pfInput2, pfInput3,
                                                   pstSig) ;
    cuError = cudaGetLastError () ;
    CHECK_CUDA_ERROR ("stressupdate kernel launch failed\n", STRESS_UPDATE_EXIT_WITH_FREE)

    if (cudaDeviceSynchronize() != cudaSuccess)
    {
        printf ("sync failed 2\n") ;
    }

#ifdef PAPI_PROF        
    papiendprofiling_(cTimerName) ;
#endif

#ifdef PRINT_DEBUG_INFO
    printf ("custressupdate: exited with no errors\n") ;
#endif

    return ;


STRESS_UPDATE_EXIT_WITH_FREE:
    cuFreeCudaMemory () ;
    return ;
}

/**
  !-----------------------------------------------------------------
  !> subroutine EquivalentBodyForce
  !! computes and updates the equivalent body-force
  !!
  !!         f = - div.( C : E^i )
  !!
  !! and the equivalent surface traction
  !!
  !!         t = n . C : E^i
  !!
  !! with n = (0,0,-1). In indicial notations
  !!
  !!         f_i = - (C_ijkl E^i_kl),j
  !!
  !! and
  !!
  !!         t_1 = n_j C_ijkl E^i_kl
  !!
  !! where f is the equivalent body-force, t is the equivalent surface
  !! traction, C is the elastic moduli tensor and E^i is the moment
  !! density tensor tensor.
  !!
  !! Divergence is computed with a mixed numerical scheme including
  !! centered finite-difference (in the vertical direction) and
  !! finite impulse response differentiator filter for derivatives
  !! estimates. see function 'stress' for further explanations.
  !-----------------------------------------------------------------
 *
 * @param   pstSig      Device pointer to the sigma. 
 * @param       dDx1[in]        Sampling size in x1(north) direction.
 * @param       dDx2[in]        Sampling size in x2(east) direction.
 * @param       dDx3[in]        Sampling size in x3(down) direction.
 * @param       iSx1[in]        The size of array in x1 direction.
 * @param       iSx2[in]        The size of array in x2 direction.
 * @param       iSx3[in]        The size of array in x3 direction.
 * @param       fData?[in]      Host pointer to the data(i.e., v? or u?) ** reserved for future ** 
 * @param       pfT?[in]        Host pointer to the data(t?)  
*/

extern "C" void cuequivalentbodyforces_ (ST_TENSOR  *pstSig,
                                         double     dDx1,
                                         double     dDx2,
                                         double     dDx3,
                                         int        iSx1,
                                         int        iSx2,
                                         int        iSx3,
                                         float      *fData1,
                                         float      *fData2,
                                         float      *fData3,
                                         float      *pfT1,
                                         float      *pfT2,
                                         float      *pfT3)
{
    int          iInd3 = 0 ;
    int          iInd3p = 0 ;
    int          iInd3m = 0 ;

    double       dPx3 = 0.0;
    cudaError_t  cuError = cudaSuccess ;

    dim3         dimGrid (iSx2, 1, 1) ;
    dim3         dimBlock (iSx1, 1, 1) ;
    dim3         dimGrid1 (iSx3, iSx2, 1) ;
    dim3         dimBlock1 (iSx1, 1, 1) ;

    cuEquivalentTraction <<<dimGrid, dimBlock>>> (pfDevTract1, pfDevTract2, pfDevTract3,
                                                  pstSig, iSx1, iSx2 ) ;

    cuError = cudaGetLastError () ;
    CHECK_CUDA_ERROR ("cuequivalentbodyforces_ : Failed in launch of cuEquivalentTraction\n", 
                      BODY_FORCES_EXIT_WITH_FREE) ;

    if (cudaSuccess != cudaDeviceSynchronize ())
    {
        printf ("cuequivalentbodyforces_ : Failed in sync 1\n") ;
        goto BODY_FORCES_EXIT_WITH_FREE ;
    }

    for (iInd3 = 0 ; iInd3 < iSx3 ; iInd3++)
    {
        if ((iInd3 >= iLen3) && (iInd3 < (iSx3-iLen3)))
        {
            continue ;
        }

        if (iInd3 == 0)
        {
            dPx3 = dDx3 ;
            iInd3p = 1 ;
            iInd3m = 0 ;
        }
        else
        {
            if (iInd3 == iSx3-1)
            {
                dPx3 = dDx3 ;
                iInd3p = iSx3-1 ;
                iInd3m = iSx3-2 ;
            }
            else
            {
                dPx3 = dDx3*2.0 ;
                iInd3m = iInd3-1 ;
                iInd3p = iInd3+1 ;
            }
        }
        cuLocalDivergenceKernel<<<dimGrid, dimBlock>>> (iInd3m, iInd3p, iLen1, iLen2, iInd3, dPx3,
                                                        iSx1, iSx2, iSx3, gpV1, gpV2, gpV3, pstSig);


        cuError = cudaGetLastError () ;
        CHECK_CUDA_ERROR ("cuequivalentbodyforces_ : cuLocalDivergenceKernel kernel launch failed\n",
                          BODY_FORCES_EXIT_WITH_FREE)

    }

    dimGrid1.x = iSx3 - 2 * iLen3 ;

    cuEquiBodyKernel <<<dimGrid1, dimBlock1>>> (pstSig, iLen1, iLen2, iLen3, 
                                                iSx1, iSx2, iSx3, gpV1, gpV2, gpV3) ;
    cuError = cudaGetLastError () ;
    CHECK_CUDA_ERROR ("cuequivalentbodyforces_ : cuEquiBodyKernel kernel launch failed\n", 
                      BODY_FORCES_EXIT_WITH_FREE)

    if (cudaSuccess != cudaDeviceSynchronize ())
    {
        printf ("cuequivalentbodyforces_ : Failed in sync 2\n") ;
        goto BODY_FORCES_EXIT_WITH_FREE ;
    }


#ifdef PRINT_DEBUG_INFO
    printf ("cubodyforces: exited with no errors\n") ;
#endif


    return ;


BODY_FORCES_EXIT_WITH_FREE :
    cuFreeCudaMemory () ;
}

/* -------------------------------------- Intermediate functions end ----------------------------------- */



/* -------------------------------- extern functions called from fortran ------------------------------- */
/**
 * This function allocates and initializes various memory required. 
 *
 * @param       iSx1[in]        The size of array in x1 direction.
 * @param       iSx2[in]        The size of array in x2 direction.
 * @param       iSx3[in]        The size of array in x3 direction.
 * @param       dDx1[in]        Sampling size in x1(north) direction.
 * @param       dDx2[in]        Sampling size in x2(east) direction.
 * @param       dDx3[in]        Sampling size in x3(down) direction.
 * @param   iRet[in,out]    Return code for any errors in allocation or initialization.
 **/

extern "C" void cuinit_ (int    iSx1,
                         int    iSx2,
                         int    iSx3,
                         double dDx1,
                         double dDx2,
                         double dDx3,
                         int    *iRet)
{
    cudaError_t cuError = cudaSuccess ;
    int         iSize = 0 ;
    int         iSize2 = 0 ;
#ifdef PRINT_DEBUG_INFO
    size_t iFreeMem ;
    size_t iTotalMem ;
#endif
    int         iDev ;

    cudaDeviceProp deviceProp;
    cuError = cudaGetDevice (&iDev) ;
    cudaGetDeviceProperties(&deviceProp, iDev);
    //printf("Device %d: \"%s\"\n", iDev, deviceProp.name);
    *iRet = 1 ;

    ihSx1 = iSx1 ;
    ihSx2 = iSx2 ;
    ihSx3 = iSx3 ;
    
    if (-1 == checkMemRequirement(iSx1,iSx2,iSx3))
    {
        printf ("********************** ERROR ******************\n") ;
        printf ("Memory required to run on GPU is insufficient\n");
        printf ("Either try reducing the grid size or run on CPU only\n") ;
        printf ("********************** ERROR ******************\n\n") ;
        return ;
    }

    iSize = sizeof (float) * (iSx1 + 2) * iSx2 * iSx3 ;

    cuError = cudaMalloc((void**)&gpV1, iSize) ;
    if (cudaSuccess != cuError)
    {
        printf ("cuinit : Failed to allocate memory 1\n") ;
        goto CUINIT_FAILURE ;
    }

    cuError = cudaMalloc((void**)&gpV2, iSize) ;
    if (cudaSuccess != cuError)
    {
        printf ("cuinit : Failed to allocate memory 2\n") ;
        goto CUINIT_FAILURE ;
    }

    cuError = cudaMalloc((void**)&gpV3, iSize) ;
    if (cudaSuccess != cuError)
    {
        printf ("cuinit : Failed to allocate memory 3\n") ;
        goto CUINIT_FAILURE ;
    }

    cuError = cudaMalloc((void**)&gpU1, iSize) ;
    if (cudaSuccess != cuError)
    {
        printf ("cuinit : Failed to allocate memory 4\n") ;
        goto CUINIT_FAILURE ;
    }

    cuError = cudaMalloc((void**)&gpU2, iSize) ;
    if (cudaSuccess != cuError)
    {
        printf ("cuinit : Failed to allocate memory 5\n") ;
        goto CUINIT_FAILURE ;
    }

    cuError = cudaMalloc((void**)&gpU3, iSize) ;
    if (cudaSuccess != cuError)
    {
        printf ("cuinit : Failed to allocate memory 6\n") ;
        goto CUINIT_FAILURE ;
    }

    iSize2 = sizeof (ST_TENSOR) * iSx1 * iSx2 * (iSx3/2) ;
    cuError = cudaMalloc((void**)&pstSig, iSize2) ;
    if (cudaSuccess != cuError)
    {
        printf ("cuinit : Failed to allocate memory 7\n") ;
        goto CUINIT_FAILURE ;
    }

    cuError = cudaMemset (pstSig, 0, iSize2) ;
    if (cudaSuccess != cuError)
    {
        printf ("cuinit : Failed in memset\n") ;
        goto CUINIT_FAILURE ;
    }

    iSize = sizeof (float) * (iSx1+2) * iSx2 ;
    cuError = cudaMalloc ((void **)&pfDevTract1, iSize) ;
    if (cudaSuccess != cuError)
    {
        printf ("cuinit : Failed to allocate memory 8\n") ;
        goto CUINIT_FAILURE ;
    }

    cuError = cudaMalloc ((void **)&pfDevTract2, iSize) ;
    if (cudaSuccess != cuError)
    {
        printf ("cuinit : Failed to allocate memory 9\n") ;
        goto CUINIT_FAILURE ;
    }
    cuError = cudaMalloc ((void **)&pfDevTract3, iSize) ;
    if (cudaSuccess != cuError)
    {
        printf ("cuinit : Failed to allocate memory 10\n") ;
        goto CUINIT_FAILURE ;
    }

    iSize = sizeof (ST_LAYER) * iSx3/2 ;
    cuError = cudaMalloc ((void **)&pstStruct, iSize) ;
    if (cudaSuccess != cuError)
    {
        printf ("cuinit : Failed to allocate memory 11\n") ;
        goto CUINIT_FAILURE ;
    }
    cuError = cudaMalloc ((void **)&pstMoment, iSize2) ;
    if (cudaSuccess != cuError)
    {
        printf ("cuinit : Failed to allocate memory 12\n") ;
        goto CUINIT_FAILURE ;
    }
    cuError = cudaMalloc((void**)&pstTau, iSize2) ;
    if (cudaSuccess != cuError)
    {
        printf ("cuinit : Failed to allocate memory 13\n") ;
        goto CUINIT_FAILURE ;
    }

    iSize = sizeof (ST_TENSOR_LAYER) * iSx3/2 ;
    cuError = cudaMalloc ((void **)&pstPrestress, iSize) ;
    if (cudaSuccess != cuError)
    {
        printf ("cuinit : Failed to allocate memory 14\n") ;
        goto CUINIT_FAILURE ;
    }

//    if (pstFlags->istransient)
    {
        cuError = cudaMalloc ((void **)&pstEpsilonik, iSize2) ;
        if (cudaSuccess != cuError)
        {
            printf ("cuinit : Failed to allocate memory 15\n") ;
            goto CUINIT_FAILURE ;
        }
        cuError = cudaMalloc ((void **)&pstEpsilonikdot, iSize2) ;
        if (cudaSuccess != cuError)
        {
            printf ("cuinit : Failed to allocate memory 16\n") ;
            goto CUINIT_FAILURE ;
        }
    }

    //memset u1, u2, u3
    iSize = (sizeof (float) * (ihSx1+2) * ihSx2 * ihSx3) ;
    cuError = cudaMemset (gpU1, 0, iSize) ;
    CHECK_CUDA_ERROR ("cuinit_ : Memset failed 1\n", CUINIT_FAILURE)
    cuError = cudaMemset (gpU2, 0, iSize) ;
    CHECK_CUDA_ERROR ("cuinit_ : Memset failed 2\n", CUINIT_FAILURE)
    cuError = cudaMemset (gpU3, 0, iSize) ;
    CHECK_CUDA_ERROR ("cuinit_ : Memset failed 3\n", CUINIT_FAILURE)

    *iRet = cuOptimalFilter (&iLen1, &iLen2, &iLen3,
                             iSx1, iSx2, iSx3, dDx1, dDx2, dDx3) ;
    if (0 != *iRet)
    {
        printf ("custressupdate_ : Something went wrong with optimal filter\n");
        goto CUINIT_FAILURE ;
    }

    createPlanForFFT (iSx1, iSx2, iSx3) ;
    *iRet = 0 ;

#ifdef PRINT_DEBUG_INFO
    cudaMemGetInfo(&iFreeMem, &iTotalMem);
    printf ("cuinit: Memory available after allocation is : %lu MB\n", iFreeMem/(1024*1024));
    printf ("cuinit: Total memory available is : %lu MB\n",iTotalMem/(1024*1024));
    printf ("cuinit: exited with no errors\n") ;
#endif

    return ;

CUINIT_FAILURE:
    cuFreeCudaMemory () ;
}

extern "C" void cuinflags_ (int istransient, 
                            )
{
    
}
/**
 * This is called from the host code. For more information check custressupdate function.
 * @param   eType[in]       Variable to indicate data to be used(i.e., u? and v?) 
 * @param   dLambda[in]     Lame's first parameter
 * @param   dMu[in]         shear modulus or Lame's second parameter. 
 * @param   dDx1[in]        Sampling size in x1(north) direction.
 * @param   dDx2[in]        Sampling size in x2(east) direction.
 * @param   dDx3[in]        Sampling size in x3(down) direction.
 * @param   iSx1[in]        The size of array in x1 direction.
 * @param   iSx2[in]        The size of array in x2 direction.
 * @param   iSx3[in]        The size of array in x3 direction.
 * @param   fData?[in]      Host pointer to the data(i.e., v? or u?) ** reserved for future ** 
 * @param   pstHostSig[in]  Host pointer to the sigma. ** legacy code. Not used. reserved for future **
 *
 **/

extern "C" void custressupdatewrapper_ (E_TYPE     eType,
                                        double     dLambda,
                                        double     dMu,
                                        double     dDx1,
                                        double     dDx2,
                                        double     dDx3,
                                        int        iSx1,
                                        int        iSx2,
                                        int        iSx3,
                                        float      *fData1,
                                        float      *fData2,
                                        float      *fData3,
                                        ST_TENSOR  *pstHostSig)
{

    switch (eType)
    {
        case E_TYPE_U :
        {
            custressupdate_ (dLambda, dMu, dDx1, dDx2, dDx3, iSx1, iSx2, iSx3,
                             pstHostSig, fData1, fData2, fData3, gpU1, gpU2, gpU3) ;
        }
        break ;
        case E_TYPE_V:
        {
            custressupdate_ (dLambda, dMu, dDx1, dDx2, dDx3, iSx1, iSx2, iSx3,
                             pstHostSig, fData1, fData2, fData3, gpV1, gpV2, gpV3) ;
        }
        break ;
        case E_INVALID_TYPE:
        {
            printf ("custressupdatewrapper_: Invalid input\n") ;
        }
    }

}

/**
 * This is called from the host code. For more information check cuequivalentbodyforces_ function.
 * @param       eType[in]       Variable to indicate data to be used(i.e., pstSig/pstMoment) 
 * @param       dDx1[in]        Sampling size in x1(north) direction.
 * @param       dDx2[in]        Sampling size in x2(east) direction.
 * @param       dDx3[in]        Sampling size in x3(down) direction.
 * @param       iSx1[in]        The size of array in x1 direction.
 * @param       iSx2[in]        The size of array in x2 direction.
 * @param       iSx3[in]        The size of array in x3 direction.
 * @param       fData?[in]      Host pointer to the data(i.e., v? or u?) ** reserved for future ** 
 * @param       pfT?[in]        Host pointer to the data(t?)   
 *
 **/


extern "C" void cubodyforceswrapper_ (E_TENSOR_TYPE  eType,
                                      double         dDx1,
                                      double         dDx2,
                                      double         dDx3,
                                      int            iSx1,
                                      int            iSx2,
                                      int            iSx3,
                                      float          *fData1,
                                      float          *fData2,
                                      float          *fData3,
                                      float          *pfT1,
                                      float          *pfT2,
                                      float          *pfT3)
{
    switch (eType)
    {
        case E_TENSOR_SIG:
        {
            cuequivalentbodyforces_ (pstSig, dDx1, dDx2, dDx3, iSx1, iSx2, iSx3,
                                     fData1, fData2, fData3, pfT1, pfT2, pfT3) ;
        }
        break ;
        case E_TENSOR_MOMENT:
        {
            cuequivalentbodyforces_ (pstMoment, dDx1, dDx2, dDx3, iSx1, iSx2, iSx3,
                                     fData1, fData2, fData3, pfT1, pfT2, pfT3) ;
        }
        break ;
        default:
        {
            printf ("cubodyforceswrapper_ : Invalid input\n") ;
        }
    }
}


extern "C" void cusource_ (double  dMu,
                           double  dS,
                           double  dX,
                           double  dY,
                           double  dZ,
                           double  dL,
                           double  dW,
                           double  dStrike,
                           double  dDip,
                           double  dRake,
                           double  dBeta,
                           int     iSx1,
                           int     iSx2,
                           int     iSx3,
                           double  dDx1,
                           double  dDx2,
                           double  dDx3,
                           float   *pfData1,
                           float   *pfData2,
                           float   *pfData3,
                           float   *pfTract1,
                           float   *pfTract2,
                           float   *pfTract3)
{
    double  dcStrike ;
    double  dsStrike ;
    double  dcDip ;
    double  dsDip ;

    double  dCr ;
    double  dSr ;
    double  dX2r ;
    double  dScale ;

    double  dXr ;
    double  dYr ;
    double  dZr ;
    double  dWp ;
    double  dLp ;

    dim3    dimGrid (iSx2, 1, 1) ;
    dim3    dimBlock (iSx1, 1, 1) ;
    dim3    dimGrid1 (iSx3/2, iSx2, 1) ;
    dim3    dimBlock1 (iSx1, 1, 1) ;

    cudaError_t  cuError = cudaSuccess ;

#ifdef PAPI_PROF        
    char    cTimerName[17] = "source          " ;
#endif

#ifdef PAPI_PROF        
    papistartprofiling_(cTimerName) ;
#endif

    dcStrike = cos (dStrike) ;
    dsStrike = sin (dStrike) ;
    dcDip = cos (dDip) ;
    dsDip = sin (dDip) ;
    dCr = cos (dRake) ;
    dSr = sin (dRake) ;
    dScale = -1.0 * (dMu * dS) ;


    dWp = dW * (1.0 + 2.0 * dBeta) / 2.0 ;
    dLp = dL * (1.0 + 2.0 * dBeta) / 2.0 ;

    dX2r = (dcStrike * dX) - (dsStrike * dY) ;
    dXr = (dcDip * dX2r) - (dsDip * dZ) ;
    dYr = (dsStrike * dX) + (dcStrike * dY) ;
    dZr = (dsDip * dX2r) + (dcDip * dZ) ;


    cuSourceTractionKernel <<<dimGrid, dimBlock>>> (iSx1, iSx2, iSx3, dDx1, dDx2, dDx3, dcStrike, 
                                                    dsStrike,
                                                    dcDip, dsDip, dCr, dSr, dScale, dWp, dLp, dX2r,
                                                    dXr, dYr, dZr, dX, dY, dW, dL, dBeta, dMu, dS,
                                                    pfDevTract1, pfDevTract2, pfDevTract3) ;

    cuError = cudaGetLastError () ;
    CHECK_CUDA_ERROR ("cusource_ : Failed in launch of cusourcetractionkernel\n", CUSOURCE_FREE_EXIT) ;

    if (cudaSuccess != cudaDeviceSynchronize())
    {
        printf ("cusource_: Synch failed\n") ;
        goto CUSOURCE_FREE_EXIT ;
    }
        
    cuSourceForceKernel <<<dimGrid1, dimBlock1>>> (iSx1, iSx2, iSx3, dDx1, dDx2, dDx3, dcStrike, 
                                                   dsStrike, dcDip, dsDip, dCr, dSr, dScale, dWp,
                                                   dLp, dX2r, dXr, dYr, dZr, dX, dY, dZ, dDip, dW,
                                                   dL, dBeta, dMu, dS, gpV1, gpV2, gpV3) ;
    cuError = cudaGetLastError () ;
    CHECK_CUDA_ERROR ("cusource_ : Failed in launch of cuSourceForceKernel\n", CUSOURCE_FREE_EXIT) ;

    if (cudaSuccess != cudaDeviceSynchronize())
    {
        printf ("cusource_: Synch failed\n") ;
        goto CUSOURCE_FREE_EXIT ;
    }


#ifdef PAPI_PROF        
    papiendprofiling_(cTimerName) ;
#endif

#ifdef PRINT_DEBUG_INFO
    printf ("cusources : exited with no errors\n") ;
#endif

    return ;

CUSOURCE_FREE_EXIT:
    cuFreeCudaMemory () ;
}

extern "C" void cutransienteigen_ (ST_LAYER          *pStruct,
                                   ST_TENSOR         *pstEpsilonik,
                                   ST_TENSOR         *pstEpsilonikdot,
                                   double            dMu,
                                   int               iSx1,
                                   int               iSx2,
                                   int               iSx3,
                                   double            dDx1,
                                   double            dDx2,
                                   double            dDx3,
                                   int               bMaxwell,
                                   double            *dMaxwell,
                                   int               bdgammadot0,
                                   float             *dgammadot0)
{
    cudaError_t  cuError = cudaSuccess ;
    int          iSize = 0 ;
    float        *devMinArray = NULL ;
    dim3         dimGrid (iSx3, iSx2, 1) ;
    dim3         dimBlock (iSx1, 1, 1) ;
    int          isdgammadot0 = false ;
    int          iRet = 0;
    double       dTemp = 0.0 ;

#ifdef PAPI_PROF
    char cTimerName[17] = "Transienteigen  " ;
    papistartprofiling_ (cTimerName) ;
#endif

    iSize = sizeof (ST_LAYER) * iSx3 ;
    if (bMaxwell)
    {
        cuError = cudaMalloc((void **) &devMinArray, sizeof (float) * iSx1 * iSx2 * iSx3) ;
        if (cudaSuccess != cuError)
        {
            printf ("cutransienteigen_: Failed to allocate 0\n") ;
            cuFreeCudaMemory () ;
        }
    }

    cuError = cudaMemcpy (pstStruct, pStruct, iSize, cudaMemcpyHostToDevice) ;
    if (cudaSuccess != cuError)
    {
        printf ("cutransienteigen_: memcpy failed 1\n") ;
        cuFreeCudaMemory () ;
    }

    /* if dgammadot0 is present then we need to add that to gammadot0 */
    if (bdgammadot0)
    {
        isdgammadot0 = true;
        copygammadot0 (iSx1, iSx2, iSx3, dgammadot0, &iRet) ;
    }

    cuTransientEigenKernel <<<dimGrid, dimBlock>>> (pstStruct, pstSig, pstMoment,
                                                    gpGammadot0, pstEpsilonik,
                                                    pstEpsilonikdot, dMu,
                                                    iSx1, iSx2, iSx3, dDx1, dDx2,
                                                    dDx3, bMaxwell, devMinArray, isdgammadot0) ;

    cuError = cudaGetLastError () ;
    if ((cudaSuccess != cuError) && (cudaSuccess != cudaDeviceSynchronize()))
    {
        printf ("cutransienteigen_: transient kernel failure \n") ;
        cuFreeCudaMemory () ;
    }

    if (bMaxwell)
    {
        thrust::device_ptr<float> dev(devMinArray);
        thrust::device_ptr<float> min = thrust::min_element(dev, dev+(iSx1 * iSx2 * iSx3)) ;

        cuError = cudaGetLastError () ;
        if (cudaSuccess != cuError)
        {
            printf ("cutransienteigen_: Thrust min element failure \n") ;
            cuFreeCudaMemory () ;
        }

        dTemp  =  *min ;
        *dMaxwell = MIN (*dMaxwell, dTemp) ;
        cudaFree (devMinArray) ;
    }

#ifdef PAPI_PROF
    papiendprofiling_ (cTimerName) ;
#endif
    return ;

}

extern "C" void cutransienteigenwrapper_ (E_TENSOR_TYPE eType,
                                          ST_LAYER   *pStruct,
                                          double     dMu,
                                          int        iSx1,
                                          int        iSx2,
                                          int        iSx3,
                                          double     dDx1,
                                          double     dDx2,
                                          double     dDx3,
                                          int        bMaxwell,
                                          double     *dMaxwell,
                                          int        bgammadot,
                                          float      *dGammadot0 = NULL)
{

    switch (eType)
    {
        case E_TENSOR_IK :
        {
            cutransienteigen_ (pStruct, pstEpsilonik, pstEpsilonikdot, dMu,
                               iSx1, iSx2, iSx3, dDx1, dDx2, dDx3,
                               bMaxwell, dMaxwell, bgammadot, dGammadot0) ;
        }
        break ;
        case E_TENSOR_IKDOT:
        {
            cutransienteigen_ (pStruct, pstEpsilonikdot, pstEpsilonikdot, dMu,
                               iSx1, iSx2, iSx3, dDx1, dDx2, dDx3,
                               bMaxwell, dMaxwell, bgammadot, dGammadot0) ;
        }
        break ;
        case E_INVALID_TYPE:
        {
            printf ("custressupdatewrapper_: Invalid input\n") ;
        }
    }
}
extern "C" void cuviscouseigen_ (ST_LAYER          *pStruct,
                                 ST_TENSOR         *pSig,
                                 ST_TENSOR         *pMoment,
                                 ST_TENSOR_LAYER   *pPrestress, 
                                 double            dMu,
                                 int               iSx1,
                                 int               iSx2,
                                 int               iSx3,
                                 double            dDx1,
                                 double            dDx2,
                                 double            dDx3,
                                 float             *dMaxwell,
                                 float             *dgammadot0,
                                 float             *pGamma)
{
    cudaError_t  cuError = cudaSuccess ;
    int          iSize = 0 ;
    int          iSize1 = 0 ;
    float        *devMinArray = NULL ;
    dim3         dimGrid (iSx3, iSx2, 1) ;
    dim3         dimBlock (iSx1, 1, 1) ;
    bool         isdgammadot0 = false ;
    bool         bMaxwell = false;
    bool         bGamma = false;
    int          iRet = 0; 

#ifdef PAPI_PROF
    char cTimerName[17] = "Eigenstress     " ;
    papistartprofiling_ (cTimerName) ;
#endif

    iSize = sizeof (ST_LAYER) * iSx3 ;
    iSize1 = sizeof (float) * (iSx1 + 2) * iSx2 * iSx3 * 2 ;

    /* Check if maxwell time is present then we need to do a reduction */
    if (cuispresent(dMaxwell))
    {
        bMaxwell = true;
        cuError = cudaMalloc((void **) &devMinArray, sizeof (float) * iSx1 * iSx2 * iSx3) ;
        CHECK_CUDA_ERROR ("cuviscouseigen_ : Failed to allocate 0\n", VISCOUS_FREE_EXIT) ;
    }

    cuError = cudaMemcpy (pstStruct, pStruct, iSize, cudaMemcpyHostToDevice) ;
    CHECK_CUDA_ERROR ("cuviscouseigen_ : memcpy failed 1\n", VISCOUS_FREE_EXIT) ;
        
    iSize = sizeof (ST_TENSOR_LAYER) * iSx3 ;
    cuError = cudaMemcpy (pstPrestress, pPrestress, iSize, cudaMemcpyHostToDevice) ;
    CHECK_CUDA_ERROR ("cuviscouseigen_ : memcpy failed 2\n", VISCOUS_FREE_EXIT) ;

    if (cuispresent(pGamma))
    {
        bGamma = true ;
        cuError = cudaMemset (gpV1, 0, iSize1) ;
        CHECK_CUDA_ERROR ("cuviscouseigen_ : memset failed 1\n", VISCOUS_FREE_EXIT) ;
    }

    /* if dgammadot0 is present then we need to add that to gammadot0 */
    if (cuispresent(dgammadot0))
    {
        isdgammadot0 = true;
        copygammadot0 (iSx1, iSx2, iSx3, dgammadot0, &iRet) ;
    }

    cuViscousEigenKernel <<<dimGrid, dimBlock>>> (pstStruct, pstSig, pstMoment, 
                                                  pstPrestress, gpGammadot0, dMu, 
                                                  iSx1, iSx2, iSx3, dDx1, dDx2, 
                                                  dDx3, devMinArray, gpV1, bMaxwell,
                                                  bGamma, isdgammadot0) ;


    cuError = cudaGetLastError () ;
    CHECK_CUDA_ERROR ("cuviscouseigen_ : Kernel launch failed\n", VISCOUS_FREE_EXIT) ;
    if (cudaSuccess != cudaDeviceSynchronize())
    {
        printf ("cuviscouseigen_: sync failed \n") ;
    }

    if (bMaxwell)
    {
        thrust::device_ptr<float> dev(devMinArray);
        thrust::device_ptr<float> min = thrust::min_element(dev, dev+(iSx1 * iSx2 * iSx3)) ;

        cuError = cudaGetLastError () ;
        CHECK_CUDA_ERROR ("cuviscouseigen_ : Thrust min element failure \n", VISCOUS_FREE_EXIT) ;

        *dMaxwell =  (float)*min ;
        cudaFree (devMinArray) ;
    }


#ifdef PAPI_PROF
    papiendprofiling_ (cTimerName) ;
#endif
    return ;

VISCOUS_FREE_EXIT:
    cuFreeCudaMemory () ;
}

extern "C" void cutensorfieldadd_ (E_TENSOR_FIELD  eField,
                                   int             iSx1,
                                   int             iSx2,
                                   int             iSx3,
                                   float           fC1,
                                   float           fC2)
{
    dim3 dimGrid (iSx3, iSx2, 1) ;
    dim3 dimBlock (iSx1, 1, 1) ;
#ifdef PAPI_PROF
    char cTimerName[17] = "tensorfieldadd  " ;
    papistartprofiling_ (cTimerName) ;
#endif

    switch (eField)
    {
        case E_TENSOR_TAU_TAU :
        {
            cuTensorFieldKernel <<<dimGrid, dimBlock>>> (pstTau, pstTau, fC1, fC2, iSx1, iSx2, iSx3) ;
            if (cudaSuccess != cudaDeviceSynchronize())
            {
                printf ("Failed to sync 1\n") ;
            }
        }
        break ;
        case E_TENSOR_SIG_TAU :
        {
            cuTensorFieldKernel <<<dimGrid, dimBlock>>> (pstSig, pstTau, fC1, fC2, iSx1, iSx2, iSx3) ;
            if (cudaSuccess != cudaDeviceSynchronize())
            {
                printf ("Failed to sync 2\n") ;
            }
        }
        break ;
        case E_TENSOR_MOM_MOM :
        {
            cuTensorFieldKernel <<<dimGrid, dimBlock>>> (pstMoment, pstMoment, fC1, fC2, iSx1, 
                                                         iSx2, iSx3) ;
            if (cudaSuccess != cudaDeviceSynchronize())
            {
                printf ("Failed to sync 3\n") ;
            }
        }
        break ;
        case E_TENSOR_SIG_MOM :
        {
            cuTensorFieldKernel <<<dimGrid, dimBlock>>> (pstSig, pstMoment, fC1, fC2, iSx1, iSx2, iSx3) ;
            if (cudaSuccess != cudaDeviceSynchronize())
            {
                printf ("Failed to sync 4\n") ;
            }
        }
        break ;
        case E_TENSOR_TAU_MOM :
        {
            cuTensorFieldKernel <<<dimGrid, dimBlock>>> (pstTau, pstMoment, fC1, fC2, iSx1, iSx2, iSx3) ;
            if (cudaSuccess != cudaDeviceSynchronize())
            {
                printf ("Failed to sync 5\n") ;
            }
        }
        break ;
        case E_TENSOR_IK_IKDOT:
        {
            cuTensorFieldKernel <<<dimGrid, dimBlock>>> (pstEpsilonik, pstEpsilonikdot, fC1, fC2, iSx1, iSx2, iSx3) ;
            if (cudaSuccess != cudaDeviceSynchronize())
            {
                printf ("Failed to sync 6\n") ;
            }
        }
        break ;
        case E_TENSOR_IKDOT_IK:
        {
            cuTensorFieldKernel <<<dimGrid, dimBlock>>> (pstEpsilonikdot, pstEpsilonik, fC1, fC2, iSx1, iSx2, iSx3) ;
            if (cudaSuccess != cudaDeviceSynchronize())
            {
                printf ("Failed to sync 6\n") ;
            }
        }
        break ;
        case E_TENSOR_IKDOT_IKDOT:
        {
            cuTensorFieldKernel <<<dimGrid, dimBlock>>> (pstEpsilonikdot, pstEpsilonikdot, fC1, fC2, iSx1, iSx2, iSx3) ;
            if (cudaSuccess != cudaDeviceSynchronize())
            {
                printf ("Failed to sync 7\n") ;
            }
        }
        break ; 
        default :
        {
            printf ("cutensorfieldadd_ : The enum sent to this function is wrong %d\n",eField) ;
            return ;
        }
    }

    if (cudaSuccess != cudaGetLastError())
    {
        printf ("cutensorfieldadd_ : Something went wrong in tensor kernel\n") ;
        cuFreeCudaMemory () ;
    }
    if (cudaSuccess != cudaDeviceSynchronize())
    {
        printf ("cutensorfieldadd_ : Failed to sync\n") ;
        cuFreeCudaMemory () ;
    }

#ifdef PAPI_PROF
    papiendprofiling_ (cTimerName) ;
#endif

    return ;
}

extern "C" void cucopytraction_ (float *pTract,
                                 int    iSx1,
                                 int    iSx2,
                                 int   *iRet)
{
    *iRet = 0 ;
    if (cudaSuccess != cudaMemcpy (pTract, pfDevTract3, sizeof (float) * (iSx1+2) * iSx2, 
                                   cudaMemcpyDeviceToHost))
    {
        printf ("Error in memcpy\n") ;
        *iRet = 1 ;
    }
}

extern "C" void cudeinit_ ()
{
    destroyPlanForFFT() ;
    cuFreeCudaMemory () ;
}

extern "C" void cutensoramp_ (E_TENSOR_AMP_TYPE  eType,
                              int                iSx1,
                              int                iSx2,
                              int                iSx3,
                              double             *dAmp)
{
    cudaError_t     cuError = cudaSuccess ;
    dim3        dimGrid (iSx3, iSx2, 1) ;
    dim3        dimBlock (iSx1, 1, 1) ;
    double      *pdTemp = NULL ;

    cuError = cudaMalloc((void **) &pdTemp, sizeof (double) * iSx1 * iSx2 * iSx3) ;
    if (cuError != cudaSuccess)
    {
        printf ("cutensoramp_ : Failed to allocate 0\n") ;
    }


    switch (eType)
    {
        case E_TENSOR_AMP_MOMENT :
        {
            cuTensorAmpKernel <<<dimGrid, dimBlock>>> (pstMoment, pdTemp, iSx1, iSx2, iSx3) ;
        }
        break ;
        case E_TENSOR_AMP_TAU :
        {
            cuTensorAmpKernel <<<dimGrid, dimBlock>>> (pstTau, pdTemp, iSx1, iSx2, iSx3) ;
        }
        break ;
        default :
        {
            printf ("The enum sent to this function is wrong\n") ;
            return ;
        }
    }
    if (cudaSuccess != cudaGetLastError())
    {
        printf ("cutensorfieldadd_ : Something went wrong in tensor kernel\n") ;
        cuFreeCudaMemory () ;
    }
    if (cudaSuccess != cudaDeviceSynchronize())
    {
        printf ("SYnch failed\n") ;
    }

    thrust::device_ptr<double> dev(pdTemp);
    double dSum = thrust::reduce(dev, dev+(iSx1 * iSx2 * iSx3), (double) 0.0, thrust::plus<double>());
    cuError = cudaGetLastError () ;
    CHECK_CUDA_ERROR ("cutensoramp_ : Thrust reduce failure \n", TENSOR_EXIT) ;

    *dAmp = dSum ;

    cudaFree (pdTemp) ;

    return ;

TENSOR_EXIT:
    cudaFree (pdTemp) ;
    cuFreeCudaMemory () ;
}

extern "C" void cufieldadd_ (E_TYPE eType,
                             float *pfVal1,
                             float *pfVal2,
                             float *pfVal3,
                             float *pfVal4,
                             float *pfVal5,
                             float *pfVal6,
                             int   iSx1,
                             int   iSx2,
                             int   iSx3,
                             float fC1,
                             float fC2)
{
    dim3 dimGrid (iSx3, iSx2, 1) ;
    dim3 dimBlock (iSx1/2, 1, 1) ;
    cudaError_t cuError = cudaSuccess ;

    switch (eType)
    {
        case E_TYPE_U :
        {
            cuFieldAddKernel <<<dimGrid, dimBlock>>> (gpU1, gpV1, fC1, fC2, iSx1, iSx2, iSx3) ;
            cuError = cudaGetLastError () ;
            CHECK_CUDA_ERROR ("cufieldadd_: Kernel launch failure 1\n", FIELD_ADD_EXIT) ;

            cuFieldAddKernel <<<dimGrid, dimBlock>>> (gpU2, gpV2, fC1, fC2, iSx1, iSx2, iSx3) ;
            cuError = cudaGetLastError () ;
            CHECK_CUDA_ERROR ("cufieldadd_: Kernel launch failure 2\n", FIELD_ADD_EXIT) ;

            cuFieldAddKernel <<<dimGrid, dimBlock>>> (gpU3, gpV3, fC1, fC2, iSx1, iSx2, iSx3) ;
            cuError = cudaGetLastError () ;
            CHECK_CUDA_ERROR ("cufieldadd_: Kernel launch failure 3\n", FIELD_ADD_EXIT) ;

            if (cudaSuccess != cudaDeviceSynchronize())
            {
                printf ("cufieldadd: failed to synchronize\n") ;
            }
        }
        break ;
        case E_TYPE_V :
        {
            cuFieldAddKernel <<<dimGrid, dimBlock>>> (gpV1, gpU1, fC1, fC2, iSx1, iSx2, iSx3) ;
            cuError = cudaGetLastError () ;
            CHECK_CUDA_ERROR ("cufieldadd_: Kernel launch failure 4\n", FIELD_ADD_EXIT) ;

            cuFieldAddKernel <<<dimGrid, dimBlock>>> (gpV2, gpU2, fC1, fC2, iSx1, iSx2, iSx3) ;
            cuError = cudaGetLastError () ;
            CHECK_CUDA_ERROR ("cufieldadd_: Kernel launch failure 5\n", FIELD_ADD_EXIT) ;

            cuFieldAddKernel <<<dimGrid, dimBlock>>> (gpV3, gpU3, fC1, fC2, iSx1, iSx2, iSx3) ;
            cuError = cudaGetLastError () ;
            CHECK_CUDA_ERROR ("cufieldadd_: Kernel launch failure 6\n", FIELD_ADD_EXIT) ;

            if (cudaSuccess != cudaDeviceSynchronize())
            {
                printf ("cufieldadd: failed to synchronize\n") ;
            }
        }
        break ;
        default :
        {
            printf ("Invalid input\n") ;
        }
    }
    return ;

FIELD_ADD_EXIT:
    printf ("Something wrong in cufieldadd\n") ;
    cuFreeCudaMemory () ;
}

extern "C" void cufieldrep_ (int iSx1,
                             int iSx2,
                             int iSx3)
{
    float *pfdata = NULL ;
    cufieldadd_ (E_TYPE_U, pfdata, pfdata, pfdata, pfdata, pfdata, pfdata, iSx1, iSx2, iSx3, 0.0, 1.0) ;
        
    return ;
}

extern "C" void curesetvectors_ ()
{
    int iSize = 0 ;
    cudaError_t cuError = cudaSuccess ;

    iSize = (sizeof (float) * (ihSx1+2) * ihSx2 * ihSx3) ;
    cuError = cudaMemset (gpV1, 0, iSize) ;
    CHECK_CUDA_ERROR ("curesetvectors_ : Memset failed 1\n", CURESET_FAILURE)
    cuError = cudaMemset (gpV2, 0, iSize) ;
    CHECK_CUDA_ERROR ("curesetvectors_ : Memset failed 2\n", CURESET_FAILURE)
    cuError = cudaMemset (gpV3, 0, iSize) ;
    CHECK_CUDA_ERROR ("curesetvectors_ : Memset failed 3\n", CURESET_FAILURE)

    iSize = (sizeof (float) * (ihSx1+2) * ihSx2) ;
    cuError = cudaMemset (pfDevTract1, 0, iSize) ;
    CHECK_CUDA_ERROR ("curesetvectors_ : Memset failed 4\n", CURESET_FAILURE)
    cuError = cudaMemset (pfDevTract2, 0, iSize) ;
    CHECK_CUDA_ERROR ("curesetvectors_ : Memset failed 5\n", CURESET_FAILURE)
    cuError = cudaMemset (pfDevTract3, 0, iSize) ;
    CHECK_CUDA_ERROR ("curesetvectors_ : Memset failed 6\n", CURESET_FAILURE)

    return  ;

CURESET_FAILURE:
    cuFreeCudaMemory () ;

}

bool cuispresent (void *pVar = NULL)
{
    int ipresent=0;
    __util_MOD_ispresent(pVar, &ipresent) ;
    return (ipresent == 1) ? true : false ;
}

void copygammadot0 (int        iSx1,
                    int        iSx2,
                    int        iSx3,
                    float      *fGammadot0,
                    int        *iRet)
{
    int         iSize = 0 ;
    cudaError_t cuError = cudaSuccess ;

    *iRet = 1 ; 
    iSize = sizeof (float) * iSx1 * iSx2 * iSx3 ;
    //allocate 
    if (NULL == gpGammadot0)
    { 
        cuError = cudaMalloc((void**)&gpGammadot0, iSize) ;
        if (cudaSuccess != cuError)
        {
            printf ("copygammadot0_ : Failed to allocate memory 1\n") ;
            goto COPY_GAMMA_DOT_0; 
        }
    }

    if (cudaSuccess != cudaMemcpy (gpGammadot0, fGammadot0, iSize, 
                                   cudaMemcpyHostToDevice))
    {
            printf ("copygammadot0_ : failed in memcpy 1\n") ;
            goto COPY_GAMMA_DOT_0;
    }

    *iRet = 0 ; 
    return ;

COPY_GAMMA_DOT_0:
    cuFreeCudaMemory () ;

}

extern "C" void copytau_ (ST_TENSOR  *pTemp,
                          int        iSx1,
                          int        iSx2,
                          int        iSx3,
                          int        iForward)
{
    if (1 == iForward)
    {
        if (cudaSuccess != cudaMemcpy (pstTau, pTemp, sizeof(ST_TENSOR)*iSx1*iSx2*iSx3, 
                                       cudaMemcpyHostToDevice))
        {
                printf ("copytau_ : failed in memcpy 1\n") ;
        }
    }
    else
    {
        if (cudaSuccess != cudaMemcpy (pTemp, pstTau, sizeof(ST_TENSOR)*iSx1*iSx2*iSx3, 
                                       cudaMemcpyDeviceToHost))
        {
            printf ("copytau_ : failed in memcpy 2\n") ;
        }
    }
}

extern "C" void cutensormemset_ (E_TENSOR_TYPE eType)
{
    cudaError_t cuError = cudaSuccess ;

    switch (eType)
    {
        case E_TENSOR_MOMENT:
        {
            cuError = cudaMemset (pstMoment, 0, sizeof (ST_TENSOR) * ihSx1 * ihSx2 * ihSx3/2) ;
        }
        break ;
        case E_TENSOR_SIG:
        {
            cuError = cudaMemset (pstSig, 0, sizeof (ST_TENSOR) * ihSx1 * ihSx2 * ihSx3/2) ;
        }
        break ;
        case E_TENSOR_TAU:
        {
            cuError = cudaMemset (pstTau, 0, sizeof (ST_TENSOR) * ihSx1 * ihSx2 * ihSx3/2) ;
        }
        break ;
        case E_INVALID_TENSOR_TYPE:
        {   
            printf ("Invalid input\n") ;
        }
        break ;
    }
    if (cudaSuccess != cuError)
    {
        printf ("cutensormemset: Error in cudamemset\n") ;
    }
}

extern "C" void cuexportspatial_ (float *pTemp1,
                                  float *pTemp2,
                                  float *pTemp3,
                                  int   iIndex3)
{
    cudaError_t cuError = cudaSuccess ;
    int iIndex = 0 ;
    iIndex = iIndex3*ihSx2*(ihSx1+2) ;

    cuError = cudaMemcpy (pTemp1, &gpU1[iIndex], (ihSx1+2)*ihSx2*sizeof (float), cudaMemcpyDeviceToHost) ;
    CHECK_CUDA_ERROR ("cuexportspatial_ : Error in memcpy 1\n", CUEXPORT_SPATIAL)

    cuError = cudaMemcpy (pTemp2, &gpU2[iIndex], (ihSx1+2)*ihSx2*sizeof (float), cudaMemcpyDeviceToHost) ;
    CHECK_CUDA_ERROR ("cuexportspatial_ : Error in memcpy 2\n", CUEXPORT_SPATIAL)

    cuError = cudaMemcpy (pTemp3, &gpU3[iIndex], (ihSx1+2)*ihSx2*sizeof (float), cudaMemcpyDeviceToHost) ;
    CHECK_CUDA_ERROR ("cuexportspatial_ : Error in memcpy 3\n", CUEXPORT_SPATIAL)

    return ;

CUEXPORT_SPATIAL :
    printf ("Error\n") ;
}
extern "C" void cuexportpoints_ (double         *pfU1,
                                 double         *pfU2,
                                 double         *pfU3,
                                 ST_TENSOR      *pstTemp,
                                 int            iInd1,
                                 int            iInd2,
                                 int            iInd3)
{
    int iIndex1 = 0 ;
    int iIndex2 = 0 ;
    cudaError_t cuError = cudaSuccess  ;
    float temp1 = 0.0 ;

    iIndex1 = (((iInd3 * ihSx2) + iInd2) * (ihSx1 + 2)) + iInd1 ;
    iIndex2 = (((iInd3 * ihSx2) + iInd2) * ihSx1) + iInd1 ;

    cuError = cudaMemcpy (&temp1, &gpU1[iIndex1], sizeof (float), cudaMemcpyDeviceToHost) ;
    CHECK_CUDA_ERROR ("cuexportpoints_ : Error in memcpy 1\n", CUEXPORT_EXIT)
    *pfU1 = (double) temp1 ;

    cuError = cudaMemcpy (&temp1, &gpU2[iIndex1], sizeof (float), cudaMemcpyDeviceToHost) ;
    CHECK_CUDA_ERROR ("cuexportpoints_ : Error in memcpy 2\n", CUEXPORT_EXIT)
    *pfU2 = (double) temp1 ;

    cuError = cudaMemcpy (&temp1, &gpU3[iIndex1], sizeof (float), cudaMemcpyDeviceToHost) ;
    CHECK_CUDA_ERROR ("cuexportpoints_ : Error in memcpy 3\n", CUEXPORT_EXIT)
    *pfU3 = (double) temp1 ;

    cuError = cudaMemcpy (pstTemp, &pstSig[iIndex2], sizeof (ST_TENSOR), cudaMemcpyDeviceToHost) ;
    CHECK_CUDA_ERROR ("cuexportpoints_ : Error in memcpy 4\n", CUEXPORT_EXIT)

    return ;

CUEXPORT_EXIT:
    printf ("Error\n") ;

}
extern "C" void cubuildgammadot_ (int     iSx1,
                                  int     iSx2,
                                  int     iSx3,
                                  double  dDx1,
                                  double  dDx2,
                                  double  dDx3,
                                  int     iNz,
                                  ST_WEAK *pDuctile,
                                  float   *pDgammadot0)
{
    cudaError_t  cuError = cudaSuccess ;
    int          iSize = 0 ;
    int          iSize1 = 0 ;
    dim3         dimGrid (iSx3, iSx2, 1) ;
    dim3         dimBlock (iSx1, 1, 1) ;
    double       dBeta = 0 ;
    float        *pfDgammadot0 = NULL ;
    ST_WEAK      *pstDuctile = NULL ;
    
    iSize = sizeof (ST_WEAK) * iNz ;
    iSize1 = sizeof (float) * iSx1 * iSx2 * iSx3 ; 

    cuError = cudaMalloc((void **) &pfDgammadot0, iSize1) ;
    if (cudaSuccess != cuError)
    {
        printf ("cubuildgammadot_ : Couldnt allocate memory 1\n") ;
        cuFreeCudaMemory() ;
    }

    cuError = cudaMalloc((void **) &pstDuctile, iSize) ;
    if (cudaSuccess != cuError)
    {
        printf ("cubuildgammadot_ : Couldnt allocate memory 2\n") ;
        cuFreeCudaMemory() ;
    }

    cuError = cudaMemcpy (pstDuctile, pDuctile, iSize, cudaMemcpyHostToDevice) ;
    if (cudaSuccess != cuError)
    {
        printf ("cubuildgammadot_ : Couldnt copy memory 3\n") ;
        cuFreeCudaMemory() ;
    }

    //call kernel
    cuBuildGammadotKernel <<<dimGrid, dimBlock>>> (iSx1, iSx2, iSx3, dDx1, dDx2,
                                                   dDx3, dBeta, pstDuctile, iNz, 
                                                   pfDgammadot0) ;

    cuError = cudaGetLastError () ;
    if (cudaSuccess != cuError)
    {
        printf ("cubuildgammadot0_ : Kernel launch failed\n") ;
    }

    if (cudaSuccess != cudaDeviceSynchronize())
    {
        printf ("cubuildgammadot0_ : sync failed \n") ;
    }

    //copy back

    cuError = cudaMemcpy (pDgammadot0, pfDgammadot0, iSize1, cudaMemcpyDeviceToHost) ;
    if (cudaSuccess != cuError)
    {
        printf ("cubuildgammadot_ : Couldnt copy memory 4\n") ;
        cuFreeCudaMemory() ;
    }

    cudaFree(pfDgammadot0);
    cudaFree(pstDuctile);

    return ; 
}

extern "C" void cufrictioneigenstress_ (double     dX,
                                        double     dY,
                                        double     dZ,
                                        double     dL,
                                        double     dW,
                                        double     dStrike,
                                        double     dDip,
                                        double     dRake,
                                        double     dBeta,
                                        double     dMu,
                                        ST_LAYER   *pStruct,
                                        int        iSx1,
                                        int        iSx2,
                                        int        iSx3,
                                        double     dDx1,
                                        double     dDx2,
                                        double     dDx3,
                                        int        bPresent,
                                        float      *dMaxwell,
                                        ST_TENSOR  *pMoment,
                                        ST_TENSOR  *pSig)
{
    double      dScaling ;
    double      dCstrike ;
    double      dSstrike ;
    double      dCdip ;
    double      dSdip ;
    double      dCr ;
    double      dSr ;
    double      dLp ;
    double      dWp ;
    double      dX2r ;
    double      dXr ;
    double      dYr ;
    double      dZr ;
    float       *devMinArray ;
    dim3        dimGrid (iSx3, iSx2, 1) ;
    dim3        dimBlock (iSx1, 1, 1) ;
    cudaError_t cuError ;
    int         iSize ;

    cuError = cudaSuccess ;
    dScaling = sqrt(PI2) * dDx1 ;

    dCstrike = cos(dStrike) ;
    dSstrike = sin(dStrike) ;
    dCdip = cos(dDip) ;
    dSdip = sin(dDip) ;
    dCr = cos(dRake) ;
    dSr = sin(dRake) ;

    dWp = dW*(1.0 + 2.0 * dBeta)/2.0 ;
    dLp = dL*(1.0 + 2.0 * dBeta)/2.0 ;

    dX2r = dCstrike * dX  - dSstrike * dY ;
    dXr = dCdip * dX2r - dSdip * dZ ;
    dYr = dSstrike * dX + dCstrike * dY ;
    dZr = dSdip * dX2r + dCdip * dZ ;

    if (1 == bPresent)
    {
        cuError = cudaMalloc((void **) &devMinArray, sizeof (float) * iSx1 * iSx2 * iSx3) ;
        CHECK_CUDA_ERROR ("cufrictioneigenstress_ : Failed to allocate 0\n", FRICTION_FREE_EXIT) ;
    }
    iSize = sizeof (ST_LAYER) * iSx3 ;
    cuError = cudaMemcpy (pstStruct, pStruct, iSize, cudaMemcpyHostToDevice) ;
    CHECK_CUDA_ERROR ("cufrictioneigenstress_ : memcpy failed 1\n", FRICTION_FREE_EXIT) ;

    cuFrictionStress <<<dimGrid, dimBlock>>> (dScaling, dCstrike, dSstrike, dCdip, dSdip, dCr, dSr,
                                              dWp, dLp, dX2r, dXr, dYr, dZr, iSx1, iSx2, iSx3, dDx1,
                                              dDx2, dDx3, dX, dY, dZ, dL, dW, dRake, dMu, bPresent,  
                                              dBeta, devMinArray, pstMoment, pstSig, pstStruct) ;

    cuError = cudaGetLastError () ;
    CHECK_CUDA_ERROR ("cufrictioneigenstress_: Kernel launch failed\n", FRICTION_FREE_EXIT) ;
    if (cudaSuccess != cudaDeviceSynchronize())
    {
        printf ("cufrictioneigenstress_ : sync failed \n") ;
    }

    if (1 == bPresent)
    {
        thrust::device_ptr<float> dev(devMinArray);
        thrust::device_ptr<float> min = thrust::min_element(dev, dev+(iSx1 * iSx2 * iSx3)) ;

        cuError = cudaGetLastError () ;
        CHECK_CUDA_ERROR ("cufrictioneigenstress_ : Thrust min element failure \n", FRICTION_FREE_EXIT) ;

        *dMaxwell =  *min ;
        cudaFree (devMinArray) ;
    }

FRICTION_FREE_EXIT:
    return ;
}
/* -------------------------------------- utility functions ----------------------------------- */

cudaError_t copyFilter1(double dDx,
                        int    iWhich) 
{
    cudaError_t cuError = cudaSuccess ;
    double ker ;
    
    memcpy (&ker, CONST_FIR_1, sizeof (double)) ;
    ker /= dDx ;

    switch (iWhich)
    {
        case 1 :
        {
            cuError = cudaMemcpyToSymbol (constdKer1, &ker, sizeof(ker)) ;
        }
        break ;
        case 2 : 
        {
            cuError = cudaMemcpyToSymbol (constdKer2, &ker, sizeof(ker)) ;
        }
        break ;
        case 3 : 
        {
            cuError = cudaMemcpyToSymbol (constdKer3, &ker, sizeof(ker)) ;
        }
        break ;
    }
    return cuError ; 
}

cudaError_t copyFilter7(double dDx,
                        int    iWhich) 
{
    cudaError_t cuError = cudaSuccess ;
    double ker[7] ;
    int     i ;
    
    memcpy (ker, CONST_FIR_7, 7 * sizeof (double)) ;
    for (i = 0 ; i < 7 ; i++)
    {
        ker[i] /= dDx ;
    }

    switch (iWhich)
    {
        case 1 :
        {
            cuError = cudaMemcpyToSymbol (constdKer1, ker, sizeof(ker)) ;
        }
        break ;
        case 2 :
        {
            cuError = cudaMemcpyToSymbol (constdKer2, ker, sizeof(ker)) ;
        }
        break ;
        case 3 :
        {
            cuError = cudaMemcpyToSymbol (constdKer3, ker, sizeof(ker)) ;
        }
        break ;
    }
return cuError ;
}

cudaError_t copyFilter14(double dDx,
                         int    iWhich)
{
    cudaError_t cuError = cudaSuccess ;
    double ker[14] ;
    int     i ;
    
    memcpy (ker, CONST_FIR_14, 14 * sizeof (double)) ;
    for (i = 0 ; i < 14 ; i++)
    {               
        ker[i] /= dDx ; 
    }               

    switch (iWhich)
    {
        case 1 :
        {
            cuError = cudaMemcpyToSymbol (constdKer1, ker, sizeof(ker)) ;
        }
        break ;
        case 2 :
        {
            cuError = cudaMemcpyToSymbol (constdKer2, ker, sizeof(ker)) ;
        }
        break ;
        case 3 :
        {
            cuError = cudaMemcpyToSymbol (constdKer3, ker, sizeof(ker)) ;
        }
        break ;
    }
    if (cudaSuccess != cuError)
    {
        printf ("copyFilter14 : %s",cudaGetErrorString(cuError)) ;
    }
    return cuError ;
}

int cuOptimalFilter (int       *iLen1,
                     int       *iLen2,
                     int       *iLen3,
                     int       iSx1,
                     int       iSx2,
                     int       iSx3, 
                     double    dDx1,
                     double    dDx2,
                     double    dDx3)
{
    cudaError_t  cuError = cudaSuccess ;
    int          iRet = 0 ;

    if ((iSx1 > 1) && (iSx1 < 5))
    {
        *iLen1 = 1 ;
        cuError = copyFilter1 (dDx1, 1) ;
    }
    else if ((iSx1 > 4) && (iSx1 < 15))
    {
         *iLen1 = 7 ;
         cuError = copyFilter7 (dDx1, 1) ;
    }
    else if (iSx1 > 15)
    {
         *iLen1 = FILTER_SIZE ;
         switch (*iLen1)
         {
             case 1 : 
                 cuError = copyFilter1 (dDx1, 1) ;
             break ;
             case 7 :
                 cuError = copyFilter7 (dDx1, 1) ;
             break ;
             case 14 : 
                 cuError = copyFilter14 (dDx1, 1) ;
             break ;
         }
    }

    if ((iSx2 > 1) && (iSx2 < 5))
    {
        *iLen2 = 1 ;
        cuError = copyFilter1 (dDx2, 2) ;
    }
    else if ((iSx2 > 4) && (iSx2 < 15))
    {
        *iLen2 = 7 ;
        cuError = copyFilter7 (dDx2, 2) ;
    }
    else if (iSx2 > 15)
    {
         *iLen2 = FILTER_SIZE ;
         switch (*iLen2)
         {
             case 1 :
                 cuError = copyFilter1 (dDx2, 2) ;
             break ;
             case 7 :
                 cuError = copyFilter7 (dDx2, 2) ;
             break ;
             case 14 :
                 cuError = copyFilter14 (dDx2, 2) ;
             break ;
         }      
    }

    if ((iSx3 > 1) && (iSx3 < 5))
    {
        *iLen3 = 1 ;   
        cuError = copyFilter1 (dDx3, 3) ;
    }
    else if ((iSx3 > 4) && (iSx3 < 15))
    {
         *iLen3 = 7 ;
         cuError = copyFilter7 (dDx3, 3) ;
    }
    else if (iSx3 > 15)
    {
         *iLen3 = FILTER_SIZE ;
         switch (*iLen3)
         {
             case 1 :
                 cuError = copyFilter1 (dDx3, 3) ;
             break ;
             case 7 :
                 cuError = copyFilter7 (dDx3, 3) ;
             break ;
             case 14 :
                 cuError = copyFilter14 (dDx3, 3) ;
             break ;
         }
    }

    if (cuError != cudaSuccess)
    {
        printf ("cuOptimalFilter: Failed in memcpy 1\n");
        iRet = 1 ;
    }

#ifdef PRINT_DEBUG_INFO
    printf ("cuOptimalFilter: exited with no errors\n") ;
#endif

    return iRet ;
}


void cuFreeCudaMemory()
{
    CUDA_FREE_MEM(pstStruct) ;
    CUDA_FREE_MEM(pstSig) ;
    CUDA_FREE_MEM(pstMoment) ;
    CUDA_FREE_MEM(pstTau) ;
    CUDA_FREE_MEM(pstEpsilonik) ;
    CUDA_FREE_MEM(pstEpsilonikdot) ;
    CUDA_FREE_MEM(pstPrestress) ;
    CUDA_FREE_MEM (gpV1) ;
    CUDA_FREE_MEM (gpV2) ;
    CUDA_FREE_MEM (gpV3) ;
    CUDA_FREE_MEM (gpU1) ;
    CUDA_FREE_MEM (gpU2) ;
    CUDA_FREE_MEM (gpU3) ;
    CUDA_FREE_MEM (pfDevTract1) ;
    CUDA_FREE_MEM (pfDevTract2) ;
    CUDA_FREE_MEM (pfDevTract3) ;
    CUDA_FREE_MEM (gpGammadot0) ;
}

int checkMemRequirement(int iSx1,
                        int iSx2,
                        int iSx3)
{
    int         liReq = 0 ;
    long int    iTemp = 0 ;
    size_t      iTotalMem = 0 ;
    size_t      iFreeMem = 0 ;

    /* Ui's, Vi's and fft's */
    iTemp=((iSx1+2)*iSx2*iSx3*sizeof(float)*8)/(1024*1024) ;
    liReq+=iTemp ;

    /* sig, moment and tau */
    iTemp=((iSx1*iSx2*iSx3/2)*sizeof(ST_TENSOR)*3)/(1024*1024) ;
    liReq+=iTemp ;  

    iTemp=((iSx1*iSx2*iSx3/2)*sizeof(ST_TENSOR)*2)/(1024*1024) ;
    liReq+=iTemp ;  
    
    /* Ti's */
    iTemp=((iSx1+2)*iSx2*sizeof(float)*3)/(1024*1024) ;
    liReq+=iTemp ;
    
    iTemp=(iSx3/2)*sizeof(ST_LAYER)/(1024*1024) ;
    liReq+=iTemp ;

    /* dMinArray */ 
    iTemp=((iSx1+2)*iSx2*iSx3*sizeof(float))/(1024*1024) ;

    cudaMemGetInfo(&iFreeMem, &iTotalMem);
    iTotalMem/=(1024*1024) ;    
    
    if ((liReq+iTemp) > iTotalMem)
    {
        printf ("\nTotal memory required is : %d MB\n", (int)(liReq+iTemp)) ;
        printf ("Total available is is : %lu MB \n", iTotalMem) ;
        return -1 ;
    }
    
    return 0;    
}


/* ------------------------------------------- utility end -------------------------------------- */


/* ------------------------------------------- Kernels ------------------------------------------ */

__global__ void cuStressUpdateKernel (int        iLen1,
                                      int        iLen2,
                                      int        iLen3,
                                      double     dLambda,
                                      double     dMu,
                                      int        iSx1,
                                      int        iSx2,
                                      int        iSx3,
                                      float      *gpV1,
                                      float      *gpV2,
                                      float      *gpV3,
                                      ST_TENSOR  *pstSig)
{
    unsigned int iInd1 ;
    unsigned int iInd2 ;
    unsigned int iInd3 ;
    unsigned int iIdx ;
    int  iIndex = 0 ;
    int  iInd1m = 0 ;
    int  iInd2m = 0 ;
    int  iInd3m = 0 ;
    int  iInd1p = 0 ;
    int  iInd2p = 0 ;
    int  iInd3p = 0 ;
    int  iTemp = 0 ;
    
    int  iOffp = 0 ;
    int  iOffm = 0 ;
    ST_TENSOR  stT = {0};

#ifdef STRESS_SHARED_MEM
    __shared__ float shdV11 [256] ;
    __shared__ float shdV21 [256] ;
    __shared__ float shdV31 [256] ;
#endif

    iInd1 = threadIdx.x ;
    iInd2 = blockIdx.y ;
    iInd3 = blockIdx.x ;

#ifdef STRESS_SHARED_MEM
    iIdx = (((iInd3 * iSx2) + iInd2)  * iSx1) + iInd1 ;

    if (0 == iInd1)
    {
        for (iIndex = 0; iIndex < iSx1 ; iIndex++)
        {
            shdV11[iIndex] = gpV1[iIdx] ;
            shdV21[iIndex] = gpV2[iIdx] ;
            shdV31[iIndex] = gpV3[iIdx] ;
            iIdx++ ;    
        }
    
    }
    __syncthreads() ;
#endif



    if ((iInd1 < iSx1) && (iInd2 < iSx2) && ((iInd3 > (iLen3-1)) && (iInd3 < (iSx3-iLen3))))
    {
        iIdx = (((iInd3 * iSx2) + iInd2)  * iSx1) + iInd1 ;

#ifdef STRESS_SHARED_MEM
        iTemp = iSx1 - 1 ;
        for (iIndex = 0; iIndex < iLen1; iIndex++)
        {
            iInd1m = ((iTemp+iInd1-iIndex) % (iSx1)) ;
            iInd1p = ((iInd1+iIndex) % (iSx1)) + 1 ;
            iInd1p = (iInd1p % iSx1) ;
            
            stT.s11 += ((shdV11[iInd1p] - shdV11[iInd1m]) * constdKer1[iIndex]) ;
            stT.s12 += ((shdV21[iInd1p] - shdV21[iInd1m]) * constdKer1[iIndex]) ;
            stT.s13 += ((shdV31[iInd1p] - shdV31[iInd1m]) * constdKer1[iIndex]) ;
        }
#else
        iTemp = iSx1 - 1 ;
        for (iIndex = 0; iIndex < iLen1; iIndex++)
        {
            iInd1m = ((iTemp+iInd1-iIndex) % (iSx1)) ;
            iInd1p = ((iInd1+iIndex) % (iSx1)) + 1 ;
            iInd1p = (iInd1p % iSx1) ;

            iOffp = (((iInd3 * iSx2) + iInd2) * (iSx1 + 2)) + iInd1p ;
            iOffm = (((iInd3 * iSx2) + iInd2) * (iSx1 + 2)) + iInd1m ;

            stT.s11 += ((gpV1[iOffp] - gpV1[iOffm]) * constdKer1[iIndex]) ;
            stT.s12 += ((gpV2[iOffp] - gpV2[iOffm]) * constdKer1[iIndex]) ;
            stT.s13 += ((gpV3[iOffp] - gpV3[iOffm]) * constdKer1[iIndex]) ;
        }
#endif
        iTemp = iSx2 - 1 ;

        for (iIndex = 0; iIndex < iLen2; iIndex++)
        {
            iInd2m = ((iTemp+iInd2-iIndex) % (iSx2)) ;
            iInd2p = ((iInd2+iIndex) % (iSx2)) + 1 ;
            iInd2p = (iInd2p % iSx2) ;

            iOffp = (((iInd3 * iSx2) + iInd2p) * (iSx1 + 2)) + iInd1 ;
            iOffm = (((iInd3 * iSx2) + iInd2m) * (iSx1 + 2)) + iInd1 ;

            stT.s12 += ((gpV1[iOffp] - gpV1[iOffm]) * constdKer2[iIndex]) ;
            stT.s22 += ((gpV2[iOffp] - gpV2[iOffm]) * constdKer2[iIndex]) ;
            stT.s23 += ((gpV3[iOffp] - gpV3[iOffm]) * constdKer2[iIndex]) ;
        }

        for (iIndex = 1; iIndex <= iLen3; iIndex++)
        {
            iInd3m = iInd3 - iIndex ;
            iInd3p = iInd3 + iIndex ;

            iOffp = (((iInd3p * iSx2) + iInd2) * (iSx1 + 2)) + iInd1 ;
            iOffm = (((iInd3m * iSx2) + iInd2) * (iSx1 + 2)) + iInd1 ;

            stT.s13 += ((gpV1[iOffp] - gpV1[iOffm]) * constdKer3[iIndex-1]) ;
            stT.s23 += ((gpV2[iOffp] - gpV2[iOffm]) * constdKer3[iIndex-1]) ;
            stT.s33 += ((gpV3[iOffp] - gpV3[iOffm]) * constdKer3[iIndex-1]) ;
        }

        stT.s12 /= 2.0 ;
        stT.s13 /= 2.0 ;
        stT.s23 /= 2.0 ;

        cuIsotrpicStressStrain (&stT, dLambda, dMu) ;
        cuTensorOperate (&pstSig[iIdx], (void *)&stT, '+') ;

    }

}


__global__ void cuLocalStressStrainKernel (int        iInd3m,
                                           int        iInd3p, 
                                           int        iLen1, 
                                           int        iLen2, 
                                           int        iInd3, 
                                           double     dPx3, 
                                           double     dLambda,
                                           double     dMu, 
                                           int        iSx1, 
                                           int        iSx2, 
                                           int        iSx3,
                                           float      *gpV1, 
                                           float      *gpV2, 
                                           float      *gpV3,
                                           ST_TENSOR  *pstSig) 
{
    ST_TENSOR stT ;
    int iInd1 = 0 ;
    int iInd2 = 0 ;
    int iIdx = 0 ;

    iInd2 = blockIdx.x ;
    iInd1 = threadIdx.x ;

    cuTensorMemset(&stT) ;
            
    iIdx = (((iInd3 * iSx2) + iInd2)  * (iSx1)) + iInd1 ;                        
    
    cuLocalStrain_ani (&stT, iInd3m, iInd3p, dPx3, iLen1, iLen2,
               iInd1, iInd2, iInd3, gpV1, gpV2, gpV3, iSx1, iSx2, iSx3) ;
        
            
    cuIsotrpicStressStrain (&stT, dLambda, dMu) ;

    cuTensorOperate (&pstSig[iIdx], (void *)&stT, '+') ;

}


__global__ void cuEquiBodyKernel (ST_TENSOR   *pstT,
                                  int         iLen1,
                                  int         iLen2,
                                  int         iLen3,
                                  int         iSx1,
                                  int         iSx2,
                                  int         iSx3,
                                  float       *fData1,
                                  float       *fData2,
                                  float       *fData3)
{
    double  f1 ;
    double  f2 ;
    double  f3 ;
    int     iInd1 = 0 ;
    int     iInd2 = 0 ;
    int     iInd3 = 0 ;
    int     iIdx = 0 ;
    int     iIndex = 0 ;
    int     iIndm = 0 ;
    int     iIndp = 0 ;
    int     iTemp = 0 ;
    int     iOffp = 0 ;
    int     iOffm = 0 ;

#ifdef ENABLE_REG_BLOCKING
    float   fTen_1_s1[14] ;
    float   fTen_1_s2[14] ;
    float   fTen_1_s3[14] ;
    float   fTen_2_s1[14] ;
    float   fTen_2_s2[14] ;
    float   fTen_2_s3[14] ;
#endif

    ST_TENSOR   stTemp1 ;
    ST_TENSOR   stTemp2 ;

    iInd1 = threadIdx.x ;
    iInd2 = blockIdx.y ;
    iInd3 = blockIdx.x  + iLen3 ;

    if ((iInd1 < iSx1) && (iInd2 < iSx2) && ((iInd3 > (iLen3-1)) && (iInd3 < (iSx3-iLen3))))
    {
#ifdef ENABLE_REG_BLOCKING
        iTemp = iSx1 - 1 ;
        for (iIndex = 0; iIndex < iLen1; iIndex++)
        {
            iIndm = ((iTemp+iInd1-iIndex) % (iSx1)) ;
            iIndp = (((iInd1+iIndex) % (iSx1)) + 1) % iSx1 ;

            iOffm = (((iInd3 * iSx2) + iInd2) * (iSx1)) + iIndm ;
            iOffp = (((iInd3 * iSx2) + iInd2) * (iSx1)) + iIndp ;

            stTemp1 = pstT[iOffp] ;
            stTemp2 = pstT[iOffm] ;
            fTen_1_s1[iIndex] = stTemp1.s11 ;
            fTen_1_s2[iIndex] = stTemp1.s12 ;
            fTen_1_s3[iIndex] = stTemp1.s13 ;

            fTen_2_s1[iIndex] = stTemp2.s11 ;
            fTen_2_s2[iIndex] = stTemp2.s12 ;
            fTen_2_s3[iIndex] = stTemp2.s13 ;
        }

#pragma unroll 
        for (iIndex = 0; iIndex < iLen1; iIndex++)
        {
            f1 += ((fTen_1_s1[iIndex] - fTen_2_s1[iIndex]) * constdKer1[iIndex]) ;
            f2 += ((fTen_1_s2[iIndex] - fTen_2_s2[iIndex]) * constdKer1[iIndex]) ;
            f3 += ((fTen_1_s3[iIndex] - fTen_2_s3[iIndex]) * constdKer1[iIndex]) ;
        }

#else
        iTemp = iSx1 - 1 ;
        for (iIndex = 0; iIndex < iLen1; iIndex++)
        {
            iIndm = ((iTemp+iInd1-iIndex) % (iSx1)) ;
            iIndp = (((iInd1+iIndex) % (iSx1)) + 1) % iSx1 ;

            iOffm = (((iInd3 * iSx2) + iInd2) * (iSx1)) + iIndm ;
            iOffp = (((iInd3 * iSx2) + iInd2) * (iSx1)) + iIndp ;

            stTemp1 = pstT[iOffp] ;
            stTemp2 = pstT[iOffm] ;

            f1 += ((stTemp1.s11 - stTemp2.s11) * constdKer1[iIndex]) ;
            f2 += ((stTemp1.s12 - stTemp2.s12) * constdKer1[iIndex]) ;
            f3 += ((stTemp1.s13 - stTemp2.s13) * constdKer1[iIndex]) ;
        }
        iTemp = iSx2 - 1 ;

        for (iIndex = 0; iIndex < iLen2; iIndex++)
        {
            iIndm = ((iTemp+iInd2-iIndex) % (iSx2)) ;
            iIndp = (((iInd2+iIndex) % (iSx2)) + 1) % iSx2 ;

            iOffm = (((iInd3 * iSx2) + iIndm) * iSx1) + iInd1 ;
            iOffp = (((iInd3 * iSx2) + iIndp) * iSx1) + iInd1 ;

            stTemp1 = pstT[iOffp] ;
            stTemp2 = pstT[iOffm] ;

            f1 += ((stTemp1.s12 - stTemp2.s12) * constdKer2[iIndex]) ;
            f2 += ((stTemp1.s22 - stTemp2.s22) * constdKer2[iIndex]) ;
            f3 += ((stTemp1.s23 - stTemp2.s23) * constdKer2[iIndex]) ;
        }

        for (iIndex = 1; iIndex <= iLen3; iIndex++)
        {
            iIndm = iInd3 - iIndex ;
            iIndp = iInd3 + iIndex ;

            iOffm = (((iIndm * iSx2) + iInd2) * iSx1) + iInd1 ;
            iOffp = (((iIndp * iSx2) + iInd2) * iSx1) + iInd1 ;

            stTemp1 = pstT[iOffp] ;
            stTemp2 = pstT[iOffm] ;

            f1 += ((stTemp1.s13 - stTemp2.s13) * constdKer3[iIndex-1]) ;
            f2 += ((stTemp1.s23 - stTemp2.s23) * constdKer3[iIndex-1]) ;
            f3 += ((stTemp1.s33 - stTemp2.s33) * constdKer3[iIndex-1]) ;
        }
#endif

        iIdx = (((iInd3 * iSx2) + iInd2) * (iSx1 + 2)) + iInd1 ;

        fData1[iIdx] -= f1 ;
        fData2[iIdx] -= f2 ;
        fData3[iIdx] -= f3 ;
    }
}


__global__ void cuLocalDivergenceKernel (int          iInd3m,
                                         int          iInd3p,
                                         int          iLen1,
                                         int          iLen2,
                                         int          iInd3,
                                         double       dPx3,
                                         int          iSx1,
                                         int          iSx2,
                                         int          iSx3,
                                         float        *gpV1,
                                         float        *gpV2,
                                         float        *gpV3,
                                         ST_TENSOR    *pstSig)
{
    int     iInd1 = 0 ;
    int     iInd2 = 0 ;
    int     iIdx = 0 ;
    double  f1 = 0.0 ;
    double  f2 = 0.0 ;
    double  f3 = 0.0 ;

    iInd2 = blockIdx.x ;
    iInd1 = threadIdx.x ;

    if ((iInd1 < iSx1) && (iInd2 < iSx2))
    {
        cuLocalDivergence_ani (pstSig, iInd3m, iInd3p, dPx3, iLen1, iLen2, iInd1, iInd2, iInd3,
                               &f1, &f2, &f3, iSx1, iSx2, iSx3) ;

        iIdx = (((iInd3 * iSx2) + iInd2) * (iSx1 + 2)) + iInd1 ;

        gpV1[iIdx] -= f1 ;
        gpV2[iIdx] -= f2 ;
        gpV3[iIdx] -= f3 ;
    }
}

__global__ void cuEquivalentTraction (float      *pfTraction1,
                                      float      *pfTraction2,
                                      float      *pfTraction3,
                                      ST_TENSOR  *pstSig,
                                      int        iSx1,
                                      int        iSx2)
{
    unsigned int iIdx1 = 0 ;
    unsigned int iIdx2 = 0 ;
    unsigned int iIdx = 0 ;

    iIdx1 = threadIdx.x ;
    iIdx2 = blockIdx.x ;

    if ((iIdx1 < iSx1) && (iIdx2 < iSx2))
    {
        iIdx = (iIdx2 * (iSx1 + 2)) + iIdx1 ;

        pfTraction1[iIdx] += pstSig[iIdx].s13 ;
        pfTraction2[iIdx] += pstSig[iIdx].s23 ;
        pfTraction3[iIdx] += pstSig[iIdx].s33 ;
    }
}

__global__ void cuSourceForceKernel (int     iSx1,
                                     int     iSx2,
                                     int     iSx3,
                                     double  dDx1,
                                     double  dDx2,
                                     double  dDx3,
                                     double  dcStrike,
                                     double  dsStrike,
                                     double  dcDip,
                                     double  dsDip,
                                     double  dCr,
                                     double  dSr,
                                     double  dScale,
                                     double  dWp,
                                     double  dLp,
                                     double  dX2r,
                                     double  dXr,
                                     double  dYr,
                                     double  dZr,
                                     double  dX,
                                     double  dY,
                                     double  dZ,
                                     double  dDip,
                                     double  dW,
                                     double  dL,
                                     double  dBeta,
                                     double  dMu,
                                     double  dS,
                                     float   *pfData1,
                                     float   *pfData2,
                                     float   *pfData3)
{
    int     iInd1 = 0 ;
    int     iInd2 = 0 ;
    int     iInd3 = 0 ;
    int     iIdx = 0 ;
    int     iTemp = 0 ;

    double  dX1 ;
    double  dX2 ;
    double  dX3 ;

    double  dX1s ;
    double  dX2s ;
    double  dX3s ;

    double  dX1i ;
    double  dX3i ;

    double  dSource ;
    double  dTemp1 ;
    double  dTemp2 ;
    double  dTemp3 ;

    double  dDblcp ;
    double  dCplei ;
//      double          dDipci ;
    double  dImage ;
    double  dDipcs ;

    iInd1 = threadIdx.x ;
    iInd2 = blockIdx.y ;
    iInd3 = blockIdx.x ;

    if ((iInd1 < iSx1) && (iInd2 < iSx2) && (iInd3 < iSx3/2))
    {
        cuShiftedCoordinates (1, 1, iInd3, iSx1, iSx2, iSx3,
                              dDx1, dDx2, dDx3, &dX1, &dX2, &dX3) ;
        iTemp = ((abs(dX3 - dZ) > dLp) && (abs(dX3 + dZ) > dLp)) ;
        if (0 == iTemp)
        {
            cuShiftedCoordinates (iInd1, iInd2, iInd3, iSx1, iSx2, iSx3,
                                  dDx1, dDx2, dDx3, &dX1, &dX2, &dX3) ;

            iTemp = ((abs(dX1 - dX) > MAX_NUM(dWp, dLp)) || (abs(dX2 - dY) > MAX_NUM(dWp, dLp))) ;
            if (0 == iTemp)
            {
                dX2r = (dcStrike * dX1) - (dsStrike * dX2) ;
                dX1s = (dcDip * dX2r) - (dsDip * dX3) ;
                dX1i = (dcDip * dX2r) + (dsDip * dX3) ;

                iTemp = ((abs(dX1s - dXr) > (7.01 * dDx1)) && (abs(dX1i - dXr) > (7.01 * dDx1))) ;
                if (0 == iTemp)
                {
                    dX2s = (dsStrike * dX1) + (dcStrike * dX2) ;
                    dX3s = (dsDip * dX2r) + (dcDip * dX3) ;
                    dX3i = (-dsDip * dX2r) + (dcDip * dX3) ;

                    dTemp1 = cuGauss (dX1s - dXr, dDx1) ;
                    dTemp2 = cuOmega ((dX2s - dYr)/dW, dBeta) ;
                    dTemp3 = cuOmega ((dX3s - dZr)/dL, dBeta) ;

                    dSource = dScale * cuGaussp (dX1s-dXr, dDx1) * dTemp2 * dTemp3 ;
                    dDblcp = (dScale / dW) * dTemp1 * cuOmegap ((dX2s-dYr)/dW, dBeta) * dTemp3 ;
                    dDipcs = (dScale / dL) * dTemp1 * dTemp2 * cuOmegap ((dX3s-dZr)/dL, dBeta) ;

                    dTemp1 = cuGauss (dX1i - dXr, dDx1) ;
                    dTemp3 = cuOmega ((dX3i + dZr)/dL, dBeta) ;
                    dImage = dScale * cuGaussp (dX1i-dXr, dDx1) * dTemp2 * dTemp3 ;
                    dCplei = (dScale / dW) * dTemp1 * cuOmegap ((dX2s-dYr)/dW, dBeta) * dTemp3 ;
//                                    dDipci = (dScale / dL) * dTemp1 * dTemp2 * cuOmegap ((dX3i+dZr)/dL, dBeta) ;

                    iIdx = ((iInd3 * iSx2) + iInd2) * (iSx1 + 2) + iInd1 ;
                    if ((2.01 * DEG2RAD) > dDip)
                    {
                        pfData1[iIdx] += ((dCr * dsStrike * (dSource + dImage)) +
                                          (dCr * dcDip * dcStrike * (dDblcp + dCplei))) ;
                        pfData2[iIdx] += ((dCr * dcStrike * (dSource + dImage)) -
                                          (dCr * dcDip * dsStrike * (dDblcp + dCplei))) ;
                        pfData3[iIdx] -= (dCr * dsDip * dDblcp) ;
                    }
                    else
                    {
                        pfData1[iIdx] += ((dCr * dsStrike * dSource) +
                                          (dCr * dcDip * dcStrike * dDblcp)) ;
                        pfData2[iIdx] += ((dCr * dcStrike * (dSource)) -
                                          (dCr * dcDip * dsStrike * (dDblcp))) ;
                        pfData3[iIdx] -= (dCr * dsDip * dDblcp) ;
                    }

                    pfData1[iIdx] += ((dcDip * dSr * dcStrike * dDipcs) +
                                      (dsDip * dSr * dcStrike * dSource)) ;
                    pfData2[iIdx] -= ((dcDip * dSr * dsStrike * dDipcs) +
                                      (dsDip * dSr * dsStrike * dSource)) ;
                    pfData3[iIdx] += ((dcDip * dSr * dSource) - (dsDip * dSr * dDipcs)) ;
                }
            }
        }
    }
}

__global__ void cuSourceTractionKernel (int     iSx1,
                                        int     iSx2,
                                        int     iSx3,
                                        double  dDx1,
                                        double  dDx2,
                                        double  dDx3,
                                        double  dcStrike,
                                        double  dsStrike,
                                        double  dcDip,
                                        double  dsDip,
                                        double  dCr,
                                        double  dSr,
                                        double  dScale,
                                        double  dWp,
                                        double  dLp,
                                        double  dX2r,
                                        double  dXr,
                                        double  dYr,
                                        double  dZr,
                                        double  dX,
                                        double  dY,
                                        double  dW,
                                        double  dL,
                                        double  dBeta,
                                        double  dMu,
                                        double  dS,
                                        float   *pfTract1,
                                        float   *pfTract2,
                                        float   *pfTract3)
{
    int     iInd1 = 0 ;
    int     iInd2 = 0 ;
    int     iInd3 = 0 ;
    int     iIdx = 0 ;
    int     iTemp = 0 ;

    double  dX1 ;
    double  dX2 ;
    double  dX3 ;

    double  dX1s ;
    double  dX2s ;
    double  dX3s ;

    double  dX1i ;
    double  dX3i ;

    double  dSource = 0.0 ;
    double  dTemp1 ;
    double  dTemp2 ;
    double  dTemp3 ;

    double  dN[3] = {0.0} ;
    double  dB[3] = {0.0} ;
    double  dBmod[3] = {0.0} ;


    ST_TENSOR  dM ;

    iInd1 = threadIdx.x ;
    iInd2 = blockIdx.x ;
    iInd3 = 0 ;

    if((iInd1 < iSx1) && (iInd2 < iSx2))
    {
        cuShiftedCoordinates (iInd1, iInd2, iInd3, iSx1, iSx2, iSx3,
                              dDx1, dDx2, dDx3, &dX1, &dX2, &dX3) ;

        iTemp = ((abs(dX1 - dX) > MAX_NUM(dWp, dLp)) || (abs(dX2 - dY) > MAX_NUM(dWp, dLp))) ;
        if (0 == iTemp)
        {
            dX2r = (dcStrike * dX1) - (dsStrike * dX2) ;
            dX1s = (dcDip * dX2r) - (dsDip * dX3) ;
            dX1i = (dcDip * dX2r) + (dsDip * dX3) ;

            iTemp = ((abs(dX1s - dXr) > (7.01 * dDx1)) && (abs(dX1i - dXr) > (7.01 * dDx1))) ;
            if (0 == iTemp)
            {
                dX2s = (dsStrike * dX1) - (dcStrike * dX2) ;
                dX3s = (dsDip * dX2r) - (dcDip * dX3) ;
                dX3i = (-dsDip * dX2r) + (dcDip * dX3) ;

                dTemp1 = cuGauss (dX1s - dXr, dDx1) ;
                dTemp2 = cuOmega ((dX2s - dYr)/dW, dBeta) ;
                dTemp3 = cuOmega ((dX3s - dZr)/dL, dBeta) ;
                dSource += (dTemp1 * dTemp2 * dTemp3) ;

                dTemp1 = cuGauss (dX1i - dXr, dDx1) ;
                dTemp3 = cuOmega ((dX3i + dZr)/dL, dBeta) ;
                dSource += (dTemp1 * dTemp2 * dTemp3) ;

                dN[0] = dcDip * dcStrike * dSource ;
                dN[1] = -dcDip * dsStrike * dSource ;
                dN[2] = -dsDip * dSource ;


                dB[0] = dsStrike * dCr ;
                dB[1] = dcStrike * dCr ;

                dB[0] += (dcStrike * dsDip * dSr) ;
                dB[1] -= (dsStrike * dsDip * dSr) ;
                dB[2] = dcDip * dSr ;

                dBmod[0] = dMu * dS * dB[0] ;
                dBmod[1] = dMu * dS * dB[1] ;
                dBmod[2] = dMu * dS * dB[2] ;

                //dyadic product
                cuTensorDyadProd (&dM, dN, dBmod) ;

                iIdx = (iInd2 * (iSx1 + 2)) + iInd1 ;

                pfTract1[iIdx] += (float)dM.s13 ;
                pfTract2[iIdx] += (float)dM.s23 ;
                pfTract3[iIdx] += (float)dM.s33 ;
            }
        }
    }
}

__global__ void cuTransientEigenKernel (ST_LAYER          *pstStruct,
                                        ST_TENSOR         *pstSig,
                                        ST_TENSOR         *pstMoment,
                                        float             *gpGammadot0,
                                        ST_TENSOR         *pstEpsilonik,
                                        ST_TENSOR         *pstEpsilonikdot,
                                        double            dMu,
                                        int               iSx1,
                                        int               iSx2,
                                        int               iSx3,
                                        double            dDx1,
                                        double            dDx2,
                                        double            dDx3,
                                        int               bMaxwell,
                                        float             *dMinArray,
                                        bool              isdgammadot0)
{
    int        iInd1 = 0 ;
    int        iInd2 = 0 ;
    int        iInd3 = 0 ;
    int        iIdx = 0 ;
    double     dPower = 0.0 ;
    double     dMuk = 0.0 ;
    double     dGammaDot0 = 0.0 ;
    double     dGammaDot = 0.0 ;
    double     dNq = 0.0 ;
    ST_TENSOR  stS = {0} ;
    ST_TENSOR  stQ = {0} ;
    ST_TENSOR  stTemp = {0} ;
    double     dTemp = 0.0 ;
    int        iCond = 0 ;

    iInd3 = blockIdx.x ;
    iInd2 = blockIdx.y ;
    iInd1 = threadIdx.x ;

    if ((iInd1 < iSx1) && (iInd2 < iSx2) && (iInd3 < iSx3))
    {
        iIdx = (((iInd3 * iSx2) + iInd2) * iSx1) + iInd1 ;
        if (bMaxwell)
        {
            dMinArray[iIdx] = 1.0e+30 ;
        }
        dPower = pstStruct[iInd3].stressexponent - 1 ;
        dMuk = pstStruct[iInd3].Gk;

        dGammaDot0 = pstStruct[iInd3].gammadot0 ;

        if (isdgammadot0) 
        {
            dGammaDot0 += gpGammadot0[iIdx] ;
        }

        iCond = (1.0e-20 > dGammaDot0) ? 1 : 0 ;
        if (0 == iCond)
        {
            cuTensorDeviatoric (&pstSig[iIdx], &stS) ;

            dTemp = 0.0 ;            
            cuTensorOperate (&stTemp, (void *)&dTemp, '*') ;

            dTemp = 2 * dMuk ;
            cuTensorOperate (&stTemp, (void *)(&pstEpsilonik[iIdx]), '=') ;

            // 2Gk*epsilonik
            cuTensorOperate (&stTemp, (void *)&dTemp, '*') ;

            cuTensorOperate (&stQ, (void *)&stS, '=') ;
            //  Q = (sigma - 2Gk*epsilonik)
            cuTensorOperate (&stQ, (void *)&stTemp, '-') ;

            // q = || Q ||
            dNq = cuTensorNorm(&stQ) ;
             
            // powerlaw viscosity
            dTemp = dNq/dMu ;
            dGammaDot = dGammaDot0 * pow (dTemp, dPower) ;
            
            cuTensorOperate (&pstEpsilonikdot[iIdx], (void *)&stQ, '=') ;

            dTemp = dGammaDot / dMu ;
            cuTensorOperate (&pstEpsilonikdot[iIdx], (void *)&dTemp, '*') ;

            dTemp = 2 * dGammaDot ;
            cuTensorOperate (&stQ, (void *)&dTemp, '*') ;

            // update moment density forcing
            cuTensorOperate (&pstMoment[iIdx], (void *)&stQ, '+') ;

            if (bMaxwell)
            {
                dTemp = 1 / dGammaDot ;
                if (0 != dTemp)
                {
                        dMinArray[iIdx] = (float)dTemp ;
                }
            }
        }
        else
        {
            dTemp=0.0;
            cuTensorOperate (&pstEpsilonikdot[iIdx], (void *)&dTemp, '*') ;
        }
    }
}

__global__ void cuViscousEigenKernel (ST_LAYER          *pstStruct,
                                      ST_TENSOR         *pstSig,
                                      ST_TENSOR         *pstMoment,
                                      ST_TENSOR_LAYER   *pstPrestress,
                                      float             *gpGammadot0,
                                      double            dMu,
                                      int               iSx1,
                                      int               iSx2,
                                      int               iSx3,
                                      double            dDx1,
                                      double            dDx2,
                                      double            dDx3,
                                      float             *dMinArray,
                                      float             *pGamma,
                                      bool              bPresent,
                                      bool              bGammaPresent,
                                      bool              isdgammadot0)
{
    int        iInd1 = 0 ;
    int        iInd2 = 0 ;
    int        iInd3 = 0 ;
    int        iIdx = 0 ;
    double     dPower = 0.0 ;
    double     dCohesion = 0.0 ;
    double     dGammaDot0 = 0.0 ;
    double     dGammaDot = 0.0 ;
    double     dGammaDotp = 0.0 ;
    double     dTau = 0.0 ;
    double     dTauc = 0.0 ;
    double     dTaup = 0.0 ;
    ST_TENSOR  stS = {0} ;
    ST_TENSOR  stP = {0} ;
    ST_TENSOR  stSP = {0} ;
    double     dTemp = 0.0 ;
    int        iCond = 0 ;

    iInd3 = blockIdx.x ;
    iInd2 = blockIdx.y ;
    iInd1 = threadIdx.x ;

    if ((iInd1 < iSx1) && (iInd2 < iSx2) && (iInd3 < iSx3))
    {
        iIdx = (((iInd3 * iSx2) + iInd2) * iSx1) + iInd1 ;
        if (bPresent)
        {
                dMinArray[iIdx] = 1.0e+30 ;
        }

        dPower = pstStruct[iInd3].stressexponent ;
        dCohesion = pstStruct[iInd3].cohesion ;

        cuTensorDeviatoric (&pstPrestress[iInd3].t, &stP) ;
        dTaup = cuTensorNorm(&stP) ;

        dGammaDot0 = pstStruct[iInd3].gammadot0 ;

        if (isdgammadot0) 
        {
            dGammaDot0 += gpGammadot0[iIdx] ;
        }

        iCond = (1.0e-20 > dGammaDot0) ? 1 : 0 ;

        if (0 == iCond)
        {
            cuTensorDeviatoric (&pstSig[iIdx], &stS) ;

            dTemp=0.0;
            cuTensorOperate (&stSP, (void *)&dTemp, '*') ;
            cuTensorOperate (&stSP, (void *)&stS, '+') ;
            cuTensorOperate (&stSP, (void *)&stP, '+') ;

            dTau = cuTensorNorm(&stSP) ;

            dTauc = MAX_NUM(0, dTau - dCohesion) ;

            iCond = (dTauc <= 1.0e-20) ? 1 : 0 ;
            if (0 == iCond)
            {
                dTemp = dTauc/dMu ;
                dGammaDot = dGammaDot0 * pow (dTemp, dPower-1) ;
                dTemp = dTaup/dMu ;
                dGammaDotp = dGammaDot0 * pow (dTaup/dMu, dPower-1) ;

                dTemp = 2 * dGammaDot ;
                cuTensorOperate (&stS, (void *)&dTemp, '*') ;
                
                dTemp = 2 * dGammaDotp ;
                cuTensorOperate (&stP, (void *)&dTemp, '*') ;

                cuTensorOperate (&pstMoment[iIdx], (void *)&stS, '+') ;
                cuTensorOperate (&pstMoment[iIdx], (void *)&stP, '-') ;

                /*if (1 == bGammaPresent)
                {
                        iIdx = (((iInd3 * iSx2) + iInd2) * (iSx1+2)) + iInd1 ;
                        pGamma[iIdx] = (float) dGammaDot ;
                }*/
                if (bPresent)
                {
                    dTemp = 1 / dGammaDot ;
                    iIdx = (((iInd3 * iSx2) + iInd2) * iSx1) + iInd1 ;
                    if (0 != dTemp)
                    {
                            dMinArray[iIdx] = (float)dTemp ;
                    }
                }
            }
        }
    }
}

__global__ void cuTensorFieldKernel (ST_TENSOR  *pstTens1,
                                     ST_TENSOR  *pstTens2,
                                     float      fC1,
                                     float      fC2,
                                     int        iSx1,
                                     int        iSx2,
                                     int        iSx3)
{
    int        iInd1 = 0 ;
    int        iInd2 = 0 ;
    int        iInd3 = 0 ;
    int        iIdx = 0 ;
    ST_TENSOR  stTemp = {0} ;
    double     dC1 = fC1 ;
    double     dC2 = fC2 ;

    iInd3 = blockIdx.x ;
    iInd2 = blockIdx.y ;
    iInd1 = threadIdx.x ;

    if ((iInd1 < iSx1) && (iInd2 < iSx2) && (iInd3 < iSx3))
    {
        iIdx = (((iInd3 * iSx2) + iInd2) * iSx1) + iInd1 ;
        cuTensorOperate (&stTemp, &pstTens2[iIdx], '+') ;
        cuTensorOperate (&pstTens1[iIdx], (void *)&dC1, '*') ;
        cuTensorOperate (&stTemp, (void *)&dC2, '*') ;
        cuTensorOperate (&pstTens1[iIdx], &stTemp, '+') ;
    }
}

__global__ void cuTensorAmpKernel (ST_TENSOR    *pstTensor,
                                   double       *pdArray,
                                   int           iSx1,
                                   int           iSx2,
                                   int           iSx3)
{
    int     iInd1 = 0 ;
    int     iInd2 = 0 ;
    int     iInd3 = 0 ;
    int     iIdx = 0 ;
    double  dTemp = 0.0 ;

    iInd3 = blockIdx.x ;
    iInd2 = blockIdx.y ;
    iInd1 = threadIdx.x ;

    if ((iInd1 < iSx1) && (iInd2 < iSx2) && (iInd3 < iSx3))
    {
        iIdx = (((iInd3 * iSx2) + iInd2) * iSx1) + iInd1 ;
        dTemp = cuTensorNorm(&pstTensor[iIdx]) ;
        pdArray[iIdx] = dTemp ;
    }
}

__global__ void cuBuildGammadotKernel (int           iSx1,
                                       int           iSx2,
                                       int           iSx3,
                                       double        dDx1,
                                       double        dDx2,
                                       double        dDx3,
                                       double        dBeta, 
                                       ST_WEAK       *pstDuctile,
                                       int           iNz,
                                       float         *pfDgammadot0)
{
    int     iInd1 = 0 ;
    int     iInd2 = 0 ;
    int     iInd3 = 0 ;
    int     iIdx = 0 ;
    
    double     dX1 ;
    double     dX2 ;
    double     dX3 ;
    double     dDum ;


    iInd3 = blockIdx.x ;
    iInd2 = blockIdx.y ;
    iInd1 = threadIdx.x ;

    if ((iInd1 < iSx1) && (iInd2 < iSx2) && (iInd3 < iSx3))
    {
        dX3 = (iInd3) * dDx3 ;
        cuShiftedCoordinates (iInd1, iInd2, iInd3, iSx1, iSx2, iSx3,
                              dDx1, dDx2, dDx3, &dX1, &dX2, &dDum) ;
    
        iIdx = (((iInd3 * iSx2) + iInd2) * iSx1) + iInd1 ;
        pfDgammadot0[iIdx] = cuDgGammaDotNot(pstDuctile, iNz, dX1, dX2, dX3, dBeta) ; 
    }
}

__global__ void cuFieldAddKernel (float  *pfData1,
                                  float  *pfData2,
                                  float  fC1,
                                  float  fC2,
                                  int    iSx1,
                                  int    iSx2,
                                  int    iSx3)
{
    int     iInd1 = 0 ;
    int     iInd2 = 0 ;
    int     iInd3 = 0 ;
    int     iIdx = 0 ;
    float   fTemp = 0 ;
    float2  *fData1 ;
    float2  *fData2 ;

    iInd3 = blockIdx.x ;
    iInd2 = blockIdx.y ;
    iInd1 = threadIdx.x ;

    if ((iInd1 < iSx1) && (iInd2 < iSx2) && (iInd3 < iSx3))
    {
        iIdx = (((iInd3 * iSx2) + iInd2) * iSx1) + 2 * iInd1 ;
        fData1 = (float2 *) &pfData1[iIdx] ;
        fData2 = (float2 *) &pfData2[iIdx] ;

        fTemp = (float) (fC2 * (*fData2).x) ;
        (*fData1).x = (float) (fC1 * (*fData1).x) ;
        (*fData1).x += fTemp ;

        fTemp = (float) (fC2 * (*fData2).y) ;
        (*fData1).y = (float) (fC1 * (*fData1).y) ;
        (*fData1).y += fTemp ;
    }
}

__global__ void cuFrictionStress (double     dScaling,
                                  double     dcStrike,
                                  double     dsStrike,
                                  double     dcDip,
                                  double     dsDip,
                                  double     dCr,
                                  double     dSr,
                                  double     dWp,
                                  double     dLp,
                                  double     dX2r,
                                  double     dXr,
                                  double     dYr,
                                  double     dZr,
                                  int        iSx1,
                                  int        iSx2,
                                  int        iSx3,
                                  double     dDx1,
                                  double     dDx2,
                                  double     dDx3,
                                  double     dX,
                                  double     dY,
                                  double     dZ,
                                  double     dL,
                                  double     dW,
                                  double     dRake,
                                  double     dMu,
                                  int        bPresent,
                                  double     dBeta,
                                  float      *dMinArray,
                                  ST_TENSOR  *pMoment,
                                  ST_TENSOR  *pSig,
                                  ST_LAYER   *pStruct)
{
    double     dN[3] ;
    double     dnTemp[3] ;
    double     dR[3] ;
    double     dT[3] ;
    double     dTs[3] ;

    int        iTemp ;
    int        iIdx ;
    int        iIdx1 ;
    int        iIdx2 ;
    int        iIdx3 ;

    double     dX1s ;
    double     dX2s ;
    double     dX3s ;

    double     dX1i ;
    double     dX3i ;

    double     dTemp1 ;
    double     dTemp2 ;
    double     dTemp3 ;

    double     dX1 ;
    double     dX2 ;
    double     dDum ;

    // make these shared 
    double     dVo ;
    double     dTauc ;
    double     dFriction ;
    double     dCohesion ;
    double     dX3 ;

    double     dSource ;
    double     dImage ;
    double     dImpulse ;

    double     dTau ;
    double     dTaun ;
    double     dTaus ;
    ST_TENSOR  stTen ;
    double     dZero = 0.0 ;

    double     dGammaDot ;

    dN[0] = dcDip * dcStrike ;
    dN[1] = -dcDip * dsStrike ;
    dN[2] = -dsDip ;

dR[0] = dsStrike * dCr + dcStrike * dsDip * dSr ;
    dR[1] = dcStrike * dCr - dsStrike * dsDip * dSr ;
    dR[2] = dcDip * dSr ;

    iIdx3 = blockIdx.x ;
    iIdx2 = blockIdx.y ;
    iIdx1 = threadIdx.x ;

    if ((iIdx1 < iSx1) && (iIdx2 < iSx2) && (iIdx3 < iSx3))
    {
        iIdx = (((iIdx3 * iSx2) + iIdx2) * iSx1) + iIdx1 ;
        if (1 == bPresent)
        {
            dMinArray[iIdx] = 1.0e+30 ;
        }

        // These should go in shared memory.
        dX3 = (iIdx3) * dDx3 ;
        iTemp = (fabs(dX3 - dZ) > dLp) && (fabs(dX3 + dZ) > dLp) ;
        if (0 == iTemp)
        {
           dVo = pStruct[iIdx3].gammadot0 ;
           dTauc = pStruct[iIdx3].stressexponent ;
           dFriction = pStruct[iIdx3].friction ;
           dCohesion = pStruct[iIdx3].cohesion ;

           cuShiftedCoordinates (iIdx1, iIdx2, iIdx3, iSx1, iSx2, iSx3,
                                 dDx1, dDx2, dDx3, &dX1, &dX2, &dDum) ;
           iTemp = ((fabs(dX1 - dX) > MAX_NUM(dWp, dLp)) || (fabs(dX2 - dY) > MAX_NUM(dWp, dLp))) ;
           if (0 == iTemp)
           {
                dX2r = (dcStrike * dX1) - (dsStrike * dX2) ;
                dX1s = (dcDip * dX2r) - (dsDip * dX3) ;
                dX1i = (dcDip * dX2r) + (dsDip * dX3) ;

                iTemp = ((fabs(dX1s - dXr) > (7.01 * dDx1)) && (fabs(dX1i - dXr) > (7.01 * dDx1))) ;
                if (0 == iTemp)
                {
                    dX2s = (dsStrike * dX1) + (dcStrike * dX2) ;
                    dX3s = (dsDip * dX2r) + (dcDip * dX3) ;
                    dX3i = (-dsDip * dX2r) + (dcDip * dX3) ;

                    dTemp1 = cuGauss (dX1s - dXr, dDx1) ;
                    dTemp2 = cuOmega ((dX2s - dYr)/dW, dBeta) ;
                    dTemp3 = cuOmega ((dX3s - dZr)/dL, dBeta) ;

                    dSource = dTemp1 * dTemp2 * dTemp3 ;
                    dTemp1 = cuGauss ((dX1i - dXr), dDx1) ;
                    dTemp3 = cuOmega ((dX3i + dZr)/dL, dBeta) ;
                    dImage = dTemp1 * dTemp2 * dTemp3 ;

                    dImpulse = dSource + dImage ;

                    stTen = pSig[iIdx] ;

                    cutDot (&stTen, dN, dT) ;
                    dTaun = cuSum (dT, dN) ;
                    cuMulSub (dTaun, dT, dN, dTs) ;

                                    //replace sqrt
                    dTemp2 = cuSum (dTs, dTs) ;

                    dTaus = sqrt (dTemp2) ;

                    dTemp1 = dTaus + (dFriction * dTaun) - dCohesion ;
                    dTau = MAX_NUM(dZero, dTemp1) ;

                    dTemp1 = cuSum (dTs, dR) ;
                    iTemp = (dTemp1 < 0.0) && (fabs(dRake) < (PI2 * 1.5)) ;
                    if (0 == iTemp)
                    {
                        dTs[0] = dTs[0] / dTaus ;
                        dTs[1] = dTs[1] / dTaus ;
                        dTs[2] = dTs[2] / dTaus ;
                        dGammaDot = dVo * 2.0 * mycuSinh (dTau/dTauc) ;

                        dTemp3 = MIN (dL, dW) ;
                        dTemp1 = ((dTau / dMu / dGammaDot) *
                        (dTemp3 / sqrt (dDx1 * dDx2))) ;

                        if (1 == bPresent)
                        {
                            dMinArray[iIdx] = (float) dTemp1 ;
                        }
                        /*if (1 == bPresent)
                        {
                        iIdx = (((iIdx3 * iSx2) + iIdx2) * (iSx1 + 2)) + iIdx1 ;
                        gpV1[iIdx] += (dGammaDot * dImpulse * dScaling) ;
                        }*/

                        dnTemp[0] = dN[0] * (2.0 * dMu * dImpulse * dGammaDot) ;
                        dnTemp[1] = dN[1] * (2.0 * dMu * dImpulse * dGammaDot) ;
                        dnTemp[2] = dN[2] * (2.0 * dMu * dImpulse * dGammaDot) ;

                        cuTensorDyadProd (&stTen, dTs, dnTemp) ;

                        cuTensorOperate (&pMoment[iIdx], (void *)&stTen, '+') ;
                    }
                }
            }
        }
    }
}

/* --------------------------------- __global__ functions end ---------------------------------- */


/* --------------------------------- host/device functions ------------------------------------- */

__device__ __host__ void print_tensor (ST_TENSOR    *pstT)
{
    printf ("%e,%e,%e,%e,%e,%e\n",pstT->s11, pstT->s12, pstT->s13, pstT->s22, pstT->s23, pstT->s33) ;
}


__host__ __device__ double cuTensorTrace (ST_TENSOR    *pstT)
{
    return (double)(pstT->s11 + pstT->s22 + pstT->s33) ;
}


__host__ __device__ void cuIsotrpicStressStrain (ST_TENSOR    *pstT,
                                                 double       dLambda,
                                                 double       dMu)
{
    double dEpskk ; 
    double dTemp = 0.0 ;
    
    dTemp = 2.0 * dMu ; 
    //dEpskk = cuTensorTrace (pstT) ;
    dEpskk = (double)(pstT->s11 + pstT->s22 + pstT->s33) ;
    cuTensorOperate (pstT, &dTemp, '*') ;
    
    pstT->s11 += (dLambda * dEpskk) ;
    pstT->s22 += (dLambda * dEpskk) ;
    pstT->s33 += (dLambda * dEpskk) ; 
}


__device__ void culocalstrain_fir2 (ST_TENSOR  *pstT,
                                    int        iLen1,
                                    int        iLen2,
                                    int        iLen3,
                                    int        iInd1,
                                    int        iInd2,
                                    int        iInd3,
                                    float      *gpV1,
                                    float      *gpV2,
                                    float      *gpV3,
                                    int        iSx1,
                                    int        iSx2,
                                    int        iSx3)
{
    int     iIndex = 0 ;
    int     iInd1m = 0 ;
    int     iInd2m = 0 ;
    int     iInd3m = 0 ;
    int     iInd1p = 0 ;
    int     iInd2p = 0 ;
    int     iInd3p = 0 ;
    int     iTemp = 0 ;

    int     iOffp = 0 ;
    int     iOffm = 0 ;

    iTemp = iSx1 - 1 ;
    for (iIndex = 0; iIndex < iLen1; iIndex++)
    {
        iInd1m = ((iTemp+iInd1-iIndex) % (iSx1)) ;
        iInd1p = ((iInd1+iIndex) % (iSx1)) + 1 ;
        iInd1p = (iInd1p % iSx1) ;

        iOffp = (((iInd3 * iSx2) + iInd2) * (iSx1 + 2)) + iInd1p ;
        iOffm = (((iInd3 * iSx2) + iInd2) * (iSx1 + 2)) + iInd1m ;

        pstT->s11 += ((gpV1[iOffp] - gpV1[iOffm]) * constdKer1[iIndex]) ;
        pstT->s12 += ((gpV2[iOffp] - gpV2[iOffm]) * constdKer1[iIndex]) ;
        pstT->s13 += ((gpV3[iOffp] - gpV3[iOffm]) * constdKer1[iIndex]) ;
    }
    
    iTemp = iSx2 - 1 ;

    for (iIndex = 0; iIndex < iLen2; iIndex++)
    {
        iInd2m = ((iTemp+iInd2-iIndex) % (iSx2)) ;
        iInd2p = ((iInd2+iIndex) % (iSx2)) + 1 ;
        iInd2p = (iInd2p % iSx2) ;

        iOffp = (((iInd3 * iSx2) + iInd2p) * (iSx1 + 2)) + iInd1 ;
        iOffm = (((iInd3 * iSx2) + iInd2m) * (iSx1 + 2)) + iInd1 ;

        pstT->s12 += ((gpV1[iOffp] - gpV1[iOffm]) * constdKer2[iIndex]) ;
        pstT->s22 += ((gpV2[iOffp] - gpV2[iOffm]) * constdKer2[iIndex]) ;
        pstT->s23 += ((gpV3[iOffp] - gpV3[iOffm]) * constdKer2[iIndex]) ;
    }

    for (iIndex = 1; iIndex <= iLen3; iIndex++)
    {
        iInd3m = iInd3 - iIndex ; 
        iInd3p = iInd3 + iIndex ; 
    
        iOffp = (((iInd3p * iSx2) + iInd2) * (iSx1 + 2)) + iInd1 ;
        iOffm = (((iInd3m * iSx2) + iInd2) * (iSx1 + 2)) + iInd1 ; 

        pstT->s13 += ((gpV1[iOffp] - gpV1[iOffm]) * constdKer3[iIndex-1]) ;
        pstT->s23 += ((gpV2[iOffp] - gpV2[iOffm]) * constdKer3[iIndex-1]) ;
        pstT->s33 += ((gpV3[iOffp] - gpV3[iOffm]) * constdKer3[iIndex-1]) ;
    }

    pstT->s12 /= 2.0 ;
    pstT->s13 /= 2.0 ;
    pstT->s23 /= 2.0 ;
}

__device__ void cuLocalStrain_ani (ST_TENSOR  *pstT, 
                                   int        iInd3m, 
                                   int        iInd3p, 
                                   double     dPx3,
                                   int        iLen1, 
                                   int        iLen2, 
                                   int        iInd1, 
                                   int        iInd2, 
                                   int        iInd3, 
                                   float      *gpV1,
                                   float      *gpV2,
                                   float      *gpV3, 
                                   int        iSx1, 
                                   int        iSx2, 
                                   int        iSx3)
{
    int     iIndex = 0 ;
    int     iInd1m = 0 ;
    int     iInd2m = 0 ;
    int     iInd1p = 0 ;
    int     iInd2p = 0 ;
    int     iTemp = 0 ;
    int     iOffp = 0 ;
    int     iOffm = 0 ;

    iTemp = iSx1 - 1 ; 

    for (iIndex = 0; iIndex < iLen1; iIndex++)
    {
        iInd1m = ((iTemp+iInd1-iIndex) % (iSx1)) ;  
        iInd1p = ((iInd1+iIndex) % (iSx1)) + 1 ;
        iInd1p = (iInd1p % iSx1) ;
        
        iOffp = (((iInd3 * iSx2) + iInd2) * (iSx1 + 2)) + iInd1p ;
        iOffm = (((iInd3 * iSx2) + iInd2) * (iSx1 + 2)) + iInd1m ;
                
        pstT->s11 += ((gpV1[iOffp] - gpV1[iOffm]) * constdKer1[iIndex]) ;   
        pstT->s12 += ((gpV2[iOffp] - gpV2[iOffm]) * constdKer1[iIndex]) ;   
        pstT->s13 += ((gpV3[iOffp] - gpV3[iOffm]) * constdKer1[iIndex]) ;   
    }
    iTemp = iSx2 - 1 ; 
    
    for (iIndex = 0; iIndex < iLen2; iIndex++)
    {
        iInd2m = ((iTemp+iInd2-iIndex) % (iSx2)) ; 
        iInd2p = ((iInd2+iIndex) % (iSx2)) + 1 ;
        iInd2p = (iInd2p % iSx2) ;

        iOffp = (((iInd3 * iSx2) + iInd2p) * (iSx1 + 2)) + iInd1 ;
        iOffm = (((iInd3 * iSx2) + iInd2m) * (iSx1 + 2)) + iInd1 ;

        pstT->s12 += ((gpV1[iOffp] - gpV1[iOffm]) * constdKer2[iIndex]) ;
        pstT->s22 += ((gpV2[iOffp] - gpV2[iOffm]) * constdKer2[iIndex]) ;
        pstT->s23 += ((gpV3[iOffp] - gpV3[iOffm]) * constdKer2[iIndex]) ;
    }

    iOffp = (((iInd3p * iSx2) + iInd2) * (iSx1 + 2)) + iInd1 ;
    iOffm = (((iInd3m * iSx2) + iInd2) * (iSx1 + 2)) + iInd1 ;

    pstT->s13 += ((gpV1[iOffp] - gpV1[iOffm]) / dPx3) ;
    pstT->s23 += ((gpV2[iOffp] - gpV2[iOffm]) / dPx3) ;
    pstT->s33 = ((gpV3[iOffp] - gpV3[iOffm]) / dPx3) ; 

    pstT->s12 /= 2.0 ;
    pstT->s13 /= 2.0 ;
    pstT->s23 /= 2.0 ; 
}

__device__ void cuLocalDivergence_ani (ST_TENSOR  *pstT,
                                       int        iInd3m,
                                       int        iInd3p,
                                       double     dPx3,
                                       int        iLen1,  
                                       int        iLen2,  
                                       int        iInd1,  
                                       int        iInd2,
                                       int        iInd3, 
                                       double     *pF1,
                                       double     *pF2,
                                       double     *pF3,
                                       int        iSx1,                   
                                       int        iSx2,   
                                       int        iSx3)   
{
    int     iIndex = 0 ;
    int     iInd1m = 0 ;
    int     iInd2m = 0 ;
    int     iInd1p = 0 ;
    int     iInd2p = 0 ;
    int     iTemp = 0 ;
    int     iOffp = 0 ;
    int     iOffm = 0 ;

    iTemp = iSx1 - 1 ;

    for (iIndex = 0; iIndex < iLen1; iIndex++)
    {
        iInd1m = ((iTemp+iInd1-iIndex) % (iSx1)) ;
        iInd1p = ((iInd1+iIndex) % (iSx1)) + 1 ;
        iInd1p = (iInd1p % iSx1) ;

        iOffp = (((iInd3 * iSx2) + iInd2) * (iSx1)) + iInd1p ;
        iOffm = (((iInd3 * iSx2) + iInd2) * (iSx1)) + iInd1m ;

        *pF1 += ((pstT[iOffp].s11 - pstT[iOffm].s11) * constdKer1[iIndex]) ;
        *pF2 += ((pstT[iOffp].s12 - pstT[iOffm].s12) * constdKer1[iIndex]) ;
        *pF3 += ((pstT[iOffp].s13 - pstT[iOffm].s13) * constdKer1[iIndex]) ;
    }
    iTemp = iSx2 - 1 ;

    for (iIndex = 0; iIndex < iLen2; iIndex++)
    {
        iInd2m = ((iTemp+iInd2-iIndex) % (iSx2)) ;
        iInd2p = ((iInd2+iIndex) % (iSx2)) + 1 ;
        iInd2p = (iInd2p % iSx2) ;

        iOffp = (((iInd3 * iSx2) + iInd2p) * iSx1) + iInd1 ;
        iOffm = (((iInd3 * iSx2) + iInd2m) * iSx1) + iInd1 ;

        *pF1 += ((pstT[iOffp].s12 - pstT[iOffm].s12) * constdKer2[iIndex]) ;
        *pF2 += ((pstT[iOffp].s22 - pstT[iOffm].s22) * constdKer2[iIndex]) ;
        *pF3 += ((pstT[iOffp].s23 - pstT[iOffm].s23) * constdKer2[iIndex]) ;
    }

    iOffp = (((iInd3p * iSx2) + iInd2) * iSx1) + iInd1 ;
    iOffm = (((iInd3m * iSx2) + iInd2) * iSx1) + iInd1 ;

    *pF1 += ((pstT[iOffp].s13 - pstT[iOffm].s13) / dPx3) ;
    *pF2 += ((pstT[iOffp].s23 - pstT[iOffm].s23) / dPx3) ;
    *pF3 += ((pstT[iOffp].s33 - pstT[iOffm].s33) / dPx3) ;
}

__device__ void cuLocalDivergence_fir (ST_TENSOR  *pstT,
                                       int        iLen1,
                                       int        iLen2,
                                       int        iLen3,
                                       int        iInd1,
                                       int        iInd2,
                                       int        iInd3,
                                       double     *pF1,
                                       double     *pF2,
                                       double     *pF3,
                                       int        iSx1,
                                       int        iSx2,
                                       int        iSx3)
{
    int     iIndex = 0 ;
    int     iIndm = 0 ;
    int     iIndp = 0 ;
    int     iTemp = 0 ;
    int     iOffp = 0 ;
    int     iOffm = 0 ;

    iTemp = iSx1 - 1 ;

    for (iIndex = 0; iIndex < iLen1; iIndex++)
    {
        iIndm = ((iTemp+iInd1-iIndex) % (iSx1)) ;
        iIndp = ((iInd1+iIndex) % (iSx1)) + 1 ;
        iIndp = (iIndp % iSx1) ;

        iOffp = (((iInd3 * iSx2) + iInd2) * (iSx1)) + iIndp ;
        iOffm = (((iInd3 * iSx2) + iInd2) * (iSx1)) + iIndm ;

        *pF1 += ((pstT[iOffp].s11 - pstT[iOffm].s11) * constdKer1[iIndex]) ;
        *pF2 += ((pstT[iOffp].s12 - pstT[iOffm].s12) * constdKer1[iIndex]) ;
        *pF3 += ((pstT[iOffp].s13 - pstT[iOffm].s13) * constdKer1[iIndex]) ;
    }
    iTemp = iSx2 - 1 ;

    for (iIndex = 0; iIndex < iLen2; iIndex++)
    {
        iIndm = ((iTemp+iInd2-iIndex) % (iSx2)) ;
        iIndp = ((iInd2+iIndex) % (iSx2)) + 1 ;
        iIndp = (iIndp % iSx2) ;

        iOffp = (((iInd3 * iSx2) + iIndp) * iSx1) + iInd1 ;
        iOffm = (((iInd3 * iSx2) + iIndm) * iSx1) + iInd1 ;

        *pF1 += ((pstT[iOffp].s12 - pstT[iOffm].s12) * constdKer2[iIndex]) ;
        *pF2 += ((pstT[iOffp].s22 - pstT[iOffm].s22) * constdKer2[iIndex]) ;
        *pF3 += ((pstT[iOffp].s23 - pstT[iOffm].s23) * constdKer2[iIndex]) ;
    }

    for (iIndex = 1; iIndex <= iLen3; iIndex++)
    {
        iIndm = iInd3 - iIndex ;
        iIndp = iInd3 + iIndex ;

        iOffp = (((iIndp * iSx2) + iInd2) * iSx1) + iInd1 ;
        iOffm = (((iIndm * iSx2) + iInd2) * iSx1) + iInd1 ;

        *pF1 += ((pstT[iOffp].s13 - pstT[iOffm].s13) * constdKer3[iIndex-1]) ;
        *pF2 += ((pstT[iOffp].s23 - pstT[iOffm].s23) * constdKer3[iIndex-1]) ;
        *pF3 += ((pstT[iOffp].s33 - pstT[iOffm].s33) * constdKer3[iIndex-1]) ;
    }
}

__host__ __device__ void cuTensorMemset (ST_TENSOR *pstT)
{
    pstT->s11 = 0 ;
    pstT->s12 = 0 ;
    pstT->s13 = 0 ;
    pstT->s22 = 0 ;
    pstT->s23 = 0 ;
    pstT->s33 = 0 ;
}

__host__ __device__ void cuTensorOperate (ST_TENSOR  *pstT,
                                          void       *pTemp,
                                          char       cOp)
{
    double      *dTemp ;
    ST_TENSOR   *pstTemp ;

    switch (cOp)
    {
        case '+':
        {
            pstTemp = (ST_TENSOR *) pTemp ;
            pstT->s11 += pstTemp->s11 ;
            pstT->s12 += pstTemp->s12 ;
            pstT->s13 += pstTemp->s13 ;
            pstT->s22 += pstTemp->s22 ;
            pstT->s23 += pstTemp->s23 ;
            pstT->s33 += pstTemp->s33 ;
        }
        break ;

        case '-':
        {
            pstTemp = (ST_TENSOR *) pTemp ;
            pstT->s11 -= pstTemp->s11 ;
            pstT->s12 -= pstTemp->s12 ;
            pstT->s13 -= pstTemp->s13 ;
            pstT->s22 -= pstTemp->s22 ;
            pstT->s23 -= pstTemp->s23 ;
            pstT->s33 -= pstTemp->s33 ;
        }
        break ;
        case '*':
        {
            dTemp = (double *) pTemp ;
            pstT->s11 *= (*dTemp) ;
            pstT->s12 *= (*dTemp) ;
            pstT->s13 *= (*dTemp) ;
            pstT->s22 *= (*dTemp) ;
            pstT->s23 *= (*dTemp) ;
            pstT->s33 *= (*dTemp) ;
        }
        break ;
        case '=':
        {
            pstTemp = (ST_TENSOR *) pTemp ;
            pstT->s11 = pstTemp->s11 ;
            pstT->s12 = pstTemp->s12 ;
            pstT->s13 = pstTemp->s13 ;
            pstT->s22 = pstTemp->s22 ;
            pstT->s23 = pstTemp->s23 ;
            pstT->s33 = pstTemp->s33 ;
        }
        break ;
    }
}

__host__ __device__ double cuGauss (double dX, 
                                    double dSigma)
{
    double dTemp = 0.0 ;
    
    dTemp = (exp (-0.5 * (dX / dSigma) * (dX / dSigma))) / (sqrt (PI2) * dSigma) ;
    
    return dTemp ;
} 

__host__ __device__ double cuOmega (double  dX,
                                    double  dBeta)
{
    double dTemp = 0.0 ; 
    double dInter = 0.0 ; 
    
    dInter = (1.0 - (2.0 * dBeta)) / (2.0 * (1.0 - dBeta)) ;
  
    if (abs (dX) <= dInter)
    {
        dTemp = 1.0 ;
    }
    else
    {
        if (abs (dX) < (1.0 / (2.0 * (1.0 - dBeta))))
        {
            dInter = cos (PI * ((1.0 - dBeta) * abs (dX) - 0.5 + dBeta) / (2.0 * dBeta)) ;
            dTemp = dInter * dInter ;
            
        }   
        else 
        {
            dTemp = 0.0 ;
        }
    }

    return dTemp ;
}

__host__ __device__ double cuGaussp (double dX,
                                     double dSigma)
{
    double dTemp = 0.0 ;

    dTemp = -dX * (exp (-0.5 * (dX / dSigma) * (dX / dSigma))) / (sqrt (PI2) * (dSigma * dSigma * dSigma)) ;
    return dTemp ;
}

__host__ __device__ double cuOmegap (double      dX,
                                     double      dBeta)
{
    double dTemp = 0.0 ;
    double dInter = 0.0 ;

    dInter = (1.0 - (2.0 * dBeta)) / (2.0 * (1.0 - dBeta)) ;

    if (abs (dX) > dInter)
    {
        if (abs (dX) < (1.0 / (2.0 * (1.0 - dBeta))))
        {
            dTemp = -DSIGN(1.0, dX) * PI * (1.0 - dBeta) / (2.0 * dBeta) * 
            sin (PI * ((1.0 - dBeta) * abs (dX) - 0.5 + dBeta) / dBeta) ;   
        }
    }

    return dTemp ;
}

__host__ __device__ void cuTensorDyadProd (ST_TENSOR *pstT,
                                           double    *pdA,
                                           double    *pdB)
{
    pstT->s11 = pdA[0] * pdB[0] ;
    pstT->s12 = ((pdA[0] * pdB[1]) + (pdA[1] * pdB[0])) / 2.0 ;
    pstT->s13 = ((pdA[0] * pdB[2]) + (pdA[2] * pdB[0])) / 2.0 ;
    pstT->s22 = pdA[1] * pdB[1] ;
    pstT->s23 = ((pdA[1] * pdB[2]) + (pdA[2] * pdB[1])) / 2.0 ; 
    pstT->s33 = pdA[2] * pdB[2] ;

    return ;
}

__host__ __device__ void cuShiftedCoordinates (int     iInd1, 
                                               int     iInd2,
                                               int     iInd3,
                                               int     iSx1,
                                               int     iSx2,
                                               int     iSx3,
                                               double  dDx1, 
                                               double  dDx2, 
                                               double  dDx3, 
                                               double  *pdX1,
                                               double  *pdX2,
                                               double  *pdX3)
{
    if (iInd1 <= iSx1/2)
    {
        *pdX1 = iInd1 * dDx1 ;
    }
    else
    {
        *pdX1 = (iInd1 - iSx1) * dDx1 ;  
    }

    if (iInd2 <= iSx2/2)
    {
        *pdX2 = iInd2 * dDx2 ; 
    }
    else
    {
        *pdX2 = (iInd2 - iSx2) * dDx2 ;
    }

    if (iInd3 <= iSx3/2)
    {
        *pdX3 = iInd3 * dDx3 ; 
    }
    else
    {
        *pdX3 = (iInd3 - iSx3) * dDx3 ;
    }
}

__host__ __device__ double cuTensorNorm (ST_TENSOR *pstTemp)
{
    double dTemp = 0.0 ;
    double dVal = 2.0 ;
    double dSqr = 0.0 ;

/* Took me 2 days to find this bug (Nightmare) */
    double s11 = (double) pstTemp->s11 ;
    double s12 = (double) pstTemp->s12 ;
    double s13 = (double) pstTemp->s13 ;
    double s22 = (double) pstTemp->s22 ;
    double s23 = (double) pstTemp->s23 ;
    double s33 = (double) pstTemp->s33 ;

    dTemp = ((s11 * s11) + 2.0 * (s12 * s12) + 2.0 * (s13 * s13) + (s22 * s22) + 2.0 * 
            (s23 * s23) + (s33 * s33)) / dVal;

    dSqr = sqrt (dTemp) ;

    return dSqr ;
}


__host__ __device__ void cuTensorDecompose (ST_TENSOR   *pstT,
                                            double      *dGamma,
                                            ST_TENSOR   *pstR)
{
    *dGamma = (float)cuTensorNorm (pstT) ;

    pstR->s11 = pstT->s11 / (*dGamma) ;
    pstR->s12 = pstT->s12 / (*dGamma) ;
    pstR->s13 = pstT->s13 / (*dGamma) ;
    pstR->s22 = pstT->s22 / (*dGamma) ;
    pstR->s23 = pstT->s23 / (*dGamma) ;
    pstR->s33 = pstT->s33 / (*dGamma) ;
}

__host__ __device__ void cuTensorDeviatoric(ST_TENSOR *pstSrc,
                                            ST_TENSOR *pstDest)
{
    float dDiag ;

    dDiag = cuTensorTrace (pstSrc) / 3.0 ;

    pstDest->s11 = pstSrc->s11 - dDiag ;
    pstDest->s12 = pstSrc->s12 ;
    pstDest->s13 = pstSrc->s13 ;
    pstDest->s22 = pstSrc->s22 - dDiag ;
    pstDest->s23 = pstSrc->s23 ;
    pstDest->s33 = pstSrc->s33 - dDiag ;
}

__host__ __device__ double cuDgGammaDotNot (ST_WEAK     *pstZones,
                                            int         iN,
                                            double      dX1,
                                            double      dX2,
                                            double      dX3,
                                            double      dBeta)
{
    double   dDg ;
    double   dX ;
    double   dY ;
    double   dZ ;
    double   dL ;
    double   dW ;
    double   dD ;
    double   dStrike ;
    double   dDip ;
    double   dLM ;
    double   dX1s ;
    double   dX2s ;
    double   dX3s ;


    double   dcStrike ;
    double   dsStrike ;
    double   dcDip ;
    double   dsDip ;

    double   dX2r ;
    double   dXr ;
    double   dYr ;
    double   dZr ;

    double   dWp ;
    double   dLp ;
    double   dDp ;

    double   dDgGammaDot0 = 0.0 ;
    int      iInd = 0 ;

    for (iInd = 0; iInd < iN; iInd++)
    {
        dDg = pstZones[iInd].dgammadot0 ;

        dX = pstZones[iInd].x ;
        dY = pstZones[iInd].y ;
        dZ = pstZones[iInd].z ;

        // ask  
        dL = pstZones[iInd].width ;
        dW = pstZones[iInd].length ;
        dD = pstZones[iInd].thickness ;
        dStrike = pstZones[iInd].strike ;
        dDip = pstZones[iInd].dip ;

        dWp = dW * (1.0 + 2.0 * dBeta) / 2.0 ;
        dLp = dL * (1.0 + 2.0 * dBeta) / 2.0 ;
        dDp = dD * (1.0 + 2.0 * dBeta) / 2.0 ;

        dLM = MAX3(dWp, dLp, dDp) ;

        if ((abs (dX3-dZ) > dLM) || (abs(dX1-dX) > dLM) || (abs(dX2-dY) > dLM))
        {
            continue ;
        }

        dcStrike = cos (dStrike) ;
        dsStrike = sin (dStrike) ;
        dcDip = cos (dDip) ;
        dsDip = sin (dDip) ;

        dX2r = (dcStrike * dX) - (dsStrike * dY) ;
        dXr = (dcDip * dX2r) - (dsDip * dZ) ;
        dYr = (dsStrike * dX) + (dcStrike * dY) ;
        dZr = (dsDip * dX2r) + (dcDip * dZ) ;

        dX2r = (dcStrike * dX1) - (dsStrike * dX2) ;
        dX1s = (dcDip * dX2r) - (dsDip * dX3) ;
        dX2s = (dsStrike * dX1) + (dcStrike * dX2) ;
        dX3s = (dsDip * dX2r) + (dcDip * dX3) ;

        dDgGammaDot0 += (cuOmega((dX1s-dXr)/dD, dBeta) * cuOmega((dX2s-dYr)/dW, dBeta) *
                         cuOmega((dX3s-dZr)/dL, dBeta) * dDg) ;
    }

    return dDgGammaDot0 ;
}

__host__ __device__ void cutDot (ST_TENSOR  *pstTemp,
                                 double     *dN,
                                 double     *dT)
{
    dT[0] = pstTemp->s11 * dN[0] + pstTemp->s12 * dN[1] + pstTemp->s13 * dN[2] ;
    dT[1] = pstTemp->s12 * dN[0] + pstTemp->s22 * dN[1] + pstTemp->s23 * dN[2] ;
    dT[2] = pstTemp->s13 * dN[0] + pstTemp->s23 * dN[1] + pstTemp->s33 * dN[2] ;
}


__host__ __device__ double cuSum (double *dN,
                                  double *dT)
{
    return (dN[0] * dT[0] + dN[1] * dT[1] + dN[2] * dT[2]) ;
}

__host__ __device__ void cuMulSub (double  dTaun,
                                   double  *dT,
                                   double  *dN,
                                   double  *dTs)
{
    dTs[0] = dT[0] - (dTaun * dN[0]) ;
    dTs[1] = dT[1] - (dTaun * dN[1]) ;
    dTs[2] = dT[2] - (dTaun * dN[2]) ;
}

__host__ __device__ double mycuSinh (double dX)
{
    double dTemp ;

    if (fabs(dX) > 11.0)
    {
        dTemp = (dX > 0 ? 1.0 : -1.0) * sinh (11.0) ;
    }
    else
    {
        dTemp = sinh(dX) ;
    }

    return dTemp ;
}

/* ---------------------------------- host/device functions end -------------------------------------- */

#endif

/* EOF */
