/*-----------------------------------------------------------------------
! Copyright 2013 Sylvain Barbot
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


#include <stdio.h>
//#define PAPI_PROF
#include "cuinclude.h"
#include <cufft.h>

#ifdef USING_CUDA

//#define ENABLE_FFTW3

#ifdef ENABLE_FFTW3
extern "C"
{
	extern void __fourier_MOD_fft3 (float *data,int *sx1,int *sx2,int *sx3,double *dx1,double *dx2, double *dx3,int *direction) ;
	extern void __fourier_MOD_fft2 (float *data,int *sx1,int *sx2,double *dx1,double *dx2, int *direction) ;
}
#endif

/* --------------------------- Forward declaration -------------------------------- */

/* ------------------------------- Global Functions ------------------------------- */ 

__global__ void scaleAndSub (float *, float *, float, int, int, int) ;

__global__ void cuAddKernel (float *, float *, float *, float2 *, float2 *, float2 *, int, int, int, 
                             int, int) ;

__global__ void scale1D (float *, float *, float *, float) ;

__global__ void cuElasticResKernel (float *, float *, float *, int, int, int, double, double, double, 
                                    double, double) ;

__global__ void cerrutiKernel (double2, double2, double2, float2 *, float2 *, float2 *, int, double, 
                               double, double, double, double, double) ;

__global__ void cuSurfaceReduction(double, double, double, double, double, double, double, float,
                                   float, float, int, int, int, float, float, float) ;

__global__ void cuSurfaceKernel (double, double, double, double, double, double, double, float [], 
                                 float [], float [], int, int, int, float *, float *, float *) ;

__global__ void cuCerrutiModKernel (double, double, double, float [], float [], float [], int, int,                                     int, double, double, double, float *, float *, float *) ;

/* ------------------------------- Global Functions ------------------------------- */ 

/* ------------------------------- Device Functions ------------------------------- */ 

__host__ __device__ void cuWaveNumber (int, int, int, int, int, int, double, double, double, 
                                       double *, double*, double *) ;

__device__ void cerrutiSolution (double, double2, double2, double2, double, double, float2 *, 
                                 float2 *, float2 *, double, double, double, double, int) ;

__device__ void cuCerrutiModDevice (double, double2, double2, double2, double, double, float2 *, 
                                    float2 *, float2 *, double, double, double) ;

__device__ void surfaceCal (double, double, double, double, double, double, double, float, float, 
                            float, int, int, int, float, float, float, unsigned int, unsigned int, 
                            unsigned int, unsigned int, unsigned int) ;

/* ------------------------------- Device Functions ------------------------------- */ 

/* ------------------------------- Host Functions --------------------------------- */ 
static cudaError scaleAll (float, float, int, int, int) ;

static int calfft (float *, float *, float *, int, int, double, double, float *, float *, float *) ;

static void OneDfft (cufftHandle, float *, float *, float *, int, double) ;

static cudaError_t doSurfaceTraction (double, double, double, int, int, int, double, double, 
                                      double) ;

static void cerruti (float *, float *, float *, int, int, int, double, double, double, double, double,                     double, double, float *, float *, float *) ;

static void cuCerrutiModified (double, double, double, double, double, double, float [], float [], 
                               float [], int, int, int, float *, float *, float *) ;

static int allocAndCopy (float *, float *, float *, int, float *, float *, float *, int, float, float, int, int,
                         int) ;

static int  inverseFFT (float *, float *, float *, int, int, int, double, double, double, float *, 
                        float *, float *) ;


/* ------------------------------- Host Functions --------------------------------- */ 
/* -------------------------------------------------------------------------------- */


/* ------------------------------- util functions --------------------------------- */

__device__ double2 comp_double_add (double2, double2) ;
__device__ double2 comp_double_mult (double2, double2) ;
__device__ double2 comp_double_sub (double2, double2) ;
__host__ __device__ double2 make_comp_double (double, double) ;
__device__ double2 comp_double_div (double2, double2) ;
__device__ float comp_float_real (float2) ;
__device__ float comp_float_img (float2) ;
__device__ float2 make_comp_float (float, float) ;
__device__ float2 comp_float_add (float2, float2) ;
__device__ float2 comp_float_sub (float2, float2) ;
__device__ float2 comp_float_mult (float2, float2) ;
__device__  float2 comp_float_div (float2, float2) ;

/* ---------------------------------------------------------------------------------*/


/* ------------------------------ Global variables --------------------------------- */
extern float    *gpV1 ;
extern float    *gpV2 ;
extern float    *gpV3 ;

extern float    *pfDevTract1 ;
extern float    *pfDevTract2 ;
extern float    *pfDevTract3 ;

extern int      ihSx1 ;
extern int      ihSx2 ;
extern int      ihSx3 ;

cufftHandle     ghPlan ;		    /* This holds the current FFT plan */
cufftComplex   	*gpCompData ;		/* Pointer to the device data */
cufftHandle 	ghThreeFftPlan ;	/* Plan which is used for 3 continuous FFT's */
cufftComplex 	*gpCompData1 ;		/* Pointer to the device data */
cufftComplex 	*gpCompData2 ;		/* Pointer to the device data */
cufftComplex 	*gpCompData3 ;		/* Pointer to the device data */
int 			giCurDirection ;	/* The current direction of FFT - > Forward or Inverse */
int 			giDim ;			    /* The dimensions of the current FFT -> 1 for 1D,2 for 2D and 3 for 3D */
cudaStream_t  	cuStreamMem ;		/* Stream in which the memory copy should take place */
//cudaStream_t  cuStreamExec ;		/* Stream in which the execution of FFT's should take place */
cufftComplex    *gpComp2dData1 ; 
cufftComplex 	*gpComp2dData2 ; 
cufftComplex 	*gpComp2dData3 ;
cufftHandle     ghTwoFftPlan ;
cufftHandle 	ghPlanOne ; 
float 			*pfB1 ;
float 			*pfB2 ;
float 			*pfB3 ;
cufftHandle     cuThreeFftPlan ;
cufftHandle     cuTwoFftPlan ;
cufftHandle     hInversePlan ;
cufftHandle     hInvTwoPlan ;
/* -----------------------------------------------------------------------------------*/

/* ------------------------------  Main Functions ----------------------------------- */
int createPlanForFFT(int iSx1, 
                     int iSx2, 
                     int iSx3)
{
    cufftResult     cuRet ;

    cuRet = cufftPlan3d(&cuThreeFftPlan, iSx3, iSx2, iSx1, CUFFT_R2C) ;
    if (CUFFT_SUCCESS != cuRet)
    {
	    printf ("createPlanForFFT : the plan creation failed 1\n") ;
        return 1 ;
    }

    cuRet = cufftPlan2d(&cuTwoFftPlan, iSx2, iSx1, CUFFT_R2C) ;
    if (CUFFT_SUCCESS != cuRet)
    {
        printf ("createPlanForFFT : the plan creation failed 2\n") ;
        cufftDestroy (cuThreeFftPlan) ;
        return 1 ;
    }
    cuRet = cufftPlan3d(&hInversePlan, iSx3, iSx2, iSx1, CUFFT_C2R) ;
    if (CUFFT_SUCCESS != cuRet)
    {
        printf ("createPlanForFFT : the plan creation failed 3\n") ;
        cufftDestroy (cuThreeFftPlan) ;
        cufftDestroy (cuTwoFftPlan) ;
        return 1 ;
    }

    cuRet = cufftPlan2d(&hInvTwoPlan, iSx2, iSx1, CUFFT_C2R) ;
    if (CUFFT_SUCCESS != cuRet)
    {
        printf ("createPlanForFFT : the plan creation failed 4\n") ;
        cufftDestroy (cuThreeFftPlan) ;
        cufftDestroy (cuTwoFftPlan) ;
        cufftDestroy (hInversePlan) ;
        return 1 ;
    }

    return 0 ;
}

int destroyPlanForFFT ()
{
    if (CUFFT_SUCCESS != cufftDestroy (cuThreeFftPlan))
    {
        printf ("destroyPlanForFFT: Destroying the plan failed 1 \n") ;
    } 
	
    if (CUFFT_SUCCESS != cufftDestroy (cuTwoFftPlan))
    {
        printf ("destroyPlanForFFT: Destroying the plan failed 2\n") ;
    }
    
    if (CUFFT_SUCCESS != cufftDestroy (hInversePlan))
    {
        printf ("destroyPlanForFFT: Destroying the plan failed 3\n") ;
    }

    if (CUFFT_SUCCESS != cufftDestroy (hInvTwoPlan))
    {
        printf ("destroyPlanForFFT: Destroying the plan failed 4 \n") ;
    }
	
    return 0 ;
}

/**
 * @breif   This function is used to create a FFT plan and if the plan creation was 
 *		    successul then allocate memory to copy the input for the FFT. 
 *
 * @param 	iSx1[in]        The size of array in x direction.
 * @param 	iSx2[in]        The size of array in y direction.
 * @param 	iSx3[in]        The size of array in z direction.
 * @param 	iDim[in]        The number of dimensions i.e., 1D 2D or 3D.
 * @param	iDirection[in]  The direction of the FFT i.e., Forward or Inverse.
 * @return 	iRet[in,out] 	Returns 1 if successful and 0 if unsuccessful.
 *		 
 */

extern "C" void initialize_cufft_ (int iSx1, 
                                   int iSx2, 
                                   int iSx3,
                                   int iDim,
                                   int iDirection,
                                   int *iRet)
{
    cufftResult cuRet ;
    cudaError_t cuError = cudaSuccess;
    *iRet = 1 ;

    switch (iDim)
    {
        case E_ONE_DIMENSION:
        break ;

        case E_TWO_DIMENSION :
        {
            if (iDirection == FFT_FORWARD)
            {
                cuRet = cufftPlan2d(&ghPlan, iSx2, iSx1, CUFFT_R2C) ;
                if (CUFFT_SUCCESS == cuRet)
                {
                    cuError = cudaMalloc((void**)&gpCompData, 
                              sizeof(cufftComplex) * ((iSx1/2) + 1) * iSx2) ;
                }
                else
                {
                    *iRet = 0 ;
                }
            }
            else
            {
                cuRet = cufftPlan2d(&ghPlan, iSx2, iSx1, CUFFT_C2R) ;
                if (CUFFT_SUCCESS == cuRet)
                {
                    cuError = cudaMalloc((void**)&gpCompData, 
                              sizeof(float) * (iSx1+2) * iSx2) ;
                }
                else
                {
                    *iRet = 0 ;
                }
            }
            if (cudaSuccess != cuError)
            {
                *iRet = 0 ;
                fprintf(stderr, "Cuda error: Failed because :  %s\n",
                cudaGetErrorString(cuError)) ;
            }
            else
            {
                giCurDirection = iDirection ;
                giDim = iDim ;
                *iRet = 1 ;
            }
        }
        break ;
                    
        case E_THREE_DIMENSION :
        {
#ifdef GPU_MEMORY_LOG
            size_t iFreeMem = 0 ;
            size_t iTotalMem = 0 ;
            cudaMemGetInfo(&iFreeMem, &iTotalMem);  
            fprintf(stdout, "Memory avaliable: Free: %lu MB, Total: %lu MB \n", 
            iFreeMem/(1024 * 1024), iTotalMem/(1024 * 1024));
#endif
                
            if (iDirection == FFT_FORWARD)
            {
                cuRet = cufftPlan3d(&ghPlan, iSx3, iSx2, iSx1, CUFFT_R2C) ;
                if (CUFFT_SUCCESS == cuRet)
                {
#ifdef GPU_MEMORY_LOG
                    cudaMemGetInfo(&iFreeMem, &iTotalMem);
                    fprintf(stdout, "Memory avaliable: Free: %lu MB, Total: %lu MB \n", 
                    iFreeMem/(1024 * 1024), iTotalMem/(1024 * 1024));
                    fprintf (stdout, "Going to allocate %lu MB\n", 
                    (sizeof(cufftComplex) * ((iSx1/2) + 1) * iSx2 * iSx3) /(1024 * 1024)) ; 
#endif
                    cuError = cudaMalloc((void**)&gpCompData, 
                              sizeof(cufftComplex) * ((iSx1/2) + 1) * iSx2 * iSx3) ;
                }
                else
                {
                    fprintf(stderr, "Plan failed\n"); 
                    *iRet = 0 ;
                }
            }
            else
            {
                cuRet = cufftPlan3d(&ghPlan, iSx3, iSx2, iSx1, CUFFT_C2R) ;
                if (CUFFT_SUCCESS == cuRet)
                {
                    cuError = cudaMalloc((void**)&gpCompData, sizeof(float) * (iSx1+2) * iSx2 * iSx3) ;
                }
                else
                {
                    fprintf(stderr, "Plan failed 2\n") ;
                    *iRet = 0 ;
                }	
            }
            if (cuError != cudaSuccess)
            {
                *iRet = 0 ;
                fprintf(stderr, "Cuda error: Failed because :  %s\n",cudaGetErrorString(cuError)) ;
            }
            else
            {
                giCurDirection = iDirection ;
                giDim = iDim ;
                *iRet = 1 ;
            }
        }
        break ;
    }
		
	return ;
}

/**
 * @breif   This function is used to execute the FFT using the plan already created.  
 *          After the execution of the FFT it also uses the scaling kernel to scale FFT output. 
 *
 * @param	fData[in,out]	The array which contains the input and we do in-place transform 
 *				            so it contains the FFT output as well.
 * @param 	fScale[in]	    This is scaling factor. 	
 * @param   iSx1[in]        The size of array in x direction.
 * @param   iSx2[in]        The size of array in y direction.
 * @param   iSx3[in]        The size of array in z direction.
 *               
 */

extern "C" void calculatefft_(float fData[],
                              float fScale,
                              int   iSx1, 
                              int   iSx2, 
                              int   iSx3) 
{
 	cufftResult_t   cuRet ;
	cudaError_t     cuError ;

    switch (giDim)
    {
        case E_ONE_DIMENSION:
        break ;
                
        case E_TWO_DIMENSION :
        {	
            if (FFT_FORWARD == giCurDirection)
            {
                cuError = cudaMemcpy (gpCompData, fData, 
                          sizeof (float) * (iSx1 + 2) * iSx2, cudaMemcpyHostToDevice) ;
                CHECK_CUDA_ERROR("calculatefft", EXIT_FUNCTION)

                cuRet = cufftExecR2C(ghPlan, (cufftReal *)gpCompData, gpCompData) ;
                if (CUFFT_SUCCESS != cuRet)
                {
                    fprintf(stderr, "CUFFT error: ExecR2C Forward failed");
                    return ;
                }
                if (cudaDeviceSynchronize() != cudaSuccess)
                {
                    fprintf(stderr, "Cuda error: Failed to synchronize\n") ;
                    return ;
                }
                else
                {
                    /*Success case : Everything went well*/
                    dim3 dimGrid(iSx2, 1, 1) ;
                    dim3 dimBlock((iSx1 + 2), 1, 1) ;
                    scaling <float> <<<dimGrid, dimBlock>>> ((float *)gpCompData, fScale, 
                                                             (iSx1 + 2), iSx2, 2) ; 
                    cuError = cudaMemcpy (fData, gpCompData, 
                    sizeof(cufftComplex) * ((iSx1/2) + 1) * iSx2, cudaMemcpyDeviceToHost) ;
                    CHECK_CUDA_ERROR("calculatefft", EXIT_FUNCTION)
                }
            }
            else
            {
                cuError = cudaMemcpy (gpCompData, fData, 
                          sizeof (cufftComplex) * (iSx1/2 + 1) * iSx2, cudaMemcpyHostToDevice) ;
                CHECK_CUDA_ERROR("calculatefft", EXIT_FUNCTION)
                                   
                cuRet = cufftExecC2R(ghPlan, gpCompData, (cufftReal *)gpCompData) ; 
                if (CUFFT_SUCCESS != cuRet)
                {
                    fprintf(stderr, "CUFFT error: ExecR2C Forward failed");
                    return ;
                }

                if (cudaDeviceSynchronize() != cudaSuccess)
                {
                    fprintf(stderr, "Cuda error: Failed to synchronize\n") ;
                    return ;
                }
                else
                {
                    /*Success case : Everything went well*/
                    dim3 dimGrid(iSx2, 1, 1) ;
                    dim3 dimBlock((iSx1 + 2), 1, 1) ;
                    scaling <float> <<<dimGrid, dimBlock>>> ((float *)gpCompData, fScale, 
                                                             (iSx1 + 2), iSx2, 2) ;
                    cuError = cudaMemcpy (fData, gpCompData, 
                    sizeof(float) * (iSx1+ 2) * iSx2, cudaMemcpyDeviceToHost) ;
                    CHECK_CUDA_ERROR("calculatefft", EXIT_FUNCTION)
                }
            }
        }
        break ;
                    
        case E_THREE_DIMENSION :
        {
            if (FFT_FORWARD == giCurDirection)
            {
                cuError = cudaMemcpy (gpCompData, fData, 
                          sizeof (float) * (iSx1 + 2) * iSx2 * iSx3, cudaMemcpyHostToDevice) ;
                CHECK_CUDA_ERROR("calculatefft", EXIT_FUNCTION)

                cuRet = cufftExecR2C(ghPlan, (cufftReal *)gpCompData, gpCompData) ; 
                if (CUFFT_SUCCESS != cuRet)
                {
                    fprintf(stderr, "CUFFT error: ExecR2C Forward failed");
                    return ;
                }

                if (cudaDeviceSynchronize() != cudaSuccess)
                {
                    fprintf(stderr, "Cuda error: Failed to synchronize\n") ;
                    return ;
                }
                else
                {
                    dim3 dimGrid(iSx3, iSx2, 1) ;
                    dim3 dimBlock((iSx1 + 2), 1, 1) ;
                    /*Success case : Everything went well*/
                    scaling <float> <<<dimGrid, dimBlock>>> ((float *)gpCompData, fScale, (iSx1 + 2), 
                                                             iSx2, iSx3) ;
                    cuError = cudaMemcpy (fData, gpCompData, 
                    sizeof(cufftComplex) * ((iSx1/2) + 1) * iSx2 * iSx3, cudaMemcpyDeviceToHost) ;
                    CHECK_CUDA_ERROR("calculatefft", EXIT_FUNCTION)
                }
            }
            else
            {
                cuError = cudaMemcpy (gpCompData, fData, 
                sizeof (cufftComplex) * (iSx1/2 + 1) * iSx2 * iSx3, cudaMemcpyHostToDevice) ;
                CHECK_CUDA_ERROR("calculatefft", EXIT_FUNCTION)

                cuRet = cufftExecC2R(ghPlan, gpCompData, (cufftReal *)gpCompData) ;
                if (CUFFT_SUCCESS != cuRet)
                {
                    fprintf(stderr, "CUFFT error: ExecR2C Forward failed");
                    return ;
                }

                if (cudaDeviceSynchronize() != cudaSuccess)
                {
                    fprintf(stderr, "Cuda error: Failed to synchronize\n") ;
                    return ;
                }
                else
                {
                    dim3 dimGrid(iSx3, iSx2, 1) ;
                    dim3 dimBlock((iSx1 + 2), 1, 1) ;
                    /*Success case : Everything went well*/
                    scaling <float> <<<dimGrid, dimBlock>>> ((float *)gpCompData, fScale, (iSx1 + 2), 
                                                             iSx2, iSx3) ;
                    cuError = cudaMemcpy (fData, gpCompData, 
                    sizeof(float) * (iSx1+ 2) * iSx2 * iSx3, cudaMemcpyDeviceToHost) ;
                    CHECK_CUDA_ERROR("calculatefft", EXIT_FUNCTION)
                }
            }
        }
        break ;
    }

EXIT_FUNCTION:	
	return ;
}

/**
 * @breif   This function is used to free the memory allocated and 
 *			destroys the FFT plan created. 
 *                               
 */

extern "C" void deinitialize_cufft_ ()
{
	cufftDestroy (ghPlan) ;
	cudaFree (gpCompData) ; 

#ifdef GPU_MEMORY_LOG
    size_t iFreeMem = 0 ; 
    size_t iTotalMem = 0 ;
    cudaMemGetInfo(&iFreeMem, &iTotalMem) ;
    fprintf(stdout, "deinitialize_cufft : Memory avaliable : Free: %lu MB, Total: %lu MB \n", 
    iFreeMem/(1024 * 1024), iTotalMem/(1024 * 1024)) ;
#endif 
    return ;
}


/* -------------------------------------------------------------------------------------------- */


/* ----------------------------------- Cuda Kernels ------------------------------------------- */

/**
 * @breif       This kernel function is used to scale the FFT output. 
 *
 * @param  pCompData[in,out] The array which contains the input to be scaled.
 * @param  fScale[in]      	 This is scaling factor.         
 * @param  iNx[in]        	 The size of array in x direction.
 * @param  iNy[in]        	 The size of array in y direction.
 * @param  iNz[in]        	 The size of array in z direction.
 *               
 
template <class T>
__global__ void scaling (T 	*pCompData,
			            T   fScale,
                        int iNx, 
                        int iNy, 
                        int iNz)
{
	int iX ; 
	int iY ;
	int iZ ; 
	unsigned int 	iIdx ;

 	iX = threadIdx.x ;
    iY = blockIdx.x ;
    iZ = blockIdx.y ;
	
	if ((iX < iNx) && (iY < iNy) && (iZ < iNz))
	{
        iIdx = (((iZ * iNy) + iY) * iNx) + iX ;
        pCompData[iIdx] = pCompData[iIdx] * fScale ;
	}
}*/

/**
 * @breif  This is a device function which calculates the wavenumber.
 *
 * @param  iInd1[in] Indice1.
 * @param  iInd2[in] Indice2.
 * @param  iInd3[in] Indice3.
 * @param  iSx1[in]  The size of array in x direction.
 * @param  iSx2[in]  The size of array in y direction.
 * @param  iSx3[in]  The size of array in z direction.
 * @param  dDx1[in]  Sampling size.
 * @param  dDx2[in]  Sampling size.
 * @param  dDx3[in]  Sampling size.
 * @param  *dK1[out] The wavenumber.
 * @param  *dK2[out] The wavenumber. 
 * @param  *dK3[out] The wavenumber.
 *               
 */

__device__ void cuWaveNumber (int  	  iInd1, 
                              int 	  iInd2, 
                              int	  iInd3,
                              int 	  iSx1, 
                              int 	  iSx2,
                              int 	  iSx3,
                              double  dDx1, 
                              double  dDx2, 
                              double  dDx3, 
                              double  *dK1, 
                              double  *dK2, 
                              double  *dK3)
{
    if (iInd3 < ((iSx3/2)+1))
    {
        *dK3 = ((double)iInd3) / (iSx3 * dDx3) ;
    }
    else
    {
        *dK3 = -((double) (iSx3 - iInd3)) / (iSx3 * dDx3) ;
    } 

    if (iInd2 < ((iSx2 /2) + 1))
    {
        *dK2 = ((double)iInd2) / (iSx2 * dDx2) ;
    }
    else
    {
        *dK2 = -((double) (iSx2 - iInd2)) / (iSx2 * dDx2) ;
    }

    *dK1 = ((double) iInd1) / (iSx1 * dDx1) ; 
	return ;
}

/**
 * @breif  This function can be called from the fortran interface, it copies the 
 *		   data onto the GPU and calculates the elastic response. 
 *
 * @param  fInData1[in,out] Pointer to the data array
 * @param  fInData2[in,out] Pointer to the data array
 * @param  fInData3[in,out] Pointer to the data array.
 * @param  iSx1[in]         The size of array in x direction.
 * @param  iSx2[in]         The size of array in y direction.
 * @param  iSx3[in]         The size of array in z direction.
 * @param  dDx1[in]         Sampling size.
 * @param  dDx2[in]         Sampling size.
 * @param  dDx3[in]         Sampling size.
 *               
 */

extern "C" void cuelasticresp_ (float   fInData1[],
                                float   fInData2[],
                                float   fInData3[],
                                int     iSx1,
                                int     iSx2,
                                int     iSx3,
                                double  dDx1,
                                double  dDx2,
                                double  dDx3,
                                double  dR1,
                                double  dRatio2)
{
	float 		 *fData1 ;
	float 		 *fData2 ;
	float 		 *fData3 ;
	cudaError_t  cuError ; 
	unsigned int iArrSize = 0 ; 

	/*Initialize the array size*/
    iArrSize = sizeof (float) * (iSx1 + 2) * iSx2 * iSx3 ; 

#ifdef PRINT_DEBUG_INFO
    fprintf (stdout, "Entered the cuelasticresp function \n") ; 
#endif

#ifdef GPU_MEMORY_LOG
    size_t iFreeMem = 0 ;
    size_t iTotalMem = 0 ;
    cudaMemGetInfo(&iFreeMem, &iTotalMem) ;
    fprintf(stdout, "Memory avaliable: Free: %lu MB, Total: %lu MB \n", iFreeMem/(1024 * 1024), iTotalMem/(1024 * 1024));

    unsigned int t ;
    t = sizeof (float) * (iSx1 + 2) * iSx2 * iSx3 ;
    fprintf (stdout, "Going to allocate:  %u MB \n ",t / (1024 * 1024)) ;
#endif	

    dim3 dimGrid(iSx3, iSx2, 1) ;
    dim3 dimBlock((iSx1/2 + 1), 1, 1) ;	
	
/*Test code */
/*   	float           *fIn ;
        float           sumBef = 0 ;
        float           sumAft = 0 ;
	int i, j, k = 0 ; 
	
	fIn = (float *)malloc (sizeof (float) * (iSx1 + 2) * iSx2 * iSx3) ;
	if (NULL != fIn)
	{ 	
		for (i = 0; i < (iSx3); i++)
		{
			for (j = 0 ; j < (iSx2); j++)
			{
				for (k = 0; k < (iSx1 + 2); k++)
				{
					fIn[(i * iSx2 + j) * (iSx1 + 2)  + k] = 1 ;
					sumBef += 1 ;
				}	
			}
		}
	}
	else 
	{
		printf ("Memory allocation failed\n") ;
	}
*/

    cuError = cudaMalloc ((void**)&fData1, iArrSize) ;
    if (cudaSuccess != cuError)
    {
        printf ("failed in memory allocation. 1 : %s\n", cudaGetErrorString(cuError)) ;
        return ;
    }

    cuError = cudaMalloc ((void**)&fData2, iArrSize) ;
    if (cudaSuccess != cuError)
    {
        printf ("failed in memory allocation. 2 :  %s\n", cudaGetErrorString (cuError)) ;
        cudaFree (fData1) ;
        return ; 
    }

    cuError = cudaMalloc ((void**)&fData3, iArrSize) ;
    if (cudaSuccess != cuError)
    {
        printf ("failed in memory allocation. 3 : %s\n", cudaGetErrorString (cuError)) ;
        cudaFree (fData1) ;
        cudaFree (fData2) ;
        return ;
    }
	
    cuError = cudaMemcpy (fData1, fInData1, iArrSize, cudaMemcpyHostToDevice) ;
	
	/* -- Test code starts -- */
	/* ---- cuError = cudaMemcpy (fData1, fIn, iArrSize, cudaMemcpyHostToDevice) ; ---  */
	/* -- Test code ends -- */
	
    CHECK_CUDA_ERROR ("cudaMemcpy1 of cuelasticresp_", EXIT_WITH_FREE)
	/* -- Test code starts -- */
	/* -- memset (fIn, 0, iArrSize) ; -- */
	/* -- Test code ends -- */

    cuError = cudaMemcpy (fData2, fInData2, iArrSize, cudaMemcpyHostToDevice) ;
    CHECK_CUDA_ERROR ("cudaMemcpy2 of cuelasticresp_", EXIT_WITH_FREE)

    cuError = cudaMemcpy (fData3, fInData3, iArrSize, cudaMemcpyHostToDevice) ;
    CHECK_CUDA_ERROR ("cudaMemcpy3 of cuelasticresp_", EXIT_WITH_FREE)


#ifdef PRINT_DEBUG_INFO 
    printf ("Memcpy done\n") ;
#endif
 
    cuElasticResKernel <<<dimGrid, dimBlock>>> (fData1, fData2, fData3, (iSx1/2 + 1), iSx2, iSx3, 
                                                dDx1, dDx2, dDx3, dR1, dRatio2) ;
    cuError = cudaGetLastError () ;
    CHECK_CUDA_ERROR ("cuElasticResKernel of cuelasticresp_", EXIT_WITH_FREE)

#ifdef PRINT_DEBUG_INFO
    printf ("Launched the kernel\n") ;	
#endif 
        
    if (cudaDeviceSynchronize() != cudaSuccess)
    {
        fprintf(stderr, "Cuda error: Failed to synchronize\n") ;
        goto EXIT_WITH_FREE ;		
    }
	
    cuError = cudaMemcpy (fInData1, fData1, iArrSize, cudaMemcpyDeviceToHost) ;
    CHECK_CUDA_ERROR ("cudaMemcpy4 of cuelasticresp_", EXIT_WITH_FREE)
	
	/* -- Test code starts -- */
	/* -- cuError = cudaMemcpy (fIn, fData1, iArrSize, cudaMemcpyDeviceToHost) ; -- */
	/* -- Test code ends -- */

	/*for (i = 0; i < (iSx3); i++)
        {
        	for (j = 0 ; j < (iSx2); j++)
                {
                	for (k = 0; k < (iSx1 + 2); k++)
                        {
                        	sumAft += fIn[i * (iSx1 + 2)  * iSx2 + j * (iSx1 + 2) + k];
                        }
                }       
        }

	if (2 * sumBef != sumAft) 
	{
		printf ("Something is wrong, The value before is %f and after is %f\n", sumBef, sumAft) ;
	}
	else 
	{
		printf ("They are equal %f and %f", sumBef, sumAft) ; 
	}
	free (fIn) ;
*/

    cuError = cudaMemcpy (fInData2, fData2, iArrSize, cudaMemcpyDeviceToHost) ;
    CHECK_CUDA_ERROR ("cudaMemcpy5 of cuelasticresp_", EXIT_WITH_FREE)
	
    cuError = cudaMemcpy (fInData3, fData3, iArrSize, cudaMemcpyDeviceToHost) ;
    CHECK_CUDA_ERROR ("cudaMemcpy6 of cuelasticresp_", EXIT_WITH_FREE)
        

EXIT_WITH_FREE :

    cudaFree (fData1) ;
    cudaFree (fData2) ;
    cudaFree (fData3) ;

    return ;
}

/**
 * @breif  This is the kernel which calculates the elastic response. 
 *
 * @param  fInData1[in,out] Pointer to equivalent body-forces in the Fourier domain
 * @param  fInData2[in,out] Pointer to equivalent body-forces in the Fourier domain
 * @param  fInData3[in,out] Pointer to equivalent body-forces in the Fourier domain
 * @param  iSx1[in]         The size of array in x direction.
 * @param  iSx2[in]         The size of array in y direction.
 * @param  iSx3[in]         The size of array in z direction.
 * @param  dDx1[in]         Sampling size.
 * @param  dDx2[in]         Sampling size.
 * @param  dDx3[in]         Sampling size.
 *               
 */
__global__ void cuElasticResKernel (float  *fData1, 
                                    float  *fData2, 
                                    float  *fData3, 
                                    int    iSx1,
                                    int    iSx2,
                                    int    iSx3,
                                    double dDx1,
                                    double dDx2,
                                    double dDx3,
                                    double dR1, 
                                    double dRatio2)
{
	/* Variable declaration */
    double 	     dR2 ; 
    double 	     dDenom ; 
    double 	     dK1 ;
    double 	     dK2 ;
    double 	     dK3 ;

    double2      cuC1 ;
    double2      cuC2 ;
    double2      cuC3 ;
    double2      cuBuf1 ; 
    double2      cuBuf2 ;
    double2      cuBuf3 ;
    unsigned int iIdx1 ; 
    unsigned int iIdx2 ; 
    unsigned int iIdx3 ; 
    unsigned int iIdx ;	
#ifdef SHARED_MEMORY_IMPL_ELASTIC	
    unsigned int i ; 
    unsigned int index ; 
    unsigned int index2 ; 

    __shared__ float pfshBuf1[258] ;
    __shared__ float pfshBuf2[258] ;
    __shared__ float pfshBuf3[258] ;
#else
    float2       *pC1 ;
    float2       *pC2 ;
    float2       *pC3 ;
#endif

	/* Get the proper index */
    iIdx1 = threadIdx.x ;
    iIdx2 = blockIdx.y ;
    iIdx3 = blockIdx.x ;

#ifdef PRINT_DEBUG_INFO
    printf ("Entered the kernel\n") ;
#endif
	
	/*Make sure the indices are within the range */
    if ((iIdx1 < iSx1) && (iIdx2 < iSx2) && (iIdx3 < iSx3))
    {
#ifdef SHARED_MEMORY_IMPL_ELASTIC
        if (iIdx1 < (2 * warpSize))
        {
            for (i = 0; i < (iSx1/warpSize); i++)
            {
                index = (iIdx1 * (iSx1 / warpSize)) + i ;
                index2 = (((iIdx3 * iSx2) + iIdx2) * 2 * iSx1) + index ;
                pfshBuf1[index] = fData1[index2] ;
                pfshBuf2[index] = fData2[index2] ;
                pfshBuf3[index] = fData3[index2] ;
            }
        }
        if (iIdx1 == (iSx1 - 1))
        {
            index = 2 * iIdx1  ;
            index2 = (((iIdx3 * iSx2) + iIdx2) * 2 * iSx1) + index ;
            pfshBuf1[index] = fData1[index2] ;
            pfshBuf1[index+1] = fData1[index2+1] ;
            pfshBuf2[index] = fData2[index2] ;
            pfshBuf2[index+1] = fData2[index2+1] ;
            pfshBuf3[index] = fData3[index2] ;
            pfshBuf3[index+1] = fData3[index2+1] ;
        }
        __syncthreads () ;	
#endif
 
#ifdef PRINT_DEBUG_INFO
        printf ("The values of indexes are iIdx1 : %u, iIdx2 : %u, iIdx3 : %u\n", iIdx1, iIdx2, iIdx3) ;
#endif 
        cuWaveNumber (iIdx1, iIdx2, iIdx3, (2 * iSx1) - 2, iSx2, iSx3, dDx1, dDx2, dDx3, &dK1, &dK2, &dK3) ;

        dR2 = dK1 * dK1 + dK2 * dK2 + dK3 * dK3 ; 
        dDenom = (dR1 / (dR2 * dR2)) ;
		
        iIdx = (((iIdx3 * iSx2) + iIdx2) * 2 * iSx1) + 2 * iIdx1 ;

#ifdef PRINT_DEBUG_INFO		
        printf ("The Idx is %u\n", iIdx) ; 	
#endif
/*		float temp = 0 ; 
		temp = fData1[iIdx] ;
		fData1[iIdx] = 2 *  temp;
		temp = fData1[iIdx + 1] ;
		fData1[iIdx + 1] = 2 *  temp;
*/		
		
/*
		//Change for 64 bit access pattern. Need to check if 128 bit is possible ?
		cuC1 = make_comp_double ((double) fData1[iIdx], (double) fData1[iIdx + 1]) ;
        cuC2 = make_comp_double ((double) fData2[iIdx], (double) fData2[iIdx + 1]) ;
        cuC3 = make_comp_double ((double) fData3[iIdx], (double) fData3[iIdx + 1]) ;
*/ 

#ifdef SHARED_MEMORY_IMPL_ELASTIC
        iIdx = 2 * iIdx1 ;
        cuC1 = make_comp_double ((double) pfshBuf1[iIdx], (double) pfshBuf1[iIdx+1]) ;
        cuC2 = make_comp_double ((double) pfshBuf2[iIdx], (double) pfshBuf2[iIdx+1]) ;
        cuC3 = make_comp_double ((double) pfshBuf3[iIdx], (double) pfshBuf3[iIdx+1]) ;
#else
        pC1 = (float2 *) &fData1[iIdx] ;
        pC2 = (float2 *) &fData2[iIdx] ;
        pC3 = (float2 *) &fData3[iIdx] ;

        cuC1 = make_comp_double ((double) ((*pC1).x), (double) ((*pC1).y)) ;
        cuC2 = make_comp_double ((double) ((*pC2).x), (double) ((*pC2).y)) ;
        cuC3 = make_comp_double ((double) ((*pC3).x), (double) ((*pC3).y)) ;
#endif		
		cuBuf1 = comp_double_mult (comp_double_sub(comp_double_mult
                 (make_comp_double(dK2 * dK2 + dK3 * dK3 + dRatio2 * dR2, 0),  cuC1) ,
                 (comp_double_mult(make_comp_double(dK1,0), 
                 comp_double_add(comp_double_mult(make_comp_double(dK2,0),  cuC2), 
                 comp_double_mult (make_comp_double(dK3,0), cuC3))))), (make_comp_double (dDenom, 0)));

        cuBuf2 = comp_double_mult (comp_double_sub(comp_double_mult	
                 (make_comp_double(dK1 * dK1 + dK3 * dK3 + dRatio2 * dR2, 0),  cuC2) ,
                 (comp_double_mult(make_comp_double(dK2,0), 
                 comp_double_add(comp_double_mult(make_comp_double(dK1,0),  cuC1), 
                 comp_double_mult (make_comp_double(dK3,0), cuC3))))), (make_comp_double(dDenom, 0)));

        cuBuf3 = comp_double_mult (comp_double_sub(comp_double_mult
                 (make_comp_double(dK1 * dK1 + dK2 * dK2 + dRatio2 * dR2, 0),  cuC3) ,
                 (comp_double_mult(make_comp_double(dK3,0), 
                 comp_double_add(comp_double_mult(make_comp_double(dK1,0),  cuC1), 
                 comp_double_mult (make_comp_double(dK2,0), cuC2))))), (make_comp_double (dDenom,0))) ;

#ifdef SHARED_MEMORY_IMPL_ELASTIC
        pfshBuf1[iIdx] = (float) cuBuf1.x ;
        pfshBuf1[iIdx+1] = (float) cuBuf1.y ;
		
        pfshBuf2[iIdx] = (float) cuBuf2.x ;
        pfshBuf2[iIdx+1] = (float) cuBuf2.y ;

        pfshBuf3[iIdx] = (float) cuBuf3.x ;
        pfshBuf3[iIdx+1] = (float) cuBuf3.y ;

        __syncthreads () ;

        if (iIdx1 < (2 * warpSize)) 
        {
            for (i = 0; i < (iSx1/warpSize); i++)
            {
                index = (iIdx1 * (iSx1 / warpSize)) + i ;
                index2 = (((iIdx3 * iSx2) + iIdx2) * 2 * iSx1) + index ;
                fData1[index2] = pfshBuf1 [index];
                fData2[index2] = pfshBuf2 [index] ;
                fData3[index2] = pfshBuf3 [index] ;
            }
        }
        if (iIdx1 == (iSx1 - 1))
        {
            index = 2 * iIdx1  ;
            index2 = (((iIdx3 * iSx2) + iIdx2) * 2 * iSx1) + index ;
            fData1[index2] = pfshBuf1 [index];
            fData1[index2+1] = pfshBuf1 [index+1] ;
            fData2[index2] = pfshBuf2 [index] ;
            fData2[index2+1] = pfshBuf2 [index+1] ;
            fData3[index2] = pfshBuf3 [index] ;
            fData3[index2+1] = pfshBuf3 [index+1] ;		
        }	
#else
        (*pC1).x = (float) cuBuf1.x ;
        (*pC1).y = (float) cuBuf1.y ;

        (*pC2).x = (float) cuBuf2.x ;
        (*pC2).y = (float) cuBuf2.y ;

        (*pC3).x = (float) cuBuf3.x ;
        (*pC3).y = (float) cuBuf3.y ;
	
#endif
		/*fData1[iIdx] = (float) cuBuf1.x ;
		fData1[iIdx + 1] = (float) cuBuf1.y ;
 
		fData2[iIdx] = (float) cuBuf2.x ;
        	fData2[iIdx + 1] = (float) cuBuf2.y ;

		fData3[iIdx] = (float) cuBuf3.x ;
        	fData3[iIdx + 1] = (float) cuBuf3.y ; */
    }	
}

/* ------------------------------- Utility functions ------------------------------------------- */

__device__ double2 comp_double_mult(double2 ab, double2 cd)
{
    return make_comp_double (ab.x * cd.x - ab.y * cd.y, ab.x * cd.y + ab.y * cd.x);
}

__device__ double2 comp_double_add(double2 a, double2 b)
{
    return make_comp_double (a.x + b.x, a.y + b.y);
}

__device__ double2 comp_double_sub(double2 a, double2 b)
{
    return make_comp_double (a.x - b.x, a.y - b.y);
}


__device__ double2 make_comp_double (double a, double b)
{
    double2 c ;
    c.x = a  ;
    c.y = b ;
    return c ;
}

__device__ double2 comp_double_div (double2 x, double2 y)
{
    double2 quot;
    double s = (fabs(y.x)) + (fabs(y.y));
    double oos = 1.0 / s;
    double ars = x.x * oos;
    double ais = x.y * oos;
    double brs = y.x * oos;
    double bis = y.y * oos;
    s = (brs * brs) + (bis * bis);
    oos = 1.0 / s;
    quot = make_comp_double (((ars * brs) + (ais * bis)) * oos,
                                 ((ais * brs) - (ars * bis)) * oos);
    return quot;
}

__device__ float comp_float_real (float2 x)
{
    return x.x;
}

__device__ float comp_float_img (float2 x)
{
    return x.y;
}

__device__ float2 make_comp_float (float r, float i)
{
    float2 res;
    res.x = r;
    res.y = i;
    return res;
}

__device__ float2 comp_float_add (float2 x, float2 y)
{
    return make_comp_float (comp_float_real(x) + comp_float_real(y),
                            comp_float_img(x) + comp_float_img(y));
}

__device__ float2 comp_float_sub (float2 x, float2 y)
{
    return make_comp_float (comp_float_real(x) - comp_float_real(y),
                            comp_float_img(x) - comp_float_img(y));
}

__device__ float2 comp_float_mult (float2 x, float2 y)
{
    float2 prod;
    prod = make_comp_float  ((comp_float_real(x) * comp_float_real(y)) -
                             (comp_float_img(x) * comp_float_img(y)),
                             (comp_float_real(x) * comp_float_img(y)) +
                             (comp_float_img(x) * comp_float_real(y)));
    return prod;
}

__device__  float2 comp_float_div (float2 x, float2 y)
{
    float2 quot;
    float s = fabsf(comp_float_real(y)) + fabsf(comp_float_img(y));
    float oos = 1.0f / s;
    float ars = comp_float_real(x) * oos;
    float ais = comp_float_img(x) * oos;
    float brs = comp_float_real(y) * oos;
    float bis = comp_float_img(y) * oos;
    s = (brs * brs) + (bis * bis);
    oos = 1.0f / s;
    quot = make_comp_float (((ars * brs) + (ais * bis)) * oos,
                                ((ais * brs) - (ars * bis)) * oos);
    return quot;
}


/* -------------------------------------------------------------------------------------------------- */



/**
 * @breif  This function allocates the device memory to hold the data
 *		   and copies the host data into it. 
 *
 * @param  fData1[in] Pointer to equivalent body-forces in the Fourier domain
 * @param  fData2[in] Pointer to equivalent body-forces in the Fourier domain
 * @param  fData3[in] Pointer to equivalent body-forces in the Fourier domain
 * @param  iSize[in]  Number of bytes to be allocated.
 *               
 */

static int allocAndCopy (float  *fData1,
                         float  *fData2,
                         float  *fData3,
                         int    iSize,
                         float  *fData4,
                         float  *fData5,
                         float  *fData6,
                         int    iSize2,
                         float  fScale,
                         float  fScale2,
                         int    iSx1,
                         int    iSx2,
                         int    iSx3)
{
	cufftResult cuRet = CUFFT_SUCCESS;
	cudaError_t cuError = cudaSuccess;

    dim3  dimGrid(iSx3, iSx2, 1) ;
    dim3  dimBlock((iSx1 + 2), 1, 1) ;

    dim3  dimGrid1(iSx2, 1, 1) ;
    dim3  dimBlock1((iSx1 + 2), 1, 1) ;

//	cufftHandle     cuThreeFftPlan ;
//        cufftHandle     cuTwoFftPlan ;
//	cudaStream_t 	cuStreamFFT ;	

/*	cudaStreamCreate (&cuStreamFFT) ;
	
	cuRet = cufftPlan3d(&cuThreeFftPlan, iSx3, iSx2, iSx1, CUFFT_R2C) ;
        if (CUFFT_SUCCESS != cuRet)
        {
                return 0 ;
        }

        cuRet = cufftPlan2d(&cuTwoFftPlan, iSx2, iSx1, CUFFT_R2C) ;
        if (CUFFT_SUCCESS != cuRet)
        {
                cufftDestroy (cuThreeFftPlan) ;
                return 0 ;
        }

	cuRet = cufftSetStream (cuThreeFftPlan, cuStreamFFT) ;
        if (CUFFT_SUCCESS != cuRet)
        {
		cufftDestroy (cuThreeFftPlan) ;
	        cufftDestroy (cuTwoFftPlan) ;	
                return 0 ;
        }
*/

	// DO NOT DELETE THIS COMMENT
	// Strange behaviour : When we create the plan here we get wrong results.
	/*cuRet = cufftPlan3d(&cuThreeFftPlan, iSx3, iSx2, iSx1, CUFFT_R2C) ;
        if (CUFFT_SUCCESS != cuRet)
        {
		goto FREE_ALL ;
        }

        cuRet = cufftPlan2d(&cuTwoFftPlan, iSx2, iSx1, CUFFT_R2C) ;
        if (CUFFT_SUCCESS != cuRet)
        {
		cufftDestroy (cuThreeFftPlan) ;
                goto FREE_ALL ;
        }*/

    cuRet = cufftExecR2C(cuThreeFftPlan, (cufftReal *)gpCompData1, gpCompData1) ;
    if (CUFFT_SUCCESS != cuRet)
    {
        printf("CUFFT error: ExecR2C Forward failed 1");
        goto FREE_ALL ;
    }

    cuRet = cufftExecR2C(cuThreeFftPlan, (cufftReal *)gpCompData2, gpCompData2) ;
    if (CUFFT_SUCCESS != cuRet)
    {
        printf("CUFFT error: ExecR2C Forward failed 2");
        goto FREE_ALL ;
    }

	// Not to remove the commented code below : 
	// Just to remember that this is already tried and doesnt improve the performance.
/*	scaling <float> <<<dimGrid, dimBlock>>> ((float *)gpCompData1, fScale, (iSx1 + 2), iSx2, iSx3) ;
	cuError = cudaGetLastError () ;	
	CHECK_CUDA_ERROR ("AllocAndCopy : scale1", FREE_ALL)
*/
        
    cuRet = cufftExecR2C(cuThreeFftPlan, (cufftReal *)gpCompData3, gpCompData3) ;
    if (CUFFT_SUCCESS != cuRet)
    {
        printf("CUFFT error: ExecR2C Forward failed 3");
        goto FREE_ALL ;
    }
	
    cuRet = cufftExecR2C(cuTwoFftPlan, (cufftReal *)gpComp2dData1, gpComp2dData1) ;
    if (CUFFT_SUCCESS != cuRet)
    {
        printf("CUFFT error: ExecR2C Forward failed 4");
        goto FREE_ALL ;
    }

    cuRet = cufftExecR2C(cuTwoFftPlan, (cufftReal *)gpComp2dData2, gpComp2dData2) ;
    if (CUFFT_SUCCESS != cuRet)
    {
        printf("CUFFT error: ExecR2C Forward failed 5");
        goto FREE_ALL ;
    }

    cuRet = cufftExecR2C(cuTwoFftPlan, (cufftReal *)gpComp2dData3, gpComp2dData3) ;
    if (CUFFT_SUCCESS != cuRet)
    {
        printf("CUFFT error: ExecR2C Forward failed 6");
        goto FREE_ALL;
    }
    if (cudaSuccess != cudaDeviceSynchronize ())
    {
        printf ("Sync failed \n") ;
    }
    cuError = scaleAll (fScale, fScale2, iSx1, iSx2, iSx3) ;
    CHECK_CUDA_ERROR ("scaleAll of allocAndCopy", FREE_ALL)
    if (cudaSuccess != cudaDeviceSynchronize ())
    {
        printf ("Sync failed \n") ;
    }

    //cufftDestroy (cuThreeFftPlan) ;
    //cufftDestroy (cuTwoFftPlan) ;
	
    return 1 ;

FREE_ALL: 
    cudaFree (gpComp2dData1) ;
    cudaFree (gpComp2dData2) ;
    cudaFree (gpComp2dData3) ;
    cudaFree (gpCompData1) ;
    cudaFree (gpCompData2) ;
    cudaFree (gpCompData3) ;
	
    return 0 ;
}


cudaError_t memcpyUsingStreams (float 		   *fDest,
                                float 		   *fSrc,
                                int 		   iBytes,
                                cudaMemcpyKind eDirection, 
                                cudaStream_t   *pCuStream) 
{
    int 		iIndex = 0 ;
    cudaError_t	cuError = cudaSuccess;
    int 		iOffset = 0 ;

    iOffset = (iBytes / NUM_STREAMS) ;

    /*Creating streams if not present */
    if (NULL == pCuStream)
    {
        // Need to destroy the streams at the end of fftnelastic function.
        pCuStream = (cudaStream_t *) malloc(NUM_STREAMS * sizeof(cudaStream_t));
        for (iIndex = 0 ; iIndex < NUM_STREAMS; iIndex++)
        {
            cuError = cudaStreamCreate (&pCuStream[iIndex]) ;
        }
    }	

    if (cuError != cudaSuccess)
    {
        cuError = cudaMemcpy (fDest, fSrc, iBytes, eDirection) ; 
    }
    else
    {
        for (iIndex = 0 ; iIndex < NUM_STREAMS; iIndex++)
        {
            iOffset = iIndex * iOffset ;
            cuError = cudaMemcpyAsync (fDest +  iOffset , fSrc + iOffset, iBytes / NUM_STREAMS , eDirection, pCuStream[iIndex]) ;
        }
    }

    if (NULL != pCuStream)
    {	
        for (iIndex = 0 ; iIndex < NUM_STREAMS; iIndex++)
        {
            cuError = cudaStreamDestroy (pCuStream[iIndex]) ;
        }
        free (pCuStream) ;
    }
	
    return cuError ;
}

/**
 * @breif	This function executes the scaling kernel for all the data.
 *  
 */

static cudaError scaleAll (float 	fScale,
                           float 	fScale2,
                           int 		iSx1,
                           int 		iSx2,
                           int 		iSx3)
{
    cudaError_t cuError = cudaSuccess; 
    dim3        dimGrid(iSx3, iSx2, 1) ;
    dim3        dimBlock((iSx1 + 2), 1, 1) ;
    dim3        dimGrid1(1, iSx2, 1) ;
    dim3        dimBlock1((iSx1 + 2), 1, 1) ;


    scaling <float> <<<dimGrid, dimBlock>>> ((float *)gpCompData1, fScale, (iSx1 + 2), iSx2, iSx3) ;
    CHECK_ERROR("scaleAll") 

    scaling <float> <<<dimGrid, dimBlock>>> ((float *)gpCompData2, fScale, (iSx1 + 2), iSx2, iSx3) ;
    CHECK_ERROR("scaleAll") 

    scaling <float> <<<dimGrid, dimBlock>>> ((float *)gpCompData3, fScale, (iSx1 + 2), iSx2, iSx3) ;
    CHECK_ERROR("scaleAll") 

    scaling <float> <<<dimGrid1, dimBlock1>>> ((float *)gpComp2dData1, fScale2, (iSx1 + 2), iSx2, 2) ;
    CHECK_ERROR("scaleAll")

    scaling <float> <<<dimGrid1, dimBlock1>>> ((float *)gpComp2dData2, fScale2, (iSx1 + 2), iSx2, 2) ;
    CHECK_ERROR("scaleAll")

    scaling <float> <<<dimGrid1, dimBlock1>>> ((float *)gpComp2dData3, fScale2, (iSx1 + 2), iSx2, 2) ;
    CHECK_ERROR("scaleAll")

    return cuError ;
}



void copyfftmemory (float *fData1,
                    float *fData2,
                    float *fData3,
                    float *fData4,
                    float *fData5,
                    float *fData6,
                    int   iSize,
                    int   iSize2,
                    int   iDirection) 
{
#ifdef ENABLE_FFTW3
cudaError_t cuError = cudaSuccess ;

    if (1 == iDirection)
    {
        cuError = cudaMemcpy (gpV1, fData1, iSize, cudaMemcpyHostToDevice) ;
        CHECK_CUDA_ERROR("calfft : memcpy error 1\n", COPY_EXIT)

        cuError = cudaMemcpy (gpV2, fData2, iSize, cudaMemcpyHostToDevice) ;
        CHECK_CUDA_ERROR("calfft : memcpy error 1\n", COPY_EXIT)

        cuError = cudaMemcpy (gpV3, fData3, iSize, cudaMemcpyHostToDevice) ;
        CHECK_CUDA_ERROR("calfft : memcpy error 1\n", COPY_EXIT)

        cuError = cudaMemcpy (pfDevTract1, fData4, iSize2, cudaMemcpyHostToDevice) ;
        CHECK_CUDA_ERROR("calfft : memcpy error 1\n", COPY_EXIT)

        cuError = cudaMemcpy (pfDevTract2, fData5, iSize2, cudaMemcpyHostToDevice) ;
        CHECK_CUDA_ERROR("calfft : memcpy error 1\n", COPY_EXIT)

        cuError = cudaMemcpy (pfDevTract3, fData6, iSize2, cudaMemcpyHostToDevice) ;
        CHECK_CUDA_ERROR("calfft : memcpy error 1\n", COPY_EXIT)

    }
    //	else
    if (0 == iDirection)
    {

        cuError = cudaMemcpy (fData1, gpV1, iSize, cudaMemcpyDeviceToHost) ;
        CHECK_CUDA_ERROR("calfft : memcpy error 1\n", COPY_EXIT)

        cuError = cudaMemcpy (fData2, gpV2, iSize, cudaMemcpyDeviceToHost) ;
        CHECK_CUDA_ERROR("calfft : memcpy error 2\n", COPY_EXIT )

        cuError = cudaMemcpy (fData3, gpV3, iSize, cudaMemcpyDeviceToHost) ;
        CHECK_CUDA_ERROR("calfft : memcpy error 3\n", COPY_EXIT)

        cuError = cudaMemcpy (fData4, pfDevTract1, iSize2, cudaMemcpyDeviceToHost) ;
        CHECK_CUDA_ERROR("calfft : memcpy error 1\n", COPY_EXIT)

        cuError = cudaMemcpy (fData5, pfDevTract2, iSize2, cudaMemcpyDeviceToHost) ;
        CHECK_CUDA_ERROR("calfft : memcpy error 1\n", COPY_EXIT)
                        
        cuError = cudaMemcpy (fData6, pfDevTract3, iSize2, cudaMemcpyDeviceToHost) ;
        CHECK_CUDA_ERROR("calfft : memcpy error 1\n", COPY_EXIT)

    }
return ;

COPY_EXIT: 
    printf ("Error\n") ;
#endif
}
/**
 * @breif       This function executes fft and elastic response in asynchronous manner with different streams.
 *  
 */


extern "C" void calfftnelasic_ (float   *fData1,
                                float   *fData2,
                                float   *fData3,
                                float   *fData4,
                                float   *fData5,
                                float   *fData6,
                                int     iSx1,
                                int     iSx2,
                                int     iSx3,
                                double  dDx1,
                                double  dDx2,
                                double  dDx3,
                                double  dLambda,
                                double  dMu,
                                double  dGamma, 
                                int     *iRet)
{
    cudaError_t	cuError ;
    size_t 		iFreeMem = 0 ;
    size_t 		iTotalMem = 0 ; 
    //	size_t 		iReqSize = 0 ;
    size_t	 	iSize = 0 ;
    size_t	 	iSize2 = 0 ;
    dim3 		dimGrid(ihSx3, ihSx2, 1) ;
    dim3 		dimBlock((ihSx1 + 2), 1, 1) ;
    dim3 		dimGrid1(ihSx3, ihSx2, 1) ;
    dim3 		dimBlock1(((ihSx1/2)+1), 1, 1) ;
    float 		fScale ;
    float 		fScale2 ;
    double 		dRatio1 ;
    double 		dRatio2 ;
#ifdef PAPI_PROF	
    char   		cTimerName[17] = "FFT             " ;
#endif

    *iRet = 1 ;
        
    cudaMemGetInfo(&iFreeMem, &iTotalMem) ;
#ifdef PAPI_PROF	
    papistartprofiling_(cTimerName) ;
#endif
    iSx1=ihSx1 ;
    iSx2=ihSx2 ;
    iSx3=ihSx3 ;

    /* Lets check the available memory before we allocate */
    iSize = (sizeof (float) * (iSx1 + 2) * (iSx2) * (iSx3)) ;
    iSize2 = (sizeof (float) * (iSx1 + 2) * (iSx2)) ;

//	copyfftmemory (fData1, fData2, fData3, fData4, fData5, fData6, iSize, iSize2, 1) ;	
//	iReqSize = iSize * 3 + (3 * iSize2) ;
	
//        printf ("Free : %lu\n", iFreeMem) ;
       /* if (iFreeMem < iReqSize)
        {
                printf ("Free : %lu, Required: %lu\n", iFreeMem, iReqSize) ;
                return ;
        }*/
    gpCompData1 = (cufftComplex *) gpV1 ;
    gpCompData2 = (cufftComplex *) gpV2 ;
    gpCompData3 = (cufftComplex *) gpV3 ;

    gpComp2dData1 = (cufftComplex *) pfDevTract1;
    gpComp2dData2 = (cufftComplex *) pfDevTract2;
    gpComp2dData3 = (cufftComplex *) pfDevTract3;

    fScale = (float) (dDx1 * dDx2 * dDx3) ;
    fScale2 = (float) (dDx1 * dDx2) ;

#ifndef ENABLE_FFTW3
	/* Allocation copying and fft */
    *iRet = allocAndCopy (fData1, fData2, fData3, iSize, fData4, fData5, fData6, 
                          iSize2, fScale, fScale2, iSx1, iSx2, iSx3) ;
    if (0 == *iRet)
    {
        printf ("Memory allocation failure 1\n") ;
        return ;        
    }
#endif

#ifdef ENABLE_FFTW3
    copyfftmemory (fData1, fData2, fData3, fData4, fData5, fData6, iSize, iSize2, 0) ;
    int iDirection ; 

#ifdef FFTW3
    iDirection = -1 ;
    __fourier_MOD_fft3 (fData1, &iSx1, &iSx2, &iSx3, &dDx1, &dDx2, &dDx3, &iDirection) ;
    __fourier_MOD_fft3 (fData2, &iSx1, &iSx2, &iSx3, &dDx1, &dDx2, &dDx3, &iDirection) ;
    __fourier_MOD_fft3 (fData3, &iSx1, &iSx2, &iSx3, &dDx1, &dDx2, &dDx3, &iDirection) ;
    iDirection = -1 ;
    __fourier_MOD_fft2 (fData4, &iSx1, &iSx2, &dDx1, &dDx2, &iDirection) ;
    __fourier_MOD_fft2 (fData5, &iSx1, &iSx2, &dDx1, &dDx2, &iDirection) ;
    __fourier_MOD_fft2 (fData6, &iSx1, &iSx2, &dDx1, &dDx2, &iDirection) ;
#endif
    copyfftmemory (fData1, fData2, fData3, fData4, fData5, fData6, iSize, iSize2, 1) ;
#endif

#ifdef PAPI_PROF
    if (cudaDeviceSynchronize() != cudaSuccess)
    {
        printf("Cuda error: Failed to synchronize\n") ;
        goto FREE_AND_EXIT ;
    }	
    papiendprofiling_(cTimerName) ;
        
    strcpy (cTimerName, "Elastic         ") ;
    papistartprofiling_(cTimerName) ;
#endif
    dRatio1 = (dLambda + dMu) / (dLambda + (2 * dMu)) / dMu / (PI2 * PI2) ;
    dRatio2 =  dMu/ (dLambda + dMu) ;
	
	/* So far so good.. lets calculate elastic response */ 
    cuElasticResKernel <<<dimGrid1, dimBlock1>>> ((float *)gpCompData1, (float *)gpCompData2, 
                        (float *)gpCompData3, (iSx1/2 + 1), iSx2, iSx3, dDx1, dDx2, dDx3, dRatio1, 
                    dRatio2) ;
    cuError = cudaGetLastError () ;
    CHECK_CUDA_ERROR ("cuElasticResKernel of calfftnelasic_", FREE_AND_EXIT)
	
    /* This is to turn the NAN to ZERO */
    cuError = cudaMemset (gpCompData1, 0,(2 * sizeof (float))) ;
    CHECK_CUDA_ERROR ("cudaMemset1 of calfftnelasic_", FREE_AND_EXIT)

    cuError = cudaMemset (gpCompData2, 0,(2 * sizeof (float))) ;
    CHECK_CUDA_ERROR ("cudaMemset1 of calfftnelasic_", FREE_AND_EXIT)

    cuError = cudaMemset (gpCompData3, 0,(2 * sizeof (float))) ;
    CHECK_CUDA_ERROR ("cudaMemset1 of calfftnelasic_", FREE_AND_EXIT)

#ifdef PAPI_PROF
    if (cudaDeviceSynchronize() != cudaSuccess)
    {
        printf("Cuda error: Failed to synchronize\n") ;
        goto FREE_AND_EXIT ;
    }	
    papiendprofiling_(cTimerName) ;

    strcpy (cTimerName, "Surface         ") ;
    papistartprofiling_(cTimerName) ;
#endif

    cuError = doSurfaceTraction (dLambda, dMu, dGamma, iSx1, iSx2, iSx3, dDx1, dDx2, dDx3) ;

	// No need to sync here, its there in doSurfacetraction function. 
    if (cuError != cudaSuccess)
    {
        printf ("Something went wrong in doSurfaceTraction check logs\n") ;
        goto FREE_AND_EXIT ;
    }

#ifdef PAPI_PROF
    strcpy (cTimerName, "cerruti         ") ;
    papiendprofiling_(cTimerName) ;
	
    strcpy (cTimerName, "inversefft      ") ;
    papistartprofiling_ (cTimerName) ;
#endif

#ifndef ENABLE_FFTW3
    *iRet = calfft ((float *)gpComp2dData1, (float *)gpComp2dData2, (float *)gpComp2dData3, 
                                             iSx1, iSx2, dDx1, dDx2, fData4, fData5, fData6)  ;
    if (0 != *iRet)
    {
        printf ("calfft failed \n") ;
        goto FREE_AND_EXIT ;
    }

    *iRet = inverseFFT ((float *)gpCompData1, (float *)gpCompData2, (float *)gpCompData3, 
                        iSx1, iSx2, iSx3, dDx1, dDx2, dDx3, fData1, fData2, fData3) ;
    if (0 != *iRet)	
    {
        printf ("inverseFFT failed \n") ;
        goto FREE_AND_EXIT ;
    }
	//copyfftmemory (fData1, fData2, fData3, fData4, fData5, fData6, iSize, iSize2, 0) ;

#endif


#ifdef ENABLE_FFTW3
    copyfftmemory (fData1, fData2, fData3, fData4, fData5, fData6, iSize, iSize2, 0) ;

#ifdef FFTW3
    iDirection = 1 ;
    __fourier_MOD_fft3 (fData1, &iSx1, &iSx2, &iSx3, &dDx1, &dDx2, &dDx3, &iDirection) ;
    __fourier_MOD_fft3 (fData2, &iSx1, &iSx2, &iSx3, &dDx1, &dDx2, &dDx3, &iDirection) ;
    __fourier_MOD_fft3 (fData3, &iSx1, &iSx2, &iSx3, &dDx1, &dDx2, &dDx3, &iDirection) ;
    iDirection = 1 ;
    __fourier_MOD_fft2 (fData4, &iSx1, &iSx2, &dDx1, &dDx2, &iDirection) ; 
    __fourier_MOD_fft2 (fData5, &iSx1, &iSx2, &dDx1, &dDx2, &iDirection) ;
    __fourier_MOD_fft2 (fData6, &iSx1, &iSx2, &dDx1, &dDx2, &iDirection) ;
#endif
    copyfftmemory (fData1, fData2, fData3, fData4, fData5, fData6, iSize, iSize2, 1) ;
#endif //ENABLE_FFTW3

#ifdef PAPI_PROF
    if (cudaDeviceSynchronize() != cudaSuccess)
    {
        printf("Cuda error: Failed to synchronize\n") ;
        goto FREE_AND_EXIT ;
    }
    papiendprofiling_(cTimerName) ;	
#endif	

    *iRet = 0 ;	
    return ;
	
FREE_AND_EXIT :
    cudaFree (gpCompData1) ;
    cudaFree (gpCompData2) ;
    cudaFree (gpCompData3) ;

    cudaFree (gpComp2dData1) ;
    cudaFree (gpComp2dData2) ;
    cudaFree (gpComp2dData3) ;
    return ;
}

static int calfft (float   *pSrcData1,
                   float   *pSrcData2,
                   float   *pSrcData3,
                   int     iSx1,
                   int     iSx2,
                   double  dDx1,
                   double  dDx2,
                   float   *pDestData1,
                   float   *pDestData2,
                   float   *pDestData3)
{
    cufftResult_t   cuRet ;
    //	cufftHandle     hInversePlan ;
    dim3            dimGrid(1, iSx2, 1) ;              
    dim3            dimBlock((iSx1 + 2), 1, 1) ; 
    cudaError_t     cuError = cudaSuccess ;
    int             iRet = 0 ;
    float           fScale = 0 ;
	
	/*cuRet = cufftPlan2d(&hInversePlan, iSx2, iSx1, CUFFT_C2R) ;
        if (CUFFT_SUCCESS != cuRet)
	{
		printf ("calfft : Failed in plan creation\n") ;
		iRet = 1 ;
              	return iRet ; 
	}*/
    cuRet = cufftExecC2R(hInvTwoPlan, (cufftComplex *)pSrcData1, (cufftReal *)pSrcData1) ;
    if (CUFFT_SUCCESS != cuRet)
    {
        printf("calfft : execC2R failed 1\n") ;
        goto CALFFT_EXIT ;
    }

    cuRet = cufftExecC2R(hInvTwoPlan, (cufftComplex *)pSrcData2, (cufftReal *)pSrcData2) ;
    if (CUFFT_SUCCESS != cuRet)
    {                       
        printf("calfft : execC2R failed 2\n") ;
        goto CALFFT_EXIT ;
    }

    cuRet = cufftExecC2R(hInvTwoPlan, (cufftComplex *)pSrcData3, (cufftReal *)pSrcData3) ;
    if (CUFFT_SUCCESS != cuRet)
    {                       
        printf("calfft : execC2R failed 3\n") ;
        goto CALFFT_EXIT ;
    }

    fScale = (float) (1.0 / (iSx1 * dDx1 * iSx2 * dDx2)) ;

    scaling <float> <<<dimGrid, dimBlock>>> ((float *)pSrcData1, fScale, (iSx1 + 2), iSx2, 2) ;
    cuError = cudaGetLastError () ;
    CHECK_CUDA_ERROR ("calfft: scaling ", CALFFT_EXIT)

    scaling <float> <<<dimGrid, dimBlock>>> ((float *)pSrcData2, fScale, (iSx1 + 2), iSx2, 2) ;
    cuError = cudaGetLastError () ;
    CHECK_CUDA_ERROR ("calfft: scaling ", CALFFT_EXIT)

    scaling <float> <<<dimGrid, dimBlock>>> ((float *)pSrcData3, fScale, (iSx1 + 2), iSx2, 2) ;
    cuError = cudaGetLastError () ;
    CHECK_CUDA_ERROR ("calfft: scaling ", CALFFT_EXIT)
        
    // cufftDestroy (hInversePlan) ;
	
    return iRet ;	
	
CALFFT_EXIT:
    //cufftDestroy (hInversePlan) ;
    printf ("Something went wrong in calfft function\n") ;
    iRet = 1 ;
    return iRet ;	
}

static int  inverseFFT (float 	*pSrcData1, 
                        float 	*pSrcData2, 
                        float 	*pSrcData3, 
                        int 	iSx1, 
                        int 	iSx2, 
                        int 	iSx3, 
                        double  dDx1, 
                        double  dDx2, 
                        double  dDx3, 
                        float   *pDestData1,
                        float   *pDestData2,
                        float   *pDestData3)
{
    cufftResult_t   cuRet ;
        //cufftHandle 	hInversePlan ; 
    cudaError_t	cuError = cudaSuccess ;
    int 		iRet = 0 ;
    dim3 		dimGrid(iSx3, iSx2, 1) ;
    dim3 		dimBlock((iSx1 + 2), 1, 1) ; 
    float 		fScale = 0 ;
	
	/*cuRet = cufftPlan3d(&hInversePlan, iSx3, iSx2, iSx1, CUFFT_C2R) ;
        if (CUFFT_SUCCESS != cuRet)
	{
		printf ("inverseFFT : Failed in plan creation\n") ;
		iRet = 1 ;
		return iRet ;
	}*/


    cuRet = cufftExecC2R (hInversePlan, (cufftComplex *)pSrcData1, (cufftReal *)pSrcData1) ;
    if (CUFFT_SUCCESS != cuRet)
    {
        printf("inverseFFT: ExecC2R failed1\n");
        goto INVERSE_EXIT ;
    }
    cuRet = cufftExecC2R (hInversePlan, (cufftComplex *)pSrcData2, (cufftReal *)pSrcData2) ;
    if (CUFFT_SUCCESS != cuRet)
    {
        printf("inverseFFT: ExecC2R failed2\n");
        goto INVERSE_EXIT ;
    }
    cuRet = cufftExecC2R (hInversePlan, (cufftComplex *)pSrcData3, (cufftReal *)pSrcData3) ;
    if (CUFFT_SUCCESS != cuRet)
    {
        printf("inverseFFT: ExecC2R failed3\n");
        goto INVERSE_EXIT ;
    }

    fScale = (float)(1.0/(iSx1*dDx1*iSx2*dDx2*iSx3*dDx3)) ;
 
    scaling <float> <<<dimGrid, dimBlock>>> ((float *)pSrcData1, fScale, (iSx1 + 2), iSx2, iSx3) ;
    cuError = cudaGetLastError () ;
    CHECK_CUDA_ERROR ("inverseFFT: scaling ", INVERSE_EXIT)

	scaling <float> <<<dimGrid, dimBlock>>> ((float *)pSrcData2, fScale, (iSx1 + 2), iSx2, iSx3) ;
    cuError = cudaGetLastError () ;
    CHECK_CUDA_ERROR ("inverseFFT: scaling ", INVERSE_EXIT)
 
    scaling <float> <<<dimGrid, dimBlock>>> ((float *)pSrcData3, fScale, (iSx1 + 2), iSx2, iSx3) ;
    cuError = cudaGetLastError () ;
    CHECK_CUDA_ERROR ("inverseFFT: scaling ", INVERSE_EXIT)
        
//        cufftDestroy (hInversePlan) ;
	
    return iRet ;	

INVERSE_EXIT :
  //      cufftDestroy (hInversePlan) ;
    printf ("Something went wrong in inverseFFT function\n") ;
    iRet = 1 ;
    return iRet ;
}


static void cuCerrutiModified (double  dGamma, 
                               double  dLambda,
                               double  dMu,
                               double  dDx1,
                               double  dDx2,
                               double  dDx3,
                               float   pfP1[],
                               float   pfP2[],
                               float   pfP3[],
                               int     iSx1,
                               int     iSx2,
                               int     iSx3,
                               float   *fData1,
                               float   *fData2,
                               float   *fData3)         
{
    dim3  dimGrid(iSx2, iSx1/2+1, 1) ;
    dim3  dimBlock(iSx3, 1, 1) ;       

	double  dAlpha ;
	
    dAlpha = (dLambda + dMu)/(dLambda + 2 * dMu) ;

    cuCerrutiModKernel <<<dimGrid, dimBlock>>>(dAlpha, dMu, dGamma, pfP1, pfP2, pfP3,
                       iSx1/2+1, iSx2, iSx3, dDx1, dDx2, dDx3, fData1, fData2, fData3) ;

    if( cudaSuccess != cudaGetLastError ())
    {
        printf ("icuCerrutiModified : Cerruti kernel failed\n");	
    }
}

__global__ void cuCerrutiModKernel (double 	dAlpha,
				    double      dMu,
				    double 	dGamma, 
				    float       pfP1[],
                        	    float       pfP2[],
                        	    float       pfP3[],
				    int         iSx1,
                        	    int         iSx2,
                        	    int         iSx3,
				    double      dDx1,
                        	    double      dDx2,
                        	    double      dDx3,
				    float       *fData1,
                        	    float       *fData2,
                        	    float       *fData3)
{
    unsigned int  iIdx1 ;
    unsigned int  iIdx2 ;
    unsigned int  iIdx3 ;
    unsigned int  iIdxSec ;
    unsigned int  iIdx ;
    double        dK1 ;
    double        dK2 ;
    double        dK3 ;
    double2       cuC1 ;
    double2       cuC2 ;
    double2       cuC3 ;
        
    float2 		  *pC1 ;
    float2 		  *pC2 ;
    float2 		  *pC3 ;
        
    float2        compSum1 ;
    float2        compSum2 ;
    float2        compSum3 ;

    iIdx3 = threadIdx.x ;
    iIdx2 = blockIdx.x ;
    iIdx1 = blockIdx.y ;

    if ((iIdx1 < iSx1) && (iIdx2 < iSx2) && (iIdx3 < iSx3))
    {
        cuWaveNumber (iIdx1, iIdx2, iIdx3, (2 * iSx1) - 2, iSx2, iSx3, dDx1, dDx2, dDx3, 
                      &dK1, &dK2, &dK3) ;

        iIdx = (iIdx2 * 2 * iSx1) + 2 * iIdx1 ; 
        iIdxSec = (((iIdx3 * iSx2) + iIdx2) * 2 * iSx1) + 2 * iIdx1 ;

        pC1 = (float2 *) &pfP1[iIdx] ;
        pC2 = (float2 *) &pfP2[iIdx] ;
        pC3 = (float2 *) &pfP3[iIdx] ;

        cuC1 = make_comp_double ((double) (*pC1).x, (double) (*pC1).y) ;
        cuC2 = make_comp_double ((double) (*pC2).x, (double) (*pC2).y) ;
        cuC3 = make_comp_double ((double) (*pC3).x, (double) (*pC3).y) ;
        /*
            cuC1 = make_comp_double ((double) pfP1[iIdx], (double) pfP1[iIdx + 1]) ;
            cuC2 = make_comp_double ((double) pfP2[iIdx], (double) pfP2[iIdx + 1]) ;
            cuC3 = make_comp_double ((double) pfP3[iIdx], (double) pfP3[iIdx + 1]) ;
        */

        cuCerrutiModDevice (dMu, cuC1, cuC2, cuC3, dAlpha, dGamma, &compSum1, &compSum2, &compSum3, 
                            dK1, dK2, dK3) ;

        atomicAdd (&fData1[iIdxSec], (float) compSum1.x) ;
        atomicAdd (&fData1[iIdxSec + 1], (float) compSum1.y) ;

        atomicAdd (&fData2[iIdxSec], (float) compSum2.x) ;
        atomicAdd (&fData2[iIdxSec + 1], (float) compSum2.y) ;

        atomicAdd (&fData3[iIdxSec], (float) compSum3.x) ;
        atomicAdd (&fData3[iIdxSec + 1], (float) compSum3.y) ;
    }
}

__device__ void cuCerrutiModDevice (double   dMu,
                                    double2  cuC1,
                                    double2  cuC2,
                                    double2  cuC3,
                                    double   dAlpha,
                                    double   dGamma,
                                    float2   *ppB1,
                                    float2   *ppB2,
                                    float2   *ppB3,
                                    double   dK1,
                                    double   dK2,
                                    double   dK3)
{
    double   dBeta ;
    double   dH ;
    double 	 dFi ;

    double2  i ;
    double2  cuB1 ;
    double2  cuB2 ;
    double2  cuB3 ;
    double2  cuTemp ;
    double2  cuTemp1 ;
    double2  cuV1 ;
    double2  cuV2 ;
    double2  cuV3 ;

    if ((0 == dK1) && (0 == dK2))
    {
        *ppB1 = make_comp_float (0, 0) ;
        *ppB2 = make_comp_float (0, 0) ;
        *ppB3 = make_comp_float (0, 0) ;
    }
    else
    {
        dBeta = PI2 * sqrt ((dK1 * dK1) + (dK2 * dK2)) ;
        i = make_comp_double (0, 1.0) ;
        dH = dGamma / dBeta ;
        dFi = (2 * dBeta) / (((PI2 * PI2) * (dK3 * dK3)) + (dBeta * dBeta)) ;

        cuTemp = make_comp_double ((1.0 / (2.0 * dMu * dBeta * dBeta * dBeta)) * dFi, 0) ;
        cuB1 = comp_double_mult (cuTemp, make_comp_double ((double)cuC1.x, (double)cuC1.y)) ;
        cuB2 = comp_double_mult (cuTemp, make_comp_double ((double)cuC2.x, (double)cuC2.y)) ; 

        //b3=(beta*p3+i*(1._8-alpha)*(pi2*k1*p1+pi2*k2*p2))/(2._8*alpha*mu*beta**4*(1+h))

        cuB3 = comp_double_div (comp_double_sub (comp_double_mult (
               make_comp_double ((double)cuC3.x, (double)cuC3.y), 
               make_comp_double (dBeta, 0)), comp_double_mult( comp_double_mult (comp_double_add (
               comp_double_mult (make_comp_double (PI2 * dK1, 0), make_comp_double (cuC1.x, cuC1.y)),   
               comp_double_mult (make_comp_double (PI2 * dK2, 0), make_comp_double (cuC2.x, cuC2.y))), 
               make_comp_double (1.0 - dAlpha, 0)), i)), 
               make_comp_double (2.0 * dAlpha * dMu * dBeta * dBeta * dBeta * dBeta * (1.0 + dH),0)) ;

    //tmp=alpha*i*beta*pi2*b3*(1._8-1._8/alpha-i*pi2*k3*fi)
    
        cuTemp = comp_double_mult (comp_double_mult (comp_double_mult (
                 make_comp_double (dAlpha * dBeta * PI2, 0), i), cuB3), 
                 comp_double_sub (make_comp_double (1.0 - (1.0 / dAlpha), 0),
                 comp_double_mult (i, make_comp_double (PI2 * dK3 * dFi, 0)))) ;
 
        cuV1 = comp_double_mult(make_comp_double(dFi, 0), 
               comp_double_mult (cuTemp, make_comp_double (dK1, 0))) ;
        cuV2 = comp_double_mult(make_comp_double(dFi, 0), 
               comp_double_mult (cuTemp, make_comp_double (dK2, 0))) ;
        cuV3 = comp_double_mult (make_comp_double(-1.0, 0), comp_double_mult (comp_double_mult (
               make_comp_double (dAlpha * dBeta * dBeta, 0), cuB3), 
               comp_double_sub (make_comp_double ((1.0 / dAlpha), 0), 
               comp_double_mult (i, make_comp_double (PI2 * dK3 * dFi, 0))))) ;

    /* cuTemp = alpha*(pi2**2)*(b1*k1+b2*k2)*(1._8-(i*pi2*k3*fi)))) */
    /* u1=CMPLX(v1+(-(2._8*beta**2*b1)+ cuTemp *k1)) */
     
        cuTemp = comp_double_mult (comp_double_mult (make_comp_double (dAlpha * PI2 * PI2, 0), 
                 comp_double_add (comp_double_mult (cuB1, make_comp_double (dK1, 0)), 
                 comp_double_mult (cuB2, make_comp_double (dK2, 0)))), 
                 comp_double_sub (make_comp_double (1.0, 0), 
                 comp_double_mult (i, make_comp_double (PI2 * dK3 * dFi, 0)))) ;
    
        cuTemp1 = comp_double_add (cuV1, comp_double_add (comp_double_mult (
                  make_comp_double (-2.0 * dBeta * dBeta, 0), cuB1), 
                  comp_double_mult (cuTemp, make_comp_double(dK1, 0)))) ;
        (*ppB1).x = (float) cuTemp1.x ;
        (*ppB1).y = (float) cuTemp1.y ;

    
        cuTemp1 = comp_double_add (cuV2, comp_double_add (comp_double_mult (
                  make_comp_double (-2.0 * dBeta * dBeta, 0), cuB2), 
                  comp_double_mult (cuTemp, make_comp_double(dK2, 0)))) ;
        (*ppB2).x = (float) cuTemp1.x ;
        (*ppB2).y = (float) cuTemp1.y ;


        cuTemp1 = comp_double_add ( comp_double_mult (cuV3, make_comp_double (dFi, 0)), 
                  comp_double_mult (make_comp_double(dAlpha * dBeta * PI2 * PI2 * dK3 * dFi, 0), 
                  comp_double_add (comp_double_mult (cuB1, make_comp_double (dK1, 0)),
                  comp_double_mult (cuB2, make_comp_double (dK2, 0))))) ;
    
        (*ppB3).x = (float) cuTemp1.x ;
        (*ppB3).y = (float) cuTemp1.y ;
    }
}


/*
 * @brief : This function to calculate the surface traction is called after the fft and elastic response.
 *          There is version of this function which is called from the surfacetractioncowling function, 
 * 	        this is to make sure that when the memory is insufficient the software just works fine.
*/

static cudaError_t doSurfaceTraction (double        dLambda,
                                      double        dMu,
                                      double        dGamma,
                                      int           iSx1,
                                      int           iSx2,
                                      int           iSx3,
                                      double        dDx1,
                                      double        dDx2,
                                      double        dDx3)
{
    cudaError_t cuError ;
    float       *pfP1 ;
    float       *pfP2 ;
    float       *pfP3 ;
    float       fScale ;

    dim3 		dimGrid(iSx2, iSx3, 1) ;
    dim3 		dimBlock((iSx1/2 + 1), 1, 1) ;

    dim3 		dimGrid1(iSx2, 1, 1) ;
    dim3 		dimBlock1((iSx1 + 2), 1, 1) ;

    double      dModulus ;
    double      dAlpha ;
    double      dGravity ;

#ifdef PAPI_PROF        
    char        cTimerName[17] = "Surface         " ;
#endif

#ifdef REDUCTION 
    int         blocks ;
    int         threads ; 

    blocks = iSx1 ;
    threads = iSx2 ;
#endif

    dModulus = dLambda + 2 * dMu ;
    dAlpha = (dLambda + dMu) / (dLambda + 2 * dMu) ;
    dGravity = 2 * dMu * dAlpha * dGamma ;

    fScale = (float)(1.0 / (iSx3 * dDx3)) ;

    cuError = cudaMalloc((void**)&pfP1, ((iSx1+2) * iSx2 * sizeof (float))) ;
    if (cudaSuccess != cuError)
    {
        printf ("doSurfaceTraction: Memory allocation failure 1\n") ;
        return cuError ;
    }

    cuError = cudaMalloc((void**)&pfP2, ((iSx1+2) * iSx2 * sizeof (float))) ;
    if (cudaSuccess != cuError)
    {
        cudaFree (pfP1) ;
        printf ("doSurfaceTraction: Memory allocation failure 2\n") ;
        return cuError ;
    }

    cuError = cudaMalloc((void**)&pfP3, ((iSx1+2) * iSx2 * sizeof (float))) ;
    if (cudaSuccess != cuError)
    {
        cudaFree (pfP1) ;
        cudaFree (pfP2) ;
        printf ("doSurfaceTraction: Memory allocation failure 3\n") ;
        return cuError ;
    }

    cuError = cudaMemset (pfP1, 0,((iSx1+2) * iSx2 * sizeof (float))) ;
    CHECK_CUDA_ERROR ("doSurfaceTraction: cudaMemset1", FREE_DATA_DO_SURFACE)

    cuError = cudaMemset (pfP2, 0,((iSx1+2) * iSx2 * sizeof (float))) ;
    CHECK_CUDA_ERROR ("doSurfaceTraction: cudaMemset2", FREE_DATA_DO_SURFACE)

    cuError = cudaMemset (pfP3, 0,((iSx1+2) * iSx2 * sizeof (float))) ;
    CHECK_CUDA_ERROR ("doSurfaceTraction: cudaMemset3", FREE_DATA_DO_SURFACE)
#ifdef REDUCTION
	
    cuSurfaceReduction <<<blocks, threads>>> (dModulus, dLambda, dMu, dGravity, dDx1, dDx2, dDx3,
                                              pfP1, pfP2, pfP3, (iSx1/2 + 1), iSx2, iSx3,
                                              (float *)gpCompData1, (float *)gpCompData2, 
                                              (float *)gpCompData3) ;
#else 
    cuSurfaceKernel <<<dimGrid, dimBlock>>> (dModulus, dLambda, dMu, dGravity, dDx1, dDx2, dDx3,
                                             pfP1, pfP2, pfP3, (iSx1/2 + 1), iSx2, iSx3,
                                             (float *)gpCompData1, (float *)gpCompData2, 
                                             (float *)gpCompData3) ;

#endif
    cuError = cudaGetLastError () ;
    CHECK_CUDA_ERROR ("doSurfaceTraction : cuSurfaceKernel ", FREE_DATA_DO_SURFACE)

    scaleAndSub <<<dimGrid1, dimBlock1>>> (pfP1, (float *) gpComp2dData1, fScale, 
                                           (iSx1 + 2), iSx2, 2) ;
    cuError = cudaGetLastError () ;
    CHECK_CUDA_ERROR ("doSurfaceTraction : scalingKernel ", FREE_DATA_DO_SURFACE)

    scaleAndSub <<<dimGrid1, dimBlock1>>> (pfP2, (float *) gpComp2dData2, fScale, 
                                           (iSx1 + 2), iSx2, 2) ;
    cuError = cudaGetLastError () ;
    CHECK_CUDA_ERROR ("doSurfaceTraction : scalingKernel ", FREE_DATA_DO_SURFACE)

    scaleAndSub <<<dimGrid1, dimBlock1>>> (pfP3, (float *) gpComp2dData3, fScale, 
                                           (iSx1 + 2), iSx2, 2) ;
    cuError = cudaGetLastError () ;
    CHECK_CUDA_ERROR ("doSurfaceTraction :	scalingKernel ", FREE_DATA_DO_SURFACE)

#ifdef PAPI_PROF        
    if (cudaDeviceSynchronize() != cudaSuccess)
    {
        printf("Cuda error: Failed to synchronize\n") ;
        goto FREE_DATA_DO_SURFACE;
    }

    papiendprofiling_(cTimerName) ;

    strcpy (cTimerName, "cerruti         ") ;
    papistartprofiling_(cTimerName) ;

#endif

    cuCerrutiModified (dGamma, dLambda, dMu, dDx1, dDx2, dDx3, (float *)pfP1, (float *)pfP2,       
                       (float *)pfP3, iSx1, iSx2, iSx3, (float *)gpCompData1, (float *)gpCompData2,
                       (float *)gpCompData3);


    if (cudaDeviceSynchronize() != cudaSuccess)
    {
        printf("Cuda error: Failed to synchronize after cuCerrutiModified call\n") ;
    }
        

FREE_DATA_DO_SURFACE :
    cudaFree (pfP1) ;
    cudaFree (pfP2) ;
    cudaFree (pfP3) ;
        
    return cuError ;	
} 

#ifdef REDUCTION

__global__ void cuSurfaceReduction(double  dModulus,
                                   double  dLambda,
                                   double  dMu,
                                   double  dGravity,
                                   double  dDx1,
                                   double  dDx2,
                                   double  dDx3,
                                   float   pfP1[],
                                   float   pfP2[],
                                   float   pfP3[],
                                   int     iSx1,
                                   int     iSx2,
                                   int     iSx3,
                                   float   *fData1,
                                   float   *fData2,
                                   float   *fData3)
{
	unsigned int iIdx ;
	unsigned int iStride ;
	unsigned int iElements ;
	unsigned int iIdxSec ;
	unsigned int iIdx1 ;
	unsigned int iIdx2 ;
	unsigned int iIdx3 ;

    iStride = 2  ;

    if ((threadIdx.x < blockDim.x) && (blockIdx.x < gridDim.x))
    {	
        iIdx = (blockIdx.x * (iSx1 * 2))  +  (threadIdx.x * (iSx2 * iSx1 * 2)) ;
        if (iIdx < (iSx2 * iSx3 * iSx1 * 2))
        {
            iElements = iIdx + (iSx1 * 2) ; 
            iIdxSec = blockIdx.x * (iSx1 * 2)  ;
        
            iIdx1 = 0 ;
            iIdx2 = blockIdx.x ;
            iIdx3 = threadIdx.x ;

            while (iIdx < iElements)
            {
                surfaceCal (dModulus, dLambda, dMu, dGravity, dDx1, dDx2, dDx3,
                            pfP1, pfP2, pfP3, iSx1, iSx2, iSx3,
                            fData1, fData2, fData3,
                            iIdx, iIdxSec, iIdx1, iIdx2, iIdx3) ;
                    
                iIdx += iStride ;
                iIdxSec += iStride ;
                iIdx1 += 1 ;	
            }
        }
    }
}

__device__ void surfaceCal (double        dModulus,
                            double        dLambda,
                            double        dMu,
                            double        dGravity,
                            double        dDx1,
                            double        dDx2,
                            double        dDx3,
                            float         pfP1[],
                            float         pfP2[],
                            float         pfP3[],
                            int           iSx1,
                            int           iSx2,
                            int           iSx3,
                            float         *fData1,
                            float         *fData2,
                            float         *fData3, 
				            unsigned int  iIdx, 
				            unsigned int  iIdxSec, 
				            unsigned int  iIdx1,
				            unsigned int  iIdx2,
				            unsigned int  iIdx3) 
{
    double   dK1 ;
    double   dK2 ;
    double   dK3 ;

    double2  cuC1 ;
    double2  cuC2 ;
    double2  cuC3 ;

    double2  compSum1 ;
    double2  compSum2 ;
    double2  compSum3 ;
        
    cuWaveNumber (iIdx1, iIdx2, iIdx3, (2 * iSx1) - 2, iSx2, iSx3, dDx1, dDx2, dDx3, 
                  &dK1, &dK2, &dK3) ;

    cuC1 = make_comp_double ((double) fData1[iIdx], (double) fData1[iIdx + 1]) ;
    cuC2 = make_comp_double ((double) fData2[iIdx], (double) fData2[iIdx + 1]) ;
    cuC3 = make_comp_double ((double) fData3[iIdx], (double) fData3[iIdx + 1]) ;

    compSum1 = comp_double_mult (comp_double_mult (make_comp_double (0, PI2), 
               make_comp_double (dMu, 0)), comp_double_add (comp_double_mult (
               make_comp_double (dK3, 0), cuC1), comp_double_mult (make_comp_double (dK1, 0), 
               cuC3))) ;

    compSum2 = comp_double_mult (comp_double_mult (make_comp_double (0, PI2), 
               make_comp_double (dMu, 0)), comp_double_add (comp_double_mult (
               make_comp_double (dK3, 0), cuC2), comp_double_mult (make_comp_double (dK2, 0), 
               cuC3))) ;

    compSum3 = comp_double_sub (comp_double_mult (make_comp_double (0, PI2),
               comp_double_add (comp_double_mult (comp_double_mult (make_comp_double 
               (dModulus, 0), make_comp_double(dK3, 0)), cuC3),
               comp_double_mult (make_comp_double (dLambda, 0),
               comp_double_add (comp_double_mult (make_comp_double (dK1, 0), cuC1),
               comp_double_mult (make_comp_double (dK2, 0), cuC2))))),
               comp_double_mult (make_comp_double (dGravity, 0), cuC3)) ;

    atomicAdd (&pfP1[iIdxSec], (float) compSum1.x) ;
    atomicAdd (&pfP1[iIdxSec + 1], (float) compSum1.y) ;

    atomicAdd (&pfP2[iIdxSec], (float) compSum2.x) ;
    atomicAdd (&pfP2[iIdxSec + 1], (float) compSum2.y) ;

    atomicAdd (&pfP3[iIdxSec], (float) compSum3.x) ;
    atomicAdd (&pfP3[iIdxSec + 1], (float) compSum3.y) ;

}
#endif 


__global__ void scaleAndSub (float  *pCompData1,
                             float 	*pCompData2,
                             float  fScale,
                             int    iNx,
                             int    iNy,
                             int    iNz)
{
    int           iX ;
    int           iY ;
    int           iZ ;
    unsigned int  iIdx ;

    iX = threadIdx.x ;
    iY = blockIdx.x ;
    iZ = blockIdx.y ;

    if ((iX < iNx) && (iY < iNy) && (iZ < iNz))
    {
        iIdx = (((iZ * iNy) + iY) * iNx) + iX ;
        pCompData1[iIdx] = pCompData2[iIdx] - (pCompData1[iIdx] * fScale) ;
    }
}

/**
 * @brief 	This function is called from the surfacetractioncowling of fortran code in case there 
 *          was insufficient memory and couldnt process all the things at once.
 */

extern "C" void cusurfacetraction_ (double 	dModulus, 
                                    double 	dAlpha, 
                                    double 	dLambda, 
                                    double 	dMu, 
                                    double 	dGravity, 
                                    float 	*fData1,
                                    float 	*fData2,
                                    float 	*fData3,
                                    double 	dDx1, 
                                    double 	dDx2, 
                                    double 	dDx3,
                                    float 	*fOutData1,
                                    float 	*fOutData2,
                                    float 	*fOutData3,
                                    int 	iSx1,
                                    int 	iSx2,
                                    int 	iSx3)
{
    cudaError_t cuError ;
    int 	    iSize ;
    float       *pfP1 ;
    float       *pfP2 ;
    float       *pfP3 ;
    float       fScale ;

    dim3 dimGrid(iSx2, iSx3, 1) ;
    dim3 dimBlock((iSx1/2 + 1), 1, 1) ;

    dim3 dimGrid1(iSx2, 1, 1) ;
    dim3 dimBlock1((iSx1 + 2), 1, 1) ;

    iSize = (sizeof (float) * (iSx1 + 2) * (iSx2) * (iSx3)) ;

    cuError = cudaMalloc((void**)&gpCompData1, iSize) ;
    if (cudaSuccess != cuError)
    {
        printf ("Memory allocation failure 1\n") ;
        return ;
    }

    cuError = cudaMalloc((void**)&gpCompData2, iSize) ;
    if (cudaSuccess != cuError)
    {
        cudaFree (gpCompData1) ;
        printf ("Memory allocation failure 2\n") ;
        return ;
    }
    cuError = cudaMalloc((void**)&gpCompData3, iSize) ;
    if (cudaSuccess != cuError)
    {
        cudaFree (gpCompData1) ;
        cudaFree (gpCompData2) ;
        printf ("Memory allocation failure 3\n") ;
        return ;
    }

    cuError = cudaMalloc((void**)&pfP1, ((iSx1+2) * iSx2 * sizeof (float))) ;
    if (cudaSuccess != cuError)
    {
        printf ("Memory allocation failure 4\n") ;
        goto FREE_DATA_SURFACE ;
    }
        
    cuError = cudaMalloc((void**)&pfP2, ((iSx1+2) * iSx2 * sizeof (float))) ;
    if (cudaSuccess != cuError)
    {
        cudaFree (pfP1) ;
        printf ("Memory allocation failure 5\n") ;
        goto FREE_DATA_SURFACE ;
    }
        
    cuError = cudaMalloc((void**)&pfP3, ((iSx1+2) * iSx2 * sizeof (float))) ;
    if (cudaSuccess != cuError)
    {
        cudaFree (pfP1) ;
        cudaFree (pfP2) ;
        printf ("Memory allocation failure 6\n") ;
        goto FREE_DATA_SURFACE ;
    }
        
    /* To make sure all are 0 since we add in the kernel */
    cudaMemset (pfP1, 0,((iSx1+2) * iSx2 * sizeof (float))) ;
    CHECK_CUDA_ERROR ("SurfaceTraction: cudaMemset1", FREE_DATA_SURFACE)

    cudaMemset (pfP2, 0,((iSx1+2) * iSx2 * sizeof (float))) ;
    CHECK_CUDA_ERROR ("SurfaceTraction: cudaMemset2", FREE_DATA_SURFACE)

    cudaMemset (pfP3, 0,((iSx1+2) * iSx2 * sizeof (float))) ;
    CHECK_CUDA_ERROR ("SurfaceTraction: cudaMemset3", FREE_DATA_SURFACE)

    cuError = cudaMemcpy (gpCompData1, fData1, iSize, cudaMemcpyHostToDevice) ;
    CHECK_CUDA_ERROR ("SurfaceTraction: cudaMemcpy1", FREE_DATA_SURFACE)

    cuError = cudaMemcpy (gpCompData2, fData2, iSize, cudaMemcpyHostToDevice) ;
    CHECK_CUDA_ERROR ("SurfaceTraction: cudaMemcpy2", FREE_DATA_SURFACE)

    cuError = cudaMemcpy (gpCompData3, fData3, iSize, cudaMemcpyHostToDevice) ;
    CHECK_CUDA_ERROR ("SurfaceTraction: cudaMemcpy3", FREE_DATA_SURFACE)

    fScale = (float)(1.0 / (iSx3 * dDx3)) ;
        
    cuSurfaceKernel <<<dimGrid, dimBlock>>> (dModulus, dLambda, dMu, dGravity, dDx1, dDx2, dDx3, 
                                             pfP1, pfP2, pfP3, (iSx1/2 + 1), iSx2, iSx3, 
                                             (float *)gpCompData1, (float *)gpCompData2, 
                                             (float *)gpCompData3) ;
    cuError = cudaGetLastError () ;
    CHECK_CUDA_ERROR ("Surfacetraction : cuSurfaceKernel ", FREE_DATA_SURFACE)
        
    scaling <float> <<<dimGrid1, dimBlock1>>> (pfP1, fScale, (iSx1 + 2), iSx2, 2) ;
    cuError = cudaGetLastError () ;
    CHECK_CUDA_ERROR ("SurfaceTraction : scalingKernel 1", FREE_DATA_SURFACE)

    scaling <float> <<<dimGrid1, dimBlock1>>> (pfP2, fScale, (iSx1 + 2), iSx2, 2) ;
    cuError = cudaGetLastError () ;
    CHECK_CUDA_ERROR ("SurfaceTraction : scalingKernel 2", FREE_DATA_SURFACE)

    scaling <float> <<<dimGrid1, dimBlock1>>> (pfP3, fScale, (iSx1 + 2), iSx2, 2) ;
    cuError = cudaGetLastError () ;
    CHECK_CUDA_ERROR ("SurfaceTraction : scalingKernel 3", FREE_DATA_SURFACE)

        
    cuError = cudaMemcpy (fOutData1, pfP1, ((iSx1+2) * iSx2 * sizeof (float)), 
                          cudaMemcpyDeviceToHost) ;
    CHECK_CUDA_ERROR ("SurfaceTraction: cudaMemcpy4", FREE_DATA_SURFACE)	

    cuError = cudaMemcpy (fOutData2, pfP2, ((iSx1+2) * iSx2 * sizeof (float)), 
                          cudaMemcpyDeviceToHost) ;
    CHECK_CUDA_ERROR ("SurfaceTraction: cudaMemcpy5", FREE_DATA_SURFACE)

    cuError = cudaMemcpy (fOutData3, pfP3, ((iSx1+2) * iSx2 * sizeof (float)), 
                          cudaMemcpyDeviceToHost) ;
    CHECK_CUDA_ERROR ("SurfaceTraction: cudaMemcpy6", FREE_DATA_SURFACE)


FREE_DATA_SURFACE :
    cudaFree (gpCompData1) ;
    cudaFree (gpCompData2) ;
    cudaFree (gpCompData3) ;
    cudaFree (pfP1) ;
    cudaFree (pfP2) ;
    cudaFree (pfP3) ;

    return ;
}
	
	
/**
 * brief :	This is kernel  funtion for surfacetraction calculations.
 *
 */

__global__ void cuSurfaceKernel (double  dModulus,
                                 double  dLambda,
                                 double  dMu,
                                 double  dGravity,
                                 double  dDx1,
                                 double  dDx2,
                                 double  dDx3,
                                 float   pfP1[],
                                 float   pfP2[],
                                 float   pfP3[],
                                 int     iSx1,
                                 int     iSx2,
                                 int     iSx3, 
                                 float 	 *fData1, 
                                 float 	 *fData2, 
                                 float 	 *fData3)
{
    unsigned int iIdx1 ;
    unsigned int iIdx2 ;
    unsigned int iIdx3 ;
    unsigned int iIdxSec ;
    unsigned int iIdx ;
        
    double       dK1 ;
    double       dK2 ;
    double       dK3 ;

    float2 		 *pC1 ;
    float2 		 *pC2 ;
    float2 		 *pC3 ;
        
    double2 	 cuC1 ;
    double2 	 cuC2 ;
    double2 	 cuC3 ;

    double2 	 compSum1 ;
    double2 	 compSum2 ;
    double2 	 compSum3 ;

    iIdx1 = threadIdx.x ;
    iIdx2 = blockIdx.x ;
    iIdx3 = blockIdx.y ;

    if ((iIdx1 < iSx1) && (iIdx2 < iSx2) && (iIdx3 < iSx3))
    {
                    
        cuWaveNumber (iIdx1, iIdx2, iIdx3, (2 * iSx1) - 2, iSx2, iSx3, dDx1, dDx2, dDx3, 
                      &dK1, &dK2, &dK3) ;

        iIdx = (((iIdx3 * iSx2) + iIdx2) * 2 * iSx1) + 2 * iIdx1 ;
                
        //printf ("iIdx = %u\n", iIdx) ;	
        pC1 = (float2 *) &fData1[iIdx] ;
        pC2 = (float2 *) &fData2[iIdx] ;
        pC3 = (float2 *) &fData3[iIdx] ;

        cuC1 = make_comp_double ((double) ((*pC1).x), (double) ((*pC1).y)) ;
        cuC2 = make_comp_double ((double) ((*pC2).x), (double) ((*pC2).y)) ;
        cuC3 = make_comp_double ((double) ((*pC3).x), (double) ((*pC3).y)) ;
            
        /*cuC1 = make_comp_double ((double) fData1[iIdx], (double) fData1[iIdx + 1]) ;
        cuC2 = make_comp_double ((double) fData2[iIdx], (double) fData2[iIdx + 1]) ;
        cuC3 = make_comp_double ((double) fData3[iIdx], (double) fData3[iIdx + 1]) ;	
        */
        compSum1 = comp_double_mult (comp_double_mult (make_comp_double (0, PI2), 
                   make_comp_double (dMu, 0)), comp_double_add (comp_double_mult (
                   make_comp_double (dK3, 0), cuC1), comp_double_mult (make_comp_double (dK1, 0), 
                   cuC3))) ;
          
        compSum2 = comp_double_mult (comp_double_mult (make_comp_double (0, PI2), 
                   make_comp_double (dMu, 0)), comp_double_add (comp_double_mult (
                   make_comp_double (dK3, 0), cuC2), comp_double_mult (make_comp_double (dK2, 0), 
                   cuC3))) ; 

        compSum3 = comp_double_sub (comp_double_mult (make_comp_double (0, PI2), 
                   comp_double_add (comp_double_mult (comp_double_mult (
                   make_comp_double (dModulus, 0), make_comp_double(dK3, 0)), cuC3), 
                   comp_double_mult (make_comp_double (dLambda, 0), 
                   comp_double_add (comp_double_mult (make_comp_double (dK1, 0), cuC1), 
                   comp_double_mult (make_comp_double (dK2, 0), cuC2))))), 
                   comp_double_mult (make_comp_double (dGravity, 0), cuC3)) ;

        iIdxSec = (iIdx2 * 2 * iSx1) + 2 * iIdx1 ;		

        atomicAdd (&pfP1[iIdxSec], (float) compSum1.x) ;
        atomicAdd (&pfP1[iIdxSec + 1], (float) compSum1.y) ;
                
        atomicAdd (&pfP2[iIdxSec], (float) compSum2.x) ;
        atomicAdd (&pfP2[iIdxSec + 1], (float) compSum2.y) ;
                
        atomicAdd (&pfP3[iIdxSec], (float) compSum3.x) ;
        atomicAdd (&pfP3[iIdxSec + 1], (float) compSum3.y) ;
        /*	
            pfP1[iIdxSec] += (float) compSum1.x ; 
            pfP1[iIdxSec + 1] += (float) compSum1.y ; 
            
            pfP2[iIdxSec] += (float) compSum2.x ;
                    pfP2[iIdxSec + 1] += (float) compSum2.y ;

            pfP3[iIdxSec] += (float) compSum3.x ;
                    pfP3[iIdxSec + 1] += (float) compSum3.y ;
     */
    }	
}

static void cuAllocate (int iSize)
{
    cudaError_t     cuError ;

    cuError = cudaMalloc((void**)&pfB1, iSize) ;
    if (cudaSuccess != cuError)
    {
        printf ("cuAllocate : Memory allocation failure 1\n") ;
        return ;
    }

    cuError = cudaMalloc((void**)&pfB2, iSize) ;
    if (cudaSuccess != cuError)
    {
        cudaFree (pfB1) ;
        printf ("cuAllocate : Memory allocation failure 2\n") ;
        return ;
    }

    cuError = cudaMalloc((void**)&pfB3, iSize) ;
    if (cudaSuccess != cuError)
    {
        cudaFree (pfB1) ;
        cudaFree (pfB2) ;
        printf ("cuAllocate : Memory allocation failure 3\n") ;
        return ;
    }
}


extern "C" void cucerruti_ (float   *p1,
                            float   *p2,
                            float   *p3,
                            float   *fData1,
                            float   *fData2,
                            float   *fData3,
                            double  dLambda,
                            double  dMu,
                            double  dGamma,
                            double  dDx1,
                            double  dDx2,
                            double  dDx3,
                            int     iSx1,
                            int     iSx2,
                            int     iSx3)
{
    double      dAlpha ;
    cudaError_t cuError ;
    int         iSize ;


    iSize = (sizeof (float) * (iSx1 + 2) * (iSx2) * (iSx3)) ;

    cuError = cudaMalloc((void**)&gpCompData1, iSize) ;
    if (cudaSuccess != cuError)
    {
        printf ("Memory allocation failure 1\n") ;
        return  ;
    }

    cuError = cudaMalloc((void**)&gpCompData2, iSize) ;
    if (cudaSuccess != cuError)
    {
        cudaFree (gpCompData1) ;
        printf ("Memory allocation failure 2\n") ;
        return ;
    }

    cuError = cudaMalloc((void**)&gpCompData3, iSize) ;
    if (cudaSuccess != cuError)
    {
        cudaFree (gpCompData1) ;
        cudaFree (gpCompData2) ;
        printf ("Memory allocation failure 3\n") ;
        return ;
    }

    cuError = cudaMemcpy (gpCompData1, fData1, iSize, cudaMemcpyHostToDevice) ;
    CHECK_CUDA_ERROR ("cucerruti_: cudaMemcpy1", FREE_DATA_CERRUTI)

    cuError = cudaMemcpy (gpCompData2, fData2, iSize, cudaMemcpyHostToDevice) ;
    CHECK_CUDA_ERROR ("cucerruti_: cudaMemcpy2", FREE_DATA_CERRUTI)

    cuError = cudaMemcpy (gpCompData3, fData3, iSize, cudaMemcpyHostToDevice) ;
    CHECK_CUDA_ERROR ("cucerruti_: cudaMemcpy3", FREE_DATA_CERRUTI)

    dAlpha = (dLambda + dMu) / (dLambda + (2 * dMu)) ;


    cerruti (p1, p2, p3, iSx1, iSx2, iSx3, dDx1, dDx2, dDx3, dMu,
             dLambda, dAlpha, dGamma, (float *)gpCompData1, 
             (float *)gpCompData2, (float *)gpCompData3) ;

    cuError = cudaMemcpy (fData1, gpCompData1, iSize, cudaMemcpyDeviceToHost) ;
    CHECK_CUDA_ERROR ("cucerruti_: cudaMemcpy1", FREE_DATA_CERRUTI)

    cuError = cudaMemcpy (fData2, gpCompData2, iSize, cudaMemcpyDeviceToHost) ;
    CHECK_CUDA_ERROR ("cucerruti_: cudaMemcpy2", FREE_DATA_CERRUTI)

    cuError = cudaMemcpy (fData3, gpCompData3, iSize, cudaMemcpyDeviceToHost) ;
    CHECK_CUDA_ERROR ("cucerruti_: cudaMemcpy3", FREE_DATA_CERRUTI)

FREE_DATA_CERRUTI:
    cudaFree (gpCompData1) ;
    cudaFree (gpCompData2) ;
    cudaFree (gpCompData3) ;

}

static void cerruti (float   *p1,
                     float   *p2,
                     float   *p3,
                     int     iSx1,
                     int     iSx2,
                     int     iSx3,
                     double  dDx1,
                     double  dDx2,
                     double  dDx3,
                     double  dMu,
                     double  dLambda,
                     double  dAlpha,
                     double  dGamma,
                     float   *pfU1,
                     float   *pfU2,
                     float   *pfU3)
{
    int          iIn1 ;
    int          iIn2 ;
    dim3         dimGrid(1, 1, 1) ;
    dim3         dimBlock (iSx3, 1, 1) ;
    int          iIdxSec ;

    cufftResult  cuRet ;
    double       dK1 ;
    double       dK2 ;
    double       dK3 ;

    double2      cuC1 ;
    double2      cuC2 ;
    double2      cuC3 ;


    cuAllocate (iSx3 * sizeof (cufftComplex)) ;

    cuRet = cufftPlan1d(&ghPlanOne, iSx3, CUFFT_C2C, 1) ;
    if (CUFFT_SUCCESS != cuRet)
    {
        printf ("Plan creation is failed\n") ;
        return ;
    }
            
    for (iIn2 = 0; iIn2 < iSx2 ; iIn2++)
    {
        for (iIn1 = 0; iIn1 < ((iSx1/2)+1) ; iIn1++)
        {
            cuWaveNumber (iIn1, iIn2, 1, iSx1, iSx2, 1, dDx1, dDx2, 1, &dK1, &dK2, &dK3) ;

            iIdxSec = (iIn2 * ((iSx1/2)+1) * 2) + 2 * iIn1 ;

            cuC1 = make_comp_double ((double)p1[iIdxSec], (double)p1[iIdxSec + 1]) ;
            cuC2 = make_comp_double ((double)p2[iIdxSec], (double)p2[iIdxSec + 1]) ;
            cuC3 = make_comp_double ((double)p3[iIdxSec], (double)p3[iIdxSec + 1]) ;


            cerrutiKernel <<<dimGrid, dimBlock>>> (cuC1, cuC2, cuC3, (float2 *)pfB1, (float2 *)pfB2,
                                                   (float2 *)pfB3, iSx3, dK1, dK2, dDx3, dMu, 
                                                   dAlpha, dGamma) ;

            if (cudaDeviceSynchronize() != cudaSuccess)
            {
                fprintf(stderr, "Cuda error: Failed to synchronize\n") ;
                return ;
            }
            OneDfft (ghPlanOne, pfB1, pfB2, pfB3, iSx3, dDx3) ;

            cuAddKernel <<<dimGrid, dimBlock>>> (pfU1, pfU2, pfU3, (float2 *)pfB1, (float2 *)pfB2, 
                                                 (float2 *)pfB3, iIn1, iIn2, iSx1, iSx2, iSx3) ;

            if (cudaDeviceSynchronize() != cudaSuccess)
            {
                fprintf(stderr, "Cuda error: Failed to synchronize\n") ;
                return ;
            }
        }
    }
    cufftDestroy (ghPlanOne) ;
    cudaFree (pfB1) ;
    cudaFree (pfB2) ;
    cudaFree (pfB3) ;

}

__global__ void cuAddKernel (float   *pfU1,
                             float   *pfU2,
                             float   *pfU3,
                             float2  *pfB1,
                             float2  *pfB2,
                             float2  *pfB3,
                             int     iIn1,
                             int     iIn2,
                             int     iSx1,
                             int     iSx2,
                             int     iSx3)
{
    int     iIdx ;
    int     iIdx1 ;

    iIdx1 = threadIdx.x ;

    if (iIdx1 < iSx3)
    {
        iIdx = ((iIdx1 * iSx2) + iIn2) * (iSx1 + 2) + 2 * iIn1 ;

        pfU1[iIdx] += pfB1[iIdx1].x ;
        pfU1[iIdx + 1] += pfB1[iIdx1].y ;

        pfU2[iIdx] += pfB2[iIdx1].x ;
        pfU2[iIdx + 1] += pfB2[iIdx1].y ;

        pfU3[iIdx] += pfB3[iIdx1].x ;
        pfU3[iIdx + 1] += pfB3[iIdx1].y ;
    }
}

__global__ void cerrutiKernel (double2  cuC1,
                               double2  cuC2,
                               double2  cuC3,
                               float2   *pfB1,
                               float2   *pfB2,
                               float2   *pfB3,
                               int      iSx3,
                               double   dK1,
                               double   dK2,
                               double   dDx3,
                               double   dMu,
                               double   dAlpha,
                               double   dGamma)
{
    unsigned int iIdx ;
    double       dX3 ;
    float2       dB1 ;
    float2       dB2 ;
    float2       dB3 ;

    iIdx = threadIdx.x ;

    if (iSx3 > iIdx)
    {
        if (iIdx < iSx3/2)
        {
            dX3 = iIdx * dDx3 ;
        }
        else
        {
            dX3 = ((iSx3 - iIdx) * dDx3) ;
        }

        cerrutiSolution (dMu, cuC1, cuC2, cuC3, dAlpha, dGamma,  &dB1,
        &dB2, &dB3, dK1, dK2, dX3, ((double) (iSx3/2) * dDx3), iIdx) ;

        pfB1[iIdx].x = dB1.x ;
        pfB1[iIdx].y = dB1.y ;

        pfB2[iIdx].x = dB2.x ;
        pfB2[iIdx].y = dB2.y ;

        pfB3[iIdx].x = dB3.x ;
        pfB3[iIdx].y = dB3.y ;
    }
}

static void OneDfft (cufftHandle  ghPlanOne,
                     float        *pB1,
                     float        *pB2,
                     float        *pB3,
                     int          iSx3,
                     double       dDx3)
{
    cufftResult     cuRet ;
    dim3 dimGrid(1,1,1) ;
    dim3 dimBlock(iSx3,1,1) ;

    cuRet = cufftExecC2C(ghPlanOne, (cufftComplex *)pB1, (cufftComplex *)pB1, CUFFT_FORWARD) ;
    if (CUFFT_SUCCESS != cuRet)
    {
        printf ("OneDfft : Failed 1\n") ;
        return ;
    }

    cuRet = cufftExecC2C(ghPlanOne, (cufftComplex *) pB2, (cufftComplex *) pB2, CUFFT_FORWARD) ;
    if (CUFFT_SUCCESS != cuRet)
    {
        printf ("OneDfft : Failed 2\n") ;
        return ;
    }
    cuRet = cufftExecC2C(ghPlanOne, (cufftComplex *)pB3, (cufftComplex *)pB3, CUFFT_FORWARD) ;
    if (CUFFT_SUCCESS != cuRet)
    {
        printf ("OneDfft : Failed 3\n") ;
    }
    scale1D <<<dimGrid, dimBlock>>> (pB1, pB2, pB3, (float)dDx3) ;

}

__global__ void scale1D (float *pB1,
                         float *pB2,
                         float *pB3,
                         float fScale)
{
    int iIdx = threadIdx.x ;

    if (iIdx < blockDim.x)
    {
        iIdx = 2 * iIdx ;

        pB1[iIdx] = pB1[iIdx] * fScale ;
        pB1[iIdx + 1] = pB1[iIdx + 1] * fScale ;
        pB2[iIdx] = pB2[iIdx] * fScale ;
        pB2[iIdx + 1] = pB2[iIdx + 1] * fScale ;
        pB3[iIdx] = pB3[iIdx] * fScale ;
        pB3[iIdx + 1] = pB3[iIdx + 1] * fScale ;
    }

}

__device__ void cerrutiSolution (double   dMu,
                                 double2  cuC1,
                                 double2  cuC2,
                                 double2  cuC3,
                                 double   dAlpha,
                                 double   dGamma,
                                 float2   *ppB1,
                                 float2   *ppB2,
                                 float2   *ppB3,
                                 double   dK1,
                                 double   dK2,
                                 double   dX3,
                                 double   dL,
				                 int 	  iIdx)
{
    double   dBeta ;
    double   dDepthDecay ;
    double   dH ;
    double   dTemp ;

    double2  i ;
    double2  cuB1 ;
    double2  cuB2 ;
    double2  cuB3 ;
    double2  cuTemp ;
    double2  cuV1 ;
    double2  cuV2 ;
    double2  cuV3 ;
            
    if ((0 == dK1) && (0 == dK2))
    {
        *ppB1 = make_comp_float ((cuC1.x / dMu) * (dX3 - dL), 0) ;
        *ppB2 = make_comp_float ((cuC2.x / dMu) * (dX3 - dL), 0) ;
        *ppB3 = make_comp_float (((cuC3.x / dMu) * (dX3 - dL) * (1.0 - dAlpha)) /
                                 (1.0 + 2.0 * (dL * dAlpha * dGamma * (1.0 - dAlpha))), 0) ;
    }
    else
    {
        dBeta = PI2 * sqrt ((dK1 * dK1) + (dK2 * dK2)) ;
        dDepthDecay = exp (-dBeta * abs (dX3)) ;
        dH = dGamma / dBeta ;
        i = make_comp_double (0, PI2) ;


        dTemp = (1.0 / (2.0 * dMu * (dBeta * dBeta * dBeta))) * dDepthDecay ;
        cuB1 = comp_double_mult (make_comp_double(dTemp, 0),
               make_comp_double ((double)cuC1.x, (double)cuC1.y)) ;
        cuB2 = comp_double_mult (make_comp_double(dTemp, 0),
               make_comp_double ((double)cuC2.x, (double)cuC2.y)) ;

        cuB3 = comp_double_div (comp_double_mult (make_comp_double (dTemp, 0),
               make_comp_double ((double) cuC3.x, (double) cuC3.y)),
               make_comp_double (1.0+dH, 0)) ;


        cuTemp = comp_double_mult (comp_double_mult (i, cuB3), make_comp_double (dBeta * (1.0 -
                                  (1.0 / dAlpha) + (dBeta * dX3)), 0)) ;

        cuV1 = comp_double_mult (cuTemp, make_comp_double (dK1, 0)) ;
        cuV2 = comp_double_mult (cuTemp, make_comp_double (dK2, 0)) ;
        cuV3 = comp_double_mult (make_comp_double ((dBeta * dBeta) *
                                 ((1.0 / dAlpha) + (dBeta * dX3)), 0), cuB3) ;
        cuV3 = comp_double_mult (cuV3, make_comp_double(-1.0,0)) ;
		

        dTemp = (PI2 * PI2) * (2.0  - (1.0/dAlpha) + (dBeta * dX3)) / (1+dH) ;
                
		cuV1 = comp_double_add (cuV1, comp_double_mult (cuB1, make_comp_double
                               (-2.0 * (dBeta * dBeta) + ((dK1 * dK1) * dTemp), 0))) ;
        cuV2 = comp_double_add (cuV2, comp_double_mult (cuB1, make_comp_double
                                (dK1 * dK2 * dTemp, 0))) ;
        cuV3 = comp_double_add (cuV3, comp_double_mult (comp_double_mult (cuB1, i), make_comp_double
                                ((dK1 * dBeta) * ((1.0/dAlpha) - 1.0 + (dBeta * dX3))/(1+dH), 0))) ;


        cuTemp = comp_double_add (cuV1, comp_double_mult (cuB2, make_comp_double
                                 ((dK1 * dK2 * dTemp), 0))) ;
        (*ppB1).x = (float) cuTemp.x ;
        (*ppB1).y = (float) cuTemp.y ;

        cuTemp = comp_double_add (cuV2, comp_double_mult (cuB2, make_comp_double
                                 ((-2.0 * (dBeta * dBeta) + (dK2 * dK2) * dTemp), 0))) ;
        (*ppB2).x = (float) cuTemp.x ;
        (*ppB2).y = (float) cuTemp.y ;
		

        cuTemp = comp_double_add (cuV3, comp_double_mult
                 (comp_double_mult (cuB2, i), make_comp_double
                 (dBeta * dK2 * ((1.0/dAlpha) - 1.0 + (dBeta * dX3)) / (1+dH), 0))) ;
        (*ppB3).x = (float) cuTemp.x ;
        (*ppB3).y = (float) cuTemp.y ;
 
    }
}

#endif

/* EOF */
