!-----------------------------------------------------------------------
! Copyright 2007, 2008, 2009 Sylvain Barbot
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
!-----------------------------------------------------------------------

#include "include.f90"

MODULE fourier

#ifdef IMKL_FFT
  USE MKL_DFTI
#endif

  IMPLICIT NONE

  PUBLIC

#ifdef FFTW3
  INCLUDE 'fftw3.f'
#endif

  INTEGER, PARAMETER :: FFT_FORWARD=-1,FFT_INVERSE=1

CONTAINS

  !---------------------------------------------------------------------
  ! subroutine wavenumbers 
  ! computes the values of the wavenumbers
  ! in the sequential order required when using subroutine FOURT
  ! to perform forward and backward inverse transforms.
  !
  ! INPUT
  ! i1 i3     running index in the discrete Fourier domain array
  ! sx1 sx3  number of elements in the 2 directions
  ! dx1 dx3  sampling interval in the 2 directions
  !
  ! OUTPUT
  ! k1 k3     wavenumbers in the 2 direction
  !
  ! sylvain barbot (04-14-07) - original form
  !---------------------------------------------------------------------
  SUBROUTINE wavenumbers(i1,i2,i3,sx1,sx2,sx3,dx1,dx2,dx3,k1,k2,k3)
    INTEGER, INTENT(IN) :: i1, i2, i3, sx1, sx2, sx3
    REAL*8, INTENT(IN) :: dx1, dx2, dx3
    REAL*8, INTENT(OUT) :: k1, k2, k3
    
    IF (i3 < sx3/2+1) THEN
       k3= (DBLE(i3)-1._8)/(sx3*dx3)
    ELSE
       k3=-(DBLE(sx3-i3)+1._8)/(sx3*dx3)
    END IF
    IF (i2 < sx2/2+1) THEN
       k2= (DBLE(i2)-1._8)/(sx2*dx2)
    ELSE
       k2=-(DBLE(sx2-i2)+1._8)/(sx2*dx2)
    END IF
    k1=(DBLE(i1)-1._8)/(sx1*dx1)
    
  END SUBROUTINE wavenumbers

  SUBROUTINE wavenumber1(i1,sx1,dx1,k1)
    INTEGER, INTENT(IN) :: i1,sx1
    REAL*8, INTENT(IN) :: dx1
    REAL*8, INTENT(OUT) :: k1

    k1=(DBLE(i1)-1._8)/(sx1*dx1)
  END SUBROUTINE wavenumber1

  SUBROUTINE wavenumber2(i2,sx2,dx2,k2)
    INTEGER, INTENT(IN) :: i2,sx2
    REAL*8, INTENT(IN) :: dx2
    REAL*8, INTENT(OUT) :: k2
    
    IF (i2 < sx2/2+1) THEN
       k2= (DBLE(i2)-1._8)/(sx2*dx2)
    ELSE
       k2=-(DBLE(sx2-i2)+1._8)/(sx2*dx2)
    END IF
  END SUBROUTINE wavenumber2

  SUBROUTINE wavenumber3(i3,sx3,dx3,k3)
    INTEGER, INTENT(IN) :: i3,sx3
    REAL*8, INTENT(IN) :: dx3
    REAL*8, INTENT(OUT) :: k3
    
    IF (i3 < sx3/2+1) THEN
       k3= (DBLE(i3)-1._8)/(sx3*dx3)
    ELSE
       k3=-(DBLE(sx3-i3)+1._8)/(sx3*dx3)
    END IF
  END SUBROUTINE wavenumber3

  !---------------------------------------------------------------------
  ! subroutine FFTshift_TF applies the transfer function 
  ! in the Fourier domain corresponding to shifting the space 
  ! domain array by sx1*dx1/2 in the 1-direction and sx3*dx3/2 
  ! in the 3-direction.
  !
  ! fftshift_tf follows the data storage convention in
  ! agreement with DFT subroutine FOURT
  !
  ! sylvain barbot (05-01-07) - original form
  !---------------------------------------------------------------------
  SUBROUTINE fftshift_tf(spec)
    REAL*4, INTENT(INOUT), DIMENSION(:,:,:) :: spec
    
    INTEGER :: sx1, sx2, sx3, i1, i2, i3
    REAL*4 :: exp1, exp2, exp3
    
    sx1=SIZE(spec, 1)-2
    sx2=SIZE(spec, 2)
    sx3=SIZE(spec, 3)
    
    DO i3=1,sx3
       IF (i3 < sx3/2+1) THEN
          exp3=-(DBLE(i3)-1._8)
       ELSE
          exp3= (DBLE(sx3-i3)+1._8)
       END IF
       DO i2=1,sx2
          IF (i2 < sx2/2+1) THEN
             exp2=-(DBLE(i2)-1._8)
          ELSE
             exp2= (DBLE(sx2-i2)+1._8)
          END IF
          DO i1=1,sx1/2+1
             exp1=(DBLE(i1)-1._8)
             spec(2*i1-1:2*i1,i2,i3) = &
                  spec(2*i1-1:2*i1,i2,i3)*((-1._4)**(exp1+exp2+exp3))
          END DO
       END DO
    END DO
  END SUBROUTINE fftshift_tf

  !----------------------------------------------------------------------
  ! subroutine FFT3 performs normalized forward and
  ! inverse fourier transforms of real 3d data
  !
  ! USES
  ! ctfft (Brenner, 1968) by default
  ! fftw3 (Frigo & Jonhson) with preproc FFTW3 flag
  ! scfft (SGI library) with preproc SGI_FFT flag
  !
  ! for real array the fourier transform returns a sx1/2+1 complex array
  ! and the enough space must be reserved
  !----------------------------------------------------------------------
#ifdef FFTW3
  !--------------------------------------------------------
  ! implementation of FFTW3
  ! must be linked with -lfftw3f (single-threaded version)
  !
  ! sylvain barbot (09-28-08) - original form
  !--------------------------------------------------------
  SUBROUTINE fft3(data,sx1,sx2,sx3,dx1,dx2,dx3,direction)
    INTEGER, INTENT(IN) :: sx1,sx2,sx3,direction
    REAL*4, DIMENSION(sx1+2,sx2,sx3), INTENT(INOUT) :: data
    REAL*8, INTENT(IN) :: dx1,dx2,dx3

    INTEGER*8 :: plan

    IF (FFT_FORWARD == direction) THEN
      CALL sfftw_plan_dft_r2c_3d(plan,sx1,sx2,sx3, &
           data(1,1,1),data(1,1,1),FFTW_ESTIMATE)
    ELSE
      CALL sfftw_plan_dft_c2r_3d(plan,sx1,sx2,sx3, &
           data(1,1,1),data(1,1,1),FFTW_ESTIMATE)
    END IF

    CALL sfftw_execute(plan)
    CALL sfftw_destroy_plan(plan)

   IF (FFT_INVERSE == direction) THEN
     data=data/(sx1*dx1*sx2*dx2*sx3*dx3)
   ELSE
     data=data*(dx1*dx2*dx3)
   END IF

  END SUBROUTINE fft3
#else
#ifdef SGI_FFT
  !--------------------------------------------------------------------
  ! implementation of SGI SCFFT
  ! must be linked with -L/usr/lib -lscs or -L/usr/lib -lscs_mp for
  ! multithread version expect up x8 performance increase compared to
  ! ctfft implementation. check out the SGI documentation at:
  !
  ! http://techpubs.sgi.com/library/tpl/cgi-bin/getdoc.cgi?coll=linux&
  !      db=man&fname=/usr/share/catman/man3/ccfft.3s.html&srch=ccfft
  !
  ! sylvain barbot (09-28-08) - original form
  !--------------------------------------------------------------------
  SUBROUTINE fft3(data,sx1,sx2,sx3,dx1,dx2,dx3,direction)
    INTEGER, INTENT(IN) :: direction,sx1,sx2,sx3
    REAL*4, DIMENSION(:,:,:), INTENT(INOUT) :: data
    REAL*8, INTENT(IN) :: dx1,dx2,dx3

    INTEGER, PARAMETER :: NF=256, NFR=256

    REAL*4, DIMENSION(sx1+NFR+(2*sx2+NF)+(2*sx3+NF)) :: table
    REAL*4, DIMENSION(sx1+4*sx3) :: work
    INTEGER, DIMENSION(2) :: isys
    REAL*4 :: scale

    isys(1)=1

    IF (FFT_FORWARD == direction) THEN
      scale=dx1*dx2*dx3
      ! initialize the sin/cos table
      CALL SCFFT3D(+0,sx1,sx2,sx3,scale,data(1,1,1),sx1+2,sx2, &
                   data(1,1,1),sx1/2+1,sx2,table,work,isys)
      CALL SCFFT3D(-1,sx1,sx2,sx3,scale,data(1,1,1),sx1+2,sx2, &
                   data(1,1,1),sx1/2+1,sx2,table,work,isys)
    ELSE
      scale=1._4/(sx1*dx1*sx2*dx2*sx3*dx3)
      ! initialize the sin/cos table
      CALL CSFFT3D(+0,sx1,sx2,sx3,scale,data(1,1,1),sx1/2+1,sx2, &
                   data(1,1,1),sx1+2,sx2,table,work,isys)
      CALL CSFFT3D(+1,sx1,sx2,sx3,scale,data(1,1,1),sx1/2+1,sx2, &
                   data(1,1,1),sx1+2,sx2,table,work,isys)
    END IF

  END SUBROUTINE fft3
#else
#ifdef IMKL_FFT
  !-------------------------------------------------------------------------
  ! implementation IMKL_FFT (Intel Math Kernel Library)
  ! for information and example calculations with the
  ! mkl FFT, see:
  !
  ! http://www.intel.com/software/products/mkl/docs/webhelp/appendices/ ...
  !                      mkl_appC_DFT.html#appC-exC-25
  !
  ! and a thread (Fortran 3-D FFT real-to-complex ...)
  ! on the intel forum
  !
  ! http://software.intel.com/en-us/forums/intel-math-kernel-library/
  !
  ! sylvain barbot (04-30-10) - original form
  !-------------------------------------------------------------------------
  SUBROUTINE fft3(data,sx1,sx2,sx3,dx1,dx2,dx3,direction)
    REAL*4, DIMENSION(0:*), INTENT(INOUT) :: data
    REAL*8, INTENT(IN) :: dx1,dx2,dx3
    INTEGER, INTENT(IN) :: sx1,sx2,sx3,direction

    INTEGER :: iret,size(3),rstrides(4),cstrides(4)
    TYPE(DFTI_DESCRIPTOR), POINTER :: desc
    REAL*4 :: scale

    rstrides=(/ 0,1,(sx1/2+1)*2,(sx1/2+1)*2*sx2 /)
    cstrides=(/ 0,1,sx1/2+1,(sx1/2+1)*sx2 /)
    size=(/ sx1,sx2,sx3 /)

    iret=DftiCreateDescriptor(desc,DFTI_SINGLE,DFTI_REAL,3,size)
    iret=DftiSetValue(desc,DFTI_CONJUGATE_EVEN_STORAGE,DFTI_COMPLEX_COMPLEX)

    WRITE_MKL_DEBUG_INFO(iret)

    IF (FFT_FORWARD == direction) THEN
       scale=dx1*dx2*dx3
       iret=DftiSetValue(desc,DFTI_FORWARD_SCALE,scale)
       iret=DftiSetValue(desc,DFTI_INPUT_STRIDES,rstrides);
       iret=DftiSetValue(desc,DFTI_OUTPUT_STRIDES,cstrides);
       iret=DftiCommitDescriptor(desc)
       iret=DftiComputeForward(desc,data)
    ELSE
       scale=1._4/(sx1*dx1*sx2*dx2*sx3*dx3)
       iret=DftiSetValue(desc,DFTI_BACKWARD_SCALE,scale)
       iret=DftiSetValue(desc,DFTI_INPUT_STRIDES,cstrides);
       iret=DftiSetValue(desc,DFTI_OUTPUT_STRIDES,rstrides);
       iret=DftiCommitDescriptor(desc)
       iret=DftiComputeBackward(desc,data)
    END IF
    iret=DftiFreeDescriptor(desc)
    WRITE_MKL_DEBUG_INFO(iret)

  END SUBROUTINE fft3
#else
  !------------------------------------------------------
  ! implementation of ctfft (N. Brenner, 1968)
  ! must be linked with ctfft.o
  !------------------------------------------------------
  SUBROUTINE fft3(data,sx1,sx2,sx3,dx1,dx2,dx3,direction)
    REAL*4, DIMENSION(:,:,:), INTENT(INOUT) :: data
    REAL*8, INTENT(IN) :: dx1,dx2,dx3
    INTEGER, INTENT(IN) :: sx1,sx2,sx3,direction

    INTEGER :: dim(3)
    INTEGER :: FOURT_DS ! data storage
    INTEGER, PARAMETER :: FOURT_NW = 128 ! extra work space size
    REAL*4, DIMENSION(FOURT_NW) :: FOURT_WORK ! extra work space

    dim=(/ sx1,sx2,sx3 /)

    IF (FFT_FORWARD == direction) THEN
       FOURT_DS=0
    ELSE
       FOURT_DS=-1
    END IF
    CALL ctfft(data,dim,3,direction,FOURT_DS,FOURT_WORK,FOURT_NW)

    IF (FFT_INVERSE == direction) THEN
       data=data/(sx1*dx1*sx2*dx2*sx3*dx3)
    ELSE
       data=data*(dx1*dx2*dx3)
    END IF

  END SUBROUTINE fft3
#endif
#endif
#endif
  !----------------------------------------------------------------------
  ! subroutine FFT2 performs normalized forward and
  ! inverse fourier transforms of real 2d data
  !
  ! USES subroutine FOURT
  ! ctfft(data,n,ndim,isign,iform,work,nwork)
  ! or
  ! fftw3
  !
  ! for real array the fourier transform returns a sx1/2+1 complex array
  ! and the enough space must be reserved
  !----------------------------------------------------------------------
#ifdef FFTW3
  SUBROUTINE fft2(data,sx1,sx2,dx1,dx2,direction)
    INTEGER, INTENT(IN) :: sx1,sx2,direction
    REAL*4, DIMENSION(sx1+2,sx2), INTENT(INOUT) :: data
    REAL*8, INTENT(IN) :: dx1,dx2

    INTEGER*8 :: plan

    IF (FFT_FORWARD == direction) THEN
      CALL sfftw_plan_dft_r2c_2d(plan,sx1,sx2, &
           data(1,1),data(1,1),FFTW_ESTIMATE)
    ELSE
      CALL sfftw_plan_dft_c2r_2d(plan,sx1,sx2, &
           data(1,1),data(1,1),FFTW_ESTIMATE)
    END IF

    CALL sfftw_execute(plan)
    CALL sfftw_destroy_plan(plan)

    IF (FFT_INVERSE == direction) THEN
      data=data/(sx1*dx1*sx2*dx2)
    ELSE
      data=data*(dx1*dx2)
    END IF

  END SUBROUTINE fft2
#else
#ifdef SGI_FFT
  SUBROUTINE fft2(data,sx1,sx2,dx1,dx2,direction)
    REAL*4, DIMENSION(:,:), INTENT(INOUT) :: data
    REAL*8, INTENT(IN) :: dx1,dx2
    INTEGER, INTENT(IN) :: sx1,sx2,direction

    INTEGER, PARAMETER :: NF=256, NFR=256

    REAL*4, DIMENSION(sx1+NFR+2*sx2+NF) :: table
    REAL*4, DIMENSION(sx1+4*sx2) :: work
    INTEGER, DIMENSION(2) :: isys
    REAL*4 :: scale

    isys(1)=1

    IF (FFT_FORWARD == direction) THEN
       scale=dx1*dx2
       ! initialize the sin/cos table
       CALL SCFFT2D(+0,sx1,sx2,scale,data(1,1),sx1+2, &
                    data(1,1),sx1/2+1,table,work,isys)
       CALL SCFFT2D(-1,sx1,sx2,scale,data(1,1),sx1+2, &
                    data(1,1),sx1/2+1,table,work,isys)
    ELSE
       scale=1._4/(sx1*dx1*sx2*dx2)
       ! initialize the sin/cos table
       CALL CSFFT2D(+0,sx1,sx2,scale,data(1,1),sx1/2+1, &
                    data(1,1),sx1+2,table,work,isys)
       CALL CSFFT2D(+1,sx1,sx2,scale,data(1,1),sx1/2+1, &
                    data(1,1),sx1+2,table,work,isys)
    END IF

  END SUBROUTINE fft2
#else
#ifdef IMKL_FFT
  !------------------------------------------------------
  ! implementation IMKL_FFT (Intel Math Kernel Library)
  ! for information and example calculations with the
  ! mkl FFT, see:
  !
  ! http://www.intel.com/software/products/mkl/ ...
  !                      docs/webhelp/appendices/ ...
  !                      mkl_appC_DFT.html#appC-exC-25
  !
  ! sylvain barbot (04-30-10) - original form
  !------------------------------------------------------
  SUBROUTINE fft2(data,sx1,sx2,dx1,dx2,direction)
    REAL*4, DIMENSION(0:*), INTENT(INOUT) :: data
    REAL*8, INTENT(IN) :: dx1,dx2
    INTEGER, INTENT(IN) :: sx1,sx2,direction

    INTEGER :: iret,size(2),rstrides(3),cstrides(3)
    TYPE(DFTI_DESCRIPTOR), POINTER :: desc
    REAL*4 :: scale

    rstrides=(/ 0,1,sx1+2 /)
    cstrides=(/ 0,1,sx1/2+1 /)
    size=(/ sx1,sx2 /)

    iret=DftiCreateDescriptor(desc,DFTI_SINGLE,DFTI_REAL,2,size);
    iret=DftiSetValue(desc,DFTI_CONJUGATE_EVEN_STORAGE,DFTI_COMPLEX_COMPLEX)

    WRITE_MKL_DEBUG_INFO(iret)

    IF (FFT_FORWARD == direction) THEN
       scale=dx1*dx2
       iret=DftiSetValue(desc,DFTI_FORWARD_SCALE,scale)
       iret=DftiSetValue(desc,DFTI_INPUT_STRIDES,rstrides);
       iret=DftiSetValue(desc,DFTI_OUTPUT_STRIDES,cstrides);
       iret=DftiCommitDescriptor(desc)
       iret=DftiComputeForward(desc,data)
    ELSE
       scale=1._4/(sx1*dx1*sx2*dx2)
       iret=DftiSetValue(desc,DFTI_BACKWARD_SCALE,scale)
       iret=DftiSetValue(desc,DFTI_INPUT_STRIDES,cstrides);
       iret=DftiSetValue(desc,DFTI_OUTPUT_STRIDES,rstrides);
       iret=DftiCommitDescriptor(desc)
       iret=DftiComputeBackward(desc,data)
    END IF
    iret=DftiFreeDescriptor(desc)
    WRITE_MKL_DEBUG_INFO(iret)

  END SUBROUTINE fft2
#else
  !------------------------------------------------------
  ! Couley-Tuckey implementation of the Fourier 
  ! transform with built-in FFT code (ctfft.f).
  !------------------------------------------------------
  SUBROUTINE fft2(data,sx1,sx2,dx1,dx2,direction)
    REAL*4, DIMENSION(:,:), INTENT(INOUT) :: data
    REAL*8, INTENT(IN) :: dx1,dx2
    INTEGER, INTENT(IN) :: sx1,sx2,direction

    INTEGER :: dim(2)
    INTEGER :: FOURT_DS ! data storage
    INTEGER, PARAMETER :: FOURT_NW = 64 ! extra work space size
    REAL*4, DIMENSION(FOURT_NW) :: FOURT_WORK ! extra work space

    dim=(/ sx1,sx2 /)

    IF (FFT_FORWARD == direction) THEN
       FOURT_DS=0
    ELSE
       FOURT_DS=-1
    END IF
    CALL ctfft(data,dim,2,direction,FOURT_DS,FOURT_WORK,FOURT_NW)

    IF (FFT_INVERSE == direction) THEN
       data=data/(sx1*dx1*sx2*dx2)
    ELSE
       data=data*(dx1*dx2)
    END IF

  END SUBROUTINE fft2
#endif
#endif
#endif

  !-----------------------------------------------------------------
  ! subroutine FFT1
  ! performs a one dimensional complex to complex Fourier
  ! transform
  !
  ! uses complex DFT ctfft (N. Brenner, 1968) by default
  ! or CCFFT (SGI library) with compile flag SGI_FFT
  !
  ! sylvain barbot (05-02-07) - original form
  !-----------------------------------------------------------------
#ifdef SGI_FFT
  !------------------------------------------------------
  ! implementation CCFFT
  !
  ! sylvain barbot (09-28-08) - original form
  !------------------------------------------------------
  SUBROUTINE fft1(data,sx,dx,direction)
    INTEGER, INTENT(IN) :: sx,direction
    COMPLEX(KIND=4), DIMENSION(:), INTENT(INOUT) :: data
    REAL*8, INTENT(IN) :: dx

    INTEGER, PARAMETER :: NF=256

    REAL*4, DIMENSION(2*sx+NF) :: table
    REAL*4, DIMENSION(2*sx) :: work
    INTEGER, DIMENSION(2) :: isys
    REAL*4 :: scale

    isys(1)=1

    IF (FFT_FORWARD == direction) THEN
       scale=dx
       ! initialize the sin/cos table
       CALL CCFFT(+0,sx,scale,data,data,table,work,isys)
       CALL CCFFT(-1,sx,scale,data,data,table,work,isys)
    ELSE
       scale=1._4/(sx*dx)
       ! initialize the sin/cos table
       CALL CCFFT(+0,sx,scale,data,data,table,work,isys)
       CALL CCFFT(+1,sx,scale,data,data,table,work,isys)
    END IF

  END SUBROUTINE fft1
#else
#ifdef IMKL_FFT
  !------------------------------------------------------
  ! implementation IMKL_FFT (Intel Math Kernel Library)
  ! evaluates a complex-to-complex Fourier transform
  !
  ! sylvain barbot (04-30-10) - original form
  !------------------------------------------------------
  SUBROUTINE fft1(data,sx,dx,direction)
    INTEGER, INTENT(IN) :: sx,direction
    COMPLEX(KIND=4), DIMENSION(0:*), INTENT(INOUT) :: data
    REAL*8, INTENT(IN) :: dx

    INTEGER :: iret
    TYPE(DFTI_DESCRIPTOR), POINTER :: desc

    REAL*4 :: scale

    iret=DftiCreateDescriptor(desc,DFTI_SINGLE,DFTI_COMPLEX,1,sx)
    WRITE_MKL_DEBUG_INFO(iret)

    IF (FFT_FORWARD == direction) THEN
       scale=dx
       iret=DftiSetValue(desc,DFTI_FORWARD_SCALE,scale)
       iret=DftiCommitDescriptor(desc)
       iret=DftiComputeForward(desc,data)
    ELSE
       scale=1._4/(sx*dx)
       iret=DftiSetValue(desc,DFTI_BACKWARD_SCALE,scale)
       iret=DftiCommitDescriptor(desc)
       iret=DftiComputeBackward(desc,data)
    END IF
    iret=DftiFreeDescriptor(desc)
    WRITE_MKL_DEBUG_INFO(iret)

  END SUBROUTINE fft1
#else
  !----------------------------------------------------
  ! implementation ctfft
  !
  ! sylvain barbot (05-02-07) - original form
  !----------------------------------------------------
  SUBROUTINE fft1(data,sx,dx,direction)
    COMPLEX(KIND=4),DIMENSION(:), INTENT(INOUT) :: data
    REAL*8, INTENT(IN) :: dx
    INTEGER, INTENT(IN) :: sx,direction

    INTEGER, PARAMETER :: FOURT_NW = 32 ! extra work space size
    REAL*4, DIMENSION(FOURT_NW) :: FOURT_WORK ! extra work space
    INTEGER :: FOURT_DS = 1

    CALL ctfft(data,sx,1,direction,FOURT_DS,FOURT_WORK,FOURT_NW)
    IF (FFT_INVERSE == direction) THEN
       data=data/(sx*dx)
    ELSE
       data=data*dx
    END IF

  END SUBROUTINE fft1
#endif
#endif

END MODULE fourier
