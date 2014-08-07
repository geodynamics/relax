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

MODULE green

  USE fourier

  IMPLICIT NONE

#include "include.f90"

  PUBLIC
  REAL*8, PRIVATE, PARAMETER :: pi   = 3.141592653589793115997963468544185161_8
  REAL*8, PRIVATE, PARAMETER :: pi2  = 6.28318530717958623199592693708837032318_8
  REAL*8, PRIVATE, PARAMETER :: pid2 = 1.57079632679489655799898173427209258079_8
    
  INTEGER, PARAMETER :: GRN_IMAGE=1,GRN_HS=0

CONTAINS

  !------------------------------------------------------------------------
  !> Subroutine ElasticResponse
  !! apply the 2d elastic (half-space) transfert function
  !! to the set of body forces.
  !!
  !! INPUT:
  !! @param mu          shear modulus
  !! @param f1,2,3      equivalent body-forces in the Fourier domain
  !! @param dx1,2,3     sampling size
  !!
  !! \author sylvain barbot (04/14/07) - original form
  !!                        (02/06/09) - parallel implementation with MPI and OpenMP
  !!                        (01/06/11) - remove implementation with MPI
  !------------------------------------------------------------------------
  SUBROUTINE elasticresponse(lambda,mu,f1,f2,f3,dx1,dx2,dx3)
    REAL*8, INTENT(IN) :: lambda,mu,dx1,dx2,dx3
    REAL*4, DIMENSION(:,:,:), INTENT(INOUT) :: f1,f2,f3
    
    REAL*8 :: k1,k2,k3,denom,r2,ratio1,ratio2
    INTEGER :: i1,i2,i3,sx1,sx2,sx3,ubound3
    COMPLEX(kind=8) :: buf1,buf2,buf3,c1,c2,c3
    
    sx1=SIZE(f2,1)-2
    sx2=SIZE(f2,2)
    sx3=SIZE(f2,3)
    
    ratio1=(lambda+mu)/(lambda+2._8*mu)/mu/(pi2**2._8)
    ratio2=mu/(lambda+mu)
    
    ubound3=sx3

    ! serial computation
!$omp parallel do private(i1,i2,k1,k2,k3,r2,denom,c1,c2,c3,buf1,buf2,buf3)
    DO i3=1,ubound3
       CALL wavenumber3(i3,sx3,dx3,k3)
       DO i2=1,sx2
          CALL wavenumber2(i2,sx2,dx2,k2)
          DO i1=1,sx1/2+1
             CALL wavenumber1(i1,sx1,dx1,k1)
             
             r2=k1**2._8+k2**2._8+k3**2._8
             denom=ratio1/r2**2
             
             c1=CMPLX(f1(2*i1-1,i2,i3),f1(2*i1,i2,i3),8)
             c2=CMPLX(f2(2*i1-1,i2,i3),f2(2*i1,i2,i3),8)
             c3=CMPLX(f3(2*i1-1,i2,i3),f3(2*i1,i2,i3),8)
             
             buf1=((k2**2._8+k3**2._8+ratio2*r2)*c1-k1*(k2*c2+k3*c3))*denom
             buf2=((k1**2._8+k3**2._8+ratio2*r2)*c2-k2*(k1*c1+k3*c3))*denom
             buf3=((k1**2._8+k2**2._8+ratio2*r2)*c3-k3*(k1*c1+k2*c2))*denom
             
             f1(2*i1-1:2*i1,i2,i3)=REAL((/ REAL(buf1),AIMAG(buf1) /))
             f2(2*i1-1:2*i1,i2,i3)=REAL((/ REAL(buf2),AIMAG(buf2) /))
             f3(2*i1-1:2*i1,i2,i3)=REAL((/ REAL(buf3),AIMAG(buf3) /))
          END DO
       END DO
    END DO
!$omp end parallel do

    ! zero wavenumber, no net body-force
    f1(1:2,1,1)=(/ 0._4, 0._4 /)
    f2(1:2,1,1)=(/ 0._4, 0._4 /)
    f3(1:2,1,1)=(/ 0._4, 0._4 /)

  END SUBROUTINE elasticresponse

  !---------------------------------------------------------------------
  !> subroutine SurfaceNormalTraction
  !! computes the two-dimensional field of surface normal stress
  !! expressed in the Fourier domain.
  !! The surface (x3=0) solution is obtained by integrating over the
  !! wavenumbers in 3-direction in the Fourier domain.
  !!
  !! \author sylvain barbot (05-01-07) - original form
  !---------------------------------------------------------------------
  SUBROUTINE surfacenormaltraction(lambda, mu, u1, u2, u3, dx1, dx2, dx3, p)
    REAL*4, INTENT(IN), DIMENSION(:,:,:) :: u1, u2, u3
    REAL*8, INTENT(IN) :: lambda, mu, dx1, dx2, dx3
    REAL*4, INTENT(OUT), DIMENSION(:,:) :: p
    
    INTEGER :: i1, i2, i3, sx1, sx2, sx3
    REAL*8 :: k1, k2, k3, modulus
    COMPLEX*8, PARAMETER :: i = CMPLX(0._8,pi2)
    COMPLEX(KIND=8) :: sum, c1, c2, c3
    
    sx1=SIZE(u1,1)-2
    sx2=SIZE(u1,2)
    sx3=SIZE(u1,3)
    
    modulus=lambda+2*mu
    
    p=0
    DO i3=1,sx3
       DO i2=1,sx2
          DO i1=1,sx1/2+1
             CALL wavenumbers(i1,i2,i3,sx1,sx2,sx3,dx1,dx2,dx3,k1,k2,k3)
             
             c1=CMPLX(u1(2*i1-1,i2,i3),u1(2*i1,i2,i3))
             c2=CMPLX(u2(2*i1-1,i2,i3),u2(2*i1,i2,i3))
             c3=CMPLX(u3(2*i1-1,i2,i3),u3(2*i1,i2,i3))
             
             sum=i*(modulus*k3*c3+lambda*(k1*c1+k2*c2))
             
             p(2*i1-1,i2)=p(2*i1-1,i2)+REAL( REAL(sum))
             p(2*i1  ,i2)=p(2*i1  ,i2)+REAL(AIMAG(sum))
          END DO
       END DO
    END DO
    p=REAL(p/(sx3*dx3),4)
    
  END SUBROUTINE surfacenormaltraction

  !---------------------------------------------------------------------
  !> subroutine Boussinesq3D
  !! computes the deformation field in the 3-dimensional grid
  !! due to a normal stress at the surface. Apply the Fourier domain
  !! solution of Steketee [1958].
  !---------------------------------------------------------------------
  SUBROUTINE boussinesq3d(p,lambda,mu,u1,u2,u3,dx1,dx2,dx3)
    REAL*4, DIMENSION(:,:), INTENT(IN) :: p
    REAL*4, DIMENSION(:,:,:), INTENT(INOUT) :: u1, u2, u3
    REAL*8, INTENT(IN) :: lambda, mu, dx1, dx2, dx3

    INTEGER :: i1, i2, i3, sx1, sx2, sx3, status
    REAL*8 :: k1, k2, k3, x3, alpha
    COMPLEX, ALLOCATABLE, DIMENSION(:) :: b1, b2, b3
    COMPLEX :: load

    sx1=SIZE(u1,1)-2
    sx2=SIZE(u1,2)
    sx3=SIZE(u1,3)
    
    ALLOCATE(b1(sx3),b2(sx3),b3(sx3),STAT=status)
    IF (0/=status) STOP "could not allocate arrays for Boussinesq3D"
    
    alpha=(lambda+mu)/(lambda+2*mu)

    DO i2=1,sx2
       DO i1=1,sx1/2+1
          CALL wavenumbers(i1,i2,1,sx1,sx2,1,dx1,dx2,1._8,k1,k2,k3)
          load=CMPLX(p(2*i1-1,i2),p(2*i1,i2))
          DO i3=1,sx3
             IF (i3<=sx3/2) THEN
                x3=DBLE(i3-1)*dx3
             ELSE
                x3=ABS(DBLE(i3-sx3-1)*dx3)
             END IF
             CALL steketeesolution(load,alpha,b1(i3),b2(i3),b3(i3),k1,k2,x3)
          END DO
          
          ! transforms the Steketee solution into a full 3-dimensional
          ! Fourier transform by 1d transforming in the 3-direction
          CALL fft1(b1,sx3,dx3,FFT_FORWARD)
          CALL fft1(b2,sx3,dx3,FFT_FORWARD)
          CALL fft1(b3,sx3,dx3,FFT_FORWARD)
          
          ! add the Boussinesq contribution to the deformation field
          DO i3=1,sx3
             u1(2*i1-1:2*i1,i2,i3)=u1(2*i1-1:2*i1,i2,i3)+ &
                  (/REAL(b1(i3)),AIMAG(b1(i3))/)
             u2(2*i1-1:2*i1,i2,i3)=u2(2*i1-1:2*i1,i2,i3)+ &
                  (/REAL(b2(i3)),AIMAG(b2(i3))/)
             u3(2*i1-1:2*i1,i2,i3)=u3(2*i1-1:2*i1,i2,i3)+ &
                  (/REAL(b3(i3)),AIMAG(b3(i3))/)
          END DO
       END DO
    END DO

    DEALLOCATE(b1,b2,b3)
    
    CONTAINS
      !-----------------------------------------------------------------
      !> subroutine SteketeeSolution
      !! computes the spectrum (two-dimensional Fourier transform)
      !! of the 3 components of the deformation field u1, u2, u3
      !! at wavenumbers k1, k2 and position x3. This is the analytical
      !! solution of [J. A. Steketee, On Volterra's dislocations in a
      !! semi-infinite elastic medium, Canadian Journal of Physics, 1958]
      !!
      !! \author sylvain barbot (05-02-07) - original form
      !-----------------------------------------------------------------
      SUBROUTINE steketeesolution(p,alpha,u1,u2,u3,k1,k2,x3)
        COMPLEX, INTENT(INOUT) :: u1, u2, u3
        REAL*8, INTENT(IN) :: alpha, k1, k2, x3
        COMPLEX, INTENT(IN) :: p
        
        REAL*8 :: beta, depthdecay
        COMPLEX, PARAMETER :: i=CMPLX(0,1)
        COMPLEX :: b
        
        beta=pi2*sqrt(k1**2._8+k2**2._8)
        depthdecay=exp(-beta*abs(x3))
        
        IF (0==k1 .AND. 0==k2) THEN
           u1=CMPLX(0.,0.)
           u2=CMPLX(0.,0.)
           u3=CMPLX(0.,0.)
        ELSE
           b=CMPLX(p/(2._8*mu*alpha*beta**3._8))
           u1=CMPLX(i*alpha*pi2*beta*b*(1._8-1._8/alpha+beta*x3)*depthdecay)
           u2=u1
           u1=CMPLX(u1*k1)
           u2=CMPLX(u2*k2)
           u3=CMPLX(-p/(2*mu*beta)*(1._8/alpha+beta*x3)*depthdecay)
        END IF
        
      END SUBROUTINE steketeesolution

  END SUBROUTINE boussinesq3d

  !---------------------------------------------------------------------
  !> subroutine SurfaceTraction
  !! computes the two-dimensional field of surface normal stress
  !! expressed in the Fourier domain.
  !! The surface (x3=0) solution is obtained by integrating over the
  !! wavenumbers in 3-direction in the Fourier domain.
  !!
  !! \author sylvain barbot (07-07-07) - original form
  !                         (02-09-09) - parallelized with mpi and openmp
  !---------------------------------------------------------------------
  SUBROUTINE surfacetraction(lambda,mu,u1,u2,u3,dx1,dx2,dx3,p1,p2,p3)
    REAL*4, INTENT(IN), DIMENSION(:,:,:) :: u1,u2,u3
    REAL*8, INTENT(IN) :: lambda,mu,dx1,dx2,dx3
    REAL*4, INTENT(OUT), DIMENSION(:,:) :: p1,p2,p3

    INTEGER :: i1,i2,i3,sx1,sx2,sx3
    REAL*8 :: k1,k2,k3,modulus
    COMPLEX(KIND=8), PARAMETER :: i=CMPLX(0._8,pi2,8)
    COMPLEX(KIND=8) :: sum1,sum2,sum3,c1,c2,c3

    sx1=SIZE(u1,1)-2
    sx2=SIZE(u1,2)
    sx3=SIZE(u1,3)

    modulus=lambda+2._8*mu

    p1=0
    p2=0
    p3=0

!$omp parallel do private(i1,i2,k1,k2,k3,c1,c2,c3,sum1,sum2,sum3), &
!$omp reduction(+:p1,p2,p3)
    DO i3=1,sx3
       DO i2=1,sx2
          DO i1=1,sx1/2+1
             CALL wavenumbers(i1,i2,i3,sx1,sx2,sx3,dx1,dx2,dx3,k1,k2,k3)

             c1=CMPLX(u1(2*i1-1,i2,i3),u1(2*i1,i2,i3),8)
             c2=CMPLX(u2(2*i1-1,i2,i3),u2(2*i1,i2,i3),8)
             c3=CMPLX(u3(2*i1-1,i2,i3),u3(2*i1,i2,i3),8)

             sum1=i*mu*(k3*c1+k1*c3)
             sum2=i*mu*(k3*c2+k2*c3)
             sum3=i*(modulus*k3*c3+lambda*(k1*c1+k2*c2))

             p1(2*i1-1:2*i1,i2)=p1(2*i1-1:2*i1,i2) &
                  +(/REAL(REAL(sum1)),REAL(AIMAG(sum1))/)
             p2(2*i1-1:2*i1,i2)=p2(2*i1-1:2*i1,i2) &
                  +(/REAL(REAL(sum2)),REAL(AIMAG(sum2))/)
             p3(2*i1-1:2*i1,i2)=p3(2*i1-1:2*i1,i2) &
                  +(/REAL(REAL(sum3)),REAL(AIMAG(sum3))/)

          END DO
       END DO
    END DO
!$omp end parallel do

    p1=p1/REAL(sx3*dx3,4)
    p2=p2/REAL(sx3*dx3,4)
    p3=p3/REAL(sx3*dx3,4)

  END SUBROUTINE surfacetraction

  !---------------------------------------------------------------------
  !> subroutine SurfaceTractionCowling
  !! computes the two-dimensional field of the resulting traction 
  !! expressed in the Fourier domain in the presence of gravity.
  !!
  !! The surface solution (x3=0) is obtained from the Fourier domain 
  !! array by integrating over the wavenumbers in 3-direction.
  !!
  !! The effective traction at x3=0 is 
  !!
  !!     t_1 = sigma_13
  !!     t_2 = sigma_23
  !!     t_3 = sigma_33 - r g u3
  !!         = sigma_33 - 2 mu alpha gamma u3
  !!
  !! \author sylvain barbot (07-07-07) - original form
  !---------------------------------------------------------------------
  SUBROUTINE surfacetractioncowling(lambda,mu,gamma,u1,u2,u3,dx1,dx2,dx3, &
       p1,p2,p3)
    REAL*4, INTENT(IN), DIMENSION(:,:,:) :: u1,u2,u3
    REAL*8, INTENT(IN) :: lambda,mu,gamma,dx1,dx2,dx3
    REAL*4, INTENT(OUT), DIMENSION(:,:) :: p1,p2,p3
    
    INTEGER :: i1,i2,i3,sx1,sx2,sx3
    REAL*8 :: k1,k2,k3,modulus,alpha,grav
    COMPLEX(KIND=8), PARAMETER :: i=CMPLX(0._8,pi2)
    COMPLEX(KIND=8) :: sum1,sum2,sum3,c1,c2,c3
    
    sx1=SIZE(u1,1)-2
    sx2=SIZE(u1,2)
    sx3=SIZE(u1,3)
    
    modulus=lambda+2._8*mu
    alpha=(lambda+mu)/(lambda+2._8*mu)
    grav=2._8*mu*alpha*gamma
    
    p1=0
    p2=0
    p3=0

!$omp parallel do private(i1,i3,k1,k2,k3,c1,c2,c3,sum1,sum2,sum3)
!!!$omp reduction(+:p1,p2,p3)
    DO i2=1,sx2
       DO i3=1,sx3
          DO i1=1,sx1/2+1
             CALL wavenumbers(i1,i2,i3,sx1,sx2,sx3,dx1,dx2,dx3,k1,k2,k3)
             
             c1=CMPLX(u1(2*i1-1,i2,i3),u1(2*i1,i2,i3))
             c2=CMPLX(u2(2*i1-1,i2,i3),u2(2*i1,i2,i3))
             c3=CMPLX(u3(2*i1-1,i2,i3),u3(2*i1,i2,i3))

             sum1=i*mu*(k3*c1+k1*c3)
             sum2=i*mu*(k3*c2+k2*c3)
             sum3=i*(modulus*k3*c3+lambda*(k1*c1+k2*c2))-grav*c3
             
             p1(2*i1-1:2*i1,i2)=p1(2*i1-1:2*i1,i2)+(/REAL(sum1,4),REAL(AIMAG(sum1),4)/)
             p2(2*i1-1:2*i1,i2)=p2(2*i1-1:2*i1,i2)+(/REAL(sum2,4),REAL(AIMAG(sum2),4)/)
             p3(2*i1-1:2*i1,i2)=p3(2*i1-1:2*i1,i2)+(/REAL(sum3,4),REAL(AIMAG(sum3),4)/)
          END DO
       END DO
    END DO
!$omp end parallel do

    p1=p1/REAL(sx3*dx3,4)
    p2=p2/REAL(sx3*dx3,4)
    p3=p3/REAL(sx3*dx3,4)
    
  END SUBROUTINE surfacetractioncowling

  !---------------------------------------------------------------------
  !> subroutine Cerruti3D
  !! computes the deformation field in the 3-dimensional grid
  !! due to an arbitrary surface traction.
  !!
  !! \author sylvain barbot (07/07/07) - original form
  !                (02/01/09) - parallelized with MPI and OpenMP
  !                (01/06/11) - remove parallelized version with MPI
  !---------------------------------------------------------------------
  SUBROUTINE cerruti3d(p1,p2,p3,lambda,mu,u1,u2,u3,dx1,dx2,dx3)
    REAL*4, DIMENSION(:,:), INTENT(IN) :: p1,p2,p3
    REAL*4, DIMENSION(:,:,:), INTENT(INOUT) :: u1,u2,u3
    REAL*8, INTENT(IN) :: lambda,mu,dx1,dx2,dx3

    INTEGER :: i1,i2,i3,ib,sx1,sx2,sx3,iostatus,buffersize
    REAL*8 :: k1,k2,k3,x3,alpha
    COMPLEX(KIND=4) :: t1,t2,t3
    INTEGER, PARAMETER :: stride=64
    COMPLEX(KIND=4), ALLOCATABLE, DIMENSION(:,:) :: b1,b2,b3

    sx1=SIZE(u1,1)-2
    sx2=SIZE(u1,2)
    sx3=SIZE(u1,3)

    alpha=(lambda+mu)/(lambda+2*mu)

    ! serial programmation implementation
!$omp parallel private(b1,b2,b3,iostatus)

    ALLOCATE(b1(sx3,stride),b2(sx3,stride),b3(sx3,stride),STAT=iostatus)
    IF (0/=iostatus) STOP "could not allocate arrays for Cerruti3D"

!$omp do private(i1,i3,ib,k1,k2,k3,t1,t2,t3,x3,buffersize)
    DO i2=1,sx2
       DO i1=1,sx1/2+1,stride

          ! buffer results
          IF (i1+stride-1 .GT. sx1/2+1) THEN
             buffersize=sx1/2+1-i1+1
          ELSE
             buffersize=stride
          END IF

          DO ib=0,buffersize-1

             CALL wavenumbers(i1+ib,i2,1,sx1,sx2,1,dx1,dx2,1._8,k1,k2,k3)
             t1=CMPLX(p1(2*(i1+ib)-1,i2),p1(2*(i1+ib),i2),4)
             t2=CMPLX(p2(2*(i1+ib)-1,i2),p2(2*(i1+ib),i2),4)
             t3=CMPLX(p3(2*(i1+ib)-1,i2),p3(2*(i1+ib),i2),4)

             DO i3=1,sx3
                IF (i3<=sx3/2) THEN
                   x3=DBLE(i3-1)*dx3
                ELSE
                   x3=ABS(DBLE(i3-sx3-1)*dx3)
                END IF
                CALL cerrutisolution(mu,t1,t2,t3,alpha,b1(i3,ib+1),b2(i3,ib+1),b3(i3,ib+1),k1,k2,x3)
             END DO

             ! transforms the Cerruti solution into a full 3-dimensional
             ! Fourier transform by 1d transforming in the 3-direction
             CALL fft1(b1(:,ib+1),sx3,dx3,FFT_FORWARD)
             CALL fft1(b2(:,ib+1),sx3,dx3,FFT_FORWARD)
             CALL fft1(b3(:,ib+1),sx3,dx3,FFT_FORWARD)

          END DO

          ! update solution displacement
          DO i3=1,sx3
             DO ib=0,buffersize-1
                u1(2*(i1+ib)-1,i2,i3)=u1(2*(i1+ib)-1,i2,i3)+REAL( REAL(b1(i3,ib+1)))
                u1(2*(i1+ib)  ,i2,i3)=u1(2*(i1+ib)  ,i2,i3)+REAL(AIMAG(b1(i3,ib+1)))
                u2(2*(i1+ib)-1,i2,i3)=u2(2*(i1+ib)-1,i2,i3)+REAL( REAL(b2(i3,ib+1)))
                u2(2*(i1+ib)  ,i2,i3)=u2(2*(i1+ib)  ,i2,i3)+REAL(AIMAG(b2(i3,ib+1)))
                u3(2*(i1+ib)-1,i2,i3)=u3(2*(i1+ib)-1,i2,i3)+REAL( REAL(b3(i3,ib+1)))
                u3(2*(i1+ib)  ,i2,i3)=u3(2*(i1+ib)  ,i2,i3)+REAL(AIMAG(b3(i3,ib+1)))
             END DO
          END DO

       END DO
    END DO

    DEALLOCATE(b1,b2,b3)
!$omp end parallel

    CONTAINS
      !-----------------------------------------------------------------
      !> subroutine CerrutiSolution
      !! computes the general solution for the deformation field in an
      !! elastic half-space due to an arbitrary surface traction.
      !! the 3 components u1, u2, u3 of the deformation field are
      !! expressed in the horizontal Fourier at depth x3.
      !! this combines the solution to the Boussinesq's and the Cerruti's
      !! problem in a half-space.
      !!
      !! \author sylvain barbot (07-07-07) - original form
      !-----------------------------------------------------------------
      SUBROUTINE cerrutisolution(mu,p1,p2,p3,alpha,u1,u2,u3,k1,k2,x3)
        COMPLEX(KIND=4), INTENT(INOUT) :: u1,u2,u3
        REAL*8, INTENT(IN) :: mu,alpha,k1,k2,x3
        COMPLEX(KIND=4), INTENT(IN) :: p1,p2,p3

        REAL*8 :: beta, depthdecay
        COMPLEX(KIND=8), PARAMETER :: i=CMPLX(0._8,pi2,8)
        REAL*8  :: temp
        COMPLEX(KIND=8) :: b1,b2,b3,tmp,v1,v2,v3

        beta=pi2*sqrt(k1**2+k2**2)
        depthdecay=exp(-beta*abs(x3))

        IF (0==k1 .AND. 0==k2) THEN
           u1=CMPLX(0._4,0._4,4)
           u2=CMPLX(0._4,0._4,4)
           u3=CMPLX(0._4,0._4,4)
        ELSE
           temp=1._8/(2._8*mu*beta**3)*depthdecay
           b1=temp*p1
           b2=temp*p2
           b3=temp*p3

           ! b3 contribution
           tmp=i*b3*(beta*(1._8-1._8/alpha+beta*x3))
           v1=tmp*k1
           v2=tmp*k2
           v3=-beta**2*b3*(1._8/alpha+beta*x3)

           ! b1 contribution
           temp=pi2**2*(2._8-1._8/alpha+beta*x3)
           v1=v1+b1*(-2._8*beta**2+k1**2*temp)
           v2=v2+b1*k1*k2*temp
           v3=v3+b1*i*k1*beta*(1._8/alpha-1._8+beta*x3)

           ! b2 contribution & switch to single-precision
           u1=CMPLX(v1+b2*k1*k2*temp)
           u2=CMPLX(v2+b2*(-2._8*beta**2+k2**2*temp))
           u3=CMPLX(v3+b2*i*k2*beta*(1._8/alpha-1._8+beta*x3))
        END IF

      END SUBROUTINE cerrutisolution
  END SUBROUTINE cerruti3d

  SUBROUTINE cerrutimodified(p1,p2,p3,lambda,mu,gamma,u1,u2,u3,dx1,dx2,dx3)
    REAL*4, DIMENSION(:,:), INTENT(IN) :: p1,p2,p3
    REAL*4, DIMENSION(:,:,:), INTENT(INOUT) :: u1,u2,u3
    REAL*8, INTENT(IN) :: lambda,mu,gamma,dx1,dx2,dx3

    INTEGER :: i,i1,i2,i3,ib,sx1,sx2,sx3,iostatus,buffersize
    REAL*8 :: k1,k2,k3,alpha
    COMPLEX(KIND=4) :: t1,t2,t3
    COMPLEX(KIND=4) :: b1,b2,b3
    REAL*8 :: taper

    sx1=SIZE(u1,1)-2
    sx2=SIZE(u1,2)
    sx3=SIZE(u1,3)

    alpha=(lambda+mu)/(lambda+2*mu)
  
!$omp parallel do private(i1,i3,k1,k2,k3,t1,t2,t3,b1,b2,b3) 
    DO i2=1,sx2
       DO i1=1,sx1/2+1
          DO i3=1, sx3

             CALL wavenumbers(i1,i2,i3,sx1,sx2,sx3,dx1,dx2,dx3,k1,k2,k3)

             t1=CMPLX(p1(2*(i1)-1,i2),p1(2*(i1),i2),4)
             t2=CMPLX(p2(2*(i1)-1,i2),p2(2*(i1),i2),4)
             t3=CMPLX(p3(2*(i1)-1,i2),p3(2*(i1),i2),4)

             CALL cerrutisolmodified(mu,t1,t2,t3,alpha,gamma, &
                  b1,b2,b3,k1,k2,k3,DBLE(sx3/2)*dx3)

             !taper=cos((pi*i3/(sx3))**(2))

             !b1=(b1*taper) 
             !b2=(b2*taper) 
             !b3=(b3*taper) 
             
             u1(2*(i1)-1,i2,i3)=u1(2*(i1)-1,i2,i3)+REAL( REAL(b1))
             u1(2*(i1)  ,i2,i3)=u1(2*(i1)  ,i2,i3)+REAL(AIMAG(b1))
             u2(2*(i1)-1,i2,i3)=u2(2*(i1)-1,i2,i3)+REAL( REAL(b2))
             u2(2*(i1)  ,i2,i3)=u2(2*(i1)  ,i2,i3)+REAL(AIMAG(b2))
             u3(2*(i1)-1,i2,i3)=u3(2*(i1)-1,i2,i3)+REAL( REAL(b3))
             u3(2*(i1)  ,i2,i3)=u3(2*(i1)  ,i2,i3)+REAL(AIMAG(b3))

          END DO
       END DO
    END DO
!$omp end parallel do

    CONTAINS 

    SUBROUTINE cerrutisolmodified(mu,p1,p2,p3,alpha,gamma,u1,u2,u3,k1,k2,k3,L)
        COMPLEX(KIND=4), INTENT(INOUT) :: u1,u2,u3
        REAL*8, INTENT(IN) :: mu,alpha,gamma,k1,k2,k3,L
        COMPLEX(KIND=4), INTENT(IN) :: p1,p2,p3

        REAL*8 :: beta, fi, h
        COMPLEX(KIND=8), PARAMETER :: i=CMPLX(0._8,1._8)
        REAL*8  :: temp
        COMPLEX(KIND=8) :: b1,b2,b3,tmp,v1,v2,v3

        beta=pi2*sqrt(k1**2+k2**2)
        fi = (2*beta)/((pi2**2)*(k3**2)+beta**2)
        h=gamma/beta

        IF (0==k1 .AND. 0==k2) THEN
           ! Check these 
           !u1=CMPLX(REAL(+p1/mu*(k3-L)),0._4)
           !u2=CMPLX(REAL(+p2/mu*(k3-L)),0._4)
           !u3=CMPLX(REAL(+p3/mu*(k3-L)*(1.d0-alpha)/(1.d0+2.d0*L*alpha*gamma*(1.d0-alpha))),0._4)
           u1=CMPLX(0._4,0._4)
           u2=CMPLX(0._4,0._4)
           u3=CMPLX(0._4,0._4)
        ELSE
           temp=(1._8/(2._8*mu*beta**3))*fi
           b1=temp*p1
           b2=temp*p2
           b3=(beta*p3-i*(1._8-alpha)*(pi2*k1*p1+pi2*k2*p2))/(2._8*alpha*mu*beta**4*(1+h))

           tmp=alpha*i*beta*pi2*b3*(1._8-1._8/alpha-i*pi2*k3*fi)
           v1=tmp*k1*fi
           v2=tmp*k2*fi
           v3=-alpha*beta**2*b3*(1._8/alpha-i*fi*pi2*k3)

           u1=CMPLX(v1+(-(2._8*beta**2*b1)+(alpha*(pi2**2)*k1*(b1*k1+b2*k2)*(1._8-(i*pi2*k3*fi)))))
           u2=CMPLX(v2+(-(2._8*beta**2*b2)+(alpha*(pi2**2)*k2*(b1*k1+b2*k2)*(1._8-(i*pi2*k3*fi)))))
           u3=CMPLX(v3*fi+alpha*beta*pi2*(k1*b1+k2*b2)*pi2*k3*fi)
        
        END IF

      END SUBROUTINE cerrutisolmodified

  END SUBROUTINE cerrutimodified


  !---------------------------------------------------------------------
  !> subroutine CerrutiCowling
  !! computes the deformation field in the 3-dimensional grid
  !! due to an arbitrary surface traction.
  !!
  !! \author sylvain barbot - (07/07/07) - original form
  !!                          (21/11/08) - gravity effect
  !!                          (02/01/09) - parallelized with MPI and OpenMP
  !!                          (01/06/11) - remove parallelized version with MPI
  !---------------------------------------------------------------------
  SUBROUTINE cerruticowling(p1,p2,p3,lambda,mu,gamma,u1,u2,u3,dx1,dx2,dx3)
    REAL*4, DIMENSION(:,:), INTENT(IN) :: p1,p2,p3
    REAL*4, DIMENSION(:,:,:), INTENT(INOUT) :: u1,u2,u3
    REAL*8, INTENT(IN) :: lambda,mu,gamma,dx1,dx2,dx3

    INTEGER :: i1,i2,i3,ib,sx1,sx2,sx3,iostatus,buffersize
    REAL*8 :: k1,k2,k3,x3,alpha
    COMPLEX(KIND=4) :: t1,t2,t3
    INTEGER, PARAMETER :: stride=64
    COMPLEX(KIND=4), ALLOCATABLE, DIMENSION(:,:) :: b1,b2,b3

    sx1=SIZE(u1,1)-2
    sx2=SIZE(u1,2)
    sx3=SIZE(u1,3)

    alpha=(lambda+mu)/(lambda+2*mu)

    ! serial programmation implementation
!$omp parallel private(b1,b2,b3,iostatus)

    ALLOCATE(b1(sx3,stride),b2(sx3,stride),b3(sx3,stride),STAT=iostatus)
    IF (0/=iostatus) STOP "could not allocate arrays for Cerruti3D"

!$omp do private(i1,i3,ib,k1,k2,k3,t1,t2,t3,x3,buffersize)
    DO i2=1,sx2
       DO i1=1,sx1/2+1,stride

          ! buffer results
          IF (i1+stride-1 .GT. sx1/2+1) THEN
             buffersize=sx1/2+1-i1+1
          ELSE
             buffersize=stride
          END IF

          DO ib=0,buffersize-1

             CALL wavenumbers(i1+ib,i2,1,sx1,sx2,1,dx1,dx2,1._8,k1,k2,k3)
             t1=CMPLX(p1(2*(i1+ib)-1,i2),p1(2*(i1+ib),i2),4)
             t2=CMPLX(p2(2*(i1+ib)-1,i2),p2(2*(i1+ib),i2),4)
             t3=CMPLX(p3(2*(i1+ib)-1,i2),p3(2*(i1+ib),i2),4)

             DO i3=1,sx3
                IF (i3<=sx3/2) THEN
                   x3=DBLE(i3-1)*dx3
                ELSE
                   x3=ABS(DBLE(i3-sx3-1)*dx3)
                END IF
                CALL cerrutisolcowling(mu,t1,t2,t3,alpha,gamma, &
                     b1(i3,ib+1),b2(i3,ib+1),b3(i3,ib+1),k1,k2,x3,DBLE(sx3/2)*dx3)
             END DO

             ! transforms the Cerruti solution into a full 3-dimensional
             ! Fourier transform by 1d transforming in the 3-direction
             CALL fft1(b1(:,ib+1),sx3,dx3,FFT_FORWARD)
             CALL fft1(b2(:,ib+1),sx3,dx3,FFT_FORWARD)
             CALL fft1(b3(:,ib+1),sx3,dx3,FFT_FORWARD)

          END DO

          ! update solution displacement
          DO i3=1,sx3
             DO ib=0,buffersize-1
                u1(2*(i1+ib)-1,i2,i3)=u1(2*(i1+ib)-1,i2,i3)+REAL( REAL(b1(i3,ib+1)))
                u1(2*(i1+ib)  ,i2,i3)=u1(2*(i1+ib)  ,i2,i3)+REAL(AIMAG(b1(i3,ib+1)))
                u2(2*(i1+ib)-1,i2,i3)=u2(2*(i1+ib)-1,i2,i3)+REAL( REAL(b2(i3,ib+1)))
                u2(2*(i1+ib)  ,i2,i3)=u2(2*(i1+ib)  ,i2,i3)+REAL(AIMAG(b2(i3,ib+1)))
                u3(2*(i1+ib)-1,i2,i3)=u3(2*(i1+ib)-1,i2,i3)+REAL( REAL(b3(i3,ib+1)))
                u3(2*(i1+ib)  ,i2,i3)=u3(2*(i1+ib)  ,i2,i3)+REAL(AIMAG(b3(i3,ib+1)))
             END DO
          END DO

       END DO
    END DO

    DEALLOCATE(b1,b2,b3)
!$omp end parallel

    CONTAINS

      !------------------------------------------------------------------------
      !> subroutine CerrutiSolCowling
      !! computes the general solution for the deformation field in an
      !! elastic half-space due to an arbitrary surface traction in the
      !! presence of gravity.
      !!
      !! The 3 components u1, u2 and u3 of the deformation field are 
      !! expressed in the horizontal Fourier domain at depth x3. 
      !!
      !! Combines the solution to the Boussinesq's and the Cerruti's 
      !! problem in a half-space with buoyancy boundary conditions.
      !!
      !! For the zero wavenumber, the solution is only depth-dependent and
      !! satisfies:
      !!
      !!                   -q_1
      !!   u1(k1, k2, x3) = --- (x3 - L)
      !!                     mu
      !!
      !!                   -q_2
      !!   u2(k1, k2, x3) = --- (x3 - L)
      !!                     mu
      !!
      !!                     1        - q_3 (1 - alpha) 
      !!   u3(k1, k2, x3) = -- -------------------------------- (x3 - L)
      !!                    mu  1 + 2 Gamma L alpha (1 - alpha)
      !!
      !! sylvain barbot (07-07-07) - original form
      !!                (08-30-10) - account for net surface traction
      !!                (04-08-12) - fix a sign error for the zero wavenumber
      !------------------------------------------------------------------------
      SUBROUTINE cerrutisolcowling(mu,p1,p2,p3,alpha,gamma,u1,u2,u3,k1,k2,x3,L)
        COMPLEX(KIND=4), INTENT(INOUT) :: u1,u2,u3
        REAL*8, INTENT(IN) :: mu,alpha,gamma,k1,k2,x3,L
        COMPLEX(KIND=4), INTENT(IN) :: p1,p2,p3
        
        REAL*8 :: beta, depthdecay, h
        COMPLEX(KIND=8), PARAMETER :: i=CMPLX(0._8,pi2)
        REAL*8  :: temp
        COMPLEX(KIND=8) :: b1,b2,b3,tmp,v1,v2,v3
        
        beta=pi2*sqrt(k1**2+k2**2)
        depthdecay=exp(-beta*abs(x3))
        h=gamma/beta
        
        IF (0==k1 .AND. 0==k2) THEN
           u1=CMPLX(REAL(+p1/mu*(x3-L)),0._4)
           u2=CMPLX(REAL(+p2/mu*(x3-L)),0._4)
           u3=CMPLX(REAL(+p3/mu*(x3-L)*(1.d0-alpha)/(1.d0+2.d0*L*alpha*gamma*(1.d0-alpha))),0._4)
           !u1=CMPLX(0._4,0._4)
           !u2=CMPLX(0._4,0._4)
           !u3=CMPLX(0._4,0._4)
        ELSE
           temp=1._8/(2._8*mu*beta**3)*depthdecay
           b1=temp*p1
           b2=temp*p2
           b3=temp*p3/(1+h)
           
           ! b3 contribution
           tmp=i*b3*(beta*(1._8-1._8/alpha+beta*x3))
           v1=tmp*k1
           v2=tmp*k2
           v3=-beta**2*b3*(1._8/alpha+beta*x3)
           
           ! b1 contribution
           temp=pi2**2*(2._8-1._8/alpha+beta*x3)/(1+h)
           v1=v1+b1*(-2._8*beta**2+k1**2*temp)
           v2=v2+b1*k1*k2*temp
           v3=v3+b1*i*k1*beta*(1._8/alpha-1._8+beta*x3)/(1+h)
           
           ! b2 contribution & switch to single-precision
           u1=CMPLX(v1+b2*k1*k2*temp)
           u2=CMPLX(v2+b2*(-2._8*beta**2+k2**2*temp))
           u3=CMPLX(v3+b2*i*k2*beta*(1._8/alpha-1._8+beta*x3)/(1+h))
        END IF

      END SUBROUTINE cerrutisolcowling

  END SUBROUTINE cerruticowling

  !---------------------------------------------------------------------
  !> subroutine CerrutiCowlingSerial
  !! computes the deformation field in the 3-dimensional grid
  !! due to an arbitrary surface traction. No parallel version.
  !
  ! sylvain barbot - 07/07/07 - original form
  !                  21/11/08 - gravity effect
  !---------------------------------------------------------------------
  SUBROUTINE cerruticowlingserial(p1,p2,p3,lambda,mu,gamma,u1,u2,u3,dx1,dx2,dx3)
    REAL*4, DIMENSION(:,:), INTENT(IN) :: p1,p2,p3
    REAL*4, DIMENSION(:,:,:), INTENT(INOUT) :: u1,u2,u3
    REAL*8, INTENT(IN) :: lambda,mu,gamma,dx1,dx2,dx3

    INTEGER :: i1,i2,i3,sx1,sx2,sx3,status
    REAL*8 :: k1,k2,k3,x3,alpha
    COMPLEX(KIND=4), ALLOCATABLE, DIMENSION(:) :: b1,b2,b3
    COMPLEX(KIND=4) :: t1,t2,t3

    sx1=SIZE(u1,1)-2
    sx2=SIZE(u1,2)
    sx3=SIZE(u1,3)
    
    ALLOCATE(b1(sx3),b2(sx3),b3(sx3),STAT=status)
    IF (0/=status) STOP "could not allocate arrays for Cerruti3D"
    
    alpha=(lambda+mu)/(lambda+2*mu)

    DO i2=1,sx2
       DO i1=1,sx1/2+1
          CALL wavenumbers(i1,i2,1,sx1,sx2,1,dx1,dx2,1._8,k1,k2,k3)
          t1=CMPLX(p1(2*i1-1,i2),p1(2*i1,i2))
          t2=CMPLX(p2(2*i1-1,i2),p2(2*i1,i2))
          t3=CMPLX(p3(2*i1-1,i2),p3(2*i1,i2))
          DO i3=1,sx3
             IF (i3<=sx3/2) THEN
                x3=DBLE(i3-1)*dx3
             ELSE
                x3=ABS(DBLE(i3-sx3-1)*dx3)
             END IF
             CALL cerrutisolcowling(t1,t2,t3,alpha,gamma, &
                  b1(i3),b2(i3),b3(i3),k1,k2,x3)
          END DO
          
          ! transforms the Cerruti solution into a full 3-dimensional
          ! Fourier transform by 1d transforming in the 3-direction
          CALL fft1(b1,sx3,dx3,FFT_FORWARD)
          CALL fft1(b2,sx3,dx3,FFT_FORWARD)
          CALL fft1(b3,sx3,dx3,FFT_FORWARD)
          
          ! add the Cerruti's contribution to the deformation field
          DO i3=1,sx3
             u1(2*i1-1,i2,i3)=u1(2*i1-1,i2,i3)+REAL( REAL(b1(i3)))
             u1(2*i1  ,i2,i3)=u1(2*i1  ,i2,i3)+REAL(AIMAG(b1(i3)))
             u2(2*i1-1,i2,i3)=u2(2*i1-1,i2,i3)+REAL( REAL(b2(i3)))
             u2(2*i1  ,i2,i3)=u2(2*i1  ,i2,i3)+REAL(AIMAG(b2(i3)))
             u3(2*i1-1,i2,i3)=u3(2*i1-1,i2,i3)+REAL( REAL(b3(i3)))
             u3(2*i1  ,i2,i3)=u3(2*i1  ,i2,i3)+REAL(AIMAG(b3(i3)))
          END DO
       END DO
    END DO
    
  CONTAINS
    !-----------------------------------------------------------------
    !> subroutine CerrutiSolCowling
    !! computes the general solution for the deformation field in an
    !! elastic half-space due to an arbitrary surface traction in the
    !! presence of gravity.
    !!
    !! The 3 components u1, u2 and u3 of the deformation field are 
    !! expressed in the horizontal Fourier at depth x3. 
    !!
    !! Combines the solution to the Boussinesq's and the Cerruti's 
    !! problem in a half-space with buoyancy boundary conditions.
    !
    ! sylvain barbot (07-07-07) - original form
    !-----------------------------------------------------------------
    SUBROUTINE cerrutisolcowling(p1,p2,p3,alpha,gamma,u1,u2,u3,k1,k2,x3)
      COMPLEX(KIND=4), INTENT(INOUT) :: u1,u2,u3
      REAL*8, INTENT(IN) :: alpha,gamma,k1,k2,x3
      COMPLEX(KIND=4), INTENT(IN) :: p1,p2,p3
        
      REAL*8 :: beta, depthdecay, h
      COMPLEX(KIND=8), PARAMETER :: i=CMPLX(0._8,pi2)
      REAL*8  :: temp
      COMPLEX(KIND=8) :: b1,b2,b3,tmp,v1,v2,v3
      
      beta=pi2*sqrt(k1**2+k2**2)
      depthdecay=exp(-beta*abs(x3))
      h=gamma/beta
      
      IF (0==k1 .AND. 0==k2) THEN
         u1=CMPLX(0._4,0._4)
         u2=CMPLX(0._4,0._4)
         u3=CMPLX(0._4,0._4)
      ELSE
         temp=1._8/(2._8*mu*beta**3)*depthdecay
         b1=temp*p1
         b2=temp*p2
         b3=temp*p3/(1+h)
           
         ! b3 contribution
         tmp=i*b3*(beta*(1._8-1._8/alpha+beta*x3))
         v1=tmp*k1
         v2=tmp*k2
         v3=-beta**2*b3*(1._8/alpha+beta*x3)
           
         ! b1 contribution
         temp=pi2**2*(2._8-1._8/alpha+beta*x3)/(1+h)
         v1=v1+b1*(-2._8*beta**2+k1**2*temp)
         v2=v2+b1*k1*k2*temp
         v3=v3+b1*i*k1*beta*(1._8/alpha-1._8+beta*x3)/(1+h)
           
         ! b2 contribution & switch to single-precision
         u1=CMPLX(v1+b2*k1*k2*temp)
         u2=CMPLX(v2+b2*(-2._8*beta**2+k2**2*temp))
         u3=CMPLX(v3+b2*i*k2*beta*(1._8/alpha-1._8+beta*x3)/(1+h))
      END IF

    END SUBROUTINE cerrutisolcowling

  END SUBROUTINE cerruticowlingserial

  !------------------------------------------------------------------
  !> subroutine GreenFunction
  !! computes (inplace) the displacement components due to a set of
  !! 3-D body-forces by application of the semi-analytic Green's
  !! function. The solution satisfies quasi-static Navier's equation
  !! including vanishing of normal traction at the surface.
  !!
  !! The surface traction can be made to vanish by application of
  !!   1) method of images + boussinesq problem (grn_method=GRN_IMAGE)
  !!   2) boussinesq's and cerruti's problems (grn_method=GRN_HS)
  !! in the first case, the body-forces are supposed by have been
  !! imaged appropriately.
  !
  ! sylvain barbot (07/07/07) - original form
  !------------------------------------------------------------------
  SUBROUTINE greenfunction(c1,c2,c3,t1,t2,t3,dx1,dx2,dx3,lambda,mu,grn_method)
    REAL*4, INTENT(INOUT), DIMENSION(:,:,:) :: c1,c2,c3
    REAL*4, INTENT(INOUT), DIMENSION(:,:) :: t1,t2,t3
    REAL*8, INTENT(IN) :: dx1,dx2,dx3
    REAL*8, INTENT(IN) :: lambda,mu
    INTEGER, INTENT(IN) :: grn_method
  
    INTEGER :: sx1,sx2,sx3,status

    REAL*4, DIMENSION(:,:), ALLOCATABLE :: p1,p2,p3

    sx1=SIZE(c1,1)-2
    sx2=SIZE(c1,2)
    sx3=SIZE(c1,3)

    ALLOCATE(p1(sx1+2,sx2),p2(sx1+2,sx2),p3(sx1+2,sx2),STAT=status)
    IF (status > 0) THEN
       WRITE_DEBUG_INFO
       WRITE(0,'("could not allocate memory for green function")')
       STOP 1
    ELSE
       p1=0;p2=0;p3=0;
    END IF

    ! forward Fourier transform equivalent body-force
    CALL fft3(c1,sx1,sx2,sx3,dx1,dx2,dx3,FFT_FORWARD)
    CALL fft3(c2,sx1,sx2,sx3,dx1,dx2,dx3,FFT_FORWARD)
    CALL fft3(c3,sx1,sx2,sx3,dx1,dx2,dx3,FFT_FORWARD)
    CALL fft2(t1,sx1,sx2,dx1,dx2,FFT_FORWARD)
    CALL fft2(t2,sx1,sx2,dx1,dx2,FFT_FORWARD)
    CALL fft2(t3,sx1,sx2,dx1,dx2,FFT_FORWARD)
   
    ! solve for displacement field
    CALL elasticresponse(lambda,mu,c1,c2,c3,dx1,dx2,dx3)
    IF (GRN_IMAGE .eq. grn_method) THEN
       CALL surfacenormaltraction(lambda,mu,c1,c2,c3,dx1,dx2,dx3,p3)
       p3=t3-p3
       CALL boussinesq3d(p3,lambda,mu,c1,c2,c3,dx1,dx2,dx3)
    ELSE
       CALL surfacetraction(lambda,mu,c1,c2,c3,dx1,dx2,dx3,p1,p2,p3)
       p1=t1-p1
       p2=t2-p2
       p3=t3-p3
       CALL cerruti3d(p1,p2,p3,lambda,mu,c1,c2,c3,dx1,dx2,dx3)
    END IF

    ! inverse Fourier transform solution displacement components
    CALL fft3(c1,sx1,sx2,sx3,dx1,dx2,dx3,FFT_INVERSE)
    CALL fft3(c2,sx1,sx2,sx3,dx1,dx2,dx3,FFT_INVERSE)
    CALL fft3(c3,sx1,sx2,sx3,dx1,dx2,dx3,FFT_INVERSE)
    CALL fft2(t1,sx1,sx2,dx1,dx2,FFT_INVERSE)
    CALL fft2(t2,sx1,sx2,dx1,dx2,FFT_INVERSE)
    CALL fft2(t3,sx1,sx2,dx1,dx2,FFT_INVERSE)

    DEALLOCATE(p1,p2,p3)
    
  END SUBROUTINE greenfunction

  !------------------------------------------------------------------
  !> subroutine GreensFunctionCowling
  !! computes (inplace) the displacement components due to a set of
  !! 3-D body-forces by application of the semi-analytic Green's
  !! function. The solution satisfies quasi-static Navier's equation
  !! with buoyancy boundary condition to simulate the effect of 
  !! gravity (the Cowling approximation).
  !!
  !! the importance of gravity depends upon the density contrast rho 
  !! at the surface, the acceleration of gravity g and the value of 
  !! shear modulus mu in the crust. effect on the displacement field
  !! is governed by the gradient
  !!
  !!            gamma = (1 - nu) rho g / mu
  !!                  = rho g / (2 mu alpha)
  !! 
  !! where nu is the Poisson's ratio. For a Poisson's solid with 
  !! nu = 1/4, with a density contrast rho = 3200 kg/m^3 and a shear
  !! modulus mu = 30 GPa, we have gamma = 0.8e-6 /m.
  !!
  !! INPUT:
  !!   @param c1,c2,c3    is a set of body forces
  !!   @param dx1,dx2,dx3 are the sampling size
  !!   @param lambda,mu   are the Lame elastic parameters
  !!   @param gamma       is the gravity coefficient
  !
  ! sylvain barbot (07/07/07) - original function greenfunction
  !                (11/21/08) - effect of gravity
  !------------------------------------------------------------------
  SUBROUTINE greenfunctioncowling(c1,c2,c3,t1,t2,t3,dx1,dx2,dx3, &
                                  lambda,mu,gamma)
    REAL*4, INTENT(INOUT), DIMENSION(:,:,:) :: c1,c2,c3
    REAL*4, INTENT(INOUT), DIMENSION(:,:) :: t1,t2,t3
    REAL*8, INTENT(IN) :: dx1,dx2,dx3
    REAL*8, INTENT(IN) :: lambda,mu,gamma
  
    INTEGER :: sx1,sx2,sx3,status,cuStatus
    INTEGER :: direction 

    REAL*4, DIMENSION(:,:), ALLOCATABLE :: p1,p2,p3

#ifdef PAPI_PROF
    CHARACTER (LEN=16) cTimerName
    cTimerName = 'FFT'
#endif

    sx1=SIZE(c1,1)-2
    sx2=SIZE(c1,2)
    sx3=SIZE(c1,3)
    cuStatus = 0

    ALLOCATE(p1(sx1+2,sx2),p2(sx1+2,sx2),p3(sx1+2,sx2),STAT=status)
    IF (status > 0) THEN
       WRITE_DEBUG_INFO
       WRITE(0,'("could not allocate memory for green function")')
       STOP 1
    ELSE
       p1=0;p2=0;p3=0;
    END IF


    ! forward Fourier transform equivalent body-force
#ifdef USING_CUDA 
    
    direction = FFT_FORWARD
    CALL calfftnelasic (c1, c2, c3, t1, t2, t3, %VAL(sx1),%VAL(sx2),%VAL(sx3), &
                         %VAL(dx1),%VAL(dx2),%VAL(dx3),%VAL(lambda),%VAL(mu),%VAL(gamma), p1, p2, p3, cuStatus)
    IF (cuStatus > 0) THEN 
       WRITE_DEBUG_INFO
       WRITE(0,'("Check the logs for the reason")')
       STOP 1
    END IF
#else
     
#ifdef PAPI_PROF
    CALL papistartprofiling(cTimerName) 
#endif
    CALL fft3(c1,sx1,sx2,sx3,dx1,dx2,dx3,FFT_FORWARD)
    CALL fft3(c2,sx1,sx2,sx3,dx1,dx2,dx3,FFT_FORWARD)
    CALL fft3(c3,sx1,sx2,sx3,dx1,dx2,dx3,FFT_FORWARD)
    CALL fft2(t1,sx1,sx2,dx1,dx2,FFT_FORWARD)
    CALL fft2(t2,sx1,sx2,dx1,dx2,FFT_FORWARD)
    CALL fft2(t3,sx1,sx2,dx1,dx2,FFT_FORWARD)

#ifdef PAPI_PROF  
    CALL papiendprofiling(cTimerName)
    cTimerName = 'Elastic'
    CALL papistartprofiling(cTimerName)
#endif


    ! solve for displacement field
    CALL elasticresponse(lambda,mu,c1,c2,c3,dx1,dx2,dx3)

#ifdef PAPI_PROF  
    CALL papiendprofiling(cTimerName)
    cTimerName = 'Surface'
    CALL papistartprofiling(cTimerName)
#endif
 
   CALL surfacetractioncowling(lambda,mu,gamma, &
         c1,c2,c3,dx1,dx2,dx3,p1,p2,p3)
    
    p1=t1-p1
    p2=t2-p2
    p3=t3-p3
#ifdef PAPI_PROF  
    CALL papiendprofiling(cTimerName)
    cTimerName = 'cerruti'
    CALL papistartprofiling(cTimerName)
#endif

#ifdef CERRUTI_FFT
    CALL cerruticowling(p1,p2,p3,lambda,mu,gamma, &
         c1,c2,c3,dx1,dx2,dx3)
#else
    
    CALL cerrutimodified(p1,p2,p3,lambda,mu,gamma, &
         c1,c2,c3,dx1,dx2,dx3)
#endif

#ifdef PAPI_PROF  
    CALL papiendprofiling(cTimerName)
#endif

#endif   
#ifndef USING_CUDA

#ifdef PAPI_PROF  
    cTimerName = 'inversefft'
    CALL papistartprofiling(cTimerName)
#endif

    CALL fft3(c1,sx1,sx2,sx3,dx1,dx2,dx3,FFT_INVERSE)   
    CALL fft3(c2,sx1,sx2,sx3,dx1,dx2,dx3,FFT_INVERSE)
    CALL fft3(c3,sx1,sx2,sx3,dx1,dx2,dx3,FFT_INVERSE)
    CALL fft2(t1,sx1,sx2,dx1,dx2,FFT_INVERSE)
    CALL fft2(t2,sx1,sx2,dx1,dx2,FFT_INVERSE)
    CALL fft2(t3,sx1,sx2,dx1,dx2,FFT_INVERSE)

#ifdef PAPI_PROF  
    CALL papiendprofiling(cTimerName)
#endif


#endif
 

    DEALLOCATE(p1,p2,p3)
    
  END SUBROUTINE greenfunctioncowling

END MODULE green
