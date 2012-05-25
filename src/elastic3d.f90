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

MODULE elastic3d

  USE types
  USE fourier

  IMPLICIT NONE

#include "include.f90"

  REAL*8, PRIVATE, PARAMETER :: pi   = 3.141592653589793115997963468544185161_8
  REAL*8, PRIVATE, PARAMETER :: pi2  = 6.28318530717958623199592693708837032318_8
  REAL*8, PRIVATE, PARAMETER :: pid2 = 1.57079632679489655799898173427209258079_8
  REAL*8, PRIVATE, PARAMETER :: DEG2RAD = 0.01745329251994329547437168059786927_8
    
  INTERFACE OPERATOR (.times.)
     MODULE PROCEDURE tensorscalarprod
  END INTERFACE

  INTERFACE OPERATOR (.minus.)
     MODULE PROCEDURE tensordiff
  END INTERFACE

  INTERFACE OPERATOR (.plus.)
     MODULE PROCEDURE tensorplus
  END INTERFACE

  INTERFACE OPERATOR (.sdyad.)
     MODULE PROCEDURE tensorsymmetricdyadprod
  END INTERFACE

  INTERFACE OPERATOR (.tdot.)
     MODULE PROCEDURE tensorvectordotprod
  END INTERFACE

CONTAINS

  !------------------------------------------------------------
  !> function fix
  !! returns the closest integer scalar
  !
  ! sylvain barbot (08/25/07) - original form
  !------------------------------------------------------------
  INTEGER FUNCTION fix(number)
    REAL*8, INTENT(IN) :: number

    INTEGER :: c,f
    f=FLOOR(number)
    c=CEILING(number)

    IF ((number-f) .gt. 0.5_8) THEN
       fix=c
    ELSE
       fix=f
    END IF

  END FUNCTION fix

  !------------------------------------------------------------
  !> function SINH
  !! computes the hyperbolic sine
  !------------------------------------------------------------
  REAL*8 FUNCTION my_sinh(x)
    REAL*8, INTENT(IN) :: x

    IF (abs(x) .GT. 11._8) THEN
       my_sinh=sign(1._8,x)*sinh(11._8)
    ELSE
       my_sinh=sinh(x)
    END IF
  END FUNCTION my_sinh

  !-----------------------------------------------------------------
  !> subroutine Neighbor
  !! computes the indices of neighbor samples (l points away)
  !! bracketing the current samples location i1,i2,i3 and
  !! assuming periodic boundary condition.
  !!
  !!           i1m < i1 < i1p
  !!           i2m < i2 < i2p
  !!           i3m < i3 < i3p
  !-----------------------------------------------------------------
  SUBROUTINE neighbor(i1,i2,i3,sx1,sx2,sx3,l,i1m,i1p,i2m,i2p,i3m,i3p)
    INTEGER, INTENT(IN) :: i1,i2,i3,sx1,sx2,sx3,l
    INTEGER, INTENT(OUT) :: i1m,i1p,i2m,i2p,i3m,i3p

    i1m=mod(sx1+i1-1-l,sx1)+1
    i1p=mod(i1-1+l,sx1)+1
    i2m=mod(sx2+i2-1-l,sx2)+1
    i2p=mod(i2-1+l,sx2)+1
    i3m=mod(sx3+i3-1-l,sx3)+1
    i3p=mod(i3-1+l,sx3)+1

  END SUBROUTINE neighbor

  !---------------------------------------------------------------
  !> subroutine IsotropicStressStrain
  !! computes in place the isotropic stress tensor from a given
  !! strain tensor using Hooke's law stress-strain relationship.
  !
  ! sylvain barbot (10/14/07) - original form
  !---------------------------------------------------------------
  SUBROUTINE isotropicstressstrain(t,lambda,mu)
    TYPE(TENSOR), INTENT(INOUT) :: t
    REAL*8, INTENT(IN) :: lambda, mu

    REAL*8 :: epskk

    epskk=tensortrace(t)

    t = REAL(2._8*mu) .times. t
    t%s11=REAL(t%s11+lambda*epskk)
    t%s22=REAL(t%s22+lambda*epskk)
    t%s33=REAL(t%s33+lambda*epskk)

  END SUBROUTINE isotropicstressstrain

  !------------------------------------------------------------
  !> function TensorDiff
  !! computes the difference between two tensors: t=t1-t2
  !
  ! sylvain barbot (07/09/07) - original form
  !------------------------------------------------------------
  TYPE(TENSOR) FUNCTION tensordiff(t1,t2)
    TYPE(TENSOR), INTENT(IN) :: t1,t2

    tensordiff=TENSOR(t1%s11-t2%s11, & ! 11
                      t1%s12-t2%s12, & ! 12
                      t1%s13-t2%s13, & ! 13
                      t1%s22-t2%s22, & ! 22
                      t1%s23-t2%s23, & ! 23
                      t1%s33-t2%s33)   ! 33

  END FUNCTION tensordiff

  !------------------------------------------------------------
  !> function TensorPlus
  !! computes the sum of two tensors: t=t1-t2
  !
  ! sylvain barbot (07/09/07) - original form
  !------------------------------------------------------------
  TYPE(TENSOR) FUNCTION tensorplus(t1,t2)
    TYPE(TENSOR), INTENT(IN) :: t1,t2

    tensorplus=TENSOR(t1%s11+t2%s11, & ! 11
                      t1%s12+t2%s12, & ! 12
                      t1%s13+t2%s13, & ! 13
                      t1%s22+t2%s22, & ! 22
                      t1%s23+t2%s23, & ! 23
                      t1%s33+t2%s33)   ! 33

  END FUNCTION tensorplus

  !------------------------------------------------------------
  !> function TensorScalarProd
  !! multiplies a tensor with a scalar
  !
  ! sylvain barbot (07/09/07) - original form
  !------------------------------------------------------------
  TYPE(TENSOR) FUNCTION tensorscalarprod(scalar,t)
    TYPE(TENSOR), INTENT(IN) :: t
    REAL*4, INTENT(IN) :: scalar

    tensorscalarprod=TENSOR(scalar*t%s11, & ! 11
                            scalar*t%s12, & ! 12
                            scalar*t%s13, & ! 13
                            scalar*t%s22, & ! 22
                            scalar*t%s23, & ! 23
                            scalar*t%s33)   ! 33

  END FUNCTION tensorscalarprod

  !------------------------------------------------------------
  !> function TensorSymmetricDyadProd
  !! computes the dyadic product of two vectors to obtain a
  !! symmetric second order tensor
  !
  ! sylvain barbot (07/09/07) - original form
  !------------------------------------------------------------
  TYPE(TENSOR) FUNCTION tensorsymmetricdyadprod(a,b)
    REAL*8, DIMENSION(3), INTENT(IN) :: a,b

    tensorsymmetricdyadprod=TENSOR( &
          REAL(a(1)*b(1)),                 & ! 11
         REAL((a(1)*b(2)+a(2)*b(1))/2._8), & ! 12
         REAL((a(1)*b(3)+a(3)*b(1))/2._8), & ! 13
          REAL(a(2)*b(2)),                 & ! 22
         REAL((a(2)*b(3)+a(3)*b(2))/2._8), & ! 23
          REAL(a(3)*b(3))                  & ! 33
          )

  END FUNCTION tensorsymmetricdyadprod

  !------------------------------------------------------------
  !> function TensorVectorDotProd
  !! compute the dot product T.v where T is a second-order
  !! tensor and v is a vector.
  !
  ! sylvain barbot (07/09/07) - original form
  !------------------------------------------------------------
  FUNCTION tensorvectordotprod(t,v)
    TYPE(TENSOR), INTENT(IN) :: t
    REAL*8, DIMENSION(3), INTENT(IN) :: v
    REAL*8, DIMENSION(3) :: tensorvectordotprod

    tensorvectordotprod= &
         (/ t%s11*v(1)+t%s12*v(2)+t%s13*v(3), &
            t%s12*v(1)+t%s22*v(2)+t%s23*v(3), &
            t%s13*v(1)+t%s23*v(2)+t%s33*v(3) /)

  END FUNCTION tensorvectordotprod

  !------------------------------------------------------------
  !> function TensorVectorDotProd
  !! compute the dot product T.v where T is a second-order
  !! tensor and v is a vector.
  !
  ! sylvain barbot (07/09/07) - original form
  !------------------------------------------------------------
  FUNCTION tensordeviatoric(t)
    TYPE(TENSOR), INTENT(IN) :: t
    TYPE(TENSOR) :: tensordeviatoric

    REAL*4 :: diag

    diag=REAL(tensortrace(t)/3._8)
    
    tensordeviatoric%s11=t%s11-diag
    tensordeviatoric%s12=t%s12
    tensordeviatoric%s13=t%s13
    tensordeviatoric%s22=t%s22-diag
    tensordeviatoric%s23=t%s23
    tensordeviatoric%s33=t%s33-diag

  END FUNCTION tensordeviatoric

  !------------------------------------------------------------
  !> function TensorTrace
  !! computes the trace of a second order tensor
  !
  ! sylvain barbot (07/09/07) - original form
  !------------------------------------------------------------
  REAL*8 FUNCTION tensortrace(t)
    TYPE(TENSOR), INTENT(IN) :: t

    tensortrace=t%s11+t%s22+t%s33

  END FUNCTION tensortrace

  !------------------------------------------------------------
  !> function TensorNorm
  !! computes the Frobenius norm of a second order tensor
  !
  ! sylvain barbot (07/09/07) - original form
  !------------------------------------------------------------
  REAL*8 FUNCTION tensornorm(t)
    TYPE(TENSOR), INTENT(IN) :: t

    tensornorm=SQRT(( &
         t%s11**2+2._8*t%s12**2+2._8*t%s13**2+ &
         t%s22**2+2._8*t%s23**2+ &
         t%s33**2)/2._8)

  END FUNCTION tensornorm

  !------------------------------------------------------------
  !> function TensorDecomposition
  !! writes a tensor t as the product of a norm and a direction
  !!
  !!         t = gamma * R
  !!
  !! where gamma is a scalar, the norm of t, and R is a unitary
  !! tensor. t is assumed to be a deviatoric tensor.
  !
  ! sylvain barbot (07/09/07) - original form
  !------------------------------------------------------------
  SUBROUTINE tensordecomposition(t,gamma,R)
    TYPE(TENSOR), INTENT(IN) :: t
    TYPE(TENSOR), INTENT(OUT) :: R
    REAL*8, INTENT(OUT) :: gamma
    
    gamma=tensornorm(t)

    R%s11=REAL(t%s11/gamma)
    R%s12=REAL(t%s12/gamma)
    R%s13=REAL(t%s13/gamma)
    R%s22=REAL(t%s22/gamma)
    R%s23=REAL(t%s23/gamma)
    R%s33=REAL(t%s33/gamma)

  END SUBROUTINE tensordecomposition


  !------------------------------------------------------------
  !> function TensorForbeniusNorm
  !! computes the Frobenius norm of a second order tensor
  !
  ! sylvain barbot (07/09/07) - original form
  !------------------------------------------------------------
  REAL*8 FUNCTION tensorfrobeniusnorm(t)
    TYPE(TENSOR), INTENT(IN) :: t

    tensorfrobeniusnorm=SQRT( &
         t%s11**2+2._8*t%s12**2+2._8*t%s13**2+ &
         t%s22**2+2._8*t%s23**2+ &
         t%s33**2)

  END FUNCTION tensorfrobeniusnorm

  !------------------------------------------------------------
  !> function VectorFieldNormMax
  !! computes the maximum value of the norm of a vector field
  !------------------------------------------------------------
  SUBROUTINE vectorfieldnormmax(v1,v2,v3,sx1,sx2,sx3,maximum,location)
    INTEGER, INTENT(IN) :: sx1,sx2,sx3
#ifdef ALIGN_DATA
    REAL*4, INTENT(IN), DIMENSION(sx1+2,sx2,sx3) :: v1,v2,v3
#else
    REAL*4, INTENT(IN), DIMENSION(sx1,sx2,sx3) :: v1,v2,v3
#endif
    REAL*8, INTENT(OUT) :: maximum
    INTEGER, INTENT(OUT), DIMENSION(3) :: location
    
    INTEGER :: i1,i2,i3
    REAL*8 :: norm

    maximum=-1._8
    DO i3=1,sx3
       DO i2=1,sx2
          DO i1=1,sx1
             norm=SQRT(v1(i1,i2,i3)**2+v2(i1,i2,i3)**2+v3(i1,i2,i3)**2)
             IF (norm .GT. maximum) THEN
                maximum=norm
                location=(/ i1,i2,i3 /)
             END IF
          END DO
       END DO
    END DO
    
  END SUBROUTINE vectorfieldnormmax

  !------------------------------------------------------------
  !> function TensorMean
  !! computesthe mean of the norm of a tensor field
  !------------------------------------------------------------
  REAL*8 FUNCTION tensormean(t)
    TYPE(TENSOR), INTENT(IN), DIMENSION(:,:,:) :: t
    
    INTEGER :: i1,i2,i3,sx1,sx2,sx3
    sx1=SIZE(t,1)
    sx2=SIZE(t,2)
    sx3=SIZE(t,3)

    DO i3=1,sx3
       DO i2=1,sx2
          DO i1=1,sx1
             tensormean=tensormean+tensornorm(t(i1,i2,i3))
          END DO
       END DO
    END DO
    tensormean=tensormean/DBLE(sx1*sx2*sx3)
    
  END FUNCTION tensormean

  !------------------------------------------------------------
  !> function TensorAmplitude
  !! computes the integral of the norm of a tensor field
  !------------------------------------------------------------
  REAL*8 FUNCTION tensoramplitude(t,dx1,dx2,dx3)
    TYPE(TENSOR), INTENT(IN), DIMENSION(:,:,:) :: t
    REAL*8, INTENT(IN) :: dx1,dx2,dx3

    INTEGER :: i1,i2,i3,sx1,sx2,sx3
    sx1=SIZE(t,1)
    sx2=SIZE(t,2)
    sx3=SIZE(t,3)

    tensoramplitude=0._8
    DO i3=1,sx3
       DO i2=1,sx2
          DO i1=1,sx1
             tensoramplitude=tensoramplitude &
                  +tensornorm(t(i1,i2,i3))
          END DO
       END DO
    END DO
    tensoramplitude=tensoramplitude*DBLE(dx1*dx2*dx3)

  END FUNCTION tensoramplitude

  !------------------------------------------------------------
  !> function TensorMeanTrace
  !! computesthe mean of the norm of a tensor field
  !------------------------------------------------------------
  REAL*8 FUNCTION tensormeantrace(t)
    TYPE(TENSOR), INTENT(IN), DIMENSION(:,:,:) :: t
    
    INTEGER :: i1,i2,i3,sx1,sx2,sx3
    sx1=SIZE(t,1)
    sx2=SIZE(t,2)
    sx3=SIZE(t,3)

    DO i3=1,sx3
       DO i2=1,sx2
          DO i1=1,sx1
             tensormeantrace= &
                  tensormeantrace+tensortrace(t(i1,i2,i3))
          END DO
       END DO
    END DO
    tensormeantrace=tensormeantrace/DBLE(sx1*sx2*sx3)
    
  END FUNCTION tensormeantrace

  !------------------------------------------------------------
  !> sinc function
  !! computes sin(pi*x)/(pi*x)
  !
  ! sylvain barbot (04-14-07) - original form
  !------------------------------------------------------------
  FUNCTION sinc(x)
    REAL*8 :: sinc
    REAL*8, INTENT(IN) :: x
    IF (x /= 0) THEN
       sinc=sin(pi*x)/(pi*x)
    ELSE
       sinc=1._8
    END IF
  END FUNCTION sinc
  
  !-------------------------------------------------------------------------
  !> function gauss computes the normalized gaussian function
  !
  ! Sylvain Barbot (06-29-07)
  !-------------------------------------------------------------------------
  FUNCTION gauss(x,sigma)
    REAL*8 :: gauss
    REAL*8, INTENT(IN) :: x,sigma
    
    gauss=exp(-0.5_8*(x/sigma)**2)/sqrt(pi2)/sigma
  END FUNCTION gauss
  
  !-------------------------------------------------------------------------
  !> function gaussp computes the normalized gaussian derivative
  !
  ! Sylvain Barbot (06-29-07)
  !-------------------------------------------------------------------------
  FUNCTION gaussp(x,sigma)
    REAL*8 :: gaussp
    REAL*8, INTENT(IN) :: x,sigma
    
    gaussp=-x*exp(-0.5_8*(x/sigma)**2)/sqrt(pi2)/sigma**3
  END FUNCTION gaussp

  !-------------------------------------------------------------------------
  !> function omega computes raised-cosine taper in the space domain
  !
  ! Sylvain Barbot (06-29-07)
  !-------------------------------------------------------------------------
  FUNCTION omega(x,beta)
    REAL*8 :: omega
    REAL*8, INTENT(IN) :: x,beta
    
    IF (abs(x) .le. (1._8-2._8*beta)/(1._8-beta)/2._8) THEN
       omega=1._8
    ELSE
       IF (abs(x) .lt. 1._8/(1-beta)/2._8) THEN
          omega=cos(pi*((1._8-beta)*abs(x)-0.5_8+beta)/2._8/beta)**2
       ELSE
          omega=0._8
       END IF
    END IF
  END FUNCTION omega

  !-------------------------------------------------------------------------
  !> function omegap computes raised-cosine taper derivative 
  !! in the space domain
  !
  ! Sylvain Barbot (06-29-07)
  !-------------------------------------------------------------------------
  FUNCTION omegap(x,beta)
    REAL*8 :: omegap
    REAL*8, INTENT(IN) :: x,beta
    
    omegap=0
    IF (abs(x) .gt. (1._8-2._8*beta)/(1._8-beta)/2._8) THEN
       IF (abs(x) .lt. 1._8/(1-beta)/2._8) THEN
          omegap=-DSIGN(1._8,x)*pi*(1._8-beta)/2._8/beta* &
               sin(pi*((1._8-beta)*abs(x)-0.5_8+beta)/beta)
       END IF
    END IF
  END FUNCTION omegap
  
  !-------------------------------------------------------------------------
  !> tapered step function (raised-cosine) of unit area in the Fourier domain
  !!
  !! INPUT
  !! @param k        wavenumber
  !! @param beta     roll-off parameter 0<beta<0.5
  !!                 no smoothing for beta close to 0
  !!                 string smoothing for beta close to 0.5
  !
  ! sylvain barbot (04-14-07) - original form
  !-------------------------------------------------------------------------
  FUNCTION omegak(k,beta)
    REAL*8 :: omegak
    REAL*8, INTENT(IN) :: k, beta
    REAL*8 :: gamma,denom,om1,om2
    
    gamma=(1._8-beta)
    denom=(gamma-(4._8*beta**2._8/gamma)*k**2._8)*2._8
    om1=sinc(k/gamma)
    om2=(1._8-2._8*beta)*sinc(((1._8-2._8*beta)/gamma)*k)
    omegak=(om1+om2)/denom

  END FUNCTION omegak

  !----------------------------------------------------------------
  !> subroutine TensorStructure
  !! constructs a vertically-stratified tensor field.
  !! The structure is defined by its interfaces: changes can be
  !! gradual or discontinuous.
  !
  ! sylvain barbot (10/25/08) - original form
  !----------------------------------------------------------------
  SUBROUTINE tensorstructure(vstruct,layers,dx3)
    TYPE(TENSOR_LAYER_STRUCT), INTENT(IN), DIMENSION(:) :: layers
    TYPE(TENSOR_LAYER_STRUCT), INTENT(OUT), DIMENSION(:) :: vstruct
    REAL*8, INTENT(IN) :: dx3

    INTEGER :: nv,k,i3s,i3e=1,i3,sx3
    REAL*8 :: z,z0,z1
    TYPE(TENSOR) :: t0,t1,t
         
    nv =SIZE(layers,1)
    sx3=SIZE(vstruct,1)

    IF (0 .ge. nv) THEN

       WRITE (0,'("error at line ",I5.5," of source file ",a)') __LINE__,__FILE__
       WRITE_DEBUG_INFO
       WRITE (0,'("invalid tensor structure. exiting.")')
       STOP 1
    END IF

    ! initialization
    vstruct(:)%z=0      ! depth is not used
    vstruct(:)%t=tensor(0._4,0._4,0._4,0._4,0._4,0._4) ! default

    z0=fix(layers(1)%z/dx3)*dx3
    DO k=1,nv
       ! project model on multiples of sampling size 'dx3'
       ! to avoid aliasing problems
       z1=fix(layers(k)%z/dx3)*dx3

       IF (z1 .lt. z0) THEN
          WRITE_DEBUG_INFO
          WRITE (0,'("invalid mechanical structure.")')
          WRITE (0,'("depths must be increasing. exiting.")')
          STOP 1
       END IF

       IF (z1 .eq. z0) THEN
          ! discontinuous interface in the elastic structure
          z0=z1
          
          t1=layers(k)%t
          
          i3e=fix(z1/dx3+1)
       ELSE
          ! interpolate linearly between current and previous value

          t1=layers(k)%t

          i3s=fix(z0/dx3)+1
          i3e=MIN(fix(z1/dx3+1),sx3)
          DO i3=i3s,i3e
             z=(i3-1._8)*dx3

             t=REAL(1._8/(z1-z0)) .times. &
                  ((REAL(z-z0) .times. t1) .plus. (REAL(z1-z) .times. t0))
             
             vstruct(i3)%t=t
 
         END DO
       END IF

       z0=z1
       t0=t1

    END DO

    ! downward-continue the last layer
    IF (fix(z1/dx3) .lt. sx3-1) THEN
       vstruct(i3e:sx3)%t=t1
    END IF

  END SUBROUTINE tensorstructure


  !----------------------------------------------------------------
  !> subroutine ViscoElasticStructure
  !! constructs a vertically-stratified viscoelastic structure.
  !! The structure is defined by its interfaces: changes can be
  !! gradual or discontinuous.
  !!
  !! EXAMPLE INPUTS:
  !!
  !! 1- elastic plate over linear viscous half-space
  !!    1
  !!    1 1.0 1.0 1.0
  !!
  !! 2- elastic plate over powerlaw viscous half-space (n=3)
  !!    1
  !!    1 1.0 1.0 3.0
  !!
  !! 3- elastic plate over viscous half-space with depth-dependent
  !!    viscosity
  !!    2
  !!    1 01.0 1.0 1.0
  !!    2 10.0 6.0 1.0
  !!
  !!    in this last example, the grid does not have to reach down
  !!    to x3=10.
  !!
  !! \author sylvain barbot (08/07/07) - original form
  !----------------------------------------------------------------
  SUBROUTINE viscoelasticstructure(vstruct,layers,dx3)
    TYPE(LAYER_STRUCT), INTENT(IN), DIMENSION(:) :: layers
    TYPE(LAYER_STRUCT), INTENT(OUT), DIMENSION(:) :: vstruct
    REAL*8, INTENT(IN) :: dx3

    INTEGER :: nv,k,i3s,i3e=1,i3,sx3
    REAL*8 :: z,z0,z1, &
         power,power0,power1, &
         gamma,gamma0,gamma1, &
         friction,friction0,friction1, &
         cohesion,cohesion0,cohesion1
         

    nv =SIZE(layers,1)
    sx3=SIZE(vstruct,1)

    IF (0 .ge. nv) THEN
       WRITE_DEBUG_INFO
       WRITE (0,'("invalid elastic structure. exiting.")')
       STOP 1
    END IF

    ! initialization
    vstruct(:)%z=0      ! depth is not used
    vstruct(:)%gammadot0=0 ! default is inviscid
    vstruct(:)%friction=0.6  ! default is friction=0.6
    vstruct(:)%cohesion=0  ! default is no cohesion
    vstruct(:)%stressexponent=layers(1)%stressexponent  ! default

    z0=fix(layers(1)%z/dx3)*dx3
    DO k=1,nv
       ! project model on multiples of sampling size 'dx3'
       ! to avoid aliasing problems
       z1=fix(layers(k)%z/dx3)*dx3

       IF (z1 .lt. z0) THEN
          WRITE_DEBUG_INFO
          WRITE (0,'("invalid mechanical structure. exiting.")')
          STOP 1
       END IF

       IF (z1 .eq. z0) THEN
          ! discontinuous interface in the elastic structure
          z0=z1
          gamma1=layers(k)%gammadot0
          power1 =layers(k)%stressexponent
          friction1=layers(k)%friction
          cohesion1=layers(k)%cohesion
          
          i3e=fix(z1/dx3+1)
       ELSE
          ! interpolate between current and previous value
          gamma1=layers(k)%gammadot0
          power1 =layers(k)%stressexponent
          friction1=layers(k)%friction
          cohesion1=layers(k)%cohesion

          i3s=fix(z0/dx3)+1
          i3e=MIN(fix(z1/dx3+1),sx3)
          DO i3=i3s,i3e
             z=(i3-1._8)*dx3
             gamma=((z-z0)*gamma1+(z1-z)*gamma0)/(z1-z0)
             power=((z-z0)*power1+(z1-z)*power0)/(z1-z0)
             friction=((z-z0)*friction1+(z1-z)*friction0)/(z1-z0)
             cohesion=((z-z0)*cohesion1+(z1-z)*cohesion0)/(z1-z0)

             vstruct(i3)%gammadot0=gamma
             vstruct(i3)%stressexponent =power
             vstruct(i3)%friction=friction
             vstruct(i3)%cohesion=cohesion
          END DO
       END IF

       z0=z1
       gamma0=gamma1
       power0=power1
       friction0=friction1
       cohesion0=cohesion1

    END DO

    ! downward-continue the last layer
    IF (fix(z1/dx3) .lt. sx3-1) THEN
       vstruct(i3e:sx3)%gammadot0=REAL(gamma1)
       vstruct(i3e:sx3)%stressexponent =REAL(power1)
       vstruct(i3e:sx3)%friction=REAL(friction1)
       vstruct(i3e:sx3)%cohesion=REAL(cohesion1)
    END IF

  END SUBROUTINE viscoelasticstructure


  !------------------------------------------------------------------
  !> function OptimalFilter
  !! load predefined Finite Impulse Response (FIR) filters of various
  !! lengths and select the most appropriate ones based on the
  !! computational grid size. result is filter kernels always smaller
  !! than available computational length.
  !! this is useful in the special cases of infinite faults where
  !! deformation is essentially two-dimensional, despite the actual
  !! three-dimensional computation. in the direction of symmetry,
  !! no strain occurs and high accuracy derivative estimates are not
  !! needed.
  !
  ! Sylvain Barbot (03/05/08) - original form
  !------------------------------------------------------------------
  SUBROUTINE optimalfilter(ker1,ker2,ker3,len1,len2,len3,sx1,sx2,sx3)
    REAL*8, DIMENSION(16), INTENT(OUT) :: ker1,ker2,ker3
    INTEGER, INTENT(OUT) :: len1,len2,len3
    INTEGER, INTENT(IN) :: sx1,sx2,sx3

    ! load FIR differentiator filter
    ! variables 'fir1', 'fir7', 'fir14'
    INCLUDE 'kernel1.inc'
    INCLUDE 'kernel7.inc'
    INCLUDE 'kernel14bis.inc'

    ! choose best differentiator kernels
    SELECT CASE(sx1)
    CASE (2:4)
       ! use centered finite difference
       len1=1
       ker1(1)=fir1(1)
    CASE (5:14)
       len1=7
       ker1(1:len1)=fir7(1:len1)
    CASE (15:)
       len1=1
       ker1(1:len1)=fir1(1:len1)
    CASE DEFAULT
       WRITE_DEBUG_INFO
       WRITE (0,'("optimalfilter: invalid dimension. exiting.")')
       STOP 2
    END SELECT

    ! choose best differentiator kernels
    SELECT CASE(sx2)
    CASE (2:4)
       ! use centered finite difference
       len2=1
       ker2(1)=fir1(1)
    CASE (5:14)
       len2=7
       ker2(1:len2)=fir7(1:len2)
    CASE (15:)
       len2=1
       ker2(1:len2)=fir1(1:len2)
    CASE DEFAULT
       WRITE_DEBUG_INFO
       WRITE (0,'("optimalfilter: invalid dimension. exiting.")')
       STOP 2
    END SELECT

    ! choose best differentiator kernels
    SELECT CASE(sx3)
    CASE (5:14)
       len3=7
       ker3(1:len3)=fir7(1:len3)
    CASE (15:)
       len3=1
       ker3(1:len3)=fir1(1:len3)
    CASE DEFAULT
       WRITE_DEBUG_INFO
       WRITE (0,'("optimalfilter: invalid dimension. exiting.")')
       STOP 2
    END SELECT

  END SUBROUTINE optimalfilter

  !-----------------------------------------------------------------
  !> subroutine StressUpdate
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
  !!\verbatim
  !! firpm(14, ...
  !!    [0 7.0e-1 8.000000e-1 8.500000e-1 9.000000e-1 1.0e+0],...
  !!    [0 7.0e-1 5.459372e-1 3.825260e-1 2.433534e-1 0.0e+0]*pi,...
  !!    'differentiator');
  !!\endverbatim
  !!
  !! The kernel is odd and antisymmetric and only half the numbers
  !! are stored in this code. Kernels of different sizes are readilly
  !! available in the 'kernelX.inc' files. Stress tensor field is
  !! obtained by application of Hooke's law
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
  !! \author sylvain barbot (10/10/07) - original form
  !!                                   - optional sample of normal stress
  !!                        (02/12/09) - OpemMP parallel implementation
  !-----------------------------------------------------------------
  SUBROUTINE stressupdate(v1,v2,v3,lambda,mu,dx1,dx2,dx3,sx1,sx2,sx3,sig)
    INTEGER, INTENT(IN) :: sx1,sx2,sx3
    REAL*8, INTENT(IN) :: dx1,dx2,dx3,lambda,mu
#ifdef ALIGN_DATA
    REAL*4, INTENT(IN), DIMENSION(sx1+2,sx2,sx3) :: v1,v2,v3
#else
    REAL*4, INTENT(IN), DIMENSION(sx1,sx2,sx3) :: v1,v2,v3
#endif
    TYPE(TENSOR), INTENT(INOUT), DIMENSION(sx1,sx2,sx3) :: sig

    TYPE(TENSOR) :: t
    INTEGER :: i1,i2,i3,i3p,i3m,len1,len2,len3
    REAL*8 :: px3
    REAL*8, DIMENSION(16) :: ker1,ker2,ker3

    ! load FIR differentiator filter
    CALL optimalfilter(ker1,ker2,ker3,len1,len2,len3,sx1,sx2,sx3)
    ker1=ker1/dx1; ker2=ker2/dx2; ker3=ker3/dx3;

    ! no periodicity in the 3rd direction
    ! use a simple finite difference scheme
    DO i3=1,sx3

       IF ((i3 .gt. len3) .and. (i3 .lt. (sx3-len3+1))) &
            CYCLE

       IF (i3 .eq. 1) THEN
          ! right-centered finite difference
          px3=dx3; i3p=2; i3m=1
       ELSE
          IF (i3 .eq. sx3) THEN
             ! left-centered finite difference
             px3=dx3; i3p=sx3; i3m=sx3-1
          ELSE
             ! centered finite difference
             px3=dx3*2._8; i3m=i3-1; i3p=i3+1
          END IF
       END IF

       DO i2=1,sx2
          DO i1=1,sx1
             CALL localstrain_ani(t,i3m,i3p,px3)
             CALL isotropicstressstrain(t,lambda,mu)
             sig(i1,i2,i3)=sig(i1,i2,i3) .plus. t
          END DO
       END DO
    END DO

    ! intermediate depth treated isotropically
!$omp parallel do private(i1,i2,t)
    DO i3=len3+1,sx3-len3
       DO i2=1,sx2
          DO i1=1,sx1
             ! Finite Impulse Response filter
             !CALL localstrain_fir(t)
             CALL localstrain_fir2(t,i1,i2,i3,ker1,ker2,ker3,len1,len2,len3,v1,v2,v3,sx1,sx2,sx3)
             CALL isotropicstressstrain(t,lambda,mu)
             sig(i1,i2,i3)=sig(i1,i2,i3) .plus. t
          END DO
       END DO
    END DO
!$omp end parallel do

  CONTAINS

    !---------------------------------------------------------------
    !> LocalStrain_FIR2
    !! implements a finite impulse response filter (FIR) to estimate
    !! derivatives and strain components. the compatibility with the
    !! OpenMP parallel execution requires that all variable be 
    !! tractable from the calling routine.
    !!
    !! \author sylvain barbot (10/10/07) - original form
    !                (03/05/08) - implements 3 filters
    !                (02/12/09) - compatibility with OpenMP (scope)
    !---------------------------------------------------------------
    SUBROUTINE localstrain_fir2(e,i1,i2,i3,ker1,ker2,ker3,len1,len2,len3,v1,v2,v3,sx1,sx2,sx3)
      TYPE(TENSOR), INTENT(OUT) :: e
      INTEGER, INTENT(IN) :: len1,len2,len3,i1,i2,i3,sx1,sx2,sx3
      REAL*8, INTENT(IN), DIMENSION(len1) :: ker1
      REAL*8, INTENT(IN), DIMENSION(len2) :: ker2
      REAL*8, INTENT(IN), DIMENSION(len3) :: ker3
      REAL*4, INTENT(IN), DIMENSION(:,:,:) :: v1,v2,v3

      INTEGER :: l,i1m,i2m,i3m,i1p,i2p,i3p

      e=TENSOR(0._4,0._4,0._4,0._4,0._4,0._4)

      DO l=1,len1
         ! neighbor samples with periodic boundary conditions
         i1m=mod(sx1+i1-1-l,sx1)+1
         i1p=mod(i1-1+l,sx1)+1

         e%s11=REAL(e%s11+(v1(i1p,i2,i3)-v1(i1m,i2,i3))*ker1(l))
         e%s12=REAL(e%s12+(v2(i1p,i2,i3)-v2(i1m,i2,i3))*ker1(l))
         e%s13=REAL(e%s13+(v3(i1p,i2,i3)-v3(i1m,i2,i3))*ker1(l))
      END DO

      DO l=1,len2
         ! neighbor samples with periodic boundary conditions
         i2m=mod(sx2+i2-1-l,sx2)+1
         i2p=mod(i2-1+l,sx2)+1

         e%s12=REAL(e%s12+(v1(i1,i2p,i3)-v1(i1,i2m,i3))*ker2(l))
         e%s22=REAL(e%s22+(v2(i1,i2p,i3)-v2(i1,i2m,i3))*ker2(l))
         e%s23=REAL(e%s23+(v3(i1,i2p,i3)-v3(i1,i2m,i3))*ker2(l))
      END DO

      DO l=1,len3
         ! neighbor samples in semi-infinite solid
         i3m=i3-l
         i3p=i3+l
         
         e%s13=REAL(e%s13+(v1(i1,i2,i3p)-v1(i1,i2,i3m))*ker3(l))
         e%s23=REAL(e%s23+(v2(i1,i2,i3p)-v2(i1,i2,i3m))*ker3(l))
         e%s33=REAL(e%s33+(v3(i1,i2,i3p)-v3(i1,i2,i3m))*ker3(l))
      END DO
      
      e%s12=REAL(e%s12/2._8)
      e%s13=REAL(e%s13/2._8)
      e%s23=REAL(e%s23/2._8)
      
    END SUBROUTINE localstrain_fir2

    !---------------------------------------------------------------
    !> LocalStrain_FIR
    !! implements a finite impulse response filter (FIR) to estimate
    !! derivatives and strain components.
    !!
    !! \author sylvain barbot (10/10/07) - original form
    !!                        (03/05/08) - implements 3 filters
    !---------------------------------------------------------------
    SUBROUTINE localstrain_fir(e)
      TYPE(TENSOR), INTENT(OUT) :: e

      INTEGER :: l,i1m,i2m,i3m,i1p,i2p,i3p

      e=TENSOR(0._4,0._4,0._4,0._4,0._4,0._4)

      DO l=1,len1
         ! neighbor samples with periodic boundary conditions
         i1m=mod(sx1+i1-1-l,sx1)+1
         i1p=mod(i1-1+l,sx1)+1

         e%s11=REAL(e%s11+(v1(i1p,i2,i3)-v1(i1m,i2,i3))*ker1(l))
         e%s12=REAL(e%s12+(v2(i1p,i2,i3)-v2(i1m,i2,i3))*ker1(l))
         e%s13=REAL(e%s13+(v3(i1p,i2,i3)-v3(i1m,i2,i3))*ker1(l))
      END DO

      DO l=1,len2
         ! neighbor samples with periodic boundary conditions
         i2m=mod(sx2+i2-1-l,sx2)+1
         i2p=mod(i2-1+l,sx2)+1

         e%s12=REAL(e%s12+(v1(i1,i2p,i3)-v1(i1,i2m,i3))*ker2(l))
         e%s22=REAL(e%s22+(v2(i1,i2p,i3)-v2(i1,i2m,i3))*ker2(l))
         e%s23=REAL(e%s23+(v3(i1,i2p,i3)-v3(i1,i2m,i3))*ker2(l))
      END DO

      DO l=1,len3
         ! neighbor samples in semi-infinite solid
         i3m=i3-l
         i3p=i3+l

         e%s13=REAL(e%s13+(v1(i1,i2,i3p)-v1(i1,i2,i3m))*ker3(l))
         e%s23=REAL(e%s23+(v2(i1,i2,i3p)-v2(i1,i2,i3m))*ker3(l))
         e%s33=REAL(e%s33+(v3(i1,i2,i3p)-v3(i1,i2,i3m))*ker3(l))
      END DO

      e%s12=REAL(e%s12/2._8)
      e%s13=REAL(e%s13/2._8)
      e%s23=REAL(e%s23/2._8)

    END SUBROUTINE localstrain_fir

    !---------------------------------------------------------------
    !> LocalStrain_ANI
    !! implements a different finite impulse response filter (FIR)
    !! in each direction (ANIsotropy) to estimate derivatives and
    !! strain components.
    !
    ! sylvain barbot (10/10/07) - original form
    !                (03/05/09) - implements 3 filters
    !---------------------------------------------------------------
    SUBROUTINE localstrain_ani(e,i3m,i3p,px3)
      TYPE(TENSOR), INTENT(OUT) :: e
      INTEGER, INTENT(IN) :: i3m, i3p
      REAL*8, INTENT(IN) :: px3

      INTEGER :: l,i1m,i2m,i1p,i2p

      e=TENSOR(0._4,0._4,0._4,0._4,0._4,0._4)

      DO l=1,len1
         ! neighbor samples with periodic boundary conditions
         i1m=mod(sx1+i1-1-l,sx1)+1
         i1p=mod(i1-1+l,sx1)+1

         e%s11=REAL(e%s11+(v1(i1p,i2,i3)-v1(i1m,i2,i3))*ker1(l))
         e%s12=REAL(e%s12+(v2(i1p,i2,i3)-v2(i1m,i2,i3))*ker1(l))
         e%s13=REAL(e%s13+(v3(i1p,i2,i3)-v3(i1m,i2,i3))*ker1(l))
      END DO

      DO l=1,len2
         ! neighbor samples with periodic boundary conditions
         i2m=mod(sx2+i2-1-l,sx2)+1
         i2p=mod(i2-1+l,sx2)+1

         e%s12=REAL(e%s12+(v1(i1,i2p,i3)-v1(i1,i2m,i3))*ker2(l))
         e%s22=REAL(e%s22+(v2(i1,i2p,i3)-v2(i1,i2m,i3))*ker2(l))
         e%s23=REAL(e%s23+(v3(i1,i2p,i3)-v3(i1,i2m,i3))*ker2(l))
      END DO

      ! finite difference in the 3rd direction
      e%s13=REAL(e%s13 + (v1(i1,i2,i3p)-v1(i1,i2,i3m))/px3)
      e%s23=REAL(e%s23 + (v2(i1,i2,i3p)-v2(i1,i2,i3m))/px3)
      e%s33=REAL((v3(i1,i2,i3p)-v3(i1,i2,i3m))/px3)

      e%s12=REAL(e%s12/2._8)
      e%s13=REAL(e%s13/2._8)
      e%s23=REAL(e%s23/2._8)

    END SUBROUTINE localstrain_ani

  END SUBROUTINE stressupdate

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
  !!
  !! \author sylvain barbot (07/09/07) - original form
  !!                        (10/09/07) - upgrade the finite difference scheme
  !!                                     to a finite impulse response filter
  !!                        (02/12/09) - OpenMP parallel implementation
  !-----------------------------------------------------------------
  SUBROUTINE equivalentbodyforce(sig,dx1,dx2,dx3,sx1,sx2,sx3, &
                                 c1,c2,c3,t1,t2,t3,mask)
    REAL*8, INTENT(IN) :: dx1,dx2,dx3
    INTEGER, INTENT(IN) :: sx1,sx2,sx3
#ifdef ALIGN_DATA
    REAL*4, INTENT(INOUT), DIMENSION(sx1+2,sx2,sx3) :: c1,c2,c3
    REAL*4, INTENT(INOUT), DIMENSION(sx1+2,sx2) :: t1,t2,t3
#else
    REAL*4, INTENT(INOUT), DIMENSION(sx1,sx2,sx3) :: c1,c2,c3
    REAL*4, INTENT(INOUT), DIMENSION(sx1,sx2) :: t1,t2,t3
#endif
    TYPE(TENSOR), INTENT(IN), DIMENSION(sx1,sx2,sx3) :: sig
    REAL*4, INTENT(IN), DIMENSION(sx3), OPTIONAL :: mask

    INTEGER :: i1,i2,i3,i3m,i3p,len1,len2,len3
    REAL*8 :: f1,f2,f3,px3
    REAL*8, DIMENSION(16) :: ker1,ker2,ker3

    CALL optimalfilter(ker1,ker2,ker3,len1,len2,len3,sx1,sx2,sx3)
    ker1=ker1/dx1; ker2=ker2/dx2; ker3=ker3/dx3

    ! equivalent surface traction
    DO i2=1,sx2
       DO i1=1,sx1
          t1(i1,i2)=t1(i1,i2)+sig(i1,i2,1)%s13
          t2(i1,i2)=t2(i1,i2)+sig(i1,i2,1)%s23
          t3(i1,i2)=t3(i1,i2)+sig(i1,i2,1)%s33
       END DO
    END DO

    ! no periodicity in the 3rd direction
    ! use a simple finite difference scheme in the 3rd direction
!$omp parallel 
!$omp do private(i1,i2,f1,f2,f3,px3,i3m,i3p)
    DO i3=1,sx3

       IF ((i3 .gt. len3) .and. (i3 .lt. (sx3-len3+1))) &
            CYCLE

       IF (PRESENT(mask)) THEN
          IF (mask(i3) .EQ. 0) THEN
             CYCLE
          END IF
       END IF

       IF (i3 .eq. 1) THEN
          ! right-centered finite difference
          px3=dx3; i3p=2; i3m=1
       ELSE
          IF (i3 .eq. sx3) THEN
             ! left-centered finite difference
             px3=dx3; i3p=sx3; i3m=sx3-1
          ELSE
             ! centered finite difference
             px3=dx3*2._8; i3m=i3-1; i3p=i3+1
          END IF
       END IF

       DO i2=1,sx2
          DO i1=1,sx1
             CALL localdivergence_ani(f1,f2,f3,i3m,i3p,px3, &
                       i1,i2,i3,ker1,ker2,ker3,len1,len2,len3,sig,sx1,sx2,sx3)

             c1(i1,i2,i3)=c1(i1,i2,i3)-REAL(f1)
             c2(i1,i2,i3)=c2(i1,i2,i3)-REAL(f2)
             c3(i1,i2,i3)=c3(i1,i2,i3)-REAL(f3)

          END DO
       END DO
    END DO
!$omp end do nowait

    ! intermediate depth treated isotropically
!$omp do private(i1,i2,f1,f2,f3)
    DO i3=len3+1,sx3-len3
       
       IF (PRESENT(mask)) THEN
          IF (mask(i3) .EQ. 0) THEN
             CYCLE
          END IF
       END IF
       
       DO i2=1,sx2
          DO i1=1,sx1
             ! Finite Impulse Response filter
             !CALL localdivergence_fir(f1,f2,f3)
             CALL localdivergence_fir2(f1,f2,f3,i1,i2,i3,ker1,ker2,ker3,len1,len2,len3,sig,sx1,sx2,sx3)

             c1(i1,i2,i3)=c1(i1,i2,i3)-REAL(f1)
             c2(i1,i2,i3)=c2(i1,i2,i3)-REAL(f2)
             c3(i1,i2,i3)=c3(i1,i2,i3)-REAL(f3)
          END DO
       END DO
    END DO
!$omp end do
!$omp end parallel

  CONTAINS

    !---------------------------------------------------------------
    ! LocalDivergence_FIR
    ! implements a finite impulse response filter (FIR) to estimate
    ! the divergence of second-order tensor.
    !
    ! ATTENTION - calls to this routine can cause memory leak.
    !
    ! sylvain barbot (10/10/07) - original form
    !                (03/05/08) - implements 3 filters
    !                (02/11/09) - compatibility with OpenMP
    !---------------------------------------------------------------
    SUBROUTINE localdivergence_fir2(f1,f2,f3,i1,i2,i3,ker1,ker2,ker3,len1,len2,len3,sig,sx1,sx2,sx3)
      REAL*8, INTENT(OUT) :: f1,f2,f3
      INTEGER, INTENT(IN) :: len1,len2,len3,i1,i2,i3,sx1,sx2,sx3
      REAL*8, INTENT(IN), DIMENSION(len1) :: ker1
      REAL*8, INTENT(IN), DIMENSION(len2) :: ker2
      REAL*8, INTENT(IN), DIMENSION(len3) :: ker3
      TYPE(TENSOR), INTENT(IN), DIMENSION(:,:,:) :: sig

      INTEGER :: l,i1m,i1p,i2m,i2p,i3m,i3p

      f1=0._8; f2=0._8; f3=0._8
      
      DO l=1,len1
         ! neighbor samples with periodic boundary conditions
         i1m=mod(sx1+i1-1-l,sx1)+1
         i1p=mod(i1-1+l,sx1)+1
         
         f1=f1+(sig(i1p,i2,i3)%s11-sig(i1m,i2,i3)%s11)*ker1(l)
         f2=f2+(sig(i1p,i2,i3)%s12-sig(i1m,i2,i3)%s12)*ker1(l)
         f3=f3+(sig(i1p,i2,i3)%s13-sig(i1m,i2,i3)%s13)*ker1(l)
      END DO
      
      DO l=1,len2
         ! neighbor samples with periodic boundary conditions
         i2m=mod(sx2+i2-1-l,sx2)+1
         i2p=mod(i2-1+l,sx2)+1
         
         f1=f1+(sig(i1,i2p,i3)%s12-sig(i1,i2m,i3)%s12)*ker2(l)
         f2=f2+(sig(i1,i2p,i3)%s22-sig(i1,i2m,i3)%s22)*ker2(l)
         f3=f3+(sig(i1,i2p,i3)%s23-sig(i1,i2m,i3)%s23)*ker2(l)
      END DO
      
      DO l=1,len3
         ! neighbor samples in semi-infinite solid
         i3m=i3-l
         i3p=i3+l
         
         f1=f1+(sig(i1,i2,i3p)%s13-sig(i1,i2,i3m)%s13)*ker3(l)
         f2=f2+(sig(i1,i2,i3p)%s23-sig(i1,i2,i3m)%s23)*ker3(l)
         f3=f3+(sig(i1,i2,i3p)%s33-sig(i1,i2,i3m)%s33)*ker3(l)
      END DO
      
    END SUBROUTINE localdivergence_fir2

    !---------------------------------------------------------------
    ! LocalDivergence_FIR
    ! implements a finite impulse response filter (FIR) to estimate
    ! the divergence of second-order tensor.
    !
    ! ATTENTION - calls to this routine can cause memory leak.
    !
    ! sylvain barbot (10/10/07) - original form
    !                (03/05/08) - implements 3 filters
    !---------------------------------------------------------------
    SUBROUTINE localdivergence_fir(f1,f2,f3)
      REAL*8, INTENT(OUT) :: f1,f2,f3

      INTEGER :: l,i1m,i1p,i2m,i2p,i3m,i3p

      f1=0._8; f2=0._8; f3=0._8

      DO l=1,len1
         ! neighbor samples with periodic boundary conditions
         i1m=mod(sx1+i1-1-l,sx1)+1
         i1p=mod(i1-1+l,sx1)+1

         f1=f1+(sig(i1p,i2,i3)%s11-sig(i1m,i2,i3)%s11)*ker1(l)
         f2=f2+(sig(i1p,i2,i3)%s12-sig(i1m,i2,i3)%s12)*ker1(l)
         f3=f3+(sig(i1p,i2,i3)%s13-sig(i1m,i2,i3)%s13)*ker1(l)
      END DO

      DO l=1,len2
         ! neighbor samples with periodic boundary conditions
         i2m=mod(sx2+i2-1-l,sx2)+1
         i2p=mod(i2-1+l,sx2)+1

         f1=f1+(sig(i1,i2p,i3)%s12-sig(i1,i2m,i3)%s12)*ker2(l)
         f2=f2+(sig(i1,i2p,i3)%s22-sig(i1,i2m,i3)%s22)*ker2(l)
         f3=f3+(sig(i1,i2p,i3)%s23-sig(i1,i2m,i3)%s23)*ker2(l)
      END DO

      DO l=1,len3
         ! neighbor samples in semi-infinite solid
         i3m=i3-l
         i3p=i3+l

         f1=f1+(sig(i1,i2,i3p)%s13-sig(i1,i2,i3m)%s13)*ker3(l)
         f2=f2+(sig(i1,i2,i3p)%s23-sig(i1,i2,i3m)%s23)*ker3(l)
         f3=f3+(sig(i1,i2,i3p)%s33-sig(i1,i2,i3m)%s33)*ker3(l)
      END DO

    END SUBROUTINE localdivergence_fir

    !---------------------------------------------------------------
    ! LocalDivergence_ANI
    ! implements a finite impulse response filter (FIR) in the
    ! horizontal direction and a finite-difference scheme in the
    ! vertical direction to estimate the divergence of second-order
    ! tensor.
    ! Finite difference scheme is left-centered, right-centered or
    ! symmetric, depending on input positions (i3m,i3p) and spacing
    ! (px3).
    !
    ! sylvain barbot (10/10/07) - original form
    !                (03/05/08) - implements 3 filters
    !                (02/12/09) - compatibility with OpenMP
    !---------------------------------------------------------------
    SUBROUTINE localdivergence_ani(f1,f2,f3,i3m,i3p,px3,i1,i2,i3,ker1,ker2,ker3,len1,len2,len3,sig,sx1,sx2,sx3)
      REAL*8, INTENT(OUT) :: f1,f2,f3
      INTEGER, INTENT(IN) :: i3m,i3p,i1,i2,i3,len1,len2,len3,sx1,sx2,sx3
      REAL*8, INTENT(IN), DIMENSION(len1) :: ker1
      REAL*8, INTENT(IN), DIMENSION(len2) :: ker2
      REAL*8, INTENT(IN), DIMENSION(len3) :: ker3
      REAL*8, INTENT(IN) :: px3
      TYPE(TENSOR), INTENT(IN), DIMENSION(:,:,:) :: sig

      INTEGER :: l,i1m,i1p,i2m,i2p

      f1=0._8; f2=0._8; f3=0._8

      ! differentiator filter in the horizontal direction
      DO l=1,len1
         ! neighbor samples with periodic boundary conditions
         i1m=mod(sx1+i1-1-l,sx1)+1
         i1p=mod(i1-1+l,sx1)+1

         f1=f1+(sig(i1p,i2,i3)%s11-sig(i1m,i2,i3)%s11)*ker1(l)
         f2=f2+(sig(i1p,i2,i3)%s12-sig(i1m,i2,i3)%s12)*ker1(l)
         f3=f3+(sig(i1p,i2,i3)%s13-sig(i1m,i2,i3)%s13)*ker1(l)
      END DO

      DO l=1,len2
         ! neighbor samples with periodic boundary conditions
         i2m=mod(sx2+i2-1-l,sx2)+1
         i2p=mod(i2-1+l,sx2)+1

         f1=f1+(sig(i1,i2p,i3)%s12-sig(i1,i2m,i3)%s12)*ker2(l)
         f2=f2+(sig(i1,i2p,i3)%s22-sig(i1,i2m,i3)%s22)*ker2(l)
         f3=f3+(sig(i1,i2p,i3)%s23-sig(i1,i2m,i3)%s23)*ker2(l)
      END DO

      ! finite difference in the 3-direction
      f1=f1+( sig(i1,i2,i3p)%s13-sig(i1,i2,i3m)%s13 )/px3
      f2=f2+( sig(i1,i2,i3p)%s23-sig(i1,i2,i3m)%s23 )/px3
      f3=f3+( sig(i1,i2,i3p)%s33-sig(i1,i2,i3m)%s33 )/px3

    END SUBROUTINE localdivergence_ani

    !-------------------------------------------------------------------
    ! subroutine LocalDivergence_CFD
    ! estimate the divergence of the stress tensor by means of simple
    ! finite difference schemes. In the horizontal direction, numerical
    ! scheme is always centered finite difference. because of the
    ! surface and bottom boundary condition, scheme in the vertical
    ! direction changes from right-centered at the top, to center in the
    ! middle, to left-centered finite difference at the bottom.
    !-------------------------------------------------------------------
    SUBROUTINE localdivergence_cfd(f1,f2,f3,i3m,i3p,px3)
      REAL*8, INTENT(OUT) :: f1,f2,f3
      REAL*8, INTENT(IN) :: px3
      INTEGER, INTENT(IN) :: i3m, i3p

      INTEGER :: i1m,i1p,i2m,i2p

      ! neighbor samples
      i1m=mod(sx1+i1-2,sx1)+1
      i1p=mod(i1,sx1)+1
      i2m=mod(sx2+i2-2,sx2)+1
      i2p=mod(i2,sx2)+1

      f1= ( sig(i1p,i2,i3)%s11-sig(i1m,i2,i3)%s11 )/dx1/2._8 &
         +( sig(i1,i2p,i3)%s12-sig(i1,i2m,i3)%s12 )/dx2/2._8 &
         +( sig(i1,i2,i3p)%s13-sig(i1,i2,i3m)%s13 )/px3
      f2= ( sig(i1p,i2,i3)%s12-sig(i1m,i2,i3)%s12 )/dx1/2._8 &
         +( sig(i1,i2p,i3)%s22-sig(i1,i2m,i3)%s22 )/dx2/2._8 &
         +( sig(i1,i2,i3p)%s23-sig(i1,i2,i3m)%s23 )/px3
      f3= ( sig(i1p,i2,i3)%s13-sig(i1m,i2,i3)%s13 )/dx1/2._8 &
         +( sig(i1,i2p,i3)%s23-sig(i1,i2m,i3)%s23 )/dx2/2._8 &
         +( sig(i1,i2,i3p)%s33-sig(i1,i2,i3m)%s33 )/px3

    END SUBROUTINE localdivergence_cfd

  END SUBROUTINE equivalentbodyforce


  !---------------------------------------------------------------------
  !> function SourceSpectrum
  !! computes the equivalent body-forces for a buried dislocation,
  !! with strike-slip and dip-slip components,
  !! slip s, width W, length L in a rigidity mu
  !!
  !! \author sylvain barbot (06-25-07) - original form
  !---------------------------------------------------------------------
  SUBROUTINE sourcespectrum(mu,s,x,y,d, &
       L,W,strike,dip,rake,beta,dx1,dx2,dx3,f1,f2,f3)
    REAL*8, INTENT(IN) :: mu,s,x,y,d,L,W,strike,dip,rake,&
         beta,dx1,dx2,dx3
    REAL*4, DIMENSION(:,:,:), INTENT(INOUT) :: f1,f2,f3

    INTEGER :: i1,i2,i3,sx1,sx2,sx3
    REAL*8 :: k1,k2,k3,k1s,k2s,k3s,k1i,k3i, &
         cstrike,sstrike,cdip,sdip,cr,sr,k2r
    COMPLEX*8 :: cbuf1,cbuf2,cbuf3,source,image,&
         shift,scale,aperture,up,down
    COMPLEX*8, PARAMETER :: i=CMPLX(0._8,pi2)

    sx1=SIZE(f2,1)-2
    sx2=SIZE(f2,2)
    sx3=SIZE(f2,3)

    cstrike=cos(strike)
    sstrike=sin(strike)
    cdip=cos(dip)
    sdip=sin(dip)
    cr=cos(rake)
    sr=sin(rake)
    scale=CMPLX(i*mu*s*L*W)

    DO i3=1,sx3
       CALL wavenumber3(i3,sx3,dx3,k3)
       down=CMPLX(exp(-i*k3*(L/2._8+d)))
       up=conjg(down)
       DO i2=1,sx2
          CALL wavenumber2(i2,sx2,dx2,k2)
          DO i1=1,sx1/2+1
             CALL wavenumber1(i1,sx1,dx1,k1)

             !rotate the wavenumbers
             k2r= cstrike*k1-sstrike*k2
             k1s= cdip*k2r-sdip*k3
             k2s= sstrike*k1+cstrike*k2
             k3s= sdip*k2r+cdip*k3
             k1i= cdip*k2r+sdip*k3
             k3i=-sdip*k2r+cdip*k3
             
             !integrate at depth and along strike with raised cosine taper
             !and shift sources to x,y,z coordinate
             shift=CMPLX(exp(-i*(x*k1+y*k2)))
             aperture=CMPLX(scale*omegak(W*k2s,beta))
             source=CMPLX(omegak(L*k3s,beta)*aperture*shift*down)
             image =CMPLX(omegak(L*k3i,beta)*aperture*shift*up)

             !convolve source and image with a 1-D gaussian
             source=CMPLX(source*exp(-(pi*dx1*k1s)**2))
             image = CMPLX(image*exp(-(pi*dx1*k1i)**2))
             
             cbuf1= CMPLX(cdip*cstrike*( &
                  -(cr*k2s+sr*k3s)*source-(cr*k2s-sr*k3i)*image) &
                  +cr*sstrike*(-k1s*source-k1i*image) &
                  -sr*sdip*cstrike*(-k1s*source-k1i*image))
             !change -sr*sdip back to +sr*sdip above and below
             cbuf2=CMPLX(-cdip*sstrike*( &
                  -(cr*k2s+sr*k3s)*source-(cr*k2s-sr*k3i)*image) &
                  +cr*cstrike*(-k1s*source-k1i*image) &
                  -sr*sdip*sstrike*(-k1s*source-k1i*image))
             !change -sdip back to +sdip here
             cbuf3=CMPLX(-sdip*((-sr*k3s-cr*k2s)*source &
                  +(-sr*k3i+cr*k2s)*image) &
                  +sr*cdip*(-k1s*source+k1i*image))

             f1(2*i1-1:2*i1,i2,i3)=&
                  f1(2*i1-1:2*i1,i2,i3)+(/REAL(cbuf1),AIMAG(cbuf1)/)
             f2(2*i1-1:2*i1,i2,i3)=&
                  f2(2*i1-1:2*i1,i2,i3)+(/REAL(cbuf2),AIMAG(cbuf2)/)
             f3(2*i1-1:2*i1,i2,i3)=&
                  f3(2*i1-1:2*i1,i2,i3)+(/REAL(cbuf3),AIMAG(cbuf3)/)
          END DO
       END DO
    END DO

  END SUBROUTINE sourcespectrum


  !---------------------------------------------------------------------
  !> function SourceSpectrumHalfSpace
  !! computes the equivalent body-forces for a buried dislocation,
  !! with strike-slip and dip-slip components,
  !! slip s, width W, length L in a rigidity mu; sources are not imaged
  !!
  !! \author sylvain barbot (06-25-07) - original form
  !---------------------------------------------------------------------
  SUBROUTINE sourcespectrumhalfspace(mu,s,x,y,d, &
       L,W,strike,dip,rake,beta,dx1,dx2,dx3,f1,f2,f3)
    REAL*8, INTENT(IN) :: mu,s,x,y,d,L,W,strike,dip,rake,&
         beta,dx1,dx2,dx3
    REAL*4, DIMENSION(:,:,:), INTENT(INOUT) :: f1,f2,f3

    INTEGER :: i1,i2,i3,sx1,sx2,sx3
    REAL*8 :: k1,k2,k3,k1s,k2s,k3s, &
         cstrike,sstrike,cdip,sdip,cr,sr,k2r
    COMPLEX*8 :: cbuf1,cbuf2,cbuf3,source,&
         shift,scale,aperture,down
    COMPLEX*8, PARAMETER :: i=CMPLX(0._8,pi2)

    sx1=SIZE(f2,1)-2
    sx2=SIZE(f2,2)
    sx3=SIZE(f2,3)

    cstrike=cos(strike)
    sstrike=sin(strike)
    cdip=cos(dip)
    sdip=sin(dip)
    cr=cos(rake)
    sr=sin(rake)
    scale=CMPLX(i*mu*s*L*W)

    DO i3=1,sx3
       CALL wavenumber3(i3,sx3,dx3,k3)
       down=CMPLX(exp(-i*k3*(L/2._8+d)))
       DO i2=1,sx2
          CALL wavenumber2(i2,sx2,dx2,k2)
          DO i1=1,sx1/2+1
             CALL wavenumber1(i1,sx1,dx1,k1)
             !rotate the wavenumbers
             k2r= cstrike*k1-sstrike*k2
             k1s= cdip*k2r-sdip*k3
             k2s= sstrike*k1+cstrike*k2
             k3s= sdip*k2r+cdip*k3
             
             !convolve source and image with a 1-D gaussian
             !integrate at depth and along strike with raised cosine taper
             !and shift sources to x,y,z coordinate
             shift=CMPLX(exp(-i*(x*k1+y*k2)))
             aperture=CMPLX(scale*omegak(W*k2s,beta)*exp(-(pi*dx1*k1s)**2))
             source=CMPLX((omegak(L*k3s,beta)*aperture)*shift*down)

             cbuf1= CMPLX(cdip*cstrike*( &
                  -(cr*k2s+sr*k3s)*source) &
                  +cr*sstrike*(-k1s*source) &
                  -sr*sdip*cstrike*(-k1s*source))
             cbuf2=CMPLX(-cdip*sstrike*( &
                  -(cr*k2s+sr*k3s)*source) &
                  +cr*cstrike*(-k1s*source) &
                  -sr*sdip*sstrike*(-k1s*source))
             cbuf3=CMPLX(-sdip*((-sr*k3s-cr*k2s)*source) &
                  +sr*cdip*(-k1s*source))

             f1(2*i1-1:2*i1,i2,i3)=&
                  f1(2*i1-1:2*i1,i2,i3)+(/REAL(cbuf1),AIMAG(cbuf1)/)
             f2(2*i1-1:2*i1,i2,i3)=&
                  f2(2*i1-1:2*i1,i2,i3)+(/REAL(cbuf2),AIMAG(cbuf2)/)
             f3(2*i1-1:2*i1,i2,i3)=&
                  f3(2*i1-1:2*i1,i2,i3)+(/REAL(cbuf3),AIMAG(cbuf3)/)
          END DO
       END DO
    END DO

  END SUBROUTINE sourcespectrumhalfspace

  !---------------------------------------------------------------------
  !> function Source computes the equivalent body-forces
  !! in the space domain for a buried dislocation with strike-slip
  !! and dip-slip components, slip s, width W, length L in a rigidity mu
  !!
  !! Default (strike=0, dip=0, rake=0) is a vertical left-lateral
  !! strike-slip fault along the x2 axis. Default fault slip is
  !! represented with the double-couple equivalent body forces:
  !!
  !!\verbatim
  !!
  !!                   x1
  !!                   |
  !!                   |   ^  f2
  !!                   |   |<-----
  !!                   +---+------+---- x2
  !!                        ----->|
  !!                              v  f1
  !!
  !!\endverbatim
  !!
  !! \author sylvain barbot (06-29-07) - original form
  !---------------------------------------------------------------------
  SUBROUTINE source(mu,s,x,y,z,L,W,strike,dip,rake, &
       beta,sx1,sx2,sx3,dx1,dx2,dx3,f1,f2,f3,t1,t2,t3)
    INTEGER, INTENT(IN) :: sx1,sx2,sx3
    REAL*8, INTENT(IN) :: mu,s,x,y,z,L,W,strike,dip,rake, &
         beta,dx1,dx2,dx3
#ifdef ALIGN_DATA
    REAL*4, DIMENSION(sx1+2,sx2,sx3), INTENT(INOUT) :: f1,f2,f3
    REAL*4, DIMENSION(sx1+2,sx2), INTENT(INOUT) :: t1,t2,t3
#else
    REAL*4, DIMENSION(sx1,sx2,sx3), INTENT(INOUT) :: f1,f2,f3
    REAL*4, DIMENSION(sx1,sx2), INTENT(INOUT) :: t1,t2,t3
#endif

    INTEGER :: i1,i2,i3
    REAL*8 :: x1,x2,x3,x1s,x2s,x3s,x1i,x3i, &
         cstrike,sstrike,cdip,sdip,cr,sr,x2r, &
         sourc,image,scale,temp1,temp2,temp3, &
         dblcp,cplei,dipcs,dipci,xr,yr,zr,Wp,Lp
    REAL(8), DIMENSION(3) :: n,b
    TYPE(TENSOR) :: m

    cstrike=cos(strike)
    sstrike=sin(strike)
    cdip=cos(dip)
    sdip=sin(dip)
    cr=cos(rake)
    sr=sin(rake)
    scale=-mu*s

    ! effective tapered dimensions
    Wp=W*(1._8+2._8*beta)/2._8
    Lp=L*(1._8+2._8*beta)/2._8

    ! rotate centre coordinates of source and images
    x2r= cstrike*x  -sstrike*y
    xr = cdip   *x2r-sdip   *z
    yr = sstrike*x  +cstrike*y
    zr = sdip   *x2r+cdip   *z
    
    ! equivalent surface traction
    i3=1
    DO i2=1,sx2
       DO i1=1,sx1
          CALL shiftedcoordinates(i1,i2,i3,sx1,sx2,sx3, &
                                  dx1,dx2,dx3,x1,x2,x3)

          IF ((ABS(x1-x).GT.MAX(Lp,Wp)).OR.(ABS(x2-y).GT.MAX(Lp,Wp))) CYCLE

          x2r= cstrike*x1-sstrike*x2
          x1s= cdip*x2r-sdip*x3
          x1i= cdip*x2r+sdip*x3
          IF ((ABS(x1s-xr).GT.7.01*dx1).AND.(ABS(x1i-xr).GT.7.01*dx1)) CYCLE
          x2s= sstrike*x1+cstrike*x2
          x3s= sdip*x2r+cdip*x3
          x3i=-sdip*x2r+cdip*x3

          ! integrate at depth and along strike with raised cosine taper
          ! and shift sources to x,y,z coordinate
          temp1=gauss(x1s-xr,dx1)
          temp2=omega((x2s-yr)/W,beta)
          temp3=omega((x3s-zr)/L,beta)
          sourc=temp1*temp2*temp3

          ! add image
          temp1=gauss(x1i-xr,dx1)
          temp3=omega((x3i+zr)/L,beta)
          sourc=sourc+temp1*temp2*temp3

          ! surface normal vector components
          n(1)=+cdip*cstrike*sourc
          n(2)=-cdip*sstrike*sourc
          n(3)=-sdip*sourc

          ! burger vector (strike-slip)
          b(1)=sstrike*cr
          b(2)=cstrike*cr

          ! burger vector (dip-slip)
          b(1)=b(1)+cstrike*sdip*sr
          b(2)=b(2)-sstrike*sdip*sr
          b(3)=    +cdip*sr

          ! principal stress (symmetric deviatoric second-order tensor)
          m=n .sdyad. (mu*s*b)

          ! surface tractions
          t1(i1,i2)=t1(i1,i2)+m%s13
          t2(i1,i2)=t2(i1,i2)+m%s23
          t3(i1,i2)=t3(i1,i2)+m%s33
             
       END DO
    END DO

    ! equivalent body-force density
!$omp parallel do private(i1,i2,x1,x2,x3,x2r,x1s,x1i,x2s,x3s,x3i,temp1,temp2,temp3), &
!$omp private(sourc,dblcp,dipcs,image,cplei,dipci)
    DO i3=1,sx3/2
       CALL shiftedcoordinates(1,1,i3,sx1,sx2,sx3,dx1,dx2,dx3,x1,x2,x3)
       IF ((abs(x3-z).gt.Lp) .and. (abs(x3+z).gt.Lp)) CYCLE

       DO i2=1,sx2
          DO i1=1,sx1
             CALL shiftedcoordinates(i1,i2,i3,sx1,sx2,sx3,dx1,dx2,dx3,x1,x2,x3)
             IF ((ABS(x1-x) .GT. MAX(Wp,Lp)) .OR.  (abs(x2-y) .GT. MAX(Wp,Lp))) CYCLE

             x2r= cstrike*x1-sstrike*x2
             x1s= cdip*x2r-sdip*x3
             x1i= cdip*x2r+sdip*x3
             IF ((ABS(x1s-xr) .GT. 7.01_8*dx1) .AND. (ABS(x1i-xr) .GT. 7.01_8*dx1)) CYCLE
             x2s= sstrike*x1+cstrike*x2
             x3s= sdip*x2r+cdip*x3
             x3i=-sdip*x2r+cdip*x3
             
             !integrate at depth and along strike with raised cosine taper
             !and shift sources to x,y,z coordinate
             temp1=gauss(x1s-xr,dx1)
             temp2=omega((x2s-yr)/W,beta)
             temp3=omega((x3s-zr)/L,beta)
             sourc=scale  *gaussp(x1s-xr,dx1) &
                          *temp2 &
                          *temp3
             dblcp=scale/W*temp1 &
                          *omegap((x2s-yr)/W,beta) &
                          *temp3
             dipcs=scale/L*temp1 &
                          *temp2 &
                          *omegap((x3s-zr)/L,beta)

             temp1=gauss(x1i-xr,dx1)
             temp3=omega((x3i+zr)/L,beta)
             image=scale  *gaussp(x1i-xr,dx1) &
                          *temp2 &
                          *temp3
             cplei=scale/W*temp1 &
                          *omegap((x2s-yr)/W,beta) &
                          *temp3
             dipci=scale/L*temp1 &
                          *temp2 &
                          *omegap((x3i+zr)/L,beta)

             ! strike-slip component

             IF (2.01_8*DEG2RAD .GT. dip) THEN
                ! use method of images for subvertical faults
                f1(i1,i2,i3)=REAL(f1(i1,i2,i3) &
                     +cr*sstrike*(sourc+image) &
                       +cr*cdip*cstrike*(dblcp+cplei))
                f2(i1,i2,i3)=REAL(f2(i1,i2,i3) &
                     +cr*cstrike*(sourc+image) &
                     -cr*cdip*sstrike*(dblcp+cplei))
                f3(i1,i2,i3)=REAL(f3(i1,i2,i3) &
                     -cr*sdip*(dblcp-cplei))
             ELSE
                ! dipping faults do not use method of image
                f1(i1,i2,i3)=REAL(f1(i1,i2,i3) &
                     +cr*sstrike*(sourc) &
                     +cr*cdip*cstrike*(dblcp))
                f2(i1,i2,i3)=REAL(f2(i1,i2,i3) &
                     +cr*cstrike*(sourc) &
                     -cr*cdip*sstrike*(dblcp))
                 f3(i1,i2,i3)=REAL(f3(i1,i2,i3) &
                     -cr*sdip*(dblcp))
             END IF

             ! dip-slip component

             f1(i1,i2,i3)=REAL(f1(i1,i2,i3) &
                  +cdip*sr*cstrike*dipcs &
                  +sdip*sr*cstrike*sourc)
             f2(i1,i2,i3)=REAL(f2(i1,i2,i3) &
                  -cdip*sr*sstrike*dipcs &
                  -sdip*sr*sstrike*sourc)
             f3(i1,i2,i3)=REAL(f3(i1,i2,i3) &
                  +cdip*sr*sourc &
                  -sdip*sr*dipcs)

          END DO
       END DO
    END DO
!$omp end parallel do

  END SUBROUTINE source

  !---------------------------------------------------------------------
  !> function TensileSource
  !! computes the equivalent body-forces in the space domain for a buried
  !! tensile crack with opening s, width W, length L and Lame parameters
  !! lambda, mu.
  !!
  !! Default (strike=0, dip=0) is a vertical opening along the x2 axis.
  !! Default fault opening is represented with the double-couple
  !! equivalent body forces:
  !!
  !!\verbatim
  !!
  !!           x1           f1
  !!           |         ^^^^^^^
  !!           |         |||||||
  !!           | -f2 <--+-------+--> f2
  !!           |         |||||||
  !!           |         vvvvvvv
  !!           |           -f1
  !!           |
  !!           +----------------------------- x2
  !!
  !!\endverbatim
  !!
  !! The eigenstrain/potency tensor for a point source is
  !!
  !!\verbatim
  !!
  !!         | 1 0 0 |
  !!   E^i = | 0 0 0 |
  !!         | 0 0 0 |
  !!
  !!\endverbatim
  !!
  !! and the corresponding moment density for a point source is
  !!
  !!\verbatim
  !!
  !!                 | lambda+2*mu    0      0   |
  !!   m = C : E^i = |      0      lambda    0   |
  !!                 |      0         0   lambda |
  !!
  !!\endverbatim
  !!
  !! Moment density is integrated along the planar surface
  !!
  !!   \f[ box(x2) \delta(x1) box(x3) \f]
  !!
  !! where box(x) and delta(x) are the boxcar and the dirac delta
  !! functions, respectively. We use a tapered boxcar, omega_beta(x) and
  !! approximate the delta function by a small gaussian function.
  !! Finally, the equivalent body force is the divergence of the moment
  !! density tensor
  !!
  !!   \f[ f_i = - ( m_{ij} )_{,j} \f]
  !!
  !! derivatives are performed analytically on the gaussian and
  !! omega_beta functions.
  !!
  !! \author sylvain barbot (05-09-08) - original form
  !---------------------------------------------------------------------
  SUBROUTINE tensilesource(lambda,mu,s,x,y,z,L,W,strike,dip, &
       beta,sx1,sx2,sx3,dx1,dx2,dx3,f1,f2,f3)
    INTEGER, INTENT(IN) :: sx1,sx2,sx3
    REAL*8, INTENT(IN) :: lambda,mu,s,x,y,z,L,W,strike,dip,&
         beta,dx1,dx2,dx3
#ifdef ALIGN_DATA
    REAL*4, DIMENSION(sx1+2,sx2,sx3), INTENT(INOUT) :: f1,f2,f3
#else
    REAL*4, DIMENSION(sx1,sx2,sx3), INTENT(INOUT) :: f1,f2,f3
#endif

    INTEGER :: i1,i2,i3
    REAL*8 :: x1,x2,x3,x1s,x2s,x3s,x1i,x3i, &
         cstrike,sstrike,cdip,sdip,x2r,&
         sourc,image,scale1,scale2,temp1,temp2,temp3, &
         dblcp,cplei,dipcs,dipci,xr,yr,zr,Wp,Lp

    cstrike=cos(strike)
    sstrike=sin(strike)
    cdip=cos(dip)
    sdip=sin(dip)

    ! effective tapered dimensions
    Wp=W*(1._8+2._8*beta)/2._8
    Lp=L*(1._8+2._8*beta)/2._8

    ! rotate centre coordinates of source and images
    x2r= cstrike*x  -sstrike*y
    xr = cdip   *x2r-sdip   *z
    yr = sstrike*x  +cstrike*y
    zr = sdip   *x2r+cdip   *z
    scale1=-s*(lambda+2._8*mu)
    scale2=-s*lambda

    DO i3=1,sx3
       CALL shiftedcoordinates(i1,i2,i3,sx1,sx2,sx3,dx1,dx2,dx3,x1,x2,x3)
       IF ((abs(x3-z).gt.Lp) .and. (abs(x3+z).gt.Lp)) CYCLE

       DO i2=1,sx2
          DO i1=1,sx1
             CALL shiftedcoordinates(i1,i2,i3,sx1,sx2,sx3,dx1,dx2,dx3,x1,x2,x3)
             IF ((abs(x1-x).gt.Wp) .or.  (abs(x2-y).gt.Wp)) CYCLE

             x2r= cstrike*x1-sstrike*x2
             x1s= cdip*x2r-sdip*x3
             x1i= cdip*x2r+sdip*x3
             IF ((abs(x1s-xr).gt.7.01*dx1).and.(abs(x1i-xr).gt.7.01*dx1)) CYCLE
             x2s= sstrike*x1+cstrike*x2
             x3s= sdip*x2r+cdip*x3
             x3i=-sdip*x2r+cdip*x3

             !integrate at depth and along strike with raised cosine taper
             !and shift sources to x,y,z coordinate
             temp1=gauss(x1s-xr,dx1)
             temp2=omega((x2s-yr)/W,beta)
             temp3=omega((x3s-zr)/L,beta)
             sourc=scale1  *gaussp(x1s-xr,dx1) &
                           *temp2 &
                           *temp3
             dblcp=scale2/W*temp1 &
                           *omegap((x2s-yr)/W,beta) &
                           *temp3
             dipcs=scale2/L*temp1 &
                           *temp2 &
                           *omegap((x3s-zr)/L,beta)

             temp1=gauss(x1i-xr,dx1)
             temp3=omega((x3i+zr)/L,beta)
             image=scale1  *gaussp(x1i-xr,dx1) &
                           *temp2 &
                           *temp3
             cplei=scale2/W*temp1 &
                           *omegap((x2s-yr)/W,beta) &
                           *temp3
             dipci=scale2/L*temp1 &
                           *temp2 &
                           *omegap((x3i+zr)/L,beta)

             ! force moments in original coordinate system

             f1(i1,i2,i3)=REAL(f1(i1,i2,i3) &
                  +cstrike*cdip*(sourc+image) &
                  +sstrike*(dblcp+cplei) &
                  +cstrike*sdip*(dipcs+dipci))
             f2(i1,i2,i3)=REAL(f2(i1,i2,i3) &
                  -sstrike*cdip*(sourc+image) &
                  +cstrike*(dblcp+cplei) &
                  -sstrike*sdip*(dipcs+dipci))
             f3(i1,i2,i3)=REAL(f3(i1,i2,i3) &
                  -sdip*(sourc-image) &
                  +cdip*(dipcs-dipci))

          END DO
       END DO
    END DO

  END SUBROUTINE tensilesource

  !---------------------------------------------------------------------
  !! function MogiSource 
  !! computes the equivalent body-forces in the space domain for a buried 
  !! dilatation point source.
  !!
  !! The point-source opening o with at position xs in the half space is
  !! associated with eigenstrain
  !!
  !!      \f[ E^i = o \frac{1}{3} I \delta(x-x_s) \f]
  !!
  !! where I is the diagonal tensor and delta is the Dirac delta function
  !! (or in index notation E^i_{ij} = o delta_{ij} / 3 delta(xs) ) and 
  !! with the moment density
  !!
  !!      \f[ m = C : E^i = K o I \delta(x-x_s) \f]
  !!
  !! The equivalent body-force density is
  !!
  !!      \f[ f = - \nabla \cdot m = K o \nabla \delta(x-x_s) \f]
  !!
  !! where nabla is the gradient operator. Default source opening is 
  !! represented with the isotropic equivalent body-force density:
  !!
  !!\verbatim
  !!
  !!                   x1
  !!                   |      f1
  !!                   |      ^
  !!                   |  f2  |  f2
  !!                   +---<--+-->---- x2
  !!                          |
  !!                          v  f1
  !!
  !!                   x3
  !!                   |      f3
  !!                   |      ^
  !!                   |  f2  |  f2
  !!                   +---<--+-->---- x2
  !!                          |
  !!                          v  f3
  !!
  !!\endverbatim
  !!
  !! \author sylvain barbot (03-24-09) - original form
  !---------------------------------------------------------------------
  SUBROUTINE mogisource(lambda,mu,o,xs,ys,zs,sx1,sx2,sx3,dx1,dx2,dx3,f1,f2,f3)
    INTEGER, INTENT(IN) :: sx1,sx2,sx3
    REAL*8, INTENT(IN) :: lambda,mu,o,xs,ys,zs,dx1,dx2,dx3
#ifdef ALIGN_DATA
    REAL*4, DIMENSION(sx1+2,sx2,sx3), INTENT(INOUT) :: f1,f2,f3
#else
    REAL*4, DIMENSION(sx1,sx2,sx3), INTENT(INOUT) :: f1,f2,f3
#endif

    INTEGER :: i1,i2,i3
    REAL*8 :: x1,x2,x3,source1,source2,source3, &
         image1,image2,image3,scale,temp1,temp2,temp3,Wp,Lp

    scale=-(lambda+2._8*mu/3._8)*o ! -kappa*o

    ! effective dimensions
    Wp=6._8*MAX(dx1,dx2,dx3)
    Lp=6._8*MAX(dx1,dx2,dx3)

    DO i3=1,sx3
       CALL shiftedcoordinates(i1,i2,i3,sx1,sx2,sx3,dx1,dx2,dx3,x1,x2,x3)
       IF ((abs(x3-zs).gt.Lp) .and. (abs(x3+zs).gt.Lp)) CYCLE
       
       DO i2=1,sx2
          DO i1=1,sx1
             CALL shiftedcoordinates(i1,i2,i3,sx1,sx2,sx3,dx1,dx2,dx3,x1,x2,x3)
             IF ((abs(x1-xs).gt.Wp) .or.  (abs(x2-ys).gt.Wp)) CYCLE

             temp1=gauss(x1-xs,dx1)
             temp2=gauss(x2-ys,dx2)
             temp3=gauss(x3-zs,dx3)

             source1=scale*gaussp(x1-xs,dx1)*temp2*temp3
             source2=scale*temp1*gaussp(x2-ys,dx2)*temp3
             source3=scale*temp1*temp2*gaussp(x3-zs,dx3)

             temp3=gauss(x3+zs,dx3)

             image1=scale*gaussp(x1-xs,dx1)*temp2*temp3
             image2=scale*temp1*gaussp(x2-ys,dx2)*temp3
             image3=scale*temp1*temp2*gaussp(x3+zs,dx3)

             ! equivalent body-force density
             f1(i1,i2,i3)=REAL(f1(i1,i2,i3)+(source1+image1))
             f2(i1,i2,i3)=REAL(f2(i1,i2,i3)+(source2+image2))
             f3(i1,i2,i3)=REAL(f3(i1,i2,i3)+(source3-image3))

          END DO
       END DO
    END DO

  END SUBROUTINE mogisource

  !---------------------------------------------------------------------
  !> subroutine Traction 
  !! assigns the traction vector at the surface.
  !!
  !! \author sylvain barbot (07-19-07) - original form
  !---------------------------------------------------------------------
  SUBROUTINE traction(mu,e,sx1,sx2,dx1,dx2,t,Dt,t3,rate)
    TYPE(EVENT_STRUC), INTENT(IN) :: e
    INTEGER, INTENT(IN) :: sx1,sx2
    REAL*8, INTENT(IN) :: mu,dx1,dx2,t,Dt
#ifdef ALIGN_DATA
    REAL*4, DIMENSION(sx1+2,sx2), INTENT(INOUT) :: t3
#else
    REAL*4, DIMENSION(sx1,sx2), INTENT(INOUT) :: t3
#endif
    LOGICAL, INTENT(IN), OPTIONAL :: rate

    INTEGER :: i,i1,i2,i3
    LOGICAL :: israte
    REAL*8 :: period,phi,amp,L,W,Lp,Wp,x1,x2,x3,x,y,beta

    REAL*8, PARAMETER :: pi=3.141592653589793115997963468544185161_8

    IF (PRESENT(rate)) THEN
       israte=rate
    ELSE
       israte=.FALSE.
    END IF

    ! loop over traction sources
    DO i=1,e%nl

       x=e%l(i)%x
       y=e%l(i)%y

       L=e%l(i)%length
       W=e%l(i)%width

       beta=e%l(i)%beta

       ! effective tapered dimensions
       Lp=L*(1._8+2._8*beta)/2._8
       Wp=W*(1._8+2._8*beta)/2._8

       i3=1
       DO i2=1,sx2
          DO i1=1,sx1
             CALL shiftedcoordinates(i1,i2,i3,sx1,sx2,1, &
                                     dx1,dx2,1.d8,x1,x2,x3)

             IF ((ABS(x1-x).GT.MAX(Lp,Wp)).OR.(ABS(x2-y).GT.MAX(Lp,Wp))) CYCLE

             amp=omega((x1-x)/L,beta)* &
                 omega((x2-y)/W,beta)* &
                 mu*e%l(i)%slip

             IF (israte) THEN
                IF (0 .NE. period) THEN
                   ! surface tractions rate
                   period=e%l(i)%period
                   phi=e%l(i)%phase

                   t3(i1,i2)=REAL(t3(i1,i2)-amp*(sin(2*pi*(t+Dt)/period+phi)-sin(2*pi*t/period+phi)))
                END IF
             ELSE
                IF (e%l(i)%period .LE. 0) THEN
                   ! surface tractions
                   t3(i1,i2)=REAL(t3(i1,i2)-amp)
                END IF
             END IF
          END DO
       END DO
    END DO

  END SUBROUTINE traction

  !---------------------------------------------------------------------
  !! function MomentDensityShear
  !! computes the inelastic irreversible moment density in the space
  !! domain corresponding to a buried dislocation with strike-slip and
  !! dip-slip components (pure shear). A fault along a surface of normal
  !! n_i with a burger vector s_i, is associated with the eigenstrain
  !!
  !!   E^i_ij = 1/2 ( n_i s_j + s_i n_j )
  !!
  !! In a heterogeneous medium of elastic moduli tensor C_ijkl, the
  !! corresponding moment density tensor is
  !!
  !!   m_ij = C_ijkl E^i_kl
  !!
  !! where C = C(x) is a function of space. Equivalent body forces
  !! representing the set of dislocations can be obtained by evaluating
  !! the divergence of the moment density tensor
  !!
  !!   f_i = - ( m_ji ),j
  !!
  !! using the function "EquivalentBodyForce" in this module.
  !!
  !! The default dislocation extends in the x2 direction, with a normal
  !! in the x1 direction. Using the following angular convention,
  !!
  !!\verbatim
  !!
  !!           x1            !           x1
  !!   n  theta |            !   n   phi  |
  !!     \  ____|            !     \  ____|
  !!       \    |            !       \    |
  !!         \  |            !         \  |
  !!      -----\+------ x2   !      -----\+------ x3
  !!        (x3 down)        !         (x2 up)
  !!
  !!\endverbatim
  !!
  !! where theta is the strike and phi is the dip (internal convention),
  !! and introducting the rotation matrices
  !!
  !!\verbatim
  !!
  !!        |  cos(theta)   sin(theta)    0 |
  !!   R1 = | -sin(theta)   cos(theta)    0 |
  !!        |      0             0        1 |
  !!
  !!        |  cos(phi)     0     sin(phi)  |
  !!   R2 = |     0         1        0      |
  !!        | -sin(phi)     0     cos(phi)  |
  !!
  !!\endverbatim
  !!
  !! a normal vector n of arbitrary orientation and the corresponding
  !! strike-slip and dip-slip vector, s and d respectively, are
  !!
  !!\verbatim
  !!
  !!             | 1 |             | 0 |             | 0 |
  !!   n = R1 R2 | 0 |,  s = R1 R2 | 1 |,  d = R1 R2 | 0 |
  !!             | 0 |             | 0 |             | 1 |
  !!
  !!\endverbatim
  !!
  !! vector n, s and d are orthogonal and the corresponding moment
  !! density second order tensor is deviatoric. The method of images is
  !! used to avoid tapering of the fault at the surface.
  !!
  !! \author sylvain barbot (03-02-08) - original form
  !---------------------------------------------------------------------
  SUBROUTINE momentdensityshear(mu,slip,x,y,z,L,W,strike,dip,rake, &
       beta,sx1,sx2,sx3,dx1,dx2,dx3,sig)
    INTEGER, INTENT(IN) :: sx1,sx2,sx3
    REAL*8, INTENT(IN) :: mu,slip,x,y,z,L,W,strike,dip,rake,&
         beta,dx1,dx2,dx3
    TYPE(TENSOR), INTENT(INOUT), DIMENSION(sx1,sx2,sx3) :: sig

    INTEGER :: i1,i2,i3
    REAL*4 :: rmu
    REAL*8 :: x1,x2,x3,x1s,x2s,x3s,x1i,x3i, &
         cstrike,sstrike,cdip,sdip,cr,sr,x2r,&
         aperture,temp1,temp2,temp3,xr,yr,zr,Wp,Lp,dum
    REAL*8, DIMENSION(3) :: n,s
    TYPE(TENSOR) :: Ei

    rmu=2._4*REAL(mu,4)

    cstrike=cos(strike)
    sstrike=sin(strike)
    cdip=cos(dip)
    sdip=sin(dip)
    cr=cos(rake)
    sr=sin(rake)

    ! effective tapered dimensions
    Wp=W*(1._8+2._8*beta)/2._8
    Lp=L*(1._8+2._8*beta)/2._8

    ! rotate centre coordinates of source and images
    x2r= cstrike*x  -sstrike*y
    xr = cdip   *x2r-sdip   *z
    yr = sstrike*x  +cstrike*y
    zr = sdip   *x2r+cdip   *z

    DO i3=1,sx3
       x3=DBLE(i3-1)*dx3
       IF (abs(x3-z) .gt. Lp) CYCLE

       DO i2=1,sx2
          DO i1=1,sx1
             CALL shiftedcoordinates(i1,i2,i3,sx1,sx2,sx3, &
                  dx1,dx2,dx3,x1,x2,dum)

             IF ((abs(x1-x).gt.Wp) .or.  (abs(x2-y).gt.Wp)) CYCLE

             x2r= cstrike*x1-sstrike*x2
             x1s= cdip*x2r-sdip*x3
             x1i= cdip*x2r+sdip*x3
             IF ((abs(x1s-xr).gt.7.01*dx1).and.(abs(x1i-xr).gt.7.01*dx1)) CYCLE
             x2s= sstrike*x1+cstrike*x2
             x3s= sdip*x2r+cdip*x3
             x3i=-sdip*x2r+cdip*x3

             ! integrate at depth and along strike with raised cosine taper
             ! and shift sources to x,y,z coordinate
             temp1=gauss(x1s-xr,dx1)
             temp2=omega((x2s-yr)/W,beta)
             temp3=omega((x3s-zr)/L,beta)
             aperture=temp1*temp2*temp3

             ! add image
             temp1=gauss(x1i-xr,dx1)
             temp3=omega((x3i+zr)/L,beta)
             aperture=aperture+temp1*temp2*temp3

             ! surface normal vector components
             n(1)=+cdip*cstrike*aperture
             n(2)=-cdip*sstrike*aperture
             n(3)=-sdip*aperture

             ! strike-slip component
             s(1)=sstrike*cr
             s(2)=cstrike*cr

             ! dip-slip component
             s(1)=s(1)+cstrike*sdip*sr
             s(2)=s(2)-sstrike*sdip*sr
             s(3)=    +cdip*sr

             ! eigenstrain (symmetric deviatoric second-order tensor)
             Ei=n .sdyad. (slip*s)

             ! moment density (pure shear)
             sig(i1,i2,i3)=sig(i1,i2,i3) .plus. (rmu .times. Ei)
             
          END DO
       END DO
    END DO

  END SUBROUTINE momentdensityshear

  !---------------------------------------------------------------------
  !> function MomentDensityTensile
  !! computes the inelastic irreversible moment density in the space
  !! domain corresponding to a buried dislocation with opening (open
  !! crack). A fault along a surface of normal n_i with a burger vector
  !! s_i, is associated with the eigenstrain
  !!
  !!   \f[ E^i_{ij} = \frac{1}{2} ( n_i s_j + s_i n_j ) \f]
  !!
  !! The eigenstrain/potency tensor for a point source opening crack is
  !!
  !!\verbatim
  !!
  !!         | 1 0 0 |
  !!   E^i = | 0 0 0 |
  !!         | 0 0 0 |
  !!
  !!\endverbatim
  !!
  !! In a heterogeneous medium of elastic moduli tensor C_ijkl, the
  !! corresponding moment density tensor is
  !!
  !!   \f[ m_{ij} = C_{ijkl} E^i_{kl} = \lambda E^i_{kk} \delta_{ij} + 2 \mu E^i_{ij} \f]
  !!
  !! where C = C(x) is a function of space. (We use isotropic elastic
  !! solid, and heterogeneous elastic moduli tensor simplifies to
  !! mu=mu(x) and lambda = lambda(x).) The moment density for a point
  !! source opening crack is
  !!
  !!\verbatim
  !!
  !!          | lambda+2*mu    0      0   |
  !!   m(x) = |      0      lambda    0   |
  !!          |      0         0   lambda |
  !!
  !!\endverbatim
  !!
  !! Moment density m(x) is integrated along the planar surface
  !!
  !!   box(x2) delta (x1) box(x3)
  !!
  !! where box(x) and delta(x) are the boxcar and the dirac delta
  !! functions, respectively. Equivalent body forces representing the
  !! set of dislocations can be obtained by evaluating the divergence
  !! of the moment density tensor
  !!
  !!   \f[ f_i = - ( m_{ji} ),j \f]
  !!
  !! The corresponding equivalent surface traction is simply
  !!
  !!   \f[ t_i = m_{ij} n_j \f]
  !!
  !! Both equivalent body forces and equivalent surface traction are
  !! computed using the function "EquivalentBodyForce" in this module.
  !!
  !! The default dislocation extends in the x2 direction, with a normal
  !! in the x1 direction. Using the following angular convention,
  !!
  !!\verbatim
  !!
  !!           x1            !           x1
  !!   n  theta |            !   n   phi  |
  !!     \  ____|            !     \  ____|
  !!       \    |            !       \    |
  !!         \  |            !         \  |
  !!      -----\+------ x2   !      -----\+------ x3
  !!        (x3 down)        !         (x2 up)
  !!
  !!\endverbatim
  !!
  !! where theta is the strike and phi is the dip, in internal
  !! convention. (Internal angular convention does not correspond to
  !! usual angular convention of geology and conversion between the two
  !! standard is necessary.) Introducting the rotation matrices,
  !!
  !!\verbatim
  !!
  !!        |  cos(theta)   sin(theta)    0 |
  !!   R1 = | -sin(theta)   cos(theta)    0 |
  !!        |      0             0        1 |
  !!
  !!        |  cos(phi)     0     sin(phi)  |
  !!   R2 = |     0         1        0      |
  !!        | -sin(phi)     0     cos(phi)  |
  !!
  !!\endverbatim
  !!
  !! a normal vector n of arbitrary orientation and the corresponding
  !! slip vector s are
  !!
  !!\verbatim
  !!
  !!             | 1 |                 | 1 |
  !!   n = R1 R2 | 0 |,  s = n = R1 R2 | 0 |
  !!             | 0 |                 | 0 |
  !!
  !!\endverbatim
  !!
  !! The method of images is used to avoid tapering of the fault at
  !! the surface.
  !!
  !! \author sylvain barbot (03-02-08) - original form
  !---------------------------------------------------------------------
  SUBROUTINE momentdensitytensile(lambda,mu,slip,x,y,z,L,W,strike,dip,rake, &
       beta,sx1,sx2,sx3,dx1,dx2,dx3,sig)
    INTEGER, INTENT(IN) :: sx1,sx2,sx3
    REAL*8, INTENT(IN) :: lambda,mu,slip,x,y,z,L,W,strike,dip,rake,&
         beta,dx1,dx2,dx3
    TYPE(TENSOR), INTENT(INOUT), DIMENSION(sx1,sx2,sx3) :: sig

    INTEGER :: i1,i2,i3
    REAL*8 :: x1,x2,x3,x1s,x2s,x3s,x1i,x3i, &
         cstrike,sstrike,cdip,sdip,cr,sr,x2r,&
         aperture,temp1,temp2,temp3,xr,yr,zr,Wp,Lp,dum
    REAL*8, DIMENSION(3) :: n
    TYPE(TENSOR) :: Ei

    cstrike=cos(strike)
    sstrike=sin(strike)
    cdip=cos(dip)
    sdip=sin(dip)
    cr=cos(rake)
    sr=sin(rake)

    ! effective tapered dimensions
    Wp=W*(1._8+2._8*beta)/2._8
    Lp=L*(1._8+2._8*beta)/2._8

    ! rotate centre coordinates of source and images
    x2r= cstrike*x  -sstrike*y
    xr = cdip   *x2r-sdip   *z
    yr = sstrike*x  +cstrike*y
    zr = sdip   *x2r+cdip   *z

    DO i3=1,sx3
       x3=DBLE(i3-1)*dx3
       IF (abs(x3-z) .gt. Lp) CYCLE

       DO i2=1,sx2
          DO i1=1,sx1
             CALL shiftedcoordinates(i1,i2,i3,sx1,sx2,sx3, &
                  dx1,dx2,dx3,x1,x2,dum)

             IF ((abs(x1-x).gt.Wp) .or.  (abs(x2-y).gt.Wp)) CYCLE

             x2r= cstrike*x1-sstrike*x2
             x1s= cdip*x2r-sdip*x3
             x1i= cdip*x2r+sdip*x3
             IF ((abs(x1s-xr).gt.7.01*dx1).and.(abs(x1i-xr).gt.7.01*dx1)) CYCLE
             x2s= sstrike*x1+cstrike*x2
             x3s= sdip*x2r+cdip*x3
             x3i=-sdip*x2r+cdip*x3

             ! integrate at depth and along strike with raised cosine taper
             ! and shift sources to x,y,z coordinate
             temp1=gauss(x1s-xr,dx1)
             temp2=omega((x2s-yr)/W,beta)
             temp3=omega((x3s-zr)/L,beta)
             aperture=temp1*temp2*temp3

             ! add image
             temp1=gauss(x1i-xr,dx1)
             temp3=omega((x3i+zr)/L,beta)
             aperture=aperture+temp1*temp2*temp3

             ! surface normal vector components
             n(1)=+cdip*cstrike*aperture
             n(2)=-cdip*sstrike*aperture
             n(3)=-sdip*aperture

             ! eigenstrain (symmetric second-order tensor)
             Ei=n .sdyad. (slip*n)

             ! moment density (isotropic Hooke's law)
             CALL isotropicstressstrain(Ei,lambda,mu)
             sig(i1,i2,i3)=sig(i1,i2,i3) .plus. Ei
             
          END DO
       END DO
    END DO

  END SUBROUTINE momentdensitytensile

  !---------------------------------------------------------------------
  !! function MomentDensityMogi
  !! computes the inelastic irreversible moment density in the space
  !! domain corresponding to a buried Mogi source. 
  !! The Mogi source is associated with the eigenstrain
  !!
  !!   \f[ E^i_{ij} = o \frac{1}{3} \delta_{ij} \f]
  !!
  !! In a heterogeneous medium of elastic moduli tensor C_ijkl, the
  !! corresponding moment density tensor is
  !!
  !!   \f[ m_{ij} = C_{ijkl} E^i_{kl} \f]
  !!
  !! where C = C(x) is a function of space. Equivalent body forces
  !! representing the set of dislocations can be obtained by evaluating
  !! the divergence of the moment density tensor
  !!
  !!   \f[ f_i = - ( m_{ji} ),j \f]
  !!
  !! using the function "EquivalentBodyForce" in this module.
  !!
  !! \author sylvain barbot (03-24-09) - original form
  !---------------------------------------------------------------------
  SUBROUTINE momentdensitymogi(lambda,mu,o,xs,ys,zs,sx1,sx2,sx3,dx1,dx2,dx3,sig)
    INTEGER, INTENT(IN) :: sx1,sx2,sx3
    REAL*8, INTENT(IN) :: lambda,mu,o,xs,ys,zs,dx1,dx2,dx3
    TYPE(TENSOR), INTENT(INOUT), DIMENSION(sx1,sx2,sx3) :: sig

    INTEGER :: i1,i2,i3
    REAL*8 :: x1,x2,x3,Wp,Lp,dum,kappa,gamma,gammai
    TYPE(TENSOR) :: m

    kappa=lambda+2._8/3._8*mu

    ! effective tapered dimensions
    Wp=6._8*MAX(dx1,dx2,dx3)
    Lp=6._8*MAX(dx1,dx2,dx3)

    DO i3=1,sx3
       x3=DBLE(i3-1)*dx3
       IF (abs(x3-zs) .gt. Lp) CYCLE

       DO i2=1,sx2
          DO i1=1,sx1
             CALL shiftedcoordinates(i1,i2,i3,sx1,sx2,sx3, &
                  dx1,dx2,dx3,x1,x2,dum)

             IF ((abs(x1-xs).gt.Wp) .or.  (abs(x2-ys).gt.Wp)) CYCLE

             ! amplitude of eigenstrain
             gamma =o*gauss(x1-xs,dx1)*gauss(x2-ys,dx2)*gauss(x3-zs,dx3)

             ! add image
             gammai=o*gauss(x1-xs,dx1)*gauss(x2-ys,dx2)*gauss(x3+zs,dx3)

             ! amplitude of moment density
             gamma=kappa*gamma
             gammai=kappa*gammai

             ! eigenstrain (diagonal second-order tensor)
             m=TENSOR(REAL(gamma),0,0,REAL(gamma),0,REAL(gamma))

             ! moment density (pure shear)
             sig(i1,i2,i3)=sig(i1,i2,i3) .plus. m
             
          END DO
       END DO
    END DO

  END SUBROUTINE momentdensitymogi

  !---------------------------------------------------------------------
  !> function Plane
  !! computes the three components, n1, n2 and n3, of the normal vector
  !! corresponding to a rectangular surface of finite size. The plane
  !! is defined by its orientation (strike and dip) and dimension.
  !!
  !!\verbatim
  !!
  !!              W
  !!       +-------------+
  !!       |             |
  !!     L |      +      | - - - > along strike direction
  !!       |   (x,y,z)   |
  !!       +-------------|
  !!              |
  !!              v
  !!      down-dip direction
  !!
  !!\endverbatim
  !!
  !! in the default orientation, for which strike=0 and dip=0, the plane
  !! is vertical along the x2 axis, such as n2(x) = n3(x) = 0 for all x.
  !! internal angular conventions are as follows:
  !!
  !!\verbatim
  !!
  !!             n   x1                          n   x1
  !!              \   |                           \   |
  !!               \  |                            \  |
  !!   90 - strike  \ |                  90 - dip   \ |
  !!               ( \|                            ( \|
  !!        ----------+------ x2            ----------+------ x3
  !!              (x3 down)                       (x2 up)
  !!
  !!\endverbatim
  !!
  !! edges of the rectangle are tapered.
  !!
  !! \author sylvain barbot (09-15-07) - original form
  !---------------------------------------------------------------------
  SUBROUTINE plane(x,y,z,L,W,strike,dip, &
       beta,sx1,sx2,sx3,dx1,dx2,dx3,n1,n2,n3)
    INTEGER, INTENT(IN) :: sx1,sx2,sx3
    REAL*8, INTENT(IN) :: x,y,z,L,W,strike,dip,beta,dx1,dx2,dx3
#ifdef ALIGN_DATA
    REAL*4, DIMENSION(sx1+2,sx2,sx3), INTENT(INOUT) :: n1,n2,n3
#else
    REAL*4, DIMENSION(sx1,sx2,sx3), INTENT(INOUT) :: n1,n2,n3
#endif

    INTEGER :: i1,i2,i3
    REAL*8 :: x1,x2,x3,x1s,x2s,x3s,x1i,x3i, &
         cstrike,sstrike,cdip,sdip,x2r,&
         temp1,temp2,temp3,sourc,image,xr,yr,zr,Wp,Lp,dum

    cstrike=cos(strike)
    sstrike=sin(strike)
    cdip=cos(dip)
    sdip=sin(dip)

    ! effective tapered dimensions
    Wp=W*(1._8+2._8*beta)/2._8
    Lp=L*(1._8+2._8*beta)/2._8

    ! rotate centre coordinates of source and images
    x2r= cstrike*x  -sstrike*y
    xr = cdip   *x2r-sdip   *z
    yr = sstrike*x  +cstrike*y
    zr = sdip   *x2r+cdip   *z

    DO i3=1,sx3
       x3=DBLE(i3-1)*dx3
       IF ((abs(x3-z).gt.Lp) .and. (abs(x3+z).gt.Lp)) CYCLE

       DO i2=1,sx2
          DO i1=1,sx1
             CALL shiftedcoordinates(i1,i2,i3,sx1,sx2,sx3, &
                  dx1,dx2,dx3,x1,x2,dum)
             IF ((abs(x1-x).gt.Wp) .or.  (abs(x2-y).gt.Wp)) CYCLE

             x2r= cstrike*x1-sstrike*x2
             x1s= cdip*x2r-sdip*x3
             x1i= cdip*x2r+sdip*x3
             IF ((abs(x1s-xr).gt.7.01*dx1).and.(abs(x1i-xr).gt.7.01*dx1)) CYCLE
             x2s= sstrike*x1+cstrike*x2
             x3s= sdip*x2r+cdip*x3
             x3i=-sdip*x2r+cdip*x3

             !integrate at depth and along strike with raised cosine taper
             !and shift sources to x,y,z coordinate
             temp1=gauss(x1s-xr,dx1)
             temp2=omega((x2s-yr)/W,beta)
             temp3=omega((x3s-zr)/L,beta)
             sourc=temp1*temp2*temp3

             temp1=gauss(x1i-xr,dx1)
             temp3=omega((x3i+zr)/L,beta)
             image=temp1*temp2*temp3

             ! surface normal vector components
             n1(i1,i2,i3)=REAL(n1(i1,i2,i3)+cdip*cstrike*(sourc+image))
             n2(i1,i2,i3)=REAL(n2(i1,i2,i3)-cdip*sstrike*(sourc+image))
             n3(i1,i2,i3)=REAL(n3(i1,i2,i3)-sdip*(sourc+image))
             
          END DO
       END DO
    END DO

  END SUBROUTINE plane

  !---------------------------------------------------------------------
  !> function MonitorStressField
  !! samples a stress field along a specified planar surface.
  !!
  !! \author sylvain barbot (10-16-07) - original form
  !---------------------------------------------------------------------
  SUBROUTINE monitorstressfield(x,y,z,L,W,strike,dip,beta, &
       sx1,sx2,sx3,dx1,dx2,dx3,sig,patch)
    INTEGER, INTENT(IN) :: sx1,sx2,sx3
    REAL*8, INTENT(IN) :: x,y,z,L,W,strike,dip,beta,dx1,dx2,dx3
    TYPE(TENSOR), INTENT(IN), DIMENSION(sx1,sx2,sx3) :: sig
    TYPE(SLIPPATCH_STRUCT), ALLOCATABLE, DIMENSION(:,:), INTENT(OUT) :: patch

    INTEGER :: px2,px3,j2,j3,status
    REAL*8 :: x1,x2,x3,xr,yr,zr,Wp,Lp, &
         cstrike,sstrike,cdip,sdip
    TYPE(TENSOR) :: lsig

    cstrike=cos(strike)
    sstrike=sin(strike)
    cdip=cos(dip)
    sdip=sin(dip)

    ! effective tapered dimensions
    Wp=W*(1._8+2._8*beta) ! horizontal dimension for vertical fault
    Lp=L*(1._8+2._8*beta) ! depth for a vertical fault

    px3=fix(Lp/dx3)
    px2=fix(Wp/dx2)

    ALLOCATE(patch(px2+1,px3+1),STAT=status)
    IF (status>0) STOP "could not allocate the slip patches for export"

    DO j3=1,px3+1
       DO j2=1,px2+1

          CALL ref2local(x,y,z,xr,yr,zr)
          
          ! no translation in out of plane direction
          yr=REAL(yr)+REAL((DBLE(j2)-DBLE(px2)/2._8-1._8)*dx2)
          zr=REAL(zr)+REAL((DBLE(j3)-DBLE(px3)/2._8-1._8)*dx3)
          
          CALL local2ref(xr,yr,zr,x1,x2,x3)
          
          ! discard out-of-bound locations
          IF (  (x1 .gt. DBLE(sx1/2-1)*dx1) .or. (x1 .lt. -DBLE(sx1/2)*dx1) &
           .or. (x2 .gt. DBLE(sx2/2-1)*dx2) .or. (x2 .lt. -DBLE(sx2/2)*dx2) &
           .or. (x3 .gt. DBLE(sx3-1)*dx3) .or. (x3 .lt. 0._8)  ) THEN
             lsig=TENSOR(0._8,0._8,0._8,0._8,0._8,0._8)
          ELSE
             CALL sampletensor(x1,x2,x3,dx1,dx2,dx3,sx1,sx2,sx3,sig,lsig)
          END IF

          patch(j2,j3)=SLIPPATCH_STRUCT(x1,x2,x3,yr,zr, &
                                        0._8,0._8,0._8,0._8,0._8,0._8,0._8,lsig)

       END DO
    END DO

  CONTAINS

    !--------------------------------------------------------------
    !> subroutine sample
    !! interpolates the value of a discretized 3-dimensional field
    !! at a subpixel location. method consists in correlating the
    !! 3D field with a delta function filter. the delta function is
    !! approximated with a narrow normalized gaussian.
    !!
    !! \author sylvain barbot (10-17-07) - original form
    !--------------------------------------------------------------
    SUBROUTINE sampletensor(x1,x2,x3,dx1,dx2,dx3,sx1,sx2,sx3,sig,lsig)
      INTEGER, INTENT(IN) :: sx1,sx2,sx3
      REAL*8, INTENT(IN) :: x1,x2,x3,dx1,dx2,dx3
      TYPE(TENSOR), INTENT(IN), DIMENSION(sx1,sx2,sx3) :: sig
      TYPE(TENSOR), INTENT(OUT) :: lsig
    
      INTEGER :: i,j,k,l1,l2,l3,i1p,i2p,i3p
      INTEGER, PARAMETER :: RANGE=2
      REAL*8 :: sum,weight,x,y,z
      REAL*8, PARAMETER :: EPS=1e-2

      sum=0._8
      lsig=TENSOR(0._8,0._8,0._8,0._8,0._8,0._8)

      ! closest sample
      CALL shiftedindex(x1,x2,x3,sx1,sx2,sx3,dx1,dx2,dx3,i,j,k)
      ! rounded coordinates of closest sample
      CALL shiftedcoordinates(i,j,k,sx1,sx2,2*sx3,dx1,dx2,dx3,x,y,z)

      ! no interpolation for node points
      IF ( (abs(x-x1) .lt. EPS*dx1) .and. &
           (abs(y-x2) .lt. EPS*dx2) .and. &
           (abs(z-x3) .lt. EPS*dx3) ) THEN
         lsig=sig(i,j,k)
         RETURN
      END IF

      DO l3=-RANGE,+RANGE
         ! no periodicity in the 3-direction
         IF ((k+l3 .le. 0) .or. (k+l3 .gt. sx3)) CYCLE

         IF (l3 .ge. 0) THEN
            i3p=mod(k-1+l3,sx3)+1
         ELSE
            i3p=mod(sx3+k-1+l3,sx3)+1
         END IF

         DO l2=-RANGE,+RANGE
            IF (l2 .ge. 0) THEN
               i2p=mod(j-1+l2,sx2)+1
            ELSE
               i2p=mod(sx2+j-1+l2,sx2)+1
            END IF

            DO l1=-RANGE,+RANGE
               IF (l1 .ge. 0) THEN
                  i1p=mod(i-1+l1,sx1)+1
               ELSE
                  i1p=mod(sx1+i-1+l1,sx1)+1
               END IF

               weight=sinc(((x+l1*dx1)-x1)/dx1)*dx1 &
                     *sinc(((y+l2*dx2)-x2)/dx2)*dx2 &
                     *sinc(((z+l3*dx3)-x3)/dx3)*dx3

               !weight=gauss((x+l1*dx1)-x1,dx1)*dx1 &
               !      *gauss((y+l2*dx2)-x2,dx2)*dx2 &
               !      *gauss((z+l3*dx3)-x3,dx3)*dx3

               lsig=lsig.plus.(REAL(weight).times.sig(i1p,i2p,i3p))
               sum  =sum  +weight

            END DO
         END DO
      END DO
      IF (sum .gt. 1e-6) lsig=REAL(1._8/sum).times.lsig

    END SUBROUTINE sampletensor

    !-----------------------------------------------
    ! subroutine ref2local
    ! convert reference Cartesian coordinates into
    ! the rotated, local fault coordinates system.
    !-----------------------------------------------
    SUBROUTINE ref2local(x,y,z,xp,yp,zp)
      REAL*8, INTENT(IN) :: x,y,z
      REAL*8, INTENT(OUT) :: xp,yp,zp

      REAL*8 :: x2

      x2 = cstrike*x  -sstrike*y
      xp = cdip   *x2 -sdip   *z
      yp = sstrike*x  +cstrike*y
      zp = sdip   *x2 +cdip   *z

    END SUBROUTINE ref2local

    !-----------------------------------------------
    ! subroutine local2ref
    ! converts a set of coordinates from the rotated
    ! fault-aligned coordinate system into the
    ! reference, Cartesian coordinates system.
    !-----------------------------------------------
    SUBROUTINE local2ref(xp,yp,zp,x,y,z)
      REAL*8, INTENT(IN) :: xp,yp,zp
      REAL*8, INTENT(OUT) :: x,y,z

      REAL*8 :: x2p

      x2p=  cdip*xp+sdip*zp
      x  =  cstrike*x2p+sstrike*yp
      y  = -sstrike*x2p+cstrike*yp
      z  = -sdip*xp    +cdip*zp

    END SUBROUTINE local2ref

  END SUBROUTINE monitorstressfield

  !---------------------------------------------------------------------
  !> function MonitorField
  !! samples a scalar field along a specified planar surface.
  !!
  !! \author sylvain barbot (10-16-07) - original form
  !---------------------------------------------------------------------
  SUBROUTINE monitorfield(x,y,z,L,W,strike,dip,beta, &
       sx1,sx2,sx3,dx1,dx2,dx3,slip,patch)
    INTEGER, INTENT(IN) :: sx1,sx2,sx3
    REAL*8, INTENT(IN) :: x,y,z,L,W,strike,dip,beta,dx1,dx2,dx3
#ifdef ALIGN_DATA
    REAL*4, DIMENSION(sx1+2,sx2,sx3), INTENT(IN) :: slip
#else
    REAL*4, DIMENSION(sx1,sx2,sx3), INTENT(IN) :: slip
#endif
    TYPE(SLIPPATCH_STRUCT), ALLOCATABLE, DIMENSION(:,:), INTENT(OUT) :: patch

    INTEGER :: px2,px3,j2,j3,status
    REAL*8 :: x1,x2,x3,xr,yr,zr,Wp,Lp, &
         cstrike,sstrike,cdip,sdip,value
    TYPE(TENSOR) :: sig0

    sig0=TENSOR(0._8,0._8,0._8,0._8,0._8,0._8)

    cstrike=cos(strike)
    sstrike=sin(strike)
    cdip=cos(dip)
    sdip=sin(dip)

    ! effective tapered dimensions
    Wp=W*(1._8+2._8*beta) ! horizontal dimension for vertical fault
    Lp=L*(1._8+2._8*beta) ! depth for a vertical fault

    px3=fix(Lp/dx3)
    px2=fix(Wp/dx2)

    ALLOCATE(patch(px2+1,px3+1),STAT=status)
    IF (status>0) STOP "could not allocate the slip patches for export"

    DO j3=1,px3+1
       DO j2=1,px2+1

          CALL ref2local(x,y,z,xr,yr,zr)
          
          ! no translation in out of plane direction
          yr=REAL(yr)+REAL((DBLE(j2)-DBLE(px2)/2._8-1._8)*dx2)
          zr=REAL(zr)+REAL((DBLE(j3)-DBLE(px3)/2._8-1._8)*dx3)
          
          CALL local2ref(xr,yr,zr,x1,x2,x3)
          
          ! discard out-of-bound locations
          IF (  (x1 .gt. DBLE(sx1/2-1)*dx1) .or. (x1 .lt. -DBLE(sx1/2)*dx1) &
           .or. (x2 .gt. DBLE(sx2/2-1)*dx2) .or. (x2 .lt. -DBLE(sx2/2)*dx2) &
           .or. (x3 .gt. DBLE(sx3-1)*dx3) .or. (x3 .lt. 0._8)  ) THEN
             value=0._8
          ELSE
             CALL sample(x1,x2,x3,dx1,dx2,dx3,sx1,sx2,sx3,slip,value)
          END IF

          patch(j2,j3)=SLIPPATCH_STRUCT(x1,x2,x3,yr,zr,value,0._8,0._8, &
                                        0._8,0._8,0._8,0._8,sig0)

       END DO
    END DO

  CONTAINS

    !--------------------------------------------------------------
    !> subroutine sample
    !! interpolates the value of a discretized 3-dimensional field
    !! at a subpixel location. method consists in correlating the
    !! 3D field with a delta function filter. the delta function is
    !! approximated with a narrow normalized gaussian.
    !!
    !! \author sylvain barbot (10-17-07) - original form
    !--------------------------------------------------------------
    SUBROUTINE sample(x1,x2,x3,dx1,dx2,dx3,sx1,sx2,sx3,field,value)
      INTEGER, INTENT(IN) :: sx1,sx2,sx3
      REAL*8, INTENT(IN) :: x1,x2,x3,dx1,dx2,dx3
      REAL*8, INTENT(OUT) :: value
#ifdef ALIGN_DATA
    REAL*4, DIMENSION(sx1+2,sx2,sx3), INTENT(IN) :: field
#else
    REAL*4, DIMENSION(sx1,sx2,sx3), INTENT(IN) :: field
#endif
    
      INTEGER :: i,j,k,l1,l2,l3,i1p,i2p,i3p
      INTEGER, PARAMETER :: RANGE=2
      REAL*8 :: sum,weight,x,y,z
      REAL*8, PARAMETER :: EPS=1e-2

      sum=0._8
      value=0._8

      ! closest sample
      CALL shiftedindex(x1,x2,x3,sx1,sx2,sx3,dx1,dx2,dx3,i,j,k)
      ! rounded coordinates of closest sample
      CALL shiftedcoordinates(i,j,k,sx1,sx2,2*sx3,dx1,dx2,dx3,x,y,z)

      ! no interpolation for node points
      IF ( (abs(x-x1) .lt. EPS*dx1) .and. &
           (abs(y-x2) .lt. EPS*dx2) .and. &
           (abs(z-x3) .lt. EPS*dx3) ) THEN
         value=field(i,j,k)
         RETURN
      END IF

      DO l3=-RANGE,+RANGE
         ! no periodicity in the 3-direction
         IF ((k+l3 .le. 0) .or. (k+l3 .gt. sx3)) CYCLE

         IF (l3 .ge. 0) THEN
            i3p=mod(k-1+l3,sx3)+1
         ELSE
            i3p=mod(sx3+k-1+l3,sx3)+1
         END IF

         DO l2=-RANGE,+RANGE
            IF (l2 .ge. 0) THEN
               i2p=mod(j-1+l2,sx2)+1
            ELSE
               i2p=mod(sx2+j-1+l2,sx2)+1
            END IF

            DO l1=-RANGE,+RANGE
               IF (l1 .ge. 0) THEN
                  i1p=mod(i-1+l1,sx1)+1
               ELSE
                  i1p=mod(sx1+i-1+l1,sx1)+1
               END IF

               weight=sinc(((x+l1*dx1)-x1)/dx1)*dx1 &
                     *sinc(((y+l2*dx2)-x2)/dx2)*dx2 &
                     *sinc(((z+l3*dx3)-x3)/dx3)*dx3

               !weight=gauss((x+l1*dx1)-x1,dx1)*dx1 &
               !      *gauss((y+l2*dx2)-x2,dx2)*dx2 &
               !      *gauss((z+l3*dx3)-x3,dx3)*dx3

               value=value+weight*field(i1p,i2p,i3p)
               sum  =sum  +weight

            END DO
         END DO
      END DO
      IF (sum .gt. 1e-6) value=value/sum

    END SUBROUTINE sample

    !-----------------------------------------------
    ! subroutine ref2local
    ! convert reference Cartesian coordinates into
    ! the rotated, local fault coordinates system.
    !-----------------------------------------------
    SUBROUTINE ref2local(x,y,z,xp,yp,zp)
      REAL*8, INTENT(IN) :: x,y,z
      REAL*8, INTENT(OUT) :: xp,yp,zp

      REAL*8 :: x2

      x2 = cstrike*x  -sstrike*y
      xp = cdip   *x2 -sdip   *z
      yp = sstrike*x  +cstrike*y
      zp = sdip   *x2 +cdip   *z

    END SUBROUTINE ref2local

    !-----------------------------------------------
    ! subroutine local2ref
    ! converts a set of coordinates from the rotated
    ! fault-aligned coordinate system into the
    ! reference, Cartesian coordinates system.
    !-----------------------------------------------
    SUBROUTINE local2ref(xp,yp,zp,x,y,z)
      REAL*8, INTENT(IN) :: xp,yp,zp
      REAL*8, INTENT(OUT) :: x,y,z

      REAL*8 :: x2p

      x2p=  cdip*xp+sdip*zp
      x  =  cstrike*x2p+sstrike*yp
      y  = -sstrike*x2p+cstrike*yp
      z  = -sdip*xp    +cdip*zp

    END SUBROUTINE local2ref

  END SUBROUTINE monitorfield

  !-----------------------------------------------------------------
  ! subroutine FieldAdd
  ! computes in place the sum of two scalar fields
  !
  !   u = c1 * u + c2 * v
  !
  ! the function is useful to add fields of different sizes.
  !
  ! sylvain barbot (07/27/07) - original form
  !-----------------------------------------------------------------
  SUBROUTINE fieldadd(u,v,sx1,sx2,sx3,c1,c2)
    INTEGER, INTENT(IN) :: sx1,sx2,sx3
    REAL*4, INTENT(INOUT), DIMENSION(sx1,sx2,sx3) :: u
    REAL*4, INTENT(IN), DIMENSION(sx1,sx2,sx3) :: v
    REAL*4, INTENT(IN), OPTIONAL :: c1,c2

    IF (PRESENT(c1)) THEN
       IF (PRESENT(c2)) THEN
          u=c1*u+c2*v
       ELSE
          u=c1*u+v
       END IF
    ELSE
       IF (PRESENT(c2)) THEN
          u=u+c2*v
       ELSE
          u=u+v
       END IF
    END IF

  END SUBROUTINE fieldadd

  !-----------------------------------------------------------------
  ! subroutine FieldRep
  !
  !   u = c1 * v
  !
  ! the function is useful to add fields of different sizes.
  !
  ! sylvain barbot (07/27/07) - original form
  !-----------------------------------------------------------------
  SUBROUTINE fieldrep(u,v,sx1,sx2,sx3,c1)
    INTEGER, INTENT(IN) :: sx1,sx2,sx3
    REAL*4, INTENT(INOUT), DIMENSION(sx1,sx2,sx3) :: u
    REAL*4, INTENT(IN), DIMENSION(sx1,sx2,sx3) :: v
    REAL*4, INTENT(IN), OPTIONAL :: c1

    IF (PRESENT(c1)) THEN
       u=u+c1*v
    ELSE
       u=v
    END IF
    
  END SUBROUTINE fieldrep

  !-----------------------------------------------------------------
  ! subroutine SliveAdd
  ! computes in place the sum of two scalar fields
  !
  !   u = c1 * u + c2 * v
  !
  ! the function is useful to add fields of different sizes.
  !
  ! sylvain barbot (10/24/08) - original form
  !-----------------------------------------------------------------
  SUBROUTINE sliceadd(u,v,sx1,sx2,sx3,index,c1,c2)
    INTEGER, INTENT(IN) :: sx1,sx2,sx3,index
    REAL*4, INTENT(INOUT), DIMENSION(sx1,sx2) :: u
    REAL*4, INTENT(IN), DIMENSION(sx1,sx2,sx3) :: v
    REAL*4, INTENT(IN), OPTIONAL :: c1,c2

    IF (PRESENT(c1)) THEN
       IF (PRESENT(c2)) THEN
          u=c1*u+c2*v(:,:,index)
       ELSE
          u=c1*u+v(:,:,index)
       END IF
    ELSE
       IF (PRESENT(c2)) THEN
          u=u+c2*v(:,:,index)
       ELSE
          u=u+v(:,:,index)
       END IF
    END IF

  END SUBROUTINE sliceadd

  !-----------------------------------------------------------------
  !> subroutine TensorFieldAdd
  !! computes the linear combination of two tensor fields
  !!
  !!     t1 = c1 * t1 + c2 * t2
  !!
  !! where t1 and t2 are two tensor fields and c1 and c2 are scalars.
  !! only tensor field t1 is modified.
  !
  ! sylvain barbot (07/27/07) - original form
  !-----------------------------------------------------------------
  SUBROUTINE tensorfieldadd(t1,t2,sx1,sx2,sx3,c1,c2)
    INTEGER, INTENT(IN) :: sx1,sx2,sx3
    TYPE(TENSOR), INTENT(INOUT), DIMENSION(sx1,sx2,sx3) :: t1
    TYPE(TENSOR), INTENT(IN), DIMENSION(sx1,sx2,sx3) :: t2
    REAL*4, INTENT(IN), OPTIONAL :: c1,c2

    INTEGER :: i1,i2,i3

    IF (PRESENT(c1)) THEN
       IF (PRESENT(c2)) THEN
          IF (0._4 .eq. c1) THEN
             IF (0._4 .eq. c2) THEN
                DO 05 i3=1,sx3; DO 05 i2=1,sx2; DO 05 i1=1,sx1
                   t1(i1,i2,i3)=TENSOR(0._4,0._4,0._4,0._4,0._4,0._4)
05                 CONTINUE
             ELSE
                DO 10 i3=1,sx3; DO 10 i2=1,sx2; DO 10 i1=1,sx1
                   t1(i1,i2,i3)=c2 .times. t2(i1,i2,i3)
10                 CONTINUE
                END IF
          ELSE
             DO 20 i3=1,sx3; DO 20 i2=1,sx2; DO 20 i1=1,sx1
                t1(i1,i2,i3)=(c1 .times. t1(i1,i2,i3)) .plus. &
                             (c2 .times. t2(i1,i2,i3))
20           CONTINUE
          END IF
       ELSE
          DO 30 i3=1,sx3; DO 30 i2=1,sx2; DO 30 i1=1,sx1
             t1(i1,i2,i3)=(c1 .times. t1(i1,i2,i3)) .plus. t2(i1,i2,i3)
30           CONTINUE
       END IF
    ELSE
       IF (PRESENT(c2)) THEN
          DO 40 i3=1,sx3; DO 40 i2=1,sx2; DO 40 i1=1,sx1
             t1(i1,i2,i3)=t1(i1,i2,i3) .plus. (c2 .times. t2(i1,i2,i3))
40        CONTINUE
       ELSE
          DO 50 i3=1,sx3; DO 50 i2=1,sx2; DO 50 i1=1,sx1
             t1(i1,i2,i3)=t2(i1,i2,i3) .plus. t2(i1,i2,i3)
50        CONTINUE
       END IF
    END IF

  END SUBROUTINE tensorfieldadd


  !-----------------------------------------------------------------
  ! subroutine TensorIntegrate
  ! computes a numercial integration with numerical viscosity
  !
  !    T^(n+1)_i = (T^n_(i-1)+T^n_(i+1))/2 + dt * S^n_i
  !
  ! instead of
  !
  !    T^(n+1)_i = T^n_i + dt * S^n_i
  !
  ! implementation is just generalized for a 3-dimensional field.
  !
  ! sylvain barbot (07/27/07) - original form
  !-----------------------------------------------------------------
  SUBROUTINE tensorintegrate(T,S,sx1,sx2,sx3,dt)
    INTEGER, INTENT(IN) :: sx1,sx2,sx3
    TYPE(TENSOR), INTENT(INOUT), DIMENSION(sx1,sx2,sx3) :: T
    TYPE(TENSOR), INTENT(IN), DIMENSION(sx1,sx2,sx3) :: S
    REAL*8, INTENT(IN) :: dt

    INTEGER :: i1,i2,i3,i1m,i2m,i3m,i1p,i2p,i3p

    DO i3=1,sx3
       i3m=mod(sx3+i3-2,sx3)+1
       i3p=mod(i3,sx3)+1
       DO i2=1,sx2
          i2m=mod(sx2+i2-2,sx2)+1
          i2p=mod(i2,sx2)+1
          DO i1=1,sx1
             i1m=mod(sx1+i1-2,sx1)+1
             i1p=mod(i1,sx1)+1
             
             T(i1,i2,i3)=( &
                  (1._4/6._4) .times. (T(i1m,i2,i3) .plus. T(i1p,i2,i3) &
                  .plus. T(i1,i2m,i3) .plus. T(i1,i2p,i3) &
                  .plus. T(i1,i2,i3m) .plus. T(i1,i2,i3p))) &
                  .plus. &
                  (REAL(dt) .times. S(i1,i2,i3))
          END DO
       END DO
    END DO

  END SUBROUTINE tensorintegrate

  !---------------------------------------------------------------------
  !> subroutine coordinates computes the xi coordinates from the
  !! array index and sampling interval
  !---------------------------------------------------------------------
  SUBROUTINE coordinates(i1,i2,i3,sx1,sx2,sx3,dx1,dx2,dx3,x1,x2,x3)
    INTEGER, INTENT(IN) :: i1,i2,i3,sx1,sx2,sx3
    REAL*8, INTENT(IN) :: dx1,dx2,dx3
    REAL*8, INTENT(OUT) :: x1,x2,x3
    
    x1=DBLE(i1-sx1/2-1)*dx1
    x2=DBLE(i2-sx2/2-1)*dx2
    x3=DBLE(i3-sx3/2-1)*dx3
  END SUBROUTINE coordinates

  !---------------------------------------------------------------------
  !> subroutine ShiftedCoordinates
  !! computes the xi coordinates from the array index and sampling
  !! interval assuming data is order like fftshift.
  !!
  !! \author sylvain barbot (07/31/07) - original form
  !---------------------------------------------------------------------
  SUBROUTINE shiftedcoordinates(i1,i2,i3,sx1,sx2,sx3,dx1,dx2,dx3,x1,x2,x3)
    INTEGER, INTENT(IN) :: i1,i2,i3,sx1,sx2,sx3
    REAL*8, INTENT(IN) :: dx1,dx2,dx3
    REAL*8, INTENT(OUT) :: x1,x2,x3

    IF (i1 .LE. sx1/2) THEN
       x1=DBLE(i1-1)*dx1
    ELSE
       x1=DBLE(i1-sx1-1)*dx1
    END IF
    IF (i2 .LE. sx2/2) THEN
       x2=DBLE(i2-1)*dx2
    ELSE
       x2=DBLE(i2-sx2-1)*dx2
    END IF
    IF (i3 .LE. sx3/2) THEN
       x3=DBLE(i3-1)*dx3
    ELSE
       x3=DBLE(i3-sx3-1)*dx3
    END IF

  END SUBROUTINE shiftedcoordinates

  !----------------------------------------------------------------------
  !> subroutine ShiftedIndex
  !! returns the integer index corresponding to the specified coordinates
  !! assuming the data are ordered following fftshift. input coordinates
  !! are assumed bounded -sx/2 <= x <= sx/2-1. out of bound input
  !! purposefully triggers a fatal error. in the x3 direction, coordinates
  !! are assumed bounded by 0 <= x3 <= (sx3-1)*dx3
  !!
  !! CALLED BY:
  !!   monitorfield/sample
  !!
  !! \author sylvain barbot (07/31/07) - original form
  !----------------------------------------------------------------------
  SUBROUTINE shiftedindex(x1,x2,x3,sx1,sx2,sx3,dx1,dx2,dx3,i1,i2,i3)
    REAL*8, INTENT(IN) :: x1,x2,x3,dx1,dx2,dx3
    INTEGER, INTENT(IN) :: sx1,sx2,sx3
    INTEGER, INTENT(OUT) :: i1,i2,i3

    IF (x1 .gt.  DBLE(sx1/2-1)*dx1) THEN
       WRITE_DEBUG_INFO
       WRITE (0,'("x1=",ES9.2E2,"; boundary at x1=",ES9.2E2)') x1, DBLE(sx1/2)*dx1
       STOP "ShiftedIndex:invalid x1 coordinates (x1 too large)"
    END IF
    IF (x1 .lt. -DBLE(sx1/2)*dx1  ) THEN
       WRITE_DEBUG_INFO
       WRITE (0,'("x1=",ES9.2E2,"; boundary at x1=",ES9.2E2)') x1, -DBLE(sx1/2)*dx1
       STOP "ShiftedIndex:coordinates out of range (-x1 too large)"
    END IF
    IF (x2 .gt.  DBLE(sx2/2-1)*dx2) THEN
       WRITE_DEBUG_INFO
       WRITE (0,'("x2=",ES9.2E2,"; boundary at x2=",ES9.2E2)') x2, DBLE(sx2/2)*dx2
       STOP "ShiftedIndex:invalid x2 coordinates (x2 too large)"
    END IF
    IF (x2 .lt. -DBLE(sx2/2)*dx2  ) THEN
       WRITE_DEBUG_INFO
       WRITE (0,'("x2=",ES9.2E2,"; boundary at x2=",ES9.2E2)') x2, -DBLE(sx2/2)*dx2
       STOP "ShiftedIndex:coordinates out of range (-x2 too large)"
    END IF
    IF (x3 .gt.  DBLE(sx3-1)*dx3) THEN
       WRITE_DEBUG_INFO
       STOP "ShiftedIndex:invalid x3 coordinates (x3 too large)"
    END IF
    IF (x3 .lt.  0              )   THEN
       WRITE (0,'("x3=",ES9.2E2)') x3
       STOP "ShiftedIndex:coordinates out of range (x3 negative)"
    END IF

    i1=MOD(sx1+fix(x1/dx1),sx1)+1
    i2=MOD(sx2+fix(x2/dx2),sx2)+1
    i3=fix(x3/dx3)+1

  END SUBROUTINE shiftedindex

  !-----------------------------------------------------------------
  ! subroutine ExportSlice
  ! computes the value of a scalar field at a horizontal plane.
  ! the field if shifted such as the (0,0) coordinate is in the 
  ! middle of the array at (sx1/2+1,sx2/2+1).
  !
  ! sylvain barbot (07/09/07) - original form
  !-----------------------------------------------------------------
  SUBROUTINE exportslice(field,odepth,dx1,dx2,dx3,s)
    REAL*4, INTENT(IN), DIMENSION(:,:,:) :: field
    REAL*8, INTENT(IN) :: dx1,dx2,dx3,odepth
    REAL*4, INTENT(OUT), DIMENSION(:,:) :: s
    
    INTEGER :: i1,i2,i3,sx1,sx2,sx3
    REAL*8 :: k3
    COMPLEX(KIND=8), PARAMETER :: i=CMPLX(0._8,pi2)
    COMPLEX(KIND=8) :: sum,exp3
    REAL*4 :: exp1,exp2
  
    sx1=SIZE(field,1)-2
    sx2=SIZE(field,2)
    sx3=SIZE(field,3)
    
    s=0
    DO i3=1,sx3
       CALL wavenumber3(i3,sx3,dx3,k3)
       exp3=exp(i*k3*odepth)
       DO i2=1,sx2
          DO i1=1,sx1/2+1
             sum=CMPLX(field(2*i1-1,i2,i3),field(2*i1,i2,i3))*exp3
             s(2*i1-1:2*i1,i2)=REAL(s(2*i1-1:2*i1,i2)+(/REAL(sum),AIMAG(sum)/))
          END DO
       END DO
    END DO
    s=REAL(s/(sx3*dx3))
    
    !fftshift
    DO i2=1,sx2
       IF (i2 < sx2/2+1) THEN
          exp2= (i2-1._4)
       ELSE
          exp2=-(sx2-i2+1._4)
       END IF
       DO i1=1,sx1/2+1
          exp1=i1-1._4
          sum=CMPLX(s(2*i1-1,i2),s(2*i1,i2))*((-1._4)**(exp1+exp2))
          s(2*i1-1:2*i1,i2)=REAL((/REAL(sum),AIMAG(sum)/))
       END DO
    END DO
    CALL fft2(s,sx1,sx2,dx1,dx2,FFT_INVERSE)
    
  END SUBROUTINE exportslice

  !-----------------------------------------------------------------
  !> subroutine ExportSpatial
  !! transfer a horizontal layer from array 'data' to smaller array
  !! 'p' and shift center position so that coordinates (0,0) are in
  !! center of array 'p'. optional parameter 'doflip' generates
  !! output compatible with grd binary format.
  !
  ! sylvain barbot (07/09/07) - original form
  !                (03/19/08) - compatibility with grd output
  !-----------------------------------------------------------------
  SUBROUTINE exportspatial(data,sx1,sx2,p,doflip)
    INTEGER, INTENT(IN) :: sx1,sx2
#ifdef ALIGN_DATA
    REAL*4, INTENT(IN), DIMENSION(sx1+2,sx2) :: data
#else
    REAL*4, INTENT(IN), DIMENSION(sx1,sx2) :: data
#endif
    REAL*4, INTENT(OUT), DIMENSION(:,:) :: p
    LOGICAL, INTENT(IN), OPTIONAL :: doflip

    INTEGER :: i1,i2,i1s,i2s
    LOGICAL :: flip

    IF (PRESENT(doflip)) THEN
       flip=doflip
    ELSE
       flip=.false.
    END IF

    DO i2=1,sx2
       IF (i2 .LE. sx2/2) THEN
          i2s=sx2/2+i2
       ELSE
          i2s=i2-sx2/2
       END IF
       DO i1=1,sx1
          IF (i1 .LE. sx1/2) THEN
             i1s=sx1/2+i1
          ELSE
             i1s=i1-sx1/2
          END IF

          IF (flip) THEN
             p(i2s,sx1-i1s+1)=data(i1,i2)
          ELSE
             p(i1s,i2s)=data(i1,i2)
          END IF

       END DO
    END DO

  END SUBROUTINE exportspatial

END MODULE elastic3d
