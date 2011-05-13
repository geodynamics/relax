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

MODULE viscoelastic3d

  USE elastic3d

  IMPLICIT NONE

#include "include.f90"

  REAL*8, PRIVATE, PARAMETER :: pi   = 3.141592653589793115997963468544185161_8
  REAL*8, PRIVATE, PARAMETER :: pi2  = 6.28318530717958623199592693708837032318_8
  REAL*8, PRIVATE, PARAMETER :: pid2 = 1.57079632679489655799898173427209258079_8
    
CONTAINS

  !-----------------------------------------------------------------
  !> subroutine ViscoElasticDeviatoricStress
  !! computes the instantaneous deviatoric stress tensor sigma_ij'
  !!
  !!  sigma_ij' = 2*mu*(-delta_ij epsilon_kk/3 + epsilon_ij) - tau_ij 
  !!
  !! such as
  !! 
  !!  sigma_kk'= 0
  !!
  !! where tau_ij is a second-order deviatoric symmetric tensor 
  !! that integrates the history of the relaxed stress. strain is
  !! estimated using a centered finite difference derivative.
  !!
  !! \author sylvain barbot (07/07/07) - original form
  !-----------------------------------------------------------------
  SUBROUTINE viscoelasticdeviatoricstress(mu,u1,u2,u3,tau,&
       dx1,dx2,dx3,sx1,sx2,sx3,sig)
    INTEGER, INTENT(IN) :: sx1,sx2,sx3
    REAL*8, INTENT(IN) :: mu,dx1,dx2,dx3
    REAL*4, INTENT(IN), DIMENSION(sx1,sx2,sx3) :: u1,u2,u3
    TYPE(TENSOR), INTENT(IN),  DIMENSION(:,:,:) :: tau
    TYPE(TENSOR), INTENT(OUT), DIMENSION(:,:,:) :: sig
    
    TYPE(TENSOR) :: s
    INTEGER :: i1,i2,i3,i1p,i2p,i3p,i1m,i2m,i3m
    REAL*8 :: epskk,px1,px2,px3

    px1=dx1*2._8
    px2=dx2*2._8
    px3=dx3*2._8
    
    ! space domain with finite difference scheme
    DO i3=1,sx3
       ! wrap around neighbor
       i3m=mod(sx3+i3-2,sx3)+1
       i3p=mod(i3,sx3)+1
       DO i2=1,sx2
          i2m=mod(sx2+i2-2,sx2)+1
          i2p=mod(i2,sx2)+1
          
          DO i1=1,sx1
             i1m=mod(sx1+i1-2,sx1)+1
             i1p=mod(i1,sx1)+1
             
             ! trace component
             epskk=((u1(i1p,i2,i3)-u1(i1m,i2,i3))/px1+&
                    (u2(i1,i2p,i3)-u2(i1,i2m,i3))/px2+&
                    (u3(i1,i2,i3p)-u3(i1,i2,i3m))/px3)/3._8
             
             s%s11=2._8*mu*( (u1(i1p,i2,i3)-u1(i1m,i2,i3))/px1-epskk )
             s%s12=     mu*( (u1(i1,i2p,i3)-u1(i1,i2m,i3))/px2+ &
                             (u2(i1p,i2,i3)-u2(i1m,i2,i3))/px1)
             s%s13=     mu*( (u1(i1,i2,i3p)-u1(i1,i2,i3m))/px3+ &
                             (u3(i1p,i2,i3)-u3(i1m,i2,i3))/px1)
             s%s22=2._8*mu*( (u2(i1,i2p,i3)-u2(i1,i2m,i3))/px2-epskk )
             s%s23=     mu*( (u2(i1,i2,i3p)-u2(i1,i2,i3m))/px3+ &
                             (u3(i1,i2p,i3)-u3(i1,i2m,i3))/px2)
             s%s33=2._8*mu*( (u3(i1,i2,i3p)-u3(i1,i2,i3m))/px3-epskk )
             
             sig(i1,i2,i3)= s .minus. tau(i1,i2,i3)
             
          END DO
       END DO
    END DO
    
    ! no normal traction at the boundary
    sig(:,:,1)%s13=0
    sig(:,:,1)%s23=0
    sig(:,:,1)%s33=0
    sig(:,:,sx3)%s13=0
    sig(:,:,sx3)%s23=0
    sig(:,:,sx3)%s33=0

  END SUBROUTINE viscoelasticdeviatoricstress

  !-----------------------------------------------------------------
  !> subroutine ViscousEigenstress
  !! computes the moment density rate due to a layered viscoelastic
  !! structure with powerlaw creep
  !!
  !!     d Ei / dt = C:F:sigma'
  !!
  !! where C is the elastic moduli tensor, F is the heterogeneous
  !! fluidity tensor and sigma' is the instantaneous deviatoric 
  !! stress. F is stress dependent (powerlaw creep.)
  !!
  !! \author sylvain barbot (08/30/08) - original form
  !-----------------------------------------------------------------
  SUBROUTINE viscouseigenstress(mu,structure,ductilezones,sig,sx1,sx2,sx3, &
       dx1,dx2,dx3,moment,beta,maxwelltime,gamma)
    REAL*8, INTENT(IN) :: mu,dx1,dx2,dx3,beta
    TYPE(LAYER_STRUCT), DIMENSION(:), INTENT(IN) :: structure
    TYPE(WEAK_STRUCT), DIMENSION(:), INTENT(IN) :: ductilezones
    INTEGER, INTENT(IN) :: sx1,sx2,sx3
    TYPE(TENSOR), INTENT(IN), DIMENSION(sx1,sx2,sx3) :: sig
    TYPE(TENSOR), INTENT(OUT), DIMENSION(sx1,sx2,sx3) :: moment
    REAL*8, OPTIONAL, INTENT(INOUT) :: maxwelltime
#ifdef ALIGN_DATA
    REAL*4, DIMENSION(sx1+2,sx2,sx3), INTENT(OUT), OPTIONAL :: gamma
#else
    REAL*4, DIMENSION(sx1,sx2,sx3), INTENT(OUT), OPTIONAL :: gamma
#endif

    INTEGER :: i1,i2,i3
    TYPE(TENSOR) :: s,R
    TYPE(TENSOR), PARAMETER :: zero = tensor(0._4,0._4,0._4,0._4,0._4,0._4)
    REAL*8 :: gammadot,tau,tauc,gammadot0,power,cohesion,x1,x2,x3,dg0,dum
    REAL*4 :: tm
    
    IF (SIZE(structure,1) .NE. sx3) RETURN

    IF (PRESENT(maxwelltime)) THEN
       tm=REAL(maxwelltime)
    ELSE
       tm=1e30
    END IF

!$omp parallel do private(i1,i2,gammadot0,power,cohesion,s,tau,R,tauc,gammadot,dg0,x1,x2,x3,dum), &
!$omp reduction(MIN:tm)
    DO i3=1,sx3
       power=structure(i3)%stressexponent
       cohesion=structure(i3)%cohesion
       x3=DBLE(i3-1)*dx3

       IF (power .LT. 0.999999_8) THEN 
          WRITE_DEBUG_INFO
          WRITE (0,'("power=",ES9.2E1)') power
          WRITE (0,'("invalid power exponent. interrupting.")')
          STOP 1
       END IF

       DO i2=1,sx2
          DO i1=1,sx1
             ! local coordinates
             CALL shiftedcoordinates(i1,i2,i3,sx1,sx2,sx3, &
                  dx1,dx2,dx3,x1,x2,dum)

             ! depth-dependent fluidity structure             
             gammadot0=structure(i3)%gammadot0

             ! perturbation from isolated viscous zones
             dg0=dgammadot0(ductilezones,x1,x2,x3,beta)

             ! local fluidity structure
             gammadot0=gammadot0+dg0

             IF (1e-9 .GT. gammadot0) CYCLE

             ! local deviatoric stress
             s=tensordeviatoric(sig(i1,i2,i3))
             
             ! s = tau * R
             CALL tensordecomposition(s,tau,R)

             ! effective stress
             tauc=tau-cohesion

             ! cohesion test
             IF (tauc .LE. 1e-9) CYCLE

             ! powerlaw viscosity
             gammadot=gammadot0*(tauc/mu)**power

             ! update moment density forcing
             moment(i1,i2,i3)=moment(i1,i2,i3) .plus. &
                  (REAL(2._8*mu*gammadot) .times. R)

             tm=MIN(tm,tauc/mu/gammadot)

             IF (PRESENT(gamma)) &
                  gamma(i1,i2,i3)=gammadot
             
          END DO
       END DO
    END DO
!$omp end parallel do

    IF (PRESENT(maxwelltime)) maxwelltime=MIN(tm,maxwelltime)

  CONTAINS

    !---------------------------------------------------------
    !> function dgammadot0
    !! evaluates the change of fluidity at position x1,x2,x3
    !! due to the presence of weak ductile zones. the extent
    !! and magnitude of ductile zones is tapered (beta).
    !!
    !! \author sylvain barbot (3/29/10) - original form
    !---------------------------------------------------------
    REAL*8 FUNCTION dgammadot0(zones,x1,x2,x3,beta)
       TYPE(WEAK_STRUCT), INTENT(IN), DIMENSION(:) :: zones
       REAL*8, INTENT(IN) :: x1,x2,x3,beta

       REAL*8 :: dg,x,y,z,L,W,D,strike,dip,LM
       REAL*8 :: cstrike,sstrike,cdip,sdip, &
                 xr,yr,zr,x2r,Wp,Lp,Dp,x1s,x2s,x3s
       INTEGER :: n,i

       ! number of ductile zones
       n=SIZE(zones,1)

       ! default is no change in fluidity
       dgammadot0=0._8

       DO i=1,n
          ! retrieve weak zone geometry
          dg=zones(i)%dgammadot0
          x=zones(i)%x;y=zones(i)%y;z=zones(i)%z
          W=zones(i)%length;L=zones(i)%width;D=zones(i)%thickness
          strike=zones(i)%strike;dip=zones(i)%dip

          ! effective tapered dimensions
          Wp=W*(1._8+2._8*beta)/2._8
          Lp=L*(1._8+2._8*beta)/2._8
          Dp=D*(1._8+2._8*beta)/2._8
          LM=MAX(Wp,Lp,Dp)

          ! check distance from weak zone
          IF ((ABS(x3-z).GT.LM) .OR. &
              (ABS(x1-x).GT.LM) .OR. &
              (ABS(x2-y).GT.LM)) CYCLE

          ! evaluate contribution from weak zone
          cstrike=cos(strike)
          sstrike=sin(strike)
          cdip=cos(dip)
          sdip=sin(dip)

          ! rotate centre coordinates of weak zone
          x2r= cstrike*x  -sstrike*y
          xr = cdip   *x2r-sdip   *z
          yr = sstrike*x  +cstrike*y
          zr = sdip   *x2r+cdip   *z

          x2r= cstrike*x1 -sstrike*x2
          x1s= cdip   *x2r-sdip   *x3
          x2s= sstrike*x1 +cstrike*x2
          x3s= sdip   *x2r+cdip   *x3

          dgammadot0=dgammadot0+omega((x1s-xr)/D,beta) &
                               *omega((x2s-yr)/W,beta) &
                               *omega((x3s-zr)/L,beta)*dg
       END DO

    END FUNCTION dgammadot0

  END SUBROUTINE viscouseigenstress

END MODULE viscoelastic3d
