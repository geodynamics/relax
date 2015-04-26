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

MODULE friction3d

  USE types
  USE elastic3d

  IMPLICIT NONE

  REAL*8, PRIVATE, PARAMETER :: pi   = 3.141592653589793115997963468544185161_8
  REAL*8, PRIVATE, PARAMETER :: pi2  = 6.28318530717958623199592693708837032318_8
  REAL*8, PRIVATE, PARAMETER :: pid2 = 1.57079632679489655799898173427209258079_8

CONTAINS

  !-----------------------------------------------------------------
  !> subroutine FrictionEigenStress
  !! compute the eigen-stress (forcing moment) to be relaxed by
  !! rate-dependent inelastic deformation in the case of a frictional
  !! surface:
  !!
  !!        sigma^i = C:F:sigma
  !!
  !! where C is the elastic moduli tensor, F is the heterogeneous
  !! fluidity moduli tensor and sigma is the instantaneous stress
  !! tensor. for a frictional surface, the eigenstrain-rate is given
  !! by
  !!
  !!        epsilon^i^dot = F:sigma = gamma^dot R
  !!
  !! where gamma^dot is the slip rate (a scalar) and R is the
  !! deviatoric, symmetric, and unitary, tensor:
  !!
  !!        R_ij = 1/2 ( t_i n_j + t_j n_i ) / sqrt( t_j t_j )
  !!
  !! where the shear traction t_i is the projection of the traction
  !! vector on the plane surface. the strain amplitude is given by
  !!
  !!        gamma^dot = H( t_j r_j ) 2 vo sinh( taus / (t_c )
  !!
  !! where taus is the effective shear on the fault plane,
  !!
  !!        taus = tau + mu*sigma
  !!
  !! where tau is the shear and sigma the normal stress. tau and sigma
  !! assumed to be the co-seismic change only, not the absolute
  !! stress. vo is a reference slip velocity, and t_c, the critical
  !! stress, corresponds to (a-b)*sigma in the framework of rate-and-
  !! state friction. the effective viscosity eta* and the fluidity
  !!
  !!        eta* = tau / gamma^dot
  !!        fluidity = 1 / eta*
  !!
  !! are used to compute the optimal time-step. H( x ) is the 
  !! Heaviside function and r_i is the rake vector. I impose
  !! gamma^dot to be zero is t_j r_j < 0. This constraint is
  !! enforced to ensure that no back slip occurs on faults.
  !!
  !! \author sylvain barbot (07/24/07) - original form
  !!                        (02/28/11) - add constraints on the direction
  !!                                     of afterslip
  !-----------------------------------------------------------------
  SUBROUTINE frictioneigenstress(x,y,z,L,W,strike,dip,rake,beta, &
       sig,prestress,mu,structure,sx1,sx2,sx3,dx1,dx2,dx3,moment,maxwelltime,vel)
    INTEGER, INTENT(IN) :: sx1,sx2,sx3
    REAL*8, INTENT(IN) :: mu,dx1,dx2,dx3,x,y,z,L,W,strike,dip,rake,beta
    TYPE(LAYER_STRUCT), DIMENSION(:), INTENT(IN) :: structure
    TYPE(TENSOR_LAYER_STRUCT), DIMENSION(:), INTENT(IN) :: prestress
    TYPE(TENSOR), INTENT(IN), DIMENSION(sx1,sx2,sx3) :: sig
    TYPE(TENSOR), INTENT(INOUT), DIMENSION(sx1,sx2,sx3) :: moment
    REAL*8, OPTIONAL, INTENT(INOUT) :: maxwelltime
#ifdef ALIGN_DATA
    REAL*4, INTENT(OUT), DIMENSION(sx1+2,sx2,sx3), OPTIONAL :: vel
#else
    REAL*4, INTENT(OUT), DIMENSION(sx1,sx2,sx3), OPTIONAL :: vel
#endif

    INTEGER :: i1,i2,i3
    TYPE(TENSOR) :: s
    REAL*8, DIMENSION(3) :: t,ts,n,r
    REAL*8 :: vo,tauc,taun,taus,gammadot,impulse, &
         friction,tau,scaling,cohesion
    REAL*8 :: x1,x2,x3,x1s,x2s,x3s,x1i,x3i, &
         cstrike,sstrike,cdip,sdip,cr,sr,x2r,&
         temp1,temp2,temp3,sourc,image,xr,yr,zr,Wp,Lp,dum
    REAL*4 :: tm

    IF (PRESENT(maxwelltime)) THEN
       tm=REAL(maxwelltime)
       i1=1
    ELSE
       tm=1e30
       i1=0
    END IF
#ifdef USING_CUDA
    CALL cufrictioneigenstress (%VAL(x), %VAL(y), %VAL(z), %VAL(L), %VAL(W),%VAL(strike), &
         %VAL(dip), %VAL(rake), %VAL(beta), %VAL(mu), structure, %VAL(sx1),%VAL(sx2), &
         %VAL(sx3), %VAL(dx1), %VAL(dx2), %VAL(dx3), %VAL(i1), tm, moment, sig)
#else
    
    ! delta function scaling
    scaling=sqrt(pi2)*dx1

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
    
    ! surface normal vector components
    n(1)=+cdip*cstrike
    n(2)=-cdip*sstrike
    n(3)=-sdip
             
    ! rake vector component
    r(1)=sstrike*cr+cstrike*sdip*sr
    r(2)=cstrike*cr-sstrike*sdip*sr
    r(3)=cdip*sr

    DO i3=1,sx3
       x3=DBLE(i3-1)*dx3
       IF ((abs(x3-z).gt.Lp) .and. (abs(x3+z).gt.Lp)) CYCLE

       vo=structure(i3)%gammadot0
       tauc=structure(i3)%stressexponent
       friction=structure(i3)%friction
       cohesion=structure(i3)%cohesion
       
       DO i2=1,sx2
          DO i1=1,sx1
             CALL shiftedcoordinates(i1,i2,i3,sx1,sx2,sx3, &
                  dx1,dx2,dx3,x1,x2,dum)
             IF ((ABS(x1-x).gt.MAX(Wp,Lp)) .OR.  (ABS(x2-y).gt.MAX(Wp,Lp))) CYCLE
             
             x2r= cstrike*x1-sstrike*x2
             x1s= cdip*x2r-sdip*x3
             x1i= cdip*x2r+sdip*x3
             IF ((ABS(x1s-xr).GT.7.01_8*dx1).AND.(ABS(x1i-xr).GT.7.01_8*dx1)) CYCLE
             x2s= sstrike*x1+cstrike*x2
             x3s= sdip*x2r+cdip*x3
             x3i=-sdip*x2r+cdip*x3

             ! integrate at depth and along strike with raised cosine taper
             ! and shift sources to x,y,z coordinate
             temp1=gauss(x1s-xr,dx1)
             temp2=omega((x2s-yr)/W,beta)
             temp3=omega((x3s-zr)/L,beta)
             sourc=temp1*temp2*temp3

             temp1=gauss(x1i-xr,dx1)
             temp3=omega((x3i+zr)/L,beta)
             image=temp1*temp2*temp3

             impulse=sourc+image

             ! traction = sigma . n
             s=sig(i1,i2,i3)
             t=s .tdot. n

             ! signed normal component
             taun=SUM(t*n)

             ! absolute value of shear component
             ts=t-taun*n
             taus=SQRT(SUM(ts*ts))

             ! effective shear stress on fault plane
             tau=MAX(0.d0,taus+friction*taun-cohesion)

             ! rake direction test only if | rake | < 3*Pi
             IF (SUM(ts*r).LT.0.d0 .AND. ABS(rake).LT.pi2*1.5d0) CYCLE

             ! warning for wrong input
             IF ((tau/tauc) .gt. 20) THEN
                WRITE_DEBUG_INFO
                WRITE (0,'("------------------------------------------")')
                WRITE (0,'("wrong value of (a-b)sigma gives rise to")')
                WRITE (0,'("(a - b) * sigma       = ",ES11.3E2)') tauc
                WRITE (0,'("tau                   = ",ES11.3E2)') tau
                WRITE (0,'("tau_s                 = ",ES11.3E2)') taus
                WRITE (0,'("tau_n                 = ",ES11.3E2)') taun
                WRITE (0,'("tau / ((a - b) sigma) = ",ES11.3E2)') tau/tauc
                WRITE (0,'("------------------------------------------")')
                STOP 5
             END IF

             ! shear traction direction
             ts=ts/taus

             ! deviatoric strain rate
             gammadot=vo*2._8*my_sinh(tau/tauc)

             tm=MIN(tm,REAL(tau/mu/gammadot*(MIN(L,W)/sqrt(dx1*dx3))))

             ! provide the strain-rate on request
             IF (PRESENT(vel)) THEN
                vel(i1,i2,i3)=REAL(vel(i1,i2,i3)+gammadot*impulse*scaling)
             END IF

             ! deviatoric strain
             moment(i1,i2,i3)=moment(i1,i2,i3) .plus. &
                  (ts .sdyad. ((2._8*mu*impulse*gammadot)*n))

          END DO
       END DO
    END DO
#endif

    IF (PRESENT(maxwelltime)) maxwelltime=MIN(tm,maxwelltime)

  END SUBROUTINE frictioneigenstress

  !---------------------------------------------------------------------
  !> function MonitorFriction
  !! samples a scalar field along a specified planar surface.
  !!
  !! input:
  !! @param x,y,z       coordinates of the creeping segment
  !! @param L           dimension of segment in the depth direction
  !! @param W           dimension of segment in the strike direction
  !! @param beta        smoothing factor
  !! @param sx1,2,3     dimension of the stress tensor array
  !! @param dx1,2,3     sampling size
  !! @param sig         stress tensor array
  !! @param structure   frictional properties as a function of depth
  !!
  !! output:
  !! @param patch       list of strike- and dip-slip as a function of position
  !!                    on the fault.     
  !! 
  !! \author sylvain barbot (10-16-07) - original form
  !---------------------------------------------------------------------
  SUBROUTINE monitorfriction(x,y,z,L,W,strike,dip,rake,beta, &
       sx1,sx2,sx3,dx1,dx2,dx3,sig,prestress,structure,patch)
    INTEGER, INTENT(IN) :: sx1,sx2,sx3
    REAL*8, INTENT(IN) :: x,y,z,L,W,strike,rake,dip,beta,dx1,dx2,dx3
    TYPE(TENSOR), DIMENSION(sx1,sx2,sx3), INTENT(IN) :: sig
    TYPE(SLIPPATCH_STRUCT), ALLOCATABLE, DIMENSION(:,:), INTENT(INOUT) :: patch
    TYPE(LAYER_STRUCT), DIMENSION(:), INTENT(IN) :: structure
    TYPE(TENSOR_LAYER_STRUCT), DIMENSION(:), INTENT(IN) :: prestress

    INTEGER :: i1,i2,i3,px2,px3,j2,j3
    REAL*8 :: cstrike,sstrike,cdip,sdip,cr,sr
    REAL*8 :: vo,tauc,taun,taus, &
         friction,tau,cohesion
    REAL*8 :: x1,x2,x3,xr,yr,zr
    TYPE(TENSOR) :: s
    REAL*8, DIMENSION(3) :: t,ts,n,sv,dv,r

    cstrike=cos(strike)
    sstrike=sin(strike)
    cdip=cos(dip)
    sdip=sin(dip)
    cr=cos(rake)
    sr=sin(rake)

    ! strike direction vector
    sv=(/ sstrike, cstrike, 0._8 /)

    ! dip direction vector
    dv=(/ -cstrike*sdip, +sstrike*sdip, -cdip /)

    ! number of samples in the dip and strike direction
    px2=SIZE(patch,1)
    px3=SIZE(patch,2)

    ! surface normal vector components
    n(1)=+cdip*cstrike
    n(2)=-cdip*sstrike
    n(3)=-sdip
             
    ! rake vector component
    r(1)=sstrike*cr+cstrike*sdip*sr
    r(2)=cstrike*cr-sstrike*sdip*sr
    r(3)=cdip*sr

    ! loop in the dip direction
    DO j3=1,px3
       ! loop in the strike direction
       DO j2=1,px2

          CALL ref2local(x,y,z,xr,yr,zr)
          
          ! no translation in out of plane direction
          yr=REAL(yr)+REAL((DBLE(j2)-DBLE(px2)/2._8-1._8)*dx2)
          zr=REAL(zr)+REAL((DBLE(j3)-DBLE(px3)/2._8-1._8)*dx3)
          
          CALL local2ref(xr,yr,zr,x1,x2,x3)

          ! initialize zero slip velocity
          patch(j2,j3)%x1=x1
          patch(j2,j3)%x2=x2
          patch(j2,j3)%x3=x3
          patch(j2,j3)%lx=yr
          patch(j2,j3)%lz=zr

          ! discard out-of-bound locations
          IF (  (x1 .GT. DBLE(sx1/2-1)*dx1) .OR. (x1 .LT. -DBLE(sx1/2)*dx1) &
           .OR. (x2 .GT. DBLE(sx2/2-1)*dx2) .OR. (x2 .LT. -DBLE(sx2/2)*dx2) &
           .OR. (x3 .GT. DBLE(sx3-1)*dx3) .OR. (x3 .LT. 0._8)  ) CYCLE

          ! evaluates instantaneous creep velocity
          CALL shiftedindex(x1,x2,x3,sx1,sx2,sx3,dx1,dx2,dx3,i1,i2,i3)

          ! retrieve friction parameters
          vo=structure(i3)%gammadot0
          tauc=structure(i3)%stressexponent
          friction=structure(i3)%friction
          cohesion=structure(i3)%cohesion
       
          ! traction = sigma . n
          s=sig(i1,i2,i3)
          patch(j2,j3)%sig=s
          t=s .tdot. n

          ! signed normal component
          taun=SUM(t*n)

          ! absolute value of shear component
          ts=t-taun*n
          taus=SQRT(SUM(ts*ts))
             
          ! effective shear stress on fault plane
          tau=MAX(0.d0,taus+friction*taun-cohesion)

          ! rake direction test only if | rake | < 3*Pi
          IF (SUM(ts*r).LT.0.d0 .AND. ABS(rake).LT.pi2*1.5d0) CYCLE

          ! shear stress
          patch(j2,j3)%taus=taus

          ! creep rate
          patch(j2,j3)%v=vo*2._8*my_sinh(tau/tauc)

          ! shear traction direction
          ts=ts/taus

          ! strike-direction creep rate
          patch(j2,j3)%vss=patch(j2,j3)%v*SUM(ts*sv)

          ! dip-direction creep rate
          patch(j2,j3)%vds=patch(j2,j3)%v*SUM(ts*dv)

       END DO
    END DO

  CONTAINS

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

  END SUBROUTINE monitorfriction

  !---------------------------------------------------------------------
  !> function FrictionAdd
  !! update the cumulative slip of a creeping segment
  !!
  !! \author sylvain barbot (04-01-12) - original form
  !---------------------------------------------------------------------
  SUBROUTINE frictionadd(np,n,dt)
    INTEGER, INTENT(IN) :: np
    TYPE(PLANE_STRUCT), INTENT(INOUT), DIMENSION(np) :: n
    REAL*8, INTENT(IN) :: dt

    INTEGER :: px2,px3,j2,j3,k

    ! number of samples in the dip and strike direction
    DO k=1,np
       px2=SIZE(n(k)%patch,1)
       px3=SIZE(n(k)%patch,2)

       ! loop in the dip direction
       DO j3=1,px3
          ! loop in the strike direction
          DO j2=1,px2
             ! cumulative strike-direction creep
             n(k)%patch(j2,j3)%ss=n(k)%patch(j2,j3)%ss+dt*n(k)%patch(j2,j3)%vss

             ! cumulative dip-direction creep
             n(k)%patch(j2,j3)%ds=n(k)%patch(j2,j3)%ds+dt*n(k)%patch(j2,j3)%vds

             ! cumulative creep
             n(k)%patch(j2,j3)%slip=(n(k)%patch(j2,j3)%ds**2+n(k)%patch(j2,j3)%ss**2)**0.5
          END DO
       END DO

    END DO

  END SUBROUTINE frictionadd

END MODULE friction3d
