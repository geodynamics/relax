MODULE green_space

  USE elastic3d

  IMPLICIT NONE

#include "../include.f90"

  PRIVATE :: okada

CONTAINS

  !----------------------------------------------------------------------------
  !> Subroutine Okada
  !! converts Aki's to Okada's fault representation and evaluates the
  !! displacement at position (xrec,yrec,zrec) due to a series of faults 
  !! with given slip, lengths, widths, rake, etc...
  !!
  !! lambda, mu = the two Lame constants in Pascal (SI unit)
  !! (xs,ys,zs) = coordinates of the start point of strike
  !! with x = north, y = east, z = downward.
  !! all angles in degree.
  !! (xrec,yrec,zrec) = cartesian coordinates of observations
  !! disp = 3 displacement components: ux,uy,uz
  !!
  !!   \author R. Wang              - Created, Potsdam, Nov, 2001
  !!   \author S. Barbot (05/07/10) - Updated to Fortran90
  !----------------------------------------------------------------------------
  SUBROUTINE okada(lambda,mu,slip,opening,x1,x2,x3, &
                   width,length,strike,dipd,rake, &
                   xrec,yrec,zrec,disp)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: lambda,mu,slip,opening
    REAL*8, INTENT(IN) :: x1,x2,x3,width,length,strike,dipd,rake
    REAL*8, INTENT(IN) :: xrec,yrec,zrec
    REAL*8, INTENT(OUT) :: disp(3)

    ! from Okada's subroutine DC3D0:
    INTEGER IRET
    REAL*4 ALPHA,X,Y,Z,DEPTH,DIP,POT3,POT4,& 
           UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ

    ! more from Okada's subroutine DC3D:
    REAL*4 AL1,AL2,AW1,AW2,DISL1,DISL2,DISL3

    ! LOCAL CONSTANTS
    REAL*8 degtorad,eps
    PARAMETER(degtorad=1.745329252E-02,eps=1.0d-06)

    REAL*8 st,di,ra
    REAL*8 csst,ssst,csra,ssra,csdi,ssdi

    ! receiver and source independent variables
    ALPHA=sngl((lambda+mu)/(lambda+2.d0*mu))
    POT3=0.0
    POT4=0.0
    DISL3=0.0
    AL1=0.0
    AW2=0.0
    Z=REAL(-zrec)

    ! initialization
    disp(:)=0.d0

    st=strike*degtorad
    csst=dcos(st)
    ssst=dsin(st)

    di=dipd*degtorad
    csdi=dcos(di)
    ssdi=dsin(di)

    ra=rake*degtorad
    csra=dcos(ra)
    ssra=dsin(ra)
 
       ! transform from Aki's to Okada's system
    X=sngl((xrec-x1)*csst+(yrec-x2)*ssst)
    Y=sngl((xrec-x1)*ssst-(yrec-x2)*csst)
    DEPTH=sngl(x3)
    DIP=sngl(dipd)

    ! finite source
    AL2=sngl(length)
    AW1=-sngl(width)
    DISL1=sngl(slip*csra)
    DISL2=sngl(slip*ssra)
    DISL3=sngl(opening)

    IRET=1
    CALL DC3D(ALPHA,X,Y,Z,DEPTH,DIP,AL1,AL2,AW1,AW2, &
              DISL1,DISL2,DISL3,UX,UY,UZ, &
              UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ,IRET)

    IF (IRET .EQ. 0) THEN
       ! transform from Okada's to Aki's system
       disp(1)=disp(1)+dble(UX)*csst+dble(UY)*ssst
       disp(2)=disp(2)+dble(UX)*ssst-dble(UY)*csst
       disp(3)=disp(3)-dble(UZ)
    ELSE
       ! singular point at fault edge
       disp=0._8
    ENDIF

  END SUBROUTINE okada

  !--------------------------------------------------------------------
  !> subroutine dislocations_okada
  !! evaluate the deformation due to the motion of dislocation (shear
  !! and opening)
  !!
  !! \author sylvain barbot (05/07/10) - original form 
  !--------------------------------------------------------------------
  SUBROUTINE grnfct_okada(lambda,mu,slip,opening,x,y,z,width,length, &
                          strike,dip,rake,sx1,sx2,sx3,dx1,dx2,dx3,u1,u2,u3)
    INTEGER, INTENT(IN) :: sx1,sx2,sx3
#ifdef ALIGN_DATA
    REAL*4, DIMENSION(sx1+2,sx2,sx3), INTENT(OUT) :: u1,u2,u3
#else
    REAL*4, DIMENSION(sx1,sx2,sx3), INTENT(OUT) :: u1,u2,u3
#endif
    REAL*8, INTENT(IN) :: lambda,mu,dx1,dx2,dx3
    REAL*8, INTENT(IN) :: slip,opening,x,y,z,width,length,strike,dip,rake

    INTEGER :: i1,i2,i3
    REAL*8 :: dum,x1,x2,x3
    REAL*8, DIMENSION(3) :: disp

    DO i3=1,sx3
       x3=DBLE(i3-1)*dx3
       DO i2=1,sx2
          DO i1=1,sx1
             CALL shiftedcoordinates(i1,i2,i3,sx1,sx2,sx3, &
                                     dx1,dx2,dx3,x1,x2,dum)
             CALL okada(lambda,mu,slip,opening,x,y,z, &
                        width,length,strike,dip,rake,x1,x2,x3,disp)
             u1(i1,i2,i3)=u1(i1,i2,i3)+REAL(disp(1),4)
             u2(i1,i2,i3)=u2(i1,i2,i3)+REAL(disp(2),4)
             u3(i1,i2,i3)=u3(i1,i2,i3)+REAL(disp(3),4)
          END DO
       END DO
    END DO

  END SUBROUTINE grnfct_okada

END MODULE green_space
