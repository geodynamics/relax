!-----------------------------------------------------------------------
! Copyright 2016 Sylvain Barbot
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
!\author Sagar Masuti
!-----------------------------------------------------------------------

#include "include.f90"

MODULE util
    
   USE elastic3d
   IMPLICIT NONE

CONTAINS 

  SUBROUTINE ispresent(var,avail)
    REAL*4, OPTIONAL :: var 
    INTEGER :: avail

    IF (PRESENT(var)) avail=1

  END SUBROUTINE
  
  INTEGER FUNCTION isallocated(var)
    TYPE(LAYER_STRUCT), DIMENSION(:), ALLOCATABLE :: var 

    isallocated=0
    IF (ALLOCATED(var)) isallocated=1

  END FUNCTION isallocated

  !------------------------------------------------------------------
  !> subroutine exportpoints
  !! sample a vector field at a series of points for export.
  !! each location is attributed a file in which the time evolution
  !! of the vector value is listed in the format:
  !!
  !!                t_0 u(t_0) v(t_0) w(t_0)
  !!                t_1 u(t_1) v(t_1) w(t_1)
  !!                ...
  !!
  !! \author sylvain barbot (11/10/07) - original form
  !------------------------------------------------------------------
  SUBROUTINE exportpoints(c1,c2,c3,sig,sx1,sx2,sx3,dx1,dx2,dx3, &
       opts,ptsname,time,wdir,isnew,x0,y0,rot)
    INTEGER, INTENT(IN) :: sx1,sx2,sx3
#ifdef ALIGN_DATA
    REAL*4, INTENT(IN), DIMENSION(sx1+2,sx2,sx3) :: c1,c2,c3
#else
    REAL*4, INTENT(IN), DIMENSION(sx1,sx2,sx3) :: c1,c2,c3
#endif
    TYPE(TENSOR), INTENT(IN), DIMENSION(sx1,sx2,sx3) :: sig
    TYPE(VECTOR_STRUCT), INTENT(IN), DIMENSION(:) :: opts
    CHARACTER(LEN=4), INTENT(IN), DIMENSION(:) :: ptsname
    REAL*8, INTENT(IN) :: dx1,dx2,dx3,time,x0,y0,rot
    CHARACTER(256), INTENT(IN) :: wdir
    LOGICAL, INTENT(IN) :: isnew

    INTEGER :: i1,i2,i3,n,k
    REAL*8 :: u1,u2,u3,v1,v2,v3,x1,x2,x3,y1,y2,y3
    TYPE(TENSOR) :: lsig
    INTEGER :: i,iostatus
    CHARACTER(256) :: file1,file2

    i=INDEX(wdir," ")
    n=SIZE(ptsname)

    DO k=1,n
       file1=wdir(1:i-1) // "/" // ptsname(k) // ".ned"
       IF (isnew) THEN
          OPEN (UNIT=15,FILE=file1,IOSTAT=iostatus,FORM="FORMATTED")
          WRITE (15,'("#          t          u1          u2          u3         ", &
                    & "s11         s12         s13         s22         s23         s33")')
       ELSE
          OPEN (UNIT=15,FILE=file1,POSITION="APPEND",&
               IOSTAT=iostatus,FORM="FORMATTED")
       END IF
       IF (iostatus>0) STOP "could not open point file for writing"

       x1=opts(k)%v1
       x2=opts(k)%v2
       x3=opts(k)%v3

       CALL shiftedindex(x1,x2,x3,sx1,sx2,sx3,dx1,dx2,dx3,i1,i2,i3)
#ifdef USING_CUDA
       CALL cuexportpoints (u1, u2, u3, lsig, %VAL(i1-1), %VAL(i2-1), %VAL(i3-1))
#else
       u1=c1(i1,i2,i3)
       u2=c2(i1,i2,i3)
       u3=c3(i1,i2,i3)
       lsig=sig(i1,i2,i3)
#endif
       ! change from computational reference frame to user reference system
       y1=x1;v1=u1
       y2=x2;v2=u2
       y3=x3;v3=u3

       CALL rotation(y1,y2,-rot)
       y1=y1+x0
       y2=y2+y0
       CALL rotation(v1,v2,-rot)

       x1=x1+x0
       x2=x2+y0

       WRITE (15,'(14ES12.4E2)') time,v1,v2,v3, &
                                 lsig%s11,lsig%s12,lsig%s13, &
                                 lsig%s22,lsig%s23,lsig%s33
       CLOSE(15)

       ! output in rotated coordinates system
       IF (0._8.NE.rot) THEN
          file2=wdir(1:i-1) // "/" // ptsname(k) // ".c.txt"
          IF (isnew) THEN
             OPEN (UNIT=16,FILE=file2,IOSTAT=iostatus,FORM="FORMATTED")
             WRITE (16,'("#         t         u1         u2         u3")')
          ELSE
             OPEN (UNIT=16,FILE=file2,POSITION="APPEND",&
                   IOSTAT=iostatus,FORM="FORMATTED")
          END IF
          IF (iostatus>0) STOP "could not open point file for writing"

          WRITE (16,'(7ES11.3E2)') time,u1,u2,u3
          CLOSE(16)
       END IF

    END DO

  CONTAINS

    !------------------------------------------------------------------
    ! subroutine Rotation
    ! rotates a point coordinate into the computational reference
    ! system.
    ! 
    ! sylvain barbot (04/16/09) - original form
    !------------------------------------------------------------------
    SUBROUTINE rotation(x,y,rot)
      REAL*8, INTENT(INOUT) :: x,y
      REAL*8, INTENT(IN) :: rot

      REAL*8 :: alpha,xx,yy
      REAL*8, PARAMETER :: DEG2RAD = 0.01745329251994329547437168059786927_8


      alpha=rot*DEG2RAD
      xx=x
      yy=y

      x=+xx*cos(alpha)+yy*sin(alpha)
      y=-xx*sin(alpha)+yy*cos(alpha)

    END SUBROUTINE rotation

  END SUBROUTINE exportpoints
  !------------------------------------------------------------------
  !> subroutine pts2series
  !! sample a vector field at a series of points for export.
  !! each location is attributed a file in which the time evolution
  !! of the vector value is listed in the format:
  !!
  !!                t_0 u(t_0) v(t_0) w(t_0)
  !!                t_1 u(t_1) v(t_1) w(t_1)
  !!                ...
  !!
  !! \author sylvain barbot (11/10/07) - original form
  !------------------------------------------------------------------
  SUBROUTINE pts2series(c1,c2,c3,sx1,sx2,sx3,dx1,dx2,dx3, &
       opts,time,index,gps)
    INTEGER, INTENT(IN) :: sx1,sx2,sx3
#ifdef ALIGN_DATA
    REAL*4, INTENT(IN), DIMENSION(sx1+2,sx2,sx3) :: c1,c2,c3
#else
    REAL*4, INTENT(IN), DIMENSION(sx1,sx2,sx3) :: c1,c2,c3
#endif
    TYPE(VECTOR_STRUCT), INTENT(IN), DIMENSION(:) :: opts
    REAL*8, INTENT(IN) :: dx1,dx2,dx3,time
    INTEGER, INTENT(IN) :: index
    TYPE(MANIFOLD_STRUCT), INTENT(INOUT), DIMENSION(:) :: gps

    INTEGER :: i1,i2,i3,k
    REAL*8 :: u1,u2,u3,x1,x2,x3
#ifdef USING_CUDA
    TYPE(TENSOR) :: lsig
#endif

    DO k=1,SIZE(opts)
       x1=opts(k)%v1
       x2=opts(k)%v2
       x3=opts(k)%v3

       CALL shiftedindex(x1,x2,x3,sx1,sx2,sx3,dx1,dx2,dx3,i1,i2,i3)

#ifdef USING_CUDA
       CALL cuexportpoints (u1, u2, u3, lsig, %VAL(i1-1), %VAL(i2-1), %VAL(i3-1))
#else
       u1=c1(i1,i2,i3)
       u2=c2(i1,i2,i3)
       u3=c3(i1,i2,i3)
#endif
       gps(k)%nepochs=index
       gps(k)%t(index)=time
       gps(k)%u1(index)=u1
       gps(k)%u2(index)=u2
       gps(k)%u3(index)=u3

   END DO
  END SUBROUTINE pts2series  
END MODULE

