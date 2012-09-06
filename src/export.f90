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

MODULE export

  USE elastic3d
  USE viscoelastic3d
  USE friction3d

  IMPLICIT NONE

  PRIVATE xyzwrite
  PRIVATE geoxyzwrite

CONTAINS

  !-------------------------------------------------------------------
  ! routine ReportTime
  ! writes the times of exports
  !
  ! sylvain barbot (04/29/09) - original form
  !-------------------------------------------------------------------
  SUBROUTINE reporttime(i,t,repfile)
    INTEGER, INTENT(IN) :: i
    CHARACTER(256), INTENT(IN) :: repfile
    REAL*8, INTENT(IN) :: t

    INTEGER :: iostatus

    IF (0 .eq. i) THEN
       OPEN (UNIT=15,FILE=repfile,IOSTAT=iostatus,FORM="FORMATTED")
    ELSE
       OPEN (UNIT=15,FILE=repfile,POSITION="APPEND",&
            IOSTAT=iostatus,FORM="FORMATTED")
    END IF
    IF (iostatus>0) THEN
       WRITE_DEBUG_INFO
       PRINT '(a)', repfile
       STOP "could not open file for export"
    END IF

    WRITE (15,'(ES11.3E2)') t

    CLOSE(15)

  END SUBROUTINE reporttime

  SUBROUTINE report(i,t,file1,file2,file3,sx1,sx2,repfile)
    INTEGER, INTENT(IN) :: i,sx1,sx2
    CHARACTER(256), INTENT(IN) :: file1,file2,file3,repfile
    REAL*8, INTENT(IN) :: t

    INTEGER :: iostatus, ind1,ind2,ind3

    IF (0 .eq. i) THEN
       OPEN (UNIT=15,FILE=repfile,IOSTAT=iostatus,FORM="FORMATTED")
    ELSE
       OPEN (UNIT=15,FILE=repfile,POSITION="APPEND",&
            IOSTAT=iostatus,FORM="FORMATTED")
    END IF
    IF (iostatus>0) THEN
       WRITE_DEBUG_INFO
       PRINT '(a)', repfile
       STOP "could not open file for export"
    END IF

    ind1=INDEX(file1," ")
    ind2=INDEX(file2," ")
    ind3=INDEX(file3," ")
    WRITE (15,'(I3.3,2I6," ",f13.4," ",a," ",a," ",a)') i,sx1,sx2,t,&
         file1(1:ind1-1),file2(1:ind2-1),file3(1:ind3-1)

    CLOSE(15)

  END SUBROUTINE report

  SUBROUTINE export2d(data,sx1,sx2,filename)
    INTEGER, INTENT(IN) :: sx1,sx2
    REAL*4, INTENT(IN), DIMENSION(sx1,sx2) :: data
    CHARACTER(256), INTENT(IN) :: filename

    INTEGER :: iostatus,i1,i2
    CHARACTER(15) :: form
    CHARACTER(5) :: digit

    WRITE (digit,'(I5.5)') sx1
    form="("//digit//"ES11.3E2)"

    OPEN (UNIT=15,FILE=filename,IOSTAT=iostatus,FORM="FORMATTED")
    IF (iostatus>0) THEN
       WRITE_DEBUG_INFO
       STOP "could not open file for export"
    END IF

    WRITE (15,form) ((data(i1,i2), i1=1,sx1), i2=1,sx2)
    CLOSE(15)

  END SUBROUTINE export2d

  !------------------------------------------------------------------
  ! subroutine geoxyzwrite
  !
  ! sylvain barbot (22/05/10) - original form
  !------------------------------------------------------------------
  SUBROUTINE geoxyzwrite(x,y,z,sx1,sx2,filename)
    INTEGER, INTENT(IN) :: sx1,sx2
    REAL*4, INTENT(IN), DIMENSION(sx1,sx2) :: z
    REAL*8, INTENT(IN), DIMENSION(sx1,sx2) :: x,y
    CHARACTER(256), INTENT(IN) :: filename

    INTEGER :: iostatus,i1,i2

    OPEN (UNIT=15,FILE=filename,IOSTAT=iostatus,FORM="FORMATTED")
    IF (iostatus>0) STOP "could not open file for proj export"

    DO i2=1,sx2
       DO i1=1,sx1
          WRITE (15,'(ES15.8E1,ES15.8E1,ES11.3E2)'), &
                 x(i1,i2),y(i1,i2),z(i1,i2)
       END DO
    END DO
    CLOSE(15)

  END SUBROUTINE geoxyzwrite

  !------------------------------------------------------------------
  ! subroutine xyzwrite
  !
  ! sylvain barbot (06/10/09) - original form
  !------------------------------------------------------------------
  SUBROUTINE xyzwrite(data,sx1,sx2,dx1,dx2,filename)
    INTEGER, INTENT(IN) :: sx1,sx2
    REAL*4, INTENT(IN), DIMENSION(sx1,sx2) :: data
    CHARACTER(256), INTENT(IN) :: filename
    REAL*8 :: dx1,dx2

    INTEGER :: iostatus,i1,i2

    OPEN (UNIT=15,FILE=filename,IOSTAT=iostatus,FORM="FORMATTED")
    IF (iostatus>0) STOP "could not open file for export"

    DO i2=1,sx2
       DO i1=1,sx1
          !x1=(mod(sx1/2+i1-1,sx1)-sx1/2)*dx1
          !x2=(mod(sx2/2+i2-1,sx2)-sx2/2)*dx2
          WRITE (15,'(ES11.3E2,ES11.3E2,ES11.3E2)'), &
                DBLE(i2-1-sx2/2)*dx2,DBLE(i1-1-sx1/2)*dx1,data(i1,i2)
       END DO
    END DO
    CLOSE(15)

  END SUBROUTINE xyzwrite

#ifdef PROJ
  !------------------------------------------------------------------
  !> subroutine ExportStressPROJ
  !! export a map view of stress with coordinates in 
  !! longitude/latitude. Text format output is the GMT-compatible
  !! .xyz file format where data in each file is organized as follows
  !!
  !! longitude latitude s11 
  !! longitude latitude s12
  !! longitude latitude s13
  !! longitude latitude s22
  !! longitude latitude s23
  !! longitude latitude s33
  !!
  !! this is an interface to exportproj.
  !!
  !! \author sylvain barbot (05/22/10) - original form
  !------------------------------------------------------------------
  SUBROUTINE exportstressproj(sig,sx1,sx2,sx3,dx1,dx2,dx3,oz, &
                              x0,y0,lon0,lat0,zone,scale,wdir,index)
    INTEGER, INTENT(IN) :: index,sx1,sx2,sx3,zone
    TYPE(TENSOR), INTENT(IN), DIMENSION(sx1,sx2,sx3) :: sig
    REAL*8, INTENT(IN) :: oz,dx1,dx2,dx3,x0,y0,lon0,lat0,scale
    CHARACTER(256), INTENT(IN) :: wdir

    REAL*4, DIMENSION(:,:), ALLOCATABLE :: t1,t2,t3
    INTEGER :: iostatus,i,j,k,l

    ALLOCATE(t1(sx1+2,sx2),t2(sx1+2,sx2),t3(sx1+2,sx2),STAT=iostatus)
    IF (iostatus>0) STOP "could not allocate memory for grid export"

    k=fix(oz/dx3)+1
    DO j=1,sx2
       DO i=1,sx1
#ifdef ALIGN_DATA
          l=(j-1)*(sx1+2)+i
#else
          l=(j-1)*sx1+i
#endif
          t1(l,1)=sig(i,j,k)%s11
          t2(l,1)=sig(i,j,k)%s12
          t3(l,1)=sig(i,j,k)%s13
       END DO
    END DO

    CALL exportproj(t1,t2,t3,sx1,sx2,1,dx1,dx2,dx3,0._8, &
                  x0,y0,lon0,lat0,zone,scale,wdir,index,convention=4)

    DO j=1,sx2
       DO i=1,sx1
#ifdef ALIGN_DATA
          l=(j-1)*(sx1+2)+i
#else
          l=(j-1)*sx1+i
#endif
          t1(l,1)=sig(i,j,k)%s22
          t2(l,1)=sig(i,j,k)%s23
          t3(l,1)=sig(i,j,k)%s33
       END DO
    END DO

    CALL exportproj(t1,t2,t3,sx1,sx2,1,dx1,dx2,dx3,0._8, &
                  x0,y0,lon0,lat0,zone,scale,wdir,index,convention=5)

    DEALLOCATE(t1,t2,t3)

  END SUBROUTINE exportstressproj

  !------------------------------------------------------------------
  !> subroutine ExportPROJ
  !! export a map view of displacements with coordinates in 
  !! longitude/latitude. Text format output is the GMT-compatible
  !! .xyz file format where data in each file is organized as follows
  !!
  !! longitude latitude u1, 
  !! longitude latitude u2 and 
  !! longitude latitude -u3
  !!
  !! for index-geo-north.xyz, 
  !!     index-geo-east.xyz and 
  !!     index-geo-up.xyz, respectively.
  !!
  !! \author sylvain barbot (05/22/10) - original form
  !------------------------------------------------------------------
  SUBROUTINE exportproj(c1,c2,c3,sx1,sx2,sx3,dx1,dx2,dx3,oz, &
                        x0,y0,lon0,lat0,zone,scale,wdir,i,convention)
    INTEGER, INTENT(IN) :: i,sx1,sx2,sx3,zone
#ifdef ALIGN_DATA
    REAL*4, INTENT(IN), DIMENSION(sx1+2,sx2,sx3) :: c1,c2,c3
#else
    REAL*4, INTENT(IN), DIMENSION(sx1,sx2,sx3) :: c1,c2,c3
#endif
    REAL*8, INTENT(IN) :: oz,dx1,dx2,dx3,x0,y0,lon0,lat0,scale
    CHARACTER(256), INTENT(IN) :: wdir
    INTEGER, INTENT(IN), OPTIONAL :: convention

    INTEGER :: iostatus,i1,i2,pos,conv
    CHARACTER(3) :: digit
    REAL*4, DIMENSION(:,:), ALLOCATABLE :: t1,t2,t3
    REAL*8, DIMENSION(:,:), ALLOCATABLE :: x,y
    CHARACTER(256) :: file1,file2,file3
    REAL*8 :: lon1,lat1

    IF (PRESENT(convention)) THEN
       conv=convention
    ELSE
       conv=1
    END IF

    lon1=lon0
    lat1=lat0

    ALLOCATE(t1(sx1,sx2),t2(sx1,sx2),t3(sx1,sx2), &
             x(sx1,sx2),y(sx1,sx2),STAT=iostatus)
    IF (iostatus>0) STOP "could not allocate memory for export"

    CALL exportspatial(c1(:,:,int(oz/dx3)+1),sx1,sx2,t1)
    CALL exportspatial(c2(:,:,int(oz/dx3)+1),sx1,sx2,t2)
    CALL exportspatial(c3(:,:,int(oz/dx3)+1),sx1,sx2,t3)
    t3=-t3

    ! grid coordinates (x=easting, y=northing)
    DO i2=1,sx2
       DO i1=1,sx1
          y(i1,i2)=(i1-sx1/2)*(dx1*scale)+x0
          x(i1,i2)=(i2-sx2/2)*(dx2*scale)+y0
       END DO
    END DO
    CALL proj(x,y,sx1*sx2,lon1,lat1,zone)

    pos=INDEX(wdir," ")
    WRITE (digit,'(I3.3)') i
    SELECT CASE(conv)
    CASE (1) ! cumulative displacement
       file1=wdir(1:pos-1) // "/" // digit // "-geo-north.xyz"
       file2=wdir(1:pos-1) // "/" // digit // "-geo-east.xyz"
       file3=wdir(1:pos-1) // "/" // digit // "-geo-up.xyz"
    CASE (2) ! postseismic displacement
       file1=wdir(1:pos-1) // "/" // digit // "-relax-geo-north.xyz"
       file2=wdir(1:pos-1) // "/" // digit // "-relax-geo-east.xyz"
       file3=wdir(1:pos-1) // "/" // digit // "-relax-geo-up.xyz"
    CASE (3) ! equivalent body forces
       file1=wdir(1:pos-1) // "/" // digit // "-eqbf-geo-north.xyz"
       file2=wdir(1:pos-1) // "/" // digit // "-eqbf-geo-east.xyz"
       file3=wdir(1:pos-1) // "/" // digit // "-eqbf-geo-up.xyz"
    CASE (4) ! equivalent body forces
       file1=wdir(1:pos-1) // "/" // digit // "-geo-s11.xyz"
       file2=wdir(1:pos-1) // "/" // digit // "-geo-s12.xyz"
       file3=wdir(1:pos-1) // "/" // digit // "-geo-s13.xyz"
    CASE (5) ! equivalent body forces
       file1=wdir(1:pos-1) // "/" // digit // "-geo-s22.xyz"
       file2=wdir(1:pos-1) // "/" // digit // "-geo-s23.xyz"
       file3=wdir(1:pos-1) // "/" // digit // "-geo-s33.xyz"
    END SELECT
    
    CALL geoxyzwrite(x,y,t1,sx1,sx2,file1)
    CALL geoxyzwrite(x,y,t2,sx1,sx2,file2)
    CALL geoxyzwrite(x,y,t3,sx1,sx2,file3)

    DEALLOCATE(t1,t2,t3)

  END SUBROUTINE exportproj
#endif

#ifdef XYZ
  !------------------------------------------------------------------
  !> subroutine ExportXYZ
  !! export a map view of surface displacement into the GMT-compatible
  !! .xyz file format where data in each file is organized as follows
  !!
  !! x1 x2 u1, x1 x2 u2 and x1 x2 -u3
  !!
  !! for index-north.xyz, index-east.xyz and index-up.xyz, 
  !! respectively.
  !!
  !! \author sylvain barbot (06/10/09) - original form
  !------------------------------------------------------------------
  SUBROUTINE exportxyz(c1,c2,c3,sx1,sx2,sx3,oz,dx1,dx2,dx3,i,wdir)
    INTEGER, INTENT(IN) :: i,sx1,sx2,sx3
#ifdef ALIGN_DATA
    REAL*4, INTENT(IN), DIMENSION(sx1+2,sx2,sx3) :: c1,c2,c3
#else
    REAL*4, INTENT(IN), DIMENSION(sx1,sx2,sx3) :: c1,c2,c3
#endif
    REAL*8, INTENT(IN) :: oz,dx1,dx2,dx3
    CHARACTER(256), INTENT(IN) :: wdir

    INTEGER :: iostatus,pos
    REAL*4, DIMENSION(:,:), ALLOCATABLE :: temp1,temp2,temp3
    CHARACTER(256) :: file1,file2,file3
    CHARACTER(3) :: digit

    ALLOCATE(temp1(sx1,sx2),temp2(sx1,sx2),temp3(sx1,sx2),STAT=iostatus)
    IF (iostatus>0) STOP "could not allocate memory for export"

    CALL exportspatial(c1(:,:,int(oz/dx3)+1),sx1,sx2,temp1)
    CALL exportspatial(c2(:,:,int(oz/dx3)+1),sx1,sx2,temp2)
    CALL exportspatial(c3(:,:,int(oz/dx3)+1),sx1,sx2,temp3)
    temp3=-temp3

    pos=INDEX(wdir," ")
    WRITE (digit,'(I3.3)') i
    file1=wdir(1:pos-1) // "/" // digit // "-north.xyz"
    file2=wdir(1:pos-1) // "/" // digit // "-east.xyz"
    file3=wdir(1:pos-1) // "/" // digit // "-up.xyz"

    CALL xyzwrite(temp1,sx1,sx2,dx1,dx2,file1)
    CALL xyzwrite(temp2,sx1,sx2,dx1,dx2,file2)
    CALL xyzwrite(temp3,sx1,sx2,dx1,dx2,file3)

    DEALLOCATE(temp1,temp2,temp3)

  END SUBROUTINE exportxyz
#endif

#ifdef TXT
  !------------------------------------------------------------------
  ! subroutine ExportTXT
  ! exports a horizontal slice of uniform depth into specified text
  ! files and adds filenames in the report file.
  ! if i is set to 0, the report file is reinitiated.
  ! input data c1,c2,c3 are in the space domain.
  !------------------------------------------------------------------
  SUBROUTINE exporttxt(c1,c2,c3,sx1,sx2,sx3,oz,dx3,i,time,wdir,reportfilename)
    INTEGER, INTENT(IN) :: i,sx1,sx2,sx3
#ifdef ALIGN_DATA
    REAL*4, INTENT(IN), DIMENSION(sx1+2,sx2,sx3) :: c1,c2,c3
#else
    REAL*4, INTENT(IN), DIMENSION(sx1,sx2,sx3) :: c1,c2,c3
#endif
    REAL*8, INTENT(IN) :: oz,dx3,time
    CHARACTER(256), INTENT(IN) :: wdir,reportfilename

    INTEGER :: iostatus,pos
    REAL*4, DIMENSION(:,:), ALLOCATABLE :: temp1,temp2,temp3
    CHARACTER(3) :: digit
    CHARACTER(256) :: file1,file2,file3
    
    ALLOCATE(temp1(sx1,sx2),temp2(sx1,sx2),temp3(sx1,sx2),STAT=iostatus)
    IF (iostatus>0) STOP "could not allocate memory for export"

    CALL exportspatial(c1(:,:,int(oz/dx3)+1),sx1,sx2,temp1)
    CALL exportspatial(c2(:,:,int(oz/dx3)+1),sx1,sx2,temp2)
    CALL exportspatial(c3(:,:,int(oz/dx3)+1),sx1,sx2,temp3)

    pos=INDEX(wdir," ")
    WRITE (digit,'(I3.3)') i
    file1=wdir(1:pos-1) // "/" // digit // "-u1.txt"
    file2=wdir(1:pos-1) // "/" // digit // "-u2.txt"
    file3=wdir(1:pos-1) // "/" // digit // "-u3.txt"
    
    CALL export2d(temp1,sx1,sx2,file1)
    CALL export2d(temp2,sx1,sx2,file2)
    CALL export2d(temp3,sx1,sx2,file3)
    
    file1=digit // "-u1.txt "
    file2=digit // "-u2.txt "
    file3=digit // "-u3.txt "
    CALL report(i,time,file1,file2,file3,sx1,sx2,reportfilename)

    DEALLOCATE(temp1,temp2,temp3)

  END SUBROUTINE exporttxt
#endif

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
       file1=wdir(1:i-1) // "/" // ptsname(k) // ".txt"
       IF (isnew) THEN
          OPEN (UNIT=15,FILE=file1,IOSTAT=iostatus,FORM="FORMATTED")
          WRITE (15,'("#         t         u1         u2         u3        ", &
                    & "s11        s12        s13        s22        s23        s33")')
       ELSE
          OPEN (UNIT=15,FILE=file1,POSITION="APPEND",&
               IOSTAT=iostatus,FORM="FORMATTED")
       END IF
       IF (iostatus>0) STOP "could not open point file for writing"

       x1=opts(k)%v1
       x2=opts(k)%v2
       x3=opts(k)%v3

       CALL shiftedindex(x1,x2,x3,sx1,sx2,sx3,dx1,dx2,dx3,i1,i2,i3)

       u1=c1(i1,i2,i3)
       u2=c2(i1,i2,i3)
       u3=c3(i1,i2,i3)
       lsig=sig(i1,i2,i3)

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

       WRITE (15,'(13ES11.3E2)') time,v1,v2,v3, &
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

  !---------------------------------------------------------------------
  !> subroutine exportoptsdat
  !! export the coordinates and name of the observation points (often
  !! coordinates of GPS instruments or such) for display with GMT in the
  !! ASCII format. The file contains a list of x1,x2,x3 coordinates and
  !! a 4-character name string.
  !!
  !! input variables
  !! @param n          - number of observation points
  !! @param opts       - coordinates of observation points
  !! @param ptsname    - name of obs. points
  !! @param filename   - output file (example: wdir/opts.xy)
  !!
  !! \author sylvain barbot (08/10/11) - original form
  !---------------------------------------------------------------------
  SUBROUTINE exportoptsdat(n,opts,ptsname,filename)
    INTEGER, INTENT(IN) :: n
    TYPE(VECTOR_STRUCT), DIMENSION(n) :: opts
    CHARACTER(LEN=4), DIMENSION(n) :: ptsname
    CHARACTER(256) :: filename

    INTEGER :: k,iostatus

    IF (n.LE.0) RETURN

    OPEN (UNIT=15,FILE=filename,IOSTAT=iostatus,FORM="FORMATTED")
    IF (iostatus>0) STOP "could not open .xy file to export observation points"
    DO k=1,n
       WRITE (15,'(3ES11.4E1,X,a)') opts(k)%v1,opts(k)%v2,opts(k)%v3,ptsname(k)
    END DO
    CLOSE(15)
    
  END SUBROUTINE exportoptsdat
    
  !---------------------------------------------------------------------
  !> subroutine exportPlaneStress
  !! samples the value of an input tensor field at the location of 
  !! defined plane (position, strike, dip, length and width).
  !!
  !! input variables
  !! @param sig        - sampled tensor array
  !! @param nop        - number of observation planes
  !! @param op         - structure of observation planes (position, orientation)
  !! @param x0, y0 - origin position of coordinate system
  !! @param dx1,2,3    - sampling size
  !! @param sx1,2,3    - size of the scalar field
  !! @param wdir       - output directory for writing
  !! @param i          - loop index to suffix file names
  !!
  !! creates files 
  !!
  !!    wdir/index.s00001.estrain.txt with TXT_EXPORTEIGENSTRAIN option
  !!
  !!    wdir/index.s00001.estrain.grd with GRD_EXPORTEIGENSTRAIN option
  !! 
  !! \author sylvain barbot (01/01/07) - original form
  !                         (02/25/10) - output in TXT and GRD formats
  !---------------------------------------------------------------------
  SUBROUTINE exportplanestress(sig,nop,op,x0,y0,dx1,dx2,dx3,sx1,sx2,sx3,wdir,i)
    INTEGER, INTENT(IN) :: nop,sx1,sx2,sx3,i
    TYPE(PLANE_STRUCT), INTENT(IN), DIMENSION(nop) :: op
    TYPE(TENSOR), INTENT(IN), DIMENSION(sx1,sx2,sx3) :: sig
    REAL*8, INTENT(IN) :: x0,y0,dx1,dx2,dx3
    CHARACTER(256), INTENT(IN) :: wdir

    INTEGER :: k,ns1,ns2
    TYPE(SLIPPATCH_STRUCT), DIMENSION(:,:), ALLOCATABLE :: slippatch
    CHARACTER(3) :: sdigit
    CHARACTER(3) :: digit
#ifdef TXT_EXPORTEIGENSTRAIN
    INTEGER :: iostatus
#endif
!#_indef GRD_EXPORTEIGENSTRAIN
    CHARACTER(256) :: fn11,fn12,fn13,fn22,fn23,fn33
    INTEGER :: j,j1,j2
    REAL*4, ALLOCATABLE, DIMENSION(:,:) :: temp11,temp12,temp13, &
                                           temp22,temp23,temp33
    REAL*8 :: rland=9998.,rdum=9999.
    REAL*8 :: xmin,ymin
    CHARACTER(256) :: title="monitor tensor field "
!#_endif

    IF (nop .le. 0) RETURN

    WRITE (digit,'(I3.3)') i

    DO k=1,nop
       CALL monitorstressfield(op(k)%x,op(k)%y,op(k)%z, &
            op(k)%width,op(k)%length,op(k)%strike,op(k)%dip, &
            0._8,sx1,sx2,sx3,dx1,dx2,dx3,sig,slippatch)

       IF (.NOT. ALLOCATED(slippatch)) THEN
          WRITE_DEBUG_INFO
          WRITE (0,'("could not monitor slip")')
          STOP 2
       END IF

       ns1=SIZE(slippatch,1)
       ns2=SIZE(slippatch,2)
          
       slippatch(:,:)%x1=slippatch(:,:)%x1+x0
       slippatch(:,:)%x2=slippatch(:,:)%x2+y0

       WRITE (sdigit,'(I3.3)') k

!#_ifdef GRD_EXPORTEIGENSTRAIN
       fn11=trim(wdir)//"/"//digit//".op"//sdigit//"-s11.grd"
       fn12=trim(wdir)//"/"//digit//".op"//sdigit//"-s12.grd"
       fn13=trim(wdir)//"/"//digit//".op"//sdigit//"-s13.grd"
       fn22=trim(wdir)//"/"//digit//".op"//sdigit//"-s22.grd"
       fn23=trim(wdir)//"/"//digit//".op"//sdigit//"-s23.grd"
       fn33=trim(wdir)//"/"//digit//".op"//sdigit//"-s33.grd"

       ! convert to c standard
       j=INDEX(fn11," ")
       fn11(j:j)=char(0)
       fn12(j:j)=char(0)
       fn13(j:j)=char(0)
       fn22(j:j)=char(0)
       fn23(j:j)=char(0)
       fn33(j:j)=char(0)

       ALLOCATE(temp11(ns1,ns2),temp12(ns1,ns2),temp13(ns1,ns2), &
                temp22(ns1,ns2),temp23(ns1,ns2),temp33(ns1,ns2),STAT=iostatus)
       IF (iostatus>0) STOP "could not allocate temporary array for GRD slip export."

       DO j2=1,ns2
          DO j1=1,ns1
             temp11(ns1+1-j1,j2)=slippatch(j1,j2)%sig%s11
             temp12(ns1+1-j1,j2)=slippatch(j1,j2)%sig%s12
             temp13(ns1+1-j1,j2)=slippatch(j1,j2)%sig%s13
             temp22(ns1+1-j1,j2)=slippatch(j1,j2)%sig%s22
             temp23(ns1+1-j1,j2)=slippatch(j1,j2)%sig%s23
             temp33(ns1+1-j1,j2)=slippatch(j1,j2)%sig%s33
          END DO
       END DO

       ! xmin is the lowest coordinates (positive eastward in GMT)
       xmin= MINVAL(slippatch(:,:)%lx)
       ! ymin is the lowest coordinates (positive northward in GMT)
       ymin=-MAXVAL(slippatch(:,:)%lz)

       ! call the c function "writegrd_"
       CALL writegrd(temp11,ns1,ns2,ymin,xmin,dx3,dx2,rland,rdum,title,fn11)
       CALL writegrd(temp12,ns1,ns2,ymin,xmin,dx3,dx2,rland,rdum,title,fn12)
       CALL writegrd(temp13,ns1,ns2,ymin,xmin,dx3,dx2,rland,rdum,title,fn13)
       CALL writegrd(temp22,ns1,ns2,ymin,xmin,dx3,dx2,rland,rdum,title,fn22)
       CALL writegrd(temp23,ns1,ns2,ymin,xmin,dx3,dx2,rland,rdum,title,fn23)
       CALL writegrd(temp33,ns1,ns2,ymin,xmin,dx3,dx2,rland,rdum,title,fn33)

       DEALLOCATE(temp11,temp12,temp13,temp22,temp23,temp33)

!#_endif

       DEALLOCATE(slippatch)
    END DO

END SUBROUTINE exportplanestress

  !---------------------------------------------------------------------
  !> subroutine exportEigenstrain
  !! samples the value of an input scalar field at the location of 
  !! defined plane (position, strike, dip, length and width).
  !!
  !! input variables
  !! @param field      - sampled scalar array
  !! @param nop        - number of observation planes
  !! @param op         - structure of observation planes (position, orientation)
  !! @param x0, y0 - origin position of coordinate system
  !! @param dx1,2,3    - sampling size
  !! @param sx1,2,3    - size of the scalar field
  !! @param wdir       - output directory for writing
  !! @param i          - loop index to suffix file names
  !!
  !! creates files 
  !!
  !!    wdir/index.s00001.estrain.txt with TXT_EXPORTEIGENSTRAIN option
  !!
  !!    wdir/index.s00001.estrain.grd with GRD_EXPORTEIGENSTRAIN option
  !! 
  !! \author sylvain barbot (01/01/07) - original form
  !                         (02/25/10) - output in TXT and GRD formats
  !---------------------------------------------------------------------
  SUBROUTINE exporteigenstrain(field,nop,op,x0,y0,dx1,dx2,dx3,sx1,sx2,sx3,wdir,i)
    INTEGER, INTENT(IN) :: nop,sx1,sx2,sx3,i
    TYPE(PLANE_STRUCT), INTENT(IN), DIMENSION(nop) :: op
#ifdef ALIGN_DATA
    REAL*4, INTENT(IN), DIMENSION(sx1+2,sx2,sx3) :: field
#else
    REAL*4, INTENT(IN), DIMENSION(sx1,sx2,sx3) :: field
#endif
    REAL*8, INTENT(IN) :: x0,y0,dx1,dx2,dx3
    CHARACTER(256), INTENT(IN) :: wdir

    INTEGER :: k,ns1,ns2,pos
    TYPE(SLIPPATCH_STRUCT), DIMENSION(:,:), ALLOCATABLE :: slippatch
    CHARACTER(5) :: sdigit
    CHARACTER(3) :: digit
#ifdef TXT_EXPORTEIGENSTRAIN
    INTEGER :: iostatus,i1,i2
    CHARACTER(256) :: outfiletxt
#endif
!#_indef GRD_EXPORTEIGENSTRAIN
    CHARACTER(256) :: outfilegrd
    INTEGER :: j,j1,j2
    REAL*4, ALLOCATABLE, DIMENSION(:,:) :: temp
    REAL*8 :: rland=9998.,rdum=9999.
    REAL*8 :: xmin,ymin
    CHARACTER(256) :: title="monitor field "
!#_endif

    IF (nop .le. 0) RETURN

    pos=INDEX(wdir," ")
    WRITE (digit,'(I3.3)') i

    DO k=1,nop
       CALL monitorfield(op(k)%x,op(k)%y,op(k)%z, &
            op(k)%width,op(k)%length,op(k)%strike,op(k)%dip, &
            0._8,sx1,sx2,sx3,dx1,dx2,dx3,field,slippatch)

       IF (.NOT. ALLOCATED(slippatch)) THEN
          WRITE_DEBUG_INFO
          WRITE (0,'("could not monitor slip")')
          STOP 2
       END IF

       ns1=SIZE(slippatch,1)
       ns2=SIZE(slippatch,2)
          
       slippatch(:,:)%x1=slippatch(:,:)%x1+x0
       slippatch(:,:)%x2=slippatch(:,:)%x2+y0

       WRITE (sdigit,'(I5.5)') k
#ifdef TXT_EXPORTEIGENSTRAIN
       outfiletxt=wdir(1:pos-1)//"/"//digit//".s"//sdigit//".estrain.txt"
       
       OPEN (UNIT=15,FILE=outfiletxt,IOSTAT=iostatus,FORM="FORMATTED")
       IF (iostatus>0) STOP "could not open file for export"
          
       WRITE (15,'(6ES11.3E2)') ((slippatch(i1,i2), i1=1,ns1), i2=1,ns2)
          
       CLOSE(15)
#endif

!#_ifdef GRD_EXPORTEIGENSTRAIN
       outfilegrd=wdir(1:pos-1)//"/"//digit//".s"//sdigit//".estrain.grd"

       ! convert to c standard
       j=INDEX(outfilegrd," ")
       outfilegrd(j:j)=char(0)

       ALLOCATE(temp(ns1,ns2),STAT=iostatus)
       IF (iostatus>0) STOP "could not allocate temporary array for GRD slip export."

       DO j2=1,ns2
          DO j1=1,ns1
             temp(ns1+1-j1,j2)=REAL(slippatch(j1,j2)%slip)
          END DO
       END DO

       ! xmin is the lowest coordinates (positive eastward in GMT)
       xmin= MINVAL(slippatch(:,:)%lx)
       ! ymin is the lowest coordinates (positive northward in GMT)
       ymin=-MAXVAL(slippatch(:,:)%lz)

       ! call the c function "writegrd_"
       CALL writegrd(temp,ns1,ns2,ymin,xmin,dx3,dx2, &
                     rland,rdum,title,outfilegrd)

       DEALLOCATE(temp)

!#_endif

       DEALLOCATE(slippatch)
    END DO

END SUBROUTINE exporteigenstrain

  !---------------------------------------------------------------------
  !> subroutine exportCreep_asc
  !! evaluates the value of creep velocity at the location of 
  !! defined plane (position, strike, dip, length and width).
  !!
  !! input variables
  !! @param np         - number of frictional planes
  !! @param n          - array of frictional planes (position, orientation)
  !! @param structure  - array of depth-dependent frictional properties
  !! @param x0, y0     - origin position of coordinate system
  !! @param dx1,2,3    - sampling size
  !! @param sx1,2,3    - size of the stress tensor field
  !! @param beta       - smoothing factor controlling the extent of planes
  !! @param wdir       - output directory for writing
  !! @param i          - loop index to suffix file names
  !!
  !! creates files 
  !!
  !!    wdir/index.s00001.creep.dat
  !!
  !! containing
  !!
  !!    x1,x2,x3,x',y',slip,ss,ds,vel,ss vel,ds vel,taus,sigij
  !!
  !! where
  !! 
  !!   x1, x2, x3          are the absolute coordinates of the fault samples
  !!   x', y'              are the local coordinates (along-strike and down-dip)
  !!   slip, ss, ds        are the total, strike- and dip- slip components
  !!   vel, ss vel, ds vel are the total, strike- and dip- velocity components
  !!   taus, sigij         are the shear stress and the 6 independent components
  !!                       of the stress tensor
  !!
  !! with TXT_EXPORTCREEP option and
  !!
  !!    wdir/index.s00001.slip-strike.grd 
  !!    wdir/index.s00001.slip-dip.grd 
  !!    wdir/index.s00001.slip.grd 
  !!
  !! with GRD_EXPORTCREEP option where the suffix -strike stands for
  !! strike-slip, -dip for dip slip.
  !!
  !! file wdir/index.s00001.creep.dat is subsampled by a factor "skip"
  !! compared to the grd files.
  !! 
  !! \author sylvain barbot (01/01/07) - original form
  !!                        (02/25/10) - output in TXT and GRD formats
  !---------------------------------------------------------------------
#define TXT_EXPORTCREEP
  SUBROUTINE exportcreep_asc(np,n,beta,sig,structure, &
                         sx1,sx2,sx3,dx1,dx2,dx3,x0,y0,wdir,i)
    INTEGER, INTENT(IN) :: np,sx1,sx2,sx3,i
    TYPE(PLANE_STRUCT), INTENT(INOUT), DIMENSION(np) :: n
    TYPE(TENSOR), INTENT(IN), DIMENSION(sx1,sx2,sx3) :: sig
    TYPE(LAYER_STRUCT), DIMENSION(:), INTENT(IN) :: structure
    REAL*8, INTENT(IN) :: x0,y0,dx1,dx2,dx3,beta
    CHARACTER(256), INTENT(IN) :: wdir

    INTEGER :: k,ns1,ns2,pos
    CHARACTER(5) :: sdigit
    CHARACTER(3) :: digit
    CHARACTER(256) :: outfile
    INTEGER :: skip=3
    INTEGER :: iostatus,i1,i2

    IF (np .le. 0) RETURN

    pos=INDEX(wdir," ")
    WRITE (digit,'(I3.3)') i

    DO k=1,np

       ns1=SIZE(n(k)%patch,1)
       ns2=SIZE(n(k)%patch,2)
          
       WRITE (sdigit,'(I5.5)') k
       outfile=wdir(1:pos-1)//"/"//digit//".s"//sdigit//".creep.dat"
       
       OPEN (UNIT=15,FILE=outfile,IOSTAT=iostatus,FORM="FORMATTED")
       IF (iostatus>0) STOP "could not open file for export"
          
       WRITE (15,'("#        x1         x2         x3          yr        yz", &
                 & "       slip strike-slip  dip-slip   velocity     ss vel", &
                 & "     ds vel       taus      sig11      sig12      sig13      sig22      sig23      sig33")')
       WRITE (15,'(18ES11.3E2)') ((n(k)%patch(i1,i2)%x1,n(k)%patch(i1,i2)%x3,n(k)%patch(i1,i2)%x3, &
                                  n(k)%patch(i1,i2)%lx,n(k)%patch(i1,i2)%lz, &
                                  n(k)%patch(i1,i2)%slip, &
                                  n(k)%patch(i1,i2)%ss, &
                                  n(k)%patch(i1,i2)%ds, &
                                  n(k)%patch(i1,i2)%v, &
                                  n(k)%patch(i1,i2)%vss, &
                                  n(k)%patch(i1,i2)%vds, &
                                  n(k)%patch(i1,i2)%taus, &
                                  n(k)%patch(i1,i2)%sig%s11, &
                                  n(k)%patch(i1,i2)%sig%s12, &
                                  n(k)%patch(i1,i2)%sig%s13, &
                                  n(k)%patch(i1,i2)%sig%s22, &
                                  n(k)%patch(i1,i2)%sig%s23, &
                                  n(k)%patch(i1,i2)%sig%s33, &
                                  i1=1,ns1,skip), i2=1,ns2,skip)
       CLOSE(15)

    END DO

END SUBROUTINE exportcreep_asc

  !---------------------------------------------------------------------
  !> subroutine exportCreep_grd
  !! evaluates the value of creep velocity at the location of 
  !! defined plane (position, strike, dip, length and width).
  !!
  !! input variables
  !! @param np         - number of frictional planes
  !! @param n          - array of frictional planes (position, orientation)
  !! @param structure  - array of depth-dependent frictional properties
  !! @param x0, y0     - origin position of coordinate system
  !! @param dx1,2,3    - sampling size
  !! @param sx1,2,3    - size of the stress tensor field
  !! @param beta       - smoothing factor controlling the extent of planes
  !! @param wdir       - output directory for writing
  !! @param i          - loop index to suffix file names
  !!
  !! creates files 
  !!
  !!    wdir/index.s00001.creep.dat
  !!
  !! containing
  !!
  !!    x1,x2,x3,x',y',slip,ss,ds,vel,ss vel,ds vel,taus,sigij
  !!
  !! where
  !! 
  !!   x1, x2, x3          are the absolute coordinates of the fault samples
  !!   x', y'              are the local coordinates (along-strike and down-dip)
  !!   slip, ss, ds        are the total, strike- and dip- slip components
  !!   vel, ss vel, ds vel are the total, strike- and dip- velocity components
  !!   taus, sigij         are the shear stress and the 6 independent components
  !!                       of the stress tensor
  !!
  !! with TXT_EXPORTCREEP option and
  !!
  !!    wdir/index.s00001.slip-strike.grd 
  !!    wdir/index.s00001.slip-dip.grd 
  !!    wdir/index.s00001.slip.grd 
  !!
  !! with GRD_EXPORTCREEP option where the suffix -strike stands for
  !! strike-slip, -dip for dip slip.
  !!
  !! file wdir/index.s00001.creep.dat is subsampled by a factor "skip"
  !! compared to the grd files.
  !! 
  !! \author sylvain barbot (01/01/07) - original form
  !!                        (02/25/10) - output in TXT and GRD formats
  !---------------------------------------------------------------------
#define TXT_EXPORTCREEP
  SUBROUTINE exportcreep_grd(np,n,beta,sig,structure, &
                         sx1,sx2,sx3,dx1,dx2,dx3,x0,y0,wdir,i)
    INTEGER, INTENT(IN) :: np,sx1,sx2,sx3,i
    TYPE(PLANE_STRUCT), INTENT(INOUT), DIMENSION(np) :: n
    TYPE(TENSOR), INTENT(IN), DIMENSION(sx1,sx2,sx3) :: sig
    TYPE(LAYER_STRUCT), DIMENSION(:), INTENT(IN) :: structure
    REAL*8, INTENT(IN) :: x0,y0,dx1,dx2,dx3,beta
    CHARACTER(256), INTENT(IN) :: wdir

    INTEGER :: k,ns1,ns2,pos
    CHARACTER(5) :: sdigit
    CHARACTER(3) :: digit
    INTEGER :: j,iostatus,i1,i2
    REAL*4, ALLOCATABLE, DIMENSION(:,:) :: temp1,temp2,temp3
    REAL*8 :: rland=9998.,rdum=9999.
    REAL*8 :: xmin,ymin
    CHARACTER(256) :: title="monitor field "
    CHARACTER(256) :: file1,file2,file3

    IF (np .le. 0) RETURN

    pos=INDEX(wdir," ")
    WRITE (digit,'(I3.3)') i

    DO k=1,np

       ns1=SIZE(n(k)%patch,1)
       ns2=SIZE(n(k)%patch,2)
          
       WRITE (sdigit,'(I5.5)') k

       ALLOCATE(temp1(ns1,ns2),temp2(ns1,ns2),temp3(ns1,ns2),STAT=iostatus)
       IF (iostatus>0) STOP "could not allocate temporary arrays for GRD slip export."

       file1=wdir(1:pos-1)//"/"//digit//".s"//sdigit//".slip-dip.grd"
       file2=wdir(1:pos-1)//"/"//digit//".s"//sdigit//".slip-strike.grd"
       file3=wdir(1:pos-1)//"/"//digit//".s"//sdigit//".slip.grd"

       ! convert to c standard
       j=INDEX(file1," ")
       file1(j:j)=char(0)
       j=INDEX(file2," ")
       file2(j:j)=char(0)
       j=INDEX(file3," ")
       file3(j:j)=char(0)

       DO i2=1,ns2
          DO i1=1,ns1
             temp1(ns1+1-i1,i2)=REAL(n(k)%patch(i1,i2)%ds)
             temp2(ns1+1-i1,i2)=REAL(n(k)%patch(i1,i2)%ss)
             temp3(ns1+1-i1,i2)=REAL(n(k)%patch(i1,i2)%slip)
          END DO
       END DO

       ! xmin is the lowest coordinates (positive eastward in GMT)
       xmin= MINVAL(n(k)%patch(:,:)%lx)
       ! ymin is the lowest coordinates (positive northward in GMT)
       ymin=-MAXVAL(n(k)%patch(:,:)%lz)

       ! call the c function "writegrd_"
       CALL writegrd(temp1,ns1,ns2,ymin,xmin,dx3,dx2, &
                     rland,rdum,title,file1)
       CALL writegrd(temp2,ns1,ns2,ymin,xmin,dx3,dx2, &
                     rland,rdum,title,file2)
       CALL writegrd(temp3,ns1,ns2,ymin,xmin,dx3,dx2, &
                     rland,rdum,title,file3)

       file1=wdir(1:pos-1)//"/"//digit//".s"//sdigit//".vel-strike.grd"
       file2=wdir(1:pos-1)//"/"//digit//".s"//sdigit//".vel-dip.grd"
       file3=wdir(1:pos-1)//"/"//digit//".s"//sdigit//".vel.grd"

       ! convert to c standard
       j=INDEX(file1," ")
       file1(j:j)=char(0)
       j=INDEX(file2," ")
       file2(j:j)=char(0)
       j=INDEX(file3," ")
       file3(j:j)=char(0)

       DO i2=1,ns2
          DO i1=1,ns1
             temp1(ns1+1-i1,i2)=REAL(n(k)%patch(i1,i2)%vds)
             temp2(ns1+1-i1,i2)=REAL(n(k)%patch(i1,i2)%vss)
             temp3(ns1+1-i1,i2)=REAL(n(k)%patch(i1,i2)%v)
          END DO
       END DO

       ! xmin is the lowest coordinates (positive eastward in GMT)
       xmin= MINVAL(n(k)%patch(:,:)%lx)
       ! ymin is the lowest coordinates (positive northward in GMT)
       ymin=-MAXVAL(n(k)%patch(:,:)%lz)

       ! call the c function "writegrd_"
       CALL writegrd(temp1,ns1,ns2,ymin,xmin,dx3,dx2, &
                     rland,rdum,title,file1)
       CALL writegrd(temp2,ns1,ns2,ymin,xmin,dx3,dx2, &
                     rland,rdum,title,file2)
       CALL writegrd(temp3,ns1,ns2,ymin,xmin,dx3,dx2, &
                     rland,rdum,title,file3)

       DEALLOCATE(temp1,temp2,temp3)

    END DO

END SUBROUTINE exportcreep_grd

  !---------------------------------------------------------------------
  !> subroutine exportCreep_vtk
  !! evaluates the value of creep velocity at the location of 
  !! defined plane (position, strike, dip, length and width).
  !!
  !! input variables
  !! @param np         - number of frictional planes
  !! @param n          - array of frictional planes (position, orientation)
  !! @param structure  - array of depth-dependent frictional properties
  !! @param x0, y0     - origin position of coordinate system
  !! @param dx1,2,3    - sampling size
  !! @param sx1,2,3    - size of the stress tensor field
  !! @param beta       - smoothing factor controlling the extent of planes
  !! @param wdir       - output directory for writing
  !! @param i          - loop index to suffix file names
  !!
  !! \author sylvain barbot (04/04/12) - original form
  !---------------------------------------------------------------------
  SUBROUTINE exportcreep_vtk(np,n,beta,sig,structure, &
                         sx1,sx2,sx3,dx1,dx2,dx3,x0,y0,wdir,i)
    INTEGER, INTENT(IN) :: np,sx1,sx2,sx3,i
    TYPE(PLANE_STRUCT), INTENT(INOUT), DIMENSION(np) :: n
    TYPE(TENSOR), INTENT(IN), DIMENSION(sx1,sx2,sx3) :: sig
    TYPE(LAYER_STRUCT), DIMENSION(:), INTENT(IN) :: structure
    REAL*8, INTENT(IN) :: x0,y0,dx1,dx2,dx3,beta
    CHARACTER(256), INTENT(IN) :: wdir

    INTEGER :: k,ns1,ns2,pos
    CHARACTER(5) :: sdigit
    CHARACTER(3) :: digit
    CHARACTER(256) :: outfile
    INTEGER :: iostatus,i1,i2

    IF (np .le. 0) RETURN

    pos=INDEX(wdir," ")
    WRITE (digit,'(I3.3)') i

    DO k=1,np

       ns1=SIZE(n(k)%patch,1)
       ns2=SIZE(n(k)%patch,2)
          
       WRITE (sdigit,'(I5.5)') k

       outfile=wdir(1:pos-1)//"/creep-"//sdigit//"s-"//digit//".vtk"
       
       OPEN (UNIT=15,FILE=outfile,IOSTAT=iostatus,FORM="FORMATTED")
       IF (iostatus>0) STOP "could not open file for export"
       
       WRITE (15,'("# vtk DataFile Version 3.0")')
       WRITE (15,'("afterslip")')
       WRITE (15,'("ASCII")')
       WRITE (15,'("DATASET POLYDATA")')
       WRITE (15,'("POINTS ",i0," float")') ns1*ns2
       WRITE (15,'(3ES11.3E2)') ((n(k)%patch(i1,i2)%x1,n(k)%patch(i1,i2)%x2,n(k)%patch(i1,i2)%x3,i1=1,ns1), i2=1,ns2)
       WRITE (15,'("")')
       WRITE (15,'("POLYGONS ",i0,X,i0)') (ns1-1)*(ns2-1)*2,(ns1-1)*(ns2-1)*2*4
       DO i2=1,ns2-1
          DO i1=1,ns1-1
             WRITE (15,'("3 ",i0,X,i0,X,i0)') (i2-1)*ns1+i1-1,(i2-1)*ns1+i1-1+1,i2*ns1+i1-1+1
             WRITE (15,'("3 ",i0,X,i0,X,i0)') (i2-1)*ns1+i1-1,(i2  )*ns1+i1-1  ,i2*ns1+i1-1+1
          END DO
       END DO
       WRITE (15,'("")')
       WRITE (15,'("POINT_DATA ",i0)') (ns1)*(ns2)
       WRITE (15,'("SCALARS afterslip float")')
       WRITE (15,'("LOOKUP_TABLE default")')
       WRITE (15,*) ((n(k)%patch(i1,i2)%slip,i1=1,ns1),i2=1,ns2)
       WRITE (15,'("SCALARS strike_slip float")')
       WRITE (15,'("LOOKUP_TABLE default")')
       WRITE (15,*) ((n(k)%patch(i1,i2)%ss,i1=1,ns1),i2=1,ns2)
       WRITE (15,'("SCALARS dip-slip float")')
       WRITE (15,'("LOOKUP_TABLE default")')
       WRITE (15,*) ((n(k)%patch(i1,i2)%ds,i1=1,ns1),i2=1,ns2)
       WRITE (15,'("SCALARS velocity float")')
       WRITE (15,'("LOOKUP_TABLE default")')
       WRITE (15,*) ((n(k)%patch(i1,i2)%v,i1=1,ns1),i2=1,ns2)
       WRITE (15,'("SCALARS velocity_strike_slip float")')
       WRITE (15,'("LOOKUP_TABLE default")')
       WRITE (15,*) ((n(k)%patch(i1,i2)%vss,i1=1,ns1),i2=1,ns2)
       WRITE (15,'("SCALARS velocity_dip_slip float")')
       WRITE (15,'("LOOKUP_TABLE default")')
       WRITE (15,*) ((n(k)%patch(i1,i2)%vds,i1=1,ns1),i2=1,ns2)
       WRITE (15,'("SCALARS shear_stress float")')
       WRITE (15,'("LOOKUP_TABLE default")')
       WRITE (15,*) ((n(k)%patch(i1,i2)%taus,i1=1,ns1),i2=1,ns2)
       CLOSE(15)

    END DO

END SUBROUTINE exportcreep_vtk

#ifdef GRD
  !------------------------------------------------------------------
  !> subroutine ExportStressGRD
  !! writes the 6 components of deformation in map view in the GMT
  !! (Generic Mapping Tools) GRD binary format. This is an interface
  !! to exportgrd.
  !!
  !! \author sylvain barbot 03/19/08 - original form
  !------------------------------------------------------------------
  SUBROUTINE exportstressgrd(sig,sx1,sx2,sx3,dx1,dx2,dx3, &
                             oz,origx,origy,wdir,index)
    INTEGER, INTENT(IN) :: sx1,sx2,sx3,index
    TYPE(TENSOR), INTENT(IN), DIMENSION(sx1,sx2,sx3) :: sig
    REAL*8, INTENT(IN) :: dx1,dx2,dx3,origx,origy,oz
    CHARACTER(256), INTENT(IN) :: wdir

    REAL*4, DIMENSION(:,:), ALLOCATABLE :: t1,t2,t3
    INTEGER :: iostatus,i,j,k,l

    ALLOCATE(t1(sx1+2,sx2),t2(sx1+2,sx2),t3(sx1+2,sx2),STAT=iostatus)
    IF (iostatus>0) STOP "could not allocate memory for grid export"

    k=fix(oz/dx3)+1
    DO j=1,sx2
       DO i=1,sx1
#ifdef ALIGN_DATA
          l=(j-1)*(sx1+2)+i
#else
          l=(j-1)*sx1+i
#endif
          t1(l,1)=sig(i,j,k)%s11
          t2(l,1)=sig(i,j,k)%s12
          t3(l,1)=sig(i,j,k)%s13
       END DO
    END DO

    CALL exportgrd(t1,t2,t3,sx1,sx2,1, &
         dx1,dx2,dx3,0._8,origx,origy,wdir,index,convention=4)

    DO j=1,sx2
       DO i=1,sx1
#ifdef ALIGN_DATA
          l=(j-1)*(sx1+2)+i
#else
          l=(j-1)*sx1+i
#endif
          t1(l,1)=sig(i,j,k)%s22
          t2(l,1)=sig(i,j,k)%s23
          t3(l,1)=sig(i,j,k)%s33
       END DO
    END DO

    CALL exportgrd(t1,t2,t3,sx1,sx2,1, &
         dx1,dx2,dx3,0._8,origx,origy,wdir,index,convention=5)

    DEALLOCATE(t1,t2,t3)

  END SUBROUTINE exportstressgrd


  !------------------------------------------------------------------
  !> subroutine ExportGRD
  !! writes the 3 components of deformation in map view in the GMT
  !! (Generic Mapping Tools) GRD binary format.
  !!
  !! \author sylvain barbot 03/19/08 - original form
  !------------------------------------------------------------------
  SUBROUTINE exportgrd(c1,c2,c3,sx1,sx2,sx3,dx1,dx2,dx3,oz,origx,origy,&
       wdir,i,convention)
    INTEGER, INTENT(IN) :: sx1,sx2,sx3,i
#ifdef ALIGN_DATA
    REAL*4, INTENT(IN), DIMENSION(sx1+2,sx2,sx3) :: c1,c2,c3
#else
    REAL*4, INTENT(IN), DIMENSION(sx1,sx2,sx3) :: c1,c2,c3
#endif
    REAL*8, INTENT(IN) :: dx1,dx2,dx3,origx,origy,oz
    CHARACTER(256), INTENT(IN) :: wdir
    INTEGER, INTENT(IN), OPTIONAL :: convention

    REAL*4, DIMENSION(:,:), ALLOCATABLE :: temp1,temp2,temp3
    REAL*8 :: rland=9998.,rdum=9999.
    INTEGER :: iostatus,k,pos,conv
    REAL*8 :: xmin,ymin
    CHARACTER(256) :: file1,file2,file3
    CHARACTER(3) :: digit

    IF (PRESENT(convention)) THEN
       conv=convention
    ELSE
       conv=1
    END IF

    ALLOCATE(temp1(sx2,sx1),temp2(sx2,sx1),temp3(sx2,sx1),STAT=iostatus)
    IF (iostatus>0) STOP "could not allocate memory for grid export"

    CALL exportspatial(c1(:,:,int(oz/dx3)+1),sx1,sx2,temp1,doflip=.true.)
    CALL exportspatial(c2(:,:,int(oz/dx3)+1),sx1,sx2,temp2,doflip=.true.)
    CALL exportspatial(c3(:,:,int(oz/dx3)+1),sx1,sx2,temp3,doflip=.true.)

    ! positive up
    temp3=-temp3
    
    pos=INDEX(wdir," ")
    WRITE (digit,'(I3.3)') i
    
    SELECT CASE(conv)
    CASE (1) ! cumulative displacement
       file1=wdir(1:pos-1) // "/" // digit // "-north.grd"
       file2=wdir(1:pos-1) // "/" // digit // "-east.grd"
       file3=wdir(1:pos-1) // "/" // digit // "-up.grd"
    CASE (2) ! postseismic displacement
       file1=wdir(1:pos-1) // "/" // digit // "-relax-north.grd"
       file2=wdir(1:pos-1) // "/" // digit // "-relax-east.grd"
       file3=wdir(1:pos-1) // "/" // digit // "-relax-up.grd"
    CASE (3) ! equivalent body forces
       file1=wdir(1:pos-1) // "/" // digit // "-eqbf-north.grd"
       file2=wdir(1:pos-1) // "/" // digit // "-eqbf-east.grd"
       file3=wdir(1:pos-1) // "/" // digit // "-eqbf-up.grd"
    CASE (4) ! equivalent body forces
       file1=wdir(1:pos-1) // "/" // digit // "-s11.grd"
       file2=wdir(1:pos-1) // "/" // digit // "-s12.grd"
       file3=wdir(1:pos-1) // "/" // digit // "-s13.grd"
    CASE (5) ! equivalent body forces
       file1=wdir(1:pos-1) // "/" // digit // "-s22.grd"
       file2=wdir(1:pos-1) // "/" // digit // "-s23.grd"
       file3=wdir(1:pos-1) // "/" // digit // "-s33.grd"
    END SELECT
    
    ! convert to c standard
    k=INDEX(file1," ")
    file1(k:k)=char(0)
    k=INDEX(file2," ")
    file2(k:k)=char(0)
    k=INDEX(file3," ")
    file3(k:k)=char(0)

    ! xmin is the lowest coordinates (positive eastward)
    xmin=origy-sx2/2*dx2
    ! ymin is the lowest coordinates (positive northward)
    ymin=origx-sx1/2*dx1

    ! call the c function "writegrd_"
    CALL writegrd(temp1,sx2,sx1,ymin,xmin,dx1,dx2, &
         rland,rdum,file1,file1)
    CALL writegrd(temp2,sx2,sx1,ymin,xmin,dx1,dx2, &
         rland,rdum,file2,file2)
    CALL writegrd(temp3,sx2,sx1,ymin,xmin,dx1,dx2, &
         rland,rdum,file3,file3)

    DEALLOCATE(temp1,temp2,temp3)

  END SUBROUTINE exportgrd
#endif

#ifdef VTK
  !------------------------------------------------------------------
  !> subroutine ExportVTK_Grid
  !! creates a .vtp file (in the VTK PolyData XML format) containing
  !! the dimension of the computational grid
  !!
  !! \author sylvain barbot 06/24/09 - original form
  !------------------------------------------------------------------
  SUBROUTINE exportvtk_grid(sx1,sx2,sx3,dx1,dx2,dx3,cgfilename)
    INTEGER, INTENT(IN) :: sx1,sx2,sx3
    REAL*8, INTENT(IN) :: dx1,dx2,dx3
    CHARACTER(256), INTENT(IN) :: cgfilename

    INTEGER :: iostatus
    CHARACTER :: q

    q=char(34)

    OPEN (UNIT=15,FILE=cgfilename,IOSTAT=iostatus,FORM="FORMATTED")
    IF (iostatus>0) THEN
       WRITE_DEBUG_INFO
       PRINT '(a)', cgfilename
       STOP "could not open file for export"
    END IF

    WRITE (15,'("<?xml version=",a,"1.0",a,"?>")') q,q
    WRITE (15,'("<VTKFile type=",a,"PolyData",a," version=",a,"0.1",a,">")') q,q,q,q
    WRITE (15,'("  <PolyData>")')
    WRITE (15,'("    <Piece NumberOfPoints=",a,"8",a," NumberOfPolys=",a,"6",a,">")'),q,q,q,q
    WRITE (15,'("      <Points>")')
    WRITE (15,'("        <DataArray type=",a,"Float32",a, &
                           & " Name=",a,"Computational Grid",a, &
                           & " NumberOfComponents=",a,"3",a, &
                           & " format=",a,"ascii",a,">")'),q,q,q,q,q,q,q,q
    WRITE (15,'(24ES9.2E1)') &
                 -sx1*dx1/2, -sx2*dx2/2, sx3*dx3/2, &
                 +sx1*dx1/2, -sx2*dx2/2, sx3*dx3/2, &
                 +sx1*dx1/2, +sx2*dx2/2, sx3*dx3/2, &   
                 -sx1*dx1/2, +sx2*dx2/2, sx3*dx3/2, &
                 -sx1*dx1/2, -sx2*dx2/2, 0.0, &
                 +sx1*dx1/2, -sx2*dx2/2, 0.0, &
                 +sx1*dx1/2, +sx2*dx2/2, 0.0, &
                 -sx1*dx1/2, +sx2*dx2/2, 0.0
    WRITE (15,'("        </DataArray>")')
    WRITE (15,'("      </Points>")')
    WRITE (15,'("      <Polys>")')
    WRITE (15,'("        <DataArray type=",a,"Int32",a, &
                            & " Name=",a,"connectivity",a, &
                            & " format=",a,"ascii",a, &
                            & " RangeMin=",a,"0",a, &
                            & " RangeMax=",a,"7",a,">")'), q,q,q,q,q,q,q,q,q,q
    WRITE (15,'("0 1 2 3 4 5 6 7 2 3 7 6 0 3 7 4 0 1 5 4 1 2 6 5")')
    WRITE (15,'("        </DataArray>")')
    WRITE (15,'("        <DataArray type=",a,"Int32",a, &
                                 & " Name=",a,"offsets",a, &
                                 & " format=",a,"ascii",a, &
                                 & " RangeMin=",a,"4",a, &
                                 & " RangeMax=",a,"24",a,">")'), q,q,q,q,q,q,q,q,q,q
    WRITE (15,'("          4 8 12 16 20 24")')
    WRITE (15,'("        </DataArray>")')
    WRITE (15,'("      </Polys>")')
    WRITE (15,'("    </Piece>")')
    WRITE (15,'("  </PolyData>")')
    WRITE (15,'("</VTKFile>")')

    CLOSE(15)

  END SUBROUTINE exportvtk_grid

  !------------------------------------------------------------------
  !> subroutine ExportXY_RFaults
  !! creates a .xy file (in the GMT closed-polygon format) containing
  !! the rectangular faults. Each fault segemnt is described by a
  !! closed polygon (rectangle) associated with a slip amplitude.
  !! use pxzy with the -Cpalette.cpt -L -M options to color rectangles 
  !! by slip.
  !!
  !! \author sylvain barbot 03/05/11 - original form
  !------------------------------------------------------------------
  SUBROUTINE exportxy_rfaults(e,x0,y0,rffilename)
    TYPE(EVENT_STRUC), INTENT(IN) :: e
    REAL*8, INTENT(IN) :: x0, y0
    CHARACTER(256), INTENT(IN) :: rffilename

    INTEGER :: iostatus,k
    CHARACTER :: q

    REAL*8 :: strike,dip,x1,x2,x3,cstrike,sstrike,cdip,sdip,L,W,slip
         
    REAL*8, DIMENSION(3) :: s,d

    ! double-quote character
    q=char(34)

    OPEN (UNIT=15,FILE=rffilename,IOSTAT=iostatus,FORM="FORMATTED")
    IF (iostatus>0) THEN
       WRITE_DEBUG_INFO
       PRINT '(a)', rffilename
       STOP "could not open file for export"
    END IF

    WRITE (15,'("> # east, north")')
    DO k=1,e%ns

       ! fault slip
       slip=e%s(k)%slip

       ! fault orientation
       strike=e%s(k)%strike
       dip=e%s(k)%dip

       ! fault center position
       x1=e%s(k)%x+x0
       x2=e%s(k)%y+y0
       x3=e%s(k)%z

       ! fault dimension
       W=e%s(k)%width
       L=e%s(k)%length

       cstrike=cos(strike)
       sstrike=sin(strike)
       cdip=cos(dip)
       sdip=sin(dip)
 
       ! strike-slip unit direction
       s(1)=sstrike
       s(2)=cstrike
       s(3)=0._8

       ! dip-slip unit direction
       d(1)=+cstrike*sdip
       d(2)=-sstrike*sdip
       d(3)=+cdip

       ! fault edge coordinates - export east (x2) and north (x1)
       WRITE (15,'("> -Z",3ES11.2)') ABS(slip)
       WRITE (15,'(3ES11.2)') x2-d(2)*W/2-s(2)*L/2, x1-d(1)*W/2-s(1)*L/2
       WRITE (15,'(3ES11.2)') x2-d(2)*W/2+s(2)*L/2, x1-d(1)*W/2+s(1)*L/2
       WRITE (15,'(3ES11.2)') x2+d(2)*W/2+s(2)*L/2, x1+d(1)*W/2+s(1)*L/2
       WRITE (15,'(3ES11.2)') x2+d(2)*W/2-s(2)*L/2, x1+d(1)*W/2-s(1)*L/2

    END DO

    CLOSE(15)

  END SUBROUTINE exportxy_rfaults

  !------------------------------------------------------------------
  !> subroutine ExportVTK_RFaults
  !! creates a .vtp file (in the VTK PolyData XML format) containing
  !! the rectangular faults. The faults are characterized with a set
  !! of subsegments (rectangles) each associated with a slip vector. 
  !!
  !! \author sylvain barbot 06/24/09 - original form
  !------------------------------------------------------------------
  SUBROUTINE exportvtk_rfaults(e,rffilename)
    TYPE(EVENT_STRUC), INTENT(IN) :: e
    CHARACTER(256), INTENT(IN) :: rffilename

    INTEGER :: iostatus,k
    CHARACTER :: q

    REAL*8 :: slip,opening
    REAL*8 :: strike,dip,x1,x2,x3,cstrike,sstrike,cdip,sdip,L,W,r
         
    REAL*8, DIMENSION(3) :: s,d,n

    ! double-quote character
    q=char(34)

    OPEN (UNIT=15,FILE=rffilename,IOSTAT=iostatus,FORM="FORMATTED")
    IF (iostatus>0) THEN
       WRITE_DEBUG_INFO
       PRINT '(a)', rffilename
       STOP "could not open file for export"
    END IF

    WRITE (15,'("<?xml version=",a,"1.0",a,"?>")') q,q
    WRITE (15,'("<VTKFile type=",a,"PolyData",a," version=",a,"0.1",a,">")') q,q,q,q
    WRITE (15,'("  <PolyData>")')

    DO k=1,e%ns

       ! dyke opening
       opening=e%s(k)%opening

       ! fault slip
       slip=e%s(k)%slip

       ! slip
       r=e%s(k)%rake

       ! fault orientation
       strike=e%s(k)%strike
       dip=e%s(k)%dip

       ! fault center position
       x1=e%s(k)%x
       x2=e%s(k)%y
       x3=e%s(k)%z

       ! fault dimension
       W=e%s(k)%width
       L=e%s(k)%length

       cstrike=cos(strike)
       sstrike=sin(strike)
       cdip=cos(dip)
       sdip=sin(dip)
 
       ! strike-slip unit direction
       s(1)=sstrike
       s(2)=cstrike
       s(3)=0._8

       ! dip-slip unit direction
       d(1)=+cstrike*sdip
       d(2)=-sstrike*sdip
       d(3)=+cdip

       ! surface normal vector components
       n(1)=+cdip*cstrike
       n(2)=-cdip*sstrike
       n(3)=-sdip

       WRITE (15,'("    <Piece NumberOfPoints=",a,"4",a," NumberOfPolys=",a,"1",a,">")'),q,q,q,q
       WRITE (15,'("      <Points>")')
       WRITE (15,'("        <DataArray type=",a,"Float32",a, &
                           & " Name=",a,"Fault Patch",a, &
                           & " NumberOfComponents=",a,"3",a, &
                           & " format=",a,"ascii",a,">")'),q,q,q,q,q,q,q,q

       ! fault edge coordinates
       WRITE (15,'(12ES11.2)') &
                     x1-d(1)*W/2-s(1)*L/2,x2-d(2)*W/2-s(2)*L/2,x3-d(3)*W/2-s(3)*L/2, &
                     x1-d(1)*W/2+s(1)*L/2,x2-d(2)*W/2+s(2)*L/2,x3-d(3)*W/2+s(3)*L/2, &
                     x1+d(1)*W/2+s(1)*L/2,x2+d(2)*W/2+s(2)*L/2,x3+d(3)*W/2+s(3)*L/2, &
                     x1+d(1)*W/2-s(1)*L/2,x2+d(2)*W/2-s(2)*L/2,x3+d(3)*W/2-s(3)*L/2

       WRITE (15,'("        </DataArray>")')
       WRITE (15,'("      </Points>")')
       WRITE (15,'("      <Polys>")')
       WRITE (15,'("        <DataArray type=",a,"Int32",a, &
                            & " Name=",a,"connectivity",a, &
                            & " format=",a,"ascii",a, &
                            & " RangeMin=",a,"0",a, &
                            & " RangeMax=",a,"3",a,">")'), q,q,q,q,q,q,q,q,q,q
       WRITE (15,'("0 1 2 3")')
       WRITE (15,'("        </DataArray>")')
       WRITE (15,'("        <DataArray type=",a,"Int32",a, &
                                 & " Name=",a,"offsets",a, &
                                 & " format=",a,"ascii",a, &
                                 & " RangeMin=",a,"4",a, &
                                 & " RangeMax=",a,"4",a,">")'), q,q,q,q,q,q,q,q,q,q
       WRITE (15,'("          4")')
       WRITE (15,'("        </DataArray>")')
       WRITE (15,'("      </Polys>")')

       WRITE (15,'("      <CellData Normals=",a,"slip",a,">")'), q,q
       WRITE (15,'("        <DataArray type=",a,"Float32",a, &
                               & " Name=",a,"slip",a, &
                               & " NumberOfComponents=",a,"3",a, &
                               & " format=",a,"ascii",a,">")'), q,q,q,q,q,q,q,q


       WRITE (15,'(3ES11.2)'), (s(1)*cos(r)+d(1)*sin(r))*slip, &
                               (s(2)*cos(r)+d(2)*sin(r))*slip, &
                               (s(3)*cos(r)+d(3)*sin(r))*slip

       WRITE (15,'("        </DataArray>")')

       WRITE (15,'("        <DataArray type=",a,"Float32",a, &
                               & " Name=",a,"opening",a, &
                               & " NumberOfComponents=",a,"3",a, &
                               & " format=",a,"ascii",a,">")'), q,q,q,q,q,q,q,q


       WRITE (15,'(3ES11.2)'), n(1)*opening, &
                               n(2)*opening, &
                               n(3)*opening

       WRITE (15,'("        </DataArray>")')
       WRITE (15,'("      </CellData>")')

       WRITE (15,'("    </Piece>")')

    END DO

    WRITE (15,'("  </PolyData>")')
    WRITE (15,'("</VTKFile>")')

    CLOSE(15)

  END SUBROUTINE exportvtk_rfaults

  !------------------------------------------------------------------
  !> subroutine ExportVTK_RFaults_Stress_Init
  !! creates a .vtp file (in the VTK PolyData XML format) containing
  !! the rectangular faults. The faults are characterized with a set
  !! of subsegments (rectangles) each associated with stress values. 
  !!
  !! \author sylvain barbot 06/06/11 - original form
  !------------------------------------------------------------------
  SUBROUTINE export_rfaults_stress_init(sig,sx1,sx2,sx3,dx1,dx2,dx3, &
                                           nsop,sop)
    INTEGER, INTENT(IN) :: sx1,sx2,sx3,nsop
    TYPE(TENSOR), INTENT(IN), DIMENSION(sx1,sx2,sx3) :: sig
    REAL*8, INTENT(IN) :: dx1,dx2,dx3
    TYPE(SEGMENT_STRUCT), INTENT(INOUT), DIMENSION(nsop) :: sop

    INTEGER :: k,i1,i2,i3
    REAL*8 :: x1,x2,x3
    ! local value of stress
    TYPE(TENSOR) :: lsig

    DO k=1,nsop
       ! fault center position
       x1=sop(k)%x
       x2=sop(k)%y
       x3=sop(k)%z

       CALL shiftedindex(x1,x2,x3,sx1,sx2,sx3,dx1,dx2,dx3,i1,i2,i3)
       lsig=sig(i1,i2,i3)

       sop(k)%sig0%s11=lsig%s11
       sop(k)%sig0%s12=lsig%s12
       sop(k)%sig0%s13=lsig%s13
       sop(k)%sig0%s22=lsig%s22
       sop(k)%sig0%s23=lsig%s23
       sop(k)%sig0%s33=lsig%s33

    END DO

  END SUBROUTINE export_rfaults_stress_init

  !------------------------------------------------------------------
  !> subroutine ExportGMT_RFaults_Stress
  !! creates a .vtp file (in the VTK PolyData XML format) containing
  !! the rectangular faults. The faults are characterized with a set
  !! of subsegments (rectangles) each associated with stress values. 
  !!
  !! \author sylvain barbot 06/06/11 - original form
  !------------------------------------------------------------------
  SUBROUTINE exportgmt_rfaults_stress(sx1,sx2,sx3,dx1,dx2,dx3, &
                          nsop,sop,rffilename,convention,sig)
    USE elastic3d
    INTEGER, INTENT(IN) :: sx1,sx2,sx3,nsop
    REAL*8, INTENT(IN) :: dx1,dx2,dx3
    TYPE(SEGMENT_STRUCT), INTENT(INOUT), DIMENSION(nsop) :: sop
    CHARACTER(256), INTENT(IN) :: rffilename
    INTEGER, INTENT(IN), OPTIONAL :: convention
    TYPE(TENSOR), INTENT(IN), DIMENSION(sx1,sx2,sx3), OPTIONAL :: sig

    INTEGER :: iostatus,k,i1,i2,i3,conv
    CHARACTER :: q
    REAL*8 :: strike,dip,x1,x2,x3,cstrike,sstrike,cdip,sdip,L,W
    ! segment normal vector, strike direction, dip direction
    REAL*8, DIMENSION(3) :: n,s,d
    ! local value of stress
    TYPE(TENSOR) :: lsig
    ! stress components
    REAL*8 :: taun,taus,taustrike,taudip,taucoulomb
    ! friction coefficient
    REAL*8 :: friction
    ! traction components
    REAL*8, DIMENSION(3) :: t,ts

    IF (0.GE.nsop) RETURN

    ! double-quote character
    q=char(34)

    IF (PRESENT(convention)) THEN
       conv=convention
    ELSE
       conv=0
    END IF

    OPEN (UNIT=15,FILE=rffilename,IOSTAT=iostatus,FORM="FORMATTED")
    IF (iostatus>0) THEN
       WRITE_DEBUG_INFO
       PRINT '(a)', rffilename
       STOP "could not open file for export"
    END IF

    DO k=1,nsop
       ! friction coefficient
       friction=sop(k)%friction

       ! fault orientation
       strike=sop(k)%strike
       dip=sop(k)%dip

       ! fault center position
       x1=sop(k)%x
       x2=sop(k)%y
       x3=sop(k)%z

       IF (PRESENT(sig)) THEN

          CALL shiftedindex(x1,x2,x3,sx1,sx2,sx3,dx1,dx2,dx3,i1,i2,i3)
          lsig=sig(i1,i2,i3)

          IF (1.EQ.conv) THEN
             lsig%s11=lsig%s11-sop(k)%sig0%s11
             lsig%s12=lsig%s12-sop(k)%sig0%s12
             lsig%s13=lsig%s13-sop(k)%sig0%s13
             lsig%s22=lsig%s22-sop(k)%sig0%s22
             lsig%s23=lsig%s23-sop(k)%sig0%s23
             lsig%s33=lsig%s33-sop(k)%sig0%s33
          END IF
       ELSE
          lsig=TENSOR(0.d0,0.d0,0.d0,0.d0,0.d0,0.d0)
       END IF

       ! fault dimension
       W=sop(k)%width
       L=sop(k)%length

       cstrike=cos(strike)
       sstrike=sin(strike)
       cdip=cos(dip)
       sdip=sin(dip)
 
       ! surface normal vector components
       n(1)=+cdip*cstrike
       n(2)=-cdip*sstrike
       n(3)=-sdip

       ! strike-slip unit direction
       s(1)=sstrike
       s(2)=cstrike
       s(3)=0._8

       ! dip-slip unit direction
       d(1)=+cstrike*sdip
       d(2)=-sstrike*sdip
       d(3)=+cdip

       ! traction vector
       t=lsig .tdot. n

       ! signed normal component
       taun=SUM(t*n)

       ! shear traction
       ts=t-taun*n

       ! absolute value of shear component
       taus=SQRT(SUM(ts*ts))

       ! strike-direction shear component
       taustrike=SUM(ts*s)

       ! dip-direction shear component
       taudip=SUM(ts*d)

       ! Coulomb stress 
       taucoulomb=taus+friction*taun

       WRITE (15,'("> -Z",5ES11.2)') taus, taun, taucoulomb, taustrike, taudip
       WRITE (15,'(3ES11.2)') x1-d(1)*W/2-s(1)*L/2, x2-d(2)*W/2-s(2)*L/2
       WRITE (15,'(3ES11.2)') x1-d(1)*W/2+s(1)*L/2, x2-d(2)*W/2+s(2)*L/2
       WRITE (15,'(3ES11.2)') x1+d(1)*W/2+s(1)*L/2, x2+d(2)*W/2+s(2)*L/2
       WRITE (15,'(3ES11.2)') x1+d(1)*W/2-s(1)*L/2, x2+d(2)*W/2-s(2)*L/2

    END DO

    CLOSE(15)

  END SUBROUTINE exportgmt_rfaults_stress

  !------------------------------------------------------------------
  !> subroutine ExportVTK_RFaults_Stress
  !! creates a .vtp file (in the VTK PolyData XML format) containing
  !! the rectangular faults. The faults are characterized with a set
  !! of subsegments (rectangles) each associated with stress values. 
  !!
  !! \author sylvain barbot 06/06/11 - original form
  !------------------------------------------------------------------
  SUBROUTINE exportvtk_rfaults_stress(sx1,sx2,sx3,dx1,dx2,dx3, &
                          nsop,sop,rffilename,convention,sig)
    USE elastic3d
    INTEGER, INTENT(IN) :: sx1,sx2,sx3,nsop
    REAL*8, INTENT(IN) :: dx1,dx2,dx3
    TYPE(SEGMENT_STRUCT), INTENT(INOUT), DIMENSION(nsop) :: sop
    CHARACTER(256), INTENT(IN) :: rffilename
    INTEGER, INTENT(IN), OPTIONAL :: convention
    TYPE(TENSOR), INTENT(IN), DIMENSION(sx1,sx2,sx3), OPTIONAL :: sig

    INTEGER :: iostatus,k,i1,i2,i3,conv
    CHARACTER :: q
    REAL*8 :: strike,dip,x1,x2,x3,cstrike,sstrike,cdip,sdip,L,W
    ! segment normal vector, strike direction, dip direction
    REAL*8, DIMENSION(3) :: n,s,d
    ! local value of stress
    TYPE(TENSOR) :: lsig
    ! stress components
    REAL*8 :: taun,taus,taustrike,taudip,taucoulomb
    ! friction coefficient
    REAL*8 :: friction
    ! traction components
    REAL*8, DIMENSION(3) :: t,ts

    IF (0.GE.nsop) RETURN

    ! double-quote character
    q=char(34)

    IF (PRESENT(convention)) THEN
       conv=convention
    ELSE
       conv=0
    END IF

    OPEN (UNIT=15,FILE=rffilename,IOSTAT=iostatus,FORM="FORMATTED")
    IF (iostatus>0) THEN
       WRITE_DEBUG_INFO
       PRINT '(a)', rffilename
       STOP "could not open file for export"
    END IF

    WRITE (15,'("<?xml version=",a,"1.0",a,"?>")') q,q
    WRITE (15,'("<VTKFile type=",a,"PolyData",a," version=",a,"0.1",a,">")') q,q,q,q
    WRITE (15,'("  <PolyData>")')

    DO k=1,nsop
       ! friction coefficient
       friction=sop(k)%friction

       ! fault orientation
       strike=sop(k)%strike
       dip=sop(k)%dip

       ! fault center position
       x1=sop(k)%x
       x2=sop(k)%y
       x3=sop(k)%z

       IF (PRESENT(sig)) THEN

          CALL shiftedindex(x1,x2,x3,sx1,sx2,sx3,dx1,dx2,dx3,i1,i2,i3)
          lsig=sig(i1,i2,i3)

          IF (1.EQ.conv) THEN
             lsig%s11=lsig%s11-sop(k)%sig0%s11
             lsig%s12=lsig%s12-sop(k)%sig0%s12
             lsig%s13=lsig%s13-sop(k)%sig0%s13
             lsig%s22=lsig%s22-sop(k)%sig0%s22
             lsig%s23=lsig%s23-sop(k)%sig0%s23
             lsig%s33=lsig%s33-sop(k)%sig0%s33
          END IF
       ELSE
          lsig=TENSOR(0.d0,0.d0,0.d0,0.d0,0.d0,0.d0)
       END IF

       ! fault dimension
       W=sop(k)%width
       L=sop(k)%length

       cstrike=cos(strike)
       sstrike=sin(strike)
       cdip=cos(dip)
       sdip=sin(dip)
 
       ! surface normal vector components
       n(1)=+cdip*cstrike
       n(2)=-cdip*sstrike
       n(3)=-sdip

       ! strike-slip unit direction
       s(1)=sstrike
       s(2)=cstrike
       s(3)=0._8

       ! dip-slip unit direction
       d(1)=+cstrike*sdip
       d(2)=-sstrike*sdip
       d(3)=+cdip

       ! traction vector
       t=lsig .tdot. n

       ! signed normal component
       taun=SUM(t*n)

       ! shear traction
       ts=t-taun*n

       ! absolute value of shear component
       taus=SQRT(SUM(ts*ts))

       ! strike-direction shear component
       taustrike=SUM(ts*s)

       ! dip-direction shear component
       taudip=SUM(ts*d)

       ! Coulomb stress 
       taucoulomb=taus+friction*taun

       WRITE (15,'("    <Piece NumberOfPoints=",a,"4",a," NumberOfPolys=",a,"1",a,">")'),q,q,q,q
       WRITE (15,'("      <Points>")')
       WRITE (15,'("        <DataArray type=",a,"Float32",a, &
                           & " Name=",a,"Fault Patch",a, &
                           & " NumberOfComponents=",a,"3",a, &
                           & " format=",a,"ascii",a,">")'),q,q,q,q,q,q,q,q
       ! fault edge coordinates
       WRITE (15,'(12ES11.2)') &
                     x1-d(1)*W/2-s(1)*L/2, x2-d(2)*W/2-s(2)*L/2, x3-d(3)*W/2-s(3)*L/2, &
                     x1-d(1)*W/2+s(1)*L/2, x2-d(2)*W/2+s(2)*L/2, x3-d(3)*W/2+s(3)*L/2, &
                     x1+d(1)*W/2+s(1)*L/2, x2+d(2)*W/2+s(2)*L/2, x3+d(3)*W/2+s(3)*L/2, &
                     x1+d(1)*W/2-s(1)*L/2, x2+d(2)*W/2-s(2)*L/2, x3+d(3)*W/2-s(3)*L/2
       WRITE (15,'("        </DataArray>")')

       WRITE (15,'("      </Points>")')
       WRITE (15,'("      <Polys>")')
       WRITE (15,'("        <DataArray type=",a,"Int32",a, &
                            & " Name=",a,"connectivity",a, &
                            & " format=",a,"ascii",a, &
                            & " RangeMin=",a,"0",a, &
                            & " RangeMax=",a,"3",a,">")'), q,q,q,q,q,q,q,q,q,q
       WRITE (15,'("0 1 2 3")')
       WRITE (15,'("        </DataArray>")')
       WRITE (15,'("        <DataArray type=",a,"Int32",a, &
                                 & " Name=",a,"offsets",a, &
                                 & " format=",a,"ascii",a, &
                                 & " RangeMin=",a,"4",a, &
                                 & " RangeMax=",a,"4",a,">")'), q,q,q,q,q,q,q,q,q,q
       WRITE (15,'("          4")')
       WRITE (15,'("        </DataArray>")')
       WRITE (15,'("      </Polys>")')

       WRITE (15,'("      <CellData Normals=",a,"stress",a,">")'), q,q

       WRITE (15,'("        <DataArray type=",a,"Float32",a, &
                               & " Name=",a,"stress tensor",a, &
                               & " NumberOfComponents=",a,"6",a, &
                               & " format=",a,"ascii",a,">")'), q,q,q,q,q,q,q,q
       WRITE (15,'(6ES11.2)'), lsig%s11,lsig%s12,lsig%s13,lsig%s22,lsig%s23,lsig%s33
       WRITE (15,'("        </DataArray>")')

       WRITE (15,'("        <DataArray type=",a,"Float32",a, &
                               & " Name=",a,"shear stress",a, &
                               & " NumberOfComponents=",a,"1",a, &
                               & " format=",a,"ascii",a,">")'), q,q,q,q,q,q,q,q
       WRITE (15,'(ES11.2)'), taus
       WRITE (15,'("        </DataArray>")')

       WRITE (15,'("        <DataArray type=",a,"Float32",a, &
                               & " Name=",a,"normal stress",a, &
                               & " NumberOfComponents=",a,"1",a, &
                               & " format=",a,"ascii",a,">")'), q,q,q,q,q,q,q,q
       WRITE (15,'(ES11.2)'), taun
       WRITE (15,'("        </DataArray>")')

       WRITE (15,'("        <DataArray type=",a,"Float32",a, &
                               & " Name=",a,"Coulomb stress",a, &
                               & " NumberOfComponents=",a,"1",a, &
                               & " format=",a,"ascii",a,">")'), q,q,q,q,q,q,q,q
       WRITE (15,'(ES11.2)'), taucoulomb
       WRITE (15,'("        </DataArray>")')

       WRITE (15,'("        <DataArray type=",a,"Float32",a, &
                               & " Name=",a,"stress in strike direction",a, &
                               & " NumberOfComponents=",a,"1",a, &
                               & " format=",a,"ascii",a,">")'), q,q,q,q,q,q,q,q
       WRITE (15,'(ES11.2)'), taustrike
       WRITE (15,'("        </DataArray>")')

       WRITE (15,'("        <DataArray type=",a,"Float32",a, &
                               & " Name=",a,"stress in dip direction",a, &
                               & " NumberOfComponents=",a,"1",a, &
                               & " format=",a,"ascii",a,">")'), q,q,q,q,q,q,q,q
       WRITE (15,'(ES11.2)'), taudip
       WRITE (15,'("        </DataArray>")')

       WRITE (15,'("      </CellData>")')

       WRITE (15,'("    </Piece>")')

    END DO

    WRITE (15,'("  </PolyData>")')
    WRITE (15,'("</VTKFile>")')

    CLOSE(15)

  END SUBROUTINE exportvtk_rfaults_stress

  !--------------------------------------------------------------------------------
  !> subroutine ExportCoulombStress
  !! sample the stress tensor, shear and normal stress and Coulomb
  !! stress at a series of locations.
  !!
  !! each fault patch is attributed to a file in which the time 
  !! evolution is listed in the following format:
  !!
  !! #t     s11     s12     s13     s22     s23     s33     taus     taud     tau     taun     Coulomb
  !! t0 s11(t0) s12(t0) s13(t0) s22(t0) s23(t0) s33(t0) taus(t0) taud(t0) tau(t0) taun(t0) Coulomb(t0)
  !! t1 s11(t1) s12(t1) s13(t1) s22(t1) s23(t1) s33(t1) taus(t1) taud(t1) tau(t1) taun(t1) Coulomb(t0)
  !!    ...
  !!
  !! where sij(t0) is the component ij of the stress tensor at time t0, taus is
  !! the component of shear in the strike direction, taud is the component of shear
  !! in the fault dip direction, tau^2=taus^2+taud^2, taun is the fault normal
  !! stress and Coulomb(t0) is the Coulomb stress tau+mu*taun. 
  !!
  !! \author sylvain barbot (10/05/11) - original form
  !--------------------------------------------------------------------------------
  SUBROUTINE exportcoulombstress(sig,sx1,sx2,sx3,dx1,dx2,dx3, &
                          nsop,sop,time,wdir,isnew)
    INTEGER, INTENT(IN) :: sx1,sx2,sx3,nsop
    TYPE(TENSOR), INTENT(IN), DIMENSION(sx1,sx2,sx3) :: sig
    TYPE(SEGMENT_STRUCT), INTENT(INOUT), DIMENSION(nsop) :: sop
    REAL*8, INTENT(IN) :: dx1,dx2,dx3,time
    CHARACTER(256), INTENT(IN) :: wdir
    LOGICAL, INTENT(IN) :: isnew

    INTEGER :: iostatus,k,i1,i2,i3
    CHARACTER :: q
    CHARACTER(4) :: digit4
    CHARACTER(256) :: file
    REAL*8 :: strike,dip,x1,x2,x3,cstrike,sstrike,cdip,sdip,L,W
    ! segment normal vector, strike direction, dip direction
    REAL*8, DIMENSION(3) :: n,s,d
    ! local value of stress
    TYPE(TENSOR) :: lsig
    ! stress components
    REAL*8 :: taun,taus,taustrike,taudip,taucoulomb
    ! friction coefficient
    REAL*8 :: friction
    ! traction components
    REAL*8, DIMENSION(3) :: t,ts

    IF (0.GE.nsop) RETURN

    ! double-quote character
    q=char(34)

    DO k=1,nsop
       WRITE (digit4,'(I4.4)') k
       file=trim(wdir)//"/cfaults-sigma-"//digit4//".txt"

       ! fault center position
       x1=sop(k)%x
       x2=sop(k)%y
       x3=sop(k)%z

       IF (isnew) THEN
          OPEN (UNIT=15,FILE=file,IOSTAT=iostatus,FORM="FORMATTED")
          WRITE (15,'("# center position (north, east, down): ",3ES9.2)') x1,x2,x3
          WRITE (15,'("#         t        s11        s12        s13        ", &
               & "s22        s23        s33       taus       taud        tau       taun    Coulomb")')
       ELSE
          OPEN (UNIT=15,FILE=file,POSITION="APPEND",&
               IOSTAT=iostatus,FORM="FORMATTED")
       END IF
       IF (iostatus>0) STOP "could not open point file for writing"

       ! friction coefficient
       friction=sop(k)%friction

       ! fault orientation
       strike=sop(k)%strike
       dip=sop(k)%dip

       CALL shiftedindex(x1,x2,x3,sx1,sx2,sx3,dx1,dx2,dx3,i1,i2,i3)
       lsig=sig(i1,i2,i3)

       ! fault dimension
       W=sop(k)%width
       L=sop(k)%length

       cstrike=cos(strike)
       sstrike=sin(strike)
       cdip=cos(dip)
       sdip=sin(dip)
 
       ! surface normal vector components
       n(1)=+cdip*cstrike
       n(2)=-cdip*sstrike
       n(3)=-sdip

       ! strike-slip unit direction
       s(1)=sstrike
       s(2)=cstrike
       s(3)=0._8

       ! dip-slip unit direction
       d(1)=+cstrike*sdip
       d(2)=-sstrike*sdip
       d(3)=+cdip

       ! traction vector
       t=lsig .tdot. n

       ! signed normal component
       taun=SUM(t*n)

       ! shear traction
       ts=t-taun*n

       ! absolute value of shear component
       taus=SQRT(SUM(ts*ts))

       ! strike-direction shear component
       taustrike=SUM(ts*s)

       ! dip-direction shear component
       taudip=SUM(ts*d)

       ! Coulomb stress 
       taucoulomb=taus+friction*taun

       WRITE (15,'(12ES11.3E2)') time, &
                                 lsig%s11,lsig%s12,lsig%s13, &
                                 lsig%s22,lsig%s23,lsig%s33, &
                                 taustrike,taudip,taus,taun,taucoulomb
       CLOSE(15)
    END DO

  END SUBROUTINE exportcoulombstress

  !------------------------------------------------------------------
  !> subroutine ExportVTK_Rectangle
  !! creates a .vtp file (in the VTK PolyData XML format) containing
  !! a rectangle.
  !!
  !! \author sylvain barbot 06/24/09 - original form
  !------------------------------------------------------------------
  SUBROUTINE exportvtk_rectangle(x1,x2,x3,L,W,strike,dip,filename)
    REAL*8 :: x1,x2,x3,L,W,strike,dip
    CHARACTER(256), INTENT(IN) :: filename

    INTEGER :: iostatus
    CHARACTER :: q

    REAL*8 :: cstrike,sstrike,cdip,sdip
    REAL*8, DIMENSION(3) :: s,d

    ! double-quote character
    q=char(34)

    OPEN (UNIT=15,FILE=filename,IOSTAT=iostatus,FORM="FORMATTED")
    IF (iostatus>0) THEN
       WRITE_DEBUG_INFO
       PRINT '(a)', filename
       STOP "could not open file for export in ExportVTK_Rectangle"
    END IF

    WRITE (15,'("<?xml version=",a,"1.0",a,"?>")') q,q
    WRITE (15,'("<VTKFile type=",a,"PolyData",a," version=",a,"0.1",a,">")') q,q,q,q
    WRITE (15,'("  <PolyData>")')

    cstrike=cos(strike)
    sstrike=sin(strike)
    cdip=cos(dip)
    sdip=sin(dip)
 
    ! strike-slip unit direction
    s(1)=sstrike
    s(2)=cstrike
    s(3)=0._8

    ! dip-slip unit direction
    d(1)=+cstrike*sdip
    d(2)=-sstrike*sdip
    d(3)=+cdip

    WRITE (15,'("    <Piece NumberOfPoints=",a,"4",a," NumberOfPolys=",a,"1",a,">")'),q,q,q,q
    WRITE (15,'("      <Points>")')
    WRITE (15,'("        <DataArray type=",a,"Float32",a, &
                        & " Name=",a,"Fault Patch",a, &
                        & " NumberOfComponents=",a,"3",a, &
                        & " format=",a,"ascii",a,">")'),q,q,q,q,q,q,q,q

    ! fault edge coordinates
    WRITE (15,'(12ES11.2)') &
                  x1-d(1)*W/2-s(1)*L/2,x2-d(2)*W/2-s(2)*L/2,x3-d(3)*W/2-s(3)*L/2, &
                  x1-d(1)*W/2+s(1)*L/2,x2-d(2)*W/2+s(2)*L/2,x3-d(3)*W/2+s(3)*L/2, &
                  x1+d(1)*W/2+s(1)*L/2,x2+d(2)*W/2+s(2)*L/2,x3+d(3)*W/2+s(3)*L/2, &
                  x1+d(1)*W/2-s(1)*L/2,x2+d(2)*W/2-s(2)*L/2,x3+d(3)*W/2-s(3)*L/2

    WRITE (15,'("        </DataArray>")')
    WRITE (15,'("      </Points>")')
    WRITE (15,'("      <Polys>")')
    WRITE (15,'("        <DataArray type=",a,"Int32",a, &
                         & " Name=",a,"connectivity",a, &
                         & " format=",a,"ascii",a, &
                         & " RangeMin=",a,"0",a, &
                         & " RangeMax=",a,"3",a,">")'), q,q,q,q,q,q,q,q,q,q
    WRITE (15,'("0 1 2 3")')
    WRITE (15,'("        </DataArray>")')
    WRITE (15,'("        <DataArray type=",a,"Int32",a, &
                              & " Name=",a,"offsets",a, &
                              & " format=",a,"ascii",a, &
                              & " RangeMin=",a,"4",a, &
                              & " RangeMax=",a,"4",a,">")'), q,q,q,q,q,q,q,q,q,q
    WRITE (15,'("          4")')
    WRITE (15,'("        </DataArray>")')
    WRITE (15,'("      </Polys>")')

    WRITE (15,'("    </Piece>")')

    WRITE (15,'("  </PolyData>")')
    WRITE (15,'("</VTKFile>")')

    CLOSE(15)

  END SUBROUTINE exportvtk_rectangle

  !------------------------------------------------------------------
  !> subroutine ExportXY_Brick
  !! creates a .xy file containing a brick (3d rectangle, cuboid).
  !!
  !! \author sylvain barbot 11/29/11 - original form
  !------------------------------------------------------------------
  SUBROUTINE exportxy_brick(x1,x2,x3,L,W,T,strike,dip,filename)
    REAL*8 :: x1,x2,x3,L,W,T,strike,dip
    CHARACTER(256), INTENT(IN) :: filename

    INTEGER :: iostatus

    REAL*8 :: cstrike,sstrike,cdip,sdip
    REAL*8, DIMENSION(3) :: s,d,n

    OPEN (UNIT=15,FILE=filename,IOSTAT=iostatus,FORM="FORMATTED")
    IF (iostatus>0) THEN
       WRITE_DEBUG_INFO
       PRINT '(a)', filename
       STOP "could not open file for export in ExportXY_Brick"
    END IF

    cstrike=cos(strike)
    sstrike=sin(strike)
    cdip=cos(dip)
    sdip=sin(dip)
 
    ! strike-slip unit direction
    s(1)=sstrike
    s(2)=cstrike
    s(3)=0._8

    ! dip-slip unit direction
    d(1)=+cstrike*sdip
    d(2)=-sstrike*sdip
    d(3)=+cdip

    ! surface normal vector components
    n(1)=+cdip*cstrike
    n(2)=-cdip*sstrike
    n(3)=-sdip

    ! fault edge coordinates
    WRITE (15,'(">")')
    WRITE (15,'(3ES11.2)') x2-d(2)*W/2-s(2)*L/2-n(2)*T/2.d0, x1-d(1)*W/2-s(1)*L/2-n(1)*T/2.d0, x3-d(3)*W/2-s(3)*L/2-n(3)*T/2.d0
    WRITE (15,'(3ES11.2)') x2-d(2)*W/2+s(2)*L/2-n(2)*T/2.d0, x1-d(1)*W/2+s(1)*L/2-n(1)*T/2.d0, x3-d(3)*W/2+s(3)*L/2-n(3)*T/2.d0
    WRITE (15,'(3ES11.2)') x2-d(2)*W/2+s(2)*L/2+n(2)*T/2.d0, x1-d(1)*W/2+s(1)*L/2+n(1)*T/2.d0, x3-d(3)*W/2+s(3)*L/2+n(3)*T/2.d0
    WRITE (15,'(3ES11.2)') x2-d(2)*W/2-s(2)*L/2+n(2)*T/2.d0, x1-d(1)*W/2-s(1)*L/2+n(1)*T/2.d0, x3-d(3)*W/2-s(3)*L/2+n(3)*T/2.d0
    WRITE (15,'(3ES11.2)') x2-d(2)*W/2-s(2)*L/2-n(2)*T/2.d0, x1-d(1)*W/2-s(1)*L/2-n(1)*T/2.d0, x3-d(3)*W/2-s(3)*L/2-n(3)*T/2.d0
    WRITE (15,'(">")')
    WRITE (15,'(3ES11.2)') x2+d(2)*W/2-s(2)*L/2-n(2)*T/2.d0, x1+d(1)*W/2-s(1)*L/2-n(1)*T/2.d0, x3+d(3)*W/2-s(3)*L/2-n(3)*T/2.d0
    WRITE (15,'(3ES11.2)') x2+d(2)*W/2+s(2)*L/2-n(2)*T/2.d0, x1+d(1)*W/2+s(1)*L/2-n(1)*T/2.d0, x3+d(3)*W/2+s(3)*L/2-n(3)*T/2.d0
    WRITE (15,'(3ES11.2)') x2+d(2)*W/2+s(2)*L/2+n(2)*T/2.d0, x1+d(1)*W/2+s(1)*L/2+n(1)*T/2.d0, x3+d(3)*W/2+s(3)*L/2+n(3)*T/2.d0
    WRITE (15,'(3ES11.2)') x2+d(2)*W/2-s(2)*L/2+n(2)*T/2.d0, x1+d(1)*W/2-s(1)*L/2+n(1)*T/2.d0, x3+d(3)*W/2-s(3)*L/2+n(3)*T/2.d0
    WRITE (15,'(3ES11.2)') x2+d(2)*W/2-s(2)*L/2-n(2)*T/2.d0, x1+d(1)*W/2-s(1)*L/2-n(1)*T/2.d0, x3+d(3)*W/2-s(3)*L/2-n(3)*T/2.d0

    CLOSE(15)

  END SUBROUTINE exportxy_brick

  !------------------------------------------------------------------
  !> subroutine ExportVTK_Brick
  !! creates a .vtp file (in the VTK PolyData XML format) containing
  !! a brick (3d rectangle, cuboid).
  !!
  !! \author sylvain barbot 06/24/09 - original form
  !------------------------------------------------------------------
  SUBROUTINE exportvtk_brick(x1,x2,x3,L,W,T,strike,dip,filename)
    REAL*8 :: x1,x2,x3,L,W,T,strike,dip
    CHARACTER(256), INTENT(IN) :: filename

    INTEGER :: iostatus
    CHARACTER :: q

    REAL*8 :: cstrike,sstrike,cdip,sdip
    REAL*8, DIMENSION(3) :: s,d,n

    ! double-quote character
    q=char(34)

    OPEN (UNIT=15,FILE=filename,IOSTAT=iostatus,FORM="FORMATTED")
    IF (iostatus>0) THEN
       WRITE_DEBUG_INFO
       PRINT '(a)', filename
       STOP "could not open file for export in ExportVTK_Brick"
    END IF

    WRITE (15,'("<?xml version=",a,"1.0",a,"?>")') q,q
    WRITE (15,'("<VTKFile type=",a,"PolyData",a," version=",a,"0.1",a,">")') q,q,q,q
    WRITE (15,'("  <PolyData>")')

    cstrike=cos(strike)
    sstrike=sin(strike)
    cdip=cos(dip)
    sdip=sin(dip)
 
    ! strike-slip unit direction
    s(1)=sstrike
    s(2)=cstrike
    s(3)=0._8

    ! dip-slip unit direction
    d(1)=+cstrike*sdip
    d(2)=-sstrike*sdip
    d(3)=+cdip

    ! surface normal vector components
    n(1)=+cdip*cstrike
    n(2)=-cdip*sstrike
    n(3)=-sdip

    WRITE (15,'("    <Piece NumberOfPoints=",a,"8",a," NumberOfPolys=",a,"1",a,">")'),q,q,q,q
    WRITE (15,'("      <Points>")')
    WRITE (15,'("        <DataArray type=",a,"Float32",a, &
                        & " Name=",a,"Weak Zone",a, &
                        & " NumberOfComponents=",a,"3",a, &
                        & " format=",a,"ascii",a,">")'),q,q,q,q,q,q,q,q

    ! fault edge coordinates
    WRITE (15,'(24ES11.2)') &
                  x1-d(1)*W/2-s(1)*L/2-n(1)*T/2.d0, x2-d(2)*W/2-s(2)*L/2-n(2)*T/2.d0, x3-d(3)*W/2-s(3)*L/2-n(3)*T/2.d0, &
                  x1-d(1)*W/2+s(1)*L/2-n(1)*T/2.d0, x2-d(2)*W/2+s(2)*L/2-n(2)*T/2.d0, x3-d(3)*W/2+s(3)*L/2-n(3)*T/2.d0, &
                  x1+d(1)*W/2+s(1)*L/2-n(1)*T/2.d0, x2+d(2)*W/2+s(2)*L/2-n(2)*T/2.d0, x3+d(3)*W/2+s(3)*L/2-n(3)*T/2.d0, &
                  x1+d(1)*W/2-s(1)*L/2-n(1)*T/2.d0, x2+d(2)*W/2-s(2)*L/2-n(2)*T/2.d0, x3+d(3)*W/2-s(3)*L/2-n(3)*T/2.d0, &
                  x1+d(1)*W/2-s(1)*L/2+n(1)*T/2.d0, x2+d(2)*W/2-s(2)*L/2+n(2)*T/2.d0, x3+d(3)*W/2-s(3)*L/2+n(3)*T/2.d0, &
                  x1-d(1)*W/2-s(1)*L/2+n(1)*T/2.d0, x2-d(2)*W/2-s(2)*L/2+n(2)*T/2.d0, x3-d(3)*W/2-s(3)*L/2+n(3)*T/2.d0, &
                  x1-d(1)*W/2+s(1)*L/2+n(1)*T/2.d0, x2-d(2)*W/2+s(2)*L/2+n(2)*T/2.d0, x3-d(3)*W/2+s(3)*L/2+n(3)*T/2.d0, &
                  x1+d(1)*W/2+s(1)*L/2+n(1)*T/2.d0, x2+d(2)*W/2+s(2)*L/2+n(2)*T/2.d0, x3+d(3)*W/2+s(3)*L/2+n(3)*T/2.d0

    WRITE (15,'("        </DataArray>")')
    WRITE (15,'("      </Points>")')
    WRITE (15,'("      <Polys>")')
    WRITE (15,'("        <DataArray type=",a,"Int32",a, &
                         & " Name=",a,"connectivity",a, &
                         & " format=",a,"ascii",a, &
                         & " RangeMin=",a,"0",a, &
                         & " RangeMax=",a,"6",a,">")'), q,q,q,q,q,q,q,q,q,q
    WRITE (15,'("7 4 5 6 7 4 3 2 7 2 1 6")')
    WRITE (15,'("        </DataArray>")')
    WRITE (15,'("        <DataArray type=",a,"Int32",a, &
                              & " Name=",a,"offsets",a, &
                              & " format=",a,"ascii",a, &
                              & " RangeMin=",a,"12",a, &
                              & " RangeMax=",a,"12",a,">")'), q,q,q,q,q,q,q,q,q,q
    WRITE (15,'("          12")')
    WRITE (15,'("        </DataArray>")')
    WRITE (15,'("      </Polys>")')
    WRITE (15,'("    </Piece>")')

    WRITE (15,'("    <Piece NumberOfPoints=",a,"8",a," NumberOfPolys=",a,"1",a,">")'),q,q,q,q
    WRITE (15,'("      <Points>")')
    WRITE (15,'("        <DataArray type=",a,"Float32",a, &
                        & " Name=",a,"Weak Zone",a, &
                        & " NumberOfComponents=",a,"3",a, &
                        & " format=",a,"ascii",a,">")'),q,q,q,q,q,q,q,q

    ! fault edge coordinates
    WRITE (15,'(24ES11.2)') &
                  x1-d(1)*W/2.d0-s(1)*L/2.d0-n(1)*T/2.d0, x2-d(2)*W/2.d0-s(2)*L/2.d0-n(2)*T/2.d0,x3-d(3)*W/2-s(3)*L/2-n(3)*T/2.d0, &
                  x1-d(1)*W/2.d0+s(1)*L/2.d0-n(1)*T/2.d0, x2-d(2)*W/2.d0+s(2)*L/2.d0-n(2)*T/2.d0,x3-d(3)*W/2+s(3)*L/2-n(3)*T/2.d0, &
                  x1+d(1)*W/2.d0+s(1)*L/2.d0-n(1)*T/2.d0, x2+d(2)*W/2.d0+s(2)*L/2.d0-n(2)*T/2.d0,x3+d(3)*W/2+s(3)*L/2-n(3)*T/2.d0, &
                  x1+d(1)*W/2.d0-s(1)*L/2.d0-n(1)*T/2.d0, x2+d(2)*W/2.d0-s(2)*L/2.d0-n(2)*T/2.d0,x3+d(3)*W/2-s(3)*L/2-n(3)*T/2.d0, &
                  x1+d(1)*W/2.d0-s(1)*L/2.d0+n(1)*T/2.d0, x2+d(2)*W/2.d0-s(2)*L/2.d0+n(2)*T/2.d0,x3+d(3)*W/2-s(3)*L/2+n(3)*T/2.d0, &
                  x1-d(1)*W/2.d0-s(1)*L/2.d0+n(1)*T/2.d0, x2-d(2)*W/2.d0-s(2)*L/2.d0+n(2)*T/2.d0,x3-d(3)*W/2-s(3)*L/2+n(3)*T/2.d0, &
                  x1-d(1)*W/2.d0+s(1)*L/2.d0+n(1)*T/2.d0, x2-d(2)*W/2.d0+s(2)*L/2.d0+n(2)*T/2.d0,x3-d(3)*W/2+s(3)*L/2+n(3)*T/2.d0, &
                  x1+d(1)*W/2.d0+s(1)*L/2.d0+n(1)*T/2.d0, x2+d(2)*W/2.d0+s(2)*L/2.d0+n(2)*T/2.d0,x3+d(3)*W/2+s(3)*L/2+n(3)*T/2.d0

    WRITE (15,'("        </DataArray>")')
    WRITE (15,'("      </Points>")')
    WRITE (15,'("      <Polys>")')
    WRITE (15,'("        <DataArray type=",a,"Int32",a, &
                         & " Name=",a,"connectivity",a, &
                         & " format=",a,"ascii",a, &
                         & " RangeMin=",a,"0",a, &
                         & " RangeMax=",a,"7",a,">")'), q,q,q,q,q,q,q,q,q,q
    WRITE (15,'("0 1 2 3 0 5 4 3 0 1 6 5")')
    WRITE (15,'("        </DataArray>")')
    WRITE (15,'("        <DataArray type=",a,"Int32",a, &
                              & " Name=",a,"offsets",a, &
                              & " format=",a,"ascii",a, &
                              & " RangeMin=",a,"12",a, &
                              & " RangeMax=",a,"12",a,">")'), q,q,q,q,q,q,q,q,q,q
    WRITE (15,'("          12")')
    WRITE (15,'("        </DataArray>")')
    WRITE (15,'("      </Polys>")')
    WRITE (15,'("    </Piece>")')
    WRITE (15,'("  </PolyData>")')
    WRITE (15,'("</VTKFile>")')

    CLOSE(15)

  END SUBROUTINE exportvtk_brick

  !------------------------------------------------------------------
  !> subroutine ExportVTK_Vectors
  !! creates a .vtr file (in the VTK Rectilinear XML format) 
  !! containing a vector field.
  !!
  !! \author sylvain barbot 06/25/09 - original form
  !------------------------------------------------------------------
  SUBROUTINE exportvtk_vectors(u1,u2,u3,sx1,sx2,sx3,dx1,dx2,dx3,j1,j2,j3,vcfilename)
    INTEGER, INTENT(IN) :: sx1,sx2,sx3,j1,j2,j3
#ifdef ALIGN_DATA
    REAL*4, INTENT(IN), DIMENSION(sx1+2,sx2,sx3) :: u1,u2,u3
#else
    REAL*4, INTENT(IN), DIMENSION(sx1,sx2,sx3) :: u1,u2,u3
#endif
    REAL*8, INTENT(IN) :: dx1,dx2,dx3
    CHARACTER(256), INTENT(IN) :: vcfilename

    INTEGER :: iostatus,idum,i1,i2
    CHARACTER :: q
    INTEGER :: k1,k2,k3
    REAL*8 :: x1,x2,x3

    ! double-quote character
    q=char(34)

    OPEN (UNIT=15,FILE=vcfilename,IOSTAT=iostatus,FORM="FORMATTED")
    IF (iostatus>0) THEN
       WRITE_DEBUG_INFO
       PRINT '(a)', vcfilename
       STOP "could not open file for export"
    END IF

    WRITE (15,'("<?xml version=",a,"0.1",a,"?>")') q,q
    WRITE (15,'("<VTKFile type=",a,"RectilinearGrid",a," version=",a,"0.1",a,">")') q,q,q,q
    WRITE (15,'("  <RectilinearGrid WholeExtent=",a,6I5.4,a,">")') q,1,sx1/j1,1,sx2/j2,1,sx3/j3,q
    WRITE (15,'("  <Piece Extent=",a,6I5.4,a,">")') q,1,sx1/j1,1,sx2/j2,1,sx3/j3,q
    WRITE (15,'("    <PointData Scalars=",a,"Vector Field",a,">")') q,q

    WRITE (15,'("    <DataArray type=",a,"Float32",a, &
                       & " Name=",a,"X Velocity",a, &
                       & " format=",a,"ascii",a,">")') q,q,q,q,q,q

    ! write first component values
    DO k3=0,sx3-1,j3
       x3=REAL(k3,8)
       DO k2=-sx2/2,sx2/2-1,j2
          x2=REAL(k2,8)
          DO k1=-sx1/2,sx1/2-1,j1
             x1=REAL(k1,8)

             CALL shiftedindex(x1,x2,1._8,sx1,sx2,sx3,1._8,1._8,1._8,i1,i2,idum)
             WRITE (15,'(ES12.2)') u1(i1,i2,k3+1)
          END DO
       END DO
    END DO
    WRITE (15,'("    </DataArray>")')

    WRITE (15,'("    <DataArray type=",a,"Float32",a, &
                       & " Name=",a,"Y Velocity",a, &
                       & " format=",a,"ascii",a,">")') q,q,q,q,q,q

    ! write second component values
    DO k3=0,sx3-1,j3
       x3=REAL(k3,8)
       DO k2=-sx2/2,sx2/2-1,j2
          x2=REAL(k2,8)
          DO k1=-sx1/2,sx1/2-1,j1
             x1=REAL(k1,8)

             CALL shiftedindex(x1,x2,1._8,sx1,sx2,sx3,1._8,1._8,1._8,i1,i2,idum)
             WRITE (15,'(ES12.2)') u2(i1,i2,k3+1)

          END DO
       END DO
    END DO
    WRITE (15,'("    </DataArray>")')

    WRITE (15,'("    <DataArray type=",a,"Float32",a, &
                       & " Name=",a,"Z Velocity",a, &
                       & " format=",a,"ascii",a,">")') q,q,q,q,q,q

    ! write third component values
    DO k3=0,sx3-1,j3
       x3=REAL(k3,8)
       DO k2=-sx2/2,sx2/2-1,j2
          x2=REAL(k2,8)
          DO k1=-sx1/2,sx1/2-1,j1
             x1=REAL(k1,8)

             CALL shiftedindex(x1,x2,1._8,sx1,sx2,sx3,1._8,1._8,1._8,i1,i2,idum)
             WRITE (15,'(ES12.2)') u3(i1,i2,k3+1)

          END DO
       END DO
    END DO
    WRITE (15,'("    </DataArray>")')

    WRITE (15,'("  </PointData>")')

    WRITE (15,'("  <Coordinates>")')

    WRITE (15,'("    <DataArray type=",a,"Float32",a, &
                       & " Name=",a,"Array 1",a, &
                       & " format=",a,"ascii",a, &
                       & " RangeMin=",a,ES12.2,a, &
                       & " RangeMax=",a,ES12.2,a,">")') q,q,q,q,q,q,q,-sx1/2*dx1,q,q,(sx1/2-1)*dx1,q
    DO k1=-sx1/2,sx1/2-1,j1
       x1=REAL(k1,8)
       WRITE (15,'(ES12.2)') x1*dx1
    END DO
    WRITE (15,'("    </DataArray>")')
    WRITE (15,'("    <DataArray type=",a,"Float32",a, &
                       & " Name=",a,"Array 2",a, &
                       & " format=",a,"ascii",a, &
                       & " RangeMin=",a,ES12.2,a, &
                       & " RangeMax=",a,ES12.2,a,">")') q,q,q,q,q,q,q,-sx2/2*dx2,q,q,(sx2/2-1)*dx2,q
    DO k2=-sx2/2,sx2/2-1,j2
       x2=REAL(k2,8)
       WRITE (15,'(ES12.2)') x2*dx2
    END DO
    WRITE (15,'("    </DataArray>")')
    WRITE (15,'("    <DataArray type=",a,"Float32",a, &
                       & " Name=",a,"Array 3",a, &
                       & " format=",a,"ascii",a, &
                       & " RangeMin=",a,ES12.2,a, &
                       & " RangeMax=",a,ES12.2,a,">")') q,q,q,q,q,q,q,0,q,q,(sx3-1)*dx3,q
    DO k3=0,sx3-1,j3
       x3=REAL(k3,8)
       WRITE (15,'(ES12.2)') x3*dx3
    END DO
    WRITE (15,'("    </DataArray>")')

    WRITE (15,'("  </Coordinates>")')
    WRITE (15,'("</Piece>")')
    WRITE (15,'("</RectilinearGrid>")')
    WRITE (15,'("</VTKFile>")')

    CLOSE(15)

  END SUBROUTINE exportvtk_vectors

  !------------------------------------------------------------------
  !> subroutine ExportVTK_Vectors_Slice
  !! creates a .vtr file (in the VTK Rectilinear XML format) 
  !! containing a vector field.
  !!
  !! \author sylvain barbot 06/25/09 - original form
  !------------------------------------------------------------------
  SUBROUTINE exportvtk_vectors_slice(u1,u2,u3,sx1,sx2,sx3,dx1,dx2,dx3,oz,j1,j2,vcfilename)
    INTEGER, INTENT(IN) :: sx1,sx2,sx3,j1,j2
#ifdef ALIGN_DATA
    REAL*4, INTENT(IN), DIMENSION(sx1+2,sx2,sx3) :: u1,u2,u3
#else
    REAL*4, INTENT(IN), DIMENSION(sx1,sx2,sx3) :: u1,u2,u3
#endif
    REAL*8, INTENT(IN) :: dx1,dx2,dx3,oz
    CHARACTER(256), INTENT(IN) :: vcfilename

    INTEGER :: iostatus,idum,i1,i2
    CHARACTER :: q
    INTEGER :: k1,k2,k3
    REAL*8 :: x1,x2,x3

    ! double-quote character
    q=char(34)

    OPEN (UNIT=15,FILE=vcfilename,IOSTAT=iostatus,FORM="FORMATTED")
    IF (iostatus>0) THEN
       WRITE_DEBUG_INFO
       PRINT '(a)', vcfilename
       STOP "could not open file for export"
    END IF

    WRITE (15,'("<?xml version=",a,"0.1",a,"?>")') q,q
    WRITE (15,'("<VTKFile type=",a,"RectilinearGrid",a," version=",a,"0.1",a,">")') q,q,q,q
    WRITE (15,'("  <RectilinearGrid WholeExtent=",a,6I5.4,a,">")') q,1,sx1/j1,1,sx2/j2,1,1,q
    WRITE (15,'("  <Piece Extent=",a,6I5.4,a,">")') q,1,sx1/j1,1,sx2/j2,1,1,q
    WRITE (15,'("    <PointData Scalars=",a,"Vector Field",a,">")') q,q

    WRITE (15,'("    <DataArray type=",a,"Float32",a, &
                       & " Name=",a,"X Velocity",a, &
                       & " format=",a,"ascii",a,">")') q,q,q,q,q,q

    ! write first component values
    x3=oz/dx3
    DO k2=-sx2/2,sx2/2-1,j2
       x2=REAL(k2,8)
       DO k1=-sx1/2,sx1/2-1,j1
          x1=REAL(k1,8)

          CALL shiftedindex(x1,x2,1._8,sx1,sx2,sx3,1._8,1._8,1._8,i1,i2,idum)
          WRITE (15,'(ES12.2)') u1(i1,i2,k3+1)
       END DO
    END DO
    WRITE (15,'("    </DataArray>")')

    WRITE (15,'("    <DataArray type=",a,"Float32",a, &
                       & " Name=",a,"Y Velocity",a, &
                       & " format=",a,"ascii",a,">")') q,q,q,q,q,q

    ! write second component values
    x3=oz/dx3
    DO k2=-sx2/2,sx2/2-1,j2
       x2=REAL(k2,8)
       DO k1=-sx1/2,sx1/2-1,j1
          x1=REAL(k1,8)

          CALL shiftedindex(x1,x2,1._8,sx1,sx2,sx3,1._8,1._8,1._8,i1,i2,idum)
          WRITE (15,'(ES12.2)') u2(i1,i2,k3+1)

       END DO
    END DO
    WRITE (15,'("    </DataArray>")')

    WRITE (15,'("    <DataArray type=",a,"Float32",a, &
                       & " Name=",a,"Z Velocity",a, &
                       & " format=",a,"ascii",a,">")') q,q,q,q,q,q

    ! write third component values
    x3=oz/dx3
    DO k2=-sx2/2,sx2/2-1,j2
       x2=REAL(k2,8)
       DO k1=-sx1/2,sx1/2-1,j1
          x1=REAL(k1,8)

          CALL shiftedindex(x1,x2,1._8,sx1,sx2,sx3,1._8,1._8,1._8,i1,i2,idum)
          WRITE (15,'(ES12.2)') u3(i1,i2,k3+1)

       END DO
    END DO
    WRITE (15,'("    </DataArray>")')

    WRITE (15,'("  </PointData>")')

    WRITE (15,'("  <Coordinates>")')

    WRITE (15,'("    <DataArray type=",a,"Float32",a, &
                       & " Name=",a,"Array 1",a, &
                       & " format=",a,"ascii",a, &
                       & " RangeMin=",a,ES12.2,a, &
                       & " RangeMax=",a,ES12.2,a,">")') q,q,q,q,q,q,q,-sx1/2*dx1,q,q,(sx1/2-1)*dx1,q
    DO k1=-sx1/2,sx1/2-1,j1
       x1=REAL(k1,8)
       WRITE (15,'(ES12.2)') x1*dx1
    END DO
    WRITE (15,'("    </DataArray>")')
    WRITE (15,'("    <DataArray type=",a,"Float32",a, &
                       & " Name=",a,"Array 2",a, &
                       & " format=",a,"ascii",a, &
                       & " RangeMin=",a,ES12.2,a, &
                       & " RangeMax=",a,ES12.2,a,">")') q,q,q,q,q,q,q,-sx2/2*dx1,q,q,(sx2/2-1)*dx2,q
    DO k2=-sx2/2,sx2/2-1,j2
       x2=REAL(k2,8)
       WRITE (15,'(ES12.2)') x2*dx2
    END DO
    WRITE (15,'("    </DataArray>")')
    WRITE (15,'("    <DataArray type=",a,"Float32",a, &
                       & " Name=",a,"Array 3",a, &
                       & " format=",a,"ascii",a, &
                       & " RangeMin=",a,ES12.2,a, &
                       & " RangeMax=",a,ES12.2,a,">")') q,q,q,q,q,q,q,oz,q,q,oz,q
    WRITE (15,'(2ES12.2)') oz
    WRITE (15,'("    </DataArray>")')

    WRITE (15,'("  </Coordinates>")')
    WRITE (15,'("</Piece>")')
    WRITE (15,'("</RectilinearGrid>")')
    WRITE (15,'("</VTKFile>")')

    CLOSE(15)

  END SUBROUTINE exportvtk_vectors_slice
#endif

END MODULE export
