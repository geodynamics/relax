!-----------------------------------------------------------------------
! Copyright 2021 Sylvain Barbot
!
! This file is part of UNICYCLE
!
! UNICYCLE is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! UNICYCLE is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with UNICYCLE.  If not, see <http://www.gnu.org/licenses/>.
!-----------------------------------------------------------------------

MODULE exportnetcdf

  USE netcdf

  IMPLICIT NONE

  PUBLIC

  INTEGER, PRIVATE, PARAMETER :: NDIMS = 2

  CHARACTER (LEN = *), PRIVATE, PARAMETER :: ACTUAL_RANGE = "actual_range"

CONTAINS

  SUBROUTINE writegrd(z,nx,ny,y0,x0,dy,dx,land,dummy,title,filename)
    REAL*4, INTENT(IN) :: z(nx,ny)
    INTEGER, INTENT(IN) :: nx,ny
    REAL*8, INTENT(IN) :: y0,x0
    REAL*8, INTENT(IN) :: dy,dx,land,dummy
    CHARACTER(LEN=256), INTENT(IN) :: title
    CHARACTER(LEN=256), INTENT(IN) :: filename

    INTEGER, PARAMETER :: registration=1
    INTEGER :: i

    REAL*8, DIMENSION(nx) :: x
    REAL*8, DIMENSION(ny) :: y

    DO i=0,nx-1
       x(i+1)=x0+REAL(i,8)*dx
    END DO

    DO i=0,ny-1
       y(i+1)=y0+REAL(i,8)*dy
    END DO

    CALL writeNetcdf(filename,nx,x,ny,y,z,title,registration)

  END SUBROUTINE writegrd

  !------------------------------------------------------------------------------
  ! subroutine writeNetcdf
  ! writes a netcdf file compatible with the GMT .grd format
  ! use 0 for gridline registration and 1 for pixel node registration
  !------------------------------------------------------------------------------
  SUBROUTINE writeNetcdf(filename,nx,x,ny,y,z,titlename,registration)
    CHARACTER(LEN=256), INTENT(IN) :: filename
    INTEGER, INTENT(IN) :: nx,ny
    REAL*8, INTENT(IN) :: x(nx),y(ny)
    REAL*4, INTENT(IN) :: z(nx,ny)
    CHARACTER(LEN=256), INTENT(IN) :: titlename
    INTEGER, INTENT(IN) :: registration

    CHARACTER (LEN = *), PARAMETER :: X_NAME = "x"
    CHARACTER (LEN = *), PARAMETER :: Y_NAME = "y"
    CHARACTER (LEN = *), PARAMETER :: Z_NAME = "z"

    CHARACTER (LEN = *), PARAMETER :: LONG_NAME = "long_name"
    CHARACTER (LEN = *), PARAMETER :: NODE_OFFSET = "node_offset"
    CHARACTER (LEN = *), PARAMETER :: TITLE = "title"
    CHARACTER (LEN = *), PARAMETER :: HISTORY = "history"
    CHARACTER (LEN = *), PARAMETER :: DESCRIPTION = "Relax output"
    CHARACTER (LEN = *), PARAMETER :: CONVENTIONS = "Conventions"
    CHARACTER (LEN = *), PARAMETER :: COARDS = "COARDS/CF-1.0"

    INTEGER :: ncid
    INTEGER :: x_dimid, y_dimid
    INTEGER :: x_varid, y_varid, z_varid
    INTEGER :: dimids(NDIMS)

    ! create the netcdf file
    CALL check(nf90_create(filename,nf90_clobber,ncid))

    ! dimensions
    CALL check(nf90_def_dim(ncid,X_NAME,nx,x_dimid))
    CALL check(nf90_def_dim(ncid,Y_NAME,ny,y_dimid))

    ! define the coordinate variables
    CALL check(nf90_def_var(ncid,X_NAME,NF90_DOUBLE,x_dimid,x_varid))
    CALL check(nf90_def_var(ncid,Y_NAME,NF90_DOUBLE,y_dimid,y_varid))

    dimids = (/ x_dimid, y_dimid /)
  
    ! prepare data variable
    CALL check(nf90_def_var(ncid,Z_NAME,NF90_REAL,dimids,z_varid))

    ! attributes
    CALL check(nf90_put_att(ncid,x_varid,LONG_NAME,X_NAME))
    CALL check(nf90_put_att(ncid,y_varid,LONG_NAME,Y_NAME))
    CALL check(nf90_put_att(ncid,z_varid,LONG_NAME,Z_NAME))

    IF (0 .EQ. registration) THEN
       CALL check(nf90_put_att(ncid,x_varid,ACTUAL_RANGE,(/ x(1), x(nx) /)))
       CALL check(nf90_put_att(ncid,y_varid,ACTUAL_RANGE,(/ y(1), y(ny) /)))
    ELSE
       CALL check(nf90_put_att(ncid,x_varid,ACTUAL_RANGE,(/ x(1)-0.5d0*(x(2)-x(1)), x(nx)+0.5d0*(x(nx)-x(nx-1)) /)))
       CALL check(nf90_put_att(ncid,y_varid,ACTUAL_RANGE,(/ y(1)-0.5d0*(y(2)-y(1)), y(ny)+0.5d0*(y(ny)-y(ny-1)) /)))
       CALL check(nf90_put_att(ncid,NF90_GLOBAL,NODE_OFFSET,1))
    END IF
    CALL check(nf90_put_att(ncid,z_varid,ACTUAL_RANGE,(/ MINVAL(z), MAXVAL(z) /)))

    ! global attributes
    CALL check(nf90_put_att(ncid,NF90_GLOBAL,TITLE,TRIM(titlename)))
    CALL check(nf90_put_att(ncid,NF90_GLOBAL,HISTORY,DESCRIPTION))
    CALL check(nf90_put_att(ncid,NF90_GLOBAL,CONVENTIONS,COARDS))

    ! end define mode
    CALL check(nf90_enddef(ncid))
 
    ! write coordinate system
    CALL check(nf90_put_var(ncid,x_varid,x))
    CALL check(nf90_put_var(ncid,y_varid,y))

    ! write the data
    CALL check(nf90_put_var(ncid,z_varid,z))

    ! flush
    CALL check(nf90_sync(ncid))

    ! close file
    CALL check(nf90_close(ncid))

  END SUBROUTINE writeNetcdf

  !------------------------------------------------------------------------------
  ! subroutine check
  ! checks the return of nf90_ functions
  !------------------------------------------------------------------------------
  SUBROUTINE check(status)
    INTEGER, INTENT(IN) :: status
          
    IF (status /= nf90_noerr) THEN 
       PRINT *, TRIM(nf90_strerror(status))
       STOP "netcdf write stopped"
    END IF
  END SUBROUTINE check

END MODULE exportnetcdf

