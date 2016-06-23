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
    
   IMPLICIT NONE

CONTAINS 

   SUBROUTINE ispresent(var,avail)
    REAL*4, OPTIONAL :: var 
    INTEGER :: avail

    IF (PRESENT(var)) avail=1

   END SUBROUTINE

   SUBROUTINE isallocated(var,avail)
    REAL*4, OPTIONAL, ALLOCATABLE :: var 
    INTEGER :: avail

    IF (ALLOCATED(var)) avail=1

   END SUBROUTINE
END MODULE

