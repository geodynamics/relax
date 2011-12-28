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

MODULE types

  TYPE SOURCE_STRUCT
     SEQUENCE
     REAL*8 :: slip,x,y,z,width,length,strike,dip,rake,period,phase,beta
  END TYPE SOURCE_STRUCT

  TYPE PLANE_STRUCT
     SEQUENCE
     REAL*8 :: x,y,z,width,length,strike,dip,rake
  END TYPE PLANE_STRUCT

  TYPE LAYER_STRUCT
     SEQUENCE
     REAL*8 :: z,gammadot0,stressexponent,cohesion,friction
  END TYPE LAYER_STRUCT

  TYPE WEAK_STRUCT
     SEQUENCE
     REAL*8 :: dgammadot0,x,y,z,width,length,thickness,strike,dip
  END TYPE WEAK_STRUCT

  TYPE VECTOR_STRUCT
     SEQUENCE
     REAL*8 :: v1,v2,v3
  END TYPE VECTOR_STRUCT

  TYPE TENSOR
     SEQUENCE
     REAL*4 :: s11,s12,s13,s22,s23,s33
  END TYPE TENSOR

  TYPE TENSOR_LAYER_STRUCT
     SEQUENCE
     REAL*4 :: z,dum
     TYPE(TENSOR) :: t
  END TYPE TENSOR_LAYER_STRUCT

  TYPE SEGMENT_STRUCT
     SEQUENCE
     REAL*8 :: x,y,z,width,length,strike,dip,friction
     TYPE(TENSOR) :: sig0
  END TYPE SEGMENT_STRUCT

  TYPE SLIPPATCH_STRUCT
     SEQUENCE
     REAL*8 :: x1,x2,x3,lx,lz,slip,ss,ds
  END TYPE SLIPPATCH_STRUCT

  TYPE EVENT_STRUC
     REAL*8 :: time
     INTEGER*4 :: i,ns,nt,nm,nl
     TYPE(SOURCE_STRUCT), DIMENSION(:), ALLOCATABLE :: s,sc,ts,tsc,m,mc,l,lc
  END TYPE EVENT_STRUC
  
  TYPE, PUBLIC :: SIMULATION_STRUC
     ! grid dimension
     INTEGER :: sx1,sx2,sx3

     ! sampling
     REAL*8 :: dx1,dx2,dx3

     ! smoothing factor
     REAL*8 :: beta

     ! filter parameter for slip models
     REAL*8 :: nyquist

     ! center coordinates and rotation
     REAL*8 :: x0,y0,rot

#ifdef PROJ
     ! geographic coordinates of center, UTM zone, length unit
     REAL*8 :: lon0,lat0,umult
     INTEGER :: zone
#endif

     ! observation depths
     REAL*8 :: oz,ozs

     ! output directory
     CHARACTER(80) :: wdir

     ! filenames
     CHARACTER(80) :: reportfilename,reporttimefilename

     ! elastic moduli and gravity parameter
     REAL*8 :: lambda,mu,gam

     ! time step parameters
     REAL*8 :: interval
     REAL*8 :: odt,tscale
     INTEGER :: skip=0

     ! number of observation planes
     INTEGER :: nop

     ! observation planes
     TYPE(PLANE_STRUCT), DIMENSION(:), ALLOCATABLE :: op

     ! number of stress observation planes
     INTEGER :: nsop

     ! stress observation planes
     TYPE(SEGMENT_STRUCT), DIMENSION(:), ALLOCATABLE :: sop

     ! number of observation points
     INTEGER :: npts

     ! observation points
     TYPE(VECTOR_STRUCT), DIMENSION(:), ALLOCATABLE :: opts

     ! observation points name
     CHARACTER(LEN=4), DIMENSION(:), ALLOCATABLE :: ptsname

     ! number of prestress interfaces
     INTEGER :: nps

     ! stress layers and stress structure
     TYPE(TENSOR_LAYER_STRUCT), DIMENSION(:), ALLOCATABLE :: stresslayer,stressstruc

     ! number of linear viscous interfaces
     INTEGER :: nv

     ! linear viscous layers and structure
     TYPE(LAYER_STRUCT), DIMENSION(:), ALLOCATABLE :: linearlayer,linearstruc

     ! number of linear weak zones
     INTEGER :: nlwz

     ! linear weak zones
     TYPE(WEAK_STRUCT), DIMENSION(:), ALLOCATABLE :: linearweakzone,linearweakzonec

     ! number of nonlinear viscous interfaces
     INTEGER :: npl

     ! nonlinear viscous layers and structure
     TYPE(LAYER_STRUCT), DIMENSION(:), ALLOCATABLE :: nonlinearlayer,nonlinearstruc

     ! number of nonlinear weak zones
     INTEGER :: nnlwz

     ! nonlinear viscous layers and structure
     TYPE(WEAK_STRUCT), DIMENSION(:), ALLOCATABLE :: nonlinearweakzone,nonlinearweakzonec

     ! number of fault creep interfaces
     INTEGER :: nfc

     ! fault creep interfaces
     TYPE(LAYER_STRUCT), DIMENSION(:), ALLOCATABLE :: faultcreeplayer,faultcreepstruc

     ! number of afterslip planes
     INTEGER :: np

     ! afterslip planes
     TYPE(PLANE_STRUCT), DIMENSION(:), ALLOCATABLE :: n

     ! interseismic event
     TYPE(EVENT_STRUC) :: inter

     ! number of coseismic events
     INTEGER :: ne

     ! coseismic events
     TYPE(EVENT_STRUC), DIMENSION(:), ALLOCATABLE :: events

     ! overrides output to formats
     LOGICAL :: isoutputproj=.TRUE.
     LOGICAL :: isoutputrelax=.TRUE.
     LOGICAL :: isoutputtxt=.TRUE.
     LOGICAL :: isoutputvtk=.TRUE.
     LOGICAL :: isoutputvtkrelax=.TRUE.
     LOGICAL :: isoutputgrd=.TRUE.
     LOGICAL :: isoutputxyz=.TRUE.
     LOGICAL :: isoutputstress=.TRUE.

     ! other options
     LOGICAL :: isdryrun=.FALSE.
     LOGICAL :: ishelp=.FALSE.

  END TYPE SIMULATION_STRUC

END MODULE types
