!-----------------------------------------------------------------------
! Copyright 2007-2012, Sylvain Barbot
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

  !-----------------------------------------------------------------------
  !> \mainpage 
  !! program relax
  !! <hr>
  !! PURPOSE:
  !!   The program RELAX computes nonlinear time-dependent viscoelastic
  !!   deformation with powerlaw rheology and rate-strengthening friction 
  !!   in a cubic, periodic grid due to coseismic stress changes, initial
  !!   stress, surface loads, and/or moving faults.
  !! 
  !! ONLINE DOCUMENTATION:
  !!   generate html documentation from the source directory with the 
  !!   doxygen (http://www.stack.nl/~dimitri/doxygen/index.html) 
  !!   program with command:
  !!
  !!     doxygen .doxygen
  !!
  !! DESCRIPTION:
  !!   Computation is done semi-analytically inside a cartesian grid.
  !!   The grid is defined by its size sx1*sx2*sx3 and the sampling
  !!   intervals dx1, dx2 and dx3. rule of thumb is to allow for at least
  !!   five samples per fault length or width, and to have the tip of any 
  !!   fault at least 10 fault widths away from any edge of the 
  !!   computational grid.
  !!
  !!   Coseismic stress changes and initial coseismic deformation results
  !!   from the presence of dislocations in the brittle layer. Fault
  !!   geometry is prescribed following Okada or Wang's convention, with the
  !!   usual slip, strike, dip and rake and is converted to a double-couple
  !!   equivalent body-force analytically. Current implementation allows 
  !!   shear fault (strike slip and dip slip), dykes, Mogi source, and
  !!   surface traction. Faults and dykes can be of arbitrary orientation 
  !!   in the half space.
  !!
  !! <hr>
  !!
  !! METHOD:
  !!   The current implementation is organized to integrate stress/strain-
  !!   rate constitutive laws (rheologies) of the form
  !! \f[
  !!       \dot{\epsilon} = f(\sigma)
  !! \f]
  !!   as opposed to epsilon^dot = f(sigma,epsilon) wich would include work-
  !!   hardening (or weakening). The time-stepping implements a second-order
  !!   Runge-Kutta numerical integration scheme with a variable time-step.
  !!   The Runge-Kutta method integrating the ODE y'=f(x,y) can be summarized
  !!   as follows:
  !! \f[
  !!          y_(n+1) = y_n + k_2
  !!              k_1 = h * f(x_n, y_n)
  !!              k_2 = h * f(x_n + h, y_n + k_1)
  !! \f]
  !!   where h is the time-step and n is the time-index. The elastic response
  !!   in the computational grid is obtained using elastic Greens functions.
  !!   The Greens functions are applied in the Fourier domain. Strain,
  !!   stress and body-forces are obtained by application of a finite impulse
  !!   response (FIR) differentiator filter in the space domain.
  !!
  !! <hr>
  !!
  !! INPUT:
  !!   Static dislocation sources are discretized into a series of planar
  !!   segments. Slip patches are defined in terms of position, orientation,
  !!   and slip, as illustrated in the following figure:
  !!\verbatim
  !!                 N (x1)
  !!                /
  !!               /| Strike
  !!   x1,x2,x3 ->@------------------------      (x2)
  !!              |\        p .            \ W
  !!              :-\      i .              \ i
  !!              |  \    l .                \ d
  !!              :90 \  S .                  \ t
  !!              |-Dip\  .                    \ h
  !!              :     \. | Rake               \
  !!              |      -------------------------
  !!              :             L e n g t h
  !!              Z (x3)
  !!\endverbatim
  !!   Dislocations are converted to double-couple equivalent body-force
  !!   analytically. Solution displacement is obtained by application of
  !!   the Greens functions in the Fourier domain.
  !!
  !!   For friction faults where slip rates are evaluated from stress and
  !!   a constitutive law, the rake corresponds to the orientation of slip. 
  !!   That is, if r_i is the rake vector and v_i is the instantaneous 
  !!   velocity vector, then r_j v_j >= 0. 
  !!
  !! <hr>
  !!
  !! OUTPUT:
  !!   The vector-valued deformation is computed everywhere in a cartesian
  !!   grid. The vector field is sampled 1) along a horizontal surface at a
  !!   specified depth and 2) at specific points. Format is always North (x1), 
  !!   East (x2) and Down (x3) components, following the right-handed reference 
  !!   system convention. North corresponds to x1-direction, East to the 
  !!   x2-direction and down to the x3-direction. The Generic Mapping Tool 
  !!   output files are labeled explicitely ???-north.grd, ???-east.grd and 
  !!   ???-up.grd (or say, ???-geo-up.grd for outputs in geographic 
  !!   coordinates), where ??? stands for an output index: 001, 002, ...
  !!
  !!   The amplitude of the inelastic (irreversible) deformation is also
  !!   tracked and can be output along a plane of arbitrary orientation.
  !!   The inelastic deformation includes the initial, constrained, slip on
  !!   fault surfaces, the time-dependent slip on frictional surfaces and
  !!   the cumulative amplitude of bulk strain in viscoelastic regions.
  !!   Slip is provided as a function of local coordinates along strike and 
  !!   dip as well as a function of the Cartesian coordinates for three-
  !!   dimensional display.
  !!
  !!   Time integration uses adaptive time steps to ensure accuracy but
  !!   results can be output either 1) at specified uniform time intervals 
  !!   or 2) at the same intervals as computed. In the later case, output 
  !!   intervals is chosen internally depending on instantaneous relaxation 
  !!   rates.
  !!
  !! <hr>
  !!
  !! TECHNICAL ASPECTS:
  !!   Most of the computational burden comes from 1) applying the elastic
  !!   Green function and 2) computing the current strain from a displacement
  !!   field. The convolution of body forces with the Green function is 
  !!   performed in the Fourier domain and the efficiency of the computation
  !!   depends essentially upon a choice of the discrete Fourier transform.
  !!   Current implementation is compatible with the Couley-Tuckey, the
  !!   Fast Fourier transform of the West (FFTW), the SGI FFT and the intel
  !!   FFT from the intel MKL library. Among these choices, the MKL FFT is
  !!   the most efficient. The FFTW, SGI FFT and MKL FFT can all be ran
  !!   in parallel on shared-memory computers.
  !!
  !!   Strain is computed using a Finite Impulse Response differentiator
  !!   filter in the space domain. Use of FIR filter give rise to very
  !!   accurate derivatives but is computationally expensive. The filter
  !!   kernels are provided in the kernel???.inc files. Use of a compact
  !!   kernel may accelerate computation significantly.
  !!
  !!   Compilation options are defined in the include.f90 file and specify
  !!   for instance the choice of DFT and the kind of output provided.
  !!
  !! MODIFICATIONS:
  !! \author Sylvain Barbot 
  !! (07-06-07) - original form                                    <br>
  !! (08-28-08) - FFTW/SGI_FFT support, FIR derivatives,
  !!              Runge-Kutta integration, tensile cracks,
  !!              GMT output, comments in input file               <br>
  !! (10-24-08) - interseismic loading, postseismic signal
  !!              output in separate files                         <br>
  !! (12-08-09) - slip distribution smoothing                      <br>
  !! (05-05-10) - lateral variations in viscous properties
  !!              Intel MKL implementation of the FFT              <br>
  !! (06-04-10) - output in geographic coordinates
  !!              and output components of stress tensor           <br>
  !! (07-19-10) - includes surface tractions initial condition
  !!              output geometry in VTK format for Paraview       <br>
  !! (02-28-11) - add constraints on the broad direction of 
  !!              afterslip, export faults to GMT xy format
  !!              and allow scaling of computed time steps.        <br>
  !! (04-26-11) - include command-line arguments
  !! (11-04-11) - compatible with gfortran                         <br>
  !!
  !! \todo 
  !!   - homogenize VTK output so that geometry of events match event index
  !!   - evaluate Green's function, stress and body forces in GPU
  !!   - write the code for MPI multi-thread
  !!   - fix the vtk export to grid for anisotropic sampling
  !!   - export position of observation points to long/lat in opts-geo.dat
  !!   - check the projected output on the south hemisphere
  !!   - check the fully-relaxed afterslip for uniform stress change
  !!   - include topography of parameter interface
  !!   - export afterslip output in VTK legacy format (binary)
  !------------------------------------------------------------------------
PROGRAM relax

  USE input
  USE green
  USE green_space
  USE elastic3d
  USE viscoelastic3d
  USE friction3d
  USE export

#include "include.f90"
  
  IMPLICIT NONE

  INTEGER, PARAMETER :: ITERATION_MAX = 99999 
  REAL*8, PARAMETER :: STEP_MAX = 1e7

  INTEGER :: i,k,e,oi,iostatus,mech(3),iTensor,iRet
  
  REAL*8 :: maxwell(3)
  TYPE(SIMULATION_STRUC) :: in
#ifdef VTK
  CHARACTER(256) :: filename,title,name
  CHARACTER(3) :: digit
#endif
  CHARACTER(4) :: digit4
  REAL*8 :: t,Dt,tm,dAmp,dAmp1
  REAL*4 :: cuC1, cuC2 
 
  ! arrays
  REAL*4, DIMENSION(:,:,:), ALLOCATABLE :: v1,v2,v3,u1,u2,u3,gamma
  REAL*4, DIMENSION(:,:,:), ALLOCATABLE :: u1r,u2r,u3r
  REAL*4, DIMENSION(:,:), ALLOCATABLE :: t1,t2,t3
  REAL*4, DIMENSION(:,:,:), ALLOCATABLE :: inter1,inter2,inter3
  TYPE(TENSOR), DIMENSION(:,:,:), ALLOCATABLE :: tau,sig,moment
 
#ifdef PAPI_PROF
    CHARACTER (LEN=16) cTimerName
    cTimerName = 'relax'
#endif

  ! read input parameters
  CALL init(in)

  ! abort calculation after help message
  ! or for dry runs
  IF (in%isdryrun) THEN
     PRINT '("dry run: abort calculation")'
  END IF
  IF (in%isdryrun .OR. in%ishelp) THEN
     ! exit program
     GOTO 100
  END IF

  CALL cuinit (%VAL(in%sx1), %VAL(in%sx2), %VAL(in%sx3), %VAL(in%dx1), %VAL(in%dx2), %VAL(in%dx3), iostatus)
  IF (iostatus>0) STOP "could not allocate memory"

  ! allocate memory
  ALLOCATE (v1(1,1,1),v2(1,1,1),v3(1,1,1), &
            u1(1,1,1),u2(1,1,1),u3(1,1,1), &
            tau(in%sx1,in%sx2,in%sx3/2),sig(1,1,1),gamma(1,1,1), &
            t1(in%sx1+2,in%sx2),t2(in%sx1+2,in%sx2),t3(in%sx1+2,in%sx2), &
            STAT=iostatus)
  IF (iostatus>0) STOP "could not allocate memory"

  CALL curesetvectors () 

  iTensor=3
  CALL cutensormemset (%VAL(iTensor))

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! -     construct pre-stress structure
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  IF (ALLOCATED(in%stresslayer)) THEN
     CALL tensorstructure(in%stressstruc,in%stresslayer,in%dx3)
     DEALLOCATE(in%stresslayer)
     
     DO k=1,in%sx3/2
        tau(:,:,k)=(-1._4) .times. in%stressstruc(k)%t
     END DO
     DEALLOCATE(in%stressstruc)
  END IF

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! -     construct linear viscoelastic structure
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  IF (ALLOCATED(in%linearlayer)) THEN
     CALL viscoelasticstructure(in%linearstruc,in%linearlayer,in%dx3)
     DEALLOCATE(in%linearlayer)
  END IF

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! -   construct nonlinear viscoelastic structure
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  IF (ALLOCATED(in%nonlinearlayer)) THEN
     CALL viscoelasticstructure(in%nonlinearstruc,in%nonlinearlayer,in%dx3)
     DEALLOCATE(in%nonlinearlayer)
  END IF

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! -   construct nonlinear fault creep structure (rate-strenghtening)
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  IF (ALLOCATED(in%faultcreeplayer)) THEN
     CALL viscoelasticstructure(in%faultcreepstruc,in%faultcreeplayer,in%dx3)
     DEALLOCATE(in%faultcreeplayer)
  END IF

  ! first event
  e=1
  ! first output
  oi=1;
  ! initial condition
  t=0

  ! sources
  CALL dislocations(in%events(e),in%lambda,in%mu,in%beta,in%sx1,in%sx2,in%sx3, &
                    in%dx1,in%dx2,in%dx3,v1,v2,v3,t1,t2,t3,tau)

  iTensor = 1  
  CALL copytau (tau, %VAL(in%sx1), %VAL(in%sx2), %VAL(in%sx3/2), %VAL(iTensor))
  
  CALL cucopytraction (t3, %VAL(in%sx1), %VAL(in%sx2), iRet)
  
  CALL traction(in%mu,in%events(e),in%sx1,in%sx2,in%dx1,in%dx2,t,0.d0,t3)

 
  PRINT '("coseismic event ",I3.3)', e
  PRINT 0990

  ! test the presence of dislocations for coseismic calculation
  IF ((in%events(e)%nt .NE. 0) .OR. &
      (in%events(e)%ns .NE. 0) .OR. &
      (in%events(e)%nm .NE. 0) .OR. &
      (in%events(e)%nl .NE. 0)) THEN

     ! apply the 3d elastic transfer function
     CALL greenfunctioncowling(v1,v2,v3,t1,t2,t3, &
                               in%dx1,in%dx2,in%dx3,in%lambda,in%mu,in%gam)

     ! add displacement from analytic solutions for small patches (avoid
     ! aliasing)
  END IF

  ! transfer solution
  CALL cufieldrep (%VAL(in%sx1+2),%VAL(in%sx2),%VAL(in%sx3/2))


  ! evaluate stress
  iTensor=2
  cuc1=0._4
  cuc2=-1._4 
  CALL cutensorfieldadd (%VAL(iTensor), %VAL(in%sx1),%VAL(in%sx2),%VAL(in%sx3/2), & 
                         %VAL(cuc1), %VAL(cuc2), sig, tau)

  iTensor = 1
  CALL custressupdatewrapper (%VAL(iTensor), %VAL(in%lambda), %VAL(in%mu), &
                              %VAL(in%dx1), %VAL(in%dx2), %VAL(in%dx3),    & 
                              %VAL(in%sx1), %VAL(in%sx2), %VAL(in%sx3/2), u1, u2, u3, sig)

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! -   export observation points, map view of displacement,
  ! -   map view of stress components, Coulomb stress on observation
  ! -   patches, and full displacement and stress field.
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#ifdef GRD
  IF (in%isoutputgrd) THEN
     CALL exportgrd(u1,u2,u3,in%sx1,in%sx2,in%sx3/2,in%dx1,in%dx2,in%dx3,in%oz,in%x0,in%y0,in%wdir,0)
  END IF
#endif

  IF (ALLOCATED(in%ptsname)) THEN
     CALL exportpoints(u1,u2,u3,sig,in%sx1,in%sx2,in%sx3/2,in%dx1,in%dx2,in%dx3, &
          in%opts,in%ptsname,0._8,in%wdir,.true.,in%x0,in%y0,in%rot)
  END IF

  iTensor=2
  CALL cutensoramp (%VAL(iTensor),%VAL(in%sx1),%VAL(in%sx2),%VAL(in%sx3/2),dAmp)
  dAmp=dAmp*DBLE(in%dx1*in%dx2*in%dx3)
  PRINT 1101,0,0._8,0._8,0._8,0._8,0._8,in%interval,0._8,dAmp
 
  IF (in%interval .LE. 0) THEN
     GOTO 100 ! no time integration
  END IF

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! -   start the relaxation
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  iTensor=2
  CALL cutensormemset (%VAL(iTensor))

  DO i=1,ITERATION_MAX
     IF (t .GE. in%interval) GOTO 100 ! proper exit

#ifdef PAPI_PROF
   ! start timer
    CALL papistartprofiling (cTimerName)
#endif
     
     ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     ! predictor
     ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

     ! initialize large time step
     tm=STEP_MAX;
     maxwell(:)=STEP_MAX;
     
     ! active mechanism flag
     mech(:)=0

     ! initialize no forcing term in tensor space
     iTensor=2
     CALL cutensormemset (%VAL(iTensor))

     ! power density from three mechanisms (linear and power-law viscosity 
     ! and fault creep)
     ! 1- linear viscosity

#ifdef PAPI_PROF
    cTimerName = 'Eigenstress'
    CALL papistartprofiling (cTimerName)
#endif 

     IF (ALLOCATED(in%linearstruc)) THEN
        CALL viscouseigenstress(in%mu,in%linearstruc,in%linearweakzone,in%nlwz, &
             sig,in%sx1,in%sx2,in%sx3/2, &
             in%dx1,in%dx2,in%dx3,moment,0.01_8,MAXWELLTIME=maxwell(1))
        mech(1)=1
     END IF
     
     ! 2- powerlaw viscosity
     IF (ALLOCATED(in%nonlinearstruc)) THEN
        CALL viscouseigenstress(in%mu,in%nonlinearstruc,in%nonlinearweakzone,in%nnlwz, &
             sig,in%sx1,in%sx2,in%sx3/2, &
             in%dx1,in%dx2,in%dx3,moment,0.01_8,MAXWELLTIME=maxwell(2))
      mech(2)=1
     END IF

     ! 3- nonlinear fault creep with rate-strengthening friction
     IF (ALLOCATED(in%faultcreepstruc)) THEN
        DO k=1,in%np
           CALL frictioneigenstress(in%n(k)%x,in%n(k)%y,in%n(k)%z, &
                in%n(k)%width,in%n(k)%length, &
                in%n(k)%strike,in%n(k)%dip,in%n(k)%rake,in%beta, &
                sig,in%mu,in%faultcreepstruc, &
                in%sx1,in%sx2,in%sx3/2,in%dx1,in%dx2,in%dx3, &
                moment,maxwelltime=maxwell(3))
        END DO
        mech(3)=1
     END IF

#ifdef PAPI_PROF
    CALL papiendprofiling (cTimerName)
#endif

     ! identify the required time step
     tm=1._8/(REAL(mech(1))/maxwell(1)+ &
              REAL(mech(2))/maxwell(2)+ &
              REAL(mech(3))/maxwell(3))
     ! force finite time step
     tm=MIN(tm,STEP_MAX)

     ! modify
     IF ((in%inter%ns .GT. 0) .OR. (in%inter%nt .GT. 0)) THEN
        IF (tm .EQ. STEP_MAX) THEN
           ! no relaxation occurs, pick a small integration time
           tm=in%interval/20._8
        END IF
     END IF

     ! choose an integration time step
     CALL integrationstep(tm,Dt,t,oi,in%odt,in%skip,in%tscale,in%events,e,in%ne)

     iTensor=4
     cuc1=0._4
     cuc2=1._4
     CALL cutensorfieldadd (%VAL(iTensor), %VAL(in%sx1),%VAL(in%sx2),%VAL(in%sx3/2), %VAL(cuc1), %VAL(cuc2), sig, moment)

     CALL curesetvectors () 

#ifdef PAPI_PROF
     cTimerName = 'bodyforce'
     CALL papistartprofiling(cTimerName)
#endif

     iTensor=1
     CALL cubodyforceswrapper (%VAL(iTensor), %VAL(in%dx1), %VAL(in%dx2), %VAL(in%dx3), &
                               %VAL(in%sx1),%VAL(in%sx2),%VAL(in%sx3/2), v1,v2,v3,t1,t2,t3)

#ifdef PAPI_PROF
    CALL papiendprofiling (cTimerName)
#endif


     CALL cucopytraction (t3, %VAL(in%sx1), %VAL(in%sx2), iRet)
     
     ! add time-dependent surface loads
     CALL traction(in%mu,in%events(e),in%sx1,in%sx2,in%dx1,in%dx2,t,Dt/2.d8,t3,rate=.TRUE.)

#ifdef PAPI_PROF
     cTimerName = 'green'
     CALL papistartprofiling(cTimerName)
#endif

     CALL greenfunctioncowling(v1,v2,v3,t1,t2,t3,in%dx1,in%dx2,in%dx3,in%lambda,in%mu,in%gam)

#ifdef PAPI_PROF
    CALL papiendprofiling (cTimerName)
#endif

#ifdef PAPI_PROF
     cTimerName = 'fieldadd'
     CALL papistartprofiling(cTimerName)
#endif 

     ! v1,v2,v3 contain the predictor displacement
     iTensor = 2
     cuC1 = REAL(Dt/2)
     cuC2 = 1._4 
     CALL cufieldadd (%VAL(iTensor), v1, v2, v3, u1, u2, u3, %VAL(in%sx1+2),%VAL(in%sx2),%VAL(in%sx3/2), %VAL(cuC1), %VAL(cuC2))

#ifdef PAPI_PROF
     CALL papiendprofiling(cTimerName)
#endif 

     iTensor=2
     cuc1=-REAL(Dt/2)
     cuc2=-1._4
     CALL cutensorfieldadd (%VAL(iTensor), %VAL(in%sx1),%VAL(in%sx2),%VAL(in%sx3/2), %VAL(cuc1), %VAL(cuc2), sig, tau)

#ifdef PAPI_PROF
     cTimerName = 'stress'
     CALL papistartprofiling(cTimerName)
#endif

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! corrector
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  iTensor = 2
  CALL custressupdatewrapper (%VAL(iTensor), %VAL(in%lambda), %VAL(in%mu), &
                              %VAL(in%dx1), %VAL(in%dx2), %VAL(in%dx3),    & 
                              %VAL(in%sx1), %VAL(in%sx2), %VAL(in%sx3/2), v1, v2, v3, sig)

#ifdef PAPI_PROF
    CALL papiendprofiling (cTimerName)
#endif

     ! reinitialize moment density tensor
     iTensor=2
     CALL cutensormemset (%VAL(iTensor)) 
     
     IF (ALLOCATED(in%linearstruc)) THEN
        ! linear viscosity
!        v1=0
        CALL viscouseigenstress(in%mu,in%linearstruc,in%linearweakzone,in%nlwz,sig, &
             in%sx1,in%sx2,in%sx3/2, &
             in%dx1,in%dx2,in%dx3,moment,0.01_8,GAMMA=v1)
        
     END IF
    
#ifdef PAPI_PROF
     cTimerName = 'Eigenstress'
     CALL papistartprofiling(cTimerName)
#endif
 
     IF (ALLOCATED(in%nonlinearstruc)) THEN
        ! powerlaw viscosity
!        v1=0
        CALL viscouseigenstress(in%mu,in%nonlinearstruc,in%nonlinearweakzone,in%nnlwz,sig, &
             in%sx1,in%sx2,in%sx3/2, &
             in%dx1,in%dx2,in%dx3,moment,0.01_8,GAMMA=v1)
 
     END IF

#ifdef PAPI_PROF
    CALL papiendprofiling (cTimerName)
#endif
        
 
     ! nonlinear fault creep with rate-strengthening friction
     IF (ALLOCATED(in%faultcreepstruc)) THEN

        ! use v1 as placeholders for the afterslip planes
        DO k=1,in%np
           ! one may use optional arguments ...,VEL=v1) to convert
           ! fault slip to eigenstrain (scalar)
           CALL frictioneigenstress(in%n(k)%x,in%n(k)%y,in%n(k)%z, &
                in%n(k)%width,in%n(k)%length, &
                in%n(k)%strike,in%n(k)%dip,in%n(k)%rake,in%beta, &
                sig,in%mu,in%faultcreepstruc,in%sx1,in%sx2,in%sx3/2, &
                in%dx1,in%dx2,in%dx3,moment)

           ! keep track if afterslip instantaneous velocity
           CALL monitorfriction(in%n(k)%x,in%n(k)%y,in%n(k)%z, &
                in%n(k)%width,in%n(k)%length, &
                in%n(k)%strike,in%n(k)%dip,in%n(k)%rake,in%beta, &
                in%sx1,in%sx2,in%sx3,in%dx1,in%dx2,in%dx3, &
                sig,in%faultcreepstruc,in%n(k)%patch)
        END DO

     END IF

     ! interseismic loading
     IF ((in%inter%ns .GT. 0) .OR. (in%inter%nt .GT. 0)) THEN
        ! vectors v1,v2,v3 are not affected.
        CALL dislocations(in%inter,in%lambda,in%mu,in%beta,in%sx1,in%sx2,in%sx3, &
             in%dx1,in%dx2,in%dx3,v1,v2,v3,t1,t2,t3,tau,eigenstress=moment)
     END IF
     

     CALL curesetvectors () 

#ifdef PAPI_PROF
     cTimerName = 'bodyforce'
     CALL papistartprofiling(cTimerName)
#endif

     iTensor=2
     CALL cubodyforceswrapper (%VAL(iTensor), %VAL(in%dx1), %VAL(in%dx2), %VAL(in%dx3), &
                               %VAL(in%sx1),%VAL(in%sx2),%VAL(in%sx3/2), v1,v2,v3,t1,t2,t3)

#ifdef PAPI_PROF
    CALL papiendprofiling (cTimerName)
#endif


     CALL cucopytraction (t3, %VAL(in%sx1), %VAL(in%sx2), iRet)

#ifdef PAPI_PROF
     cTimerName = 'green'
     CALL papistartprofiling(cTimerName)
#endif

     CALL greenfunctioncowling(v1,v2,v3,t1,t2,t3,in%dx1,in%dx2,in%dx3,in%lambda,in%mu,in%gam)

#ifdef PAPI_PROF
    CALL papiendprofiling (cTimerName)
#endif

#ifdef PAPI_PROF
     cTimerName = 'fieldadd'
     CALL papistartprofiling(cTimerName)
#endif 

     iTensor = 1
     cuC1 = 1._4
     cuC2 = REAL(Dt) 
     CALL cufieldadd (%VAL(iTensor), u1, u2, u3, v1, v2, v3, %VAL(in%sx1+2),%VAL(in%sx2),%VAL(in%sx3/2), %VAL(cuC1), %VAL(cuC2))

#ifdef PAPI_PROF
     CALL papiendprofiling(cTimerName)
#endif 


     iTensor=5
     cuc1=1._4
     cuc2=REAL(Dt)
     CALL cutensorfieldadd (%VAL(iTensor), %VAL(in%sx1),%VAL(in%sx2),%VAL(in%sx3/2), %VAL(cuc1), %VAL(cuc2), tau, moment)
     
     CALL frictionadd(in%np,in%n,Dt)
     
     ! time increment
     t=t+Dt
    
     ! next event
     IF (e .LT. in%ne) THEN
        IF (abs(t-in%events(e+1)%time) .LT. 1e-6) THEN
           e=e+1
           in%events(e)%i=i

           PRINT '("coseismic event ",I3.3)', e
           PRINT 0990
          
           CALL curesetvectors ()
           iTensor = 0
           CALL copytau (tau, %VAL(in%sx1), %VAL(in%sx2), %VAL(in%sx3/2), %VAL(iTensor)) 
           
           CALL dislocations(in%events(e),in%lambda,in%mu, &
                in%beta,in%sx1,in%sx2,in%sx3, &
                in%dx1,in%dx2,in%dx3,v1,v2,v3,t1,t2,t3,tau)

           iTensor = 1
           CALL copytau (tau, %VAL(in%sx1), %VAL(in%sx2), %VAL(in%sx3/2), %VAL(iTensor))
           CALL cucopytraction (t3, %VAL(in%sx1), %VAL(in%sx2), iRet)
           
           CALL traction(in%mu,in%events(e),in%sx1,in%sx2,in%dx1,in%dx2,t,0.d0,t3)

           ! apply the 3d elastic transfert function
           CALL greenfunctioncowling(v1,v2,v3,t1,t2,t3, &
                in%dx1,in%dx2,in%dx3,in%lambda,in%mu,in%gam)

     iTensor = 1
     cuC1 = 1._4
     cuC2 = 1._4 
     CALL cufieldadd (%VAL(iTensor), u1, u2, u3, v1, v2, v3, %VAL(in%sx1+2),%VAL(in%sx2),%VAL(in%sx3/2), %VAL(cuC1), %VAL(cuC2))
        END IF
     END IF

     iTensor=2
     cuc1=0._4
     cuc2=-1._4
     CALL cutensorfieldadd (%VAL(iTensor), %VAL(in%sx1),%VAL(in%sx2),%VAL(in%sx3/2), %VAL(cuc1), %VAL(cuc2), sig, tau)


#ifdef PAPI_PROF
     cTimerName = 'stress'
     CALL papistartprofiling(cTimerName)
#endif

  iTensor = 1
  CALL custressupdatewrapper (%VAL(iTensor), %VAL(in%lambda), %VAL(in%mu), &
                              %VAL(in%dx1), %VAL(in%dx2), %VAL(in%dx3),    & 
                              %VAL(in%sx1), %VAL(in%sx2), %VAL(in%sx3/2), u1, u2, u3, sig)

#ifdef PAPI_PROF
     CALL papiendprofiling(cTimerName) 
#endif
     ! points are exported at all time steps
     IF (ALLOCATED(in%ptsname)) THEN
        CALL exportpoints(u1,u2,u3,sig,in%sx1,in%sx2,in%sx3/2,in%dx1,in%dx2,in%dx3, &
             in%opts,in%ptsname,t,in%wdir,.FALSE.,in%x0,in%y0,in%rot)
     END IF

     iTensor=1
     CALL cutensoramp (%VAL(iTensor),%VAL(in%sx1),%VAL(in%sx2),%VAL(in%sx3/2),dAmp, moment)
     iTensor=2
     CALL cutensoramp (%VAL(iTensor),%VAL(in%sx1),%VAL(in%sx2),%VAL(in%sx3/2),dAmp1, tau)
     dAmp=dAmp*DBLE(in%dx1*in%dx2*in%dx3)
     dAmp1=dAmp1*DBLE(in%dx1*in%dx2*in%dx3)
     PRINT 1101,i,Dt,maxwell,t,in%interval, dAmp, dAmp1
        ! update output counter
     oi=oi+1

#ifdef PAPI_PROF
    cTimerName = 'relax'
    CALL papiendprofiling (cTimerName)
#endif

  END DO

100 CONTINUE

  DO i=1,in%ne
     IF (ALLOCATED(in%events(i)%s))  DEALLOCATE(in%events(i)%s,in%events(i)%sc)
     IF (ALLOCATED(in%events(i)%ts)) DEALLOCATE(in%events(i)%ts,in%events(i)%tsc)
  END DO
  IF (ALLOCATED(in%events)) DEALLOCATE(in%events)

  CALL cudeinit()

  ! free memory
  IF (ALLOCATED(gamma)) DEALLOCATE(gamma)
  IF (ALLOCATED(in%opts)) DEALLOCATE(in%opts)
  IF (ALLOCATED(in%ptsname)) DEALLOCATE(in%ptsname)
  IF (ALLOCATED(in%op)) DEALLOCATE(in%op)
  IF (ALLOCATED(in%sop)) DEALLOCATE(in%sop)
  IF (ALLOCATED(in%n)) DEALLOCATE(in%n)
  IF (ALLOCATED(in%stressstruc)) DEALLOCATE(in%stressstruc)
  IF (ALLOCATED(in%stresslayer)) DEALLOCATE(in%stresslayer)
  IF (ALLOCATED(in%linearstruc)) DEALLOCATE(in%linearstruc)
  IF (ALLOCATED(in%linearlayer)) DEALLOCATE(in%linearlayer)
  IF (ALLOCATED(in%linearweakzone)) DEALLOCATE(in%linearweakzone)
  IF (ALLOCATED(in%nonlinearstruc)) DEALLOCATE(in%nonlinearstruc)
  IF (ALLOCATED(in%nonlinearlayer)) DEALLOCATE(in%nonlinearlayer)
  IF (ALLOCATED(in%nonlinearweakzone)) DEALLOCATE(in%nonlinearweakzone)
  IF (ALLOCATED(in%faultcreepstruc)) DEALLOCATE(in%faultcreepstruc)
  IF (ALLOCATED(in%faultcreeplayer)) DEALLOCATE(in%faultcreeplayer)
  IF (ALLOCATED(sig)) DEALLOCATE(sig)
  IF (ALLOCATED(tau)) DEALLOCATE(tau)
  IF (ALLOCATED(moment)) DEALLOCATE(moment)
  IF (ALLOCATED(in%stresslayer)) DEALLOCATE(in%stresslayer)
  IF (ALLOCATED(in%linearlayer)) DEALLOCATE(in%linearlayer)
  IF (ALLOCATED(in%nonlinearlayer)) DEALLOCATE(in%nonlinearlayer)
  IF (ALLOCATED(in%faultcreeplayer)) DEALLOCATE(in%faultcreeplayer)
  IF (ALLOCATED(v1)) DEALLOCATE(v1,v2,v3,t1,t2,t3)
  IF (ALLOCATED(u1)) DEALLOCATE(u1,u2,u3)
  IF (ALLOCATED(inter1)) DEALLOCATE(inter1,inter2,inter3)


0990 FORMAT (" I  |   Dt   | tm(ve) | tm(pl) | tm(as) |     t/tmax     | power  |  C:E^i | ")
1100 FORMAT (I3.3," ",ES9.2E2,3ES9.2E2,ES9.2E2,"/",ES7.2E1,2ES9.2E2)
1101 FORMAT (I3.3,"*",ES9.2E2,3ES9.2E2,ES9.2E2,"/",ES7.2E1,2ES9.2E2)

CONTAINS

  !--------------------------------------------------------------------
  !> subroutine dislocations
  !! assigns equivalent body forces or moment density to simulate
  !! shear dislocations and fault opening. add the corresponding moment
  !! density in the cumulative relaxed moment so that fault slip does
  !! not reverse in the postseismic time.
  !--------------------------------------------------------------------
  SUBROUTINE dislocations(event,lambda,mu,beta,sx1,sx2,sx3,dx1,dx2,dx3, &
                          v1,v2,v3,t1,t2,t3,tau,factor,eigenstress)
    TYPE(EVENT_STRUC), INTENT(IN) :: event
    INTEGER, INTENT(IN) :: sx1,sx2,sx3
    REAL*8, INTENT(IN) :: lambda,mu,beta,dx1,dx2,dx3
    REAL*4, DIMENSION(:,:,:), INTENT(INOUT) :: v1,v2,v3
    REAL*4, DIMENSION(:,:), INTENT(INOUT) :: t1,t2,t3
    TYPE(TENSOR), DIMENSION(:,:,:), INTENT(INOUT) :: tau
    REAL*8, INTENT(IN), OPTIONAL :: factor
    TYPE(TENSOR), DIMENSION(:,:,:), INTENT(INOUT), OPTIONAL :: eigenstress
    
    INTEGER :: i
    REAL*8 :: slip_factor
    
    IF (PRESENT(factor)) THEN
       slip_factor=factor
    ELSE
       slip_factor=1._8
    END IF
    
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! -             load shear dislocations
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    IF (.NOT. (PRESENT(eigenstress))) THEN
       ! forcing term in equivalent body force
       DO i=1,event%ns
          IF (in%nyquist*MIN(in%dx1,in%dx2,in%dx3).LT.event%s(i)%length .OR. &
              in%nyquist*MIN(in%dx1,in%dx2,in%dx3).LT.event%s(i)%width) THEN
             ! adding sources in the space domain
             CALL source(mu,slip_factor*event%s(i)%slip, &
                  event%s(i)%x,event%s(i)%y,event%s(i)%z, &
                  event%s(i)%width,event%s(i)%length, &
                  event%s(i)%strike,event%s(i)%dip,event%s(i)%rake, &
                  event%s(i)%beta,sx1,sx2,sx3,dx1,dx2,dx3,v1,v2,v3,t1,t2,t3)
          END IF
       END DO
    ELSE
       ! forcing term in moment density
       DO i=1,event%ns
          CALL momentdensityshear(mu,slip_factor*event%s(i)%slip, &
               event%s(i)%x,event%s(i)%y,event%s(i)%z, &
               event%s(i)%width,event%s(i)%length, &
               event%s(i)%strike,event%s(i)%dip,event%s(i)%rake, &
               event%s(i)%beta,sx1,sx2,sx3/2,dx1,dx2,dx3,eigenstress)
       END DO
    END IF

    DO i=1,event%ns
       ! remove corresponding eigenmoment
       CALL momentdensityshear(mu,slip_factor*event%s(i)%slip, &
            event%s(i)%x,event%s(i)%y,event%s(i)%z, &
            event%s(i)%width,event%s(i)%length, &
            event%s(i)%strike,event%s(i)%dip,event%s(i)%rake, &
            event%s(i)%beta,sx1,sx2,sx3/2,dx1,dx2,dx3,tau)
    END DO
    
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! -             load tensile cracks
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    IF (.NOT. (PRESENT(eigenstress))) THEN
       ! forcing term in equivalent body force
       DO i=1,event%nt
          IF (in%nyquist*MIN(in%dx1,in%dx2,in%dx3).LT.event%tsc(i)%length .OR. &
              in%nyquist*MIN(in%dx1,in%dx2,in%dx3).LT.event%tsc(i)%width) THEN
             ! adding sources in the space domain
             CALL tensilesource(lambda,mu,slip_factor*event%ts(i)%opening, &
                  event%ts(i)%x,event%ts(i)%y,event%ts(i)%z, &
                  event%ts(i)%width,event%ts(i)%length, &
                  event%ts(i)%strike,event%ts(i)%dip, &
                  beta,sx1,sx2,sx3,dx1,dx2,dx3,v1,v2,v3)
          END IF
       END DO
    ELSE
       ! forcing term in moment density
       DO i=1,event%nt
          CALL momentdensitytensile(lambda,mu,slip_factor*event%ts(i)%opening, &
               event%ts(i)%x,event%ts(i)%y,event%ts(i)%z,&
               event%ts(i)%width,event%ts(i)%length, &
               event%ts(i)%strike,event%ts(i)%dip,event%ts(i)%rake, &
               beta,sx1,sx2,sx3/2,dx1,dx2,dx3,eigenstress)
       END DO
    END IF

    DO i=1,event%nt
       ! removing corresponding eigenmoment
       CALL momentdensitytensile(lambda,mu,slip_factor*event%ts(i)%opening, &
            event%ts(i)%x,event%ts(i)%y,event%ts(i)%z,&
            event%ts(i)%width,event%ts(i)%length, &
            event%ts(i)%strike,event%ts(i)%dip,event%ts(i)%rake, &
            beta,sx1,sx2,sx3/2,dx1,dx2,dx3,tau)
    END DO

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! -             load point dilatation sources
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    IF (.NOT. (PRESENT(eigenstress))) THEN
       ! forcing term in equivalent body force
       DO i=1,event%nm
          ! adding sources in the space domain
          CALL mogisource(lambda,mu,slip_factor*event%m(i)%slip, &
               event%m(i)%x,event%m(i)%y,event%m(i)%z, &
               sx1,sx2,sx3,dx1,dx2,dx3,v1,v2,v3)
       END DO
    ELSE
       ! forcing term in moment density
       DO i=1,event%nm
          CALL momentdensitymogi(lambda,mu,slip_factor*event%m(i)%slip, &
               event%m(i)%x,event%m(i)%y,event%m(i)%z, &
               sx1,sx2,sx3/2,dx1,dx2,dx3,eigenstress)
       END DO
    END IF

    DO i=1,event%nm
       ! remove corresponding eigenmoment
       CALL momentdensitymogi(lambda,mu,slip_factor*event%m(i)%slip, &
            event%m(i)%x,event%m(i)%y,event%m(i)%z, &
            sx1,sx2,sx3/2,dx1,dx2,dx3,tau)
    END DO
    
  END SUBROUTINE dislocations

  !--------------------------------------------------------------------
  !> subroutine dislocations_disp
  !! evaluate the displacement due to the motion of dislocation (shear
  !! and opening)
  !!
  !! \author sylvain barbot (01/01/08) - original form 
  !--------------------------------------------------------------------
  SUBROUTINE dislocations_disp(event,lambda,mu,sx1,sx2,sx3,dx1,dx2,dx3, &
                          v1,v2,v3,factor)
    TYPE(EVENT_STRUC), INTENT(IN) :: event
    INTEGER, INTENT(IN) :: sx1,sx2,sx3
    REAL*8, INTENT(IN) :: lambda,mu,dx1,dx2,dx3
    REAL*4, DIMENSION(:,:,:), INTENT(INOUT) :: v1,v2,v3
    REAL*8, INTENT(IN), OPTIONAL :: factor
    
    INTEGER :: i
    REAL*8 :: slip_factor
    
    IF (PRESENT(factor)) THEN
       slip_factor=factor
    ELSE
       slip_factor=1._8
    END IF
    
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! -             shear dislocations
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    ! forcing term in equivalent body force
    DO i=1,event%ns
       IF (in%nyquist*MIN(in%dx1,in%dx2,in%dx3).GE.event%s(i)%length .OR. &
           in%nyquist*MIN(in%dx1,in%dx2,in%dx3).GE.event%s(i)%width) THEN
          ! adding sources in the space domain
          CALL grnfct_okada(lambda,mu,slip_factor*event%sc(i)%slip,0._8, &
               event%sc(i)%x,event%sc(i)%y,event%sc(i)%z, &
               event%sc(i)%width,event%sc(i)%length, &
               event%sc(i)%strike,event%sc(i)%dip,event%sc(i)%rake, &
               sx1,sx2,sx3/2,dx1,dx2,dx3,v1,v2,v3)
       END IF
    END DO

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! -             load tensile cracks
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    ! forcing term in equivalent body force
    DO i=1,event%nt
       IF (in%nyquist*MIN(in%dx1,in%dx2,in%dx3).GE.event%tsc(i)%length .OR. &
           in%nyquist*MIN(in%dx1,in%dx2,in%dx3).GE.event%tsc(i)%width) THEN
          ! adding sources in the space domain
          CALL grnfct_okada(lambda,mu,0._8,slip_factor*event%tsc(i)%opening, &
               event%tsc(i)%x,event%tsc(i)%y,event%tsc(i)%z, &
               event%tsc(i)%width,event%tsc(i)%length, &
               event%tsc(i)%strike,event%tsc(i)%dip,0._8, &
               sx1,sx2,sx3/2,dx1,dx2,dx3,v1,v2,v3)
       END IF
    END DO

  END SUBROUTINE dislocations_disp

  !--------------------------------------------------------------------
  !> function IsOutput
  !! checks if output should be written based on user choices: if output
  !! time interval (odt) is positive, output is written only if time
  !! is an integer of odt. If odt is negative output is written at times
  !! corresponding to internally chosen time steps.
  !!
  !! @return IsOutput is true only at discrete intervals (skip=0,odt>0),
  !! or at every "skip" computational steps (skip>0,odt<0),
  !! or anytime a coseismic event occurs
  !
  ! Sylvain Barbot (07/06/09) - original form
  !--------------------------------------------------------------------
  LOGICAL FUNCTION isoutput(skip,t,i,odt,oi,etime)
    INTEGER, INTENT(IN) :: skip,i,oi
    REAL*8, INTENT(IN) :: t,odt,etime

    IF (((0 .EQ. skip) .AND. (abs(t-oi*odt) .LT. 1e-6*odt)) .OR. &
        ((0 .LT. skip) .AND. (MOD(i-1,skip) .EQ. 0)) .OR. &
         (abs(t-etime) .LT. 1e-6)) THEN
       isoutput=.TRUE.
    ELSE
       isoutput=.FALSE.
    END IF

  END FUNCTION isoutput

  !--------------------------------------------------------------------
  !> subroutine IntegrationStep
  !! find the time-integration forward step for the predictor-corrector
  !! scheme.
  !!
  !! input file line
  !!
  !!    time interval, (positive dt step) or (negative skip and scaling)
  !!
  !! can be filled by either 1)
  !!
  !!   T, dt
  !!
  !! where T is the time interval of the simulation and dt is the
  !! output time step, or 2)
  !!
  !!   T, -n, t_s
  !!
  !! where n indicates the number of computational steps before 
  !! outputing results, t_s is a scaling applied to internally
  !! computed time step.
  !!
  !! for case 1), an optimal time step is evaluated internally to
  !! ensure stability (t_m/10) of time integration. The actual
  !! time step Dt is chosen as
  !!
  !!    Dt = min( t_m/10, ((t%odt)+1)*odt-t )
  !!
  !! where t is the current time in the simulation. regardless of 
  !! time step Dt, results are output if t is a multiple of dt.
  !!
  !! for case 2), the time step is chosen internally based on an 
  !! estimate of the relaxation time (t_m/10). Results are output
  !! every n steps. The actual time step is chosen as
  !!
  !!    Dt = min( t_m/10*t_s, t(next event)-t )
  !!
  !! where index is the number of computational steps after a coseismic
  !! event and t(next event) is the time of the next coseismic event.
  !!
  !! \author sylvain barbot (01/01/08) - original form 
  !--------------------------------------------------------------------
  SUBROUTINE integrationstep(tm,Dt,t,oi,odt,skip,tscale,events,e,ne)
    REAL*8, INTENT(INOUT) :: tm,Dt,odt
    REAL*8, INTENT(IN) :: t,tscale
    INTEGER, INTENT(IN) :: oi,e,ne,skip
    TYPE(EVENT_STRUC), INTENT(IN), DIMENSION(:) :: events

    ! output at optimal computational intervals
    Dt=tm/10._8

    ! reduce time in case something happens in [ t, t+Dt ]
    IF (0 .EQ. skip) THEN
       ! reduce time step so that t+Dt is time at next 
       ! user-required output time
       IF ((t+Dt) .GE. (dble(oi)*odt)-Dt*0.04d0) THEN
          ! pick a smaller time step to reach :
          ! integers of odt
          Dt=dble(oi)*odt-t
       END IF
    ELSE
       ! scale the estimate of optimal time step
       Dt=Dt*tscale

       ! reduce time step so that t+Dt is time to next event
       IF (e .LT. ne) THEN
          IF ((t+Dt-events(e+1)%time) .GE. 0._8) THEN
             ! pick a smaller time step to reach 
             ! next event time
             Dt=events(e+1)%time-t
          END IF
       END IF
    END IF

  END SUBROUTINE integrationstep


END PROGRAM relax
