!-----------------------------------------------------------------------
! Copyright 2007-2013, Sylvain Barbot
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
  !! (07-25-13) - include cylindrical and spherical ductile zones  <br>
  !! (01-03-14) - accelerate processing of ductile zones           <br>
  !!
  !! \todo 
  !!   - homogenize VTK output so that geometry of events match event index
  !!   - write the code for MPI multi-thread
  !!   - fix the vtk export to grid for anisotropic sampling
  !!   - export position of observation points to long/lat in opts-geo.dat
  !!   - check the projected output on the south hemisphere
  !!   - check the fully-relaxed afterslip for uniform stress change
  !!   - include topography of parameter interface
  !!   - export afterslip output in VTK legacy format (binary)
  !!   - export ductile zones for cylindrical and spherical geometries
  !------------------------------------------------------------------------

#include "include.f90"
  
SUBROUTINE relaxlite(in,gps,isverbose)

  USE types
  USE green
  USE green_space
  USE elastic3d
  USE viscoelastic3d
  USE friction3d
  USE export

  IMPLICIT NONE

  TYPE(SIMULATION_STRUCT), INTENT(INOUT) :: in
  TYPE(MANIFOLD_STRUCT), INTENT(OUT) :: gps(in%npts)
  LOGICAL, INTENT(IN) :: isverbose
  
  INTEGER, PARAMETER :: ITERATION_MAX = 9999
  REAL*8, PARAMETER :: STEP_MAX = 1e7

  INTEGER :: i,k,e,oi,iostatus,mech(3)
#ifdef FFTW3_THREADS
  INTEGER :: iret
!$  INTEGER :: omp_get_max_threads
#endif
  REAL*8 :: maxwell(3)
  CHARACTER(256) :: filename,title,name

  CHARACTER(4) :: digit4
  REAL*8 :: t,Dt,tm
  
  ! arrays
  REAL*4, DIMENSION(:,:,:), ALLOCATABLE :: v1,v2,v3,u1,u2,u3,gamma
  REAL*4, DIMENSION(:,:,:), ALLOCATABLE :: lineardgammadot0,nonlineardgammadot0
  REAL*4, DIMENSION(:,:), ALLOCATABLE :: t1,t2,t3
  TYPE(TENSOR), DIMENSION(:,:,:), ALLOCATABLE :: tau,sig,moment

#ifdef FFTW3_THREADS
  CALL sfftw_init_threads(iret)
#ifdef _OPENMP
  CALL sfftw_plan_with_nthreads(omp_get_max_threads())
#else
  CALL sfftw_plan_with_nthreads(4)
#endif
#endif

  ! allocate space for time series simulations
  DO k=1,in%npts
     gps(k)%nepochs=0
     ALLOCATE(gps(k)%t(1+ITERATION_MAX), &
              gps(k)%u1(1+ITERATION_MAX), &
              gps(k)%u2(1+ITERATION_MAX), &
              gps(k)%u3(1+ITERATION_MAX), STAT=iostatus)
     IF (iostatus>0) STOP "could not allocate prediction array"
  END DO

  ! allocate memory
  ALLOCATE (v1(in%sx1+2,in%sx2,in%sx3),v2(in%sx1+2,in%sx2,in%sx3),v3(in%sx1+2,in%sx2,in%sx3), &
            u1(in%sx1+2,in%sx2,in%sx3/2),u2(in%sx1+2,in%sx2,in%sx3/2),u3(in%sx1+2,in%sx2,in%sx3/2), &
            tau(in%sx1,in%sx2,in%sx3/2),sig(in%sx1,in%sx2,in%sx3/2),gamma(in%sx1+2,in%sx2,in%sx3/2), &
            t1(in%sx1+2,in%sx2),t2(in%sx1+2,in%sx2),t3(in%sx1+2,in%sx2),STAT=iostatus)
  IF (iostatus>0) STOP "could not allocate memory"

  v1=0;v2=0;v3=0;u1=0;u2=0;u3=0;gamma=0;t1=0;t2=0;t3=0
  i=0
  CALL tensorfieldadd(tau,tau,in%sx1,in%sx2,in%sx3/2,c1=0._4,c2=0._4)

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! -
  ! -     construct pre-stress structure
  ! -
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  IF (ALLOCATED(in%stresslayer)) THEN
     ! depth-dependent background stress
     CALL tensorstructure(in%stressstruc,in%stresslayer,in%dx3)
  ELSE
     ! background stress is zero
     in%stressstruc(:)%t=tensor(0._4,0._4,0._4,0._4,0._4,0._4)
  END IF
  DO k=1,in%sx3/2
     tau(:,:,k)=(-1._4) .times. in%stressstruc(k)%t
  END DO

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! -
  ! -     first event
  ! -
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  e=1
  ! first output
  oi=1;
  ! initial condition
  t=0

  ! sources
  CALL dislocations(in%events(e),in%lambda,in%mu,in%beta,in%sx1,in%sx2,in%sx3, &
                    in%dx1,in%dx2,in%dx3,v1,v2,v3,t1,t2,t3,tau)
  CALL traction(in%mu,in%events(e),in%sx1,in%sx2,in%dx1,in%dx2,t,0.d0,t3)
  
  IF (isverbose) THEN
     PRINT '("# event ",I3.3)', e
     PRINT 0990
  END IF

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
     CALL dislocations_disp(in%events(e),in%lambda,in%mu, &
                            in%sx1,in%sx2,in%sx3, &
                            in%dx1,in%dx2,in%dx3,v1,v2,v3)
  END IF
  
  ! transfer solution
  CALL fieldrep(u1,v1,in%sx1+2,in%sx2,in%sx3/2)
  CALL fieldrep(u2,v2,in%sx1+2,in%sx2,in%sx3/2)
  CALL fieldrep(u3,v3,in%sx1+2,in%sx2,in%sx3/2)

  ! evaluate stress
  CALL tensorfieldadd(sig,tau,in%sx1,in%sx2,in%sx3/2,c1=0._4,c2=-1._4)
  CALL stressupdate(u1,u2,u3,in%lambda,in%mu, &
                    in%dx1,in%dx2,in%dx3,in%sx1,in%sx2,in%sx3/2,sig)

  ! export the initial displacements
  CALL pts2series(u1,u2,u3,in%sx1,in%sx2,in%sx3,in%dx1,in%dx2,in%dx3,in%opts,0._8,1,gps)

  IF (ALLOCATED(in%ptsname)) THEN
     CALL exportpoints(u1,u2,u3,sig,in%sx1,in%sx2,in%sx3/2,in%dx1,in%dx2,in%dx3, &
                      in%opts,in%ptsname,0._8,in%wdir,.true.,in%x0,in%y0,in%rot)
  END IF

  WRITE (digit4,'(I4.4)') 0

  IF (isverbose) THEN
     PRINT 1101,0,0._8,0._8,0._8,0._8,0._8,in%interval,0._8,tensoramplitude(tau,in%dx1,in%dx2,in%dx3)
  END IF
  IF (in%interval .LE. 0) THEN
     GOTO 100 ! no time integration
  END IF

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! -     construct linear viscoelastic structure
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  IF (ALLOCATED(in%linearlayer)) THEN
     CALL viscoelasticstructure(in%linearstruc,in%linearlayer,in%dx3)
     !The caller has to take care of this memory release
     !DEALLOCATE(in%linearlayer)

     IF (0 .LT. in%nlwz) THEN
        ALLOCATE(lineardgammadot0(in%sx1,in%sx2,in%sx3/2),STAT=iostatus)
        IF (iostatus.GT.0) STOP "could not allocate lineardgammadot0"
        CALL builddgammadot0(in%sx1,in%sx2,in%sx3/2,in%dx1,in%dx2,in%dx3,0._8, &
                             in%nlwz,in%linearweakzone,lineardgammadot0)
     END IF
  END IF

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! -   construct nonlinear viscoelastic structure
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  IF (ALLOCATED(in%nonlinearlayer)) THEN
     CALL viscoelasticstructure(in%nonlinearstruc,in%nonlinearlayer,in%dx3)
     !The caller has to take care of this memory release
     !DEALLOCATE(in%nonlinearlayer)

     IF (0 .LT. in%nnlwz) THEN
        ALLOCATE(nonlineardgammadot0(in%sx1,in%sx2,in%sx3/2),STAT=iostatus)
        IF (iostatus.GT.0) STOP "could not allocate nonlineardgammadot0"
        CALL builddgammadot0(in%sx1,in%sx2,in%sx3/2,in%dx1,in%dx2,in%dx3,0._8, &
                             in%nnlwz,in%nonlinearweakzone,nonlineardgammadot0)
     END IF
  END IF

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! -   construct nonlinear fault creep structure (rate-strenghtening)
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  IF (ALLOCATED(in%faultcreeplayer)) THEN
     CALL viscoelasticstructure(in%faultcreepstruc,in%faultcreeplayer,in%dx3)
     !The caller has to take care of this memory release
     !DEALLOCATE(in%faultcreeplayer)
  END IF

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! -   start the relaxation
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ALLOCATE(moment(in%sx1,in%sx2,in%sx3/2),STAT=iostatus)
  IF (iostatus>0) STOP "could not allocate the mechanical structure"

  CALL tensorfieldadd(moment,moment,in%sx1,in%sx2,in%sx3/2,c1=0._4,c2=0._4)  

  DO i=1,ITERATION_MAX
     IF (t .GE. in%interval) GOTO 100 ! proper exit
     
     ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     ! predictor
     ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

     ! initialize large time step
     tm=STEP_MAX;
     maxwell(:)=STEP_MAX;
     
     ! active mechanism flag
     mech(:)=0

     ! initialize no forcing term in tensor space
     CALL tensorfieldadd(moment,moment,in%sx1,in%sx2,in%sx3/2,0._4,0._4)

     ! power density from three mechanisms (linear and power-law viscosity 
     ! and fault creep)
     ! 1- linear viscosity
     IF (ALLOCATED(in%linearstruc)) THEN
        IF (0 .LT. in%nlwz) THEN
           CALL viscouseigenstress(in%mu,in%linearstruc, &
                sig,in%stressstruc,in%sx1,in%sx2,in%sx3/2, &
                in%dx1,in%dx2,in%dx3,moment,DGAMMADOT0=lineardgammadot0,MAXWELLTIME=maxwell(1))
        ELSE
           CALL viscouseigenstress(in%mu,in%linearstruc, &
                sig,in%stressstruc,in%sx1,in%sx2,in%sx3/2, &
                in%dx1,in%dx2,in%dx3,moment,MAXWELLTIME=maxwell(1))
        END IF
        mech(1)=1
     END IF
     
     ! 2- powerlaw viscosity
     IF (ALLOCATED(in%nonlinearstruc)) THEN
        IF (0 .LT. in%nnlwz) THEN
           CALL viscouseigenstress(in%mu,in%nonlinearstruc, &
                sig,in%stressstruc,in%sx1,in%sx2,in%sx3/2, &
                in%dx1,in%dx2,in%dx3,moment,DGAMMADOT0=nonlineardgammadot0,MAXWELLTIME=maxwell(2))
        ELSE
           CALL viscouseigenstress(in%mu,in%nonlinearstruc, &
                sig,in%stressstruc,in%sx1,in%sx2,in%sx3/2, &
                in%dx1,in%dx2,in%dx3,moment,MAXWELLTIME=maxwell(2))
        END IF
        mech(2)=1
     END IF
     
     ! 3- nonlinear fault creep with rate-strengthening friction
     IF (ALLOCATED(in%faultcreepstruc)) THEN
        DO k=1,in%np
           CALL frictioneigenstress(in%n(k)%x,in%n(k)%y,in%n(k)%z, &
                in%n(k)%width,in%n(k)%length, &
                in%n(k)%strike,in%n(k)%dip,in%n(k)%rake,in%beta, &
                sig,in%stressstruc,in%mu,in%faultcreepstruc, &
                in%sx1,in%sx2,in%sx3/2,in%dx1,in%dx2,in%dx3, &
                moment,maxwelltime=maxwell(3))
        END DO
        mech(3)=1
     END IF

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

     CALL tensorfieldadd(sig,moment,in%sx1,in%sx2,in%sx3/2,c1=0.0_4,c2=1._4)
     
     v1=0;v2=0;v3=0;t1=0;t2=0;t3=0;
     CALL equivalentbodyforce(sig,in%dx1,in%dx2,in%dx3,in%sx1,in%sx2,in%sx3/2,v1,v2,v3,t1,t2,t3)

     ! add time-dependent surface loads
     CALL traction(in%mu,in%events(e),in%sx1,in%sx2,in%dx1,in%dx2,t,Dt/2.d8,t3,rate=.TRUE.)

     CALL greenfunctioncowling(v1,v2,v3,t1,t2,t3,in%dx1,in%dx2,in%dx3,in%lambda,in%mu,in%gam)
     
     ! v1,v2,v3 contain the predictor displacement
     CALL fieldadd(v1,u1,in%sx1+2,in%sx2,in%sx3/2,c1=REAL(Dt/2))
     CALL fieldadd(v2,u2,in%sx1+2,in%sx2,in%sx3/2,c1=REAL(Dt/2))
     CALL fieldadd(v3,u3,in%sx1+2,in%sx2,in%sx3/2,c1=REAL(Dt/2))
     CALL tensorfieldadd(sig,tau,in%sx1,in%sx2,in%sx3/2,c1=-REAL(Dt/2),c2=-1._4)

     ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     ! corrector
     ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     CALL stressupdate(v1,v2,v3,in%lambda,in%mu, &
                       in%dx1,in%dx2,in%dx3,in%sx1,in%sx2,in%sx3/2,sig)

     ! reinitialize moment density tensor
     CALL tensorfieldadd(moment,moment,in%sx1,in%sx2,in%sx3/2,0._4,0._4)
     
     IF (ALLOCATED(in%linearstruc)) THEN
        ! linear viscosity
        v1=0
        IF (0 .LT. in%nlwz) THEN
           CALL viscouseigenstress(in%mu,in%linearstruc, &
                sig,in%stressstruc,in%sx1,in%sx2,in%sx3/2, &
                in%dx1,in%dx2,in%dx3,moment,DGAMMADOT0=lineardgammadot0,GAMMA=v1)
        ELSE
           CALL viscouseigenstress(in%mu,in%linearstruc, &
                sig,in%stressstruc,in%sx1,in%sx2,in%sx3/2, &
                in%dx1,in%dx2,in%dx3,moment,GAMMA=v1)
        END IF
        
        ! update slip history
        CALL fieldadd(gamma,v1,in%sx1+2,in%sx2,in%sx3/2,c2=REAL(Dt))
     END IF
     
     IF (ALLOCATED(in%nonlinearstruc)) THEN
        ! powerlaw viscosity
        v1=0
        IF (0 .LT. in%nnlwz) THEN
           CALL viscouseigenstress(in%mu,in%nonlinearstruc, &
                sig,in%stressstruc,in%sx1,in%sx2,in%sx3/2, &
                in%dx1,in%dx2,in%dx3,moment,DGAMMADOT0=nonlineardgammadot0,GAMMA=v1)
        ELSE
           CALL viscouseigenstress(in%mu,in%nonlinearstruc, &
                sig,in%stressstruc,in%sx1,in%sx2,in%sx3/2, &
                in%dx1,in%dx2,in%dx3,moment,GAMMA=v1)
        END IF
        
        ! update slip history
        CALL fieldadd(gamma,v1,in%sx1+2,in%sx2,in%sx3/2,c2=REAL(Dt))
     END IF
     
     ! nonlinear fault creep with rate-strengthening friction
     IF (ALLOCATED(in%faultcreepstruc)) THEN

        ! use v1 as placeholders for the afterslip planes
        DO k=1,in%np
           ! one may use optional arguments ...,VEL=v1) to convert
           ! fault slip to eigenstrain (scalar)
           CALL frictioneigenstress(in%n(k)%x,in%n(k)%y,in%n(k)%z, &
                in%n(k)%width,in%n(k)%length, &
                in%n(k)%strike,in%n(k)%dip,in%n(k)%rake,in%beta, &
                sig,in%stressstruc,in%mu,in%faultcreepstruc,in%sx1,in%sx2,in%sx3/2, &
                in%dx1,in%dx2,in%dx3,moment)

           ! keep track if afterslip instantaneous velocity
           CALL monitorfriction(in%n(k)%x,in%n(k)%y,in%n(k)%z, &
                in%n(k)%width,in%n(k)%length, &
                in%n(k)%strike,in%n(k)%dip,in%n(k)%rake,in%beta, &
                in%sx1,in%sx2,in%sx3,in%dx1,in%dx2,in%dx3, &
                sig,in%stressstruc,in%faultcreepstruc,in%n(k)%patch)
        END DO

     END IF

     ! interseismic loading
     IF ((in%inter%ns .GT. 0) .OR. (in%inter%nt .GT. 0)) THEN
        ! vectors v1,v2,v3 are not affected.
        CALL dislocations(in%inter,in%lambda,in%mu,in%beta,in%sx1,in%sx2,in%sx3, &
             in%dx1,in%dx2,in%dx3,v1,v2,v3,t1,t2,t3,tau,eigenstress=moment)
     END IF
     
     v1=0;v2=0;v3=0;t1=0;t2=0;t3=0;
     CALL equivalentbodyforce(moment,in%dx1,in%dx2,in%dx3,in%sx1,in%sx2,in%sx3/2,v1,v2,v3,t1,t2,t3)

     ! add time-dependent surface loads
     CALL traction(in%mu,in%events(e),in%sx1,in%sx2,in%dx1,in%dx2,t,Dt,t3,rate=.TRUE.)

     CALL greenfunctioncowling(v1,v2,v3,t1,t2,t3,in%dx1,in%dx2,in%dx3,in%lambda,in%mu,in%gam)

     ! update deformation field
     CALL fieldadd(u1,v1,in%sx1+2,in%sx2,in%sx3/2,c2=REAL(Dt))
     CALL fieldadd(u2,v2,in%sx1+2,in%sx2,in%sx3/2,c2=REAL(Dt))
     CALL fieldadd(u3,v3,in%sx1+2,in%sx2,in%sx3/2,c2=REAL(Dt))
     CALL tensorfieldadd(tau,moment,in%sx1,in%sx2,in%sx3/2,c2=REAL(Dt))
     CALL frictionadd(in%np,in%n,Dt)
     
     ! time increment
     t=t+Dt
    
     ! next event
     IF (e .LT. in%ne) THEN
        IF (abs(t-in%events(e+1)%time) .LT. 1e-6) THEN
           e=e+1
           in%events(e)%i=i

           IF (isverbose) THEN
              PRINT '("coseismic event ",I3.3)', e
              PRINT 0990
           END IF

           v1=0;v2=0;v3=0;t1=0;t2=0;t3=0;
           CALL dislocations(in%events(e),in%lambda,in%mu, &
                in%beta,in%sx1,in%sx2,in%sx3, &
                in%dx1,in%dx2,in%dx3,v1,v2,v3,t1,t2,t3,tau)
           CALL traction(in%mu,in%events(e),in%sx1,in%sx2,in%dx1,in%dx2,t,0.d0,t3)

           ! apply the 3d elastic transfert function
           CALL greenfunctioncowling(v1,v2,v3,t1,t2,t3, &
                in%dx1,in%dx2,in%dx3,in%lambda,in%mu,in%gam)
           
           ! transfer solution
           CALL fieldadd(u1,v1,in%sx1+2,in%sx2,in%sx3/2)
           CALL fieldadd(u2,v2,in%sx1+2,in%sx2,in%sx3/2)
           CALL fieldadd(u3,v3,in%sx1+2,in%sx2,in%sx3/2)

        END IF
     END IF

     CALL tensorfieldadd(sig,tau,in%sx1,in%sx2,in%sx3/2,c1=0._4,c2=-1._4)
     CALL stressupdate(u1,u2,u3,in%lambda,in%mu, &
                       in%dx1,in%dx2,in%dx3,in%sx1,in%sx2,in%sx3/2,sig)

     ! export displacements
     CALL pts2series(u1,u2,u3,in%sx1,in%sx2,in%sx3,in%dx1,in%dx2,in%dx3,in%opts,t,1+i,gps)

     IF (ALLOCATED(in%ptsname)) THEN
        CALL exportpoints(u1,u2,u3,sig,in%sx1,in%sx2,in%sx3/2,in%dx1,in%dx2,in%dx3, &
                          in%opts,in%ptsname,t,in%wdir,.FALSE.,in%x0,in%y0,in%rot)
     END IF

     ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     ! -   export displacement and stress
     ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

     ! output only at discrete intervals (skip=0, odt>0),
     ! or every "skip" computational steps (skip>0, odt<0),
     ! or anytime a coseismic event occurs
     IF (isverbose) THEN
        IF (isoutput(in%skip,t,i,in%odt,oi,in%events(e)%time)) THEN
        
           WRITE (digit4,'(I4.4)') oi
           PRINT 1101,i,Dt,maxwell,t,in%interval, &
                tensoramplitude(moment,in%dx1,in%dx2,in%dx3), &
                tensoramplitude(tau,in%dx1,in%dx2,in%dx3)

           ! update output counter
           oi=oi+1
        ELSE
           PRINT 1100,i,Dt,maxwell,t,in%interval, &
                tensoramplitude(moment,in%dx1,in%dx2,in%dx3), &
                tensoramplitude(tau,in%dx1,in%dx2,in%dx3)
        END IF
     END IF

  END DO

100 CONTINUE

!  DO i=1,in%ne
!     IF (ALLOCATED(in%events(i)%s))  DEALLOCATE(in%events(i)%s,in%events(i)%sc)
!     IF (ALLOCATED(in%events(i)%ts)) DEALLOCATE(in%events(i)%ts,in%events(i)%tsc)
!  END DO
!  IF (ALLOCATED(in%events)) DEALLOCATE(in%events)

  ! free memory
  IF (ALLOCATED(gamma)) DEALLOCATE(gamma)
  IF (ALLOCATED(sig)) DEALLOCATE(sig)
  IF (ALLOCATED(tau)) DEALLOCATE(tau)
  IF (ALLOCATED(moment)) DEALLOCATE(moment)
  IF (ALLOCATED(v1)) DEALLOCATE(v1,v2,v3,t1,t2,t3)
  IF (ALLOCATED(u1)) DEALLOCATE(u1,u2,u3)
  IF (ALLOCATED(in%stressstruc)) DEALLOCATE(in%stressstruc)

#ifdef FFTW3_THREADS
  CALL sfftw_cleanup_threads()
#endif

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

END SUBROUTINE relaxlite 
