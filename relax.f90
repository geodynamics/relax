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

PROGRAM relax
  !-----------------------------------------------------------------------
  ! PURPOSE:
  !   The program RELAX computes nonlinear time-dependent viscoelastic
  !   deformation with powerlaw rheology and rate-strengthening friction 
  !   in a cubic, periodic grid due to coseismic stress changes, initial
  !   stress, surface loads, and/or moving faults.
  !
  ! DESCRIPTION:
  !   Computation is done semi-analytically inside a cartesian grid.
  !   The grid is defined by its size sx1*sx2*sx3 and the sampling
  !   intervals dx1, dx2 and dx3. rule of thumb is to allow for at least
  !   five samples per fault length or width, and to have the tip of any 
  !   fault at least 10 fault widths away from any edge of the 
  !   computational grid.
  !
  !   Coseismic stress changes and initial coseismic deformation results
  !   from the presence of dislocations in the brittle layer. Fault
  !   geometry is prescribed following Okada or Wang's convention, with the
  !   usual slip, strike, dip and rake and is converted to a double-couple
  !   equivalent body-force analytically. Current implementation allows 
  !   shear fault (strike slip and dip slip), dykes, Mogi source, and
  !   surface traction. Faults and dykes can be of arbitrary orientation 
  !   in the half space.
  !
  ! METHOD:
  !   The current implementation is organized to integrate stress/strain-
  !   rate constitutive laws (rheologies) of the form
  !
  !       epsilon^dot = f(sigma)
  !
  !   as opposed to epsilon^dot = f(sigma,epsilon) wich would include work-
  !   hardening (or weakening). The time-stepping implements a second-order
  !   Runge-Kutta numerical integration scheme with a variable time-step.
  !   The Runge-Kutta method integrating the ODE y'=f(x,y) can be summarized
  !   as follows:
  !
  !          y_(n+1) = y_n + k_2
  !              k_1 = h * f(x_n, y_n)
  !              k_2 = h * f(x_n + h, y_n + k_1)
  !
  !   where h is the time-step and n is the time-index. The elastic response
  !   in the computational grid is obtained using elastic Greens functions.
  !   The Greens functions are applied in the Fourier domain. Strain,
  !   stress and body-forces are obtained by application of a finite impulse
  !   response (FIR) differentiator filter in the space domain.
  !
  ! INPUT:
  !   Static dislocation sources are discretized into a series of planar
  !   segments. Slip patches are defined in terms of position, orientation,
  !   and slip, as illustrated in the following figure:
  !
  !                 N (x1)
  !                /
  !               /| Strike
  !   x1,x2,x3 ->@------------------------      (x2)
  !              |\        p .            \ W
  !              :-\      i .              \ i
  !              |  \    l .                \ d
  !              :90 \  S .                  \ t
  !              |-Dip\  .                    \ h
  !              :     \. | Rake               \
  !              |      -------------------------
  !              :             L e n g t h
  !              Z (x3)
  !
  !   Dislocations are converted to double-couple equivalent body-force
  !   analytically. Solution displacement is obtained by application of
  !   the Greens functions in the Fourier domain.
  !
  !   For friction faults where slip rates are evaluated from stress and
  !   a constitutive law, the rake corresponds to the orientation of slip. 
  !   That is, if r_i is the rake vector and v_i is the instantaneous 
  !   velocity vector, then r_j v_j >= 0. 
  !
  ! OUTPUT:
  !   The vector-valued deformation is computed everywhere in a cartesian
  !   grid. The vector field is sampled 1) along a horizontal surface at a
  !   specified depth and 2) at specific points. Format is always North (x1), 
  !   East (x2) and Down (x3) components, following the right-handed reference 
  !   system convention. North corresponds to x1-direction, East to the 
  !   x2-direction and down to the x3-direction. The Generic Mapping Tool 
  !   output files are labeled explicitely ???-north.grd, ???-east.grd and 
  !   ???-up.grd (or say, ???-geo-up.grd for outputs in geographic 
  !   coordinates), where ??? stands for an output index: 001, 002, ...
  !
  !   The amplitude of the inelastic (irreversible) deformation is also
  !   tracked and can be output along a plane of arbitrary orientation.
  !   The inelastic deformation includes the initial, constrained, slip on
  !   fault surfaces, the time-dependent slip on frictional surfaces and
  !   the cumulative amplitude of bulk strain in viscoelastic regions.
  !   Slip is provided as a function of local coordinates along strike and 
  !   dip as well as a function of the Cartesian coordinates for three-
  !   dimensional display.
  !
  !   Time integration uses adaptive time steps to ensure accuracy but
  !   results can be output either 1) at specified uniform time intervals 
  !   or 2) at the same intervals as computed. In the later case, output 
  !   intervals is chosen internally depending on instantaneous relaxation 
  !   rates.
  !
  ! TECHNICAL ASPECTS:
  !   Most of the computational burden comes from 1) applying the elastic
  !   Green function and 2) computing the current strain from a displacement
  !   field. The convolution of body forces with the Green function is 
  !   performed in the Fourier domain and the efficiency of the computation
  !   depends essentially upon a choice of the discrete Fourier transform.
  !   Current implementation is compatible with the Couley-Tuckey, the
  !   Fast Fourier transform of the West (FFTW), the SGI FFT and the intel
  !   FFT from the intel MKL library. Among these choices, the MKL FFT is
  !   the most efficient. The FFTW, SGI FFT and MKL FFT can all be ran
  !   in parallel on shared-memory computers.
  !
  !   Strain is computed using a Finite Impulse Response differentiator
  !   filter in the space domain. Use of FIR filter give rise to very
  !   accurate derivatives but is computationally expensive. The filter
  !   kernels are provided in the kernel???.inc files. Use of a compact
  !   kernel may accelerate computation significantly.
  !
  !   Compilation options are defined in the include.f90 file and specify
  !   for instance the choice of DFT and the kind of output provided.
  !
  ! MODIFICATIONS:
  !   sylvain barbot (07-06-07) - original form
  !                  (08-28-08) - FFTW/SGI_FFT support, FIR derivatives,
  !                               Runge-Kutta integration, tensile cracks,
  !                               GMT output, comments in input file
  !                  (10-24-08) - interseismic loading, postseismic signal
  !                               output in separate files
  !                  (12-08-09) - slip distribution smoothing
  !                  (05-05-10) - lateral variations in viscous properties
  !                               Intel MKL implementation of the FFT
  !                  (06-04-10) - output in geographic coordinates
  !                               and output components of stress tensor
  !                  (07-19-10) - includes surface tractions initial condition
  !                               output geometry in VTK format for Paraview
  !                  (02-28-11) - add constraints on the broad direction of 
  !                               afterslip, export faults to GMT xy format
  !                               and allow scaling of computed time steps.
  !-----------------------------------------------------------------------

  USE green
  USE elastic3d
  USE viscoelastic3d
  USE friction3d
  USE export

#include "include.f90"
  
  IMPLICIT NONE
  
  REAL*8, PARAMETER :: DEG2RAD = 0.01745329251994329547437168059786927_8
  INTEGER, PARAMETER :: ITERATION_MAX = 9900
  REAL*8, PARAMETER :: STEP_MAX = 1e7

  INTEGER :: i,k,sx1,sx2,sx3,e,ne,nv,np,nop,npl,nps,oi,nfc, &
       unit,iostatus,iargc,npts,skip=0,mech(3),nlwz,nnlwz
#ifdef FFTW3_THREADS
  INTEGER :: iret
!$  INTEGER :: omp_get_max_threads
#endif
  REAL*8 :: beta,lambda,mu,gam,x0,y0,interval, &
       minlength,minwidth,rot,maxwell(3),nyquist
#ifdef PROJ
  REAL*8 :: lon0,lat0,umult
  INTEGER :: zone
#endif
  CHARACTER(80) :: wdir,reporttimefilename,reportfilename, &
                   inputfile,logfilename,inputfilename
#ifdef VTK
  INTEGER :: j
  CHARACTER(80) :: rffilename,vcfilename,cgfilename
  CHARACTER(3) :: digit
#endif
  REAL*8 :: dx1,dx2,dx3,oz,ozs,t,Dt,tm,odt,tscale
  ! coseismic events
  TYPE(EVENT_STRUC), DIMENSION(:), ALLOCATABLE :: events
  TYPE(EVENT_STRUC) :: inter
  
  ! input dislocation (shear and tensile cracks)
  TYPE(PLANE_STRUCT), DIMENSION(:), ALLOCATABLE :: n, op
  TYPE(LAYER_STRUCT), DIMENSION(:), ALLOCATABLE :: linearlayer,nonlinearlayer
  TYPE(LAYER_STRUCT), DIMENSION(:), ALLOCATABLE :: faultcreeplayer
  TYPE(LAYER_STRUCT), DIMENSION(:), ALLOCATABLE :: linearstruc,nonlinearstruc
  TYPE(LAYER_STRUCT), DIMENSION(:), ALLOCATABLE :: faultcreepstruc
  TYPE(TENSOR_LAYER_STRUCT), DIMENSION(:), ALLOCATABLE :: stresslayer,stressstruc
  TYPE(WEAK_STRUCT), DIMENSION(:), ALLOCATABLE :: linearweakzone,linearweakzonec, &
                                            nonlinearweakzone,nonlinearweakzonec
  
  ! arrays
  REAL*4, DIMENSION(:,:,:), ALLOCATABLE :: v1,v2,v3,u1,u2,u3,gamma
  REAL*4, DIMENSION(:,:), ALLOCATABLE :: t1,t2,t3
  REAL*4, DIMENSION(:,:,:), ALLOCATABLE :: inter1,inter2,inter3
  TYPE(TENSOR), DIMENSION(:,:,:), ALLOCATABLE :: tau,sig,moment
  TYPE(VECTOR_STRUCT), DIMENSION(:), ALLOCATABLE :: opts
  CHARACTER(LEN=4), DIMENSION(:), ALLOCATABLE :: ptsname
  
#ifdef FFTW3_THREADS
  CALL sfftw_init_threads(iret)
#ifdef _OPENMP
  CALL sfftw_plan_with_nthreads(omp_get_max_threads())
#else
  CALL sfftw_plan_with_nthreads(4)
#endif
#endif

  ! read standard input or filename given in argument
  IF (0 .EQ. iargc()) THEN
     ! standard input
     unit=5
  ELSE
     ! open input file
     CALL getarg(1,inputfile)

     OPEN (UNIT=15,FILE=inputfile,IOSTAT=iostatus,FORM="FORMATTED")
     IF (iostatus .GT. 0) THEN
        WRITE_DEBUG_INFO
        WRITE (0,'("unable to access input file ",a)') inputfile
        STOP 1
     END IF
     ! input file
     unit=15
  END IF

  CALL init(UNIT=unit)

  ! close input file
  IF (iargc() .GT. 0) CLOSE(15)

  ALLOCATE (v1(sx1+2,sx2,sx3),v2(sx1+2,sx2,sx3),v3(sx1+2,sx2,sx3), &
            u1(sx1+2,sx2,sx3/2),u2(sx1+2,sx2,sx3/2),u3(sx1+2,sx2,sx3/2), &
            inter1(sx1+2,sx2,2),inter2(sx1+2,sx2,2),inter3(sx1+2,sx2,2), &
            tau(sx1,sx2,sx3/2),gamma(sx1+2,sx2,sx3/2), &
            t1(sx1+2,sx2),t2(sx1+2,sx2),t3(sx1+2,sx2), &
            STAT=iostatus)
  IF (iostatus>0) STOP "could not allocate memory"
  v1=0;v2=0;v3=0;u1=0;u2=0;u3=0;gamma=0;t1=0;t2=0;t3=0
  CALL tensorfieldadd(tau,tau,sx1,sx2,sx3/2,c1=0._4,c2=0._4)

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! -     construct pre-stress structure
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  IF (ALLOCATED(stresslayer)) THEN
     CALL tensorstructure(stressstruc,stresslayer,dx3)
     DEALLOCATE(stresslayer)
     
     DO k=1,sx3/2
        tau(:,:,k)=(-1._4) .times. stressstruc(k)%t
     END DO
     DEALLOCATE(stressstruc)
  END IF

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! -     construct linear viscoelastic structure
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  IF (ALLOCATED(linearlayer)) THEN
     CALL viscoelasticstructure(linearstruc,linearlayer,dx3)
     DEALLOCATE(linearlayer)
  END IF

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! -   construct nonlinear viscoelastic structure
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  IF (ALLOCATED(nonlinearlayer)) THEN
     CALL viscoelasticstructure(nonlinearstruc,nonlinearlayer,dx3)
     DEALLOCATE(nonlinearlayer)
  END IF

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! -   construct nonlinear fault creep structure (rate-strenghtening)
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  IF (ALLOCATED(faultcreeplayer)) THEN
     CALL viscoelasticstructure(faultcreepstruc,faultcreeplayer,dx3)
     DEALLOCATE(faultcreeplayer)
  END IF

  ! first event
  e=1
  ! first output
  oi=1;

  ! sources
  CALL dislocations(events(e),lambda,mu,beta,sx1,sx2,sx3, &
       dx1,dx2,dx3,v1,v2,v3,t1,t2,t3,tau)
  CALL traction(mu,events(e),sx1,sx2,dx1,dx2,t3)
  
  PRINT '("coseismic event ",I3.3)', e
  PRINT 0990

  ! export the amplitude of eigenstrain
  CALL exporteigenstrain(gamma,nop,op,x0,y0,dx1,dx2,dx3,sx1,sx2,sx3/2,wdir,0)
  
  ! export equivalent body forces
  IF (isoutput(skip,t,i,odt,oi,events(e)%time)) THEN
#ifdef GRD_EQBF
     CALL exportgrd(v1,v2,v3,sx1,sx2,sx3/2,dx1,dx2,dx3,0.7_8,x0,y0,wdir,0,convention=3)
#endif
  END IF

  ! test the presence of dislocations for coseismic calculation
  IF ((events(e)%nt .NE. 0) .OR. &
      (events(e)%ns .NE. 0) .OR. &
      (events(e)%nm .NE. 0)) THEN

     ! apply the 3d elastic transfer function
     CALL greenfunctioncowling(v1,v2,v3,t1,t2,t3,dx1,dx2,dx3,lambda,mu,gam)
  END IF
  
  ! transfer solution
  CALL fieldrep(u1,v1,sx1+2,sx2,sx3/2)
  CALL fieldrep(u2,v2,sx1+2,sx2,sx3/2)
  CALL fieldrep(u3,v3,sx1+2,sx2,sx3/2)

  ! export
#ifdef TXT
  CALL exporttxt(u1,u2,u3,sx1,sx2,sx3/2,oz,dx3,0,0._8,wdir,reportfilename)
#endif
#ifdef XYZ
  CALL exportxyz(u1,u2,u3,sx1,sx2,sx3/2,oz,dx1,dx2,dx3,0,wdir)
#endif
#ifdef GRD
  CALL exportgrd(u1,u2,u3,sx1,sx2,sx3/2,dx1,dx2,dx3,oz,x0,y0,wdir,0)
  CALL exportgrd(inter1,inter2,inter3,sx1,sx2,sx3/2, &
       dx1,dx2,dx3,0._8,x0,y0,wdir,0,convention=2)
#endif
#ifdef PROJ
  CALL exportproj(u1,u2,u3,sx1,sx2,sx3/2,dx1,dx2,dx3,oz, &
                  x0,y0,lon0,lat0,zone,umult,wdir,0)
#endif
#ifdef VTK
  j=INDEX(wdir," ")
  vcfilename=wdir(1:j-1)//"/disp-000.vtr"
  CALL exportvtk_vectors(u1,u2,u3,sx1,sx2,sx3/4,dx1,dx2,dx3,8,8,8,vcfilename)
  !CALL exportvtk_vectors_slice(u1,u2,u3,sx1,sx2,sx3/2,dx1,dx2,dx3,oz,8,8,vcfilename)
#endif
  IF (ALLOCATED(ptsname)) THEN
     CALL exportpoints(u1,u2,u3,sx1,sx2,sx3/2,dx1,dx2,dx3, &
          opts,ptsname,0._8,wdir,.true.,x0,y0,rot)
  END IF
  CALL reporttime(0,0._8,reporttimefilename)

  PRINT 1101,0,0._8,0._8,0._8,0._8,0._8,interval,0._8,tensoramplitude(tau,dx1,dx2,dx3)
  IF (interval .LE. 0) THEN
     GOTO 100 ! no time integration
  END IF

  ALLOCATE(moment(sx1,sx2,sx3/2),sig(sx1,sx2,sx3/2),STAT=iostatus)
  IF (iostatus>0) STOP "could not allocate the mechanical structure"

  CALL tensorfieldadd(sig,sig,sx1,sx2,sx3/2,c1=0._4,c2=0._4)
  CALL tensorfieldadd(moment,moment,sx1,sx2,sx3/2,c1=0._4,c2=0._4)  

  t=0
  DO i=1,ITERATION_MAX
     IF (t > (interval+1e-6)) GOTO 100 ! proper exit
     
     ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     ! predictor
     ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     CALL tensorfieldadd(sig,tau,sx1,sx2,sx3/2,c1=0._4,c2=-1._4)
     CALL stressupdate(u1,u2,u3,lambda,mu,dx1,dx2,dx3,sx1,sx2,sx3/2,sig)

     IF (isoutput(skip,t,i,odt,oi,events(e)%time)) THEN
        ! export stress
#ifdef GRD
        CALL exportstressgrd(sig,sx1,sx2,sx3/2,dx1,dx2,dx3,ozs,x0,y0,wdir,i-1)
#endif
#ifdef PROJ
        CALL exportstressproj(sig,sx1,sx2,sx3/2,dx1,dx2,dx3,ozs, &
                              x0,y0,lon0,lat0,zone,umult,wdir,i-1)
#endif
     END IF

     ! initialize large time step
     tm=STEP_MAX;
     maxwell(:)=STEP_MAX;
     
     ! active mechanism flag
     mech(:)=0

     ! initialize no forcing term in tensor space
     CALL tensorfieldadd(moment,moment,sx1,sx2,sx3/2,0._4,0._4)

     ! power density from three mechanisms (linear and power-law viscosity 
     ! and fault creep)
     ! 1- linear viscosity
     IF (ALLOCATED(linearstruc)) THEN
        CALL viscouseigenstress(mu,linearstruc,linearweakzone,sig,sx1,sx2,sx3/2, &
             dx1,dx2,dx3,moment,0.01_8,MAXWELLTIME=maxwell(1))
        mech(1)=1
     END IF
     
     ! 2- powerlaw viscosity
     IF (ALLOCATED(nonlinearstruc)) THEN
        CALL viscouseigenstress(mu,nonlinearstruc,nonlinearweakzone,sig,sx1,sx2,sx3/2, &
             dx1,dx2,dx3,moment,0.01_8,MAXWELLTIME=maxwell(2))
        mech(2)=1
     END IF
     
     ! 3- nonlinear fault creep with rate-strengthening friction
     IF (ALLOCATED(faultcreepstruc)) THEN
        DO k=1,np
           CALL frictioneigenstress(n(k)%x,n(k)%y,n(k)%z, &
                n(k)%width,n(k)%length,n(k)%strike,n(k)%dip,n(k)%rake,beta, &
                sig,mu,faultcreepstruc,sx1,sx2,sx3/2,dx1,dx2,dx3,moment, &
                maxwelltime=maxwell(3))
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
     IF ((inter%ns .GT. 0) .OR. (inter%nt .GT. 0)) THEN
        IF (tm .EQ. STEP_MAX) THEN
           ! no relaxation occurs, pick a small integration time
           tm=interval/20._8
        END IF
     END IF
     
     ! choose an integration time step
     CALL integrationstep(tm,Dt,t,oi,odt,skip,tscale,events,e,ne)

     CALL tensorfieldadd(sig,moment,sx1,sx2,sx3/2,c1=0.0_4,c2=1._4)
     
     v1=0;v2=0;v3=0;t1=0;t2=0;t3=0;
     CALL equivalentbodyforce(sig,dx1,dx2,dx3,sx1,sx2,sx3/2,v1,v2,v3,t1,t2,t3)
     CALL greenfunctioncowling(v1,v2,v3,t1,t2,t3,dx1,dx2,dx3,lambda,mu,gam)
     
     ! v1,v2,v3 contain the predictor displacement
     CALL fieldadd(v1,u1,sx1+2,sx2,sx3/2,c1=REAL(Dt/2))
     CALL fieldadd(v2,u2,sx1+2,sx2,sx3/2,c1=REAL(Dt/2))
     CALL fieldadd(v3,u3,sx1+2,sx2,sx3/2,c1=REAL(Dt/2))
     CALL tensorfieldadd(sig,tau,sx1,sx2,sx3/2,c1=-REAL(Dt/2),c2=-1._4)

     ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     ! corrector
     ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     CALL stressupdate(v1,v2,v3,lambda,mu,dx1,dx2,dx3,sx1,sx2,sx3/2,sig)

     ! reinitialize moment density tensor
     CALL tensorfieldadd(moment,moment,sx1,sx2,sx3/2,0._4,0._4)
     
     IF (ALLOCATED(linearstruc)) THEN
        ! linear viscosity
        v1=0
        CALL viscouseigenstress(mu,linearstruc,linearweakzone,sig,sx1,sx2,sx3/2, &
             dx1,dx2,dx3,moment,0.01_8,GAMMA=v1)
        
        ! update slip history
        CALL fieldadd(gamma,v1,sx1+2,sx2,sx3/2,c2=REAL(Dt))
     END IF
     
     IF (ALLOCATED(nonlinearstruc)) THEN
        ! powerlaw viscosity
        v1=0
        CALL viscouseigenstress(mu,nonlinearstruc,nonlinearweakzone,sig,sx1,sx2,sx3/2, &
             dx1,dx2,dx3,moment,0.01_8,GAMMA=v1)
        
        ! update slip history
        CALL fieldadd(gamma,v1,sx1+2,sx2,sx3/2,c2=REAL(Dt))
     END IF
     
     ! nonlinear fault creep with rate-strengthening friction
     IF (ALLOCATED(faultcreepstruc)) THEN
        ! use v1 as placeholders for the afterslip planes
        v1=0
        DO k=1,np
           ! one may use optional arguments ...,VEL=v1) to convert
           ! fault slip to eigenstrain (scalar)
           CALL frictioneigenstress(n(k)%x,n(k)%y,n(k)%z, &
                n(k)%width,n(k)%length,n(k)%strike,n(k)%dip,n(k)%rake,beta, &
                sig,mu,faultcreepstruc,sx1,sx2,sx3/2,dx1,dx2,dx3,moment)
        END DO
        
        ! update slip history
        CALL fieldadd(gamma,v1,sx1+2,sx2,sx3/2,c2=REAL(Dt))

        ! export strike and dip creep velocity
        IF (isoutput(skip,t,i,odt,oi,events(e)%time)) THEN
           CALL exportcreep(np,n,beta,sig,faultcreepstruc, &
                            sx1,sx2,sx3/2,dx1,dx2,dx3,x0,y0,wdir,oi)
        END IF
     END IF
     
     ! interseismic loading
     IF ((inter%ns .GT. 0) .OR. (inter%nt .GT. 0)) THEN
        ! vectors v1,v2,v3 are not affected.
        CALL dislocations(inter,lambda,mu,beta,sx1,sx2,sx3, &
             dx1,dx2,dx3,v1,v2,v3,t1,t2,t3,tau,factor=Dt,eigenstress=moment)
     END IF
     
     v1=0;v2=0;v3=0;t1=0;t2=0;t3=0;
     CALL equivalentbodyforce(moment,dx1,dx2,dx3,sx1,sx2,sx3/2,v1,v2,v3,t1,t2,t3)

     ! export equivalent body forces
     IF (isoutput(skip,t,i,odt,oi,events(e)%time)) THEN
#ifdef VTK_EQBF
        WRITE (digit,'(I3.3)') oi
        j=INDEX(wdir," ")
        vcfilename=wdir(1:j-1)//"/eqbf-"//digit//".vtr"
        CALL exportvtk_vectors(v1,v2,v3,sx1,sx2,sx3/4,dx1,dx2,dx3,8,8,8,vcfilename)
#endif
#ifdef GRD_EQBF
        CALL exportgrd(v1,v2,v3,sx1,sx2,sx3/2,dx1,dx2,dx3,30.7_8,x0,y0,wdir,oi,convention=3)
#endif
     END IF

     CALL greenfunctioncowling(v1,v2,v3,t1,t2,t3,dx1,dx2,dx3,lambda,mu,gam)

     ! update deformation field
     CALL fieldadd(u1,v1,sx1+2,sx2,sx3/2,c2=REAL(Dt))
     CALL fieldadd(u2,v2,sx1+2,sx2,sx3/2,c2=REAL(Dt))
     CALL fieldadd(u3,v3,sx1+2,sx2,sx3/2,c2=REAL(Dt))
     CALL tensorfieldadd(tau,moment,sx1,sx2,sx3/2,c2=REAL(Dt))
     
     ! keep track of the viscoelastic contribution alone
     CALL sliceadd(inter1(:,:,1),v1,sx1+2,sx2,sx3,int(oz/dx3)+1,c2=REAL(Dt))
     CALL sliceadd(inter2(:,:,1),v2,sx1+2,sx2,sx3,int(oz/dx3)+1,c2=REAL(Dt))
     CALL sliceadd(inter3(:,:,1),v3,sx1+2,sx2,sx3,int(oz/dx3)+1,c2=REAL(Dt))

     ! time increment
     t=t+Dt
     
     ! next event
     IF (e .LT. ne) THEN
        IF (abs(t-events(e+1)%time) .LT. 1e-6) THEN
           e=e+1
           events(e)%i=i
           PRINT '("coseismic event ",I3.3)', e
           PRINT 0990
           
           v1=0;v2=0;v3=0;t1=0;t2=0;t3=0;
           CALL dislocations(events(e),lambda,mu,beta,sx1,sx2,sx3, &
                dx1,dx2,dx3,v1,v2,v3,t1,t2,t3,tau)
           CALL traction(mu,events(e),sx1,sx2,dx1,dx2,t3)

           ! apply the 3d elastic transfert function
           CALL greenfunctioncowling(v1,v2,v3,t1,t2,t3,dx1,dx2,dx3,lambda,mu,gam)
           
           ! transfer solution
           CALL fieldadd(u1,v1,sx1+2,sx2,sx3/2)
           CALL fieldadd(u2,v2,sx1+2,sx2,sx3/2)
           CALL fieldadd(u3,v3,sx1+2,sx2,sx3/2)

        END IF
     END IF

     ! points are exported at all time steps
     IF (ALLOCATED(ptsname)) THEN
        CALL exportpoints(u1,u2,u3,sx1,sx2,sx3/2,dx1,dx2,dx3, &
             opts,ptsname,t,wdir,.false.,x0,y0,rot)
     END IF

     ! output only at discrete intervals (skip=0, odt>0),
     ! or every "skip" computational steps (skip>0, odt<0),
     ! or anytime a coseismic event occurs
     IF (isoutput(skip,t,i,odt,oi,events(e)%time)) THEN
        
        CALL reporttime(1,t,reporttimefilename)

        ! export
#ifdef TXT
        CALL exporttxt(u1,u2,u3,sx1,sx2,sx3/2,oz,dx3,oi,t,wdir,reportfilename)
#endif  
#ifdef XYZ
        CALL exportxyz(u1,u2,u3,sx1,sx2,sx3/2,oz,dx1,dx2,dx3,i,wdir)
        !CALL exportxyz(inter1,inter2,inter3,sx1,sx2,sx3/2,0.0_8,dx1,dx2,dx3,i,wdir)
#endif
        CALL exporteigenstrain(gamma,nop,op,x0,y0,dx1,dx2,dx3,sx1,sx2,sx3/2,wdir,oi)
#ifdef GRD
        CALL exportgrd(inter1,inter2,inter3,sx1,sx2,sx3/2, &
                       dx1,dx2,dx3,0._8,x0,y0,wdir,oi,convention=2)
        CALL exportgrd(u1,u2,u3,sx1,sx2,sx3/2,dx1,dx2,dx3,oz,x0,y0,wdir,oi)
#endif
#ifdef PROJ
        CALL exportproj(inter1,inter2,inter3,sx1,sx2,sx3/2, &
                        dx1,dx2,dx3,oz,x0,y0, &
                        lon0,lat0,zone,umult,wdir,oi,convention=2)
        CALL exportproj(u1,u2,u3,sx1,sx2,sx3/2,dx1,dx2,dx3,oz,x0,y0, &
                        lon0,lat0,zone,umult,wdir,oi)
#endif
#ifdef VTK
        WRITE (digit,'(I3.3)') oi
        j=INDEX(wdir," ")
        ! export total displacement in VTK XML format
        vcfilename=wdir(1:j-1)//"/disp-"//digit//".vtr"
        CALL exportvtk_vectors(u1,u2,u3,sx1,sx2,sx3/4,dx1,dx2,dx3,8,8,8,vcfilename)
        !CALL exportvtk_vectors_slice(u1,u2,u3,sx1,sx2,sx3/2,dx1,dx2,dx3,oz,8,8,vcfilename)

        ! export instantaneous velocity in VTK XML format
        vcfilename=wdir(1:j-1)//"/vel-"//digit//".vtr"
        CALL exportvtk_vectors(v1,v2,v3,sx1,sx2,sx3/4,dx1,dx2,dx3,8,8,8,vcfilename)
        !CALL exportvtk_vectors_slice(v1,v2,v3,sx1,sx2,sx3/2,dx1,dx2,dx3,oz,8,8,vcfilename)
#endif

        PRINT 1101,i,Dt,maxwell,t,interval, &
             tensoramplitude(moment,dx1,dx2,dx3), &
             tensoramplitude(tau,dx1,dx2,dx3)

        ! update output counter
        oi=oi+1
     ELSE
        PRINT 1100,i,Dt,maxwell,t,interval, &
             tensoramplitude(moment,dx1,dx2,dx3), &
             tensoramplitude(tau,dx1,dx2,dx3)
     END IF

  END DO

100 CONTINUE

  DO i=1,ne
     IF (ALLOCATED(events(i)%s))  DEALLOCATE(events(i)%s,events(i)%sc)
     IF (ALLOCATED(events(i)%ts)) DEALLOCATE(events(i)%ts,events(i)%tsc)
  END DO
  IF (ALLOCATED(events)) DEALLOCATE(events)

  ! free memory
  IF (ALLOCATED(gamma)) DEALLOCATE(gamma)
  IF (ALLOCATED(opts)) DEALLOCATE(opts)
  IF (ALLOCATED(op)) DEALLOCATE(op)
  IF (ALLOCATED(n)) DEALLOCATE(n)
  IF (ALLOCATED(stressstruc)) DEALLOCATE(stressstruc)
  IF (ALLOCATED(linearstruc)) DEALLOCATE(linearstruc)
  IF (ALLOCATED(nonlinearstruc)) DEALLOCATE(nonlinearstruc)
  IF (ALLOCATED(faultcreepstruc)) DEALLOCATE(faultcreepstruc)
  IF (ALLOCATED(sig)) DEALLOCATE(sig)
  IF (ALLOCATED(tau)) DEALLOCATE(tau)
  IF (ALLOCATED(moment)) DEALLOCATE(moment)
  IF (ALLOCATED(stresslayer)) DEALLOCATE(stresslayer)
  IF (ALLOCATED(linearlayer)) DEALLOCATE(linearlayer)
  IF (ALLOCATED(nonlinearlayer)) DEALLOCATE(nonlinearlayer)
  IF (ALLOCATED(faultcreeplayer)) DEALLOCATE(faultcreeplayer)
  DEALLOCATE(v1,v2,v3,t1,t2,t3)
  DEALLOCATE(u1,u2,u3)
  DEALLOCATE(inter1,inter2,inter3)


#ifdef FFTW3_THREADS
  CALL sfftw_cleanup_threads()
#endif

0990 FORMAT (" I  |   Dt   | tm(ve) | tm(pl) | tm(as) |     t/tmax     | power  |  C:E^i | ")
1000 FORMAT (I3.3,"*",ES9.2E2,"                            ",ES9.2E2,"/",ES7.2E1)
1100 FORMAT (I3.3," ",ES9.2E2,3ES9.2E2,ES9.2E2,"/",ES7.2E1,2ES9.2E2)
1101 FORMAT (I3.3,"*",ES9.2E2,3ES9.2E2,ES9.2E2,"/",ES7.2E1,2ES9.2E2)
1200 FORMAT ("----------------------------------------------------------------------------")

CONTAINS

  !--------------------------------------------------------------------
  ! subroutine eqbf_mask
  ! fills an array with positive values if some linear/nonlinear/creep
  ! is expected at the corresponding depth, zero otherwise.
  !
  ! the mask can be given to the routine "equivalentBodyForce" to skip
  ! these depths where no creep happens.
  !--------------------------------------------------------------------
  SUBROUTINE eqbf_mask(mask,sx)
    INTEGER, INTENT(IN) :: sx
    REAL*4, DIMENSION(sx), INTENT(OUT) :: mask
    
    IF (ALLOCATED(linearstruc)) THEN
       DO k=1,sx
          mask(k)=MAX(mask(k),REAL(linearstruc(k)%gammadot0,4))
       END DO
    END IF
    IF (ALLOCATED(nonlinearstruc)) THEN
       DO k=1,sx
          mask(k)=MAX(mask(k),REAL(nonlinearstruc(k)%gammadot0,4))
       END DO
    END IF
    IF (ALLOCATED(faultcreepstruc)) THEN
       DO k=1,sx
          mask(k)=MAX(mask(k),REAL(faultcreepstruc(k)%gammadot0,4))
       END DO
    END IF

    ! smooth the mask in the depth direction
    mask(1:sx-2)=(mask(1:sx-2)+mask(2:sx-1)+mask(3:sx))/3._4

  END SUBROUTINE eqbf_mask

  !---------------------------------------------------------------------
  ! subroutine Traction 
  ! assigns the traction vector at the surface
  !
  ! sylvain barbot (07-19-07) - original form
  !---------------------------------------------------------------------
  SUBROUTINE traction(mu,e,sx1,sx2,dx1,dx2,t3)
    TYPE(EVENT_STRUC), INTENT(IN) :: e
    INTEGER, INTENT(IN) :: sx1,sx2
    REAL*8, INTENT(IN) :: mu,dx1,dx2
#ifdef ALIGN_DATA
    REAL*4, DIMENSION(sx1+2,sx2), INTENT(INOUT) :: t3
#else
    REAL*4, DIMENSION(sx1,sx2), INTENT(INOUT) :: t3
#endif

    INTEGER :: i1,i2,i3

    DO i=1,e%nl
       CALL shiftedindex(e%l(i)%x,e%l(i)%y,0._8,sx1,sx2,1,dx1,dx2,1._8,i1,i2,i3)

       ! surface tractions
       t3(i1,i2)=t3(i1,i2)-e%l(i)%slip*mu
    END DO
             
  END SUBROUTINE traction

  !--------------------------------------------------------------------
  ! subroutine dislocations
  ! assigns equivalent body forces or moment density to simulate
  ! shear dislocations and fault opening. add the corresponding moment
  ! density in the cumulative relaxed moment so that fault slip does
  ! not reverse in the postseismic time.
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
    REAL*8 :: slip_factor=1._8
    
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
          ! adding sources in the space domain
          CALL source(mu,slip_factor*event%s(i)%slip, &
               event%s(i)%x,event%s(i)%y,event%s(i)%z, &
               event%s(i)%width,event%s(i)%length, &
               event%s(i)%strike,event%s(i)%dip,event%s(i)%rake, &
               beta,sx1,sx2,sx3,dx1,dx2,dx3,v1,v2,v3,t1,t2,t3)
       END DO
    ELSE
       ! forcing term in moment density
       DO i=1,event%ns
          CALL momentdensityshear(mu,slip_factor*event%s(i)%slip, &
               event%s(i)%x,event%s(i)%y,event%s(i)%z, &
               event%s(i)%width,event%s(i)%length, &
               event%s(i)%strike,event%s(i)%dip,event%s(i)%rake, &
               beta,sx1,sx2,sx3/2,dx1,dx2,dx3,eigenstress)
       END DO
    END IF

    DO i=1,event%ns
       ! remove corresponding eigenmoment
       CALL momentdensityshear(mu,slip_factor*event%s(i)%slip, &
            event%s(i)%x,event%s(i)%y,event%s(i)%z, &
            event%s(i)%width,event%s(i)%length, &
            event%s(i)%strike,event%s(i)%dip,event%s(i)%rake, &
            beta,sx1,sx2,sx3/2,dx1,dx2,dx3,tau)
    END DO
    
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! -             load tensile cracks
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    IF (.NOT. (PRESENT(eigenstress))) THEN
       ! forcing term in equivalent body force
       DO i=1,event%nt
          ! adding sources in the space domain
          CALL tensilesource(lambda,mu,slip_factor*event%ts(i)%slip, &
               event%ts(i)%x,event%ts(i)%y,event%ts(i)%z, &
               event%ts(i)%width,event%ts(i)%length, &
               event%ts(i)%strike,event%ts(i)%dip, &
               beta,sx1,sx2,sx3,dx1,dx2,dx3,v1,v2,v3)
       END DO
    ELSE
       ! forcing term in moment density
       DO i=1,event%nt
          CALL momentdensitytensile(lambda,mu,slip_factor*event%ts(i)%slip, &
               event%ts(i)%x,event%ts(i)%y,event%ts(i)%z,&
               event%ts(i)%width,event%ts(i)%length, &
               event%ts(i)%strike,event%ts(i)%dip,event%ts(i)%rake, &
               beta,sx1,sx2,sx3/2,dx1,dx2,dx3,eigenstress)
       END DO
    END IF

    DO i=1,event%nt
       ! removing corresponding eigenmoment
       CALL momentdensitytensile(lambda,mu,slip_factor*event%ts(i)%slip, &
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

  SUBROUTINE init(unit)
    INTEGER, OPTIONAL, INTENT(INOUT) :: unit

    INTEGER :: k,iostatus,i,e
    CHARACTER(180) :: dataline
#ifdef VTK
    INTEGER :: j
    CHARACTER(3) :: digit
#endif
    INTEGER :: iunit
!$  INTEGER :: omp_get_num_procs,omp_get_max_threads
    REAL*8 :: dummy,dum1,dum2

    ! default is standard input
    IF (.NOT. PRESENT(unit)) THEN
       iunit=5
    ELSE
       iunit=unit
    END IF

    PRINT 2000
    PRINT '("     nonlinear viscoelastic postseismic relaxation")'
#ifdef FFTW3
#ifdef FFTW3_THREADS
    PRINT '("     * FFTW3 (multi-threaded) implementation of the FFT")'
#else
    PRINT '("     * FFTW3 implementation of the FFT")'
#endif
#else
#ifdef SGI_FFT
    PRINT '("     * SGI_FFT implementation of the FFT")'
#else
#ifdef IMKL_FFT
    PRINT '("     * Intel MKL implementation of the FFT")'
#else
    PRINT '("     * fourt implementation of the FFT")'
#endif
#endif
#endif
!$  PRINT '("     * parallel OpenMP implementation with ",I3.3,"/",I3.3," threads")', &
!$                  omp_get_max_threads(),omp_get_num_procs()
#ifdef GRD
    PRINT '("     * export to GRD format")'
#endif
#ifdef TXT
    PRINT '("     * export to TXT format")'
#endif
#ifdef VTK
    PRINT '("     * export to VTK format")'
#endif
#ifdef PROJ
    PRINT '("     * export to longitude/latitude text format")'
#endif
    PRINT 2000

    PRINT '(a)', "grid dimension (sx1,sx2,sx3)"
    CALL getdata(iunit,dataline)
    READ (dataline,*) sx1,sx2,sx3
    PRINT '(3I5)', sx1,sx2,sx3

    PRINT '(a)', "sampling (dx1,dx2,dx3), smoothing (beta, nyquist)"
    CALL getdata(iunit,dataline)
    READ  (dataline,*) dx1,dx2,dx3,beta,nyquist
    PRINT '(5ES9.2E1)', dx1,dx2,dx3,beta,nyquist

    PRINT '(a)', "origin position (x0,y0) and rotation"
    CALL getdata(iunit,dataline)
    READ  (dataline,*) x0, y0, rot
    PRINT '(3ES9.2E1)', x0, y0, rot

#ifdef PROJ
    PRINT '(a)', "geographic origin (longitude, latitude, UTM zone, unit)"
    CALL getdata(iunit,dataline)
    READ  (dataline,*) lon0,lat0,zone,umult
    PRINT '(2ES9.2E1,I3.2,ES9.2E1)',lon0,lat0,zone,umult
    IF (zone.GT.60 .OR. zone.LT.1) THEN
       WRITE_DEBUG_INFO
       WRITE (0,'("invalid UTM zone ",I," (1<=zone<=60. exiting.)")') zone
       STOP 1
    END IF
#endif

    PRINT '(a)', "observation depth (displacement and stress)"
    CALL getdata(iunit,dataline)
    READ  (dataline,*) oz,ozs
    PRINT '(2ES9.2E1)', oz,ozs

    PRINT '(a)', "output directory"
    CALL getdata(iunit,dataline)
    READ (dataline,'(a)') wdir
    i=INDEX(wdir," ")
    reporttimefilename=wdir(1:i-1)//"/time.txt"
    reportfilename=wdir(1:i-1)//"/report.txt"
    logfilename=wdir(1:i-1)//"/relax.log"
    inputfilename=wdir(1:i-1)//"/relax.inp"
#ifdef TXT
    PRINT '(" ",a," (report: ",a,")")', wdir(1:i-1),reportfilename(1:i+10)
#else
    PRINT '(" ",a," (time report: ",a,")")', wdir(1:i-1),reporttimefilename(1:i+8)
#endif

    ! test write permissions on output directory
    OPEN (UNIT=14,FILE=reportfilename,POSITION="APPEND",&
            IOSTAT=iostatus,FORM="FORMATTED")
    IF (iostatus>0) THEN
       WRITE_DEBUG_INFO
       WRITE (0,'("unable to access ",a)') reporttimefilename(1:i+10)
       STOP 1
    END IF
    CLOSE(14)
    ! end test

#ifdef VTK
    cgfilename=wdir(1:i-1)//"/cgrid.vtp"
    CALL exportvtk_grid(sx1,sx2,sx3,dx1,dx2,dx3,x0,y0,cgfilename)
#endif

    PRINT '(a)', "lambda, mu, gamma (gamma = (1 - nu) rho g / mu)"
    CALL getdata(iunit,dataline)
    READ (dataline,*) lambda,mu,gam
    PRINT '(3ES10.2E2)',lambda,mu,gam

    PRINT '(a)', "time interval, (positive time step) or (negative skip, scaling)"
    CALL getdata(unit,dataline)
    READ  (dataline,*) interval, odt
    IF (odt .LT. 0.) THEN
       READ  (dataline,*) dum1, dum2, tscale
       skip=ceiling(-odt)
       PRINT '(ES9.2E1," (output every ",I3.3," steps, dt scaled by ",ES7.2E1,")")', &
             interval,skip,tscale
    ELSE
       PRINT '(ES9.2E1," (output every ",ES9.2E1," time unit)")', interval,odt
    END IF

    
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !         O B S E R V A T I O N       P L A N E S
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    PRINT '(a)', "number of observation planes"
    CALL getdata(unit,dataline)
    READ  (dataline,*) nop
    PRINT '(I5)', nop
    IF (nop .gt. 0) THEN
       ALLOCATE(op(nop),STAT=iostatus)
       IF (iostatus>0) STOP "could not allocate the observation plane list"
       PRINT 2000
       PRINT 2100
       PRINT 2000
       DO k=1,nop
          CALL getdata(unit,dataline)
          READ  (dataline,*) i,op(k)%x,op(k)%y,op(k)%z,&
               op(k)%length,op(k)%width,op(k)%strike,op(k)%dip

          PRINT '(I3.3," ",5ES9.2E1,2f7.1)', &
               k,op(k)%x,op(k)%y,op(k)%z, &
               op(k)%length,op(k)%width,op(k)%strike,op(k)%dip

          IF (i .ne. k) THEN
             WRITE_DEBUG_INFO
             WRITE (0,*) "error in input file: plane index misfit", k,"<>",i
             WRITE (0,*) op(k)
             STOP 1
          END IF

          ! comply to Wang's convention
          CALL wangconvention(dummy,op(k)%x,op(k)%y,op(k)%z,&
               op(k)%length,op(k)%width,op(k)%strike,op(k)%dip,dummy,rot)

       END DO
    END IF


    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !         O B S E R V A T I O N       P O I N T S
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    PRINT '(a)', "number of observation points"
    CALL getdata(iunit,dataline)
    READ  (dataline,*) npts
    PRINT '(I5)', npts
    IF (npts .gt. 0) THEN
       ALLOCATE(opts(npts),STAT=iostatus)
       IF (iostatus>0) STOP "could not allocate the observation point list"
       ALLOCATE(ptsname(npts),STAT=iostatus)
       IF (iostatus>0) STOP "could not allocate the list of point name"

       PRINT 2000
       PRINT 2300
       PRINT 2000
       DO k=1,npts
          CALL getdata(iunit,dataline)
          READ (dataline,*) i,ptsname(k),opts(k)%v1,opts(k)%v2,opts(k)%v3

          PRINT '(I3.3," ",A4,3ES9.2E1)', i,ptsname(k), &
               opts(k)%v1,opts(k)%v2,opts(k)%v3

          ! shift and rotate coordinates
          opts(k)%v1=opts(k)%v1-x0
          opts(k)%v2=opts(k)%v2-y0
          CALL rotation(opts(k)%v1,opts(k)%v2,rot)

          IF (i .ne. k) THEN
             WRITE_DEBUG_INFO
             WRITE (0,'("error in input file: points index misfit")')
             STOP 1
          END IF
       END DO

    END IF

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !                     P R E S T R E S S
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    PRINT '(a)', "number of prestress interfaces"
    CALL getdata(unit,dataline)
    READ  (dataline,*) nps
    PRINT '(I5)', nps

    IF (nps .GT. 0) THEN
       ALLOCATE(stresslayer(nps),stressstruc(sx3/2),STAT=iostatus)
       IF (iostatus>0) STOP "could not allocate the stress layer structure"
       
       PRINT 2000
       PRINT '(a)', "no.    depth  sigma11  sigma12  sigma13  sigma22  sigma23  sigma33"
       PRINT 2000
       DO k=1,nps
          CALL getdata(unit,dataline)
          READ  (dataline,*) i,stresslayer(k)%z, &
               stresslayer(k)%t%s11, stresslayer(k)%t%s12, &
               stresslayer(k)%t%s13, stresslayer(k)%t%s22, &
               stresslayer(k)%t%s23, stresslayer(k)%t%s33
          
          PRINT '(I3.3,7ES9.2E1)', i, &
               stresslayer(k)%z, &
               stresslayer(k)%t%s11, stresslayer(k)%t%s12, &
               stresslayer(k)%t%s13, stresslayer(k)%t%s22, &
               stresslayer(k)%t%s23, stresslayer(k)%t%s33
          
          IF (i .ne. k) THEN
             WRITE_DEBUG_INFO
             WRITE (0,'("error in input file: index misfit")')
             STOP 1
          END IF
       END DO
    END IF



    ! - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  L I N E A R    V I S C O U S    I N T E R F A C E
    ! - - - - - - - - - - - - - - - - - - - - - - - - - -
    PRINT '(a)', "number of linear viscous interfaces"
    CALL getdata(unit,dataline)
    READ  (dataline,*) nv
    PRINT '(I5)', nv
    
    IF (nv .GT. 0) THEN
       ALLOCATE(linearlayer(nv),linearstruc(sx3/2),STAT=iostatus)
       IF (iostatus>0) STOP "could not allocate the layer structure"
       
       PRINT 2000
       PRINT '(a)', "no.     depth    gamma0  cohesion"
       PRINT 2000
       DO k=1,nv
          CALL getdata(unit,dataline)
          READ  (dataline,*) i,linearlayer(k)%z, &
               linearlayer(k)%gammadot0, linearlayer(k)%cohesion

          linearlayer(k)%stressexponent=1

          PRINT '(I3.3,3ES10.2E2)', i, &
               linearlayer(k)%z, &
               linearlayer(k)%gammadot0, &
               linearlayer(k)%cohesion
          
          ! check positive strain rates
          IF (linearlayer(k)%gammadot0 .LT. 0) THEN
             WRITE_DEBUG_INFO
             WRITE (0,'("error in input file: strain rates must be positive")')
             STOP 1
          END IF

          IF (i .ne. k) THEN
             WRITE_DEBUG_INFO
             WRITE (0,'("error in input file: index misfit")')
             STOP 1
          END IF
#ifdef VTK
          ! export the viscous layer in VTK format
          j=INDEX(wdir," ")
          WRITE (digit,'(I3.3)') k

          rffilename=wdir(1:j-1)//"/linearlayer-"//digit//".vtp"
          CALL exportvtk_rectangle(0.d0,0.d0,linearlayer(k)%z, &
                                   DBLE(sx1)*dx1,DBLE(sx2)*dx2, &
                                   0._8,1.57d0,rffilename)
#endif
       END DO

       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !                 L I N E A R   W E A K   Z O N E S
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       PRINT '(a)', "number of linear weak zones (nlwz)"
       CALL getdata(iunit,dataline)
       READ  (dataline,*) nlwz
       PRINT '(I5)', nlwz
       IF (nlwz .GT. 0) THEN
          ALLOCATE(linearweakzone(nlwz),linearweakzonec(nlwz),STAT=iostatus)
          IF (iostatus>0) STOP "could not allocate the linear weak zones"
          PRINT 2000
          PRINT '(a)', "no. dgammadot0     x1       x2       x3  length   width thickn. strike   dip"
          PRINT 2000
          DO k=1,nlwz
             CALL getdata(iunit,dataline)
             READ  (dataline,*) i, &
                  linearweakzone(k)%dgammadot0, &
                  linearweakzone(k)%x,linearweakzone(k)%y,linearweakzone(k)%z,&
                  linearweakzone(k)%length,linearweakzone(k)%width,linearweakzone(k)%thickness, &
                  linearweakzone(k)%strike,linearweakzone(k)%dip
          
             linearweakzonec(k)=linearweakzone(k)
             
             PRINT '(I3.3,4ES9.2E1,3ES8.2E1,f7.1,f6.1)',k,&
                  linearweakzone(k)%dgammadot0, &
                  linearweakzone(k)%x,linearweakzone(k)%y,linearweakzone(k)%z, &
                  linearweakzone(k)%length,linearweakzone(k)%width, &
                  linearweakzone(k)%thickness, &
                  linearweakzone(k)%strike,linearweakzone(k)%dip
             
             IF (i .ne. k) THEN
                WRITE_DEBUG_INFO
                WRITE (0,'("error in input file: source index misfit")')
                STOP 1
             END IF
             ! comply to Wang's convention
             CALL wangconvention( &
                  dummy, & 
                  linearweakzone(k)%x,linearweakzone(k)%y,linearweakzone(k)%z, &
                  linearweakzone(k)%length,linearweakzone(k)%width, &
                  linearweakzone(k)%strike,linearweakzone(k)%dip,dummy,rot)
#ifdef VTK
                  ! export the ductile zone in VTK format
                  j=INDEX(wdir," ")-1
                  WRITE (digit,'(I3.3)') k

                  rffilename=wdir(1:j)//"/weakzone-"//digit//".vtp"
                  CALL exportvtk_brick(linearweakzone(k)%x,linearweakzone(k)%y,linearweakzone(k)%z, &
                                       linearweakzone(k)%length,linearweakzone(k)%width,linearweakzone(k)%thickness, &
                                       linearweakzone(k)%strike,linearweakzone(k)%dip,rffilename)
#endif
          END DO
       END IF
    END IF ! end linear viscous
       
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  N O N L I N E A R    V I S C O U S    I N T E R F A C E
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    PRINT '(a)', "number of nonlinear viscous interfaces"
    CALL getdata(unit,dataline)
    READ  (dataline,*) npl
    PRINT '(I5)', npl

    IF (npl .GT. 0) THEN
       ALLOCATE(nonlinearlayer(npl),nonlinearstruc(sx3/2),STAT=iostatus)
       IF (iostatus>0) STOP "could not allocate the layer structure"
       
       PRINT 2000
       PRINT '(a)', "no.     depth    gamma0     power  cohesion"
       PRINT 2000
       DO k=1,npl
          CALL getdata(unit,dataline)

          READ  (dataline,*) i,nonlinearlayer(k)%z, &
               nonlinearlayer(k)%gammadot0, &
               nonlinearlayer(k)%stressexponent, &
               nonlinearlayer(k)%cohesion

          PRINT '(I3.3,4ES10.2E2)', i, &
               nonlinearlayer(k)%z, &
               nonlinearlayer(k)%gammadot0, &
               nonlinearlayer(k)%stressexponent, &
               nonlinearlayer(k)%cohesion
          
          ! check positive strain rates
          IF (nonlinearlayer(k)%gammadot0 .LT. 0) THEN
             WRITE_DEBUG_INFO
             WRITE (0,'("error in input file: strain rates must be positive")')
             STOP 1
          END IF

          IF (i .ne. k) THEN
             WRITE_DEBUG_INFO
             WRITE (0,'("error in input file: index misfit")')
             STOP 1
          END IF
          
       END DO

       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !           N O N L I N E A R   W E A K   Z O N E S
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       PRINT '(a)', "number of nonlinear weak zones (nnlwz)"
       CALL getdata(iunit,dataline)
       READ  (dataline,*) nnlwz
       PRINT '(I5)', nnlwz
       IF (nnlwz .GT. 0) THEN
          ALLOCATE(nonlinearweakzone(nnlwz),nonlinearweakzonec(nnlwz),STAT=iostatus)
          IF (iostatus>0) STOP "could not allocate the nonlinear weak zones"
          PRINT 2000
          PRINT '(a)', "no. dgammadot0     x1       x2       x3  length   width thickn. strike   dip"
          PRINT 2000
          DO k=1,nnlwz
             CALL getdata(iunit,dataline)
             READ  (dataline,*) i, &
                  nonlinearweakzone(k)%dgammadot0, &
                  nonlinearweakzone(k)%x,nonlinearweakzone(k)%y,nonlinearweakzone(k)%z,&
                  nonlinearweakzone(k)%length,nonlinearweakzone(k)%width,nonlinearweakzone(k)%thickness, &
                  nonlinearweakzone(k)%strike,nonlinearweakzone(k)%dip
          
             nonlinearweakzonec(k)=nonlinearweakzone(k)
             
             PRINT '(I3.3,4ES9.2E1,3ES8.2E1,f7.1,f6.1)',k,&
                  nonlinearweakzone(k)%dgammadot0, &
                  nonlinearweakzone(k)%x,nonlinearweakzone(k)%y,nonlinearweakzone(k)%z, &
                  nonlinearweakzone(k)%length,nonlinearweakzone(k)%width, &
                  nonlinearweakzone(k)%thickness, &
                  nonlinearweakzone(k)%strike,nonlinearweakzone(k)%dip
             
             IF (i .ne. k) THEN
                WRITE_DEBUG_INFO
                WRITE (0,'("error in input file: source index misfit")')
                STOP 1
             END IF
             ! comply to Wang's convention
             CALL wangconvention( &
                  dummy, & 
                  nonlinearweakzone(k)%x,nonlinearweakzone(k)%y,nonlinearweakzone(k)%z, &
                  nonlinearweakzone(k)%length,nonlinearweakzone(k)%width, &
                  nonlinearweakzone(k)%strike,nonlinearweakzone(k)%dip,dummy,rot)
          END DO
       END IF
    END IF ! end nonlinear viscous

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !                 F A U L T    C R E E P
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    PRINT '(a)', "number of fault creep interfaces"
    CALL getdata(unit,dataline)
    READ  (dataline,*) nfc
    PRINT '(I5)', nfc

    IF (nfc .GT. 0) THEN
       ALLOCATE(faultcreeplayer(nfc),faultcreepstruc(sx3/2),STAT=iostatus)
       IF (iostatus>0) STOP "could not allocate the layer structure"

       PRINT 2000
       PRINT '(a)', "no.    depth   gamma0 (a-b)sig friction cohesion"
       PRINT 2000
       DO k=1,nfc
          CALL getdata(unit,dataline)
          READ  (dataline,*) i,faultcreeplayer(k)%z, &
               faultcreeplayer(k)%gammadot0, &
               faultcreeplayer(k)%stressexponent, &
               faultcreeplayer(k)%friction, &
               faultcreeplayer(k)%cohesion

          PRINT '(I3.3,5ES9.2E1)', i, &
               faultcreeplayer(k)%z, &
               faultcreeplayer(k)%gammadot0, &
               faultcreeplayer(k)%stressexponent, &
               faultcreeplayer(k)%friction, &
               faultcreeplayer(k)%cohesion

          ! check positive strain rates
          IF (faultcreeplayer(k)%gammadot0 .LT. 0) THEN
             WRITE_DEBUG_INFO
             WRITE (0,'("error in input file: slip rates must be positive")')
             STOP 1
          END IF

          IF (i .ne. k) THEN
             WRITE_DEBUG_INFO
             WRITE (0,'("error in input file: index misfit")')
             STOP 1
          END IF

       END DO

       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !             A F T E R S L I P       P L A N E S
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       PRINT '(a)', "number of afterslip planes"
       CALL getdata(unit,dataline)
       READ  (dataline,*) np
       PRINT '(I5)', np
       
       IF (np .gt. 0) THEN
          ALLOCATE(n(np),STAT=iostatus)
          IF (iostatus>0) STOP "could not allocate the plane list"
       
          PRINT 2000
          PRINT 2500
          PRINT 2000
          
          DO k=1,np
             CALL getdata(unit,dataline)
             READ  (dataline,*) i,n(k)%x,n(k)%y,n(k)%z,&
                  n(k)%length,n(k)%width,n(k)%strike,n(k)%dip,n(k)%rake
             
             PRINT '(I3.3," ",5ES9.2E1,3f7.1)',i, &
                  n(k)%x,n(k)%y,n(k)%z, &
                  n(k)%length,n(k)%width,n(k)%strike,n(k)%dip,n(k)%rake

             IF (i .ne. k) THEN
                WRITE_DEBUG_INFO
                WRITE (0,'("error in input file: plane index misfit")')
                STOP 1
             END IF

             ! modify rake for consistency with slip model
             IF (n(k)%rake .GE. 0.d0) THEN
                n(k)%rake=n(k)%rake-180.d0
             ELSE             
                n(k)%rake=n(k)%rake+180.d0
             END IF

             ! comply to Wang's convention
             CALL wangconvention(dummy,n(k)%x,n(k)%y,n(k)%z,&
                  n(k)%length,n(k)%width,n(k)%strike,n(k)%dip,n(k)%rake,rot)

#ifdef VTK
             ! export the afterslip segment in VTK format
             j=INDEX(wdir," ")
             WRITE (digit,'(I3.3)') k

             rffilename=wdir(1:j-1)//"/aplane-"//digit//".vtp"
             CALL exportvtk_rectangle(n(k)%x,n(k)%y,n(k)%z,n(k)%length,n(k)%width, &
                                      n(k)%strike,n(k)%dip,rffilename)
#endif

          END DO
       END IF
       
    END IF

    ! - - - - - - - - - - - - - - - - - - - - - - - - - -
    !     I N T E R - S E I S M I C    L O A D I N G
    ! - - - - - - - - - - - - - - - - - - - - - - - - - -
    minlength=sx1*dx1+sx2*dx2
    minwidth=sx3*dx3
    
    ! - - - - - - - - - - - - - - - - - - - - - - - - - -
    !        S H E A R     S O U R C E S   R A T E
    ! - - - - - - - - - - - - - - - - - - - - - - - - - -
    PRINT '(a)', "number of inter-seismic strike-slip segments"
    CALL getdata(iunit,dataline)
    READ  (dataline,*) inter%ns
    PRINT '(I5)', inter%ns
    IF (inter%ns .GT. 0) THEN
       ALLOCATE(inter%s(inter%ns),inter%sc(inter%ns),STAT=iostatus)
       IF (iostatus>0) STOP "could not allocate the source list"
       PRINT 2000
       PRINT '(a)',"no.  slip  xs ys zs  length width  strike dip rake"
       PRINT 2000
       DO k=1,inter%ns
          CALL getdata(iunit,dataline)
          READ (dataline,*) i,inter%s(k)%slip, &
               inter%s(k)%x,inter%s(k)%y,inter%s(k)%z, &
               inter%s(k)%length,inter%s(k)%width, &
               inter%s(k)%strike,inter%s(k)%dip,inter%s(k)%rake
          ! copy the input format for display
          inter%sc(k)=inter%s(k)
             
          PRINT '(I3.3,4ES9.2E1,2ES8.2E1,f7.1,f6.1,f7.1)',i, &
               inter%sc(k)%slip,&
               inter%sc(k)%x,inter%sc(k)%y,inter%sc(k)%z, &
               inter%sc(k)%length,inter%sc(k)%width, &
               inter%sc(k)%strike,inter%sc(k)%dip, &
               inter%sc(k)%rake
          
          IF (i .ne. k) THEN
             WRITE_DEBUG_INFO
             WRITE (0,'("error in input file: source index misfit")')
             STOP 1
          END IF
          IF (MAX(inter%s(k)%length,inter%s(k)%width) .LE. 0._8) THEN
             WRITE_DEBUG_INFO
             WRITE (0,'("error in input file: lengths must be positive.")')
             STOP 1
          END IF
          IF (inter%s(k)%length .lt. minlength) THEN
             minlength=inter%s(k)%length
          END IF
          IF (inter%s(k)%width  .lt. minwidth ) THEN
             minwidth =inter%s(k)%width
          END IF
          
          ! smooth out the slip distribution
          CALL antialiasingfilter(inter%s(k)%slip, &
                      inter%s(k)%length,inter%s(k)%width, &
                      dx1,dx2,dx3,nyquist)

          ! comply to Wang's convention
          CALL wangconvention(inter%s(k)%slip, &
               inter%s(k)%x,inter%s(k)%y,inter%s(k)%z, &
               inter%s(k)%length,inter%s(k)%width, &
               inter%s(k)%strike,inter%s(k)%dip, &
               inter%s(k)%rake,rot)

       END DO
       PRINT 2000
    END IF
    
    ! - - - - - - - - - - - - - - - - - - - - - - - - - -
    !       T E N S I L E   S O U R C E S   R A T E
    ! - - - - - - - - - - - - - - - - - - - - - - - - - -
    PRINT '(a)', "number of inter-seismic tensile segments"
    CALL getdata(iunit,dataline)
    READ  (dataline,*) inter%nt
    PRINT '(I5)', inter%nt
    IF (inter%nt .GT. 0) THEN
       ALLOCATE(inter%ts(inter%nt),inter%tsc(inter%nt),STAT=iostatus)
       IF (iostatus>0) STOP "could not allocate the tensile source list"
       PRINT 2000
       PRINT '(a)',"no. opening xs ys zs  length width  strike dip"
       PRINT 2000
       DO k=1,inter%nt
          CALL getdata(iunit,dataline)
          READ  (dataline,*) i,inter%ts(k)%slip, &
               inter%ts(k)%x,inter%ts(k)%y,inter%ts(k)%z, &
               inter%ts(k)%length,inter%ts(k)%width, &
               inter%ts(k)%strike,inter%ts(k)%dip
          ! copy the input format for display
          inter%tsc(k)=inter%ts(k)
          
          PRINT '(I3.3,4ES9.2E1,2ES8.2E1,f7.1,f6.1)', i, &
               inter%tsc(k)%slip,&
               inter%tsc(k)%x,inter%tsc(k)%y,inter%tsc(k)%z, &
               inter%tsc(k)%length,inter%tsc(k)%width, &
               inter%tsc(k)%strike,inter%tsc(k)%dip
          
          IF (i .ne. k) THEN
             WRITE_DEBUG_INFO
             WRITE (0,'("error in input file: tensile source index misfit")')
             STOP 1
          END IF
          IF (MAX(inter%ts(k)%length,inter%ts(k)%width) .LE. 0._8) THEN
             WRITE_DEBUG_INFO
             WRITE (0,'("error in input file: lengths must be positive.")')
             STOP 1
          END IF
          IF (inter%ts(k)%length .lt. minlength) THEN
             minlength=inter%ts(k)%length
          END IF
          IF (inter%ts(k)%width  .lt. minwidth) THEN
             minwidth =inter%ts(k)%width
          END IF
          
          ! smooth out the slip distribution
          CALL antialiasingfilter(inter%ts(k)%slip, &
                           inter%ts(k)%length,inter%ts(k)%width, &
                           dx1,dx2,dx3,nyquist)

          ! comply to Wang's convention
          CALL wangconvention(inter%ts(k)%slip, &
               inter%ts(k)%x,inter%ts(k)%y,inter%ts(k)%z, &
               inter%ts(k)%length,inter%ts(k)%width, &
               inter%ts(k)%strike,inter%ts(k)%dip,dummy,rot)

       END DO
       PRINT 2000
    END IF
       
    ! - - - - - - - - - - - - - - - - - - - - - - - - - -
    !       C 0 - S E I S M I C     E V E N T S
    ! - - - - - - - - - - - - - - - - - - - - - - - - - -
    PRINT '(a)', "number of events"
    CALL getdata(iunit,dataline)
    READ (dataline,*) ne
    PRINT '(I5)', ne
    IF (ne .GT. 0) ALLOCATE(events(ne),STAT=iostatus)
    IF (iostatus>0) STOP "could not allocate the event list"
    
    DO e=1,ne
       IF (1 .NE. e) THEN
          PRINT '("time of next coseismic event")'
          CALL getdata(iunit,dataline)
          READ (dataline,*) events(e)%time
          
          IF (0 .EQ. skip) THEN
             ! change event time to multiples of output time step
             events(e)%time=fix(events(e)%time/odt)*odt
          END IF

          PRINT '(ES9.2E1," (multiple of ",ES9.2E1,")")', &
               events(e)%time,odt

          IF (events(e)%time .LE. events(e-1)%time) THEN
             WRITE_DEBUG_INFO
             WRITE (0,'(a,a)') "input file error. ", &
                  "coseismic source time must increase. interrupting."
             STOP 1
          END IF
       ELSE
          events(1)%time=0._8
          events(1)%i=0
       END IF

       ! - - - - - - - - - - - - - - - - - - - - - - - - - -
       !           S H E A R     S O U R C E S
       ! - - - - - - - - - - - - - - - - - - - - - - - - - -
       PRINT '(a)', "number of coseismic strike-slip segments (ns)"
       CALL getdata(iunit,dataline)
       READ  (dataline,*) events(e)%ns
       PRINT '(I5)', events(e)%ns
       IF (events(e)%ns .GT. 0) THEN
          ALLOCATE(events(e)%s(events(e)%ns),events(e)%sc(events(e)%ns), &
               STAT=iostatus)
          IF (iostatus>0) STOP "could not allocate the source list"
          PRINT 2000
          PRINT '(a)',"no.     slip       xs       ys       zs  length   width strike   dip   rake"
          PRINT 2000
          DO k=1,events(e)%ns
             CALL getdata(iunit,dataline)
             READ (dataline,*) i,events(e)%s(k)%slip, &
                  events(e)%s(k)%x,events(e)%s(k)%y,events(e)%s(k)%z, &
                  events(e)%s(k)%length,events(e)%s(k)%width, &
                  events(e)%s(k)%strike,events(e)%s(k)%dip,events(e)%s(k)%rake
             ! copy the input format for display
             events(e)%sc(k)=events(e)%s(k)
             
             PRINT '(I3.3,4ES9.2E1,2ES8.2E1,f7.1,f6.1,f7.1)',i, &
                  events(e)%sc(k)%slip,&
                  events(e)%sc(k)%x,events(e)%sc(k)%y,events(e)%sc(k)%z, &
                  events(e)%sc(k)%length,events(e)%sc(k)%width, &
                  events(e)%sc(k)%strike,events(e)%sc(k)%dip, &
                  events(e)%sc(k)%rake
             
             IF (i .ne. k) THEN
                WRITE_DEBUG_INFO
                WRITE (0,'("invalid shear source definition ")')
                WRITE (0,'("error in input file: source index misfit")')
                STOP 1
             END IF
             IF (MAX(events(e)%s(k)%length,events(e)%s(k)%width) .LE. 0._8) THEN
                WRITE_DEBUG_INFO
                WRITE (0,'("error in input file: lengths must be positive.")')
                STOP 1
             END IF
             IF (events(e)%s(k)%length .lt. minlength) THEN
                minlength=events(e)%s(k)%length
             END IF
             IF (events(e)%s(k)%width  .lt. minwidth ) THEN
                minwidth =events(e)%s(k)%width
             END IF
             
             ! smooth out the slip distribution
             CALL antialiasingfilter(events(e)%s(k)%slip, &
                              events(e)%s(k)%length,events(e)%s(k)%width, &
                              dx1,dx2,dx3,nyquist)

             ! comply to Wang's convention
             CALL wangconvention(events(e)%s(k)%slip, &
                  events(e)%s(k)%x,events(e)%s(k)%y,events(e)%s(k)%z, &
                  events(e)%s(k)%length,events(e)%s(k)%width, &
                  events(e)%s(k)%strike,events(e)%s(k)%dip, &
                  events(e)%s(k)%rake,rot)

          END DO

#ifdef VTK
          ! export the fault segments in VTK format for the current event
          j=INDEX(wdir," ")
          WRITE (digit,'(I3.3)') e

          rffilename=wdir(1:j-1)//"/rfaults-"//digit//".vtp"
          CALL exportvtk_rfaults(events(e),rffilename)
#endif
          rffilename=wdir(1:j-1)//"/rfaults-"//digit//".xy"
          CALL exportxy_rfaults(events(e),rffilename)

          PRINT 2000
       END IF
       
       ! - - - - - - - - - - - - - - - - - - - - - - - - - -
       !          T E N S I L E      S O U R C E S
       ! - - - - - - - - - - - - - - - - - - - - - - - - - -
       PRINT '(a)', "number of coseismic tensile segments (nt)"
       CALL getdata(iunit,dataline)
       READ  (dataline,*) events(e)%nt
       PRINT '(I5)', events(e)%nt
       IF (events(e)%nt .GT. 0) THEN
          ALLOCATE(events(e)%ts(events(e)%nt),events(e)%tsc(events(e)%nt), &
               STAT=iostatus)
          IF (iostatus>0) STOP "could not allocate the tensile source list"
          PRINT 2000
          PRINT '(a)',"no. opening xs ys zs  length width  strike dip"
          PRINT 2000
          DO k=1,events(e)%nt
             CALL getdata(iunit,dataline)
             READ  (dataline,*) i,events(e)%ts(k)%slip, &
                  events(e)%ts(k)%x,events(e)%ts(k)%y,events(e)%ts(k)%z, &
                  events(e)%ts(k)%length,events(e)%ts(k)%width, &
                  events(e)%ts(k)%strike,events(e)%ts(k)%dip
             ! copy the input format for display
             events(e)%tsc(k)=events(e)%ts(k)
             
             PRINT '(I3.3,4ES9.2E1,2ES8.2E1,f7.1,f6.1)',k, &
                  events(e)%tsc(k)%slip,&
                  events(e)%tsc(k)%x,events(e)%tsc(k)%y,events(e)%tsc(k)%z, &
                  events(e)%tsc(k)%length,events(e)%tsc(k)%width, &
                  events(e)%tsc(k)%strike,events(e)%tsc(k)%dip
             
             IF (i .ne. k) THEN
                PRINT *, "error in input file: source index misfit"
                STOP 1
             END IF
             IF (events(e)%ts(k)%length .lt. minlength) THEN
                minlength=events(e)%ts(k)%length
             END IF
             IF (events(e)%ts(k)%width  .lt. minwidth) THEN
                minwidth =events(e)%ts(k)%width
             END IF
             
             ! smooth out the slip distribution
             CALL antialiasingfilter(events(e)%ts(k)%slip, &
                              events(e)%ts(k)%length,events(e)%ts(k)%width, &
                              dx1,dx2,dx3,nyquist)

             ! comply to Wang's convention
             CALL wangconvention(events(e)%ts(k)%slip, &
                  events(e)%ts(k)%x,events(e)%ts(k)%y,events(e)%ts(k)%z, &
                  events(e)%ts(k)%length,events(e)%ts(k)%width, &
                  events(e)%ts(k)%strike,events(e)%ts(k)%dip,dummy,rot)

          END DO
          PRINT 2000
       END IF
       
       ! - - - - - - - - - - - - - - - - - - - - - - - - - -
       !                M O G I      S O U R C E S
       ! - - - - - - - - - - - - - - - - - - - - - - - - - -
       PRINT '(a)', "number of coseismic dilatation point sources"
       CALL getdata(iunit,dataline)
       READ  (dataline,*) events(e)%nm
       PRINT '(I5)', events(e)%nm
       IF (events(e)%nm .GT. 0) THEN
          ALLOCATE(events(e)%m(events(e)%nm),events(e)%mc(events(e)%nm), &
               STAT=iostatus)
          IF (iostatus>0) STOP "could not allocate the tensile source list"
          PRINT 2000
          PRINT '(a)',"no. strain (positive for extension) xs ys zs"
          PRINT 2000
          DO k=1,events(e)%nm
             CALL getdata(iunit,dataline)
             READ  (dataline,*) i,events(e)%m(k)%slip, &
                  events(e)%m(k)%x,events(e)%m(k)%y,events(e)%m(k)%z
             ! copy the input format for display
             events(e)%mc(k)=events(e)%m(k)
             
             PRINT '(I3.3,4ES9.2E1)',k, &
                  events(e)%mc(k)%slip,&
                  events(e)%mc(k)%x,events(e)%mc(k)%y,events(e)%mc(k)%z
             
             IF (i .ne. k) THEN
                PRINT *, "error in input file: source index misfit"
                STOP 1
             END IF
             
             ! rotate the source in the computational reference frame
             CALL rotation(events(e)%m(k)%x,events(e)%m(k)%y,rot)
          END DO
          PRINT 2000
       END IF

       ! - - - - - - - - - - - - - - - - - - - - - - - - - -
       !             S U R F A C E   L O A D S
       ! - - - - - - - - - - - - - - - - - - - - - - - - - -
       PRINT '(a)', "number of surface loads"
       CALL getdata(iunit,dataline)
       READ  (dataline,*) events(e)%nl
       PRINT '(I5)', events(e)%nl
       IF (events(e)%nl .GT. 0) THEN
          ALLOCATE(events(e)%l(events(e)%nl),events(e)%lc(events(e)%nl), &
               STAT=iostatus)
          IF (iostatus>0) STOP "could not allocate the load list"
          PRINT 2000
          PRINT '(a)',"no. xs ys t3 (force/surface/rigidity, positive down)"
          PRINT 2000
          DO k=1,events(e)%nl
             CALL getdata(iunit,dataline)
             READ  (dataline,*) i, &
                  events(e)%l(k)%x,events(e)%l(k)%y,events(e)%l(k)%slip
             ! copy the input format for display
             events(e)%lc(k)=events(e)%l(k)
             
             PRINT '(I3.3,4ES9.2E1)',k, &
                  events(e)%lc(k)%x,events(e)%lc(k)%y,events(e)%lc(k)%slip
             
             IF (i .NE. k) THEN
                PRINT *, "error in input file: source index misfit"
                STOP 1
             END IF
             
             ! rotate the source in the computational reference frame
             CALL rotation(events(e)%l(k)%x,events(e)%l(k)%y,rot)
          END DO
          PRINT 2000
       END IF
       
    END DO

    ! test the presence of dislocations for coseismic calculation
    IF ((events(1)%nt .EQ. 0) .AND. &
        (events(1)%ns .EQ. 0) .AND. &
        (events(1)%nm .EQ. 0) .AND. &
        (events(1)%nl .EQ. 0) .AND. &
        (interval .LE. 0._8)) THEN

       WRITE_DEBUG_INFO
       WRITE (0,'("**** error **** ")')
       WRITE (0,'("no input dislocations or dilatation point sources")')
       WRITE (0,'("or surface tractions for first event . exiting.")')
       STOP 1
    END IF

    ! maximum recommended sampling size
    PRINT '(a,2ES8.2E1)', &
         "max sampling size (hor.,vert.):", minlength/2.5_8,minwidth/2.5_8

    PRINT 2000

2000 FORMAT ("----------------------------------------------------------------------------")
2100 FORMAT ("no.        x1       x2       x3   length    width strike    dip")
2200 FORMAT ("no. slip        x1         x2         x3    length   width strike  dip  rake")
2300 FORMAT ("no. name       x1       x2       x3 (name is a 4-character string)")
2400 FORMAT ("no. strain       x1       x2       x3 (positive for extension)")
2500 FORMAT ("no.        x1       x2       x3   length    width strike    dip   rake")

  END SUBROUTINE init

  !--------------------------------------------------------------------
  ! function IsOutput
  ! checks if output should be written based on user choices: if output
  ! time interval (odt) is positive, output is written only if time
  ! is an integer of odt. If odt is negative output is written at times
  ! corresponding to internally chosen time steps.
  !
  ! IsOutput is true only at discrete intervals (skip=0,odt>0),
  ! or at every "skip" computational steps (skip>0,odt<0),
  ! or anytime a coseismic event occurs
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
  ! subroutine IntegrationStep
  ! find the time-integration forward step for the predictor-corrector
  ! scheme.
  !
  ! input file line
  !
  !    time interval, (positive dt step) or (negative skip and scaling)
  !
  ! can be filled by either 1)
  !
  !   T, dt
  !
  ! where T is the time interval of the simulation and dt is the
  ! output time step, or 2)
  !
  !   T, -n, t_s
  !
  ! where n indicates the number of computational steps before 
  ! outputing results, t_s is a scaling applied to internally
  ! computed time step.
  !
  ! for case 1), an optimal time step is evaluated internally to
  ! ensure stability (t_m/10) of time integration. The actual
  ! time step Dt is chosen as
  !
  !    Dt = min( t_m/10, ((t%odt)+1)*odt-t )
  !
  ! where t is the current time in the simulation. regardless of 
  ! time step Dt, results are output if t is a multiple of dt.
  !
  ! for case 2), the time step is chosen internally based on an 
  ! estimate of the relaxation time (t_m/10). Results are output
  ! every n steps. The actual time step is chosen as
  !
  !    Dt = min( t_m/10*t_s, t(next event)-t )
  !
  ! where index is the number of computational steps after a coseismic
  ! event and t(next event) is the time of the next coseismic event.
  !
  ! sylvain barbot (01/01/08) - original form 
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

    alpha=rot*DEG2RAD
    xx=x
    yy=y

    x=+xx*cos(alpha)+yy*sin(alpha)
    y=-xx*sin(alpha)+yy*cos(alpha)

  END SUBROUTINE rotation

  !-------------------------------------------------------------------
  ! subroutine AntiAliasingFilter
  ! smoothes a slip distribution model to avoid aliasing of
  ! the source geometry. Aliasing occurs is a slip patch has 
  ! dimensions (width or length) smaller than the grid sampling.
  !
  ! if a patch length is smaller than a critical size L=dx*nyquist, it 
  ! is increased to L and the slip (or opening) is scaled accordingly
  ! so that the moment M = s*L*W is conserved.
  !
  ! sylvain barbot (12/08/09) - original form
  !-------------------------------------------------------------------
  SUBROUTINE antialiasingfilter(slip,length,width,dx1,dx2,dx3,nyquist)
    REAL*8, INTENT(INOUT) :: slip,length,width
    REAL*8, INTENT(IN) :: dx1,dx2,dx3,nyquist

    REAL*8 :: dx

    ! minimum slip patch dimension
    dx=MIN(dx1,dx2,dx3)*nyquist

    ! update length
    IF (length .LT. dx) THEN
       slip=slip*length/dx
       length=dx
    END IF
    ! update width
    IF (width .LT. dx) THEN
       slip=slip*width/dx
       width=dx
    END IF

  END SUBROUTINE antialiasingfilter

  !------------------------------------------------------------------
  ! subroutine WangConvention
  ! converts a fault slip model from a geologic description including
  ! fault length, width, strike, dip and rake into a description
  ! compatible with internal convention of the program.
  !
  ! Internal convention describes a fault patch by the location of
  ! its center, instead of an upper corner and its orientation by
  ! the deviation from the vertical, instead of the angle from the
  ! horizontal and by the angle from the x2 axis (East-West)
  !------------------------------------------------------------------
  SUBROUTINE wangconvention(slip,x,y,z,length,width,strike,dip,rake,rot)
    REAL*8, INTENT(OUT) :: slip, x,y,z,strike,dip,rake
    REAL*8, INTENT(IN) :: length,width,rot

    slip=-slip
    strike=-90._8-strike
    dip   = 90._8-dip

    strike=strike*DEG2RAD
    dip=dip*DEG2RAD
    rake=rake*DEG2RAD

    x=x-x0-length/2._8*sin(strike)+width /2._8*sin(dip)*cos(strike)
    y=y-y0-length/2._8*cos(strike)-width /2._8*sin(dip)*sin(strike)
    z=z+width /2._8*cos(dip)

    CALL rotation(x,y,rot)

    strike=strike+rot*DEG2RAD

  END SUBROUTINE wangconvention
  
END PROGRAM relax
