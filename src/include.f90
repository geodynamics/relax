#include "config.h"

! implement SGI Fast Fourier Transforms library
!#define SGI_FFT 1

! export data to GMT XYZ text format
!#define XYZ 1

! export data to GMT GRD binary format
#define GRD 1

! export equivalent body forces in GRD format
!#define GRD_EQBF 1

! export amplitude of scalar fields 
! along a plane in GRD binary format
#define GRD_EXPORTEIGENSTRAIN 1

! export creep velocity along a frictional 
! plane in GRD binary format
#define GRD_EXPORTCREEP 1

! export data to the TXT format
!#define TXT 1

! export amplitude of scalar fields along 
! an observation plane in text format
#define TXT_EXPORTEIGENSTRAIN 1

! export creep velocity along a frictional 
! plane in text format
!#define TXT_EXPORTCREEP 1

! export data to VTK format (for visualization in Paraview)
#define VTK 1
!#define VTK_EQBF 1

#define WRITE_DEBUG_INFO WRITE (0,'("error at line ",I5.5," of source file ",a)') __LINE__,__FILE__


#ifdef IMKL_FFT
#define WRITE_MKL_DEBUG_INFO(i) IF(i.NE.0)THEN;IF(.NOT.DftiErrorClass(i,DFTI_NO_ERROR))THEN;WRITE_DEBUG_INFO;WRITE (0,*) DftiErrorMessage(i);STOP 1;END IF;END IF
#endif

! adjust data alignment for the Fourier transform
#ifdef FFTW3
#define ALIGN_DATA 1
#else
#ifdef SGI_FFT
#define ALIGN_DATA 1
#else
#ifdef IMKL_FFT
#define ALIGN_DATA 1
#endif
#endif
#endif
