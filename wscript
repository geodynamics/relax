import os
def options(opt):
    opt.load('compiler_c compiler_fc')
    opt.load('compiler_cxx')

    mkl=opt.add_option_group('MKL Options')
    mkl.add_option('--mkl-dir',
                   help='Base directory where mkl is installed')
    mkl.add_option('--mkl-incdir',
                   help='Directory where mkl include files are installed')
    mkl.add_option('--mkl-libdir',
                   help='Directory where mkl library files are installed')
    mkl.add_option('--mkl-libs',
                   help='Names of the mkl libraries without prefix or suffix\n'
                   '(e.g. "mkl_intel_lp64 mkl_intel_thread mkl_core"')

    fftw=opt.add_option_group('FFTW Options')
    fftw.add_option('--use-fftw', action='store_true', default=False,
                    help='Use FFTW instead of MKL')
    fftw.add_option('--fftw-dir',
                    help='Base directory where fftw is installed')
    fftw.add_option('--fftw-incdir',
                    help='Directory where fftw include files are installed')
    fftw.add_option('--fftw-libdir',
                    help='Directory where fftw library files are installed')

    cuda=opt.add_option_group('CUDA Options')
    cuda.add_option('--use-cuda', action='store_true', default=False,
                    help='Uses GPU for computation')
    cuda.add_option('--cuda-dir',
                    help='Base Directory where the cuda is installed')
    cuda.add_option('--cuda-incdir',
                    help='Directory where cuda include files are installed.')
    cuda.add_option('--cuda-libdir',
                    help='Directory where cuda library files are installed')

    papi=opt.add_option_group('PAPI Options')
    papi.add_option('--use-papi', action='store_true', default=False,
                    help='Use PAPI for profiling')
    papi.add_option('--papi-dir',
                    help='Base directory where papi is installed')
    papi.add_option('--papi-incdir',
                    help='Directory where papi include files are installed')
    papi.add_option('--papi-libdir',
                    help='Directory where papi library files are installed')

    proj=opt.add_option_group('Proj Options')
    proj.add_option('--proj-dir',
                    help='Base directory where proj is installed')
    proj.add_option('--proj-incdir',
                    help='Directory where proj include files are installed')
    proj.add_option('--proj-libdir',
                    help='Directory where proj library files are installed')

    gmt=opt.add_option_group('Gmt Options')
    gmt.add_option('--gmt-dir',
                    help='Base directory where gmt is installed')
    gmt.add_option('--gmt-incdir',
                    help='Directory where gmt include files are installed')
    gmt.add_option('--gmt-libdir',
                    help='Directory where gmt library files are installed')

    other=opt.add_option_group('Other Options')
    other.add_option('--openmp-flag',
                     help="C and Fortran compiler flag for OpenMP")
    other.add_option('--cpp-flag',
                     help="Fortran compiler flag to run the C preprocessor")
    other.add_option('--use-ctfft', action='store_true', default=False,
                     help='Use slow internal CTFFT instead of MKL or FFTW')
    other.add_option('--length-flag',
                     help='Fortran compiler option to allow unlimited line length')
    other.add_option('--relax-lite', action='store_true', default=False,
                     help="Generating the relax library")

def configure(cnf):
    cnf.load('compiler_c compiler_fc')
    cnf.load('compiler_cxx')

 # We set the flags here 
    if not cnf.env.CFLAGS:
        cnf.env.CFLAGS=['-O3','-fPIC']
    if not cnf.env.FCFLAGS:
        cnf.env.FCFLAGS=['-O3']

    cnf.check_fortran()

    # Find Proj
    if cnf.options.proj_dir:
        if not cnf.options.proj_incdir:
            cnf.options.proj_incdir=cnf.options.proj_dir + "/include"
        if not cnf.options.proj_libdir:
            cnf.options.proj_libdir=cnf.options.proj_dir + "/lib"
    cnf.check_cc(header_name='proj_api.h',uselib_store='proj',
                 includes=[cnf.options.proj_incdir],
                 libpath=[cnf.options.proj_libdir],
                 rpath=[cnf.options.proj_libdir],
                 lib='proj',define_name="PROJ")

    # Find GMT
    if not cnf.options.relax_lite:
        if cnf.options.gmt_dir:
            if not cnf.options.gmt_incdir:
                cnf.options.gmt_incdir=cnf.options.gmt_dir + "/include"
            if not cnf.options.gmt_libdir:
                cnf.options.gmt_libdir=cnf.options.gmt_dir + "/lib"
        if cnf.options.gmt_incdir:
            includedirs=[cnf.options.gmt_incdir]
        else:
            includedirs=['','/usr/include/gmt']
        found_gmt=False
        for inc in includedirs:
            try:
                cnf.check_cc(msg="Checking for gmt.h in '" + inc + "'",
                             header_name='gmt.h', includes=inc,
                             libpath=[cnf.options.gmt_libdir],
                             rpath=[cnf.options.gmt_libdir],
                             lib=['gmt','netcdf'], uselib_store='gmt')
            except cnf.errors.ConfigurationError:
                pass
            else:
                found_gmt=True
                break
            if not found_gmt:
                cnf.fatal('Could not find gmt')

    #Find Cuda
    if cnf.options.use_cuda:
        cnf.env.CUDA=cnf.options.use_cuda
        cnf.load('cuda',tooldir='.')   
        if not cnf.env.CUDAFLAGS:
            cnf.env.CUDAFLAGS = ['-gencode','arch=compute_35,code=sm_35','-Xcompiler','-fPIC']
    #       cnf.env.CUDAFLAGS += ['-Xptxas', '-dlcm=cg']
    #       cnf.env.CUDAFLAGS += ['--maxrregcount=32']
    #       cnf.env.CUDAFLAGS = ['-gencode','arch=compute_30,code=sm_30']
    #       cnf.env.CUDAFLAGS = ['-gencode','arch=compute_20,code=sm_20']
            cnf.env.CXXFLAGS=['-m64']
        if cnf.options.cuda_dir:
            if not cnf.options.cuda_incdir:
                cnf.options.cuda_incdir=cnf.options.cuda_dir + "/include"
            if not cnf.options.cuda_libdir:
                cnf.options.cuda_libdir=cnf.options.cuda_dir + "/lib64"
        if cnf.options.cuda_incdir:
            includedirs=[cnf.options.cuda_incdir]
        else:
            includedirs=['','/usr/local/cuda/include']
        found_cuda=False
        for inc in includedirs:
            try:
                cnf.check_cc(msg="Checking for cuda.h", 
                             header_name='cuda.h', includes=inc,
                             libpath=[cnf.options.cuda_libdir], 
                             rpath=[cnf.options.cuda_libdir], 
                             lib=['cudart', 'cufft','stdc++'], uselib_store='cuda',
                             define_name="USING_CUDA") 
            except cnf.errors.ConfigurationError:
                pass
            else:
                found_cuda=True
                break
        if not found_cuda:
            cnf.fatal('Could not find cuda')

    # Find PAPI
    if cnf.options.use_papi:
        if cnf.options.papi_dir:
            if not cnf.options.papi_incdir:
                cnf.options.papi_incdir=cnf.options.papi_dir + "/include"
            if not cnf.options.papi_libdir:
                cnf.options.papi_libdir=cnf.options.papi_dir + "/lib"
        if cnf.options.papi_incdir:
            includedirs=[cnf.options.papi_incdir]
        else:
            includedirs=['','/usr/include/papi']
        found_papi=False
        for inc in includedirs:
            try:
                cnf.check_cc(msg="Checking for papi.h in '" + inc + "'",
                             header_name='papi.h', includes=inc,
                             libpath=[cnf.options.papi_libdir],
                             rpath=[cnf.options.papi_libdir],
                             lib=['papi'], uselib_store='papi',define_name="PAPI_PROF")
            except cnf.errors.ConfigurationError:
                pass
            else:
                found_papi=True
                break
        if not found_papi:
            cnf.fatal('Could not find papi')

    # Find OpenMP
    openmp_msg="Checking for openmp flag "
    openmp_fragment="program main\n  call omp_get_num_threads()\nend program main"
    found_openmp=False
    openmp_flags=['-fopenmp','-openmp','-mp','-xopenmp','-omp','-qsmp=omp']
    if cnf.options.openmp_flag:
        openmp_flags=[cnf.options.openmp_flag]
    for flag in openmp_flags:
        try:
            cnf.check_fc(msg=openmp_msg+flag, fragment=openmp_fragment,
                         fcflags=flag, linkflags=flag, uselib_store='openmp')
        except cnf.errors.ConfigurationError:
            continue
        else:
            found_openmp=True
            break
    if not found_openmp:
        cnf.fatal('Could not find OpenMP')

    # Find FFTW
    #if not cnf.options.use_cuda:
    if cnf.options.use_fftw:
        if cnf.options.fftw_dir:
            if not cnf.options.fftw_incdir:
                cnf.options.fftw_incdir=cnf.options.fftw_dir + "/include"
            if not cnf.options.fftw_libdir:
                cnf.options.fftw_libdir=cnf.options.fftw_dir + "/lib"
        frag="program main\n" + 'include "fftw3.f"\n' \
            + "end program main\n"
        if cnf.options.fftw_incdir:
            inc=cnf.options.fftw_incdir
        else:
            inc='/usr/include'
        cnf.check_fc(msg="Checking for fftw",
                     includes=[inc], fragment=frag, uselib_store='fftw',
                     libpath=[cnf.options.fftw_libdir],
                     rpath=[cnf.options.fftw_libdir],
                     lib=['fftw3f','fftw3f_threads'], define_name="FFTW3")
    # Find MKL
    elif not cnf.options.use_ctfft:
        if cnf.options.mkl_dir:
            if not cnf.options.mkl_incdir:
                cnf.options.mkl_incdir=cnf.options.mkl_dir + "/include"
            if not cnf.options.mkl_libdir:
                cnf.options.mkl_libdir=cnf.options.mkl_dir + "/lib"
        if cnf.options.mkl_libs:
            libs=cnf.options.mkl_libs.split()
        else:
            libs=['mkl_intel_lp64', 'mkl_intel_thread',
                  'mkl_core']
        cnf.check_fc(msg="Checking for MKL", lib=libs,
                     uselib_store='imkl',
                     includes=[cnf.options.mkl_incdir],
                     libpath=[cnf.options.mkl_libdir],
                     rpath=[cnf.options.mkl_libdir],
                     use='openmp', define_name='IMKL_FFT')

    # Check for C preprocessor option
    frag="program main\n  INTEGER :: foo\n" + "end program main\n"
    cpp_flags=['-cpp','-Mpreprocess','-fpic']
    if cnf.options.cpp_flag:
        cpp_flags=[cnf.options.cpp_flag]
    found_cpp=False
    for flag in cpp_flags:
        try:
            cnf.check_fc(fragment=frag,msg="Checking preprocessor option " + flag,
                         fcflags=flag,uselib_store='cpp')
        except:
            continue
        else:
            found_cpp=True
            break
    if not found_cpp:
        cnf.fatal("Could not find an option for running the C preprocessor")

    # Check for line length option
    frag="program main\n" +"  WRITE(*,*) \"                                                                                                                                                      \"\n" + "end program main\n"
    length_flags=['','-ffree-line-length-none']
    if cnf.options.length_flag:
        length_flags=[cnf.options.length_flag]
    found_length=False
    for flag in length_flags:
        try:
            cnf.check_fc(fragment=frag,msg="Checking length option " + flag,
                         fcflags=flag,uselib_store='length')
        except:
            continue
        else:
            found_length=True
            break
    if not found_length:
        cnf.fatal("Could not find an option for allowing long lines")


    cnf.write_config_header('config.h')


def lite(ctx) :
    if(ctx.env.CUDA):
        ctx.shlib(features='c fc fcprogram cxx',
                source=['src/curelaxlite.f90',
                        'src/ctfft.f',
                        'src/types.f90',
                        'src/fourier.f90',
                        'src/green.f90',
                        'src/okada/green_space.f90',
                        'src/okada/dc3d.f',
                        'src/elastic3d.f90',
                        'src/friction3d.f90',
                        'src/viscoelastic3d.f90',
                        'src/proj.c',
                        'src/getdata.f',
                        'src/getopt_m.f90',
                        'src/util.f90',
                        'src/mkl_dfti.f90',
                        'src/papi_prof.c',
                        'src/cugreen.cu',
                        'src/cuelastic.cu'],
                install_path='${PREFIX}/bin',
                includes=['build'],
                use=['proj','openmp','fftw','imkl','cpp','length','cuda','papi','stdc++'],
                target='librelax.so'
                )
    else:
        ctx.shlib(features='c fc fcprogram',
                source=['src/relaxlite.f90',
                        'src/ctfft.f',
                        'src/types.f90',
                        'src/fourier.f90',
                        'src/green.f90',
                        'src/okada/green_space.f90',
                        'src/okada/dc3d.f',
                        'src/elastic3d.f90',
                        'src/friction3d.f90',
                        'src/viscoelastic3d.f90',
                        'src/proj.c',
                        'src/getdata.f',
                        'src/getopt_m.f90',
                        'src/util.f90',
                        'src/mkl_dfti.f90',
                        'src/papi_prof.c'],
                install_path='${PREFIX}/bin',
                includes=['build'],
                use=['proj','openmp','fftw','imkl','cpp','length','papi','stdc++'],
                target='librelax.so'
                )

from waflib.Build import BuildContext

class miracle(BuildContext):
    cmd = 'lite'
    fun = 'lite'

def build(bld):
    if bld.env.CUDA:    
        bld.program(features='c fc fcprogram cxx',
                source=['src/curelax.f90',
                        'src/ctfft.f',
                        'src/types.f90',
                        'src/fourier.f90',
                        'src/green.f90',
                        'src/okada/green_space.f90',
                        'src/okada/dc3d.f',
                        'src/elastic3d.f90',
                        'src/friction3d.f90',
                        'src/viscoelastic3d.f90',
                        'src/writevtk.c',
                        'src/writegrd4.2.c',
                        'src/proj.c',
                        'src/export.f90',
                        'src/getdata.f',
                        'src/getopt_m.f90',
                        'src/input.f90',
                        'src/util.f90',
                        'src/mkl_dfti.f90',
                        'src/papi_prof.c',
                        'src/cugreen.cu',
                        'src/cuelastic.cu'],
                install_path='${PREFIX}/bin',
                includes=['build'],
                use=['gmt','proj','openmp','fftw','imkl','cpp','length','cuda','papi','stdc++'],
                target='relax'
                )
    else:
        bld.program(features='c fc fcprogram',
                source=['src/relax.f90',
                        'src/ctfft.f',
						'src/types.f90',
                        'src/fourier.f90',
                        'src/green.f90',
                        'src/okada/green_space.f90',
                        'src/okada/dc3d.f',
                        'src/elastic3d.f90',
                        'src/friction3d.f90',
                        'src/viscoelastic3d.f90',
                        'src/writevtk.c',
                        'src/writegrd4.2.c',
                        'src/proj.c',
                        'src/export.f90',
                        'src/getdata.f',
                        'src/getopt_m.f90',
                        'src/input.f90',
                        'src/util.f90',
                        'src/mkl_dfti.f90',
                        'src/papi_prof.c'],
                install_path='${PREFIX}/bin',
                includes=['build'],
                use=['gmt','proj','openmp','fftw','imkl','cpp','length','papi','stdc++'],
                target='relax'
                )
