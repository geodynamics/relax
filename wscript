def options(opt):
    opt.load('compiler_c compiler_fc')

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
    other.add_option('--openmp-flag', help="Compiler flag for OpenMP")
    other.add_option('--use-ctfft', action='store_true', default=False,
                     help='Use slow internal CTFFT instead of MKL or FFTW')


def configure(cnf):
    cnf.load('compiler_c compiler_fc')

    # Find Proj
    if cnf.options.proj_dir:
        if not cnf.options.proj_incdir:
            cnf.options.proj_incdir=cnf.options.proj_dir + "/include"
        if not cnf.options.proj_libdir:
            cnf.options.proj_libdir=cnf.options.proj_dir + "/lib"
    cnf.check_cc(header_name='proj_api.h',uselib_store='proj',
                 includes=[cnf.options.proj_incdir],
                 libpath=[cnf.options.proj_libdir],
                 lib='proj')

    # Find GMT
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
                         lib=['gmt','netcdf'], uselib_store='gmt')
        except cnf.errors.ConfigurationError:
            pass
        else:
            found_gmt=True
            break
    if not found_gmt:
        cnf.fatal('Could not find gmt')

    # Find OpenMP
    openmp_msg="Checking for openmp flag "
    openmp_fragment="program main\n  call omp_get_num_threads()\nend program main"
    found_openmp=False
    openmp_flags=['-fopenmp','-openmp','-mp','-xopenmp','-omp','-qsmp=omp']
    if cnf.options.openmp_flag:
        openmp_flag=cnf.options.openmp_flag
    for flag in openmp_flags:
        try:
            cnf.check_fc(msg=openmp_msg+flag, fragment=openmp_fragment, fcflags=flag,
                         linkflags=flag, uselib_store='openmp')
        except cnf.errors.ConfigurationError:
            continue
        else:
            found_openmp=True
            break
    if not found_openmp:
        cnf.fatal('Could not find OpenMP')

    # Find FFTW
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
                     use='openmp', define_name='IMKL_FFT')

    cnf.write_config_header('config.h')

def build(bld):
    uses=['gmt','proj','openmp','fftw','imkl']

    bld.program(features='c fc fcprogram',
                source=['relax.f90',
                        'types.f90',
                        'ctfft.f',
                        'fourier.f90',
                        'green.f90',
                        'elastic3d.f90',
                        'friction3d.f90',
                        'viscoelastic3d.f90',
                        'writegrd4.2.c',
                        'proj.c',
                        'export.f90',
                        'getdata.f',
                        'getopt_m.f90',
                        'input.f90',
                        'mkl_dfti.f90'],
                includes=['build'],
                use=uses,
                fcflags=['-cpp','-zero'],
                target='relax'
                )
