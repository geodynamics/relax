def options(opt):
    opt.load('compiler_c compiler_fc')
    opt.add_option('--use-fftw', action='store_true', default=False,
                   help='use fftw instead of mkl')
    opt.add_option('--use-ctfft', action='store_true', default=False,
                   help='use internal ctfft instead of mkl')

def configure(cnf):
    cnf.load('compiler_c compiler_fc')

    # Find Proj
    cnf.check_cc(header_name='proj_api.h',uselib_store='proj', lib='proj')

    # Find GMT
    includedirs=['','/usr/include/gmt']
    found_gmt=False
    for inc in includedirs:
        try:
            cnf.check_cc(msg="Checking for gmt.h in '" + inc + "'",
                         header_name='gmt.h', includes=inc,
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
    for flag in ['-fopenmp','-openmp','-mp','-xopenmp','-omp','-qsmp=omp']:
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

    # Find FFTW or IMKL
    if cnf.options.use_fftw:
        frag="program main\n" + 'include "fftw3.f"\n' \
            + "end program main\n"
        cnf.check_fc(msg="Checking for fftw in /usr/include",
                     includes=['/usr/include'], fragment=frag, uselib_store='fftw',
                     lib=['fftw3f','fftw3f_threads'], define_name="FFTW3")
    elif not cnf.options.use_ctfft:
        cnf.check_fc(lib=['mkl_intel_lp64', 'mkl_intel_thread',
                          'mkl_core'], uselib_store='imkl',
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
