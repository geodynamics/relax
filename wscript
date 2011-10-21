def check_fc_header(cnf,header_name,incs,libname,link_libs, defname):
    frag="program main\n" + 'include "' + header_name + '"\n' \
        + "end program main\n"
    cnf.check_fc(msg="Checking for " + header_name + " in " + incs,
                 includes=incs, fragment=frag, uselib_store=libname,
                 lib=link_libs, define_name=defname)

def options(opt):
    opt.load('compiler_c compiler_fc')
    opt.add_option('--use-fftw', action='store_true', default=False,
                   help='use fftw instead of mkl')
    opt.add_option('--use-ctfft', action='store_true', default=False,
                   help='use internal ctfft instead of mkl')

def configure(cnf):
    cnf.load('compiler_c compiler_fc')
    cnf.check_cc(header_name='proj_api.h',uselib_store='proj', lib='proj')
    includedirs=['','/usr/include/gmt']
    found_gmt=False
    for inc in includedirs:
        try:
            cnf.check_cc(msg="Checking for gmt.h in '" + inc + "'",
                         header_name='gmt.h', includes=inc,
                         lib=['gmt','netcdf'], uselib_store='gmt')
        except:
            pass
        else:
            found_gmt=True
            break
    if not found_gmt:
        cnf.fatal('Could not find gmt')
    try:
        cnf.check_fc(fcflags='-openmp', linkflags='-openmp', uselib_store='openmp')
    except:
        cnf.check_fc(fcflags='-fopenmp', linkflags='-fopenmp', uselib_store='openmp')

    if cnf.options.use_fftw:
        check_fc_header(cnf,"fftw3.f",'/usr/include','fftw',
                        ['fftw3f','fftw3f_threads'],"FFTW3")
    elif not cnf.options.use_ctfft:
        cnf.check_fc(lib=['mkl_intel_lp64', 'mkl_intel_thread',
                          'mkl_core'], uselib_store='imkl',
                     use='openmp', define_name='IMKL_FFT')

    cnf.write_config_header('config.h')

def build(bld):
    sources=['relax.f90',
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
             'input.f90']
    uses=['gmt','proj','openmp']
    if bld.options.use_fftw:
        uses.append('fftw')
    elif not bld.options.use_ctfft:
        uses.append('imkl')
        sources.append('mkl_dfti.f90')
    bld.program(features='c fc fcprogram',
                source=sources,
                includes=['build'],
                use=uses,
                fcflags=['-cpp','-zero'],
                target='relax'
                )
