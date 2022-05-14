project('SU2', 'c', 'cpp',
        version: '7.3.1 "Blackbird"',
        license: 'LGPL2',
        default_options: ['buildtype=release',
                          'warning_level=3',
                          'c_std=c99',
                          'cpp_std=c++11'])


pymod = import('python')
python = pymod.find_installation()

su2_cpp_args = []
su2_deps     = [declare_dependency(include_directories: 'externals/CLI11')]

default_warning_flags = []
if build_machine.system() != 'windows'
  if meson.get_compiler('cpp').get_id() != 'intel'
    default_warning_flags += ['-Wno-empty-body']
  endif
  default_warning_flags += ['-Wno-unused-parameter',
                            '-Wno-deprecated-declarations']

endif

# meson script path
script_path = meson.source_root() / 'meson_scripts'

# mpi flag
mpi = false

# If custom mpi mode is enabled, include and library pathes for MPI have to be set manually to env variables
if get_option('custom-mpi')
  mpi_dep = []
  mpi = true
# Otherwise they are automatically determined
else
  mpi_dep = [dependency('mpi', language:'c', required : get_option('with-mpi')),
             dependency('mpi', language:'cpp', required : get_option('with-mpi'))]
  if mpi_dep[0].found() or mpi_dep[1].found()
    mpi = true
  endif
endif

# OpenMP
omp = get_option('with-omp')

if omp
  omp_dep = dependency('openmp', language:'cpp')
endif

# check for non-debug build
if get_option('buildtype')!='debug'
  su2_cpp_args += '-DNDEBUG'
endif

# check if MPI dependencies are found and add them
if mpi

  # add MPI dependency
  su2_deps += mpi_dep
  su2_cpp_args += '-DHAVE_MPI'

  # compile metis
  subdir('externals/metis')

  # add metis dependency
  su2_deps     += metis_dep
  su2_cpp_args += '-DHAVE_METIS'

  # compile parmetis
  subdir('externals/parmetis')

  # add parmetis dependency
  su2_deps     += parmetis_dep
  su2_cpp_args += '-DHAVE_PARMETIS'

endif

if omp
  # add OpenMP dependency
  su2_deps += omp_dep

endif

# compile common library
subdir('Common/src')
# compile SU2_CFD executable
subdir('SU2_CFD/src')

message('''-------------------------------------------------------------------------
         |    ___ _   _ ___                                                      |
         |   / __| | | |_  )   Release 7.3.1 "Blackbird"                         |
         |   \__ \ |_| |/ /                                                      |
         |   |___/\___//___|   Meson Configuration Summary                       |
         |                                                                       |
         -------------------------------------------------------------------------
         
         Please be sure to add the $SU2_HOME and $SU2_RUN environment variables,
         and update your $PATH (and $PYTHONPATH if applicable) with $SU2_RUN

         Based on the input to this configuration, add these lines to your .bashrc file:

         export SU2_RUN=@0@
         export SU2_HOME=@1@
         export PATH=$PATH:$SU2_RUN
         export PYTHONPATH=$PYTHONPATH:$SU2_RUN

         Use './ninja -C @2@ install' to compile and install SU2
'''.format(get_option('prefix')+'/bin', meson.source_root(),
           meson.build_root().split('/')[-1]))
