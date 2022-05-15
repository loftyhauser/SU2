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

# check for non-debug build
if get_option('buildtype')!='debug'
  su2_cpp_args += '-DNDEBUG'
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
