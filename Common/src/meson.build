common_include = include_directories('./')
common_src =files(['graph_coloring_structure.cpp',
           'CConfig.cpp',
           '../include/parallelization/mpi_structure.cpp'])

subdir('linear_algebra')
subdir('toolboxes')
subdir('geometry')
subdir('geometry/elements')
subdir('geometry/dual_grid')
subdir('geometry/primal_grid')
subdir('geometry/meshreader')
subdir('adt')

if get_option('enable-normal')
  common = static_library('SU2Common', 
                          common_src,
                          install : false, 
                          dependencies : su2_deps, 
                          cpp_args: [default_warning_flags, su2_cpp_args])

  common_dep = declare_dependency(link_with: common,
                                  include_directories : common_include)
endif
