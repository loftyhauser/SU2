su2_cfd_include = include_directories('./')
su2_cfd_src = files(['python_wrapper_structure.cpp',
                     'CMarkerProfileReaderFVM.cpp'])

su2_cfd_src += files(['fluid/CFluidModel.cpp',
                      'fluid/CIdealGas.cpp',
                      'fluid/CPengRobinson.cpp',
                      'fluid/CVanDerWaalsGas.cpp'])

su2_cfd_src += files(['output/COutputFactory.cpp',
                      'output/CBaselineOutput.cpp',
                      'output/CFlowCompOutput.cpp',
                      'output/CFVMOutput.cpp',
                      'output/CFlowOutput.cpp',
                      'output/CMeshOutput.cpp',
                      'output/COutput.cpp',
                      'output/filewriter/CParallelDataSorter.cpp',
                      'output/filewriter/CFVMDataSorter.cpp',
                      'output/filewriter/CSurfaceFVMDataSorter.cpp',
                      'output/filewriter/CParallelFileWriter.cpp',
                      'output/filewriter/CParaviewFileWriter.cpp',
                      'output/filewriter/CParaviewBinaryFileWriter.cpp',
                      'output/filewriter/CTecplotFileWriter.cpp',
                      'output/filewriter/CTecplotBinaryFileWriter.cpp',
                      'output/filewriter/CCSVFileWriter.cpp',
                      'output/filewriter/CSTLFileWriter.cpp',
                      'output/filewriter/CSU2FileWriter.cpp',
                      'output/filewriter/CSU2BinaryFileWriter.cpp',
                      'output/filewriter/CParaviewXMLFileWriter.cpp',
                      'output/filewriter/CParaviewVTMFileWriter.cpp',
                      'output/filewriter/CSU2MeshFileWriter.cpp',
                      'output/tools/CWindowingTools.cpp'])

su2_cfd_src += files(['variables/CTransLMVariable.cpp',
                      'variables/CTurbVariable.cpp',
                      'variables/CScalarVariable.cpp',
                      'variables/CBaselineVariable.cpp',
                      'variables/CTurbSAVariable.cpp',
                      'variables/CMeshVariable.cpp',
                      'variables/CMeshBoundVariable.cpp',
                      'variables/CMeshElement.cpp',
                      'variables/CTurbSSTVariable.cpp',
                      'variables/CVariable.cpp',
                      'variables/CNSVariable.cpp',
                      'variables/CFlowVariable.cpp',
                      'variables/CEulerVariable.cpp',
                      'variables/CSobolevSmoothingVariable.cpp'])

su2_cfd_src += files(['solvers/CSolverFactory.cpp',
                      'solvers/CBaselineSolver.cpp',
                      'solvers/CEulerSolver.cpp',
                      'solvers/CNSSolver.cpp',
                      'solvers/CSolver.cpp',
                      'solvers/CTemplateSolver.cpp',
                      'solvers/CTransLMSolver.cpp',
                      'solvers/CTurbSolver.cpp',
                      'solvers/CTurbSASolver.cpp',
                      'solvers/CTurbSSTSolver.cpp'])

su2_cfd_src += files(['numerics/CNumerics.cpp',
                      'numerics/template.cpp',
                      'numerics/transition.cpp',
                      'numerics/flow/convection/roe.cpp',
                      'numerics/flow/convection/fvs.cpp',
                      'numerics/flow/convection/cusp.cpp',
                      'numerics/flow/convection/hllc.cpp',
                      'numerics/flow/convection/ausm_slau.cpp',
                      'numerics/flow/flow_diffusion.cpp',
                      'numerics/flow/flow_sources.cpp',
                      'numerics/continuous_adjoint/adj_convection.cpp',
                      'numerics/continuous_adjoint/adj_diffusion.cpp',
                      'numerics/continuous_adjoint/adj_sources.cpp',
                      'numerics/scalar/scalar_sources.cpp',
                      'numerics/CGradSmoothing.cpp'])

su2_cfd_src += files(['../include/numerics_simd/CNumericsSIMD.cpp'])

su2_cfd_src += files(['drivers/CDriver.cpp',
                      'drivers/CSinglezoneDriver.cpp',
                      'drivers/CDummyDriver.cpp'])

su2_cfd_src += files(['integration/CIntegration.cpp',
                      'integration/CIntegrationFactory.cpp',
                      'integration/CSingleGridIntegration.cpp',
                      'integration/CMultiGridIntegration.cpp',
                      'integration/CNewtonIntegration.cpp'])

su2_cfd_src += files(['iteration/CIteration.cpp',
                      'iteration/CIterationFactory.cpp',
                      'iteration/CFluidIteration.cpp'])

su2_cfd_src += files(['limiters/CLimiterDetails.cpp'])

if get_option('enable-normal')
  su2_cfd_lib = static_library('SU2core',
                               su2_cfd_src,
                               install : false,
		               dependencies : [su2_deps, common_dep],
		               cpp_args:  [default_warning_flags, su2_cpp_args])
  su2_cfd_dep = declare_dependency(link_with: su2_cfd_lib,
                                   include_directories: su2_cfd_include)
  su2_cfd = executable('SU2_CFD',
                       'SU2_CFD.cpp',
                       install : true,
                       dependencies : [su2_cfd_dep, su2_deps, common_dep],
                       cpp_args:  ['-fPIC'] + [default_warning_flags, su2_cpp_args])
endif

