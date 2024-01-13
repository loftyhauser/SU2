# SU^2 autoconf file
# Stanford University Aerospace Design Lab
#
# Created: 		10/26/12
# Last updated: 4/4/2013
#
# Created by: 		M. Colonno  10/26/12
# Updated by: 		M. Colonno  3/16/13
#                   T. Lukaczyk 4/4/13
#
##########################

AC_INIT([SU^2], [2.0.7], [susquared-dev@lists.stanford.edu],
        [SU^2], [http://su2.stanford.edu])
AC_PREREQ([2.59])
AM_INIT_AUTOMAKE([subdir-objects])  # [SU^2], [2.0]

# optional bindings
AC_ARG_WITH(MPI,
    	AS_HELP_STRING([--with-MPI[=ARG]], [Build parallel version with MPI tools, ARG = MPI C++ compiler]), 
	[with_MPI=$withval], [with_MPI="no"])
AC_ARG_ENABLE(CFD,
   	AS_HELP_STRING([--disable-CFD], [build the SU2_CFD executable (default = yes)]), 
	[build_CFD=$enableval], [build_CFD="yes"])

##########################

have_MPI="no"
if test "$with_MPI" != "no"
then
	CXX="$with_MPI"	
	have_MPI="yes"
else
	CXXFLAGS=""$CXXFLAGS" -DNO_MPI"
fi

###########################

AM_CONDITIONAL([BUILD_CFD],[test $build_CFD != "no"])

###########################

# AC_CONFIG_HEADERS([config.h])
AC_PROG_CXX
# AC_PROG_CC
# AM_PROG_CC_C_O
AC_PROG_RANLIB
AC_OUTPUT(Makefile \
	Common/lib/Makefile \
	SU2_CFD/obj/Makefile \
	SU2_PY/Makefile)

abs_srcdir=`(cd $srcdir && pwd)`
abs_prefix=`(cd $prefix && pwd)`


AC_MSG_RESULT([

Build Configuration Summary:

    Source code location: ${abs_srcdir}
    Install location: 	  ${abs_prefix}
    Version:              ${VERSION}
    Compiler:             ${CXX}
    Preprocessor flags:   ${CPPFLAGS}
    Compiler flags:       ${CXXFLAGS}
    Linker flags:         ${LDFLAGS}
    MPI support:          $have_MPI

    Build SU2_CFD:        $build_CFD

Please be sure to add the [\$]SU2_HOME and [\$]SU2_RUN environment variables, and update your [\$]PATH with [\$]SU2_RUN

Based on the input to this configuration, add these lines to your .bashrc file -

export SU2_RUN="${abs_prefix}/bin"
export SU2_HOME="${abs_srcdir}"
export PATH=[\$]PATH:[\$]SU2_RUN


])
