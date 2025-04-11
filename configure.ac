AC_INIT

dnl -----------------------------------------------
dnl Package name and version number (user defined)
dnl -----------------------------------------------

GENERIC_LIBRARY_NAME=golem
GENERIC_PACKAGE_NAME=golem95

#release versioning
GENERIC_MAJOR_VERSION=1
GENERIC_MINOR_VERSION=3
GENERIC_MICRO_VERSION=4

#API version (often = GENERIC_MAJOR_VERSION.GENERIC_MINOR_VERSION)
GENERIC_API_VERSION=$GENERIC_MAJOR_VERSION.$GENERIC_MINOR_VERSION
AC_SUBST(GENERIC_API_VERSION)

#shared library versioning
GENERIC_LIBRARY_VERSION=2:0:1
#                       | | |
#                +------+ | +---+
#                |        |     |
#             current:revision:age
#                |        |     |
#                |        |     +- increment if interfaces have been added
#                |        |        set to zero if interfaces have been removed
#                                  or changed
#                |        +- increment if source code has changed
#                |           set to zero if current is incremented
#                +- increment if interfaces have been added, removed or changed

dnl --------------------------------
dnl Package name and version number
dnl --------------------------------

AC_SUBST(GENERIC_LIBRARY_VERSION)

PACKAGE=$GENERIC_PACKAGE_NAME
AC_SUBST(GENERIC_LIBRARY_NAME)

GENERIC_VERSION=$GENERIC_MAJOR_VERSION.$GENERIC_MINOR_VERSION.$GENERIC_MICRO_VERSION
AC_SUBST(GENERIC_VERSION)

VERSION=$GENERIC_VERSION

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR([config.aux])

AM_INIT_AUTOMAKE($PACKAGE, $VERSION, no-define, no-dependencies, tar-ustar)
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES])

dnl AC_DISABLE_STATIC([])
LT_INIT

dnl -----------------------------------------------
dnl Checks for programs.
dnl -----------------------------------------------

AC_LANG_PUSH([Fortran])
AC_PROG_F77([gfortran ifort g77 f77 xlf frt pgf77 fort77 fl32 af77])
AC_PROG_FC([gfortran ifort g77 f77 xlf frt pgf77 fort77 fl32 af77])

dnl When configuring for `make dist' purposes, skip checks that may yield fatal
dnl errors when there is no working F90 compiler.
if test -z "$with_dist"; then
   dnl Check the flags needed to link fc programs with ld (i.e. cc)
   AC_FC_LIBRARY_LDFLAGS
   dnl Check for underscoring of external names
   AC_FC_WRAPPERS
   dnl We need to use .f90 and not .f to enable Automake FC support
   dnl Some Intel ifc/ifort do not understand .f95.  :-/
   AC_FC_SRCEXT([f90])
   AC_FC_FREEFORM
fi
AC_LANG_POP

LT_INIT

AM_SANITY_CHECK

dnl AC_PATH_PROG(PERL,perl,perl)

dnl ----------------------------------------------
dnl tensorial reconstruction interface
dnl ----------------------------------------------
dnl user should be able to switch off tensorial reconstruction
AC_ARG_ENABLE([tensrec],
  [AS_HELP_STRING([--disable-tensrec],
             [switch off tensorial reconstruction interface])],
  [case "${enableval}" in
     yes) with_tensrec=yes ;;
     no)  with_tensrec=no ;;
     *) AC_MSG_ERROR(bad value ${enableval} for --enable-tensrec) ;;
   esac],
  [with_tensrec=yes])

AM_CONDITIONAL(COMPILE_TENSREC, [test "x$with_tensrec" != xno])

dnl -----------------------------------------------
dnl link LoopTools
dnl -----------------------------------------------
AC_ARG_WITH([looptools],
  [AS_HELP_STRING([--with-looptools],
  [enable linking to LoopTools])],
  [],
  [with_looptools=no])

LIBLOOPTOOLS=
AS_IF(
  [test "x$with_looptools" == xyes],
  [AC_CHECK_LIB([ooptools -lgfortran], [ltexi_],
     [AC_SUBST([LIBLOOPTOOLS], ["-looptools"])
      AC_DEFINE([HAVE_LT], [1], [Define if you have Looptools])
     ],
     [AC_MSG_FAILURE(
        [--with-looptools was given, but test for -looptools failed. \
 Consider using --with-looptools=path/libooptools.a .])],
     [-looptools])],
   [test "x$with_looptools" != xno],
   [AC_CHECK_FILE([$with_looptools], 
     [AC_SUBST([LIBLOOPTOOLS], ["$with_looptools"])
      AC_DEFINE([HAVE_LT], [1], [Define if you have Looptools])
     ],
     [AC_MSG_FAILURE(
        [--with-looptools was given, but location '$with_looptools' \
 is wrong.])])]
  )

AS_IF([test "x$with_looptools" != xno],
		[AC_SUBST([case_with_lt], ["    "])],
		[AC_SUBST([case_with_lt], ["!AC!"])])
AS_IF([test "x$with_looptools" != xno],
		[AC_SUBST([case_wout_lt], ["!AC!"])],
		[AC_SUBST([case_wout_lt], ["    "])])

dnl -----------------------------------------------
dnl Set the precision used by LoopTools
dnl -----------------------------------------------
AC_ARG_WITH([lt-precision], [AS_HELP_STRING([--with-lt-precision],
	    [set the precision used by LoopTools to either
	     'double' or 'quadruple'. @<:@default=double@:>@])],
	     [], [with_lt_precision=double])

AS_IF(
	[test "x$with_lt_precision" == xquadruple],
		[AC_SUBST([lt_real_kind], ["selected_real_kind(32,50)"])],
	[test "x$with_lt_precision" == xquad],
		[AC_SUBST([lt_real_kind], ["selected_real_kind(32,50)"])],
	[test "x$with_lt_precision" == xdouble],
		[AC_SUBST([lt_real_kind], ["kind(1.0d0)"])],
	[AC_MSG_FAILURE([--with-lt-precision was given with an
	 unrecognized parameter])])

dnl -----------------------------------------------
dnl link AvH's OneLOop
dnl -----------------------------------------------
dnl AC_ARG_WITH([avh_olo],
dnl   [AS_HELP_STRING([--without-avh_olo],
dnl   [disable integrals from OneLOop])],
dnl   [],
dnl   [with_avh_olo=auto])
dnl             
dnl LIBAVHOLO=
dnl AS_IF([test "x$with_avh_olo" == xauto],
dnl   [AC_CHECK_LIB([avh_olo], [avh_olo_d0_],
dnl      [AC_SUBST([LIBAVHOLO], ["-avh_olo"])
dnl       AC_DEFINE([HAVE_AVH], [1], [Define if you have AvHOneLOop])
dnl      ],
dnl      [AC_MSG_FAILURE(
dnl         [No path specified for avh_olo and test for avh_olo failed. \
dnl  Consider using --with-avh_olo=path/libavh_olo.a \
dnl  or --without-avh_olo.])],
dnl      [-avh_olo])],
dnl   [test "x$with_avh_olo" == xyes],
dnl   [AC_CHECK_LIB([avh_olo], [avh_olo_d0_],
dnl      [AC_SUBST([LIBAVHOLO], ["-avh_olo"])
dnl       AC_DEFINE([HAVE_AVH], [1], [Define if you have AvHOneLOop])
dnl      ],
dnl      [AC_MSG_FAILURE(
dnl         [--with-avh_olo was given, but test for avh_olo failed. \
dnl  Consider using --with-avh_olo=path/libavh_olo.a .])],
dnl      [-avh_olo])],
dnl 
dnl   [test "x$with_avh_olo" != xno],
dnl   [AC_CHECK_FILE([$with_avh_olo], 
dnl      [AC_SUBST([LIBAVHOLO], ["$with_avh_olo"])
dnl       AC_DEFINE([HAVE_AVH], [1], [Define if you have AvHOneLOop])
dnl      ],
dnl      [AC_MSG_FAILURE(
dnl         [--with-avh_olo was given, but location '$with_avh_olo' \
dnl  is wrong.])])]
dnl   )
dnl 
dnl AS_IF([test "x$with_avh_olo" != xno],
dnl 		[AC_SUBST([case_with_avh], ["    "])],
dnl 		[AC_SUBST([case_with_avh], ["!AC!"])])
dnl AS_IF([test "x$with_avh_olo" != xno],
dnl 		[AC_SUBST([case_wout_avh], ["!AC!"])],
dnl 		[AC_SUBST([case_wout_avh], ["    "])])
dnl 
dnl -----------------------------------------------
dnl Set the precision
dnl -----------------------------------------------
AC_ARG_WITH([precision], [AS_HELP_STRING([--with-precision],
	    [set the precision of the library to either
	     'double' or 'quadruple'. @<:@default=double@:>@])],
	     [], [with_precision=double])

AS_IF(
	[test "x$with_precision" == xquadruple],
		[AC_SUBST([fortran_real_kind], ["selected_real_kind(32,50)"])],
	[test "x$with_precision" == xquad],
		[AC_SUBST([fortran_real_kind], ["selected_real_kind(32,50)"])],
	[test "x$with_precision" == xdouble],
		[AC_SUBST([fortran_real_kind], ["kind(1.0d0)"])],
	[test "x$with_precision" == xintermediate],
		[AC_SUBST([fortran_real_kind], ["selected_real_kind(18,4931)"])],
	[AC_MSG_FAILURE([--with-precision was given with an unrecognized
	 parameter])])

dnl -----------------------------------------------
dnl Set the precision of AVH's OneLOop
dnl -----------------------------------------------
AC_ARG_WITH([avh_olo_precision], [AS_HELP_STRING([--with-avh_olo-precision],
	    [set the precision of OneLOop to either
	     'double' or 'quadruple'. @<:@default=double@:>@])],
	     [], [with_avh_olo_precision=double])

AS_IF(
	[test "x$with_avh_olo_precision" == xquadruple],
		[AC_SUBST([avh_olo_real_kind],["selected_real_kind(32,50)"])],
	[test "x$with_avh_olo_precision" == xquad],
		[AC_SUBST([avh_olo_real_kind],["selected_real_kind(32,50)"])],
	[test "x$with_avh_olo_precision" == xdouble],
		[AC_SUBST([avh_olo_real_kind],["kind(1.0d0)"])],
	[test "x$with_precision" == xintermediate],
		[AC_SUBST([avh_olo_real_kind],["selected_real_kind(18,4931)"])],
	[AC_MSG_FAILURE([--with-avh_olo-precision
	 was given with an unrecognized parameter])])

dnl -----------------------------------------------
dnl Generates Makefile's, configuration files and scripts
dnl -----------------------------------------------

AC_CONFIG_FILES([Makefile \
   avh_olo-3.6.1/Makefile \
   avh_olo-3.6.1/avh_olo_xkind.f90 \
   src/Makefile \
   src/module/Makefile \
   src/module/precision_golem.f90 \
   src/integrals/Makefile \
   src/integrals/two_point/Makefile \
   src/integrals/one_point/Makefile \
   src/integrals/three_point/Makefile \
   src/integrals/four_point/Makefile \
   src/integrals/four_point/generic_function_4p.f90 \
   src/interface/Makefile \
   src/numerical/Makefile \
   src/kinematic/Makefile \
   src/form_factor/Makefile \
   src/higherrank/Makefile \
   demos/Makefile \
   golem95.pc])
dnl   test/Makefile_proto

AC_OUTPUT
