noinst_LTLIBRARIES=libgolem95_integrals_one_point.la

AM_FCFLAGS= \
	-I$(top_builddir)/src/module \
	-I$(top_builddir)/src/kinematic \
	-I$(top_builddir)/src/numerical \
	-I$(top_builddir)/avh_olo-3.6.1

libgolem95_integrals_one_point_la_SOURCES= generic_function_1p.f90
libgolem95_integrals_one_point_la_FCFLAGS= $(AM_FCFLAGS)

nodist_pkginclude_HEADERS= generic_function_1p.mod
CLEANFILES=*.mod

include Makefile.dep
