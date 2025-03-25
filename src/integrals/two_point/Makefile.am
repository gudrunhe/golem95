noinst_LTLIBRARIES=libgolem95_integrals_two_point.la

AM_FCFLAGS= \
	-I$(top_builddir)/src/module \
	-I$(top_builddir)/src/kinematic \
	-I$(top_builddir)/src/numerical \
	-I$(top_builddir)/src/integrals/one_point \
	-I$(top_builddir)/avh_olo-3.6.1

libgolem95_integrals_two_point_la_SOURCES= \
	function_2p0m_1mi.f90 function_2p_m1m2.f90 generic_function_2p.f90

nodist_pkginclude_HEADERS= function_2p_m1m2.mod function_2p0m_1mi.mod \
	generic_function_2p.mod
CLEANFILES=*.mod

include Makefile.dep
