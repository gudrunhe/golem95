noinst_LTLIBRARIES=libgolem95_integrals_three_point.la

AM_FCFLAGS= \
	-I$(top_builddir)/src/module \
	-I$(top_builddir)/src/kinematic \
	-I$(top_builddir)/src/numerical \
	-I$(top_builddir)/src/integrals/one_point \
	-I$(top_builddir)/src/integrals/two_point \
	-I$(top_builddir)/avh_olo-3.6.1

libgolem95_integrals_three_point_la_SOURCES= \
	mod_gn.f90 mod_h0.f90 mod_he.f90 mod_hf.f90 \
	function_3p0m_1mi.f90 function_3p1m_1mi.f90 function_3p1m_2mi.f90 \
	function_3p1m.f90 function_3p2m_1mi.f90 function_3p2m.f90 \
	function_3p3m.f90 function_3p_finite.f90 \
	generic_function_3p.f90

nodist_pkginclude_HEADERS= function_3p3m.mod generic_function_3p.mod \
	function_3p1m_2mi.mod func_h0.mod func_he.mod \
	function_3p2m_1mi.mod func_hf.mod function_3p_finite.mod \
	func_gn.mod function_3p0m_1mi.mod function_3p1m_1mi.mod \
	function_3p2m.mod function_3p1m.mod
CLEANFILES=*.mod

include Makefile.dep
