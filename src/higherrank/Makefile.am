noinst_LTLIBRARIES=libgolem95_higherrank.la

AM_FCFLAGS=\
	-I$(top_builddir)/src/module \
	-I$(top_builddir)/src/kinematic \
	-I$(top_builddir)/src/numerical \
	-I$(top_builddir)/src/form_factor \
	-I$(top_builddir)/src/integrals/one_point \
	-I$(top_builddir)/src/integrals/two_point \
	-I$(top_builddir)/src/integrals/three_point \
	-I$(top_builddir)/src/integrals/four_point \
	-I$(top_builddir)/avh_olo-3.6.1

libgolem95_higherrank_la_SOURCES=generic_function_np.f90 form_factor_higher_ranks.f90

nodist_pkginclude_HEADERS= generic_function_np.mod form_factor_higher_ranks.mod

CLEANFILES=*.mod

include Makefile.dep
