noinst_LTLIBRARIES=libgolem95_formfactor.la

AM_FCFLAGS= \
	-I $(top_builddir)/src/module \
	-I $(top_builddir)/src/kinematic \
	-I $(top_builddir)/src/numerical \
	-I $(top_builddir)/src/interface \
	-I $(top_builddir)/src/integrals/one_point \
	-I $(top_builddir)/src/integrals/two_point \
	-I $(top_builddir)/src/integrals/three_point \
	-I $(top_builddir)/src/integrals/four_point \
	-I $(top_builddir)/avh_olo-3.6.1

libgolem95_formfactor_la_SOURCES= form_factor_1p.f90 form_factor_2p.f90 \
	form_factor_3p.f90 form_factor_4p.f90 form_factor_5p.f90 \
	form_factor_6p.f90
nodist_pkginclude_HEADERS= form_factor_5p.mod form_factor_6p.mod form_factor_4p.mod \
	form_factor_3p.mod form_factor_2p.mod form_factor_1p.mod
CLEANFILES=*.mod

include Makefile.dep
