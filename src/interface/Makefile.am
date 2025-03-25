noinst_LTLIBRARIES=libgolem95_interface.la

AM_FCFLAGS= \
	-I$(top_builddir)/src/module \
	-I$(top_builddir)/src/kinematic \
	-I$(top_builddir)/src/form_factor \
	-I$(top_builddir)/src/integrals/one_point \
	-I$(top_builddir)/src/integrals/two_point \
	-I$(top_builddir)/src/integrals/three_point \
	-I$(top_builddir)/src/integrals/four_point \
	-I$(top_builddir)/src/higherrank \
	-I$(top_builddir)/src/numerical \
	-I$(top_builddir)/avh_olo-3.6.1

libgolem95_interface_la_SOURCES= \
	tool_lt_to_golem.f90 tensor_integrals.f90 \
	gb0.f90 gc0.f90 gd0.f90 ge0.f90 gf0.f90 ga0.f90

nodist_pkginclude_HEADERS= tool_lt_to_golem.mod tensor_integrals.mod
if COMPILE_TENSREC
libgolem95_interface_la_SOURCES+=tens_rec.f90 tens_comb.f90
nodist_pkginclude_HEADERS+=tens_rec.mod tens_comb.mod
endif

CLEANFILES=*.mod

include Makefile.dep
