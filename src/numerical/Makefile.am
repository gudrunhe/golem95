noinst_LTLIBRARIES=libgolem95_numerical.la

AM_FCFLAGS=-I$(top_builddir)/src/module -I$(top_builddir)/avh_olo-3.6.1

libgolem95_numerical_la_SOURCES= mod_adapt_gauss.f90 mod_numeric.f90

nodist_pkginclude_HEADERS= adapt_gauss.mod numerical_evaluation.mod
CLEANFILES=*.mod

include Makefile.dep
