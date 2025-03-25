noinst_LTLIBRARIES=libgolem95_kinematics.la

AM_FCFLAGS=-I$(top_builddir)/src/module -I$(top_builddir)/avh_olo-3.6.1

libgolem95_kinematics_la_SOURCES= inverse_matrice.f90 matrice_s.f90

nodist_pkginclude_HEADERS= matrice_s.mod inverse_matrice.mod
CLEANFILES=*.mod

include Makefile.dep
