noinst_LTLIBRARIES=libgolem95_module.la
libgolem95_module_la_SOURCES= precision_golem.f90 \
	tri.f90 array.f90 parametre.f90 constante.f90 \
	s_matrix_type.f90 equal.f90 multiply_div.f90 \
	sortie_erreur.f90 cache.f90 \
	form_factor_type.f90 kronecker.f90 z_log.f90\
	spinor.f90 translate.f90 zdilog.f90 cache_generic.f90
nodist_pkginclude_HEADERS= precision_golem.mod \
	tri_croissant.mod array.mod \
	cache.mod constante.mod s_matrix_type.mod dilogarithme.mod \
	equal.mod form_factor_type.mod kronecker.mod logarithme.mod \
	multiply_div.mod parametre.mod \
	sortie_erreur.mod spinor.mod translate.mod  cache_generic.mod 
AM_FCFLAGS= \
	-I$(top_builddir)/avh_olo-3.6.1

CLEANFILES=$(nodist_pkginclude_HEADERS)

include Makefile.dep
