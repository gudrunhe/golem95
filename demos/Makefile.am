EXTRA_PROGRAMS= demo_3point demo_4photon demo_4point demo_5point \
	demo_6point demo_a55_dets_sing demo_detg demo_LT demo_cmplx_masses \
	SusyLandau

if COMPILE_TENSREC
EXTRA_PROGRAMS+=demo_tens_rec
demo_tens_rec_SOURCES=demo_tens_rec.f90 demo_tens_mod.f90
demo_tens_rec_LDADD=${link_to_golem95}
endif

LIBAVHOLO=${libdir}/libavh_olo.a
link_to_golem95=${libdir}/lib@GENERIC_LIBRARY_NAME@.a ${LIBAVHOLO} ${LIBLOOPTOOLS}
AM_FCFLAGS=-I${pkgincludedir}

demo_3point_SOURCES=demo_3point.f90
demo_4photon_SOURCES=demo_4photon.f90
demo_4point_SOURCES=demo_4point.f90
demo_5point_SOURCES=demo_5point.f90
demo_6point_SOURCES=demo_6point.f90
demo_a55_dets_sing_SOURCES=demo_a55_dets_sing.f90
demo_detg_SOURCES=demo_detg.f90
demo_LT_SOURCES=demo_LT.f
demo_cmplx_masses_SOURCES=demo_cmplx_masses.f90
SusyLandau_SOURCES=SusyLandau.f90

demo_3point_LDADD=${link_to_golem95}
demo_4photon_LDADD=${link_to_golem95}
demo_4point_LDADD=${link_to_golem95}
demo_5point_LDADD=${link_to_golem95}
demo_6point_LDADD=${link_to_golem95}
demo_a55_dets_sing_LDADD=${link_to_golem95}
demo_detg_LDADD=${link_to_golem95}
demo_LT_LDADD=${link_to_golem95}
demo_cmplx_masses_LDADD=${link_to_golem95}
SusyLandau_LDADD=${link_to_golem95}

all-local: ${EXTRA_PROGRAMS}

demo_tens_rec.o: demo_tens_mod.o
demo_tens_rec.lo: demo_tens_mod.lo
demo_tens_rec.obj: demo_tens_mod.obj
