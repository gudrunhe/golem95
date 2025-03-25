# vim: ts=3:sw=3
DIST_SUBDIRS = demos src avh_olo-3.6.1
SUBDIRS = avh_olo-3.6.1 src

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = golem95.pc

lib_LTLIBRARIES         = libgolem.la libavh_olo.la
libgolem_la_SOURCES     =
libavh_olo_la_SOURCES   =
libavh_olo_la_LIBADD    = \
       $(builddir)/avh_olo-3.6.1/libavh_olo.la
libgolem_la_LIBADD      = \
	$(builddir)/src/form_factor/libgolem95_formfactor.la \
	$(builddir)/src/integrals/four_point/libgolem95_integrals_four_point.la \
	$(builddir)/src/integrals/three_point/libgolem95_integrals_three_point.la \
	$(builddir)/src/integrals/two_point/libgolem95_integrals_two_point.la \
	$(builddir)/src/integrals/one_point/libgolem95_integrals_one_point.la \
	$(builddir)/src/numerical/libgolem95_numerical.la \
	$(builddir)/src/kinematic/libgolem95_kinematics.la \
	$(builddir)/src/module/libgolem95_module.la \
	$(builddir)/src/interface/libgolem95_interface.la \
   $(builddir)/avh_olo-3.6.1/libavh_olo.la \
   $(builddir)/src/higherrank/libgolem95_higherrank.la

EXTRA_DIST = \
	$(srcdir)/demos/DemoContents \
	$(srcdir)/demos/plot_demo_A55.gp \
	$(srcdir)/demos/plotdetG_Im.gp \
	$(srcdir)/demos/plotdetG_Re.gp \
	$(srcdir)/demos/README \
	$(srcdir)/demos/table_of_results_3point_option1.txt \
	$(srcdir)/demos/table_of_results_3point_option2.txt \
	$(srcdir)/demos/table_of_results_3point_option3.txt \
	$(srcdir)/demos/table_of_results_3point_option4.txt \
	$(srcdir)/demos/table_of_results_3point_option5.txt \
	$(srcdir)/demos/table_of_results_3point_option6.txt \
	$(srcdir)/demos/table_of_results_3point_option7.txt \
	$(srcdir)/demos/table_of_results_4point_option1.txt \
	$(srcdir)/demos/table_of_results_4point_option2.txt \
	$(srcdir)/demos/table_of_results_4point_option3.txt \
	$(srcdir)/demos/table_of_results_4point_option4.txt \
	$(srcdir)/demos/table_of_results_4point_option5.txt \
	$(srcdir)/demos/table_of_results_4point_option6.txt \
	$(srcdir)/demos/table_of_results_4point_option7.txt \
	$(srcdir)/demos/table_of_results_4point_option8.txt \
	$(srcdir)/demos/table_of_results_5point.txt \
	$(srcdir)/demos/table_of_results_6point.txt \
	$(srcdir)/demos/table_of_results_LT.txt \
	$(srcdir)/demos/table_of_results_4photon.txt \
	$(srcdir)/demos/demo_a55_dets_sing.txt \
	$(srcdir)/demos/demo_detG.txt \
	$(srcdir)/doc/robodoc.rc \
	$(srcdir)/doc/html/numerical/*.html \
	$(srcdir)/doc/html/integrals/one_point/*.html \
	$(srcdir)/doc/html/integrals/two_point/*.html \
	$(srcdir)/doc/html/integrals/three_point/*.html \
	$(srcdir)/doc/html/integrals/four_point/*.html \
	$(srcdir)/doc/html/form_factor/*.html \
	$(srcdir)/doc/html/kinematic/*.html \
	$(srcdir)/doc/html/module/*.html \
	$(srcdir)/doc/html/interface/*.html \
	$(srcdir)/doc/html/higherrank/*.html \
	$(srcdir)/doc/html/*.html \
	$(srcdir)/doc/html/robodoc.css \
	$(srcdir)/tool/highrank/gen_form_factors.py \
	$(srcdir)/tool/tens_rec/robodoc.txt \
	$(srcdir)/tool/tens_rec/tens.py \
	makedepend.perl \
	autogen.sh \
	#$(srcdir)/test/cleanup.pl \
	#$(srcdir)/test/Makefile_proto.in \
	#$(srcdir)/test/makenumerator.pl \
	#$(srcdir)/test/maketest.pl \
	#$(srcdir)/test/mask_3point \
	#$(srcdir)/test/mask_4point \
	#$(srcdir)/test/mask_5point \
	#$(srcdir)/test/mask_6point \
	#$(srcdir)/test/momenta.dat \
	#$(srcdir)/test/mom_rambo.f \
	#$(srcdir)/test/param.input \
	#$(srcdir)/test/Readme


ACLOCAL_AMFLAGS = -I m4
EXTRA_DIST+=m4

dist-hook:
	find $(distdir) -name '*.mod' -exec rm \{} \;
