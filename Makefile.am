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
	meson.build \
	meson_options.txt \
	install_mod_files.py \
	$(srcdir)/src/meson.build \
	$(srcdir)/src/form_factor/meson.build \
	$(srcdir)/src/higherrank/meson.build \
	$(srcdir)/src/integrals/meson.build \
	$(srcdir)/src/integrals/four_point/meson.build \
	$(srcdir)/src/integrals/one_point/meson.build \
	$(srcdir)/src/integrals/three_point/meson.build \
	$(srcdir)/src/integrals/two_point/meson.build \
	$(srcdir)/src/interface/meson.build \
	$(srcdir)/src/kinematic/meson.build \
	$(srcdir)/src/module/meson.build \
	$(srcdir)/src/numerical/meson.build \
	$(srcdir)/subprojects \
	$(srcdir)/subprojects/avh_olo.wrap \
	$(srcdir)/subprojects/packagefiles/avh_olo/create.py.patch \
	$(srcdir)/subprojects/packagefiles/avh_olo/install_mod_files.py \
	$(srcdir)/subprojects/packagefiles/avh_olo/meson.build 

ACLOCAL_AMFLAGS = -I m4
EXTRA_DIST+=m4

dist-hook:
	find $(distdir) -name '*.mod' -exec rm \{} \;
