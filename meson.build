project(
    'golem95', 
    'fortran', 
    version: '1.3.4',
    default_options: ['buildtype=release']
    )

pkg = import('pkgconfig')

add_project_arguments('-w', language: 'fortran')

deps = [dependency('avh_olo')]

golem95_src = files([])

config_data = configuration_data()

if get_option('precision') == 'intermediate'
    config_data.set('fortran_real_kind', 'selected_real_kind(18,4931)')
elif get_option('precision') == 'quadruple'
    config_data.set('fortran_real_kind', 'selected_real_kind(32,50)')
else
    config_data.set('fortran_real_kind', 'kind(1.0d0)')
endif

config_data.set('case_with_lt', '!')
config_data.set('lt_real_kind', 'kind(1.0d0)')

subdir('src')

golem95 = both_libraries('golem95',
    golem95_src,
    dependencies: deps,
    install: true,
    install_rpath:  get_option('prefix') / get_option('libdir')
    )

golem95_dep = declare_dependency(link_with: golem95)

meson.add_install_script(
  find_program(files('install_mod_files.py')),
  meson.current_build_dir(),
  get_option('includedir') / 'golem95',
)

meson.override_dependency('golem95', golem95_dep)

pkg.generate(golem95,
    subdirs: 'golem95',
    description: 'A package for the numerical evaluation of integrals and tensor form factors with up to six external legs'
    )