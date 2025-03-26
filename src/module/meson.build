golem95_src += files(['tri.f90', 'array.f90', 'parametre.f90', 'constante.f90',
    's_matrix_type.f90', 'equal.f90', 'multiply_div.f90', 'sortie_erreur.f90', 'cache.f90',
    'form_factor_type.f90', 'kronecker.f90', 'z_log.f90', 'spinor.f90', 'translate.f90', 'zdilog.f90',
    'cache_generic.f90'])

golem95_src += configure_file(
    input:'precision_golem.f90.in', 
    output: 'precision_golem.f90',
    configuration: config_data
    )