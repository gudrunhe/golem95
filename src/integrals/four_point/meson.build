golem95_src += files(['function_4p1m.f90', 'function_4p2m_3mi_onshell.f90', 'function_4p2m_adj.f90',
    'function_4p2m_opp.f90', 'function_4p3m.f90', 'function_4p4m.f90', 'function_4p_ql10.f90',
    'function_4p_ql11.f90', 'function_4p_ql12.f90', 'function_4p_ql13.f90', 'function_4p_ql14.f90',
    'function_4p_ql15.f90', 'function_4p_ql16.f90', 'function_4p_ql6.f90', 'function_4p_ql7.f90',
    'function_4p_ql8.f90', 'function_4p_ql9.f90'])

golem95_src += configure_file(
    input: 'generic_function_4p.f90.in',
    output: 'generic_function_4p.f90',
    configuration: config_data
)