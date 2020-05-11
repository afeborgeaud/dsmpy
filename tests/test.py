from pydsm import dsm, rootdsm_psv

parameter_file = rootdsm_psv + 'test1.inf'
print('Read inputs')
inputs = dsm.PyDSMInput.input_from_file(parameter_file)
print('Compute')
outputs = dsm.compute(inputs)
