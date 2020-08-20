from enum import IntEnum

class Component(IntEnum):
    '''Seismic component computed by DSM.'''
    Z = 0
    R = 1
    T = 2

    @staticmethod
    def parse_component(str):
        if (
            str == 'Z'
            or str == 'vertical'
            or str.endswith('Z')):
            return Component.Z
        elif (
            str == 'R'
            or str == 'radial'
            or str.endswith('R')):
            return Component.R
        elif (
            str == 'T'
            or str == 'trnsvers'
            or str.endswith('T')):
            return Component.T