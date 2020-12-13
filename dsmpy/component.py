from enum import IntEnum

class Component(IntEnum):
    '''Seismic record geographical components (R T Z).
    
    '''
    Z = 0
    R = 1
    T = 2

    @staticmethod
    def parse_component(str):
        '''Parse component from a str.

            Args:
                str (str): str representation of a component
                    (e.g., 'Z', 'vertical' -> Z)

        '''
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
        else:
            return None