from pydsm.spc.stf import SourceTimeFunction
from pydsm import root_resources
import os
import pickle

class STFCatalog:
    '''Utility class to read and save catalog (dict) of
    source time functions with entries:
        event_id (str): SourceTimeFunction
    '''

    @staticmethod
    def read_from_file(path):
        '''Read source time functions from a file.
        Args:
            path (str): path of a file containing source time function
            informations. The first 3 columns of each row are:
                (1) event_id (GCMT id); (2) duration (s);
                (3) amplitude correction
            remaining columns will be ignored. Columns starting with '#'
            will be ignored.
        Returns:
            catalog (dict): source time function catalog
        '''
        with open(path, 'r') as f:
            lines = f.readlines()
        catalog = {
            ss[0].strip():SourceTimeFunction(
                'triangle', float(ss[1]), float(ss[2]))
            for line in lines
            for ss in line.split()}
        return catalog

    @staticmethod
    def _format_line(event_id, duration, amp_corr, *args):
        line = '{} {} {}'.format(event_id, duration, amp_corr)
        if args:
            for arg in args:
                line += ' {}'.format(arg)
        return line

    @staticmethod
    def read_scardec():
        '''Returns the scardec source time function catalog.
        Returns:
            catalog (dict): source time function catalog
        '''
        path = os.path.join(root_resources, 'scardec.pkl')
        catalog = STFCatalog.load(path)
        return catalog

    @staticmethod
    def save(path, catalog):
        '''Save catalog using pickle.dump().
        Args:
            path (str): name of the output file
            catalog (dict): source time function catalog
        '''
        with open(path, 'wb') as f:
            pickle.dump(catalog, f)

    @staticmethod
    def load(path):
        '''Read path into a catalog using pickle.load().
        Args:
            path (str): name of the file that contains a catalog
        Return:
            output (PyDSMOutput)
        '''
        with open(path, 'rb') as f:
            output = pickle.load(f)
        return output