""" 3. HDF5 Matrix Reader """
import h5py

class H5MatrixReader:
    """ Handles reading matrix data from HDF5 files """
    def __init__(self, h5_file):
        self.h5_file = h5_file
    
    def read_matrix_data(self):
        """Extract matrix data and meta data from HDF5 file """
        with h5py.File(self.h5_file, 'r') as f:
            matrix = f['matrix']

            return {
                'data': matrix['data'][:],
                'indices': matrix['indices'][:],
                'indptr': matrix['indptr'][:],
                'shape': matrix['shape'][:],
                'barcodes': matrix['barcodes'][:].astype(str),
                'genes':self._extract_geneid(matrix)
            }

    def _extract_geneid(self, matrix):
        """ Extrac ensembl ID given gene name """
        if 'name' in matrix['features']:
            return matrix['features']['id'][:].astype(str)
        else:
            raise ValueError("Can't find gene names in 'features'.")
