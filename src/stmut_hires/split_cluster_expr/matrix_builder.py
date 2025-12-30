"""4. Expression Matrix Builder  """
import scipy.sparse

class ExpressionMatrixBuilder:
    """ Handles building and managing expression matrix """

    @staticmethod
    def build_sparse_matrix(matrix_data):
        """Construct sparse matrix from h5_file"""
        return scipy.sparse.csc_matrix(
            (matrix_data['data'], matrix_data['indices'], matrix_data['indptr']),
            shape = matrix_data['shape']
        )
    
    @staticmethod
    def build_barcodes_index_map(matrix_data):
        """Create barcode to index mapping """
        return {bc: i for i, bc in enumerate(matrix_data['barcodes'])}

    @staticmethod
    def extract_cluster_expression(barcode_index_map, expr_matrix, all_barcodes, cluster_barcodes):
        """Extract expression data for specific cluster barcodes """
        valid_indices = [
            barcode_index_map[bc] for bc in cluster_barcodes 
            if bc in barcode_index_map
        ]

        if not valid_indices:
            return None, None
        sub_expr = expr_matrix[:, valid_indices].toarray()  # shape: [genes x barcodes]
        sub_barcodes = [all_barcodes[i] for i in valid_indices]
        return sub_expr, sub_barcodes
