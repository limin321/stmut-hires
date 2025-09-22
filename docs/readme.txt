merge_gene_mem_efficient_csv_chunks.py
   input gene-expression is csv file;
   separate data to chunks so the computer can handle large dataset, but very slow.


merge_gene_parquet_optimized.py
    pros:
        support parallel grouping genes for clusters;
        Automate update window-size if total genes of all given windows < cutoff 1k genes, for example.
        The higher cutoff, the faster.
        If set num-threads = num-clusters, all clusters run at the same time.




group_barcodes.py
    Improved version of merge_gene_mem_efficient.py. 
    pros:
        Gene expression matrix is h5 file, no need to convert to csv.
    cons:
        Run for loop, very time consuming.
        Not automatically update window-size when the given (ex 100) window-size used up, the total gene counts still not reach the min gene count = 1000.



    (scvi-env) [lchen@localhost stmut-hires]$ python group_barcodes.py --help
    usage: group_barcodes.py [-h] --expression_file EXPRESSION_FILE --cluster_file CLUSTER_FILE --spatial_file SPATIAL_FILE --output_dir
                            OUTPUT_DIR

    Grouping barcodes within each cluster.

    optional arguments:
    -h, --help            show this help message and exit
    --expression_file EXPRESSION_FILE
                            Gene expression matrix
    --cluster_file CLUSTER_FILE
                            Barcodes cluster info.
    --spatial_file SPATIAL_FILE
                            Barcodes spatial coordinates.
    --output_dir OUTPUT_DIR
                            Output directory.
    
legacy/module/cluster_gene_exp.py
    input gene exp is csv.
    generate cluster gene expression file. 
    cons -- very slow.

