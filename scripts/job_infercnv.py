from call_cnv.callCNV_workflow import CNVCallPlotWorkflow

def main():
    dir1 = "/stomics_data/liminData/Visium/stmut_python"
    output_dir = "/stomics_data/liminData/Visium/stmut_python/two_steps_out"

    graph_based_csv = f"{dir1}/BD17_bin8_inputs/Graph-Based.csv"
    annotate_csv = f"{dir1}/BD17_bin8_inputs/annotate.csv"
    bulk_csv = f"{output_dir}/tables/bulkCNV.csv"
    
    ncluster = 6
    pmtimes = 5

    # Instantiate workflow
    workflow = CNVCallPlotWorkflow(graph_based_csv, annotate_csv, output_dir)

    # Run the full workflow
    workflow.cnv_workflow(
        ncluster=ncluster,
        bulk_csv=bulk_csv,
        pmtimes=pmtimes
    )

if __name__ == '__main__':
    main()


