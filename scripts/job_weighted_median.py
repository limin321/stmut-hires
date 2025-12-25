import os
import logging

from weighted_median.parallel_wtmedian import WorkflowOrchestrator

def main():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    logger = logging.getLogger(__name__)

    output_dir = "/stomics_data/liminData/Visium/stmut_python/two_steps_out/test"
    cores = 4

    orchestrator = WorkflowOrchestrator(output_dir, cores)
    orchestrator.run_parallel_workflow()

if __name__ == '__main__':
    main()
