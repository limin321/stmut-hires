import json
import os
from datetime import datetime
from typing import Dict, List, Optional


class PipelineMetadataManager:
    """ 
    save metadata
    load metadata
    validate metadata
    manage pipeline metadata for reproducibility, provenance tracking, and downstream workflow recovery.
    """
    METADATA_FILENAME = "metadata.json"
    PIPELINE_VERSION = "1.0.0"


    def __init__(self, output_dir: str):
        self.output_dir = output_dir
        self.metadata_path = os.path.join(output_dir, self.METADATA_FILENAME)

    def initialize_metadata(
        self,
        platform: str,
        canonical_inputs: Dict,
        raw_inputs: Dict,
        runtime_params=None,
        completed_steps: Optional[List[str]] = None
    ) -> Dict:
        """ 
        Create initial metadata structure.
        """
        metadata = {
            "pipeline_version": self.PIPELINE_VERSION,
            "platform": platform,
            "created_at": datetime.utcnow().isoformat(),
            "canonical_inputs": canonical_inputs,
            "raw_inputs":raw_inputs,
            "runtime_params": runtime_params or {},
            "completed_steps": completed_steps or []
        }
        self.save_metadata(metadata)

        return metadata



    def save_metadata(self, metadata: Dict):
        """
        Save metadata JSON to disk.
        """
        os.makedirs(self.output_dir, exist_ok=True)

        with open(self.metadata_path, "w") as f:
            json.dump(
                metadata,
                f,
                indent=4
            )
    
    def load_metadata(self) -> Dict:
        """ 
        Load metadata JSON from disk
        """
        if not os.path.exists(self.metadata_path):
            raise FileNotFoundError(
                f"Metadata file not found: {self.metadata_path}"
            )
        
        with open(self.metadata_path, "r") as f:
            metadata = json.load(f)

        return metadata

    def update_completed_steps(
        self,
        step_name: str
    ):
        """ 
        Append completed pipeline step.
        """
        metadata = self.load_metadata()

        if step_name not in metadata["completed_steps"]:
            metadata["completed_steps"].append(step_name)

        self.save_metadata(metadata)

    def update_canonical_inputs(
        self,
        new_inputs: Dict
    ):
        """ 
        Update canonical input records.
        """
        metadata = self.load_metadata()
        metadata["canonical_inputs"].update(new_inputs)
        self.save_metadata(metadata)

    def update_raw_inputs(
        self,
        new_inputs:Dict
    ):
        """ 
        Update raw input records.
        """
        metadata = self.load_metadata()
        metadata["raw_inputs"].update(new_inputs)
        self.save_metadata(metadata)

    def metadata_exists(self) -> bool:
        """ 
        check whether metadata.json exists.
        """

        return os.path.exists(self.metadata_path)