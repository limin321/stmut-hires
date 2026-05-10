"""
OutputCleaner — removes intermediate output folders, keeping only `figures` and `tables`.
 
Integration:
    from data_io.output_cleaner import OutputCleaner
 
    cleaner = OutputCleaner(output_dir)
    cleaner.clean()                      # always clean
    cleaner.clean(enabled=args.clean)    # respect CLI flag from cli_parser.py
"""
 
import shutil
from pathlib import Path
 
 
class OutputCleaner:
    """
    Cleans an output directory by removing all subdirectories except
    a configurable keep-list (default: figures, tables).
 
    Parameters
    ----------
    output_dir : str | Path
        Root output directory to clean.
    keep : list[str], optional
        Folder names to preserve. Defaults to ["figures", "tables"].
    verbose : bool
        Print what is removed / kept.
    """
 
    DEFAULT_KEEP = {"figures", "tables"}
 
    def __init__(
        self,
        output_dir: str | Path,
        keep: list[str] | None = None,
        verbose: bool = True,
    ):
        self.output_dir = Path(output_dir)
        self.keep = set(keep) if keep is not None else self.DEFAULT_KEEP
        self.verbose = verbose
 
    def clean(self, enabled: bool = True) -> None:
        """
        Run the cleaning step.
 
        Parameters
        ----------
        enabled : bool
            When False, returns immediately without touching anything.
            Pass args.clean directly from cli_parser.py.
        """
        if not enabled:
            if self.verbose:
                print("[OutputCleaner] Cleaning skipped.")
            return
 
        if not self.output_dir.exists():
            raise FileNotFoundError(f"Output directory not found: {self.output_dir}")
 
        removed, kept, missing = [], [], []
 
        for item in sorted(self.output_dir.iterdir()):
            if not item.is_dir():
                continue  # leave loose files alone
            if item.name in self.keep:
                kept.append(item.name)
            else:
                shutil.rmtree(item)
                removed.append(item.name)
 
        for name in self.keep:
            if name not in kept:
                missing.append(name)
 
        if self.verbose:
            self._report(removed, kept, missing)
 
    def _report(self, removed: list, kept: list, missing: list) -> None:
        print(f"[OutputCleaner] Directory : {self.output_dir}")
        if kept:
            print(f"  Kept    : {', '.join(sorted(kept))}")
        if removed:
            print(f"  Removed : {', '.join(sorted(removed))}")
        if missing:
            print(f"  Missing : {', '.join(sorted(missing))} (not present, nothing to keep)")

            