import tempfile
import subprocess
import os
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from .Utils import AnalysisUtils as au

logger = logging.getLogger(__name__)
if not logger.hasHandlers():
    logging.basicConfig(level=logging.DEBUG)


class FemtoMerger:
    """
    Merge femto outputs using rootcp + hadd.
    """

    def __init__(self, config):
        """
        config: dict
            {
                "inputs": ["file.root:SE", "file.root:ME"],
                "output": "merged.root:femto",
                "overwrite": False,      # optional
                "mergeJobs": 4,          # optional, defaults to os.cpu_count()
                "compression": 505       # optional
            }
        """
        self.inputs = config.get("inputs", None)
        self.output = config.get("output", None)

        if not self.inputs:
            raise ValueError("'inputs' must be provided and non-empty.")
        if not self.output:
            raise ValueError("'output' must be provided.")

        self.overwrite = config.get("overwrite", False)
        self.mergeJobs = config.get("mergeJobs", os.cpu_count())
        self.compression = config.get("compression", 505)

    @staticmethod
    def _split(path):
        if ":" not in path:
            raise ValueError(f"Invalid ROOT path (expected file.root:dir): {path}")
        return path.split(":", 1)

    def _copy_single_file(self, entry, i, tmpdir, out_dir):
        """Copy a single file using rootcp."""
        src_file, src_dir = self._split(entry)
        if not os.path.exists(src_file):
            raise RuntimeError(f"Input ROOT file does not exist: {src_file}")
        tmp_file = os.path.join(tmpdir, f"input_{i}.root")
        logger.debug("Copying %s -> %s:%s", entry, tmp_file, out_dir)
        subprocess.check_call(
            [
                "rootcp",
                "--compress",
                str(self.compression),
                "--recursive",
                "--recreate",
                f"{src_file}:{src_dir}",
                f"{tmp_file}:{out_dir}",
            ]
        )
        return tmp_file

    def merge(self):
        out_file, out_dir = self._split(self.output)
        if not out_dir:
            raise ValueError("Output directory name must not be empty.")
        if os.path.exists(out_file) and not self.overwrite:
            logger.warning(
                "Output file %s exists and overwrite=False. Skipping merge.", out_file
            )
            return

        with tempfile.TemporaryDirectory(prefix="femto_merge_") as tmpdir:
            tmp_inputs = [None] * len(self.inputs)

            with ThreadPoolExecutor(max_workers=self.mergeJobs) as executor:
                futures = {
                    executor.submit(
                        self._copy_single_file, entry, i, tmpdir, out_dir
                    ): i
                    for i, entry in enumerate(self.inputs)
                }
                for future in as_completed(futures):
                    i = futures[future]
                    tmp_inputs[i] = future.result()

            logger.debug("Merging %d files -> %s", len(tmp_inputs), out_file)
            au.CreateOutputDir(out_file)
            subprocess.check_call(
                [
                    "hadd",
                    f"-f{self.compression}",
                    "-j",
                    str(self.mergeJobs),
                    out_file,
                ]
                + tmp_inputs
            )
            logger.debug("Merge complete")
