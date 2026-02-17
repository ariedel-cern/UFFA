import tempfile
import subprocess
import os
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from .Utils import AnalysisUtils as au

# Global logger
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
                "overwrite": False  # optional
                "mergeJobs": #cores  # optional
                "compression": 505  # optional
            }
        """
        self.inputs = config.get("inputs", None)
        self.output = config.get("output", None)
        self.overwrite = config.get("overwrite", False)
        self.mergeJobs = config.get("mergeJobs", os.cpu_count())
        self.compression = config.get("compression", 505)

    @staticmethod
    def _split(path):
        if ":" not in path:
            raise ValueError(f"Invalid ROOT path (expected file.root:dir): {path}")
        return path.split(":", 1)

    def _file_exists(self, filepath):
        return os.path.exists(filepath)

    def _copy_single_file(self, entry, i, tmpdir, out_dir):
        """Copy a single file using rootcp."""
        src_file, src_dir = self._split(entry)
        if not self._file_exists(src_file):
            raise RuntimeError(f"Input ROOT file does not exist: {src_file}")

        tmp_file = os.path.join(tmpdir, f"input_{i}.root")
        logger.debug(f"Copying {entry} -> {tmp_file}:{out_dir}")
        subprocess.check_call(
            [
                "rootcp",
                "--compress",
                f"{self.compression}",
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
            raise ValueError("Output directory name must not be empty")

        if os.path.exists(out_file) and not self.overwrite:
            logger.warning(
                f"Output file {out_file} exists and overwrite=False. Skipping merge."
            )
            return

        with tempfile.TemporaryDirectory(prefix="femto_merge_") as tmpdir:
            tmp_inputs = [None] * len(self.inputs)

            # Parallel copy using ThreadPoolExecutor
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

            logger.debug(f"Merging {len(tmp_inputs)} files -> {out_file}")
            au.CreateOutputDir(out_file)
            subprocess.check_call(
                [
                    "hadd",
                    f"-f{self.compression}",
                    "-j",
                    f"{self.mergeJobs}",
                    "-f",
                    out_file,
                ]
                + tmp_inputs
            )
            logger.debug("Merge complete")
