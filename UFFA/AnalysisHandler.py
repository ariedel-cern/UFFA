import ROOT as rt

rt.gROOT.SetBatch(True)
rt.TH1.AddDirectory(False)

import itertools
import os
from functools import partial
import logging
from concurrent.futures import ProcessPoolExecutor
import tempfile

from . import CorrelationHandler as ch
from .Utils import AnalysisUtils as au

logger = logging.getLogger(__name__)
if not logger.hasHandlers():
    logging.basicConfig(level=logging.DEBUG)


class AnalysisHandler:
    """
    AnalysisHandler class
    """

    def __init__(self, analysis_dict):
        """
        AnalysisHandler constructor

        Args:
            analysis_dict (dict): Dictionary to configure AnalysisHandler
        """

        # # declare class variables
        self._config_list = []
        self._output_path_list = []

        # get input file
        self._input_file_name = analysis_dict.get("Input_File", None)

        # by default it is assumed that SE and ME are in the same file
        # adding option to specify the file separately
        self._input_file_name_se = analysis_dict.get("Input_File_SE", None)
        self._input_file_name_me = analysis_dict.get("Input_File_ME", None)

        # paths to SE and ME distribution
        self._path_se = analysis_dict.get("Path_SE", "SE")
        self._path_me = analysis_dict.get("Path_ME", "ME")

        self._output_file_name = analysis_dict.get("Output_File", "./output.root")
        logger.debug("Output file name: %s", self._output_file_name)
        au.CreateOutputDir(self._output_file_name)
        self._output_dir_name = analysis_dict.get("Output_Dir", "output")
        logger.debug(
            "Name for TDirectoryFile in output file: %s", self._output_dir_name
        )

        # these are checked and logged later in CorrelationHandler class
        self._normalization_range = analysis_dict.get(
            "Normalization_Range", (0.24, 0.34)
        )

        if (
            not isinstance(self._normalization_range, (tuple, list))
            or len(self._normalization_range) != 2
        ):
            raise ValueError(
                "Normalization_Range must be a tuple/list of two floats (low, high)"
            )

        self._rebin_list = analysis_dict.get("Rebin_Factor_Kstar_Axis", [1])
        self._axis_kstar = analysis_dict.get("Index_Kstar_Axis", 0)
        self._rescale_kstar = analysis_dict.get("Rescale_Factor_Kstar_Axis", -1)
        self._axis_reweight = analysis_dict.get("Index_Reweight_Axis", -1)
        self._range_list = analysis_dict.get("Bins", [])

        # init same and mixed event distribution to none
        self._Se = None
        self._Me = None

    def _GetHistograms(self):
        """
        Get histograms from input file, if defined
        """

        # if same and mixed event are set externally, break out early
        if self._Se is not None and self._Me is not None:
            return

        if self._input_file_name:
            input_file = rt.TFile(self._input_file_name, "READ")

            self._Se = au.GetObjectFromFile(input_file, self._path_se)
            logger.debug(
                "Get Same event distribution from %s at path: %s",
                self._input_file_name,
                self._path_se,
            )

            self._Me = au.GetObjectFromFile(input_file, self._path_me)
            logger.debug(
                "Mixed event distribution from %s at path: %s",
                self._input_file_name,
                self._path_me,
            )
            input_file.Close()
        elif self._input_file_name_se and self._input_file_name_me:
            input_file_se = rt.TFile(self._input_file_name_se, "READ")

            self._Se = au.GetObjectFromFile(input_file_se, self._path_se)
            logger.debug(
                "Get Same event distribution from %s at path: %s",
                self._input_file_name_se,
                self._path_se,
            )
            input_file_se.Close()

            input_file_me = rt.TFile(self._input_file_name_me, "READ")
            self._Me = au.GetObjectFromFile(input_file_me, self._path_me)
            logger.debug(
                "Get Mixed event distribution from %s at path: %s",
                self._input_file_name_me,
                self._path_me,
            )
            input_file_me.Close()
        else:
            raise ValueError(
                "Either Input_File or both Input_File_SE and Input_File_ME must be provided."
            )

    def SetHistograms(self, Se, Me):
        """
        Set histograms directly without retrieving them from a file

        Args:
            Se (THX): Same event distribtuion
            Me (THX): Mixed event distribution
        """
        self._Se = Se.Clone("SE")
        self._Me = Me.Clone("ME")
        logger.debug("Same and mixed event distribution passed directly")

    def _ProcessHandler(self, index, tmp_dir):
        """
        Process CorrelationHandler
        Args:
            index (int): index of configuration

        Returns:
            CorrelationHandler object that has been processed (i.e. correlation function has been computed)
        """
        Handler = ch.CorrelationHandler(
            self._Se,
            self._Me,
            self._normalization_range,
            self._config_list[index]["Rebin"],
            self._rescale_kstar,
            self._axis_kstar,
            self._axis_reweight,
            self._config_list[index]["Ranges"],
        )
        Handler.FinalTouch()

        handler_outputfile_path = os.path.join(
            tmp_dir, f"{index}_" + os.path.basename(self._output_file_name)
        )

        OutputFile = rt.TFile(handler_outputfile_path, "RECREATE")
        output_dir = rt.TDirectoryFile(
            self._output_dir_name, self._output_dir_name, "", OutputFile
        )
        HandlerOutputDir = rt.TDirectoryFile(
            self._config_list[index]["Name"],
            self._config_list[index]["Name"],
            "",
            output_dir,
        )
        Handler.SaveOutput(HandlerOutputDir)
        OutputFile.Write()
        OutputFile.Close()

        return handler_outputfile_path

    def SteerAnalysis(self, parallel=False, workers=0):
        """
        Steer analysis
        Args:
            parallel (bool): Steer analysis in parallel if true
            workers (int): Number of launched processes. If number is less then 0, use all available cores
        """

        self._GetHistograms()

        # generate all configurations, i.e. all rebins and all ranges
        self._GenerateConfigurations()

        with tempfile.TemporaryDirectory(prefix="UFFA_") as tmp_dir:
            logger.debug("Partial outputs will be stored in %s", tmp_dir)
            if parallel:
                if workers <= 0:
                    workers = os.cpu_count()
                with ProcessPoolExecutor(max_workers=workers) as executor:
                    func = partial(self._ProcessHandler, tmp_dir=tmp_dir)
                    self._output_path_list = list(
                        executor.map(func, range(len(self._config_list)))
                    )
            else:
                for i, tmp in [(i, tmp_dir) for i in range(len(self._config_list))]:
                    self._output_path_list.append(self._ProcessHandler(i, tmp))
            self.MergeOutputs()

    def _GenerateConfigurations(self):
        """
        Generate all configurations
        """

        # split all given ranges
        range_list = []
        if self._range_list:
            ranges = []
            # loop over all dimensions
            for r in self._range_list:
                temp = []
                # if only one range is defined, take it
                if len(r) <= 2:
                    temp.append(tuple(r))
                # if we have more than 1, split it into upper and lower limits
                else:
                    for i in range(0, len(r) - 1):
                        temp.append((r[i], r[i + 1]))
                ranges.append(temp)
            # get all combinations
            range_list = list(itertools.product(*ranges))
        else:
            range_list = [[]]

        # now get all combinations of ranges and rebins
        for rebin in self._rebin_list:
            for ranges in range_list:
                d = {
                    "Name": self._GetConfigName(rebin, ranges),
                    "Rebin": rebin,
                    "Ranges": ranges,
                }
                self._config_list.append(d)

    def _GetConfigName(self, rebin, ranges):
        """
        Generate a name for a configuration
        Name will be given to the TDirectoryFile used to store the output of configuration with given index
        Args:
            rebin (int): Rebin factor
            ranges (list): list of ranges
        """
        # rebin 1 is a special case
        if rebin == 1:
            ConfigName = "NoRebin"
        else:
            ConfigName = f"Rebin_{rebin}"

        if ranges is not None:
            for i, e in enumerate(ranges):
                # only add parts to a name if there is a range defined for a dimension
                if e:
                    # set name to the edges of the bin
                    ConfigName = ConfigName + f"_Dim_{i}-{str(e[0])}-{str(e[1])}"
        return ConfigName

    def MergeOutputs(self):
        """
        Merge output of all CorrelationHandler
        """

        logger.debug("Merging final output into %s", self._output_file_name)

        merger = rt.TFileMerger()
        merger.OutputFile(self._output_file_name, "RECREATE")
        for file in self._output_path_list:
            merger.AddFile(file)

        merger.Merge()
