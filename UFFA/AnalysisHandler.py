import ROOT as rt
import itertools
import os
import logging
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor

from . import CorrelationHandler as ch
from .Utils import AnalysisUtils as au

logger = logging.getLogger(__name__)


class AnalysisHandler:
    """
    AnalysisHandler class
    """

    def __init__(self, analysis_dict):
        """
        AnalysisHandler constructor

        Args:
            AnalysisDict (dict): Dictionary to configure AnalysisHandler
        """

        # declare class variables
        self._input_file_name = None
        self._input_file = None
        self._output_file_name = None
        self._output_file = None
        self._output_dir_name = None
        self._path_se = None
        self._path_me = None
        self._normalization_range = None
        self._rebin_kstar = None
        self._axis_kstar = None
        self._axis_reweight = None
        self._ranges = None

        self._number_configs = 1
        self._config_list = []
        self._handler_list = []

        # Check input file
        au.CheckDictEntry(analysis_dict, "Input_File", str)
        self._input_file_name = analysis_dict["Input_File"]

        if not rt.gSystem.AccessPathName(self._input_file_name):
            self._input_file = rt.TFile(self._input_file_name, "READ")
            logger.debug("Opened input file: %s", self._input_file_name)
        else:
            raise ValueError(f"{self._input_file_name} does not exist")

        au.CheckDictEntry(analysis_dict, "Output_File", str)
        self._output_file_name = analysis_dict["Output_File"]
        # logger.debug("(Re)created output file: %s", self._OutputFileName)

        # check name of output TDirectoryFile
        au.CheckDictEntry(analysis_dict, "Output_Dir", str)
        self._output_dir_name = analysis_dict["Output_Dir"]
        logger.debug("Name for TDirectoryFile in output file: %s", self._output_dir_name)

        # check same event distribution
        au.CheckDictEntry(analysis_dict, "Path_SE", str)
        self._path_se = analysis_dict["Path_SE"]
        self._Se = self._input_file.Get(self._path_se)
        if self._Se == None:
            raise ValueError(
                f"Same event distribution not found. Is '{self._path_se}' the correct path?"
            )
        logger.debug(
            "Same event distribution retrieved from input file: %s", self._path_se
        )

        # check mixed event distribution
        au.CheckDictEntry(analysis_dict, "Path_ME", str)
        self._path_me = analysis_dict["Path_ME"]
        self._Me = self._input_file.Get(self._path_me)
        if self._Me == None:
            raise ValueError(
                f"Mixed event distribution not found. Is '{self._path_me}' the correct path?"
            )
        logger.debug(
            "Mixed event distribution retrieved from input file: %s", self._path_me
        )

        # these are checked and logged later in CorrelationHandler class
        au.CheckDictEntry(analysis_dict, "Normalization_Range", tuple)
        self._normalization_range = analysis_dict["Normalization_Range"]

        au.CheckDictEntry(analysis_dict, "Rebin_Factor_Kstar_Axis", int)
        self._rebin_kstar = analysis_dict["Rebin_Factor_Kstar_Axis"]

        au.CheckDictEntry(analysis_dict, "Index_Kstar_Axis", int)
        self._axis_kstar = analysis_dict["Index_Kstar_Axis"]

        au.CheckDictEntry(analysis_dict, "Index_Reweight_Axis", int)
        self._axis_reweight = analysis_dict["Index_Reweight_Axis"]

        au.CheckDictEntry(analysis_dict, "Bins", list)
        self._ranges = analysis_dict["Bins"]

    def ProcessHandler(self, Index):
        """
        Process CorrelationHandler
        Args:
            Index (int): Index of configuration

        Returns:
            CorrelationHandler object that has been processed (i.e. correlation function has been computed)
        """
        Handler = ch.CorrelationHandler(
            self._Se,
            self._Me,
            self._normalization_range,
            self._rebin_kstar,
            self._axis_kstar,
            self._axis_reweight,
            self._config_list[Index],
        )
        Handler.FinalTouch()
        return Handler

    def SteerAnalysis(self, parallel=False, workers=-1):
        """
        Steer analysis
        Args:
            parallel (bool): Steer analysis in parallel if true
            workers (int): Number of launched processes. If number is less then 0, use all avaiable cores
        """

        # split bins to get all configurations
        self.SplitBins()

        if parallel == True:
            if workers <= 0:
                workers = os.cpu_count()
            with ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
                self._handler_list = list(
                    executor.map(
                        self.ProcessHandler, [i for i in range(self._number_configs)]
                    )
                )
        else:
            for i in range(self._number_configs):
                self._handler_list.append(self.ProcessHandler(i))

        self.SaveToOutput(parallel, workers)

    def SplitBins(self):
        """
        Split list of bins and generate a list of all configurations
        """

        def split_into_pairs(lst):
            """(Helper function) Split a list into consecutive pairs"""
            return [lst[i : i + 2] for i in range(len(lst) - 1)]

        processed_sublists = []
        for sublist in self._ranges:
            # Process bins for each dimension
            # If more than 2 edges are defined, there is more than 1 bin
            if len(sublist) > 2:
                # Split list into pairs
                processed_sublists.append(split_into_pairs(sublist))
            else:
                # If there is only one bin keep the sublist as it is
                processed_sublists.append([sublist])

        # Generate all combinations using Cartesian product
        combinations = list(itertools.product(*processed_sublists))
        # CorrelationHandler expects list of tuples, so transfrom lists to tuples
        self._config_list = [
            [tuple(inner_list) for inner_list in sublist] for sublist in combinations
        ]
        self._number_configs = len(self._config_list)

    def GetConfigName(self, index):
        """
        Generate a name for a configuration
        Name will be given to the TDirectoryFile used to store the output of configuration with given index
        Args:
            index (int): Index of the configuration
        """
        # get configuration at passed index
        Config = self._config_list[index]
        ConfigName = ""
        for i, e in enumerate(Config):
            # only add parts to a name if there is a bin defined for a dimension
            if e:
                # set name to the edged of the bin
                ConfigName = ConfigName + f"{str(e[0])}-{str(e[1])}"
                # add underscore if there follows another bin
                if i < len(Config) - 1:
                    ConfigName = ConfigName + "_"
        # catch corner case where there are no cuts applied
        if ConfigName == "":
            ConfigName = "Analysis"
        return ConfigName

    def SaveHandler(self, index):
        """
        Save CorrelationHandler

        Args:
            index (int): Save output of configuration with given index to TDirectoryFile
        """
        HandlerOutputDir = rt.TDirectoryFile(
            self.GetConfigName(index), self.GetConfigName(index), "", self._OutputDir
        )
        self._handler_list[index].SaveOutput(HandlerOutputDir)

    def SaveToOutput(self, parallel=False, workers=-1):
        """
        Save analysis output to file
        Args:
            parallel (bool): Save analysis output in parallel if true
            workers (int): Number of launched threads(!). If number is less then 0, use all avaiable cores.
        """
        OutputFile = rt.TFile(self._output_file_name, "RECREATE")
        self._OutputDir = rt.TDirectoryFile(
            self._output_dir_name, self._output_dir_name, "", OutputFile
        )

        if parallel == True:
            if workers <= 0:
                workers = os.cpu_count()
            with ThreadPoolExecutor(max_workers=workers) as executor:
                executor.map(
                    self.SaveHandler, [index for index in range(self._number_configs)]
                )
        else:
            for index in range(self._number_configs):
                HandlerOutputDir = rt.TDirectoryFile(
                    self.GetConfigName(index),
                    self.GetConfigName(index),
                    "",
                    self._OutputDir,
                )
                self._handler_list[index].SaveOutput(HandlerOutputDir)

        self._OutputDir.Write(self._output_dir_name, rt.TObject.kSingleKey)

        OutputFile.Save()
        OutputFile.Close()
