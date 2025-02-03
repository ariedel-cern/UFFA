import ROOT as rt

# rt.EnableThreadSafety()

import itertools
import os
import logging
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor

from . import CorrelationHandler as ch
from .Utils import AnalysisUtils as au
from .Utils import CorrelationUtils as cu

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
        self._axis_kstar = None
        self._rescale_kstar = -1
        self._axis_reweight = None

        self._rebin_list = []
        self._range_list = []
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

        # check name of output TDirectoryFile
        au.CheckDictEntry(analysis_dict, "Output_Dir", str)
        self._output_dir_name = analysis_dict["Output_Dir"]
        logger.debug(
            "Name for TDirectoryFile in output file: %s", self._output_dir_name
        )

        # check same event distribution
        au.CheckDictEntry(analysis_dict, "Path_SE", str)
        self._path_se = analysis_dict["Path_SE"]
        self._Se = au.GetObjectFromFile(self._input_file, self._path_se)
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
        self._Me = au.GetObjectFromFile(self._input_file, self._path_me)
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

        # au.CheckDictEntry(analysis_dict, "Rebin_Factor_Kstar_Axis", int)
        self._rebin_list = analysis_dict["Rebin_Factor_Kstar_Axis"]
        if not self._rebin_list:
            self._rebin_list = [1]

        au.CheckDictEntry(analysis_dict, "Index_Kstar_Axis", int)
        self._axis_kstar = analysis_dict["Index_Kstar_Axis"]

        au.CheckDictEntry(analysis_dict, "Rescale_Factor_Kstar_Axis", int)
        self._rescale_kstar = analysis_dict["Rescale_Factor_Kstar_Axis"]

        au.CheckDictEntry(analysis_dict, "Index_Reweight_Axis", int)
        self._axis_reweight = analysis_dict["Index_Reweight_Axis"]

        au.CheckDictEntry(analysis_dict, "Bins", list)
        self._range_list = analysis_dict["Bins"]

    def ProcessHandler(self, index):
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

        return Handler

    def SteerAnalysis(self, parallel=False, workers=0):
        """
        Steer analysis
        Args:
            parallel (bool): Steer analysis in parallel if true
            workers (int): Number of launched processes. If number is less then 0, use all avaiable cores
        """

        # generate all configurations, i.e. all rebins and all ranges
        self.GenerateConfigurations()

        if parallel == True:
            if workers <= 0:
                workers = os.cpu_count()
            with ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
                self._handler_list = list(
                    executor.map(
                        self.ProcessHandler, list(range(len(self._config_list)))
                    )
                )
        else:
            for i in range(len(self._config_list)):
                self._handler_list.append(self.ProcessHandler(i))

        self.SaveToOutput(parallel, workers)

        self._input_file.Close()

    def GenerateConfigurations(self):
        """
        Generate all configurations
        """

        # split all given ranges
        range_list = []
        if self._range_list is not None:
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
                    "Name": self.GetConfigName(rebin, ranges),
                    "Rebin": rebin,
                    "Ranges": ranges,
                }
                self._config_list.append(d)

    def GetConfigName(self, rebin, ranges):
        """
        Generate a name for a configuration
        Name will be given to the TDirectoryFile used to store the output of configuration with given index
        Args:
            rebin (int): Rebin factor
            range (list): list of ranges
        """
        ConfigName = f"Rebin_{rebin}"

        if ranges is not None:
            for i, e in enumerate(ranges):
                # only add parts to a name if there is a range defined for a dimension
                if e:
                    # set name to the edges of the bin
                    ConfigName = ConfigName + f"_Dim_{i}-{str(e[0])}-{str(e[1])}"
        return ConfigName

    def SaveHandler(self, index):
        """
        Save CorrelationHandler

        Args:
            dict (dictionary): Save output of configuration for given index pair. First element is the range and second index is the rebin
        """

        HandlerOutputDir = rt.TDirectoryFile(
            self._config_list[index]["Name"],
            self._config_list[index]["Name"],
            "",
            self._OutputDir,
        )
        self._handler_list[index].SaveOutput(HandlerOutputDir)

    def SaveToOutput(self, parallel=False, workers=-1):
        """
        Save analysis output to file
        Args:
            parallel (bool): Save analysis output in parallel if true
            workers (int): Number of launched threads(!). If number is less then 0, use all avaiable cores.
        """
        OutputFile = rt.TFile(self._output_file_name, "UPDATE")
        # check first if a tdirectoryfile with the same name already exists
        temp_dir = OutputFile.Get(self._output_dir_name)
        if temp_dir is not None:
            OutputFile.Delete(f"{self._output_dir_name};*")

        self._OutputDir = rt.TDirectoryFile(
            self._output_dir_name, self._output_dir_name, "", OutputFile
        )

        if parallel == True:
            if workers <= 0:
                workers = os.cpu_count()
            with ThreadPoolExecutor(max_workers=workers) as executor:
                executor.map(self.SaveHandler, list(range(len(self._handler_list))))
        else:
            for i in range(len(self._handler_list)):
                self.SaveHandler(i)

        self._OutputDir.Write(
            self._output_dir_name, rt.TObject.kSingleKey + rt.TObject.kWriteDelete
        )

        OutputFile.Save()
        OutputFile.Close()
