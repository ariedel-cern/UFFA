import ROOT as rt

rt.gROOT.SetBatch(True)
rt.EnableThreadSafety()
rt.TH1.AddDirectory(False)

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
        self.__input_file_name = None
        self.__output_file_name = None
        self.__output_file = None
        self.__output_dir_name = None
        self.__path_se = None
        self.__path_me = None
        self.__normalization_range = None
        self.__axis_kstar = None
        self.__rescale_kstar = -1
        self.__axis_reweight = None

        self.__rebin_list = []
        self.__range_list = []
        self.__config_list = []
        self.__handler_list = []

        # Check input file
        self.__input_file_name = analysis_dict.get("Input_File", None)
        logger.debug("Input file: name %s", self.__input_file_name)

        # check same event distribution
        self.__path_se = analysis_dict.get("Path_SE", "SE")
        self.__path_me = analysis_dict.get("Path_ME", "ME")

        self.__output_file_name = analysis_dict.get("Output_File", "./output.root")
        logger.debug("Output file name: %s", self.__output_file_name)
        au.CreateOutputDir(self.__output_file_name)
        self.__output_dir_name = analysis_dict.get("Output_Dir", "analysis")
        logger.debug(
            "Name for TDirectoryFile in output file: %s", self.__output_dir_name
        )

        # these are checked and logged later in CorrelationHandler class
        self.__normalization_range = analysis_dict.get(
            "Normalization_Range", (0.24, 0.34)
        )

        self.__rebin_list = analysis_dict.get("Rebin_Factor_Kstar_Axis", [1])
        self.__axis_kstar = analysis_dict.get("Index_Kstar_Axis", 0)
        self.__rescale_kstar = analysis_dict.get("Rescale_Factor_Kstar_Axis", -1)
        self.__axis_reweight = analysis_dict.get("Index_Reweight_Axis", -1)
        self.__range_list = analysis_dict.get("Bins", [])

    def GetHistograms(self):
        """
        Get histograms from input file, if defined
        """
        if self.__input_file_name == None:
            return

        input_file = rt.TFile(self.__input_file_name, "READ")

        self._Se = au.GetObjectFromFile(input_file, self.__path_se)
        logger.debug("Same event distribution at path: %s", self.__path_se)

        self._Me = au.GetObjectFromFile(input_file, self.__path_me)
        logger.debug("Mixed event distribution at path: %s", self.__path_me)

    def SetHistograms(self, Se, Me):
        """
        Set histograms directly without retrieving them from a file

        Args:
            Se (THX): Same event distribtuion
            Me (THX): Mixed event distribution
        """
        self._Se = Se.Clone("SE")
        self._Me = Se.Clone("ME")
        logger.debug("Same and mixed event distribution passed directly")

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
            self.__normalization_range,
            self.__config_list[index]["Rebin"],
            self.__rescale_kstar,
            self.__axis_kstar,
            self.__axis_reweight,
            self.__config_list[index]["Ranges"],
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

        self.GetHistograms()

        # generate all configurations, i.e. all rebins and all ranges
        self.GenerateConfigurations()

        if parallel == True:
            if workers <= 0:
                workers = os.cpu_count()
            with ProcessPoolExecutor(max_workers=workers) as executor:
                self.__handler_list = list(
                    executor.map(
                        self.ProcessHandler, list(range(len(self.__config_list)))
                    )
                )
        else:
            for i in range(len(self.__config_list)):
                self.__handler_list.append(self.ProcessHandler(i))

        self.SaveToOutput(parallel, workers)

    def GenerateConfigurations(self):
        """
        Generate all configurations
        """

        # split all given ranges
        range_list = []
        if self.__range_list is not None:
            ranges = []
            # loop over all dimensions
            for r in self.__range_list:
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
        for rebin in self.__rebin_list:
            for ranges in range_list:
                d = {
                    "Name": self.GetConfigName(rebin, ranges),
                    "Rebin": rebin,
                    "Ranges": ranges,
                }
                self.__config_list.append(d)

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
            self.__config_list[index]["Name"],
            self.__config_list[index]["Name"],
            "",
            self._OutputDir,
        )
        self.__handler_list[index].SaveOutput(HandlerOutputDir)

    def SaveToOutput(self, parallel=False, workers=-1):
        """
        Save analysis output to file
        Args:
            parallel (bool): Save analysis output in parallel if true
            workers (int): Number of launched threads(!). If number is less then 0, use all avaiable cores.
        """
        OutputFile = rt.TFile(self.__output_file_name, "UPDATE")
        # check first if a tdirectoryfile with the same name already exists
        temp_dir = OutputFile.Get(self.__output_dir_name)
        if temp_dir is not None:
            OutputFile.Delete(f"{self.__output_dir_name};*")

        self._OutputDir = rt.TDirectoryFile(
            self.__output_dir_name, self.__output_dir_name, "", OutputFile
        )

        if parallel == True:
            if workers <= 0:
                workers = os.cpu_count()
            with ThreadPoolExecutor(max_workers=workers) as executor:
                executor.map(self.SaveHandler, list(range(len(self.__handler_list))))
        else:
            for i in range(len(self.__handler_list)):
                self.SaveHandler(i)

        self._OutputDir.Write(
            self.__output_dir_name, rt.TObject.kSingleKey + rt.TObject.kWriteDelete
        )

        OutputFile.Save()
        OutputFile.Close()
