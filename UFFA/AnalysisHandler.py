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

    def __init__(self, AnalysisDict):
        """
        AnalysisHandler constructor

        Args:
            AnalysisDict (dict): Dictionary to configure AnalysisHandler
        """

        # declare class variables
        self._InputFileName = None
        self._InputFile = None
        self._OutputFileName = None
        self._OutputFile = None
        self._OutputDirName = None
        self._PathSe = None
        self._PathMe = None
        self._NormRange = None
        self._RebinFaktorKstarAxis = None
        self._IndexKstarAxis = None
        self._IndexReweightAxis = None
        self._BinList = None

        self._NConfigs = 1
        self._ConfigList = []
        self._HandlerList = []

        # Check input file
        au.CheckDictEntry(AnalysisDict, "Input_File", str)
        self._InputFileName = AnalysisDict["Input_File"]

        if not rt.gSystem.AccessPathName(self._InputFileName):
            self._InputFile = rt.TFile(self._InputFileName, "READ")
            logger.debug("Opened input file: %s", self._InputFileName)
        else:
            raise ValueError(f"{self._InputFileName} does not exist")

        au.CheckDictEntry(AnalysisDict, "Output_File", str)
        self._OutputFileName = AnalysisDict["Output_File"]
        # logger.debug("(Re)created output file: %s", self._OutputFileName)

        # check name of output TDirectoryFile
        au.CheckDictEntry(AnalysisDict, "Output_Dir", str)
        self._OutputDirName = AnalysisDict["Output_Dir"]
        logger.debug("Name for TDirectoryFile in output file: %s", self._OutputDirName)

        # check same event distribution
        au.CheckDictEntry(AnalysisDict, "Path_SE", str)
        self._PathSe = AnalysisDict["Path_SE"]
        self._Se = self._InputFile.Get(self._PathSe)
        if self._Se == None:
            raise ValueError(
                f"Same event distribution not found. Is '{self._PathSe}' the correct path?"
            )
        logger.debug(
            "Same event distribution retrieved from input file: %s", self._PathSe
        )

        # check mixed event distribution
        au.CheckDictEntry(AnalysisDict, "Path_ME", str)
        self._PathMe = AnalysisDict["Path_ME"]
        self._Me = self._InputFile.Get(self._PathMe)
        if self._Me == None:
            raise ValueError(
                f"Mixed event distribution not found. Is '{self._PathMe}' the correct path?"
            )
        logger.debug(
            "Mixed event distribution retrieved from input file: %s", self._PathMe
        )

        # these are checked and logged later in CorrelationHandler class
        au.CheckDictEntry(AnalysisDict, "Normalization_Range", tuple)
        self._NormRange = AnalysisDict["Normalization_Range"]

        au.CheckDictEntry(AnalysisDict, "Rebin_Factor_Kstar_Axis", int)
        self._RebinFaktorKstarAxis = AnalysisDict["Rebin_Factor_Kstar_Axis"]

        au.CheckDictEntry(AnalysisDict, "Index_Kstar_Axis", int)
        self._IndexKstarAxis = AnalysisDict["Index_Kstar_Axis"]

        au.CheckDictEntry(AnalysisDict, "Index_Reweight_Axis", int)
        self._IndexReweightAxis = AnalysisDict["Index_Reweight_Axis"]

        au.CheckDictEntry(AnalysisDict, "Bins", list)
        self._BinList = AnalysisDict["Bins"]

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
            self._NormRange,
            self._RebinFaktorKstarAxis,
            self._IndexKstarAxis,
            self._IndexReweightAxis,
            self._ConfigList[Index],
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
                self._HandlerList = list(
                    executor.map(
                        self.ProcessHandler, [i for i in range(self._NConfigs)]
                    )
                )
        else:
            for i in range(self._NConfigs):
                self._HandlerList.append(self.ProcessHandler(i))

        self.SaveToOutput(parallel, workers)

    def SplitBins(self):
        """
        Split list of bins and generate a list of all configurations
        """

        def split_into_pairs(lst):
            """(Helper function) Split a list into consecutive pairs"""
            return [lst[i : i + 2] for i in range(len(lst) - 1)]

        processed_sublists = []
        for sublist in self._BinList:
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
        self._ConfigList = [
            [tuple(inner_list) for inner_list in sublist] for sublist in combinations
        ]
        self._NConfigs = len(self._ConfigList)

    def GetConfigName(self, index):
        """
        Generate a name for a configuration
        Name will be given to the TDirectoryFile used to store the output of configuration with given index
        Args:
            index (int): Index of the configuration
        """
        # get configuration at passed index
        Config = self._ConfigList[index]
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
        self._HandlerList[index].SaveOutput(HandlerOutputDir)

    def SaveToOutput(self, parallel=False, workers=-1):
        """
        Save analysis output to file
        Args:
            parallel (bool): Save analysis output in parallel if true
            workers (int): Number of launched threads(!). If number is less then 0, use all avaiable cores.
        """
        OutputFile = rt.TFile(self._OutputFileName, "RECREATE")
        self._OutputDir = rt.TDirectoryFile( self._OutputDirName, self._OutputDirName, "", OutputFile)

        if parallel == True:
            if workers <= 0:
                workers = os.cpu_count()
            with ThreadPoolExecutor(max_workers=workers) as executor:
                executor.map(
                    self.SaveHandler, [index for index in range(self._NConfigs)]
                )
        else:
            for index in range(self._NConfigs):
                HandlerOutputDir = rt.TDirectoryFile(
                    self.GetConfigName(index),
                    self.GetConfigName(index),
                    "",
                    self._OutputDir,
                )
                self._HandlerList[index].SaveOutput(HandlerOutputDir)

        self._OutputDir.Write(self._OutputDirName, rt.TObject.kSingleKey)

        OutputFile.Save()
        OutputFile.Close()
