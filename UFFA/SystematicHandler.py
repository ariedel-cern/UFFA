import ROOT as rt

rt.gROOT.SetBatch(True)
rt.EnableThreadSafety()
rt.TH1.AddDirectory(False)

import numpy as np
import logging
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor

from . import CorrelationHandler as ch
from .Utils import AnalysisUtils as au

logger = logging.getLogger(__name__)


class SystematicHandler:
    """
    SystematicHandler class
    """

    def __init__(self, systematic_dict):
        """
        SystematicHandler constructor

        Args:
            AnalysisDict (dict): Dictionary to configure AnalysisHandler
        """

        # declare class variables
        self.__output_file_name = systematic_dict.get("Output_File", "systematics.root")
        self.__output_dir_name = systematic_dict.get("Output_Dir", "Systematics")

        self.__default_dict = systematic_dict.get("Default", {})
        self.__variation_dicts = systematic_dict.get("Variations", [])
        self.__stats_dict = systematic_dict.get("Stats", {})

        self.__Binning = systematic_dict.get("Binning")
        self.__signal_region = systematic_dict.get("SignalRegion", (0, 200))
        self.__max_signal_variation = systematic_dict.get("MaxSignalVariation", 0.1)

    def __GetObjects(self):
        if "Histogram" in self.__default_dict:
            self.__default = self.__default_dict.get("Object").Clone("Default")
            logger.debug("Get default object %s", self.__default.GetName())
        else:
            InputFile = rt.TFile(self.__default_dict.get("File"), "READ")
            self.__default = au.GetObjectFromFile(
                InputFile, self.__default_dict.get("Path")
            ).Clone("Default")
            InputFile.Close()
            logger.debug(
                "Get default object %s from file %s at path %s",
                self.__default.GetName(),
                self.__default_dict.get("File"),
                self.__default_dict.get("Path"),
            )

        Signal = self.__default.Integral(
            self.__signal_region[0], self.__signal_region[1]
        )
        logger.debug(
            "Integral %.2f in signal region (%.2f,%.2f)",
            Signal,
            self.__signal_region[0],
            self.__signal_region[1],
        )

        self.__variations = []

        self.__systematics_signal = rt.TH1F(
            "SignalVariation",
            "SignalVariation",
            len(self.__variation_dicts),
            0.5,
            len(self.__variation_dicts) + 0.5,
        )

        for i, variation_dict in enumerate(self.__variation_dicts):
            if "Histogram" in variation_dict:
                variation = variation_dict.get("Object").Clone(f"Variation_{i:02d}")
                logger.debug("Get Variation %d object %s", i, self.__default.GetName())
            else:
                InputFile = rt.TFile(variation_dict.get("File"), "READ")
                variation = au.GetObjectFromFile(
                    InputFile,
                    variation_dict.get("Path"),
                ).Clone(f"Variation_{i:02d}")
                InputFile.Close()
                logger.debug(
                    "Get variation %d object %s from file %s at path %s",
                    i,
                    variation.GetName(),
                    variation_dict.get("File"),
                    variation_dict.get("Path"),
                )

            Signal_Variation = variation.Integral(
                self.__signal_region[0], self.__signal_region[1]
            )
            diff = (Signal - Signal_Variation) / Signal
            logger.debug(
                "Signal %.2f with relative deviation from the signal %.2f",
                Signal_Variation,
                diff,
            )

            self.__systematics_signal.SetBinContent(i + 1, diff)

            if np.abs(diff) <= self.__max_signal_variation:
                self.__variations.append(variation)
                logger.debug("Accept variation %d", i)
            else:
                logger.debug("Reject variation %d", i)
            InputFile.Close()

    def __ComputeSystemtatics(self):

        self.__systematics_hist = self.__default.Clone("SystematicsHist")
        self.__systematics_hist.Reset()
        self.__systematics_dist = rt.TH2D(
            "SystematicDist",
            "SystematicDist",
            self.__systematics_hist.GetNbinsX(),
            self.__systematics_hist.GetXaxis().GetXmin(),
            self.__systematics_hist.GetXaxis().GetXmax(),
            self.__Binning[0],
            self.__Binning[1],
            self.__Binning[2],
        )

        for bin in range(1, self.__default.GetNbinsX()):
            var_min = 0
            var_max = 0
            for variation in self.__variations:
                diff = self.__default.GetBinContent(bin) - variation.GetBinContent(bin)
                center = self.__default.GetBinCenter(bin)
                self.__systematics_dist.Fill(center, diff)

                if diff < var_min:
                    var_min = diff
                if diff > var_max:
                    var_max = diff
            sys = np.abs(var_max - var_min) / np.sqrt(12.0)
            self.__systematics_hist.SetBinContent(bin, sys)

    def __GenerateObjects(self):
        if "Graph" in self.__stats_dict:
            self.__stats_graph = self.__stats_dict.get("Graph").Clone("Graph_Stat")
            logger.debug(
                "Get graph with stat. uncertainty %s",
                self.__default.GetName(),
            )
        else:
            InputFile = rt.TFile(self.__stats_dict.get("File"), "READ")
            self.__stats_graph = au.GetObjectFromFile(
                InputFile, self.__stats_dict.get("Path")
            ).Clone("Graph_Stat")
            InputFile.Close()
            logger.debug(
                "Get graph with stat. uncertainty %s from file %s at path %s",
                self.__stats_graph.GetName(),
                self.__stats_dict.get("File"),
                self.__stats_dict.get("Path"),
            )

            self.__systs_graph = self.__stats_graph.Clone("Graph_Syst")
            self.__bins_graph = self.__stats_graph.Clone("Graph_BinWidth")

            for bin in range(self.__stats_graph.GetN()):
                self.__systs_graph.SetPointError(
                    bin,
                    self.__stats_graph.GetErrorX(bin),
                    self.__systematics_hist.GetBinContent(bin + 1),
                )

            for bin in range(self.__stats_graph.GetN()):
                self.__bins_graph.SetPoint(
                    bin,
                    self.__default.GetBinCenter(bin + 1),
                    self.__default.GetBinContent(bin + 1),
                )
                self.__bins_graph.SetPointError(
                    bin, self.__default.GetBinWidth(bin) / 2, 0
                )

            self.__systs_hist = self.__default.Clone("Hist_Syst")
            self.__stats_hist = self.__default.Clone("Hist_Stat")

            for bin in range(self.__stats_graph.GetN()):
                self.__systs_hist.SetBinError(
                    bin + 1, self.__systematics_hist.GetBinContent(bin + 1)
                )
                self.__stats_hist.SetBinError(
                    bin + 1, self.__stats_graph.GetErrorY(bin)
                )

    def SteerSystematics(self):
        if "Histogram" in self.__default_dict:
            self.__default = self.__default_dict.get("Object").Clone()
            logger.debug("Get default object %s", self.__default.GetName())
        else:
            InputFile = rt.TFile(self.__default_dict.get("File"), "READ")
            self.__default = au.GetObjectFromFile(
                InputFile, self.__default_dict.get("Path")
            ).Clone()
            InputFile.Close()
            logger.debug(
                "Get default object %s from file %s at path %s",
                self.__default.GetName(),
                self.__default_dict.get("File"),
                self.__default_dict.get("Path"),
            )
        self.__GetObjects()

        self.__ComputeSystemtatics()

        self.__GenerateObjects()

        OutputFile = rt.TFile(self.__output_file_name, "RECREATE")

        OutputDir = rt.TDirectoryFile(
            self.__output_dir_name, self.__output_dir_name, "", OutputFile
        )

        systematics_dir = rt.TDirectoryFile("Systematics", "Systematics", "", OutputDir)
        systematics_dir.Add(self.__systematics_hist)
        systematics_dir.Add(self.__systematics_dist)
        systematics_dir.Add(self.__systematics_signal)
        systematics_dir.Add(self.__stats_graph)
        systematics_dir.Add(self.__systs_graph)
        systematics_dir.Add(self.__bins_graph)
        systematics_dir.Add(self.__systs_hist)
        systematics_dir.Add(self.__stats_hist)
        systematics_dir.Write()

        variations_dir = rt.TDirectoryFile("Variations", "Variations", "", OutputDir)
        variations_dir.Add(self.__default)
        for variation in self.__variations:
            variations_dir.Add(variation)

        variations_dir.Write()

        OutputDir.Write(
            self.__output_dir_name, rt.TObject.kSingleKey + rt.TObject.kWriteDelete
        )

        OutputFile.Save()
        # OutputFile.Close()
