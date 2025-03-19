import ROOT as rt
import logging
import numpy as np
import ctypes as ct
from .Utils import CorrelationUtils as cu
from .Utils import HistUtils as hu

logger = logging.getLogger(__name__)


class PurityHandler:
    """
    CorrelationHandler class
    """

    def __init__(self, hist_dict):

        if "Histogram" in hist_dict:
            self.__hist = hist_dict.get("Histogram").Clone()
            self.__hist.SetDirectory(0)
            logger.debug("Set histogram %s from dict", self.__hist.GetName())
        else:
            with rt.TFile(hist_dict.get("File", ""), "READ") as file:
                self.__hist = au.GetObjectFromFile(
                    file, hist_dict.get("Path", "path")
                ).Clone()
                self.__hist.SetDirectory(0)
                logger.debug(
                    "Open file %s to retrieve histogram at path %s",
                    hist_dict.get("File", ""),
                    hist_dict.get("Path", "path"),
                )


    :q

