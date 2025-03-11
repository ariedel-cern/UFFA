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


    def __init__(self):
