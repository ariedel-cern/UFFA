# UFFA/__init__.py

# check if ROOT is present
try:
    import ROOT
except ImportError as e:
    raise ImportError(
        "UFFA requires PyROOT (ROOT's Python bindings). "
        "Please install ROOT (e.g. via aliBuild, alisw, or from https://root.cern)."
    ) from e

# main packages
from .AnalysisHandler import AnalysisHandler
from .CorrelationHandler import CorrelationHandler
from .SystematicHandler import SystematicHandler
from .FemtoMerger import FemtoMerger

__all__ = [
    "AnalysisHandler",
    "CorrelationHandler",
    "SystematicHandler",
]

# sub packages
from .Utils import AnalysisUtils, CorrelationUtils, HistUtils

__all__.extend(["AnalysisUtils", "CorrelationUtils", "HistUtils"])
