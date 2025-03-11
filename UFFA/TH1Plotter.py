import ROOT as rt
import logging

from .Utils import AnalysisUtils as au

logger = logging.getLogger(__name__)


class TH1Plotter:
    """
    TH1Plotter
    """

    def __init__(self, hist_dict):

        if "Histogram" in hist_dict:
            self._hist = hist_dict.get("Histogram").Clone()
            self._hist.SetDirectory(0)
            print(self._hist)
            logger.debug("Set histogram %s from dict", self._hist.GetName())
        else:
            with rt.TFile(hist_dict.get("File", ""), "READ") as file:
                self._hist = au.GetObjectFromFile(
                    file, hist_dict.get("Path", "path")
                ).Clone()
                self._hist.SetDirectory(0)
                logger.debug(
                    "Open file %s to retrieve histogram at path %s",
                    hist_dict.get("File", ""),
                    hist_dict.get("Path", "path"),
                )

        # name and title
        self._histName = hist_dict.get("Name", "histogram")
        self._histTitle = hist_dict.get("Title", "histogram")

        # axis settings
        self._xaxisTitle = hist_dict.get("XTitle", "xaxis")
        self._xaxisTitleSize = hist_dict.get("XTitleSize", 0.04)
        self._xaxisLabelSize = hist_dict.get("XLabelSize", 0.04)

        self._yaxisTitle = hist_dict.get("YTitle", "yaxis")
        self._yaxisTitleSize = hist_dict.get("YTitleSize", 0.04)
        self._yaxisLabelSize = hist_dict.get("YLabelSize", 0.04)

        # style settings
        self._LineColor = hist_dict.get("LineColor", rt.kBlack)
        self._LineStyle = hist_dict.get("LineStyle", 1)
        self._LineWidth = hist_dict.get("LineWidth", 1)

        self._MarkerColor = hist_dict.get("MarkerColor", rt.kBlack)
        self._MarkerStyle = hist_dict.get("MarkerStyle", 1)
        self._MarkerSize = hist_dict.get("MarkerSize", 1.0)

        self._FillColor = hist_dict.get("FillColor", 0)
        self._FillStyle = hist_dict.get("FillStyle", 1001)

        self._Norm = hist_dict.get("Norm", -1)

        self._draw_options = hist_dict.get("DrawOption", "SAME")

        self._legend_defined = False
        # check for legend entry, this is the most important key
        if "LegendEntry" in hist_dict:
            self._legend_defined = True
            self._legend_entry = hist_dict.get("LegendEntry", "")
            self._legend_option = hist_dict.get("LegendOption", "lepf")

    def Style(self):
        """
        Apply the stored style settings to the histogram
        """

        if self._Norm > 0:
            self._hist.Scale(self._Norm, "nosw2")

        self._hist.SetName(self._histName)
        self._hist.SetTitle(self._histTitle)
        logger.debug("Histogram name %s", self._histName)
        logger.debug("Histogram title %s", self._histTitle)

        self._hist.GetXaxis().SetTitle(self._xaxisTitle)
        self._hist.GetXaxis().SetTitleSize(self._xaxisTitleSize)
        self._hist.GetXaxis().SetLabelSize(self._xaxisLabelSize)
        logger.debug("X axis tile %s", self._xaxisTitle)
        logger.debug("X axis title size %.2f", self._xaxisTitleSize)
        logger.debug("X axis label size %.2f", self._xaxisLabelSize)

        self._hist.GetYaxis().SetTitle(self._yaxisTitle)
        self._hist.GetYaxis().SetTitleSize(self._yaxisTitleSize)
        self._hist.GetYaxis().SetLabelSize(self._yaxisLabelSize)
        logger.debug("Y axis tile %s", self._yaxisTitle)
        logger.debug("Y axis title size %.2f", self._yaxisTitleSize)
        logger.debug("Y axis label size %.2f", self._yaxisLabelSize)

        self._hist.SetLineColor(self._LineColor)
        self._hist.SetLineStyle(self._LineStyle)
        self._hist.SetLineWidth(self._LineWidth)
        logger.debug("Line color: %d", self._LineColor)
        logger.debug("Line style: %d", self._LineStyle)
        logger.debug("Line width: %d", self._LineWidth)

        self._hist.SetMarkerColor(self._MarkerColor)
        self._hist.SetMarkerStyle(self._MarkerStyle)
        self._hist.SetMarkerSize(self._MarkerSize)
        logger.debug("Marker color: %d", self._MarkerColor)
        logger.debug("Marker style: %d", self._MarkerStyle)
        logger.debug("Marker size: %.2f", self._MarkerSize)

        self._hist.SetFillColor(self._FillColor)
        self._hist.SetFillStyle(self._FillStyle)
        logger.debug("Fill color: %d", self._FillColor)
        logger.debug("Fill style: %d", self._FillStyle)

    def Draw(self, style=True):
        """
        Draw styled histogram
        """
        logger.debug("Style and draw histogram")
        if style:
            self.Style()
        self._hist.Draw(self._draw_options)
        logger.debug("Draw options: %s", self._draw_options)

    def GetHistogram(self):
        """
        Return histogram
        """
        return self._hist

    def GetLegendEntry(self):
        """
        Return legend entry
        """
        logger.debug(
            "Histogram with %s with legend Entry: %s",
            self._histName,
            self._legend_entry,
        )
        return self._legend_entry

    def GetLegendOption(self):
        """
        Return legend entry
        """
        logger.debug(
            "Histogram with %s with legend opton: %s",
            self._histName,
            self._legend_option,
        )
        return self._legend_option

    def IsLegendDefined(self):
        """
        Whether legend entry should be added
        """
        return self._legend_defined
