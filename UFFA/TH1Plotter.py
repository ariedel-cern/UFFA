import ROOT as rt

rt.TH1.AddDirectory(False)

import logging

from .Utils import AnalysisUtils as au

logger = logging.getLogger(__name__)


class TH1Plotter:
    """
    TH1Plotter
    """

    def __init__(self, hist_dict):

        if "Histogram" in hist_dict:
            self.__hist = hist_dict.get("Histogram").Clone()
            self.__hist.SetDirectory(0)
            logger.debug("Set histogram %s from dict", self.__hist.GetName())
        else:
            au.GetObject(hist_dict.get("File", ""), hist_dict.get("Path", "path"))
            logger.debug(
                "Open file %s to retrieve histogram at path %s",
                hist_dict.get("File", ""),
                hist_dict.get("Path", "path"),
            )

        # name and title
        self.__histName = hist_dict.get("Name", "histogram")
        self.__histTitle = hist_dict.get("Title", "histogram")

        # axis settings
        self.__xaxisTitle = hist_dict.get("XTitle", "xaxis")
        self.__xaxisTitleSize = hist_dict.get("XTitleSize", 0.04)
        self.__xaxisLabelSize = hist_dict.get("XLabelSize", 0.04)

        self.__yaxisTitle = hist_dict.get("YTitle", "yaxis")
        self.__yaxisTitleSize = hist_dict.get("YTitleSize", 0.04)
        self.__yaxisLabelSize = hist_dict.get("YLabelSize", 0.04)

        # style settings
        self.__LineColor = hist_dict.get("LineColor", rt.kBlack)
        self.__LineStyle = hist_dict.get("LineStyle", 1)
        self.__LineWidth = hist_dict.get("LineWidth", 1)

        self.__MarkerColor = hist_dict.get("MarkerColor", rt.kBlack)
        self.__MarkerStyle = hist_dict.get("MarkerStyle", 1)
        self.__MarkerSize = hist_dict.get("MarkerSize", 1.0)

        self.__FillColor = hist_dict.get("FillColor", 0)
        self.__FillAlpha = hist_dict.get("FillAlpha", -1)
        self.__FillStyle = hist_dict.get("FillStyle", 1001)

        self.__Norm = hist_dict.get("Norm", -1)

        self.__draw_options = hist_dict.get("DrawOption", "SAME")

        self.__legend_defined = False
        # check for legend entry, this is the most important key
        if "LegendEntry" in hist_dict:
            self.__legend_defined = True
            self.__legend_entry = hist_dict.get("LegendEntry", "")
            self.__legend_option = hist_dict.get("LegendOption", "lepf")

    def Style(self):
        """
        Apply the stored style settings to the histogram
        """

        if self.__Norm > 0:
            self.__hist.Scale(self.__Norm, "nosw2")

        self.__hist.SetName(self.__histName)
        self.__hist.SetTitle(self.__histTitle)
        logger.debug("Histogram name %s", self.__histName)
        logger.debug("Histogram title %s", self.__histTitle)

        self.__hist.GetXaxis().SetTitle(self.__xaxisTitle)
        self.__hist.GetXaxis().SetTitleSize(self.__xaxisTitleSize)
        self.__hist.GetXaxis().SetLabelSize(self.__xaxisLabelSize)
        logger.debug("X axis tile %s", self.__xaxisTitle)
        logger.debug("X axis title size %.2f", self.__xaxisTitleSize)
        logger.debug("X axis label size %.2f", self.__xaxisLabelSize)

        self.__hist.GetYaxis().SetTitle(self.__yaxisTitle)
        self.__hist.GetYaxis().SetTitleSize(self.__yaxisTitleSize)
        self.__hist.GetYaxis().SetLabelSize(self.__yaxisLabelSize)
        logger.debug("Y axis tile %s", self.__yaxisTitle)
        logger.debug("Y axis title size %.2f", self.__yaxisTitleSize)
        logger.debug("Y axis label size %.2f", self.__yaxisLabelSize)

        self.__hist.SetLineColor(self.__LineColor)
        self.__hist.SetLineStyle(self.__LineStyle)
        self.__hist.SetLineWidth(self.__LineWidth)
        logger.debug("Line color: %d", self.__LineColor)
        logger.debug("Line style: %d", self.__LineStyle)
        logger.debug("Line width: %d", self.__LineWidth)

        self.__hist.SetMarkerColor(self.__MarkerColor)
        self.__hist.SetMarkerStyle(self.__MarkerStyle)
        self.__hist.SetMarkerSize(self.__MarkerSize)
        logger.debug("Marker color: %d", self.__MarkerColor)
        logger.debug("Marker style: %d", self.__MarkerStyle)
        logger.debug("Marker size: %.2f", self.__MarkerSize)

        if self.__FillAlpha > 0:
            self.__hist.SetFillColorAlpha(self.__FillColor, self.__FillAlpha)
            logger.debug(
                "Fill color: %d with Alpha: %.2f", self.__FillColor, self.__FillAlpha
            )
        else:
            self.__hist.SetFillColor(self.__FillColor)
            logger.debug("Fill color: %d", self.__FillColor)
        self.__hist.SetFillStyle(self.__FillStyle)
        logger.debug("Fill style: %d", self.__FillStyle)

    def Draw(self, style=True):
        """
        Draw styled histogram
        """
        logger.debug("Style and draw histogram")
        if style:
            self.Style()
        self.__hist.Draw(self.__draw_options)
        logger.debug("Draw options: %s", self.__draw_options)

    def GetHistogram(self):
        """
        Return histogram
        """
        return self.__hist

    def GetLegendEntry(self):
        """
        Return legend entry
        """
        logger.debug(
            "Histogram with %s with legend Entry: %s",
            self.__histName,
            self.__legend_entry,
        )
        return self.__legend_entry

    def GetLegendOption(self):
        """
        Return legend entry
        """
        logger.debug(
            "Histogram with %s with legend opton: %s",
            self.__histName,
            self.__legend_option,
        )
        return self.__legend_option

    def IsLegendDefined(self):
        """
        Whether legend entry should be added
        """
        return self.__legend_defined
