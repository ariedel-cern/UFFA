import ROOT as rt
import logging

from .Utils import AnalysisUtils as au

logger = logging.getLogger(__name__)


class TGraphPlotter:
    """
    TGraphPlotter
    """

    def __init__(self, graph_dict):
        if "Graph" in graph_dict:
            self.__graph = graph_dict.get("Graph").Clone()
            logger.debug("Set graph %s from dict", graph_dict["Graph"].GetName())
        else:
            with rt.TFile(graph_dict.get("File"), "READ") as file:
                self.__graph = au.GetObjectFromFile(
                    file, graph_dict.get("Path")
                ).Clone()
                logger.debug(
                    "Open file %s to retrieve histogram at path %s",
                    graph_dict.get("File"),
                    graph_dict.get("Path"),
                )

        # file to fetch histogram (if needed)
        self.__filename = graph_dict.get("File", "")
        self.__path = graph_dict.get("Path", "")

        # name and title
        self.__histName = graph_dict.get("Name", "graph")
        self.__histTitle = graph_dict.get("Title", "graph")

        # axis settings
        self.__xaxisTitle = graph_dict.get("XTitle", "xaxis")
        self.__xaxisTitleSize = graph_dict.get("XTitleSize", 0.04)
        self.__xaxisLabelSize = graph_dict.get("XLabelSize", 0.04)

        self.__yaxisTitle = graph_dict.get("YTitle", "yaxis")
        self.__yaxisTitleSize = graph_dict.get("YTitleSize", 0.04)
        self.__yaxisLabelSize = graph_dict.get("YLabelSize", 0.04)

        # style settings
        self.__LineColor = graph_dict.get("LineColor", rt.kBlack)
        self.__LineAlpha = graph_dict.get("LineAlpha", -1)
        self.__LineStyle = graph_dict.get("LineStyle", 1)
        self.__LineWidth = graph_dict.get("LineWidth", 1)

        self.__MarkerColor = graph_dict.get("MarkerColor", rt.kBlack)
        self.__MarkerStyle = graph_dict.get("MarkerStyle", 1)
        self.__MarkerSize = graph_dict.get("MarkerSize", 1.0)

        self.__FillColor = graph_dict.get("FillColor", 0)
        self.__FillAlpha = graph_dict.get("FillAlpha", -1)
        self.__FillStyle = graph_dict.get("FillStyle", 1001)

        self.__draw_options = graph_dict.get("DrawOption", "P")

        self.__legend_defined = False
        # check for legend entry, this is the most important key
        if "LegendEntry" in graph_dict:
            self.__legend_defined = True
            self.__legend_entry = graph_dict.get("LegendEntry", "")
            self.__legend_option = graph_dict.get("LegendOption", "lepf")

    def Style(self):
        """
        Apply the stored style settings to the histogram
        """

        self.__graph.SetName(self.__histName)
        self.__graph.SetTitle(self.__histTitle)
        logger.debug("Histogram name %s", self.__histName)
        logger.debug("Histogram title %s", self.__histTitle)

        self.__graph.GetXaxis().SetTitle(self.__xaxisTitle)
        self.__graph.GetXaxis().SetTitleSize(self.__xaxisTitleSize)
        self.__graph.GetXaxis().SetLabelSize(self.__xaxisLabelSize)
        logger.debug("X axis tile %s", self.__xaxisTitle)
        logger.debug("X axis title size %.2f", self.__xaxisTitleSize)
        logger.debug("X axis label size %.2f", self.__xaxisLabelSize)

        self.__graph.GetYaxis().SetTitle(self.__yaxisTitle)
        self.__graph.GetYaxis().SetTitleSize(self.__yaxisTitleSize)
        self.__graph.GetYaxis().SetLabelSize(self.__yaxisLabelSize)
        logger.debug("Y axis tile %s", self.__yaxisTitle)
        logger.debug("Y axis title size %.2f", self.__yaxisTitleSize)
        logger.debug("Y axis label size %.2f", self.__yaxisLabelSize)

        if self.__LineAlpha > 0:
            self.__graph.SetLineColorAlpha(self.__LineColor, self.__LineAlpha)
        else:
            self.__graph.SetLineColor(self.__LineColor)
        self.__graph.SetLineStyle(self.__LineStyle)
        self.__graph.SetLineWidth(self.__LineWidth)
        logger.debug("Line color: %d", self.__LineColor)
        logger.debug("Line style: %d", self.__LineStyle)
        logger.debug("Line width: %d", self.__LineWidth)

        self.__graph.SetMarkerColor(self.__MarkerColor)
        self.__graph.SetMarkerStyle(self.__MarkerStyle)
        self.__graph.SetMarkerSize(self.__MarkerSize)
        logger.debug("Marker color: %d", self.__MarkerColor)
        logger.debug("Marker style: %d", self.__MarkerStyle)
        logger.debug("Marker size: %.2f", self.__MarkerSize)

        if self.__FillAlpha > 0:
            self.__graph.SetFillColorAlpha(self.__FillColor, self.__FillAlpha)
            logger.debug(
                "Fill color: %d with Alpha: %.2f", self.__FillColor, self.__FillAlpha
            )
        else:
            self.__graph.SetFillColor(self.__FillColor)
            logger.debug("Fill color: %d", self.__FillColor)
        self.__graph.SetFillStyle(self.__FillStyle)
        logger.debug("Fill color: %d", self.__FillColor)
        logger.debug("Fill style: %d", self.__FillStyle)

    def Draw(self, style=True):
        """
        Draw styled histogram
        """
        logger.debug("Style and draw histogram")
        if style:
            self.Style()
        self.__graph.Draw(self.__draw_options)
        logger.debug("Draw options: %s", self.__draw_options)

    def GetGraph(self):
        """
        Return histogram
        """
        return self.__graph

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
