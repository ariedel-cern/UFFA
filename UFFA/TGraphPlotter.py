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
            logger.debug("Set graph %s from dict", graph_dict["Graph"].GetName())
            self._graph = graph_dict.get("Graph").Clone()
        else:
            with rt.TFile(graph_dict.get("File"), "READ") as file:
                self._graph = au.GetObjectFromFile(file, graph_dict.get("Path")).Clone()
                logger.debug(
                    "Open file %s to retrieve histogram at path %s",
                    graph_dict.get("File"),
                    graph_dict.get("Path"),
                )

        # file to fetch histogram (if needed)
        self._filename = graph_dict.get("File", "")
        self._path = graph_dict.get("Path", "")

        # name and title
        self._histName = graph_dict.get("Name", "graph")
        self._histTitle = graph_dict.get("Title", "graph")

        # axis settings
        self._xaxisTitle = graph_dict.get("XTitle", "xaxis")
        self._xaxisTitleSize = graph_dict.get("XTitleSize", 0.04)
        self._xaxisLabelSize = graph_dict.get("XLabelSize", 0.04)

        self._yaxisTitle = graph_dict.get("YTitle", "yaxis")
        self._yaxisTitleSize = graph_dict.get("YTitleSize", 0.04)
        self._yaxisLabelSize = graph_dict.get("YLabelSize", 0.04)

        # style settings
        self._LineColor = graph_dict.get("LineColor", rt.kBlack)
        self._LineStyle = graph_dict.get("LineStyle", 1)
        self._LineWidth = graph_dict.get("LineWidth", 1)

        self._MarkerColor = graph_dict.get("MarkerColor", rt.kBlack)
        self._MarkerStyle = graph_dict.get("MarkerStyle", 1)
        self._MarkerSize = graph_dict.get("MarkerSize", 1.0)

        self._FillColor = graph_dict.get("FillColor", 0)
        self._FillStyle = graph_dict.get("FillStyle", 1001)

        self._draw_options = graph_dict.get("DrawOption", "P")

        self._legend_defined = False
        # check for legend entry, this is the most important key
        if "LegendEntry" in graph_dict:
            self._legend_defined = True
            self._legend_entry = graph_dict.get("LegendEntry", "")
            self._legend_option = graph_dict.get("LegendOption", "lepf")

    def Style(self):
        """
        Apply the stored style settings to the histogram
        """

        self._graph.SetName(self._histName)
        self._graph.SetTitle(self._histTitle)
        logger.debug("Histogram name %s", self._histName)
        logger.debug("Histogram title %s", self._histTitle)

        self._graph.GetXaxis().SetTitle(self._xaxisTitle)
        self._graph.GetXaxis().SetTitleSize(self._xaxisTitleSize)
        self._graph.GetXaxis().SetLabelSize(self._xaxisLabelSize)
        logger.debug("X axis tile %s", self._xaxisTitle)
        logger.debug("X axis title size %.2f", self._xaxisTitleSize)
        logger.debug("X axis label size %.2f", self._xaxisLabelSize)

        self._graph.GetYaxis().SetTitle(self._yaxisTitle)
        self._graph.GetYaxis().SetTitleSize(self._yaxisTitleSize)
        self._graph.GetYaxis().SetLabelSize(self._yaxisLabelSize)
        logger.debug("Y axis tile %s", self._yaxisTitle)
        logger.debug("Y axis title size %.2f", self._yaxisTitleSize)
        logger.debug("Y axis label size %.2f", self._yaxisLabelSize)

        self._graph.SetLineColor(self._LineColor)
        self._graph.SetLineStyle(self._LineStyle)
        self._graph.SetLineWidth(self._LineWidth)
        logger.debug("Line color: %d", self._LineColor)
        logger.debug("Line style: %d", self._LineStyle)
        logger.debug("Line width: %d", self._LineWidth)

        self._graph.SetMarkerColor(self._MarkerColor)
        self._graph.SetMarkerStyle(self._MarkerStyle)
        self._graph.SetMarkerSize(self._MarkerSize)
        logger.debug("Marker color: %d", self._MarkerColor)
        logger.debug("Marker style: %d", self._MarkerStyle)
        logger.debug("Marker size: %.2f", self._MarkerSize)

        self._graph.SetFillColor(self._FillColor)
        self._graph.SetFillStyle(self._FillStyle)
        logger.debug("Fill color: %d", self._FillColor)
        logger.debug("Fill style: %d", self._FillStyle)

    def Draw(self, style=True):
        """
        Draw styled histogram
        """
        logger.debug("Style and draw histogram")
        if style:
            self.Style()
        self._graph.Draw(self._draw_options)
        logger.debug("Draw options: %s", self._draw_options)

    def GetGraph(self):
        """
        Return histogram
        """
        return self._graph

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
