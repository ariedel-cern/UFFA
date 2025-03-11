import ROOT as rt
import logging
import copy
import numpy as np
import ctypes as ct


from . import TH1Plotter as th1p
from . import TGraphPlotter as tgp
from .Utils import AnalysisUtils as au

logger = logging.getLogger(__name__)


class PlotHandler:
    """
    PlotHandler
    """

    output_formats = [".pdf", ".png", ".svg"]

    def __init__(self, canvas_dict, legend_dict=None):
        """
        PlotHandler constructor

        Args:
            canvas_dict (dict): Dictionary to configure canvas
            hist_dict (dict): Dictionary to configure histograms
            legend_dict (dict): Dictionary to configure legends
        """
        # canvas style
        self._canvas_name = canvas_dict.get("Name", "canvas")
        self._canvas_title = canvas_dict.get("Title", "title")
        self._canvas_width = canvas_dict.get("Width", 700)
        self._canvas_height = canvas_dict.get("Height", 500)

        self._canvas_xmin = canvas_dict.get("Xmin", 0)
        self._canvas_xmax = canvas_dict.get("Xmax", 0)
        self._canvas_ymin = canvas_dict.get("Ymin", 0)
        self._canvas_ymax = canvas_dict.get("Ymax", 0)

        self._canvas_left_margin = canvas_dict.get("MarginLeft", 0.1)
        self._canvas_right_margin = canvas_dict.get("MarginRight", 0.1)
        self._canvas_bottom_margin = canvas_dict.get("MarginBottom", 0.1)
        self._canvas_top_margin = canvas_dict.get("MarginTop", 0.1)

        self._canvas_log = canvas_dict.get("Log", "")
        self._canvas_max_digits = canvas_dict.get("MaxDigits", 2)

        # axis style
        self._xaxis_title = canvas_dict.get("XTitle", "xaxis")
        self._yaxis_title = canvas_dict.get("YTitle", "yaxis")

        self._xaxis_title_size = canvas_dict.get("XTitleSize", 0.04)
        self._yaxis_title_size = canvas_dict.get("YTitleSize", 0.04)

        self._xaxis_label_size = canvas_dict.get("XLabelSize", 0.04)
        self._yaxis_label_size = canvas_dict.get("YLabelSize", 0.04)

        self._xaxis_bin_labels = canvas_dict.get("XBinLabels", [])
        self._yaxis_bin_labels = canvas_dict.get("YBinLabels", [])

        self._xaxis_title_offset = canvas_dict.get("XTitleOffset", 1)
        self._yaxis_title_offset = canvas_dict.get("YTitleOffset", 1)

        self._hist_dicts = []
        self._graph_dicts = []

        # legend style
        self._create_legend = False
        if legend_dict is not None:
            self._create_legend = True
            self._legend_xmin = legend_dict.get("Xmin", 0)
            self._legend_xmax = legend_dict.get("Xmax", 0)
            self._legend_ymin = legend_dict.get("Ymin", 0)
            self._legend_ymax = legend_dict.get("Ymax", 0)

            self._legend_bordersize = legend_dict.get("BorderSize", 0)
            self._legend_fillstyle = legend_dict.get("FillStyle", 0)
            self._legend_header = legend_dict.get("Header", "")

        self._output_name = canvas_dict.get("OutputName", "plot")
        self._canvas = None
        self._histograms = []
        self._graphs = []
        self._legend = None

    def SetHistograms(self, hist_dicts):
        self._hist_dicts = copy.deepcopy(hist_dicts)

    def SetGraphs(self, graph_dicts):
        self._graph_dicts = copy.deepcopy(graph_dicts)

    def CreateCanvas(self):
        """
        Create canvas for the plot
        """
        # set max digits per axis
        rt.TGaxis.SetMaxDigits(self._canvas_max_digits)

        # create canvas
        self._canvas = rt.TCanvas(
            self._canvas_name,
            self._canvas_title,
            self._canvas_width,
            self._canvas_height,
        )
        logger.debug("Create canvas with name: %s", self._canvas_name)
        logger.debug("Create canvas with title: %s", self._canvas_title)
        logger.debug(
            "Create canvas with dimensions (width,height): (%d,%d)",
            self._canvas_width,
            self._canvas_height,
        )

        # set margins
        self._canvas.SetMargin(
            self._canvas_left_margin,
            self._canvas_right_margin,
            self._canvas_bottom_margin,
            self._canvas_top_margin,
        )
        logger.debug("Set left Margin: %.2f", self._canvas_left_margin)
        logger.debug("Set right Margin: %.2f", self._canvas_right_margin)
        logger.debug("Set bottom Margin: %.2f", self._canvas_bottom_margin)
        logger.debug("Set top Margin: %.2f", self._canvas_top_margin)

        self._frame = self._canvas.DrawFrame(
            self._canvas_xmin,
            self._canvas_ymin,
            self._canvas_xmax,
            self._canvas_ymax,
        )
        logger.debug(
            "Create canvas with Frame Xmin: %.2f -> Xmax: %.2f and Ymin: %.2f -> Ymax: %.2f",
            self._canvas_xmin,
            self._canvas_xmax,
            self._canvas_ymin,
            self._canvas_ymax,
        )

        self._frame.SetTitle(self._canvas_title)
        self._frame.SetName("canvas_frame")

        if self._xaxis_bin_labels:
            self._frame.SetBins(
                len(self._xaxis_bin_labels), self._canvas_xmin, self._canvas_xmax
            )
            for i, name in enumerate(self._xaxis_bin_labels):
                self._frame.GetXaxis().SetBinLabel(i + 1, name)
                logger.debug("Set label %s on bin %d", name, i + 1)
                # self._hist.GetXaxis().ChangeLabel(i + 1, 320, -1, 13, -1, -1, name)

        self._frame.GetXaxis().SetLabelSize(self._xaxis_label_size)
        self._frame.GetXaxis().SetTitleSize(self._xaxis_title_size)
        self._frame.GetXaxis().SetTitleOffset(self._xaxis_title_offset)
        self._frame.GetXaxis().SetTitle(self._xaxis_title)
        logger.debug("X axis title: %s", self._xaxis_title)
        logger.debug("X axis title size: %.2f", self._xaxis_title_size)
        logger.debug("X axis label size: %.2f", self._xaxis_label_size)
        self._frame.GetYaxis().SetLabelSize(self._yaxis_label_size)
        self._frame.GetYaxis().SetTitleSize(self._yaxis_title_size)
        self._frame.GetYaxis().SetTitleOffset(self._yaxis_title_offset)
        self._frame.GetYaxis().SetTitle(self._yaxis_title)
        logger.debug("Y axis title: %s", self._yaxis_title)
        logger.debug("Y axis title size: %.2f", self._yaxis_title_size)
        logger.debug("Y axis label size: %.2f", self._yaxis_label_size)

        if "z" in self._canvas_log:
            self._canvas.SetLogz()
            logger.debug("Use log scale for z axis")
        if "y" in self._canvas_log:
            self._canvas.SetLogy()
            logger.debug("Use log scale for y axis")
        if "x" in self._canvas_log:
            self._canvas.SetLogx()
            logger.debug("Use log scale for x axis")

    def CreateLegend(self):
        """
        Add legend to the plot
        """
        # normalized coordinated
        if self._create_legend is False:
            return

        self._legend = rt.TLegend(
            self._legend_xmin,
            self._legend_ymin,
            self._legend_xmax,
            self._legend_ymax,
        )
        logger.debug(
            "Create legend at position Xmin: %.2f -> Xmax: %.2f and Ymin: %.2f -> Ymax: %.2f",
            self._legend_xmin,
            self._legend_xmax,
            self._legend_ymin,
            self._legend_ymax,
        )

        # defaults for legend
        self._legend.SetBorderSize(self._legend_bordersize)
        self._legend.SetFillStyle(self._legend_fillstyle)
        if self._legend_header != "":
            self._legend.SetHeader(self._legend_header)
            logger.debug("Add legend header: %s", self._legend_header)
        for hist in self._histograms:
            if hist.IsLegendDefined():
                logger.debug(
                    "Add histogram %s to legend with entry %s and option %s",
                    hist.GetHistogram().GetName(),
                    hist.GetLegendEntry(),
                    hist.GetLegendOption(),
                )
                self._legend.AddEntry(
                    hist.GetHistogram(),
                    hist.GetLegendEntry(),
                    hist.GetLegendOption(),
                )
        for graph in self._graphs:
            if graph.IsLegendDefined():
                logger.debug(
                    "Add histogram %s to legend with entry %s and option %s",
                    graph.GetGraph().GetName(),
                    graph.GetLegendEntry(),
                    graph.GetLegendOption(),
                )
                self._legend.AddEntry(
                    graph.GetGraph(),
                    graph.GetLegendEntry(),
                    graph.GetLegendOption(),
                )

    def FetchObjects(self):
        """
        Fetch histograms
        """
        logger.debug("Fetch histograms")
        # override titles and names
        for i, d in enumerate(self._hist_dicts):
            d.update(
                {
                    "Name": f"Histogram_{i}",
                }
            )
        # create wrappers
        self._histograms = [th1p.TH1Plotter(d) for d in self._hist_dicts]

        logger.debug("Fetch graphs")
        # override titles and names
        for i, d in enumerate(self._graph_dicts):
            d.update(
                {
                    "Name": f"Graph_{i}",
                }
            )
        self._graphs = [tgp.TGraphPlotter(d) for d in self._graph_dicts]

    def Plot(self):
        """
        Create plot overlaying all pass objects
        """
        logger.debug("Draw plot")
        # fetch histogram and
        self.FetchObjects()
        # create canvas and legend
        self.CreateCanvas()
        self.CreateLegend()
        # draw histograms with configured style
        for hist in self._histograms:
            hist.Draw(style=True)
        # draw graphs with configured style
        for graph in self._graphs:
            graph.Draw(style=True)
        # draw legend
        self._legend.Draw()
        # save output
        for format in self.output_formats:
            output_name = self._output_name + format
            au.CreateOutputDir(output_name, rename_old=True)
            self._canvas.SaveAs(output_name)
