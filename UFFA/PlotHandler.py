import ROOT as rt

rt.gROOT.SetBatch(True)
rt.EnableThreadSafety()
rt.TH1.AddDirectory(False)

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
        self.__canvas_name = canvas_dict.get("Name", "canvas")
        self.__canvas_title = canvas_dict.get("Title", "title")
        self.__canvas_width = canvas_dict.get("Width", 700)
        self.__canvas_height = canvas_dict.get("Height", 500)

        self.__canvas_xmin = canvas_dict.get("Xmin", 0)
        self.__canvas_xmax = canvas_dict.get("Xmax", 0)
        self.__canvas_ymin = canvas_dict.get("Ymin", 0)
        self.__canvas_ymax = canvas_dict.get("Ymax", 0)

        self.__canvas_left_margin = canvas_dict.get("MarginLeft", 0.1)
        self.__canvas_right_margin = canvas_dict.get("MarginRight", 0.1)
        self.__canvas_bottom_margin = canvas_dict.get("MarginBottom", 0.1)
        self.__canvas_top_margin = canvas_dict.get("MarginTop", 0.1)

        self.__canvas_log = canvas_dict.get("Log", "")
        self.__canvas_max_digits = canvas_dict.get("MaxDigits", 2)

        # axis style
        self.__xaxis_title = canvas_dict.get("XTitle", "xaxis")
        self.__yaxis_title = canvas_dict.get("YTitle", "yaxis")

        self.__xaxis_title_size = canvas_dict.get("XTitleSize", 0.07)
        self.__yaxis_title_size = canvas_dict.get("YTitleSize", 0.07)

        self.__xaxis_label_size = canvas_dict.get("XLabelSize", 0.05)
        self.__yaxis_label_size = canvas_dict.get("YLabelSize", 0.05)

        self.__xaxis_bin_labels = canvas_dict.get("XBinLabels", [])
        # self.__yaxis_bin_labels = canvas_dict.get("YBinLabels", [])

        self.__xaxis_title_offset = canvas_dict.get("XTitleOffset", 1)
        self.__yaxis_title_offset = canvas_dict.get("YTitleOffset", 1)
        self.__ticks_x = canvas_dict.get("TicksX", 1)
        self.__ticks_y = canvas_dict.get("TicksY", 1)

        self.__end_error_size = canvas_dict.get("EndErrorSize", 1)

        self.__hist_dicts = []
        self.__graph_dicts = []
        self.__line_dicts = []
        self.__text_dicts = []

        # legend style
        self.__create_legend = False
        if legend_dict is not None:
            self.__create_legend = True
            self.__legend_xmin = legend_dict.get("Xmin", 0)
            self.__legend_xmax = legend_dict.get("Xmax", 0)
            self.__legend_ymin = legend_dict.get("Ymin", 0)
            self.__legend_ymax = legend_dict.get("Ymax", 0)

            self.__legend_font = legend_dict.get("Font", 42)
            self.__legend_textsize = legend_dict.get("TextSize", 0)
            self.__legend_bordersize = legend_dict.get("BorderSize", 0)
            self.__legend_fillcolor = legend_dict.get("FillColor", 0)
            self.__legend_fillstyle = legend_dict.get("FillStyle", 0)
            self.__legend_columns = legend_dict.get("NColumns", 1)
            self.__legend_header = legend_dict.get("Header", "")

        self.__output_name = canvas_dict.get("OutputName", "plot")
        self.__canvas = None
        self.__histograms = []
        self.__graphs = []
        self.__lines = []
        self.__texts = []
        self.__legend = None

    def SetHistograms(self, hist_dicts):
        self.__hist_dicts = copy.deepcopy(hist_dicts)

    def SetHistogramsWithRatio(self, ratio_dict, hist_dicts):
        ratioHistDict = copy.deepcopy(ratio_dict)

        if "Histogram" in ratioHistDict:
            RatioHist = ratio_dict.get("Histogram").Clone()
            logger.debug("Get histogram %s from dict", RatioHist.GetName())
        else:
            with rt.TFile(ratio_dict.get("File", ""), "READ") as file:
                RatioHist = au.GetObjectFromFile(
                    file, ratio_dict.get("Path", "path")
                ).Clone()
                logger.debug(
                    "Open file %s to retrieve histogram at path %s",
                    ratio_dict.get("File", ""),
                    ratio_dict.get("Path", "path"),
                )

        ratioHistDicts = copy.deepcopy(hist_dicts)
        for hist_dict in ratioHistDicts:
            if "Histogram" in ratioHistDict:
                TempHist = hist_dict.get("Histogram").Clone("Ratio")
            else:
                with rt.TFile(hist_dict.get("File", ""), "READ") as file:
                    TempHist = au.GetObjectFromFile(
                        file, hist_dict.get("Path", "Path")
                    ).Clone("Ratio")
                    TempHist.SetDirectory(0)
                del hist_dict["File"]
                del hist_dict["Path"]

            TempHist.Divide(RatioHist)
            hist_dict["Histogram"] = TempHist
            self.__hist_dicts.append(hist_dict)

    def SetGraphs(self, graph_dicts):
        self.__graph_dicts = copy.deepcopy(graph_dicts)

    def SetLines(self, line_dicts):
        self.__line_dicts = copy.deepcopy(line_dicts)

    def SetTexts(self, text_dicts):
        self.__text_dicts = copy.deepcopy(text_dicts)

    def CreateCanvas(self):
        """
        Create canvas for the plot
        """
        # set max digits per axis
        rt.TGaxis.SetMaxDigits(self.__canvas_max_digits)

        # create canvas
        self.__canvas = rt.TCanvas(
            self.__canvas_name,
            self.__canvas_title,
            self.__canvas_width,
            self.__canvas_height,
        )
        logger.debug("Create canvas with name: %s", self.__canvas_name)
        logger.debug("Create canvas with title: %s", self.__canvas_title)
        logger.debug(
            "Create canvas with dimensions (width,height): (%d,%d)",
            self.__canvas_width,
            self.__canvas_height,
        )

        # set margins
        self.__canvas.SetMargin(
            self.__canvas_left_margin,
            self.__canvas_right_margin,
            self.__canvas_bottom_margin,
            self.__canvas_top_margin,
        )
        logger.debug("Set left Margin: %.2f", self.__canvas_left_margin)
        logger.debug("Set right Margin: %.2f", self.__canvas_right_margin)
        logger.debug("Set bottom Margin: %.2f", self.__canvas_bottom_margin)
        logger.debug("Set top Margin: %.2f", self.__canvas_top_margin)

        self.__frame = self.__canvas.DrawFrame(
            self.__canvas_xmin,
            self.__canvas_ymin,
            self.__canvas_xmax,
            self.__canvas_ymax,
        )
        logger.debug(
            "Create canvas with Frame Xmin: %.2f -> Xmax: %.2f and Ymin: %.2f -> Ymax: %.2f",
            self.__canvas_xmin,
            self.__canvas_xmax,
            self.__canvas_ymin,
            self.__canvas_ymax,
        )

        self.__frame.SetTitle(self.__canvas_title)
        self.__frame.SetName("canvas_frame")

        if self.__xaxis_bin_labels:
            self.__frame.SetBins(
                len(self.__xaxis_bin_labels), self.__canvas_xmin, self.__canvas_xmax
            )
            for i, name in enumerate(self.__xaxis_bin_labels):
                self.__frame.GetXaxis().SetBinLabel(i + 1, name)
                logger.debug("Set label %s on bin %d", name, i + 1)
                # self.__hist.GetXaxis().ChangeLabel(i + 1, 320, -1, 13, -1, -1, name)

        self.__frame.GetXaxis().SetLabelSize(self.__xaxis_label_size)
        self.__frame.GetXaxis().SetTitleSize(self.__xaxis_title_size)
        self.__frame.GetXaxis().SetTitleOffset(self.__xaxis_title_offset)
        self.__frame.GetXaxis().SetTitle(self.__xaxis_title)
        logger.debug("X axis title: %s", self.__xaxis_title)
        logger.debug("X axis title size: %.2f", self.__xaxis_title_size)
        logger.debug("X axis label size: %.2f", self.__xaxis_label_size)
        self.__frame.GetYaxis().SetLabelSize(self.__yaxis_label_size)
        self.__frame.GetYaxis().SetTitleSize(self.__yaxis_title_size)
        self.__frame.GetYaxis().SetTitleOffset(self.__yaxis_title_offset)
        self.__frame.GetYaxis().SetTitle(self.__yaxis_title)
        logger.debug("Y axis title: %s", self.__yaxis_title)
        logger.debug("Y axis title size: %.2f", self.__yaxis_title_size)
        logger.debug("Y axis label size: %.2f", self.__yaxis_label_size)

        if "z" in self.__canvas_log:
            self.__canvas.SetLogz()
            logger.debug("Use log scale for z axis")
        if "y" in self.__canvas_log:
            self.__canvas.SetLogy()
            logger.debug("Use log scale for y axis")
        if "x" in self.__canvas_log:
            self.__canvas.SetLogx()
            logger.debug("Use log scale for x axis")

    def CreateLegend(self):
        """
        Add legend to the plot
        """
        # normalized coordinated
        if self.__create_legend is False:
            return

        self.__legend = rt.TLegend(
            self.__legend_xmin,
            self.__legend_ymin,
            self.__legend_xmax,
            self.__legend_ymax,
        )
        logger.debug(
            "Create legend at position Xmin: %.2f -> Xmax: %.2f and Ymin: %.2f -> Ymax: %.2f",
            self.__legend_xmin,
            self.__legend_xmax,
            self.__legend_ymin,
            self.__legend_ymax,
        )

        # defaults for legend
        rt.gStyle.SetLegendFont(self.__legend_font)
        rt.gStyle.SetLegendTextSize(self.__legend_textsize)
        rt.gStyle.SetLegendFillColor(self.__legend_fillcolor)

        # self.__texts[-1].SetTextAlign(text_dict.get("Align", 11))
        # self.__legend.SetTextFont(self.__legend_font)
        # self.__legend.SetTextSize(self.__legend_textsize)

        self.__legend.SetBorderSize(self.__legend_bordersize)
        self.__legend.SetFillStyle(self.__legend_fillstyle)
        self.__legend.SetNColumns(self.__legend_columns)
        if self.__legend_header != "":
            self.__legend.SetHeader(self.__legend_header)
            logger.debug("Add legend header: %s", self.__legend_header)
        for hist in self.__histograms:
            if hist.IsLegendDefined():
                logger.debug(
                    "Add histogram %s to legend with entry %s and option %s",
                    hist.GetHistogram().GetName(),
                    hist.GetLegendEntry(),
                    hist.GetLegendOption(),
                )
                self.__legend.AddEntry(
                    hist.GetHistogram(),
                    hist.GetLegendEntry(),
                    hist.GetLegendOption(),
                )
        for graph in self.__graphs:
            if graph.IsLegendDefined():
                logger.debug(
                    "Add histogram %s to legend with entry %s and option %s",
                    graph.GetGraph().GetName(),
                    graph.GetLegendEntry(),
                    graph.GetLegendOption(),
                )
                self.__legend.AddEntry(
                    graph.GetGraph(),
                    graph.GetLegendEntry(),
                    graph.GetLegendOption(),
                )

    def DrawLines(self):
        """
        Draw TLineObjects objects to plot text in the plot
        """
        for line_dict in self.__line_dicts:
            self.__lines.append(rt.TLine())
            self.__lines[-1].SetLineWidth(line_dict.get("Width", 3))
            self.__lines[-1].SetLineStyle(line_dict.get("Style", 1))
            if "Alpha" in line_dict:
                self.__lines[-1].SetLineColorAlpha(
                    line_dict.get("Color", 1), line_dict.get("Alpha", 1)
                )
            else:
                self.__lines[-1].SetLineColor(line_dict.get("Color", 1))
            self.__lines[-1].DrawLine(
                line_dict.get("Xmin", 0),
                line_dict.get("Ymin", 1),
                line_dict.get("Xmax", 0),
                line_dict.get("Ymax", 1),
            )

    def DrawTexts(self):
        """
        Draw TLatex objects on plot
        """
        for text_dict in self.__text_dicts:
            self.__texts.append(rt.TLatex())
            self.__texts[-1].SetTextAlign(text_dict.get("Align", 11))
            self.__texts[-1].SetTextFont(text_dict.get("Font", 42))
            self.__texts[-1].SetTextSize(text_dict.get("Size", 0.08))
            self.__texts[-1].DrawLatexNDC(
                text_dict.get("X", 0.5),
                text_dict.get("Y", 0.5),
                text_dict.get("Text", "text"),
            )

    def FetchObjects(self):
        """
        Fetch histograms
        """
        logger.debug("Fetch histograms")
        self.__histograms = [th1p.TH1Plotter(d) for d in self.__hist_dicts]

        logger.debug("Fetch graphs")
        self.__graphs = [tgp.TGraphPlotter(d) for d in self.__graph_dicts]

    def Plot(self):
        """
        Create plot overlaying all pass objects
        """
        logger.debug("Draw plot")
        # global style options
        rt.gStyle.SetEndErrorSize(self.__end_error_size)
        rt.gStyle.SetPadTickX(self.__ticks_x)
        rt.gStyle.SetPadTickY(self.__ticks_y)

        # fetch histogram and graphs
        self.FetchObjects()
        # create canvas
        self.CreateCanvas()
        # create legend
        self.CreateLegend()
        # create texts
        self.DrawTexts()
        # draw lines
        self.DrawLines()
        # draw histograms with configured style
        for hist in self.__histograms:
            hist.Draw(style=True)
        # draw graphs with configured style
        for graph in self.__graphs:
            graph.Draw(style=True)
        # draw text
        # for text in self.__texts:
        #     text.Draw()
        # draw legend
        if self.__create_legend:
            self.__legend.Draw()
        # save output
        for format in self.output_formats:
            output_name = self.__output_name + format
            au.CreateOutputDir(output_name, rename_old=True)
            self.__canvas.SaveAs(output_name)
