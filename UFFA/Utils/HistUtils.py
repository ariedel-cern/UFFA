import ROOT as rt
import numpy as np


def FindBinWithUpperEdgeDetection(Axis, Value):
    """
    Find bin of a TAxis for a given value with upper edge detection
    Args:
        Axis (TAxis): Axis with bins
        Value (float): Value whose bin is being searched

    Returns:

    """
    FoundBin = Axis.FindBin(Value)
    edge = Axis.GetBinLowEdge(FoundBin)
    if abs(Value - edge) < 1e-7 and FoundBin != 1:
        FoundBin -= 1
    return FoundBin


def GetHistDimension(hist):
    """
    Return Dimensions of histogram
    """
    if hist.InheritsFrom(rt.THnBase.Class()):
        return hist.GetNdimensions()
    elif hist.InheritsFrom(rt.TH3.Class()):
        return 3
    elif hist.InheritsFrom(rt.TH2.Class()):
        return 2
    elif hist.InheritsFrom(rt.TH1.Class()):
        return 1
    else:
        return 0


def SetHistRanges(hist, ranges):
    """
    Apply range limits to a histogram
    """
    dimension = GetHistDimension(hist)

    if dimension == 1:
        if ranges[0]:
            hist.GetXaxis().SetRangeUser(ranges[0][0], ranges[0][1])
    elif dimension == 2:
        if ranges[0]:
            hist.GetXaxis().SetRangeUser(ranges[0][0], ranges[0][1])
        if ranges[1]:
            hist.GetYaxis().SetRangeUser(ranges[1][0], ranges[1][1])
    elif dimension == 3:
        if ranges[0]:
            hist.GetXaxis().SetRangeUser(ranges[0][0], ranges[0][1])
        if ranges[1]:
            hist.GetYaxis().SetRangeUser(ranges[1][0], ranges[1][1])
        if ranges[2]:
            hist.GetZaxis().SetRangeUser(ranges[2][0], ranges[2][1])
    else:
        for dim, cut in enumerate(ranges):
            if cut:
                bin_low = hist.GetAxis(dim).FindBin(cut[0])
                bin_high = FindBinWithUpperEdgeDetection(hist.GetAxis(dim), cut[1])
                hist.GetAxis(dim).SetRange(bin_low, bin_high)


def RescaleHist(hist, scale, axis, suffix="_Rescaled"):
    """
    Rescale axis of a histogram
    """

    dimension = GetHistDimension(hist)
    name = hist.GetName() + f"{suffix}"
    hist_rescaled = None

    if axis < 0 or axis >= dimension:
        raise ValueError(f"Invalid axis {axis} for {dimension}D histogram")

    if dimension == 1:
        bins = hist.GetNbinsX()
        x_min = hist.GetXaxis().GetXmin() * scale
        x_max = hist.GetXaxis().GetXmax() * scale
        hist_rescaled = rt.TH1F(name, name, bins, x_min, x_max)
        for i in range(1, bins + 1):
            bin_content = hist.GetBinContent(i)
            bin_error = hist.GetBinError(i)
            hist_rescaled.SetBinContent(i, bin_content)
            hist_rescaled.SetBinError(i, bin_error)

    elif dimension == 2:
        # Get the old axis range
        x_min = hist.GetXaxis().GetXmin()
        x_max = hist.GetXaxis().GetXmax()
        x_bins = hist.GetNbinsX()
        y_min = hist.GetYaxis().GetXmin()
        y_max = hist.GetYaxis().GetXmax()
        y_bins = hist.GetNbinsY()
        if axis == 0:
            x_min *= scale
            x_max *= scale
        elif axis == 1:
            y_min *= scale
            y_max *= scale
        hist_rescaled = rt.TH2F(name, name, x_bins, x_min, x_max, y_bins, y_min, y_max)
        for i in range(1, x_bins + 1):
            for j in range(1, y_bins + 1):
                bin_content = hist.GetBinContent(i, j)
                bin_error = hist.GetBinError(i, j)
                hist_rescaled.SetBinContent(i, j, bin_content)
                hist_rescaled.SetBinError(i, j, bin_error)

    elif dimension == 3:
        # Get the old axis range
        x_min = hist.GetXaxis().GetXmin()
        x_max = hist.GetXaxis().GetXmax()
        x_bins = hist.GetNbinsX()
        y_min = hist.GetYaxis().GetXmin()
        y_max = hist.GetYaxis().GetXmax()
        y_bins = hist.GetNbinsY()
        z_min = hist.GetZaxis().GetXmin()
        z_max = hist.GetZaxis().GetXmax()
        z_bins = hist.GetNbinsZ()
        if axis == 0:
            x_min *= scale
            x_max *= scale
        elif axis == 1:
            y_min *= scale
            y_max *= scale
        elif axis == 2:
            z_min *= scale
            z_max *= scale
        hist_rescaled = rt.TH3F(
            name, name, x_bins, x_min, x_max, y_bins, y_min, y_max, z_bins, z_min, z_max
        )
        for i in range(1, x_bins + 1):
            for j in range(1, y_bins + 1):
                for k in range(1, z_bins + 1):
                    bin_content = hist.GetBinContent(i, j, k)
                    bin_error = hist.GetBinError(i, j, k)
                    hist_rescaled.SetBinContent(i, j, k, bin_content)
                    hist_rescaled.SetBinError(i, j, k, bin_error)

    # keep some meta data
    hist_rescaled.SetTitle(hist.GetTitle())
    if hist.GetSumw2N() == 0:
        hist_rescaled.Sumw2()

    return hist_rescaled


def NormalizeHistogramIntegral(hist, integral=1, name="normalized"):
    """
    Normalized integral of histogram to a given value
    Args:
        hist (TH1): histogram to be normalized
        integral (flaot): Value of the integral, by default equal to 1
        name (string): name of the normalized histogram

    Returns:
        Scaled histogram
    """
    hist_clone = hist.Clone(name)
    hist_integral = hist.Integral()
    if hist_integral != 0:
        hist_clone.Scale(integral / hist_integral)
    return hist_clone
