import ROOT as rt
import numpy as np
import ctypes as ct


def Proj2dTo1d(hist_2d, axis, name):
    """
    Project a TH2 to a TH1 histogram

    Args:
        hist_2d (TH2): 2d histogram to project down into 1d
        axis (int): Index of the axis to project down (0->x, 1->y)
        name (string): Name of the projected histogram

    Returns:
       Projected TH1 Histogram
    """
    Proj = None
    # use option "e" to trigger calculation of errors
    if axis == 0:
        Proj = hist_2d.ProjectionX(name, option="e")
    elif axis == 1:
        Proj = hist_2d.ProjectionY(name, option="e")
    return Proj


def Proj3dTo2d(hist_3d, axis_x, axis_y, name):
    """
    Project a TH3 to a TH2 histogram

    Args:
        Hist3d (TH3): 3d histogram to project down into 2d
        AxisIndexX (int): Index of the axis to project into x-axis
        AxisIndexY (int): Index of the axis to project into y-axis
        Name (string): Name of projected histogram

    Example: If AxisIndexX=2 and AxisIndexY=0 is provided, then
    the z-axis (AxisIndexX=2) of the 3d histogram is projected onto the x-axis and
    the x-axis (AxisIndexY=0) of the 3d histogram is projected onto the y-axis
    of the projected 2d histogram

    Returns:
       Projected TH2 Histogram
    """
    Proj = None
    proj = ""
    # ROOT convention, the y axis has to be specified first
    if axis_y == 0:
        proj = proj + "x"
    elif axis_y == 1:
        proj = proj + "y"
    elif axis_y == 2:
        proj = proj + "z"
    if axis_x == 0:
        proj = proj + "x"
    elif axis_x == 1:
        proj = proj + "y"
    elif axis_x == 2:
        proj = proj + "z"
    Proj = hist_3d.Project3D(proj)
    Proj.SetName(name)
    return Proj


def Proj3dTo1d(hist_3d, axis, name):
    """
    Project TH3 to a TH1 histogram
    Args:
        Hist3d (TH3): 3d histogram to project down into 2d
        AxisIndex (int): Index of axis to project down on
        Name (string): Name of the new histogram

    Returns:
       Projected TH1 histogram
    """
    Proj = None
    proj = ""
    if axis == 0:
        proj = "x"
    elif axis == 1:
        proj = "y"
    elif axis == 2:
        proj = "z"
    Proj = hist_3d.Project3D(proj)
    Proj.SetName(name)
    return Proj


def ProjNdTo2d(hist_Nd, axis_x, axis_y, name):
    """
    Project THn histogram onto a TH2 histogram
    Args:
        HistNd (THn): nd histogram to project down into 2d
        AxisIndexX (int): Index of the axis to project into x-axis
        AxisIndexY (int): Index of the axis to project into y-axis
        Name (string): Name of projected histogram
    Returns:
       Projected TH2 Histogram
    """
    Proj = hist_Nd.Projection(axis_y, axis_x, "e")
    Proj.SetName(name)
    return Proj


def ProjNdTo1d(hist_Nd, axis, name):
    """
    Project THn histogram onto a TH1 histogram
    Args:
        HistNd (THn): nd histogram to project down into 2d
        AxisIndex (int): Index of axis to project down on
        Name (string): Name of projected histogram
    Returns:
       Projected TH1 Histogram
    """
    Proj = hist_Nd.Projection(axis, "e")
    Proj.SetName(name)
    return Proj


def Proj2dTo2d(hist_2d, axis_x, axis_y, name):
    """
    Project TH2 histogram onto a TH2 histogram
    This is a convient wrapper to swap x and y axis of a TH2 histogram
    If the order of the axis are as desired, just return a clone of the original histogram with a new name
    Args:
        Hist2d (THn): 2d histogram to project down into 2d
        AxisIndexX (int): Index of the axis to project into x-axis
        AxisIndexY (int): Index of the axis to project into y-axis
        Name (string): Name of projected histogram
    Returns:
       Projected TH2 Histogram
    """
    Hist = None
    if axis_x == 0 and axis_y == 1:
        Hist = hist_2d.Clone(name)
    else:
        Hist = rt.TH2D(
            name,
            name,
            hist_2d.GetNbinsY(),
            hist_2d.GetYaxis().GetXmin(),
            hist_2d.GetYaxis().GetXmax(),
            hist_2d.GetNbinsX(),
            hist_2d.GetXaxis().GetXmin(),
            hist_2d.GetXaxis().GetXmax(),
        )
        for i in range(hist_2d.GetNbinsX()):
            for j in range(hist_2d.GetNbinsY()):
                Hist.SetBinContent(j, i, hist_2d.GetBinContent(i, j))
    return Hist


def Reweight2d(se, me, name):
    """
    Args:
        Se (TH2): Same event distribution (kstar vs reweight axis)
        Me (TH2): Mixed event distribution (kstar vs reweight axis)
        Name (string): Name of reweighted mixed event distribution

    Returns:
        Reweighted mixed event distribtion
    """
    MeRw = me.Clone(name)
    MeRw.Reset("ICESM")

    for ybin in range(1, se.GetNbinsY()):
        SeBin = se.ProjectionX(f"se_bin_{ybin}", ybin, ybin)
        MeBin = me.ProjectionX(f"me_bin_{ybin}", ybin, ybin)
        SeInt = SeBin.Integral()
        MeInt = MeBin.Integral()
        if MeInt > 0.0 and SeInt > 0.0:
            MeBin.Scale(SeInt / MeInt)
            for xbin in range(1, MeBin.GetNbinsX() + 1):
                MeRw.SetBinContent(xbin, ybin, MeBin.GetBinContent(xbin))
    return MeRw


def Normalize(distribution, range):
    """
    Normalize distribution in given range
    Args:
        Dist (TH1): 1d histogram
        range (tuple): Tuple with lower and upper bound

    Returns:
       Normalized histogram
    """
    LowerBound = range[0]
    UpperBound = range[1]
    LowerBin = distribution.FindBin(LowerBound)
    LowerEdge = distribution.GetBinLowEdge(LowerBin)
    UpperBin = distribution.FindBin(UpperBound)
    UpperEdge = distribution.GetBinLowEdge(UpperBin + 1)
    Dist_Integral = distribution.Integral(LowerBin, UpperBin, "width")
    try:
        distribution.Scale((UpperEdge - LowerEdge) / Dist_Integral)
    except:
        print("Integral in normalization range is zero")
    return distribution


def RecenterCorrelationFunction(cf, me):
    """
    Recenter bins of the correlation function according to mixed event distribution
    """

    rebin = int(me.GetNbinsX() / cf.GetNbinsX())
    Nbins = cf.GetNbinsX()

    BinCenterX = []
    BinErrorX = []
    BinCenterY = []
    BinErrorY = []

    for bin in range(1, Nbins + 1):
        value = 0
        norm = 0
        error = 0
        for bins in range(1 + (bin - 1) * rebin, bin * rebin):
            value = value + me.GetBinCenter(bins) * me.GetBinContent(bins)
            norm = norm + me.GetBinContent(bins)
        value = value / norm
        for bins in range(1 + (bin - 1) * rebin, bin * rebin):
            error = error + me.GetBinContent(bins) * np.power(
                me.GetBinCenter(bins) - value, 2
            )
        BinCenterX.append(value)
        BinErrorX.append(np.sqrt(error / norm))
        BinCenterY.append(cf.GetBinContent(bin))
        BinErrorY.append(cf.GetBinError(bin))

    BinCenterX = np.array(BinCenterX, dtype=ct.c_double)
    BinErrorX = np.array(BinErrorX, dtype=ct.c_double)
    BinCenterY = np.array(BinCenterY, dtype=ct.c_double)
    BinErrorY = np.array(BinErrorY, dtype=ct.c_double)

    return rt.TGraphErrors(Nbins, BinCenterX, BinCenterY, BinErrorX, BinErrorY)
