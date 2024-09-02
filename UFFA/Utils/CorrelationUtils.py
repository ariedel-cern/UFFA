import ROOT as rt


def Proj2dTo1d(Hist2d, AxisIndex, Name):
    """
    Project a TH2 to a TH1 histogram

    Args:
        Hist2d (TH2): 2d histogram to project down into 1d
        AxisIndex (int): Index of the axis to project down (0->x, 1->y)
        Name (string): Name of the projected histogram

    Returns:
       Projected TH1 Histogram
    """
    Proj = None
    # use option "e" to trigger calculation of errors
    if AxisIndex == 0:
        Proj = Hist2d.ProjectionX(Name, option="e")
    elif AxisIndex == 1:
        Proj = Hist2d.ProjectionY(Name, option="e")
    return Proj


def Proj3dTo2d(Hist3d, AxisIndexX, AxisIndexY, Name):
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
    if AxisIndexY == 0:
        proj = proj + "x"
    elif AxisIndexY == 1:
        proj = proj + "y"
    elif AxisIndexY == 2:
        proj = proj + "z"
    if AxisIndexX == 0:
        proj = proj + "x"
    elif AxisIndexX == 1:
        proj = proj + "y"
    elif AxisIndexX == 2:
        proj = proj + "z"
    Proj = Hist3d.Project3D(proj)
    Proj.SetName(Name)
    return Proj


def Proj3dTo1d(Hist3d, AxisIndex, Name):
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
    if AxisIndex == 0:
        proj = "x"
    elif AxisIndex == 1:
        proj = "y"
    elif AxisIndex == 2:
        proj = "z"
    Proj = Hist3d.Project3D(proj)
    Proj.SetName(Name)
    return Proj


def ProjNdTo2d(HistNd, AxisIndexX, AxisIndexY, Name):
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
    Proj = HistNd.Projection(AxisIndexY, AxisIndexX, "e")
    Proj.SetName(Name)
    return Proj


def ProjNdTo1d(HistNd, AxisIndex, Name):
    """
    Project THn histogram onto a TH1 histogram
    Args:
        HistNd (THn): nd histogram to project down into 2d
        AxisIndex (int): Index of axis to project down on
        Name (string): Name of projected histogram
    Returns:
       Projected TH1 Histogram
    """
    Proj = HistNd.Projection(AxisIndex, "e")
    Proj.SetName(Name)
    return Proj


def Proj2dTo2d(Hist2d, AxisIndexX, AxisIndexY, Name):
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
    if AxisIndexX == 0 and AxisIndexY == 1:
        Hist = Hist2d.Clone(Name)
    else:
        Hist = rt.TH2D(
            Name,
            Name,
            Hist2d.GetNbinsY(),
            Hist2d.GetYaxis().GetXmin(),
            Hist2d.GetYaxis().GetXmax(),
            Hist2d.GetNbinsX(),
            Hist2d.GetXaxis().GetXmin(),
            Hist2d.GetXaxis().GetXmax(),
        )
        for i in range(Hist2d.GetNbinsX()):
            for j in range(Hist2d.GetNbinsY()):
                Hist.SetBinContent(j, i, Hist2d.GetBinContent(i, j))
    return Hist


def Reweight2d(Se, Me, Name):
    """
    Args:
        Se (TH2): Same event distribution (kstar vs reweight axis)
        Me (TH2): Mixed event distribution (kstar vs reweight axis)
        Name (string): Name of reweighted mixed event distribution

    Returns:
        Reweighted mixed event distribtion
    """
    MeRw = Me.Clone(Name)
    MeRw.Reset("ICESM")

    for ybin in range(1, Se.GetNbinsY()):
        SeBin = Se.ProjectionX(f"se_bin_{ybin}", ybin, ybin)
        MeBin = Me.ProjectionX(f"me_bin_{ybin}", ybin, ybin)
        SeInt = SeBin.Integral()
        MeInt = MeBin.Integral()
        if MeInt > 0.0 and SeInt > 0.0:
            MeBin.Scale(SeInt / MeInt)
            for xbin in range(1, MeBin.GetNbinsX() + 1):
                MeRw.SetBinContent(xbin, ybin, MeBin.GetBinContent(xbin))
    return MeRw


def Normalize(Dist, range):
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
    LowerBin = Dist.FindBin(LowerBound)
    LowerEdge = Dist.GetBinLowEdge(LowerBin)
    UpperBin = Dist.FindBin(UpperBound)
    UpperEdge = Dist.GetBinLowEdge(UpperBin + 1)
    Dist_Integral = Dist.Integral(LowerBin, UpperBin, "width")
    try:
        Dist.Scale((UpperEdge - LowerEdge) / Dist_Integral)
    except:
        print("Integral in normalization range is zero")
    return Dist
