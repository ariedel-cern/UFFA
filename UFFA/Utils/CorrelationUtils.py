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
        Proj = hist_2d.ProjectionX(name, 0, -1, "e")
    elif axis == 1:
        Proj = hist_2d.ProjectionY(name, 0, -1, "e")
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
    MeRw.Reset()

    for ybin in range(1, se.GetNbinsY()):
        SeBin = se.ProjectionX(f"se_bin_{ybin}", ybin, ybin)
        MeBin = me.ProjectionX(f"me_bin_{ybin}", ybin, ybin)
        SeInt = SeBin.Integral()
        MeInt = MeBin.Integral()
        if MeInt > 0.0 and SeInt > 0.0:
            MeBin.Scale(SeInt / MeInt)
            for xbin in range(1, MeBin.GetNbinsX() + 1):
                MeRw.SetBinContent(xbin, ybin, MeBin.GetBinContent(xbin))
                MeRw.SetBinError(xbin, ybin, MeBin.GetBinError(xbin))
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


def Recenter(cf, me, skipEmptyBins=False):
    """
    Recenter bins of the correlation function according to mixed event distribution with a finer binning
    Args:
        cf (TH1): correlation function to be recentered
        me (TH1): fine binned mixed event distribution
        skipEmptyBins (bool): skip empty bins, so there is not data point at 0

    Returns: Rescaled correlation functions as TGraphErrors
    """

    # get rebin factor
    rebin = int(me.GetNbinsX() / cf.GetNbinsX())
    Nbins = cf.GetNbinsX()

    BinCenterX = []
    BinErrorX = []
    BinCenterY = []
    BinErrorY = []

    # loop ever all bins in the correlation function
    for bin in range(1, Nbins + 1):

        value = 0
        norm = 0
        error = 0

        # loop over underlying bins in fine binned mixed event distribution
        for bins in range(1 + (bin - 1) * rebin, bin * rebin + 1):
            value = value + me.GetBinCenter(bins) * me.GetBinContent(bins)
            norm = norm + me.GetBinContent(bins)

        # if requested, we skip empty bins
        if skipEmptyBins == True and norm == 0:
            continue

        if norm > 0:
            value = value / norm
        else:
            # if we do not skip empty bins, then use the original bin center
            value = cf.GetBinCenter(bin)

        for bins in range(1 + (bin - 1) * rebin, bin * rebin + 1):
            error = error + me.GetBinContent(bins) * np.power(
                me.GetBinCenter(bins) - value, 2
            )

        BinCenterX.append(value)
        if norm > 0:
            BinErrorX.append(np.sqrt(error / norm))
        else:
            BinErrorX.append(cf.GetBinWidth(bin) / 2.0)
        BinCenterY.append(cf.GetBinContent(bin))
        BinErrorY.append(cf.GetBinError(bin))

    BinCenterX = np.array(BinCenterX, dtype=ct.c_double)
    BinErrorX = np.array(BinErrorX, dtype=ct.c_double)
    BinCenterY = np.array(BinCenterY, dtype=ct.c_double)
    BinErrorY = np.array(BinErrorY, dtype=ct.c_double)

    return rt.TGraphErrors(
        len(BinCenterX), BinCenterX, BinCenterY, BinErrorX, BinErrorY
    )


def RescaleGraph(graph, scale, suffix="_Rescaled"):
    """
    Rescale x axis of TGraphErrors
    """
    BinCenterX = []
    BinErrorX = []
    BinCenterY = []
    BinErrorY = []
    n_bins = graph.GetN()

    for i in range(n_bins):
        BinCenterX.append(graph.GetPointX(i) * scale)
        BinErrorX.append(graph.GetErrorX(i) * scale)
        BinCenterY.append(graph.GetPointY(i))
        BinErrorY.append(graph.GetErrorY(i))
    BinCenterX = np.array(BinCenterX, dtype=ct.c_double)
    BinErrorX = np.array(BinErrorX, dtype=ct.c_double)
    BinCenterY = np.array(BinCenterY, dtype=ct.c_double)
    BinErrorY = np.array(BinErrorY, dtype=ct.c_double)
    g = rt.TGraphErrors(n_bins, BinCenterX, BinCenterY, BinErrorX, BinErrorY)
    g.SetNameTitle(f"{graph.GetName()}{suffix}", f"{graph.GetTitle()}{suffix}")
    return g


def GetRelativeUncertainties(object):
    """
    Return relative uncertainties bin by bin or point by point
    """
    rel_unc = []
    if object.InheritsFrom(rt.TH1.Class()):
        for i in range(object.GetNbinsX()):
            y = object.GetBinContent(i + 1)
            ey = object.GetBinError(i + 1)

            if y != 0:
                rel = ey / y
            else:
                rel = 0
            rel_unc.append(rel)
        return np.array(rel_unc)
    if object.InheritsFrom(rt.TGraph.Class()):
        for i in range(object.GetN()):
            y = object.GetPointY(i)
            ey = object.GetErrorY(i)
            if y != 0:
                rel = ey / y
            else:
                rel = 0
            rel_unc.append(rel)
        return np.array(rel_unc)


def SetMinimalUncertainty(graph, threshold):
    """
    Adjusts the uncertainties in the TGraph if they are lower than the threshold.

    Parameters:
    graph (ROOT.TGraph): The TGraph object to process.
    threshold (float): The minimum value for the uncertainties.
    """

    g = graph.Clone()

    # Get the number of points in the graph
    n_points = g.GetN()

    if isinstance(g, rt.TGraphAsymmErrors):
        # Get the number of points in the graph
        n_points = g.GetN()

        # Loop over each point
        for i in range(n_points):
            # Get the current Y error (uncertainty) at the i-th point
            y_error_low = g.GetErrorYlow(i)
            y_error_up = g.GetErrorYhigh(i)
            # If the current error is lower than the threshold, adjust it
            if y_error_low < threshold:
                g.SetPointError(
                    i, g.GetErrorXlow(i), g.GetErrorXhigh(i), threshold, y_error_up
                )
            if y_error_up < threshold:
                g.SetPointError(
                    i, g.GetErrorXlow(i), g.GetErrorXhigh(i), y_error_low, threshold
                )

    elif isinstance(g, rt.TGraph):
        # For TGraph (symmetric errors), adjust only the Y error
        n_points = g.GetN()

        # Loop over each point
        for i in range(n_points):
            # Get the current Y error (uncertainty) at the i-th point
            y_error = g.GetErrorY(i)

            # If the current error is lower than the threshold, adjust it
            if y_error < threshold:
                g.SetPointError(
                    i, g.GetErrorX(i), threshold
                )  # Adjust Y error to threshold

    else:
        raise TypeError("Input graph must be a TGraph or TGraphAsymmErrors.")

    return g


def GetNsigma(stat, syst, theory):
    """
    Computes the deviation between the data and the theory in terms of n-sigma using interpolation,
    and returns the result as a TGraphErrors.

    Parameters:
    stat_graph (ROOT.TGraphErrors): Graph with data points and statistical uncertainties.
    sys_graph (ROOT.TGraphErrors): Graph with data points and systematic uncertainties.
    theory_graph (ROOT.TGraphErrors): Graph with theory predictions (may have more points).

    Returns:
    ROOT.TGraphErrors: A TGraphErrors object containing the n-sigma deviations for each data point.
    """
    # Extract points from the data graph (x and y) and their uncertainties
    n_data = stat.GetN()
    stat_x = np.array([stat.GetX()[i] for i in range(n_data)])
    stat_y = np.array([stat.GetY()[i] for i in range(n_data)])
    stat_error = np.array([stat.GetErrorY(i) for i in range(n_data)])

    # Extract systematic uncertainties from the sys graph
    sys_error = np.array([syst.GetErrorY(i) for i in range(n_data)])

    # Calculate total error for the data points
    total_error = np.sqrt(stat_error**2 + sys_error**2)

    # Extract theory points (x and y) and their uncertainties
    n_theory = theory.GetN()
    theory_x = np.array([theory.GetX()[i] for i in range(n_theory)])
    theory_y = np.array([theory.GetY()[i] for i in range(n_theory)])
    theory_error = np.array([theory.GetErrorY(i) for i in range(n_theory)])

    # Create a list to store the deviations (n-sigma)
    n_sigma_deviations = []

    # Loop over the data points and compare with the closest theory point
    for i in range(n_data):
        # Find the nearest theory point by minimizing the difference in x-values
        min_diff = np.abs(theory_x - stat_x[i])
        nearest_theory_idx = np.argmin(min_diff)

        # Compute the absolute difference between the data and theory at the (approximately) same x
        y_diff = stat_y[i] - theory_y[nearest_theory_idx]

        # Compute the combined error (data + theory error) for the n-sigma
        combined_error = np.sqrt(
            total_error[i] ** 2 + theory_error[nearest_theory_idx] ** 2
        )

        # Compute the deviation in n-sigma
        n_sigma = y_diff / combined_error
        n_sigma_deviations.append(n_sigma)

    # Create the  object to store the n-sigma deviations
    n_sigma_graph = rt.TGraph(n_data)

    # Fill the TGraphErrors with the n-sigma values
    for i in range(n_data):
        n_sigma_graph.SetPoint(i, stat_x[i], n_sigma_deviations[i])

    return n_sigma_graph
