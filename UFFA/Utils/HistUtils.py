import ROOT as rt


def FindBinWithUpperEdgeDetection(Axis, Value):
    """
    Find bin of a TAxis for a given value with edge detection
    Args:
        Axis (TAxis): Axis with bins
        Value (float): Value whose bin is being searched

    Returns:

    """
    FoundBin = Axis.FindBin(Value)
    if Value == Axis.GetBinLowEdge(FoundBin) and FoundBin != 1:
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
