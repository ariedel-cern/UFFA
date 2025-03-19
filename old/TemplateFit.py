import ROOT
import FemtoAnalysis as FA
import array as arr

class ftotal():
    def __init__(self, data, mcdata):
        self.data = data
        self.histos = mcdata
        self.axis = data.GetXaxis()
    def __call__(self, arr, par):
        total = 0
        nbin = self.axis.FindBin(arr[0])
        for n in range(len(self.histos)):
            total += par[n]*self.histos[n].GetBinContent(nbin)

        return total

def find_bin_reduce_on_lower_edge(axis, value):
    found_bin = axis.FindBin(value)
    if value == axis.GetBinLowEdge(found_bin):
        found_bin -= 1
    return found_bin

def TemplateFit(fname, dca_data, dca_templates, dcacpa, dca_names, fit_range, pt_bins, pt_rebin, dirOut, temp_init, temp_limits, temp_fraction):
    xAxis = dca_data.GetXaxis()
    yAxis = dca_data.GetYaxis()

    # Output file
    ofile = fname
    if type(fname) == str:
        ofile = ROOT.TFile(dirOut + "TemplateFit_" + fname, "recreate")

    # initialize list of values if single rebin value was provided
    if pt_rebin and (len(pt_rebin) < len(pt_bins)):
        pt_rebin = pt_rebin*len(pt_bins)

    # setup for the bin ranges for the projections
    # output: [(bin1, bin2), (bin2, bin3), ...]
    if type(pt_bins) == int:            # int that describes how many bins to individually fit
        pt_ent = pt_bins
        pt_range = []
        for n in range(pt_ent):
            pt_range.append((n + 1, n + 1))
    elif type(pt_bins) == list:         # list of pt edges
        pt_ent = len(pt_bins) - 1
        pt_range = []
        for i in range(pt_ent):
            bin_value1 = xAxis.FindBin(pt_bins[i])
            bin_value2 = find_bin_reduce_on_lower_edge(xAxis, pt_bins[i + 1])
            pt_range.append((bin_value1, bin_value2))
    else:
        print("TemplateFit.py: pt_bins not an int or list of ranges: " + pt_bins)
        exit()

    # setup for the fraction fitting parameters
    if type(temp_fraction) == dict:
        temp_fraction = [temp_fraction]

    # titles of the graphs
    pt_names = []
    for n in range(pt_ent):
        #name = "p_{T} range: [%.3f-%.3f) GeV" % (xAxis.GetBinLowEdge(pt_range[n][0]), xAxis.GetBinLowEdge(pt_range[n][1]))
        name = "p_{T} range: " + f"[{pt_bins[n]}-{pt_bins[n + 1]}) GeV"
        if pt_rebin:
            name += " rebin: " + str(pt_rebin[n])
        pt_names.append(name)

    dca_ent = len(dca_templates)
    if type(dcacpa) == str:
        if dcacpa.lower() == "dca":
            dcacpa = 'dca'
        elif dcacpa.lower() == "cpa":
            dcacpa = 'cpa'
        else:
            print("Error in plot type: \"dca\" or \"cpa\"")

    # dca cut
    CPAcut = 0.99
    fitmax = []
    fitmin = []
    if dcacpa == 'dca':
        if type(fit_range) == list:
            for value in fit_range:
                fitmax.append(value)
                fitmin.append(-value)
            if len(fitmax) < dca_ent:
                fitmax += fitmax[-1]*(pt_ent - len(fitmax))    # append the last fitrange to be used for the rest of the pt bins
        else:
            fitmax = [fit_range]*pt_ent
            fitmin = [-fit_range]*pt_ent
    else:
        fitmin = 0.9
        fitmax = 1.0

    # graph initialization
    gChi = ROOT.TGraph(pt_ent - 1)
    m2ent = [0]*dca_ent
    dataEntries = []
    mcEntries = []
    parDCA_mc = []
    gDCA_mc = []
    for n in range(dca_ent):
        parDCA_mc.append([])
        gDCA_mc.append(ROOT.TGraph(pt_ent - 1))
        gDCA_mc[n].SetName(dca_names[n])
        gDCA_mc[n].SetTitle(dca_names[n])
        gDCA_mc[n].SetLineWidth(2)
        gDCA_mc[n].SetLineColor(n + 3)
        gDCA_mc[n].SetMarkerStyle(21)
        gDCA_mc[n].SetMarkerSize(1)
        gDCA_mc[n].SetMarkerColor(n + 3)
        gDCA_mc[n].GetXaxis().SetLabelSize(0.05)
        gDCA_mc[n].GetYaxis().SetLabelSize(0.05)
        gDCA_mc[n].GetXaxis().SetTitleSize(0.05)
        gDCA_mc[n].GetXaxis().SetTitle("<p_{T}> (GeV)")

    # main loop for fitting
    for n in range(pt_ent):
        # canvas for fitting
        canvas = ROOT.TCanvas("canvas_" + str(n + 1), "canvas_" + str(n + 1))
        canvas.SetCanvasSize(1024, 768)
        canvas.cd()
        ROOT.gPad.SetLogy()

        # data
        data = dca_data.ProjectionY("hDCAxy_" + str(n + 1), pt_range[n][0], pt_range[n][1])
        if pt_rebin and pt_rebin[n] != 1:
            data.Rebin(pt_rebin[n])
        data.SetAxisRange(fitmin[n], fitmax[n])
        if dcacpa == "cpa":
            data_int = data.Integral(data.FindBin(CPAcut), data.FindBin(1))
        else:
            data_int = data.Integral(data.FindBin(fitmin[n]), data.FindBin(fitmax[n]))
        dataEntries.append(data_int)
        if data.Integral():
            data.Scale(1. / data.Integral())

        # MC templates
        hDCA_mc = []
        for i in range(dca_ent):
            hDCA_mc.append(dca_templates[i].ProjectionY(dca_names[i] + '_' + str(n + 1), pt_range[n][0], pt_range[n][1]))
            if pt_rebin and pt_rebin[n] != 1:
                hDCA_mc[i].Rebin(pt_rebin[n])
            hDCA_mc[i].SetAxisRange(fitmin[n], fitmax[n])
            hDCA_mc[i].SetTitle(pt_names[n])
            if hDCA_mc[i].Integral():
                hDCA_mc[i].Scale(1. / hDCA_mc[i].Integral())

        # fit function
        adj = ftotal(data, hDCA_mc)
        ftot = ROOT.TF1("ftot", adj, fitmin[n], fitmax[n], dca_ent)

        for i in range(dca_ent):
            ftot.SetParameter(i, 1. / dca_ent)
            ftot.SetParLimits(i, 0., 1.)

            if temp_init:
                ftot.SetParameter(i, temp_init[i])
            if temp_limits:
                ftot.SetParLimits(i, temp_limits[i][0], temp_limits[i][1])
            if temp_fraction:
                for entry in temp_fraction:
                    if entry['temp_name'] == dca_names[i]:
                        if entry['temp_init']:
                            ftot.SetParameter(i, entry['temp_init'][n])
                            ftot.SetParLimits(i, entry['temp_init'][n], entry['temp_init'][n])
                        if entry['temp_limits']:
                            if type(entry['temp_limits'][n]) == list:
                                ftot.SetParLimits(i, entry['temp_limits'][n][0], entry['temp_limits'][n][1])
                            else:
                                if entry['temp_limits'][n] == 0 and not entry['temp_init']:
                                    if not entry['temp_init']:
                                        ftot.FixParameter(i, 0)

        # draw histograms
        data.Fit("ftot", "S, N, R, M", "", fitmin[n], fitmax[n])
        data.SetLineColor(1)
        data.SetLineWidth(2)
        data.Draw("hist")

        for i in range(dca_ent):
            hDCA_mc[i].Scale(ftot.GetParameter(i))
            hDCA_mc[i].SetLineColor(i + 3)
            hDCA_mc[i].SetLineWidth(2)
            hDCA_mc[i].Draw("same")

        htot = hDCA_mc[0].Clone("htot_" + str(n + 1))
        for i in range(1, dca_ent):
            htot.Add(hDCA_mc[i])
        htot.SetLineColor(2)
        htot.SetLineWidth(2)
        htot.Draw("same")

        # legend
        legend = ROOT.TLegend(0.15, 0.40, 0.35, 0.8)
        legend.SetTextSize(0.05)
        legend.SetBorderSize(0)
        legend.SetLineColor(1)
        legend.SetLineWidth(1)
        legend.SetFillColor(0)
        legend.SetFillStyle(1001)

        legend.AddEntry(data, "total", "l")
        legend.AddEntry(htot, "fit", "l")
        for i in range(dca_ent):
            legend.AddEntry(hDCA_mc[i], dca_names[i], "l")
        legend.Draw()

        # write to file
        ofile.cd()
        for i in range(dca_ent):
            hDCA_mc[i].Write()
        canvas.Write()

        # fractions
        if dcacpa == "cpa":
            htot_int = htot.Integral(htot.FindBin(CPAcut), htot.FindBin(1))
        else:
            htot_int = htot.Integral(htot.FindBin(fitmin[n]), htot.FindBin(fitmax[n]))
        mcEntries.append(htot_int)
        if dcacpa == "cpa":
            for i in range(dca_ent):
                parDCA_mc[i].append(hDCA_mc[i].Integral(hDCA_mc[i].FindBin(CPAcut), hDCA_mc[i].FindBin(1)) / mcEntries[n])
        else:
            for i in range(dca_ent):
                tmp = hDCA_mc[i].Integral(hDCA_mc[i].FindBin(fitmin[n]), hDCA_mc[i].FindBin(fitmax[n]))
                if mcEntries[n]:
                    parDCA_mc[i].append(tmp / mcEntries[n])
                else:
                    parDCA_mc[i].append(tmp)

        # fill x and y lists for graphs
        pt_avg = (xAxis.GetBinLowEdge(pt_range[n][0]) + xAxis.GetBinLowEdge(pt_range[n][1])) / 2
        for i in range(dca_ent):
            gDCA_mc[i].SetPoint(n, pt_avg, parDCA_mc[i][n])
        gChi.SetPoint(n, pt_avg, ftot.GetChisquare() / (ftot.GetNDF() + dca_ent))

        canvas.Close()
        del canvas

    # canvas for chi^2
    c2 = ROOT.TCanvas("c2", "c2", 0, 0, 1024, 768)
    c2.cd()
    gChi.SetLineColor(1)
    gChi.Draw()

    # calculate final pT weighted result
    m2tot = 0
    for i in range(pt_ent):
        m2tot += dataEntries[i]
        for j in range(dca_ent):
            m2ent[j] += dataEntries[i]*parDCA_mc[j][i]

    data_bins = dca_data.GetXaxis().GetNbins()
    BinRanges_pT = arr.array('d', [0]*(data_bins + 1))
    BinRanges_pT[data_bins] = dca_data.GetXaxis().GetBinUpEdge(data_bins)
    pT_Weights = ROOT.TH1F("pT_Weights", "pT_Weights", pt_ent, BinRanges_pT)
    for i in range(pt_ent - 1):
        pT_Weights.SetBinContent(i, dataEntries[i] / m2tot)

    # chi2 graph
    gChi.SetLineWidth(2)
    gChi.GetXaxis().SetLabelSize(0.05)
    gChi.GetYaxis().SetLabelSize(0.05)
    gChi.GetYaxis().SetTitle("chi2/ndf")
    gChi.GetXaxis().SetTitle("<p_{T}> (GeV)")
    gChi.GetXaxis().SetTitleSize(0.05)
    gChi.GetYaxis().SetTitleSize(0.05)

    # canvas for all fractions
    fractions = ROOT.TCanvas("fractions", "fractions", 1024, 768)
    ROOT.gPad.SetGridx()
    ROOT.gPad.SetGridy()
    fractions.cd()
    fractions.SetFillColor(0)
    fractions.SetBorderMode(0)
    fractions.SetBorderSize(2)
    fractions.SetFrameBorderMode(0)
    fractions.SetMargin(0.15, 0.05, 0.2, 0.05)

    # fraction graph
    gEmpty = gDCA_mc[0]
    gEmpty.SetTitle("")
    gEmpty.SetFillStyle(1000)
    gEmpty.SetMinimum(0.)
    gEmpty.SetMaximum(1.0)
    gx = gEmpty.GetXaxis()
    gy = gEmpty.GetYaxis()
    gx.SetLimits(xAxis.GetBinLowEdge(pt_range[0][0] - 1), xAxis.GetBinLowEdge(pt_range[-1][1] + 1))
    gx.SetTitle("<p_{T}> (GeV)")
    gx.SetLabelFont(42)
    gx.SetLabelSize(0.06)
    gx.SetTitleSize(0.07)
    gx.SetTitleOffset(1)
    gx.SetTitleFont(42)
    gy.SetTitle("fractions")
    gy.SetLabelFont(42)
    gy.SetLabelSize(0.06)
    gy.SetTitleSize(0.07)
    gy.SetTitleFont(42)
    gEmpty.Draw("alp")
    for i in range(1, dca_ent):
        gDCA_mc[i].Draw("lp same")
    for i in range(pt_ent):
        line = ROOT.TLine(xAxis.GetBinLowEdge(pt_range[i][0]), parDCA_mc[0][i], xAxis.GetBinLowEdge(pt_range[i][1]), parDCA_mc[0][i])
        line.DrawClone("same")

    # primary average
    primAvg = ROOT.TF1("primAvg", "[0]", 0, 5)
    primAvg.SetParameter(0, m2ent[0] / m2tot)
    primAvg.SetLineColor(12)
    primAvg.SetLineStyle(2)
    primAvg.Draw("same")

    for i in range(dca_ent):
        print("%-12s %.6f" % (dca_names[i], m2ent[i] / m2tot))

    ofile.cd()
    for i in range(dca_ent):
        gDCA_mc[i].Write()
    gChi.Write("chi2/ndf")
    fractions.Write("all")
    primAvg.Write()
    pT_Weights.Write()

    c2.Close()
    fractions.Close()
    del c2, fractions, htot, ftot, adj

