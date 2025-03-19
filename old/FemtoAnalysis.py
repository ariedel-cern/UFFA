import ROOT
import FileUtils as FU
import FemtoDreamSaver as FDS
import FemtoDreamReader as FDR
import CorrelationHandler as CH
import CombinedTemplateFit as TF

def UFFA(settings):
    conf = config(settings)
    if conf['function'] == 'cf':
        UFFA_cf(conf)
    elif conf['function'] == 'tf':
        UFFA_tf(conf)
    elif conf['function'] == 'tf2d':
        UFFA_tf2d(conf)
    elif conf['function'] == 'ctf':
        UFFA_ctf(conf)
    elif conf['function'] == 'syst':
        if conf['htype'] in ['mtmult', 'rew3d', '4d', 'rew4d']:
            UFFA_syst_3d(conf)
        else:
            UFFA_syst(conf)

# correlation function
def UFFA_cf(settings):
    conf = config(settings)
    fdr = FDR.FemtoDreamReader(conf['fullpath'], conf['fileTDir'])
    ch = cf_handler(fdr, conf)
    fds = FDS.FemtoDreamSaver(conf, ch.get_histos())

# template fits
def UFFA_tf(settings):
    conf = config(settings)
    if conf['file']:
        fdr1 = FDR.FemtoDreamReader(conf['fullpath'], conf['fileTDir'])
        dca_data = fdr1.get_dca()
    elif conf['data']:
        dca_data = conf['data']
    else:
        print('UFFA_tf: Missing input data!')
    if conf['templates']:
        if type(conf['templates']) == str:
            fdr2 = FDR.FemtoDreamReader(conf['templates'], conf['mcTDir'])
            dca_mcplots = fdr2.get_dca_mc()
        else:
            dca_mcplots = conf['templates']
    else:
        dca_mcplots = fdr1.get_dca_mc()

    fds = FDS.FemtoDreamSaver(settings)
    ofile = fds.getFile()

    TF.TemplateFit(ofile, dca_data, dca_mcplots, conf['tftype'], conf['namelist'], conf['fitrange'], conf['signalrange'], conf['bins'], conf['rebin'], conf['outDir'], conf['temp_init'], conf['temp_limits'], conf['temp_fraction'], conf['print'])

# template fits 2d
def UFFA_tf2d(settings):
    conf = config(settings)
    dca_data = conf['data']
    dca_mcplots = conf['templates']

    fds = FDS.FemtoDreamSaver(settings)
    ofile = fds.getFile()

    TF.TemplateFit2D(ofile, dca_data, dca_mcplots, conf['namelist'], conf['fitrange'], conf['signalrange'], conf['bins'], conf['rebin'], conf['outDir'], conf['temp_init'], conf['temp_limits'], conf['temp_fraction'], conf['print'], conf['debug'])

# combined template fits
def UFFA_ctf(settings):
    conf = config(settings)
    if conf['file']:
        fdr1 = FDR.FemtoDreamReader(conf['fullpath'], conf['fileTDir'])
        dca_data = fdr1.get_dca()
    elif conf['data']:
        dca_data = conf['data']
    else:
        print('UFFA_tf: Missing input data!')
    if conf['templates']:
        if type(conf['templates']) == str:
            fdr2 = FDR.FemtoDreamReader(conf['templates'], conf['mcTDir'])
            dca_mcplots = fdr2.get_dca_mc()
        else:
            dca_mcplots = conf['templates']
    else:
        dca_mcplots = fdr1.get_dca_mc()

    fds = FDS.FemtoDreamSaver(settings)
    ofile = fds.getFile()

    TF.CombinedFit(ofile, conf['outDir'], dca_data, dca_mcplots, conf['namelist'], conf['fitrange'], conf['signalrange'], conf['bins'], conf['rebin'], conf['temp_init'], conf['temp_limits'], conf['temp_fraction'], conf['print'])

# systematics
def UFFA_syst(settings):
    conf = config(settings)
    fdr = FDR.FemtoDreamReader(conf['fullpath'], conf['fileTDir'])

    # default cf
    ch = cf_handler(fdr, conf)
    cf, cf_unw = ch.get_cf()                                # [[cf, [rebins]], [bin2...], ...], [[cf unw, [rebins]], [bin2...], ...]

    # input same event for yield filtering
    if conf['yield']:
        se = fdr.get_se()
        pair_num_se = se.Integral(se.FindBin(0), se.FindBin(conf['yield'][0]))
    if conf['debug']:
        se_all = ch.get_se()

    cf_list = []
    if conf['rebin']:
        len_rebin = len(conf['rebin'])

    if conf['atype'] == 'int':                              # integrated
        ck, ck_rebin = cf[0]
        cf_list.append([ck, ck_rebin])
        syst = [[Systematics(ck), []]]                          # [[syst cf, [rebins]]]
        if conf['rebin']:
            for i in range(len_rebin):
                syst[0][1].append(Systematics(ck_rebin[i]))
    elif conf['atype'] == 'dif':                            # differential
        syst = []
        for n, [ck, ck_rebin] in enumerate(cf):
            cf_list.append([ck, ck_rebin])
            syst.append([Systematics(ck), []])                  # [[syst cf, [rebins]], [bin2...], ...]
            if conf['rebin']:
                for i in range(len_rebin):
                    syst[n][1].append(Systematics(ck_rebin[i]))

    # loop over data variations in file and calculate the cf for each
    # which is then saved in a th2 from which the systematic error is computed and saved in a th1
    file_dir = fdr.get_dir();
    fdr.cd(0)           # class method of FileSaver to return to root of file
    folders = fdr.get_folder_names()
    for folder in folders:
        fdr.cd(folder)

        # allows to include/exclude specific variations
        if conf['exclude'] and folder in conf['exclude']:
            continue
        elif conf['include']:
            if folder in conf['include']:
                pass
            else:
                continue
        elif folder.rsplit('_')[-1][:3] != "Var":
            continue

        ch_var = cf_handler(fdr, conf)
        cf_var, cf_var_unw = ch_var.get_cf()

        if conf['debug']:
            print("Variation: \"" + folder + "\"")
        if conf['yield']:
            se_var = fdr.get_se()
            pair_num_var = se_var.Integral(se_var.FindBin(0), se_var.FindBin(conf['yield'][0]))
            deviation = abs(pair_num_se - pair_num_var) / pair_num_se
            if deviation > conf['yield'][1]:
                if conf['debug']:
                    dev = deviation * 100
                    print("Integrated yield k*: [0, " + str(conf['yield'][0]) + ") differs by " + f"{dev:.1f} %")
                    if deviation > conf['yield'][1]:
                        print("Variation: Excluded!\n")
                        continue
        if conf['debug'] and conf['htype'] != 'k':
            se_var_all = ch_var.get_se()
            tab = '\t'
            print("Differential yield:")
            for n, bin1 in enumerate(se_var_all):
                yield_all = se_all[n][0].Integral()
                yield_all_var = se_var_all[n][0].Integral()
                deviation = (abs(yield_all - yield_all_var) / yield_all) * 100
                print(f"{tab}{conf['htype']:s}:  [{conf['bins'][n]:.2f}, {conf['bins'][n + 1]:.2f}) {tab} {deviation:5.2f} %")
            print()

        for n, [ck_var, ck_var_rebin] in enumerate(cf_var):
            syst[n][0].AddVar(ck_var)
            if conf['rebin']:
                for i in range(len_rebin):
                    syst[n][1][i].AddVar(ck_var_rebin[i])
        del ch_var

    # generate th2 plots for systematics
    for n in range(len(syst)):
        syst[n][0].GenSyst()
        if conf['rebin']:
            for i in range(len_rebin):
                syst[n][1][i].GenSyst()

    syst_plots = []                                         # [[[cf, diff, syst, dev], [rebins]], [bin2...], ...]
    for n in range(len(syst)):
        syst_plots.append([syst[n][0].GetAll(), []])
        if conf['rebin']:
            for i in range(len_rebin):
                syst_plots[n][1].append(syst[n][1][i].GetAll())

    # generates the graphs with the systematic errors for the cf and the rebinned entries
    tgraphs = []
    for n, (hist, hist_rebin) in enumerate(cf_list):
        tgraphs.append([ROOT.TGraphErrors(), []])
        for i in range(1, hist.GetNbinsX() + 1):
            tgraphs[n][0].SetName("CF_syst_graph")
            tgraphs[n][0].SetPoint(i - 1, hist.GetBinCenter(i), hist.GetBinContent(i))
            tgraphs[n][0].SetPointError(i - 1, 0, syst_plots[n][0][2].GetBinContent(i))
        if conf['rebin']:
            for i in range(len_rebin):
                tgraphs[n][1].append(ROOT.TGraphErrors())
                for j in range(1, hist.GetNbinsX() + 1):
                    tgraphs[n][1][i].SetName("CF_syst_graph")
                    tgraphs[n][1][i].SetPoint(j - 1, hist_rebin[i].GetBinCenter(j), hist_rebin[i].GetBinContent(j))
                    tgraphs[n][1][i].SetPointError(j - 1, 0, syst_plots[n][1][i][2].GetBinContent(j))

    histos = (cf_list, syst_plots, tgraphs)
    fds = FDS.FemtoDreamSaver(conf, histos)

# systematics
def UFFA_syst_3d(settings):
    conf = config(settings)
    fdr = FDR.FemtoDreamReader(conf['fullpath'], conf['fileTDir'])

    # default cf
    ch = cf_handler(fdr, conf)
    histos = ch.get_cf_3d()                                # [[cf, [rebins]], [bin2...], ...], [[cf unw, [rebins]], [bin2...], ...]

    # input same event for yield filtering
    if conf['yield']:
        se = fdr.get_se()
        pair_num_se = se.Integral(se.FindBin(0), se.FindBin(conf['yield'][0]))
    if conf['debug']:
        se_all = ch.get_se_3d()

    syst = []
    syst_plots = []
    cf_raw = []

    if conf['rebin']:
        len_rebin = len(conf['rebin'])

    # create systematic object for all entries
    for n, bin1 in enumerate(histos):
        syst.append([])
        for nn, [cf, cf_rebin] in enumerate(bin1):
            syst[n].append([Systematics(cf), []])
            if conf['rebin']:
                for nnn in range(len_rebin):
                    syst[n][nn][1].append(Systematics(cf_rebin[nnn]))

    # loop over data variations in file and calculate the cf for each
    # which is then saved in a th2 from which the systematic error is computed and saved in a th1
    file_dir = fdr.get_dir();
    fdr.cd(0)                               # class method of FileSaver to return to root of file
    folders = fdr.get_folder_names()
    folder_counter = -1
    for folder in folders:
        fdr.cd(folder)

        # include/exclude specific variations
        if conf['exclude'] and folder in conf['exclude']:
            continue
        elif conf['include']:
            if folder in conf['include']:
                pass
            else:
                continue
        elif folder.rsplit('_')[-1][:3] != "Var":
            continue

        ch_var = cf_handler(fdr, conf)
        histos_var = ch_var.get_cf_3d()

        if conf['debug']:
            print("Variation: \"" + folder + "\"")
        # compare integrated yields in given range
        if conf['yield']:
            se_var = fdr.get_se()
            pair_num_var = se_var.Integral(se_var.FindBin(0), se_var.FindBin(conf['yield'][0]))
            deviation = abs(pair_num_se - pair_num_var) / pair_num_se
            if deviation > conf['yield'][1]:
                if conf['debug']:
                    dev = deviation * 100
                    print("Integrated yield k*: [0, " + str(conf['yield'][0]) + ") differs by " + f"{dev:.1f} %")
                    if deviation > conf['yield'][1]:
                        print("Variation: Excluded!\n")
                        continue
        if conf['debug']:
            se_var_all = ch_var.get_se_3d()
            tab = '\t'
            for n, bin1 in enumerate(se_var_all):
                print(f"Differential yield {conf['diff3d']:s}: [{conf['bins3d'][n]:.2f}, {conf['bins3d'][n + 1]:.2f})")
                for nn, bin2 in enumerate(bin1):
                    yield_all = se_all[n][nn][0].Integral()
                    yield_all_var = se_var_all[n][nn][0].Integral()
                    deviation = (abs(yield_all - yield_all_var) / yield_all) * 100
                    print(f"{tab}{conf['diff3d2']:s}:  [{conf['bins'][nn]:.2f}, {conf['bins'][nn + 1]:.2f}) {tab} {deviation:5.2f} %")
                print()
            if conf['interactive']:
                option = input("Include [Y/n] ")
                if option and option.lower()[0] == 'n':
                    print("\"" + folder + "\" excluded!\n")
                    continue
        folder_counter += 1

        cf_raw.append([])   # add entry for folder
        # add rebinned variations
        for n, bin1 in enumerate(histos_var):
            cf_raw[folder_counter].append([])
            for nn, [cf, cf_rebin] in enumerate(bin1):
                cf_raw[folder_counter][n].append([cf.Clone("CF_" + folder.rsplit('_')[-1]), []])
                syst[n][nn][0].AddVar(cf)
                if conf['rebin']:
                    for nnn in range(len_rebin):
                        cf_raw[folder_counter][n][nn][1].append(cf_rebin[nnn].Clone("CF_" + folder.rsplit('_')[-1]))
                        syst[n][nn][1][nnn].AddVar(cf_rebin[nnn])
        del ch_var

    # generate th2 plots for systematics
    for n, bin1 in enumerate(syst):
        syst_plots.append([])
        for nn, bin2 in enumerate(bin1):
            syst[n][nn][0].GenSyst()
            syst_plots[n].append([syst[n][nn][0].GetAll(), []])
            if conf['rebin']:
                for nnn in range(len_rebin):
                    syst[n][nn][1][nnn].GenSyst()
                    syst_plots[n][nn][1].append(syst[n][nn][1][nnn].GetAll())

    # generates the graphs with the systematic errors for the cf and the rebinned entries
    tgraphs = []
    for n, bin1 in enumerate(histos):
        tgraphs.append([])
        for nn, [hist, hist_rebin] in enumerate(bin1):
            tgraphs[n].append([ROOT.TGraphErrors(), []])
            for nnn in range(1, hist.GetNbinsX() + 1):
                tgraphs[n][nn][0].SetName("CF syst graph")
                tgraphs[n][nn][0].SetPoint(nnn - 1, hist.GetBinCenter(nnn), hist.GetBinContent(nnn))
                tgraphs[n][nn][0].SetPointError(nnn - 1, 0, syst_plots[n][nn][0][2].GetBinContent(nnn))
            if conf['rebin']:
                for nnn in range(len_rebin):
                    tgraphs[n][nn][1].append(ROOT.TGraphErrors())
                    for nnnn in range(1, hist.GetNbinsX() + 1):
                        tgraphs[n][nn][1][nnn].SetName("CF syst graph")
                        tgraphs[n][nn][1][nnn].SetPoint(nnnn - 1, hist_rebin[nnn].GetBinCenter(nnnn), hist_rebin[nnn].GetBinContent(nnnn))
                        tgraphs[n][nn][1][nnn].SetPointError(nnnn - 1, 0, syst_plots[n][nn][1][nnn][2].GetBinContent(nnnn))

    all_histos = (histos, syst_plots, tgraphs, cf_raw)
    fds = FDS.FemtoDreamSaver(conf, all_histos)

# class that returns the systematics of a cf
# add variations with AddVar(var) before calling GenSyst()
# GetAll() returns [th2 cf, th2 difference, th1 systematics, th1 std dev]
class Systematics():
    counter = 0
    ybins = 1200
    def __init__(self, cf):
        self._cf = cf
        self._xaxis = cf.GetXaxis()
        self._xbins = cf.GetNbinsX()

        self._var = ROOT.TH2D("CF_th2", "CF_th2", self._xbins, self._xaxis.GetXmin(), self._xaxis.GetXmax(), Systematics.ybins, 0, 10)
        self._dif = ROOT.TH2D("diff_th2", "diff_th2", self._xbins, self._xaxis.GetXmin(), self._xaxis.GetXmax(), Systematics.ybins, -5, 5)
        self._sys = ROOT.TH1D("syst_th1", "syst_th1", self._xbins, self._xaxis.GetXmin(), self._xaxis.GetXmax())
        self._dev = ROOT.TH1D("dev_th1", "dev_th1", self._xbins, self._xaxis.GetXmin(), self._xaxis.GetXmax())
        Systematics.counter = Systematics.counter + 1

    def AddVar(self, cf_var):
        for i in range(1, self._xbins + 1):
            if self._cf.GetBinCenter(i) > 3:    # break over 3GeV
                break
            self._var.Fill(cf_var.GetBinCenter(i), cf_var.GetBinContent(i))     # fill th2 cf histo with variation

    def GenSyst(self):
        for i in range(1, self._xbins + 1):
            if self._cf.GetBinCenter(i) > 3:    # break over 3GeV
                break
            cf_proj = self._var.ProjectionY("cf_xbin" + str(i), i, i)
            dev = cf_proj.GetStdDev()
            self._dev.SetBinContent(i, dev)
            cont_def = self._cf.GetBinContent(i)
            var_min = cf_proj.GetBinCenter(cf_proj.FindFirstBinAbove(0))
            var_max = cf_proj.GetBinCenter(cf_proj.FindLastBinAbove(0))
            self._dif.SetBinContent(i, self._dif.GetYaxis().FindBin(var_min), 1)
            self._dif.SetBinContent(i, self._dif.GetYaxis().FindBin(var_max), 1)
            dif_proj = self._dif.ProjectionY("diff_xbin" + str(i), i, i)
            proj_min = dif_proj.GetBinCenter(dif_proj.FindFirstBinAbove(0))     # minimum value of difference
            proj_max = dif_proj.GetBinCenter(dif_proj.FindLastBinAbove(0))      # maximum value of difference
            self._sys.SetBinContent(i, (proj_max - proj_min) / (12**0.5))       # assume a square distribution
        self._var.SetDirectory(0)
        self._dif.SetDirectory(0)
        self._sys.SetDirectory(0)
        self._dev.SetDirectory(0)

    def SetBinning(self, n):
        Systematics.ybins = n

    def GetVar(self):
        return self._var

    def GetDiff(self):
        return self._dif

    def GetSyst(self):
        return self._sys

    def GetDev(self):
        return self._dev

    def GetAll(self):
        return [self._var, self._dif, self._sys, self._dev]

# class that handles the retrieving of histos and computing of correlation functions
class cf_handler():
    def __init__(self, FileReader, conf):
        self._file  = FileReader
        self._pair  = conf['pair']
        self._atype = conf['atype']             # analysis type
        self._htype = conf['htype']             # histo type
        self._mc    = conf['mc']                # bool monte carlo data
        self._bins  = conf['bins']              # bin range for differential
        self._diff3d = conf['diff3d']           # which axis to split first in a 3D analysis
        self._bins3d = conf['bins3d']           # bin range for the first differential split in case of a 3D analysis
        self._rebin = conf['rebin']             # rebin factors for all se, me, cf plots
        self._norm  = conf['normalize']         # normalization range
        self._perc  = conf['percentile']        # percentile range
        self._rew_range = conf['rewrange']      # reweighting range
        self._name_se = conf['nameSE']
        self._name_me = conf['nameME']
        self._se = None
        self._me = None
        self._se_mc = None
        self._me_mc = None
        self._event = None
        self._tracks = None
        self._tracks_mc = None
        self._v0 = None
        self._get_histos()

    # retrieves histos from the provided file reader
    def _get_histos(self):
        if self._name_se and self._name_me:
            self._se = self._file.get_histo(self._name_se)
            self._me = self._file.get_histo(self._name_me)
        elif self._htype == 'k':        # TH1 kstar
            self._se, self._me = self._file.get_kstar()
            if self._mc:
                self._se_mc, self._me_mc = self._file.get_kstar_mc()
        elif self._htype == 'mult':      # TH2 k-mult
            self._se, self._me = self._file.get_kmult()
            if self._mc:
                self._se_mc, self._me_mc = self._file.get_kmult_mc()
        elif self._htype == 'mult3d':      # TH3 k-mult
            se, me = self._file.get_kmtmult()
            self._se = se.Project3D("zx").Clone()
            self._me = me.Project3D("zx").Clone()
            if self._mc:
                self._se_mc, self._me_mc = self._file.get_kmult_mc()
        elif self._htype == 'mt':      # TH2 k-mt
            self._se, self._me = self._file.get_kmt()
            if self._mc:
                self._se_mc, self._me_mc = self._file.get_kmt_mc()
        elif self._htype == 'mt3d':
            se, me = self._file.get_kmtmult()
            self._se = se.Project3D("yx").Clone()
            self._me = me.Project3D("yx").Clone()
            if self._mc:
                self._se_mc, self._me_mc = self._file.get_kmt_mc()
        elif self._htype in ['mtmult', 'rew3d']:
            self._se, self._me = self._file.get_kmtmult()
            if self._mc:
                self._se_mc, self._me_mc = self._file.get_kmtmult_mc()
        elif self._htype in ['4d', 'rew4d']:
            self._se, self._me = self._file.get_4d()

        if self._name_se and not self._name_me:
            self._se = self._file.get_histo(self._name_se)
        if self._name_me and not self._name_se:
            self._me = self._file.get_histo(self._name_me)

        if self._mc:
            self._tracks_mc = self._file.get_tracks_mc()
        if self._pair == 'pp':
            self._event = self._file.get_event()
            self._tracks = self._file.get_tracks()
        elif self._pair == 'pl':
            self._v0 = self._file.get_v0()

    # computes the cf for integrated or differential analysis and for mc data
    # and returns the histos for all the different options
    # [histos, histos_unw, histos_mc, histos_unw_mc, self._event, self._tracks, self._tracks_mc]
    def get_histos(self):
        histos = []
        histos_mc = []
        histos_unw = []
        histos_unw_mc = []
        if self._atype == 'int':        # integrated analysis
            histos, histos_unw = getIntegrated(self._se, self._me, self._htype, self._rebin, self._norm, self._rew_range)
            if self._mc:
                histos_mc, histos_unw_mc = getIntegrated(self._se_mc, self._me_mc, self._htype, self._rebin, self._norm, self._rew_range)
        elif self._atype == 'dif':      # differential analysis
            if self._htype == 'mtmult': # 3D differantial analysis
                histos = getDifferential3D(self._se, self._me, self._diff3d, self._bins3d, self._bins, self._rebin, self._norm)
            elif self._htype == 'rew3d':
                histos, histos_unw = getDiffReweight3D(self._se, self._me, self._diff3d, self._bins3d, self._bins, self._rebin, self._norm, self._rew_range)
            elif self._htype in ['4d', 'rew4d']:
                se3d, me3d = getProj4d(self._se, self._me, self._perc)
                if self._htype == '4d':
                    histos = getDifferential3D(se3d, me3d, self._diff3d, self._bins3d, self._bins, self._rebin, self._norm)
                else:
                    histos, histos_unw = getDiffReweight3D(se3d, me3d, self._diff3d, self._bins3d, self._bins, self._rebin, self._norm, self._rew_range)
            else:
                histos = getDifferential(self._se, self._me, self._htype, self._bins, self._rebin, self._norm)
                if self._mc:
                    histos_mc = getDifferential(self._se_mc, self._me_mc, self._htype, self._bins, self._rebin, self._norm)

        return [histos, histos_unw, histos_mc, histos_unw_mc, self._event, self._tracks, self._tracks_mc, self._v0]

    # returns a list of cf and their rebinned version
    # [[cf, [rebin 1, rebin 2, ...]], [bin2...], ...] same for unweighted if integrated analysis
    def get_cf(self):
        histos = []
        histos_unw = []
        cf_list = []
        cf_list_unw = []

        # integrated analysis
        if self._atype == 'int':
            histos, histos_unw = getIntegrated(self._se, self._me, self._htype, self._rebin, self._norm)
            if self._htype == 'mult':
                cf_list_unw.append(histos_unw[1])
                cf_list_unw.append([])
        # differential analysis
        elif self._atype == 'dif':
            histos = getDifferential(self._se, self._me, self._htype, self._bins, self._rebin, self._norm)
        cf_list.append([histos[1][2], []])                          # cf, for differential 1st bin

        # rebinned entries appended to the empty list for the first bin
        if self._rebin:
            for n in range(len(self._rebin)):
                cf_list[0][1].append(histos[1][3][n][2])            # rebinned cf
                if self._atype == 'int' and self._htype == 'mult':
                    cf_list_unw[1].append(histos_unw[2][n][1])      # rebinned unw cf for integrated

        # repeat for the rest of the bins in case of differential analysis
        if self._atype == 'dif':
            for n in range(2, len(self._bins)):
                cf_list.append([histos[n][2], []])
                if self._rebin:
                    for nn in range(len(self._rebin)):
                        cf_list[n][1].append(histos[n][3][nn][2])    # rebinned cf appended to rebin list
        return [cf_list, cf_list_unw]

    # returns a list of se and their rebinned version
    def get_se(self):
        histos = []
        se_list = []

        # integrated analysis
        if self._atype == 'int':
            histos, histos_unw = getIntegrated(self._se, self._me, self._htype, self._rebin, self._norm)
        # differential analysis
        elif self._atype == 'dif':
            histos = getDifferential(self._se, self._me, self._htype, self._bins, self._rebin, self._norm)
        se_list.append([histos[1][0], []])                          # se for differential 1st bin

        # rebinned entries appended to the empty list for the first bin
        if self._rebin:
            for n in range(len(self._rebin)):
                se_list[0][1].append(histos[1][3][n][0])            # rebinned se

        # repeat for the rest of the bins in case of differential analysis
        if self._atype == 'dif':
            for n in range(1, len(self._bins) - 1):
                se_list.append([histos[n][0], []])
                if self._rebin:
                    for nn in range(len(self._rebin)):
                        se_list[n][1].append(histos[n][3][nn][0])    # rebinned se appended to rebin list
        return se_list

    # returns all the cf's for a 3D mt/mult histo
    # [[[bin1-1 cf, [rebin cf]], [bin1-2 cf, [rebin cf]], ...], [[bin2-1 cf, [rebin cf]], [bin2-2 cf, [rebin cf]], ...], ...]
    def get_cf_3d(self):
        cf_list = []

        se = self._se
        me = self._me
        if self._htype in ['4d', 'rew4d']:
            se, me = getProj4d(self._se, self._me, self._perc)

        if self._htype in ['rew3d', 'rew4d']:
            histos, histos_unw = getDiffReweight3D(se, me, self._diff3d, self._bins3d, self._bins, self._rebin, self._norm, self._rew_range)
        else:
            histos = getDifferential3D(se, me, self._diff3d, self._bins3d, self._bins, self._rebin, self._norm)

        histos = histos[1:]     # remove TH3 histos
        for n, bin1 in enumerate(histos):
            cf_list.append([])
            bin1 = bin1[1:]     # remove TH2 histos
            for nn, th1 in enumerate(bin1):
                cf_list[n].append([th1[2], []])
                if self._rebin:
                    for nnn in range(len(self._rebin)):
                        cf_list[n][nn][1].append(th1[3][nnn][2])

        return cf_list

    # returns all the cf's for a 3D mt/mult histo
    def get_se_3d(self):
        se_list = []

        se = self._se
        me = self._me
        if self._htype in ['4d', 'rew4d']:
            se, me = getProj4d(self._se, self._me, self._perc)

        if self._htype in ['rew3d', 'rew4d']:
            histos, histos_unw = getDiffReweight3D(se, me, self._diff3d, self._bins3d, self.bins, self._rebin, self._nrm, self._rew_range)
        else:
            histos = getDifferential3D(se, me, self._diff3d, self._bins3d, self._bins, self._rebin, self._norm)

        histos = histos[1:]     # remove TH3 histos
        for n, bin1 in enumerate(histos):
            se_list.append([])
            bin1 = bin1[1:]     # remove TH2 histos
            for nn, th1 in enumerate(bin1):
                se_list[n].append([th1[0], []])
                if self._rebin:
                    for nnn in range(len(self._rebin)):
                        se_list[n][nn][1].append(th1[3][nnn][0])

        return se_list

# splits th2 in section based on provided bins
def getBinRangeHistos(iSE, iME, bins):
    """
    This function splits 2D histograms in the ranges
    defined in the option 'bins'.

    The output is a list of ["range", SE, ME] for each bin:
        [[name, SE, ME], [bin2], ...]
    where the name is a string containing the limits.
    """
    yAxis = iSE.GetYaxis()

    if type(bins) == list:
        limits = []
        for n in range(len(bins) - 1):
            bin_low = yAxis.FindBin(bins[n])
            bin_up  = FU.find_bin_reduce_on_lower_edge(yAxis, bins[n + 1])
            limits.append((bin_low, bin_up))
    else:
        print("Error in getBinRangeHistos: bin input \"" + str(bins) + "\" not a list of ranges!")
        exit()

    histos = []
    for n in range(len(limits)):
        name = "[%.2f-%.2f)" % (bins[n], bins[n + 1])
        se = iSE.ProjectionX("se_k", limits[n][0], limits[n][1])
        me = iME.ProjectionX("me_k", limits[n][0], limits[n][1])
        histos.append([name, se.Clone(), me.Clone()])

    return histos

# splits th3 in section based on provided bins
def getBinRangeHistos3D(iSE, iME, diff3d, bins3d):
    """
    This function takes as input 3D SE and ME plots
    and splits them into 2D plots in mt/mult according to 'diff3d'
    in the ranges defined in 'bins3d'.

    The output is a list of ["range", SE mt/mult vs k*, ME mt/mult vs k*] for each bin:
        [[name, SE, ME], [bin2], ...]
    where the name is a string containing the limits and SE, ME are 2D plots.
    """
    if diff3d == 'mt':
        diffAxisSE = iSE.GetYaxis()
        diffAxisME = iME.GetYaxis()
        projOpt = "zx"
    elif diff3d == 'mult':
        diffAxisSE = iSE.GetZaxis()
        diffAxisME = iME.GetZaxis()
        projOpt = "yx"
    else:
        print("Error in getBinRangeHistos: diff3d axis not known. Please choose either 'mt' or 'mult'")
        exit()

    if type(bins3d) == list:
        limits = []
        for n in range(len(bins3d) - 1):
            bin_low = diffAxisSE.FindBin(bins3d[n])
            bin_up  = FU.find_bin_reduce_on_lower_edge(diffAxisSE, bins3d[n + 1])
            limits.append((bin_low, bin_up))
    else:
        print("Error in getBinRangeHistos: bin input \"" + str(bins3d) + "\" not a list of ranges!")
        exit()

    histos = []
    for n in range(len(limits)):
        name = diff3d + ": [%.2f-%.2f)" % (bins3d[n], bins3d[n + 1])
        diffAxisSE.SetRange(limits[n][0], limits[n][1])
        diffAxisME.SetRange(limits[n][0], limits[n][1])
        se = iSE.Project3D(projOpt)
        me = iME.Project3D(projOpt)
        histos.append([name, se.Clone(f"SE_{projOpt}_{name}"), me.Clone(f"ME_{projOpt}_{name}")])

    return histos

# helper function for the CF
# output: [se, me, cf]
def getCorrelation(se, me, name, conf, norm = None):
    ch = CH.CorrelationHandler(name, se, me)
    #ch.normalize()
    ch.make_cf()
    minmax = norm if norm else [0.24, 0.34]
    ch.normalize_cf(minmax[0], minmax[1])
    se = ch.get_se().Clone("SE")
    me = ch.get_me().Clone("ME")
    cf = ch.get_cf().Clone("CF")
    se.SetTitle(conf)
    me.SetTitle(conf)
    cf.SetTitle(conf)
    del ch
    return [se, me, cf]

# returns [[iSE, iME], [se, me, cf]] for a list of mt or mult ranges
# [[iSE, iME], [se, me, cf, [rebin: [...], [...], ...], [bin 2 [rebin]], ...]
def getDifferential(iSE, iME, htype, bins, rebin, norm, title = None):
    histos = []
    conf = "" if not title else title + " "         # append to given name
    if htype == 'mult':
        conf += "mult: "
        histos.append([iSE.Clone("SE_kmult"), iME.Clone("ME_kmult")])
    elif htype == 'mt':
        conf += "mt: "
        histos.append([iSE.Clone("SE_kmT"), iME.Clone("ME_kmT")])
    else:
        print("getDifferential: no kmT or kmult input!")
        exit()

    mt_histos = getBinRangeHistos(iSE, iME, bins)
    for n, [name, se, me] in enumerate(mt_histos, 1):
        histos.append(getCorrelation(se, me, name, conf + name, norm))
        histos[n].append([])
        if rebin:       # append a list of rebinned [se, me, cf] in the original [se, me, cf, []]
            for factor in rebin:
                se_rebin = rebin_hist(se, factor)
                me_rebin = rebin_hist(me, factor)
                rebin_conf = " rebin: " + str(factor)
                histos[n][3].append(getCorrelation(se_rebin, me_rebin, name, conf + name + rebin_conf, norm))
    return histos

# [[iSE, iME], [[1st proj SE, 1st proj ME], [se, me, cf, [rebin], [bin 2 [rebin]]]], ...]
def getDifferential3D(iSE, iME, diff3d, bins3d, bins, rebin, norm):
    """
    This function takes as input 3D mult-mt-k* plots
    and splits them first in mt/mult according to 'diff3d' in the limits defined in 'bins3d'.
    The 2D mt/mult-k* plots are then projected in k* in the limits defined in 'bins'.

    The output is a list of the mt bins, each bin is then a list of the mult bins
    which include the SE, ME, CF and the rebinned plots:
        [ [ [SE, ME, CF, [rebinned SE, ME, CF]], [mult/mt bin2], ... ], [mt/mult bin2], ... ]
    """
    histos = []
    histos.append([iSE.Clone("SE_kmTmult"), iME.Clone("ME_kmTmult")])

    diff3d_histos = getBinRangeHistos3D(iSE, iME, diff3d, bins3d)

    diff_2 = 'mult'
    if (diff3d == "mult"):
        diff_2 = 'mt'

    for title, se, me in diff3d_histos:
        histos.append(getDifferential(se, me, diff_2, bins, rebin, norm, title))

    return histos

# [[iSE, iME], [[1st proj SE, 1st proj ME], [se, me, cf, [rebin], [bin 2 [rebin]]]], ...]
def getDiffReweight3D(iSE, iME, diff3d, bins3d, bins, rebin, norm, rew_range):
    """
    This function takes as input 3D mult-mt-k* plots
    and splits them first in mt according to the limits defined in 'bins3d'.

    The 2D mult-k* plots are then reweighted in mult and splits according
    to the limits defined in 'bins'.

    The output is a list of the mt bins, each bin is then a list of the mult bins
    which include the SE, ME, CF and the rebinned plots:
        [ [ [SE, ME, CF, [rebinned SE, ME, CF]], [mult bin2], ... ], [mt bin2], ... ]
    """
    histos = []
    histos.append([iSE.Clone("SE_kmTmult"), iME.Clone("ME_kmTmult")])

    diff_2 = 'mult'
    if (diff3d == "mult"):
        diff_2 = 'mt'

    histos_unw = getDifferential3D(iSE, iME, diff3d, bins3d, bins, rebin, norm)
    histos_diff3d = reweight3D(iSE, iME, diff3d, bins3d, rew_range)

    for title, se, me in histos_diff3d:
        histos.append(getDifferential(se, me, diff_2, bins, rebin, norm, title))

    return histos, histos_unw

# returns a list of [[iSE, iME], [se, me, cf]] for rel pair k* input
# or reweights and returns ([[iSE, iME], [se, me, cf, [rebin]]], [me_unw, cf_unw, [rebin]]) for kmult
# and [[iSE, iME], [se, me, cf, [rebin]]] for kmT
def getIntegrated(iSE, iME, htype, rebin, norm, rew_range):
    histos = []
    histos_unw = []
    if htype == 'k':      # k* input
        histos.append([iSE.Clone("SE_kstar"), iME.Clone("ME_kstar")])
        se = iSE
        me = iME
    elif htype == 'mult':    # kmult input
        histos.append([iSE.Clone("SE_kmult"), iME.Clone("ME_kmult")])
        hReweight = reweight(iSE, iME, rew_range)
        se          = hReweight[0]
        me          = hReweight[1]
        me_unw      = hReweight[2]
        se_mult     = hReweight[4]
        me_mult     = hReweight[5]
        me_mult_unw = hReweight[6]
        histos[0].append(se_mult.Clone("SE_mult"))
        histos[0].append(me_mult.Clone("ME_mult"))
        histos[0].append(me_mult_unw.Clone("ME_mult unw"))
    elif htype == 'mt':    # kmT input
        histos.append([iSE.Clone("SE_kmT"), iME.Clone("ME_kmT")])
        hReweight = reweight(iSE, iME, rew_range)
        se = hReweight[0]
        me = hReweight[2]   # unweighted me, i.e. normal me projection of the kmT histo

    histos.append(getCorrelation(se, me, "cf", "", norm))
    if rebin:               # append rebinned histos to list of histos
        histos_rebin = []
        for factor in rebin:
            se_rebin = rebin_hist(se, factor)
            me_rebin = rebin_hist(me, factor)
            rebin_conf = " rebin: " + str(factor)
            histos_rebin.append(getCorrelation(se_rebin, me_rebin, "rebin: " + str(factor), rebin_conf, norm))
        histos[1].append(histos_rebin)
    if htype == 'mult':      # 2nd list with unweighted histos
        se, me, cf = getCorrelation(se, me_unw, "cf_unw", "unweighted", norm)
        histos_unw.append(me.Clone("ME_unw"))
        histos_unw.append(cf.Clone("CF unw"))
        histos_unw.append([])
        if rebin:           # append rebinned histos to list of histos
            histos_rebin = []
            for factor in rebin:
                se_rebin = rebin_hist(se, factor)
                me_rebin = rebin_hist(me, factor)
                rebin_conf = " rebin: " + str(factor)
                se_rebin, me_rebin, cf_rebin = getCorrelation(se_rebin, me_rebin, "rebin: " + str(factor), rebin_conf, norm)
                histos_unw[2].append([me_rebin.Clone("ME_unw"), cf_rebin.Clone("CF_unw")])

    return histos, histos_unw

# projects and reweights se and me from kmult histos
# returns [se, me, me unweighted, se mult, me mult, me mult unweighted]
def reweight(iSE, iME, rew_range):
    """
    This function takes as input a 2D mult vs k* SE and ME distribution
    and reweights the ME distribution in each bin projection of mt/mult.

    The output is a list that includes all plots that can be generated:
        [0] SE 1D k*
        [1] ME 1D k* reweighted
        [2] ME 1D k* unweighted
        [3] ME 2D mt/mult vs k* reweighted
        [4] SE 1D mt/mult
        [5] ME 1D mt/mult reweighted
        [6] ME 1D mt/mult unweighted
    """

    me = iME.Clone("ME_kmult_reweighted")
    me.Reset("ICESM")
    me_axis = me.GetYaxis()

    se_k = iSE.ProjectionX("se_k")
    me_k = iME.ProjectionX("me_k")

    int_min = 0
    int_max = se_k.GetNbinsX()
    if rew_range:
        int_min = se_k.FindBin(rew_range[0])
        int_max = se_k.FindBin(FU.get_bin_reduce_on_lower_edge(iME.GetXaxis(), rew_range[1]))

    se_mult = iSE.ProjectionY("se_mult")
    me_mult = iME.ProjectionY("me_mult")

    me_k_unw = iME.ProjectionX("me_k_unw")
    me_mult_unw = iME.ProjectionY("me_mult_unw")

    me_k.Reset("ICESM")
    me_mult.Reset("ICESM")

    # loop for the projection of each multiplicity slice
    for ybin in range(1, iSE.GetNbinsY()):
        se_n = iSE.ProjectionX("se_bin", ybin, ybin)
        me_n = iME.ProjectionX("me_bin", ybin, ybin)

        se_int = se_n.Integral(int_min, int_max)
        me_int = me_n.Integral(int_min, int_max)

        if me_int > 0. and se_int > 0.:
            me_n.Scale(se_int / me_int)
            me_mult.SetBinContent(ybin, me_n.Integral(int_min, int_max))
            me_k.Add(me_n)
            for xbin in range(1, me_n.GetNbinsX() + 1):        # fill th2 reweighted ME
                #me.Fill(me_n.GetBinContent(xbin), me_axis.GetBinCenter(ybin))
                #me.Fill(me_axis.GetBinCenter(ybin), me_n.GetBinContent(xbin))
                me.SetBinContent(xbin, ybin, me_n.GetBinContent(xbin))

    return [se_k, me_k, me_k_unw, me, se_mult, me_mult, me_mult_unw]

# split th3 in mt range and reweight each slice in multiplicity
# output: [[name, mult-k SE, reweighted mult-k ME], [bin 2], ...]
def reweight3D(iSE, iME, diff3d, bins3d, rew_range):
    """
    This function takes as input a 3D SE and ME distribution
    and splits them in mt by the provided binning.
    The resulting 2D mult/k* plots are then reweighted.

    The output is a list for the individual mt bins with the name (mt limits), SE, ME:
        [[name, SE mult/k*, ME mult/k* reweighted], [bin 2], ...]
    """
    histos = getBinRangeHistos3D(iSE, iME, diff3d, bins3d)        # split th3 in mt and get a list of mult-k histos
    out = []

    for hist in histos:
        out.append([hist[0], hist[1], reweight(hist[1], hist[2], rew_range)[3].Clone(hist[0])])        # append the reweighted th2 ME distribution

    return out

# 4d percentile histos
def getProj4d(iSE, iME, perc_range):
    se4d = iSE.Clone("4d_perc_se")
    me4d = iME.Clone("4d_perc_me")

    axis = se4d.GetAxis(3)
    bin_low = axis.FindBin(perc_range[0])
    bin_up  = FU.find_bin_reduce_on_lower_edge(axis, perc_range[1])

    se4d.GetAxis(3).SetRange(bin_low, bin_up)
    me4d.GetAxis(3).SetRange(bin_low, bin_up)

    se = se4d.Projection(0, 1, 2)
    me = me4d.Projection(0, 1, 2)

    return [se, me]

# returns rebinned copy of histo
def rebin_hist(input_histo, binning):
    histo = input_histo.Clone()
    histo = histo.Rebin(binning)
    return histo

# generates list with rebin factors
def bin2list(rebin):
    rebin_list = []
    if type(rebin) == int or type(rebin) == str:
        rebin = [rebin]
    elif type(rebin) != list:
        return None
    rebin_list.extend(rebin)
    return rebin_list

# generates the proper settings dictionary
def config(dic_conf):
    """
    This function sets up all the configurable options
    so that it is consitent and expandable.

    Current options include:
            "function":     'cf', 'syst', 'tf' -> for correlation function, systematics, template fits
            "pair":         'pp', 'pl' -> for q&a plots relevant for individual analyses
            "path":         "string" -> full path to the root file, might include ~/ for home directory
            "file":         "string" -> name of the root file
            "fullpath":     "string" -> full path and file name equal to "path" + "file"
            "outDir":       "string" -> output directory
            "rename":       "string" -> rename output file
            "fileTDir":     "string" -> root file directory: path to directory inside the root file
            "nameSE":       "string" -> path + name of the se plot inside the provided "fileTDir" if given
            "nameME":       "string" -> same as nameSE but for the ME distribution
            "newfile":      'new', 'recreate', 'update' -> same option as in ROOT, 'new' will rename if file already exists
            "mc":           'true', 'false' -> save monte carlo data from provided root file
            "mcTDir":       "string" -> root file directory for the monte carlo data
            "bins":         [list of floats] -> binning for differential analysis
            "diff3d":       'mt', 'mult' -> project 3D plots first in mt/mult 2D and after in mult/mt 1D or vice versa
            "bins3d":       [list of floats] -> binning for 3D plots to 2D plots
            "yield":        [GeV, Deviation] -> integrated analysis: include systematics inside deviation for the GeV range
            "rebin":        int or [list of ints] -> rebin output plots
                                for tf: int -> all dca/cpa rebinned with int
                                or      [list of ints] -> each int will correspond to one range of the binning if provided
            "atype":        'int', 'dif' -> integrated analysis or differential analysis
            "htype":        'k', 'mt', 'mult', 'mt3d', 'mult3d', 'mtmult', 'rew3d'
                                'k'     -> k* - relative pair momentum distribution
                                'mt'  -> mt vs k* distribution
                                'mult'  -> multiplicity vs k* distribution
                                'mt3d' -> mt vs k* from 3D distribution but integrated in mult
                                'mult3d' -> mult vs k* from 3D distribution but integrated in mt
                                'mtmult' -> mult vs mt vs k* 3D distribution
                                'rew3d' -> mult vs mt vs k* 3D differentially in mt and reweighted in mult
            "tftype":       'dca', 'cpa' -> option for the template fit plots
            "templates":    [list of th1 plots] -> list of dca/cpa plots for fitting
            "temp_init":   list of values to initialize fitting parameters
            "temp_limits": list of limits of the fitting parameters
            "namelist":     [list of strings] -> names of dca/cpa plots for fitting
            "fitrange":     float -> fitrange for the template fitter
            "normalize":    [float, float] -> normalization range for the correlation function
            "percentile":   [low, upper] or int -> lower and upper edges or 0 to int for percentile cut
            "include":      "string" or [list of strings] -> include these variations in the systematics
            "exclude":      "string" or [list of strings] -> exclude these variations in the systematics
            "interactive":  'True', 'False' -> include/exclude interactively variations in terminal
            "debug":        'True', 'False' -> debug information in console
            "print":        'True', 'False' -> print canvas as png
    """

    # settings dictionary skeleton
    dic = {
            "function":         None,
            "pair":             None,
            "path":             "",
            "file":             None,
            "fullpath":         None,
            "fileTDir":         "",
            "nameSE":           "",
            "nameME":           "",
            "newfile":          None,
            "mc":               None,
            "mcTDir":           "",
            "outDir":           "",
            "rename":           None,
            "bins":             None,
            "bins3d":           None,
            "diff3d":           "",
            "diff3d2":          "",
            "yield":            None,
            "rebin":            None,
            "atype":            None,
            "htype":            None,
            "tftype":           None,
            "data":             None,
            "templates":        None,
            "temp_init":        None,
            "temp_limits":      None,
            "temp_fraction":    None,
            "namelist":         None,
            "fitrange":         None,
            "percentile":       None,
            "rewrange":         None,
            "normalize":        None,
            "include":          None,
            "exclude":          None,
            "debug":            False,
            "print":            False,
            "interactive":      False,
        }

    # keys to set values
    keys_k      = ['k', 'kstar']
    keys_mult   = ['mult', 'kmult']
    keys_mult3d = ['mult3d', 'kmult3d']
    keys_rew3d  = ['rew3d', 'rewmult']
    keys_mt     = ['mt', 'kmt']
    keys_mt3d   = ['mt3d', 'kmt3d']
    keys_mtmult = ['mtmult','kmtmult']
    keys_4d     = ['perc', '4d', '4dim', '4dims']
    keys_rew4d  = ['rew4d', 'rewperc']

    keys_int    = ['int', 'integrated']
    keys_dif    = ['diff', 'dif', 'differential']

    # initialize values
    entries = ['function',      # function to be used
               'path',          # path to input file
               'fileTDir',      # root folder where getter functions are used
               'mc',            # path to mc root file
               'mcTDir',        # mc folder inside root file
               'rename',        # rename output file
               'templates',     # list of template histos
               'namelist',      # list of template histos names
               'temp_init',     # initialize template fitting values
               'temp_limits',   # limits for template fitting parameters
               'temp_fraction', # set fraction with dictionary {'name', 'temp_init', 'temp_limits'}
               'fitrange',      # fit range for templates
               'signalrange',   # signal range to evaluate template fractions
               'normalize',     # range to normalize cf
               'data',          # not used
               'bins3d',        # bins to split 3d histo
               'bins',          # bins to split 2d histo
               'nameSE',        # root folder with SE histo
               'nameME',        # root folder with ME histo
               ]
    for entry in entries:
        if entry in dic_conf:
            dic[entry] = dic_conf[entry]

    # type of particle pair
    if 'pair' in dic_conf:
        if dic_conf['pair']:
            dic['pair'] = dic_conf['pair'].lower()

    # file name, file directory
    if 'file' in dic_conf:
        if dic_conf['file']:
            path_name = dic_conf['file'].rsplit('/', 1)
            if len(path_name) == 1:
                dic['file']  = path_name[0]
            else:
                dic['path'] = FU.path_expand(path_name[0]) + '/'
                dic['file']  = path_name[1]
            dic['fullpath'] = dic['path'] + dic['file']

    # output directory
    if 'outDir' in dic_conf:
        if dic_conf['outDir'] != "" and dic_conf['outDir']:
            dic['outDir'] = FU.path_expand(dic_conf['outDir'])
            if ROOT.gSystem.AccessPathName(dic['outDir']):
                print("output directory \"" + dic['outDir'] + "\" does not exist!")
                exit()
    else:
        dic['outDir'] = dic['path']

    # create file
    if 'newfile' in dic_conf:
        if dic_conf['newfile'] in [1, "new"]:
            dic['newfile'] = "new"
        if dic_conf['newfile'] in [2, "recreate"]:
            dic['newfile'] = "recreate"
        if dic_conf['newfile'] in [3, "update"]:
            dic['newfile'] = "update"

    # histogram and analysis type
    if 'type' in dic_conf:
        dic_conf['atype'] = dic_conf['type'][0]
        dic_conf['htype'] = dic_conf['type'][1]
        if len(dic_conf['type']) > 2:
            dic_conf['diff3d'] = dic_conf['type'][2]

    # analysis type
    if 'atype' in dic_conf:
        atype = dic_conf['atype']
        if type(atype) == str:
            atype = atype.lower()
        if atype in keys_int:
            dic['atype'] = 'int'
        elif atype in keys_dif:
            dic['atype'] = 'dif'

    # histogram type
    if 'htype' in dic_conf:
        htype = dic_conf['htype']
        if type(htype) == str:
            htype = htype.lower()
        if htype in keys_4d:
            dic['htype'] = '4d'
        elif htype in keys_rew4d:
            dic['htype'] = 'rew4d'
            dic['diff3d'] = 'mt'
        elif htype in keys_k:
            dic['htype'] = 'k'
        elif htype in keys_mult:
            dic['htype'] = 'mult'
        elif htype in keys_mult3d:
            dic['htype'] = 'mult3d'
        elif htype in keys_rew3d:
            dic['htype'] = 'rew3d'
            dic['diff3d'] = 'mt'
        elif htype in keys_mt:
            dic['htype'] = 'mt'
        elif htype in keys_mt3d:
            dic['htype'] = 'mt3d'
        elif htype in keys_mtmult:
            dic['htype'] = 'mtmult'

    # template fit type
    if 'tftype' in dic_conf:
        if dic_conf['tftype']:
            tftype = dic_conf['tftype'].lower()
            if tftype == 'dca':
                dic['tftype'] = 'dca'
            elif tftype == 'cpa':
                dic['tftype'] = 'cpa'

    # which axis to be used for the first split in a 3D analysis
    if 'diff3d' in dic_conf:
        diff3d = dic_conf['diff3d']
        if type(diff3d) == str:
            diff3d = diff3d.lower()
        if diff3d in keys_mult:
            dic['diff3d'] = 'mult'
            dic['diff3d2'] = 'mt'
        elif diff3d in keys_mt:
            dic['diff3d'] = 'mt'
            dic['diff3d2'] = 'mult'

    # yield setting to exclude systematic variations below a value of GeV that vary by a given percentage
    # input: [GeV, %]
    if 'yield' in dic_conf:
        dic['yield'] = dic_conf['yield']
        if dic['yield'] and len(dic['yield']) != 2:
            print("'yield' accepts [GeV, deviation], where GeV defines [0, GeV) and deviation the percentage!")

    # rebin factor/s
    if 'rebin' in dic_conf:
        dic['rebin'] = bin2list(dic_conf['rebin'])

    # rewrange
    if 'rewrange' in dic_conf:
        if dic_conf['rewrange'] and type(dic_conf['rewrange']) != list:
            print("'rewrange' accepts [min, max]!")
        dic['rewrange'] = dic_conf['rewrange']

    # percentile range
    if 'percentile' in dic_conf:
        if type(dic_conf['percentile']) == int:
            dic['percentile'] = [0, dic_conf['percentile']]
        elif type(dic_conf['percentile']) == list:
            if len(dic_conf['percentile']) == 1:
                dic['percentile'] = [0, dic_conf['percentile'][0]]
            else:
                dic['percentile'] = dic_conf['percentile']

    # include variations
    if 'include' in dic_conf:
        dic['include'] = bin2list(dic_conf['include'])

    # exclude variations
    if 'exclude' in dic_conf:
        dic['exclude'] = bin2list(dic_conf['exclude'])

    if 'debug' in dic_conf:
        dic['debug'] = bool(dic_conf['debug'])

    if 'print' in dic_conf:
        dic['print'] = bool(dic_conf['print'])

    if 'interactive' in dic_conf:
        dic['interactive'] = bool(dic_conf['interactive'])

    return dic

