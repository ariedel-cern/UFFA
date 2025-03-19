import ROOT
import FileUtils as FU
import FileSaver as FS

# class that handles the saving of the histograms in the correct file structure
class FemtoDreamSaver(FS.FileSaver):
    last_edit = None

    def __init__(self, settings, histos = None):
        FS.FileSaver.__init__(self)
        self._histos = histos
        self._func  = settings['function']          # save structure for 'function'
        self._pair  = settings['pair']
        self._oname = settings['rename'] if settings['rename'] \
                else  settings['file']              # output file name
        self._nfile = settings['newfile']           # "new", "recreate" or "update"
        self._ipath = settings['path']              # input directory
        self._opath = settings['outDir']            # output directory
        self._idir  = settings['fileTDir']          # file TDirectory
        self._atype = settings['atype']             # analysis type
        self._htype = settings['htype']             # histo type
        self._mc    = settings['mc']                # bool monte carlo data
        self._bins  = settings['bins']              # bin range for differential
        self._diff3d = settings['diff3d']
        self._bins3d = settings['bins3d']
        self._rebin = settings['rebin']             # rebin factors for all se, me, cf plots
        self._debug = settings['debug']
        FS.FileSaver.setDebug(self._debug)
        if self._func == 'cf':
            self._save_cf()
        elif self._func in ['tf']:
            self._save_tf()
        elif self._func in ['ctf']:
            self._save_ctf()
        elif self._func in ['tf2d']:
            self._save_tf2d()
        elif self._func == 'syst':
            if self._htype in ['mtmult', 'rew3d', '4d', 'rew4d']:
                self._save_syst_3d()
            else:
                self._save_syst()

    # create the output file
    def _create_output(self, prepend):
        if self._nfile == 'new':
            new_name = FU.file_exists(self._opath + prepend + self._oname)
            FemtoDreamSaver.last_edit = new_name
        elif self._nfile == 'recreate':
            new_name = self._opath + prepend + self._oname
            FemtoDreamSaver.last_edit = None
        elif self._nfile == 'update':
            if FemtoDreamSaver.last_edit:
                new_name = FemtoDreamSaver.last_edit
            else:
                new_name = self._opath + prepend + self._oname
        self.touch(new_name, self._nfile)

    # function that saves all the histos in the correct file format
    # saver for UFFA_cf()
    def _save_cf(self):
        self._create_output("UFFA_")

        # histos setup
        hist_std = self._histos[0]
        hist_in = hist_std[0]               # [iSE, iME]
        hist_smc = hist_std[1:]             # [[se, me, cf, [rebin]], ...]
        if self._mc:
            hist_std_mc = self._histos[2]
            hist_in_mc  = hist_std_mc[0]    # [iSE, iME] <- MC
            hist_smc_mc = hist_std_mc[1:]   # [[se, me, cf, [rebin]], ...] <- MC
        if self._htype == 'mult':
            hist_unw = self._histos[1]
            if self._mc:
                hist_unw_mc = self._histos[3]
        if self._htype in ['rew3d', 'rew4d']:
            hist_unw_smc = self._histos[1][1:]

        hist_event = self._histos[4]
        hist_track = self._histos[5]
        hist_v0    = self._histos[7]
        if self._mc:
            hist_track_mc = self._histos[6]
            hist_pur = getPurity(hist_track[0], hist_track_mc[0])

        default = "femto-dream-pair-task-track-track"
        if self._idir == default:
            dir_root = self._file.mkdir(default + "_std")
        else:
            dir_root = self._file.mkdir(self._idir)
        dir_root.cd()

        if self._atype == 'int':                        # integrated
            self.write(hist_in)                             # iSE, iME
            self.write(hist_smc[0][:3])                     # se, me, cf
            if self._htype == 'mult':
                self.write(hist_unw[:2])                    # me unw, cf unw
            if self._rebin:                                 # all rebinned se, me, cf, me unw, cf unw, in individual directories
                for n in range(len(self._rebin)):
                    dir_rebin = dir_root.mkdir(f"rebin_{self._rebin[n]}")
                    dir_rebin.cd()
                    self.write(hist_smc[0][3][n])           # [[se, me, cf, [rebins]]]
                    if self._htype == 'mult':
                        self.write(hist_unw[2][n])          # [me unw, cf unw, [rebins]]
                    dir_root.cd()
                    del dir_rebin
            self.writeInDir(dir_root, "Event", hist_event)
            self.writeInDir(dir_root, "Tracks_one", hist_track)
            if self._pair == 'pl':
                self.writeInDir(dir_root, "V0_two", hist_v0[0])
                self.writeInDir(dir_root, "V0Child_pos", hist_v0[1])
                self.writeInDir(dir_root, "V0Child_neg", hist_v0[2])

            if self._mc:
                dir_mc = dir_root.mkdir("mc")
                dir_mc.cd()
                self.write(hist_in_mc)                      # [iSE, iME]
                self.write(hist_smc_mc[0][:3])              # [[se, me, cf, [rebins]]
                if self._htype == 'mult':
                    self.write(hist_unw_mc[:2])             # [me unw, cf unw, [rebins]]
                if self._rebin:                             # rebin directories in main directory
                    for n in range(len(self._rebin)):
                        dir_rebin = dir_mc.mkdir("rebin_{self._rebin[n]}")
                        dir_rebin.cd()
                        self.write(hist_smc_mc[0][3][n])
                        if self._htype == 'mult':
                            self.write(hist_unw_mc[2][n])
                        dir_mc.cd()
                        del dir_rebin
                self.writeInDir(dir_root, "Tracks_one_MC", hist_track_mc)
                hist_pur.Write()
        elif self._atype == 'dif':                      # differential
            if self._htype in ['mtmult', 'rew3d', '4d', 'rew4d']:                 # 3D histogramms
                self.write(hist_in)                                   # [iSE, iME]
                for n in range(len(self._bins3d) - 1):            # list: [1, 2, 3, 4] -> ranges: [1-2, 2-3, 3-4]
                    dir_bindiff3d = dir_root.mkdir(f"bin_{self._diff3d}_" + str(n + 1))
                    dir_bindiff3d.cd()

                    self.write(hist_smc[n][0])                        # 2D projection in the first bin of the 3D histogram
                    if self._htype in ['rew3d', 'rew4d']:
                        hist_unw_smc[n][0][1].Write("ME_kmult_unw")
                    for nm in range(1, len(self._bins)):
                        dir_bin = dir_bindiff3d.mkdir(f"bin_{nm}")
                        dir_bin.cd()                                # [[se, me, cf, [rebins]], ...]
                        self.write(hist_smc[n][nm][:3])             #   ^^  ^^  ^^
                                                                    #  rebin directories inside each mt/mult bin directory
                        if self._htype in ['rew3d', 'rew4d']:
                            hist_unw_smc[n][nm][1].Write("ME_unw")
                            hist_unw_smc[n][nm][2].Write("CF_unw")
                        if self._rebin:
                            for m in range(len(self._rebin)):
                                dir_rebin = dir_bin.mkdir(f"rebin_{self._rebin[m]}")
                                dir_rebin.cd()
                                self.write(hist_smc[n][nm][3][m])   # [[se, me, cf, [rebins]], ...]
                                if self._htype in ['rew3d', 'rew4d']:
                                    hist_unw_smc[n][nm][3][m][1].Write("ME_unw")
                                    hist_unw_smc[n][nm][3][m][2].Write("CF_unw")
                                dir_bin.cd()
                                del dir_rebin
                        dir_bindiff3d.cd()
                        del dir_bin

                    dir_root.cd()
                    del dir_bindiff3d
            else:
                self.write(hist_in)                             # [iSE, iME]
                for n in range(len(self._bins) - 1):            # list: [1, 2, 3, 4] -> ranges: [1-2, 2-3, 3-4]
                    dir_bin = dir_root.mkdir("bin_" + str(n + 1))
                    dir_bin.cd()
                    self.write(hist_smc[n][:3])                 # [[se, me, cf, [rebins]], ...]
                    if self._rebin:                             # rebin directories inside each mt/mult bin directory
                        for m in range(len(self._rebin)):
                            dir_rebin = dir_bin.mkdir(f"rebin_{self._rebin[m]}")
                            dir_rebin.cd()
                            self.write(hist_smc[n][3][m])       # [[se, me, cf, [rebins]], ...]
                            dir_bin.cd()
                            del dir_rebin
                    dir_root.cd()
                    del dir_bin

            self.writeInDir(dir_root, "Event", hist_event)
            self.writeInDir(dir_root, "Tracks_one", hist_track)
            if self._pair == 'pl':
                self.writeInDir(dir_root, "V0_two", hist_v0[0])
                self.writeInDir(dir_root, "V0Child_pos", hist_v0[1])
                self.writeInDir(dir_root, "V0Child_neg", hist_v0[2])

            if self._mc:                            # monte carlo directory
                dir_mc = dir_root.mkdir("mc")
                dir_mc.cd()
                self.write(hist_in_mc)                      # [iSE, iME]
                for n in range(len(self._bins) - 1):        # list: [1, 2, 3, 4] -> ranges: [1-2, 2-3, 3-4]
                    dir_bin = dir_mc.mkdir("bin_" + str(n + 1))
                    dir_bin.cd()
                    if self._rebin:                         # rebin directories inside each mt/mult bin directory
                        self.write(hist_smc_mc[n][:3])
                        for m in range(len(self._rebin)):
                            dir_rebin = dir_bin.mkdir(f"rebin_{self._rebin[m]}")
                            dir_rebin.cd()
                            self.write(hist_smc_mc[n][3][m])
                            dir_bin.cd()
                            del dir_rebin
                    del dir_bin
                del dir_mc
                self.writeInDir(dir_root, "Tracks_one_MC", hist_track_mc)

    # saver for UFFA_tf()
    def _save_tf(self):
        self._create_output("TemplateFit_")

    # saver for UFFA_ctf()
    def _save_tf2d(self):
        self._create_output("TemplateFit2D_")

    # saver for UFFA_ctf()
    def _save_ctf(self):
        self._create_output("CombinedFit_")

    # saver for UFFA_cf()
    def _save_syst(self):
        # create the output file
        self._create_output("UFFA_syst_")

        # [[syst cf, [rebins]], [bin2...], ...], [[[cf, diff, syst, dev], [rebins]], [bin2...], ...]
        syst, syst_plots, tgraphs = self._histos

        default = "femto-dream-pair-task-track-track"
        if self._idir == default:
            dir_root = self._file.mkdir(default + "_std")
        else:
            dir_root = self._file.mkdir(self._idir)
        dir_root.cd()

        if self._atype == 'int':
            syst[0][0].Write()
            tgraphs[0][0].Write()
            self.write(syst_plots[0][0])
            if self._rebin:
                for i, factor in enumerate(self._rebin):
                    dir_rebin = dir_root.mkdir(f"rebin_{factor}")
                    dir_rebin.cd()
                    syst[0][1][i].Write()
                    tgraphs[0][1][i].Write()
                    self.write(syst_plots[0][1][i])
                    dir_root.cd()
                    del dir_rebin
        elif self._atype == 'dif':
            for n in range(len(syst)):
                dir_bin = dir_root.mkdir("bin_" + str(n + 1))
                dir_bin.cd()
                syst[n][0].Write()
                tgraphs[n][0].Write()
                self.write(syst_plots[n][0])
                if self._rebin:
                    for i, factor in enumerate(self._rebin):
                        dir_rebin = dir_bin.mkdir(f"rebin_{factor}")
                        dir_rebin.cd()
                        syst[n][1][i].Write()
                        tgraphs[n][1][i].Write()
                        self.write(syst_plots[n][1][i])
                        dir_bin.cd()
                        del dir_rebin
                dir_root.cd()
                del dir_bin

    # [[[bin1-1 cf, [rebin cf]], [bin1-2 cf, [rebin cf]], ...], [[bin2-1 cf, [rebin cf]], [bin2-2 cf, [rebin cf]], ...], ...]
    # [[syst cf, [rebins]], [bin2...], ...], [[[cf, diff, syst, dev], [rebins]], [bin2...], ...] <---- not actually the real layout
    def _save_syst_3d(self):
        # create the output file
        self._create_output("UFFA_syst_")

        syst, syst_plots, tgraphs, cf_raw = self._histos

        default = "femto-dream-pair-task-track-track"
        if self._idir == default:
            dir_root = self._file.mkdir(default + "_std")
        else:
            dir_root = self._file.mkdir(self._idir)
        dir_root.cd()

        cf_title = syst[0][0][0].GetTitle().replace(":", "_").rsplit(' ')
        bin1_title = cf_title[0]
        bin2_title = cf_title[2]
        for n, bin1 in enumerate(syst):
            dir_bin1 = dir_root.mkdir(bin1_title + str(n + 1))
            dir_bin1.cd()
            for nn, bin2 in enumerate(bin1):
                dir_bin2 = dir_bin1.mkdir(bin2_title + str(nn + 1))
                dir_bin2.cd()
                syst[n][nn][0].Write()
                dir_cf = dir_bin2.mkdir("cf_var")
                dir_cf.cd()
                for folder in cf_raw:
                    folder[n][nn][0].Write()
                dir_bin2.cd()
                del dir_cf
                tgraphs[n][nn][0].Write()
                self.write(syst_plots[n][nn][0])
                if self._rebin:
                    for nnn, factor in enumerate(self._rebin):
                        dir_rebin = dir_bin2.mkdir(f"rebin_{factor}")
                        dir_rebin.cd()
                        syst[n][nn][1][nnn].Write()
                        dir_cf = dir_rebin.mkdir("cf_var")
                        dir_cf.cd()
                        for folder in cf_raw:
                            folder[n][nn][1][nnn].Write()
                        dir_rebin.cd()
                        tgraphs[n][nn][1][nnn].Write()
                        self.write(syst_plots[n][nn][1][nnn])
                        dir_bin2.cd()
                        del dir_rebin
                dir_bin1.cd()
                del dir_bin2
            del dir_bin1

