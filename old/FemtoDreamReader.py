import ROOT
import FileReader as FR

class FemtoDreamReader(FR.FileReader):
    def __init__(self, ifile, directory = None):
        if directory:
            if directory == "":
                directory = "femto-dream-pair-task-track-track"
            elif directory[0] == '_':
                directory = "femto-dream-pair-task-track-track" + directory
        FR.FileReader.__init__(self, ifile, directory)

    ### Getter Functions ###
    def get_pt(self):
        return self.GetHisto("Tracks_one/hPt")

    def get_pt_mc(self):
        return self.GetHisto("Tracks_one_MC/hPt_ReconNoFake")

    def get_dca(self):
        return self.GetHisto("Tracks_one/hDCAxy")

    def get_dca_mc(self):
        return [self.GetHisto("Tracks_one_MC/hDCAxy_Primary"), \
                self.GetHisto("Tracks_one_MC/hDCAxy_Secondary"), \
                self.GetHisto("Tracks_one_MC/hDCAxy_SecondaryDaughterLambda"), \
                self.GetHisto("Tracks_one_MC/hDCAxy_SecondaryDaughterSigmaplus"), \
                self.GetHisto("Tracks_one_MC/hDCAxy_Fake"), \
                self.GetHisto("Tracks_one_MC/hDCAxy_Material")]

    def get_zvxt(self):
        return self.GetHisto("Event/zvtxhist")

    def get_event(self):
        return self.GetHistos("Event")

    def get_tracks(self):
        return self.GetHistos("Tracks_one")

    def get_tracks_mc(self):
        return self.GetHistos("Tracks_one_MC")

    def get_v0(self):
        return [self.GetHistos("V0_two"), \
                self.GetHistos("V0Child_pos"), \
                self.GetHistos("V0Child_neg")]

    def get_se(self):
        return self.GetHisto("SameEvent/relPairDist")

    def get_me(self):
        return self.GetHisto("MixedEvent/relPairDist")

    def get_kstar(self):
        return self.GetHisto("SameEvent/relPairDist"), self.GetHisto("MixedEvent/relPairDist")

    def get_kstar_mc(self):
        return self.GetHisto("SameEvent_MC/relPairDist"), self.GetHisto("MixedEvent_MC/relPairDist")

    def get_kmt(self):
        return self.GetHisto("SameEvent/relPairkstarmT"), self.GetHisto("MixedEvent/relPairkstarmT")

    def get_kmt_mc(self):
        return self.GetHisto("SameEvent_MC/relPairkstarmT"), self.GetHisto("MixedEvent_MC/relPairkstarmT")

    def get_kmtmult(self):
        return self.GetHisto("SameEvent/relPairkstarmTMult"), self.GetHisto("MixedEvent/relPairkstarmTMult")

    def get_kmtmult_mc(self):
        return self.GetHisto("SameEvent_MC/relPairkstarmTMult"), self.GetHisto("MixedEvent_MC/relPairkstarmTMult")

    def get_kmult(self):
        return self.GetHisto("SameEvent/relPairkstarMult"), self.GetHisto("MixedEvent/relPairkstarMult")

    def get_kmult_mc(self):
        return self.GetHisto("SameEvent_MC/relPairkstarMult"), self.GetHisto("MixedEvent_MC/relPairkstarMult")

    def get_4d(self):
        return self.GetHisto("SameEvent/relPairkstarmTMultMultPercentile"), self.GetHisto("MixedEvent/relPairkstarmTMultMultPercentile")

