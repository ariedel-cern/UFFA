import ROOT
import numpy as np


class CorrelationHandler:
    """
    ---------------------------------------------------
    Class used to compute the correlation function
    from same-event and mixed-event k* distributions.

    Parameters:
    __se = same-event k* distribution (TH1F)
    __me = mixed-event k* distribution (TH1F)
    ---------------------------------------------------
    """

    def __init__(self, name, iSE, iME):

        self.name = name

        self.__se = iSE.Clone(f'hSE_{self.name}')
        self.__me = iME.Clone(f'hME_{self.name}')

        self.__cf = None

        if not iSE.InheritsFrom(ROOT.TH1.Class()) and iSE.InheritsFrom(ROOT.TH2.Class()):
            print("Input SE is not a TH1!")
        if not iME.InheritsFrom(ROOT.TH1.Class()) and iSE.InheritsFrom(ROOT.TH2.Class()):
            print("Input ME is not a TH1!")

    def _normalize(self, hist):
        nBins = hist.GetNbinsX()
        hist.Scale(1 / hist.Integral(1, nBins + 1))

    def normalize(self):
        self.__se = self._normalize(self.__se)
        self.__me = self._normalize(self.__me)

    def make_cf(self):
        self.__cf = self.__se.Clone(f'hCF_{self.name}')
        self.__cf.Reset()
        self.__cf.GetYaxis().SetTitle('C(k*)')
        for i in range(1, self.__cf.GetNbinsX()):
            vSame = self.__se.GetBinContent(i)
            eSame = self.__se.GetBinError(i)
            vMixed = self.__me.GetBinContent(i)
            eMixed = self.__me.GetBinError(i)
            if vSame < 1.e-12 or vMixed < 1.e-12:
                continue
            else:
                vCF = vSame / vMixed
                eCF = np.sqrt((eMixed/vMixed)*(eMixed/vMixed) +
                              (eSame/vSame)*(eSame/vSame)) * vCF
                self.__cf.SetBinContent(i, vCF)
                self.__cf.SetBinError(i, eCF)

    def normalize_cf(self, low=0.6, high=1.):
        if not self.__cf:
            print('Compute the correlation function first')
            return
        low_bin = self.__cf.FindBin(low)
        low_edge = self.__cf.GetBinLowEdge(low_bin)
        high_bin = self.__cf.FindBin(high)
        high_edge = self.__cf.GetBinLowEdge(high_bin + 1)
        cf_integral = self.__cf.Integral(low_bin, high_bin, 'width')
        if cf_integral > 0:
            self.__cf.Scale((high_edge - low_edge) / cf_integral)

    ###   Getter functions   ###
    # Histogramms for same event distribution
    def get_se(self):
        return self.__se

    # Histograms for the mixed event distribution
    def get_me(self):
        return self.__me

    # Correlation function
    def get_cf(self):
        return self.__cf

