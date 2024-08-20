import ROOT as rt
import logging
import numpy as np
import ctypes as ct
from .Utils import CorrelationUtils as cu

logger = logging.getLogger(__name__)


class CorrelationHandler:
    """
    CorrelationHandler class
    """

    # names of histograms
    SeName = "SE_Original"
    MeName = "ME_Original"
    SeWithCutName = "SE_with_Cuts"
    MeWithCutName = "ME_with_Cuts"
    Se2dName = "SE_2d"
    Me2dName = "ME_2d"
    Me2dRwName = "ME_2d_Reweighted"
    Se1dName = "SE"
    Se1dNormName = "SE_Normalized"
    Me1dName = "ME"
    Me1dNormName = "ME_Normalized"
    Me1dRwName = "ME_Reweighted"
    Me1dRwNormName = "ME_Reweighted_Normalized"
    CfName = "CF"
    CfRwName = "CF_Reweighted"

    def __init__(
        self,
        Se,
        Me,
        NormRange,
        RebinFactorKstarAxis=-1,
        IndexKstarAxis=0,
        IndexReweightAxis=-1,
        Bins=[],
    ):
        """
        CorrelationHanderl constructor
        Read in all information necessary to compute the correlation function
        Args:
            SE (TH1/TH2/TH3/THnSparse): Same Event distribution
            ME (same as SE): Mixed Event distribution
            NormRange (tuple): Tuple containg upper and lower bound for normalization
            IndexKstarAxis (int): Index of kstar axis (axis start counting from 0)
            IndexReweightAxis (int): Index of axis to use for reweighting (axis start counting from 0, set to -1 to deactivate)
            RebinFactorKstarAxis (int): Rebin factor of kstar axis (set to -1 to deactivate)
            Bins (list of tuples): List of tuples to indicate bin for each axis in a tuple (Min,Max)
        """

        # declaring all class variables
        self._Se = None  # original same event distribution
        self._Me = None  # original mixed event distribution
        self._SeInBin = None  # same event distribution in bin
        self._MeInBin = None  # mixed event distribution in bin
        self._HistDimension = None  # dimensions of SE/ME histograms
        self._IndexAxisKstar = None  # index of kstar axis, axis are counted from 0
        self._IndexAxisReweight = (
            None  # index of axis used for reweighting, axis are counted from 0
        )
        self._Se2d = (
            None  # same event distribution projected down to kstar vs reweight axis
        )
        self._Me2d = (
            None  # mixed event distribution projected down to kstar vs reweight axis
        )
        self._Se1d = None  # same event distribution projected down to kstar
        self._Se1dNorm = None  # same event distribution projected down to kstar, normalized in NormRange
        self._Me1d = None  # mixed event distribution projected down to kstar
        self._Me1dNorm = None  # mixed event distribution projected down to kstar, normalized in NormRange
        self._Cf = None  # correlation function
        self._Me2dRw = None  # mixed event distribution projected down to kstar vs reweight axis (reweighted)
        self._Me1dRw = (
            None  # mixed event distribution projected down to kstar (reweighted)
        )
        self._Me1dRwNorm = None  # mixed event distribution projected down to kstar (reweighted), normalized in NormRange
        self._CfRw = None  # correlation function (reweighted)
        self._NormRange = None  # normalization range
        self._RebinFactorKstarAxis = None
        self._Bins = None  # list of bins

        ### From here: Argument santiszation and parsing

        # check that SE and ME are not null pointers
        if Se == None or Me == None:
            raise TypeError("Se or Me are null pointers")
        # check that SE and ME have the same class
        if not (Se.IsA() == Me.IsA()):
            raise TypeError("Histogram type of SE and ME do not match")

        # bool to check if Se/Me inherit from TH1 (i.e. do not inherit from THn)
        # if this is true, call SetDirectory(0) so CorrelationHandler can manage histograms (otherwise the might be name clashes)
        self._InheritsTH1 = True

        # Get histogram type of SE and ME to determine their dimensions
        if Se.InheritsFrom(rt.THnBase.Class()):
            self._HistDimension = Se.GetNdimensions()
            self._InheritsTH1 = False
            logger.debug(
                "Detected histogram type: THn with %d dimensions", self._HistDimension
            )
        elif Se.InheritsFrom(rt.TH3.Class()):
            self._HistDimension = 3
            logger.debug("Detected histogram type: %s", Se.ClassName())
        elif Se.InheritsFrom(rt.TH2.Class()):
            self._HistDimension = 2
            logger.debug("Detected histogram type: %s", Se.ClassName())
        elif Se.InheritsFrom(rt.TH1.Class()):
            self._HistDimension = 1
            logger.debug("Detected histogram type: %s", Se.ClassName())
        else:
            raise TypeError(f"Histogram type not supported: {Se.ClassName()}")

        # Clone original SE/ME
        self._Se = Se.Clone(self.SeName)
        self._Me = Me.Clone(self.MeName)

        # check that normalization range is passed as a tuple of floats
        assert isinstance(NormRange, tuple), "NormRange is not an tuple"
        assert isinstance(NormRange[0], float), "Lower edge of NormRange is not a float"
        assert isinstance(NormRange[1], float), "Upper edge of NormRange is not a float"
        if NormRange[0] > NormRange[1]:
            raise ValueError(
                f"Lower edge of NormRange is larger than the upper edge (%{NormRange[0]}, {NormRange[1]})"
            )
        self._NormRange = NormRange
        logger.debug(
            "Normalization range: (%f,%f)",
            NormRange[0],
            NormRange[1],
        )

        # assert index of kstar axis is an integer
        assert isinstance(IndexKstarAxis, int), "IndexKstarAxis is not an integer"

        # Check if index of kstar axis is smaller than number of dimensions
        # Axis indices are counted from 0, so the index can go from 0,..,N-1
        if self._HistDimension - 1 < IndexKstarAxis:
            raise ValueError(
                f"Index of kstar axis (={IndexKstarAxis}) is larger or equal the histogram dimension (={self._HistDimension})"
            )

        # Set kstar axis
        if self._HistDimension == 1:
            self._IndexAxisKstar = 0
            logger.debug("TH1 histogram detected, set index of kstar axis to 0")
        else:
            self._IndexAxisKstar = IndexKstarAxis
            logger.debug("Provided index of kstar axis: %d", self._IndexAxisKstar)

        # Check if reweighting is activated
        if IndexReweightAxis >= 0:
            # assert index of reweight axis is an integer
            assert isinstance(
                IndexReweightAxis, int
            ), "IndexReweightAxis is not an integer"

            # check that dimensions of histograms is larger than 1, otherwise reweighting is not possible
            if self._HistDimension <= 1:
                raise ValueError("1D Histogram provided, but reweighting is activated")

            # Check if index of reweight axis is smaller than number of dimensions
            # Axis indices are counted from 0, so the index can go from 0,..,N-1
            if self._HistDimension - 1 < IndexReweightAxis:
                raise ValueError(
                    f"Index of reweight axis (={IndexReweightAxis}) is larger or equal the histograms dimension (={self._HistDimension})"
                )

            # check that index of kstar axis is different from index of reweight axis
            if IndexKstarAxis == IndexReweightAxis:
                raise ValueError(
                    f"Index of kstar axis (={IndexKstarAxis}) and reweight axis (={IndexReweightAxis}) are the same"
                )

            # now set the reweight index
            self._IndexAxisReweight = IndexReweightAxis
            logger.debug("Provided index of reweight axis: %d", self._IndexAxisReweight)

        # check type of rebin fator
        if RebinFactorKstarAxis > 0:
            assert isinstance(
                RebinFactorKstarAxis, int
            ), "Rebin factor is not an integer"
            self._RebinFactorKstarAxis = RebinFactorKstarAxis
            logger.debug("Provided rebin factor: %d", self._RebinFactorKstarAxis)

        if Bins:
            # check that provided list has the correct number of dimensions
            if len(Bins) != self._HistDimension:
                logger.critical(
                    "Number of bins (=%d) and number of dimension (#=%d) do not match",
                    len(Bins),
                    self._HistDimension,
                )
                raise ValueError(
                    f"Number of bins (={len(Bins)}) and number of dimension (={self._HistDimension}) do not match"
                )
            for dim, bin in enumerate(Bins):
                assert isinstance(bin, tuple), "bin range is not an tuple"
                if bin:
                    assert isinstance(
                        bin[0], float
                    ), "Lower edge of bin range is not a float"
                    assert isinstance(
                        bin[1], float
                    ), "Upper edge of bin range is not a float"
                    if bin[0] > bin[1]:
                        raise ValueError(
                            f"Lower edge of bin is larger than the upper edge ({bin[0]},{bin[1]})"
                        )
                    logger.debug(
                        "A bin (%f,%f) is configured in dimension %d",
                        bin[0],
                        bin[1],
                        dim,
                    )
                else:
                    logger.debug("No bin is configured in dimension %d", dim)
            self._Bins = Bins

    def DoRebin(self):
        """
        Rebin Se/Me histogram
        """
        logger.debug("Rebin SE/ME")
        if self._HistDimension == 1:
            self._Se.Rebin(self._RebinFactorKstarAxis)
            self._Me.Rebin(self._RebinFactorKstarAxis)
        elif self._HistDimension == 2:
            if self._IndexAxisKstar == 0:
                self._Se.RebinX(self._RebinFactorKstarAxis)
                self._Me.RebinX(self._RebinFactorKstarAxis)
            else:
                self._Se.RebinY(self._RebinFactorKstarAxis)
                self._Me.RebinY(self._RebinFactorKstarAxis)
        elif self._HistDimension == 3:
            if self._IndexAxisKstar == 0:
                self._Se.RebinX(self._RebinFactorKstarAxis)
                self._Me.RebinX(self._RebinFactorKstarAxis)
            elif self._IndexAxisKstar == 1:
                self._Se.RebinY(self._RebinFactorKstarAxis)
                self._Me.RebinY(self._RebinFactorKstarAxis)
            else:
                self._Se.RebinZ(self._RebinFactorKstarAxis)
                self._Me.RebinZ(self._RebinFactorKstarAxis)
        else:
            RebinFactors = np.ones((self._HistDimension))
            RebinFactors[self._IndexAxisKstar] = self._RebinFactorKstarAxis
            RebinFactors = RebinFactors.astype(dtype=ct.c_int32)
            self._Se = self._Se.Rebin(RebinFactors)
            self._Se.SetName(self.SeName)
            self._Me = self._Me.Rebin(RebinFactors)
            self._Me.SetName(self.MeName)

    def DoCuts(self):
        """
        Apply cuts to same and mixed event distribution
        """
        logger.debug("Apply cuts to SE/ME")

        self._SeInBin = self._Se.Clone(self.SeWithCutName)
        self._MeInBin = self._Me.Clone(self.MeWithCutName)

        if self._HistDimension == 1:
            if self._Bins[0]:
                self._SeInBin.GetXaxis().SetRangeUser(
                    self._Bins[0][0], self._Bins[0][1]
                )
                self._MeInBin.GetXaxis().SetRangeUser(
                    self._Bins[0][0], self._Bins[0][1]
                )
        elif self._HistDimension == 2:
            if self._Bins[0]:
                self._SeInBin.GetXaxis().SetRangeUser(
                    self._Bins[0][0], self._Bins[0][1]
                )
                self._MeInBin.GetXaxis().SetRangeUser(
                    self._Bins[0][0], self._Bins[0][1]
                )
            if self._Bins[1]:
                self._SeInBin.GetYaxis().SetRangeUser(
                    self._Bins[1][0], self._Bins[1][1]
                )
                self._MeInBin.GetYaxis().SetRangeUser(
                    self._Bins[1][0], self._Bins[1][1]
                )
        elif self._HistDimension == 3:
            if self._Bins[0]:
                self._SeInBin.GetXaxis().SetRangeUser(
                    self._Bins[0][0], self._Bins[0][1]
                )
                self._MeInBin.GetXaxis().SetRangeUser(
                    self._Bins[0][0], self._Bins[0][1]
                )
            if self._Bins[1]:
                self._SeInBin.GetYaxis().SetRangeUser(
                    self._Bins[1][0], self._Bins[1][1]
                )
                self._MeInBin.GetYaxis().SetRangeUser(
                    self._Bins[1][0], self._Bins[1][1]
                )
            if self._Bins[2]:
                self._SeInBin.GetZaxis().SetRangeUser(
                    self._Bins[2][0], self._Bins[2][1]
                )
                self._MeInBin.GetZaxis().SetRangeUser(
                    self._Bins[2][0], self._Bins[2][1]
                )
        else:
            for dim, cut in enumerate(self._Bins):
                if cut:
                    BinLow = self._SeInBin.GetAxis(dim).FindBin(cut[0])
                    BinHigh = cu.FindBinWithUpperEdgeDetection(
                        self._SeInBin.GetAxis(dim), cut[1]
                    )
                    self._SeInBin.GetAxis(dim).SetRange(BinLow, BinHigh)
                    self._MeInBin.GetAxis(dim).SetRange(BinLow, BinHigh)

    def DoProjections(self):
        """
        Project SE and ME onto 1d Histogram (kstar axis)
        """
        logger.debug(
            "Project SE/ME (Dim=%d) onto kstar axis (Index=%d)",
            self._HistDimension,
            self._IndexAxisKstar,
        )

        # use cuts if they are defined
        if self._Bins is None:
            Se = self._Se
            Me = self._Me
        else:
            Se = self._SeInBin
            Me = self._MeInBin

        if self._HistDimension == 1:
            self._Se1d = Se.Clone(self.Se1dName)
            self._Me1d = Me.Clone(self.Me1dName)

        elif self._HistDimension == 2:
            self._Se1d = cu.Proj2dTo1d(Se, self._IndexAxisKstar, self.Se1dName)
            self._Me1d = cu.Proj2dTo1d(Me, self._IndexAxisKstar, self.Me1dName)

        elif self._HistDimension == 3:
            self._Se1d = cu.Proj3dTo1d(Se, self._IndexAxisKstar, self.Se1dName)
            self._Me1d = cu.Proj3dTo1d(Me, self._IndexAxisKstar, self.Me1dName)
        else:
            self._Se1d = cu.ProjNdTo1d(Se, self._IndexAxisKstar, self.Se1dName)
            self._Me1d = cu.ProjNdTo1d(Me, self._IndexAxisKstar, self.Me1dName)

    def DoProjectionsWithRw(self):
        """
        Project SE and ME onto 2d (kstar vs reweight) Histogram and 1d (kstar) axis Histogram
        """
        logger.debug(
            "Project SE/ME (Dim=%d) onto kstar axis (Index=%d) vs reweight axis (Index=%d)",
            self._HistDimension,
            self._IndexAxisKstar,
            self._IndexAxisReweight,
        )

        # use cuts if they are defined
        if self._Bins is None:
            Se = self._Se
            Me = self._Me
        else:
            Se = self._SeInBin
            Me = self._MeInBin

        # if reweighting is activated, SE/ME need to hava at least 2 dimensions
        if self._HistDimension == 2:
            # Swap x and y axis in case kstar axis is not index 0
            self._Se2d = cu.Proj2dTo2d(
                Se, self._IndexAxisKstar, self._IndexAxisReweight, self.Se2dName
            )
            self._Me2d = cu.Proj2dTo2d(
                Me, self._IndexAxisKstar, self._IndexAxisReweight, self.Me2dName
            )

        elif self._HistDimension == 3:
            self._Se2d = cu.Proj3dTo2d(
                Se, self._IndexAxisKstar, self._IndexAxisReweight, self.Se2dName
            )
            self._Me2d = cu.Proj3dTo2d(
                Me, self._IndexAxisKstar, self._IndexAxisReweight, self.Me2dName
            )
        else:
            self._Se2d = cu.ProjNdTo2d(
                Se, self._IndexAxisKstar, self._IndexAxisReweight, self.Se2dName
            )
            self._Me2d = cu.ProjNdTo2d(
                Me, self._IndexAxisKstar, self._IndexAxisReweight, self.Me2dName
            )

        # reweight kstar vs reweight 2d histograms (kstar axis is x axis)
        self._Me2dRw = cu.Reweight2d(self._Se2d, self._Me2d, self.Me2dRwName)
        # project onto kstar on x axis (-> index 0)
        self._Me1dRw = cu.Proj2dTo1d(self._Me2dRw, 0, self.Me1dRwName)

    def DoCorrelations(self):
        """
        Compute correlation functions
        """
        logger.debug("Compute correlation functions")
        self._Cf = self._Se1d.Clone(self.CfName)
        # active Sumw2 if not already for error computation
        if self._Cf.GetSumw2N() == 0:
            self._Cf.Sumw2(True)
        self._Cf.Divide(self._Me1d)
        if self._IndexAxisReweight is not None:
            self._CfRw = self._Se1d.Clone(self.CfRwName)
            # active Sumw2 if not already for error computation
            if self._CfRw.GetSumw2N() == 0:
                self._CfRw.Sumw2(True)
            self._CfRw.Divide(self._Me1dRw)

    def DoNormalizations(self):
        """
        Normalize correlation function and SE/ME in normalization range
        """
        logger.debug("Normalize correlation functions")
        self._Se1dNorm = cu.Normalize(
            self._Se1d.Clone(self.Se1dNormName), self._NormRange
        )
        self._Me1dNorm = cu.Normalize(
            self._Me1d.Clone(self.Me1dNormName), self._NormRange
        )
        self._Cf = cu.Normalize(self._Cf, self._NormRange)

        if self._IndexAxisReweight is not None:
            self._Me1dRwNorm = cu.Normalize(
                self._Me1dRw.Clone(self.Me1dRwNormName), self._NormRange
            )
            self._CfRw = cu.Normalize(self._CfRw, self._NormRange)


    def ManageHistograms(self):
        """
        Call SetDirectory(0) on all histograms to obtain ownership
        """
        logger.debug("Gain ownership of histograms")
        if self._InheritsTH1:
            self._Se.SetDirectory(0)
            self._Me.SetDirectory(0)
            if self._Bins is not None:
                self._SeInBin.SetDirectory(0)
                self._MeInBin.SetDirectory(0)

        self._Se2d.SetDirectory(0)
        self._Me2d.SetDirectory(0)
        self._Se1d.SetDirectory(0)
        self._Se1dNorm.SetDirectory(0)
        self._Me1d.SetDirectory(0)
        self._Me1dNorm.SetDirectory(0)
        self._Cf.SetDirectory(0)

        if self._IndexAxisReweight is not None:
            self._Me2dRw.SetDirectory(0)
            self._Me1dRw.SetDirectory(0)
            self._Me1dRwNorm.SetDirectory(0)
            self._CfRw.SetDirectory(0)

    def FinalTouch(self):
        """
        Run all steps to generate correlation function
        """
        logger.debug("Final Touch")
        # rebin kstar axis
        if self._RebinFactorKstarAxis is not None:
            self.DoRebin()
        # apply cuts if defined
        if self._Bins is not None:
            self.DoCuts()
        # do projections onto kstar axis
        self.DoProjections()
        # do projections onto kstar vs reweight axis if reweight axis is defined
        if self._IndexAxisReweight is not None:
            self.DoProjectionsWithRw()
        # compute correlation function
        self.DoCorrelations()
        # normalize correlation function/distributions
        self.DoNormalizations()
        # gain ownership of all histograms
        self.ManageHistograms()

    def SaveOutput(self, TdirFile):
        """
        Save output to a TDirectoryFile
        Args:
            TdirFile (TDirectoryFile): TDirectoryFile for saving objects
        """
        logger.debug("Save output to TDirectoryFile: %s", TdirFile.GetName())
        # create subdirectory for same event
        DirSe = rt.TDirectoryFile("SE", "SE", "", TdirFile)
        DirSe.Add(self._Se, True)
        if self._Bins is not None:
            DirSe.Add(self._SeInBin, True)
        if self._IndexAxisReweight is not None:
            DirSe.Add(self._Se2d, True)
        DirSe.Add(self._Se1d, True)
        DirSe.Add(self._Se1dNorm, True)

        # create subdirectory for mixed event
        DirMe = rt.TDirectoryFile("ME", "ME", "", TdirFile)
        DirMe.Add(self._Me, True)
        if self._Bins is not None:
            DirMe.Add(self._MeInBin, True)
        if self._IndexAxisReweight is not None:
            DirMe.Add(self._Me2d, True)
            DirMe.Add(self._Me2dRw, True)
        DirMe.Add(self._Me1d, True)
        DirMe.Add(self._Me1dNorm, True)
        if self._IndexAxisReweight is not None:
            DirMe.Add(self._Me1dRw, True)
            DirMe.Add(self._Me1dRwNorm, True)

        # create subdirectory for correlation function
        DirCf = rt.TDirectoryFile("CF", "CF", "", TdirFile)
        DirCf.Add(self._Cf, True)
        if self._IndexAxisReweight is not None:
            DirCf.Add(self._CfRw, True)

        DirSe.Write()
        DirMe.Write()
        DirCf.Write()
