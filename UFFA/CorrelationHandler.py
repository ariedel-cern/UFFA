import ROOT as rt
import logging
import numpy as np
import ctypes as ct
from .Utils import CorrelationUtils as cu
from .Utils import HistUtils as hu

logger = logging.getLogger(__name__)


class CorrelationHandler:
    """
    CorrelationHandler class
    """

    # names of histograms
    se_name = "SE_Original"
    me_name = "ME_Original"
    se_inRange_name = "SE_inRange"
    me_inRange_name = "ME_inRange"
    se_2d_name = "SE_2d"
    se_2d_noRebin_name = "SE_2d_NoRebin"
    me_2d_name = "ME_2d"
    me_2d_noRebin_name = "ME_2d_NoRebin"
    me_2d_reweighted_name = "ME_2d_Reweighted"
    me_2d_reweighted_noRebin_name = "ME_2d_Reweighted_NoRebin"
    se_1d_name = "SE"
    s1_1d_normalized_name = "SE_Normalized"
    me_1d_name = "ME"
    me_1d_noRebin_name = "ME_NoRebin"
    me_1d_normalized_name = "ME_Normalized"
    me_1d_reweighted_name = "ME_Reweighted"
    me_1d_reweighted_noRebin_name = "ME_Reweighted_NoRebin"
    me_1d_reweighted_normalized_name = "ME_Reweighted_Normalized"
    cf_name = "CF"
    cf_recentered_name = "Cf_Recentered"
    cf_reweighted_name = "CF_Reweighted"

    def __init__(
        self,
        se,
        me,
        normalization_range,
        rebin_kstar=-1,
        index_kstar=0,
        index_reweight=-1,
        ranges=[],
    ):
        """
        CorrelationHanderl constructor
        Read in all information necessary to compute the correlation function
        Args:
            se (TH1/TH2/TH3/THnSparse): Same Event distribution
            me (same as SE): Mixed Event distribution
            normalization_range (tuple): Tuple containg upper and lower bound for normalization
            rebin_kstar (int): Rebin factor of kstar axis (set to -1 to deactivate)
            index_kstar (int): Index of kstar axis (axis start counting from 0)
            index_reweight (int): Index of axis to use for reweighting (axis start counting from 0, set to -1 to deactivate)
            ranges (list of tuples): List of tuples to indicate bin for each axis in a tuple (Min,Max). If no range should be applied, the value should be None
        """

        # declaring all class variables
        self._se = None  # original same event distribution
        self._me = None  # original mixed event distribution
        self._se_inRange = None  # same event distribution in bin
        self._me_inRange = None  # mixed event distribution in bin
        self._dimension = None  # dimensions of SE/ME histograms
        self._index_kstar = None  # index of kstar axis, axis are counted from 0
        self._index_reweight = (
            None  # index of axis used for reweighting, axis are counted from 0
        )
        self._se_2d = (
            None  # same event distribution projected down to kstar vs reweight axis
        )
        self._me_2d = (
            None  # mixed event distribution projected down to kstar vs reweight axis
        )
        self._me_2d_noRebin = None  # mixed event distribution projected down to kstar vs reweight axis (not rebinned)
        self._se_1d = None  # same event distribution projected down to kstar
        self._se_1d_normalized = None  # same event distribution projected down to kstar, normalized in NormRange
        self._me_1d = None  # mixed event distribution projected down to kstar
        self._me_1d_normalized = None  # mixed event distribution projected down to kstar, normalized in NormRange
        self._me_1d_noRebin = (
            None  # mixed event distribution projected down to kstar (not rebinned)
        )
        self._cf = None  # correlation function
        self._cf_recentered = None  # correlation function (recentered)
        self._me_2d_reweighted = None  # mixed event distribution projected down to kstar vs reweight axis (reweighted)
        self._me_2d_noRebin_reweighted = None  # mixed event distribution projected down to kstar vs reweight axis (not rebinned)
        self._me_1d_reweighted = (
            None  # mixed event distribution projected down to kstar (reweighted)
        )
        self._me_1d_noRebin_reweighted = None  # mixed event distribution projected down to kstar (reweighted and not rebinned)
        self._me_1d_reweighed_normalized = None  # mixed event distribution projected down to kstar (reweighted), normalized in NormRange
        self._cf_reweighted = None  # correlation function (reweighted)
        self._normalization_range = None  # normalization range
        self._rebin_kstar = None
        self._ranges = None  # list of ranges

        # bools
        self._setRanges = False
        self._reweighting = False
        self._rebinning = False

        ### From here: Argument santiszation and parsing

        # check that SE and ME are not null pointers
        if se == None or me == None:
            raise TypeError("Se or Me are null pointers")
        # check that SE and ME have the same class
        if not (se.IsA() == me.IsA()):
            raise TypeError("Histogram type of SE and ME do not match")

        # bool to check if Se/Me inherit from TH1 (i.e. do not inherit from THn)
        # if this is true, call SetDirectory(0) so CorrelationHandler can manage histograms (otherwise the might be name clashes)
        if se.InheritsFrom(rt.TH1.Class()):
            self._inheritsTH1 = True

        # Get histogram type of SE and ME to determine their dimensions
        self._dimension = hu.GetHistDimension(se)
        logger.debug("Number of dimensions of SE/ME: %d", self._dimension)

        # Clone original SE/ME
        self._se = se.Clone(self.se_name)
        self._me = me.Clone(self.me_name)

        # check that normalization range is passed as a tuple of floats
        assert isinstance(normalization_range, tuple), "NormRange is not an tuple"
        assert isinstance(
            normalization_range[0], float
        ), "Lower edge of NormRange is not a float"
        assert isinstance(
            normalization_range[1], float
        ), "Upper edge of NormRange is not a float"
        if normalization_range[0] > normalization_range[1]:
            raise ValueError(
                f"Lower edge of NormRange is larger than the upper edge (%{normalization_range[0]}, {normalization_range[1]})"
            )
        self._normalization_range = normalization_range
        logger.debug(
            "Normalization range: (%f,%f)",
            normalization_range[0],
            normalization_range[1],
        )

        # assert index of kstar axis is an integer
        assert isinstance(index_kstar, int), "IndexKstarAxis is not an integer"

        # Check if index of kstar axis is smaller than number of dimensions
        # Axis indices are counted from 0, so the index can go from 0,..,N-1
        if self._dimension - 1 < index_kstar:
            raise ValueError(
                f"Index of kstar axis (={index_kstar}) is larger or equal the histogram dimension (={self._dimension})"
            )

        # Set kstar axis
        if self._dimension == 1:
            self._index_kstar = 0
            logger.debug("TH1 histogram detected, set index of kstar axis to 0")
        else:
            self._index_kstar = index_kstar
            logger.debug("Provided index of kstar axis: %d", self._index_kstar)

        # Check if reweighting is activated
        if index_reweight >= 0:
            # assert index of reweight axis is an integer
            assert isinstance(
                index_reweight, int
            ), "IndexReweightAxis is not an integer"

            # check that dimensions of histograms is larger than 1, otherwise reweighting is not possible
            if self._dimension <= 1:
                raise ValueError("1D Histogram provided, but reweighting is activated")

            # Check if index of reweight axis is smaller than number of dimensions
            # Axis indices are counted from 0, so the index can go from 0,..,N-1
            if self._dimension - 1 < index_reweight:
                raise ValueError(
                    f"Index of reweight axis (={index_reweight}) is larger or equal the histograms dimension (={self._dimension})"
                )

            # check that index of kstar axis is different from index of reweight axis
            if index_kstar == index_reweight:
                raise ValueError(
                    f"Index of kstar axis (={index_kstar}) and reweight axis (={index_reweight}) are the same"
                )

            # now set the reweight index
            self._index_reweight = index_reweight
            self._reweighting = True
            logger.debug("Provided index of reweight axis: %d", self._index_reweight)

        # check type of rebin fator
        if rebin_kstar > 0:
            assert isinstance(rebin_kstar, int), "Rebin factor is not an integer"
            self._rebin_kstar = rebin_kstar
            self._rebinning = True
            logger.debug("Provided rebin factor: %d", self._rebin_kstar)

        if ranges:
            # check that provided list has the correct number of dimensions
            if len(ranges) != self._dimension:
                logger.critical(
                    "Number of bins (=%d) and number of dimension (#=%d) do not match",
                    len(ranges),
                    self._dimension,
                )
                raise ValueError(
                    f"Number of bins (={len(ranges)}) and number of dimension (={self._dimension}) do not match"
                )
            for dim, bin in enumerate(ranges):
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
            self._setRanges = True
            self._ranges = ranges

    def DoRangeLimits(self):
        """
        Apply range limits to same and mixed event distribution
        """
        logger.debug("Applying range limits to Se/Me distribution")

        self._se_inRange = self._se.Clone(self.se_inRange_name)
        self._me_inRange = self._me.Clone(self.me_inRange_name)

        hu.SetHistRanges(self._se_inRange, self._ranges)
        hu.SetHistRanges(self._me_inRange, self._ranges)

    def DoProjections1d(self):
        """
        Project SE and ME onto 1d Histogram (kstar axis)
        """
        logger.debug(
            "Project SE/ME (Dim=%d) onto kstar axis (Index=%d)",
            self._dimension,
            self._index_kstar,
        )

        # use cuts if they are defined
        if self._setRanges:
            se = self._se_inRange
            me = self._me_inRange
        else:
            se = self._se
            me = self._me

        if self._dimension == 1:
            self._se_1d = se.Clone(self.se_1d_name)
            self._me_1d = me.Clone(self.me_1d_name)

        elif self._dimension == 2:
            self._se_1d = cu.Proj2dTo1d(se, self._index_kstar, self.se_1d_name)
            self._me_1d = cu.Proj2dTo1d(me, self._index_kstar, self.me_1d_name)

        elif self._dimension == 3:
            self._se_1d = cu.Proj3dTo1d(se, self._index_kstar, self.se_1d_name)
            self._me_1d = cu.Proj3dTo1d(me, self._index_kstar, self.me_1d_name)
        else:
            self._se_1d = cu.ProjNdTo1d(se, self._index_kstar, self.se_1d_name)
            self._me_1d = cu.ProjNdTo1d(me, self._index_kstar, self.me_1d_name)

    def DoProjections2d(self):
        """
        Project SE and ME onto 2d (kstar vs reweight) Histogram and 1d (kstar) axis Histogram
        """
        logger.debug(
            "Project SE/ME (Dim=%d) onto kstar axis (Index=%d) vs reweight axis (Index=%d)",
            self._dimension,
            self._index_kstar,
            self._index_reweight,
        )

        # use cuts if they are defined
        if self._setRanges:
            se = self._se_inRange
            me = self._me_inRange
        else:
            se = self._se
            me = self._me

        if self._dimension == 2:
            # Swap x and y axis in case kstar axis is not index 0, otherwise just return histogram again
            self._se_2d = cu.Proj2dTo2d(
                se, self._index_kstar, self._index_reweight, self.se_2d_name
            )
            self._me_2d = cu.Proj2dTo2d(
                me, self._index_kstar, self._index_reweight, self.me_2d_name
            )

        elif self._dimension == 3:
            self._se_2d = cu.Proj3dTo2d(
                se, self._index_kstar, self._index_reweight, self.se_2d_name
            )
            self._me_2d = cu.Proj3dTo2d(
                me, self._index_kstar, self._index_reweight, self.me_2d_name
            )
        else:
            self._se_2d = cu.ProjNdTo2d(
                se, self._index_kstar, self._index_reweight, self.se_2d_name
            )
            self._me_2d = cu.ProjNdTo2d(
                me, self._index_kstar, self._index_reweight, self.me_2d_name
            )

    def DoRebin(self):
        """
        Rebin Se/Me histogram
        """
        logger.debug("Rebin SE/ME")

        self._se_rebin = self._se.Clone()
        self._me_rebin = self._me.Clone()

        rebi

        # rebin 1d histograms, keep copy of original mixed event distribution
        self._me_1d_noRebin = self._me_1d.Clone(self.me_1d_noRebin_name)
        self._me_1d.Rebin(self._rebin_kstar)
        self._se_1d.Rebin(self._rebin_kstar)

        # rebin 2d histograms, keep copy of original mixed event distribution
        if self._reweighting:
            self._se_2d_noRebin = self._se_2d.Clone(self.se_2d_noRebin_name)
            print("asdf")
            print(self._se_2d.GetMean(2))
            self._se_2d.RebinX(self._rebin_kstar)
            print(self._se_2d.GetMean(2))
            self._me_2d_noRebin = self._me_2d.Clone(self.me_2d_noRebin_name)
            self._me_2d.RebinX(self._rebin_kstar)

    def DoReweight(self):
        """
        Reweight Mixed event distributions
        """
        logger.debug("Reweight ME distribution")

        # reweight kstar vs reweight 2d histograms (kstar axis is x axis)
        self._me_2d_reweighted = cu.Reweight2d(
            self._se_2d, self._me_2d, self.me_2d_reweighted_name
        )
        # project onto kstar on x axis (-> index 0)
        self._me_1d_reweighted = cu.Proj2dTo1d(
            self._me_2d_reweighted, 0, self.me_1d_reweighted_name
        )

        if self._rebinning:
            self._me_2d_noRebin_reweighted = cu.Reweight2d(
                self._se_2d_noRebin, self._me_2d_noRebin, self.me_2d_reweighted_noRebin_name
            )
            self._me_1d_noRebin_reweighted = cu.Proj2dTo1d(
                self._me_2d_reweighted, 0, self.me_1d_reweighted_noRebin_name
            )

    def DoCorrelations(self):
        """
        Compute correlation functions
        """
        logger.debug("Compute correlation functions")
        self._cf = self._se_1d.Clone(self.cf_name)
        # active Sumw2 if not already for error computation
        if self._cf.GetSumw2N() == 0:
            self._cf.Sumw2(True)
        self._cf.Divide(self._me_1d)

        if self._reweighting:
            self._cf_reweighted = self._se_1d.Clone(self.cf_reweighted_name)
            # active Sumw2 if not already for error computation
            if self._cf_reweighted.GetSumw2N() == 0:
                self._cf_reweighted.Sumw2(True)
            self._cf_reweighted.Divide(self._me_1d_reweighted)

    def DoNormalizations(self):
        """
        Normalize correlation function and SE/ME in normalization range
        """
        logger.debug("Normalize correlation functions")
        self._se_1d_normalized = cu.Normalize(
            self._se_1d.Clone(self.s1_1d_normalized_name), self._normalization_range
        )
        self._me_1d_normalized = cu.Normalize(
            self._me_1d.Clone(self.me_1d_normalized_name), self._normalization_range
        )
        self._cf = cu.Normalize(self._cf, self._normalization_range)

        if self._reweighting:
            self._me_1d_reweighed_normalized = cu.Normalize(
                self._me_1d_reweighted.Clone(self.me_1d_reweighted_normalized_name),
                self._normalization_range,
            )
            self._cf_reweighted = cu.Normalize(
                self._cf_reweighted, self._normalization_range
            )

    def DoRecenter(self):
        """
        Generate TGraphError for correlation function, recentering bins with ME
        """
        logger.debug("Compute recentered correlation function")

        Nbins = self._cf.GetNbinsX()

        BinCenterX = []
        BinErrorX = []
        BinCenterY = []
        BinErrorY = []

        for bin in range(1, Nbins + 1):
            value = 0
            norm = 0
            error = 0
            for bins in range(
                1 + (bin - 1) * self._rebin_kstar, bin * self._rebin_kstar
            ):
                value = value + self._me_1d_noRebin.GetBinCenter(
                    bins
                ) * self._me_1d_noRebin.GetBinContent(bins)
                norm = norm + self._me_1d_noRebin.GetBinContent(bins)
            value = value / norm
            for bins in range(
                1 + (bin - 1) * self._rebin_kstar, bin * self._rebin_kstar
            ):
                error = error + self._me_1d_noRebin.GetBinContent(bins) * np.power(
                    self._me_1d_noRebin.GetBinCenter(bins) - value, 2
                )
            BinCenterX.append(value)
            BinErrorX.append(np.sqrt(error / norm))
            BinCenterY.append(self._cf.GetBinContent(bin))
            BinErrorY.append(self._cf.GetBinError(bin))

        BinCenterX = np.array(BinCenterX, dtype=ct.c_double)
        BinErrorX = np.array(BinErrorX, dtype=ct.c_double)
        BinCenterY = np.array(BinCenterY, dtype=ct.c_double)
        BinErrorY = np.array(BinErrorY, dtype=ct.c_double)

        self._cf_recentered = rt.TGraphErrors(
            Nbins, BinCenterX, BinCenterY, BinErrorX, BinErrorY
        )

    def ManageHistograms(self):
        """
        Call SetDirectory(0) on all histograms to obtain ownership
        """
        logger.debug("Gain explicit ownership of all histograms")

        if self._inheritsTH1:
            self._se.SetDirectory(0)
            self._me.SetDirectory(0)
            if self._setRanges:
                self._se_inRange.SetDirectory(0)
                self._me_inRange.SetDirectory(0)

        self._se_1d.SetDirectory(0)
        self._se_1d_normalized.SetDirectory(0)
        self._me_1d.SetDirectory(0)
        self._me_1d_normalized.SetDirectory(0)
        self._cf.SetDirectory(0)

        if self._reweighting:
            self._se_2d.SetDirectory(0)
            self._me_2d.SetDirectory(0)

            self._me_2d_reweighted.SetDirectory(0)
            self._me_1d_reweighted.SetDirectory(0)
            self._me_1d_reweighed_normalized.SetDirectory(0)
            self._cf_reweighted.SetDirectory(0)

    def FinalTouch(self):
        """
        Run all steps to generate correlation function
        """
        logger.debug("Running Final Touch")
        # rebin se/me distribution
        if self._rebinning:
            self.DoRebin()
        # apply cuts if defined
        if self._setRanges:
            self.DoRangeLimits()
        # do projections onto kstar axis
        self.DoProjections1d()
        # do projections onto kstar vs reweight axis if reweight axis is defined
        if self._reweighting:
            self.DoProjections2d()
        # rebin projected histograms
        if self._rebinning:
            self.DoRebin()
        # compute reweighted mixed event distributions
        if self._reweighting:
            self.DoReweight()
        # compute correlation function
        self.DoCorrelations()
        # normalize correlation function/distributions
        self.DoNormalizations()
        # recenter correlation function if rebinned
        if self._rebinning:
            self.DoRecenter()
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
        DirSe.Add(self._se, True)
        if self._setRanges:
            DirSe.Add(self._se_inRange, True)
        if self._reweighting:
            print("yash")
            print(self._se_2d.GetMean(2))
            DirSe.Add(self._se_2d, True)
        DirSe.Add(self._se_1d, True)
        DirSe.Add(self._se_1d_normalized, True)

        # create subdirectory for mixed event
        DirMe = rt.TDirectoryFile("ME", "ME", "", TdirFile)
        DirMe.Add(self._me, True)
        if self._setRanges:
            DirMe.Add(self._me_inRange, True)
        if self._reweighting:
            DirMe.Add(self._me_2d, True)
            DirMe.Add(self._me_2d_reweighted, True)
        DirMe.Add(self._me_1d, True)
        DirMe.Add(self._me_1d_normalized, True)
        if self._reweighting:
            DirMe.Add(self._me_1d_reweighted, True)
            DirMe.Add(self._me_1d_reweighed_normalized, True)

        # create subdirectory for correlation function
        DirCf = rt.TDirectoryFile("CF", "CF", "", TdirFile)
        DirCf.Add(self._cf, True)
        if self._index_reweight is not None:
            DirCf.Add(self._cf_reweighted, True)
        if self._rebinning:
            DirCf.Add(self._cf_recentered, True)

        DirSe.Write()
        DirMe.Write()
        DirCf.Write()
