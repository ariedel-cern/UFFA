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
    me_2d_name = "ME_2d"
    me_2d_reweighted_name = "ME_2d_Reweighted"
    se_1d_name = "SE"
    se_1d_normalized_name = "SE_Normalized"
    me_1d_name = "ME"
    me_1d_normalized_name = "ME_Normalized"
    me_1d_reweighted_name = "ME_Reweighted"
    me_1d_reweighted_normalized_name = "ME_Reweighted_Normalized"
    cf_name = "CF"
    cf_reweighted_name = "CF_Reweighted"

    def __init__(
        self,
        se,
        me,
        normalization_range,
        rebin_kstar=-1,
        axis_kstar=0,
        axis_reweight=-1,
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
            axis_kstar (int): Index of kstar axis (axis start counting from 0)
            axis_reweight (int): Index of axis to use for reweighting (axis start counting from 0, set to -1 to deactivate)
            ranges (list of tuples): List of tuples to indicate bin for each axis in a tuple (Min,Max). Empty tuple () means no cut
        """

        # declaring all class variables
        self._se = None  # original same event distribution
        self._me = None  # original mixed event distribution
        self._se_inRange = None  # same event distribution in range
        self._me_inRange = None  # mixed event distribution in range
        self._dim = None  # dimensions of SE/ME histograms
        self._axis_kstar = None  # kstar axis, axis are counted from 0
        self._axis_reweight = (
            None  # axis used for reweighting, axis are counted from 0
        )
        self._se_2d = (
            None  # same event distribution projected down to kstar vs reweight axis
        )
        self._me_2d = (
            None  # mixed event distribution projected down to kstar vs reweight axis
        )
        self._se_1d = None  # same event distribution projected down to kstar
        self._se_1d_normalized = None  # same event distribution projected down to kstar (normalized)
        self._me_1d = None  # mixed event distribution projected down to kstar
        self._me_1d_normalized = None  # mixed event distribution projected down to kstar, (normalized)
        self._cf = None  # correlation function
        self._me_2d_reweighted = None  # mixed event distribution projected down to kstar vs reweight axis (reweighted)
        self._me_1d_reweighted = (
            None  # mixed event distribution projected down to kstar (reweighted)
        )
        self._me_1d_reweighted_normalized = None  # mixed event distribution projected down to kstar (reweighted, normalized)
        self._cf_reweighted = None  # correlation function (reweighted)
        self._normalization_range = None  # normalization range
        self._rebin_kstar = None # rebin factor for kstar
        self._ranges = None  # list of ranges in each axis as tuples. If the tuple is empty (), no range is applied
        self._inheritsTH1 = False

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

        # get dimensions
        self._dim = hu.GetHistDimension(se)

        # Clone original SE/ME
        self._se = se.Clone(self.se_name)
        self._me = me.Clone(self.me_name)

        # check that normalization range is passed as a tuple of floats
        assert isinstance(normalization_range, tuple), "NormRange is not an tuple"
        assert isinstance(normalization_range[0], float), "Lower edge of NormRange is not a float"
        assert isinstance(normalization_range[1], float), "Upper edge of NormRange is not a float"
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
        assert isinstance(axis_kstar, int), "IndexKstarAxis is not an integer"

        # Check if index of kstar axis is smaller than number of dimensions
        # Axis indices are counted from 0, so the index can go from 0,..,N-1
        if self._dim - 1 < axis_kstar:
            raise ValueError(
                f"Index of kstar axis (={axis_kstar}) is larger or equal the histogram dimension (={self._dim})"
            )

        # Set kstar axis
        if self._dim == 1:
            self._axis_kstar = 0
            logger.debug("TH1 histogram detected, set index of kstar axis to 0")
        else:
            self._axis_kstar = axis_kstar
            logger.debug("Provided index of kstar axis: %d", self._axis_kstar)

        # Check if reweighting is activated
        if axis_reweight >= 0:
            # assert index of reweight axis is an integer
            assert isinstance(
                axis_reweight, int
            ), "IndexReweightAxis is not an integer"

            # check that dimensions of histograms is larger than 1, otherwise reweighting is not possible
            if self._dim <= 1:
                raise ValueError("1D Histogram provided, but reweighting is activated")

            # Check if index of reweight axis is smaller than number of dimensions
            # Axis indices are counted from 0, so the index can go from 0,..,N-1
            if self._dim - 1 < axis_reweight:
                raise ValueError(
                    f"Index of reweight axis (={axis_reweight}) is larger or equal the histograms dimension (={self._dim})"
                )

            # check that index of kstar axis is different from index of reweight axis
            if axis_kstar == axis_reweight:
                raise ValueError(
                    f"Index of kstar axis (={axis_kstar}) and reweight axis (={axis_reweight}) are the same"
                )

            # now set the reweight index
            self._axis_reweight = axis_reweight
            logger.debug("Provided index of reweight axis: %d", self._axis_reweight)

        # check type of rebin fator
        if rebin_kstar > 0:
            assert isinstance(
                rebin_kstar, int
            ), "Rebin factor is not an integer"
            self._rebin_kstar = rebin_kstar
            logger.debug("Provided rebin factor: %d", self._rebin_kstar)

        if ranges:
            # check that provided list has the correct number of dimensions
            if len(ranges) != self._dim:
                logger.critical(
                    "Number of bins (=%d) and number of dimension (#=%d) do not match",
                    len(ranges),
                    self._dim,
                )
                raise ValueError(
                    f"Number of bins (={len(ranges)}) and number of dimension (={self._dim}) do not match"
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
                        "Range (%f,%f) is configured in dimension %d",
                        bin[0],
                        bin[1],
                        dim,
                    )
                else:
                    logger.debug("No bin is configured in dimension %d", dim)
            self._ranges = ranges

    def DoRebin(self):
        """
        Rebin Se/Me histogram
        """
        logger.debug("Rebin SE/ME")
        if self._dim == 1:
            self._se.Rebin(self._rebin_kstar)
            self._me.Rebin(self._rebin_kstar)
        elif self._dim == 2:
            if self._axis_kstar == 0:
                self._se.RebinX(self._rebin_kstar)
                self._me.RebinX(self._rebin_kstar)
            else:
                self._se.RebinY(self._rebin_kstar)
                self._me.RebinY(self._rebin_kstar)
        elif self._dim == 3:
            if self._axis_kstar == 0:
                self._se.RebinX(self._rebin_kstar)
                self._me.RebinX(self._rebin_kstar)
            elif self._axis_kstar == 1:
                self._se.RebinY(self._rebin_kstar)
                self._me.RebinY(self._rebin_kstar)
            else:
                self._se.RebinZ(self._rebin_kstar)
                self._me.RebinZ(self._rebin_kstar)
        else:
            RebinFactors = np.ones((self._dim))
            RebinFactors[self._axis_kstar] = self._rebin_kstar
            RebinFactors = RebinFactors.astype(dtype=ct.c_int32)
            self._se = self._se.Rebin(RebinFactors)
            self._se.SetName(self.se_name)
            self._me = self._me.Rebin(RebinFactors)
            self._me.SetName(self.me_name)

    def DoRanges(self):
        """
        Apply ranges to same and mixed event distribution
        """
        logger.debug("Apply ranges to SE/ME")

        self._se_inRange = self._se.Clone(self.se_inRange_name)
        self._me_inRange = self._me.Clone(self.me_inRange_name)

        hu.SetHistRanges(self._se_inRange,self._ranges)
        hu.SetHistRanges(self._me_inRange,self._ranges)

    def DoProjections(self):
        """
        Project SE and ME onto 1d Histogram (kstar axis)
        """
        logger.debug(
            "Project SE/ME (Dim=%d) onto kstar axis (Index=%d)",
            self._dim,
            self._axis_kstar,
        )

        # use cuts if they are defined
        if self._ranges is None:
            se = self._se
            me = self._me
        else:
            se = self._se_inRange
            me = self._me_inRange

        if self._dim == 1:
            self._se_1d = se.Clone(self.se_1d_name)
            self._me_1d = me.Clone(self.me_1d_name)

        elif self._dim == 2:
            print(se)
            print(me)
            self._se_1d = cu.Proj2dTo1d(se, self._axis_kstar, self.se_1d_name)
            self._me_1d = cu.Proj2dTo1d(me, self._axis_kstar, self.me_1d_name)

        elif self._dim == 3:
            self._se_1d = cu.Proj3dTo1d(se, self._axis_kstar, self.se_1d_name)
            self._me_1d = cu.Proj3dTo1d(me, self._axis_kstar, self.me_1d_name)
        else:
            self._se_1d = cu.ProjNdTo1d(se, self._axis_kstar, self.se_1d_name)
            self._me_1d = cu.ProjNdTo1d(me, self._axis_kstar, self.me_1d_name)

    def DoProjectionsWithRw(self):
        """
        Project SE and ME onto 2d (kstar vs reweight) Histogram and 1d (kstar) axis Histogram
        """
        logger.debug(
            "Project SE/ME (Dim=%d) onto kstar axis (Index=%d) vs reweight axis (Index=%d)",
            self._dim,
            self._axis_kstar,
            self._axis_reweight,
        )

        # use cuts if they are defined
        if self._ranges is None:
            Se = self._se
            Me = self._me
        else:
            Se = self._se_inRange
            Me = self._me_inRange

        # if reweighting is activated, SE/ME need to hava at least 2 dimensions
        if self._dim == 2:
            # Swap x and y axis in case kstar axis is not index 0
            self._se_2d = cu.Proj2dTo2d(
                Se, self._axis_kstar, self._axis_reweight, self.se_2d_name
            )
            self._me_2d = cu.Proj2dTo2d(
                Me, self._axis_kstar, self._axis_reweight, self.me_2d_name
            )

        elif self._dim == 3:
            self._se_2d = cu.Proj3dTo2d(
                Se, self._axis_kstar, self._axis_reweight, self.se_2d_name
            )
            self._me_2d = cu.Proj3dTo2d(
                Me, self._axis_kstar, self._axis_reweight, self.me_2d_name
            )
        else:
            self._se_2d = cu.ProjNdTo2d(
                Se, self._axis_kstar, self._axis_reweight, self.se_2d_name
            )
            self._me_2d = cu.ProjNdTo2d(
                Me, self._axis_kstar, self._axis_reweight, self.me_2d_name
            )

        # reweight kstar vs reweight 2d histograms (kstar axis is x axis)
        self._me_2d_reweighted = cu.Reweight2d(self._se_2d, self._me_2d, self.me_2d_reweighted_name)
        # project onto kstar on x axis (-> index 0)
        self._me_1d_reweighted = cu.Proj2dTo1d(self._me_2d_reweighted, 0, self.me_1d_reweighted_name)

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
        if self._axis_reweight is not None:
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
            self._se_1d.Clone(self.se_1d_normalized_name), self._normalization_range
        )
        self._me_1d_normalized = cu.Normalize(
            self._me_1d.Clone(self.me_1d_normalized_name), self._normalization_range
        )
        self._cf = cu.Normalize(self._cf, self._normalization_range)

        if self._axis_reweight is not None:
            self._me_1d_reweighted_normalized = cu.Normalize(
                self._me_1d_reweighted.Clone(self.me_1d_reweighted_normalized_name), self._normalization_range
            )
            self._cf_reweighted = cu.Normalize(self._cf_reweighted, self._normalization_range)


    def ManageHistograms(self):
        """
        Call SetDirectory(0) on all histograms to obtain ownership
        """
        logger.debug("Gain ownership of histograms")
        if self._inheritsTH1:
            self._se.SetDirectory(0)
            self._me.SetDirectory(0)
            if self._ranges is not None:
                self._se_inRange.SetDirectory(0)
                self._me_inRange.SetDirectory(0)

        self._se_2d.SetDirectory(0)
        self._me_2d.SetDirectory(0)
        self._se_1d.SetDirectory(0)
        self._se_1d_normalized.SetDirectory(0)
        self._me_1d.SetDirectory(0)
        self._me_1d_normalized.SetDirectory(0)
        self._cf.SetDirectory(0)

        if self._axis_reweight is not None:
            self._me_2d_reweighted.SetDirectory(0)
            self._me_1d_reweighted.SetDirectory(0)
            self._me_1d_reweighted_normalized.SetDirectory(0)
            self._cf_reweighted.SetDirectory(0)

    def FinalTouch(self):
        """
        Run all steps to generate correlation function
        """
        logger.debug("Final Touch")
        # rebin kstar axis
        if self._rebin_kstar is not None:
            self.DoRebin()
        # apply cuts if defined
        if self._ranges is not None:
            self.DoRanges()
        # do projections onto kstar axis
        self.DoProjections()
        # do projections onto kstar vs reweight axis if reweight axis is defined
        if self._axis_reweight is not None:
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
        DirSe.Add(self._se, True)
        if self._ranges is not None:
            DirSe.Add(self._se_inRange, True)
        if self._axis_reweight is not None:
            DirSe.Add(self._se_2d, True)
        DirSe.Add(self._se_1d, True)
        DirSe.Add(self._se_1d_normalized, True)

        # create subdirectory for mixed event
        DirMe = rt.TDirectoryFile("ME", "ME", "", TdirFile)
        DirMe.Add(self._me, True)
        if self._ranges is not None:
            DirMe.Add(self._me_inRange, True)
        if self._axis_reweight is not None:
            DirMe.Add(self._me_2d, True)
            DirMe.Add(self._me_2d_reweighted, True)
        DirMe.Add(self._me_1d, True)
        DirMe.Add(self._me_1d_normalized, True)
        if self._axis_reweight is not None:
            DirMe.Add(self._me_1d_reweighted, True)
            DirMe.Add(self._me_1d_reweighted_normalized, True)

        # create subdirectory for correlation function
        DirCf = rt.TDirectoryFile("CF", "CF", "", TdirFile)
        DirCf.Add(self._cf, True)
        if self._axis_reweight is not None:
            DirCf.Add(self._cf_reweighted, True)

        DirSe.Write()
        DirMe.Write()
        DirCf.Write()
