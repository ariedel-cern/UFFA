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
    cf_recentered_name = "CF_Recentered"
    cf_reweighted_recentered_name = "CF_Reweighted_Recentered"

    def __init__(
        self,
        se,
        me,
        normalization_range=(),
        rebin_kstar=0,
        rescale_kstar=-1,
        axis_kstar=0,
        axis_reweight=-1,
        ranges=[],
    ):
        """
        CorrelationHandler constructor
        Read in all information necessary to compute the correlation function
        Args:
            se (TH1/TH2/TH3/THnSparse): Same Event distribution
            me (same as SE): Mixed Event distribution
            normalization_range (tuple): Tuple containg upper and lower bound for normalization
            rebin_kstar (int): Rebin factor of kstar axis (set to 1 to deactivate)
            rescale_kstar (float): Rescaling factor for kstar axis. It rescales the axis for the correlation functions only (set to 0 to deactivate)
            axis_kstar (int): Index of kstar axis. Default is 0 (axis start counting from 0)
            axis_reweight (int): Index of axis to use for reweighting (axis start counting from 0, set to -1 to deactivate)
            ranges (list of tuples): List of tuples to indicate range for each axis in a tuple (Min,Max). Empty tuple () means full range
        """

        # declaring all class variables
        self.__handler_noRebin = None  # in case the user wants to rebin, generate a handler which is not rebinned
        self.__se = None  # original same event distribution
        self.__me = None  # original mixed event distribution
        self.__se_inRange = None  # same event distribution in range
        self.__me_inRange = None  # mixed event distribution in range
        self.__se_2d = (
            None  # same event distribution projected down to kstar vs reweight axis
        )
        self.__me_2d = (
            None  # mixed event distribution projected down to kstar vs reweight axis
        )
        self.__se_1d = None  # same event distribution projected down to kstar
        self.__se_1d_normalized = (
            None  # same event distribution projected down to kstar (normalized)
        )
        self.__me_1d = None  # mixed event distribution projected down to kstar
        self.__me_1d_normalized = (
            None  # mixed event distribution projected down to kstar, (normalized)
        )
        self.__cf = None  # correlation function
        self.__cf_reweighted = None  # correlation function (reweighted)
        self.__cf_recentered = None  # correlation function (recentered)
        self.__cf_reweighted_recentered = (
            None  # correlation function (reweighted and recentered)
        )
        self.__me_2d_reweighted = None  # mixed event distribution projected down to kstar vs reweight axis (reweighted)
        self.__me_1d_reweighted = (
            None  # mixed event distribution projected down to kstar (reweighted)
        )
        self.__me_1d_reweighted_normalized = None  # mixed event distribution projected down to kstar (reweighted, normalized)

        self.__dim = 1  # dimensions of SE/ME histograms
        self.__normalization_range = ()  # normalization range
        self.__rebin_kstar = 1  # rebin factor for kstar
        self.__rescale_kstar = 0  # rescale factor for kstar
        self.__axis_kstar = 1  # kstar axis, axis are counted from 0
        self.__axis_reweight = -1  # axis used for reweighting, axis are counted from 0
        self.__ranges = (
            []
        )  # list of ranges in each axis as tuples. If the tuple is empty (), no range is applied

        # booleans for easier configuration
        self.__do_normalization = False
        self.__do_rebin = False
        self.__do_rescale = False
        self.__do_reweight = False
        self.__do_ranges = False
        self.__inheritsTH1 = False

        ### From here: Argument santiszation and parsing

        # check that SE and ME are not null pointers
        if se == None or me == None:
            raise TypeError("Se or Me are null pointers")
        # check that SE and ME have the same class
        if not (se.IsA() == me.IsA()):
            raise TypeError("Histogram type of SE and ME do not match")

        # bool to check if Se/Me inherit from TH1 (i.e. do not inherit from THn)
        # if this is true, we call SetDirectory(0) at the end so CorrelationHandler can manage histograms
        if se.InheritsFrom(rt.TH1.Class()):
            self.__inheritsTH1 = True

        # get dimensions
        self.__dim = hu.GetHistDimension(se)

        # Clone original SE/ME
        self.__se = se.Clone(self.se_name)
        self.__me = me.Clone(self.me_name)

        # Check if index of kstar axis is smaller than number of dimensions
        # Axis indices are counted from 0, so the index can go from 0,..,N-1
        if self.__dim - 1 < axis_kstar:
            raise ValueError(
                f"Index of kstar axis (={axis_kstar}) is larger or equal the histogram dimension (={self.__dim})"
            )

        # Set kstar axis
        if self.__dim == 1:
            self.__axis_kstar = 0
            logger.debug("TH1 histogram detected, set index of kstar axis to 0")
        else:
            self.__axis_kstar = axis_kstar
            logger.debug("Provided index of kstar axis: %d", self.__axis_kstar)

        # set normalization range
        if len(normalization_range) == 2:
            if normalization_range[0] > normalization_range[1]:
                raise ValueError(
                    f"Lower edge of NormRange is larger than the upper edge ({normalization_range[0]}, {normalization_range[1]})"
                )
            self.__normalization_range = normalization_range
            self.__do_normalization = True
            logger.debug(
                "Normalization range: (%f,%f)",
                normalization_range[0],
                normalization_range[1],
            )

        # Check if reweighting is activated
        if axis_reweight >= 0:
            # assert index of reweight axis is an integer
            assert isinstance(axis_reweight, int), "IndexReweightAxis is not an integer"

            # check that dimensions of histograms is larger than 1, otherwise reweighting is not possible
            if self.__dim <= 1:
                raise ValueError("1D Histogram provided, but reweighting is activated")

            # Check if index of reweight axis is smaller than number of dimensions
            # Axis indices are counted from 0, so the index can go from 0,..,N-1
            if self.__dim - 1 < axis_reweight:
                raise ValueError(
                    f"Index of reweight axis (={axis_reweight}) is larger or equal to the histograms dimension (={self.__dim})"
                )

            # check that index of kstar axis is different from index of reweight axis
            if axis_kstar == axis_reweight:
                raise ValueError(
                    f"Index of kstar axis (={axis_kstar}) and reweight axis (={axis_reweight}) are the same"
                )

            # now set the reweight index
            self.__axis_reweight = axis_reweight
            self.__do_reweight = True
            logger.debug("Provided index of reweight axis: %d", self.__axis_reweight)

        # check rebin fator
        if rebin_kstar > 1:
            assert isinstance(rebin_kstar, int), "rebin_kstar is not an integer"
            self.__rebin_kstar = rebin_kstar
            self.__do_rebin = True
            logger.debug("Provided rebin factor: %d", self.__rebin_kstar)

        # check rescale factor
        if rescale_kstar > 0:
            self.__rescale_kstar = rescale_kstar
            self.__do_rescale = True
            logger.debug(
                "Provided rescaling factor for kstar axis: %d", self.__rescale_kstar
            )

        # check supplied ranges
        if ranges:
            # check that provided list has the correct number of dimensions
            if len(ranges) != self.__dim:
                raise ValueError(
                    f"Number of ranges (={len(ranges)}) and number of dimension (={self.__dim}) do not match"
                )
            for dim, bin in enumerate(ranges):
                if bin:
                    if bin[0] > bin[1]:
                        raise ValueError(
                            f"Lower edge of range is larger than the upper edge ({bin[0]},{bin[1]}) in dimension {self.__dim}"
                        )
                    logger.debug(
                        "Range (%f,%f) is configured in dimension %d",
                        bin[0],
                        bin[1],
                        dim,
                    )
                else:
                    logger.debug(
                        "No range is configured in dimension %d. Use full range", dim
                    )
            self.__do_ranges = True
            self.__ranges = ranges

    def DoNoRebin(self):
        """
        Generate CorrelationHandler object where rebin is turned off.
        With this object we can recentered the correlation function according to the ME distribution.
        """
        logger.debug("Generated CorrelationHandler without rebinning")
        self.__handler_noRebin = CorrelationHandler(
            self.__se,
            self.__me,
            self.__normalization_range,
            0,  # dont rebin
            self.__rescale_kstar,
            self.__axis_kstar,
            self.__axis_reweight,
            self.__ranges,
        )
        self.__handler_noRebin.FinalTouch()

    def DoRebin(self):
        """
        Rebin Se/Me histogram
        """
        logger.debug("Rebin SE/ME with %d", self.__rebin_kstar)
        if self.__inheritsTH1:
            if self.__dim == 1:
                self.__se.Rebin(self.__rebin_kstar)
                self.__me.Rebin(self.__rebin_kstar)
            elif self.__dim == 2:
                if self.__axis_kstar == 0:
                    self.__se.RebinX(self.__rebin_kstar)
                    self.__me.RebinX(self.__rebin_kstar)
                else:
                    self.__se.RebinY(self.__rebin_kstar)
                    self.__me.RebinY(self.__rebin_kstar)
            elif self.__dim == 3:
                if self.__axis_kstar == 0:
                    self.__se.RebinX(self.__rebin_kstar)
                    self.__me.RebinX(self.__rebin_kstar)
                elif self.__axis_kstar == 1:
                    self.__se.RebinY(self.__rebin_kstar)
                    self.__me.RebinY(self.__rebin_kstar)
                else:
                    self.__se.RebinZ(self.__rebin_kstar)
                    self.__me.RebinZ(self.__rebin_kstar)
        else:
            RebinFactors = np.ones((self.__dim))
            RebinFactors[self.__axis_kstar] = self.__rebin_kstar
            RebinFactors = RebinFactors.astype(dtype=ct.c_int32)
            self.__se = self.__se.Rebin(RebinFactors)
            self.__se.SetName(self.se_name)
            self.__me = self.__me.Rebin(RebinFactors)
            self.__me.SetName(self.me_name)

    def DoRanges(self):
        """
        Apply ranges to same and mixed event distribution
        """
        logger.debug("Apply ranges to SE/ME: %s", self.__ranges)

        self.__se_inRange = self.__se.Clone(self.se_inRange_name)
        self.__me_inRange = self.__me.Clone(self.me_inRange_name)

        hu.SetHistRanges(self.__se_inRange, self.__ranges)
        hu.SetHistRanges(self.__me_inRange, self.__ranges)

    def DoProjections(self):
        """
        Project SE and ME onto 1d Histogram (kstar axis)
        """
        logger.debug(
            "Project SE/ME (Dim=%d) onto kstar axis (Index=%d)",
            self.__dim,
            self.__axis_kstar,
        )

        # use ranged se/me if they are defined
        if self.__do_ranges:
            se = self.__se_inRange
            me = self.__me_inRange
        else:
            se = self.__se
            me = self.__me

        if self.__inheritsTH1:
            if self.__dim == 1:
                self.__se_1d = se.Clone(self.se_1d_name)
                self.__me_1d = me.Clone(self.me_1d_name)

            elif self.__dim == 2:
                self.__se_1d = cu.Proj2dTo1d(se, self.__axis_kstar, self.se_1d_name)
                self.__me_1d = cu.Proj2dTo1d(me, self.__axis_kstar, self.me_1d_name)

            elif self.__dim == 3:
                self.__se_1d = cu.Proj3dTo1d(se, self.__axis_kstar, self.se_1d_name)
                self.__me_1d = cu.Proj3dTo1d(me, self.__axis_kstar, self.me_1d_name)
        else:
            self.__se_1d = cu.ProjNdTo1d(se, self.__axis_kstar, self.se_1d_name)
            self.__me_1d = cu.ProjNdTo1d(me, self.__axis_kstar, self.me_1d_name)

    def DoProjectionsWithRw(self):
        """
        Project SE and ME onto 2d (kstar vs reweight) Histogram and 1d (kstar) axis Histogram
        """
        logger.debug(
            "Project SE/ME (Dim=%d) onto kstar axis (Index=%d) vs reweight axis (Index=%d)",
            self.__dim,
            self.__axis_kstar,
            self.__axis_reweight,
        )

        # use cuts if they are defined
        if self.__do_ranges:
            Se = self.__se_inRange
            Me = self.__me_inRange
        else:
            Se = self.__se
            Me = self.__me

        # if reweighting is activated, SE/ME need to hava at least 2 dimensions
        if self.__inheritsTH1:
            if self.__dim == 2:
                # Swap x and y axis in case kstar axis is not index 0
                self.__se_2d = cu.Proj2dTo2d(
                    Se, self.__axis_kstar, self.__axis_reweight, self.se_2d_name
                )
                self.__me_2d = cu.Proj2dTo2d(
                    Me, self.__axis_kstar, self.__axis_reweight, self.me_2d_name
                )

            elif self.__dim == 3:
                self.__se_2d = cu.Proj3dTo2d(
                    Se, self.__axis_kstar, self.__axis_reweight, self.se_2d_name
                )
                self.__me_2d = cu.Proj3dTo2d(
                    Me, self.__axis_kstar, self.__axis_reweight, self.me_2d_name
                )
        else:
            self.__se_2d = cu.ProjNdTo2d(
                Se, self.__axis_kstar, self.__axis_reweight, self.se_2d_name
            )
            self.__me_2d = cu.ProjNdTo2d(
                Me, self.__axis_kstar, self.__axis_reweight, self.me_2d_name
            )

        # reweight kstar vs reweight 2d histograms (kstar axis is x axis)
        self.__me_2d_reweighted = cu.Reweight2d(
            self.__se_2d, self.__me_2d, self.me_2d_reweighted_name
        )
        # project onto kstar on x axis (-> index 0)
        self.__me_1d_reweighted = cu.Proj2dTo1d(
            self.__me_2d_reweighted, 0, self.me_1d_reweighted_name
        )

    def DoCorrelations(self):
        """
        Compute correlation functions
        """
        logger.debug("Compute correlation functions")
        self.__cf = self.__se_1d.Clone(self.cf_name)
        # active Sumw2 if not already for error computation
        if self.__cf.GetSumw2N() == 0:
            self.__cf.Sumw2(True)
        self.__cf.Divide(self.__me_1d)
        if self.__do_reweight:
            self.__cf_reweighted = self.__se_1d.Clone(self.cf_reweighted_name)
            # active Sumw2 if not already for error computation
            if self.__cf_reweighted.GetSumw2N() == 0:
                self.__cf_reweighted.Sumw2(True)
            self.__cf_reweighted.Divide(self.__me_1d_reweighted)

    def DoNormalizations(self):
        """
        Normalize correlation function and SE/ME in normalization range
        """
        logger.debug(
            "Normalize correlation functions in range (%f,%f)",
            self.__normalization_range[0],
            self.__normalization_range[1],
        )

        self.__se_1d_normalized = cu.Normalize(
            self.__se_1d.Clone(self.se_1d_normalized_name), self.__normalization_range
        )
        self.__me_1d_normalized = cu.Normalize(
            self.__me_1d.Clone(self.me_1d_normalized_name), self.__normalization_range
        )
        self.__cf = cu.Normalize(self.__cf, self.__normalization_range)

        if self.__do_reweight:
            self.__me_1d_reweighted_normalized = cu.Normalize(
                self.__me_1d_reweighted.Clone(self.me_1d_reweighted_normalized_name),
                self.__normalization_range,
            )
            self.__cf_reweighted = cu.Normalize(
                self.__cf_reweighted, self.__normalization_range
            )

    def DoRescale(self):
        """
        Rescale kstar axis
        """
        logger.debug(
            "Rescale correlation functions and distributions with %f",
            self.__rescale_kstar,
        )

        self.__se_1d = hu.RescaleHist(self.__se_1d, self.__rescale_kstar, 0)
        self.__se_1d_normalized = hu.RescaleHist(
            self.__se_1d_normalized, self.__rescale_kstar, 0
        )
        self.__me_1d = hu.RescaleHist(self.__me_1d, self.__rescale_kstar, 0)
        self.__me_1d_normalized = hu.RescaleHist(
            self.__me_1d_normalized, self.__rescale_kstar, 0
        )
        self.__cf = hu.RescaleHist(self.__cf, self.__rescale_kstar, 0)

        # if we reweight, then we also need the 2d histograms
        if self.__do_reweight:
            self.__se_2d = hu.RescaleHist(self.__se_2d, self.__rescale_kstar, 0)
            self.__me_2d = hu.RescaleHist(self.__me_2d, self.__rescale_kstar, 0)
            self.__me_1d_reweighted = hu.RescaleHist(
                self.__me_1d_reweighted, self.__rescale_kstar, 0
            )
            self.__cf_reweighted = hu.RescaleHist(
                self.__cf_reweighted, self.__rescale_kstar, 0
            )

        if self.__do_rebin:
            self.__cf_recentered = cu.RescaleGraph(
                self.__cf_recentered, self.__rescale_kstar
            )
            if self.__do_reweight:
                self.__cf_reweighted_recentered = cu.RescaleGraph(
                    self.__cf_reweighted_recentered, self.__rescale_kstar
                )

    def DoRecenter(self):
        """
        Generate TGraphError for correlation function, recentering bins with ME
        """
        logger.debug("Compute recentered correlation function")

        self.__cf_recentered = cu.Recenter(
            self.__cf, self.__handler_noRebin.GetMe1D(), True
        )
        self.__cf_recentered.SetNameTitle(
            self.cf_recentered_name, self.cf_recentered_name
        )
        if self.__do_reweight:
            self.__cf_reweighted_recentered = cu.Recenter(
                self.__cf_reweighted,
                self.__handler_noRebin.GetMe1DRw(),
                True,
            )
            self.__cf_reweighted_recentered.SetNameTitle(
                self.cf_reweighted_recentered_name, self.cf_reweighted_recentered_name
            )

    def ManageHistograms(self):
        """
        Call SetDirectory(0) on all histograms to obtain ownership
        """
        logger.debug("Gain ownership of histograms")
        if self.__inheritsTH1:
            self.__se.SetDirectory(0)
            self.__me.SetDirectory(0)
            if self.__do_ranges:
                self.__se_inRange.SetDirectory(0)
                self.__me_inRange.SetDirectory(0)

        self.__se_2d.SetDirectory(0)
        self.__me_2d.SetDirectory(0)
        self.__se_1d.SetDirectory(0)
        self.__se_1d_normalized.SetDirectory(0)
        self.__me_1d.SetDirectory(0)
        self.__me_1d_normalized.SetDirectory(0)
        self.__cf.SetDirectory(0)

        if self.__do_reweight:
            self.__me_2d_reweighted.SetDirectory(0)
            self.__me_1d_reweighted.SetDirectory(0)
            self.__me_1d_reweighted_normalized.SetDirectory(0)
            self.__cf_reweighted.SetDirectory(0)

    def FinalTouch(self):
        """
        Run all steps to generate correlation function
        """
        logger.debug("Final Touch")
        # rebin kstar axis
        if self.__do_rebin:
            self.DoNoRebin()
            self.DoRebin()
        # apply cuts if defined
        if self.__do_ranges:
            self.DoRanges()
        # do projections onto kstar axis
        self.DoProjections()
        # do projections onto kstar vs reweight axis if reweight axis is defined
        if self.__do_reweight:
            self.DoProjectionsWithRw()
        # if rescale factor is given, rescale kstar axis of all histograms
        # compute correlation function
        self.DoCorrelations()
        # normalize correlation function/distributions
        if self.__do_normalization:
            self.DoNormalizations()
        # if rebin is applied, compute bin center according to mixed event
        if self.__do_rebin:
            self.DoRecenter()
        # rescale kstar axis of 1D/2D histograms and correlation function
        if self.__do_rescale:
            self.DoRescale()
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
        DirSe.Add(self.__se, True)
        if self.__do_ranges:
            DirSe.Add(self.__se_inRange, True)
        if self.__do_reweight:
            DirSe.Add(self.__se_2d, True)
        DirSe.Add(self.__se_1d, True)
        DirSe.Add(self.__se_1d_normalized, True)

        # create subdirectory for mixed event
        DirMe = rt.TDirectoryFile("ME", "ME", "", TdirFile)
        DirMe.Add(self.__me, True)
        if self.__do_ranges:
            DirMe.Add(self.__me_inRange, True)
        if self.__do_reweight:
            DirMe.Add(self.__me_2d, True)
            DirMe.Add(self.__me_2d_reweighted, True)
        DirMe.Add(self.__me_1d, True)
        DirMe.Add(self.__me_1d_normalized, True)
        if self.__do_reweight:
            DirMe.Add(self.__me_1d_reweighted, True)
            DirMe.Add(self.__me_1d_reweighted_normalized, True)

        # create subdirectory for correlation function
        DirCf = rt.TDirectoryFile("CF", "CF", "", TdirFile)
        DirCf.Add(self.__cf, True)
        if self.__do_reweight:
            DirCf.Add(self.__cf_reweighted, True)
        if self.__do_rebin:
            DirCf.Add(self.__cf_recentered, True)
            if self.__do_reweight:
                DirCf.Add(self.__cf_reweighted_recentered, True)

        DirSe.Write()
        DirMe.Write()
        DirCf.Write()

        if self.__do_rebin:
            DirNoRebin = rt.TDirectoryFile("No_Rebin", "No_Rebin", "", TdirFile)
            self.__handler_noRebin.SaveOutputSmall(DirNoRebin)
            DirNoRebin.Write()

    def SaveOutputSmall(self, TdirFile):
        """
        Save output to a TDirectoryFile
        Args:
            TdirFile (TDirectoryFile): TDirectoryFile for saving objects
        """
        logger.debug("Save small output to TDirectoryFile: %s", TdirFile.GetName())
        # create subdirectory for same event
        TdirFile.Add(self.__se_1d, True)
        TdirFile.Add(self.__me_1d, True)
        if self.__do_reweight:
            TdirFile.Add(self.__me_1d_reweighted, True)
        TdirFile.Add(self.__cf, True)

    def GetMe1D(self):
        return self.__me_1d

    def GetMe1DRw(self):
        return self.__me_1d_reweighted
