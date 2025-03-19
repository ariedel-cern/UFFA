import FemtoAnalysis as FA
import FemtoDreamReader as FDR
#import sys

#   Current options include:
#           "function":     'cf', 'syst', 'tf' -> for correlation function, systematics, template fits
#           "pair":         'pp', 'pl' -> for q&a plots relevant for individual analyses
#           "path":         "string" -> full path to the root file, might include ~/ for home directory
#           "file":         "string" -> name of the root file
#           "fullpath":     "string" -> full path and file name equal to "path" + "file"
#           "outDir":       "string" -> output directory
#           "rename":       "string" -> rename output file
#           "fileTDir":     "string" -> root file directory: path to directory inside the root file
#           "nameSE":       "string" -> path + name of the se plot inside the provided "fileTDir" if given
#           "nameME":       "string" -> same as nameSE but for the ME distribution
#           "newfile":      'new', 'recreate', 'update' -> same option as in ROOT, 'new' will rename if file already exists
#           "mc":           'true', 'false' -> save monte carlo data from provided root file
#           "mcTDir":       "string" -> root file directory for the monte carlo data
#           "bins":         [list of floats] -> binning for differential analysis
#           "diff3d":       'mt', 'mult' -> project 3D plots first in mt/mult 2D and after in mult/mt 1D or vice versa
#           "bins3d":       [list of floats] -> binning for 3D plots to 2D plots
#           "yield":        [GeV, Deviation] -> integrated analysis: include systematics inside deviation for the GeV range
#           "rebin":        int or [list of ints] -> rebin output plots
#                               for tf: int -> all dca/cpa rebinned with int
#                               or      [list of ints] -> each int will correspond to one range of the binning if provided
#           "atype":        'int', 'dif' -> integrated analysis or differential analysis
#           "htype":        'k', 'mt', 'mult', 'mt3d', 'mult3d', 'mtmult', 'rew3d'
#                               'k'     -> k* - relative pair momentum distribution
#                               'mt'  -> mt vs k* distribution
#                               'mult'  -> multiplicity vs k* distribution
#                               'mt3d' -> mt vs k* from 3D distribution but integrated in mult
#                               'mult3d' -> mult vs k* from 3D distribution but integrated in mt
#                               'mtmult' -> mult vs mt vs k* 3D distribution
#                               'rew3d' -> mult vs mt vs k* 3D differentially in mt and reweighted in mult
#           "tftype":       'dca', 'cpa' -> option for the template fit plots
#           "templates":    [list of th1 plots] -> list of dca/cpa plots for fitting
#           "namelist":     [list of strings] -> names of dca/cpa plots for fitting
#           "fitrange":     float -> fitrange for the template fitter
#           "normalize":    [float, float] -> normalization range for the correlation function
#           "include":      "string" or [list of strings] -> include these variations in the systematics
#           "exclude":      "string" or [list of strings] -> exclude these variations in the systematics
#           "interactive":  'True', 'False' -> include/exclude interactively variations in terminal
#           "debug":        'True', 'False' -> debug information in console

ipath = ""
opath = ""

Yield = [0.2,0.3]
Rebin = [2, 5, 10]

mtBins = [0, 1.2, 1.56, 4.5]
mtBins = [0,1.02,1.14,1.2,1.26,1.38,1.56,1.86,4.5]
multBins = [0, 11, 20, 200]
multBins = [0, 200]

#filename = sys.argv[1]
filename = "root_input/pL_22all.root"
filename = "AnalysisResults_sqm.root"
filename = "maps.root"

templates = ""

# f1 = FDR.FemtoDreamReader(templates, "_base-tempFit")
# mc_plots = []
# mc_plots.append(f1.getHisto("Tracks_one_MC/hDCAxy_Primary"))
# mc_plots.append(f1.getHisto("Tracks_one_MC/hDCAxy_DaughterLambda"))
# mc_plots.append(f1.getHisto("Tracks_one_MC/hDCAxy_DaughterSigmaplus"))
# mc_plots.append(f1.getHisto("Tracks_one_MC/hDCAxy_Fake"))

namelist = ['prim', 'lam', 'sig', 'fake']

settings_cf = {
        "function":     'cf',
        "pair":         'pp',
        "file":         filename,
        "fileTDir":     "femto-dream-pair-task-track-track",
        #"nameSE":       "SameEvent/relPairDist",
        #"nameME":       "MixedEvent/relPairDist",
        "newfile":      "new",
        "outDir":       opath,
        "mc":           None,
        "mcTDir":       "",
        "rename":       "test-limits",
        #"type":         ['dif', 'rew4d', 'mt'],
        "type":         ['dif', '4d', 'mult'],
        #"atype":        'dif',
        #"htype":        'rew4d',
        #"diff3d":       'mt',
        "bins3d":       multBins,
        "bins":         mtBins,
        #"bins3d":       mtBins,
        #"bins":         multBins,
        "rebin":        [2, 5, 10],
        "rewrange":     [0, 1.0],
        "percentile":   [0, 20],
        "normalize":    [0.24, 0.34],
        "debug":        True
    }

settings_tf = {
        "function":     'tf',
        "file":         ipath + filename,
        "fileTDir":     "_base-tempFit",
        "newfile":      "new",
        "templates":    templates,
        #"templates":    mc_plots,
        "mcTDir":       "_base-tempFit",
        "outDir":       opath,
        "rename":       None,
        "bins":         [0.500, 0.678, 0.856, 1.033, 1.211, 1.388, 1.567, 1.743, 1.921, 2.099, 2.453, 2.986, 4.051],
        "rebin":        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 4],
        "tftype":       'dca',
        "namelist":     namelist,
        "fitrange":     0.3,
        "debug":        False
        }

settings_syst = {
        "function":     'syst',
        "file":         filename,
        "fileTDir":     "",
        "newfile":      "new",
        "outDir":       opath,
        "mc":           None,
        "mcTDir":       "",
        #"rename":       "test.root",
        "atype":        'dif',
        "htype":        'mtmult',
        "diff3d":       'mt',
        "bins3d":       mtBins,
        "bins":         multBins,
        #"yield":        [0.3, 0.2],
        #"atype":        "dif",
        #"htype":        "kmult",
        #"bins":         [0, 20, 40, 60],
        #"bins":         [0.5, 1.2, 2.5, 4],
        #"rebin":        2,
        "rebin":        [2, 5, 10],
        "yield":        [0.2, 0.6],
        "normalize":    [0.24, 0.34],
        #"exclude":      '',
        #"interactive":  True,
        "debug":        True
    }

settings_syst2 = {
        "function":     'syst',
        "file":         filename,
        "fileTDir":     "",
        "newfile":      "recreate",
        "outDir":       opath,
        "rename":       "pp_mtDiff.root",
        "atype":        'dif',
        "htype":        'mt',
        "bins":          [0,1.02,1.26,1.86,4.5],
        "rebin":         2,
        "debug":         True
    }

FA.UFFA(settings_cf)

