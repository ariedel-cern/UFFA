# 🌀 UFFA — Unified Framework for Femtoscopy Analysis

UFFA (Unified Framework for Femtoscopy Analysis) is a modular Python-based framework designed for **femtoscopic analysis**.
It automates the workflow from histogram input to correlation function output, supporting rebinning, range selection, normalization, and result merging — all output is stored in a structured ROOT file.

UFFA uses **ROOT (PyROOT)** for histogram manipulation and supports both **sequential** and **parallel** execution.

---

## 📦 Installation

UFFA requires **Python ≥ 3.9** and **ROOT with PyROOT** enabled. You have two ways to install it:

### Install via pip (easiest)

If you have a Python environment with ROOT installed (e.g., an O2Physics environment), you can install UFFA directly:

```bash
# Activate your o2physics environment first
alienv enter O2Physics/latest

# Clone the repository
git clone https://github.com/ariedel-cern/UFFA.git
cd UFFA

# install via pip
pip install .
```
Then you can simply import it in your scripts, for example

```python
from UFFA import AnalysisHandler
```

### Use sys.path (Recommended for development)

If you cloned the repository locally, you can add it to Python path without installing:

```python
import sys
sys.path.append("/path/to/your/UFFA")
from UFFA import AnalysisHandler
```

Use this method if you want to modify the source code directly.

## ⚙️ Example Usage

```python
from UFFA import AnalysisHandler

analysis_config = {
    "Input_File": "./input.root", # path to input file
    "Path_SE": "femto/SE", # path to same event distribution in Input_File
    "Path_ME": "femto/ME", # path to mixed event distribution in Input_File
    "Output_File": "output/output.root", # path to output file (will be overwritten)
    "Output_Dir": "analysis", # name of TDirecotryFile to stored results in Output_File
    "Normalization_Range": (0.24, 0.34),
    "Rebin_Factor_Kstar_Axis": [1, 2, 4, 8],
    "Index_Kstar_Axis": 0, # index of kstar axis in same and mixed event distribution
    "Index_Reweight_Axis": 2, # index of reweight axis (should be multiplicity) in same and mixed event distribution
    "Bins": [[], [1.4, 1.65, 2.7], [], []],
    # in this example,
    # the 0th dimension is kstar, no cuts are applied so we pass an empty list
    # the 1st dimension is transverse mass, here we want to 2 bins so we specify three limits, -> [1.4,1.65] & [1.65,2.7]
    # the 2st dimension is multiplicity, no cuts are applied so we pass empty list
    # the 3st dimension is multiplicity percentile, no cuts are applied so we pass empty list
}

Handler = AnalysisHandler(analysis_config)
Handler.SteerAnalysis(parallel=True)
```

📁 Output Structure

UFFA stores all analysis results in the directory defined by Output_Dir:

Example output:

```bash
> rootls output/output.root:analysis
Rebin_1_Dim_1-1.4-1.65
Rebin_1_Dim_1-1.65-2.7
Rebin_2_Dim_1-1.4-1.65
Rebin_2_Dim_1-1.65-2.7
Rebin_4_Dim_1-1.4-1.65
Rebin_4_Dim_1-1.65-2.7
Rebin_8_Dim_1-1.4-1.65
Rebin_8_Dim_1-1.65-2.7
```

Each directory corresponds to a single configuration, identified by:

- the rebin factor (Rebin_\<Rebin factor\>)
- the range (Dim_\<ith dimension\>_\<Lower Limit\>-\<Upper Limit\>)

The ranges correspond typically to bins in transverse mass, invariant mass, multiplicity, multiplicity percentile, ...

Each configuration directory contains:
```
SE/   → same-event histograms
ME/   → mixed-event histograms
CF/   → correlation functions (SE/ME)
```
The user can inspect the output thoroughly after execution.
