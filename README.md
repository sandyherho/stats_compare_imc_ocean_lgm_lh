# Supplementary Materials: "Statistical Characterization of LGM and Late Holocene Oceanography in the Indonesian Maritime Continent"


[![No Maintenance Intended](http://unmaintained.tech/badge.svg)](http://unmaintained.tech/)
[![License: WTFPL](https://img.shields.io/badge/License-WTFPL-brightgreen.svg)](http://www.wtfpl.net/about/)
[![Climate Dynamics](https://img.shields.io/badge/Climate-Dynamics-orange.svg)](https://link.springer.com/journal/382)

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![NumPy](https://img.shields.io/badge/numpy-%23013243.svg?logo=numpy&logoColor=white)](https://numpy.org/)
[![SciPy](https://img.shields.io/badge/SciPy-%230C55A5.svg?logo=scipy&logoColor=white)](https://scipy.org/)
[![Pandas](https://img.shields.io/badge/pandas-%23150458.svg?logo=pandas&logoColor=white)](https://pandas.pydata.org/)
[![Matplotlib](https://img.shields.io/badge/Matplotlib-%23ffffff.svg?logo=Matplotlib&logoColor=black)](https://matplotlib.org/)
[![xarray](https://img.shields.io/badge/xarray-blue.svg)](https://xarray.dev/)

This repository contains the supplementary materials, analysis code, and statistical results for the paper:

**Mazdhania, D.Z., Herho, S.H.S., Anwar, I.P., Khadami, F., and Rachmayani, R. (2025).** Statistical Characterization of LGM and Late Holocene Oceanography in the Indonesian Maritime Continent. *[xxxxx]*, *[xxxx]*, xxx.


## Overview

This study provides a comprehensive statistical characterization of surface oceanography in the Indonesian Maritime Continent (IMC) during two distinct paleoclimate epochs:
- **Last Glacial Maximum (LGM)**: 23–19 ka
- **Late Holocene (LH)**: 4–0 ka

Using the lgmDA v2.1 paleoclimate data assimilation product, we analyze:
- Sea surface temperature (SST)
- Sea surface salinity (SSS)
- Seawater oxygen isotope composition ( $\delta^{18}$ Osw)

## Repository Structure

```
.
├── script/
│   ├── map_plot.py          # Spatial distribution visualization
│   ├── stats_plot.py        # KDE and boxplot generation
│   └── stats_test.py        # Comprehensive statistical testing
├── stats/
│   ├── complete_statistical_test_report.txt
│   ├── sst_statistical_analysis.txt
│   ├── sss_statistical_analysis.txt
│   ├── d18osw_statistical_analysis.txt
│   └── map_stats.txt
├── figs/                    # Generated figures
├── raw_data/               # lgmDA v2.1 NetCDF files
├── processed_data/         # CSV outputs
└── LICENSE                 # WTFPL License
```

## Key Findings

### Temperature (SST)
- **LGM cooling**: 2.8°C cooler than Late Holocene
- **Effect size**: Cliff's δ = -0.993 (large)
- **Spatial pattern**: Strongest cooling in eastern warm pool

### Salinity (SSS)
- **LGM salinity**: 0.88 PSU higher than Late Holocene
- **Effect size**: Cliff's $\delta$ = 0.469 (medium)
- **Spatial variability**: LGM shows reduced heterogeneity

### Oxygen Isotopes ( $\delta^{18}$ Osw)
- **LGM enrichment**: 1.0‰ higher than Late Holocene
- **Effect size**: Cliff's δ = 0.994 (large)
- **Pattern**: Ubiquitous enrichment across domain

## Statistical Methods

### Normality Assessment (5 tests)
- Shapiro-Wilk Test
- Anderson-Darling Test
- Kolmogorov-Smirnov Test
- D'Agostino-Pearson Omnibus Test
- Jarque-Bera Test

**Result**: All variables in both epochs show strong non-normality

### Nonparametric Hypothesis Testing (6 tests)
- Mann-Whitney U Test
- Wilcoxon Rank-Sum Test
- Kruskal-Wallis H Test
- Mood's Median Test
- Kolmogorov-Smirnov Two-Sample Test
- Epps-Singleton Test

**Result**: All tests show $p < 0.001$, indicating strong evidence for distributional differences

### Information-Theoretic Measures
- Shannon Entropy (distributional complexity)
- Approximate Entropy (regularity)
- Sample Entropy (pattern regularity)
- Coefficient of Variation (relative variability)

## Usage

### Prerequisites

```bash
pip install numpy scipy pandas xarray matplotlib
```

### Running the Analysis

1. **Generate spatial statistics and maps**:
```bash
python script/map_plot.py
```

2. **Create KDE and boxplot visualizations**:
```bash
python script/stats_plot.py
```

3. **Perform comprehensive statistical testing**:
```bash
python script/stats_test.py
```


## Data Source

This analysis uses the **Last Glacial Maximum Data Assimilation version 2.1 (lgmDA v2.1)**:

- **Reference**: Tierney et al. (2020a, 2020b); Tierney (2022)
- **DOI**: [10.5281/zenodo.5171432](https://doi.org/10.5281/zenodo.5171432)
- **License**: GNU General Public License

The lgmDA product synthesizes marine geochemical proxy data with isotope-enabled climate model simulations through offline ensemble data assimilation.

## Authors

- **Daffa Z. Mazdhania**
- **Sandy H. S. Herho**
- **Iwan P. Anwar**
- **Faruq Khadami**
- **Rima Rachmayani**


## Citation

```bibtex
@article{mazdhania202xstatistical,
  title={Statistical Characterization of LGM and Late Holocene Oceanography in the Indonesian Maritime Continent},
  author={Mazdhania, Daffa Z. and Herho, Sandy H. S. and Anwar, Iwan P. and Khadami, Faruq and Rachmayani, Rima},
  journal={xxxx},
  volume={xxxx},
  pages={xxx},
  year={202x},
  publisher={xxxxx}
}
```


## Acknowledgements

This study was supported by:
- Dean's Distinguished Fellowship from the College of Natural and Agricultural Sciences, University of California, Riverside (2023) awarded to S.H.S.H.
- Bandung Institute of Technology Research, Community Service and Innovation Program (PPMI-ITB 2025) awarded to F.K. and I.P.A.

## License

This project is licensed under the **WTFPL (Do What The Fuck You Want To Public License)** - see the [LICENSE](LICENSE) file for details.

**Note**: This repository contains supplementary materials for a published research article. The code is provided as-is for reproducibility and transparency purposes. No active maintenance is intended after publication.
