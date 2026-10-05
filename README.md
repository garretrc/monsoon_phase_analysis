# Monsoon phase analysis

This example reproduces the Indian Summer Monsoon phase analysis in
"Sliced Elastic Distance for Evaluating Amplitude and Phase Differences in
Precipitation Models," including Figures 6–8 in the manuscript.

## Data

`data/` contains unsmoothed daily precipitation climatologies for
GPCP and six CMIP6 models: NorESM2-MM, AWI-CM-1-1-MR, CESM2-WACCM,
ICON-ESM-LR, IPSL-CM5A2-INCA and NESM3. These seven NetCDF files are the only
required input data. The native spatial grid is provided for each data product.

## Running the example

The scripts require R and the packages `tidyverse`, `Matrix`, `ncdf4`, `glmgen`,
`FNN`, `fdasrvf`, `viridis`, `ggrepel`, `patchwork`, `sf`, `maps`, `rworldmap`
and `rworldxtra`.

The input climatologies are included in data/ using Git LFS. Install Git LFS before cloning this repository:

```sh
git lfs install
git clone https://github.com/garretrc/monsoon_phase_analysis.git
cd monsoon_phase_analysis
```

If you already cloned the repository, run `git lfs install` followed by `git lfs pull` from the repository directory.

From the `monsoon_phase_analysis/` directory, run:

```sh
Rscript scripts/run.R
```

## Scripts

- `scripts/run.R` sets the parameters and runs the analysis.
- `scripts/slicing.R` smooths the native climatologies, remaps by nearest
  neighbour and computes convolutional slices at the monsoon locations.
- `scripts/monsoon.R` determines GPCP onset and retreat dates, aligns the model
  slices to GPCP and calculates onset and retreat timing biases.
- `scripts/plotting.R` creates the monsoon figures.

The default parameters reproduce the manuscript analysis.

## Output

PDF and PNG figures are saved in `output/`:

| Figure | Files |
|---|---|
| 6: GPCP onset and retreat dates | `onset_retreat.pdf`, `onset_retreat.png` |
| 7: Phase alignment at 22.5°N, 78°E | `model_slices.pdf`, `model_slices.png` |
| 8: Onset and retreat timing biases | `onset_bias.pdf`, `onset_bias.png`, `retreat_bias.pdf`, `retreat_bias.png` |

## License and acknowledgments

The analysis code is distributed under the GNU General Public License,
version 3 (GPL-3.0-only); see [LICENSE](LICENSE). It is provided without warranty.
The NetCDF data retain their separate source terms; see
[DATA_LICENSE.md](DATA_LICENSE.md) and [CC BY 4.0](LICENSE-CC-BY-4.0.txt).

The `ndims` and `trapz` functions in `scripts/monsoon.R` are copied from
**J. Derek Tucker's fdasrvf 2.2.0**. Citation: Tucker, J. D. (2024).
*fdasrvf: Elastic Functional Data Analysis*, R package version 2.2.0.
[Package and source](https://CRAN.R-project.org/package=fdasrvf).

We acknowledge the World Climate Research Programme's Working Group on Coupled
Modelling, which coordinated CMIP6, the modeling groups listed below for producing
and sharing their output, and the Earth System Grid Federation for archiving and
providing access. The six CMIP6 datasets are licensed under CC BY 4.0.

| Model | Source credit and CMIP6 dataset citation |
|---|---|
| NorESM2-MM | Bentsen, M., et al. (2019), NCC. [CMIP6.506](https://doi.org/10.22033/ESGF/CMIP6.506) |
| AWI-CM-1-1-MR | Semmler, T., et al. (2018), AWI. [CMIP6.359](https://doi.org/10.22033/ESGF/CMIP6.359) |
| CESM2-WACCM | Danabasoglu, G. (2019), NCAR. [CMIP6.10024](https://doi.org/10.22033/ESGF/CMIP6.10024) |
| ICON-ESM-LR | Lorenz, S., et al. (2021), MPI-M. [CMIP6.743](https://doi.org/10.22033/ESGF/CMIP6.743) |
| IPSL-CM5A2-INCA | Boucher, O., et al. (2020), IPSL. [CMIP6.13642](https://doi.org/10.22033/ESGF/CMIP6.13642) |
| NESM3 | Cao, J., and Wang, B. (2019), NUIST. [CMIP6.2021](https://doi.org/10.22033/ESGF/CMIP6.2021) |

GPCP Daily Version 1.3 was obtained from NOAA NCEI and developed by Robert Adler,
Jian-Jian Wang and Mathew Sapiano at the University of Maryland. Please cite:

- Adler, R., Wang, J.-J., Sapiano, M., Huffman, G., Bolvin, D., Nelkin, E., and
  NOAA CDR Program (2017). *Global Precipitation Climatology Project (GPCP)
  Climate Data Record, Version 1.3 (Daily)*. NOAA NCEI.
  [doi:10.7289/V5RX998Z](https://doi.org/10.7289/V5RX998Z).
  The subset used here covers 1997–2014.
- Huffman, G. J., et al. (2001). Global Precipitation at One-Degree Daily
  Resolution from Multisatellite Observations. *Journal of Hydrometeorology*,
  2(1), 36–50.
