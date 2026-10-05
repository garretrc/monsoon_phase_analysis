# Data licensing

The code license in `LICENSE` does not apply to the NetCDF data.

## CMIP6

The six CMIP6 climatologies derive from NorESM2-MM, AWI-CM-1-1-MR,
CESM2-WACCM, ICON-ESM-LR, IPSL-CM5A2-INCA and NESM3 historical precipitation
simulations. All six sources are currently licensed under
[Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/),
as recorded in the [WCRP model license register](https://wcrp-cmip.github.io/CMIP6_CVs/docs/CMIP6_source_id_licenses.html).
The full license is included in `LICENSE-CC-BY-4.0.txt`.

Retain the source attribution, license link and processing notice when sharing
these data. Source credits and dataset citations are at the end of `README.md`.

## GPCP

The GPCP climatology derives from the NOAA NCEI GPCP Daily Version 1.3 Climate
Data Record. NOAA's [GPCP use agreement](https://www.ncei.noaa.gov/pub/data/sds/cdr/CDRs/Precipitation_GPCP-Daily/UseAgreement_01B-35.pdf)
states that the CDR data are non-proprietary, publicly available and unrestricted
in use. It requests acknowledgment and citation of the data producers and NOAA
CDR Program. Those credits are included in `README.md`. We do not assign a new
Creative Commons license to the underlying GPCP data.

## Processing notice

The supplied files are author-prepared derivatives, not the original daily
archives: 365-day precipitation climatologies for 1997–2014, excluding February
29, expressed in mm/day and retained on native spatial grids. They contain no
smoothing or spatial slicing. The arrays were exported to NetCDF with lossless
compression. The original providers do not endorse this processing or analysis.
Data are supplied without warranty; see the source terms and licenses above.
