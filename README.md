# StableClim

## Code used for creating and validating the *StableClim* database.

This repository contains the code used to create and validate the *StableClim* database. The code used in creating the database is presented in a series of PDF documents (generated using R markdown), whilst validation code is available in the .R documents. Specific details regarding the climatology windows etc. are available in the main manuscript.

The CMIP5 data used in the database can be downloaded from the ESGF system using the scripts available [here](https://github.com/GlobalEcologyLab/ESGF_ClimateDownloads). TraCE-21ka data is available through [PaleoView](https://github.com/GlobalEcologyLab/PaleoView).

The code presented here can be used to generate, for example, ensemble median (ensemble mean is used in *StableClim*) estimates of climate under four different RCP scenarios. The code can also be altered to generate the same data at a finer spatial resolution.

Files are structured as:

01 = Pre-processing code

02 = Regression code

03 = Validation code

04 = custom functions used in validation

*StableClim* outputs are available through [figshare](https://doi.org/10.25909/5ea59831121bc)

*Notes*
The 03_Validation_PRIMER_SNR_Comparison_StableClim.tar.gz file contains a [PRIMER-6](https://www.primer-e.com/) workspace containing validation data and PERMDISP and PERMANOVA analysis results.
