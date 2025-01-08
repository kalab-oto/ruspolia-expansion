**Estimation of the recent expansion rate of _Ruspolia nitidula_ (Orthoptera) on a regional and landscape scale**

Data and scripts for analysis in https://doi.org/10.3390/insects12070639

Repository includes code for the data downloads ([`scripts/data_download`](https://github.com/kalab-oto/ruspolia-wip/tree/master/scripts/data_download)), complete analysis ([`scripts/`](https://github.com/kalab-oto/ruspolia-wip/tree/master/scripts)), generating plots [`scripts/6_plot.R`](https://github.com/kalab-oto/ruspolia-wip/blob/b838c3853a68eef03b9e6c311973f4df238f36b5/scripts/6_plot.R) and maps from QGIS [`scripts/7_qgis_maps.py`](https://github.com/kalab-oto/ruspolia-wip/blob/b838c3853a68eef03b9e6c311973f4df238f36b5/scripts/7_qgis_maps) (QGIS project is included in [`qgis/`](https://github.com/kalab-oto/ruspolia-wip/tree/master/qgis) directory).

### Maps
To generate complete maps, it is necessary to download and preprocess data that cannot be provided in this repository:

1. Download SRTM data (part of the [`scripts/data_download/env_download.R`](https://github.com/kalab-oto/ruspolia-wip/blob/master/scripts/data_download/env_download.R))

2. Download NDOP data ([`scripts/data_download/occ_download.R`](https://github.com/kalab-oto/ruspolia-wip/blob/master/scripts/data_download/occ_download.R))
[`scripts/1_occ_preprocess.R`](https://github.com/kalab-oto/ruspolia-wip/blob/master/scripts/1_occ_preprocess.R))

3. Preprocess occurrence data ([`scripts/2_env_data_preprocess.R`](https://github.com/kalab-oto/ruspolia-wip/blob/master/scripts/2_env_data_preprocess.R))
