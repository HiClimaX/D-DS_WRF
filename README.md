# D-DS_WRF

Direct downscaling using WRF

## Authors
* Do Ngoc Khanh - Shibaura Institute of Technology

## Installation

```bash
conda env create -f environment.yml
```

## Workflow

1. Get URLs of CMIP6 6-hour output file.
```bash
./get_6h_url.sh
```

2. Download the 6-hour output files.
```bash
cd 6hr_data
wget -cN -i urls.txt
```

3. Convert to WPS interim format. This script's output is equivalent to ungrib.exe output.
```bash
./to_wps_interim.sh
```

4. Run other WRF steps (`geogrid.exe`, `metgrid.exe`, `real.exe`, and `wrf.exe`) using the data in `wps_interm/`
