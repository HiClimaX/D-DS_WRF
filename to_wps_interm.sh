#!/bin/bash

output_dir=./wps_interm
mkdir -p "$output_dir"

python to_wps_interm.py \
	--model_conf_path model_conf/MIROC6.csv \
	--start_date 2010-01-05 \
	--end_date 2010-01-12 \
	--interval 6h \
	--prefix "$output_dir/MIROC6_HIST" \
	--nc_files 6hr_data/*MIROC6_historical*.nc

python to_wps_interm.py \
	--model_conf_path model_conf/MIROC6.csv \
	--start_date 2050-01-05 \
	--end_date 2050-01-12 \
	--interval 6h \
	--prefix "$output_dir/MIROC6_SSP370" \
	--nc_files 6hr_data/*MIROC6_ssp370*.nc
