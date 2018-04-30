%=====================================================================
% Purpose: This script provides the required lat_lon file for the OTIS
% software package in order to compare the DG tide simulation with TPXO7.2 
% dataset afterwards.
% Usage: for a given TIDG output grid, one may build this file only once.
% Written by: Hesam Salehipour (2012, June 2016).
%=====================================================================
% Running OTISC software afterwards, follow these instructions:
% 1) cp lat_lon ~/tide-workspace/3rd_party/OTIS_netcdf/OTPSnc/
% 2) cd ~/workspace/OTIS_netcdf/OTPSnc/
% 3) Edit "setup.inp" as required based on the description in "README" file
% 4) ./extract_HC<setup.inp 
% 5) rename the sample.out file to your desired name and move it to data/
% 6) use tpxo_scanf.m to read the files and use them
% ...
output = [LAT(:)'; LON(:)'];
filename = 'lat_lon';
fid = fopen(filename, 'wt');
fprintf(fid, '%18.6f \t %18.6f \n', output);
fclose(fid);
