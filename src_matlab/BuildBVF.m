% BuildBVF(output_flag, src_filename, bvf_filename)
% Purpose: Build global depth-averaged Brunt–Väisälä frequency matrix from
%          CESM1 computations

clear all; clc;

src_filename = 'bf/cesm1_fv0.9x1.25_T_S_PD_Q_RHO.nc';   % pre-industrial
bvf_filename = 'bf/bvf_pi.nc';                          % pre-industrial

% src_filename = 'bf/cesm1lgm_yr1590_T85_T_S_PD_Q_RHO.nc';    % lgm
% bvf_filename = 'bf/bvf_lgm.nc';                             % lgm

% BVF: Brunt–Väisälä frequency (N)% -- Read data : QSS: Static Stability
[outBVF,~,lat,lon] = BuoyancyFreq(src_filename);

lon2   = [lon(lon>180)-360  ; lon(lon<=180)];
outBVF2= [outBVF(lon>180,:) ; outBVF(lon<=180,:)];

nlat = length(lat);
nlon = length(lon2);
icell = ones(nlon,nlat);    
icell(isnan(outBVF2)) = 0;     

% Write the BF into a NetCDF file
% Open netCDF file.
ncid = netcdf.create(bvf_filename,'NC_WRITE');

% Define the dimensions of the variable.
x_dimid = netcdf.defDim(ncid,'nlat',nlat);
y_dimid = netcdf.defDim(ncid,'nlon',nlon);
dimid = [y_dimid, x_dimid];

% Define a new variable in the file.
bf_varID  = netcdf.defVar(ncid,'bf','double',dimid);
lat_varID = netcdf.defVar(ncid,'lat','double',x_dimid);
lon_varID = netcdf.defVar(ncid,'lon','double',y_dimid);
 ic_varID = netcdf.defVar(ncid,'icell','double',dimid);

% Leave define mode and enter data mode to write data.
netcdf.endDef(ncid);

% Write data to variable.
netcdf.putVar(ncid,lat_varID,lat);
netcdf.putVar(ncid,lon_varID,lon2);
netcdf.putVar(ncid, bf_varID,outBVF2);
netcdf.putVar(ncid, ic_varID,icell );

% Verify that the variable was created.
[varname xtype dimid natts ] = netcdf.inqVar(ncid,0)

% close the file
netcdf.close(ncid);


