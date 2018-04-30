function [PIT, PBL, amp, phs, lat, lon, dep] = read_ppTide(filename)

% function [PIT, PBL, amp, phs, lat, lon, dep] = read_ppTide(filename)
% Purpose  : Read the pp'ed NetCDF output file ppTide.nc

% open netCDF file
ncid = netcdf.open(filename,'NC_NOWRITE');

% inquire the variable info within .nc file
% [numdims, numvars, numglobalatts, unlimdimID] = netcdf.inq(ncid);
% [varname, xtype, varDimIDs, varAtts] = netcdf.inqVar(ncid,4);

% Read data
lat = netcdf.getVar(ncid,0,'double')';
lon = netcdf.getVar(ncid,1,'double')';
dep = netcdf.getVar(ncid,2,'double')';
amp = netcdf.getVar(ncid,3,'double')';
phs = netcdf.getVar(ncid,4,'double')';
PIT = netcdf.getVar(ncid,5,'double')';
PBL = netcdf.getVar(ncid,6,'double')';

end



