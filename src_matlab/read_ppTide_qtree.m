function data = read_ppTide_qtree(filename)

% function [PIT, PBL, amp, phs, lat, lon, dep] = read_ppTide(filename)
% Purpose  : Read the pp'ed NetCDF output file ppTide.nc

% open netCDF file
ncid = netcdf.open(filename,'NC_NOWRITE');

% inquire the variable info within .nc file
% [numdims, numvars, numglobalatts, unlimdimID] = netcdf.inq(ncid);
% [varname, xtype, varDimIDs, varAtts] = netcdf.inqVar(ncid,4);

% Read data
data.lat = netcdf.getVar(ncid,0,'double');
data.lon = netcdf.getVar(ncid,1,'double');
data.dep = netcdf.getVar(ncid,2,'double');
data.amp = netcdf.getVar(ncid,3,'double');
data.phs = netcdf.getVar(ncid,4,'double');
data.amp_uu = netcdf.getVar(ncid,5,'double');
data.phs_uu = netcdf.getVar(ncid,6,'double');
data.amp_vv = netcdf.getVar(ncid,7,'double');
data.phs_vv = netcdf.getVar(ncid,8,'double');
data.PIT = netcdf.getVar(ncid,9,'double');
data.PBL = netcdf.getVar(ncid,10,'double');

end



