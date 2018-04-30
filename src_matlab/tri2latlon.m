function [phi_latlon, lat, lon] = tri2latlon(phi, interp_file)


% open netCDF file
ncid = netcdf.open(interp_file,'NC_NOWRITE');

lat = netcdf.getVar(ncid,0,'double');
lon = netcdf.getVar(ncid,1,'double');
icell = netcdf.getVar(ncid,2);
t1 = netcdf.getVar(ncid,3);
t2 = netcdf.getVar(ncid,4);
t3 = 1-t1-t2;

Np = size(phi,1);       % N = (N+1)*(N+2)/2;
N = (sqrt(8*Np+1)-3)/2;
Nfp = N+1;

phi_va = phi(1,:)';
phi_vb = phi(Nfp,:)';
phi_vc = phi(Np,:)';

phi_latlon = NaN(size(t1,1), size(t1,2));
ids = find(icell~=0);
phi_latlon(ids) =   t1(ids).*phi_va(icell(ids)) + ...
                    t2(ids).*phi_vb(icell(ids)) + ...
                    t3(ids).*phi_vc(icell(ids));
phi_latlon = double(phi_latlon);

end