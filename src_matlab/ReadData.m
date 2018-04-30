%% Read the ppTide.nc data
fname = [fadrs,'ppTide.nc'];

% % Deprecated: Read the ppTide.nc file on a gaussian lat-lon mesh 
% fname = [fadrs,'ppTide_gauss.nc'];
% gauss = read_ppTide_gauss(fname);
% PIT  = gauss.PIT;       PBL     = gauss.PBL;
% amp    = gauss.amp;     phs     = gauss.phs;
% amp_uu = gauss.amp_uu;  phs_uu  = gauss.phs_uu;
% amp_vv = gauss.amp_vv;  phs_vv  = gauss.phs_vv;
% lat    = gauss.lat;     lon     = gauss.lon;
% depth  = gauss.dep;

% % Read the ppTide.nc file on qtree mesh 
qtree = read_ppTide_qtree(fname);

% Interpolate the results from their native qtree mesh into a gaussian
% lat-lon grid as specified in the interp_fname
[PIT, lat, lon] = tri2latlon(qtree.PIT, interp_fname);
PBL = tri2latlon(qtree.PBL, interp_fname);
amp    = tri2latlon(qtree.amp, interp_fname);  
phs    = tri2latlon(qtree.phs, interp_fname);
amp_uu = tri2latlon(qtree.amp_uu, interp_fname);  
phs_uu = tri2latlon(qtree.phs_uu, interp_fname);
amp_vv = tri2latlon(qtree.amp_vv, interp_fname);  
phs_vv = tri2latlon(qtree.phs_vv, interp_fname);
depth  = tri2latlon(qtree.dep   , interp_fname);

nlat  = length(lat);
nlon  = length(lon);

% In TiDg, the calculated drag terms are already multiplied by g as they
% appear in the Eq (2) of Salehipour et al 2013
PIT = rho/g*PIT;
PBL = rho/g*PBL;

%% Read N2 from CESM1 file
% Read data : QSS: Static Stability, PD: Potential density;

% bf_fname = 'bf/cesm1lgm_yr1590_T85_T_S_PD_Q_RHO.nc';    % lgm
% bf_fname = 'bf/cesm1_fv0.9x1.25_T_S_PD_Q_RHO.nc';       % pre-industrial

ncid = netcdf.open(bf_fname,'NC_NOWRITE');
zlev = netcdf.getVar(ncid,2)/100; 
clat = netcdf.getVar(ncid,3);    
clon = netcdf.getVar(ncid,4);   
PD   = netcdf.getVar(ncid,6);       % Potential density in g/cm^3   
QSS  = netcdf.getVar(ncid,7);       % static stability (drho/dz) in g/cm^4
netcdf.close(ncid);

% Change units into SI
PD = PD*10^3;       % kg/m^3
QSS = QSS*10^5;     % kg/m^4


% compute the buoyancy frequency
cN2 = -g*QSS./PD;    

% adjust the formatting of the data
nzlev = length(zlev);
ind = clon>180;
clon = [clon(ind)-360  ; clon(~ind)];
cN2  = [cN2(ind,:,:)     ;  cN2(~ind,:,:)];

% interpolate from cesm1 grid into our tide model lat-lon grid (zlev)
[cLAT, cLON]=meshgrid(clat, clon);
[ LAT,  LON]=meshgrid( lat,  lon);
N2 = zeros(nlon,nlat,nzlev);
for iz=1:nzlev
    N2(:,:,iz) = interp2(cLAT,cLON,cN2(:,:,iz),LAT,LON);
end

%% Read TPX data
% Note that the files should be already prepared as per instructions in
% "tpxo_input.m".

% tpxo_fadrs = '../data/tpxo7.2/'; 
% 
% % read the elevation, uu and vv files separately
% [lat_tpxo, lon_tpxo, amp_tpxo , phs_tpxo] = ...
%     tpxo_scanf([tpxo_fadrs,'m2_elev.out'],nlon,nlat);
% [~,~, Auu_tpxo, phsuu_tpxo] = tpxo_scanf([tpxo_fadrs,'m2_uu.out'],nlon,nlat);
% [~,~, Avv_tpxo, phsvv_tpxo] = tpxo_scanf([tpxo_fadrs,'m2_vv.out'],nlon,nlat);
% 
% uu_tpxo  = build_fld(Auu_tpxo(:)*1e-2,phsuu_tpxo(:),omega,time);
% vv_tpxo  = build_fld(Avv_tpxo(:)*1e-2,phsvv_tpxo(:),omega,time);
% eta_tpxo  = build_fld(amp_tpxo(:),phs_tpxo(:),omega,time);
% Htot_tpxo = eta_tpxo + depth(:)*ones(1,ntime);
% 
% PE_tpxo = reshape(mean(0.5*rho*g*eta_tpxo.^2,2),nlon,nlat);
% KE_tpxo = reshape(mean(0.5*rho.*Htot_tpxo.*(uu_tpxo.^2+vv_tpxo.^2),2),nlon,nlat);

