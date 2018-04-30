function [avgBVF,BVF,lat,lon] = BuoyancyFreq(filename)

% function [avgBVF,BVF,lat,lon] = BuoyancyFreq(filename)
% Purpose: Calculate the Brunt–Väisälä frequency (or N buoyancy frequency)

% Read data : QSS: Static Stability, PD: Potential density;
ncid = netcdf.open(filename,'NC_NOWRITE');
z_t  = netcdf.getVar(ncid,2); 
lat  = netcdf.getVar(ncid,3);    
lon  = netcdf.getVar(ncid,4);   
PD   = netcdf.getVar(ncid,6);       % Potential density in g/cm^3   
QSS  = netcdf.getVar(ncid,7);       % static stability (drho/dz) in g/cm^4

% Change units into SI
PD = PD*10^3;       % kg/m^3
QSS = QSS*10^5;     % kg/m^4


% compute the frequency
g  = 9.80616;
N2 = -g*QSS./PD;    
BVF = sqrt(abs(N2));


% Compute the vertical average avgBVF = 1/H * int_{-H}^0 BVF dz
nlon = size(PD,1);       
nlat = size(PD,2);
avgBVF = nan(nlon,nlat);
for ilat=1:nlat;
    for ilon=1:nlon;
        tBVF = squeeze(BVF(ilon,ilat,:));
        ind  = find(~isnan(tBVF));
        if ~isempty(ind)
            %plot(tBVF(ids),-z_t(ids),'-ob'); drawnow;
            tBVF = trapz(z_t(ind),tBVF(ind));
            maxH = max(z_t(ind))-min(z_t(ind));            
            avgBVF(ilon,ilat) = tBVF/maxH;
        end
    end
end
% contourf(lon,lat,avgBVF');
