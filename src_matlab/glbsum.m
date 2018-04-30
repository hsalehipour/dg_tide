function fsum = glbsum(fld,lat,lon,ind_filter)
% Purpose: Calculate the global sum of a filed fld(lat,lon) on the global 
% sphere accounting for latuitudinal difference in area.
% Note that lat and lon are assumed to be in degrees.
% Optional: also computes global sums by filtering out the data specified
% by ind_filter

% initial setup
earthR = 6.37122e6;
if size(fld,1)==length(lat);
    fld = fld';
end
if nargin>3 
    if size(ind_filter,1)==length(lat);
        ind_filter = ind_filter';
    end
    fld(~ind_filter)=NaN;
end


% Ensure the f field is on a regular grid
nlat = 512;
nlon = 2*nlat;
lat_reg = linspace(90,-90,nlat);
lon_reg = linspace(min(lon),max(lon),nlon);

[LAT    , LON    ]=meshgrid(lat    ,  lon    );
[LAT_reg, LON_reg]=meshgrid(lat_reg,  lon_reg);

% interpolate into regularly spaced grids
fld_reg = interp2(LAT,LON,fld,LAT_reg,LON_reg,'nearest');    


% LAT = repmat(lat',[1 nlon]);
dlat = abs(lat_reg(1)-lat_reg(2));
dlon = abs(lon_reg(2)-lon_reg(1));


% Now compute the surface integral over sphere
W   = cos(pi/180*LAT_reg);
dA = earthR^2*dlat*dlon*(pi/180)^2;           % dA = grid surface area 

ind = ~isnan(fld_reg);
fsum = sum(fld_reg(ind).*W(ind).*dA);



end