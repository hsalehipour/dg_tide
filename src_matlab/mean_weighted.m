function favg = mean_weighted(lat,f)
% Purpose: Calculate the weighted averaged on the global sphere accounting
% for latuitudinal difference in area.

nlon = size(f,2);
LAT = repmat(lat',[1 nlon]);
W   = cos(pi/180*LAT);
ind = ~isnan(f);
favg = sum(f(ind).*W(ind)) / sum(W(ind));

end