%% PlotData

figure;
pcolor(lon,lat,amp'); shading flat
% pcolor(lon,lat,sqrt(amp')); shading flat
% caxis([0 sqrt(3)]);
% colorbar('peer',gca,'YTick',sqrt([0 0.5 1 1.8 3]),...
%     'YTickLabel',{'0', '0.5','1','1.8','3'});

% figure;
% phs2 = phs;
% ind = find(phs2<0);
% phs2(ind) = phs2(ind)+360;
% pcolor(lon,lat,phs2'); shading flat

figure;
pcolor(lon,lat,log10(PIT')); shading flat
% caxis([1 4]);
caxis([-4 0]);


%% vertical profile of globally mean quantities
% figure(10); plot(squeeze(nanmean(nanmean(visc_diss))),-zlev); hold all;
% figure(10); plot(squeeze(nanmean(nanmean(Gamma_param_ub))),-zlev); hold all;
% figure(10); plot(squeeze(nanmean(nanmean(N2))),-zlev); hold all;

%% depth-lat contours of Zonally averaged quantities
% figure; pcolor(lat,-zlev,squeeze(log10(nanmean(Reb)))'); shading flat
% figure; pcolor(lat,-zlev,squeeze((nanmean(Gamma_param_lb)))'); shading flat
% figure; pcolor(lat,-zlev,squeeze((nanmean(Gamma_param_ub)))'); shading flat


%% Interpolate the results at a given set of lat lon locations.
% A = loadtxt('lat_lon',2,0);
% lat_in = A(1,:);        lon_in = A(2,:);
% lim_box = 1;
% lat_tri = qtree.lat;
% lon_tri = qtree.lon;
% interp_ind = lat_tri<=max(lat_in) + lim_box & ...
%              lat_tri>=min(lat_in) - lim_box & ...
%              lon_tri<=max(lon_in) + lim_box & ...
%              lon_tri>=min(lon_in) - lim_box;
% % interp_ind1 = lat_tri<=max(lat_in) + lim_box & ...
% %              lat_tri>=min(lat_in) - lim_box & ...
% %              lon_tri<=max(lon_in) + lim_box & ...
% %              lon_tri>=min(lon_in) - lim_box;
% % interp_ind2 = lat_tri<=max(lat_in) - lim_box & ...
% %              lat_tri>=min(lat_in) + lim_box & ...
% %              lon_tri<=max(lon_in) - lim_box & ...
% %              lon_tri>=min(lon_in) + lim_box;
% % interp_ind = interp_ind1 & ~interp_ind2;
% 
% % interpolation method
% method = 'natural';
% 
% % lat-lon
% tidg.lat = lat_in';      tidg.lon = lon_in';
%  
% % elev-amp
% F = scatteredInterpolant(lat_tri(interp_ind),lon_tri(interp_ind),qtree.amp(interp_ind),method);
% tidg.amp = F(lat_in,lon_in)';
% 
% % elev-phase
% F = scatteredInterpolant(lat_tri(interp_ind),lon_tri(interp_ind),qtree.phs(interp_ind),method);       
% tidg.phs = F(lat_in,lon_in)';
% 
% % uu-amp
% F = scatteredInterpolant(lat_tri(interp_ind),lon_tri(interp_ind),qtree.amp_uu(interp_ind),method);    
% tidg.amp_uu = F(lat_in,lon_in)';
% 
% % uu-phase
% F = scatteredInterpolant(lat_tri(interp_ind),lon_tri(interp_ind),qtree.phs_uu(interp_ind),method);    
% tidg.phs_uu = F(lat_in,lon_in)';
% % tidg.phs_uu = smooth(F(lat_in,lon_in),'rlowess')';
% 
% % vv-amp
% F = scatteredInterpolant(lat_tri(interp_ind),lon_tri(interp_ind),qtree.amp_vv(interp_ind),method);    
% tidg.amp_vv = F(lat_in,lon_in)';
% 
% % vv-phase
% F = scatteredInterpolant(lat_tri(interp_ind),lon_tri(interp_ind),qtree.phs_vv(interp_ind),method);    
% tidg.phs_vv = F(lat_in,lon_in)';
% % tidg.phs_vv = smooth(F(lat_in,lon_in),'rlowess')';

% 

%% Outputing E(x,y) field for use in CESM1 on a special grid
% fname_cesm = '../guido/cesm1_tidalmix_fv1.nc';
% ncid = netcdf.open(fname_cesm,'NC_NOWRITE');
% lat_in = netcdf.getVar(ncid,1); 
% lon_in = netcdf.getVar(ncid,2); 
% [LAT_cesm, LON_cesm]=meshgrid(lat_in,  lon_in);
% netcdf.close(ncid);
% 
% % interpolate into cesm grid
% ind = find(LON_cesm>=180);
% LON_cesm(ind) = LON_cesm(ind)-360;
% PITtot_cesm = interp2(LAT,LON,PIT_tot,LAT_cesm,LON_cesm,'nearest');    
% % 
% 
% ncid = netcdf.create('../guido/fromHS/Exy_21ka.nc','NC_WRITE');
% % ncid = netcdf.create('../guido/fromHS/Exy_LGM.nc','NC_WRITE');
% x_dimid = netcdf.defDim(ncid,'nlat',length(lat_in));
% y_dimid = netcdf.defDim(ncid,'nlon',length(lon_in));
% dimid = [y_dimid, x_dimid];
% lat_varID = netcdf.defVar(ncid,'lat','double',x_dimid);
% lon_varID = netcdf.defVar(ncid,'lon','double',y_dimid);
% Exy_varID = netcdf.defVar(ncid,'Exy','double',dimid);
% netcdf.endDef(ncid);
% netcdf.putVar(ncid,lat_varID,lat_in);
% netcdf.putVar(ncid,lon_varID,lon_in);
% netcdf.putVar(ncid,Exy_varID,PITtot_cesm);
% netcdf.close(ncid);

%% test
% figure; pcolor(lon,lat,amp'); shading flat

% tmp00 = amp00; tmp00(isnan(amp00)) = -100;
% tmp120 = amp120; tmp120(isnan(amp120)) = -200;
% lndmsk0_120 = tmp00-tmp120;
% figure; 
% set(gca, 'color', [0.8 0.8 0.8]); 
% hold on; 
% pcolor(lon,lat,amp120');    shading flat; axis equal;

% avg_zonal(topography) = depth(latitude)
% topo_yz = zeros(nlat,1);
% bb = zeros(nlat,1);
% for i=1:nlat
%     ind = ~isnan(depth(:,i)); 
%     topo_yz(i) = mean(depth(ind,i),1); 
%     if ~isempty(find(ind, 1)); 
%         bb(i) = trapz(LON(ind,i),depth(ind,i),1)/(max(LON(ind,i))-min(LON(ind,i))); 
%     end
% end