%% Read and plot open boundary conditions around SCS.
clear all; close all; clc;
load('../Pu/tidl_force_1h_m2.mat');
% lat = [ss(:,1); nn(:,1); ww(:,1); ee(:,1)];
% lon = [ss(:,2); nn(:,2); ww(:,2); ee(:,2)];
lat = linspace(min(ss(:,1)),max(nn(:,1)),size(ww,1));
lon = linspace(min(ww(:,2)),max(ee(:,2)),size(ss,1));
[LAT, LON] = meshgrid(lat,lon);

ufld = LAT*NaN; vfld=LAT*NaN;
ntime = size(time_1h,1);

for indt=1:ntime;

ufld(:,1)   = usouth(:,indt);       vfld(:,1)   = vsouth(:,indt);
ufld(:,end) = unorth(:,indt);       vfld(:,end) = vnorth(:,indt);
ufld(1,:)   = uwest(:,indt);        vfld(1,:)   = vwest(:,indt);
ufld(end,:) = ueast(:,indt);        vfld(end,:) = veast(:,indt);

indnan = ~isnan(ufld) & ~isnan(vfld);
quiver(LON(indnan)',LAT(indnan)',ufld(indnan)',vfld(indnan)',10,'k','filled')
title(['Time = ', datestr(time_1h(indt,:))])
ylim([round(min(lat-1)) round(max(lat+1))]);
xlim([round(min(lon-1)) round(max(lon+1))]);
% datenum(time_1h)
pause(0.5)

end

%% Compare tidg results with  those obtained directly from OTSC model
load('../Pu/tpxo_z_amp_pha.mat');
tpxo.m2.amp = [APsouth(:,1); APnorth(:,1); APwest(:,1); APeast(:,1)];
tpxo.m2.phs = [APsouth(:,2); APnorth(:,2); APwest(:,2); APeast(:,2)];
load('../Pu/tpxo_u_amp_pha.mat');
tpxo.m2.amp_uu = [APsouth(:,1); APnorth(:,1); APwest(:,1); APeast(:,1)];
tpxo.m2.phs_uu = [APsouth(:,2); APnorth(:,2); APwest(:,2); APeast(:,2)];
load('../Pu/tpxo_v_amp_pha.mat');
tpxo.m2.amp_vv = [APsouth(:,1); APnorth(:,1); APwest(:,1); APeast(:,1)];
tpxo.m2.phs_vv = [APsouth(:,2); APnorth(:,2); APwest(:,2); APeast(:,2)];

tpnlat = length(APwest(:,1));
tpnlon = length(APsouth(:,1));
bndry_ind = [1 tpnlon; tpnlon+1 2*tpnlon;
             2*tpnlon+1  2*tpnlon+tpnlat; ...
             2*tpnlon+tpnlat+1 2*(tpnlon+tpnlat)];
for i=1:4;
    istart = bndry_ind(i,1);
    iend   = bndry_ind(i,2);
    figure(1);  hold all;
    plot(tpxo.m2.amp(istart:iend), tidg.amp(istart:iend),'o');
    figure(2);  hold all;
    plot(tpxo.m2.phs(istart:iend), tidg.phs(istart:iend),'o');
    figure(3);  hold all;
    plot(tpxo.m2.amp_uu(istart:iend), tidg.amp_uu(istart:iend)*100,'o');
    figure(4);  hold all;
    plot(tpxo.m2.phs_uu(istart:iend), tidg.phs_uu(istart:iend),'o');
    figure(5);  hold all;
    plot(tpxo.m2.amp_vv(istart:iend), tidg.amp_vv(istart:iend)*100,'o');
    figure(6);  hold all;
    plot(tpxo.m2.phs_vv(istart:iend), tidg.phs_vv(istart:iend),'o');

end
grid on; box on;
xlabel('$A_\eta$ (m), OTIS', 'interpreter','latex')
ylabel('$A_\eta$ (m), TIDG', 'interpreter','latex')
hndl = legend('south','north','west','east');

%%
 xpatch = [min(lon_in), max(lon_in), max(lon_in), min(lon_in)];
 ypatch = [min(lat_in), min(lat_in), max(lat_in), max(lat_in)];
 patch(xpatch,ypatch,[1 0.4 0.4]);