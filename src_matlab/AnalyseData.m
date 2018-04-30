%% Analyse Data



%% Compute PE and KE
% Note: pptide subroutine calculates the AMP and PHS based on 
%  e.g. fld=amp.*cos(omega*time-phs);

ntime = 100;
time = linspace(0,2*pi/omega,ntime);    % one tidal cycle

eta = build_fld(amp(:),phs(:),omega,time);
uu  = build_fld(amp_uu(:),phs_uu(:),omega,time);
vv  = build_fld(amp_vv(:),phs_vv(:),omega,time);

Htot = eta + depth(:)*ones(1,ntime);

PE=reshape(mean(0.5*rho*g*eta.^2,2),nlon,nlat);
KE=reshape(mean(0.5*rho.*Htot.*(uu.^2+vv.^2),2),nlon,nlat);


%% Analyse PIT : global sum, etc

lim_depth = 1000;   
lim_lat   = 66;

% indeces of areas to be filtered out
ind_dw = depth>=lim_depth;
ind_sw = ~ind_dw;
ind_lat= LAT<=lim_lat & LAT>=-lim_lat;

% Dissipation due to internal tide and bottom drag
PIT_sum(ic) = glbsum(PIT,lat,lon)/1e12;    % Terra Watts
PBL_sum(ic) = glbsum(PBL,lat,lon)/1e12;    % Terra Watts

% Shallow water and deep ocean dissipation
Ptot = PIT+PBL;
Psw_sum(ic) = glbsum(Ptot,lat,lon, ind_sw)/1e12;    % Terra Watts
Pdw_sum(ic) = glbsum(Ptot,lat,lon, ind_dw)/1e12;    % Terra Watts
Ptot_sum = Pdw_sum + Psw_sum;


% add dissipation due to multiple constituents
% if ic==1; PIT_tot=PIT*0; PBL_tot=PBL*0; end;
% PIT_tot = PIT_tot + PIT;
% PBL_tot = PBL_tot + PBL;
% Pdw = PIT*0;        Psw = PIT*0;
% Pdw(ind_dw) = PIT_tot(ind_dw)+PBL_tot(ind_dw);
% Psw(ind_sw) = PIT_tot(ind_sw)+PBL_tot(ind_sw);


% Total kinetic energy and potential energy
KE_sum(ic) = glbsum(KE,lat,lon)/1e17;    % 10^17 Watts
PE_sum(ic) = glbsum(PE,lat,lon)/1e17;    % 10^17 Watts



%% Find diapycnal diffusivity 

ksi = 500;
q   = 1/3;
flux_coeff = 0.2;

Fz  = NaN(nlon,nlat,nzlev);
Exy = Fz*0;

for i=1:nzlev;  
    Exy(:,:,i) = PIT;   
end


for i=1:nlon;
    for j=1:nlat;
        depthij = depth(i,j);
        ind = zlev<=depthij;
        Fz(i,j,ind) =   exp((zlev(ind)-depthij)/ksi)...
                    /(1-exp(-depthij/ksi))/ksi;
    end
end

visc_diss = (q/rho)*Exy.*Fz;      % turbulent dissipation
Reb = visc_diss./N2/nu;           % buoyancy Reynolds number
Reb(Reb==inf)=NaN;
Reb(Reb<0)=NaN;

Ri = 0.4;              % Ri^start
[eff_DNS_par, Reb_DNS_par,Ri_DNS_par] = mixeff_param(Reb(:),Ri);
Gamma = reshape(eff_DNS_par./(1-eff_DNS_par),[nlon,nlat,nzlev]);

param = @(x,x0,y0)  2*y0*(x/x0).^(0.5)./(1+x/x0);
Gamma_param_lb = param(Reb,100,0.2);
Gamma_param_ub = param(Reb,300,0.5);

kappa_osb = nu*flux_coeff*Reb;     % Osborn    diapycnal diffusivity
kappa_par = nu*Gamma.*Reb;         % DNS-based diapycnal diffusivity


%% Find discrepency in tidal elevation by comparison to tpxo data

% t1 = reshape(nanmean((eta-eta_tpxo).^2,2),nlon,nlat);
% t2 = reshape(nanmean((eta_tpxo).^2,2),nlon,nlat);
% 
% % t1 = (amp-amp_tpxo).^2;
% % t2 = amp_tpxo.^2;
% % 
% % t1 = (phs-phs_tpxo).^2;
% % t2 = phs_tpxo.^2;
% 
% ocn_area1 = glbsum(t1*0+1,lat,lon,ind_dw & ind_lat);
% ocn_area2 = glbsum(t2*0+1,lat,lon,ind_dw & ind_lat);
% misfit_eta(ic) = sqrt(glbsum(t1,lat,lon,ind_dw & ind_lat)./ocn_area1); 
% misfit_tpxo(ic) = sqrt(glbsum(t2,lat,lon,ind_dw & ind_lat)./ocn_area2);
% misfit_variance =1-(misfit_eta./misfit_tpxo).^2;
% 
