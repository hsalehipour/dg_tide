function [eff_par, Reb_par,Ri_par] = mixeff_param(reb,ri)
% =========================================================================
% ========== multiparameter Parameterization of mixing efficiency==========
% =========================================================================

% Written by: Hesam Salehipour,     Aug 14, 2015
% Purpose:
% multiparameter parametrization of irreversible mixing
% efficiency based on the DNS analyses of KH-induced stratified turbulence.
% 
% inputs: 
% ri :   gradient Richardson number, Ri  = N^2/<S>^2
% reb:   buoyancy Reynolds number,   Reb = eps/(nu*N^2)
%
% Citation: 
% [1] Salehipour et al (2015), GRL, submitted.
% [2] Salehipour, H., and W. R. Peltier (2015), Diapycnal diffusivity,
% turbulent prandtl number and mixing efficiency in boussinesq stratified
% turbulence, J. Fluid Mech.
%
%

[Reb_par,Ri_par]=meshgrid(reb,ri);
eff_par = Reb_par*0;   

%===================
% Using Ri = <N^2>/<S>^2
%===================
% eff_peak = Ri_par/6.*(1-Ri_par)./(Ri_par.^3-1.95*Ri_par.^2+Ri_par+0.02);
% Psi    = 3.9*Ri_par+1.1;
% Reb_peak = (Psi./eff_peak).^2;
% 
% ind_lf = Reb_par<Reb_peak;
% ind_rf = Reb_par>=Reb_peak;
% 
% % eff_par(ind_lf) = -1/6*Psi(ind_lf).*Reb_par(ind_lf).^(-0.9)+eff_peak(ind_lf); 
% eff_par(ind_lf) = (-2*Reb_par(ind_lf).^(-0.9)+1).*eff_peak(ind_lf); 
% 
% eff_par(ind_rf) = Psi(ind_rf).*Reb_par(ind_rf).^(-0.5);
% eff_par(eff_par<0) = 1e-6;


% %===================
% % Using Ri=<N^2>/<S^2> : BUT discontinous
% %===================
% Psi = mixeff_param_Psi(Ri_par);
% eff_peak = mixeff_param_effpeak(Ri_par);
% Reb_peak = (Psi./eff_peak).^2;
% 
% ind_lf = Reb_par<Reb_peak;
% ind_rf = Reb_par>Reb_peak;
% 
% eff_par(ind_lf) = (-2*Reb_par(ind_lf).^(-0.8)+1).*eff_peak(ind_lf); 
% eff_par(ind_rf) = Psi(ind_rf).*Reb_par(ind_rf).^(-0.5);
% eff_par(eff_par<0) = 1e-6;

% %===================
% % Using Ri=<N^2>/<S^2>  : continuous
% %===================
% Psi = mixeff_param_Psi(Ri_par);
% eff_peak = mixeff_param_effpeak(Ri_par);
% Reb_peak = (Psi./eff_peak).^2;
% 
% ind_lf = Reb_par<=Reb_peak;
% ind_rf = Reb_par>Reb_peak;
% 
% eff_par(ind_lf) = (Reb_par(ind_lf)./Reb_peak(ind_lf)).^(0.1).*eff_peak(ind_lf); 
% eff_par(ind_rf) = Psi(ind_rf).*Reb_par(ind_rf).^(-0.5);
% eff_par(eff_par<0) = 1e-6;


%===================
% Using Ri=<N^2>/<S^2>  : tst attempt to get rid of spike in histograms
%===================
Psi = mixeff_param_Psi(Ri_par);
eff_peak = mixeff_param_effpeak(Ri_par);
Reb_peak = 4/9*(Psi./eff_peak).^2;

ind_lf = Reb_par<=Reb_peak;
ind_rf = Reb_par>Reb_peak;

param = @(x,p)  (1+2*p)*x.^p./(1+2*p*x.^(p+0.5));

eff_par(ind_lf) = eff_peak(ind_lf).*param((Reb_par(ind_lf)./Reb_peak(ind_lf)),0.55); 
eff_par(ind_rf) = eff_peak(ind_rf).*param((Reb_par(ind_rf)./Reb_peak(ind_rf)),1); 
eff_par(eff_par<0) = 1e-6;
end
