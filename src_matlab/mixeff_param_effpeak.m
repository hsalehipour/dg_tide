function eff_peak = mixeff_param_effpeak(ri)
% =========================================================================
% ========== multiparameter Parameterization of mixing efficiency==========
% =========================================================================

% Written by: Hesam Salehipour,     Aug 14, 2015
% Purpose:
% multiparameter parametrization of irreversible mixing
% efficiency based on the DNS analyses of KH-induced stratified turbulence.
% 
% inputs: 
% ri :   gradient Richardson number, Ri  = N^2/S^2
% reb:   buoyancy Reynolds number,   Reb = eps/(nu*N^2)
%
% Citation: 
% [1] Salehipour et al (2015), GRL, submitted.
% [2] Salehipour, H., and W. R. Peltier (2015), Diapycnal diffusivity,
% turbulent prandtl number and mixing efficiency in boussinesq stratified
% turbulence, J. Fluid Mech.
%
%

% original fit (commented to be "unwiledly")/ Let's make it simpler
% eff_peak = 0.04*ri.*(1-2*ri)./...
%     (ri.^3-1.05*ri.^2+0.26*ri+0.01);

% % FOR LATER: A refined, simpler eqn:
Ri_peak = 0.4;
x = ri./Ri_peak;
eff_peak = (1/3)*9*(x.^1)./(8+x.^9);
end
