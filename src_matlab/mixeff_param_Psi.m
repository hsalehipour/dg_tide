function Psi = mixeff_param_Psi(ri)
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

% original fit to obtain right flank directly;
% Psi      = 0.02*exp(13*ri)+1.3;

% a revised fit to obtain right flank indirectly using (eff_peak, Reb_peak)
% as Psi = 3/2*sqrt(Reb_peak)*eff_peak
Psi      = 0.04*exp(12*ri)+1.5;

end
