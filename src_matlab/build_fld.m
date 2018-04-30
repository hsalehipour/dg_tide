function fld = build_fld(amp,phs,omega,time)
% Purpose: Build a time and space dependent field fld(x,t) using the given
% amplitude, phase, frequency and time
% Note: AMP and PHS have been calculated based on
% fld=amp.*cos(omega*time-phs);

i = sqrt(-1);
fld_r = amp.*cos(phs*pi/180);     % real part
fld_c = amp.*sin(phs*pi/180);     % complex part

fld=real((fld_r + i*fld_c)*exp(-i*omega*time));

end
