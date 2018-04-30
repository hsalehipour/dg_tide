%% 0) Simulation startup

% Problem Setup
fname = [fadrs,'siminfo.dat'];
fid = fopen(fname);  
fgetl(fid);
bf_fname = fgetl(fid);    
fgetl(fid);
interp_fname = fgetl(fid);    
fgetl(fid);
omega    = fscanf(fid, '%f');
fclose(fid); 

% define constants
g   = 9.80616;
rho = 1035;
earthR = 6.37122e6;
nu = 1e-6;