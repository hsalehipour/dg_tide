function A = loadtxt(fname,ncol,nskip)
% Purpose: Read any columnar tabluated .dat file
% INPUTS: 
%       fname: file name
%       ncol : # of data columns
%       nskip: # of rows to skip at the beginning


fid = fopen(fname);
for i=1:nskip
    fgetl(fid)
end    
A = fscanf(fid, '%g', [ncol inf]);
fclose(fid);
end

