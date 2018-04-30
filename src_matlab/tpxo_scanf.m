function [tpxlat, tpxlon, tpxamp, tpxphs] = tpxo_scanf(filename,nlon,nlat)

% Purpose: Read the output data of Tpx 7.2 dataset by filtering out the 
% text within the file.

clc;
SearchedString = '************* Site is out of model grid OR land ***************';
fid = fopen(filename);
i = 1;
fgetl(fid);
fgetl(fid);
fgetl(fid);
tline = fgetl(fid);
tpxlat = [];   tpxlon = [];  tpxamp = [];  tpxphs = [];
while ischar(tline)
   data = str2num(tline);
   if isempty(data)
       indx = strfind(tline, SearchedString);
       data = sscanf(tline(1:indx-1), '%f %f');
       tpxlat = [tpxlat; data(1)];
       tpxlon = [tpxlon; data(2)];
       tpxamp = [tpxamp; NaN];
       tpxphs  = [tpxphs;  NaN];
   else
       tpxlat = [tpxlat; data(1)];
       tpxlon = [tpxlon; data(2)];
       tpxamp = [tpxamp; data(3)];
       tpxphs  = [tpxphs;  data(4)];
   end
   tline = fgetl(fid);
   i = i + 1;
end
fclose(fid);

tpxlat = reshape(tpxlat, nlon,nlat);
tpxlon = reshape(tpxlon, nlon,nlat);
tpxamp = reshape(tpxamp, nlon,nlat);
tpxphs = reshape(tpxphs, nlon,nlat);

end