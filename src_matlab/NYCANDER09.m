%% NYCANDER09.m
%  MCB, NRL, 2013-05-03
%  modified to run on UTIL in parallel mode
%  create NYCANDER tensor for topography
%  *using GEBCO averaged to 120/60 arc seconds and
%  filled GDEM for stratification
%  *use filter wave length a2 for semidiurnal frequencies
%  *limit radius to 5a when computing dga
%  *split the world in smaller subgrids to do the analyiss
%  make sure there is an overlap of x
%  updated for higher-mode filtering
%
%% grid layaout
%            y
%           |
%           |
% __________.___________.-corner point
%           |           |
%           |    h11    |
%     h01   |   (i,j)-main point
%           |           |
%           |           |
% _____dJ/dx_(0,0)______.________ x
%           |           |
%           |           |
%    h00    |    h10    |DY
%           |           |
%           |           |
% __________.____DX_____.
%
%
%            y
% __________.___________.  INDICES relative to main grid
%           |           |          from which corner points are computed
%           |           |
%  (j+1,i)  | (j+1,i+1) |
%           |           |
%           |           |
% ________(j,i)_________._ x
%           |           |
%           |           |
%   (j,i)   |  (j,i+1)  |
%           |           |
%           |           |
% __________.___________.
%
% depth limit in hycom is 7200 m


clear all

matlabpool close force local
myCluster = parcluster();
matlabpool(myCluster)

addpath	/u/home/mbui/matcodes/NRL
addpath	/u/home/mbui/matcodes/NRL/funcs

aaa = tic;

%runname = 'lowres';
%runname = 'gebco120_01';
%runname = 'gebco120_angelique1f';  %old N
%runname = 'gebco120_angelique1g'; %better,deeper N

%runname = 'gebco120_hawaii2';  %old N
%runname = 'gebco120_hawaii3'; %better,deeper N

%runname = 'gebco60_tst1';
%runname = 'gebco60_tst2';
%runname = 'gebco60_tst7'; %global experiment, every 80 cells
%runname = 'gebco60_tst8';  %hires experiment, every 1 cell, 2,1-7
%runname = 'gebco60_tst9';  %hires experiment, every 1 cell, 2,1-7, no parfor, 100,10
runname = 'gebco60_r1';  %hires experiment, every 1 cell, 2,1-7, no parfor, 100,10

%% set switches/constants ========================================
savflg = 1; % save stuff
areafl = 0; % select small area
Ndeep  = 1; % use deep N, gives lower values in deep water!

% range over which bessel filter works
numr = 5; %
modenum = 3; %filter over shorther wave lengths

% number of boxes
nbx = 80;
nby = 10;

% indexes of the OUTER boxes ixo iyo
% find max a to determine range of overlap
%maxa = max(a2r(:));
%maxa = 80e3;  %use 80km as a limit
%maxa = 40e3;  %use 40km as a limit
maxa = 60e3;  %use 60km as a limit
maxD = maxa*numr;

%dirsave = '/p/CWFS5/mbui/data/TENSOROUT/';
%dirsave = '/home/mbui/data/matout/hycom_tensor/';
dirsave = '/net/pargo/export/data/matout/hycom_tensor/';


%% load bottom topo gebco_08_120sec, made in subsample_GEBCO_01.m
%dirin  = '/p/CWFS5/mbui/data/TENSORIN/';
%dirin  = '/home/mbui/data/projects/GEBCO/'

%load([dirin 'gebco_08_120sec.mat']) % z lon lat
load([dirin 'gebco_08_60sec.mat']) % 60 secs
disp(['finished loading grid file ....' ]);

sx = 1; sy = 1; %subsample to make it faster hires
lonS  = lon(1:sx:end);
latS  = lat(1:sy:end)';
zbot  = z(1:sy:end,1:sx:end);

ddxx1 = lon(2)-lon(1);   %% degrees
ddxx2 = lonS(2)-lonS(1); %% degrees
clear lon lat z

% check
if abs(360-ddxx2*length(lonS)) > 1e-3
    disp(['potential error; ddxx2 * ' num2str(length(lonS)) ' not equal to 360'])
    disp(['difference is ' num2str(abs(360-ddxx2*length(lonS)))])
else
    disp(['dx fills 360 circle; difference is ' num2str(abs(360-ddxx2*length(lonS)))])
end

% limit depths
mindepth = -11999;
zbot(zbot<mindepth) = mindepth;

disp(['max depth is ' num2str(min(zbot(:)))])

[nny,nnx] = size(zbot);
plon  = repmat(lonS,[nny 1]);
plat  = repmat(latS,[1 nnx]);

%% fill in depth SOUTH of 84
zbot(plat<-84 & zbot<0)=0.1;

%% do a small area
if areafl
    %xran = [-110 -100]+360; yran = [5 15]; % angelique bathy range
    %xran = [-180-2 -155+2]+360; yran = [15-2 40+2];  % hawaii
    xran = [195 207]; yran = [15 25];  % hawaii2
    
    Ixs = find(plon(1,:)>xran(1) & plon(1,:)<xran(2));
    Iys = find(plat(:,1)>yran(1) & plat(:,1)<yran(2));
    plon = plon(Iys,Ixs);
    plat = plat(Iys,Ixs);
    zbot = zbot(Iys,Ixs);
    
    %    figure
    %    pcolor(plon,plat,zbot); shading flat
    
    [nny,nnx] = size(zbot)
end


%% Split world in boxes ==============================================
nix = floor(nnx/nbx);
niy = floor(nny/nby);

% indexes of the INNER boxes ix iy
BX=[];

isy=0;
for jj=1:nby
    isx=0;
    for ii=1:nbx
        BX(jj,ii).ix = [isx+1:isx+nix];
        isx = isx+nix;
        
        BX(jj,ii).iy = [isy+1:isy+niy];
    end
    
    % remainder x
    BX(jj,ii).ix = BX(jj,ii).ix(1):nnx;
    
    isy = isy+niy;
end

% remainder y
for ii=1:nbx
    BX(jj,ii).iy = BX(jj,ii).iy(1):nny;
end

k=0;
for jj=1:nby
    for ii=1:nbx
        % find maximum dx and dy
        DXmax = max( spheric_dist(plat(BX(jj,ii).iy([1 end]),BX(jj,ii).ix(2)),plat(BX(jj,ii).iy([1 end]),BX(jj,ii).ix(1)),...
            plon(BX(jj,ii).iy([1 end]),BX(jj,ii).ix(2)),plon(BX(jj,ii).iy([1 end]),BX(jj,ii).ix(1))) );
        
        DYmax =      spheric_dist(plat(BX(jj,ii).iy(1),      BX(jj,ii).ix(1)),plat(BX(jj,ii).iy(2),      BX(jj,ii).ix(1)),...
            plon(BX(jj,ii).iy(1),      BX(jj,ii).ix(1)),plon(BX(jj,ii).iy(2),      BX(jj,ii).ix(1))) ;
        
        % find indices on main grid that border the box within
        numx = round(maxD/DXmax);
        numy = round(maxD/DYmax);
        
        ix = BX(jj,ii).ix;
        iy = BX(jj,ii).iy;
        
        % works only if indices are x,left to right and y,south to north
        % along x
        if     ix(1) == 1 % first point
            BX(jj,ii).ixo = [nnx-numx+1:nnx ix ix(end)+1:ix(end)+1+numx-1]; % add western boundary
        elseif ix(end) == nnx % last point
            BX(jj,ii).ixo = [ix(1)-numx:ix(1)-1 ix 1:numx];                    % add eastern boundary
        elseif ix(1)-numx<1
            BX(jj,ii).ixo = [nnx+(ix(1)-numx):nnx 1:ix(end)+1+numx-1];        % add western boundary
        elseif ix(end)+numx>nnx
            BX(jj,ii).ixo = [ix(1)-numx:nnx 1:ix(end)+numx-nnx];                % add eastern boundary
        else
            BX(jj,ii).ixo = [ix(1)-numx:ix(1)-1 ix ix(end)+1:ix(end)+1+numx-1];
        end
        
        % along y
        if     iy(1) == 1 % first point
            BX(jj,ii).iyo = [                   iy iy(end)+1:iy(end)+1+numy-1]; % at southern boundary
        elseif iy(end) == nny % last point
            BX(jj,ii).iyo = [iy(1)-numy:iy(1)-1 iy];                            % at northern boundary
        else
            BX(jj,ii).iyo = [iy(1)-numy:iy(1)-1 iy iy(end)+1:iy(end)+1+numy-1]; % add eastern boundary
        end
        
        %         k=k+1;
        %         if k>4; k=1; end
        %         plot(plon(BX(jj,ii).iyo,BX(jj,ii).ixo),plat(BX(jj,ii).iyo,BX(jj,ii).ixo),char(clr(k)))
        %         drawnow
    end
end


% figure
% hist(a2r(:)/1e3,100)

%% load stratification N^2:NN and depth-integrated stratifcation:NI
%dirout = '/p/CWFS5/mbui/data/TENSORIN/';
%dirout = '/home/mbui/data/matout/GDEM4/'
dirout = '/net/pargo/export/data/matout/GDEM4/'


if Ndeep
    load([dirout 'GDEM4_NN_NIdeep.mat']); % longv4 latgv4 depgv4  NN NI
else
    load([dirout 'GDEM4_NN_NI.mat']); % longv4 latgv4 depgv4  NN NI
end
disp(['finished loading stratification file ....' ]);

disp(['max depth from stratifcation file is ' num2str(max(depgv4(:)))])

% to be used for interpolating Nbot at zmid levels
Nb = sqrt(NN); clear NN

% add to the eastern boundaries for gridding purposes
longv4(end+1) = longv4(1)+360;
Nb(:,end+1,:) = Nb(:,1,:);
NI(:,end+1,:) = NI(:,1,:);

% coefficient to be multiplied with NI
lata  = repmat(latgv4,[1 length(longv4)]);

bet   = 1.455;
om2   = 12.1408331767/3600/24;
fcor  = coriolis(lata);
%acoef = bet./(pi*sqrt(om2^2-fcor.^2));
acoef = bet./(modenum*pi*sqrt(om2^2-fcor.^2));
acoef(abs(imag(acoef))>0)=0; %set to zero beyond turning latitude
NI    = NI.*repmat(acoef,[1 1 size(NI,3)]);


%% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%  LOOP over all boxes
%  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

%% input from screen
inputvaly1 = input(['type JJJ start: (max=' num2str(nby) '): '])
inputvaly2 = input(['type JJJ end    (max=' num2str(nby) '): '])
inputvalx1 = input(['type III start  (max=' num2str(nbx) '): '])
inputvalx2 = input(['type III end    (max=' num2str(nbx) '): '])

if inputvalx2>=inputvalx1
    xsteps = [inputvalx1:1:inputvalx2]
else
    xsteps = [inputvalx1:-1:inputvalx2]
end

if inputvaly2>=inputvaly1
    ysteps = [inputvaly1:1:inputvaly2]
else
    ysteps = [inputvaly1:-1:inputvaly2]
end


for JJJ=ysteps
    %    disp(['row ' num2str(JJJ) '  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'])
    
    for III=xsteps
        disp(['row    ' num2str(JJJ) '  ==================================================='])
        disp(['column ' num2str(III) '  ==================================================='])
        
        TT = tic;
        
        %% select outerbox area ==============================================
        if areafl == 0
            %  this box overlaps in x and y
            Iselx = BX(JJJ,III).ixo;
            Isely = BX(JJJ,III).iyo;
        else
            % angelique, outer bounds
            Iselx = find(plon(1,:)>xran(1) & plon(1,:)<xran(2));
            Isely = find(plat(:,1)>yran(1) & plat(:,1)<yran(2));
            
            % fudge it
            ddxx = BX(JJJ,III).ix(1)-BX(JJJ,III).ixo(1);
            ddyy = BX(JJJ,III).iy(1)-BX(JJJ,III).iyo(1);
            BX(JJJ,III).ix  = Iselx(1)+ddxx:Iselx(end)-ddxx;
            BX(JJJ,III).iy  = Isely(1)+ddyy:Isely(end)-ddyy;
            BX(JJJ,III).ixo = Iselx;
            BX(JJJ,III).iyo = Isely;
        end
        
        lon2 = plon(Isely,Iselx); % main coordinates of outer grid
        lat2 = plat(Isely,Iselx);
        zbo2 = zbot(Isely,Iselx);
        
        [ny,nx]=size(zbo2);
        %% correct for overlapping grid points at begin/end of grid
        %  make continuous lattitude
        if III==1
            Idiff = find(abs(diff(lon2(1,:)))>10); %difference larger than 10 degree
            lon2(:,1:Idiff) = lon2(:,1:Idiff)-360;
        elseif III==nbx
            Idiff = find(abs(diff(lon2(1,:)))>10); %difference larger than 10 degree
            lon2(:,Idiff+1:end) = lon2(:,Idiff+1:end)+360;
        end
        
        %% compute h-grid corner points and DX and DY ==============================================
        
        % mean depth at corner points
        zboc = zbo2(1:end-1,1:end-1)/4 + zbo2(1:end-1,2:end,:)/4 + ...
            zbo2(2:end,1:end-1)/4   + zbo2(2:end,2:end,:)/4;
        
        % mean positions at corner points
        lonc = lon2(1:end-1,1:end-1)/4 + lon2(1:end-1,2:end,:)/4 + ...
            lon2(2:end,1:end-1)/4   + lon2(2:end,2:end,:)/4;
        
        latc = lat2(1:end-1,1:end-1)/4 + lat2(1:end-1,2:end,:)/4 + ...
            lat2(2:end,1:end-1)/4   + lat2(2:end,2:end,:)/4;
        
        [nyc,nxc]=size(zboc);
        
        %         ff2=figure
        %         contour(lonc,latc,zboc,[0 0],'k-');
        %         holder
        
        % get dx and dy and map to main topographic (NOT corner) points
        % works for equidistant lat lon only
        dx1 = spheric_dist(lat2(:,1),lat2(:,1),lon2(:,1),lon2(:,2));
        dx  = repmat(dx1,[1 size(lat2,2)]);
        
        dy1 = spheric_dist(lat2(1,1),lat2(2,1),lon2(1,1),lon2(1,1));
        dy  = repmat(dy1,size(lat2));
        
        
        %% compute dhdx-y at corner points
        dx2  = dx(1:end-1,1:end-1)/2 + dx(2:end,1:end-1)/2;  %at
        dhdx = ( zbo2(1:end-1,2:end)   + zbo2(2:end,2:end) ...
            - zbo2(1:end-1,1:end-1) - zbo2(2:end,1:end-1) )./(2*dx2);
        
        dhdy = ( zbo2(2:end,  1:end-1) + zbo2(2:end,  2:end) ...
            - zbo2(1:end-1,1:end-1) - zbo2(1:end-1,2:end) )./(2*dy(1,1));
        
        %         figure
        %         pcolor((dhdx-dhdx2)./dhdx2*100); shading flat
        %         caxis([-1 1])
        %         colorbar
        
        %% ======================================================================
        % compute dJ/dx at corner points, BIG loop over zboc
        % ======================================================================
        %        flg=1; % for plotting below
        
        % prepare indexes referenced on main grid for inner grids
        % first cell inner grid may start later than outer grid
        Ipx = find( BX(JJJ,III).ix(1)==BX(JJJ,III).ixo );
        Ipy = find( BX(JJJ,III).iy(1)==BX(JJJ,III).iyo );
        
        INx = [1:length(BX(JJJ,III).ix)]+Ipx-1;
        INy = [1:length(BX(JJJ,III).iy)]+Ipy-1;
        
        % omit top corner point
        if JJJ == nby; INy(end) = []; end
        
        % declare matrices at corner points
        %        dJdxc = zeros(length(INy),length(INx));
        dJdxc = zeros(1,length(INy)*length(INx));
        dJdyc = dJdxc;
        dJdx  = dJdxc; dJdy = dJdxc;
        Nboci = dJdxc; Aci  = dJdxc;
        
        %if III==1;  ffff = figure; end
        %         k=k+1;
        %         if k>4; k=1; end
        %         clr = {'ro' 'b+' 'y*' 'g<'}
        %         figure(ff2)
        %         plot(lonc(INy,INx),latc(INy,INx),char(clr(k)))
        %         holder
        
        %        pause
        % loop over all, get indices from square matrix
        
        lenmi = repmat(1:length(INx),[length(INy) 1])';  %flip to have x in columns
        lenmj = repmat([1:length(INy)]',[1 length(INx)])';
        
        lenmi = reshape(lenmi,1,length(INy)*length(INx)); %along x, stacked
        lenmj = reshape(lenmj,1,length(INy)*length(INx));
        
        lenmx = repmat(INx,[length(INy) 1])';
        lenmy = repmat(INy',[1 length(INx)])';
        
        lenmx = reshape(lenmx,1,length(INy)*length(INx));
        lenmy = reshape(lenmy,1,length(INy)*length(INx));
        
        lenlen = length(INy)*length(INx);
        whose
        parfor lll=1:lenlen
            %         for lll=1:lenlen
            %         for lll=floor(lenlen*3/4):lenlen
            ll = lenmj(lll); %y index inner grid
            kk = lenmi(lll); %x index inner grid
            jj = lenmy(lll); %y index main grid
            ii = lenmx(lll); %x index main grid
            
            maxdep = -30;
            
            if zboc(jj,ii)<maxdep & sqrt(om2^2-coriolis(latc(jj,ii)).^2)>0
                
                if rem(lll,1000)==0; disp(['% ' num2str(lll/lenlen*100)]); end
                
                %%height at origin
                horg = zboc(jj,ii);
                
                %% 3D interpolate NI and NB on corner points of zboc
                
                % find nearest lona and lonb points
                % alway use subsequent points "<=" sign
                Ie = find(lonc(jj,ii)-longv4<=0); Ie = Ie(1);
                Iw = find(lonc(jj,ii)-longv4>0);  Iw = Iw(end);
                
                In = find(latc(jj,ii)-latgv4<=0); In = In(1);
                Is = find(latc(jj,ii)-latgv4>0);  Is = Is(end);
                
                zg   = -depgv4;
                zmid = zg(1:end-1)/2+zg(2:end)/2;
                Iug = find(horg-zg<0); Iug = Iug(end);
                Idg = find(horg-zg>0); Idg = Idg(1);
                
                Ium = find(horg-zmid<0); Ium = Ium(end);
                Idm = find(horg-zmid>0); Idm = Idm(1);
                
                % make matrices
                lona = repmat([longv4([Iw Ie]) ; longv4([Iw Ie])],[1 1 2]);
                lata = repmat([latgv4([Is In])   latgv4([Is In])],[1 1 2]);
                
                zzg=zeros(2,2,2); zzm=zzg;
                zzg(:,:,1) = zg([Iug Iug; Iug Iug]);
                zzg(:,:,2) = zg([Idg Idg; Idg Idg]);
                
                zzm(:,:,1) = zmid([Ium Ium; Ium Ium]);
                zzm(:,:,2) = zmid([Idm Idm; Idm Idm]);
                
                % interpolate Nb(zmid)
                Nboci(lll) = griddata(lona,lata,zzm,Nb([Is In],[Iw Ie],[Ium Idm]),...
                    lonc(jj,ii),latc(jj,ii),horg);
                % interpolate NI(zg)
                Aci(lll)   = griddata(lona,lata,zzg,NI([Is In],[Iw Ie],[Iug Idg]),...
                    lonc(jj,ii),latc(jj,ii),horg);
                
                if isnan(Nboci(lll)) | isnan(Aci(lll)); disp(['warning ' num2str(lll)]); end
                
                %% distance from origin to i,j points (r-r_ij)
                DS = spheric_dist(latc(jj,ii),lat2,lonc(jj,ii),lon2);
                %                 figure
                %                 pcolor(DS); shading flat
                
                % limit matrix to numr*Aci (Green and Nycander use 5*Aci)
                % differences are small
                Iyd = find(DS(:,ii)<numr*Aci(lll));
                Ixd = find(DS(jj,:)<numr*Aci(lll));
                
                %% obtain circle
                DS = DS(Iyd,Ixd);
                DS(DS>=numr*Aci(lll))=0;
                
                
                %             figure
                %             pcolor(lon2(Iyd,Ixd),lat2(Iyd,Ixd),DS);
                %             hold
                
                jjj = find(Iyd==jj);  % disp(['jjj = ' num2str([jjj jj])])
                iii = find(Ixd==ii);  %  disp(['iii = ' num2str([iii ii])])
                
                %% make sure there is data within 5*A (not in shallow water)
                if ~isempty(Iyd) & ~isempty(Ixd)
                    %% Green and Nycander (2012) page 7
                    % compute once as a matrix over main points
                    % local A is zero at topography
                    % BES also goes to inf for large values of x (large distance and shallow depth (is small A))
                    % derivative of bessel function ga'(r), comments redbook page 59
                    
                    nu   = 0;
                    a    = Aci(lll);
                    rr   = (DS/a).^2/8;
                    bes  = besseli(nu,rr);
                    expf = exp(-rr);
                    dr   = 2*DS/(a^2*8);
                    dexp = expf.*-dr;
                    %dbes = ( nu./rr.*besseli(nu,rr) + besseli(nu+1,rr) ).*dr;
                    dbes = ( 0 + besseli(nu+1,rr) ).*dr; %is faster
                    
                    dga = 1/a * ( ...
                        -a./DS.^2 ...
                        -sqrt(pi)/2*dexp.*bes ...
                        -sqrt(pi)/2*expf.*dbes );
                    
                    %% ============================================================
                    %  now along x
                    %  ============================================================
                    
                    %% find horizontal distance between xorg-xij
                    
                    % compute dx at level of lon2(JJ,:)
                    % mcb, this is much faster
                    dxp1      = -spheric_dist(lat2(Iyd,Ixd),lat2(Iyd,Ixd),lonc(jj,ii),lon2(Iyd,Ixd));
                    iw        = find(lon2(1,Ixd)<lonc(jj,ii)); %points to the west of lonc(jj,ii)
                    dxp       =  dxp1;
                    dxp(:,iw) = -dxp1(:,iw);
                    
                    
                    %% combine the stuff
                    % has infinity values at some topography, because A=ZERO
                    % Duh!
                    dJdxi = (zbo2(Iyd,Ixd)-horg).*dga.*dxp./DS.*dx(Iyd,Ixd).*dy(Iyd,Ixd); %at main points
                    dJdxi(DS==0)=0;
                    
                    %% add correction; h is on the main grid and the first corner grid is to
                    % the northeast of the the main grid origin
                    % h10 + h11 -h00 -h01 = h(jj,ii+1) + h(jj+1,ii+1) - h(jj,ii) - h(jj+1,ii)
                    hsumx = zbo2(jj,ii+1) + zbo2(jj+1,ii+1) - zbo2(jj,ii) - zbo2(jj+1,ii);
                    dxa  = dx(jj,ii+1)/4 + dx(jj+1,ii+1)/4 + dx(jj,ii)/4 + dx(jj+1,ii)/4;
                    dya  = dy(jj,ii+1)/4 + dy(jj+1,ii+1)/4 + dy(jj,ii)/4 + dy(jj+1,ii)/4;
                    s = dxa/dya;
                    
                    corrx = ( 2/s*log(sqrt(s^2+1)+s)-4*s^2/(1+s^2)^(3/2) )*hsumx;
                    corX  = ( 2/s*log(sqrt(s^2+1)+s) )*hsumx;
                    
                    %% add correction
                    % what is effect of the sum
                    % sum over a limited area
                    %                 dJdx(ll,kk)       = corrx + sum(dJdxi(~isnan(dJdxi)));
                    %                 dJdx_nocor(ll,kk) =         sum(dJdxi(~isnan(dJdxi)));
                    
                    % exlude the error correction, just remove these cells around corner point jj,ii
                    % include check if not enough points - like in shallow water
                    if length(Iyd)==1 & length(Ixd)>1
                        dJdxc(lll)      = corX  + sum(dJdxi(:)) - ...
                            (dJdxi(jjj,iii) + dJdxi(jjj,iii+1));
                    elseif length(Iyd)>1 & length(Ixd)==1
                        dJdxc(lll)      = corX  + sum(dJdxi(:)) - ...
                            (dJdxi(jjj,iii) + dJdxi(jjj+1,iii));
                    elseif length(Iyd)==1 & length(Ixd)==1
                        dJdxc(lll)      = corX  + sum(dJdxi(:)) - ...
                            (dJdxi(jjj,iii) + dJdxi(jjj,iii));
                    else
                        dJdxc(lll)      = corX  + sum(dJdxi(:)) - ...
                            (dJdxi(jjj,iii) + dJdxi(jjj,iii+1) + dJdxi(jjj+1,iii) + dJdxi(jjj+1,iii+1));
                    end
                    
                    % uncorrected
                    dJdx(lll)      = sum(dJdxi(:));
                    
                    %% ============================================================
                    %  now along y
                    %  ============================================================
                    
                    %% find horizontal distance between yorg-yij
                    % ofcourse the distances are the same along y axis as for method 1 !!!
                    dyp1  = -spheric_dist(latc(jj,ii),lat2(Iyd,Ixd),lon2(Iyd,Ixd),lon2(Iyd,Ixd));
                    
                    iw        =  find(lat2(Iyd,ii)<latc(jj,ii)); %points to the south of latc(jj,ii)
                    dyp       =  dyp1;
                    dyp(iw,:) = -dyp1(iw,:);
                    
                    
                    %% combine the stuff, all points, including the bad ones
                    dJdyi = (zbo2(Iyd,Ixd)-horg).*dga.*dyp./DS.*dx(Iyd,Ixd).*dy(Iyd,Ixd); %at main points
                    dJdyi(DS==0)=0;
                    
                    %% add correction; h is on the main grid
                    % jj,ii on corner grid is to the north-east of main grid
                    hsumy = zbo2(jj+1,ii+1) + zbo2(jj+1,ii) - zbo2(jj,ii+1) - zbo2(jj,ii);
                    
                    t = dya/dxa;
                    corry = ( 2/t*log(sqrt(t^2+1)+t)-4*t^2/(1+t^2)^(3/2) )*hsumy;
                    corY  = ( 2/t*log(sqrt(t^2+1)+t) )*hsumy;
                    
                    %% add correction
                    % exlude the error correction, just remove these cells around corner point jj,ii
                    % include check if not enough points - like in shallow water
                    if length(Iyd)==1 & length(Ixd)>1
                        dJdyc(lll)      = corY  + sum(dJdyi(:)) - ...
                            (dJdyi(jjj,iii) + dJdyi(jjj,iii+1));
                    elseif length(Iyd)>1 & length(Ixd)==1
                        dJdyc(lll)      = corY  + sum(dJdyi(:)) - ...
                            (dJdyi(jjj,iii) + dJdyi(jjj+1,iii));
                    elseif length(Iyd)==1 & length(Ixd)==1
                        dJdyc(lll)      = corY  + sum(dJdyi(:)) - ...
                            (dJdyi(jjj,iii) + dJdyi(jjj,iii));
                    else
                        dJdyc(lll)      = corY + sum(dJdyi(:)) - ...
                            (dJdyi(jjj,iii) + dJdyi(jjj,iii+1) + dJdyi(jjj+1,iii) + dJdyi(jjj+1,iii+1));
                    end
                    
                    % uncorrected
                    dJdy(lll)      = sum(dJdyi(:));
                    
                    
                    %         disp([ num2str(lll) '; ' num2str(dJdy(lll))])
                    %                else
                    %                    flgsw=1;
                    %                    isw = isw+1;
                    %                end
                end %if enough datapoints
            end  %if south of critical latitude
            
            %            if flgsw==1; disp([num2str(isw) ' shallow water points']); end
            
            %            avtime=toc;
            %            timeleft = avtime*(length(INy)-ll);
            %            disp(['time in this step: ' num2str(avtime) ' seconds']);
            %            disp(['time left:         ' num2str(timeleft/3600) ' hours']);
            
            %            disp([' ']);
            %                    disp(['end lll=' num2str(lll)]) ;
            
        end  % loop over x*x
        
        %reshape to square matrix
        dJdxc = reshape(dJdxc,length(INx),length(INy))';
        dJdyc = reshape(dJdyc,length(INx),length(INy))';
        dJdx  = reshape(dJdx,length(INx),length(INy))';
        dJdy  = reshape(dJdy,length(INx),length(INy))';
        Aci   = reshape(Aci,length(INx),length(INy))';
        Nboci = reshape(Nboci,length(INx),length(INy))';
        %        lenmi = reshape(lenmi,length(INx),length(INy))';
        %        lenmi(1,[1 end]),lenmi(2,[1 end]),whos INx INy
        %        lenmj = reshape(lenmj,length(INx),length(INy))';
        %        lenmj(1,[1 end]),lenmj(2,[1 end]),whos INx INy
        %        pause
        
        %% save inner grid only
        lonci  = lonc(INy,INx);
        latci  = latc(INy,INx);
        zboci  = zboc(INy,INx);
        dhdxi  = dhdx(INy,INx);
        dhdyi  = dhdy(INy,INx);
        
        if savflg == 1
            eval(['save ' dirsave 'tensor_' runname '_' num2str(JJJ) '_' num2str(III) '.mat JJJ III BX dJdxc dJdx dJdyc dJdy Aci Nboci lonci latci zboci dhdxi dhdyi'])
            disp(['tensor_' runname '_' num2str(JJJ) '_' num2str(III) '.mat saved ....'])
            disp([' '])
        end
        
        disp(['-------------------------------------------------------------------']);
        disp(['time in this tile: ' num2str(toc(TT)/3600) ' h and ' num2str(toc(TT)) ' s']);
        disp(['-------------------------------------------------------------------']);
        pause(0.25)
        
        %          figure(ff2)
        %          plot(lonci([1 end],:),latci([1 end],:),'ro')
        %          plot(lonci(:,[1 end]),latci(:,[1 end]),'ro')
        %          drawnow
        
        clear lonci latci zboci dhdxi dhdyi
        
    end
end

toc(aaa)




