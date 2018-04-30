%==========================================================================
%                 Reading/ Analysing the output files of the
%                    TiDg Global Barotropic Tidal model 
% Ref:              Salehipour et al 2013 (Ocean Modell) 
%
% Written by :              Hesam Salehipour
%                           March 25th, 2015
%==========================================================================
% This is a rather automated code to analyse the data;
% Multiple notes to which one should pay attention:
%   1) Make sure the data.in file contains the correct path to your desired
%      cases 
%   2) All the neccessary analysis should be performed inside AnalyseData.m
%   3) All the ploting commands should be invoked inside PlotData.m
%==========================================================================


clear all;  close all;  clc;
% close all;  clc;

% Read number of cases to be compared (default = 1)
fname = 'data.in';
fidm  = fopen(fname);     
ncase = fscanf(fidm, '%f');



for ic=1:ncase
    % Reads the address for the current case
    fadrs = fgetl(fidm);
    
    % 1) Startup
    StartUp;
    
    % 1) Reading data
    ReadData;
        
    % 2) Analysis
    AnalyseData;

    % 3) Plot data
    PlotData;

end
fclose(fidm); 




