% This script prepares chlorophyll climatology
% for may and jun only...
addpath ~/'OneDrive - University of Tasmania'/MatlabData/HPdata/IOVoyage/
clear;clc
lon = ncread('MODISmclimchla0319.nc', 'LON6241_7440');
lat = ncread('MODISmclimchla0319.nc', 'LAT961_2161');
times = ncread('MODISmclimchla0319.nc', 'TIME');
chla = ncread('MODISmclimchla0319.nc', 'CHLOR_A'); % mg/m3

%% Calculating May-June climatology

edates = datenum(0001,01,01, 00, 00, 00) + times;
t = datetime(datestr(edates));
edaymonth = month(t);
%
tInd = edaymonth >=5 & edaymonth <=6;

% checking dates
clc;
datestr(edates(tInd))

%% get corresponding chlorophyll-a data

chlamj = chla(:, :, tInd);
chlamj = squeeze(nanmean(chlamj, 3));

%% Checking data

figure(1);clf
pcolor(lon, lat, log10(chlamj)');
shading flat
cmocean('algae')

%% saving the data
lonmodis = lon;
latmodis = lat;

save chlamodisaquaWA lonmodis latmodis chlamj