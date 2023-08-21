% This script I read Hydrographic data of DM1963_V02 voyage along 110E line


clear; clc; close all
CTD_dir = './netcdf/DM1963V02/';
% [ctddmv02, ind_deepv02, ind_shalv02] = ra_1960DMread(CTD_dir);
% 
% CTD_dir = './netcdf/DM1963V01/';
% [ctddmv01, ind_deepv01, ind_shalv01] = ra_1960DMread(CTD_dir);

files_nc=dir([CTD_dir,'*.nc']);
nfiles=length(files_nc);

hy63v02 = struct();

for j=1:nfiles
    name = files_nc(j).name;
    [var, att, ~] = nc2mat([CTD_dir, name]);
    hy63v02(j).cruise = att.SURVEY_NAME; % string 
    hy63v02(j).station = att.STATION_NO; % string 
    hy63v02(j).date = att.START_TIME; % string in yyyy-mm-dd hh:mm:ss
    hy63v02(j).lat =  att.START_LAT;
    hy63v02(j).lon =  att.START_LON;
    hy63v02(j).bot_depth =  att.BOTTOM_DEPTH; 
    hy63v02(j).pres =  var.PRESSURE; %dbar
    hy63v02(j).max_depth = max(var.PRESSURE);
    hy63v02(j).temp =  var.TEMPERATURE; %C
    hy63v02(j).psal =  var.SALINITY; %psu
    hy63v02(j).oxy =  var.OXYGEN; %uM
    hy63v02(j).no3 =  var.NITRATE; %uM
    hy63v02(j).nit =  var.NITRITE;
    hy63v02(j).phos =  var.PHOSPHATE; %uM
    hy63v02(j).nh3 =  var.AMMONIA;
    % Flags: 'Good = 0; Below detection limit = 63; Suspect = 69; Bad = 133; Missing = 141; Unprocessed = 192'
    clear var name dind att
end
clear CTD_dir files_nc nfiles j

% Checking Station positions and pressure
lon = [hy63v02.lon]; lat = [hy63v02.lat]; max_press = [hy63v02.max_depth];
%
figure(1)
clf
scatter(lon, lat, 140, max_press)
colorbar
caxis([400, 5500])

%% Separate deep and shallow stations
% The station is considered in the vicinity of 110E if the actual
% measurement of a longitude deviate less than 1 degree longitudes. This is
% to accounting for the drift of ships during the package deployment. 
deployment = [hy63v02.station];

% This would gives us stations numbers that I care to seperate both shallow
% stations and deep stations
ind_deep = find(abs(lon-110)<1 & max_press > 2000);
ind_shal = find(abs(lon-110)<1 & max_press < 600);

%
figure(4); clf
subplot(1,2,1) % deep
plot(lon,lat,'k.','markersize',8)
hold on 
plot(lon(ind_deep),lat(ind_deep),'ro','markersize',8)
text(lon(ind_deep)+0.2,lat(ind_deep), num2str(deployment(ind_deep)'), 'color', 'b')
text(lon(ind_deep)-2, lat(ind_deep), num2str(max_press(ind_deep)'), 'color', 'k')
title('110^{o}E Deep casts DM1963V02')
hold off

subplot(1,2,2) % shallow
plot(lon,lat,'k.','markersize',8)
hold on 
plot(lon(ind_shal),lat(ind_shal),'ro','markersize',8)
text(lon(ind_shal)+0.2,lat(ind_shal), num2str(deployment(ind_shal)'), 'color', 'k')
text(lon(ind_shal)-2, lat(ind_shal), num2str(max_press(ind_shal)'), 'color', 'k')
hold off
title('110^{o}E Shallow casts DM1963V02')

% print(gcf, '-dpng', '-r100', '-painters', 'DM1963_V02SamplingHydro')
%% Comparing stations position with 2019 voyage
load deepshalhydrocast.mat

figure(5);clf
subplot(1,2,1)
plot([ctd19v03shal.longitude], [ctd19v03shal.latitude], '.k', 'MarkerSize', 8)
hold on
plot(lon(ind_shal),lat(ind_shal),'ko','markersize',8)
hold off
legend('2019 shal', '1963V02 shal')
title('Shallow Casts: 2019 and 1963V02')

subplot(1,2,2)
plot([ctd19v03deep.longitude], [ctd19v03deep.latitude], '.r', 'MarkerSize', 8)
hold on
% plot([ctd19v03deep.longitude], [ctd19v03deep.latitude], '.r', 'MarkerSize', 8)
% plot(lon(ind_shal),lat(ind_shal),'ko','markersize',8)
plot(lon(ind_deep),lat(ind_deep),'ro','markersize',8)
hold off
text(lon(ind_deep)-0.1,lat(ind_deep), num2str(deployment(ind_deep)'), 'color', 'r')
% legend('2019 shal', '2019 deep', '1963V02shal', '1963V02deep')
legend('2019 deep', '1963V02deep')
title('Selected stations 110E: 2019 & 1963V02')
% % print(gcf, '-dpng', '-r100', '-painters', 'DM1963V02nIN2019V03SelectHydro')

% figure(6);clf
% plot([ctd19v03shal.longitude], [ctd19v03shal.latitude], '.k', 'MarkerSize', 8)
% hold on
% plot([ctd19v03deep.longitude], [ctd19v03deep.latitude], '.r', 'MarkerSize', 8)
% plot(lon(ind_shal),lat(ind_shal),'ko','markersize',8)
% plot(lon(ind_deep),lat(ind_deep),'ro','markersize',8)
% for stn = 1:length(hy63v02)
%     ox = hy63v02(stn).oxy;
%     if all(isnan(ox))
%         text(lon(stn)-0.2, lat(stn), num2str(deployment(stn)), 'color', 'b')
%         plot(lon(stn)-0.22, lat(stn), 'xr', 'linewi', 1.5)
%     end
% end
% hold off
% text(lon(ind_shal)+0.1,lat(ind_shal), num2str(deployment(ind_shal)'), 'color', 'k', 'linewi', 2)
% text(lon(ind_deep)-0.1,lat(ind_deep), num2str(deployment(ind_deep)'), 'color', 'r')
% legend('2019V03 Shal', '2019V03 Deep', '1963V02 Shal', '1963V02 Deep')
% xlim([109.6, 110.6])
% title('Selecting 1963V02 Stations')
% set(gca, 'linewi', 2, 'fontsize', 12, 'fontweigh', 'bold')

% print(gcf, '-dpng', '-r200', '-painters', 'DM1963V02cIN2019V03SelectHydro')

clear ctd*
%% Rearranging latitudes from 40s to 10s
% I am going for only deep cast because, no nutrient information for
% shallow cast in this voyage - as per the .nc files recording.. 

[~, ind] = sort(lat(ind_deep));
ind_deep1 = ind_deep(ind);

%
hy63v02deep = hy63v02(ind_deep1);

%% Fixing nitrate data as pointed out by Peter Thompson
% Thompson words "I note there is an error in the data that gives that anomalous nitrate point at 27.5S and 73m depth in 1963.  
% It seems likely there was a transcription error with the 73m depth reported as 32.3µM NO3 and 1930m depth reported
% as 0.400 µM NO3 on this cast.  Clearly the latter is impossible.  
% I have reported it to CSIRO.  You can either delete the points or change the numbers.  
% I think I have done the analysis both ways and it makes little difference."

% So, I decided to change the numbers, rather than deleting them..
a = hy63v02deep(4).no3;
b = a; % Here 4 position = 32.3uM while 14 is 0.400uM. So I swithed them..
b(4) = a(14); b(14) = a(4);

hy63v02deep(4).no3 = b;

%% Rewriting data into proper formate at uniform depths
% Creating new structure for data

ctdHy63v02D500 = struct('Cruise', {}, 'Stations', {}, 'time', {}, 'lat', {}, 'lon', {}, 'sdepth', {}, 't', {}, 'SP', {}, ...
    'upres', {}, 'ct', {}, 'SA', {}, 'oxy', {}, 'no3', {}, 'no2', {}, 'phos', {}, 'ammonia', {});

% The idea here is to project data on regular pressure grids intervalled at
% 1dbar from 0dbar to 600 dbar to accommodate the deepest cast
% upres = [0:5:100, 110:10:200, 250:50:500];
upres = [0:5:500];
%----
for stn = 1:length(hy63v02deep)
    disp(['station', ' ', num2str(hy63v02deep(stn).station)])
    ctdHy63v02D500(stn).Cruise = hy63v02deep(stn).cruise;
    ctdHy63v02D500(stn).Stations = hy63v02deep(stn).station; 
    ctdHy63v02D500(stn).lat = hy63v02deep(stn).lat; 
    ctdHy63v02D500(stn).lon = hy63v02deep(stn).lon;
    ctdHy63v02D500(stn).time = hy63v02deep(stn).date;
    ctdHy63v02D500(stn).sdepth = double(hy63v02deep(stn).pres);
    ctdHy63v02D500(stn).upres = upres';
    p = double(hy63v02deep(stn).pres);
    % writing biogeochemical properties at uniform depths
    ctdHy63v02D500(stn).oxy = gsw_tracer_interp(hy63v02deep(stn).oxy, p, upres); % micromol l-1 uM
    ctdHy63v02D500(stn).no3 = gsw_tracer_interp(hy63v02deep(stn).no3, p, upres); % micromol l-1 uM
    ctdHy63v02D500(stn).no2 = gsw_tracer_interp(hy63v02deep(stn).nit, p, upres); % micromol l-1 uM
    ctdHy63v02D500(stn).phos = gsw_tracer_interp(hy63v02deep(stn).phos, p, upres); % micromol l-1 uM
    ctdHy63v02D500(stn).ammonia = gsw_tracer_interp(hy63v02deep(stn).nh3, p, upres); % micromol l-1 uM

    % reading original data for physical properties
    ctdHy63v02D500(stn).t = hy63v02deep(stn).temp; % degree C
    ctdHy63v02D500(stn).SP = hy63v02deep(stn).psal; % psu
    % Converting to CT and SA and interpolating at uniform depths
    sa = gsw_SA_from_SP(hy63v02deep(stn).psal, p, hy63v02deep(stn).lon, hy63v02deep(stn).lat); % Absolute Salinity g/kg
    ct = gsw_CT_from_t(sa, [hy63v02deep(stn).temp], p); % Conservative Temperature
    % projecting all on the uniform depth
    [SA, CT] = gsw_SA_CT_interp(sa, ct, p, upres);
    %
    ctdHy63v02D500(stn).ct = CT; % degree C
    ctdHy63v02D500(stn).SA = SA; % g/kg
   
    clear p sa ct SA CT 
end
clear stn

%% Checking original and interpolated values 

figure(5);clf
for stn = 1:15
    plot(hy63v02deep(stn).no3, ctdHy63v02D500(stn).sdepth, '*k')
    hold on
    plot(ctdHy63v02D500(stn).no3, ctdHy63v02D500(stn).upres, 'or'); hold off; axis ij
    ylim([0, 600])
    pause()
end

%% Saving final product for further analysis 
% save DM63V02Hy ctdHy63v02D500


%% Preparing deep cast for other analysis.. 
% By rewriting data into proper formate at uniform depths
% Creating new structure for data

ctdHy63v02D = struct('Cruise', {}, 'Stations', {}, 'time', {}, 'lat', {}, 'lon', {}, 'sdepth', {}, 't', {}, 'SP', {}, ...
    'upres', {}, 'ct', {}, 'SA', {}, 'dynhgt', {}, 'oxy', {}, 'no3', {}, 'no2', {}, 'phos', {}, 'ammonia', {});

% The idea here is to project data on regular pressure grids intervalled at
% 1dbar from 0dbar to 600 dbar to accommodate the deepest cast
upres = [0:5:100, 110:10:200, 250:50:500, 600:100:2000, 2500:500:4500]';

%----
for stn = 1:length(hy63v02deep)
    disp(['station', ' ', num2str(hy63v02deep(stn).station)])
    ctdHy63v02D(stn).Cruise = hy63v02deep(stn).cruise;
    ctdHy63v02D(stn).Stations = hy63v02deep(stn).station; 
    ctdHy63v02D(stn).lat = hy63v02deep(stn).lat; 
    ctdHy63v02D(stn).lon = hy63v02deep(stn).lon;
    ctdHy63v02D(stn).time = hy63v02deep(stn).date;
    ctdHy63v02D(stn).sdepth = double(hy63v02deep(stn).pres);
    ctdHy63v02D(stn).upres = upres;
    p = double(hy63v02deep(stn).pres);
    % writing biogeochemical properties at uniform depths
    ctdHy63v02D(stn).oxy = gsw_tracer_interp(hy63v02deep(stn).oxy, p, upres); % micromol l-1 uM
    ctdHy63v02D(stn).no3 = gsw_tracer_interp(hy63v02deep(stn).no3, p, upres); % micromol l-1 uM
    ctdHy63v02D(stn).no2 = gsw_tracer_interp(hy63v02deep(stn).nit, p, upres); % micromol l-1 uM
    ctdHy63v02D(stn).phos = gsw_tracer_interp(hy63v02deep(stn).phos, p, upres); % micromol l-1 uM
    ctdHy63v02D(stn).ammonia = gsw_tracer_interp(hy63v02deep(stn).nh3, p, upres); % micromol l-1 uM

    % reading original data for physical properties
    ctdHy63v02D(stn).t = hy63v02deep(stn).temp; % degree C
    ctdHy63v02D(stn).SP = hy63v02deep(stn).psal; % psu
    % Converting to CT and SA and interpolating at uniform depths
    sa = gsw_SA_from_SP(hy63v02deep(stn).psal, p, hy63v02deep(stn).lon, hy63v02deep(stn).lat); % Absolute Salinity g/kg
    ct = gsw_CT_from_t(sa, [hy63v02deep(stn).temp], p); % Conservative Temperature
    % projecting all on the uniform depth
    [SA, CT] = gsw_SA_CT_interp(sa, ct, p, upres);
    %
    ctdHy63v02D(stn).ct = CT; % degree C
    ctdHy63v02D(stn).SA = SA; % g/kg  
    % adding Dynamic Height
    dyheight = gsw_geo_strf_dyn_height(SA, CT, upres, 2000); % in $m^2/s^2$
    ctdHy63v02D(stn).dynhgt = dyheight;
   
    clear p sa ct SA CT dyheight
end
clear stn

%% Checking original and interpolated values 

figure(5);clf
for stn = 1:15
    plot(hy63v02deep(stn).no3, ctdHy63v02D(stn).sdepth, '*k')
    hold on
    plot(ctdHy63v02D(stn).no3, ctdHy63v02D(stn).upres, 'or'); hold off; axis ij
    ylim([0, 4500])
    pause()
end

%% Saving final product for further analysis 
% save DM63V02HyD ctdHy63v02D