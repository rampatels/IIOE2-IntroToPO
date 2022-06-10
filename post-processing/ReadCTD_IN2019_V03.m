% This script is to read CTD data in a proper structure for further
% analysis

% clear;clc; close all
% 
% CTD_dir = './in2019_v03CtdAvg_nc/';
% files_nc=dir([CTD_dir,'*.nc']);
% nfiles=length(files_nc);
% 
% for j=1:nfiles
%     [var, ~, ~]=nc2mat([CTD_dir,files_nc(j).name]);
%     CTD(j) = var;
%     lon(j) = var.longitude;
%     lat(j) = var.latitude;
%     CTD(j).time = var.time/60/24+datenum(1900,1,1,0,0,0); % convert to datenum format
%     max_press(j) = var.pressure(end);
%     clear var
% end
% clear CTD_dir files_nc nfiles j
% %   struct with fields:
% % 
% %                    time: 6.2781e+07
% %               woce_date: 20190515
% %               woce_time: 1.5317e+04
% %                pressure: [494×1 double]
% %                latitude: -34.7073
% %               longitude: 113.4239
% %             temperature: [1×494 double]
% %         temperatureFlag: [1×494 double]
% %            conductivity: [1×494 double]
% %        conductivityFlag: [1×494 double]
% %                salinity: [1×494 double]
% %            salinityFlag: [1×494 double]
% %           temperature_2: [1×494 double]
% %       temperature_2Flag: [1×494 double]
% %          conductivity_2: [1×494 double]
% %      conductivity_2Flag: [1×494 double]
% %              salinity_2: [1×494 double]
% %          salinity_2Flag: [1×494 double]
% %                  oxygen: [1×494 double]
% %              oxygenFlag: [1×494 double]
% %                oxygen_2: [1×494 double]
% %            oxygen_2Flag: [1×494 double]
% %             fluorometer: [1×494 double]
% %         fluorometerFlag: [1×494 double]
% %                     par: [1×494 double]
% %                 parFlag: [1×494 double]
% %               altimeter: [1×494 double]
% %           altimeterFlag: [1×494 double]
% %         transmissometer: [1×494 double]
% %     transmissometerFlag: [1×494 double]

% %% Checking Station positions and pressure
% 
% figure(1)
% clf
% scatter(lon, lat, 140, max_press)
% colorbar
% caxis([400, 5500])
% % So this figures shows that there are some stations away from 110E lines
% % which for now I would like to omit from the analysis. 
% % Therefore, it is important to account for the longitudes of the stations
% % and I am interested in looking at full column stations only. This leads
% % to the subsequent steps where, I seperate stations

% %% Separate deep and shallow stations
% % The station is considered in the vicinity of 110E if the actual
% % measurement of a longitude deviate less than 1 degree longitudes. This is
% % to accounting for the drift of ships during the package deployment. 
% 
% % This would gives us stations numbers that I care to seperate both shallow
% % stations and deep stations
% ind_deep = find(abs(lon-110)<1 & max_press>2000);
% ind_shal = find(abs(lon-110)<1 & max_press<600);
% 
% % Checking the casts
% figure(3);clf
% plot(lat(ind_deep), ind_deep, 'sb')
% hold on
% text(lat(ind_deep), ind_deep, num2str(ind_deep'), 'color', 'b')
% plot(lat(ind_shal), ind_shal, '.k')
% text(lat(ind_shal), ind_shal, num2str(ind_shal'), 'color', 'k')
% hold off
% 
% % The below notes talks about number of profiles not the deployments
% %----- Here 19 corresponds to 49th deployment and 20 corresponds to 47th--
% 
% %% NB - Deep casts: Notes from Helen
% % Stations 19 and 20 were switched, 19 is at 11.5S and 20 is at 12.5S so
% % need to sort in latitude before plotting
% [~, ind] = sort(lat(ind_deep));
% ind_deep = ind_deep(ind);
% 
% %NB - Shallow casts
% % Stations 19-21 are all at the same station at 12.5S because of CTD
% % problem trouble shooting.
% % Plot data to see which one to use. Last one had missing salinity above
% % 30m. Last one is the one that water samples came from. Should use this
% % one but use #20 since S is good to the surface
% ind_shal = [ind_shal(1:18) ind_shal(20) ind_shal(22)];
% 
% % checking stations order for the data recording
% figure(4);clf
% plot(ind_deep, 'sb')
% text(1:20, ind_deep, num2str(ind_deep'), 'color', 'r')
% hold on
% plot(ind_shal, '.k')
% text(1:20, ind_shal, num2str(ind_shal'), 'color', 'k')
% hold off
% %%
% station_110E=[1:20]; 
% n_stations = length(station_110E);
% station_lat=[-39.5:1.5:-12.5 -11.5];
% 
% 
% % figure(3)
% % clf
% % for j=1:n_stations
% %     subplot(1,3,1)
% %     plot(CTD(ind_shal(j)).temperature,CTD(ind_shal(j)).pressure,'-','linew',1.5)
% %     axis ij
% %     xlabel('T (C)')
% %     set(gca,'ylim',[0 500])
% %     subplot(1,3,2)
% %     plot(CTD(ind_shal(j)).salinity,CTD(ind_shal(j)).pressure,'-','linew',1.5)
% %     axis ij
% %     xlabel('S (psu)')
% %     set(gca,'ylim',[0 500])
% %     title(['Station ',int2str(station_110E(j)),' Lat ',num2str(CTD(ind_shal(j)).latitude)])
% %     subplot(1,3,3)
% %     plot(CTD(ind_shal(j)).oxygen,CTD(ind_shal(j)).pressure,'-','linew',1.5)
% %     axis ij
% %     xlabel('O2 (uM)')
% %     set(gca,'ylim',[0 500])
% %     pause
% % end
% 
% 
% figure(4)
% clf
% subplot(1,2,1) % deep
% plot(lon,lat,'k.','markersize',8)
% hold on 
% plot(lon(ind_deep),lat(ind_deep),'ro','markersize',8)
% text(lon(ind_deep)+0.05,lat(ind_deep), num2str(ind_deep'), 'color', 'b')
% title('110^{o}E Deep casts')
% hold off
% 
% subplot(1,2,2) % deep
% plot(lon,lat,'k.','markersize',8)
% hold on 
% plot(lon(ind_shal),lat(ind_shal),'ro','markersize',8)
% text(lon(ind_shal)+0.05,lat(ind_shal), num2str(ind_shal'), 'color', 'k')
% hold off
% title('110^{o}E Shallow casts')

% print(gcf, '-dpng', '-r100', '-painters', 'IN2019_V03Sampling')
% %% Seperating data first and saving them temporarily
% 
% % deep casts
% ctd19v03deep = CTD(ind_deep);
% for stn = 1:length(ind_deep)
%     ctd19v03deep(stn).deployment = ind_deep(stn);
% end
% 
% % shallow casts
% ctd19v03shal = CTD(ind_shal);
% for stn = 1:length(ind_shal)
%     ctd19v03shal(stn).deployment = ind_shal(stn);
% end
% 
% % Reorder the fields
% % I first identify how many fields are there in a structure
% % Then by looking the field name find that I need to move last added field
% % to first 
% 
% aa = fieldnames(ctd19v03deep);
% stlen = length(aa);
% ctd19v03deep = orderfields(ctd19v03deep, [stlen,1:stlen-1]);
% 
% % likewise for shallow casts
% aa = fieldnames(ctd19v03shal);
% stlen = length(aa);
% ctd19v03shal = orderfields(ctd19v03shal, [stlen,1:stlen-1]); clear stlen aa

%
% save deepshalcast ctd19v03deep ctd19v03shal

%% This part is based on commented portion
clear;clc; close all
addpath ~/'OneDrive - University of Tasmania'/MatlabData/HPdata/IOVoyage/

load deepshalcast % This contains all the stations of interest

% The idea here is to project data on regular pressure grids intervalled at
% 1dbar from 0dbar to 5850 dbar to accommodate the deepest cast

% Creating new structure for data
%
upres = [0:5850]'; %#ok

ctd19v03D = struct('Deployment', {}, 'lat', {}, 'lon', {}, 'time', {}, 'upres', {}, 'p', {}, 't', {}, 'SP', {}, 'oxy', {}, 't2', {}, 'SP2', {},  'oxy2', {});

for stn = 1:length(ctd19v03deep)
    ctd19v03D(stn).Deployment = ctd19v03deep(stn).deployment; 
    ctd19v03D(stn).lat = ctd19v03deep(stn).latitude; 
    ctd19v03D(stn).lon = ctd19v03deep(stn).longitude;
    ctd19v03D(stn).time = ctd19v03deep(stn).time;
    ctd19v03D(stn).p = ctd19v03deep(stn).pressure;
    ctd19v03D(stn).upres = upres;
    %
    p = ctd19v03deep(stn).pressure; 
    ctd19v03D(stn).t = ra_unipresvar(p, [ctd19v03deep(stn).temperature]', upres); % degree C
    ctd19v03D(stn).t2 = ra_unipresvar(p, [ctd19v03deep(stn).temperature_2]', upres);% degree C
    ctd19v03D(stn).SP = ra_unipresvar(p, [ctd19v03deep(stn).salinity]', upres); % psu
    ctd19v03D(stn).SP2 = ra_unipresvar(p, [ctd19v03deep(stn).salinity_2]', upres); % psu
    ctd19v03D(stn).oxy = ra_unipresvar(p, [ctd19v03deep(stn).oxygen]', upres); % micromol l-1 uM
    ctd19v03D(stn).oxy2 = ra_unipresvar(p, [ctd19v03deep(stn).oxygen_2]', upres); % micromol l-1 uM

%     ctd19v03D(ii).par = ra_unipresvar(p,ctd19v03deep.par); % microE/m^2/sec
%     ctd19v03D(ii).flu = ra_unipresvar(p,ctd19v03deep.fluorometer); % mg/m^3
%     ctd19v03D(ii).tr = ra_unipresvar(p,ctd19v03deep.transmissometer); % in %
%     ctd19v03D(ii).scat = ra_unipresvar(p,ctd19v03deep.obs); % m-1/sr %
 
    clear p
end
clear stn

% save IN19V03D ctd19v03D % This is what you will use for your analysis.. 

%% For shallow casts
clear;clc; close all

load deepshalcast

clear ctd19v03deep

% Creating new structure for data
ctd19v03S = struct('Deployment', {}, 'lat', {}, 'lon', {}, 'time', {}, 'upres', {}, 'p', {}, 't', {}, 'SP', {}, 'oxy', {}, 't2', {}, 'SP2', {},  'oxy2', {});

% The idea here is to project data on regular pressure grids intervalled at
% 1dbar from 0dbar to 600 dbar to accommodate the deepest cast
upres = [1:600]'; %#ok


for stn = 1:length(ctd19v03shal)
    disp(['station', ' ', num2str(ctd19v03shal(stn).deployment)])
    ctd19v03S(stn).Deployment = ctd19v03shal(stn).deployment; 
    ctd19v03S(stn).lat = ctd19v03shal(stn).latitude; 
    ctd19v03S(stn).lon = ctd19v03shal(stn).longitude;
    ctd19v03S(stn).time = ctd19v03shal(stn).time;
    ctd19v03S(stn).p = ctd19v03shal(stn).pressure;
    ctd19v03S(stn).upres = upres;
    %
    p = ctd19v03shal(stn).pressure; 
    ctd19v03S(stn).t = ra_unipresvar(p, [ctd19v03shal(stn).temperature]', upres); % degree C
    ctd19v03S(stn).t2 = ra_unipresvar(p, [ctd19v03shal(stn).temperature_2]', upres);% degree C
    ctd19v03S(stn).SP = ra_unipresvar(p, [ctd19v03shal(stn).salinity]', upres); % psu
    ctd19v03S(stn).SP2 = ra_unipresvar(p, [ctd19v03shal(stn).salinity_2]', upres); % psu
    ctd19v03S(stn).oxy = ra_unipresvar(p, [ctd19v03shal(stn).oxygen]', upres); % micromol l-1 uM
    ctd19v03S(stn).oxy2 = ra_unipresvar(p, [ctd19v03shal(stn).oxygen_2]', upres); % micromol l-1 uM

%     ctd19v03S(ii).par = ra_unipresvar(p,ctd19v03shal.par); % microE/m^2/sec
%     ctd19v03S(ii).flu = ra_unipresvar(p,ctd19v03shal.fluorometer); % mg/m^3
%     ctd19v03S(ii).tr = ra_unipresvar(p,ctd19v03shal.transmissometer); % in %
%     ctd19v03S(ii).scat = ra_unipresvar(p,ctd19v03shal.obs); % m-1/sr %
    
    clear p
end
clear stn

% save IN19V03S ctd19v03S % This is what you will use for your analysis.. 
