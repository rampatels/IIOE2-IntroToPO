% This script I read Hydrographic data of IN2019_V03 voyage along 110E line

% The data that used in this script were obtained from Helen, currently at
% OneDrive>don't sync>IN2019_V03 Analysis>processed Hydrography

clear; clc; close all
CTD_dir = './in2019_v03_HydroDep_NetCDF/';
files_nc=dir([CTD_dir,'*.nc']);
nfiles=length(files_nc);

for j=1:nfiles
    name = files_nc(j).name;
    [var, ~, ~]=nc2mat([CTD_dir, name]);
    dind = regexp(name, '\d'); % gets all the digits for me.
    deployment(j) = str2num(name(dind(end-2:end)));
    CTD(j) = var;
    lon(j) = var.longitude;
    lat(j) = var.latitude;
    CTD(j).time = var.time/60/24+datenum(1900,1,1,0,0,0); % convert to datenum format
    max_press(j) = var.pressure(end);
    clear var name dind
end
clear CTD_dir files_nc nfiles j
%   struct with fields:
% 
%                    time: 6.2781e+07
%               woce_date: 20190515
%               woce_time: 1.5317e+04
%                pressure: [494×1 double]
%                latitude: -34.7073
%               longitude: 113.4239
%             temperature: [1×494 double]
%         temperatureFlag: [1×494 double]
%            conductivity: [1×494 double]
%        conductivityFlag: [1×494 double]
%                salinity: [1×494 double]
%            salinityFlag: [1×494 double]
%           temperature_2: [1×494 double]
%       temperature_2Flag: [1×494 double]
%          conductivity_2: [1×494 double]
%      conductivity_2Flag: [1×494 double]
%              salinity_2: [1×494 double]
%          salinity_2Flag: [1×494 double]
%                  oxygen: [1×494 double]
%              oxygenFlag: [1×494 double]
%                oxygen_2: [1×494 double]
%            oxygen_2Flag: [1×494 double]
%             fluorometer: [1×494 double]
%         fluorometerFlag: [1×494 double]
%                     par: [1×494 double]
%                 parFlag: [1×494 double]
%               altimeter: [1×494 double]
%           altimeterFlag: [1×494 double]
%         transmissometer: [1×494 double]
%     transmissometerFlag: [1×494 double]

%% Checking Station positions and pressure

figure(1)
clf
scatter(lon, lat, 140, max_press)
colorbar
caxis([400, 5500])
% So this figures shows that there are some stations away from 110E lines
% which for now I would like to omit from the analysis. 
% Therefore, it is important to account for the longitudes of the stations
% and I am interested in looking at full column stations only. This leads
% to the subsequent steps where, I seperate stations

%% Separate deep and shallow stations
% The station is considered in the vicinity of 110E if the actual
% measurement of a longitude deviate less than 1 degree longitudes. This is
% to accounting for the drift of ships during the package deployment. 

% This would gives us stations numbers that I care to seperate both shallow
% stations and deep stations
ind_deep = find(abs(lon-110)<1 & max_press>2000);
ind_shal = find(abs(lon-110)<1 & max_press<600);

% Checking the casts
figure(3);clf
plot(lat(ind_deep), deployment(ind_deep), 'sb')
hold on
text(lat(ind_deep), deployment(ind_deep), num2str(deployment(ind_deep)'), 'color', 'b')
plot(lat(ind_shal), deployment(ind_shal), '.k')
text(lat(ind_shal), deployment(ind_shal), num2str(deployment(ind_shal)'), 'color', 'k')
legend('Deep Casts', 'Shallow Casts', 'location', 'northwest')
hold off

% The below notes talks about number of profiles not the deployments
%----- Here 19 corresponds to 49th deployment and 20 corresponds to 48th--

%% NB - Deep casts: Notes from Helen and additional notes from Me
% Stations 19 and 20 were switched, 19 is at 11.5S and 20 is at 12.5S so
% need to sort in latitude before plotting
[~, ind] = sort(lat(ind_deep));
ind_deep = ind_deep(ind);

%NB - Shallow casts
% From the Above figure, I can see Stations 20-21 are at the same latitude - 12.5S because of CTD
% problem trouble shooting.
% Plot data to see which one to use. Last one had missing salinity above
% 30m. Last one is the one that water samples came from. Should use this
% one but use #20 since S is good to the surface
ind_shal = [ind_shal(1:19) ind_shal(21)]; % NOTE: I have used deployment 43 here, contrary to 42 for CTD.

% checking stations order for the data recording
figure(4);clf
plot(ind_deep, 'sb')
text(1:20, ind_deep, num2str(deployment(ind_deep)'), 'color', 'r')
hold on
plot(ind_shal, '.k')
text(1:20, ind_shal, num2str(deployment(ind_shal)'), 'color', 'k')
legend('Deep Casts', 'Shallow Casts', 'location', 'northwest')
hold off
%%
% station_110E=[1:20]; 
% n_stations = length(station_110E);
% station_lat=[-39.5:1.5:-12.5 -11.5];


% figure(3)
% clf
% for j=1:n_stations
%     subplot(1,3,1)
%     plot(CTD(ind_shal(j)).temperature,CTD(ind_shal(j)).pressure,'-','linew',1.5)
%     axis ij
%     xlabel('T (C)')
%     set(gca,'ylim',[0 500])
%     subplot(1,3,2)
%     plot(CTD(ind_shal(j)).salinity,CTD(ind_shal(j)).pressure,'-','linew',1.5)
%     axis ij
%     xlabel('S (psu)')
%     set(gca,'ylim',[0 500])
%     title(['Station ',int2str(station_110E(j)),' Lat ',num2str(CTD(ind_shal(j)).latitude)])
%     subplot(1,3,3)
%     plot(CTD(ind_shal(j)).oxygen,CTD(ind_shal(j)).pressure,'-','linew',1.5)
%     axis ij
%     xlabel('O2 (uM)')
%     set(gca,'ylim',[0 500])
%     pause
% end


figure(4)
clf
subplot(1,2,1) % deep
plot(lon,lat,'k.','markersize',8)
hold on 
plot(lon(ind_deep),lat(ind_deep),'ro','markersize',8)
text(lon(ind_deep)+0.05,lat(ind_deep), num2str(deployment(ind_deep)'), 'color', 'b')
title('110^{o}E Deep casts 4 Hydro')
hold off

subplot(1,2,2) % deep
plot(lon,lat,'k.','markersize',8)
hold on 
plot(lon(ind_shal),lat(ind_shal),'ro','markersize',8)
text(lon(ind_shal)+0.05,lat(ind_shal), num2str(deployment(ind_shal)'), 'color', 'k')
hold off
title('110^{o}E Shallow casts 4 Hydro')

% print(gcf, '-dpng', '-r100', '-painters', 'IN2019_V03SamplingHydro')
%% Seperating data first and saving them temporarily
deepdep = deployment(ind_deep);
% deep casts
ctd19v03deep = CTD(ind_deep);
for stn = 1:length(ind_deep)
    ctd19v03deep(stn).deployment = deepdep(stn);
end

% shallow casts
shaldep = deployment(ind_shal);
ctd19v03shal = CTD(ind_shal);
for stn = 1:length(ind_shal)
    ctd19v03shal(stn).deployment = shaldep(stn);
end

% Reorder the fields
% I first identify how many fields are there in a structure
% Then by looking the field name find that I need to move last added field
% to first 

aa = fieldnames(ctd19v03deep);
stlen = length(aa);
ctd19v03deep = orderfields(ctd19v03deep, [stlen,1:stlen-1]);

% likewise for shallow casts
aa = fieldnames(ctd19v03shal);
stlen = length(aa);
ctd19v03shal = orderfields(ctd19v03shal, [stlen,1:stlen-1]); clear stlen aa



% save deepshalhydrocast ctd19v03deep ctd19v03shal

% %% Cleaning shallow casts for further analysis
% clear; clc; close all
% load deepshalhydrocast.mat
% clear ctd19v03deep
% 
% % Reading data 
% ctd.lat = [ctd19v03shal.latitude]; 
% ctd.lon = [ctd19v03shal.longitude];
% ctd.p = [ctd19v03shal.pressure]; % dbar
% ctd.t = cat(1, ctd19v03shal.temperature)'; % degree C
% ctd.SP = cat(1, ctd19v03shal.salinity)'; % psu
% ctd.oxy = cat(1, ctd19v03shal.oxygen)'; % uM - micro mol/l
% ctd.nox = cat(1, ctd19v03shal.nox)'; % micro M where M = mol/l therefore umol/l
% ctd.nit = cat(1, ctd19v03shal.nitrite)'; % micro M
% ctd.pho = cat(1, ctd19v03shal.phosphate)'; % micro M
% ctd.sili = cat(1, ctd19v03shal.silicate)'; % micro M
% ctd.nh3 = cat(1, ctd19v03shal.ammonia)'; % micro M
% 
% % Retrieving sampling depths based on NaN at recorded property
% % saving sampling depths
% sdepth = ctd.p;
% sdepth(isnan(ctd.nox)) = NaN; % It is understood that all nutrients must be sampled at same depth. 
% ctd.sdepth = sdepth;
% 
% %% Checking how many bottles were sampled at each stations
% for stn = 1:length(ctd19v03shal)
%     dp = sdepth(:, stn);
%     dp = dp(~isnan(dp));
%     bott(1, stn) = length(dp); %#ok
%     clear dp
% end
% %
% figure(1);
% bar(bott)
% 
% % This figure indicates that bottle samples were not collected
% % consistently.
% %% Let's identify the pattern for bottle sampling locations
% for stn = 1:length(ctd19v03shal)
%     dp = sdepth(:, stn);
%     dp = dp(~isnan(dp));
%     sampledp{1,stn} = dp; %#ok
%     sampdiff{1, stn} = diff(round(dp));
%     clear dp
% end
% % creating column
% samdep = cell2col(sampledp);
% samdiff = cell2col(sampdiff);
% 
% % removing nans
% samdep(isnan(samdep)) = [];
% samdiff(isnan(samdiff)) = [];
% 
% % rounding and finding unique depths
% uqdepth = unique(round(samdep));
% 
% % plotting both to check the results
% figure(2);clf
% plot(round(samdep), '.k')
% hold on
% plot(uqdepth, 'or')
% hold off
% axis ij
% 
% % 
% figure(3);clf
% histogram(samdiff,'BinWidth', 5)
% %
% figure(4);clf
% plot(samdiff, 'o')
% 
% %% Interpolation grid levels for uniform representation
% udepth = [0:5:100, 110:10:200, 250:50:500];
% 
% % %% Interpolating all the stations to uniform depth intervals
% % % In this case I'm using GSW package to do so..
% aa = gsw_tracer_interp([ctd.nox],round([ctd.sdepth]), udepth);
% 
% figure(5);clf
% for stn = 1:20
%     plot(ctd.nox(:,stn), sdepth(:, stn), '-*k')
%     hold on
%     plot(aa(:,stn), udepth, 'or'); hold off; axis ij
%     pause()
% end

%% 
clear; clc; close all
load deepshalhydrocast.mat
clear ctd19v03deep

%
% Creating new structure for data
ctdHy19v03S = struct('Stations', {}, 'time', {}, 'lat', {}, 'lon', {}, 'sdepth', {}, 't', {}, 'SP', {}, ...
    'upres', {}, 'ct', {}, 'SA', {}, 'oxy', {}, 'nox', {}, 'no2', {}, 'phos', {}, 'silicate', {}, 'ammonia', {});

% Interpolation grid levels for uniform representation
% upres = [0:5:100, 110:10:200, 250:50:500]';
upres = [0:5:500]';

%----
for stn = 1:length(ctd19v03shal)
    disp(['Station', ' ', num2str(ctd19v03shal(stn).deployment)])
    ctdHy19v03S(stn).Stations = ctd19v03shal(stn).deployment; 
    ctdHy19v03S(stn).lat = ctd19v03shal(stn).latitude; 
    ctdHy19v03S(stn).lon = ctd19v03shal(stn).longitude;
    ctdHy19v03S(stn).time = ctd19v03shal(stn).time;
    % Bottle must be collected at same sampling depth
    p = ctd19v03shal(stn).pressure; n = ctd19v03shal(stn).nox';
    p(isnan(n)) = NaN;
    ctdHy19v03S(stn).sdepth = p(~isnan(p));
    ctdHy19v03S(stn).upres = upres;
    % writing biogeochemical properties at uniform depths
    ctdHy19v03S(stn).oxy = gsw_tracer_interp(ctd19v03shal(stn).oxygen', p, upres); % micromol l-1 uM
    ctdHy19v03S(stn).nox = gsw_tracer_interp(ctd19v03shal(stn).nox', p, upres); % micromol l-1 uM
    ctdHy19v03S(stn).no2 = gsw_tracer_interp(ctd19v03shal(stn).nitrite', p, upres); % micromol l-1 uM
    ctdHy19v03S(stn).phos = gsw_tracer_interp(ctd19v03shal(stn).phosphate', p, upres); % micromol l-1 uM
    ctdHy19v03S(stn).silicate = gsw_tracer_interp(ctd19v03shal(stn).silicate', p, upres); % micromol l-1 uM
    ctdHy19v03S(stn).ammonia = gsw_tracer_interp(ctd19v03shal(stn).ammonia', p, upres); % micromol l-1 uM

    % reading original data for physical properties
    ctdHy19v03S(stn).t = ctd19v03shal(stn).temperature'; % degree C
    ctdHy19v03S(stn).SP = ctd19v03shal(stn).salinity'; % psu
    % Converting to CT and SA and interpolating at uniform depths
    sa = gsw_SA_from_SP(ctd19v03shal(stn).salinity', p, ctd19v03shal(stn).longitude, ctd19v03shal(stn).latitude); % Absolute Salinity g/kg
    ct = gsw_CT_from_t(sa, ctd19v03shal(stn).temperature', p); % Conservative Temperature
    % projecting all on the uniform depth
    [SA, CT] = gsw_SA_CT_interp(sa, ct, p, upres);
    %
    ctdHy19v03S(stn).ct = CT; % degree C
    ctdHy19v03S(stn).SA = SA; % g/kg
   
    clear p sa ct SA CT 
end
clear stn

%% checking data
figure(5);clf
for stn = 1:15
    plot(ctd19v03shal(stn).nox, ctd19v03shal(stn).pressure, '*k')
    hold on
    plot(ctdHy19v03S(stn).nox, ctdHy19v03S(stn).upres, 'or'); hold off; axis ij
    ylim([0, 600])
    pause()
end

%% Saving data
% save IN19V03HyS ctdHy19v03S

%% Preparing Deep Hydro Cast for the analysis

clear; clc; close all
load deepshalhydrocast.mat
clear ctd19v03shal

%
% Creating new structure for data
ctdHy19v03D = struct('Stations', {}, 'time', {}, 'lat', {}, 'lon', {}, 'sdepth', {}, 't', {}, 'SP', {}, ...
    'upres', {}, 'ct', {}, 'SA', {}, 'dynhgt', {},'oxy', {}, 'nox', {}, 'no2', {}, 'phos', {}, 'silicate', {}, 'ammonia', {});

% Interpolation grid levels for uniform representation
upres = [0:5:100, 110:10:200, 250:50:500, 600:100:2000, 2500:500:4500]';

%----
for stn = 1:length(ctd19v03deep)
    disp(['Station', ' ', num2str(ctd19v03deep(stn).deployment)])
    ctdHy19v03D(stn).Stations = ctd19v03deep(stn).deployment; 
    ctdHy19v03D(stn).lat = ctd19v03deep(stn).latitude; 
    ctdHy19v03D(stn).lon = ctd19v03deep(stn).longitude;
    ctdHy19v03D(stn).time = ctd19v03deep(stn).time;
    % Bottle must be collected at same sampling depth
    p = ctd19v03deep(stn).pressure; n = ctd19v03deep(stn).nox';
    p(isnan(n)) = NaN;
    ctdHy19v03D(stn).sdepth = p(~isnan(p));
    ctdHy19v03D(stn).upres = upres;
    % writing biogeochemical properties at uniform depths
    ctdHy19v03D(stn).oxy = gsw_tracer_interp(ctd19v03deep(stn).oxygen', p, upres); % micromol l-1 uM
    ctdHy19v03D(stn).nox = gsw_tracer_interp(ctd19v03deep(stn).nox', p, upres); % micromol l-1 uM
    ctdHy19v03D(stn).no2 = gsw_tracer_interp(ctd19v03deep(stn).nitrite', p, upres); % micromol l-1 uM
    ctdHy19v03D(stn).phos = gsw_tracer_interp(ctd19v03deep(stn).phosphate', p, upres); % micromol l-1 uM
    ctdHy19v03D(stn).silicate = gsw_tracer_interp(ctd19v03deep(stn).silicate', p, upres); % micromol l-1 uM
    ctdHy19v03D(stn).ammonia = gsw_tracer_interp(ctd19v03deep(stn).ammonia', p, upres); % micromol l-1 uM

    % reading original data for physical properties
    ctdHy19v03D(stn).t = ctd19v03deep(stn).temperature'; % degree C
    ctdHy19v03D(stn).SP = ctd19v03deep(stn).salinity'; % psu
    % Converting to CT and SA and interpolating at uniform depths
    sa = gsw_SA_from_SP(ctd19v03deep(stn).salinity', p, ctd19v03deep(stn).longitude, ctd19v03deep(stn).latitude); % Absolute Salinity g/kg
    ct = gsw_CT_from_t(sa, ctd19v03deep(stn).temperature', p); % Conservative Temperature
    % projecting all on the uniform depth
    [SA, CT] = gsw_SA_CT_interp(sa, ct, p, upres);
    %
    ctdHy19v03D(stn).ct = CT; % degree C
    ctdHy19v03D(stn).SA = SA; % g/kg
    % adding Dynamic Height
    dyheight = gsw_geo_strf_dyn_height(SA, CT, upres, 2000); % in $m^2/s^2$
    ctdHy19v03D(stn).dynhgt = dyheight;
    
    clear p sa ct SA CT dyheight 
end
clear stn

%% checking data
figure(5);clf
for stn = 1:15
    plot(ctd19v03deep(stn).nox, ctd19v03deep(stn).pressure, '*k')
    hold on
    plot(ctdHy19v03D(stn).nox, ctdHy19v03D(stn).upres, 'or'); hold off; axis ij
    ylim([0, 600])
    pause()
end

%% Saving data
% save IN19V03HyD ctdHy19v03D
