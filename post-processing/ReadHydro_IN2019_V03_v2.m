% This script I read Hydrographic data of IN2019_V03 voyage along 110E line
% This is analogous to ReadyHydro_IN2019_V03.m with subtle difference in
% vertical interpolation. 
% NOTE: Here I will use gsw_tracer_CT_interp instead of gsw_tracer_interp function.

% The data that used in this script were obtained from Helen, currently at
% OneDrive>don't sync>IN2019_V03 Analysis>processed Hydrography

%% 

% Recreating Shallow cast interpolated at uniform pressure levels using
% gsw_tracer_CT_interp function. The default factor to scale the slop is 9.
% I am using this function because it preserve the shape of the tracer-CT
% diagram, thereby, removes anomalous watermass properties. 

% NOTE FROM FUNCTION DISCRIPTION
% The interpolation method is is the 
% MRST-PCHIP method described in Barker and McDougall (2020), with the
% tracer data being used in place of Absoluate Salinity data. 
% 
% This function requires scaling the tracer and temperature data so that
% the tracer-CT diagram reflects the relative variation of the tracer and 
% temperature in the world ocean.  Specifically, "factor" should be chosen
% to be the ratio of the global range of CT to that of the tracer variable
% in the world ocean.  A list of suitable values of "factor" for various
% tracers is given here. 
% 
%       TRACER              UNITS             FACTOR
%     Absolute Salinity     g/kg                9
%     dissolved oxygen        ?                 ?
%     AOU                     ?                 ?
%     silicic acid            ?                 ?
%     nitrate                 ?                 ?
%     phosphate               ?                 ?
%     carbon 14               ?                 ?
%     tritium                 ?                 ?
%     eastward velocity      m/s               100
%     westward velocity      m/s               100
% 
% If an input value of "factor" is not given in the function call, it is 
% set equal to 9.
% 
% Any interpolated bottles that have pressures shallower than the 
% shallowest observed bottle are set equal to the shallowest observed 
% bottle.

% NOTE FROM gsw_tracer_interp
% When the tracer is in-situ temperature we have found a suitable value 
% for the scale_factor is 0.33, so that the final scaling factor is 0.33
% times the maximum magnitude (over all data pairs) of the slope on the 
% [tracer - bottle_number] diagram.  We expect 0.33 will be a suitable 
% scale_factor for other tracer data.

%% Reading hydrographic data from ReadHydro_IN2019_V03.m commented section

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
    
    % writing biogeochemical properties at uniform depths
    scale_factor = 0.33; % see NOTE 
    ctdHy19v03S(stn).oxy = gsw_tracer_CT_interp(ctd19v03shal(stn).oxygen', ct, p, upres, scale_factor); % micromol l-1 uM
    ctdHy19v03S(stn).nox = gsw_tracer_CT_interp(ctd19v03shal(stn).nox', ct, p, upres, scale_factor); % micromol l-1 uM
    ctdHy19v03S(stn).no2 = gsw_tracer_CT_interp(ctd19v03shal(stn).nitrite', ct, p, upres, scale_factor); % micromol l-1 uM
    ctdHy19v03S(stn).phos = gsw_tracer_CT_interp(ctd19v03shal(stn).phosphate', ct, p, upres, scale_factor); % micromol l-1 uM
    ctdHy19v03S(stn).silicate = gsw_tracer_CT_interp(ctd19v03shal(stn).silicate', ct, p, upres, scale_factor); % micromol l-1 uM
    ctdHy19v03S(stn).ammonia = gsw_tracer_CT_interp(ctd19v03shal(stn).ammonia', ct, p, upres, scale_factor); % micromol l-1 uM
   
    clear p sa ct SA CT 
end
clear stn

%% checking data
figure(5);clf
for stn = 1:15
    plot(ctd19v03shal(stn).temperature, ctd19v03shal(stn).pressure, '*k')
    hold on
    plot(ctdHy19v03S(stn).ct, ctdHy19v03S(stn).upres, 'or'); hold off; axis ij
    ylim([0, 600])
    pause()
end

%% Saving data
ctdHy19v03Sv2 = ctdHy19v03S;
% save IN19V03HySv2 ctdHy19v03Sv2
clear ctdHy19v03S
%% Comparing performance of old version and v2 
load IN19V03HyS.mat

figure(5);clf
for stn = 1:15
    plot(ctdHy19v03Sv2(stn).nox, ctdHy19v03Sv2(stn).upres, '*k')
    hold on
    plot(ctdHy19v03S(stn).nox, ctdHy19v03S(stn).upres, 'or'); hold off; axis ij
    ylim([0, 600])
    legend('with CT', 'without CT')
    pause()
end
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
% Check notes at the begining of the script
scale_factor = 0.33;% 

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
    
    % writing biogeochemical properties at uniform depths
    ctdHy19v03D(stn).oxy = gsw_tracer_CT_interp(ctd19v03deep(stn).oxygen', ct, p, upres, scale_factor);% micromol l-1 uM
    ctdHy19v03D(stn).nox = gsw_tracer_CT_interp(ctd19v03deep(stn).nox', ct, p, upres, scale_factor); % micromol l-1 uM
    ctdHy19v03D(stn).no2 = gsw_tracer_CT_interp(ctd19v03deep(stn).nitrite', ct, p, upres, scale_factor); % micromol l-1 uM
    ctdHy19v03D(stn).phos = gsw_tracer_CT_interp(ctd19v03deep(stn).phosphate', ct, p, upres, scale_factor); % micromol l-1 uM
    ctdHy19v03D(stn).silicate = gsw_tracer_CT_interp(ctd19v03deep(stn).silicate', ct, p, upres, scale_factor); % micromol l-1 uM
    ctdHy19v03D(stn).ammonia = gsw_tracer_CT_interp(ctd19v03deep(stn).ammonia', ct, p, upres, scale_factor); % micromol l-1 uM

    % adding Dynamic Height
    dyheight = gsw_geo_strf_dyn_height(SA, CT, upres, 2000); % in $m^2/s^2$
    ctdHy19v03D(stn).dynhgt = dyheight;
    
    clear p sa ct SA CT dyheight n
end
clear stn

%% Saving data
ctdHy19v03Dv2 = ctdHy19v03D; 
% save IN19V03HyDv2 ctdHy19v03Dv2

%% comparing performance of deeper level cast with ct and without ct
load IN19V03HyD.mat

figure(5);clf
for stn = 1:15
    plot(ctdHy19v03Dv2(stn).nox, ctdHy19v03Dv2(stn).upres, '*k')
    hold on
    plot(ctdHy19v03D(stn).nox, ctdHy19v03D(stn).upres, 'or'); hold off; axis ij
    ylim([0, 6000])
    legend('with CT', 'without CT')
    pause()
end