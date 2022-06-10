% Reading data
clear; clc; %close all
%
addpath ~/'OneDrive - University of Tasmania'/MatlabData/HPdata/IOVoyage/

load IN19V03D.mat
%
lon = [ctd19v03D.lon];
lat = [ctd19v03D.lat];
stn = [ctd19v03D.Deployment];
% Computing depths
pres = [ctd19v03D.upres]; % as all the data are projected on uniform pressure
z = gsw_z_from_p(pres, lat);
depth = gsw_depth_from_z(z); clear z
%
pres = pres(:, 1); pind10 = pres >= 10;

% For now I am focussing on first sensor only.. 
% computing physical properties
sa = gsw_SA_from_SP([ctd19v03D.SP], pres, lon, lat); % Absolute Salinity g/kg
ct = gsw_CT_from_t(sa, [ctd19v03D.t], pres); % Conservative Temperature
IN19pden = gsw_rho(sa, ct, 10.1325); % Potential Density 10.1325 dbar considered as reference level

% Computing geostrophic velocity from CTD data 
% stations are arranged in South to north. Thus, +ve values would suggest
% westward flow and -ve would suggest eastward flow.
dyheight = gsw_geo_strf_dyn_height(sa, ct, pres, 2000); % in $m^2/s^2$
[IN19GeoVel, IN19mid_lat, IN19mid_long] = gsw_geostrophic_velocity(dyheight, lon, lat, pres);

% Brunt-Vaisala frequency
[lats, ps] = meshgrid(lat, pres);
%
[IN19N2, IN19p_mid] = gsw_Nsquared(sa, ct, ps, lats);
% PV Helen suggested...
[IPV_vs_fN2, p_mid] = gsw_IPV_vs_fNsquared_ratio(sa, ct, ps);
[lats, ~] = meshgrid(lat, p_mid(:, 1));
f = gsw_f(lats);
IN19PV2 = f.*IN19N2.*IPV_vs_fN2;

for ii = 1:size(IN19PV2, 2)
    IN19PV2sm(:, ii) = smooth(IN19PV2(:, ii), 50); %#ok
end

% Neutral density
nctd = size(ct, 2);
for ii = 1:nctd
    [g_n, ~, ~] = eos80_legacy_gamma_n(sa(:, ii), ct(:, ii), pres, lon(ii), lat(ii));
    IN19gm_n(:, ii) = g_n; %#ok
    % Computing PV
    Prho = IN19pden(:, ii);
    depth = gsw_depth_from_z(gsw_z_from_p(pres, lat(ii)));
    [pv, z] = ra_pVorticityCTD(Prho, depth, lat(ii));
%     pv = smooth(pv, 100);
    IN19PV1(:, ii) = pv;%#ok
    IN19Z(:, ii) = z; %#ok
    % mixed layer depth
    mld = ra_CTDmld(sa(pind10, ii), ct(pind10, ii), pres(pind10), 0.2);
    IN19mld(:, ii) = mld; %#ok
    clear g_n mld Prho pv z depth
end
% PV computation following SW package methodology..
[lats, ps] = meshgrid(lat, IN19p_mid(:,  1));
g = gsw_grav(lats, ps);
IN19PV = (f.*IN19N2)./g;

clear ii lat lon pind10 ii lats ps f p_mid nctd g

%% Mixed layer depth and potential vorticity

close all
FX = 23.4; FY = 15; % page size A4
FIG1 = figure('units', 'centimeters', 'Position', [0 0 FX FY]);
set(FIG1, 'PaperUnits', 'centimeters')
set(FIG1, 'Papersize', [FX FY])
set(FIG1, 'PaperPositionMode', 'manual');
set(FIG1, 'PaperPosition', [0 0 FX FY]);

% Axis set up for the plots
dy = 0.8; 
width = 18; height = 12;
left = [2.2, width+2.2+0.5];
bottom = [1.2, 1.2+dy+height, 1.2+dy*2+height*2];

%
xtick = 110; ytick = -35:5:-10; 
[~, yticklabel] = XYTickLabel(xtick, ytick);
ylabelname = 'Pressure [dbar]';
yst = 0; yint = 100; yen = 1500;

% y = repmat([ctd19v03D.lat], size(IN19PV, 1), 1);
figure(1);clf
salint = 0:0.5:4.5;
axes('units', 'centimeters', 'position', [left(1) bottom(1) width height]);
contourf([ctd19v03D.lat], IN19p_mid(:, 1), -1*IN19PV2sm*10^9, salint, 'linest', 'none');
hold on
[c, h] = contour([ctd19v03D.lat], pres, IN19gm_n, [ 26.6, 26.9, 27, 27.2, 27.4, 27.6, 27.8], 'linest', '-', 'color', 'b', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'b', 'fontweigh', 'bold'); clear c h
% mld 
plot([ctd19v03D.lat], IN19mld, '-r', 'linewi', 2)
hold off
cmocean('turbid', length(salint)-1); caxis([salint(1), salint(end)])
xlim([-35, -10]); ylim([0, 1500])
ra_subsurfembellish('xy', yst, yint, yen, ytick, yticklabel, ylabelname)

% colorbar for two figures.. 
ra_colorbarV2('vertical', 'centimeters', [left(2) bottom(1) 0.45 height], 'Potential Vorticity [ x 10^{-9} m^{-1}s^{-1}]', 'in')

% print(gcf, '-dpng', '-r300', '-painters', 'IN19mNPVmld_PFig')

%% Plotting Geostrophic Velocity

close all
FX = 24.5; FY = 15; % page size A4
FIG1 = figure('units', 'centimeters', 'Position', [0 0 FX FY]);
set(FIG1, 'PaperUnits', 'centimeters')
set(FIG1, 'Papersize', [FX FY])
set(FIG1, 'PaperPositionMode', 'manual');
set(FIG1, 'PaperPosition', [0 0 FX FY]);
%
xtick = 110; ytick = -35:5:-10; 
[xticklabel, yticklabel] = XYTickLabel(xtick, ytick);
ylabelname = 'Pressure [dbar]';
yst = 0; yint = 100; yen = 1500;
%
% Axis set up for the plots
dy = 1.5; 
width = 18; height = 12;
left = [2.5, width+2.5+0.5];
bottom = [1.2, 1.2+dy+height, 1.2+dy*2+height*2];
velint = -0.2:0.02:0.2;

% Remapping geostrophic velocity to original latitudes..
% [lat, p] = meshgrid([ctd19v03D.lat], pres);
% govel = interp2(IN19mid_lat,repmat(pres, 1, 19), IN19GeoVel, lat, p, 'extrap');
figure(1);clf
axes('units', 'centimeters', 'position', [left(1) bottom(1) width height]);
contourf(IN19mid_lat,repmat(pres, 1, 19), -IN19GeoVel, velint, 'linest', 'none');
hold on
% mld 
%plot([ctd19v03D.lat], IN19mld, '-r', 'linewi', 2)
[c, h] = contour([ctd19v03D.lat], pres, IN19gm_n, [ 26.6, 26.9, 27, 27.2, 27.4, 27.6, 27.8], 'linest', '-', 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 16, 'color', 'k', 'fontweigh', 'bold'); clear c h
hold off
cmocean('balance', length(velint)-1); caxis([velint(1), velint(end)])
xlim([-40, -10]); ylim([0, 1500])
% title('IN2019\_V03 CTD observations ({\gamma}_n)')
ra_subsurfembellish('xy', yst, yint, yen, ytick, yticklabel, ylabelname)
set(gca, 'XMinorTick', 'on','TickLength', [0.02, 0.02], 'Xcolor', [0.3, 0.3, 0.3], 'Ycolor', [0.3, 0.3, 0.3], 'box', 'on')
ra_colorbarV2('vertical', 'centimeters', [left(2)-0.4, bottom(1) 0.45 height], 'Geostrophic Volocity [ ms^{-1}]', 'in')

% axes('units', 'centimeters', 'position', [left(1) bottom(2) width 5]);
% plot([ctd19v03D.lat], dyheight(100, :), '-o', 'color', rgb('orange'), 'linewi', 1.5)
% grid on
% xlim([-40, -10]);
% set(gca, 'XMinorTick', 'on', 'TickDir', 'out', 'TickLength', [0.02, 0.02], 'Xcolor', [0.3, 0.3, 0.3], 'Ycolor', [0.3, 0.3, 0.3], 'box', 'on')
% ntitle({'','100 dbar/2000 dbar'}, 'FontSize', 16, 'Fontweigh', 'bold')
% ra_subsurfembellish('y', 16, 1, 20, ytick, yticklabel, 'D_h [m^2/s^2]')

% print(gcf, '-dpng', '-r300', '-painters', 'IN19V03_GeovelCTDmldGN')
% jprint('./', 'IN19V03_GeovelCTDmldGN','-dpng', '-r300','-painters')
%% Conservative Temperature for the entire water column
left = 0.08; bottom = [0.55, 0.05];
width = 0.8; height = [0.35, 0.5];
%
xtick = 110; ytick = -35:5:-10; 
[~, yticklabel] = XYTickLabel(xtick, ytick);
ylabelname = 'Pressure [dbar]';
%
x = [ctd19v03D.lat];
clear z*
% surface to 500 m
lev = pres <= 500;
z02 = ct(lev, :);
y02 = pres(lev); %r02 = prho(lev, :);
% 500m to 5000 m
lev215 = pres > 500 & pres <= 5000;
z215 = ct(lev215, :); %z215 = z215(:, ix);
y215 = pres(lev215); %r215 = prho(lev215, :); %r215 = r215(:, ix);

%
varint = 0:1:30; cval = [1:3:30];
figure(3);clf
axes('position', [left, bottom(1), width, height(1)]);
contourf(x, y02, z02, varint, 'LineStyle','none');
hold on
[c, h] = contour(x, y02, z02, cval, 'linest', '-', 'color', 'w', 'linewi', 2);
clabel(c, h, 'fontsize', 16, 'color', 'w', 'fontweigh', 'bold'); clear c h
% 
hold off
cmocean('thermal', length(varint)-1); caxis([varint(1), varint(end)])
%
% title('IN2019\_V03, 110^{\circ}E')
%
yst = 0; yint= 100; yen = 495;
xlim([-40, -10]); ylim([yst, yen])
axis tight
ra_subsurfembellish('y', yst, yint, yen, ytick, yticklabel, ylabelname)

% 500 to 5000 dbar
axes('position', [left, bottom(2), width, height(2)]);
contourf(x, y215, z215, varint, 'LineStyle','none');
hold on
[c, h] = contour(x, y215, z215, cval, 'linest', '-', 'color', 'w', 'linewi', 2);
clabel(c, h, 'fontsize', 16, 'color', 'w', 'fontweigh', 'bold'); clear c h
hold off
cmocean('thermal', length(varint)-1); caxis([varint(1), varint(end)])
%
ra_colorbarV2('vertical', 'normalized', [0.885, bottom(2)+0.02, 0.015, 0.8], 'Conservative Temperature [^{\circ}C]', 'in')
yst = 500; yint= 500; yen = 5000;
xlim([-40, -10]); ylim([yst, yen])
axis tight
ra_subsurfembellish('xy', yst, yint, yen, ytick, yticklabel, ylabelname)

jprint('./', 'IN2019V03D_ct_ctd', '-r300', '-painters')
clear varint cval
%% Absolute Salinity
clear z*
% surface to 500 m
lev = pres <= 500;
z02 = sa(lev, :);
y02 = pres(lev); %r02 = prho(lev, :);
% 500m to 5000 m
lev215 = pres > 500 & pres <= 5000;
z215 = sa(lev215, :); %z215 = z215(:, ix);
y215 = pres(lev215); %r215 = prho(lev215, :); %r215 = r215(:, ix);

%
varint = 34.3:0.1:35.9; cval = [34.45, 34.75, 34.9, 35, 35.5, 35.8];
figure(3);clf
axes('position', [left, bottom(1), width, height(1)]);
contourf(x, y02, z02, varint, 'LineStyle','none');
hold on
[c, h] = contour(x, y02, z02, cval([1:2, 4:end]), 'linest', '-', 'color', 'w', 'linewi', 2);
clabel(c, h, 'fontsize', 16, 'color', 'w', 'fontweigh', 'bold', 'LabelSpacing', 400); clear c h
% 
hold off
cmocean('haline', length(varint)-1); caxis([varint(1), varint(end)])
%
% title('IN2019\_V03, 110^{\circ}E')
%
yst = 0; yint= 100; yen = 495;
xlim([-35, -10]); ylim([yst, yen])
axis tight
ra_subsurfembellish('y', yst, yint, yen, ytick, yticklabel, ylabelname)

% 500 to 5000 dbar
axes('position', [left, bottom(2), width, height(2)]);
contourf(x, y215, z215, varint, 'LineStyle','none');
hold on
[c, h] = contour(x, y215, z215, cval, 'linest', '-', 'color', 'w', 'linewi', 2);
clabel(c, h, 'fontsize', 16, 'color', 'w', 'fontweigh', 'bold', 'LabelSpacing', 500); clear c h
hold off
cmocean('haline', length(varint)-1); caxis([varint(1), varint(end)])
%
ra_colorbarV2('vertical', 'normalized', [0.885, bottom(2)+0.02, 0.015, 0.8], 'Absolute Salinity [g/kg]', 'in')
yst = 500; yint= 500; yen = 5000;
xlim([-35, -10]); ylim([yst, yen])
axis tight
ra_subsurfembellish('xy', yst, yint, yen, ytick, yticklabel, ylabelname)

 jprint('./', 'IN2019V03D_sa_ctd', '-r300', '-painters')
clear varint cval
%% Potential density 
clear z*
% surface to 500 m
lev = pres <= 500;
z02 = IN19pden(lev, :)-1000;
y02 = pres(lev); %r02 = prho(lev, :);
% 500m to 5000 m
lev215 = pres > 500 & pres <= 5000;
z215 = IN19pden(lev215, :)-1000; %z215 = z215(:, ix);
y215 = pres(lev215); %r215 = prho(lev215, :); %r215 = r215(:, ix);

%
varint =  22:0.2:28; cval = [22, 24, 25, 26, 26.6, 26.8, 26.9, 27.2:0.2:28];
figure(3);clf
axes('position', [left, bottom(1), width, height(1)]);
contourf(x, y02, z02, varint, 'LineStyle','none');
hold on
[c, h] = contour(x, y02, z02, cval, 'linest', '-', 'color', 'w', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'w', 'fontweigh', 'bold', 'LabelSpacing', 400); clear c h
% 
hold off
cmocean('dense', length(varint)-1); caxis([varint(1), varint(end)])
%
title('IN2019\_V03, 110^{\circ}E')
%
yst = 0; yint= 100; yen = 495;
xlim([-35, -10]); ylim([yst, yen])
axis tight
ra_subsurfembellish('y', yst, yint, yen, ytick, yticklabel, ylabelname)

% 500 to 5000 dbar
axes('position', [left, bottom(2), width, height(2)]);
contourf(x, y215, z215, varint, 'LineStyle','none');
hold on
[c, h] = contour(x, y215, z215, cval, 'linest', '-', 'color', 'w', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'w', 'fontweigh', 'bold', 'LabelSpacing', 500); clear c h
hold off
cmocean('dense', length(varint)-1); caxis([varint(1), varint(end)])
%
ra_colorbarV2('vertical', 'normalized', [0.885, bottom(2)+0.02, 0.015, 0.8], 'Potential Density [kg/m^3]', 'in')
yst = 500; yint= 500; yen = 5000;
xlim([-35, -10]); ylim([yst, yen])
axis tight
ra_subsurfembellish('xy', yst, yint, yen, ytick, yticklabel, ylabelname)

% jprint('./', 'IN2019V03D_pden_ctd', '-r300', '-painters')
clear varint cval
%% Plotting hydrography Nutrient data.. 
clear;clc
load IN19V03HyDv2.mat
ctdHy19v03D = ctdHy19v03Dv2;
ctdHy19v03D(18)=[]; % Removed station 35 at index 18 due to sampling below 250 m. 
%
% Interpolating stations between 11S and 39.5;
ctdHy19v03Dlat = ra_projstations(ctdHy19v03D, -39.5, -11, 0.1, 1903);
%
% computing potential density
prho = gsw_rho(ctdHy19v03Dlat.SA, ctdHy19v03Dlat.ct, 10.1325) - 1000; % Potential Density 10.1325 dbar considered as reference level
%
sa = ctdHy19v03Dlat.SA;
ct = ctdHy19v03Dlat.ct;
lat = ctdHy19v03Dlat.lat;
lon = ctdHy19v03Dlat.lon;
pres = ctdHy19v03Dlat.upres;

% Neutral density
nctd = length(sa);
for ii = 1:nctd
    [g_n, ~, ~] = eos80_legacy_gamma_n(sa(:, ii), ct(:, ii), pres, lon(ii), lat(ii));
    gamma_n(:, ii) = g_n; %#ok
    clear g_n
end
clear ii

%% Oxygen
left = 0.08; bottom = [0.55, 0.05];
width = 0.8; height = [0.35, 0.5];
%
xtick = 110; ytick = -35:5:-10; 
[~, yticklabel] = XYTickLabel(xtick, ytick);
ylabelname = 'Pressure [dbar]';
%
x = lat;
var = [ctdHy19v03Dlat.oxy];
%
clear z*
% surface to 500 m
lev = pres <= 500;
z02 = var(lev, :);
y02 = pres(lev); %r02 = prho(lev, :);
% 500m to 5000 m
lev215 = pres > 500 & pres <= 5000;
z215 = var(lev215, :); %z215 = z215(:, ix);
y215 = pres(lev215); %r215 = prho(lev215, :); %r215 = r215(:, ix);

%
varint = 60:5:280; cval = [100, 150, 180, 240, 260];
figure(3);clf
axes('position', [left, bottom(1), width, height(1)]);
contourf(x, y02, z02, varint, 'LineStyle','none');
hold on
for stn = 1:length(ctdHy19v03D)
    ra_plotsampleloc(ctdHy19v03D(stn).lat, ctdHy19v03D(stn).sdepth, 12, [0.6, 0.6, 0.6])
end
[c, h] = contour(x, y02, z02, cval, 'linest', '-', 'color', 'w', 'linewi', 2);
clabel(c, h, 'fontsize', 16, 'color', 'w', 'fontweigh', 'bold'); clear c h
% 
hold off
cmocean('oxy', length(varint)-1); caxis([varint(1), varint(end)])
%
% title('IN2019\_V03, 110^{\circ}E')
%
yst = 0; yint= 100; yen = 495;
xlim([-35, -10]); ylim([yst, yen])
axis tight
ra_subsurfembellish('y', yst, yint, yen, ytick, yticklabel, ylabelname)

% 500 to 5000 dbar
axes('position', [left, bottom(2), width, height(2)]);
contourf(x, y215, z215, varint, 'LineStyle','none');
hold on
for stn = 1:length(ctdHy19v03D)
    ra_plotsampleloc(ctdHy19v03D(stn).lat, ctdHy19v03D(stn).sdepth, 12, [0.6, 0.6, 0.6])
end
[c, h] = contour(x, y215, z215, cval, 'linest', '-', 'color', 'w', 'linewi', 2);
clabel(c, h, 'fontsize', 16, 'color', 'w', 'fontweigh', 'bold'); clear c h
hold off
cmocean('oxy', length(varint)-1); caxis([varint(1), varint(end)])
%
ra_colorbarV2('vertical', 'normalized', [0.885, bottom(2)+0.02, 0.015, 0.8], 'Dissolved Oxygen [{\mu}mol L^{-1}]', 'in')
yst = 500; yint= 500; yen = 5000;
xlim([-35, -10]); ylim([yst, yen])
axis tight
ra_subsurfembellish('xy', yst, yint, yen, ytick, yticklabel, ylabelname)

 jprint('./', 'IN2019V03D_oxy_hydro', '-r300', '-painters')
clear varint cval
%% Nitrate
x = lat;
var = [ctdHy19v03Dlat.nox];
%
clear z*
% surface to 500 m
lev = pres <= 500;
z02 = var(lev, :);
y02 = pres(lev); %r02 = prho(lev, :);
% 500m to 5000 m
lev215 = pres > 500 & pres <= 5000;
z215 = var(lev215, :); %z215 = z215(:, ix);
y215 = pres(lev215); %r215 = prho(lev215, :); %r215 = r215(:, ix);
%
varint = 0:1:40; cval = [1, 5:5:40, 34, 36];
figure(3);clf
axes('position', [left, bottom(1), width, height(1)]);
contourf(x, y02, z02, varint, 'LineStyle','none');
hold on
for stn = 1:length(ctdHy19v03D)
    ra_plotsampleloc(ctdHy19v03D(stn).lat, ctdHy19v03D(stn).sdepth, 12, [0.6, 0.6, 0.6])
end
[c, h] = contour(x, y02, z02, cval, 'linest', '-', 'color', 'w', 'linewi', 2);
clabel(c, h, 'fontsize', 16, 'color', 'w', 'fontweigh', 'bold'); clear c h
% 
hold off
cmocean('matter', length(varint)-1); caxis([varint(1), varint(end)])
%
% title('IN2019\_V03, 110^{\circ}E')
%
yst = 0; yint= 100; yen = 495;
xlim([-35, -10]); ylim([yst, yen])
axis tight
ra_subsurfembellish('y', yst, yint, yen, ytick, yticklabel, ylabelname)

% 500 to 5000 dbar
axes('position', [left, bottom(2), width, height(2)]);
contourf(x, y215, z215, varint, 'LineStyle','none');
hold on
for stn = 1:length(ctdHy19v03D)
    ra_plotsampleloc(ctdHy19v03D(stn).lat, ctdHy19v03D(stn).sdepth, 12, [0.6, 0.6, 0.6])
end
[c, h] = contour(x, y215, z215, cval, 'linest', '-', 'color', 'w', 'linewi', 2);
clabel(c, h, 'fontsize', 16, 'color', 'w', 'fontweigh', 'bold'); clear c h
hold off
cmocean('matter', length(varint)-1); caxis([varint(1), varint(end)])
%
ra_colorbarV2('vertical', 'normalized', [0.885, bottom(2)+0.02, 0.015, 0.8], 'Nitrate [{\mu}mol L^{-1}]', 'in')
yst = 500; yint= 500; yen = 5000;
xlim([-35, -10]); ylim([yst, yen])
axis tight
ra_subsurfembellish('xy', yst, yint, yen, ytick, yticklabel, ylabelname)
 jprint('./', 'IN2019V03D_nox_hydro', '-r300', '-painters')
%% Silicate
x = lat;
var = [ctdHy19v03Dlat.silicate];
%
clear z*
% surface to 500 m
lev = pres <= 500;
z02 = var(lev, :);
y02 = pres(lev); %r02 = prho(lev, :);
% 500m to 5000 m
lev215 = pres > 500 & pres <= 5000;
z215 = var(lev215, :); %z215 = z215(:, ix);
y215 = pres(lev215); %r215 = prho(lev215, :); %r215 = r215(:, ix);
%
varint = 0:2:135; cval = [5:10:40, 50:20:130];
figure(3);clf
axes('position', [left, bottom(1), width, height(1)]);
contourf(x, y02, z02, varint, 'LineStyle','none');
hold on
for stn = 1:length(ctdHy19v03D)
    ra_plotsampleloc(ctdHy19v03D(stn).lat, ctdHy19v03D(stn).sdepth, 12, [0.6, 0.6, 0.6])
end
[c, h] = contour(x, y02, z02, cval, 'linest', '-', 'color', 'w', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'w', 'fontweigh', 'bold'); clear c h
% 
hold off
cmocean('matter', length(varint)-1); caxis([varint(1), varint(end)])
%
title('IN2019\_V03, 110^{\circ}E')
%
yst = 0; yint= 100; yen = 495;
xlim([-35, -10]); ylim([yst, yen])
axis tight
ra_subsurfembellish('y', yst, yint, yen, ytick, yticklabel, ylabelname)

% 500 to 5000 dbar
axes('position', [left, bottom(2), width, height(2)]);
contourf(x, y215, z215, varint, 'LineStyle','none');
hold on
for stn = 1:length(ctdHy19v03D)
    ra_plotsampleloc(ctdHy19v03D(stn).lat, ctdHy19v03D(stn).sdepth, 12, [0.6, 0.6, 0.6])
end
[c, h] = contour(x, y215, z215, cval, 'linest', '-', 'color', 'w', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'w', 'fontweigh', 'bold'); clear c h
hold off
cmocean('matter', length(varint)-1); caxis([varint(1), varint(end)])
%
ra_colorbarV2('vertical', 'normalized', [0.885, bottom(2)+0.02, 0.015, 0.8], 'Silicate [{\mu}mol L^{-1}]', 'in')
yst = 500; yint= 500; yen = 5000;
xlim([-35, -10]); ylim([yst, yen])
axis tight
ra_subsurfembellish('xy', yst, yint, yen, ytick, yticklabel, ylabelname)
% jprint('./', 'IN2019V03D_sil_hydro', '-r300', '-painters')

%% Phosphate
x = lat;
var = [ctdHy19v03Dlat.phos];
%
clear z*
% surface to 500 m
lev = pres <= 500;
z02 = var(lev, :);
y02 = pres(lev); %r02 = prho(lev, :);
% 500m to 5000 m
lev215 = pres > 500 & pres <= 5000;
z215 = var(lev215, :); %z215 = z215(:, ix);
y215 = pres(lev215); %r215 = prho(lev215, :); %r215 = r215(:, ix);
%
varint = 0:0.1:3; cval = [0.1, 0.5, 1, 1.5, 2, 2.5, 2.3, 2.4];
figure(3);clf
axes('position', [left, bottom(1), width, height(1)]);
contourf(x, y02, z02, varint, 'LineStyle','none');
hold on
for stn = 1:length(ctdHy19v03D)
    ra_plotsampleloc(ctdHy19v03D(stn).lat, ctdHy19v03D(stn).sdepth, 12, [0.6, 0.6, 0.6])
end
[c, h] = contour(x, y02, z02, cval, 'linest', '-', 'color', 'w', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'w', 'fontweigh', 'bold'); clear c h
% 
hold off
cmocean('matter', length(varint)-1); caxis([varint(1), varint(end)])
%
title('IN2019\_V03, 110^{\circ}E')
%
yst = 0; yint= 100; yen = 495;
xlim([-35, -10]); ylim([yst, yen])
axis tight
ra_subsurfembellish('y', yst, yint, yen, ytick, yticklabel, ylabelname)

% 500 to 5000 dbar
axes('position', [left, bottom(2), width, height(2)]);
contourf(x, y215, z215, varint, 'LineStyle','none');
hold on
for stn = 1:length(ctdHy19v03D)
    ra_plotsampleloc(ctdHy19v03D(stn).lat, ctdHy19v03D(stn).sdepth, 12, [0.6, 0.6, 0.6])
end
[c, h] = contour(x, y215, z215, cval, 'linest', '-', 'color', 'w', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'w', 'fontweigh', 'bold'); clear c h
hold off
cmocean('matter', length(varint)-1); caxis([varint(1), varint(end)])
%
ra_colorbarV2('vertical', 'normalized', [0.885, bottom(2)+0.02, 0.015, 0.8], 'Phosphate [{\mu}mol L^{-1}]', 'in')
yst = 500; yint= 500; yen = 5000;
xlim([-35, -10]); ylim([yst, yen])
axis tight
ra_subsurfembellish('xy', yst, yint, yen, ytick, yticklabel, ylabelname)
% jprint('./', 'IN2019V03D_phos_hydro', '-r300', '-painters')

%% Si* = Si-Nox

x = lat;
var = [ctdHy19v03Dlat.silicate] - [ctdHy19v03Dlat.nox];
%
clear z*
% surface to 500 m
lev = pres <= 500;
z02 = var(lev, :);
y02 = pres(lev); %r02 = prho(lev, :);
% 500m to 5000 m
lev215 = pres > 500 & pres <= 5000;
z215 = var(lev215, :); %z215 = z215(:, ix);
y215 = pres(lev215); %r215 = prho(lev215, :); %r215 = r215(:, ix);
%
varint = [-100:2.5:100]; cval = [-10, 0, 10, 30, 50, 80, 90];
figure(3);clf
axes('position', [left, bottom(1), width, height(1)]);
contourf(x, y02, z02, varint, 'LineStyle','none');
hold on
for stn = 1:length(ctdHy19v03D)
    ra_plotsampleloc(ctdHy19v03D(stn).lat, ctdHy19v03D(stn).sdepth, 12, [0.6, 0.6, 0.6])
end
[c, h] = contour(x, y02, z02, cval, 'linest', '-', 'color', 'w', 'linewi', 2);
clabel(c, h, 'fontsize', 16, 'color', 'w', 'fontweigh', 'bold'); clear c h
% 
hold off
cmocean('curl', length(varint)-1); caxis([varint(1), varint(end)])
%
% title('IN2019\_V03, 110^{\circ}E')
%
yst = 0; yint= 100; yen = 495;
xlim([-35, -10]); ylim([yst, yen])
axis tight
ra_subsurfembellish('y', yst, yint, yen, ytick, yticklabel, ylabelname)

% 500 to 5000 dbar
axes('position', [left, bottom(2), width, height(2)]);
contourf(x, y215, z215, varint, 'LineStyle','none');
hold on
for stn = 1:length(ctdHy19v03D)
    ra_plotsampleloc(ctdHy19v03D(stn).lat, ctdHy19v03D(stn).sdepth, 12, [0.6, 0.6, 0.6])
end
[c, h] = contour(x, y215, z215, cval, 'linest', '-', 'color', 'w', 'linewi', 2);
clabel(c, h, 'fontsize', 16, 'color', 'w', 'fontweigh', 'bold'); clear c h
hold off
cmocean('curl', length(varint)-1); caxis([varint(1), varint(end)])
%
ra_colorbarV2('vertical', 'normalized', [0.885, bottom(2)+0.02, 0.015, 0.8], 'Si^* (Silicate - Nox) [{\mu}mol L^{-1}]', 'in')
yst = 500; yint= 500; yen = 5000;
xlim([-35, -10]); ylim([yst, yen])
axis tight
ra_subsurfembellish('xy', yst, yint, yen, ytick, yticklabel, ylabelname)
 jprint('./', 'IN2019V03D_sistar_hydro', '-r300', '-painters')

%% N* = Nox-16Po4

x = lat;
var = [ctdHy19v03Dlat.nox] - 16*[ctdHy19v03Dlat.phos];
%
clear z*
% surface to 500 m
lev = pres <= 500;
z02 = var(lev, :);
y02 = pres(lev); %r02 = prho(lev, :);
% 500m to 5000 m
lev215 = pres > 500 & pres <= 5000;
z215 = var(lev215, :); %z215 = z215(:, ix);
y215 = pres(lev215); %r215 = prho(lev215, :); %r215 = r215(:, ix);
%
varint = [-7:0.2:7]; cval = [-5.5, -4, -2, -1];
figure(3);clf
axes('position', [left, bottom(1), width, height(1)]);
contourf(x, y02, z02, varint, 'LineStyle','none');
hold on
for stn = 1:length(ctdHy19v03D)
    ra_plotsampleloc(ctdHy19v03D(stn).lat, ctdHy19v03D(stn).sdepth, 12, [0.6, 0.6, 0.6])
end
[c, h] = contour(x, y02, z02, cval, 'linest', '-', 'color', 'w', 'linewi', 2);
clabel(c, h, 'fontsize', 16, 'color', 'w', 'fontweigh', 'bold'); clear c h
% 
hold off
cmocean('balance', length(varint)-1); caxis([varint(1), varint(end)])
%
% title('IN2019\_V03, 110^{\circ}E')
%
yst = 0; yint= 100; yen = 495;
xlim([-35, -10]); ylim([yst, yen])
axis tight
ra_subsurfembellish('y', yst, yint, yen, ytick, yticklabel, ylabelname)

% 500 to 5000 dbar
axes('position', [left, bottom(2), width, height(2)]);
contourf(x, y215, z215, varint, 'LineStyle','none');
hold on
for stn = 1:length(ctdHy19v03D)
    ra_plotsampleloc(ctdHy19v03D(stn).lat, ctdHy19v03D(stn).sdepth, 12, [0.6, 0.6, 0.6])
end
[c, h] = contour(x, y215, z215, cval, 'linest', '-', 'color', 'w', 'linewi', 2);
clabel(c, h, 'fontsize', 16, 'color', 'w', 'fontweigh', 'bold'); clear c h
hold off
cmocean('balance', length(varint)-1); caxis([varint(1), varint(end)])
%
ra_colorbarV2('vertical', 'normalized', [0.885, bottom(2)+0.02, 0.015, 0.8], 'N^* (Nox - 16xPO4) [{\mu}mol L^{-1}]', 'in')
yst = 500; yint= 500; yen = 5000;
xlim([-35, -10]); ylim([yst, yen])
axis tight
ra_subsurfembellish('xy', yst, yint, yen, ytick, yticklabel, ylabelname)
 jprint('./', 'IN2019V03D_Nstar_hydro', '-r300', '-painters')