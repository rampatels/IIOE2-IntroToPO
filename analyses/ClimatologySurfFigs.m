% This script plots climatological properties over the east indian ocean
% region for Indian Ocean voyage first paper. 
% I am focussing on only Salinity and Chlorophyll data given they are
% important for the paper.. So this script will prepare the final version of the
% Figures..

% Reading data
clear; clc; close all

% CARS climatology of salinity
load salcarsWA.mat
% monthly climatology
msalmj = squeeze(nanmean(msal(5:6, :, :, :), 1));
msalmj = squeeze(msalmj(1, :, :)); % only surface value

% computing physical properties
sa = gsw_SA_from_SP(msalmj, depthsal(1), lonsal, latsal); % Absolute Salinity g/kg

figure();clf
subplot(211)
contourf(lonsal, latsal, msalmj', 32:0.1:36, 'linest', 'none')
cmocean haline
caxis([32, 36])

subplot(212)
contourf(lonsal, latsal, sa', 32:0.1:36, 'linest', 'none')
cmocean haline
caxis([32, 36])

% Loading mean dynamic topography
load adtmjwa

figure();clf
contourf(lonmdt, latmdt, adt)

% Chlorophyll-a from Modis-aqua
load chlamodisaquaWA.mat

figure();clf
contourf(lonmodis, latmodis, log10(chlamj)')
cmocean algae

% Loading stations 
load IN19V03HyS.mat
load DM63V02Hy.mat

%% Surface Salinity and SSH contours for Paper...
close all
% FX = 16; FY = 13.7; % page size A4
% FIG1 = figure('units', 'centimeters', 'Position', [0 0 FX FY]);
% set(FIG1, 'PaperUnits', 'centimeters')
% set(FIG1, 'Papersize', [FX FY])
% set(FIG1, 'PaperPositionMode', 'manual');
% set(FIG1, 'PaperPosition', [0 0 FX FY]);
% width = 14; height = 12;
% left = [1.5, 1.5+width+0.1];
% bottom = [1.2, 2];

%
m_proj('mercator', 'lat', [-40, -5], 'lon', [100, 130]);
set(gcf, 'color', 'w');  % otherwise 'print' turns lakes black
% salinity and vectors
salInt = 32:0.2:36;
%
% axes('units', 'centimeters', 'position', [left(1) bottom(1) width height]);
m_contourf(lonsal, latsal, sa', salInt, 'linest', 'none')
cmocean('haline', length(salInt) -1)
caxis([salInt(1), salInt(end)])
% ra_colorbarV2('vertical', 'centimeters', [left(2) bottom(1) 0.45 height], 'Salinity [PSU]', 'in')
h = colorbar;
ylabel(h, 'Absolute Salinity [g/kg]', 'fontsize', 16, 'fontweight', 'normal');
set(h, 'fontsize', 16, 'fontweight', 'normal'); clear h
hold on
% coastline
m_usercoast('WAcoastf', 'patch', [0.7, 0.7, 0.7], 'linewidth', 2, 'edgecolor', [0, 0, 0]);
% mdt
[c, h] = m_contour(lonmdt, latmdt, adt, 50:5:120, 'linest', '-', 'color', [1, 1, 1], 'linewi', 1.5);
v = 50:10:120; %[60, 80, 100, 120];
clabel(c, h, v, 'fontsize', 12, 'color', [1, 1, 1], 'fontweigh', 'bold'); clear c h
% stations
m_plot([ctdHy63v02D500.lon], [ctdHy63v02D500.lat], '*m', 'Markersize', 9, 'linewidth', 1.2)
m_plot([ctdHy19v03S.lon], [ctdHy19v03S.lat], 'ok', 'MarkerSize', 6, 'linewi', 1.5)
hold off
m_grid('box','fancy','fontsize',16,'fontweigh','bold','linest','none','tickdir','out', 'xtick', 100:10:130, 'ytick', -40:5:-5);

% print(gcf, '-dpng', '-r300', '-painters', 'ClimatologySSSmdtCTD')
%jprint('./', 'ClimatologySSSmdtCTD','-dpng', '-r300','-painters')
%% chlorophyll climatology
figure();clf
m_proj('mercator', 'lat', [-40, -5], 'lon', [100, 130]);
set(gcf, 'color', 'w');  % otherwise 'print' turns lakes black
chlInt = 0:0.025:0.8;
m_contourf(lonmodis, latmodis, chlamj', chlInt, 'linest', 'none')
% cmocean('algae', length(chlInt) -1)
colormap(m_colmap('jet', length(chlInt) - 1))
caxis([chlInt(1), chlInt(end)])
h = colorbar;
ylabel(h, 'Chlorophyll-a [mg/m^{3}]', 'fontsize', 16, 'fontweight', 'normal');
set(h, 'fontsize', 14, 'fontweight', 'normal')
% yticks = chlInt(1):0.2:chlInt(end);
% set(h, 'ytick', yticks , 'yticklabel', round(10.^yticks, 1))
hold on
% [c, h] = m_contour(lonmdt, latmdt, adt, 40:10:120, 'linest', '-', 'color', [0.6, 0.6, 0.6], 'linewi', 1.2);
% v = [60, 80, 100, 120];
% clabel(c, h, v, 'fontsize', 12, 'color', [0.6, 0.6, 0.6], 'fontweigh', 'bold'); clear c h
m_plot([ctdHy63v02D500.lon], [ctdHy63v02D500.lat], '*m', 'Markersize', 9, 'linewidth', 1.2)
m_plot([ctdHy19v03S.lon], [ctdHy19v03S.lat], 'ok', 'MarkerSize', 6, 'linewi', 1.5)
m_text([ctdHy19v03S.lon]-2.4, [ctdHy19v03S.lat], num2str([1:20]'), 'fontsize', 14)
hold off
m_usercoast('WAcoastf', 'patch', [0.7, 0.7, 0.7], 'linewidth', 2, 'edgecolor', [0, 0, 0]);
m_grid('box','fancy','fontsize',14,'fontweigh','bold','linest','none','tickdir','out', 'xtick', 100:10:130, 'ytick', -40:5:-5);

% print(gcf, '-dpng', '-r300', '-painters', 'ClimatologyChlaCTD')
%jprint('./', 'ClimatologyChlaCTD','-dpng', '-r300','-painters')
%% 1) Removing small and very large velocity 
% I will keep velocities from 1 cm/s to 100 cm/s 

mag = sqrt(umoscar.^2 + vmoscar.^2);
keep = mag >= 2 & mag <= 90;
u = umoscar; v = vmoscar;
u(~keep) = NaN; v(~keep) = NaN;

% %%
% % cint = [-350, -200:25:0, 50, 100]; scale = 0.01;
% % %
% % [cs, ch] = m_contourf(lon1, lat1, mthf', cint, 'linest', 'none');
% % cmocean('-tarn', 'pivot', 0); caxis([cint(1), cint(end)])
% % ax = m_contfbar([0.1, 0.55], 0.15,  cs, ch, 'edgecolor', 'none', 'endpiece', 'yes' ,'axfrac', .03, 'levels', 'match'); %0.85, [0.2, 0.8]
% % title(ax, 'MTHF [MW m^{-1}]', 'fontsize', 14, 'fontweigh', 'bold'); 
% % set(ax, 'fontsize', 12, 'fontweigh', 'bold')


% Figures for spatial distribution at surface 
% close all
m_proj('mercator', 'lat', [-40, -5], 'lon', [100, 130]);
sc = 0.2;

figure(1);clf
set(gcf, 'color', 'w');  % otherwise 'print' turns lakes black
% salinity and vectors
salInt = 32:0.1:36;
subplot_tight(2, 1, 1, [0.06, 0.2])
m_contourf(lonsal, latsal, squeeze(msalmj(1, :, :))', salInt, 'linest', 'none')
cmocean('haline', length(salInt) -1)
caxis([salInt(1), salInt(end)])
colorbar
h = colorbar;
ylabel(h, 'Salinity [PSU]', 'fontsize', 14, 'fontweight', 'bold');
set(h, 'fontsize', 14, 'fontweight', 'bold')
hold on
m_usercoast('WAcoastf', 'patch', [1, 1, 1], 'linewidth', 2, 'edgecolor', [0, 0, 0]);
m_quiver(x(1:2:end, 1:2:end), y(1:2:end, 1:2:end), u(1:2:end, 1:2:end)'.*sc, v(1:2:end, 1:2:end)'.*sc,  0, 'color', [0.4, 0.4, 0.4],'linewidth',0.75,'MaxHeadSize',0.1)
m_text(114.15, -23.5, '10 cm s^{-1}', 'fontsize', 12, 'fontweigh', 'bold', 'color',[0.4, 0.4, 0.4])
m_quiver(115.2, -22.5, 10.*sc, 0, 0, 'color', [0.4, 0.4, 0.4],'MaxHeadSize',1,'linewidth', 2)
m_plot([ctdHy63v02D500.lon], [ctdHy63v02D500.lat], '.r', 'Markersize', 15)
m_plot([ctdHy19v03S.lon], [ctdHy19v03S.lat], 'ow', 'MarkerSize', 10, 'linewi', 1.5)
hold off
m_grid('box','fancy','fontsize',14,'fontweigh','bold','linest','none','tickdir','out', 'xticklabel', [], 'xtick', 100:5:130, 'ytick', -40:5:-5);




%% Geostrophic velocity overlaid plot that Helen gave me..
clear;clc

% CARS climatology of salinity
load salcarsWA.mat
% monthly climatology
msalmj = squeeze(nanmean(msal(5:6, :, :, :), 1));

% Chlorophyll-a from Modis-aqua
load chlamodisaquaWA.mat

% mean dynamic topography
load adtmjwa

% New velocity that Helen gave
var = nc2mat('cars-ug-vg-RL2000-adjusted.nc');
%
%% Reading velocity data
[X, Y] = meshgrid(var.XAXM, var.YAXM);

u = var.UG; v = var.VG;
U = squeeze(nanmean(u(:, :, 1, 5:6), 4))'; V = squeeze(nanmean(v(:, :, 1, 5:6), 4))';
U = U(:, 1:240)*100; V = V(1:280, :)*100; 
%
m_proj('mercator', 'lat', [-40, -10], 'lon', [90, 120]);
sc = 0.05;
figure(1);clf
m_quiver(X(1:3:end, 1:3:end), Y(1:3:end, 1:3:end), U(1:3:end, 1:3:end).*sc, V(1:3:end, 1:3:end).*sc,  0, 'color', 'k','linewidth',0.75,'MaxHeadSize',0.1) % ,'ShowArrowHead','off');
hold off
m_usercoast('WAcoastf', 'patch', [0.7, 0.7, 0.7], 'linewidth', 2.5, 'edgecolor', [0, 0, 0]);
m_grid('box','fancy','fontsize',14,'fontweigh','bold','linest','none','tickdir','out');
% 
%% 1) Removing small and very large velocity 
% I will keep velocities from 1 cm/s to 100 cm/s 

mag = sqrt(U.^2 + V.^2);
keep = mag >= 1 & mag <= 40;
u = U; v = V;
u(~keep) = NaN; v(~keep) = NaN;

% Figures for spatial distribution at surface 
% close all
m_proj('mercator', 'lat', [-40, -10], 'lon', [90, 120]);
sc = 0.08;

figure(1);clf
set(gcf, 'color', 'w');  % otherwise 'print' turns lakes black
% salinity and vectors
salInt = 33.6:0.1:35.9;
subplot_tight(2, 1, 1, [0.01, 0.15])
m_contourf(lonsal, latsal, squeeze(msalmj(1, :, :))', salInt, 'linest', 'none')
cmocean('haline', length(salInt) -1)
caxis([salInt(1), salInt(end)])
colorbar
h = colorbar;
ylabel(h, 'Salinity [PSU]', 'fontsize', 14, 'fontweight', 'bold');
set(h, 'fontsize', 14, 'fontweight', 'bold')
hold on
m_usercoast('WAcoastf', 'patch', [1, 1, 1], 'linewidth', 2, 'edgecolor', [0, 0, 0]);
m_quiver(X(1:3:end, 1:3:end), Y(1:3:end, 1:3:end), u(1:3:end, 1:3:end).*sc, v(1:3:end, 1:3:end).*sc,  0,  'color', [0.4, 0.4, 0.4],'linewidth',0.75,'MaxHeadSize',0.1)
m_text(114.15, -23.5, '10 cm s^{-1}', 'fontsize', 12, 'fontweigh', 'bold', 'color',[0.4, 0.4, 0.4])
m_quiver(115.2, -22.5, 10.*sc, 0, 0, 'color', [0.4, 0.4, 0.4],'MaxHeadSize',1,'linewidth', 2)
hold off
m_grid('box','fancy','fontsize',14,'fontweigh','bold','linest','none','tickdir','out', 'xticklabel', [], 'xtick', 90:5:120, 'ytick', -40:5:-10);

% chlorophyll and mdt
chlInt = -0.1:0.05:0.8;
subplot_tight(2, 1, 2, [0.01, 0.15])
m_contourf(lonmodis, latmodis, chlamj', chlInt, 'linest', 'none')
cmocean('algae', length(chlInt) -1)
caxis([chlInt(1), chlInt(end)])
colorbar
h = colorbar;
ylabel(h, 'Chlorophyll-a [mg/m^{3}]', 'fontsize', 14, 'fontweight', 'bold');
set(h, 'fontsize', 14, 'fontweight', 'bold')
yticks = chlInt(1):0.2:chlInt(end);
set(h, 'ytick', yticks , 'yticklabel', round(10.^yticks, 1))
hold on
[c, h] = m_contour(lonmdt, latmdt, adt, 40:10:120, 'linest', '-', 'color', [0.5, 0.5, 0.5], 'linewi', 1.2);
v = [60, 80, 100, 120];
clabel(c, h, v, 'fontsize', 12, 'color', [0.5, 0.5, 0.5], 'fontweigh', 'bold'); clear c h
hold off
m_usercoast('WAcoastf', 'patch', [1, 1, 1], 'linewidth', 1.2, 'edgecolor', [0, 0, 0]);
m_grid('box','fancy','fontsize',14,'fontweigh','bold','linest','none','tickdir','out','xtick', 90:5:120, 'ytick', -40:5:-10);

% print(gcf, '-dpng', '-r200', '-painters', 'ClimatologySurface240geoCTD')