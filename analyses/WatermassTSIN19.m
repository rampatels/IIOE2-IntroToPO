% This script prepares Theta-SA plot coloured by various variables to
% identify different watermasses along 110E line for the entire column
% stations only...

% Reading data
clear; clc; close all
%
load IN19V03D.mat
%
lon = [ctd19v03D.lon];
lat = [ctd19v03D.lat];
stn = [ctd19v03D.Deployment];
% Computing depths
pres = [ctd19v03D.upres]; % as all the data are projected on uniform pressure
z = gsw_z_from_p(pres, lat);
depth = gsw_depth_from_z(z); clear z

pres = pres(:, 1);
% For now I am focussing on first sensor only.. 
% computing physical properties
sa = gsw_SA_from_SP([ctd19v03D.SP], pres, lon, lat); % Absolute Salinity g/kg
ct = gsw_CT_from_t(sa, [ctd19v03D.t], pres); % Conservative Temperature
prho = gsw_rho(sa, ct, 0) - 1000; % Potential Density 10.1325 dbar considered as reference level

% 
oxy = [ctd19v03D.oxy];

%% Potential vorticity computing from the CTD station
% Theoratically, it is vertical density gradient times f/rho. 
% For the system where realtive worticity is negligible.
%
Prho = prho + 1000;
[PV, Z] = deal(NaN(size(Prho)));
%
figure(5);clf
for ii =  1:size(Prho, 2)
    [pv, z] = ra_pVorticityCTD(Prho(:, ii), depth(:, ii), lat(ii));
    pv = pv * 10^9;
    pv = smooth(pv, 100); % 100 data points
    PV(:, ii) = pv; % Multiplication to remove power so when you write axis use x10^-9
    Z(:, ii) = z;
    plot(pv, z); hold on
    clear pv z
end
clear ii
hold off
axis ij

%% plotting CT SA diagram colored by latitudes and then by oxygen 
% Potential Density 10.1325 dbar considered as reference level
smin = 34; smax = 36.3;
tmin = 0; tmax = 30;
ss = smin:0.1:smax;
tt = tmin:0.1:tmax;
[sgrid, tgrid] = meshgrid(ss, tt);
pden = gsw_rho(sgrid, tgrid, 0) - 1000;
clear tt ss
% color to choose from
eddiescol = cbrewer('qual', 'Paired', 10);

% Setting up axis
left = 0.1; bottom = 0.2;
width = 0.75; height = 0.7;
%% 1) each station is coloured by latitudes to identify their position

xtick = 110; ytick = -35:5:-10;
[~, yticklabel] = XYTickLabel(xtick, ytick);

cpmin = min(lat); cprng = max(lat) - min(lat); msz = 1; % marker size
figure(3);clf
set(gcf, 'color', 'w');
axes('position',[left, bottom, width, height])
[c, h] = contour(sgrid, tgrid, pden, 20:0.25:30, 'color', [.7 .7 .7]);
clabel(c , h, 22:30, 'fontsize', 14); hold on
for ii = 1:length(stn)
    colourplot(sa(:, ii), ct(:, ii), lat(ii).*ones(size(sa(:, ii))), '.-', cpmin, cprng, msz, length(lat)-1)
    hold on
end
% [c, h] = contour(sgrid, tgrid, pden, 20:0.25:30, 'color', [.7 .7 .7]);
% clabel(c , h, 22:30)
% contour(sgrid, tgrid, pden, [26, 26], 'color', eddiescol(2, :), 'linest', '--', 'linewi', 2.5); % STUW
% text(34.1, 12, 'STUW', 'fontsize', 12, 'fontweigh', 'bold', 'color', eddiescol(2, :))
% contour(sgrid, tgrid, pden, [26.8, 26.9], 'color', eddiescol(1, :), 'linest', '-.', 'linewi', 2.5); % SAMW
% text(34.1, 7, 'SAMW', 'fontsize', 12, 'fontweigh', 'bold', 'color', eddiescol(1, :))
% contour(sgrid, tgrid, pden, [27, 27.3], 'color', eddiescol(8, :), 'linest', ':', 'linewi', 2.5); % AAIW
% text(34.1, 2, 'AAIW', 'fontsize', 12, 'fontweigh', 'bold', 'color', eddiescol(8, :))

% % Water type
% plot([33.8, 33.8], [-1.85, 5],'.','MarkerSize',20, 'color', 'k')
% text(33.82, 5, 'Antarctic Surface', 'fontsize', 12, 'fontweigh', 'bold')
% plot(34.35, 4, '.','MarkerSize',20, 'color', 'k')
% text(34.37, 4, 'Modified Intermediate water', 'fontsize', 12, 'fontweigh', 'bold')
% plot(34.56, 8.8,'.k','MarkerSize',20, 'color', 'k')
% text(34.58, 8.8, 'Mode Water', 'fontsize', 12, 'fontweigh', 'bold')
% plot(34.62, 2.5,'.k','MarkerSize',20, 'color', 'k')
% text(34.64, 2.5, 'UCDW', 'fontsize', 12, 'fontweigh', 'bold')
hold off
% hold off
axis([smin smax tmin tmax])
h = colorbar('vertical','position',[0.86 0.2 0.025 0.7]);
set(h,'fontsize',16,'fontweigh','bold')
set(h,'ytick', ytick, 'yticklabel', yticklabel)
% ylabel(h,'Latitudes [\circ]','fontsize',14,'fontweigh','bold')
xlabel('Absolute Salinity [g/kg]','fontsize',16,'fontweight','bold')
ylabel('Conservative Temperature [^oC]','fontsize',16,'fontweight','bold')
set(gca,'box','on','linewidth',2,'fontsize',16,'fontweight','bold','tickdir','in')

% print(gcf, '-dpng', '-r200', '-painters', 'CTvsSAcbyLat19V03sensor1')
% jprint('./', 'CTvsSAcbyLat19V03sensor1','-dpng', '-r300','-painters')
%% 2) each station is coloured by dissolved oxygen to identify old waters

cpmin = min(oxy(:)); cpmax = max(oxy(:));
figure(3);clf
axes('position',[left, bottom, width, height])
for ii = 1:length(stn)
    scatter(sa(:, ii), ct(:, ii), 8, oxy(:, ii), 'o', 'filled')
    hold on
end
cmocean('oxy')
caxis([round(cpmin), round(cpmax)])
[c, h] = contour(sgrid, tgrid, pden, 20:0.25:30, 'color', [.7 .7 .7]);
clabel(c , h, 22:30, 'fontsize', 14)
% contour(sgrid, tgrid, pden, [26, 26], 'color', eddiescol(2, :), 'linest', '--', 'linewi', 2.5); % STUW
% text(34.1, 12, 'STUW', 'fontsize', 12, 'fontweigh', 'bold', 'color', eddiescol(2, :))
% contour(sgrid, tgrid, pden, [26.8, 26.9], 'color', eddiescol(1, :), 'linest', '-.', 'linewi', 2.5); % SAMW
% text(34.1, 7, 'SAMW', 'fontsize', 12, 'fontweigh', 'bold', 'color', eddiescol(1, :))
% contour(sgrid, tgrid, pden, [27, 27.3], 'color', eddiescol(8, :), 'linest', ':', 'linewi', 2.5); % AAIW
% text(34.1, 2, 'AAIW', 'fontsize', 12, 'fontweigh', 'bold', 'color', eddiescol(8, :))
% Water type
% plot([33.8, 33.8], [-1.85, 5],'.','MarkerSize',20, 'color', 'k')
% text(33.82, 5, 'Antarctic Surface', 'fontsize', 12, 'fontweigh', 'bold')
% plot(34.35, 4, '.','MarkerSize',20, 'color', 'k')
% text(34.37, 4, 'Modified Intermediate water', 'fontsize', 12, 'fontweigh', 'bold')
% plot(34.56, 8.8,'.k','MarkerSize',20, 'color', 'k')
% text(34.58, 8.8, 'Mode Water', 'fontsize', 12, 'fontweigh', 'bold')
% plot(34.62, 2.5,'.k','MarkerSize',20, 'color', 'k')
% text(34.64, 2.5, 'UCDW', 'fontsize', 12, 'fontweigh', 'bold')
% @ 14 august for the review added reference profiles
% plot(Refs, Reft, 'linewidth', 2, 'color', ra_color('brown')) % obsRef
% plot(aRefs, aReft ,'color', ra_color('maroon'), 'linewidth', 2) % annual climatology
% plot(mRefs(3,:), mReft(3,:), 'color',[ 1, 0.3, 1],'linewidth', 2) % mar-climatology
% plot(mRefs(4,:), mReft(4,:), 'color', 'b','linewidth', 2) % april-climatology
hold off
axis([smin smax tmin tmax])
h = colorbar('vertical','position',[0.86 0.2 0.025 0.7]);
set(h,'fontsize',16,'fontweigh','bold')
ylabel(h,'Oxygen [{\mu}mol L^{-1}]','fontsize',16,'fontweigh','bold')
% xlabel('Absolute Salinity [g/kg]')
% ylabel('Conservative Temperature [^oC]')
% xlabel('Salinity [psu]','fontsize',12,'fontweight','bold')
% ylabel(' \Theta [^oC]','fontsize',12,'fontweight','bold')
xlabel('Absolute Salinity [g/kg]','fontsize',16,'fontweight','bold')
ylabel('Conservative Temperature [^oC]','fontsize',16,'fontweight','bold')
set(gca,'box','on','linewidth',2,'fontsize',16,'fontweight','bold','tickdir','in')
% print(gcf,'-dpng','-r300','-painters','CTvsSAcbyOxy19V03sensor1')
% jprint('./', 'CTvsSAcbyOxy19V03sensor1','-dpng', '-r300','-painters')
%% 3) each station is coloured by depth

cpmin = min(depth(:)); cpmax = max(depth(:));
figure(3);clf
axes('position',[left, bottom, width, height])
for ii = 1:length(stn)
    scatter(sa(:, ii), ct(:, ii), 8, depth(:, ii), 'o', 'filled')
    hold on
end
cmocean('deep', 50)
% colormap(parula(6))
caxis([round(cpmin), 2000])
[c, h] = contour(sgrid, tgrid, pden, 20:0.25:30, 'color', [.7 .7 .7]);
clabel(c , h, 22:30)
contour(sgrid, tgrid, pden, [26, 26], 'color', eddiescol(2, :), 'linest', '--', 'linewi', 2.5); % STUW
text(34.1, 12, 'STUW', 'fontsize', 12, 'fontweigh', 'bold', 'color', eddiescol(2, :))
contour(sgrid, tgrid, pden, [26.8, 26.9], 'color', eddiescol(1, :), 'linest', '-.', 'linewi', 2.5); % SAMW
text(34.1, 7, 'SAMW', 'fontsize', 12, 'fontweigh', 'bold', 'color', eddiescol(1, :))
contour(sgrid, tgrid, pden, [27, 27.3], 'color', eddiescol(8, :), 'linest', ':', 'linewi', 2.5); % AAIW
text(34.1, 2, 'AAIW', 'fontsize', 12, 'fontweigh', 'bold', 'color', eddiescol(8, :))
% Water type
% plot([33.8, 33.8], [-1.85, 5],'.','MarkerSize',20, 'color', 'k')
% text(33.82, 5, 'Antarctic Surface', 'fontsize', 12, 'fontweigh', 'bold')
% plot(34.35, 4, '.','MarkerSize',20, 'color', 'k')
% text(34.37, 4, 'Modified Intermediate water', 'fontsize', 12, 'fontweigh', 'bold')
% plot(34.56, 8.8,'.k','MarkerSize',20, 'color', 'k')
% text(34.58, 8.8, 'Mode Water', 'fontsize', 12, 'fontweigh', 'bold')
% plot(34.62, 2.5,'.k','MarkerSize',20, 'color', 'k')
% text(34.64, 2.5, 'UCDW', 'fontsize', 12, 'fontweigh', 'bold')
% @ 14 august for the review added reference profiles
% plot(Refs, Reft, 'linewidth', 2, 'color', ra_color('brown')) % obsRef
% plot(aRefs, aReft ,'color', ra_color('maroon'), 'linewidth', 2) % annual climatology
% plot(mRefs(3,:), mReft(3,:), 'color',[ 1, 0.3, 1],'linewidth', 2) % mar-climatology
% plot(mRefs(4,:), mReft(4,:), 'color', 'b','linewidth', 2) % april-climatology
hold off
axis([smin smax tmin tmax])
h = colorbar('vertical','position',[0.87 0.2 0.025 0.7]);
set(h,'fontsize',11,'fontweigh','bold')
ylabel(h,'Depth [m]','fontsize',14,'fontweigh','bold')
% xlabel('Absolute Salinity [g/kg]')
% ylabel('Conservative Temperature [^oC]')
% xlabel('Salinity [psu]','fontsize',12,'fontweight','bold')
% ylabel(' \Theta [^oC]','fontsize',12,'fontweight','bold')
xlabel('Absolute Salinity [g/kg]','fontsize',14,'fontweight','bold')
ylabel('Conservative Temperature [^oC]','fontsize',14,'fontweight','bold')
set(gca,'box','on','linewidth',2,'fontsize',14,'fontweight','bold','tickdir','in')
% print(gcf,'-dpng','-r300','-painters','CTvsSAcbydepth19V03sensor1')


%% 4) each station is coloured by Potential vorticity to identify SAMW

cpmin = min(PV(:)); cpmax = max(PV(:));
figure(3);clf
axes('position',[left, bottom, width, height])
for ii = 1:length(stn)
    scatter(sa(:, ii), ct(:, ii), 8, PV(:, ii), 'o', 'filled')
    hold on
end
% cmocean('-oxy', 25)
colormap(parula(25))
caxis([-0.1, round(cpmax)])
[c, h] = contour(sgrid, tgrid, pden, 20:0.25:30, 'color', [.7 .7 .7]);
clabel(c , h, 22:30)
contour(sgrid, tgrid, pden, [26, 26], 'color', eddiescol(2, :), 'linest', '--', 'linewi', 2.5); % STUW
text(34.1, 12, 'STUW', 'fontsize', 12, 'fontweigh', 'bold', 'color', eddiescol(2, :))
contour(sgrid, tgrid, pden, [26.8, 26.9], 'color', eddiescol(1, :), 'linest', '-.', 'linewi', 2.5); % SAMW
text(34.1, 7, 'SAMW', 'fontsize', 12, 'fontweigh', 'bold', 'color', eddiescol(1, :))
contour(sgrid, tgrid, pden, [27, 27.3], 'color', eddiescol(8, :), 'linest', ':', 'linewi', 2.5); % AAIW
text(34.1, 2, 'AAIW', 'fontsize', 12, 'fontweigh', 'bold', 'color', eddiescol(8, :))
% Water type
% plot([33.8, 33.8], [-1.85, 5],'.','MarkerSize',20, 'color', 'k')
% text(33.82, 5, 'Antarctic Surface', 'fontsize', 12, 'fontweigh', 'bold')
% plot(34.35, 4, '.','MarkerSize',20, 'color', 'k')
% text(34.37, 4, 'Modified Intermediate water', 'fontsize', 12, 'fontweigh', 'bold')
% plot(34.56, 8.8,'.k','MarkerSize',20, 'color', 'k')
% text(34.58, 8.8, 'Mode Water', 'fontsize', 12, 'fontweigh', 'bold')
% plot(34.62, 2.5,'.k','MarkerSize',20, 'color', 'k')
% text(34.64, 2.5, 'UCDW', 'fontsize', 12, 'fontweigh', 'bold')
% @ 14 august for the review added reference profiles
% plot(Refs, Reft, 'linewidth', 2, 'color', ra_color('brown')) % obsRef
% plot(aRefs, aReft ,'color', ra_color('maroon'), 'linewidth', 2) % annual climatology
% plot(mRefs(3,:), mReft(3,:), 'color',[ 1, 0.3, 1],'linewidth', 2) % mar-climatology
% plot(mRefs(4,:), mReft(4,:), 'color', 'b','linewidth', 2) % april-climatology
hold off
axis([smin smax tmin tmax])
h = colorbar('vertical','position',[0.87 0.2 0.025 0.7]);
set(h,'fontsize',11,'fontweigh','bold')
ylabel(h,'Potential Vorticity [ x 10^{-9} m^{-1}s^{-1}]','fontsize',14,'fontweigh','bold')
% xlabel('Absolute Salinity [g/kg]')
% ylabel('Conservative Temperature [^oC]')
% xlabel('Salinity [psu]','fontsize',12,'fontweight','bold')
% ylabel(' \Theta [^oC]','fontsize',12,'fontweight','bold')
xlabel('Absolute Salinity [g/kg]','fontsize',14,'fontweight','bold')
ylabel('Conservative Temperature [^oC]','fontsize',14,'fontweight','bold')
set(gca,'box','on','linewidth',2,'fontsize',14,'fontweight','bold','tickdir','in')
% print(gcf,'-dpng','-r300','-painters','CTvsSAcbyPV19V03sensor1')
