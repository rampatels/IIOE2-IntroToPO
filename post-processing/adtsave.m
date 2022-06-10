% This script prepares absolute dynamic topography climatology for the
% Indian Ocean region.


% mean dynamic topography
[var, ~, ~] = nc2mat('./mdt9312APDRC.nc');
[varmsla, ~, ~] = nc2mat('./AVISOupdmclimDTMSLA9312.nc');

% Reading both data

lonmdt = var.LON321_521; latmdt = var.LAT117_317;
lonmsla = varmsla.LON241_391; latmsla = varmsla.LAT285_458;

mdt = var.MDOT; % 'mean dynamic ocean topography (cm)'
msla = varmsla.SLA;
msla = squeeze(nanmean(msla, 3));
%% Projecting data
[xq, yq] = meshgrid(lonmdt, latmdt);
[x, y] = meshgrid(lonmsla, latmsla);
%
msla2amdt = interp2(x, y, msla', xq, yq);

%%
figure(1);clf
subplot(1,2,1)
pcolor(lonmsla, latmsla, msla'); shading interp
title('original')

subplot(1,2,2)
pcolor(lonmdt, latmdt, msla2amdt); shading interp
title('projected')

%% 
adt = mdt + msla2amdt';

figure(2);clf
pcolor(lonmdt, latmdt, adt'); shading interp

%% 
adt = adt';
save adtmjwa lonmdt latmdt adt