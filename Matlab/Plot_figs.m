% Plot_figs.m
% This script plots results reported in Cacho and Hester (2022a and 2022b)
%
%% =================================== Plot phase diagram in Figure 6  
clear;
path(path, '.\Data\');
load('vic_2010.mat'); % load file to plot
Wmean = permute(sum(ev.W .* garea)./sum(garea),[3,2,1]);
Bmean = permute(sum(ev.B .* garea)./sum(garea),[3,2,1]);
figure(1);
plot(Bmean', Wmean','.-');
axis([0,7,0,1.1]);
xlabel('Biocontrol density (adults/km^2)');
ylabel('Wasp density (nests/km^2)');
title('VIC ({\itx_c}=20, {\it x_p}=10)')
grid
%% ================================== Plot wasp density maps in Figure 7
clear;
path(path, '.\Data\');
load('vic_2010.mat'); % load file to plot
load('W0_prm_vic.mat', 'gpoly'); % load map polyshape
Wmap = permute(mean(ev.W,2), [1,3,2]);
wmax = 1; % max value for mapping to colormap 
figure(2)
subplot(1,2,1);
  Plot_map(gpoly,Wmap(:,1),wmax);
  title('VIC treatment 1');
  xlabel('Longitude');
  ylabel('Latitude');
  grid;
subplot(1,2,2);
  Plot_map(gpoly,Wmap(:,9),wmax);
  title('VIC treatment 9');
  xlabel('Longitude');
  ylabel('Latitude');
  grid;
%% ======================== Plot wasp GA parameter comparison in Figure 8
clear;
path(path, '.\Data\');
%
load base_data_vic;
wp_mat_vic = wp_mat;
load base_data_nsw;
wp_mat_nsw = wp_mat;
%
figure(3)
subplot(2,2,1);
    ksdensity(wp_mat_vic(:,1))
    hold on
    ksdensity(wp_mat_nsw(:,1))
    grid
    xlabel('x');
    ylabel('f(x)');
    title('A (\alpha_W)')
    hold off
subplot(2,2,2);
    ksdensity(wp_mat_vic(:,2))
    hold on
    ksdensity(wp_mat_nsw(:,2))
    grid
    xlabel('x');
    ylabel('f(x)');
    title('B (\kappa_W)')
    hold off
subplot(2,2,3);
    ksdensity(wp_mat_vic(:,3), 'Support', 'positive')
    hold on
    ksdensity(wp_mat_nsw(:,3))
    grid
    xlabel('x');
    ylabel('f(x)');
    title('C (\gamma_W)')
    legend('VIC','NSW');
    hold off
subplot(2,2,4);
    ksdensity(wp_mat_vic(:,4))
    hold on
    ksdensity(wp_mat_nsw(:,4))
    grid
    xlabel('x');
    ylabel('f(x)');
    title('D (\delta_W)')
    hold off
%
%% =================== Plot wasp density vs household density
clear;
path(path, '.\Data\');
load('vic_2010.mat'); % load file to plot
Wmap = permute(mean(ev.W,2), [1,3,2]);
figure(4)
  plot(log(hh_dens),(Wmap(:,1)), '.');
  hold on;
  plot(log(hh_dens),(Wmap(:,9)), '.');
  title('VIC');
  xlabel('Household density per km^2 (log_e)');
  ylabel('Mean wasp density (nests/km^2');
  legend('trt 1', 'trt 9');
  grid;
  hold off
%% =================================================== embedded functions
% --------------------------------------------------- Plot_map funtion
function [pg] = Plot_map(m, x, x_max)
    % Plots the values of x on map m using x_max as the upper bound for the 
    %   colormap, where m is the map expressed as a polyshape type.  
    %
    pg = plot(m);
    cmap = colormap;
    nmap = size(cmap,1);
    ip = (nmap-1) * (x./x_max); 
    ip = round(ip) + 1;
    for i = 1 : length(pg)
        pg(i).FaceColor = cmap(ip(i),:);
        pg(i).FaceAlpha = 0.8;
    end
    colorbar;
end % Plot_map function

