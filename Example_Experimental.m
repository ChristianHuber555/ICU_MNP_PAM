% clear all;
addpath("data")
addpath("functions")

load("RF_MNP_3mgml_MI_1.mat"); % directly received by Verasonics Vantage 64 with L11-5v 
fs_meas = 31.25e6; % System sampling rate
dataRF = double(dataRF(:,1+32:128-32,:)); % only middle 64 elements

t_meas = 0:1/fs_meas:size(dataRF,1)/fs_meas - 1/fs_meas;

param.xarray = (-31.5:1:31.5)*0.3e-3;
param.c = 1540;
%% interpolate to 55.11 MHz

param.fs = 55.11e6;
param.t = t_meas(1):1/param.fs:t_meas(end);
RF = zeros(size(param.t,2),size(dataRF,2),size(dataRF,3));
for i = 1:10
    RF(:,:,i) = interp1(t_meas,squeeze(dataRF(:,:,i)),param.t,'spline');
end

%% Cavitation map geometry
param.xdim = 101;
param.zdim = 201;
start_x = -5e-3;
end_x = 5e-3;
start_z = 35e-3;
end_z = 45e-3;
x = linspace(start_x,end_x,param.xdim);
z = linspace(start_z,end_z,param.zdim);

%%
tic
Cav_Map_DAS = DAS(RF, param,x ,z);
time_DAS = toc/10;
tic;
Cav_Map_DMAS = DMAS(RF, param,x ,z);
time_DMAS = toc/10;
tic
Cav_Map_RCB = RCB(RF, param,x ,z,0.1);
time_RCB = toc/10;

%% averaging and normalization in dB
Cav_Map_DAS_avg = squeeze(mean(Cav_Map_DAS,1));
Cav_Map_DAS_avg = Cav_Map_DAS_avg./max(Cav_Map_DAS_avg(:));
[CR_DAS,in_DAS,out_DAS] = compute_CR(Cav_Map_DAS_avg,x,z);
Cav_Map_DAS_avg = 10*log10(Cav_Map_DAS_avg);

Cav_Map_DMAS_avg = squeeze(mean(Cav_Map_DMAS,1));
Cav_Map_DMAS_avg = Cav_Map_DMAS_avg./max(Cav_Map_DMAS_avg(:));
[CR_DMAS,in_DMAS,out_DMAS] = compute_CR(Cav_Map_DMAS_avg,x,z);
Cav_Map_DMAS_avg = 10*log10(Cav_Map_DMAS_avg);

Cav_Map_RCB_avg = squeeze(mean(Cav_Map_RCB,1));
Cav_Map_RCB_avg = Cav_Map_RCB_avg./max(Cav_Map_RCB_avg(:));
[CR_RCB, in_RCB, out_RCB] = compute_CR(Cav_Map_RCB_avg,x,z);
Cav_Map_RCB_avg = 10*log10(Cav_Map_RCB_avg);
%% Imaging for -10 dB to 0 dB
clim = [-10, 0];

subplot(1,3,1)
imagesc(Cav_Map_DAS_avg, clim)
subplot(1,3,2)
imagesc(Cav_Map_DMAS_avg, clim)
subplot(1,3,3)
imagesc(Cav_Map_RCB_avg, clim)


function [cr,in,out] = compute_CR(Cav_MAP,x,z)
    [X,Z] = meshgrid(x,z);
    % position of flow channel
    xf = 0;
    zf = 40e-3;
    r = 1e-3;

    % Data of Cavitation Map within flow channel
    Cav_Map_in = Cav_MAP;
    Cav_Map_in(sqrt((X-xf).^2 + (Z-zf).^2)>r) = 0;
    % Data of Cavitation Map outside flow channel in 5 mm radius
    Cav_Map_out = Cav_MAP;
    Cav_Map_out(sqrt((X-xf).^2 + (Z-zf).^2)<=r) = 0;
    Cav_Map_out(sqrt((X-xf).^2 + (Z-zf).^2)>5e-3) = 0;
    
    in = sum(Cav_Map_in,'all');
    out = sum(Cav_Map_out,'all');
    cr = 10*log10((sum(Cav_Map_in,'all')/sum(Cav_Map_out,'all')));
end
