% This code is contains an example on how to use the displSNR function to
% calculate ARF-related mean displacement SNR (SNRdmean) for a range of
% mechanical and MR tissue constants as well as a range or MEG durations
% and HIFU offset times

% More information in the following paper
% Dadakova T, et al. Magn Reson Med (2017). Optimization of acoustic radiation force imaging: Influence of timing parameters on sensitivity

% Contact: Tetiana Dadakova, tetiana.dadakova@uniklinik-freiburg.de or tetiana.d@gmail.com
% 1. Department of Radiology - Medical Physics, Medical Center - University of Freiburg, Germany
% 2. Faculty of Medicine, University of Freiburg, Germany

clear all
% close all

%% Tissue constants
% Ranges of the values
xmax = 1e-3 : 5e-3 : 40e-3; %mm, maximum ARF-related tissue displacement (would be achieved if the HIFU duration is infinite) 
taur = 1 : 1 : 15; %ms, rise time constant of the overdamped responce model of the displacement (Eq. 2)
taud = 1 : 1 : 15; %ms, decay time constant of the overdamped responce model of the displacement (Eq. 2)
T2s = 20 : 5 : 50; %ms, T2star transverse relaxation time
D = 0.2e-6 :  0.2e-6 : 2e-6; %mm^2/ms, diffusion coefficient

% % Phantom values
% xmax = 6.8e-3; %mm
% taur = 4; %ms
% taud = 2; %ms
% T2s = 65; %ms
% D = 1.0e-6; %mm^2/ms

%% Experimentally controlled (changeable) variables
delta = 2 : 2 : 60; %ms, TOTAL duration of MEG (Fig. 1)
toffset = -60 : 2 : 30; %ms, time when the HIFU starts relative to the MEG start (Fig. 1)

%% Call function displSNR, which calculates the displacement SNR for a mentione d range of paramenters
SNRdmean = displSNR (delta, toffset, xmax, taur, taud, T2s, D, 1);

%% Plot the displacement SNR as a function of MEG duration delta and HIFU offset time 
figure; contourf(delta, toffset, SNRdmean)
colorbar
xlabel('MEG duration (ms)', 'FontSize',14)
ylabel('HIFU offset time (ms)', 'FontSize',14)
caxis([0 0.6])
