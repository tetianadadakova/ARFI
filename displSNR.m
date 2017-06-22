function [SNRdmean] = displSNR (delta, toffset, xmax, taur, taud, T2s, D, MEGtype)
% displSNR calculated ARF-related mean displacement SNR (SNRdmean) for a range of
% mechanical and MR tissue constants as well as a range or MEG durations
% and HIFU offset times
%
% --- Function input ---
% delta    -TOTAL duration of MEG (Fig. 1), in ms
% toffset  -time when the HIFU starts relative to the MEG start (Fig. 1), ms
% xmax     -maximum ARF-related tissue displacement (would be achieved if the HIFU duration is infinite), in mm
% taur     -rise  time constant of the overdamped responce model of the displacement (Eq. 2), ms
% taud     -decay time constant of the overdamped responce model of the displacement (Eq. 2), ms
% T2s      -2star transverse relaxation time, ms
% D        -diffusion coefficient, mm^2/ms
% MEGtype  -MEG type: 1 corresponds to bipolar MEG and 2 correspond to tripolar MEG

% Contact: Tetiana Dadakova, tetiana.dadakova@uniklinik-freiburg.de or tetiana.d@gmail.com
% 1. Department of Radiology - Medical Physics, Medical Center - University of Freiburg, Germany
% 2. Faculty of Medicine, University of Freiburg, Germany
%

%% Time samples
dt = 0.05; % sampling step of the simulation
t = 0:dt:100; % ms, simulation os done over 100 ms

%% MR related constants
gamma = 42.576e3;  %kHz/T or 1/(ms*T), gyromagnetic ratio for the hydrogen protons
gammarad = gamma * 2 * pi; % in rad/ms*T, gyromagnetic ratio for the hydrogen protons in radians
G0 = 40e-6; %T/mm, MEG amplitude
risetime = 0.25; %ms, MEG rise time

%% Variables that depend on delta
TE = 4 + delta;

%% Displacement
disp('Calculating displacement x(t) - this takes some time ...')
tic
[tg, xmaxg, taurg, taudg, toffsetg, deltag] = ndgrid(t, xmax, taur, taud, toffset, delta); % replicates vectors t, xmax, taur, taud, toffset, and delta to crease a full grid. This is done in order to calculatre the displacement x for all combinations of values of these vectors

% switch MEGtype
% case 1 % Bipolar MEG
%     toffsetg = toffsetg + deltag./2; % the displacement starts with the second lobe of bipolar MEG
% case 2 % Three-lobed MEG
%     toffsetg = toffsetg + ((deltag - 6*risetime - 2*risetime) ./ 4 + 2*risetime); % the displacement starts with the second lobe of three-lobed MEG        
% end

toffg = toffsetg + deltag/2; % time when the HIFU ends. HIFU lasts for 1/2 of MEG duration

% Displacement calculation according to Eq 2:
x = (heaviside(tg-toffsetg) .* heaviside(toffg-tg)) .*     (-xmaxg .* (1-exp(-((tg-toffsetg) ./ taurg)))) + ...
    (heaviside(tg-toffg).* heaviside(tg(end)+1-tg)) .*     (-xmaxg .* (1-exp(-((toffg-toffsetg) ./ taurg)))) ...
    .* exp(-((tg-toffg) ./ taudg));
toc

%% MEG
switch MEGtype
case 1 % Bipolar MEG
    disp('Calculating bipolar MEG ...')
    tic
    [tgG, deltagG] = ndgrid(t, delta);
    
    flattopgG = deltagG./2 - 2*risetime; % ms, flat top duration of each lobe of MEG
    G = (heaviside(tgG - 0) .* heaviside(risetime - tgG)) ...
        .* G0 / risetime .* (tgG - 0) + ...
        (heaviside(tgG - risetime) .* heaviside((risetime + flattopgG) - tgG)) ...
        .* G0 + ...
        (heaviside(tgG - (risetime + flattopgG)) .* heaviside((3*risetime + flattopgG) - tgG)) ...
        .* (-G0 / risetime .* (tgG - (risetime + flattopgG)) + G0) + ...
        (heaviside(tgG - (3*risetime + flattopgG)) .* heaviside((3*risetime + 2*flattopgG) - tgG)) ...
        .* (-G0) + ...
        (heaviside(tgG - (3*risetime + 2*flattopgG)) .* heaviside((4*risetime + 2*flattopgG) - tgG)) ...
        .* (G0/risetime .* (tgG - (3*risetime + 2*flattopgG) ) - G0);
    toc

case 2 % Tripolar MEG
    disp('Calculating tripolar MEG ...')
    tic
    [tgG, deltagG] = ndgrid(t, delta);
    
    flattopS = (deltagG - 6*risetime - 2*risetime) / 4; % ms, flat top duration of first and third lobe of MEG
    flattopB = 2*flattopS + 2*risetime; % ms, flat top duration of middle lobe of MEG
    G = (heaviside(tgG - 0) .* heaviside(risetime-tgG)) ...
        .* G0 / risetime .* (tgG - 0) + ...
        (heaviside(tgG - risetime) .* heaviside((risetime + flattopS) - tgG)) ...
        .* G0 + ...
        (heaviside(tgG - (risetime + flattopS)) .* heaviside((3*risetime + flattopS) - tgG)) ...
        .* (-G0 / risetime .* (tgG - (risetime + flattopS) ) + G0) + ...
        (heaviside(tgG - (3*risetime + flattopS)) .* heaviside((3*risetime + flattopS + flattopB) - tgG)) ...
        .* (-G0) + ...
        (heaviside(tgG - (3*risetime + flattopS + flattopB)) .* heaviside((5*risetime + flattopS + flattopB) - tgG)) ...
        .* (G0 / risetime .* (tgG - (3*risetime + flattopS + flattopB)) - G0) + ...
        (heaviside(tgG - (5*risetime + flattopS + flattopB)) .* heaviside((5*risetime + 2*flattopS + flattopB) - tgG)) ...
        .* G0 + ...
        (heaviside(tgG - (5*risetime + 2*flattopS + flattopB)) .* heaviside((6*risetime + 2*flattopS + flattopB) - tgG)) ...
        .* (-G0 / risetime .* (tgG - (5*risetime + 2*flattopS + flattopB) ) + G0);
    toc
end

% % Optional: Plot MEG amplitude vs time
% figure
% plot(t(1:31/dt),G(1:31/dt, 8),'o-')
% ylim([-45e-6 45e-6])

%% Phase
disp('Calculating displacement phase ...')
tic
tempG = repmat(reshape(G, [size(t,2) 1 1 1 1 size(delta,2)]), ...
    [1 size(xmax,2) size(taur,2) size(taud,2) size(toffset,2) 1]); % reshape and replicate G matrix to make it the same size as x matrix, because they will be multilied element-wise for phase calculation

% Phase calculation according to Eq 1:
Phase = (gammarad) * dt * sum(tempG .* x, 1);
Phase = squeeze(Phase);
toc

%% b-value
disp('Calculating b-value ...')

% b-value calculation according to Eq 4:
tic
b = (gammarad)^2 * (dt)^3 * sum(cumsum(G,1).^2);
toc

%% SNR
disp('Calculating SNRd ...')

tic
% reshape and replicate matrixes to make all of them the same dimensions to be used in SNRd calculations:
tempPhase = repmat(Phase, [1 1 1 1 1 size(D,2) size(T2s,2)]);
tempb = repmat(reshape(b, [1 1 1 1 size(delta,2) 1 1]),   ... % b-value
    [size(xmax,2) size(taur,2) size(taud,2) size(toffset,2) 1             size(D,2) size(T2s,2)  ]); 
tempD = repmat(reshape(D, [1 1 1 1 1 size(D,2) 1]),       ... % diffusion coefficient
    [size(xmax,2) size(taur,2) size(taud,2) size(toffset,2) size(delta,2) 1         size(T2s,2)  ]);
tempTE = repmat(reshape(TE, [1 1 1 1 size(delta,2) 1 1]), ... % TE
    [size(xmax,2) size(taur,2) size(taud,2) size(toffset,2) 1             size(D,2) size(T2s,2)  ]);
tempT2s = repmat(reshape(T2s, [1 1 1 1 1 1 size(T2s,2)]), ... % T2 star
    [size(xmax,2) size(taur,2) size(taud,2) size(toffset,2) size(delta,2) size(D,2) 1            ]);

% SNRd calculation according to Eq 3:
SNRd = tempPhase .* squeeze(exp(-tempb .* tempD)) .* squeeze(exp(-tempTE ./ tempT2s));
SNRd = abs(SNRd);
toc

% mean SNR over xmax, taur, taud, T2s, and D:
SNRdmean = squeeze(mean (mean (mean (mean (mean (SNRd, 7), 6), 3), 2), 1)  );

disp('Simulations of SNRdmean is completed.')

end





