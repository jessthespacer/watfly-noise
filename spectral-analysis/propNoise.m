clear all; clc; close all;

% Background noise spectrum
% fbacg = 'Spectra/bacg_1.txt'

% Motor noise spectrum
fmotor = 'Spectra/theta=0,r=1ft-bacg+motor-b5-2500rpm.txt'

% Overall noise spectrum
ftotal = 'Spectra/theta=0,r=1ft-bacg+motor+prop-b5-2500rpm.txt'

% Load spectra
% bacg = readSpec(fbacg);
motor = readSpec(fmotor);
total = readSpec(ftotal);

% Subtract bacg + motor noise
correctedSpec = removeNoise(total, motor, true);

% % Plot original recording
% semilogx(total(:, 1), total(:, 2));
% hold on;

% % Plot bacg + motor noise
% semilogx(motor(:, 1), motor(:, 2));

% % Plot corrected noise
% semilogx(correctedSpec(:, 1), correctedSpec(:, 2));

% % Label axes
% legend('Total', 'Background + Motor', 'Propeller');
% xlabel('Frequency [Hz]');
% ylabel('SPL [dB]');