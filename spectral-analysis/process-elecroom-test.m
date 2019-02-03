clear all; clc; close all;

% Background noise spectrum
fbacg = 'Spectra/test-elecroom-12-jan-2019/elecroom-bacgnoise.csv'

% Motor noise spectrum
fmotor = 'Spectra/test-elecroom-12-jan-2019/theta=0,r=1ft-bacg+motor-b5-2500rpm.txt'

% Overall noise spectrum
ftotal = 'Spectra/test-elecroom-12-jan-2019/theta=0,r=1ft-bacg+motor+prop-b5-2500rpm.txt'

% Load spectra
bacg = readSpec(fbacg, 'rew');
motor = readSpec(fmotor, 'rew');
total = readSpec(ftotal, 'rew');

% Get propeller noise
% propeller = removeNoise(total, motor, true);

% motor = abs(removeNoise(motor, bacg, true));

% Plot original recording
semilogx(total(:, 1), total(:, 2)); hold on;

% % Plot bacg + motor noise
% semilogx(motor(:, 1), motor(:, 2));

% % Plot background noise
% semilogx(bacg(:, 1), bacg(:, 2));

% semilogx(motor(:, 1), motor(:, 2));

% Plot corrected noise
% semilogx(propeller(:, 1), propeller(:, 2));

% Label axes
% legend('Background + Motor + Propeller', 'Background + Motor', 'Background');
% xlabel('Frequency [Hz]');
% ylabel('SPL [dB]');