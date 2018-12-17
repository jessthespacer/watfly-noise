clear all; clc; close all;

% Motor noise spectrum
fmotor = 'Spectra/theta=0,r=1ft-bacg+motor-b5-2500rpm.txt'

% Overall noise spectrum
ftotal = 'Spectra/theta=0,r=1ft-bacg+motor+prop-b5-2500rpm.txt'

% Convert to matrix
motor = readSpec(fmotor);
total = readSpec(ftotal);

% Add +100 dB
motor(:, 2) += 100;
total(:, 2) += 100;

% Display sample vals
disp('   freq		ampl')
disp('motor')
motor(1:10, :)
disp('total')
total(1:10, :)

% Sample calculation with sample row
signal = zeros(10, 2);
signal(:, 1) = total(1:10, 1);
signal(:, 2) = 10 .* log10(10.^(total(1:10, 2) ./ 10) - 10.^(motor(1:10, 2) ./ 10));
signal