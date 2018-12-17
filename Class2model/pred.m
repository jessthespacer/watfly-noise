clear all; clc; close all;

% Jin Sing Sia (jssia@uwaterloo.ca), 2018/12/14
% Based on code by Brooks and Marcolini

% -- SOURCE --:
% See NASA Technical Report 19890016302: "Airfoil self-noise and prediction",
% Brooks, Thomas F., Pope, D. Stuart, and Marcolini, Michael A.

% This code is based on provided FORTRAN code for predicting propeller noise in
% Appendix D. Some variable names have changed and code has been vectorized
% where possible.

% Num. blade segments (default)
Nseg = 10;

% Create 1/3 octave centered frequency array
f = [100 125 160 200 250 315 400 500 630 800 1000 1250 1600 2000 2500 ...
	 3150 4000 5000 6300 8000 10000 12500 16000 20000 25000 31500 40000];

% Initialize airfoil parameter arrays for each segment and default vals.

% Segment chord length [m]
c = ones(1, Nseg) * 1.0;

% Segment span length [m]
L = ones(1, Nseg) * 0.1;

% Segment to observer distance [m]
r =  ones(1, Nseg) * 1.0;

% Directivity angle [deg]
Theta = ones(1, Nseg) * 90.0;

% Directivity angle [deg]
Phi = ones(1, Nseg) * 90.0;

% Segment angle of attack [deg]
alpha_star = ones(1, Nseg) * 0.0;

% Segment tip bluntness [m]
h = ones(1, Nseg) * 0.0005;

% Bluntness angle [deg]
Psi = ones(1, Nseg) * 14.0;

% Segment freestream velocity [m/s]
U = ones(1, Nseg) * 100.0;

% Tip angle of attack [deg]
alpha_tip = 0;

% Tip lift curve slope
alprat = 1.0;

% Rounded tip
tip_round = false;

% Fluid transport properties
% Kinematic viscosity (nu) [m^2/s]
visc = 1.4529E-5;

% Speed of sound [m/s]
c0 = 340.46;

% Simulation settings
% Trip boundary layer
Itrip = false;

% Compute laminar boundary layer (LBL) noise
Ilam = false;

% Compute turbulent boundary layer (TBL-TE) noise
Iturb = false;

% Compute bluntness noise
Iblunt = false;

% Compute tip noise
Itip = false;

% --- Initialize your own settings here ---

% Initialize aeroacoustic parameter/result arrays for each frequency
% Strouhal number
St = zeros(size(f));

% SPL from laminar BL vortex shedding
SPL_LBL = zeros(size(f));

% SPL from turbulent BL vortex shedding
SPL_TBL = zeros(size(f));

% SPL from TBL-TE
SPL_P = zeros(size(f));

% SPL from TBL-TE
SPL_S = zeros(size(f));

% SPL from TBL-TE
SPL_ALPH = zeros(size(f));

% SPL from tip noise
SPL_TIP = zeros(size(f));

% SPL from blunt noise
SPL_BLUNT = zeros(size(f));

% Total SPL matrix
SPL = zeros(7, size(f)(2));

% Pressures for TBL-TE
P = zeros(7, size(f)(2));

% Prediction code
% Iterate over all segments to predict SPL
for i = 1:Nseg
	if Ilam
		SPL_LBL(i, :) = LBL_VS(alpha_star(i), c(i), U(i), f, Theta(i), ...
			Phi(i), L(i), r(i), visc, c0);
	end

	if Iturb
	% TODO
		[SPL_P(i, :), SPL_S(i, :), SPL_ALPH(i, :), SPL_TBL(i, :)] = ...
		TBL_TE(alpha_star(i), c(i), U(i), f, Itrip, Theta(i), Phi(i), ...
			L(i), r(i), visc, c0);
	end

	if Iblunt
	% TODO, TODO, TODO TODO TODO TODO TODOOOOOOOO TODO DO DO
		SPL_BLUNT(i, :) = BLUNT(alpha_star(i), c(i), U(i), f, Itrip, ...
			Theta(i), 	Phi(i), L(i), r(i), h(i), Psi(i), visc, c0);
	end
end

if Itip
	% TODO, TODO, TODO TODO TODO TODO TODOOOOOOOO!
	% Must ensure that only last segment is used
	SPL_TIP(end) = TIP_NOISE(alpha_tip, alprat, c(end), U(end), f, Theta, ...
		Phi, r(end), visc, c0, tip_round);
	end

% Pressure contribution, use mean-square pressure basis (vectorized)
if Ilam
	P(5, :) += 10.^(SPL_LBL / 10);
end

if Iturb
	P(1, :) += 10.^(SPL_P / 10);
	P(2, :) += 10.^(SPL_S / 10);
	P(3, :) += 10.^(SPL_ALPH / 10);
end

if Iblunt
	P(6, :) += 10.^(SPL_BLUNT / 10);
end

if Itip 
	% Just the tip...
	P(7, end) += 10.^(SPL_TIP / 10);
end

% Compute total pressure for all mechanisms for each segment; this basically
% flattens the matrix
P(4, :) = sum([P([1:3 5:6], :)]);

% For segments/freqs where P != 0, convert to SPL (vectorized)
SPL(P ~= 0) = 10 .* log10(P(P ~= 0));

% DISCLAIMER: Limitations of this code include not knowing the mean airspeed
% velocity of an unladen swallow, whether African or European.