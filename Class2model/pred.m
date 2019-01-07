clear all; clc; close all;

% Jin Sing Sia (jssia@uwaterloo.ca), 2018/12/14
% Based on code by Brooks and Marcolini

% -- SOURCE --:
% See NASA Technical Report 19890016302: "Airfoil self-noise and prediction",
% Brooks, Thomas F., Pope, D. Stuart, and Marcolini, Michael A.

% This code is based on provided FORTRAN code for predicting propeller noise in
% Appendix D. Some variable names have changed and code has been vectorized
% where possible.

% Num. blade segments
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

Nseg = 1
c = 0.3048;
L = 0.4572;
r = 1.22;
Theta = 90;
Phi = 90;
alpha_star = 1.516;
U = 71.3;
Itrip = 0;
Ilam = 1;
Iturb = 1;

Theta *= pi / 180;
Phi *= pi / 180;

% Initialize aeroacoustic parameter/result arrays for each frequency
size_parray = [Nseg, size(f, 2)];

% SPL from laminar BL vortex shedding
SPL_LBL = zeros(size_parray);

% SPL from turbulent BL vortex shedding
SPL_TBL = zeros(size_parray);

% SPL from TBL-TE
SPL_P = zeros(size_parray);

% SPL from TBL-TE
SPL_S = zeros(size_parray);

% SPL from TBL-TE
SPL_ALPH = zeros(size_parray);

% SPL from tip noise
SPL_TIP = zeros(size_parray);

% SPL from blunt noise
SPL_BLUNT = zeros(size_parray);

% Total SPL matrix
SPL = zeros(7, size(f, 2));

% Pressures for TBL-TE
P = zeros(7, size(f, 2));

% Prediction code
% Iterate over all segments to predict SPL
for I = 1:Nseg
	if Ilam
		SPL_LBL(I, :) = LBL_VS(alpha_star(I), c(I), U(I), f, Theta(I), ...
			Phi(I), L(I), r(I), visc, c0, Itrip);
	end

	if Iturb
		[SPL_P(I, :), SPL_S(I, :), SPL_ALPH(I, :), SPL_TBL(I, :)] = ...
		TBL_TE(alpha_star(I), c(I), U(I), f, Itrip, Theta(I), Phi(I), ...
			L(I), r(I), visc, c0);
	end

	if Iblunt
		SPL_BLUNT(I, :) = blunt(alpha_star(I), c(I), U(I), f, Itrip, ...
			Theta(I), Phi(I), L(I), r(I), h(I), Psi(I), visc, c0);
	end

	% Iterate over all frequencies
	for J = 1:size(f, 2)
		% Add segment's contribution on mean-square pressure basis

		if Ilam
			P(5, J) += 10^(SPL_LBL(I, J) / 10);
		end

		if Iturb
			P(1, J) += 10^(SPL_P(I, J) / 10);
			P(2, J) += 10^(SPL_S(I, J) / 10);
			P(3, J) += 10^(SPL_ALPH(I, J) / 10);
		end

		if Iblunt
			P(6, J) += 10^(SPL_BLUNT(I, J) / 10);
		end

		% Tip noise for last segment only
		if Itip && I == Nseg
			P(7, J) += 10^(SPL_TIP(I, J) / 10);
		end

		% Compute total pressure for the segment for all mechanisms
		P(4, J) = sum(P([1 2 3 5 6 7], J));
	end
end

% For segments/freqs where P != 0, convert to SPL (vectorized)
SPL(P ~= 0) = 10 .* log10(P(P ~= 0));

% DISCLAIMER: Limitations of this code include not knowing the airspeed
% velocity of an unladen swallow, whether African or European.