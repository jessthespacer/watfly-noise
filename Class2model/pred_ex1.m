% Define simulation inputs
Nseg = 1			% Number of segments
c = 0.3048			% Chord length [m]
L = 0.4572			% Span [m]
r = 1.22			% Observer distance [m]
Theta =	90.0		% Observer angle from x-axis [deg]
Phi = 90.0			% Observer angle from y-axis [deg]
alpha_star = 1.516	% Aerodynamic angle of attack [deg]
alpha_tip = 0.0		% Tip flow angle [deg]
U = 71.3			% Free-stream velocity [m/s]
h = 0				% Trailing edge bluntness [m]
Psi = 0				% Trailing edge angle [deg]

% Program parameters and settings
Itrip =	0			% 0: Use untripped BL condition
					% 1: Use tripped BL condition
					% 2: Use lightly tripped BL condition

Ilam = true			% Compute LBL-VS noise
Iturb =	true		% Dompute turbulent TBL-TE noise
Iblunt = false 		% Compute TE bluntness noise
Itip = false 		% Compute tip noise
Iround = false 		% Use rounded tip; if false, use square tip

% Define fluid transport properties
visc = 1.4529E-5	% Kinematic viscosity nu [m^2/s]
c0 = 340.46			% Speed of sound [m/s]