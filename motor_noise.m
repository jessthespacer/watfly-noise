% See https://ntrs.nasa.gov/search.jsp?R=20160014910
% See http://www.sengpielaudio.com/calculator-soundpower.htm
% kW: Consumed power by motor
% RPM: Motor RPM
% A: Area of motor bounding cylinder	[m^2]
% h: Motor height						[m]
% D: Motor diameter						[m]
% r: Observer distance					[m]
% Q: Directivity function:				[dimless]
% 	Q = 1 for full sphere
%	Q = 2 for half-sphere
% 	Q = 4 for quarter-sphere
% 	Q = 8 for eighth-sphere

% Area of cylinder
A_cyl = @(h, D) pi / 4 * D.^2 * 2 + pi * D * h;

% Sound power level (SWL)
SWL = @(kW, RPM, A) 27 + 10 * log10(kW) + 15 * log10(RPM) + 10 * log10(A);

% Convert SWL to SPL
SPL = @(SWL, r, Q) SWL - abs(10 .* log10(Q / (4 * pi * r^2)));

% EMRAX 208 characteristics:
h = 0.085
D = 0.208
A = A_cyl(h, D)

% Operational characteristics:
P = 55			% [kW]
N = 5000		% [RPM]
Q = 2
r = 15.24		% [m]

% SWL calculation
PWL = SWL(P, N, A);

% Convert to SPL
OASPL_r = SPL(PWL, r, Q);
OASPL_0 = SPL(PWL, 1, Q);

disp(['OASPL @ 1 m = ' num2str(OASPL_0) ' dBA'])
disp(['OASPL @ ' num2str(r) ' m = ' num2str(OASPL_r) ' dBA']);