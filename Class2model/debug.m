f = [100 125 160 200 250 315 400 500 630 800 1000 1250 1600 2000 2500 ...
	 3150 4000 5000 6300 8000 10000 12500 16000 20000 25000 31500 40000];

Nseg = 1;
c = ones(1, Nseg) * 1.0;
L = ones(1, Nseg) * 0.1;
r =  ones(1, Nseg) * 1.0;
Theta = ones(1, Nseg) * 90.0;
Phi = ones(1, Nseg) * 90.0;
alpha_star = ones(1, Nseg) * 0.0;
h = ones(1, Nseg) * 0.0005;
Psi = ones(1, Nseg) * 14.0;
U = ones(1, Nseg) * 100.0;
alpha_tip = 0;
alprat = 1.0;
tip_round = false;
visc = 1.4529E-5;
c0 = 340.46;
Itrip = false;
Ilam = false;
Iturb = false;
Iblunt = false;
Itip = false;

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

size_parray = [Nseg, size(f)];
SPL_LBL = zeros(size_parray);
SPL = zeros(7, size(f));
P = zeros(7, size(f));

% Code for segment 1, compute laminar
I = 1;
SPL_LBL(I, :) = LBL_VS(alpha_star(I), c(I), U(I), f, Theta(I), ...
			Phi(I), L(I), r(I), visc, c0, Itrip);

% Correct SPL_LBL vector
SPL_LBL_CORRECT = [-17.142 -13.285 -9.018 -5.161 -1.304 -2.690 6.82 ...
	10.677 14.671 18.801 22.658 26.515 30.782 37.725 47.262 48.959 ...
	41.796 32.428 28.433 24.304 20.477 16.590 12.323 8.466 4.609 ...
	0.614 -3.515];

result = [SPL_LBL(I, :); SPL_LBL_CORRECT]';
[f' result]
residual = norm(result(:, 1) - result(:, 2), 2)

semilogx(f, result(:, 1), "b-o", f, result(:, 2), "r-o", "linewidth", 3)