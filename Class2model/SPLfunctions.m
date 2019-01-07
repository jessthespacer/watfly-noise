% Jin Sing Sia (jssia@uwaterloo.ca), 2018/12/14
% Based on code by Brooks and Marcolini

% -- SOURCE --:
% See NASA Technical Report 19890016302: "Airfoil self-noise and prediction",
% Brooks, Thomas F., Pope, D. Stuart, and Marcolini, Michael A.

% This code is based on provided FORTRAN code for predicting propeller noise in
% Appendix D. Some variable names have changed and code has been vectorized
% where possible.

% Subroutine AMIN
function amin_a = amin(a)
	% Defines curve fit corresponding to a-curve for minimum allowed Re
	x1 = abs(a);
	if x1 <= 0.204
		amin_a = sqrt(67.552 - 886.788 * x1^2) - 8.219;
	elseif x1 <= 0.244
		amin_a = -32.665 * x1 + 3.981;
	else
		amin_a = -142.795 * x1^3 + 103.656 * x1^2 - 57.757 * x1 + 6.006;
	end
end

% Suburoutine AMAX
function amax_a = amax(a)
	% Defines curve fit corresponding to a-curve for maximum allowed Re
	x1 = abs(a);
	if x1 <= 0.13
		amax_a = sqrt(67.552 - 886.788 * x1^2) - 8.219;
	elseif x1 <= 0.321
		amax_a = -15.901 * x1 + 1.098;
	else
		amax_a = -4.669 * x1^3 + 3.491 * x1^2 - 16.699 * x1 + 1.149;
	end
end

% Subroutine BMIN
function bmin_b = bmin(b)
	% Defines curve fit corresponding to b-curve for minimum allowed Re
	x1 = abs(b);
	if x1 <= 0.13
		bmin_b = sqrt(16.888 - 886.788 * x1^2) - 4.109;
	elseif x1 <= 0.145
		bmin_b = -83.607 * x1 + 8.138;
	else
		bmin_b = -817.81 * x1^3 + 355.21 * x1^2 - 135.024 * x1 + 10.619;
	end
end

% Subroutine BMAX
function bmax_b = bmax(b)
	% Defines curve fit corresponding to b-curve for maximum allowed Re
	x1 = abs(b);
	if x1 <= 0.1
		bmax_b = sqrt(16.888 - 886.788 * x1^2) - 4.109;
	elseif x1 <= 0.187
		bmax_b = -31.313 * x1 + 1.854;
	else
		bmax_b = -80.541 * x1^3 + 44.174 * x1^2 - 39.381 * x1 + 2.344;
	end
end

% Subroutine A0COMP
function a0 = a0comp(Re)
	% Determines where the a-curve is -20 dB
	if Re < 9.52E+04
		a0 = 0.57;
	elseif Re < 8.57E+05
		a0 = -9.57E-13 * (Re - 8.57E+05)^2 + 1.13;
	else
		a0 = 1.13;
	end
end

% Subroutine DIRECTH
function Dbar_h = direct_h(M, Theta, Phi)
	% Computes high frequency directivity function for a given observer
	degrad = 0.017453;

	Mc = 0.8 * M;
	ThetaR = Theta * degrad;
	PhiR = Phi * degrad;

	Dbar_h = 2 * sin(ThetaR / 2)^2 * sin(PhiR)^2 / ...
		((1 + M * cos(ThetaR)) * (1 + (M - Mc) * cos(ThetaR))^2);
end

% Subroutine DIRECTL
function Dbar_l = direct_l(M, Theta, Phi)
	% Computes low frequency directivity function for a given observer
	degrad = 0.017453;

	Mc = 0.8 * M;
	ThetaR = Theta * degrad;
	PhiR = Phi * degrad;

	Dbar_l = (sin(ThetaR) * sin(PhiR))^2 / (1 + M * cos(ThetaR))^4;
end

% Subroutine blunt
function G5 = G5comp(hdstar, eta)
	if hdstar < 0.25
		mu = 0.1211;
	elseif hdstar <= 0.62
		mu = -0.2175 * hdstar + 0.1755;
	elseif hdstar < 1.15
		mu = -0.0308 * hdstar + 0.0596;
	else
		mu = 0.0242;
	end

	if hdstar <= 0.02
		M = 0.0;
	elseif hdstar < 0.5
		M = 68.724 * hdstar - 1.35;
	elseif hdstar <= 0.62
		M = 308.475 * hdstar - 121.23;
	elseif hdstar <= 1.15
		M = 224.811 * hdstar - 69.354;
	elseif hdstar < 1.2
		M = 1583.28 * hdstar - 1631.592;
	else
		M = 268.344;
	end

	if M < 0
		M = 0;
	end

	eta0 = sqrt((M * M * mu^4) / (6.25 + M * M * mu * mu));

	K = 2.5 * sqrt(1 - (eta0 / mu)^2) - 2.5 - M * eta0;

	if eta <= eta0
		G5 = M * eta + K;
	elseif eta > eta0 && eta <= 0
		G5 = 2.6 * sqrt(1 - (eta / mu)^2) - 2.5;
	elseif eta > 0 && eta <= 0.03616
		G5 = sqrt(1.5625 - 1194.99 * eta^2) - 1.25;
	else
		G5 = -155.543 * eta + 4.375;
	end
end

% Subroutine THICK
function [deltaP, dstrp, dstrs] = BL_thick(c, U, alpha_star, Itrip, c0, visc)
	% Compute boundary layer thickness
	M = U / c0;
	Re = U * c / visc;

	delta0 = 10^(1.6569 - 0.9045 * log10(Re) + 0.0596 * (log10(Re))^2) * c;

	if Itrip == 2
		delta0 *= 0.6;
	end

	% Compute pressure-side boundary layer thickness
	deltaP = 10^(-0.04175 * alpha_star + 0.00106 * alpha_star^2) * delta0;

	% Compute zero AOA displacement thickness
	if Itrip == 1 || Itrip == 2
		if Re <= 0.3E+06
			dstr0 = 0.0601 * Re^(-0.114 * c);
		else
			dstr0 = 10^(3.411 - 1.5397 * log10(Re) + 0.1059 * (log10(Re))^2)*c;
		end
		if Itrip = 2
			dstr0 *= 0.6;
		end
	else
		dstr0 = 10^(3.0187 - 1.5397 * log10(Re) + 0.1059 * (log10(Re))^2) * c;
	end

	% Pressure side displacement thickness
	dstrp = 10^(-0.0432 * alpha_star + 0.00113 * alpha_star^2) * dstr0;
	if Itrip == 3
		dstrp *= 1.48;
	end

	% Suction side displcement thickness
	if Itrip == 1
		if alpha_star <= 5
			dstrs = 10^(0.0679 * alpha_star) * dstr0;
		elseif alpha_star <= 12.5
			dstrs = 0.381 * 10^(0.1516 * alpha_star) * dstr0;
		else
			dstrs = 14.296 * 10^(0.0258 * alpha_star) * dstr0;
		end
	else
		if alpha_star <= 7.5
			dstrs = 10^(0.0679 * alpha_star) * dstr0;
		elseif alpha_star <= 12.5
			dstrs = 0.0162 * 10^(0.3066 * alpha_star) * dstr0;
		else
			dstrs = 52.42 * 10^(0.0258 * alpha_star) * dstr0;
		end
	end
end

function SPL_BLUNT = blunt(alpha_star, c, U, f, Itrip, Theta, Phi, L, r, ...
	h, Psi, visc, c0)
	% Calculates SPL for blunt noise
	SPL_BLUNT = zeros(size(f));
	stppp = zeros(size(f));

	M = U / c0;
	Re = U * c / visc;

	[deltaP, dstrp, dstrs] = BL_thick(c, U, alpha_star, Itrip, c0, visc);

	% Compute average displacement thickness
	dstravg = (dstrs + dstrp) / 2;
	hdstar = h / dstravg;
	dstarh = 1 / hdstar;

	% Compute directivity function
	Dbar_h = direct_h(M, Theta, Phi);

	% Compute peak Strouhal number
	aterm = 0.212 - 0.0045 * Psi;

	if hdstar >= 0.2
		stpeak = aterm / (1 + 0.235 * dstarh - 0.0132 * dstarh^2);
	else
		stpeak = 0.1 * hdstar + 0.095 - 0.00243 * Psi;
	end

	% Compute scaled spectrum level
	if hdstar <= 5
		G4 = 17.5 * log10(hdstar) + 157.5 - 1.114 * Psi;
	else
		G4 = 169.7 - 1.114 * Psi;
	end

	% For each frequency, compute spectrum shape referenced to 0 dB
	for i = 1:size(f, 2)
		% I think NASA just gave up trying to explain the code a long, long
		% time ago...
		stppp(i) = f(i) * h / U;
		eta = log10(stppp(i) / stpeak);

		hdstarl = hdstar;
		G514 = G5comp(hdstarl, eta);
		hdstarp = 6.724 * hdstar^2 - 4.019 * hdstar + 1.107;
		G50 = G5comp(hdstarp, eta);

		G5 = G50 + 0.0714 * Psi * (G514 - G50);

		if G5 > 0
			G5 = 0;
		end

		F4temp = G5comp(0.25, eta);

		if G5 > F4temp
			G5 = F4temp;
		end

		Scale = 10 * log10(M^5.5 * h * Dbar_h * L / r^2);

		SPL_BLUNT(i) = G4 + G5 + Scale;
	end
end

% Subroutine LBLVS
function SPL_LBL = LBL_VS(alpha_star, c, U, f, Theta, Phi, L, r, visc, c0)
	% Calculates SPL for LBL-VS noise
	% Taken from p. 118
	SPL_LBL = zeros(size(f));

	% Strouhal number based on pressure side BL thickness
	stprim = zeros(size(f));

	% Calculate Reynolds and Mach numbers
	M = U / c0;
	Re = U * c / visc;

	% Calculate boundary layer thicknesses
	% deltaP: Pressure side BL thickness
	% dstrp: Pressure side BL displacement thickness
	% dstrs: Suction side BL displacement thickness
	% TODO: Implement
	[deltaP, dstrp, dstrs] = BL_thick(c, U, alpha_star, Itrip, c0, visc);

	% Compute directivity function
	% Dbar_h: High frequency directivity
	% TODO: Implement
	Dbar_h = direct_h(M, Theta, Phi)

	% Compute reference Strouhal number
	% Why the hell is it called st1prim anyway?
	if Re <= 1.3E+5
		st1prim = 0.18;
	elseif Re < 4.0E+5
		st1prim = 0.001756 * Re.^0.3931;
	else
		st1prim = 0.28;
	end

	% Peak Strouhal number
	stpkprm = 10^(-0.04 * alpha_star) * st1prim;

	% Compute reference Reynolds number
	if alpha_star <= 3
		Re0 = 10.^(0.215 .* alpha_star(alpha_star <= 3) + 4.978);
	else
		Re0 = 10.^(0.12 .* alpha_star(alpha_star > 3) + 5.263);
	end

	% Compute peak scaled spectrum level
	D = Re / Re0;

	% Compute G2 directivity function
	if D <= 0.3237
		G2 = 77.852 * log10(D) + 15.328;
	elseif D <= 0.5689
		G2 = 65.188 * log10(D) + 9.125;
	elseif D <= 1.7579
		G2 = -114.052 * (log10(D))^2;
	elseif D <= 3.0889
		G2 = -65.188 * log10(D) + 9.125;
	else
		G2 = -77.852 * log10(D) + 15.328;
	end

	G3 = 171.04 - 3.03 * alpha_star;

	Scale = 10 * log10(deltaP * M^5 * Dbar_h * L / r^2);

	% Compute scaled SPLs for each Strouhal number (vectorized)
	stprim = f * deltaP / U;

	% Compute Strouhal number ratio
	E = stprim ./ stpkprm;

	% Compute G1 directivity function
	for i = 1:size(f, 2)
		if E < 0.5974
			G1 = 39.8 * log10(E) - 11.12;
		elseif E < 0.8545
			G1 = 98.409 * log10(E) + 2.0;
		elseif E < 1.17
			G1 = -5.076 + sqrt(2.484 - 506.25 * (log10(E))^2);
		elseif E < 1.674
			G1 = -98.409 * log10(E) + 2.0;
		else
			G1 = -39.8 * log10(E) - 11.12;
		end

		% Compute scaled SPL for each Strouhal number
		SPL_LBL(i) = G1 + G2 + G3 + Scale;
	end
end

% Subroutine TBLTE
% Look, it's a long boi...
function [SPL_P, SPL_S, SPL_ALPH, SPL_TBL] = TBL_TE(alpha_star, c, U, f, ... Itrip, Theta, Phi, L, r, visc, c0)
	% Calculates SPL for TBL-TE noise

	% SPL vectors for each mechanism
	SPL_TBL = zeros(size(f));
	SPL_P = zeros(size(f));
	SPL_S = zeros(size(f));
	SPL_ALPH = zeros(size(f));

	% Pressure side Strouhal number
	stp = zeros(size(f));

	% Suction side Strouhal number
	sts = zeros(size(f));

	% Compute angle of attack contribution
	computeAOAcont = false;

	% Compute Reynold's number
	Re = U * c / visc;

	% Compute Mach number
	M = U / c0;

	% Compute boundary layer thicknesses
	[deltaP, dstrp, dstrs] = BL_thick(c, U, alpha_star, Itrip, c0, visc);

	% Compute directivity functions
	Dbar_h = direct_h(M, Theta, Phi);
	Dbar_l = direct_l(M, Theta, Phi);	% TODO: Implement

	% Calculate Re based on pressure and suction displacement thicknesses
	rdstrs = dstrs * U / visc;
	rdstrp = dstrp * U / visc;

	% Find peak Strouhal numbers to be used for 'A' and 'B' curve calcs.
	st1 = 0.02 * M^(-0.6);

	if alpha_star <= 1.333
		st2 = st1;
	elseif alpha_star <= 12.5
		st2 = st1 * 10^(0.0054 * (alpha_star - 1.333)^2);
	else
		st2 = 4.72 * st1;
	end

	st1prim = (st1 + st2) / 2;

	a0 = a0comp(Re);		% TODO: Implement
	a02 = a0comp(3 * Re);	% TODO: Implement

	% Evaluate minimum and maximum 'A' curves at A0
	% TODO: Implement
	amin_a0 = amin(a0);
	amax_a0 = amax(a0);
	amin_a02 = amin(a02);
	amax_a02 = amax(a02);

	% Compute 'A' max/min ratio
	ara0 = (20 + amin_a0) / (amin_a0 - amax_a0);
	ara02 = (20 + amin_a02) / (amin_a02 - amax_a02);

	% Compute B0 to be used in 'B' curve calculations
	if Re < 9.52E+04
		b0 = 0.3;
	elseif Re < 8.57E+05
		b0 = -4.48E-13 * (Re - 8.57E+05)^2 + 0.56;
	else
		b0 = 0.56;
	end

	% Evaluate minimum and maximum 'B' curves at B0
	% TODO: Implement
	bmin_b0 = bmin(b0);
	bmax_b0 = bmax(b0);

	% Compute 'B' max/min ratio
	brb0 = (20 + bmin_b0) / (bmin_b0 - bmax_b0);

	% For each center frequency, compute 'A' prediction for pressure side
	stpeak = st1;

	for i = 1:size(f, 2)
		stp(i) = f(i) * dstrp / U;
		a = log10(stp(i) / stpeak);
		amin_a = amin(a);
		amax_a = amax(a);
		aa = amin_a + ara0 * (amax_a - amin_a);

		if Re < 2.47E+05
			K1 = -4.31 * log10(Re) + 156.3;
		elseif Re < 8.0E+05
			K1 = -9.0 * log10(Re) + 181.6;
		else
			K1 = 128.5;
		end

		if rdstrp <= 5000
			delk1 = -alpha_star * (5.29 - 1.43 * log10(rdstrp));
		else
			delk1 = 0.0;
		end
		
		SPL_P(i) = aa + K1 - 3 + ...
			10 * log10(dstrp * M^5 * Dbar_h * L / r^2) + delk1;

		Gamma = 27.094 * M + 3.31;
		Beta = 72.65 * M + 10.74;
		Gamma0 = 23.43 * M + 4.651;
		Beta0 = -34.19 * M - 13.82;

		if alpha_star <= (Gamma0 - Gamma)
			K2 = -1000;
		elseif alpha_star <= (Gamma0 + Gamma)
			K2 = sqrt(Beta^2 - (Beta / Gamma)^2 * (alpha_star - Gamma0)^2) ...
				+ Beta0;
		else
			K2 = -12;
		end

		K2 = K2 + K1;

		sts(i) = f(i) * dstrs / U;
		
		% Check for 'A' computation for suction side
		xcheck = Gamma0;
		computeAOAcont = false;

		if alpha_star >= xcheck || alpha_star > 12.5
			computeAOAcont = true;
		end
		if !computeAOAcont
			a = log10(sts(i) / st1prim);
			amin_a = amin(a);
			amax_a = amax(a);
			aa = amin_a + ara0 * (amax_a - amin_a);

			SPL_S(i) = aa + K1 - 3 + 10^log10(dstrs * M^5 * Dbar_h * L / r^2);

			% Check for 'B' computation for suction side
			b = abs(log10(sts(i) / st2));
			bmin_b = bmin(b);
			bmax_b = bmax(b);
			bb = bmin_b + brb0 * (bmax_b - bmin_b);
			SPL_ALPH(i) = bb + K2 + 10 * log10(dstrs * M^5 * Dbar_h * L / r^2);
		else
			% Drop 'A' computation
			SPL_S(i) = 10 * log10(dstrs * M^5 * Dbar_l * L / r^2);
			SPL_P(i) = 10 * log10(dstrs * M^5 * Dbar_l * L / r^2);

			b = abs(log10(sts(i) / st2));
			amin_b = amin(b);
			amax_b = amax(b);
			bb = amin_b + ara02 * (amax_b - amin_b);

			SPL_ALPH(i) = bb + K2 + 10 * log10(dstrs * M^5 * Dbar_l * L / r^2);
		end

		% Sum all contributions from 'A' and 'B' on both pressure and suction
		% sides on a mean-square pressure basis
		if SPL_P(i) < -100
			SPL_P(i) = -100;
		end
		if SPL_S(i) < -100
			SPL_S(i) = -100;
		end
		if SPL_ALPH(i) < -100
			SPL_ALPH(i) = -100;
		end

		P1 = 10^(SPL_P(i) / 10);
		P2 = 10^(SPL_S(i) / 10);
		P4 = 10^(SPL_ALPH(i) / 10);

		SPL_TBL(i) = 10 * log10(P1 + P2 + P4);
	end
end

% Subroutine TIPNOIS
function SPL_TIP = tip_noise(alpha_tip, alprat, c, U, f, Theta, Phi, r, ...
	visc, c0, tip_round)
	% Calculates tip noise
	SPL_TIP = zeros(size(f));
	alptipp = alpha_tip * alprat;
	M = U / c0;

	Dbar_h = direct_h(M, Theta, Phi);

	if tip_round
		L = 0.008 * alptipp * c;
	else
		if abs(alptipp) <= 2
			L = (0.023 + 0.0169 * alptipp) * c;
		else
			L = (0.0378 + 0.0095 * alptipp) * c;
		end
	end

	MM = (1 + 0.036 * alptipp) * M;
	UM = MM * c0;

	term = M * M * MM^3 * L^2 * Dbar_h / r^2;

	if term ~= 0
		Scale = 10 * log10(term);
	else
		Scale = 0;
	end

	for i = 1:size(f)
		stpp = f(i) * L / UM;
		SPL_TIP(i) = 126 - 30.5 * (log10(stpp) + 0.3)^2 + Scale;
	end
end