% Jin Sing Sia (jssia@uwaterloo.ca), 2018/12/14
% Based on code by Brooks and Marcolini

% -- SOURCE --:
% See NASA Technical Report 19890016302: "Airfoil self-noise and prediction",
% Brooks, Thomas F., Pope, D. Stuart, and Marcolini, Michael A.

% This code is based on provided FORTRAN code for predicting propeller noise in
% Appendix D. Some variable names have changed and code has been vectorized
% where possible.

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