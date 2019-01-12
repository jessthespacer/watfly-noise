% Subroutine TBLTE
% Look, it's a long boi...
function [SPL_P, SPL_S, SPL_ALPH, SPL_TBL] = TBL_TE(alpha_star, c, U, f, Itrip, Theta, Phi, L, r, visc, c0)
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
	Dbar_l = direct_l(M, Theta, Phi);

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

	a0 = a0comp(Re);
	a02 = a0comp(3 * Re);

	% Evaluate minimum and maximum 'A' curves at A0
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
	bmin_b0 = bmin(b0);
	bmax_b0 = bmax(b0);

	% Compute 'B' max/min ratio
	brb0 = (20 + bmin_b0) / (bmin_b0 - bmax_b0);

	% For each center frequency, compute 'A' prediction for pressure side
	stpeak = st1;

	for I = 1:size(f, 2)
		stp(I) = f(I) * dstrp / U;
		a = log10(stp(I) / stpeak);
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
		
		SPL_P(I) = aa + K1 - 3 + ...
			10 * log10(dstrp * M^5 * Dbar_h * L / r^2) + delk1;

		Gamma 	=  27.094 * M + 3.31;
		Beta 	=  72.650 * M + 10.74;
		Gamma0 	=  23.430 * M + 4.651;
		Beta0 	= -34.190 * M - 13.82;

		if alpha_star <= (Gamma0 - Gamma)
			K2 = -1000;
		elseif alpha_star <= (Gamma0 + Gamma)
			K2 = sqrt(Beta^2 - (Beta / Gamma)^2 * (alpha_star - Gamma0)^2) ...
				+ Beta0;
		else
			K2 = -12.0;
		end

		K2 += K1;

		sts(I) = f(I) * dstrs / U;
		
		% Check for 'A' computation for suction side
		xcheck = Gamma0;
		computeAOAcont = false;

		if ((alpha_star >= xcheck) || (alpha_star > 12.5))
			computeAOAcont = true;
		end
		if !computeAOAcont
			a = log10(sts(I) / st1prim);
			amin_a = amin(a);
			amax_a = amax(a);
			aa = amin_a + ara0 * (amax_a - amin_a);

			SPL_S(I) = aa + K1 - 3 + ...
				10 * log10(dstrs * M^5 * Dbar_h * L / r^2);

			% Check for 'B' computation for suction side
			b = abs(log10(sts(I) / st2));
			bmin_b = bmin(b);
			bmax_b = bmax(b);
			bb = bmin_b + brb0 * (bmax_b - bmin_b);
			SPL_ALPH(I) = bb + K2 + 10 * log10(dstrs * M^5 * Dbar_h * L / r^2);
		else
			% Drop 'A' computation
			SPL_S(I) = 10 * log10(dstrs * M^5 * Dbar_l * L / r^2);
			SPL_P(I) = 10 * log10(dstrs * M^5 * Dbar_l * L / r^2);

			b = abs(log10(sts(I) / st2));
			amin_b = amin(b);
			amax_b = amax(b);
			bb = amin_b + ara02 * (amax_b - amin_b);

			SPL_ALPH(I) = bb + K2 + 10 * log10(dstrs * M^5 * Dbar_l * L / r^2);
		end

		% Sum all contributions from 'A' and 'B' on both pressure and suction
		% sides on a mean-square pressure basis
		if SPL_P(I) < -100
			SPL_P(I) = -100;
		end
		if SPL_S(I) < -100
			SPL_S(I) = -100;
		end
		if SPL_ALPH(I) < -100
			SPL_ALPH(I) = -100;
		end

		P1 = 10^(SPL_P(I) / 10);
		P2 = 10^(SPL_S(I) / 10);
		P4 = 10^(SPL_ALPH(I) / 10);

		SPL_TBL(I) = 10 * log10(P1 + P2 + P4);
	end
end
