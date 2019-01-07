% Subroutine LBLVS
function SPL_LBL = LBL_VS(alpha_star, c, U, f, Theta, Phi, L, r, visc, c0, ...
	Itrip)
	% Calculates SPL for LBL-VS noise for 1 segment
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
	Dbar_h = direct_h(M, Theta, Phi);

	% Compute reference Strouhal number
	if Re <= 1.3E+5
		st1prim = 0.18;
	elseif Re < 4.0E+5
		st1prim = 0.001756 * Re^0.3931;
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

	% Compute G1 directivity function
	for i = 1:size(f, 2)
		% Compute Strouhal number ratio
		E = stprim(i) / stpkprm;

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
