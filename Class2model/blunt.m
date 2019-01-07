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