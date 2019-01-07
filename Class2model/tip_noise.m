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