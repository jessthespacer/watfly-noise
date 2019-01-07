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