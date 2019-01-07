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