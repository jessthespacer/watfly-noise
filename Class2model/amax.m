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