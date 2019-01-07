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