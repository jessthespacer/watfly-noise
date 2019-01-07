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