% Subroutine THICK
function [deltaP, dstrp, dstrs] = BL_thick(c, U, alpha_star, Itrip, c0, visc)
	% Compute boundary layer thickness
	M = U / c0;
	Re = U * c / visc;

	delta0 = 10^(1.6569 - 0.9045 * log10(Re) + 0.0596 * (log10(Re))^2) * c;

	if Itrip == 2
		delta0 *= 0.6;
	end

	% Compute pressure-side boundary layer thickness
	deltaP = 10^(-0.04175 * alpha_star + 0.00106 * alpha_star^2) * delta0;

	% Compute zero AOA displacement thickness
	if Itrip == 1 || Itrip == 2
		if Re <= 0.3E+06
			dstr0 = 0.0601 * Re^(-0.114 * c);
		else
			dstr0 = 10^(3.411 - 1.5397 * log10(Re) + 0.1059 * (log10(Re))^2)*c;
		end
		if Itrip == 2
			dstr0 *= 0.6;
		end
	else
		dstr0 = 10^(3.0187 - 1.5397 * log10(Re) + 0.1059 * (log10(Re))^2) * c;
	end

	% Pressure side displacement thickness
	dstrp = 10^(-0.0432 * alpha_star + 0.00113 * alpha_star^2) * dstr0;
	if Itrip == 3
		dstrp *= 1.48;
	end

	% Suction side displcement thickness
	if Itrip == 1
		if alpha_star <= 5
			dstrs = 10^(0.0679 * alpha_star) * dstr0;
		elseif alpha_star <= 12.5
			dstrs = 0.381 * 10^(0.1516 * alpha_star) * dstr0;
		else
			dstrs = 14.296 * 10^(0.0258 * alpha_star) * dstr0;
		end
	else
		if alpha_star <= 7.5
			dstrs = 10^(0.0679 * alpha_star) * dstr0;
		elseif alpha_star <= 12.5
			dstrs = 0.0162 * 10^(0.3066 * alpha_star) * dstr0;
		else
			dstrs = 52.42 * 10^(0.0258 * alpha_star) * dstr0;
		end
	end
end