% Subroutine DIRECTH
function Dbar_h = direct_h(M, Theta, Phi)
	% Computes high frequency directivity function for a given observer
	degrad = 0.017453;

	Mc = 0.8 * M;
	ThetaR = Theta * degrad;
	PhiR = Phi * degrad;

	Dbar_h = 2 * sin(ThetaR / 2)^2 * sin(PhiR)^2 / ...
		((1 + M * cos(ThetaR)) * (1 + (M - Mc) * cos(ThetaR))^2);
end