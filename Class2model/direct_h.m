% Subroutine DIRECTH
function Dbar_h = direct_h(M, Theta, Phi)
	% Computes high frequency directivity function for a given observer
	Mc = 0.8 * M;
	ThetaR = Theta * pi / 180;
	PhiR = Phi * pi / 180;

	Dbar_h = 2 * (sin(ThetaR / 2))^2 * (sin(PhiR))^2 / ...
		((1 + M * cos(ThetaR)) * (1 + (M - Mc) * cos(ThetaR))^2);
end