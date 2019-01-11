% Subroutine DIRECTL
function Dbar_l = direct_l(M, Theta, Phi)
	% Computes low frequency directivity function for a given observer
	Mc = 0.8 * M;
	ThetaR = Theta * pi / 180;
	PhiR = Phi * pi / 180;

	Dbar_l = (sin(ThetaR) * sin(PhiR))^2 / (1 + M * cos(ThetaR))^4;
end