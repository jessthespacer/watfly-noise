% Subroutine DIRECTL
function Dbar_l = direct_l(M, Theta, Phi)
	% Computes low frequency directivity function for a given observer
	degrad = 0.017453;

	Mc = 0.8 * M;
	ThetaR = Theta * degrad;
	PhiR = Phi * degrad;

	Dbar_l = (sin(ThetaR) * sin(PhiR))^2 / (1 + M * cos(ThetaR))^4;
end