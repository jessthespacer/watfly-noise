% Subroutine blunt
function G5 = G5comp(hdstar, eta)
	if hdstar < 0.25
		mu = 0.1211;
	elseif hdstar <= 0.62
		mu = -0.2175 * hdstar + 0.1755;
	elseif hdstar < 1.15
		mu = -0.0308 * hdstar + 0.0596;
	else
		mu = 0.0242;
	end

	if hdstar <= 0.02
		M = 0.0;
	elseif hdstar < 0.5
		M = 68.724 * hdstar - 1.35;
	elseif hdstar <= 0.62
		M = 308.475 * hdstar - 121.23;
	elseif hdstar <= 1.15
		M = 224.811 * hdstar - 69.354;
	elseif hdstar < 1.2
		M = 1583.28 * hdstar - 1631.592;
	else
		M = 268.344;
	end

	if M < 0
		M = 0;
	end

	eta0 = sqrt((M * M * mu^4) / (6.25 + M * M * mu * mu));

	K = 2.5 * sqrt(1 - (eta0 / mu)^2) - 2.5 - M * eta0;

	if eta <= eta0
		G5 = M * eta + K;
	elseif eta > eta0 && eta <= 0
		G5 = 2.6 * sqrt(1 - (eta / mu)^2) - 2.5;
	elseif eta > 0 && eta <= 0.03616
		G5 = sqrt(1.5625 - 1194.99 * eta^2) - 1.25;
	else
		G5 = -155.543 * eta + 4.375;
	end
end