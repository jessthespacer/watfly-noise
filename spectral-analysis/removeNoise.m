function signal = removeNoise(combined, noise, ignoreLowSNR=false)
	% removeNoise 	Remove noise from a signal
	%	signal = removeNoise(combined, noise) removes noise from combined,
	%	given both in dB
	% 
	% Combined > noise. Otherwise, something is wrong.
	%
	% If combined < 5 dB louder than noise (SNR too low)
	% Requires both frequency vectors to be identical.

	if !ignoreLowSNR && any(combined - noise < 5)
		error("error: SNR too low. Combined signal and noise should be at least 5 dB greater than noise. To ignore this error, set ignoreLowSNR to true.");
	end

	if combined(:, 1) ~= noise(:, 1)
		error("error: Frequency vectors do not match.")
	end

	f = combined(:, 1);
	combinedAmp = combined(:, 2);
	noiseAmp = noise(:, 2);
	signalAmp = zeros(size(combinedAmp));

	% Remove noise level from signal
	valid = combinedAmp > noiseAmp;
	invalid = not(valid);
	signalAmp(valid) = ...
		10 .* log10(10.^(combinedAmp(valid) ./ 10) - 10.^(noiseAmp(valid) ./ 10));

	signal = [f signalAmp];
	signal(invalid, :) = 0;
end