function spectrum = readSpec(f, source)
	% readSpec 	Read spectrum to matrix from FFT
	%	spectrum = readSpec(f, source) reads FFT from f and stores data
	%	in spectrum
	%
	%	Set source to 'aud' to read Audacity FFTs
	%	Set source to 'rew' to read REW FFTs
	%
	%	spectrum(:, 1) returns frequency vector
	%	spectrum(:, 2) returns amplitude vector

	if source == 'aud'
		spectrum = dlmread(f, '\t', 1, 0);
	elseif source == 'rew'
		spectrum = dlmread(f, ',', 14, 0);
end