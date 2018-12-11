function spectrum = readSpec(f)
	% readSpec 	Read spectrum to matrix from Audacity FFT
	%	spectrum = readSpec(f) reads Audacity FFT from f and stores data
	%	in spectrum
	%
	%	spectrum(:, 1) returns frequency vector
	%	spectrum(:, 2) returns amplitude vector

	spectrum = dlmread(f, '\t', 1, 0);
end