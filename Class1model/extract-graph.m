clear all; close all; clc;

format long;

% To get line eqn. from y = m log10(x) + b graphs
p1 = [10.256668331442935, -11.786864645472527];
p2 = [0.05436930797862856, 8.657323819565827];

m = (p2(2) - p1(2))/(log10(p2(1)) - log10(p1(1)));
b = p1(2) - m * log10(p1(1));

if b >= 0
	sgn = "+";
else
	sgn = "-";
	b = -b;
end

disp(["y = " mat2str(m) " * log10(D) " sgn " " mat2str(b)])