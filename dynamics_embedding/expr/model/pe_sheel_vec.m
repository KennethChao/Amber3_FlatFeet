function [x_pe_sheel_vec] = pe_sheel_vec(x)
x_pe_sheel_vec=[(-0.2032).*cos(x(3)) + x(1),(-0.2032).*sin(x(3)) + x(2)];
