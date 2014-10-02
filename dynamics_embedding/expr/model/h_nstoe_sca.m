function [x_h_nstoe_sca] = h_nstoe_sca(x)
x_h_nstoe_sca=0.0004.*(127..*(cos(x(3)) + 8..*cos(x(3) + x(4)) + 8..*cos(x(3) +  ...
  x(4) + x(5)) + (-8.).*cos(x(3) + x(4) + x(5) + x(6) + (-1.).*x(7)) ...
   + (-8.).*cos(x(3) + x(4) + x(5) + x(6) + (-1.).*x(7) + (-1.).*x( ...
  8)) + (-1.).*cos(x(3) + x(4) + x(5) + x(6) + (-1.).*x(7) + (-1.).* ...
  x(8) + (-1.).*x(9)) + (-3.).*sin(x(3)) + 3..*sin(x(3) + x(4) + x( ...
  5) + x(6) + (-1.).*x(7) + (-1.).*x(8) + (-1.).*x(9))) + 2500..*x( ...
  2));