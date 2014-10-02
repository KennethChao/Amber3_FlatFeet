function [x_hdot_nstoe_sca] = hdot_nstoe_sca(x)
x_hdot_nstoe_sca=0.0004.*(2500..*x(11) + 127..*((-3.).*cos(x(3)).*x(12) + (-1.).* ...
  sin(x(3)).*x(12) + (-8.).*sin(x(3) + x(4)).*(x(12) + x(13)) + ( ...
  -8.).*sin(x(3) + x(4) + x(5)).*(x(12) + x(13) + x(14)) + 8..*sin( ...
  x(3) + x(4) + x(5) + x(6) + (-1.).*x(7)).*(x(12) + x(13) + x(14) + ...
   x(15) + (-1.).*x(16)) + 8..*sin(x(3) + x(4) + x(5) + x(6) + (-1.) ...
  .*x(7) + (-1.).*x(8)).*(x(12) + x(13) + x(14) + x(15) + (-1.).*x( ...
  16) + (-1.).*x(17)) + 3..*cos(x(3) + x(4) + x(5) + x(6) + (-1.).* ...
  x(7) + (-1.).*x(8) + (-1.).*x(9)).*(x(12) + x(13) + x(14) + x(15)  ...
  + (-1.).*x(16) + (-1.).*x(17) + (-1.).*x(18)) + sin(x(3) + x(4) +  ...
  x(5) + x(6) + (-1.).*x(7) + (-1.).*x(8) + (-1.).*x(9)).*(x(12) +  ...
  x(13) + x(14) + x(15) + (-1.).*x(16) + (-1.).*x(17) + (-1.).*x(18) ...
  )));