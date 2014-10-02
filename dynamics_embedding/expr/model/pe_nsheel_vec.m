function [x_pe_nsheel_vec] = pe_nsheel_vec(x)
x_pe_nsheel_vec=[0.0004.*((-127.).*(3..*cos(x(3)) + cos(x(3) + x(4) + x(5) + x(6)  ...
  + (-1.).*x(7) + (-1.).*x(8) + (-1.).*x(9)) + sin(x(3)) + 8..*sin( ...
  x(3) + x(4)) + 8..*sin(x(3) + x(4) + x(5)) + (-8.).*sin(x(3) + x( ...
  4) + x(5) + x(6) + (-1.).*x(7)) + (-8.).*sin(x(3) + x(4) + x(5) +  ...
  x(6) + (-1.).*x(7) + (-1.).*x(8)) + (-1.).*sin(x(3) + x(4) + x(5)  ...
  + x(6) + (-1.).*x(7) + (-1.).*x(8) + (-1.).*x(9))) + 2500..*x(1)), ...
  0.0004.*((-127.).*((-1.).*cos(x(3)) + (-8.).*cos(x(3) + x(4)) + ( ...
  -8.).*cos(x(3) + x(4) + x(5)) + 8..*cos(x(3) + x(4) + x(5) + x(6)  ...
  + (-1.).*x(7)) + 8..*cos(x(3) + x(4) + x(5) + x(6) + (-1.).*x(7) + ...
   (-1.).*x(8)) + cos(x(3) + x(4) + x(5) + x(6) + (-1.).*x(7) + ( ...
  -1.).*x(8) + (-1.).*x(9)) + 3..*sin(x(3)) + sin(x(3) + x(4) + x(5) ...
   + x(6) + (-1.).*x(7) + (-1.).*x(8) + (-1.).*x(9))) + 2500..*x(2)) ...
  ];
