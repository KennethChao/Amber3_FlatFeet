function [x_y_ComX] = y_ComX(x)
x_y_ComX=1.4193960469820091e-11.*((-9.968927475e9).*cos(x(3)) +  ...
  3.29502525e8.*cos(x(3) + x(4) + x(5) + x(6) + (-1.).*x(7) + (-1.) ...
  .*x(8) + (-1.).*x(9)) + (-3.452578425e9).*sin(x(3)) + ( ...
  -2.5658269992e10).*sin(x(3) + x(4)) + (-2.1731369032e10).*sin(x(3) ...
   + x(4) + x(5)) + (-6.54037659e9).*sin(x(3) + x(4) + x(5) + x(6))  ...
  + 6.900526968e9.*sin(x(3) + x(4) + x(5) + x(6) + (-1.).*x(7)) +  ...
  2.973626008e9.*sin(x(3) + x(4) + x(5) + x(6) + (-1.).*x(7) + (-1.) ...
  .*x(8)) + 1.26408575e8.*sin(x(3) + x(4) + x(5) + x(6) + (-1.).*x( ...
  7) + (-1.).*x(8) + (-1.).*x(9)) + 7.04525e10.*x(1));
