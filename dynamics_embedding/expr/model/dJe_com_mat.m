function [x_dJe_com_mat] = dJe_com_mat(x)
x_dJe_com_mat=[0.,0.,1.4193960469820091e-11.*(9.968927475e9.*cos(x(3)).*x(12) +  ...
  3.452578425e9.*sin(x(3)).*x(12) + 2.5658269992e10.*sin(x(3) + x(4) ...
  ).*(x(12) + x(13)) + 2.1731369032e10.*sin(x(3) + x(4) + x(5)).*(x( ...
  12) + x(13) + x(14)) + 6.54037659e9.*sin(x(3) + x(4) + x(5) + x(6) ...
  ).*(x(12) + x(13) + x(14) + x(15)) + (-6.900526968e9).*sin(x(3) +  ...
  x(4) + x(5) + x(6) + (-1.).*x(7)).*(x(12) + x(13) + x(14) + x(15)  ...
  + (-1.).*x(16)) + (-2.973626008e9).*sin(x(3) + x(4) + x(5) + x(6)  ...
  + (-1.).*x(7) + (-1.).*x(8)).*(x(12) + x(13) + x(14) + x(15) + ( ...
  -1.).*x(16) + (-1.).*x(17)) + (-3.29502525e8).*cos(x(3) + x(4) +  ...
  x(5) + x(6) + (-1.).*x(7) + (-1.).*x(8) + (-1.).*x(9)).*(x(12) +  ...
  x(13) + x(14) + x(15) + (-1.).*x(16) + (-1.).*x(17) + (-1.).*x(18) ...
  ) + (-1.26408575e8).*sin(x(3) + x(4) + x(5) + x(6) + (-1.).*x(7) + ...
   (-1.).*x(8) + (-1.).*x(9)).*(x(12) + x(13) + x(14) + x(15) + ( ...
  -1.).*x(16) + (-1.).*x(17) + (-1.).*x(18))), ...
  1.4193960469820091e-11.*(2.5658269992e10.*sin(x(3) + x(4)).*(x(12) ...
   + x(13)) + 2.1731369032e10.*sin(x(3) + x(4) + x(5)).*(x(12) + x( ...
  13) + x(14)) + 6.54037659e9.*sin(x(3) + x(4) + x(5) + x(6)).*(x( ...
  12) + x(13) + x(14) + x(15)) + (-6.900526968e9).*sin(x(3) + x(4) + ...
   x(5) + x(6) + (-1.).*x(7)).*(x(12) + x(13) + x(14) + x(15) + ( ...
  -1.).*x(16)) + (-2.973626008e9).*sin(x(3) + x(4) + x(5) + x(6) + ( ...
  -1.).*x(7) + (-1.).*x(8)).*(x(12) + x(13) + x(14) + x(15) + (-1.) ...
  .*x(16) + (-1.).*x(17)) + (-3.29502525e8).*cos(x(3) + x(4) + x(5)  ...
  + x(6) + (-1.).*x(7) + (-1.).*x(8) + (-1.).*x(9)).*(x(12) + x(13)  ...
  + x(14) + x(15) + (-1.).*x(16) + (-1.).*x(17) + (-1.).*x(18)) + ( ...
  -1.26408575e8).*sin(x(3) + x(4) + x(5) + x(6) + (-1.).*x(7) + ( ...
  -1.).*x(8) + (-1.).*x(9)).*(x(12) + x(13) + x(14) + x(15) + (-1.) ...
  .*x(16) + (-1.).*x(17) + (-1.).*x(18))),1.4193960469820091e-11.*( ...
  2.1731369032e10.*sin(x(3) + x(4) + x(5)).*(x(12) + x(13) + x(14))  ...
  + 6.54037659e9.*sin(x(3) + x(4) + x(5) + x(6)).*(x(12) + x(13) +  ...
  x(14) + x(15)) + (-6.900526968e9).*sin(x(3) + x(4) + x(5) + x(6) + ...
   (-1.).*x(7)).*(x(12) + x(13) + x(14) + x(15) + (-1.).*x(16)) + ( ...
  -2.973626008e9).*sin(x(3) + x(4) + x(5) + x(6) + (-1.).*x(7) + ( ...
  -1.).*x(8)).*(x(12) + x(13) + x(14) + x(15) + (-1.).*x(16) + (-1.) ...
  .*x(17)) + (-3.29502525e8).*cos(x(3) + x(4) + x(5) + x(6) + (-1.) ...
  .*x(7) + (-1.).*x(8) + (-1.).*x(9)).*(x(12) + x(13) + x(14) + x( ...
  15) + (-1.).*x(16) + (-1.).*x(17) + (-1.).*x(18)) + ( ...
  -1.26408575e8).*sin(x(3) + x(4) + x(5) + x(6) + (-1.).*x(7) + ( ...
  -1.).*x(8) + (-1.).*x(9)).*(x(12) + x(13) + x(14) + x(15) + (-1.) ...
  .*x(16) + (-1.).*x(17) + (-1.).*x(18))),1.4193960469820091e-11.*( ...
  6.54037659e9.*sin(x(3) + x(4) + x(5) + x(6)).*(x(12) + x(13) + x( ...
  14) + x(15)) + (-6.900526968e9).*sin(x(3) + x(4) + x(5) + x(6) + ( ...
  -1.).*x(7)).*(x(12) + x(13) + x(14) + x(15) + (-1.).*x(16)) + ( ...
  -2.973626008e9).*sin(x(3) + x(4) + x(5) + x(6) + (-1.).*x(7) + ( ...
  -1.).*x(8)).*(x(12) + x(13) + x(14) + x(15) + (-1.).*x(16) + (-1.) ...
  .*x(17)) + (-3.29502525e8).*cos(x(3) + x(4) + x(5) + x(6) + (-1.) ...
  .*x(7) + (-1.).*x(8) + (-1.).*x(9)).*(x(12) + x(13) + x(14) + x( ...
  15) + (-1.).*x(16) + (-1.).*x(17) + (-1.).*x(18)) + ( ...
  -1.26408575e8).*sin(x(3) + x(4) + x(5) + x(6) + (-1.).*x(7) + ( ...
  -1.).*x(8) + (-1.).*x(9)).*(x(12) + x(13) + x(14) + x(15) + (-1.) ...
  .*x(16) + (-1.).*x(17) + (-1.).*x(18))),1.4193960469820091e-11.*( ...
  6.900526968e9.*sin(x(3) + x(4) + x(5) + x(6) + (-1.).*x(7)).*(x( ...
  12) + x(13) + x(14) + x(15) + (-1.).*x(16)) + 2.973626008e9.*sin( ...
  x(3) + x(4) + x(5) + x(6) + (-1.).*x(7) + (-1.).*x(8)).*(x(12) +  ...
  x(13) + x(14) + x(15) + (-1.).*x(16) + (-1.).*x(17)) +  ...
  3.29502525e8.*cos(x(3) + x(4) + x(5) + x(6) + (-1.).*x(7) + (-1.) ...
  .*x(8) + (-1.).*x(9)).*(x(12) + x(13) + x(14) + x(15) + (-1.).*x( ...
  16) + (-1.).*x(17) + (-1.).*x(18)) + 1.26408575e8.*sin(x(3) + x(4) ...
   + x(5) + x(6) + (-1.).*x(7) + (-1.).*x(8) + (-1.).*x(9)).*(x(12)  ...
  + x(13) + x(14) + x(15) + (-1.).*x(16) + (-1.).*x(17) + (-1.).*x( ...
  18))),1.4193960469820091e-11.*(2.973626008e9.*sin(x(3) + x(4) + x( ...
  5) + x(6) + (-1.).*x(7) + (-1.).*x(8)).*(x(12) + x(13) + x(14) +  ...
  x(15) + (-1.).*x(16) + (-1.).*x(17)) + 3.29502525e8.*cos(x(3) + x( ...
  4) + x(5) + x(6) + (-1.).*x(7) + (-1.).*x(8) + (-1.).*x(9)).*(x( ...
  12) + x(13) + x(14) + x(15) + (-1.).*x(16) + (-1.).*x(17) + (-1.) ...
  .*x(18)) + 1.26408575e8.*sin(x(3) + x(4) + x(5) + x(6) + (-1.).*x( ...
  7) + (-1.).*x(8) + (-1.).*x(9)).*(x(12) + x(13) + x(14) + x(15) +  ...
  (-1.).*x(16) + (-1.).*x(17) + (-1.).*x(18))), ...
  1.4193960469820091e-11.*(3.29502525e8.*cos(x(3) + x(4) + x(5) + x( ...
  6) + (-1.).*x(7) + (-1.).*x(8) + (-1.).*x(9)).*(x(12) + x(13) + x( ...
  14) + x(15) + (-1.).*x(16) + (-1.).*x(17) + (-1.).*x(18)) +  ...
  1.26408575e8.*sin(x(3) + x(4) + x(5) + x(6) + (-1.).*x(7) + (-1.) ...
  .*x(8) + (-1.).*x(9)).*(x(12) + x(13) + x(14) + x(15) + (-1.).*x( ...
  16) + (-1.).*x(17) + (-1.).*x(18)));0.,0.,1.4193960469820091e-11.* ...
  ((-3.452578425e9).*cos(x(3)).*x(12) + 9.968927475e9.*sin(x(3)).*x( ...
  12) + (-2.5658269992e10).*cos(x(3) + x(4)).*(x(12) + x(13)) + ( ...
  -2.1731369032e10).*cos(x(3) + x(4) + x(5)).*(x(12) + x(13) + x(14) ...
  ) + (-6.54037659e9).*cos(x(3) + x(4) + x(5) + x(6)).*(x(12) + x( ...
  13) + x(14) + x(15)) + 6.900526968e9.*cos(x(3) + x(4) + x(5) + x( ...
  6) + (-1.).*x(7)).*(x(12) + x(13) + x(14) + x(15) + (-1.).*x(16))  ...
  + 2.973626008e9.*cos(x(3) + x(4) + x(5) + x(6) + (-1.).*x(7) + ( ...
  -1.).*x(8)).*(x(12) + x(13) + x(14) + x(15) + (-1.).*x(16) + (-1.) ...
  .*x(17)) + 1.26408575e8.*cos(x(3) + x(4) + x(5) + x(6) + (-1.).*x( ...
  7) + (-1.).*x(8) + (-1.).*x(9)).*(x(12) + x(13) + x(14) + x(15) +  ...
  (-1.).*x(16) + (-1.).*x(17) + (-1.).*x(18)) + (-3.29502525e8).* ...
  sin(x(3) + x(4) + x(5) + x(6) + (-1.).*x(7) + (-1.).*x(8) + (-1.) ...
  .*x(9)).*(x(12) + x(13) + x(14) + x(15) + (-1.).*x(16) + (-1.).*x( ...
  17) + (-1.).*x(18))),1.4193960469820091e-11.*((-2.5658269992e10).* ...
  cos(x(3) + x(4)).*(x(12) + x(13)) + (-2.1731369032e10).*cos(x(3) + ...
   x(4) + x(5)).*(x(12) + x(13) + x(14)) + (-6.54037659e9).*cos(x(3) ...
   + x(4) + x(5) + x(6)).*(x(12) + x(13) + x(14) + x(15)) +  ...
  6.900526968e9.*cos(x(3) + x(4) + x(5) + x(6) + (-1.).*x(7)).*(x( ...
  12) + x(13) + x(14) + x(15) + (-1.).*x(16)) + 2.973626008e9.*cos( ...
  x(3) + x(4) + x(5) + x(6) + (-1.).*x(7) + (-1.).*x(8)).*(x(12) +  ...
  x(13) + x(14) + x(15) + (-1.).*x(16) + (-1.).*x(17)) +  ...
  1.26408575e8.*cos(x(3) + x(4) + x(5) + x(6) + (-1.).*x(7) + (-1.) ...
  .*x(8) + (-1.).*x(9)).*(x(12) + x(13) + x(14) + x(15) + (-1.).*x( ...
  16) + (-1.).*x(17) + (-1.).*x(18)) + (-3.29502525e8).*sin(x(3) +  ...
  x(4) + x(5) + x(6) + (-1.).*x(7) + (-1.).*x(8) + (-1.).*x(9)).*(x( ...
  12) + x(13) + x(14) + x(15) + (-1.).*x(16) + (-1.).*x(17) + (-1.) ...
  .*x(18))),1.4193960469820091e-11.*((-2.1731369032e10).*cos(x(3) +  ...
  x(4) + x(5)).*(x(12) + x(13) + x(14)) + (-6.54037659e9).*cos(x(3)  ...
  + x(4) + x(5) + x(6)).*(x(12) + x(13) + x(14) + x(15)) +  ...
  6.900526968e9.*cos(x(3) + x(4) + x(5) + x(6) + (-1.).*x(7)).*(x( ...
  12) + x(13) + x(14) + x(15) + (-1.).*x(16)) + 2.973626008e9.*cos( ...
  x(3) + x(4) + x(5) + x(6) + (-1.).*x(7) + (-1.).*x(8)).*(x(12) +  ...
  x(13) + x(14) + x(15) + (-1.).*x(16) + (-1.).*x(17)) +  ...
  1.26408575e8.*cos(x(3) + x(4) + x(5) + x(6) + (-1.).*x(7) + (-1.) ...
  .*x(8) + (-1.).*x(9)).*(x(12) + x(13) + x(14) + x(15) + (-1.).*x( ...
  16) + (-1.).*x(17) + (-1.).*x(18)) + (-3.29502525e8).*sin(x(3) +  ...
  x(4) + x(5) + x(6) + (-1.).*x(7) + (-1.).*x(8) + (-1.).*x(9)).*(x( ...
  12) + x(13) + x(14) + x(15) + (-1.).*x(16) + (-1.).*x(17) + (-1.) ...
  .*x(18))),1.4193960469820091e-11.*((-6.54037659e9).*cos(x(3) + x( ...
  4) + x(5) + x(6)).*(x(12) + x(13) + x(14) + x(15)) +  ...
  6.900526968e9.*cos(x(3) + x(4) + x(5) + x(6) + (-1.).*x(7)).*(x( ...
  12) + x(13) + x(14) + x(15) + (-1.).*x(16)) + 2.973626008e9.*cos( ...
  x(3) + x(4) + x(5) + x(6) + (-1.).*x(7) + (-1.).*x(8)).*(x(12) +  ...
  x(13) + x(14) + x(15) + (-1.).*x(16) + (-1.).*x(17)) +  ...
  1.26408575e8.*cos(x(3) + x(4) + x(5) + x(6) + (-1.).*x(7) + (-1.) ...
  .*x(8) + (-1.).*x(9)).*(x(12) + x(13) + x(14) + x(15) + (-1.).*x( ...
  16) + (-1.).*x(17) + (-1.).*x(18)) + (-3.29502525e8).*sin(x(3) +  ...
  x(4) + x(5) + x(6) + (-1.).*x(7) + (-1.).*x(8) + (-1.).*x(9)).*(x( ...
  12) + x(13) + x(14) + x(15) + (-1.).*x(16) + (-1.).*x(17) + (-1.) ...
  .*x(18))),1.4193960469820091e-11.*((-6.900526968e9).*cos(x(3) + x( ...
  4) + x(5) + x(6) + (-1.).*x(7)).*(x(12) + x(13) + x(14) + x(15) +  ...
  (-1.).*x(16)) + (-2.973626008e9).*cos(x(3) + x(4) + x(5) + x(6) +  ...
  (-1.).*x(7) + (-1.).*x(8)).*(x(12) + x(13) + x(14) + x(15) + (-1.) ...
  .*x(16) + (-1.).*x(17)) + (-1.26408575e8).*cos(x(3) + x(4) + x(5)  ...
  + x(6) + (-1.).*x(7) + (-1.).*x(8) + (-1.).*x(9)).*(x(12) + x(13)  ...
  + x(14) + x(15) + (-1.).*x(16) + (-1.).*x(17) + (-1.).*x(18)) +  ...
  3.29502525e8.*sin(x(3) + x(4) + x(5) + x(6) + (-1.).*x(7) + (-1.) ...
  .*x(8) + (-1.).*x(9)).*(x(12) + x(13) + x(14) + x(15) + (-1.).*x( ...
  16) + (-1.).*x(17) + (-1.).*x(18))),1.4193960469820091e-11.*(( ...
  -2.973626008e9).*cos(x(3) + x(4) + x(5) + x(6) + (-1.).*x(7) + ( ...
  -1.).*x(8)).*(x(12) + x(13) + x(14) + x(15) + (-1.).*x(16) + (-1.) ...
  .*x(17)) + (-1.26408575e8).*cos(x(3) + x(4) + x(5) + x(6) + (-1.) ...
  .*x(7) + (-1.).*x(8) + (-1.).*x(9)).*(x(12) + x(13) + x(14) + x( ...
  15) + (-1.).*x(16) + (-1.).*x(17) + (-1.).*x(18)) + 3.29502525e8.* ...
  sin(x(3) + x(4) + x(5) + x(6) + (-1.).*x(7) + (-1.).*x(8) + (-1.) ...
  .*x(9)).*(x(12) + x(13) + x(14) + x(15) + (-1.).*x(16) + (-1.).*x( ...
  17) + (-1.).*x(18))),1.4193960469820091e-11.*((-1.26408575e8).* ...
  cos(x(3) + x(4) + x(5) + x(6) + (-1.).*x(7) + (-1.).*x(8) + (-1.) ...
  .*x(9)).*(x(12) + x(13) + x(14) + x(15) + (-1.).*x(16) + (-1.).*x( ...
  17) + (-1.).*x(18)) + 3.29502525e8.*sin(x(3) + x(4) + x(5) + x(6)  ...
  + (-1.).*x(7) + (-1.).*x(8) + (-1.).*x(9)).*(x(12) + x(13) + x(14) ...
   + x(15) + (-1.).*x(16) + (-1.).*x(17) + (-1.).*x(18)))];