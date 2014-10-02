function [x_dJh_nstoe_mat] = dJh_nstoe_mat(x)
x_dJh_nstoe_mat=[0.,0.,0.0508.*(3..*cos(x(3)).*x(12) + sin(x(3)).*x(12) + 8..*sin( ...
  x(3) + x(4)).*(x(12) + x(13)) + 8..*sin(x(3) + x(4) + x(5)).*(x( ...
  12) + x(13) + x(14)) + (-8.).*sin(x(3) + x(4) + x(5) + x(6) + ( ...
  -1.).*x(7)).*(x(12) + x(13) + x(14) + x(15) + (-1.).*x(16)) + ( ...
  -8.).*sin(x(3) + x(4) + x(5) + x(6) + (-1.).*x(7) + (-1.).*x(8)).* ...
  (x(12) + x(13) + x(14) + x(15) + (-1.).*x(16) + (-1.).*x(17)) + ( ...
  -3.).*cos(x(3) + x(4) + x(5) + x(6) + (-1.).*x(7) + (-1.).*x(8) +  ...
  (-1.).*x(9)).*(x(12) + x(13) + x(14) + x(15) + (-1.).*x(16) + ( ...
  -1.).*x(17) + (-1.).*x(18)) + (-1.).*sin(x(3) + x(4) + x(5) + x(6) ...
   + (-1.).*x(7) + (-1.).*x(8) + (-1.).*x(9)).*(x(12) + x(13) + x( ...
  14) + x(15) + (-1.).*x(16) + (-1.).*x(17) + (-1.).*x(18))), ...
  0.0508.*(8..*sin(x(3) + x(4)).*(x(12) + x(13)) + 8..*sin(x(3) + x( ...
  4) + x(5)).*(x(12) + x(13) + x(14)) + (-8.).*sin(x(3) + x(4) + x( ...
  5) + x(6) + (-1.).*x(7)).*(x(12) + x(13) + x(14) + x(15) + (-1.).* ...
  x(16)) + (-8.).*sin(x(3) + x(4) + x(5) + x(6) + (-1.).*x(7) + ( ...
  -1.).*x(8)).*(x(12) + x(13) + x(14) + x(15) + (-1.).*x(16) + (-1.) ...
  .*x(17)) + (-3.).*cos(x(3) + x(4) + x(5) + x(6) + (-1.).*x(7) + ( ...
  -1.).*x(8) + (-1.).*x(9)).*(x(12) + x(13) + x(14) + x(15) + (-1.) ...
  .*x(16) + (-1.).*x(17) + (-1.).*x(18)) + (-1.).*sin(x(3) + x(4) +  ...
  x(5) + x(6) + (-1.).*x(7) + (-1.).*x(8) + (-1.).*x(9)).*(x(12) +  ...
  x(13) + x(14) + x(15) + (-1.).*x(16) + (-1.).*x(17) + (-1.).*x(18) ...
  )),0.0508.*(8..*sin(x(3) + x(4) + x(5)).*(x(12) + x(13) + x(14)) + ...
   (-8.).*sin(x(3) + x(4) + x(5) + x(6) + (-1.).*x(7)).*(x(12) + x( ...
  13) + x(14) + x(15) + (-1.).*x(16)) + (-8.).*sin(x(3) + x(4) + x( ...
  5) + x(6) + (-1.).*x(7) + (-1.).*x(8)).*(x(12) + x(13) + x(14) +  ...
  x(15) + (-1.).*x(16) + (-1.).*x(17)) + (-3.).*cos(x(3) + x(4) + x( ...
  5) + x(6) + (-1.).*x(7) + (-1.).*x(8) + (-1.).*x(9)).*(x(12) + x( ...
  13) + x(14) + x(15) + (-1.).*x(16) + (-1.).*x(17) + (-1.).*x(18))  ...
  + (-1.).*sin(x(3) + x(4) + x(5) + x(6) + (-1.).*x(7) + (-1.).*x(8) ...
   + (-1.).*x(9)).*(x(12) + x(13) + x(14) + x(15) + (-1.).*x(16) + ( ...
  -1.).*x(17) + (-1.).*x(18))),0.0508.*((-8.).*sin(x(3) + x(4) + x( ...
  5) + x(6) + (-1.).*x(7)).*(x(12) + x(13) + x(14) + x(15) + (-1.).* ...
  x(16)) + (-8.).*sin(x(3) + x(4) + x(5) + x(6) + (-1.).*x(7) + ( ...
  -1.).*x(8)).*(x(12) + x(13) + x(14) + x(15) + (-1.).*x(16) + (-1.) ...
  .*x(17)) + (-3.).*cos(x(3) + x(4) + x(5) + x(6) + (-1.).*x(7) + ( ...
  -1.).*x(8) + (-1.).*x(9)).*(x(12) + x(13) + x(14) + x(15) + (-1.) ...
  .*x(16) + (-1.).*x(17) + (-1.).*x(18)) + (-1.).*sin(x(3) + x(4) +  ...
  x(5) + x(6) + (-1.).*x(7) + (-1.).*x(8) + (-1.).*x(9)).*(x(12) +  ...
  x(13) + x(14) + x(15) + (-1.).*x(16) + (-1.).*x(17) + (-1.).*x(18) ...
  )),(-0.0508).*((-8.).*sin(x(3) + x(4) + x(5) + x(6) + (-1.).*x(7)) ...
  .*(x(12) + x(13) + x(14) + x(15) + (-1.).*x(16)) + (-8.).*sin(x(3) ...
   + x(4) + x(5) + x(6) + (-1.).*x(7) + (-1.).*x(8)).*(x(12) + x(13) ...
   + x(14) + x(15) + (-1.).*x(16) + (-1.).*x(17)) + (-3.).*cos(x(3)  ...
  + x(4) + x(5) + x(6) + (-1.).*x(7) + (-1.).*x(8) + (-1.).*x(9)).*( ...
  x(12) + x(13) + x(14) + x(15) + (-1.).*x(16) + (-1.).*x(17) + ( ...
  -1.).*x(18)) + (-1.).*sin(x(3) + x(4) + x(5) + x(6) + (-1.).*x(7)  ...
  + (-1.).*x(8) + (-1.).*x(9)).*(x(12) + x(13) + x(14) + x(15) + ( ...
  -1.).*x(16) + (-1.).*x(17) + (-1.).*x(18))),(-0.0508).*((-8.).* ...
  sin(x(3) + x(4) + x(5) + x(6) + (-1.).*x(7) + (-1.).*x(8)).*(x(12) ...
   + x(13) + x(14) + x(15) + (-1.).*x(16) + (-1.).*x(17)) + (-3.).* ...
  cos(x(3) + x(4) + x(5) + x(6) + (-1.).*x(7) + (-1.).*x(8) + (-1.) ...
  .*x(9)).*(x(12) + x(13) + x(14) + x(15) + (-1.).*x(16) + (-1.).*x( ...
  17) + (-1.).*x(18)) + (-1.).*sin(x(3) + x(4) + x(5) + x(6) + (-1.) ...
  .*x(7) + (-1.).*x(8) + (-1.).*x(9)).*(x(12) + x(13) + x(14) + x( ...
  15) + (-1.).*x(16) + (-1.).*x(17) + (-1.).*x(18))),0.0508.*(3..* ...
  cos(x(3) + x(4) + x(5) + x(6) + (-1.).*x(7) + (-1.).*x(8) + (-1.) ...
  .*x(9)) + sin(x(3) + x(4) + x(5) + x(6) + (-1.).*x(7) + (-1.).*x( ...
  8) + (-1.).*x(9))).*(x(12) + x(13) + x(14) + x(15) + (-1.).*x(16)  ...
  + (-1.).*x(17) + (-1.).*x(18));0.,0.,0.0508.*((-1.).*cos(x(3)).*x( ...
  12) + 3..*sin(x(3)).*x(12) + (-8.).*cos(x(3) + x(4)).*(x(12) + x( ...
  13)) + (-8.).*cos(x(3) + x(4) + x(5)).*(x(12) + x(13) + x(14)) +  ...
  8..*cos(x(3) + x(4) + x(5) + x(6) + (-1.).*x(7)).*(x(12) + x(13) + ...
   x(14) + x(15) + (-1.).*x(16)) + 8..*cos(x(3) + x(4) + x(5) + x(6) ...
   + (-1.).*x(7) + (-1.).*x(8)).*(x(12) + x(13) + x(14) + x(15) + ( ...
  -1.).*x(16) + (-1.).*x(17)) + cos(x(3) + x(4) + x(5) + x(6) + ( ...
  -1.).*x(7) + (-1.).*x(8) + (-1.).*x(9)).*(x(12) + x(13) + x(14) +  ...
  x(15) + (-1.).*x(16) + (-1.).*x(17) + (-1.).*x(18)) + (-3.).*sin( ...
  x(3) + x(4) + x(5) + x(6) + (-1.).*x(7) + (-1.).*x(8) + (-1.).*x( ...
  9)).*(x(12) + x(13) + x(14) + x(15) + (-1.).*x(16) + (-1.).*x(17)  ...
  + (-1.).*x(18))),0.0508.*((-8.).*cos(x(3) + x(4)).*(x(12) + x(13)) ...
   + (-8.).*cos(x(3) + x(4) + x(5)).*(x(12) + x(13) + x(14)) + 8..* ...
  cos(x(3) + x(4) + x(5) + x(6) + (-1.).*x(7)).*(x(12) + x(13) + x( ...
  14) + x(15) + (-1.).*x(16)) + 8..*cos(x(3) + x(4) + x(5) + x(6) +  ...
  (-1.).*x(7) + (-1.).*x(8)).*(x(12) + x(13) + x(14) + x(15) + (-1.) ...
  .*x(16) + (-1.).*x(17)) + cos(x(3) + x(4) + x(5) + x(6) + (-1.).* ...
  x(7) + (-1.).*x(8) + (-1.).*x(9)).*(x(12) + x(13) + x(14) + x(15)  ...
  + (-1.).*x(16) + (-1.).*x(17) + (-1.).*x(18)) + (-3.).*sin(x(3) +  ...
  x(4) + x(5) + x(6) + (-1.).*x(7) + (-1.).*x(8) + (-1.).*x(9)).*(x( ...
  12) + x(13) + x(14) + x(15) + (-1.).*x(16) + (-1.).*x(17) + (-1.) ...
  .*x(18))),0.0508.*((-8.).*cos(x(3) + x(4) + x(5)).*(x(12) + x(13)  ...
  + x(14)) + 8..*cos(x(3) + x(4) + x(5) + x(6) + (-1.).*x(7)).*(x( ...
  12) + x(13) + x(14) + x(15) + (-1.).*x(16)) + 8..*cos(x(3) + x(4)  ...
  + x(5) + x(6) + (-1.).*x(7) + (-1.).*x(8)).*(x(12) + x(13) + x(14) ...
   + x(15) + (-1.).*x(16) + (-1.).*x(17)) + cos(x(3) + x(4) + x(5) + ...
   x(6) + (-1.).*x(7) + (-1.).*x(8) + (-1.).*x(9)).*(x(12) + x(13) + ...
   x(14) + x(15) + (-1.).*x(16) + (-1.).*x(17) + (-1.).*x(18)) + ( ...
  -3.).*sin(x(3) + x(4) + x(5) + x(6) + (-1.).*x(7) + (-1.).*x(8) +  ...
  (-1.).*x(9)).*(x(12) + x(13) + x(14) + x(15) + (-1.).*x(16) + ( ...
  -1.).*x(17) + (-1.).*x(18))),0.0508.*(8..*cos(x(3) + x(4) + x(5) + ...
   x(6) + (-1.).*x(7)).*(x(12) + x(13) + x(14) + x(15) + (-1.).*x( ...
  16)) + 8..*cos(x(3) + x(4) + x(5) + x(6) + (-1.).*x(7) + (-1.).*x( ...
  8)).*(x(12) + x(13) + x(14) + x(15) + (-1.).*x(16) + (-1.).*x(17)) ...
   + cos(x(3) + x(4) + x(5) + x(6) + (-1.).*x(7) + (-1.).*x(8) + ( ...
  -1.).*x(9)).*(x(12) + x(13) + x(14) + x(15) + (-1.).*x(16) + (-1.) ...
  .*x(17) + (-1.).*x(18)) + (-3.).*sin(x(3) + x(4) + x(5) + x(6) + ( ...
  -1.).*x(7) + (-1.).*x(8) + (-1.).*x(9)).*(x(12) + x(13) + x(14) +  ...
  x(15) + (-1.).*x(16) + (-1.).*x(17) + (-1.).*x(18))),(-0.0508).*( ...
  8..*cos(x(3) + x(4) + x(5) + x(6) + (-1.).*x(7)).*(x(12) + x(13) + ...
   x(14) + x(15) + (-1.).*x(16)) + 8..*cos(x(3) + x(4) + x(5) + x(6) ...
   + (-1.).*x(7) + (-1.).*x(8)).*(x(12) + x(13) + x(14) + x(15) + ( ...
  -1.).*x(16) + (-1.).*x(17)) + cos(x(3) + x(4) + x(5) + x(6) + ( ...
  -1.).*x(7) + (-1.).*x(8) + (-1.).*x(9)).*(x(12) + x(13) + x(14) +  ...
  x(15) + (-1.).*x(16) + (-1.).*x(17) + (-1.).*x(18)) + (-3.).*sin( ...
  x(3) + x(4) + x(5) + x(6) + (-1.).*x(7) + (-1.).*x(8) + (-1.).*x( ...
  9)).*(x(12) + x(13) + x(14) + x(15) + (-1.).*x(16) + (-1.).*x(17)  ...
  + (-1.).*x(18))),(-0.0508).*(8..*cos(x(3) + x(4) + x(5) + x(6) + ( ...
  -1.).*x(7) + (-1.).*x(8)).*(x(12) + x(13) + x(14) + x(15) + (-1.) ...
  .*x(16) + (-1.).*x(17)) + cos(x(3) + x(4) + x(5) + x(6) + (-1.).* ...
  x(7) + (-1.).*x(8) + (-1.).*x(9)).*(x(12) + x(13) + x(14) + x(15)  ...
  + (-1.).*x(16) + (-1.).*x(17) + (-1.).*x(18)) + (-3.).*sin(x(3) +  ...
  x(4) + x(5) + x(6) + (-1.).*x(7) + (-1.).*x(8) + (-1.).*x(9)).*(x( ...
  12) + x(13) + x(14) + x(15) + (-1.).*x(16) + (-1.).*x(17) + (-1.) ...
  .*x(18))),(-0.0508).*(cos(x(3) + x(4) + x(5) + x(6) + (-1.).*x(7)  ...
  + (-1.).*x(8) + (-1.).*x(9)) + (-3.).*sin(x(3) + x(4) + x(5) + x( ...
  6) + (-1.).*x(7) + (-1.).*x(8) + (-1.).*x(9))).*(x(12) + x(13) +  ...
  x(14) + x(15) + (-1.).*x(16) + (-1.).*x(17) + (-1.).*x(18));0.,0., ...
  0.,0.,0.,0.,0.,0.,0.];
