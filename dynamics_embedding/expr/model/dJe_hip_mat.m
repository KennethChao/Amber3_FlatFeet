function [x_dJe_hip_mat] = dJe_hip_mat(x)
x_dJe_hip_mat=[0.,0.,(-0.0508).*((-3.).*cos(x(3)).*x(12) + (-1.).*sin(x(3)).*x( ...
  12) + 8..*((-1.).*sin(x(3) + x(4)).*(x(12) + x(13)) + (-1.).*sin( ...
  x(3) + x(4) + x(5)).*(x(12) + x(13) + x(14)))),(-0.4064).*((-1.).* ...
  sin(x(3) + x(4)).*(x(12) + x(13)) + (-1.).*sin(x(3) + x(4) + x(5)) ...
  .*(x(12) + x(13) + x(14))),0.4064.*sin(x(3) + x(4) + x(5)).*(x(12) ...
   + x(13) + x(14)),0.,0.,0.,0.;0.,0.,0.0508.*((-1.).*cos(x(3)).*x( ...
  12) + 3..*sin(x(3)).*x(12) + (-8.).*cos(x(3) + x(4)).*(x(12) + x( ...
  13)) + (-8.).*cos(x(3) + x(4) + x(5)).*(x(12) + x(13) + x(14))), ...
  0.0508.*((-8.).*cos(x(3) + x(4)).*(x(12) + x(13)) + (-8.).*cos(x( ...
  3) + x(4) + x(5)).*(x(12) + x(13) + x(14))),(-0.4064).*cos(x(3) +  ...
  x(4) + x(5)).*(x(12) + x(13) + x(14)),0.,0.,0.,0.];
