function [x_Dya2_mat] = Dya2_mat(x)
x_Dya2_mat=[(27030/28181),0,(1/70452500000).*((-3306401425).*cos(x(3)) + ( ...
  -24488853992).*cos(x(3) + x(4)) + (-20561953032).*cos(x(3) + x(4)  ...
  + x(5)) + (-6540376590).*cos(x(3) + x(4) + x(5) + x(6)) +  ...
  7853801920.*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7)) +  ...
  1804210008.*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x( ...
  8)) + (-2122690952).*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) +  ...
  (-1).*x(8) + (-1).*x(9)) + (-329502525).*sin(x(3))),( ...
  1/35226250000).*((-12244426996).*cos(x(3) + x(4)) + (-10280976516) ...
  .*cos(x(3) + x(4) + x(5)) + (-3270188295).*cos(x(3) + x(4) + x(5)  ...
  + x(6)) + 3926900960.*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7))  ...
  + 902105004.*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x( ...
  8)) + (-1061345476).*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) +  ...
  (-1).*x(8) + (-1).*x(9))),(1/35226250000).*((-10280976516).*cos(x( ...
  3) + x(4) + x(5)) + (-3270188295).*cos(x(3) + x(4) + x(5) + x(6))  ...
  + 773012.*(5080.*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7)) +  ...
  1167.*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8)) + ( ...
  -1373).*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8) +  ...
  (-1).*x(9)))),(1/35226250000).*((-3270188295).*cos(x(3) + x(4) +  ...
  x(5) + x(6)) + 773012.*(5080.*cos(x(3) + x(4) + x(5) + x(6) + (-1) ...
  .*x(7)) + 1167.*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1) ...
  .*x(8)) + (-1373).*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + ( ...
  -1).*x(8) + (-1).*x(9)))),(193253/8806562500).*((-5080).*cos(x(3)  ...
  + x(4) + x(5) + x(6) + (-1).*x(7)) + (-1167).*cos(x(3) + x(4) + x( ...
  5) + x(6) + (-1).*x(7) + (-1).*x(8)) + 1373.*cos(x(3) + x(4) + x( ...
  5) + x(6) + (-1).*x(7) + (-1).*x(8) + (-1).*x(9))),( ...
  193253/8806562500).*((-1167).*cos(x(3) + x(4) + x(5) + x(6) + (-1) ...
  .*x(7) + (-1).*x(8)) + 1373.*cos(x(3) + x(4) + x(5) + x(6) + (-1) ...
  .*x(7) + (-1).*x(8) + (-1).*x(9))),(265336369/8806562500).*cos(x( ...
  3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8) + (-1).*x(9)),0, ...
  0,0,0,0,0,0,0,0;0,(27030/28181),(1/70452500000).*(329502525.*cos( ...
  x(3)) + (-3306401425).*sin(x(3)) + (-24488853992).*sin(x(3) + x(4) ...
  ) + (-20561953032).*sin(x(3) + x(4) + x(5)) + (-6540376590).*sin( ...
  x(3) + x(4) + x(5) + x(6)) + 7853801920.*sin(x(3) + x(4) + x(5) +  ...
  x(6) + (-1).*x(7)) + 1804210008.*sin(x(3) + x(4) + x(5) + x(6) + ( ...
  -1).*x(7) + (-1).*x(8)) + (-2122690952).*sin(x(3) + x(4) + x(5) +  ...
  x(6) + (-1).*x(7) + (-1).*x(8) + (-1).*x(9))),(1/35226250000).*(( ...
  -12244426996).*sin(x(3) + x(4)) + (-10280976516).*sin(x(3) + x(4)  ...
  + x(5)) + (-3270188295).*sin(x(3) + x(4) + x(5) + x(6)) +  ...
  3926900960.*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7)) +  ...
  902105004.*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8) ...
  ) + (-1061345476).*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + ( ...
  -1).*x(8) + (-1).*x(9))),(1/35226250000).*((-10280976516).*sin(x( ...
  3) + x(4) + x(5)) + (-3270188295).*sin(x(3) + x(4) + x(5) + x(6))  ...
  + 773012.*(5080.*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7)) +  ...
  1167.*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8)) + ( ...
  -1373).*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8) +  ...
  (-1).*x(9)))),(1/35226250000).*((-3270188295).*sin(x(3) + x(4) +  ...
  x(5) + x(6)) + 773012.*(5080.*sin(x(3) + x(4) + x(5) + x(6) + (-1) ...
  .*x(7)) + 1167.*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1) ...
  .*x(8)) + (-1373).*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + ( ...
  -1).*x(8) + (-1).*x(9)))),(-193253/8806562500).*(5080.*sin(x(3) +  ...
  x(4) + x(5) + x(6) + (-1).*x(7)) + 1167.*sin(x(3) + x(4) + x(5) +  ...
  x(6) + (-1).*x(7) + (-1).*x(8)) + (-1373).*sin(x(3) + x(4) + x(5)  ...
  + x(6) + (-1).*x(7) + (-1).*x(8) + (-1).*x(9))),( ...
  -193253/8806562500).*(1167.*sin(x(3) + x(4) + x(5) + x(6) + (-1).* ...
  x(7) + (-1).*x(8)) + (-1373).*sin(x(3) + x(4) + x(5) + x(6) + (-1) ...
  .*x(7) + (-1).*x(8) + (-1).*x(9))),(265336369/8806562500).*sin(x( ...
  3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8) + (-1).*x(9)),0, ...
  0,0,0,0,0,0,0,0;0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,( ...
  -254/625).*(cos(x(4)) + cos(x(4) + x(5)) + (-2).*cos(x(4) + x(5) + ...
   x(6) + (-1).*x(7) + (-1/2).*x(8)).*cos((1/2).*x(8))),(-254/625).* ...
  (cos(x(4) + x(5)) + (-2).*cos(x(4) + x(5) + x(6) + (-1).*x(7) + ( ...
  -1/2).*x(8)).*cos((1/2).*x(8))),(508/625).*cos(x(4) + x(5) + x(6)  ...
  + (-1).*x(7) + (-1/2).*x(8)).*cos((1/2).*x(8)),(-508/625).*cos(x( ...
  4) + x(5) + x(6) + (-1).*x(7) + (-1/2).*x(8)).*cos((1/2).*x(8)),( ...
  -254/625).*cos(x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8)),0,0, ...
  0,0,0,0,0,0,0,0;0,0,0,(-254/625).*(sin(x(4)) + sin(x(4) + x(5)) +  ...
  (-2).*cos((1/2).*x(8)).*sin(x(4) + x(5) + x(6) + (-1).*x(7) + ( ...
  -1/2).*x(8))),(-254/625).*(sin(x(4) + x(5)) + (-2).*cos((1/2).*x( ...
  8)).*sin(x(4) + x(5) + x(6) + (-1).*x(7) + (-1/2).*x(8))),( ...
  508/625).*cos((1/2).*x(8)).*sin(x(4) + x(5) + x(6) + (-1).*x(7) +  ...
  (-1/2).*x(8)),(-508/625).*cos((1/2).*x(8)).*sin(x(4) + x(5) + x(6) ...
   + (-1).*x(7) + (-1/2).*x(8)),(-254/625).*sin(x(4) + x(5) + x(6) + ...
   (-1).*x(7) + (-1).*x(8)),0,0,0,0,0,0,0,0,0,0;0,0,0,1,1,1,(-1),( ...
  -1),(-1),0,0,0,0,0,0,0,0,0];