function [x_DLfya2_mat] = DLfya2_mat(x)
x_DLfya2_mat=[0,0,(1/70452500000).*(((-329502525).*cos(x(3)) + 3306401425.*sin( ...
  x(3)) + 24488853992.*sin(x(3) + x(4)) + 20561953032.*sin(x(3) + x( ...
  4) + x(5)) + 6540376590.*sin(x(3) + x(4) + x(5) + x(6)) + ( ...
  -7853801920).*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7)) + ( ...
  -1804210008).*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).* ...
  x(8)) + 2122690952.*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + ( ...
  -1).*x(8) + (-1).*x(9))).*x(12) + 2.*(12244426996.*sin(x(3) + x(4) ...
  ) + 10280976516.*sin(x(3) + x(4) + x(5)) + 3270188295.*sin(x(3) +  ...
  x(4) + x(5) + x(6)) + (-3926900960).*sin(x(3) + x(4) + x(5) + x(6) ...
   + (-1).*x(7)) + (-902105004).*sin(x(3) + x(4) + x(5) + x(6) + ( ...
  -1).*x(7) + (-1).*x(8)) + 1061345476.*sin(x(3) + x(4) + x(5) + x( ...
  6) + (-1).*x(7) + (-1).*x(8) + (-1).*x(9))).*x(13) + 2.*( ...
  10280976516.*sin(x(3) + x(4) + x(5)) + 3270188295.*sin(x(3) + x(4) ...
   + x(5) + x(6)) + (-773012).*(5080.*sin(x(3) + x(4) + x(5) + x(6)  ...
  + (-1).*x(7)) + 1167.*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + ...
   (-1).*x(8)) + (-1373).*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) ...
   + (-1).*x(8) + (-1).*x(9)))).*x(14) + 2.*(3270188295.*sin(x(3) +  ...
  x(4) + x(5) + x(6)) + (-773012).*(5080.*sin(x(3) + x(4) + x(5) +  ...
  x(6) + (-1).*x(7)) + 1167.*sin(x(3) + x(4) + x(5) + x(6) + (-1).* ...
  x(7) + (-1).*x(8)) + (-1373).*sin(x(3) + x(4) + x(5) + x(6) + (-1) ...
  .*x(7) + (-1).*x(8) + (-1).*x(9)))).*x(15) + 1546024.*(5080.*sin( ...
  x(3) + x(4) + x(5) + x(6) + (-1).*x(7)) + 1167.*sin(x(3) + x(4) +  ...
  x(5) + x(6) + (-1).*x(7) + (-1).*x(8)) + (-1373).*sin(x(3) + x(4)  ...
  + x(5) + x(6) + (-1).*x(7) + (-1).*x(8) + (-1).*x(9))).*x(16) +  ...
  1546024.*(1167.*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1) ...
  .*x(8)) + (-1373).*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + ( ...
  -1).*x(8) + (-1).*x(9))).*x(17) + (-2122690952).*sin(x(3) + x(4) + ...
   x(5) + x(6) + (-1).*x(7) + (-1).*x(8) + (-1).*x(9)).*x(18)),( ...
  1/35226250000).*((12244426996.*sin(x(3) + x(4)) + 10280976516.* ...
  sin(x(3) + x(4) + x(5)) + 3270188295.*sin(x(3) + x(4) + x(5) + x( ...
  6)) + (-3926900960).*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7)) + ...
   (-902105004).*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).* ...
  x(8)) + 1061345476.*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + ( ...
  -1).*x(8) + (-1).*x(9))).*x(12) + (12244426996.*sin(x(3) + x(4)) + ...
   10280976516.*sin(x(3) + x(4) + x(5)) + 3270188295.*sin(x(3) + x( ...
  4) + x(5) + x(6)) + (-3926900960).*sin(x(3) + x(4) + x(5) + x(6) + ...
   (-1).*x(7)) + (-902105004).*sin(x(3) + x(4) + x(5) + x(6) + (-1) ...
  .*x(7) + (-1).*x(8)) + 1061345476.*sin(x(3) + x(4) + x(5) + x(6) + ...
   (-1).*x(7) + (-1).*x(8) + (-1).*x(9))).*x(13) + (10280976516.* ...
  sin(x(3) + x(4) + x(5)) + 3270188295.*sin(x(3) + x(4) + x(5) + x( ...
  6)) + (-773012).*(5080.*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) ...
  ) + 1167.*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8)) ...
   + (-1373).*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x( ...
  8) + (-1).*x(9)))).*x(14) + (3270188295.*sin(x(3) + x(4) + x(5) +  ...
  x(6)) + (-773012).*(5080.*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x( ...
  7)) + 1167.*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x( ...
  8)) + (-1373).*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).* ...
  x(8) + (-1).*x(9)))).*x(15) + 773012.*(5080.*sin(x(3) + x(4) + x( ...
  5) + x(6) + (-1).*x(7)) + 1167.*sin(x(3) + x(4) + x(5) + x(6) + ( ...
  -1).*x(7) + (-1).*x(8)) + (-1373).*sin(x(3) + x(4) + x(5) + x(6) + ...
   (-1).*x(7) + (-1).*x(8) + (-1).*x(9))).*x(16) + 773012.*(1167.* ...
  sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8)) + (-1373) ...
  .*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8) + (-1).* ...
  x(9))).*x(17) + (-1061345476).*sin(x(3) + x(4) + x(5) + x(6) + ( ...
  -1).*x(7) + (-1).*x(8) + (-1).*x(9)).*x(18)),(1/35226250000).*(( ...
  10280976516.*sin(x(3) + x(4) + x(5)) + 3270188295.*sin(x(3) + x(4) ...
   + x(5) + x(6)) + (-773012).*(5080.*sin(x(3) + x(4) + x(5) + x(6)  ...
  + (-1).*x(7)) + 1167.*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + ...
   (-1).*x(8)) + (-1373).*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) ...
   + (-1).*x(8) + (-1).*x(9)))).*x(12) + (10280976516.*sin(x(3) + x( ...
  4) + x(5)) + 3270188295.*sin(x(3) + x(4) + x(5) + x(6)) + ( ...
  -773012).*(5080.*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7)) +  ...
  1167.*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8)) + ( ...
  -1373).*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8) +  ...
  (-1).*x(9)))).*x(13) + (10280976516.*sin(x(3) + x(4) + x(5)) +  ...
  3270188295.*sin(x(3) + x(4) + x(5) + x(6)) + (-773012).*(5080.* ...
  sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7)) + 1167.*sin(x(3) + x( ...
  4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8)) + (-1373).*sin(x(3) +  ...
  x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8) + (-1).*x(9)))).*x( ...
  14) + (3270188295.*sin(x(3) + x(4) + x(5) + x(6)) + (-773012).*( ...
  5080.*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7)) + 1167.*sin(x(3) ...
   + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8)) + (-1373).*sin(x( ...
  3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8) + (-1).*x(9)))) ...
  .*x(15) + 773012.*(5080.*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x( ...
  7)) + 1167.*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x( ...
  8)) + (-1373).*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).* ...
  x(8) + (-1).*x(9))).*x(16) + 773012.*(1167.*sin(x(3) + x(4) + x(5) ...
   + x(6) + (-1).*x(7) + (-1).*x(8)) + (-1373).*sin(x(3) + x(4) + x( ...
  5) + x(6) + (-1).*x(7) + (-1).*x(8) + (-1).*x(9))).*x(17) + ( ...
  -1061345476).*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).* ...
  x(8) + (-1).*x(9)).*x(18)),(1/35226250000).*((3270188295.*sin(x(3) ...
   + x(4) + x(5) + x(6)) + (-773012).*(5080.*sin(x(3) + x(4) + x(5)  ...
  + x(6) + (-1).*x(7)) + 1167.*sin(x(3) + x(4) + x(5) + x(6) + (-1) ...
  .*x(7) + (-1).*x(8)) + (-1373).*sin(x(3) + x(4) + x(5) + x(6) + ( ...
  -1).*x(7) + (-1).*x(8) + (-1).*x(9)))).*x(12) + (3270188295.*sin( ...
  x(3) + x(4) + x(5) + x(6)) + (-773012).*(5080.*sin(x(3) + x(4) +  ...
  x(5) + x(6) + (-1).*x(7)) + 1167.*sin(x(3) + x(4) + x(5) + x(6) +  ...
  (-1).*x(7) + (-1).*x(8)) + (-1373).*sin(x(3) + x(4) + x(5) + x(6)  ...
  + (-1).*x(7) + (-1).*x(8) + (-1).*x(9)))).*x(13) + (3270188295.* ...
  sin(x(3) + x(4) + x(5) + x(6)) + (-773012).*(5080.*sin(x(3) + x(4) ...
   + x(5) + x(6) + (-1).*x(7)) + 1167.*sin(x(3) + x(4) + x(5) + x(6) ...
   + (-1).*x(7) + (-1).*x(8)) + (-1373).*sin(x(3) + x(4) + x(5) + x( ...
  6) + (-1).*x(7) + (-1).*x(8) + (-1).*x(9)))).*x(14) + ( ...
  3270188295.*sin(x(3) + x(4) + x(5) + x(6)) + (-773012).*(5080.* ...
  sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7)) + 1167.*sin(x(3) + x( ...
  4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8)) + (-1373).*sin(x(3) +  ...
  x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8) + (-1).*x(9)))).*x( ...
  15) + 773012.*(5080.*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7)) + ...
   1167.*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8)) +  ...
  (-1373).*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8) + ...
   (-1).*x(9))).*x(16) + 773012.*(1167.*sin(x(3) + x(4) + x(5) + x( ...
  6) + (-1).*x(7) + (-1).*x(8)) + (-1373).*sin(x(3) + x(4) + x(5) +  ...
  x(6) + (-1).*x(7) + (-1).*x(8) + (-1).*x(9))).*x(17) + ( ...
  -1061345476).*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).* ...
  x(8) + (-1).*x(9)).*x(18)),(193253/8806562500).*((5080.*sin(x(3) + ...
   x(4) + x(5) + x(6) + (-1).*x(7)) + 1167.*sin(x(3) + x(4) + x(5) + ...
   x(6) + (-1).*x(7) + (-1).*x(8)) + (-1373).*sin(x(3) + x(4) + x(5) ...
   + x(6) + (-1).*x(7) + (-1).*x(8) + (-1).*x(9))).*x(12) + (5080.* ...
  sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7)) + 1167.*sin(x(3) + x( ...
  4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8)) + (-1373).*sin(x(3) +  ...
  x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8) + (-1).*x(9))).*x(13) ...
   + 5080.*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7)).*x(14) +  ...
  1167.*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8)).*x( ...
  14) + (-1373).*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).* ...
  x(8) + (-1).*x(9)).*x(14) + 5080.*sin(x(3) + x(4) + x(5) + x(6) +  ...
  (-1).*x(7)).*x(15) + 1167.*sin(x(3) + x(4) + x(5) + x(6) + (-1).* ...
  x(7) + (-1).*x(8)).*x(15) + (-1373).*sin(x(3) + x(4) + x(5) + x(6) ...
   + (-1).*x(7) + (-1).*x(8) + (-1).*x(9)).*x(15) + (-5080).*sin(x( ...
  3) + x(4) + x(5) + x(6) + (-1).*x(7)).*x(16) + (-1167).*sin(x(3) + ...
   x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8)).*x(16) + 1373.*sin( ...
  x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8) + (-1).*x(9)) ...
  .*x(16) + (-1167).*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + ( ...
  -1).*x(8)).*x(17) + 1373.*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x( ...
  7) + (-1).*x(8) + (-1).*x(9)).*x(17) + 1373.*sin(x(3) + x(4) + x( ...
  5) + x(6) + (-1).*x(7) + (-1).*x(8) + (-1).*x(9)).*x(18)),( ...
  193253/8806562500).*((1167.*sin(x(3) + x(4) + x(5) + x(6) + (-1).* ...
  x(7) + (-1).*x(8)) + (-1373).*sin(x(3) + x(4) + x(5) + x(6) + (-1) ...
  .*x(7) + (-1).*x(8) + (-1).*x(9))).*x(12) + (1167.*sin(x(3) + x(4) ...
   + x(5) + x(6) + (-1).*x(7) + (-1).*x(8)) + (-1373).*sin(x(3) + x( ...
  4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8) + (-1).*x(9))).*x(13) + ...
   1167.*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8)).* ...
  x(14) + (-1373).*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1) ...
  .*x(8) + (-1).*x(9)).*x(14) + 1167.*sin(x(3) + x(4) + x(5) + x(6)  ...
  + (-1).*x(7) + (-1).*x(8)).*x(15) + (-1373).*sin(x(3) + x(4) + x( ...
  5) + x(6) + (-1).*x(7) + (-1).*x(8) + (-1).*x(9)).*x(15) + (-1167) ...
  .*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8)).*x(16)  ...
  + 1373.*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8) +  ...
  (-1).*x(9)).*x(16) + (-1167).*sin(x(3) + x(4) + x(5) + x(6) + (-1) ...
  .*x(7) + (-1).*x(8)).*x(17) + 1373.*sin(x(3) + x(4) + x(5) + x(6)  ...
  + (-1).*x(7) + (-1).*x(8) + (-1).*x(9)).*x(17) + 1373.*sin(x(3) +  ...
  x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8) + (-1).*x(9)).*x(18)) ...
  ,(-265336369/8806562500).*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x( ...
  7) + (-1).*x(8) + (-1).*x(9)).*(x(12) + x(13) + x(14) + x(15) + ( ...
  -1).*x(16) + (-1).*x(17) + (-1).*x(18)),(27030/28181),0,( ...
  1/70452500000).*((-3306401425).*cos(x(3)) + (-24488853992).*cos(x( ...
  3) + x(4)) + (-20561953032).*cos(x(3) + x(4) + x(5)) + ( ...
  -6540376590).*cos(x(3) + x(4) + x(5) + x(6)) + 7853801920.*cos(x( ...
  3) + x(4) + x(5) + x(6) + (-1).*x(7)) + 1804210008.*cos(x(3) + x( ...
  4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8)) + (-2122690952).*cos( ...
  x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8) + (-1).*x(9))  ...
  + (-329502525).*sin(x(3))),(1/35226250000).*((-12244426996).*cos( ...
  x(3) + x(4)) + (-10280976516).*cos(x(3) + x(4) + x(5)) + ( ...
  -3270188295).*cos(x(3) + x(4) + x(5) + x(6)) + 3926900960.*cos(x( ...
  3) + x(4) + x(5) + x(6) + (-1).*x(7)) + 902105004.*cos(x(3) + x(4) ...
   + x(5) + x(6) + (-1).*x(7) + (-1).*x(8)) + (-1061345476).*cos(x( ...
  3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8) + (-1).*x(9))),( ...
  1/35226250000).*((-10280976516).*cos(x(3) + x(4) + x(5)) + ( ...
  -3270188295).*cos(x(3) + x(4) + x(5) + x(6)) + 773012.*(5080.*cos( ...
  x(3) + x(4) + x(5) + x(6) + (-1).*x(7)) + 1167.*cos(x(3) + x(4) +  ...
  x(5) + x(6) + (-1).*x(7) + (-1).*x(8)) + (-1373).*cos(x(3) + x(4)  ...
  + x(5) + x(6) + (-1).*x(7) + (-1).*x(8) + (-1).*x(9)))),( ...
  1/35226250000).*((-3270188295).*cos(x(3) + x(4) + x(5) + x(6)) +  ...
  773012.*(5080.*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7)) +  ...
  1167.*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8)) + ( ...
  -1373).*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8) +  ...
  (-1).*x(9)))),(193253/8806562500).*((-5080).*cos(x(3) + x(4) + x( ...
  5) + x(6) + (-1).*x(7)) + (-1167).*cos(x(3) + x(4) + x(5) + x(6) + ...
   (-1).*x(7) + (-1).*x(8)) + 1373.*cos(x(3) + x(4) + x(5) + x(6) +  ...
  (-1).*x(7) + (-1).*x(8) + (-1).*x(9))),(193253/8806562500).*(( ...
  -1167).*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8)) + ...
   1373.*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8) + ( ...
  -1).*x(9))),(265336369/8806562500).*cos(x(3) + x(4) + x(5) + x(6)  ...
  + (-1).*x(7) + (-1).*x(8) + (-1).*x(9));0,0,(1/70452500000).*((-1) ...
  .*(3306401425.*cos(x(3)) + 24488853992.*cos(x(3) + x(4)) +  ...
  20561953032.*cos(x(3) + x(4) + x(5)) + 6540376590.*cos(x(3) + x(4) ...
   + x(5) + x(6)) + (-7853801920).*cos(x(3) + x(4) + x(5) + x(6) + ( ...
  -1).*x(7)) + (-1804210008).*cos(x(3) + x(4) + x(5) + x(6) + (-1).* ...
  x(7) + (-1).*x(8)) + 2122690952.*cos(x(3) + x(4) + x(5) + x(6) + ( ...
  -1).*x(7) + (-1).*x(8) + (-1).*x(9)) + 329502525.*sin(x(3))).*x( ...
  12) + (-2).*(12244426996.*cos(x(3) + x(4)) + 10280976516.*cos(x(3) ...
   + x(4) + x(5)) + 3270188295.*cos(x(3) + x(4) + x(5) + x(6)) + ( ...
  -3926900960).*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7)) + ( ...
  -902105004).*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x( ...
  8)) + 1061345476.*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + ( ...
  -1).*x(8) + (-1).*x(9))).*x(13) + 2.*((-10280976516).*cos(x(3) +  ...
  x(4) + x(5)) + (-3270188295).*cos(x(3) + x(4) + x(5) + x(6)) +  ...
  773012.*(5080.*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7)) +  ...
  1167.*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8)) + ( ...
  -1373).*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8) +  ...
  (-1).*x(9)))).*x(14) + 2.*((-3270188295).*cos(x(3) + x(4) + x(5) + ...
   x(6)) + 773012.*(5080.*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) ...
  ) + 1167.*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8)) ...
   + (-1373).*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x( ...
  8) + (-1).*x(9)))).*x(15) + 1546024.*((-5080).*cos(x(3) + x(4) +  ...
  x(5) + x(6) + (-1).*x(7)) + (-1167).*cos(x(3) + x(4) + x(5) + x(6) ...
   + (-1).*x(7) + (-1).*x(8)) + 1373.*cos(x(3) + x(4) + x(5) + x(6)  ...
  + (-1).*x(7) + (-1).*x(8) + (-1).*x(9))).*x(16) + 1546024.*(( ...
  -1167).*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8)) + ...
   1373.*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8) + ( ...
  -1).*x(9))).*x(17) + 2122690952.*cos(x(3) + x(4) + x(5) + x(6) + ( ...
  -1).*x(7) + (-1).*x(8) + (-1).*x(9)).*x(18)),(1/35226250000).*(( ...
  -1).*(12244426996.*cos(x(3) + x(4)) + 10280976516.*cos(x(3) + x(4) ...
   + x(5)) + 3270188295.*cos(x(3) + x(4) + x(5) + x(6)) + ( ...
  -3926900960).*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7)) + ( ...
  -902105004).*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x( ...
  8)) + 1061345476.*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + ( ...
  -1).*x(8) + (-1).*x(9))).*x(12) + (-1).*(12244426996.*cos(x(3) +  ...
  x(4)) + 10280976516.*cos(x(3) + x(4) + x(5)) + 3270188295.*cos(x( ...
  3) + x(4) + x(5) + x(6)) + (-3926900960).*cos(x(3) + x(4) + x(5) + ...
   x(6) + (-1).*x(7)) + (-902105004).*cos(x(3) + x(4) + x(5) + x(6)  ...
  + (-1).*x(7) + (-1).*x(8)) + 1061345476.*cos(x(3) + x(4) + x(5) +  ...
  x(6) + (-1).*x(7) + (-1).*x(8) + (-1).*x(9))).*x(13) + (( ...
  -10280976516).*cos(x(3) + x(4) + x(5)) + (-3270188295).*cos(x(3) + ...
   x(4) + x(5) + x(6)) + 773012.*(5080.*cos(x(3) + x(4) + x(5) + x( ...
  6) + (-1).*x(7)) + 1167.*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x( ...
  7) + (-1).*x(8)) + (-1373).*cos(x(3) + x(4) + x(5) + x(6) + (-1).* ...
  x(7) + (-1).*x(8) + (-1).*x(9)))).*x(14) + ((-3270188295).*cos(x( ...
  3) + x(4) + x(5) + x(6)) + 773012.*(5080.*cos(x(3) + x(4) + x(5) + ...
   x(6) + (-1).*x(7)) + 1167.*cos(x(3) + x(4) + x(5) + x(6) + (-1).* ...
  x(7) + (-1).*x(8)) + (-1373).*cos(x(3) + x(4) + x(5) + x(6) + (-1) ...
  .*x(7) + (-1).*x(8) + (-1).*x(9)))).*x(15) + 773012.*((-5080).* ...
  cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7)) + (-1167).*cos(x(3) +  ...
  x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8)) + 1373.*cos(x(3) +  ...
  x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8) + (-1).*x(9))).*x(16) ...
   + 773012.*((-1167).*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) +  ...
  (-1).*x(8)) + 1373.*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + ( ...
  -1).*x(8) + (-1).*x(9))).*x(17) + 1061345476.*cos(x(3) + x(4) + x( ...
  5) + x(6) + (-1).*x(7) + (-1).*x(8) + (-1).*x(9)).*x(18)),( ...
  1/35226250000).*((-1).*(10280976516.*cos(x(3) + x(4) + x(5)) +  ...
  3270188295.*cos(x(3) + x(4) + x(5) + x(6)) + 773012.*((-5080).* ...
  cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7)) + (-1167).*cos(x(3) +  ...
  x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8)) + 1373.*cos(x(3) +  ...
  x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8) + (-1).*x(9)))).*x( ...
  12) + (-1).*(10280976516.*cos(x(3) + x(4) + x(5)) + 3270188295.* ...
  cos(x(3) + x(4) + x(5) + x(6)) + 773012.*((-5080).*cos(x(3) + x(4) ...
   + x(5) + x(6) + (-1).*x(7)) + (-1167).*cos(x(3) + x(4) + x(5) +  ...
  x(6) + (-1).*x(7) + (-1).*x(8)) + 1373.*cos(x(3) + x(4) + x(5) +  ...
  x(6) + (-1).*x(7) + (-1).*x(8) + (-1).*x(9)))).*x(13) + (( ...
  -10280976516).*cos(x(3) + x(4) + x(5)) + (-3270188295).*cos(x(3) + ...
   x(4) + x(5) + x(6)) + 773012.*(5080.*cos(x(3) + x(4) + x(5) + x( ...
  6) + (-1).*x(7)) + 1167.*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x( ...
  7) + (-1).*x(8)) + (-1373).*cos(x(3) + x(4) + x(5) + x(6) + (-1).* ...
  x(7) + (-1).*x(8) + (-1).*x(9)))).*x(14) + ((-3270188295).*cos(x( ...
  3) + x(4) + x(5) + x(6)) + 773012.*(5080.*cos(x(3) + x(4) + x(5) + ...
   x(6) + (-1).*x(7)) + 1167.*cos(x(3) + x(4) + x(5) + x(6) + (-1).* ...
  x(7) + (-1).*x(8)) + (-1373).*cos(x(3) + x(4) + x(5) + x(6) + (-1) ...
  .*x(7) + (-1).*x(8) + (-1).*x(9)))).*x(15) + 773012.*((-5080).* ...
  cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7)) + (-1167).*cos(x(3) +  ...
  x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8)) + 1373.*cos(x(3) +  ...
  x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8) + (-1).*x(9))).*x(16) ...
   + 773012.*((-1167).*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) +  ...
  (-1).*x(8)) + 1373.*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + ( ...
  -1).*x(8) + (-1).*x(9))).*x(17) + 1061345476.*cos(x(3) + x(4) + x( ...
  5) + x(6) + (-1).*x(7) + (-1).*x(8) + (-1).*x(9)).*x(18)),( ...
  1/70452500000).*(((-6540376590).*cos(x(3) + x(4) + x(5) + x(6)) +  ...
  1546024.*(5080.*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7)) +  ...
  1167.*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8)) + ( ...
  -1373).*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8) +  ...
  (-1).*x(9)))).*x(12) + 2.*((-3270188295).*cos(x(3) + x(4) + x(5) + ...
   x(6)) + 773012.*(5080.*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) ...
  ) + 1167.*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8)) ...
   + (-1373).*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x( ...
  8) + (-1).*x(9)))).*x(13) + 2.*((-3270188295).*cos(x(3) + x(4) +  ...
  x(5) + x(6)) + 773012.*(5080.*cos(x(3) + x(4) + x(5) + x(6) + (-1) ...
  .*x(7)) + 1167.*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1) ...
  .*x(8)) + (-1373).*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + ( ...
  -1).*x(8) + (-1).*x(9)))).*x(14) + 2.*((-3270188295).*cos(x(3) +  ...
  x(4) + x(5) + x(6)) + 773012.*(5080.*cos(x(3) + x(4) + x(5) + x(6) ...
   + (-1).*x(7)) + 1167.*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7)  ...
  + (-1).*x(8)) + (-1373).*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x( ...
  7) + (-1).*x(8) + (-1).*x(9)))).*x(15) + 1546024.*((-5080).*cos(x( ...
  3) + x(4) + x(5) + x(6) + (-1).*x(7)) + (-1167).*cos(x(3) + x(4) + ...
   x(5) + x(6) + (-1).*x(7) + (-1).*x(8)) + 1373.*cos(x(3) + x(4) +  ...
  x(5) + x(6) + (-1).*x(7) + (-1).*x(8) + (-1).*x(9))).*x(16) +  ...
  1546024.*((-1167).*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + ( ...
  -1).*x(8)) + 1373.*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + ( ...
  -1).*x(8) + (-1).*x(9))).*x(17) + 2122690952.*cos(x(3) + x(4) + x( ...
  5) + x(6) + (-1).*x(7) + (-1).*x(8) + (-1).*x(9)).*x(18)),( ...
  -193253/8806562500).*((5080.*cos(x(3) + x(4) + x(5) + x(6) + (-1) ...
  .*x(7)) + 1167.*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1) ...
  .*x(8)) + (-1373).*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + ( ...
  -1).*x(8) + (-1).*x(9))).*x(12) + (5080.*cos(x(3) + x(4) + x(5) +  ...
  x(6) + (-1).*x(7)) + 1167.*cos(x(3) + x(4) + x(5) + x(6) + (-1).* ...
  x(7) + (-1).*x(8)) + (-1373).*cos(x(3) + x(4) + x(5) + x(6) + (-1) ...
  .*x(7) + (-1).*x(8) + (-1).*x(9))).*x(13) + 5080.*cos(x(3) + x(4)  ...
  + x(5) + x(6) + (-1).*x(7)).*x(14) + 1167.*cos(x(3) + x(4) + x(5)  ...
  + x(6) + (-1).*x(7) + (-1).*x(8)).*x(14) + (-1373).*cos(x(3) + x( ...
  4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8) + (-1).*x(9)).*x(14) +  ...
  5080.*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7)).*x(15) + 1167.* ...
  cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8)).*x(15) +  ...
  (-1373).*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8) + ...
   (-1).*x(9)).*x(15) + (-5080).*cos(x(3) + x(4) + x(5) + x(6) + ( ...
  -1).*x(7)).*x(16) + (-1167).*cos(x(3) + x(4) + x(5) + x(6) + (-1) ...
  .*x(7) + (-1).*x(8)).*x(16) + 1373.*cos(x(3) + x(4) + x(5) + x(6)  ...
  + (-1).*x(7) + (-1).*x(8) + (-1).*x(9)).*x(16) + (-1167).*cos(x(3) ...
   + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8)).*x(17) + 1373.* ...
  cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8) + (-1).*x( ...
  9)).*x(17) + 1373.*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + ( ...
  -1).*x(8) + (-1).*x(9)).*x(18)),(193253/8806562500).*(((-1167).* ...
  cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8)) + 1373.* ...
  cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8) + (-1).*x( ...
  9))).*x(12) + ((-1167).*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) ...
   + (-1).*x(8)) + 1373.*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7)  ...
  + (-1).*x(8) + (-1).*x(9))).*x(13) + (-1167).*cos(x(3) + x(4) + x( ...
  5) + x(6) + (-1).*x(7) + (-1).*x(8)).*x(14) + 1373.*cos(x(3) + x( ...
  4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8) + (-1).*x(9)).*x(14) +  ...
  (-1167).*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8)) ...
  .*x(15) + 1373.*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1) ...
  .*x(8) + (-1).*x(9)).*x(15) + 1167.*cos(x(3) + x(4) + x(5) + x(6)  ...
  + (-1).*x(7) + (-1).*x(8)).*x(16) + (-1373).*cos(x(3) + x(4) + x( ...
  5) + x(6) + (-1).*x(7) + (-1).*x(8) + (-1).*x(9)).*x(16) + 1167.* ...
  cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8)).*x(17) +  ...
  (-1373).*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8) + ...
   (-1).*x(9)).*x(17) + (-1373).*cos(x(3) + x(4) + x(5) + x(6) + ( ...
  -1).*x(7) + (-1).*x(8) + (-1).*x(9)).*x(18)),( ...
  265336369/8806562500).*cos(x(3) + x(4) + x(5) + x(6) + (-1).*x(7)  ...
  + (-1).*x(8) + (-1).*x(9)).*(x(12) + x(13) + x(14) + x(15) + (-1) ...
  .*x(16) + (-1).*x(17) + (-1).*x(18)),0,(27030/28181),( ...
  1/70452500000).*(329502525.*cos(x(3)) + (-3306401425).*sin(x(3)) + ...
   (-24488853992).*sin(x(3) + x(4)) + (-20561953032).*sin(x(3) + x( ...
  4) + x(5)) + (-6540376590).*sin(x(3) + x(4) + x(5) + x(6)) +  ...
  7853801920.*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7)) +  ...
  1804210008.*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x( ...
  8)) + (-2122690952).*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) +  ...
  (-1).*x(8) + (-1).*x(9))),(1/35226250000).*((-12244426996).*sin(x( ...
  3) + x(4)) + (-10280976516).*sin(x(3) + x(4) + x(5)) + ( ...
  -3270188295).*sin(x(3) + x(4) + x(5) + x(6)) + 3926900960.*sin(x( ...
  3) + x(4) + x(5) + x(6) + (-1).*x(7)) + 902105004.*sin(x(3) + x(4) ...
   + x(5) + x(6) + (-1).*x(7) + (-1).*x(8)) + (-1061345476).*sin(x( ...
  3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8) + (-1).*x(9))),( ...
  1/35226250000).*((-10280976516).*sin(x(3) + x(4) + x(5)) + ( ...
  -3270188295).*sin(x(3) + x(4) + x(5) + x(6)) + 773012.*(5080.*sin( ...
  x(3) + x(4) + x(5) + x(6) + (-1).*x(7)) + 1167.*sin(x(3) + x(4) +  ...
  x(5) + x(6) + (-1).*x(7) + (-1).*x(8)) + (-1373).*sin(x(3) + x(4)  ...
  + x(5) + x(6) + (-1).*x(7) + (-1).*x(8) + (-1).*x(9)))),( ...
  1/35226250000).*((-3270188295).*sin(x(3) + x(4) + x(5) + x(6)) +  ...
  773012.*(5080.*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7)) +  ...
  1167.*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8)) + ( ...
  -1373).*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8) +  ...
  (-1).*x(9)))),(-193253/8806562500).*(5080.*sin(x(3) + x(4) + x(5)  ...
  + x(6) + (-1).*x(7)) + 1167.*sin(x(3) + x(4) + x(5) + x(6) + (-1) ...
  .*x(7) + (-1).*x(8)) + (-1373).*sin(x(3) + x(4) + x(5) + x(6) + ( ...
  -1).*x(7) + (-1).*x(8) + (-1).*x(9))),(-193253/8806562500).*( ...
  1167.*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8)) + ( ...
  -1373).*sin(x(3) + x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8) +  ...
  (-1).*x(9))),(265336369/8806562500).*sin(x(3) + x(4) + x(5) + x(6) ...
   + (-1).*x(7) + (-1).*x(8) + (-1).*x(9));0,0,0,0,0,0,0,0,0,0,0,0, ...
  1,1,1,0,0,0;0,0,0,(254/625).*((sin(x(4)) + sin(x(4) + x(5)) + (-2) ...
  .*cos((1/2).*x(8)).*sin(x(4) + x(5) + x(6) + (-1).*x(7) + (-1/2).* ...
  x(8))).*x(13) + (sin(x(4) + x(5)) + (-2).*cos((1/2).*x(8)).*sin(x( ...
  4) + x(5) + x(6) + (-1).*x(7) + (-1/2).*x(8))).*x(14) + (-2).*cos( ...
  (1/2).*x(8)).*sin(x(4) + x(5) + x(6) + (-1).*x(7) + (-1/2).*x(8)) ...
  .*x(15) + 2.*cos((1/2).*x(8)).*sin(x(4) + x(5) + x(6) + (-1).*x(7) ...
   + (-1/2).*x(8)).*x(16) + sin(x(4) + x(5) + x(6) + (-1).*x(7) + ( ...
  -1).*x(8)).*x(17)),(254/625).*((sin(x(4) + x(5)) + (-2).*cos((1/2) ...
  .*x(8)).*sin(x(4) + x(5) + x(6) + (-1).*x(7) + (-1/2).*x(8))).*x( ...
  13) + (sin(x(4) + x(5)) + (-2).*cos((1/2).*x(8)).*sin(x(4) + x(5)  ...
  + x(6) + (-1).*x(7) + (-1/2).*x(8))).*x(14) + (-2).*cos((1/2).*x( ...
  8)).*sin(x(4) + x(5) + x(6) + (-1).*x(7) + (-1/2).*x(8)).*x(15) +  ...
  2.*cos((1/2).*x(8)).*sin(x(4) + x(5) + x(6) + (-1).*x(7) + (-1/2) ...
  .*x(8)).*x(16) + sin(x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8)) ...
  .*x(17)),(-254/625).*(2.*cos((1/2).*x(8)).*sin(x(4) + x(5) + x(6)  ...
  + (-1).*x(7) + (-1/2).*x(8)).*x(13) + 2.*cos((1/2).*x(8)).*sin(x( ...
  4) + x(5) + x(6) + (-1).*x(7) + (-1/2).*x(8)).*x(14) + 2.*cos(( ...
  1/2).*x(8)).*sin(x(4) + x(5) + x(6) + (-1).*x(7) + (-1/2).*x(8)).* ...
  x(15) + (-2).*cos((1/2).*x(8)).*sin(x(4) + x(5) + x(6) + (-1).*x( ...
  7) + (-1/2).*x(8)).*x(16) + (-1).*sin(x(4) + x(5) + x(6) + (-1).* ...
  x(7) + (-1).*x(8)).*x(17)),(254/625).*(2.*cos((1/2).*x(8)).*sin(x( ...
  4) + x(5) + x(6) + (-1).*x(7) + (-1/2).*x(8)).*x(13) + 2.*cos(( ...
  1/2).*x(8)).*sin(x(4) + x(5) + x(6) + (-1).*x(7) + (-1/2).*x(8)).* ...
  x(14) + 2.*cos((1/2).*x(8)).*sin(x(4) + x(5) + x(6) + (-1).*x(7) + ...
   (-1/2).*x(8)).*x(15) + (-2).*cos((1/2).*x(8)).*sin(x(4) + x(5) +  ...
  x(6) + (-1).*x(7) + (-1/2).*x(8)).*x(16) + (-1).*sin(x(4) + x(5) + ...
   x(6) + (-1).*x(7) + (-1).*x(8)).*x(17)),(254/625).*sin(x(4) + x( ...
  5) + x(6) + (-1).*x(7) + (-1).*x(8)).*(x(13) + x(14) + x(15) + ( ...
  -1).*x(16) + (-1).*x(17)),0,0,0,0,(-254/625).*(cos(x(4)) + cos(x( ...
  4) + x(5)) + (-2).*cos(x(4) + x(5) + x(6) + (-1).*x(7) + (-1/2).* ...
  x(8)).*cos((1/2).*x(8))),(-254/625).*(cos(x(4) + x(5)) + (-2).* ...
  cos(x(4) + x(5) + x(6) + (-1).*x(7) + (-1/2).*x(8)).*cos((1/2).*x( ...
  8))),(508/625).*cos(x(4) + x(5) + x(6) + (-1).*x(7) + (-1/2).*x(8) ...
  ).*cos((1/2).*x(8)),(-508/625).*cos(x(4) + x(5) + x(6) + (-1).*x( ...
  7) + (-1/2).*x(8)).*cos((1/2).*x(8)),(-254/625).*cos(x(4) + x(5) + ...
   x(6) + (-1).*x(7) + (-1).*x(8)),0;0,0,0,(-254/625).*((cos(x(4)) + ...
   cos(x(4) + x(5)) + (-2).*cos(x(4) + x(5) + x(6) + (-1).*x(7) + ( ...
  -1/2).*x(8)).*cos((1/2).*x(8))).*x(13) + (cos(x(4) + x(5)) + (-2) ...
  .*cos(x(4) + x(5) + x(6) + (-1).*x(7) + (-1/2).*x(8)).*cos((1/2).* ...
  x(8))).*x(14) + (-2).*cos(x(4) + x(5) + x(6) + (-1).*x(7) + (-1/2) ...
  .*x(8)).*cos((1/2).*x(8)).*x(15) + 2.*cos(x(4) + x(5) + x(6) + ( ...
  -1).*x(7) + (-1/2).*x(8)).*cos((1/2).*x(8)).*x(16) + cos(x(4) + x( ...
  5) + x(6) + (-1).*x(7) + (-1).*x(8)).*x(17)),(-254/625).*((cos(x( ...
  4) + x(5)) + (-2).*cos(x(4) + x(5) + x(6) + (-1).*x(7) + (-1/2).* ...
  x(8)).*cos((1/2).*x(8))).*x(13) + (cos(x(4) + x(5)) + (-2).*cos(x( ...
  4) + x(5) + x(6) + (-1).*x(7) + (-1/2).*x(8)).*cos((1/2).*x(8))).* ...
  x(14) + (-2).*cos(x(4) + x(5) + x(6) + (-1).*x(7) + (-1/2).*x(8)) ...
  .*cos((1/2).*x(8)).*x(15) + 2.*cos(x(4) + x(5) + x(6) + (-1).*x(7) ...
   + (-1/2).*x(8)).*cos((1/2).*x(8)).*x(16) + cos(x(4) + x(5) + x(6) ...
   + (-1).*x(7) + (-1).*x(8)).*x(17)),(254/625).*(2.*cos(x(4) + x(5) ...
   + x(6) + (-1).*x(7) + (-1/2).*x(8)).*cos((1/2).*x(8)).*x(13) +  ...
  2.*cos(x(4) + x(5) + x(6) + (-1).*x(7) + (-1/2).*x(8)).*cos((1/2) ...
  .*x(8)).*x(14) + 2.*cos(x(4) + x(5) + x(6) + (-1).*x(7) + (-1/2).* ...
  x(8)).*cos((1/2).*x(8)).*x(15) + (-2).*cos(x(4) + x(5) + x(6) + ( ...
  -1).*x(7) + (-1/2).*x(8)).*cos((1/2).*x(8)).*x(16) + (-1).*cos(x( ...
  4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8)).*x(17)),(-254/625).*( ...
  2.*cos(x(4) + x(5) + x(6) + (-1).*x(7) + (-1/2).*x(8)).*cos((1/2) ...
  .*x(8)).*x(13) + 2.*cos(x(4) + x(5) + x(6) + (-1).*x(7) + (-1/2).* ...
  x(8)).*cos((1/2).*x(8)).*x(14) + 2.*cos(x(4) + x(5) + x(6) + (-1) ...
  .*x(7) + (-1/2).*x(8)).*cos((1/2).*x(8)).*x(15) + (-2).*cos(x(4) + ...
   x(5) + x(6) + (-1).*x(7) + (-1/2).*x(8)).*cos((1/2).*x(8)).*x(16) ...
   + (-1).*cos(x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8)).*x(17)) ...
  ,(-254/625).*cos(x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8)).*( ...
  x(13) + x(14) + x(15) + (-1).*x(16) + (-1).*x(17)),0,0,0,0,( ...
  -254/625).*(sin(x(4)) + sin(x(4) + x(5)) + (-2).*cos((1/2).*x(8)) ...
  .*sin(x(4) + x(5) + x(6) + (-1).*x(7) + (-1/2).*x(8))),(-254/625) ...
  .*(sin(x(4) + x(5)) + (-2).*cos((1/2).*x(8)).*sin(x(4) + x(5) + x( ...
  6) + (-1).*x(7) + (-1/2).*x(8))),(508/625).*cos((1/2).*x(8)).*sin( ...
  x(4) + x(5) + x(6) + (-1).*x(7) + (-1/2).*x(8)),(-508/625).*cos(( ...
  1/2).*x(8)).*sin(x(4) + x(5) + x(6) + (-1).*x(7) + (-1/2).*x(8)),( ...
  -254/625).*sin(x(4) + x(5) + x(6) + (-1).*x(7) + (-1).*x(8)),0;0, ...
  0,0,0,0,0,0,0,0,0,0,0,1,1,1,(-1),(-1),(-1)];
