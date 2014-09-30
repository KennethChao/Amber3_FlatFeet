function [ output_args ] = FootTrajGen( step_height, stpe_length)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% Foot Traj X
xstart=0;
PolyOrder = 5,

A=zeros(1,PolyOrder+1);

T=[ 1 0 0 0 0 0 ;
    0 1 0 0 0 0 ;
    0 0 2 0 0 0 ;
    1 1 1 1 1 1 ;
    0 1 2 3 4 5 ;
    0 0 2 6 12 20 ;    
  ];

A=T^-1*[xstart;0;0;stpe_length;0;0];
tPoly=0;

syms x t real
a = sym('a', [1,PolyOrder+1]);

for i=1:PolyOrder+1
    if i==1
    x = a(i);
    tPoly = t;
    else
    x = x+a(i)*tPoly;
    tPoly = tPoly*t;
    end
    
end

% x=0; (a(1)+a(2)*(t)+a(3)*(t)^2+a(4)*(t)^3+a(5)*(t)^4+a(6)*(t)^5);
for i=1:PolyOrder+1
 x=subs(x,a(i),A(i));   
end
% x=subs(x,a,A)
x = simplify(x);
dxdt = simplify(diff(x,t));
ddxdt = simplify(diff(dxdt,t));

matlabFunction(x, 'file', 'FootTrajGenX');
matlabFunction(dxdt, 'file', 'FootTrajGenXder');
matlabFunction(ddxdt, 'file', 'FootTrajGenXder2');
% % test
% i=1;
% for t = linspace(0,1,100)
% xt(i)=FootTragGenXTest(t);
% i=i+1;
% end
% plot(diff(xt))
% % test
% 
% % test
% hold on
% i=1;
% for t = linspace(0,1,100)
% xt(i)=FootTragGenXTestder(t)/100;
% i=i+1;
% end
% plot(xt,'r')
% % test
%% Foot Traj Y
xstart=0;
PolyOrder = 6,

A=zeros(1,PolyOrder+1);

T=[ 1 0 0 0 0 0 0;
    0 1 0 0 0 0 0;
    0 0 2 0 0 0 0;
    1 1 1 1 1 1 1;
    0 1 2 3 4 5 6;
    0 0 2 6 12 20 30;    
    1 0.5 0.5^2 0.5^3 0.5^4 0.5^5 0.5^6;
  ];

A=T^-1*[xstart;0;0;xstart;0;0;step_height];
tPoly=0;

syms y real
a = sym('a', [1,PolyOrder+1]);

for i=1:PolyOrder+1
    if i==1
    y = a(i);
    tPoly = t;
    else
    y = y+a(i)*tPoly;
    tPoly = tPoly*t;
    end
    
end

% x=0; (a(1)+a(2)*(t)+a(3)*(t)^2+a(4)*(t)^3+a(5)*(t)^4+a(6)*(t)^5);
for i=1:PolyOrder+1
 y=subs(y,a(i),A(i));   
end
% x=subs(x,a,A)
y=simplify(y);
dydt = simplify(diff(y,t));
ddydt = simplify(diff(dydt,t));

matlabFunction(y, 'file', 'FootTrajGenZ');
matlabFunction(dydt, 'file', 'FootTrajGenZder');
matlabFunction(ddydt, 'file', 'FootTrajGenZder2');




end

