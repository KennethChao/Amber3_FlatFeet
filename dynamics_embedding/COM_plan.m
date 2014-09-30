function [ output_args ] = COM_plan( cog_h, hstep_leng,x0,domain_type, domains,ndomains)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


if (nargin<2)
  % Half step lenth 
  close all
  hstep_leng = 0.1;  
  cog_h = 0.9;
  %Inital condition: COM Position, Velocity and ZMP position/COM
  %acceleration
  x0=[0;0;0];
  domain_type = 2;
  domain_time = 2.5;
  Sstime = 1.5;%domains{1,1}.time;
  Dstime = 1.0;%domains{2,1}.time;
  StepT = 2.5;
  domains{1,1}.time=1.5;
  domains{2,1}.time=1.0;
  ndomains=2;
end

% Time Specs
ST = 0.1; % sampling time / time step
StepT = 0;
for i = 1:ndomains
    StepT=StepT+domains{i,1}.time;
end
LfowardT = StepT; % look forward time
StepCount = StepT/ST;
% desire terminal com velocity
desireV = hstep_leng/StepT/4;
desireV = 0;
Sstime = domains{1,1}.time;
Dstime = domains{2,1}.time;
N1 = floor(Sstime/StepT*StepCount);
N2 = floor(Dstime/StepT*StepCount);
StepSelection = blkdiag(ones(N1,1),ones(N2,1));
Step1 = StepSelection(:,1);
Step2 = StepSelection(:,2);
%%% LIPM
% Specs
T=ST;
N=LfowardT/ST;
g=9.8;
hcom=cog_h;
omega2= g/hcom;

% Discrete State Space
    % State Variable: COM COM_Vel, ZMP ( or COM_Acel)
%     syms omega2
%     syms T
%     syms g
%     syms hcom
%     A = [0 1 0 ;
%          omega2 0 -omega2 ;
%          0 0 0]
%     B = [0; 0; 1];
%     C = [1 0 -hcom/g];
%     sys = ss(A,B,C,[],T);

A = [1        T     0;
     omega2*T 1 -omega2*T;
     0        0     1;
     ];

B=[0;0;T];
C=[0 0 1 ];

AB = [A,B,zeros(3,N-1)];
CC=C*AB;
for i = 1:N-1
    ABtemp = [A*AB(end-2:end,1:3+i),B,zeros(3,N-1-i)];
    AB=[AB;ABtemp];
    CC=[CC;C*ABtemp];
end

Acom = AB(1:3:end,1:3);
Bcom = AB(1:3:end,4:end);
AcomV = AB(2:3:end,1:3);
BcomV = AB(2:3:end,4:end);
Azmp = CC(:,1:3);
Bzmp = CC(:,4:end);

%%% QP
w1=1; %c input
w2=0.1; % com 
w3=0.1; % comV

H =  w1*eye(N)+w2*Bcom'*Bcom+w3*BcomV'*BcomV ;
%%Insert Terminal COM position (hstep_leng) & Vel (desireV)
f = w2*Bcom'*(Acom*x0-hstep_leng)+w3*BcomV'*(AcomV*x0-desireV);

% ineqaulity constraint
Pzmp1 = 0.15;
Nzmp1 = 0.05;
Pzmp2 = 0.15;
Nzmp2 = 0.05;

Ax_ineq = [];
bx_ineq = [];

Azmp_Piq = Bzmp; bzmp_Piq = -Azmp*x0+Pzmp1*Step1+Pzmp2*Step2;
Azmp_Niq = -Bzmp; bzmp_Niq = Azmp*x0+Nzmp1*Step1+Nzmp2*Step2;

Ax_ineq = [Azmp_Piq;Azmp_Niq];
bx_ineq = [bzmp_Piq;bzmp_Niq];

Ax_eq = [];
bx_eq = [];

temp = Acom*x0;
Acom_eq = Bcom(end,:); bcom_eq = hstep_leng-temp(end);
temp = AcomV*x0;
AcomV_eq = BcomV(end,:); bcomV_eq = -temp(end);

Ax_eq = [Acom_eq;AcomV_eq];
bx_eq = [bcom_eq;bcomV_eq];

options = optimset('Algorithm','interior-point-convex','Display','final');

u0 = -H\f;

save('test.mat','H','f')

[u1,fval,exitflag] = quadprog(H,f,Ax_ineq,bx_ineq,Ax_eq,bx_eq,[],[],u0,options);
exitflag
%  Ucalc = -H^-1*f;
Z = Azmp*x0+Bzmp*u1;
C = Acom*x0+Bcom*u1;
CV = AcomV*x0+BcomV*u1;

%% Curve fitting
xdata = linspace(0,1,StepCount)';
ydata = C;
PolyOrder=7;
p = mmpolyfit(xdata,ydata,PolyOrder,'Point',[xdata(1) ydata(1);xdata(end) ydata(end)],'Slope',[xdata(1:2) zeros(2,1);xdata(end-1:end) zeros(2,1)]);
ydataLearned = polyval(p,xdata);
h=figure()
% set(h,'position',[100,100,1024*1.5,768*1.5]);
hold on
% plot(Zmprefx(1:N),'r','Linewidth',3,'LineStyle','-.')
plot(Z,'Linewidth',2)
plot(C,'g','Linewidth',2)
plot(ydataLearned,'m','Linewidth',2)
plot(CV,'kx','Markersize',3)
hold off
legend('ZMP from QP','COG from QP', 'Polyfit-7 order','COG Velocity')
%% Save as a function

syms x t real

for i=1:PolyOrder+1
    if i==1
    x = p(PolyOrder+2-i);
    tPoly = t/StepT;
    else
    x = x+p(PolyOrder+2-i)*tPoly;
    tPoly = tPoly*t/StepT;
    end
    
end

x = simplify(x);
dxdt = simplify(diff(x,t));
ddxdt = simplify(diff(dxdt,t));

matlabFunction(x, 'file', 'ZMPTrajGenX_v2');
matlabFunction(dxdt, 'file', 'ZMPTrajGenXder_v2');
matlabFunction(ddxdt, 'file', 'ZMPTrajGenXder2_v2');

end

