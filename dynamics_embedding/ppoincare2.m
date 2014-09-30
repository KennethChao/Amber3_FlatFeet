function [x_plus] = ppoincare2()

% load robot model

addpath('../matlab_utilities/matlab');
matlab_utilities_depends('general', 'sim');

addpath('./expr/model');
addpath('./expr/control');

% nsteps = 1;
Energy0 = 497.3859;

ndof = 9;
nb = 2;

% xb_minus = convert_ic();
xb_minus = convert_ic_com();
% load('A_fit_single','a_fit_single');
load('A_fit_cwf','a_fit');
% a_fit = a_fit_double;
load('p_opt','p_opt','vhip');
% a_fit_double = zeros(4,7);
% a_fit_double(6:7,:) = [a_fit_double(5,:);a_fit_double(5,:)];
% a_fit_double = a_fit_single;
% a_fit_double(2,:) = a_fit_single(1,:);
ds_domain = domainConfig('double',a_fit,p_opt,vhip);
% 
% load('A_fit_single','a_fit_single');
% load('p_opt_single','p_opt','vhip');
ss_domain = domainConfig('single',a_fit,p_opt,vhip);
domains = {ds_domain,ss_domain};
q_minus = xb_minus(nb+1:nb+ndof);
dq_minus = xb_minus(2*nb+ndof+1:end);

x0 = [q_minus;dq_minus];
tol     = 1e-20;
opts    = optimset('display', 'iter', 'tolfun', tol);
ref = Ref();
ref.h = struct();
ref.h.calcs = {};
ref.h.calc = {};
[x,~,exitflag,~,DP] = fsolve(@(x) x-pacre2(x, domains,Energy0,ref),x0,opts);


abs(eig(eye(18)-DP))

% q_minus = [-0.05;x(1:model.dim_q-1)];
% dq_minus = x(model.dim_q:end);
q_minus = x(1:model.dim_q);
dq_minus = x(1+model.dim_q:end);
x_minus = [0;0;q_minus;0;0;dq_minus];
save('x_minus');


    
end