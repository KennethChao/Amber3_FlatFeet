function [] = walking_slip(nsteps)

global use_one_ecwf use_cwf

if nargin < 1
    nsteps = 1;
end

% options
use_spatial  = true;
use_extra    = true;
use_one_ecwf = true;
use_cwf      = true;

addpath('../matlab_utilities/matlab');
matlab_utilities_depends('general', 'sim');

addpath('./expr/model');
addpath('./expr/control');
addpath('./expr/outputs');
addpath('./slip/');

% nsteps = 1;
global Energy0 kspring grav mTotal L0 kdamping
load('slip_orbit_ic.mat','sliprobot');

L0 = sliprobot.L0;
mTotal = sliprobot.m;
grav = 9.81;
kspring = sliprobot.k;
Energy0 = sliprobot.E;
kdamping = 2*sqrt(mTotal*kspring);

ndof = 9;
nb = 2;

domain_time = 10; % s


x_minus = convert_ic_com_case4(false);
% 
if use_one_ecwf
    %     load('slip/case4/A_fit_full_case4.mat','a_opt');
    %     a_fit = a_opt;
    if use_cwf
        load('slip/case4/A_fit_full_case4.mat','a_fit');
    else
        load('slip/case4/A_bezier_full_case4','a_fit');
    end
end
load('slip/case4/p_opt_full_case4');

ndomains = 3;
domains = cell(ndomains,1);
for i = 1:ndomains
    if ~use_one_ecwf
        load(['slip/case4/A_fit_domain',num2str(i),'_case4.mat'],'a_fit');
    end
    domains{i} = domainConfig_case4(i,a_fit,p_opt,vhip);
end


ref = Ref();
ref.h = struct();
ref.h.calcs = {};
ref.h.calc = {};

ref.h.useSpatial = use_spatial;
ref.h.use_extra  = use_extra;
model = get_spatial_model(ref.h.use_extra);
ref.h.model = model;


calcs_set = cell(nsteps*ndomains,1);

t_plus = 0;
% save('xb_minus','xb_minus')
x_plus = apply_reset(x_minus,domains{end},ref);
% xb_plus = convert_ic_plus();
x0 = x_plus;
legId = 1;
for k = 1:nsteps
    pcm_i = pe_com_vec(x0);
    t_i = t_plus;
    tic    
    for j = 1:ndomains
        
        odeopts = odeset('MaxStep', 1e-2, 'OutputFcn', @outputfcn, ...
            'Events', @(t, x, ref) eventfcn(t, x,domains{j},ref),'RelTol',1e-6,'AbsTol',1e-6);
        ref.h.calcs = {};
        
        
        sol = ode45(@(t, x, ref) calcEOM(t,x,domains{j},legId,ref), [0, domain_time], x0, ...
            odeopts, ref);
        
        % Store calculations
        calcs = horzcat_fields(cell2mat(ref.h.calcs));
        calcs.t = t_plus + calcs.t;
        calcs_set{(k-1)*ndomains+j} = calcs;
        
        % Reset for next integration
        t_plus = t_plus + sol.x(end);
        x_minus = sol.y(:, end);
        x_plus  = apply_reset(x_minus,domains{j},ref);
        
        %         if j == 1
        %             no = size(actuals, 1);
        %             a_opt = domains{j+1}.a;
        %             for i = 1:no
        %                 cur = actuals{i,1};
        %                 Jcur = actuals{i,2};
        %                 a0 = 0.1.*a_opt(i,[2 4 6]);
        %                 a_opt(i,:) = find_transition(a0,xb_plus,xb_minus,cur,Jcur);
        %             end
        %             domains{j+1}.a = a_opt;
        %         end
        x0 = x_plus;
    end
    
    toc
    
    
    legId = 1 - legId;
    %     pe_nsm = pe_nsm_vec(x_minus);
    pcm_f = pe_com_vec(x_minus);
    te = t_plus-t_i;
    speed = (pcm_f(1)-pcm_i(1))/te;
    %     q_ns = atan2(pcm_f(2)-pe_nsm(2),pcm_f(1)-pe_nsm(1));
    disp(['Step ',num2str(k),...
        ' completed at time t = ', num2str(te),...
        ' Speed = ', num2str(speed)]);%,...
        %         ' Strike length = ', num2str(pe_nsm(1))
        %         ' Touch-down angle = ', num2str((pi-q_ns)*180/pi)]);
    
    
end


calcs_full = [calcs_set{:}];
out = horzcat_fields_domains(calcs_full);
save('figs/simulation_data','out','calcs_set','Energy0','kspring','mTotal','L0');


% Animate
aniplot();


%% Plot Data
plotData();




end




