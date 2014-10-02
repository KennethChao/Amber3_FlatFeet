% function [] = walking_Case8(nsteps)
close all;

first_step = 1;
last_step = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Options
    opt = struct();
    nsteps =2;                  % Domains number (1SSP/1DSP: nsteps=2)
    opt.ctrl_method = 'IO'; % Controller mode 'IO' 'QP-MinNorm' 'QP-CLF'
    opt.dt1 = 2;                % 1st domain time (SSP)
    opt.dt2 = 1;                % 2nd domain time (DSP)
    step_h = 0.1;               % half step length
    hstep_l = 0.1;              % Step height
%     use_replan_estep = true;
% options
    use_spatial  = true  
    use_extra    = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Options    

addpath('../matlab_utilities/matlab');
matlab_utilities_depends('general', 'sim');

addpath('./expr/model');
addpath('./expr/control');
% addpath('./build_torso');


%%
ndof = 7;   %%%%%%%%%%%%%%%%%%%%%%% For Case8
nb = 2;

%%
x_minus = [zeros(ndof+nb,1);zeros(nb,1);ones(ndof,1)*0.0];
% Initial Configuration (Squad)
ang = 14/360*2*pi;
    x_minus(4)=-ang;
    x_minus(5)=2*ang;
    x_minus(6)= -ang;
    x_minus(7)=-ang;
    x_minus(8)=2*ang;
    x_minus(9)=-ang;
    
% Plot_Configuration
    Data = jpos_mat(x_minus(1:ndof+nb,1));

    plot(Data(1,:),Data(2,:),'-o')
%     hold on
% %     plot(y_ComX(x_minus(1:ndof+nb,1)),y_ComZ(x_minus(1:ndof+nb,1)),'r*')
%     hold off
    
%     y_ComX(x_minus(1:ndof+nb,1))
%     y_ComZ(x_minus(1:ndof+nb,1))

ndomains = 2; %%%%%%%%%%%%%%%%%%%%%%% For Case8: One step
domains = cell(ndomains,1);
for i = 1:ndomains
    domains{i} = domainConfig_Case8(i,opt);
end

FootTrajGen(hstep_l,step_h);
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

y0=zeros(6,1);
    y0(1)= y_ComX(x0);
    y0(2)= y_ComZ(x0);
    y0(3)= y_TorsoAngle(x0);
    y0(4)= y_NonStanceFootX(x0);
    y0(5)= y_NonStanceFootZ(x0);
    y0(6)= y_NonStanceFootAng(x0);
%     y0 = ya2_vec(x0);

    ComX=y0(1);% y_ComX(x0);
    ComVx=0; %y_ComVX(x0);        
    ComZ= y0(2);%y_ComZ(x0);       
    Zmp =y0(1);% y_ComX(x0);
%     ComAx = 9.81/ComZ*(ComX-Zmp);
    Xinitial = ComX;    
    X0 =[ComX-Xinitial;ComVx;Zmp-Xinitial];
Index =1;              
     %% Com terminal position/ velocity setting inside
%      COM_plan( ComZ, hstep_l,X0,domains{Index,1}.type, domains,ndomains)
legId = 1;

profile on
for k = 1:nsteps
%     pcm_i = pe_com_vec(x0);
    t_i = t_plus;
    tic    
%     COMTrajGen(0.15,0.1);
   
    odeopts = odeset('MaxStep', 1e-2, 'OutputFcn', @outputfcn, ...
         'RelTol',1e-6,'AbsTol',1e-6);
    ref.h.calcs = {};
    if(k>1)
        first_step=0;
    end
    if(k==nsteps&&k~=1)
        last_step=1;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Options           
%     if(k==2) % For debugging
%         Stime=0.03;
%     else
%         Stime =domains{Index,1}.time;
%     end        
  Stime =domains{Index,1}.time; % For normal domain time, please uncomment it

%   Replan for each step
%     if(mod(k,2)==1)
%         COM_plan( ComZ, hstep_l,X0,domains{Index,1}.type, domains,ndomains)
%     end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Options   
        sol = ode45(@(t, x, ref) calcEOM_Case8_1(t,x,domains{Index},legId,ref,y0,first_step,last_step,domains,ndomains), [0, Stime], x0, ...domains{1,1}.time
            odeopts, ref);
        display('nsteps = '  );
        display(k);
        display('ndomains = ' );
        display(Index) ;        
        % Store calculations
        calcs = horzcat_fields(cell2mat(ref.h.calcs));
        calcs.t = t_plus + calcs.t;
 
        
        % Reset for next integration
        t_plus = t_plus + sol.x(end);
        x_minus = sol.y(:, end);
        x_plus  = apply_reset(x_minus,domains{Index},ref);
        x0 = x_plus;
        Data1=jpos_mat(x_minus(1:ndof+nb,1));
        plot(Data1(1,:),Data1(2,:),'-o')
        hold on
        Data2=jpos_mat(x_plus(1:ndof+nb,1));
        plot(Data2(1,:)+0.1,Data2(2,:),'r-o')        
        hold off
        calcs_set{k} = calcs;
            temp(1)= y_ComX(x0);
            temp(2)= y_ComZ(x0);
            temp(3)= 0;%y_TorsoAngle(x0);
            temp(4)= y_NonStanceFootX(x0);
            temp(5)= y_NonStanceFootZ(x0);
            temp(6)= 0;%y_NonStanceFootAng(x0);
            
%             ComX= y_ComX(x0);
%             ComVx= y_ComVX(x0);        
%             ComZ= y_ComZ(x0);       
%             Zmp = calcs.ZMP(end);
%             ComAx = 9.81/ComZ*(-(Zmp-Xinitial)+(ComX-Xinitial));
%             ComAx2 = calcs.ComAx(end);
%             X0 =[ComX-Xinitial;ComVx;Zmp-Xinitial-0.1];
%               temp = ya2_vec(x0); 
                 
      if Index==1      
          Index=2;
          for q=2:6
           y0(q) = temp(q);
          end 
      else
          Index=1; 
          for q=1:6
           y0(q) = temp(q);
          end           
      end
    toc
      
    legId = 1 - legId;
    %     pe_nsm = pe_nsm_vec(x_minus);
%     pcm_f = pe_com_vec(x_minus);
    te = t_plus-t_i;
%     speed = (pcm_f(1)-pcm_i(1))/te;
    %     q_ns = atan2(pcm_f(2)-pe_nsm(2),pcm_f(1)-pe_nsm(1));
    disp(['Step ',num2str(k),...
        ' completed at time t = ', num2str(te)] );%,...
        %         ' Speed = ', num2str(speed)]);%,...
        %         ' Strike length = ', num2str(pe_nsm(1))
        %         ' Touch-down angle = ', num2str((pi-q_ns)*180/pi)]);
        if(mod(k,5)==0) % For multiple steps checking
            pause();
        end
     
end

profile viewer
profile resume
profile off

calcs_full = [calcs_set{:}];
out = horzcat_fields_domains(calcs_full);
save('figs\simulation_data','out','calcs_set');

% Animate
aniplot();

%% Plot Data
plotData();

% end




