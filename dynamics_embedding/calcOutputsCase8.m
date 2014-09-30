function [ya, Jya, dJya, yd, dydt,ddydt, ...
          output_indices, slip_indices, force_indices,tau] = calcOutputsCase8(domain,x,t,f_firststep,f_laststep,domains)
 
% global use_cwf not used in case 8
nb = 2;
ndof = 7;
qe = x(1:nb+ndof);
dqe = x(nb+ndof+1:end);
% a  = domain.a;
% p  = domain.p;
% vhip = domain.vhip;

%% outputs
% Actual outputs
switch domain.type
    case 1
        if domain.qp.useSprings
            slip_indices = 2;
            output_indices = [1 2 3 4 5 6];%[1 2 3 4];
            force_indices = [];
        else
            slip_indices = [];
            output_indices = [1 2 3 4 5 6];%[1 2 3 4 7];
            force_indices = [];
        end
    case 2
        if domain.qp.useSprings
            slip_indices = 2;
            output_indices = [1 2 3 4 6];
            force_indices = [];
        else
            slip_indices = [];
            output_indices = [1 2 3];
            force_indices = [];
        end
    case 3
        if domain.qp.useSprings
            slip_indices = 2;
            output_indices = [1 2 3 4];
            force_indices = [];
        else
            slip_indices = [];
            output_indices = [1 2 3 4 5 7];
            force_indices = [];
        end
    otherwise
        error('Unknown domain Type.');
end
    

% actuals = {
%         @y_StanceKnee, @Jy_StanceKnee;
%         @y_NonStanceKnee, @Jy_NonStanceKnee;     
%         @y_TorsoAngle, @Jy_TorsoAngle;        
%         @y_NonStanceAnkle, @Jy_NonStanceAnkle;
%         @y_NonStanceSpring, @Jy_NonStanceSpring;
%         @y_DeltaNsSlope, @Jy_DeltaNsSlope;
%         @y_StanceSpring, @Jy_StanceSpring;
%         };

% outputs = {ComX, ComZ, TorsoAngle, NonStanceFootX, NonStanceFootZ,   NonStanceFootAng}
actuals = {
        @y_ComX, @Jy_ComX;
        @y_ComZ, @Jy_ComZ;     
        @y_TorsoAngle, @Jy_TorsoAngle;        
        @y_NonStanceFootX, @Jy_NonStanceFootX;
        @y_NonStanceFootZ, @Jy_NonStanceFootZ;
        @y_NonStanceFootAng, @Jy_NonStanceFootAng;
        };
% actuals = {
%         @y_ComX, @Jy_ComX, @dJy_ComX;
%         @y_ComZ, @Jy_ComZ, @dJy_ComZ;
%         @y_TorsoAngle, @Jy_TorsoAngle, @dJy_TorsoAngle;      
%         @y_NonStanceFootX, @Jy_NonStanceFootX, @dJy_NonStanceFootX;
%         @y_NonStanceFootZ, @Jy_NonStanceFootZ,@dJy_NonStanceFootZ;
%         @y_NonStanceFootAng, @Jy_NonStanceFootAng, @dJy_NonStanceFootAng;
%         };    
% 
 no = size(actuals, 1);
ya = zeros(no, 1);
Jya = zeros(no, 2*(nb+ndof));
dJya = zeros(no, 2*(nb+ndof));
% % 
% for i = 1:no
%     % Compute actuals
%          cur = actuals{i, 1};
%     Jcur = actuals{i, 2};
%     %     dJcur = actuals{i, 3};
%      dJcur = actuals{i, 3};
%      ya(i) = cur(x);
%     Jya(i, :) = Jcur(x);
% %     ya(i) = Jya(i,1:nb+ndof)*qe;    
% %        dJya(i,1+nb+ndof:end) = Jya(i,1:nb+ndof);
%        dJya(i,:) =dJcur(x);
% end

%  ya = yout(x);
% % sum(ya2-ya)
%  Jya = Jyout(x);
% sum(sum(Jya2-Jya))
%  dJya = dJy_out(x);

 ya = ya2_vec(x);
 Jya = Dya2_mat(x);
 dJya = DLfya2_mat(x);


% sum(sum(dJya2-dJya))
% ya2 = y_out(x);
% sum(ya2-ya)
% Purely state-based desired outputs
yd = zeros(no, 1);
dydt = zeros(no, 1);
ddydt= zeros(no, 1);
Jyd = zeros(no, 2*(nb+ndof));
dJyd = zeros(no, 2*(nb+ndof));
% phip = pe_com_vec(x);
% Jphip = Je_com_mat(x);
% dJphip = dJe_com_mat(x);

% phip = delta_phip_sca(x);
% Jphip = Jdelta_phip_sca(x);
% dJphip = zeros(1,nb+ndof);
% dtdp = 1/vhip;
 tau   = 0;%(phip(1) - p(1))/vhip;


%% COMZ: 0.9072 Total Mass:79.7761

% if use_cwf
%     yd = cwf_time(a,tau);
%     dydt = dcwf_time(a,tau);
%     ddydt = d2cwf_time(a,tau);
% else
%     yd = bezier(a,tau);
%     dydt = dbezier(a,tau);
%     ddydt = d2bezier(a,tau);
% end
% global first_step last_step
 FirstStepDSP=1;  
if f_firststep==1
   StrideGain=1;
%    FirstStepDSP=0;
%    first_step=false;
else
   StrideGain=2; 
%    FirstStepDSP=1;   
end
if f_laststep==8
   MovGain=0;
   StrideGain=1;   
%    last_step=false;
else
   MovGain=1;    
end        
       
switch domain.type
    case 1

yd(1) = ZMPTrajGenX_v2(t);%com_onestep_time(t);
yd(2) = 0;
yd(3) = 0;
yd(4) = FootTrajGenX(t/domain.time)*StrideGain;
yd(5) = FootTrajGenZ(t/domain.time);
yd(6) = 0;

dydt(1) = ZMPTrajGenXder_v2(t);%com_onestep_time_der(t);
dydt(2) = 0;
dydt(3) = 0;
dydt(4) = FootTrajGenXder(t/domain.time)/domain.time*StrideGain;
dydt(5) = FootTrajGenZder(t/domain.time)/domain.time;
dydt(6) = 0;


ddydt(1) = ZMPTrajGenXder2_v2(t);%com_onestep_time_der2(t);
ddydt(2) = 0;
ddydt(3) = 0;
ddydt(4) = FootTrajGenXder2(t/domain.time)/domain.time^2*StrideGain;
ddydt(5) = FootTrajGenZder2(t/domain.time)/domain.time^2;
ddydt(6) = 0;
    case 2
 
        
yd(1) = ZMPTrajGenX_v2(t+domains{1,1}.time)*MovGain*FirstStepDSP;
yd(2) = 0;
yd(3) = 0;
yd(4) = 0;
yd(5) = 0;
yd(6) = 0;

dydt(1) = ZMPTrajGenXder_v2(t+domains{1,1}.time)*MovGain*FirstStepDSP;
dydt(2) = 0;
dydt(3) = 0;
dydt(4) = 0;
dydt(5) = 0;
dydt(6) = 0;


ddydt(1) = ZMPTrajGenXder2_v2(t+domains{1,1}.time)*MovGain*FirstStepDSP;
ddydt(2) = 0;
ddydt(3) = 0;
ddydt(4) = 0;
ddydt(5) = 0;
ddydt(6) = 0;        

    otherwise
        error('Unknown domain Type.');
end

% for i = 1:no
% %     Jyd(i, :) = [dydt(i) * Jphip(1,:)*dtdp,zeros(1,nb+ndof)];
% %     dJyd(i, :) = [ddydt(i) * ((Jphip(1,:)*dtdp)*dqe) * (Jphip(1,:)*dtdp) + ...
% %         dydt(i) * dJphip(1,:),...
% %         dydt(i) * (Jphip(1,:)*dtdp)];
% % %     Jyd(i, :) = [zeros(1,nb+ndof),zeros(1,nb+ndof)];
% % %     dJyd(i, :) = [zeros(1,nb+ndof),zeros(1,nb+ndof)];    
%     
% end

end