function [dx] = calcEOM_Case8_1(t,x,domain,legId,ref,y0,f_firststep,f_laststep,domains,ndomains)

% global kspring kdamping

%% dynamics
ndof = 7;
nb = 2;
qe  = x(1:nb+ndof);
dqe = x(nb+ndof+1:end);
qb = qe(1:nb);
dqb = dqe(1:nb);
q  = qe(nb+1:end);
dq = dqe(nb+1:end);
% x = [q;dq];

% Lagrange Model
if ref.h.useSpatial
    model = ref.h.model; % Will be used later on
    [He, Ce] = HandC_extra(model, qe, dqe, ref.h.use_extra);
else
    He    = De_mat(x);
    Ce_o  = Ce_mat(x);
    Ge    = Ge_vec(x);
    Ce  = Ce_o*dqe + Ge;
end

% Constraints
[Je,Jedot,Be] = calcConstraintJacobian(domain,x);
% spring_joint_indices = [4 11]; % not used in case 8
% 
n_act = size(Be,2);% not used in case 8
n_ext = size(Je,1);% not used in case 8

% control parameters
ep = domain.qp.ep;

% Outputs
[ya, Jya, dJya, yd, dydt, ddydt, ...
    output_indices, ~, ~,~] = calcOutputsCase8(domain,x,t,f_firststep,f_laststep,domains);
% t
% not used in case 8
    % ya1 = Jdelta_phip_sca(x)*dqe;
    % yd1 = domain.vhip;
    % Jya1 = [zeros(1, nb+ndof),Jdelta_phip_sca(x)];
    % Jyd1 = zeros(1, 2*(nb+ndof));
    % y1 = ya1 - yd1;
    % Dy1 = Jya1 - Jyd1;
% not used in case 8

y2 = ya - (yd+y0);

Dy2 = Jya;% - Jyd;
DLfy2= dJya;% - dJyd;

% SLIP Dynamics
% [slip_dynamics,slip_force,pcom,pcom_dot,Energy] = calcSlipDynamics(domain,x);% not used in case 8

%% Controller
switch domain.controllerType
    case 'IO'
        Ie = eye(nb+ndof);
        XiInv = Je * (He \ Je');
%                 vfc = [dqe;-He\Ce];
%                 gfc = [zeros(size(Be));He\Be];  
        vfc = [
            dqe;
            He \ ((Je' * (XiInv \ (Je / He)) - Ie) * Ce - Je' * (XiInv \ Jedot * dqe))];
%         
        gfc = [
            zeros(size(Be));
            He \ (Ie - Je'* (XiInv \ (Je / He))) * Be];
        % Lie Derivatives
        Lfy2 = Dy2(output_indices,:)*vfc;
        % Second Lie Derivatives
        LfLfy2 = DLfy2(output_indices,:)*vfc;
        LgLfy2 = DLfy2(output_indices,:)*gfc;
        
        
        % feedback controller
        A = LgLfy2;
        u = -A \ (LfLfy2-ddydt(output_indices) + 2*ep*(Lfy2-dydt(output_indices)) + ep^2*y2(output_indices));
%         %         u = zeros(4,1);
%         u

        dx = vfc + gfc * u;
        Fe = -XiInv \ (Jedot * dqe + Je / He * (Be * u - Ce));

    case 'QP-MinNorm'

        Ie = eye(nb+ndof);
        XiInv = Je * (He \ Je');

        vfc = [
            dqe;
            He \ ((Je' * (XiInv \ (Je / He)) - Ie) * Ce - Je' * (XiInv \ Jedot * dqe))];
        gfc = [
            zeros(size(Be));
            He \ (Ie - Je'* (XiInv \ (Je / He))) * Be];        

        % first Lie Derivatives
        %         Lfy1 = Dy1*vf;
        Lfy2 = Dy2(output_indices,:)*vfc;
        %         % Second Lie Derivatives
        LfLfy2 = DLfy2(output_indices,:)*vfc;

        LgLfy2 = DLfy2(output_indices,:)*gfc;
    %% CLF with Min-Norm ctrl
         y2(1)=y2(1)*100; 
%         y2(2)=y2(2)*100;        
        etac = [ 
        y2(output_indices); 
        Lfy2(output_indices)-dydt(output_indices)];
        
        outputNumb=size(output_indices,2);
        F=[zeros(outputNumb,outputNumb) eye(outputNumb);zeros(outputNumb,outputNumb*2)];    
        G=[zeros(outputNumb,outputNumb); eye(outputNumb)];
        Q=eye(outputNumb*2);
        P=care(F,G,Q);  
        Pe = blkdiag(eye(outputNumb)*ep,eye(outputNumb))*P*blkdiag(eye(outputNumb)*ep,eye(outputNumb));
        Pe(1,1)=1000;
        V=etac'*Pe*etac;
        LfV = etac'*(F'*Pe+Pe*F)*etac;
        LgV = 2*etac'*Pe*G;

        [~, eigvl]=eig(Pe);
        lambdaPemax = max(max(eigvl));
%         [eigvc eigvl]=eig(Q);
%         lambdaQemin = min(min(eigvl));        
        gamma =1/lambdaPemax ;
        psi0= LfV+gamma*ep*V;
        psi1=LgV';
        if(psi0>0)
        mu =( -1*psi0*psi1/(psi1'*psi1));
        else
        mu = zeros(outputNumb,1);
        end
        
        % feedback controller
        A = LgLfy2;
        u = A\ (-LfLfy2+ddydt(output_indices)+mu );
        
%         u = u_qp(1:n_act);
%         Fe = u_qp(1+n_act:n_act+n_ext);

        Fe = -XiInv \ (Jedot * dqe + Je / He * (Be * u - Ce));
%         u
        dx = vfc + gfc * u;
    case 'QP-CLF'
        Ie = eye(nb+ndof);
        XiInv = Je * (He \ Je');
        vfc = [
            dqe;
            He \ ((Je' * (XiInv \ (Je / He)) - Ie) * Ce - Je' * (XiInv \ Jedot * dqe))];
        gfc = [
            zeros(size(Be));
            He \ (Ie - Je'* (XiInv \ (Je / He))) * Be];
        % First Lie Derivatives
        %         Lfy1 = Dy1*vf;
        Lfy2 = Dy2(output_indices,:)*vfc;
        % Second Lie Derivatives
        LfLfy2 = DLfy2(output_indices,:)*vfc;

        LgLfy2 = DLfy2(output_indices,:)*gfc;
%% CLF
        etac = [ 
        y2(output_indices); 
        Lfy2-dydt(output_indices)];
%         y2(1)=y2(1)*100;
        outputNumb=size(output_indices,2);
        F=[zeros(outputNumb,outputNumb) eye(outputNumb);zeros(outputNumb,outputNumb*2)];    
        G=[zeros(outputNumb,outputNumb); eye(outputNumb)];

        Q=eye(outputNumb*2);
        P=care(F,G,Q);  
        
        Pe = blkdiag(eye(outputNumb)*ep,eye(outputNumb))*P*blkdiag(eye(outputNumb)*ep,eye(outputNumb));

        V=etac'*Pe*etac;
        LfV = etac'*(F'*Pe+Pe*F)*etac;
        LgV = 2*etac'*Pe*G;

        [~, eigvl]=eig(Pe);
        lambdaPemax = max(max(eigvl));
       
        gamma =1/lambdaPemax ;
        psi0= LfV+gamma*ep*V;
        psi1=LgV';
        
        % feedback controller
        A = LgLfy2;        

%% u-based        
        if(det(A'*A)<10^-12)       
           try
           H=A'*A+eye(size(A,2))*10^-20;
%            det(A'*A)
           catch
           ppp=1;
           end

        else
        H=A'*A;
        end 
        f=2*(LfLfy2-ddydt(output_indices))'*A;
        %         AeqLagrangec = Je * (He \ [Be, Je']);
        %         beqLagrangec = -1 * (Jedot * dqe - Je * (He \ Ce));

        Aiq_CLF = psi1'*A;
        biq_CLF = -1*(psi0+psi1'*(LfLfy2-ddydt(output_indices)));
        %         Aiq=[-1,Aiq];

        % Quadprog        
        %     opts=optimset('Display','off','Algorithm', 'interior-point-convex',...
        %     'TolX',1e-6,'TolFun',1e-6)      

%         [u,~,EXITFLAG] = quadprog(2*H,f',Aiq,biq,[],[],[],[],[],domain.qp.opts);
%         if (domain.type==1)
%         EXITFLAG;
%         else
%         EXITFLAG            
%         end    
%         mu = A*u+(LfLfy2-ddydt(output_indices));
%% u-based w/ relaxation

%% w/ torqe constraint & ZMP constaint

% ZMP inequality constraints
        pos = jpos_mat(qe);
        pos(2,:) = [];  
       
%         xzmp = 1.0;
%         nxzmp =  -0.03;%-0.014-0.2;
        uMax=30;
        ActuatorNumb = 6;    
        
        Hp=blkdiag(10 ,H);
        fp=[0,f];       
        
        Aiq_CLFr=[-1,Aiq_CLF]; %w/ relaxation & CLF constrains
        
% torque inequality constraints added        
        Aiq_torqueH =[zeros(ActuatorNumb,1) eye(ActuatorNumb) ]; % torque Aiq: Higher bound
        Aiq_torqueL =[zeros(ActuatorNumb,1) eye(ActuatorNumb).*(-1) ]; % torque Aiq: Lower bound
       
        biq_torque = uMax.*ones(ActuatorNumb,1); % negative bound & positive bound share the same biq_torque
        
        Aiq = [Aiq_CLFr;Aiq_torqueH;Aiq_torqueL ]  ; %Gathered Aiq, biq        
        biq = [biq_CLF;biq_torque;biq_torque];            

        LfLfZMPx =  (-XiInv \ (Jedot * dqe + Je / He * ( - Ce)));
        LfLgZMPx =  (-XiInv\(Je/He*Be ));

 if domain.type ==1
 % calc.ZMP =-Fe(3)/Fe(2);  
 % Fe = -XiInv \ (Jedot * dqe + Je / He * (Be * u - Ce));

    if f_firststep==1
%         xzmp = 0.12;
%         nxzmp =  -0.02;
        xzmp = -0.03;
        nxzmp =  -0.12;
    else
        xzmp = -0.03;
        nxzmp =  -0.12;  
%         xzmp = 1.2+pos(11);
%         nxzmp = -0.03+pos(11);%-0.014-0.2;        
    end
        
        Aiq_ZMPH = [zeros(1,1), [0 -xzmp -1]*LfLgZMPx];
        Aiq_ZMPL = [zeros(1,1), [0 nxzmp 1]*LfLgZMPx]; 
        biq_ZMPH = [0 -xzmp -1]*(-LfLfZMPx);           
        biq_ZMPL = [0 nxzmp 1]*(-LfLfZMPx);   
        
        Aiq = [Aiq;Aiq_ZMPH ;Aiq_ZMPL];
        biq = [biq;biq_ZMPH ;biq_ZMPL];
 else
 % calc.ZMP =( Fe(3)+Fe(6)-pos(2)*Fe(2)-pos(12)*Fe(5) )/(-Fe(2)-Fe(5));
 
%         LfLfZMPx =  (-XiInv \ (Jedot * dqe + Je / He * ( - Ce)));
%         LfLgZMPx =  (-XiInv\(Je/He*Be ));
        xzmp = -0.03+0.1;
        nxzmp =  -0.12;
        pos = jpos_mat(qe);
        pos(2,:) = []; 
        
        Aiq_ZMPH = [zeros(1,1), [0 (-pos(2)+xzmp) 1 0 (-pos(12)+xzmp) 1]*(-LfLgZMPx)];
        Aiq_ZMPL = [zeros(1,1), [0 (-pos(2)+nxzmp) 1 0 (-pos(12)+nxzmp) 1]*(LfLgZMPx)]; 
        biq_ZMPH = [0 (-pos(2)+xzmp) 1 0 (-pos(12)+xzmp) 1]*(LfLfZMPx);           
        biq_ZMPL = [0 (-pos(2)+nxzmp) 1 0 (-pos(12)+nxzmp) 1]*(-LfLfZMPx); 
        
        Aiq_GRF_R = [zeros(1,1), [0 0 0 0 -1 0 ]*(LfLgZMPx)];
        biq_GRF_R = [0 0 0 0 -1 0 ]*(-LfLfZMPx)+5;           
        Aiq_GRF_L = [zeros(1,1), [0 -1 0 0 0 0 ]*(LfLgZMPx)];
        biq_GRF_L = [0 -1 0 0 0 0 ]*(-LfLfZMPx)+5;           
        
        Aiq = [Aiq;Aiq_GRF_R ;Aiq_GRF_L;Aiq_ZMPH;Aiq_ZMPL];
        biq = [biq;biq_GRF_R ;biq_GRF_L;biq_ZMPH;biq_ZMPL];  
%         Aiq = [Aiq;Aiq_ZMPH;Aiq_ZMPL];
%         biq = [biq;biq_ZMPH;biq_ZMPL];          
%         Aiq = [Aiq ];
%         biq = [biq ];     
   
%         Aiq = [Aiq;[zeros(ActuatorNumb,1) eye(ActuatorNumb)*-1 ];[zeros(ActuatorNumb,1) eye(6) ];[zeros(1,1),  [0 pos(2)-xzmp -1 0 pos(12)-xzmp -1]*(-XiInv\(Je/He*Be ))];[zeros(1,1), [0 -pos(2)+xzmp 1 0 -pos(12)+xzmp 1]*(-XiInv\(Je/He*Be ))] ]  ;
%         biq = [biq;uMax*ones(ActuatorNumb,1);uMax*ones(ActuatorNumb,1); [0 pos(2)-xzmp -1 0 pos(12)-xzmp -1]*(XiInv \ (Jedot * dqe + Je / He * ( - Ce)));[0 -pos(2)+xzmp 1 0 -pos(12)+xzmp 1]*(XiInv \ (Jedot * dqe + Je / He * ( - Ce))) ]; 

 end    
     [   up,xxx,existflag] = quadprog(2*Hp,fp',Aiq,biq,[],[],[],[],[],domain.qp.opts);   
      fprintf('\ntime = %e, existflag = %u',t,existflag)  ;     %existflag
     try
          u=up(2:7);
     catch
     [   up,xxx,existflag] = quadprog(2*Hp,fp',Aiq,biq,[],[],[],[],[],domain.qp.opts);          
         u
     end

        mu = A*u+(LfLfy2-ddydt(output_indices));
%% mu-based
%         Aiq = psi1';
%         biq = -psi0;
%         H=eye(outputNumb);
%         f = [];
% % Quadprog        
% %     opts=optimset('Display','off','Algorithm', 'interior-point-convex',...
% %     'TolX',1e-6,'TolFun',1e-6)      
%         [mu,FVAL,EXITFLAG] = quadprog(2*H,f,Aiq,biq,[],[],[],[],[],domain.qp.opts);
%         u = A\ (-LfLfy2+ddydt(output_indices)+mu );
%         if (domain.type==1)
%         EXITFLAG;
%         else
%         EXITFLAG            
%         end    
%         
%% Update Fe and dx
        Fe = -XiInv \ (Jedot * dqe + Je / He * (Be * u - Ce));
        dx = vfc + gfc * u;
      
   
    otherwise
        error('Invalid control: %s', domain.controller );
end

if nargin > 4
    
    switch domain.type
        case 1
%             Fc = Fe(5);
%             Fr = Fe([1 2 3 4]);
%             Fv_ns = Fe(4);
            Fc = 0;
            Fr = [Fe([1 2]);0;0];
            Fv_ns = 0;            
        case 2
%             Fc = 0;
%             Fr = [Fe([1 2]);0;0];
%             Fv_ns = 0;
            Fc = 0;
            Fr = Fe([1 2 3 4]);
            Fv_ns = 0;           
        case 3
            Fc = 0;
            Fr = Fe([1 2 3 4]);
            Fv_ns = 0;
        otherwise
            error('Unknown Domain Type.')
    end
    
%     ur = u(2:7);
%     us = u([1 8]);
    
    if legId == 0
        Fr = [Fr(3:4);Fr(1:2)];
    end
    

    calc = struct();
    calc.t =  t;
    calc.x = x;
%         calc.qe = qe;
%         calc.dqe = dqe;
%         calc.ddqe = dx(1:ndof+nb);
%     calc.ya1 = ya1;
%     calc.yd1 = yd1;
    calc.y2 = y2;
    calc.ya = ya;
    calc.yd = yd+y0;
    calc.dydt = dydt;
    calc.ddydt = ddydt;
    
    calc.ComAx = DLfy2(1,:)*vfc+DLfy2(1,:)*gfc*u;
%     calc.V = V_sca(x);
    calc.T = 1/2 * dqe' * He * dqe;
    calc.ur = u;
    if(domain.type==1)
    calc.Fe = [Fe;0;0;0];
    else
    calc.Fe = Fe;       
    end
    calc.Fr = Fr;
    calc.Fc = Fc;
    calc.q_anim = qe;

%     calc.ComVx= y_ComVX(x);        
    pos = jpos_mat(qe);
    pos(2,:) = [];      
    
    if domain.type ==1
    calc.ZMP = -Fe(3)/Fe(2);
    else
       

    calc.ZMP =( Fe(3)+Fe(6)-pos(2)*Fe(2)-pos(12)*Fe(5) )/(-Fe(2)-Fe(5));

    end
    
%         calc.Reg = calc.ComVx/(9.81/y0(2))^0.5+ya(1);        
    
    
    if legId == 1
        
%     calc.us = us;
%         calc.qs = qe(spring_joint_indices);
        calc.qe = qe;
        calc.dqe = dqe;
        calc.ddqe = dx(1+ndof+nb:end);
    else
%         calc.us = -flipud(us);        
        calc.qe = [qb;apply_relabel(q)];
        calc.dqe = [dqb;apply_relabel(dq)];
        calc.ddqe = [dqb;apply_relabel(dx(1+ndof+2*nb:end))];
%         calc.qs = calc.qe(spring_joint_indices);
    end
    calc.Fs = sqrt(Fr(1)^2+Fr(2)^2);
    calc.Fns = sqrt(Fr(3)^2+Fr(4)^2);
%     calc.pcom = pcom';
    if strcmp(domain.controllerType ,'QP-CLF')
    calc.Vclf=V;
    calc.Vdot=LfV+LgV*mu;   
%     calc.Vdot2=etac'*(Pe)*[Lfy2-dydt(output_indices);LfLfy2-ddydt(output_indices)+LgLfy2*mu]+[Lfy2-dydt(output_indices);LfLfy2-ddydt(output_indices)+LgLfy2*mu]'*(Pe)*etac;   
    calc.geV=-gamma*ep*V;  
    end    
%     calc.pcomdot = pcom_dot;
%     calc.Energy = Energy;
%     calc.Fv_ns = Fv_ns;
%     calc.tau = tau;
    ref.h.calc = calc;
end


end


