function [dx] = calcEOM(t,x,domain,legId,ref)


global kspring kdamping


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
% spring_joint_indices = [4 11];

% n_act = size(Be,2);
% n_ext = size(Je,1);


% control parameters
ep = domain.qp.ep;

% Outputs
[ya, Jya, dJya, yd, Jyd, dJyd, ...
    output_indices, slip_indices, force_indices,tau,dydt,ddydt] = calcOutputsCase8(domain,x,t);

% ya1 = Jdelta_phip_sca(x)*dqe;
% yd1 = domain.vhip;
% Jya1 = [zeros(1, nb+ndof),Jdelta_phip_sca(x)];
% Jyd1 = zeros(1, 2*(nb+ndof));
% y1 = ya1 - yd1;
% Dy1 = Jya1 - Jyd1;

y2 = ya - yd;
y2(1) =y2(1) +0.1881;
y2(4) =y2(4) +0.2;
y2(5) =y2(5) -0.05;

Dy2 = Jya;% - Jyd;
DLfy2= dJya;% - dJyd;

% SLIP Dynamics
% [slip_dynamics,slip_force,pcom,pcom_dot,Energy] = calcSlipDynamics(domain,x);

%% controller
switch domain.controllerType
    case 'IO'
        Ie = eye(nb+ndof);
        XiInv = Je * (He \ Je');
        vfc = [
            dqe;
            He \ ((Je' * (XiInv \ (Je / He)) - Ie) * Ce - Je' * (XiInv \ Jedot * dqe))];
        
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
        u = -A \ (LfLfy2-ddydt + 2*ep*(Lfy2-dydt) + ep^2*y2(output_indices));
%         %         u = zeros(4,1);
%         u
        y2;
        dx = vfc + gfc * u;
        Fe = -XiInv \ (Jedot * dqe + Je / He * (Be * u - Ce));
            
        
        
        

    case 'QP'
        vf = [dqe;-He\Ce];
        %         gf = [zeros(size(Be));He\Be];     
        
        
        AeqLagrangec = Je * (He \ [Be, Je']);
        beqLagrangec = -1 * (Jedot * dqe - Je * (He \ Ce));
        
        
        
        Bbarc = [Be, Je'];
        gbar_fc = [zeros(size(Bbarc)); He\Bbarc];
        
        
        
        % first Lie Derivatives
        Lfy1 = Dy1*vf;
        Lfy2 = Dy2*vf;
        % Second Lie Derivatives
        LfLfy2 = DLfy2*vf;
        
        Lf_gammac = [Lfy1;
            LfLfy2(output_indices)];
        
        etac = [y1; 
            y2(output_indices); 
            Lfy2(output_indices)];
        
        Lgbar_Lf_y1 = Dy1 * gbar_fc;
        Lgbar_Lf_y2 = DLfy2(output_indices,:) * gbar_fc;
        Abarc =[Lgbar_Lf_y1; 
            Lgbar_Lf_y2];
        
        if ~isempty(slip_indices)
            dJy_slip = dJy_NonLinearOutputs(x);
            AeqSlipDynamics = dJy_slip(slip_indices,:)*gbar_fc;
            beqSlipDynamics = slip_dynamics(slip_indices) - dJy_slip(slip_indices,:)*vf;
        else
            AeqSlipDynamics = [];
            beqSlipDynamics = [];
        end
        
        if isempty(force_indices)
            AeqSpringForce = zeros(1,n_act + n_ext);
            AeqSpringForce(1,1) = 1; % stance spring force
            AeqSpringForce(2,8) = 1; % non-stance spring force
            
            if domain.type == 2 || domain.type==1
                beqSpringForce = -kspring.*qe(spring_joint_indices);
            else
                beqSpringForce = -kspring.*qe(spring_joint_indices);
            end
            [u_qp] = clfMultiLagrange(domain.qp, etac, Abarc,Lf_gammac,...
                AeqLagrangec,beqLagrangec,AeqSlipDynamics,beqSlipDynamics,...
                AeqSpringForce,beqSpringForce);
        else
            Id_f = eye(n_ext);
            force_num = length(force_indices);
            nClf = length(domain.qp.clfs);
            AiqForces = [zeros(force_num,nClf), zeros(force_num,n_act),...
                Id_f(force_indices,:),-eye(force_num);...
                zeros(force_num,nClf), zeros(force_num,n_act),...
                -Id_f(force_indices,:),-eye(force_num);];
            
            biqForces = [slip_force(force_indices);-slip_force(force_indices)];
            
            AeqSpringForce = zeros(2,n_act + n_ext);
            AeqSpringForce(1,1) = 1; % stance spring force
            AeqSpringForce(2,8) = 1; % non-stance spring force
            
            if domain.type == 2
                beqSpringForce = -kspring.*qe(spring_joint_indices);
            else
                beqSpringForce = -kspring.*qe(spring_joint_indices);
            end
            [u_qp] = clfMultiLagrangeForce(domain.qp, etac, Abarc,Lf_gammac,...
                AeqLagrangec,beqLagrangec,AeqSlipDynamics,beqSlipDynamics,AiqForces,biqForces,...
                AeqSpringForce, beqSpringForce);
        end
        
        
        
        u = u_qp(1:n_act);
        Fe = u_qp(1+n_act:n_act+n_ext);

        
        Ie = eye(nb+ndof);
        XiInv = Je * (He \ Je');
        vfc = [
            dqe;
            He \ ((Je' * (XiInv \ (Je / He)) - Ie) * Ce - Je' * (XiInv \ Jedot * dqe))];
        
        gfc = [
            zeros(size(Be));
            He \ (Ie - Je'* (XiInv \ (Je / He))) * Be];
%         u
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
            Fc = 0;
            Fr = [Fe([1 2]);0;0];
            Fv_ns = 0;
        case 3
            Fc = 0;
            Fr = Fe([1 2 3 4]);
            Fv_ns = 0;
        otherwise
            error('Unknown Domain Type.')
    end
    
    ur = u(2:5);
    us = u([1 6]);
    
    if legId == 0
        Fr = [Fr(3:4);Fr(1:2)];
    end
    

    calc = struct();
    calc.t =  t;
    calc.x = x;
        calc.qe = qe;
        calc.dqe = dqe;
        calc.ddqe = dx(1:ndof+nb);
%     calc.ya1 = ya1;
%     calc.yd1 = yd1;
    calc.y2 = y2;
    calc.ya = ya;
    calc.yd = yd;
    calc.dydt = dydt;
    calc.ddydt = ddydt;
    
    calc.V = V_sca(x);
    calc.T = 1/2 * dqe' * He * dqe;
    calc.ur = ur;
    %     calc.Fe = Fe;
    calc.Fr = Fr;
    calc.Fc = Fc;
    calc.q_anim = qe;
    if legId == 1
        
%         calc.us = us;
%         calc.qs = qe(spring_joint_indices);
        calc.qe = qe;
        calc.dqe = dqe;
        calc.ddqe = dx(1+ndof+nb:end);
    else
        calc.us = -flipud(us);        
        calc.qe = [qb;apply_relabel(q)];
        calc.dqe = [dqb;apply_relabel(dq)];
        calc.ddqe = [dqb;apply_relabel(dx(1+ndof+2*nb:end))];
        calc.qs = calc.qe(spring_joint_indices);
    end
    calc.Fs = sqrt(Fr(1)^2+Fr(2)^2);
    calc.Fns = sqrt(Fr(3)^2+Fr(4)^2);
%     calc.pcom = pcom';
%     calc.pcomdot = pcom_dot;
%     calc.Energy = Energy;
%     calc.Fv_ns = Fv_ns;
%     calc.tau = tau;


        
    ref.h.calc = calc;
end


end


