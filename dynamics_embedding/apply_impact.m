function [dq_post] = apply_impact(x_pre,domain,ref)

nb = 2;
ndof = 7;
% 

qe_pre    = x_pre(1:nb+ndof);
dqe_pre   = x_pre(nb+ndof+1:end);

switch domain.type
%     case 1
%         % there is no impact, identity map
%         no_impact = true;
%     case 2
%         
%         Jh_stoe = Jh_stoe_mat(x_pre);
%         Jh_nsheel = Jh_nsheel_mat(x_pre);
%         %         Je = [Jh_stoe(1:2,:);Jh_nsheel(1:2,:)];
%         Je = Jh_nsheel(1:2,:);
%         no_impact = false;
    case 1
        % there is no impact, identity map
        Jh_stoe = Jh_stoe_mat(x_pre);
        Jh_nstoe = Jh_nstoe_mat(x_pre);
        %         Je = [Jh_stoe(1:2,:);Jh_nsheel(1:2,:)];
        Je = Jh_nstoe(1:2,:);
        no_impact = false;           

    case 2
          no_impact = true;      
     
%     case 3
%         
%         Jh_stoe = Jh_stoe_mat(x_pre);
%         Jh_nsheel = Jh_nsheel_mat(x_pre);
%         %         Je = [Jh_stoe(1:2,:);Jh_nsheel(1:2,:)];
%         Je = Jh_nsheel(1:2,:);
%         no_impact = false;     
%     case 4
%           no_impact = true;      
             
        
    otherwise
        error('Unknown domain Type.');
end

if no_impact
    dq_post = dqe_pre(nb+1:end);
else
    if ref.h.useSpatial
        model = ref.h.model; % Will be used later on
        De = HandC_extra(model, qe_pre, dqe_pre, ref.h.use_extra);
    else
        De    = De_mat(x_pre);
    end
    
    % Compute Dynamically Consistent Contact Null-Space from Lagrange
    % Multiplier Formulation
    %     DbarInv = Je * (De \ Je');
    %     I = eye(nb + ndof);
    %     Jbar = De \ (Je' / DbarInv);
    %     Nbar = I - Jbar * Je;
    %
    %     % Apply null-space projection
    %     dqe_post = Nbar * dqe_pre;
    nConstraints = size(Je,1);
    A = [De -Je'; Je zeros(nConstraints)];
    b = [De*dqe_pre; zeros(nConstraints,1)];
    y = A\b;
    VelFoot = y(1:2);
    ImpF = y(nb+ndof+1:end);
    dqe_post = y(1:nb+ndof);
    dq_post = dqe_post(nb+1:end);
end

end
