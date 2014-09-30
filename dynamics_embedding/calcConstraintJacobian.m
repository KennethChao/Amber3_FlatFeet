function [Je,Jedot,Be] = calcConstraintJacobian(domain,x)


nb = 2;
ndof = 7;

Id_b = eye(nb+ndof);
% Constraints


act_indices = 4:9;
Be = Id_b(:,act_indices);

switch domain.type
    case 1
%         Jh_stoe = Jh_stoe_mat(x);
%         Jh_nstoe = Jh_nstoe_mat(x);
%         Je = [Jh_stoe(1:2,:);Jh_nstoe(1:2,:);Jh_stoe(3,:)];
%         
%         dJh_stoe = dJh_stoe_mat(x);
%         dJh_nstoe = dJh_nstoe_mat(x);
%         Jedot = [dJh_stoe(1:2,:);dJh_nstoe(1:2,:);dJh_nstoe(3,:)];
        Jh_stoe = Jh_stoe_mat(x);
        Je = Jh_stoe;
        
        dJh_stoe = dJh_stoe_mat(x);
        Jedot = dJh_stoe;
    case 2
%         Jh_stoe = Jh_stoe_mat(x);
%         Je = Jh_stoe;
%         
%         dJh_stoe = dJh_stoe_mat(x);
%         Jedot = dJh_stoe;

         Jh_stoe = Jh_stoe_mat(x);
         Jh_nstoe = Jh_nstoe_mat(x);
         Je = [Jh_stoe;Jh_nstoe];
%         
         dJh_stoe = dJh_stoe_mat(x);
         dJh_nstoe = dJh_nstoe_mat(x);
         Jedot = [dJh_stoe;dJh_nstoe];
         
%     case 3
%         Jh_stoe = Jh_stoe_mat(x);
%         Jh_nsheel = Jh_nsheel_mat(x);
%         Je = [Jh_stoe(1:2,:);Jh_nsheel(1:2,:)];
%         
%         dJh_stoe = dJh_stoe_mat(x);
%         dJh_nsheel = dJh_nsheel_mat(x);
%         Jedot = [dJh_stoe(1:2,:);dJh_nsheel(1:2,:)];
    otherwise
end





end