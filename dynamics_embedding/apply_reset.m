function [x_plus] = apply_reset(x_pre,domain,ref)

ndof = 7;
nb = 2;

q_post = x_pre(nb+1:nb+ndof);
dq_post = apply_impact(x_pre,domain,ref);
qb = zeros(nb,1);
dqb = zeros(nb,1);

if domain.SwapLegs
    q_plus  = apply_relabel(q_post);
    dq_plus = apply_relabel(dq_post);
    qe_post = [qb;q_plus];
    dqe_post = [dqb;dq_plus];
    x_plus = [
        qe_post;
        dqe_post
        ];
else
    qe_post = [qb;q_post];
    dqe_post = [dqb;dq_post];
    x_plus = [
        qe_post;
        dqe_post
        ];
end

% vcom_pre = Je_com_mat(x_pre)*x_pre(1+nb+ndof:end)
% 
% vcom_plus = Je_com_mat(x_plus)*x_plus(1+nb+ndof:end)
end
