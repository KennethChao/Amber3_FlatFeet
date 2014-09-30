%> @brief HandC_custom Perform HandC as normal, but add in model.H_const to
%> inertia matrix (for motor inertias)
%> @brief use_extra Incorporate boom and motor inertia dynamics
function [H, C, Xbase] = HandC_extra(model, q, dq, use_extra, fext)

%> @todo This can be made much quicker for the boom dynamics by making a
%> custom modification:
%> Use kinematics computed in the normal routine, and simply and the
%> inertia and bias terms due to the boom dynamics as a separate inertial
%> element.
%> CHALLENGE: How to exclude gravity from the 'avp' variable?

if nargin < 4
    use_extra = true;
end
if nargin < 5
    fext = {};
end

[H, C, Xbase] = HandC_kinematic(model, q, dq, fext);

if use_extra
    % Add extra dynamics
    %     [H_extra, C_extra] = HandC(model.extra, q, dq); % fext?
    H = H +  model.H_const;
    %     C = C + C_extra;
end

end
