%> @brief get_spatial_model Get an extended model to use with spatial_v2 in
%> addition to customized spatial_v2 functions. Also loads dependencies.
%> @note To get additional motor inertias, be sure to use HandC_custom()
%> @author Eric Cousineau
function [model] = get_spatial_model(use_extra)

if nargin < 1 || isempty(use_extra)
    use_extra = true;
end

addpath('../matlab_utilities/matlab');
matlab_utilities_depends('sim', 'general', 'spatial_v2', 'yaml');

% Build sva model
body_struct = cell_to_matrix_scan(yaml_read_file('body_struct.yaml'));
model = planar_body_struct_to_spatial_model(body_struct);

% actuator types
% 0 - no actuator
% 1 - small actuator
% 2 - large actuator
motor_types = [0,0,0,2,2,2,2,2,2];


% Add boom and motor submodel
if use_extra
    % Add boom and motors
    % Boom - this is the effective linear inertia (mass) that stores energy
    % at the stance hip.
    %     m_boom = body_struct.boom.I(3, 3) / body_struct.boom.L^2;
    %     Isp_boom = mcI(m_boom, [0; 0], 0);
    %
    %     % Create minimal structure with boom
    %     % No potential energy, just inertia at the hip
    %     extra = struct();
    %     extra.NB = model.NB;
    %     extra.parent = model.parent;
    %     extra.jtype = model.jtype;
    %     extra.Xtree = model.Xtree;
    %     extra.I = model.I;
    %     extra.gravity = zeros(2, 1);
    %     for k = 1:extra.NB
    %         extra.I{k} = zeros(3);
    %     end
    %     % Add boom inertia to hip (dof: 3)
    %     dof_hip = 3;
    %     extra.I{dof_hip} = Isp_boom;
    
    % Add motor inertias - constant addition to inertia matrix
    %> @todo Handle floating-base case
    D_motor = zeros(model.NB);
    for i=1:model.NB
        if motor_types(i) == 1
            D_motor(i,i) = body_struct.motors.smallMotors.I;
        elseif motor_types(i) == 2
            D_motor(i,i) = body_struct.motors.largeMotors.I;
        end
    end
    % Using H in terms of spatial_v2 - inertia matrix
    %     extra.H_const = D_motor;
    model.H_const = D_motor;
end

end
