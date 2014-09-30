% configure domain type by its number
% domain #1: Both feet are flat on the ground
% domain #2: non-stance heel leaves the ground
% domain #3: non-stance toe leaves the ground
% domain #4: non-stance heel hits the ground
function [domain] = domainConfig_case3(type,a_opt,p_opt,vhip)


domain = struct();
domain.type = type; 

domain.a  = a_opt; % extend_a_mat(a_opt);
domain.p  = p_opt;
domain.vhip = vhip;

ep = 50;


controller = struct();
controller.useMuDirect = false;
controller.useUTransposeU = false;
controller.useUPrev = false;
controller.useSVD = false;
controller.uMax = 50000;
controller.extForceHACK = false;
controller.MaxConvergence = false;
controller.uRateLimit = nan;
controller.posDefSafety =1e-9;
controller.ZMP = false;
controller.frictionCoeff= [];

controller.opts = optimset('Display','off','Algorithm', 'interior-point-convex',...
    'TolX',1e-6,'TolFun',1e-6);


controller.useSprings = true;
controller.ep = ep;

controller.extForceMat = 1;



switch type
    case 1
        domain.kp_energy = 100;
        domain.controllerType = 'QP';
        domain.SwapLegs  = false;
        if controller.useSprings
            nClf = 5;
        else
            nClf = 5;
        end
        penalty = 5*ones(1,nClf);
        relaxation = -1*ones(1,nClf);        
        relDegree = 2*ones(1,nClf);
        controller.ExtForceDirections = [0,1,0,0,1,0];
    case 2 
        domain.kp_energy = 100;
        domain.controllerType = 'QP';
        domain.SwapLegs  = false;
        if controller.useSprings
            nClf = 5;
        else
            nClf = 6;
        end
        penalty = 5*ones(1,nClf);
        relaxation = -1*ones(1,nClf);        
        relDegree = 2*ones(1,nClf);
        controller.ExtForceDirections = [0,1,0,0,1];
    case 3
        domain.kp_energy = 100;
        domain.controllerType = 'QP';
        domain.SwapLegs  = false;
        if controller.useSprings
            nClf = 7;
        else
            nClf = 8;
        end
        penalty = 5*ones(1,nClf);
        relaxation = -1*ones(1,nClf);        
        relDegree = 2*ones(1,nClf);
        controller.ExtForceDirections = [0,1,0];
    case 4
        domain.kp_energy = 100;
        domain.controllerType = 'QP';
        domain.SwapLegs  = true;
        if controller.useSprings
            nClf = 5;
        else
            nClf = 6;
        end
        penalty = 5*ones(1,nClf);
        relaxation = -1*ones(1,nClf);        
        relDegree = 2*ones(1,nClf);
        controller.ExtForceDirections = [0,1,0,0,1];
    otherwise
        error('Unkown domain type.');
end

clfs = clf_construct(nClf,penalty,relaxation,relDegree,ep);

controller.clfs = clfs;

domain.qp = controller;
end


function [clfs] = clf_construct(nClf,penalty,relaxation,relDegree,ep)


clfs = struct([]);  
    
for i = 1:nClf
    nClfOutputs = 1;
    outputIndices = i;
    
    yIndices = i;
    dyIndices = i + nClf;
    clfs(i).penalty = penalty(i);
    clfs(i).relaxation = relaxation(i);
    
    
    % See if it's D1 (only one D1 output)
    etaIndices = [yIndices, dyIndices];
    
    
    
    qcare = care_gen(0, nClfOutputs, ep);
    
    clfs(i).care = qcare;
    clfs(i).nOutputs = nClfOutputs;
    clfs(i).etaIndices = etaIndices;
    clfs(i).outputIndices = outputIndices;
    clfs(i).relDegree = relDegree(i);
    
    
    
    
end

end

function [qcare] = care_gen(nD1, nD2, ep)

qcare = struct();
qcare.nD1 = nD1;
qcare.nD2 = nD2;
qcare.G = [
    eye(nD1)         zeros(nD1, nD2);
    zeros(nD2, nD1) zeros(nD2, nD2);
    zeros(nD2, nD1) eye(nD2)
    ];
qcare.F = [
    zeros(nD1, nD1)  zeros(nD1, 2*nD2);
    zeros(nD2, nD1) zeros(nD2, nD2)   eye(nD2);
    zeros(nD2, nD1) zeros(nD2, 2*nD2)
    ];
qcare.P = care(qcare.F, qcare.G, eye(nD1 + nD2 * 2));
qcare.C3 = min(eig(eye(nD1 + nD2*2))) /max(eig(qcare.P));

if nargin >= 3
    e = 1 / ep;
    qcare.emat = blkdiag(eye(nD1), 1 / e .* eye(nD2), eye(nD2));
    qcare.Pe = qcare.emat * qcare.P * qcare.emat;
end

end