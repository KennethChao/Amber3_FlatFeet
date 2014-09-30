% configure domain type {'single','double'}
function [domain] = domainConfig_case2(type,a_opt,p_opt,vhip)



domain = struct();
domain.type = type;

% load('aFitSpline');

domain.ep = 10;
domain.a  = a_opt; % extend_a_mat(a_opt);
domain.p  = p_opt;
domain.vhip = vhip;

domain.umax = 30;
domain.paramTimeType = 'velocity';


qp = qp_construct(type);
domain.qp = qp;

if strcmp(type,'single')
    domain.kp_energy = 0;
    % construct single support domain    
    domain.controller = 'QP';
    domain.constraintType = 'SingleSupport';
    domain.impactType = 'discrete';
    domain.SwapLegs = true;
    domain.GroundHeight = 0;
    domain.HipPosition = NaN;
    domain.eventType = 'footHeight';
elseif strcmp(type, 'double')
    domain.kp_energy = 50;
    % construct double support domain
    domain.controller = 'QP';
    domain.constraintType = 'DoubleSupport';
    domain.impactType = 'continuous';
    domain.SwapLegs = false;
    domain.GroundHeight = NaN;
    domain.HipPosition = -0.0732;
    domain.LegLength = 0.85;
    domain.eventType = 'reactForce';
else
    error('Unknown domain type!!!');
end

end


function [qp] = qp_construct(type)

controller = struct();
controller.useMuDirect = false;
controller.useUTransposeU = false;
controller.useUPrev = false;
controller.useSVD = false;
controller.uMax = 3000;
controller.extForceHACK = false;
controller.useSprings = true;
controller.uRateLimit = nan;

controller.posDefSafety =1e-9;


controller.ZMP = false;
controller.frictionCoeff= [];

controller.MaxConvergence = false;
controller.extForceMat = 1;
controller.SSextForceDirections = [0,0,0,1];
controller.SSextForceMat = diag(controller.SSextForceDirections);
controller.DSextForceDirections = [0,0,0,1,0,1];
controller.DSextForceMat = diag(controller.DSextForceDirections);
    
if strcmp(type,'single')
    nClf = 5;
    clfs = struct([]);  
    penalty = [50 100 50 50 100];
    relaxation = [-1 -1 -1 -1 -1];
    controller.ep = 50;
    % for clf = clfs(:)'
    for i = 1:nClf
        nClfOutputs = 1;
        outputIndices = i;
        
        yIndices = i;
        dyIndices = i + nClf;
        clfs(i).penalty = penalty(i);
        clfs(i).relaxation = relaxation(i);
        
        
        % See if it's D1 (only one D1 output)
        etaIndices = [yIndices, dyIndices];
        
        
        relDegree = 2;
        qcare = care_gen(0, nClfOutputs, controller.ep);
        
        clfs(i).care = qcare;
        clfs(i).nOutputs = nClfOutputs;
        clfs(i).etaIndices = etaIndices;
        clfs(i).outputIndices = outputIndices;
        clfs(i).relDegree = relDegree;
        
        
        
        
    end
else
    nClf = 1;
    clfs = struct([]);
    controller.ep = 50;
    penalty = [100 100 50];
    relaxation = [-1 -1 -1];
    for i = 1:nClf
        
        nClfOutputs = 1;
        outputIndices = i;
        
        yIndices = i;
        dyIndices = i + nClf;
        clfs(i).penalty = penalty(i);
        clfs(i).relaxation = relaxation(i);
        % See if it's D1 (only one D1 output)
        etaIndices = [yIndices, dyIndices];
        
        
        relDegree = 2;
        qcare = care_gen(0, nClfOutputs, controller.ep);
        
        clfs(i).care = qcare;
        clfs(i).nOutputs = nClfOutputs;
        clfs(i).etaIndices = etaIndices;
        clfs(i).outputIndices = outputIndices;
        clfs(i).relDegree = relDegree;
        
    end
end

controller.clfs = clfs;

qp = controller;
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