function [ya, Jya, dJya, yd, Jyd, dJyd, ...
    output_indices, slip_indices, force_indices,tau] = calcOutputs(domain,x)
 
global use_cwf
nb = 2;
ndof = 7;
qe = x(1:nb+ndof);
dqe = x(nb+ndof+1:end);
a  = domain.a;
p  = domain.p;
vhip = domain.vhip;

%% outputs
% Actual outputs
switch domain.type
    case 1
        if domain.qp.useSprings
            slip_indices = 2;
            output_indices = [1 2 3 4 5 6 ];%[1 2 3 4];
            force_indices = [];
        else
            slip_indices = [];
            output_indices = [1 2 3 4 5 6 ];%[1 2 3 4 7];
            force_indices = [];
        end
    case 2
        if domain.qp.useSprings
            slip_indices = 2;
            output_indices = [1 2 3 4 6];
            force_indices = [];
        else
            slip_indices = [];
            output_indices = [1 2 3 4 5 6 7];
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


no = size(actuals, 1);
ya = zeros(no, 1);
Jya = zeros(no, 2*(nb+ndof));
dJya = zeros(no, 2*(nb+ndof));

for i = 1:no
    % Compute actuals
    %     cur = actuals{i, 1};
    Jcur = actuals{i, 2};
    %     dJcur = actuals{i, 3};
    Jya(i, :) = Jcur(x);
    ya(i) = Jya(i,1:nb+ndof)*qe;    
    dJya(i,1+nb+ndof:end) = Jya(i,1:nb+ndof);
end


% Purely state-based desired outputs
yd = zeros(no, 1);
Jyd = zeros(no, 2*(nb+ndof));
dJyd = zeros(no, 2*(nb+ndof));
% phip = pe_com_vec(x);
% Jphip = Je_com_mat(x);
% dJphip = dJe_com_mat(x);

phip = delta_phip_sca(x);
Jphip = Jdelta_phip_sca(x);
dJphip = zeros(1,nb+ndof);
dtdp = 1/vhip;
tau   = (phip(1) - p(1))/vhip;


%% COMZ: 0.9072 Total Mass:79.7761


if use_cwf
    yd = cwf_time(a,tau);
    dydt = dcwf_time(a,tau);
    ddydt = d2cwf_time(a,tau);
else
    yd = bezier(a,tau);
    dydt = dbezier(a,tau);
    ddydt = d2bezier(a,tau);
end


yd(1) = com_onestep_time();
for i = 1:no
    Jyd(i, :) = [dydt(i) * Jphip(1,:)*dtdp,zeros(1,nb+ndof)];
    dJyd(i, :) = [ddydt(i) * ((Jphip(1,:)*dtdp)*dqe) * (Jphip(1,:)*dtdp) + ...
        dydt(i) * dJphip(1,:),...
        dydt(i) * (Jphip(1,:)*dtdp)];
    
end

end