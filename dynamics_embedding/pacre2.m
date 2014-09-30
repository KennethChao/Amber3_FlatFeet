function ret = pacre2(x,domains,Energy0,ref)

ndomains = size(domains,2);
ndof = 9;
nb = 2;

q_minus = x(1:ndof);
dq_minus = x(1+ndof:end);
xb_minus = [0;0;q_minus;0;0;dq_minus];


x0 = apply_reset(xb_minus,domains{end});

% [PE,KE] = robot.atriasEnergy(model,x0);
% Energy0 = PE+KE;

t_start = 0;
for j = 1:ndomains
        %         t_end = t_start + Tp(j);
        ref.h.calcs = {};
        odeopts = odeset('MaxStep', 1e-2, 'OutputFcn', @outputfcn, ...
            'Events', @(t, x, ref) eventfcn(t, x,domains{j}));
        
        
        sol = ode45(@(t, x,ref) calcEOM(t,x,domains{j},Energy0,ref,1), [t_start, 2], x0, ...
            odeopts, ref);
        t_start = sol.x(end);
        xb_minus = sol.y(:,end);
        x0 = apply_reset(xb_minus,domains{end});
        
        
end

q_minus = xb_minus(nb+1:nb+ndof);
dq_minus = xb_minus(nb+ndof+1:end);
ret = [q_minus(1:end);dq_minus];
    

end