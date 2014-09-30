close all

load('figs/simulation_data');
nSteps = size(calcs_set,1)/3;
nb = 2;
ndof = 7;
% slip_data = get_slip_data(nSteps);

labels = {'sa', 'sk', 'sh', 'nsh', 'nsk','nsa'};

% rspring_indices = [2 9];
rigid_indices = [1 :6];
act_indices = [1:6];

% Es = out.T + out.V;

qbs = out.qe(1:nb, :);
qrs = out.qe(nb + (1:ndof), :);
dqrs = out.dqe(nb + (1:ndof), :);
ddqrs = out.ddqe(nb + (1:ndof), :);

% Torques
h = figure(1);
clf();
% set(h,'position',[100,100,480,360]);
% % subplot(2,1,1);
% plot(out.t, out.us);
% legend(labels(spring_indices));
% title('Spring Force (Stance / Nonstance)');
% xlabel('t [s]');
% ylabel('u [N]');
% subplot(2,1,2);
plot(out.t, out.ur);
legend(labels(act_indices));
title('Joint Torque (Stance / Nonstance)');
xlabel('t [s]');
ylabel('u [N m]');
print(h,'-dpng','-r100','figs/torques');


% angles
h =figure(2);
clf();
% hold('all');
% subplot(2,1,1);
% plot(out.t, qrs(spring_indices,:));
% legend(labels(spring_indices));
% title('Spring Deflection (Stance / Nonstance)');
% xlabel('t [s]');
% ylabel('q [rad]');
% subplot(2,1,2);
plot(out.t, qrs(rigid_indices,:));
legend(labels(rigid_indices));
title('Joint Position (Stance / Nonstance)');
xlabel('t [s]');
ylabel('q [rad]');
print(h,'-dpng','-r100','figs/angles');


% velocities
h =figure(3);
clf();
% hold('all');
% subplot(2,1,1);
% plot(out.t, dqrs(spring_indices,:));
% legend(labels(spring_indices));
% title('Spring Velocity (Stance / Nonstance)');
% xlabel('t [s]');
% ylabel('dq [rad/s]');
% subplot(2,1,2);
plot(out.t, dqrs(rigid_indices,:));
legend(labels(rigid_indices));
title('Joint Velocity (Stance / Nonstance)');
xlabel('t [s]');
ylabel('dq [rad/s]');
print(h,'-dpng','-r100','figs/velocities');


% accelerations

h =figure(30);
clf();
% hold('all');
% subplot(2,1,1);
% plot(out.t, ddqrs(spring_indices,:));
% legend(labels(spring_indices));
% title('Spring Accelerations (Stance / Nonstance)');
% xlabel('t [s]');
% ylabel('ddq [rad/s^2]');
% subplot(2,1,2);
plot(out.t, ddqrs(rigid_indices,:));
legend(labels(rigid_indices));
title('Joint Accelerations (Stance / Nonstance)');
xlabel('t [s]');
ylabel('ddq [rad/s^2]');
print(h,'-dpng','-r100','figs/accelerations');


% actual energy
    % h = figure(4);
    % clf();
    % hold('all');
    % plot(out.t, out.T);
    % plot(out.t, out.V);
    % plot(out.t, Es);
    % legend({'Kinetic', 'Potential', 'Total'});
    % xlabel('t [s]');
    % ylabel('E [J]');
    % print(h,'-dpng','-r100','figs/energies');

% output comparison

h =figure(4); 
clf();
set(h,'position',[100,100,1024,768]);
hold('all');

% subplot(2,4,1);
% plot(out.t,out.ya1,out.t,out.yd1,'r--');
% xlim([out.t(1) out.t(end)]);
% title('Non-stance Spring');
% xlabel('t [s]');
% ylabel('r [m]');
subplot(2,3,1);
plot(out.t,out.ya(1,:),out.t,out.yd(1,:),'r--');
xlim([out.t(1) out.t(end)]);
title('Com X');
xlabel('t [s]');
ylabel('r [m]');
subplot(2,3,2);
plot(out.t,out.ya(2,:),out.t,out.yd(2,:),'r--');
xlim([out.t(1) out.t(end)]);
title('Com Z');
xlabel('t [s]');
ylabel('r [m]');
subplot(2,3,3);
plot(out.t,out.ya(3,:),out.t,out.yd(3,:),'r--');
xlim([out.t(1) out.t(end)]);
title('Torso Angle');
xlabel('t [s]');
ylabel('q [rad]');
subplot(2,3,4);
plot(out.t,out.ya(4,:),out.t,out.yd(4,:),'r--');
xlim([out.t(1) out.t(end)]);
title('Swing Foot X');
xlabel('t [s]');
ylabel('r [m]');
subplot(2,3,5);
plot(out.t,out.ya(5,:),out.t,out.yd(5,:),'r--');
xlim([out.t(1) out.t(end)]);
title('Swing Foot Z');
xlabel('t [s]');
ylabel('r [m]');
subplot(2,3,6);
plot(out.t,out.ya(6,:),out.t,out.yd(6,:),'r--');
xlim([out.t(1) out.t(end)]);
title('Swing Foot Angle');
xlabel('t [s]');
ylabel('q [rad]');
% subplot(2,4,7);
% plot(out.t,out.ya(7,:),out.t,out.yd(7,:),'r--');
% xlim([out.t(1) out.t(end)]);
% title('Non-stance Spring');
% xlabel('t [s]');
% ylabel('q [rad]');
 print(h,'-dpng','-r100','figs/outputs');


 
% output comparison error terms

h =figure(5); 
clf();
set(h,'position',[100,100,1024,768]);
hold('all');

% subplot(2,4,1);
% plot(out.t,out.ya1,out.t,out.yd1,'r--');
% xlim([out.t(1) out.t(end)]);
% title('Non-stance Spring');
% xlabel('t [s]');
% ylabel('r [m]');
subplot(2,3,1);
plot(out.t,out.ya(1,:)-out.yd(1,:));
xlim([out.t(1) out.t(end)]);
title('Com X');
xlabel('t [s]');
ylabel('r [m]');
subplot(2,3,2);
plot(out.t,out.ya(2,:)-out.yd(2,:));
xlim([out.t(1) out.t(end)]);
title('Com Z');
xlabel('t [s]');
ylabel('r [m]');
subplot(2,3,3);
plot(out.t,out.ya(3,:)-out.yd(3,:));
xlim([out.t(1) out.t(end)]);
title('Torso Angle');
xlabel('t [s]');
ylabel('q [rad]');
subplot(2,3,4);
plot(out.t,out.ya(4,:)-out.yd(4,:));
xlim([out.t(1) out.t(end)]);
title('Swing Foot X');
xlabel('t [s]');
ylabel('r [m]');
subplot(2,3,5);
plot(out.t,out.ya(5,:)-out.yd(5,:));
xlim([out.t(1) out.t(end)]);
title('Swing Foot Z');
xlabel('t [s]');
ylabel('r [m]');
subplot(2,3,6);
plot(out.t,out.ya(6,:)-out.yd(6,:));
xlim([out.t(1) out.t(end)]);
title('Swing Foot Angle');
xlabel('t [s]');
ylabel('q [rad]');
% subplot(2,4,7);
% plot(out.t,out.ya(7,:),out.t,out.yd(7,:),'r--');
% xlim([out.t(1) out.t(end)]);
% title('Non-stance Spring');
% xlabel('t [s]');
% ylabel('q [rad]');
 print(h,'-dpng','-r100','figs/outputs'); 
 

%%
h =figure(6); clf;
hold('all');
xlim([out.t(1) out.t(end)]);
% plot(out.t,out.pcom(1,:)); hold on;

% plot((-calcs.Fe(3,:)./calcs.Fe(2,:))')
% ylim('auto');


% line('XData', [0 size(out.t,2)], 'YData', [0.014-0.2 0.014-0.2], 'LineWidth', 2, ...
%     'LineStyle', '-.', 'Color', [1 0.6 0]);
% line('XData', [0 size(out.t,2)], 'YData', [-0.189 -0.189], 'LineWidth', 2, ...
%     'LineStyle', '-.', 'Color', [0 0.6 0]);
plot(out.t, out.ZMP); 
% legend('UpperBound', 'LowerBound' ,'Zmp X');
legend('Zmp X');
xlabel('t [s]');
ylabel('ZMPx [m]');
title('ZMP x'); 
print(h,'-dpng','-r100','figs/CoM_position');

%%
h =figure(7); clf;
hold('all');
xlim([out.t(1) out.t(end)]);
% plot(out.t,out.pcom(1,:)); hold on;

% plot((-calcs.Fe(3,:)./calcs.Fe(2,:))')
% ylim('auto');


% line('XData', [0 size(out.t,2)], 'YData', [0.014-0.2 0.014-0.2], 'LineWidth', 2, ...
%     'LineStyle', '-.', 'Color', [1 0.6 0]);
% line('XData', [0 size(out.t,2)], 'YData', [-0.189 -0.189], 'LineWidth', 2, ...
%     'LineStyle', '-.', 'Color', [0 0.6 0]);
plot(out.t, out.Fe); 
% legend('UpperBound', 'LowerBound' ,'Zmp X');
legend('Fsx','Fsz', 'Tsy', 'Fnsx','Fnsz', 'Tnsy');
xlabel('t [s]');
ylabel('Force[N] or Torque [Nm]');
title('Ground Reaction Force'); 
print(h,'-dpng','-r100','figs/CoM_position');

%%

 % hip position
% h =figure(7); clf;
% subplot(2,1,1);
% plot(out.t,out.pcom(1,:)); hold on;
% plot(slip_data.ts, slip_data.xcomc,'r--'); 
% legend('Proxi','SLIP');
% xlabel('t [s]');
% ylabel('p [m]');
% title('Horizontal Position');
% subplot(2,1,2);
% plot(out.t,out.pcom(2,:)); hold on;
% plot(slip_data.ts, slip_data.ycom,'r--');
% legend('Proxi','SLIP');
% xlabel('t [s]');
% ylabel('p [m]');
% title('Vertical Position');
% print(h,'-dpng','-r100','figs/CoM_position');
% 
% % hip position
% h =figure(8);clf;
% subplot(2,1,1);
% plot(out.t,out.pcomdot(1,:)); hold on;
% plot(slip_data.ts, slip_data.dxcom,'r--'); 
% legend('Proxi','SLIP');
% xlabel('t [s]');
% ylabel('v [m/s]');
% title('Horizontal Velocity');
% subplot(2,1,2);
% plot(out.t,out.pcomdot(2,:)); hold on;
% plot(slip_data.ts, slip_data.dycom,'r--');
% legend('Proxi','SLIP');
% xlabel('t [s]');
% ylabel('v [m/s]');
% title('Vertical Velocity');
% print(h,'-dpng','-r100','figs/CoM_velocity');
% 
% h =figure(9);clf;
% subplot(2,1,1);
% plot(out.pcom(1,:),out.pcomdot(1,:)); hold on;
% plot(slip_data.xcomc, slip_data.dxcom,'r--'); 
% legend('Proxi','SLIP');
% xlabel('t [s]');
% ylabel('v [m/s]');
% title('Horizontal Velocity');
% subplot(2,1,2);
% plot(out.pcom(2,:),out.pcomdot(2,:)); hold on;
% plot(slip_data.ycom, slip_data.dycom,'r--');
% legend('Proxi','SLIP');
% xlabel('t [s]');
% ylabel('v [m/s]');
% title('Vertical Velocity');
% print(h,'-dpng','-r100','figs/limit_cycles');
% 
% % energy
% h =figure(10); clf;
% plot(out.t,out.Energy); hold on;
% plot(out.t,Energy0*ones(size(out.t)));
% legend('Actual','Desired');
% xlabel('t [s]');
% ylabel('E [J]');
% title('SLIP-like Energy')
% print(h,'-dpng','-r100','figs/energy');
% 
% 
% % ground reaction forces
% h =figure(11); clf;
% subplot(2,1,1);
% plot(out.t,out.Fr(1,:)); hold on;
% plot(slip_data.ts,slip_data.StanceForceX,'r--');
% legend('Proxi','SLIP');
% xlabel('t [s]');
% ylabel('F [N]');
% title('Horizontal Force');
% subplot(2,1,2);
% plot(out.t,out.Fr(2,:));
% hold on;
% plot(slip_data.ts,slip_data.StanceForceY,'r--');
% legend('Proxi','SLIP');
% xlabel('t [s]');
% ylabel('F [N]');
% title('Vertical Force');
% print(h,'-dpng','-r100','figs/reaction_forces_left');
% 
% 
% 
% h =figure(12); clf;
% subplot(2,1,1);
% plot(out.t,out.Fr(3,:)); hold on;
% plot(slip_data.ts,slip_data.NonStanceForceX,'r--');
% legend('Proxi','SLIP');
% xlabel('t [s]');
% ylabel('F [N]');
% title('Horizontal Force');
% subplot(2,1,2);
% plot(out.t,out.Fr(4,:));
% hold on;
% plot(slip_data.ts,slip_data.NonStanceForceY,'r--');
% legend('Proxi','SLIP');
% xlabel('t [s]');
% ylabel('F [N]');
% title('Vertical Force');
% print(h,'-dpng','-r100','figs/reaction_forces_right');
% 
% 
% 
% h =figure(13); clf;
% plot(out.qs(1,:),out.us(1,:),'*');
% xlabel('delta [m]');
% ylabel('force [N]');
% title('Spring deflection vs. Spring Force (Left)');
% print(h,'-dpng','-r100','figs/spring_deflection_vs_force_left');
% 
% 
% h =figure(14);clf;
% plot(out.qs(2,:),out.us(2,:),'*');
% 
% xlabel('delta [m]');
% ylabel('force [N]');
% title('Spring deflection vs. Spring Force (Right)');
% print(h,'-dpng','-r100','figs/spring_deflection_vs_force_right');
% 
% 
% h = figure(15); clf;
% subplot(2,1,1);
% plot(out.t,out.ur(1,:));
% subplot(2,1,2);
% plot(out.t,out.dqe(4,:));


% output comparison

% h =figure(4); clf();
% set(h,'position',[100,100,1024,768]);
% hold('all');
