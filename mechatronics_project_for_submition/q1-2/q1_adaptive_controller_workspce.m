%% system constans
L = 1;
m = 1;
k = 1;
c = 1;
J = 2*L^2/3*m;
g = 9.8;

alfa_1 = L*m*g;
alfa_2 = k*L^2;
alfa_3 = c*L^2;
%% controller constans
Ts = 0.001;
K = 10;
lamda_0 = 5;
gama =5;
%% run simulation
out = sim('q1_adaptive_controller_simulink');
%% make graphs
% Angle V.S Time
figure(1)
set(gcf,'color','w');
hold on;
plot(out.r.Time, out.r.Data(1,:));
plot(out.theta.Time, out.theta.Data);
grid on;
title('Angle V.S Time');
legend('Theta design', 'Theta');
xlabel('Time [sec]');
ylabel('Angle [rad]');
theta_q2 = out.theta.Data;
%% Eror V.S Time
figure(2)
set(gcf,'color','w');
hold on;
plot(out.e.Time, out.e.Data);
grid on;
title('Eror V.S Time');
xlabel('Time [sec]');
ylabel('Eror [rad]');
%% Z V.S Time
figure(3)
set(gcf,'color','w');
set(gcf,'color','w');
hold on;
plot(out.Z.Time, out.Z.Data);
grid on;
title('Z V.S Time');
xlabel('Time [sec]');
ylabel('Z [rad/sec]');
%% J-Estimation V.S Time
figure(4)
set(gcf,'color','w');
hold on;
plot(out.J_est.Time, J*ones(1, length(out.J_est.Data)));
plot(out.J_est.Time, out.J_est.Data)
grid on;
title('J Estimation V.S Time');
xlabel('Time [sec]');
ylabel('Moment of Inertia [N*m]');
legend('J Real', 'J Estimation');
%% Alfa_1 Estimation V.S Time
figure(5)
set(gcf,'color','w');
hold on;
plot(out.alfa_1_est.Time, alfa_1*ones(1, length(out.alfa_1_est.Data)));
plot(out.alfa_1_est.Time, out.alfa_1_est.Data)
grid on;
title('Alfa_1 Estimation V.S Time');
xlabel('Time [sec]');
ylabel('Alfa_1 [N*m]');
legend('Alfa_1 Real', 'Alfa_1 Estimation');
%% Alfa_2 Estimation V.S Time
figure(6)
set(gcf,'color','w');
hold on;
plot(out.alfa_2_est.Time, alfa_2*ones(1, length(out.alfa_2_est.Data)));
plot(out.alfa_2_est.Time, out.alfa_2_est.Data)
grid on;
title('Alfa_2 Estimation V.S Time');
xlabel('Time [sec]');
ylabel('Alfa_2 [N*m]');
legend('Alfa_2 Real', 'Alfa_2 Estimation');
%% Alfa_3 Estimation V.S Time
figure(7)
set(gcf,'color','w');
hold on;
plot(out.alfa_3_est.Time, alfa_3*ones(1, length(out.alfa_3_est.Data)));
plot(out.alfa_3_est.Time, out.alfa_3_est.Data)
grid on;
title('Alfa_3 Estimation V.S Time');
xlabel('Time [sec]');
ylabel('Alfa_3 [N*m]');
legend('Alfa_3 Real', 'Alfa_3 Estimation');
%% u(t) V.S Time
figure(8)
set(gcf,'color','w');
hold on;
plot(out.u.Time, out.u.Data)
grid on;
title('Controller Signal V.S Time');
xlabel('Time [sec]');
ylabel('u(t) [N*m]');