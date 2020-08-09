%% configure plant
%system constans
L = 1;
m = 1;
k = 1;
c = 1;
J = 2*(L^2/3)*m;
g = 9.8;

alfa_1 = L*m*g/J;
alfa_2 = k*L^2/J;
alfa_3 = c*L^2/J;
%% init RBF
%build NN
n = 5; %viden layer size
%domaim of angle and velocity
theta_max = pi;
theta_min = -pi;
theta_dot_max = 2;
theta_dot_min = -2;
delta_theta = (theta_max - theta_min)/n;
delta_theta_dot = (theta_dot_max - theta_dot_min)/n;
%difdelta_thetaine the rbf
%c = [linspace(theta_min, theta_max, n) ; linspace(theta_dot_min, theta_dot_max, n)]; 
c1 = linspace(theta_min, theta_max ,n);
c2 = linspace(theta_dot_min, theta_dot_max, n);
[C1,C2] = meshgrid(c1,c2);
c = [C1(:),C2(:)]';
b = (theta_max - theta_min)/n;
W0 = (rand(1, n^2) - 0.5)/10;%initial waits
n = n^2;
%% define the adaptation law
gama = 40;
eta = 0.1;
m_min = 0.5; %the inverse of the maximum evaluation of J
%% define the contriller
Ts = 0.001;
kp = 100;
kd = 10;
K = [kp;kd];
%% build the eror state
Fai = [0 1; -kp -kd];
A = Fai';
B = [0; 1];
%% solve lyaponov equsaion
Q = [500 0; 0 500];
P = lyap(A, Q);
%% run simulation
out = sim('q3_rbf_adaptive_controller_sim_at2019a');
%% build graphics 
%%Angle V.S Time
figure(1)
set(gcf,'color','w');
hold on;
plot(out.r.Time, out.r.Data(1,:), 'b-.', 'LineWidth', 1.25);
plot(out.theta.Time, out.theta.Data, 'r-');
grid on;
title('Angle V.S Time');
legend('Theta design', 'Theta');
xlabel('Time [sec]');
ylabel('Angle [rad]');
%% Eror V.S Time
figure(2)
hold on;
set(gcf,'color','w');
plot(out.e.Time, out.e.Data);
grid on;
title('Eror V.S Time');
xlabel('Time [sec]');
ylabel('Eror [rad]');
%% f real V.S f est
figure(3)
set(gcf,'color','w');
hold on;
plot(out.f_real.Time, out.f_real.Data);
plot(out.f_est.Time, out.f_est.Data);
grid on;
title('f real V.S f est');
legend('f real', 'f estimation');
xlabel('Time [sec]');
ylabel('f [N/kg*m]');
%% J-Estimation V.S Time
figure(4)
set(gcf,'color','w');
hold on;
plot(out.m_est.Time, J*ones(1, length(out.m_est.Time)));
plot(out.m_est.Time, 1/out.m_est.Data);
grid on;
title('J Estimation V.S Time');
xlabel('Time [sec]');
ylabel('Moment of Inertia [N*m]');
legend('J Real', 'J Estimation');
%% u(t) V.S Time
figure(5)
set(gcf,'color','w');
hold on;
plot(out.u.Time, out.u.Data)
grid on;
title('Controller Signal V.S Time');
xlabel('Time [sec]');
ylabel('u(t) [N*m]');
%% W V.S Time
figure(6)
set(gcf,'color','w');
hold on;
plot(out.W.Time, out.W.Data(:,1));
plot(out.W.Time, out.W.Data(:,2));
plot(out.W.Time, out.W.Data(:,3));
plot(out.W.Time, out.W.Data(:,4));
plot(out.W.Time, out.W.Data(:,5));
grid on;
title('W V.S Time');
xlabel('Time [sec]');
ylabel('W [-]');
legend('W_1', 'W_2', 'W_3', 'W_4', 'W_5');