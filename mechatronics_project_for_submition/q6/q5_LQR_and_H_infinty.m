%% build system
%syms k_12 k_23 m1 m2 m3 gama1 gama2;
k_12 = 800;
k_23 = 100;
m1 = 2;
m2 = 1.5;
m3 = 1;
T = 20;
%% build matrix
A = [0 0 0 1 0 0;
    0 0 0 0 1 0;
    0 0 0 0 0 1;
    -k_12/m1 k_12/m1 0 0 0 0;
    k_12/m2 -(k_12+k_23)/m2 k_23/m2 0 0 0;
    0 k_23/m3 -k_23/m3 0 0 0];
B1 = [0 0 0 1/m1 0 0]';
B2 = [0 0 0 0 1/m2 0]';
B = [B1, B2];
C = [0 0 1 0 0 0];
D = [0,0];
%% check if the system is contrillabilit
S1 = [B1 A*B1 (A.^2)*B1 (A.^3)*B1 (A.^4)*B1 (A.^5)*B1];
S2 = [B2 A*B2 (A.^2)*B2 (A.^3)*B2 (A.^4)*B2 (A.^5)*B2];
S = [B A*B (A.^2)*B (A.^3)*B (A.^4)*B (A.^5)*B];
ranks_LQR_1 = rank(S1);
ranks_LQR_2 = rank(S2);
ranks_LQR = rank(S);
eigen_open = eig(A);

%% build contriller
gama1 = 0.0001;
gama2 = 0.0001;
Q = C'*C;
R = [gama1 0; 0 gama2];
[K] = lqr(A, B, Q, R);
eigen_LQR_close = eig(A - B*K);
rank_LQR_close = rank(A - B*K);
%% H-infinity
B2h = [0 0 0 1/m1 0 0]';
B1h = [0 0 0 0 1/m2 0]';
Bh = [B1h, B2h];
gama2_h = 5;
gama1_h = 0.1;
mh1 = size(B1h,6);
mh2 = size(B2h,6);
Rh = [-gama2_h^2*eye(mh1) zeros(mh1,mh2) ; zeros(mh2,mh1) eye(mh2)];
Ph = care(A,Bh,C'*C,Rh);
Kh = inv(gama1_h)*B2h'*Ph;
eigen_Hinf_close = eig(A - B2h*Kh);
rank_Hinf_close = rank(A - B2h*Kh);
%% run simulation of LQR controller

gama1 = 0.01;
gama2 = 0.01;
R = [gama1 0; 0 gama2];
[K] = lqr(A, B, Q, R);
LQR_G1_001_G2_001 = sim('q5_LQR_sim');
%% LQR for divergent gamma 2
gama1 = 0.01;
gama2 = 0.1;
R = [gama1 0; 0 gama2];
[K] = lqr(A, B, Q, R);
LQR_G1_001_G2_01 = sim('q5_LQR_sim');
%% LQR for divergent gamma 2
gama1 = 0.01;
gama2 = 0.001;
R = [gama1 0; 0 gama2];
[K] = lqr(A, B, Q, R);
LQR_G1_001_G2_0001 = sim('q5_LQR_sim');
%% LQR for divergent gamma 2
gama1 = 0.01;
gama2 = 0.0001;
R = [gama1 0; 0 gama2];
[K] = lqr(A, B, Q, R);
LQR_G1_001_G2_00001 = sim('q5_LQR_sim');
%% LQR for divergent gamma 2
gama1 = 0.01;
gama2 = 1;
R = [gama1 0; 0 gama2];
[K] = lqr(A, B, Q, R);
LQR_G1_001_G2_1 = sim('q5_LQR_sim');
%% LQR for divergent gamma 1
gama1 = 0.001;
gama2 = 0.01;
R = [gama1 0; 0 gama2];
[K] = lqr(A, B, Q, R);
LQR_G1_0001_G2_001 = sim('q5_LQR_sim');
%% LQR for divergent gamma 1
gama1 = 1;
gama2 = 0.01;
R = [gama1 0; 0 gama2];
[K] = lqr(A, B, Q, R);
LQR_G1_1_G2_001 = sim('q5_LQR_sim');
%% LQR for divergent gamma 1
gama1 = 0.0001;
gama2 = 0.01;
R = [gama1 0; 0 gama2];
[K] = lqr(A, B, Q, R);
LQR_G1_00001_G2_001 = sim('q5_LQR_sim');
%% LQR for divergent gamma 1
gama1 = 0.1;
gama2 = 0.01;
R = [gama1 0; 0 gama2];
[K] = lqr(A, B, Q, R);
LQR_G1_01_G2_001 = sim('q5_LQR_sim');
%% gamma1 = 1 and gamma2  =1
gama1 = 1;
gama2 = 1;
R = [gama1 0; 0 gama2];
[K] = lqr(A, B, Q, R);
LQR_G1_1_G2_1 = sim('q5_LQR_sim');
%% makes graphs for LQR controller

%% LQR y V.S Time - gamma 2 divergent
figure(1)
set(gcf,'color','w');
hold on;
plot(LQR_G1_001_G2_01.y.Time, LQR_G1_001_G2_01.y.Data, ...
LQR_G1_001_G2_0001.y.Time, LQR_G1_001_G2_0001.y.Data, ...
LQR_G1_001_G2_001.y.Time, LQR_G1_001_G2_001.y.Data, ...
LQR_G1_001_G2_1.y.Time, LQR_G1_001_G2_1.y.Data, ...
LQR_G1_001_G2_00001.y.Time, LQR_G1_001_G2_00001.y.Data);
grid on;
xlabel("Time [sec]"); ylabel("Y");
title("LQR - \gamma1 = 0.01 , \gamma2 = 0.0001 - 1");
legend('\gamma2 = 1', '\gamma2 = 0.1', '\gamma2 = 0.01', '\gamma2 = 0.001', '\gamma2 = 0.0001');
axis([0 5 -1 2.5]);
%% y V.S Time - gamma 1 divergent
figure(2)
set(gcf,'color','w');
hold on;
plot(LQR_G1_1_G2_001.y.Time, LQR_G1_1_G2_001.y.Data, ...
LQR_G1_01_G2_001.y.Time, LQR_G1_01_G2_001.y.Data, ...
LQR_G1_001_G2_001.y.Time, LQR_G1_001_G2_001.y.Data, ...
LQR_G1_0001_G2_001.y.Time, LQR_G1_0001_G2_001.y.Data, ...
LQR_G1_00001_G2_001.y.Time, LQR_G1_00001_G2_001.y.Data);
grid on;
xlabel("Time [sec]"); ylabel("Y");
title("LQR - \gamma1 = (0.0001 - 1) , \gamma2 = 0.01");
legend('\gamma1 = 1', '\gamma1 = 0.1', '\gamma1 = 0.01', '\gamma1 = 0.001', '\gamma1 = 0.0001');
axis([0 5 -1 2.5]);
%% state space V.S Time - big gammas
figure(3)
set(gcf,'color','w');
subplot(2,1,1)
hold on;
plot(LQR_G1_1_G2_1.x.Time, LQR_G1_1_G2_1.x.Data(:,1),...
    LQR_G1_1_G2_1.x.Time, LQR_G1_1_G2_1.x.Data(:,2));
plot(LQR_G1_1_G2_1.x.Time, LQR_G1_1_G2_1.x.Data(:,3), 'k-', 'LineWidth', 1.25);
grid on;
axis([0 5 min([min(LQR_G1_1_G2_1.x.Data(:,1)) min(LQR_G1_1_G2_1.x.Data(:,2)) min(LQR_G1_1_G2_1.x.Data(:,3))]) max([max(LQR_G1_1_G2_1.x.Data(:,1)) max(LQR_G1_1_G2_1.x.Data(:,2)) max(LQR_G1_1_G2_1.x.Data(:,3))])]);
xlabel("Time [sec]"); ylabel("state");
legend ('x1', 'x2', 'x3');
title("LQR - position states V.S Time  \gamma1 = 1 , \gamma2 = 1");

subplot(2,1,2)
hold on;
plot(LQR_G1_001_G2_001.x.Time, LQR_G1_001_G2_001.x.Data(:,1),...
    LQR_G1_001_G2_001.x.Time, LQR_G1_001_G2_001.x.Data(:,2));
plot(LQR_G1_001_G2_001.x.Time, LQR_G1_001_G2_001.x.Data(:,3), 'k-', 'LineWidth', 1.25);
grid on;
axis([0 5 min([min(LQR_G1_001_G2_001.x.Data(:,1)) min(LQR_G1_001_G2_001.x.Data(:,2)) min(LQR_G1_001_G2_001.x.Data(:,3))]) max([max(LQR_G1_001_G2_001.x.Data(:,1)) max(LQR_G1_001_G2_001.x.Data(:,2)) max(LQR_G1_001_G2_001.x.Data(:,3))])]);
xlabel("Time [sec]"); ylabel("state");
legend ('x1', 'x2', 'x3');
title("LQR - position states V.S Time \gamma1 = 0.01 , \gamma2 = 0.01");
%% LQR - the state response for different gamma
figure(4)
set(gcf,'color','w');
subplot(2,1,1)
hold on;
plot(LQR_G1_1_G2_1.x.Time, LQR_G1_1_G2_1.x.Data(:,4),...
    LQR_G1_1_G2_1.x.Time, LQR_G1_1_G2_1.x.Data(:,5));
plot(LQR_G1_1_G2_1.x.Time, LQR_G1_1_G2_1.x.Data(:,6), 'k-', 'LineWidth', 1.25);
grid on;
axis([0 5 min([min(LQR_G1_1_G2_1.x.Data(:,4)) min(LQR_G1_1_G2_1.x.Data(:,5)) min(LQR_G1_1_G2_1.x.Data(:,6))]) max([max(LQR_G1_1_G2_1.x.Data(:,4)) max(LQR_G1_1_G2_1.x.Data(:,5)) max(LQR_G1_1_G2_1.x.Data(:,6))])]);
xlabel("Time [sec]"); ylabel("state");
legend ('x1', 'x2', 'x3');
title("LQR - speed states V.S Time  \gamma1 = 1 , \gamma2 = 1");

subplot(2,1,2)
hold on;
plot(LQR_G1_001_G2_001.x.Time, LQR_G1_001_G2_001.x.Data(:,4),...
    LQR_G1_001_G2_001.x.Time, LQR_G1_001_G2_001.x.Data(:,5));
plot(LQR_G1_001_G2_001.x.Time, LQR_G1_001_G2_001.x.Data(:,6), 'k-', 'LineWidth', 1.25);
grid on;
axis([0 5 min([min(LQR_G1_001_G2_001.x.Data(:,4)) min(LQR_G1_001_G2_001.x.Data(:,5)) min(LQR_G1_001_G2_001.x.Data(:,6))]) max([max(LQR_G1_001_G2_001.x.Data(:,4)) max(LQR_G1_001_G2_001.x.Data(:,5)) max(LQR_G1_001_G2_001.x.Data(:,6))])]);
xlabel("Time [sec]"); ylabel("state");
legend ('x1', 'x2', 'x3');
title("LQR - speed states V.S Time \gamma1 = 0.01 , \gamma2 = 0.01");
%% controll signal VS states
figure(5)
set(gcf,'color','w');
subplot(2,1,1)
plot(LQR_G1_1_G2_001.u.Time, LQR_G1_1_G2_001.u.Data(:,1),...
    LQR_G1_1_G2_001.u.Time, LQR_G1_1_G2_001.u.Data(:,2));
grid on;
%axis([0 10 -25 25]);
xlabel("Time [sec]"); ylabel("control signal");
legend ('u1', 'u2');
title("LQR - control signal V.S Time  \gamma1 = 1 , \gamma2 = 0.01");

subplot(2,1,2)
plot(LQR_G1_001_G2_1.u.Time, LQR_G1_001_G2_1.u.Data(:,1),...
    LQR_G1_001_G2_1.u.Time, LQR_G1_001_G2_1.u.Data(:,2));
grid on;
%axis([0 10 -25 25]);
xlabel("Time [sec]"); ylabel("control signal");
legend ('u1', 'u2');
title("LQR - control signal V.S Time  \gamma1 = 0.01 , \gamma2 = 1");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% run simulation of H-inf controller
%% H - inf for divergent gamma 2 
gama2_h = 10;
gama1_h = 0.01;
Rh = [-gama2_h^2*eye(mh1) zeros(mh1,mh2) ; zeros(mh2,mh1) eye(mh2)];
Ph = care(A,Bh,C'*C,Rh);
Kh = inv(gama1_h)*B2h'*Ph;
Hinf_G1_001_G2_10 = sim('q6_H_inf_sim');
%% H - inf for divergent gamma 2 
gama2_h = 100;
gama1_h = 0.01;
Rh = [-gama2_h^2*eye(mh1) zeros(mh1,mh2) ; zeros(mh2,mh1) eye(mh2)];
Ph = care(A,Bh,C'*C,Rh);
Kh = inv(gama1_h)*B2h'*Ph;
Hinf_G1_001_G2_100 = sim('q6_H_inf_sim');
%% H - inf for divergent gamma 2 
gama2_h = 1.45;
gama1_h = 0.01;
Rh = [-gama2_h^2*eye(mh1) zeros(mh1,mh2) ; zeros(mh2,mh1) eye(mh2)];
Ph = care(A,Bh,C'*C,Rh);
Kh = inv(gama1_h)*B2h'*Ph;
Hinf_G1_001_G2_145 = sim('q6_H_inf_sim');
%% H - inf for divergent gamma 2 
gama2_h = 1000;
gama1_h = 0.01;
Rh = [-gama2_h^2*eye(mh1) zeros(mh1,mh2) ; zeros(mh2,mh1) eye(mh2)];
Ph = care(A,Bh,C'*C,Rh);
Kh = inv(gama1_h)*B2h'*Ph;
Hinf_G1_001_G2_1000 = sim('q6_H_inf_sim');
%% H - inf for divergent gamma 2 
gama2_h = 10000;
gama1_h = 0.01;
Rh = [-gama2_h^2*eye(mh1) zeros(mh1,mh2) ; zeros(mh2,mh1) eye(mh2)];
Ph = care(A,Bh,C'*C,Rh);
Kh = inv(gama1_h)*B2h'*Ph;
Hinf_G1_001_G2_10000 = sim('q6_H_inf_sim');
%% H - inf for divergent gamma 1 and gamma 2 is H-inf gain
gama2_h = 1.45;
gama1_h = 1;
Rh = [-gama2_h^2*eye(mh1) zeros(mh1,mh2) ; zeros(mh2,mh1) eye(mh2)];
Ph = care(A,Bh,C'*C,Rh);
Kh = inv(gama1_h)*B2h'*Ph;
Hinf_G1_1_G2_145 = sim('q6_H_inf_sim');
%% H - inf for divergent gamma 1 and gamma 2 is H-inf gain
gama2_h = 1.45;
gama1_h = 0.1;
Rh = [-gama2_h^2*eye(mh1) zeros(mh1,mh2) ; zeros(mh2,mh1) eye(mh2)];
Ph = care(A,Bh,C'*C,Rh);
Kh = inv(gama1_h)*B2h'*Ph;
Hinf_G1_01_G2_145 = sim('q6_H_inf_sim');
%% H - inf for divergent gamma 1 and gamma 2 is H-inf gain
gama2_h = 1.45;
gama1_h = 0.01;
Rh = [-gama2_h^2*eye(mh1) zeros(mh1,mh2) ; zeros(mh2,mh1) eye(mh2)];
Ph = care(A,Bh,C'*C,Rh);
Kh = inv(gama1_h)*B2h'*Ph;
Hinf_G1_001_G2_145 = sim('q6_H_inf_sim');
%% H - inf for divergent gamma 1 and gamma 2 is H-inf gain
gama2_h = 1.45;
gama1_h = 10;
Rh = [-gama2_h^2*eye(mh1) zeros(mh1,mh2) ; zeros(mh2,mh1) eye(mh2)];
Ph = care(A,Bh,C'*C,Rh);
Kh = inv(gama1_h)*B2h'*Ph;
Hinf_G1_10_G2_145 = sim('q6_H_inf_sim');
%% H - inf for divergent gamma 1 and gamma 2 is H-inf gain
gama2_h = 1.45;
gama1_h = 100;
Rh = [-gama2_h^2*eye(mh1) zeros(mh1,mh2) ; zeros(mh2,mh1) eye(mh2)];
Ph = care(A,Bh,C'*C,Rh);
Kh = inv(gama1_h)*B2h'*Ph;
Hinf_G1_100_G2_145 = sim('q6_H_inf_sim');
%% H - inf for divergent gamma 1 = 1 and gamma 2 = 1000
gama2_h = 1000;
gama1_h = 1;
Rh = [-gama2_h^2*eye(mh1) zeros(mh1,mh2) ; zeros(mh2,mh1) eye(mh2)];
Ph = care(A,Bh,C'*C,Rh);
Kh = inv(gama1_h)*B2h'*Ph;
Hinf_G1_1_G2_1000 = sim('q6_H_inf_sim');
%% H - inf for divergent gamma 1 = 1 and gamma 2 = 5
gama2_h = 5;
gama1_h = 1;
Rh = [-gama2_h^2*eye(mh1) zeros(mh1,mh2) ; zeros(mh2,mh1) eye(mh2)];
Ph = care(A,Bh,C'*C,Rh);
Kh = inv(gama1_h)*B2h'*Ph;
Hinf_G1_1_G2_5 = sim('q6_H_inf_sim');
%% H - inf for divergent gamma 1 = 0.01 and gamma 2 = 5
gama2_h = 5;
gama1_h = 0.01;
Rh = [-gama2_h^2*eye(mh1) zeros(mh1,mh2) ; zeros(mh2,mh1) eye(mh2)];
Ph = care(A,Bh,C'*C,Rh);
Kh = inv(gama1_h)*B2h'*Ph;
Hinf_G1_001_G2_5 = sim('q6_H_inf_sim');
%% H - inf for divergent gamma 1 = 0.1 and gamma 2 = 5
gama2_h = 5;
gama1_h = 0.1;
Rh = [-gama2_h^2*eye(mh1) zeros(mh1,mh2) ; zeros(mh2,mh1) eye(mh2)];
Ph = care(A,Bh,C'*C,Rh);
Kh = inv(gama1_h)*B2h'*Ph;
Hinf_G1_01_G2_5 = sim('q6_H_inf_sim');
%% H - inf for divergent gamma 1 = 100 and gamma 2 = 10
gama2_h = 10;
gama1_h = 10;
Rh = [-gama2_h^2*eye(mh1) zeros(mh1,mh2) ; zeros(mh2,mh1) eye(mh2)];
Ph = care(A,Bh,C'*C,Rh);
Kh = inv(gama1_h)*B2h'*Ph;
Hinf_G1_100_G2_10 = sim('q6_H_inf_sim');
%% makes graphs for H-inf controller

%% H-inf y V.S Time - gamma 2 divergent
figure(6)
set(gcf,'color','w');
hold on;
plot(Hinf_G1_001_G2_145.y.Time, Hinf_G1_001_G2_145.y.Data, ...
Hinf_G1_001_G2_10000.y.Time, Hinf_G1_001_G2_10000.y.Data);
grid on;
xlabel("Time [sec]"); ylabel("Y");
title("H-inf, y V.S Time for \gamma1 = 0.01 , \gamma2 = 1.45 and 1000");
legend('\gamma2 = 1.45', '\gamma2 = 10000');
%% H-inf y V.S Time - gamma 1 divergent
figure(7)
set(gcf,'color','w');
hold on;
subplot(2,1,1)
plot(Hinf_G1_001_G2_145.y.Time, Hinf_G1_001_G2_145.y.Data, ...
Hinf_G1_01_G2_145.y.Time, Hinf_G1_01_G2_145.y.Data,'k')
legend('\gamma1 = 0.01', '\gamma1 = 0.1');
title("H-inf, y for \gamma1 = 0.01', '\gamma1 = 0.1");
xlabel("Time [sec]"); ylabel("Y");
axis([0 20 min([min(Hinf_G1_001_G2_145.y.Data) min(Hinf_G1_01_G2_145.y.Data)]) max([max(Hinf_G1_001_G2_145.y.Data) max(Hinf_G1_01_G2_145.y.Data)])]);

grid on
title("H-inf, y for \gamma1 = 0.01 and 0.1 , \gamma2 = 1.45");
subplot(2,1,2)
plot(Hinf_G1_1_G2_145.y.Time, Hinf_G1_1_G2_145.y.Data, ...
Hinf_G1_10_G2_145.y.Time, Hinf_G1_10_G2_145.y.Data, ...
Hinf_G1_100_G2_145.y.Time, Hinf_G1_100_G2_145.y.Data);
legend('\gamma1 = 1', '\gamma1 = 10', '\gamma1 = 100');
title("H-inf, y for \gamma1 = 1, 10 and 100 , \gamma2 = 1.45");
grid on;
xlabel("Time [sec]"); ylabel("Y");
axis([0 20 min([min(Hinf_G1_1_G2_145.y.Data) min(Hinf_G1_10_G2_145.y.Data) min(Hinf_G1_100_G2_145.y.Data)]) max([max(Hinf_G1_1_G2_145.y.Data) max(Hinf_G1_10_G2_145.y.Data) max(Hinf_G1_100_G2_145.y.Data)])]);
%% state space - positions V.S Time - big gammas
figure(8)
set(gcf,'color','w');
subplot(2,1,1)
plot(Hinf_G1_1_G2_1000.x.Time, Hinf_G1_1_G2_1000.x.Data(:,1),...
    Hinf_G1_1_G2_1000.x.Time, Hinf_G1_1_G2_1000.x.Data(:,2),...
    Hinf_G1_1_G2_1000.x.Time, Hinf_G1_1_G2_1000.x.Data(:,3));
grid on;
%axis([0 10 -5 5]);
xlabel("Time [sec]"); ylabel("state");
legend ('x1', 'x2', 'x3');
title("H-inf, position states V.S Time  \gamma1 = 1 , \gamma2 = 1000");

subplot(2,1,2)
plot(Hinf_G1_1_G2_5.x.Time, Hinf_G1_1_G2_5.x.Data(:,1),...
    Hinf_G1_1_G2_5.x.Time, Hinf_G1_1_G2_5.x.Data(:,2),...
    Hinf_G1_1_G2_5.x.Time, Hinf_G1_1_G2_5.x.Data(:,3));
grid on;
%axis([0 10 -5 5]);
xlabel("Time [sec]"); ylabel("state");
title("H-inf, position states V.S Time  \gamma1 = 1 , \gamma2 = 5");
legend ('x1', 'x2', 'x3');
%% H-inf, state space - speeds V.S Time - big gammas
figure(9)
set(gcf,'color','w');
subplot(2,1,1)
plot(Hinf_G1_1_G2_1000.x.Time, Hinf_G1_1_G2_1000.x.Data(:,4),...
    Hinf_G1_1_G2_1000.x.Time, Hinf_G1_1_G2_1000.x.Data(:,5),...
    Hinf_G1_1_G2_1000.x.Time, Hinf_G1_1_G2_1000.x.Data(:,6));
grid on;
% axis([0 10 -5 5]);
xlabel("Time [sec]"); ylabel("state");
legend ('x4', 'x5', 'x6');
title("H-inf, speed states V.S Time  \gamma1 = 1 , \gamma2 = 5");

subplot(2,1,2)
plot(Hinf_G1_1_G2_5.x.Time, Hinf_G1_1_G2_5.x.Data(:,4),...
    Hinf_G1_1_G2_5.x.Time, Hinf_G1_1_G2_5.x.Data(:,5),...
    Hinf_G1_1_G2_5.x.Time, Hinf_G1_1_G2_5.x.Data(:,6));
grid on;
% axis([0 10 -5 5]);
xlabel("Time [sec]"); ylabel("state");
legend ('x4', 'x5', 'x6');
title("H-inf, speed states V.S Time \gamma1 = 1 , \gamma2 = 1000");
%% state space - positions V.S Time - gamma1 - small
figure(10)
set(gcf,'color','w');
subplot(2,1,1)
plot(Hinf_G1_001_G2_5.x.Time, Hinf_G1_001_G2_5.x.Data(:,1),...
    Hinf_G1_001_G2_5.x.Time, Hinf_G1_001_G2_5.x.Data(:,2),...
    Hinf_G1_001_G2_5.x.Time, Hinf_G1_001_G2_5.x.Data(:,3));
grid on;
% axis([0 10 -5 5]);
xlabel("Time [sec]"); ylabel("state");
legend ('x1', 'x2', 'y');
title("H-inf, position states V.S Time  \gamma1 = 0.01 , \gamma2 = 5")
subplot(2,1,2)
plot(Hinf_G1_01_G2_5.x.Time, Hinf_G1_01_G2_5.x.Data(:,1),...
    Hinf_G1_01_G2_5.x.Time, Hinf_G1_01_G2_5.x.Data(:,2),...
    Hinf_G1_01_G2_5.x.Time, Hinf_G1_01_G2_5.x.Data(:,3));
grid on;
% axis([0 10 -5 5]);
xlabel("Time [sec]"); ylabel("state");
% legend ('x1', 'x2', 'x3');
title("H-inf, position states V.S Time \gamma1 = 0.1 , \gamma2 = 5");
%% H-inf, contntrrol signal
figure(11)
plot(Hinf_G1_1_G2_1000.u.Time, Hinf_G1_1_G2_1000.u.Data,...
    Hinf_G1_1_G2_5.u.Time, Hinf_G1_1_G2_5.u.Data);
grid on;
%axis([0 10 -10 10]);
xlabel("Time [sec]"); ylabel("control signal");
legend ('\gamma2 = 1000 ', '\gamma2 = 5');
title("control signal V.S Time  \gamma1 = 1");
%% H-inf, contntrrol signal
figure(12)
plot(Hinf_G1_001_G2_10.u.Time, Hinf_G1_001_G2_10.u.Data,...
    Hinf_G1_100_G2_10.u.Time, Hinf_G1_100_G2_10.u.Data);
grid on;
%axis([0 100 -40 40]);
xlabel("Time [sec]"); ylabel("control signal");
legend ('\gamma1 = 0.01 ', '\gamma1 = 100');
title("control signal V.S Time  \gamma2 = 5");
