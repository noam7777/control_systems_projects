[K_opt, K_ADP, A, B, Ts, C]=ADP();
%% run simulation in simulink and the input is dirak delta function
T = 10;
k_sim = K_opt;
x0 = [0 0 0 0 0 0];
out_optimal = sim('q7_LQR_based_data_set_sim');
k_sim = K_ADP;
out_ADP = sim('q7_LQR_based_data_set_sim');
%% y for K optimal and K ADP 
figure(3)
set(gcf,'color','w');
hold on;
plot(out_optimal.y.Time, out_optimal.y.Data, 'r-');
plot(out_ADP.y.Time, out_ADP.y.Data,'b-.', 'LineWidth', 1.25);
grid on;
title('y for K optimal and K ADP');
legend('K optimal', 'K ADP');
xlabel('Time [sec]'); ylabel('y [m]');
%% control signal for K optimal and K ADP
figure(4)
set(gcf,'color','w');
hold on;
plot(out_optimal.u.Time, out_optimal.u.Data(:,1));
plot(out_optimal.u.Time, out_optimal.u.Data(:,2));
plot(out_ADP.u.Time, out_ADP.u.Data(:,1),'b-.', 'LineWidth', 1.25);
plot(out_ADP.u.Time, out_ADP.u.Data(:,2),'g-.', 'LineWidth', 1.25);
grid on;
title(' control signal for K optimal and K ADP');
legend('u(1) optimal', 'u(2) optimal', 'u(1) ADP', 'u(2) ADP');
xlabel('Time [sec]'); ylabel('u [N]');

%% state parameter with optimal controller and optimal ADP
figure(5)
set(gcf,'color','w');
subplot(3,2,1)
hold on;
plot(out_optimal.x.Time, out_optimal.x.Data(:,1), 'r-');
plot(out_ADP.x.Time, out_ADP.x.Data(:,1),'b-.', 'LineWidth', 1.25);
grid on;
legend('K optimal', 'K ADP');
title(' X1 optmal V.S X1 ADP');
xlabel('Time [sec]'); ylabel('X_1 [m]');

subplot(3,2,2)
hold on;
plot(out_optimal.x.Time, out_optimal.x.Data(:,2), 'r-');
plot(out_ADP.x.Time, out_ADP.x.Data(:,2),'b-.', 'LineWidth', 1.25);
grid on;
legend('K optimal', 'K ADP');
title(' X2 optmal V.S X2 ADP');
xlabel('Time [sec]'); ylabel('X_2 [m]');

subplot(3,2,3)
hold on;
plot(out_optimal.x.Time, out_optimal.x.Data(:,3), 'r-');
plot(out_ADP.x.Time, out_ADP.x.Data(:,3),'b-.', 'LineWidth', 1.25);
grid on;
legend('K optimal', 'K ADP');
title(' X3 optmal V.S X3 ADP');
xlabel('Time [sec]'); ylabel('X_3 [m]');

subplot(3,2,4)
hold on;
plot(out_optimal.x.Time, out_optimal.x.Data(:,4), 'r-');
plot(out_ADP.x.Time, out_ADP.x.Data(:,4),'b-.', 'LineWidth', 1.25);
grid on;
legend('K optimal', 'K ADP');
title(' X4 optmal V.S X4 ADP');
xlabel('Time [sec]'); ylabel('X_4 [m/sec]');

subplot(3,2,5)
hold on;
plot(out_optimal.x.Time, out_optimal.x.Data(:,5), 'r-');
plot(out_ADP.x.Time, out_ADP.x.Data(:,5),'b-.', 'LineWidth', 1.25);
grid on;
legend('K optimal', 'K ADP');
title(' X5 optmal V.S X5 ADP');
xlabel('Time [sec]'); ylabel('X_5 [m/sec]');

subplot(3,2,6)
hold on;
plot(out_optimal.x.Time, out_optimal.x.Data(:,6), 'r-');
plot(out_ADP.x.Time, out_ADP.x.Data(:,6),'b-.', 'LineWidth', 1.25);
grid on;
legend('K optimal', 'K ADP');
title(' X6 optmal V.S X6 ADP');
xlabel('Time [sec]'); ylabel('X_6 [m/sec]');
%% destruction
figure(6)
set(gcf,'color','w');
hold on;
plot(out_ADP.input);
title('destruction V.S Time');
xlabel('Time [sec]'); ylabel('destruction');
grid on;
%% eror between LQR based model end LQR based data set
figure(7)
set(gcf,'color','w');
hold on;
plot(out_optimal.y.Time, out_optimal.y.Data - out_ADP.y.Data);
xlabel('Time [sec]'); ylabel('error [m]');
grid on;
title('eror between LQR based model end LQR based data set');
%%
function [K0, K, A, B, Ts, C]=ADP()
clc;
x_save=[];
t_save=[];
flag=1;  % 1: learning is on. 0: learning is off.
%% build system
%syms k_12 k_23 m1 m2 m3 gama1 gama2;
k_12 = 800;
k_23 = 100;
m1 = 2;
m2 = 1.5;
m3 = 1;
Ts = 0.01;
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

%% build contriller
gama1 = 0.001;
gama2 = 0.001;

Q = C'*C;
R = [gama1 0; 0 gama2];
[K0, P0] = lqr(A, B, Q, R); % Calculate the ideal solution for comparion purpose
eigen_LQR_close = eig(A - B*K0);
rank_LQR_close = rank(A - B*K0);

[xn,un]=size(B);%size of B. un-column #, xn row #

% Initialize the feedback gain matrix
K=0.5 * K0; %K = zeros(un,xn) Only if A is Hurwitz, K can be set as zero.
N=200;  %Length of the window, should be at least greater than xn^2
NN=10;  %Max iteration times
T=.01;  %Duration of time for each integration

%x0=[10;2;100;2;-1;-2]; %Initial condition
x0=[0;0;0;0;0;0];
i1=(rand(1,100)-.5)*5;
i2=(rand(1,100)-.5)*5;

Dxx=[];XX=[];XU=[];  % Data matrices

X=[x0;kron(x0',x0')';kron(x0,zeros(un,1))]';

% Run the simulation and obtain the data matrices \delta_{xx}, I_{xx},
% and I_{xu}

for i=1:N
    % Simulation the system and at the same time collect online info.
    [t,X]=ode45(@mysys, [(i-1)*T,i*T],X(end,:));

    %Append new data to the data matrices
    Dxx=[Dxx;kron(X(end,1:xn),X(end,1:xn))-kron(X(1,1:xn),X(1,1:xn))];
    XX=[XX;X(end,xn+1:xn+xn^2)-X(1,xn+1:xn+xn^2)];
    XU=[XU;X(end,xn+xn^2+1:end)-X(1,xn+xn^2+1:end)];

    % Keep track of the system trajectories
    x_save=[x_save;X];
    t_save=[t_save;t];
end

Dxx=processing_Dxx(Dxx); % Only the distinct columns left
% K=zeros(un,xn);  % Initial stabilizing feedback gain matrix
K = 5*K0;
P_old=zeros(xn);P=eye(xn)*10; % Initialize the previous cost matrix
it=0;            % Counter for iterations
iter=[it];
p_save=[];       % Track the cost matrices in all the iterations
k_save=[];       % Track the feedback gain matrix in each iterations
k_save=[norm(K-K0)];
norm_p_i = [norm(P)];
K11=[K(1,1)]; K12=[K(1,2)]; K13=[K(1,3)]; K14=[K(1,4)]; K15=[K(1,5)]; K16=[K(1,6)];
K21=[K(2,1)]; K22=[K(2,2)]; K23=[K(2,3)]; K24=[K(2,4)]; K25=[K(2,5)]; K26=[K(2,6)];
while norm(P-P_old)>0.001 && it < 20    % Stopping criterion for learning
    it=it+1;                        % Update and display the # of iters
    iter=[iter it];
    P_old=P;                         % Update the previous cost matrix
    QK=Q+K'*R*K;                     % Update the Qk matrix
    X2=XX*kron(eye(xn),K');          %
    X1=[Dxx,-X2-XU];                 % Left-hand side of the key equation
    Y=-XX*QK(:);                    % Right-hand side of the key equation
    pp=X1\Y;                         % Solve the equations in the LS sense
    P=reshape_p(pp);                 % Reconstruct the symmetric matrix
    p_save=[p_save,norm(P-P0)];      % Keep track of the cost matrix
    norm_p_i = [norm_p_i norm(P)];
    BPv=pp(end-(xn*un-1):end);
    K=inv(R)*reshape(BPv,un,xn)/2;   % Get the improved gain matrix
    K11=[K11 K(1,1)]; K12=[K12 K(1,2)]; K13=[K13 K(1,3)]; K14=[K14 K(1,4)]; K15=[K15 K(1,5)]; K16=[K16 K(1,6)];
    K21=[K21 K(2,1)]; K22=[K22 K(2,2)]; K23=[K23 K(2,3)]; K24=[K24 K(2,4)]; K25=[K25 K(2,5)]; K26=[K26 K(2,6)];
    k_save=[k_save,norm(K-K0)];     % Keep track of the control gains
end
% The following nested function gives the dynamics of the sytem. Also,
% integraters are included for the purpose of data collection.
    function dX=mysys(t,X)
        %global A B xn un i1 i2 K flag
        x=X(1:xn);
        if t>=2   % See if learning is stopped
            flag=0;
        end
        if flag==1
            u=zeros(un,1);
            for i=i1
                u(1)=u(1)+sin(i*t)/length(i1); % constructing the
                % exploration noise
            end
            for i=i2
                u(2)=u(2)+sin(i*t)/length(i2);
            end
            u=u;
        else
            u=-K*x;
        end
        dx=A*x+B*u;
        dxx=kron(x',x')';
        dux=kron(x',u')';
        dX=[dx;dxx;dux];
    end
%% This nested function reconstruct the P matrix from its distinct elements
    function P=reshape_p(p)
        P=zeros(xn);
        ij=0;
        for i=1:xn
            for j=1:i
                ij=ij+1;
                P(i,j)=p(ij);
                P(j,i)=P(i,j);
            end
        end
    end

function [ A,RHS ] = Normalizing( A,RHS)
% Normalize A,RHS by diag(A) in order to lower cond(A) number.
RHS=RHS./diag(A);
A=bsxfun(@rdivide, A, diag(A));
end

%% The following nested function removes the repeated columns from Dxx
function Dxx=processing_Dxx(Dxx)
    ij=[]; ii=[];
    for i=1:xn
        ii=[ii (i-1)*xn+i];
    end
     for i=1:xn-1
        for j=i+1:xn
            ij=[ij (i-1)*xn+j];
        end
     end
     Dxx(:,ii)=Dxx(:,ii)/2;
     Dxx(:,ij)=[];
     Dxx=Dxx*2;
end
%% control gain V.S iteration
figure(1)
set(gcf,'color','w');
subplot(2,1,1)
plot(iter, K11, 'o-', iter, K0(1,1)*ones(size(K11)),...
    iter, K12, 'o-', iter, K0(1,2)*ones(size(K12)),...
    iter, K13, 'o-', iter, K0(1,3)*ones(size(K13)),...
    iter, K14, 'o-', iter, K0(1,4)*ones(size(K14)),...
    iter, K15, 'o-', iter, K0(1,5)*ones(size(K15)),...
    iter, K16, 'o-', iter, K0(1,6)*ones(size(K16)));
title('K(1,:), control gain V.S iteration');
legend('K11', 'K11 optimal',...
    'K12', 'K12 optimal',...
    'K13', 'K13 optimal',...
    'K14', 'K14 optimal',...
    'K15', 'K15 optimal',...
    'K16', 'K16 optimal');
grid on;

subplot(2,1,2)
plot(iter, K21, 'o-', iter, K0(2,1)*ones(size(K21)),...
    iter, K22, 'o-', iter, K0(2,2)*ones(size(K22)),...
    iter, K23, 'o-', iter, K0(2,3)*ones(size(K23)),...
    iter, K24, 'o-', iter, K0(2,4)*ones(size(K24)),...
    iter, K25, 'o-', iter, K0(2,5)*ones(size(K25)),...
    iter, K26, 'o-', iter, K0(2,6)*ones(size(K26)));
title('K(2,:), control gain V.S iteration');
legend('K21', 'K21 optimal',...
    'K22', 'K22 optimal',...
    'K23', 'K23 optimal',...
    'K24', 'K24 optimal',...
    'K25', 'K25 optimal',...
    'K26', 'K26 optimal');
grid on;
%% norm P_i V.S iteration
figure(2)
hold on;
set(gcf,'color','w');
plot(iter, norm_p_i, 'r-o');
plot(iter, norm(P0)*ones(size(norm_p_i)), 'b-');
legend('||P_i ADP||', '||P_i optimal||');
xlabel('number of iteration'); ylabel('||P||');
title('||P_i|| V.S iteration');
grid on;
end