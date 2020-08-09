lamda_2 = 8.4;
lamda_1 = 4;
lamda_0 = 7;

TT = 0.05;

k1 = 10;
k2 = 10;

m = 10;
K = 10;
gama = 15;

% k1_11= 10; k2_11= 10; m_11= 10; x_11= 0.1; xm_11 = 0.1; e_11= 0.1;
% yr_dt_11= 0.1;

h = tf(lamda_2, [1 lamda_1 lamda_2]);
hd = c2d(h,TT);
num = hd.Numerator{1,1};
dnum = hd.Denominator{1,1};

p_lead = -1*lamda_1/2;

z_lead = -2/TT;

k_lead = (TT/2)*((((11*lamda_1/2)^2 + abs(lamda_2 - 0.25*lamda_1^2))^0.5*(0.25*lamda_1^2 + abs(lamda_2 - 0.25*lamda_1^2))^0.5)/((0.5*lamda_1 + 2/TT)^2 + abs(lamda_2 - 0.25*lamda_1^2)));

d_lead_c = zpk(z_lead, p_lead, k_lead);

d_lead_d = c2d(d_lead_c, TT);

num_l = d_lead_d.Z{1};
dnum_l = d_lead_d.P{1};
k_lead_d = d_lead_d.K;



t = 0:20;
x_05 = out.x;
x_m_05 = out.x_m; 
e_005 = out.e;
z_1 = out.z;
m_l_1 = out.m_l;
k1_l_1 = out.k1_l;
k2_l_1 = out.k2_l;
input_05 = out.input;
u_05 = out.u; 
x_mc=out.x_mc;
% 
%  figure 
%  plot(input_001);
%  hold on
%  grid on;
%  plot(x_05);
%  hold on
%  plot(x_m_05);
%  %plot(x_01);
%  hold on
%  %plot(x_m_01);
%  %plot(x_mc);
%  xlabel('t [sec]');
%  ylabel('x [m]');
%  
%  legend('input_001', 'x_05', 'x_m_05');
% hold off;
% 
% figure 
% plot(k1_l);
% hold on
% plot(t, ones(length(t))*k1);
% hold on
% grid on;
% xlabel('t [sec]');
% ylabel('k1 [N/m]');
% legend('k1', 'k1_real');
% hold off;
% 
% figure 
% plot(k2_l);
% hold on;
% plot(t, ones(length(t))*k2);
% hold on
% grid on
% xlabel('t [sec]');
% ylabel('k2 [N/m^3]');
% legend('k2', 'k2_real');
% hold off;
% 
% figure 
% plot(m_l);
% hold on
% plot(t, ones(length(t))*m);
% hold on
% grid on;
% xlabel('t [sec]');
% ylabel('m [kg]');
% legend('m', 'm_real');
% hold off;
% 
% figure 
% plot(z);
% grid on;
% xlabel('t [sec]');
% ylabel('Z');
% hold off;
% 
figure 
hold on
plot(e_001);
plot(e_005);
grid on;
xlabel('t [sec]');
ylabel('e [m]');

% 
% figure 
% plot(u);
% grid on;
% xlabel('t [sec]');
% ylabel('u [N]');
% hold off;