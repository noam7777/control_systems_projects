lamda_2 = 8.4;
lamda_1 = 4;
lamda_0 = 7;

TT = 1;

k1 = 10;
k2 = 10;

m = 10;
K = 4;
gama = 2;

k1_11= 10;
k2_11= 10;
m_11= 10;
x_11= 0.1;
xm_11 = 0.1;
e_11= 0.1;
yr_dt_11= 0.1;

h = tf(lamda_2, [1 lamda_1 lamda_2]);
hd = c2d(h,TT);
num = hd.Numerator{1,1};
dnum = hd.Denominator{1,1};

t = 0:20;
x = out.x;
x_m = out.x_m; 
e = out.e;
z = out.z;
m_l = out.m_l;
k1_l = out.k1_l;
k2_l = out.k2_l;
input = out.input;
u = out.u; 

figure 
plot(input);
hold on
grid on;
plot(x);
hold on
plot(x_m);
xlabel('t [sec]');
ylabel('x [m]');
legend('input', 'x', 'x_m');
hold off;

figure 
plot(k1_l);
hold on
plot(t, ones(length(t))*k1);
hold on
grid on;
xlabel('t [sec]');
ylabel('k1 [N/m]');
legend('k1', 'k1_real');
hold off;

figure 
plot(k2_l);
hold on;
plot(t, ones(length(t))*k2);
hold on
grid on
xlabel('t [sec]');
ylabel('k2 [N/m^3]');
legend('k2', 'k2_real');
hold off;

figure 
plot(m_l);
hold on
plot(t, ones(length(t))*m);
hold on
grid on;
xlabel('t [sec]');
ylabel('m [kg]');
legend('m', 'm_real');
hold off;

figure 
plot(z);
grid on;
xlabel('t [sec]');
ylabel('Z');
hold off;

figure 
plot(e);
grid on;
xlabel('t [sec]');
ylabel('e [m]');
hold off;

figure 
plot(u);
grid on;
xlabel('t [sec]');
ylabel('u [N]');
hold off;