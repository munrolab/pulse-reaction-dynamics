% 1D fit params 
% ka = 0.1;
% kn = 0.16;
% pa = 0.16;
% n = 2.7;
% kp = 0.12;
% kr = 0.1;
% m=1;

% 2D fit params (produces oscillator)
ka = 0.13;
kn = 0.17;
pa = 0.2;
n = 1.3;
kp = 0.12;
kr = 0.1;
m=1;

% adjustments to create pulses (still has a stable high state)
% ka = 0.11;
% kn = 0.17;
% pa = 0.18;
% n = 1.7;
% kp = 0.12;
% kr = 0.1;
% m=1;

% focus on early stage and add nonlinear p->r
% ka = 0.095;
% kn = 0.14;
% pa = 0.14;
% n = 2.7;
% kp = 0.22;
% kr = 0.18;
% m=1.6;

% set free parameters
bet = kp*kn/kr/kr*pa^(m-1);
alph = ka/pa/kr;

% determine rescalings
p_sc = pa;
r_sc = pa*kr/kn;
t_sc = 1/kr;

xq = linspace(-1,2/p_sc,100);
xs = linspace(-1,2/r_sc,100);
s_q = alph*(xq.^n)./(1+xq.^n);
s_s = bet*xq.^m;

figure;
plot(xq*p_sc,s_q*r_sc,'g');
hold on
plot(xq*p_sc,s_s*r_sc,'r');
[gq, gs] = meshgrid(linspace(-1,2/p_sc,20),linspace(-1,2/r_sc,20));
dq = alph*gq.^n./(1+gq.^n) - gs;
ds = (bet*gq.^m - gs) ;
% 
% quiver(gq*p_sc,gs*r_sc,dq*p_sc,ds*r_sc);
t = 0:0.01:100;
[~, sol] = ode45(@pulseon_simple_fun, t,[0.5,0],[],alph,bet,n,m);
% plot(sol(:,1)*p_sc,sol(:,2)*r_sc,'k');

% xlim([0 1]);
% ylim([0 1]);
% 
% linkaxes

% plot raw data
% plot(p(1:end-3),r(1:end-3));
% xlabel('Rho')
% ylabel('RGA')


plot(t*t_sc,sol(:,1)*p_sc)
hold on
plot(t*t_sc,sol(:,2)*r_sc);

% dt = find(max(p)==p)-t_sc*t(sol(:,1)==max(sol(:,1)))+1;
% plot(t*t_sc+dt,sol(:,1)*p_sc,'g');
% hold on
% plot(t*t_sc+dt,sol(:,2)*r_sc,'r');
% plot(p,'g-')
% plot(r,'r-')
% xlim([0,length(p)])