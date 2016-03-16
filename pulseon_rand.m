ka = 0.16;
kn = 0.19;
n = 2;
pa = 0.31;
kp = 0.16;
kr = 0.14;
m = 1;
pc = 0.1;

bet = kp*kn/kr/kr/pa;
alph = ka/kr/pa;
q0 = pc/pa;

xq = linspace(0,10,30);
xs = linspace(0,10,30);
s_q = alph*(xq.^n)./(1+xq.^n);
s_s = bet*(xq.^m)./(q0^m+xq.^m);

figure;
plot(xq,s_q,'g');
hold on
plot(xq,s_s,'r');
[gq, gs] = meshgrid(xq,xs);
dq = alph*gq.^n./(1+gq.^n) - gs;
ds = bet*gq.^m./(q0^m+gq.^m) - gs ;

quiver(gq,gs,dq,ds);
t = 0:0.01:1000;
[~, sol] = ode45(@pulseon_fun, t,[0.8 0],[],alph,bet,q0,n,m);
plot(sol(:,1),sol(:,2),'k');

xlim([0 10]);
ylim([0 10]);

linkaxes

q_p = p/pa;
s_r = r*kn/kr/pa;
plot(q_p,s_r);
xlabel('Rho')
ylabel('RGA')