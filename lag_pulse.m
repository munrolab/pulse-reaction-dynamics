alph = 2;
bet = 5;
n = 2;
m = 1.35;

gam = 1.1*bet;

t = 0:0.1:100;
[t, c] = ode45(@lag_pulse_fun,t,[0.8 0 0],[],alph,bet,gam,n,m);


%plot(t,p(:,1),t,p(:,2))
figure;
plot(c(:,1),c(:,3),'b','LineWidth',2)
hold on

bnd=4;
cc = 0:0.25:bnd;
nq = alph*cc.^n./(1+cc.^n);
ns = (gam*cc).^m/bet;
plot(cc,nq,'r','LineWidth',2);
plot(cc,ns,'g','LineWidth',2);

[q, s] = meshgrid(cc,cc);
dq = alph*q.^n./(1+q.^n) - s;
dr = (gam*q).^m - bet*s;
quiver(q,s,dr,dq,'LineWidth',1)
xlim([0,bnd])
ylim([0,bnd])