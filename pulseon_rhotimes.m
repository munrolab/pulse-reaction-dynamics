load('average_rho_rga.mat')
p_data = A(1:51,1)-A(1,1);
r_data = A(1:51,2)-A(1,2);
dp_data = diff(p_data);  %time happens to be unit 1, lucky us
dr_data = diff(r_data);
p_data = p_data(1:end-1);
r_data = r_data(1:end-1);


stp = 50;
spt = 20;
figure


n=2.4;
m=1.4;

% 1D fit for dp rising
dpf_1 = @(b,x)(b(1)*x(:,1).^n./(b(2)^n+x(:,1).^n));
beta = nlinfit(p_data(1:spt),dp_data(1:spt),dpf_1,[0.1,0.1]);
kqq = beta(1); p0 = beta(2); 
subplot(2,2,1)
plot(p_data(1:spt),dp_data(1:spt),'.')
hold on
plot(p_data(1:stp),dpf_1(beta,p_data(1:stp)))
xlabel('[Rho]')
ylabel('d[Rho]/dt')

% 1D fit for dp falling
dpf_2 = @(b,x)(b(1)-b(2).*x(:,1).*x(:,2));
beta = nlinfit([p_data(spt:stp),r_data(spt:stp)],dp_data(spt:stp),dpf_2,[0.1,0.1]);
koff = beta(2);
subplot(2,2,2)
plot(p_data(spt:stp).*r_data(spt:stp),dp_data(spt:stp),'.')
hold on
plot(p_data(spt:stp).*r_data(spt:stp),dpf_2(beta,[p_data(spt:stp),r_data(spt:stp)]))
xlabel('[Rho]\cdot[RGA]')
ylabel('d[Rho]/dt')

% 2D fit for dp
dpf = @(b,x)(b(1)*(x(:,1).^n)./(b(2)^n+x(:,1).^n)-b(3).*x(:,1).*x(:,2));
beta = nlinfit([p_data(1:stp),r_data(1:stp)],dp_data(1:stp),dpf,[0.1,0.1,0.1]);
kqq = beta(1); p0 = beta(2); koff = beta(3);
subplot(2,2,3)
plot3(p_data(1:stp),r_data(1:stp),dp_data(1:stp),'k.')
hold on
[X,Y] = meshgrid(0:0.1:1,0:0.1:1);
sol = dpf(beta,[X(:) Y(:)]);
s=surf(X,Y,reshape(sol,size(X)));
alpha(s,0.4)
xlabel('[Rho]')
ylabel('[RGA]')
zlabel('d[RGA]/dt')

%2D fit for dr
drf = @(b,x)(b(1)*x(:,1).^m-b(2)*x(:,2));
beta = nlinfit([p_data(spt:stp),r_data(spt:stp)],dr_data(spt:stp),drf,[0.1,0.1]);
ks = beta(1); kt = beta(2); 
subplot(2,2,4)
plot3(p_data(1:stp),r_data(1:stp),dr_data(1:stp),'k.')
hold on
[X,Y] = meshgrid(0:0.1:1,0:0.1:1);
sol = drf(beta,[X(:) Y(:)]);
s=surf(X,Y,reshape(sol,size(X)));
alpha(s,0.4)
xlabel('[Rho]')
ylabel('[RGA]')
zlabel('d[RGA]/dt')




% set free parameters
alph = kqq/p0/kt;
bet = ks*koff/kt^2*p0^m;

% determine rescalings
p_sc = p0;
r_sc = kt/koff;
t_sc = 1/kt;


xq = linspace(0,2/p_sc,100);
xs = linspace(0,2/r_sc,100);
s_q = alph*(xq.^(n-1))./(1+xq.^n);
s_s = bet*xq.^m;

figure;
plot(xq*p_sc,s_q*r_sc,'g');
hold on
plot(xq*p_sc,s_s*r_sc,'r');
[gq, gs] = meshgrid(linspace(-1,2/p_sc,20),linspace(-1,2/r_sc,20));
dq = alph*gq.^n./(1+gq.^n) - gq.*gs;
ds = (bet*gq.^m - gs);
% 
quiver(gq*p_sc,gs*r_sc,dq*p_sc,ds*r_sc);
t = 0:0.01:100;
[~, sol] = ode45(@pulseon_rhotimes_fun, t,[0.0,0.0],odeset('MaxStep',0.1),alph,bet,n,m);
plot(sol(:,1)*p_sc,sol(:,2)*r_sc,'k');

xlim([-1*p_sc,2]);
ylim([-1*r_sc,2]);
% 
% linkaxes

%plot raw data
plot(p_data(1:stp),r_data(1:stp));
xlabel('Rho')
ylabel('RGA')

figure
plot(t*t_sc,sol(:,1)*p_sc)
hold on
plot(t*t_sc,sol(:,2)*r_sc);

figure
dt = find(max(p_data)==p_data)-t_sc*t(sol(:,1)==max(sol(:,1)))+1;
plot(t*t_sc+dt,sol(:,1)*p_sc,'g');
hold on
plot(t*t_sc+dt,sol(:,2)*r_sc,'r');
plot(p_data,'go')
plot(r_data,'ro')
xlim([0,40])