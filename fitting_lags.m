load('average_rho_rga.mat')

p = A(:,1);
r = A(:,2);
p = p-min(p);
r = r-mean(r(1:10));
% r(1:10)=0*r(1:10);
dp = diff(p);
dr = diff(r);
p = p(1:end-1);
r = r(1:end-1);


dp = dp(1:end-15);
dr = dr(1:end-15);
p = p(1:end-15);
r = r(1:end-15);


t = (0:1: length(p)-1)';

% cftool(p,r,dp);


x0 = [0.1 10 1 2];
xdata = [t; p; r];
ydata = dr;
x = lsqcurvefit(@lag_pulse_dr_fit,x0,xdata,ydata,[0 0 0 0]);

figure
plot3(p,r,dr)
hold on
plot3(p,r,lag_pulse_dr_fit(x,xdata),'.')

a = exp(-t/x(2)).*cumtrapz(t,exp(t/x(2)).*p);
a = a/max(a)*max(p);

figure;
plot(p);
hold on
plot(a);
plot(r);