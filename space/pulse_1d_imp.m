
datalen = 101;
Dc = 1;
Dp = 1;
L = 100;
l = 5;

m0 = 50;

c0 = 0.7971;
p0 = 0.1328;
kp = 0.1;
kpp = 8;
koff = 1;
kcp = 6;
np = 2;

x = linspace(0,L,200)';
c = c0*(ones(size(x))-0.01*(rand(size(x))-0.5)-0.01*cos(2*pi*x/L));
p = p0*ones(size(x));
dt = 1;
[t, y] = ode23tb(@pulse_1d_ode,0:dt:dt*(datalen-1),[c; p],odeset('NonNegative',1:length(c),'RelTol',1e-3),x,m0,Dc,Dp,l,L,kp,kpp,koff,kcp,np);
c = y(:,1:end/2);
p = y(:,end/2+1:end);
x = x';

figure;
movind = 1;
for i = 1:datalen
    m = m0*c(i,:);
    Gr = (cosh((L+x(1)-x)/l)-cosh((x(1)-x)/l))/2/l^2/(cosh(L/l)-1);
    v(1) = trapz(x,Gr.*m);
    Gl = (cosh((L-x(end)+x)/l)-cosh((x(end)-x)/l))/2/l^2/(cosh(L/l)-1);
    v(length(x)) = -trapz(x,Gl.*m);
    for ind = 2:length(x)-1
        Gr = (cosh((L+x(ind)-x)/l)-cosh((x(ind)-x)/l))/2/l^2/(cosh(L/l)-1);
        Gl = (cosh((L-x(ind)+x)/l)-cosh((x(ind)-x)/l))/2/l^2/(cosh(L/l)-1);
        v(ind) = trapz(x(ind:end),Gr(ind:end).*m(ind:end)) -  trapz(x(1:ind),Gl(1:ind).*m(1:ind));
    end
    
%     subplot(2,1,2);
%     plot(x,v,'r.');
%     subplot(2,1,1);
    plot(x,c(i,:)/mean(mean(c)),'r');
    hold on
    plot(x,p(i,:)/mean(mean(p)),'g');
    ylim([0 1.2*max(max(max(c))/mean(mean(c)),max(max(p))/mean(mean(p)))]);
    
    drawnow
    mov(movind,:) = getframe;
    movind = movind+1;
    clf
end

% figure;
% movind = 1;
% for i = 1:datalen
%     [xx, yy ]=meshgrid(x,x);
%     
%     plot(x,c(i,:)/mean(mean(c)),'r');
%     hold on
%     plot(x,p(i,:)/mean(mean(p)),'g');
%     ylim([0 1.2*max(max(max(c))/mean(mean(c)),max(max(p))/mean(mean(p)))]);
%     
%     drawnow
%     mov(movind,:) = getframe;
%     movind = movind+1;
%     clf
% end