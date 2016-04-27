load('average_rho_rga.mat')
p_data = A(1:51,1)-A(1,1);
r_data = A(1:51,2)-A(1,2);
dp_data = diff(p_data);  %time happens to be unit 1, lucky us
dr_data = diff(r_data);
p_data = p_data(1:end-1);
r_data = r_data(1:end-1);


stp = 50;
spt = 20;

nlist = [1.4 1.7 2];
mlist = [1.4 1.7 2];
% compare rhotimes
figure
ind=1;
for n = nlist
    for m = mlist
        ind
        subplot(length(nlist),length(mlist),ind)

        % 2D fit for dp
        dpf = @(b,x)(b(1)*(x(:,1).^n)./(b(2)^n+x(:,1).^n)-b(3).*x(:,1).*x(:,2));
        beta = nlinfit([p_data(1:stp),r_data(1:stp)],dp_data(1:stp),dpf,[0.1,0.1,0.1]);
        kqq = beta(1); p0 = beta(2); koff = beta(3);

        %2D fit for dr
        drf = @(b,x)(b(1)*x(:,1).^m-b(2)*x(:,2));
        beta = nlinfit([p_data(spt:stp),r_data(spt:stp)],dr_data(spt:stp),drf,[0.1,0.1]);
        ks = beta(1); kt = beta(2); 
        
        % set free parameters
        alph = kqq/p0/kt;
        bet = ks*koff/kt^2*p0^m;

        % determine rescalings
        p_sc = p0;
        r_sc = kt/koff;
        t_sc = 1/kt;

        
        [~, sol] = ode45(@pulseon_rhotimes_fun, t,[0.0,0.0],odeset('MaxStep',0.1),alph,bet,n,m);
        dt = find(max(p_data)==p_data)-t_sc*t(sol(:,1)==max(sol(:,1)))+1;
        plot(t*t_sc+dt,sol(:,1)*p_sc);
        hold on
        plot(t*t_sc+dt,sol(:,2)*r_sc);
        plot(p_data,'.')
        plot(r_data,'.')
        xlim([0,200])
        
        ind=ind+1;
    end
end


