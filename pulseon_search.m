global force;
force = false;
% kq = 0.5;
% kqq = 1;
% koff = 3;
% ks = 2;
n = 2;
stoP = [];
stoF = [];
stoR = [];
stoE = [];
for kN = 1%logspace(log10(0.3),log10(2),4)
    for ka = linspace(1,50,4)
        for k0 = linspace(2,5,3)
            h1 = figure;
%             h2 = figure;
            ind = 1;
            for kp = tan(linspace(pi/16,pi/2-pi/16,9))
%                 koff = kq+kq*kqq/2 - koff_set;
                xp = linspace(0,10,200);
                s_q = xp;
                s_s = kN*xp;
                figure(h1)
                subplot(3,3,ind)
                plot(xq,s_q,'g');
                hold on
                plot(xq,s_s,'r');
                x = [];
                y = [];
                u = [];
                v = [];
                i = 1;
                for q = linspace(0,5,20);
                    for s = linspace(-koff,kq*(1+kqq/2)*1.5,20);
                        x(i) = q;
                        y(i) = s;
                        dx = pulseon_fun(0,[q; s],kq,kqq,koff,ks,n,0);
                        u(i) = dx(1);
                        v(i) = dx(2);
                        i = i+1;
                    end
                end
                hold on
                quiver(x,(y+koff)/kq,u,v);
                t = 0:0.01:10;
                [~, osol] = ode45(@pulseon_fun, t,[0.01 0.01],[],kq,kqq,koff,ks,n,0);
                e = [0.1 0.2 0.4 0.8];
                holdit = [];
                for i = 1:4
                    force = true;
                    [~, sol] = ode45(@pulseon_fun, t,osol(end,:),odeset('MaxStep',0.01),kq,kqq,koff,ks,n,e(i)+osol(end,1));
                    holdit = [holdit e(i)+osol(end,1) max(sol(:,1)) sol(end,1)];
                        
                    figure(h1);
                    plot([osol(end,1);sol(:,1)],(koff+[osol(end,2);sol(:,2)])/kq,'k');
                    xlim([0 5]);
                    ylim([0 (1+kqq/2)*1.5]);
%                     figure(h2);
%                     subplot(2,2,i);
%                     plot(t,sol(:,1)/(osol(end,1)),'g');
%                     hold on
%                     plot(t,sol(:,2)/(osol(end,2)),'r');
                end
                stoF = [stoF; holdit(1:3:end)];
                stoR = [stoR; holdit(2:3:end)];
                stoE = [stoE; holdit(3:3:end)];
                stoP = [stoP; kq kqq koff ks];
                linkaxes
                ind = ind+1;
            end
        end
    end
end