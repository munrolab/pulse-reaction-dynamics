function dy = superpulse_1d_ode(t,y,x,m0,Dc,Dp,l,L,kp,koff,kc,n)
    t
    c = y(1:end/2);
    p = y(end/2+1:end);
    dx = x(2)-x(1);
    dc = zeros(size(c));
    dp = zeros(size(p));
    v = zeros(size(c));
%     m = m0*p;
%     Gr = (cosh((L+x(1)-x)/l)-cosh((x(1)-x)/l))/2/l^2/(cosh(L/l)-1);
%     v(1) = trapz(x,Gr.*m);
%     Gl = (cosh((L-x(end)+x)/l)-cosh((x(end)-x)/l))/2/l^2/(cosh(L/l)-1);
%     v(length(x)) = -trapz(x,Gl.*m);
%     for ind = 2:length(x)-1
%         Gr = (cosh((L+x(ind)-x)/l)-cosh((x(ind)-x)/l))/2/l^2/(cosh(L/l)-1);
%         Gl = (cosh((L-x(ind)+x)/l)-cosh((x(ind)-x)/l))/2/l^2/(cosh(L/l)-1);
%         v(ind) = trapz(x(ind:end),Gr(ind:end).*m(ind:end)) -  trapz(x(1:ind),Gl(1:ind).*m(1:ind));
%     end
    for ind = 1:length(v)
        left = ind-1;
        if(left<1)
            left = length(v);
        end
        right = ind+1;
        if(right>length(v))
            right = 1;
        end
        jL = (c(left)+c(ind))*(v(left)+v(ind))/4 - Dc*(c(ind)-c(left))/dx;
        jR = (c(right)+c(ind))*(v(right)+v(ind))/4 - Dc*(c(right)-c(ind))/dx;
        dc(left) = dc(left) - jL/2;
        dc(right) = dc(right) + jR/2;
        dc(ind) = dc(ind) + jL/2 - jR/2;
    end
    
    for ind = 1:length(v)
        left = ind-1;
        if(left<1)
            left = length(v);
        end
        right = ind+1;
        if(right>length(v))
            right = 1;
        end
        jL = (p(left)+p(ind))*(v(left)+v(ind))/4 - Dp*(p(ind)-p(left))/dx;
        jR = (p(right)+p(ind))*(v(right)+v(ind))/4 - Dp*(p(right)-p(ind))/dx;
        dp(left) = dp(left) - jL/2;
        dp(right) = dp(right) + jR/2;
        dp(ind) = dp(ind) + jL/2 - jR/2;
    end
    boost = 0;
    if(t>15&&t<20)
        boost = 1*double(x>L/2-2&x<L/2+2);
    end
    dc = dc + (kc*p - c);
    dp = dp + kp + boost  - (koff./(1+p.^n)+c).*p;
    dy = [dc;dp];
end