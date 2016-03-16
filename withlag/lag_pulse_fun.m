function dc = lag_pulse_fun(t,c,alph,bet,gam,n,m)
    q = c(1);
    b = c(2);
    s = c(3);
    dc = [alph*q.^n./(1+q.^n) - s ;
        gam*q - b;
        b.^m - bet*s;];
end