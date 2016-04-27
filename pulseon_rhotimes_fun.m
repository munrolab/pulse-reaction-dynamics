function dx = pulseon_rhotimes_fun(t,x,alph,bet,n,m)
    q = x(1);
    s = x(2);
    dx = [5*(t<2.1&&t>2)+alph*q^n/(1+q^n)  - s*q ;
        (bet*q^m - s) ] ;
end