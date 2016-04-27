function dx = pulseon_simple_fun(t,x,alph,bet,n,m)
    q = x(1);
    s = x(2);
    dx = [2*(t<2.1&&t>2)+alph*q^n/(1+q^n)  - s ;
        (bet*q^m - s) ] ;
end