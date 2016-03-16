function dx = pulseon_fun(t,x,alph,bet,q0,n,m)
    q = x(1);
    s = x(2);
    dx = [alph*q^n/(1+q^n)  - s;
        (bet*q^m/(q0^m+q^m) - s)];
    
end