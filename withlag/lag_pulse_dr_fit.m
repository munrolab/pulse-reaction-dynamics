function dr = lag_pulse_dr_fit(x,xdata)
    kx = x(1);
    tau = x(2);
    kr = x(3);
    m= x(4);
    t = xdata(1:length(xdata)/3);
    p = xdata(length(xdata)/3+1:2*length(xdata)/3);
    r = xdata(2*length(xdata)/3+1:end);
    dr = kx*(exp(-t/tau).*cumtrapz(t,exp(t/tau).*p)).^m-kr*r;
end