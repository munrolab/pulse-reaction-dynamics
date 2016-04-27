# pulse-reaction-dynamics

This is the fitting and modeling code used for paper X.  We used two fitting schemes (1d and 2d) and two model equations (simple and rhotimes).  pulseon_simple visualizes the fitting functions and the simulation output for the case where 

$$ dp/dt = k_{qq}\frac{p^n}{p0^n + p^n} - k_{off}r$$

$$ dr/dt = k_{s}p^m - k_{t}r$$

pulseon_rhotimes using the slightly more conventional, but slightly less accurate form

$$ dp/dt = k_{qq}\frac{p^n}{p0^n + p^n} - k_{off}pr$$

$$ dr/dt = k_{s}p^m - k_{t}r$$


The set of scripts comparor_* looks at the behavior of each model and fitting scheme when utilizing different values of n and m.

More documentation will be forthcoming as we find out how much is going to be used.
