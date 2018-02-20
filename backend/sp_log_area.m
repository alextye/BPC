function y = sp_log_area(knots, coefs, N_pts)
    %approximate area under a log-probability plot (when translated to
    %linear space) using the trapezoid method.
    
    %PARAMETERS
    %knots - knot locations for b-spline
    %coefs - spline coefficients
    %N_pts - number of points to use when discretizing the b-spline
    %   function
    
    %OUTPUT
    %y - log area underneath the candidate PDF constructed from the spline
    %coefficients in coefs.

    x = min(knots):range(knots)/N_pts:max(knots);
    
    %exponentiate the b-spline curve--the b-spline curve is in log
    %probability density - log age space.
    p = exp(fnval(spmak(knots, coefs), x));
    y = log(trapz(x, p));
    
    
end