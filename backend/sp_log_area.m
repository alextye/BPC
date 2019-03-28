function y = sp_log_area(coefs, x, area_basis)
    %approximate area under a log-probability plot (when translated to
    %linear space) using the trapezoid method.
    
    %PARAMETERS
    %coefs - spline coefficients
    %x - the set of x values for which function value is queried to
    %calculate area, is constant throughout a given run of pchain_gen (aka 
    %a single MCMC routine).
    %area_basis - a matrix consisting of the impulse response of each
    %value in x to the coefficients, i.e. what would the y value be at each
    %x location if a given coefficient were set to 1.  This value can be
    %multiplied by the coefficient value in the current model to determine
    %the x value given the current model.
    
    %OUTPUT
    %y - log area underneath the candidate PDF constructed from the spline
    %coefficients in coefs.

    %figure out each the function value at each x location by multiplying
    %the impulse responses of each x to each coefficient by the value of
    %that coefficient and summing across all coefficients
    %coefsmat = repmat(reshape(coefs,1,[]),size(area_basis,1),1);
    %fnvals = sum(area_basis.*coefsmat,2);
    fnvals = area_basis*coefs';
    %keyboard
    
    %trapezoidal approximation of the integral under the exponentiated
    %curve
%    y = sum(exp([fnvals(1) fnvals(2:end-1)' fnvals(2:end-1)' fnvals(end)]))/2*(x(2)-x(1))-1;
    y = log(trapz(x,exp(fnvals)));
    %y = log(sum(exp(fnvals))*(x(2)-x(1)));
end