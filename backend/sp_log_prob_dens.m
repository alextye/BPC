function y = sp_log_prob_dens(coefs, logsigfrac, x_eval, marg_basis, lnorm_densrep)
    
%function sp_log_prob calculates the log probability of observing the
%given data in obs_data given the candidate PDF described by knots and
%coefs.
    
%PARAMETERS
%coefs - coefficients of the spline basis functions.
%N_sig - the interval at which to discretize the evaluation of the
%   candidate PDF and Gaussian analytical uncertainty, in terms of 
%   the analytical uncertainty of each grain (i.e. if N_sig = 1, then
%   the discrete evaluation points for each grain will be age-4*unc, 
%   age-3*unc, age-2*unc, age-unc, age, age+unc, age+2*unc, age+3*unc, 
%   and age+4*unc; if N_sig = 2, then the discrete evaluation points 
%   will be age-4*unc, age-3.5*unc, age-3*unc, age-2.5*unc, etc., where
%   age is the measured age for each grain and unc is the 1 sigma 
%   analytical uncertainty.
%x_eval - the matrix of discrete x-values at which to evaluate
%   analytical uncertainty and candidate PDF function value, generated
%   by the function x_eval.
%marg_basis - matrix that contains the impulse response of each x location
%to be queried to each coefficient
    
%OUTPUT
%y is the log likelihood of observing the set of zircon ages contained in
%obs_data given the candidate PDF constructed from knots and coefs.

    %retrieve the function values (log likelihood values) of the candidate
    %PDF at each x value that corresponds to one of the discretization
    %points for the Gaussians associated with each individual zircon grain.
    
    %coefmat = repmat(reshape(coefs,1,1,[]),size(x_eval,1),size(x_eval,2));
    %vals = sum(marg_basis.*coefmat,3);
    vals = reshape(marg_basis,[],length(coefs))*coefs';
    vals = reshape(vals, size(marg_basis,1),[]);
    
    
    
    
    %the likelihood of observing each individual grain age given the
    %candidate PDF is the sum of the products of the Gaussian
    %function values (analytical uncertainty of each grain age) and
    %candidate PDF function values (in linear space) at a set of discrete evaluation points. This
    %line performs that operation (which involves exponentiating and then
    %taking the log) in a way that preserves accuracy of values in linear
    %form. Once the likelihoods of observing each grain age are again in
    %logarithmic form, log(sig_frac) is added because we use a column
    %method to approximate the integral of 
    %(Gaussian analytical uncertainty * candidate PDF), and the width of
    %the columns is sig_frac, so the discretized results must be multiplied
    %by sig_frac to be a true area, which is what is achieved by
    %integration.
    if size(vals,2)>1
        %find the maximum log likelihood value associated with each zircon
        %grain (each row of x_eval and vals matrices correspond to one zircon).
        maxval = max(vals, [], 2);

        %replicate the above column vector to be a matrix of the same
        %dimensions as vals.
        maxvalrep = repmat(maxval, 1, size(vals,2));
        
        sp_test_y = maxval + logsigfrac + log(sum(exp(vals + lnorm_densrep - maxvalrep), 2));
    else
        sp_test_y = vals;
    end
    
    %keyboard
    %add together the log likelihood of observing each individual grain age
    %to obtain the likelihood of observing the entire sample.
    y = sum(sp_test_y);
    
end