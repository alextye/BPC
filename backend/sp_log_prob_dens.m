function y = sp_log_prob_dens(knots, coefs, obs_data, N_sig, x_eval)
    
%function sp_log_prob calculates the log probability of observing the
%given data in obs_data given the candidate PDF described by knots and
%coefs.
    
%PARAMETERS
%knots - knots of the b-spline curve.
%coefs - coefficients of the spline basis functions.
%obs_data - observed zircon ages and analytical uncertainties.
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
    
%OUTPUT
%y is the log likelihood of observing the set of zircon ages contained in
%obs_data given the candidate PDF constructed from knots and coefs.
    
    
    %parse observed zircon data into ages (m) and uncertainties (s).
    
    obs_data_m = obs_data(:,1);
    obs_data_s = obs_data(:,2);
    
    %define the fraction of 1 sigma over which to discretize the normal
    %distribution for each point.
    sig_frac = 1/N_sig;
    
    sp_test = spmak(knots, coefs);
    
    %calculate the log of the function value of a Gaussian distribution
    %from -4 sigma to 4 sigma, at interval sig_frac.
    lnorm_dens = [-4:sig_frac:4];
    lnorm_dens = log(1/sqrt(2*pi)) - 0.5 * lnorm_dens.^2;

    %retrieve the function values (log likelihood values) of the candidate
    %PDF at each x value that corresponds to one of the discretization
    %points for the Gaussians associated with each individual zircon grain.
    vals = fnval(sp_test, x_eval);
        
    %find the maximum log likelihood value associated with each zircon
    %grain (each row of x_eval and vals matrices correspond to one zircon).
    maxval = max(vals, [], 2);
    
    %replicate the above column vector to be a matrix of the same
    %dimensions as vals.
    maxvalrep = repmat(maxval, 1, size(vals,2));
    
    %replicate the log Gaussian function values found above to have the
    %same dimensions as vals.
    lnorm_densrep = repmat(lnorm_dens, size(vals,1), 1);
    
    %initialize vector containing the log likelihood of each observed grain
    %age given the candidate PDF.
    sp_test_y = zeros(length(obs_data_m),1);
    
    %the likelihood of observing each individual grain age given the
    %candidate PDF is the (linear) sum of the products of the Gaussian
    %function values (analytical uncertainty of each grain age) and
    %candidate PDF function values at each discrete evaluation point. This
    %line performs that operation (which involves exponentiating and then
    %taking the log) in a way that preserves accuracy of values in linear
    %form. Once the likelihoods of observing each grain age are again in
    %logarithmic form, log(sig_frac) is added because we use a column
    %method to approximate the integral of 
    %(Gaussian analytical uncertainty * candidate PDF), and the width of
    %the columns is sig_frac, so the discretized results must be multiplied
    %by sig_frac to be a true area, which is what is achieved by
    %integration.
    sp_test_y = maxval + log(sig_frac) + log(sum(exp(vals + lnorm_densrep - maxvalrep), 2));
    
    %add together the log likelihood of observing each individual grain age
    %to obtain the likelihood of observing the entire sample.
    y = sum(sp_test_y);
    
end