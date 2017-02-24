%Alex Tye
%Dec 2016
function [pchain, logLk] = chain_const(N_coefs, obs_data, x_min, x_max, p_sig, delta, N_pts,fileID)

%function constructs expected values and covariance matrix that are 
%parameters for the prior probability distribution and call function that
%samples the posterior.

%OUTPUT
%[pchain, logLk] are the spline coefficients and natural log likelihood
%values of the models of the Markov chain--a representative sample of the
%posterior

%PARAMETERS
%N_coefs is the number of spline coefficients to be used,
%obs_data is the set of observed ages and uncertainties
%x_min and x_max are the minimum and maximum x-values to be covered by the
%   spline function
%p_sig controls the variance of the spline coefficients enforced by the
%   prior
%delta controls the influence of each spline coefficient on neighboring
%   coefficients
%N_pts is the number of points by which to discretize candidate probability
%   models for likelihood calculations.
%fileID is the reference to the open log file.


    %Define knots for the cubic spline based on N_coefs and the range of 
    %values represented by the data.
    knots = [x_min x_min x_min:(x_max-x_min)/(N_coefs-2):x_max x_max x_max];

    %set expected coefficient values for a uniform distribution.
    priors_mu = zeros(1, N_coefs) -log(x_max-x_min);

    %construct a prior covariance matrix--diagonal elements are variance of
    %each spline parameter; off diagonal elements reflect influence of one
    %spline parameter on another. Magnitude of the off-diagonal elements is
    %controlled by parameter delta.

    priors_sigma = p_sig * ones(1, N_coefs);
    priors_sigma2 = priors_sigma.^2;

    cov_prior = diag(priors_sigma2);
    for i = 1:length(priors_sigma2)
        for j = i+1:length(priors_sigma2)
            cov_prior(i,j) = priors_sigma(i) * priors_sigma(j) * exp(-abs(i-j)/delta);
            cov_prior(j,i) = cov_prior(i,j);
        end
    end

    %call function to sample the posterior distribution that corresponds to
    %the current sample

    [best_coefs, bestLk, hessian, pchain, logLk, accRatio] = pchain_gen(obs_data, knots, priors_mu, cov_prior, N_pts, 1,fileID);

end