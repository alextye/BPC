%Alex Tye
%Dec 2016
function [pchain, logLk, logPost, fileID] = chain_const_joint(N_coefs, knots, obs_data1, obs_data2, x_min, x_max, p_sig, delta, N_pts,filename,QUICKMFLAG)

%function constructs expected values and covariance matrix that are 
%parameters for the prior probability distribution and call function that
%samples the posterior.

%OUTPUT
%[pchain, logLk] are the spline coefficients and natural log likelihood
%values of the models of the Markov chain--a representative sample of the
%posterior

%PARAMETERS
%N_coefs is the number of spline coefficients to be used,
%obs_data1 is the set of observed ages and uncertainties of zircon sample 1
%obs_data2 is the set of observed ages and uncertainties of zircon sample 2
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
    %knots = [x_min x_min x_min:(x_max-x_min)/(N_coefs-2):x_max x_max x_max];

    %set expected coefficient values for a uniform distribution.
    priors_mu = zeros(1, N_coefs) -log(x_max-x_min);

    %find mean value of each basis function in order to construct prior.
    x = [[0:.01:6.9] log([1000:10:4000])];
    bas_mean = zeros(1,N_coefs);
    
    for i = 1:N_coefs
        
        %iterate thru each possible basis function coefficient, and
        %construct the spline curve that consists of only that basis
        %function (coefficient=1).  Then, determine the overall function
        %values at both the x locations that will be queried for likelihood
        %(depends on the zircon age observations) and the evenly spaced
        %locations that will be used to estimate area under the PDF.
        
        coefs = zeros(1,N_coefs);
        coefs(i) = 1;
        sp0 = spmak(knots, coefs);
        
        %calculate mean value by integration, being careful to normalize by
        %area
        bas_mean(i) = trapz(exp(x),fnval(sp0,x))/trapz(exp(x),fnval(sp0,x)./exp(x));
    end
    
    mmeandif = mean(bas_mean(2:end)-bas_mean(1:end-1));
    
    %construct a prior covariance matrix--diagonal elements are variance of
    %each spline parameter; off diagonal elements reflect influence of one
    %spline parameter on another. Magnitude of the off-diagonal elements is
    %controlled by parameter delta.

    priors_sigma = p_sig * ones(1, N_coefs);
    priors_sigma2 = priors_sigma.^2;

    cov_prior = diag(priors_sigma2);
    for i = 1:length(priors_sigma2)
        for j = i+1:length(priors_sigma2)
            %cov_prior(i,j) = priors_sigma(i) * priors_sigma(j) * exp(-abs(bas_mean(i)-bas_mean(j))/(mmeandif*delta));
            cov_prior(i,j) = priors_sigma(i) * priors_sigma(j) * exp(-abs(i-j)/delta);
            cov_prior(j,i) = cov_prior(i,j);
        end
    end

    obs_data_comp = [obs_data1; obs_data2];

    fileID = fopen(filename,'w');
    
    %call function to sample the posterior distribution that corresponds to
    %the 'joint' zircon sample, which includes the ages of both zircon
    %samples

    [best_coefs, bestLk, hessian, pchain, logLk, accRatio, pre_norm_area, logPost] = pchain_gen(obs_data_comp, knots, priors_mu, cov_prior, N_pts, 1, fileID,QUICKMFLAG);

    fclose(fileID);
    
end