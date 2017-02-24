function y = x_eval(obs_data, N_sig)

%function generates a matrix of x-values at which to evaluate a candidate
%probability model (a PDF) in order to calculate the likelihood of
%observing a given zircon sample. This set of x-values permits marginal
%analysis of the observed ages and their uncertainties--the likelihood of a
%specific grain can be understood as the integral of the candidate PDF
%times the Gaussian distribution constructed from the grain age and
%analytical uncertainty. This integral is calculated by summing the
%function values of the candidate PDF times the Gaussian that represents
%the grain age and analytical uncertainty at a discrete interval.

%PARAMETERS
%obs_data is the set of zircon ages and analytical uncertainties for the
%   dataset
%N_sig is the interval at which to discretize the evaluation of the
%   candidate PDF and Gaussian analytical uncertainty, in terms of 
%   the analytical uncertainty of each grain (i.e. if N_sig = 1, then the
%   discrete evaluation points for each grain will be age-4*unc, age-3*unc,
%   age-2*unc, age-unc, age, age+unc, age+2*unc, age+3*unc, and age+4*unc;
%   if N_sig = 2, then the discrete evaluation points will be age-4*unc,
%   age-3.5*unc, age-3*unc, age-2.5*unc, etc., where age is the measured
%   age for each grain and unc is the 1 sigma analytical uncertainty.

%OUTPUT
%y is an N*(8*N_sig + 1) matrix consisting of x values at which to evaluate
%   the candidate PDF and Gaussians that correspond to the measured grain
%   ages and analytical uncertainties.  Each row corresponds to a different
%   grain and the column values are the x-values to evaluate for that
%   grain.

    %define the fraction of 1 sigma over which to discretize the normal
    %distribution for each point.
    sig_frac = 1/N_sig;

    %parse the 2 column matrix obs_data into mu (first column) and sigma 
    %values (second column).
    if(size(obs_data,2)==2)
        obs_data_m = obs_data(:,1);
        obs_data_s = obs_data(:,2);
    elseif(size(obs_data,2)==3)
        obs_data_m = obs_data(:,2);
        obs_data_s = obs_data(:,3);
    else
        'need 2 or 3 column matrix'
    end
    
    %generate and fill matrix x_eval2, which contains the x locations for
    %function evaluation, as described above.
    x_eval2 = zeros(length(obs_data_m), 8/sig_frac+1);
    for i = 1:length(obs_data_m)
        try
            x_eval2(i,:) = [obs_data_m(i) - 4 * obs_data_s(i) : obs_data_s(i) * sig_frac : obs_data_m(i) + 4 * obs_data_s(i)];
        catch
            keyboard;
        end
    end

    y = x_eval2;
end