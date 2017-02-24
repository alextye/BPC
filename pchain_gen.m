function [y,fval,hessian,pchain,logLk,accRatio] = pchain_gen(obs_data, knots, priors_mu, priors_sigma, N_pts, MCMC_FLAG,varargin)
%function solves for the maximum likelihood b-spline probability model to
%describe a set of observed zircon ages. The function then optimally
%executes a Metropolis-Hastings Markov Chain Monte Carlo (MCMC) search to
%describe the posterior probability distribution in model space that
%describes the set of models (possible parent populations) likely to have
%generated the observed data.

%PARAMETERS
%obs_data - observed zircon ages (column 1) and analytical uncertainties
%   (column 2)
%knots - knots for b-spline curve.
%priors_mu - expected value of each spline coefficient.
%priors_sigma - covariance matrix for the spline coefficients.
%N_pts - number of points to use in discretization of b-spline curves when
%   area underneath a candidate PDF is calculated.

%OUTPUT
%y is the set of spline coefficients representing the maximum likelihood
%   model, which is the candidate PDF that maximizes the likelihood of
%   observing the zircon age data in obs_data
%fval is the log likelihood * (-1) of the maximum likelihood model
%hessian is the Hessian derived from solving for the maximum likelihood
%   model by minimization.
%pchain is the Markov chain containing the sets of spline coefficients
%   accpeted by the Metropolis-Hastings algorithm. pchain is a NxM matrix
%   where N is the length of the Markov chain (equivalently, the number of
%   samples of the posterior) and M is the number of spline parameters used
%   (typically 50).
%logLk is the log likelihood values of the models contained in pchain.
%   logLk is a column vector where element P is the log likelihood of the
%   model whose coefficients are listed in row P of pchain.
%accRatio is the acceptance ratio of candidate models into the Markov
%   chain, as determined by the Metropolis-Hastings method.
%varargin can be fileID, the handle of an open log file.
    
    if min(size(varargin))>0
        fileID = varargin{1};
    else
        fileID = fopen('dump.txt','w');
    end

    %define degrees of freedom and number of parameters. 5 degrees of
    %freedom has been shown to be optimal for Cauchy distributions in
    %Markov processes.
    dof = 5;
    Np = length(priors_mu);

    %calculate x-values for function evaluation of candidate probability
    %models (PDFs) and the Gaussians that correspond to the measured age
    %and analytical uncertainty of each grain
    x_evalu = x_eval(obs_data, 10);

    %execute a constrained minimization to find the maximum likelihood PDF
    %for the given zircon sample. The minimized function is
    %-log(likelihood). The minimized function has an additional constraint,
    %which is that the area underneath the exponentiated b-spline curve 
    %must integrate to 1 because it is a PDF.
    [y,fval,exitflag,output,lambda,grad,hessian] = fmincon(@log_eval_opt,priors_mu,[],[],[],[],[],[],@log_area_loc);

    opt_area = sp_log_area(knots,y,N_pts);
    fprintf(fileID, '%s\n\n', strcat('opt_area = ',mat2str(opt_area)));

    %use the nearest symmetric, positive-definite matrix to the inverse of
    %the Hessian returned by the minimization as the covariance matrix for
    %the model parameters during the subsequent Monte Carlo simulations.
    opt_cov = nearestSPD(inv(hessian));
    fprintf(fileID, '%s\n\n', strcat('MLM = ',mat2str(y)));
    fprintf(fileID, '%s\n\n', strcat('cov matrix = ',mat2str(opt_cov)));

    %check and modify the covariance matrix so that it really reflects the
    %distances required to decrease the maximum likelihood by half in each
    %direction.
    [opt_cov, Dcheck] = mvncovcheck(y,opt_cov,@log_eval_opt);

    opt_cov = nearestSPD(opt_cov);
    fprintf(fileID, '%s\n\n', strcat('corrected cov matrix = ',mat2str(opt_cov)));
    fprintf(fileID, '%s\n\n', strcat('Dcheck = ',mat2str(Dcheck)));

    if(MCMC_FLAG==1)
        %only run if MCMC is desired, otherwise return only the maximum
        %likelihood model
        
        nc_flag = 1;
        fedup = 1;
        while(nc_flag==1|fedup==1)
            %find ideal step size for MCMC algorithm--the acceptance ratio of a
            %Metro-Hastings MCMC algorithm can be modeled as a logistic
            %function of the log of the step size for the MCMC algorithm.
        
            %find 2 points of reasonable acceptance ratio (>0.01) in order to
            %estimate the logistic function that relates acceptance ratio to
            %log step size.
            testt1 = 5/(10^(1/4));
            testt2 = 5;
            accRatio1 = 0;
            accRatio2 = 0;
            while(~(accRatio1 > 0.1 & accRatio2 > 0.1) | max(accRatio1,accRatio2)<0.24)
                testt1 = testt1/(10^(1/4));
                testt2 = testt2/(10^(1/4));
                [pchain1, logLk1, accRatio1] = mcmcmetro_bpdf(2000, y, @log_eval, testt1*opt_cov, knots, N_pts);
                [pchain2, logLk2, accRatio2] = mcmcmetro_bpdf(2000, y, @log_eval, testt2*opt_cov, knots, N_pts);
                fprintf(fileID, '%s\n\n', strcat('testt1 = ',mat2str(testt1)));
                fprintf(fileID, '%s\n\n', strcat('testt2 = ',mat2str(testt2)));
                fprintf(fileID, '%s\n\n', strcat('accRatio1 = ',mat2str(accRatio1)));
                fprintf(fileID, '%s\n\n', strcat('accRatio2 = ',mat2str(accRatio2)));
            end

            %use the step sizes and resulting acceptance ratios found above to
            %find the relation between acceptance ratio and log step size, then
            %use the relation to solve for the step size necessary to achieve
            %the desired acceptance ratio (0.24)
            x1 = log(testt1);
            x2 = log(testt2);
            solt = 0;
            if accRatio1>0.24 & accRatio2>0.24
                solt = max(x1,x2);
            else
                accRatiofctn = @(f)(((accRatio1 - (1-exp(f(1)*x1+f(2))))^2 + (accRatio2 - (1-exp(f(1)*x2+f(2))))^2)^(1/2));
                sol = patternsearch(accRatiofctn,[0 0]);
                targetAccRatio = @(t)((0.24 - (1-exp(sol(1)*t+sol(2))))^2);
                solt = patternsearch(targetAccRatio,0);
            end
            fprintf(fileID, '%s\n\n', strcat('exp(solt) = ',mat2str(exp(solt))));

            %find the autocorrelation and resulting 'effective sample size'
            %from initial MCMC run of 2000 steps, then initiate the main part
            %of the MCMC algorithm, which will run until the number of
            %'effective samples' is judged to be 30 in the direction of each
            %vector of the basis of the model covariance matrix.
            [pchain, logLk, accRatio] = mcmcmetro_bpdf(2000, y, @log_eval, exp(solt)*opt_cov, knots, N_pts);
            fprintf(fileID, '%s\n\n', strcat('test accRatio = ',mat2str(accRatio)));

            [effSampN,effSampNErr,logEffSampN,logEffSampNErr,AR1coeff,nc_flag,fedup] = mcmcautocorr_bpdf(pchain);
            fprintf(fileID, '%s\n\n', strcat('test [effSampN] = ',mat2str([effSampN])));
            fprintf(fileID, '%s\n\n', strcat('test [effSampNErr] = ',mat2str([effSampNErr])));
            fprintf(fileID, '%s\n\n', strcat('test [logEffSampN] = ',mat2str([logEffSampN])));
            fprintf(fileID, '%s\n\n', strcat('test [logEffSampNErr] = ',mat2str([logEffSampNErr])));
            fprintf(fileID, '%s\n\n', strcat('test [AR1coeff] = ',mat2str([AR1coeff])));
            fprintf(fileID, '%s\n\n', strcat('test [nc_flag] = ',mat2str([nc_flag])));
            fprintf(fileID, '%s\n\n', strcat('test [fedup] = ',mat2str([fedup])));
        end
        [pchain,logLk,effSampN,accVec] = mcmcsampsiz_bpdf(100,@log_eval,exp(solt)*opt_cov,pchain,logLk,effSampN,accRatio,'metro',knots,N_pts,fileID);
    else
        pchain = 0;
        logLk = 0;
        accRatio = 0;
    end
    
    
    %subtract log(prior) values from each log(posterior) value in order to
    %return log(likelihood) values.
    for i = 1:length(logLk)
        coef = pchain(i,:);
        logLk(i) = logLk(i) - sum(log(gamma((dof+Np)/2)) - log(gamma(dof/2)) - Np/2*log(pi) - Np/2*log(dof) - 1/2*log(det(priors_sigma)) - (dof+Np)/2*log(1+1/dof.*(coef-priors_mu)*inv(priors_sigma)*(coef-priors_mu)'));
    end
    
    function z = log_eval(coefs)

    %function which will be minimized during MCMC procedure. As with 
    %log_eval_opt below, this function generates a b-spline
    %function using the input coefficients and adds the function values at
    %each of the observed datapoints in obs_data, using marginal analysis
    %(procedure described in comments in function x_eval) to account for
    %analytical uncertainty of grain ages. These function values are the
    %log likelihoods of observing each individual zircon age given the
    %candidate PDF described by the coefficients (coefs). These log
    %likelihoods, after summation, are multiplied by -1 because the MATLAB
    %solver is a minimizer. Unlike log_eval_opt below, this function always
    %receives a set of coefficients that results in a candidate PDF with an
    %area of 1 underneath the curve, so no correction for area under the
    %candidate PDF is needed.

        %calculate the log likelihood of observing the given zircon sample
        %given the set of coefficients (coefs)
        log_prob = sp_log_prob_dens(knots, coefs, obs_data, 10, x_evalu);
        
        %add the log prior (a Cauchy distribution) to the log likelihood,
        %resulting in the log posterior; return the log posterior.

        z = (-1)*(log_prob + sum(log(gamma((dof+Np)/2)) - log(gamma(dof/2)) - Np/2*log(pi) - Np/2*log(dof) - 1/2*log(det(priors_sigma)) - (dof+Np)/2*log(1+1/dof.*(coefs-priors_mu)*inv(priors_sigma)*(coefs-priors_mu)')));
    
    end

    function z = log_eval_opt(coefs)

    %function which will be minimized during initial solution for the 
    %maximum likelihood model. This function generates a b-spline
    %function using the input coefficients and adds the function values at
    %each of the observed datapoints in obs_data, using marginal analysis
    %(procedure described in comments in function x_eval) to account for
    %analytical uncertainty of grain ages. These function values are the
    %log likelihoods of observing each individual zircon age given the
    %candidate PDF described by the coefficients (coefs). These log
    %likelihoods, after summation, are multiplied by -1 because the MATLAB
    %solver is a minimizer.

        %find the log of the area underneath the candidate PDF and subtract
        %it from each coefficient (results in a candidate PDF with area
        %underneath the curve of 1)
        
        %Note that the the log_area value is typically on the order of
        %10^(-6) or less, making this somewhat of a formality.
        log_area = sp_log_area(knots, coefs, N_pts);
        coefs = coefs - log_area;
        
        %calculate the log likelihood of observing the given zircon sample
        %given the set of coefficients (coefs)
        log_prob = sp_log_prob_dens(knots, coefs, obs_data, 10, x_evalu);
        
        %add the log prior (a Cauchy distribution) to the log likelihood,
        %resulting in the log posterior; return the log posterior.
        z = (-1)*(log_prob + sum(log(gamma((dof+Np)/2)) - log(gamma(dof/2)) - Np/2*log(pi) - Np/2*log(dof) - 1/2*log(det(priors_sigma)) - (dof+Np)/2*log(1+1/dof.*(coefs-priors_mu)*inv(priors_sigma)*(coefs-priors_mu)')));
        
    end

    function [c, ceq] = log_area_loc(coefs)
        
        %function returns the log of the area underneath the exponentiated
        %spline curve (a PDF). The output of log_area_loc is required to be
        %zero by the constrained minimizer for the maximum likelihood
        %solution spline coefficient set.
        c = [];
        ceq = (sp_log_area(knots, coefs, N_pts))^2;
        
    end
end