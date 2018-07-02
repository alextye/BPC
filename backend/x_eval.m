function y = x_eval(obs_data, N_sig)

    %function generates a set of x values that correspond with age
    %observations to be queried during the MCMC process for likelihood
    %calculation
    
    %PARAMETERS
    %obs_data - a 2 column matrix with age in the first column and 1-sigma
    %analytical uncertainty in the second
    
    %N_sig - the number of points to query for each standard deviation of
    %analytical uncertainty in a measured age, i.e. if an age is 100 +/- 1,
    %then an N_sig value of 10 would result in evaluation every 0.1 from -4
    %sigma to +4 sigma, in this case 96 to 104. An N_sig value of 1 in the
    %same situation would result in 1 query point per standard deviation
    %(i.e. the set of values 96, 97, 98, 99, ...).
    
    %OUTPUT
    
    %y is the x locations to be queried later on.

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
    end
    %generate and fill matrix x_eval, which contains the x locations for
    %marginal analysis of each datapoint.

    if sum(obs_data_s)==0
        y = obs_data_m;
    else
        %if analytical uncertainties are included to be modeled, each row of
        %the x_eval matrix represents a given zircon grain age, and each
        %column represents an x value to query for the purposes of using
        %marginal analysis to account for analytical uncertainty.
        %The locations of these query points are determined by the analytical 
        %uncertainty of each observation.
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
end
