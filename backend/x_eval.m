function y = x_eval(obs_data, N_sig)

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
    %generate and fill matrix x_eval, which contains the x locations for
    %marginal analysis of each datapoint.

    if sum(obs_data_s)==0
        y = obs_data_m;
    else

    %DIM1 = datapoint
    %DIM2 = locations for marginal analysis
    %i.e. each row represents one grain age observation and each column is
    %a location about that measurement to evaluate for marginal analysis.
    %The exact locations of these points are determined by the analytical 
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
