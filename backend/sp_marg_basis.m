function [marg_basis, areax, area_basis] = sp_marg_basis(knots, N_coefs, x_eval, N_pts_area)

    %function generates several matrices that make it much more
    %computationally efficient to calculate PDF function values later on.
    
    %PARAMETERS
    
    %knots - knots of the spline
    
    %N_coefs - number of coefficients used for spline modeling
    
    %x_eval - the x locations which will need to be evaluated for
    %likelihood--these points include measured zircon ages and possibly
    %sets of nearby ages if the user wants to account for analytical
    %uncertainty thru marginalization
    
    %N_pts_area - number of points to spread evenly across the modeled x
    %scale for the purposes of evaluating the area under the PDF curve.
    %Typically 1000.
    
    %OUTPUT
    
    %marg_basis - matrix containing the impulse response of each location
    %in x_eval to each coefficient, i.e. if a coefficient's value were 1
    %and all other coefficients were 0, what would the overall function
    %value be at each x_eval location.  Dimensions and size depend on the
    %dimesions of x_eval, which will be 1d if no modeling of analytical
    %uncertainties is desired, or 2d if such modeling is desired.

    %areax - evenly spaced x values throughout the modeled domain queried
    %to calculate area under the PDF curves.
    
    %area_basis - impulse response of each x location in areax to each
    %coefficient, similar to marg_basis
    
    
    
    marg_basis = zeros(size(x_eval,1),size(x_eval,2),N_coefs);
    
    %areax = (min(knots)+.5*range(knots)/N_pts_area):range(knots)/N_pts_area:(max(knots)-.5*range(knots)/N_pts_area);
    
    areabreak = floor(5/7*(N_pts_area));
    areaxlog = [min(knots):(log(1000)-min(knots))/areabreak:log(1000)];
    areaxlin = [1000:(exp(max(knots))-1000)/(N_pts_area-areabreak):exp(max(knots))];
    areax = [areaxlog log(areaxlin(2:end))];
    %keyboard
        
    area_basis = zeros(length(areax),N_coefs);

    %figure
    %hold on
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
        
        marg_basis(:,:,i) = fnval(sp0,x_eval);
        area_basis(:,i) = fnval(sp0,areax);
        %plot(areax,fnval(sp0,areax));
    end
    %keyboard
end